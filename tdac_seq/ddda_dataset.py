import mappy
import re
import numpy as np
from tqdm import tqdm
import os
from scipy.stats import norm, binom
import scipy
import multiprocessing as mp
import pickle
import tempfile
import logging
from collections.abc import Iterable
from .utils import count_reads_in_fastq, get_fastq_iter
from .refseq import extract_ref_seq, extract_ref_seq_from_fasta
from .dedup import dedup_read_inds
import pandas as pd

def get_footprint_single_scale(edits, bias, fp_radius=50, combine_sides="max_pval", pseudo_count_factor=0.01):
    '''
    Compute the footprint of a single read at a single scale (i.e., a single footprint radius)
    :param edits: Vector containing edits (1 is edited, 0 is unedited) at each position
    :param bias: Vector of the same length as edits containing bias of the corresponding strand at each position
    :param fp_radius: Int, radius of the footprint window
    :param combine_sides: One of "max_pval" (default), "mean_pval", and "min_pval". Determines how to combine center-vs-left and center-vs-right test statistics.
    :param pseudo_count_factor: Float. Used to determin the pseudo count when calculating center-vs-flank bias ratio.
    :return pvals: Vector of the same length as edits containing the footprint pvalues at each position
    '''

    # Create the left, center, and right flank convolution kernels
    # We will later run a sliding window sum of the edit and bias counts using these kernels
    # This is mathematically equivalent to convolution
    left_flank_kernel = np.concatenate([np.ones(fp_radius * 3), np.zeros(fp_radius)])
    center_kernel = np.concatenate([np.zeros(fp_radius), np.ones(fp_radius * 2), np.zeros(fp_radius)])
    right_flank_kernel = np.concatenate([np.zeros(fp_radius), np.ones(fp_radius * 3)])

    # Sliding window sum of left + center, center, and center + right edit counts
    sliding_sums = {}
    sliding_sums["left_edit"] = np.convolve(edits, left_flank_kernel)
    sliding_sums["center_edit"] = np.convolve(edits, center_kernel)
    sliding_sums["right_edit"] = np.convolve(edits, right_flank_kernel)

    # Sliding window sum of left + center, center, and center + right bias
    sliding_sums["left_bias"] = np.convolve(bias, left_flank_kernel)
    sliding_sums["center_bias"] = np.convolve(bias, center_kernel)
    sliding_sums["right_bias"] = np.convolve(bias, right_flank_kernel)

    # Compare center-to-(center + flank) ratio of edits versus the same ratio of biases
    # Use binomial test to calculate p-values
    pval_dict = {}
    for side in ["left", "right"]:
    
        # Choose an appropriate pseudo-count to prevent division by zero
        # We can't just use a fixed number becauese different bias models have different ranges
        side_sum = sliding_sums[side + "_bias"]
        pseudo_count = np.min(side_sum[side_sum > 0]) * pseudo_count_factor
    
        # Probability p used for binomial test
        p = sliding_sums["center_bias"] / (sliding_sums[side + "_bias"] + pseudo_count)  # + pseudo_count to prevent division by zero
    
        # Calculate pvals using binomial test
        pval_dict[side] = binom.cdf(
            k=sliding_sums["center_edit"],
            n=sliding_sums[side + "_edit"],
            p=p
        )

    # Keep the average p-value on the two sides
    if combine_sides == "min_pval":
        pvals_padded = np.minimum(pval_dict["left"], pval_dict["right"])
    elif combine_sides == "max_pval":
        pvals_padded = np.maximum(pval_dict["left"], pval_dict["right"])
    elif combine_sides == "mean_pval":
        pvals_padded = (pval_dict["left"] + pval_dict["right"]) / 2

    # Remove padded edges
    pvals = pvals_padded[(2 * fp_radius):(len(pvals_padded) - 2 * fp_radius + 1)]

    return pvals

def get_footprint_multi_scale(edits, bias, fp_radii, combine_sides="max_pval"):
    '''
    Compute the footprint of a single read at multiple scales (i.e., multiple footprint radii)
    :param edits: Vector containing edits (1 is edited, 0 is unedited) at each position
    :param bias: Vector of the same length as edits containing bias of the corresponding strand at each position
    :param fp_radii: List or array of footprint radii to use
    :return: Array of shape (len(fp_radii), length_seq, 1) containing the footprint scores at each position using each footprint radius
    '''

    footprints = []
    for fp_radius in fp_radii:
        footprints.append(get_footprint_single_scale(edits, bias, fp_radius, combine_sides))
    return np.array(footprints)

def align_helper(args):

    read_id, read_seq, ref_fasta_path, start_gap_threshold, \
        end_gap_threshold, del_size_threshold, \
        ref_seq_dict, mask_CpG = args

    aligner = mappy.Aligner(ref_fasta_path, preset='map-ont')
    
    # Align the current read
    alignments = [a for a in aligner.map(read_seq, cs=True)]
    
    # Remove multi-mappers and unmapped reads
    if len(alignments) == 0 or len(alignments) > 1:
        return None
    else:
        alignment = alignments[0]
    
    # Retrieve the locus that this read is aligned to
    locus = alignment.ctg
    seq_len = len(ref_seq_dict[locus])
    
    # Remove reads with large gaps at the beginning or end
    if alignment.r_st >= start_gap_threshold:
        return None
    if alignment.r_en <= seq_len - end_gap_threshold:
        return None
    
    # Parse information from the cs string (for details see https://lh3.github.io/minimap2/minimap2.html)
    cs_info = re.findall("([=:*+-~][\\dATCGNatcgn]+)", alignment.cs)
    
    # Encode C-to-T mutations and G-to-A mutations into a vector
    edit_vec = np.zeros((seq_len, 2), dtype=bool)
    del_vec = np.zeros(seq_len, dtype=bool)
    ABE_edit_vec = np.zeros(seq_len, dtype=bool)
    current_pos = alignment.r_st
    ref_seq = ref_seq_dict[locus]
    for i in cs_info:
        match i[0]:
            case ":":  # Match compared to the reference
                current_pos = current_pos + int(i[1:])
            case "-":  # Deletion compared to the reference
                del_size = len(i[1:])
                if del_size > del_size_threshold:  # Record deletions passing size threshold
                    del_vec[current_pos:(current_pos + del_size)] = True
                current_pos = current_pos + del_size
            case "+":  # Insertion compared to the reference
                pass
            case "*":  # Mismatch compared to the reference
                if mask_CpG:
                    is_CpG = (ref_seq[current_pos:(current_pos + 2)] == "CG") or \
                    (ref_seq[(current_pos - 1):(current_pos + 1)] == "CG")
                else:
                    is_CpG = False # IF we don't mask CpG, just treat every CpG as non-CpG
                if i[1] == "c" and i[2] == "t" and not is_CpG:
                    edit_vec[current_pos, 0] = True
                if i[1] == "g" and i[2] == "a" and not is_CpG:
                    edit_vec[current_pos, 1] = True
                if (i[1] == "a" and i[2] == "g") or (i[1] == "t" and i[2] == "c"):
                    ABE_edit_vec[current_pos] = True
                current_pos = current_pos + 1
    
    # Determine if the current read is a C-to-T strand or G-to-A strand
    strand_ind = np.argmax(np.sum(edit_vec, axis=0))

    # Also mask read edges because different reads have different edge gaps
    edit_vec[:start_gap_threshold] = False
    edit_vec[-end_gap_threshold:] = False
    ABE_edit_vec[:start_gap_threshold] = False
    ABE_edit_vec[-end_gap_threshold:] = False
        
    return {
        "del_vec" : scipy.sparse.csr_array(del_vec[np.newaxis]),
        "edit_vec" : scipy.sparse.csr_array(edit_vec[:, strand_ind][np.newaxis]),
        "locus" : locus,
        "strand" : strand_ind,
        "read_id" : read_id,
        "ABE_edit_vec" : scipy.sparse.csr_array(ABE_edit_vec[np.newaxis])
    }

class ddda_dataset:
    '''
    Class for storing DddA editing data
    '''

    def __init__(self, ID, region_dict: dict[str, tuple[str, int, int]], fastq_file, working_dir, genome_file: str | None = None, ref_fasta_path: str | None = None):
        '''
        Initialize the dda_dataset object
        :param ID: str, unique identifier for the dataset
        :param region_dict: dict, dictionary of regions to analyze. Format: {region_name: (chromosome, start, end)}
            Coordinates are 1-based and inclusive.
        :param genome_file: str, path to the reference genome file (e.g. hg38.fa)
        :param fastq_file: str, path to the fastq file containing nanopore sequencing reads
        :param working_dir: str, path to the working directory
        '''

        logging.info("Initializing ddda_dataset object")
        self.ID = ID
        self.region_dict = region_dict
        self.fastq_file = fastq_file
        self.working_dir = working_dir
        os.makedirs(self.working_dir, exist_ok=True)
        self.bias_dict = None
        self.bias_model = None

        if genome_file is not None:
            assert ref_fasta_path is None, "Either genome_file or ref_fasta_path should be provided, not both"
            self.ref_seq_dict = extract_ref_seq(genome_file, region_dict, os.path.join(self.working_dir, "ref.fa"))
        elif ref_fasta_path is not None:
            assert genome_file is None, "Either genome_file or ref_fasta_path should be provided, not both"
            # Build the reference sequence dictionary from file
            self.ref_seq_dict = extract_ref_seq_from_fasta(ref_fasta_path, region_dict.keys())
        else:
            raise ValueError("Either genome_file or ref_fasta_path should be provided")
        assert all(" " not in locus for locus in self.ref_seq_dict.keys()), "Region names cannot contain spaces"

    def load_bias_model(self, bias_model_fname: str):
        '''
        Load DddA bias model. We will retrieve bias using bias_model["forward"][seq_context] or bias_model["reverse"][seq_context]
        '''
        # Load DddA bias model
        with open(bias_model_fname, "rb") as f:
            self.bias_model = pickle.load(f)

    def get_bias_single_seq(self, seq, seq_len_threshold=5):
        '''
        Compute the sequence bias for a given sequence
        :param seq: String. Sequence context
        :param seq_len_threshold: Positions at the edge of the locus will have a shorter context (we don't do padding).
            For these edge positions we set a super high bias to prevent false positive footprints.
        :return: List of [C-to-T bias,  G-to-A bias]
        '''
        assert self.bias_model is not None, "Need to run load_bias_model() to load DddA bias model first"
    
        # If the sequence is too short, return a high bias
        if len(seq) < seq_len_threshold:
            return [100, 100]
    
        # Extract the center base being edited
        seq_radius = int((len(seq) - 1) / 2)
        center_base = seq[seq_radius]
    
        # A and T can't be edited so they have a bias of 0
        if center_base in ["A", "T", "N"]:
            return [0, 0]
    
        elif center_base == "C":
            return [self.bias_model["forward"][seq], 0]
    
        elif center_base == "G":
            return [0, self.bias_model["reverse"][seq]]

    def get_bias(self, context_radius=1):
        '''
        Pre-compute the sequence bias for each position in the reference sequences
        :param context_radius: Radius of the context to consider for each position. For example, if context_radius=2, we will consider the 2 bases to the left and the 2 bases to the right of the position being edited. Total length is then 2 * context_radius + 1
        :return None: the self.bias_dict attribute is updated. Keys are the loci and values are numpy arrays of shape
            (len(ref_seq), 2) containing the C-to-T and G-to-A biases for each position
        '''

        # Pre-compute sequence bias for each position
        logging.info("Computing sequence bias")
        self.bias_dict = {}
        for locus in self.ref_seq_dict.keys():
            ref_seq = self.ref_seq_dict[locus]
            len_seq = len(ref_seq)
            context = [ref_seq[max(0, i - context_radius):min(len_seq, i + context_radius + 1)] for i in range(len_seq)]
            bias = [self.get_bias_single_seq(c, 2 * context_radius + 1) for c in context]
            self.bias_dict[locus] = np.array(bias)

    def align_reads(self, del_size_threshold=5, start_gap_threshold=500, 
                    end_gap_threshold=500, mask_CpG=True):
        '''
        Align nanopore reads to the reference sequences and extract the editing events
        :param del_size_threshold: Only record deletions of at least this size
        :param start_gap_threshold: Remove reads with a gap at the start of the alignment larger than this threshold
        :param end_gap_threshold: Remove reads with a gap at the end of the alignment larger than this threshold
        :param mask_CpG: Whether to mask CpG positions from footprinting. Methylated CpG can confound DddA footprints.
        :return: None, the self.edit_dict and self.del_dict attributes are updated. Keys are the loci and values are lists
            of numpy arrays containing the editing events and deletions for each read. Additionally, the self.read_ids attribute
            is updated. Keys are the loci and values are lists of read IDs. Order of the reads is  the same in these attributes.
        '''

        self.start_gap_threshold = start_gap_threshold
        self.end_gap_threshold = end_gap_threshold

        # First examine whether the bias has been computed. If not, compute it
        if self.bias_dict is None:
            self.get_bias()

        datasets = pd.read_csv('data/datasets.tsv', sep='\t', dtype=str, keep_default_na=False)
        dataset_mask = datasets['fastq_file'] == self.fastq_file
        if (datasets.loc[dataset_mask, 'num_reads'] != '').any():
            total_reads = next(int(x) for x in datasets.loc[dataset_mask, 'num_reads'] if x != '')
        else:
            logging.info("Estimating number of reads to be processed")
            total_reads = count_reads_in_fastq(self.fastq_file)

        # Align reads to the reference sequences and extract editing events
        logging.info(f"Aligning reads in {self.fastq_file}")
        n_aligned, n_unaligned = 0, 0
        self.edit_dict = {locus: [] for locus in self.ref_seq_dict.keys()}
        self.del_dict = {locus: [] for locus in self.ref_seq_dict.keys()}
        self.read_ids = {locus: [] for locus in self.ref_seq_dict.keys()}
        self.read_strands = {locus: [] for locus in self.ref_seq_dict.keys()}
        self.ABE_edit_dict = {locus: [] for locus in self.ref_seq_dict.keys()}
        with mp.Pool(mp.cpu_count() + 2) as pool, tempfile.NamedTemporaryFile("w+t") as ref_fasta:
            ref_fasta.writelines([f">{locus}\n{seq}\n" for locus, seq in self.ref_seq_dict.items() if locus in self.region_dict])
            ref_fasta.flush()
            for result in pool.imap(align_helper, (
                    [read_id, read_seq, ref_fasta.name, start_gap_threshold, end_gap_threshold, del_size_threshold, self.ref_seq_dict, mask_CpG]
                    for read_id, read_seq, qual in tqdm(get_fastq_iter(self.fastq_file), total=total_reads)
                )):

                if result is None:
                    n_unaligned += 1
                    continue
                n_aligned += 1

                locus = result["locus"]
                self.edit_dict[locus].append(result["edit_vec"])
                self.del_dict[locus].append(result["del_vec"])
                self.read_ids[locus].append(result["read_id"])
                self.read_strands[locus].append(result["strand"])
                self.ABE_edit_dict[locus].append(result["ABE_edit_vec"])

        logging.info("{n_aligned} reads aligned. {n_unaligned} reads unaligned".format(
            n_aligned=n_aligned, n_unaligned=n_unaligned)
        )
        # reload the file so that updates from other parallel processes don't get overridden
        datasets = pd.read_csv('data/datasets.tsv', sep='\t', dtype=str, keep_default_na=False)
        datasets.loc[dataset_mask, 'num_reads'] = str(n_aligned + n_unaligned)
        datasets.to_csv('data/datasets.tsv', sep='\t', index=False)
        logging.info(f"Updated {dataset_mask.sum()} entries in datasets.tsv: {self.fastq_file} has {n_aligned + n_unaligned:,} reads")

        # Positions with high editing frequency are likely to be SNPs. Remove them from the bias_dict and edit_dict
        for locus in self.edit_dict.keys():
            if len(self.edit_dict[locus]) == 0:
                self.edit_dict[locus] = None
                self.del_dict[locus] = None
                self.read_strands[locus] = None
                self.edit_dict[locus] = None
            else:
                self.edit_dict[locus] = scipy.sparse.vstack(self.edit_dict[locus])
                self.del_dict[locus] = scipy.sparse.vstack(self.del_dict[locus])
                self.read_strands[locus] = np.array(self.read_strands[locus])                
                self.ABE_edit_dict[locus] = scipy.sparse.vstack(self.ABE_edit_dict[locus])
                    
                if mask_CpG:
                    # We already masked CpG positions in DddA edits during alignment. 
                    # Now we mask CpG positions in the bias vector
                    ref_seq = self.ref_seq_dict[locus]
                    len_seq = len(ref_seq)
                    # Extract 3-mer context
                    context = [ref_seq[max(0, i - 1):min(len_seq, i + 2)] for i in range(len_seq)] 
                    CpG_mask = np.array([i for i in range(len(context))\
                        if context[i][:2] == "CG" or context[i][1:] == "CG"])
                    self.bias_dict[locus][CpG_mask, :] = 0
                    
            logging.info(f"{locus}: n_reads, len_read are {np.shape(self.edit_dict[locus])}")

    def check_valid_alignment(self, start_gap_threshold: int | None = None, end_gap_threshold: int | None = None) -> bool:
        logging.info("Checking alignment data")
        for attr in ['edit_dict', 'del_dict', 'read_ids', 'read_strands']:
            if not hasattr(self, attr):
                logging.info(f"Missing {attr}")
                return False
            if missing_locus := next((locus for locus in self.region_dict.keys() if locus not in getattr(self, attr)), None):
                logging.info(f"Missing {attr} for locus {missing_locus}")
                return False
        if start_gap_threshold is not None and (old_start_thresh := getattr(self, 'start_gap_threshold', None)) != start_gap_threshold:
            logging.info(f"Start gap threshold ({old_start_thresh}) does not agree with desired setting ({start_gap_threshold})")
            return False
        if end_gap_threshold is not None and (old_end_thresh := getattr(self, 'end_gap_threshold', None)) != end_gap_threshold:
            logging.info(f"End gap threshold ({old_end_thresh}) does not agree with desired setting ({end_gap_threshold})")
            return False
        return True

    def get_footprints(self, selection_dict, footprint_radii, combine_sides="mean_pval"):
        '''
        Compute footprints
        :param selection_dict: Dictionary that specifies the reads to footprint for each locus. Keys are the loci and
            values are lists of read_ids
        :param footprint_radii: List or array of radii to compute footprints for
        :return footprint_dict: Dictionary of computed footprints. Keys are the loci and values are the footprints for the specified reads
        '''

        footprint_dict = {locus: {} for locus in self.ref_seq_dict.keys()}
        for locus in selection_dict.keys():
            print("Computing footprints for", locus)
            for read_id in tqdm(selection_dict[locus]):
                read_ind = self.read_ids[locus].index(read_id)
                read_edits = np.array(self.edit_dict[locus][read_ind, :].todense())[0, :]
                read_strand = self.read_strands[locus][read_ind]
                read_bias = self.bias_dict[locus][:, read_strand]
                del_mask = 1 - np.array(self.del_dict[locus][read_ind, :].todense())[0, :]
                read_bias = np.multiply(read_bias, del_mask) # For deleted positions,set bias to be zero
                footprint = [get_footprint_single_scale(read_edits, read_bias, 
                                                        fp_radius=r, combine_sides=combine_sides) for r in footprint_radii]
                footprint_dict[locus][read_id] = footprint
        return footprint_dict

    def dedup_reads(self, locus, read_ids: Iterable[str], threshold_factor=0.01) -> np.ndarray[str]:
        '''
        :param selection_dict: Dictionary that specifies the reads to dedup for each locus. Keys are the loci and
            values are lists of read_ids
        :param read_ids: IDs of reads to dedup. Must be a subset of self.read_ids[locus]
        :param threshold_factor: Used to calculate threshold of edit number differences for calling dup reads. Threshold = number_of_C_and_G * threshold_factor
        :return deduped_ids: IDs of reads passing deduplication
        '''
        logging.warn("Prefer `dedup_all` which saves results to object")
        deduped_inds = self.dedup_read_inds(locus=locus, read_ids=read_ids, threshold_factor=threshold_factor)
        return np.array(self.read_ids[locus])[deduped_inds]

    def dedup_read_inds(self, locus, read_ids: Iterable[str], threshold_factor=0.01, **kwargs) -> np.ndarray[int]:
        # For each read ID, get its index in the full read ID list
        locus_ids = self.read_ids[locus]
        locus_id_dict = dict(zip(locus_ids, np.arange(len(locus_ids))))
        read_inds = [*map(locus_id_dict.get, read_ids)]

        # Determine threshold for calling duplicate reads
        ref_seq = self.ref_seq_dict[locus]
        CG_mask = np.arange(len(ref_seq))[self.start_gap_threshold:-self.end_gap_threshold]
        C_pos = np.array([i for i in CG_mask if ref_seq[i] == "C"])
        G_pos = np.array([i for i in CG_mask if ref_seq[i] == "G"])
        threshold = int((len(C_pos) + len(G_pos)) * threshold_factor)
        self.threshold = threshold
        logging.info(f"Threshold for calling duplicate reads: {threshold}")

        return dedup_read_inds(self.edit_dict[locus], read_inds, threshold, **kwargs)

    def dedup_all(self, threshold_factor: float = 0.01, **kwargs):
        self.deduped_inds = {}
        for locus, read_ids in self.read_ids.items():
            self.deduped_inds[locus] = self.dedup_read_inds(
                locus=locus, read_ids=read_ids, threshold_factor=threshold_factor, **kwargs
            )
