import operator
from collections import defaultdict

from ragoo_utilities.utilities import summarize_planesweep, binary_search


class ContigAlignment:
    """
    This object will represent a contig's alignment to a reference chromosome.
    Specifically, this class will organize aligments in a one-to-many fashion, where
    we will only have one object per contig, and each object will store all of the alignments
    of that contig.

    The ultimate goal is to uniquely assign every contig to a chromosome. In the case where the
    contig only aligns to one chromosome, the answer is already solved.

    If one contig maps to many chromosomes, we will have to find a way to find the 'best'
    chromosome for that contig given the nucmer alignment metrics. This object should store
    all such alignment metrics which might be helpful for this.
    """

    def __init__(self, in_contig):
        """
        Each object will refer to a single contig via its header name.

        All other alignment metrics will be stored in lists, where the list offsets
        correspond to the alignment number.

        Don't actually add any alignments through object instantiation. Only use the add_alignment method.
        """
        self.contig = in_contig

        self.query_lens = []
        self.query_starts = []
        self.query_ends = []
        self.strands = []
        self.ref_headers = []
        self.ref_lens = []
        self.ref_starts = []
        self.ref_ends = []
        self.num_matches = []
        self.aln_lens = []
        self.mapqs = []

    def __repr__(self):
        return '<ContigAlignment' + self.contig + '>'

    def __str__(self):
        """ Return a tab delimited string. Number of lines = number of alignments."""
        return ','.join(list(set(self.ref_headers)))

    def get_attr_lens(self):
        all_lens = [
            len(self.query_lens),
            len(self.query_starts),
            len(self.query_ends),
            len(self.strands),
            len(self.ref_headers),
            len(self.ref_lens),
            len(self.ref_starts),
            len(self.ref_ends),
            len(self.num_matches),
            len(self.aln_lens),
            len(self.mapqs)
        ]
        return all_lens

    def add_alignment(self, paf_line):
        """
        The only way to add alignments for this contig. Only accepts a PAFLine type object which represents
        as single line of a paf file.
        """
        if not paf_line.contig == self.contig:
            raise ValueError('Only can add alignments with query header %s' % self.contig)

        self.query_lens.append(paf_line.query_len)
        self.query_starts.append(paf_line.query_start)
        self.query_ends.append(paf_line.query_end)
        self.strands.append(paf_line.strand)
        self.ref_headers.append(paf_line.ref_header)
        self.ref_lens.append(paf_line.ref_len)
        self.ref_starts.append(paf_line.ref_start)
        self.ref_ends.append(paf_line.ref_end)
        self.num_matches.append(paf_line.num_match)
        self.aln_lens.append(paf_line.aln_len)
        self.mapqs.append(paf_line.mapq)

        # Check that all attribute lists have same length
        all_lens = self.get_attr_lens()
        assert len(set(all_lens)) == 1

    def has_unique_chr_match(self):
        s1 = set(self.ref_headers)
        if len(s1) == 1:
            return True
        return False

    def count_chr_matches(self):
        s1 = set(self.ref_headers)
        return len(s1)

    def filter_ref_chroms(self, in_chroms):
        """
        Given a list of chromosomes, mutate this object so as to only include alignments to those chromosomes.
        """
        hits = []
        for i in range(len(self.ref_headers)):
            if self.ref_headers[i] in in_chroms:
                hits.append(i)

        self.query_lens = [self.query_lens[i] for i in hits]
        self.query_starts = [self.query_starts[i] for i in hits]
        self.query_ends = [self.query_ends[i] for i in hits]
        self.strands = [self.strands[i] for i in hits]
        self.ref_headers = [self.ref_headers[i] for i in hits]
        self.ref_lens = [self.ref_lens[i] for i in hits]
        self.ref_starts = [self.ref_starts[i] for i in hits]
        self.ref_ends = [self.ref_ends[i] for i in hits]
        self.num_matches = [self.num_matches[i] for i in hits]
        self.aln_lens = [self.aln_lens[i] for i in hits]
        self.mapqs = [self.mapqs[i] for i in hits]

    def sort_by_ref(self):
        ref_pos = []
        for i in range(len(self.ref_headers)):
            ref_pos.append((self.ref_starts[i], self.ref_ends[i], i))
        hits = [i[2] for i in sorted(ref_pos)]

        self.query_lens = [self.query_lens[i] for i in hits]
        self.query_starts = [self.query_starts[i] for i in hits]
        self.query_ends = [self.query_ends[i] for i in hits]
        self.strands = [self.strands[i] for i in hits]
        self.ref_headers = [self.ref_headers[i] for i in hits]
        self.ref_lens = [self.ref_lens[i] for i in hits]
        self.ref_starts = [self.ref_starts[i] for i in hits]
        self.ref_ends = [self.ref_ends[i] for i in hits]
        self.num_matches = [self.num_matches[i] for i in hits]
        self.aln_lens = [self.aln_lens[i] for i in hits]
        self.mapqs = [self.mapqs[i] for i in hits]

    def exclude_ref_chroms(self, exclude_list):
        hits = []
        for i in range(len(self.ref_headers)):
            if self.ref_headers[i] not in exclude_list:
                hits.append(i)

        self.query_lens = [self.query_lens[i] for i in hits]
        self.query_starts = [self.query_starts[i] for i in hits]
        self.query_ends = [self.query_ends[i] for i in hits]
        self.strands = [self.strands[i] for i in hits]
        self.ref_headers = [self.ref_headers[i] for i in hits]
        self.ref_lens = [self.ref_lens[i] for i in hits]
        self.ref_starts = [self.ref_starts[i] for i in hits]
        self.ref_ends = [self.ref_ends[i] for i in hits]
        self.num_matches = [self.num_matches[i] for i in hits]
        self.aln_lens = [self.aln_lens[i] for i in hits]
        self.mapqs = [self.mapqs[i] for i in hits]

    def filter_lengths(self, l):
        hits = []
        for i in range(len(self.ref_headers)):
            if self.aln_lens[i] > l:
                hits.append(i)

        self.query_lens = [self.query_lens[i] for i in hits]
        self.query_starts = [self.query_starts[i] for i in hits]
        self.query_ends = [self.query_ends[i] for i in hits]
        self.strands = [self.strands[i] for i in hits]
        self.ref_headers = [self.ref_headers[i] for i in hits]
        self.ref_lens = [self.ref_lens[i] for i in hits]
        self.ref_starts = [self.ref_starts[i] for i in hits]
        self.ref_ends = [self.ref_ends[i] for i in hits]
        self.num_matches = [self.num_matches[i] for i in hits]
        self.aln_lens = [self.aln_lens[i] for i in hits]
        self.mapqs = [self.mapqs[i] for i in hits]

    def unique_anchor_filter(self):
        """
        The contents of this method are either influenced by or directly copied from "Assemblytics_uniq_anchor.py"
        written by Maria Nattestad. The original script can be found here:

        https://github.com/MariaNattestad/Assemblytics

        And the publication associated with Maria's work is here:

        Nattestad, Maria, and Michael C. Schatz. "Assemblytics: a
        web analytics tool for the detection of variants from an
        assembly." Bioinformatics 32.19 (2016): 3021-3023.
        """
        lines_by_query = []
        for i, j in zip(self.query_starts, self.query_ends):
            lines_by_query.append((i, j))

        hits = summarize_planesweep(lines_by_query, 10000)

        self.query_lens = [self.query_lens[i] for i in hits]
        self.query_starts = [self.query_starts[i] for i in hits]
        self.query_ends = [self.query_ends[i] for i in hits]
        self.strands = [self.strands[i] for i in hits]
        self.ref_headers = [self.ref_headers[i] for i in hits]
        self.ref_lens = [self.ref_lens[i] for i in hits]
        self.ref_starts = [self.ref_starts[i] for i in hits]
        self.ref_ends = [self.ref_ends[i] for i in hits]
        self.num_matches = [self.num_matches[i] for i in hits]
        self.aln_lens = [self.aln_lens[i] for i in hits]
        self.mapqs = [self.mapqs[i] for i in hits]


class UniqueContigAlignment:
    """
    A class representing the reference chromosome to which a particular contig will be assigned to.
    At first, the incoming alignments may be to multiple chromosomes, but this class will assign a single
    reference chromosome to an input contig.
    """

    def __init__(self, in_contig_aln):
        if not isinstance(in_contig_aln, ContigAlignment):
            message = 'Can only make a unique alignment from a ContigAlignment object. Got type %s instead.' % str(type(in_contig_aln))
            raise ValueError(message)

        self.contig = in_contig_aln.contig

        # Find the chromosome which was covered the most in these alignments.
        self.ref_chrom = None
        self.confidence = 0.0
        self._get_best_ref_cov(in_contig_aln)

    def __str__(self):
        return '<UniqueContigAlignment:' + self.contig + ' - ' + self.ref_chrom + '>'

    def _get_best_ref_cov(self, alns):
        # Get list of all references chromosomes
        all_chroms = set(alns.ref_headers)
        if len(all_chroms) == 1:
            self.ref_chrom = alns.ref_headers[0]
            self.confidence = 1
            return

        # Initialize coverage counts for each chromosome
        ranges = dict()
        for i in all_chroms:
            ranges[i] = 0

        # Get all the ranges in reference to all chromosomes
        all_intervals = defaultdict(list)
        for i in range(len(alns.ref_headers)):
            this_range = (alns.ref_starts[i], alns.ref_ends[i])
            this_chrom = alns.ref_headers[i]
            all_intervals[this_chrom].append(this_range)

        # For each chromosome, sort the range and get the union interval length.
        # Algorithm and pseudocode given by Michael Kirsche
        for i in all_intervals.keys():
            sorted_intervals = sorted(all_intervals[i], key=lambda tup: tup[0])
            max_end = -1
            for j in sorted_intervals:
                start_new_terr = max(j[0], max_end)
                ranges[i] += max(0, j[1] - start_new_terr)
                max_end = max(max_end, j[1])

        assert ranges
        max_chrom = max(ranges.items(), key=operator.itemgetter(1))[0]
        self.ref_chrom = max_chrom

        # Now get the confidence of this chromosome assignment
        # Equal to the max range over all ranges
        self.confidence = ranges[max_chrom]/sum(ranges.values())
        assert self.confidence >= 0
        assert self.confidence <= 1

class LongestContigAlignment:
    """
    Given a ContigAlignment object, find the longest alignment with respect to a given references chromosome.

    1. Filter out any alignments that are not aligned to specified chromosome.
    2. Find longest alignment remaining.
    """

    def __init__(self, in_contig_aln):
        self.best_index = self._get_best_index(in_contig_aln)
        self.contig = in_contig_aln.contig
        self.ref_start = in_contig_aln.ref_starts[self.best_index]
        self.ref_end = in_contig_aln.ref_ends[self.best_index]
        self.aln_len = in_contig_aln.aln_lens[self.best_index]
        self.strand = in_contig_aln.strands[self.best_index]
        self.interval = (self.ref_start, self.ref_end)

    def _get_best_index(self, aln):
        max_index = -1
        max_len = -1
        for i in range(len(aln.aln_lens)):
            if aln.aln_lens[i] > max_len:
                max_len = aln.aln_lens[i]
                max_index = i
        return max_index