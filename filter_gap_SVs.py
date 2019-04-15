#!/usr/bin/env python
import re

from intervaltree import IntervalTree

from ragoo_utilities.SeqReader import SeqReader


class BaseSequence(object):

    def __init__(self, in_sequence):
        if not isinstance(in_sequence, str):
            raise AttributeError('Only a string can be used to instantiate this class.')
        self.sequence = in_sequence.upper()


class GapSequence(BaseSequence):
    def get_gap_coords(self):
        """ Find all of the gap string indices for this sequence. """
        return re.finditer(r'N+', self.sequence)


def make_gaps_tree(in_file):
    # A dictionary to store an interval tree for each chromosome header.
    all_trees = dict()
    x = SeqReader(in_file)
    for header, sequence in x.parse_fasta():
        # Remove the greater than sign and only get first token if delimited by spaces
        header = header[1:].split(' ')[0]
        all_trees[header] = IntervalTree()
        gap_sequence = GapSequence(sequence)
        all_coordinates = [(m.start(0), m.end(0)) for m in gap_sequence.get_gap_coords()]
        for i in all_coordinates:
            all_trees[header][i[0]:i[1]] = i
    return all_trees


def get_query_bed_coords(in_line):
    L1 = in_line.split(':')
    header = L1[0]
    L2 = L1[1].split('-')
    start, stop = int(L2[0]), int(L2[1])
    return header, start, stop


def make_svs_bed(in_trees_q, in_trees_r):
    final_lines = []
    with open('assemblytics_out.Assemblytics_structural_variants.bed', 'r') as f:

        # Make a new header column for the field I am about to create
        header = f.readline().rstrip() + '\tquery_gap_ovlp'
        final_lines.append(header)

        # Calculate the percentage overlap for each structural variant.
        for line in f:
            pct_ovlp = 0
            ovlp = 0
            L1 = line.rstrip().split('\t')
            head, start, end = get_query_bed_coords(L1[9])
            query = in_trees_q[head][start:end]
            if query:
                for interval in query:
                    gap_start, gap_end= interval.begin, interval.end

                    # Calculate the amount these to intervals overlap
                    ovlp += max(0, min(end, gap_end) - max(start, gap_start))

                pct_ovlp = ovlp/(end-start)

            # Add the new value to the original line
            L1.append(pct_ovlp)
            L1 = [str(i) for i in L1]
            final_lines.append('\t'.join(L1))

    with open('assemblytics_out.Assemblytics_structural_variants.bed', 'w') as out_file:
        out_file.write('\n'.join(final_lines))

    # Now repeat but with respect to the reference
    final_lines = []
    with open('assemblytics_out.Assemblytics_structural_variants.bed', 'r') as f:
        # Make a new header column for the field I am about to create
        header = f.readline().rstrip() + '\tref_gap_ovlp'
        final_lines.append(header)

        # Calculate the percentage overlap for each structural variant.
        for line in f:
            pct_ovlp = 0
            ovlp = 0
            L1 = line.rstrip().split('\t')
            head, start, end = L1[0], int(L1[1]), int(L1[2])
            query = in_trees_r[head][start:end]
            if query:
                for interval in query:
                    gap_start, gap_end = interval.begin, interval.end

                    # Calculate the amount these to intervals overlap
                    ovlp += max(0, min(end, gap_end) - max(start, gap_start))

                pct_ovlp = ovlp / (end - start)

            # Add the new value to the original line
            L1.append(pct_ovlp)
            L1 = [str(i) for i in L1]
            final_lines.append('\t'.join(L1))

        with open('assemblytics_out.Assemblytics_structural_variants.bed', 'w') as out_file:
            out_file.write('\n'.join(final_lines))


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Annotate the SV bed file with the amount degree in which the SV overlaps with a gap')
    parser.add_argument("ref_file", metavar="<ref.fasta>", type=str,help="reference fasta file")
    args = parser.parse_args()
    ref_file = args.ref_file

    # Make a new bed file w.r.t the query
    all_trees_q = make_gaps_tree('../ragoo.fasta')
    all_trees_r = make_gaps_tree(ref_file)
    make_svs_bed(all_trees_q, all_trees_r)


if __name__ == "__main__":
    main()
