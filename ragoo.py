#!/usr/bin/env python
from collections import defaultdict
from collections import OrderedDict
import copy

from intervaltree import IntervalTree

from ragoo_utilities.PAFReader import PAFReader
from ragoo_utilities.SeqReader import SeqReader
from ragoo_utilities.ReadCoverage import ReadCoverage
from ragoo_utilities.ContigAlignment import ContigAlignment
from ragoo_utilities.ContigAlignment import UniqueContigAlignment
from ragoo_utilities.ContigAlignment import LongestContigAlignment
from ragoo_utilities.GFFReader import GFFReader
from ragoo_utilities.utilities import run, log, reverse_complement, read_contigs, read_gz_contigs
from ragoo_utilities.break_chimera import get_ref_parts, cluster_contig_alns, avoid_gff_intervals, update_gff, break_contig, get_intra_contigs


def update_misasm_features(features, breaks, contig, ctg_len):

    # Get ctg len from ReadCoverage object
    break_list = [0] + sorted(breaks) + [ctg_len]
    borders = []
    for i in range(len(break_list) - 1):
        borders.append((break_list[i], break_list[i+1]))

    # Pop the features to be updated
    contig_feats = features.pop(contig)

    # Initialize lists for new broken contig headers
    for i in range(len(borders)):
        features[contig + '_misasm_break:' + str(borders[i][0]) + '-' + str(borders[i][1])] = []

    t = IntervalTree()
    for i in borders:
        t[i[0]:i[1]] = i

    for i in contig_feats:
        query = t[i.start]
        assert len(query) == 1
        break_start = list(query)[0].begin
        break_end = list(query)[0].end
        query_border = (break_start, break_end)
        break_number = borders.index(query_border)
        i.seqname = contig + '_misasm_break:' + str(borders[break_number][0]) + '-' + str(borders[break_number][1])
        i.start = i.start - break_start
        i.end = i.end - break_start
        features[
            contig + '_misasm_break:' + str(borders[break_number][0]) + '-' + str(borders[break_number][1])].append(i)

    return features


def remove_gff_breaks(gff_ins, breaks):
    """
    Given a list of candidate breakpoints proposed by misassembly correction, remove any such break points that
    fall within the interval of a gff feature. This should be called once per contig.
    :param gff_ins: List of GFFLines
    :param breaks: candidate break points
    :return:
    """
    # Make an interval tree from the intervals of the gff lines
    t = IntervalTree()
    for line in gff_ins:
        # If the interval is one bp long, skip
        if line.start == line.end:
            continue
        t[line.start:line.end] = (line.start, line.end)

    return [i for i in breaks if not t[i]]


def write_misasm_broken_ctgs(contigs_file, breaks, out_prefix, in_gff=None, in_gff_name=None):
    current_path = os.getcwd()
    os.chdir('ctg_alignments')

    if in_gff and in_gff_name:
        with open(in_gff_name, 'w') as f:
            for i in in_gff.keys():
                for j in in_gff[i]:
                    f.write(str(j) + '\n')

    x = SeqReader("../../" + contigs_file)
    f = open(out_prefix + ".misasm.break.fa", 'w')
    for header, seq in x.parse_fasta():
        header = header[1:]
        if header not in breaks:
            f.write(">" + header + "\n")
            f.write(seq + "\n")
        else:
            # Break the contig
            ctg_len = len(seq)
            break_list = [0] + sorted(breaks[header]) + [ctg_len]
            for i in range(len(break_list) - 1):
                f.write(">" + header + "_misasm_break:" + str(break_list[i]) + "-" + str(break_list[i+1]) + "\n")
                f.write(seq[break_list[i]:break_list[i+1]] + "\n")
    os.chdir(current_path)


def align_misasm_broken(out_prefix):
    current_path = os.getcwd()
    os.chdir('ctg_alignments')

    ctgs_file = out_prefix + ".misasm.break.fa"
    cmd = '{} -k19 -w19 -t{} ../../{}  {} ' \
          '> contigs_brk_against_ref.paf 2> contigs_brk_against_ref.paf.log'.format(minimap_path, t, reference_file,
                                                                            ctgs_file)
    if not os.path.isfile('contigs_brk_against_ref.paf'):
        run(cmd)
    os.chdir(current_path)


def write_contig_clusters(unique_dict, thresh, skip_list):
    # Get a list of all chromosomes
    all_chroms = set([unique_dict[i].ref_chrom for i in unique_dict.keys()])
    current_path = os.getcwd()
    output_path = current_path + '/groupings'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    os.chdir('groupings')
    for i in all_chroms:
        open(i + '_contigs.txt', 'w').close()

    for i in unique_dict.keys():
        this_chr = unique_dict[i].ref_chrom
        this_confidence = unique_dict[i].confidence
        if this_confidence > thresh:
            if not i in skip_list:
                file_name = str(this_chr) + '_contigs.txt'
                with open(file_name, 'a') as f:
                    f.write(i + '\t' + str(this_confidence) + '\n')
    os.chdir(current_path)


def clean_alignments(in_alns, l=10000, in_exclude_file='', uniq_anchor_filter=False, merge=False, quality=0):
    # Exclude alignments to undesired reference headers and filter alignment lengths.
    exclude_list = []
    if in_exclude_file:
        with open('../' + in_exclude_file) as f:
            for line in f:
                exclude_list.append(line.rstrip().replace('>', '').split()[0])

    empty_headers = []
    for header in in_alns.keys():
        in_alns[header].exclude_ref_chroms(exclude_list)
        in_alns[header].filter_lengths(l)
        if uniq_anchor_filter:
            in_alns[header].unique_anchor_filter()
        if quality:
            in_alns[header].filter_quality(quality)

        if merge:
            in_alns[header].merge_alns()

        # Check if our filtering has removed all alignments for a contig
        if len(in_alns[header].ref_headers) == 0:
            empty_headers.append(header)

    for header in empty_headers:
        in_alns.pop(header)
    return in_alns


def read_paf_alignments(in_paf):
    # Read in PAF alignments
    # Initialize a dictionary where key is contig header, and value is ContigAlignment.
    alns = dict()
    x = PAFReader(in_paf)
    for paf_line in x.parse_paf():
        if paf_line.contig in alns:
            alns[paf_line.contig].add_alignment(paf_line)
        else:
            alns[paf_line.contig] = ContigAlignment(paf_line.contig)
            alns[paf_line.contig].add_alignment(paf_line)
    return alns


def get_contigs_from_groupings(in_file):
    contigs = []
    with open(in_file) as f:
        for line in f:
            contigs.append(line.split('\t')[0])
    return contigs


def get_location_confidence(in_ctg_alns):
    # Use interval tree to get all alignments with the reference span
    # Go through each of them and if any start is less than the min_pos or any end is greater than
    # the max_pos, change the borders to those values. Then use the algorithm that Mike gave me.
    min_pos = min(in_ctg_alns.ref_starts)
    max_pos = max(in_ctg_alns.ref_ends)
    t = IntervalTree()

    # Put the reference start and end position for every alignment into the tree
    for i in range(len(in_ctg_alns.ref_headers)):
        t[in_ctg_alns.ref_starts[i]:in_ctg_alns.ref_ends[i]] = (in_ctg_alns.ref_starts[i], in_ctg_alns.ref_ends[i])

    overlaps = t[min_pos:max_pos]
    if not overlaps:
        return 0

    # If any intervals fall beyond the boundaries, replace the start/end with the boundary it exceeds
    ovlp_list = [i.data for i in overlaps]
    bounded_list = []
    for i in ovlp_list:
        if i[0] < min_pos:
            i[0] = min_pos
        if i[1] > max_pos:
            i[1] = max_pos
        bounded_list.append(i)

    # Now can just calculate the total range covered by the intervals
    ovlp_range = 0
    sorted_intervals = sorted(bounded_list, key=lambda tup: tup[0])
    max_end = -1
    for j in sorted_intervals:
        start_new_terr = max(j[0], max_end)
        ovlp_range += max(0, j[1] - start_new_terr)
        max_end = max(max_end, j[1])

    return ovlp_range / (max_pos - min_pos)


def order_orient_contigs(in_unique_contigs, in_alns):
    current_path = os.getcwd()
    output_path = current_path + '/orderings'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Get longest alignments
    longest_contigs = dict()
    for i in in_alns.keys():
        # Only consider alignments to the assigned chromosome
        uniq_aln = UniqueContigAlignment(in_alns[i])
        best_header = uniq_aln.ref_chrom
        ctg_alns = copy.deepcopy(in_alns[i])
        ctg_alns.filter_ref_chroms([best_header])
        longest_contigs[i] = LongestContigAlignment(ctg_alns)

    # Save the orientations
    final_orientations = dict()
    for i in longest_contigs.keys():
        final_orientations[i] = longest_contigs[i].strand

    # Get the location and orientation confidence scores
    orientation_confidence = dict()
    location_confidence = dict()
    forward_bp = 0
    reverse_bp = 0
    for i in in_alns.keys():
        uniq_aln = UniqueContigAlignment(in_alns[i])
        best_header = uniq_aln.ref_chrom
        ctg_alns = copy.deepcopy(in_alns[i])
        ctg_alns.filter_ref_chroms([best_header])

        # Orientation confidence scores
        # Every base pair votes for the orientation of the alignment in which it belongs
        # Score is # votes for the assigned orientation over all votes
        for j in range(len(ctg_alns.ref_headers)):
            if ctg_alns.strands[j] == '+':
                forward_bp += ctg_alns.aln_lens[j]
            else:
                reverse_bp += ctg_alns.aln_lens[j]

        if final_orientations[i] == '+':
            orientation_confidence[i] = forward_bp / (forward_bp + reverse_bp)
        else:
            orientation_confidence[i] = reverse_bp / (forward_bp + reverse_bp)

        forward_bp = 0
        reverse_bp = 0

        # Location confidence
        location_confidence[i] = get_location_confidence(ctg_alns)

    all_chroms = set([in_unique_contigs[i].ref_chrom for i in in_unique_contigs.keys()])

    for this_chrom in all_chroms:

        # Intialize the list of start and end positions w.r.t the query
        ref_pos = []

        groupings_file = 'groupings/' + this_chrom + '_contigs.txt'
        contigs_list = get_contigs_from_groupings(groupings_file)

        for i in range(len(contigs_list)):
            # There is a scope issue here. Pass this (longest_contigs) to the method explicitly.
            ref_pos.append((longest_contigs[contigs_list[i]].ref_start, longest_contigs[contigs_list[i]].ref_end, i))

        final_order = [contigs_list[i[2]] for i in sorted(ref_pos)]

        # Get ordering confidence
        # To do this, get the max and min alignments to this reference chromosome
        # Then within that region, what percent of bp are covered

        with open('orderings/' + this_chrom + '_orderings.txt', 'w') as out_file:
            for i in final_order:
                # Also have a scope issue here.
                out_file.write(i + '\t' + final_orientations[i] + '\t' + str(location_confidence[i]) + '\t' + str(orientation_confidence[i]) + '\n')


def get_orderings(in_orderings_file):
    all_orderings = []
    with open(in_orderings_file) as f:
        for line in f:
            L1 = line.split('\t')
            all_orderings.append((L1[0], L1[1].rstrip()))
    return all_orderings


def create_pseudomolecules(in_contigs_file, in_unique_contigs, gap_size, chr0=True):
    """
    Need to make a translation table for easy lift-over.
    :param in_contigs_file:
    :param in_unique_contigs:
    :param gap_size:
    :return:
    """
    # First, read all of the contigs into memory
    remaining_contig_headers = []
    all_seqs = OrderedDict()
    x = SeqReader('../' + in_contigs_file)
    if in_contigs_file.endswith(".gz"):
        for header, seq in x.parse_gzip_fasta():
            remaining_contig_headers.append(header.split(' ')[0])
            all_seqs[header.split(' ')[0]] = seq
    else:
        for header, seq in x.parse_fasta():
            remaining_contig_headers.append(header.split(' ')[0])
            all_seqs[header.split(' ')[0]] = seq

    # Get all reference chromosomes
    all_chroms = sorted(list(set([in_unique_contigs[i].ref_chrom for i in in_unique_contigs.keys()])))

    # Iterate through each orderings file and store sequence in a dictionary
    all_pms = dict()
    pad = ''.join('N' for i in range(gap_size))
    for this_chrom in all_chroms:
        orderings_file = 'orderings/' + this_chrom + '_orderings.txt'
        orderings = get_orderings(orderings_file)
        if orderings:
            seq_list = []
            for line in orderings:
                # Mark that we have seen this contig
                remaining_contig_headers.pop(remaining_contig_headers.index('>' + line[0]))
                if line[1] == '+':
                    seq_list.append(all_seqs['>' + line[0]])
                else:
                    assert line[1] == '-'
                    seq_list.append(reverse_complement(all_seqs['>' + line[0]]))
            all_pms[this_chrom] = pad.join(seq_list)
            all_pms[this_chrom] += '\n'

    # Get unincorporated sequences and place them in Chr0
    if remaining_contig_headers:
        if chr0:
            chr0_headers = []
            chr0_seq_list = []
            for header in remaining_contig_headers:
                chr0_headers.append(header)
                chr0_seq_list.append(all_seqs[header])
            all_pms['Chr0'] = pad.join(chr0_seq_list)
            all_pms['Chr0'] += '\n'

            # Write out the list of chr0 headers
            f_chr0_g = open('groupings/Chr0_contigs.txt', 'w')
            f_chr0_o = open('orderings/Chr0_orderings.txt', 'w')
            for i in chr0_headers:
                f_chr0_g.write(i[1:] + "\t" + "0" + '\n')
                f_chr0_o.write(i[1:] + '\t' + "+" + '\t' + "0" + '\t' + "0" + '\n')
            f_chr0_g.close()
            f_chr0_o.close()
        else:
            # Instead of making a chromosome 0, add the unplaced sequences as is.
            for header in remaining_contig_headers:
                all_pms[header[1:]] = all_seqs[header] + "\n"
                f_chr0_g = open('groupings/' + header[1:] + '_contigs.txt', 'w')
                f_chr0_o = open('orderings/' + header[1:] + '_orderings.txt', 'w')
                f_chr0_g.write(header[1:] + "\t" + "0" + '\n')
                f_chr0_o.write(header[1:] + '\t' + "+" + '\t' + "0" + '\t' + "0" + '\n')
                f_chr0_g.close()
                f_chr0_o.close()

    # Write the final sequences out to a file
    with open('ragoo.fasta', 'w') as f:
        for out_header in all_pms:
            f.write(">" + out_header + "_RaGOO\n")
            f.write(all_pms[out_header])


def write_broken_files(in_contigs, in_contigs_name, in_gff=None, in_gff_name=None):
    current_path = os.getcwd()
    output_path = current_path + '/chimera_break'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    os.chdir('chimera_break')
    if in_gff and in_gff_name:
        with open(in_gff_name, 'w') as f:
            for i in in_gff.keys():
                for j in in_gff[i]:
                    f.write(str(j) + '\n')

    with open(in_contigs_name, 'w') as f:
        for i in in_contigs.keys():
            f.write('>' + i + '\n')
            f.write(in_contigs[i] + '\n')

    os.chdir(current_path)


def align_breaks(break_type, m_path, in_reference_file, in_contigs_file, in_num_threads):
    current_path = os.getcwd()
    os.chdir('chimera_break')
    if break_type == 'inter':
        cmd = '{} -k19 -w19 -t{} ../../{} {} ' \
          '> inter_contigs_against_ref.paf 2> inter_contigs_against_ref.paf.log'.format(m_path, in_num_threads, in_reference_file, in_contigs_file)
        if not os.path.isfile('inter_contigs_against_ref.paf'):
            run(cmd)
    else:
        cmd = '{} -k19 -w19 -t{} ../../{} {} ' \
              '> intra_contigs_against_ref.paf 2> intra_contigs_against_ref.paf.log'.format(m_path, in_num_threads, in_reference_file, in_contigs_file)
        if not os.path.isfile('intra_contigs_against_ref.paf'):
            run(cmd)

    os.chdir(current_path)


def align_pms(m_path, num_threads, in_reference_file):
    current_path = os.getcwd()
    output_path = current_path + '/pm_alignments'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir('pm_alignments')

    cmd = '{} -ax asm5 --cs -t{} -I20G  ../../{} {} ' \
          '> pm_against_ref.sam 2> pm_contigs_against_ref.sam.log'.format(m_path, num_threads,
                                                                                        in_reference_file, '../ragoo.fasta')
    if not os.path.isfile('pm_against_ref.sam'):
        run(cmd)

    os.chdir(current_path)


def get_SVs(sv_min, sv_max, in_ref_file):
    current_path = os.getcwd()
    os.chdir('pm_alignments')
    # Change this when setup.py is ready. Just call script directly
    cmd = 'sam2delta.py pm_against_ref.sam'
    if not os.path.isfile('pm_against_ref.sam.delta'):
        run(cmd)

    cmd_2 = 'Assemblytics_uniq_anchor.py --delta pm_against_ref.sam.delta --unique-length 10000 --out assemblytics_out --keep-small-uniques'
    if not os.path.isfile('assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta'):
        run(cmd_2)

    cmd_3 = 'Assemblytics_between_alignments.pl assemblytics_out.coords.tab %r %r all-chromosomes exclude-longrange bed > assemblytics_out.variants_between_alignments.bed' %(sv_min, sv_max)
    if not os.path.isfile('assemblytics_out.variants_between_alignments.bed'):
        run(cmd_3)

    cmd_4 = 'Assemblytics_within_alignment.py --delta assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta --min %r > assemblytics_out.variants_within_alignments.bed' %(sv_min)
    if not os.path.isfile('assemblytics_out.variants_within_alignments.bed'):
        run(cmd_4)

    header = "reference\tref_start\tref_stop\tID\tsize\tstrand\ttype\tref_gap_size\tquery_gap_size\tquery_coordinates\tmethod\n"

    with open('assemblytics_out.variants_between_alignments.bed', 'r')as f1:
        b1 = f1.read()

    with open('assemblytics_out.variants_within_alignments.bed', 'r') as f2:
        b2 = f2.read()

    with open('assemblytics_out.Assemblytics_structural_variants.bed', 'w') as f:
        f.write(header)
        # Might need to add newlines here
        f.write(b1)
        f.write(b2)

    # Filter out SVs caused by gaps
    cmd_5 = 'filter_gap_SVs.py ../../%s' %(in_ref_file)
    run(cmd_5)

    os.chdir(current_path)


def align_reads(m_path, num_threads, in_ctg_file, reads, tech='ont'):
    current_path = os.getcwd()
    output_path = current_path + '/ctg_alignments'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir('ctg_alignments')

    if tech == 'sr':
        cmd = '{} -x sr -t{} ../../{} ../../{} ' \
              '> reads_against_ctg.paf 2> reads_against_ctg.paf.log'.format(m_path, num_threads, in_ctg_file, reads)
    elif tech == 'corr':
        cmd = '{} -x asm10 -t{} ../../{} ../../{} ' \
              '> reads_against_ctg.paf 2> reads_against_ctg.paf.log'.format(m_path, num_threads, in_ctg_file, reads)
    else:
        raise ValueError("Only 'sr' or 'corr' are accepted for read type.")

    if not os.path.isfile('reads_against_ctg.paf'):
        run(cmd)

    os.chdir(current_path)


if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description='order and orient contigs according to minimap2 alignments to a reference (v1.1)')
    parser.add_argument("contigs", metavar="<contigs.fasta>", type=str, help="fasta file with contigs to be ordered and oriented (gzipped allowed)")
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file (gzipped allowed)")
    #parser.add_argument("-o", metavar="PATH", type=str, default="ragoo_output", help="output directory name")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-gff", metavar="<annotations.gff>", type=str, default='', help="lift-over gff features to chimera-broken contigs")
    parser.add_argument("-m", metavar="PATH", type=str, default="minimap2", help='path to minimap2 executable')
    parser.add_argument("-b", action='store_true', default=False, help="Break chimeric contigs")
    parser.add_argument("-R", metavar="<reads.fasta>", type=str, default="", help="Turns on misassembly correction. Align provided reads to the contigs to aid misassembly correction. fastq or fasta allowed. Gzipped files allowed. Turns off '-b'.")
    parser.add_argument("-T", metavar="sr", type=str, default="", help="Type of reads provided by '-R'. 'sr' and 'corr' accepted for short reads and error corrected long reads respectively.")
    parser.add_argument("-I", type=str, metavar="index_split", default="4G", help="Minimap2: split index for every ~NUM input bases [4G]")
    parser.add_argument("--mini-extra", type=str, help="Extra flags to pass-through to MiniMap2.")
    parser.add_argument("-p", metavar="5", type=int, default=5, help=argparse.SUPPRESS)
    parser.add_argument("-l", metavar="10000", type=int, default=10000, help=argparse.SUPPRESS)
    parser.add_argument("-r", metavar="100000", type=int, default=100000, help=argparse.SUPPRESS)
    parser.add_argument("-c", metavar="1000000", type=int, default=1000000, help=argparse.SUPPRESS)
    parser.add_argument("-d", metavar="2000000", type=int, default=2000000, help=argparse.SUPPRESS)
    parser.add_argument("-t", metavar="3", type=int, default=3, help="Number of threads when running minimap.")
    parser.add_argument("-g", metavar="100", type=int, default=100, help="Gap size for padding in pseudomolecules.")
    parser.add_argument("-s", action='store_true', default=False, help="Call structural variants")
    parser.add_argument("-a", metavar="50", type=int, default=50, help=argparse.SUPPRESS)
    parser.add_argument("-f", metavar="10000", type=int, default=10000, help=argparse.SUPPRESS)
    parser.add_argument("-i", metavar="0.2", type=float, default=0.2, help="Minimum grouping confidence score needed to be localized.")
    parser.add_argument("-j", metavar="<skip.txt>", type=str, default="", help="List of contigs to automatically put in chr0.")
    parser.add_argument("-C", action='store_true', default=False, help="Write unplaced contigs individually instead of making a chr0")
    parser.add_argument("-q", "--quality", default=0, type=int, metavar="min. mapping quality",
                        help="Minimum mapping quality for Minimap2 alignments. Default: 0 (no filtering)")

    # Get the command line arguments
    args = parser.parse_args()
    contigs_file = args.contigs
    reference_file = args.reference
    #output_path = args.o
    exclude_file = args.e
    minimap_path = args.m
    break_chimeras = args.b
    gff_file = args.gff
    min_break_pct = args.p
    min_len = args.l
    min_range = args.r
    intra_wrt_ref_min = args.d
    intra_wrt_ctg_min = args.c
    t = args.t
    g = args.g
    call_svs = args.s
    a = args.a
    f = args.f
    group_score_thresh = args.i
    skip_file = args.j
    corr_reads = args.R
    corr_reads_tech = args.T
    make_chr0 = not args.C

    if corr_reads:
        log("Misassembly correction has been turned on. This automatically inactivates chimeric contig correction.")
        break_chimeras = False

    # Make sure that if -R, -T has been specified
    if corr_reads and not corr_reads_tech:
        raise ValueError("'-T' must be provided when using -R.")

    skip_ctg = []
    if skip_file:
        with open(skip_file) as f:
            for line in f:
                skip_ctg.append(line.rstrip())

    current_path = os.getcwd()
    output_path = current_path + '/ragoo_output'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir(output_path)

    # Run minimap2
    cmd = '{} -k19 -I {} {} -w19 -t{} ../{} ../{} ' \
          '> contigs_against_ref.paf 2> contigs_against_ref.paf.log'.format(minimap_path, args.I, args.mini_extra,
                                                                            t, reference_file, contigs_file)

    if not os.path.isfile('contigs_against_ref.paf'):
        run(cmd)

    # Read in the minimap2 alignments just generated
    log('Reading alignments')
    alns = read_paf_alignments('contigs_against_ref.paf')
    alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file, quality=args.quality)

    # Process the gff file
    if gff_file:
        log('Getting gff features')
        features = defaultdict(list)
        z = GFFReader('../' + gff_file)
        for i in z.parse_gff():
            features[i.seqname].append(i)

    # Break chimeras if desired
    if break_chimeras:
        # Record how many contigs are broken
        total_inter_broken = 0
        total_intra_broken = 0

        alns = clean_alignments(alns, l=10000, in_exclude_file=exclude_file, uniq_anchor_filter=True)
        # Process contigs
        log('Getting contigs')
        if contigs_file.endswith(".gz"):
            contigs_dict = read_gz_contigs('../' + contigs_file)
        else:
            contigs_dict = read_contigs('../' + contigs_file)

        log('Finding interchromosomally chimeric contigs')
        all_chimeras = dict()
        for i in alns.keys():
            ref_parts = get_ref_parts(alns[i], min_len, min_break_pct, min_range)
            if len(ref_parts) > 1:
                all_chimeras[i] = ref_parts

        log('Finding break points and breaking interchromosomally chimeric contigs')
        break_intervals = dict()
        for i in all_chimeras.keys():
            break_intervals[i] = cluster_contig_alns(i, alns, all_chimeras[i], min_len)

            # If its just going to break it into the same thing, skip it.
            if len(break_intervals[i]) <= 1:
                continue

            if gff_file:
                # If desired, ensure that breakpoints don't disrupt any gff intervals
                break_intervals[i] = avoid_gff_intervals(break_intervals[i], features[i])
                features = update_gff(features, break_intervals[i], i)

            # Break contigs according to the final break points
            contigs_dict = break_contig(contigs_dict, i, break_intervals[i])
            total_inter_broken += 1

        # Next, need to re-align before finding intrachromosomal chimeras
        # First, write out the interchromosomal chimera broken fasta
        out_inter_fasta = contigs_file[:contigs_file.rfind('.')] + '.inter.chimera.broken.fa'
        if gff_file:
            out_gff = gff_file[:gff_file.rfind('.')] + '.inter.chimera_broken.gff'
            write_broken_files(contigs_dict, out_inter_fasta, features, out_gff)
        else:
            write_broken_files(contigs_dict, out_inter_fasta)

        # Next, realign the chimera broken contigs
        align_breaks('inter', minimap_path, reference_file, out_inter_fasta, t)

        # Now, use those new alignments for intrachromosomal chimeras
        log('Reading interchromosomal chimera broken alignments')
        inter_alns = read_paf_alignments('chimera_break/inter_contigs_against_ref.paf')
        inter_alns = clean_alignments(inter_alns, l=1000, in_exclude_file=exclude_file)

        log('Finding intrachromosomally chimeric contigs')
        # Find intrachromosomally chimeric contigs
        for i in inter_alns.keys():
            intra = get_intra_contigs(inter_alns[i], 15000, intra_wrt_ref_min, intra_wrt_ctg_min)
            if intra:
                if gff_file:
                    intra_break_intervals = avoid_gff_intervals(intra[1], features[intra[0]])
                else:
                    intra_break_intervals = intra[1]
                # Check if the avoidance of gff intervals pushed the break point to the end of the contig.
                if intra_break_intervals[-1][0] == intra_break_intervals[-1][1]:
                    continue

                # break the contigs and update features if desired
                contigs_dict = break_contig(contigs_dict, intra[0], intra_break_intervals)
                total_intra_broken += 1

                if gff_file:
                    features = update_gff(features, intra_break_intervals, intra[0])

        # Write out the intrachromosomal information
        out_intra_fasta = contigs_file[:contigs_file.rfind('.')] + '.intra.chimera.broken.fa'
        if gff_file:
            out_intra_gff = gff_file[:gff_file.rfind('.')] + '.intra.chimera_broken.gff'
            write_broken_files(contigs_dict, out_intra_fasta, features, out_intra_gff)
        else:
            write_broken_files(contigs_dict, out_intra_fasta)

        # Re align the contigs
        # Next, realign the chimera broken contigs
        align_breaks('intra', minimap_path, reference_file, out_intra_fasta, t)

        # Read in alignments of intrachromosomal chimeras and proceed with ordering and orientation
        log('Reading intrachromosomal chimera broken alignments')
        alns = read_paf_alignments('chimera_break/intra_contigs_against_ref.paf')
        alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file)
        contigs_file = '/ragoo_output/chimera_break/' + out_intra_fasta
        log('The total number of interchromasomally chimeric contigs broken is %r' % total_inter_broken)
        log('The total number of intrachromasomally chimeric contigs broken is %r' % total_intra_broken)

    # Check if misassembly correction is turned on. This is mutually exclusive with chimeric contig correction
    if corr_reads:
        # Align the raw reads to the assembly.
        log('Aligning raw reads to contigs')
        align_reads(minimap_path, t, contigs_file, corr_reads, corr_reads_tech)
        log('Computing contig coverage')
        cov_map = ReadCoverage('ctg_alignments/reads_against_ctg.paf')
        alns = clean_alignments(alns, l=10000, in_exclude_file=exclude_file, uniq_anchor_filter=True, merge=True)

        # Get the initial candidate break points.
        candidate_breaks = dict()
        for i in alns:
            candidates = alns[i].get_break_candidates()
            if candidates:
                candidate_breaks[i] = candidates

        # Validate each breakpoint by checking for excessively high or low coverage
        # Also, if a gff is provided, check to ensure that we don't break within a gff feature interval
        val_candidate_breaks = dict()
        for i in candidate_breaks:
            candidates = cov_map.check_break_cov(i, candidate_breaks[i])
            if gff_file:
                candidates = remove_gff_breaks(features[i], candidates)
            if candidates:
                val_candidate_breaks[i] = list(set(candidates))
                if gff_file:
                    features = update_misasm_features(features, val_candidate_breaks[i], i, cov_map.ctg_lens[i])

        # Break the contigs
        if gff_file:
            out_misasm_gff = gff_file[:gff_file.rfind('.')] + '.misasm.broken.gff'
            write_misasm_broken_ctgs(contigs_file, val_candidate_breaks, contigs_file[:contigs_file.rfind('.')], in_gff=features, in_gff_name=out_misasm_gff)
        else:
            write_misasm_broken_ctgs(contigs_file, val_candidate_breaks, contigs_file[:contigs_file.rfind('.')])

        # Align the broken contigs back to the reference
        align_misasm_broken(contigs_file[:contigs_file.rfind('.')])
        alns = read_paf_alignments('ctg_alignments/contigs_brk_against_ref.paf')
        alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file)
        contigs_file = '/ragoo_output/ctg_alignments/' + contigs_file[:contigs_file.rfind('.')] + ".misasm.break.fa"

    # Assign each contig to a corresponding reference chromosome.
    log('Assigning contigs')
    all_unique_contigs = dict()
    for i in alns.keys():
        all_unique_contigs[i] = UniqueContigAlignment(alns[i])

    # Add to this the list of headers that did not make it
    write_contig_clusters(all_unique_contigs, group_score_thresh, skip_ctg)

    log('Ordering and orienting contigs')
    order_orient_contigs(all_unique_contigs, alns)

    log('Creating pseudomolecules')
    create_pseudomolecules(contigs_file, all_unique_contigs, g, make_chr0)

    if call_svs:
        log('Aligning pseudomolecules to reference')
        align_pms(minimap_path, t, reference_file)

        log('Getting structural variants')
        get_SVs(a, f, reference_file)

    log('goodbye')
