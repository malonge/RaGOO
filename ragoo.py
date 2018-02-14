from collections import defaultdict

from utilities.PAFReader import PAFReader
from utilities.SeqReader import SeqReader
from utilities.ContigAlignment import ContigAlignment
from utilities.ContigAlignment import UniqueContigAlignment
from utilities.ContigAlignment import LongestContigAlignment
from utilities.GFFReader import GFFReader
from utilities.utilities import run, log, reverse_complement, read_contigs
from utilities.break_chimera import get_ref_parts, cluster_contig_alns, avoid_gff_intervals, update_gff, break_contig, get_intra_contigs
import shutil


def write_contig_clusters(unique_dict):
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
        file_name = str(this_chr) + '_contigs.txt'
        with open(file_name, 'a') as f:
            f.write(i + '\n')
    os.chdir(current_path)


def clean_alignments(in_alns, l=10000, in_exclude_file='', uniq_anchor_filter=False):
    # Exclude alignments to undesired reference headers and filter alignment lengths.
    exclude_list = []
    if in_exclude_file:
        with open('../' + in_exclude_file) as f:
            for line in f:
                exclude_list.append(line.rstrip())

    empty_headers = []
    for header in in_alns.keys():
        in_alns[header].exclude_ref_chroms(exclude_list)
        in_alns[header].filter_lengths(l)
        if uniq_anchor_filter:
            in_alns[header].unique_anchor_filter()
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


def get_contigs(in_file):
    contigs = []
    with open(in_file) as f:
        for line in f:
            contigs.append(line.rstrip())
    return contigs


def order_contigs(in_unique_contigs):
    current_path = os.getcwd()
    output_path = current_path + '/orderings'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    all_chroms = set([in_unique_contigs[i].ref_chrom for i in in_unique_contigs.keys()])

    # For chromosome 0, this reports all contigs in alignments not placed.
    # But it still does not report all contigs not even in alignments in the first place.
    for this_chrom in all_chroms:

        # Intialize the list of start and end positions w.r.t the query
        ref_pos = []

        groupings_file = 'groupings/' + this_chrom + '_contigs.txt'
        contigs_list = get_contigs(groupings_file)

        # Add edges for each of n choose 2 pairs of alignments
        for i in range(len(contigs_list)):
            ref_pos.append((longest_contigs[contigs_list[i]].ref_start, longest_contigs[contigs_list[i]].ref_end, i))

        final_order = [contigs_list[i[2]] for i in sorted(ref_pos)]

        with open('orderings/' + this_chrom + '_orderings.txt', 'w') as out_file:
            for i in final_order:
                out_file.write(i + '\t' + final_orientations[i] + '\n')


def get_orderings(in_orderings_file):
    all_orderings = []
    with open(in_orderings_file) as f:
        for line in f:
            L1 = line.split('\t')
            all_orderings.append((L1[0], L1[1].rstrip()))
    return all_orderings


def create_pseudomolecules(in_contigs_file, in_unique_contigs, gap_size):
    """
    Need to make a translation table for easy lift-over.
    :param in_contigs_file:
    :param in_unique_contigs:
    :param gap_size:
    :return:
    """
    # First, read all of the contigs into memory
    remaining_contig_headers = []
    all_seqs = dict()
    x = SeqReader('../' + in_contigs_file)
    for header, seq in x.parse_fasta():
        remaining_contig_headers.append(header.split(' ')[0])
        all_seqs[header.split(' ')[0]] = seq

    # Get all reference chromosomes
    all_chroms = sorted(list(set([in_unique_contigs[i].ref_chrom for i in in_unique_contigs.keys()])))

    # Iterate through each orderings file and store sequence in a dictionary
    all_pms = dict()
    for this_chrom in all_chroms:
        all_pms[this_chrom] = ''
        orderings_file = 'orderings/' + this_chrom + '_orderings.txt'
        orderings = get_orderings(orderings_file)
        for line in orderings:
            # Mark that we have seen this contig
            remaining_contig_headers.pop(remaining_contig_headers.index('>' + line[0]))
            if line[1] == '+':
                all_pms[this_chrom] += all_seqs['>' + line[0]]
                all_pms[this_chrom] += ''.join('N' for i in range(gap_size))
            else:
                assert line[1] == '-'
                all_pms[this_chrom] += reverse_complement(all_seqs['>' + line[0]])
                all_pms[this_chrom] += ''.join('N' for i in range(gap_size))
        all_pms[this_chrom] += '\n'

    # Get unincorporated sequences and place them in Chr0
    all_pms['Chr0'] = ''
    for header in remaining_contig_headers:
        all_pms['Chr0'] += all_seqs[header]
        all_pms['Chr0'] += ''.join('N' for i in range(gap_size))
    all_pms['Chr0'] += '\n'

    # Write the final sequences out to a file
    with open('ragoo.fasta', 'w') as f:
        f.write('>Chr0_RaGOO\n')
        f.write(all_pms['Chr0'])
        for header in all_chroms:
            f.write('>' + header + '_RaGOO\n')
            f.write(all_pms[header])


def write_broken_files(in_contigs, in_contigs_name, in_gff=None, in_gff_name=None):
    current_path = os.getcwd()
    output_path = current_path + '/chimera_break'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    os.chdir('chimera_break')
    if in_gff:
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

    cmd = '{} -ax asm5 -t{} ../../{} {} ' \
          '> pm_against_ref.sam 2> pm_contigs_against_ref.sam.log'.format(m_path, num_threads,
                                                                                        in_reference_file, '../ragoo.fasta')
    if not os.path.isfile('pm_against_ref.sam'):
        run(cmd)

    os.chdir(current_path)


def get_SVs():
    current_path = os.getcwd()
    os.chdir('pm_alignments')
    # Change this when setup.py is ready. Just call script directly
    cmd = 'python3 ~/Projects/RaGOO/sam2delta.py pm_against_ref.sam'
    if not os.path.isfile('pm_against_ref.sam.delta'):
        run(cmd)

    cmd_2 = 'python3 ~/Projects/RaGOO/Assemblytics_uniq_anchor.py --delta pm_against_ref.sam.delta --unique-length 10000 --out assemblytics_out --keep-small-uniques'
    if not os.path.isfile('assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta'):
        run(cmd_2)

    cmd_3 = '~/Projects/RaGOO/Assemblytics_between_alignments.pl assemblytics_out.coords.tab 50 10000 all-chromosomes exclude-longrange bed > assemblytics_out.variants_between_alignments.bed'
    if not os.path.isfile('assemblytics_out.variants_between_alignments.bed'):
        run(cmd_3)

    cmd_4 = 'python3 ~/Projects/RaGOO/Assemblytics_within_alignment.py --delta assemblytics_out.Assemblytics.unique_length_filtered_l10000.delta --min 50 > assemblytics_out.variants_within_alignments.bed'
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

    os.chdir(current_path)


if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description='order and orient contigs according to minimap2 alignments to a reference')
    parser.add_argument("contigs", metavar="<contigs.fasta>", type=str, help="fasta file with contigs to be ordered and oriented")
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file")
    #parser.add_argument("-o", metavar="PATH", type=str, default="ragoo_output", help="output directory name")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-gff", metavar="<annotations.gff>", type=str, default='', help="Annotations for lift over")
    parser.add_argument("-m", metavar="PATH", type=str, default="", help='path to minimap2 executable')
    parser.add_argument("-b", action='store_true', default=False, help="Break chimeric contigs")
    parser.add_argument("-p", metavar="5", type=int, default=5, help="(for chimera breaking) if this percent or more aligns to additional chromosomes, break.")
    parser.add_argument("-l", metavar="10000", type=int, default=10000, help="(for chimera breaking) minimum alignment lengths to consider.")
    parser.add_argument("-r", metavar="100000", type=int, default=100000, help='(for chimera breaking) minimum ranges to consider')
    parser.add_argument("-c", metavar="1000000", type=int, default=1000000, help="When findng intrachromosomal chimeras, minimum gap length with respect to the query.")
    parser.add_argument("-d", metavar="20000000", type=int, default=20000000, help="When findng intrachromosomal chimeras, minimum gap length with respect to the reference.")
    parser.add_argument("-t", metavar="3", type=int, default=3, help="Number of threads when running minimap.")

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

    current_path = os.getcwd()
    output_path = current_path + '/ragoo_output'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir(output_path)

    # Run minimap2
    cmd = '{} -k19 -w19 -t{} ../{} ../{} ' \
          '> contigs_against_ref.paf 2> contigs_against_ref.paf.log'.format(minimap_path, t, reference_file, contigs_file)

    if not os.path.isfile('contigs_against_ref.paf'):
        run(cmd)

    # Read in the minimap2 alignments just generated
    log('-- Reading alignments')
    alns = read_paf_alignments('contigs_against_ref.paf')
    alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file)

    # Process the gff file
    if gff_file:
        log('-- Getting gff features')
        features = defaultdict(list)
        z = GFFReader('../' + gff_file)
        for i in z.parse_gff():
            features[i.seqname].append(i)

    # Break chimeras if desired
    if break_chimeras:
        alns = clean_alignments(alns, l=15000, in_exclude_file=exclude_file, uniq_anchor_filter=True)
        # Process contigs
        log('-- Getting contigs')
        contigs_dict = read_contigs('../' + contigs_file)

        log('-- Finding interchromosomally chimeric contigs')
        all_chimeras = dict()
        for i in alns.keys():
            ref_parts = get_ref_parts(alns[i], min_len, min_break_pct, min_range)
            if len(ref_parts) > 1:
                all_chimeras[i] = ref_parts

        log('-- Finding break points and breaking interchromosomally chimeric contigs')
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
        log('-- Reading interchromosomal chimera broken alignments')
        inter_alns = read_paf_alignments('chimera_break/inter_contigs_against_ref.paf')
        inter_alns = clean_alignments(inter_alns, l=1000, in_exclude_file=exclude_file)

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
        log('-- Reading intrachromosomal chimera broken alignments')
        alns = read_paf_alignments('chimera_break/intra_contigs_against_ref.paf')
        alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file)
        contigs_file = '/ragoo_output/chimera_break/' + out_intra_fasta



    # Assign each contig to a corresponding reference chromosome.
    log('-- Assigning contigs')
    all_unique_contigs = dict()
    for i in alns.keys():
        all_unique_contigs[i] = UniqueContigAlignment(alns[i])

    # Add to this the list of headers that did not make it
    write_contig_clusters(all_unique_contigs)

    # Get longest alignments
    longest_contigs = dict()
    for i in alns.keys():
        longest_contigs[i] = LongestContigAlignment(alns[i])

    log('-- Orienting Contigs')
    # Save the orientations
    final_orientations = dict()
    for i in longest_contigs.keys():
        final_orientations[i] = longest_contigs[i].strand

    log('-- Ordering Contigs')
    order_contigs(all_unique_contigs)

    log('-- Creating Pseudomolecules')
    create_pseudomolecules(contigs_file, all_unique_contigs, 500)

    log('-- Aligning pseudomolecules to reference')
    align_pms(minimap_path, t, reference_file)

    log('-- Getting structural variants')
    get_SVs()

    log('-- goodbye')