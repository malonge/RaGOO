from utilities.PAFReader import PAFReader
from utilities.ContigAlignment import ContigAlignment
from utilities.ContigAlignment import UniqueContigAlignment
from utilities.ContigAlignment import LongestContigAlignment
from utilities.utilities import run, log
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


def clean_alignments(in_alns, l=10000, in_exclude_file=''):
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
        if len(in_alns[header].ref_headers) == 0:
            empty_headers.append(header)

    for header in empty_headers:
        in_alns.pop(header)
    return in_alns


def read_paf_alignments():
    # Read in PAF alignments
    # Initialize a dictionary where key is contig header, and value is ContigAlignment.
    alns = dict()
    x = PAFReader('contigs_against_ref.paf')
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


if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description='order and orient contigs according to minimap2 alignments to a reference')
    parser.add_argument("contigs", metavar="<contigs.fasta>", type=str, help="fasta file with contigs to be ordered and oriented")
    parser.add_argument("reference", metavar="<reference.fasta>", type=str, help="reference fasta file")
    #parser.add_argument("-o", metavar="PATH", type=str, default="ragoo_output", help="output directory name")
    parser.add_argument("-e", metavar="<exclude.txt>", type=str, default="", help="single column text file of reference headers to ignore")
    parser.add_argument("-m", metavar="PATH", type=str, default="", help='path to minimap2 executable')

    args = parser.parse_args()
    args = parser.parse_args()
    contigs_file = args.contigs
    reference_file = args.reference
    #output_path = args.o
    exclude_file = args.e
    minimap_path = args.m

    current_path = os.getcwd()
    output_path = current_path + '/ragoo_output'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir(output_path)

    # Run minimap2
    cmd = '/usr/local/src/minimap2/minimap2 -k19 -w19 -t3 ../{} ../{} ' \
          '> contigs_against_ref.paf 2> contigs_against_ref.paf.log'.format(reference_file, contigs_file)

    if not os.path.isfile('contigs_against_ref.paf'):
        run(cmd)

    log('-- Reading alignments')
    alns = read_paf_alignments()
    alns = clean_alignments(alns, l=13000, in_exclude_file=exclude_file)

    # Assign each contig to a corresponding reference chromosome.
    log('-- Assigning contigs')
    all_unique_contigs = dict()
    for i in alns.keys():
        all_unique_contigs[i] = UniqueContigAlignment(alns[i])

    # Add to this the list of headers that did not make it
    write_contig_clusters(all_unique_contigs)

    log('-- Getting longest alignments')
    longest_contigs = dict()
    for i in alns.keys():
        longest_contigs[i] = LongestContigAlignment(alns[i])

    log('-- Orienting Contigs')
    # Save the orientations
    final_orientations = dict()
    for i in longest_contigs.keys():
        final_orientations[i] = longest_contigs[i].strand

    log('-- Ordering Contigs')
    current_path = os.getcwd()
    output_path = current_path + '/orderings'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    all_chroms = set([all_unique_contigs[i].ref_chrom for i in all_unique_contigs.keys()])
    for this_chrom in all_chroms:

        # Intialize the list of start and end positions w.r.t the query
        ref_pos = []

        contigs_file = 'groupings/' + this_chrom + '_contigs.txt'
        contigs_list = get_contigs(contigs_file)

        # Add edges for each of n choose 2 pairs of alignments
        for i in range(len(contigs_list)):
            ref_pos.append((longest_contigs[contigs_list[i]].ref_start, longest_contigs[contigs_list[i]].ref_end, i))

        final_order = [contigs_list[i[2]] for i in sorted(ref_pos)]

        with open('orderings/' + this_chrom + '_orderings.txt', 'w') as out_file:
            for i in final_order:
                out_file.write(i + '\t' + final_orientations[i] + '\n')

    log('goodbye')