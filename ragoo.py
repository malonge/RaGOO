from utilities.PAFReader import PAFReader
from utilities.SeqReader import SeqReader
from utilities.ContigAlignment import ContigAlignment
from utilities.ContigAlignment import UniqueContigAlignment
from utilities.ContigAlignment import LongestContigAlignment
from utilities.utilities import run, log, reverse_complement
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

    # Read in the minimap2 alignments just generated
    log('-- Reading alignments')
    alns = read_paf_alignments()
    alns = clean_alignments(alns, l=1000, in_exclude_file=exclude_file)

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

    log('goodbye')