
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='given and orderings file and a contigs fasta index, print a bed file of contig placements in the pseudomolecules.')
    parser.add_argument("orderings", metavar="<orderings.txt>", type=str, help="orderings file from RaGOO")
    parser.add_argument("fai", metavar="<contigs.fasta.fai>", type=str, help="index file for contigs (samtools faidx contigs.fasta)")
    parser.add_argument("gap_len", metavar="100", type=int, help="Gap size used for pseudomolecule padding.")

    # Get the command line arguments
    args = parser.parse_args()
    orderings_file = args.orderings
    fai_file = args.fai
    gap_len = args.gap_len

    # Save the contig orderings
    ctgs = []
    with open(orderings_file, 'r') as f:
        for line in f:
            ctgs.append(line.rstrip().split('\t')[0])

    # Get contig lengths
    ctg_lens = dict()
    with open(fai_file, 'r') as f:
        for line in f:
            L1 = line.split('\t')
            ctg_lens[L1[0]] = int(L1[1])

    # Get contig borders
    final_bed = []
    current_pos = 0

    for ctg in ctgs:
        start = current_pos
        end = current_pos + ctg_lens[ctg]
        current_pos += ctg_lens[ctg]
        current_pos += gap_len
        pm_header = orderings_file[orderings_file.rfind('/')+1:orderings_file.rfind('_')] + '_RaGOO'
        final_bed.append('%s\t%r\t%r' % (pm_header, start, end))

    print('\n'.join(final_bed))