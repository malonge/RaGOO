#!/usr/bin/env python
def get_contig_lengths(in_fai):
    """ Get contig lengths from a fasta index. """
    lens = dict()
    with open(in_fai, 'r') as f:
        for line in f:
            L1 = line.rstrip().split('\t')
            lens[L1[0]] = int(L1[1])
    return lens


def get_contig_orderings(in_groupings):
    """
    From the orderings files, return the following:

    1. Dictionary associating a reference header with a list of ordered contig headers
    2. Dictionary associating a contig with an orientation
    3. Dicitonary associating a contig with its assigned reference header
    """
    orderings = dict()
    orientations = dict()
    ref_chr = dict()

    # Iterate through the orderings files
    with open(in_groupings, 'r') as f1:
        for line1 in f1:
            header = line1.rstrip().replace('_orderings.txt', '_RaGOO')
            header = header[header.rfind('/') + 1:]

            # Initialize the list for the orderings
            orderings[header] = []
            with open(line1.rstrip(), 'r') as f2:
                for line2 in f2:
                    L1 = line2.rstrip().split('\t')
                    orderings[header].append(L1[0])
                    orientations[L1[0]] = L1[1]
                    ref_chr[L1[0]] = header
    return orderings, orientations, ref_chr


def get_reverse_coords(start_coord, end_coord, seq_length):
    """
    Returns new genomic coordinates of a region that has undergone reverse complementation.
    new start coordinate = seqLength - endCoord
    new end coordinate = seqLength - startCoord
    """
    return seq_length - (end_coord - 1), seq_length - (start_coord - 1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Lift-over gff coordinates to from contigs to RaGOO pseudomolecules.'
                                                 'Please make sure that you are using the updated gff file and set of contigs if chimera correction was used.'
                                                 'Also make sure that the gap size (-g) matches that which was used during ordering and orienting.')
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="Gff file to be lifted-over")
    parser.add_argument("orderings", metavar="<orderings.fofn>", type=str, help="List of RaGOO 'orderings' files. 'ls ragoo_output/orderings/* > orderings.fofn'")
    parser.add_argument("fai", metavar="<contigs.fasta.fai>", type=str, help="Fasta index for contigs.'")
    parser.add_argument("-g", metavar="100", type=int, default=100, help="Gap size for padding in pseudomolecules (must match what was used for 'ragoo.py'.")


    # Get the command line arguments
    args = parser.parse_args()
    gff_file = args.gff
    orderings_file = args.orderings
    fai_file = args.fai
    gap_size = args.g

    # Get the contig lengths
    ctg_lens = get_contig_lengths(fai_file)

    # Get the contig orderings and orientations
    ctg_orderings, ctg_orientations, ctg_chr = get_contig_orderings(orderings_file)

    #Iterate through the GFF features and update
    offsets = dict()
    with open(gff_file) as f:
        for gff_line in f:
            if gff_line.startswith("#"):
                print(gff_line.rstrip())
            else:
                gff_fields = gff_line.rstrip().split('\t')
                gff_header, start, end, strand = gff_fields[0], int(gff_fields[3]), int(gff_fields[4]), gff_fields[6]
                new_header = ctg_chr[gff_header]

                # Check that this contig header is in our list of headers
                if gff_header not in ctg_lens.keys():
                    err_msg = """ %s was not found in the list of orderings files provided.
                    Please check that you are using the correct files. If chimeric contig correction was
                    used, please use the corrected gff and fai file. 
                    """
                    raise ValueError(err_msg)

                # Check if the contig has been reverse complemented. Update accordingly
                if ctg_orientations[gff_header] == '-':
                    start, end = get_reverse_coords(start, end, ctg_lens[gff_header])

                    # Set the opposite strand
                    if strand == '+':
                        strand = '-'
                    else:
                        strand = '+'

                # Check if the offset for this contig has already been calculated
                if gff_header in offsets:
                    offset = offsets[gff_header]
                else:
                    ctg_idx = ctg_orderings[new_header].index(gff_header)
                    offset = 0

                    for ctg in ctg_orderings[new_header][:ctg_idx]:
                        offset += ctg_lens[ctg]
                        offset += gap_size

                    # memoize the offset
                    offsets[gff_header] = offset

                new_start = start + offset
                new_end = end + offset

                gff_fields[0] = new_header
                gff_fields[3] = str(new_start)
                gff_fields[4] = str(new_end)
                gff_fields[6] = strand
                print('\t'.join(gff_fields))
