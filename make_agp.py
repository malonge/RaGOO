#!/usr/bin/env python
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce an AGP v2.0 file describing the ordering and orienting performed by RaGOO.')
    parser.add_argument("orderings", metavar="<orderings.fofn>", type=str, help="file of orderings files ($ls ragoo_output/orderings/* > orderings.fofn)")
    parser.add_argument("fai", metavar="<contigs.fasta.fai>", type=str, help="index file for contigs ($samtools faidx contigs.fasta)")
    parser.add_argument("gap_len", metavar="gap_size", type=int, help="Gap size used for pseudomolecule padding.")

    # Get the command line arguments
    args = parser.parse_args()
    orderings_fofn = args.orderings
    fai_file = args.fai
    gap_len = args.gap_len

    # Save the contig orderings
    orderings = dict()
    with open(orderings_fofn, 'r') as f1:
        for line1 in f1:
            chrom = line1[line1.rfind("/")+1:line1.rfind("_orderings.txt")] + "_RaGOO"
            orderings[chrom] = []
            with open(line1.rstrip(), 'r') as f2:
                for line2 in f2:
                    L1 = line2.rstrip().split('\t')
                    # Append (ctg_header, orientation)
                    orderings[chrom].append((L1[0], L1[1]))

    # Get contig lengths
    ctg_lens = dict()
    with open(fai_file, 'r') as f:
        for line in f:
            L1 = line.split('\t')
            ctg_lens[L1[0]] = int(L1[1])

    # Write the final AGP file
    current_pos = 1
    idx = 1
    chrom_idx = 1
    out_buff = []

    sys.stdout.write("## AGP-version 2.0\n")
    sys.stdout.write("## AGP constructed by RaGOO\n")
    for chrom in orderings:
        for ctg, strand in orderings[chrom]:
            start = current_pos
            end = current_pos + ctg_lens[ctg] - 1
            line_buff = list()

            # Write the line for the contig
            line_buff.append(chrom)
            line_buff.append(str(start))
            line_buff.append(str(end))
            line_buff.append(str(idx))
            line_buff.append("W")
            line_buff.append(ctg)
            line_buff.append("1")
            line_buff.append(str(end - start + 1))
            line_buff.append(strand)
            out_buff.append("\t".join(line_buff))

            # Write the line for the gap
            idx += 1
            current_pos += ctg_lens[ctg]
            start = current_pos
            end = current_pos + gap_len - 1
            line_buff = list()

            line_buff.append(chrom)
            line_buff.append(str(start))
            line_buff.append(str(end))
            line_buff.append(str(idx))
            line_buff.append("N")
            line_buff.append(str(gap_len))
            line_buff.append("scaffold")
            line_buff.append("yes")
            line_buff.append("align_genus")
            out_buff.append("\t".join(line_buff))

            idx += 1
            current_pos += gap_len

        # Pop the last gap since we don't end with padding
        out_buff.pop()
        idx =1
        chrom_idx += 1
        current_pos = 1

    sys.stdout.write("\n".join(out_buff) + "\n")