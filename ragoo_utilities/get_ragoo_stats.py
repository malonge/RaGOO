import argparse

parser = argparse.ArgumentParser(description='Summary stats about contig scaffolding with RaGOO.')
parser.add_argument("index", metavar="<contigs.fasta.fai>", type=str, help="Samtools fasta index file for input contigs. If chimera breaking mode was used, this must be"
                                                                           "the index file of the chimera broken contigs, which can be found in ragoo_output/chimera_break."
                                                                           "The correct file to use is the file with the .intra.chimera.broken.fa suffix.")
parser.add_argument("groupings", metavar="<groupings.fofn>", type=str, help="file of file names for all *_groupings.txt produced by RaGOO. Single column with full path to each grouping file.")

args = parser.parse_args()
contigs_index = args.index
grouping_fofn = args.groupings

remaining_ctg = []
all_ctg_len = dict()
with open(contigs_index) as f:
    for line in f:
        L1 = line.split('\t')
        all_ctg_len[L1[0]] = int(L1[1])
        remaining_ctg.append(L1[0])

grouping_files = []
with open(grouping_fofn) as f:
    for line in f:
        grouping_files.append(line.rstrip())

num_ctg_localized = 0
num_bp_localized = 0

for group_file in grouping_files:
    with open(group_file) as f:
        for line in f:
            L1 = line.split('\t')
            header = L1[0].rstrip()
            num_ctg_localized += 1
            num_bp_localized += all_ctg_len[header]
            assert header in remaining_ctg
            remaining_ctg.pop(remaining_ctg.index(header))

num_ctg_unlocalized = 0
num_bp_unlocalized = 0
for ctg in remaining_ctg:
    num_ctg_unlocalized += 1
    num_bp_unlocalized += all_ctg_len[ctg]

print('%r contigs were localized by RaGOO' %(num_ctg_localized))
print('%r bp were localized by RaGOO' %(num_bp_localized))
print('%r contigs were unlocalized by RaGOO' %(num_ctg_unlocalized))
print('%r bp were unlocalized by RaGOO' %(num_bp_unlocalized))

print('%r %% of contigs were localized by RaGOO' %((num_ctg_localized/(num_ctg_localized + num_ctg_unlocalized))*100))
print('%r %% of bp were localized by RaGOO' %((num_bp_localized/(num_bp_localized + num_bp_unlocalized))*100))




