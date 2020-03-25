import os
import sys

pbs_header = "# This is a bedgraph normalization script\n#PBS -N normalize_bedgraphsdg\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o normal.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

out_dir = sys.argv[1]
folder_name = sys.argv[2]
allele_specific = sys.argv[3]
allele_1 = sys.argv[4]
allele_2 = sys.argv[5]
ref = sys.argv[6]
library = sys.argv[7]

pbs_file_name = out_dir + folder_name + "/bam_normalize.pbs"

fixmate = ""

if (library == "paired"):
	fixmate = "_fixmate"

idxstats = open(out_dir + folder_name + "/" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_idxstats.txt", 'r')

total_chrom = 0
total_mapped_reads = 0 

cmd = "cd " + out_dir + folder_name + "\n\n"

cmd = cmd + "module load bedtools\n\n"

for line in idxstats:
	total_chrom += 1
	if (total_chrom <= 22):
		fields = line.rstrip().split('\t')
		reads = fields[2].rstrip()
		total_mapped_reads = total_mapped_reads + int(reads)

if (allele_specific == "yes"):
	if (library == "paired"):
		cmd = cmd + "bedtools genomecov -pc -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bam -g " + ref + ".fa.fai > " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bdg\n"
		cmd = cmd + "bedtools genomecov -pc -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bam -g " + ref + ".fa.fai > " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bdg\n"		
	elif (library == "single"):
		cmd = cmd + "bedtools genomecov -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" +  fixmate + "_CS_dedup_filtered.genome1.bam -g " + ref + ".fa.fai > " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bdg\n"
		cmd = cmd + "bedtools genomecov -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" +  fixmate + "_CS_dedup_filtered.genome2.bam -g " + ref + ".fa.fai > " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bdg\n"
	
	cmd = cmd + "sort -k1,1 -k2,2n -k3,3n " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bdg > " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1_sorted.bdg\n"
	cmd = cmd + "rm " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bdg\n"
	cmd = cmd + "bedGraphToBigWig " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1_sorted.bdg " + ref + ".fa.fai " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1_sorted.bw\n\n"
	
	cmd = cmd + "sort -k1,1 -k2,2n -k3,3n " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bdg > " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2_sorted.bdg\n"
	cmd = cmd + "rm " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bdg\n"
	cmd = cmd + "bedGraphToBigWig " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2_sorted.bdg " + ref + ".fa.fai " + allele_2 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2_sorted.bw\n\n"

if (library == "paired"):
	cmd = cmd + "bedtools genomecov -pc -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam -g " + ref + ".fa.fai > comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bdg\n\n"
else:
	cmd = cmd + "bedtools genomecov -scale " + str((1000000 * 2)/total_mapped_reads) +  " -bga -ibam " +  folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam -g " + ref + ".fa.fai > comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bdg\n\n"

cmd = cmd + "sort -k1,1 -k2,2n -k3,3n comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bdg > comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_sorted.bdg\n"
cmd = cmd + "rm comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bdg\n"
cmd = cmd + "bedGraphToBigWig comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_sorted.bdg " + ref + ".fa.fai comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_sorted.bw\n\n"

idxstats.close()
pbs_script = open(pbs_file_name,'w+')
pbs_script.write('%s\n'%(pbs_header + cmd))
pbs_script.close()
