import os
import sys
import argparse
import itertools 

parser = argparse.ArgumentParser(description='Process some inputs')

parser.add_argument('-r', '--ref', type=str, help='Path to reference genome (with prefix)')
parser.add_argument('-f', '--fastq_dir', type=str, help='Path to fastq_files')
parser.add_argument('-o', '--out_dir', type=str, help='output directory')
parser.add_argument('-s', '--samples_file', type=str, help='samples.txt file (having sample names and metadata)')
parser.add_argument('-c', '--controls_file', type=str, help='controls.txt file (having control names and metadata)')
parser.add_argument('-v', '--snp_file', type=str, help='snp file required for SNPsplit')
args = parser.parse_args()

ref = args.ref
fastq_dir = args.fastq_dir
out_dir = args.out_dir
samples_file = args.samples_file
controls_file = args.controls_file
snp_file = args.snp_file

sample_names = open(samples_file,'r')
control_names = open(controls_file, 'r')

sample_control = []
controls = []
trimming = []
encodings = []
libraries = []
peak_types = []
effective_genomes = []
allele_specifics = []
alleles_1 = []
alleles_2 = []

total_samples = 0

for line in sample_names:
	total_samples += 1
	fields = line.rstrip().split("\t")
	sample_control.append(fields[0].rstrip())
	controls.append(fields[1].rstrip())
	trimming.append(fields[2].rstrip())
	encodings.append(fields[3].rstrip())
	libraries.append(fields[4].rstrip())
	effective_genomes.append(fields[5].rstrip())
	allele_specifics.append(fields[6].rstrip())
	alleles_1.append(fields[7].rstrip())
	alleles_2.append(fields[8].rstrip())
	peak_types.append(fields[9].rstrip())


for line in control_names:
	fields = line.rstrip().split("\t")
	sample_control.append(fields[0].rstrip())
	trimming.append(fields[1].rstrip())
	encodings.append(fields[2].rstrip())
	libraries.append(fields[3].rstrip())
	effective_genomes.append(fields[4].rstrip())
	allele_specifics.append(fields[5].rstrip())
	alleles_1.append(fields[6].rstrip())
	alleles_2.append(fields[7].rstrip())
	peak_types.append(fields[8].rstrip())


sample_names.close()
control_names.close()

##########################################################################################################################################################################################################
################################################################################## bowtie alignment pbs script creation ##################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a bowtie2-align script\n#PBS -N bowtie2-align\n#PBS -l nodes=1:ppn=20\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o bowtie2-align_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	qual_encoding = encodings[i]
	library = libraries[i]
	trimmed = trimming[i]

	pbs_file_name = out_dir + folder_name + "/bowtie2_align.pbs"

	os.system("mkdir -p " + out_dir + folder_name)

	cmd_1 = "module load bowtie2\nmodule load samtools\n\n"
	cmd_2 = "cd " + out_dir + folder_name + "\n\n"

	if (library == "paired"):
		if (trimmed == "no"):
			cmd_3 = "( time bowtie2 -p 6 --" + qual_encoding + " -x " + ref + " -1 " + fastq_dir + "raw/" + folder_name + "_1.fastq.gz -2 " + fastq_dir + "raw/" + folder_name + "_2.fastq.gz|samtools view -q 20 -@ 6 -O BAM -o " + folder_name + ".bam ) > bowtie2_align_out.log 2>&1\n\n"
		elif (trimmed == "yes"):
			cmd_3 = "( time bowtie2 -p 6 --" + qual_encoding + " -x " + ref + " -1 " + fastq_dir + "trimmed/" + folder_name + "_1_val_1.fq.gz -2 " + fastq_dir + "trimmed/" + folder_name + "_2_val_2.fq.gz|samtools view -q 20 -@ 6 -O BAM -o " + folder_name + ".bam ) > bowtie2_align_out.log 2>&1\n\n"	
	elif (library == "single"):
		if (trimmed == "no"):
			cmd_3 = "( time bowtie2 -p 6 --" + qual_encoding + " -x " + ref + " " + fastq_dir + "raw/" + folder_name + "_trimmed.fq.gz " + "|samtools view -q 20 -@ 6 -O BAM -o " + folder_name + ".bam ) > bowtie2_align_out.log 2>&1\n\n"
		elif (trimmed == "yes"):
			cmd_3 = "( time bowtie2 -p 6 --" + qual_encoding + " -x " + ref + " " + fastq_dir + "trimmed/" + folder_name + "_trimmed.fq.gz " + "|samtools view -q 20 -@ 6 -O BAM -o " + folder_name + ".bam ) > bowtie2_align_out.log 2>&1\n\n"

	cmd_4 = "qsub samtools_markdup.pbs\n"
	pbs_script = open(pbs_file_name,'w+')
	pbs_script.write('%s\n'%(pbs_header + cmd_1 + cmd_2 + cmd_3 + cmd_4))
	pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## samtools markdup pbs script creation ##################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a samtools duplicate marking pipeline script\n#PBS -N mark_dup\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o samtools_markdup_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

cmd_4 = ""

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	allele_specific = allele_specifics[i]
	library = libraries[i]
	allele_1 = alleles_1[i]
	allele_2 = alleles_2[i]
	pbs_file_name = out_dir + folder_name + "/samtools_markdup.pbs"
	fixmate = ""

	cmd_1 = "module load samtools\n\n"
	cmd_2 = "cd " + out_dir + folder_name + "\n\n"
	cmd_7 = "( time samtools view -q 20 -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_view_filter_out.log 2>&1\n\n"

	cmd_3 = "( time samtools sort -@ 13 -n -O BAM -o " + folder_name + "_NS.bam " + folder_name + ".bam ) > samtools_name_sort_out.log 2>&1\n"
	if (library == "paired"):
		fixmate = "_fixmate"
		cmd_4 = "( time samtools fixmate -@ 13 -m -O BAM " + folder_name + "_NS.bam " + folder_name + "_NS" + fixmate + ".bam ) > samtools_fixmate_out.log 2>&1\n"
		cmd_7 = "( time samtools view -q 20 -f 3 -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_view_filter_out.log 2>&1\n\n"
	cmd_5 = "( time samtools sort -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS" + fixmate + ".bam ) > samtools_coord_sort_out.log 2>&1\n"
	cmd_6 = "( time samtools markdup -@ 13 -r -O BAM " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_markdup_out.log 2>&1\n"
	cmd_8 = ("rm " + folder_name + "_NS.bam " + folder_name + "_NS_fixmate.bam " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS_fixmate_CS_dedup.bam\n\n"  )
	cmd_9 = ("samtools index " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam\n")
	cmd_10 = ("samtools idxstats " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam > " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_idxstats.txt\n\n")

	if (allele_specific == "yes"):
		cmd_11 = "qsub SNPsplit.pbs\n"
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4 + cmd_5 + cmd_6 + cmd_7 + cmd_8 + cmd_9 + cmd_10 + cmd_11
	else:
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_5 + cmd_6 + cmd_7 + cmd_8 + cmd_9 + cmd_10
		cmd = cmd + "cd " + out_dir + "\ncd ..\n\n"
		cmd = cmd + "python scripts/bam_normalize.py " + out_dir + " " + folder_name + " " + allele_specific + " " + allele_1 + " " + allele_2 + " " + ref + " " + library + "\n\n"
		cmd = cmd + "cd " + out_dir + folder_name + "\n\n"
		cmd = cmd + "qsub bam_normalize.pbs\n"

	pbs_script = open(pbs_file_name,'w+')
	pbs_script.write('%s\n'%(pbs_header + cmd))
	pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## SNP_split pbs script creation #########################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a SNPsplit script\n#PBS -N SNPsplit\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o SNPsplit_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	library = libraries[i]
	allele_specific = allele_specifics[i]
	allele_1 = alleles_1[i]
	allele_2 = alleles_2[i]

	if (allele_specific == "yes"):
		pbs_file_name = out_dir + folder_name + "/SNPsplit.pbs"
		cmd_1 = "module load SNPsplit\nmodule load samtools\n\n"
		cmd_2 = "cd " + out_dir + folder_name + "\n\n"

		if (library == "paired"):
			cmd_3 = "( time SNPsplit --paired --no_sort --snp_file "  + snp_file + " " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam ) > SNPsplit_out.log 2>&1\n\n"
		elif (library == "single"):
			cmd_3 = "( time SNPsplit --bam --no_sort --snp_file "  + snp_file + " " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam ) > SNPsplit_out.log 2>&1\n\n"

		cmd_4 = "rm " + folder_name + "_NS_fixmate_CS_dedup_filtered.allele_flagged.bam\n\n"
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4
		cmd = cmd + "cd " + out_dir + "\ncd ..\n\n"
		cmd = cmd + "python scripts/bam_normalize.py " + out_dir + " " + folder_name + " " + allele_specific + " " + allele_1 + " " + allele_2 + " " + ref + " " + library + "\n\n"
		cmd = cmd + "cd " + out_dir + folder_name + "\n\n"
		cmd = cmd + "qsub bam_normalize.pbs\n"

		pbs_script = open(pbs_file_name,'w+')
		pbs_script.write('%s\n'%(pbs_header + cmd))
		pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## macs2 callpeak pbs script creation ####################################################################################
##########################################################################################################################################################################################################

# pbs_header = "# This is a macs2_callpeak script\n#PBS -N macs2_callpeak\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=48:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o macs2_callpeak.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

# for i in range(total_samples):

# 	folder_name = sample_control[i]
# 	control_name = out_dir + controls[i] + "/" + controls[i]
# 	peak_type = peak_types[i]
# 	control = controls[i]
# 	effective_genome = effective_genomes[i]
# 	allele_specific = allele_specifics[i]
# 	allele_1 = alleles_1[i]
# 	allele_2 = alleles_2[i]
# 	library =libraries[i]

# 	pbs_file_name = out_dir + folder_name + "/macs2_callpeak.pbs"

# 	file_format = ""

# 	if (library == "paired"):
# 		file_format = file_format + " -f BAMPE "

# 	cmd_1 = "module load macs2\n\n"
# 	cmd_2 = "cd " + out_dir + folder_name + "\n\n"
	

# 	if (allele_specific == "yes"):
# 		if (control != "NA"):
# 			cmd_3 = "mkdir -p " + allele_1 + "\n"
# 			cmd_4 = "mkdir -p " + allele_2 + "\n"
# 			cmd_5 = "mkdir -p comp\n\n"
# 			cmd_6 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.genome1.bam -c " + control_name + "_NS_fixmate_CS_dedup_filtered.genome1.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n " + allele_1 + "_" + folder_name + " --outdir " + allele_1 + " )  > macs2_allele_1_out.log 2>&1\n"
# 			cmd_7 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.genome2.bam -c " + control_name + "_NS_fixmate_CS_dedup_filtered.genome2.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n " + allele_2 + "_" + folder_name + " --outdir " + allele_2 + " )  > macs2_allele_2_out.log 2>&1\n"
# 			cmd_8 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam -c " + control_name + "_NS_fixmate_CS_dedup_filtered.bam" + file_format +  "--" + peak_type + " --SPMR -B -g " + effective_genome + " -n comp_" + folder_name + " --outdir comp )  > macs2_comp_out.log 2>&1\n"
# 			cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4 + cmd_5 + cmd_6 + cmd_7 + cmd_8 + "\n\n"

# 		elif (control == "NA"):
# 			cmd_3 = "mkdir " + allele_1 + "\n"
# 			cmd_4 = "mkdir " + allele_2 + "\n"
# 			cmd_5 = "mkdir comp\n"  
# 			cmd_6 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.genome1.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n " + allele_1 + "_" + folder_name + " --outdir " + allele_1 + " ) > macs2_allele_1_out.log 2>&1\n"
# 			cmd_7 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.genome2.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n " + allele_2 + "_" + folder_name + " --outdir " + allele_2 + " ) > macs2_allele_2_out.log 2>&1\n"
# 			cmd_8 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n comp_" + folder_name + " --outdir comp )  > macs2_comp_out.log 2>&1\n"
# 			cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4 + cmd_5 + cmd_6 + cmd_7 + cmd_8 


# 	elif (allele_specific == "no"):
# 		if (control != "NA"):
# 			cmd_3 = "mkdir comp\n"
# 			cmd_4 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam -c " + control_name + "_NS_fixmate_CS_dedup_filtered.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n comp_" + folder_name + " --outdir comp )  > macs2_comp_out.log 2>&1\n"
# 			cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4
		
# 		elif (control == "NA"):
# 			cmd_3 = "mkdir comp\n"
# 			cmd_4 = "( time macs2 callpeak -t " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam" + file_format +  "--" + peak_type + " -B -g " + effective_genome + " -n comp_" + folder_name + " --outdir comp )  > macs2_comp_out.log 2>&1\n"
# 			cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4

# 	cmd = cmd + "cd " + out_dir + "\ncd ..\n\n"
# 	cmd = cmd + "python scripts/bedgraph_normalize.py " + out_dir + " " + folder_name + " " + allele_specific + " " + allele_1 + " " + allele_2 + " " + ref + "\n\n"
# 	cmd = cmd + "cd " + out_dir + folder_name + "\n\n"
# 	cmd = cmd + "qsub bedgraph_normalize.pbs\n"

# 	pbs_script = open(pbs_file_name,'w+')
# 	pbs_script.write('%s\n'%(pbs_header + cmd))
# 	pbs_script.close()
