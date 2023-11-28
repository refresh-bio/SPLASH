# This script generates 3 files for annotated exon boundaries, splice sites, and gene coordinates from a gtf file 

library(data.table)

#sbatch -p quake,normal,horence --time=24:00:00 --mem=60000 --wrap="/oak/stanford/groups/horence/Roozbeh/software/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_files --genomeFastaFiles GCF_021976095.1_XeniaSp_v1_genomic.fna --sjdbGTFfile GCF_021976095.1_XeniaSp_v1_genomic.gtf"
#sbatch -p horence,owners,quake --time=24:00:00 --mem=60000 --wrap="/home/groups/horence/applications/bowtie2-2.2.1/bowtie2-build GCF_021976095.1_XeniaSp_v1_genomic.fna Bowtie_index_files/GCF_021976095.1_XeniaSp_v1_bt2"



args <- commandArgs(trailingOnly = TRUE)
gtf_file = args[1]
hisat2_directory = args[2]
outfile_name = args[3]


##############################################################################
############# making the file for annotated splice sites file ################
##############################################################################
system(paste(hisat2_directory,"/extract_splice_sites.py ", gtf_file, " > ", outfile_name, "_known_splice_sites.txt", sep=""))

##############################################################################
##############################################################################
##############################################################################


###############################################################################
###############  making the file for annotated exon coordinates ###############
##############################################################################
system(paste(hisat2_directory,"/hisat2_extract_exons.py ", gtf_file, " > ", outfile_name, "_exon_coordinates.bed", sep=""))
f = fread(paste(outfile_name, "_exon_coordinates.bed", sep = ""))
f[,V2_1:=V2+1,by=V2]
f[,V2_2:=V2-1,by=V2]
f[,V3_2:=V3-1,by=V3]
f[,V3_1:=V3+1,by=V3]
f[,chr_V2:=paste(V1,V2,sep=""),by=1:nrow(f)]
f[,chr_V2_1:=paste(V1,V2_1,sep=""),by=1:nrow(f)]
f[,chr_V2_2:=paste(V1,V2_2,sep=""),by=1:nrow(f)]
f[,chr_V3:=paste(V1,V3,sep=""),by=1:nrow(f)]
f[,chr_V3_1:=paste(V1,V3_1,sep=""),by=1:nrow(f)]
f[,chr_V3_2:=paste(V1,V3_2,sep=""),by=1:nrow(f)]
write.table(f, paste(outfile_name, "_exon_coordinates.bed", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
################################################################################
################################################################################
##############################################################################

###############################################################################
################ commands for making the gene coordinates files ###############
##############################################################################
genes = fread(gtf_file, skip = ";")
genes = genes[V3 == "gene"]
genes[, gene_name:=strsplit(V9,split=";")[[1]][1], by = V9]
genes[, gene_name:=gsub("gene_id ","",gene_name), by = gene_name]
genes[, gene_name:=gsub("\"","",gene_name), by = gene_name]
genes[, V9:=NULL]
genes = genes[,list(V1,V4,V5,gene_name,V6,V7)]
write.table(genes, paste(outfile_name, "_genes.bed", sep = ""), row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
################################################################################
################################################################################
##############################################################################
