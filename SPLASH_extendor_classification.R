if (!require("stringdist")) {
  install.packages("stringdist", dependencies = TRUE)
  library(stringdist)
}

if (!require("Biostrings")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
  library(Biostrings)
}

if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

if (!require("GenomicAlignments")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GenomicAlignments")
  library(GenomicAlignments)
}

if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if (!require("tictoc")) {
  install.packages("tictoc", dependencies = TRUE)
  library(tictoc)
}


riffle <- function(a, b) {   # this function interleaves the elements of two vectors into a vector
  seqmlab <- seq(length=length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}

##################################################
############ input arguments #####################
args <- commandArgs(trailingOnly = TRUE)
SPLASH_directory = args[1]                 # the output directory used for the SPLASH run
which_anchors_file = args[2]               # specify which anchor file (after correction or all anchors) to use, could be "after_correction" or "all" 
effect_size_cutoff = args[3]               # the effect size cutoff for significant anchors (suggested 0.2) 
num_samples_cutoff = args[4]               # the minimum number of samples for an anchor to be called (suggested 10)
STAR_executable = args[5]                  # path to STAR executable file
Samtools_executable = args[6]              # path to Samtools executable file
bedtools_executable = args[7]              # path to Bedtools executable file
Bowtie2_executable = args[8]               # path to Bowtie2 executable file
STAR_reference = args[9]                   # path to STAR index files for the reference genome
Bowtie2_reference = args[10]               # path to the Bowtie2 index for the reference genome
Bowtie2_univec_index = args[11]            # path to the Bowtie2 index for univec
annotated_splice_juncs_file = args[12]     # path to the file containing annotated splice junctions
annotated_exon_boundaries_file = args[13]  # path to the file containing annotated exon boundaries
gene_coords_file = args[14]                # path to the file containing gene coordinates
repeats_annotation_file = args[15]         # path to the file containing repeats coordinates
##################################################
##################################################

setwd(SPLASH_directory)
print(SPLASH_directory)

tic("reading anchors file")
if (which_anchors_file == "all"){
  anchors_file = paste(SPLASH_directory,"result.after_correction.all_anchors.tsv",sep="")
} else if (which_anchors_file=="after_correction"){
  anchors_file = paste(SPLASH_directory,"result.after_correction.scores.tsv",sep="")
}

anchors = fread(anchors_file, sep = "\t")
anchors = anchors[effect_size_bin > effect_size_cutoff & number_nonzero_samples > num_samples_cutoff ]

anchors = anchors[!(most_freq_target_1=="-" | most_freq_target_2=="-")]
anchors_target_2 = anchors
anchors_target_2$help=1
anchors_target_2[,help:=NULL]
anchors_target_2[,c("most_freq_target_1","cnt_most_freq_target_1"):=NULL]
anchors[,c("most_freq_target_2","cnt_most_freq_target_2"):=NULL]
setnames(anchors_target_2,c("most_freq_target_2","cnt_most_freq_target_2"),c("target","target_count"))
setnames(anchors,c("most_freq_target_1","cnt_most_freq_target_1"),c("target","target_count"))
anchors = rbind(anchors,anchors_target_2)

setnames(anchors,c("M","anch_uniqTargs"),c("anchor_count","num_extendor_per_anchor"))
anchors = anchors[num_extendor_per_anchor > 1]  
anchors$extendor=paste(anchors$anchor,anchors$target,sep="")
toc()

anchors = setorder(anchors, -number_nonzero_samples, anchor,-target_count)
anchors[,extendor_order:=1:2,by=anchor]
anchor_rank_dt = data.table(unique(anchors$anchor),1:length(unique(anchors$anchor)))
names(anchor_rank_dt) = c("anchor","anchor_index")
anchors = merge(anchors,anchor_rank_dt,all.x=TRUE,all.y=FALSE,by.x="anchor",by.y="anchor")

##### removing those anchors with polyAs or those with one of the top two targets with polyAs ##################
anchors = anchors[! ((extendor %like% "AAAAAAAAA") | (extendor %like% "TTTTTTTTT") | (extendor %like% "GGGGGGGGG") | (extendor %like% "CCCCCCCCC"))]

anchors[,num_target_per_anchor:=length(unique(target)),by=anchor]
anchors = anchors[num_target_per_anchor > 1]  
anchors[,num_target_per_anchor:=NULL]
################################################################################################################

anchors_high_rank = anchors[extendor_order<3]
anchors_high_rank_shifted = anchors_high_rank[, data.table::shift(.SD, 1, NA, "lead", TRUE), .SDcols=1:ncol(anchors_high_rank)]
anchors_high_rank = cbind(anchors_high_rank,anchors_high_rank_shifted[,list(anchor_lead_1,target_lead_1)])
anchors_high_rank = anchors_high_rank[anchor_lead_1==anchor]
tic("compute lev and hamming distance")
anchors_high_rank[,lev_dist:=stringdist(target,target_lead_1, method = "lv"),by=1:nrow(anchors_high_rank)]
anchors_high_rank[,ham_dist:=stringdist(target,target_lead_1, method = "hamming"),by=1:nrow(anchors_high_rank)]
toc()
tic("compute run_D and run_I")
anchors_high_rank[,lev_operations:=as.character(attributes(adist(target,target_lead_1, count=T))$trafos),by=1:nrow(anchors_high_rank)]
anchors_high_rank[,lcs_operations:=as.character(attributes(adist(target,target_lead_1, count=T, costs = list(ins=1, del=1, sub=100)))$trafos),by=1:nrow(anchors_high_rank)]
anchors_high_rank[lev_operations%like%"D",run_length_D:=max(rle( strsplit(lev_operations,"")[[1]] == 'D')$lengths[which( rle( strsplit(lev_operations,"")[[1]] == 'D')$values=="TRUE")]),by=lev_operations]
anchors_high_rank[lev_operations%like%"I",run_length_I:=max(rle( strsplit(lev_operations,"")[[1]] == 'I')$lengths[which( rle( strsplit(lev_operations,"")[[1]] == 'I')$values=="TRUE")]),by=lev_operations]
anchors_high_rank[is.na(run_length_I),run_length_I:=0]
anchors_high_rank[is.na(run_length_D),run_length_D:=0]
anchors = merge(anchors,unique(anchors_high_rank[,list(anchor,ham_dist,lev_dist,lev_operations,lcs_operations,run_length_D,run_length_I)]),all.x=TRUE,all.y=FALSE,by.x="anchor",by.y="anchor")
toc()

#################################################################
##### Finding anchors whose target diversity can be explained by a unique base pair change mapping ##########################
#################################################################
tic("RNA editing")
RNA_editing = anchors[,list(anchor,target,anchor_count,target_count,most_freq_target_3,most_freq_target_4,cnt_most_freq_target_3,cnt_most_freq_target_4)]
RNA_editing[,frac_most_freq_target_3:=cnt_most_freq_target_3/anchor_count]
RNA_editing[,frac_most_freq_target_4:=cnt_most_freq_target_4/anchor_count]
RNA_editing_3 = RNA_editing[frac_most_freq_target_3>0.05 & !duplicated(anchor),list(anchor,anchor_count,most_freq_target_3,cnt_most_freq_target_3)]
setnames(RNA_editing_3,c("most_freq_target_3","cnt_most_freq_target_3"),c("target","target_count"))
RNA_editing_4 = RNA_editing[frac_most_freq_target_4>0.05 & !duplicated(anchor),list(anchor,anchor_count,most_freq_target_4,cnt_most_freq_target_4)]
setnames(RNA_editing_4,c("most_freq_target_4","cnt_most_freq_target_4"),c("target","target_count"))
RNA_editing[,c("most_freq_target_3","most_freq_target_4","cnt_most_freq_target_3","cnt_most_freq_target_4","frac_most_freq_target_3","frac_most_freq_target_4"):=NULL]
RNA_editing = rbind(RNA_editing,RNA_editing_3,RNA_editing_4)


for (nt1 in c("A","C","G","T")){
  for (nt2 in c("A","C","G","T")){
    if(nt1!=nt2){
      RNA_editing[,edited_target:=gsub(nt1,nt2,paste(target))] # changing all nt1 to nt2 in each target 
      print(paste("subs", nt1, nt2))
      
      RNA_editing[,num_unique_edit:=length(unique(edited_target)),by="anchor"]
      RNA_editing[num_unique_edit==1,RNA_editing:=paste(nt1,"_to_",nt2,sep="")]
    }
  }
}

anchors = merge(anchors,RNA_editing[!duplicated(anchor),list(anchor,RNA_editing)],all.x=TRUE,all.y=FALSE,by.x="anchor",by.y="anchor")
remove(RNA_editing)
remove(RNA_editing_3)
remove(RNA_editing_4)
toc()
#################################################################
#################################################################
#################################################################



###########################################################################
############# STAR alignment of the extendors to the genome ##############
###########################################################################

### I first add ranks for the anchors
tic("make extendor index")
anchors$extendor_index = paste(">",anchors$anchor_index,"_",anchors$extendor_order,sep="")
toc()

tic("riffle to make fasta")
extendors_fasta = riffle(anchors$extendor_index,anchors$extendor)
extendors_fasta  = data.table(extendors_fasta) # the fasta file containing all extendor sequences
toc()
tic("STAR alignment")
write.table(extendors_fasta,paste(SPLASH_directory,"extendors_fasta.fa",sep = ""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep = "\t")
system(paste(STAR_executable," --runThreadN 4 --genomeDir  ", STAR_reference," --readFilesIn ", SPLASH_directory,"extendors_fasta.fa", " --outFileNamePrefix ", SPLASH_directory,"STAR_alignment/T2T_extendors"," --twopassMode Basic --alignIntronMax 1000000 --limitOutSJcollapsed 3000000 --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 0 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10 --outSAMtype SAM --chimOutType SeparateSAMold --outSAMunmapped None --clip3pAdapterSeq AAAAAAAAA --outSAMattributes NH HI AS nM NM ",sep = ""))

alignment_info_extendors = fread(paste(SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.sam",sep=""),header=FALSE,skip="NH:",select = c("V1","V2","V3","V4","V6","V16")) ## now grabbing alignment information for extendor sequences after running STAR
alignment_info_extendors$V1=paste(">",alignment_info_extendors$V1, sep = "")
alignment_info_extendors[,num_alignments:=.N,by=V1]
alignment_info_extendors[,STAR_T2T_num_mismatches:=as.numeric(strsplit(V16, split = ":")[[1]][3]), by = V16]

anchors[,c("STAR_T2T_flag","STAR_T2T_chr","STAR_T2T_coord","STAR_T2T_CIGAR","STAR_T2T_num_alignments","STAR_T2T_num_mismatches"):=NULL]
anchors = merge(anchors,alignment_info_extendors[!duplicated(V1),list(V1,V2,V3,V4,V6,num_alignments,STAR_T2T_num_mismatches)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="V1")
setnames(anchors,c("V2","V3","V4","V6","num_alignments"),c("STAR_T2T_flag","STAR_T2T_chr","STAR_T2T_coord","STAR_T2T_CIGAR","STAR_T2T_num_alignments"))

### below I add the flags for the STAR alignment status of each extendor
anchors[,extendor_index:=gsub(">","",extendor_index), by = extendor_index]
anchors[, is.aligned_STAR:=0]   # whether extendor has been mapped by STAR
anchors[, is.STAR_chimeric:=0]  # whether STAR reports a chimeric alignment for the extendor
anchors[, is.STAR_SJ:=0]        # whether STAR reports a gapped alignment (splice junction for the extendor)

num_chimeric_alignments = data.table()
num_chimeric_alignments = fread(paste(SPLASH_directory,"STAR_alignment/T2T_extendorsLog.final.out",sep=""),sep="|",skip=35)
num_chimeric_alignments[,V2:=gsub("\t","",V2),by=V2]
if (num_chimeric_alignments$V2[1]!=0){
  chimeric_alignment_info_extendors = fread(paste(SPLASH_directory,"STAR_alignment/T2T_extendorsChimeric.out.sam",sep=""),header=FALSE,skip="NH:")
  anchors[extendor_index%in%chimeric_alignment_info_extendors$V1,is.STAR_chimeric:=1]
}
anchors[STAR_T2T_CIGAR%like%"N",is.STAR_SJ:=1] # the extendors with non NA STAR alignment are flagged as mapped by STAR
anchors[!is.na(STAR_T2T_chr) & (is.STAR_chimeric==0),is.aligned_STAR:=1] # the extendors with non NA STAR alignment are flagged as mapped by STAR

##########################################################################################
##########################################################################################
##########################################################################################
toc()

tic("add gene name")
########################################################################################
################ adding gene names based on extendor alignments to the genome #########
########################################################################################
anchors[,extendor_gene:=NULL]
system(paste(Samtools_executable, " view -S -b ",SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.sam > ",SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.bam",sep=""))
system(paste(bedtools_executable, " bamtobed -split -i ",SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.bam | sed '/^chr/!d' | sort -k1,1 -k2,2n > ", SPLASH_directory,"STAR_alignment/T2T_called_exons.bed",sep=""))
system(paste(bedtools_executable, " intersect -a ",SPLASH_directory,"STAR_alignment/T2T_called_exons.bed -b ", gene_coords_file, " -wb -loj | cut -f 4,10   | ",bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory,"STAR_alignment/T2T_extendor_genes.txt",sep=""))
extendor_genes = fread(paste(SPLASH_directory,"STAR_alignment/T2T_extendor_genes.txt",sep=""),sep="\t",header=FALSE)
names(extendor_genes) = c("extendor_index","extendor_gene")
anchors = merge(anchors,extendor_genes[!duplicated(extendor_index)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="extendor_index")
anchors[extendor_gene==".",extendor_gene:=NA]
anchors[,num_extendor_gene_anchor:=length(unique(extendor_gene)),by=anchor]  #number of unique extendor gene names for each anchor (is useful for distinguishing TE-like events in which different extendors might map to multiple genes)
########################################################################################
toc()

tic("extract and annotate splice junctions")
##################################################################################################################################
############################ extracting splice junctions for those with reported STAR splice alignment ###########################
##################################################################################################################################
# I first write a sam file that has only the best splice alignment for each extendor
alignment_info_extendors = fread(paste(SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.sam",sep=""),header=FALSE,skip="NH:")
alignment_info_extendors = alignment_info_extendors[!duplicated(V1)][V6%like%"N"]
write.table(alignment_info_extendors,paste(SPLASH_directory,"STAR_alignment/T2T_only_top_splice_alignments.out.sam",sep=""),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
system(paste(Samtools_executable, " view -H ",SPLASH_directory,"STAR_alignment/T2T_extendorsAligned.out.bam > ", SPLASH_directory, "STAR_alignment/sam_header.txt", sep=""))

# now I need to concatenate the header to the new sam file and then convert the sam file to a bam file and then run bamtobed split to get the extracted splice junctions
system(paste("cat ",SPLASH_directory,"STAR_alignment/sam_header.txt ",SPLASH_directory,"STAR_alignment/T2T_only_top_splice_alignments.out.sam > ", SPLASH_directory, "STAR_alignment/T2T_only_top_splice_alignments_with_header.out.sam" ,sep=""))
system(paste(Samtools_executable, " view -S -b ",SPLASH_directory,"STAR_alignment/T2T_only_top_splice_alignments_with_header.out.sam > ",SPLASH_directory,"STAR_alignment/T2T_only_top_splice_alignments_with_header.out.bam",sep=""))
system(paste(bedtools_executable, " bamtobed -split -i ",SPLASH_directory,"STAR_alignment/T2T_only_top_splice_alignments_with_header.out.bam > ",SPLASH_directory,"STAR_alignment/extracted_splice_junction.bed",sep=""))
splice_junctions = fread(paste(SPLASH_directory,"STAR_alignment/extracted_splice_junction.bed",sep=""),header=FALSE) # now I read in the extracted spliec junctions and then by creating a shifted version of them and then appending them to the original data table columnwise I extract spliec junctions as the 2nd coord from the first line and 1st coord from the second line + 1
splice_junctions_shifted = splice_junctions[, data.table::shift(.SD, 1, NA, "lead", TRUE), .SDcols=1:6]
splice_junctions = cbind(splice_junctions,splice_junctions_shifted)
splice_junctions = splice_junctions[V4==V4_lead_1] # to subset to only those with the same extendor_index
splice_junctions$splice_junc=paste(splice_junctions$V1,":",splice_junctions$V3,":",splice_junctions$V2_lead_1+1,sep="")
splice_junctions[,all_splice_juncs:=paste(splice_junc,collapse = "--"),by=V4]
remove(alignment_info_extendors)

#######################################
##### annotating AS splice sites ######
known_splice_sites = fread(annotated_splice_juncs_file)

# now we want to do the same analysis by looking at the splice sites that are involved in at least two distinct junctions
known_splice_sites$chr_V2=paste(known_splice_sites$V1,known_splice_sites$V2,sep=":")
known_splice_sites$chr_V3=paste(known_splice_sites$V1,known_splice_sites$V3,sep=":")
known_splice_sites[,num_uniq_V2_for_V3:=length(unique(chr_V2)),by = chr_V3] # number of partner splice cites for each V2
known_splice_sites[,num_uniq_V3_for_V2:=length(unique(chr_V3)),by = chr_V2] # number of partner splice cites for each V3
alt_v2 = known_splice_sites[num_uniq_V3_for_V2>1]$V2 # the V2 coordinates that have more than one splice site partner
alt_v3 = known_splice_sites[num_uniq_V2_for_V3>1]$V3 # the V3 coordinates that have more than one splice site partner
total = c(alt_v2, alt_v3, alt_v2-1, alt_v3-1, alt_v2+1, alt_v3+1) # I concatenate all coordinates with their +-1 counterparts
splice_junctions[,SSA_AS_annot:=0]
splice_junctions[,SSB_AS_annot:=0]
splice_junctions[V3%in%total,SSA_AS_annot:=1]
splice_junctions[V2_lead_1%in%total,SSB_AS_annot:=1]
splice_junctions$SS_AS_annot=paste(splice_junctions$SSA_AS_annot,":",splice_junctions$SSB_AS_annot,sep="") # the flag that shows whether the 5' and 3' SS of each splice alignment is annotated as being involved in AS
splice_junctions[,all_SS_AS_annot:=paste(SS_AS_annot,collapse = "--"),by=V4]
#######################################
#######################################

##############################################################################
### now check to see if splice site is an annotated exon boundary ############
known_exon_boundaries = fread(annotated_exon_boundaries_file, header = TRUE, sep = "\t")
known_exon_boundaries = known_exon_boundaries[!duplicated(paste(V1,V2,V3))]
total = c(known_exon_boundaries$chr_V2,known_exon_boundaries$chr_V2_1,known_exon_boundaries$chr_V2_2,known_exon_boundaries$chr_V3,known_exon_boundaries$chr_V3_1,known_exon_boundaries$chr_V3_2)
splice_junctions[,SSA_annot:=0]
splice_junctions[,SSB_annot:=0]
splice_junctions[paste(V1,V3,sep="")%in%total,SSA_annot:=1]
splice_junctions[paste(V1,V2_lead_1,sep="")%in%total,SSB_annot:=1]
splice_junctions$SS_annot=paste(splice_junctions$SSA_annot,":",splice_junctions$SSB_annot,sep="") # the flag that shows whether the 5' and 3' SS of each splice alignment is annotated as an exon boundary
splice_junctions[,all_SS_annot:=paste(SS_annot,collapse = "--"),by=V4]
###############################################################################
##############################################################################

anchors[,c("all_splice_juncs","all_SS_AS_annot","all_SS_annot"):=NULL]
splice_junctions = splice_junctions[!duplicated(V4)]
anchors = merge(anchors,splice_junctions[,list(V4,all_splice_juncs,all_SS_AS_annot,all_SS_annot)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="V4")
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
toc()

tic("bowtie alignment of anchors")
#########################################################################################################
######## Bowtie2 alignment of anchors to extendors to remove duplicate anchors ##########################
#########################################################################################################
anchors_fasta = riffle(paste(">",anchors[!duplicated(anchor)]$anchor_index,sep=""),anchors[!duplicated(anchor)]$anchor)
write.table(anchors_fasta,paste(SPLASH_directory,"anchors_fasta.fa",sep = ""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep = "\t")
system(paste("mkdir ",SPLASH_directory,"Bowtie_alignment",sep="")) # this creates the directory for BOWTIE alignment and index files
system(paste(Bowtie2_executable, "-build ", SPLASH_directory,"extendors_fasta.fa ",SPLASH_directory,"Bowtie_alignment/extendors_Bowtie2_index",sep="")) # making the index files for extendor sequences
#aligning anchors to extendors
system(paste(Bowtie2_executable, " -f -k 10 -x ",SPLASH_directory,"Bowtie_alignment/extendors_Bowtie2_index ","-U ", SPLASH_directory,"anchors_fasta.fa -S ",SPLASH_directory,"Bowtie_alignment/anchors_to_extendors_bowtie.sam",sep=""))

Bowtie_alignment_info_extendors = fread(paste(SPLASH_directory,"Bowtie_alignment/anchors_to_extendors_bowtie.sam",sep=""), header=FALSE, fill = TRUE, skip = "AS:",select = c("V1","V3","V12"))
Bowtie_alignment_info_extendors = Bowtie_alignment_info_extendors[V12=="AS:i:0"]
# now I find the corresponding anchor index for all the extendors to which each anchor has been mapped
# and then those that have been mapped to a extendor with higher rank will be flagged as duplicate
Bowtie_alignment_info_extendors[,aligned_extendor_anchor_index:=as.numeric(strsplit(V3,split="_")[[1]][1]),by=V3]
Bowtie_alignment_info_extendors$aligned_extendor_anchor_index = as.numeric(Bowtie_alignment_info_extendors$aligned_extendor_anchor_index)
Bowtie_alignment_info_extendors[,min_aligned_extendor_anchor_index := min(aligned_extendor_anchor_index),by=V1]
Bowtie_alignment_info_extendors[,duplicate_anchor:=0]
Bowtie_alignment_info_extendors[min_aligned_extendor_anchor_index<V1,duplicate_anchor:=1]
Bowtie_alignment_info_extendors[duplicate_anchor==1 & aligned_extendor_anchor_index==min_aligned_extendor_anchor_index,high_ranked_extendors:=paste(V3,collapse = ":"),by=V1]
Bowtie_alignment_info_extendors = Bowtie_alignment_info_extendors[duplicate_anchor==1 & !is.na(high_ranked_extendors)]
Bowtie_alignment_info_extendors = Bowtie_alignment_info_extendors[!duplicated(V1)]
anchors = merge(anchors,Bowtie_alignment_info_extendors[,list(V1,duplicate_anchor,high_ranked_extendors)],all.x=TRUE,all.y=FALSE,by.x="anchor_index",by.y="V1")
anchors[is.na(duplicate_anchor),duplicate_anchor:=0]
##########################################################################################
##########################################################################################
##########################################################################################
toc()

tic("bowtie to univec")
######################################################################################################
################ Bowtie2 alignment of anchors and targets to univec to remove barcodes ###############
#####################################################################################################
system(paste(Bowtie2_executable, " -f -x ", Bowtie2_univec_index, " -U ", SPLASH_directory,"anchors_fasta.fa -p 4 --quiet | sed '/^@/d' | cut -f1,3 | sort | uniq > ",SPLASH_directory,"Bowtie_alignment/anchor_hits_illumina_adapters.tsv",sep=""))
Bowtie_alignment_info_anchors = fread(paste(SPLASH_directory,"Bowtie_alignment/anchor_hits_illumina_adapters.tsv",sep=""), header=FALSE, fill=TRUE)
Bowtie_alignment_info_anchors = Bowtie_alignment_info_anchors[V2!="*"] # finding anchors that have been mapped to a barcode
anchors = anchors[!anchor_index%in%Bowtie_alignment_info_anchors$V1] # removing anchors that have been mapped to a barcode

targets_fasta = riffle(paste(">",anchors$extendor_index,sep=""),anchors$target)
write.table(targets_fasta,paste(SPLASH_directory,"targets_fasta.fa",sep = ""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep = "\t")
system(paste(Bowtie2_executable, " -f -x ", Bowtie2_univec_index, " -U ", SPLASH_directory,"targets_fasta.fa -p 4 --quiet | sed '/^@/d' | cut -f1,3 | sort | uniq > ",SPLASH_directory,"Bowtie_alignment/target_hits_illumina_adapters.tsv",sep=""))

Bowtie_alignment_info_targets = fread(paste(SPLASH_directory,"Bowtie_alignment/target_hits_illumina_adapters.tsv",sep=""), header=FALSE, fill=TRUE)
Bowtie_alignment_info_targets = Bowtie_alignment_info_targets[V2!="*"] # finding targets that have been mapped to a barcode
anchors = anchors[!extendor_index%in%Bowtie_alignment_info_targets$V1] # removing anchors that have at least one target mapped to a barcode

anchors[,num_target_per_anchor:=length(unique(target)),by=anchor] # making sure that only anchors with at least two targets are remaining
anchors = anchors[num_target_per_anchor > 1]  
anchors[,num_target_per_anchor:=NULL]
####################################################################################################
####################################################################################################
####################################################################################################
toc()



tic("soft clip realignment")
#######################################################################################
################## split mapping of those with softclipped bases ######################
#######################################################################################
## I select those for realignment that have at least 10 softclipped bases
anchors[,c("Bowtie_T2T_softclip_flag","Bowtie_T2T_softclip_chr","Bowtie_T2T_softclip_coord","Bowtie_T2T_softclip_CIGAR","Bowtie_T2T_softclip_num_alignments","soft_clipped_bases"):=NULL]
anchors[!is.na(STAR_T2T_CIGAR) & (STAR_T2T_CIGAR%like%"S"),num_STAR_soft_clipped:=max(explodeCigarOpLengths(STAR_T2T_CIGAR, ops=c("S"))[[1]]),by=STAR_T2T_CIGAR]
soft_clipped_extendors_for_realignment = anchors[num_STAR_soft_clipped > 9]
#soft_clipped_extendors_for_realignment[,num_CIGAR_S:=length(explodeCigarOpLengths(STAR_T2T_CIGAR, ops=c("S"))[[1]]),by=STAR_T2T_CIGAR]

soft_clipped_extendors_for_realignment[,first_CIGAR_num:=explodeCigarOpLengths(STAR_T2T_CIGAR)[[1]][1],by=STAR_T2T_CIGAR] # I need this (to know which end of the read softclipping occurs) for extracting softclipped bases
soft_clipped_extendors_for_realignment[first_CIGAR_num==num_STAR_soft_clipped & (STAR_T2T_flag == 0 | STAR_T2T_flag==256),soft_clipped_bases:=substr(extendor,1,num_STAR_soft_clipped)]
soft_clipped_extendors_for_realignment[first_CIGAR_num!=num_STAR_soft_clipped & (STAR_T2T_flag == 16 | STAR_T2T_flag==272),soft_clipped_bases:=substr(extendor,1,num_STAR_soft_clipped)]
soft_clipped_extendors_for_realignment[first_CIGAR_num==num_STAR_soft_clipped & (STAR_T2T_flag == 16 | STAR_T2T_flag==272),soft_clipped_bases:=substr(extendor,nchar(extendor)-num_STAR_soft_clipped+1,nchar(extendor))]
soft_clipped_extendors_for_realignment[first_CIGAR_num!=num_STAR_soft_clipped & (STAR_T2T_flag == 0 | STAR_T2T_flag==256),soft_clipped_bases:=substr(extendor,nchar(extendor)-num_STAR_soft_clipped+1,nchar(extendor))]

soft_clipped_extendors_for_realignment[,extendor_index:=paste(">",anchor_index,"_",extendor_order,sep=""),by=1:nrow(soft_clipped_extendors_for_realignment)] # I use extendor_index as the unique identifier for each extendor sequence for the alignment step
soft_clipped_extendors_for_realignment_fasta = riffle(soft_clipped_extendors_for_realignment$extendor_index,soft_clipped_extendors_for_realignment$soft_clipped_bases)
soft_clipped_extendors_for_realignment_fasta = data.table(soft_clipped_extendors_for_realignment_fasta ) # the fasta file containing all extendor sequences
write.table(soft_clipped_extendors_for_realignment_fasta,paste(SPLASH_directory,"Bowtie_alignment/soft_clipped_extendors_for_realignment_fasta.fa",sep = ""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep = "\t")
system(paste(Bowtie2_executable, " -f -x ", Bowtie2_reference, " -U ", SPLASH_directory,"Bowtie_alignment/soft_clipped_extendors_for_realignment_fasta.fa -S " ,SPLASH_directory,"Bowtie_alignment/soft_clipped_extendors_bowtie.sam --no-unal -k 10 --sam-nohead",sep=""))

alignment_info_extendors = fread(paste(SPLASH_directory,"Bowtie_alignment/soft_clipped_extendors_bowtie.sam",sep=""), header=FALSE, fill=TRUE)
alignment_info_extendors = alignment_info_extendors[V12%like%"AS:i:0"]
alignment_info_extendors[,num_alignments:=.N,by=V1]
alignment_info_extendors = alignment_info_extendors[num_alignments==1]

soft_clipped_extendors_for_realignment[,extendor_index:=gsub(">","",extendor_index),by=extendor_index]
soft_clipped_extendors_for_realignment = merge(soft_clipped_extendors_for_realignment,alignment_info_extendors[!duplicated(V1),list(V1,V2,V3,V4,V6,num_alignments)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="V1")
setnames(soft_clipped_extendors_for_realignment,c("V2","V3","V4","V6","num_alignments"),c("Bowtie_T2T_softclip_flag","Bowtie_T2T_softclip_chr","Bowtie_T2T_softclip_coord","Bowtie_T2T_softclip_CIGAR","Bowtie_T2T_softclip_num_alignments"))


anchors = merge(anchors,soft_clipped_extendors_for_realignment[,list(extendor_index,first_CIGAR_num,soft_clipped_bases,Bowtie_T2T_softclip_flag,Bowtie_T2T_softclip_chr,Bowtie_T2T_softclip_coord,Bowtie_T2T_softclip_CIGAR,Bowtie_T2T_softclip_num_alignments)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="extendor_index")
## below I try to find those softclipped realignment that are evidence for circular RNA
anchors[,is.circle:=0]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==0 | STAR_T2T_flag==256) &(Bowtie_T2T_softclip_flag==0 | Bowtie_T2T_softclip_flag==256) & (first_CIGAR_num==num_STAR_soft_clipped)  & (Bowtie_T2T_softclip_coord > STAR_T2T_coord), is.circle:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==0 | STAR_T2T_flag==256) &(Bowtie_T2T_softclip_flag==0 | Bowtie_T2T_softclip_flag==256) & (first_CIGAR_num!=num_STAR_soft_clipped)  & (Bowtie_T2T_softclip_coord < STAR_T2T_coord), is.circle:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==16 | STAR_T2T_flag==272) &(Bowtie_T2T_softclip_flag==16 | Bowtie_T2T_softclip_flag==272) & (first_CIGAR_num==num_STAR_soft_clipped)  & (Bowtie_T2T_softclip_coord > STAR_T2T_coord), is.circle:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==16 | STAR_T2T_flag==272) &(Bowtie_T2T_softclip_flag==16 | Bowtie_T2T_softclip_flag==272) & (first_CIGAR_num!=num_STAR_soft_clipped)  & (Bowtie_T2T_softclip_coord < STAR_T2T_coord), is.circle:=1 ]
anchors[,is.strandcross:=0]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==0 | STAR_T2T_flag==256) &(Bowtie_T2T_softclip_flag==16 | Bowtie_T2T_softclip_flag==272), is.strandcross:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==16 | STAR_T2T_flag==272) &(Bowtie_T2T_softclip_flag==0 | Bowtie_T2T_softclip_flag==256), is.strandcross:=1 ]
anchors[,is.fusion:=0]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr!=STAR_T2T_chr), is.fusion:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==0 | STAR_T2T_flag==256) &(Bowtie_T2T_softclip_flag==0 | Bowtie_T2T_softclip_flag==256)& abs(Bowtie_T2T_softclip_coord -STAR_T2T_coord)>1000000, is.fusion:=1 ]
anchors[!is.na(soft_clipped_bases) & STAR_T2T_num_alignments==1 & Bowtie_T2T_softclip_num_alignments==1 & (Bowtie_T2T_softclip_chr==STAR_T2T_chr) & (STAR_T2T_flag==16 | STAR_T2T_flag==272) &(Bowtie_T2T_softclip_flag==16 | Bowtie_T2T_softclip_flag==272)& abs(Bowtie_T2T_softclip_coord -STAR_T2T_coord)>1000000, is.fusion:=1 ]
##############################################################################################################
##############################################################################################################
toc()


tic("annotate repeats")
#############################################################################
########################## annotate repeats #################################
#############################################################################
anchors[,repeat_annotation:=NULL]
system(paste(bedtools_executable, " intersect -a ", SPLASH_directory, "STAR_alignment/T2T_called_exons.bed -b ", repeats_annotation_file, " -wb -loj | cut -f 4,10 | sort -nrk1 | ",bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory,"STAR_alignment/T2T_repeat_annotation.txt",sep=""))
repeat_annotation = fread(paste(SPLASH_directory, "STAR_alignment/T2T_repeat_annotation.txt", sep = ""), sep = "\t", header = FALSE)
names(repeat_annotation) = c("extendor_index", "repeat_annotation")
anchors = merge(anchors,repeat_annotation[!duplicated(extendor_index)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="extendor_index")
anchors[,repeat_annotation:=gsub(".,", "", repeat_annotation, fixed = TRUE), by = repeat_annotation]
##############################################################################################################
##############################################################################################################
toc()

tic("assigning classes")
##############################################################################################################################################################
########## below the anchors are classified to splicing or based_pair_change anchors based on their ham_dist, lev_dist and STAR alignment info  ##############
anchors[,anchor_event:=""] # anchor event defines the corresponding event for the anchor
anchors[ham_dist==lev_dist,anchor_event:=paste("Base_pair_change_",ham_dist,sep="")]
anchors[,is_STAR_SJ_in_top_two:=0]
anchors[(extendor_order<3),is_STAR_SJ_in_top_two:=sum(is.STAR_SJ),by=anchor]
anchors[,is_STAR_SJ_in_top_two:=max(is_STAR_SJ_in_top_two),by=anchor]
anchors[,is_aligned_STAR_in_top_two:=0]
anchors[(extendor_order<3),is_aligned_STAR_in_top_two:=sum(is.aligned_STAR),by=anchor]
anchors[,is_aligned_STAR_in_top_two:=max(is_aligned_STAR_in_top_two),by=anchor]
anchors[ is_STAR_SJ_in_top_two>0 & (ham_dist!=lev_dist | ham_dist>5) & num_extendor_gene_anchor==1, anchor_event:="splicing"] # if both top extendors are mapped as SJ by STAR we keep event as splicing
anchors[,vdj:=0]
anchors[((extendor_gene%like%"IGH") |(extendor_gene%like%"IGK") | (extendor_gene%like%"IGV")) & !extendor_gene%like%"PIGK", vdj:=1 ]
anchors[,vdj:=max(vdj),by=anchor]

repeat_anchors = unique(anchors[anchor_event == ""  & repeat_annotation != "." & (extendor_order == 1 | extendor_order == 2)]$anchor_index)
anchors[anchor_index %in% repeat_anchors, anchor_event := "Repeat"] 
################################################################################################################################################################
toc()

write.table(anchors,paste(SPLASH_directory,"classified_anchors.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
