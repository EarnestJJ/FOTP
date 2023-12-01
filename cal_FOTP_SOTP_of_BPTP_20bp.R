#########提取左4右9的fragment的FOTP########

library(data.table)
library(tidyr)
library(GenomicRanges)
library(SummarizedExperiment)
library(stringr)
setwd("fragment_path") ###The path to extract 20bp breakpoint sequences information
sapply(c("FOTP_20BPTP","SOTP_20BPTP","single_base_freq_20BPTP"),dir.create)

heal_file <- list.files(path = ".", recursive = F, pattern = "fragment")
CN <- unlist(str_split_fixed(heal_file,"[.]",2))[,1]
file_path <- paste0("full_path", heal_file)
bins_5MB <- read.table("./bins_5MB_clean.bed")
colnames(bins_5MB) <- c("chrom","start","end")
bins_5MB.gr <- makeGRangesFromDataFrame(bins_5MB,keep.extra.columns = T)


for (i in c(2:length(heal_file))){    
    heal_sample <- fread(paste0("full_path",heal_file[i]))
    heal_sample <- heal_sample[,c(1,2,9,5,10,6,11)]
    colnames(heal_sample) <- c("chrom","start","end","EDM1","EDM2","BPM1","BPM2")               

    wrong <- which(heal_sample$start > heal_sample$end)    #Partial inversion of chromosome coordinates needs to be corrected
    new_end <- heal_sample[wrong,]$start
    new_start <- heal_sample[wrong,]$end
    heal_sample[wrong,]$start <- new_start
    heal_sample[wrong,]$end <- new_end

    heal_sample.gr <- makeGRangesFromDataFrame(heal_sample,keep.extra.columns=T)
    heal_sample <- as.data.frame(subsetByOverlaps(heal_sample.gr,bins_5MB.gr))

    matrix_2mer_trans <- as.numeric()
    table_single_base <- prop.table(table(c(str_sub(heal_sample$BPM1,17,17),str_sub(heal_sample$BPM2,17,17))))
    table_initial <- table(c(str_sub(heal_sample$BPM1,1,1),str_sub(heal_sample$BPM2,1,1)))
    k=1
    while(k<=39){
        table_m_base <- table(c(str_sub(heal_sample$BPM1,k+1,k+1),str_sub(heal_sample$BPM2,k+1,k+1)))
        table_initial <- c(table_initial,table_m_base)
        k=k+1
    }
    table_initial_matrix <- matrix(table_initial,nrow=4,ncol=40)
    row.names(table_initial_matrix) <- c("A","C","G","T")
    colnames(table_initial_matrix) <- paste0(rep("pos",40),c(1:40))    
    write.table(table_initial_matrix,file=paste0("./single_base_freq_20BPTP/",CN[i],".txt"))    
    
    m=17
    while(m<=28){
        matrix_2mer_T_name <- matrix(data = rep(0, 16), nrow = 16, ncol = 1)
        rownames(matrix_2mer_T_name) <- c("AA", "AC", "AG", "AT","CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT","TA", "TC", "TG", "TT")
        
        #table_m_base <- table(c(str_sub(heal_sample$BPM1,k+1,k+1),str_sub(heal_sample$BPM2,k+1,k+1)))
        table_2mer_base <- table(c(str_sub(heal_sample$BPM1,m,m+1),str_sub(heal_sample$BPM2,m,m+1)))
        m_num_2mer <- match(names(table_2mer_base), rownames(matrix_2mer_T_name))
        matrix_2mer_T_name[m_num_2mer,1] <- table_2mer_base
        matrix_2mer_trans <- c(matrix_2mer_trans,matrix_2mer_T_name[1:4,1]/sum(matrix_2mer_T_name[1:4,1]),matrix_2mer_T_name[5:8,1]/sum(matrix_2mer_T_name[5:8,1]),
                        matrix_2mer_T_name[9:12,1]/sum(matrix_2mer_T_name[9:12,1]),matrix_2mer_T_name[13:16,1]/sum(matrix_2mer_T_name[13:16,1]))
        #table_initial <- c(table_initial,table_m_base)
        m = m+1
    }
    #print(dim(table_single_base))
    print(length(matrix_2mer_trans))
    matrix_2mer_trans <- c(table_single_base,matrix_2mer_trans)   
    write.table(matrix_2mer_trans,file=paste0("./FOTP_20BPTP/",CN[i],".txt"))
		
    matrix_3mer_trans <- as.numeric()
	table_12_p <- as.numeric(prop.table(table(c(str_sub(heal_sample$BPM1,17,18),str_sub(heal_sample$BPM2,17,18)))))
    j=17
    while(j<=27){
        matrix_3mer <- matrix(data = rep(0, 64), nrow = 64, ncol = 1)
        RN_12_3_64D <- read.table("./trans_12_3_64D_name.txt") #rownames of SOTP (See the attached file "trans_12_3_64D_name.txt" for details)
        rownames(matrix_3mer) <- RN_12_3_64D$V1
        table_3mer_base <- table(c(str_sub(heal_sample$BPM1,j,j+2),str_sub(heal_sample$BPM2,j,j+2)))

        m_num_3mer <- match(names(table_3mer_base), rownames(matrix_3mer))
        matrix_3mer[m_num_3mer,1] <- table_3mer_base
        revised_trans <- c(matrix_3mer[1:4,1]/sum(matrix_3mer[1:4,1]),matrix_3mer[5:8,1]/sum(matrix_3mer[5:8,1]),matrix_3mer[9:12,1]/sum(matrix_3mer[9:12,1]),
                            matrix_3mer[13:16,1]/sum(matrix_3mer[13:16,1]),matrix_3mer[17:20,1]/sum(matrix_3mer[17:20,1]),matrix_3mer[21:24,1]/sum(matrix_3mer[21:24,1]),
                            matrix_3mer[25:28,1]/sum(matrix_3mer[25:28,1]),matrix_3mer[29:32,1]/sum(matrix_3mer[29:32,1]),matrix_3mer[33:36,1]/sum(matrix_3mer[33:36,1]),
                            matrix_3mer[37:40,1]/sum(matrix_3mer[37:40,1]),matrix_3mer[41:44,1]/sum(matrix_3mer[41:44,1]),matrix_3mer[45:48,1]/sum(matrix_3mer[45:48,1]),
                            matrix_3mer[49:52,1]/sum(matrix_3mer[49:52,1]),matrix_3mer[53:56,1]/sum(matrix_3mer[53:56,1]),matrix_3mer[57:60,1]/sum(matrix_3mer[57:60,1]), 
                            matrix_3mer[61:64,1]/sum(matrix_3mer[61:64,1]))


        matrix_3mer_trans <- c(matrix_3mer_trans,revised_trans)
        j = j+1
    }
    print(length(matrix_3mer_trans))
	matrix_3mer_trans <- c(table_12_p,matrix_3mer_trans)
    write.table(matrix_3mer_trans,file=paste0("./SOTP_BPTP/",CN[i],".txt"))    
}
