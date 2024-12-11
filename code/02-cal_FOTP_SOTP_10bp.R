###calculate FOTP and SOTP (after 10bp sequence extraction)
library(stringr)
library(data.table)
setwd("fragment_detail_path")  ###The file path for 10bp sequences information (extraction from bam file)
sapply(c("FOTP_2mer","SOTP_3mer","single_base_freq"),dir.create)
heal_file <- list.files(".",pattern="detail")
CN <- unlist(str_split_fixed(heal_file,"_fragment_",2))[,1]


for (i in c(1:length(heal_file))){    
    heal_sample <- fread(paste0("./",heal_file[i]))
    heal_sample <- na.omit(heal_sample)                

    matrix_2mer_trans <- as.numeric()
    table_single_base <- prop.table(table(c(str_sub(heal_sample$read1,1,1),str_sub(heal_sample$read2,1,1))))
    table_initial <- table(c(str_sub(heal_sample$read1,1,1),str_sub(heal_sample$read2,1,1)))
    m=1
    while(m<=49){
        matrix_2mer_T_name <- matrix(data = rep(0, 16), nrow = 16, ncol = 1)
        rownames(matrix_2mer_T_name) <- c("AA", "AC", "AG", "AT","CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT","TA", "TC", "TG", "TT")
        
        table_m_base <- table(c(str_sub(heal_sample$read1,m+1,m+1),str_sub(heal_sample$read2,m+1,m+1)))
        table_2mer_base <- table(c(str_sub(heal_sample$read1,m,m+1),str_sub(heal_sample$read2,m,m+1)))
        m_num_2mer <- match(names(table_2mer_base), rownames(matrix_2mer_T_name))
        matrix_2mer_T_name[m_num_2mer,1] <- table_2mer_base
        matrix_2mer_trans <- c(matrix_2mer_trans,matrix_2mer_T_name[1:4,1]/sum(matrix_2mer_T_name[1:4,1]),matrix_2mer_T_name[5:8,1]/sum(matrix_2mer_T_name[5:8,1]),
                        matrix_2mer_T_name[9:12,1]/sum(matrix_2mer_T_name[9:12,1]),matrix_2mer_T_name[13:16,1]/sum(matrix_2mer_T_name[13:16,1]))
        table_initial <- c(table_initial,table_m_base)
        m = m+1
    }
    #print(dim(table_single_base))
    print(length(matrix_2mer_trans))
    matrix_2mer_trans <- c(table_single_base,matrix_2mer_trans)
    write.table(table_initial,file=paste0("./single_base_freq/",CN[i],".txt"))
    write.table(matrix_2mer_trans,file=paste0("./FOTP_2mer/",CN[i],".txt")) ###save result of FOTP
	
	
    matrix_3mer_trans <- as.numeric()
	table_12_p <- as.numeric(prop.table(table(c(str_sub(heal_sample$read1,1,2),str_sub(heal_sample$read2,1,2)))))
    j=1
    while(j<=8){
        matrix_3mer <- matrix(data = rep(0, 64), nrow = 64, ncol = 1)
        RN_12_3_64D <- read.table("/FOTP/code/trans_64D_name.txt")  #rownames of SOTP (See the attached file "trans_64D_name.txt" for details)
        rownames(matrix_3mer) <- RN_12_3_64D$V1
        table_3mer_base <- table(c(str_sub(heal_sample$read1,j,j+2),str_sub(heal_sample$read2,j,j+2)))

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
    write.table(matrix_3mer_trans,file=paste0("./SOTP_3mer/",CN[i],".txt"))  ###save results of SOTP    
}
