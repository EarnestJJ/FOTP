# FOTP
The First-order Transition Probability (FOTP) based on the transition probability of adjacent nucleotides in a N bp length at each 5' end of cfDNA fragments.

#step1: extract 10bp end sequences 
python3 reads_extract_10bp.py -v bins_5MB_clean.bed -i bam_file -s txt

#step2: calculate FOTP and SOTP feature of 10 bp end sequences
Rscript cal_FOTP_SOTP_10bp.R

#step3: extract 20bp breakpoint suquences  (including other features)
python fragment_BPTP20bp.py -b bam_file -f hg19.fa -o bam_file.fragment

#step4: calculate FOTP feature of 20bp breakpoint transition probability (BPTP)
Rscript cal_FOTP_SOTP_of_BPTP_20bp.R
