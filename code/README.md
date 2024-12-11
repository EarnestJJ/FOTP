# code
Detailed description of the relevant code

#step1: extract 10bp end sequences in target regions
python3 01-reads_extract_10bp.py -v bins_5MB_clean.bed -i bam_file -s txt

#step2: calculate FOTP and SOTP feature of 10 bp end sequences
Rscript 02-cal_FOTP_SOTP_10bp.R

#step3: extract 20bp breakpoint suquences  (including other features)
python 03-fragment_BPTP20bp.py -b bam_file -f hg19.fa -o bam_file.fragment

#step4: calculate FOTP feature of 20bp breakpoint transition probability (BPTP)
Rscript 04-cal_FOTP_SOTP_of_BPTP_20bp.R
