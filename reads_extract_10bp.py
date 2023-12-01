import os
import re
import argparse
import pysam
from collections import defaultdict as dd


def get_region(rf):
    rd = {}
    for line in open(rf, "r"):
        fields = line.rstrip().split()
        rd[fields[0], fields[1], fields[2]] = 1
    return rd


def fragment_features(rd, bf, suffix):
    name = os.path.basename(bf).split(".")[0]

    #bam = pysam.AlignmentFile(bf, "rb")
    bam = pysam.AlignmentFile(bf, "rb", ignore_truncation=True)
    read_dict = dd(list)
    for i in rd.keys():
        reads = bam.fetch(i[0], int(i[1]), int(i[2]))
        for read in reads:
            if read.is_proper_pair:
                if read.mapq <= 30:
                    continue
                if re.search("S|H", read.cigarstring):
                    continue
                if read.is_qcfail:
                    continue
                if read.is_duplicate:
                    continue
                if read.is_secondary:
                    continue
                if read.is_unmapped:
                    continue
                if read.is_supplementary:
                    continue
                end_motif = ""
                if read.flag == 99 or read.flag == 163:          # forward strand
                    end_motif = read.get_reference_sequence().upper()[:10]
                elif read.flag == 83 or read.flag == 147:      # reverse strand
                    end_motif = read.get_reference_sequence().upper()[-10:][::-1]
                    end_motif = end_motif.translate(str.maketrans('ATCG', 'TAGC'))
                read_dict[i[0], i[1], i[2], read.query_name].append((read.pos+1, abs(read.tlen), end_motif))

    of_name1 = name + "_fragment_feature." + suffix + ".detail"
    #of_name2 = name + "_end_motif." + suffix + ".summary"
    of1 = open(of_name1, "w")
    #of2 = open(of_name2, "w")

    print("read_name", "read_chr", "read_start", "read_end", "read1", "read2", "fragment_size",
          "region_chr", "region_start", "region_end", sep="\t", file=of1)

    motif_dict = dd(int)
    motif_num = 0

    for k in read_dict.keys():
        if len(read_dict[k]) != 2:
            continue

        mate1 = read_dict[k][0]
        mate2 = read_dict[k][1]
        distance = mate1[0] - mate2[0]
        motif_dict[mate1[2]] += 1
        motif_dict[mate2[2]] += 1
        motif_num += 2
        if distance < 0:
            start = mate1[0]
        else:
            start = mate2[0]
        end = start + mate1[1]
        print(k[3], k[0], start, end, mate1[2], mate2[2], mate1[1], k[0], k[1], k[2], sep="\t", file=of1)

    #for m in motif_dict.keys():
        #motif_ratio = motif_dict[m] / motif_num
        #print(m, motif_ratio, sep="\t", file=of2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, dest="input",
                        help='Input bam file.')
    parser.add_argument("-v", "--interval", type=str, required=True, dest="interval",
                        help="Region used for extracting reads.")
    parser.add_argument("-s", "--suffix", required=True, type=str, dest="suffix",
                        help='suffix of output files')
    args = parser.parse_args()

    regions = args.interval.split(",")

    region_dict = get_region(args.interval)
    fragment_features(region_dict, args.input, args.suffix)
