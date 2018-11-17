#!/usr/bin/env python
# count reads to bins with cell barcode

import argparse

parser = argparse.ArgumentParser(description='count reads to bins with bin/cell coordinate')
parser.add_argument('-i', type=str, dest="input", help='input bam file')
parser.add_argument('-g', '--genomeSize', type=str, dest="genomeSize", help='genome size')
parser.add_argument('-b', '--binSize', type=int, dest="binSize", default=5000, help='bin size')
parser.add_argument('-o', type=str, dest="output", help='output prefix')

args = parser.parse_args()

import sys
import pysam
import copy
import numpy as np
from collections import Counter
import pybedtools

def is_sorted_queryname(header):
    """
    Check if bam fiel is sorted by read name.
    """
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                return True
    return False

def bins_maker(genomeSize):
    genomeSizes = pybedtools.BedTool(genomeSize)
    bins = genomeSizes.window_maker(g=genomeSize, w=bin_size)
    bin_dict = {}
    idx = 0
    for i in bins:
        chr = i.chrom
        start = i.start
        end = i.start
        key = (chr+"\t"+str(start)+"\t"+str(end))
        bin_dict[key] = idx
        idx += 1
    return bin_dict

def get_genSize(headerSQ):
    genSize = {}
    for i in range(len(headerSQ)):
        key, values = headerSQ[i]["SN"], headerSQ[i]["LN"]
        genSize[key]= values
    return genSize
    
def main():
    bamf = args.input
    bin_size = int(args.binSize)
    genomeSize = args.genomeSize
    outf = args.output
    
    # start reading the bam
    samfile = pysam.AlignmentFile(input_bam, "rb")
    # check if bam file is sorted by name or not
    if not is_sorted_queryname(samfile.header):
        sys.exit("Error: bam needs to be sorted by read name")
    
    # get chormosome size from header 
    headerSQ = samfile.header["SQ"]
    genSizes = get_genSize(headerSQ)
    bins_dict = bins_maker(genomeSize)
    idx = 0 
    counts = {}
    barcodes = {}
    barcode_tmp = "AAA"

    while True:
        try:
            read1 = next(samfile)
            read2 = next(samfile)
        except: # to the end
            break
        if read1.is_proper_pair:
	    if read1.is_read1:
                pass
            else:
                tmp = copy.deepcopy(read1)
                read1 = read2
                read2 = tmp
            ref1 = read1.reference_name
            ref2 = read2.reference_name
            barcode = read1.qname.strip().split(":")[0]
            if barcode not in barcodes:
                barcodes[barcode] = idx
                idx += 1
            length = read1.reference_length
            if not read1.is_reverse:
                pos = read1.reference_start
            else:
                pos = read2.reference_start
            chrSize1 = int(genSizes[ref1])
            chrSize2 = int(genSizes[ref2])
            if(ref1 == ref2 and abs(read1.reference_start//bin_size - read1.next_reference_start//bin_size) <= 1):
                if (pos//bin_size + 1) * bin_size < chrSize1:
                    key = (ref1+"\t"+str((pos//bin_size) * bin_size)+"\t"+str((pos//bin_size + 1) * bin_size)+"\t"+barcode)
                    if key in counts:
                        counts[key] += 1
                    else:
                        counts[key] = 1
                else:
                    key = (ref1+"\t"+str((pos//bin_size) * bin_size)+"\t"+str(chrSize1)+"\t"+barcode)
                    if key in counts:
                        counts[key] += 1
                    else:
                        counts[key] = 1
            if(ref1 == ref2 and abs(read1.reference_start//bin_size - read1.next_reference_start//bin_size) > 1):
                midpos = pos + length * 0.5
                if (midpos//bin_size + 1) * bin_size < chrSize1:
                   key = (ref1+"\t"+str((midpos//bin_size) * bin_size)+"\t"+str((midpos//bin_size + 1) * bin_size)+"\t"+barcode)
                    if key in counts:
                        counts[key] += 1
                    else:
                        counts[key] = 1
                else:
                    key = (ref1+"\t"+str((midpos//bin_size) * bin_size)+"\t"+str(chrSize1)+"\t"+barcode)
                    if key in counts:
                        counts[key] += 1
                    else:
                        counts[key] = 1

    samfile.close()
    outmat = ".".join(outf,"mat")
    outxgi = ".".join(outf,"xgi")
    np.savetxt(outxgi, barcodes)
    for key,val in counts.items():
        items = key.strip().split("\t")
        xgi = items[3]
        ygi = items[0]+"\t"+items[1]+"\t"+items[2]
        xgi_idx = barcode[xgi]
        ygi_idx = bins_dict[ygi]
        outmat.write(str(xgi_idx) + "\t" + str(ygi_idx) + "\t" + str(val) + "\n")
    outmat.close()
 
if __name__ == "__main__":
    main()

