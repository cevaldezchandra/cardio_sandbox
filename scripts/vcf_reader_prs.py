### Goal: to read in a VCF file into a pandas frame,
### get a subset of rsIDs for analysis, and send those
### vectors to the risk prediction program

from pysam import VariantFile
import numpy as np
import sys
import os
import re
import pandas as pd

# Goal is to take a patient's VCF and to read
# a set of rsIDs from it
# Iterator that returns the next line in
# genotype counts form from the vcf file
# skips multi-allelic snps
def flip_alleles(records, to_flip):
    """
    Input: list of genotype records, in format:
    1. [[rsid, chrom, pos, ref, alt, array of geno counts]]
    2. name of rsids to flip
    """
    for i, k in enumerate(records):
        if k[0] in to_flip:
            records[i] = flip_one(k)
    return records

def flip_one(k):
    # Flips single record
    new_ref = k[4]
    new_alt = k[3]
    new_counts = -1*k[5] + 2.0
    return [k[0], k[1], k[2], new_ref, new_alt, new_counts]

class genoVCFreader:
    def __init__(self, vcfname, sampleids=None, rsids=None):
        self.vcfhandle = VariantFile(vcfname)
        self.rsids = rsids
        if sampleids == None:
            self.samples = self.vcfhandle.header.samples
        else:
            self.samples = [s for s in self.vcfhandle.header.samples if s in sampleids]


    def __iter__(self):
        return self

    def get_sample_list(self):
        # Returns a list of the samples
        return [s for s in self.samples]

    def next(self):
        """
        Returns a list of the following elements for the next included rsID:
        1. rsid
        2. chrom
        3. position
        4. reference allele
        5. alternate allele
        6. genotype array in counts for requested sample id
        """
        # Read the next line in the vcf file
        rec = self.vcfhandle.next()
        # Go to next with id in rsid list
        if self.rsids:
            while (not rec.id in self.rsids) or len(rec.alts)>1:
                rec = self.vcfhandle.next()
        gts = [ s['GT'] for i,s in enumerate(rec.samples.values()) if rec.header.samples[i] in self.samples]
        gts = np.array(gts, dtype=np.float)
        gt_sums = np.sum(gts, axis=1)
        print rec, gts
        return [rec.id, rec.chrom, rec.pos, rec.ref, rec.alts[0], gt_sums]
