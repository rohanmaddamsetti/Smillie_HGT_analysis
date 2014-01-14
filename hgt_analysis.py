#!/usr/bin/python

##hgt_analysis.py by Rohan Maddamsetti.

from subprocess import call
from os import listdir, environ
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

##put the path from bash into the path in this script.
mypath = environ['PATH']
sys.path.append(mypath)

def makeBLASTinputs(fasta_db):
    hgt_db = open(fasta_db)
    header = hgt_db.readline()
    seq = hgt_db.readline()
    while (header and seq):
        filename = header.replace(" ", "__")[1:].strip()
        cur_output = open("input/"+filename+".fasta",'w')
        cur_output.write(header)
        cur_output.write(seq)
        cur_output.close()
        header = hgt_db.readline()
        seq = hgt_db.readline()

def runBLAST():
    queries = [x for x in listdir("input") if x.endswith(".fasta")]
    for cur_query in queries:
        prefix = cur_query.split('.fasta')[0]
        blast_cmd = ["blastn", "-task", "megablast", "-db", "REL606blastdb",
                     "-outfmt", "5", "-max_target_seqs", "1",
                     "-query", "input/"+cur_query, "-out",
                     "output/"+prefix+".xml"]
        blast_cmd_string = ' '.join(blast_cmd)
        call(blast_cmd_string, shell=True)

def parseBLASTresults(grep_patterns=None):
    '''grep_patterns is a list of patterns used to select subsets of the 
	hgt.fst sequences
    based on matching a pattern in the name. This is not a regex, but a
    substring in the current implementation.'''
    HGT_hits = {}
    ##preprocess the reference genome.
    gene_starts = {}
    gene_ends = {}
    ref_genome = SeqIO.read("../REL606.gbk", "genbank")
    for feature in ref_genome.features:
        if feature.type != "CDS":
            continue
        else:
            cur_locus_tag = feature.qualifiers['locus_tag'][0]
            cur_start = feature.location.start
            cur_end = feature.location.end
            ##print cur_start, cur_end, cur_locus_tag
            gene_starts[cur_locus_tag] = cur_start
            gene_ends[cur_locus_tag] = cur_end

    all_results = [x for x in listdir("output") if x.endswith(".xml")]
    results = []
    if grep_patterns:
        for pattern in grep_patterns:
            grepped_results = [x for x in all_results if pattern in x]
            results = results + grepped_results
    for cur_blast_output in results:
        result_handle = open("output/" + cur_blast_output)
        blast_record = NCBIXML.parse(result_handle).next()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.0000000001:
                    e_value = hsp.expect
                    subject_start = hsp.sbjct_start
                    subject_end = hsp.sbjct_end
                    ##enforce that subject_start < subject_end.
                    if subject_end < subject_start:
                        temp = subject_end
                        subject_end = subject_start
                        subject_start = temp
                    subject_start = ExactPosition(subject_start)
                    subject_end = ExactPosition(subject_end)
                    ##print e_value, subject_start, subject_end
                    ##extract the locus_tags of all genes present in this range.
                    starts_in_range = [x for x in gene_starts.keys() if gene_starts[x] > subject_start and gene_starts[x] < subject_end]
                    ends_in_range = [y for y in gene_ends.keys() if gene_ends[y] > subject_start and gene_ends[y] < subject_end]
                    loci_in_range = [z for z in starts_in_range if z in ends_in_range]
                    ##print loci_in_range
                    for locus in loci_in_range:
                        try:
                            HGT_hits[locus] = HGT_hits[locus] + 1
                        except KeyError:
                            HGT_hits[locus] = 1
    ##Now print the results to a csv file.
    ##pipe to output as so: 'python hgt_analysis.py > out.csv
    print "locus_tag,HGT.hits"
    for tag in sorted(HGT_hits.keys()):
        print tag + "," + str(HGT_hits[tag])

def main():
    
    ## The following two lines are commented out because they only need
    ## to be run once for the analysis.
    ##makeBLASTinputs("hgt_seqs/hgt.fst")
    ##runBLAST()
    parseBLASTresults(grep_patterns=["Escherichia", "Shigella"])

main()
