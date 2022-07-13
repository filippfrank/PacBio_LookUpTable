#!/usr/bin/python
from Bio.Seq import Seq
#import csv 
import pandas as pd
import matplotlib.pyplot as plt
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
import numpy as np
import os

os.chdir('/Users/filipp/OneDrive - Emory University/Documents - Ortlund Lab - RADx VTF/RADx VTF/HTS/LibraryPrep/PacBio/N/Libraries_04_05/Lib04/')

def process(lines):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}
    
fn = "N-Lib04.hifi_reads.fastq"


n = 4
reads = {}

read_lengths = []

with open(fn, 'r') as fh:
    lines= []
    for line in fh:
        lines.append(line.strip('\n'))
        if len(lines) == n:
            record = process(lines)
            #sys.stderr.write("Record: %s\n" % (str(record)))
            #print(record)
            lines = []
            read_lengths.append(len(record["sequence"]))
            reads[record["name"]] = record

number_of_reads = 0
number_of_reads_with_anchor = 0
number_of_reads_with_stop_plus_one = 0 
number_of_reads_with_start = 0
number_of_reads_with_correct_length = 0
number_of_reads_with_stop_minus_one = 0
number_of_16ntbc_reads_with_correct_length = 0
barcodes = {}

#Define some pieces of data to extract
seq_lengths = []
barcode_seqs = []
DNA_seqs = []
aa_seqs = []
read_IDs = []
Mut_Name = []
AA_WT = []
AA_Mut = []
AA_number = []
number_of_mutations = []


reference_aa = "SDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQVILLNKHIDAYKTFPPTEPKKDKKKKADETQALPQRQKKQQTVTLLPAADLDDFSKQLQQSMSSADSTQA"

#for read in list(reads.keys())[0:100]:
for read in reads.keys():
    number_of_reads += 1
    seq_ = Seq(reads[read]["sequence"])
    #Search for expected 3' end sequence
    barcode_start = (seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA") +
                         len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
    #print(barcode_start)
    if barcode_start == 44:
        #Search for 3' end sequence with AAGtoAAA mutation (both code for lysine)
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAACCACGTTAA")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #Search for TAG stop codon and AAG
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAG")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #Search for TAG stop codon and AAA
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAACCACGTTAG")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #If the expected 3' sequence was still not found, the sequence may be the reverse
        seq_ = seq_.reverse_complement()
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #Search the reverse complement for AAGtoAAA mutation
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAACCACGTTAA")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #Search the reverse complement for TAG stop codon
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAG")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start == 44:
        #Search the reverse complement for TAG stop codon and AAGtoAAA
        barcode_start = (int(
            seq_.find("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAACCACGTTAG")) +
            len("TCCCTTATTATTCTCATCATGCTGTGGCAGAAGAAGCCACGTTAA"))
        #print(barcode_start)
    if barcode_start != 44:
        #print(barcode_start)
        number_of_reads_with_anchor += 1
        #Get the barcode
        barcode = seq_[barcode_start:barcode_start+15]
        barcodes[barcode] = True
        #Find the starting sequence
        if len(barcode) == 15:
            is_5prime_present = (seq_.find(
                                "ATGGAGTTCGGGCTCAGCTGGGTGTTTTTGGTAGCTCTCTTTAGAGGAGTGCAGTGCGGAGGGTCCGAGCAGAAACTGATATCAGAAGAGGATCTTGGGGGGTCAGGAGGGAGC"))
            start_of_n = (seq_.find(
                                "ATGGAGTTCGGGCTCAGCTGGGTGTTTTTGGTAGCTCTCTTTAGAGGAGTGCAGTGCGGAGGGTCCGAGCAGAAACTGATATCAGAAGAGGATCTTGGGGGGTCAGGAGGGAGC") +
                                      len("ATGGAGTTCGGGCTCAGCTGGGTGTTTTTGGTAGCTCTCTTTAGAGGAGTGCAGTGCGGAGGGTCCGAGCAGAAACTGATATCAGAAGAGGATCTTGGGGGGTCAGGAGGGAGC"))
            N_seq_ = Seq(seq_[start_of_n:barcode_start-195])
            if is_5prime_present != -1:
                #print(start_of_n)
                #Count reads with correct start (IgG4, Myc) sequence
                number_of_reads_with_start += 1
                seq_lengths.append(len(N_seq_))
                if len(N_seq_) == 1254:
                    number_of_reads_with_correct_length += 1
                    mut_ = []
                    AA_WT_tmp = []
                    AA_Mut_tmp = []
                    A_number_tmp = []
                    Mut_Name_tmp = []
                    for i, a in enumerate(N_seq_.translate()):
                        if a != reference_aa[i]:
                            mut_.append(1)
                            AA_WT_tmp = reference_aa[i]
                            AA_number_tmp = i+2
                            AA_Mut_tmp = a
                            Mut_Name_tmp = "%c%i%c" % (reference_aa[i], i+2, a) 
                            #print(a)
                    number_of_mutations.append(len(mut_))
                    if len(mut_) == 1:
                        AA_WT.append(AA_WT_tmp)
                        AA_number.append(AA_number_tmp)
                        AA_Mut.append(AA_Mut_tmp)
                        Mut_Name.append(Mut_Name_tmp)
                        #Extract data and sequences
                        barcode_seqs.append(barcode)
                        #DNA_seqs.append(N_seq_)
                        aa_seqs.append(N_seq_.translate())
                        read_IDs.append(list(read.split("/"))[1])

####################################
#Combine all data
####################################
d = {'Names': Mut_Name, 
     'Barcode': barcode_seqs, 
     'AA_WT': AA_WT,
     'AA_Number': AA_number,
     'AA_Mut': AA_Mut,
     #'Number of Mutations': number_of_mutations,
     #'Sequence Length': seq_lengths, 
     #'DNA': DNA_seqs, 
     #'AA_Seq': aa_seqs
     }
#print(d)
df = pd.DataFrame(d)
#Save to a .csv file
df.to_csv('data.csv',index=False)
#print(df)

#############################
###Map Mutations to Barcodes
#############################

from collections import defaultdict

barcodes = {}

with open("data.csv", "rU") as f: #csv; AA name, barcode, wt AA, AA number, mut AA
	next(f)
	for line in f:
		name_ = line.split(",")[0]
		barcode_ = line.split(",")[1]
		if barcodes.get(barcode_): # if we've encountered this barcode before
			if barcodes[barcode_].get(name_): # if we've seen this mutation before
				barcodes[barcode_][name_] += 1
			else:
				barcodes[barcode_][name_] = 1
		else:
			barcodes[barcode_] = defaultdict(int)
			barcodes[barcode_][name_] = 1
		
with open("mapped_barcodes_out.txt", "w") as f:
	f.write("barcode\tnumber of unique mutations\tmutations\tdepth\n")
	for key, value in barcodes.items():
		f.write("\t".join([key,str(len(value)),",".join(value.keys()),",".join([str(x) for x in value.values()])]) + "\n")










