#!/bin/python

from sys import argv
import argparse
import re
from subprocess import call
import sys
import gffutils
from collections import defaultdict

# These lines assign the arguments given at command line to variables that can be used within the script. You can type 'python proteinFromGFF.py -h' at the command line to see this.
parser = argparse.ArgumentParser(description='Parse a gff file (supplied by augustus) and retrieve the protein sequence prediction from that file.')
parser.add_argument('inputFile',
                    help='The gff file you wish to get protein sequence from.')
parser.add_argument('--outputDirectory', default=".",
                    help='The directory in which to save the extracted sequence. Default is current working directory.')
parser.add_argument('--geneNameOn', dest='geneNameOn', default=False,
                    help='Toggle flag telling the program whether the strain name is preceeded by a gene name in the input filename. Default is False.')
parser.add_argument('--fullGenome', dest='fullGenome', default=False,
                    help='Toggle flag telling the program whether the gff file contains a whole genome annotation. Default is False.')
parser.add_argument('--multiSeq', dest='multiSeq', default=False,
                    help='Toggle flag telling the program whether the gff file contains multiple sequences. Default is False.')

args = parser.parse_args()
inputFile = args.inputFile
outputDirectory = args.outputDirectory
geneNameOn = args.geneNameOn
fullGenome = args.fullGenome
multiSeq = args.multiSeq

# This bit is specific to my file structure to a certain extent. It assumes that the input file will be named something like 'NCYC1_contig123_pos456to789+.fasta.gff' (because I have annotated a fasta file containing a hmmer match).
# I intend to make this part a bit more flexible. Here I just take the last string in the filepath (the filename) and split it up until I have the various bits of location information from the filename as separate variables.
splitInput = inputFile.split("/")
fileName = splitInput[-1]
if "*" in fileName:
    sys.exit(0)

if fullGenome == "True" or fullGenome == "true":
    geneName = ""
    strainName = fileName.split(".")[0]
    contig = ""
    start = ""
    end = ""
    #print "FULL GENOME"
    #db = gffutils.create_db(inputFile, dbfn='getProtein.db', force=True, merge_strategy='merge')
elif geneNameOn == "True" or geneNameOn == "true":
    splitFilename = re.split("[_]", fileName)
    geneName = splitFilename[0] + "_"
    strainName = splitFilename[1]
    contig = re.split("[a-z]+", splitFilename[2])[1]
    start = re.split("[.+a-z]+", splitFilename[3])[1]
    end = re.split("[.+a-z]+", splitFilename[3])[2]
    #print splitFilename
elif multiSeq == "True" or multiSeq == "true":
    splitFilename = re.split("[_]", fileName)
    geneName = ""
    strainName = splitFilename[0]
    contig = ""
    start = ""
    end = ""
else:
    splitFilename = re.split("[._+a-z]+", fileName)
    geneName = ""
    strainName = splitFilename[0]
    contig = splitFilename[1]
    start = splitFilename[2]
    end = splitFilename[3]
#print strainName + "_" + contig + ":" + start + "-" + end

# open the input file (for reading by default)
fh = open(inputFile)

labelLine = False
geneSeq = False
genes = []
genes2 = defaultdict(list)
prot = ""

# This 'for-loop' runs through each line in the file and checks whether it is a gene sequence (it toggles the 'geneSeq' flag to be on if we hit a line saying 'protein sequence' and off if we hit one saying 'end gene').
# The loop takes everything in between and then filters out the actual sequence (by getting rid of whitespaces, # symbols, brackets, and lower case letters). It then saves each sequence 'prot' to a dictionary 'genes'.
for line in fh:
    if line.startswith("# start gene"):
        labelLine = True
        geneNo = re.split("[ \n]", line)[3]
        prot = ""
        continue
    if labelLine == True:
        labels = line.split("\t")
        #print labels
        contig = labels[0]
        start = labels[3]
        end = labels[4]
        labelLine = False
    if line.startswith("# end gene"):
        geneSeq = False
        prot2 = re.findall("[A-Z]+", prot)
        prot3 = ""
        for seq in prot2:
            prot3 = prot3 + seq
        #print prot3
        #genes[geneNo] = prot3
        lst = [geneNo, contig, start, end, prot3]
        genes.append(lst)
    if line.startswith("# protein sequence"):
        geneSeq = True
    if geneSeq == True:
        prot = prot + line
        #print prot

fh.close()
#print genes

for i, j, k, l, m in genes:
    genes2[i].append([j, k, l, m])

# if one doesn't exist, make a directory to put the sequences into (if there are any genes to store).
if genes2 and multiSeq:
    command1 = "mkdir -p " + outputDirectory
    #print command1
    call(command1, shell=True)
elif genes2:
    command1 = "mkdir -p " + outputDirectory + "/" + strainName
    #print command1
    call(command1, shell=True)

# This last loop iterates through the 'genes' dictionary and constructs a fasta file containing all the location information from before, plus the gene number (g1, g2, g3, etc - from the gff annotation), and the protein sequence.
# The file is saved to the directory created above (i.e. saved by strain), in the output directory specified at the command line.
#for gene, sequence in genes.iteritems():
#    #print gene, sequence
#    outputFile = outputDirectory + "/" + strainName + "/" + geneName + strainName + "_" + contig + "_" + start + "-" + end + "_" + gene + "_prot.fasta"
#    fastaHeader = ">" + geneName + strainName + "_" + contig + ":" + start + "-" + end + "_" + gene + "\n"
#    file = open(outputFile, 'w')
#    file.write(fastaHeader)
#    file.write(sequence + "\n")
#    file.close()

if multiSeq:
    outputFile = outputDirectory + "/" + strainName + "_prot.fasta"
    with open(outputFile, 'w') as ch:
        ch.write("")
    
    for gene, attributes in genes2.iteritems():
        #print gene, attributes
        contig = attributes[0][0].split(":")[0] + "_" + attributes[0][0].split(":")[1]
        start = attributes[0][1]
        end = attributes[0][2]
        sequence = attributes[0][3]
        fastaHeader = ">" + geneName + contig + ":" + start + "-" + end + "_" + gene + "\n"
        file = open(outputFile, 'a')
        file.write(fastaHeader)
        file.write(sequence + "\n")
        file.close()
else:
    for gene, attributes in genes2.iteritems():
        #print gene, attributes
        contig = attributes[0][0]
        start = attributes[0][1]
        end = attributes[0][2]
        sequence = attributes[0][3]
        outputFile = outputDirectory + "/" + strainName + "/" + geneName + strainName + "_" + contig + "_" + start + "-" + end + "_" + gene + "_prot.fasta"
        fastaHeader = ">" + geneName + strainName + "_" + contig + ":" + start + "-" + end + "_" + gene + "\n"
        file = open(outputFile, 'w')
        file.write(fastaHeader)
        file.write(sequence + "\n")
        file.close()




