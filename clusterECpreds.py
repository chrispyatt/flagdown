#!/bin/python

import sys
from sys import argv
import argparse
import re
from collections import defaultdict
from Bio import ExPASy
from Bio.ExPASy import Enzyme
from decimal import *
getcontext().prec = 3


# These lines assign the arguments given at command line to variables that can be used within the script. You can type 'python proteinFromGFF.py -h' at the command line to see this.
parser = argparse.ArgumentParser(description='Parse a file containing EC number predictions (from DETECT) and determine whether there are any clustered enzymes in the genome of interest')
parser.add_argument('inputFile',
                    help='The DETECT output file (*.top).')
parser.add_argument('enzymeDB',
                    help='The path to the ExPASy Enzyme database.')
parser.add_argument('--outputDir', default=".",
                    help='The directory in which to save the output (cluster info). Default is current working directory.')
parser.add_argument('--sepGenes', default="6",
                    help='The distance allowed between clustered enzymes (measured in number of genes). Default is six. Lowering this parameter will more strictly define each cluster but will also risk splitting up more widely spaced clusters.')
parser.add_argument('--sepDist', default="10000",
                    help='The distance allowed between clustered enzymes (measured in number of bases). Default is 10000 (a very wide net to cast). Lowering this parameter will have the same effect as lowering --sepGenes')
parser.add_argument('--minClusterSize', default="2",
                    help='The minimum number of genes together to be regarded as a cluster. Default is two.')
parser.add_argument('--strainName', default=None,
                    help='The name of the strain in the input file. This will be used to name the output file. The default behaviour is to take the input filename minus the ".top" part.')


args = parser.parse_args()
inputFile = args.inputFile
enzymeDB = args.enzymeDB
outputDir = args.outputDir
sepGenes = args.sepGenes
sepDist = args.sepDist
minClusterSize = args.minClusterSize
strainName = args.strainName

enzymeDB_dict = {}
#db = open(enzymeDB)
with open(enzymeDB) as db:
    for record in Enzyme.parse(db):
        id_ec = record["ID"]
        de = record["DE"]
        enzymeDB_dict[id_ec] = de
#db.close()

fileName = inputFile.split("/")[-1]
if not (fileName.split(".")[-1] == "top"):
    sys.exit('ERROR! Wrong filetype! Input should be a ".top" file!')
if not strainName:
    strainName = fileName.split(".")[0]

# open the input file (for reading by default)
#fh = open(inputFile)
# initialise dictionary to hold enzyme data for each contig
group_enzymes = defaultdict(list)

with open(inputFile) as fh:
    for line in fh:
        if line.startswith("ID"):
            #print "ID line is " + line
            continue
        else:
            #print "Not ID line: " + line
            id = line.split("\t")[0]
            ec = line.split("\t")[1]
            ils = line.split("\t")[2]
            num_pos_hits = line.split("\t")[3]
            num_neg_hits = re.split("[\t\n]", line)[4]
            
            strain = re.split("[-_:]", id)[0]
            contig = re.split("[-_:]", id)[1]
            start = re.split("[-_:]", id)[2]
            end = re.split("[-_:]", id)[3]
            genename = re.split("[-_:]", id)[4]
            
            if len(ec.split(";")) > 1:
                multi_ec = ec.split(";")
                enzyme_name = enzymeDB_dict[multi_ec[0]]
                for ec_i in multi_ec[1:]:
                    alt_name = enzymeDB_dict[ec_i]
                    enzyme_name = enzyme_name + ";" + alt_name
            else:
                enzyme_name = enzymeDB_dict[ec]
            
            attributes = [strain, contig, start, end, genename, ec, ils, num_pos_hits, num_neg_hits, enzyme_name]
            
            group_enzymes[contig].append(attributes)


#fh.close()
#print group_enzymes

# initialise new dictionary to hold lists of enzyme clusters for each contig
clusters_dict = defaultdict(list)
num_clusters = 0
min_len_cluster = None
max_len_cluster = None
total_len_clusters = 0

for contig, enzymes in group_enzymes.iteritems():
    if len(enzymes) >= 2:
        enzymes.sort(key=lambda x: x[4])
        previous = None
        prev_end = 0
        cluster = []
        for enzyme in enzymes:
            numeric = int(enzyme[4].split("g")[1])
            curr_start = int(enzyme[2])
            if not previous or (numeric - previous <= int(sepGenes) and curr_start - prev_end <= int(sepDist)):
                cluster.append(enzyme)
            else:
                if len(cluster) >= int(minClusterSize):
                    clusters_dict[contig].append(cluster)
                cluster = [enzyme]
            previous = int(enzyme[4].split("g")[1])
            prev_end = int(enzyme[3])
        if len(cluster) >= int(minClusterSize):
            clusters_dict[contig].append(cluster)


#print clusters_dict

for contig, clusters in clusters_dict.iteritems():
    num_clusters = num_clusters + len(clusters)
    for cluster in clusters:
        total_len_clusters = total_len_clusters + len(cluster)
        if not min_len_cluster or len(cluster) < min_len_cluster:
            min_len_cluster = len(cluster)
        if not max_len_cluster or len(cluster) > max_len_cluster:
            max_len_cluster = len(cluster)

av_len_cluster = Decimal(total_len_clusters) / Decimal(num_clusters)
outFileName = outputDir + "/" + strainName + ".EC_clusters"
#outFile = open(outFileName, 'w')
with open(outFileName, 'w') as outFile:
    
    header = ("# This output was generated by clusterECpreds.py \n"
        + "# This program takes the output from DETECT (http://www.compsysbio.org/projects/DETECT) and finds clusters\n# of enzymes within contigs (cluster = within x genes of each other where x is specified by the user). \n"
        + "# Author: Chris Pyatt \n"
        + "# Year: 2018 \n"
        + "#\n-- Input parameters --\n"
        + "# The input file was: " + inputFile + "\n"
        + "# The output directory was: " + outputDir + "\n"
        + "# The specified distance between clustered enzymes was: " + sepGenes + " genes and " + sepDist + " bases.\n"
        + "# The specified minimum cluster size was: " + minClusterSize + " genes.\n"
        + "#\n-- Summary Statistics --\n" + "# There were " + str(num_clusters) + " enzyme clusters found in " + str(len(clusters_dict.keys())) + " contigs.\n"
        + "# Average cluster length was " + str(av_len_cluster) + " genes (min=" + str(min_len_cluster) + ", max=" + str(max_len_cluster) + ").\n")
    outFile.write(header)

    for contig, clusters in clusters_dict.iteritems():
        contig_record = "#\n#\n###\n-------------------- Predicted enzyme clusters on contig number: " + contig + " --------------------\n###\n"
        for cluster in clusters:
            cluster_record = "# Start cluster \n# Gene_Number Start End EC_Number Recommended_Name ILS NumPosHits NumNegHits \n"
            for gene in cluster:
                cluster_record = cluster_record + gene[4] + "\t" + gene[2] + "\t" + gene[3] + "\t" + gene[5] + "\t" + gene[9] + "\t" + gene[6] + "\t" + gene[7] + "\t" + gene[8] + "\n"
            
            cluster_record = cluster_record + "# End cluster \n###\n"
            contig_record = contig_record + cluster_record
        outFile.write(contig_record)

#outFile.close()










