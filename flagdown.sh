#!/usr/bin/bash

echo ""
echo "You are running Flagdown. This will (eventually) check a single assembled genome for clusters of enzymes in proximity to user specified Flag gene types."
echo "If you have any problems, consult the README. If that doesn't help, tweet me."
echo "Author: Chris Pyatt, NCYC, Quadram Institute."
echo ""

if (($# == 0)); then echo "No options provided. Use -h for help."; exit 1; fi

# Install python dependencies (only needs doing once - comment out once installed).
echo "Installing python dependencies if necessary."
pip install pyfaidx
pip install biopython
pip install Pillow

# Check hmmer installation.
hash nhmmer 2>/dev/null || { echo >&2 "Requirement: HMMER doesn't appear to be installed. Make sure it's added to your PATH.  Aborting."; exit 1; }


while getopts ":g:a:o:i:he" opt; do
    case ${opt} in
        h )
            echo ""
            echo "Usage:"
            echo "Run the following commands from the FindClustersPipeline directory (containing the scripts)."
            echo ""
            echo "    bash runClusterPipeline.sh -h                  Display this help message."
            echo ""
            echo "    bash runClusterPipeline.sh -e                  Run an example analysis to see what the results should look like."
            echo ""
            echo "    bash runClusterPipeline.sh -g [assemblies] -a [alignments] -o [output] -i [image]"
            echo ""
            echo "        -g    genome assemblies: the list of genome assemblies to be searched. These must be in FASTA format (with .fa or .fasta extensions) and the filename must be the strain/organism (e.g. NCYC123.fa)."
            echo "        -a    alignments: the list of stockholm alignments of target genes. These must be in STOCKHOLM format (with .stockholm, .sto, or .stk extensions) and the filename must be the gene name (e.g. EMT1.sto)."
            echo "        -o    output: your desired output directory, to which results will be written. Specify the full or relative path."
            echo "        -i    image: the name to give the final gene cluster image (no extension). Currently only PNG images are produced."
            echo ""
            exit 0
            ;;
        e )
            echo "Running an example analysis."
            assemblies="./example/NCYC1384.FA ./example/NCYC3267.fasta ./example/NCYC3431.FASTA"
            alignments="./example/EMT1.stockholm ./example/MAC1.stockholm ./example/MMF1.STO"
            output="./example"
            image="example"
            ;;
        g ) # assemblies
            assemblies=$OPTARG
            ;;
        a ) # alignments
            alignments=$OPTARG
            ;;
        o ) # output
            output=$OPTARG
            ;;
        i ) # image
            image=$OPTARG
            ;;
        : )
            echo "Option -$OPTARG requires an argument. Use -h for help." 1>&2
            exit 1
            ;;
        \? )
            echo "Invalid option. Use -h for help." 1>&2
            exit 1
            ;;
    esac
done

shift $((OPTIND -1))

# check whether provided variables actually contain files
for f in $assemblies; do [ -e "$f" ] && ass=0 || ass=1; break; done
for f in $alignments; do [ -e "$f" ] && aln=0 || aln=1; break; done

# exit or assign defaults if options not specified
if [ -z "$assemblies" ] || [ $ass == 1 ]; then echo "No genome assemblies provided. Use -g [assemblies]."; echo ""; exit 1; fi
if [ -z "$alignments" ] || [ $aln == 1 ]; then echo "No gene alignments provided. Use -a [alignments]."; echo ""; exit 1; fi
if [ -z "$output" ]; then echo "No output directory provided. Defaulting to current directory"; output="."; fi
if [ -z "$image" ]; then echo "No image name provided. Defaulting to 'image'"; image="image"; fi

echo ""
echo "Assemblies provided:           " $assemblies
echo "Alignments provided:           " $alignments
echo "Output directory provided:     " $output
echo "Image name chosen:             " $image

echo ""
echo "Running HMMER search."
python runHMMERsearch.py $output $assemblies $alignments

echo ""
echo "Getting top hit FASTA files for matching sequences."

for assembly in $assemblies; do
    strain=$(basename -- "$assembly")
    strain="${strain%.*}"
    echo "Processing ..... $strain"
    files="$output"/hmmer_outputs/"$strain"/*.out
    for file in $files; do
        echo "Parsing ..... $file"
        python mapCoordinates.py $file $assembly --outputDirectory=$output --numHits=1 --geneNameOn=True
    done
done

echo ""
echo "Finding clusters of matching sequences (same contig & within 10,000 nt)"

# get number of genes (length of cluster)
arr=($alignments)
len=${#arr[@]}

python findClusters.py "$output"/found_cluster_files "$output"/clusters.txt $len

echo ""
echo "Drawing PNG image of clusters."
python drawClusters.py "$output"/clusters.txt "$output"/"$image"
echo ""
