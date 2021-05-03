# Deafult
MERGE=false

# Get the flags
while getopts e:f:g:o:m flag
do
    case "${flag}" in
        e) EXPERIMENT_NAME=${OPTARG};;
        f) FASTA=${OPTARG};;
        g) GENBANK=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        m) MERGE=true
    esac
done

echo "EXPERIMENT_NAME: $EXPERIMENT_NAME";
echo "FASTA: $FASTA";
echo "GENBANK: $GENBANK";
echo "OUTPUT: $OUTPUT";
echo "MERGE: $MERGE";

# Create the index
bowtie-build data/$EXPERIMENT_NAME/references/$FASTA.fasta data/$EXPERIMENT_NAME/indexes/$OUTPUT

# Run the python script
python scripts/createTAlist.py --experiment=$EXPERIMENT_NAME --fasta=$FASTA --genbank=$GENBANK --output=$OUTPUT

# Optional merge


$