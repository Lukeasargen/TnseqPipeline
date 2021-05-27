
# Get the flags
while getopts :e:f:g:o: flag
do
    case "${flag}" in
        e) EXPERIMENT_NAME=${OPTARG};;
        f) FASTA=${OPTARG};;
        g) GENBANK=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

echo "EXPERIMENT_NAME: $EXPERIMENT_NAME";
echo "FASTA: $FASTA";
echo "GENBANK: $GENBANK";
echo "OUTPUT: $OUTPUT";


index_dir=data/$EXPERIMENT_NAME/indexes/$OUTPUT

shopt -s extglob nullglob globstar
indexes=($index_dir*.ebwt)
if [[ -f $indexes ]]
then
    echo "Indexes already exist on your filesystem."
else
    # Create the bowtie index
    bowtie-build data/$EXPERIMENT_NAME/references/$FASTA data/$EXPERIMENT_NAME/indexes/$OUTPUT
fi


# Run the python script
python3 scripts/referenceTAlist.py --experiment=$EXPERIMENT_NAME --fasta=$FASTA --genbank=$GENBANK --output=$OUTPUT


# Generate the next possible command
printf "\nCommand for processing reads:\n"

cmd_str="./scripts/reads1.sh -e $EXPERIMENT_NAME -i $OUTPUT";

adapter_dir=data/$EXPERIMENT_NAME/adapters/
reads_dir=data/$EXPERIMENT_NAME/reads/

adapters=($adapter_dir*.fasta)
adapters+=($adapter_dir*.fa)
reads=($reads_dir*.fastq)
reads+=($reads_dir*.fq)

delim=""
adapter_str=""
for item in "${adapters[@]}"; do
    filename=$(basename -- "${item}");
    adapter_str="$adapter_str$delim$filename"
    delim=","
done
cmd_str+=" -a $adapter_str"


delim=""
reads_str=""
for item in "${reads[@]}"; do
    filename=$(basename -- "${item}");
    reads_str="$reads_str$delim$filename"
    delim=","
done
cmd_str+=" -r $reads_str"

echo "$cmd_str";
