
IFS=,  # comma delimiter for args

while getopts :e:i:a:r: flag
do
    case $flag in
        e ) EXPERIMENT_NAME=${OPTARG};;
        i ) INDEX=${OPTARG};;
        a ) adapters=($OPTARG) ;;
        r ) reads=($OPTARG) ;;
    esac
done

echo "EXPERIMENT : $EXPERIMENT_NAME";
echo "INDEX : $INDEX";
echo "ADAPTERS :";
for i in "${adapters[@]}";
do
    echo "    ${i}";
done
echo "READS :";
for i in "${reads[@]}";
do
    echo "    ${i}";
done

# Make the adpater list into a string for trimmomatic
adapter_str="ILLUMINACLIP"  # Always is illuminaclip adpater
for i in "${adapters[@]}";
do
    adapter_str+=":data/$EXPERIMENT_NAME/adapters/${i}";
done
# From the manual
# seedMismatches:palindromeClipThreshold:simpleClipThreshold
# seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
# palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
# simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
adapter_str+=":2:30:10";
# echo $adapter_str;


# Process each read
for i in "${reads[@]}";
do

    echo " * Begin Aligning ${i} with $INDEX"

    # Make the file paths as strings
    input_str=data/$EXPERIMENT_NAME/reads/${i};
    filename=$(basename -- "${i}");
    extension="${filename##*.}";
    filename="${filename%.*}";
    trimmed_str=data/$EXPERIMENT_NAME/reads/processed/${filename}_trimmed.$extension;
    mapped_str=data/$EXPERIMENT_NAME/reads/processed/${filename}_trimmed_mapped;
    index_str=data/$EXPERIMENT_NAME/indexes/$INDEX;
    # echo "input_str : $input_str";
    # echo "trimmed_str : $trimmed_str";
    # echo "mapped_str : $mapped_str";
    # echo "index_str : $index_str";

    # If the files exist, then skip the processing
    if [[ -f $trimmed_str ]]
    then
        echo "$trimmed_str exists on your filesystem."
    else
        # Trim, crop, and output file
        echo " * Begin Trimmomatic for ${i}...";
        # LEADING: Cut bases off the start of a read, if below a threshold quality
        # TRAILING: Cut bases off the end of a read, if below a threshold quality
        # MINLEN: Drop the read if it is below a specified length
        # CROP: Cut the read to a specified length by removing bases from the end
        # LEADING:20 TRAILING:20 MINLEN:48 CROP:20

        # Important: MINLEN:48 REMOVES READS WITH 2 BAD BP TRAILING
        java -jar tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $input_str $trimmed_str $adapter_str TRAILING:20 MINLEN:48 CROP:20;
    fi

    if [[ -f $mapped_str ]]
    then
        echo "$mapped_str exists on your filesystem."
    else
        # Run the bowtie alignment
        echo " * Begin Bowtie1 for ${i} and $INDEX...";
        # -m 1
        bowtie -t -v 3 -a --best --strata $index_str $trimmed_str $mapped_str;
    fi

    # Create TA map for the read to the index
    # This will also try to map to a combined TAmap if one exists
    echo " * Creating TAmap for ${i} and $INDEX...";
    python3 scripts/readTAmap.py --experiment=$EXPERIMENT_NAME --index=$INDEX --map=${i};

done

echo "Finished processing the reads.";
# Generate the next possible command

cmd_str="python3 scripts/analysis.py --experiment $EXPERIMENT_NAME --index $INDEX";

reads_dir=data/$EXPERIMENT_NAME/reads/

shopt -s extglob nullglob globstar
reads=($reads_dir*.fastq)
reads+=($reads_dir*.fq)

cmd_str+=" --controls read_name";
cmd_str+=" --samples read_name";

printf "\nExample command for analysis:\n";
echo "$cmd_str";

printf "\nThese are all the reads in the experiment...\n"
for item in "${reads[@]}"; do
    filename=$(basename -- "${item}");
    printf '%s\n' "${filename%.*}"
done
