
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

# Save all terminal outputs here
out_filename=data/$EXPERIMENT_NAME/reads/processed/terminal_capture.txt

cli=".scripts/reads1.sh -e "
cli+=$@

printf "\nSTART NEW COMMAND\n%s\n" "$cli" >> $out_filename

echo "EXPERIMENT : $EXPERIMENT_NAME" |& tee -a $out_filename;
echo "INDEX : $INDEX" |& tee -a $out_filename;
echo "ADAPTERS :" |& tee -a $out_filename;
for i in "${adapters[@]}";
do
    echo "    ${i}" |& tee -a $out_filename;
done
echo "READS :" |& tee -a $out_filename;
for i in "${reads[@]}";
do
    echo "    ${i}" |& tee -a $out_filename;
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

    echo " * Begin Aligning ${i} with $INDEX" |& tee -a $out_filename;

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
        echo "$trimmed_str exists on your filesystem." |& tee -a $out_filename;
    else
        # Trim, crop, and output file
        echo " * Begin Trimmomatic for ${i}..." |& tee -a $out_filename;
        # LEADING: Cut bases off the start of a read, if below a threshold quality
        # TRAILING: Cut bases off the end of a read, if below a threshold quality
        # MINLEN: Drop the read if it is below a specified length
        # CROP: Cut the read to a specified length by removing bases from the end
        # LEADING:20 TRAILING:20 MINLEN:48 CROP:20

        # Important: MINLEN:48 REMOVES READS WITH 2 BAD BP TRAILING
        java -jar tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 $input_str $trimmed_str $adapter_str LEADING:20 TRAILING:20 MINLEN:48 CROP:20 |& tee -a $out_filename;
    fi

    if [[ -f $mapped_str ]]
    then
        echo "$mapped_str exists on your filesystem." |& tee -a $out_filename;
    else
        # Run the bowtie alignment
        echo " * Begin Bowtie1 for ${i} and $INDEX..." |& tee -a $out_filename;
        # -t = print time in terminal
        # -v 3 = v-alignment(no quality checks), 3 reportable alignments allowed
        # -a = reportable all alignments
        # -m 1 = only report 1 unique alignment
        # --best --strata = tbh, not sure behavior, output in best-to-worst order
        bowtie -t -v 3 -a -m 1 --best --strata $index_str $trimmed_str $mapped_str |& tee -a $out_filename;
    fi

    # Create TA map for the read to the index
    # This will also try to map to a combined TAmap if one exists
    echo " * Creating TAmap for $INDEX and ${i}..." |& tee -a $out_filename;
    python3 scripts/readTAmap.py --experiment=$EXPERIMENT_NAME --index=$INDEX --map=${i} |& tee -a $out_filename;

done

echo "Finished processing the reads." |& tee -a $out_filename;
# Generate the next possible command

cmd_str="python3 scripts/analysis.py --experiment $EXPERIMENT_NAME --index $INDEX";

reads_dir=data/$EXPERIMENT_NAME/reads/

shopt -s extglob nullglob globstar
reads=($reads_dir*.fastq)
reads+=($reads_dir*.fq)

cmd_str+=" --controls read_name";
cmd_str+=" --samples read_name";

printf "\nExample command for analysis:\n" |& tee -a $out_filename;
echo "$cmd_str" |& tee -a $out_filename;

printf "\nThese are all the reads in the experiment...\n" |& tee -a $out_filename;
for item in "${reads[@]}"; do
    filename=$(basename -- "${item}");
    printf '%s\n' "${filename%.*}" |& tee -a $out_filename;
done
