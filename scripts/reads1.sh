

# Get the positinoal arguments

EXP=$1
INDEX=$2
READ=$3

echo "EXPERIMENT : $EXP";
echo "INDEX : $INDEX";
echo "READ : $READ";

# Adapter argument is not required
if [ "$4" != "" ]; then
    ADAPTER=$4
    echo "Using ILLUMINACLIP and supplied adapter..."
    echo "ADAPTER : $ADAPTER";
else
    echo "Using only ILLUMINACLIP adapter"
fi



#java -jar Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 Data/$CONTROL.fastq Data/temp.fastq ILLUMINACLIP:Trimmomatic-0.36/PolyC_Adapter.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:48;
#java -jar Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 Data/temp.fastq Data/$CONTROL-trimmed.fastq CROP:20;
# bowtie-1.0.0/bowtie -v 3 -a --best --strata -m 1 -q bowtie-1.0.0/indexes/$BowtieName Data/$CONTROL-trimmed.fastq Data/$CONTROL-trimmed-bowtieMap;
# TODO - create TA map using the genome index

