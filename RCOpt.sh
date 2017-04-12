#!/bin/bash

####################### define executables
FASTAtoRF="./bin/FASTAtoRF"
rndForest="./src/_R/_predict.RF.R"
RCOpt="./bin/RCOpt"
memebin="<MEME Suite bin folder>"

####################### identify the input arguments
jobid=$1
proteins=$2
peaks=$3

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash RCOpt.sh <jobID> <C2H2_ZFP.fasta> <ChIP_seq.fasta>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input FASTA file for the target protein(s): "$proteins
echo "Input FASTA file for the peaks: "$peaks

if [ -e "$proteins" ]; then
	echo "Protein sequence file found."
else
	echo "ERROR: Protein sequence file was not found."
	exit
fi

if [ -e "$peaks" ]; then
	echo "Peak sequence file found."
else
	echo "ERROR: Peak sequence file not found."
	exit
fi

if [ -d "$memebin" ]; then
	echo "MEME bin directory found."
else
	echo "ERROR: MEME bin directory was not found (RCOpt.sh line 7)."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder
RF_in=$tmp_folder"/_predict.in"
RF_out=$tmp_folder"/_predict.RF.out"

centered=$tmp_folder"/centered.fa"
shuffled=$tmp_folder"/shuffled.fa"
all=$tmp_folder"/all.fa"

####################### prepare the peak sequences
# get the central region of the sequences, and also dinucleotide-shuffled sequences	
$memebin/fasta-center -len 100 < $peaks 1> $centered
$memebin/fasta-dinucleotide-shuffle -f $centered -t -dinuc 1> $shuffled
cat $centered $shuffled > $all

####################### define the output path
out_folder="./out/"$jobid
mkdir -p $out_folder
rm -f $out_folder/log.step1.txt
rm -f $out_folder/log.step2.txt
out_file=$out_folder"/results"

####################### convert the input FASTA file to a covariate matrix file for the RF script
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
	$FASTAtoRF -minl 2 -maxl 8 -span $i -fasta $proteins -out $RF_in.span$i >>$out_folder/log.step1.txt
	if [ $i == 2 ]; then
		cat $RF_in.span$i > $RF_in
	else
		cat $RF_in.span$i | sed 1d >> $RF_in
	fi
done

####################### run the RF script, and reformat it for the next step
Rscript $rndForest $jobid
sed 's/"//g' $RF_out > $out_file.RF_out.txt

####################### run the RCOpt script
$RCOpt -rf $out_file.RF_out.txt -fasta $all -out $out_file -mode 3  >>$out_folder/log.step2.txt









#*****************************************************************************************
# The following lines check the input/output, and produce appropriate messages
# If no error was detected in either input or output, the info messages will be written in
# ./out/<jobID>/log.info.txt
# Otherwise, the error messages will be written in
# ./out/<jobID>/log.error.txt
#*****************************************************************************************



####################### identify the input arguments
err=""
info=""

####################### define log files
step1=$out_folder/log.step1.txt
step2=$out_folder/log.step2.txt
report=$out_folder/results.report.txt

####################### check if any of the B1H-RC motifs were optimized

optimized=`cat $report | cut -f3 | grep '1' | wc -l`

if [ "$optimized" -le 0 ]; then
	err="ERROR: None of the predicted motifs from the provided C2H2-ZF proteins were enriched in the ChIP-seq peaks.\n"
else
	info="The predicted motifs of $optimized possible C2H2-ZF arrays were enriched in the ChIP-seq peaks.\n"$info
fi


####################### check if the peak sequence file had any valid sequences

numPeaks=`cat $step2 | grep 'sequences were read,' | head -n 1 | cut -d ' ' -f1`

if [ "$numPeaks" = "" ]; then
	err="ERROR: No sequences were found in the input FASTA for ChIP-seq peaks. Please check the input format.\n"
elif [ "$numPeaks" = "ERROR:" ]; then
	err="ERROR: No sequences were found in the input FASTA for ChIP-seq peaks. Please check the input format.\n"
elif [ "$numPeaks" -le 0 ]; then
	err="ERROR: No sequences were found in the input FASTA for ChIP-seq peaks. Please check the input format.\n"
else
	numPeaks=$((numPeaks/2))
	info="$numPeaks sequences were found in the input FASTA for ChIP-seq peaks.\n"$info
fi


####################### check if the C2H2-ZF sequences have had any ZF arrays

numArrays=`cat $step2 | grep 'motifs were read.' | head -n 1 | cut -d ' ' -f1`

if [ "$numArrays" = "" ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
elif [ "$numArrays" = "ERROR:" ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
elif [ "$numArrays" -le 0 ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
else
	info="$numArrays possible C2H2-ZF arrays were tested.\n"$info
fi


####################### check if the C2H2-ZF file had any valid sequences

numC2H2=`cat $step1 | grep 'sequences were read.' | head -n 1 | cut -d ' ' -f1`

if [ "$numC2H2" = "" ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
elif [ "$numC2H2" = "ERROR:" ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
elif [ "$numC2H2" -le 0 ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
else
	info="$numC2H2 sequences were found in the input FASTA for C2H2-ZF proteins.\n"$info
fi


####################### write the appropriate messages to the output

if [ "$err" = "" ]; then
	echo -e -n $info > $out_folder/log.info.txt
else
	echo -e -n $err > $out_folder/log.error.txt
fi
