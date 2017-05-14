#!/bin/bash

####################### define executables
FASTAtoRF="./bin/FASTAtoRF"
rndForest="./src/_R/_predict.RF.R"
RC="./bin/RC"

####################### identify the input arguments
jobid=$1
proteins=$2

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash RCOpt.sh <jobID> <C2H2_ZFP.fasta>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input FASTA file for the target protein(s): "$proteins

if [ -e "$proteins" ]; then
	echo "Protein sequence file found."
else
	echo "ERROR: Protein sequence file was not found."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder
RF_in=$tmp_folder"/_predict.in"
RF_out=$tmp_folder"/_predict.RF.out"

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
$RC -rf $out_file.RF_out.txt -out $out_file  >>$out_folder/log.step2.txt









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
