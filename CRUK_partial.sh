#!/bin/bash
set -euo pipefail

#Description: CRUK Basespace app pipeline- draft version
#Author: Sara Rey
#Status: DEVELOPMENT/TESTING
Version=0.0


# How to use
# bash CRUK_draft.sh <path_to_sample_sheet> <path_to_results_location> <config_file_name> <Sample Pairs text file>
# /Users/sararey/Documents/cruk_test_data/rawFQs/ # for reference- path to sample sheet and fastqs
# Requires bash version 4 or above

# Usage checking
if [ "$#" -lt 2 ]
	then
		echo "Commandline args incorrect. Usage: $0 <path_to_sample_sheet> <path_to_results_location>." 
		exit 0
fi

if [ "$#" -lt 3 ]
	then
		SAMPLEPAIRS="SamplePairs.txt"
	else
		SAMPLEPAIRS="$3"
fi

# Variables
APPNAME="SMP2 v2"
INPUTFOLDER="$1"
RESULTSFOLDER="$2"


# Declare an array to store the patient name and sample id
declare -A samplePatient


# Parse SampleSheet
function parseSampleSheet {

	echo "Parsing sample sheet"
	
	# Obtain project name from sample sheet
	#projectName=$(grep "Experiment Name" "$INPUTFOLDER"SampleSheet.csv | cut -d, -f2 | tr -d " ")
	projectName="sr2" #temp var	

	# Obtain list of samples from sample sheet
	for line in $(sed "1,/Sample_ID/d" "$INPUTFOLDER"SampleSheet.csv | tr -d " ")
	do 
		
		# Obtain sample name and patient name		
		samplename=$(printf "$line" | cut -d, -f1 | sed 's/[^a-zA-Z0-9]+/-/g')
		patientname=$(printf "$line" | cut -d, -f2 | sed 's/[^a-zA-Z0-9]+/-/g')

	 	# Skip any empty sample ids- both empty and whitespace characters (but not tabs at present)
	 	if [[ "${#samplename}" = 0 ]] || [[ "$samplename" =~ [" "] ]]
		then
			continue
	 	fi

	 	# Append information to array
		samplePatient["$samplename"]="$patientname"
	done

}

function pairSamples {

	echo "Pairing samples"

	# Create/clear file which holds the sample name and the patient identifiers
	>"$SAMPLEPAIRS"

	count=0

	awk -v var="${!samplePatient[*]}" 'NR % 2 {print var[1]}'	
	
	for sample in ${!samplePatient[@]}
	do

		echo $sample	
		#if (( $count % 2 == 0 ))
	 	#then
	 	 	#tumour="$sample"
	 	#else
	 	 	#normal="$sample"
	 	 	# Add paired samples to a file- tumour sample first, normal sample second
	 	 	#printf "$tumour">>"$SAMPLEPAIRS"
	 	 	#printf "%s\t">>"$SAMPLEPAIRS"
	 	 	#printf "$normal">>"$SAMPLEPAIRS"
	 	 	#printf "%s\n">>"$SAMPLEPAIRS"
		#fi
		#((count++))
	done

}


function locateFastqs {

	echo "Uploading fastqs"

	for fastq in $( printf -- '%s\n' "${samplePatient[@]}" | grep -f "not_bs_samples.txt" -v )
	do
		f1=$INPUTFOLDER${fastq}*_R1_*.fastq.gz
		f2=$INPUTFOLDER${fastq}*_R2_*.fastq.gz

	done
}



# Call the functions

# Parse sample sheet to obtain required information
parseSampleSheet $INPUTFOLDER


# Pair samples according to order in sample sheet- make a command line argument optional to manually create
# for NEQAS samples etc.
pairSamples


# Get fastqs
locateFastqs $INPUTFOLDER


# Launch app for each pair of samples in turn as tumour normal pairs then download analysis files
echo "Launching app"
while read pair
do
	printf $pair	
	tum=$(printf "$pair" | cut -d" " -f1)
	nor=$(printf "$pair" | cut -d" " -f2)

	echo $tum
	echo $nor

done <"$SAMPLEPAIRS"


