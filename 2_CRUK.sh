#!/bin/bash
set -euo pipefail

# Description: CRUK BaseSpace app pipeline
# Author: Sara Rey and Chris Medway
# Status: RELEASE
Version="1.1.3"

# location of basespace CLI v2 binary
BS=/home/transfer/bs

# How to use
# bash 2_CRUK.sh <path_to_sample_sheet> <name_of_negative_control_sample> <sample_pairs_text_file (optional)>

# variables
CONFIG="pmg-euc1"

# load variables
. variables

# usage checking
if [ "$#" -lt 2 ]
    then
    echo "Commandline args incorrect. Usage: $0 <path_to_sample_sheet> <name_of_negative_control_sample> <sample_pairs_text_file (optional)>." 
    exit -1
fi

# variables dependent on command line arguments
INPUTFOLDER="$1"
NEGATIVE="$2"
NOTBASESPACE="$INPUTFOLDER""not_bs.txt"
FASTQFOLDER="$INPUTFOLDER""/*/trimmed/"


# check if the sample sheet indicates that a manual pairs file should be created
if [ $pairs == 0 ]
then
    SAMPLEPAIRS="$INPUTFOLDER""SamplePairs.txt"
    makePairs=1
elif [ $pairs == 1 ] && [ "$#" -lt 3 ]
    then
    echo "SamplePairs file requires manual generation. Create in script directory and relaunch" \
    "2_CRUK.sh passing pairs file as the third command line argument."
    exit 1
elif [ $pairs == 1 ] && [ "$#" -eq 3 ]
    then
    SAMPLEPAIRS="$3"
    # Skip generation of a SamplePairs.txt file
    makePairs=-1
fi


# check for the presence of the file with samples not to upload to BaseSpace in the same directory as the script
if [[ -e $NOTBASESPACE ]]
then
    samples_to_skip=1
    # check that the provided file is not empty
    if ! [[ -s $NOTBASESPACE ]]
    then
        echo "The file "$NOTBASESPACE" is empty. When this file exists, it must contain the names of samples that are in the SampleSheet.csv, but should not be uploaded to BaseSpace."
        exit -1
    fi
else
    samples_to_skip=-1
    # notify the user that all samples in the sample sheet will be uploaded
    echo "No "$NOTBASESPACE" file found in the same directory as the script. All samples on the SampleSheet.csv will be uploaded to BaseSpace."
fi


# declare an array to store the sample ids in order
declare -a samplesArr
# initial entry created to avoid downstream error when appending to array
samplesArr+=1 


# parse SampleSheet
function parseSampleSheet {

    echo "Parsing sample sheet"
	
    # obtain project name from sample sheet
    projectName=$(grep "Experiment Name" "$INPUTFOLDER""SampleSheet.csv" | cut -d, -f2 | tr -d " ")

    echo $projectName

    # obtain list of samples from sample sheet
    for line in $(sed "1,/Sample_ID/d" "$INPUTFOLDER""SampleSheet.csv" | tr -d " ")
    do
        # obtain sample name and patient name		
        samplename=$(printf "$line" | cut -d, -f1 | sed 's/[^a-zA-Z0-9]+/-/g')

        # skip any empty sample ids- both empty and whitespace characters (but not tabs at present)
        if [[ "${#samplename}" = 0 ]] || [[ "$samplename" =~ [" "] ]]
        then
	    continue
        fi

        # append information to list array- to retain order for sample pairing
        samplesArr=("${samplesArr[@]}" "$samplename")
    done
}


function pairSamples {

    echo "Pairing samples"

    # create/clear file which holds the sample name and the patient identifiers
    > "$SAMPLEPAIRS"
	
    # iterate through the samples and exclude any samples that are not for basespace
    # pair the samples assuming the order tumour then normal and create a file of these pairs
    # create array containing the samples that are not tumour-normal pairs
    # check if there are any samples on the run that are not for BaseSpace and so should not be paired
    if [[ -e $NOTBASESPACE ]]
    then
        mapfile -t notPairs < $NOTBASESPACE
        notPairs=("${notPairs[@]}" "$NEGATIVE") 
    else
        notPairs+=("$NEGATIVE")
    fi	
	
    # exclude non tumour-normal pairs from pair file creation		
    grep -f <(printf -- '%s\n' "${notPairs[@]}") -v <(printf '%s\n' "${samplesArr[@]:1}") | awk -F '\t' 'NR % 2 {printf "%s\t", $1;} !(NR % 2) {printf "%s\n", $1;}' >"$SAMPLEPAIRS"
}


function locateFastqs {

    echo "Uploading fastqs"

    if [[ "$samples_to_skip" == 1 ]]
    then
        fastqlist=$( printf -- '%s\n' "${samplesArr[@]:1}" | grep -f "$NOTBASESPACE" -v )
    else
        fastqlist=$(printf -- '%s\n' "${samplesArr[@]:1}")
    fi
	
    for fastq in $(printf -- '%s\n' "$fastqlist")
    do
        f1=$FASTQFOLDER${fastq}*_R1_*.fastq.gz
        f2=$FASTQFOLDER${fastq}*_R2_*.fastq.gz

        # added in version 1.1.3. bscli v2 requires sample unicity
        # therefore sample names are prefixed with project name
        cp $f1 ./"$projectName"-`basename $f1`
        cp $f2 ./"$projectName"-`basename $f2`

        f1=./"$projectName"-`basename $f1`
        f2=./"$projectName"-`basename $f2`
                        
        # upload fastq to biosample
        $BS upload dataset --config "$CONFIG" --project $projectId $f1 $f2
    done
}


function launchApp {

    # launch app for each pair of samples in turn as tumour normal pairs then download analysis files
	
    # obtain basespace ID of negative control- this is not an optional input through the commandline app launch
    negId=$($BS list biosample --config "$CONFIG" --filter-field BioSampleName --filter-term "$projectName"-"$NEGATIVE" --terse)
	
    while read pair
    do
        # stop iteration on first empty line of SamplePairs.txt file in case EOF marker is absent for any reason
        if [[ -z $pair ]]
        then
            return 0
        fi

        echo "Launching app for ""$pair"
			
        tum=$(printf "$pair" | cut -d$'\t' -f1)
        nor=$(printf "$pair" | cut -d$'\t' -f2)

        # obtain sample ids from basespace
        tumId=$($BS list biosample --config "$CONFIG" --filter-field BioSampleName --filter-term "$projectName"-"$tum" --terse)
        norId=$($BS list biosample --config "$CONFIG" --filter-field BioSampleName --filter-term "$projectName"-"$nor" --terse)

        # launch app and store the appsession ID	
       appSessionId=$($BS launch application \
               --config "$CONFIG" \
               --name "SMP2 v2" \
               --app-version "1.1.2" \
               --option tumour-sample-id:$tumId \
               --option normal-sample-id:$norId \
               --option negative-sample-id:$negId \
               --option project-id:$projectId \
               --option basespace-labs:1 \
               --terse )

        # save file that will track the appsession IDs for each sample pair
        echo -e $appSessionId $tum $nor $projectName >> ./appsessions.txt

    done <"$SAMPLEPAIRS"
}


# call the functions
# check sample sheet exists at location provided
if ! [[ -e "$INPUTFOLDER""SampleSheet.csv" ]]
then
    echo "Sample Sheet not found at input folder location"
    exit -1
fi

# parse sample sheet to obtain required information
parseSampleSheet

# pair samples according to order in sample sheet if manually created pairs file has not been supplied
if [[ "$makePairs" == 1 ]]
then
    pairSamples
fi

# count number of paired samples
numPairs=$(cat "$SAMPLEPAIRS" | cut -f2 | sed '/^\s*$/d' | wc -l)

# read out the sample pairs in the order tumour blood with each pair on a new line 
echo "Displaying sample pairs:" 
cat "$SAMPLEPAIRS"
printf $'\n'
echo "Abort the script if the samples are paired incorrectly and create a file of the pairs (see README.MD for details about this file)." 
printf $'\n'

# create project in basespace
echo "Creating project"
$BS create project --name "$projectName" --config "$CONFIG"

# get project ID
projectId=$($BS get project --name $projectName --config $CONFIG --terse)

# Get fastqs and upload to basespace
locateFastqs

# Kick off the app for each pair in turn
if [ -e "appsessions.txt" ]; then rm appsessions.txt; fi
launchApp

# queue next script in the pipeline for half an hours time
at now +50 minutes -f ./3_CRUK.sh >3_CRUK.out 2>3_CRUK.err
