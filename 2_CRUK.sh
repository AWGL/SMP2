#!/bin/bash
set -euo pipefail

# Description: CRUK BaseSpace app pipeline
# Author: Sara Rey and Chris Medway
# Status: RELEASE
Version="1.1.4"

# location of basespace CLI v2 binary
BS=/home/transfer/bs

# How to use
# bash 2_CRUK.sh <path_to_sample_sheet> <sample_pairs_text_file (optional)>

# load sample variables
. *.variables

# loaded in from sample variables
NEGATIVE="$negative"

# variables
CONFIG="pmg-euc1"

# usage checking
if [ "$#" -lt 1 ]
    then
    echo "Commandline args incorrect. Usage: $0 <path_to_sample_sheet> <sample_pairs_text_file (optional)>." 
    exit -1
fi

# variables dependent on command line arguments
INPUTFOLDER="$1"
FASTQFOLDER="$INPUTFOLDER""/*/trimmed/"


# check if the sample sheet indicates that a manual pairs file should be created
if [ $pairs == 0 ]
then
    SAMPLEPAIRS="$INPUTFOLDER""SamplePairs.txt"
    makePairs=1
elif [ $pairs == 1 ] && [ "$#" -lt 2 ]
    then
    echo "SamplePairs file requires manual generation. Create in script directory and relaunch" \
    "2_CRUK.sh passing a properly formatted pairs file as the third command line argument. See README.MD for details."
    exit 1
elif [ $pairs == 1 ] && [ "$#" -eq 2 ]
    then
    SAMPLEPAIRS="$2"
    # Skip generation of a SamplePairs.txt file as it is already created and passed in as a command line argument
    makePairs=-1
fi


# notify the user that all samples in the sample sheet will be uploaded
echo "All samples on the SampleSheet.csv will be uploaded to BaseSpace."


# declare an array to store the sample ids in order
declare -a samplesArr


# parse SampleSheet
function parseSampleSheet {

    echo "Parsing sample sheet"
	
    # obtain project name from sample sheet
    projectName=$(grep "Experiment Name" "$INPUTFOLDER""SampleSheet.csv" | cut -d, -f2 | tr -d " ")

    echo "Project name is "$projectName

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

    # create/clear file which holds the sample name
    # pair the samples assuming the order tumour then normal and create a file of these pairs
    > "$SAMPLEPAIRS"
	
    
    # create array containing the samples that are not tumour-normal pairs (supports the option for more than one excluded sample)
    notPairs+=("$NEGATIVE")

	
    # exclude non tumour-normal pairs from pair file creation (e.g. negative control)	
    grep -f <(printf -- '%s\n' "${notPairs[@]}") -v <(printf '%s\n' "${samplesArr[@]}") | awk -F '\t' 'NR % 2 {printf "%s\t", $1;} !(NR % 2) {printf "%s\n", $1;}' >"$SAMPLEPAIRS"

}


function locateFastqs {

    # create array of all files for fq upload- use sample pairs file to support option for manual creation of sample pairs file
    mapfile -t pairsArr < "$SAMPLEPAIRS"
    fastqArr=("${notPairs[@]}" "${pairsArr[@]}")

    echo "Upload fastqs"
	
    for fastq in $(printf -- '%s\n' "${fastqArr[@]}")
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

# loop over app sessions and check to see if all app sessions are complete
numberOfAppSessions=$(wc -l  < appsessions.txt)
numberOfCompleteAppSessions=0

until [ $numberOfCompleteAppSessions -eq $numberOfAppSessions ]
do

    numberOfCompleteAppSessions=0

    while read session
    do
      	appId=$(echo $session | cut -d' ' -f1)
        status=$($BS list appsessions --config "pmg-euc1" --filter-field Id --filter-term $appId -f csv -F ExecutionStatus)

        if [ $status == "Complete" ]; then numberOfCompleteAppSessions=$[$numberOfCompleteAppSessions +1]; fi

    done < ./appsessions.txt
done

echo "$numberOfCompleteAppSessions app sessions are complete. Downloading data..."

# kick off next script
bash ./3_CRUK.sh >3_CRUK.out 2>3_CRUK.err
