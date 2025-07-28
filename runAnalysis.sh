#!/bin/bash

# Set the username for which to check the jobs
runningTimeLimitHours=12                          #hours
sleepTime=600                                      #seconds
username=$(whoami)
condorQFile='condor.out'
runningTime=0
runningTimeLimit=$((runningTimeLimitHours * 3600)) #seconds

# Set the output directory and file name
filename="productionList.txt"
folder=$( tail -n 1 $filename ) # reading the last line
outDate=`date +%d-%m`
details="_newEtaCorr"
outputDir="output/production/${folder}"




#Function to clean the output directory
clean(){
    echo "Deleting production..."
    rm -r output/production
    
    echo "Deleting csh..."
    rm -r output/csh
    
    echo "Deleting Local Libraries..."
    rm -r LocalLibraries.package
    rm LocalLibraries.zip
    
    echo "Deleting report..."
    rm -r output/report
    
    echo "Deleting list..."
    rm -r output/list
    
    echo "Deleting jobs..."
    rm -r output/jobs
    
    echo "Deleting scheduler files..."
    rm sched*
    rm *.session.xml
    
    echo "Done!"
}


# Function to check if all jobs for the user are finished
check_jobs() {
    condor_q >$condorQFile
    
    # Extract the line containing the user's job information
    userJobsLine=$(grep "Total for $username" $condorQFile)
    
    #example userJobLine:
    #Total for username: 10 jobs; 0 completed, 2 removed, 3 idle, 1 running, 4 held, 0 suspended
    
    # Extract the number of jobs from the user's line
    leftJobs=$(echo "$userJobsLine" | awk '{print $4}')
    idleJobs=$(echo "$userJobsLine" | awk '{print $10}')
    runningJobs=$(echo "$userJobsLine" | awk '{print $12}')
    heldJobs=$(echo "$userJobsLine" | awk '{print $14}')
    
    echo "$userJobsLine"
    # Check if all jobs are finished
    if [ "$leftJobs" = "0" ]; then
        echo "Jobs for user $username are finished! -------------------> Let's merge the root files!"
        return 1
        # Check if the jobs are running for too long
        elif [ $runningTime -gt $runningTimeLimit ]; then
        echo "Jobs for user $username are running for too long! Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1
        # Check if there are held jobs
        elif [ "$heldJobs" != "0" ] && [ "$idleJobs" = "0" ] && [ "$runningJobs" = "0" ]; then
        echo "There are held jobs for user $username. Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1
        
    else
        # If there are still running jobs, wait for 10 seconds
        echo "Some jobs are still running for user $username. Please wait."
        echo "Sleep $sleepTime seconds..."
        
        runningTime=$((runningTime + sleepTime)) # Increment runningTime by 10 seconds
        # output percentage of time passed
        echo "Percentage of time passed: $((runningTime * 100 / runningTimeLimit))%"
        echo " "
        sleep $sleepTime
    fi
}

# Function to send and monitor job batches
process_jobs() {
    local start=$1
    local end=$2
    echo "Sending jobs from $start to $end..."
    star-submit -r "$start-$end" sched*.session.xml
    echo "Jobs from $start to $end are sent!"
    echo "========================================================================================================================"
    
    echo "Checking the job status..."
    while :; do
        check_jobs || break
    done
    
    rm "$condorQFile"
    runningTime=0
}

# Function to merge the root files
merge_files() {
    echo " "
    echo "Merging root files..."

    setup 64b || { echo "Failed to setup 64b"; exit 1; }
    setup root 6.20.08 || { echo "Failed to setup ROOT"; exit 1; }

    # Ensure output directory exists
    mkdir -p ${outputDir}/merged

    # Merge in batches if file count is too high
    file_count=$(ls ${outputDir}/*.root | wc -l)
    
    if [[ $file_count -gt 5000 ]]; then
        echo "Too many files (${file_count}). Merging in smaller batches..."
        
        batch_size=5000
        batch_number=1
        temp_files=""
        
        for batch in $(seq 0 $batch_size $file_count); do
            batch_output="${outputDir}/merged/batch_${batch_number}.root"
            hadd -f -k -j 4 ${batch_output} $(ls ${outputDir}/*.root | head -n $batch_size)
            temp_files+="${batch_output} "
            ((batch_number++))
        done
        
        # Merge final batches
        echo "Merging all batches into final output..."
        hadd -f -k -j 4 ${outputFile} ${temp_files}
        
    else
        echo "Merging all files at once..."
        hadd -f -k -j 4 ${outputFile} ${outputDir}/*.root
    fi
}



merge() {
    local filename='productionList.txt'
    local folder=$(tail -n 1 "$filename")
    local outDate=$(date +%d-%m)
    local details="corrnHitsFit2Poss"
    local fileName="outInclE_${outDate}_${details}.root"

    echo "Merging root files into ${fileName}..."
    ./haddMulti.pl "$fileName" output/production/"$folder"/*.root
}

# Main loop to do the analysis
clean
echo "Output is cleaned!"
echo "========================================================================================================================"

echo "simulating jobs..."
./runJob.sh
echo "simulation is finished!"
echo "========================================================================================================================"

# Process job batches
process_jobs 0 1499
process_jobs 1500 2999
process_jobs 3000 4499
process_jobs 4500 5999
process_jobs 6000 7143

echo "All jobs are finished!"
echo "Checking which jobs are missing..."
./checkMissingFiles.sh
echo "All missing jobs are resubmitted!"
echo "========================================================================================================================"
check_jobs
echo "All jobs are finished!"
echo "========================================================================================================================"

echo "Merging all output files..."
merge
echo "All output files are merged!"
