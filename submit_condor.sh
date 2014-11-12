###Condor submitting template for plain root jobs. Run on all the files separately in a given folder, the DATASET folder below
###TA-WEI WANG, 02/20/2014 created
###Your plain root file needs to be modified , see example loop.C
###Please compile the .C before submission to make sure your code is working.
###Checking condor jobs status: condor_q <username> 

###Plain root .C to be run
CONFIGFILE="loop.C"

###All the header/related files needed
TRANSFERFILE="loop.C,loop.h"

###Folder location within which files are to be run
DATASET=/mnt/hadoop/cms/store/user/wangj/HI_Btuple/20140218_PAMuon_HIRun2013_PromptReco_v1/*

###Output file location
DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/bntuple/root_condor/fout

###Output file name
OUTFILE="test_output"

###Maximum number of files to be run
MAXFILES=3

###Log file location and it's name
LOGDIR=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/bntuple/root_condor/logout
LOGNAME=testrootcondor

########################## Create subfile ###############################
dateTime=$(date +%Y%m%d%H%M)
fileCounter=0
INFILE=""
mkdir -p $DESTINATION
mkdir -p $LOGDIR

for file in $DATASET
do
    if [ $fileCounter -ge $MAXFILES ]
    then
	break
    fi

    INFILE="$file"
    fileCounter=$((fileCounter+1))

# make the condor file
cat > subfile <<EOF

Universe = vanilla
Initialdir = .
Executable = exec_condor.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ${OUTFILE}_${fileCounter}.root $INFILE
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.out
Error        = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.err
Log          = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $TRANSFERFILE

Queue
EOF

############################ Submit ###############################

#cat subfile
condor_submit subfile
mv subfile $LOGDIR/$LOGNAME-$dateTime-$fileCounter.subfile
done
echo "Submitted $fileCounter jobs to Condor."
