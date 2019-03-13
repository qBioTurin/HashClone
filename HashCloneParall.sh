#!/bin/bash

kmer=$1
hash=$2
coll=$3
threashold=$4
type=$5
output=$(perl -e "use File::Spec; print(File::Spec->rel2abs(\"${6}\"),\"\n\")")
mail=$7
spike_in=$8
path=${PWD}

declare -i numf
array=( "$@" )
arraylength=${#array[@]}
numf=arraylength-8




echo
echo "======================================================="
echo "		HashClone Execution"
echo "======================================================="

if [ ! -d "$output" ]; then
  echo "\nWarning: Output folder does not exit --> It will be created"
  mkdir $output
  if [ $? -ne 0 ]; then
    echo "Error creating the output folder"
    exit 1
fi
fi

if [ $numf -lt 2 ]
then
    echo
    echo "Error: you have to specify at least two input file"
    echo
    echo "Command line:"
    echo
    echo "./HashCloneR <k-mer size>  <hash size>  <collision size> <threashold> <IGH/IGK> <output folder> <mail> <spike_in> <input_file_1> ...<input_file_n>"
    echo
    exit 0
fi

if [ "$type" != "IGH" ] && [ "$type" != "IGK" ]
then
    echo
    echo "Error: you have to specify IGH or IGK"
    echo
    echo "Command line:"
    echo
    echo "./HashCloneR <k-mer size>  <hash size>  <collision size> <threashold> <IGH/IGK> <output folder> <mail> <spike_in> <input_file_1> ...<input_file_n>"
    echo
    exit 0
fi

echo
echo "-------------------------------------------------------"
echo "		INPUT PARAMETERS"
echo "-------------------------------------------------------"
echo
echo "Kmer size: $kmer"
echo "Hash table size: $hash"
echo "Collision list size: $coll"
echo "Threashold: $threashold"
echo "Type of IG: $type"
echo "Output folder: $output"
 if [ "$mail" != "null" ]
 then
      echo "mail: $mail"  
 else
      echo "mail: DISABLE"
 fi

if [ "$spike_in" != "null" ]
 then
      echo "spike_in: $spike_in"
 else
      echo "spike_in: DISABLE"
 fi
   

inputpath=( "$@" )
inputpathlength=${#inputpath[@]}
inputname=( "$@" )

echo "Files:"
for (( i=9; i<${arraylength}+1; i++ ));
do
        IFS='/' read -a inputpath <<< "${array[$i-1]}" 					
        IFS='.' read -a inputname <<< "${inputpath[-1]}"	
	echo "    ${inputname[0]}"
	
done

echo "Num. of files: $numf"
echo "Current folder: ${PWD}"
echo 
echo "-------------------------------------------------------"

echo 

echo "Removing  *.hash *.freq *.kmer *.cvs *.cou* *.col* *.sing*"
rm ${output}/*.hash 2> /dev/null
rm ${output}/*.freq 2> /dev/null
rm ${output}/*.kmer 2> /dev/null
rm ${output}/*.cvs  2> /dev/null
rm ${output}/*.cou* 2> /dev/null
rm ${output}/*.col* 2> /dev/null
rm ${output}/*.sign* 2> /dev/null




echo
for (( i=9; i<${arraylength}+1; i++ ));
do
        IFS='/' read -a inputpath <<< "${array[$i-1]}" 					
        IFS='.' read -a inputname <<< "${inputpath[-1]}"	
	 					
   rm  ${output}/input-${inputname[0]}-0.fastq
#   tmp_name=$(readlink -f ${array[$i-1]})
   tmp_name=$(perl -e "use File::Spec; print(File::Spec->rel2abs(\"${array[$i-1]}\"),\"\n\")")
#  ln -s ${array[$i-1]} ${output}/input-${inputname[0]}-0.fastq
   ln -s $tmp_name ${output}/input-${inputname[0]}-0.fastq
#   echo " ln -s ${array[$i-1]} ${output}/input-${inputname[0]}-0.fastq"
   echo "ln -s $tmp_name ${output}/input-${inputname[0]}-0.fastq"
done
echo

freq=""
for (( i=9; i<${arraylength}+1; i++ ));
do
   IFS='/' read -a inputpath <<< "${array[$i-1]}"	
   IFS='.' read -a inputname <<< "${inputpath[-1]}"
	
  echo "./HashCheckerFreq  ${output}/input-${inputname[0]}- fastq 1 ${output}/input-${inputname[0]}- fastq 1 ${output}/output-${inputname[0]} $kmer 1 $hash $coll" 
  ${path}/HashCheckerFreq  ${output}/input-${inputname[0]}- fastq 1 ${output}/input-${inputname[0]}- fastq 1 ${output}/output-${inputname[0]} $kmer 1 $hash $coll &
  pids="$pids $!"
  freq+="${output}/output-${inputname[0]}.freq "
done

#Parall waiting 
wait $pids

echo

echo "./CompCheckerKmer    ${output}/KmerSignificativi0.kmer $threashold $type $numf $freq"
${path}/CompCheckerKmer    ${output}/KmerSignificativi0.kmer $threashold $type $numf $freq
if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: CompCheckerKmer: output is save in ${output} "
    echo 
    if [ "$mail" != "null" ]
    then
        mail -s "Program terminated with an error" -r $mail $mail <<< "Error output is save in ${output}"
    fi    
    exit 0
fi

echo 

sig=""
pids=""

for (( i=9; i<${arraylength}+1; i++ ));
do
  IFS='/' read -a inputpath <<< "${array[$i-1]}"
  IFS='.' read -a inputname <<< "${inputpath[-1]}"
  
  echo "./HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-${inputname[0]}- fastq 1 ${output}/signature-in-${inputname[0]} $kmer 1 $hash $coll"
  ${path}/HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-${inputname[0]}- fastq 1 ${output}/signature-in-${inputname[0]} $kmer 1 $hash $coll &
  pids="$pids $!" 
  sig+="${output}/signature-in-${inputname[0]}.signature "
done

#Parall waiting 
wait $pids

if [ "$spike_in" != "null" ]
 then
  
 echo "./CompCheckerRead ${output}/Clone.csv $threashold $type $numf  $sig $spike_in "  
 ${path}/CompCheckerRead ${output}/Clone.csv $threashold $type $numf  $sig $spike_in
 else
 
   echo "./CompCheckerRead ${output}/Clone.csv $threashold $type $numf  $sig "   
   ${path}/CompCheckerRead ${output}/Clone.csv $threashold $type $numf  $sig
  fi

if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: CompCheckerRead"
    echo
     if [ "$mail" != "null" ]
     then
        mail -s "Program terminated with an error" -r $mail $mail <<< "Error output is save in ${output}"
     fi

    exit 0
fi


echo
echo "Removing temprorary files..."
for (( i=9; i<${arraylength}+1; i++ ));
do
   IFS='/' read -a inputpath <<< "${array[$i-1]}"	
   IFS='.' read -a inputname <<< "${inputpath[-1]}"

   rm ${output}/input-${inputname[0]}-0.fastq
 #  echo  rm ${output}/input-${inputname[0]}-0.fastq
done
echo
echo "Program  finished successfully: output is save in ${output}"
echo
if [ "$mail" != "null" ]
then
    mail -s "Program  finished successfully" -r $mail $mail <<< "Output is save in ${output}"
fi


echo
echo
echo "======================================================="
echo "	              END EXECUTION"
echo "======================================================="
echo
echo
