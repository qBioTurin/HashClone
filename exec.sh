#!/bin/sh

kmer=$1
hash=$2
coll=$3
threashold=$4
output=$5
mail=$6
path=${PWD}
echo
echo "======================================================="
echo "			INPUT PARAMETERS"
echo "======================================================="
echo
echo "Kmer size: $kmer"
echo "Hash table size: $hash"
echo "Collision list size: $coll"
echo "Threashold: $threashold"
echo "Output folder: $output"
 if [ "$mail" != "null" ]
 then
      echo "mail: $mail"  
 else
      echo "mail: DISABLE"
 fi  
array=( "$@" )
arraylength=${#array[@]}
echo "Files:"
for (( i=7; i<${arraylength}+1; i++ ));
do
   echo "     ${array[$i-1]}"
done
let numf=arraylength-6;
echo "Num. of files: $numf"
echo "PATH: ${PWD}"
echo "$@"
echo 
echo "======================================================="

echo 
rm ${output}/*.hash
rm ${output}/*.freq
rm ${output}/*.kmer
rm ${output}/*.cvs
rm ${output}/*.cou*
rm ${output}/*.col*
rm ${output}/*.sign*




echo
for (( i=7; i<${arraylength}+1; i++ ));
do
   ln -s ${array[$i-1]} ${output}/input-$i-0.fastq
   echo " ln -s ${array[$i-1]} ${output}/input-$i-0.fastq"
done
echo

freq=""
for (( i=7; i<${arraylength}+1; i++ ));
do
  echo "./HashCheckerFreq  ${output}/input-$i- fastq 1 ${output}/input-$i- fastq 1 ${output}/output-$i $kmer 1 $hash $coll"
  ${path}/HashCheckerFreq  ${output}/input-$i- fastq 1 ${output}/input-$i- fastq 1 ${output}/output-$i $kmer 1 $hash $coll
  if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: HashCheckerFreq: output is save in ${output}"
    echo
    if [ "$mail" != "null" ] 
    then
        mail -s "Program terminated with an error" -r $mail $mail <<< "Error output is save in ${output}"
    fi 
    
    exit 0
  fi
  
  freq+="${output}/output-$i.freq "
done
 
echo

echo "./CompCheckerKmer    ${output}/KmerSignificativi0.kmer $threashold $numf $freq"
${path}/CompCheckerKmer    ${output}/KmerSignificativi0.kmer $threashold $numf $freq
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
for (( i=7; i<${arraylength}+1; i++ ));
do
  echo "./HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-$i- fastq 1 ${output}/signature-in-$i $kmer 1 $hash $coll"
  ${path}/HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-$i- fastq 1 ${output}/signature-in-$i $kmer 1 $hash $coll
    if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: HashCheckerSignature:  output is save in ${output}"
    echo
     if [ "$mail" != "null" ]
     then
        mail -s "Program terminated with an error" -r $mail $mail <<< "Error output is save in ${output}"
     fi   
    exit 0
  fi
 
  sig+="${output}/signature-in-$i.signature "
done
echo

echo "./CompCheckerRead ${output}/Clone.cvs $threashold $numf  $sig"
${path}/CompCheckerRead ${output}/Clone.cvs $threashold $numf  $sig
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


for (( i=7; i<${arraylength}+1; i++ ));
do
   rm ${output}/input-$i-0.fastq
   echo  rm ${output}/input-$i-0.fastq
done
echo
echo "Program  finished successfully: output is save in ${output}"
echo
if [ "$mail" != "null" ]
then
    mail -s "Program  finished successfully" -r $mail $mail <<< "Output is save in ${output}"
fi
