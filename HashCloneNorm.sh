#!/bin/bash

kmer=$1
hash=$2
coll=$3
threashold=$4
output=$(perl -e "use File::Spec; print(File::Spec->rel2abs(\"${5}\"),\"\n\")")
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
inputpath=( "$@" )
inputpathlength=${#inputpath[@]}
inputname=( "$@" )

echo "Files:"
for (( i=7; i<${arraylength}+1; i++ ));
do
        IFS='/' read -a inputpath <<< "${array[$i-1]}" 					
        IFS='.' read -a inputname <<< "${inputpath[-1]}"	
	echo ${inputname[0]}
	
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
for (( i=7; i<${arraylength}+1; i++ ));
do
   IFS='/' read -a inputpath <<< "${array[$i-1]}"	
   IFS='.' read -a inputname <<< "${inputpath[-1]}"
	
  echo "./HashCheckerFreqRPM  ${output}/input-${inputname[0]}- fastq 1 ${output}/input-${inputname[0]}- fastq 1 ${output}/output-${inputname[0]} $kmer 1 $hash $coll"
  ${path}/HashCheckerFreqRPM  ${output}/input-${inputname[0]}- fastq 1 ${output}/input-${inputname[0]}- fastq 1 ${output}/output-${inputname[0]} $kmer 1 $hash $coll
  if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: HashCheckerFreqRPM: output is save in ${output}"
    echo
    if [ "$mail" != "null" ] 
    then
        mail -s "Program terminated with an error" -r $mail $mail <<< "Error output is save in ${output}"
    fi 
    
    exit 0
  fi
  
  freq+="${output}/output-${inputname[0]}.freq "
done

echo

echo "./CompCheckerKmerRPM    ${output}/KmerSignificativi0.kmer $threashold $numf $freq"
${path}/CompCheckerKmerRPM    ${output}/KmerSignificativi0.kmer $threashold $numf $freq
if [ $? -ne 0 ]
  then
    echo
    echo "Failed command: CompCheckerKmerRPM: output is save in ${output} "
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
  IFS='/' read -a inputpath <<< "${array[$i-1]}"
  IFS='.' read -a inputname <<< "${inputpath[-1]}"
  
  echo "./HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-${inputname[0]}- fastq 1 ${output}/signature-in-${inputname[0]} $kmer 1 $hash $coll"
  ${path}/HashCheckerSignature ${output}/KmerSignificativi  kmer  1 ${output}/input-${inputname[0]}- fastq 1 ${output}/signature-in-${inputname[0]} $kmer 1 $hash $coll
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
 
  sig+="${output}/signature-in-${inputname[0]}.signature "
done
echo

echo "./CompCheckerRead ${output}/Clone.csv $threashold $numf  $sig"
${path}/CompCheckerRead ${output}/Clone.csv $threashold $numf  $sig
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
   IFS='/' read -a inputpath <<< "${array[$i-1]}"	
   IFS='.' read -a inputname <<< "${inputpath[-1]}"

   rm ${output}/input-${inputname[0]}-0.fastq
   echo  rm ${output}/input-${inputname[0]}-0.fastq
done
echo
echo "Program  finished successfully: output is save in ${output}"
echo
if [ "$mail" != "null" ]
then
    mail -s "Program  finished successfully" -r $mail $mail <<< "Output is save in ${output}"
fi
