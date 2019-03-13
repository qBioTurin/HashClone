INSTALLING AND USING HashClone SUITE
======================================================

1.Build Requirements
===================================

To compile/install the package you need:
1) GNU make utility;
2) GNU Compiler Collection (version 4.8 or later);
3) Internet connection to download IGMT database (http://www.imgt.org). To use the IMGT germline databases (IMGT/GENE-DB), you have to agree to IMGT license: academic research only, provided that it is referred to IMGT®, and cited as "IMGT®, the international  ImMunoGeneTics information system®".
4) Bash shell 
5) perl
6) Java SE Runtime Environment 8 (http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
   
To compile and install the package go through the following steps:

 - get the zipped archive hashclone.XXX-src.tar.gz
 
 - create a new directory where you want the tool to reside, for example:  
    mkdir ./HashClone
 
 
 - extract the archive into the installation directory
     mv ./hashclone.XXX-src.tar.gz  ./HashClone
     cd ./HashClone
     tar xvfz hashclone.XXX-src.tar.gz 
     cd ./HashCheck
 - type
      make 
    
      
At the end of compiling process the following programs are created:

 		1. HashChecker: it is the standard version to compare reads: all the reads with at least a k-mer not present in the Hash Table are reported in output.
 		
 		2. HashCheckerDebugOutput: it is the debug output version of HashChecker.

 		3. HashCheckerFreq: it computes and returns k-mer frequencies.

 		4. HashCheckerSimilarRead: it  improves the standard version in the output: the reads for which all k-mers are in the hash table are also reported in output (.sim file). 

 		5. HashCheckerSignature: it  computes and returns the signatures of reads. For each reads it computes the signatures as the largest sequence of consecutive k-mers found in the Hash Table.

 		6. HashCheckerPairedEnd: it  takes as input paired-end reads (e.g. <file>.R1.<num>.<ext>,  <file>.R2.<num>.<ext>) and returns the paired-end in which at least one mate has a k-mer not presented in the Hash Table.

 		7. HashCheckerPairedEndSimilar: it  takes as input paired-end reads (e.g. <file>.R1.<num>.<ext>,  <file>.R2.<num>.<ext>) and it  improves the standard version in the output: the reads for which all k-mer of all mathes are in the hash table are also reported in output (.sim). 

		 8. FreqChecker: it  takes in input a set of reads with their positions in terms of chromosome and relative position, and it returns a vector showing the reads distribution (and related peaks).

		 9. Parser: it divides the reads in the input files into two files so that reads in the same pair-end are separated. It checks the read name to discovered error in the pair-end (i.e. pair-end with more han 2 reads.)
		 
 		10. CreateKmer: it  creates all the k-mers from a set of reads. It is able to manage pair-end reads.
		
		11. CompCheckerKmer: it compares k-mers among different set of experiments and returns those are different frequencies.
 
 		12. CompCheckerRead: it compares reads  among different set of experiments  and returns CVS and GUI files with a frequency table. 
 
 		13. FindRead: it  takes in input two files. The first contains a list of read names, while the second a set of reads in ``fasq" format. It stores in the output file the reads with a name contained in  the first  file.
 		
 		14. Alignment: it  computes similarity between two input strings using Smith and Waterman algorithm or Levenshtein Distance.
 
 		15. AlignmentBlast: it  computes similarity between two input strings using Smith and Waterman or LevenshteinDistance. It visualizes the alignment in Blast format.
 
		16. AlignmentIGH  computes similarity between two input string and IGH references
 
 		17. FASTA2FASTQ.pl converts a file with FASTA format into a new file with FASTQ format (use:  perl FASTA2FASTQ.pl file.fa > file.fastq)

Observe that folder "./HackClone/Data" contains the instructions on how to execute these  tools on  simple (provided) datasets  and the corresponding outputs.

2.Getting Started
===================================

- To run HashClone GUI you must type the following commands from a terminal:

    bash HashCloneGUI.sh

Observe that HashCloneGUI requires   java  installed at system level. 


- To run HashClone without GUI you must type the following commands from a terminal:

    bash HashCloneParall.sh <k-mer>  <hash_size> <collision_list_size> <threshold> <type of IG> <output_folder> <mail> <spike-in> <input_file_1> ....<input_file_n>
  

where:

	<k-mer> is the size of k-mers encoded in the hash table.  This must be a value between 1 and 32. 


    <hash_size> is a (prime) number indicating the size of the hash table. Increasing this value reduces the execution time but increases the memory utilization.  Ideally, this value should be close to the number of different k-mers stored in the hash table;

    <collision_list_size> the maximum number of different k-mers that the tool might need to store in the hast table. This parameter is required to optimize the memory utilization.  

    <threshold> this value is the  threshold used to select significant k-mers;
    
    <type of IG>  Immunoglobulin Heavy Chain (IGH) or Immunoglobulin Light chain (IGK);
    
    <output_folder> the folder where the output will be saved;
    
    <mail>  null or a valid mail address where to send information on  termination status of the run;
    
    <spike-in> the fasta file containing the list of the spike in sequences that you want to research in the samples
    
    <input_file_1> ....<input_file_n> the list of input patient files (one for each follow-up in fastq format)


3. An example
===================================
    
To analyze PAT A by terminal assuming "~/output"" as the output folder path, and ~/input as the input folder path where PAT A data are copied:     
    
   bash  HashCloneParall.sh 26   10999997 10999997  1 IGH  ~/output  null null  ~/input/A-S1.fastq  ~/input/AFU1.fastq ~/input/AFU2.fastq ~/input/AFU3.fastq 
  
    
