ifeq ($(shell uname -s),Darwin)
OPTION =   -O3 -Wall  -fmessage-length=0  -std=c++11
else
OPTION =   -O3 -Wall  -fmessage-length=0 -fopenmp -std=c++11
endif
CC = gcc
C++ = g++


all:	HashChecker HashCheckerDebugOutput HashCheckerFreq HashCheckerFreqRPM HashCheckerSimilarRead  HashCheckerSignature  HashCheckerPairedEnd HashCheckerPairedEndSimilar CreateKmer TestFreq  FreqChecker Parser CompCheckerKmer CompCheckerKmerRPM CompCheckerRead FindRead Alignment AlignmentBlast AlignmentIGH UpdateFreq CompCheckerReadNoGerm 

clean:	
	rm  HashChecker  HashCheckerDebugOutput HashCheckerFreq HashCheckerFreqRPM HashCheckerSimilarRead HashCheckerSignature HashCheckerPairedEnd HashCheckerPairedEndSimilar CreateKmer TestFreq FreqChecker Parser CompCheckerKmer CompCheckerRead FindRead Alignment AlignmentBlast AlignmentIGH UpdateFreq CompCheckerReadNoGerm


HashChecker:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -o  HashChecker classIO.cpp  main.cpp  general.cpp     $(OPTION)
	
HashCheckerDebugOutput:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D DEBUGOUTPUT=1 -o  HashCheckerDebugOutput classIO.cpp  main.cpp  general.cpp     $(OPTION)
	

HashCheckerFreq:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D DEBUGFREQ=2 -D RPM=0 -o  HashCheckerFreq classIO.cpp  main.cpp  general.cpp     $(OPTION)

HashCheckerFreqRPM:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D DEBUGFREQ=2 -D RPM=1 -o  HashCheckerFreqRPM classIO.cpp  main.cpp  general.cpp     $(OPTION)

	
HashCheckerSimilarRead:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D SIMILARITY=1 -o  HashCheckerSimilarRead classIO.cpp  main.cpp  general.cpp     $(OPTION)	
	
HashCheckerSignature:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h 
	$(C++) -D SIGNATURE=1 -o  HashCheckerSignature classIO.cpp  main.cpp  general.cpp  ssw_cpp.cpp ssw.c   $(OPTION)	

HashCheckerPairedEnd:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D SIGNATURE=0 -D  PAIREDEND=1 -o  HashCheckerPairedEnd classIO.cpp  main.cpp  general.cpp    $(OPTION)

HashCheckerPairedEndSimilar:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++) -D SIGNATURE=0 -D  PAIREDEND=1 -D  SIMILARITY=1 -o  HashCheckerPairedEndSimilar classIO.cpp  main.cpp  general.cpp    $(OPTION)	
	
FreqChecker:	mainfreq.cpp
	$(C++) -o FreqChecker mainfreq.cpp   general.cpp  $(OPTION)
	
Parser:		mainparser.cpp	
	$(C++) -o Parser mainparser.cpp    $(OPTION)
	
TestFreq:	test_frequency.cpp
	 $(C++) -o TestFreq test_frequency.cpp   general.cpp  $(OPTION)	
	
CreateKmer:	mainKmer.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp 
	$(C++)  -o CreateKmer  classIO.cpp  mainKmer.cpp  general.cpp     $(OPTION)
	
CompCheckerKmer:	  CompareChecker.cpp general.cpp conf.h	general.h
	$(C++)  -D KMER=1 -D RPM=0 -o CompCheckerKmer  CompareChecker.cpp  general.cpp  $(OPTION)

CompCheckerKmerRPM:	  CompareChecker.cpp general.cpp conf.h	general.h
	$(C++)  -D KMER=1 -D RPM=1 -o CompCheckerKmerRPM  CompareChecker.cpp  general.cpp  $(OPTION)
	
CompCheckerRead:	  CompareChecker.cpp general.cpp conf.h	general.h ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h 
	$(C++) -D KMER=0 -D GERM=1 -o CompCheckerRead  CompareChecker.cpp  general.cpp ssw_cpp.cpp ssw.c  $(OPTION)	
	
CompCheckerReadNoGerm:	  CompareChecker.cpp general.cpp conf.h	general.h ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h 
	$(C++) -D KMER=0 -D GERM=0 -o CompCheckerReadNoGerm  CompareChecker.cpp  general.cpp ssw_cpp.cpp ssw.c  $(OPTION)		
	
FindRead:	  findRead.cpp general.cpp	general.h
	$(C++) -o FindRead  findRead.cpp  general.cpp  $(OPTION)
	
Alignment:	  EditDistance.cpp ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h general.cpp general.h
	$(C++) -D BLAST_OUTPUT=0 -o Alignment EditDistance.cpp ssw_cpp.cpp ssw.c general.cpp $(OPTION)  
	
AlignmentBlast:	  EditDistance.cpp ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h general.cpp general.h
	$(C++) -D BLAST_OUTPUT=1 -o AlignmentBlast EditDistance.cpp ssw_cpp.cpp ssw.c general.cpp $(OPTION) 
	
AlignmentIGH:	  EditDistance.cpp ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h general.cpp general.h
	$(C++) -D BLAST_OUTPUT=0 -D  IGH=1 -o AlignmentIGH EditDistance.cpp ssw_cpp.cpp ssw.c general.cpp $(OPTION)  	
	
UpdateFreq:	UpdateFreq.cpp general.cpp conf.h general.h
	$(C++) -o UpdateFreq  UpdateFreq.cpp  general.cpp  $(OPTION)
	
DataBase:	
	cd Germline; ./get_germline.sh
	
doc:     
	doxygen  DoxyFile

	
