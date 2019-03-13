






using namespace std;

#ifndef __FSTREAM__
	#define __FSTREAM__
	#include <fstream>
#endif

#ifndef __IOS_H__
	#define __IOS_H__
	#include <iostream>
#endif 


#ifndef __STDL__
	#define __STDL__
	#include <stdlib.h>
#endif

#ifndef __STR_H__
	#define __STR_H__
	#include <string.h>

#endif

#ifndef __GEN_H__
	#define __GEN_H__
	#include "general.h"
#endif
 
#ifndef __TIM_H__
	#define __TIM_H__
	#include <time.h>

#endif

#ifndef __MAP_H__
	#define __MAP_H__
	#include <map>

#endif

#ifndef RES_H_
  #define RES_H
  #include <sys/resource.h>
#endif

#ifndef MTH_H_
  #define MTH_H_
  #include <math.h>   
#endif


int main(int argc, char **argv) {
  clock_t startGlobal,endGlobal;
  startGlobal=clock();
  
  char delim[]=" \t";
  class general::Parser parser;
 
  cout<<"\n\n =========================================================\n";
  cout<<"|	      	        UpdateFreq.       	          |\n";
  cout<<"|	      	       Kmers Frequency    	          |\n";
  cout<<" =========================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n\n";

  if (argc<4){
    std::cerr<<"\n\nUSE: UpdateFreq  <path/InputFile> <path/Outputfile>   <value>\n\n";
    std::cerr<<"\t<path/Outputfile>:\t output file\n";
    std::cerr<<"\t<path/Inputfile>:\t input file\n\n\n";
      std::cerr<<"\t<value>:\t a value which will be used to moltiply the current frequencies\n\n\n";
    exit(EXIT_FAILURE);
  }
  cout<<"\n\n___________________________________________________________\n\n\t\t\tInput parameters \n___________________________________________________________\n\n";
   cout<<"\n\t\t Input name: "<<argv[1];
   cout<<"\n\t\t Output name: "<<argv[2];
   cout<<"\n\t\t value: "<<argv[3]<<endl<<endl;
  cout<<"___________________________________________________________\n";
   double val=atof(argv[3]);
   ifstream in(argv[1],ifstream::in);
   if(!in){
      cerr << "\n*****Error opening input file "<< argv[1] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    ofstream out(argv[2],ifstream::out);
    if(!out){
      cerr << "\n*****Error opening output file "<< argv[2] <<" *****" << endl;
      exit(EXIT_FAILURE);
    } 
    
    int proc=0;    
    while (!in.eof()){
	std::string buffer("");
        getline(in,buffer);
	if (buffer!=""){
	  if (++proc%1000000==1)	  
	    cout<<"Processed seq. "<<proc<<endl;
	  parser.update(delim,buffer);	
	  if (parser.size()<2){  
	    //comment to do format wrong
	    cerr<<"\n*****Error wrong file format"<< argv[1] <<" *****" << endl;
	  exit(EXIT_FAILURE); 
	  }
	int tmp=  (int)(atoi(parser.get(1).c_str())*val);
	if (tmp<=0)
	  tmp=1;
	out<<parser.get(0)<<"\t"<<tmp<<endl;
	}
    }
    out.close();
    in.close();

   
   endGlobal=clock();

   cout<<"\n\nEND EXECUTION"<<endl;

cout<<"\n=========================== TIME ===========================\n\n\t";
cout<<"Total time required: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"\n=========================== TIME ===========================\n\n";
}