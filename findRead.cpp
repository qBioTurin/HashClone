#ifndef __SET_H__
#define __SET_H__
#include <set>
#endif

#ifndef __FST_H__
#define __FST_H__
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

#ifndef __GNR_H__
	#define __GNR_H__
	#include "general.h"

#endif

#define MAXSIZE 2048


using namespace std;
using namespace general;

//It takes in input two files. The first contains a list of read names, while the second a set of reads in fasq format.
//It stores in the output file the reads with a name contained in the first file. 

int main(int argc, char **argv) {

time_t time_1,time_4;

if (argc<4)
        {
        std::cerr<<"\n\nUSE: FindRead <input_file1> <input_file2> <delimitator_name> <out_file> \n\n";
        exit(EXIT_FAILURE);
        }

time(&time_1);

ifstream in(argv[1],ifstream::in);
if(!in)
        {
          cout<<"*****Error opening the first input file *****\n\n";
        exit(EXIT_FAILURE);
        }
        
ifstream in1(argv[2],ifstream::in);
if(!in1)
        {
        cout<<"*****Error opening the second input file *****\n\n";
        exit(EXIT_FAILURE);
        }
               
ofstream out(argv[4],ofstream::out);
if(!out) 
        {
          cout<<"*****Error opening the  output file *****\n\n"; 
          exit(EXIT_FAILURE); 
        }



cout<<"\n\nSTART EXECUTION..."<<endl;

char buffer[MAXSIZE];
unsigned int line=0;
char* delimC = argv[3];
Parser parser;


set <string> read_name;

while (!in.eof()){

      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
            int num=in.gcount();
      if (buffer[num-1]!='\0'){
        buffer[num]='\0';
        num++;
      }
      if(buffer[0]!='\0'){
        line++;
	if (line%1001==0)
	  cout<<"Processing line:"<<line<<endl;
	parser.update(delimC,buffer);
	read_name.insert(parser.get(0));
      }
}
in.close();
cout<<"Total number of read names:"<<read_name.size()<<endl;


while (!in1.eof()){

      buffer[0]='\0';
      in1.getline(buffer,MAXSIZE);
            int num=in.gcount();
      if (buffer[num-1]!='\0'){
        buffer[num]='\0';
        num++;
      }
      if(buffer[0]!='\0'){
	if (read_name.find(buffer)!=read_name.end()){
	  out<<buffer<<endl;
	  in1.getline(buffer,MAXSIZE);//nucleotides
	  out<<buffer<<endl;
	  in1.getline(buffer,MAXSIZE);//+
	  out<<buffer<<endl;
	  in1.getline(buffer,MAXSIZE);//qualities
	  out<<buffer<<endl;
	}
	else
	{
	    in1.getline(buffer,MAXSIZE);//nucleotides
	    in1.getline(buffer,MAXSIZE);//+
	    in1.getline(buffer,MAXSIZE);//qualities
	}
        line++;
	if (line%100001==0)
	  cout<<"Processing reads:"<<line<<endl;

      }
}
in1.close();
out.close();

time(&time_4);
cout<<"\n\nEND EXECUTION"<<endl;
cout<<"\nMatching reads are saved in: "<<argv[3]<<endl;
cout<<"\n=========================== TIME ===========================\n\n\t";
cout<<"Total time required: "<<(time_4-time_1)<<"s."<<endl;
cout<<"\n=========================== TIME ===========================\n\n";
return EXIT_SUCCESS;
}
