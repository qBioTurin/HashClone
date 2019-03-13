#ifndef __MAP_H__
#define __MAP_H__
#include <map>
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

#define MAXSIZE 1048
#define OUTPUT 0

using namespace std;


//It divides the reads in the input files in two files so that reads in the same pair-end are separated.
//It check the read name to discovered error in the pair-end. I.E. pair-end with more then 2 reads.
int main(int argc, char **argv) {

time_t time_1,time_4;

if (argc<4)
        {
        std::cerr<<"\n\nUSE: Parser <input_file> <out_file1> <out_file2>\n\n";
        exit(EXIT_FAILURE);
        }

time(&time_1);

ifstream in(argv[1],ifstream::in);
if(!in)
        {
          cout<<"*****Error opening the first input file *****\n\n";
        exit(EXIT_FAILURE);
        }

ofstream out1(argv[2],ofstream::out);
if(!out1) 
        {
          cout<<"*****Error opening the first output file *****\n\n"; 
          exit(EXIT_FAILURE); 
        }

ofstream out2(argv[3],ofstream::out);
if(!out2)
        {
        cout<<"*****Error opening the second  output file *****\n\n";
        exit(EXIT_FAILURE);
        }

cout<<"\n\nSTART EXECUTION..."<<endl;

char buffer[MAXSIZE];
unsigned int line=0;

map <string,int> reads;

while (!in.eof()){

      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
            int num=in.gcount();
      if (buffer[num-1]!='\0')
      {
        buffer[num]='\0';
        num++;
      }
     if(buffer[0]!='\0')
        {
        line++;
      if (line%1000001==0)
        cout<<"Processing line:"<<line<<endl;
      if ((buffer[0]!='>')&&(buffer[0]!='@'))
        {
        cout<<"Error at line: "<<line<<endl;
        cout<<buffer<<endl;
         exit(EXIT_FAILURE);
        }
        if (reads.find(std::string(buffer))!=reads.end())
        {
        if (reads[std::string(buffer)]!=1)
                {
                cout<<"Error pair end at line:"<<line<<" "<<buffer<< " "<<reads[std::string(buffer)] <<endl;
                //exit(EXIT_FAILURE);
                }
        reads[std::string(buffer)]++;
#if OUTPUT
        out2<<buffer<<"/2"<<endl;
#endif
        in.getline(buffer,MAXSIZE);
#if OUTPUT
        out2<<buffer<<endl;
#endif
        in.getline(buffer,MAXSIZE);
#if OUTPUT
        out2<<buffer<<endl;
#endif
        in.getline(buffer,MAXSIZE);
#if OUTPUT
        out2<<buffer<<endl;
#endif
        line=line+3;
        }
        else
        {
         reads[buffer]=1;
#if OUTPUT
         out1<<buffer<<"/1"<<endl;
#endif
         in.getline(buffer,MAXSIZE);
#if OUTPUT
         out1<<buffer<<endl;
#endif
         in.getline(buffer,MAXSIZE);
#if OUTPUT
         out1<<buffer<<endl;
#endif
         in.getline(buffer,MAXSIZE);
#if OUTPUT
         out1<<buffer<<endl;
#endif
         line=line+3;
        }
}
}
    out1.close();
    out2.close();
    in.close();
        map<std::string,int>::iterator it=reads.begin();
        int count=0;
    while (it!=reads.end())
        {
        if (it->second!=2 )
                {
                cout<< it->first<<" "<<it->second<<endl;
                count++;
                }
        it++;
        }
        cout<<"reads rep. "<<count<<endl;
    time(&time_4);
    cout<<"\n\nEND EXECUTION"<<endl;
    cout<<"\nResults are saved in: "<<argv[2]<<" "<<argv[3]<<endl;
    cout<<"\n=========================== TIME ===========================\n\n\t";
    cout<<"Total time required: "<<(time_4-time_1)<<"s."<<endl;
    cout<<"\n=========================== TIME ===========================\n\n";
    return EXIT_SUCCESS;
}
