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

#ifndef __SET_H__
	#define __SET_H__
	#include <set>

#endif

#ifndef RES_H_
  #define RES_H
  #include <sys/resource.h>
#endif

#define DNASIZE 2846426792
#define MAXSIZE 1024
#define READSIZE 100

using namespace std;

unsigned short* FreqVec;



int main(int argc, char **argv) {
  clock_t startGlobal,endGlobal;
  startGlobal=clock();
  bool pos=false;
  cout<<"\n\n =========================================================\n";
  cout<<"|	      	        FreqChecker       	          |\n";
  cout<<"|	      	       READs Frequency    	          |\n";
  cout<<" =========================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

    if ((argc!=6)&&(argc!=8))
  {
    std::cerr<<"\n\nUSE: FreqChecker <path/InputFile>  <path/OutputFile> <id_chromosome>  <threshold> <window_size> [Option]\n\n\tOption\n\t\t-f <file_position>\t-> chromosome pos_1 \\n ..\\n pos_n\n\n";
    exit(EXIT_FAILURE);
  }
 
 ifstream fpos;
 
  if (argc==8)
  {
    pos=true;
    fpos.open(argv[7],ifstream::in);
    if(!fpos) 
    {
      cerr << "\n*****Error opening input file "<< argv[7] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    
  }
ifstream in(argv[1],ifstream::in);
    if(!in) 
    {
      cerr << "\n*****Error opening input file "<< argv[1] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
ofstream out(argv[2],ofstream::out);
    if(!out) 
    {
      cerr << "\n*****Error opening output file "<< argv[2] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    

std::string filter(argv[3]);
set<int> lpos;

if (pos)
{
   getline(fpos,filter);
   cout<<filter<<endl;
   while (!fpos.eof())
    {
	std::string buffer("");
        getline(fpos,buffer);
	if (buffer!="")
	lpos.insert(atoi(buffer.c_str()));
    }
}


unsigned int threshold = (unsigned int) (atoi(argv[4]));

unsigned int window_size= (unsigned int)(atoi(argv[5]));

cout<<"\n++++++++++++++++++++++++++++\n";
cout<<"Window size is: "<<window_size<<endl;
cout<<"++++++++++++++++++++++++++++\n\n";

FreqVec =(unsigned short*)malloc(sizeof(unsigned short)*DNASIZE);
memset (FreqVec,0,sizeof(unsigned short)*DNASIZE);

char buffer[MAXSIZE];
char delim[]=" \t";
class general::Parser parser;
int count=0;
 //For each read in the pool file
    while (!in.eof())
    {
      buffer[0]='\0';
      count++;
      if (count%1000000==1)
      {
	cout<<"Reads: "<<count<<endl;
      }
      //Sequence
      in.getline(buffer,MAXSIZE);
      int num=in.gcount();
      if (buffer[num-1]!='\0')
      {
	buffer[num]='\0';
	num++;
      }
      //cout<<buffer<<endl;
      if ((buffer[0]!='\0')&&(buffer[0]!='@'))
      {
       
       parser.update(delim,buffer);
       if ((parser.get(2)==filter)&&((atoi(parser.get(1).c_str())==0)||(atoi(parser.get(1).c_str())==16)))
       {
 	unsigned int start= (unsigned int)atoi(parser.get(3).c_str()); 
	//cout<<start<<endl;
	//exit(1);
	if (start+READSIZE>DNASIZE)
	{
	 cerr << "\n*****Error vector size  is lower than "<< start+READSIZE <<" *****" << endl;
         exit(EXIT_FAILURE);
	}
	else
	{
	  if ((atoi(parser.get(1).c_str())==0))
	    for (unsigned int i=start; i<=start+READSIZE;i++)
	    {
	    FreqVec[i]++;
	    }
	  else
	  {
	    for (unsigned int i=start-READSIZE; i<=start;i++)
	    {
	    FreqVec[i]++;
	    }
	  }
	}
       }
      }

    }
    
if(!pos)
{
//print frequency greater than threshold
unsigned int sum=0;
for (unsigned int i=0;i<DNASIZE;i++)
{
 // if (FreqVec[i]>threshold)
 //   out<<i<<"\t"<<FreqVec[i]<<endl;
  if (i % window_size != 0)
  {
    sum = sum +FreqVec[i];
  }
  else
  {
    sum = sum +FreqVec[i];
    if (sum != 0)
    {
      if (sum>threshold)
      {
	out<<i<<"\t"<<sum<<endl;
      }
      else
	out<<i<<"\t"<<0<<endl;
    sum=0;  
    }
  }

}
}
else
{
  set<int>::iterator it=lpos.begin();
  while (it!=lpos.end())
  {
    out<<*it<<"\t"<<FreqVec[*it];
   /* int sum =0;
    for (int i=*it-50;i<=*it+50;i++){
      if ((i>=0)&&(i<DNASIZE))
      {
	sum+=FreqVec[i];
	//cout<<i<<"\t"<<FreqVec[i];
      }
    }
      out<<"\t("<<sum<<")"<<endl;
      */
    ++it;
  }
}

 endGlobal=clock();
 cout<<"\n=========================== TIME ===========================\n\n";
 cout<<"\tTime to compute frequencies: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n============================================================\n\n";
 free(FreqVec);
 out.close();
 in.close();
}
    