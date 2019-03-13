  #ifndef __IOS_H__
  #define __IOS_H__
  #include <iostream>
  #endif 
  
  #ifndef __FST_H__
  #define __FST_H__
  #include <fstream>
  #endif 

  #ifndef __STDL__
  #define __STDL__
  #include <stdlib.h>
  #endif

  #ifndef __GEN_H__
	#define __GEN_H__
	#include "general.h"
  #endif

  #define MAXSIZE 1048


  using namespace std;

  int main(int argc, char **argv) {

  time_t time_1,time_4;
  int **freq;

  if (argc<2)
  {
    std::cerr<<"\n\nUSE: test_freq <input_file> <output_file> <number_probability_distribution>\n\n";
    exit(EXIT_FAILURE);
  }

  time(&time_1);

  ifstream in(argv[1],ifstream::in);
  if(!in)
  {
    cout<<"*****Error opening the input file "<<argv[1]<<"*****\n\n";
    exit(EXIT_FAILURE);
  }

  ofstream out(argv[2],ofstream::out);
  if(!out) 
  {
    cout<<"*****Error opening the  output file "<<argv[2]<<" *****\n\n"; 
    exit(EXIT_FAILURE); 
  }

  int num_prob=atoi(argv[3]);

  if (num_prob<=0)
  {
    cout<<"*****Error the number of probability distributions must be greater than 0 *****\n\n";
    exit(EXIT_FAILURE); 
  }

  //create the frequency vector 4 X num_prop
  freq = (int **) malloc(sizeof(int*)*num_prob+1);


  for (int i=0;i<num_prob+1;i++)
  {
    freq[i]=(int*) malloc(sizeof(int)*4);
  }

  char buffer[MAXSIZE];
  unsigned int line=0;

  while (!in.eof())
  {
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
      if ((buffer[0]=='>')|| (buffer[0]=='@'))
      {
	class general::Parser parser("*-",buffer);
	if (parser.size()<3)
	{
	  freq[num_prob][0]++;
	  
	}
	else
	{
	  unsigned int tmp=2;
	  while(tmp<parser.size())
	  {
	    int index = atoi (parser.get(tmp).c_str());
	    if ((index<0)||(index>num_prob+1))
	    {
	      cout<<"*****Error the index of read "<<buffer<<" is wrong ("<<index<<")*****\n\n";
	      exit(EXIT_FAILURE); 
	    }
	    char c= parser.get(tmp+1)[0];
	    switch(c)
	    {
	      case 'A':
	      case 'a':
		freq[index][0]++;
		break;
	      case 'C':
	      case 'c':
		freq[index][1]++;
		break;
	      case 'G':
	      case 'g':
		freq[index][2]++;
		break;
	      case 'T':
	      case 't':
		freq[index][3]++;
		break; 
	      default :
		cout<<"*****Error the SNP of read "<<buffer<<" is wrong ("<< c <<")*****\n\n";
		exit(EXIT_FAILURE); 
	      break;
	    }
	    tmp=tmp+3;
	  }
	
	}
       
      }
    }
  }
  //free memory
  for (int i=0;i<num_prob;i++)
  {
    out<<i<<";"<<freq[i][0]<<";"<<freq[i][1]<<";"<<freq[i][2]<<";"<<freq[i][3]<<endl;  
    free(freq[i]);
  }

  free(freq);
  out.close();
  time(&time_4);
  cout<<"\n\nEND EXECUTION"<<endl;
  cout<<"\nResults are saved in: "<<argv[2]<<endl;
  cout<<"\n=========================== TIME ===========================\n\n\t";
  cout<<"Total time required: "<<(time_4-time_1)<<"s."<<endl;
  cout<<"\n=========================== TIME ===========================\n\n";
  return EXIT_SUCCESS;
  }