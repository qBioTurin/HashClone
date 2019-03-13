/***************************************************************************
 *   Copyright (C) 2013 by Marco Beccuti                                   *
 *   beccuti@di.unito.it                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef THR_H_
        #define THR_H_
        #include <pthread.h>
#endif

#ifndef CIO_H_
  #define CIO_H
  #include "classIO.h"
#endif

#ifndef TIM_H_
  #define TIM_H
  #include <sys/time.h>
#endif

#ifndef RES_H_
  #define RES_H
  #include <sys/resource.h>
#endif

#include <stdlib.h>

using namespace IOFn;

bool ContSeach=false;

//int mybitset::ksize=0;

class IOF* p_ref; 

bool pairend=false;


unsigned int* NotFReads;

int main(int argc, char **argv) {
  clock_t startGlobal,endGlobal;
  time_t time1,time2;
  time1=time(NULL);
  startGlobal=clock();
  
  cout<<"\n\n =========================================================\n";
  cout<<"|                     Create K-mer                       |\n";
  cout<<" =========================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

  if ((argc!=8)&&(argc!=5))
  {
    std::cerr<<"\n\nUSE x Pair-end: CreateKmer <path/ReadFile1> <file_extension1> <input_file_number1> <path/ReadFile2>  <file_extension2>  <input_file_number2>  <k-mer_size>\n";
     std::cerr<<"\n\nUSE x Single-end: CreateKmer <path/ReadFile1> <file_extension1> <input_file_number1>  <k-mer_size>\n\n";
    exit(EXIT_FAILURE);
  }
  
 //Initialize input data 
 unsigned int qgram=0;

 if (argc==8)
    pairend=true; 
 
 if (pairend)
    qgram=atoi(argv[7]);
 else
    qgram=atoi(argv[4]); 
 //mybitset::setKmer(qgram*2);


 unsigned int files1=atoi(argv[3]);
 unsigned int files2=0;
  //Initialize input/output files
 
 string ReadFname1(argv[1]);
 string ExtFile1(argv[2]);


 string ReadFname2="";
 string ExtFile2="";
 if  (pairend){
    files2=atoi(argv[6]); 
    ReadFname2=string(argv[4]);
    ExtFile2=string(argv[5]);
 }
 //Initialize input data 
 
 //Initialize input/output files

 
 string IOread="";
  //Initialize input/output files
 
 cout<<"\n\n___________________________________________________________\n\n\t\t\tInput parameters \n___________________________________________________________\n\n";
 
 cout<<"\tInput 1:";
  cout<<"\n\t\tfiles name: "<<ReadFname1;
  cout<<"\n\t\tInput files extension: "<<ExtFile1;
  cout<<"\n\t\t#Input files: "<<files1;
  cout<<"\n\t\t#Output files:" <<ReadFname1<<"Output";
  //cout<<"\n\t\tfiles name: "<<IOread<<endl;
  if (pairend){
    cout<<"\n\n\tInput 2:";
    cout<<"\n\t\tfiles name: "<<ReadFname2;
    cout<<"\n\t\tfiles extension: "<<ExtFile2;
    cout<<"\n\t\t#Input files2: "<<files2;
    cout<<"\n\t\t#Output files:"<<ReadFname2<<"Output";
  }
 cout<<"\n\tWindow's size: "<<qgram<<"\n";
 cout<<"___________________________________________________________\n";

 cout<<"\n\nSTART FORMATTING...\n"<<endl;
 
 

 
  
 
 
 
 //!create variable for deriving the memory utilization  
 int who = RUSAGE_SELF;  
 struct rusage usage;
 //!create variable for deriving the memory utilization  
 
 class IOF ref(ReadFname1,ExtFile1,ReadFname2,ExtFile2,IOread,qgram,files1,files2);
 

 time1=time(NULL);
 startGlobal=clock();
 cout<<"\n\nSTART MAPPING...\n"<<endl;
 ref.ReadXKerm(0,files1-1, p1);
 if (pairend)
 ref.ReadXKerm(0,files2-1, p2);
 //compute memory utilization
 getrusage(who,&usage); 
 //compute memory utilization
 
 cout<<"\n\n\tTotal memory used: "<<usage.ru_maxrss<<"KB"<<endl;


 cout<<"\n\nEND FORMATTING\n"<<endl;

 
 endGlobal=clock();
 time2=time(NULL);
 cout<<"\n=========================== TIME ===========================\n\n";
  cout<<"\tReal Time to  create K-mer: "<<time2-time1<<"s."<<endl;
 cout<<"\tTime to create k-mer: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n============================================================\n\n";
 

 exit(EXIT_SUCCESS);
}
