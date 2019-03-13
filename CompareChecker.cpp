#include "CompareChercker.h"
#include <fstream>
#include <string>




using namespace std;

double  TAU=1.0;
std::string type;
bool    spike_in_condition=false;

#if !KMER
multimap<string, class Info > freq;
#else
 #if RPM
  map<string,double*> freq;
 #else 
  map<string,int*> freq;
 #endif
#endif
int num_input_file=0;
//! It is used for edit distance.
double  d[2][MAXSIZE];
#if GUI
ofstream outGUI;
#endif

#if !KMER
inline bool Edistance(const string& T,const string& S,string& r_T){
  int missD=0,missR=0;
  int size=min3(S.size(),T.size(),r_T.size());
  for(int i=0;i<size&&(missD<THEDIT||missR<THEDIT);++i){
    if (T[i]!=S[i]) missD++;
    if (r_T[i]!=S[i]) missR++;
  }
  //cout<<missD<<" "<<missR<<endl;
  if ((missD<THEDIT)||(missR<THEDIT))
    return true;
  else
    return false;
}



bool AminoDistance(const string& T,const string& S,const string& r_T){
 
  int size=min3(S.size(),T.size(),r_T.size());
  int amini_freq[3][5] {0};
  for(int i=0;i<size;++i){
       switch(S[i]){
	case 'A':
	 ++amini_freq[0][0];
	break;
	case 'T':
	 ++amini_freq[0][1];
	break;
	case 'C':
	 ++amini_freq[0][2];
	break;     
	case 'G':
	 ++amini_freq[0][3];
	break;
	case 'N':
	 ++amini_freq[0][4];
	break;
      }
      switch(T[i]){
	case 'A':
	 ++amini_freq[1][0];
	break;
	case 'T':
	 ++amini_freq[1][1];
	break;
	case 'C':
	 ++amini_freq[1][2];
	break;     
	case 'G':
	 ++amini_freq[1][3];
	break;
	case 'N':
	 ++amini_freq[1][4];
	break;
      }
      switch(r_T[i]){
	case 'A':
	 ++amini_freq[2][0];
	break;
	case 'T':
	 ++amini_freq[2][1];
	break;
	case 'C':
	 ++amini_freq[2][2];
	break;     
	case 'G':
	 ++amini_freq[2][3];
	break;
	case 'N':
	 ++amini_freq[2][4];
	break;
      }
  }
  int max_1=0,max_2=0;
  for(int i=0;i<5;++i){
    if ( amini_freq[0][i]>amini_freq[1][i])
      max_1+=amini_freq[0][i]-amini_freq[1][i];
    if ( amini_freq[0][i]>amini_freq[2][i])
      max_2+=amini_freq[0][i]-amini_freq[2][i];
  }
  return (max_1<THEDIT||max_2<THEDIT)?true:false;
}    
   




bool LevenshteinDistance(const string& s1,const string& s2){
  const unsigned int m(s1.size());
  const unsigned int n(s2.size());
 
  if (( m==0 )||( n==0 ))  exit(EXIT_FAILURE);
 
  size_t costs[MAXSIZE];

 
  for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
  size_t i = 0;
  for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
  {
    costs[0] = i+1;
    size_t corner = i;
 
    size_t j = 0;
    for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
    {
      size_t upper = costs[j+1];
      if( *it1 == *it2 )
      {
		  costs[j+1] = corner;
	  }
      else
	  {
		size_t t(upper<corner?upper:corner);
        costs[j+1] = (costs[j]<t?costs[j]:t)+1;
	  }
 
      corner = upper;
    }
  }
  int diff=abs((int)m-(int)n);
  if (diff > MINSIGNATURE)
    diff=0;
  return (costs[n]-diff)<THEDIT?true:false;
}
 
bool SmithWaterman(const string& S,const string& T,int& score){

    // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  int size=0;
  if (S.size()>T.size()){
    aligner.Align(T.c_str(), S.c_str(), S.size(), filter, &alignment);
    size=T.size();
  }
  else{
     aligner.Align(S.c_str(), T.c_str(), T.size(), filter, &alignment);
     size=S.size();
  }

   return ((score=alignment.sw_score)>=(size*THEDIT)?true:false);
}


double SmithWaterman(const string& S,const string& T,double& matches,int& positionBclone, int& positionEclone,int& positionBigh, int& positionEigh,string& cigar){
    // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner(2,4,10,1);
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  //Qui è stata switchata la S con la T
  aligner.Align(S.c_str(), T.c_str(), T.size(), filter, &alignment);

  //compute matches
  int csize= alignment.cigar.size();
  matches=0;
  double tot=0;
  for (int i=0;i<csize;++i){
    if(cigar_int_to_op(alignment.cigar[i])=='='){
      matches+=(double)cigar_int_to_len(alignment.cigar[i]);
      tot+=(double)cigar_int_to_len(alignment.cigar[i]);
      }
   // tot+=(double)cigar_int_to_len(alignment.cigar[i]);
    else
      if (cigar_int_to_op(alignment.cigar[i])=='X')
	tot+=(double)cigar_int_to_len(alignment.cigar[i]);
	
  }
  matches=(matches/tot)*100;
  //matches=(matches/S.size())*100;
  positionBclone=alignment.ref_begin;
  positionEclone=alignment.ref_end;
  positionBigh=alignment.query_begin;
  positionEigh=alignment.query_end;
  cigar=alignment.cigar_string;
  return alignment.sw_score;
}

#if RAF
void search(const string& S, const string& T,int freqInc,int nfile){
    SmithWaterman
    string r_T,r_S; 
    int size=T.size();
    
    for(int i=size-1;i>=0;--i){
      switch(T[i]){
	case 'A':
	 r_T.push_back('T');
	break;
	case 'T':
	  r_T.push_back('A');
	break;
	case 'C':
	  r_T.push_back('G');
	break;     
	case 'G':
	  r_T.push_back('C');
	break;
	case 'N':
	  r_T.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }
    
    size=S.size();
    for(int i=size-1;i>=0;--i){
      switch(S[i]){
	case 'A':
	 r_S.push_back('T');
	break;
	case 'T':
	  r_S.push_back('A');
	break;
	case 'C':
	  r_S.push_back('G');
	break;     
	case 'G':
	  r_S.push_back('C');
	break;
	case 'N':
	  r_S.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }
    
   bool found=false;
   //auto it=freq.end();
   for (auto it=freq.begin();it!=freq.end();++it){
     if (S.size()>it->first.size()){
       if (((S.find(it->first)!=std::string::npos)||(r_S.find(it->first)!=std::string::npos))){//&& (LevenshteinDistance(T,it->second.first)||LevenshteinDistance(it->second.first,r_T))){
	 found=true;
	 it->second.second[nfile]+=freqInc;
       }
     }
     else
     {
       if (((it->first.find(S)!=std::string::npos)||(it->first.find(r_S)!=std::string::npos))){//&&(LevenshteinDistance(T,it->second.first)||LevenshteinDistance(it->second.first,r_T))){
	found=true;
	it->second.second[nfile]+=freqInc;
       }
     }	 
   }
   if (!found){
     int *tmp=(int*)malloc(sizeof(int)*num_input_file); 
	 for (int j=0;j<num_input_file;++j)
	    tmp[j]=0;
     tmp[nfile]=freqInc;
     freq.insert(std::pair<string,pair <string,int*>>(S,std::pair<string,int*>(T,tmp)));
     
   }
     
}
#else
void search(const string& S, const string& T,int freqInc,int nfile){
    
    string r_T,r_S; 
    int size=T.size();
    
    for(int i=size-1;i>=0;--i){
      switch(T[i]){
	case 'A':
	 r_T.push_back('T');
	break;
	case 'T':
	  r_T.push_back('A');
	break;
	case 'C':
	  r_T.push_back('G');
	break;     
	case 'G':
	  r_T.push_back('C');
	break;
	case 'N':
	  r_T.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }
    
    size=S.size();
    for(int i=size-1;i>=0;--i){
      switch(S[i]){
	case 'A':
	 r_S.push_back('T');
	break;
	case 'T':
	  r_S.push_back('A');
	break;
	case 'C':
	  r_S.push_back('G');
	break;     
	case 'G':
	  r_S.push_back('C');
	break;
	case 'N':
	  r_S.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }
   int score=0; 
   auto fit=freq.end();
   //auto it=freq.end();
   for (auto it=freq.begin();it!=freq.end();++it){
      int score1=0,score2=0,score3=0,score4=0;
      if ((SmithWaterman(S,it->first,score1)||SmithWaterman(r_S,it->first,score2))&&(SmithWaterman(T,it->second.clone,score3)||SmithWaterman(r_T,it->second.clone,score4))){
	  if ((score1>score||score2>score)){
	    fit=it;
	    score =max(score1,score2);
	  }
      }
   }
   if (fit==freq.end()){
     vector <int>tmp (num_input_file,0);
     tmp[nfile]=freqInc;
     freq.insert(std::pair<string,class Info>(S, Info(T,tmp)));
   }
   else{
     fit->second.freq[nfile]+=freqInc;
   }
     
}


#endif

#endif



void  save_data(ofstream &out, char **argv){

   
#if !KMER
  char delim[]="/";
  class general::Parser parser; 
  //Column header in textual file
  out<<"Signature\tClone\t";
  for (int i=0;i<num_input_file;i++)
	      {
		parser.update(delim,string(argv[i]));
		out<<parser.get(parser.size()-1)<<"\t";
	      }
  if(spike_in_condition != true){
    if(type=="IGH"){
        out<<"V-GENE\tV-G-ID%\tJ-GENE\tJ-G-ID%\tD-GENE\tD-G-ID%\tReverse\tV-ReadPosB\tV-ReadPosE\tV-IghPosB\tV-IghPosE\tV-Cigar\tJ-ReadPosB\tJ-ReadPosE\tJ-IghPosB\tJ-IghPosE\tJ-Cigar\tD-ReadPosB\tD-ReadPosE\tD-IghPosB\tD-IghPosE\tD-Cigar\n";
    }
    if(type=="IGK"){
        out<<"V-GENE\tV-G-ID%\tJ-GENE\tJ-G-ID%\tReverse\tV-ReadPosB\tV-ReadPosE\tV-IghPosB\tV-IghPosE\tV-Cigar\tJ-ReadPosB\tJ-ReadPosE\tJ-IghPosB\tJ-IghPosE\tJ-Cigar\n";
    }    
  }
  
  else{
     if(type=="IGH"){ 
        out<<"V-GENE\tV-G-ID%\tJ-GENE\tJ-G-ID%\tD-GENE\tD-G-ID%\tSPIKE-IN\tReverse\tV-ReadPosB\tV-ReadPosE\tV-IghPosB\tV-IghPosE\tV-Cigar\tJ-ReadPosB\tJ-ReadPosE\tJ-IghPosB\tJ-IghPosE\tJ-Cigar\tD-ReadPosB\tD-ReadPosE\tD-IghPosB\tD-IghPosE\tD-Cigar\n";
     }
     if(type=="IGK"){
        out<<"V-GENE\tV-G-ID%\tJ-GENE\tJ-G-ID%\tSPIKE-IN\tReverse\tV-ReadPosB\tV-ReadPosE\tV-IghPosB\tV-IghPosE\tV-Cigar\tJ-ReadPosB\tJ-ReadPosE\tJ-IghPosB\tJ-IghPosE\tJ-Cigar\n";
     }    
    }
  
  
  #if GUI
//creating columns  
 outGUI<<"[\n { \"columns\": [\n\t{\n\t\t\"text\": \"ID\",\n\t\t\"datafield\": \"id\",\n\t\t\"width\": \"10\"\n\t},\n\t{\n\t\t\"text\": \"Clone\",\n\t\t\"datafield\": \"clone\",\n\t\t\"width\": \"100\"\n\t},\n\t{\n\t\t\"text\": \"Signature\",\n\t\t\"datafield\": \"sig\",\n\t\t\"width\": \"100\"\n\t}";

//for frequencies
for (int i=0;i<num_input_file;i++)
	      {
	       parser.update(delim,string(argv[i]));
	       outGUI<<",\n\t{\n\t\t\"text\": \""<<parser.get(parser.size()-1)<<"\",\n\t\t\"datafield\": \"f"<<i+1<<"\",\n\t\t\"width\": \"80\"\n\t}";
	      }
outGUI<<",\n\t{\n\t\t\"text\": \"V-GENE\",\n\t\t\"datafield\": \"vgene\",\n\t\t\"width\": \"110\"\n\t}";    
outGUI<<",\n\t{\n\t\t\"text\": \"V-G ID%\",\n\t\t\"datafield\": \"vid\",\n\t\t\"width\": \"75\"\n\t}";
outGUI<<",\n\t{\n\t\t\"text\": \"J-GENE\",\n\t\t\"datafield\": \"jgene\",\n\t\t\"width\": \"110\"\n\t}";
outGUI<<",\n\t{\n\t\t\"text\": \"J-G ID%\",\n\t\t\"datafield\": \"jid\",\n\t\t\"width\": \"75\"\n\t}";
if(type=="IGH"){
outGUI<<",\n\t{\n\t\t\"text\": \"D-GENE\",\n\t\t\"datafield\": \"dgene\",\n\t\t\"width\": \"110\"\n\t}";
outGUI<<",\n\t{\n\t\t\"text\": \"D-G ID%\",\n\t\t\"datafield\": \"did\",\n\t\t\"width\": \"110\"\n\t}";
}

if(spike_in_condition == true){  
outGUI<<",\n\t{\n\t\t\"text\": \"SPIKE-IN\",\n\t\t\"datafield\": \"spike_in\",\n\t\t\"width\": \"40\"\n\t}";
}

outGUI<<",\n\t{\n\t\t\"text\": \"Reverse\",\n\t\t\"datafield\": \"reverse\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"V-ReadPosB\",\n\t\t\"datafield\": \"vreadposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"V-ReadPosE\",\n\t\t\"datafield\": \"vreadposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"V-IghPosB\",\n\t\t\"datafield\": \"vighposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"V-IghPosE\",\n\t\t\"datafield\": \"vighposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"V-Cigar\",\n\t\t\"datafield\": \"vcigar\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"J-ReadPosB\",\n\t\t\"datafield\": \"jreadposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"J-ReadPosE\",\n\t\t\"datafield\": \"jreadposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"J-IghPosB\",\n\t\t\"datafield\": \"jighposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"J-IghPosE\",\n\t\t\"datafield\": \"jighposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"J-Cigar\",\n\t\t\"datafield\": \"jcigar\",\n\t\t\"width\": \"110\"\n\t}"; 
if(type=="IGH"){
outGUI<<",\n\t{\n\t\t\"text\": \"D-ReadPosB\",\n\t\t\"datafield\": \"dreadposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"D-ReadPosE\",\n\t\t\"datafield\": \"dreadposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"D-IghPosB\",\n\t\t\"datafield\": \"dighposB\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"D-IghPosE\",\n\t\t\"datafield\": \"dighposE\",\n\t\t\"width\": \"110\"\n\t}"; 
outGUI<<",\n\t{\n\t\t\"text\": \"D-Cigar\",\n\t\t\"datafield\": \"dcigar\",\n\t\t\"width\": \"110\"\n\t}";
}
outGUI<<"\n\t]},\n { \"rows\": [\n";	

int count=1;      
  #endif 
#endif
  
   for (auto it=freq.begin();it!=freq.end();++it)
   {    
     bool end_search=false;
     for (int i=1;(i<num_input_file)&&!(end_search);i++){
       double log1=0;
#if KMER       
       if (it->second[i-1]!=0)
	 log1=log10(it->second[i-1]);
#else
       if (it->second.freq[i-1]!=0)
	 log1=log10(it->second.freq[i-1]);
#endif       
       double log2=0;
#if KMER          
       if (it->second[i]!=0)
       log2=log10(it->second[i]);
#else 
       if (it->second.freq[i]!=0)
       log2=log10(it->second.freq[i]);
#endif       
       double difflog=fabs(log1-log2);
       if (((difflog>=TAU)||isinf(difflog))){
#if KMER	 
	// if (((difflog>=TAU)||isinf(difflog))){
        out<<">Signature_"<<i<<endl;
        out<<it->first<<endl<<"+"<<endl;
        for (int j=0;j<num_input_file;++j)
        out<<(double)it->second[j]<<" ";
        out<<endl;
        
#else	
   
  //  if ((it->second.tmp[0] > 0) || (it->second.tmp[1] > 0) || (it->second.tmp[2] > 0) ){ 
     
	 out<<it->first<<"\t"<<it->second.clone;
	 out<<"\t"<<it->second.freq[0];
	 for (int i=1;i<num_input_file;i++){
        out<<"\t"<<it->second.freq[i];
	  }
	  
	  //VDJ POSITION 
     
     if(type=="IGH"){
     for (int i=0;i<3;i++){
        out<<"\t"<<it->second.info_identity[i];
        out<<"\t"<<it->second.info_pos[i];
       //    out<<"\t"<<it->second.tmp[i];
     }
     }
     
     if(type=="IGK"){
     for (int i=0;i<2;i++){
        out<<"\t"<<it->second.info_identity[i];
        out<<"\t"<<it->second.info_pos[i];
       //    out<<"\t"<<it->second.tmp[i];
     }
     }
         
    
       
            
     if(spike_in_condition == true){
           
        if(it->second.spike_in == true){
            out<<"\t"<< it->second.spike_in_name;
             
        }      
        else{
            out<<"\t"<<"Null";
        }
            
      }
      
      
    out<<"\t"<< it->second.align[1].reverse_align;  
    
    if(type=="IGH"){
      for (int i=0;i<3;i++){
        out<<"\t"<<it->second.align[i].positionBclone<<"\t"<<it->second.align[i].positionEclone;
        out<<"\t"<<it->second.align[i].positionBigh<<"\t"<<it->second.align[i].positionEigh;
        out<<"\t"<<it->second.align[i].cigar;
     }
    }
    
    if(type=="IGK"){
      for (int i=0;i<2;i++){
        out<<"\t"<<it->second.align[i].positionBclone<<"\t"<<it->second.align[i].positionEclone;
        out<<"\t"<<it->second.align[i].positionBigh<<"\t"<<it->second.align[i].positionEigh;
        out<<"\t"<<it->second.align[i].cigar;
     }
    }
         
    
	out<<endl;
 
//}  
  
#if GUI
	if (count>1)
	  outGUI<<",\n";
//	if ((it->second.info_identity[i] != null_position) && (it->second.info_pos[i] != 0)){    
      outGUI<<"\t{\"id\" : \""<<count<<"\", \"clone\" : \""<<it->second.clone<<"\", \"sig\" :\""<<it->first<<"\"";
	//for frequencies
	for (int i=0;i<num_input_file;i++){
        outGUI<<", \"f"<<i+1<<"\":"<<it->second.freq[i]; 
	}
	
	//IGH information
    outGUI<<", \"vgene\" : \""<<it->second.info_identity[0]<<"\""; 
    outGUI<<", \"vid\" : "<<it->second.info_pos[0];
    outGUI<<", \"jgene\" : \""<<it->second.info_identity[1]<<"\""; 
    outGUI<<", \"jid\" : "<<it->second.info_pos[1];
    if(type=="IGH"){
    outGUI<<", \"dgene\" : \""<<it->second.info_identity[2]<<"\""; 
    outGUI<<", \"did\" : "<<it->second.info_pos[2];
    }

    if(spike_in_condition == true){   
        if(it->second.spike_in ==true)
            outGUI<<", \"spike_in\" : \""<<it->second.spike_in_name<<"\"";
        else
            outGUI<<", \"spike_in\" : "<<"\"Null"<<"\"";
    }
   // }
    
        outGUI<<", \"reverse\" : \""<< it->second.align[1].reverse_align<<"\"";  
    
      
        outGUI<<", \"vreadposB\" : \""<<it->second.align[0].positionBclone<<"\"";
        outGUI<<", \"vreadposE\" : \""<<it->second.align[0].positionEclone<<"\"";
        outGUI<<", \"vighposB\" : \""<<it->second.align[0].positionBigh<<"\"";
        outGUI<<", \"vighposB\" : \""<<it->second.align[0].positionEigh<<"\"";
        outGUI<<", \"vcigar\" : \""<<it->second.align[0].cigar<<"\"";
        
        outGUI<<", \"jreadposB\" : \""<<it->second.align[1].positionBclone<<"\"";
        outGUI<<", \"jreadposE\" : \""<<it->second.align[1].positionEclone<<"\"";
        outGUI<<", \"jighposB\" : \""<<it->second.align[1].positionBigh<<"\"";
        outGUI<<", \"jighposB\" : \""<<it->second.align[1].positionEigh<<"\"";
        outGUI<<", \"jcigar\" : \""<<it->second.align[1].cigar<<"\"";
        if(type=="IGH"){
        outGUI<<", \"dreadposB\" : \""<<it->second.align[2].positionBclone<<"\"";
        outGUI<<", \"dreadposE\" : \""<<it->second.align[2].positionEclone<<"\"";
        outGUI<<", \"dighposB\" : \""<<it->second.align[2].positionBigh<<"\"";
        outGUI<<", \"dighposB\" : \""<<it->second.align[2].positionEigh<<"\"";
        outGUI<<", \"dcigar\" : \""<<it->second.align[2].cigar<<"\"";
        }
    
	double val =(it->second.freq[0]!=0)?log10(it->second.freq[0]):0.0;
	outGUI<<", \"points\": { \"data\" : [ [1,"<<val<<"]";
	for (int i=1;i<num_input_file;i++){ 
	 val =(it->second.freq[i]!=0)?log10(it->second.freq[i]):0.0; 
	outGUI<<", ["<<i+1<<", "<<val<<"]"; 
	}
	outGUI<<"], \"label\": \"Clone "<<count<<"\"} }";
	++count;	
#endif	
	
#endif
	 end_search=true;
      }
     }
   }  

#if GUI
outGUI<<"\n\t]}\n]\n";
#endif
}


#if !KMER
 inline void reverse(string& S,string& r_S){
 int size=S.size();   
 for(int i=size-1;i>=0;--i){
      switch(S[i]){
	case 'A':  
	 r_S.push_back('T');
	break;
	case 'T':
	  r_S.push_back('A');
	break;
	case 'C':
	  r_S.push_back('G');
	break;     
	case 'G':
	  r_S.push_back('C');
	break;
	case 'N':
	  r_S.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }

}
    

void compute_position(void){
  
    
  std::string germfile[3] { "./Germline/IGHV.fa", "./Germline/IGHJ.fa",  "./Germline/IGHD.fa"};
  int numGerm=3;
  
  
  if(type=="IGK"){
    germfile[0] = "./Germline/IGKV.fa";
    germfile[1] = "./Germline/IGKJ.fa";
    germfile[2] = "";
    numGerm=2; 
  }
  
  clock_t startGlobal,endGlobal;

  set< pair<string,string> > IGH[numGerm];
 
  startGlobal=clock();

  for (int i=0;i<numGerm;++i){
   
        cout<<"Reading file: "<< germfile[i]<<endl; 
    
    ifstream f_IGH(germfile[i],ifstream::in);
    if(!f_IGH) {
      cerr << "\n*****Error opening input file "<< germfile[i] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    std::string id("");
    std::string seq("");
    char delim[]="|";
    class general::Parser parser;
    bool pending=false;
    while (!f_IGH.eof()){ 
	std::string buffer("");
        getline(f_IGH,buffer);
	if ((buffer[0]=='>')||(buffer=="")){
	  if (seq!=""){
	   IGH[i].insert(make_pair(id,seq)); 
	   seq="";
	   pending=false;
	  }
	  if (buffer!=""){
	    parser.update(delim,buffer);
	    if (parser.size()<2){
	      cerr << "\n*****Error parsing ID: "<<buffer<< "*****" << endl;
	      exit(EXIT_FAILURE);
	    }
	    id=parser.get(1)+"_"+parser.get(3);
	  }
	  }
	else{
	  seq+=buffer;
	  pending=true;
	  }
    }
    if (pending){
       IGH[i].insert(make_pair(id,seq)); 
    }
    f_IGH.close();
  }
  cout<<"Mapping putative clone: "<<endl; 
  for (auto itc=freq.begin();itc!=freq.end(); ++itc){//for all clone 
  //  if (itc->second.spike_in == false){ 
   
    double Maxscore[numGerm] {-999999, -999999, -999999};
    double MaxscoreR[numGerm] {-999999, -999999, -999999};
    
    //For alignement information
    align_info MaxAlign[numGerm];
    align_info MaxAlignR[numGerm];
    std::string cigar;
    std::string Rcigar;
    
    double Max_matches[numGerm] {0,0,0}, Max_matchesR[numGerm] {0,0,0};
    set< pair<string,string> >::iterator Maxit[numGerm]  {IGH[0].end()};
    set< pair<string,string> >::iterator MaxitR[numGerm]  {IGH[0].end()}; 
    
    
    int max_position[numGerm]  {0,0,0};
    int max_positionR[numGerm] {0,0,0};
    
    std::string cloneR("");
    std::string clone(itc->second.clone);
    
   
    reverse(itc->second.clone,cloneR); 
    for (int i=0;i<numGerm;i++){//for all IGH 
      for (auto it=IGH[i].begin();it!=IGH[i].end();++it){//for all germline
	double matches=0,matchesR=0, score=0, scoreR=0;
	int positionBclone=0, positionBRclone=0, positionEclone=0, positionERclone=0;
	int positionBigh=0, positionBRigh=0, positionEigh=0, positionERigh=0;
	if (clone.size()>MINSIZEALIGN)
	  score=SmithWaterman(clone,it->second,matches,positionBclone, positionEclone,positionBigh, positionEigh,cigar);
	if (cloneR.size()>MINSIZEALIGN)
	  scoreR=SmithWaterman(cloneR,it->second,matchesR,positionBRclone, positionERclone,positionBRigh, positionERigh,Rcigar);
	if (score>Maxscore[i]){
	  Maxscore[i]=score;
	  Maxit[i]=it;
	  Max_matches[i]=matches;
      
      //For alignement information
      MaxAlign[i].cigar=cigar;
      MaxAlign[i].positionBclone=i!=0?positionBclone+MaxAlign[0].positionEclone:positionBclone;
      MaxAlign[i].positionEclone=i!=0?positionEclone+MaxAlign[0].positionEclone:positionEclone;
      MaxAlign[i].positionBigh=positionBigh;
      MaxAlign[i].positionEigh= positionEigh;

      
      if (i==0)//case IGHV
	    max_position[i]=positionEclone;
	  else //IGHJ
	    max_position[i]=positionBclone; 
	} 
	if (scoreR>MaxscoreR[i]){
	  MaxscoreR[i]=scoreR;
	  MaxitR[i]=it;
	  Max_matchesR[i]=matchesR;
      
      //For alignement information
      MaxAlignR[i].cigar=Rcigar;
      MaxAlignR[i].positionBclone=i!=0?positionBRclone+MaxAlignR[0].positionEclone:positionBRclone;
      MaxAlignR[i].positionEclone=i!=0?positionERclone+MaxAlignR[0].positionEclone:positionERclone;
      MaxAlignR[i].positionBigh=positionBRigh;
      MaxAlignR[i].positionEigh= positionERigh;
      
      
	  if (i==0)//case IGHV
	    max_positionR[i]=positionERclone;
	  else //IGHJ
	    max_positionR[i]=positionBRclone;  
	} 
      }
     if (i==0){//next IGHJ
       clone=clone.erase(0, max_position[i]);
       cloneR=cloneR.erase(0,max_positionR[i]);
     }
     else{//next IGHD
       clone=clone.substr(0,max_position[i]);
       cloneR=cloneR.substr(0,max_positionR[i]); 
     }
    }
    //select the best mapping 
     double sum=0.0,sumR=0.0;
     for (int i=0;i<numGerm;i++){
       sum+=Maxscore[i];
       sumR+=MaxscoreR[i];
     }
     if (sum>=sumR){
      for (int i=0;i<numGerm;i++){
          
	
        //For alignement information  
        itc->second.align.push_back(MaxAlign[i]);
        //modified by Greta
        if (Maxscore[i] <= 0){ 
            itc->second.info_pos.push_back(0.0);  
            itc->second.info_identity.push_back(Maxit[i]->first);
            itc->second.tmp.push_back(Maxscore[i]);
            
        }
        else{
            itc->second.tmp.push_back(Maxscore[i]);
            itc->second.info_pos.push_back(Max_matches[i]);
            itc->second.info_identity.push_back(Maxit[i]->first);
	    }      
      }
    }
    else{
      for (int i=0;i<numGerm;i++){  
        //For alignement information
        MaxAlignR[i].reverse_align=true;  
        itc->second.align.push_back(MaxAlignR[i]);
        //modified by Greta      
		if (MaxscoreR[i] <= 0){
           itc->second.info_pos.push_back(0.0);  
           itc->second.info_identity.push_back(MaxitR[i]->first);
           itc->second.tmp.push_back(MaxscoreR[i]);
        }
        else{
            itc->second.tmp.push_back(MaxscoreR[i]);
            itc->second.info_identity.push_back(MaxitR[i]->first);
            itc->second.info_pos.push_back(Max_matchesR[i]);  
        }    
        }
      }

   // }
  }
 endGlobal=clock();   
 cout<<"\n=========================== TIME ===========================\n\n";
 cout<<"\tTime to identify alignment: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n============================================================\n\n";

  }



  
void Spike_in_Research(ifstream &in){
    

     int count=0;
     std::cout<<"\n\n =========================================================\n";
     std::cout<<"|	      	        CompChecker       	          |\n";
     std::cout<<"|	      	       Spike In Research  	          |\n";
     std::cout<<" =========================================================\n\n";
     std::cout << "Starting the research of Spike In" << std::endl;  
     
     std::string  spike_in_id="";
     while( !in.eof()){
            std::string buffer2("");
            getline(in,buffer2);
            if((buffer2!="") && (buffer2[0]=='>')){
                    spike_in_id=std::string(buffer2.begin()+1, buffer2.end());
                    }            
            else{ 
               if(((buffer2!="") || buffer2[0]!='>' || buffer2[0]!='@')&& (buffer2[0]=='A' || buffer2[0]=='a' || buffer2[0]=='C' || buffer2[0]=='c' || buffer2[0]=='G' || buffer2[0]=='g' || buffer2[0]=='T' ||     buffer2[0]=='t' || buffer2[0]=='N' || buffer2[0]=='n' )){
                   for (auto k=freq.begin();k!=freq.end(); ++k){
                       if (k->second.spike_in != true){ //it is TRUE iff I have already found it
                            int score=0, scoreR=0;
                            std::string cloneR("");
                            reverse(k->second.clone,cloneR); 
                            if (SmithWaterman(buffer2,k->second.clone,score)||(SmithWaterman(buffer2,cloneR,scoreR)))
                                    if (score > MINSCORE || scoreR > MINSCORE ){
                                            k->second.spike_in = true;
                                            k->second.spike_in_name=spike_in_id;
                                     //       std::cout <<  "Spike:" << k->second.spike_in_name << std::endl;
                                    //      std::cout <<  "Score:" << score << "\nScore R:" << scoreR << std::endl;
                                            count++;
                                        }
                          }
                    }
               
                }
          
            }
     }
     std::cout<<"\n Found:"<< count << "\t Spike In" <<"\n" <<std::endl; 
    std::cout<<" =========================================================\n"<<std::endl;
}
#endif


int main(int argc, char **argv) {
  clock_t startGlobal,endGlobal;
  startGlobal=clock();
  
  char delim[]=" \t";
  class general::Parser parser;
 
  cout<<"\n\n =========================================================\n";
  cout<<"|	      	        CompChecker       	          |\n";
#if !KMER 
  cout<<"|	      	       READs Frequency    	          |\n";
#else  
  cout<<"|	      	       Kmers Frequency    	          |\n";
#endif
  cout<<" =========================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n\n";

  if (argc<7){      
    std::cerr<<"\n\nUSE: CompChecker <path/Outputfile> <TAU> <n> <path/InputFile2> ... <path/InputFileN> <path/spike_in_file> \n\n";
    std::cerr<<"\t<path/Outputfile>:\t output file\n";
    std::cerr<<"\t<TAU>:\t\t\t threshold used to select significant k-mer identification (default 1.0)\n";
    std::cerr<<"\t<type>:\t\t\t IGH or IGK\n";
    std::cerr<<"\t<n>:\t\t\t number of input files (Must be at least 2)\n";
    std::cerr<<"\t<path/InputfileI>:\t ith input file\n";
    std::cerr<<"\t<path/Spike_In_File>:\t path/Spike_In_File\n\n\n";
    exit(EXIT_FAILURE);
  }
  
  TAU=atof(argv[2]);
  type = argv[3];
  num_input_file = atoi(argv[4]);
  
  if((type!="IGH") && (type!="IGK")){ 
    std::cerr<<"\tError: you must specify IGH or IGK\n";
     exit(EXIT_FAILURE);
  }    
      
  if (num_input_file<2){
     std::cerr<<"\t<n>Error: number of input files must be at least 2)\n";
     exit(EXIT_FAILURE);
  }
	
  for (int i=0;i<num_input_file;i++){
   ifstream in(argv[i+5],ifstream::in);
   if(!in){
      cerr << "\n*****Error opening input file "<< argv[i+5] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    int proc=0;    
    while (!in.eof()){
	std::string buffer("");
        getline(in,buffer);
	if (buffer!=""){
#if KMER	  
	if (++proc%1000000==1)
#else
	if (++proc%1000==1)  
#endif		  
	  cout<<"Processed seq. "<<proc<<" "<<"Size: "<<freq.size()<<endl;
	parser.update(delim,buffer);
#if KMER	
	if (parser.size()<2){
#else
	 if (parser.size()<3){ 
#endif	  
	  //comment to do format wrong
	  cerr<<"\n*****Error wrong file format"<< argv[i+5] <<" *****" << endl;
	 exit(EXIT_FAILURE); 
	}
#if KMER
	auto it=freq.find(parser.get(0));
	if (it==freq.end()){
	 #if RPM
   	 freq[parser.get(0)]=(double*)malloc(sizeof(double)*num_input_file);
	 #else
	 freq[parser.get(0)]=(int*)malloc(sizeof(int)*num_input_file); 
	 #endif
         for (int j=0;j<num_input_file;++j)
	    freq[parser.get(0)][j]=0;
	}
	#if RPM
	freq[parser.get(0)][i]=atof(parser.get(1).c_str());
	#else	
	freq[parser.get(0)][i]=atoi(parser.get(1).c_str());
	#endif	
#else
      	
      search(parser.get(0),parser.get(2),atoi(parser.get(1).c_str()),i);
#endif	
	}
      }
      
  in.close();
  cout<<"Size:"<<freq.size()<<endl;
  }

    //open output file 
    ofstream out(argv[1],ofstream::out);
    if(!out) 
    {  
      cerr << "\n*****Error opening output file "<< argv[1] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
#if GUI   
    string filename=string(argv[1])+".gui";
    outGUI.open(filename.c_str(),ofstream::out);
     if(!outGUI) 
    {
      cerr << "\n*****Error opening output file "<< argv[1] <<" for GUI *****" << endl;
      exit(EXIT_FAILURE);
    } 
#endif  

#if !KMER  
   //perform research of spike-in
      if(argc>5+num_input_file){
        spike_in_condition = true; 
        std::ifstream in;
        in.open (argv[5+num_input_file],std::ifstream::in);
        if(!in){
            cerr << "\n*****Error opening input file "<< argv[5+num_input_file] <<" *****" << endl;
            exit(EXIT_FAILURE);
        }
      
   
     Spike_in_Research(in); 
     in.close();
    }
#endif 

#if !KMER && GERM
   //perform mapping with germlines

    compute_position();

#endif
   //save read or kmer with their frequencies 
   save_data(out,(argv+5)); 
   out.close();
#if GUI
    outGUI.close();
#endif   
   
   endGlobal=clock();

   cout<<"\n\nEND EXECUTION"<<endl;

cout<<"\n=========================== TIME ===========================\n\n\t";
cout<<"Total time required: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
cout<<"\n=========================== TIME ===========================\n\n";
}
