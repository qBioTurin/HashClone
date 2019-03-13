#include <string>
#include <iostream>
#include <vector>
#include <set>    
#include <fstream>
#include <cstring>
#include "ssw_cpp.h"
#include "ssw.h"

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

#ifndef __GEN_H__
	#define __GEN_H__
	#include "general.h"
#endif


using namespace std;
using std::string;
using std::cout;
using std::endl;
 
//const std::string* idRef {nullptr};
//const std::string* idRead  {nullptr};

 
double  d[2][MAXSIZE];
int input_matches {0};

class align_info{
public:    
    int positionBclone {0};
    int positionEclone {0};
	int positionBigh   {0}; 
    int positionEigh   {0};
    std::string cigar    {""};
    bool reverse_align {false};
};


#if BLAST_OUTPUT
static void ssw_write (StripedSmithWaterman::Alignment& a,
			const string& ref_seq, const string& read_seq,
            const string& idRef, const string& idRead) 
{
	static int count=0;
	cout<<"\n\n\n________________________________________________________________________________\n";
	cout<<"===================================MAPPING "<<++count<<"====================================\n"<<endl;
	
	cout<<"Ref:\n";
	cout<<idRef;
	for (unsigned int i=0;i<ref_seq.size();++i){
	  if (i%80==0)
	    cout<<"\n";
	  cout<<ref_seq[i];
	}
	cout<<"\n\nRead:\n";
	cout<<idRead;
	for (unsigned int i=0;i<read_seq.size();++i){
	  if (i%80==0)
	    cout<<"\n";
	  cout<<read_seq[i];
	}
	       
	cout<<"\n\n================================Info alignment==================================\n"<<endl;
	fprintf(stdout, "\toptimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t\n", a.sw_score, a.sw_score_next_best);
	if (a.ref_begin + 1) fprintf(stdout, "\ttarget_begin: %d\t", a.ref_begin + 1);
	fprintf(stdout, "target_end: %d\t\n", a.ref_end + 1);
	if (a.query_begin + 1) fprintf(stdout, "\tquery_begin: %d\t", a.query_begin + 1);
	fprintf(stdout, "\tquery_end: %d\n\n", a.query_end + 1);
	cout<<"\n===================================Alignment====================================\n"<<endl;
	if (a.cigar.size()>0) {
		int32_t c = 0, left = 0, e = 0, qb = a.ref_begin, pb = a.query_begin;
		uint32_t i;
		int cigarLen=a.cigar.size();
		while (e < cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if ((letter == 'I')||(letter == 'S')) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", ref_seq[q]);
						++ q;
					}
					++ count;
					if (count == 60) goto step2;
				}
			}
step2:
			fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i){
					if (letter == '=') {
						fprintf(stdout, "|");
						++q;
						++p;
					} else {
						
						fprintf(stdout, "*");
						if ((letter == 'I')||(letter == 'S')) ++p;
						else ++q;
					}
					++ count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
step3:
			p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if ((letter == 'D')||(letter == 'S')) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", read_seq[p]);
						++p;
					}
					++ count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
end:
			fprintf(stdout, "    %d\n\n", p);
		}
	}
        cout<<"________________________________________________________________________________\n";
	cout<<"================================================================================"<<endl;
}

#endif 

//#if  IGH
//double SmithWaterman(const string& S,const string& T,double& matches,int& position_b, int& position_e){
//#else
//double SmithWaterman(const string& S,const string& T,double& matches){
//#endif  
    // Declares a default Aligner
//  StripedSmithWaterman::Aligner aligner(2,4,10,1);
  // Declares a default filter
//   StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
 //  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  
   //aligner.Align(T.c_str(), S.c_str(), S.size(), filter, &alignment);

  //compute matches
 //  int csize= alignment.cigar.size();
 //  matches=0;
 //  double tot=0;
 //  for (int i=0;i<csize;++i){

   //  if(cigar_int_to_op(alignment.cigar[i])=='='){
    //   matches+=(double)cigar_int_to_len(alignment.cigar[i]);
    //   tot+=(double)cigar_int_to_len(alignment.cigar[i]);
      // }
   // tot+=(double)cigar_int_to_len(alignment.cigar[i]);
   //  else
   //    if(cigar_int_to_op(alignment.cigar[i])=='X')
// 	tot+=(double)cigar_int_to_len(alignment.cigar[i]);
	
  // }

  
 //#if BLAST_OUTPUT
  // if (matches>input_matches)  
  //   ssw_write(alignment,S,T);
 //#endif 
// matches=(matches/tot)*100;
//#if IGH 
  //  position_b=alignment.ref_begin;
  // position_e=alignment.ref_end;
//#endif  
//   return alignment.sw_score;
 //}
//#endif 

double LevenshteinDistance(const string& S,const string& T){
  
int sizeS=S.size();
int sizeT=T.size();

d[0][0]=0.0;
//double maxW=0;

for (int i=1;i<=sizeS;i++){
	d[0][i]=i;
}

for (int j=0;j<sizeT;j++){
	d[1][0]=j+1;
	for (int i=0;i<sizeS;i++){
	  if (S[i] == T[j]){
	    d[1][i+1] = d[0][i];
	  }
	  else{
	    d[1] [i+1] =min3(d[1][i]+INS,d[0][i+1]+DEL,d[0][i]+SUB);
	  }
	}
	for(int i=0; i<sizeS+1;i++){
	  d[0][i]=d[1][i];
	}
}
for(int i=0; i<sizeS+1;i++){
  d[0][i]=d[1][i];
}

return d[0][sizeS];
//maxW=d[0][sizeS];
//return maxW<THEDIT?true:false;
}

size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2)
{
  const size_t m(s1.size());
  const size_t n(s2.size());
  static int y=0;
  if( m==0 ) return n;
  if( n==0 ) return m;
 //  cout<<m<<" "<<n<<endl;
  size_t costs[10000];
  ++y;
 if (y%100000==1)
    cout<<"@@@:"<<y<<endl;
 
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
  cout<<"DIFF:"<<abs((int)m-(int)n)<<endl;
  cout<<"ED:"<<(costs[n]-abs((int)m-(int)n))<<endl;
  size_t result = costs[n];
  //delete [] costs;
 
  return result;
}
 
 inline void reverse(const string& S,string& r_S){
 int size=S.size();   
 for(int i=size-1;i>=0;--i){
      switch(S[i]){
	case 'A':
	case 'a':  
	 r_S.push_back('T');
	break;
	case 'T':
	case 't':  
	  r_S.push_back('A');
	break;
	case 'C':
	case 'c': 
	  r_S.push_back('G');
	break;     
	case 'G':
	case 'g':  
	  r_S.push_back('C');
	break;
	case 'N':
	case 'n':  
	  r_S.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }

}
#if  IGH
double SmithWaterman(const string& S,const string& T, const string& idS, const string& idT,
                     double& matches,int& positionBclone, int& positionEclone){
#else    
double SmithWaterman(const string& S,const string& T, const string& idS, const string& idT,
                     double& matches,int& positionBclone, int& positionEclone,int& positionBigh, int& positionEigh,string& cigar){
#endif   
    // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner(2,4,10,1);
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  
  aligner.Align(T.c_str(), S.c_str(), S.size(), filter, &alignment);
  
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
  #if BLAST_OUTPUT
   if (matches>input_matches)  
     ssw_write(alignment,S,T,idS,idT);
 #endif 
  matches=(matches/tot)*100;
  //matches=(matches/S.size())*100;
#if IGH 
   positionBclone=alignment.ref_begin;
   positionEclone=alignment.ref_end;
#else    
  positionBclone=alignment.ref_begin;
  positionEclone=alignment.ref_end;
  positionBigh=alignment.query_begin;
  positionEigh=alignment.query_end;
  cigar=alignment.cigar_string;
#endif    
  
  return alignment.sw_score;
}  

 
    
int main(int argc, char **argv)
{
  cout<<"\n\n ================================================================================\n";
  cout<<"|	      	        	   ALIGNMENT TOOL            		         |\n";
  cout<<"|	      	      		                                                 |\n";
  cout<<" ================================================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

#if IGH
 if (argc<2)
  {
    std::cerr<<"\n\nUSE:  AlignmentIGH <input>\n\n"<<endl;
#else  
 if (argc<4)
  {
#if  BLAST_OUTPUT    //aggiungere interactive
    std::cerr<<"\n\nUSE: AlignmentBlast <input> n <file_1> <file_2> .... <file_n> [<number_matches>]\n\n"<<endl;
    std::cerr<<"\t <number_matches> is a threshold (in terms of matches) to decide when the blast output is reported."<<endl;
#else
     std::cerr<<"\n\nUSE: Alignment <input> n <file_1> <file_2> .... <file_n> \n\n"<<endl;
#endif
     
#endif
    
    
    exit(EXIT_FAILURE);
  }
 /*#if !IGH //interactive 
if (string(argv[1])=="-I"){
  string in1(argv[2]);
  string in2(argv[3]);	
  string in2Rev;
  reverse(in2,in2Rev);
  double matches=0,matchesR=0, score=0, scoreR=0;
  idRef=new string("REF1");
  idRead=new string("REF1R");			
  score=SmithWaterman(in1,in2,matches);
  scoreR=SmithWaterman(in1,in2Rev,matchesR);
  cout<<score<<" "<<scoreR<<" "<<matches<<" "<<matchesR<<endl;
  return 1;
}  
#endif */
  
#if IGH
 std::string germfile[3] { "./Germline/IGHV.fa", "./Germline/IGHJ.fa",  "./Germline/IGHD.fa"};
 int num_file=3;
#else   
 int num_file=atoi(argv[2]);
 if (argc>num_file+3){
    input_matches=atoi(argv[num_file+3]);
    cout<<"\nOnly the alignment with matches >= "<< input_matches<<" are reported."<<endl;
 }
#endif 
  clock_t startGlobal,endGlobal;
  startGlobal=clock();

  set< pair<string,string> > ref[num_file];
 
  startGlobal=clock();

  

  for (int i=0;i<num_file;++i){
 #if IGH
    cout<<"\nReading  file: "<< germfile[i]<<endl; 
    ifstream f_ref(germfile[i],ifstream::in);
 #else    
    cout<<"\nReading  file: "<< argv[i+3]<<endl; 
    ifstream f_ref(argv[i+3],ifstream::in);
 #endif    
    if(!f_ref) {
      cerr << "\n*****Error opening input file "<< argv[i+3] <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    std::string id("");
    std::string seq("");
#if IGH	   
    char delim[]="|";
    class general::Parser parser;
#endif 
    bool pending=false;
    while (!f_ref.eof()){ 
	std::string buffer("");
        getline(f_ref,buffer);
	//std::cout<<buffer<<endl;
	if ((buffer[0]=='>')||(buffer[0]=='@')||(buffer=="")){
	  if (seq!=""){
	   ref[i].insert(make_pair(id,seq)); 
	   seq="";
	   pending=false;
       
	  }
	  if (buffer!=""){
#if IGH	
	  parser.update(delim,buffer);
	  if (parser.size()<2){
	    cerr << "\n*****Error parsing ID: "<<buffer<< "*****" << endl;
	    exit(EXIT_FAILURE);
	  }
  
	  id=parser.get(1)+" "+parser.get(3);	
#else
	  id=buffer;
#endif	  
	   }
	  }
	else{
	  seq+=buffer;
	  pending=true;
	  }
    }
  
    if (pending){
       ref[i].insert(make_pair(id,seq)); 
    }
    f_ref.close();
  }
  cout<<"\nReading clone file: "<< argv[1]<<endl; 
  ifstream f_clone(argv[1],ifstream::in);
  if(!f_clone) {
      cerr << "\n*****Error opening input file "<< argv[1] <<" *****" << endl;
      exit(EXIT_FAILURE);
  } 
  

  cout<<"\nStart Mapping: "<<endl; 
  while (!f_clone.eof()){//for all reads
   std::string id("");
   getline(f_clone,id);
   std::string clone("");
   getline(f_clone,clone);
   if ((id!="")&&(clone!="")){
   std::string cloneR(""); 
   reverse(clone,cloneR); 
   //scores
   vector<double> Maxscore(num_file,-999999),  MaxscoreR(num_file,-999999);
   //matches
   vector<double> max_matches(num_file,0) , max_matchesR(num_file,0);
   //positions
   vector<int> max_position(num_file,0), max_positionR(num_file,0);
   
#if !IGH
   // for alignement information 
    align_info MaxAlign[3];
    align_info MaxAlignR[3];
    std::string cigar;
    std::string Rcigar;
#endif
    
   set< pair<string,string> >::iterator Maxit[num_file];
   set< pair<string,string> >::iterator MaxitR[num_file];
   vector <string> similar_aligm(num_file,""), similar_aligmR(num_file,""); //to store other alignments with score equal to the max score 
   for (int i=0;i<num_file;i++){//for all files

//#if IGH    
     Maxit[i]=MaxitR[i]=ref[i].end();
     similar_aligm[i].clear();
     similar_aligmR[i].clear(); //to store other alignments with score equal to the max score 
//#endif  
     
     for (auto it=ref[i].begin();it!=ref[i].end();++it){   
        double matches=0,matchesR=0, score=0, scoreR=0;
        int positionB=0, positionBR=0, positionE=0, positionER=0;
        
#if IGH
	if (clone.size()>MINSIZEALIGN)
	   score=SmithWaterman(clone,it->second,string(),string(),matches,positionB, positionE);
	if (cloneR.size()>MINSIZEALIGN)
	   scoreR=SmithWaterman(cloneR,it->second,string(),string(),matchesR,positionBR, positionER);
#else
    
    int positionBigh=0, positionBRigh=0, positionEigh=0, positionERigh=0;

	//const string& idRef = id;
	//const string& idRead = it->second;
    std::string seqR;
    reverse(it->second,seqR); //PERCHÈ NON GLI PIACE
   
    score=SmithWaterman(clone,it->second,id,it->first,matches,positionB, positionE,positionBigh, positionEigh,cigar);  
	scoreR=SmithWaterman(clone,seqR,id,it->first,matchesR,positionBR, positionER,positionBRigh, positionERigh,Rcigar);
    std::cout << cigar << "\n";
    std::cout<< Rcigar << "\n";
     ofstream out(argv[num_file+3],ofstream::app);
      if(!out) 
      {  
      cerr << "\n*****Error opening output file "<< argv[num_file+3] <<" *****" << endl;
      exit(EXIT_FAILURE);
      }
      
      out << cigar << " " << score << " " << positionB << " " << positionE 
          << " " << positionBigh << " " << positionEigh << "    ";
      out << Rcigar << " " << scoreR << " " << positionBR << " " << positionER 
          << " " << positionBRigh << " " << positionERigh << "\n";
      
      out.close();
      
#endif 	
        
	if (score>Maxscore[i]){
	  Maxscore[i]=score;//*100/(2*(min(clone.size(),it->second.size())));
	  Maxit[i]=it;
	  max_matches[i]=matches;
#if !IGH      
      //For alignement information
      MaxAlign[i].cigar=cigar;
      MaxAlign[i].positionBclone=i!=0?positionB+MaxAlign[0].positionEclone:positionB;
      MaxAlign[i].positionEclone=i!=0?positionE+MaxAlign[0].positionEclone:positionE;
      MaxAlign[i].positionBigh=positionBigh;
      MaxAlign[i].positionEigh= positionEigh;
      
     
#endif
      
#if IGH	  
	  if (i==0)//case IGHV
	    max_position[i]=positionE;
	  else //IGHJ
	    max_position[i]=positionB; 
      similar_aligm[i].clear();
#endif 	  
	}
#if IGH
    else
    {
    if (score==Maxscore[i])
        {
         similar_aligm[i]+=string(string(" ")+it->first);
        }     
    }
#endif
	if (scoreR>MaxscoreR[i]){
	  MaxscoreR[i]=scoreR;//*100/(2*(min(clone.size(),it->second.size())));
	  MaxitR[i]=it;
	  max_matchesR[i]=matchesR;
#if IGH	  
	  if (i==0)//case IGHV
	    max_positionR[i]=positionER;
	  else //IGHJ
	    max_positionR[i]=positionBR; 
     similar_aligmR[i].clear(); 
#endif 	  
	}
#if IGH	
    else
    {
     if (scoreR==MaxscoreR[i])
     {
         similar_aligmR[i]+=string(" ")+it->first;
     }
         
    }
#endif 	
     }
#if IGH   
     if (i==0){//next IGHJ
       clone=clone.erase(0, max_position[i]);
       cloneR=cloneR.erase(0,max_positionR[i]);
     }
     else{//next IGHJD
       clone=clone.substr(0,max_position[i]);
       cloneR=cloneR.substr(0,max_positionR[i]); 
     }
#endif
     
   }
 cout<<"\n\n ===============================================================================";  
 cout<<"\n|				BEST ALIGNMENT	                                |\n"; 
 cout<<" ===============================================================================\n\n";  
      cout<<"Ref:  "<<id<<endl<<endl;
      double sum=0.0,sumR=0.0;
      for (int i=0;i<num_file;i++){
	sum+=Maxscore[i];
	sumR+=MaxscoreR[i];
	//cout<<Maxit[i]->first<<endl<<"\t"<<Maxscore[i]<<"\t"<< max_matches[i]<<endl;
      }
      if (sum>sumR){
	for (int i=0;i<num_file;i++){ 
#if IGH
       
#else        
	  cout<<"\t======================================\n\t File: "<< argv[i+3]<<endl;
#endif      
	  if (Maxit[i]==ref[i].end())
	    cout<<"----"<<endl;
	  else  
	  {
#if IGH
          cout<<"\t Read: "<<Maxit[i]->first;
          (similar_aligm[i]=="")?cout<<endl:cout<<" ("<<similar_aligm[i]<<")"<<endl;
#else
           cout<<"\t Read: "<<Maxit[i]->first<<endl;
#endif    
	  }
	  cout<<"\t Score SW:"<<Maxscore[i]<<"\t Matches:"<< max_matches[i]<<endl;
	} 
	cout<<"\t======================================\n";
	cout<<"Tot. score:"<<sum<<endl;
      }
      else{
	cout<<"\t(reverse complement)"<<endl;
	for (int i=0;i<num_file;i++){
#if IGH
        cout<<"\t======================================\n\t File: "<< germfile[i]<<endl; 
#else
        cout<<"\t======================================\n\t File: "<< argv[i+3]<<endl;
#endif        
	 if (MaxitR[i]==ref[i].end()) 
	   cout<<"----"<<endl;
	 else 
 #if IGH
       cout<<"\t Read: "<<MaxitR[i]->first;
       (similar_aligmR[i]=="")?cout<<endl:cout<<" ("<<similar_aligmR[i]<<")"<<endl;
#else         
	   cout<<"\t Read: "<<MaxitR[i]->first<<endl;
#endif         
	   cout<<"\t Score SW:"<<MaxscoreR[i]<<"\t Matches:"<< max_matchesR[i]<<endl;
	}
	 cout<<"\t======================================\n";
	 cout<<"Tot. score:"<<sumR<<endl;
      }
   }
 }
 f_clone.close();
cout<<"\n================================================================================\n\n";

 endGlobal=clock();
         
 cout<<"\n=====================================TIME=======================================\n\n";
 cout<<"\tTime to identify alignment: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n================================================================================\n\n";

}
