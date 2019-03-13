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


#ifndef __CON_H__
	#define __CON_H__
	#include "conf.h"
	#include "const.h"
	#include "ssw_cpp.h"
	#include "ssw.h"
	
#endif

#ifndef __VCT_H__
	#define __VCT_H__
	#include <vector>
#endif

#ifndef __SET_H__
	#define __SET_H__
	#include <set>
#endif


class align_info{
public:    
    int positionBclone {0};
    int positionEclone {0};
	int positionBigh   {0}; 
    int positionEigh   {0};
    std::string cigar    {""};
    bool reverse_align {false};
};


class Info{
  public:
    
  //! nucleotide clone string   
  std::string clone {""};
  //! clone frequencies in any the sample
  std::vector <int> freq;
  //! genomic position IHGV,IHGJ,IHGD
  std::vector <std::string> info_identity;
  //! genomic identity 
  std::vector <double> info_pos;
  

  std::vector <align_info> align; 

  //! spike in presence
  std::string spike_in_name {""}; 
  bool spike_in {false};
  
  //! Max score for V,D,J alignemnt
  std::vector <int> tmp;  
  
  //!Constructor. 
 Info(const std::string& clone,  std::vector <int>& freq,  std::vector <std::string>&  info_identity,  std::vector <double>& info_pos, bool & spike_in, std::string& spike_in_name){
   this->clone=clone;
   this->freq=freq;
   this->info_pos=info_pos;
   this->info_identity=info_identity;
   this->spike_in=spike_in;
   this->spike_in_name=spike_in_name;
    
 }
 //!Constructor.
 Info(const std::string& clone,  std::vector <int>& freq){
   this->clone=clone;
   this->freq=freq;
 }
};
