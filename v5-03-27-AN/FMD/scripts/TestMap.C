//____________________________________________________________________
//
// $Id$
//
// Test I/O of ALiFMDMap
//
/** @file    TestMap.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sat Dec 16 01:29:03 2006
    @brief   Test of uniquenss of map
    @ingroup FMD_script     
*/
#include "STEER/AliFMDFloatMap.h"
#include <TArrayI.h>
#include <iostream>
#include <iomanip>
#include <map>

struct Check 
{
  UShort_t d;
  Char_t   r;
  UShort_t s;
  UShort_t t;
  Float_t  v;
};

std::ostream& operator<<(std::ostream& o, const Check& c) 
{
  UShort_t v = UShort_t(c.v);
  o << "FMD" << c.d << c.r << '[' 
    << std::setw(2) << c.s << ',' 
    << std::setw(3) << c.t << "] (" 
    << std::setw(6) << v << ')';
   return o;
}

Check Name(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v) 
{
  Check c;
  c.d = d;
  c.r = r;
  c.s = s;
  c.t = t;
  c.v = v;
  return c;
}


void
TestMap() 
{
  AliFMDFloatMap m;
  // m.SetBit(AliFMDMap::kNeedUShort);
  typedef std::map<Int_t, Check> CheckMap;
  CheckMap c;
  
  for (UShort_t d = 1; d <= 3; d++) {
    Char_t rings[] = { 'I', 'O', '\0' };
    for (Char_t* rp = rings; *rp; rp++) {
      Char_t r = *rp;
      std::cout << "FMD" << d << r << " " << std::flush;
      for (UShort_t s = 0; s < 40; s++) {
	std::cout << "." << std::flush;
	for (UShort_t t = 0; t < 512; t++) {
	  Int_t i = m.CheckIndex(d, r, s, t);
	  CheckMap::iterator z = c.find(i);
	  Float_t v = (t + 512 * (s + 40 * ((r=='I'? 0 : 1) + 2 * d)));
	  if (z != c.end()) 
	    std::cout << '\n' << Name(d,r,s,t,v)  << " but aleady seen " 
		      << z->second  << std::flush;
	  else {
	    m(d,r,s,t) = v;
	    c[i] = Name(d,r,s,t,v);
	  } // else
	} // for t
      } // for s 
      std::cout << "done" << std::endl;
    }
  }
}

    
	
