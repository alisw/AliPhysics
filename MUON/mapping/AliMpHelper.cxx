/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpHelper.cxx,v 1.3 2005/09/19 19:01:31 ivana Exp $

#include "AliMpHelper.h"

#include "AliLog.h"
#include "TArrayI.h"
#include "TClass.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

ClassImp(AliMpHelper)

//_____________________________________________________________________________
AliMpHelper::AliMpHelper() : TObject()
{
  //
  // Default (empty) ctor.
  // 
} 
 
//_____________________________________________________________________________
AliMpHelper::~AliMpHelper()
{
  //
  // Dtor.
  //
}

//_____________________________________________________________________________
void AliMpHelper::DecodeName(const char* name, char sep, TArrayI& theList)
{
  //
  // From a string of the form "i-j;k;l;m-n" returns an integer array
  // containing all the integers from i to j, then k, l and then from m to
  // n.
  //
  theList.Set(0);
  
  TString str(name);
  
  if ( str.Length() == 0 )
  {
    // protection against empty input string.
    return;
  }
  
  // Get substrings separated by 'sep'
  TObjArray* ranges = str.Tokenize(sep);
  
  // Finally takes each substring (which ought to be a range of the form
  // x-y), and decode it into the theList integer vector.
  for ( Int_t i = 0; i < ranges->GetEntriesFast(); ++i )
  {
    int m1;
    int m2;
    int n;
    int incr;
    TString& s = ((TObjString*)ranges->At(i))->String();
    GetRange(s.Data(),m1,m2,incr,n);
    int m = m1;
    while ( n > 0 )
    {
      theList.Set(theList.GetSize()+1);
      theList[theList.GetSize()-1] = m;
      m += incr;
      --n;
    }
  }
  
  delete ranges;
}

//_____________________________________________________________________________
void 
AliMpHelper::GetRange(const char* cstr, Int_t& begin, Int_t& end, 
                      Int_t& incr, Int_t& n)
{
  //
  // From a string of the form "m-n" returns a range (begin,end),
  // its ordering (incr=+-1) and its size (abs(begin-end)+1)
  //
  TString str(cstr);
  
  incr = 1;
  Ssiz_t pos = str.First('-');
  if ( pos < 0 )
  {
    begin = str.Atoi();
    end = -1;
    n = 1;
  }
  else
  {
    begin = str.Atoi();
    end = TString(str(pos+1,str.Length()-pos)).Atoi();
    if ( begin > end )
    {
      incr = -1;
      n = begin-end+1;
    }
    else
    {
      n = end-begin+1;
    }    
  }
}

//_____________________________________________________________________________
TString AliMpHelper::Normalize(const char* line)
{
  //
  // Remove multiple blanks, and blanks in the begining/end.
  //
  TString rv(line);
  
  if ( rv.Length() <= 0 ) return TString();
  
  while ( rv[0] == ' ' )
  {
    rv.Remove(0,1);
  }
  while ( rv[rv.Length()-1] == ' ' )
  {
    rv.Remove(rv.Length()-1,1);
  }
  Ssiz_t i(0);
  bool kill = false;
  for ( i = 0; i < rv.Length(); ++i )
  {
    if ( rv[i] == ' ' )
    {
      if (kill)
	    {
	      rv.Remove(i,1);
	      --i;
	    }
      else
	    {
	      kill = true;
	    }
    }
    else
    {
      kill = false;
    }
  }
  return rv;
}
