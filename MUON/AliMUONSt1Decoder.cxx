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
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1Decoder
// -----------------------
// A generic set of functions (defined in the <decoder> namespace).
// Used to decode formatted strings, eg. a list of integer ranges,
// or a list of sub-strings separated by delimiters such as '(','{','[', ... .
// Example: 
//   (string 1) (string 2) [string3] {string4} [ (string5.1) (string5.2) ] 
//   note :                                      |_____________________|
//                                                         |
//                                             this is just ONE substring.
// Included in AliRoot 2003/01/28

#include "AliMUONSt1Decoder.h"

//_____________________________________________________________________
StringVector decoder::SplitNtuples(const string& s,
                                   const string& leftSep,
                                   const string& rightSep)

{
// separate substrings in <s>, by using the first level of delimiters
// given in the argument
// Example: 
// (string 1) (string 2) [string3] {string4} [ (string5.1) (string5.2) ]
// returns a list of 5 substrings
// "string 1"
// "string 2"
// "string 3"
// "string 4" and
// " (string5.1) (string5.2) "
// --
  StringVector ans;
  string::size_type idx = 0;
  do {
    idx = s.find_first_of(leftSep,idx);
    if (idx != string::npos) {
      string::size_type sepNum = leftSep.find_first_of(s[idx],0);
      if (sepNum>=rightSep.length()){
      	idx++;
	continue;
      }
      Int_t count=1;
      string::size_type idx2 = idx+1;
      while ((count>0) && (idx2<s.length())){
      	if (s[idx2] == leftSep[sepNum])  count++;
      	if (s[idx2] == rightSep[sepNum]) count--;
	idx2++;
      }
      if (count != 0) return ans; // bad format...stop here
      ans.push_back(s.substr(idx+1,idx2-idx-2));
      idx=idx2;
    }
  } while (idx != string::npos);
  return ans;
}

//_____________________________________________________________________
StringVector decoder::SplitList(const string& s,const string& sep)
{
// split <s> into several substrings, by using any of the delimters in <sep> 
// Example : str1 ; str2 , str3
// gives a list of 3 substrings ("str1", "str2" and "str3") if
// the delimiter parameter is ",;"
// and a list of 2 substrings ("str1","str2,str3") if
// the delimiter parameter is ";"
// --
  StringVector ans;
  string::size_type i=0,j;
  while ((j=s.find_first_of(sep,i))!=string::npos){
    ans.push_back(s.substr(i,j-i));
    i=j+1;
  }
  ans.push_back(s.substr(i,s.length()-i));
  return ans;
}

//_____________________________________________________________________
IntVector 
decoder::DecodeListRanges(const string& s, const string& sep,
                          const string& rangeSep)
{
// decode <s> as a list of integers
// Example: 192/199 ; -10/-7
// gives a list of 12 integers:
// 192, 193, 194, 195, 196, 197, 198, 199, -10, -9, -8, -7
// --
  IntVector ans;
  IntPairVector rangeList = DecodeListOfIntRanges(s,sep,rangeSep);
  for (UInt_t i = 0 ;i<rangeList.size();i++){
    for (Int_t j=rangeList[i].first;j<=rangeList[i].second;j++){
      ans.push_back(j);
    }
  }
  return ans;
}

//_____________________________________________________________________
IntPairVector
decoder::DecodeListOfIntRanges(const string& s, const string& sep,
                               const string& rangeSep)
{
// decodes <s> as a list of int ranges
// Example: 192/303 ; -10/-1
// gives a list of two int pairs:
// pair(192,303) and another pair (-10,-1)
// --
  string::size_type i=0;
  IntPairVector ans;
  if (s.empty()) return ans;

  StringVector parts = decoder::SplitList(s,sep);

  for (UInt_t k=0;k<parts.size();k++){
    i=parts[k].find_first_of(rangeSep);
    Int_t from,to;
    if (i == string::npos) {
      from = atoi(parts[k].c_str());
      to = from;
    } else {
      from=atoi(parts[k].substr(0,i).c_str());
      to  =atoi(parts[k].substr(i+1,parts[k].length()-i-1).c_str());
    }
    ans.push_back(IntPair(from,to));
  }
  return ans;
}

//_____________________________________________________________________
DoublePairVector 
decoder::DecodeListOfFloatRanges(const string& s, const string& sep,
                                 const string& rangeSep)
{
// decodes <s> as a list of double (floating point) ranges
// Example : 192.33/303.26 ; -10.2/-1.41
// gives a list of two double precision floating point (phew!) pairs:
// pair (192.33,303.26) and another pair (-10.2,-1.41)
// --
  string::size_type i=0;
  DoublePairVector ans;
  if (s.empty()) return ans;

  StringVector parts = decoder::SplitList(s,sep);

  for (UInt_t k=0;k<parts.size();k++){
    i=parts[k].find_first_of(rangeSep);
    Double_t from,to;
    if (i == string::npos) {
      from = atof(parts[k].c_str());
      to = from;
    } else {
      from=atof(parts[k].substr(0,i).c_str());
      to  =atof(parts[k].substr(i+1,parts[k].length()-i-1).c_str());
    }
    ans.push_back(DoublePair(from,to));
  }
  return ans;
}
