#ifndef ALI_MUON_ST1_DECODER_H
#define ALI_MUON_ST1_DECODER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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


#include <vector>
#include <utility>
#include <string>

#include "AliMUONSt1Types.h"

namespace decoder
{
  vector<string> SplitNtuples(const string& s,
                              const string& leftSep ="({[\"'/",
                              const string& rightSep=")}]\"'/");
  vector<string> SplitList(const string& s,const string& sep=";,");
  vector<int> DecodeListRanges(const string& s,const string& sep=";,",const string& rangeSep="/");
  vector< pair<int,int> > DecodeListOfIntRanges(const string& s,const string& sep=";,",const string& rangeSep="/");
  vector< pair<double,double> > DecodeListOfFloatRanges(const string& s,const string& sep=";,",const string& rangeSep="/");
}

#endif //ALI_MUON_ST1_DECODER_H 
