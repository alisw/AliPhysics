#ifndef ALI_MUON_ST1_DECODER_H
#define ALI_MUON_ST1_DECODER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

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

#include <string>
#include <utility>
#include <vector>

#ifndef __HP_aCC
  using std::string;
  using std::pair;
  using std::vector;
#endif  

#include <Rtypes.h>

typedef vector<string>  StringVector; 
typedef vector<Int_t>   IntVector; 
typedef pair<Int_t, Int_t> IntPair;
typedef vector<IntPair>    IntPairVector;
typedef pair<Double_t, Double_t>  DoublePair;
typedef vector<DoublePair>        DoublePairVector;

namespace decoder
{
  StringVector     SplitNtuples(const string& s,
                                const string& leftSep ="({[\"'/",
                                const string& rightSep=")}]\"'/");
  StringVector     SplitList(const string& s, const string& sep=";,");

  IntVector        DecodeListRanges(const string& s, const string& sep=";,",
                                    const string& rangeSep="/");
  IntPairVector    DecodeListOfIntRanges(const string& s, const string& sep=";,",
                                    const string& rangeSep="/");
  DoublePairVector DecodeListOfFloatRanges(const string& s, const string& sep=";,",
                                    const string& rangeSep="/");
}

#endif //ALI_MUON_ST1_DECODER_H 
