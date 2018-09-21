#ifndef ALIROOT_ALIPARSER_H
#define ALIROOT_ALIPARSER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/// \ingroup STAT
/// \class AliParser
/*!
* \brief Class provides flexibility static methods for tokenizing strings.
*  You can find explanation, details and examples in description of each function.
* \author  <a href="mailto:marian.ivanov@cern.ch">Marian Ivanov</a>, <a href="mailto:boris.rumyantsev@cern.ch">Boris Rumyantsev</a>
*/

#include "TObject.h"
#include "TString.h"
#include "TMatrixD.h"
#include <map>
#include <vector>

class AliParser {
  public:
    static std::vector<TString> ExtractBetween(const char *inputString, const char *startStr, const char *endStr, Int_t verbose=0);
    static std::vector<TString> Split(const char *inputString, const char delimiter=',', Int_t verbose=0);
    static std::vector<TString> ExtractSurroundingBy(const char *inputString, const char begin='(', const char end = ')', Int_t verbose=0);
    static std::map<TString, TString> Parse(const char *inputString, Int_t verbose=0, std::vector<TString> defKeys = std::vector<TString>{});
    static  std::vector<Int_t> Slice2IArray(const char *inputString);
    static TMatrixD Slice2Matrix(const char *inputString, Int_t verbose=0);
  //private:
    static void Map2Array (std::vector<Double_t> &array, std::map<Int_t, std::vector<Double_t> > iRanges, Int_t *indexes, Int_t cnt);
    static void FillFloatFlagArray(std::vector<Double_t> &, TString);
  ClassDef(AliParser,1);
};
#endif
