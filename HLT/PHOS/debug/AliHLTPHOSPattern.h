//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSPATTERN_H
#define ALIHLTPHOSPATTERN_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 


//#include "AliHLTPHOSConstants.h"
//using  namespace PhosHLTConst;
#include "AliHLTPHOSBase.h"
class AliHLTPHOSUtilities;


class AliHLTPHOSPattern : public AliHLTPHOSBase
{
 public:
  AliHLTPHOSPattern(const int *pattern, const int length);
  ~AliHLTPHOSPattern();
//   const int AddPattern(const int *readbackpattern,  const int nsamples =  ALTRO_MAX_SAMPLES, const int nPresamples = 0);
//   const bool CheckDoExistPattern(const int *readbackpattern,  const int nsamples =  ALTRO_MAX_SAMPLES, const int nPresamples = 0);
//   const AliHLTPHOSPattern *GetNextPtr() const {return fPattern;}; 
//   const int  GetPattern(int *pattern,  const int maxlengths =  ALTRO_MAX_SAMPLES) const;
//   const int GetPatternLength() const {return  fPatternLength;};
//   const int ValidatePattern(const int *readbackpattern,  const int nsamples =  ALTRO_MAX_SAMPLES, const int nPresamples = 0) const;
  int AddPattern(const int *readbackpattern,  const int nsamples =  ALTROMAXSAMPLES, const int nPresamples = 0);
  bool CheckDoExistPattern(const int *readbackpattern,  const int nsamples =  ALTROMAXSAMPLES, const int nPresamples = 0);
  AliHLTPHOSPattern *GetNextPtr() const {return fPattern;}; 
  int  GetPattern(int *pattern,  const int maxlengths =  ALTROMAXSAMPLES) const;
  int GetPatternLength() const {return  fPatternLength;};
  int ValidatePattern(const int *readbackpattern,  const int nsamples =  ALTROMAXSAMPLES, const int nPresamples = 0) const;
  void PrintPattern(const int nPerLine = ALTROMAXSAMPLES);
 
 
 private:
  AliHLTPHOSPattern(const AliHLTPHOSPattern &);
  AliHLTPHOSPattern & operator = (const AliHLTPHOSPattern &); 
  AliHLTPHOSPattern();

  // const int DoComparePattern(const int *inputPattern1, const int *inputPattern2, const int nsamples) const;
//   const bool CheckPatternLength(const int nsamples, const int nPresamples) const;
  int DoComparePattern(const int *inputPattern1, const int *inputPattern2, const int nsamples) const;
  bool CheckPatternLength(const int nsamples, const int nPresamples) const;
  void SetPattern(const int *pattern, const int length);
  int fVal[ALTROMAXSAMPLES];
  int fPatternLength;
  int fN; /**<The number of detected patterns equal to fPattern*/ 
  int fCnt; 
  AliHLTPHOSPattern *fPattern;
   
  AliHLTPHOSUtilities *fUtilitiesPtr;

};

#endif

