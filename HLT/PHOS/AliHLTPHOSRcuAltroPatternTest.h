//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUALTROPATTERNTEST_H
#define ALIHLTPHOSRCUALTROPATTERNTEST_H

// 1
// 2
// 3
// 4
// 5

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSBase.h"

//#include "AliHLTPHOSConstants.h"
//#include "AliHLTDataTypes.h"

//using  namespace PhosHLTConst;


class AliHLTPHOSPattern;

class AliHLTPHOSRcuAltroPatternTest : public AliHLTPHOSBase
//class AliHLTPHOSRcuAltroPatternTest
{
 public:
  AliHLTPHOSRcuAltroPatternTest(const AliHLTUInt8_t moduleID, const AliHLTUInt8_t rcuX, const AliHLTUInt8_t rcuZ, const int *pattern, const int length);
  virtual ~AliHLTPHOSRcuAltroPatternTest();
//   const int ValidateAltroPattern(const int *inputPattern,  const int length =  ALTRO_MAX_SAMPLES, const int presamples = 0) const;
//   const int AddPattern(const int *inputPattern,  const int z, const int x, const int gain, const int length =  ALTRO_MAX_SAMPLES, const int presamples = 0); 
//   const int countPatterns(const AliHLTPHOSPattern *pattern) const;
//   const int countAllPatterns(const int length, const bool printpatterns = true);
  int ValidateAltroPattern(const int *inputPattern,  const int length =  ALTROMAXSAMPLES, const int presamples = 0) const;
  int AddPattern(const int *inputPattern,  const int z, const int x, const int gain, const int length =  ALTROMAXSAMPLES, const int presamples = 0); 
  int countPatterns(const AliHLTPHOSPattern *pattern) const;
  int countAllPatterns(const int length, const bool printpatterns = true);
  void PrintStatistics() const;
 
 //  void PrintPatterns(const AliHLTPHOSPattern *pattern) const;
  
 private:
  void PrintPatterns(AliHLTPHOSPattern *pattern);
  AliHLTPHOSRcuAltroPatternTest();
  AliHLTPHOSRcuAltroPatternTest(const AliHLTPHOSRcuAltroPatternTest & );
  AliHLTPHOSRcuAltroPatternTest & operator = (const AliHLTPHOSRcuAltroPatternTest &);
 
  //  void PrintPattern() const;
  //  bool Compare(const AliHLTPHOSPattern);

  unsigned long  fNEqual[NZROWSRCU][NXCOLUMNSRCU][NGAINS];
  unsigned long  fNNotEqual[NZROWSRCU][NXCOLUMNSRCU][NGAINS];
  const AliHLTUInt8_t fModuleID;    /**<ID of the module this component read data from (0-4)*/
  const AliHLTUInt8_t fRcuX;        /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  const AliHLTUInt8_t fRcuZ;        /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTPHOSPattern *fReferenceAltroPattern; /**<The pattern stored in the altro*/
  AliHLTPHOSPattern *fPerChannelPatterns[NZROWSRCU][NXCOLUMNSRCU][NGAINS]; /**<Pattern actually read back from the electronics*/
  
  unsigned long fCnt; //REMOVE !!

};

#endif
