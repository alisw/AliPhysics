// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTTPCDEDXDATA_H
#define ALIHLTTPCDEDXDATA_H

struct AliHLTTPCdEdxInfo
{
  float fdEdxTotIROC; //Do not change the order, elements are accessed as array float[10]
  float fdEdxTotOROC1;
  float fdEdxTotOROC2;
  float fdEdxTotOROC;
  float fdEdxTotTPC;
  float fdEdxMaxIROC;
  float fdEdxMaxOROC1;
  float fdEdxMaxOROC2;
  float fdEdxMaxOROC;
  float fdEdxMaxTPC;
  short nHitsIROC;
  short nHitsOROC1;
  short nHitsOROC2;
  short nHitsSubThresholdIROC;
  short nHitsSubThresholdOROC1;
  short nHitsSubThresholdOROC2;
};

struct AliHLTTPCdEdxData
{
  int fCount; // Number of tracks we have dedx info for
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCdEdxInfo fdEdxInfo[1]; // array of dedx info 
#else
  AliHLTTPCdEdxInfo fdEdxInfo[0]; // array of dedx info
#endif

};

#endif
