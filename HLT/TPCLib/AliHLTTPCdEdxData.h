// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTTPCDEDXDATA_H
#define ALIHLTTPCDEDXDATA_H

struct AliHLTTPCdEdxData
{
  int fVersion; //Currently, only version 1 present
  int fValuesPerTrack; //Version 1: 10 values: qTot - IROC OROC OROCLONG OROCALL TPCALL ; qMax - IROC OROC OROCLONG OROCALL TPCALL
  int fCount; // Number of tracks for which we have information
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  float fdEdxInfo[1]; // array of clusters  
#else
  float fdEdxInfo[0]; // array of clusters 
#endif

};

#endif
