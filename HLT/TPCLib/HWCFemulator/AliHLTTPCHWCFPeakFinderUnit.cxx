// $Id: AliHLTTPCHWCFPeakFinderUnit.cxx 51236 2011-08-22 16:01:48Z sgorbuno $
//****************************************************************************
//* This file is property of and copyright by the ALICE HLT Project          * 
//* ALICE Experiment at CERN, All rights reserved.                           *
//*                                                                          *
//* Primary Authors: Sergey Gorbunov, Torsten Alt                            *
//* Developers:      Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de> *
//*                  Torsten Alt <talt@cern.ch>                              *
//*                  for The ALICE HLT Project.                              *
//*                                                                          *
//* Permission to use, copy, modify and distribute this software and its     *
//* documentation strictly for non-commercial purposes is hereby granted     *
//* without fee, provided that the above copyright notice appears in all     *
//* copies and that both the copyright notice and this permission notice     *
//* appear in the supporting documentation. The authors make no claims       *
//* about the suitability of this software for any purpose. It is            *
//* provided "as is" without express or implied warranty.                    *
//****************************************************************************

//  @file   AliHLTTPCHWCFPeakFinderUnit.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @date   
//  @brief  Channel Processor unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note 

#include "AliHLTTPCHWCFPeakFinderUnit.h"
#include <iostream>
#include <cstdio>


AliHLTTPCHWCFPeakFinderUnit::AliHLTTPCHWCFPeakFinderUnit()
  :
  fOutput(),
  fkBunch(0),
  fChargeFluctuation(0),
  fDebug(0)
{
  //constructor 
  Init();
}


AliHLTTPCHWCFPeakFinderUnit::~AliHLTTPCHWCFPeakFinderUnit()
{   
  //destructor 
}

AliHLTTPCHWCFPeakFinderUnit::AliHLTTPCHWCFPeakFinderUnit(const AliHLTTPCHWCFPeakFinderUnit&)
  :
  fOutput(),
  fkBunch(0),
  fChargeFluctuation(0),
  fDebug(0)
{
  // dummy
  Init();
}

AliHLTTPCHWCFPeakFinderUnit& AliHLTTPCHWCFPeakFinderUnit::operator=(const AliHLTTPCHWCFPeakFinderUnit&)
{
  // dummy  
  return *this;
}

int AliHLTTPCHWCFPeakFinderUnit::Init()
{
  // initialise  

  fkBunch = 0;
  return 0;
}

int AliHLTTPCHWCFPeakFinderUnit::InputStream( const AliHLTTPCHWCFBunch *bunch )
{
  // input stream of data 

  if( bunch && fDebug ){
    printf("\nHWCF Peak Finder: input bunch F %1d R %3d P %3d  NS %2ld:\n",
	   bunch->fFlag, bunch->fRow, bunch->fPad, bunch->fData.size());
    for( unsigned int i=0; i<bunch->fData.size(); i++ ){
      const AliHLTTPCHWCFDigit &d = bunch->fData[i];
      printf("   q %2d t %3d ", d.fQ, d.fTime);
      printf("(");
      for( int j=0; j<3; j++ ) printf(" {%d,%2.0f}",d.fMC.fClusterID[j].fMCID, d.fMC.fClusterID[j].fWeight );
      printf(" )\n");      
    }
  }

  fkBunch = bunch;
  return 0;
}

const AliHLTTPCHWCFBunch *AliHLTTPCHWCFPeakFinderUnit::OutputStream()
{ 
  // output stream of data 

  if( !fkBunch ) return 0;

  fOutput.fFlag = fkBunch->fFlag;
  fOutput.fRow = fkBunch->fRow;
  fOutput.fPad = fkBunch->fPad;
  fOutput.fBranch = fkBunch->fBranch;
  fOutput.fBorder = fkBunch->fBorder;
  fOutput.fGain = fkBunch->fGain;
  fOutput.fData.clear();
  fOutput.fData.insert(fOutput.fData.end(),fkBunch->fData.begin(), fkBunch->fData.end());
  fkBunch = 0;

  if( fOutput.fFlag !=1 ){ // rcu trailer word, forward it 
    return &fOutput;
  }

  bool slope = 0;
  AliHLTUInt32_t qLast = 0;
  AliHLTUInt32_t n = fOutput.fData.size();
  
  for( AliHLTUInt32_t i=0; i<n; i++ ){
    AliHLTUInt32_t q = fOutput.fData[i].fQ;    
    if( !slope ){
      if(q + fChargeFluctuation < qLast   ){ // peak
	slope = 1;
	if( i>0 ) fOutput.fData[i-1].fPeak = 1;
      }
    }else if( q > qLast + fChargeFluctuation ){ // minimum
      slope = 0;
      if( i>0 ) fOutput.fData[i-1].fPeak = 2;
    }
    qLast = q;
  }
  
  if( n>0 ){
    if( !slope ) fOutput.fData[n-1].fPeak = 1;
    else fOutput.fData[n-1].fPeak = 2;
  }
  
  return &fOutput;
}
