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

//  @file   AliHLTTPCHWCFDivisionUnit.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @date   
//  @brief  Division unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note 

#include "AliHLTTPCHWCFDivisionUnit.h"
#include <iostream>
#include <algorithm>


AliHLTTPCHWCFDivisionUnit::AliHLTTPCHWCFDivisionUnit()
  : 
  fSinglePadSuppression(1), fClusterLowerLimit(0), fkInput(0),fOutput()
{
  //constructor 
}


AliHLTTPCHWCFDivisionUnit::~AliHLTTPCHWCFDivisionUnit()
{   
  //destructor 
}

AliHLTTPCHWCFDivisionUnit::AliHLTTPCHWCFDivisionUnit(const AliHLTTPCHWCFDivisionUnit&)
  : 
  fSinglePadSuppression(1),fClusterLowerLimit(0),fkInput(0),fOutput()
{
}


AliHLTTPCHWCFDivisionUnit& AliHLTTPCHWCFDivisionUnit::operator=(const AliHLTTPCHWCFDivisionUnit&)
{
  // dummy  
  return *this;
}

int AliHLTTPCHWCFDivisionUnit::Init()
{
  // initialise
  fkInput = 0;
  return 0;
}
  

int AliHLTTPCHWCFDivisionUnit::InputStream( const AliHLTTPCHWCFClusterFragment *fragment )
{
  // input stream of data
  fkInput = fragment;
  if( fkInput ){
    //if( fInput->fFlag==1 && fInput->fQ>0 ){
    //std::cout<<"Division: input F "<<fInput->fFlag<<" R "<<fInput->fRow
    //       <<" C "<<(fInput->fQ>>12)
    //       <<" P "<<fInput->fPad
    //       <<" T "<<fInput->fTMean
    //       <<" CoGPad "<<((float)fInput->fP)/fInput->fQ
    //       <<" CoGTime "<<((float)fInput->fT)/fInput->fQ
    //       <<std::endl;
    //} else std::cout<<"Division: input F "<<fInput->fFlag<<" R "<<fInput->fRow<<" P "<<fInput->fPad
    //	    <<" Q "<<(fInput->fQ>>12)<<std::endl;
  }
  return 0;
}

const AliHLTTPCHWCFCluster *AliHLTTPCHWCFDivisionUnit::OutputStream()
{ 
  // output stream of data

  if( !fkInput ) return 0;

  if( fkInput->fFlag==2 ){ // RCU trailer word
    fOutput.fFlag = 2;
    fOutput.fRowQ = fkInput->fRow; // rcu word
    fkInput = 0;
    return &fOutput;;
  }

  if( fkInput->fFlag!=1 ) return 0;
  
  if( fkInput->fQ==0 ) return 0;
  if( fSinglePadSuppression && fkInput->fQ==fkInput->fLastQ && !fkInput->fBorder ) return 0;
  if( fkInput->fQ < fClusterLowerLimit ) return 0;

  AliHLTFloat32_t q = fkInput->fQ;
  
  fOutput.fFlag = 1;
  fOutput.fRowQ = (((AliHLTUInt32_t) 0x3)<<30) + ((fkInput->fRow &0x3f)<<24) + ((fkInput->fQ>>(AliHLTTPCHWCFDefinitions::kFixedPoint-6))&0xFFFFFF);		  
  *((AliHLTFloat32_t*)&fOutput.fP) = (float)fkInput->fP/q;
  *((AliHLTFloat32_t*)&fOutput.fT) = (float)fkInput->fT/q;
  *((AliHLTFloat32_t*)&fOutput.fP2) = (float)fkInput->fP2/q;
  *((AliHLTFloat32_t*)&fOutput.fT2) = (float)fkInput->fT2/q;
 
  // MC part

  AliHLTTPCClusterMCWeight emptyWeight = {-1,0};

  fOutput.fMC.fClusterID[0] = emptyWeight;
  fOutput.fMC.fClusterID[1] = emptyWeight;
  fOutput.fMC.fClusterID[2] = emptyWeight;
  
  vector<AliHLTTPCClusterMCWeight> labels = fkInput->fMC;
  sort(labels.begin(), labels.end(), CompareMCLabels);
  for( unsigned int i=1; i<labels.size(); i++ ){
    if(labels[i-1].fMCID==labels[i].fMCID ){
      labels[i-1].fWeight+=labels[i].fWeight;
      labels[i].fWeight = 0;
    }
  }

  sort(labels.begin(), labels.end(), CompareMCWeights );
    
  for( unsigned int i=0; i<3 && i<labels.size(); i++ ){
    if( labels[i].fMCID <0 ) continue;
    fOutput.fMC.fClusterID[i] = labels[i];
  }

  fkInput = 0;

  return &fOutput;
}

