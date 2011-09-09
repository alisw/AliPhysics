// $Id$
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
  fSinglePadSuppression(1), fClusterLowerLimit(0), fkInput(0),fOutput(), fDebug(0)
{
  //constructor 
}


AliHLTTPCHWCFDivisionUnit::~AliHLTTPCHWCFDivisionUnit()
{   
  //destructor 
}

AliHLTTPCHWCFDivisionUnit::AliHLTTPCHWCFDivisionUnit(const AliHLTTPCHWCFDivisionUnit&)
  : 
  fSinglePadSuppression(1),fClusterLowerLimit(0),fkInput(0),fOutput(), fDebug(0)
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
  if( fkInput && fDebug ){
    std::cout<<"HWCF Division: input Br: "<<fragment->fBranch<<" F: "<<fragment->fFlag<<" R: "<<fragment->fRow
	     <<" Q: "<<(fragment->fQ>>AliHLTTPCHWCFDefinitions::kFixedPoint)
	     <<" P: "<<fragment->fPad<<" Tmean: "<<fragment->fTMean;	        
    if( fragment->fFlag==1 && fragment->fQ > 0 ){
      std::cout<<" Pw: "<<((float)fragment->fP)/fragment->fQ
	       <<" Tw: "<<((float)fragment->fT)/fragment->fQ;
      std::cout<<"   MC: ";
      for( unsigned int j=0; j<fragment->fMC.size(); j++ ){
	for( int k=0; k<3; k++ ){
	  std::cout<<"("<<fragment->fMC[j].fClusterID[k].fMCID<<" "<<fragment->fMC[j].fClusterID[k].fWeight<<") ";
	}
      }
      std::cout<<std::endl;
    }
    else std::cout<<std::endl;      
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
  fOutput.fRowQ = (((AliHLTUInt32_t) 0x3)<<30) + ((fkInput->fRow &0x3f)<<24) + ((fkInput->fQmax)&0xFFFFFF);
  fOutput.fQ = fkInput->fQ;
  *((AliHLTFloat32_t*)&fOutput.fP) = (float)fkInput->fP/q;
  *((AliHLTFloat32_t*)&fOutput.fT) = (float)fkInput->fT/q;
  *((AliHLTFloat32_t*)&fOutput.fP2) = (float)fkInput->fP2/q;
  *((AliHLTFloat32_t*)&fOutput.fT2) = (float)fkInput->fT2/q;
 
  // MC part

  AliHLTTPCClusterMCWeight emptyWeight;

  fOutput.fMC.fClusterID[0] = emptyWeight;
  fOutput.fMC.fClusterID[1] = emptyWeight;
  fOutput.fMC.fClusterID[2] = emptyWeight;
  
  vector<AliHLTTPCClusterMCWeight> labels;
  for( unsigned i=0; i<fkInput->fMC.size(); i++){
    labels.push_back(fkInput->fMC[i].fClusterID[0]);
    labels.push_back(fkInput->fMC[i].fClusterID[1]);
    labels.push_back(fkInput->fMC[i].fClusterID[2]);
  }
  sort(labels.begin(), labels.end(), CompareMCLabels);
  for( unsigned int i=1; i<labels.size(); i++ ){
    if(labels[i-1].fMCID==labels[i].fMCID ){
      labels[i].fWeight+=labels[i-1].fWeight;
      labels[i-1].fWeight = 0;
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

