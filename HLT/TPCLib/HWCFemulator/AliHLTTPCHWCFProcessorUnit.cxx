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

//  @file   AliHLTTPCHWCFProcessorUnit.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @date   
//  @brief  Channel Processor unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note 

#include "AliHLTTPCHWCFProcessorUnit.h"
#include <iostream>


AliHLTTPCHWCFProcessorUnit::AliHLTTPCHWCFProcessorUnit()
  :
  fOutput(),
  fkBunch(0),
  fBunchIndex(0),
  fDeconvolute(0),
  fSingleSeqLimit(0)
{
  //constructor 
  Init();
}


AliHLTTPCHWCFProcessorUnit::~AliHLTTPCHWCFProcessorUnit()
{   
  //destructor 
}

AliHLTTPCHWCFProcessorUnit::AliHLTTPCHWCFProcessorUnit(const AliHLTTPCHWCFProcessorUnit&)
  :
  fOutput(),
  fkBunch(0),
  fBunchIndex(0),
  fDeconvolute(0),
  fSingleSeqLimit(0)
{
  // dummy
  Init();
}

AliHLTTPCHWCFProcessorUnit& AliHLTTPCHWCFProcessorUnit::operator=(const AliHLTTPCHWCFProcessorUnit&)
{
  // dummy  
  return *this;
}

int AliHLTTPCHWCFProcessorUnit::Init()
{
  // initialise  

  fkBunch = 0;
  return 0;
}

int AliHLTTPCHWCFProcessorUnit::InputStream( const AliHLTTPCHWCFBunch *bunch )
{
  // input stream of data 
  
  //if( bunch && bunch->fRow==0 ){
  //std::cout<<"Processor: input bunch F "<<bunch->fFlag<<" R "<<bunch->fRow<<" P "<<bunch->fPad <<" Br "<<bunch->fBranch<<" bord "<<bunch->fBorder
  //<<" Time "<<bunch->fTime<<" NS= "<<bunch->fData.size()<<std::endl;
  //}
  fkBunch = bunch;
  fBunchIndex = 0;
  return 0;
}

const AliHLTTPCHWCFClusterFragment *AliHLTTPCHWCFProcessorUnit::OutputStream()
{ 
  // output stream of data 

  if( !fkBunch ) return 0;
  
  fOutput.fFlag = fkBunch->fFlag;
  fOutput.fRow = fkBunch->fRow;
  fOutput.fPad = fkBunch->fPad;
  fOutput.fBranch = fkBunch->fBranch;
  fOutput.fBorder = fkBunch->fBorder;  
  fOutput.fQ = 0;
  fOutput.fT = 0;
  fOutput.fT2 = 0;
  fOutput.fP = 0;
  fOutput.fP2 = 0;
  fOutput.fTMean =  fkBunch->fTime;

  fOutput.fMC.clear();
  
  if( fkBunch->fFlag==2 && fkBunch->fData.size()==1 ){ // rcu trailer word, forward it 
    fOutput.fRow = fkBunch->fData[0];
  }	
  
  if( fkBunch->fFlag >1 ){
    fkBunch = 0;
    fBunchIndex = 0;
    return &fOutput;
  }

  if( fkBunch->fFlag < 1 ) return 0;


  if( fBunchIndex >= fkBunch->fData.size() || fkBunch->fTime < fBunchIndex ) return 0;  

  AliHLTInt32_t bunchTime0 = fkBunch->fTime - fBunchIndex;
  AliHLTInt32_t bunchTime = bunchTime0;

  AliHLTUInt64_t qLast = 0;
  bool slope = 0;
  AliHLTUInt32_t length = 0;
  for( ; fBunchIndex<fkBunch->fData.size() && bunchTime>=0; fBunchIndex++, bunchTime--, length++ ){
    AliHLTUInt64_t q = fkBunch->fData[fBunchIndex]*fkBunch->fGain;
    if( fDeconvolute && slope && q>qLast ){
      cout<<"deconvolution time!!!"<<endl;
      if( length==1 && fOutput.fQ<fSingleSeqLimit ){
	fOutput.fQ = 0;
	fOutput.fT = 0;
	fOutput.fT2 = 0;
	fOutput.fP = 0;
	fOutput.fP2 = 0;
	bunchTime0 = fkBunch->fTime - fBunchIndex;
	qLast = 0;
	slope = 0;
	length = 0;
      } else {      
	break;
      }
    }
    if( q<qLast ) slope = 1;
    qLast = q;
    fOutput.fQ += q;
    fOutput.fT += q*bunchTime;
    fOutput.fT2+= q*bunchTime*bunchTime;
    fOutput.fP += q*fkBunch->fPad;
    fOutput.fP2+= q*fkBunch->fPad*fkBunch->fPad;
    fOutput.fMC.insert(fOutput.fMC.end(),fkBunch->fMC.begin(), fkBunch->fMC.end() );
  }
  
  fOutput.fTMean = (AliHLTUInt64_t)( (bunchTime0 + bunchTime + 1)/2 );

  if( length==1 && fOutput.fQ < fSingleSeqLimit ) return 0;

  return &fOutput;
}
