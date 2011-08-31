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

//  @file   AliHLTTPCHWCFExtractorUnit.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @date  
//  @brief  Channel Extractor unit of FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note 

#include "AliHLTTPCHWCFExtractorUnit.h"

AliHLTTPCHWCFExtractorUnit::AliHLTTPCHWCFExtractorUnit()
  :
  fStatus(kReadingData),
  fInput(0),
  fInputStatus(kEmpty),
  fkMapping(0),  
  fBunch(fMemory),
  fPendingOutput(0),
  fChannelNumWordsLeft(0),
  fBunchNumWordsLeft(0),  
  fkMCLabels(0),
  fNMCLabels(0),
  fCurrentMCLabel(0)
{
  //constructor 
  fBunch->fFlag = 0; // wait for the next channel
}
 

AliHLTTPCHWCFExtractorUnit::~AliHLTTPCHWCFExtractorUnit()
{   
  //destructor  
}

AliHLTTPCHWCFExtractorUnit::AliHLTTPCHWCFExtractorUnit(const AliHLTTPCHWCFExtractorUnit&)
  :
  fStatus(kReadingData),
  fInput(0),
  fInputStatus(kEmpty),
  fkMapping(0),
  fBunch(fMemory),
  fPendingOutput(0),
  fChannelNumWordsLeft(0),
  fBunchNumWordsLeft(0),  
  fkMCLabels(0),
  fNMCLabels(0),
  fCurrentMCLabel(0)
{
  // dummy
}

AliHLTTPCHWCFExtractorUnit& AliHLTTPCHWCFExtractorUnit::operator=(const AliHLTTPCHWCFExtractorUnit&)
{
  // dummy  
  return *this;
}

int AliHLTTPCHWCFExtractorUnit::Init( const AliHLTUInt32_t *mapping, const AliHLTTPCClusterMCLabel *mcLabels, AliHLTUInt32_t nMCLables )
{
  // initialisation

  fStatus = kReadingData;
  fInput = 0;
  fInputStatus = kEmpty;
  fkMapping = mapping;
  fBunch->fFlag = 0; // wait for the next channel
  fBunch->fData.clear();
  fBunch->fMC.clear();
  fPendingOutput = 0;  
  fkMCLabels = mcLabels;
  fNMCLabels = nMCLables;
  fCurrentMCLabel = 0;
  if( !fkMapping ) return  -1;
  return 0;
}

int AliHLTTPCHWCFExtractorUnit::InputEndOfData()
{
  // input "end of data" signal 

  int ret = 0;
  if( fInputStatus != kEmpty ) ret = -1;
  fInputStatus = kEndOfData;
  return ret;
}

int AliHLTTPCHWCFExtractorUnit::InputStream( AliHLTUInt32_t word )
{
  // input stream of data

  int ret = 0;
  if( fInputStatus != kEmpty ) ret = -1;
  fInputStatus = kData;
  fInput = word;
  return ret;
}

const AliHLTTPCHWCFBunch *AliHLTTPCHWCFExtractorUnit::OutputStream()
{
  // output stream of data
  
  AliHLTTPCHWCFExtractorInputStatus inpStatus = fInputStatus;
  fInputStatus = kEmpty;

  if( fPendingOutput ){
    fPendingOutput = 0;
    return fBunch;
  }

  if( fStatus == kStop ) return 0;  
 
  if( !fkMapping ) inpStatus = kEndOfData;

  AliHLTTPCHWCFBunch *oldBunch = fBunch;
  AliHLTTPCHWCFBunch *newBunch = ( fBunch == fMemory ) ?fMemory+1 :fMemory;

  if( fStatus == kFinishing || inpStatus == kEndOfData){
    if( fBunch->fFlag == 1 ){
      fBunch = newBunch;
      fPendingOutput = 1;
    }
    fBunch->fData.clear();
    fBunch->fMC.clear();
    fBunch->fFlag = 3; // end of data
    fStatus = kStop;    
    return oldBunch;   
  }

  if( inpStatus == kEmpty ) return 0;  

  AliHLTUInt32_t flag = (fInput >> 30); 
  
  if( flag >= 0x2 ){ //  RCU trailer  
    if( fBunch->fFlag == 1 ){
      fBunch = newBunch;
      fPendingOutput = 1;
    }
    fBunch->fData.clear();
    fBunch->fMC.clear();
    fBunch->fFlag = 2; // rcu
    fBunch->fData.push_back(fInput);
    fStatus = ( flag == 0x2 ) ?kReadingRCU :kFinishing;
    return oldBunch;   
  }
  
  if( fStatus!=kReadingData ) return 0;  
  
  if( flag==0x1 ){ // header of a new channel
        
    if( fBunch->fFlag==1 ){
      //cout<<"Extractor: Bunch finished: "<<fBunch->fFlag
      //<<" R "<<fBunch->fRow<<" P "<<fBunch->fPad<<" T "<<fBunch->fTime<<" NS "<<fBunch->fNSignals<<endl;
      fBunch = newBunch;
    }
    AliHLTUInt32_t  hwAddress = fInput & 0xFFF;

    fBunch->fFlag = 1;
    fBunch->fBranch = (hwAddress >> 11) & 0x1;
    if( hwAddress>=fkMapping[0] ) fBunch->fFlag = 0; //readout errors
    else{
      AliHLTUInt32_t configWord = fkMapping[hwAddress+1];
      fBunch->fRow = (configWord>>8) & 0x3F;
      fBunch->fPad =  configWord & 0xFF;
      fBunch->fBorder = (configWord>>14) & 0x1;
      if( !( (configWord>>15) & 0x1 ) ) fBunch->fFlag = 0;// channel not active
      fBunch->fGain = (configWord>>16 ) & 0x1FFF;
    }
    fBunch->fData.clear();
    fBunch->fMC.clear();
    fBunch->fTime = 0xFFFFFFFF;
    fChannelNumWordsLeft= (fInput >> 16) & 0x3FF; // payload size in 10-bit words
    fBunchNumWordsLeft = 0;
   
    if( (fInput >> 29) & 0x1 ) fBunch->fFlag = 0; // there were readout errors

    //cout<<"Extractor: Header of new channel F "<<fBunch->fFlag
    //<<" R "<<fBunch->fRow<<" P "<<fBunch->fPad<<" NWords10 "<<fChannelNumWordsLeft<<endl;
  
  } else if( flag==0x0 ){ 
    
    // bunch data, read three 10-bit words
    
    for( int ishift=20; fBunch->fFlag==1 && ishift>=0; ishift-=10 ){
      
      AliHLTUInt32_t word10 = (fInput >> ishift) & 0x3FF;
      
      if( fChannelNumWordsLeft <= 0 || fBunchNumWordsLeft <= 0 ){ // bunch finished
	//cout<<"Extractor: Bunch finished: "<<fBunch->fFlag
	//<<" R "<<fBunch->fRow<<" P "<<fBunch->fPad<<" T "<<fBunch->fTime<<" NS "<<fBunch->fNSignals<<endl;	
	fBunch = newBunch; // push to the output
	
	if( fChannelNumWordsLeft <= 0 ){ // wait for the next channel
	  fBunch->fFlag = 0;
	} else {
	  fBunch->fFlag = 1;
	  fBunch->fRow  = oldBunch->fRow;
	  fBunch->fPad = oldBunch->fPad;
	  fBunch->fBranch = oldBunch->fBranch;
	  fBunch->fBorder = oldBunch->fBorder;
	  fBunch->fGain = oldBunch->fGain;
	  fBunch->fData.clear();
	  fBunch->fMC.clear();	  
	  fBunch->fTime = 0xFFFFFFFF;
	  fBunchNumWordsLeft = word10;
	}
      } else { // continue the brunch
	if( fBunch->fTime > AliHLTTPCHWCFDefinitions::kMaxNTimeBins ){ // time has not been read so far
	  fBunch->fTime = word10;
	  //cout<<"Extractor: Bunch time: "<<fBunch->fTime<<endl;
	} else { // read the signal
	  fBunch->fData.push_back(word10);
	  if( fkMCLabels && fCurrentMCLabel<=fNMCLabels ){
	    fBunch->fMC.push_back( fkMCLabels[fCurrentMCLabel] );
	    fCurrentMCLabel++;
	  }
	  //cout<<"Extractor: Bunch signal["<<fBunch->fNSignals<<"]: "<<word10<<endl;
	}
      }
      fChannelNumWordsLeft--;
      fBunchNumWordsLeft--;   
    }
  }

  if( fBunch==newBunch && oldBunch->fFlag==1 && oldBunch->fData.size()>0 ){
    return oldBunch;
  }
  return 0;
}
