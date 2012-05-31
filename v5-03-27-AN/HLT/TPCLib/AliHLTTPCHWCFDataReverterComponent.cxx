// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt <Kenneth.Aamodt@student.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCHWCFDataReverterComponent.cxx
/// @author Kenneth Aamodt
/// @date   
/// @brief  Component for reverting data for the HW clusterfinder
///

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHWCFDataReverterComponent.h"
#include "AliHLTTPCDigitReader32Bit.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCMapping.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include "TString.h"
#include <sys/time.h>
#include "AliHLTAltroEncoder.h"
#include "AliRawDataHeader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCHWCFDataReverterComponent)

  AliHLTTPCHWCFDataReverterComponent::AliHLTTPCHWCFDataReverterComponent()
    :
    fDigitReader(NULL),
    fRowPadVector(),
    fNumberOfPadsInRow(NULL),
    fFirstPadHigh(NULL),
    fNumberOfRows(0),
    fCurrentPatch(0),
    fFirstRow(0),
    fLastRow(0),
    fNTimeBins(0),
    fVectorInitialized(kFALSE),
    fMapping(NULL),
    fInterleave(kTRUE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCHWCFDataReverterComponent::~AliHLTTPCHWCFDataReverterComponent()
{
  // see header file for class documentation
  if(fVectorInitialized){
    DeInitializePadArray();
  }
  if(fNumberOfPadsInRow){
    delete [] fNumberOfPadsInRow;
    fNumberOfPadsInRow=NULL;
  }
  if(fFirstPadHigh){
    delete [] fFirstPadHigh;
    fFirstPadHigh=NULL;
  }
  if(fDigitReader){
    delete fDigitReader;
    fDigitReader=NULL;
  }
  if(fMapping){
    delete fMapping;
    fMapping = NULL;
  }
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCHWCFDataReverterComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCHWCFDataReverter";
}

void AliHLTTPCHWCFDataReverterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCHWCFDataReverterComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeDDLRaw;
}

int AliHLTTPCHWCFDataReverterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeDDLRaw);
  return tgtList.size();
}

void AliHLTTPCHWCFDataReverterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTTPCHWCFDataReverterComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCHWCFDataReverterComponent();
}
	
int AliHLTTPCHWCFDataReverterComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- number of timebins
    if ( !strcmp( argv[i], "-timebins" )) {
      fNTimeBins = strtoul( argv[i+1], &cpErr ,0);
      AliHLTTPCTransform::SetNTimeBins(fNTimeBins);
      if ( *cpErr ) {
	HLTError("Cannot convert ntimebins specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
    // -- number of timebins
    if ( !strcmp( argv[i], "-interleave-off" )) {
      fInterleave = kFALSE;
      i++;
      continue;
    }
     
    HLTError("HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );

    return EINVAL;
  }

  fDigitReader = new AliHLTTPCDigitReader32Bit();

  return 0;
}

int AliHLTTPCHWCFDataReverterComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

Int_t AliHLTTPCHWCFDataReverterComponent::DeInitializePadArray()
{
  // see header file for class documentation
  if(fVectorInitialized){
    for(Int_t i=0;i<fNumberOfRows;i++){
      for(Int_t j=0;j<fNumberOfPadsInRow[i];j++){
	delete fRowPadVector[i][j];
	fRowPadVector[i][j]=NULL;
      }
      fRowPadVector[i].clear();
    }
    fRowPadVector.clear();
  }
  return 1;
} 

void AliHLTTPCHWCFDataReverterComponent::InitializePadArray(){
  // see header file for class documentation
  if(fCurrentPatch>5){
    HLTFatal("Patch is not set");
    return;
  }

  fFirstRow = AliHLTTPCTransform::GetFirstRow(fCurrentPatch);
  fLastRow = AliHLTTPCTransform::GetLastRow(fCurrentPatch);

  fNumberOfRows=fLastRow-fFirstRow+1;
  fNumberOfPadsInRow= new Int_t[fNumberOfRows];
  fFirstPadHigh= new Int_t[fNumberOfRows];

  memset( fNumberOfPadsInRow, 0, sizeof(Int_t)*(fNumberOfRows));
  memset( fFirstPadHigh, 0, sizeof(Int_t)*(fNumberOfRows));

  for(Int_t i=0;i<fNumberOfRows;i++){
    fNumberOfPadsInRow[i]=AliHLTTPCTransform::GetNPads(i+fFirstRow);
    AliHLTTPCPadVector tmpRow;
    for(Int_t j=0;j<fNumberOfPadsInRow[i];j++){
      AliHLTTPCPad *tmpPad = new AliHLTTPCPad();
      if(fFirstPadHigh[i] == 0){
	if(fMapping->GetHwAddress(i,j) > 2047){
	  fFirstPadHigh[i]=j;
	}
      }
      tmpPad->SetID(i,j);
      tmpPad->SetDataToDefault();
      tmpRow.push_back(tmpPad);
    }
    fRowPadVector.push_back(tmpRow);
  }
  fVectorInitialized=kTRUE;
}


int AliHLTTPCHWCFDataReverterComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						 const AliHLTComponentBlockData* blocks, 
						 AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
						 AliHLTUInt32_t& size, 
						 vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;
  if (!fDigitReader) return -ENODEV;

  AliHLTUInt32_t capacity=size;
  size=0;
  if (!IsDataEvent()) return 0;

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  //reading the data
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ){

    iter = blocks+ndx;
      
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	     evtData.fEventID, evtData.fEventID, 
	     DataType2Text( iter->fDataType).c_str(), 
	     DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)){
      continue;
    }

    if (iter->fSize<=sizeof(AliRawDataHeader)) {
      // forward empty DDLs
      outputBlocks.push_back(*iter);
      continue;
    }

    UInt_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    UInt_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

    if(fDigitReader->InitBlock(iter->fPtr,iter->fSize,patch,slice)<0){
      HLTWarning("Decoder failed to initialize, event aborted.");
      continue;
    }
      
    fMapping = new AliHLTTPCMapping(patch);

    if(!fVectorInitialized){
      fCurrentPatch=patch;
      InitializePadArray();
    }

    //Here the reading of the data and the zerosuppression takes place
    while(fDigitReader->NextChannel()){//Pad

      Int_t row=fDigitReader->GetRow();
      Int_t pad=fDigitReader->GetPad();

      if(row >= fNumberOfRows || row < 0){
	continue;
      }
      else if(pad >= fNumberOfPadsInRow[row] || pad < 0){
	continue;
      }  
	
      AliHLTTPCPad *tmpPad = fRowPadVector[row][pad];
      if (tmpPad){
	tmpPad->SetDataToDefault();
      }
	
      //reading data to pad
      while(fDigitReader->NextBunch()){
	const UInt_t *bunchData= fDigitReader->GetSignals();
	Int_t time=fDigitReader->GetTime();
	for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){
	  if(bunchData[i]>0){// disregarding 0 data.
	    if(time+i >= 0 && time+i < AliHLTTPCTransform::GetNTimeBins()){
	      if (tmpPad){
		tmpPad->SetDataSignal(time+i,bunchData[i]);
	      }
	    }
	  }
	}
      }
    }

    if( iter->fSize > sizeof(AliRawDataHeader )){
  
      AliHLTAltroEncoder *altroEncoder = new AliHLTAltroEncoder;
      altroEncoder->SetUse32BitFormat(kTRUE);
      Int_t ddlno=768;
      if (patch>1) ddlno+=72+4*slice+(patch-2);
      else ddlno+=2*slice+patch;
      altroEncoder->SetDDLid(ddlno);
      altroEncoder->SetSlice(slice);
      altroEncoder->SetPartition(patch);
      
      altroEncoder->SetBuffer(outputPtr,capacity); //tests if one overwrite the buffer is done in the encoder

      // set CDH from the beginning of buffer
      altroEncoder->SetCDH((AliHLTUInt8_t*)iter->fPtr,sizeof(AliRawDataHeader));

      UChar_t *RCUTrailer=NULL;
      Int_t RCUTrailerSize=fDigitReader->GetRCUTrailerSize();
      if (RCUTrailerSize<=0 || !fDigitReader->GetRCUTrailerData( RCUTrailer )) {
	if(RCUTrailer==NULL){
	  HLTWarning("can not find RCU trailer for data block %s 0x%08x: skipping data block",
		     DataType2Text(iter->fDataType).c_str(), iter->fSpecification);
	  continue;
	}
      }
      altroEncoder->SetRCUTrailer(RCUTrailer, RCUTrailerSize);

      for(Int_t row = 0; row< fNumberOfRows;row++){

	if(fInterleave == kFALSE){
	  Int_t currentTime = 0;
	  Int_t bunchSize = 0;
	  for(Int_t ipad=0;ipad<fNumberOfPadsInRow[row];ipad++){
	    AliHLTTPCPad * pad = fRowPadVector[row][ipad];
	    if(pad->GetNAddedSignals() > 0){
	      while(pad->GetNextGoodSignal(currentTime, bunchSize)){
		for(Int_t i=bunchSize-1;i>=0;i--){
		  if (altroEncoder->AddSignal((AliHLTUInt16_t)(pad->GetDataSignal(currentTime+i)),(AliHLTUInt16_t)(currentTime+i))<0) {
		    HLTWarning("can not add channel: slice %d, partition %d, hw address %d, row %d, pad %d, time %d, bunch size %d Charge %d",
			       slice, patch, fMapping->GetHwAddress(row,ipad), row, ipad, currentTime+i, bunchSize,pad->GetDataSignal(currentTime+i));
		    break;
		  }
		}
	      }
	      altroEncoder->SetChannel(fMapping->GetHwAddress(row,ipad));
	      currentTime = 0;
	      bunchSize = 0;
	    }
	  }
	}
	else{
	  Int_t padHigh=fFirstPadHigh[row];
	  
	  Int_t padLowIndex=0;
	  Int_t padHighIndex= padHigh;
	  
	  while(padLowIndex < padHigh || padHighIndex < fNumberOfPadsInRow[row]){
	    Int_t currentTime = 0;
	    Int_t bunchSize = 0;
	    //add the data from low side
	    if(padLowIndex < padHigh){
	      AliHLTTPCPad * lowPad= fRowPadVector[row][padLowIndex];
	      if(lowPad->GetNAddedSignals()>0){
		while(lowPad->GetNextGoodSignal(currentTime, bunchSize)){
		  for(Int_t i=bunchSize-1;i>=0;i--){
		    if (altroEncoder->AddSignal((AliHLTUInt16_t)(lowPad->GetDataSignal(currentTime+i)),(AliHLTUInt16_t)(currentTime+i))<0) {
		      HLTWarning("can not add channel: slice %d, partition %d, hw address %d, row %d, pad %d, time %d, bunch size %d",
				 slice, patch, fMapping->GetHwAddress(row,padLowIndex), row, padLowIndex, currentTime+i, bunchSize);
		      break;
		    }
		  }
		}
		altroEncoder->SetChannel(fMapping->GetHwAddress(row,padLowIndex));
	      }
	    }
	    currentTime = 0;
	    bunchSize = 0;
	    //add the data from the high side
	    if(padHighIndex < fNumberOfPadsInRow[row]){
	      AliHLTTPCPad * highPad= fRowPadVector[row][padHighIndex];
	      if(highPad->GetNAddedSignals()>0){
		while(highPad->GetNextGoodSignal(currentTime, bunchSize)){
		  for(Int_t i=bunchSize-1;i>=0;i--){
		    if (altroEncoder->AddSignal((AliHLTUInt16_t)(highPad->GetDataSignal(currentTime+i)),(AliHLTUInt16_t)(currentTime+i))<0) {
		      HLTWarning("can not add channel: slice %d, partition %d, hw address %d, row %d, pad %d, time %d, bunch size %d",
				 slice, patch, fMapping->GetHwAddress(row,padHighIndex), row, padHighIndex, currentTime+i, bunchSize);
		      break;
		    }
		  }
		}
		altroEncoder->SetChannel(fMapping->GetHwAddress(row,padHighIndex));
	      }
	    }
	    padLowIndex++;
	    padHighIndex++;
	  }
	}
      }	
      
      int sizeOfData=altroEncoder->SetLength();
      

      if (sizeOfData<0) {
	HLTError("data encoding failed");
	iResult=sizeOfData;
	break;
      }
      if(sizeOfData>(int)capacity){
	HLTWarning("Buffer too small too add the altrodata: %d of %d byte(s) already used", sizeOfData, capacity);
	iResult=-ENOSPC;
	break;
      }

      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = 0;
      bd.fSize = sizeOfData;
      bd.fDataType = kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC;
      bd.fSpecification = iter->fSpecification;
      outputBlocks.push_back( bd );

      size+=sizeOfData;
    }
    fDigitReader->Reset();
  }
  
  if(iResult < 0) {
    fDigitReader->Reset();
    size=0;
  }
  HLTDebug("Total size of output is: %d ",size);
  return iResult;
}
