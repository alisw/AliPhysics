// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kenneth Aamodt <Kenneth.Aamodt@student.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCZeroSuppressionComponent.cxx
    @author Kenneth Aamodt
    @date   
    @brief  The TPC ZeroSuppression component
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCZeroSuppressionComponent.h"
#include "AliHLTTPCDigitReaderDecoder.h"
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
ClassImp(AliHLTTPCZeroSuppressionComponent)

AliHLTTPCZeroSuppressionComponent::AliHLTTPCZeroSuppressionComponent()
    :
    fDigitReader(NULL),
    fRowPadVector(),
    fNumberOfPadsInRow(NULL),
    fNumberOfRows(0),
    fCurrentPatch(0),
    fFirstRow(0),
    fLastRow(0),
    fStartTimeBin(0),
    fEndTimeBin(AliHLTTPCTransform::GetNTimeBins()),
    fNTimeBins(0),
    fNRMSThreshold(0),
    fSignalThreshold(0),
    fMinimumNumberOfSignals(AliHLTTPCTransform::GetNTimeBins()/2),
    fOldRCUFormat(0),
    fSortPads(0),
    fVectorInitialized(kFALSE),
    fValueBelowAverage(5),
    fLeftTimeBin(5),
    fRightTimeBin(5),
    fGetActivePads(kFALSE),
    fSkipSendingZSData(kFALSE),
    fSendHWList(kFALSE),
    fHwAddressList()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCZeroSuppressionComponent::~AliHLTTPCZeroSuppressionComponent()
{
  // see header file for class documentation
  if(fVectorInitialized){
    DeInitializePadArray();
  }
  if(fNumberOfPadsInRow){
    delete [] fNumberOfPadsInRow;
    fNumberOfPadsInRow=NULL;
  }
  if(fDigitReader){
    delete fDigitReader;
    fDigitReader=NULL;
  }
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCZeroSuppressionComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCZeroSuppression";
}

void AliHLTTPCZeroSuppressionComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCZeroSuppressionComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
  //return kAliHLTDataTypeDDLRaw;
  //  return AliHLTTPCDefinitions::fgkUnpackedRawDataType;
}

int AliHLTTPCZeroSuppressionComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeDDLRaw);
  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  return tgtList.size();
}

void AliHLTTPCZeroSuppressionComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=2.0;
}

AliHLTComponent* AliHLTTPCZeroSuppressionComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCZeroSuppressionComponent();
}
	
int AliHLTTPCZeroSuppressionComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- zero suppression threshold
    if ( !strcmp( argv[i], "-signal-threshold" ) || !strcmp( argv[i], "signal-threshold" ) ) {
      fSignalThreshold = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert signal-threshold specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for nsigma-threshold, used in 2007 December run in ZeroSuppression
    if ( !strcmp( argv[i], "-rms-threshold" ) ||  !strcmp( argv[i], "rms-threshold" ) ) {
      fNRMSThreshold = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert rms-threshold specifier '%s'. Must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- number of timebins
    if ( !strcmp( argv[i], "-timebins" ) || !strcmp( argv[i], "ntimebins" ) || !strcmp( argv[i], "-ntimebins" )) {
      fNTimeBins = strtoul( argv[i+1], &cpErr ,0);
      AliHLTTPCTransform::SetNTimeBins(fNTimeBins);
      if(fEndTimeBin>AliHLTTPCTransform::GetNTimeBins()){
	fEndTimeBin = AliHLTTPCTransform::GetNTimeBins();
      }
      if ( *cpErr ) {
	HLTError("Cannot convert ntimebins specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- first timebin
    if ( !strcmp( argv[i], "-start-timebin" ) || !strcmp( argv[i], "start-timebin" ) ) {
      fStartTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert start-timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- last timebin
    if ( !strcmp( argv[i], "-end-timebin" ) || !strcmp( argv[i], "end-timebin" ) ) {
      if(strtoul( argv[i+1], &cpErr ,0)<=(UInt_t)AliHLTTPCTransform::GetNTimeBins()){
	fEndTimeBin = strtoul( argv[i+1], &cpErr ,0);
      }
      if ( *cpErr ) {
	HLTError("Cannot convert end-timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- timebins to keep left of signal
    if ( !strcmp( argv[i], "-timebin-left" ) || !strcmp( argv[i], "timebin-left" ) ) {
      fLeftTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert timebin-left specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- timebin to keep right of signal
    if ( !strcmp( argv[i], "-timebin-right" ) || !strcmp( argv[i], "timebin-right" ) ) {
      fRightTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert timebin-right specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- value below average to subtract
    if ( !strcmp( argv[i], "-value-below-average" ) || !strcmp( argv[i], "value-below-average" ) ) {
      fValueBelowAverage = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert value-below-average specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[i], "-occupancy-limit" ) || !strcmp( argv[i], "occupancy-limit" ) ) {
      fMinimumNumberOfSignals = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert occupancy-limit specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for rcu format
    if ( !strcmp( argv[i], "oldrcuformat" ) ) {
      fOldRCUFormat = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert oldrcuformat specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for pad sorting
    if ( !strcmp( argv[i], "-sort-pads" ) || !strcmp( argv[i], "sort-pads" ) ) {
      fSortPads = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert sort-pads specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for skipZSdatashipping
    if ( !strcmp( argv[i], "-skip-sending-data" ) ) {
      fSkipSendingZSData = kTRUE;
      i++;
      continue;
    }

    // -- checking for hw address shipping
    if ( !strcmp( argv[i], "-send-hw-list" ) ) {
      fSendHWList = kTRUE;
      i++;
      continue;
    }
      
    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  if(fSkipSendingZSData == kTRUE && fSendHWList == kFALSE){
    HLTError("Component will have no output, check your configuration.");
  }
  

  HLTDebug("using AliHLTTPCDigitReaderDecoder");
  fDigitReader = new AliHLTTPCDigitReaderDecoder();

  fHwAddressList.clear();

  return 0;
}

int AliHLTTPCZeroSuppressionComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

Int_t AliHLTTPCZeroSuppressionComponent::DeInitializePadArray()
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

void AliHLTTPCZeroSuppressionComponent::InitializePadArray(){
  // see header file for class documentation
  if(fCurrentPatch>5){
    HLTFatal("Patch is not set");
    return;
  }

  fFirstRow = AliHLTTPCTransform::GetFirstRow(fCurrentPatch);
  fLastRow = AliHLTTPCTransform::GetLastRow(fCurrentPatch);

  fNumberOfRows=fLastRow-fFirstRow+1;
  fNumberOfPadsInRow= new Int_t[fNumberOfRows];

  memset( fNumberOfPadsInRow, 0, sizeof(Int_t)*(fNumberOfRows));

  for(Int_t i=0;i<fNumberOfRows;i++){
    fNumberOfPadsInRow[i]=AliHLTTPCTransform::GetNPads(i+fFirstRow);
    AliHLTTPCPadVector tmpRow;
    for(Int_t j=0;j<fNumberOfPadsInRow[i];j++){
      AliHLTTPCPad *tmpPad = new AliHLTTPCPad();
      tmpPad->SetID(i,j);
      tmpRow.push_back(tmpPad);
    }
    fRowPadVector.push_back(tmpRow);
  }
  fVectorInitialized=kTRUE;
}


int AliHLTTPCZeroSuppressionComponent::DoEvent( const AliHLTComponentEventData& evtData, 
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

  //  HLTInfo("Entered DoEvent in AliHLTTPCZeroSuppressionComponent");

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  //  HLTInfo("Number of blocks: ",evtData.fBlockCnt);

  fHwAddressList.clear();
  //reading the data
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text( iter->fDataType).c_str(), 
	       DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

      if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)){
	continue;
      }

      UInt_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      UInt_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

      if(!fVectorInitialized){
	fCurrentPatch=patch;
	InitializePadArray();
      }
      
      if(fDigitReader->InitBlock(iter->fPtr,iter->fSize,patch,slice)<0){
	HLTWarning("Decoder failed to initialize, event aborted.");
	continue;
      }
      AliHLTTPCMapping mapping(patch);

      int nTotalChannels=0;
      int nSkippedChannels=0;
      short lowestOccupancy=-1;
      float sumOccupancy=0;

      //Here the reading of the data and the zerosuppression takes place
      while(fDigitReader->NextChannel()){//Pad
	int sumSignals=0;
	int maxSignal=0;
	int nofSignals=0;
	Int_t row=fDigitReader->GetRow();
	Int_t pad=fDigitReader->GetPad();
	if(row==1000 || pad==1000){
	  continue;
	}
	if(row>=fNumberOfRows||row<0){
	  continue;
	}
	else if(pad>=fNumberOfPadsInRow[row]||pad<0){
	  continue;
	}  
	
	AliHLTTPCPad *tmpPad = NULL;
	if (!fSkipSendingZSData) tmpPad=fRowPadVector[row][pad];
	if (tmpPad) tmpPad->SetDataToDefault();

	//reading data to pad
	while(fDigitReader->NextBunch()){
	  const UInt_t *bunchData= fDigitReader->GetSignals();
	  Int_t time=fDigitReader->GetTime();
	  for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){
	    if(bunchData[i]>0){// disregarding 0 data.
	      if(time+i>=fStartTimeBin && time+i<=fEndTimeBin){
		if (tmpPad) tmpPad->SetDataSignal(time+i,bunchData[i]);
		sumSignals+=bunchData[i];
		if (maxSignal<(int)bunchData[i]) maxSignal=bunchData[i];
		nofSignals++;
	      }
	    }
	  }
	}

	nTotalChannels++;
	if (lowestOccupancy<0 || lowestOccupancy>nofSignals)
	  lowestOccupancy=nofSignals;
	sumOccupancy+=nofSignals;

	if(nofSignals>=fMinimumNumberOfSignals){
	  if (tmpPad) {
	  tmpPad->ZeroSuppress(fNRMSThreshold, fSignalThreshold, fMinimumNumberOfSignals, fStartTimeBin, fEndTimeBin, fLeftTimeBin, fRightTimeBin, fValueBelowAverage, fSkipSendingZSData);
	  if(tmpPad->GetNAddedSignals()>0){
	    assert((int)mapping.GetRow(fDigitReader->GetAltroBlockHWaddr())==row);
	    assert((int)mapping.GetPad(fDigitReader->GetAltroBlockHWaddr())==pad);
	    fHwAddressList.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr());
	  }
	  } else {
	    assert(fSkipSendingZSData);
	    if (nofSignals>0 && maxSignal>(sumSignals/nofSignals)+fSignalThreshold) {
	      fHwAddressList.push_back((AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr());
	    }
	  }
	} else {
	  nSkippedChannels++;
	}
      }
      if (nSkippedChannels>0) {
	HLTWarning("skipped %d of %d channels because of low occupancy: average %.2f, lowest %d, threshold %d",
		   nSkippedChannels, nTotalChannels, sumOccupancy/nTotalChannels, lowestOccupancy, fMinimumNumberOfSignals);
      }

      AliHLTUInt32_t dataOffsetBeforeHW=0;

      if(fSkipSendingZSData == kFALSE && iter->fSize>sizeof(AliRawDataHeader)) {
  
	AliHLTAltroEncoder *altroEncoder = new AliHLTAltroEncoder;
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

	for(unsigned int channel=0; channel<fHwAddressList.size(); channel++){
	  int row=mapping.GetRow(fHwAddressList[channel]);
	  int pad=mapping.GetPad(fHwAddressList[channel]);
	  if (true) {
	    AliHLTTPCPad * zeroSuppressedPad= fRowPadVector[row][pad];
	    Int_t currentTime=0;
	    Int_t bunchSize=0;
	    if(zeroSuppressedPad->GetNAddedSignals()>0){
	      while(zeroSuppressedPad->GetNextGoodSignal(currentTime, bunchSize)){
		for(Int_t i=0;i<bunchSize;i++){
		  if (altroEncoder->AddSignal((AliHLTUInt16_t)(zeroSuppressedPad->GetDataSignal(currentTime+i)),(AliHLTUInt16_t)(currentTime+i))<0) {
		    // Matthias 01.10.2008: there is a problem with certain real data which produces the same
		    // bunch multiple times, needs investigation. I found an examplary case in run 53465
		    // (08000053465011.450.root) equipment 981.
		    // needs to be followed up.
		    // addon 02.10.2008 I just corrected a bug concerning the loop over the active channels.
		    // A double loop over all rows and pads also considered pads which had not even been
		    // filled, and thus not properly cleaned. Maybe the bug above is related to that.
		    HLTWarning("can not add channel: slice %d, partition %d, hw address %d, row %d, pad %d, time %d, bunch size %d",
			       slice, patch, fHwAddressList[channel], row, pad, currentTime+i, bunchSize);
		    break;
		  }
		}
	      }
	      altroEncoder->SetChannel(fHwAddressList[channel]);
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

	//Push back the zerosuppressed altro data to the output
	if(true/* condition was deprecated but keep for the sake of formatting*/){
	  AliHLTComponentBlockData bd;
	  FillBlockData( bd );
	  bd.fOffset = 0;
	  bd.fSize = sizeOfData;
	  bd.fDataType = kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC;
	  bd.fSpecification = iter->fSpecification;
	  outputBlocks.push_back( bd );
    
	  //Push back the list of hardware addresses to the output
	  dataOffsetBeforeHW=sizeOfData;
	}
      }

      AliHLTUInt32_t sizeOfHWArray = 0;
      if(fSendHWList == kTRUE){

	if(true/* condition was deprecated but keep for the sake of formatting*/){
	  sizeOfHWArray = fHwAddressList.size()*sizeof(AliHLTUInt16_t);
      
	  if(dataOffsetBeforeHW+sizeOfHWArray>capacity){
	    HLTWarning("Buffer too small too add the active channels: %d of %d byte(s) already used", dataOffsetBeforeHW + sizeOfHWArray, capacity);
	    iResult=-ENOSPC;
	    break;
	  }
      
	  AliHLTUInt16_t*outputHWPtr=(AliHLTUInt16_t*)(outputPtr+dataOffsetBeforeHW);
	  memcpy(outputHWPtr,&fHwAddressList[0],sizeOfHWArray);
	  AliHLTComponentBlockData bdHW;
	  FillBlockData( bdHW );
	  bdHW.fOffset = dataOffsetBeforeHW;
	  bdHW.fSize = sizeOfHWArray;
	  bdHW.fDataType = kAliHLTDataTypeHwAddr16|kAliHLTDataOriginTPC;
	  bdHW.fSpecification = iter->fSpecification;
	  outputBlocks.push_back( bdHW );
	}
      }
      size = dataOffsetBeforeHW+sizeOfHWArray;
      fDigitReader->Reset();
    }

  if (iResult<0) {
    fDigitReader->Reset();
    size=0;
  }

  return iResult;
}
