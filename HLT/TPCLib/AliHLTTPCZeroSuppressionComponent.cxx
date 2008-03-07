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
#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include <sys/time.h>


AliHLTTPCZeroSuppressionComponent gAliHLTTPCZeroSuppressionComponent;

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
  return AliHLTTPCDefinitions::fgkUnpackedRawDataType;
}

int AliHLTTPCZeroSuppressionComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkUnpackedRawDataType);
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
    if ( !strcmp( argv[i], "signal-threshold" ) ) {
      fSignalThreshold = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert signal-threshold specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for nsigma-threshold, used in 2007 December run in ZeroSuppression
    if ( !strcmp( argv[i], "rms-threshold" ) ) {
      fNRMSThreshold = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert rms-threshold specifier '%s'. Must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- number of timebins
    if ( !strcmp( argv[i], "ntimebins" ) ) {
      fNTimeBins = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert ntimebins specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- first timebin
    if ( !strcmp( argv[i], "start-timebin" ) ) {
      fStartTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert start-timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- last timebin
    if ( !strcmp( argv[i], "end-timebin" ) ) {
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
    if ( !strcmp( argv[i], "timebin-left" ) ) {
      fLeftTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert timebin-left specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- timebin to keep right of signal
    if ( !strcmp( argv[i], "timebin-right" ) ) {
      fRightTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert timebin-right specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- value below average to subtract
    if ( !strcmp( argv[i], "value-below-average" ) ) {
      fValueBelowAverage = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert value-below-average specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[i], "occupancy-limit" ) ) {
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

    // -- checking for rcu format
    if ( !strcmp( argv[i], "sort-pads" ) ) {
      fSortPads = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert sort-pads specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
      
    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

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

  //  HLTInfo("Entered DoEvent in AliHLTTPCZeroSuppressionComponent");

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  //  HLTInfo("Number of blocks: ",evtData.fBlockCnt);

  Bool_t wasInput = 0;

  fHwAddressList.clear();
  //reading the data
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text( iter->fDataType).c_str(), 
	       DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

      if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType &&
	  GetEventCount()<2) {
	HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!",
		   DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
		   DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
      }
      
      if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) &&
	   iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ){
	continue;
      }

      wasInput = 1;

      UInt_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      UInt_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

      if(!fVectorInitialized){
	fCurrentPatch=patch;
	InitializePadArray();
      }
      
      fDigitReader->InitBlock(iter->fPtr,iter->fSize,patch,slice);

      //Here the reading of the data and the zerosuppression takes place
      while(fDigitReader->NextChannel()){//Pad
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
	
	AliHLTTPCPad *tmpPad = fRowPadVector[row][pad];
	tmpPad->SetDataToDefault();
	
	//reading data to pad
	while(fDigitReader->NextBunch()){
	  const UInt_t *bunchData= fDigitReader->GetSignals();
	  Int_t time=fDigitReader->GetTime();
	  for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){
	    if(bunchData[i]>0){// disregarding 0 data.
	      if(time+i>=fStartTimeBin && time+i<=fEndTimeBin){
		tmpPad->SetDataSignal(time+i,bunchData[i]);
	      }
	    }
	  }
	}
	if(tmpPad->GetNAddedSignals()>=(UInt_t)fMinimumNumberOfSignals){
	  fHwAddressList.push_back(fDigitReader->GetAltroBlockHWaddr());
	  tmpPad->ZeroSuppress(fNRMSThreshold, fSignalThreshold, fMinimumNumberOfSignals, fStartTimeBin, fEndTimeBin, fLeftTimeBin, fRightTimeBin, fValueBelowAverage);
	}
      }
    }
  Int_t nAdded=0;

  HLTDebug("Max numberof signals: %d",size/sizeof(Int_t));

  //writing to output
  AliHLTUInt8_t* outBPtr;
  UInt_t* outPtr;
  unsigned long long outputSize = 0;
  unsigned long blockOutputSize = 0;
  outBPtr = outputPtr;
  outPtr = (UInt_t*)outputPtr;

  outPtr[nAdded]=0;
  UInt_t *numberOfChannels=&outPtr[nAdded];
  nAdded++;
  outputSize      += sizeof(UInt_t);
  blockOutputSize += sizeof(UInt_t);
  for(Int_t row=0;row<fNumberOfRows;row++){
    for(Int_t pad=0;pad<fNumberOfPadsInRow[row];pad++){
      AliHLTTPCPad * zeroSuppressedPad= fRowPadVector[row][pad];
      Int_t currentTime=0;
      Int_t bunchSize=0;
      Bool_t newPad=kTRUE;
      UInt_t *nBunches=NULL;
      while(zeroSuppressedPad->GetNextGoodSignal(currentTime, bunchSize,0)){
	if(newPad){
   	  (*numberOfChannels)++;
	  
	  outPtr[nAdded]=(UInt_t)row;
	  nAdded++;
	  outputSize      += sizeof(UInt_t);
	  blockOutputSize += sizeof(UInt_t);

	  outPtr[nAdded]=(UInt_t)pad;
	  nAdded++;
	  outputSize      += sizeof(UInt_t);
	  blockOutputSize += sizeof(UInt_t);

	  
	  outPtr[nAdded]=0;
	  nBunches=&outPtr[nAdded];
	  nAdded++;
	  outputSize      += sizeof(UInt_t);
	  blockOutputSize += sizeof(UInt_t);
	  
	  newPad=kFALSE;
	}
	(*nBunches)++;

	outPtr[nAdded]=(UInt_t)currentTime;
	nAdded++;
	outputSize      += sizeof(UInt_t);
	blockOutputSize += sizeof(UInt_t);

	outPtr[nAdded]=(UInt_t)bunchSize;
	nAdded++;
	outputSize      += sizeof(UInt_t);
	blockOutputSize += sizeof(UInt_t);
	
	for(Int_t t=0;t<bunchSize;t++){
	  outPtr[nAdded]=(UInt_t)zeroSuppressedPad->GetDataSignal(currentTime+t);
	  nAdded++;
	  outputSize      += sizeof(UInt_t);
	  blockOutputSize += sizeof(UInt_t);
	}
      }
      zeroSuppressedPad->SetDataToDefault();
    }
  }

  if( wasInput ){

    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = outputSize-blockOutputSize;
    bd.fSize = blockOutputSize;
    bd.fSpecification = iter->fSpecification;
    Logging( kHLTLogDebug, "HLT::TPCZeroSuppressionComponent::DoEvent", "Event received", 
	     "Event 0x%08LX (%Lu) output data block %lu of %lu bytes at offset %lu",
	     evtData.fEventID, evtData.fEventID, ndx, blockOutputSize, outputSize-blockOutputSize );
    outputBlocks.push_back( bd );
  }
  else{
    size=0;
  }
  
  return 0;
}
