// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck, Matthias Richter                      *
//* Developers:      Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
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

/** @file   AliHLTTPCClusterFinderComponent.cxx
    @author Kenneth Aamodt <kenneth.aamodt@student.uib.no>
    @date   
    @brief  The TPC cluster finder processing component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCClusterFinderComponent.h"
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCDigitReaderDecoder.h"
#include "AliHLTTPCClusterFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include "AliHLTTPCDefinitions.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterFinderComponent)

AliHLTTPCClusterFinderComponent::AliHLTTPCClusterFinderComponent(int mode)
  :
  fClusterFinder(NULL),
  fReader(NULL),
  fDeconvTime(kFALSE),
  fDeconvPad(kFALSE),
  fClusterDeconv(false),
  fXYClusterError(-1),
  fZClusterError(-1),
  fModeSwitch(mode),
  fUnsorted(1),
  fPatch(0),
  fGetActivePads(0),
  fFirstTimeBin(-1),
  fLastTimeBin(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusterFinderComponent::~AliHLTTPCClusterFinderComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCClusterFinderComponent::GetComponentID()
{
  // see header file for class documentation
  switch(fModeSwitch){
  case kClusterFinderPacked:
    return "TPCClusterFinderPacked";
    break;
  case kClusterFinderUnpacked: 	 
    return "TPCClusterFinderUnpacked"; 	 
    break;
  case kClusterFinderDecoder:
    return "TPCClusterFinderDecoder";
    break;
  }
  HLTFatal("unknown digit reader type");
  return "";
}

void AliHLTTPCClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  switch(fModeSwitch){
  case kClusterFinderPacked:
    list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
    break;
  case kClusterFinderUnpacked: 	 
    list.push_back( AliHLTTPCDefinitions::fgkUnpackedRawDataType ); 	 
    break;
  case kClusterFinderDecoder:
    list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
    break;
  }
}

AliHLTComponentDataType AliHLTTPCClusterFinderComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTTPCClusterFinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  return tgtList.size();
}

void AliHLTTPCClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  switch(fModeSwitch){
  case 0:
    inputMultiplier = (6 * 0.4);
    break;
  case 1:
    inputMultiplier = 0.4;
    break;
  case 2:
    inputMultiplier = (6 * 0.4);
    break;
  }
}

AliHLTComponent* AliHLTTPCClusterFinderComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCClusterFinderComponent(fModeSwitch);
}
	
int AliHLTTPCClusterFinderComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  if ( fClusterFinder )
    return EINPROGRESS;

  fClusterFinder = new AliHLTTPCClusterFinder();

  Float_t occulimit = 1.0;

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- deconvolute-time option
    if ( !strcmp( argv[i], "-deconvolute-time" ) ) {
      fDeconvTime = kTRUE;
      i++;
      continue;
    }

    // -- deconvolute-pad option
    if ( !strcmp( argv[i], "-deconvolute-pad" ) ) {
      fDeconvPad = kTRUE;
      i++;
      continue;
    }

    // -- number of timebins (default 1024)
    if (!strcmp( argv[i], "-timebins") || !strcmp( argv[i], "timebins" )){
      TString parameter(argv[i+1]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	AliHLTTPCTransform::SetNTimeBins(parameter.Atoi());
	HLTInfo("number of timebins set to %d, zbin=%f", AliHLTTPCTransform::GetNTimeBins(), AliHLTTPCTransform::GetZWidth());
	fClusterFinder->UpdateLastTimeBin();
      } else {
	HLTError("Cannot timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      if(!strcmp( argv[i], "timebins")){
	HLTWarning("Argument 'timebins' is old, please switch to new argument naming convention (-timebins). The timebins argument will still work, but please change anyway.");
      }
      i+=2;
      continue;
    }

    // -first-timebin (default 0)
    if ( !strcmp( argv[i], "-first-timebin" ) ) {
      TString parameter(argv[i+1]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()){
	fFirstTimeBin=parameter.Atoi();
	HLTDebug("fFirstTimeBin set to %d",fFirstTimeBin);
      } 
      else {
	HLTError("Cannot -first-timebin specifier '%s'. Not a number.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -last-timebin (default 1024)
    if ( !strcmp( argv[i], "-last-timebin" ) ) {
      TString parameter(argv[i+1]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()){
	fLastTimeBin=parameter.Atoi();
	HLTDebug("fLastTimeBin set to %d",fLastTimeBin);
      } 
      else {
	HLTError("Cannot -last-timebin specifier '%s'. Not a number.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // --  unsorted option
    if ( !strcmp( argv[i], "-sorted" ) ) {
      fUnsorted=0;
      i++;
      continue;
    }

      
    // -- checking for active pads, used in 2007 December run
    if ( !strcmp( argv[i], "-active-pads" ) || !strcmp( argv[i], "activepads" ) ) {
      if(!strcmp( argv[i], "activepads" )){
	HLTWarning("Please change to new component argument naming scheme and use '-active-pads' instead of 'active-pads'");
      }
      fGetActivePads = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert activepads specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[i], "-occupancy-limit" ) || !strcmp( argv[i], "occupancy-limit" ) ) {
      if(!strcmp( argv[i], "occupancy-limit" )){
	HLTWarning("Please switch to new component argument naming convention, use '-occupancy-limit' instead of 'occupancy-limit'");
      }
      occulimit = strtod( argv[i+1], &cpErr);
      if ( *cpErr ) {
	HLTError("Cannot convert occupancy specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      if(fModeSwitch!=kClusterFinderPacked){
	HLTWarning("Argument '-occupancy-limit' is only used with -sorted set and with the TPCClusterFinderPacked , argument is deprecated");
      }
      i+=2;
      continue;
    }

    // -- raw reader mode option
    if ( !strcmp( argv[i], "rawreadermode" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Missing rawreadermode", "Raw Reader Mode not specified. rawreadermode is no longer a valid argument and will be deprecated even if rawreadermode is specified." );
	return ENOTSUP;
      }

      HLTWarning("Argument 'rawreadermode' is deprecated");      

      i += 2;
      continue;
    }

    // -- pp-run option
    if ( !strcmp( argv[i], "pp-run") ) {
      HLTWarning("Argument 'pp-run' is obsolete, deconvolution is swiched off in both time and pad directions by default.");
      fClusterDeconv = false;
      i++;
      continue;
    }

    // -- zero suppression threshold
    if ( !strcmp( argv[i], "adc-threshold" ) ) {
      strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert threshold specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      HLTWarning("'adc-threshold' is no longer a valid argument, please use TPCZeroSuppression component if you want to zerosuppress data.");
      i+=2;
      continue;
    }

    // -- checking for rcu format
    if ( !strcmp( argv[i], "oldrcuformat" ) ) {
      strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert oldrcuformat specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      HLTWarning("Argument 'oldrcuformat' is deprecated.");
      i+=2;
      continue;
    }
      
    // -- checking for unsorted clusterfinding (default 1)
    if ( !strcmp( argv[i], "unsorted" ) ) {
      fUnsorted = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert unsorted specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      HLTWarning("Argument 'unsorted' is old and does not follow the new argument naming convention. A change has been made, and the clusterfinder will read the data unsorted by default. For sorted reading, please use '-sorted' as argument. (unsorted 0 will do the same job, but please change anyway.)");
      i+=2;
      continue;
    }

    // -- checking for nsigma-threshold, used in 2007 December run in ZeroSuppression
    if ( !strcmp( argv[i], "nsigma-threshold" ) ) {
      strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert nsigma-threshold specifier '%s'. Must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      HLTWarning("Argument 'nsigma-threshold' argument is obsolete.");
      continue;
    }

    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  //Checking for conflicting arguments
  if(fClusterDeconv){
    if(fDeconvPad==kTRUE || fDeconvTime==kTRUE){
      HLTWarning("Conflicting arguments: argument 'pp-run' will be ignored.");
    }
  }
  if(occulimit!=1.0 && fUnsorted){
    HLTWarning("Argument 'occupancy-limit' is deprecated when doing unsorted data reading.");
  }
  if(fGetActivePads==kTRUE && fUnsorted==kFALSE){
    HLTWarning("Argument '-active-pads' only work with unsorted data reading. Active pads list will not be produced.");
  }
  

  // Choose reader
  if (fModeSwitch==kClusterFinderPacked) {
      HLTDebug("using AliHLTTPCDigitReaderPacked");
      fReader = new AliHLTTPCDigitReaderPacked();
      if(fUnsorted==1){	fReader->SetUnsorted(kTRUE); }
      fClusterFinder->SetReader(fReader);
  }
  else if(fModeSwitch==kClusterFinderUnpacked){ 	 
    HLTDebug("using AliHLTTPCDigitReaderUnpacked"); 	 
    fReader = new AliHLTTPCDigitReaderUnpacked(); 	 
    fClusterFinder->SetReader(fReader);
  }
  else if(fModeSwitch==kClusterFinderDecoder){
    HLTDebug("using AliHLTTPCDigitReaderDecoder");
    fReader = new AliHLTTPCDigitReaderDecoder();
    fClusterFinder->SetReader(fReader);
  }
  else{
    HLTFatal("No mode set for clusterfindercomponent");
  }
  // if pp-run use occupancy limit else set to 1. ==> use all 
  if ( !fClusterDeconv )
    fClusterFinder->SetOccupancyLimit(occulimit);
  else 
    fClusterFinder->SetOccupancyLimit(1.0);
      
  
  fClusterFinder->SetDeconv(fClusterDeconv);
  fClusterFinder->SetDeconvPad(fDeconvPad);
  fClusterFinder->SetDeconvTime(fDeconvPad);
  fClusterFinder->SetXYError( fXYClusterError );
  fClusterFinder->SetZError( fZClusterError );
  if ( (fXYClusterError>0) && (fZClusterError>0) ){
    fClusterFinder->SetCalcErr( false );
  }

  if(fFirstTimeBin>0){
    fClusterFinder->SetFirstTimeBin(fFirstTimeBin);
  }
  if(fLastTimeBin>0 && fLastTimeBin>fFirstTimeBin && fLastTimeBin<=AliHLTTPCTransform::GetNTimeBins()){
    fClusterFinder->SetLastTimeBin(fLastTimeBin);
  }

  return 0;
}

int AliHLTTPCClusterFinderComponent::DoDeinit()
{
  // see header file for class documentation

  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;
 
  if ( fReader )
    delete fReader;
  fReader = NULL;
    
  return 0;
}

int AliHLTTPCClusterFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, 
					      const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  if(fReader == NULL){
    HLTFatal("Digit reader not initialized, aborting event.");
    size=0;
    return 0;    
  }

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
    size=0;
    return 0;
  }

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  //  == OUTdatatype pointer
  AliHLTTPCClusterData* outPtr;

  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, nSize, tSize = 0;

  outBPtr = outputPtr;
  outPtr = (AliHLTTPCClusterData*)outBPtr;

  Int_t slice, patch;
  unsigned long maxPoints, realPoints = 0;

  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;


      if (fModeSwitch==0 || fModeSwitch==2) {
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
	     iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ) continue;

      }
      else if(fModeSwitch==1){
	HLTWarning("Unpacked data type depreciated, block aborted");
	HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
		 evtData.fEventID, evtData.fEventID, 
		 DataType2Text( iter->fDataType).c_str(), 
		 DataType2Text(AliHLTTPCDefinitions::fgkUnpackedRawDataType).c_str());
	continue;
      }

      slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

      if(fUnsorted){
	fClusterFinder->SetUnsorted(fUnsorted);
	fClusterFinder->SetPatch(patch);
      }

      outPtr = (AliHLTTPCClusterData*)outBPtr;

      maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);

      fClusterFinder->InitSlice( slice, patch, maxPoints );
      fClusterFinder->SetOutputArray( (AliHLTTPCSpacePointData*)outPtr->fSpacePoints );
	
      if(fUnsorted){
      	if(fGetActivePads){
	  fClusterFinder->SetDoPadSelection(kTRUE);
	}	
	if(fDeconvTime){
	  fClusterFinder->ReadDataUnsortedDeconvoluteTime(iter->fPtr, iter->fSize);
	}
	else{
	  fClusterFinder->ReadDataUnsorted(iter->fPtr, iter->fSize);
	}

	fClusterFinder->FindClusters();
      }
      else{
	fClusterFinder->Read(iter->fPtr, iter->fSize );
	fClusterFinder->ProcessDigits();
      }

      realPoints = fClusterFinder->GetNumberOfClusters();
	
      outPtr->fSpacePointCnt = realPoints;
      nSize = sizeof(AliHLTTPCSpacePointData)*realPoints;
      mysize += nSize+sizeof(AliHLTTPCClusterData);

      Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Spacepoints", 
	       "Number of spacepoints: %lu Slice/Patch/RowMin/RowMax: %d/%d/%d/%d.",
	       realPoints, slice, patch,AliHLTTPCTransform::GetFirstRow( patch ) , AliHLTTPCTransform::GetLastRow( patch ) );
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fSpecification = iter->fSpecification;
      bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
      outputBlocks.push_back( bd );
	
      tSize += mysize;
      outBPtr += mysize;
      outPtr = (AliHLTTPCClusterData*)outBPtr;
	

      if ( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::TPCClusterFinder::DoEvent", "Too much data", 
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		   tSize, size );
	  return -ENOSPC;
	}
	
      if(fUnsorted && fGetActivePads){
       Int_t maxNumberOfHW=(Int_t)((size-tSize)/sizeof(AliHLTUInt16_t)-1);
       AliHLTUInt16_t* outputHWPtr= (AliHLTUInt16_t*)(outputPtr+tSize);
       Int_t nHWAdd = fClusterFinder->FillHWAddressList(outputHWPtr, maxNumberOfHW);
      
       //cout<<"Number of hardwareaddresses: "<<nHWAdd<<endl;
       for(AliHLTUInt16_t test=0;test<nHWAdd;test++){
	 //cout<<"The HW address is: "<<(AliHLTUInt16_t)outputHWPtr[test]<<endl;
       }
       AliHLTComponentBlockData bdHW;
       FillBlockData( bdHW );
       bdHW.fOffset = tSize ;
       bdHW.fSize = nHWAdd*sizeof(AliHLTUInt16_t);
       bdHW.fSpecification = iter->fSpecification;
       bdHW.fDataType = kAliHLTDataTypeHwAddr16;
       outputBlocks.push_back( bdHW );
       
       tSize+=nHWAdd*sizeof(AliHLTUInt16_t);
      }
      fReader->Reset();
    }
    
  size = tSize;

  return 0;
}

int AliHLTTPCClusterFinderComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  const char* path="HLT/ConfigTPC";
  if (cdbEntry) path=cdbEntry;
  if (path) {
    HLTInfo("reconfigure from entry %s, chain id %s", path, (chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object: %s", pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  return 0;
}
