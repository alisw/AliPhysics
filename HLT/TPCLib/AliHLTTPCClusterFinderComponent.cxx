// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTPCClusterFinderComponent.cxx
    @author Timm Steinbeck, Matthias Richter, Jochen Thaeder, Kenneth Aamodt
    @date   
    @brief  The TPC cluster finder processing component
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCClusterFinderComponent.h"
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCDigitReaderRaw.h"
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

// this is a global object used for automatic component registration, do not use this
// use fPackedSwitch = true for packed inputtype "gkDDLPackedRawDataType"
// use fPackedSwitch = false for unpacked inputtype "gkUnpackedRawDataType"
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponentPacked(AliHLTTPCClusterFinderComponent::kClusterFinderPacked);
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponentUnpacked(AliHLTTPCClusterFinderComponent::kClusterFinderUnpacked);
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponentDecoder(AliHLTTPCClusterFinderComponent::kClusterFinderDecoder);

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterFinderComponent)

AliHLTTPCClusterFinderComponent::AliHLTTPCClusterFinderComponent(int mode)
  :
  fClusterFinder(NULL),
  fReader(NULL),
  fClusterDeconv(true),
  fXYClusterError(-1),
  fZClusterError(-1),
  fModeSwitch(mode),
  fUnsorted(0),
  fPatch(0),
  fPadArray(NULL),
  fGetActivePads(0)
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
  tgtList.push_back(AliHLTTPCDefinitions::fgkActivePadsDataType);
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

  Int_t rawreadermode =  -1;
  Int_t sigthresh = -1;
  Double_t sigmathresh= -1;
  Float_t occulimit = 1.0;
  Int_t oldRCUFormat=0;
  // Data Format version numbers:
  // 0: RCU Data format as delivered during TPC commissioning, pads/padrows are sorted, RCU trailer is one 32 bit word.
  // 1: As 0, but pads/padrows are delivered "as is", without sorting
  // 2: As 0, but RCU trailer is 3 32 bit words.
  // 3: As 1, but RCU trailer is 3 32 bit words.
  // -1: use offline raw reader

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- raw reader mode option
    if ( !strcmp( argv[i], "rawreadermode" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Missing rawreadermode", "Raw Reader Mode not specified" );
	return ENOTSUP;
      }

      // Decodes the rawreader mode: either number or string and returns the rawreadermode
      // -1 on failure, -2 for offline
      rawreadermode = AliHLTTPCDigitReaderRaw::DecodeMode( argv[i+1] );

      if (rawreadermode == -1 ) {
	Logging( kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Missing rawreadermode", "Cannot convert rawreadermode specifier '%s'.", argv[i+1] );
	return EINVAL;
      }

      i += 2;
      continue;
    }

    // -- pp run option
    if ( !strcmp( argv[i], "pp-run" ) ) {
      fClusterDeconv = false;
      i++;
      continue;
    }

    // -- zero suppression threshold
    if ( !strcmp( argv[i], "adc-threshold" ) ) {
      sigthresh = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert threshold specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[i], "occupancy-limit" ) ) {
      occulimit = strtod( argv[i+1], &cpErr);
      if ( *cpErr ) {
	HLTError("Cannot convert occupancy specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- number of timebins (default 1024)
    if ( !strcmp( argv[i], "timebins" ) ) {
      TString parameter(argv[i+1]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	AliHLTTPCTransform::SetNTimeBins(parameter.Atoi());
	HLTInfo("number of timebins set to %d, zbin=%f", AliHLTTPCTransform::GetNTimeBins(), AliHLTTPCTransform::GetZWidth());
      } else {
	HLTError("Cannot timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for rcu format
    if ( !strcmp( argv[i], "oldrcuformat" ) ) {
      oldRCUFormat = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert oldrcuformat specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
      
    // -- checking for unsorted clusterfinding
    if ( !strcmp( argv[i], "unsorted" ) ) {
      fUnsorted = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert unsorted specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
      
    // -- checking for active pads, used in 2007 December run
    if ( !strcmp( argv[i], "activepads" ) ) {
      fGetActivePads = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert activepads specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for nsigma-threshold, used in 2007 December run in ZeroSuppression
    if ( !strcmp( argv[i], "nsigma-threshold" ) ) {
       sigmathresh = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert nsigma-threshold specifier '%s'. Must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  // Choose reader
  if (fModeSwitch==kClusterFinderPacked) { 
    if (rawreadermode == -2) {
      HLTDebug("using AliHLTTPCDigitReaderPacked");
      fReader = new AliHLTTPCDigitReaderPacked();
      if(oldRCUFormat==1){
	fReader->SetOldRCUFormat(kTRUE);
      }
      else if(oldRCUFormat!=0){
	HLTWarning("Wrong oldrcuformat specifier %d; oldrcuformat set to default(kFALSE)",oldRCUFormat);
      }
      if(fUnsorted==1){
	fReader->SetUnsorted(kTRUE);
      }
      fClusterFinder->SetReader(fReader);
    } 
    else {
#if defined(HAVE_TPC_MAPPING)
      HLTDebug("using AliHLTTPCDigitReaderRaw mode %d", rawreadermode);
      fReader = new AliHLTTPCDigitReaderRaw(rawreadermode);
      fClusterFinder->SetReader(fReader);
#else //! defined(HAVE_TPC_MAPPING)
      HLTFatal("DigitReaderRaw not available - check your build");
      return -ENODEV;
#endif //defined(HAVE_TPC_MAPPING)
    }
  }
  else if(fModeSwitch==kClusterFinderUnpacked){
    HLTDebug("using AliHLTTPCDigitReaderUnpacked");
    fReader = new AliHLTTPCDigitReaderUnpacked();
    fClusterFinder->SetReader(fReader);
  }
  else if(fModeSwitch==kClusterFinderDecoder){
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
      
  // Variables to setup the Clusterfinder
  // TODO: this sounds strange and has to be verified; is the cluster finder not working when
  // fClusterDeconv = false ?
  fClusterDeconv = true;
  fXYClusterError = -1;
  fZClusterError = -1;

 
  fClusterFinder->SetDeconv( fClusterDeconv );
  fClusterFinder->SetXYError( fXYClusterError );
  fClusterFinder->SetZError( fZClusterError );
  if ( (fXYClusterError>0) && (fZClusterError>0) )
    fClusterFinder->SetCalcErr( false );
  fClusterFinder->SetSignalThreshold(sigthresh);
  fClusterFinder->SetNSigmaThreshold(sigmathresh);

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

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  //  == OUTdatatype pointer
  AliHLTTPCClusterData* outPtr;

  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, nSize, tSize = 0;

  outBPtr = outputPtr;
  outPtr = (AliHLTTPCClusterData*)outBPtr;

  Int_t slice, patch, row[2];
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
	HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
		 evtData.fEventID, evtData.fEventID, 
		 DataType2Text( iter->fDataType).c_str(), 
		 DataType2Text(AliHLTTPCDefinitions::fgkUnpackedRawDataType).c_str());

	if ( iter->fDataType != AliHLTTPCDefinitions::fgkUnpackedRawDataType ) continue;

      }
    	
      slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
      row[0] = AliHLTTPCTransform::GetFirstRow( patch );
      row[1] = AliHLTTPCTransform::GetLastRow( patch );


      if(fUnsorted){
	fClusterFinder->SetUnsorted(fUnsorted);
	fClusterFinder->SetPatch(patch);
      }
      /*      if(fUnsorted){
	      if(fPadArray==NULL){
	      fClusterFinder->SetUnsorted(fUnsorted);
	      //fPadArray = new AliHLTTPCPadArray(patch);
	      //fPadArray->InitializeVector(fModeSwitch);
	      }
	      else if(fPadArray->GetPatch()!=patch||fPadArray->GetPatch()==-1){
	      if (GetEventCount()<3) {
	      HLTWarning("pad array not initialized for data of specification 0x%08x, block skipped", iter->fSpecification);
	      } else if ((GetEventCount()%5000)==0) { // assuming 0.5 to 1kHz this gives a message rate of 0.1 to 0.5 Hz
	      HLTWarning("reminder: pad array not initialized for data of specification 0x%08x", iter->fSpecification);
	      }
	      continue;
	      }
	      }
      */
      

      outPtr = (AliHLTTPCClusterData*)outBPtr;

      maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);

      fClusterFinder->InitSlice( slice, patch, row[0], row[1], maxPoints );
      fClusterFinder->SetOutputArray( (AliHLTTPCSpacePointData*)outPtr->fSpacePoints );
	
      if(fUnsorted){
	fClusterFinder->ReadDataUnsorted(iter->fPtr, iter->fSize, fModeSwitch);

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
	       realPoints, slice, patch, row[0], row[1] );
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fSpecification = iter->fSpecification;
      bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
      //AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
      outputBlocks.push_back( bd );
	
      tSize += mysize;
      outBPtr += mysize;
      outPtr = (AliHLTTPCClusterData*)outBPtr;
	
      /*      if(fGetActivePads){
	      AliHLTTPCPadArray::AliHLTTPCActivePads* outPtrActive;
	      UInt_t activePadsSize, activePadsN = 0;
	      outPtrActive = (AliHLTTPCPadArray::AliHLTTPCActivePads*)outBPtr;
	      offset=tSize;
	      Int_t maxActivePads = (size-tSize)/sizeof(AliHLTTPCPadArray::AliHLTTPCActivePads);
	      activePadsSize= fClusterFinder->GetActivePads((AliHLTTPCPadArray::AliHLTTPCActivePads*)outPtrActive,maxActivePads)*sizeof(AliHLTTPCPadArray::AliHLTTPCActivePads);
 	
	      AliHLTComponentBlockData bdActive;
	      FillBlockData( bdActive );
	      bdActive.fOffset = offset;
	      bdActive.fSize = activePadsSize;
	      bdActive.fSpecification = iter->fSpecification;
	      bdActive.fDataType = AliHLTTPCDefinitions::fgkActivePadsDataType;
	      outputBlocks.push_back( bdActive );
 	
	      tSize+=activePadsSize;
	      outBPtr += activePadsSize;
	      outPtrActive = (AliHLTTPCPadArray::AliHLTTPCActivePads*)outBPtr;
	      }
      */
 

      if ( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::TPCClusterFinder::DoEvent", "Too much data", 
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		   tSize, size );
	  return EMSGSIZE;
	}
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
