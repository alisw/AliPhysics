//-*- Mode: C++ -*-
// $Id: AliHLTGlobalEsdToFlatConverterComponent.cxx $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Steffen Weber                                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTGlobalEsdToFlatConverterComponent.cxx
    @author  Steffen Weber
    @brief   Component to convert ESD objects to flatESD objects
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "TList.h"
#include "AliESDEvent.h"
#include "AliFlatESDEvent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTGlobalEsdToFlatConverterComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "TTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliSysInfo.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalEsdToFlatConverterComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTGlobalEsdToFlatConverterComponent::AliHLTGlobalEsdToFlatConverterComponent() :
  AliHLTProcessor()
  {
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTGlobalEsdToFlatConverterComponent::~AliHLTGlobalEsdToFlatConverterComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTGlobalEsdToFlatConverterComponent::GetComponentID() { 
  // see header file for class documentation
  return "GlobalEsdToFlatConverter";
}

// #################################################################################
void AliHLTGlobalEsdToFlatConverterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
	list.push_back(kAliHLTDataTypeESDTree|kAliHLTDataOriginOut );
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut );
  list.push_back( kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginOut );
}


AliHLTComponentDataType AliHLTGlobalEsdToFlatConverterComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTGlobalEsdToFlatConverterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){ 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut );
  tgtList.push_back( kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut );
  return tgtList.size();
}

// #################################################################################
void AliHLTGlobalEsdToFlatConverterComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 10000;
  inputMultiplier = 10.0;
}


// #################################################################################
AliHLTComponent* AliHLTGlobalEsdToFlatConverterComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTGlobalEsdToFlatConverterComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTGlobalEsdToFlatConverterComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
  printf("AliHLTGlobalEsdToFlatConverterComponent::DoInit\n");
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;

  // default list of skiped ESD objects
  TString skipObjects=
    // "AliESDRun,"
    // "AliESDHeader,"
    // "AliESDZDC,"
    "AliESDFMD,"
    // "AliESDVZERO,"
    // "AliESDTZERO,"
    // "TPCVertex,"
    // "SPDVertex,"
    // "PrimaryVertex,"
    // "AliMultiplicity,"
    // "PHOSTrigger,"
    // "EMCALTrigger,"
    // "SPDPileupVertices,"
    // "TrkPileupVertices,"
    "Cascades,"
    "Kinks,"
    "AliRawDataErrorLogs,"
    "AliESDACORDE";

  iResult=Reconfigure(NULL, NULL);
  TString allArgs = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !allArgs.IsNull() ) allArgs += " ";
    allArgs += argv[i];
  }

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->String();	
      if (argument.IsNull()) continue;
 if (argument.Contains("-skipobject=")) {
	argument.ReplaceAll("-skipobject=", "");
	skipObjects=argument;
      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }


  if (iResult>=0) {
    SetupCTPData();
  }

  return iResult;
}



// #################################################################################
Int_t AliHLTGlobalEsdToFlatConverterComponent::DoDeinit() {
  // see header file for class documentation

	
  return 0;
}

// #################################################################################
Int_t AliHLTGlobalEsdToFlatConverterComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
						    const AliHLTComponentBlockData* /*blocks*/, 
						    AliHLTComponentTriggerData& /*trigData*/,
						    AliHLTUInt8_t* outputPtr, 
						    AliHLTUInt32_t& size,
						    AliHLTComponentBlockDataList& outputBlocks) {
  // see header file for class documentation

	
	
	
  AliSysInfo::AddStamp("AliHLTGlobalEsdToFlatConverterComponent::DoEvent.Start");
  Int_t iResult=0;

  Bool_t err = kFALSE;
	
	
  size_t maxOutputSize = size;
  size = 0;
	
  // -- Only use data event
  if( !IsDataEvent() ) return 0;

  AliESDEvent *esd =0;
  AliESDfriend *esdFriend=0;
	
	
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject | kAliHLTDataOriginOut); iter != NULL; iter = GetNextInputObject() ) {
    esd =dynamic_cast<AliESDEvent*>(const_cast<TObject*>(iter));
    if( esd ){
      esd->GetStdContent();
      iResult=1;
    } else {
    }
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject | kAliHLTDataOriginOut); iter != NULL; iter = GetNextInputObject() ) {     
    esdFriend = dynamic_cast<AliESDfriend*>(const_cast<TObject*>(iter));
    if( esdFriend ){
    } else {
    }
  }
  AliFlatESDEvent *flatEsd = reinterpret_cast<AliFlatESDEvent*>(outputPtr);
  new (flatEsd) AliFlatESDEvent;
  
  Int_t converted1 = flatEsd->SetFromESD(AliFlatESDEvent::EstimateSize(esd),esd, kTRUE); 
  
  if( converted1>=0 && maxOutputSize > flatEsd->GetSize() ){
    AliHLTComponentBlockData outBlock;
    FillBlockData( outBlock );
    outBlock.fOffset = size;
    outBlock.fSize = flatEsd->GetSize();
    outBlock.fDataType = kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut;
    outputBlocks.push_back( outBlock );
    size += outBlock.fSize;
  }	
  else {    
    return 0;	
  }

  if(esdFriend){
    AliFlatESDFriend *flatFriend = reinterpret_cast<AliFlatESDFriend*>(outputPtr + size); 
    size_t freeSpaceTotal = maxOutputSize - size;
    size_t freeSpace = freeSpaceTotal;

    err = ( freeSpace < sizeof( AliFlatESDEvent ) );    
    if( err ) return 0;
    
    new (flatFriend) AliFlatESDFriend;
    freeSpace = freeSpaceTotal - flatFriend->GetSize();
    
    Int_t converted = flatFriend->SetFromESDfriend(AliFlatESDFriend::EstimateSize(esdFriend),esdFriend); 
    
    Int_t insizeEvent, insizeFriend, outsizeEvent, outsizeFriend =  0 ;
    
    insizeEvent = AliSysInfo::EstimateObjectSize(esd);
    
    insizeFriend = AliSysInfo::EstimateObjectSize(esdFriend);
    
    outsizeEvent = (int) flatEsd->GetSize();
    outsizeFriend = (int) flatFriend->GetSize();
    AliSysInfo::AddStamp("AliHLTGlobalEsdToFlatConverterComponent::DoEvent.Stop", insizeEvent, outsizeEvent, outsizeFriend, insizeFriend);
    
    if( converted>=0 ){ // set up the output block description
      
      AliHLTComponentBlockData outBlock;
      FillBlockData( outBlock );
      outBlock.fOffset = size;
      outBlock.fSize = flatFriend->GetSize();
      outBlock.fDataType = kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut;
      outputBlocks.push_back( outBlock );
      size += outBlock.fSize;
    }
  }
  return iResult;
}


// #################################################################################
Int_t AliHLTGlobalEsdToFlatConverterComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}


int AliHLTGlobalEsdToFlatConverterComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->String();	
      if (argument.IsNull()) continue;      
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
      break;
    }  
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  return iResult;
}

int AliHLTGlobalEsdToFlatConverterComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=NULL;
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->String().Data());
	iResult=Configure(pString->String().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;
}
