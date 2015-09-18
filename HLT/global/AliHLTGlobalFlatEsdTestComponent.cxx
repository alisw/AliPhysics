// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

//  @file   AliHLTGlobalFlatEsdTestComponent.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Global ESD converter component.
// 

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTGlobalFlatEsdTestComponent.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDTrack.h"
#include "AliFlatExternalTrackParam.h"
#include "AliExternalTrackParam.h"

#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTrackMCLabel.h"
#include "AliHLTCTPData.h"
#include "AliHLTErrorGuard.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliPID.h"
#include "TTree.h"
#include "TList.h"
#include "TClonesArray.h"
//#include "AliHLTESDCaloClusterMaker.h"
//#include "AliHLTCaloClusterDataStruct.h"
//#include "AliHLTCaloClusterReader.h"
//#include "AliESDCaloCluster.h"
//#include "AliESDVZERO.h"
#include "AliHLTGlobalVertexerComponent.h"
#include "AliHLTVertexFinderBase.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterMCData.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalFlatEsdTestComponent)

AliHLTGlobalFlatEsdTestComponent::AliHLTGlobalFlatEsdTestComponent()
  : AliHLTProcessor()
  , fWriteClusters(0)
  , fVerbosity(0)  
  , fSolenoidBz(-5.00668)
  , fBenchmark("FlatEsdTest")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalFlatEsdTestComponent::~AliHLTGlobalFlatEsdTestComponent()
{
  // see header file for class documentation
}

int AliHLTGlobalFlatEsdTestComponent::Configure(const char* arguments)
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

int AliHLTGlobalFlatEsdTestComponent::Reconfigure(const char* cdbEntry, const char* chainId)
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

void AliHLTGlobalFlatEsdTestComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginOut);
}

AliHLTComponentDataType AliHLTGlobalFlatEsdTestComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut;
}

void AliHLTGlobalFlatEsdTestComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=2000000;
  inputMultiplier=10.0;
}

int AliHLTGlobalFlatEsdTestComponent::DoInit(int argc, const char** argv)
{
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

      // -noclusters
      if (argument.CompareTo("-noclusters")==0) {
	fWriteClusters=0;	
	// -clusters
      } else if (argument.CompareTo("-clusters")==0) {
	fWriteClusters=1;
      } else if (argument.Contains("-skipobject=")) {
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

  fSolenoidBz=GetBz();

  if (iResult>=0) {
    SetupCTPData();
  }

  fBenchmark.SetTimer(0,"total");

  return iResult;
}

int AliHLTGlobalFlatEsdTestComponent::DoDeinit()
{
  // see header file for class documentation

  return 0;
}

int AliHLTGlobalFlatEsdTestComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
						    const AliHLTComponentBlockData* /*blocks*/, 
						    AliHLTComponentTriggerData& /*trigData*/,
					            AliHLTUInt8_t* /*outputPtr*/, 
						    AliHLTUInt32_t& size,
						    AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  // see header file for class documentation
  int iResult=0;

  if (!IsDataEvent()) return iResult;

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  size = 0;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject | kAliHLTDataOriginOut); iter != NULL; iter = GetNextInputObject() ) {
    //fBenchmark.AddInput(pBlock->fSize);
    cout<<"Found ESD in esd test component !!!"<<endl;
    AliESDEvent *esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(iter));//dynamic_cast<const AliESDEvent*>(iter);
    if( esd ){
      esd->GetStdContent();
      cout<<"N ESD tracks: "<<esd->GetNumberOfTracks()<<endl;
    } else {
      cout<<"ESD pointer is NULL "<<endl;
    }
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject | kAliHLTDataOriginOut); iter != NULL; iter = GetNextInputObject() ) {     
    //fBenchmark.AddInput(pBlock->fSize);
    cout<<"Found ESD friend in esd test component !!!"<<endl;
    const AliESDfriend *esdFriend = dynamic_cast<const AliESDfriend*>(iter);
    if( esdFriend ){
      cout<<"N friend tracks: "<<esdFriend->GetNumberOfTracks()<<endl;
    } else {
      cout<<"ESD friend pointer is NULL "<<endl;
    }
  }

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
       pBlock!=NULL; pBlock=NULL/*GetNextInputBlock() there is only one block*/) {
    fBenchmark.AddInput(pBlock->fSize);
    AliFlatESDEvent *flatEvent = reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
    cout<<"N flat tracks: "<<flatEvent->GetNumberOfTracks()<<endl;
  }

  
  fBenchmark.Stop(0);
  HLTWarning( fBenchmark.GetStatistics() );

  return iResult;
}

