//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: S.Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>    *
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

/** @file   AliHLTdNdPtAnalysisMergerComponent.cxx
    @author Sergey Gorbunov
    @brief  Component for ploting charge in clusters
*/

#include "AliHLTdNdPtAnalysisMergerComponent.h"
#include "AlidNdPtAnalysisPbPb.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliHLTMessage.h"
#include "TTimeStamp.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

using namespace std;

AliHLTdNdPtAnalysisMergerComponent::AliHLTdNdPtAnalysisMergerComponent() :
  fUID(0),
  fBenchmark("dNdPtMerger"),
  fAnalysis(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTdNdPtAnalysisMergerComponent::~AliHLTdNdPtAnalysisMergerComponent()
{
  // see header file for class documentation

  if( fAnalysis ) delete fAnalysis;
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTdNdPtAnalysisMergerComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "dNdPtAnalysisMerger";
}

void AliHLTdNdPtAnalysisMergerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypedNdPt|kAliHLTDataOriginAny );
}

AliHLTComponentDataType AliHLTdNdPtAnalysisMergerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypedNdPt  | kAliHLTDataOriginOut;
}

void AliHLTdNdPtAnalysisMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTdNdPtAnalysisMergerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTdNdPtAnalysisMergerComponent;
}

void AliHLTdNdPtAnalysisMergerComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA merger component
  // Some parameters can be later overwritten from the OCDB

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");  
}



int AliHLTdNdPtAnalysisMergerComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA merger component from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;
    HLTError( "Unknown option \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }
  delete pTokens;
  
  if ( bMissingParam ) {
    HLTError( "Specifier missed for parameter \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }

  return iResult;
}

int AliHLTdNdPtAnalysisMergerComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
{
  // see header file for class documentation
  
  const char* defaultNotify = "";
  
  if ( !cdbEntry ) {
    return 0;
    cdbEntry = "HLT/ConfigAnalysis/dNdPtAnalysisMerger";
    defaultNotify = " (default)";
    chainId = 0;
  }

  HLTInfo( "configure from entry \"%s\"%s, chain id %s", cdbEntry, defaultNotify, ( chainId != NULL && chainId[0] != 0 ) ? chainId : "<none>" );
  AliCDBEntry *pEntry = AliCDBManager::Instance()->Get( cdbEntry );//,GetRunNo());

  if ( !pEntry ) {
    HLTError( "cannot fetch object \"%s\" from CDB", cdbEntry );
    return -EINVAL;
  }

  TObjString* pString = dynamic_cast<TObjString*>( pEntry->GetObject() );

  if ( !pString ) {
    HLTError( "configuration object \"%s\" has wrong type, required TObjString", cdbEntry );
    return -EINVAL;
  }

  HLTInfo( "received configuration object string: \"%s\"", pString->GetString().Data() );

  return  ReadConfigurationString( pString->GetString().Data() );
}

int AliHLTdNdPtAnalysisMergerComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step


  //* read hard-coded values
  
  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}

int AliHLTdNdPtAnalysisMergerComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB

  return Configure( cdbEntry, chainId, NULL );
}


int AliHLTdNdPtAnalysisMergerComponent::DoInit( int argc, const char** argv ) 
{
  // init

  fUID = 0;
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");

  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  int iResult = Configure( NULL, NULL, arguments.Data()  );

  return iResult; 
}
  
int AliHLTdNdPtAnalysisMergerComponent::DoDeinit()
{
  // see header file for class documentation
  fUID = 0;
  if( fAnalysis ) delete fAnalysis;
  fAnalysis = 0;
  return 0;
}


int AliHLTdNdPtAnalysisMergerComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) ) {
    return 0;
  }

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }

  TList *list = new TList();

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypedNdPt); iter != NULL; iter = GetNextInputObject() ) {

    AlidNdPtAnalysisPbPb *instance = dynamic_cast<AlidNdPtAnalysisPbPb*>(const_cast<TObject*>( iter ) );
    if( !instance ){ 
      HLTWarning("Wrong AlidNdPtAnalysisPbPb object received");
      continue;
    }
    list->Add(instance); 
  }

  fBenchmark.Start(1);
  if( !list->IsEmpty() ){
    delete fAnalysis;
    fAnalysis = new AlidNdPtAnalysisPbPb;
    if( !fAnalysis ){
      HLTError("Out of memory");
      return -ENOMEM;
    }
    fAnalysis->Init();  
    fAnalysis->Merge(list);
  }
  fBenchmark.Stop(1);

  if( fAnalysis ) PushBack( (TObject*) fAnalysis, kAliHLTDataTypedNdPt,fUID);
  delete list;

  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
  return 0;
}


