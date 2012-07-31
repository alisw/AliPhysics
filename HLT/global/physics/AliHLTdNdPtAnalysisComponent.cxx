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

/** @file   AliHLTdNdPtAnalysisComponent.cxx
    @author Sergey Gorbunov
    @brief  Component for ploting charge in clusters
*/

#include "AliHLTdNdPtAnalysisComponent.h"
#include "AlidNdPtAnalysisPbPb.h"
#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AliESDtrackCuts.h"

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

AliHLTdNdPtAnalysisComponent::AliHLTdNdPtAnalysisComponent() :
  fUID(0),
  fBenchmark("dNdPt"),
  fAnalysis(0),
  fEventCuts(0),
  fAcceptanceCuts(0),
  fESDtrackCuts(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTdNdPtAnalysisComponent::~AliHLTdNdPtAnalysisComponent()
{
  // see header file for class documentation
  if( fAnalysis ){
    fAnalysis->SetEventCuts(0);
    fAnalysis->SetAcceptanceCuts(0);
    fAnalysis->SetTrackCuts(0);
  }
  delete fAnalysis;
  delete fEventCuts;
  delete fAcceptanceCuts;
  delete fESDtrackCuts;
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTdNdPtAnalysisComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "dNdPtAnalysis";
}

void AliHLTdNdPtAnalysisComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut );
}

AliHLTComponentDataType AliHLTdNdPtAnalysisComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypedNdPt  | kAliHLTDataOriginOut;
}

void AliHLTdNdPtAnalysisComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTdNdPtAnalysisComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTdNdPtAnalysisComponent;
}

void AliHLTdNdPtAnalysisComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA merger component
  // Some parameters can be later overwritten from the OCDB

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");    

  if( !fEventCuts || !fAcceptanceCuts || !fESDtrackCuts){
    HLTError("Out of memory");
    return ;
  }
  
  fEventCuts->SetTriggerRequired(kFALSE);

  fESDtrackCuts->SetRequireTPCRefit(kFALSE);
  fESDtrackCuts->SetAcceptKinkDaughters(kTRUE);
  fESDtrackCuts->DefineHistograms();

  ReadConfigurationString("-vertexZRange 20. -meanVertexXYZ 0. 0. 0. -meanVertexXYZSigma 1. 1. 10. -etaRange 0.9 -ptRange 0.15 1.e10 -maxDCAr 3.0 -maxDCAz 30.0 -maxDCAToVertexXY 3.0 -maxDCAToVertexZ 3.0 -requireSigmaToVertex 1");

}



int AliHLTdNdPtAnalysisComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA merger component from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  if( !fEventCuts || !fAcceptanceCuts || !fESDtrackCuts){
    HLTError("Out of memory");
    return -ENOMEM;
  }

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;

    if ( argument.CompareTo( "-vertexZRange" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float z = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fEventCuts->SetZvRange(-z, z );
      HLTInfo( "Z vertex range is set to [%f,%f]", -z,z);
      continue;
    }
 
    if ( argument.CompareTo( "-meanVertexXYZ" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float y = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float z = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fEventCuts->SetMeanXYZv(x,y,z);
      HLTInfo( "mean vertex is set to (%f,%f,%f)", x,y,z);
      continue;
    }
 
    if ( argument.CompareTo( "-meanVertexXYZSigma" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float y = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float z = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fEventCuts->SetSigmaMeanXYZv(x,y,z);
      HLTInfo( "mean vertex sigma is set to (%f,%f,%f)", x,y,z);
      continue;
    }

    if ( argument.CompareTo( "-etaRange" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float eta = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fAcceptanceCuts->SetEtaRange(-eta,eta);     
      HLTInfo( "eta range is set to [%f,%f]", -eta,eta);
      continue;
    }

    if ( argument.CompareTo( "-ptRange" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float ptMin = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float ptMax = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fAcceptanceCuts->SetPtRange(ptMin,ptMax);
      HLTInfo( "pt range is set to [%f,%f]", -ptMin,ptMax);
      continue;
    }
 
    if ( argument.CompareTo( "-maxDCAr" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fAcceptanceCuts->SetMaxDCAr(x);
      HLTInfo( "max DCA R is set to %f", x);
      continue;
    }
 
    if ( argument.CompareTo( "-maxDCAz" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fAcceptanceCuts->SetMaxDCAz(x);
      HLTInfo( "max DCA Z is set to %f", x);
      continue;
    }

    if ( argument.CompareTo( "-maxDCAToVertexXY" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fESDtrackCuts->SetMaxDCAToVertexXY(x);    
      HLTInfo( "max DCA to XY vertex  is set to %f", x);
      continue;
    }

    if ( argument.CompareTo( "-maxDCAToVertexZ" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      float x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      fESDtrackCuts->SetMaxDCAToVertexZ(x);    
      HLTInfo( "max DCA Z to vertex  is set to %f", x);
      continue;
    }

    if ( argument.CompareTo( "-requireSigmaToVertex" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      bool x = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
      fESDtrackCuts->SetRequireSigmaToVertex(x);
      HLTInfo( "sigma to vertex requirement is set to %i",x);
      continue;
    }

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

int AliHLTdNdPtAnalysisComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
{
  // see header file for class documentation

  const char* defaultNotify = "";
  
  if ( !cdbEntry ) {
    cdbEntry = "HLT/ConfigAnalysis/dNdPtAnalysis"; 
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

int AliHLTdNdPtAnalysisComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  if( !fEventCuts ) fEventCuts = new AlidNdPtEventCuts;
  if( !fAcceptanceCuts ) fAcceptanceCuts = new AlidNdPtAcceptanceCuts;
  if( !fESDtrackCuts ) fESDtrackCuts = new AliESDtrackCuts;

  if( !fEventCuts || !fAcceptanceCuts || !fESDtrackCuts){
    HLTError("Out of memory");
    return -ENOMEM;
  }

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

int AliHLTdNdPtAnalysisComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB

  return Configure( cdbEntry, chainId, NULL );
}


int AliHLTdNdPtAnalysisComponent::DoInit( int argc, const char** argv ) 
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

  if( fAnalysis ){
    fAnalysis->SetEventCuts(0);
    fAnalysis->SetAcceptanceCuts(0);
    fAnalysis->SetTrackCuts(0);
  }
  delete fAnalysis;  

  fAnalysis = new AlidNdPtAnalysisPbPb;

  if( !fAnalysis || !fEventCuts || !fAcceptanceCuts || !fESDtrackCuts){
    HLTError("Out of memory");
    return -ENOMEM;
  }

  fAnalysis->SetUseMCInfo(kFALSE);
  fAnalysis->SetAnalysisMode(AlidNdPtHelper::kTPC);
  fAnalysis->SetTrigger(AliTriggerAnalysis::kMB1);
  fAnalysis->SetTriggerClass(0) ;
  fAnalysis->SetParticleMode(AlidNdPtHelper::kAllPart);
  fAnalysis->SetPhysicsTriggerSelection(0);

  fAnalysis->SetEventCuts(fEventCuts);
  fAnalysis->SetAcceptanceCuts(fAcceptanceCuts);
  fAnalysis->SetTrackCuts(fESDtrackCuts);

  fAnalysis->SetBackgroundCuts(0);
  
  fAnalysis->Init();

  return iResult; 
}
  
int AliHLTdNdPtAnalysisComponent::DoDeinit()
{
  // see header file for class documentation
  fUID = 0;
  if( fAnalysis ){
    fAnalysis->SetEventCuts(0);
    fAnalysis->SetAcceptanceCuts(0);
    fAnalysis->SetTrackCuts(0);
  }
  delete fAnalysis;
  fAnalysis = 0;
  return 0;
}


int AliHLTdNdPtAnalysisComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) ) {
    return 0;
  }

  if( !fAnalysis ){
    HLTWarning("Analysis objects has not been created!!!");
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

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {

    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !esdEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      continue;
    }
    esdEvent->GetStdContent();
    fBenchmark.Start(1);
    fAnalysis->Process(esdEvent);
    fBenchmark.Stop(1);
  }

  if( fAnalysis ) PushBack( (TObject*) fAnalysis, kAliHLTDataTypedNdPt,fUID);

  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
  return 0;
}


