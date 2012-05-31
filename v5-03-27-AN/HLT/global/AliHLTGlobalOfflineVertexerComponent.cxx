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

//  @file   AliHLTGlobalOfflineVertexerComponent.cxx
//  @author Matthias Richter
//  @date   2010-04-19
//  @brief  Component wrapping the offline vertexer
//  @ingroup alihlt_global

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTGlobalOfflineVertexerComponent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliGRPRecoParam.h"
#include "AliLog.h"
#include "TMap.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalOfflineVertexerComponent)

AliHLTGlobalOfflineVertexerComponent::AliHLTGlobalOfflineVertexerComponent()
  : AliHLTProcessor()
  , fVertexer(NULL)
  , fBenchmark("GlobalOfflineVertexer")
{
  // The component subscribes to the HLT ESD and calculates the vertex
  // using the offline vertexer.
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

AliHLTGlobalOfflineVertexerComponent::~AliHLTGlobalOfflineVertexerComponent()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTGlobalOfflineVertexerComponent::GetComponentID()
{ 
  // component property: id
  return "GlobalOfflineVertexer";
}

void AliHLTGlobalOfflineVertexerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
    list.push_back(kAliHLTDataTypeESDObject);
}

AliHLTComponentDataType AliHLTGlobalOfflineVertexerComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTAnyDataType;
}

void AliHLTGlobalOfflineVertexerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 0;
  inputMultiplier = 1.0;
}

void AliHLTGlobalOfflineVertexerComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  // targetMap->Add(new TObjString("HLT/ConfigHLT/GlobalOfflineVertexer"),
  // 		new TObjString("configuration object"));
  targetMap->Add(new TObjString("GRP/GRP/Data"),
		new TObjString("GRP object"));
}

AliHLTComponent* AliHLTGlobalOfflineVertexerComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTGlobalOfflineVertexerComponent;
}

int AliHLTGlobalOfflineVertexerComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // init stage 1: default values for all data members

  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  TString cdbPath="HLT/ConfigHLT/";
  cdbPath+=GetComponentID();

  // TODO: activate later when the object has been copied to the
  // Grid OCDB
  //iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  if (iResult>=0) {
    // implement the component initialization
    fVertexer=new AliVertexerTracks;
    if (fVertexer) {
      fVertexer->SetFieldkG(GetBz());

      // initialize for TPC + ITS primary vertex
      fVertexer->SetITSMode();
      fVertexer->SetConstraintOff();
      // get cuts for vertexer from AliGRPRecoParam
      AliGRPRecoParam *grpRecoParam=NULL;
      TObject* obj=LoadAndExtractOCDBObject("GRP/Calib/RecoParam");
      if (obj) grpRecoParam=dynamic_cast<AliGRPRecoParam*>(obj);
      if (grpRecoParam) {
	Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
	Double_t *cutsVertexer = new Double_t[nCutsVertexer];
	grpRecoParam->GetVertexerTracksCutsITS(cutsVertexer,nCutsVertexer);
	fVertexer->SetCuts(cutsVertexer,nCutsVertexer);
	delete [] cutsVertexer; cutsVertexer = NULL;
      }
      // to run on HLT events we do not require ITS refit
      fVertexer->SetITSrefitNotRequired();
    }
  }

  if (iResult<0) {
    // implement cleanup
  }

  fBenchmark.SetTimer(0,"total");

  return iResult;
}

int AliHLTGlobalOfflineVertexerComponent::ScanConfigurationArgument(int /*argc*/, const char** /*argv*/)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  return 0;
}

int AliHLTGlobalOfflineVertexerComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fVertexer) delete fVertexer;
  fVertexer=NULL;

  return 0;
}

int AliHLTGlobalOfflineVertexerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  int iResult=0;

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!fVertexer) return -ENODEV;
  if (!IsDataEvent()) return 0;

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  for (const TObject* obj = GetFirstInputObject(kAliHLTDataTypeESDObject, "AliESDEvent");
       obj; obj = GetNextInputObject()) {

  // input objects are not supposed to be changed by the component, so they
  // are defined const. However, the implementation of AliESDEvent does not
  // support this and we need the const_cast
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  if (esd != NULL) {
    esd->GetStdContent();
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      AliESDtrack* track = esd->GetTrack(i);
      // TODO: this is just a quick hack, the vertexer requires kITSrefit and
      // at least 4 points
      // the HLT reconstruction needs to be adapted according to that requirement
      // we have to check which flag needs to be set at what stage
      track->SetStatus(AliESDtrack::kITSrefit);
    }
    AliESDVertex* vertex=fVertexer->FindPrimaryVertex(esd);
    if (vertex) {
      AliInfoClass("Offline Vertexer using the HLT ESD:");
      vertex->Print();
      iResult = PushBack( vertex, kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut);
    }
    break;
  }
  }

  fBenchmark.Stop(0);
  AliInfoClass( fBenchmark.GetStatistics() );

  return iResult;
}

int AliHLTGlobalOfflineVertexerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  HLTInfo("reconfigure '%s' from entry %s", chainId, cdbEntry);

  return 0;
}

int AliHLTGlobalOfflineVertexerComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
