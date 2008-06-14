// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCOfflineTrackerComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline tracker
*/

#include "AliHLTTPCOfflineTrackerComponent.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCtrackerMI.h"
#include "AliTPCClustersRow.h"
#include "AliESDEvent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliTracker.h"
#include "AliMagFMaps.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCOfflineTrackerComponent)

AliHLTTPCOfflineTrackerComponent::AliHLTTPCOfflineTrackerComponent() : AliHLTProcessor(),
fOutputPercentage(100),
fTPCGeomParam(0),
fTracker(0),
fESD(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCOfflineTrackerComponent::~AliHLTTPCOfflineTrackerComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCOfflineTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCOfflineTracker";
}

void AliHLTTPCOfflineTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // get input data type
  list.push_back(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/);
}

AliHLTComponentDataType AliHLTTPCOfflineTrackerComponent::GetOutputDataType()
{
  // create output data type
  return kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineTrackSegmentsDataType*/;
}

void AliHLTTPCOfflineTrackerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // get output data size
  constBase = 2000000;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTTPCOfflineTrackerComponent::Spawn()
{
  // create instance of the component
  return new AliHLTTPCOfflineTrackerComponent;
}

int AliHLTTPCOfflineTrackerComponent::DoInit( int argc, const char** argv )
{
  // init configuration 
  //
  int iResult=0;
#ifdef HAVE_NOT_TPC_LOAD_CLUSTERS
  HLTError("AliRoot version > v4-13-Release required");
  return -EFAULT;
#endif

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0 && !configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  // TPC geometry parameters
  fTPCGeomParam = new AliTPCParamSR;
  if (fTPCGeomParam) {
    fTPCGeomParam->ReadGeoMatrices();
  }

  // Init clusterer
  fTracker = new AliTPCtrackerMI(fTPCGeomParam);

  // AliESDEvent event needed by AliTPCtrackerMI
  // output of the component
  fESD = new AliESDEvent();
  if (fESD) {
    fESD->CreateStdContent();
  }

  // TODO: set the magnetic field correctly
  // the tracker needs the field map correctly initialized in AliTracker.
  // init from HLT/ConfigHLT/SolenoidBz or other appropriate CDB entry.
  // temporarily set to 5kG
  if (!AliTracker::GetFieldMap()) {
    // this instance must never be deleted, the AliRoot framework and the design
    // of AliTracker just does not support this. That's why we do not keep the
    // pointer. The memory leak is relativly small.
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
    AliTracker::SetFieldMap(field,kTRUE);
  }

  if (!fTracker || !fESD || !fTPCGeomParam) {
    HLTError("failed creating internal objects");
    iResult=-ENOMEM;
  }

  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::DoDeinit()
{
  // deinit configuration

  if(fTPCGeomParam) delete fTPCGeomParam; fTPCGeomParam = 0; 
  if(fTracker) delete fTracker; fTracker = 0; 
  if(fESD) delete fESD; fESD = 0;

  return 0;
}

int AliHLTTPCOfflineTrackerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // tracker function
  HLTInfo("DoEvent processing data");

//   Logging(kHLTLogDebug, "AliHLTTPCOfflineTrackerComponent::DoEvent", "Trigger data received",
//          "Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  TObjArray *clusterArray = 0;

  int iResult=0;

  if (fTracker && fESD) {
    // loop over input data blocks: TObjArrays of clusters
    for (TObject *pObj = (TObject *)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/,"TObjArray",0);
	 pObj !=0 && iResult>=0;
	 pObj = (TObject *)GetNextInputObject(0)) {
      clusterArray = dynamic_cast<TObjArray*>(pObj);
      if (!clusterArray) continue;
//       int lower=clusterArray->LowerBound();
//       int entries=clusterArray->GetEntries();
//       if (entries<=lower) continue;
//       if (clusterArray->At(lower)==NULL) continue; 
//       if (dynamic_cast<AliTPCClustersRow*>(clusterArray->At(lower))==NULL) continue;

      HLTInfo("load %d cluster rows from block %s 0x%08x", clusterArray->GetEntries(), DataType2Text(GetDataType(pObj)).c_str(), GetSpecification(pObj));
#ifndef HAVE_NOT_TPC_LOAD_CLUSTERS
      fTracker->LoadClusters(clusterArray);
#endif //HAVE_NOT_TPC_LOAD_CLUSTERS
    }// end loop over input objects

    // set magnetic field for the ESD, assumes correct initialization of
    // the field map
    fESD->SetMagneticField(AliTracker::GetBz());

    // run tracker
    fTracker->Clusters2Tracks(fESD);
    fTracker->UnloadClusters();

    Int_t nTracks = fESD->GetNumberOfTracks();
    HLTInfo("Number of tracks %d", nTracks);

    // TODO: calculate specification from the specification of input data blocks
    PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC, 0);

    // Alternatively: Push back tracks
//     for (Int_t it = 0; it < nTracks; it++) {
//       AliESDtrack* track = fESD->GetTrack(it);
//       PushBack(track, AliHLTTPCDefinitions::fgkOfflineTrackSegmentsDataType, 0);
//     }

    // is this necessary? If yes, we have to keep all the created TObjArrays
    // from the loop above
    // clear clusters
    //clusterArray->Clear();
    //clusterArray->Delete();

    // reset ESDs
    fESD->Reset();
  } else {
    HLTError("component not initialized");
    iResult=-ENOMEM;
  }

  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::Configure(const char* arguments)
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
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-something")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;

      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}
