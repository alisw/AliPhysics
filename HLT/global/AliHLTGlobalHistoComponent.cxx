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

/// @file   AliHLTGlobalHistoComponent.cxx
/// @author Matthias Richter
/// @date   2010-09-16
/// @brief  A histogramming component for global ESD properties based
///         on the AliHLTTTreeProcessor

#include "AliHLTGlobalHistoComponent.h"
#include "AliESDEvent.h"
#include "TTree.h"
#include "TString.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalHistoComponent)

AliHLTGlobalHistoComponent::AliHLTGlobalHistoComponent()
  : AliHLTTTreeProcessor()
  , fEvent(0)
  , fNofTracks(0)
  , fNofV0s(0)
  , fNofContributors(0)
  , fVertexX(-99)
  , fVertexY(-99)
  , fVertexZ(-99)
  , fVertexStatus(kFALSE)
  , fTrackVariables()
  , fTrackVariablesInt()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTGlobalHistoComponent::~AliHLTGlobalHistoComponent()
{
  // see header file for class documentation
  fTrackVariables.Reset();
  fTrackVariablesInt.Reset();
}

void AliHLTGlobalHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTAllDataTypes);
}

TTree* AliHLTGlobalHistoComponent::CreateTree(int /*argc*/, const char** /*argv*/)
{
  // create the tree and branches
  int iResult=0;
  TTree* pTree=new TTree("ESDproperties", "HLT ESD properties");
  if (!pTree) return NULL;

  const char* trackVariableNames = {
    // Note the black at the end of each name!
    "Track_pt "
    "Track_phi "
    "Track_eta "
    "Track_p "
    "Track_theta "
    "Track_Nclusters "
    "Track_status "
    "Track_charge "
    "Track_DCAr "
    "Track_DCAz "
    "Track_dEdx "
  };
  
  const char* trackIntVariableNames = {
    // Note the black at the end of each name!
    "Track_status "
  };
  
  int maxTrackCount=20000; // FIXME: make configurable
  
  if ((iResult=fTrackVariables.Init(maxTrackCount, trackVariableNames))<0) {
    HLTError("failed to initialize internal structure for track properties (float)");
  }
  if ((iResult=fTrackVariablesInt.Init(maxTrackCount, trackIntVariableNames))<0) {
    HLTError("failed to initialize internal structure for track properties (int)");
  }
  
  if (iResult>=0) {
    pTree->Branch("event",        &fEvent,           "event/I");
    pTree->Branch("trackcount",   &fNofTracks,       "trackcount/I");
    pTree->Branch("vertexX",      &fVertexX,         "vertexX/F");
    pTree->Branch("vertexY",      &fVertexY,         "vertexY/F");
    pTree->Branch("vertexZ",      &fVertexZ,         "vertexZ/F");
    pTree->Branch("nV0",          &fNofV0s,          "nV0/I");
    pTree->Branch("nContributors",&fNofContributors, "nContributors/I");
    pTree->Branch("vertexStatus", &fVertexStatus,    "vertexStatus/I");

    int i=0;
    // FIXME: this is a bit ugly since type 'f' and 'i' are specified
    // explicitely. Would be better to use a function like
    // AliHLTGlobalHistoVariables::GetType but could not get this working
    for (i=0; i<fTrackVariables.Variables(); i++) {
      TString specifier=fTrackVariables.GetKey(i);
      float* pArray=fTrackVariables.GetArray(specifier);
      specifier+="[trackcount]/f";
      pTree->Branch(fTrackVariables.GetKey(i), pArray, specifier.Data());
    }
    for (i=0; i<fTrackVariablesInt.Variables(); i++) {
      TString specifier=fTrackVariablesInt.GetKey(i);
      int* pArray=fTrackVariablesInt.GetArray(specifier);
      specifier+="[trackcount]/i";
      pTree->Branch(fTrackVariablesInt.GetKey(i), pArray, specifier.Data());
    }
  } else {
    delete pTree;
    pTree=NULL;
  }

  return pTree;
}

void AliHLTGlobalHistoComponent::FillHistogramDefinitions()
{
  /// default histogram definitions
}

int AliHLTGlobalHistoComponent::FillTree(TTree* pTree, const AliHLTComponentEventData& /*evtData*/, 
					 AliHLTComponentTriggerData& /*trigData*/ )
{
  /// fill the tree from the ESD
  int iResult=0;
  if (!IsDataEvent()) return 0;

  ResetVariables();

  // fetch ESD from input stream
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  esd->GetStdContent();

  // fill track variables
  fNofTracks       = esd->GetNumberOfTracks();
  fVertexX         = esd->GetPrimaryVertexTracks()->GetX();
  fVertexY         = esd->GetPrimaryVertexTracks()->GetY();
  fVertexZ         = esd->GetPrimaryVertexTracks()->GetZ();
  fNofV0s          = esd->GetNumberOfV0s();
  fNofContributors = esd->GetPrimaryVertexTracks()->GetNContributors();
  fVertexStatus    = esd->GetPrimaryVertexTracks()->GetStatus();
  
  for (int i=0; i<fNofTracks; i++) {
    AliESDtrack *esdTrack = esd->GetTrack(i);
    if (!esdTrack) continue;
    
    Float_t DCAr, DCAz = -99;
    esdTrack->GetImpactParametersTPC(DCAr, DCAz);
    
    fTrackVariables.Fill("Track_pt"        , esdTrack->Pt()                      );
    fTrackVariables.Fill("Track_phi"       , esdTrack->Phi()*TMath::RadToDeg()   );
    fTrackVariables.Fill("Track_eta"       , esdTrack->Theta()                   );
    fTrackVariables.Fill("Track_p"         , esdTrack->P()                       );
    fTrackVariables.Fill("Track_theta"     , esdTrack->Theta()*TMath::RadToDeg() );
    fTrackVariables.Fill("Track_Nclusters" , esdTrack->GetTPCNcls()              );
    fTrackVariables.Fill("Track_status"    , esdTrack->GetStatus()               );
    fTrackVariables.Fill("Track_charge"    , esdTrack->Charge()                  );
    fTrackVariables.Fill("Track_DCAr"      , DCAr                        	 );
    fTrackVariables.Fill("Track_DCAz"      , DCAz                                );   
    fTrackVariables.Fill("Track_dEdx"      , esdTrack->GetTPCsignal()            );   

    fTrackVariablesInt.Fill("Track_status" , esdTrack->GetStatus()               );   
   }
  //HLTInfo("added parameters for %d tracks", fNofTracks);

  if (iResult<0) {
    // fill an empty event
    ResetVariables();
  }
  fEvent++;

  pTree->Fill();
  return iResult;
}

int AliHLTGlobalHistoComponent::ResetVariables()
{
  /// reset all filling variables
  fNofTracks=0;
  fTrackVariables.ResetCount();
  fTrackVariablesInt.ResetCount();
  return 0;
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOriginDataType() const
{
  // get the origin of the output data
  return kAliHLTVoidDataType;
}

// AliHLTUInt32_t AliHLTGlobalHistoComponent::GetDataSpec() const;
// {
//   // get specifications of the output data
//   return 0;
// }
