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
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

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
  , fMaxTrackCount(20000)
  , fMaxV0Count(1000)
  , fFillV0(kFALSE)
  , fTrackVariables()
  , fTrackVariablesInt()
  , fV0Variables()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalHistoComponent::~AliHLTGlobalHistoComponent(){
  // see header file for class documentation
  fTrackVariables.Reset();
  fTrackVariablesInt.Reset();
  fV0Variables.Reset();
}

void AliHLTGlobalHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list){
  // see header file for class documentation
  list.push_back(kAliHLTAllDataTypes);
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOutputDataType(){
  // see header file for class documentation
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginOut;
}

TTree* AliHLTGlobalHistoComponent::CreateTree(int /*argc*/, const char** /*argv*/){
// create the tree and branches
 
  int iResult=0;
  TTree* pTree = new TTree("ESDproperties", "HLT ESD properties");
  if (!pTree) return NULL;

  const char* trackVariableNames = {
    // Note the black at the end of each name!
    "Track_pt "
    "Track_phi "
    "Track_eta "
    "Track_p "
    "Track_theta "
    "Track_TPCclus "
    "Track_ITSclus "
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

  const char* V0VariableNames = {
    // Note the black at the end of each name!
    "V0_AP "
    "V0_pt "  
    "clust1 "
    "clust2 "
    "dev1 "
    "dev2 "
    "devPrim "
    "length "
    "sigmaLength "
    "r "
  };
     
  if((iResult=fTrackVariables.Init(fMaxTrackCount, trackVariableNames))<0){
    HLTError("failed to initialize internal structure for track properties (float)");
  }
  if((iResult=fTrackVariablesInt.Init(fMaxTrackCount, trackIntVariableNames))<0){
    HLTError("failed to initialize internal structure for track properties (int)");
  }
  if((iResult=fV0Variables.Init(fMaxV0Count, V0VariableNames))<0){
      HLTError("failed to initialize internal structure for V0 properties (float)");
  }
  
  if(iResult>=0){     
     pTree->Branch("event",	   &fEvent,	      "event/I");
     pTree->Branch("trackcount",   &fNofTracks,       "trackcount/I");
     pTree->Branch("vertexX",	   &fVertexX,	      "vertexX/F");
     pTree->Branch("vertexY",	   &fVertexY,	      "vertexY/F");
     pTree->Branch("vertexZ",	   &fVertexZ,	      "vertexZ/F");
     pTree->Branch("nContributors",&fNofContributors, "nContributors/I");
     pTree->Branch("vertexStatus", &fVertexStatus,    "vertexStatus/I");
     if(fFillV0==kTRUE) pTree->Branch("V0", &fNofV0s, "V0/I");

     int i=0;
     // FIXME: this is a bit ugly since type 'f' and 'i' are specified
     // explicitely. Would be better to use a function like
     // AliHLTGlobalHistoVariables::GetType but could not get this working
    
     for(i=0; i<fTrackVariables.Variables(); i++){
         TString specifier=fTrackVariables.GetKey(i);
         float* pArray=fTrackVariables.GetArray(specifier);
         specifier+="[trackcount]/f";
         pTree->Branch(fTrackVariables.GetKey(i), pArray, specifier.Data());
     }
     for(i=0; i<fTrackVariablesInt.Variables(); i++){
         TString specifier=fTrackVariablesInt.GetKey(i);
         int* pArray=fTrackVariablesInt.GetArray(specifier);
         specifier+="[trackcount]/i";
         pTree->Branch(fTrackVariablesInt.GetKey(i), pArray, specifier.Data());
     }    
     if(fFillV0==kTRUE){
        for(i=0; i<fV0Variables.Variables(); i++){
            TString specifier=fV0Variables.GetKey(i);
            float* pArray=fV0Variables.GetArray(specifier);
            specifier+="[V0]/f";
            pTree->Branch(fV0Variables.GetKey(i), pArray, specifier.Data());
        }
     }
  } else {
    delete pTree;
    pTree = NULL;
  }
  
  return pTree;
}

void AliHLTGlobalHistoComponent::FillHistogramDefinitions(){
  /// default histogram definitions
}

int AliHLTGlobalHistoComponent::FillTree(TTree* pTree, const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ){

  /// fill the tree from the ESD
  int iResult=0;
  if (!IsDataEvent()) return 0;

  ResetVariables();

  // fetch ESD from input stream
  const TObject *obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  esd->GetStdContent();

  // fill track variables
  fNofTracks       = esd->GetNumberOfTracks();
  fVertexX         = esd->GetPrimaryVertexTracks()->GetX();
  fVertexY         = esd->GetPrimaryVertexTracks()->GetY();
  fVertexZ         = esd->GetPrimaryVertexTracks()->GetZ();
  fNofContributors = esd->GetPrimaryVertexTracks()->GetNContributors();
  fVertexStatus    = esd->GetPrimaryVertexTracks()->GetStatus();
  fNofV0s          = esd->GetNumberOfV0s();
  
  if(fFillV0==kTRUE && (fNofV0s > fMaxV0Count)){
     HLTWarning("Found V0s are %d, while component argument is %d, the respective TTree branch is not filled properly. Need to reconfigure.\n", fNofV0s, fMaxV0Count);
  }
  
  for(int i=0; i<fNofTracks; i++){    
      AliESDtrack *esdTrack = esd->GetTrack(i);
      if (!esdTrack) continue;
      
      Float_t DCAr, DCAz = -99;
      esdTrack->GetImpactParametersTPC(DCAr, DCAz);
      
      fTrackVariables.Fill("Track_pt"	     , esdTrack->Pt()			   );
      fTrackVariables.Fill("Track_phi"       , esdTrack->Phi()*TMath::RadToDeg()   );
      fTrackVariables.Fill("Track_eta"       , esdTrack->Eta()   		   );
      fTrackVariables.Fill("Track_p"	     , esdTrack->P()			   );
      fTrackVariables.Fill("Track_theta"     , esdTrack->Theta()*TMath::RadToDeg() );
      fTrackVariables.Fill("Track_TPCclus"   , esdTrack->GetTPCNcls()		   );
      fTrackVariables.Fill("Track_ITSclus"   , esdTrack->GetNcls(0)		   );
      fTrackVariables.Fill("Track_status"    , esdTrack->GetStatus()		   );
      fTrackVariables.Fill("Track_charge"    , esdTrack->Charge()		   );
      fTrackVariables.Fill("Track_DCAr"      , DCAr				   );
      fTrackVariables.Fill("Track_DCAz"      , DCAz				   );	
      fTrackVariables.Fill("Track_dEdx"      , esdTrack->GetTPCsignal() 	   );	
      fTrackVariablesInt.Fill("Track_status" , esdTrack->GetStatus()		   );	   
  }
  
  if(fFillV0==kTRUE){
     for(int i=0; i<fNofV0s; i++){     
   	 AliESDv0 *esdV0 = esd->GetV0(i);
   	 if(!esdV0) continue;
   	    
   	//fV0Variables.Fill("V0_AP", ap);
   	//fV0Variables.Fill("V0_pt", pt); 
   	//fV0Variables.Fill("clust1", t1->GetTPCNcls()); 
   	//fV0Variables.Fill("clust2", t2->GetTPCNcls()); 
   	//fV0Variables.Fill("dev1", dev1); 
   	//fV0Variables.Fill("dev2", dev2); 
   	//fV0Variables.Fill("devPrim", devPrim); 
   	//fV0Variables.Fill("length", length); 
   	//fV0Variables.Fill("sigmaLength", sigmaLength); 
   	//fV0Variables.Fill("r", r); 
   	  
     } // end of loop over V0s
  }
  
  if(iResult<0){
    // fill an empty event
    ResetVariables();
  }
  
  fEvent++;
  pTree->Fill();
  return iResult;
}

int AliHLTGlobalHistoComponent::ResetVariables(){
/// reset all filling variables
  fNofTracks=0;
  fNofV0s=0;
  fTrackVariables.ResetCount();
  fTrackVariablesInt.ResetCount();
  fV0Variables.ResetCount();
  return 0;
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOriginDataType() const{
// get the origin of the output data
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT;
}

int AliHLTGlobalHistoComponent::ScanConfigurationArgument(int argc, const char** argv){
/// inherited from AliHLTComponent, scan argument

  if (argv==NULL || argc<1) return 0;

  int i=0;
  TString argument = argv[i];

  // -max-track-count
  if(argument.CompareTo("-max-track-count")==0){    
     if (++i>=argc) return -EPROTO; // missing parameter
     argument = argv[i];
     fMaxTrackCount = argument.Atoi();    
     HLTInfo("got %s with parameter %s", argument.Data(), argv[i]);
     return ++i; 
  }
 
  // -max-V0-count 
  if(argument.CompareTo("-max-V0-count")==0){    
     if (++i>=argc) return -EPROTO; // missing parameter
     argument = argv[i];
     fMaxV0Count = argument.Atoi();    
     HLTInfo("got %s with parameter %s", argument.Data(), argv[i]);
     return ++i;
  }
  
  // -fill-V0
  if(argument.CompareTo("-fill-V0")==0){
     fFillV0 = kTRUE;
     HLTInfo("got %s", argument.Data());
     return ++i;
  }
  // no recognized argument, forward to base class
  return AliHLTTTreeProcessor::ScanConfigurationArgument(argc, argv);
}

int AliHLTGlobalHistoComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/){  
// see header file for class documentation

  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=GetComponentID();
  }

  return ConfigureFromCDBTObjString(cdbPath.Data());
}


