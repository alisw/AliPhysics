// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelectionMCEl.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  MC El selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "THnSparse.h"
#include "AliAODMCParticle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "AliAODTrack.h"
#include "AliReducedParticle.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;
using std::vector;


/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionMCEl)

AliDxHFEParticleSelectionMCEl::AliDxHFEParticleSelectionMCEl(const char* opt)
  : AliDxHFEParticleSelectionEl(opt)
  , fMCTools()
  , fPDGnotMCElectron(NULL)
  , fOriginMother(0)
  , fResultMC(0)
  , fUseKine(kFALSE)
  , fMotherPDGs(0)
  , fUseMCReco(kFALSE)
  , fSelectionStep(AliDxHFEParticleSelectionEl::kNoCuts)
{
  // constructor
  // 
  // 
  // 
  // 

  fMCTools.~AliDxHFEToolsMC();

  ParseArguments(opt);

  // This is also checked in AddTasks, but just for safety! 
  if(fUseMCReco && fUseKine) AliFatal("CAN'T SET BOTH usekine AND elmcreco AT THE SAME TIME");

  // TODO: argument scan, build tool options accordingly
  // e.g. set mc mode first/last, skip control histograms
  TString toolopt("pdg=11 mc-last");
  if(fUseKine) toolopt+=" usekine";
  new (&fMCTools) AliDxHFEToolsMC(toolopt);
}

//at the moment not used, keep for now (to be used when plotting)
const char* AliDxHFEParticleSelectionMCEl::fgkPDGMotherBinLabels[]={
  "d",
  "u",
  "s",
  "c",
  "b",
  "gluon",
  "gamma",
  "#pi^{0}",
  "#eta",
  "proton",
  "others"
};

                             
const char*  AliDxHFEParticleSelectionMCEl::fgkPDGBinLabels[]={
  "positron",
  "electron",
  "#mu+",
  "#mu-",
  "#pi+",
  "#pi-",
  "K+",
  "K-",
  "proton",
  "antiproton",
  "others"
};

AliDxHFEParticleSelectionMCEl::~AliDxHFEParticleSelectionMCEl()
{
  // destructor

  if(fPDGnotMCElectron){
    delete fPDGnotMCElectron;
    fPDGnotMCElectron=NULL;
  }

}

int AliDxHFEParticleSelectionMCEl::Init()
{
  //
  // init function
  // 

  // Particles considered HFE background. can be expanded
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGpi0); 
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGeta);
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGgamma);
  fMotherPDGs.push_back(AliDxHFEToolsMC::kPDGJpsi);

  int iResult=0;
  iResult=AliDxHFEParticleSelectionEl::Init();
  if (iResult<0) return iResult;

  // Histo containing PDG of track which was not MC truth electron
  // TODO: Add to list in baseclass
  fPDGnotMCElectron=CreateControlHistogram("fPDGnotMCElectron","PDG of track not MC truth electron",AliDxHFEToolsMC::kNofPDGLabels,fgkPDGBinLabels);  
  AddControlObject(fPDGnotMCElectron);
  return 0;
}

THnSparse* AliDxHFEParticleSelectionMCEl::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 

  const int thnSize = 4;
  InitTHnSparseArray(thnSize);
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  //     		       0    1      2     3     
  // 	 	               Pt   Phi   Eta   mother 
  int    thnBins[thnSize] = { 1000,  200, 500,    15  };
  double thnMin [thnSize] = {    0,    0, -1.,  -1.5  };
  double thnMax [thnSize] = {  100, 2*Pi,  1.,  13.5  };
  const char* thnNames[thnSize]={
    "Pt",
    "Phi",
    "Eta", 
    "Mother", //bin==-1: Not MC truth electron
  };
 
  return CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);
}

int AliDxHFEParticleSelectionMCEl::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliAODTrack *track=(AliAODTrack*)p;
  if (!track) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=track->Eta();
  data[i++]=fOriginMother; 
  
  return i;
}

int AliDxHFEParticleSelectionMCEl::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// overloaded from AliDxHFEParticleSelection: check particle
  /// H: Have changed function. Now doing particle selection first, then run MC over 
  /// selected tracks. Could configure it to be configurable, but not sure if it
  /// is needed.  
  /// Result from normal track selection is returned, result from MC is stored in
  /// THnSparse. 

  int iResult=0;
  if (!p || !pEvent){
    return -EINVAL;
  }
  fOriginMother=-1;

  if(!fUseKine && !fUseMCReco){
  // step 1:
  // optional MC selection before the particle selection
  if (fMCTools.MCFirst() && (iResult=CheckMC(p, pEvent))==0) {
    // histograming?
    return iResult;
  }

  // step 2 or 1, depending on sequence:
  // normal particle selection
  iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
  if (fMCTools.MCFirst() || iResult==0) return iResult;


  }

  // step 2, only executed if MC check is last
  // optional MC selection after the particle selection
  // result stored to be filled into THnSparse
  // TODO: strictly speaken the particles should be rejected
  // if not mc selected, however skip this for the moment, because of
  // the logic outside
  // This line will run always when running directly over the stack or over MC reconstructed tracks
  // For MC reconstructed tracks, should consider adding more constrictions on the tracks,
  // maybe do the track selection first (which means could call AliDxHFEParticleSelectionEl::IsSelected()
  // and don't do PID
  
  if(fUseMCReco && ( fSelectionStep > AliDxHFEParticleSelectionEl::kNoCuts )){
    iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
    if(iResult ==0) return iResult;
  }
  fResultMC=CheckMC(p, pEvent);
  
  if(fUseKine || fUseMCReco)
    return fResultMC;

  return iResult;
   
}

int AliDxHFEParticleSelectionMCEl::CheckMC(AliVParticle* p, const AliVEvent* pEvent)
{
  /// check if MC criteria are fulfilled
  if (!p || !pEvent){
    return -EINVAL;
  }
  int iResult=0;

  if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
    // TODO: message? but has to be filtered in order to avoid message flood
    return 0; // no meaningful filtering on mc possible
  }

  if(fUseKine){
    // If run on kinematical level, check if particle is electron (only ones who are relevant) 
    // and if pass IsPhysicalPrimary()
    Int_t test = fMCTools.CheckMCParticle(p);
    if(test==1) return 0;
    // Add additional contraints on the electrons? 
    // remove delta e? also remove E<300MeV?
    // Should also mark dalitz decay and gamma conversion..
    AliAODMCParticle *mcp=dynamic_cast<AliAODMCParticle*>(p);
    if(!mcp->IsPhysicalPrimary()) return 0; 

  }

  int pdgParticle=-1;
  if (!fUseKine && fMCTools.RejectByPDG(p,false, &pdgParticle)) {
    // rejected by pdg
    // TODO: Move this to fMCTools???? Can this be part of the statistics in the MC class?
    fPDGnotMCElectron->Fill(fMCTools.MapPDGLabel(pdgParticle));
    return 0;
  }

  int pdgMother=0;
  // Find PDG of first mother
  pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother);

  // Check if first mother is counted as background
  Bool_t isNotBackground=fMCTools.RejectByPDG(pdgMother,fMotherPDGs);

  if(!isNotBackground){
    // Set fOriginMother if mother is counted as background
    // TODO: Could this be done in a more elegant way?
    switch(pdgMother){
    case(AliDxHFEToolsMC::kPDGpi0): fOriginMother=AliDxHFEToolsMC::kNrOrginMother; break;
    case(AliDxHFEToolsMC::kPDGeta): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+1; break;
    case(AliDxHFEToolsMC::kPDGgamma): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+2;break;
    case(AliDxHFEToolsMC::kPDGJpsi): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+3;break;
    }
  }
  else{
    // If loop over Stack, also checks if First mother is a HF meson
    Bool_t isHFmeson =fMCTools.TestMotherHFMeson(TMath::Abs(pdgMother));

    if(isHFmeson){
      // If first mother is a HF meson, loops back to find 
      // original quark + if there was a gluon. Result is 
      // stored in fOriginMother
      pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
      fOriginMother=fMCTools.GetOriginMother();
    }
    else{
      fOriginMother=AliDxHFEToolsMC::kNrOrginMother+4;
    }

  }

  return 1;
}

void AliDxHFEParticleSelectionMCEl::Clear(const char* option)
{
  /// clear internal memory
  fMCTools.Clear(option);
}

AliVParticle *AliDxHFEParticleSelectionMCEl::CreateParticle(AliVParticle* track)
{

  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),track->Charge(),fOriginMother);

  return part;

}

int AliDxHFEParticleSelectionMCEl::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;

  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if (argument.BeginsWith("usekine") ){
      fUseKine=kTRUE;
      continue;
    }
    if (argument.BeginsWith("elmcreco")){
      fUseMCReco=kTRUE;
      if(argument.BeginsWith("elmcreco=")){
	argument.ReplaceAll("elmcreco=", "");
	if(argument.CompareTo("aftertrackcuts")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTPC;
	if(argument.CompareTo("aftertofpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOF;
	if(argument.CompareTo("afterfullpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOFTPC;

	AliDxHFEParticleSelectionEl::SetFinalCutStep(fSelectionStep);
    
      }
      continue;
    }
    // forwarding of single argument works, unless key-option pairs separated
    // by blanks are introduced
    AliDxHFEParticleSelection::ParseArguments(argument);
  }
  
  return 0;
}
