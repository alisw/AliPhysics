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
#include "TList.h"
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
  , fHistoList(NULL)
  , fOriginMother(0)
  , fResultMC(0)
  , fUseKine(kFALSE)
  , fMotherPDGs(0)
  , fUseMCReco(kFALSE)
  , fSelectionStep(AliDxHFEParticleSelectionEl::kNoCuts)
  , fStoreCutStepInfo(kFALSE)
  , fElSelection(kAllPassingSelection)
  , fStoreOnlyMCElectrons(kFALSE)
  , fMCInfo(kMCLast)
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
  TString toolopt("pdg=11");
  if(fMCInfo==kMCLast) toolopt+=" mc-last";
  if(fMCInfo==kMCOnly) toolopt+=" mc-first";
  if(fMCInfo==kMCFirst) toolopt+=" mc-first";
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
  if(fHistoList){
    delete fHistoList;
    fHistoList=NULL;
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
  fHistoList=new TList;
  fHistoList->SetName("HFE MC Histograms");
  fHistoList->SetOwner();

  fHistoList->Add(CreateControlHistogram("fPDGnotMCElectron","PDG of track not MC truth electron",AliDxHFEToolsMC::kNofPDGLabels,fgkPDGBinLabels));
  fHistoList->Add(CreateControlHistogram("fPDGNotHFMother","PDG of mother not HF",5000));  
  fHistoList->Add(CreateControlHistogram("fPDGHadronMother","PDG of mother of hadrons",5000));  

  AddControlObject(fHistoList);

  return 0;
}

THnSparse* AliDxHFEParticleSelectionMCEl::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 

  const double Pi=TMath::Pi();
  TString name;
  THnSparse* thn=NULL;
  name.Form("%s info", GetName());

  if(fStoreCutStepInfo){
    const int thnSizeExt =5;

    InitTHnSparseArray(thnSizeExt);

    // TODO: Redo binning of distributions more?
    //     		              0    1      2     3     
    // 	 	                      Pt   Phi   Eta   mother 
    int    thnBinsExt[thnSizeExt] = { 100,  100, 100,    15, kNCutLabels-1 };
    double thnMinExt [thnSizeExt] = {   0,    0, -1.,  -1.5, kRecKineITSTPC-0.5};
    double thnMaxExt [thnSizeExt] = {  10, 2*Pi,  1.,  13.5, kSelected-0.5};
    const char* thnNamesExt[thnSizeExt]={
      "Pt",
      "Phi",
      "Eta", 
      "Mother", //bin==-1: Not MC truth electron
      "Last survived cut step"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,thnSizeExt,thnBinsExt,thnMinExt,thnMaxExt,thnNamesExt);
  }
  else{

    const int thnSize =4;
    InitTHnSparseArray(thnSize);

    // TODO: Redo binning of distributions more?
    //     		       0       1      2     3     
    // 	 	               Pt     Phi   Eta   mother 
    int    thnBins[thnSize] = { 100,  50, 100,    15  };
    double thnMin [thnSize] = {   0,    0, -1.,  -1.5  };
    double thnMax [thnSize] = {  10, 2*Pi,  1.,  13.5  };
    const char* thnNames[thnSize]={
      "Pt",
      "Phi",
      "Eta", 
      "Mother", //bin==-1: Not MC truth electron
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);

  }
  return thn;
}

int AliDxHFEParticleSelectionMCEl::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  if (!p) return -ENODATA;
  // handle different types of tracks, can be extended
  AliReducedParticle *trRP=dynamic_cast<AliReducedParticle*>(p);
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  memset(data, 0, dimension*sizeof(data[0]));
  data[i++]=p->Pt();
  data[i++]=p->Phi();
  data[i++]=p->Eta();
  if (trRP) data[i]=trRP->GetOriginMother();
  else data[i]=fOriginMother;
  i++; // take out of conditionals to be save
  if (i<dimension) {
    if(fStoreCutStepInfo) data[i]=GetLastSurvivedCutsStep();
    i++; // take out of conditionals to be save
  }
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
  fResultMC=0;
  if (!p || !pEvent){
    return -EINVAL;
  }
  fOriginMother=-1;

  if(!fUseKine && !fUseMCReco){
  // step 1:
  // optional MC selection before the particle selection
    if (fMCTools.MCFirst()){
      fResultMC=CheckMC(p, pEvent);
      // histograming?
      if(fMCInfo==kMCOnly) return fResultMC;
      if(fResultMC==0) return fResultMC;
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
  
  if(fUseMCReco) {// && ( fSelectionStep > AliDxHFEParticleSelectionEl::kNoCuts )){
    // always check the base class method in mode fUseMCReco
    iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
    if(iResult == 0) return iResult;
  }
  fResultMC=CheckMC(p, pEvent);
  
  // Return only result of MC checks (rather than result from reconstruction) for kinematical studies, 
  // when looking at the various selection steps (want only electrons) or when we only want to store specific
  // sources from the reconstruction (hadrons/nonHFE/HFE/onlyc/onlyb)
  if(fUseKine || fUseMCReco || fElSelection>kAllPassingSelection || fStoreOnlyMCElectrons)
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
    // TODO: remove delta e? also remove E<300MeV?
    // Should also mark dalitz decay and gamma conversion..
    AliAODMCParticle *mcp=dynamic_cast<AliAODMCParticle*>(p);
    if(!mcp || !mcp->IsPhysicalPrimary()) return 0; 
    if(TMath::Abs(mcp->Eta())>0.8){ return 0;}

  }

  int pdgMother=0;
  int pdgParticle=-1;
  if (!fUseKine && fMCTools.RejectByPDG(p,false, &pdgParticle)) {
    // rejected by pdg
    // TODO: Move this to fMCTools???? Can this be part of the statistics in the MC class?
    ((TH1D*)fHistoList->FindObject("fPDGnotMCElectron"))->Fill(fMCTools.MapPDGLabel(pdgParticle));
    ((TH1D*)fHistoList->FindObject("fPDGHadronMother"))->Fill(fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother));

    // If the hadrons (and only hadrons) are going to be stored
    if(fElSelection==kHadron) return 1;
    else return 0;
  }

  // Find PDG of first mother
  pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother);

  // Check if first mother is counted as background
  Bool_t isNotBackground=fMCTools.RejectByPDG(pdgMother,fMotherPDGs);
  Int_t selection=-1;

  if(!isNotBackground){
    // Set fOriginMother if mother is counted as background
    // TODO: Could this be done in a more elegant way?
    switch(pdgMother){
    case(AliDxHFEToolsMC::kPDGpi0): fOriginMother=AliDxHFEToolsMC::kNrOrginMother; break;
    case(AliDxHFEToolsMC::kPDGeta): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+1; break;
    case(AliDxHFEToolsMC::kPDGgamma): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+2;break;
    case(AliDxHFEToolsMC::kPDGJpsi): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+3;break;
    }
    selection=kNonHFE;
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
      //NotHFmother
      ((TH1D*)fHistoList->FindObject("fPDGNotHFMother"))->Fill(pdgMother);
      fOriginMother=AliDxHFEToolsMC::kNrOrginMother+4;
      selection=kNonHFE;
    }

  }

  // Checks on the electrons, to return only specific sources
  if(fOriginMother >= AliDxHFEToolsMC::kOriginNone && fOriginMother <= AliDxHFEToolsMC::kOriginStrange)
    selection=kNonHFE;

  if(fOriginMother == AliDxHFEToolsMC::kOriginCharm || fOriginMother == AliDxHFEToolsMC::kOriginGluonCharm)
    selection=kOnlyc;

  if(fOriginMother == AliDxHFEToolsMC::kOriginBeauty || fOriginMother == AliDxHFEToolsMC::kOriginGluonBeauty)
    selection=kOnlyb;

  if(fElSelection==kNonHFE && ( selection!= kNonHFE)){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kHFE && (selection != kOnlyc && selection != kOnlyb)) {
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kOnlyc && selection !=kOnlyc){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
  }
  if(fElSelection==kOnlyb && selection !=kOnlyb){
    AliDebug(2,Form("Particle selected as: %d, want only to select %d",selection, fElSelection));
    return 0;
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
	if(argument.CompareTo("alltracks")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kNoCuts;
	else if(argument.CompareTo("afterreckineitstpc")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRecKineITSTPC;
	else if(argument.CompareTo("afterrecprim")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kRecPrim;
	else if(argument.CompareTo("afterhfeits")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsITS;
	else if(argument.CompareTo("afterhfetof")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTOF;
	else if(argument.CompareTo("afterhfetpc")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTPC;
	else if(argument.CompareTo("aftertrackcuts")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kHFEcutsTPC;
	else if(argument.CompareTo("aftertofpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOF;
	else if(argument.CompareTo("aftertpcpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTPC;
	else if(argument.CompareTo("afterfullpid")==0) fSelectionStep=AliDxHFEParticleSelectionEl::kPIDTOFTPC;
	else AliFatal(Form("unknown argument '%s'", argument.Data()));

	AliDxHFEParticleSelectionEl::SetFinalCutStep(fSelectionStep);
      }
      continue;
    }
    if(argument.BeginsWith("storeonlyMCelectrons")){
      AliInfo("Store only MC truth electrons");
      fStoreOnlyMCElectrons=kTRUE;
      continue;
    }
    if(argument.BeginsWith("storelastcutstep")){
      AliInfo("Stores the last cut step");
      fUseMCReco=kTRUE;
      fStoreCutStepInfo=kTRUE;
      AliDxHFEParticleSelectionEl::SetStoreLastCutStep(kTRUE);
      continue;
    }
    if(argument.BeginsWith("ElSelection=")){
      argument.ReplaceAll("ElSelection=","");
      if(argument.CompareTo("hadron")==0){ fElSelection=kHadron;}
      else if(argument.CompareTo("nonHFE")==0) {fElSelection=kNonHFE;}
      else if(argument.CompareTo("HFE")==0){ fElSelection=kHFE;}
      else if(argument.CompareTo("Onlyc")==0){ fElSelection=kOnlyc;}
      else if(argument.CompareTo("Onlyb")==0){ fElSelection=kOnlyb;}
      else AliFatal(Form("unknown argument '%s'", argument.Data()));
      AliInfo(Form("Selecting only source %d",fElSelection));
      continue;
    }
    if(argument.BeginsWith("mc-only")){
      AliInfo("Do only test on MC info");
      fMCInfo=kMCOnly;
      continue;
    }
    if(argument.BeginsWith("mc-first")){
      AliInfo("Do test on MC info first");
      fMCInfo=kMCFirst;
      continue;
    }
    if(argument.BeginsWith("mc-last")){
      AliInfo("Do test on MC info last");
      fMCInfo=kMCLast;
      continue;
    }
    // forwarding of single argument works, unless key-option pairs separated
    // by blanks are introduced
    AliDxHFEParticleSelection::ParseArguments(argument);
  }
  
  return 0;
}
