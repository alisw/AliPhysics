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

/// @file   AliDxHFEParticleSelectionMCD0.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  MC D0 selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCuts.h"
#include "AliVParticle.h"
#include "AliReducedParticle.h"
#include "AliRDHFCuts.h"
#include "THnSparse.h"
#include "TH1F.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionMCD0)

AliDxHFEParticleSelectionMCD0::AliDxHFEParticleSelectionMCD0(const char* opt)
  : AliDxHFEParticleSelectionD0(opt)
  , fMCTools()
  , fPDGnotMCD0(NULL)
  , fResultMC(0)
  , fOriginMother(0)
  , fUseKine(kFALSE)  
  , fD0PropertiesKine(NULL)
  , fStoreOnlyMCD0(kFALSE)
  , fMCInfo(kMCLast)
  , fRequireD0toKpi(kFALSE)
  , fSystem(0)
  , fUseCentrality(0)
{
  // constructor
  // 
  // 
  // 
  ParseArguments(opt);

  // TODO: argument scan, pass only relevant arguments to tools
  fMCTools.~AliDxHFEToolsMC();
  TString toolopt("pdg=421");
  if(fMCInfo==kMCLast) toolopt+=" mc-last";
  if(fMCInfo==kMCOnly) toolopt+=" mc-first";
  if(fMCInfo==kMCFirst) toolopt+=" mc-first";

  if(fUseKine) toolopt+=" usekine";  
  new (&fMCTools) AliDxHFEToolsMC(toolopt);
}

AliDxHFEParticleSelectionMCD0::~AliDxHFEParticleSelectionMCD0()
{
  // destructor
  if (fD0PropertiesKine) {
    delete fD0PropertiesKine;
    fD0PropertiesKine=NULL;
  }
  if(fPDGnotMCD0){
    delete fPDGnotMCD0;
    fPDGnotMCD0=NULL;
  }
}

int AliDxHFEParticleSelectionMCD0::InitControlObjects()
{
  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly
  AliInfo("Setting up control objects");

  if(fUseKine) {
    fD0PropertiesKine=DefineTHnSparse();
    AddControlObject(fD0PropertiesKine);
    return AliDxHFEParticleSelection::InitControlObjects();
  }
  else{
    return AliDxHFEParticleSelectionD0::InitControlObjects();
  }
}

THnSparse* AliDxHFEParticleSelectionMCD0::DefineTHnSparse()
{
  //
  // Defines the THnSparse.

  // here is the only place to change the dimension
  const int thnSize2 = 6;
  const int thnSize3 = 5;
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  if(2==fSystem){//Reduced binning for p-Pb
    InitTHnSparseArray(thnSize3);
    // 			             0     1        3        4     5
    // 	 	                     Pt   Phi   D0InvMass  Eta  mother 
    int         thnBins [thnSize3] = {28, 100,    150,     100,  10  };
    double      thnMin  [thnSize3] = {   2,  0,  1.5848,   -1., -1.5  };
    double      thnMax  [thnSize3] = { 16, 2*Pi, 2.1848,    1.,  8.5  };
    const char* thnNames[thnSize3] = {
      "Pt",
      "Phi",
      "D0InvMass", 
      "Eta",
      "Mother of D0"  // Bin -1 = not MC truth D0, rest OK
    };
    // Add Histo displaying pdg of D0 candidates not passing MatchToMC()
    // TODO: Add it to the TList of D0 main class
    fPDGnotMCD0= new TH1F("fPDGnotMCD0","PDG of track not MC truth D0",1002,-2.5,999.5);
    AddControlObject(fPDGnotMCD0);
    return CreateControlTHnSparse(name,thnSize3,thnBins,thnMin,thnMax,thnNames);
  }
  else{
    InitTHnSparseArray(thnSize2);
    // 			             0     1      2       3        4     5
    // 	 	                     Pt   Phi   Ptbin  D0InvMass  Eta  mother 
    int         thnBins [thnSize2] = {1000, 200,   15,     200,     500,  10  };
    double      thnMin  [thnSize2] = {   0,  0,     0,    1.5648,   -1., -1.5  };
    double      thnMax  [thnSize2] = { 100, 2*Pi,  14,    2.1648,    1.,  8.5  };
    const char* thnNames[thnSize2] = {
      "Pt",
      "Phi",
      "Ptbin", 
      "D0InvMass", 
      "Eta",
      "Mother of D0"  // Bin -1 = not MC truth D0, rest OK
    };
    // Add Histo displaying pdg of D0 candidates not passing MatchToMC()
    // TODO: Add it to the TList of D0 main class
    fPDGnotMCD0= new TH1F("fPDGnotMCD0","PDG of track not MC truth D0",1002,-2.5,999.5);
    AddControlObject(fPDGnotMCD0);
    return CreateControlTHnSparse(name,thnSize2,thnBins,thnMin,thnMax,thnNames);
  }
}

int AliDxHFEParticleSelectionMCD0::HistogramParticleProperties(AliVParticle* p, int selectionCode)
{

  // When looping on kinematical level, need a different HistogramParticleProperties than on reconstructed tracks
  if(fUseKine){
    /// histogram particle properties
    if (!p) return -EINVAL;

    // fill the common histograms
    AliDxHFEParticleSelection::HistogramParticleProperties(p, selectionCode);

    // no daughters to fill if 0 (= no candidate)
    if (selectionCode==0){
      return 0;
    }
    AliAODMCParticle* partMC=dynamic_cast<AliAODMCParticle*>(p);

    if(!partMC) {
      return 0;
    }

    SetInvMass(partMC->GetCalcMass());
    AliRDHFCuts *cuts=GetHFCuts();
    int ptbin=cuts->PtBin(partMC->Pt());
    SetPtBin(ptbin);

    // Fills only for D0 or both.. 
    if ((selectionCode==1 || selectionCode==3) && GetFillOnlyD0D0bar()<2) {
      if(fD0PropertiesKine && ParticleProperties()) {
	memset(ParticleProperties(), 0, GetDimTHnSparse()*sizeof(ParticleProperties()[0]));
	FillParticleProperties(p, ParticleProperties(), GetDimTHnSparse());
	fD0PropertiesKine->Fill(ParticleProperties());
      }
    }

    // Fills for D0bar or both
    if ((selectionCode==1 || selectionCode==2) && (GetFillOnlyD0D0bar()==0 || GetFillOnlyD0D0bar()==2)) {
      if(fD0PropertiesKine && ParticleProperties()) {
	memset(ParticleProperties(), 0, GetDimTHnSparse()*sizeof(ParticleProperties()[0]));
	FillParticleProperties(p, ParticleProperties(), GetDimTHnSparse());
	fD0PropertiesKine->Fill(ParticleProperties());
      }

    }
    return 0;
  }
  else {
    return AliDxHFEParticleSelectionD0::HistogramParticleProperties(p,selectionCode);
  }
}

int AliDxHFEParticleSelectionMCD0::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
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
  if(fSystem!=2) data[i++]=AliDxHFEParticleSelectionMCD0::GetPtBin(); 
  data[i++]=AliDxHFEParticleSelectionMCD0::GetInvMass();
  data[i++]=track->Eta();
  data[i++]=fOriginMother; // at the moment not included background. Should expand

  return i;
}

int AliDxHFEParticleSelectionMCD0::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// overloaded from AliDxHFEParticleSelection: check particle
  /// H: Have changed function. Now doing particle selection first, then run MC over 
  /// selected tracks. Could configure it to be configurable, but not sure if it
  /// is needed.  
  /// result from normal track selection is returned, result from MC is stored in
  /// THnSparse. 

  int iResult=0;
  fResultMC=0;
  fOriginMother=-1;
  if(fUseKine){
    // Will here loop on all tracks in the stack, and checks whether they are D0s (through CheckMCParticle())
    if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
      return 0; // no meaningful filtering on mc possible
    }
    bool result = fMCTools.CheckMCParticle(p);
    if(!result) return 0;
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(p);
    if (!mcPart) {
      AliWarning("Particle not found in tree, skipping"); 
      return 0;
    } 
    AliRDHFCuts* cuts=const_cast<AliRDHFCuts*>(GetHFCuts());
    if (!cuts) {
      return 0;
    } 
    else if(cuts->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y()) ) {
      if(TMath::Abs(mcPart->Eta()) < 0.8 ){
	fMCTools.FindMotherPDG(p);
	fOriginMother=fMCTools.GetOriginMother();
	if(fRequireD0toKpi){

	  TClonesArray* fMCArray = dynamic_cast<TClonesArray*>(fMCTools.GetMCArray());
	  if(!fMCArray) {cout << "no array" << endl; return 0;}

	  //Make sure only two daughters
	  int nrDaughters=mcPart->GetNDaughters();
	  if(nrDaughters!= 2)
	    return 0;

	  Int_t label0 = mcPart->GetDaughterLabel(0);
	  Int_t label1 = mcPart->GetDaughterLabel(1);
	  AliDebug(2,Form("label0 = %d, label1 = %d",label0,label1));

	  if (label1<=0 || label0 <= 0){
	    AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
	    return 0;  
	  }

	  //fetch daughters
	  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label0));
	  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label1));
	  int pdg0 = TMath::Abs(mcPartDaughter0->PdgCode());
	  int pdg1 = TMath::Abs(mcPartDaughter1->PdgCode());

	  //Check pdg of daughters
	  if(!(pdg0==AliDxHFEToolsMC::kPDGkaon || pdg0==AliDxHFEToolsMC::kPDGpion))
	    return 0;
	  if(!(pdg1==AliDxHFEToolsMC::kPDGkaon || pdg1==AliDxHFEToolsMC::kPDGpion))
	    return 0;
	}
	iResult=1;
      }
    }

    //TODO: Should also return whether D0 or D0bar... at the moment only care of absolute value of D0
    return iResult;
  }
  else{
  // step 1:
  // MC selection
  if (fMCTools.MCFirst()){
    fResultMC=CheckMC(p, pEvent);
    // histograming?
    if(fMCInfo==kMCOnly) return fResultMC;
    if(fResultMC==0) return fResultMC;
  }

  // step 2 or 1, depending on sequence:
  // normal particle selection
  iResult=AliDxHFEParticleSelectionD0::IsSelected(p, pEvent);
  if (fMCTools.MCFirst() || iResult==0) return iResult;

  // step 2, only executed if MC check is last
  // MC selection  - > Should maybe also distinguish between D0 and D0bar
  // result stored to be filled into THnSparse
  // TODO: strictly speaken the particles should be rejected
  // if not mc selected, however skip this for the moment, because of
  // the logic outside
  fResultMC=CheckMC(p, pEvent);

  if(fStoreOnlyMCD0) return fResultMC;
  
  return iResult;
  }

  return 0;
}

int AliDxHFEParticleSelectionMCD0::CheckMC(AliVParticle* p, const AliVEvent* pEvent)
{
  /// check if MC criteria are fulfilled
  // Check both D0 and D0bar (for now only D0)

  if (!p || !pEvent){
    return -EINVAL;
  }
  int iResult=0;

  if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
    // TODO: message? but has to be filtered in order to avoid message flood
    return 0; // no meaningful filtering on mc possible
  }

  AliAODRecoDecayHF2Prong *particle = dynamic_cast<AliAODRecoDecayHF2Prong*>(p);

  if(!particle) return 0;

  Int_t pdgDgD0toKpi[2]={AliDxHFEToolsMC::kPDGkaon,AliDxHFEToolsMC::kPDGpion};

  TClonesArray* fMCArray = dynamic_cast<TClonesArray*>(fMCTools.GetMCArray());
  if(!fMCArray) {cout << "no array" << endl; return -1;}

  // find associated MC particle for D0->Kpi
  Int_t MClabel=-9999;

  //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx). Checks both D0s and daughters
  MClabel=particle->MatchToMC(AliDxHFEToolsMC::kPDGD0,fMCArray,2,pdgDgD0toKpi); 
  
  // TODO: Need a different strategy!!!
  // ALSO: look at AliAnalysisTaskSED0Mass for tips
  if(MClabel<0){
    // Checking PDG of particle if not MC truth D0
    // TODO: done the right way??
    Int_t MCl = p->GetLabel();
    if(MCl<0) {
      fPDGnotMCD0->Fill(-2);
      return 0;
    }
    int pdgPart=-1;
    AliAODMCParticle* aodmcp=0;
    aodmcp=dynamic_cast<AliAODMCParticle*>(fMCArray->At(MCl));
    if (aodmcp)
      pdgPart=TMath::Abs(aodmcp->GetPdgCode());
    if (pdgPart<0){
      fPDGnotMCD0->Fill(-1);
      return 0;
    }
    else{
      fPDGnotMCD0->Fill(pdgPart);
    }
    fOriginMother=-1;
    return 0;
  }

  fMCTools.SetMClabel(MClabel);
  fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
  fOriginMother=fMCTools.GetOriginMother();

  return 1;
}

void AliDxHFEParticleSelectionMCD0::Clear(const char* option)
{
  /// clear internal memory
  fMCTools.Clear(option);
}

AliVParticle *AliDxHFEParticleSelectionMCD0::CreateParticle(AliVParticle* track)
{
  //
  //Created object which contain variables needed for correlation. 
  //

  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),AliDxHFEParticleSelectionMCD0::GetInvMass(),AliDxHFEParticleSelectionMCD0::GetPtBin(), fOriginMother);

  return part;

}

int AliDxHFEParticleSelectionMCD0::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
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
    if(argument.BeginsWith("storeonlyMCD0")){
      AliInfo("Store only MC truth D0");
      fStoreOnlyMCD0=kTRUE;
      continue;
    }
    if(argument.BeginsWith("RequireD0toKpi")|| argument.BeginsWith("requireD0toKpi")){
      AliInfo("Store only D0 to Kpi");
      fRequireD0toKpi=kTRUE;
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
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) {fSystem=0; fUseCentrality=0;}
      else if (argument.CompareTo("Pb-Pb")==0) {fSystem=1; fUseCentrality=1;}
      else if (argument.CompareTo("p-Pb")==0) {fSystem=2; fUseCentrality=0;}
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    AliDxHFEParticleSelection::ParseArguments(argument);
  }
  return 0;
}

double AliDxHFEParticleSelectionMCD0::GetD0Eff(AliVParticle* tr){


  Double_t D0eff=1;
  AliReducedParticle *track=(AliReducedParticle*)tr;
  if (!track) return -ENODATA;
  Double_t pt=track->Pt();
  //  Double_t origin=track->GetOriginMother();

  //  Bool_t isCharm=(fOriginMother==AliDxHFEToolsMC::kOriginCharm || 
  //		  fOriginMother==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(fOriginMother==AliDxHFEToolsMC::kOriginBeauty || 
		   fOriginMother==AliDxHFEToolsMC::kOriginGluonBeauty);

  AliHFAssociatedTrackCuts* cuts=dynamic_cast<AliHFAssociatedTrackCuts*>(GetEffCutObject());
  if (!cuts) {
    if (GetEffCutObject())
      AliError(Form("cuts object of wrong type %s, required AliHFAssociatedTrackCuts", cuts->ClassName()));
    else
      AliError("mandatory cuts object missing");
    return -EINVAL;
  }

  if(isBeauty)
    D0eff=cuts->GetTrigWeightB(pt,GetEventMult());
  else
    D0eff=cuts->GetTrigWeight(pt,GetEventMult());

  return D0eff;

}
