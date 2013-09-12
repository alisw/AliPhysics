
//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
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

// $Id$

/// @file   AliDxHFECorrelationMC.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-11-02
/// @brief  Worker class for DxHFE correlation
///

#include "AliDxHFECorrelationMC.h"
#include "AliVParticle.h"
#include "AliLog.h"
//#include "AliAnalysisCuts.h"         // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"      // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h"   // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"          // required dependency libCORRFW.so
//#include "AliAODRecoDecayHF2Prong.h"   // libPWGHFvertexingHF
#include "TObjArray.h"
#include "AliDxHFECorrelation.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "AliReducedParticle.h"
#include "AliDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliDxHFEToolsMC.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelationMC)

AliDxHFECorrelationMC::AliDxHFECorrelationMC(const char* name)
  : AliDxHFECorrelation(name?name:"AliDxHFECorrelationMC")
  , fMCEventType(0)
  , fStoreOriginEl(kAll)
  , fStoreOriginD(kAll)
  , fD0EffMapP(NULL)
  , fD0EffMapFD(NULL)
{
  // default constructor
  // 
  //

}


AliDxHFECorrelationMC::~AliDxHFECorrelationMC()
{
  // destructor
  //
  //
  fD0EffMapFD=NULL;
  fD0EffMapP=NULL;

}

THnSparse* AliDxHFECorrelationMC::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 
  if(RunFullMode())
    AliDebug(1, "Creating Full Corr THnSparse");
  else
    AliDebug(1, "Creating Reduced Corr THnSparse");
  // here is the only place to change the dimension
  static const int sizeEventdphi = 10;  
  static const int sizeEventdphiReduced = 5;  
  InitTHnSparseArray(sizeEventdphi);
  const double pi=TMath::Pi();
  Double_t minPhi=GetMinPhi();
  Double_t maxPhi=GetMaxPhi();
  THnSparse* thn=NULL;

  // TODO: Everything here needed for eventmixing? 
  TString name;
  name.Form("%s info", GetName());

  if(RunFullMode()){
    // 			                        0        1     2      3      4     5       6      7    8      9
    // 			                      D0invmass  PtD0 PhiD0 PtbinD0 Pte  dphi    dEta OrigD0 origEl process
    int         binsEventdphi[sizeEventdphi] = {   200,    100,  100,  21,   100, 100,     100,   10,     14,  100 };
    double      minEventdphi [sizeEventdphi] = { 1.5648,     0,    0,   0,     0 , minPhi, -0.9, -1.5,   -1.5, -0.5 };
    double      maxEventdphi [sizeEventdphi] = { 2.1648,    50, 2*pi,  20,    10, maxPhi,  0.9,  8.5,   12.5, 99.5 };
    const char* nameEventdphi[sizeEventdphi] = {
      "D0InvMass",
      "PtD0",
      "PhiD0",
      "PtBinD0",
      "PtEl",
      "#Delta#Phi",
      "#Delta#eta", 
      "Origin D0", 
      "Origin Electron",
      "Original Process"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphi,binsEventdphi,minEventdphi,maxEventdphi,nameEventdphi);
  }
  else{//Could also here consider removing process
    //           			                    0        1      2    3      4  
    // 		                	                 D0invmass  PtD0  Pte   dphi    dEta
    int         binsEventdphiRed[sizeEventdphiReduced] = {   200,    100,  21, 100,     100 };
    double      minEventdphiRed [sizeEventdphiReduced] = { 1.5648,     0,   0, minPhi, -0.9 };
    double      maxEventdphiRed [sizeEventdphiReduced] = { 2.1648,    50,  20, maxPhi,  0.9 };
    const char* nameEventdphiRed[sizeEventdphiReduced] = {
      "D0InvMass",
      "PtD0",
      "PtEl",
      "#Delta#Phi",
      "#Delta#eta" 
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiReduced,binsEventdphiRed,minEventdphiRed,maxEventdphiRed,nameEventdphiRed);

  }
  return thn;


}

Bool_t AliDxHFECorrelationMC::TestParticle(AliVParticle* p, Int_t id){

  AliReducedParticle *part=(AliReducedParticle*)p;
  if (!part) return -ENODATA;

  Bool_t selected = kTRUE;
  Bool_t isCharm=(part->GetOriginMother()==AliDxHFEToolsMC::kOriginCharm || 
		  part->GetOriginMother()==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(part->GetOriginMother()==AliDxHFEToolsMC::kOriginBeauty || 
		   part->GetOriginMother()==AliDxHFEToolsMC::kOriginGluonBeauty);

  if(id==kD){
    if(!fStoreOriginD==kAll){
      if(fStoreOriginD==kC ){
	// Test if particle really is from C
	if(!isCharm)
	  selected =kFALSE;
      }
      else if(fStoreOriginD==kB ){
	// Test if particle really is from B
	if(!isBeauty)
	  selected =kFALSE;
      }
      else if(fStoreOriginD==kHF ){
	// Test if particle really is HF
	if(!(isCharm || isBeauty))
	  selected=kFALSE;
      }
    }
  }

  // Test to see if test for D/el and whether to do further selection
  if(id==kElectron){
    if(!fStoreOriginEl==kAll){
      if(fStoreOriginEl==kC){ // Test if particle really is from C
	if(!isCharm)
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kB){ // Test if particle really is from B

	if(!isBeauty)
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kHF){ // Test if particle really is HF	
	if((!isCharm) || (!isBeauty))
	  selected=kFALSE;
      }
  
      // Test extra for source of el
      if(fStoreOriginEl==kNonHF){
	// Test if particle really is from NonHF
	if(!((part->GetOriginMother() >= AliDxHFEToolsMC::kOriginNone && part->GetOriginMother()<=AliDxHFEToolsMC::kOriginStrange)
	     || (part->GetOriginMother() > AliDxHFEToolsMC::kOriginGluonBeauty )))
	  selected =kFALSE;
      }
      else if(fStoreOriginEl==kHadrons){
	// Test if particle really is from hadrons
	if(! (part->GetOriginMother()<AliDxHFEToolsMC::kOriginNone))
	  selected =kFALSE;
      }
    }
  }

  return selected;

}

int AliDxHFECorrelationMC::FillParticleProperties(AliVParticle* tr, AliVParticle *as, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliReducedParticle *ptrigger=(AliReducedParticle*)tr;
  AliReducedParticle *assoc=(AliReducedParticle*)as;
  if (!ptrigger || !assoc) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
    data[i++]=ptrigger->GetInvMass();
    data[i++]=ptrigger->Pt();
    if(RunFullMode())data[i++]=ptrigger->Phi();
    if(RunFullMode()) {data[i++]=ptrigger->GetPtBin(); }
    data[i++]=assoc->Pt();
  } 
  else{
    data[i++]=assoc->GetInvMass();
    data[i++]=assoc->Pt();
    data[i++]=assoc->Phi();
    if(RunFullMode()) data[i++]=assoc->GetPtBin(); 
    data[i++]=ptrigger->Pt();
  }
  data[i++]=AliDxHFECorrelation::GetDeltaPhi();
  data[i++]=AliDxHFECorrelation::GetDeltaEta();
  if(RunFullMode()){
    if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
      data[i++]=ptrigger->GetOriginMother();
      data[i++]=assoc->GetOriginMother();
    }
    else {
      data[i++]=assoc->GetOriginMother();
      data[i++]=ptrigger->GetOriginMother();
    }
  }
  if(RunFullMode()) data[i++]=fMCEventType;
  return i;
}

int AliDxHFECorrelationMC::Fill(const TObjArray* triggerCandidates, TObjArray* associatedTracks, const AliVEvent* pEvent)
{
  // TODO: Implement more on MC?? (Needed?)
  return AliDxHFECorrelation::Fill(triggerCandidates,associatedTracks,pEvent);
}

double AliDxHFECorrelationMC::GetD0Eff(AliVParticle* tr){

  Double_t D0eff=1;
  AliReducedParticle *track=(AliReducedParticle*)tr;
  if (!track) return -ENODATA;
  Double_t pt=track->Pt();
  Double_t origin=track->GetOriginMother();

  Bool_t isCharm=(origin==AliDxHFEToolsMC::kOriginCharm || 
		  origin==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(origin==AliDxHFEToolsMC::kOriginBeauty || 
		   origin==AliDxHFEToolsMC::kOriginGluonBeauty);

  TH1F* effMap=NULL;
  // TODO: redefine how D0 origin is set, so for now only say that on uses 
  // Feeddown correction when it's defined as beauty, for the rest apply
  // Prompt correction
  if(isBeauty)
    effMap=(TH1F*)fD0EffMapFD;
  else
    effMap=(TH1F*)fD0EffMapP;

  if(isCharm)  AliDebug(2, "Correcting for Prompt D0");
  else AliDebug(2, "Correcting for Feeddown D0");

  Int_t bin=effMap->FindBin(pt);
  if(effMap->IsBinUnderflow(bin)|| effMap->IsBinOverflow(bin)) 
    D0eff = 1.;
  else 
    D0eff = effMap->GetBinContent(bin);
  return D0eff;
}

void AliDxHFECorrelationMC::SetD0EffMap(TH1* eff, int whichMap){

  if(whichMap==AliDxHFECorrelation::kPrompt) {fD0EffMapP=(TH1F*)eff; }
  if(whichMap==AliDxHFECorrelation::kFeedDown) { fD0EffMapFD=(TH1F*)eff;}
}


int AliDxHFECorrelationMC::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if (argument.BeginsWith("storeoriginD=")){
      argument.ReplaceAll("storeoriginD=", "");
      if (argument.CompareTo("all")==0) { fStoreOriginD=kAll; AliInfo("Store all Ds"); }
      else if (argument.CompareTo("charm")==0) { fStoreOriginD=kC; AliInfo("Store only D from charm");}
      else if (argument.CompareTo("beauty")==0){ fStoreOriginD=kB; AliInfo("Store only D from beauty");}
      else if (argument.CompareTo("HF")==0){ fStoreOriginD=kHF; AliInfo("Store D from HF");}
      continue;
    }  
    if (argument.BeginsWith("storeoriginEl=")){
      argument.ReplaceAll("storeoriginEl=", "");
      if (argument.CompareTo("all")==0) { fStoreOriginEl=kAll; AliInfo("Store all electrons"); }
      else if (argument.CompareTo("charm")==0) { fStoreOriginEl=kC; AliInfo("Store only electrons from charm");}
      else if (argument.CompareTo("beauty")==0){ fStoreOriginEl=kB; AliInfo("Store only electrons from beauty");}
      else if (argument.CompareTo("HF")==0){ fStoreOriginEl=kHF; AliInfo("Store electrons from HF");}
      else if (argument.CompareTo("nonHF")==0){ fStoreOriginEl=kNonHF; AliInfo("Store electrons from nonHF");}
      else if (argument.CompareTo("hadrons")==0){ fStoreOriginEl=kHadrons; AliInfo("Store electrons candidates from hadrons");}
      continue;
    }  

    //    AliWarning(Form("unknown argument '%s'", argument.Data()));
    AliDxHFECorrelation::ParseArguments(argument);      
  }

  return 0;
}
