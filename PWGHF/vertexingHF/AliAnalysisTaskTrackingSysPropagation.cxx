/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
// AliAnalysisTask for Tracking Systematics (matching efficiency          //
// + tracking efficiency) porpagation at the D-meson level                //
// (including daughter's kinematics)                                      //
////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include <Rtypes.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskTrackingSysPropagation.h"
#include "AliMultiplicity.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"


ClassImp(AliAnalysisTaskTrackingSysPropagation)
/* $Id$ */

//________________________________________________________________________
AliAnalysisTaskTrackingSysPropagation::AliAnalysisTaskTrackingSysPropagation():
AliAnalysisTaskSE("taskTrackingSysProp"),
  fPartName(""),
  fOutput(0x0),
  fAnalysisCuts(0),
  fHistNEvents(0x0),
  fHistMESyst(0x0),
  fHistTrEffSyst(0x0),
  fhPtDauVsD(0x0),
  fhSystMatchEffD(0x0),
  fDecayChannel(AliAnalysisTaskTrackingSysPropagation::kDplustoKpipi),
  fPDGcode(411),
  fAODProtection(1),
  fMaxPt(60.)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskTrackingSysPropagation::AliAnalysisTaskTrackingSysPropagation(AliAnalysisTaskTrackingSysPropagation::DecChannel ch, AliRDHFCuts* cuts, TH1F *HistMESys, TH1F *HistTrEffSys):
  AliAnalysisTaskSE("taskTrackingSysProp"),
  fPartName(""),
  fOutput(0x0),
  fAnalysisCuts(cuts),
  fHistNEvents(0x0),
  fHistMESyst(0x0),
  fHistTrEffSyst(0x0),
  fhPtDauVsD(0x0),
  fhSystMatchEffD(0x0),
  fDecayChannel(ch),
  fPDGcode(411),
  fAODProtection(1),
  fMaxPt(60.)
{
  fHistMESyst = new TH1F(*(static_cast<TH1F*>(HistMESys)));
  fHistTrEffSyst = new TH1F(*(static_cast<TH1F*>(HistTrEffSys)));
    
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskTrackingSysPropagation::~AliAnalysisTaskTrackingSysPropagation(){
    
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistMESyst;
    delete fHistTrEffSyst;
    delete fhPtDauVsD;
    delete fhSystMatchEffD;
  }        
  delete fOutput;
  delete fAnalysisCuts;
  fOutput = 0;
}

//________________________________________________________________________
void AliAnalysisTaskTrackingSysPropagation::UserCreateOutputObjects(){
    
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
    
  fHistNEvents = new TH1F("hNEvents", "number of events ",15,-0.5,14.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents Mismatched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"no. of 3 prong candidates");
  fHistNEvents->GetXaxis()->SetBinLabel(13,"no. of D after filtering cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(14,"no. of D after selection cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(15,"no. of not on-the-fly rec D");
    
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
    
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
   
  int nbins = (int)(10*fMaxPt);
  fhSystMatchEffD = new TH2F("fhSystMatchEffD","Matching Efficiency; p_{T} D; syst. (\%)",nbins,0.,fMaxPt,50,0.,25.);
  fhPtDauVsD      = new TH2F("fhPtDauVsD","Pt Dau vs D; p_{T} D; p_{T} daugh",nbins,0.,fMaxPt,nbins,0.,fMaxPt);
    
  fOutput->Add(fHistMESyst);
  fOutput->Add(fHistTrEffSyst);
  fOutput->Add(fhSystMatchEffD);
  fOutput->Add(fhPtDauVsD);
    
  PostData(1,fOutput);
}
//________________________________________________________________________
void AliAnalysisTaskTrackingSysPropagation::Init(){
    
  if(fDecayChannel == kDplustoKpipi) {
    fAnalysisCuts = new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kDstoKKpi) {
    fAnalysisCuts = new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kD0toKpi) {
    fAnalysisCuts = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kDstartoKpipi) {
    fAnalysisCuts = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fAnalysisCuts)));
  }
  else {
    AliFatal("The decay channel MUST be defined according to AliCFVertexing::DecayChannel - Exiting...");
  }
}
//________________________________________________________________________
void AliAnalysisTaskTrackingSysPropagation::UserExec(Option_t *){
    
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(fInputEvent);
  fHistNEvents->Fill(0); // all events
    
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(2);
      PostData(1,fOutput);
      return;
    }
  }
  TClonesArray *arrayBranch=0;
    
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi)
	arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      else if(fDecayChannel == kD0toKpi)
	arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      else if(fDecayChannel == kDstartoKpipi)
	arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  }
  else {
    if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi)
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    else if(fDecayChannel == kD0toKpi)
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    else if(fDecayChannel == kDstartoKpipi)
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("Dstar");
  }
    
  if (!arrayBranch) {
    AliError("Could not find array of HF vertices");
    PostData(1,fOutput);
    return;
  }
    
  AliAODVertex *aodVtx = (AliAODVertex*)aod->GetPrimaryVertex();
  if (!aodVtx || TMath::Abs(aod->GetMagneticField())<0.001) {
    AliDebug(3, "The event was skipped due to missing vertex or magnetic field issue");
    PostData(1,fOutput);
    return;
  }
  fHistNEvents->Fill(3); // count event
  if(fDecayChannel == kDplustoKpipi){
    fPDGcode = 411;
  }else if(fDecayChannel == kDstoKKpi){
    fPDGcode = 431;
  }else if(fDecayChannel == kD0toKpi){
    fPDGcode = 421;
  }else if(fDecayChannel == kDstartoKpipi){
    fPDGcode = 413;
  }
  Bool_t isEvSel  = fAnalysisCuts->IsEventSelected(aod);
  Float_t ntracks = aod->GetNumberOfTracks();
    
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()) fHistNEvents->Fill(5);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) fHistNEvents->Fill(6);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) fHistNEvents->Fill(7);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) fHistNEvents->Fill(8);
  if(fAnalysisCuts->IsEventRejectedDueToPileup()) fHistNEvents->Fill(9);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()) fHistNEvents->Fill(10);
    
  Int_t runNumber = aod->GetRunNumber();
    
  TClonesArray *arrayMC    =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  AliAODMCHeader *mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    
  if(!arrayMC) {
    AliError("AliAnalysisTaskTrackingSysPropagation::UserExec: MC particles branch not found!\n");
    PostData(1,fOutput);
    return;
  }
  if(!mcHeader) {
    AliError("AliAnalysisTaskTrackingSysPropagation::UserExec: MC header branch not found!\n");
    PostData(1,fOutput);
    return;
  }
    
  // arrayMC->Print();
  if(aod->GetTriggerMask()==0 && (runNumber>=195344 && runNumber<=195677)) {
    // protection for events with empty trigger mask in p-Pb
    PostData(1,fOutput);
    return;
  }
  if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0) {
    // events not passing the centrality selection can be removed immediately.
    PostData(1,fOutput);
    return;
  }
  Double_t zMCVertex = mcHeader->GetVtxZ();
  if (TMath::Abs(zMCVertex) > fAnalysisCuts->GetMaxVtxZ()) {
    PostData(1,fOutput);
    return;
  }
  if(!isEvSel){
    PostData(1,fOutput);
    return;
  }
  fHistNEvents->Fill(4);
    
  Int_t nCand = arrayBranch->GetEntriesFast();
    
  Int_t nprongs = -1;
  Int_t pdgDaughter[3];
  Int_t pdg2Daughter[2];
  if(fDecayChannel == kDplustoKpipi){
    fPDGcode = 411;
    nprongs  = 3;
    fPartName="Dplus";
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
    pdgDaughter[2]=211;
  }else if(fDecayChannel == kDstoKKpi){
    fPDGcode = 431;
    nprongs  = 3;
    fPartName="Ds";
    pdgDaughter[0]=321;
    pdgDaughter[1]=321;
    pdgDaughter[2]=211;
  }else if(fDecayChannel == kD0toKpi){
    fPDGcode = 421;
    nprongs  = 2;
    fPartName="D0";
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
  }else if(fDecayChannel == kDstartoKpipi){
    fPDGcode = 413;
    nprongs  = 2;
    fPartName="Dstar";
    pdgDaughter[0]=421;
    pdgDaughter[1]=211;
    pdg2Daughter[0]=321;
    pdg2Daughter[1]=211;
  }else{
    AliError("WRONG DECAY SETTING");
    PostData(1,fOutput);
    return;
  }
    
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
    
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
        
    AliAODRecoDecayHF* d = 0x0;
        
    if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi) {
      d = (AliAODRecoDecayHF3Prong*)arrayBranch->UncheckedAt(iCand);
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d)) {
	fHistNEvents->Fill(14);
	continue;
      }
    }
    else if(fDecayChannel == kD0toKpi) {
      d = (AliAODRecoDecayHF2Prong*)arrayBranch->UncheckedAt(iCand);
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d)) {
	fHistNEvents->Fill(14);
	continue;
      }
    }
    else if(fDecayChannel == kDstartoKpipi) {
      d = (AliAODRecoCascadeHF*)arrayBranch->At(iCand);
      if(!d) continue;
      Bool_t isDStarCand =kTRUE;
      if(!vHF->FillRecoCasc(aod,(AliAODRecoCascadeHF*)d,isDStarCand)) {
	fHistNEvents->Fill(14);
	continue;
      }
      if(!d->GetSecondaryVtx()) continue;
    }
        
    fHistNEvents->Fill(11);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(aodVtx);
      unsetvtx=kTRUE;
    }
        
    Bool_t recVtx=kFALSE;
    AliAODVertex *origownvtx=0x0;
        
    Double_t ptD = d->Pt();
    Double_t rapid  = d->Y(fPDGcode);
    Bool_t isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(ptD,rapid);
        
    if(isFidAcc){
      Int_t retCodeAnalysisCuts = fAnalysisCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
      if(retCodeAnalysisCuts > 0) {
	fHistNEvents->Fill(13);
	if(fDecayChannel == kDstoKKpi){
	  Int_t isPhiKKpi = retCodeAnalysisCuts&4; //for Ds
	  Int_t isPhipiKK = retCodeAnalysisCuts&8; //for Ds
	  if(!(isPhiKKpi || isPhipiKK)) continue;
	}
                
	double syst = 0.;
	Int_t nDau = d->GetNDaughters();
	AliAODTrack *track;
                
	Int_t mcLabel=-1;
	int nbinsME = fHistMESyst->GetNbinsX();
	int bin = 0;
	if(!(fDecayChannel == kDstartoKpipi)) {
	  mcLabel = d->MatchToMC(fPDGcode,arrayMC,nprongs,pdgDaughter);
	  if (mcLabel < 0) continue;
	  for(Int_t iDau=0; iDau < nDau; iDau++){
	    track = (AliAODTrack*)d->GetDaughter(iDau);
	    if(!track){
	      AliError("Daughter particle track not found");
	      PostData(1,fOutput);
	      return;
	    }
	    Int_t labDau = track->GetLabel();
	    AliAODMCParticle* p = (AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau));
	    double ptdau = track->Pt();
               
	    bin = fHistMESyst->FindBin(ptdau);
	    if(bin > 0 && bin <= nbinsME) syst += fHistMESyst->GetBinContent(bin);
	    else if(bin > nbinsME) syst += fHistMESyst->GetBinContent(nbinsME);
	    else if(bin == 0) {
	      AliError("Check input histo at low pt!");
	      PostData(1,fOutput);
	      return;
	    }
	    fhPtDauVsD->Fill(ptD,ptdau);
	  }
	}
	else {  //D*
	  mcLabel = ((AliAODRecoCascadeHF*)d)->MatchToMC(413,421,pdgDaughter,pdg2Daughter,arrayMC,kFALSE);
	  if (mcLabel < 0) continue;
                    
	  //soft pion
	  track = (AliAODTrack*)(((AliAODRecoCascadeHF*)d)->GetBachelor());
	  if(!track) {
	    PostData(1,fOutput);
	    return;
	  }
	  Int_t labDau = track->GetLabel();
	  AliAODMCParticle* p = (AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau));
	  double ptdau = track->Pt();
                   
	  bin = fHistMESyst->FindBin(ptdau);
	  if(bin > 0 && bin <= nbinsME) syst += fHistMESyst->GetBinContent(bin);
	  else if(bin > nbinsME) syst += fHistMESyst->GetBinContent(nbinsME);
	  else if(bin == 0) {
	    AliError("Check input histo at low pt!");
	    PostData(1,fOutput);
	    return;
	  }
	  fhPtDauVsD->Fill(ptD,ptdau);
	  //D0
	  AliAODRecoDecayHF2Prong *D0 = ((AliAODRecoCascadeHF*)d)->Get2Prong();
	  for(Int_t iDau=0; iDau < 2; iDau++){
	    track = 0x0;
	    track = (AliAODTrack*)D0->GetDaughter(iDau);
	    if(!track) {
	      PostData(1,fOutput);
	      return;
	    }
	    labDau = track->GetLabel();
	    p=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau));
	    ptdau = track->Pt();
                        
	    bin = fHistMESyst->FindBin(ptdau);
	    if(bin > 0 && bin <= nbinsME) syst += fHistMESyst->GetBinContent(bin);
	    else if(bin > nbinsME) syst += fHistMESyst->GetBinContent(nbinsME);
	    else if(bin == 0) {
	      AliError("Check input histo at low pt!");
	      PostData(1,fOutput);
	      return;
	    }
	    fhPtDauVsD->Fill(ptD,ptdau);
	  }
	}
	int nbinsTE = fHistTrEffSyst->GetNbinsX();
	bin = fHistTrEffSyst->FindBin(ptD);
	double trackEffSys = 0.;
	if(bin > 0 && bin <= nbinsTE) trackEffSys = fHistTrEffSyst->GetBinContent(bin);
	else if(bin > nbinsTE) trackEffSys = fHistTrEffSyst->GetBinContent(nbinsTE);
	else if(bin == 0) {
	  AliError("Check input histo at low pt!");
	  PostData(1,fOutput);
	  return;
	}
	syst = TMath::Sqrt(syst*syst+trackEffSys*trackEffSys);
	fhSystMatchEffD->Fill(ptD,syst);
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  PostData(1,fOutput);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskTrackingSysPropagation::Terminate(Option_t *) {
    
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    AliError("ERROR: fOutput not available\n");
    return;
  }
  if(fHistNEvents){
    Printf("Number of Analyzed Events = %f",fHistNEvents->GetBinContent(1));
  }else{
    AliError("ERROR: fHistNEvents not available\n");
  }
    
  Printf("end of Terminate");
  return;
}
