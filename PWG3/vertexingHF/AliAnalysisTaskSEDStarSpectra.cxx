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
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//                  Base class for DStar Analysis
//
//
//  The D* spectra study is done in pt bins:
//  [0,0.5] [0.5,1] [1,2] [2,3] [3,4] [4,5] [5,6] [6,7] [7,8],
//  [8,10],[10,12], [12,16], [16,20] and [20,24]
//
//  Cuts arew centralized in AliRDHFCutsDStartoKpipi
//  Side Band and like sign background are implemented in the macro
//
//-----------------------------------------------------------------------
//
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl,
//                         Author Y.Wang
//        University of Heidelberg - yifei@physi.uni-heidelberg.de
//                         Author C.Ivan 
//             ERC-QGP Utrecht University - c.ivan@uu.nl,
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDStarSpectra.h"
#include "AliNormalizationCounter.h"

ClassImp(AliAnalysisTaskSEDStarSpectra)

//__________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra():  
  AliAnalysisTaskSE(),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE),
  fDoSearch(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPID(0),
  fNSigma(3),
  fCuts(0),
  fCEvents(0),     
  fTrueDiff2(0),
  fDeltaMassD1(0),
  fCounter(0),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const Char_t* name, AliRDHFCutsDStartoKpipi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fAnalysis(0),
  fD0Window(0),
  fPeakWindow(0),
  fUseMCInfo(kFALSE),
  fDoSearch(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPID(0),
  fNSigma(3),
  fCuts(0),
  fCEvents(0),     
  fTrueDiff2(0),
  fDeltaMassD1(0),
  fCounter(0),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarSpectra","Calling Constructor");

  fCuts=cuts;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());  //All Entries output
  DefineOutput(3,TList::Class());  //3sigma PID output
  DefineOutput(4,AliRDHFCutsDStartoKpipi::Class());   //My private output
  DefineOutput(5,AliNormalizationCounter::Class());   // normalization
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::~AliAnalysisTaskSEDStarSpectra() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEDStarSpectra","Calling Destructor");
  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }
  if (fOutputPID) {
    delete fOutputPID;
    fOutputPID = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if(fCEvents){
    delete fCEvents;
    fCEvents =0;
  }
  if(fDeltaMassD1){
    delete fDeltaMassD1;
    fDeltaMassD1 =0;
  }
  for(Int_t i=0; i<5; i++){
    delete fHistMassPtImpParTCDs[i];
  }
}
//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::Init(){
  //
  // Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");
   AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
  // Post the data
  PostData(4,copyfCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  fEvents++;

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi=0;

  fCEvents->Fill(1);

  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  } else {
    arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
  }

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
  fCEvents->Fill(2);

  fCounter->StoreEvent(aodEvent,fCuts,fUseMCInfo);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass=aodEvent->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  if(!fCuts->IsEventSelected(aodEvent)) {
    if(fCuts->GetWhyRejection()==6) // rejected for Z vertex
      fCEvents->Fill(6);
    return;
  }

  Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
  if(!isEvSel) return;

  // Load the event
  //  AliInfo(Form("Event %d",fEvents));
  //if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));

  // counters for efficiencies
  Int_t icountReco = 0;
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if(!vtx1) return;
  if(vtx1->GetNContributors()<1) return;

  if (!arrayDStartoD0pi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast())); 

  Int_t nSelectedAna =0;
  Int_t nSelectedProd =0;

  // loop over the tracks to search for candidates soft pion 
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {

    // D* candidates and D0 from D*
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    if(!dstarD0pi->GetSecondaryVtx()) continue;
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    if (!theD0particle) continue;
    
    Int_t isDStar = 0;   
    TClonesArray *mcArray = 0;
    AliAODMCHeader *mcHeader=0;

    Bool_t isPrimary=kTRUE;
    Float_t truePt=0.;
    Float_t xDecay=0.;
    Float_t yDecay=0.;
    Float_t zDecay=0.;
    Float_t pdgCode=-2;
    Float_t trueImpParXY=0.;

    // mc analysis 
    if(fUseMCInfo){
    //MC array need for maching
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
	AliError("Could not find Monte-Carlo in AOD");
	return;
      }
      // load MC header
      mcHeader =  (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if(!mcHeader) {
	printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
	return;
      }
      // find associated MC particle for D* ->D0toKpi
      Int_t mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      if(mcLabel>=0){

	AliAODMCParticle *partDSt = (AliAODMCParticle*)mcArray->At(mcLabel);
	Int_t checkOrigin = CheckOrigin(mcArray,partDSt);
       	if(checkOrigin==5) isPrimary=kFALSE;
	AliAODMCParticle *dg0 = (AliAODMCParticle*)mcArray->At(partDSt->GetDaughter(0));
	AliAODMCParticle *dg0_1 = (AliAODMCParticle*)mcArray->At(dg0->GetDaughter(0));
	truePt=dg0->Pt();
	xDecay=dg0_1->Xv();	  
	yDecay=dg0_1->Yv();	  
	zDecay=dg0_1->Zv();
	pdgCode=TMath::Abs(partDSt->GetPdgCode());
	if(!isPrimary){
	  trueImpParXY=GetTrueImpactParameterD0(mcHeader,mcArray,dg0)*1000.;
	}
	isDStar = 1;
      }else{
	pdgCode=-1;
      }
    }
    
    Int_t ptbin=fCuts->PtBin(dstarD0pi->Pt());
    
    // quality selction on tracks and region of interest
    Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
    if(!isTkSelected) continue;

    if(!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(),dstarD0pi->YDstar())) continue;    


    //histos for impact par studies - D0!!!
    Double_t ptCand = dstarD0pi->Get2Prong()->Pt();
    Double_t invMass=dstarD0pi->InvMassD0();
    Double_t impparXY=dstarD0pi->Get2Prong()->ImpParXY()*10000.;

    Double_t arrayForSparse[3]={invMass,ptCand,impparXY};
    Double_t arrayForSparseTrue[3]={invMass,ptCand,trueImpParXY};
   
  // set the D0 search window bin by bin - useful to calculate side band bkg
    if (ptbin==0){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==1){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==2){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (ptbin==3){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.022;
	fPeakWindow=0.0016;
      }
    }
    if (ptbin==4){ 
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    }
    if (ptbin==5){
      if(fAnalysis==1){
	fD0Window=0.045;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    } 
   if (ptbin==6){
      if(fAnalysis==1){
	fD0Window=0.045;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.006;
      }
    } 
    if (ptbin==7){
      if(fAnalysis==1){
	fD0Window=0.055;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.006;
      }
    }
    if (ptbin>7){
      if(fAnalysis==1){
	fD0Window=0.074;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.006;
      }
    }
    
    nSelectedProd++;
    nSelectedAna++;
    
    // check that we are close to signal in the DeltaM - here to save time for PbPb
    Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();  
    Double_t mPDGDstar=TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t invmassDelta = dstarD0pi->DeltaInvMass();
    
    if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))>fPeakWindow) continue;
    
    
    Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
    
    // after cuts
    if(fDoImpParDstar && isSelected){
      fHistMassPtImpParTCDs[0]->Fill(arrayForSparse);
      if(isPrimary) fHistMassPtImpParTCDs[1]->Fill(arrayForSparse);
      else{
	fHistMassPtImpParTCDs[2]->Fill(arrayForSparse);
	fHistMassPtImpParTCDs[3]->Fill(arrayForSparseTrue);
      }
    }
    
    // fill PID
    FillSpectrum(dstarD0pi,isDStar,fCuts,isSelected,fOutputPID);
    SideBandBackground(dstarD0pi,fCuts,isSelected, fOutputPID);
    //WrongSignForDStar(dstarD0pi,fCuts,fOutputPID);

    //swich off the PID selection
    fCuts->SetUsePID(kFALSE);
    Int_t isSelectedNoPID=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
    fCuts->SetUsePID(kTRUE);

    FillSpectrum(dstarD0pi,isDStar,fCuts,isSelectedNoPID,fOutputPID);
    SideBandBackground(dstarD0pi,fCuts,isSelectedNoPID, fOutputPID);

    // rare D search ------ 
    if(fDoSearch){
      TLorentzVector LorentzTrack1(0,0,0,0); // lorentz 4 vector
      TLorentzVector LorentzTrack2(0,0,0,0); // lorentz 4 vector
      
      for (Int_t i=0; i<aodEvent->GetNTracks(); i++){ 
	
	AliAODTrack* aodTrack = aodEvent->GetTrack(i);
	
	if(dstarD0pi->Charge() == aodTrack->Charge()) continue;
	if((!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
	if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))>0.02) continue;
	
	//build the D1 mass
	Double_t mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();

	LorentzTrack1.SetPxPyPzE( dstarD0pi->Px(),dstarD0pi->Py(), dstarD0pi->Pz(), dstarD0pi->E(413) );
	LorentzTrack2.SetPxPyPzE( aodTrack->Px(),aodTrack->Py(), aodTrack->Pz(),aodTrack->E(mass) );
	
	//D1 mass
	Double_t d1mass = ((LorentzTrack1+LorentzTrack2).M());
	//mass difference - at 0.4117 and 0.4566
	fDeltaMassD1->Fill(d1mass-dstarD0pi->InvMassDstarKpipi());
      }
    }

    if(isDStar == 1) { 
      fTrueDiff2->Fill(dstarD0pi->Pt(),dstarD0pi->DeltaInvMass());
    }
    
  }
  
  fCounter->StoreCandidates(aodEvent,nSelectedProd,kTRUE);  
  fCounter->StoreCandidates(aodEvent,nSelectedAna,kFALSE); 

  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));
  
  PostData(1,fOutput);
  PostData(2,fOutputAll);
  PostData(3,fOutputPID);
  PostData(5,fCounter);

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarSpectra::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fDeltaMassD1     = dynamic_cast<TH1F*>(fOutput->FindObject("fDeltaMassD1"));
  fTrueDiff2      = dynamic_cast<TH2F*>(fOutput->FindObject("fTrueDiff2"));

  fOutputAll = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputAll) {
    printf("ERROR: fOutputAll not available\n");
    return;
  }
  fOutputPID = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputPID) {
    printf("ERROR: fOutputPID not available\n");
    return;
  }
 

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserCreateOutputObjects() { 
 // output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("listAll");

  fOutputPID = new TList();
  fOutputPID->SetOwner();
  fOutputPID->SetName("listPID");
    
  // define histograms
  DefineHistograms();

 //Counter for Normalization
 fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
 fCounter->Init();

 if(fDoImpParDstar) CreateImpactParameterHistos();

  PostData(1,fOutput);
  PostData(2,fOutputAll);
  PostData(3,fOutputPID);

  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSEDStarSpectra::DefineHistograms(){

  fCEvents = new TH1F("fCEvents","conter",11,0,11);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");
  fOutput->Add(fCEvents);

  fTrueDiff2 = new TH2F("DiffDstar_pt","True Reco diff vs pt",200,0,15,900,0,0.3);
  fOutput->Add(fTrueDiff2);

  fDeltaMassD1 = new TH1F("DeltaMassD1","delta mass d1",600,0,0.8);
  fOutput->Add(fDeltaMassD1);

  const Int_t nhist=14;
  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  for(Int_t i=-2;i<nhist;i++){
    nameMass="histDeltaMass_";
    nameMass+=i+1;
    nameSgn="histDeltaSgn_";
    nameSgn+=i+1;
    nameBkg="histDeltaBkg_";
    nameBkg+=i+1; 
    
    if (i==-2) {
      nameMass="histDeltaMass";
      nameSgn="histDeltaSgn";
      nameBkg="histDeltaBkg";
    }
    
    TH1F* spectrumMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} invariant mass; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    TH1F* spectrumSgn = new TH1F(nameSgn.Data(), "D^{*}-D^{0} Signal invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    TH1F* spectrumBkg = new TH1F(nameBkg.Data(), "D^{*}-D^{0} Background invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",200,0.1,0.2);
    
    nameMass="histD0Mass_";
    nameMass+=i+1;
    nameSgn="histD0Sgn_";
    nameSgn+=i+1;
    nameBkg="histD0Bkg_";
    nameBkg+=i+1;
    
    if (i==-2) {
      nameMass="histD0Mass";
      nameSgn="histD0Sgn";
      nameBkg="histD0Bkg";
    }

    TH1F* spectrumD0Mass = new TH1F(nameMass.Data(),"D^{0} invariant mass; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);
    TH1F* spectrumD0Sgn = new TH1F(nameSgn.Data(), "D^{0} Signal invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);
    TH1F* spectrumD0Bkg = new TH1F(nameBkg.Data(), "D^{0} Background invariant mass - MC; M(D^{0}) [GeV/c^{2}]; Entries",200,1.75,1.95);

    nameMass="histDstarMass_";
    nameMass+=i+1;
    nameSgn="histDstarSgn_";
    nameSgn+=i+1;
    nameBkg="histDstarBkg_";
    nameBkg+=i+1;

    if (i==-2) {
      nameMass="histDstarMass";
      nameSgn="histDstarSgn";
      nameBkg="histDstarBkg";
    }

    TH1F* spectrumDstarMass = new TH1F(nameMass.Data(),"D^{*} invariant mass; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarSgn = new TH1F(nameSgn.Data(), "D^{*} Signal invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);
    TH1F* spectrumDstarBkg = new TH1F(nameBkg.Data(), "D^{*} Background invariant mass - MC; M(D^{*}) [GeV/c^{2}]; Entries",200,1.9,2.1);

    nameMass="histSideBandMass_";
    nameMass+=i+1;
    if (i==-2) { 
      nameMass="histSideBandMass";
    }
    
    TH1F* spectrumSideBandMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} sideband mass; M(D^{*}) [GeV/c^{2}]; Entries",200,0.1,0.2);

    nameMass="histWrongSignMass_";
    nameMass+=i+1;
    if (i==-2) { 
      nameMass="histWrongSignMass";
    }
    
    TH1F* spectrumWrongSignMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} wrongsign mass; M(D^{*}) [GeV/c^{2}]; Entries",200,0.1,0.2);


    spectrumMass->Sumw2();
    spectrumSgn->Sumw2();
    spectrumBkg->Sumw2();
    
    spectrumMass->SetLineColor(6);
    spectrumSgn->SetLineColor(2);
    spectrumBkg->SetLineColor(4);
    
    spectrumMass->SetMarkerStyle(20);
    spectrumSgn->SetMarkerStyle(20);
    spectrumBkg->SetMarkerStyle(20);
    spectrumMass->SetMarkerSize(0.6);
    spectrumSgn->SetMarkerSize(0.6);
    spectrumBkg->SetMarkerSize(0.6);
    spectrumMass->SetMarkerColor(6);
    spectrumSgn->SetMarkerColor(2);
    spectrumBkg->SetMarkerColor(4);

    spectrumD0Mass->Sumw2();
    spectrumD0Sgn->Sumw2();
    spectrumD0Bkg->Sumw2();

    spectrumD0Mass->SetLineColor(6);
    spectrumD0Sgn->SetLineColor(2);
    spectrumD0Bkg->SetLineColor(4);

    spectrumD0Mass->SetMarkerStyle(20);
    spectrumD0Sgn->SetMarkerStyle(20);
    spectrumD0Bkg->SetMarkerStyle(20);
    spectrumD0Mass->SetMarkerSize(0.6);
    spectrumD0Sgn->SetMarkerSize(0.6);
    spectrumD0Bkg->SetMarkerSize(0.6);
    spectrumD0Mass->SetMarkerColor(6);
    spectrumD0Sgn->SetMarkerColor(2);
    spectrumD0Bkg->SetMarkerColor(4);

    spectrumDstarMass->Sumw2();
    spectrumDstarSgn->Sumw2();
    spectrumDstarBkg->Sumw2();

    spectrumDstarMass->SetLineColor(6);
    spectrumDstarSgn->SetLineColor(2);
    spectrumDstarBkg->SetLineColor(4);

    spectrumDstarMass->SetMarkerStyle(20);
    spectrumDstarSgn->SetMarkerStyle(20);
    spectrumDstarBkg->SetMarkerStyle(20);
    spectrumDstarMass->SetMarkerSize(0.6);
    spectrumDstarSgn->SetMarkerSize(0.6);
    spectrumDstarBkg->SetMarkerSize(0.6);
    spectrumDstarMass->SetMarkerColor(6);
    spectrumDstarSgn->SetMarkerColor(2);
    spectrumDstarBkg->SetMarkerColor(4);

    spectrumSideBandMass->Sumw2();
    spectrumSideBandMass->SetLineColor(4);
    spectrumSideBandMass->SetMarkerStyle(20);
    spectrumSideBandMass->SetMarkerSize(0.6);
    spectrumSideBandMass->SetMarkerColor(4);

    spectrumWrongSignMass->Sumw2();
    spectrumWrongSignMass->SetLineColor(4);
    spectrumWrongSignMass->SetMarkerStyle(20);
    spectrumWrongSignMass->SetMarkerSize(0.6);
    spectrumWrongSignMass->SetMarkerColor(4);

    TH1F* allMass = (TH1F*)spectrumMass->Clone();
    TH1F* allSgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* allBkg  = (TH1F*)spectrumBkg->Clone();

    TH1F* pidMass = (TH1F*)spectrumMass->Clone();
    TH1F* pidSgn  = (TH1F*)spectrumSgn->Clone();
    TH1F* pidBkg  = (TH1F*)spectrumBkg->Clone();

    fOutputAll->Add(allMass);
    fOutputAll->Add(allSgn);
    fOutputAll->Add(allBkg);

    fOutputPID->Add(pidMass);
    fOutputPID->Add(pidSgn);
    fOutputPID->Add(pidBkg);

    TH1F* allD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* allD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* allD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pidD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pidD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pidD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    fOutputAll->Add(allD0Mass);
    fOutputAll->Add(allD0Sgn);
    fOutputAll->Add(allD0Bkg);

    fOutputPID->Add(pidD0Mass);
    fOutputPID->Add(pidD0Sgn);
    fOutputPID->Add(pidD0Bkg);
  
    TH1F* allDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* allDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* allDstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    TH1F* pidDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pidDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pidDstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    fOutputAll->Add(allDstarMass);
    fOutputAll->Add(allDstarSgn);
    fOutputAll->Add(allDstarBkg);

    fOutputPID->Add(pidDstarMass);
    fOutputPID->Add(pidDstarSgn);
    fOutputPID->Add(pidDstarBkg);

    TH1F* allSideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pidSideBandMass = (TH1F*)spectrumSideBandMass->Clone();

    fOutputAll->Add(allSideBandMass);
    fOutputPID->Add(pidSideBandMass);
   
    TH1F* allWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pidWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();

    fOutputAll->Add(allWrongSignMass);
    fOutputPID->Add(pidWrongSignMass);
   
  }
  
  // pt spectra
  nameMass="ptMass";
  nameSgn="ptSgn";
  nameBkg="ptBkg";
  
  TH1F* ptspectrumMass = new TH1F(nameMass.Data(),"D^{*} p_{T}; p_{T} [GeV]; Entries",200,0,10);
  TH1F* ptspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal p_{T} - MC; p_{T} [GeV]; Entries",200,0,10);
  TH1F* ptspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background p_{T} - MC; p_{T} [GeV]; Entries",200,0,10);
  
  ptspectrumMass->Sumw2();
  ptspectrumSgn->Sumw2();
  ptspectrumBkg->Sumw2();
  
  ptspectrumMass->SetLineColor(6);
  ptspectrumSgn->SetLineColor(2);
  ptspectrumBkg->SetLineColor(4);
  
  ptspectrumMass->SetMarkerStyle(20);
  ptspectrumSgn->SetMarkerStyle(20);
  ptspectrumBkg->SetMarkerStyle(20);
  ptspectrumMass->SetMarkerSize(0.6);
  ptspectrumSgn->SetMarkerSize(0.6);
  ptspectrumBkg->SetMarkerSize(0.6);
  ptspectrumMass->SetMarkerColor(6);
  ptspectrumSgn->SetMarkerColor(2);
  ptspectrumBkg->SetMarkerColor(4);
  
  TH1F* ptallMass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptallSgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptallBkg = (TH1F*)ptspectrumBkg->Clone();
  
  TH1F* ptpidMass = (TH1F*)ptspectrumMass->Clone();
  TH1F* ptpidSgn = (TH1F*)ptspectrumSgn->Clone();
  TH1F* ptpidBkg = (TH1F*)ptspectrumBkg->Clone();
    
  fOutputAll->Add(ptallMass);
  fOutputAll->Add(ptallSgn);
  fOutputAll->Add(ptallBkg);
  
  fOutputPID->Add(ptpidMass);
  fOutputPID->Add(ptpidSgn);
  fOutputPID->Add(ptpidBkg);
  
  // eta spectra
  nameMass="etaMass";
  nameSgn="etaSgn";
  nameBkg="etaBkg";
  
  TH1F* etaspectrumMass = new TH1F(nameMass.Data(),"D^{*} #eta; #eta; Entries",200,-1,1);
  TH1F* etaspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal #eta - MC; #eta; Entries",200,-1,1);
  TH1F* etaspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background #eta - MC; #eta; Entries",200,-1,1);
  
  etaspectrumMass->Sumw2();
  etaspectrumSgn->Sumw2();
  etaspectrumBkg->Sumw2();
  
  etaspectrumMass->SetLineColor(6);
  etaspectrumSgn->SetLineColor(2);
  etaspectrumBkg->SetLineColor(4);
  
  etaspectrumMass->SetMarkerStyle(20);
  etaspectrumSgn->SetMarkerStyle(20);
  etaspectrumBkg->SetMarkerStyle(20);
  etaspectrumMass->SetMarkerSize(0.6);
  etaspectrumSgn->SetMarkerSize(0.6);
  etaspectrumBkg->SetMarkerSize(0.6);
  etaspectrumMass->SetMarkerColor(6);
  etaspectrumSgn->SetMarkerColor(2);
  etaspectrumBkg->SetMarkerColor(4);
  
  TH1F* etaallMass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etaallSgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etaallBkg = (TH1F*)etaspectrumBkg->Clone();
  
  TH1F* etapidMass = (TH1F*)etaspectrumMass->Clone();
  TH1F* etapidSgn = (TH1F*)etaspectrumSgn->Clone();
  TH1F* etapidBkg = (TH1F*)etaspectrumBkg->Clone();
    
  fOutputAll->Add(etaallMass);
  fOutputAll->Add(etaallSgn);
  fOutputAll->Add(etaallBkg);
  
  fOutputPID->Add(etapidMass);
  fOutputPID->Add(etapidSgn);
  fOutputPID->Add(etapidBkg);
  
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, AliRDHFCutsDStartoKpipi *cuts,Int_t isSel, TList *listout){
  //
  // Fill histos for D* spectrum
  //
  
  if(!isSel) return;

  // D0 window
  Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();  
  if (TMath::Abs(invmassD0-mPDGD0)>fD0Window) return; 


  Int_t ptbin=cuts->PtBin(part->Pt());  
  Double_t pt = part->Pt();
  Double_t eta = part->Eta();  
  
  Double_t invmassDelta = part->DeltaInvMass();
  Double_t invmassDstar = part->InvMassDstarKpipi();
  
  TString fillthis="";
  Bool_t massInRange=kFALSE;
  
  Double_t mPDGDstar=TDatabasePDG::Instance()->GetParticle(413)->Mass();
  
  // delta M(Kpipi)-M(Kpi)
  if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<fPeakWindow) massInRange=kTRUE;  
  
  if(fUseMCInfo) {
    if(isDStar==1) {
      fillthis="histD0Sgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histD0Sgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histDstarSgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDstarSgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDeltaSgn_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      fillthis="histDeltaSgn";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
    if (massInRange) {
	fillthis="ptSgn";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
	fillthis="etaSgn";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
      }
    }
    else {//background
      fillthis="histD0Bkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histD0Bkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      fillthis="histDstarBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDstarBkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
      fillthis="histDeltaBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
      fillthis="histDeltaBkg";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
     if (massInRange) {
	fillthis="ptBkg";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
	fillthis="etaBkg";
	((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
      }
    }
  }
  //no MC info, just cut selection
  fillthis="histD0Mass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
  fillthis="histD0Mass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
  fillthis="histDstarMass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
  fillthis="histDstarMass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDstar);
  fillthis="histDeltaMass_";
  fillthis+=ptbin;
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
  fillthis="histDeltaMass";
  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
  
  if (massInRange) {
    fillthis="ptMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(pt);
    fillthis="etaMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(eta);
  }
 
  return;
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarSpectra::SideBandBackground(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   

  if(!isSel) return;

  Int_t ptbin=cuts->PtBin(part->Pt());
  
  Bool_t massInRange=kFALSE;
  
  // select the side bands intervall
  Double_t invmassD0    = part->InvMassD0();
  if(TMath::Abs(invmassD0-1.865)>4*fD0Window && TMath::Abs(invmassD0-1.865)<8*fD0Window){
    
    // for pt and eta
    Double_t invmassDelta = part->DeltaInvMass();
    if (TMath::Abs(invmassDelta-0.14557)<fPeakWindow) massInRange=kTRUE;
    
    TString fillthis="";
    fillthis="histSideBandMass_";
    fillthis+=ptbin;
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
    fillthis="histSideBandMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassDelta);
    
  }
}
//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::WrongSignForDStar(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, TList *listout){
  //
  // assign the wrong charge to the soft pion to create background
  //
  Int_t ptbin=cuts->PtBin(part->Pt());
  
  Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();  
  if (TMath::Abs(invmassD0-mPDGD0)>fD0Window) return; 

  AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)part->Get2Prong();

  Int_t okD0WrongSign,okD0barWrongSign;
  Double_t wrongMassD0=0.;
  
  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
   if (!isSelected){
    return;
  }

  okD0WrongSign =  1;
  okD0barWrongSign = 1;
  
  //if is D*+ than assume D0bar
  if(part->Charge()>0 && (isSelected ==1)) { 
    okD0WrongSign = 0;
  }
  if(part->Charge()<0 && (isSelected ==2)){
    okD0barWrongSign = 0;
  }
  
  // assign the wrong mass in case the cuts return both D0 and D0bar
  if(part->Charge()>0 && (isSelected ==3)) { 
    okD0WrongSign = 0;
  } else if(part->Charge()<0 && (isSelected ==3)){
    okD0barWrongSign = 0;
  }
  
  //wrong D0 inv mass
  if(okD0WrongSign!=0){
    wrongMassD0 = theD0particle->InvMassD0();
  }else if(okD0WrongSign==0){
    wrongMassD0 = theD0particle->InvMassD0bar();
  }
  
  if(TMath::Abs(wrongMassD0-1.865)<fD0Window){
    
    // wrong D* inv mass   
    Double_t e[3];
    if (part->Charge()>0){
      e[0]=theD0particle->EProng(0,321);
      e[1]=theD0particle->EProng(1,211);
    }else{
      e[0]=theD0particle->EProng(0,211);
      e[1]=theD0particle->EProng(1,321);
    }
    e[2]=part->EProng(0,211);
    
    Double_t esum = e[0]+e[1]+e[2];
    Double_t pds = part->P();

    Double_t   wrongMassDstar = TMath::Sqrt(esum*esum-pds*pds);

    TString fillthis="";
    fillthis="histWrongSignMass_";
    fillthis+=ptbin;
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
    fillthis="histWrongSignMass";
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(wrongMassDstar-wrongMassD0);
    
  }
}
//-------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEDStarSpectra::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isFromB) return 5;
  else return 4;
}
//-------------------------------------------------------------------------------------
Float_t AliAnalysisTaskSEDStarSpectra::GetTrueImpactParameterD0(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp) const {
  // true impact parameter calculation

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partDp->XvYvZv(origD);	  
  Short_t charge=partDp->Charge();
  Double_t pXdauTrue[3],pYdauTrue[3],pZdauTrue[3];
  Int_t labelFirstDau = partDp->GetDaughter(0); 

  Int_t nDau=partDp->GetNDaughters();

  Int_t theDau=0;
  if(nDau==2){
    for(Int_t iDau=0; iDau<2; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	AliError("Daughter particle not found in MC array");
	return 99999.;
      } 
      Int_t pdgCode=TMath::Abs(part->GetPdgCode());
      if(pdgCode==211 || pdgCode==321){
	pXdauTrue[theDau]=part->Px();
	pYdauTrue[theDau]=part->Py();
	pZdauTrue[theDau]=part->Pz();
	++theDau;
      }
    }
  }
  if(theDau!=2){
    AliError("Wrong number of decay prongs");
    return 99999.;
  }

  Double_t d0dummy[3]={0.,0.,0.};
  AliAODRecoDecayHF aodD0MC(vtxTrue,origD,3,charge,pXdauTrue,pYdauTrue,pZdauTrue,d0dummy);
  return aodD0MC.ImpParXY();

}
//______________________________________________________-
void AliAnalysisTaskSEDStarSpectra::CreateImpactParameterHistos(){
  // Histos for impact paramter study

  Int_t nbins[3]={400,200,fNImpParBins};
  Double_t xmin[3]={1.75,0.,fLowerImpPar};
  Double_t xmax[3]={1.98,20.,fHigherImpPar};

  fHistMassPtImpParTCDs[0]=new THnSparseF("hMassPtImpParAll",
					"Mass vs. pt vs.imppar - All",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTCDs[1]=new THnSparseF("hMassPtImpParPrompt",
					"Mass vs. pt vs.imppar - promptD",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTCDs[2]=new THnSparseF("hMassPtImpParBfeed",
					"Mass vs. pt vs.imppar - DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTCDs[3]=new THnSparseF("hMassPtImpParTrueBfeed",
					"Mass vs. pt vs.true imppar -DfromB",
					3,nbins,xmin,xmax);
  fHistMassPtImpParTCDs[4]=new THnSparseF("hMassPtImpParBkg",
				        "Mass vs. pt vs.imppar - backgr.",
					3,nbins,xmin,xmax);

  for(Int_t i=0; i<5;i++){
    fOutput->Add(fHistMassPtImpParTCDs[i]);
  }
}

