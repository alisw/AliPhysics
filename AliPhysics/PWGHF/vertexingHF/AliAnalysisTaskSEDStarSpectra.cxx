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
#include "AliNormalizationCounter.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEDStarSpectra.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDStarSpectra);
/// \endcond

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
  fAODProtection(1),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fNPtBins(0),
  fAllhist(0x0),
  fPIDhist(0x0),
  fDoDStarVsY(kFALSE)
{
  //
  /// Default ctor
  //
  for(Int_t i=0;i<5;i++) fHistMassPtImpParTCDs[i]=0;
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
  fAODProtection(1),
  fDoImpParDstar(kFALSE),
  fNImpParBins(400),
  fLowerImpPar(-2000.),
  fHigherImpPar(2000.),
  fNPtBins(0),
  fAllhist(0x0),
  fPIDhist(0x0),
    fDoDStarVsY(kFALSE)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarSpectra","Calling Constructor");

  fCuts=cuts;
  for(Int_t i=0;i<5;i++) fHistMassPtImpParTCDs[i]=0;


  DefineOutput(1,TList::Class());  //counters
  DefineOutput(2,TList::Class());  //All Entries output
  DefineOutput(3,TList::Class());  //3sigma PID output
  DefineOutput(4,AliRDHFCutsDStartoKpipi::Class());   //My private output
  DefineOutput(5,AliNormalizationCounter::Class());   // normalization
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::~AliAnalysisTaskSEDStarSpectra() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEDStarSpectra","Calling Destructor");
  
  delete fOutput;
  delete fOutputAll;
  delete fOutputPID;
  delete fCuts;
  delete fCEvents;
  delete fDeltaMassD1;
  for(Int_t i=0; i<5; i++){
    delete fHistMassPtImpParTCDs[i];
  }
  for(Int_t i=0; i<((fNPtBins+2)*18); i++){
    delete fAllhist[i];
    delete fPIDhist[i];
  }
  delete [] fAllhist;
  delete [] fPIDhist;

}
//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::Init(){
  //
  /// Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");
   AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
 fNPtBins=fCuts->GetNPtBins();
  // Post the data
  PostData(4,copyfCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserExec(Option_t *)
{
  /// user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  fCEvents->Fill(0);//all events
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fCEvents->Fill(8);
      return;
    }
    fCEvents->Fill(1);
  }

  fEvents++;

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi=0;
  TClonesArray *arrayD0toKpi=0;

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
      arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else {
    arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
    arrayD0toKpi=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
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
  fCEvents->Fill(3);
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
  fCEvents->Fill(4);

  if (!arrayDStartoD0pi || !arrayD0toKpi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast())); 

  Int_t nSelectedAna =0;
  Int_t nSelectedProd =0;
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  // loop over the tracks to search for candidates soft pion 
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {

    // D* candidates and D0 from D*
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    AliAODRecoDecayHF2Prong *trackD0;
    if(dstarD0pi->GetIsFilled()<1){
      trackD0 = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->At(dstarD0pi->GetProngID(1));
    } else {
      trackD0 = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    }

    fCEvents->Fill(10);
    TObjArray arrTracks(3);
    for(Int_t ipr=0;ipr<3;ipr++){
      AliAODTrack *tr;
      if(ipr == 0) tr=vHF->GetProng(aodEvent,dstarD0pi,ipr); //soft pion
      else         tr=vHF->GetProng(aodEvent,trackD0,ipr-1); //D0 daughters
      arrTracks.AddAt(tr,ipr);
    }
    if(!fCuts->PreSelect(arrTracks)){
      fCEvents->Fill(13);
      continue;
    }
    
    Bool_t isDStarCand =kTRUE;
    if(!(vHF->FillRecoCasc(aodEvent,dstarD0pi,isDStarCand))) {//Fill the data members of the candidate only if they are empty.
      fCEvents->Fill(12); //monitor how often this fails
      continue;
    }
    if(!dstarD0pi->GetSecondaryVtx()) continue;
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    if (!theD0particle) continue;
    
    Int_t isDStar = 0;
    TClonesArray *mcArray = 0;
    AliAODMCHeader *mcHeader=0;

    Bool_t isPrimary=kTRUE;
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
	AliAODMCParticle *dg0 = (AliAODMCParticle*)mcArray->At(partDSt->GetDaughterLabel(0));
	//	AliAODMCParticle *dg01 = (AliAODMCParticle*)mcArray->At(dg0->GetDaughterLabel(0));
	pdgCode=TMath::Abs(partDSt->GetPdgCode());
	if(!isPrimary){
	  trueImpParXY=GetTrueImpactParameterD0(mcHeader,mcArray,dg0)*1000.;
	}
	isDStar = 1;
      }else{
	pdgCode=-1;
      }
    }
   
    if(pdgCode==-1) AliDebug(2,"No particle assigned! check\n");

    Double_t Dstarpt = dstarD0pi->Pt();
    
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
   
  // set the D0 and D* search window  bin by bin - D* window useful to speed up the reconstruction and D0 window used *ONLY* to calculate side band bkg for the background subtraction methods, for the standard analysis the value in the cut file is considered
    
    if (0<=Dstarpt && Dstarpt<0.5){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (0.5<=Dstarpt && Dstarpt<1.0){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (1.0<=Dstarpt && Dstarpt<2.0){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.020;
	fPeakWindow=0.0018;
      }
    }
    if (2.0<=Dstarpt && Dstarpt<3.0){
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.022;
	fPeakWindow=0.0016;
      }
    }
    if (3.0<=Dstarpt && Dstarpt<4.0){ 
      if(fAnalysis==1){
	fD0Window=0.035;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    }
    if (4.0<=Dstarpt && Dstarpt<5.0){
      if(fAnalysis==1){
	fD0Window=0.045;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.0014;
      }
    } 
   if (5.0<=Dstarpt && Dstarpt<6.0){
      if(fAnalysis==1){
	fD0Window=0.045;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.006;
      }
    } 
    if (6.0<=Dstarpt && Dstarpt<7.0){
      if(fAnalysis==1){
	fD0Window=0.055;
	fPeakWindow=0.03;
      }else{
	fD0Window=0.026;
	fPeakWindow=0.006;
      }
    }
    if (Dstarpt>=7.0){
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
    Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate,aodEvent); //selected
    if(isSelected>0) fCEvents->Fill(11);

    // after cuts
    if(fDoImpParDstar && isSelected){
      fHistMassPtImpParTCDs[0]->Fill(arrayForSparse);
      if(isPrimary) fHistMassPtImpParTCDs[1]->Fill(arrayForSparse);
      else{
	fHistMassPtImpParTCDs[2]->Fill(arrayForSparse);
	fHistMassPtImpParTCDs[3]->Fill(arrayForSparseTrue);
      }
    }
    
    if (fDoDStarVsY && isSelected){
      ((TH3F*) (fOutputPID->FindObject("deltamassVsyVsPt")))->Fill(dstarD0pi->DeltaInvMass(),dstarD0pi->YDstar(),dstarD0pi->Pt() );  
    }

    
    // fill PID
    FillSpectrum(dstarD0pi,isDStar,fCuts,isSelected,fOutputPID, fPIDhist);
    SideBandBackground(dstarD0pi,fCuts,isSelected, fOutputPID, fPIDhist);
    //WrongSignForDStar(dstarD0pi,fCuts,fOutputPID);

    //swich off the PID selection
    fCuts->SetUsePID(kFALSE);
    Int_t isSelectedNoPID=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate, aodEvent); //selected
    fCuts->SetUsePID(kTRUE);

    FillSpectrum(dstarD0pi,isDStar,fCuts,isSelectedNoPID,fOutputAll, fAllhist);
    //    SideBandBackground(dstarD0pi,fCuts,isSelectedNoPID, fOutputAll);

    // rare D search ------ 
    if(fDoSearch){
      TLorentzVector lorentzTrack1(0,0,0,0); // lorentz 4 vector
      TLorentzVector lorentzTrack2(0,0,0,0); // lorentz 4 vector
      
      for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){ 
	
	AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
	if(!aodTrack) AliFatal("Not a standard AOD");
	
	if(dstarD0pi->Charge() == aodTrack->Charge()) continue;
	if((!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
	if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))>0.02) continue;
	
	//build the D1 mass
	Double_t mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();

	lorentzTrack1.SetPxPyPzE( dstarD0pi->Px(),dstarD0pi->Py(), dstarD0pi->Pz(), dstarD0pi->E(413) );
	lorentzTrack2.SetPxPyPzE( aodTrack->Px(),aodTrack->Py(), aodTrack->Pz(),aodTrack->E(mass) );
	
	//D1 mass
	Double_t d1mass = ((lorentzTrack1+lorentzTrack2).M());
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
  
  delete vHF;

  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));
  
  PostData(1,fOutput);
  PostData(2,fOutputAll);
  PostData(3,fOutputPID);
  PostData(5,fCounter);

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarSpectra::Terminate(Option_t*)
{    
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.
  
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
 /// output
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
  /// Create histograms

  fCEvents = new TH1F("fCEvents","counter",14,0,14);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");
  fCEvents->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fCEvents->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fCEvents->GetXaxis()->SetBinLabel(3,"good prim vtx and B field");
  fCEvents->GetXaxis()->SetBinLabel(4,"no event selected");
  fCEvents->GetXaxis()->SetBinLabel(5,"no vtx contributors");
  fCEvents->GetXaxis()->SetBinLabel(6,"trigger for PbPb");
  fCEvents->GetXaxis()->SetBinLabel(7,"no z vtx");
  fCEvents->GetXaxis()->SetBinLabel(9,"nEvents Mismatched dAOD");
  fCEvents->GetXaxis()->SetBinLabel(11, "no. of cascade candidates");
  fCEvents->GetXaxis()->SetBinLabel(12, "no. of Dstar after selection cuts");
  fCEvents->GetXaxis()->SetBinLabel(13, "no. of not on-the-fly rec Dstar");
  fCEvents->GetXaxis()->SetBinLabel(14, "no. of Dstar rejected by preselect"); //toadd

  fOutput->Add(fCEvents);

  fTrueDiff2 = new TH2F("DiffDstar_pt","True Reco diff vs pt",200,0,15,900,0,0.3);
  fOutput->Add(fTrueDiff2);

  fDeltaMassD1 = new TH1F("DeltaMassD1","delta mass d1",600,0,0.8);
  fOutput->Add(fDeltaMassD1);
  //temp a

  fAllhist = new TH1F*[(fNPtBins+2)*18];
  fPIDhist = new TH1F*[(fNPtBins+2)*18];

  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  for(Int_t i=-2;i<fNPtBins;i++){
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
    
    TH1F* spectrumMass = new TH1F(nameMass.Data(),"D^{*}-D^{0} invariant mass; #DeltaM [GeV/c^{2}]; Entries",700,0.13,0.2);
    TH1F* spectrumSgn = new TH1F(nameSgn.Data(), "D^{*}-D^{0} Signal invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",700,0.13,0.2);
    TH1F* spectrumBkg = new TH1F(nameBkg.Data(), "D^{*}-D^{0} Background invariant mass - MC; #DeltaM [GeV/c^{2}]; Entries",700,0.13,0.2);
    
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
    fAllhist[i+2+((fNPtBins+2)*kDeltaMass)]=allMass;
    fAllhist[i+2+((fNPtBins+2)*kDeltaSgn)]=allSgn;
    fAllhist[i+2+((fNPtBins+2)*kDeltaBkg)]=allBkg;

    fOutputPID->Add(pidMass);
    fOutputPID->Add(pidSgn);
    fOutputPID->Add(pidBkg);
    fPIDhist[i+2+((fNPtBins+2)*kDeltaMass)]=pidMass;
    fPIDhist[i+2+((fNPtBins+2)*kDeltaSgn)]=pidSgn;
    fPIDhist[i+2+((fNPtBins+2)*kDeltaBkg)]=pidBkg;

    TH1F* allD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* allD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* allD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    TH1F* pidD0Mass = (TH1F*)spectrumD0Mass->Clone();
    TH1F* pidD0Sgn  = (TH1F*)spectrumD0Sgn->Clone();
    TH1F* pidD0Bkg  = (TH1F*)spectrumD0Bkg->Clone();

    fOutputAll->Add(allD0Mass);
    fOutputAll->Add(allD0Sgn);
    fOutputAll->Add(allD0Bkg);
    fAllhist[i+2+((fNPtBins+2)*kDzMass)]=allD0Mass;
    fAllhist[i+2+((fNPtBins+2)*kDzSgn)]=allD0Sgn;
    fAllhist[i+2+((fNPtBins+2)*kDzBkg)]=allD0Bkg;

    fOutputPID->Add(pidD0Mass);
    fOutputPID->Add(pidD0Sgn);
    fOutputPID->Add(pidD0Bkg);
    fPIDhist[i+2+((fNPtBins+2)*kDzMass)]=pidD0Mass;
    fPIDhist[i+2+((fNPtBins+2)*kDzSgn)]=pidD0Sgn;
    fPIDhist[i+2+((fNPtBins+2)*kDzBkg)]=pidD0Bkg;
  
    TH1F* allDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* allDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* allDstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    TH1F* pidDstarMass = (TH1F*)spectrumDstarMass->Clone();
    TH1F* pidDstarSgn = (TH1F*)spectrumDstarSgn->Clone();
    TH1F* pidDstarBkg = (TH1F*)spectrumDstarBkg->Clone();
    
    fOutputAll->Add(allDstarMass);
    fOutputAll->Add(allDstarSgn);
    fOutputAll->Add(allDstarBkg);
    fAllhist[i+2+((fNPtBins+2)*kDstarMass)]=allDstarMass;
    fAllhist[i+2+((fNPtBins+2)*kDstarSgn)]=allDstarSgn;
    fAllhist[i+2+((fNPtBins+2)*kDstarBkg)]=allDstarBkg;

    fOutputPID->Add(pidDstarMass);
    fOutputPID->Add(pidDstarSgn);
    fOutputPID->Add(pidDstarBkg);
    fPIDhist[i+2+((fNPtBins+2)*kDstarMass)]=pidDstarMass;
    fPIDhist[i+2+((fNPtBins+2)*kDstarSgn)]=pidDstarSgn;
    fPIDhist[i+2+((fNPtBins+2)*kDstarBkg)]=pidDstarBkg;

    TH1F* allSideBandMass = (TH1F*)spectrumSideBandMass->Clone();
    TH1F* pidSideBandMass = (TH1F*)spectrumSideBandMass->Clone();

    fOutputAll->Add(allSideBandMass);
    fOutputPID->Add(pidSideBandMass);
    fAllhist[i+2+((fNPtBins+2)*kSideBandMass)]=allSideBandMass;
    fPIDhist[i+2+((fNPtBins+2)*kSideBandMass)]=pidSideBandMass;
   
    TH1F* allWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();
    TH1F* pidWrongSignMass = (TH1F*)spectrumWrongSignMass->Clone();

    fOutputAll->Add(allWrongSignMass);
    fOutputPID->Add(pidWrongSignMass);
    fAllhist[i+2+((fNPtBins+2)*kWrongSignMass)]=allWrongSignMass;
    fPIDhist[i+2+((fNPtBins+2)*kWrongSignMass)]=pidWrongSignMass;
   
  }

  // pt spectra
  nameMass="ptMass";
  nameSgn="ptSgn";
  nameBkg="ptBkg";
  
  TH1F* ptspectrumMass = new TH1F(nameMass.Data(),"D^{*} p_{T}; p_{T} [GeV]; Entries",400,0,50);
  TH1F* ptspectrumSgn = new TH1F(nameSgn.Data(), "D^{*} Signal p_{T} - MC; p_{T} [GeV]; Entries",400,0,50);
  TH1F* ptspectrumBkg = new TH1F(nameBkg.Data(), "D^{*} Background p_{T} - MC; p_{T} [GeV]; Entries",400,0,50);
  
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
  fAllhist[((fNPtBins+2)*kptMass)]=ptallMass;
  fAllhist[((fNPtBins+2)*kptSgn)]=ptallSgn;
  fAllhist[((fNPtBins+2)*kptBkg)]=ptallBkg;

  fOutputPID->Add(ptpidMass);
  fOutputPID->Add(ptpidSgn);
  fOutputPID->Add(ptpidBkg);
  fPIDhist[(fNPtBins+2)*kptMass]=ptpidMass;
  fPIDhist[(fNPtBins+2)*kptSgn]=ptpidSgn;
  fPIDhist[(fNPtBins+2)*kptBkg]=ptpidBkg;
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
  fAllhist[(fNPtBins+2)*ketaMass]=etaallMass;
  fAllhist[(fNPtBins+2)*ketaSgn]=etaallSgn;
  fAllhist[(fNPtBins+2)*ketaBkg]=etaallBkg;
  
  fOutputPID->Add(etapidMass);
  fOutputPID->Add(etapidSgn);
  fOutputPID->Add(etapidBkg);
  fPIDhist[(fNPtBins+2)*ketaMass]=etapidMass;
  fPIDhist[(fNPtBins+2)*ketaSgn]=etapidSgn;
  fPIDhist[(fNPtBins+2)*ketaBkg]=etapidBkg;

  if (fDoDStarVsY){  
    TH3F* deltamassVsyVsPtPID = new TH3F("deltamassVsyVsPt", "delta mass Vs y Vs pT;  #DeltaM [GeV/c^{2}]; y; p_{T} [GeV/c]", 700,0.13,0.2, 40, -1, 1, 36, 0., 36.);
    fOutputPID->Add(deltamassVsyVsPtPID);
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, AliRDHFCutsDStartoKpipi *cuts,Int_t isSel, TList *listout, TH1F** histlist){
  //
  /// Fill histos for D* spectrum
  //
  
  if(!isSel) return;

  // D0 window
  Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();  
 

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
      histlist[ptbin+1+((fNPtBins+2)*kDzSgn)]->Fill(invmassD0);
      histlist[(fNPtBins+2)*kDzSgn]->Fill(invmassD0);
      histlist[ptbin+1+((fNPtBins+2)*kDstarSgn)]->Fill(invmassDstar);
      histlist[(fNPtBins+2)*kDstarSgn]->Fill(invmassDstar);
      histlist[ptbin+1+((fNPtBins+2)*kDeltaSgn)]->Fill(invmassDelta);
      histlist[(fNPtBins+2)*kDeltaSgn]->Fill(invmassDelta);
    if (massInRange) {
	histlist[(fNPtBins+2)*kptSgn]->Fill(pt);
        histlist[(fNPtBins+2)*ketaSgn]->Fill(eta);
      }
    }
    else {//background
      histlist[ptbin+1+((fNPtBins+2)*kDzBkg)]->Fill(invmassD0);
      histlist[(fNPtBins+2)*kDzBkg]->Fill(invmassD0);
      histlist[ptbin+1+((fNPtBins+2)*kDstarBkg)]->Fill(invmassDstar);
      histlist[(fNPtBins+2)*kDstarBkg]->Fill(invmassDstar);
      histlist[ptbin+1+((fNPtBins+2)*kDeltaBkg)]->Fill(invmassDelta);
      histlist[(fNPtBins+2)*kDeltaBkg]->Fill(invmassDelta);
     if (massInRange) {
        histlist[(fNPtBins+2)*kptBkg]->Fill(pt);
        histlist[(fNPtBins+2)*ketaBkg]->Fill(eta);
      }
    }
  }
  //no MC info, just cut selection
  histlist[ptbin+1+((fNPtBins+2)*kDzMass)]->Fill(invmassD0);
  histlist[(fNPtBins+2)*kDzMass]->Fill(invmassD0);
  histlist[ptbin+1+((fNPtBins+2)*kDstarMass)]->Fill(invmassDstar);
  histlist[(fNPtBins+2)*kDstarMass]->Fill(invmassDstar);
  histlist[ptbin+1+((fNPtBins+2)*kDeltaMass)]->Fill(invmassDelta);
  histlist[(fNPtBins+2)*kDeltaMass]->Fill(invmassDelta);
  
  if (massInRange) {
    histlist[(fNPtBins+2)*kptMass]->Fill(pt);
    histlist[(fNPtBins+2)*ketaMass]->Fill(eta);
  }
 
  return;
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarSpectra::SideBandBackground(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist){

  ///  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas
  /// (expected detector resolution) on the left and right frm the D0 mass. Each band
  ///  has a width of ~5 sigmas. Two band needed  for opening angle considerations

  if(!isSel) return;

  Int_t ptbin=cuts->PtBin(part->Pt());
    
  // select the side bands intervall
  Double_t invmassD0    = part->InvMassD0();
  if(TMath::Abs(invmassD0-1.865)>4*fD0Window && TMath::Abs(invmassD0-1.865)<8*fD0Window){
    
    // for pt and eta
    Double_t invmassDelta = part->DeltaInvMass();
    
    histlist[ptbin+1+((fNPtBins+2)*kSideBandMass)]->Fill(invmassDelta);
    histlist[(fNPtBins+2)*kSideBandMass]->Fill(invmassDelta);
    
    
  }
}
//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::WrongSignForDStar(AliAODRecoCascadeHF *part,  AliRDHFCutsDStartoKpipi *cuts, TList *listout){
  //
  /// assign the wrong charge to the soft pion to create background
  //
  Int_t ptbin=cuts->PtBin(part->Pt());
  
  Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t invmassD0   = part->InvMassD0();  
  if (TMath::Abs(invmassD0-mPDGD0)>fD0Window) return; 

  AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)part->Get2Prong();

  Int_t okDzWrongSign;
  Double_t wrongMassD0=0.;
  
  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
   if (!isSelected){
    return;
  }

  okDzWrongSign =  1;
  
  //if is D*+ than assume D0bar
  if(part->Charge()>0 && (isSelected ==1)) { 
    okDzWrongSign = 0;
  }
  
  // assign the wrong mass in case the cuts return both D0 and D0bar
  if(part->Charge()>0 && (isSelected ==3)) { 
    okDzWrongSign = 0;
  }
  
  //wrong D0 inv mass
  if(okDzWrongSign!=0){
    wrongMassD0 = theD0particle->InvMassD0();
  }else if(okDzWrongSign==0){
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
Int_t AliAnalysisTaskSEDStarSpectra::CheckOrigin(TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
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
Float_t AliAnalysisTaskSEDStarSpectra::GetTrueImpactParameterD0(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const {
  /// true impact parameter calculation

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partDp->XvYvZv(origD);	  
  Short_t charge=partDp->Charge();
  Double_t pXdauTrue[3],pYdauTrue[3],pZdauTrue[3];
  Int_t labelFirstDau = partDp->GetDaughterLabel(0); 

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
  /// Histos for impact paramter study

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

