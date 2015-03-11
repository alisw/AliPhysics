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
//
//
// Task for selecting D mesons to be used as an input for D-Jet correlations
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL) XMZhang@lbl.gov
// S. Aiola (Yale University) salvatore.aiola@cern.ch
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TH3F.h>

#include "AliLog.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSEDmesonsFilterCJ.h"

ClassImp(AliAnalysisTaskSEDmesonsFilterCJ)

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ() :
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fUseReco(kTRUE),
  fCandidateType(0),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fCuts(0),
  fMinMass(0.),
  fMaxMass(0.),
  fInhibitTask(kFALSE),
  fInitOk(kFALSE),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fMCarray(0),
  fCandidateArray(0),
  fSideBandArray(0),
  fOutput(0),
  fHistStat(0),
  fHistNSBCandEv(0),
  fHistNCandEv(0),
  fHistImpParS(0),
  fHistImpParB(0),
  fHistPtPion(0),
  fHistInvMassPtD(0),
  fHistInvMassS(0),
  fHistInvMassB(0),
  fHistAlphaDDS(0),
  fHistAlphaDpisS(0),
  fHistAlphaDpiS(0),
  fHistAlphaDKS(0),
  fHistAlphaDDB(0),
  fHistAlphaDpisB(0),
  fHistAlphaDpiB(0),
  fHistAlphaDKB(0),
  fHistDeltaRDDS(0),
  fHistDeltaRDpisS(0),
  fHistDeltaRDpiS(0),
  fHistDeltaRDKS(0),
  fHistDeltaRDDB(0),
  fHistDeltaRDpisB(0),
  fHistDeltaRDpiB(0),
  fHistDeltaRDKB(0),
  fHistAlphaDpiR(0),
  fHistAlphaDKR(0),
  fHistDeltaRDpiR(0),
  fHistDeltaRDKR(0)
{
   //
   // Default constructor
   //
   
   for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
   for (Int_t i=30; i--;) fSigmaD0[i] = 0.;
}

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::AliAnalysisTaskSEDmesonsFilterCJ(const char *name, AliRDHFCuts *cuts, ECandidateType candtype) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fUseReco(kTRUE),
  fCandidateType(candtype),
  fCandidateName(""),
  fPDGmother(0),
  fNProngs(0),
  fBranchName(""),
  fCuts(cuts),
  fMinMass(0.),
  fMaxMass(0.),
  fInhibitTask(kFALSE),
  fInitOk(kFALSE),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fMCarray(0),
  fCandidateArray(0),
  fSideBandArray(0),
  fOutput(0),
  fHistStat(0),
  fHistNSBCandEv(0),
  fHistNCandEv(0),
  fHistImpParS(0),
  fHistImpParB(0),
  fHistPtPion(0),
  fHistInvMassPtD(0),
  fHistInvMassS(0),
  fHistInvMassB(0),
  fHistAlphaDDS(0),
  fHistAlphaDpisS(0),
  fHistAlphaDpiS(0),
  fHistAlphaDKS(0),
  fHistAlphaDDB(0),
  fHistAlphaDpisB(0),
  fHistAlphaDpiB(0),
  fHistAlphaDKB(0),
  fHistDeltaRDDS(0),
  fHistDeltaRDpisS(0),
  fHistDeltaRDpiS(0),
  fHistDeltaRDKS(0),
  fHistDeltaRDDB(0),
  fHistDeltaRDpisB(0),
  fHistDeltaRDpiB(0),
  fHistDeltaRDKB(0),
  fHistAlphaDpiR(0),
  fHistAlphaDKR(0),
  fHistDeltaRDpiR(0),
  fHistDeltaRDKR(0)
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //
   
   Info("AliAnalysisTaskSEDmesonsFilterCJ","Calling Constructor");
   
   for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
   for (Int_t i=30; i--;) fSigmaD0[i] = 0.;
   
   const Int_t nptbins = fCuts->GetNPtBins();
   Float_t defaultSigmaD013[13] = { 0.012, 0.012, 0.012, 0.015, 0.015, 0.018, 0.018, 0.020, 0.020, 0.030, 0.030, 0.037, 0.040 };
   
   switch (fCandidateType) {
   case 0 :
      fCandidateName = "D0";
      fPDGmother = 421;
      fNProngs = 2;
      fPDGdaughters[0] = 211;  // pi 
      fPDGdaughters[1] = 321;  // K
      fPDGdaughters[2] = 0;    // empty
      fPDGdaughters[3] = 0;    // empty
      fBranchName = "D0toKpi";
      break;
   case 1 :
      fCandidateName = "Dstar";
      fPDGmother = 413;
      fNProngs = 3;
      fPDGdaughters[1] = 211; // pi soft
      fPDGdaughters[0] = 421; // D0
      fPDGdaughters[2] = 211; // pi fromD0
      fPDGdaughters[3] = 321; // K from D0
      fBranchName = "Dstar";
      
      if (nptbins<=13) {
      	 for (Int_t ipt=0; ipt<nptbins;ipt++) fSigmaD0[ipt] = defaultSigmaD013[ipt];
      }
      else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
   default :
     Warning("AliAnalysisTaskSEDmesonsFilterCJ", "%d not accepted!!", fCandidateType);
     break;
   }
   
   if (fCandidateType == kD0toKpi) SetMassLimits(0.15, fPDGmother);
   if (fCandidateType == kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
   
   DefineOutput(1, TList::Class());       // histos
   DefineOutput(2, AliRDHFCuts::Class()); // my cuts
   DefineOutput(3, TClonesArray::Class()); //array of candidates
   DefineOutput(4, TClonesArray::Class()); //array of SB candidates
}

//_______________________________________________________________________________
AliAnalysisTaskSEDmesonsFilterCJ::~AliAnalysisTaskSEDmesonsFilterCJ()
{
   //
   // Destructor
   //
   
   Info("~AliAnalysisTaskSEDmesonsFilterCJ","Calling Destructor");  
   
   if (fOutput) { delete fOutput; fOutput = 0; }
   if (fCuts)   { delete fCuts;   fCuts   = 0; }
   if (fCandidateArray)  { delete fCandidateArray;  fCandidateArray = 0; }
   delete fSideBandArray;
   
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Init()
{
   //
   // Initialization
   //
   
  Info("AnalysisTaskSEDmesonsForJetCorrelations::Init()", "Entering method");
   
   switch (fCandidateType) {
   case 0: 
      {
   	 AliRDHFCutsD0toKpi* copyfCutsDzero = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
   	 copyfCutsDzero->SetName("AnalysisCutsDzero");
   	 PostData(2, copyfCutsDzero);  // Post the data
      } break;
   case 1: 
      {
      	 AliRDHFCutsDStartoKpipi* copyfCutsDstar = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      	 copyfCutsDstar->SetName("AnalysisCutsDStar");
      	 PostData(2, copyfCutsDstar); // Post the cuts
      } break;
   default:
     return;
   }
   
   return;
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects()
{ 
   //
   // Create output objects 
   //
   
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s", GetName());
   
   fOutput = new TList();
   fOutput->SetOwner();
   DefineHistoForAnalysis(); // define histograms
   
   if (fCandidateType == kD0toKpi){
      fCandidateArray = new TClonesArray("AliAODRecoDecayHF",0);
      fSideBandArray = new TClonesArray("AliAODRecoDecayHF",0); 
   }
   else if (fCandidateType == kDstartoKpipi) {
      fCandidateArray = new TClonesArray("AliAODRecoCascadeHF",0);
      fSideBandArray = new TClonesArray("AliAODRecoCascadeHF",0); 
   }
   else {
     AliWarning(Form("Candidate type %d not recognized!", fCandidateType));
     return;
   }
   
   fCandidateArray->SetOwner();
   fCandidateArray->SetName(Form("fCandidateArray%s%s",fCandidateName.Data(),fUseReco ? "rec" : "gen"));
   
   //this is used for the DStar side bands and MC!
   fSideBandArray->SetOwner();
   fSideBandArray->SetName(Form("fSideBandArray%s%s",fCandidateName.Data(),fUseReco ? "rec" : "gen"));
  
   PostData(1, fOutput);
   PostData(3, fCandidateArray);
   PostData(4, fSideBandArray);
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ExecOnce()
{
  //
  // To be executed only once, for the first event
  //

  if (fInhibitTask) return;

  AliDebug(2, "Entering ExecOnce()");

  // Load the event
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  if (fAodEvent) {
    fArrayDStartoD0pi = (TClonesArray*)fAodEvent->GetList()->FindObject(fBranchName.Data());
  }
  else {
    if (AODEvent() && IsStandardAOD()) {
      
      // In case there is an AOD handler writing a standard AOD, use the AOD 
      // event in memory rather than the input (ESD) event.    
      fAodEvent = dynamic_cast<AliAODEvent*>(AODEvent());
      
      // in this case the branches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
      if(aodHandler->GetExtensions()) {
        AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
        AliAODEvent *aodFromExt = ext->GetAOD();
        fArrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      }
      else {
        AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
        fInhibitTask = kTRUE;
        return;
      }
    }
  }
   
  if (!fArrayDStartoD0pi) {
    AliError(Form("Could not find array %s, skipping the event. Task '%s' will be disabled!", fBranchName.Data(), GetName()));
    fInhibitTask = kTRUE;
    return;
  }

  if (fUseMCInfo) {
    fMCarray = dynamic_cast<TClonesArray*>(fAodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMCarray) {
      AliError(Form("MC particles not found! Task '%s' will be disabled!", GetName()));
      fInhibitTask = kTRUE;
      return;
    }
  }

  fInitOk = kTRUE;
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ProcessD0(AliAODRecoDecayHF* charmCand, Int_t isSelected, Int_t mcLabel)
{
  //
  // Process the d0 candidate
  //

  Double_t masses[2] = {0.};

  //D0 from D0 bar
  Int_t pdgdaughtersD0[2] = { 211, 321 };     // pi,K 
  Int_t pdgdaughtersD0bar[2] = { 321, 211 };  // K,pi 

  //needed quantities
  masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
  masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
  fHistStat->Fill(3);

  Double_t ptD = charmCand->Pt();
      	 
  // mass vs pt
  if (isSelected == 1 || isSelected == 3) fHistInvMassPtD->Fill(masses[0], ptD);
  if (isSelected >= 2) fHistInvMassPtD->Fill(masses[1], ptD);
      	         	
  if (fUseMCInfo) {  //fill histograms of kinematics, using MC truth
    Double_t aD = charmCand->Phi();
    Double_t adaugh[2] = {charmCand->PhiProng(0), charmCand->PhiProng(1)};
    AliAODTrack* p0 = (AliAODTrack*)charmCand->GetDaughter(0); 
    AliAODTrack* p1 = (AliAODTrack*)charmCand->GetDaughter(1);
    Float_t dR0 = DeltaR(charmCand, p0);
    Float_t dR1 = DeltaR(charmCand, p1);

    Bool_t isMCBkg = (mcLabel <= 0);
    
    if (isMCBkg) { //background
      if (isSelected == 1 || isSelected == 3) { // selected as D0
        fHistAlphaDKB->Fill(aD-adaugh[0],ptD);
        fHistAlphaDpiB->Fill(aD-adaugh[1],ptD);

        fHistDeltaRDKB->Fill(dR0,ptD);
        fHistDeltaRDpiB->Fill(dR1,ptD);   
      }
      if (isSelected >= 2) { //selected as D0bar
        fHistAlphaDpiB->Fill(aD-adaugh[0],ptD);
        fHistAlphaDKB->Fill(aD-adaugh[1],ptD);
            
        fHistDeltaRDpiB->Fill(dR0,ptD);
        fHistDeltaRDKB->Fill(dR1,ptD);
      }
    }
    else { //signal and reflections
      AliAODMCParticle* charmPart = (AliAODMCParticle*)fMCarray->At(TMath::Abs(mcLabel));
      if (!charmPart) {
        AliError("Could not find MC particle associated with charm candidate!");
        return;
      }
      
      Int_t pdgCode = charmPart->GetPdgCode();
        
      Bool_t isD0 = kFALSE;
      if (pdgCode ==  421)  {
        isD0 = kTRUE;
      }
      else if (pdgCode == -421)  {
        isD0 = kFALSE;
      }
      else {
        AliDebug(2, "Not a D0/D0bar meson!");
        return;
      }
      
      if (isD0) { //D0
        fHistAlphaDKS->Fill(aD-adaugh[0],ptD);
        fHistAlphaDpiS->Fill(aD-adaugh[1],ptD);

        fHistDeltaRDKS->Fill(dR0,ptD);
        fHistDeltaRDpiS->Fill(dR1,ptD);
            
        if (isSelected == 3) { // selected as both D0bar/D0
          fHistAlphaDpiR->Fill(aD-adaugh[0],ptD);
          fHistAlphaDKR->Fill(aD-adaugh[1],ptD);

          fHistDeltaRDpiR->Fill(dR0,ptD);
          fHistDeltaRDKR->Fill(dR1,ptD);
        }
      }
      else { //D0bar
        fHistAlphaDKS->Fill(aD-adaugh[1],ptD);
        fHistAlphaDpiS->Fill(aD-adaugh[0],ptD);

        fHistDeltaRDKS->Fill(dR1,ptD);
        fHistDeltaRDpiS->Fill(dR0,ptD);

        if (isSelected == 3) { // selected as both D0bar/D0
          fHistAlphaDpiR->Fill(aD-adaugh[1],ptD);
          fHistAlphaDKR->Fill(aD-adaugh[0],ptD);

          fHistDeltaRDpiR->Fill(dR1, ptD);
          fHistDeltaRDKR->Fill(dR0, ptD);
        }
      }
    } //end signal and reflections 
  }// end MC
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::ProcessDstar(AliAODRecoCascadeHF* dstar, Int_t mcLabel)
{
  //
  // Process the d star candidate
  //

  AliDebug(2, "Entering method");
  
  //softpion from D* decay
  AliAODTrack *track2 = (AliAODTrack*)dstar->GetBachelor();

  Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

  Double_t ptD = dstar->Pt();
      
  // select D* in the D0 window.
  // In the cut object window is loose to allow for side bands
      
  // retrieve the corresponding pt bin for the candidate
  // and set the expected D0 width (x3)

  Int_t bin = fCuts->PtBin(ptD);
  if (bin < 0 || bin >= fCuts->GetNPtBins()) {
    AliError(Form("Pt %.3f out of bounds", ptD));
    return;
  }
      	 
  AliDebug(1, Form("Pt bin %d and sigma D0 %.4f", bin, fSigmaD0[bin]));
  //consider the Dstar candidates only if the mass of the D0 is in 3 sigma wrt the PDG value
  if ((dstar->InvMassD0()>=(mPDGD0-3.*fSigmaD0[bin])) && (dstar->InvMassD0()<=(mPDGD0+3.*fSigmaD0[bin]))) {
    AliDebug(2, "D0 mass within 3 sigma");
      	    
    // D* delta mass
    AliDebug(2, "Filling invariant mass vs pt histogram");
    fHistInvMassPtD->Fill(dstar->DeltaInvMass(), ptD); // 2 D slice for pt bins
      	    
    // D* pt and soft pion pt for good candidates  	      	
    Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t invmassDelta = dstar->DeltaInvMass();
    if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<0.0021) {
      AliDebug(2, "Filling pion pt histogram");
      fHistPtPion->Fill(track2->Pt());
    }
  }

  if (fUseMCInfo) { //fill histograms of kinematics, using MC truth
    Double_t aD  = dstar->Phi();
    Double_t apis= track2->Phi();
      	             
    AliAODRecoDecayHF2Prong* D0fromDstar = dstar->Get2Prong();
    Double_t aD0    = D0fromDstar->Phi();
    Int_t    isD0   = D0fromDstar->Charge()>0 ? kTRUE : kFALSE;
    Double_t aK     = isD0 ? D0fromDstar->PhiProng(0) : D0fromDstar->PhiProng(1);
    Double_t api    = isD0 ? D0fromDstar->PhiProng(1) : D0fromDstar->PhiProng(0);
    Double_t dRDD0  = DeltaR(dstar,D0fromDstar);
    Double_t dRDpis = DeltaR(dstar,track2);
    Double_t dRDpi  = DeltaR(dstar, isD0 ? (AliVParticle*)D0fromDstar->GetDaughter(1) : (AliVParticle*)D0fromDstar->GetDaughter(0));
    Double_t dRDK   = DeltaR(dstar, isD0 ? (AliVParticle*)D0fromDstar->GetDaughter(0) : (AliVParticle*)D0fromDstar->GetDaughter(1));

    Bool_t isMCBkg = (mcLabel <= 0);
      	    
    if (isMCBkg) {
      fHistAlphaDDB  ->Fill(aD-aD0, ptD);
      fHistAlphaDpisB->Fill(aD-apis, ptD);
      fHistAlphaDpiB ->Fill(aD-api, ptD);
      fHistAlphaDKB  ->Fill(aD-aK, ptD);
      	    
      fHistDeltaRDDB  ->Fill(dRDD0, ptD);
      fHistDeltaRDpisB->Fill(dRDpis, ptD);
      fHistDeltaRDpiB ->Fill(dRDpi, ptD);
      fHistDeltaRDKB  ->Fill(dRDK, ptD);
    
      fHistImpParB->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
      fHistImpParB->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
      fHistImpParB->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));
    }
    else {
      fHistAlphaDDS  ->Fill(aD-aD0, ptD);
      fHistAlphaDpisS->Fill(aD-apis, ptD);
      fHistAlphaDpiS ->Fill(aD-api, ptD);
      fHistAlphaDKS  ->Fill(aD-aK, ptD);
      	    
      fHistDeltaRDDS  ->Fill(dRDD0, ptD);
      fHistDeltaRDpisS->Fill(dRDpis, ptD);
      fHistDeltaRDpiS ->Fill(dRDpi, ptD);
      fHistDeltaRDKS  ->Fill(dRDK, ptD);
      
      fHistImpParS->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
      fHistImpParS->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
      fHistImpParS->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));   
    }
  }
  else {
    fHistImpParS->Fill(dstar->Getd0Prong(0), dstar->PtProng(0)); //bachelor
    AliAODRecoDecayHF2Prong* D0fromDstar = dstar->Get2Prong();
    fHistImpParS->Fill(D0fromDstar->Getd0Prong(0), D0fromDstar->PtProng(0));
    fHistImpParS->Fill(D0fromDstar->Getd0Prong(1), D0fromDstar->PtProng(1));
  }

  AliDebug(2, "Exiting method");
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::UserExec(Option_t *)
{
  //
  // Analysis execution
  //
 
  if (!fInitOk) ExecOnce();
  if (fInhibitTask) return;

  AliDebug(2, "Entering UserExec()");
  
  fHistStat->Fill(0);
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return;
   
  //Event selection
  Bool_t iseventselected = fCuts->IsEventSelected(fAodEvent);
  if (!iseventselected) return;
  fHistStat->Fill(1);

  AliDebug(2, "Event selected");
  
  const Int_t nD = fArrayDStartoD0pi->GetEntriesFast();
  AliDebug(2, Form("Found %d vertices", nD));
  if (!fUseMCInfo) fHistStat->Fill(2, nD);

  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
  Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi
   
  //D0 from D0 bar
  Int_t pdgdaughtersD0[2] = { 211, 321 };     // pi,K 
  Int_t pdgdaughtersD0bar[2] = { 321, 211 };  // K,pi 

  const Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
  Int_t pdgMeson = 413;
  if (fCandidateType == kD0toKpi) pdgMeson = 421;

  //clear the TClonesArray from the previous event
  fCandidateArray->Clear();
  fSideBandArray->Clear();
  AliDebug(2, "TClonesArray cleared");

  Int_t iCand = 0;
  Int_t iSBCand = 0;
  
  for (Int_t icharm=0; icharm<nD; icharm++) {   //loop over D candidates
    Bool_t isMCBkg = kFALSE;
    Int_t mcLabel = -9999;
    Int_t isSelected = 0;

    AliAODRecoCascadeHF* dstar = 0;
    AliAODMCParticle* charmPart = 0;
    
    AliAODRecoDecayHF* charmCand = (AliAODRecoDecayHF*)fArrayDStartoD0pi->At(icharm); // D candidates
    if (!charmCand) continue;

    Double_t ptD = charmCand->Pt();

    // region of interest + cuts
    if (!fCuts->IsInFiducialAcceptance(ptD, charmCand->Y(pdgMeson))) continue;
    
    if (fCandidateType == kDstartoKpipi) dstar = (AliAODRecoCascadeHF*)charmCand;
      
    if (fUseMCInfo) { // Look in MC, try to simulate the z
      if (fCandidateType == kDstartoKpipi) {
        mcLabel = dstar->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, fMCarray);
      }
      else if (fCandidateType == kD0toKpi) {
        mcLabel = charmCand->MatchToMC(421, fMCarray, fNProngs, fPDGdaughters);
      }
      
      if (mcLabel <= 0) {
        AliDebug(2, Form("Charm candidate %d is background", icharm));
        isMCBkg = kTRUE;
        fHistStat->Fill(5);
      }
      else {
        AliDebug(2, Form("Charm candidate %d is *NOT* background", icharm));
        isMCBkg = kFALSE;  // this is redundant, but better keep it!
        charmPart = (AliAODMCParticle*)fMCarray->At(mcLabel);
        fHistStat->Fill(2);  
      }
    } // if fUseMCInfo
      
    if (!fUseMCInfo && fCandidateType == kDstartoKpipi) {
      //select by track cuts the side band candidates (don't want mass cut)
      isSelected = fCuts->IsSelected(charmCand, AliRDHFCuts::kTracks, fAodEvent); 
      if (!isSelected) continue;
      //add a reasonable cut on the invariant mass (e.g. (+-2\sigma, +-10 \sigma), with \sigma = fSigmaD0[bin])
      Int_t bin = fCuts->PtBin(ptD);
      if (bin < 0 || bin >= fCuts->GetNPtBins()) {
        AliError(Form("Pt %.3f out of bounds", ptD));
        continue;
      }
      
      //if data and Dstar from D0 side band
      if (((dstar->InvMassD0()<=(mPDGD0-3.*fSigmaD0[bin])) && (dstar->InvMassD0()>(mPDGD0-10.*fSigmaD0[bin]))) /*left side band*/   ||
          ((dstar->InvMassD0()>=(mPDGD0+3.*fSigmaD0[bin])) && (dstar->InvMassD0()<(mPDGD0+10.*fSigmaD0[bin]))) /*right side band*/) {	
      	    
        new ((*fSideBandArray)[iSBCand]) AliAODRecoCascadeHF(*dstar);
        iSBCand++;
      }
    }
    
    //candidate selected by cuts and PID
    isSelected = fCuts->IsSelected(charmCand, AliRDHFCuts::kAll, fAodEvent); //selected
    if (!isSelected) continue;
      
    Int_t nprongs = charmCand->GetNProngs();
    AliDebug(2, Form("Candidate is %d, and nprongs = %d", fCandidateType, nprongs));
    //for MC background fill fSideBandArray (which is instead filled above for DStar in case of data for the side bands candidates)
    if (fUseMCInfo && isMCBkg) { 
      if (fUseReco) {
        if (fCandidateType == kDstartoKpipi) {
          new ((*fSideBandArray)[iSBCand]) AliAODRecoCascadeHF(*dstar);
          fHistInvMassB->Fill(dstar->DeltaInvMass());
          AliDebug(2, Form("Now filling background inv mass histogram for charm candidate %d", icharm));
        }
        else if (fCandidateType == kD0toKpi) {
          new ((*fSideBandArray)[iSBCand]) AliAODRecoDecayHF(*charmCand);
          for (Int_t kd = 0; kd < nprongs; kd++) {
            fHistImpParB->Fill(charmCand->Getd0Prong(kd), charmCand->PtProng(kd));
            AliDebug(2, Form("Prong %d: Getd0Prong = %.3f, PtProng = %.3f", kd, charmCand->Getd0Prong(kd), charmCand->PtProng(kd)));
          }
          Double_t masses[2];
          masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
          masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
          fHistInvMassB->Fill(masses[0]);
          fHistInvMassB->Fill(masses[1]);
        }
        iSBCand++;
      }
    }
    //for data and MC signal fill fCandidateArray
    else {
      // for data or MC with the requirement fUseReco fill with candidates	       
      if (fUseReco) {
        if (fCandidateType == kDstartoKpipi) {  // D*
          new ((*fCandidateArray)[iCand]) AliAODRecoCascadeHF(*dstar);
          fHistInvMassS->Fill(dstar->DeltaInvMass());
          AliDebug(2, Form("Now filling signal inv mass histogram for charm candidate %d (mcLabel = %d, isMCBkg = %d", icharm, mcLabel, isMCBkg));
        }
        else {  //D0
          new ((*fCandidateArray)[iCand]) AliAODRecoDecayHF(*charmCand);
          for (Int_t kd = 0; kd < nprongs; kd++) {
            fHistImpParS->Fill(charmCand->Getd0Prong(kd), charmCand->PtProng(kd));
            AliDebug(2, Form("Prong %d: Getd0Prong = %.3f, PtProng = %.3f", kd, charmCand->Getd0Prong(kd), charmCand->PtProng(kd)));
          }
          Double_t masses[2];
          masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
          masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
          fHistInvMassS->Fill(masses[0]);
          fHistInvMassS->Fill(masses[1]);
        } // if (fCandidateType == kDstartoKpipi)
        
        fHistStat->Fill(3);
      } 
      // for MC with requirement particle level fill with AliAODMCParticle
      else { 
        new ((*fCandidateArray)[iCand]) AliAODMCParticle(*charmPart);
        fHistStat->Fill(3);
      } // if (fUseReco)
      	 
      iCand++;	
    } // if (fUseMCInfo && isMCBkg)

    
    if (fCandidateType == kDstartoKpipi) { //D*->D0pi->Kpipi
      ProcessDstar(dstar, mcLabel);
    }
    else if (fCandidateType==kD0toKpi) { //D0->Kpi
      ProcessD0(charmCand, isSelected, mcLabel);
    }
  } // end of D cand loop

  AliDebug(2, "Loop done");
  
  fHistNCandEv->Fill(fCandidateArray->GetEntriesFast());
  if (fCandidateType == kDstartoKpipi || fUseMCInfo) {
    Int_t nsbcand = fSideBandArray->GetEntriesFast();
    fHistStat->Fill(4, nsbcand);
    fHistNSBCandEv->Fill(nsbcand);
  }

  PostData(1, fOutput);
  PostData(3, fCandidateArray);
  PostData(4, fSideBandArray);

  AliDebug(2, "Exiting method");
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
   
  Info("Terminate"," terminate");
  AliAnalysisTaskSE::Terminate();
   
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    AliError("fOutput not available");
    return;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t range, Int_t pdg)
{
  //
  // AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits
  //
   
  Float_t mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();
   
  // compute the Delta mass for the D*
  if (fCandidateType==kDstartoKpipi) mass -= TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
  fMinMass = mass - range;
  fMaxMass = mass + range;
   
  AliInfo(Form("Setting mass limits to %f, %f", fMinMass, fMaxMass));
  if ((fMinMass<0.) || (fMaxMass<=0.) || (fMaxMass<fMinMass)) {
    AliError(Form("Wrong mass limits! Task '%s' will be disabled!", GetName()));
    fInhibitTask = kTRUE;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits(Double_t lowlimit, Double_t uplimit)
{
  //
  // AliAnalysisTaskSEDmesonsFilterCJ::SetMassLimits
  //
   
  if (uplimit>lowlimit) {
    fMinMass = lowlimit;
    fMaxMass = uplimit;
  }
  else {
    printf("Error! Lower limit larger than upper limit!\n");
    fMinMass = uplimit - uplimit*0.2;
    fMaxMass = uplimit;
  }
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar(Int_t nptbins, Float_t *width)
{
  //
  // AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar
  //
   
  if (nptbins > 30) {
    AliWarning("Maximum number of bins allowed is 30!");
    return kFALSE;
  }
   
  if (!width) return kFALSE;
  for (Int_t ipt=0; ipt<nptbins; ipt++) fSigmaD0[ipt] = width[ipt];
   
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis()
{
  //
  // AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis
  //
   
  // Statistics 
  fHistStat = new TH1I("fHistStat", "Statistics",6,-0.5,5.5);
  fHistStat->GetXaxis()->SetBinLabel(1, "N ev anal");
  fHistStat->GetXaxis()->SetBinLabel(2, "N ev sel");
  if (fUseMCInfo) {
    if (fUseReco) {
      fHistStat->GetXaxis()->SetBinLabel(3, "N D");
    }
    else {
      fHistStat->GetXaxis()->SetBinLabel(3, "N Gen D");
    }
  }
  else {
    fHistStat->GetXaxis()->SetBinLabel(3, "N Cand");
  }
  fHistStat->GetXaxis()->SetBinLabel(4, "N cand sel cuts");
  if (fCandidateType == kDstartoKpipi) {
    fHistStat->GetXaxis()->SetBinLabel(5, "N side band cand");
  }
  if (fUseMCInfo) {
    fHistStat->GetXaxis()->SetBinLabel(6, "N Background");
  }
  fHistStat->SetNdivisions(1);
  fOutput->Add(fHistStat);
   
  fHistNCandEv = new TH1F("fHistNCandEv", "Number of candidates per event (after cuts);# cand/ev", 100, 0., 100.);
  fOutput->Add(fHistNCandEv);

  if ((fCandidateType == kDstartoKpipi) || fUseMCInfo) {
    fHistNSBCandEv = new TH1F("hnSBCandEv", "Number of side bands candidates per event (after cuts);# cand/ev", 100, 0.,100.);
    fOutput->Add(fHistNSBCandEv);
  }

  if (fCandidateType == kDstartoKpipi) {
    fHistPtPion = new TH1F("fHistPtPion", "Primary pions candidates pt", 500, 0., 10.);
    fHistPtPion->SetStats(kTRUE);
    fHistPtPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPtPion->GetYaxis()->SetTitle("entries");
    fOutput->Add(fHistPtPion);
  }
   
  // Invariant mass related histograms
  const Int_t nbinsmass = 200;
  const Int_t ptbinsD = 100;
  Float_t ptmin =0.;
  Float_t ptmax =50.;
  fHistInvMassPtD = new TH2F("fHistInvMassPtD", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass, ptbinsD, ptmin, ptmax);
  fHistInvMassPtD->SetStats(kTRUE);
  fHistInvMassPtD->GetXaxis()->SetTitle("mass (GeV/c)");
  fHistInvMassPtD->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fOutput->Add(fHistInvMassPtD);

  fHistImpParS = new TH2F("fHistImpParS", "Impact parameter of daughter tracks; Getd0Prong();#it{p}^{daugh}_{T} (GeV/c)",200, -0.1,0.1,ptbinsD, ptmin, ptmax); //same range of pt of the D, but pt daughter used
  fOutput->Add(fHistImpParS);

  fHistInvMassS = new TH1F("fHistInvMassS", "D invariant mass distribution (filled with fCandidateArray)", nbinsmass, fMinMass, fMaxMass);
  fHistInvMassS->SetStats(kTRUE);
  fHistInvMassS->GetXaxis()->SetTitle("mass (GeV/c)");
  fHistInvMassS->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fOutput->Add(fHistInvMassS);

  const Int_t nbinsalpha = 200;
  Float_t minalpha = -TMath::Pi();
  Float_t maxalpha = TMath::Pi();
  const Int_t nbinsdeltaR = 200;
  Float_t mindeltaR = 0.;
  Float_t maxdeltaR = 10.;
  
  if (fUseMCInfo) {
    fHistImpParB = new TH2F("fHistImpParB", "Impact parameter of daughter tracks (Background); Getd0Prong();#it{p}^{daugh}_{T} (GeV/c)",200, -0.1,0.1,ptbinsD, ptmin, ptmax); //same range of pt of the D, but pt daughter used
    fOutput->Add(fHistImpParB);
    
    fHistInvMassB = new TH1F("fHistInvMassB", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass);
    fHistInvMassB->SetStats(kTRUE);
    fHistInvMassB->GetXaxis()->SetTitle("mass (GeV/c)");
    fHistInvMassB->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fOutput->Add(fHistInvMassB);
    
    if (fCandidateType == kDstartoKpipi) {
      
      fHistAlphaDDS    = new TH2F("fHistAlphaDDS","Angle D^{*}-D^{0} (Signal);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpisS  = new TH2F("fHistAlphaDpisS","Angle D^{*}-#pi_{soft} (Signal);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiS   = new TH2F("fHistAlphaDpiS","Angle D^{*}-#pi (Signal);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKS    = new TH2F("fHistAlphaDKS","Angle D^{*}-K (Signal);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistAlphaDDB    = new TH2F("fHistAlphaDDB","Angle D^{*}-D^{0} (Background);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpisB  = new TH2F("fHistAlphaDpisB","Angle D^{*}-#pi_{soft} (Background);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiB   = new TH2F("fHistAlphaDpiB","Angle D^{*}-#pi (Background);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKB    = new TH2F("fHistAlphaDKB","Angle D^{*}-K (Background);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      	 
      fHistDeltaRDDS   = new TH2F("fHistDeltaRDDS","Angle D^{*}-D^{0} (Signal);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpisS = new TH2F("fHistDeltaRDpisS","Angle D^{*}-#pi_{soft} (Signal);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiS  = new TH2F("fHistDeltaRDpiS","Angle D^{*}-#pi (Signal);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKS   = new TH2F("fHistDeltaRDKS","Angle D^{*}-K (Signal);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fHistDeltaRDDB   = new TH2F("fHistDeltaRDDB","Angle D^{*}-D^{0} (Background);#varphi (D^{*}) - #varphi (D0);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpisB = new TH2F("fHistDeltaRDpisB","Angle D^{*}-#pi_{soft} (Background);#varphi (D^{*}) - #varphi (#pi_{soft});p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiB  = new TH2F("fHistDeltaRDpiB","Angle D^{*}-#pi (Background);#varphi (D^{*}) - #varphi (#pi);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKB   = new TH2F("fHistDeltaRDKB","Angle D^{*}-K (Background);#varphi (D^{*}) - #varphi (K);p_{T}^{D*}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fOutput->Add(fHistAlphaDDS);
      fOutput->Add(fHistAlphaDpisS);
      fOutput->Add(fHistAlphaDpiS);
      fOutput->Add(fHistAlphaDKS);
      fOutput->Add(fHistAlphaDDB);
      fOutput->Add(fHistAlphaDpisB);
      fOutput->Add(fHistAlphaDpiB);
      fOutput->Add(fHistAlphaDKB);

      fOutput->Add(fHistDeltaRDDS);
      fOutput->Add(fHistDeltaRDpisS);
      fOutput->Add(fHistDeltaRDpiS);
      fOutput->Add(fHistDeltaRDKS);
      fOutput->Add(fHistDeltaRDDB);
      fOutput->Add(fHistDeltaRDpisB);
      fOutput->Add(fHistDeltaRDpiB);
      fOutput->Add(fHistDeltaRDKB);
    }
    else if (fCandidateType == kD0toKpi) {
      	 
      fHistAlphaDpiS  = new TH2F("fHistAlphaDpiS","Angle D^{0}-#pi (Signal);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKS   = new TH2F("fHistAlphaDKS","Angle D^{0}-K (Signal);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDpiR  = new TH2F("fHistAlphaDpiR","Angle D^{0}-#pi (Reflections);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKR   = new TH2F("fHistAlphaDKR","Angle D^{0}-K (Reflections);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      	 
      fHistAlphaDpiB  = new TH2F("fHistAlphaDpiB","Angle D^{0}-#pi (Background);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);
      fHistAlphaDKB   = new TH2F("fHistAlphaDKB","Angle D^{0}-K (Background);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsalpha, minalpha, maxalpha, ptbinsD, ptmin, ptmax);

      fHistDeltaRDpiS = new TH2F("fHistDeltaRDpiS","Angle D^{0}-#pi (Signal);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKS  = new TH2F("fHistDeltaRDKS","Angle D^{0}-K (Signal);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDpiR = new TH2F("fHistDeltaRDpiR","Angle D^{0}-#pi (Reflections);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKR  = new TH2F("fHistDeltaRDKR","Angle D^{0}-K (Reflections);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      	 
      fHistDeltaRDpiB = new TH2F("fHistDeltaRDpiB","Angle D^{0}-#pi (Background);#varphi (D^{0}) - #varphi (#pi);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);
      fHistDeltaRDKB  = new TH2F("fHistDeltaRDKB","Angle D^{0}-K (Background);#varphi (D^{0}) - #varphi (K);p_{T}^{D0}",nbinsdeltaR, mindeltaR, maxdeltaR, ptbinsD, ptmin, ptmax);

      fOutput->Add(fHistAlphaDpiS);
      fOutput->Add(fHistAlphaDKS);
      fOutput->Add(fHistAlphaDpiR);
      fOutput->Add(fHistAlphaDKR);
      fOutput->Add(fHistAlphaDpiB);
      fOutput->Add(fHistAlphaDKB);

      fOutput->Add(fHistDeltaRDpiS);
      fOutput->Add(fHistDeltaRDKS);
      fOutput->Add(fHistDeltaRDpiR);
      fOutput->Add(fHistDeltaRDKR);
      fOutput->Add(fHistDeltaRDpiB);
      fOutput->Add(fHistDeltaRDKB);
    }
  }
  
  return kTRUE;
}

//_______________________________________________________________________________
Float_t AliAnalysisTaskSEDmesonsFilterCJ::DeltaR(AliVParticle *p1, AliVParticle *p2) const
{
  //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
   
  if(!p1 || !p2) return -1;
  Double_t phi1=p1->Phi(),eta1=p1->Eta();
  Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;
   
  Double_t dPhi=phi1-phi2;
  if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
  if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());
   
  Double_t dEta=eta1-eta2;
  Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
  return deltaR;
}
