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
// A.Grelli (Utrecht University) a.grelli@uu.nl
// Xiaoming Zhang (LBNL)  XMZhang@lbl.gov
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include "TROOT.h"
#include <TH3F.h>

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
fOutput(0),
fCuts(0),
fMinMass(0.),
fMaxMass(0.),
fCandidateArray(0),
fSideBandArray(0)

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
fOutput(0),
fCuts(cuts),
fMinMass(0.),
fMaxMass(0.),
fCandidateArray(0),
fSideBandArray(0)
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
      } else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
   default :
      printf("%d not accepted!!\n",fCandidateType);
      break;
   }
   
   if (fCandidateType==kD0toKpi) SetMassLimits(0.15, fPDGmother);
   if (fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
   
   DefineOutput(1, TList::Class());       // histos
   DefineOutput(2, AliRDHFCuts::Class()); // my cuts
   DefineOutput(3, TClonesArray::Class()); //array of candidates
   DefineOutput(4, TClonesArray::Class()); //array of SB candidates
}

//_______________________________________________________________________________

AliAnalysisTaskSEDmesonsFilterCJ::~AliAnalysisTaskSEDmesonsFilterCJ()
{
   //
   // destructor
   //
   
   Info("~AliAnalysisTaskSEDmesonsFilterCJ","Calling Destructor");  
   
   if (fOutput) { delete fOutput; fOutput = 0; }
   if (fCuts)   { delete fCuts;   fCuts   = 0; }
   if (fCandidateArray)  { delete fCandidateArray;  fCandidateArray  = 0; }
   delete fSideBandArray;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskSEDmesonsFilterCJ::Init()
{
   //
   // Initialization
   //
   
   if(fDebug>1) printf("AnalysisTaskSEDmesonsForJetCorrelations::Init() \n");
   
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
   default: return;
   }
   
   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskSEDmesonsFilterCJ::UserCreateOutputObjects()
{ 
   //
   // output 
   //
   
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   
   fOutput = new TList(); fOutput->SetOwner();
   DefineHistoForAnalysis(); // define histograms
   
   if (fCandidateType==kD0toKpi){
      fCandidateArray = new TClonesArray("AliAODRecoDecayHF",0);
      fSideBandArray = new TClonesArray("AliAODRecoDecayHF",0); 
   }
   
   if (fCandidateType==kDstartoKpipi) {
      fCandidateArray = new TClonesArray("AliAODRecoCascadeHF",0);
      fSideBandArray = new TClonesArray("AliAODRecoCascadeHF",0); 
   }
   
   fCandidateArray->SetOwner();
   fCandidateArray->SetName(Form("fCandidateArray%s%s",fCandidateName.Data(),fUseReco ? "rec" : "gen"));
   
   //this is used for the DStar side bands and MC!
   fSideBandArray->SetOwner();
   fSideBandArray->SetName(Form("fSideBandArray%s%s",fCandidateName.Data(),fUseReco ? "rec" : "gen"));
  
   PostData(1, fOutput);
   PostData(3, fCandidateArray);
   PostData(4, fSideBandArray);
 
   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskSEDmesonsFilterCJ::UserExec(Option_t *){
   //
   // user exec
   //
   // Load the event
   AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
   
   TClonesArray *arrayDStartoD0pi = 0;
   if (!aodEvent && AODEvent() && IsStandardAOD()) {
      
      // In case there is an AOD handler writing a standard AOD, use the AOD 
      // event in memory rather than the input (ESD) event.    
      aodEvent = dynamic_cast<AliAODEvent*>(AODEvent());
      
      // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
      if(aodHandler->GetExtensions()) {
      	 AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      	 AliAODEvent *aodFromExt = ext->GetAOD();
      	 arrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      }
   } else {
      arrayDStartoD0pi = (TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
   }
   
   if (!arrayDStartoD0pi) {
      AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
      return;
   } else {
      AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   
   }
   
   TClonesArray* mcArray = 0x0;
   if (fUseMCInfo) {
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
      	 printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      	 return;
      }
   }
   
   //Histograms
   TH1I* hstat = (TH1I*)fOutput->FindObject("hstat");
   TH1F* hnSBCandEv=(TH1F*)fOutput->FindObject("hnSBCandEv");
   TH1F* hnCandEv=(TH1F*)fOutput->FindObject("hnCandEv");
   TH2F* hInvMassptD = (TH2F*)fOutput->FindObject("hInvMassptD");
   
   TH1F* hPtPion=0x0;
   if (fCandidateType==kDstartoKpipi) hPtPion = (TH1F*)fOutput->FindObject("hPtPion");
   hstat->Fill(0);
   
   // fix for temporary bug in ESDfilter 
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
   
   //Event selection
   Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   //TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if (!iseventselected) return;
   hstat->Fill(1);
   
   const Int_t nD = arrayDStartoD0pi->GetEntriesFast();
   if(fUseReco) hstat->Fill(2, nD);
   
   //D* and D0 prongs needed to MatchToMC method
   Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
   Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi
   
   //D0 from D0 bar
   Int_t pdgdaughtersD0[2] = { 211, 321 };     // pi,K 
   Int_t pdgdaughtersD0bar[2] = { 321, 211 };  // K,pi 
   
   Int_t iCand =0;
   Int_t iSBCand=0;
   Int_t isSelected = 0;
   AliAODRecoDecayHF *charmCand = 0;
   AliAODRecoCascadeHF *dstar = 0;
   AliAODMCParticle *charmPart = 0;
   Bool_t isMCBkg=kFALSE;
   
   Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
   
   Int_t mcLabel = -9999;
   Int_t pdgMeson = 413;
   if (fCandidateType==kD0toKpi) pdgMeson = 421;
   
   //clear the TClonesArray from the previous event
   fCandidateArray->Clear();
   fSideBandArray->Clear();
   
   for (Int_t icharm=0; icharm<nD; icharm++) {   //loop over D candidates
      charmCand = (AliAODRecoDecayHF*)arrayDStartoD0pi->At(icharm); // D candidates
      if (!charmCand) continue;
      
      if (fCandidateType==kDstartoKpipi) dstar = (AliAODRecoCascadeHF*)charmCand;
      
      if (fUseMCInfo) { // Look in MC, try to simulate the z
      	 if (fCandidateType==kDstartoKpipi) {
      	    mcLabel = dstar->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      	 }
      	 
      	 if (fCandidateType==kD0toKpi) 
      	    mcLabel = charmCand->MatchToMC(421,mcArray,fNProngs,fPDGdaughters);
      	 
      	 if (mcLabel<=0) isMCBkg=kTRUE;
      	 else hstat->Fill(2);
      	 if (!isMCBkg) charmPart=(AliAODMCParticle*)mcArray->At(mcLabel);
      }
      
      Double_t ptD = charmCand->Pt();
      
      // region of interest + cuts
      if (!fCuts->IsInFiducialAcceptance(ptD,charmCand->Y(pdgMeson))) continue;    
      
      if(!fUseMCInfo && fCandidateType==kDstartoKpipi){
      	 //select by track cuts the side band candidates (don't want mass cut)
      	 isSelected = fCuts->IsSelected(charmCand,AliRDHFCuts::kTracks,aodEvent); 
      	 if (!isSelected) continue;
      	 //add a reasonable cut on the invariant mass (e.g. (+-2\sigma, +-10 \sigma), with \sigma = fSigmaD0[bin])
      	 Int_t bin = fCuts->PtBin(ptD);
      	 if(bin<0 || bin>=fCuts->GetNPtBins()) {
      	    AliError(Form("Pt %.3f out of bounds",ptD));
      	    continue;
      	 }
      	 //if data and Dstar from D0 side band
      	 if (((dstar->InvMassD0()<=(mPDGD0-3.*fSigmaD0[bin])) && (dstar->InvMassD0()>(mPDGD0-10.*fSigmaD0[bin]))) /*left side band*/||  ((dstar->InvMassD0()>=(mPDGD0+3.*fSigmaD0[bin])) && (dstar->InvMassD0()<(mPDGD0+10.*fSigmaD0[bin])))/*right side band*/){	
      	    
      	    new ((*fSideBandArray)[iSBCand]) AliAODRecoCascadeHF(*dstar);
      	    iSBCand++;
      	 }
      }
      //candidate selected by cuts and PID
      isSelected = fCuts->IsSelected(charmCand,AliRDHFCuts::kAll,aodEvent); //selected
      if (!isSelected) continue;
      
      //for data and MC signal fill fCandidateArray
      if(!fUseMCInfo || (fUseMCInfo && !isMCBkg)){
      	 // for data or MC with the requirement fUseReco fill with candidates
      	 if(fUseReco) {
      	    if (fCandidateType==kDstartoKpipi){
      	       new ((*fCandidateArray)[iCand]) AliAODRecoCascadeHF(*dstar);
      	       AliInfo(Form("Dstar delta mass = %f",dstar->DeltaInvMass()));
      	    } else{
      	       new ((*fCandidateArray)[iCand]) AliAODRecoDecayHF(*charmCand);
      	       //Printf("Filling reco");
      	    }      	    
      	    
      	    hstat->Fill(3);
      	 }
      	 // for MC with requirement particle level fill with AliAODMCParticle
      	 else if (fUseMCInfo) {
      	    new ((*fCandidateArray)[iCand]) AliAODMCParticle(*charmPart);
      	    //Printf("Filling gen");
      	    hstat->Fill(3);
      	 }
      	 
      	 iCand++;	
      }
      //for MC background fill fSideBandArray (which is instead filled above for DStar in case of data for the side bands candidates)
      else if(fUseReco){
      	 if (fCandidateType==kDstartoKpipi){
      	    new ((*fSideBandArray)[iSBCand]) AliAODRecoCascadeHF(*dstar);
      	 }
      	 if (fCandidateType==kD0toKpi){
      	    new ((*fSideBandArray)[iSBCand]) AliAODRecoDecayHF(*charmCand);
      	 }
      	 iSBCand++;
      }
      
      
      Double_t masses[2];
      if (fCandidateType==kDstartoKpipi) { //D*->D0pi->Kpipi
      	 
      	 //softpion from D* decay
      	 AliAODTrack *track2 = (AliAODTrack*)dstar->GetBachelor();  
      	 
      	 // select D* in the D0 window.
      	 // In the cut object window is loose to allow for side bands
      	 
      	 
      	 // retrieve the corresponding pt bin for the candidate
      	 // and set the expected D0 width (x3)
      	 //    static const Int_t n = fCuts->GetNPtBins();
      	 Int_t bin = fCuts->PtBin(ptD);
      	 if(bin<0 || bin>=fCuts->GetNPtBins()) {
      	    AliError(Form("Pt %.3f out of bounds",ptD));
      	    continue;
      	 }
      	 
      	 AliInfo(Form("Pt bin %d and sigma D0 %.4f",bin,fSigmaD0[bin]));
      	 //consider the Dstar candidates only if the mass of the D0 is in 3 sigma wrt the PDG value
      	 if ((dstar->InvMassD0()>=(mPDGD0-3.*fSigmaD0[bin])) && (dstar->InvMassD0()<=(mPDGD0+3.*fSigmaD0[bin]))) {	
      	    masses[0] = dstar->DeltaInvMass(); //D*
      	    masses[1] = 0.; //dummy for D*
      	    
      	    //D*  delta mass
      	    hInvMassptD->Fill(masses[0], ptD); // 2 D slice for pt bins
      	    
      	    // D* pt and soft pion pt for good candidates  	      	
      	    Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
      	    Double_t invmassDelta = dstar->DeltaInvMass();
      	    if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<0.0021) hPtPion->Fill(track2->Pt());
      	 }
      } //Dstar specific
      
      if (fCandidateType==kD0toKpi) { //D0->Kpi
      	 
      	 //needed quantities
      	 masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
      	 masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
      	 hstat->Fill(3);
      	 
      	 // mass vs pt
      	 if (isSelected==1 || isSelected==3) hInvMassptD->Fill(masses[0],ptD);
      	 if (isSelected>=2) hInvMassptD->Fill(masses[1],ptD);
      } //D0 specific
      
      charmCand = 0;
      isMCBkg=kFALSE;
   } // end of D cand loop
   
   hnCandEv->Fill(fCandidateArray->GetEntriesFast());
   if (fCandidateType==kDstartoKpipi || fUseMCInfo) {
      Int_t nsbcand=fSideBandArray->GetEntriesFast();
      hstat->Fill(4,nsbcand);
      hnSBCandEv->Fill(nsbcand);
   }
   Printf("N candidates selected %d, counter = %d",fCandidateArray->GetEntries(), iCand);
   PostData(1, fOutput);
   PostData(3, fCandidateArray);
   PostData(4, fSideBandArray);
  
   return;
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
      printf("ERROR: fOutput not available\n");
      return;
   }
   
   return;
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
   
   AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
   if ((fMinMass<0.) || (fMaxMass<=0.) || (fMaxMass<fMinMass)) AliFatal("Wrong mass limits!\n");
   
   return;
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
   } else {
      printf("Error! Lower limit larger than upper limit!\n");
      fMinMass = uplimit - uplimit*0.2;
      fMaxMass = uplimit;
   }
   
   return;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar(Int_t nptbins, Float_t *width)
{
   //
   // AliAnalysisTaskSEDmesonsFilterCJ::SetD0WidthForDStar
   //
   
   if (nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   
   if (!width) return kFALSE;
   for (Int_t ipt=0; ipt<nptbins; ipt++) fSigmaD0[ipt]=width[ipt];
   
   return kTRUE;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis()
{
   //
   // AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis
   //
   
   // Statistics 
   TH1I* hstat = new TH1I("hstat","Statistics",5,-0.5,4.5);
   hstat->GetXaxis()->SetBinLabel(1, "N ev anal");
   hstat->GetXaxis()->SetBinLabel(2, "N ev sel");
   if(fUseReco) hstat->GetXaxis()->SetBinLabel(3, "N cand");
   else hstat->GetXaxis()->SetBinLabel(3, "N Gen D");
   hstat->GetXaxis()->SetBinLabel(4, "N cand sel cuts");
   if (fCandidateType==kDstartoKpipi) hstat->GetXaxis()->SetBinLabel(5, "N side band cand");  
   hstat->SetNdivisions(1);
   fOutput->Add(hstat);
   
   TH1F* hnCandEv=new TH1F("hnCandEv", "Number of candidates per event (after cuts);# cand/ev", 100, 0.,100.);
   fOutput->Add(hnCandEv);
   
   // Invariant mass related histograms
   const Int_t nbinsmass = 200;
   TH2F *hInvMass = new TH2F("hInvMassptD", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass, 100, 0., 50.);
   hInvMass->SetStats(kTRUE);
   hInvMass->GetXaxis()->SetTitle("mass (GeV/c)");
   hInvMass->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   fOutput->Add(hInvMass);
   
   if (fCandidateType==kDstartoKpipi) {
      TH1F* hnSBCandEv=new TH1F("hnSBCandEv", "Number of side bands candidates per event (after cuts);# cand/ev", 100, 0.,100.);
      fOutput->Add(hnSBCandEv);
      
      TH1F* hPtPion = new TH1F("hPtPion", "Primary pions candidates pt", 500, 0., 10.);
      hPtPion->SetStats(kTRUE);
      hPtPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtPion->GetYaxis()->SetTitle("entries");
      fOutput->Add(hPtPion);
   }
   
   return kTRUE; 
}
