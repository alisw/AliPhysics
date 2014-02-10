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
//  Analysis Taks for reconstructed particle correlation 
//  (first implementation done for D mesons) with jets 
//  (use the so called Emcal framework)
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include "TROOT.h"
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskFlavourFilterAndJetCorrelations.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliPicoTrack.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"

ClassImp(AliAnalysisTaskFlavourFilterAndJetCorrelations)


//_______________________________________________________________________________

AliAnalysisTaskFlavourFilterAndJetCorrelations::AliAnalysisTaskFlavourFilterAndJetCorrelations() :
AliAnalysisTaskEmcalJet("",kFALSE),
fUseMCInfo(kTRUE), 
fUseReco(kTRUE),
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fmyOutput(0),
fmyOutputF(0),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fCandArrName(0),
fLeadingJetOnly(kFALSE),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0)
{
   //
   // Default ctor
   //
   //SetMakeGeneralHistograms(kTRUE);
   
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourFilterAndJetCorrelations::AliAnalysisTaskFlavourFilterAndJetCorrelations(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype) :
AliAnalysisTaskEmcalJet(name,kFALSE),
fUseMCInfo(kTRUE),
fUseReco(kTRUE),  
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fmyOutput(0),
fmyOutputF(0),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fCandArrName(0),
fLeadingJetOnly(kFALSE),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0)
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //
   
   Info("AliAnalysisTaskFlavourFilterAndJetCorrelations","Calling Constructor");
   fCuts=cuts;
   fCandidateType=candtype;
   const Int_t nptbins=fCuts->GetNPtBins();
   Float_t defaultSigmaD013[13]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040};
   
   switch(fCandidateType){
   case 0:
      fPDGmother=421;
      fNProngs=2;
      fPDGdaughters[0]=211;//pi 
      fPDGdaughters[1]=321;//K
      fPDGdaughters[2]=0; //empty
      fPDGdaughters[3]=0; //empty
      fBranchName="D0toKpi";
      fCandArrName="D0";
      break;
   case 1: 
      fPDGmother=413;
      fNProngs=3;
      fPDGdaughters[1]=211;//pi soft
      fPDGdaughters[0]=421;//D0
      fPDGdaughters[2]=211;//pi fromD0
      fPDGdaughters[3]=321; // K from D0
      fBranchName="Dstar";
      fCandArrName="Dstar";
      
      if(nptbins<=13){
      	 for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]= defaultSigmaD013[ipt];
      } else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
   default:
      printf("%d not accepted!!\n",fCandidateType);
      break;
   }
   
   if(fCandidateType==kD0toKpi)SetMassLimits(0.15,fPDGmother);
   if(fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
   
   DefineOutput(1,TList::Class()); // histos
   DefineOutput(2,AliRDHFCuts::Class()); // my cuts
   DefineOutput(3,TList::Class()); //histos filter
   
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourFilterAndJetCorrelations::~AliAnalysisTaskFlavourFilterAndJetCorrelations() {
   //
   // destructor
   //
   
   Info("~AliAnalysisTaskFlavourFilterAndJetCorrelations","Calling Destructor");  
   
   delete fmyOutput;
   delete fmyOutputF;
   delete fCuts;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::Init(){
   //
   // Initialization
   //
   
   if(fDebug > 1) printf("AnalysisTaskRecoJetCorrelations::Init() \n");
   switch(fCandidateType){
   case 0:
      {
      	 AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDzero");
      	 // Post the data
      	 PostData(2,copyfCuts);
      }
      break;
   case 1:
      {
      	 AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDStar");
      	 // Post the cuts
      	 PostData(2,copyfCuts);
      }
      break;
   default:
      return;
   }
   
   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::UserCreateOutputObjects() { 
   
   fCandidateArray = new TClonesArray("AliAODRecoDecayHF",0);
   fCandidateArray->SetName(Form("fCandidateArray%s%s",fCandArrName.Data(),fUseReco ? "rec" : "gen"));
   
   if (fCandidateType==kDstartoKpipi){
      fSideBandArray = new TClonesArray("AliAODRecoCascadeHF",0); //this is for the DStar only!
      fSideBandArray->SetName(Form("fSideBandArray%s%s",fCandArrName.Data(),fUseReco ? "rec" : "gen"));
   }
   
   // output 
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   fmyOutput = new TList();
   fmyOutput->SetOwner();
   fmyOutput->SetName("pippo");
   // define histograms
   DefineHistoForAnalysis();
   PostData(1,fmyOutput);

   fmyOutputF = new TList();
   fmyOutputF->SetOwner();
   fmyOutputF->SetName("pluto");
   // define histograms
   DefineHistoForFiltering();
   PostData(3,fmyOutputF);
   
   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::UserExec(Option_t *)
{
   // user exec
   if (!fInitialized){
      AliAnalysisTaskEmcalJet::ExecOnce();
   }
   
   // Load the event
   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
   
   TClonesArray *arrayDStartoD0pi=0;
   TClonesArray *trackArr = 0;
   
   if (!aodEvent && AODEvent() && IsStandardAOD()) {
      
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
      	 arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      }
   } else {
      arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
   }
   
   if (!arrayDStartoD0pi) {
      AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
      //  return;
   } else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   
   
   TClonesArray* mcArray = 0x0;
   if (fUseMCInfo) {
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
      	 printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      	 return;
      }
   }
   
   //retrieve jets
   trackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoTracks"));
   //clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClustersCorr"));
   //jetArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrName));
   fJetRadius=GetJetContainer(0)->GetJetRadius();
   
   if(!trackArr){
      AliInfo(Form("Could not find the track array\n"));
      return;
   }
   
   
   //FILTER D MESONS

   fCandidateArray->Clear();
   fSideBandArray->Clear();
   
   //Histograms
   TH1I* hstatF = (TH1I*)fmyOutputF->FindObject("hstatF");
   TH1F* hnSBCandEv=(TH1F*)fmyOutputF->FindObject("hnSBCandEv");
   TH1F* hnCandEv=(TH1F*)fmyOutputF->FindObject("hnCandEv");
   TH2F* hInvMassptDF = (TH2F*)fmyOutputF->FindObject("hInvMassptDF");
   
   TH1F* hPtPion=0x0;
   if (fCandidateType==kDstartoKpipi) hPtPion = (TH1F*)fmyOutputF->FindObject("hPtPion");
   hstatF->Fill(0);
   
   // fix for temporary bug in ESDfilter 
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
   
   //Event selection
   Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   //TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if (!iseventselected) return;
   hstatF->Fill(1);
   
   const Int_t nD = arrayDStartoD0pi->GetEntriesFast();
   if(fUseReco) hstatF->Fill(2, nD);
   
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
   AliAODMCParticle *charmPart = 0;
   Bool_t isMCBkg=kFALSE;
   
   Double_t mPDGD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
   
   Int_t mcLabel = -9999;
   Int_t pdgMeson = 413;
   if (fCandidateType==kD0toKpi) pdgMeson = 421;
   
   for (Int_t icharm=0; icharm<nD; icharm++) {   //loop over D candidates
      charmCand = (AliAODRecoDecayHF*)arrayDStartoD0pi->At(icharm); // D candidates
      if (!charmCand) continue;
      
      
      if (fUseMCInfo) { // Look in MC, try to simulate the z
      	 if (fCandidateType==kDstartoKpipi) {
      	    AliAODRecoCascadeHF *temp = (AliAODRecoCascadeHF*)charmCand;
      	    mcLabel = temp->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      	 }
      	 
      	 if (fCandidateType==kD0toKpi) 
      	    mcLabel = charmCand->MatchToMC(421,mcArray,fNProngs,fPDGdaughters);
      	 
      	 if (mcLabel<=0) isMCBkg=kTRUE;
      	 else hstatF->Fill(2);
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
      	 AliAODRecoCascadeHF *temp = (AliAODRecoCascadeHF*)charmCand;
      	 //if data and Dstar from D0 side band
      	 if (((temp->InvMassD0()<=(mPDGD0-3.*fSigmaD0[bin])) && (temp->InvMassD0()>(mPDGD0-10.*fSigmaD0[bin]))) /*left side band*/||  ((temp->InvMassD0()>=(mPDGD0+3.*fSigmaD0[bin])) && (temp->InvMassD0()<(mPDGD0+10.*fSigmaD0[bin])))/*right side band*/){	
      	    
      	    new ((*fSideBandArray)[iSBCand]) AliAODRecoCascadeHF(*temp);
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
      	    new ((*fCandidateArray)[iCand]) AliAODRecoDecayHF(*charmCand);
      	    //Printf("Filling reco");
      	    hstatF->Fill(3);
      	 }
      	 // for MC with requirement particle level fill with AliAODMCParticle
      	 else if (fUseMCInfo) {
      	    new ((*fCandidateArray)[iCand]) AliAODMCParticle(*charmPart);
      	    //Printf("Filling gen");
      	    hstatF->Fill(3);
      	 }
      	 
      	 iCand++;	
      }
      //for MC background fill fSideBandArray (which is instead filled above for DStar in case of data for the side bands candidates)
      else if(fUseReco){
      	 new ((*fSideBandArray)[iSBCand]) AliAODRecoDecayHF(*charmCand);
      	 iSBCand++;
      }
      
      
      Double_t masses[2];
      if (fCandidateType==kDstartoKpipi) { //D*->D0pi->Kpipi
      	 
      	 //softpion from D* decay
      	 AliAODRecoCascadeHF *temp = (AliAODRecoCascadeHF*)charmCand;
      	 AliAODTrack *track2 = (AliAODTrack*)temp->GetBachelor();  
      	 
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
      	 if ((temp->InvMassD0()>=(mPDGD0-3.*fSigmaD0[bin])) && (temp->InvMassD0()<=(mPDGD0+3.*fSigmaD0[bin]))) {	
      	    masses[0] = temp->DeltaInvMass(); //D*
      	    masses[1] = 0.; //dummy for D*
      	    
      	    //D*  delta mass
      	    hInvMassptDF->Fill(masses[0], ptD); // 2 D slice for pt bins
      	    
      	    // D* pt and soft pion pt for good candidates  	      	
      	    Double_t mPDGDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
      	    Double_t invmassDelta = temp->DeltaInvMass();
      	    if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<0.0021) hPtPion->Fill(track2->Pt());
      	 }
      } //Dstar specific
      
      if (fCandidateType==kD0toKpi) { //D0->Kpi
      	 
      	 //needed quantities
      	 masses[0] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0); //D0
      	 masses[1] = charmCand->InvMass(fNProngs, (UInt_t*)pdgdaughtersD0bar); //D0bar
      	 hstatF->Fill(3);
      	 
      	 // mass vs pt
      	 if (isSelected==1 || isSelected==3) hInvMassptDF->Fill(masses[0],ptD);
      	 if (isSelected>=2) hInvMassptDF->Fill(masses[1],ptD);
      } //D0 specific
      
      charmCand = 0;
      isMCBkg=kFALSE;
   } // end of D cand loop
   
   hnCandEv->Fill(fCandidateArray->GetEntriesFast());
   if (fCandidateType==kDstartoKpipi || fUseMCInfo) {
      Int_t nsbcand=fSideBandArray->GetEntriesFast();
      hstatF->Fill(4,nsbcand);
      hnSBCandEv->Fill(nsbcand);
   }
   
   //CORRELATION WITH JETS
   /*
   TString arrname="fCandidateArray";
   TString candarrname=Form("%s%s%s",arrname.Data(),fCandArrName.Data(),fUseReco ? "rec" : "gen");
   fCandidateArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(candarrname));
   if (!fCandidateArray) {
      Printf("%s array not found",candarrname.Data());
      InputEvent()->GetList()->ls();
      return;
   }
   //Printf("ncandidates found %d",fCandidateArray->GetEntriesFast());
   
   TString arrSBname="fSideBandArray";
   TString sbarrname=Form("%s%s%s",arrSBname.Data(),fCandArrName.Data(),fUseReco ? "rec" : "gen");
   fSideBandArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(sbarrname));
   if (fCandidateType==1 && !fSideBandArray) {
      Printf("%s array not found",sbarrname.Data());
      InputEvent()->GetList()->ls();
      return;
   }
   //Printf("nSBCand or Bkg found %d",fSideBandArray->GetEntriesFast());
   */
   
   //Histograms
   TH1I* hstat=(TH1I*)fmyOutput->FindObject("hstat");
   TH1F* hPtJetTrks=(TH1F*)fmyOutput->FindObject("hPtJetTrks");
   TH1F* hPhiJetTrks=(TH1F*)fmyOutput->FindObject("hPhiJetTrks");
   TH1F* hEtaJetTrks=(TH1F*)fmyOutput->FindObject("hEtaJetTrks");
   TH1F* hEjetTrks=(TH1F*)fmyOutput->FindObject("hEjetTrks");
   TH1F* hPtJet=(TH1F*)fmyOutput->FindObject("hPtJet");
   TH1F* hPhiJet=(TH1F*)fmyOutput->FindObject("hPhiJet");
   TH1F* hEtaJet=(TH1F*)fmyOutput->FindObject("hEtaJet");
   TH1F* hEjet=(TH1F*)fmyOutput->FindObject("hEjet");
   TH1F* hNtrArr=(TH1F*)fmyOutput->FindObject("hNtrArr");
   TH1F* hNJetPerEv=(TH1F*)fmyOutput->FindObject("hNJetPerEv");
   TH1F* hdeltaRJetTracks=(TH1F*)fmyOutput->FindObject("hdeltaRJetTracks");
   TH1F* hNDPerEvNoJet=(TH1F*)fmyOutput->FindObject("hNDPerEvNoJet");
   TH1F* hptDPerEvNoJet=(TH1F*)fmyOutput->FindObject("hptDPerEvNoJet");
   TH1F* hNJetPerEvNoD=(TH1F*)fmyOutput->FindObject("hNJetPerEvNoD");
   TH1F* hPtJetPerEvNoD=(TH1F*)fmyOutput->FindObject("hPtJetPerEvNoD");
   THnSparseF* hnspDstandalone=(THnSparseF*)fmyOutput->FindObject("hsDstandalone");
       
   hstat->Fill(0);
   
   // fix for temporary bug in ESDfilter 
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
   
   //Event selection
   //Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if(!iseventselected) return;
   
   hstat->Fill(1);

   //retrieve charm candidates selected
   Int_t candidates = fCandidateArray->GetEntriesFast();
  
   //trigger on jets
   
   Int_t njets=GetJetContainer()->GetNJets();
   if(njets == 0) {
      hstat->Fill(6, candidates);
      hNDPerEvNoJet->Fill(candidates);
      for(Int_t iD=0;iD<candidates;iD++){
      	 AliVParticle* cand=(AliVParticle*)fCandidateArray->At(iD);
      	 if(!cand) continue;
      	 hptDPerEvNoJet->Fill(cand->Pt());
      
      }
      return;
      
   }
   
   //loop on candidates standalone (checking the candidates are there and their phi-eta distributions)
   
   for(Int_t ic = 0; ic < candidates; ic++) {
      
      // D* candidates
      AliVParticle* charm=0x0;
      charm=(AliVParticle*)fCandidateArray->At(ic);
      if(!charm) continue;
      hstat->Fill(2);
      
      Double_t candsparse[4]={charm->Eta(), charm->Phi(), charm->Pt(), 0};
      
      if(fCandidateType==kDstartoKpipi) {
      	 AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)charm;
      	 Double_t deltamass= dstar->DeltaInvMass();
      	 candsparse[3]=deltamass;
      	 hnspDstandalone->Fill(candsparse);
      }
      if(fCandidateType==kD0toKpi){
      	 AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)charm;
      	 Int_t isselected=fCuts->IsSelected(dzero,AliRDHFCuts::kAll,aodEvent);
      	 
      	 Double_t masses[2];
      	 //Int_t pdgdaughtersD0[2]={211,321};//pi,K 
      	 //Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      	 
      	 masses[0]=dzero->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
      	 masses[1]=dzero->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
      	 if(isselected==1 || isselected==3) {
      	    candsparse[3]=masses[0];
      	    hnspDstandalone->Fill(candsparse);
      	 }
      	 if(isselected>=2){
      	    candsparse[3]=masses[1];
      	    hnspDstandalone->Fill(candsparse);
    	    
      	 }
      }
   }

   // we start with jets
   Double_t ejet   = 0;
   Double_t phiJet = 0;
   Double_t etaJet = 0;
   Double_t ptjet = 0;
   Double_t leadingJet =0;
   
   Int_t ntrarr=trackArr->GetEntriesFast();
   hNtrArr->Fill(ntrarr);
   
   for(Int_t i=0;i<ntrarr;i++){
      AliVTrack *jtrack=static_cast<AliVTrack*>(trackArr->At(i));
      //reusing histograms
      hPtJetTrks->Fill(jtrack->Pt());
      hPhiJetTrks->Fill(jtrack->Phi());
      hEtaJetTrks->Fill(jtrack->Eta());
      hEjetTrks->Fill(jtrack->E());
   }
   
   
   //option to use only the leading jet
   if(fLeadingJetOnly){
      for (Int_t iJetsL = 0; iJetsL<njets; iJetsL++) {    
      	 AliEmcalJet* jetL = (AliEmcalJet*)GetJetFromArray(iJetsL);
      	 ptjet   = jetL->Pt();
      	 if(ptjet>leadingJet ) leadingJet = ptjet;
      	 
      }
   }
   
   Int_t cntjet=0;
   //loop over jets and consider only the leading jet in the event
   for (Int_t iJets = 0; iJets<njets; iJets++) {    
      //Printf("Jet N %d",iJets);
      AliEmcalJet* jet = (AliEmcalJet*)GetJetFromArray(iJets);
      //Printf("Corr task Accept Jet");
      if(!AcceptJet(jet)) {
      	 hstat->Fill(5);
      	 continue;
      }
      
      //jets variables
      ejet   = jet->E();
      phiJet = jet->Phi();
      etaJet = jet->Eta();
      ptjet = jet->Pt();
      
      // choose the leading jet
      if(fLeadingJetOnly && (ejet<leadingJet)) continue;
      //used for normalization
      hstat->Fill(3);
      cntjet++;
      // fill energy, eta and phi of the jet
      hEjet   ->Fill(ejet);
      hPhiJet ->Fill(phiJet);
      hEtaJet ->Fill(etaJet);
      hPtJet  ->Fill(ptjet);
      
      //loop on jet particles
      Int_t ntrjet=  jet->GetNumberOfTracks();    
      for(Int_t itrk=0;itrk<ntrjet;itrk++){
      	 
      	 AliPicoTrack* jetTrk=(AliPicoTrack*)jet->TrackAt(itrk,trackArr);     
      	 hdeltaRJetTracks->Fill(DeltaR(jet,jetTrk));
      	 
      }//end loop on jet tracks
      
      if(candidates==0){
      	 hstat->Fill(7);
      	 hPtJetPerEvNoD->Fill(jet->Pt());
      }
      
      //Printf("N candidates %d ", candidates);
      for(Int_t ic = 0; ic < candidates; ic++) {
      	 
      	 // D* candidates
      	 AliVParticle* charm=0x0;
      	 charm=(AliVParticle*)fCandidateArray->At(ic);
      	 if(!charm) continue;
      	 
      	 
      	 FlagFlavour(charm, jet);
      	 if (jet->TestFlavourTag(AliEmcalJet::kDStar)) hstat->Fill(4);
      	 
      	 FillHistogramsRecoJetCorr(charm, jet);
      	 
      }
      //retrieve side band background candidates for Dstar
      Int_t nsbcand = 0; if(fSideBandArray) nsbcand=fSideBandArray->GetEntriesFast();
      
      for(Int_t ib=0;ib<nsbcand;ib++){
      	 if(fCandidateType==kDstartoKpipi && !fUseMCInfo){
      	    AliAODRecoCascadeHF *sbcand=(AliAODRecoCascadeHF*)fSideBandArray->At(ib);
      	    if(!sbcand) continue;
      	    SideBandBackground(sbcand,jet);
      	 }
      	 if(fUseMCInfo){
      	    AliAODRecoDecayHF* charmbg = 0x0;
      	    charmbg=(AliAODRecoDecayHF*)fCandidateArray->At(ib);
      	    if(!charmbg) continue;
      	    MCBackground(charmbg,jet);
      	 }
      }
   } // end of jet loop
   
   hNJetPerEv->Fill(cntjet);
   if(candidates==0) hNJetPerEvNoD->Fill(cntjet);
   PostData(1,fmyOutput);
   PostData(3,fmyOutputF);
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::Terminate(Option_t*)
{    
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
   Info("Terminate"," terminate");
   AliAnalysisTaskSE::Terminate();
   
   fmyOutput = dynamic_cast<TList*> (GetOutputData(1));
   if (!fmyOutput) {     
      printf("ERROR: fmyOutput not available\n");
      return;
   }
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourFilterAndJetCorrelations::SetMassLimits(Double_t range, Int_t pdg){
   Float_t mass=0;
   Int_t abspdg=TMath::Abs(pdg);
   
   mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
   // compute the Delta mass for the D*
   if(fCandidateType==kDstartoKpipi){
      Float_t mass1=0;
      mass1=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      mass = mass-mass1;
   }
   
   fMinMass = mass-range;
   fMaxMass = mass+range;
   
   AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
   if (fMinMass<0 || fMaxMass<=0 || fMaxMass<fMinMass) AliFatal("Wrong mass limits!\n");
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourFilterAndJetCorrelations::SetMassLimits(Double_t lowlimit, Double_t uplimit){
   if(uplimit>lowlimit)
   {
      fMinMass = lowlimit;
      fMaxMass = uplimit;
   }
   else{
      printf("Error! Lower limit larger than upper limit!\n");
      fMinMass = uplimit - uplimit*0.2;
      fMaxMass = uplimit;
   }
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourFilterAndJetCorrelations::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________

Double_t AliAnalysisTaskFlavourFilterAndJetCorrelations::Z(AliVParticle* part,AliEmcalJet* jet) const{
   if(!part || !jet){
      printf("Particle or jet do not exist!\n");
      return -99;
   }
   Double_t p[3],pj[3];
   Bool_t okpp=part->PxPyPz(p);
   Bool_t okpj=jet->PxPyPz(pj);
   if(!okpp || !okpj){
      printf("Problems getting momenta\n");
      return -999;
   }
   Double_t pjet=jet->P();
   Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet*pjet);
   return z;
}

//_______________________________________________________________________________

Bool_t  AliAnalysisTaskFlavourFilterAndJetCorrelations::DefineHistoForAnalysis(){
   
   // Statistics 
   TH1I* hstat=new TH1I("hstat","Statistics",8,-0.5,7.5);
   hstat->GetXaxis()->SetBinLabel(1,"N ev anal");
   hstat->GetXaxis()->SetBinLabel(2,"N ev sel");
   hstat->GetXaxis()->SetBinLabel(3,"N cand sel & jet");
   hstat->GetXaxis()->SetBinLabel(4,"N jets");
   hstat->GetXaxis()->SetBinLabel(5,"N cand in jet");
   hstat->GetXaxis()->SetBinLabel(6,"N jet rej");
   hstat->GetXaxis()->SetBinLabel(7,"N cand sel & !jet");
   hstat->GetXaxis()->SetBinLabel(8,"N jets & !D");
   hstat->SetNdivisions(1);
   fmyOutput->Add(hstat);
   
   const Int_t nbinsmass=200;
   const Int_t nbinsptjet=500;
   const Int_t nbinsptD=100;
   const Int_t nbinsz=100;
   const Int_t nbinsphi=200;
   const Int_t nbinseta=100;
   
   const Float_t ptjetlims[2]={0.,200.};
   const Float_t ptDlims[2]={0.,50.};
   const Float_t zlims[2]={0.,1.2};
   const Float_t philims[2]={0.,6.3};
   const Float_t etalims[2]={-1.5,1.5};
   
   if(fCandidateType==kDstartoKpipi) 
   {
      
      TH2F* hDiffSideBand = new TH2F("hDiffSideBand","M(kpipi)-M(kpi) Side Band Background",nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      hDiffSideBand->SetStats(kTRUE);
      hDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV");
      hDiffSideBand->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      hDiffSideBand->Sumw2();
      fmyOutput->Add(hDiffSideBand); 
      
      
      TH1F* hPtPion = new TH1F("hPtPion","Primary pions candidates pt ",500,0,10);
      hPtPion->SetStats(kTRUE);
      hPtPion->GetXaxis()->SetTitle("GeV/c");
      hPtPion->GetYaxis()->SetTitle("Entries");
      hPtPion->Sumw2();
      fmyOutput->Add(hPtPion);
      
   }
   // D related histograms
   TH1F *hNDPerEvNoJet=new TH1F("hNDPerEvNoJet","Number of candidates per event with no jets; N candidate/ev with no jet", 20, 0., 20.);
   hNDPerEvNoJet->Sumw2();
   fmyOutput->Add(hNDPerEvNoJet);
   
   TH1F *hptDPerEvNoJet=new TH1F("hptDPerEvNoJet","pt distribution of candidates per events with no jets; p_{t}^{D} (GeV/c)",nbinsptD, ptDlims[0],ptDlims[1]);
   hptDPerEvNoJet->Sumw2();
   fmyOutput->Add(hptDPerEvNoJet);
   
   const Int_t    nAxisD=4;
   const Int_t    nbinsSparseD[nAxisD]={nbinsphi,nbinseta,nbinsptD,nbinsmass};
   const Double_t minSparseD[nAxisD]  ={philims[0],etalims[0],ptDlims[0],fMinMass};
   const Double_t maxSparseD[nAxisD]  ={philims[1],etalims[1],ptDlims[1],fMaxMass};
   THnSparseF *hsDstandalone=new THnSparseF("hsDstandalone","#phi, #eta, p_{T}^{D}, and mass", nAxisD, nbinsSparseD, minSparseD, maxSparseD);
   hsDstandalone->Sumw2();
   
   fmyOutput->Add(hsDstandalone);
   
   // jet related fistograms
   
   TH1F* hEjetTrks      = new TH1F("hEjetTrks",  "Jet tracks energy distribution;Energy (GeV)",500,0,200);
   hEjetTrks->Sumw2();
   TH1F* hPhiJetTrks    = new TH1F("hPhiJetTrks","Jet tracks #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   hPhiJetTrks->Sumw2();
   TH1F* hEtaJetTrks    = new TH1F("hEtaJetTrks","Jet tracks #eta distribution; #eta",  nbinseta,etalims[0],etalims[1]);
   hEtaJetTrks->Sumw2();
   TH1F* hPtJetTrks     = new TH1F("hPtJetTrks",  "Jet tracks Pt distribution; p_{T} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
   hPtJetTrks->Sumw2();
   
   TH1F* hEjet      = new TH1F("hEjet",  "Jet energy distribution;Energy (GeV)",500,0,200);
   hEjet->Sumw2();
   TH1F* hPhiJet    = new TH1F("hPhiJet","Jet #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   hPhiJet->Sumw2();
   TH1F* hEtaJet    = new TH1F("hEtaJet","Jet #eta distribution; #eta", nbinseta,etalims[0],etalims[1]);
   hEtaJet->Sumw2();
   TH1F* hPtJet      = new TH1F("hPtJet",  "Jet Pt distribution; p_{T} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
   hPtJet->Sumw2();
   
   TH3F* hPtJetWithD=new TH3F("hPtJetWithD","D-Jet Pt distribution; p_{T} (GeV/c);delta mass (GeV/c^{2})",nbinsptjet,ptjetlims[0],ptjetlims[1],nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
   hPtJetWithD->Sumw2();
   //for the MC this histogram is filled with the real background
   TH3F* hPtJetWithDsb=new TH3F("hPtJetWithDsb","D(background)-Jet Pt distribution; p_{T} (GeV/c);delta mass (GeV/c^{2});p_{T}^{D} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1],nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
   hPtJetWithDsb->Sumw2();
   
   TH1F* hdeltaRJetTracks=new TH1F("hdeltaRJetTracks","Delta R of tracks in the jets",200, 0.,10.);
   hdeltaRJetTracks->Sumw2();
   
   TH1F* hNtrArr= new TH1F("hNtrArr", "Number of tracks in the array of jets; number of tracks",500,0,1000);
   hNtrArr->Sumw2();
   
   TH1F *hNJetPerEv=new TH1F("hNJetPerEv","Number of jets used per event; number of jets/ev",10,-0.5,9.5);
   hNJetPerEv->Sumw2();
   
   TH1F *hNJetPerEvNoD=new TH1F("hNJetPerEvNoD","Number of jets per event with no D; number of jets/ev with no D",10,-0.5,9.5);
   hNJetPerEvNoD->Sumw2();
   
   TH1F *hPtJetPerEvNoD=new TH1F("hPtJetPerEvNoD","pt distribution of jets per event with no D; p_{T}^{jet} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
   hPtJetPerEvNoD->Sumw2();
    
   fmyOutput->Add(hEjetTrks);
   fmyOutput->Add(hPhiJetTrks);
   fmyOutput->Add(hEtaJetTrks);
   fmyOutput->Add(hPtJetTrks);
   fmyOutput->Add(hEjet);
   fmyOutput->Add(hPhiJet);
   fmyOutput->Add(hEtaJet);
   fmyOutput->Add(hPtJet);
   fmyOutput->Add(hPtJetWithD);
   fmyOutput->Add(hPtJetWithDsb);
   fmyOutput->Add(hdeltaRJetTracks);
   fmyOutput->Add(hNtrArr);
   fmyOutput->Add(hNJetPerEv);
   fmyOutput->Add(hNJetPerEvNoD);
   fmyOutput->Add(hPtJetPerEvNoD);
   
   TH1F* hDeltaRD=new TH1F("hDeltaRD","#Delta R distribution of D candidates selected;#Delta R",200, 0.,10.);
   hDeltaRD->Sumw2();
   fmyOutput->Add(hDeltaRD);
   
   //background (side bands for the Dstar and like sign for D0)
   fJetRadius=GetJetContainer(0)->GetJetRadius();
   TH2F* hInvMassptD = new TH2F("hInvMassptD",Form("D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
   hInvMassptD->SetStats(kTRUE);
   hInvMassptD->GetXaxis()->SetTitle("mass (GeV)");
   hInvMassptD->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
   hInvMassptD->Sumw2();
   
   fmyOutput->Add(hInvMassptD);
   
   if(fUseMCInfo){
      TH2F* hInvMassptDbg = new TH2F("hInvMassptDbg",Form("Background D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      hInvMassptDbg->GetXaxis()->SetTitle(Form("%s (GeV)",(fCandidateType==kDstartoKpipi) ? "M(Kpipi)" : "M(Kpi)"));
      hInvMassptDbg->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      hInvMassptDbg->Sumw2();
      fmyOutput->Add(hInvMassptDbg);
      
   }
   Double_t pi=TMath::Pi(), philims2[2];
   philims2[0]=-pi*1./2.;
   philims2[1]=pi*3./2.;
   const Int_t nAxis=6;   
   const Int_t nbinsSparse[nAxis]={nbinsz,nbinsphi,nbinsptjet,nbinsptD,nbinsmass,2};
   const Double_t minSparse[nAxis]={zlims[0],philims2[0],ptjetlims[0],ptDlims[0],fMinMass,-0.5};
   const Double_t maxSparse[nAxis]={zlims[1],philims2[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5};
   THnSparseF *hsDphiz=new THnSparseF("hsDphiz","Z and #Delta#phi vs p_{T}^{jet}, p_{T}^{D}, and mass", nAxis, nbinsSparse, minSparse, maxSparse);
   hsDphiz->Sumw2();
   
   fmyOutput->Add(hsDphiz);

   PostData(1, fmyOutput);
   
   return kTRUE; 
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourFilterAndJetCorrelations::DefineHistoForFiltering()
{
   //
   // AliAnalysisTaskSEDmesonsFilterCJ::DefineHistoForAnalysis
   //
   
   // Statistics 
   TH1I* hstat = new TH1I("hstatF","Statistics",5,-0.5,4.5);
   hstat->GetXaxis()->SetBinLabel(1, "N ev anal");
   hstat->GetXaxis()->SetBinLabel(2, "N ev sel");
   if(fUseReco) hstat->GetXaxis()->SetBinLabel(3, "N cand");
   else hstat->GetXaxis()->SetBinLabel(3, "N Gen D");
   hstat->GetXaxis()->SetBinLabel(4, "N cand sel cuts");
   if (fCandidateType==kDstartoKpipi) hstat->GetXaxis()->SetBinLabel(5, "N side band cand");  
   hstat->SetNdivisions(1);
   fmyOutputF->Add(hstat);
   
   TH1F* hnCandEv=new TH1F("hnCandEv", "Number of candidates per event (after cuts);# cand/ev", 100, 0.,100.);
   fmyOutputF->Add(hnCandEv);
   
   // Invariant mass related histograms
   const Int_t nbinsmass = 200;
   TH2F *hInvMass = new TH2F("hInvMassptDF", "D invariant mass distribution", nbinsmass, fMinMass, fMaxMass, 100, 0., 50.);
   hInvMass->SetStats(kTRUE);
   hInvMass->GetXaxis()->SetTitle("mass (GeV/c)");
   hInvMass->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   fmyOutputF->Add(hInvMass);
   
   if (fCandidateType==kDstartoKpipi) {
      TH1F* hnSBCandEv=new TH1F("hnSBCandEv", "Number of side bands candidates per event (after cuts);# cand/ev", 100, 0.,100.);
      fmyOutputF->Add(hnSBCandEv);
      
      TH1F* hPtPion = new TH1F("hPtPion", "Primary pions candidates pt", 500, 0., 10.);
      hPtPion->SetStats(kTRUE);
      hPtPion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtPion->GetYaxis()->SetTitle("entries");
      fmyOutputF->Add(hPtPion);
   }
   PostData(3, fmyOutputF);
   return kTRUE; 
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FillHistogramsRecoJetCorr(AliVParticle* candidate, AliEmcalJet *jet){
   
   Double_t ptD=candidate->Pt();
   Double_t ptjet=jet->Pt();
   Double_t deltaR=DeltaR(candidate,jet);
   Double_t deltaphi = jet->Phi()-candidate->Phi();
   if(deltaphi<=-(TMath::Pi())/2) deltaphi = deltaphi+2*(TMath::Pi());
   if(deltaphi>(3*(TMath::Pi()))/2) deltaphi = deltaphi-2*(TMath::Pi());
   Double_t z=Z(candidate,jet);
   
   TH1F* hDeltaRD=(TH1F*)fmyOutput->FindObject("hDeltaRD");
   hDeltaRD->Fill(deltaR); 
   if(fUseReco){
      if(fCandidateType==kD0toKpi) {
      	 AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)candidate;
      	 FillHistogramsD0JetCorr(dzero,deltaphi,z,ptD,ptjet,deltaR, AODEvent());
      	 
      }
      
      if(fCandidateType==kDstartoKpipi) {
      	 AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)candidate;
      	 FillHistogramsDstarJetCorr(dstar,deltaphi,z,ptD,ptjet,deltaR);
      	 
      }
   } else{ //generated level
      //AliInfo("Non reco");
      FillHistogramsMCGenDJetCorr(deltaphi,z,ptD,ptjet,deltaR);
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t ptj,Double_t deltaR, AliAODEvent* aodEvent){


   Double_t masses[2]={0.,0.};
   Int_t pdgdaughtersD0[2]={211,321};//pi,K 
   Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
   
   masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
   masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
   
   TH3F* hPtJetWithD=(TH3F*)fmyOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fmyOutput->FindObject("hsDphiz");
   Double_t point[5]={z,dPhi,ptj,ptD,masses[0]};
   
   Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);
   if(isselected==1 || isselected==3) {
      
      if(deltaR<fJetRadius) hPtJetWithD->Fill(ptj,masses[0],ptD);
      
      FillMassHistograms(masses[0], ptD, deltaR);
      hsDphiz->Fill(point,1.);
   }
   if(isselected>=2) {
      if(deltaR<fJetRadius) hPtJetWithD->Fill(ptj,masses[1],ptD);
      
      FillMassHistograms(masses[1], ptD, deltaR);
      point[4]=masses[1];
      hsDphiz->Fill(point,1.);
   }
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi, Double_t z, Double_t ptD, Double_t ptj, Double_t deltaR){
  //dPhi and z not used at the moment,but will be (re)added

   AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
   Double_t deltamass= dstar->DeltaInvMass(); 
   //Double_t massD0= dstar->InvMassD0();
   
   
   TH1F* hPtPion=(TH1F*)fmyOutput->FindObject("hPtPion");
   hPtPion->Fill(softpi->Pt());
   
   TH3F* hPtJetWithD=(TH3F*)fmyOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fmyOutput->FindObject("hsDphiz");
   Int_t isSB=IsDzeroSideBand(dstar);
   
   Double_t point[6]={z,dPhi,ptj,ptD,deltamass,isSB};

   if(deltaR<fJetRadius) hPtJetWithD->Fill(ptj,deltamass,ptD);
   
   FillMassHistograms(deltamass, ptD, deltaR);
   hsDphiz->Fill(point,1.);
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FillHistogramsMCGenDJetCorr(Double_t dPhi,Double_t z,Double_t ptD,Double_t ptjet,Double_t deltaR){
   
   Double_t pdgmass=0;
   TH3F* hPtJetWithD=(TH3F*)fmyOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fmyOutput->FindObject("hsDphiz");
   Double_t point[6]={z,dPhi,ptjet,ptD,pdgmass,0};

   if(fCandidateType==kD0toKpi) pdgmass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
   if(fCandidateType==kDstartoKpipi) pdgmass=TDatabasePDG::Instance()->GetParticle(413)->Mass()-TDatabasePDG::Instance()->GetParticle(421)->Mass();
   point[4]=pdgmass;

   if(deltaR<fJetRadius) {
     hPtJetWithD->Fill(ptjet,pdgmass,ptD); // candidates within a jet
     hsDphiz->Fill(point,1.);
   }
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FillMassHistograms(Double_t mass,Double_t ptD, Double_t deltaR){
   
   if(deltaR<fJetRadius) {
      TH2F* hInvMassptD=(TH2F*)fmyOutput->FindObject("hInvMassptD");
      hInvMassptD->Fill(mass,ptD);
   }
}

//________________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::FlagFlavour(AliVParticle *charm, AliEmcalJet *jet){
   Double_t deltaR=DeltaR(charm, jet);
   AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
   if (fCandidateType==kD0toKpi) tag=AliEmcalJet::kD0;
   if (deltaR<fJetRadius) jet->AddFlavourTag(tag);
   
   return;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::SideBandBackground(AliAODRecoCascadeHF *candDstar, AliEmcalJet *jet){
   
   //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
   // (expected detector resolution) on the left and right frm the D0 mass. Each band
   //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
   TH2F* hDiffSideBand=(TH2F*)fmyOutput->FindObject("hDiffSideBand");
   TH3F* hPtJetWithDsb=(TH3F*)fmyOutput->FindObject("hPtJetWithDsb");
   
   Double_t deltaM=candDstar->DeltaInvMass(); 
   //Printf("Inv mass = %f between %f and %f or %f and %f?",invM, sixSigmal,fourSigmal,fourSigmar,sixSigmar);
   Double_t ptD=candDstar->Pt();
   Double_t ptjet=jet->Pt();
   Double_t dPhi=jet->Phi()-candDstar->Phi();
   Double_t deltaR=DeltaR(candDstar,jet);
   if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
   if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();
   
   Int_t isSideBand=IsDzeroSideBand(candDstar);
   //fill the background histogram with the side bands or when looking at MC with the real background
   if(isSideBand==1){
      hDiffSideBand->Fill(deltaM,ptD); // M(Kpipi)-M(Kpi) side band background    
      //hdeltaPhiDjaB->Fill(deltaM,ptD,dPhi);
      
      if(deltaR<fJetRadius){  // evaluate in the near side	
      	 //hzptDB->Fill(Z(candDstar,jet),deltaM,ptD);
      	 hPtJetWithDsb->Fill(ptjet,deltaM,ptD);
      }
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourFilterAndJetCorrelations::MCBackground(AliAODRecoDecayHF *candbg, AliEmcalJet *jet){
   
   //need mass, deltaR, pt jet, ptD
   //two cases: D0 and Dstar
   
   Int_t isselected=fCuts->IsSelected(candbg,AliRDHFCuts::kAll, AODEvent());
   if(!isselected) return;
   
   Double_t ptD=candbg->Pt();
   Double_t ptjet=jet->Pt();
   Double_t deltaR=DeltaR(candbg,jet);
   
   TH2F* hInvMassptDbg=(TH2F*)fmyOutput->FindObject("hInvMassptDbg");
   TH3F* hPtJetWithDsb=(TH3F*)fmyOutput->FindObject("hPtJetWithDsb");
   
   
   if(fCandidateType==kDstartoKpipi){
      AliAODRecoCascadeHF* dstarbg = (AliAODRecoCascadeHF*)candbg;
      Double_t deltaM=dstarbg->DeltaInvMass();
      hInvMassptDbg->Fill(deltaM,ptD);
      if(deltaR<fJetRadius) hPtJetWithDsb->Fill(ptjet,deltaM,ptD);
   }
   
   if(fCandidateType==kD0toKpi){
      Double_t masses[2]={0.,0.};
      Int_t pdgdaughtersD0[2]={211,321};//pi,K 
      Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      
      masses[0]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
      masses[1]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
      
      
      if(isselected==1 || isselected==3) {
      	 if(deltaR<fJetRadius) hPtJetWithDsb->Fill(ptjet,masses[0],ptD);
      	 hInvMassptDbg->Fill(masses[0],ptD);
      }
      if(isselected>=2) {
      	 if(deltaR<fJetRadius) hPtJetWithDsb->Fill(ptjet,masses[1],ptD);
      	 hInvMassptDbg->Fill(masses[1],ptD);
      }
      
      
   }
}

//_______________________________________________________________________________

Float_t AliAnalysisTaskFlavourFilterAndJetCorrelations::DeltaR(AliVParticle *p1, AliVParticle *p2) const {
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

//_______________________________________________________________________________

Int_t AliAnalysisTaskFlavourFilterAndJetCorrelations::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){
   
   Double_t ptD=candDstar->Pt();
   Int_t bin = fCuts->PtBin(ptD);
   if (bin < 0){
      // /PWGHF/vertexingHF/AliRDHFCuts::PtBin(Double_t) const may return values below zero depending on config.
      bin = 9999; // void result code for coverity (bin later in the code non-zero) - will coverity pick up on this?      
      return -1;  // out of bounds
   }
   
   Double_t invM=candDstar->InvMassD0();
   Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();

   Float_t fourSigmal= mPDGD0-4.*fSigmaD0[bin] , sixSigmal= mPDGD0-8.*fSigmaD0[bin];
   Float_t fourSigmar= mPDGD0+4.*fSigmaD0[bin] , sixSigmar= mPDGD0+8.*fSigmaD0[bin];
   
   if((invM>=sixSigmal && invM<fourSigmal) || (invM>fourSigmar && invM<=sixSigmar)) return 1;
   else return 0;   

}
