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
#include <TSystem.h>
#include <TObjectTable.h>

#include "AliAnalysisTaskFlavourJetCorrelations.h"
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
#include "AliRhoParameter.h"

ClassImp(AliAnalysisTaskFlavourJetCorrelations)


//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations() :
AliAnalysisTaskEmcalJet("",kTRUE),
fUseMCInfo(kTRUE), 
fUseReco(kTRUE),
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fCandArrName(0),
fLeadingJetOnly(kFALSE),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0),
fJetOnlyMode(0),
fPmissing(),
fPtJet(0),
fIsDInJet(0),
fTypeDInJet(2),
fTrackArr(0),
fSwitchOnSB(0),
fSwitchOnPhiAxis(0),
fSwitchOnOutOfConeAxis(0),
fSwitchOnSparses(1),
fNAxesBigSparse(9)
{
   //
   // Default ctor
   //
   //SetMakeGeneralHistograms(kTRUE);
   
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype, Bool_t jetOnly) :
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseMCInfo(kTRUE),
fUseReco(kTRUE),  
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fCandArrName(0),
fLeadingJetOnly(kFALSE),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0),
fJetOnlyMode(jetOnly),
fPmissing(),
fPtJet(0),
fIsDInJet(0),
fTypeDInJet(2),
fTrackArr(0),
fSwitchOnSB(0),
fSwitchOnPhiAxis(0),
fSwitchOnOutOfConeAxis(0),
fSwitchOnSparses(1),
fNAxesBigSparse(9)
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //
   
   Info("AliAnalysisTaskFlavourJetCorrelations","Calling Constructor");
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
   
   if(fJetOnlyMode){
      DefineOutput(1,TList::Class()); //histos with jet info
      DefineOutput(2,AliRDHFCuts::Class()); // my cuts
   }
   else{
      DefineInput(1, TClonesArray::Class());
      DefineInput(2, TClonesArray::Class());
      
      DefineOutput(1,TList::Class()); // histos
      DefineOutput(2,AliRDHFCuts::Class()); // my cuts
   }
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::~AliAnalysisTaskFlavourJetCorrelations() {
   //
   // destructor
   //
   
   Info("~AliAnalysisTaskFlavourJetCorrelations","Calling Destructor");  
   
   delete fCuts;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::Init(){
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

void AliAnalysisTaskFlavourJetCorrelations::UserCreateOutputObjects() { 
   // output 
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   AliAnalysisTaskEmcal::UserCreateOutputObjects();
   // define histograms
   // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()
   DefineHistoForAnalysis();
   PostData(1,fOutput);
   
   return;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::Run()
{
   // user exec from AliAnalysisTaskEmcal is used
    
   // Load the event
   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
   
   TClonesArray *arrayDStartoD0pi=0;
   
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
      	 return kFALSE;
      }
   }
   
   //retrieve jets
   fTrackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoTracks"));
   //clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClustersCorr"));
   //jetArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrName));
   fJetRadius=GetJetContainer(0)->GetJetRadius();
   
   if(!fTrackArr){
      AliInfo(Form("Could not find the track array\n"));
      return kFALSE;
   }
   
    
   fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
   if (!fCandidateArray) return kFALSE;
   if (fCandidateType==1 && fSwitchOnSB) {
      fSideBandArray = dynamic_cast<TClonesArray*>(GetInputData(2));
      if (!fSideBandArray) return kFALSE;
   }
   //Printf("ncandidates found %d",fCandidateArray->GetEntriesFast());
   
   //Histograms
   TH1I* hstat=(TH1I*)fOutput->FindObject("hstat");
   TH1F* hPtJetTrks=(TH1F*)fOutput->FindObject("hPtJetTrks");
   TH1F* hPhiJetTrks=(TH1F*)fOutput->FindObject("hPhiJetTrks");
   TH1F* hEtaJetTrks=(TH1F*)fOutput->FindObject("hEtaJetTrks");
   TH1F* hEjetTrks=(TH1F*)fOutput->FindObject("hEjetTrks");
   TH1F* hPtJet=(TH1F*)fOutput->FindObject("hPtJet");
   TH1F* hPhiJet=(TH1F*)fOutput->FindObject("hPhiJet");
   TH1F* hEtaJet=(TH1F*)fOutput->FindObject("hEtaJet");
   TH1F* hEjet=(TH1F*)fOutput->FindObject("hEjet");
   TH1F* hNtrArr=(TH1F*)fOutput->FindObject("hNtrArr");
   TH1F* hNJetPerEv=(TH1F*)fOutput->FindObject("hNJetPerEv");
   TH1F* hdeltaRJetTracks=(TH1F*)fOutput->FindObject("hdeltaRJetTracks");
   TH1F* hNDPerEvNoJet=(TH1F*)fOutput->FindObject("hNDPerEvNoJet");
   TH1F* hptDPerEvNoJet=(TH1F*)fOutput->FindObject("hptDPerEvNoJet");
   TH1F* hNJetPerEvNoD=(TH1F*)fOutput->FindObject("hNJetPerEvNoD");
   TH1F* hPtJetPerEvNoD=(TH1F*)fOutput->FindObject("hPtJetPerEvNoD");
   THnSparseF* hnspDstandalone=(THnSparseF*)fOutput->FindObject("hsDstandalone");
   THnSparseF* hsJet=(THnSparseF*)fOutput->FindObject("hsJet");
   
   TH1F* hztracksinjet=(TH1F*)fOutput->FindObject("hztracksinjet");
   TH1F* hDiffPtTrPtJzNok=(TH1F*)fOutput->FindObject("hDiffPtTrPtJzNok");
   TH1F* hDiffPtTrPtJzok=(TH1F*)fOutput->FindObject("hDiffPtTrPtJzok");
   TH1F* hDiffPzTrPtJzok=(TH1F*)fOutput->FindObject("hDiffPzTrPtJzok");
   TH1F* hDiffPzTrPtJzNok=(TH1F*)fOutput->FindObject("hDiffPzTrPtJzNok");
   TH1I* hNtrkjzNok=(TH1I*)fOutput->FindObject("hNtrkjzNok");
   TH1F* hztracksinjetT=(TH1F*)fOutput->FindObject("hztracksinjetT");
   
   
   hstat->Fill(0);
   
   // fix for temporary bug in ESDfilter 
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return kFALSE;
   
   //Event selection
   Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if(!iseventselected) return kFALSE;
   
   hstat->Fill(1);

   //retrieve charm candidates selected
   Int_t candidates = 0;
   Int_t njets=GetJetContainer()->GetNJets();
   
   if(!fJetOnlyMode) {
      candidates = fCandidateArray->GetEntriesFast();
  
   //trigger on jets  
   if(njets == 0) {
      hstat->Fill(6, candidates);
      hNDPerEvNoJet->Fill(candidates);
      for(Int_t iD=0;iD<candidates;iD++){
      	 AliVParticle* cand=(AliVParticle*)fCandidateArray->At(iD);
      	 if(!cand) continue;
      	 hptDPerEvNoJet->Fill(cand->Pt());
      
      }
      return kFALSE;
      
   }
   
   //loop on candidates standalone (checking the candidates are there and their phi-eta distributions)
   
   for(Int_t ic = 0; ic < candidates; ic++) {
      
      // D* candidates
      AliAODRecoDecayHF* charm=0x0;
      AliAODRecoCascadeHF* dstar=0x0;
      
      
      charm=(AliAODRecoDecayHF*)fCandidateArray->At(ic);
      if(!charm) continue;
      
      if(fCandidateType==kDstartoKpipi){ 
      	 dstar=(AliAODRecoCascadeHF*)fCandidateArray->At(ic);
      	 if(!dstar) continue;
      }
      
      hstat->Fill(2);
      
      Double_t candsparse[4]={charm->Eta(), charm->Phi(), charm->Pt(), 0};
      
      if(fCandidateType==kDstartoKpipi) {
      	 if(fUseReco){
      	    Double_t deltamass= dstar->DeltaInvMass();
      	    candsparse[3]=deltamass;
      	 } else candsparse[3] = 0.145; //for generated
      	 if(fSwitchOnSparses) hnspDstandalone->Fill(candsparse);
      }
      if(fCandidateType==kD0toKpi){
      	 if(fUseReco){
      	    AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)charm;
      	    Int_t isselected=fCuts->IsSelected(dzero,AliRDHFCuts::kAll,aodEvent);
      	    //not a further selection,just to get the value of isselected for the D0
      	    Double_t masses[2];
      	    Int_t pdgdaughtersD0[2]={211,321};//pi,K 
      	    Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      	    
      	    masses[0]=dzero->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
      	    masses[1]=dzero->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
      	    if(isselected==1 || isselected==3) {
      	       candsparse[3]=masses[0];
      	       if(fSwitchOnSparses) hnspDstandalone->Fill(candsparse);
      	    }
      	    if(isselected>=2){
      	       candsparse[3]=masses[1];
      	       if(fSwitchOnSparses) hnspDstandalone->Fill(candsparse);
      	       
      	    }
      	 } else { //generated
      	    Int_t pdg=((AliAODMCParticle*)charm)->GetPdgCode();
      	    candsparse[3]=TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
      	    if(fSwitchOnSparses) hnspDstandalone->Fill(candsparse);
      	 }
      }
   }
   }
    
    //Background Subtraction for jets
    AliRhoParameter *rho = 0;
    Double_t rhoval = 0;
    TString sname("Rho");
    if (!sname.IsNull()) {
        rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(sname));
        if(rho) rhoval = rho->GetVal();
    }

   
   // we start with jets
   Double_t ejet   = 0;
   Double_t phiJet = 0;
   Double_t etaJet = 0;
   Double_t ptjet = 0;
   Double_t leadingJet =0;
   Double_t pointJ[6];
   
   Int_t ntrarr=fTrackArr->GetEntriesFast();
   hNtrArr->Fill(ntrarr);
   
   for(Int_t i=0;i<ntrarr;i++){
      AliVTrack *jtrack=static_cast<AliVTrack*>(fTrackArr->At(i));
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
      	 ptjet   = jetL->Pt() - jetL->Area()*rhoval; //background subtraction
      	 if(ptjet>leadingJet ) leadingJet = ptjet;
      	 
      }
   }
   
   Int_t cntjet=0;
   //loop over jets and consider only the leading jet in the event
   for (Int_t iJets = 0; iJets<njets; iJets++) {
      fPmissing[0]=0;
      fPmissing[1]=0;
      fPmissing[2]=0;
      
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
      fPtJet = jet->Pt() - jet->Area()*rhoval; //background subtraction
      Double_t origPtJet=fPtJet;
      
      pointJ[0] = phiJet;
      pointJ[1] = etaJet;
      pointJ[2] = ptjet - jet->Area()*rhoval; //background subtraction
      pointJ[3] = ejet;
      pointJ[4] = jet->GetNumberOfConstituents();
      pointJ[5] = jet->Area();
      
      // choose the leading jet
      if(fLeadingJetOnly && (ejet<leadingJet)) continue;
      //used for normalization
      hstat->Fill(3);
      cntjet++;
      // fill energy, eta and phi of the jet
      hEjet   ->Fill(ejet);
      hPhiJet ->Fill(phiJet);
      hEtaJet ->Fill(etaJet);
      hPtJet  ->Fill(fPtJet);
      if(fJetOnlyMode && fSwitchOnSparses) hsJet->Fill(pointJ,1);
      //loop on jet particles
      Int_t ntrjet=  jet->GetNumberOfTracks(); 
      Double_t sumPtTracks=0,sumPzTracks=0;
      Int_t zg1jtrk=0;
      for(Int_t itrk=0;itrk<ntrjet;itrk++){
      	 
      	 AliPicoTrack* jetTrk=(AliPicoTrack*)jet->TrackAt(itrk,fTrackArr);     
      	 hdeltaRJetTracks->Fill(DeltaR(jet,jetTrk));
      	 Double_t ztrackj=Z(jetTrk,jet);
     	 hztracksinjet->Fill(ztrackj);
     	 sumPtTracks+=jetTrk->Pt(); 
     	 sumPzTracks+=jetTrk->Pz(); 
     	 if(ztrackj>1){
     	    zg1jtrk++;
     	 }
     	 
     	 Double_t ztrackjTr=Z(jetTrk,jet,kTRUE);
     	 hztracksinjetT->Fill(ztrackjTr);
     	 
      }//end loop on jet tracks
      
      hNtrkjzNok->Fill(zg1jtrk);
      Double_t diffpt=origPtJet-sumPtTracks;
      Double_t diffpz=jet->Pz()-sumPzTracks;
      if(zg1jtrk>0){
      	 hDiffPtTrPtJzNok->Fill(diffpt);
      	 hDiffPzTrPtJzNok->Fill(diffpz);
      
      }else{
      	 hDiffPtTrPtJzok->Fill(diffpt);
      	 hDiffPzTrPtJzok->Fill(diffpz);
      }
      
      if(candidates==0){
      	 hstat->Fill(7);
      	 hPtJetPerEvNoD->Fill(fPtJet);
      }
      if(!fJetOnlyMode) {
      	 //Printf("N candidates %d ", candidates);
      	 for(Int_t ic = 0; ic < candidates; ic++) {
      	    
      	    // D* candidates
      	    AliVParticle* charm=0x0;
      	    charm=(AliVParticle*)fCandidateArray->At(ic);
      	    if(!charm) continue;
      	    AliAODRecoDecayHF *charmdecay=(AliAODRecoDecayHF*) charm;
      	    fIsDInJet=IsDInJet(jet, charmdecay, kTRUE);
      	    if (fIsDInJet) FlagFlavour(jet);
      	    if (jet->TestFlavourTag(AliEmcalJet::kDStar)) hstat->Fill(4);
      	    
      	    //Note: the z component of the jet momentum comes from the eta-phi direction of the jet particles, it is not calculated from the z component of the tracks since, as default, the scheme used for jet reco is the pt-scheme which sums the scalar component, not the vectors. Addind the D daughter momentum component by componet as done here is not 100% correct, but the difference is small, for fairly collimated particles.

      	    Double_t pjet[3];
      	    jet->PxPyPz(pjet);
             //background subtraction
             pjet[0] = jet->Px() - jet->Area()*(rhoval*TMath::Cos(jet->AreaPhi()));
             pjet[1] = jet->Py() - jet->Area()*(rhoval*TMath::Sin(jet->AreaPhi()));
             pjet[2] = jet->Pz() - jet->Area()*(rhoval*TMath::SinH(jet->AreaEta()));
      	    RecalculateMomentum(pjet,fPmissing);      	          	    
      	    fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]);
      	    
      	    
      	    //debugging histograms
      	    Double_t pmissing=TMath::Sqrt(fPmissing[0]*fPmissing[0]+fPmissing[1]*fPmissing[1]+fPmissing[2]*fPmissing[2]);
      	    for(Int_t k=0;k<3;k++) ((TH1F*)fOutput->FindObject(Form("hMissP%d",k)))->Fill(fPmissing[k]);
      	    
      	    ((TH1F*)fOutput->FindObject("hmissingp"))->Fill(pmissing);
      	    Double_t ptdiff=fPtJet-origPtJet;
      	    ((TH1F*)fOutput->FindObject("hDeltaPtJet"))->Fill(ptdiff);
      	    ((TH1F*)fOutput->FindObject("hRelDeltaPtJet"))->Fill(ptdiff/origPtJet);
      	    
      	    FillHistogramsRecoJetCorr(charm, jet, aodEvent);
      	    
      	 }//end loop on candidates
      	 
      	 //retrieve side band background candidates for Dstar
      	 Int_t nsbcand = 0; if(fSideBandArray) nsbcand=fSideBandArray->GetEntriesFast();
      	 
      	 for(Int_t ib=0;ib<nsbcand;ib++){
      	    if(fCandidateType==kDstartoKpipi && !fUseMCInfo){
      	       AliAODRecoCascadeHF *sbcand=(AliAODRecoCascadeHF*)fSideBandArray->At(ib);
      	       if(!sbcand) continue;
      	        
      	       fIsDInJet=IsDInJet(jet, sbcand,kFALSE);
      	       Double_t pjet[3];
      	       jet->PxPyPz(pjet);
                //background subtraction
                pjet[0] = jet->Px() - jet->Area()*(rhoval*TMath::Cos(jet->AreaPhi()));
                pjet[1] = jet->Py() - jet->Area()*(rhoval*TMath::Sin(jet->AreaPhi()));
                pjet[2] = jet->Pz() - jet->Area()*(rhoval*TMath::SinH(jet->AreaEta()));
      	       RecalculateMomentum(pjet,fPmissing);      	          	    
      	       fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]);
      	       
     	       SideBandBackground(sbcand,jet);
     	       
      	    }
      	    if(fUseMCInfo){
      	       AliAODRecoDecayHF* charmbg = 0x0;
      	       charmbg=(AliAODRecoDecayHF*)fSideBandArray->At(ib);
      	       if(!charmbg) continue;
      	       fIsDInJet=IsDInJet(jet, charmbg,kFALSE);
      	       if (fIsDInJet) FlagFlavour(jet); //this are backgroud HF jets, but flagged as signal at the moment. Can use the bkg flavour flag in the future. This info is not stored now a part in the jet
      	       Double_t pjet[3];
      	       jet->PxPyPz(pjet);
      	       //background subtraction
      	       pjet[0] = jet->Px() - jet->Area()*(rhoval*TMath::Cos(jet->AreaPhi()));
      	       pjet[1] = jet->Py() - jet->Area()*(rhoval*TMath::Sin(jet->AreaPhi()));
      	       pjet[2] = jet->Pz() - jet->Area()*(rhoval*TMath::SinH(jet->AreaEta()));
      	       RecalculateMomentum(pjet,fPmissing);      	          	    
      	       fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]);
      	       
      	       MCBackground(charmbg,jet);
      	    }
      	 }
      }
   } // end of jet loop
   
   hNJetPerEv->Fill(cntjet);
   if(candidates==0) hNJetPerEvNoD->Fill(cntjet);
   PostData(1,fOutput);
   return kTRUE;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::Terminate(Option_t*)
{    
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
   Info("Terminate"," terminate");
   AliAnalysisTaskSE::Terminate();
   
   fOutput = dynamic_cast<TList*> (GetOutputData(1));
   if (!fOutput) {     
      printf("ERROR: fOutput not available\n");
      return;
   }
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t range, Int_t pdg){
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

void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t lowlimit, Double_t uplimit){
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

Bool_t AliAnalysisTaskFlavourJetCorrelations::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________

Double_t AliAnalysisTaskFlavourJetCorrelations::Z(AliVParticle* part,AliEmcalJet* jet, Bool_t transverse) const{
   if(!part || !jet){
      printf("Particle or jet do not exist!\n");
      return -99;
   }
   Double_t p[3],pj[3];
   Bool_t okpp=part->PxPyPz(p);
   Bool_t okpj=jet->PxPyPz(pj);
    
    //Background Subtraction
    AliRhoParameter *rho = 0;
    Double_t rhoval = 0;
    TString sname("Rho");
    if (!sname.IsNull()) {
        rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(sname));
        if(rho){
            rhoval = rho->GetVal();
            pj[0] = jet->Px() - jet->Area()*(rhoval*TMath::Cos(jet->AreaPhi()));
            pj[1] = jet->Py() - jet->Area()*(rhoval*TMath::Sin(jet->AreaPhi()));
            pj[2] = jet->Pz() - jet->Area()*(rhoval*TMath::SinH(jet->AreaEta()));
        }
    }

    
    
   if(!okpp || !okpj){
      printf("Problems getting momenta\n");
      return -999;
   }
   
   RecalculateMomentum(pj, fPmissing);
   if(transverse) return ZT(p,pj);
   else
   return Z(p,pj);
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelations::Z(Double_t* p, Double_t *pj) const{
   
   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1]+pj[2]*pj[2];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet2);
   return z;
}


//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelations::ZT(Double_t* p, Double_t *pj) const{
   
   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1])/(pjet2);
   return z;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::RecalculateMomentum(Double_t* pj, const Double_t* pmissing) const {

   pj[0]+=pmissing[0];
   pj[1]+=pmissing[1];
   pj[2]+=pmissing[2];

}

//_______________________________________________________________________________

Bool_t  AliAnalysisTaskFlavourJetCorrelations::DefineHistoForAnalysis(){
   
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
   fOutput->Add(hstat);
   
   const Int_t nbinsmass=200;
   const Int_t nbinsptjet=500;
   const Int_t nbinsptD=100;
   const Int_t nbinsz=100;
   const Int_t nbinsphi=200;
   const Int_t nbinseta=100;
   
   //binning for THnSparse
   const Int_t nbinsSpsmass=50;
   const Int_t nbinsSpsptjet=100;
   const Int_t nbinsSpsptD=50;
   const Int_t nbinsSpsz=100;
   const Int_t nbinsSpsphi=100;
   const Int_t nbinsSpseta=60;
   const Int_t nbinsSpsContrib=100;
   const Int_t nbinsSpsA=100;
    
   const Float_t ptjetlims[2]={0.,200.};
   const Float_t ptDlims[2]={0.,50.};
   const Float_t zlims[2]={0.,1.2};
   const Float_t philims[2]={0.,6.3};
   const Float_t etalims[2]={-1.5,1.5};
   const Int_t   nContriblims[2]={0,100};
   const Float_t arealims[2]={0.,2};
   
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
   
   TH1F* hdeltaRJetTracks=new TH1F("hdeltaRJetTracks","Delta R of tracks in the jets",200, 0.,10.);
   hdeltaRJetTracks->Sumw2();
   
   TH1F* hNtrArr= new TH1F("hNtrArr", "Number of tracks in the array of jets; number of tracks",500,0,1000);
   hNtrArr->Sumw2();
   
   TH1F *hNJetPerEv=new TH1F("hNJetPerEv","Number of jets used per event; number of jets/ev",10,-0.5,9.5);
   hNJetPerEv->Sumw2();
   
   
   fOutput->Add(hEjetTrks);
   fOutput->Add(hPhiJetTrks);
   fOutput->Add(hEtaJetTrks);
   fOutput->Add(hPtJetTrks);
   fOutput->Add(hEjet);
   fOutput->Add(hPhiJet);
   fOutput->Add(hEtaJet);
   fOutput->Add(hPtJet);
   fOutput->Add(hdeltaRJetTracks);
   fOutput->Add(hNtrArr);
   fOutput->Add(hNJetPerEv);
   
   if(fJetOnlyMode && fSwitchOnSparses){
      //thnsparse for jets
      const Int_t nAxis=6;   
      const Int_t nbinsSparse[nAxis]={nbinsSpsphi,nbinsSpseta, nbinsSpsptjet, nbinsSpsptjet,nbinsSpsContrib, nbinsSpsA};
      const Double_t minSparse[nAxis]={philims[0],etalims[0],ptjetlims[0],ptjetlims[0],static_cast<Double_t>(nContriblims[0]),arealims[0]};
      const Double_t 
	maxSparse[nAxis]={philims[1],etalims[1],ptjetlims[1],ptjetlims[1],static_cast<Double_t>(nContriblims[1]),arealims[1]};
      THnSparseF *hsJet=new THnSparseF("hsJet","#phi, #eta, p_{T}^{jet}, E^{jet}, N contrib, Area", nAxis, nbinsSparse, minSparse, maxSparse);
      hsJet->Sumw2();
      
      fOutput->Add(hsJet);
   
   }

   if(!fJetOnlyMode){
      
      //debugging histograms
      TH1I* hControlDInJ=new TH1I("hControlDInJ","Checks D in Jet",8, -0.5,7.5);
      hControlDInJ->GetXaxis()->SetBinLabel(1,"DR In,1 daugh out");
      hControlDInJ->GetXaxis()->SetBinLabel(2,"DR In,2 daugh out");
      hControlDInJ->GetXaxis()->SetBinLabel(3,"DR In,3 daugh out");
      hControlDInJ->GetXaxis()->SetBinLabel(4,"Tot tracks non matched");
      hControlDInJ->GetXaxis()->SetBinLabel(5,"D0 daug missing");
      hControlDInJ->GetXaxis()->SetBinLabel(6,"soft pi missing");
      hControlDInJ->GetXaxis()->SetBinLabel(7,"IDprong diff GetDau");
      hControlDInJ->GetXaxis()->SetBinLabel(8,"still z>1 (cand)");
      
      hControlDInJ->SetNdivisions(1);
      hControlDInJ->GetXaxis()->SetLabelSize(0.05);
      fOutput->Add(hControlDInJ);
      
      TH1F *hmissingp=new TH1F("hmissingp", "Distribution of missing momentum (modulo);|p_{missing}|", 200, 0.,20.);
      fOutput->Add(hmissingp);
      
      for(Int_t k=0;k<3;k++) {
      	 TH1F *hMissPi=new TH1F(Form("hMissP%d",k), Form("Missing p component {%d};p_{%d}",k,k), 400, -10.,10.);
      	 fOutput->Add(hMissPi);
      }
      TH1F *hDeltaPtJet=new TH1F("hDeltaPtJet", "Difference between the jet pt before and after correction;p_{T}^{jet,recalc}-p_{T}^{jet,orig}", 200, 0.,20.);
      
      fOutput->Add(hDeltaPtJet);
      TH1F *hRelDeltaPtJet=new TH1F("hRelDeltaPtJet", "Difference between the jet pt before and after correction/ original pt;(p_{T}^{jet,recalc}-p_{T}^{jet,orig})/p_{T}^{jet,orig}", 200, 0.,1.);
      fOutput->Add(hRelDeltaPtJet);
      
      TH1F* hzDinjet=new TH1F("hzDinjet","Z of candidates with daughters in jet;z",nbinsz,zlims[0],zlims[1]);
      fOutput->Add(hzDinjet);
      //frag func of tracks in the jet
      TH1F* hztracksinjet=new TH1F("hztracksinjet","Z of tracks in jet;z",nbinsz,zlims[0],zlims[1]);
      fOutput->Add(hztracksinjet);
      
      //check jet and tracks
      TH1F* hDiffPtTrPtJzok=new TH1F("hDiffPtTrPtJzok","Sum p_{T}^{trks} - p_{T}^{jet}, for z<1;(#Sigma p_{T}^{trks}) - p_{T}^{jet}", 100,-0.2,0.2);
      fOutput->Add(hDiffPtTrPtJzok);
      TH1F* hDiffPtTrPtJzNok=new TH1F("hDiffPtTrPtJzNok","Sum p_{T}^{trks} - p_{T}^{jet}, for z>1;(#Sigma p_{T}^{trks}) - p_{T}^{jet}", 100,-0.2,0.2);
      fOutput->Add(hDiffPtTrPtJzNok);
      TH1F* hDiffPzTrPtJzok=new TH1F("hDiffPzTrPtJzok","Sum p_{z}^{trks} - p_{z}^{jet}, for z<1;(#Sigma p_{z}^{trks}) - p_{z}^{jet}", 100,-0.2,0.2);
      fOutput->Add(hDiffPzTrPtJzok);
      TH1F* hDiffPzTrPtJzNok=new TH1F("hDiffPzTrPtJzNok","Sum p_{z}^{trks} - p_{z}^{jet}, for z>1;(#Sigma p_{z}^{trks}) - p_{z}^{jet}", 100,-0.2,0.2);
      fOutput->Add(hDiffPzTrPtJzNok);
      TH1I* hNtrkjzNok=new TH1I("hNtrkjzNok", "Number of tracks in a jet with z>1;N^{tracks} (z>1)",5,0,5);
      fOutput->Add(hNtrkjzNok);
      
      //calculate frag func with pt (simply ptD(or track)\cdot pt jet /ptjet^2)
      TH1F* hzDT=new TH1F("hzDT", "Z of D in jet in transverse components;(p_{T}^{D} dot p_{T}^{jet})/p_{T}^{jet}^{2} ",nbinsz,zlims[0],zlims[1]);
      fOutput->Add(hzDT);
      TH1F* hztracksinjetT=new TH1F("hztracksinjetT", "Z of jet tracks in transverse components;(p_{T}^{trks} dot p_{T}^{jet})/p_{T}^{jet}^{2}",nbinsz,zlims[0],zlims[1]);
      fOutput->Add(hztracksinjetT);
      
      TH1I* hIDddaugh=new TH1I("hIDddaugh", "ID of daughters;ID", 2001,-1000,1000);
      fOutput->Add(hIDddaugh);
      TH1I* hIDddaughOut=new TH1I("hIDddaughOut", "ID of daughters not found in jet;ID", 2001,-1000,1000);
      fOutput->Add(hIDddaughOut);
      TH1I* hIDjetTracks=new TH1I("hIDjetTracks", "ID of jet tracks;ID", 2001,-1000,1000);
      fOutput->Add(hIDjetTracks);
      
      TH1F* hDRdaughOut=new TH1F("hDRdaughOut", "#Delta R of daughters not belonging to the jet tracks (D in jet);#Delta R",200, 0.,2.);
      fOutput->Add(hDRdaughOut);
      
      
      if(fCandidateType==kDstartoKpipi) 
      {
      	 
      	 TH2F* hDiffSideBand = new TH2F("hDiffSideBand","M(kpipi)-M(kpi) Side Band Background",nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      	 hDiffSideBand->SetStats(kTRUE);
      	 hDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV");
      	 hDiffSideBand->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	 hDiffSideBand->Sumw2();
      	 fOutput->Add(hDiffSideBand); 
      	 
      	 
      	 TH1F* hPtPion = new TH1F("hPtPion","Primary pions candidates pt ",500,0,10);
      	 hPtPion->SetStats(kTRUE);
      	 hPtPion->GetXaxis()->SetTitle("GeV/c");
      	 hPtPion->GetYaxis()->SetTitle("Entries");
      	 hPtPion->Sumw2();
      	 fOutput->Add(hPtPion);
      	 
      }
      // D related histograms
      TH1F *hNDPerEvNoJet=new TH1F("hNDPerEvNoJet","Number of candidates per event with no jets; N candidate/ev with no jet", 20, 0., 20.);
      hNDPerEvNoJet->Sumw2();
      fOutput->Add(hNDPerEvNoJet);
      
      TH1F *hptDPerEvNoJet=new TH1F("hptDPerEvNoJet","pt distribution of candidates per events with no jets; p_{t}^{D} (GeV/c)",nbinsptD, ptDlims[0],ptDlims[1]);
      hptDPerEvNoJet->Sumw2();
      fOutput->Add(hptDPerEvNoJet);
      
      if(fSwitchOnSparses){
      	 const Int_t    nAxisD=4;
      	 const Int_t    nbinsSparseD[nAxisD]={nbinsSpseta,nbinsSpsphi,nbinsSpsptD,nbinsSpsmass};
      	 const Double_t minSparseD[nAxisD]  ={etalims[0],philims[0],ptDlims[0],fMinMass};
      	 const Double_t maxSparseD[nAxisD]  ={etalims[1],philims[1],ptDlims[1],fMaxMass};
      	 THnSparseF *hsDstandalone=new THnSparseF("hsDstandalone","#phi, #eta, p_{T}^{D}, and mass", nAxisD, nbinsSparseD, minSparseD, maxSparseD);
      	 hsDstandalone->Sumw2();
      	 
      	 fOutput->Add(hsDstandalone);
      }
      
      /*
      TH3F* hPtJetWithD=new TH3F("hPtJetWithD","D-Jet Pt distribution; p_{T} (GeV/c);delta mass (GeV/c^{2})",nbinsptjet,ptjetlims[0],ptjetlims[1],nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      hPtJetWithD->Sumw2();
      //for the MC this histogram is filled with the real background
      TH3F* hPtJetWithDsb=new TH3F("hPtJetWithDsb","D(background)-Jet Pt distribution; p_{T} (GeV/c);delta mass (GeV/c^{2});p_{T}^{D} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1],nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      hPtJetWithDsb->Sumw2();
      
      fOutput->Add(hPtJetWithD);
      fOutput->Add(hPtJetWithDsb);

      */
      TH1F *hNJetPerEvNoD=new TH1F("hNJetPerEvNoD","Number of jets per event with no D; number of jets/ev with no D",10,-0.5,9.5);
      hNJetPerEvNoD->Sumw2();
      
      TH1F *hPtJetPerEvNoD=new TH1F("hPtJetPerEvNoD","pt distribution of jets per event with no D; p_{T}^{jet} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
      hPtJetPerEvNoD->Sumw2();
      
      fOutput->Add(hNJetPerEvNoD);
      fOutput->Add(hPtJetPerEvNoD);
      
      TH1F* hDeltaRD=new TH1F("hDeltaRD","#Delta R distribution of D candidates selected;#Delta R",200, 0.,10.);
      hDeltaRD->Sumw2();
      fOutput->Add(hDeltaRD);
      
      //background (side bands for the Dstar and like sign for D0)
      fJetRadius=GetJetContainer(0)->GetJetRadius();
      TH2F* hInvMassptD = new TH2F("hInvMassptD",Form("D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      hInvMassptD->SetStats(kTRUE);
      hInvMassptD->GetXaxis()->SetTitle("mass (GeV)");
      hInvMassptD->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      hInvMassptD->Sumw2();
      
      fOutput->Add(hInvMassptD);
      
      if(fUseMCInfo){
      	 TH2F* hInvMassptDbg = new TH2F("hInvMassptDbg",Form("Background D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      	 hInvMassptDbg->GetXaxis()->SetTitle(Form("%s (GeV)",(fCandidateType==kDstartoKpipi) ? "M(Kpipi)" : "M(Kpi)"));
      	 hInvMassptDbg->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	 hInvMassptDbg->Sumw2();
      	 fOutput->Add(hInvMassptDbg);
      	 
      }
      
      if(fSwitchOnSparses){
      	 Double_t pi=TMath::Pi(), philims2[2];
      	 philims2[0]=-pi*1./2.;
      	 philims2[1]=pi*3./2.;
      	 THnSparseF *hsDphiz=0x0; //definition below according to the switches
      	 
      	 if(fSwitchOnSB && fSwitchOnPhiAxis && fSwitchOnOutOfConeAxis){
      	    AliInfo("Creating a 9 axes container: This might cause large memory usage");
      	    const Int_t nAxis=9;   
      	    const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsphi,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,2, 2, 2, 2};
      	    const Double_t minSparse[nAxis]={zlims[0],philims2[0],ptjetlims[0],ptDlims[0],fMinMass,-0.5, -0.5,-0.5,-0.5};
      	    const Double_t maxSparse[nAxis]={zlims[1],philims2[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5, 1.5, 1.5 , 1.5};
      	    fNAxesBigSparse=nAxis;
      	    
      	    hsDphiz=new THnSparseF("hsDphiz","Z and #Delta#phi vs p_{T}^{jet}, p_{T}^{D}, mass. SB? D within the jet cone?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
      	 }
      	 
      	 if(!fSwitchOnPhiAxis || !fSwitchOnOutOfConeAxis || !fSwitchOnSB){
      	    fSwitchOnPhiAxis=0;
      	    fSwitchOnOutOfConeAxis=0;
      	    fSwitchOnSB=0;
      	    if(fUseMCInfo){
      	       AliInfo("Creating a 7 axes container (MB background candidates)");
      	       const Int_t nAxis=9;   
      	       const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,2, 2, 2};
      	       const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass, -0.5,-0.5,-0.5};
      	       const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5, 1.5 , 1.5};
      	       fNAxesBigSparse=nAxis;
      	       hsDphiz=new THnSparseF("hsDphiz","Z vs p_{T}^{jet}, p_{T}^{D}, mass. Bkg?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
      	       
      	    } else{
      	       AliInfo("Creating a 6 axes container");
      	       const Int_t nAxis=6;
      	       const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass, 2, 2};
      	       const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,-0.5,-0.5};
      	       const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5, 1.5};
      	       fNAxesBigSparse=nAxis;      	 
      	       
      	       hsDphiz=new THnSparseF("hsDphiz","Z vs p_{T}^{jet}, p_{T}^{D}, mass., D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
      	    }
      	 }
      	 if(!hsDphiz) AliFatal("No THnSparse created");
      	 hsDphiz->Sumw2();
      	 
      	 fOutput->Add(hsDphiz);
      }
   }
   PostData(1, fOutput);
   
   return kTRUE; 
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsRecoJetCorr(AliVParticle* candidate, AliEmcalJet *jet,  AliAODEvent* aodEvent){
   
   Double_t ptD=candidate->Pt();
   Double_t deltaR=DeltaR(candidate,jet);
   Double_t phiD=candidate->Phi();
   Double_t deltaphi = jet->Phi()-phiD;
   if(deltaphi<=-(TMath::Pi())/2.) deltaphi = deltaphi+2.*(TMath::Pi());
   if(deltaphi>(3.*(TMath::Pi()))/2.) deltaphi = deltaphi-2.*(TMath::Pi());
   Double_t z=Z(candidate,jet);
   /*
   if(z>1) {
      ((TH1I*)fOutput->FindObject("hControlDInJ"))->Fill(7);
      Double_t pmissing=TMath::Sqrt(fPmissing[0]*fPmissing[0]+fPmissing[1]*fPmissing[1]+fPmissing[2]*fPmissing[2]);
      
      for(Int_t k=0;k<3;k++) ((TH1F*)fOutput->FindObject(Form("hMissP%d",k)))->Fill(fPmissing[k]);
      
      ((TH1F*)fOutput->FindObject("hmissingp"))->Fill(pmissing);
      Double_t ptdiff=fPtJet-jet->Pt();
      ((TH1F*)fOutput->FindObject("hDeltaPtJet"))->Fill(ptdiff);
      ((TH1F*)fOutput->FindObject("hRelDeltaPtJet"))->Fill(ptdiff/jet->Pt());

      
   }
   */
   if(fIsDInJet)((TH1F*)fOutput->FindObject("hzDT"))->Fill(Z(candidate,jet,kTRUE));
   
   TH1F* hDeltaRD=(TH1F*)fOutput->FindObject("hDeltaRD");
   hDeltaRD->Fill(deltaR);
   Bool_t bDInEMCalAcc=InEMCalAcceptance(candidate);
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);
   if(fUseReco){
      if(fCandidateType==kD0toKpi) {
      	 AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)candidate;
      	 
      	 FillHistogramsD0JetCorr(dzero,deltaphi,z,ptD,fPtJet,bDInEMCalAcc,bJetInEMCalAcc, aodEvent);
      	 
      }
      
      if(fCandidateType==kDstartoKpipi) {
      	 AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)candidate;
      	 FillHistogramsDstarJetCorr(dstar,deltaphi,z,ptD,fPtJet,bDInEMCalAcc,bJetInEMCalAcc);
      	 
      }
   } else{ //generated level
      //AliInfo("Non reco");
      FillHistogramsMCGenDJetCorr(deltaphi,z,ptD,fPtJet,bDInEMCalAcc,bJetInEMCalAcc);
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t ptj, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent){


   Double_t masses[2]={0.,0.};
   Int_t pdgdaughtersD0[2]={211,321};//pi,K 
   Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
   
   masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
   masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
   
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fOutput->FindObject("hsDphiz");
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptj;
      point[3]=ptD;
      point[4]=masses[0];
      point[5]=0;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==6){
      point=new Double_t[6];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[5]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
}
  if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=0;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);

   }
   
   
   Printf("Candidate in FillHistogramsD0JetCorr IsA %s", (candidate->IsA())->GetName());   
   Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);
   if(isselected==1 || isselected==3) {
      
      //if(fIsDInJet) hPtJetWithD->Fill(ptj,masses[0],ptD);
      
      FillMassHistograms(masses[0], ptD);
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
   }
   if(isselected>=2) {
      //if(fIsDInJet) hPtJetWithD->Fill(ptj,masses[1],ptD);
      
      FillMassHistograms(masses[1], ptD);
      if(fNAxesBigSparse==9) point[4]=masses[1];
      if(fNAxesBigSparse==6 || fNAxesBigSparse==7) point[3]=masses[1];
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
   }
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi,  Double_t z, Double_t ptD, Double_t ptj, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){
  //dPhi and z not used at the moment,but will be (re)added

   AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
   Double_t deltamass= dstar->DeltaInvMass(); 
   //Double_t massD0= dstar->InvMassD0();
   
   
   TH1F* hPtPion=(TH1F*)fOutput->FindObject("hPtPion");
   hPtPion->Fill(softpi->Pt());
   
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fOutput->FindObject("hsDphiz");
   Int_t isSB=0;//IsDzeroSideBand(dstar);
   
   //Double_t point[]={z,dPhi,ptj,ptD,deltamass,static_cast<Double_t>(isSB), static_cast<Double_t>(fIsDInJet ? 1 : 0),bDInEMCalAcc,bJetInEMCalAcc};
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptj;
      point[3]=ptD;
      point[4]=deltamass;
      point[5]=static_cast<Double_t>(isSB);
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==6){
      point=new Double_t[6];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=deltamass;
      point[4]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[5]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=deltamass;
      point[4]=0;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   //if(fIsDInJet) hPtJetWithD->Fill(ptj,deltamass,ptD);
   
   FillMassHistograms(deltamass, ptD);
   if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsMCGenDJetCorr(Double_t dPhi,Double_t z,Double_t ptD,Double_t ptjet, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){
   
   Double_t pdgmass=0;
   if(fCandidateType==kD0toKpi) pdgmass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
   if(fCandidateType==kDstartoKpipi) pdgmass=TDatabasePDG::Instance()->GetParticle(413)->Mass()-TDatabasePDG::Instance()->GetParticle(421)->Mass();
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   THnSparseF* hsDphiz=(THnSparseF*)fOutput->FindObject("hsDphiz");
   //Double_t point[9]={z,dPhi,ptjet,ptD,pdgmass,0, static_cast<Double_t>(fIsDInJet ? 1 : 0),bDInEMCalAcc,bJetInEMCalAcc};
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptjet;
      point[3]=ptD;
      point[4]=pdgmass;
      point[5]=0;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==6){
      point=new Double_t[6];
      point[0]=z;
      point[1]=ptjet;
      point[2]=ptD;
      point[3]=pdgmass;
      point[4]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[5]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
      if(fNAxesBigSparse==7){
      point=new Double_t[6];
      point[0]=z;
      point[1]=ptjet;
      point[2]=ptD;
      point[3]=pdgmass;
      point[4]=1;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }


   
   if(fNAxesBigSparse==9) point[4]=pdgmass;
   if(fNAxesBigSparse==6 || fNAxesBigSparse==7) point[3]=pdgmass;
   if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
   //if(fIsDInJet) {
   //  hPtJetWithD->Fill(ptjet,pdgmass,ptD); // candidates within a jet
   //}
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillMassHistograms(Double_t mass,Double_t ptD){
   
   if(fIsDInJet) {
      TH2F* hInvMassptD=(TH2F*)fOutput->FindObject("hInvMassptD");
      hInvMassptD->Fill(mass,ptD);
   }
}

//________________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FlagFlavour(AliEmcalJet *jet){
   
   AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
   if (fCandidateType==kD0toKpi) tag=AliEmcalJet::kD0;
   if (fIsDInJet) jet->AddFlavourTag(tag);
   
   return;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::SideBandBackground(AliAODRecoCascadeHF *candDstar, AliEmcalJet *jet){
   
   //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
   // (expected detector resolution) on the left and right frm the D0 mass. Each band
   //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
   TH2F* hDiffSideBand=(TH2F*)fOutput->FindObject("hDiffSideBand");
   //TH3F* hPtJetWithDsb=(TH3F*)fOutput->FindObject("hPtJetWithDsb");
   THnSparseF* hsDphiz=(THnSparseF*)fOutput->FindObject("hsDphiz");
   
   Bool_t bDInEMCalAcc=InEMCalAcceptance(candDstar);  
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);
   
   Double_t deltaM=candDstar->DeltaInvMass(); 
   //Printf("Inv mass = %f between %f and %f or %f and %f?",invM, sixSigmal,fourSigmal,fourSigmar,sixSigmar);
   Double_t z=Z(candDstar,jet);
   Double_t ptD=candDstar->Pt();

   Double_t dPhi=jet->Phi()-candDstar->Phi();

   if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
   if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());
   
   Int_t isSideBand=1;//IsDzeroSideBand(candDstar);
   //if no SB analysis we should not enter here, so no need of checking the number of axes
   Double_t point[9]={z,dPhi,fPtJet,ptD,deltaM,static_cast<Double_t>(isSideBand), static_cast<Double_t>(fIsDInJet ? 1 : 0),static_cast<Double_t>(bDInEMCalAcc? 1 : 0),static_cast<Double_t>(bJetInEMCalAcc? 1 : 0)};
   //fill the background histogram with the side bands or when looking at MC with the real background
   if(isSideBand==1){
      hDiffSideBand->Fill(deltaM,ptD); // M(Kpipi)-M(Kpi) side band background    
      //hdeltaPhiDjaB->Fill(deltaM,ptD,dPhi);
      if(fSwitchOnSparses) hsDphiz->Fill(point,1.);
      //if(fIsDInJet){  // evaluate in the near side	
      //	 hPtJetWithDsb->Fill(fPtJet,deltaM,ptD);
      //}
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::MCBackground(AliAODRecoDecayHF *candbg,AliEmcalJet* jet){
   
   //need mass, deltaR, pt jet, ptD
   //two cases: D0 and Dstar
   
   Int_t isselected=fCuts->IsSelected(candbg,AliRDHFCuts::kAll, AODEvent());
   if(!isselected) return;
   
   Double_t ptD=candbg->Pt();
   Double_t phiD=candbg->Phi();
   Double_t deltaphi = jet->Phi()-phiD;
   if(deltaphi<=-(TMath::Pi())/2.) deltaphi = deltaphi+2.*(TMath::Pi());
   if(deltaphi>(3.*(TMath::Pi()))/2.) deltaphi = deltaphi-2.*(TMath::Pi());
   Double_t z=Z(candbg,jet);

   Bool_t bDInEMCalAcc=InEMCalAcceptance(candbg);
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

   TH2F* hInvMassptDbg=(TH2F*)fOutput->FindObject("hInvMassptDbg");
   //TH3F* hPtJetWithDsb=(TH3F*)fOutput->FindObject("hPtJetWithDsb");

   THnSparseF* hsDphiz=(THnSparseF*)fOutput->FindObject("hsDphiz");
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=deltaphi;
      point[2]=fPtJet;
      point[3]=ptD;
      point[4]=-999; //set below
      point[5]=1;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=fPtJet;
      point[2]=ptD;
      point[3]=-999; //set below
      point[4]=1;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   if(fCandidateType==kDstartoKpipi){
      AliAODRecoCascadeHF* dstarbg = (AliAODRecoCascadeHF*)candbg;
      Double_t deltaM=dstarbg->DeltaInvMass();
      hInvMassptDbg->Fill(deltaM,ptD);
      //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,deltaM,ptD);
      if(fNAxesBigSparse==9) point[4]=deltaM;
      if(fNAxesBigSparse==6 || fNAxesBigSparse==7) point[3]=deltaM;
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);      
   }
   
   if(fCandidateType==kD0toKpi){
      Double_t masses[2]={0.,0.};
      Int_t pdgdaughtersD0[2]={211,321};//pi,K 
      Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      
      masses[0]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
      masses[1]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
      
      
      if(isselected==1 || isselected==3) {
      	 //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,masses[0],ptD);
      	 hInvMassptDbg->Fill(masses[0],ptD);
      	 if(fNAxesBigSparse==9) point[4]=masses[0];
      	 if(fNAxesBigSparse==6 || fNAxesBigSparse==7) point[3]=masses[0];
      	 if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
     }
      if(isselected>=2) {
      	 //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,masses[1],ptD);
      	 hInvMassptDbg->Fill(masses[1],ptD);
      	 if(fNAxesBigSparse==9) point[4]=masses[1];
      	 if(fNAxesBigSparse==6 || fNAxesBigSparse==7) point[3]=masses[1];
      	 if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) hsDphiz->Fill(point,1.);
      	 
      }
      
      
   }
}

//_______________________________________________________________________________

Float_t AliAnalysisTaskFlavourJetCorrelations::DeltaR(AliVParticle *p1, AliVParticle *p2) const {
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

Int_t AliAnalysisTaskFlavourJetCorrelations::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){
   
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

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::AreDaughtersInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Int_t*& daughOutOfJetID, AliAODTrack**& dtrks, Bool_t fillH){
   //returns true/false and the array of particles not found among jet constituents
   
   TH1I* hControlDInJ=(TH1I*)fOutput->FindObject("hControlDInJ");
   TH1I* hIDddaugh   =(TH1I*)fOutput->FindObject("hIDddaugh");
   TH1I* hIDddaughOut=(TH1I*)fOutput->FindObject("hIDddaughOut");
   TH1I* hIDjetTracks=(TH1I*)fOutput->FindObject("hIDjetTracks");
   
   Int_t daughtersID[3];
   Int_t ndaugh=0;
   Int_t dmatchedID[3]={0,0,0};
   Int_t countmatches=0;
   if (fCandidateType==kDstartoKpipi){
      AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)charm;
      AliAODRecoDecayHF2Prong* dzero=dstar->Get2Prong();
      daughtersID[0]=dzero->GetProngID(0);
      daughtersID[1]=dzero->GetProngID(1);
      daughtersID[2]=dstar->GetBachelor()->GetID();
      ndaugh=3;
     
      dtrks=new AliAODTrack*[3];
      dtrks[0]=(AliAODTrack*)dzero->GetDaughter(0);
      dtrks[1]=(AliAODTrack*)dzero->GetDaughter(1);
      dtrks[2]=(AliAODTrack*)dstar->GetBachelor();
  
      //check
      if(fillH){
      	 if(daughtersID[0]!=((AliAODTrack*)dtrks[0])->GetID() || daughtersID[1]!=((AliAODTrack*)dtrks[1])->GetID())  hControlDInJ->Fill(6);
      	 
      	 hIDddaugh->Fill(daughtersID[0]);
      	 hIDddaugh->Fill(daughtersID[1]);
      	 hIDddaugh->Fill(daughtersID[2]);
      	 
      }
      //Printf("ID daughters %d, %d, %d",daughtersID[0], daughtersID[1], daughtersID[2]);
   }
   
   if (fCandidateType==kD0toKpi){
      daughtersID[0]=charm->GetProngID(0);
      daughtersID[1]=charm->GetProngID(1);
      ndaugh=2;
      if(fillH){
      	 hIDddaugh->Fill(daughtersID[0]);
      	 hIDddaugh->Fill(daughtersID[1]);
      }
      dtrks=new AliAODTrack*[2];
      dtrks[0]=(AliAODTrack*)charm->GetDaughter(0);
      dtrks[1]=(AliAODTrack*)charm->GetDaughter(1);

   }
   
   const Int_t cndaugh=ndaugh;
   daughOutOfJetID=new Int_t[cndaugh];

   Int_t nchtrk=thejet->GetNumberOfTracks();
   for(Int_t itk=0;itk<nchtrk;itk++){
      AliVTrack* tkjet=((AliPicoTrack*)thejet->TrackAt(itk,fTrackArr))->GetTrack();
      Int_t idtkjet=tkjet->GetID();
      if(fillH) hIDjetTracks->Fill(idtkjet);
      for(Int_t id=0;id<ndaugh;id++){
      	 if(idtkjet==daughtersID[id]) {
      	    dmatchedID[id]=daughtersID[id]; //daughter which matches a track in the jet
      	    countmatches++;
      	    
      	 }
      	 
      	 if(countmatches==ndaugh) break;
      }
   }
   //reverse: include in the array the ID of daughters not matching

   for(Int_t id=0;id<ndaugh;id++){
      if(dmatchedID[id]==0) {
      	 daughOutOfJetID[id]=daughtersID[id];
      	 if(fillH) hIDddaughOut->Fill(daughOutOfJetID[id]);
      }
      else daughOutOfJetID[id]=0;
   }
   if(fillH){
      if((ndaugh-countmatches) == 1) hControlDInJ->Fill(0);
      if((ndaugh-countmatches) == 2) hControlDInJ->Fill(1);
      if((ndaugh-countmatches) == 3) hControlDInJ->Fill(2);
   }
   if(ndaugh!=countmatches){
      return kFALSE;
   }
   
   return kTRUE;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::IsDInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Bool_t fillH){
   
   //check the conditions type:
   //type 0 : DeltaR < jet radius, don't check daughters (and don't correct pt jet) 
   //type 1 : DeltaR < jet radius and check for all daughters among jet tracks, don't correct ptjet
   //type 2 (default) : DeltaR < jet radius and check for all daughters among jet tracks, if not present, correct the ptjet
    
   TH1I* hControlDInJ=(TH1I*)fOutput->FindObject("hControlDInJ");
   TH1F* hDRdaughOut=(TH1F*)fOutput->FindObject("hDRdaughOut");
   
   fPmissing[0]=0; 
   fPmissing[1]=0;
   fPmissing[2]=0;
   
   Bool_t testDeltaR=kFALSE;
   if(DeltaR(thejet,charm)<fJetRadius) testDeltaR=kTRUE;
   
   Int_t* daughOutOfJet;
   AliAODTrack** charmDaugh;
   Bool_t testDaugh=AreDaughtersInJet(thejet, charm, daughOutOfJet,charmDaugh,fillH);
   if(testDaugh && testDeltaR) {
      //Z of candidates with daughters in the jet
      ((TH1F*)fOutput->FindObject("hzDinjet"))->Fill(Z(charm,thejet));
      
   }
   if(!testDaugh && testDeltaR && fTypeDInJet==2){
      
      Int_t ndaugh=3;
      if(fCandidateType==kD0toKpi) ndaugh=2;
      Int_t nOut=ndaugh;
      
      for(Int_t id=0;id<ndaugh;id++){
      	 if(daughOutOfJet[id]!=0){ //non-matched daughter
      	    //get track and its momentum
      	    nOut--;
      	    if(fillH) {
      	       hControlDInJ->Fill(3);
      	       if(id<2) hControlDInJ->Fill(4);
      	       if(id==2)hControlDInJ->Fill(5);
      	       hDRdaughOut->Fill(DeltaR(thejet, charmDaugh[id]));
      	    }
      	    fPmissing[0]+=charmDaugh[id]->Px(); 
      	    fPmissing[1]+=charmDaugh[id]->Py();
      	    fPmissing[2]+=charmDaugh[id]->Pz();
      	 }
      
      }
      
      //now the D in within the jet
      testDaugh=kTRUE;
   
   }
   
   delete[] charmDaugh;
   
   Bool_t result=0;
   switch(fTypeDInJet){
   case 0:
      result=testDeltaR;
      break;
   case 1:
      result=testDeltaR && testDaugh;
      break;
   case 2:
      result=testDeltaR && testDaugh;
      break;
   default:
      AliInfo("Selection type not specified, use 1");
      result=testDeltaR && testDaugh;
      break;
   }
 return result;
}

Bool_t AliAnalysisTaskFlavourJetCorrelations::InEMCalAcceptance(AliVParticle *vpart){
   //check eta phi of a VParticle: return true if it is in the EMCal acceptance, false otherwise
   
   Double_t phiEMCal[2]={1.405,3.135},etaEMCal[2]={-0.7,0.7};
   Bool_t binEMCal=kTRUE;
   Double_t phi=vpart->Phi(), eta=vpart->Eta();
   if(phi<phiEMCal[0] || phi>phiEMCal[1]) binEMCal=kFALSE;
   if(eta<etaEMCal[0] || eta>etaEMCal[1]) binEMCal=kFALSE;
   return binEMCal;


}
