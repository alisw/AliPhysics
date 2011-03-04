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
//             Base class for DStar in Jets Analysis
//
//-----------------------------------------------------------------------
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include "TROOT.h"

#include "AliAnalysisTaskSEDStarJets.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskSEDStarJets)

//__________________________________________________________________________
AliAnalysisTaskSEDStarJets::AliAnalysisTaskSEDStarJets() :
  AliAnalysisTaskSE(),
  fEvents(0),
  fchargeFrCorr(0),
  fUseMCInfo(kTRUE), 
  fRequireNormalization(kTRUE),
  fOutput(0),
  fCuts(0),          
  ftrigger(0),   
  fPtPion(0),        
  fInvMass(0),       
  fRECOPtDStar(0), 
  fRECOPtBkg(0),   
  fDStar(0),          
  fDiff(0),           
  fDiffSideBand(0),  
  fDStarMass(0),    
  fPhi(0),       
  fPhiBkg(0),        
  fTrueDiff(0),       
  fResZ(0),        
  fResZBkg(0),               
  fEjet(0),        
  fPhijet(0),        
  fEtaJet(0),         
  theMCFF(0),
  fDphiD0Dstar(0),
  fPtJet(0)           
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarJets::AliAnalysisTaskSEDStarJets(const Char_t* name, AliRDHFCutsDStartoKpipi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fchargeFrCorr(0),
  fUseMCInfo(kTRUE),
  fRequireNormalization(kTRUE),
  fOutput(0),
  fCuts(0),          
  ftrigger(0),   
  fPtPion(0),        
  fInvMass(0),       
  fRECOPtDStar(0),
  fRECOPtBkg(0),     
  fDStar(0),          
  fDiff(0),           
  fDiffSideBand(0),  
  fDStarMass(0),    
  fPhi(0),       
  fPhiBkg(0),        
  fTrueDiff(0),       
  fResZ(0),        
  fResZBkg(0),       
  fEjet(0),        
  fPhijet(0),        
  fEtaJet(0),         
  theMCFF(0),
  fDphiD0Dstar(0),
  fPtJet(0)               
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  fCuts=cuts;
  Info("AliAnalysisTaskSEDStarJets","Calling Constructor");
 
  DefineOutput(1,TList::Class()); // histos
  DefineOutput(2,AliRDHFCutsDStartoKpipi::Class()); // my cuts
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarJets::~AliAnalysisTaskSEDStarJets() {
  //
  // destructor
  //

  Info("~AliAnalysisTaskSEDStarJets","Calling Destructor");  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  
  if (ftrigger) { delete ftrigger; ftrigger = 0;} 
  if (fPtPion)  { delete fPtPion;  fPtPion = 0;} 
  if (fInvMass) { delete fInvMass; fInvMass = 0;} 
  if (fRECOPtDStar) { delete fRECOPtDStar; fRECOPtDStar = 0;} 
  if (fRECOPtBkg)   { delete fRECOPtBkg; fRECOPtBkg = 0;} 
  if (fDStar) { delete fDStar; fDStar = 0;} 
  if (fDiff)  { delete fDiff; fDiff = 0;} 
  if (fDiffSideBand) { delete fDiffSideBand; fDiffSideBand = 0;} 
  if (fDStarMass)    { delete fDStarMass; fDStarMass = 0;} 
  if (fPhi)     { delete fPhi; fPhi = 0;} 
  if (fPhiBkg)  { delete fPhiBkg; fPhiBkg = 0;} 
  if (fTrueDiff){ delete fTrueDiff; fTrueDiff = 0;} 
  if (fResZ)    { delete fResZ;  fResZ = 0;} 
  if (fResZBkg) { delete fResZBkg; fResZBkg = 0;}
  if (fEjet)    { delete fEjet; fEjet = 0;}
  if (fPhijet)  { delete fPhijet; fPhijet = 0;}
  if (fEtaJet)  { delete fEtaJet; fEtaJet = 0;} 
  if (theMCFF)  { delete theMCFF; theMCFF = 0;} 
  if (fDphiD0Dstar) { delete fDphiD0Dstar; fDphiD0Dstar = 0;} 
  if (fPtJet) { delete fPtJet; fPtJet = 0;} 
}

//___________________________________________________________
void AliAnalysisTaskSEDStarJets::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskSEDStarJets::Init() \n");
  AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
  // Post the cuts
  PostData(2,copyfCuts);
  
  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarJets::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  fEvents++;
  AliInfo(Form("Event %d",fEvents));
  if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));

  // Load the event
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
 
  TClonesArray *arrayDStartoD0pi=0;

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
  
  if (!arrayDStartoD0pi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;

  // Simulate a jet triggered sample
  TClonesArray *arrayofJets = (TClonesArray*)aodEvent->GetJets();
  if(aodEvent->GetNJets()<=0) return;

  // counters for efficiencies
  Int_t icountReco = 0;
  
  // Normalization factor
  if(fRequireNormalization){       
    ftrigger->Fill(1);
  }
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  Double_t max =0;
  Double_t ejet   = 0;
  Double_t phiJet = 0;
  Double_t etaJet = 0;
  Double_t ptjet = 0;

  //loop over jets and consider only the leading jey in the event
  for (Int_t iJets = 0; iJets<arrayofJets->GetEntriesFast(); iJets++) {    
    AliAODJet* jet = (AliAODJet*)arrayofJets->At(iJets);
     
    //jets variables
    ejet   = jet->E();

    if(ejet>max){
      max = jet->E();
      phiJet = jet->Phi();
      etaJet = jet->Eta();
      ptjet = jet->Pt();
    }
    
    // fill energy, eta and phi of the jet
    fEjet   ->Fill(ejet);
    fPhijet ->Fill(phiJet);
    fEtaJet ->Fill(etaJet);
    fPtJet->Fill(ptjet);
  }

  //loop over D* candidates
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {
    
    // D* candidates
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();

    Double_t finvM =-999;
    Double_t finvMDStar =-999;
    Double_t dPhi =-999;
    Bool_t isDStar =0;
    Int_t mcLabel = -9999;

    // find associated MC particle for D* ->D0toKpi
    if(fUseMCInfo){
      TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!mcArray) {
	printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
	return;
      }
      mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);

      if(mcLabel>=0) isDStar = 1;
      if(mcLabel>0){
	Double_t zMC =-999;
	AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
      //fragmentation function in mc
	zMC= FillMCFF(mcPart,mcArray,mcLabel);
	if(zMC>0) theMCFF->Fill(zMC);
      }	      
    }

    // soft pion
    AliAODTrack *track2 = (AliAODTrack*)dstarD0pi->GetBachelor();              
    Double_t pt = dstarD0pi->Pt();

    // track quality cuts
    Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
    if(!isTkSelected) continue;

    // region of interest + topological cuts + PID
    if(!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(),dstarD0pi->YDstar())) continue;    
    Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
    if (!isSelected) continue;
    
    // fill histos
    finvM = dstarD0pi->InvMassD0();
    fInvMass->Fill(finvM); 

    //DStar invariant mass
    finvMDStar = dstarD0pi->InvMassDstarKpipi();
    
    Double_t EGjet = 0;
    Double_t dStarMom = dstarD0pi->P();
    Double_t phiDStar = dstarD0pi->Phi(); 
    Double_t phiD0    = theD0particle->Phi();
    //check suggested by Federico
    Double_t dPhiD0Dstar = phiD0 - phiDStar;
    
    dPhi = phiJet - phiDStar;

    //plot right dphi
    if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
    if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());
    
    Double_t corrFactorCharge = (ejet/100)*20;
    EGjet =  ejet + corrFactorCharge;	
    
    // fill D* candidates
    Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    if(finvM >= (mPDGD0-0.05) && finvM <=(mPDGD0+0.05)){ // ~3 sigma (sigma=17MeV, conservative)

      if(isDStar == 1) { 
	fDphiD0Dstar->Fill(dPhiD0Dstar);
	fDStarMass->Fill(finvMDStar); 
	fTrueDiff->Fill(finvMDStar-finvM);
      }
      if(isDStar == 0) fDphiD0Dstar->Fill(dPhiD0Dstar); // angle between D0 and D*

      fDStar->Fill(finvMDStar);
      fDiff ->Fill(finvMDStar-finvM);  	      
      
      Double_t mPDGDstar=TDatabasePDG::Instance()->GetParticle(413)->Mass();
      Double_t invmassDelta = dstarD0pi->DeltaInvMass();
      
      // now the dphi signal and the fragmentation function 
      if (TMath::Abs(invmassDelta-(mPDGDstar-mPDGD0))<0.0019){
	//fill candidates D* and soft pion reco pt
	
	fRECOPtDStar->Fill(pt);
	fPtPion->Fill(track2->Pt());
	
	fPhi ->Fill(dPhi);

	Double_t jetCone = 0.4;
	if(dPhi>=-jetCone && dPhi<=jetCone){  // evaluate in the near side inside UA1 radius		      
	  Double_t zFrag = (TMath::Cos(dPhi)*dStarMom)/EGjet;                               
	  fResZ->Fill(TMath::Abs(zFrag));		      
	}		    
      }
    }    
    // evaluate side band background
    SideBandBackground(finvM, finvMDStar, dStarMom, EGjet, dPhi);
    
  } // D* cand

AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));

PostData(1,fOutput);
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarJets::Terminate(Option_t*)
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
  
  fDStarMass    = dynamic_cast<TH1F*>(fOutput->FindObject("fDStarMass"));
  fTrueDiff     = dynamic_cast<TH1F*>(fOutput->FindObject("fTrueDiff"));
  fInvMass      = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass"));
  fPtPion       = dynamic_cast<TH1F*>(fOutput->FindObject("fPtPion "));
  fDStar        = dynamic_cast<TH1F*>(fOutput->FindObject("fDStar"));
  fDiff         = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff"));
  fDiffSideBand = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand"));
  ftrigger      = dynamic_cast<TH1F*>(fOutput->FindObject("ftrigger"));
  fRECOPtDStar  = dynamic_cast<TH1F*>(fOutput->FindObject("fRECOPtDStar"));
  fRECOPtBkg    = dynamic_cast<TH1F*>(fOutput->FindObject("fRECOPtBkg"));
  fEjet         = dynamic_cast<TH1F*>(fOutput->FindObject("fEjet"));
  fPhijet       = dynamic_cast<TH1F*>(fOutput->FindObject("fPhijet"));
  fEtaJet       = dynamic_cast<TH1F*>(fOutput->FindObject("fEtaJet"));
  fPhi          = dynamic_cast<TH1F*>(fOutput->FindObject("fPhi"));
  fResZ         = dynamic_cast<TH1F*>(fOutput->FindObject("fResZ"));
  fResZBkg      = dynamic_cast<TH1F*>(fOutput->FindObject("fResZBkg"));
  fPhiBkg       = dynamic_cast<TH1F*>(fOutput->FindObject("fPhiBkg"));
  theMCFF       = dynamic_cast<TH1F*>(fOutput->FindObject("theMCFF"));
  fDphiD0Dstar  = dynamic_cast<TH1F*>(fOutput->FindObject("fDphiD0Dstar"));
  fPtJet  = dynamic_cast<TH1F*>(fOutput->FindObject("fPtJet"));

}
//___________________________________________________________________________

void AliAnalysisTaskSEDStarJets::UserCreateOutputObjects() { 
 // output 
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  // define histograms
  DefineHistoFroAnalysis();

  return;
}
//___________________________________ hiostograms _______________________________________

Bool_t  AliAnalysisTaskSEDStarJets::DefineHistoFroAnalysis(){
  
  // Invariant mass related histograms
  fInvMass = new TH1F("invMass","Kpi invariant mass distribution",1500,.5,3.5);
  fInvMass->SetStats(kTRUE);
  fInvMass->GetXaxis()->SetTitle("GeV/c");
  fInvMass->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fInvMass);
  
  fDStar = new TH1F("invMassDStar","DStar invariant mass after D0 cuts ",600,1.8,2.4);
  fDStar->SetStats(kTRUE);
  fDStar->GetXaxis()->SetTitle("GeV/c");
  fDStar->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDStar);

  fDiff = new TH1F("Diff","M(kpipi)-M(kpi)",750,0.1,0.2);
  fDiff->SetStats(kTRUE);
  fDiff->GetXaxis()->SetTitle("M(kpipi)-M(kpi) GeV/c^2");
  fDiff->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDiff);
  
  fDiffSideBand = new TH1F("DiffSide","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand->SetStats(kTRUE);
  fDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDiffSideBand); 
 
  fDStarMass = new TH1F("RECODStar2","RECO DStar invariant mass distribution",750,1.5,2.5);
  fOutput->Add(fDStarMass);

  fTrueDiff  = new TH1F("dstar","True Reco diff",750,0,0.2);
  fOutput->Add(fTrueDiff);

  // trigger normalization
  ftrigger = new TH1F("Normalization","Normalization factor for correlations",1,0,10);
  ftrigger->SetStats(kTRUE);
  fOutput->Add(ftrigger);

  //correlation fistograms
  fPhi = new TH1F("phi","Delta phi between Jet axis and DStar ",25,-1.57,4.72);
  fPhi->SetStats(kTRUE);
  fPhi->GetXaxis()->SetTitle("#Delta #phi (rad)");
  fPhi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPhi);

  fDphiD0Dstar = new TH1F("phiD0Dstar","Delta phi between D0 and DStar ",1000,-6.5,6.5);
  fOutput->Add(fDphiD0Dstar);

  fPhiBkg = new TH1F("phiBkg","Delta phi between Jet axis and DStar background ",25,-1.57,4.72);
  fPhiBkg->SetStats(kTRUE);
  fPhiBkg->GetXaxis()->SetTitle("#Delta #phi (rad)");
  fPhiBkg->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPhiBkg);

  fRECOPtDStar = new TH1F("RECODStar1","RECO DStar pt distribution",30,0,30);
  fRECOPtDStar->SetStats(kTRUE);
  fRECOPtDStar->SetLineColor(2);
  fRECOPtDStar->GetXaxis()->SetTitle("GeV/c");
  fRECOPtDStar->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fRECOPtDStar);
  
  fRECOPtBkg = new TH1F("RECOptBkg","RECO pt distribution side bands",30,0,30);
  fRECOPtBkg->SetStats(kTRUE);
  fRECOPtBkg->SetLineColor(2);
  fRECOPtBkg->GetXaxis()->SetTitle("GeV/c");
  fRECOPtBkg->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fRECOPtBkg);

  fPtPion = new TH1F("pionpt","Primary pions candidates pt ",500,0,10);
  fPtPion->SetStats(kTRUE);
  fPtPion->GetXaxis()->SetTitle("GeV/c");
  fPtPion->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPtPion);
    
  // jet related fistograms
  fEjet      = new TH1F("ejet",  "UA1 algorithm jet energy distribution",1000,0,500);
  fPhijet    = new TH1F("Phijet","UA1 algorithm jet #phi distribution",  200,-7,7);
  fEtaJet    = new TH1F("Etajet","UA1 algorithm jet #eta distribution",  200,-7,7);
  fPtJet      = new TH1F("PtJet",  "UA1 algorithm jet Pt distribution",1000,0,500);
  fOutput->Add(fEjet);
  fOutput->Add(fPhijet);
  fOutput->Add(fEtaJet);
  fOutput->Add(fPtJet);

  theMCFF    = new TH1F("FragFuncMC","Fragmentation function in MC for FC ",100,0,10);
  fResZ      = new TH1F("FragFunc","Fragmentation function ",50,0,1);
  fResZBkg   = new TH1F("FragFuncBkg","Fragmentation function background",50,0,1);  
  fOutput->Add(theMCFF);
  fOutput->Add(fResZ);
  fOutput->Add(fResZBkg);

  return kTRUE; 
}

//______________________________ side band background for D*___________________________________

void AliAnalysisTaskSEDStarJets::SideBandBackground(Double_t invM, Double_t invMDStar, Double_t dStarMomBkg, Double_t EGjet, Double_t dPhi){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
  
  if((invM>=1.7 && invM<=1.8) || (invM>=1.92 && invM<=2.19)){    
    fDiffSideBand->Fill(invMDStar-invM); // M(Kpipi)-M(Kpi) side band background    
    if ((invMDStar-invM)<=0.14732 && (invMDStar-invM)>=0.14352) {                                                  
      fPhiBkg->Fill(dPhi);
      fRECOPtBkg->Fill(dStarMomBkg);
      if(dPhi>=-0.4 && dPhi<=0.4){  // evaluate in the near side	
	Double_t zFragBkg = (TMath::Cos(dPhi)*dStarMomBkg)/EGjet;                               
	fResZBkg->Fill(TMath::Abs(zFragBkg));	
      }
    }
  }
}

//_____________________________________________________________________________________________-
double AliAnalysisTaskSEDStarJets::FillMCFF(AliAODMCParticle* mcPart, TClonesArray* mcArray, Int_t mcLabel){
  //
  // GS from MC
  // UA1 jet algorithm reproduced in MC
  //
  Double_t zMC2 =-999;
  
  Double_t leading =0;
  Double_t PartE = 0;
  Double_t PhiLeading = -999;
  Double_t EtaLeading = -999;
  Double_t PtLeading = -999;
  Int_t counter =-999;

  //find leading particle
  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
    AliAODMCParticle* Part = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!Part) {
      AliWarning("MC Particle not found in tree, skipping"); 
      continue;
    } 
   
    // remove quarks and the leading particle (it will be counted later)
    if(iPart == mcLabel) continue;
    if(iPart <= 8) continue;
    
    //remove resonances not directly detected in detector
    Int_t PDGCode = Part->GetPdgCode();

    // be sure the particle reach the detector
    Double_t x = Part->Xv();
    Double_t y = Part->Yv();
    Double_t z = Part->Zv();
    
    if(TMath::Abs(PDGCode)== 2212 && x<3 && y<3) continue;
    if(TMath::Abs(x)>30 || TMath::Abs(y)>30 || TMath::Abs(z)>30 ) continue; 
    if(TMath::Abs(PDGCode)!=211 && TMath::Abs(PDGCode)!=321 &&  TMath::Abs(PDGCode)!=11 &&  TMath::Abs(PDGCode)!=13 &&  TMath::Abs(PDGCode)!=2212) continue;

    Int_t daug0 = -999;
    Double_t xd =-999;
    Double_t yd =-999;
    Double_t zd =-999;
    
    daug0 = Part->GetDaughter(0);
    
    if(daug0>=0){
      AliAODMCParticle* tdaug = dynamic_cast<AliAODMCParticle*>(mcArray->At(daug0));
      if(tdaug){ 
      xd = tdaug->Xv();
      yd = tdaug->Yv();
      zd = tdaug->Zv();
      }
    }
    if(TMath::Abs(xd)<3 || TMath::Abs(yd)<3) continue;

    Bool_t AliceAcc = (TMath::Abs(Part->Eta()) <= 0.9);
    if(!AliceAcc) continue; 

    PartE  = Part->E();

    if(PartE>leading){
      leading = Part->E();
      PhiLeading = Part->Phi();
      EtaLeading = Part->Eta();
      PtLeading = Part->Pt();
      counter = iPart;
    } 
 
  }

  Double_t jetEnergy = 0;

  //reconstruct the jet
  for (Int_t iiPart=0; iiPart<mcArray->GetEntriesFast(); iiPart++) { 
    AliAODMCParticle* tPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iiPart));
    if (!tPart) {
      AliWarning("MC Particle not found in tree, skipping"); 
      continue;
    } 
    // remove quarks and the leading particle (it will be counted later)
    if(iiPart == counter) continue; // do not count again the leading particle
    if(iiPart == mcLabel) continue;
    if(iiPart <= 8) continue;

    //remove resonances not directly detected in detector
    Int_t PDGCode = tPart->GetPdgCode();

    // be sure the particle reach the detector
    Double_t x = tPart->Xv();
    Double_t y = tPart->Yv();
    Double_t z = tPart->Zv();
    
    if(TMath::Abs(PDGCode)== 2212 && (x<3 && y<3)) continue;    
    if(TMath::Abs(x)>30 || TMath::Abs(y)>30 || TMath::Abs(z)>30 ) continue; // has to be generated at least in the silicon tracker or beam pipe
    if(TMath::Abs(PDGCode)!=211 && TMath::Abs(PDGCode)!=321 &&  TMath::Abs(PDGCode)!=11 &&  TMath::Abs(PDGCode)!=13 &&  TMath::Abs(PDGCode)!=2212) continue;


    Int_t daug0 = -999;
    Double_t xd =-999;
    Double_t yd =-999;
    Double_t zd =-999;

    daug0 = tPart->GetDaughter(0);

    if(daug0>=0){
      AliAODMCParticle* tdaug = dynamic_cast<AliAODMCParticle*>(mcArray->At(daug0));
      if(tdaug){ 
      xd = tdaug->Xv();
      yd = tdaug->Yv();
      zd = tdaug->Zv();
      }
    }
    if(TMath::Abs(xd)<3 && TMath::Abs(yd)<3) continue;
    //remove particles not in ALICE acceptance
    if(tPart->Pt()<0.07) continue;
    Bool_t AliceAcc = (TMath::Abs(tPart->Eta()) <= 0.9); 
    if(!AliceAcc) continue; 

    Double_t EtaMCp = tPart->Eta();
    Double_t PhiMCp = tPart->Phi();

    Double_t DphiMClead = PhiLeading-PhiMCp;

    if(DphiMClead<=-(TMath::Pi())/2) DphiMClead = DphiMClead+2*(TMath::Pi());
    if(DphiMClead>(3*(TMath::Pi()))/2) DphiMClead = DphiMClead-2*(TMath::Pi());

    Double_t deta = (EtaLeading-EtaMCp);
    //cone radius
    Double_t radius = TMath::Sqrt((DphiMClead*DphiMClead)+(deta*deta));

    if(radius>0.4) continue; // in the jet cone
    if(tPart->Charge()==0) continue; // only charged fraction

    jetEnergy = jetEnergy+(tPart->E());
  }

  jetEnergy = jetEnergy + leading;

  // delta phi D*, jet axis
  Double_t dPhi = PhiLeading - (mcPart->Phi());

  //plot right dphi
  if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
  if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());

  zMC2 = (TMath::Cos(dPhi)*(mcPart->P()))/jetEnergy;

  return zMC2; 
}
