/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <fstream>
#include <map>

// root headers
#include <TMath.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TObjString.h>
#include <TList.h>
#include <TChain.h>

// aliroot headers
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>
#include <AliAODInputHandler.h>
#include <AliMuonTrackCuts.h>

// my headers
#include "AliAnalysisTaskNanoMUON.h"
// ----------------------------------------------------------------------------------------------------------------------------------
class AliAnalysisTaskNanoMUON;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskNanoMUON) // classimp: necessary for root

// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoMUON::AliAnalysisTaskNanoMUON() : AliAnalysisTaskSE(), 
  fMuonTrackCuts(0x0), fPeriod(0), fTrigger(0), fIsMC(0), fIsScalingOn(0), fAOD(0), fMC(0), fOutputList(0),fCounterH(0), fNumberMuonsH(0), fNumberMCMuonsH(0),  
  fMapRunAndLumi(), fMapAnalysedMC(),
  // fRAbsMuonH(0), fMuMuMassPtH(0),  
  // fZNAEnergyTimingH(0), fZNCEnergyTimingH(0), fZNATDCTimingH(0), fZNCTDCTimingH(0),
  // fZNAEnergyTimingAllH(0), fZNCEnergyTimingAllH(0), fZNATDCTimingAllH(0), fZNCTDCTimingAllH(0),
  // fZNAEnergy0NH(0), fZNCEnergy0NH(0)
  fRecTree(0), fRunNum(0), fL0inputs(0),
  // fTracklets(0),
  fZNCEnergy(-999), fZNAEnergy(-999), 
  // fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10),  fV0AFiredCells(-10), fV0CFiredCells(-10), fADADecision(-10), fADCDecision(-10), fIsZNAFired(-10), fIsZNCFired(-10),
  // fIR1Map(0), fIR2Map(0), 
  fMuMuPt(0), 
  // fMuMuPhi(0), 
  fMuMuY(0), fMuMuM(0), 
  // fMuPt1(0), fMuPt2(0), fMuEta1(0), fMuEta2(0), fMuPhi1(0), fMuPhi2(0), fMuQ1(0), fMuQ2(0),
  fGenPart(0), fGenTree(0), fMCRunNum(0), fMCMuMuPt(0),
  // fMCMuMuPhi(0), 
  fMCMuMuY(0), fMCMuMuM(0),
  // fMCMuPt1(0), fMCMuPt2(0), fMCMuEta1(0), fMCMuEta2(0), fMCMuPhi1(0), fMCMuPhi2(0), fMCMuPDG1(0), fMCMuPDG2(0),
  fTrgTree(0), fTrgRunNum(0), fCMUP6Decision(-10), fCMUP10Decision(-10), fCMUP11Decision(-10)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoMUON::AliAnalysisTaskNanoMUON(const char* name) : AliAnalysisTaskSE(name),
  fMuonTrackCuts(0x0), fPeriod(0), fTrigger(0), fIsMC(0), fIsScalingOn(0), fAOD(0), fMC(0), fOutputList(0),fCounterH(0), fNumberMuonsH(0), fNumberMCMuonsH(0), 
  fMapRunAndLumi(), fMapAnalysedMC(),
  // fRAbsMuonH(0), fMuMuMassPtH(0),  
  // fZNAEnergyTimingH(0), fZNCEnergyTimingH(0), fZNATDCTimingH(0), fZNCTDCTimingH(0),
  // fZNAEnergyTimingAllH(0), fZNCEnergyTimingAllH(0), fZNATDCTimingAllH(0), fZNCTDCTimingAllH(0),
  // fZNAEnergy0NH(0), fZNCEnergy0NH(0)
  fRecTree(0), fRunNum(0), fL0inputs(0),
  // fTracklets(0),
  fZNCEnergy(-999), fZNAEnergy(-999), 
  // fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10),  fV0AFiredCells(-10), fV0CFiredCells(-10), fADADecision(-10), fADCDecision(-10), fIsZNAFired(-10), fIsZNCFired(-10),
  // fIR1Map(0), fIR2Map(0), 
  fMuMuPt(0), 
  // fMuMuPhi(0), 
  fMuMuY(0), fMuMuM(0), 
  // fMuPt1(0), fMuPt2(0), fMuEta1(0), fMuEta2(0), fMuPhi1(0), fMuPhi2(0), fMuQ1(0), fMuQ2(0),
  fGenPart(0), fGenTree(0), fMCRunNum(0), fMCMuMuPt(0),
  // fMCMuMuPhi(0), 
  fMCMuMuY(0), fMCMuMuM(0),
  // fMCMuPt1(0), fMCMuPt2(0), fMCMuEta1(0), fMCMuEta2(0), fMCMuPhi1(0), fMCMuPhi2(0), fMCMuPDG1(0), fMCMuPDG2(0),
  fTrgTree(0), fTrgRunNum(0), fCMUP6Decision(-10), fCMUP10Decision(-10), fCMUP11Decision(-10)
{
  // constructor
  DefineInput(0, TChain::Class());   
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());   
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class()); 
}
// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoMUON::~AliAnalysisTaskNanoMUON()
{
  // destructor
  // liberate all allocated memory
  if(fOutputList) {delete fOutputList;}     	
  if(fMuonTrackCuts) {delete fMuonTrackCuts;}
  if(fRecTree) {delete fRecTree;}
  if(fGenTree) {delete fGenTree;}
  if(fTrgTree) {delete fTrgTree;}
  if(fCounterH) {delete fCounterH;}
  if(fNumberMuonsH) {delete fNumberMuonsH;}
  if(fNumberMCMuonsH) {delete fNumberMCMuonsH;}
  // if(fRAbsMuonH) {delete fRAbsMuonH;}
  // if(fMuMuMassPtH) {delete fMuMuMassPtH;}
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)

  ////////////////////////////////////////
  //Muon track cuts
  ////////////////////////////////////////
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchLpt);	
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");

  ////////////////////////////////////////
  //Analysed output tree
  ////////////////////////////////////////
  fRecTree = new TTree("fRecTree", "fRecTree");
  fRecTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fRecTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  // fRecTree ->Branch("fTracklets", &fTracklets, "fTracklets/I");	
  fRecTree ->Branch("fZNCEnergy", &fZNCEnergy, "fZNCEnergy/F");  
  fRecTree ->Branch("fZNAEnergy", &fZNAEnergy, "fZNAEnergy/F");
  // fRecTree ->Branch("fZPCEnergy", &fZPCEnergy, "fZPCEnergy/D");
  // fRecTree ->Branch("fZPAEnergy", &fZPAEnergy, "fZPAEnergy/D");  
  fRecTree ->Branch("fZNATDC", &fZNATDC[0], "fZNATDC[4]/F");
  fRecTree ->Branch("fZNCTDC", &fZNCTDC[0], "fZNCTDC[4]/F");  
  // fRecTree ->Branch("fZPATDC", &fZPATDC[0], "fZPATDC[4]/D");
  // fRecTree ->Branch("fZPCTDC", &fZPCTDC[0], "fZPCTDC[4]/D"); 
  fRecTree ->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
  fRecTree ->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");
  fRecTree ->Branch("fV0AFiredCells", &fV0AFiredCells, "fV0AFiredCells/I");
  fRecTree ->Branch("fV0CFiredCells", &fV0CFiredCells, "fV0CFiredCells/I");
  fRecTree ->Branch("fADADecision", &fADADecision, "fADADecision/I");
  fRecTree ->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
  fRecTree ->Branch("fIsZNAFired", &fIsZNAFired, "fIsZNAFired/I");
  fRecTree ->Branch("fIsZNCFired", &fIsZNCFired, "fIsZNCFired/I");
  // fRecTree ->Branch("fIR1Map", &fIR1Map);
  // fRecTree ->Branch("fIR2Map", &fIR2Map);  
  fRecTree ->Branch("fMuMuPt", &fMuMuPt, "fMuMuPt/F");
  // fRecTree ->Branch("fMuMuPhi", &fMuMuPhi, "fMuMuPhi/D");
  fRecTree ->Branch("fMuMuY", &fMuMuY, "fMuMuY/F");
  fRecTree ->Branch("fMuMuM", &fMuMuM, "fMuMuM/F");
  // fRecTree ->Branch("fMuPt1", &fMuPt1, "fMuPt1/D");
  // fRecTree ->Branch("fMuPt2", &fMuPt2, "fMuPt2/D");
  // fRecTree ->Branch("fMuEta1", &fMuEta1, "fMuEta1/D");
  // fRecTree ->Branch("fMuEta2", &fMuEta2, "fMuEta2/D");
  // fRecTree ->Branch("fMuPhi1", &fMuPhi1, "fMuPhi1/D");
  // fRecTree ->Branch("fMuPhi2", &fMuPhi2, "fMuPhi2/D");
  // fRecTree ->Branch("fMuQ1", &fMuQ1, "fMuQ1/D");
  // fRecTree ->Branch("fMuQ2", &fMuQ2, "fMuQ2/D");
  // post data
  PostData(1, fRecTree);

  ////////////////////////////////////////
  //MC generated output tree
  ////////////////////////////////////////
  fGenTree = new TTree("fGenTree", "fGenTree");
  if(fIsMC){
    fGenPart = new TClonesArray("TParticle", 1000);
    fGenTree ->Branch("fMCRunNum", &fMCRunNum, "fMCRunNum/I");
    fGenTree ->Branch("fMCMuMuPt", &fMCMuMuPt, "fMCMuMuPt/F");
    // fGenTree ->Branch("fMCMuMuPhi", &fMCMuMuPhi, "fMCMuMuPhi/D");
    fGenTree ->Branch("fMCMuMuY", &fMCMuMuY, "fMCMuMuY/F");
    fGenTree ->Branch("fMCMuMuM", &fMCMuMuM, "fMCMuMuM/F");
    // fGenTree ->Branch("fMCMuPt1", &fMCMuPt1, "fMCMuPt1/D");
    // fGenTree ->Branch("fMCMuPt2", &fMCMuPt2, "fMCMuPt2/D");
    // fGenTree ->Branch("fMCMuEta1", &fMCMuEta1, "fMCMuEta1/D");
    // fGenTree ->Branch("fMCMuEta2", &fMCMuEta2, "fMCMuEta2/D");
    // fGenTree ->Branch("fMCMuPhi1", &fMCMuPhi1, "fMCMuPhi1/D");
    // fGenTree ->Branch("fMCMuPhi2", &fMCMuPhi2, "fMCMuPhi2/D");
    // fGenTree ->Branch("fMCMuPDG1", &fMCMuPDG1, "fMCMuPDG1/D");
    // fGenTree ->Branch("fMCMuPDG2", &fMCMuPDG2, "fMCMuPDG2/D");
    // post data
  }  
  PostData(3, fGenTree);

  ////////////////////////////////////////
  //Trigger information tree
  ////////////////////////////////////////
  fTrgTree = new TTree("fTrgTree", "fTrgTree");
  if(!fIsMC){
    fTrgTree ->Branch("fTrgRunNum", &fTrgRunNum, "fTrgRunNum/I");
    fTrgTree ->Branch("fCMUP6Decision", &fCMUP6Decision, "fCMUP6Decision/I");
    fTrgTree ->Branch("fCMUP10Decision", &fCMUP10Decision, "fCMUP10Decision/I");
    fTrgTree ->Branch("fCMUP11Decision", &fCMUP11Decision, "fCMUP11Decision/I");
    // post data
  }  
  PostData(4, fTrgTree);

  ////////////////////////////////////////
  //output histograms
  ////////////////////////////////////////
  fOutputList = new TList();          // this is a list which will contain all  histograms
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  //  counter for events passing each cut    
  fCounterH = new TH1F("fCounterH", "fCounterH", 25, 0., 25.);
  fOutputList->Add(fCounterH);
  // number of positive and negative muons passing the muon selection
  fNumberMuonsH = new TH2F("fNumberMuonsH", "fNumberMuonsH", 12, 0., 12.,12, 0., 12.); 
  fOutputList->Add(fNumberMuonsH);        // don't forget to add it to the list!
  // number of positive and negative MC muons passing the muon selection
  if(fIsMC){
    fNumberMCMuonsH = new TH2F("fNumberMCMuonsH", "fNumberMCMuonsH", 12, 0., 12.,12, 0., 12.); 
    fOutputList->Add(fNumberMCMuonsH); 
  }
  // // Rabs of positive and negative muons passing the muon selection
  // fRAbsMuonH = new TH2F("fRAbsMuonH", "fRAbsMuonH", 100, 0, 100, 100, 0, 100);
  // fOutputList->Add(fRAbsMuonH);
  // // kinematics of dimouns	
  // fMuMuMassPtH = new TH2F("fMuMuMassPtH", "fMuMuMassPtH", 1500, 0, 150, 150, 0, 15);
  // fOutputList->Add(fMuMuMassPtH);

  // post data
  PostData(2, fOutputList);           
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::PostAllData()
{
  // Post data
  PostData(1, fRecTree);
  PostData(2, fOutputList);
  PostData(3, fGenTree);
  PostData(4, fTrgTree);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::TwoMuonAna(Int_t *idxPosMuons, Int_t *idxNegMuons)
{
  // Get muon masss fromn PDG
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  TParticlePDG *partMuon = pdgdat->GetParticle(13);
  Double_t MuonMass = partMuon->Mass();

  // create all four vectors
  // --  positive muon
  TLorentzVector PosMuon1;
  AliAODTrack *PosTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosMuons[0]));
  PosMuon1.SetPtEtaPhiM(PosTrack->Pt(), PosTrack->Eta(), PosTrack->Phi(), MuonMass);
  // --  negative muon
  TLorentzVector NegMuon1;
  AliAODTrack *NegTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegMuons[0]));
  NegMuon1.SetPtEtaPhiM(NegTrack->Pt(), NegTrack->Eta(), NegTrack->Phi(), MuonMass);

  // fill Rabs histo
  // fRAbsMuonH->Fill(PosTrack->GetRAtAbsorberEnd(),NegTrack->GetRAtAbsorberEnd());
  
  // fill in dimuon kinematics 
  TLorentzVector MuMu = NegMuon1+PosMuon1;
  // fMuMuMassPtH->Fill(MuMu.M(),MuMu.Pt());

  // set tree variables
  fMuMuPt = MuMu.Pt(); 
  // fMuMuPhi = MuMu.Phi();
  fMuMuY = MuMu.Rapidity(); 
  fMuMuM = MuMu.M();

  // fMuPt1 = PosTrack->Pt(); 
  // fMuEta1 = PosTrack->Eta(); 
  // fMuPhi1 = PosTrack->Phi();
  // fMuQ1 = PosTrack->Charge(); 

  // fMuPhi2 = NegTrack->Phi();
  // fMuEta2 = NegTrack->Eta();
  // fMuPt2 = NegTrack->Pt();
  // fMuQ2 = NegTrack->Charge();
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::TwoMCMuonAna(Int_t *idxMCPosMuons, Int_t *idxMCNegMuons)
{
  // get tracks  
  AliMCParticle *PosMCPart = static_cast<AliMCParticle*>(fMC->GetTrack(idxMCPosMuons[0])); 
  AliMCParticle *NegMCPart = static_cast<AliMCParticle*>(fMC->GetTrack(idxMCNegMuons[0]));

  // create lorentz vectors   
  TLorentzVector PosMCMuon,NegMCMuon,MCMuMu;

  //fill fourvectors
  PosMCMuon.SetPxPyPzE(PosMCPart->Px(), PosMCPart->Py(), PosMCPart->Pz(), PosMCPart->E());
  NegMCMuon.SetPxPyPzE(NegMCPart->Px(), NegMCPart->Py(), NegMCPart->Pz(), NegMCPart->E());
  MCMuMu = PosMCMuon+NegMCMuon;

  //connect variables
  fMCRunNum = fAOD->GetRunNumber();
  fMCMuMuM = MCMuMu.M();
  fMCMuMuPt = MCMuMu.Pt();
  fMCMuMuY = MCMuMu.Rapidity();
  // fMCMuMuPhi = MCMuMu.Phi();

  // fMCMuPt1 = PosMCPart->Pt(); 
  // fMCMuEta1 = PosMCPart->Eta(); 
  // fMCMuPhi1 = PosMCPart->Phi();
  // fMCMuPDG1 = PosMCPart->PdgCode();

  // fMCMuPhi2 = NegMCPart->Phi();
  // fMCMuEta2 = NegMCPart->Eta();
  // fMCMuPt2 = NegMCPart->Pt(); 
  // fMCMuPDG2 = NegMCPart->PdgCode();
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::UserExec(Option_t *)
{
  Int_t iSelectionCounter = 0; // no selection applied yet 
  fCounterH->Fill(iSelectionCounter); // entering UserExec 1/1 (data/MC)
  iSelectionCounter++;

  ////////////////////////////////////////////
  // Geting the AOD event
  ////////////////////////////////////////////
   // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
  if(!fAOD) {
    PostAllData();
    return;
  }                                  
  fCounterH->Fill(iSelectionCounter); // AOD event found 2/2
  iSelectionCounter++;

  ////////////////////////////////////////////
  // Selecting good runs
  ////////////////////////////////////////////
  Bool_t IsGoodRun = kFALSE;

  for(auto i : fMapRunAndLumi){
    if( fAOD->GetRunNumber() == i.first ) IsGoodRun = kTRUE;
  }
  if(!IsGoodRun) {
    PostAllData();
    return;
  }                                  
  fCounterH->Fill(iSelectionCounter); // Good run selected 3/3
  iSelectionCounter++;
  ////////////////////////////////////////////
  //  MC Luminosity scaling
  ////////////////////////////////////////////
  const char *InputName = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  Int_t scaling = 20000;

  if( strstr(InputName,"/kCohJpsiToMu/") ){
    scaling = 40000;
  } 

  if(fIsMC && fIsScalingOn && (fMapAnalysedMC[fAOD->GetRunNumber()] > Int_t(fMapRunAndLumi[fAOD->GetRunNumber()]*scaling) ) )  {
    PostAllData();
    return;
  }                 
  fMapAnalysedMC[fAOD->GetRunNumber()]++;

  ////////////////////////////////////////////
  //  MC generated particles analysis
  ////////////////////////////////////////////
  if(fIsMC){
    fGenPart->Clear("C");
    fMC = dynamic_cast<AliMCEvent*>(MCEvent()); 
    if(!fMC){
      PostAllData();
      return;
    }  
    fCounterH->Fill(iSelectionCounter); // MC generated event found -/4
    iSelectionCounter++;

    //are there particles at all?
    Int_t nMCParticles(fMC->GetNumberOfTracks()); 
    if(nMCParticles<1) {
      PostAllData();
      return;
    } 
    fCounterH->Fill(iSelectionCounter); // At least one MC generated particle -/5
    iSelectionCounter++;

    // loop over MC tracks and select muons
    Int_t nGoodMCPosMuons = 0;
    Int_t nGoodMCNegMuons = 0;  
    Int_t *idxMCPosMuons = new Int_t[nMCParticles];
    Int_t *idxMCNegMuons = new Int_t[nMCParticles];  
    for(Int_t iMCParticle = 0; iMCParticle < nMCParticles; iMCParticle++) {
      // get track
      AliMCParticle *MCPart = static_cast<AliMCParticle*>(fMC->GetTrack(iMCParticle)); 
      if(!MCPart) return;
      // Particle is primary (for LHC16b2) or is not primary and it is coming from J/Psi or Psi' decay (for LHC18l7)
      if(MCPart->GetMother() == -1){
        // if muons increase counter and store indices
        if(MCPart->PdgCode() == 13){
          idxMCPosMuons[nGoodMCPosMuons] = iMCParticle;
          nGoodMCPosMuons++;
        } else  if(MCPart->PdgCode() == -13){
          idxMCNegMuons[nGoodMCNegMuons] = iMCParticle;
          nGoodMCNegMuons++;
        }
      } else {
        AliMCParticle *MCMother = static_cast<AliMCParticle*>(fMC->GetTrack(MCPart->GetMother()));
        if(MCMother->PdgCode() != 443 && MCMother->PdgCode() != 100443) continue;
        // if muons increase counter and store indices
        if(MCPart->PdgCode() == 13){
          idxMCPosMuons[nGoodMCPosMuons] = iMCParticle;
          nGoodMCPosMuons++;
        } else  if(MCPart->PdgCode() == -13){
          idxMCNegMuons[nGoodMCNegMuons] = iMCParticle;
          nGoodMCNegMuons++;
        }
      } 
    }
    // store number of muons
    fNumberMCMuonsH->Fill(nGoodMCPosMuons,nGoodMCNegMuons);

    ////////////////////////////////////////////
    // two MC muon analysis
    ////////////////////////////////////////////
    if (!(nGoodMCPosMuons == 1 && nGoodMCNegMuons == 1)) {
      PostAllData();
      return;
    } 
    fCounterH->Fill(iSelectionCounter); // exactly one positive and one negative MC generated muons -/6
    iSelectionCounter++;
    TwoMCMuonAna(idxMCPosMuons,idxMCNegMuons);
    fGenTree->Fill();
    }
   // end of MC generated particles

  ////////////////////////////////////////////
  //  Trigger information
  ////////////////////////////////////////////
  // in 2018 q,r : CMUP6-B-NOPF-MUFAST = *0VBA 0MUL ,  
  // in 2018 q,r and 2015 o:  CMUP11-B-NOPF-MUFAST = *0VBA *0UBA *0UBC 0MUL,
  // in 2015 o : CMUP10-B-NOPF-MUFAST = = *0VBA *0UBA *0UBC 0MSL , 
  TString trigger = fAOD->GetFiredTriggerClasses();
  
  Bool_t isTriggered = kFALSE;

  if (fIsMC) {
    isTriggered = kTRUE; // No trigger required for MC
  } else {
    fTrgRunNum = fAOD->GetRunNumber();
    // ###### CMUP10+CMUP11 triggers    
    if (fTrigger.Contains("CMUP11")){
      if ((fPeriod.Contains("15o") || fPeriod.Contains("18q") || fPeriod.Contains("18r")) && trigger.Contains("CMUP11-B-NOPF-MUFAST")) {
        isTriggered = kTRUE;
        fCMUP11Decision = 1;
      } else {
        fCMUP11Decision = 0;
      }

      if ((fPeriod.Contains("15o")) && trigger.Contains("CMUP10-B-NOPF-MUFAST")) {
        isTriggered = kTRUE;
        fCMUP10Decision = 1;
      } else {
        fCMUP10Decision = 0;
      }
    }
    // ###### CMUP6 trigger  
    if (fTrigger.Contains("CMUP6")){
      if ((fPeriod.Contains("18q") || fPeriod.Contains("18r")) && trigger.Contains("CMUP6-B-NOPF-MUFAST")) {
        isTriggered = kTRUE;
        fCMUP6Decision = 1;
      } else {
        fCMUP6Decision = 0;
      }
    }
    // Fill trigger tree
    fTrgTree->Fill();
  }

  if (!isTriggered) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // right trigger found 4/7
  iSelectionCounter++;

  // get the run number and trigger inputs
  fRunNum = fAOD->GetRunNumber();
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  
  ////////////////////////////////////////////
  //  find muons
  ////////////////////////////////////////////
  //are there tracks at all?
  Int_t nTracks(fAOD->GetNumberOfTracks()); 
  if(nTracks<1) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // At least one track 5/8
  iSelectionCounter++;

  // loop over tracks and select good muons
  Int_t nGoodPosMuons = 0;
  Int_t nGoodNegMuons = 0;  
  Int_t *idxPosMuons = new Int_t[nTracks];
  Int_t *idxNegMuons = new Int_t[nTracks];  
  for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    // get track
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack)); 
    if(!track) return;

    // is it a good muon track?
    if(!track->IsMuonTrack()) continue;
    if(!fMuonTrackCuts->IsSelected(track)) continue;
    if( (track->GetRAtAbsorberEnd() < 17.5) || (track->GetRAtAbsorberEnd() > 89.5) ) continue;

    // increase counter and store indices
    if (track->Charge() > 0) {
      idxPosMuons[nGoodPosMuons] = iTrack;
      nGoodPosMuons++;
    } else  if (track->Charge() < 0) {
      idxNegMuons[nGoodNegMuons] = iTrack;
      nGoodNegMuons++;
    } 
  }
  // store number of muons
  fNumberMuonsH->Fill(nGoodPosMuons,nGoodNegMuons);

  ////////////////////////////////////////////
  // two muon analysis
  ////////////////////////////////////////////
  if (!(nGoodPosMuons == 1 && nGoodNegMuons == 1)) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // exactly one positive and one negative muons 6/9
  iSelectionCounter++;
  TwoMuonAna(idxPosMuons,idxNegMuons);

  ////////////////////////////////////////////
  // info to determine exclusivity
  ////////////////////////////////////////////
  // ---SPD
  // fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  // ---ZDC 
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // ZDC info is present 7/10
  iSelectionCounter++;

  fZNAEnergy = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy = dataZDC->GetZNCTowerEnergy()[0];
  // fZPAEnergy = dataZDC->GetZPATowerEnergy()[0];
  // fZPCEnergy = dataZDC->GetZPCTowerEnergy()[0];
  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  // for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  // for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  // at least one ZDC hit in the timing window
  fIsZNAFired = 0;
  fIsZNCFired = 0;
  for (Int_t i=0;i<4;i++){
    if ( (fZNATDC[i]>-2.) && (fZNATDC[i]<2.) ) fIsZNAFired = 1;
    if ( (fZNCTDC[i]>-2.) && (fZNCTDC[i]<2.) ) fIsZNCFired = 1;  
  }

  // ---Checks for ZDC informations
  // all ZDC hits in the timing window
  //   if (fZNATDC[iZDC] > -999.){
  //     if ((fZNATDC[iZDC]<-2.) || (fZNATDC[iZDC]>2.)) isZNAfiredAll = kFALSE;
  //   }
  //   if (fZNCTDC[iZDC] > -999.){
  //     if ((fZNCTDC[iZDC]<-2.) || (fZNCTDC[iZDC]>2.)) isZNCfiredAll = kFALSE;
  //   }    
  // }
  // // events with no hits in the ZDC
  // if ( (fZNATDC[0]<-998.5)  && (fZNATDC[1]<-998.5) && (fZNATDC[2]<-998.5) && (fZNATDC[3]<-998.5) ){
  //   isZNAfiredAll = kFALSE;
  //   fZNAEnergy0NH->Fill(fZNAEnergy);
  // }
  // if ( (fZNCTDC[0]<-998.5)  && (fZNCTDC[1]<-998.5) && (fZNCTDC[2]<-998.5) && (fZNCTDC[3]<-998.5) ){
  //   isZNCfiredAll = kFALSE;
  //  fZNCEnergy0NH->Fill(fZNCEnergy);
  // }

  // if(isZNAfired) fZNAEnergyTimingH->Fill(fZNAEnergy);  
  // if(isZNCfired) fZNCEnergyTimingH->Fill(fZNCEnergy);

  // if(isZNAfiredAll) fZNAEnergyTimingAllH->Fill(fZNAEnergy);  
  // if(isZNCfiredAll) fZNCEnergyTimingAllH->Fill(fZNCEnergy);

  // // check for timing of all hits in the triggered events
  // for(Int_t iZDC = 0; iZDC < 4 ; iZDC++) {
  //   if(isZNAfired) fZNATDCTimingH->Fill(fZNATDC[iZDC]);
  //   if(isZNCfired) fZNCTDCTimingH->Fill(fZNATDC[iZDC]);

  //   if(isZNAfiredAll) fZNATDCTimingAllH->Fill(fZNATDC[iZDC]);
  //   if(isZNCfiredAll) fZNCTDCTimingAllH->Fill(fZNATDC[iZDC]);
  // }

  // ---V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); //  V0 info 8/11
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  Int_t nV0CFiredCells = 0;
  Int_t nV0AFiredCells = 0;

  for(Int_t i = 0; i < 64; i++) {
    if(dataVZERO->GetBBFlag(i) == kTRUE) {
      if(i < 32) {
        nV0CFiredCells += 1;
      } else {
        nV0AFiredCells += 1;
      }
    }
  }

  fV0CFiredCells = nV0CFiredCells;
  fV0AFiredCells = nV0AFiredCells;

  // ---AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(!dataAD){
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); //  AD info 9/12
  iSelectionCounter++;

  fADADecision = dataAD->GetADADecision();
  fADCDecision = dataAD->GetADCDecision();

  // //Past-future protection maps
  // fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  // fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

  // fill the tree
  fRecTree->Fill();

  // post the data
  PostAllData();

  // clean up
  delete [] idxPosMuons;
  delete [] idxNegMuons;

}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoMUON::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
// ----------------------------------------------------------------------------------------------------------------------------------


