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

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TList.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"

#include "AliAnalysisTaskPWG4PidDetEx.h"

// STL includes
#include <iostream>
using namespace std;

//
// Analysis class example for PID using detector signals.
// Works on ESD and AOD; has a flag for MC analysis 
//
//    Alexandru.Dobrin@hep.lu.se
//    Peter.Christiansen@hep.lu.se
// 

ClassImp(AliAnalysisTaskPWG4PidDetEx)

const Double_t AliAnalysisTaskPWG4PidDetEx::fgkCtau = 370;                //  distance for kaon decay
const Double_t AliAnalysisTaskPWG4PidDetEx::fgkPionMass = 1.39570000000000000e-01;
const Double_t AliAnalysisTaskPWG4PidDetEx::fgkKaonMass = 4.93599999999999983e-01; 
const Double_t AliAnalysisTaskPWG4PidDetEx::fgkProtonMass = 9.38270000000000048e-01;
const Double_t AliAnalysisTaskPWG4PidDetEx::fgkC = 2.99792458e-2;





//_____________________________________________________________________________
AliAnalysisTaskPWG4PidDetEx::AliAnalysisTaskPWG4PidDetEx():
  AliAnalysisTaskSE(),
  fESD(0),
  fAOD(0), 
  fListOfHists(0), 
  fEtaCut(0.9), 
  fPtCut(0.5), 
  fXbins(20), 
  fXmin(0.),  
  fTOFCutP(2.), 
  fTrackFilter(0x0),
  fAnalysisMC(kTRUE),
  fAnalysisType("ESD"),
  fTriggerMode(kMB2),
  fEvents(0), fEffTot(0), fEffPID(0), fAccP(0), fAccPt(0), fKaonDecayCorr(0),
  fdNdPt(0), fMassAll(0),
  fdNdPtPion(0), fMassPion(0), fdEdxTPCPion(0), fbgTPCPion(0),
  fdNdPtKaon(0), fMassKaon(0), fdEdxTPCKaon(0), fbgTPCKaon(0),
  fdNdPtProton(0), fMassProton(0), fdEdxTPCProton(0), fbgTPCProton(0),
  fdNdPtMC(0), fdNdPtMCPion(0), fdNdPtMCKaon(0), fdNdPtMCProton(0) 
{
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskPWG4PidDetEx::AliAnalysisTaskPWG4PidDetEx(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0),
  fAOD(0), 
  fListOfHists(0), 
  fEtaCut(0.9), 
  fPtCut(0.5), 
  fXbins(20), 
  fXmin(0.), 
  fTOFCutP(2.), 
  fTrackFilter(0x0),
  fAnalysisMC(kTRUE),
  fAnalysisType("ESD"),
  fTriggerMode(kMB2),
  fEvents(0), fEffTot(0), fEffPID(0), fAccP(0), fAccPt(0), fKaonDecayCorr(0),
  fdNdPt(0), fMassAll(0),
  fdNdPtPion(0), fMassPion(0), fdEdxTPCPion(0), fbgTPCPion(0),
  fdNdPtKaon(0), fMassKaon(0), fdEdxTPCKaon(0), fbgTPCKaon(0),
  fdNdPtProton(0), fMassProton(0), fdEdxTPCProton(0), fbgTPCProton(0),
  fdNdPtMC(0), fdNdPtMCPion(0), fdNdPtMCKaon(0), fdNdPtMCProton(0) 
{
  // Output slot #0 writes into a TList
  DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliAnalysisTaskPWG4PidDetEx::~AliAnalysisTaskPWG4PidDetEx()
{
  // Destructor
}

// //______________________________________________________________________________
// void AliAnalysisTaskPWG4PidDetEx::ConnectInputData(Option_t *)
// {
//   // Connect AOD here
//   // Called once
//   if (fDebug > 1)  AliInfo("ConnectInputData() \n");
  
//   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
//   if (!tree) {
//     Printf("ERROR: Could not read chain from input slot 0");
//   } 
//   else {
//     if(fAnalysisType == "ESD") {
//       AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());    
//       if (!esdH) {
// 	Printf("ERROR: Could not get ESDInputHandler");
//       } else
// 	fESD = esdH->GetEvent();
//     }
//     else if(fAnalysisType == "AOD") {
//       AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
//       if (!aodH) {
// 	Printf("ERROR: Could not get AODInputHandler");
//       } else
// 	fAOD = aodH->GetEvent();
//     }
//   }

// }

//______________________________________________________________________________
void AliAnalysisTaskPWG4PidDetEx::UserCreateOutputObjects()
{ 
  OpenFile(1);
  fListOfHists = new TList();

  fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 1, 0, 1);
  fListOfHists->Add(fEvents);

  fEffTot = new TH1F ("fEffTot","Efficiency; p_{T} [GeV/c]; Counts", fXbins, fXmin, fTOFCutP);
  fEffTot->Sumw2();
  fListOfHists->Add(fEffTot);

  fEffPID = new TH1F ("fEffPID","Efficiency; p_{T} [GeV/c]; Counts", fXbins, fXmin, fTOFCutP);
  fEffPID->Sumw2();
  fListOfHists->Add(fEffPID);

  fAccP = new TH1F ("fAccP","Acceptance; p_{T} [GeV/c]; Counts", fXbins, fXmin, fTOFCutP);
  fAccP->Sumw2();
  fListOfHists->Add(fAccP);

  fAccPt = new TH1F ("fAccPt","Acceptance; p_{T} [GeV/c]; Counts", fXbins, fXmin, fTOFCutP);
  fAccPt->Sumw2();
  fListOfHists->Add(fAccPt);

  fKaonDecayCorr = new TProfile("fKaonDecayCorr","Kaon decay correction vs p_{T}; p_{T} [GeV/c]; Kaon decay correction", fXbins, fXmin, fTOFCutP);
  fListOfHists->Add(fKaonDecayCorr);

  fdNdPt = new TH1F("fdNdPt","p_{T} distribution; p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", 100 , 0, 10);
  fdNdPt->Sumw2();
  fListOfHists->Add(fdNdPt);

  fMassAll = new TH2F("fMassAll","Mass^{2} vs momentum (ToF); p [GeV/c]; M^{2} [GeV^{2}/c^{4}]", 100, -5, 5, 100, -0.2, 1.2);
  fListOfHists->Add(fMassAll);

  //Pion
  fdNdPtPion = new TH1F("fdNdPtPion","p_{T} distribution (Pions); p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
  fdNdPtPion->Sumw2();
  fListOfHists->Add(fdNdPtPion);

  fMassPion = new TH2F("fMassPion","Mass squared for pions after cuts vs pmean; p [GeV/c]; Mass^{2} [GeV^{2}/c^{4}]", 100, -5, 5, 100, -0.2, 1.2);
  fListOfHists->Add(fMassPion);

  fdEdxTPCPion = new TH2F("fdEdxTPCPion","Energy loss vs momentum; p [GeV/c]; dE/dx [a. u.]", 40, -2, 2, 100, 0, 200);
  fListOfHists->Add(fdEdxTPCPion);

  fbgTPCPion = new TH2F("fbgTPCPion", "Energy loss vs #beta#gamma (Pions);#beta#gamma; dE/dx [a. u.]", 100, 0, 15, 100, 0, 200);
  fListOfHists->Add(fbgTPCPion);

  //Kaon
  fdNdPtKaon = new TH1F("fdNdPtKaon","p_{T} distribution (Kaons);p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
  fdNdPtKaon->Sumw2();
  fListOfHists->Add(fdNdPtKaon);

  fMassKaon = new TH2F("fMassKaon","Mass squared for kaons after cuts vs pmean; p [GeV/c]; Mass^{2} [GeV^{2}/c^{4}]", 100, -5, 5, 100, -0.2, 1.2);
  fListOfHists->Add(fMassKaon);

  fdEdxTPCKaon = new TH2F("fdEdxTPCKaon","Energy loss vs momentum; p [GeV/c]; dE/dx [a. u.]", 40, -2, 2, 100, 0, 200);
  fListOfHists->Add(fdEdxTPCKaon);

  fbgTPCKaon = new TH2F("fbgTPCKaon", "Energy loss vs #beta#gamma (Kaons);#beta#gamma; dE/dx [a. u.]", 100, 0, 15, 100, 0, 200);
  fListOfHists->Add(fbgTPCKaon);

  //Proton
  fdNdPtProton = new TH1F("fdNdPtProton","p_{T} distribution (Protons);p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
  fdNdPtProton->Sumw2();
  fListOfHists->Add(fdNdPtProton); 
  
  fMassProton = new TH2F("fMassProton","Mass squared for protons after cuts vs pmean; p [GeV/c]; Mass^{2} [GeV^{2}/c^{4}]", 100, -5, 5, 100, -0.2, 1.2);
  fListOfHists->Add(fMassProton);

  fdEdxTPCProton = new TH2F("fdEdxTPCProton","Energy loss vs momentum; p [GeV/c]; dE/dx [a. u.]", 40, -2, 2, 100, 0, 200);
  fListOfHists->Add(fdEdxTPCProton);

  fbgTPCProton = new TH2F("fbgTPCProton", "Energy loss vs #beta#gamma (Protons);#beta#gamma; dE/dx [a. u.]", 100, 0, 15, 100, 0, 200);
  fListOfHists->Add(fbgTPCProton);


  //MC
  if(fAnalysisMC){
    fdNdPtMC = new TH1F("fdNdPtMC","p_{T} distribution;p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", 100, 0, 10);
    fdNdPtMC->SetLineColor(2);
    fdNdPtMC->Sumw2();
    fListOfHists->Add(fdNdPtMC);

    fdNdPtMCPion = new TH1F("fdNdPtMCPion","p_{T} distribution;p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
    fdNdPtMCPion->SetLineColor(2);
    fdNdPtMCPion->Sumw2();
    fListOfHists->Add(fdNdPtMCPion);
  
    fdNdPtMCKaon = new TH1F("fdNdPtMCKaon","p_{T} distribution;p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
    fdNdPtMCKaon->SetLineColor(2);
    fdNdPtMCKaon->Sumw2();
    fListOfHists->Add(fdNdPtMCKaon);

    fdNdPtMCProton = new TH1F("fdNdPtMCProton","p_{T} distribution;p_{T} [GeV/c]; #frac{1}{2#pip_{T}} #frac{1}{N_{ev}} #frac{dN}{dp_{T}}", fXbins, fXmin, fTOFCutP);
    fdNdPtMCProton->SetLineColor(2);
    fdNdPtMCProton->Sumw2();
    fListOfHists->Add(fdNdPtMCProton);
  }

}

//______________________________________________________________________________
void AliAnalysisTaskPWG4PidDetEx::UserExec(Option_t *) 
{

  // Create histograms
  // Called once
  if (fAnalysisType == "AOD") {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      return;
    }
    else{
      //  assume that the AOD is in the general output...
      //  fAOD  = AODEvent();
      //  if(!fAOD){
      //	Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
    }
  }
  else if  (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
  }

  // Main loop
  // Called for each event
  if (fDebug > 1) AliInfo("Exec() \n" );

  if(fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    if(IsEventTriggered(fESD, fTriggerMode)) {
      Printf("PWG4 PID ESD analysis");
      AnalyzeESD(fESD);
    }//triggered event
  }//ESD analysis  
              
  else if(fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    } 
    if(IsEventTriggered(fAOD, fTriggerMode)) {
      Printf("PWG4 PID AOD analysis");
      AnalyzeAOD(fAOD);
    }//triggered event
  }//AOD analysis  

  // Post output data.
  PostData(1, fListOfHists);
}

//________________________________________________________________________
void AliAnalysisTaskPWG4PidDetEx::AnalyzeESD(AliESDEvent* esd)
{
  // Get vertex 
  const AliESDVertex* primaryVertex = esd->GetPrimaryVertex();   
  const Double_t vertexResZ = primaryVertex->GetZRes();

  // Select only events with a reconstructed vertex
  if(vertexResZ<5.0) {
    const Int_t nESDTracks = esd->GetNumberOfTracks();
    for(Int_t i = 0; i < nESDTracks; i++) {
      AliESDtrack* trackESD = esd->GetTrack(i);

      //Cuts
      UInt_t selectInfo = 0;
      if (fTrackFilter) {
	selectInfo = fTrackFilter->IsSelected(trackESD);
	if (!selectInfo) continue;
      }

      if ((trackESD->Pt() < fPtCut) || (TMath::Abs(trackESD->Eta()) > fEtaCut ))
	continue;

      //Corrections
      fEffTot->Fill(trackESD->Pt());
      if (trackESD->GetTOFsignal() !=0)
	fEffPID->Fill(trackESD->Pt());

      if (trackESD->P() < fTOFCutP) fAccP->Fill(trackESD->Pt());
      if (trackESD->Pt() < fTOFCutP) fAccPt->Fill(trackESD->Pt());

      //Analysis
      fdNdPt->Fill(trackESD->Pt(), 1.0/trackESD->Pt());

      //TOF 
      if ((trackESD->GetIntegratedLength() == 0) || (trackESD->GetTOFsignal() == 0))
	continue;

      fMassAll->Fill(trackESD->Charge()*trackESD->P(), MassSquared(trackESD));

      if ((MassSquared(trackESD) < 0.15) && (MassSquared(trackESD) > -0.15) && (trackESD->P() < fTOFCutP)){
	fdNdPtPion->Fill(trackESD->Pt(), 1.0/trackESD->Pt());
	fMassPion->Fill(trackESD->Charge()*trackESD->P(), MassSquared(trackESD));
	fdEdxTPCPion->Fill(trackESD->Charge()*trackESD->P(),trackESD->GetTPCsignal());
	fbgTPCPion->Fill(trackESD->P()/fgkPionMass,trackESD->GetTPCsignal());
      }

      if ((MassSquared(trackESD) > 0.15) && (MassSquared(trackESD) < 0.45) && (trackESD->P() < fTOFCutP)){
	fdNdPtKaon->Fill(trackESD->Pt(), 1.0/trackESD->Pt());
	fMassKaon->Fill(trackESD->Charge()*trackESD->P(),  MassSquared(trackESD));
	fdEdxTPCKaon->Fill(trackESD->Charge()*trackESD->P(), trackESD->GetTPCsignal());
	fbgTPCKaon->Fill(trackESD->P()/fgkKaonMass, trackESD->GetTPCsignal());
	//Kaon decay correction
	fKaonDecayCorr->Fill(trackESD->Pt(), TMath::Exp(KaonDecay(trackESD)));
      }

      if ((MassSquared(trackESD) > 0.75) && (MassSquared(trackESD) < 1.05) && (trackESD->P() < fTOFCutP)){
	fdNdPtProton->Fill(trackESD->Pt(), 1.0/trackESD->Pt());
	fMassProton->Fill(trackESD->Charge()*trackESD->P(), MassSquared(trackESD));
	fdEdxTPCProton->Fill(trackESD->Charge()*trackESD->P(),trackESD->GetTPCsignal());
	fbgTPCProton->Fill(trackESD->P()/fgkProtonMass,trackESD->GetTPCsignal());
      }
    }//ESD track loop

    //MC
    if(fAnalysisMC){      
      AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!mcHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
	return;
      }

      AliMCEvent* mcEvent = mcHandler->MCEvent();
      if (!mcEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }

      AliStack* mcStack = mcEvent->Stack();
      if (!mcStack) {
	Printf("ERROR: Could not retrieve MC stack");
	return;
      }

      const Int_t nTracksMC = mcStack->GetNtrack();      
      for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
	//Cuts
	if(!(mcStack->IsPhysicalPrimary(iTracks)))
	  continue;

	TParticle* mcTrack = mcStack->Particle(iTracks);
     
	Double_t charge = mcTrack->GetPDG()->Charge();
	if (charge == 0)
	  continue;
	
	if ((mcTrack->Pt() < fPtCut) || (TMath::Abs(mcTrack->Eta()) > fEtaCut ))
	  continue;

	//Analysis
	fdNdPtMC->Fill(mcTrack->Pt(), 1.0/mcTrack->Pt());

	if ((mcTrack->GetPdgCode() == 211) || (mcTrack->GetPdgCode() == -211))
	  fdNdPtMCPion->Fill(mcTrack->Pt(),1.0/mcTrack->Pt());

	if ((mcTrack->GetPdgCode() == 321) || (mcTrack->GetPdgCode() == -321))
	  fdNdPtMCKaon->Fill(mcTrack->Pt(),1.0/mcTrack->Pt());

	if ((mcTrack->GetPdgCode() == 2212) || (mcTrack->GetPdgCode() == -2212))
	  fdNdPtMCProton->Fill(mcTrack->Pt(),1.0/mcTrack->Pt());
      }//MC track loop 
    }//if MC
    fEvents->Fill(0);
  }//if vertex

}

//________________________________________________________________________
void AliAnalysisTaskPWG4PidDetEx::AnalyzeAOD(AliAODEvent* aod)
{   
  // Get vertex 
  AliAODVertex* primaryVertex = aod->GetPrimaryVertex();    
  const Double_t vertexResZ = TMath::Sqrt(primaryVertex->RotatedCovMatrixZZ());

  // Select only events with a reconstructed vertex
  if (vertexResZ < 5) {  
    //AOD
    Int_t nGoodTracks = aod->GetNumberOfTracks();

    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {     
      AliAODTrack* trackAOD = aod->GetTrack(iTracks);
      //Cuts
      if (!(trackAOD->IsPrimaryCandidate()))
	continue;

      if ((trackAOD->Pt() < fPtCut) || (TMath::Abs(trackAOD->Eta()) > fEtaCut ))
	continue;

      //Corrections
      fEffTot->Fill(trackAOD->Pt());
      if (trackAOD->GetDetPid()->GetTOFsignal() !=0 )
	fEffPID->Fill(trackAOD->Pt());

      if (trackAOD->P() < fTOFCutP) fAccP->Fill(trackAOD->Pt());
      if (trackAOD->Pt() < fTOFCutP) fAccPt->Fill(trackAOD->Pt());

      //Analysis
      fdNdPt->Fill(trackAOD->Pt(), 1.0/trackAOD->Pt());

      //TOF
      if ((IntegratedLength(trackAOD) == 0) || (trackAOD->GetDetPid()->GetTOFsignal() == 0))
	continue;

      fMassAll->Fill(trackAOD->Charge()*trackAOD->P(), MassSquared(trackAOD));

      if ((MassSquared(trackAOD) < 0.15) && (MassSquared(trackAOD) > -0.15) && (trackAOD->P() < fTOFCutP)){
	fdNdPtPion->Fill(trackAOD->Pt(), 1.0/trackAOD->Pt());
	fMassPion->Fill(trackAOD->Charge()*trackAOD->P(), MassSquared(trackAOD));
	fdEdxTPCPion->Fill(trackAOD->Charge()*trackAOD->P(),trackAOD->GetDetPid()->GetTPCsignal());
	fbgTPCPion->Fill(trackAOD->P()/fgkPionMass,trackAOD->GetDetPid()->GetTPCsignal());
      }

      if ((MassSquared(trackAOD) > 0.15) && (MassSquared(trackAOD) < 0.45) && (trackAOD->P() < fTOFCutP)){
	fdNdPtKaon->Fill(trackAOD->Pt(), 1.0/trackAOD->Pt());
	fMassKaon->Fill(trackAOD->Charge()*trackAOD->P(), MassSquared(trackAOD));
	fdEdxTPCKaon->Fill(trackAOD->Charge()*trackAOD->P(),trackAOD->GetDetPid()->GetTPCsignal());
	fbgTPCKaon->Fill(trackAOD->P()/fgkKaonMass,trackAOD->GetDetPid()->GetTPCsignal());
	//Kaon decay correction
	fKaonDecayCorr->Fill(trackAOD->Pt(), TMath::Exp(KaonDecay(trackAOD)));
      }

      if ((MassSquared(trackAOD) > 0.75) && (MassSquared(trackAOD) < 1.05) && (trackAOD->P() < fTOFCutP)){
	fdNdPtProton->Fill(trackAOD->Pt(), 1.0/trackAOD->Pt());
	fMassProton->Fill(trackAOD->Charge()*trackAOD->P(), MassSquared(trackAOD));
	fdEdxTPCProton->Fill(trackAOD->Charge()*trackAOD->P(),trackAOD->GetDetPid()->GetTPCsignal());
	fbgTPCProton->Fill(trackAOD->P()/fgkProtonMass,trackAOD->GetDetPid()->GetTPCsignal());
      }
    }//AOD track loop

    //MC
    if(fAnalysisMC){
      TClonesArray* farray = (TClonesArray*)aod->FindListObject("mcparticles");
      Int_t ntrks = farray->GetEntries();
      for(Int_t i =0 ; i < ntrks; i++){   
	AliAODMCParticle* trk = (AliAODMCParticle*)farray->At(i);
	//Cuts
	if (!(trk->IsPhysicalPrimary()))
	  continue;

	if (trk->Charge() == 0)
	  continue;

	if ((trk->Pt() < fPtCut) || (TMath::Abs(trk->Eta()) > fEtaCut ))
	  continue;
	
	//Analysis
	fdNdPtMC->Fill(trk->Pt(), 1.0/trk->Pt());

	if ((trk->GetPdgCode() == 211) || (trk->GetPdgCode() == -211))
	  fdNdPtMCPion->Fill(trk->Pt(),1.0/trk->Pt());

	if ((trk->GetPdgCode() == 321) || (trk->GetPdgCode() == -321))
	  fdNdPtMCKaon->Fill(trk->Pt(),1.0/trk->Pt());

	if ((trk->GetPdgCode() == 2212) || (trk->GetPdgCode() == -2212))
	  fdNdPtMCProton->Fill(trk->Pt(),1.0/trk->Pt());
      }//MC track loop 
    }//if MC 
    fEvents->Fill(0);   
  }//if vertex resolution

}

//________________________________________________________________________
Bool_t AliAnalysisTaskPWG4PidDetEx::IsEventTriggered(AliVEvent* ev, TriggerMode trigger) 
{
  //adapted from PWG2 (AliAnalysisTaskProton)

  ULong64_t triggerMask = ev->GetTriggerMask();
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  switch (trigger) {
  case kMB1: 
    if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  
  case kMB2: 
    if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  
  case kSPDFASTOR: 
    if (triggerMask & spdFO)
      return kTRUE;
    break;
  
  }//switch

  return kFALSE;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskPWG4PidDetEx::IntegratedLength(AliVTrack* track) const
{
  Double_t intTime [5];
  for (Int_t i = 0; i < 5; i++) intTime[i] = -100.;
  Double_t timeElectron = 0, intLength = 0;

  AliAODTrack* trackAOD = (AliAODTrack*)track;
  trackAOD->GetDetPid()->GetIntegratedTimes(intTime);
  timeElectron = intTime[0];
  intLength = fgkC*timeElectron;

  return intLength;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskPWG4PidDetEx::MassSquared(AliVTrack* track) const
{
  Double_t beta = -10, mass = -10;

  if(fAnalysisType == "ESD"){
    AliESDtrack* trackESD = (AliESDtrack*)track;
    beta = trackESD->GetIntegratedLength()/trackESD->GetTOFsignal()/fgkC;
    mass = trackESD->P()*trackESD->P()*(1./(beta*beta) - 1.0);
  }

  if(fAnalysisType == "AOD"){
    AliAODTrack* trackAOD = (AliAODTrack*)track;
    beta =IntegratedLength(trackAOD)/trackAOD->GetDetPid()->GetTOFsignal()/fgkC;
    mass = trackAOD->P()*trackAOD->P()*(1./(beta*beta) - 1.0);
  }
  
  return mass;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskPWG4PidDetEx::KaonDecay(AliVTrack* track) const
{
  Double_t decay = -10;

  if(fAnalysisType == "ESD"){
    AliESDtrack* trackESD = (AliESDtrack*)track;
    decay = trackESD->GetIntegratedLength()*fgkKaonMass/fgkCtau/trackESD->P();
  }

  if (fAnalysisType == "AOD"){
    AliAODTrack* trackAOD = (AliAODTrack*)track;
    decay = IntegratedLength(trackAOD)*fgkKaonMass/fgkCtau/trackAOD->P();
  }

  return decay;
}

//_____________________________________________________________________________
void AliAnalysisTaskPWG4PidDetEx::Terminate(Option_t *)
{ 
  // Terminate loop
  if (fDebug > 1) Printf("Terminate()");
 
  //
  // The followig code has now been moved to drawPid.C
  //

//   fListOfHists = dynamic_cast<TList*> (GetOutputData(0));
//   if (!fListOfHists) {
//     Printf("ERROR: fListOfHists not available");
//     return;
//   }
 
//   TH1I* hevents = dynamic_cast<TH1I*> (fListOfHists->FindObject("fEvents"));
//   Int_t nEvents = (Int_t) hevents->GetBinContent(1);

//   const Float_t normalization = 2.0*fEtaCut*2.0*TMath::Pi();

//   //-----------------
//   TCanvas* c1 = new TCanvas("c1", "c1");
//   c1->cd();

//   TH2F* hMassAll =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fMassAll"));
//   hMassAll->SetLineColor(1);
//   hMassAll->SetMarkerColor(1);
//   hMassAll->DrawCopy();

//   TH2F* hMassPion =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fMassPion"));
//   hMassPion->SetLineColor(2);
//   hMassPion->SetMarkerColor(2);
//   hMassPion->SetMarkerStyle(7);
//   hMassPion->DrawCopy("same");

//   TH2F* hMassKaon =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fMassKaon"));
//   hMassKaon->SetLineColor(3);
//   hMassKaon->SetMarkerColor(3);
//   hMassKaon->SetMarkerStyle(7);
//   hMassKaon->DrawCopy("same");

//   TH2F* hMassProton =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fMassProton"));
//   hMassProton->SetLineColor(4);
//   hMassProton->SetMarkerColor(4);
//   hMassProton->SetMarkerStyle(7);
//   hMassProton->DrawCopy("same");

//   TLegend* legend1 = new TLegend(0.8, 0.8, 1, 1);    
//   legend1->SetBorderSize(0);
//   legend1->SetFillColor(0);
//   legend1->AddEntry(hMassAll, "All","LP");
//   legend1->AddEntry(hMassPion, "Pions","LP");
//   legend1->AddEntry(hMassKaon, "Kaons","LP");
//   legend1->AddEntry(hMassProton, "Protons","LP");  
//   legend1->Draw();

//   c1->Update();

//   //-----------------
//   TCanvas* c2 = new TCanvas("c2", "c2");
//   c2->cd();

//   TH2F* hdEdxTPCPion =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fdEdxTPCPion"));
//   hdEdxTPCPion->SetTitle("dE/dx vs p (TPC)");
//   hdEdxTPCPion->SetLineColor(2);
//   hdEdxTPCPion->SetMarkerColor(2);
//   hdEdxTPCPion->SetMarkerStyle(7);
//   hdEdxTPCPion->DrawCopy();

//   TH2F* hdEdxTPCKaon =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fdEdxTPCKaon"));
//   hdEdxTPCKaon->SetLineColor(3);
//   hdEdxTPCKaon->SetMarkerColor(3);
//   hdEdxTPCKaon->SetMarkerStyle(7);
//   hdEdxTPCKaon->DrawCopy("same");

//   TH2F* hdEdxTPCProton =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fdEdxTPCProton"));
//   hdEdxTPCProton->SetLineColor(4);
//   hdEdxTPCProton->SetMarkerColor(4);
//   hdEdxTPCProton->SetMarkerStyle(7);
//   hdEdxTPCProton->DrawCopy("same");

//   TLegend* legend2 = new TLegend(0.66, 0.66, 0.88, 0.88);    
//   legend2->SetBorderSize(0);
//   legend2->SetFillColor(0);
//   legend2->AddEntry(hdEdxTPCPion, "Pions","LP");
//   legend2->AddEntry(hdEdxTPCKaon, "Kaons","LP");
//   legend2->AddEntry(hdEdxTPCProton, "Protons","LP");
//   legend2->Draw();

//   c2->Update();

//   //-----------------
//   TCanvas* c3 = new TCanvas("c3", "c3");
//   c3->cd();

//   TH2F* hdEdxTPCbgPion =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fbgTPCPion"));
//   hdEdxTPCbgPion->SetTitle("dE/dx vs #beta#gamma (TPC)");
//   hdEdxTPCbgPion->SetLineColor(2);
//   hdEdxTPCbgPion->SetMarkerColor(2);
//   hdEdxTPCbgPion->SetMarkerStyle(7);
//   hdEdxTPCbgPion->DrawCopy();

//   TH2F* hdEdxTPCbgKaon =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fbgTPCKaon"));
//   hdEdxTPCbgKaon->SetLineColor(3);
//   hdEdxTPCbgKaon->SetMarkerColor(3);
//   hdEdxTPCbgKaon->SetMarkerStyle(7);
//   hdEdxTPCbgKaon->DrawCopy("same");

//   TH2F* hdEdxTPCbgProton =  dynamic_cast<TH2F*> (fListOfHists->FindObject("fbgTPCProton"));
//   hdEdxTPCbgProton->SetLineColor(4);
//   hdEdxTPCbgProton->SetMarkerColor(4);
//   hdEdxTPCbgProton->SetMarkerStyle(7);
//   hdEdxTPCbgProton->DrawCopy("same");

//   TLegend* legend3 = new TLegend(0.66, 0.66, 0.88, 0.88);    
//   legend3->SetBorderSize(0);
//   legend3->SetFillColor(0);
//   legend3->AddEntry(hdEdxTPCbgPion, "Pions","LP");
//   legend3->AddEntry(hdEdxTPCbgKaon, "Kaons","LP");
//   legend3->AddEntry(hdEdxTPCbgProton, "Protons","LP");
//   legend3->Draw();

//   c3->Update();

//   //-----------------
//   TCanvas* c4 = new TCanvas("c4", "c4", 100, 100, 500, 900);
//   c4->Divide(1,3);

//   c4->cd(1);
//   TH1F* hAccepatncePt = dynamic_cast<TH1F*> (fListOfHists->FindObject("fAccPt"));
//   TH1F* hAcceptanceP = dynamic_cast<TH1F*> (fListOfHists->FindObject("fAccP"));
//   hAccepatncePt->Divide(hAcceptanceP);
//   hAccepatncePt->SetTitle("Acceptance correction");
//   hAccepatncePt->GetYaxis()->SetTitle("Acceptance correction");
//   hAccepatncePt->DrawCopy();

//   c4->cd(2);
//   TH1F* hEfficiencyPID = dynamic_cast<TH1F*> (fListOfHists->FindObject("fEffPID"));
//   TH1F* hEfficiencyTot = dynamic_cast<TH1F*> (fListOfHists->FindObject("fEffTot"));
//   hEfficiencyPID->Divide(hEfficiencyTot);
//   hEfficiencyPID->SetTitle("Efficiency correction");
//   hEfficiencyPID->GetYaxis()->SetTitle("Efficiency correction");
//   hEfficiencyPID->DrawCopy();

//   c4->cd(3);
//   TProfile* hKDecayCorr = dynamic_cast<TProfile*> (fListOfHists->FindObject("fKaonDecayCorr"));
//   hKDecayCorr->SetTitle("Kaon decay correction");
//   hKDecayCorr->GetYaxis()->SetTitle("Kaon decay correction");
//   hKDecayCorr->GetXaxis()->SetTitle(" p_{T} [GeV/c]");
//   hKDecayCorr->DrawCopy();

//   c4->Update();

//   //---------------------
//   TCanvas* c5 = new TCanvas("c5", "c5", 100, 100, 600, 900);
//   c5->Divide(1,2);

//   c5->cd(1);
//   TH1F* hPtPionNoCorr = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtPion"));
//   hPtPionNoCorr->Scale(1.0/nEvents/normalization/hPtPionNoCorr->GetBinWidth(1));
//   hPtPionNoCorr->SetTitle("p_{T} distribution (no corrections)");
//   hPtPionNoCorr->SetLineColor(2);
//   hPtPionNoCorr->GetYaxis()->SetRangeUser(1E-3,1E1);
//   hPtPionNoCorr->DrawCopy();

//   TH1F* hPtKaonNoCorr = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtKaon"));
//   hPtKaonNoCorr->Scale(1.0/nEvents/normalization/hPtKaonNoCorr->GetBinWidth(1));
//   hPtKaonNoCorr->SetLineColor(3);
//   hPtKaonNoCorr->GetYaxis()->SetRangeUser(1E-3,1E1);
//   hPtKaonNoCorr->DrawCopy("same");

//   TH1F* hPtProtonNoCorr = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtProton"));
//   hPtProtonNoCorr->Scale(1.0/nEvents/normalization/hPtProtonNoCorr->GetBinWidth(1));
//   hPtProtonNoCorr->SetLineColor(4);
//   hPtProtonNoCorr->GetYaxis()->SetRangeUser(1E-3,1E1);
//   hPtProtonNoCorr->DrawCopy("same");

//   TH1F* hPt = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPt"));
//   hPt->Scale(1.0/nEvents/normalization/hPt->GetBinWidth(1));
//   hPt->GetYaxis()->SetRangeUser(1E-3,1E1);
//   hPt->DrawCopy("same");

//   TLegend* legend4 = new TLegend(0.63, 0.63, 0.88, 0.88);    
//   legend4->SetBorderSize(0);
//   legend4->SetFillColor(0);
//   legend4->AddEntry(hPt, "p_{T} dist (all)","L");
//   legend4->AddEntry(hPtPionNoCorr, "p_{T} dist (Pions)","L");
//   legend4->AddEntry(hPtKaonNoCorr, "p_{T} dist (Kaons)","L");
//   legend4->AddEntry(hPtProtonNoCorr, "p_{T} dist (Protons)","L");
//   legend4->Draw();

//   gPad->SetLogy();


//   c5->cd(2);
//   TH1F* hPtPionCorr = static_cast<TH1F*>(hPtPionNoCorr->Clone());
//   hPtPionCorr->SetTitle("p_{T} distribution (with corrections)");
//   hPtPionCorr->Multiply(hAccepatncePt);
//   hPtPionCorr->Divide(hEfficiencyPID);
//   hPtPionCorr->GetYaxis()->SetRangeUser(1E-2,1E1);
//   hPtPionCorr->DrawCopy();

//   TH1F* hPtKaonCorr = static_cast<TH1F*>(hPtKaonNoCorr->Clone());
//   hPtKaonCorr->Multiply(hAccepatncePt);
//   hPtKaonCorr->Divide(hEfficiencyPID);
//   hPtKaonCorr->Multiply(hKDecayCorr);
//   hPtKaonCorr->GetYaxis()->SetRangeUser(1E-2,1E1);
//   hPtKaonCorr->DrawCopy("same");

//   TH1F* hPtProtonCorr = static_cast<TH1F*>(hPtProtonNoCorr->Clone());
//   hPtProtonCorr->Multiply(hAccepatncePt);
//   hPtProtonCorr->Divide(hEfficiencyPID);
//   hPtProtonCorr->GetYaxis()->SetRangeUser(1E-2,1E1);
//   hPtProtonCorr->DrawCopy("same");

//   hPt->GetYaxis()->SetRangeUser(1E-2,1E1);
//   hPt->DrawCopy("same");

//   TLegend* legend5 = new TLegend(0.63, 0.63, 0.88, 0.88);    
//   legend5->SetBorderSize(0);
//   legend5->SetFillColor(0);
//   legend5->AddEntry(hPt, "p_{T} dist (all)","L");
//   legend5->AddEntry(hPtPionCorr, "p_{T} dist (Pions)","L");
//   legend5->AddEntry(hPtKaonCorr, "p_{T} dist (Kaons)","L");
//   legend5->AddEntry(hPtProtonCorr, "p_{T} dist (Protons)","L");
//   legend5->Draw();

//   gPad->SetLogy();

//   c5->Update();

//   //-----------------
//   if (fAnalysisMC){
//     TCanvas* c6 = new TCanvas("c6", "c6", 100, 100, 1200, 800); 
//     c6->Divide(2,2);

//     c6->cd(1);
//     TH1F* hPt_clone = static_cast<TH1F*>(hPt->Clone()); 
//     hPt_clone->GetYaxis()->SetRangeUser(1E-5,1E1);
//     hPt_clone->SetTitle("p_{T} distribution (all)");
//     hPt_clone->SetLineColor(1);
//     hPt_clone->DrawCopy();
 
//     TH1F* hPtMC = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtMC"));
//     hPtMC->Scale(1.0/nEvents/normalization/hPtMC->GetBinWidth(1));
//     hPtMC->GetYaxis()->SetRangeUser(1E-5,1E1);
//     hPtMC->DrawCopy("same");

//     TLegend* legend6 = new TLegend(0.57, 0.57, 0.87, 0.87);    
//     legend6->SetBorderSize(0);
//     legend6->SetFillColor(0);
//     legend6->AddEntry(hPt_clone, "p_{T} dist (Rec)", "L");
//     legend6->AddEntry(hPtMC, "p_{T} dist (MC)", "L");
//     legend6->Draw();

//     gPad->SetLogy();

//     c6->cd(2);
//     TH1F* hPtPion_clone = static_cast<TH1F*>(hPtPionCorr->Clone()); 
//     hPtPion_clone->GetYaxis()->SetRangeUser(1E-2,1E1);
//     hPtPion_clone->SetTitle("p_{T} distribution (Pions)");
//     hPtPion_clone->SetLineColor(1);
//     hPtPion_clone->DrawCopy();

//     TH1F* hPtMCPion = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtMCPion"));
//     hPtMCPion->Scale(1.0/nEvents/normalization/hPtMCPion->GetBinWidth(1));
//     hPtMCPion->GetYaxis()->SetRangeUser(1E-2,1E1);
//     hPtMCPion->DrawCopy("same");

//     TLegend* legend7 = new TLegend(0.57, 0.57, 0.87, 0.87);    
//     legend7->SetBorderSize(0);
//     legend7->SetFillColor(0);
//     legend7->AddEntry(hPtPion_clone, "p_{T} dist (Rec)", "L");
//     legend7->AddEntry(hPtMCPion, "p_{T} dist (MC)", "L");
//     legend7->Draw();

//     gPad->SetLogy();

//     c6->cd(3);
//     TH1F* hPtKaon_clone = static_cast<TH1F*>(hPtKaonCorr->Clone()); 
//     hPtKaon_clone->GetYaxis()->SetRangeUser(1E-2,1E0);
//     hPtKaon_clone->SetLineColor(1);
//     hPtKaon_clone->DrawCopy();

//     TH1F* hPtMCKaon = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtMCKaon"));
//     hPtMCKaon->Scale(1.0/nEvents/normalization/hPtMCKaon->GetBinWidth(1));
//     hPtMCKaon->GetYaxis()->SetRangeUser(1E-2,1E0);
//     hPtMCKaon->DrawCopy("same");

//     TLegend* legend8 = new TLegend(0.57, 0.57, 0.87, 0.87);    
//     legend8->SetBorderSize(0);
//     legend8->SetFillColor(0);
//     legend8->AddEntry(hPtKaon_clone, "p_{T} dist (Rec)", "L");
//     legend8->AddEntry(hPtMCKaon, "p_{T} dist (MC)", "L");
//     legend8->Draw();

//     gPad->SetLogy(); 

//     c6->cd(4);
//     TH1F* hPtProton_clone = static_cast<TH1F*>(hPtProtonCorr->Clone()); 
//     hPtProton_clone->GetYaxis()->SetRangeUser(1E-2,1E-1);
//     hPtProton_clone->SetLineColor(1);
//     hPtProton_clone->DrawCopy();

//     TH1F* hPtMCProton = dynamic_cast<TH1F*> (fListOfHists->FindObject("fdNdPtMCProton"));
//     hPtMCProton->Scale(1.0/nEvents/normalization/hPtMCProton->GetBinWidth(1));
//     hPtMCProton->GetYaxis()->SetRangeUser(1E-2,1E-1);
//     hPtMCProton->DrawCopy("same");

//     TLegend* legend9 = new TLegend(0.2, 0.25, 0.5, 0.55);    
//     legend9->SetBorderSize(0);
//     legend9->SetFillColor(0);
//     legend9->AddEntry(hPtProton_clone, "p_{T} dist (Rec)", "L");
//     legend9->AddEntry(hPtMCProton, "p_{T} dist (MC)", "L");
//     legend9->Draw();

//     gPad->SetLogy(); 

//     c6->Update();
//   }

}
