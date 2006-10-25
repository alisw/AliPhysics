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

/* $Id$ */

// Macro (upgraded version of MUONmassPlot_ESD.C, better handling of Jpsi) to make : 
// 1) Ntuple (Ktuple) containing Upsilon kinematics variables (from kinematics.root files) 
// 2) Ntuple (ESDtuple) containing Upsilon kinematics variables from reconstruction and 
// combinations of 2 muons with opposite charges (ESDtupleBck will be used later)
// 3) Some QA histograms
// Ntuple are stored in the file MUONefficiency.root and  ESD tree and QA histograms in AliESDs.root

// Christophe Suire, IPN Orsay



// Arguments:
//   FirstEvent (default 0)
//   LastEvent (default 1.e6)
//   ResType (default 553)
//      553 for Upsilon, 443 for J/Psi
//   Chi2Cut (default 100)
//      to keep only tracks with chi2 per d.o.f. < Chi2Cut



#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TTree.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TParticle.h"
#include "TTree.h"
#include "TString.h"
#include <Riostream.h>

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliMagFMaps.h"
#include "AliESD.h"
#include "AliTracker.h"

// MUON includes
#include "AliMUONTrackParam.h"
#include "AliESDMuonTrack.h"
#endif


Bool_t MUONefficiency( Int_t ResType = 553, Int_t FirstEvent = 0, Int_t LastEvent = 1000000,
		       char* esdFileName = "AliESDs.root", char* filename = "galice.root")
{ // MUONefficiency starts

  Double_t MUON_MASS = 0.105658369;
  Double_t UPSILON_MASS = 9.4603 ;
  Double_t JPSI_MASS = 3.097;

  // Upper and lower bound for counting entries in the mass peak
  // +/- 300 MeV/c^2 in this case.
  Float_t countingRange = 0.300 ;  
  
  Float_t massResonance = 5.;
  Float_t invMassMinInPeak = 0. ; 
  Float_t invMassMaxInPeak = 0. ; 
  
  Float_t nBinsPerGev = 40 ; 
  Float_t invMassMin = 0;   Float_t invMassMax = 20; 
  Float_t ptMinResonance = 0 ; Float_t ptMaxResonance = 20 ; Int_t ptBinsResonance = 100; 

  if (ResType==443) {
    massResonance = JPSI_MASS ;
    invMassMinInPeak =  JPSI_MASS - countingRange  ; invMassMaxInPeak = JPSI_MASS + countingRange ; 
    //limits for histograms
    invMassMin = 0 ; invMassMax = 6.;
    ptMinResonance = 0 ; ptMaxResonance = 20 ; ptBinsResonance = 100; 
  }
  if (ResType==553) {
    massResonance = UPSILON_MASS;
    invMassMinInPeak = UPSILON_MASS - countingRange ; invMassMaxInPeak = UPSILON_MASS + countingRange; 
    //limits for histograms 
    invMassMin = 0 ; invMassMax = 12.;
    ptMinResonance = 0 ; ptMaxResonance = 20 ; ptBinsResonance = 100; 
  }
  
  // Single Tracks muon cuts
  Float_t Chi2Cut = 100.;
  Float_t PtCutMin = 0. ;
  Float_t PtCutMax = 10000. ; 


  // Limits for histograms 
  Float_t ptMinMuon = 0. ; Float_t ptMaxMuon = 20.; Int_t ptBinsMuon = 100 ;
  Float_t pMinMuon = 0.  ; Float_t pMaxMuon = 200.; Int_t pBinsMuon = 100 ;
 

  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // Printing Level 
  Int_t PRINTLEVEL = 0 ;
  
  //for kinematic, i.e. reference tracks
  TNtuple *Ktuple = new TNtuple("Ktuple","Kinematics NTuple","ev:npart:id:idmo:idgdmo:p:pt:y:theta:pseudorap:vx:vy:vz");

  //for reconstruction  
  TH1F *hPtMuon = new TH1F("hPtMuon", "Muon Pt (GeV/c)", ptBinsMuon, ptMinMuon, ptMaxMuon);
  TH1F *hPtMuonPlus = new TH1F("hPtMuonPlus", "Muon+ Pt (GeV/c)", ptBinsMuon, ptMinMuon, ptMaxMuon);
  TH1F *hPtMuonMinus = new TH1F("hPtMuonMinus", "Muon- Pt (GeV/c)", ptBinsMuon, ptMinMuon, ptMaxMuon);
  TH1F *hPMuon = new TH1F("hPMuon", "Muon P (GeV/c)", pBinsMuon, pMinMuon, pMaxMuon);
  
  TH1F *hInvMassAll;
  TH1F *hInvMassBg;
  TH2F *hInvMassAll_vs_Pt;
  TH2F *hInvMassBgk_vs_Pt;
  TH1F *hInvMassRes;


  hInvMassAll = new TH1F("hInvMassAll", "Mu+Mu- invariant mass (GeV/c2)", (Int_t) (nBinsPerGev*(invMassMax - invMassMin)), invMassMin, invMassMax);
  hInvMassBg = new TH1F("hInvMassBg", "Mu+Mu- invariant mass BG(GeV/c2)", (Int_t) (nBinsPerGev*(invMassMax- invMassMin)), invMassMin, invMassMax);
  hInvMassAll_vs_Pt = new TH2F("hInvMassAll_vs_Pt","hInvMassAll_vs_Pt",(Int_t) (nBinsPerGev*(invMassMax- invMassMin)), invMassMin, invMassMax,ptBinsResonance,ptMinResonance,ptMaxResonance);
  hInvMassBgk_vs_Pt = new TH2F("hInvMassBgk_vs_Pt","hInvMassBgk_vs_Pt",(Int_t) (nBinsPerGev*(invMassMax- invMassMin)), invMassMin, invMassMax,ptBinsResonance,ptMinResonance,ptMaxResonance);
  
  hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around Resonance",(Int_t) (nBinsPerGev*3*countingRange*2),massResonance-3*countingRange,massResonance+3*countingRange);

  TH1F *hPrimaryVertex = new TH1F("hPrimaryVertex","SPD reconstructed Z vertex",150,-15,15);
  TH1F *hChi2PerDof = new TH1F("hChi2PerDof", "Muon track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hNumberOfTrack = new TH1F("hNumberOfTrack","nb of track /evt ",20,-0.5,19.5);
  TH1F *hRapMuon = new TH1F("hRapMuon"," Muon Rapidity",50,-4.5,-2);
  TH1F *hRapResonance = new TH1F("hRapResonance"," Resonance Rapidity",50,-4.5,-2);
  TH1F *hPtResonance = new TH1F("hPtResonance", "Resonance Pt (GeV/c)", 100, 0., 20.);
  TH2F *hThetaPhiPlus = new TH2F("hThetaPhiPlus", "Theta vs Phi +", 760, -190., 190., 400, 160., 180.);
  TH2F *hThetaPhiMinus = new TH2F("hThetaPhiMinus", "Theta vs Phi -", 760, -190., 190., 400, 160., 180.);
  
  TNtuple *ESDtuple = new TNtuple("ESDtuple","Reconstructed Mu+Mu- pairs and Upsilon","ev:tw:pt:y:theta:minv:pt1:y1:theta1:q1:trig1:pt2:y2:theta2:q2:trig2");
  TNtuple *ESDtupleBck = new TNtuple("ESDtupleBck","Reconstructed Mu+Mu- pairs for Background","ev:pt:y:theta:minv:pt1:y1:theta1:pt2:y2:theta2");


  // Variables
  Int_t EventInMass = 0;
  Int_t EventInMassMatch = 0;
  Int_t NbTrigger = 0;
  Int_t ptTrig = 0;

  Double_t fXVertex=0;
  Double_t fYVertex=0;
  Double_t fZVertex=0;

  Double_t thetaX, thetaY, pYZ;
  Double_t fPxRec1, fPyRec1, fPzRec1, fE1;
  Double_t fPxRec2, fPyRec2, fPzRec2, fE2;
  Int_t fCharge1, fCharge2;

  Int_t ntrackhits, nevents; 
  Int_t nprocessedevents = 0 ;
  Double_t fitfmin;

  TLorentzVector fV1, fV2, fVtot;

  // set  mag field 
  // waiting for mag field in CDB 
  printf("Loading field map...\n");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open(filename);
  if (!runLoader) {
    Error("MUONefficiency", "getting run loader from file %s failed", filename);
    return kFALSE;
  }

  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("MUONefficiency", "no galice object found");
    return kFALSE;
  }
  
  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("MUONefficiency", "opening ESD file %s failed", esdFileName);
    return kFALSE;
  }
  
  AliESD* esd = new AliESD();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  tree->SetBranchAddress("ESD", &esd);

  runLoader->LoadHeader();
  nevents = runLoader->GetNumberOfEvents();
  AliMUONTrackParam trackParam;

  // to access the particle  Stack
  runLoader->LoadKinematics("READ");

  Int_t numberOfGeneratedResonances = 0 ;

  TParticle *particle; 
  
  Int_t track1Trigger = 0 ;
  Float_t track1TriggerChi2 = 0 ;
  Int_t track2Trigger = 0 ;
  Float_t track2TriggerChi2 = 0 ;


  // Loop over events
  for (Int_t iEvent = FirstEvent; iEvent <= TMath::Min(LastEvent, nevents - 1); iEvent++) { // Start event loop

    if (iEvent%1000 == 0 )
      printf("\n Nb of events analysed: %d \n",iEvent);

    // get current event
    runLoader->GetEvent(iEvent);
    nprocessedevents++;

    // get the stack and fill the kine tree
    AliStack *theStack = runLoader->Stack();
    if (PRINTLEVEL > 0) theStack->DumpPStack ();    
    
    Int_t nparticles = (Int_t)runLoader->TreeK()->GetEntries();
    Int_t nprimarypart = theStack->GetNprimary();
    Int_t ntracks = theStack->GetNtrack();
  
    if (PRINTLEVEL || (iEvent%100==0)) printf("\n  >>> Event %d \n",iEvent);
    if (PRINTLEVEL) cout << nprimarypart << " Particles generated (total is " << ntracks << ")"<< endl ;    
    
    for(Int_t iparticle=0; iparticle<nparticles; iparticle++) { // Start loop over particles
      particle = theStack->Particle(iparticle);
      
      Int_t   muId = particle->GetPdgCode(); 
      Int_t   muM  = particle->GetFirstMother();
      Int_t   muGM = 0;
      Float_t muP  = particle->P();
      Float_t muPt = TMath::Sqrt(particle->Px()*particle->Px()+particle->Py()*particle->Py());
      Float_t muY  = 0.5*TMath::Log((particle->Energy()+particle->Pz()+1.e-13)/(particle->Energy()-particle->Pz()+1.e-13));
      if (muM >= 0) {
	TParticle *theMum = theStack->Particle(muM);
	muM  =  theMum->GetPdgCode();
	muGM  = theMum->GetFirstMother() ;
	if (muGM >= 0){
	  TParticle *grandMa = theStack->Particle(muGM);
	  muGM = grandMa->GetPdgCode();
	}
	else muGM=0;
      }
      else muM=0;
    
      if (muId==ResType) numberOfGeneratedResonances++;

  
      Float_t muT  = particle->Theta()*180/TMath::Pi();
      Float_t muE  = particle->Eta();
      
      Float_t muVx = particle->Vx();
      Float_t muVy = particle->Vy();
      Float_t muVz = particle->Vz();
      
      // If a write error occurs, the number of bytes returned is -1.
      // If no data are written, because e.g. the branch is disabled,
      // the number of bytes returned is 0.
      Int_t errCode = Ktuple->Fill(iEvent,nparticles,muId,muM,muGM,muP,muPt,muY,muT,muE,muVx,muVy,muVz);
      if (PRINTLEVEL || errCode < 1) printf("iEvent %d,nparticles %d,muId %d,muM %d,muGM %d,muP %.2f,muPt  %.2f,muY  %.2f,muT  %.2f,muE  %.2f,muVx  %.2f,muVy  %.2f,muVz  %.2f \n", iEvent,nparticles,muId,muM,muGM,muP,muPt,muY,muT,muE,muVx,muVy,muVz);

    } // End loop over particles
    

    
    // get the event summary data
    tree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return kFALSE;
    }
    
    // get the SPD reconstructed vertex (vertexer) and fill the histogram
    AliESDVertex* Vertex = (AliESDVertex*) esd->GetVertex();
    if (Vertex->GetNContributors()) {
      fZVertex = Vertex->GetZv();
      fYVertex = Vertex->GetYv();
      fXVertex = Vertex->GetXv();      
    }
    hPrimaryVertex->Fill(fZVertex);
    
    Int_t triggerWord = esd->GetTriggerMask();
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ; 

    if (PRINTLEVEL > 0){
      printf("\n Nb of events analysed: %d \n",iEvent);
      cout << " number of tracks: " << nTracks  <<endl;
    }

    // loop over all reconstructed tracks (also first track of combination)
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {

      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);

      if (!Vertex->GetNContributors()) {
	//re-extrapolate to vertex, if not kown before.
	trackParam.GetParamFrom(*muonTrack);
	trackParam.ExtrapToVertex(fXVertex, fYVertex, fZVertex);
	trackParam.SetParamFor(*muonTrack);
      }

      // Trigger
      if (PRINTLEVEL > 5) cout << "MatchTrigger " << muonTrack->GetMatchTrigger() << " and Chi2 of matching tracks " << track1TriggerChi2 <<   endl ;
      track1Trigger = muonTrack->GetMatchTrigger();
      if (track1Trigger)
	track1TriggerChi2 = muonTrack->GetChi2MatchTrigger();
      else 
	track1TriggerChi2 = 0. ;

      thetaX = muonTrack->GetThetaX();
      thetaY = muonTrack->GetThetaY();

      pYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      fPzRec1  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaY));
      fPxRec1  = fPzRec1 * TMath::Tan(thetaX);
      fPyRec1  = fPzRec1 * TMath::Tan(thetaY);
      fCharge1 = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
      
      fE1 = TMath::Sqrt(MUON_MASS * MUON_MASS + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
      fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);

      ntrackhits = muonTrack->GetNHit();
      fitfmin    = muonTrack->GetChi2();

      // transverse momentum
      Float_t pt1 = fV1.Pt();
      
      // total momentum
      Float_t p1 = fV1.P();
      
      // Rapidity
      Float_t rapMuon1 = fV1.Rapidity();
      
      // chi2 per d.o.f.
      
      Float_t ch1 =  fitfmin / (2.0 * ntrackhits - 5);
      if (PRINTLEVEL > 5 ) printf(" px %f py %f pz %f pt %f NHits %d  Norm.chi2 %f charge %d\n",fPxRec1, fPyRec1, fPzRec1, pt1, ntrackhits, ch1, fCharge1);
      
      
      if ((ch1 < Chi2Cut) && (pt1 > PtCutMin) && (pt1 < PtCutMax)) { // condition for good track (Chi2Cut and PtCut)
	if (PRINTLEVEL > 8) cout << "inside pt and chi2 cuts " << endl ; 
	
	// fill histos hPtMuon and hChi2PerDof
	hPtMuon->Fill(pt1);
	hPMuon->Fill(p1);
	hChi2PerDof->Fill(ch1);
	hRapMuon->Fill(rapMuon1);

	if (fCharge1 > 0) {
	  hPtMuonPlus->Fill(pt1);
	  hThetaPhiPlus->Fill(TMath::ATan2(fPyRec1,fPxRec1)*180./TMath::Pi(),TMath::ATan2(pt1,fPzRec1)*180./3.1415);
	} else {
	  hPtMuonMinus->Fill(pt1);
	  hThetaPhiMinus->Fill(TMath::ATan2(fPyRec1,fPxRec1)*180./TMath::Pi(),TMath::ATan2(pt1,fPzRec1)*180./3.1415);
	}

	// loop over second track of combination
	for (Int_t iTrack2 = iTrack + 1; iTrack2 < nTracks; iTrack2++) {
	  
	  AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack2);
          
	  if (!Vertex->GetNContributors()) {
	    trackParam.GetParamFrom(*muonTrack);
	    trackParam.ExtrapToVertex(fXVertex, fYVertex, fZVertex);
	    trackParam.SetParamFor(*muonTrack);
	  }

	  track2Trigger = muonTrack->GetMatchTrigger();
	  if (track2Trigger) 
	    track2TriggerChi2 = muonTrack->GetChi2MatchTrigger();
	  else 
	    track2TriggerChi2 = 0. ;

	  thetaX = muonTrack->GetThetaX();
	  thetaY = muonTrack->GetThetaY();

	  pYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
	  fPzRec2  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaY));
	  fPxRec2  = fPzRec2 * TMath::Tan(thetaX);
	  fPyRec2  = fPzRec2 * TMath::Tan(thetaY);
	  fCharge2 = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));

	  fE2 = TMath::Sqrt(MUON_MASS * MUON_MASS + fPxRec2 * fPxRec2 + fPyRec2 * fPyRec2 + fPzRec2 * fPzRec2);
	  fV2.SetPxPyPzE(fPxRec2, fPyRec2, fPzRec2, fE2);

	  ntrackhits = muonTrack->GetNHit();
	  fitfmin    = muonTrack->GetChi2();

	  // transverse momentum
	  Float_t pt2 = fV2.Pt();

	  // chi2 per d.o.f.
	  Float_t ch2 = fitfmin  / (2.0 * ntrackhits - 5);


	  // condition for good track (Chi2Cut and PtCut)
	  if ((ch2 < Chi2Cut) && (pt2 > PtCutMin)  && (pt2 < PtCutMax)) {
	    
	    // condition for opposite charges
	    if ((fCharge1 * fCharge2) == -1) {
	      
	      if (PRINTLEVEL > 8) cout << "---------> Now filling the Ntuple "  <<  endl ;
	      
	      // invariant mass
	      fVtot = fV1 + fV2;
	      Float_t invMass = fVtot.M();
	      
	      if (fCharge1 < 0){ //mu_minus is index 1 in the ntuple
		Float_t ESDFill[16] = {iEvent,triggerWord,fVtot.Pt(),fVtot.Rapidity(),fVtot.Theta()/TMath::Pi()*180,invMass,fV1.Pt(),fV1.Rapidity(),fV1.Theta()/TMath::Pi()*180,fCharge1,track1TriggerChi2,fV2.Pt(),fV2.Rapidity(),fV2.Theta()/TMath::Pi()*180,fCharge2,track2TriggerChi2};
		ESDtuple->Fill(ESDFill);
	      }
	      else{
		Float_t ESDFill[16] = {iEvent,triggerWord,fVtot.Pt(),fVtot.Rapidity(),fVtot.Theta()/TMath::Pi()*180,invMass,fV2.Pt(),fV2.Rapidity(),fV2.Theta()/TMath::Pi()*180,fCharge2,track2TriggerChi2,fV1.Pt(),fV1.Rapidity(),fV1.Theta()/TMath::Pi()*180,fCharge1,track1TriggerChi2};
		ESDtuple->Fill(ESDFill);
	      }
	      
	      // fill histos hInvMassAll and hInvMassRes
	      hInvMassAll->Fill(invMass);
	      hInvMassRes->Fill(invMass);
	      hInvMassAll_vs_Pt->Fill(invMass,fVtot.Pt());
	      
	      //trigger info 
	      if (ResType == 553)
		ptTrig = 0x20;// mask for Hpt unlike sign pair
	      else if (ResType == 443)
		ptTrig = 0x10;// mask for Lpt unlike sign pair
	      
	      
	      if (esd->GetTriggerMask() &  ptTrig) NbTrigger++;
	      
	      if (invMass > invMassMinInPeak && invMass < invMassMaxInPeak) {
		EventInMass++;
		hRapResonance->Fill(fVtot.Rapidity());
		hPtResonance->Fill(fVtot.Pt());
		
		// match with trigger
		if (muonTrack->GetMatchTrigger() && (esd->GetTriggerMask() & ptTrig))  EventInMassMatch++;
		
	      }
	      
	    } // if (fCharge1 * fCharge2) == -1)
	  } // if ((ch2 < Chi2Cut) && (pt2 > PtCutMin) && (pt2 < PtCutMax))
	} //  for (Int_t iTrack2 = iTrack + 1; iTrack2 < iTrack; iTrack2++)
      } // if (ch1 < Chi2Cut) && (pt1 > PtCutMin)&& (pt1 < PtCutMax) )
    } // for (Int_t iTrack = 0; iTrack < nrectracks; iTrack++)
    
    hNumberOfTrack->Fill(nTracks);
    //    esdFile->Delete();
  
  } // End of event loop


  // Loop over events for bg event

  Double_t thetaPlus,  phiPlus;
  Double_t thetaMinus, phiMinus;
  Float_t PtMinus, PtPlus;
  
  for (Int_t iEvent = 0; iEvent < hInvMassAll->Integral(); iEvent++) {  // Loop over events for bg event
    // according to Christian a 3d phi-theta-pt random pick  would take better care 
    // of all correlations 

    hThetaPhiPlus->GetRandom2(phiPlus, thetaPlus);
    hThetaPhiMinus->GetRandom2(phiMinus,thetaMinus);
    PtPlus = hPtMuonPlus->GetRandom();
    PtMinus = hPtMuonMinus->GetRandom();

    fPxRec1  = PtPlus * TMath::Cos(TMath::Pi()/180.*phiPlus);
    fPyRec1  = PtPlus * TMath::Sin(TMath::Pi()/180.*phiPlus);
    fPzRec1  = PtPlus / TMath::Tan(TMath::Pi()/180.*thetaPlus);

    fE1 = TMath::Sqrt(MUON_MASS * MUON_MASS + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
    fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);

    fPxRec2  = PtMinus * TMath::Cos(TMath::Pi()/180.*phiMinus);
    fPyRec2  = PtMinus * TMath::Sin(TMath::Pi()/180.*phiMinus);
    fPzRec2  = PtMinus / TMath::Tan(TMath::Pi()/180.*thetaMinus);

    fE2 = TMath::Sqrt(MUON_MASS * MUON_MASS + fPxRec2 * fPxRec2 + fPyRec2 * fPyRec2 + fPzRec2 * fPzRec2);
    fV2.SetPxPyPzE(fPxRec2, fPyRec2, fPzRec2, fE2);

    // invariant mass
    fVtot = fV1 + fV2;
      
    // fill histos hInvMassAll and hInvMassRes
    hInvMassBg->Fill(fVtot.M());
    hInvMassBgk_vs_Pt->Fill( fVtot.M(), fVtot.Pt() );

    // Ntuple for background... more convenient
    ESDtupleBck->Fill(iEvent,fVtot.Pt(),fVtot.Rapidity(),fVtot.Theta()/TMath::Pi()*180,fVtot.M(),fV2.Pt(),fV2.Rapidity(),fV2.Theta()/TMath::Pi()*180,fV1.Pt(),fV1.Rapidity(),fV1.Theta()/TMath::Pi()*180);

  } // End loop over events for background


  // File for histograms and histogram booking
  TString outfilename = "MUONefficiency.root";
  TFile *ntupleFile = new TFile(outfilename.Data(), "RECREATE");  
  
  Ktuple->Write();
  ESDtuple->Write();
  ESDtupleBck->Write();

  ntupleFile->Close();
  
  TFile *histoFile = new TFile("MUONhistos.root", "RECREATE");  
  hPrimaryVertex->Write();
  hPtMuon->Write();
  hPtMuonPlus->Write();
  hPtMuonMinus->Write();
  hPMuon->Write();
  hChi2PerDof->Write();
  hInvMassAll->Write();
  hInvMassBg->Write();
  hInvMassAll_vs_Pt ->Write();
  hInvMassBgk_vs_Pt->Write();
  hInvMassRes->Write();
  hNumberOfTrack->Write();
  hRapMuon ->Write();
  hRapResonance ->Write();
  hPtResonance ->Write();
  hThetaPhiPlus ->Write();
  hThetaPhiMinus ->Write();
  histoFile->Close();

  cout << "" << endl ;
  cout << "*************************************************" << endl;
 
  cout << "MUONefficiency : " << nprocessedevents  << " events processed" << endl;
  if (ResType==443)
    cout << "Number of generated J/Psi (443)  : " <<  numberOfGeneratedResonances  << endl ;
  if (ResType==553)
    cout << "Number of generated Upsilon (553)  :" <<  numberOfGeneratedResonances  << endl ;
  cout << "Chi2Cut for muon tracks = " << Chi2Cut << endl;
  cout << "PtCutMin for muon tracks = " << PtCutMin << endl;
  cout << "PtCutMax for muon tracks = " << PtCutMax << endl;
  cout << "Entries (unlike sign dimuons) in the mass range  ["<<invMassMinInPeak<<";"<<invMassMaxInPeak<<"] : " << EventInMass <<endl;
  if (ptTrig==0x800) cout << "Unlike Pair - All Pt" ;   
  if (ptTrig==0x400) cout << "Unlike Pair - High Pt" ;   
  if (ptTrig==0x200) cout << "Unlike Pair - Low Pt" ; 
  cout << " triggers : " << NbTrigger << endl;
  
  cout << "Entries in the mass range with matching between reconstructed tracks and trigger tracks " << EventInMassMatch << endl;


  return kTRUE;
}
