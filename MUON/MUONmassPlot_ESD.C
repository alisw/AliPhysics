#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TParticle.h"
#include "TTree.h"
#include <Riostream.h>

// STEER includes
#include "AliRun.h"
#include "AliLog.h"
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
//
// Macro MUONmassPlot.C for ESD
// Ch. Finck, Subatech, April. 2004
//

// macro to make invariant mass plots
// for combinations of 2 muons with opposite charges,
// from root file "MUON.tracks.root" containing the result of track reconstruction.
// Histograms are stored on the "MUONmassPlot.root" file.
// introducing TLorentzVector for parameter calculations (Pt, P,rap,etc...)
// using Invariant Mass for rapidity.

// Arguments:
//   FirstEvent (default 0)
//   LastEvent (default 0)
//   ResType (default 553)
//      553 for Upsilon, anything else for J/Psi
//   Chi2Cut (default 100)
//      to keep only tracks with chi2 per d.o.f. < Chi2Cut
//   PtCutMin (default 1)
//      to keep only tracks with transverse momentum > PtCutMin
//   PtCutMax (default 10000)
//      to keep only tracks with transverse momentum < PtCutMax
//   massMin (default 9.17 for Upsilon) 
//      &  massMax (default 9.77 for Upsilon) 
//         to calculate the reconstruction efficiency for resonances with invariant mass
//         massMin < mass < massMax.

// Add parameters and histograms for analysis 

Bool_t MUONmassPlot(char* filename = "galice.root", Int_t FirstEvent = 0, Int_t LastEvent = 10000,
		  char* esdFileName = "AliESDs.root", Int_t ResType = 553, 
                  Float_t Chi2Cut = 100., Float_t PtCutMin = 1., Float_t PtCutMax = 10000.,
                  Float_t massMin = 9.17,Float_t massMax = 9.77)
{
  cout << "MUONmassPlot " << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "ResType " << ResType << endl;
  cout << "Chi2Cut " << Chi2Cut << endl;
  cout << "PtCutMin " << PtCutMin << endl;
  cout << "PtCutMax " << PtCutMax << endl;
  cout << "massMin " << massMin << endl;
  cout << "massMax " << massMax << endl;

 
  //Reset ROOT and connect tree file
  gROOT->Reset();

  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONmassPlot.root", "RECREATE");
  TH1F *hPtMuon = new TH1F("hPtMuon", "Muon Pt (GeV/c)", 100, 0., 20.);
  TH1F *hPtMuonPlus = new TH1F("hPtMuonPlus", "Muon+ Pt (GeV/c)", 100, 0., 20.);
  TH1F *hPtMuonMinus = new TH1F("hPtMuonMinus", "Muon- Pt (GeV/c)", 100, 0., 20.);
  TH1F *hPMuon = new TH1F("hPMuon", "Muon P (GeV/c)", 100, 0., 200.);
  TH1F *hChi2PerDof = new TH1F("hChi2PerDof", "Muon track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hInvMassAll = new TH1F("hInvMassAll", "Mu+Mu- invariant mass (GeV/c2)", 480, 0., 12.);
  TH1F *hInvMassBg = new TH1F("hInvMassBg", "Mu+Mu- invariant mass BG(GeV/c2)", 480, 0., 12.);
  TH2F *hInvMassAll_vs_Pt = new TH2F("hInvMassAll_vs_Pt","hInvMassAll_vs_Pt",480,0.,12.,80,0.,20.);
  TH2F *hInvMassBgk_vs_Pt = new TH2F("hInvMassBgk_vs_Pt","hInvMassBgk_vs_Pt",480,0.,12.,80,0.,20.);
  TH1F *hInvMassRes;
  TH1F *hPrimaryVertex = new TH1F("hPrimaryVertex","SPD reconstructed Z vertex",150,-15,15);

  if (ResType == 553) {
    hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around Upsilon", 60, 8., 11.);
  } else {
    hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around J/Psi", 80, 0., 5.);
  }

  TH1F *hNumberOfTrack = new TH1F("hNumberOfTrack","nb of track /evt ",20,-0.5,19.5);
  TH1F *hRapMuon = new TH1F("hRapMuon"," Muon Rapidity",50,-4.5,-2);
  TH1F *hRapResonance = new TH1F("hRapResonance"," Resonance Rapidity",50,-4.5,-2);
  TH1F *hPtResonance = new TH1F("hPtResonance", "Resonance Pt (GeV/c)", 100, 0., 20.);
  TH2F *hThetaPhiPlus = new TH2F("hThetaPhiPlus", "Theta vs Phi +", 760, -190., 190., 400, 160., 180.);
  TH2F *hThetaPhiMinus = new TH2F("hThetaPhiMinus", "Theta vs Phi -", 760, -190., 190., 400, 160., 180.);


  // settings
  Int_t EventInMass = 0;
  Int_t EventInMassMatch = 0;
  Int_t NbTrigger = 0;

  Float_t muonMass = 0.105658389;
//   Float_t UpsilonMass = 9.46037;
//   Float_t JPsiMass = 3.097;

  Double_t thetaX, thetaY, pYZ;
  Double_t fPxRec1, fPyRec1, fPzRec1, fE1;
  Double_t fPxRec2, fPyRec2, fPzRec2, fE2;
  Int_t fCharge, fCharge2;

  Int_t ntrackhits, nevents;
  Double_t fitfmin;
  Double_t fZVertex=0;
  Double_t fYVertex=0;
  Double_t fXVertex=0;
 
  TLorentzVector fV1, fV2, fVtot;

  // set mag field
  // waiting for mag field in CDB 
  printf("Loading field map...\n");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k4kG);
  AliTracker::SetFieldMap(field, kFALSE);

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open(filename);
  if (!runLoader) {
    Error("MUONmass_ESD", "getting run loader from file %s failed", 
	    filename);
    return kFALSE;
  }

  if (!gAlice) {
    Error("MUONmass_ESD", "no galice object found");
    return kFALSE;
  }
  

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("MUONmass_ESD", "opening ESD file %s failed", esdFileName);
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

  // Loop over events
  for (Int_t iEvent = FirstEvent; iEvent <= TMath::Min(LastEvent, nevents - 1); iEvent++) {

    // get current event
    runLoader->GetEvent(iEvent);
   
    // get the event summary data
    tree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return kFALSE;
    }

    // get the SPD reconstructed vertex (vertexer) and fill the histogram
    AliESDVertex* Vertex = (AliESDVertex*) esd->AliESD::GetVertex();

    if (Vertex) {
      fZVertex = Vertex->GetZv();
      fYVertex = Vertex->GetYv();
      fXVertex = Vertex->GetXv();

    }
    hPrimaryVertex->Fill(fZVertex);

    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ; 

    //    printf("\n Nb of events analysed: %d\r",iEvent);
    //      cout << " number of tracks: " << nTracks  <<endl;
  
    // loop over all reconstructed tracks (also first track of combination)
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {

      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);

      if (!Vertex) {
	//re-extrapolate to vertex, if not kown before.
	trackParam.GetParamFrom(*muonTrack);
	trackParam.ExtrapToVertex(fXVertex, fYVertex, fZVertex);
	trackParam.SetParamFor(*muonTrack);
      }
      thetaX = muonTrack->GetThetaX();
      thetaY = muonTrack->GetThetaY();

      pYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      fPzRec1  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaY));
      fPxRec1  = fPzRec1 * TMath::Tan(thetaX);
      fPyRec1  = fPzRec1 * TMath::Tan(thetaY);
      fCharge = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));

      fE1 = TMath::Sqrt(muonMass * muonMass + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
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
//      printf(" px %f py %f pz %f NHits %d  Norm.chi2 %f charge %d\n", 
// 	     fPxRec1, fPyRec1, fPzRec1, ntrackhits, ch1, fCharge);

      // condition for good track (Chi2Cut and PtCut)

      if ((ch1 < Chi2Cut) && (pt1 > PtCutMin) && (pt1 < PtCutMax)) {

	// fill histos hPtMuon and hChi2PerDof
	hPtMuon->Fill(pt1);
	hPMuon->Fill(p1);
	hChi2PerDof->Fill(ch1);
	hRapMuon->Fill(rapMuon1);
	if (fCharge > 0) {
	  hPtMuonPlus->Fill(pt1);
	  hThetaPhiPlus->Fill(TMath::ATan2(fPyRec1,fPxRec1)*180./TMath::Pi(),TMath::ATan2(pt1,fPzRec1)*180./3.1415);
	} else {
	  hPtMuonMinus->Fill(pt1);
	  hThetaPhiMinus->Fill(TMath::ATan2(fPyRec1,fPxRec1)*180./TMath::Pi(),TMath::ATan2(pt1,fPzRec1)*180./3.1415);
	}
	// loop over second track of combination
	for (Int_t iTrack2 = iTrack + 1; iTrack2 < nTracks; iTrack2++) {
	  
	  AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack2);

	  if (!Vertex) {
	    trackParam.GetParamFrom(*muonTrack);
	    trackParam.ExtrapToVertex(fXVertex, fYVertex, fZVertex);
	    trackParam.SetParamFor(*muonTrack);
	  }

	  thetaX = muonTrack->GetThetaX();
	  thetaY = muonTrack->GetThetaY();

	  pYZ      =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
	  fPzRec2  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaY));
	  fPxRec2  = fPzRec2 * TMath::Tan(thetaX);
	  fPyRec2  = fPzRec2 * TMath::Tan(thetaY);
	  fCharge2 = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));

	  fE2 = TMath::Sqrt(muonMass * muonMass + fPxRec2 * fPxRec2 + fPyRec2 * fPyRec2 + fPzRec2 * fPzRec2);
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
	    if ((fCharge * fCharge2) == -1) {

	      // invariant mass
	      fVtot = fV1 + fV2;
	      Float_t invMass = fVtot.M();
	            
	      // fill histos hInvMassAll and hInvMassRes
	      hInvMassAll->Fill(invMass);
	      hInvMassRes->Fill(invMass);
	      hInvMassAll_vs_Pt->Fill(invMass,fVtot.Pt());
	      Int_t ptTrig;
	      if (ResType == 553) 
		ptTrig =  0x400;// mask for Hpt unlike sign pair
	      else 
		ptTrig =  0x200;// mask for Lpt unlike sign pair

	      if (esd->GetTriggerMask() &  ptTrig) NbTrigger++; 
	      if (invMass > massMin && invMass < massMax) {
		EventInMass++;
		if (muonTrack->GetMatchTrigger() && (esd->GetTriggerMask() & ptTrig))// match with trigger
		  EventInMassMatch++;

  		hRapResonance->Fill(fVtot.Rapidity());
  		hPtResonance->Fill(fVtot.Pt());
	      }

	    } // if (fCharge * fCharge2) == -1)
	  } // if ((ch2 < Chi2Cut) && (pt2 > PtCutMin) && (pt2 < PtCutMax))
	} //  for (Int_t iTrack2 = iTrack + 1; iTrack2 < iTrack; iTrack2++)
      } // if (ch1 < Chi2Cut) && (pt1 > PtCutMin)&& (pt1 < PtCutMax) )
    } // for (Int_t iTrack = 0; iTrack < nrectracks; iTrack++)

    hNumberOfTrack->Fill(nTracks);
    //    esdFile->Delete();
  } // for (Int_t iEvent = FirstEvent;

// Loop over events for bg event

  Double_t thetaPlus,  phiPlus;
  Double_t thetaMinus, phiMinus;
  Float_t PtMinus, PtPlus;
  
  for (Int_t iEvent = 0; iEvent < hInvMassAll->Integral(); iEvent++) {

    hThetaPhiPlus->GetRandom2(phiPlus, thetaPlus);
    hThetaPhiMinus->GetRandom2(phiMinus,thetaMinus);
    PtPlus = hPtMuonPlus->GetRandom();
    PtMinus = hPtMuonMinus->GetRandom();

    fPxRec1  = PtPlus * TMath::Cos(TMath::Pi()/180.*phiPlus);
    fPyRec1  = PtPlus * TMath::Sin(TMath::Pi()/180.*phiPlus);
    fPzRec1  = PtPlus / TMath::Tan(TMath::Pi()/180.*thetaPlus);

    fE1 = TMath::Sqrt(muonMass * muonMass + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
    fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);

    fPxRec2  = PtMinus * TMath::Cos(TMath::Pi()/180.*phiMinus);
    fPyRec2  = PtMinus * TMath::Sin(TMath::Pi()/180.*phiMinus);
    fPzRec2  = PtMinus / TMath::Tan(TMath::Pi()/180.*thetaMinus);

    fE2 = TMath::Sqrt(muonMass * muonMass + fPxRec2 * fPxRec2 + fPyRec2 * fPyRec2 + fPzRec2 * fPzRec2);
    fV2.SetPxPyPzE(fPxRec2, fPyRec2, fPzRec2, fE2);

    // invariant mass
    fVtot = fV1 + fV2;
      
    // fill histos hInvMassAll and hInvMassRes
    hInvMassBg->Fill(fVtot.M());
    hInvMassBgk_vs_Pt->Fill( fVtot.M(), fVtot.Pt() );
  }

  histoFile->Write();
  histoFile->Close();

  cout << endl;
  cout << "EventInMass " << EventInMass << endl;
  cout << "NbTrigger " << NbTrigger << endl;
  cout << "EventInMass match with trigger " << EventInMassMatch << endl;

  return kTRUE;
}

