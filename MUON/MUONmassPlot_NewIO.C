#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TParticle.h"
#include "TTree.h"

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliESDMuonTrack.h"
#endif
//
// Macro MUONmassPlot.C for new I/O
// Ch. Finck, Subatech, Jan. 2004
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

void MUONmassPlot(char* filename="galice.root", Int_t FirstEvent = 0, Int_t LastEvent = 0, Int_t ResType = 553, 
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
  TH1F *hPMuon = new TH1F("hPMuon", "Muon P (GeV/c)", 100, 0., 200.);
  TH1F *hChi2PerDof = new TH1F("hChi2PerDof", "Muon track chi2/d.o.f.", 100, 0., 20.);
  TH1F *hInvMassAll = new TH1F("hInvMassAll", "Mu+Mu- invariant mass (GeV/c2)", 480, 0., 12.);
  TH1F *hInvMassRes;

  if (ResType == 553) {
    hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around Upsilon", 60, 8., 11.);
  } else {
    hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around J/Psi", 80, 0., 5.);
  }

  TH1F *hNumberOfTrack = new TH1F("hNumberOfTrack","nb of track /evt ",20,-0.5,19.5);
  TH1F *hRapMuon = new TH1F("hRapMuon"," Muon Rapidity",50,-4.5,-2);
  TH1F *hRapResonance = new TH1F("hRapResonance"," Resonance Rapidity",50,-4.5,-2);
  TH1F *hPtResonance = new TH1F("hPtResonance", "Resonance Pt (GeV/c)", 100, 0., 20.);


  // settings
  Int_t EventInMass = 0;
  Float_t muonMass = 0.105658389;
//   Float_t UpsilonMass = 9.46037;
//   Float_t JPsiMass = 3.097;

  Double_t bendingSlope, nonBendingSlope, pYZ;
  Double_t fPxRec1, fPyRec1, fPzRec1, fZRec1, fE1;
  Double_t fPxRec2, fPyRec2, fPzRec2, fZRec2, fE2;
  Int_t fCharge, fCharge2;

  Int_t ntrackhits, nevents;
  Double_t fitfmin;

  TClonesArray * recTracksArray;
  TLorentzVector fV1, fV2, fVtot;
  
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadTracks("READ");

  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
  nevents = RunLoader->GetNumberOfEvents();
  
  AliMUONTrack * rectrack;
  AliMUONTrackParam *trackParam;
        
  // Loop over events
  for (Int_t ievent = FirstEvent; ievent <= TMath::Min(LastEvent, nevents - 1); ievent++) {

    // get current event
    RunLoader->GetEvent(ievent);
    
    muondata.SetTreeAddress("RT");
    muondata.GetRecTracks();
    recTracksArray = muondata.RecTracks();

    Int_t nrectracks = (Int_t) recTracksArray->GetEntriesFast(); //

    printf("\n Nb of events analysed: %d\r",ievent);
    //   cout << " number of tracks: " << nrectracks  <<endl;
  
    // loop over all reconstructed tracks (also first track of combination)
     for (Int_t irectracks = 0; irectracks <  nrectracks;  irectracks++) {

      rectrack = (AliMUONTrack*) recTracksArray->At(irectracks);

      trackParam = rectrack->GetTrackParamAtVertex();
      bendingSlope   = trackParam->GetBendingSlope();
      nonBendingSlope = trackParam->GetNonBendingSlope();

      pYZ     = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
      fPzRec1  = - pYZ / TMath::Sqrt(1.0 + bendingSlope * bendingSlope); // spectro. (z<0)
      fPxRec1  = fPzRec1 * nonBendingSlope;
      fPyRec1  = fPzRec1 * bendingSlope;
      fZRec1   = trackParam->GetZ();
      fCharge = Int_t(TMath::Sign(1., trackParam->GetInverseBendingMomentum()));

      fE1 = TMath::Sqrt(muonMass * muonMass + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
      fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);

      ntrackhits = rectrack->GetNTrackHits();
      fitfmin = rectrack->GetFitFMin();

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

	// loop over second track of combination
	for (Int_t irectracks2 = irectracks + 1; irectracks2 < nrectracks; irectracks2++) {

	  rectrack = (AliMUONTrack*) recTracksArray->At(irectracks2);
	 
	  trackParam = rectrack->GetTrackParamAtVertex();
	  bendingSlope    = trackParam->GetBendingSlope();
	  nonBendingSlope = trackParam->GetNonBendingSlope();
	  
	  pYZ      = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
	  fPzRec2  = - pYZ / TMath::Sqrt(1.0 + bendingSlope * bendingSlope); // spectro. (z<0)
	  fPxRec2  = fPzRec2 * nonBendingSlope;
	  fPyRec2  = fPzRec2 * bendingSlope;
	  fZRec2   = trackParam->GetZ();
	  fCharge2 = Int_t(TMath::Sign(1., trackParam->GetInverseBendingMomentum()));

	  fE2 = TMath::Sqrt(muonMass * muonMass + fPxRec2 * fPxRec2 + fPyRec2 * fPyRec2 + fPzRec2 * fPzRec2);
	  fV2.SetPxPyPzE(fPxRec2, fPyRec2, fPzRec2, fE2);

	  ntrackhits = rectrack->GetNTrackHits();
	  fitfmin = rectrack->GetFitFMin();

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

	      if (invMass > massMin && invMass < massMax) {
		EventInMass++;
  		hRapResonance->Fill(fVtot.Rapidity());
  		hPtResonance->Fill(fVtot.Pt());
	      }

	    } // if (fCharge * fCharge2) == -1)
	  } // if ((ch2 < Chi2Cut) && (pt2 > PtCutMin) && (pt2 < PtCutMax))
	} //  for (Int_t irectracks2 = irectracks + 1; irectracks2 < irectracks; irectracks2++)
      } // if (ch1 < Chi2Cut) && (pt1 > PtCutMin)&& (pt1 < PtCutMax) )
    } // for (Int_t irectracks = 0; irectracks < nrectracks; irectracks++)

    hNumberOfTrack->Fill(nrectracks);
  } // for (Int_t ievent = FirstEvent;

  histoFile->Write();
  histoFile->Close();

  cout << "MUONmassPlot " << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "ResType " << ResType << endl;
  cout << "Chi2Cut " << Chi2Cut << endl;
  cout << "PtCutMin " << PtCutMin << endl;
  cout << "PtCutMax " << PtCutMax << endl;
  cout << "massMin " << massMin << endl;
  cout << "massMax " << massMax << endl;
  cout << "EventInMass " << EventInMass << endl;
}

