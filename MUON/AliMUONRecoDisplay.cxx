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

/*
$Log$
Revision 1.3  2000/12/21 17:51:54  morsch
RN3 violations corrected

Revision 1.2  2000/11/23 10:09:39  gosset
Bug correction in AliMUONRecoDisplay.
Copyright, $Log$
Copyright, Revision 1.3  2000/12/21 17:51:54  morsch
Copyright, RN3 violations corrected
Copyright,, $Id$, comments at the right place for automatic documentation,
in AliMUONRecoEvent and AliMUONRecoDisplay

*/

//Authors: Mihaela Gheata, Andrei Gheata 09/10/00
//////////////////////////////////////////////////////////////////////
//                                                                  //
// AliMUONRecoDisplay						    //
//								    //
// This class subclasses AliDisplay and provides display of         //
// reconstructed tracks with following functionality : 		    //
//	- front/top/side/3D display of MUON reconstructed tracks    //
//        as polylines ;                                            //
//	- context menu activated when the main pad is right-clicked //
//	The context menu contains following functions :		    //
//	* SetDrawHits()	- switches on or off Geant hits ;	    //
//	* CutMomentum()	- displays only tracks within Pmin - Pmax   //
//	* ListTracks()	- prints ID and momentum info. for all	    //
//	tracks within momentum range Pmin,Pmax ;		    //
//	* Highlight()	- shows only one selected reco. track	    //
//	and its best matching Geant track;			    //
//	* UnHighlight()	- self explaining;			    //
//	* RecoEfficiency() - compute reco. efficiency for all events//
//        from galice.root file; also fake track percentage; make   //
//        plots for momentum precision                              //
//      * XYPlot()      - make X-Y plots of reconstructed and       //
//        generated tracks in all chambers                          //
//								    //
//      Starting : generate and reconstruct events, then use the    //
//                 MUONrecodisplay.C macro                          //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <AliRun.h>
#include <TClonesArray.h>
#include "AliMUONRecoEvent.h"
#include "AliMUONRecoDisplay.h"
#include <TROOT.h>
#include <AliPoints.h>
#include <TSlider.h>
#include <TView.h>
#include <TGeometry.h>

ClassImp(AliMUONRecoDisplay)

//-------------------------------------------------------------------
AliMUONRecoDisplay::AliMUONRecoDisplay(Int_t nevent)
                  :AliDisplay(750)
{
//************ Constructor of the reco. event display**********
   // get reconstructed event from file
   fFile = new TFile("tree_reco.root");
   if (!fFile) {
      cout << "File tree_reco.root not found\n";
      gApplication->Terminate(0);
   }
   fEvReco = 0;
   fTree = (TTree *) fFile->Get("TreeRecoEvent");
   if (!fTree) {
      cout << "Tree of reconstructed events not found on file. Abort.\n";
      gApplication->Terminate(0);
   }
   fEvent = nevent;
   fEmpty = kFALSE;
   TBranch *branch = fTree->GetBranch("Event");
   branch->SetAddress(&fEvReco);
   
   fEvGen  = new AliMUONRecoEvent();

   TFile *galiceFile = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   galiceFile->cd();
   if (nevent > gAlice->TreeE()->GetEntries() - 1) {
      cout << "Event number out of range !\n";
      gApplication->Terminate(0);
   }
   gAlice->GetEvent(nevent);


   fRecoTracks = 0;
   fGenTracks  = 0;
   fPolyRecoList = 0;
   fPolyGenList  = 0;
   fHighlited   = -1;
   fMinMomentum = 0;
   fMaxMomentum = 999;
   // map this event
   MapEvent(nevent);
}

//-------------------------------------------------------------------
AliMUONRecoDisplay::~AliMUONRecoDisplay()
{
// Destructor of display object
   if (fPolyRecoList) {
      fPolyRecoList->Delete();
      delete fPolyRecoList;
   }
   if (fPolyGenList) {
      fPolyGenList->Delete();
      delete fPolyGenList;
   }
   delete fEvGen;
}
//-------------------------------------------------------------------
Bool_t AliMUONRecoDisplay::Event(Int_t nevent)
{
// Go to event nevent
   fEvent = nevent;
   for (Int_t entry=0; entry<fTree->GetEntries(); entry++) {
      fTree->GetEntry(entry);
      if (fEvReco->GetNoEvent() == nevent) return kTRUE;
   }
   cout << "Event number " << nevent << " empty\n";
   fRecoTracks = 0;
   fPolyRecoList = 0;
   fGenTracks = 0;
   fPolyGenList = 0;
   fHighlited   = -1;
   fMinMomentum = 0;
   fMaxMomentum = 999;
   
   return kFALSE;
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::MapEvent(Int_t nevent)
{
// get generated event (nevent) from galice.root and corresponding
// reconstructed event from tree_reco.root; 
   cout << "mapping event " << nevent << endl;
   fEvent = nevent;
   fEmpty = kFALSE;
   fFile->cd();
   // check if the event is not empty and make fEvReco to point to it
   fEmpty = !Event(nevent);
   // just testing
//   if (!fEvReco) {
//      cout << "failed to load reco event ! Correct this.\n";
//      gApplication->Terminate(0);
//   }
   if (fEmpty) return;
   fRecoTracks   = fEvReco->TracksPtr();
   fPolyRecoList = MakePolyLines3D(fRecoTracks);   
   // clear previous event
   if (fEvGen) fEvGen->Clear();
   fEvGen->SetNoEvent(nevent);
   // get list of particles
   TClonesArray *particles = gAlice->Particles();
   // connect MUON module
   AliDetector *pMUON = gAlice->GetDetector("MUON");
   if (!pMUON) {
      cout << "MUON module not present.\n";
      gApplication->Terminate(0);
   }
   // get the number of generated tracks
   Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
   // Fill the fEvGen object
   AliMUONRecoTrack *gtrack = 0;
   AliMUONHit *hit = 0;
   TParticle *particle;
   Int_t ch;
   // loop all tracks
   for (Int_t track=0; track<ntracks; track++) {
      hit = (AliMUONHit *) pMUON->FirstHit(track);
      if (!hit) continue;
      particle = (TParticle *) particles->UncheckedAt(hit->Track());
      if (IsReconstructible(track) && TMath::Abs(particle->GetPdgCode())==13) {
         gtrack = fEvGen->AddEmptyTrack();
	 gtrack->SetSign(TMath::Sign((Int_t)1, -particle->GetPdgCode()));
	 // reset hits
	 for (ch=0; ch<10; ch++) gtrack->SetHitPosition(ch,0,0,0);
	 // loop all hits
	 for (AliMUONHit *muonHit=(AliMUONHit*)pMUON->FirstHit(track);
	      muonHit;
	      muonHit=(AliMUONHit*)pMUON->NextHit()) {
	    ch = muonHit->fChamber - 1;
	    if (ch<0 || ch>9) continue;
	    gtrack->SetHitPosition(ch, muonHit->X(),  muonHit->Y(),  muonHit->Z());
	    gtrack->SetMomReconstr(particle->Px(), particle->Py(), particle->Pz());
	    gtrack->SetVertexPos(particle->Vz());
	    gtrack->SetChi2r(0.0);
	 } 
      }
   }
   fGenTracks   = fEvGen->TracksPtr();
   fPolyGenList = MakePolyLines3D(fGenTracks);
}
//-------------------------------------------------------------------
void AliMUONRecoDisplay::XYPlot()
{
// Plot reco. tracks hits in all chambers for current event:
//  - open blue squares : generated muons that were reconstructed accurate
//  - open cyan squares : generated muons that were not reconstructed
//  - filled green circles : reco. tracks (accurate)
//  - filled red squares   : fake tracks 
   if (fEmpty) return;
   Double_t kMaxRadius[10];
   kMaxRadius[0] =  kMaxRadius[1] = 91.5;
   kMaxRadius[2] =  kMaxRadius[3] = 122.5;
   kMaxRadius[4] =  kMaxRadius[5] = 158.3;
   kMaxRadius[6] =  kMaxRadius[7] = 260.0;
   kMaxRadius[8] =  kMaxRadius[9] = 260.0;

   TH2F *xygenFound[10]; 
   TH2F *xygenLost[10]; 
   TH2F *xyrecoGood[10];
   TH2F *xyrecoFake[10];
   Double_t x,y,r;
   Int_t matches[500];
   Int_t index, ch;

   TPad *pad = (TPad*)gROOT->GetSelectedPad();
   TCanvas *canvas = new TCanvas("xy", "Reconstruction efficiency");
   canvas->Clear();
   canvas->cd();
   canvas->Divide(3,4);
   canvas->cd(11);
   gPad->Delete();
   canvas->cd(12);
   gPad->Delete();
   canvas->cd(1);
   
   // Define histograms for x-y plots
   for (ch=0; ch<10; ch++) {
      xygenFound[ch] = new TH2F("xygen_found","",50,-kMaxRadius[ch],kMaxRadius[ch],50,-kMaxRadius[ch],kMaxRadius[ch]);
      xygenLost[ch] = new TH2F("xygen_lost","",50,-kMaxRadius[ch],kMaxRadius[ch],50,-kMaxRadius[ch],kMaxRadius[ch]);
      xyrecoGood[ch] = new TH2F("xyreco_good","",50,-kMaxRadius[ch],kMaxRadius[ch],50,-kMaxRadius[ch],kMaxRadius[ch]);
      xyrecoFake[ch] = new TH2F("xyreco_fake","",50,-kMaxRadius[ch],kMaxRadius[ch],50,-kMaxRadius[ch],kMaxRadius[ch]);
   }
   // find list of matching tracks
   fPrinted = kTRUE; // no need to print
   for (index=0; index<fRecoTracks->GetEntriesFast(); index++) {
      matches[index] = GetBestMatch(index);  
   }
   // and fill them
   for (index=0; index<fGenTracks->GetEntriesFast(); index++) {	// GEANT
      Bool_t wasreconst = kFALSE;
      for (Int_t i=0; i<fRecoTracks->GetEntriesFast(); i++) {
         if (matches[i] == index) {
	    wasreconst = kTRUE;
	    break;
	 }
      }
      AliMUONRecoTrack *current = (AliMUONRecoTrack*)fGenTracks->UncheckedAt(index);
      for (ch=0; ch<10; ch++) {
         x = current->GetPosX(ch);
         y = current->GetPosY(ch);
         r = TMath::Sqrt(x*x +y*y);
         if (r >= 10) {
            if (wasreconst) {
               xygenFound[ch]->Fill(x,y);
            } else {xygenLost[ch]->Fill(x,y);}
         }
      }
   }
   for (index=0; index<fRecoTracks->GetEntriesFast(); index++) {	// reco
      AliMUONRecoTrack *current = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(index);
      for (ch=0; ch<10; ch++) {
         x = current->GetPosX(ch);
         y = current->GetPosY(ch);
         r = TMath::Sqrt(x*x +y*y);
         if (r >= 10) {
            if (matches[index] >= 0) {
               xyrecoGood[ch]->Fill(x,y);
            } else {xyrecoFake[ch]->Fill(x,y);}
         }
      }
   }
   // finally plot them
   for (ch=0; ch<10; ch++) {
      canvas->cd(ch+1);   
      xygenFound[ch]->SetMarkerColor(kBlue);
      xygenFound[ch]->SetMarkerStyle(4);
      xygenFound[ch]->SetMarkerSize(0.5);
      xygenFound[ch]->SetStats(kFALSE);
      xygenFound[ch]->Draw();
      xygenLost[ch]->SetMarkerColor(kCyan);
      xygenLost[ch]->SetMarkerStyle(4);
      xygenLost[ch]->SetMarkerSize(0.5);
      xygenLost[ch]->SetStats(kFALSE);
      xygenLost[ch]->Draw("SAME");
      xyrecoGood[ch]->SetMarkerColor(kGreen);
      xyrecoGood[ch]->SetMarkerStyle(20);
      xyrecoGood[ch]->SetMarkerSize(0.4);
      xyrecoGood[ch]->SetStats(kFALSE);
      xyrecoGood[ch]->Draw("SAME");
      xyrecoFake[ch]->SetMarkerColor(kRed);
      xyrecoFake[ch]->SetMarkerStyle(20);
      xyrecoFake[ch]->SetMarkerSize(0.5);
      xyrecoFake[ch]->SetStats(kFALSE);
      xyrecoFake[ch]->Draw("SAME");
   }
   canvas->SetTitle("y vs. x for simulated and reconstructed tracks");
   pad->cd();
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::RecoEfficiency(Int_t first, Int_t last)
{
// Loop selected reconstructed events, compute efficiency as total number of 
// reconstructed tracks (accurate) over total number of generated
// reconstructible muons. Make histogram for momentum precision and profile
// of mean momentum precision versus total momentum (generated)
   Int_t nevents = (Int_t)gAlice->TreeE()->GetEntries();
   if (last > nevents) last = nevents - 1;
   if (first < 0) first = 0;
   nevents = last - first + 1;
   Int_t track;
   Float_t efficiency;
   Float_t generated=0, found=0, fake=0; // number of generated/found/fake tracks
   Double_t pgen, preco, dP;             // dP=preco-pgen
   AliMUONRecoTrack *rtrack, *gtrack;    // generated/reco. current tracks
   fPrinted = kTRUE;                     // no need to print
   Int_t currentEvent = gAlice->GetHeader()->GetEvent();
   
   TH1F *pReso = new TH1F("Momentum precision", "dP = Prec - Pgen", 100, -5.0, 5.0);
   pReso->SetXTitle("dP [GeV/c]");
   pReso->SetYTitle("dN/dP");
   
   TProfile *dPp = new TProfile("", "dP vs. P", 50, 0, 100, 0, 5);
   dPp->SetXTitle("P [GeV/c]");
   dPp->SetYTitle("<dP> [GeV/c]"); 

   TPad *pad = (TPad*)gROOT->GetSelectedPad();
   TCanvas *canvas = new TCanvas("c", "Reconstruction efficiency");
   canvas->Divide(1,2);
   canvas->cd(1);

   // loop events
   for (Int_t event=first; event<=last; event++) {
      // get the reco. & gen. events
      TFile *galiceFile = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
      galiceFile->cd();
      gAlice->GetEvent(event);
      MapEvent(event);
      if (fEmpty) {
      // skip empty events
         fEmpty = kFALSE;
	 continue;
      }
      generated += fGenTracks->GetEntriesFast();
      // loop reco. tracks
      for (track=0; track<fRecoTracks->GetEntriesFast(); track++) {
         rtrack = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(track);
         preco  = rtrack->P();
	 // find correspondent generated track
	 Int_t ind = GetBestMatch(track);
	 if (ind >=0) {
	    gtrack = (AliMUONRecoTrack*)fGenTracks->UncheckedAt(ind);
	    pgen  = gtrack->P();
	    found += 1;
	    dP = preco - pgen;
	    pReso->Fill(dP);
	    dPp->Fill(pgen, TMath::Abs(dP));
	 } else {
	    fake += 1;
	 }
      }
   }
   efficiency = found/generated;
   cout << "=================================================================\n";
   cout << "||       Reconstruction efficiency and momentum precision      ||\n";
   cout << "=================================================================\n";
   cout << "Number of events processed              " << nevents << endl;
   cout << "Total number of reconstructible muons : " << (Int_t)generated << endl;
   cout << "RECONSTRUCTION EFFICIENCY :             " << efficiency*100 << " %" << endl;
   cout << "Fake track rate :                       " << 100*fake/generated << " %" << endl;
   cout << "Momentum precision fit : \n" << endl;
   pReso->Fit("gaus");
   cout << "=================================================================\n";    
   
   canvas->cd(1);
   pReso->Draw();
   canvas->cd(2);
   dPp->Draw("HIST");
   pad->cd();
   ShowNextEvent(currentEvent-last);
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::ShowNextEvent(Int_t delta)
{
// overwritten from AliDisplay in order to get also mapping of next event
   TFile *galiceFile = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   galiceFile->cd();
   if (delta) {
      gAlice->Clear();
      Int_t currentEvent = gAlice->GetHeader()->GetEvent();
      Int_t newEvent     = currentEvent + delta;
      if (newEvent<0 || newEvent>(gAlice->TreeE()->GetEntries() - 1)) return;
      Int_t nparticles = gAlice->GetEvent(newEvent);
      cout << "Event : " << newEvent << " with " << nparticles << " particles\n";
      if (!gAlice->TreeH()) return;
      MapEvent(newEvent);
      fHighlited = -1;
   }
   LoadPoints();
   fPad->cd();
   Draw();
   if (gROOT->GetListOfCanvases()->FindObject("xy")) XYPlot();
}
//-------------------------------------------------------------------
Bool_t AliMUONRecoDisplay::IsReconstructible(Int_t track)
{
// true if at least three hits in first 2 stations, 3 in last 2 stations
// and one in station 3
   if (fEmpty) return kFALSE;
   AliDetector *pMUON = gAlice->GetDetector("MUON");
   Bool_t chHit[10];
   Int_t ch;
   for (ch=0; ch<10; ch++) chHit[ch] = kFALSE;
   //loop hits
   for (AliMUONHit *muonHit=(AliMUONHit*)pMUON->FirstHit(track);
        muonHit;
	muonHit=(AliMUONHit*)pMUON->NextHit()) {
      ch = muonHit->fChamber - 1;
      if (ch<0 || ch>9) continue;
      chHit[ch] = kTRUE;
   }
   Int_t nhits = 0;
   for (ch=0; ch<4; ch++) nhits += (chHit[ch])?1:0;
   if (nhits < 3) return kFALSE;
   nhits = 0;
   for (ch=4; ch<6; ch++) nhits+= (chHit[ch])?1:0;
   if (nhits < 1) return kFALSE;
   
   for (ch=7; ch<10; ch++) nhits+= (chHit[ch])?1:0;
   if (nhits < 3) return kFALSE;
   return kTRUE;
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::DrawView(Float_t theta, Float_t phi, Float_t psi)
{
// ovewritten from base class to change the range for MUON
   gPad->SetCursor(kWatch);
   gPad->SetFillColor(1);
   gPad->Clear();
   
   Int_t iret;
   TView *view = new TView(1);
   Float_t range = fRrange*fRangeSlider->GetMaximum();
   view->SetRange(-range, -range, 0, range, range, 5*range);
   fZoomX0[0] = -1;
   fZoomY0[0] = -1;
   fZoomX1[0] =  1;
   fZoomY1[0] =  1;
   // Display Alice geometry
   gAlice->GetGeometry()->Draw("same");
   //Loop on all detectors to add their products to the pad
   DrawHits();
   // add itself to the list (last)
   AppendPad();
   
   view->SetView(phi, theta, psi, iret);
}
//-------------------------------------------------------------------
void AliMUONRecoDisplay::SetDrawHits(Bool_t hits)
{
// Turns on/off Geant hits drawing
   fDrawHits = hits;
   DrawView(0,-90);
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::CutMomentum(Double_t min, Double_t max)
{
// Define momentum cut for displayed reconstructed tracks
   fMinMomentum = min;
   fMaxMomentum = max;
   if (fHighlited >= 0) UnHighlight();
   DrawView(0,-90);
}

//-------------------------------------------------------------------
Int_t AliMUONRecoDisplay::GetBestMatch(Int_t indr, Float_t tolerance)
{
// Find the index of best Geant track matching a reconstructed track : track
//	with maximum number of compatible hits (within tolerance*sigma bending and
//	non-bending resolution) and minimum number of fake hits;
// If no match is found within a given tolerance, the method is called recursively
//      with increasing tolerance, until tolerance = 10;
   if (fEmpty) return -1;
   if (indr<0 || indr>=fRecoTracks->GetEntriesFast()) return -1;
   AliMUONRecoTrack *rtrack = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(indr);
   AliMUONRecoTrack *gtrack = 0;
   Int_t bestMatch = -1;
   Int_t maxNcompat = 0;
   Int_t minNfakes = 10;
   Double_t xrhit,yrhit,radius,xghit,yghit,dX,dY;
// loop over all Geant tracks
   for (Int_t indg=0; indg<fGenTracks->GetEntriesFast(); indg++) {
      gtrack = (AliMUONRecoTrack*)fGenTracks->UncheckedAt(indg);
      if (!gtrack) continue;
      Int_t ncompat = 0;      // number of compat. hits for this track
      Int_t nfake = 0;        // number of fakes
      // loop chambers to find compatible hits
      for (Int_t ch=0; ch<10; ch++) {
         xrhit = rtrack->GetPosX(ch);
         yrhit = rtrack->GetPosY(ch);
         radius = TMath::Sqrt(xrhit*xrhit + yrhit*yrhit);
         if (radius<10) continue; // skip null hits
         xghit = gtrack->GetPosX(ch);
         yghit = gtrack->GetPosY(ch);
         dX = TMath::Abs(xghit-xrhit);
         dY = TMath::Abs(yghit-yrhit);
         if (dX<tolerance*0.144 && dY<tolerance*0.01) {// within tol*sigma resolution
            ncompat++;
            continue;      // compatible hit
         } else nfake++;      // fake hit
      }
      if (ncompat && ncompat>=maxNcompat && nfake<minNfakes) { // this is best matching
         maxNcompat = ncompat;
         minNfakes = nfake;
         bestMatch = indg;
      }
   }
   if (bestMatch<0 && tolerance<=9.) bestMatch = GetBestMatch(indr, tolerance+=1);
   if (!fPrinted) {
      rtrack = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(indr);
      Int_t sign = rtrack->GetSign();
      cout << "Reconstructed track : " << indr << "(" << sign << ")" << endl;
      rtrack->TrackInfo();
      printf("Best matching Geant track within %i*sgm : %i\n", (Int_t)tolerance, bestMatch);
      if (bestMatch >=0) {
         gtrack = (AliMUONRecoTrack*)fGenTracks->UncheckedAt(bestMatch);
         gtrack->TrackInfo();
      }
      cout << "-----------------------------------------------------------------\n";
   }
   fPrinted = kTRUE;

   return bestMatch;
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::Highlight(Int_t track)
{
// Highlight the specified track
   if (fEmpty) return; 
   if (fHighlited >=0) UnHighlight();
   if (track<0 || track>fPolyRecoList->GetEntries()) return;
   TPolyLine3D *line = (TPolyLine3D*)fPolyRecoList->UncheckedAt(track);
   line->SetLineColor(kYellow);
   line->SetLineWidth(1);
   fHighlited = track;
//   DrawView(15,-45,135);
   fPad->cd();
   Draw();
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::UnHighlight()
{
// Unhighlight a previous highlighted track
   if (fHighlited < 0 || fEmpty) return;      // nothing to do
   TPolyLine3D *line = (TPolyLine3D*)fPolyRecoList->UncheckedAt(fHighlited);
   line->SetLineColor(kRed);
   line->SetLineWidth(1);
   fHighlited = -1;
//   DrawView(0,-90);
   fPad->cd();
   Draw();
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::DrawHits()
{
//    Draw hits for all ALICE detectors. Overwrites the DrawHits() method of the
//	base class for reco. track drawing

   Float_t cutmin, cutmax, etamin, etamax, pmom, smin, smax, eta, theta, r;
   const Float_t kptcutmax = 2;
   const Float_t ketacutmax = 1.5;
   Float_t *pxyz;
   Int_t ntracks,track;
   TParticle *particle;
   TObjArray *points;
   AliPoints *pm;
      
   //Get cut slider
   smax   = fCutSlider->GetMaximum();
   smin   = fCutSlider->GetMinimum();
   cutmin = kptcutmax*smin;
   if (smax < 0.98) cutmax = kptcutmax*smax;
   else             cutmax = 100000;
   
   //Get eta slider
   smax   = fEtaSlider->GetMaximum();
   smin   = fEtaSlider->GetMinimum();
   etamin = ketacutmax*(2*smin-1);
   etamax = ketacutmax*(2*smax-1);
   if (smin < 0.02) etamin = -1000;
   if (smax > 0.98) etamax =  1000;
      
   TIter next(gAlice->Modules());
   AliModule *module;
   fHitsCuts = 0;
   if (fDrawHits) {
      // draw hits in all modules
      while((module = (AliModule*)next())) {
         if (!module->IsActive()) continue;
         points = module->Points();
         if (!points) continue;
         ntracks = points->GetEntriesFast();
         for (track=0;track<ntracks;track++) {
            pm = (AliPoints*)points->UncheckedAt(track);
            if (!pm) continue;
            particle = pm->GetParticle();
            if (!particle) continue;
            pmom = particle->P();
            if (pmom < cutmin) continue;
            if (pmom > cutmax) continue;
            // as a first approximation, take eta of first point
            pxyz  = pm->GetP();
            r     = TMath::Sqrt(pxyz[0]*pxyz[0] + pxyz[1]*pxyz[1]);
            theta = TMath::ATan2(r,TMath::Abs(pxyz[2]));
            if(theta) eta = -TMath::Log(TMath::Tan(0.5*theta)); else eta = 1e10;
            if (pxyz[2] < 0) eta = -eta;
            if (eta < etamin || eta > etamax) continue;
            pm->Draw();
            fHitsCuts += pm->GetN();
         }
      }
   }
   // draw reconstructed tracks
   if (fEmpty) return;
   TPolyLine3D *line, *gline;
   Int_t bestMatch;
   Double_t px,py,pz,p;
   AliMUONRecoTrack *rtrack;

   if (fHighlited >= 0) {
      line = (TPolyLine3D*)fPolyRecoList->UncheckedAt(fHighlited);
      line->Draw();
      fPrinted = kFALSE;
      bestMatch = GetBestMatch(fHighlited);
      if (bestMatch>=0) {
         gline = (TPolyLine3D*)fPolyGenList->UncheckedAt(bestMatch);
         gline->SetLineColor(kRed);
         gline->SetLineWidth(2);
         gline->SetLineStyle(2);
         gline->Draw();
      }
   } else {
      for (track=0; track<fPolyRecoList->GetEntries(); track++) {
         rtrack = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(track);
         px = rtrack->GetMomReconstr(0);
         py = rtrack->GetMomReconstr(1);
         pz = rtrack->GetMomReconstr(2);
         p  = rtrack->P();
         if (p>fMinMomentum && p<fMaxMomentum) {
            line = (TPolyLine3D*)fPolyRecoList->UncheckedAt(track);
            line->Draw();
         }
      }
   }
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::ListTracks()
{
// List momentum information of all reconstructed traccks within fPmin and fPmax
//	cuts, as well as their best matching Geant tracks
   if (fEmpty) return;
   cout << "================================================================\n";
   printf("Reconstructed tracks with momentum in range : %g , %g [GeV/c]\n",
         fMinMomentum, fMaxMomentum);
   cout << "----------------------------------------------------------------\n";
   AliMUONRecoTrack *rtrack;
   Double_t p;
   Int_t sign;
   for (Int_t ind=0; ind<fRecoTracks->GetEntries(); ind++) {
      rtrack = (AliMUONRecoTrack*)fRecoTracks->UncheckedAt(ind);
      p  = rtrack->P();
      if (p>fMinMomentum && p<fMaxMomentum) {
         fPrinted = kFALSE;
         GetBestMatch(ind);
         sign = rtrack->GetSign();
      }
   }      
   cout << "================================================================\n";
}

//-------------------------------------------------------------------
TClonesArray* AliMUONRecoDisplay::MakePolyLines3D(TClonesArray *tracklist)
{
// Makes the list of polylines3D corresponding to the list of tracks
   if (fEmpty) return 0;
   if (tracklist!=fRecoTracks && tracklist!=fGenTracks) return 0;
   Bool_t reco = (tracklist==fRecoTracks)?kTRUE:kFALSE;
   // make sure there is no other list in memory
   if (reco) {
      if (fPolyRecoList) {
         fPolyRecoList->Delete();
         delete fPolyRecoList;
         fPolyRecoList = 0;
      }
   } else {
      if (fPolyGenList) {
         fPolyGenList->Delete();
         delete fPolyGenList;
         fPolyGenList = 0;
      }      
   }
   if (!tracklist->GetEntries()) return 0;
   
   AliMUONRecoTrack* track = 0;
   TClonesArray *polyLines3D = new TClonesArray("TPolyLine3D",1000);
   TClonesArray &polylist = *polyLines3D;
   TPolyLine3D *polyline = 0;
   
   // loop all tracks
   for (Int_t i=0; i<tracklist->GetEntries(); i++) {
      track = (AliMUONRecoTrack*)tracklist->UncheckedAt(i);
      Int_t ch = 0;
      Double_t x,y,z,r;
      polyline = new(polylist[i]) TPolyLine3D(2,"");
      polyline->SetLineColor(kRed);
      polyline->SetLineWidth(1);
      polyline->SetNextPoint(0,0,track->GetVertexPos()); // vertex point
      // loop chambers to fill TPolyLine3D objects
      for (ch=0; ch<10; ch++) {
         x = track->GetPosX(ch);
         y = track->GetPosY(ch);
         r = TMath::Sqrt(x*x + y*y);
         if (r < 10) continue;
         z = track->GetPosZ(ch);
         polyline->SetNextPoint(x,y,z);
      }      
   }
   return polyLines3D;
}

//-------------------------------------------------------------------
void AliMUONRecoDisplay::PolyLineInfo(TClonesArray *line3Dlist)
{
// Prints information (x, y, z coordinates) for all constructed polylines
   if (fEmpty) return;
   if (line3Dlist) {
      TPolyLine3D *polyline = 0;
      for(Int_t trackIndex=0; trackIndex<line3Dlist->GetEntries(); trackIndex++) {
         polyline = (TPolyLine3D*)line3Dlist->UncheckedAt(trackIndex);
         polyline->ls();
         Float_t *pl = polyline->GetP();
         for (Int_t i=0; i<polyline->GetN() ;i++) {
            printf(" x[%d]=%g, y[%d]=%g, z[%d]=%g\n",i,pl[3*i],i,pl[3*i+1],i,pl[3*i+2]);
         }
      }
   }
}


