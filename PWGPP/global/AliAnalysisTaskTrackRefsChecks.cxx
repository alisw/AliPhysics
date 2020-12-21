#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskTrackRefsChecks.h"
#include <TTree.h>


/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskTrackRefsChecks
// Task to access the TrackRefs
// 
//
// Authors: C. Zampolli
//          
//          
//*************************************************************************

ClassImp(AliAnalysisTaskTrackRefsChecks)
//______________________________________________________________________________
AliAnalysisTaskTrackRefsChecks::AliAnalysisTaskTrackRefsChecks() : 
  AliAnalysisTaskSE("TrackRefsChecks"), 
  fOutput(0x0),

  fhdtgl_pos(0x0),
  fhdtgl_prod_pos(0x0),

  fhdpt_pos(0x0),
  fhdpt_prod_pos(0x0),
  fhratiopt_pos(0x0),
  fhratiopt_prod_pos(0x0),

  fhdp_pos(0x0),
  fhdp_prod_pos(0x0),
  fhratiop_pos(0x0),
  fhratiop_prod_pos(0x0),

  fhdpz_pos(0x0),
  fhdpz_prod_pos(0x0),
  fhratiopz_pos(0x0),
  fhratiopz_prod_pos(0x0),

  fhdtgl_neg(0x0),
  fhdtgl_prod_neg(0x0),

  fhdpt_neg(0x0),
  fhdpt_prod_neg(0x0),
  fhratiopt_neg(0x0),
  fhratiopt_prod_neg(0x0),

  fhdp_neg(0x0),
  fhdp_prod_neg(0x0),
  fhratiop_neg(0x0),
  fhratiop_prod_neg(0x0),

  fhdpz_neg(0x0),
  fhdpz_prod_neg(0x0),
  fhratiopz_neg(0x0),
  fhratiopz_prod_neg(0x0),

  fhNTrackRefsInITS(0x0)
{
  //

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}


//___________________________________________________________________________
AliAnalysisTaskTrackRefsChecks::~AliAnalysisTaskTrackRefsChecks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;

  if(fOutput && !fOutput->IsOwner()){
    for (Int_t i = 0; i < 6; i++){

      delete fhdtgl_pos[i];
      delete fhdtgl_prod_pos[i];

      delete fhdpt_pos[i];
      delete fhdpt_prod_pos[i];
      delete fhratiopt_pos[i];
      delete fhratiopt_prod_pos[i];

      delete fhdp_pos[i];
      delete fhdp_prod_pos[i];
      delete fhratiop_pos[i];
      delete fhratiop_prod_pos[i];

      delete fhdpz_pos[i];
      delete fhdpz_prod_pos[i];
      delete fhratiopz_pos[i];
      delete fhratiopz_prod_pos[i];

      delete fhdtgl_neg[i];
      delete fhdtgl_prod_neg[i];

      delete fhdpt_neg[i];
      delete fhdpt_prod_neg[i];
      delete fhratiopt_neg[i];
      delete fhratiopt_prod_neg[i];

      delete fhdp_neg[i];
      delete fhdp_prod_neg[i];
      delete fhratiop_neg[i];
      delete fhratiop_prod_neg[i];

      delete fhdpz_neg[i];
      delete fhdpz_prod_neg[i];
      delete fhratiopz_neg[i];
      delete fhratiopz_prod_neg[i];
    }
    
    delete [] fhdtgl_pos;
    delete [] fhdtgl_prod_pos;

    delete [] fhdpt_pos;
    delete [] fhdpt_prod_pos;
    delete [] fhratiopt_pos;
    delete [] fhratiopt_prod_pos;

    delete [] fhdp_pos;
    delete [] fhdp_prod_pos;
    delete [] fhratiop_pos;
    delete [] fhratiop_prod_pos;

    delete [] fhdpz_pos;
    delete [] fhdpz_prod_pos;
    delete [] fhratiopz_pos;
    delete [] fhratiopz_prod_pos;

    delete [] fhdtgl_neg;
    delete [] fhdtgl_prod_neg;

    delete [] fhdpt_neg;
    delete [] fhdpt_prod_neg;
    delete [] fhratiopt_neg;
    delete [] fhratiopt_prod_neg;

    delete [] fhdp_neg;
    delete [] fhdp_prod_neg;
    delete [] fhratiop_neg;
    delete [] fhratiop_prod_neg;

    delete [] fhdpz_neg;
    delete [] fhdpz_prod_neg;
    delete [] fhratiopz_neg;
    delete [] fhratiopz_prod_neg;

    delete fhNTrackRefsInITS;
    
  }
    
  delete fOutput;

}
 
//___________________________________________________________________________
void AliAnalysisTaskTrackRefsChecks::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fhdtgl_pos = new TH3F*[6];
  fhdtgl_prod_pos = new TH3F*[6];

  fhdpt_pos = new TH3F*[6];
  fhdpt_prod_pos = new TH3F*[6];
  fhratiopt_pos = new TH3F*[6];
  fhratiopt_prod_pos = new TH3F*[6];

  fhdp_pos = new TH3F*[6];
  fhdp_prod_pos = new TH3F*[6];
  fhratiop_pos = new TH3F*[6];
  fhratiop_prod_pos = new TH3F*[6];
  
  fhdpz_pos = new TH3F*[6];
  fhdpz_prod_pos = new TH3F*[6];
  fhratiopz_pos = new TH3F*[6];
  fhratiopz_prod_pos = new TH3F*[6];

  fhdtgl_neg = new TH3F*[6];
  fhdtgl_prod_neg = new TH3F*[6];

  fhdpt_neg = new TH3F*[6];
  fhdpt_prod_neg = new TH3F*[6];
  fhratiopt_neg = new TH3F*[6];
  fhratiopt_prod_neg = new TH3F*[6];

  fhdp_neg = new TH3F*[6];
  fhdp_prod_neg = new TH3F*[6];
  fhratiop_neg = new TH3F*[6];
  fhratiop_prod_neg = new TH3F*[6];
  
  fhdpz_neg = new TH3F*[6];
  fhdpz_prod_neg = new TH3F*[6];
  fhratiopz_neg = new TH3F*[6];
  fhratiopz_prod_neg = new TH3F*[6];

  for (Int_t i = 0; i < 6; i++){
    fhdtgl_pos[i] = new TH3F(Form("hdtgl%d_pos", i), Form("Tangent Lambda change %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.05, 0.05);
    fOutput->Add(fhdtgl_pos[i]);

    fhdtgl_prod_pos[i] = new TH3F(Form("hdtgl%d_prod_pos", i), Form("Tangent Lambda change %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.05, 0.05);
    fOutput->Add(fhdtgl_prod_pos[i]);

    fhdpt_pos[i] = new TH3F(Form("hdpt%d_pos", i), Form("Pt difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpt_pos[i]);

    fhdpt_prod_pos[i] = new TH3F(Form("hdpt%d_prod_pos", i), Form("Pt difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpt_prod_pos[i]);

    fhratiopt_pos[i] = new TH3F(Form("hratiopt%d_pos", i), Form("Pt ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopt_pos[i]);

    fhratiopt_prod_pos[i] = new TH3F(Form("hratiopt%d_prod_pos", i), Form("Pt ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopt_prod_pos[i]);

    fhdp_pos[i] = new TH3F(Form("hdp%d_pos", i), Form("P difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdp_pos[i]);

    fhdp_prod_pos[i] = new TH3F(Form("hdp%d_prod_pos", i), Form("P difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdp_prod_pos[i]);

    fhratiop_pos[i] = new TH3F(Form("hratiop%d_pos", i), Form("P ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiop_pos[i]);

    fhratiop_prod_pos[i] = new TH3F(Form("hratiop%d_prod_pos", i), Form("P ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiop_prod_pos[i]);

    fhdpz_pos[i] = new TH3F(Form("hdpz%d_pos", i), Form("Pz difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpz_pos[i]);

    fhdpz_prod_pos[i] = new TH3F(Form("hdpz%d_prod_pos", i), Form("Pz difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpz_prod_pos[i]);

    fhratiopz_pos[i] = new TH3F(Form("hratiopz%d_pos", i), Form("Pz ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopz_pos[i]);

    fhratiopz_prod_pos[i] = new TH3F(Form("hratiopz%d_prod_pos", i), Form("Pz ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopz_prod_pos[i]);

    // negative tracks

    fhdtgl_neg[i] = new TH3F(Form("hdtgl%d_neg", i), Form("Tangent Lambda change %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.05, 0.05);
    fOutput->Add(fhdtgl_neg[i]);

    fhdtgl_prod_neg[i] = new TH3F(Form("hdtgl%d_prod_neg", i), Form("Tangent Lambda change %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.05, 0.05);
    fOutput->Add(fhdtgl_prod_neg[i]);

    fhdpt_neg[i] = new TH3F(Form("hdpt%d_neg", i), Form("Pt difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpt_neg[i]);

    fhdpt_prod_neg[i] = new TH3F(Form("hdpt%d_prod_neg", i), Form("Pt difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpt_prod_neg[i]);

    fhratiopt_neg[i] = new TH3F(Form("hratiopt%d_neg", i), Form("Pt ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopt_neg[i]);

    fhratiopt_prod_neg[i] = new TH3F(Form("hratiopt%d_prod_neg", i), Form("Pt ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopt_prod_neg[i]);

    fhdp_neg[i] = new TH3F(Form("hdp%d_neg", i), Form("P difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdp_neg[i]);

    fhdp_prod_neg[i] = new TH3F(Form("hdp%d_prod_neg", i), Form("P difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdp_prod_neg[i]);

    fhratiop_neg[i] = new TH3F(Form("hratiop%d_neg", i), Form("P ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiop_neg[i]);

    fhratiop_prod_neg[i] = new TH3F(Form("hratiop%d_prod_neg", i), Form("P ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiop_prod_neg[i]);

    fhdpz_neg[i] = new TH3F(Form("hdpz%d_neg", i), Form("Pz difference %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpz_neg[i]);

    fhdpz_prod_neg[i] = new TH3F(Form("hdpz%d_prod_neg", i), Form("Pz difference %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 40, -0.02, 0.02);
    fOutput->Add(fhdpz_prod_neg[i]);

    fhratiopz_neg[i] = new TH3F(Form("hratiopz%d_neg", i), Form("Pz ratio %d vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopz_neg[i]);

    fhratiopz_prod_neg[i] = new TH3F(Form("hratiopz%d_prod_neg", i), Form("Pz ratio %d wrt production value vs 1/pT vs phi", i), 400, 0., 2*TMath::Pi(), 10, 0, 4, 20, 0.95, 1.05);
    fOutput->Add(fhratiopz_prod_neg[i]);

  }

  fhNTrackRefsInITS = new TH1I("fhNTrackRefsInITS", "Number of track refs in ITS", 20, -0.5, 19.5);
  fOutput->Add(fhNTrackRefsInITS);
  
  PostData(1, fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskTrackRefsChecks::UserExec(Option_t *)
{
  //

  //Float_t r[6] = { 2.94, 3.9, 7.6, 11.04, 15.0, 23.9, 29.44 ,38.0, 43.0}; // beam pipe, spd, spd, sdd shield, sdd, sdd, sdd shield, ssd, ssd
  Float_t rmin[6] = { 3.2, 6., 12., 22., 35., 40.};
  Float_t rmax[6] = {4.5, 8., 18., 28., 39., 45.};
  
  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskTrackRefsChecks::UserExec(): bad ESD\n");
    return;
  } 

  /*
if (!fMCEvent) {
    AliFatal("NO MC INFO FOUND");
    return;
  }
  */
  
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  	     
  AliMCEvent* mcevent = mcinfo->MCEvent();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  
  // loop on the MC Generated Particles
  //Printf("There are %d particles", mcevent->GetNumberOfTracks());
  for (Int_t ipart=0; ipart < mcevent->GetNumberOfTracks(); ipart++) {

    AliMCParticle *mcPart = (AliMCParticle*)mcevent->GetTrack(ipart);
    if (!mcevent->IsPhysicalPrimary(TMath::Abs(mcPart->GetLabel()))) continue;
    if (TMath::Abs(mcPart->PdgCode()) != 211) continue;
    Bool_t isPositive = kTRUE;
    Short_t charge = mcPart->Charge();
    if (charge < 0) isPositive = kFALSE;
    Float_t ptMC = mcPart->Pt();
    Float_t pzMC = mcPart->Pz();
    Float_t pMC = mcPart->P();
    Float_t tglMC = pzMC/ptMC;
    Int_t nTrackRefs = mcPart->GetNumberOfTrackReferences();
    //Printf("\n\n\n=====================");
    //Printf("\n\n\nparticle %d is physical primary and pion and has %d track refs; original momentum pt = %f, pz = %f", ipart, nTrackRefs, ptMC, pzMC);
    Int_t nTrackRefsInITS=0;
    for(Int_t iTrackRef = 0; iTrackRef < nTrackRefs; iTrackRef++) {
      AliTrackReference * trackRef0 = mcPart->GetTrackReference(iTrackRef);
      if (trackRef0) {
	Int_t detID0 = trackRef0->DetectorId();
	Float_t r0 = trackRef0->R();
	Float_t pt0 = trackRef0->Pt(); 
	Float_t pz0 = trackRef0->Pz(); 
	Float_t p0 = trackRef0->P(); 
	Float_t phi0 = trackRef0->Phi();
	//Printf("\nTrack ref %d, R = %f, pt0 = %f, pz0 = %f, detID: %d", iTrackRef, r0, pt0, pz0, detID0);
	if (detID0 == -1){
	  //Printf("Track ref #1 with DetID = -1 --> skipping");
	  continue;
	}
	if (detID0 > 0){
	  //Printf("Track ref #1 not in ITS --> breaking");
	  break;
	}
	if (pt0 == 0.){
	  //Printf("Track ref #1 with pt = 0 --> skipping");
	  continue;
	}
	nTrackRefsInITS++;
	Float_t tgl0 = pz0/pt0;
	for (Int_t ilayer = 0; ilayer < 6; ilayer++){
	  if (r0 > rmin[ilayer] && r0 < rmax[ilayer]){ // track ref in layer ilayer
	    //Printf("Track ref is in   %d-th layer between %f and %f", ilayer, rmin[ilayer], rmax[ilayer]);
	    //if (ilayer < 5) //Printf("We will check the %d-th layer between %f and %f", ilayer+1, rmin[ilayer+1], rmax[ilayer+1]);
	    //else Printf("We will not check for further layers, but only fill the histo wrt MC values");
	    //Printf("Filling first histogram wrt production values");
	    Float_t dtglMC = tgl0 - tglMC;
	    if (isPositive){
	      fhdtgl_prod_pos[ilayer]->Fill(phi0, 1./pt0, dtglMC);
 	      fhdpt_prod_pos[ilayer]->Fill(phi0, 1./pt0, pt0 - ptMC);
 	      fhratiopt_prod_pos[ilayer]->Fill(phi0, 1./pt0, pt0/ptMC);
 	      fhdp_prod_pos[ilayer]->Fill(phi0, 1./pt0, p0 - pMC);
 	      fhratiop_prod_pos[ilayer]->Fill(phi0, 1./pt0, p0/pMC);
 	      fhdpz_prod_pos[ilayer]->Fill(phi0, 1./pt0, pz0 - pzMC);
 	      fhratiopz_prod_pos[ilayer]->Fill(phi0, 1./pt0, pz0/pzMC);
	    }
	    else {
	      fhdtgl_prod_neg[ilayer]->Fill(phi0, 1./pt0, dtglMC);
	      fhdpt_prod_neg[ilayer]->Fill(phi0, 1./pt0, pt0 - ptMC);
	      fhratiopt_prod_neg[ilayer]->Fill(phi0, 1./pt0, pt0/ptMC);
	      fhdp_prod_neg[ilayer]->Fill(phi0, 1./pt0, p0 - pMC);
	      fhratiop_prod_neg[ilayer]->Fill(phi0, 1./pt0, p0/pMC);
 	      fhdpz_prod_neg[ilayer]->Fill(phi0, 1./pt0, pz0 - pzMC);
 	      fhratiopz_prod_neg[ilayer]->Fill(phi0, 1./pt0, pz0/pzMC);
	    }
	      // compare to following track ref only if you have layers left
	    if (ilayer < 5){
	      for (Int_t iTrackRef1 = iTrackRef+1; iTrackRef1 < nTrackRefs; iTrackRef1++){
		AliTrackReference * trackRef1 = mcPart->GetTrackReference(iTrackRef1);
		if (trackRef1) {
		  // check that it is a trackRef on the subsequent layer
		  Float_t r1 = trackRef1->R();
		  Int_t detID1 = trackRef1->DetectorId();
		  Float_t pt1 = trackRef1->Pt(); 
		  Float_t pz1 = trackRef1->Pz(); 
		  Float_t p1 = trackRef1->P(); 
		  //Printf("r1 = %f, detID1 = %d", r1, detID1);
		  if (detID1 == -1){
		    //Printf("Track ref #2 with DetID = -1 --> breaking");
		    break;
		  }
		  if (detID1 > 0){
		    //Printf("Track ref #2 outside ITS --> breaking");
		    break;
		  }
		  if (pt1 == 0.){
		    //Printf("Track ref #2 with pt = 0 --> skipping");
		    continue;
		  }
		  if (r1 > rmin[ilayer+1] && r1 < rmax[ilayer+1]){ // consider only subsequent layers
		    //Printf("Second track ref is in the layer we are looking at");
		    //Printf("pt0 = %f, pz0 = %f, pt1 = %f, pz1 = %f", pt0, pz0, pt1, pz1);
		    Float_t tgl1 = pz1/pt1;
		    Float_t dtgl = tgl1-tgl0;
		    //Printf("FILLING with %f", dtgl);
		    if (isPositive){
		      fhdtgl_pos[ilayer]->Fill(phi0, 1./pt0, dtgl);
		      fhdpt_pos[ilayer]->Fill(phi0, 1./pt0, pt1-pt0);
		      fhratiopt_pos[ilayer]->Fill(phi0, 1./pt0, pt1/pt0);
		      fhdp_pos[ilayer]->Fill(phi0, 1./pt0, p1-p0);
		      fhratiop_pos[ilayer]->Fill(phi0, 1./pt0, p1/p0);
		      fhdpz_pos[ilayer]->Fill(phi0, 1./pt0, pz1-pz0);
		      fhratiopz_pos[ilayer]->Fill(phi0, 1./pt0, pz1/pz0);
		    }
		    else {
		      fhdtgl_neg[ilayer]->Fill(phi0, 1./pt0, dtgl);
		      fhdpt_neg[ilayer]->Fill(phi0, 1./pt0, pt1-pt0);
		      fhratiopt_neg[ilayer]->Fill(phi0, 1./pt0, pt1/pt0);
		      fhdp_neg[ilayer]->Fill(phi0, 1./pt0, p1-p0);
		      fhratiop_neg[ilayer]->Fill(phi0, 1./pt0, p1/p0);
		      fhdpz_neg[ilayer]->Fill(phi0, 1./pt0, pz1-pz0);
		      fhratiopz_neg[ilayer]->Fill(phi0, 1./pt0, pz1/pz0);
		    }
		    break;
		  }
		  else {
		    //Printf("Second track ref is NOT in the layer we are looking at");
		    if (r1 > rmax[ilayer+1]){
		      //Printf("It is not worth to check further layers, as they are even more to the outside");
		      break;
		    }
		    if (r1 > rmin[ilayer] && r1 < rmax[ilayer]){
		      //Printf("We are still in the same layer, not using this track ref");
		      continue;
		    }
		  }
		}
	      }
	    }
	    break;
	  }
	}
      }
    }
    fhNTrackRefsInITS->Fill(nTrackRefsInITS);
  } 

      
    /*
  Int_t ntracks = esd->GetNumberOfTracks();

  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    
    AliESDtrack * track = esd->GetTrack(iTrack);
    if (!track) continue;
    fITSncls->Fill(track->GetITSNcls());
    fITSClusterMap->Fill(track->GetITSClusterMap());
  }
    */

    PostData(1,fOutput);

   
}
//______________________________________________________________________________
void AliAnalysisTaskTrackRefsChecks::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  return;
}





