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

/// \ingroup macros
/// \file fastMUONSim.C
/// \brief An example how to do the fast simulation and storing
/// the muons that survive the reconstruction
///
/// \author A. De Falco, H. Woehri, INFN Cagliari, April 2007
///
/// An example how to do the fast simulation and storing
/// the muons that survive the reconstruction in 
/// AliMUONTrackLight and AliMUONPairLight objects for
/// further analysis. Should be used together with the
/// macro fastMUONGen.C

#if !defined(__CINT__) || defined(__MAKECINT__)
//STEER includes
#include "AliStack.h"
#include "AliRun.h"
#include "AliHeader.h"

//FASTSIM includes
#include "AliFastMuonTriggerEff.h"
#include "AliFastMuonTrackingAcc.h"
#include "AliFastMuonTrackingEff.h"
#include "AliFastMuonTrackingRes.h"
#include "AliFastDetector.h"

//MUON includes
#include "AliMUONTrackLight.h"
#include "AliMUONPairLight.h"

//ROOT includes
#include "TParticle.h"
#include "TRandom.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"

#endif

Float_t fThetaMin = 171., fThetaMax = 178.;
//====================================================
void fastMUONSim(Float_t bkg=0, char* outFileName = "fastSim_pp.root"){

  //setting up the fast simulator
  AliFastMuonTriggerEff *trigeff = new AliFastMuonTriggerEff();
  AliFastMuonTrackingAcc *acc = new AliFastMuonTrackingAcc();
  AliFastMuonTrackingEff *eff = new AliFastMuonTrackingEff();
  AliFastMuonTrackingRes *res = new AliFastMuonTrackingRes();
  acc->SetBackground(bkg);
  eff->SetBackground(bkg);
  res->SetBackground(bkg);  
  acc->Init(); 
  eff->Init(); 
  res->Init(); 
  trigeff->Init();

  trigeff->SetCut(0); // 0... trigger pt cut at 1 GeV
                      // 1... trigger pt cut at 2 GeV

  //prepare the arrays to store the generated/simulated muons
  TClonesArray *muonArray   = new TClonesArray("AliMUONTrackLight",100); 
  TClonesArray *dimuonArray = new TClonesArray("AliMUONPairLight",100); 
  TTree *treeOut = new TTree("tree","tree"); 
  TFile *fout = new TFile(outFileName,"recreate"); 
  fout->cd(); 
  treeOut->Branch("muons",&muonArray); 
  treeOut->Branch("dimuons",&dimuonArray); 
  //
  // Stack of particle for each event
  AliStack* stack;
  // Creating Run Loader and opening file containing Hits
  // AliRunLoader * RunLoader = AliRunLoader::Open("galice.root","MUONFolder","READ");
  AliRunLoader *RunLoader = AliRunLoader::Open("galice.root");
  if (RunLoader == 0x0) {
    printf(">>> Error : Error Opening %s file \n", "galice.root");
    return;
  }

  RunLoader->LoadKinematics("READ");
  Int_t nevents = RunLoader->GetNumberOfEvents();
  printf("nevents %d\n", nevents);

  TParticle *part = 0, *partMother=0;
  AliMUONTrackLight *muon = 0;
  Double_t radeg =  180./TMath::Pi(); 
  for(int iev = 0; iev < nevents; iev++){

    if(iev%1000==0) printf("Event %d\n", iev);
    RunLoader->GetEvent(iev); 
    stack = RunLoader->Stack();
    //stack->DumpPStack();

    //reset muon info
    Int_t nMuons = 0;
    muonArray->Clear();     // clean muon and dimuon arrays 
    dimuonArray->Clear(); 
    //
    Int_t nMuonsStored = 0; 
    Int_t nPart = stack->GetNtrack();
    for(int iPart = 0; iPart < nPart; iPart++){

      part = stack->Particle(iPart);
      if(fabs(part->GetPdgCode()) == 13){ //muon found

	nMuons++;
	Double_t xyz[3] = {part->Vx(),part->Vy(),part->Vz()};
	// do not take muons that were decayed after the absorber
	if(xyz[2] < -(90.+ 40.)) continue; //decay muons after 90cm + 1 int. length
	Int_t charge = (part->GetPdgCode() > 0) ? -1 : +1;
	Double_t px = part->Px();
	Double_t py = part->Py();
	Double_t pz = part->Pz();
	Double_t E = part->Energy();
	TLorentzVector pGen(px,py,pz,E);
	//
	//fast simulation code
	res->SetCharge(charge);
	eff->SetCharge(charge);
	acc->SetCharge(charge);
	Float_t p = (Float_t) pGen.P();
	Float_t pt = (Float_t) pGen.Pt();
	Float_t theta = (Float_t) radeg*pGen.Theta();
	Float_t phi = (Float_t) radeg*pGen.Phi(); 
	
	//continue only if we have a muon within the MUON acceptance
	if(theta < fThetaMin || theta > fThetaMax) continue;

	theta = 180. - theta; // correct by hand the 'new' coordinates frame
	phi   = 180. - phi; // correct by hand the 'new' coordinates frame
	Float_t prec = 0, thetarec = 0, phirec = 0; 
	res->Evaluate(p, theta ,phi, prec, thetarec, phirec);
	Float_t precx = prec * TMath::Sin(thetarec/radeg) * TMath::Cos(phirec/radeg); 
	Float_t precy = prec * TMath::Sin(thetarec/radeg) * TMath::Sin(phirec/radeg); 
	Float_t precz = prec * TMath::Cos(thetarec/radeg);
	precz = -precz; // correct by hand the 'new' coordinates frame
	precx   = -precx; // correct by hand the 'new' coordinates frame

  	Float_t trkeff = eff->Evaluate(charge, pt, theta, phi);
 	Float_t accep  = acc->Evaluate(charge, pt, theta, phi);
 	Float_t trgeff = trigeff->Evaluate(charge,pt,theta,phi);
	if(trkeff > 1.) printf("tracking efficiency > 1\n");
	if(accep > 1.) printf("acceptance efficiency > 1\n");
	if(trgeff > 1.) printf("trigger efficiency > 1\n");
	if (gRandom->Rndm() > trkeff || gRandom->Rndm() > accep) continue; 

	//only if we have a muon in the acceptance, store it as an
	//AliMUONTrackLight object and use it to form dimuons...
	muon = new AliMUONTrackLight();
	muon->SetCharge(charge);
	muon->SetTrackPDGCode(part->GetPdgCode()); //must be set, otherwise -999
	muon->SetTrackPythiaLine(iPart);//must be set, otherwise -999
	muon->SetPGen(pGen);
  	muon->FillMuonHistory(stack, part);
	// set vertex to mother's vx,vy,vz to have the primary vertex

	Int_t idMother = -1; 
	Int_t id2 = 0; 
	while (idMother < 0) idMother = muon->GetParentPythiaLine(id2++); 
	partMother = stack->Particle(idMother);
	Double_t xyzMother[3] = {partMother->Vx(),
				 partMother->Vy(),
				 partMother->Vz()};
	muon->SetVertex(xyzMother);
	if (gRandom->Rndm() > trgeff) muon->SetTriggered(kFALSE); 
	else muon->SetTriggered(kTRUE);
	muon->SetPxPyPz(precx,precy,precz);

	//store the referenced track in the muonArray:
	TClonesArray &muons = *muonArray;
	new (muons[nMuonsStored++]) AliMUONTrackLight(*muon);
	delete muon;
      }

    }//part

    Int_t nmuons = muonArray->GetEntriesFast(); 
    Int_t ndimuons = 0; 
    for(Int_t nMu1 = 0; nMu1 < nmuons-1; nMu1++) { 
      AliMUONTrackLight* mu1 = (AliMUONTrackLight*) muonArray->At(nMu1); 
      for(Int_t nMu2 = nMu1+1; nMu2 < nmuons; nMu2++){
	AliMUONTrackLight* mu2 = (AliMUONTrackLight*) muonArray->At(nMu2); 
	AliMUONPairLight *dimuLight = new AliMUONPairLight(); 
	dimuLight->SetMuons(*mu1, *mu2);
	TClonesArray &dimuons = *dimuonArray;
	new (dimuons[ndimuons++]) AliMUONPairLight(*dimuLight);
      }
    }// end dimuons
    treeOut->Fill(); 
    stack->Reset();
  }//end of events

  RunLoader->UnloadKinematics();
  fout->cd(); 
  treeOut->Write();
}
