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

//-----------------------------------------------------------------
//         AliAnalysisTaskPhiEffMc class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskPhiEffMc.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODTrackCuts.h"
#include "PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODEventCuts.h"
#include "AliStack.h"
#include <TMCProcess.h>
#include "AliAnalyseLeadingTrackUE.h"

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskPhiEffMc) 

//________________________________________________________________________
AliAnalysisTaskPhiEffMc::AliAnalysisTaskPhiEffMc(const char *name) : AliAnalysisTaskSE(name), 
  fAOD(0x0), 
  fIsMC(0),
  fOutput(0x0),
  fHelperPID(0x0),
  fTrackCuts(0x0),
  fEventCuts(0x0),
  fPtCut(0.)
{
  // Default constructor
  

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliHelperPID::Class());
  DefineOutput(3, AliSpectraAODEventCuts::Class());
  DefineOutput(4, AliSpectraAODTrackCuts::Class());
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskPhiEffMc::UserCreateOutputObjects()
{
  Printf("\n\n\n\n\n\n In CreateOutput Object:");
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  if (!fHelperPID)  AliFatal("HelperPID should be set in the steering macro");

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("BlaBla");

  // event histogram
  const Int_t nEventDim = 4;
  //                      cent  vz MCmult RecoMult
  Int_t nEventBins[4] =   { 10,  41,  50,   50};
  Double_t nEventMin[4] = {  0, -10,   0,    0};
  Double_t nEventMax[4] = {100,  10, 100,  100};

  THnSparseF* hEvent = new THnSparseF("hEvent", "hEvent", nEventDim, nEventBins, nEventMin, nEventMax);
  hEvent->GetAxis(0)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hEvent->GetAxis(1)->SetTitle("v_{z} (cm)");
  hEvent->GetAxis(2)->SetTitle("# MC particles");
  hEvent->GetAxis(3)->SetTitle("# reconstructed tracks");
  fOutput->Add(hEvent);

  TH3F* hKcorr = new TH3F("hKcorr","hKcorr",100,0,10,100,0,10,100,0,10);
  hKcorr->GetXaxis()->SetTitle("#phi p_{T} (GeV/c)");
  hKcorr->GetYaxis()->SetTitle("k+ p_{T} (GeV/c)");
  hKcorr->GetZaxis()->SetTitle("k- p_{T} (GeV/c)");
  fOutput->Add(hKcorr);

  // single charged particles -- bins for THnSparse histograms
  const Int_t nTrackDim = 6;
  //                                PID ID, pt, y, cent, eta, phi
  Int_t    nTrackBins[nTrackDim] =   {   7, 50, 20,  10,  20, 30};
  Double_t nTrackMin[nTrackDim] =    {-3.5,  0, -1,   0,  -1,  0};
  Double_t nTrackMax[nTrackDim] =    { 3.5,  5,  1, 100,   1, 2*TMath::Pi()};

  // single charged particles -- Monte Carlo particles
  THnSparseF* hTrackMc = new THnSparseF("hTrackMc", "hTrackMc", nTrackDim, nTrackBins, nTrackMin, nTrackMax);
  hTrackMc->GetAxis(0)->SetTitle("PID ID * Charge");
  hTrackMc->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hTrackMc->GetAxis(2)->SetTitle("y");
  hTrackMc->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hTrackMc->GetAxis(4)->SetTitle("#eta");
  hTrackMc->GetAxis(5)->SetTitle("#phi (rad)");
  fOutput->Add(hTrackMc);

  // single charged particles -- reconstructed tracks
  THnSparseF* hTrackReco = new THnSparseF("hTrackReco", "hTrackReco", nTrackDim, nTrackBins, nTrackMin, nTrackMax);
  hTrackReco->GetAxis(0)->SetTitle("PID ID * Charge");
  hTrackReco->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hTrackReco->GetAxis(2)->SetTitle("y");
  hTrackReco->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hTrackReco->GetAxis(4)->SetTitle("#eta");
  hTrackReco->GetAxis(5)->SetTitle("#phi (rad)");
  fOutput->Add(hTrackReco);

  // single charged particles -- Monte Carlo particles matched to reconstructed tracks
  THnSparseF* hTrackMatch = new THnSparseF("hTrackMatch", "hTrackMatch", nTrackDim, nTrackBins, nTrackMin, nTrackMax);
  hTrackMatch->GetAxis(0)->SetTitle("PID ID * Charge");
  hTrackMatch->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hTrackMatch->GetAxis(2)->SetTitle("y");
  hTrackMatch->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hTrackMatch->GetAxis(4)->SetTitle("#eta");
  hTrackMatch->GetAxis(5)->SetTitle("#phi (rad)");
  fOutput->Add(hTrackMatch);

  // kaon pairs -- bins for THnSparse histograms
  const Int_t nPairDim = 7;
  //                             InvMass, pt, y, cent, eta, phi,          pair ID
  Int_t    nPairBins[nPairDim] =   { 100, 50, 20,  10,  20, 30,              3};
  Double_t nPairMin[nPairDim] =    {0.98,  0, -1,   0,  -1,  0,           -0.5};
  Double_t nPairMax[nPairDim] =    { 1.1,  5,  1, 100,   1, 2*TMath::Pi(), 2.5};
  // pair ID = 0 unlike sign k+k-
  //           1 like sign pos k+k+
  //           2 like sign neg k-k-

  // kaon pairs -- Monte Carlo particles -- real phi's
  THnSparseF* hPairMcPhi = new THnSparseF("hPairMcPhi", "hPairMcPhi", nPairDim, nPairBins, nPairMin, nPairMax);
  hPairMcPhi->GetAxis(0)->SetTitle("Invariant Mass (GeV/c/c)");
  hPairMcPhi->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hPairMcPhi->GetAxis(2)->SetTitle("y");
  hPairMcPhi->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hPairMcPhi->GetAxis(4)->SetTitle("#eta");
  hPairMcPhi->GetAxis(5)->SetTitle("#phi (rad)");
  hPairMcPhi->GetAxis(6)->SetTitle("pair ID");
  fOutput->Add(hPairMcPhi);

  // kaon pairs -- Monte Carlo particles -- real phi's with kaon daughters that pass pt and eta cuts
  THnSparseF* hPairMcPhiCuts = new THnSparseF("hPairMcPhiCuts", "hPairMcPhiCuts", nPairDim, nPairBins, nPairMin, nPairMax);
  hPairMcPhiCuts->GetAxis(0)->SetTitle("Invariant Mass (GeV/c/c)");
  hPairMcPhiCuts->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hPairMcPhiCuts->GetAxis(2)->SetTitle("y");
  hPairMcPhiCuts->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hPairMcPhiCuts->GetAxis(4)->SetTitle("#eta");
  hPairMcPhiCuts->GetAxis(5)->SetTitle("#phi (rad)");
  hPairMcPhiCuts->GetAxis(6)->SetTitle("pair ID");
  fOutput->Add(hPairMcPhiCuts);

  // kaon pairs -- Monte Carlo particles -- phi's reconstructed from MC kaons
  THnSparseF* hPairMc = new THnSparseF("hPairMc", "hPairMc", nPairDim, nPairBins, nPairMin, nPairMax);
  hPairMc->GetAxis(0)->SetTitle("Invariant Mass (GeV/c/c)");
  hPairMc->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hPairMc->GetAxis(2)->SetTitle("y");
  hPairMc->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hPairMc->GetAxis(4)->SetTitle("#eta");
  hPairMc->GetAxis(5)->SetTitle("#phi (rad)");
  hPairMc->GetAxis(6)->SetTitle("pair ID");
  fOutput->Add(hPairMc);

  // kaon pairs -- reconstructed tracks
  THnSparseF* hPairReco = new THnSparseF("hPairReco", "hPairReco", nPairDim, nPairBins, nPairMin, nPairMax);
  hPairReco->GetAxis(0)->SetTitle("Invariant Mass (GeV/c/c)");
  hPairReco->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hPairReco->GetAxis(2)->SetTitle("y");
  hPairReco->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hPairReco->GetAxis(4)->SetTitle("#eta");
  hPairReco->GetAxis(5)->SetTitle("#phi (rad)");
  hPairReco->GetAxis(6)->SetTitle("pair ID");
  fOutput->Add(hPairReco);

  // kaon pairs -- Monte Carlo particles matched to reconstructed tracks
  THnSparseF* hPairMatch = new THnSparseF("hPairMatch", "hPairMatch", nPairDim, nPairBins, nPairMin, nPairMax);
  hPairMatch->GetAxis(0)->SetTitle("Invariant Mass (GeV/c/c)");
  hPairMatch->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  hPairMatch->GetAxis(2)->SetTitle("y");
  hPairMatch->GetAxis(3)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  hPairMatch->GetAxis(4)->SetTitle("#eta");
  hPairMatch->GetAxis(5)->SetTitle("#phi (rad)");
  hPairMatch->GetAxis(6)->SetTitle("pair ID");
  fOutput->Add(hPairMatch);

  PostData(1, fOutput);
  PostData(2, fHelperPID);
  PostData(3, fEventCuts);
  PostData(4, fTrackCuts);
}

//________________________________________________________________________
void AliAnalysisTaskPhiEffMc::UserExec(Option_t *)
{
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }

  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection

  Double_t cent = fEventCuts->GetCent();

  THnSparseF* hEvent = (THnSparseF*)fOutput->FindObject("hEvent");
  THnSparseF* hTrackMc = (THnSparseF*)fOutput->FindObject("hTrackMc");
  THnSparseF* hTrackReco = (THnSparseF*)fOutput->FindObject("hTrackReco");
  THnSparseF* hTrackMatch = (THnSparseF*)fOutput->FindObject("hTrackMatch");
  THnSparseF* hPairMcPhi = (THnSparseF*)fOutput->FindObject("hPairMcPhi");
  THnSparseF* hPairMcPhiCuts = (THnSparseF*)fOutput->FindObject("hPairMcPhiCuts");
  THnSparseF* hPairMc = (THnSparseF*)fOutput->FindObject("hPairMc");
  THnSparseF* hPairReco = (THnSparseF*)fOutput->FindObject("hPairReco");
  THnSparseF* hPairMatch = (THnSparseF*)fOutput->FindObject("hPairMatch");
  TH3F* hKcorr = (TH3F*)fOutput->FindObject("hKcorr");

  TObjArray* kaonsPosMc = new TObjArray();
  TObjArray* kaonsNegMc = new TObjArray();
  TObjArray* kaonsPosData = new TObjArray();
  TObjArray* kaonsNegData = new TObjArray();
  TObjArray* kaonsPosGen = new TObjArray();
  TObjArray* kaonsNegGen = new TObjArray();

  Int_t nMc = 0; // number of charged MC particles which satisfy cuts (except PID)
  Int_t nReco = 0; // number of reco tracks which satisfy cuts (except PID)

  //MC Loop
  TClonesArray *arrayMC = 0;
  if (fIsMC)
    {
      arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!arrayMC) {
	AliFatal("Error: MC particles branch not found!\n");
      }


      Int_t nMC = arrayMC->GetEntries();
      for (Int_t iMC = 0; iMC < nMC; iMC++)
	{
	  AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);

	  if(partMC->Pt()<fPtCut) continue;
	  if(partMC->Eta()>fTrackCuts->GetEtaMax() || partMC->Eta()<fTrackCuts->GetEtaMin() ) continue;

	  // PDG codes: pion(211), kaon(+/-321), proton(2212), phi(333)

	  // make pt spectrum of only MC phis that decay to k+k-
	  if(partMC->GetPdgCode()==333)
	    {
	      if(partMC->GetNDaughters()==2)
		{
		  if(TMath::Abs(((AliAODMCParticle*) arrayMC->At(TMath::Abs(partMC->GetDaughter(0))))->GetPdgCode()) == 321 && TMath::Abs(((AliAODMCParticle*) arrayMC->At(TMath::Abs(partMC->GetDaughter(1))))->GetPdgCode()) == 321)
		    {
		      AliAODMCParticle *d1 = (AliAODMCParticle*) arrayMC->At(TMath::Abs(partMC->GetDaughter(0)));
		      AliAODMCParticle *d2 = (AliAODMCParticle*) arrayMC->At(TMath::Abs(partMC->GetDaughter(1)));
		      if(d1->Charge() > 0) hKcorr->Fill(partMC->Pt(),d1->Pt(),d2->Pt());
		      else hKcorr->Fill(partMC->Pt(),d2->Pt(),d1->Pt());
		      TLorentzVector *d3 = (TLorentzVector*)makePhi(d1,d2);
		      // InvMass, pt, y, eta, phi, pair ID
		      Double_t invMass = (d3->M2() > 0 ? sqrt(d3->M2()) : 0);
		      if(invMass < 1.11)
			{
			  Double_t varfill1[7] = {invMass, d3->Pt(), d3->Rapidity(), cent, d3->Eta(), (d3->Phi() > 0 ? d3->Phi() : d3->Phi()+2*TMath::Pi()), 0};
			  hPairMcPhi->Fill(varfill1);
			  if(d1->Pt()>fPtCut && d2->Pt()>fPtCut && d1->Eta()<fTrackCuts->GetEtaMax() && d1->Eta()>fTrackCuts->GetEtaMin() && d2->Eta()<fTrackCuts->GetEtaMax() && d2->Eta()>fTrackCuts->GetEtaMin())
			    hPairMcPhiCuts->Fill(varfill1);
			}
		      delete d3;
		    }
		}
	    }

	  if(partMC->Charge()==0) continue; // remove neutrals
	  if(!partMC->IsPhysicalPrimary()) continue;  // need to remove all non-final state particles

	  nMc++;

	  Int_t mcID = fHelperPID->GetParticleSpecies(partMC);
	  if(mcID>3 || mcID < -3) continue;
	  // PID ID, pt, y, eta, phi
	  Double_t varfill2[6] = {(partMC->Charge() > 0 ? mcID+1 : -1*(mcID+1)), partMC->Pt(), partMC->Y(), cent, partMC->Eta(), partMC->Phi()};
	  hTrackMc->Fill(varfill2);

	  if(mcID == 1)
	    {
	      if(partMC->Charge() > 0) kaonsPosMc->Add(partMC);
	      else if(partMC->Charge() < 0) kaonsNegMc->Add(partMC);
	    }
	  
	} // end MC track loop
      
      UnlikeSign(kaonsPosMc,kaonsNegMc,hPairMc,cent);
      LikeSign(kaonsPosMc,kaonsNegMc,hPairMc,cent);
    }


  //track loop
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track,kTRUE)) continue; //track selection (rapidity selection NOT in the standard cuts?)
    if(track->Charge()==0) continue;
    if(track->Pt()<fPtCut) continue;
    if(track->Eta()>fTrackCuts->GetEtaMax() || track->Eta()<fTrackCuts->GetEtaMin()) continue;

    nReco++;

    Int_t dataID=fHelperPID->GetParticleSpecies(track,kTRUE);
    if(dataID>3 || dataID < -3) continue;

    // PID ID, pt, y, eta, phi
    Double_t varfill3[6] = {(track->Charge() > 0 ? dataID+1 : -1*(dataID+1)), track->Pt(), track->Y(), cent, track->Eta(), track->Phi()};
    hTrackReco->Fill(varfill3);

    if(dataID==1)
      {
	if(track->Charge() > 0) kaonsPosData->Add(track);
	else if(track->Charge() < 0) kaonsNegData->Add(track);
      }

    if (fIsMC)
      {
	Int_t mcLabel = track->GetLabel();
	AliAODMCParticle* tempMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(mcLabel));
	if(!tempMC)
	  {
	    AliError("Cannot get MC particle");
	    continue; 
	  }
	
	Int_t genID = fHelperPID->GetParticleSpecies(tempMC);
	if(genID != dataID) continue;
	if(genID>3 || genID < -3) continue;

	// PID ID, pt, y, eta, phi
	Double_t varfill4[6] = {(track->Charge() > 0 ? genID+1 : -1*(genID+1)), track->Pt(), track->Y(), cent, track->Eta(), track->Phi()};
	hTrackMatch->Fill(varfill4);

	if(genID==1)
	  {
	    Int_t motherId = TMath::Abs(tempMC->GetMother());
	    
	    AliAODMCParticle* genMother = (AliAODMCParticle*)arrayMC->At(motherId);
	    if(!genMother) continue;
	    
	    if(genMother->GetPdgCode() != 333) continue;
	    
	    if(track->Charge() > 0)
	      {
		for(Int_t k = 0; k < kaonsNegGen->GetEntriesFast(); k++)
		  {
		    AliVParticle* tempGenMatch = (AliVParticle*)kaonsNegGen->UncheckedAt(k);
		    AliAODMCParticle* tempMatchMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(tempGenMatch->GetLabel()));
		    if(!tempMatchMC) continue;
		    if(tempMatchMC->GetMother() != motherId) continue;
		    
		    TLorentzVector*c = (TLorentzVector*)makePhi(track,tempGenMatch);
		    Double_t invMass = (c->M2() > 0 ? sqrt(c->M2()) : 0);
		    if(invMass < 1.11)
		      {
			Double_t fillgen[7] = {invMass, c->Pt(), c->Rapidity(), cent, c->Eta(), (c->Phi() > 0 ? c->Phi() : c->Phi()+2*TMath::Pi()), 0};
			hPairMatch->Fill(fillgen);
		      }
		    delete c;
		  } 
	      }
	    else if(track->Charge() < 0)
	      {
		for(Int_t k = 0; k < kaonsPosGen->GetEntriesFast(); k++)
		  {
		    AliVParticle* tempGenMatch = (AliVParticle*)kaonsPosGen->UncheckedAt(k);
		    AliAODMCParticle* tempMatchMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(tempGenMatch->GetLabel()));
		    if(!tempMatchMC) continue;
		    if(tempMatchMC->GetMother() != motherId) continue;
		    
		    TLorentzVector*c = (TLorentzVector*)makePhi(track,tempGenMatch);
		    Double_t invMass = (c->M2() > 0 ? sqrt(c->M2()) : 0);
		    if(invMass < 1.11)
		      {
			Double_t fillgen2[7] = {invMass, c->Pt(), c->Rapidity(), cent, c->Eta(), (c->Phi() > 0 ? c->Phi() : c->Phi()+2*TMath::Pi()), 0};
			hPairMatch->Fill(fillgen2);
		      }
		    delete c;
		  }
	      }
	    
	    if(track->Charge() > 0) kaonsPosGen->Add(track);
	    else if(track->Charge() < 0) kaonsNegGen->Add(track);
	  }
      }
	
  } // end loop on tracks
  
  UnlikeSign(kaonsPosData,kaonsNegData,hPairReco,cent);
  LikeSign(kaonsPosData,kaonsNegData,hPairReco,cent);
  
  Double_t varfill5[4] = {cent, fAOD->GetPrimaryVertex()->GetZ(), nMc, nReco};
  hEvent->Fill(varfill5);

  kaonsPosMc->Clear();
  kaonsNegMc->Clear();
  kaonsPosData->Clear();
  kaonsNegData->Clear();
  kaonsPosGen->Clear();
  kaonsNegGen->Clear();

  if(kaonsPosMc) delete kaonsPosMc;
  if(kaonsNegMc) delete kaonsNegMc;
  if(kaonsPosData) delete kaonsPosData;
  if(kaonsNegData) delete kaonsNegData;
  if(kaonsPosGen) delete kaonsPosGen;
  if(kaonsNegGen) delete kaonsNegGen;

  PostData(1,fOutput);
  PostData(2, fHelperPID);
  PostData(3, fEventCuts);
  PostData(4, fTrackCuts);
  //Printf("............. end of Exec");
  
}

//_________________________________________________________________
void   AliAnalysisTaskPhiEffMc::Terminate(Option_t *)
{
  // Terminate analysis
  //
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  } 
  
    printf("AliAnalysisTaskPhiEffMc: Terminate() \n");
  
}

//_________________________________________________________________
void AliAnalysisTaskPhiEffMc::UnlikeSign(TObjArray* kaonsPos, TObjArray* kaonsNeg, THnSparseF* h, Double_t cent)
{
  Int_t countPos = kaonsPos->GetEntriesFast();
  Int_t countNeg = kaonsNeg->GetEntriesFast();

  if(countPos+countNeg < 2) return;

  for(Int_t cp = 0; cp < countPos; cp++) // form k+k- pairs and compute invariant mass
    {
      AliVParticle* temp1 = (AliVParticle*)kaonsPos->UncheckedAt(cp);
      for(Int_t cn = 0; cn < countNeg; cn++)
	{
	  AliVParticle* temp2 = (AliVParticle*)kaonsNeg->UncheckedAt(cn);
	  TLorentzVector*c = (TLorentzVector*)makePhi(temp1,temp2);
	  // InvMass, pt, y, eta, phi, pair ID
	  Double_t invMass = (c->M2() > 0 ? sqrt(c->M2()) : 0);
	  if(invMass < 1.11)
	    {
	      Double_t varfill[7] = {invMass, c->Pt(), c->Rapidity(), cent, c->Eta(), (c->Phi() > 0 ? c->Phi() : c->Phi()+2*TMath::Pi()), 0};
	      h->Fill(varfill);
	    }
	  delete c;
	}
    }
}

void AliAnalysisTaskPhiEffMc::LikeSign(TObjArray* kaonsPos, TObjArray* kaonsNeg, THnSparseF* h, Double_t cent)
{
  Int_t countPos = kaonsPos->GetEntriesFast();
  Int_t countNeg = kaonsNeg->GetEntriesFast();

  if(countPos < 2 && countNeg < 2) return;

  for(Int_t cp = 0; cp < countPos; cp++) // form k+k+ pairs and compute invariant mass
    {
      AliVParticle* temp1 = (AliVParticle*)kaonsPos->UncheckedAt(cp);
      for(Int_t cn = cp+1; cn < countPos; cn++)
	{
	  AliVParticle* temp2 = (AliVParticle*)kaonsPos->UncheckedAt(cn);
	  TLorentzVector*c = (TLorentzVector*)makePhi(temp1,temp2);
	  // InvMass, pt, y, eta, phi, pair ID
	  Double_t invMass = (c->M2() > 0 ? sqrt(c->M2()) : 0);
	  if(invMass < 1.11)
	    {
	      //cout << "y = " << 0.5*log((c->E()+c->Pz())/(c->E()-c->Pz())) << "   " << c->Rapidity() << endl;
	      //cout << "eta = " << -1.*log(tan(c->Theta()/2.)) << "   " << c->Eta() << endl;
	      Double_t varfill1[7] = {invMass, c->Pt(), c->Rapidity(), cent, c->Eta(), (c->Phi() > 0 ? c->Phi() : c->Phi()+2*TMath::Pi()), 1};
	      h->Fill(varfill1);
	    }
	  delete c;
	}
    }
  for(Int_t cp = 0; cp < countNeg; cp++) // form k-k- pairs and compute invariant mass
    {
      AliVParticle* temp1 = (AliVParticle*)kaonsNeg->UncheckedAt(cp);
      for(Int_t cn = cp+1; cn < countNeg; cn++)
	{
	  AliVParticle* temp2 = (AliVParticle*)kaonsNeg->UncheckedAt(cn);
	  TLorentzVector*c = (TLorentzVector*)makePhi(temp1,temp2);
	  Double_t invMass = (c->M2() > 0 ? sqrt(c->M2()) : 0);
	  if(invMass < 1.11)
	    {
	      Double_t varfill2[7] = {invMass, c->Pt(), c->Rapidity(), cent, c->Eta(), (c->Phi() > 0 ? c->Phi() : c->Phi()+2*TMath::Pi()), 2};
	      h->Fill(varfill2);
	    }
	  delete c;
	} 
    }
}

TLorentzVector* AliAnalysisTaskPhiEffMc::makePhi(AliVParticle* p1, AliVParticle* p2)
{
  TLorentzVector* a = new TLorentzVector(p1->Px(),p1->Py(),p1->Pz(),sqrt(pow(p1->P(),2)+pow(4.93676999999999977e-01, 2)));
  TLorentzVector* b = new TLorentzVector(p2->Px(),p2->Py(),p2->Pz(),sqrt(pow(p2->P(),2)+pow(4.93676999999999977e-01, 2)));
  TLorentzVector* c = new TLorentzVector((*a)+(*b));
  delete a;
  delete b;
  return c;
}


