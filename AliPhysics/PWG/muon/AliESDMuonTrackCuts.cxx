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

#include "AliESDMuonTrackCuts.h"

#include <AliESDMuonTrack.h>
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliLog.h>

#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectory.h>

//____________________________________________________________________
ClassImp(AliESDMuonTrackCuts)

// Cut names
const Char_t* AliESDMuonTrackCuts::fgkCutNames[kNCuts] = {
 "p",
 "p_{T}",
 "p_{x}",
 "p_{y}",
 "p_{z}",
 "y",
 "eta"
};

//____________________________________________________________________
AliESDMuonTrackCuts::AliESDMuonTrackCuts(const Char_t* name, const Char_t* title) : AliAnalysisCuts(name,title),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapMin(0),
  fRapMax(0),
  fHistogramsOn(0),
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // constructor
  //
  Init();

  //##############################################################################
  // setting default cuts
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapRange();

  SetHistogramsOn();
}

//_____________________________________________________________________________
AliESDMuonTrackCuts::AliESDMuonTrackCuts(const AliESDMuonTrackCuts &c) : AliAnalysisCuts(c),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapMin(0),
  fRapMax(0),
  fHistogramsOn(0),
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // copy constructor
  //
  ((AliESDMuonTrackCuts &) c).Copy(*this);
}

AliESDMuonTrackCuts::~AliESDMuonTrackCuts()
{
  //
  // destructor
  //
  for (Int_t i=0; i<2; i++) {    
    if (fhPt[i])
      delete fhPt[i];
    if (fhEta[i])
      delete fhEta[i];
  }

  if (fhCutStatistics)
    delete fhCutStatistics;             
  if (fhCutCorrelation)
    delete fhCutCorrelation;            
}

void AliESDMuonTrackCuts::Init()
{
  //
  // sets everything to zero
  //
  fPMin = 0;
  fPMax = 0;
  fPtMin = 0;
  fPtMax = 0;
  fPxMin = 0;
  fPxMax = 0;
  fPyMin = 0;
  fPyMax = 0;
  fPzMin = 0;
  fPzMax = 0;
  fEtaMin = 0;
  fEtaMax = 0;
  fRapMin = 0;
  fRapMax = 0;

  fHistogramsOn = kFALSE;

  for (Int_t i=0; i<2; ++i)
  {
    fhPt[i] = 0;
    fhEta[i] = 0;
  }
  fhCutStatistics = 0;
  fhCutCorrelation = 0;
}

//_____________________________________________________________________________
AliESDMuonTrackCuts &AliESDMuonTrackCuts::operator=(const AliESDMuonTrackCuts &c)
{
  //
  // Assignment operator
  //
  if (this != &c) ((AliESDMuonTrackCuts &) c).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliESDMuonTrackCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //

   AliESDMuonTrackCuts& target = (AliESDMuonTrackCuts &) c;
 
  target.Init();
// 
  target.fPMin = fPMin;
  target.fPMax = fPMax;
  target.fPtMin = fPtMin;
  target.fPtMax = fPtMax;
  target.fPxMin = fPxMin;
  target.fPxMax = fPxMax;
  target.fPyMin = fPyMin;
  target.fPyMax = fPyMax;
  target.fPzMin = fPzMin;
  target.fPzMax = fPzMax;
  target.fEtaMin = fEtaMin;
  target.fEtaMax = fEtaMax;
  target.fRapMin = fRapMin;
  target.fRapMax = fRapMax;

  target.fHistogramsOn = fHistogramsOn;

  for (Int_t i=0; i<2; ++i)
  {     
    if (fhPt[i]) target.fhPt[i] = (TH1F*) fhPt[i]->Clone();
    if (fhEta[i]) target.fhEta[i] = (TH1F*) fhEta[i]->Clone();
  }

  if (fhCutStatistics) target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  TNamed::Copy(c);
}

//_____________________________________________________________________________
Long64_t AliESDMuonTrackCuts::Merge(TCollection* list) {
  // Merge a list of AliESDMuonTrackCuts objects with this (needed for PROOF)
  // Returns the number of merged objects (including this)

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  if (!fHistogramsOn)
    return 0;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collection of measured and generated histograms
  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliESDMuonTrackCuts* entry = dynamic_cast<AliESDMuonTrackCuts*>(obj);
    if (entry == 0)
      continue;

    if (!entry->fHistogramsOn)
      continue;
    
    for (Int_t i=0; i<2; i++) {
      fhPt[i]->Add(entry->fhPt[i]); 
      fhEta[i]->Add(entry->fhEta[i]); 
    }      

    fhCutStatistics->Add(entry->fhCutStatistics);        
    fhCutCorrelation ->Add(entry->fhCutCorrelation);      

    count++;
  }

  return count+1;
}

void AliESDMuonTrackCuts::EnableNeededBranches(TTree* tree)
{
  // enables the branches needed by AcceptTrack, for a list see comment of AcceptTrack

  tree->SetBranchStatus("fTracks.fFlags", 1);
  tree->SetBranchStatus("fTracks.fP*", 1);
  tree->SetBranchStatus("fTracks.fR*", 1);  //detector response probability
}

//____________________________________________________________________
Bool_t
AliESDMuonTrackCuts::AcceptTrack(AliESDMuonTrack* esdMuTrack) {
  // 
  // figure out if the tracks survives all the track cuts defined
  //
  // the different kinematic values are first
  // retrieved from the track. then it is found out what cuts the
  // track did not survive and finally the cuts are imposed.

  // getting the kinematic variables of the track
  // (assuming the mass is known)
  Double_t p[3];
  esdMuTrack->PxPyPz(p);
  Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
  Float_t pt       = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2));
  Float_t energy   = TMath::Sqrt(TMath::Power(esdMuTrack->M(),2) + TMath::Power(momentum,2));


  //y-eta related calculations
  Float_t eta = -100.;
  Float_t y   = -100.;
  if((momentum != TMath::Abs(p[2]))&&(momentum != 0))
    eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  if((energy != TMath::Abs(p[2]))&&(momentum != 0))
    y = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
     
  //########################################################################
  // cut the track?
  
  Bool_t cuts[kNCuts];
  for (Int_t i=0; i<kNCuts; i++) cuts[i]=kFALSE;
  
  // track kinematics cut
  if((momentum < fPMin) || (momentum > fPMax)) 
    cuts[0]=kTRUE;
  if((pt < fPtMin) || (pt > fPtMax)) 
    cuts[1] = kTRUE;
  if((p[0] < fPxMin) || (p[0] > fPxMax)) 
    cuts[2] = kTRUE;
  if((p[1] < fPyMin) || (p[1] > fPyMax)) 
    cuts[3] = kTRUE;
  if((p[2] < fPzMin) || (p[2] > fPzMax))
    cuts[4] = kTRUE;
  if((eta < fEtaMin) || (eta > fEtaMax)) 
    cuts[5] = kTRUE;
  if((y < fRapMin) || (y > fRapMax)) 
    cuts[6] = kTRUE;

  Bool_t cut=kFALSE;
  for (Int_t i=0; i<kNCuts; i++) 
    if (cuts[i]) cut = kTRUE;
  
  //########################################################################
  // filling histograms
  if (fHistogramsOn) {
    fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n tracks")));
    
    if (cut)
      fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n cut tracks")));
    
    for (Int_t i=0; i<kNCuts; i++) {
      if (cuts[i])
 	fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin(fgkCutNames[i])));
      
      for (Int_t j=i; j<kNCuts; j++) {
 	if (cuts[i] && cuts[j]) {
 	  Float_t xC = fhCutCorrelation->GetXaxis()->GetBinCenter(fhCutCorrelation->GetXaxis()->FindBin(fgkCutNames[i]));
 	  Float_t yC = fhCutCorrelation->GetYaxis()->GetBinCenter(fhCutCorrelation->GetYaxis()->FindBin(fgkCutNames[j]));
 	  fhCutCorrelation->Fill(xC, yC);
 	}
      }
    }
    
    fhPt[0]->Fill(pt);
    fhEta[0]->Fill(eta);

  }

  //########################################################################
  // cut the track!
  if (cut) return kFALSE;

  //########################################################################
  // filling histograms after cut
  if (fHistogramsOn) {
// 
    fhPt[1]->Fill(pt);
    fhEta[1]->Fill(eta);
    
  }

  return kTRUE;
}

//____________________________________________________________________
TObjArray* AliESDMuonTrackCuts::GetAcceptedTracks(AliESD* esd)
{
  //
  // returns an array of all tracks that pass the cuts
  //

  TObjArray* acceptedTracks = new TObjArray();

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
    AliESDMuonTrack* track = esd->GetMuonTrack(iTrack);

    if (AcceptTrack(track))
      acceptedTracks->Add(track);
  }

  return acceptedTracks;
}

//____________________________________________________________________
Int_t AliESDMuonTrackCuts::CountAcceptedTracks(AliESD* esd)
{
  //
  // returns an the number of tracks that pass the cuts
  //

  Int_t count = 0;

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
    AliESDMuonTrack* track = esd->GetMuonTrack(iTrack);

    if (AcceptTrack(track))
      count++;
  }

  return count;
}

//____________________________________________________________________
TObjArray* AliESDMuonTrackCuts::GetAcceptedTracks(AliESDEvent* esd)
{
  //
  // returns an array of all tracks that pass the cuts
  //

  TObjArray* acceptedTracks = new TObjArray();

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
    AliESDMuonTrack* track = esd->GetMuonTrack(iTrack);

    if (AcceptTrack(track))
      acceptedTracks->Add(track);
  }

  return acceptedTracks;
}

//____________________________________________________________________
Int_t AliESDMuonTrackCuts::CountAcceptedTracks(AliESDEvent* esd)
{
  //
  // returns an the number of tracks that pass the cuts
  //

  Int_t count = 0;

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
    AliESDMuonTrack* track = esd->GetMuonTrack(iTrack);

    if (AcceptTrack(track))
      count++;
  }

  return count;
}

//____________________________________________________________________
 void AliESDMuonTrackCuts::DefineHistograms(Int_t color) {
   // 
   // diagnostics histograms are defined
   // 

   fHistogramsOn=kTRUE;
   
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   
   //###################################################################################
   // defining histograms

   fhCutStatistics = new TH1F("cut_statistics","cut statistics",kNCuts+4,-0.5,kNCuts+3.5);

   fhCutStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
   fhCutStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");

   fhCutCorrelation = new TH2F("cut_correlation","cut correlation",kNCuts,-0.5,kNCuts-0.5,kNCuts,-0.5,kNCuts-0.5);;
  
   for (Int_t i=0; i<kNCuts; i++) {
     fhCutStatistics->GetXaxis()->SetBinLabel(i+4,fgkCutNames[i]);
     fhCutCorrelation->GetXaxis()->SetBinLabel(i+1,fgkCutNames[i]);
     fhCutCorrelation->GetYaxis()->SetBinLabel(i+1,fgkCutNames[i]);
   } 

  fhCutStatistics  ->SetLineColor(color);
  fhCutCorrelation ->SetLineColor(color);
  fhCutStatistics  ->SetLineWidth(2);
  fhCutCorrelation ->SetLineWidth(2);

  Char_t str[256];
  for (Int_t i=0; i<2; i++) {
    if (i==0) snprintf(str,256," ");
    else snprintf(str,256,"_cut");

    fhPt[i]                  = new TH1F(Form("pt%s",str)     ,"p_{T} distribution;p_{T} (GeV/c)",500,0.0,100.0);
    fhEta[i]                 = new TH1F(Form("eta%s",str)     ,"#eta distribution;#eta",40,-2.0,2.0);
  }

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
Bool_t AliESDMuonTrackCuts::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  fhCutStatistics = dynamic_cast<TH1F*> (gDirectory->Get("cut_statistics"));
  fhCutCorrelation = dynamic_cast<TH2F*> (gDirectory->Get("cut_correlation"));

  Char_t str[5];
  for (Int_t i=0; i<2; i++) {
    if (i==0)
    {
      gDirectory->cd("before_cuts");
      str[0] = 0;
    }
    else
    {
      gDirectory->cd("after_cuts");
      snprintf(str,5,"_cut");
    }

    fhPt[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("pt%s",str)));
    fhEta[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("eta%s",str)));

    gDirectory->cd("../");
  }

  gDirectory->cd("..");

  return kTRUE;
}

//____________________________________________________________________
void AliESDMuonTrackCuts::SaveHistograms(const Char_t* dir) {
  //
  // saves the histograms in a directory (dir)
  //

  if (!fHistogramsOn) {
    AliDebug(0, "Histograms not on - cannot save histograms!!!");
    return;
  }

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");

  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t i=0; i<2; i++) {
    if (i==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");

    fhPt[i]                  ->Write();
    fhEta[i]                 ->Write();
    
    gDirectory->cd("../");
  }

  gDirectory->cd("../");
}

//____________________________________________________________________
void AliESDMuonTrackCuts::DrawHistograms()
{
  gStyle->SetPalette(1);
  gStyle->SetFrameFillColor(10);
  gStyle->SetCanvasColor(10);
  
  TCanvas* canvas1 = new TCanvas(Form("%s_1", GetName()), "Track Cut Results", 800, 500);
  canvas1->Divide(2, 1);

  canvas1->cd(1);
  fhCutStatistics->SetStats(kFALSE);
  fhCutStatistics->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  fhCutStatistics->Draw();

  canvas1->cd(2);
  fhCutCorrelation->SetStats(kFALSE);
  fhCutCorrelation->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  gPad->SetLeftMargin(0.3);
  fhCutCorrelation->Draw("COLZ");
  canvas1->Update();
  canvas1->SaveAs(Form("%s_%s.gif", GetName(), canvas1->GetName()));
  
}

