/* $Id$ */

#include "AlidNdEtaCorrection.h"

#include <TCanvas.h>
#include <TH3F.h>
#include <TH1D.h>

//____________________________________________________________________
ClassImp(AlidNdEtaCorrection)

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection(Char_t* name) 
  : TNamed(name, name),
  fNEvents(0),
  fNTriggeredEvents(0)
{  
  // constructor
  //

  fTrack2ParticleCorrection = new AliCorrectionMatrix3D("nTrackToNPart", "nTrackToNPart",80,-20,20,120,-6,6, 100, 0, 10);

  Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 
			    10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
  Float_t binLimitsVtx[] = {-20,-15,-10,-6,-3,0,3,6,10,15,20};
  
  fVertexRecoCorrection        = new AliCorrectionMatrix2D("vtxReco",       "vtxReco",10,binLimitsVtx ,22,binLimitsN);

  fTriggerBiasCorrection       = new AliCorrectionMatrix2D("triggerBias",   "triggerBias",120,-6,6,100, 0, 10);

  fTrack2ParticleCorrection ->SetAxisTitles("vtx z [cm]", "#eta", "p_{T}");
  fVertexRecoCorrection        ->SetAxisTitles("vtx z [cm]", "n particles/tracks/tracklets?");

  fTriggerBiasCorrection       ->SetAxisTitles("#eta", "p_{T} [GeV/c]");
}

//____________________________________________________________________
void
AlidNdEtaCorrection::Finish() {
  //
  // finish method
  //
  // divide the histograms in the AliCorrectionMatrix2D objects to get the corrections


  fTrack2ParticleCorrection->Divide();

  fVertexRecoCorrection->Divide();

  fTriggerBiasCorrection->GetMeasuredHistogram()->Scale(Double_t(fNTriggeredEvents)/Double_t(fNEvents));
  fTriggerBiasCorrection->Divide();

}

//____________________________________________________________________
Long64_t 
AlidNdEtaCorrection::Merge(TCollection* list) {
  // Merge a list of dNdEtaCorrection objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionNtrackToNparticle = new TList;
  TList* collectionVertexReco        = new TList;
  TList* collectionTriggerBias       = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AlidNdEtaCorrection* entry = dynamic_cast<AlidNdEtaCorrection*> (obj);
    if (entry == 0) 
      continue;

    collectionNtrackToNparticle ->Add(entry->GetTrack2ParticleCorrection());
    collectionVertexReco        ->Add(entry->GetVertexRecoCorrection());
    collectionTriggerBias        ->Add(entry->GetTriggerBiasCorrection());

    count++;
  }
  fTrack2ParticleCorrection ->Merge(collectionNtrackToNparticle);
  fVertexRecoCorrection        ->Merge(collectionVertexReco);
  fTriggerBiasCorrection        ->Merge(collectionTriggerBias);
  
  delete collectionNtrackToNparticle;
  delete collectionVertexReco;
  delete collectionTriggerBias;

  return count+1;
}


//____________________________________________________________________
Bool_t
AlidNdEtaCorrection::LoadHistograms(Char_t* fileName, Char_t* dir) {
  //
  // loads the histograms
  //

  fTrack2ParticleCorrection ->LoadHistograms(fileName, dir);
  fVertexRecoCorrection        ->LoadHistograms(fileName, dir);
  fTriggerBiasCorrection       ->LoadHistograms(fileName, dir);
  
  return kTRUE;
}


//____________________________________________________________________
void
AlidNdEtaCorrection::SaveHistograms() {
  //
  // save the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fTrack2ParticleCorrection ->SaveHistograms();
  fVertexRecoCorrection        ->SaveHistograms();
  fTriggerBiasCorrection       ->SaveHistograms();

  gDirectory->cd("../");
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the two AliCorrectionMatrix2D objects

  fTrack2ParticleCorrection ->DrawHistograms();
  fVertexRecoCorrection        ->DrawHistograms();
  fTriggerBiasCorrection       ->DrawHistograms();

}

//____________________________________________________________________
Float_t AlidNdEtaCorrection::GetMeasuredFraction(Float_t ptCutOff, Float_t eta, Bool_t debug)
{
  // calculates the fraction of particles measured (some are missed due to the pt cut off)
  // uses the generated particle histogram from fTrack2ParticleCorrection

  TH3F* generated = fTrack2ParticleCorrection->GetGeneratedHistogram();

  // find eta borders, if eta is negative assume -0.8 ... 0.8
  Int_t etaBegin = 0;
  Int_t etaEnd = 0;
  if (eta < 0)
  {
    etaBegin = generated->GetYaxis()->FindBin(-0.8);
    etaEnd = generated->GetYaxis()->FindBin(0.8);
  }
  else
  {
    etaBegin = generated->GetYaxis()->FindBin(eta);
    etaEnd = etaBegin;
  }

  Int_t vertexBegin = generated->GetXaxis()->FindBin(-10);
  Int_t vertexEnd = generated->GetXaxis()->FindBin(10);

  TH1D* ptProj = dynamic_cast<TH1D*> (generated->ProjectionZ(Form("%s_pt", GetName()), vertexBegin, vertexEnd, etaBegin, etaEnd));

  Int_t ptBin = ptProj->FindBin(ptCutOff);
  Float_t abovePtCut = ptProj->Integral(ptBin, ptProj->GetNbinsX());
  Float_t all = ptProj->Integral();

  if (all == 0)
    return -1;

  Float_t fraction = abovePtCut / all;

  if (debug)
  {
    new TCanvas;
    ptProj->Draw();
  }

  return fraction;
}

