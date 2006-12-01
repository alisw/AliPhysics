/* $Id$ */

#include "AlidNdEtaCorrection.h"

#include <AliLog.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TH1D.h>
#include <AliCorrection.h>
#include <AliCorrectionMatrix2D.h>
#include <AliCorrectionMatrix3D.h>

//____________________________________________________________________
ClassImp(AlidNdEtaCorrection)

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection()
  : TNamed(),
  fTrack2ParticleCorrection(0),
  fVertexRecoCorrection(0),
  fTriggerBiasCorrectionMBToINEL(0),
  fTriggerBiasCorrectionMBToNSD(0),
  fTriggerBiasCorrectionMBToND(0)
{
  // default constructor
}

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection(const Char_t* name, const Char_t* title)
  : TNamed(name, title),
  fTrack2ParticleCorrection(0),
  fVertexRecoCorrection(0),
  fTriggerBiasCorrectionMBToINEL(0),
  fTriggerBiasCorrectionMBToNSD(0),
  fTriggerBiasCorrectionMBToND(0)
{
  //
  // constructor
  //

  fTrack2ParticleCorrection = new AliCorrection("Track2Particle", "Track2Particle");
  fVertexRecoCorrection     = new AliCorrection("VertexReconstruction", "VertexReconstruction");

  fTriggerBiasCorrectionMBToINEL = new AliCorrection("TriggerBias_MBToINEL", "TriggerBias_MBToINEL");
  fTriggerBiasCorrectionMBToNSD  = new AliCorrection("TriggerBias_MBToNSD", "TriggerBias_MBToNSD");
  fTriggerBiasCorrectionMBToND   = new AliCorrection("TriggerBias_MBToND", "TriggerBias_MBToND");
}

//____________________________________________________________________
AlidNdEtaCorrection::~AlidNdEtaCorrection()
{
  // destructor

  if (fTrack2ParticleCorrection) {
    delete fTrack2ParticleCorrection;
    fTrack2ParticleCorrection = 0;
  }

  if (fVertexRecoCorrection) {
    delete fVertexRecoCorrection;
    fVertexRecoCorrection = 0;
  }

  if (fTriggerBiasCorrectionMBToINEL) {
    delete fTriggerBiasCorrectionMBToINEL;
    fTriggerBiasCorrectionMBToINEL = 0;
  }

  if (fTriggerBiasCorrectionMBToNSD) {
    delete fTriggerBiasCorrectionMBToNSD;
    fTriggerBiasCorrectionMBToNSD = 0;
  }

  if (fTriggerBiasCorrectionMBToND) {
    delete fTriggerBiasCorrectionMBToND;
    fTriggerBiasCorrectionMBToND = 0;
  }
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
  fTriggerBiasCorrectionMBToINEL->Divide();
  fTriggerBiasCorrectionMBToNSD->Divide();
  fTriggerBiasCorrectionMBToND->Divide();
}

//____________________________________________________________________
Long64_t AlidNdEtaCorrection::Merge(TCollection* list)
{
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
  TList* collectionNtrackToNparticle    = new TList;
  TList* collectionVertexReco           = new TList;
  TList* collectionTriggerBiasMBToINEL  = new TList;
  TList* collectionTriggerBiasMBToNSD   = new TList;
  TList* collectionTriggerBiasMBToND    = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AlidNdEtaCorrection* entry = dynamic_cast<AlidNdEtaCorrection*> (obj);
    if (entry == 0)
      continue;

    collectionNtrackToNparticle  ->Add(entry->fTrack2ParticleCorrection);
    collectionVertexReco         ->Add(entry->fVertexRecoCorrection);
    collectionTriggerBiasMBToINEL->Add(entry->fTriggerBiasCorrectionMBToINEL);
    collectionTriggerBiasMBToNSD ->Add(entry->fTriggerBiasCorrectionMBToNSD);
    collectionTriggerBiasMBToND  ->Add(entry->fTriggerBiasCorrectionMBToND);

    count++;
  }
  fTrack2ParticleCorrection      ->Merge(collectionNtrackToNparticle);
  fVertexRecoCorrection          ->Merge(collectionVertexReco);
  fTriggerBiasCorrectionMBToINEL ->Merge(collectionTriggerBiasMBToINEL);
  fTriggerBiasCorrectionMBToNSD  ->Merge(collectionTriggerBiasMBToNSD);
  fTriggerBiasCorrectionMBToND   ->Merge(collectionTriggerBiasMBToND);

  delete collectionNtrackToNparticle;
  delete collectionVertexReco;
  delete collectionTriggerBiasMBToINEL;
  delete collectionTriggerBiasMBToNSD;
  delete collectionTriggerBiasMBToND;

  return count+1;
}

//____________________________________________________________________
Bool_t AlidNdEtaCorrection::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  fTrack2ParticleCorrection      ->LoadHistograms();
  fVertexRecoCorrection          ->LoadHistograms();
  fTriggerBiasCorrectionMBToINEL ->LoadHistograms();
  fTriggerBiasCorrectionMBToNSD  ->LoadHistograms();
  fTriggerBiasCorrectionMBToND   ->LoadHistograms();

  gDirectory->cd("..");

  return kTRUE;
}

//____________________________________________________________________
void AlidNdEtaCorrection::SaveHistograms()
{
  //
  // save the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fTrack2ParticleCorrection     ->SaveHistograms();
  fVertexRecoCorrection         ->SaveHistograms();
  fTriggerBiasCorrectionMBToINEL->SaveHistograms();
  fTriggerBiasCorrectionMBToNSD ->SaveHistograms();
  fTriggerBiasCorrectionMBToND  ->SaveHistograms();

  gDirectory->cd("..");
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the correction
  //

  fTrack2ParticleCorrection     ->DrawHistograms();
  fVertexRecoCorrection         ->DrawHistograms();
  fTriggerBiasCorrectionMBToINEL->DrawHistograms();
  fTriggerBiasCorrectionMBToNSD ->DrawHistograms();
  fTriggerBiasCorrectionMBToND  ->DrawHistograms();
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillMCParticle(Float_t vtx, Float_t eta, Float_t pt, Bool_t trigger, Bool_t vertex, Int_t processType)
{
  // fills a particle in the corrections
  // it is filled in generated or measured depending of the flags

  fTriggerBiasCorrectionMBToINEL->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (processType != 92 && processType != 93)
    fTriggerBiasCorrectionMBToNSD->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (processType!=92 && processType!=93 && processType!=94)
    fTriggerBiasCorrectionMBToND->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (!trigger)
    return;

  fTriggerBiasCorrectionMBToINEL->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTriggerBiasCorrectionMBToNSD->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTriggerBiasCorrectionMBToND->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fVertexRecoCorrection->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (!vertex)
    return;

  fVertexRecoCorrection->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTrack2ParticleCorrection->GetTrackCorrection()->FillGene(vtx, eta, pt);
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillTrackedParticle(Float_t vtx, Float_t eta, Float_t pt)
{
  // fills a tracked particle in the corrections

  fTrack2ParticleCorrection->GetTrackCorrection()->FillMeas(vtx, eta, pt);
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillEvent(Float_t vtx, Float_t n, Bool_t trigger, Bool_t vertex, Int_t processType)
{
  // fills an event int he correction
  // it is filled in generated or measured depending of the flags

  fTriggerBiasCorrectionMBToINEL->GetEventCorrection()->FillGene(vtx, n);

  if (processType != 92 && processType != 93)
    fTriggerBiasCorrectionMBToNSD->GetEventCorrection()->FillGene(vtx, n);

  if (processType!=92 && processType!=93 && processType!=94)
    fTriggerBiasCorrectionMBToND->GetEventCorrection()->FillGene(vtx, n);

  if (!trigger)
    return;

  fTriggerBiasCorrectionMBToINEL->GetEventCorrection()->FillMeas(vtx, n);
  fTriggerBiasCorrectionMBToNSD->GetEventCorrection()->FillMeas(vtx, n);
  fTriggerBiasCorrectionMBToND->GetEventCorrection()->FillMeas(vtx, n);
  fVertexRecoCorrection->GetEventCorrection()->FillGene(vtx, n);

  if (!vertex)
    return;

  fVertexRecoCorrection->GetEventCorrection()->FillMeas(vtx, n);
}

//____________________________________________________________________
Float_t AlidNdEtaCorrection::GetMeasuredFraction(CorrectionType correctionType, Float_t ptCutOff, Float_t eta, Bool_t debug)
{
  // calculates the fraction of particles measured (some are missed due to the pt cut off)
  //
  // uses the generated particle histogram from the correction passed, e.g. pass GetTrack2ParticleCorrection()

  const TH3F* generated = 0;

  switch (correctionType)
  {
    case kNone : return -1;
    case kTrack2Particle : generated = fTrack2ParticleCorrection->GetTrackCorrection()->GetGeneratedHistogram(); break;
    case kVertexReco : generated = fVertexRecoCorrection->GetTrackCorrection()->GetGeneratedHistogram(); break;
    case kINEL : generated = fTriggerBiasCorrectionMBToINEL->GetTrackCorrection()->GetGeneratedHistogram(); break;
    case kNSD: generated = fTriggerBiasCorrectionMBToNSD->GetTrackCorrection()->GetGeneratedHistogram(); break;
    case kND: generated = fTriggerBiasCorrectionMBToND->GetTrackCorrection()->GetGeneratedHistogram(); break;
  }
  
  // find eta borders, if eta is negative assume -0.8 ... 0.8
  Int_t etaBegin = 0;
  Int_t etaEnd = 0;
  if (eta < -99)
  {
    etaBegin = generated->GetYaxis()->FindBin(-0.8);
    etaEnd = generated->GetYaxis()->FindBin(0.8);
  }
  else
  {
    etaBegin = generated->GetYaxis()->FindBin(eta);
    etaEnd = etaBegin;
  }

  Int_t vertexBegin = generated->GetXaxis()->FindBin(-4.99);
  Int_t vertexEnd = generated->GetXaxis()->FindBin(4.99);

  TH1D* ptProj = dynamic_cast<TH1D*> (generated->ProjectionZ(Form("%s_pt", generated->GetName()), vertexBegin, vertexEnd, etaBegin, etaEnd));
  //printf("GetMeasuredFraction: bin range %d %d %d %d\n", vertexBegin, vertexEnd, etaBegin, etaEnd);
  ptProj->GetXaxis()->SetTitle(generated->GetZaxis()->GetTitle());

  Int_t ptBin = ptProj->FindBin(ptCutOff);
  //printf("GetMeasuredFraction: bin range %d %d\n", ptBin, ptProj->GetNbinsX());
  Float_t abovePtCut = ptProj->Integral(ptBin, ptProj->GetNbinsX());
  Float_t all = ptProj->Integral();

  if (all == 0)
    return -1;

  Float_t fraction = abovePtCut / all;

  //printf("GetMeasuredFraction: all %f above %f fraction %f\n", all, abovePtCut, fraction);

  if (debug)
  {
    new TCanvas;
    ptProj->Draw();
  }
  else
    delete ptProj;

  if (debug)
    printf("AlidNdEtaCorrection::GetMeasuredFraction: pt cut off = %f, eta = %f, => fraction = %f\n", ptCutOff, eta, fraction);

  return fraction;
}

//____________________________________________________________________
void AlidNdEtaCorrection::ReduceInformation()
{
  // this function deletes the measured and generated histograms from the corrections to reduce the amount of data
  // in memory

  // these are needed for GetMeasuredFraction(): fTrack2ParticleCorrection->ReduceInformation();
  fVertexRecoCorrection          ->ReduceInformation();
  fTriggerBiasCorrectionMBToINEL ->ReduceInformation();
  fTriggerBiasCorrectionMBToNSD  ->ReduceInformation();
  fTriggerBiasCorrectionMBToND   ->ReduceInformation();
}

