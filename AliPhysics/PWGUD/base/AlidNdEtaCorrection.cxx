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

#include "AlidNdEtaCorrection.h"

#include <AliLog.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TDirectory.h>
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
  fTriggerBiasCorrectionMBToND(0),
  fTriggerBiasCorrectionMBToOnePart(0)
{
  // default constructor
}

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysis)
  : TNamed(name, title),
  fTrack2ParticleCorrection(0),
  fVertexRecoCorrection(0),
  fTriggerBiasCorrectionMBToINEL(0),
  fTriggerBiasCorrectionMBToNSD(0),
  fTriggerBiasCorrectionMBToND(0),
  fTriggerBiasCorrectionMBToOnePart(0)
{
  //
  // constructor
  //

  fTrack2ParticleCorrection = new AliCorrection("Track2Particle", "Track2Particle", analysis);
  fVertexRecoCorrection     = new AliCorrection("VertexReconstruction", "VertexReconstruction", analysis);

  fTriggerBiasCorrectionMBToINEL = new AliCorrection("TriggerBias_MBToINEL", "TriggerBias_MBToINEL", analysis);
  fTriggerBiasCorrectionMBToNSD  = new AliCorrection("TriggerBias_MBToNSD", "TriggerBias_MBToNSD", analysis);
  fTriggerBiasCorrectionMBToND   = new AliCorrection("TriggerBias_MBToND", "TriggerBias_MBToND", analysis);
  fTriggerBiasCorrectionMBToOnePart   = new AliCorrection("TriggerBias_MBToOnePart", "TriggerBias_MBToOnePart", analysis);
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

  if (fTriggerBiasCorrectionMBToOnePart) {
    delete fTriggerBiasCorrectionMBToOnePart;
    fTriggerBiasCorrectionMBToOnePart = 0;
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
  fTriggerBiasCorrectionMBToOnePart->Divide();
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
  TList* collectionTriggerBiasMBToOnePart    = new TList;

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
    collectionTriggerBiasMBToOnePart  ->Add(entry->fTriggerBiasCorrectionMBToOnePart);

    count++;
  }
  fTrack2ParticleCorrection      ->Merge(collectionNtrackToNparticle);
  fVertexRecoCorrection          ->Merge(collectionVertexReco);
  fTriggerBiasCorrectionMBToINEL ->Merge(collectionTriggerBiasMBToINEL);
  fTriggerBiasCorrectionMBToNSD  ->Merge(collectionTriggerBiasMBToNSD);
  fTriggerBiasCorrectionMBToND   ->Merge(collectionTriggerBiasMBToND);
  fTriggerBiasCorrectionMBToOnePart->Merge(collectionTriggerBiasMBToOnePart);

  delete collectionNtrackToNparticle;
  delete collectionVertexReco;
  delete collectionTriggerBiasMBToINEL;
  delete collectionTriggerBiasMBToNSD;
  delete collectionTriggerBiasMBToND;
  delete collectionTriggerBiasMBToOnePart;

  return count+1;
}

//____________________________________________________________________
void AlidNdEtaCorrection::Add(AlidNdEtaCorrection* aCorrectionsToAdd, Float_t c) {
  //
  // adds the measured and generated of aCorrectionsToAdd to measured and generated
  // of all corrections in this

  fTrack2ParticleCorrection      ->Add(aCorrectionsToAdd->GetTrack2ParticleCorrection() ,c);
  fVertexRecoCorrection          ->Add(aCorrectionsToAdd->GetVertexRecoCorrection()     ,c);
  fTriggerBiasCorrectionMBToINEL ->Add(aCorrectionsToAdd->GetTriggerBiasCorrectionINEL(),c);
  fTriggerBiasCorrectionMBToNSD  ->Add(aCorrectionsToAdd->GetTriggerBiasCorrectionNSD() ,c);
  fTriggerBiasCorrectionMBToND   ->Add(aCorrectionsToAdd->GetTriggerBiasCorrectionND()  ,c);
  fTriggerBiasCorrectionMBToOnePart   ->Add(aCorrectionsToAdd->GetTriggerBiasCorrectionOnePart()  ,c);
}

//____________________________________________________________________
void AlidNdEtaCorrection::Scale(Float_t c) 
{
  //
  // scales all contained corrections
  // 

  fTrack2ParticleCorrection      ->Scale(c);
  fVertexRecoCorrection          ->Scale(c);
  fTriggerBiasCorrectionMBToINEL ->Scale(c);
  fTriggerBiasCorrectionMBToNSD  ->Scale(c);
  fTriggerBiasCorrectionMBToND   ->Scale(c);
  fTriggerBiasCorrectionMBToOnePart   ->Scale(c);
}

//____________________________________________________________________
void AlidNdEtaCorrection::Reset(void) {
  //
  // reset all corrections
  // 

  fTrack2ParticleCorrection      ->Reset();
  fVertexRecoCorrection          ->Reset();
  fTriggerBiasCorrectionMBToINEL ->Reset();
  fTriggerBiasCorrectionMBToNSD  ->Reset();
  fTriggerBiasCorrectionMBToND   ->Reset();
  fTriggerBiasCorrectionMBToOnePart   ->Reset();
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
  fTriggerBiasCorrectionMBToOnePart   ->LoadHistograms();

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
  fTriggerBiasCorrectionMBToOnePart->SaveHistograms();

  gDirectory->cd("..");
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the corrections
  //

  fTrack2ParticleCorrection     ->DrawHistograms();
  fVertexRecoCorrection         ->DrawHistograms();
  fTriggerBiasCorrectionMBToINEL->DrawHistograms();
  fTriggerBiasCorrectionMBToNSD ->DrawHistograms();
  fTriggerBiasCorrectionMBToND  ->DrawHistograms();
  fTriggerBiasCorrectionMBToOnePart  ->DrawHistograms();
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawOverview(const char* canvasName)
{
  //
  // call the DrawOverview histogram method of the corrections
  //

  fTrack2ParticleCorrection     ->DrawOverview(canvasName);
  fVertexRecoCorrection         ->DrawOverview(canvasName);
  fTriggerBiasCorrectionMBToINEL->DrawOverview(canvasName);
  fTriggerBiasCorrectionMBToNSD ->DrawOverview(canvasName);
  fTriggerBiasCorrectionMBToND  ->DrawOverview(canvasName);
  fTriggerBiasCorrectionMBToOnePart  ->DrawOverview(canvasName);
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillMCParticle(Float_t vtx, Float_t eta, Float_t pt, Bool_t trigger, Bool_t vertex, Int_t processType)
{
  // fills a particle in the corrections
  // it is filled in generated or measured depending of the flags

  fTriggerBiasCorrectionMBToINEL->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if ((processType & AliPWG0Helper::kSD) == 0)
    fTriggerBiasCorrectionMBToNSD->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (processType & AliPWG0Helper::kND )
    fTriggerBiasCorrectionMBToND->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (processType & AliPWG0Helper::kOnePart)
    fTriggerBiasCorrectionMBToOnePart->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (!trigger)
    return;

  fTriggerBiasCorrectionMBToINEL->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTriggerBiasCorrectionMBToNSD->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTriggerBiasCorrectionMBToND->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTriggerBiasCorrectionMBToOnePart->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fVertexRecoCorrection->GetTrackCorrection()->FillGene(vtx, eta, pt);

  if (!vertex)
    return;

  fVertexRecoCorrection->GetTrackCorrection()->FillMeas(vtx, eta, pt);
  fTrack2ParticleCorrection->GetTrackCorrection()->FillGene(vtx, eta, pt);
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillTrackedParticle(Float_t vtx, Float_t eta, Float_t pt, Double_t weight)
{
  // fills a tracked particle in the corrections

  fTrack2ParticleCorrection->GetTrackCorrection()->FillMeas(vtx, eta, pt, weight);
}

//____________________________________________________________________
void AlidNdEtaCorrection::FillEvent(Float_t vtx, Float_t n, Bool_t trigger, Bool_t vertex, Int_t processType)
{
  // fills an event int he correction
  // it is filled in generated or measured depending of the flags

  fTriggerBiasCorrectionMBToINEL->GetEventCorrection()->FillGene(vtx, n);

  if ((processType & AliPWG0Helper::kSD) == 0)
    fTriggerBiasCorrectionMBToNSD->GetEventCorrection()->FillGene(vtx, n);

  if (processType & AliPWG0Helper::kND )
    fTriggerBiasCorrectionMBToND->GetEventCorrection()->FillGene(vtx, n);

  if (processType & AliPWG0Helper::kOnePart)
    fTriggerBiasCorrectionMBToOnePart->GetEventCorrection()->FillGene(vtx, n);

  if (!trigger)
    return;

  fTriggerBiasCorrectionMBToINEL->GetEventCorrection()->FillMeas(vtx, n);
  fTriggerBiasCorrectionMBToNSD->GetEventCorrection()->FillMeas(vtx, n);
  fTriggerBiasCorrectionMBToND->GetEventCorrection()->FillMeas(vtx, n);
  fTriggerBiasCorrectionMBToOnePart->GetEventCorrection()->FillMeas(vtx, n);
  fVertexRecoCorrection->GetEventCorrection()->FillGene(vtx, n);

  if (!vertex)
    return;

  fVertexRecoCorrection->GetEventCorrection()->FillMeas(vtx, n);
}

//____________________________________________________________________
Float_t AlidNdEtaCorrection::GetMeasuredFraction(CorrectionType correctionType, Float_t ptCutOff, Float_t eta, Int_t vertexBegin, Int_t vertexEnd, Bool_t debug)
{
  // calculates the fraction of particles measured (some are missed due to the pt cut off)
  //
  // uses the generated particle histogram from the correction passed, e.g. pass GetTrack2ParticleCorrection()

  if (!GetCorrection(correctionType))
    return -1;

  const TH3* generated = GetCorrection(correctionType)->GetTrackCorrection()->GetGeneratedHistogram();

  // find eta borders, if eta is negative assume -0.8 ... 0.8
  Int_t etaBegin = 0;
  Int_t etaEnd = 0;
  const TAxis * xax = generated->GetXaxis();
  const TAxis * yax = generated->GetYaxis();
  if (eta < -99)
  {
    etaBegin = yax->FindFixBin(-0.8);
    etaEnd = yax->FindFixBin(0.8);
  }
  else
  {
    etaBegin = yax->FindFixBin(eta);
    etaEnd = etaBegin;
  }

  if (vertexBegin == -1)
    vertexBegin = xax->FindFixBin(-9.99);

  if (vertexEnd == -1)
    vertexEnd = xax->FindFixBin(9.99);

  TH1D* ptProj = dynamic_cast<TH1D*> (generated->ProjectionZ(Form("%s_pt", generated->GetName()), vertexBegin, vertexEnd, etaBegin, etaEnd));
  //printf("GetMeasuredFraction: bin range %d %d %d %d\n", vertexBegin, vertexEnd, etaBegin, etaEnd);
  ptProj->GetXaxis()->SetTitle(generated->GetZaxis()->GetTitle());

  Int_t ptBin = ptProj->FindBin(ptCutOff);
  //printf("GetMeasuredFraction: bin range %d %d\n", ptBin, ptProj->GetNbinsX());
  Float_t abovePtCut = ptProj->Integral(ptBin, ptProj->GetNbinsX()+1);
  Float_t all = ptProj->Integral(1, ptProj->GetNbinsX()+1);

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
TH1* AlidNdEtaCorrection::GetMeasuredEventFraction(CorrectionType correctionType, Int_t multCut)
{
  // calculates the fraction of events above multCut (but including it)
  //
  // uses the generated event histogram from the correction passed, e.g. pass GetTrack2ParticleCorrection()

  if (!GetCorrection(correctionType))
    return 0;

  const TH2* generated = GetCorrection(correctionType)->GetEventCorrection()->GetGeneratedHistogram();

  TH1* allEvents = generated->ProjectionX(Form("%s_all", generated->GetName()), 1, generated->GetNbinsY());
  TH1* aboveEvents = generated->ProjectionX(Form("%s_above", generated->GetName()), generated->GetYaxis()->FindFixBin(multCut), generated->GetNbinsY());
  
  aboveEvents->Divide(aboveEvents, allEvents, 1, 1, "B");

  return aboveEvents;  
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
  fTriggerBiasCorrectionMBToOnePart   ->ReduceInformation();
}

//____________________________________________________________________
AliCorrection* AlidNdEtaCorrection::GetCorrection(CorrectionType correctionType)
{
  // returns the given correction

  switch (correctionType)
  {
    case kNone : return 0;
    case kTrack2Particle : return fTrack2ParticleCorrection;
    case kVertexReco :     return fVertexRecoCorrection;
    case kINEL :           return fTriggerBiasCorrectionMBToINEL;
    case kNSD :            return fTriggerBiasCorrectionMBToNSD;
    case kND :             return fTriggerBiasCorrectionMBToND;
    case kOnePart :             return fTriggerBiasCorrectionMBToOnePart;
  }

  return 0;
}
