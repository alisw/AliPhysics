/* $Id$ */

#include "AlidNdEtaCorrection.h"

#include <TCanvas.h>
#include <TH3F.h>
#include <TH1D.h>

//____________________________________________________________________
ClassImp(AlidNdEtaCorrection)

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection()
  : TNamed(),
  fTrack2ParticleCorrection(0),
  fVertexRecoCorrection(0),
  fTriggerCorrection(0),
  fTriggerBiasCorrection(0)
  //fNEvents(0),
  //fNTriggeredEvents(0)
{
  // default constructor
}

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection(const Char_t* name, const Char_t* title)
  : TNamed(name, title),
  fTrack2ParticleCorrection(0),
  fVertexRecoCorrection(0),
  fTriggerCorrection(0),
  fTriggerBiasCorrection(0)
  //fNEvents(0),
  //fNTriggeredEvents(0)
{
  // constructor
  //

  Float_t binLimitsPt[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 10.0, 100.0};

  TString matrixName;
  matrixName.Form("%s_nTrackToNPart", name);

  fTrack2ParticleCorrection = new AliCorrectionMatrix3D(matrixName, matrixName, 40, -20, 20, 60, -6, 6, 14, binLimitsPt);

  Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
			    10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
  Float_t binLimitsVtx[] = {-20,-15,-10,-6,-3,0,3,6,10,15,20};

  matrixName.Form("%s_vtxReco", name);
  fVertexRecoCorrection        = new AliCorrectionMatrix2D(matrixName, matrixName, 10,binLimitsVtx ,22,binLimitsN);

  matrixName.Form("%s_trigger", name);
  fTriggerCorrection           = new AliCorrectionMatrix2D(matrixName, matrixName, 10,binLimitsVtx ,22,binLimitsN);

  matrixName.Form("%s_triggerBias", name);
  fTriggerBiasCorrection       = new AliCorrectionMatrix2D(matrixName, matrixName, 120,-6,6,100, 0, 10);

  fTrack2ParticleCorrection ->SetAxisTitles("vtx z [cm]", "#eta", "p_{T} [GeV/c]");
  fVertexRecoCorrection        ->SetAxisTitles("vtx z [cm]", "Ntracks");
  fTriggerCorrection        ->SetAxisTitles("vtx z [cm]", "Ntracks");

  fTriggerBiasCorrection       ->SetAxisTitles("#eta", "p_{T} [GeV/c]");
}

//____________________________________________________________________
AlidNdEtaCorrection::~AlidNdEtaCorrection()
{
  // destructor

  if (fTrack2ParticleCorrection)
  {
    delete fTrack2ParticleCorrection;
    fTrack2ParticleCorrection = 0;
  }

  if (fVertexRecoCorrection)
  {
    delete fVertexRecoCorrection;
    fVertexRecoCorrection = 0;
  }

  if (fTriggerCorrection)
  {
    delete fTriggerCorrection;
    fTriggerCorrection = 0;
  }

  if (fTriggerBiasCorrection)
  {
    delete fTriggerBiasCorrection;
    fTriggerBiasCorrection = 0;
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

  TH3F* hist = fTrack2ParticleCorrection->GetCorrectionHistogram();
  Int_t emptyBins = 0;
  for (Int_t x=hist->GetXaxis()->FindBin(-10); x<=hist->GetXaxis()->FindBin(10); ++x)
    for (Int_t y=hist->GetYaxis()->FindBin(-0.8); y<=hist->GetYaxis()->FindBin(0.8); ++y)
      for (Int_t z=hist->GetZaxis()->FindBin(0.3); z<=hist->GetZaxis()->FindBin(9.9); ++z)
        if (hist->GetBinContent(x, y, z) == 0)
        {
          printf("Empty bin in fTrack2ParticleCorrection at vtx = %f, eta = %f, pt = %f\n", hist->GetXaxis()->GetBinCenter(x), hist->GetYaxis()->GetBinCenter(y), hist->GetZaxis()->GetBinCenter(z));
          ++emptyBins;
        }

  printf("INFO: In the central region fTrack2ParticleCorrection has %d empty bins\n", emptyBins);

  fVertexRecoCorrection->Divide();
  fTriggerCorrection->Divide();

  if (GetNevents() <= 0)
  {
    printf("ERROR: GetNevents() returns 0. Cannot scale histogram. Skipping processing of fTriggerBiasCorrection\n");
    return;
  }
  fTriggerBiasCorrection->GetMeasuredHistogram()->Scale(Double_t(GetNtriggeredEvents())/Double_t(GetNevents()));
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
  TList* collectionTrigger           = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AlidNdEtaCorrection* entry = dynamic_cast<AlidNdEtaCorrection*> (obj);
    if (entry == 0)
      continue;

    collectionNtrackToNparticle ->Add(entry->GetTrack2ParticleCorrection());
    collectionVertexReco        ->Add(entry->GetVertexRecoCorrection());
    collectionTriggerBias        ->Add(entry->GetTriggerBiasCorrection());
    collectionTrigger        ->Add(entry->GetTriggerCorrection());

    count++;

    //fNEvents += entry->fNEvents;
    //fNTriggeredEvents += entry->fNTriggeredEvents;
  }
  fTrack2ParticleCorrection ->Merge(collectionNtrackToNparticle);
  fVertexRecoCorrection        ->Merge(collectionVertexReco);
  fTriggerBiasCorrection        ->Merge(collectionTriggerBias);
  fTriggerCorrection        ->Merge(collectionTrigger);

  delete collectionNtrackToNparticle;
  delete collectionVertexReco;
  delete collectionTriggerBias;
  delete collectionTrigger;

  return count+1;
}

//____________________________________________________________________
Bool_t
AlidNdEtaCorrection::LoadHistograms(const Char_t* fileName, const Char_t* dir) {
  //
  // loads the histograms
  //

  fTrack2ParticleCorrection ->LoadHistograms(fileName, dir);
  fVertexRecoCorrection        ->LoadHistograms(fileName, dir);
  fTriggerCorrection        ->LoadHistograms(fileName, dir);
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

  fTrack2ParticleCorrection->SaveHistograms();
  fVertexRecoCorrection->SaveHistograms();
  fTriggerCorrection->SaveHistograms();
  fTriggerBiasCorrection->SaveHistograms();

  gDirectory->cd("../");
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the two AliCorrectionMatrix2D objects

  fTrack2ParticleCorrection ->DrawHistograms();
  fVertexRecoCorrection        ->DrawHistograms();
  fTriggerCorrection        ->DrawHistograms();
  fTriggerBiasCorrection       ->DrawHistograms();

}

//____________________________________________________________________
Float_t AlidNdEtaCorrection::GetMeasuredFraction(Float_t ptCutOff, Float_t eta, Bool_t debug)
{
  // calculates the fraction of particles measured (some are missed due to the pt cut off)
  // uses the generated particle histogram from fTrack2ParticleCorrection

  const TH3F* generated = fTrack2ParticleCorrection->GetGeneratedHistogram();

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

  Int_t vertexBegin = generated->GetXaxis()->FindBin(-10);
  Int_t vertexEnd = generated->GetXaxis()->FindBin(10);

  TH1D* ptProj = dynamic_cast<TH1D*> (generated->ProjectionZ(Form("%s_pt", generated->GetName()), vertexBegin, vertexEnd, etaBegin, etaEnd));
  ptProj->GetXaxis()->SetTitle(generated->GetZaxis()->GetTitle());

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
  else
    delete ptProj;

  if (debug)
    printf("AlidNdEtaCorrection::GetMeasuredFraction: pt cut off = %f, eta = %f, => fraction = %f\n", ptCutOff, eta, fraction);

  return fraction;
}

void AlidNdEtaCorrection::ReduceInformation()
{
  // this function deletes the measured and generated histograms from the corrections to reduce the amount of data
  // in memory

  // these are needed for GetMeasuredFraction(): fTrack2ParticleCorrection->ReduceInformation();
  fVertexRecoCorrection->ReduceInformation();
  fTriggerCorrection->ReduceInformation();
  fTriggerBiasCorrection->ReduceInformation();
}

Int_t AlidNdEtaCorrection::GetNevents()
{
  // Return total number of events used for this correction map (taken from fTriggerCorrection)

  if (!fTriggerCorrection)
    return -1;

  return (Int_t) fTriggerCorrection->GetGeneratedHistogram()->Integral(0, fTriggerCorrection->GetGeneratedHistogram()->GetNbinsX()+1);
}

Int_t AlidNdEtaCorrection::GetNtriggeredEvents()
{
  // Return total number of triggered events used for this correction map (taken from fTriggerCorrection)

  if (!fTriggerCorrection)
    return -1;

  return (Int_t) fTriggerCorrection->GetMeasuredHistogram()->Integral(0, fTriggerCorrection->GetMeasuredHistogram()->GetNbinsX()+1);
}
