/* $Id$ */

// ------------------------------------------------------
//
// Class to handle corrections.
//
// ------------------------------------------------------
//

#include <TFile.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TH2F.h>

#include <AliLog.h>
#include "AliCorrectionMatrix2D.h"
#include "AliCorrectionMatrix3D.h"

#include "AliCorrection.h"

//____________________________________________________________________
ClassImp(AliCorrection)

//____________________________________________________________________
AliCorrection::AliCorrection() : TNamed(),
  fEventCorr(0),
  fTrackCorr(0)
{
  // default constructor
}

//____________________________________________________________________
AliCorrection::AliCorrection(const Char_t* name, const Char_t* title) : TNamed(name, title),
  fEventCorr(0),
  fTrackCorr(0)
{
  // constructor initializing tnamed

  Float_t binLimitsPt[] = {0.0, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0, 100.0};
  Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
  //Float_t binLimitsVtx[] = {-20,-15,-10,-6,-3,0,3,6,10,15,20};
  //Float_t binLimitsVtx[] = {-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  Float_t binLimitsVtx[] = {-20,-15,-10,-8,-6,-4,-2,0,2,4,6,8,10,15,20};
  Float_t binLimitsEta[] = {-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
			    0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};

  TH3F* dummyBinning = new TH3F("","",14, binLimitsVtx, 40, binLimitsEta , 28, binLimitsPt);

  fEventCorr = new AliCorrectionMatrix2D("EventCorrection", Form("%s EventCorrection", title), 14, binLimitsVtx, 22, binLimitsN);
  fTrackCorr = new AliCorrectionMatrix3D("TrackCorrection", Form("%s TrackCorrection", title), dummyBinning);

  delete dummyBinning;

  fEventCorr->SetAxisTitles("vtx z [cm]", "Ntracks");
  fTrackCorr->SetAxisTitles("vtx z [cm]", "#eta", "p_{T} [GeV/c]");
}

//____________________________________________________________________
AliCorrection::AliCorrection(const AliCorrection& c) : TNamed(c),
  fEventCorr(0),
  fTrackCorr(0)
{
  // copy constructor
  ((AliCorrection &)c).Copy(*this);
}

//____________________________________________________________________
AliCorrection::~AliCorrection()
{
  //
  // destructor
  //

  if (fEventCorr)
  {
    delete fEventCorr;
    fEventCorr = 0;
  }

  if (fTrackCorr)
  {
    delete fTrackCorr;
    fTrackCorr = 0;
  }
}

//____________________________________________________________________
AliCorrection &AliCorrection::operator=(const AliCorrection &c)
{
  // assigment operator

  if (this != &c)
    ((AliCorrection &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliCorrection::Copy(TObject& c) const
{
  // copy function

  AliCorrection& target = (AliCorrection &) c;

  if (fEventCorr)
    target.fEventCorr = dynamic_cast<AliCorrectionMatrix2D*> (fEventCorr->Clone());

  if (fTrackCorr)
    target.fTrackCorr = dynamic_cast<AliCorrectionMatrix3D*> (fTrackCorr->Clone());
}

//____________________________________________________________________
Long64_t AliCorrection::Merge(TCollection* list)
{
  // Merge a list of AliCorrection objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionEvent = new TList;
  TList* collectionTrack = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliCorrection* entry = dynamic_cast<AliCorrection*> (obj);
    if (entry == 0) 
      continue;

    collectionEvent->Add(entry->fEventCorr);
    collectionTrack->Add(entry->fTrackCorr);

    count++;
  }
  fEventCorr->Merge(collectionEvent);
  fTrackCorr->Merge(collectionTrack);

  delete collectionEvent;
  delete collectionTrack;

  return count+1;
}

//____________________________________________________________________
void AliCorrection::Divide()
{
  //
  // divide the histograms to get the correction
  //
  
  if (!fEventCorr || !fTrackCorr)
    return;
    
  fEventCorr->Divide();
  fTrackCorr->Divide();

  Int_t emptyBins = fTrackCorr->CheckEmptyBins(-9.99, 9.99, -0.79, 0.79, 0.3, 9.9);
  printf("INFO: In the central region the track correction of %s has %d empty bins\n", GetName(), emptyBins);
}

//____________________________________________________________________
void AliCorrection::Add(AliCorrection* aCorrectionToAdd, Float_t c)
{
  //
  // add to measured and generated the measured and generated of aCorrectionToAdd
  // with the weight c

  fEventCorr->Add(aCorrectionToAdd->GetEventCorrection(),c);
  fTrackCorr->Add(aCorrectionToAdd->GetTrackCorrection(),c);
}


//____________________________________________________________________
Bool_t AliCorrection::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!fEventCorr || !fTrackCorr)
    return kFALSE;

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  Bool_t success = fEventCorr->LoadHistograms();
  success &= fTrackCorr->LoadHistograms();

  gDirectory->cd("..");

  return success;
}

//____________________________________________________________________
void AliCorrection::SaveHistograms()
{
  //
  // saves the histograms in a directory with the name of this object (GetName)
  //
  
  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());

  if (fEventCorr)
    fEventCorr->SaveHistograms();

  if (fTrackCorr)
    fTrackCorr->SaveHistograms();
    
  gDirectory->cd("..");
}

//____________________________________________________________________
void AliCorrection::ReduceInformation()
{
  // this function deletes the measured and generated histograms to reduce the amount of data
  // in memory

  if (!fEventCorr || !fTrackCorr)
    return;

  fEventCorr->ReduceInformation();
  fTrackCorr->ReduceInformation();
}

//____________________________________________________________________
void AliCorrection::Reset(Option_t* option)
{
  // resets the histograms

  if (fEventCorr)
    fEventCorr->Reset(option);

  if (fTrackCorr)
    fTrackCorr->Reset(option);
}

//____________________________________________________________________
void AliCorrection::DrawHistograms(const Char_t* name)
{
  // draws the corrections

  if (!name)
    name = GetName();

  if (fEventCorr)
    fEventCorr->DrawHistograms(Form("%s event", name));

  if (fTrackCorr)
    fTrackCorr->DrawHistograms(Form("%s track", name));
}

//____________________________________________________________________
void AliCorrection::SetCorrectionToUnity()
{
  // set the corrections to unity

  if (fEventCorr)
    fEventCorr->SetCorrectionToUnity();

  if (fTrackCorr)
    fTrackCorr->SetCorrectionToUnity();
}

//____________________________________________________________________
void AliCorrection::Multiply()
{
  // call Multiply

  if (fEventCorr)
  {
    fEventCorr->Multiply();
    // now we manually copy the overflow bin of the y axis (multiplicity) over. This is important to get the event count correct
    TH2F* hist = fEventCorr->GetMeasuredHistogram();
    for (Int_t x = 1; x <= hist->GetNbinsX(); ++x)
      fEventCorr->GetGeneratedHistogram()->SetBinContent(x, hist->GetNbinsY() + 1, hist->GetBinContent(x, hist->GetNbinsY() + 1));
  }

  if (fTrackCorr)
    fTrackCorr->Multiply();
}
