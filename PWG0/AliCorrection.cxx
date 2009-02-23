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
#include <TMath.h>

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
AliCorrection::AliCorrection(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysisMode) : TNamed(name, title),
  fEventCorr(0),
  fTrackCorr(0)
{
  // constructor initializing tnamed

  Float_t* binLimitsPt = 0;
  Int_t nBinsPt = 0;

  // different binnings, better solution could be anticipated...
  if (analysisMode == AliPWG0Helper::kTPC || analysisMode == AliPWG0Helper::kTPCITS)
  {
    static Float_t binLimitsPtTmp[] = {0.0, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0, 100.0};
    binLimitsPt = (Float_t*) binLimitsPtTmp;
    nBinsPt = 28;
  }
  else if (analysisMode == AliPWG0Helper::kSPD)
  {
    static Float_t binLimitsPtTmp[] = {-0.5, 100.0};
    binLimitsPt = (Float_t*) binLimitsPtTmp;
    nBinsPt = 1;
  }

  if (!binLimitsPt)
  {
    Printf("FATAL: Invalid binning");
    return;
  }

  // third axis for track histogram
  Int_t nBinsN2 = 1;
  Float_t binLimitsN2[]   = {-0.5, 1000};
  const char* title3 = "Ntracks";
  //Int_t nBinsN2 = 21;
  //Float_t binLimitsN2[]   = {-0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 1000.5};
  // phi
  //Int_t nBinsN2 = 36;
  //Float_t binLimitsN2[]   = {0.000, 0.175, 0.349, 0.524, 0.698, 0.873, 1.047, 1.222, 1.396, 1.571, 1.745, 1.920, 2.094, 2.269, 2.443, 2.618, 2.793, 2.967, 3.142, 3.316, 3.491, 3.665, 3.840, 4.014, 4.189, 4.363, 4.538, 4.712, 4.887, 5.061, 5.236, 5.411, 5.585, 5.760, 5.934, 6.109, TMath::Pi() * 2};
  //const char* title3 = "#phi";

  // mult axis for event histogram
  Int_t nBinsN = 22;
  Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
  //Int_t nBinsN = 52;
  //Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5, 100.5, 300.5};

  //Float_t binLimitsVtx[] = {-20,-15,-10,-6,-3,0,3,6,10,15,20};
  //Float_t binLimitsVtx[] = {-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  Float_t binLimitsVtx[] = {-30,-25,-20,-15,-10,-8,-6,-4,-2,0,2,4,6,8,10,15,20,25,30};
  /*Float_t binLimitsEta[] = {-3.0,-2.8,-2.6,-2.4,-2.2,
                            -2.0,-1.8,-1.6,-1.4,-1.2,
                            -1.0,-0.8,-0.6,-0.4,-0.2, 0.0,
                             0.2, 0.4, 0.6, 0.8, 1.0,
                             1.2, 1.4, 1.6, 1.8, 2.0,
                             2.2, 2.4, 2.6, 2.8, 3.0};*/
  Float_t binLimitsEta[] = { -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0 };
//  Float_t binLimitsEta[] = {-2,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
//  Float_t binLimitsEta[] = {-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
//  Float_t binLimitsEta[] = { -7.0, -6.9, -6.8, -6.7, -6.6, -6.5, -6.4, -6.3, -6.2, -6.1, -6.0, -5.9, -5.8, -5.7, -5.6, -5.5, -5.4, -5.3, -5.2, -5.1, -5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0 };

  fEventCorr = new AliCorrectionMatrix2D("EventCorrection", Form("%s EventCorrection", fTitle.Data()), 18, binLimitsVtx, nBinsN, binLimitsN);

  if (analysisMode == AliPWG0Helper::kSPD)
  {
    TH3F* dummyBinning = new TH3F("dummyBinning","dummyBinning",18, binLimitsVtx, 60, binLimitsEta, nBinsN2, binLimitsN2);
    fTrackCorr = new AliCorrectionMatrix3D("TrackCorrection", Form("%s TrackCorrection", fTitle.Data()), dummyBinning);
    fTrackCorr->SetAxisTitles("vtx-z (cm)", "#eta", title3);
    delete dummyBinning;
  }
  else
  {
    TH3F* dummyBinning = new TH3F("dummyBinning","dummyBinning",18, binLimitsVtx, 60, binLimitsEta , nBinsPt, binLimitsPt);
    fTrackCorr = new AliCorrectionMatrix3D("TrackCorrection", Form("%s TrackCorrection", fTitle.Data()), dummyBinning);
    fTrackCorr->SetAxisTitles("vtx-z (cm)", "#eta", "p_{T} (GeV/c)");
    delete dummyBinning;
  }


  fEventCorr->SetAxisTitles("vtx-z (cm)", "Ntracks");
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

  Int_t emptyBins = fTrackCorr->CheckEmptyBins(-9.99, 9.99, -0.79, 0.79, 0.2, 4.9);
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

void AliCorrection::DrawOverview(const char* canvasName)
{
  // draw projection of the corrections
  //   to the 3 axis of fTrackCorr and 2 axis of fEventCorr

  TString canvasNameTmp(GetName());
  if (canvasName)
    canvasNameTmp = canvasName;

  TCanvas* canvas = new TCanvas(canvasNameTmp, canvasNameTmp, 1200, 800);
  canvas->Divide(3, 2);

  if (fTrackCorr) {
    canvas->cd(1);
    fTrackCorr->Get1DCorrectionHistogram("x", 0.3, 5, -1, 1)->DrawCopy()->GetYaxis()->SetRangeUser(0, 10);

    canvas->cd(2);
    fTrackCorr->Get1DCorrectionHistogram("y", 0.3, 5, 0, 0)->DrawCopy()->GetYaxis()->SetRangeUser(0, 10);

    canvas->cd(3);
    fTrackCorr->Get1DCorrectionHistogram("z", 0, -1, -1, 1)->DrawCopy()->GetYaxis()->SetRangeUser(0, 10);
  }

  if (fEventCorr)
  {
    canvas->cd(4);
    fEventCorr->Get1DCorrectionHistogram("x")->DrawCopy();

    canvas->cd(5);
    fEventCorr->Get1DCorrectionHistogram("y")->DrawCopy()->GetXaxis()->SetRangeUser(0, 30);
  }
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
    TH2* hist = fEventCorr->GetMeasuredHistogram();
    for (Int_t x = 1; x <= hist->GetNbinsX(); ++x)
      fEventCorr->GetGeneratedHistogram()->SetBinContent(x, hist->GetNbinsY() + 1, hist->GetBinContent(x, hist->GetNbinsY() + 1));
  }

  if (fTrackCorr)
    fTrackCorr->Multiply();
}

//____________________________________________________________________
void AliCorrection::Scale(Double_t factor)
{
  // scales the two contained corrections

  fEventCorr->Scale(factor);
  fTrackCorr->Scale(factor);
}

//____________________________________________________________________
void AliCorrection::PrintStats(Float_t zRange, Float_t etaRange, Float_t ptCut)
{
  // prints statistics and effective correction factors

  Printf("AliCorrection::PrintInfo: Values in |eta| < %.2f, |vtx-z| < %.2f and pt > %.2f:", etaRange, zRange, ptCut);

  // prevent to be on bin border
  zRange -= 0.1;
  etaRange -= 0.1;

  TH3* measured = GetTrackCorrection()->GetMeasuredHistogram();
  TH3* generated = GetTrackCorrection()->GetGeneratedHistogram();

  TH2* measuredEvents = GetEventCorrection()->GetMeasuredHistogram();
  TH2* generatedEvents = GetEventCorrection()->GetGeneratedHistogram();

  Float_t tracksM = measured->Integral(measured->GetXaxis()->FindBin(-zRange), measured->GetXaxis()->FindBin(zRange), measured->GetYaxis()->FindBin(-etaRange), measured->GetYaxis()->FindBin(etaRange), measured->GetZaxis()->FindBin(ptCut), measured->GetZaxis()->GetNbins());
  Float_t tracksG = generated->Integral(generated->GetXaxis()->FindBin(-zRange), generated->GetXaxis()->FindBin(zRange), generated->GetYaxis()->FindBin(-etaRange), generated->GetYaxis()->FindBin(etaRange), generated->GetZaxis()->FindBin(ptCut), generated->GetZaxis()->GetNbins());

  Float_t eventsM = measuredEvents->Integral(measuredEvents->GetXaxis()->FindBin(-zRange), measuredEvents->GetXaxis()->FindBin(zRange), 1, measuredEvents->GetNbinsY());
  Float_t eventsG = generatedEvents->Integral(generatedEvents->GetXaxis()->FindBin(-zRange), generatedEvents->GetXaxis()->FindBin(zRange), 1, generatedEvents->GetNbinsY());

  Printf("tracks measured: %.1f tracks generated: %.1f, events measured: %.1f, events generated %.1f", tracksM, tracksG, eventsM, eventsG);

  if (tracksM > 0)
    Printf("Effective average correction factor for TRACKS: %.3f", tracksG / tracksM);
  if (eventsM > 0)
    Printf("Effective average correction factor for EVENTS: %.3f", eventsG / eventsM);

  if (eventsM > 0 && eventsG > 0)
  {
    // normalize to number of events;
    tracksM /= eventsM;
    tracksG /= eventsG;

    Printf("%.2f tracks/event measured, %.2f tracks/event after correction --> effective average correction factor is %.3f (pt cut %.2f GeV/c)", tracksM, tracksG, (tracksM > 0) ? (tracksG / tracksM) : -1, ptCut);
  }
}

//____________________________________________________________________
void AliCorrection::PrintInfo(Float_t ptCut)
{
  // prints some stats

  TH3* measured = GetTrackCorrection()->GetMeasuredHistogram();
  TH3* generated = GetTrackCorrection()->GetGeneratedHistogram();

  TH2* measuredEvents = GetEventCorrection()->GetMeasuredHistogram();
  TH2* generatedEvents = GetEventCorrection()->GetGeneratedHistogram();

  Printf("AliCorrection::PrintInfo: Whole phasespace:");

  Printf("tracks measured: %.1f tracks generated: %.1f, events measured: %.1f, events generated %.1f", measured->Integral(), generated->Integral(), measuredEvents->Integral(), generatedEvents->Integral());

  Printf("Example centered bin: tracks measured: %.1f tracks generated: %.1f, events measured: %.1f, events generated %.1f", measured->GetBinContent(measured->GetNbinsX() / 2, measured->GetNbinsY() / 2, measured->GetNbinsZ() / 2), generated->GetBinContent(measured->GetNbinsX() / 2, measured->GetNbinsY() / 2, measured->GetNbinsZ() / 2), measuredEvents->GetBinContent(measuredEvents->GetNbinsX() / 2, measuredEvents->GetNbinsY() / 2), generatedEvents->GetBinContent(measuredEvents->GetNbinsX() / 2, measuredEvents->GetNbinsY() / 2));

  PrintStats(10, 1.0, ptCut);
  PrintStats(10, 1.5, ptCut);
}
