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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class maintains a histogram that is used to monitor the quality     //
//  of the recorded data.                                                    //
//  The histogram is created and filled by a sub class of AliMonitor.        //
//  It can be compared to a reference histogram. For the comparison a        //
//  maximal deviation (in standard deviations) can be specified.             //
//  The bins where the maximal deviation is exceeded are drawn in red.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TProfile.h>
#include <TH2.h>
#include <TVirtualPad.h>

#include "AliLog.h"

#include "AliMonitorHisto.h"


ClassImp(AliMonitorHisto) 


Int_t   AliMonitorHisto::fgNHistosMax = 10;


//_____________________________________________________________________________
AliMonitorHisto::AliMonitorHisto() :
  AliMonitorPlot(),
  fHisto(NULL),
  fHistoList(),
  fNHistos(0),
  fHistoRun(NULL),
  fHistoDraw(NULL),
  fHistoRef(NULL),
  fHistoCompare(NULL),
  fNorm(kNormNone)
{
// default contructor

}

//_____________________________________________________________________________
AliMonitorHisto::AliMonitorHisto(const AliMonitorHisto& histo) :
  AliMonitorPlot(histo),
  fHisto(NULL),
  fHistoList(),
  fNHistos(histo.fNHistos),
  fHistoRun(NULL),
  fHistoDraw(NULL),
  fHistoRef(NULL),
  fHistoCompare(NULL),
  fNorm(histo.fNorm)
{
// copy constructor

  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  if (histo.fHisto) fHisto = (TH1*) histo.fHisto->Clone();
  TObjLink* link = histo.fHistoList.FirstLink();
  for (Int_t i = 0; i < fNHistos; i++) {
    fHistoList.Add(link->GetObject()->Clone());
    link = link->Next();
  }
  if (histo.fHistoRun) fHistoRun = (TH1*) histo.fHistoRun->Clone();
  if (histo.fHistoRef) fHistoRef = (TH1*) histo.fHistoRef->Clone();
  TH1::AddDirectory(addStatus);
}

//_____________________________________________________________________________
AliMonitorHisto& AliMonitorHisto::operator =(const AliMonitorHisto& histo)
{
// assignment operator

  AliMonitorPlot::operator =(histo);

  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHisto = NULL;
  if (histo.fHisto) fHisto = (TH1*) histo.fHisto->Clone();
  fNHistos = histo.fNHistos;
  TObjLink* link = histo.fHistoList.FirstLink();
  for (Int_t i = 0; i < fNHistos; i++) {
    fHistoList.Add(link->GetObject()->Clone());
    link = link->Next();
  }
  fHistoRun = NULL;
  if (histo.fHistoRun) fHistoRun = (TH1*) histo.fHistoRun->Clone();
  fHistoDraw = NULL;
  fHistoRef = NULL;
  if (histo.fHistoRef) fHistoRef = (TH1*) histo.fHistoRef->Clone();
  fHistoCompare = NULL;
  fNorm = histo.fNorm;
  TH1::AddDirectory(addStatus);

  return *this;
}

//_____________________________________________________________________________
AliMonitorHisto::AliMonitorHisto(TH1* histo, ENorm norm) :
  AliMonitorPlot(histo->GetName(), histo->GetTitle()),
  fHisto(histo),
  fHistoList(),
  fNHistos(0),
  fHistoRun(NULL),
  fHistoDraw(NULL),
  fHistoRef(NULL),
  fHistoCompare(NULL),
  fNorm(norm)
{
// create a monitor histogram from the given histogram

  if (histo->GetDimension() > 2) {
    AliFatal("3 dimensional histograms are not supported");
  }

  histo->SetDirectory(NULL);
  histo->Reset();
  fHisto->Sumw2();
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoRun = (TH1*) histo->Clone();
  TH1::AddDirectory(addStatus);
}

//_____________________________________________________________________________
AliMonitorHisto::~AliMonitorHisto()
{
// delete all histograms

  if (fHisto) delete fHisto;
  fHistoList.Delete();
  if (fHistoRun) delete fHistoRun;
  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
}


//_____________________________________________________________________________
void AliMonitorHisto::SetReference(TH1* ref)
{
// set the reference histogram for comparison

  if (fHistoRef) delete fHistoRef;
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(addStatus);
  fHistoRef = (TH1*) ref->Clone();
  if (fHistoCompare) fHistoCompare->Reset();
}

//_____________________________________________________________________________
void AliMonitorHisto::SetReference(AliMonitorPlot* ref)
{
// set the reference histogram for comparison

  if (!ref->InheritsFrom(AliMonitorHisto::Class())) return;
  ((AliMonitorHisto*)ref)->GetRun();
  SetReference(((AliMonitorHisto*)ref)->fHistoDraw);
}


//_____________________________________________________________________________
void AliMonitorHisto::Fill(Axis_t x)
{
// fill the monitor histogram

  fHisto->Fill(x);
}

//_____________________________________________________________________________
void AliMonitorHisto::Fill(Axis_t x, Axis_t y)
{
// fill the monitor histogram

  fHisto->Fill(x, y);
}

//_____________________________________________________________________________
void AliMonitorHisto::Fill(Axis_t x, Axis_t y, Stat_t w)
{
// fill the monitor histogram

  if (fHisto->InheritsFrom(TH2::Class())) {
    ((TH2*)fHisto)->Fill(x, y, w);
  } else if (fHisto->InheritsFrom(TProfile::Class())) {
    ((TProfile*)fHisto)->Fill(x, y, w);
  } else {
    AliError("trying to fill x and y of a 1 dimensinal histogram");
    return;
  }
}

//_____________________________________________________________________________
void AliMonitorHisto::ScaleErrorBy(Double_t factor)
{
// multiply the error of each bin by the given factor

  Int_t yMax = 1; 
  if (fHisto->GetDimension() > 1) 
    yMax = fHisto->GetYaxis()->GetNbins();
  Int_t xMax = fHisto->GetXaxis()->GetNbins();
  for (Int_t iY = 1; iY <= yMax; iY++) {
    for (Int_t iX = 1; iX <= xMax; iX++) {
      Int_t iBin = fHisto->GetBin(iX, iY);
      fHisto->SetBinError(iBin, factor * fHisto->GetBinError(iBin));
    }
  }
}


//_____________________________________________________________________________
void AliMonitorHisto::Update()
{
// update the normalized data histogram

  Update(fHisto);
  fHisto->Reset();
}

//_____________________________________________________________________________
void AliMonitorHisto::Update(TH1* histo)
{
// update the normalized data histogram using the given histo instead of fHisto

  fNumberOfEvents++;
  while (fNHistos >= fgNHistosMax) {
    fHistoList.Remove(fHistoList.LastLink());
    fNHistos--;
  }
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoList.AddFirst(histo->Clone());
  TH1::AddDirectory(addStatus);
  fNHistos++;
  fHistoRun->Add(histo);
  fHisto->Reset();
}

//_____________________________________________________________________________
void AliMonitorHisto::Add(AliMonitorPlot* plot)
{
// merge the given histo to this one

  if (!plot->InheritsFrom(AliMonitorHisto::Class())) return;
  AliMonitorHisto* histo = (AliMonitorHisto*) plot;

  fNumberOfEvents += histo->fNumberOfEvents;
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TObjLink* link = histo->fHistoList.LastLink();
  while (link) {
    fHistoList.AddFirst(link->GetObject()->Clone());
    link = link->Prev();
  }
  TH1::AddDirectory(addStatus);
  fNHistos += histo->fNHistos;
  while (fNHistos > fgNHistosMax) {
    fHistoList.Remove(fHistoList.LastLink());
    fNHistos--;
  }
  fHistoRun->Add(histo->fHistoRun);
}

//_____________________________________________________________________________
void AliMonitorHisto::Reset()
{
// reset the monitor histogram for a new run

  fHisto->Reset();
  fHistoList.Delete();
  fNHistos = 0;
  fHistoRun->Reset();
  if (fHistoDraw) delete fHistoDraw;
  fHistoDraw = NULL;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;
  fNumberOfEvents = 0;
}

//_____________________________________________________________________________
void AliMonitorHisto::ResetList()
{
// reset the the list of monitor histograms

  fHistoList.Delete();
  fNHistos = 0;
}


//_____________________________________________________________________________
void AliMonitorHisto::Scale(Int_t nEvents)
{
// scale the histogram to the correct normalization

  Double_t scale = 1.;
  switch (fNorm) {
  case kNormNone    : scale = 1.; break;
  case kNormEvents  : scale = 1./nEvents; break;
  case kNormEntries : scale = ((fHistoDraw->GetEntries() > 0) ? 
			       1./fHistoDraw->GetEntries() : 1.); break;
  case kNormIntegral: scale = ((fHistoDraw->Integral() > 0) ? 
			       1./fHistoDraw->Integral() : 1.); break;
  }
  fHistoDraw->Scale(scale);
}

//_____________________________________________________________________________
Bool_t AliMonitorHisto::ComparePlot()
{
// compare the data histogram to the reference histogram
// if they deviate by more than fgThreshold standard deviations in a bin,
// this bin is set in fHistoCompare and kFALSE is returned
  
  if (!fHistoRef) return kTRUE;
  if (fgThreshold <= 0) return kTRUE;
  if (!fHistoDraw) {
    AliWarning("no data histogram available for comparison\ncall DrawEvent, DrawSum or DrawRaw before calling Compare");
    return kTRUE;
  }
  if (fHistoCompare) delete fHistoCompare;
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoCompare = (TH1*) fHistoDraw->Clone();
  TH1::AddDirectory(addStatus);
  fHistoCompare->Reset();
  Bool_t result = kTRUE;

  Int_t yMax = 1; 
  if (fHistoDraw->GetDimension() > 1) 
    yMax = fHistoDraw->GetYaxis()->GetNbins();
  Int_t xMax = fHistoDraw->GetXaxis()->GetNbins();
  for (Int_t iY = 1; iY <= yMax; iY++) {
    for (Int_t iX = 1; iX <= xMax; iX++) {
      Int_t iBin = fHistoDraw->GetBin(iX, iY);
      Double_t delta = TMath::Abs(fHistoDraw->GetBinContent(iBin) -
				  fHistoRef->GetBinContent(iBin));
      Double_t errorData = fHistoDraw->GetBinError(iBin);
      Double_t errorRef = fHistoRef->GetBinError(iBin);
      Double_t sigma = TMath::Sqrt(errorData*errorData + errorRef*errorRef);
      if (delta > fgThreshold*sigma) {
	fHistoCompare->SetBinContent(iBin, fHistoDraw->GetBinContent(iBin));
	fHistoCompare->SetBinError(iBin, errorData);
	result = kFALSE;
      }
    }
  }

  return result;
}

//_____________________________________________________________________________
Bool_t AliMonitorHisto::GetEvent(Int_t number)
{
// get the normalized monitor histogram for the "number"th last event

  if (fNHistos == 0) {
    AliWarning("there are no histograms for single events available");
    return kFALSE;
  }
  if (number > fNHistos) {
    AliError(Form("requested event number (%d) exceeds range of available events (%d)", 
		  number, fNHistos));
    return kFALSE;
  }
  if (number <= 0) return kFALSE;

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;

  TObjLink* link = fHistoList.FirstLink();
  for (Int_t i = 1; i < number; i++) link = link->Next();
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoDraw = (TH1*) link->GetObject()->Clone();
  TH1::AddDirectory(addStatus);

  Scale(1);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorHisto::GetSum(Int_t number)
{
// get the normalized monitor histogram for the sum of the last 
// "number" events

  if (fNHistos == 0) {
    AliWarning("there are no histograms for single events available");
    return kFALSE;
  }
  if (number > fNHistos) {
    AliError(Form("requested number of events (%d) exceeds range of available events (%d)\nusing last %d event(s)", 
		  number, fNHistos, fNHistos));
    number = fNHistos;
  }
  if (number <= 0) return kFALSE;

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;

  TObjLink* link = fHistoList.FirstLink();
  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoDraw = (TH1*) link->GetObject()->Clone();
  TH1::AddDirectory(addStatus);
  for (Int_t i = 1; i < number; i++) {
    link = link->Next();
    fHistoDraw->Add((TH1*) link->GetObject());
  }

  Scale(number);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorHisto::GetRun()
{
// get the normalized monitor histogram for all monitored events 
// of the current run

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;

  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistoDraw = (TH1*) fHistoRun->Clone();
  TH1::AddDirectory(addStatus);

  Scale(fNumberOfEvents);
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorHisto::DrawPlot()
{
// draw the histograms

  fHistoDraw->SetMarkerColor(fgColorData);
  fHistoDraw->SetLineColor(fgColorData);
  fHistoDraw->SetFillColor(fgColorData);
  fHistoDraw->DrawCopy();

  if (fHistoRef && fgDrawRef) {
    char option[256];
    sprintf(option, "%sSAME", fHistoDraw->GetOption());

    if (fHistoRef->GetMaximum() > fHistoDraw->GetMaximum()) {
      fHistoDraw->SetMaximum(fHistoRef->GetMaximum() * 1.1);
    }

    fHistoRef->SetMarkerColor(fgColorRef);
    fHistoRef->SetLineColor(fgColorRef);
    fHistoRef->SetFillColor(fgColorRef);
    fHistoRef->DrawCopy(option);

    fHistoDraw->DrawCopy(option);

    if (fHistoCompare && (fgThreshold > 0)) {
      fHistoCompare->SetMarkerColor(fgColorCompare);
      fHistoCompare->SetLineColor(fgColorCompare);
      fHistoCompare->SetFillColor(fgColorCompare);
      fHistoCompare->DrawCopy(option);
    }
  }

  gPad->Update();
}
