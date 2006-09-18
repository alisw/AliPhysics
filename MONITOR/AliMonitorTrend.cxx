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
//  This class maintains a plot for the evolution of a value that is used    //
//  tomonitor the quality of the recorded data.                              //
//  The trendgram is created and filled by a sub class of AliMonitor.        //
//  It can be compared to a reference trendgram. For the comparison a        //
//  maximal deviation (in standard deviations) can be specified.             //
//  The bins where the maximal deviation is exceeded are drawn in red.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TH1.h>
#include <TVirtualPad.h>
#include <TLine.h>

#include "AliLog.h"

#include "AliMonitorTrend.h"


ClassImp(AliMonitorTrend) 


Int_t AliMonitorTrend::fgIncSize = 10;


//_____________________________________________________________________________
AliMonitorTrend::AliMonitorTrend():
  AliMonitorPlot(),
  fLabel(),
  fMin(0),
  fMax(0),
  fData(),
  fHistoDraw(NULL),
  fRefMean(0),
  fRefSigma(-1),
  fHistoCompare(NULL)
{
// default contructor

}

//_____________________________________________________________________________
AliMonitorTrend::AliMonitorTrend(const AliMonitorTrend& trend) :
  AliMonitorPlot(trend),
  fLabel(trend.fLabel),
  fMin(trend.fMin),
  fMax(trend.fMax),
  fData(trend.fData),
  fHistoDraw(NULL),
  fRefMean(trend.fRefMean),
  fRefSigma(trend.fRefSigma),
  fHistoCompare(NULL)
{
// copy constructor

}

//_____________________________________________________________________________
AliMonitorTrend& AliMonitorTrend::operator =(const AliMonitorTrend& trend)
{
// assignment operator

  AliMonitorPlot::operator =(trend);

  fLabel = trend.fLabel;
  fMin = trend.fMin;
  fMax = trend.fMax;
  trend.fData.Copy(fData);

  fHistoDraw = NULL;
  fRefMean = trend.fRefMean;
  fRefSigma = trend.fRefSigma;
  fHistoCompare = NULL;

  return *this;
}

//_____________________________________________________________________________
AliMonitorTrend::AliMonitorTrend(const char* name, const char* title,
		  const char* label, Double_t min, Double_t max) :
  AliMonitorPlot(name, title),
  fLabel(label),
  fMin(min),
  fMax(max),
  fData(),
  fHistoDraw(NULL),
  fRefMean(0),
  fRefSigma(0),
  fHistoCompare(NULL)
{
// create a monitor trend

}

//_____________________________________________________________________________
AliMonitorTrend::~AliMonitorTrend()
{
// delete all histograms

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
}


//_____________________________________________________________________________
void AliMonitorTrend::SetReference(TH1* ref)
{
// set the reference trend for comparison

  Int_t n = ref->GetXaxis()->GetNbins();
  if (n <= 0) return;

  Double_t sum = 0;
  Double_t sum2 = 0;
  for (Int_t i = 1; i <= n; i++) {
    sum += ref->GetBinContent(i);
    sum2 += ref->GetBinContent(i) * ref->GetBinContent(i);
  }

  fRefMean = sum / n;
  fRefSigma = TMath::Sqrt(sum2 - sum*sum/n) / n;
}

//_____________________________________________________________________________
void AliMonitorTrend::SetReference(AliMonitorPlot* ref)
{
// set the reference trendgram for comparison

  if (!ref->InheritsFrom(AliMonitorTrend::Class())) return;
  fRefMean = ((AliMonitorTrend*)ref)->GetMean();
  fRefSigma = ((AliMonitorTrend*)ref)->GetSigma();
}


//_____________________________________________________________________________
void AliMonitorTrend::Fill(Double_t x)
{
// add a value to the monitor trend

  if (fNumberOfEvents >= fData.GetSize()) {
    fData.Set(fNumberOfEvents + fgIncSize);
  }
  fData[fNumberOfEvents] = x;
}


//_____________________________________________________________________________
void AliMonitorTrend::Update()
{
// update

  fNumberOfEvents++;
}

//_____________________________________________________________________________
void AliMonitorTrend::Add(AliMonitorPlot* plot)
{
// merge the given trend to this one

  if (!plot->InheritsFrom(AliMonitorTrend::Class())) return;
  AliMonitorTrend* trend = (AliMonitorTrend*) plot;

  Int_t numberOfEvents = fNumberOfEvents + trend->fNumberOfEvents;
  if (numberOfEvents >=  fData.GetSize()) {
    fData.Set(numberOfEvents + fgIncSize);
  }
  for (Int_t i = 0; i < trend->fNumberOfEvents; i++) {
    fData[fNumberOfEvents + i] = trend->fData[i];
  }
  fNumberOfEvents = numberOfEvents;
}

//_____________________________________________________________________________
void AliMonitorTrend::Reset()
{
// reset the monitor trend for a new run

  fData.Set(fgIncSize);
  if (fHistoDraw) delete fHistoDraw;
  fHistoDraw = NULL;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;
  fNumberOfEvents = 0;
}

//_____________________________________________________________________________
void AliMonitorTrend::ResetList()
{
// reset the the list of monitor histograms
// (not applicable for trend)

}


//_____________________________________________________________________________
Bool_t AliMonitorTrend::ComparePlot()
{
// compare the data trend to the reference
// if they deviate by more than fgThreshold standard deviations in a bin,
// this bin is set in fHistoCompare and kFALSE is returned
  
  if (fRefSigma < 0) return kTRUE;
  if (fgThreshold <= 0) return kTRUE;
  if (!fHistoDraw) {
    AliWarning("no data trend available for comparison\ncall DrawSum or DrawRaw before calling Compare");
    return kTRUE;
  }

  Int_t nBins = fHistoDraw->GetXaxis()->GetNbins();
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = CreateHisto(nBins);
  fHistoCompare->Reset();
  fHistoCompare->SetOption("P");
  fHistoCompare->SetMarkerStyle(kFullCircle);
  fHistoCompare->SetFillStyle(0);
  Bool_t result = kTRUE;

  for (Int_t iBin = 1; iBin <= nBins; iBin++) {
    Double_t delta = TMath::Abs(fHistoDraw->GetBinContent(iBin) - fRefMean);
    if (delta > fgThreshold*fRefSigma) {
      fHistoCompare->SetBinContent(iBin, fHistoDraw->GetBinContent(iBin));
      result = kFALSE;
    }
  }

  return result;
}

//_____________________________________________________________________________
Bool_t AliMonitorTrend::GetEvent(Int_t)
{
// there is no single event trend

//  Info("GetEvent", "there is no trend for single events available");
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMonitorTrend::GetSum(Int_t number)
{
// get the monitor trend for the last  "number" events

  if (number > fNumberOfEvents) {
    AliError(Form("requested number of events (%d) exceeds range of available events (%d)\nusing last %d event(s)", 
	  number, fNumberOfEvents, fNumberOfEvents));
    number = fNumberOfEvents;
  }
  if (number <= 0) return kFALSE;

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;

  fHistoDraw = CreateHisto(number);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorTrend::GetRun()
{
// get the monitor trend for all monitored events of the current run

  if (fHistoDraw) delete fHistoDraw;
  if (fHistoCompare) delete fHistoCompare;
  fHistoCompare = NULL;
  if (fNumberOfEvents <= 0) return kFALSE;

  fHistoDraw = CreateHisto(fNumberOfEvents);
  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorTrend::DrawPlot()
{
// draw the trendgrams

  fHistoDraw->SetMarkerColor(fgColorData);
  fHistoDraw->SetLineColor(fgColorData);
  fHistoDraw->SetLineWidth(2);
  fHistoDraw->DrawCopy();

  if ((fRefSigma > 0) && fgDrawRef) {
    if ((fRefMean+fRefSigma > fHistoDraw->GetMaximum()) && !(fMax > fMin)) {
      fHistoDraw->SetMaximum(fRefMean+fRefSigma * 1.1);
    }

    Double_t xMin = fHistoDraw->GetXaxis()->GetXmin();
    Double_t xMax = fHistoDraw->GetXaxis()->GetXmax();
    TLine* mean = new TLine(xMin, fRefMean, xMax, fRefMean);
    mean->SetLineColor(fgColorRef);
    mean->SetLineWidth(2);
    mean->Draw();
    TLine* high = new TLine(xMin, fRefMean+fRefSigma, 
			    xMax, fRefMean+fRefSigma);
    high->SetLineColor(fgColorRef);
    high->SetLineWidth(2);
    high->SetLineStyle(2);
    high->Draw();
    TLine* low = new TLine(xMin, fRefMean-fRefSigma, 
			   xMax, fRefMean-fRefSigma);
    low->SetLineColor(fgColorRef);
    low->SetLineWidth(2);
    low->SetLineStyle(2);
    low->Draw();

//    char option[256];
//    sprintf(option, "%sSAME", fHistoDraw->GetOption());
//    fHistoDraw->DrawCopy(option);

    if (fHistoCompare && (fgThreshold > 0)) {
      char option[256];
      sprintf(option, "%sSAME", fHistoCompare->GetOption());
      fHistoCompare->SetMarkerColor(fgColorCompare);
      fHistoCompare->SetLineColor(fgColorCompare);
      fHistoCompare->SetLineWidth(2);
      fHistoCompare->DrawCopy(option);
    }
  }

  gPad->Update();
}


//_____________________________________________________________________________
TH1* AliMonitorTrend::CreateHisto(Int_t nBins)
{
// create a histogram for a trend plot with the last nBin entries

  TH1* result = new TH1D(GetName(), GetTitle(), nBins, -nBins-0.5, -0.5);
  result->GetXaxis()->SetTitle("N_{event}");
  result->GetYaxis()->SetTitle(fLabel.Data());
  if (fMax > fMin) {
    result->SetMinimum(fMin);
    result->SetMaximum(fMax);
  }
  result->SetOption("L");

  Double_t sum = 0;
  Double_t sum2 = 0;
  for (Int_t i = 0; i < nBins; i++) {
    Double_t data = fData[fNumberOfEvents-1-i];
    sum += data;
    sum2 += data * data;
    result->SetBinContent(nBins-i, data);
  }
  Stat_t stats[4];
  stats[0] = nBins;
  stats[1] = nBins * nBins;
  stats[2] = sum;
  stats[3] = sum2;
  result->PutStats(stats);

  return result;
}

//_____________________________________________________________________________
Double_t AliMonitorTrend::GetMean() const
{
// get the mean value

  if (fNumberOfEvents <= 0) return 0;

  Double_t sum = 0;
  for (Int_t i = 0; i < fNumberOfEvents; i++) {
    sum += fData[i];
  }
  return sum / fNumberOfEvents;
}

//_____________________________________________________________________________
Double_t AliMonitorTrend::GetSigma() const
{
// get the rms value

  if (fNumberOfEvents <= 0) return 0;

  Double_t sum = 0;
  Double_t sum2 = 0;
  for (Int_t i = 0; i < fNumberOfEvents; i++) {
    sum += fData[i];
    sum2 += fData[i] * fData[i];
  }
  return TMath::Sqrt(sum2 - sum*sum/fNumberOfEvents) / fNumberOfEvents;
}

