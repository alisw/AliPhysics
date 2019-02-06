#include "UtilHisto.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "UtilMessages.h"
#include <iostream>

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Utilities for histogram management et similia                         //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

//_________________________________________________________________________________________________
//Function to get the starting point of the histogram with a minimum request on the count number
Bool_t SameBinning(const TH1* h1, const TH1* h2, const Bool_t check)
{
  Infomsg("SameBinning", "Checking same binning");
  if (!h1 || !h2) {
    Errormsg("SameBinning", "Input histogram is not found!");
    return 1;
  }
  const Int_t n[2] = { h1->GetNbinsX(), h2->GetNbinsX() };
  if (n[0] != n[1]) {
    Errormsg("SameBinning", Form("Different number of bins: %s %i and %s %i", h1->GetName(), n[0], h2->GetName(), n[1]));
    return kFALSE;
  }

  if (check) {
    for (Int_t i = 1; i <= n[0]; i++) {
      const Double_t x[4] = {
        h1->GetXaxis()->GetBinLowEdge(i),
        h2->GetXaxis()->GetBinLowEdge(i),
        h1->GetXaxis()->GetBinUpEdge(i),
        h2->GetXaxis()->GetBinUpEdge(i)
      };
      for (Int_t j = 0; j < 2; j++) {
        if (x[0 + 2 * j] != x[1 + 2 * j]) {
          Errormsg("SameBinning", Form("Binning is different %i [%f, %f]", i, x[0 + 2 * j], x[1 + 2 * j]));
          return kFALSE;
        }
      }
    }
  }
  return kTRUE;
}

//_________________________________________________________________________________________________
//Function to get the starting point of the histogram with a minimum request on the count number
Int_t GetHistoLowRange(TH1F* h, const Int_t nmin, const Int_t ncounts, const Int_t secondmin)
{
  if (!h || h->GetEffectiveEntries() < nmin)
    return 1;
  const Int_t firstbin = h->FindFirstBinAbove(nmin);
  Int_t combo = 0;
  for (Int_t i = firstbin; i < h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) > secondmin)
      combo++;
    else
      combo = 0;
    if (combo > ncounts)
      return i;
  }
  return 1;
}

//_________________________________________________________________________________________________
//Function to get the first histogram of a series
Int_t GetFirstHistogram(const TObjArray* templates, const Int_t nmin)
{
  const Int_t entries = templates->GetEntries();
  Int_t binindex[entries];
  Int_t sorted[entries];
  Int_t N = 0;
  for (Int_t i = 0; i < entries; i++) {
    TH1F *H = dynamic_cast<TH1F*>(templates->At(i));
    if(!H) {
      Infomsg("GetFirstHistogram", Form("While looking for the first Histogram I skipped the entry %i as it is not a TH1F: %s is a %s", i, templates->At(i)->GetName(), templates->At(i)->ClassName()));
      continue;
    }
    binindex[i] = GetHistoLowRange(H, nmin);
    N++;
  }
  TMath::Sort(N, binindex, sorted, kFALSE); //First shall be the one with the lowest value!
  return sorted[0];
}

//_________________________________________________________________________________________________
//Function to get the starting point after a certain value
Int_t GetHistoLowRangeAfter(TH1F* h, const Int_t nmin, const Int_t after)
{
  if (h == 0x0) {
    Errormsg("GetHistoLowRangeAfter", "Input histogram is not found!");
    return 1;
  }
  if (h->GetEffectiveEntries() < nmin) {
    Errormsg("GetHistoLowRangeAfter", "Histogram has not enough entries");
    return 1;
  }
  for (Int_t i = after + 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > nmin)
      return i;
  Warningmsg("GetHistoLowRangeAfter", "Cannot find more ranges");
  return after;
}

//_________________________________________________________________________________________________
//Function to get the starting point after a certain value, without any holes in the histogram
Int_t GetHistoNoHolesAfter(TH1F* h, const Int_t after, const Int_t check)
{
  if (!h) {
    Errormsg("GetHistoNoHolesAfter", "Input histogram is not found!");
    return 1;
  }
  Infomsg("GetHistoNoHolesAfter", Form("Checking that there are no holes in %s", h->GetName()));
  for (Int_t i = after + check; i > after; i--) {
    //     cout<<i<<" "<<h->GetBinContent(i)<<endl;
    if (h->GetBinContent(i) < 1)
      return i + 1;
  }

  return after;
}

//_________________________________________________________________________________________________
//Function to get the starting point before a certain value, without any holes in the histogram
Int_t GetHistoNoHolesBefore(TH1F* h, const Double_t before)
{
  Infomsg("GetHistoNoHolesBefore", Form("Checking that there are no holes in %s", h->GetName()));
  const Int_t start = h->FindFirstBinAbove(0);
  const Int_t stop = h->FindBin(before);
  for (Int_t i = stop; i > start; i--) {
    if (h->GetBinContent(i) < 1)
      return i + 1;
  }
  return start;
}

//_________________________________________________________________________________________________
//Function to get the last histogram of a series
Int_t GetHistoNoHolesBefore(TObjArray* h, const Double_t before)
{
  Infomsg("GetHistoNoHolesBefore", "Checking that there are no holes in array of histograms");
  Int_t limit = 1;
  for (Int_t i = 0; i < h->GetEntries(); i++) {
    const TString objclass = h->At(i)->ClassName();
    if (!objclass.Contains("TH1"))
      continue;
    Int_t l = GetHistoNoHolesBefore(static_cast<TH1F*>(h->At(i)));
    if (l > limit)
      limit = l;
  }
  return limit;
}

//_________________________________________________________________________________________________
//Function to check if the histogram is in range
Bool_t IsHistogramInRange(TH1* h, const Double_t rangelow, const Double_t rangehigh, const Double_t threshold, const Bool_t verbose)
{
  const Int_t binlow = h->FindBin(rangelow);
  const Int_t binhigh = h->FindBin(rangehigh);
  if (verbose)
    Printf("For %s searching for content in interval [%i (%f), %i (%f)] with threshold %f", h->GetName(), binlow, rangelow, binhigh, rangehigh, threshold);
  //
  Double_t integral = 0;
  for (Int_t i = binlow; i <= binhigh; i++) {
    const Double_t value = h->GetBinContent(i);
    if (verbose)
      Printf("Bin %i [%f, %f, %f] has content %f", i, h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinCenter(i), h->GetXaxis()->GetBinUpEdge(i), value);
    if (value <= 0)
      continue;
    integral += value;
  }
  if (integral > threshold) {
    if (verbose)
      Printf("Intgral is above threshold: %f > %f", integral, threshold);
    return kTRUE;
  } else {
    if (verbose)
      Printf("Intgral is below threshold: %f < %f", integral, threshold);
    return kFALSE;
  }
}

//_________________________________________________________________________________________________
//Function to get the residual yield from a fit, within a certain range
Double_t GetResidualYield(TH1* hdata, TH1* hfit, const Double_t rangelow, const Double_t rangehigh)
{
  if (!SameBinning(hdata, hfit))
    Fatalmsg("GetResidualYield", "Different binnings");
  const Int_t binlow = hdata->FindBin(rangelow);
  const Int_t binhigh = hdata->FindBin(rangehigh);
  Double_t integraldata = 0;
  Double_t integralfit = 0;
  for (Int_t i = binlow; i <= binhigh; i++) {
    const Double_t valuedata = hdata->GetBinContent(i);
    const Double_t valuefit = hfit->GetBinContent(i);
    if (valuedata > 0)
      integraldata += valuedata;
    if (valuefit > 0)
      integralfit += valuefit;
  }
  const Double_t diff = integraldata - integralfit;
  Infomsg("GetResidualYield", Form("Residual Yield is %f", diff));
  // if(diff > 0) return diff;
  // else
  return diff;
}

//_________________________________________________________________________________________________
//Function to get the overlap fraction of two histograms
Double_t GetOverlapFraction(const TH1F* h, const TH1F* bkg, Double_t& error, const Bool_t show)
{
  Int_t firstbin = h->FindFirstBinAbove(2);
  Int_t lastbin = h->FindLastBinAbove(2);
  if (firstbin == lastbin) {
    firstbin = 1;
    lastbin = h->GetNbinsX();
  }
  const Double_t first = h->GetXaxis()->GetBinCenter(firstbin);
  const Double_t last = h->GetXaxis()->GetBinCenter(lastbin);

  const Int_t firstbinbkg = h->GetXaxis()->FindBin(first);
  const Int_t lastbinbkg = h->GetXaxis()->FindBin(last);
  Double_t yielderror = 0;
  const Double_t yield = h->IntegralAndError(firstbin, lastbin, yielderror);
  Double_t bkgerror = 0;
  const Double_t bkgyield = bkg->IntegralAndError(firstbinbkg, lastbinbkg, bkgerror);

  Infomsg("GetOverlapFraction", Form("Integrating signal %s with %.0f entries [%f (%i), %f (%i)] \nand background %s with %.0f entries [%f (%i),%f (%i)]", h->GetName(), h->GetEntries(), first, firstbin, last, lastbin, bkg->GetName(), bkg->GetEntries(), first, firstbinbkg, last, lastbinbkg));

  if (show) {
    TVirtualPad* pad = gPad;
    TCanvas* cCanvas = new TCanvas("cOverlapFraction", "Overlap Fraction");
    cCanvas->cd();
    h->DrawCopy()->SetName("signal");
    bkg->DrawCopy("same")->SetName("background");
    TH1F* ranges = (TH1F*)h->Clone("overlapranges");
    ranges->Reset();
    ranges->SetLineColor(kCyan);
    ranges->SetMarkerColor(kCyan);
    ranges->Fill(first, h->GetMaximum());
    ranges->Fill(last, h->GetMaximum());
    ranges->DrawCopy("same");
    pad->cd();
  }

  const Double_t sum = yield + bkgyield;
  const Double_t sumerror = TMath::Sqrt(yielderror * yielderror + bkgerror * bkgerror);
  const Double_t fraction = sum > 0 ? bkgyield / (sum) : -10;
  const Double_t fractionerror = TMath::Abs(yield * bkgyield) > 0 ? fraction * TMath::Sqrt(TMath::Power(bkgerror / bkgyield, 2) + TMath::Power(sumerror / sum, 2)) : 0;
  error = sum > 0 ? fractionerror : 0;

  if (yield + bkgyield > 0)
    return bkgyield / (yield + bkgyield);
  else
    return -10;
}

//_________________________________________________________________________________________________
//Function to clean a histogram, i.e. remove the holes
TH1* Clean(TH1* h, const Float_t center, const Bool_t add, const Int_t nchain, const Float_t thr)
{
  TH1* hdiff = (TH1*)h->Clone(Form("clean%s", h->GetName()));
  hdiff->SetMarkerStyle(4);
  hdiff->SetLineColor(kGray);
  hdiff->Reset();
  Int_t bin0 = h->GetXaxis()->FindBin(center);
  Int_t nconsec = 0;
  for (Int_t i = bin0; i < h->GetNbinsX(); i++) {
    if (nconsec == nchain) {
      for (Int_t j = 2 + i - nchain; j < h->GetNbinsX(); j++)
        hdiff->SetBinContent(j, h->GetBinContent(j));
      //
      break;
    }
    if (h->GetBinContent(i) <= thr)
      nconsec++;
    else
      nconsec = 0;
  }
  //
  nconsec = 0;
  for (Int_t i = bin0; i >= 1; i--) {
    if (nconsec == nchain) {
      for (Int_t j = i + nchain - 2; j >= 1; j--)
        hdiff->SetBinContent(j, h->GetBinContent(j));
      //
      break;
    }
    if (h->GetBinContent(i) <= thr)
      nconsec++;
    else
      nconsec = 0;
  }
  if (add)
    h->Add(hdiff, -1);
  return hdiff;
}
