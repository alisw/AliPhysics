#ifndef UtilHisto_h
#define UtilHisto_h

#include <iostream>
#include "UtilMessages.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

//_________________________________________________________________________________________________
Bool_t SameBinning(const TH1 *h1, const TH1 *h2, const Bool_t check = kFALSE){
  Infomsg("SameBinning", "Checking same binning");
  const Int_t n[2] = {h1->GetNbinsX(), h2->GetNbinsX()};
  if(n[0] != n[1]){
    Errormsg("SameBinning", Form("Different number of bins"));
    return kFALSE;
  } 
  
  if(check){
    for(Int_t i = 1; i <= n[0]; i++){
      const Double_t x[4] = {
        h1->GetXaxis()->GetBinLowEdge(i),
        h2->GetXaxis()->GetBinLowEdge(i),
        h1->GetXaxis()->GetBinUpEdge(i),
        h2->GetXaxis()->GetBinUpEdge(i)
      };
      for(Int_t j = 0; j < 2; j++){
        if(x[0+2*j] != x[1+2*j]){
          Errormsg("SameBinning", Form("Binning is different %i [%f, %f]", i, x[0+2*j], x[1+2*j]));
          return kFALSE;
        }
      }
    }
  }
  return kTRUE;
  
}


//_________________________________________________________________________________________________
Int_t GetHistoLowRange(TH1F *h, const Int_t nmin, const Int_t ncounts = 3, const Int_t secondmin = 1){
  if(h->GetEffectiveEntries() < nmin) return 1;
  const Int_t firstbin = h->FindFirstBinAbove(nmin);
  Int_t combo = 0;
  for(Int_t i = firstbin; i < h->GetNbinsX(); i++){
    if(h->GetBinContent(i) > secondmin) combo++;
    else combo = 0;
    if(combo > ncounts) return i;
  }
  return 1;
}

//_________________________________________________________________________________________________
Int_t GetFirstHistogram(const TObjArray* templates, const Int_t nmin = 10){
  const Int_t entries = templates->GetEntries();
  Int_t binindex[entries];
  Int_t sorted[entries];
  for(Int_t i = 0; i < entries; i++){
    binindex[i] = GetHistoLowRange(static_cast<TH1F*>(templates->At(i)), nmin);
  }
  TMath::Sort(entries, binindex, sorted, kFALSE);//First shall be the one with the lowest value!
  
  return sorted[0];
}

//_________________________________________________________________________________________________
Int_t GetHistoLowRangeAfter(TH1F *h, const Int_t nmin, const Int_t after = 1){
  if(h->GetEffectiveEntries() < nmin){
    Errormsg("GetHistoLowRangeAfter", "Histogram has not enough entries");
    return 1;
  }
  for(Int_t i = after +1; i <= h->GetNbinsX(); i++) if(h->GetBinContent(i) > nmin) return i;
  Warningmsg("GetHistoLowRangeAfter", "Cannot find more ranges");
  return after;
  
}

//_________________________________________________________________________________________________
Int_t GetHistoNoHolesAfter(TH1F *h, const Int_t after = 1, const Int_t check = 20){
  Infomsg("GetHistoNoHolesAfter", Form("Checking that there are no holes in %s", h->GetName()));
  for(Int_t i = after + check; i > after; i--){
    //     cout<<i<<" "<<h->GetBinContent(i)<<endl;
    if(h->GetBinContent(i) < 1) return i+1;
  }
  
  return after;
  
}

//_________________________________________________________________________________________________
Int_t GetHistoNoHolesBefore(TH1F *h, const Double_t before = 0){
  Infomsg("GetHistoNoHolesBefore", Form("Checking that there are no holes in %s", h->GetName()));
  const Int_t start = h->FindFirstBinAbove(0);
  const Int_t stop = h->FindBin(before);
  for(Int_t i = stop; i > start; i--){
    if(h->GetBinContent(i) < 1) return i+1;
  }
  return start;
}

//_________________________________________________________________________________________________
Int_t GetHistoNoHolesBefore(TObjArray *h, const Double_t before = 0){
  Infomsg("GetHistoNoHolesBefore", "Checking that there are no holes in array of histograms");
  Int_t limit = 1;
  for(Int_t i = 0; i < h->GetEntries(); i++){
    const TString objclass = h->At(i)->ClassName();
    if(!objclass.Contains("TH1")) continue;
    Int_t l = GetHistoNoHolesBefore(static_cast<TH1F*>(h->At(i)));
    if(l > limit) limit = l;
  }
  return limit;
}

//_________________________________________________________________________________________________
Bool_t IsHistogramInRange(TH1 *h, const Double_t rangelow, const Double_t rangehigh, const Double_t threshold = 1., const Bool_t verbose = kFALSE){
  const Int_t binlow = h->FindBin(rangelow);
  const Int_t binhigh = h->FindBin(rangehigh);
  if(verbose) cout<<"Search for content in interval ["<<binlow<<" ("<<rangelow<<"), "<<binhigh<<" ("<<rangehigh<<")]"<<endl;
  Double_t integral = 0;
  for(Int_t i = binlow; i <= binhigh; i++){
    const Double_t value = h->GetBinContent(i);
    if(value <= 0) continue;
    if(verbose) cout<<"Bin "<<i<<"("<<h->GetXaxis()->GetBinCenter(i)<<") has content "<<value<<endl;
    integral += value;
  }
  if(integral > threshold) return kTRUE;
  else return kFALSE;
}

//_________________________________________________________________________________________________
Double_t GetResidualYield(TH1 *hdata, TH1 *hfit, const Double_t rangelow, const Double_t rangehigh){
  if(!SameBinning(hdata, hfit)) Fatalmsg("GetResidualYield", "Different binnings");
  const Int_t binlow = hdata->FindBin(rangelow);
  const Int_t binhigh = hdata->FindBin(rangehigh);
  Double_t integraldata = 0;
  Double_t integralfit = 0;
  for(Int_t i = binlow; i <= binhigh; i++){
    const Double_t valuedata = hdata->GetBinContent(i);
    const Double_t valuefit = hfit->GetBinContent(i);
    if(valuedata > 0) integraldata += valuedata;
    if(valuefit > 0) integralfit += valuefit;
  }
  const Double_t diff = integraldata - integralfit;
  Infomsg("GetResidualYield", Form("Residual Yield is %f", diff));
  // if(diff > 0) return diff;
  // else 
  return diff;
}

//_________________________________________________________________________________________________
Double_t GetOverlapFraction(const TH1F *h, const TH1F *bkg, Double_t &error, const Bool_t show = kFALSE){
  Int_t firstbin = h->FindFirstBinAbove(2);
  Int_t lastbin = h->FindLastBinAbove(2);
  if(firstbin == lastbin){
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
  
  if(show){
    TVirtualPad *pad = gPad;
    TCanvas *cCanvas = new TCanvas("cOverlapFraction", "Overlap Fraction");
    cCanvas->cd();
    h->DrawCopy()->SetName("signal");
    bkg->DrawCopy("same")->SetName("background");
    TH1F *ranges = (TH1F*) h->Clone("overlapranges");
    ranges->Reset();
    ranges->SetLineColor(kCyan);
    ranges->SetMarkerColor(kCyan);
    ranges->Fill(first, h->GetMaximum());
    ranges->Fill(last, h->GetMaximum());
    ranges->DrawCopy("same");
    pad->cd();
  }
  
  const Double_t sum = yield + bkgyield;
  const Double_t sumerror = TMath::Sqrt(yielderror*yielderror + bkgerror*bkgerror);
  const Double_t fraction = sum > 0 ? bkgyield/(sum) : -10;
  const Double_t fractionerror = TMath::Abs(yield*bkgyield) > 0 ? fraction*TMath::Sqrt(TMath::Power(bkgerror/bkgyield, 2) + TMath::Power(sumerror/sum, 2)) : 0;
  error = sum > 0 ? fractionerror : 0;
  
  if(yield + bkgyield > 0) return bkgyield/(yield + bkgyield);
  else return -10;
}

#endif
