#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TArrayI.h>
#include <TRandom.h>
#include <TLegend.h>
#include <iostream>
#include <TH1.h>

void
convert(UShort_t rate, Int_t adc, Int_t last, TArrayI& counts)
{
  if (rate == 1) { 
    Float_t a = adc;
    if (a < 0) a = 0;
    counts[0] = UShort_t(TMath::Min(a, 1023.F));
    return;
  }
  Float_t b = 6;
  for (Ssiz_t i = 0; i < rate; i++) { 
    Float_t t = Float_t(i) / rate + 1./rate;
    Float_t s = adc + (last - adc) * TMath::Exp(-b * t);
    Float_t a = s;
    if (a < 0) a = 0;
    counts[i] = UShort_t(TMath::Min(a, 1023.F));
    std::cout << " rate=" << rate << "\tadc=" << adc 
	      << "\tcount[" << i << "]=" << counts[i] 
	      << "\ts=" << adc << " + (" << last << " - " 
	      << adc << ") * exp(-" << b << " * " << t 
	      << ")=" << s << std::endl;
  }
  return;
}

TGraph* 
makeGraph(const TArrayI& adcs, Int_t rate) 
{
  Int_t    last = adcs.fArray[0];
  TArrayI  counts(4);
  TGraph*  graph = new TGraph(rate * adcs.fN);
  graph->SetLineColor(rate);
  graph->SetMarkerColor(rate);
  graph->SetMarkerStyle(20+rate);
  graph->SetLineStyle(rate);
  graph->SetName(Form("rate%d", rate));
  graph->SetTitle(Form("Rate %d", rate));
  for (Int_t i = 0; i < adcs.fN; i++) { 
    counts.Reset(-1);
    convert(rate, adcs.fArray[i], last, counts);
    
    for (Int_t j = 0; j < rate; j++) { 
      Int_t    idx = (i * rate + j);
      Double_t x   = (i + (rate > 1 ? Float_t(j+1) / rate-1 : 0));
      graph->SetPoint(idx, x, counts[j]);
    }
    last = counts[rate - 1];
  }
  return graph;
}

void
TestShaping(int max=4)
{
  TArrayI adcs(10);
  TGraph* orig = new TGraph(adcs.fN);
  orig->SetName("Original");
  orig->SetTitle("Original");
  orig->SetMarkerStyle(25);
  orig->SetMarkerColor(1);
  orig->SetMarkerSize(2);
  orig->SetLineColor(1);
  for (Int_t i = 0; i < adcs.fN; i++) { 
    adcs.fArray[i] = Int_t(gRandom->Uniform(0, 1023));
    orig->SetPoint(i, i, adcs.fArray[i]);
  }

  TCanvas* c = new TCanvas("c", "c");
  c->SetFillColor(0);
  c->SetTopMargin(.02);
  c->SetRightMargin(.02);
  
  TH1* h = new TH1F("frame","frame", adcs.fN+1, -2, adcs.fN);
  h->SetMinimum(0);
  h->SetMaximum(1300);
  h->SetStats(0);
  h->Draw("");
  orig->Draw("pl same");

  TLegend* l = new TLegend(adcs.fN*3./4, 1023, adcs.fN, 1300, "", "");
  l->SetFillColor(0);
  l->SetBorderSize(1);
  l->AddEntry(orig, orig->GetTitle(), "lp");
  
  for (int i = 1; i <= max; i++) {
    TGraph* g = makeGraph(adcs, i);
    g->Draw("pl same");
    l->AddEntry(g, g->GetTitle(), "lp");
  }
  l->Draw();
  
  c->Modified();
  c->Update();
  c->cd();
}

  
  
    

     
