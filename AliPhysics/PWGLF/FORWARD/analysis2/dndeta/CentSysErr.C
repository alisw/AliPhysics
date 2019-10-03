#include <TFile.h>
#include <TCollection.h>
#include <TList.h>
#include <TClass.h>
#include <THStack.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TError.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TParameter.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TMath.h>

//____________________________________________________________________
TObject*
GetObject(TDirectory* d, const char* name, TClass* cls=0, Bool_t verb=true)
{
  if (!d) {
    Error("GetOject", "No directory");
    return 0;
  }
  
  TObject* o = d->Get(name);
  if (!o) {
    if (verb) 
      Error("GetObject", "No object %s in directory %s", name, d->GetName());
    return 0;
  }

  if (!cls) return o;

  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb) 
      Error("GetObject", "%s from %s is not a %s, but a %s",
	    o->GetName(), d->GetName(), cls->GetName(), o->ClassName());
    return 0;
  }

  return o;
}

//____________________________________________________________________
TObject*
GetObject(TCollection* d, const char* name, TClass* cls=0, Bool_t verb=true)
{
  if (!d) {
    Error("GetOject", "No collection");
    return 0;
  }
  
  TObject* o = d->FindObject(name);
  if (!o) {
    if (verb)
      Error("GetObject", "No object %s in collection %s", name, d->GetName());
    return 0;
  }

  if (!cls) return o;

  if (!o->IsA()->InheritsFrom(cls)) {
    if (verb)
      Error("GetObject", "%s from %s is not a %s, but a %s",
	    o->GetName(), d->GetName(), cls->GetName(), o->ClassName());
    return 0;
  }

  return o;
}
//____________________________________________________________________
TCollection*
GetCollection(TDirectory* d, const char* name)
{
  return static_cast<TCollection*>(GetObject(d, name, TCollection::Class()));
}
//____________________________________________________________________
TCollection*
GetCollection(TCollection* d, const char* name)
{
  return static_cast<TCollection*>(GetObject(d, name, TCollection::Class()));
}
//____________________________________________________________________
TH2*
GetH2(TCollection* c, const char* name)
{
  return static_cast<TH2*>(GetObject(c, name, TH2::Class()));
}
//____________________________________________________________________
TH1*
GetH1(TCollection* c, const char* name)
{
  return static_cast<TH1*>(GetObject(c, name, TH1::Class()));
}
//____________________________________________________________________
THStack*
GetStack(TCollection* c, const char* name, Bool_t verb=true)
{
  return static_cast<THStack*>(GetObject(c, name, THStack::Class(),verb));
}
//____________________________________________________________________
TParameter<double>*
GetParam(TCollection* c, const char* name, Bool_t verb=false)
{
  return static_cast<TParameter<double>*>(GetObject(c, name,
						    TParameter<double>::Class(),
						    verb));
}


//____________________________________________________________________
void
ProcessOne(const char* meth,
	   Int_t       style,
	   TList*      stacks,
	   TList*      ratios,
	   TList*      mins,
	   TList*      maxs,
	   TList*      avgs,
	   TAxis*&     centAxis)
{
  TFile* file = TFile::Open(Form("PbPb_2760_dndeta_%s_nosec_20140512_1628/"
				 "forward_dndeta.root", meth), "READ");
  TCollection* res = GetCollection(file, "ForwarddNdetaResults");
  Bool_t first = false;
  if (!centAxis) {
    centAxis = static_cast<TAxis*>(GetObject(res, "centAxis", TAxis::Class()));
    Info("ProcessOne", "Got centrality axis %p", centAxis);
    first = true;
  }

  Info("ProcessOne", "Doing stuff for %s (%sfirst)",
       meth, (first ? "" : "not "));
  for (Int_t i = 1; i <= centAxis->GetNbins(); i++) {
    Double_t c1 = centAxis->GetBinLowEdge(i);
    Double_t c2 = centAxis->GetBinUpEdge(i);
    if (meth[0] == 'Z' && c1 > 30) break;
    
    TString binName(Form("cent%03dd%02d_%03dd%02d",
			 Int_t(c1), Int_t(c1 *100) % 100,
			 Int_t(c2), Int_t(c2 *100) % 100));
    TCollection* bin   = GetCollection(res, binName);
    if (!bin) continue;
    
    TH1*         here  = GetH1(bin, "dndetaForward");
    if (!here) continue;
    
    THStack*             stack = GetStack(stacks, binName, false);
    THStack*             ratio = GetStack(ratios, binName, false);
    TParameter<double>*  pmin  = GetParam(mins, binName, false);
    TParameter<double>*  pmax  = GetParam(maxs, binName, false);
    TParameter<double>*  pavg  = GetParam(avgs, binName, false);
    if (!stack) {
      Info("ProcessOne", "Creating stacks %s", binName.Data());
      stack = new THStack(binName, Form("%3d-%3d%%", Int_t(c1), Int_t(c2)));
      stack->SetUniqueID(here->GetMarkerColor());
      stacks->AddAt(stack, i-1);
      ratio = static_cast<THStack*>(stack->Clone());
      ratios->AddAt(ratio, i-1);
      pmin = new TParameter<double>(binName, +1000.);
      pmax = new TParameter<double>(binName, -1000.);
      pavg = new TParameter<double>(binName, -1000.);
      mins->AddAt(pmin, i-1);
      maxs->AddAt(pmax, i-1);
      avgs->AddAt(pavg, i-1);
    }
    here->SetTitle(meth);
    here->SetMarkerStyle(style);
    stack->Add(here);
    

    if (first) continue;

    TH1* denom = static_cast<TH1*>(stack->GetHists()->At(0));
    TH1* num   = static_cast<TH1*>(here->Clone(Form("%s/%s",
						    meth, denom->GetName())));
    num->SetDirectory(0);
    num->Divide(denom);

    Double_t min = pmin->GetVal();
    Double_t max = pmax->GetVal();
    Double_t sum = 0;
    Int_t    cnt = 0;
    for (Int_t j = 1; j <= num->GetNbinsX(); j++) {
      Double_t c = num->GetBinContent(j);
      if (c < 1e-6) continue;
      min        =  TMath::Min(c, min);
      max        =  TMath::Max(c, max);
      sum        += TMath::Abs(1-c);
      cnt++;
      num->SetBinContent(j, c+(centAxis->GetNbins()-i));
    }
    ratio->Add(num);
    pmin->SetVal(min);
    pmax->SetVal(max);
    sum /= cnt;
    pavg->SetVal(TMath::Max(pavg->GetVal(), sum));

    Printf("Method %10s %16s: min=%7.5f max=%7.5f avg=%7.5f",
	   meth,stack->GetTitle(), min, max, sum);
			
  }
}

    
void
CentSysErr(Bool_t div=false, Bool_t zem=false)
{
  const char* dndeta = "1/#it{N} d#it{N}_{ch}/d#it{#eta}";
  TList* stacks   = new TList;
  TList* ratios   = new TList;
  TList* mins     = new TList;
  TList* maxs     = new TList;
  TList* avgs     = new TList;
  TAxis* centAxis = 0;
  TString all     = "";
  
  const char* meths[]  = { "V0M", "CL1", "TRK", (zem ? "ZEMVSZDC" : 0), 0 };
  Int_t       styles[] = { 20,    21,    22,    23,         0 };

  const char** pmeth   = meths;
  Int_t*       pstyle  = styles;

  TCanvas* c = new TCanvas("c", "C", 1200, 800);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->Divide(1,2,0,0);
  c->GetPad(1)->SetRightMargin(0.01);
  c->GetPad(2)->SetRightMargin(0.01);
  gStyle->SetOptTitle(0);

  TLegend* mleg = new TLegend(0.35, (meths[3] ? 0.05 : 0.3), 0.65, 0.95);
  mleg->SetFillStyle(0);
  mleg->SetBorderSize(0);
  mleg->SetTextFont(42);
  TLegend* leg = new TLegend(0.35, 0.15, 0.55, 0.95);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetNColumns(2);
  
 
  while (*pmeth) {
    ProcessOne(*pmeth, *pstyle, stacks, ratios, mins, maxs, avgs, centAxis);
    TLegendEntry* e = mleg->AddEntry("", *pmeth, "p");
    e->SetMarkerStyle(*pstyle);
    all.Append(Form("_%s", *pmeth));
    pmeth++;
    pstyle++;
    
  }

  TGraphAsymmErrors*  trend = new TGraphAsymmErrors(centAxis->GetNbins());
  trend->SetMarkerStyle(20);
  trend->SetFillColor(kBlue-9);
  trend->SetFillStyle(3001);
  trend->SetMarkerColor(kBlue-7);
  trend->SetLineColor(kBlue-9);
  TIter               nextS(stacks);
  TIter               nextR(ratios);
  TIter               nextN(mins);
  TIter               nextM(maxs);
  TIter               nextA(avgs);
  Bool_t              first   = true;
  THStack*            stack   = 0;
  THStack*            ratio   = 0;
  TParameter<double>* pmin    = 0;
  TParameter<double>* pmax    = 0;
  TParameter<double>* pavg    = 0;
  Int_t               cnt     = 0;
  Double_t            least   = +1000;
  Double_t            largest = -1000;
  while ((stack = static_cast<THStack*>(nextS())) &&
	 (ratio = static_cast<THStack*>(nextR())) &&
	 (pmin = static_cast<TParameter<double>*>(nextN())) &&
	 (pmax = static_cast<TParameter<double>*>(nextM())) &&
	 (pavg = static_cast<TParameter<double>*>(nextA()))
	 ) {
    c->cd(1);
    stack->Draw(first ? "nostack" : "nostack same");

    Int_t  y   = centAxis->GetNbins()-cnt++;
    UInt_t col = stack->GetUniqueID();
    TLegendEntry* e = leg->AddEntry("",
				    Form("%10s (+%d)",
					 stack->GetTitle(), y-1),
				    "f");
    e->SetFillColor(col);
    e->SetFillStyle(1001);
    e->SetLineColor(col);
    c->cd(2);
    ratio->Draw(first ? "nostack" : "nostack same");

    Double_t min = 100*(1 - pmin->GetVal());
    Double_t max = 100*(pmax->GetVal()-1);
    Double_t avg = 100*pavg->GetVal();
    Double_t c1  = centAxis->GetBinLowEdge(cnt);
    Double_t c2  = centAxis->GetBinUpEdge(cnt);
    Double_t cc  = (c2+c1)/2;
    trend->SetPoint(cnt-1, cc, avg);
    trend->SetPointError(cnt-1, cc-c1, c2-cc,
			 (min > 0 ? TMath::Abs(avg-min) : 0), max-avg);
    least = TMath::Min(avg, least);
    largest = TMath::Max(avg, largest);
    
    TString  errs = Form("#pm%4.1f%%", avg);
    if (div) {
      if (min > 0) errs = Form("{}^{+%4.1f}_{-%4.1f}%%", max, min);
      else         errs = Form("+%4.1f%%", TMath::Max(min,max));
    }
    if (leg->GetNColumns() > 1) {
      e = leg->AddEntry("", errs, "");
      e->SetFillColor(kWhite);
      e->SetFillStyle(0);
      e->SetLineColor(kWhite);
      e->SetMarkerColor(kWhite);
    }
    else {
      TLatex*  ltx = new TLatex(5.9, y, errs);
      ltx->SetTextAlign(32);
      ltx->SetTextFont(42);
      ltx->SetTextSize(0.04);
      ltx->SetTextColor(col+1);
      ltx->Draw();
    }
    
    if (!first) continue;
    
    ratio->SetMaximum(centAxis->GetNbins()+1.3);
    ratio->SetMinimum(.1);

    Double_t ts = 0.06;
    TH1* hr = ratio->GetHistogram();
    hr->SetXTitle("#it{#eta}");
    hr->SetYTitle("Ratio to V0M");
    hr->GetXaxis()->SetTitleSize(ts);
    hr->GetYaxis()->SetTitleSize(ts);
    hr->GetXaxis()->SetLabelSize(ts);
    hr->GetYaxis()->SetLabelSize(ts);
    hr->GetXaxis()->SetTitleOffset(1-4*ts);
    hr->GetYaxis()->SetTitleOffset(1-5*ts);
    hr->GetXaxis()->SetNdivisions(10);
    hr->GetYaxis()->SetNdivisions(10);
    
    TH1* hs = stack->GetHistogram();
    hs->SetXTitle("#it{#eta}");
    hs->SetYTitle(dndeta);
    hs->GetXaxis()->SetTitleSize(ts);
    hs->GetYaxis()->SetTitleSize(ts);
    hs->GetXaxis()->SetLabelSize(ts);
    hs->GetYaxis()->SetLabelSize(ts);
    hs->GetXaxis()->SetTitleOffset(1-4*ts);
    hs->GetYaxis()->SetTitleOffset(1-5*ts);
    hs->GetXaxis()->SetNdivisions(10);
    hs->GetYaxis()->SetNdivisions(10);
    

    first = false;
  }

  c->cd(1);
  mleg->Draw();
  c->cd(2);
  leg->Draw();
    
  c->Print(Form("cent_syserr_%s%s.pdf", (div ? "div" : "avg"), all.Data()));
  c->SaveAs(Form("cent_syserr_%s%s.root", (div ? "div" : "avg"), all.Data()));

  TCanvas* aux = new TCanvas("aux", "aix");
  aux->SetTopMargin(0.01);
  aux->SetRightMargin(0.01);
  
  trend->SetTitle();
  trend->Draw("a3p");
  trend->GetHistogram()->SetXTitle("Centrality [%]");
  trend->GetHistogram()->SetYTitle(Form("#delta#left[%s#right] [%%]",dndeta));

  TF1* f = new TF1("f", "[0]*x*x+[1]", 0, 100);
  f->SetLineColor(kMagenta+2);
  Double_t cl = centAxis->GetBinCenter(1);
  Double_t ch = centAxis->GetBinCenter(centAxis->GetNbins());
  f->SetParameters(largest/ch/ch, least);
  f->Draw("same");
  
  TLatex* ll = new TLatex(0.98, 0.98,
			  Form("Sys. Uncertainty on %s from Centrality",
			       dndeta));
  ll->SetTextAlign(33);
  ll->SetTextSize(0.05);
  ll->SetTextFont(42);
  ll->SetNDC();
  ll->Draw();

  TLatex* fl = new TLatex(0.2, 0.8,
			  Form("#delta#approx%5.1fc^{2}+%5.1f",
			       largest/ch/ch*100*100,
			       least));
  fl->SetTextAlign(13);
  fl->SetTextSize(0.04);
  fl->SetTextFont(42);
  fl->SetNDC();
  fl->Draw();
  
  aux->Print(Form("cent_syserr_%s%s_trend.pdf",
		  (div ? "div" : "avg"), all.Data()));
}

      

  
