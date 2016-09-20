#ifndef __CINT__
#include <TH1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TList.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TError.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TApplication.h>
#include <TSystem.h>
#else
class TH1;
class TArrayI;
class TLegend;
class THStack;
class TCanvas;
class TList;
class TVirtualPad;
#endif

struct DetailsComparer
{
  TCanvas*     fCanvas;
  TVirtualPad* fTop;
  TVirtualPad* fBody;
  TVirtualPad* fBottom;
  
  void MakeCanvas()
  {
    fCanvas = new TCanvas("c","summary.pdf", 800, 1000);
    fTop    = new TPad("top",   "",0,.9,1,1);
    fBody   = new TPad("body",  "",0,.2,1,.9);
    fBottom = new TPad("bottom","",0,0, 1,.2);
    fTop->SetFillColor(kGray);
    fBody->SetFillColor(0);
    fBottom->SetFillColor(0);

    fCanvas->cd();
    fTop->Draw();
    fCanvas->cd();
    fBody->Draw();
    fCanvas->cd();
    fBottom->Draw();
    fCanvas->cd();
    gSystem->Exec("rm -rf compare");
    gSystem->Exec("mkdir -p compare");
    fCanvas->SaveAs(Form("compare/%s[", fCanvas->GetTitle()), "PDF");
  }
  void ClearCanvas()
  {
    fTop->Clear();
    fBody->Clear();
    fBottom->Clear();
  }
  void PrintCanvas(const char* title, const char* shortTitle=0)
  {
    TLatex* ltx = new TLatex(.5,.5,title);
    ltx->SetNDC();
    ltx->SetTextSize(.6);
    ltx->SetTextAlign(22);
    fTop->cd();
    ltx->Draw();

    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
    fCanvas->SaveAs(Form("compare/%s", fCanvas->GetTitle()),
		    Form("PDF title=%s", title));
    if (shortTitle && shortTitle[0] != '\0')
      fCanvas->SaveAs(Form("compare/%s.png", shortTitle));
    fCanvas->WaitPrimitive();
  }
  void CloseCanvas()
  {
    fCanvas->SaveAs(Form("compare/%s]", fCanvas->GetTitle()), "PDF");
    fCanvas->Delete();
  }  
  Color_t IndexColor(Int_t cnt)
  {
    static Color_t colors[] = { kRed+2,
				kOrange+2,
				kYellow+2,
				kSpring+2,
				kGreen+2,
				kTeal+2,
				kCyan+2,
				kBlue+2,
				kViolet+2,
				kMagenta+2,
				kPink+2,
				kBlack };
    return colors[(cnt/2+(cnt%2)*6) % 12];
  }
  void MakeTitlePage(const TList* files)
  {
    ClearCanvas();
    fBody->cd();
    
    TLatex* l1 = new TLatex(.1,.95, "Files");
    l1->SetTextFont(62);
    l1->SetTextAlign(13);
    l1->SetTextSize(0.06);
    l1->SetNDC();
    l1->Draw();

    TLatex* l2 = new TLatex(.15,.90,"");
    l2->SetTextFont(82);
    l2->SetTextAlign(13);
    l2->SetTextSize(0.05);
    l2->SetNDC();
    TLatex* l3 = new TLatex(.18,.90,"");
    l3->SetTextFont(42);
    l3->SetTextAlign(13);
    l3->SetTextSize(0.04);
    l3->SetNDC();
    TIter    next(files);
    TFile*   file = 0;
    Double_t y    = l2->GetY();
    while ((file = static_cast<TFile*>(next()))) {
      TString tit = file->GetTitle();
      TString nam = file->GetName();
      if (file == files->First()) nam.Append(" *");
      l2->DrawLatex(l2->GetX(), y, nam);      
      if (!nam.Contains(tit)) {
	y -= 1.1*l3->GetTextSize();
	l3->DrawLatex(l3->GetX(), y, tit);
      }
      y -= 1.1*l2->GetTextSize();
    }
    y -= 1.1*l2->GetTextSize();
    TLatex* l4 = new TLatex(.1, y, "* Other files compared to this file");
    l4->SetTextFont(42);
    l4->SetTextAlign(13);
    l4->SetTextSize(0.05);
    l4->SetNDC();
    l4->Draw();
    
    PrintCanvas("Comparison of details", "");
  }
    
	       
  THStack* RatioStack(THStack* stack)
  {
    TH1*     h0     = static_cast<TH1*>(stack->GetHists()->First());
    THStack* ratios = new THStack("ratios", Form("to %s", h0->GetName()));
    TIter    nextH (stack->GetHists());
    TH1*     hr     = 0;
    Double_t min    = +1e9;
    Double_t max    = -1e9;
    while ((hr = static_cast<TH1*>(nextH()))) {
      if (hr == h0) continue;
      hr = static_cast<TH1*>(hr->Clone());
      hr->Divide(h0);
      hr->SetDirectory(0);
      ratios->Add(hr);
      for (Int_t i = 1; i <= hr->GetNbinsX(); i++) {
	Double_t c = hr->GetBinContent(i);
	Double_t e = hr->GetBinError(i);
	if (c < 1e-6 || e < 1e-9) continue;
	min = TMath::Min(min, c-2*e);
	max = TMath::Max(max, c+2*e);
      }
    }
    ratios->SetMinimum(min);
    ratios->SetMaximum(max);
    return ratios;
  }
  void DrawStack(TVirtualPad* mother,
		 Int_t        sub,
		 THStack*     stack,
		 Bool_t       log,
		 Double_t     tbase=0.05)
  {
    if (!stack || !stack->GetHists()) {
      Warning("DrawStack", "No stack passed");
      return;
    }
    if (TString(stack->GetName()).EqualTo("scaleProj")) log = false;

    TVirtualPad* p = mother->cd(sub);
    gStyle->SetTitleFontSize(2*tbase/p->GetHNDC());
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetLeftMargin(0.15);
    p->SetBottomMargin(0.15);
    p->Divide(1,2,0,0);

    TVirtualPad* q = p->cd(1);
    q->SetRightMargin(0.01);
    q->SetLogx(log);
    q->SetLogy(log);
    q->SetGridx();
    q->SetGridy();
    q->SetTicks();
    TString tit(stack->GetTitle());
    stack->SetTitle("");
    stack->Draw("nostack");
    if (stack->GetHistogram()) {
      stack->GetHistogram()->GetXaxis()->SetTitleSize(tbase/p->GetHNDC());
      stack->GetHistogram()->GetXaxis()->SetLabelSize(tbase/p->GetHNDC());
      stack->GetHistogram()->GetXaxis()->SetTitleOffset(0.4);
      stack->GetHistogram()->GetXaxis()->SetNdivisions(207);
      stack->GetHistogram()->GetYaxis()->SetTitleSize(tbase/p->GetHNDC());
      stack->GetHistogram()->GetYaxis()->SetLabelSize(tbase/p->GetHNDC());
      stack->GetHistogram()->GetYaxis()->SetNdivisions(207);
      stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
      stack->GetHistogram()->GetYaxis()->CenterTitle(true);
      stack->GetHistogram()->SetYTitle(tit);
    }
    q->Modified();

    q = p->cd(2);
    q->SetRightMargin(0.01);
    q->SetLogx(log);
    q->SetTicks();
    q->SetGridx();
    q->SetGridy();
    
    THStack* ratios = RatioStack(stack);
    if (!ratios) {
      Warning("DrawStack", "Failed to make ratio");
      return;
    }
    tit = ratios->GetTitle();
    ratios->SetTitle("");
    ratios->Draw("nostack");
    if (ratios->GetHistogram()) {
      ratios->GetHistogram()->GetXaxis()->SetTitleSize(tbase/p->GetHNDC());
      ratios->GetHistogram()->GetXaxis()->SetLabelSize(tbase/p->GetHNDC());
      ratios->GetHistogram()->GetXaxis()->SetTitleOffset(0.4);
      ratios->GetHistogram()->GetXaxis()->SetNdivisions(207);
      ratios->GetHistogram()->GetYaxis()->SetTitleSize(tbase/p->GetHNDC());
      ratios->GetHistogram()->GetYaxis()->SetLabelSize(tbase/p->GetHNDC());
      ratios->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
      ratios->GetHistogram()->GetYaxis()->SetNdivisions(207);
      // ratios->GetHistogram()->SetYTitle(tit);
    }
    q->Modified();
    p->Modified();
    
  }
  TH1* AddHisto(TDirectory*    dir,
		Int_t          dimen,
		const TString& binName,
		const char*    sub,
		THStack*       stack,
		Color_t        col)
  {
    TString hname;
    hname.Form("%s/%s%dd/%s",binName.Data(),sub, dimen, stack->GetName());
    TH1* h = static_cast<TH1*>(dir->Get(hname));
    if (!h) {
      Warning("Addhisto", "No histogram %s in %s", hname.Data(),dir->GetName());
      return 0;;
    }
    stack->Add(h);
    TString fn(dir->GetName());
    TString dn(gSystem->DirName(fn));
    h->SetName(dn.IsNull() || dn == "." ? fn : dn);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    return h;
  }
  Int_t FillStacks(TList*         files,
		   const TArrayI& dimens,
		   const TString& binName,
		   const char*    sub,
		   TList*         stacks,
		   TLegend*       l)
  {
    TIter       nextD(files);
    TDirectory* dir = 0;
    Int_t       cnt = 0;
    while ((dir = static_cast<TDirectory*>(nextD()))) {
      TIter         nextS(stacks);
      THStack*      stack = 0;
      Color_t       col   = IndexColor(cnt);
      TString       fn    = dir->GetTitle();
      TLegendEntry* e     = l->AddEntry("dummy", fn, "f");
      e->SetFillColor(col);
      e->SetFillStyle(1001);
      Bool_t        ok    = true;
      while ((stack = static_cast<THStack*>(nextS()))) 
	if (!AddHisto(dir, dimens[cnt], binName, sub, stack, col)) ok = false;
      if (!ok) continue;
      cnt++;
    }
    return cnt;
  }
  void DrawStacks(TVirtualPad* mother, TList* stacks, Bool_t log)
  {
    Int_t cnt = 1;
    TIter nextS(stacks);
    THStack* stack = 0;
    while ((stack = static_cast<THStack*>(nextS()))) {
      DrawStack(mother, cnt, stack, log);
      cnt++;
    }
  }    
  void ProcessOne(TList*         files,
		  const TArrayI& dimens,
		  const TString& binName,
		  const char*    sub,
		  const TString& title,
		  const char**   names,
		  const char**   titles,
		  Bool_t         log=false)
  {
    TList* stacks = new TList;
    const char** pname = names;
    const char** ptit  = titles;
    while (*pname) {
      TString n(*pname);
      TString t(*ptit);
      pname++;
      ptit++;
      THStack* stack = new THStack(n, t);
      stacks->Add(stack);
    }
    if (stacks->GetEntries() <= 0) return;
    ClearCanvas();
    
    fBody->SetTopMargin(0.01);
    fBody->SetRightMargin(0.01);
    Int_t nCol = 2;
    Int_t nRow = (stacks->GetEntries()+1)/2;
    fBody->Divide(nCol, nRow, 0, 0);
    
    
    TLegend*    l = new TLegend(.05, .05, .95, .95);
    l->SetFillStyle(0);
    l->SetBorderSize(0);

    if (FillStacks(files, dimens, binName, sub, stacks, l) <= 0) return;
    DrawStacks(fBody, stacks, log);

    fBottom->cd();
    l->Draw();

    PrintCanvas(title,Form("%s_%s", sub, binName.Data()));
  }
    
  void ProcessCentDelta(TList*         files,
			const TArrayI& dimens,
			const TString& binName,
			const TString& centTitle)
  {
    const char* names[] = { "realDeltaM", "realDeltaI",
			    "simDeltaM" , "simDeltaI" ,
			    "simDeltaC" , "simDeltaP" ,
			    "simDeltaS" , "scaleProj",
			    0 };
    const char* titles[] = { "#Delta_{M}", "#Delta_{I}",
			     "#Delta_{M'}", "#Delta_{I'}",
			     "#Delta_{C'}", "#Delta_{P'}",
			     "#Delta_{S'}", "#LTk#it{k}#GT",
			     0 };

    ProcessOne(files, dimens, binName, "delta",
	       Form("#Delta #minus %s", centTitle.Data()),
	       names, titles, true);
  }
  void ProcessCentResult(TList*         files,
			 const TArrayI& dimens,
			 const TString& binName,
			 const TString& centTitle)
  {
    const char* names[] = { "result", "simG", "realS", "simS",  
    			    "realM",  "simM", "realC", "simC",  0 };
    const char* titles[] = { "dN/d#eta", "G'", "S", "S'",
			     "M", "M'", "C", "C'", 0 };
    ProcessOne(files, dimens, binName, "results",
	       Form("Parts #minus %s", centTitle.Data()),
	       names, titles);
  }
      
				 
  void ProcessCentBin(TList*         files,
		      const TArrayI& dimens,
		      const TString& binName,
		      const TString& centTitle)
  {
    ProcessCentDelta (files, dimens, binName, centTitle);
    ProcessCentResult(files, dimens, binName, centTitle);
  }
    
  void Run(const char* f1,   const char* f2,
	   const char* f3=0, const char* f4=0,
	   const char* f5=0, const char* f6=0,
	   const char* f7=0, const char* f8=0)
  {
    const char* filenames[] = { f1, f2, f3, f4,
				f5, f6, f7, f8,
				0 };
    Run(filenames);
  }
  void Run()
  {
    Int_t        argc = gApplication->Argc();
    char**       argv = gApplication->Argv();
    char**       args = new char*[argc];
    Int_t        j    = 0;
    for (Int_t i = 1; i < argc; i++) {
      if (!TString(argv[i]).EndsWith(".root")) continue;
      args[j] = StrDup(argv[i]);
      Printf("Got file %s", args[j]);
      j++;
    }
    args[j] = 0;
    // gApplication->ClearInputFiles(); 
    this->Run(const_cast<const char**>(args));
  }
  void Run(const char** filenames)
  {
    MakeCanvas();
    TList*  files = new TList;
    TArrayI dimens;
    Int_t   cnt   = 0;
    const char** fptr = filenames;    
    while (*fptr) {
      TString fn(*fptr);
      TString dn(gSystem->DirName(fn));
      TString tit(fn); tit.ReplaceAll(".root","");
      if (!dn.IsNull() && !dn.EqualTo(".")) tit = dn;
      Int_t   col = tit.Index(":");
      if (col != kNPOS) {
	tit.Remove(col, tit.Length()-col);
	fn.Remove(0,fn.Index(":")+1);
      }
      Printf("Got file %s with title %s (%d)", fn.Data(), tit.Data(), col);
      TFile* file = TFile::Open(fn, "READ");
      file->SetTitle(tit);
      fptr++;
      if (!file) continue;
      cnt++;
      Int_t   dim = 0;
      if      (fn.EndsWith("etaipz.root") || dn.EndsWith("etaipz")) dim = 3;
      else if (fn.EndsWith("eta.root")    || dn.EndsWith("eta"))    dim = 2;
      else if (fn.EndsWith("const.root")  || dn.EndsWith("const"))  dim = 1;
      
      dimens.Set(cnt);
      dimens[cnt-1] = dim;
      files->Add(file);
    }
    if (cnt <= 0) {
      Warning("", "No files opened");
      return;
    }
    TFile* f0 = static_cast<TFile*>(files->At(0));
    TH1* cent = static_cast<TH1*>(f0->Get("realCent"));
    Printf("Got centrality histogram %p", cent);

    MakeTitlePage(files);
    
    for (Int_t i = 1; i <= cent->GetNbinsX(); i++) {
      Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
      Double_t c2 = cent->GetXaxis()->GetBinUpEdge(i);
      TString  cn, ct;
      cn.Form("cent%06.2f_%06.2f", c1, c2); cn.ReplaceAll(".","d");
      ct.Form("%5.1f-%5.1f%%", c1, c2);
      Printf("Processing centrality bin %6.2f-%6.2f%%", c1, c2);
      ProcessCentBin(files, dimens, cn, ct);
    }
    CloseCanvas();
  }
};

void
DetailsCompare()
{
  DetailsComparer d;
  d.Run(); // "hijing_etaipz.root", "eposlhc_etaipz.root", "whijing_const.root");
}
