/**
 * @file   DetailsCompare.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 20:57:27 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
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
#include <TPaletteAxis.h>
#else
class TH1;
class TArrayI;
class TLegend;
class THStack;
class TCanvas;
class TList;
class TVirtualPad;
class TFile;
#endif

//====================================================================
namespace {
  /** 
   * A guard to suppress messages 
   */
  struct SuppressGuard3
  {
    /** The previous message level */
    Int_t save = 0;
    /** 
     * Constructor 
     * 
     * @param lvl Level to suppress to 
     */
    SuppressGuard3(Int_t lvl=2000)
    {
      save = gErrorIgnoreLevel;
      gErrorIgnoreLevel = lvl;
    }
    /** 
     * Destructor 
     */
    ~SuppressGuard3()
    {
      gErrorIgnoreLevel = save;
    }
  };
}

/**
 * Compare details
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct DetailsComparer
{
  TCanvas*     fCanvas;
  TVirtualPad* fTop;
  TVirtualPad* fBody;
  TVirtualPad* fBottom;
  TFile*       fOutput;
  TLegend*     fLegend;
  Bool_t       f2D;
  
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
    SuppressGuard3 g;
    fCanvas->SaveAs(Form("compare/%s[", fCanvas->GetTitle()), "PDF");

    fOutput = TFile::Open("compare/result.root", "RECREATE");
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

    SuppressGuard3 g;
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
    fCanvas->Close();
    fOutput->Write();
    // fOutput->Close();
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
    fOutput->cd();
    fLegend->Write("legend");
    files->Write("files", TObject::kSingleKey);
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
	min = TMath::Min(min, c-1.1*e);
	max = TMath::Max(max, c+1.1*e);
      }
    }
    ratios->SetMinimum(TMath::Max(0.5,min));
    ratios->SetMaximum(TMath::Min(1.5,max));
    return ratios;
  }
  TH1* AddHisto(TDirectory*    dir,
		Int_t          dimen,
		const TString& binName,
		const char*    sub,
		THStack*       stack,
		Color_t        col)
  {
    TString hname;
    hname.Form("%s/%s%dd/%s%s",binName.Data(),sub, dimen,
	       (f2D && TString(sub).EqualTo("results") ? "full/" : ""),
	       stack->GetName());
    // Printf("Will get histogram %s", hname.Data());
    TH1* h = static_cast<TH1*>(dir->Get(hname));
    if (!h) {
      Warning("Addhisto", "No histogram %s in %s", hname.Data(),dir->GetName());
      return 0;;
    }
    stack->Add(h);
    Double_t min = TMath::Min(stack->GetMinimum("nostack"), h->GetMinimum());
    Double_t max = TMath::Max(stack->GetMaximum("nostack"), h->GetMaximum());
    stack->SetMaximum(max);
    stack->SetMinimum(min);
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
		   TList*         stacks)
  {
    TIter       nextD(files);
    TDirectory* dir = 0;
    Int_t       cnt = 0;
    while ((dir = static_cast<TDirectory*>(nextD()))) {
      TIter         nextS(stacks);
      THStack*      stack = 0;
      Color_t       col   = IndexColor(cnt);
      Bool_t        ok    = true;
      while ((stack = static_cast<THStack*>(nextS()))) 
	if (!AddHisto(dir, dimens[cnt], binName, sub, stack, col)) ok = false;
      if (!ok) continue;
      cnt++;
    }
    return cnt;
  }
  void DrawStack(TVirtualPad* mother,
		 Int_t        sub,
		 THStack*     stack,
		 TDirectory*  out, 
		 Bool_t       log,
		 Double_t     tbase=0.05)
  {
    if (!stack || !stack->GetHists()) {
      Warning("DrawStack", "No stack passed");
      return;
    }
    if (TString(stack->GetName()).EqualTo("scaleProj")) log = false;

    TDirectory* dir = out->mkdir(stack->GetName());
    dir->cd();
    stack->Write("full");
    
    TVirtualPad* p = mother->cd(sub);
    // gStyle->SetTitleFontSize(2*tbase/p->GetHNDC());
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetLeftMargin(0.15);
    p->SetBottomMargin(0.15);
    p->Divide(1,2,0,0);

    TVirtualPad* q = p->cd(1);
    q->SetRightMargin(f2D ? 0.1 : 0.01);
    q->SetLogx(log);
    q->SetLogy(log);
    q->SetGridx();
    q->SetGridy();
    q->SetTicks();
    TString tit(stack->GetTitle());
    stack->SetTitle("");
    stack->Draw(f2D ? "nostack colz" : "nostack");
    if (stack->GetHistogram()) {
      TH1* h = stack->GetHistogram();
      h->GetXaxis()->SetTitleSize(tbase/p->GetHNDC());
      h->GetXaxis()->SetLabelSize(tbase/p->GetHNDC());
      h->GetXaxis()->SetTitleOffset(0.4);
      h->GetXaxis()->SetNdivisions(207);
      h->GetYaxis()->SetTitleSize(tbase/p->GetHNDC());
      h->GetYaxis()->SetLabelSize(tbase/p->GetHNDC());
      h->GetYaxis()->SetNdivisions(207);
      h->GetYaxis()->SetTitleOffset(0.4);
      h->GetYaxis()->CenterTitle(true);
      h->SetYTitle(tit);
      if (f2D) {
	q->Update();
	h = static_cast<TH1*>(stack->GetHists()->Last());
	TPaletteAxis* pal = static_cast<TPaletteAxis*>(h->GetListOfFunctions()
						       ->FindObject("palette"));
	if (pal) {
	  pal->GetAxis()->SetTitleSize(tbase/p->GetHNDC());
	  pal->GetAxis()->SetLabelSize(tbase/p->GetHNDC());
	  pal->GetAxis()->SetNdivisions(207);
	  pal->GetAxis()->SetTitleOffset(0.4);
	}
      }
    }
    q->Modified();
    
    q = p->cd(2);
    q->SetRightMargin(f2D ? 0.1 : 0.01);
    q->SetLogx(log);
    q->SetTicks();
    q->SetGridx();
    q->SetGridy();
    
    THStack* ratios = RatioStack(stack);
    if (!ratios) {
      Warning("DrawStack", "Failed to make ratio");
      return;
    }
    // if (ratios->GetMaximum("nostack") > 1.5) ratios->SetMaximum(1.5);
    // if (ratios->GetMinimum("nostack") < 0.5) ratios->SetMinimum(0.5);
    ratios->Write("ratios");
    tit = ratios->GetTitle();
    ratios->SetTitle("");
    ratios->Draw(f2D ? "nostack colz" : "nostack");
    if (ratios->GetHistogram()) {
      TH1* h = ratios->GetHistogram();
      h->GetXaxis()->SetTitleSize(tbase/p->GetHNDC());
      h->GetXaxis()->SetLabelSize(tbase/p->GetHNDC());
      h->GetXaxis()->SetTitleOffset(0.4);
      h->GetXaxis()->SetNdivisions(207);
      h->GetYaxis()->SetTitleSize(tbase/p->GetHNDC());
      h->GetYaxis()->SetLabelSize(tbase/p->GetHNDC());
      h->GetYaxis()->SetTitleOffset(0.4);
      h->GetYaxis()->SetNdivisions(207);
      h->GetYaxis()->CenterTitle(true);
      h->SetYTitle("Ratio");
      if (f2D) {
	q->Update();
	h = static_cast<TH1*>(ratios->GetHists()->Last());
	TPaletteAxis* pal = static_cast<TPaletteAxis*>(h->GetListOfFunctions()
						       ->FindObject("palette"));
	if (pal) {
	  pal->GetAxis()->SetTitleSize(tbase/p->GetHNDC());
	  pal->GetAxis()->SetLabelSize(tbase/p->GetHNDC());
	  pal->GetAxis()->SetNdivisions(207);
	  pal->GetAxis()->SetTitleOffset(0.4);
	}
      }
    }
    q->Modified();
    p->Modified();
  }
  void DrawStacks(TVirtualPad* mother,
		  TList*       stacks,
		  TDirectory*  out,
		  Bool_t       log)
  {
    Int_t cnt = 1;
    TIter nextS(stacks);
    THStack* stack = 0;
    while ((stack = static_cast<THStack*>(nextS()))) {
      DrawStack(mother, cnt, stack, out, log);
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
		  TDirectory*    out, 
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
    
    
    if (FillStacks(files, dimens, binName, sub, stacks) <= 0) return;
    DrawStacks(fBody, stacks,out,  log);

    fBottom->cd();
    fLegend->Draw();

    PrintCanvas(title,Form("%s_%s", sub, binName.Data()));
  }
    
  void ProcessCentDelta(TList*         files,
			const TArrayI& dimens,
			const TString& binName,
			const TString& centTitle,
			TDirectory*    out)
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
    TDirectory* dir = out->mkdir("delta");
    
    ProcessOne(files, dimens, binName, "delta",
	       Form("#Delta #minus %s", centTitle.Data()),
	       names, titles, dir, true);
  }
  void ProcessCentResult(TList*         files,
			 const TArrayI& dimens,
			 const TString& binName,
			 const TString& centTitle,
			 TDirectory*    out)
  {
    const char* names[] = { "result", "simG", "realS", "simS",  
    			    "realM",  "simM", "realC", "simC",  0 };
    const char* titles[] = { "dN/d#eta", "G'", "S", "S'",
			     "M", "M'", "C", "C'", 0 };
    TDirectory* dir = out->mkdir("parts");
    ProcessOne(files, dimens, binName, "results",
	       Form("Parts #minus %s", centTitle.Data()),
	       names, titles, dir);
  }
      
				 
  void ProcessCentBin(TList*         files,
		      const TArrayI& dimens,
		      const TString& binName,
		      const TString& centTitle,
		      TDirectory*    out)
  {
    TDirectory* dir = out->mkdir(binName);
    ProcessCentDelta (files, dimens, binName, centTitle, dir);
    ProcessCentResult(files, dimens, binName, centTitle, dir);
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
      TString tmp(argv[i]);
      if (tmp.EndsWith(".C") ||
	  tmp.EndsWith(".C+") ||
	  tmp.BeginsWith("-")) continue;
      args[j] = StrDup(argv[i]);
      Printf("Argument %s", args[j]);
      j++;
    }
    args[j] = 0;
    // gApplication->ClearInputFiles(); 
    this->Run(const_cast<const char**>(args));
  }
  void Run(const char** filenames)
  {
    MakeCanvas();
    fLegend = new TLegend(.05, .05, .95, .95);
    fLegend->SetFillStyle(0);
    fLegend->SetBorderSize(0);

    TList*  files = new TList;
    TArrayI dimens;
    Int_t   cnt   = 0;
    TString last  = "";
    const char** fptr = filenames;    
    while (*fptr) {
      // Extract information 
      TString fn(*fptr);
      if (!fn.EndsWith(".root")) {
	last = fn;
	fptr++;
	continue;
      }
      TString dn(gSystem->DirName(fn));
      TString tit(last);
      if (tit.IsNull()) {
	tit = fn;
	tit.ReplaceAll(".root",""); 
	if (!dn.IsNull() && !dn.EqualTo(".")) tit = dn;
      }
      last = "";
      
      // Try o open the file 
      TFile* file = TFile::Open(fn, "READ");
      fptr++;
      if (!file) continue;

      // Set the title and increase count 
      file->SetTitle(tit);
      cnt++;

      // Create legend entry 
      Color_t       color = IndexColor(cnt-1);
      TLegendEntry* e     = fLegend->AddEntry("dummy", tit, "f");
      e->SetFillColor(color);
      e->SetLineColor(color);
      e->SetFillStyle(1001);

      // Figure out dimension 
      Int_t   dim = 0;
      if      (fn.EndsWith("etaipz.root") || dn.EndsWith("etaipz")) dim = 3;
      else if (fn.EndsWith("eta.root")    || dn.EndsWith("eta"))    dim = 2;
      else if (fn.EndsWith("const.root")  || dn.EndsWith("const"))  dim = 1;      
      dimens.Set(cnt);
      dimens[cnt-1] = dim;

      // Add to list of files 
      files->Add(file);
    }
    if (cnt <= 0) {
      Warning("", "No files opened");
      return;
    }
    f2D = false; // cnt == 2;
    TFile* f0 = static_cast<TFile*>(files->At(0));
    TH1* cent = static_cast<TH1*>(f0->Get("realCent"));
    // Printf("Got centrality histogram %p", cent);

    MakeTitlePage(files);
    
    for (Int_t i = 1; i <= cent->GetNbinsX(); i++) {
      Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
      Double_t c2 = cent->GetXaxis()->GetBinUpEdge(i);
      TString  cn, ct;
      cn.Form("cent%06.2f_%06.2f", c1, c2); cn.ReplaceAll(".","d");
      ct.Form("%5.1f-%5.1f%%", c1, c2);
      Printf("Processing centrality bin %6.2f-%6.2f%%", c1, c2);
      ProcessCentBin(files, dimens, cn, ct, fOutput);
    }
    CloseCanvas();
  }
};

/** 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
void
DetailsCompare()
{
  DetailsComparer d;
  d.Run(); // "hijing_etaipz.root", "eposlhc_etaipz.root", "whijing_const.root");
}
//
// EOF
//
