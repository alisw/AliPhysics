/**
 * @file   CompareResults.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:51:29 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
#ifndef __CINT__
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TCanvas.h>
# include <TH1.h>
# include <TFile.h>
# include <THStack.h>
# include <TError.h>
# include <TApplication.h>
#else
class TH1;
class TDirectory;
class TLegend;
#endif

//____________________________________________________________________
struct Compare
{
  Double_t fMin;
  Double_t fMax;
  TLegend* fLegend;
  Compare()
    : fMin(+1e9),
      fMax(-1e9)
  {}

  //____________________________________________________________________
  static TObject* GetO(TDirectory* dir, const char* name, TClass* cls=0)
  {
    if (!dir) {
      Warning("GetO", "No directory passed");
      return 0;
    }

    TObject* o = dir->Get(name);
    if (!o) {
      Warning("GetO", "object %s not found in %s",
	      name, dir->GetPath());
      return 0;
    }
    if (!cls) return o;
    if (!o->IsA()->InheritsFrom(cls)) {
      Warning("GetO", "Object %s in %s is not a %s, but a %s",
	      name, dir->GetPath(), cls->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  //____________________________________________________________________
  static TDirectory* GetD(TDirectory* dir, const char* name)
  {
    return static_cast<TDirectory*>(GetO(dir,name,TDirectory::Class()));
  }

  //____________________________________________________________________
  static TH1* GetH1(TDirectory* dir, const char* name)
  {
    return static_cast<TH1*>(GetO(dir,name,TH1::Class()));
  }
  //____________________________________________________________________
  TH1* One(TDirectory* newDir, TDirectory* oldDir, Double_t c1, Double_t c2)
  {
    TString name;
    name.Form("cent%03dd%02d_%03dd%02d",
	      Int_t(c1), Int_t(c1*100)%100,
	      Int_t(c2), Int_t(c2*100)%100);
    TDirectory* newSubDir = GetD(newDir, name);
    TDirectory* oldSubDir = GetD(oldDir, name);
    if (!newSubDir || !oldSubDir) return 0;
    Int_t newDim = 0;
    if      (TString(newDir->GetName()).Contains("etaipz")) newDim = 3;
    else if (TString(newDir->GetName()).Contains("eta"))    newDim = 2;
    else if (TString(newDir->GetName()).Contains("const"))  newDim = 1;
    Int_t oldDim = 0;
    if      (TString(oldDir->GetName()).Contains("etaipz")) oldDim = 3;
    else if (TString(oldDir->GetName()).Contains("eta"))    oldDim = 2;
    else if (TString(oldDir->GetName()).Contains("const"))  oldDim = 1;

    TDirectory* newSubSubDir = GetD(newSubDir, Form("results%dd",newDim));
    TDirectory* oldSubSubDir = GetD(oldSubDir, Form("results%dd",oldDim));
    if (!newSubSubDir || !oldSubSubDir) return 0;

    TH1* newRes = GetH1(newSubSubDir, "result");
    TH1* oldRes = GetH1(oldSubSubDir, "result");
    if (!newRes || !oldRes) return 0;

    TH1* ratio = static_cast<TH1*>(newRes->Clone(name));
    ratio->SetDirectory(0);
    ratio->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
    ratio->SetYTitle("New / Old");
    ratio->Divide(oldRes);
    fMin = TMath::Min(fMin, ratio->GetMinimum());
    fMax = TMath::Max(fMax, ratio->GetMaximum());

    Printf("Calculated %s/%s", newDir->GetName(), oldDir->GetName());
    if (!fLegend) return ratio;

    
    TLegendEntry* e =
      fLegend->AddEntry("", Form("%3.0f - %3.0f%%", c1, c2), "f");
    e->SetFillStyle(1001);
    e->SetFillColor(ratio->GetMarkerColor());
  
    return ratio;
  }
  //____________________________________________________________________
  void Run(const char* newName,        const char* oldName,
	   const char* newTitle="New", const char* oldTitle="Old")
  {
    TFile* newFile = TFile::Open(newName,"READ");
    TFile* oldFile = TFile::Open(oldName,"READ");
    if (!newFile || !oldFile) return;

    TH1* newCent = GetH1(newFile, "realCent");
    TH1* oldCent = GetH1(oldFile, "realCent");
    if (!newCent || !oldCent) return;

    TString  t; t.Form("#it{R}=#frac{%s}{%s}", newTitle, oldTitle);
    TCanvas* c     = new TCanvas("c", t, 1200, 800);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.20);
    fLegend = new TLegend(1-c->GetRightMargin(),
			  c->GetBottomMargin(),
			  1, 1-c->GetTopMargin(),
			  t);
    fLegend->SetFillStyle(0);
    fLegend->SetBorderSize(0);
    THStack* stack = new THStack("ratios","");
			       
    fMin = +1e6;
    fMax = -1e6;
    TH1* one = 0;
    for (Int_t i = newCent->GetNbinsX(); i--;) {
      Double_t c1 = newCent->GetXaxis()->GetBinLowEdge(i+1);
      Double_t c2 = newCent->GetXaxis()->GetBinUpEdge(i+1);
      Info("", "c1=%f c2=%f", c1, c2);
      TH1*     r  = One(newFile, oldFile, c1, c2);    
      if (!r) continue;
      if (!one) {
	one = static_cast<TH1*>(r->Clone("one"));
	one->SetDirectory(0);
	one->Reset();
	for (Int_t j = 1; j <= one->GetNbinsX(); j++) {
	  one->SetBinContent(j,1);
	  one->SetBinError  (j,0);
	}
      }
      // r->Add(one, i-1);
      // r->Scale(TMath::Power(10,i));
      stack->Add(r);
    }
    stack->Draw("nostack");
    stack->SetMinimum(0.95*fMin);
    stack->SetMaximum(1.05*fMax);
    stack->GetHistogram()->SetXTitle("#eta");
    stack->GetHistogram()->SetYTitle("#it{R}");
    fLegend->Draw();
    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs(Form("%sover%s.png", newTitle, oldTitle));
  }  
};

void CompareResults(const char* newFile, const char* oldFile,
		    const char* newTit,  const char* oldTit)
{
  Compare* c = new Compare;
  c->Run(newFile,oldFile,newTit,oldTit);
}

void CompareResults(const char** argv)
{
  TString newFile;
  TString oldFile;
  TString newTit("");
  TString oldTit("");
  const char** ptr = argv;
  while ((*ptr)) {
    TString argi = *ptr;
    ptr++;
    if (argi.Contains("help")) {
      Printf("Usage: CompareResults AFILE BFILE [ATITLTE [BTITLE]]");
      return;
    }
    if (argi.Contains("CompareResults.C")) continue;
    if (argi.BeginsWith("-")) continue;
    if (argi.EndsWith(".root")) {
      if (newFile.IsNull()) newFile = argi;
      else                  oldFile = argi;
    }
    else {
      if (newTit.IsNull())  newTit = argi;
      else  	            oldTit = argi;
    }
  }
  if (newTit.IsNull()) newTit = "New";
  if (oldTit.IsNull()) oldTit = "Old";
  
      
  CompareResults(newFile, oldFile, newTit, oldTit);
}

void CompareResults()
{
  CompareResults(const_cast<const char**>(&(gApplication->Argv()[1])));
  gApplication->ClearInputFiles();
}



//
// EOF
//

