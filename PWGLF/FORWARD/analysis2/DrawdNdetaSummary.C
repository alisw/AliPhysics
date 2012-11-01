/**
 * @file   DrawAODSummary.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 30 09:47:30 2012
 * 
 * @brief  Script to draw summary of AOD pass into a PDF 
 * 
 * 
 */
#ifndef __CINT__
# include <TCollection.h>
# include <TH1.h>
# include <TH2.h>
# include <THStack.h>
# include <TParameter.h>
# include <TFile.h>
# include <TCanvas.h>
# include <TStyle.h>
# include <TError.h>
# include <TLegend.h>
# include <TLatex.h>
# include <TLegendEntry.h>
# include <TPaveText.h>
# include <TMath.h>
# include <TSystem.h>
#else
class TCollection;
class TDirectory;
class TH1;
class TH2;
class THStack;
class TCanvas;
class TLatex;
class TVirtualPad;
class TLegend;
class TAxis;
#endif

//____________________________________________________________________
/** 
 * Find an object in a collection
 * 
 * @param parent Parent list
 * @param name   Name of object
 * 
 * @return Pointer to object or null 
 */
TObject* GetObject(const TCollection* parent, const TString& name)
{
  // Info("GetObject", "Getting object %s from %p", name.Data(), parent);
  // --- Check parent ------------------------------------------------
  if (!parent) {
    Warning("GetObject", "No parent list");
    return 0;
  }
  // --- Check name --------------------------------------------------
  if (name.IsNull()) { 
    Warning("GetObject", "No name specified");
    return 0;
  }
  // --- Find the object ---------------------------------------------
  TObject* o = parent->FindObject(name);
  if (!o) {
    Warning("GetObject", "Object \"%s\" not found in parent \"%s\"",
	    name.Data(), parent->GetName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
/** 
 * Find an object in a directory
 * 
 * @param parent Parent directory
 * @param name   Name of object
 * 
 * @return Pointer to object or null 
 */
TObject* GetObject(TDirectory* parent, const TString& name)
{
  // Info("GetObject", "Getting object %s from %p", name.Data(), parent);
  // --- Check parent ------------------------------------------------
  if (!parent) {
    Warning("GetObject", "No parent directory");
    return 0;
  }
  // --- Check name --------------------------------------------------
  if (name.IsNull()) { 
    Warning("GetObject", "No name specified");
    return 0;
  }
  // --- Find the object ---------------------------------------------
  TObject* o = parent->Get(name);
  if (!o) {
    Warning("GetObject", "Object \"%s\" not found in parent \"%s\"",
	    name.Data(), parent->GetName());
    return 0;
  }
  return o;
}

//____________________________________________________________________
/** 
 * Check the type of a found object 
 * 
 * @param o   Object 
 * @param cl  Class 
 * @param src Source of object
 * 
 * @return true on success, false otherwise 
 */
Bool_t CheckType(const TObject* o, const TClass* cl, const TString& src)
{
  // Info("CheckType", "Checking type of %s vs %s", o->GetName(), cl->GetName());
  if (!o->IsA()->InheritsFrom(cl)) { 
    Warning("CheckType", "Object \"%s\" retrieved from \"%s\" is not a "
	    "%s but a %s", o->GetName(), src.Data(), cl->GetName(), 
	    o->ClassName());
    return false;
  }
  return true;
}

//_____________________________________________________________________
void GetParameter(const TCollection* c, const TString& name, UShort_t& value)
{
  // Info("GetParameter", "Getting parameter of %s from %p", name.Data(), c);
  TObject* o = GetObject(c, name);
  if (!o) return;
  value = o->GetUniqueID();
}
//_____________________________________________________________________
void GetParameter(const TCollection* c, const TString& name, Int_t& value)
{
  // Info("GetParameter", "Getting parameter of %s from %p", name.Data(), c);
  TObject* o = GetObject(c, name);
  if (!o) return;
  value = o->GetUniqueID();
}
//_____________________________________________________________________
void GetParameter(const TCollection* c, const TString& name, Double_t& value)
{
  // Info("GetParameter", "Getting parameter of %s from %p", name.Data(), c);
  TObject* o = GetObject(c, name);
  if (!o) return;
  UInt_t  i = o->GetUniqueID();
  Float_t v = *reinterpret_cast<Float_t*>(&i);
  value = v;
}
//_____________________________________________________________________
void GetParameter(const TCollection* c, const TString& name, Bool_t& value)
{
  // Info("GetParameter", "Getting parameter of %s from %p", name.Data(), c);
  TObject* o = GetObject(c, name);
  if (!o) return;
  value = o->GetUniqueID();
}

//____________________________________________________________________
/** 
 * Find a collection in another collection 
 * 
 * @param parent Parent collection 
 * @param name   Name of the collection 
 * 
 * @return pointer to collection on success, otherwise null 
 */
TCollection* GetCollection(const TCollection* parent, const TString& name)
{
  // Info("GetCollection", "Getting collection of %s from %p", name.Data(), c);
  // --- Find the object ---------------------------------------------
  TObject* o = GetObject(parent, name);
  if (!o) return 0;
  
  // --- Check type of found object ----------------------------------
  if (!CheckType(o, TCollection::Class(), parent->GetName())) return 0;
  
  // --- Return the collection ---------------------------------------
  return static_cast<TCollection*>(o);
}

//____________________________________________________________________
/** 
 * Find a collection in a directory
 * 
 * @param parent Parent directory
 * @param name   Name of the collection 
 * 
 * @return pointer to collection on success, otherwise null 
 */
TCollection* GetCollection(TDirectory* parent, const TString& name)
{
  // Info("GetCollection", "Getting collection of %s from %p", 
  //       name.Data(), parent);
  // --- Find the object ---------------------------------------------
  TObject* o = GetObject(parent, name);
  if (!o) return 0;

  // --- Check the type of object ------------------------------------
  if (!CheckType(o, TCollection::Class(), parent->GetName())) return 0;
  
  // --- Return the collection ---------------------------------------
  return static_cast<TCollection*>(o);
}

//____________________________________________________________________
/** 
 * Get a 1D histogram from a collection
 * 
 * @param parent Parent collection 
 * @param name   Name of histogram 
 * 
 * @return pointer or null
 */
TH1* GetH1(const TCollection* parent, const TString& name)
{
  // Info("GetH1", "Getting 1D histogram of %s from %p", name.Data(), c);
  // --- Find the object ---------------------------------------------
  TObject* o = GetObject(parent, name);
  if (!o) return 0;

  // --- Check the type of object ------------------------------------
  if (!CheckType(o, TH1::Class(), parent->GetName())) return 0;
  
  // --- Return the collection ---------------------------------------
  return static_cast<TH1*>(o);
}
//____________________________________________________________________
/** 
 * Get a 2D histogram from a collection
 * 
 * @param parent Parent collection 
 * @param name   Name of histogram 
 * 
 * @return pointer or null
 */
TH2* GetH2(const TCollection* parent, const TString& name)
{
  // Info("GetH2", "Getting 2D histogram of %s from %p", name.Data(), c);
  // --- Find the object ---------------------------------------------
  TObject* o = GetObject(parent, name);
  if (!o) return 0;

  // --- Check the type of object ------------------------------------
  if (!CheckType(o, TH2::Class(), parent->GetName())) return 0;
  
  // --- Return the collection ---------------------------------------
  return static_cast<TH2*>(o);
}
//____________________________________________________________________
/** 
 * Get a histogram stack from a collection
 * 
 * @param parent Parent collection 
 * @param name   Name of histogram 
 * 
 * @return pointer or null
 */
THStack* GetStack(const TCollection* parent, const TString& name,
		  const char* sub=0)
{
  // Info("GetStack", "Getting histogram stack %s from %p", name.Data(), parent);
  // --- Find the object ---------------------------------------------
  TObject* o = GetObject(parent, name);
  if (!o) return 0;

  // --- Check the type of object ------------------------------------
  if (!CheckType(o, THStack::Class(), parent->GetName())) return 0;
  
  THStack* stack = static_cast<THStack*>(o);
  if (sub == 0) return stack;
  
  if (stack->GetHists()->GetEntries() <= 0 ||stack->GetMaximum() < 1) { 
    stack->GetHists()->Delete();
    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(parent, *ptr);
      if (!sc) { ptr++; continue; }

      TH2* h = GetH2(sc, sub);
      if (!h) continue;
      TH1* p = h->ProjectionX(*ptr, 1, h->GetNbinsY(), "e");
      p->Scale(1., "width");
      p->SetTitle(*ptr);
      p->SetDirectory(0);
      stack->Add(p);
      ptr++;
    }
  }
  // --- Return the collection ---------------------------------------
  return stack;
}

//____________________________________________________________________
void Pause()
{
  printf("Press enter to continue");
  std::cin.get();
}
  
//____________________________________________________________________
/** 
 * Clear canvas 
 * 
 * @param c Canvas to clear 
 *
 * @ingroup pwglf_forward_scripts_corr
 */
void
ClearCanvas(TCanvas* c)
{
  // Info("ClearCanvas", "Clearing canvas");
  c->SetLeftMargin(.1);
  c->SetRightMargin(.05);
  c->SetBottomMargin(.1);
  c->SetTopMargin(.05);
  c->Clear();

  Float_t dy = .05;
  TPad* p1 = new TPad("top", "Top", 0, 1-dy, 1, 1, 0, 0);
  p1->SetNumber(1);
  p1->SetFillColor(kBlue-5);
  p1->SetBorderSize(0);
  p1->SetBorderMode(0);
  c->cd();
  p1->Draw();

  TPad* p2 = new TPad("body", "Body", 0, 0, 1, 1-dy, 0, 0);
  p2->SetNumber(2);
  p2->SetFillColor(0);
  p2->SetFillStyle(0);
  p2->SetBorderSize(0);
  p2->SetBorderMode(0);
  c->cd();
  p2->Draw();
  p2->cd();

}
//____________________________________________________________________
/** 
 * Create a canvas 
 * 
 * @param pname Name of PDF file to make 
 * 
 * @return Created canvas 
 */
TCanvas* CreateCanvas(const TString& pname)
{
  // Info("CreateCanvas", "Creating canvas");
  Int_t size = 1000;
  TCanvas* c = new TCanvas("c", pname.Data(), size / TMath::Sqrt(2), size);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->Print(Form("%s[", pname.Data()));
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.5);
  gStyle->SetTitleY(1);
  gStyle->SetTitleW(.8);
  gStyle->SetTitleH(.09);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(1);
  gStyle->SetPalette(1);

  ClearCanvas(c);

  return c;
}

//____________________________________________________________________
/** 
 * Close the PDF
 * 
 * @param c Canvas 
 */
void CloseCanvas(TCanvas* c)
{
  // Info("CloseCanvas", "Closing canvas");
  ClearCanvas(c);
  c->Print(Form("%s]", c->GetTitle()));
}

//____________________________________________________________________
/** 
 * Print the canvas 
 * 
 * @param c      Canvas 
 * @param title  Title 
 */
void PrintCanvas(TCanvas* c, const TString& title, 
		 Float_t size=.7, Bool_t pause=false)
{
  // Info("PrintCanvas", "Printing page %s", title.Data());
  TString tit;
  tit.Form("Title:%s", title.Data());

  c->cd(1);
  TLatex* ltx = new TLatex(.5, .5, title);
  ltx->SetNDC();
  ltx->SetTextAlign(22);
  ltx->SetTextSize(size);
  ltx->SetTextColor(kWhite);
  ltx->SetTextFont(62);
  ltx->Draw();
  
  c->Modified();
  c->Update();
  c->cd();
  gSystem->RedirectOutput("/dev/null");
  c->Print(c->GetTitle(), tit);
  gSystem->RedirectOutput(0);
  
  // std::cin.peek();
  if (pause) Pause();

  ClearCanvas(c);
}
//____________________________________________________________________
/** 
 * Make a chapter page 
 * 
 * @param c     Canvas 
 * @param title Title 
 */
void MakeChapter(TCanvas* c, const TString& title)
{
  c->cd(2);
  
  // Info("MakeChapter", "Making chapter %s", title.Data());
  TLatex* ltx = new TLatex(.5, .5, title);
  ltx->SetNDC();
  ltx->SetTextAlign(22);
  ltx->Draw();

  PrintCanvas(c, title);
}
//____________________________________________________________________
void DrawInPad(TVirtualPad* c, Int_t padNo, TObject* h, Option_t* opts="",
	       UShort_t flags=0x0)
{
  // Info("DrawInPad", "Drawing %p in pad # %d of %p w/options %s, flags 0x%x", 
  //      h, padNo, c, opts, flags);
  TVirtualPad* p = c->cd(padNo);
  if (!p) { 
    Warning("DrawInPad", "Pad # %d not found in %s", padNo, c->GetName());
    return;
  }
  if (flags & 0x1) p->SetLogx();
  if (flags & 0x2) p->SetLogy();
  if (flags & 0x4) p->SetLogz();
  p->SetFillColor(0);
  TString o(opts);
  if (o.Contains("colz", TString::kIgnoreCase)) 
    p->SetRightMargin(0.15);
  if (o.Contains("text", TString::kIgnoreCase)) { 
    TH1* hh = static_cast<TH1*>(h);
    hh->SetMaximum(1.1*hh->GetMaximum());
    hh->SetMarkerSize(2);
    o.Append("30");
  }
  
  if (!h) {
    Warning("DrawInPad", "Nothing to draw in pad # %d", padNo);
    return;
  }
  h->Draw(o);

  if (flags& 0x10) { 
    TLegend* l = p->BuildLegend();
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
  }
  p->Modified();
  p->Update();
  p->cd();
}
//____________________________________________________________________
void CreateTemplates(TLatex*& name, TLatex*& value, Float_t size=.03)
{
  Double_t x1 = .1;
  Double_t x2 = .6;
  Double_t y  = .8;
  name = new TLatex(x1, y, "");
  name->SetTextAlign(13);
  name->SetNDC();
  name->SetTextSize(size);

  value = new TLatex(x2, y, "");
  value->SetTextAlign(13);
  value->SetNDC();
  value->SetTextSize(size);
}
  
//____________________________________________________________________
void DrawParameter(TLatex* name, TLatex* value, Double_t& y, 
		   const TString& sName, const TString& sValue)
{
  if (sName.IsNull() && sValue.IsNull()) return;
  if (!sName.IsNull())
    name->DrawLatex(name->GetX(), y, Form("%s:", sName.Data()));
  if (!sValue.IsNull())
    value->DrawLatex(value->GetX(), y, sValue.Data());
  y -= name->GetTextSize() + .02;
}  


//____________________________________________________________________
void DrawCentSum(const TCollection* sums, TCanvas* can, const TString& base, 
		 Int_t cLow, Int_t cHigh)
{
  TString folder; 
  if (cLow < 0 || cHigh < 0 || cLow >= cHigh) {
    folder = "all";
    cHigh  *= -1;
  }
  else folder.Form("cent%03d_%03d", cLow, cHigh);
  
  TCollection* c = GetCollection(sums, folder);
  if (!c) return;

  TVirtualPad* body = can->cd(2);
  body->Divide(2, 2);

  DrawInPad(body, 1, GetH1(c, "triggers"),                   "HIST TEXT");
  DrawInPad(body, 2, GetH1(c, Form("%sEvents",base.Data())), "HIST TEXT");
  DrawInPad(body, 3, GetH2(c, base.Data()),                  "colz");
  DrawInPad(body, 4, GetH2(c, Form("%s0", base.Data())),     "colz");
  
  PrintCanvas(can, Form("%s sums: %3d%% - %3d%%", base.Data(), cLow, cHigh));
}

//____________________________________________________________________
void DrawSums(TDirectory* top, const TString& base, TCanvas* can, bool onlyMB)
{
  TCollection* c = GetCollection(top, Form("%sSums", base.Data()));
  if (!c) return;

  TAxis* centAxis = static_cast<TAxis*>(GetObject(c, "centAxis"));

  TVirtualPad* body = can->cd(2);
  body->Divide(1, 2);
  
  body->cd(1);
  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value);
  for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
    DrawParameter(name, value, y, (i == 1 ? "Centrality classes" : ""),
		  Form("%3d%% - %3d%%", 
		       Int_t(centAxis->GetBinLowEdge(i)), 
		       Int_t(centAxis->GetBinUpEdge(i))));
  }
  Int_t sys, sNN, scheme, trigger;
  GetParameter(c, "sNN",     sNN); 
  GetParameter(c, "sys",     sys); 
  GetParameter(c, "scheme",  scheme); 
  GetParameter(c, "trigger", trigger); 
  DrawParameter(name, value, y, "Collision system", 
		(sys == 1 ? "pp" : (sys == 2 ? "PbPb" : (sys == 3 ? "pPb" : 
							 "unknown"))));
  DrawParameter(name, value, y, "#sqrt{s_{NN}}", Form("%4dGeV", sNN));
  DrawParameter(name, value, y, "Normalization scheme", Form("0x%x", scheme));
  DrawParameter(name, value, y, "Triggers",      Form("0x%x", trigger));
  
  DrawInPad(body, 2, GetH1(c, "cent"));

  PrintCanvas(can, Form("%s sums", base.Data()));

  DrawCentSum(c, can, base, centAxis->GetXmin(), -centAxis->GetXmax());
  if (onlyMB) return;

  for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
    DrawCentSum(c, can, base, centAxis->GetBinLowEdge(i), 
		centAxis->GetBinUpEdge(i));
}

//____________________________________________________________________
void CleanStack(THStack* stack, TLegend* l, const TAxis* axis)
{
  TList* hists = stack->GetHists();
  // Clean up list of histogram.  Histograms with no entries or 
  // no functions are deleted.  We have to do this using the TObjLink 
  // objects stored in the list since ROOT cannot guaranty the validity 
  // of iterators when removing from a list - tsck.  Should just implement
  // TIter::Remove(). 
  TObjLink* lnk = hists->FirstLink();
  Int_t j = 0;
  while (lnk) {
    TH1* h = static_cast<TH1*>(lnk->GetObject());
    bool remove = false;
    TString name(h->GetName());
    if (name.Contains("_mirror")) {
      // AliWarning(Form("No entries in %s - removing", h->GetName()));
      remove = true;
    }
    else if (l) { 
      j++;
      name.Form("%3d%% - %3d%%", 
		Int_t(axis->GetBinLowEdge(j)), 
		Int_t(axis->GetBinUpEdge(j)));
      TLegendEntry* e = l->AddEntry("dummy", name, "f");
      e->SetFillStyle(1001);
      e->SetFillColor(h->GetMarkerColor());
    }
    if (remove) {
      TObjLink* keep = lnk->Next();
      hists->Remove(lnk);
      lnk = keep;
      continue;
    }
    lnk = lnk->Next();
  }
}

//____________________________________________________________________
void DrawCentRes(const TCollection* sums, TCanvas* can, const TString& base, 
		 Int_t cLow, Int_t cHigh)
{
  TString folder; 
  if (cLow < 0 || cHigh < 0 || cLow >= cHigh) {
    folder =  "all";
    cHigh  *= -1;
  }
  else folder.Form("cent%03d_%03d", cLow, cHigh);
  
  TCollection* c = GetCollection(sums, folder);
  if (!c) return;

  TVirtualPad* body = can->cd(2);
  body->Divide(2, 3);

  DrawInPad(body, 1, GetH1(c, "triggers"), "HIST TEXT");
  DrawInPad(body, 2, GetH1(c, Form("norm%s",base.Data())));
  DrawInPad(body, 3, GetH1(c, Form("dndeta%s",base.Data())));
  DrawInPad(body, 4, GetH1(c, Form("dndeta%s_rebin05",base.Data())));
  DrawInPad(body, 5, GetH2(c, Form("d2Ndetadphi%s", base.Data())),"colz");
  
  TObject*   normCalc = GetObject(c, "normCalc");
  TString    calc     = normCalc->GetTitle();
  TObjArray* lines    = calc.Tokenize("\n");
  TPaveText* disp     = new TPaveText(.1,.1,.9,.9, "NDC");
  TIter      next(lines);
  TObject*   line     = 0;
  while ((line = next())) 
    disp->AddText(line->GetName());
  disp->SetBorderSize(0);
  disp->SetBorderSize(0);
  disp->SetFillStyle(0);
  DrawInPad(body, 6, disp);

  PrintCanvas(can, Form("%s result: %3d%% - %3d%%", base.Data(), cLow, cHigh));
}  

//____________________________________________________________________
THStack* 
DrawRes(TDirectory* top, const TString& base, TCanvas* can, Bool_t onlyMB)
{
  TCollection* c = GetCollection(top, Form("%sResults", base.Data()));
  if (!c) return 0;

  TAxis* centAxis = static_cast<TAxis*>(GetObject(c, "centAxis"));

  can->cd(2);

  TLatex*  name;
  TLatex*  value;
  Double_t y = .9;
  CreateTemplates(name, value, .02);
  for (Int_t i = 1; i <= centAxis->GetNbins(); i++) { 
    DrawParameter(name, value, y, (i == 1 ? "Centrality classes" : ""),
		  Form("%3d%% - %3d%%", 
		       Int_t(centAxis->GetBinLowEdge(i)), 
		       Int_t(centAxis->GetBinUpEdge(i))));
  }

  DrawParameter(name, value, y, "Collision system", 
		GetObject(c, "sys")->GetTitle());
  DrawParameter(name, value, y, "#sqrt{s_{NN}}",GetObject(c,"sNN")->GetTitle());
  DrawParameter(name, value, y, "trigger",GetObject(c,"trigger")->GetTitle());
  DrawParameter(name, value, y, "scheme", GetObject(c,"scheme")->GetTitle());

  Double_t epsT, epsT0;
  GetParameter(c, "triggerEff",  epsT);
  GetParameter(c, "triggerEff0", epsT0);
  DrawParameter(name, value, y, "#epsilon_{T}", Form("%f", epsT));
  DrawParameter(name, value, y, "#epsilon_{T,zero bin}", Form("%f", epsT0));

  PrintCanvas(can, Form("%s results", base.Data()));
  
  TVirtualPad* body = can->cd(2);
  body->Divide(1, 3);

  TLegend* l = new TLegend(0.1, 0.1, 0.9, 0.9, "Centralities");
  l->SetNColumns(2);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  THStack* dndeta = GetStack(c, "dndeta");
  CleanStack(dndeta, l, centAxis);
  THStack* dndeta5 = GetStack(c, "dndeta_rebin05");
  CleanStack(dndeta5, 0, 0);

  DrawInPad(body, 1, dndeta,  "nostack");
  DrawInPad(body, 2, dndeta5, "nostack");
  DrawInPad(body, 3, l, "");

  PrintCanvas(can, Form("%s results - stacks", base.Data()));
  
  DrawCentRes(c, can, base, centAxis->GetXmin(), -centAxis->GetXmax());
  if (onlyMB) return dndeta;

  for (Int_t i = 1; i <= centAxis->GetNbins(); i++) 
    DrawCentRes(c, can, base, centAxis->GetBinLowEdge(i), 
		centAxis->GetBinUpEdge(i));

  return dndeta;
}

  
//____________________________________________________________________
void DrawdNdetaSummary(const char* fname="forward_dndeta.root",
		       bool onlyMB=true)
{
  // --- Open the file -----------------------------------------------
  TString filename(fname);
  TFile* file = TFile::Open(filename.Data(), "READ");
  if (!file) { 
    Error("DrawAODSummary", "Failed to open \"%s\"", filename.Data());
    return;
  }

  // --- Make our canvas ---------------------------------------------
  TString pdfName(filename);
  pdfName.ReplaceAll(".root", ".pdf");

  TCanvas* c = CreateCanvas(pdfName);

  // --- Do each sub-algorithm ---------------------------------------
  DrawSums(file, "Forward", c, onlyMB);
  THStack* rF = DrawRes(file, "Forward", c, onlyMB);

  DrawSums(file, "Central", c, onlyMB);
  THStack* rC = DrawRes(file, "Central", c, onlyMB);

  TIter next(rF->GetHists());
  TH1*  h  = 0;
  while ((h = static_cast<TH1*>(next()))) rC->Add(h);
  
  rC->Draw("nostack");
  PrintCanvas(c, "Both");
    
  CloseCanvas(c);
}
  

  
