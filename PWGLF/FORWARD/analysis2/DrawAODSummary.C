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
# include <THStack.h>
# include <TH1.h>
# include <TH2.h>
# include <TParameter.h>
# include <TCanvas.h>
# include <TList.h>
# include <TFile.h>
# include <TError.h>
# include <TLatex.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TMath.h>
# include <TString.h>
# include <TStyle.h>
# include <TSystem.h>
# include <TProfile.h>
# include <TGaxis.h>
#else
class THStack;
class TH1;
class TH2;
class TCollection;
class TCanvas;
class TVirtualPad;
class TLatex;
#endif

bool fgPause = false;
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
TObject* GetObject(const TDirectory* parent, const TString& name)
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
  TObject* o = const_cast<TDirectory*>(parent)->Get(name);
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
TCollection* GetCollection(const TDirectory* parent, const TString& name)
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
void PrintCanvas(TCanvas* c, const TString& title, Float_t size=.7)
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
  if (fgPause) Pause();

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
  
  if (!h) {
    Warning("DrawInPad", "Nothing to draw in pad # %d", padNo);
    return;
  }
  h->Draw(opts);

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
void DrawTwoInPad(TVirtualPad* c, Int_t padNo, TH1* h1, TH1* h2,
		  Option_t* opts="", UShort_t flags=0x0)
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
  o.ToLower();
  TString fopts(o);
  Bool_t e3 = o.Contains("e3");
  if (e3) {
    fopts.ReplaceAll("e3", " same");
  }

  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTicks("");
  h1->GetYaxis()->SetNdivisions(0);
  h1->DrawCopy(o); // First draw with opts 
  if (e3) h1->DrawCopy(fopts);
  p->Update();

  // Make axis 
  Double_t m1 = 1.05 * h1->GetMaximum();
  TGaxis*  a1 = new TGaxis(p->GetUxmin(), p->GetUymin(), 
			   p->GetUxmin(), p->GetUymax(), 
			   0, m1, 510);
  a1->SetLineColor(h1->GetLineColor());
  a1->Draw();
  
  o.Append(" same");
  Double_t m2    = 1.1 * h2->GetMaximum();
  Double_t scale = m1 / m2;
  h2->Scale(scale);
  h2->DrawCopy(o);
  if (e3) h2->DrawCopy(fopts);

  TGaxis*  a2 = new TGaxis(p->GetUxmax(), p->GetUymin(), 
			   p->GetUxmax(), p->GetUymax(), 
			   0, m2, 510, "+L");
  a2->SetLineColor(h2->GetLineColor());
  a2->Draw();
  
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
  name->DrawLatex(name->GetX(), y, Form("%s:", sName.Data()));
  value->DrawLatex(value->GetX(), y, sValue.Data());
  y -= name->GetTextSize() + .02;
}  

  
//____________________________________________________________________
void DrawEventInspector(const TCollection* forward, TCanvas* can)
{
  Info("DrawEventInspector", "Drawing event inspector from %p", 
       forward);
  TCollection* c = GetCollection(forward, "fmdEventInspector");
  if (!c) return;

  can->cd(2);
  
  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value);

  Int_t sys, sNN, field, runNo, lowFlux, nPileUp;
  Int_t aliRev, aliBra;
  Bool_t fpVtx, v0and;
  Double_t dPileUp;
  
  GetParameter(c, "sys",           sys);
  GetParameter(c, "sNN",           sNN);
  GetParameter(c, "field",         field);
  GetParameter(c, "runNo",         runNo);
  GetParameter(c, "lowFlux",       lowFlux);
  GetParameter(c, "fpVtx",         fpVtx);
  GetParameter(c, "v0and",         v0and);
  GetParameter(c, "nPileUp",       nPileUp);
  GetParameter(c, "dPileup",       dPileUp);
  GetParameter(c, "alirootRev",    aliRev);
  GetParameter(c, "alirootBranch", aliBra);

  DrawParameter(name, value, y, "System", (sys == 1 ? "pp" : 
					   sys == 2 ? "PbPb" : 
					   sys == 3 ? "pPb" : "unknown"));
  DrawParameter(name, value, y, "#sqrt{s_{NN}}", Form("%5dGeV", sNN));
  DrawParameter(name, value, y, "L3 B field", Form("%+2dkG", field));
  DrawParameter(name, value, y, "Run #", Form("%6d", runNo));
  DrawParameter(name, value, y, "Low flux cut", Form("%6d", lowFlux));
  DrawParameter(name, value, y, "Use PWG-UD vertex", (fpVtx ? "yes" : "no"));
  DrawParameter(name, value, y, "Use V0AND for NSD", (v0and ? "yes" : "no"));
  DrawParameter(name, value, y, "Least # of pile-up vertex", 
		Form("%d", nPileUp));
  DrawParameter(name, value, y, "Least distance of pile-up vertex", 
		Form("%fcm", dPileUp));
  DrawParameter(name, value, y, "AliROOT", 
		Form("%7lu/0x%8lx", ULong_t(aliRev), ULong_t(aliBra)));

  PrintCanvas(can, "Event Inspector");

  TVirtualPad* body = can->cd(2);
  body->Divide(2,4);

  TH1*    nEventsTr    = GetH1(c, "nEventsTr");
  TH1*    nEventsTrVtx = GetH1(c, "nEventsTrVtx");
  TH1*    vertex       = GetH1(c, "vertex");
  Bool_t  mc           = (vertex != 0);
  if (nEventsTr)    nEventsTr->Rebin(2);
  if (nEventsTrVtx) nEventsTrVtx->Rebin(2);
  if (vertex) {
    // vertex->Rebin(2);
    vertex->SetFillColor(kMagenta+2);
  }
  DrawInPad(body, 1, nEventsTr);
  DrawInPad(body, 1, vertex, "same");
  DrawInPad(body, 1, nEventsTrVtx, "same"); 
  DrawInPad(body, 1, GetH1(c, "nEventsAccepted"), "same", 0x10);


  DrawInPad(body, 2, GetH2(c, "nEventsAcceptedXY"), "colz", 0x4);
  DrawInPad(body, 3, GetH1(c, "triggers"),          "hist text");
  DrawInPad(body, 4, GetH2(c, "triggerCorr"),       "colz", 0x4);
  DrawInPad(body, 5, GetH1(c, "status"),            "hist text");
  DrawInPad(body, 6, GetH1(c, "type"),              "hist text");
  DrawInPad(body, 7, GetH1(c, "cent"));
  DrawInPad(body, 8, GetH2(c, "centVsQuality"), "colz", 0x4);

  PrintCanvas(can, "EventInspector - Histograms");  

  if (!mc) return; // not MC 
  
  TH1* phiR         = GetH1(c, "phiR");
  TH1* b            = GetH1(c, "b");
  TH2* bVsNpart     = GetH2(c, "bVsParticipants");
  TH2* bVsNbin      = GetH2(c, "bVsBinary");
  TH2* bVsCent      = GetH2(c, "bVsCentrality");
  TH2* vzComparison = GetH2(c, "vzComparison");
  TH2* centVsNpart  = GetH2(c, "centralityVsParticipans");// Spelling!
  TH2* centVsNbin   = GetH2(c, "centralityVsBinary");
  
  body = can->cd(2);
  body->Divide(2,3);

  DrawInPad(body, 1, phiR);
  DrawInPad(body, 2, vzComparison, "colz", 0x4);
  DrawInPad(body, 3, b);

  TProfile* nPartB = bVsNpart->ProfileX("nPartB",1,-1,"s");
  TProfile* nBinB  = bVsNbin->ProfileX("nBinB",1,-1,"s");
  nPartB->SetMarkerColor(kBlue+2);
  nPartB->SetMarkerStyle(20);
  nPartB->SetLineColor(kBlue+2);
  nPartB->SetFillColor(kBlue-10);
  nPartB->SetFillStyle(1001);
  nPartB->SetMarkerSize(0.7);
  nBinB->SetMarkerColor(kRed+2);
  nBinB->SetMarkerStyle(21);
  nBinB->SetLineColor(kRed+2);
  nBinB->SetFillColor(kRed-10);
  nBinB->SetMarkerSize(0.7);
  nBinB->SetFillStyle(1001);

  DrawTwoInPad(body, 4, nPartB, nBinB, "e3 p", 0x10);

  DrawInPad(body, 5, bVsCent, "colz", 0x4);

  TProfile* nPartC = centVsNpart->ProfileY("nPartC",1,-1,"s");
  TProfile* nBinC  = centVsNbin->ProfileY("nBinC",1,-1,"s");
  nPartC->SetMarkerColor(kBlue+2);
  nPartC->SetMarkerStyle(20);
  nPartC->SetLineColor(kBlue+2);
  nPartC->SetFillColor(kBlue-10);
  nPartC->SetFillStyle(1001);
  nPartC->SetMarkerSize(0.7);
  nBinC->SetMarkerColor(kRed+2);
  nBinC->SetMarkerStyle(21);
  nBinC->SetLineColor(kRed+2);
  nBinC->SetFillColor(kRed-10);
  nBinC->SetMarkerSize(0.7);
  nBinC->SetFillStyle(1001);

  DrawTwoInPad(body, 6, nPartC, nBinC, "e3 p", 0x10);

  PrintCanvas(can, "EventInspector - Monte-Carlo");  
  
}

//____________________________________________________________________
void DrawSharingFilter(const TCollection* forward, TCanvas* can)
{
  Info("DrawEventInspector", "Drawing sharing filter from %p", 
       forward);
  TCollection* c = GetCollection(forward, "fmdSharingFilter");
  if (!c) return;

  TVirtualPad* body = can->cd(2);
  body->Divide(1, 3);
  body->cd(1);
  
  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value, .05);

  Bool_t angle, lowSignal, simple;
  
  GetParameter(c, "angle",     angle);
  GetParameter(c, "lowSignal", lowSignal);
  GetParameter(c, "simple",    simple);

  DrawParameter(name, value, y, "Angle correct", (angle ? "yes" : "no")); 
  DrawParameter(name, value, y, "Lower signal",  (lowSignal ? "yes" : "no"));
  DrawParameter(name, value, y, "Simple method", (simple ? "yes" : "no"));

  DrawInPad(body, 2, GetH2(c, "lowCuts"), "colz");
  DrawInPad(body, 3, GetH2(c, "highCuts"), "colz");
  
  PrintCanvas(can, "Sharing filter");

  const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
  const char** ptr   = subs;
  while (*ptr) { 
    TCollection* sc = GetCollection(c, *ptr);
    if (!sc) { ptr++; continue; }
    
    body = can->cd(2);
    body->Divide(2,3);
    DrawInPad(body, 1, GetH1(sc, "esdEloss"),       "",     0x2);
    DrawInPad(body, 1, GetH1(sc, "anaEloss"),       "same", 0x12);
    DrawInPad(body, 2, GetH1(sc, "singleEloss"),    "",     0x2);
    DrawInPad(body, 2, GetH1(sc, "doubleEloss"),    "same", 0x2);
    DrawInPad(body, 2, GetH1(sc, "tripleEloss"),    "same", 0x12);  
    DrawInPad(body, 3, GetH2(sc, "singlePerStrip"), "colz", 0x4);
    DrawInPad(body, 4, GetH1(sc, "distanceBefore"), "",     0x2);
    DrawInPad(body, 4, GetH1(sc, "distanceAfter"),  "same", 0x12);

    TH2* nB = GetH2(sc, "neighborsBefore");
    if (nB) { 
      nB->GetXaxis()->SetRangeUser(0,8); 
      nB->GetYaxis()->SetRangeUser(0,8); 
    }
    DrawInPad(body, 5, nB, "colz", 0x4);
    DrawInPad(body, 5, GetH2(sc, "neighborsAfter"), "p same", 0x4);
    DrawInPad(body, 6, GetH2(sc, "beforeAfter"),    "colz",   0x4);

    PrintCanvas(can, Form("Sharing filter - %s", *ptr));
    ptr++;
  }

  TCollection* cc = GetCollection(c, "esd_mc_comparion"); // Spelling!
  if (!cc) return; // Not MC 

  body = can->cd(2);
  body->Divide(2,3);
  DrawInPad(body, 1, GetH2(cc, "FMD1i_corr"), "colz", 0x4);
  DrawInPad(body, 3, GetH2(cc, "FMD2i_corr"), "colz", 0x4);
  DrawInPad(body, 4, GetH2(cc, "FMD2o_corr"), "colz", 0x4);
  DrawInPad(body, 5, GetH2(cc, "FMD3i_corr"), "colz", 0x4);
  DrawInPad(body, 6, GetH2(cc, "FMD3o_corr"), "colz", 0x4);

  PrintCanvas(can, "Sharing filter - MC vs Reco");

  TCollection* mc = GetCollection(c, "mcTrackDensity"); // Spelling!
  if (!mc) return; // Not MC 

  body = can->cd(2);
  body->Divide(2,3);
  DrawInPad(body, 1, GetH2(mc, "binFlow"),    "colz", 0x4);
  DrawInPad(body, 2, GetH2(mc, "binFlowEta"), "colz", 0x4);
  DrawInPad(body, 3, GetH2(mc, "binFlowPhi"), "colz", 0x4);
  DrawInPad(body, 4, GetH1(mc, "nRefs"),       "",    0x2);
  DrawInPad(body, 4, GetH1(mc, "clusterRefs"), "same");
  DrawInPad(body, 4, GetH1(mc, "clusterSize"), "same");
  DrawInPad(body, 4, GetH1(mc, "nClusters"),    "same", 0x10);
  DrawInPad(body, 5, GetH2(mc, "clusterVsRefs"),"colz", 0x4);


  PrintCanvas(can, "Sharing filter - MC");
  
}

//____________________________________________________________________
void DrawDensityCalculator(const TCollection* forward, TCanvas* can)
{
  Info("DrawDensityCalculator", "Drawing density calculator from %p", 
       forward);
  TCollection* c = GetCollection(forward, "fmdDensityCalculator");
  if (!c) return;

  TVirtualPad* body = can->cd(2);
  body->Divide(2, 2);
  body->cd(1);
  
  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value, .05);

  Int_t maxParticles, phiAcceptance, etaLumping, phiLumping;
  Bool_t method, recalcEta, recalcPhi;
  
  GetParameter(c, "maxParticle",     maxParticles);
  GetParameter(c, "phiAcceptance",   phiAcceptance);
  GetParameter(c, "etaLumping",      etaLumping);
  GetParameter(c, "phiLumping",      phiLumping);
  GetParameter(c, "method",          method);
  GetParameter(c, "recalcEta",       recalcEta);
  GetParameter(c, "recalcPhi",       recalcPhi);

  DrawParameter(name, value, y, "Method", (method ? "Poisson" : "#DeltaE")); 
  DrawParameter(name, value, y, "Recalculate #eta",(recalcEta ? "yes" : "no")); 
  DrawParameter(name, value, y, "Recalculate #phi",(recalcPhi ? "yes" : "no")); 
  DrawParameter(name, value, y, "#phi acceptance method", 
		(phiAcceptance == 1 ? "N_{ch}" : 
		 phiAcceptance == 2 ? "#DeltaE" : "none"));
  DrawParameter(name, value, y, "Region size (sector#timesstrip)", 
		Form("%2d #times %2d", phiLumping, etaLumping));

  TVirtualPad* p = body; // body->cd(2);
  // p->Divide(3,1);

  TH1* accI = GetH1(c, "accI");
  TH1* accO = GetH1(c, "accO");
  if (accI) { 
    Double_t scale = 1./accI->GetMaximum();
    accI->Scale(scale); 
    accO->Scale(scale);
    accI->SetMinimum(0); 
  }
  DrawInPad(p, 2, accI); 
  DrawInPad(p, 2, accO, "same", 0x10); 
  DrawInPad(p, 3, GetH2(c, "lowCuts"), "colz");
  DrawInPad(p, 4, GetH2(c, "maxWeights"), "colz");
  
  PrintCanvas(can, "Density calculator");

  const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
  const char** ptr   = subs;
  while (*ptr) { 
    TCollection* sc = GetCollection(c, *ptr);
    if (!sc) { ptr++; continue; }
    
    body = can->cd(2);
    body->Divide(2,3);
    
    DrawInPad(body, 1, GetH2(sc, "elossVsPoisson"), "colz",   0x4);
    DrawInPad(body, 2, GetH1(sc, "diffElossPoisson"), "",     0x2);
    DrawInPad(body, 3, GetH1(sc, "occupancy"),        "",     0x2);
    DrawInPad(body, 4, GetH1(sc, "eloss"),            "",     0x2);
    DrawInPad(body, 4, GetH1(sc, "elossUsed"),        "same", 0x12);
    TH1* phiB = GetH1(sc, "phiBefore");
    TH1* phiA = GetH1(sc, "phiAfter");
    if (phiB && phiA) { 
      phiA->Add(phiB, -1);
      phiA->Divide(phiB);
      phiA->SetTitle("#Delta#phi from Ip (x,y) correction");
      phiA->SetYTitle("(#phi_{after}-#phi_{before})/#phi_{before}");
    }
    DrawInPad(body, 5, phiA);
    DrawInPad(body, 6, GetH2(sc, "phiAcc"), "colz",   0x4);
    
    PrintCanvas(can, Form("Density calculator - %s", *ptr));
    ptr++;    
  }

  TCollection* cc = GetCollection(c, "esd_mc_comparison"); 
  if (!cc) return; // Not MC 

  body = can->cd(2);
  body->Divide(2,5);
  DrawInPad(body, 1, GetH2(cc, "FMD1I_corr_mc_esd"), "colz", 0x4);
  DrawInPad(body, 3, GetH2(cc, "FMD2I_corr_mc_esd"), "colz", 0x4);
  DrawInPad(body, 5, GetH2(cc, "FMD2O_corr_mc_esd"), "colz", 0x4);
  DrawInPad(body, 7, GetH2(cc, "FMD3O_corr_mc_esd"), "colz", 0x4);
  DrawInPad(body, 9, GetH2(cc, "FMD3I_corr_mc_esd"), "colz", 0x4);
  DrawInPad(body, 2,  GetH1(cc, "FMD1I_diff_mc_esd"), "", 0x2);
  DrawInPad(body, 4,  GetH1(cc, "FMD2I_diff_mc_esd"), "", 0x2);
  DrawInPad(body, 6,  GetH1(cc, "FMD2O_diff_mc_esd"), "", 0x2);
  DrawInPad(body, 8,  GetH1(cc, "FMD3O_diff_mc_esd"), "", 0x2);
  DrawInPad(body, 10, GetH1(cc, "FMD3I_diff_mc_esd"), "", 0x2);

  PrintCanvas(can, "Density calculator - MC vs Reco");
}

//____________________________________________________________________
void DrawCorrector(const TCollection* forward, TCanvas* can)
{
  Info("DrawCorrector", "Drawing corrector from %p", forward);
  TCollection* c = GetCollection(forward, "fmdCorrector");
  if (!c) return;
  
  TVirtualPad* body = can->cd(2);
  body->cd();
  
  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value, .05);
  
  Bool_t secondary, vertexBias, acceptance, merging;  
  GetParameter(c, "secondary",    secondary);
  GetParameter(c, "acceptance",   acceptance);
  GetParameter(c, "vertexBias",   vertexBias);
  GetParameter(c, "merging",      merging);
  
  DrawParameter(name, value, y, "Secondary corr.", secondary ? "yes" : "no");
  DrawParameter(name, value, y, "Acceptance corr.", acceptance ? "yes" : "no");
  DrawParameter(name, value, y, "Vertex bias corr.", vertexBias ? "yes" : "no");
  DrawParameter(name, value, y, "Merging eff.", merging ? "yes" : "no");
    
  PrintCanvas(can, "Corrector");

  TCollection* cc = GetCollection(c, "esd_mc_comparison"); 
  if (!cc) return; // Not MC 

  body = can->cd(2);
  body->Divide(2,3);
  DrawInPad(body, 1, GetH2(cc, "FMD1I_esd_vs_mc"), "colz", 0x0);
  DrawInPad(body, 3, GetH2(cc, "FMD2I_esd_vs_mc"), "colz", 0x0);
  DrawInPad(body, 4, GetH2(cc, "FMD2O_esd_vs_mc"), "colz", 0x0);
  DrawInPad(body, 5, GetH2(cc, "FMD3O_esd_vs_mc"), "colz", 0x0);
  DrawInPad(body, 6, GetH2(cc, "FMD3I_esd_vs_mc"), "colz", 0x0);

  PrintCanvas(can, "Corrector - MC vs Reco");
}

//____________________________________________________________________
void DrawHistCollector(const TCollection* forward, TCanvas* can)
{
  Info("DrawHistCollector", "Drawing histogram collector from %p", forward);
  TCollection* c = GetCollection(forward, "fmdHistCollector");
  if (!c) return;

  TVirtualPad* body = can->cd(2);
  body->Divide(1, 3);
  body->cd(1);

  TLatex*  name;
  TLatex*  value;
  Double_t y = .8;
  CreateTemplates(name, value, .05);
  
  Int_t nCutBins, fiducial, merge, skipRings;
  Double_t fiducialCut;
  Bool_t  bgAndHits;

  GetParameter(c, "nCutBins",       nCutBins);
  GetParameter(c, "skipRings",      skipRings);
  GetParameter(c, "bgAndHits",      bgAndHits);
  GetParameter(c, "merge",          merge);
  GetParameter(c, "fiducial",       fiducial);
  // GetParameter(c, "correctionCut",  fiducialCut);
  GetParameter(c, "fiducialCut",  fiducialCut);

  DrawParameter(name, value, y, "# of bins to cut",      Form("%d", nCutBins));
  DrawParameter(name, value, y, "Bg & hit maps stored.", bgAndHits?"yes":"no");
  DrawParameter(name, value, y, "Fiducial method.", 
		fiducial == 0 ? "cut" : "distance");
  DrawParameter(name, value, y, "Fiducial cut.", Form("%f", fiducialCut));
  DrawParameter(name, value, y, "Merge method", 
		(merge == 0 ? "straight mean" :
		 merge == 1 ? "straight mean, no zeroes" : 
		 merge == 2 ? "weighted mean" : 
		 merge == 3 ? "least error" : 
		 merge == 4 ? "sum" : "unknown"));
  TString skipped;
  if (skipRings & 0x11) skipped.Append("FMD1i ");
  if (skipRings & 0x21) skipped.Append("FMD2i ");
  if (skipRings & 0x22) skipped.Append("FMD2o ");
  if (skipRings & 0x31) skipped.Append("FMD3i ");
  if (skipRings & 0x32) skipped.Append("FMD3o ");
  DrawParameter(name, value, y, "Skipped rings", skipped);
		 
  DrawInPad(body, 2, GetH2(c, "sumRings"), "colz"); 
  DrawInPad(body, 3, GetH2(c, "coverage"), "colz");

  body->cd(1)->Modified();
  body->cd(2)->Modified();
  body->cd(3)->Modified();
  body->cd(1)->Update();
  body->cd(2)->Update();
  body->cd(3)->Update();
  PrintCanvas(can, "Histogram collector");
}
  
//____________________________________________________________________
void AddToAll(THStack* all, const THStack* stack, Int_t curr, Int_t step)
{
  if (!stack) return;

  TIter   next(stack->GetHists());
  TH1*    h = 0;
  while ((h = static_cast<TH1*>(next()))) {
    TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
    copy->SetDirectory(0);
    if (curr != step) {
      copy->SetMarkerColor(kGray);
      copy->SetLineColor(kGray);
    }
    all->Add(copy);
  }
}

//____________________________________________________________________
void DrawStep(THStack*     delta, 
	      THStack*     nchs, 
	      THStack*     prims, 
	      THStack*     rings,
	      TH1*         dndeta,
	      TVirtualPad* can, 
	      Int_t        step)
{
  THStack* all = 0;
  TH1* res = 0;
  if (step < 6) {
    all = new THStack;
    if (step != 1) AddToAll(all, delta,  step, 1);
    if (step != 2) AddToAll(all, nchs,   step, 2);
    if (step != 3) AddToAll(all, prims,  step, 3);
    if (step != 4) AddToAll(all, rings,  step, 4);

    res = static_cast<TH1*>(dndeta->Clone("dNdeta"));
    res->SetTitle("dN/d#eta");
    res->SetMarkerColor(step == 5 ? kBlack : kGray);
    res->SetLineColor(step == 5 ? kBlack : kGray);  
    res->SetDirectory(0);
    all->Add(res);
    if (step == 1) AddToAll(all, delta,  step, 1);
    if (step == 2) AddToAll(all, nchs,   step, 2);
    if (step == 3) AddToAll(all, prims,  step, 3);
    if (step == 4) AddToAll(all, rings,  step, 4);
  }

  TVirtualPad* p = can->cd(step);
  p->SetFillColor(kWhite);
  p->SetRightMargin(0.02);
  p->SetTopMargin(0.02);

  if (all) {
    all->Draw("nostack");
    all->GetHistogram()->SetXTitle("#eta");
    all->GetHistogram()->SetYTitle("signal");
    // all->GetHistogram()->GetXaxis()->SetLabelFont(132);
    // all->GetHistogram()->GetXaxis()->SetTitleFont(132);
    // all->GetHistogram()->GetYaxis()->SetLabelFont(132);
    // all->GetHistogram()->GetYaxis()->SetTitleFont(132);
  }

  TLegend* l = 0;
  if (step < 6)
    l = new TLegend(.33, .2, .53, .9);
  else 
    l = new TLegend(0.1, 0.1, .9, .9);
  l->SetFillColor(kWhite);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  TLegendEntry* e = 0;

  TIter next(delta->GetHists());
  TH1*  h = 0;
  if (step == 6) {
    while ((h = static_cast<TH1*>(next()))) {
      e = l->AddEntry("dummy", h->GetTitle(), "pl");
      e->SetMarkerStyle(20);
      e->SetMarkerColor(h->GetMarkerColor());
    }
  }
  else {
    h = static_cast<TH1*>(delta->GetHists()->At(0));
    e = l->AddEntry("dummy", delta->GetTitle(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(step != 1 ? kGray : kBlack);
    e->SetLineColor(e->GetMarkerColor());
    
    h = static_cast<TH1*>(nchs->GetHists()->At(0));
    e = l->AddEntry("dummy", nchs->GetTitle(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(step != 2 ? kGray : kBlack);
    e->SetLineColor(e->GetMarkerColor());
    
    h = static_cast<TH1*>(prims->GetHists()->At(0));
    e = l->AddEntry("dummy", prims->GetTitle(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(step != 3 ? kGray : kBlack);
    e->SetLineColor(e->GetMarkerColor());
    
    h = static_cast<TH1*>(rings->GetHists()->At(0));
    e = l->AddEntry("dummy", rings->GetTitle(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(step != 4 ? kGray : kBlack);
    e->SetLineColor(e->GetMarkerColor());

    h = res;
    e = l->AddEntry("dummy", h->GetTitle(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(step != 5 ? kGray : kBlack);
    e->SetLineColor(e->GetMarkerColor());
  }

  l->Draw("nostack");
  
  TString what;
  if (step > 0) {
    switch (step) { 
    case 1: 
      what = "After merging";
      break;
    case 2: 
      what = "After particle counting";
      break;
    case 3: 
      what = "After corrections";
      break;
    case 4: 
      what = "After normalisation";
      break;
    case 5: 
      what = "Result";
      break;
    case 6: 
      break;
    default: 
      Error("DrawStep", "Unknown step: %d (must be in 1-6)", step);
      break;
    }
  }
  TLatex* ltx = new TLatex(.97, .97, what);
  ltx->SetNDC();
  ltx->SetTextSize(.06);
  ltx->SetTextAlign(33);
  // ltx->SetTextFont(132);
  ltx->Draw();
  
  if (step < 6)
    ltx->DrawLatex(.32, .97, Form("Step %d", step));

 
}

//____________________________________________________________________
void FixStack(THStack* stack, const TString& title, Int_t marker)
{
  if (!stack) return;
  stack->SetTitle(title);
  TIter next(stack->GetHists());
  TH1*  h = 0;
  while ((h = static_cast<TH1*>(next())))  h->SetMarkerStyle(marker);
}
//____________________________________________________________________
void DrawSteps(const TCollection* forward, TCanvas* can)
{
  // MakeChapter(can, "Steps");

  THStack* deltas = GetStack(GetCollection(forward, "fmdSharingFilter"), 
			     "sums", "summed");
  THStack* nchs   = GetStack(GetCollection(forward, "fmdDensityCalculator"), 
			     "sums", "inclDensity");
  THStack* prims   = GetStack(GetCollection(forward, "fmdCorrector"), 
			      "sums", "primaryDensity");
  THStack* rings   = GetStack(GetCollection(forward, "ringResults"), "all");
  TH1*     dndeta  = GetH1(forward, "dNdeta");

  FixStack(deltas,	"#sum_{} #Delta/#Delta_{mip}",  25);
  FixStack(nchs,	"#sum_{} N_{ch,incl}", 		21);
  FixStack(prims,	"#sum_{} N_{ch,primary}",       22);
  FixStack(rings,	"dN/d#eta per ring",            23);
  
  TVirtualPad* body = can->cd(2);
  body->Divide(2,3,0,0);
  for (Int_t step = 1; step <= 6; step++) { 
    DrawStep(deltas, nchs, prims, rings, dndeta, body, step);
  }
  PrintCanvas(can, "Steps");
}


//____________________________________________________________________
void DrawResults(const TCollection* forward, 
		 const TCollection* forwardRes, TCanvas* can)
{
  // MakeChapter(can, "Results");

  TVirtualPad* body = can->cd(2);
  body->Divide(2,3, .1, 0);

  TCollection* c = GetCollection(forwardRes, "ringResults");
  if (!c) return;
  
  DrawInPad(body, 1, GetStack(c, "all"), "nostack");

  DrawInPad(body, 2, GetH1(forwardRes, "dNdeta"));
  DrawInPad(body, 3, GetH1(forwardRes, "dNdeta_"));
  DrawInPad(body, 4, GetH1(forwardRes, "norm"));
  DrawInPad(body, 4, GetH1(forwardRes, "phi"), "same", 0x10);
  DrawInPad(body, 5, GetH1(forward,    "d2Ndetadphi"), "colz");

  TCollection* mc = GetCollection(forwardRes, "mcRingResults");
  if (mc) {
    THStack* sMc = new THStack;
    sMc->SetTitle("MC dN/d#eta per Ring");
    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(mc, *ptr);
      if (!sc) { ptr++; continue; }
      
      TH1* h = GetH1(sc, "dndeta_eta");
      if (!h) { ptr++; continue; }
      sMc->Add(h);
      ptr++;
    }
    DrawInPad(body, 1, sMc, "nostack same", 0x10);
  }

  
  PrintCanvas(can, "Results");

}
  
  
//____________________________________________________________________
void DrawAODSummary(const char* fname="forward.root",
		    UShort_t what=0x7F)
{
  // --- Open the file -----------------------------------------------
  TString filename(fname);
  TFile* file = TFile::Open(filename.Data(), "READ");
  if (!file) { 
    Error("DrawAODSummary", "Failed to open \"%s\"", filename.Data());
    return;
  }

  // --- Get top-level collection ------------------------------------
  TCollection* forward = GetCollection(file, "Forward");
  if (!forward) return;

  // --- Make our canvas ---------------------------------------------
  TString pdfName(filename);
  pdfName.ReplaceAll(".root", ".pdf");

  TCanvas* c = CreateCanvas(pdfName);

  fgPause = what & 0x80;

  // --- Do each sub-algorithm ---------------------------------------
  if (what & 0x01) DrawEventInspector(forward,c);
  if (what & 0x02) DrawSharingFilter(forward,c);
  if (what & 0x04) DrawDensityCalculator(forward, c);
  if (what & 0x08) DrawCorrector(forward, c);
  if (what & 0x10) DrawHistCollector(forward, c);
  
  // --- Do the results ----------------------------------------------
  TCollection* forwardRes = GetCollection(file, "ForwardResults");
  if (!forwardRes) return;

  if (what & 0x20) DrawSteps(forwardRes, c);
  if (what & 0x40) DrawResults(forward, forwardRes, c);
  
  CloseCanvas(c);
}
  

  
