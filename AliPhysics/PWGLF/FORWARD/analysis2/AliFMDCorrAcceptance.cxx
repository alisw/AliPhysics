//
// This class contains the acceptance correction due to dead channels 
// 
//
#include "AliFMDCorrAcceptance.h"
#include <TBrowser.h>
#include <TH2D.h>
#include <AliLog.h>
#include <AliForwardUtil.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMath.h>
#include <THStack.h>
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrAcceptance::AliFMDCorrAcceptance()
  : fRingArray(), 
    fCache(0),
    fVertexAxis(0,0,0),
    fHasOverflow(false)
{
  // 
  // Default constructor 
  //
  fRingArray.SetOwner(kTRUE);
  fRingArray.SetName("rings");
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
  
}
//____________________________________________________________________
AliFMDCorrAcceptance::AliFMDCorrAcceptance(const 
					       AliFMDCorrAcceptance& o)
  : TObject(o), 
    fRingArray(o.fRingArray), 
    fCache(o.fCache),
    fVertexAxis(o.fVertexAxis.GetNbins(), o.fVertexAxis.GetXmin(), 
		o.fVertexAxis.GetXmax()),
    fHasOverflow(o.fHasOverflow)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  fVertexAxis.SetName("vtxAxis");
  fVertexAxis.SetTitle("v_{z} [cm]");
}
//____________________________________________________________________
AliFMDCorrAcceptance::~AliFMDCorrAcceptance()
{
  //
  // Destructor 
  // 
  //
  fRingArray.Clear();
  if (fCache) fCache->Clear();
}
//____________________________________________________________________
AliFMDCorrAcceptance&
AliFMDCorrAcceptance::operator=(const AliFMDCorrAcceptance& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this;

  fRingArray        = o.fRingArray;
  fCache            = o.fCache;
  fHasOverflow      = o.fHasOverflow;
  SetVertexAxis(o.fVertexAxis);

  return *this;
}
//____________________________________________________________________
TH2D*
AliFMDCorrAcceptance::GetCorrection(UShort_t d, Char_t r, Double_t v) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetCorrection(d, r, UShort_t(b));
}
//____________________________________________________________________
TH2D*
AliFMDCorrAcceptance::GetCorrection(UShort_t d, Char_t r, UShort_t b) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  return static_cast<TH2D*>(GetObject(fRingArray, d, r, b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrAcceptance::GetPhiAcceptance(UShort_t d, Char_t r, Double_t v) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    v  Primary interaction point @f$z@f$ coordinate
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0) return 0;
  return GetPhiAcceptance(d, r, UShort_t(b));
}
//____________________________________________________________________
TH1D*
AliFMDCorrAcceptance::GetPhiAcceptance(UShort_t d, Char_t r, UShort_t b) const
{
  // 
  // Get the acceptance correction @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  if (!fHasOverflow) return 0;
  if (!fCache) FillCache();
  return static_cast<TH1D*>(GetObject(*fCache, d, r, b));
}
  
//____________________________________________________________________
Int_t
AliFMDCorrAcceptance::FindVertexBin(Double_t v) const
{
  // 
  // Find the vertex bin that corresponds to the passed vertex 
  // 
  // Parameters:
  //    vertex The interaction points @f$z@f$-coordinate 
  // 
  // Return:
  //    Vertex bin in @f$[1,N_{\mbox{vertex}}]@f$ or negative if 
  // out of range 
  //
  if (fVertexAxis.GetNbins() <= 0) { 
    AliWarning("No vertex array defined");
    return 0;
  }
  Int_t bin = const_cast<TAxis&>(fVertexAxis).FindBin(v);
  if (bin <= 0 || bin > fVertexAxis.GetNbins()) { 
    AliWarning(Form("vertex %+8.4f out of range [%+8.4f,%+8.4f]",
		    v, fVertexAxis.GetXmin(), fVertexAxis.GetXmax()));
    return 0;
  }
  return bin;
}
//____________________________________________________________________
Int_t
AliFMDCorrAcceptance::GetRingIndex(UShort_t d, Char_t r) const
{
  // Get the index corresponding to the given ring 
  // 
  // Parameters:
  //    d Detector
  //    r Ring 
  // 
  // Return:
  //    Index (0 based) or negative in case of errors
  //
  switch (d) {
  case 1:  return 0;
  case 2:  return (r == 'I' || r == 'i' ? 1 : 2); break;  
  case 3:  return (r == 'I' || r == 'i' ? 3 : 4); break;  

  }
  AliWarning(Form("Index for FMD%d%c not found", d, r));
  return -1;
}
//____________________________________________________________________
TObject*
AliFMDCorrAcceptance::GetObject(const TObjArray& m, UShort_t d, 
				Char_t r, UShort_t b) const
{
  // 
  // Get the object @f$ a_{r,v}@f$ 
  // 
  // Parameters:
  //    m  Mother list
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    b  Bin corresponding to the primary interaction point 
  //           @f$z@f$ coordinate (1 based)
  // 
  // Return:
  //    The correction @f$ a_{r,v}@f$ 
  //
  TObjArray* ringArray = GetRingArray(m, d, r);
  if (!ringArray) return 0;

  if (b <= 0 || b > ringArray->GetEntriesFast()) {
    AliWarning(Form("vertex bin %d out of range [1,%d]", 
		    b, ringArray->GetEntriesFast()));
    return 0;
  }

  TObject* o = ringArray->At(b-1);
  if (o) return o;

  AliWarning(Form("No dead channels map found for FMD%d%c in vertex bin %d",
		    d,r,b));
  return 0;
}
//____________________________________________________________________
TObjArray*
AliFMDCorrAcceptance::GetRingArray(const TObjArray& m, 
				   UShort_t d, Char_t r) const
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or null in case of problems
  //
  Int_t idx = GetRingIndex(d,r);
  if (idx < 0) return 0;
  
  TObject* o = m.At(idx);
  if (!o) { 
    AliWarning(Form("No array found for FMD%d%c", d, r));
    return 0;
  }

  return static_cast<TObjArray*>(o);
}
//____________________________________________________________________
TObjArray*
AliFMDCorrAcceptance::GetOrMakeRingArray(TObjArray& m, 
					 UShort_t d, Char_t r) const
{
  // 
  // Get the ring array corresponding to the specified ring
  // 
  // Parameters:
  //    d Detector 
  //    r Ring 
  // 
  // Return:
  //    Pointer to ring array, or newly created container 
  //
  Int_t idx = GetRingIndex(d,r);
  if (idx < 0) return 0;
  
  TObject* o = m.At(idx);
  if (!o) { 
    TObjArray* a = new TObjArray(fVertexAxis.GetNbins());
    a->SetName(Form("FMD%d%c", d, r));
    a->SetOwner(kTRUE);
    m.AddAtAndExpand(a, idx);
    return a;
  }

  return static_cast<TObjArray*>(m.At(idx));
}

//____________________________________________________________________
Bool_t
AliFMDCorrAcceptance::SetCorrection(UShort_t d, Char_t r, 
				    UShort_t b, TH2D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$ 
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    b    Bin corresponding to the primary interaction point 
  //             @f$z@f$ coordinate  (1 based)
  //    h    @f$ a_{r,v}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  TObjArray* ringArray = GetOrMakeRingArray(fRingArray, d, r);
  if (!ringArray) return false;
  
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex bin %3d out of range [1,%3d]", 
		    b, fVertexAxis.GetNbins()));
    return false;
  }
  h->SetName(Form("FMD%d%c_vtxbin%03d", d, r, b));
  h->SetTitle(Form("Acceptance correction for FMD%d%c "
		   "in vertex bin %d [%+8.4f,%+8.4f]", 
		   d, r, b, fVertexAxis.GetBinLowEdge(b), 
		   fVertexAxis.GetBinUpEdge(b)));
  h->SetXTitle("#eta");
  h->SetYTitle("N_{strips,OK}/N_{strips}");
  h->SetFillStyle(3001);
  h->SetDirectory(0);
  h->SetStats(0);
  ringArray->AddAtAndExpand(h, b-1);
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDCorrAcceptance::SetCorrection(UShort_t d, Char_t r, 
				    Double_t v, TH2D*  h) 
{
  // 
  // Set the acceptance correction @f$ a_{r,v}(\eta)@f$.
  // Note, that the object takes ownership of the passed pointer.
  // 
  // Parameters:
  //    d    Detector number (1-3)
  //    r    Ring identifier (I or O)
  //    v    Primary interaction point @f$z@f$ coordinate  
  //    h    @f$ a_{r,v}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  Int_t b = FindVertexBin(v);
  if (b <= 0 || b > fVertexAxis.GetNbins()) { 
    AliWarning(Form("Vertex %+8.4f out of range [%+8.4f,%+8.4f]", 
		    v, fVertexAxis.GetXmin(), fVertexAxis.GetXmax()));
    return false;
  }
  return SetCorrection(d, r, UShort_t(b), h);
}

//____________________________________________________________________
void
AliFMDCorrAcceptance::FillCache() const
{
  if (fCache) return;

  fCache = new TObjArray;
  fCache->SetOwner(kTRUE);
  fCache->SetName("cache");

  Int_t nV = fVertexAxis.GetNbins();
  for (UShort_t v = 1; v <= nV; v++) {
    for(UShort_t d = 1; d <= 3;d++) { 
      UShort_t nR = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nR; q++) { 
	Char_t   r  = (q == 0 ? 'I' : 'O');

	TObjArray* a = GetOrMakeRingArray(*fCache, d, r);

	TH2D* corr = GetCorrection(d, r, v);
	if (!corr) continue;

	Int_t nY = corr->GetNbinsY();
	TH1D* h  = corr->ProjectionX("tmp", nY+1, nY+1, "");
	h->SetName(Form("FMD%d%c_vtxbin%03d", d, r, v));
	h->SetTitle(Form("#phi acceptance correction for FMD%d%c "
			   "in vertex bin %d [%+8.4f,%+8.4f]", 
			   d, r, v, fVertexAxis.GetBinLowEdge(v), 
			   fVertexAxis.GetBinUpEdge(v)));
	h->SetXTitle("#eta");
	h->SetYTitle("N_{strips}/N_{strips,OK}");
	h->SetFillStyle(3001);
	h->SetDirectory(0);
	h->SetStats(0);
	a->AddAtAndExpand(h,v-1);

	if (fHasOverflow) continue;

	// Construct the overflow bin from 
	Int_t nX = corr->GetNbinsX();
	for (Int_t eta = 1; eta <= nX; eta++) { 
	  Double_t sum = 0;
	  for (Int_t phi = 1; phi <= nY; phi++) 
	    sum += corr->GetBinContent(eta, phi);
	  if (nY <= 0) continue;
	  h->SetBinContent(eta, nY/sum);
	} // for eta 
      } // for q 
    } // for d 
  } // for v 
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::Browse(TBrowser* b)
{
  // 
  // Browse this object in the browser
  // 
  // Parameters:
  //    b 
  //
  b->Add(&fRingArray);
  if (fCache) b->Add(fCache);
  b->Add(&fVertexAxis);
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::Print(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  std::cout << "Acceptance correction due to dead channels" << std::endl;
  fRingArray.Print(option);
  fVertexAxis.Print(option);
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::ls(Option_t* option) const
{
  // 
  // Print this object 
  // 
  // Parameters:
  //    option 
  //  
  TObject::ls(option);
  gROOT->IncreaseDirLevel();
  fVertexAxis.ls(option);
  fRingArray.ls(option);
  gROOT->DecreaseDirLevel();
}

#if 0
namespace {
  void ClearCanvas(TVirtualPad* c)
  {
    c->SetLeftMargin(.1);
    c->SetRightMargin(.05);
    c->SetBottomMargin(.1);
    c->SetTopMargin(.05);
    c->Clear();
  }
}

//____________________________________________________________________
void
AliFMDCorrAcceptance::SaveAs(const Char_t* filename, Option_t* option) const
{
  // 
  // Override to allow saving to a PDF 
  // 
  TString fileName(filename);
  if (!fileName.EndsWith(".pdf")) {
    TObject::SaveAs(fileName, option);
    return;
  }
  
  TVirtualPad* c = new TCanvas(filename, GetTitle(), 800/TMath::Sqrt(2), 800);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->Print(Form("%s[", filename));

  //__________________________________________________________________
  // Create a title page 
  TLatex* ll = new TLatex(.5,.8, filename);
  ll->SetTextAlign(22);
  ll->SetTextSize(0.03);
  ll->SetNDC();
  ll->Draw();

  TLatex* l = new TLatex(.5,.8, filename);
  l->SetNDC();
  l->SetTextSize(0.03);
  l->SetTextFont(132);
  l->SetTextAlign(12);
  l->DrawLatex(0.2, 0.70, "Acceptance due to dead channels");
  l->SetTextAlign(22);
  l->DrawLatex(0.5, 0.60, "c_{v,r}(#eta,#phi)=#frac{"
	       "#sum active strips#in(#eta,#phi)}{"
	       "#sum strips#in(#eta,#phi)}");
  
  c->Print(filename, "Title:Title page");
  
  //__________________________________________________________________
  // Draw all corrections
  const TAxis& vtxAxis = GetVertexAxis();
  Int_t        nVtx    = vtxAxis.GetNbins();

  // --- Loop over detectors -----------------------------------------
  for (UShort_t d = 1; d <= 3; d++) {
    UShort_t     nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');

      ClearCanvas(c);
      c->Divide(2, (nVtx+1)/2);
      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p = c->cd(v);
	p->SetFillColor(kWhite);
      
	TH2* h2 = GetCorrection(d, r, v);
	if (!h2) { 
	  Warning("DrawCorrAcc", "No correction for r=%c, v=%d", r, v);
	  continue;
	}
	h2->Draw(option);
      }
      c->Print(filename, Form("Title:FMD%d%c", d, r));
    }
  }
  if (HasOverflow()){
    const_cast<AliFMDCorrAcceptance*>(this)->Draw(Form("%s phi", option));
    c->Print(filename, "Title:Phi Acceptance");
  }
  const_cast<AliFMDCorrAcceptance*>(this)->Draw(option);
  c->Print(filename, "Title:Summary");
  c->Print(Form("%s]", filename));
}
//____________________________________________________________________
void
AliFMDCorrAcceptance::Draw(Option_t* option)
{
  //
  // Draw this object 
  // 
  // Parameters: 
  //   option 
  // 
  TString opt(option);
  opt.ToLower();
  Bool_t over = opt.Contains("phi");
  opt.ReplaceAll("phi", "");

  TVirtualPad* c = gPad;
  if (!c) c = new TCanvas(GetName(), GetTitle());
  
  const TAxis& vtxAxis = fVertexAxis;
  Int_t        nVtx    = vtxAxis.GetNbins();
  Int_t        ipad    = 0;
  c->SetLeftMargin(.1);
  c->SetRightMargin(.05);
  c->SetBottomMargin(.1);
  c->SetTopMargin(.05);
  c->Clear();
  c->Divide((nVtx+2)/3, 3, 0, 0);

  // Draw all corrections
  for (UShort_t v = 1; v <= nVtx; v++) { 
    ipad++;
    if (ipad == 1 || ipad == 12) ipad++;

    TVirtualPad* p = c->cd(ipad);
    p->SetFillColor(kWhite);
        
    THStack* stack = new THStack(Form("vtx%02d", v),
				 Form("%+5.1f<v_{z}<%+5.1f",
				      vtxAxis.GetBinLowEdge(v),
				      vtxAxis.GetBinUpEdge(v)));
    for (UShort_t d = 1; d <= 3; d++) {
      UShort_t     nQ = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nQ; q++) { 
	Char_t r = (q == 0 ? 'I' : 'O');
	
	if (over) { 
	  TH1* hp = GetPhiAcceptance(d, r, v);
	  if (!hp) { 
	    Error("", "No phi acceptance at v=%d", v-1);
	    continue;
	  }
	  hp->SetDirectory(0);
	  hp->SetMarkerColor(AliForwardUtil::RingColor(d, r));
	  hp->SetLineColor(AliForwardUtil::RingColor(d, r));
	  hp->SetFillColor(AliForwardUtil::RingColor(d, r));
	  hp->SetFillStyle(3001);
	  // Info("", "Adding phi acceptance plot %d", int(hp->GetEntries()));
	  stack->Add(hp);
	  continue;
	}
	  
	TH2* h1 = GetCorrection(d, r, v);
	if (!h1) { 
	  Warning("Draw", "No correction for r=%c, v=%d", r, v);
	  continue;
	}
	Int_t nY = h1->GetNbinsY();
	TH1* hh = h1->ProjectionX(Form("FMD%d%c", d, r), 1, nY);
	hh->Scale(1. / nY);
	hh->SetDirectory(0);
	hh->SetMarkerColor(AliForwardUtil::RingColor(d, r));
	hh->SetLineColor(AliForwardUtil::RingColor(d, r));
	hh->SetFillColor(AliForwardUtil::RingColor(d, r));
	hh->SetFillStyle(3004);

	stack->Add(hh);
      }
    }
    stack->SetMaximum(1.2);
    stack->Draw(Form("nostack %s", opt.Data()));
  }
}
#endif

//____________________________________________________________________
//
// EOF
//
