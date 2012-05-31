#ifdef SPECTRA_BUILD
#include <TAxis.h>
#include <TNamed.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TBrowser.h>
#include <TFile.h>

namespace Spectra 
{
  //__________________________________________________________________
  /**
   * An Element 
   * 
   */
  struct Element : public TNamed
  {
    /** 
     * Destructor 
     */
    virtual ~Element() 
    { 
      if (fFull) delete fFull;
      if (fFull) delete fCut; 
      fChildren.Delete(); 
    } 
    /** 
     * Draw this element.  Draws the full histogram and the selected
     * range on top. 
     * 
     * @param option Drawing option  
     * @param l      Low cut
     * @param h      High cut
     */
    virtual void DrawIt(Option_t* option="", Double_t l=-1, Double_t h=-1) //*MENU*
    {
      if (!fFull) return;
      if (fFull->GetMinimum() < 0) gPad->SetLogy(kFALSE);
      fFull->Draw(option);
      if (!fCut || l == h) return;
      Double_t lx = fFull->GetXaxis()->GetXmin();
      Double_t hx = fFull->GetXaxis()->GetXmax();
      Double_t rr = (hx-lx);
      Double_t ll = rr * l + lx;
      Double_t hh = rr * h + lx;
      for (Int_t i = 1; i <= fFull->GetNbinsX(); i++) { 
	if (fFull->GetBinCenter(i) <= ll || 
	    fFull->GetBinCenter(i) >  hh) { 
	  fCut->SetBinContent(i, 0); 
	  continue;
	}
	fCut->SetBinContent(i, fFull->GetBinContent(i));
      }
      fCut->Draw(Form("%s same", option));
      gPad->Modified();
      gPad->Update();
      gPad->cd();
    }//*MENU*
    /** 
     * Draw 
     * 
     * @param option 
     */
    virtual void Draw(Option_t* option="") 
    {
      if (!fFull) return;
      if (fFull->GetMinimum() < 0) gPad->SetLogy(kFALSE);
      fFull->Draw(option);
      gPad->Modified();
      gPad->Update();
      gPad->cd();
    } //*MENU*
    /** 
     * Get the top-level node
     * 
     * @return Top-level node 
     */
    virtual Element& GetParent() { return *fParent; }
    /** 
     * Make the histograms 
     * 
     * @param axis Axis to use 
     */
    virtual void MakeHistograms(const TAxis& axis)
    {
      if (fFull) return;
      if (axis.GetNbins() <= 1) return;
      if (axis.IsVariableBinSize()) 
	fFull = new TH1F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 axis.GetXbins()->fN, axis.GetXbins()->fArray);
      else 
	fFull = new TH1F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
      fCut = static_cast<TH1F*>(fFull->Clone(Form("c_%s", GetName())));
      fCut->SetTitle(Form("%s restricted spectra", GetName()));
      fFull->SetDirectory(0);
      fFull->SetFillColor(kRed);
      fFull->SetFillStyle(3001);
      fCut->SetDirectory(0);
      fCut->SetFillColor(kBlue);
      fCut->SetFillStyle(3001);
    }
    /** 
     * Get Children
     *
     * @return Array of children
     */
    TObjArray& Children() { return fChildren; }
    /** 
     * Get Children
     *
     * @return Array of children
     */
    const TObjArray& Children() const { return fChildren; }
    /** 
     * Write to a directory 
     * 
     * @param d Directory to write to
     */
    void WriteOut(TDirectory* d) 
    {
      TDirectory* dd = d->mkdir(GetName());
      dd->cd();
      if (fFull) fFull->Write();
      if (fCut)  fCut->Write();
    
      TIter next(&fChildren); 
      Element* e = 0;
      while ((e = static_cast<Element*>(next()))) 
	e->WriteOut(dd);
      d->cd();
    }
    /** 
     * This is a folder
     * 
     * @return Always true 
     */
    virtual Bool_t IsFolder() const { return kTRUE; }
    /** 
     * Browse this element 
     * 
     * @param b Brower to use 
     */
    virtual void Browse(TBrowser* b) 
    {
      if (fFull) b->Add(fFull);
      if (fCut)  b->Add(fCut);
      
      TIter next(&fChildren);
      TObject* o = 0;
      while ((o = next())) b->Add(o);
      // if (fChildren.GetEntriesFast() > 0) b->Add(&fChildren);
    }
  protected:
    Element() : TNamed(), fFull(0), fCut(0), fParent(0), fChildren(0) {}
    /** 
     * Constructor 
     * 
     * @param name   Name 
     * @param title  Title
     */
    Element(const char* name, const char* title,  Element& parent) 
      : TNamed(name, title), fFull(0), fCut(0),
	fParent(&parent), fChildren(0)
    {
      fChildren.SetOwner();
      fChildren.SetName("children");
    }
    /** 
     * Fill value 
     * 
     * @param v Value to fill
     */
    void DoFill(Double_t v) { if (fFull) fFull->Fill(v); }
    /** 
     * Get a child or null 
     * 
     * @param id Id of child
     * 
     * @return Pointer to child, or null
     */
    Element* GetChild(UShort_t id) const 
    {
      if (id >= fChildren.GetEntriesFast()) return 0;
      return static_cast<Element*>(fChildren.At(id));
    }
    TH1*      fFull; // Full histogram
    TH1*      fCut;  // Selected data 
    Element*  fParent;
    TObjArray fChildren; 

    ClassDef(Element,1); 
  };

  //__________________________________________________________________
  struct Detector;
  /**
   * Top-level node 
   * 
   */
  struct Top : public Element
  {
    /** 
     * Constructor 
     */
    Top() : Element("top", "Top", *this), fAxis(1,0,1) {}
    /** 
     * Get or add a detector
     * 
     * @param d Detector (1-3)
     * 
     * @return Reference to detector 
     */
    Detector& GetOrAdd(UShort_t d);
    /** 
     * Fill in values 
     * 
     * @param d  Detector
     * @param r  Ring 
     * @param s  Sector
     * @param t  Strip 
     * @param v  Value 
     */
    void Fill(UShort_t d, Char_t r, UShort_t s, UShort_t t, Double_t v);
    /** 
     * Get the axis to use 
     *
     * @return Reference to axis 
     */
    const TAxis& GetAxis() const { return fAxis; }
    /** 
     * Set the axis to use 
     * 
     * @param a Axis used 
     */
    void SetAxis(const TAxis& a) 
    { 
      fAxis.Set(a.GetNbins(),a.GetXmin(),a.GetXmax()); 
    }
    void SetAxis(Int_t n, Double_t l, Double_t h) { fAxis.Set(n, l, h); }
  protected:
    TAxis fAxis; // The axis 

    ClassDef(Top,1); 
  };

  //__________________________________________________________________
  struct Ring;
  /**
   * Detector
   * 
   */
  struct Detector : public Element
  {
    Detector() : Element(), fId(0) {}
    /** 
     * Constrictor 
     * 
     * @param d    Id
     * @param top  PArent node
     * @param a    Axis 
     */
    Detector(UShort_t d, Top& top)
      : Element(Form("FMD%d", d), Form("FMD%d",d), top), fId(d)
    {}
    /** 
     * Get ID
     * 
     * @return The ID
     */
    UShort_t Id() const { return fId; } 
    /** 
     * Get or add a sub element
     * 
     * @param id Id of sub-element
     * 
     * @return Sub element
     */
    Ring& GetOrAdd(Char_t id);
    /** 
     * Fill in value 
     * 
     * @param r  Ring
     * @param s  Sector
     * @param t  Strip
     * @param v  Value
     */
    void Fill(Char_t r, UShort_t s, UShort_t t, Double_t v);

    /** 
     * Get Top
     * 
     * @return top
     */
    Top& M() { return static_cast<Top&>(*fParent); }
    /** 
     * Get Top
     * 
     * @return top
     */
    const Top& M() const { return static_cast<Top&>(*fParent); }
  protected:
    UShort_t R2Id(Char_t r) { return (r == 'I' || r == 'i') ? 0 : 1; }
    UShort_t fId;

    ClassDef(Detector,1); 
  };

  //__________________________________________________________________
  struct Sector;
  /**
   * A ring
   * 
   */
  struct Ring : public Element
  {
    Ring() : Element(), fId('\0') {}
    Ring(UShort_t r, Detector& d)
      : Element(Form("%s%c", d.GetName(), r), 
		Form("%s%c", d.GetName(), r), d), fId(r)
    {}
    /** 
     * Get ID
     * 
     * @return The ID
     */
    Char_t Id() const { return fId; } 
    /** 
     * Get or add a sub element
     * 
     * @param id Id of sub-element
     * 
     * @return Sub element
     */
    Sector& GetOrAdd(UShort_t id);
    /** 
     * Fill in value 
     * 
     * @param s  Sector
     * @param t  Strip
     * @param v  Value
     */
    void Fill(UShort_t s, UShort_t t, Double_t v);

    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    Detector& D() { return static_cast<Detector&>(*fParent); }
    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    const Detector& D() const { return static_cast<Detector&>(*fParent); }
    /** 
     * Get Top
     * 
     * @return top
     */
    Top& M() { return D().M(); }
    /** 
     * Get Top
     * 
     * @return top
     */
    const Top& M() const { return D().M(); }

    /** 
     * Draw 
     * 
     * @param option 
     */
    virtual void Draw(Option_t* option="lego2z") 
    {
      if (!fFull) return;
      gPad->SetLogy(fFull->GetMinimum() > 0);
      fFull->Draw(option);
      gPad->cd();
    } //*MENU*
    void MakeHistograms(const TAxis& axis)
    {
      if (fFull) return;
      if (axis.GetNbins() <= 1) return;
      Int_t nSec = (fId == 'I' || fId == 'i' ? 20 : 40);
      if (axis.IsVariableBinSize()) 
	fFull = new TH2F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 nSec, -.5, nSec-.5,
			 axis.GetXbins()->fN, axis.GetXbins()->fArray); 
      else 
	fFull = new TH2F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 nSec, -.5, nSec-.5,
			 axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
    }
  protected:
    Char_t fId;

    ClassDef(Ring,1); 
  };


  //__________________________________________________________________
  struct Strip;
  /**
   * A sector
   * 
   */
  struct Sector : public Element
  {
    Sector() : Element(), fId(1024) {}
    Sector(UShort_t s, Ring& r)
      : Element(Form("%s[%02d]", r.GetName(), s),
		Form("%s[%02d]", r.GetName(), s), r), fId(s)
    {
    }
    /** 
     * Get ID
     * 
     * @return The ID
     */
    UShort_t Id() const { return fId; } 
    /** 
     * Get or add a sub element
     * 
     * @param id Id of sub-element
     * 
     * @return Sub element
     */
    Strip& GetOrAdd(UShort_t id);
    /** 
     * Fill in value 
     * 
     * @param t  Strip
     * @param v  Value
     */
    void Fill(UShort_t t, Double_t v); 

    /** 
     * Draw 
     * 
     * @param option 
     */
    virtual void Draw(Option_t* option="lego2z") 
    {
      if (!fFull) return;
      gPad->SetLogy(fFull->GetMinimum() > 0);
      fFull->Draw(option);
      gPad->cd();
    } //*MENU*
    void MakeHistograms(const TAxis& axis)
    {
      if (fFull) return;
      if (axis.GetNbins() <= 1) return;
      Int_t nStr = (R().Id() == 'I' || fId == 'i' ? 512 : 256);
      if (axis.IsVariableBinSize()) 
	fFull = new TH2F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 nStr, -.5, nStr-.5,
			 axis.GetXbins()->fN, axis.GetXbins()->fArray); 
      else 
	fFull = new TH2F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
			 nStr, -.5, nStr-.5,
			 axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
    }

    /** 
     * Get parent ring
     * 
     * @return Parent ring
     */
    Ring& R() { return static_cast<Ring&>(*fParent); }
    /** 
     * Get parent ring
     * 
     * @return Parent ring
     */
    const Ring& R() const { return static_cast<Ring&>(*fParent); }
    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    Detector& D() { return R().D(); }
    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    const Detector& D() const { return R().D(); }
    /** 
     * Get Top
     * 
     * @return top
     */
    Top& M() { return D().M(); }
    /** 
     * Get Top
     * 
     * @return top
     */
    const Top& M() const { return D().M(); }
  protected:
    UShort_t fId;

    ClassDef(Sector,1); 
  };

  //__________________________________________________________________
  /**
   * A stripo
   * 
   */
  struct Strip : public Element
  {
    Strip() : Element(), fId(1024) {}
    Strip(UShort_t t, Sector& s)
      : Element(Form("%s[%03d]", s.GetName(), t),
		Form("%s[%03d]", s.GetName(), t), s), fId(t) {}
    /** 
     * Get ID
     * 
     * @return The ID
     */
    UShort_t Id() const { return fId; } 
    /** 
     * Fill in value 
     * 
     * @param v  Value
     */
    void Fill(Double_t v) 
    { 
      // Info("Fill", "%s: Filling with %f", GetName(), v);
      DoFill(v); 
    }
    Bool_t IsFolder() const { return 0; } 
    void Browse(TBrowser* b) { Draw(b->GetDrawOption()); }
    /** 
     * Get parent sector
     * 
     * @return Parent sector
     */
    Sector& S() { return static_cast<Sector&>(*fParent); }
    /** 
     * Get parent sector
     * 
     * @return Parent sector
     */
    const Sector& S() const { return static_cast<Sector&>(*fParent); }
    /** 
     * Get parent ring
     * 
     * @return Parent ring
     */
    Ring& R() { return S().R(); }
    /** 
     * Get parent ring
     * 
     * @return Parent ring
     */
    const Ring& R() const { return S().R(); }
    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    Detector& D() { return R().D(); }
    /** 
     * Get parent detector
     * 
     * @return Parent detector
     */
    const Detector& D() const { return R().D(); }
    /** 
     * Get Top
     * 
     * @return top
     */
    Top& M() { return D().M(); }
    /** 
     * Get Top
     * 
     * @return top
     */
    const Top& M() const { return D().M(); }
  protected:
    UShort_t fId;

    ClassDef(Strip,1); 
  };

  //==================================================================
  Strip& Sector::GetOrAdd(UShort_t id)
  {
    Strip* t = static_cast<Strip*>(GetChild(id));
    if (!t) {
      t = new Strip(id, *this);
      t->MakeHistograms(M().GetAxis());
      // Info("GetOrAdd", "%s: Adding strip @ %d", GetName(), id);
      fChildren.AddAtAndExpand(t, id);
    }
    return *t;  
  }
  //__________________________________________________________________
  Sector& Ring::GetOrAdd(UShort_t id) 
  { 
    Sector* s = static_cast<Sector*>(GetChild(id));
    if (!s) {
      s = new Sector(id, *this);
      s->MakeHistograms(M().GetAxis());
      // Info("GetOrAdd", "%s: Adding sector @ %d", GetName(), id);
      fChildren.AddAtAndExpand(s, id);
    }
    return *s;
  }
  //__________________________________________________________________
  Ring& Detector::GetOrAdd(Char_t id) 
  { 
    UShort_t idd = R2Id(id);
    Ring*    r   = static_cast<Ring*>(GetChild(idd));
    if (!r) {
      r = new Ring(id, *this);
      r->MakeHistograms(M().GetAxis());
      // Info("GetOrAdd", "%s: Adding ring @ %d", GetName(), idd);
      fChildren.AddAtAndExpand(r, idd);
    }
    return *r;
  }
  //__________________________________________________________________
  Detector& Top::GetOrAdd(UShort_t id) 
  { 
    Int_t     idx = id - 1;
    Detector* d   = static_cast<Detector*>(fChildren.At(idx));
    if (!d) { 
      d = new Detector(id, *this);
      d->MakeHistograms(fAxis);
      // Info("GetOrAdd", "%s: Adding detector @ %d", GetName(), id);
      fChildren.AddAtAndExpand(d, idx);
    }
    return *d;
  }
  //__________________________________________________________________
  void Top::Fill(UShort_t d, Char_t r, UShort_t s, UShort_t t, Double_t v)
  {
    DoFill(v);
    Detector& sub = GetOrAdd(d);
    //Info("Fill", "%s: Filling %d,%c,%d,%d with %f", GetName(), d, r, s, t, v);
    sub.Fill(r, s, t, v);
  }
  //__________________________________________________________________
  void Detector::Fill(Char_t r, UShort_t s, UShort_t t, Double_t v)
  {
    DoFill(v);
    Ring& sub = GetOrAdd(r);
    //Info("Fill", "%s: Filling %c,%d,%d with %f", GetName(), r, s, t, v);
    sub.Fill(s, t, v);
  }
  //__________________________________________________________________
  void Ring::Fill(UShort_t s, UShort_t t, Double_t v)
  {
    fFull->Fill(s, v);
    Sector& sub = GetOrAdd(s);
    //Info("Fill", "%s: Filling %d,%d with %f", GetName(), s, t, v);
    sub.Fill(t, v);
  }
  //__________________________________________________________________
  void Sector::Fill(UShort_t t, Double_t v)
  {
    fFull->Fill(t, v);
    Strip& sub = GetOrAdd(t);
    //Info("Fill", "%s: Filling %d with %f", GetName(), t, v);
    sub.Fill(v);
  }
}
//======================================================================
#include <AliFMDDigit.h>
#include <AliFMDSDigit.h>
#include <AliFMDRecPoint.h>
#include <AliFMDHit.h>
#include <AliESDFMD.h>
#include <AliFMDInput.h>
#include <AliFMDRawReader.h>
#include <AliFMDParameters.h>

struct SpectraCollector : public AliFMDInput
{
  //__________________________________________________________________
  SpectraCollector(UShort_t what, const char* src=0)
    : AliFMDInput("galice.root")
  {
    Setup(what, src);
  }
  //__________________________________________________________________
  SpectraCollector(const char* what, const char* src=0)
    : AliFMDInput("galice.root")
  {
    Setup(ParseLoad(what), src);
  }
  //__________________________________________________________________
  void Setup(UShort_t what, const char* src)
  { 
    switch (what) { 
    case kHits:		AddLoad(kHits); 	break;
    case kDigits:	AddLoad(kDigits); 	break;
    case kSDigits:	AddLoad(kSDigits); 	break;
    case kRaw:		AddLoad(kRaw); 		break;
    case kRawCalib:	AddLoad(kRawCalib); 	break;
    case kRecPoints:	AddLoad(kRecPoints); 	break;
    case kESD:		AddLoad(kESD); 		break;
    default: 
      Fatal("Spectra", "Unknown type %d reguested", what);
      return;
    }
    if ((what == kRaw || what == kRawCalib)) {
      if (!src || src[0] == '\0')
	Fatal("Spectra", "No raw input specified!");
      SetRawFile(src);
    }
    
    switch (what) { 
    case kHits:     fTop.SetAxis(500, 0, 1000); break;
    case kDigits:   // Fall-through
    case kSDigits:  // Fall-through
    case kRaw:      // Fall-through
    case kRawCalib: // This is the same for all kinds of digits
      fTop.SetAxis(1024, -.5, 1023.5); 
      break;
    case kRecPoints:  // Fall-through
    case kESD:        // This is the same for ESD and rec-points
      fTop.SetAxis(200, -.05, 19.95);
      break;
    }
  }
  //____________________________________________________________________
  Bool_t ProcessHit(AliFMDHit* h, TParticle* /* p */) 
  {
    if (h)
      fTop.Fill(h->Detector(), h->Ring(), h->Sector(), h->Strip(), h->Edep());
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessSDigit(AliFMDSDigit* d)
  {
    if (d)
      fTop.Fill(d->Detector(), d->Ring(), d->Sector(), d->Strip(), d->Counts());
    return kTRUE;
  }  
  //__________________________________________________________________
  Bool_t ProcessRawDigit(AliFMDDigit* digit)
  {
    return ProcessDigit(digit);
  }
  //__________________________________________________________________
  Bool_t ProcessDigit(AliFMDDigit* d)
  {
    if (d)
      fTop.Fill(d->Detector(), d->Ring(), d->Sector(),d->Strip(), d->Counts());
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessRawCalibDigit(AliFMDDigit* digit)
  {
    AliFMDParameters* parm = AliFMDParameters::Instance();
    UShort_t d          =  digit->Detector();
    Char_t   r          =  digit->Ring();
    UShort_t s          =  digit->Sector();
    UShort_t t          =  digit->Strip();
    Double_t g          =  parm->GetPulseGain(d, r, s, t);
    Double_t p          =  parm->GetPedestal(d, r, s, t);
    Double_t w          =  parm->GetPedestalWidth(d, r, s, t);
    UShort_t c          =  digit->Counts();
    Double_t x          =  0;
    if (fFMDReader && fFMDReader->IsZeroSuppressed(d-1))
      x = c + fFMDReader->NoiseFactor(d-1) * w;
    else 
      x = c - p;

    Double_t m = x / (g * parm->GetDACPerMIP());
    if (g < 0.1 || g > 10) m = 0;

    fTop.Fill(d, r, s, t, m);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessRecPoint(AliFMDRecPoint* r)
  {
    if (r)
      fTop.Fill(r->Detector(),r->Ring(),r->Sector(),r->Strip(),r->Particles());
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessESD(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
		    Float_t, Float_t m)  
  {
    if (m != AliESDFMD::kInvalidMult) 
      fTop.Fill(d,r,s,t,m);
    return kTRUE;
  }
    
  //__________________________________________________________________
  Bool_t Finish()
  {
    TBrowser* b = new TBrowser("b", "b");
    b->Add(&fTop);

#if 1
    TFile* file = TFile::Open("spectra.root", "RECREATE");
    // fTop.WriteOut(file);
    fTop.Write();
    file->Close();
#endif

    return kTRUE;
  }
  Spectra::Top fTop;
  ClassDef(SpectraCollector,0);
};

#else
//======================================================================
void
MakeSpectra(UShort_t what, Int_t n=1000, const char* src=0)
{
  gROOT->LoadMacro("$ALICE_ROOT.trunk/FMD/scripts/Compile.C");
  gSystem->AddIncludePath("-DSPECTRA_BUILD");
  const char* script = "$ALICE_ROOT.trunk/FMD/scripts/MakeSpectra.C";
  const char* here   = gSystem->BaseName(script);
  Int_t ret = 0;
  if ((ret = gSystem->CopyFile(gSystem->ExpandPathName(script), here, true))) { 
    Error("MakeSpectra", "Failed to copy %s to %s: %d", script, here, ret);
    return;
  }
  Compile(here, "+g");
  SpectraCollector* sd = new SpectraCollector(what, src);
  sd->Run(n);
}

#endif // SPECTRA_BUILD
//
// EOF
//
