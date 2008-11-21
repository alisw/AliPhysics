// 
//____________________________________________________________________
//
// $Id: DrawHits.C 22496 2007-11-26 13:50:44Z cholm $
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//

void
SpectraMonitor(const char* file="", 
	       Int_t       runno=0, 
	       const char* cdbSrc="local://$ALICE_ROOT", 
	       UShort_t    over=0)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libFMDanalysis.so");
  gSystem->Load("libFMDutil.so");
  // AliLog::SetModuleDebugLevel("FMD", 8);
  TString fname(file);
  if (fname.CompareTo("help", TString::kIgnoreCase) == 0) { 
    std::cout << "Usage: RunSpectraMonitor(<src>[,<runno>[,<cdb>[,<over>]]])\n"
	      << "\n"
	      << "Where:\n\n"
	      << "   <src>         Is a data source (online, file)\n"
	      << "   <runno>       Is the (optional) run number\n"
	      << "   <cdb>         Is the (optional) CDB storage\n"
	      << "   <over>        Is the (optional) over sampling rate\n\n"
	      << "Defaults are <runno>=0 and cdb=\"local://$ALICE_ROOT\"\n" 
	      << "<over> allows one to override the CDB setting.  Default\n"
	      << "is to use the CDB setting.\n\n"
	      << "Note: This script _must_ be compiled with ACLic"
	      << std::endl;
    return;
  }
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbSrc);
  cdb->SetRun(runno);
  UInt_t what = (AliFMDParameters::kPulseGain      |
		 AliFMDParameters::kPedestal       |
		 AliFMDParameters::kDeadMap        |
		 AliFMDParameters::kZeroSuppression|
		 AliFMDParameters::kAltroMap       |
		 AliFMDParameters::kStripRange);
  if (over != 0) what |= AliFMDParameters::kSampleRate;
  AliFMDParameters::Instance()->Init(kFALSE, what);
  if (over != 0) AliFMDParameters::Instance()->SetSampleRate(over);
  
  AliFMDSpectraDisplay* d = new AliFMDSpectraDisplay;
  // d->AddLoad(AliFMDInput::kRaw);
  d->AddLoad(AliFMDInput::kDigits);
  // d->SetRawFile(file);
  d->Run();
}

#if 0
#include <AliCDBManager.h>
#include <AliFMDParameters.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDSDigit.h>
#include <AliFMDRecPoint.h>
#include <AliESDFMD.h>
#include <AliFMDPattern.h>
#include <TBrowser.h>
#include <TDirectory.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TF1.h>
#include <iostream>
#include <TStyle.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TGFrame.h>
#include <TGCanvas.h>
#include <TGListTree.h>
#include <TGClient.h>
#include <TSystem.h>
#include <KeySymbols.h>
#include <TClass.h>
#include <RQ_OBJECT.h>
#include <TSlider.h>

namespace Spectra { 
  struct Top;
  struct Detector;
  struct Ring;
  struct Sector;
  struct Strip;
  
  //__________________________________________________________________
  struct Element : public TNamed 
  {
    virtual ~Element() {}
    virtual void Draw(Option_t* option, Double_t l, Double_t h);
    virtual Top& GetTop() = 0;
    virtual void MakeHistograms(TAxis* axis);
  protected:
    Element(const char* name, const char* title) 
      : TNamed(name, title), fFull(0), fCut(0)
    {}
    void DoFill(Double_t v);
    TH1* fFull;
    TH1* fCut;
  };

  //__________________________________________________________________
  struct Top : public Element
  {
    RQ_OBJECT("Spectra::Top")
    public:
    Top(TGCompositeFrame& frame, TCanvas* canvas);
    TGListTree& GetList() { return fList; }
    Detector& GetOrAdd(UShort_t id);
    void      Fill(UShort_t d, Char_t ring, 
		   UShort_t sec, UShort_t str, Double_t v);
    TAxis* GetAxis() { return fAxis; }
    void SetAxis(TAxis* a);
    Top&   GetTop() { return *this; }
       
    /** Handle entries 
	@param e  selected entry, if any 
	@param id Id of entry */
    virtual void HandleEntry(TGListTreeItem* e, Int_t id);
    /** Handle key strokes 
	@param f      Item selected, if any 
	@param keysym Key symbol 
	@param mask   Modifier mask */
    virtual void HandleKey(TGListTreeItem* f, UInt_t keysym, UInt_t mask);
    /** Handle Return 
	@param f Selected item, if any */
    virtual void HandleReturn(TGListTreeItem* f);
    /** Clear the list */
    virtual void ClearList();
    /** Clear the canvas */ 
    virtual void ClearCanvas();
    /** Update canvas */ 
    virtual void UpdateCanvas();
    /** Update canvas */ 
    virtual void UpdateList();
    /** Return the currently selected entry */ 
    TGListTreeItem* CurrentEntry() const { return fCurrentEntry; }
    /** @return the currently selected user data (possibly 0) */
    TObject* Current() const;
    /** Selection changed signal */
    void SelectionChanged() { Emit("SelectionChanged()"); }//*SIGNAL*
    /** Get Picture for 1D histogram */
    const TGPicture* GetH1Pic() { return fHist1DIcon; }
    /** 2D Histogram Icon */
    const TGPicture* GetH2Pic() { return fHist2DIcon; }
    /** 3D Histogram Icon */
    const TGPicture* GetH3Pic() { return fHist3DIcon; }
    /** Graph Icon */
    const TGPicture* GetGPic() { return fGraphIcon; }
    TGListTreeItem&  GetEntry() { return fEntry; }
  protected:
    TGLayoutHints fHints;
    TGCanvas      fContainer;
    TGListTree    fList;
    TObjArray     fChildren;
    TGListTreeItem* fCurrentEntry;
    TCanvas*        fCanvas;
    /** 1D Histogram Icon */
    const TGPicture*    fHist1DIcon;
    /** 2D Histogram Icon */
    const TGPicture*    fHist2DIcon;
    /** 3D Histogram Icon */
    const TGPicture*    fHist3DIcon;
    /** Graph Icon */
    const TGPicture*    fGraphIcon;
    /** The axis to use */ 
    TAxis* fAxis;
    /** Top entry */ 
    TGListTreeItem& fEntry;
  };

  //__________________________________________________________________
  struct Detector : public Element 
  {
    Detector(UShort_t det, Top& top);
    UShort_t Id() const { return fId; }
    Top&     GetTop() { return fParent; }
    Top&     GetParent() { return fParent; }
    Ring&    GetOrAdd(Char_t id);
    void     Fill(Char_t ring, UShort_t sec, UShort_t str, Double_t v);
    TGListTreeItem&  GetEntry() { return fEntry; }
  protected:
    UShort_t        fId;
    Top&            fParent;
    TObjArray       fChildren;
    TGListTreeItem& fEntry;
  };

  //__________________________________________________________________
  struct Ring : public Element 
  {
    Ring(Char_t id, Detector& d);
    Char_t    Id() const { return fId; }
    UShort_t  DetectorId() { return fParent.Id(); }
    Top&      GetTop() { return fParent.GetTop(); }
    Detector& GetDetector() { return GetParent(); }
    Detector& GetParent() { return fParent; }
    Sector&   GetOrAdd(UShort_t id);
    void      Fill(UShort_t sec, UShort_t str, Double_t v);
    TGListTreeItem&  GetEntry() { return fEntry; }
  protected:
    Detector&       fParent;
    Char_t          fId;
    TObjArray       fChildren;
    TGListTreeItem& fEntry;
  };
  
  //__________________________________________________________________
  struct Sector : public Element 
  {
    Sector(UShort_t id, Ring& r);
    UShort_t  Id() const    { return fId; }
    UShort_t  DetectorId()  { return fParent.DetectorId(); }
    Char_t    RingId()      { return fParent.Id(); }
    Top&      GetTop()      { return fParent.GetTop(); }
    Detector& GetDetector() { return fParent.GetDetector(); }
    Ring&     GetRing()     { return fParent; }
    Ring&     GetParent()   { return fParent; }
    Strip&    GetOrAdd(UShort_t id);
    void      Fill(UShort_t str, Double_t v);
    TGListTreeItem&  GetEntry() { return fEntry; }
  protected:
    Ring&           fParent;
    UShort_t        fId;
    TObjArray       fChildren;
    TGListTreeItem& fEntry;
  };

  //__________________________________________________________________
  struct Strip : public Element 
  {
    Strip(UShort_t id, Sector& s);
    UShort_t  Id() const    { return fId; }
    UShort_t  DetectorId()  { return fParent.DetectorId(); }
    Char_t    RingId()      { return fParent.RingId(); }
    UShort_t  SectorId()    { return fParent.Id(); }
    Top&      GetTop()      { return fParent.GetTop(); }
    Detector& GetDetector() { return fParent.GetDetector(); }
    Ring&     GetRing()     { return fParent.GetRing(); }
    Sector&   GetSector()   { return fParent; }
    Sector&   GetParent()   { return fParent; }
    void      Fill(Double_t v);
    TGListTreeItem&  GetEntry() { return fEntry; }
  protected:
    Sector&         fParent;
    UShort_t        fId;
    TGListTreeItem& fEntry;
  };

  //==================================================================
  void Element::MakeHistograms(TAxis* axis) 
  {
    if (fFull) return;
    if (axis->IsVariableBinSize()) {
      fFull = new TH1F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
		       axis->GetXbins()->fN, axis->GetXbins()->fArray);
      fCut = new TH1F(Form("c_%s", GetName()), 
		      Form("%s restricted spectra", GetName()),
		      axis->GetXbins()->fN, axis->GetXbins()->fArray);
    }
    else { 
      fFull = new TH1F(Form("f_%s", GetName()), Form("%s spectra", GetName()),
		       axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
      fCut = new TH1F(Form("c_%s", GetName()), 
		      Form("%s restricted spectra", GetName()),
		      axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
    }
    fFull->SetFillColor(kRed);
    fFull->SetFillStyle(3001);
    fCut->SetFillColor(kBlue);
    fCut->SetFillStyle(3001);
  }
  //__________________________________________________________________
  void Element::DoFill(Double_t v) 
  {
    if (fFull) fFull->Fill(v);
  }
  //__________________________________________________________________
  void Element::Draw(Option_t* option, Double_t l, Double_t h) 
  {
    if (!fFull) return;
    if (fFull->GetMinimum() < 0) gPad->SetLogy(kFALSE);
    fFull->Draw(option);
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
  }
  //==================================================================
  Top::Top(TGCompositeFrame& frame, TCanvas* canvas) 
    : Element("All", "Everything"), 
      fHints(kLHintsExpandX|kLHintsExpandY,3, 3, 3, 3),
      fContainer(&frame, 200, 350), 
      fList(&fContainer, kHorizontalFrame),
      fChildren(0),
      fCurrentEntry(0),
      fCanvas(canvas),
      fHist1DIcon(gClient->GetPicture("h1_t.xpm")),
      fHist2DIcon(gClient->GetPicture("h2_t.xpm")),
      fHist3DIcon(gClient->GetPicture("h3_t.xpm")),
      fGraphIcon(gClient->GetPicture("graph.xpm")), 
      fEntry(*(fList.AddItem(0, GetName(), this)))
  {
    fContainer.AddFrame(&fList, &fHints);
    frame.AddFrame(&fContainer, &fHints);

    fList.Connect("Clicked(TGListTreeItem*,Int_t)", "Spectra::Top", this, 
		  "HandleEntry(TGListTreeItem*,Int_t)");
    fList.Connect("KeyPressed(TGListTreeItem*,ULong_t,ULong_t)", 
		  "Spectra::Top", this, 
		  "HandleKey(TGListTreeItem*,UInt_t,UInt_t)");
    fList.Connect("ReturnPressed(TGListTreeItem*)", "Spectra::Top", this, 
		  "HandleReturn(TGListTreeItem*)");
  }
  //____________________________________________________________________
  void
  Top::SetAxis(TAxis* axis) 
  {
    fAxis = axis;
    MakeHistograms(axis);
  }
  //____________________________________________________________________
  void
  Top::ClearCanvas()
  {
    if (!fCanvas) return;
    fCanvas->Clear();
  }
  
  //____________________________________________________________________
  void
  Top::ClearList()
  {
    fList.DeleteItem(fList.GetFirstItem());
    UpdateList();
  }
  
  
  //____________________________________________________________________
  void
  Top::HandleReturn(TGListTreeItem * f)
  {
    
    if (!f) { 
      fList.UnselectAll(kFALSE);
      fList.SetSelected(0);
      return;
    }
    fList.ToggleItem(f);
    UpdateList();
  }
  
  
  //____________________________________________________________________
  void
  Top::HandleKey(TGListTreeItem * f, UInt_t keysym, UInt_t /*mask*/)
  {
    if (!f) { 
      fList.UnselectAll(kFALSE);
      fList.SetSelected(0);
      return;
    }
    TGListTreeItem* next = 0;
    switch (keysym) {
    case kKey_Up:
      next = f->GetPrevSibling();
    if (!next) { 
      next = f->GetParent();
      if (next) fList.CloseItem(next);
    }
    break;
    case kKey_Down:
      next = f->GetNextSibling();
      if (!next && f->GetParent()) {
	next = f->GetParent()->GetNextSibling();
	fList.CloseItem(f->GetParent());
      }
      break;
    case kKey_Left:
      next = f->GetParent();
      if (next) fList.CloseItem(next);
      break;
    case kKey_Right:
      next = f->GetFirstChild();
      if (next) fList.OpenItem(f);
      break;
    case kKey_PageUp:
      fList.PageUp(kTRUE);
      next = fList.GetSelected();
      break;
    case kKey_PageDown:
      fList.PageDown(kTRUE);
      next = fList.GetSelected();
      break;
    }
    if (next) gClient->NeedRedraw(&fList);
    if (next && next != f) {
      fList.ClearHighlighted();
      fList.SetSelected(next);
      HandleEntry(next,0);
    }
  }

  //____________________________________________________________________
  void
  Top::HandleEntry(TGListTreeItem* entry, Int_t /*id*/) 
  {
    TGListTreeItem* old = fCurrentEntry;
    if (entry) {
      if (!entry->GetUserData()) return;
      fCurrentEntry = entry;
    }
    else {
      fCurrentEntry = 0;
      ClearCanvas();
    }
    if (old != fCurrentEntry && fCanvas) fCanvas->cd();
    SelectionChanged();
  }

  //____________________________________________________________________
  void
  Top::UpdateList() 
  {
    gClient->NeedRedraw(&fList);
  }

  //____________________________________________________________________
  void
  Top::UpdateCanvas() 
  {
    if (!fCanvas) return;
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
  }
  //____________________________________________________________________
  TObject* Top::Current() const
  {
    if (!fCurrentEntry) return 0;
    if (!fCurrentEntry->GetUserData()) return 0;
    return static_cast<TObject*>(fCurrentEntry->GetUserData());
  }
  
  //__________________________________________________________________
  Detector& Top::GetOrAdd(UShort_t id) 
  { 
    Int_t     idx = id - 1;
    Detector* d   = 0;
    if (fChildren.GetEntriesFast() <= idx ||
	!(d = static_cast<Detector*>(fChildren.At(idx))))  { 
      d = new Detector(id, *this);
      fChildren.AddAtAndExpand(d, idx);
      // fList.SortChildren(&fEntry);
    }
    return *d;
  }
  //__________________________________________________________________
  void Top::Fill(UShort_t det, Char_t ring, 
		 UShort_t sec, UShort_t str, Double_t v) 
  { 
    Detector& d = GetOrAdd(det);
    d.Fill(ring, sec, str, v);
    DoFill(v);
  }
  //==================================================================
  Detector::Detector(UShort_t det, Top& tree) 
    : Element(Form("FMD%d", det), "FMD Sub-detector"), 
      fId(det), 
      fParent(tree),
      fChildren(0),
      fEntry(*(tree.GetList().AddItem(&(tree.GetEntry()), GetName())))
  {
    fEntry.SetUserData(this);
    if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
  }
  //__________________________________________________________________
  Ring& Detector::GetOrAdd(Char_t id) 
  { 
    Int_t idx = (id == 'I' || id == 'i') ? 0 : 1;
    Ring* r   = 0;;
    if (fChildren.GetEntriesFast() <= idx ||
	!(r = static_cast<Ring*>(fChildren.At(idx))))  { 
      r = new Ring(id, *this);
      fChildren.AddAtAndExpand(r, idx);
      // GetTop().GetList().SortChildren(&fEntry);
    }
    return *r;
  }
  //__________________________________________________________________
  void Detector::Fill(Char_t ring, UShort_t sec, UShort_t str, Double_t v) 
  { 
    Ring& r = GetOrAdd(ring);
    r.Fill(sec, str, v);
    DoFill(v);
  }  
  //==================================================================
  Ring::Ring(Char_t id, Detector& d) 
    : Element(Form("FMD%d%c", d.Id(), id), "FMD Ring"), 
      fParent(d),
      fId(id),
      fChildren(0),
      fEntry(*(GetTop().GetList().AddItem(&(d.GetEntry()), GetName(), this)))
  {
    if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
  }
  //__________________________________________________________________
  Sector& Ring::GetOrAdd(UShort_t id) 
  { 
    Sector* s = 0;
    if (fChildren.GetEntriesFast() <= id ||
	!(s = static_cast<Sector*>(fChildren.At(id)))) {
      s = new Sector(id, *this);
      fChildren.AddAtAndExpand(s, id);
      // GetTop().GetList().SortChildren(&fEntry);
    }
    return *s;
  }
  //__________________________________________________________________
  void Ring::Fill(UShort_t sec, UShort_t str, Double_t v) 
  { 
    Sector& s = GetOrAdd(sec);
    s.Fill(str, v);
    DoFill(v);
  }
  //==================================================================
  Sector::Sector(UShort_t id, Ring& r) 
    : Element(Form("FMD%d%c_%02d", r.DetectorId(), r.Id(), id), "FMD Sector"), 
      fParent(r),
      fId(id),
      fChildren(0),
      fEntry(*(GetTop().GetList().AddItem(&(r.GetEntry()), GetName(), this)))
  {
    if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
  }
  //__________________________________________________________________
  Strip& Sector::GetOrAdd(UShort_t id) 
  { 
    Strip* s = 0;
    if (fChildren.GetEntriesFast() <= id || 
	!(s = static_cast<Strip*>(fChildren.At(id)))) {
      s = new Strip(id, *this);
      fChildren.AddAtAndExpand(s, id);
      // GetTop().GetList().SortChildren(&fEntry);
    }
    return *s;
  }
  //__________________________________________________________________
  void Sector::Fill(UShort_t str, Double_t v) 
  { 
    Strip& s = GetOrAdd(str);
    s.Fill(v);
    DoFill(v);
  }
  //==================================================================
  Strip::Strip(UShort_t id, Sector& s) 
    : Element(Form("FMD%d%c_%02d_%03d", s.DetectorId(), s.RingId(), 
		   s.Id(), id), "FMD Strip"), 
      fParent(s),
      fId(id),
      fEntry(*(GetTop().GetList().AddItem(&(s.GetEntry()), GetName(), this)))
  {
    fEntry.SetPictures(GetTop().GetH1Pic(), GetTop().GetH1Pic());
    if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
  }
  //__________________________________________________________________
  void Strip::Fill(Double_t v) 
  { 
    DoFill(v);
  }
}


//====================================================================
class SpectraMonitor : public AliFMDPattern
{
private:
  TGMainFrame  fSelector;
  Spectra::Top fTop;
public:
  SpectraMonitor()
    : AliFMDPattern(), 
      fSelector(gClient->GetRoot(), 100, 100), 
      fTop(fSelector, fAux)
  {
    AddLoad(AliFMDInput::kRaw);
    SetName("RAW");
    SetTitle("RAW");

    fTop.Connect("SelectionChanged()", "SpectraMonitor", this, "HandleDraw()");

    fSelector.MapSubwindows();
    fSelector.Resize(fSelector.GetDefaultSize());
    fSelector.MapWindow();
  }
  // void SetAxis(TAxis* axis) { fTop.SetAxis(axis); }
  Bool_t HandleDraw()
  {
    TObject* user = fTop.Current();
    if (!user) return kFALSE;
    if (!user->InheritsFrom(Spectra::Element::Class())) { 
      Warning("HandleDraw", "%s does not inherit from Spectra::Element", 
	      user->GetName());
      return kFALSE;
    }
    fAux->cd();
    Spectra::Element* e = static_cast<Spectra::Element*>(user);
    e->Draw("hist", fSlider->GetMinimum(), fSlider->GetMaximum());
    fAux->Modified();
    fAux->Update();
    fAux->cd();
    return kTRUE;
  }
  void MakeAux()
  {
    AliFMDPattern::MakeAux();
    if (!fAux) return;
    fTop.SetAxis(fSpec->GetXaxis());
  }
  void DrawAux()
  {
    // Draw in the Aux the canvas 
    // For example draw the spectra 
    // or such stuff 
    if (fTop.Current() != &fTop && HandleDraw()) return;
    AliFMDPattern::DrawAux();
  }
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p) 
  {
    fTop.Fill(hit->Detector(), 
	      hit->Ring(), 
	      hit->Sector(),
	      hit->Strip(), 
	      hit->Edep());
    return AliFMDPattern::ProcessHit(hit, p);
  }
  Bool_t ProcessDigit(AliFMDDigit* digit)
  {
    fTop.Fill(digit->Detector(), 
	      digit->Ring(), 
	      digit->Sector(),
	      digit->Strip(), 
	      digit->Counts());
    return AliFMDDisplay::ProcessDigit(digit);
  }
  Bool_t ProcessSDigit(AliFMDSDigit* sdigit)
  {
    fTop.Fill(sdigit->Detector(), 
	      sdigit->Ring(), 
	      sdigit->Sector(),
	      sdigit->Strip(), 
	      sdigit->Counts());
    return AliFMDDisplay::ProcessSDigit(sdigit);
  }
  Bool_t ProcessRawDigit(AliFMDDigit* digit)
  {
    return ProcessDigit(digit);
  }
  Bool_t ProcessRecPoint(AliFMDRecPoint* recpoint)
  {
    fTop.Fill(recpoint->Detector(), 
	      recpoint->Ring(), 
	      recpoint->Sector(),
	      recpoint->Strip(), 
	      recpoint->Particles());
    return AliFMDDisplay::ProcessRecPoint(recpoint);
  }
  Bool_t ProcessESD(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
		    Float_t x, Float_t mult)
  {
    fTop.Fill(det, rng, sec, str, mult);
    return AliFMDDisplay::ProcessESD(det, rng, sec, str, x, mult);
  }
};

void
RunSpectraMonitor(const char* file="", 
		  Int_t       runno=0, 
		  const char* cdbSrc="local://$ALICE_ROOT", 
		  UShort_t    over=0)
{
  // AliLog::SetModuleDebugLevel("FMD", 8);
  TString fname(file);
  if (fname.CompareTo("help", TString::kIgnoreCase) == 0) { 
    std::cout << "Usage: RunSpectraMonitor(<src>[,<runno>[,<cdb>[,<over>]]])\n"
	      << "\n"
	      << "Where:\n\n"
	      << "   <src>         Is a data source (online, file)\n"
	      << "   <runno>       Is the (optional) run number\n"
	      << "   <cdb>         Is the (optional) CDB storage\n"
	      << "   <over>        Is the (optional) over sampling rate\n\n"
	      << "Defaults are <runno>=0 and cdb=\"local://$ALICE_ROOT\"\n" 
	      << "<over> allows one to override the CDB setting.  Default\n"
	      << "is to use the CDB setting.\n\n"
	      << "Note: This script _must_ be compiled with ACLic"
	      << std::endl;
    return;
  }
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbSrc);
  cdb->SetRun(runno);
  UInt_t what = (AliFMDParameters::kPulseGain      |
		 AliFMDParameters::kPedestal       |
		 AliFMDParameters::kDeadMap        |
		 AliFMDParameters::kZeroSuppression|
		 AliFMDParameters::kAltroMap       |
		 AliFMDParameters::kStripRange);
  if (over != 0) what |= AliFMDParameters::kSampleRate;
  AliFMDParameters::Instance()->Init(kFALSE, what);
  if (over != 0) AliFMDParameters::Instance()->SetSampleRate(over);
  
  SpectraMonitor* d = new SpectraMonitor;
  d->AddLoad(AliFMDInput::kRaw);
  d->SetRawFile(file);
  d->Run();
}
#endif

  

//____________________________________________________________________
//
// EOF
//
