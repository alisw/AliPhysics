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
#include <AliCDBManager.h>
#include <AliFMDParameters.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDSDigit.h>
#include <AliFMDRecPoint.h>
#include <AliESDFMD.h>
#include <AliFMDPattern.h>
#include <AliFMDSpectraDisplay.h>
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

//==================================================================
void AliFMDSpectraDisplayElement::MakeHistograms(TAxis* axis) 
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
void AliFMDSpectraDisplayElement::DoFill(Double_t v) 
{
  if (fFull) fFull->Fill(v);
}
//__________________________________________________________________
void AliFMDSpectraDisplayElement::Show(Option_t* option, 
				       Double_t l, Double_t h) 
{
  if (!fFull) return;
  gPad->SetLogy(fFull->GetMaximum() > 10);
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

//__________________________________________________________________
Int_t
AliFMDSpectraDisplayElement::Compare(const TObject*) const
{
  return -1;
}


//==================================================================
AliFMDSpectraDisplayTop::AliFMDSpectraDisplayTop(TGCompositeFrame& frame, 
						 TCanvas* canvas) 
  : AliFMDSpectraDisplayElement("All", "Everything"), 
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
    fAxis(0),
    fEntry(*(fList.AddItem(0, GetName(), this)))
{
  fContainer.AddFrame(&fList, &fHints);
  frame.AddFrame(&fContainer, &fHints);

  fList.Connect("Clicked(TGListTreeItem*,Int_t)", 
		"AliFMDSpectraDisplayTop", this, 
		"HandleEntry(TGListTreeItem*,Int_t)");
  fList.Connect("KeyPressed(TGListTreeItem*,ULong_t,ULong_t)", 
		"AliFMDSpectraDisplayTop", this, 
		"HandleKey(TGListTreeItem*,UInt_t,UInt_t)");
  fList.Connect("ReturnPressed(TGListTreeItem*)", 
		"AliFMDSpectraDisplayTop", this, 
		"HandleReturn(TGListTreeItem*)");
}
//____________________________________________________________________
void
AliFMDSpectraDisplayTop::SetAxis(TAxis* axis) 
{
  fAxis = axis;
  MakeHistograms(axis);
}
//____________________________________________________________________
void
AliFMDSpectraDisplayTop::ClearCanvas()
{
  if (!fCanvas) return;
  fCanvas->Clear();
}
  
//____________________________________________________________________
void
AliFMDSpectraDisplayTop::ClearList()
{
  fList.DeleteItem(fList.GetFirstItem());
  UpdateList();
}
  
  
//____________________________________________________________________
void
AliFMDSpectraDisplayTop::HandleReturn(TGListTreeItem * f)
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
AliFMDSpectraDisplayTop::HandleKey(TGListTreeItem * f, UInt_t keysym, UInt_t)
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
  // if (next) gClient->NeedRedraw(&fList);
  if (next && next != f) {
    fList.ClearHighlighted();
    fList.SetSelected(next);
    HandleEntry(next,0);
  }
  if (next) gClient->NeedRedraw(&fList);
}

//____________________________________________________________________
void
AliFMDSpectraDisplayTop::HandleEntry(TGListTreeItem* entry, Int_t /*id*/) 
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
AliFMDSpectraDisplayTop::UpdateList() 
{
  gClient->NeedRedraw(&fList);
}

//____________________________________________________________________
void
AliFMDSpectraDisplayTop::UpdateCanvas() 
{
  if (!fCanvas) return;
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();
}
//____________________________________________________________________
TObject* AliFMDSpectraDisplayTop::Current() const
{
  if (!fCurrentEntry) return 0;
  if (!fCurrentEntry->GetUserData()) return 0;
  return static_cast<TObject*>(fCurrentEntry->GetUserData());
}
  
//__________________________________________________________________
AliFMDSpectraDisplayDetector& AliFMDSpectraDisplayTop::GetOrAdd(UShort_t id) 
{ 
  Int_t     idx = id - 1;
  AliFMDSpectraDisplayDetector* d   = 0;
  if (fChildren.GetEntriesFast() <= idx ||
      !(d = static_cast<AliFMDSpectraDisplayDetector*>(fChildren.At(idx))))  { 
    d = new AliFMDSpectraDisplayDetector(id, *this);
    fChildren.AddAtAndExpand(d, idx);
    // GetTop().GetList().SortChildren(&fEntry);
    // GetTop().GetList().Sort(&(d->GetEntry()));
  }
  return *d;
}
//__________________________________________________________________
void AliFMDSpectraDisplayTop::Fill(UShort_t det, Char_t ring, 
	       UShort_t sec, UShort_t str, Double_t v) 
{ 
  AliFMDSpectraDisplayDetector& d = GetOrAdd(det);
  d.Fill(ring, sec, str, v);
  DoFill(v);
}
//__________________________________________________________________
Int_t
AliFMDSpectraDisplayTop::Compare(const TObject*) const
{
  return -1;
}
//==================================================================
AliFMDSpectraDisplayDetector::AliFMDSpectraDisplayDetector(UShort_t det, 
				           AliFMDSpectraDisplayTop& tree) 
  : AliFMDSpectraDisplayElement(Form("FMD%d", det), "FMD Sub-detector"), 
    fId(det), 
    fParent(tree),
    fChildren(0),
    fEntry(*(tree.GetList().AddItem(&(tree.GetEntry()), GetName())))
{
  fEntry.SetUserData(this);
  fEntry.SetText(GetName());
  if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
}
//__________________________________________________________________
AliFMDSpectraDisplayRing& AliFMDSpectraDisplayDetector::GetOrAdd(Char_t id) 
{ 
  Int_t idx = (id == 'I' || id == 'i') ? 0 : 1;
  AliFMDSpectraDisplayRing* r   = 0;;
  if (fChildren.GetEntriesFast() <= idx ||
      !(r = static_cast<AliFMDSpectraDisplayRing*>(fChildren.At(idx))))  { 
    r = new AliFMDSpectraDisplayRing(id, *this);
    fChildren.AddAtAndExpand(r, idx);
    // GetTop().GetList().SortChildren(&fEntry);
    // GetTop().GetList().Sort(&(r->GetEntry()));
  }
  return *r;
}
//__________________________________________________________________
void AliFMDSpectraDisplayDetector::Fill(Char_t ring, UShort_t sec, 
					UShort_t str, Double_t v) 
{ 
  AliFMDSpectraDisplayRing& r = GetOrAdd(ring);
  r.Fill(sec, str, v);
  DoFill(v);
}  
//__________________________________________________________________
Int_t
AliFMDSpectraDisplayDetector::Compare(const TObject* o) const
{
  std::cout << "Comparing detector to a " << o->ClassName() << std::endl;
  if (o->IsA() == AliFMDSpectraDisplayDetector::Class()) { 
    const AliFMDSpectraDisplayDetector* ro = 
      static_cast<const AliFMDSpectraDisplayDetector*>(o);
    return (Id() <  ro->Id() ? -1 : 
	    Id() == ro->Id() ? 0 : 1);
  }
  return -1;
}
//==================================================================
AliFMDSpectraDisplayRing::AliFMDSpectraDisplayRing(Char_t id,
                            AliFMDSpectraDisplayDetector& d) 
  : AliFMDSpectraDisplayElement(Form("FMD%d%c", d.Id(), id), "FMD Ring"), 
    fParent(d),
    fId(id),
    fChildren(0),
    fEntry(*(GetTop().GetList().AddItem(&(d.GetEntry()), GetName(), this)))
{
  fEntry.SetText(GetName());
  if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
}
//__________________________________________________________________
AliFMDSpectraDisplaySector& AliFMDSpectraDisplayRing::GetOrAdd(UShort_t id) 
{ 
  AliFMDSpectraDisplaySector* s = 0;
  if (fChildren.GetEntriesFast() <= id ||
      !(s = static_cast<AliFMDSpectraDisplaySector*>(fChildren.At(id)))) {
    s = new AliFMDSpectraDisplaySector(id, *this);
    fChildren.AddAtAndExpand(s, id);
    // GetTop().GetList().SortChildren(&fEntry);
    // GetTop().GetList().Sort(&(s->GetEntry()));
  }
  return *s;
}
//__________________________________________________________________
void AliFMDSpectraDisplayRing::Fill(UShort_t sec, UShort_t str, Double_t v) 
{ 
  AliFMDSpectraDisplaySector& s = GetOrAdd(sec);
  s.Fill(str, v);
  DoFill(v);
}
//__________________________________________________________________
Int_t
AliFMDSpectraDisplayRing::Compare(const TObject* o) const
{
  std::cout << "Comparing ring to a " << o->ClassName() << std::endl;
  if (o->IsA() == AliFMDSpectraDisplayRing::Class()) { 
    const AliFMDSpectraDisplayRing* ro = 
      static_cast<const AliFMDSpectraDisplayRing*>(o);
    return (Id() <  ro->Id() ? -1 : 
	    Id() == ro->Id() ? 0 : 1);
  }
  return -1;
}
//==================================================================
AliFMDSpectraDisplaySector::AliFMDSpectraDisplaySector(UShort_t id, 
                                      AliFMDSpectraDisplayRing& r) 
  : AliFMDSpectraDisplayElement(Form("FMD%d%c_%02d",r.DetectorId(),r.Id(),id), 
				"FMD Sector"), 
    fParent(r),
    fId(id),
    fChildren(0),
    fEntry(*(GetTop().GetList().AddItem(&(r.GetEntry()), GetName(), this)))
{
  fEntry.SetText(GetName());
  if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
}
//__________________________________________________________________
AliFMDSpectraDisplayStrip& AliFMDSpectraDisplaySector::GetOrAdd(UShort_t id) 
{ 
  AliFMDSpectraDisplayStrip* s = 0;
  if (fChildren.GetEntriesFast() <= id || 
      !(s = static_cast<AliFMDSpectraDisplayStrip*>(fChildren.At(id)))) {
    s = new AliFMDSpectraDisplayStrip(id, *this);
    fChildren.AddAtAndExpand(s, id);
    // GetTop().GetList().SortChildren(&fEntry);
    // GetTop().GetList().Sort(&(s->GetEntry()));
  }
  return *s;
}
//__________________________________________________________________
void AliFMDSpectraDisplaySector::Fill(UShort_t str, Double_t v) 
{ 
  AliFMDSpectraDisplayStrip& s = GetOrAdd(str);
  s.Fill(v);
  DoFill(v);
}
//__________________________________________________________________
Int_t
AliFMDSpectraDisplaySector::Compare(const TObject* o) const
{
  std::cout << "Comparing sector to a " << o->ClassName() << std::endl;
  if (o->IsA() == AliFMDSpectraDisplaySector::Class()) { 
    const AliFMDSpectraDisplaySector* ro = 
      static_cast<const AliFMDSpectraDisplaySector*>(o);
    return (Id() <  ro->Id() ? -1 : 
	    Id() == ro->Id() ? 0 : 1);
  }
  return -1;
}
//==================================================================
AliFMDSpectraDisplayStrip::AliFMDSpectraDisplayStrip(UShort_t id, 
                                  AliFMDSpectraDisplaySector& s) 
  : AliFMDSpectraDisplayElement(Form("FMD%d%c_%02d_%03d", 
				     s.DetectorId(), s.RingId(), 
				     s.Id(), id), "FMD Strip"), 
    fParent(s),
    fId(id),
    fEntry(*(GetTop().GetList().AddItem(&(s.GetEntry()), GetName(), this)))
{
  fEntry.SetText(GetName());
  fEntry.SetPictures(GetTop().GetH1Pic(), GetTop().GetH1Pic());
  if (GetTop().GetAxis()) MakeHistograms(GetTop().GetAxis());
}
//__________________________________________________________________
void AliFMDSpectraDisplayStrip::Fill(Double_t v) 
{ 
  DoFill(v);
}
//__________________________________________________________________
Int_t
AliFMDSpectraDisplayStrip::Compare(const TObject* o) const
{
  std::cout << "Comparing strip to a " << o->ClassName() << std::endl;
  if (o->IsA() == AliFMDSpectraDisplayStrip::Class()) { 
    const AliFMDSpectraDisplayStrip* ro = 
      static_cast<const AliFMDSpectraDisplayStrip*>(o);
    return (Id() <  ro->Id() ? -1 : 
	    Id() == ro->Id() ? 0 : 1);
  }
  return -1;
}

//==================================================================
AliFMDSpectraDisplay::AliFMDSpectraDisplay()
  : AliFMDPattern(), 
    fSelector(gClient->GetRoot(), 100, 100), 
    fTop(fSelector, fAux)
{
  AddLoad(AliFMDInput::kRaw);
  SetName("RAW");
  SetTitle("RAW");
  
  fTop.Connect("SelectionChanged()", 
	       "AliFMDSpectraDisplay", this, "HandleDraw()");
  
  fSelector.MapSubwindows();
  fSelector.Resize(fSelector.GetDefaultSize());
  fSelector.MapWindow();
}

//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::HandleDraw()
{
  TObject* user = fTop.Current();
  if (!user) return kFALSE;
  if (!user->InheritsFrom(AliFMDSpectraDisplayElement::Class())) { 
    Warning("HandleDraw", "%s does not inherit from Spectra::Element", 
	    user->GetName());
    return kFALSE;
  }
  fAux->cd();
  AliFMDSpectraDisplayElement* e 
    = static_cast<AliFMDSpectraDisplayElement*>(user);
  e->Show("hist", fSlider->GetMinimum(), fSlider->GetMaximum());
  fAux->Modified();
  fAux->Update();
  fAux->cd();
  return kTRUE;
}
//__________________________________________________________________
void 
AliFMDSpectraDisplay::MakeAux()
{
  AliFMDPattern::MakeAux();
  if (!fAux) return;
  fTop.SetAxis(fSpec->GetXaxis());
}
//__________________________________________________________________
void 
AliFMDSpectraDisplay::DrawAux()
{
  // Draw in the Aux the canvas 
  // For example draw the spectra 
  // or such stuff 
  if (fTop.Current() != &fTop && HandleDraw()) return;
  AliFMDPattern::DrawAux();
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessHit(AliFMDHit* hit, TParticle* p) 
{
  fTop.Fill(hit->Detector(), 
	    hit->Ring(), 
	    hit->Sector(),
	    hit->Strip(), 
	    hit->Edep());
  return AliFMDPattern::ProcessHit(hit, p);
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessDigit(AliFMDDigit* digit)
{
  fTop.Fill(digit->Detector(), 
	    digit->Ring(), 
	    digit->Sector(),
	    digit->Strip(), 
	    digit->Counts());
  return AliFMDDisplay::ProcessDigit(digit);
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessSDigit(AliFMDSDigit* sdigit)
{
  fTop.Fill(sdigit->Detector(), 
	    sdigit->Ring(), 
	    sdigit->Sector(),
	    sdigit->Strip(), 
	    sdigit->Counts());
  return AliFMDDisplay::ProcessSDigit(sdigit);
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessRawDigit(AliFMDDigit* digit)
{
  return ProcessDigit(digit);
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessRecPoint(AliFMDRecPoint* recpoint)
{
  fTop.Fill(recpoint->Detector(), 
	    recpoint->Ring(), 
	    recpoint->Sector(),
	    recpoint->Strip(), 
	    recpoint->Particles());
  return AliFMDDisplay::ProcessRecPoint(recpoint);
}
//__________________________________________________________________
Bool_t 
AliFMDSpectraDisplay::ProcessESD(UShort_t det, Char_t rng, UShort_t sec, 
				 UShort_t str, Float_t x, Float_t mult)
{
  fTop.Fill(det, rng, sec, str, mult);
  return AliFMDDisplay::ProcessESD(det, rng, sec, str, x, mult);
}
//__________________________________________________________________
//
// EOF
//
