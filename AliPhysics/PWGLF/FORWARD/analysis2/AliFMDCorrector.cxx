//
// This class calculates the exclusive charged particle density
// in each for the 5 FMD rings. 
//
#include "AliFMDCorrector.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliFMDCorrVertexBias.h"
#include "AliFMDCorrMergingEfficiency.h"
#include "AliFMDCorrAcceptance.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TROOT.h>
#include <THStack.h>
#include <iostream>
#include <iomanip>

ClassImp(AliFMDCorrector)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDCorrector::AliFMDCorrector()
  : TNamed(), 
    fRingHistos(),
    fUseSecondaryMap(true),
    fUseVertexBias(false),
    fUseAcceptance(false),
    fUseMergingEfficiency(false),
    fDebug(0)
{
  // Constructor
  DGUARD(fDebug, 3, "Default CTOR of AliFMDCorrector");
}

//____________________________________________________________________
AliFMDCorrector::AliFMDCorrector(const char* title)
  : TNamed("fmdCorrector", title), 
    fRingHistos(), 
    fUseSecondaryMap(true),
    fUseVertexBias(false),
    fUseAcceptance(false),
    fUseMergingEfficiency(false),
    fDebug(0)
{
  // Constructor 
  // 
  // Parameters: 
  //   title   Title
  DGUARD(fDebug, 3, "Named CTOR of AliFMDCorrector: %s", title);
  fRingHistos.SetName(GetName());
}

//____________________________________________________________________
AliFMDCorrector::AliFMDCorrector(const AliFMDCorrector& o)
  : TNamed(o), 
    fRingHistos(), 
    fUseSecondaryMap(o.fUseSecondaryMap),
    fUseVertexBias(o.fUseVertexBias),
    fUseAcceptance(o.fUseAcceptance),
    fUseMergingEfficiency(o.fUseMergingEfficiency),
    fDebug(o.fDebug)
{
  // Copy constructor 
  // 
  // Parameters: 
  //  o  Object to copy from 
  DGUARD(fDebug, 3, "Copy CTOR of AliFMDCorrector");
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDCorrector::~AliFMDCorrector()
{
  // Destructor 
  // 
  // 
  DGUARD(fDebug, 3, "DTOR of AliFMDCorrector");
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDCorrector&
AliFMDCorrector::operator=(const AliFMDCorrector& o)
{
  // Assignment operator 
  // 
  // Parameters:
  //   o   Object to assign from 
  DGUARD(fDebug, 3, "Assignment of AliFMDCorrector");
  if (&o == this) return *this; 
  TNamed::operator=(o);

  fDebug   = o.fDebug;
  fRingHistos.Delete();
  fUseSecondaryMap = o.fUseSecondaryMap;
  fUseVertexBias = o.fUseVertexBias;
  fUseAcceptance = o.fUseAcceptance;
  fUseMergingEfficiency = o.fUseMergingEfficiency;
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
void
AliFMDCorrector::SetupForData(const TAxis&)
{
  //
  // Initialize this object
  //
  // Parameters:
  //    etaAxis Eta axis to use
  //
  DGUARD(fDebug, 1, "Initialization of AliFMDCorrector");
  if (!fUseSecondaryMap)
    AliWarning("Secondary maps not used - BE CAREFUL");
  if (!fUseVertexBias)
    AliWarning("Vertex bias not used");
  if (!fUseAcceptance)
    AliWarning("Acceptance from dead-channels not used");
}

//____________________________________________________________________
AliFMDCorrector::RingHistos*
AliFMDCorrector::GetRingHistos(UShort_t d, Char_t r) const
{
  // 
  // Get the ring histogram container 
  // Parameters:
  //    d Detector
  //    r Ring 
  // 
  // Return:
  //    Ring histogram container 
  //
  Int_t idx = -1;
  switch (d) { 
  case 1: idx = 0; break;
  case 2: idx = 1 + (r == 'I' || r == 'i' ? 0 : 1); break;
  case 3: idx = 3 + (r == 'I' || r == 'i' ? 0 : 1); break;
  }
  if (idx < 0 || idx >= fRingHistos.GetEntries()) return 0;
  
  return static_cast<RingHistos*>(fRingHistos.At(idx));
}
    
//____________________________________________________________________
void
AliFMDCorrector::DivideMap(TH2* num, const TH2* denom,
			   Bool_t alsoUnderOver) const
{
  // 
  // Implement TH1::Divide but 
  // - Assume compatible histograms 
  // - Unless third argument is true, do not divide over/under flow bins
  // 
  if (!num || !denom) return;

  Int_t first = (alsoUnderOver ? 0 : 1);
  Int_t lastX = num->GetNbinsX() + (alsoUnderOver ? 1 : 0);
  Int_t lastY = num->GetNbinsY() + (alsoUnderOver ? 1 : 0);
  
  for (Int_t ix = first; ix <= lastX; ix++) {
    for (Int_t iy = first; iy <= lastY; iy++) { 
      Int_t    bin = num->GetBin(ix,iy);
      Double_t c0  = num->GetBinContent(bin);
      Double_t c1  = denom->GetBinContent(bin);
      if (!c1) { 
	num->SetBinContent(bin,0);
	num->SetBinError(bin, 0);
	continue;
      }
      Double_t w   = c0 / c1;
      Double_t e0  = num->GetBinError(bin);
      Double_t e1  = denom->GetBinError(bin);
      Double_t c12 = c1*c1;
      Double_t e2  = (e0*e0*c1*c1 + e1*e1*c0*c0)/(c12*c12);
      
      num->SetBinContent(bin, w);
      num->SetBinError(bin, TMath::Sqrt(e2));
    }
  }
}
//____________________________________________________________________
Bool_t
AliFMDCorrector::Correct(AliForwardUtil::Histos& hists,
			 UShort_t                vtxbin)
{
  // 
  // Do the calculations 
  // Parameters:
  //    hists    Cache of histograms 
  //    vtxBin   Vertex bin 
  // 
  // Return:
  //    true on successs 
  //
  DGUARD(fDebug, 1, "Correct histograms of AliFMDCorrector");
  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

  UShort_t uvb = vtxbin;
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r  = (q == 0 ? 'I' : 'O');
      TH2D*       h  = hists.Get(d,r);
      RingHistos* rh = GetRingHistos(d,r);

      if (fUseSecondaryMap) {
        TH2D*  bg = fcm.GetSecondaryMap()->GetCorrection(d, r, uvb);
        if (!bg) {
          AliWarning(Form("No secondary correction for FMDM%d%c "
			  "in vertex bin %d", d, r, uvb));
          continue;
        }
        // Divide by primary/total ratio
	DivideMap(h, bg, false);
      }
      if (fUseVertexBias) {
        TH2D*  ef = fcm.GetVertexBias()->GetCorrection(r, uvb);
        if (!ef) {
          AliWarning(Form("No event %s vertex bias correction in vertex bin %d",
                          (r == 'I' || r == 'i' ? "inner" : "outer"), uvb));
          continue;
        }
        // Divide by the event selection efficiency
	DivideMap(h, ef, false);
      }
      if (fUseAcceptance) {
        TH2D*  ac = fcm.GetAcceptance()->GetCorrection(d, r, uvb);
        if (!ac) {
          AliWarning(Form("No acceptance correction for FMD%d%c in "
			  "vertex bin %d", d, r, uvb));
          continue;
        }
	// Fill overflow bin with ones 
	for (Int_t i = 1; i <= h->GetNbinsX(); i++) 
	  h->SetBinContent(i, h->GetNbinsY()+1, 1);

        // Divide by the acceptance correction - 
	DivideMap(h, ac, fcm.GetAcceptance()->HasOverflow());
      }

      if (fUseMergingEfficiency) {
	if (!fcm.GetMergingEfficiency()) { 
	  AliWarning("No merging efficiencies");
	  continue;
	}
	TH1D* sf = fcm.GetMergingEfficiency()->GetCorrection(d,r,uvb);
	if (!fcm.GetMergingEfficiency()->GetCorrection(d,r,uvb)) { 
	  AliWarning(Form("No merging efficiency for FMD%d%c at vertex bin %d",
			  d, r, uvb));
	  continue;
	}
	
      
	for (Int_t ieta = 1; ieta <= h->GetNbinsX(); ieta++) {
	  Float_t c  = sf->GetBinContent(ieta);
	  Float_t ec = sf->GetBinError(ieta);
	  	  
	  for (Int_t iphi = 1; iphi <= h->GetNbinsY(); iphi++) { 
	    if (c == 0) {
	      h->SetBinContent(ieta,iphi,0);
	      h->SetBinError(ieta,iphi,0);
	      continue;
	    }

	    Double_t m  = h->GetBinContent(ieta, iphi) / c;
	    Double_t em = h->GetBinError(ieta, iphi);
	  
	    Double_t e  = TMath::Sqrt(em * em + (m * ec) * (m * ec));
	    
	    h->SetBinContent(ieta,iphi,m);
	    h->SetBinError(ieta,iphi,e);
	  }
	}
      }
      
      rh->fDensity->Add(h);
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDCorrector::Terminate(const TList* dir, TList* output, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // Parameters:
  //    dir     Where the output is stored
  //    nEvents Number of events 
  //
  DGUARD(fDebug, 1, "Scale histograms of AliFMDCorrector");
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TList* out = new TList;
  out->SetName(d->GetName());
  out->SetOwner();

  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  THStack* sums = new THStack("sums", "Sums of ring results");
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Terminate(d, nEvents);
    if (!o->fDensity) { 
      Warning("Terminate", "No density from %s", o->GetName());
      continue;
    }
    TH1D* sum = o->fDensity->ProjectionX(o->GetName(), 1, 
					 o->fDensity->GetNbinsY(),"e");
    sum->Scale(1., "width");
    sum->SetTitle(o->GetName());
    sum->SetDirectory(0);
    sum->SetYTitle("#sum N_{ch,primary}");
    sums->Add(sum);
  }
  out->Add(sums);
  output->Add(out);
}
//____________________________________________________________________
void
AliFMDCorrector::CreateOutputObjects(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  DGUARD(fDebug, 1, "Define output of AliFMDCorrector");
  TList* d = new TList;
  d->SetOwner();
  d->SetName(GetName());
  dir->Add(d);

  d->Add(AliForwardUtil::MakeParameter("secondary", fUseSecondaryMap));
  d->Add(AliForwardUtil::MakeParameter("vertexBias", fUseVertexBias));
  d->Add(AliForwardUtil::MakeParameter("acceptance", fUseAcceptance));
  d->Add(AliForwardUtil::MakeParameter("merging", fUseMergingEfficiency));
  

  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->CreateOutputObjects(d);
  }
}

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliFMDCorrector::Print(Option_t* /* option */) const
{
  // 
  // Print information
  // Parameters:
  //    option Not used 
  //
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFB("Use secondary maps",		fUseSecondaryMap );
  PFB("Use vertex bias",		fUseVertexBias );
  PFB("Use acceptance",			fUseAcceptance );
  PFB("Use merging efficiency",		fUseMergingEfficiency);
  gROOT->DecreaseDirLevel();  
}

//====================================================================
AliFMDCorrector::RingHistos::RingHistos()
  : AliForwardUtil::RingHistos(), 
    fDensity(0)
{
  // Constructor 
  //
  // 
}

//____________________________________________________________________
AliFMDCorrector::RingHistos::RingHistos(UShort_t d, Char_t r)
  : AliForwardUtil::RingHistos(d,r), 
    fDensity(0)
{
  // 
  // Constructor
  // Parameters:
  //    d detector
  //    r ring 
  //
  fDensity = new TH2D("primaryDensity", 
		      "#sum N_{ch,primary}/(#Delta#eta#Delta#phi)", 
		      200, -4, 6, (r == 'I' || r == 'i' ? 20 : 40), 
		      0, 2*TMath::Pi());
  fDensity->SetDirectory(0);
  fDensity->SetMarkerColor(Color());
  fDensity->Sumw2();
  fDensity->SetXTitle("#eta");
  fDensity->SetYTitle("#phi [radians]");
  fDensity->SetZTitle("Primary N_{ch} density");
}
//____________________________________________________________________
AliFMDCorrector::RingHistos::RingHistos(const RingHistos& o)
  : AliForwardUtil::RingHistos(o), 
    fDensity(o.fDensity)
{
  // 
  // Copy constructor 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDCorrector::RingHistos&
AliFMDCorrector::RingHistos::operator=(const RingHistos& o)
{
  // 
  // Assignment operator 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this 
  //
  if (&o == this) return *this; 
  AliForwardUtil::RingHistos::operator=(o);
  
  if (fDensity) delete fDensity;
  
  fDensity = static_cast<TH2D*>(o.fDensity->Clone());
  
  return *this;
}
//____________________________________________________________________
AliFMDCorrector::RingHistos::~RingHistos()
{
  // 
  // Destructor 
  //
  // if (fDensity) delete fDensity;
}

//____________________________________________________________________
void
AliFMDCorrector::RingHistos::CreateOutputObjects(TList* dir)
{
  // 
  // Make output 
  // Parameters:
  //    dir Where to put it 
  //
  TList* d = DefineOutputList(dir);
  d->Add(fDensity);
}

//____________________________________________________________________
void
AliFMDCorrector::RingHistos::Terminate(TList* dir, Int_t nEvents)
{ 
  // 
  // Scale the histograms to the total number of events 
  // Parameters:
  //    dir     where the output is stored
  //    nEvents Number of events 
  //
  TList* l = GetOutputList(dir);
  if (!l) return; 

  TH2D* density = static_cast<TH2D*>(GetOutputHist(l,"primaryDensity"));
  if (density) density->Scale(1./nEvents);
  fDensity = density;
}

//____________________________________________________________________
//
// EOF
//
	  


