#include "AliFMDCorrections.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliFMDAnaParameters.h"
#include "AliLog.h"
#include <TH2D.h>

ClassImp(AliFMDCorrections)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDCorrections::AliFMDCorrections()
  : TNamed(), 
    fRingHistos(),
    fMultCut(0.3)
{}

//____________________________________________________________________
AliFMDCorrections::AliFMDCorrections(const char* title)
  : TNamed("fmdCorrections", title), 
    fRingHistos(), 
    fMultCut(0.3)
{
  fRingHistos.SetName(GetName());
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
}

//____________________________________________________________________
AliFMDCorrections::AliFMDCorrections(const AliFMDCorrections& o)
  : TNamed(o), 
    fRingHistos(), 
    fMultCut(o.fMultCut)
{
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
}

//____________________________________________________________________
AliFMDCorrections::~AliFMDCorrections()
{
  fRingHistos.Delete();
}

//____________________________________________________________________
AliFMDCorrections&
AliFMDCorrections::operator=(const AliFMDCorrections& o)
{
  SetName(o.GetName());
  SetTitle(o.GetTitle());

  fMultCut = o.fMultCut;

  fRingHistos.Delete();
  TIter    next(&o.fRingHistos);
  TObject* obj = 0;
  while ((obj = next())) fRingHistos.Add(obj);
  
  return *this;
}

//____________________________________________________________________
AliFMDCorrections::RingHistos*
AliFMDCorrections::GetRingHistos(UShort_t d, Char_t r) const
{
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
Bool_t
AliFMDCorrections::Correct(AliForwardUtil::Histos& hists,
			   Int_t                   vtxbin)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();

  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      TH2D*       h = hists.Get(d,r);
      RingHistos* rh= GetRingHistos(d,r);
      TH2F*       bg= pars->GetBackgroundCorrection(d, r, vtxbin);
      TH2F*       ef= pars->GetEventSelectionEfficiency("INEL",vtxbin,r);
      if (!bg) { 
	AliWarning(Form("No secondary correction for FMDM%d%c in vertex bin %d",
			d, r, vtxbin));
	continue;
      }
      if (!ef) { 
	AliWarning(Form("No event %s selection efficiency in vertex bin %d",
			"INEL", vtxbin));
	continue;
      }

      // Divide by primary/total ratio
      h->Divide(bg);
      
      // Divide by the event selection efficiency 
      h->Divide(ef);

      if(!pars->SharingEffPresent()) { 
	AliWarning("No sharing efficiencies");
	continue;
      }

      TH1F* sf = pars->GetSharingEfficiencyTrVtx(d,r,vtxbin); 
      
      for (Int_t ieta = 1; ieta <= h->GetNbinsX(); ieta++) {
	Float_t c  = sf->GetBinContent(ieta);
	Float_t ec = sf->GetBinError(ieta);
	
	if (c == 0) continue;

	for (Int_t iphi = 1; iphi <= h->GetNbinsY(); iphi++) { 
	  Double_t m  = h->GetBinContent(ieta, iphi) / c;
	  Double_t em = h->GetBinError(ieta, iphi);
	  
	  Double_t e  = TMath::Sqrt(em * em + (m * ec) * (m * ec));
	  
	  h->SetBinContent(ieta,iphi,m);
	  h->SetBinError(ieta,iphi,e);
	}
      }

      rh->fDensity->Add(h);
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDCorrections::ScaleHistograms(Int_t nEvents)
{
  if (nEvents <= 0) return;

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->fDensity->Scale(1. / nEvents);
  }
}

//____________________________________________________________________
void
AliFMDCorrections::Output(TList* dir)
{
  TList* d = new TList;
  d->SetName(GetName());
  dir->Add(d);
  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next()))) {
    o->Output(d);
  }
}

//====================================================================
AliFMDCorrections::RingHistos::RingHistos()
  : fDet(0),
    fRing('\0'),
    fDensity(0)
{}

//____________________________________________________________________
AliFMDCorrections::RingHistos::RingHistos(UShort_t d, Char_t r)
  : fDet(d), 
    fRing(r),
    fDensity(0)
{
  fDensity = new TH2D(Form("FMD%d%c_Primary_Density", d, r), 
		      Form("in FMD%d%c", d, r), 
		      200, -4, 6, (r == 'I' || r == 'i' ? 20 : 40), 
		      0, 2*TMath::Pi());
  fDensity->SetDirectory(0);
  fDensity->SetXTitle("#eta");
  fDensity->SetYTitle("#phi [radians]");
  fDensity->SetZTitle("Primary N_{ch} density");
}
//____________________________________________________________________
AliFMDCorrections::RingHistos::RingHistos(const RingHistos& o)
  : TObject(o), 
    fDet(o.fDet), 
    fRing(o.fRing),
    fDensity(o.fDensity)
{}

//____________________________________________________________________
AliFMDCorrections::RingHistos&
AliFMDCorrections::RingHistos::operator=(const RingHistos& o)
{
  fDet  = o.fDet;
  fRing = o.fRing;
  
  if (fDensity) delete fDensity;
  
  fDensity = static_cast<TH2D*>(o.fDensity->Clone());
  
  return *this;
}
//____________________________________________________________________
AliFMDCorrections::RingHistos::~RingHistos()
{
  if (fDensity) delete fDensity;
}

//____________________________________________________________________
void
AliFMDCorrections::RingHistos::Output(TList* dir)
{
  TList* d = new TList;
  d->SetName(Form("FMD%d%c", fDet, fRing)); 
  d->Add(fDensity);
  dir->Add(d);
}

//____________________________________________________________________
//
// EOF
//
	  


