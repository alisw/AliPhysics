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
#include "AliLog.h"
#include <TH2D.h>
#include <TROOT.h>
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
    fUseVertexBias(true),
    fUseAcceptance(true),
    fUseMergingEfficiency(true),
    fDebug(0)
{
  // Constructor
}

//____________________________________________________________________
AliFMDCorrector::AliFMDCorrector(const char* title)
  : TNamed("fmdCorrector", title), 
    fRingHistos(), 
    fUseSecondaryMap(true),
    fUseVertexBias(true),
    fUseAcceptance(true),
    fUseMergingEfficiency(true),
    fDebug(0)
{
  // Constructor 
  // 
  // Parameters: 
  //   title   Title
  fRingHistos.SetName(GetName());
  fRingHistos.Add(new RingHistos(1, 'I'));
  fRingHistos.Add(new RingHistos(2, 'I'));
  fRingHistos.Add(new RingHistos(2, 'O'));
  fRingHistos.Add(new RingHistos(3, 'I'));
  fRingHistos.Add(new RingHistos(3, 'O'));
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
AliFMDCorrector::Init(const TAxis&)
{
  //
  // Initialize this object
  //
  // Parameters:
  //    etaAxis Eta axis to use
  //
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
          AliWarning(Form("No secondary correction for FMDM%d%c in vertex bin %d",
                          d, r, uvb));
          continue;
        }
        // Divide by primary/total ratio
        h->Divide(bg);
      }
      if (fUseVertexBias) {
        TH2D*  ef = fcm.GetVertexBias()->GetCorrection(r, uvb);
        if (!ef) {
          AliWarning(Form("No event %s vertex bias correction in vertex bin %d",
                          (r == 'I' || r == 'i' ? "inner" : "outer"), uvb));
          continue;
        }
        // Divide by the event selection efficiency
        h->Divide(ef);
      }
      if (fUseAcceptance) {
        TH2D*  ac = fcm.GetAcceptance()->GetCorrection(d, r, uvb);
        if (!ac) {
          AliWarning(Form("No acceptance correction for FMD%d%c in vertex bin %d",
			d, r, uvb));
          continue;
        }
        // Divide by the acceptance correction
        h->Divide(ac);
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
	  
	  if (c == 0) continue;
	  
	  for (Int_t iphi = 1; iphi <= h->GetNbinsY(); iphi++) { 
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
AliFMDCorrector::ScaleHistograms(TList* dir, Int_t nEvents)
{
  // 
  // Scale the histograms to the total number of events 
  // Parameters:
  //    dir     Where the output is stored
  //    nEvents Number of events 
  //
  if (nEvents <= 0) return;
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));
  if (!d) return;

  TIter    next(&fRingHistos);
  RingHistos* o = 0;
  while ((o = static_cast<RingHistos*>(next())))
    o->ScaleHistograms(d, nEvents);
}
//____________________________________________________________________
void
AliFMDCorrector::DefineOutput(TList* dir)
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

//____________________________________________________________________
void
AliFMDCorrector::Print(Option_t* /* option */) const
{
  // 
  // Print information
  // Parameters:
  //    option Not used 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << "AliFMDCorrector: " << GetName() <<  "\n"
            << std::boolalpha
            << ind << " Use secondary maps:     " << fUseSecondaryMap << "\n"
            << ind << " Use vertex bias:        " << fUseVertexBias << "\n"
            << ind << " Use acceptance:         " << fUseAcceptance << "\n"
            << ind << " Use merging efficiency: " << fUseMergingEfficiency
            << std::noboolalpha << std::endl;
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
  if (fDensity) delete fDensity;
}

//____________________________________________________________________
void
AliFMDCorrector::RingHistos::Output(TList* dir)
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
AliFMDCorrector::RingHistos::ScaleHistograms(TList* dir, Int_t nEvents)
{ 
  // 
  // Scale the histograms to the total number of events 
  // Parameters:
  //    dir     where the output is stored
  //    nEvents Number of events 
  //
  TList* l = GetOutputList(dir);
  if (!l) return; 

  TH1* density = GetOutputHist(l,Form("%s_Primary_Density", fName.Data()));
  if (density) density->Scale(1./nEvents);
}

//____________________________________________________________________
//
// EOF
//
	  


