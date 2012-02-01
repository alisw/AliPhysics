// 
// This class calculates the exclusive charged particle density
// in each for the 5 FMD rings. 
//
// Input:
//   - 5 RingHistos objects - each with a number of vertex dependent 
//     2D histograms of the inclusive charge particle density 
//
// Output:
//   - 5 RingHistos objects - each with a number of vertex dependent 
//     2D histograms of the exclusive charge particle density 
// 
// Corrections used: 
//   - AliFMDCorrSecondaryMap;
//   - AliFMDCorrVertexBias
//   - AliFMDCorrMergingEfficiency
//
#include "AliFMDMCCorrector.h"
#include <AliESDFMD.h>
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrVertexBias.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TROOT.h>
#include <TProfile2D.h>
#include <iostream>
ClassImp(AliFMDMCCorrector)
#if 0
; // For Emacs
#endif 


//____________________________________________________________________
AliFMDMCCorrector::~AliFMDMCCorrector()
{
  // 
  // Destructor 
  //
  if (fComps) fComps->Clear();
  if (fFMD1i) delete fFMD1i;
  if (fFMD2i) delete fFMD2i;
  if (fFMD2o) delete fFMD2o;
  if (fFMD3i) delete fFMD3i;
  if (fFMD3o) delete fFMD3o;
}

//____________________________________________________________________
AliFMDMCCorrector&
AliFMDMCCorrector::operator=(const AliFMDMCCorrector& o)
{
  // 
  // Assignement operator
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object
  //
  if (&o == this) return *this; 
  AliFMDCorrector::operator=(o);
  fSecondaryForMC = o.fSecondaryForMC;
  return *this;
}

//____________________________________________________________________
Bool_t
AliFMDMCCorrector::CorrectMC(AliForwardUtil::Histos& hists,
			       UShort_t                vtxbin)
{
  // 
  // Do the calculations 
  // 
  // Parameters:
  //    hists    Cache of histograms 
  //    vtxBin   Vertex bin 
  // 
  // Return:
  //    true on successs 
  //
  if ((!fUseSecondaryMap || !fSecondaryForMC) && fUseVertexBias) 
    return kTRUE;

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

  UShort_t uvb = vtxbin;
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r  = (q == 0 ? 'I' : 'O');
      TH2D*       h  = hists.Get(d,r);


      if (fUseSecondaryMap && fSecondaryForMC) {
        TH2D*       bg = fcm.GetSecondaryMap()->GetCorrection(d,r,uvb);
        if (!bg) {
          AliWarning(Form("No secondary correction for FMDM%d%c in vertex bin %d",
              d, r, uvb));
          continue;
        }
        // Divide by primary/total ratio
        h->Divide(bg);
      }
      if (fUseVertexBias) {
        TH2D*       ef = fcm.GetVertexBias()->GetCorrection(r, uvb);
	if (!ef) {
          AliWarning(Form("No event vertex bias correction in vertex bin %d",
			uvb));
          continue;
        }
        // Divide by the event selection efficiency
        h->Divide(ef);
      }
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDMCCorrector::Init(const TAxis& eAxis)
{
  // 
  // Initialize this object 
  // 
  // Parameters:
  //    etaAxis Eta axis to use 
  //
  AliFMDCorrector::Init(eAxis);

  fFMD1i = Make(1,'I',eAxis);
  fFMD2i = Make(2,'I',eAxis);
  fFMD2o = Make(2,'O',eAxis);
  fFMD3i = Make(3,'I',eAxis);
  fFMD3o = Make(3,'O',eAxis);

  fComps->Add(fFMD1i);
  fComps->Add(fFMD2i);
  fComps->Add(fFMD2o);
  fComps->Add(fFMD3i);
  fComps->Add(fFMD3o);
}

//____________________________________________________________________
TProfile2D*
AliFMDMCCorrector::Make(UShort_t d, Char_t r, 
				const TAxis& axis) const
{
  // 
  // MAke comparison profiles
  // 
  // Parameters:
  //    d     Detector 
  //    r     Ring 
  //    axis  Eta axis 
  // 
  // Return:
  //    Newly allocated profile object
  //
  TProfile2D* ret = new TProfile2D(Form("FMD%d%c_esd_vs_mc", d, r),
				   Form("ESD/MC signal for FMD%d%c", d, r),
				   axis.GetNbins(), 
				   axis.GetXmin(),
				   axis.GetXmax(), 
				   (r == 'I' || r == 'i') ? 20 : 40,
				   0, 2*TMath::Pi());
  ret->GetXaxis()->SetTitle("#eta");
  ret->GetYaxis()->SetTitle("#varphi [degrees]");
  ret->GetZaxis()->SetTitle("#LT primary density ESD/MC#GT");
  ret->SetDirectory(0);
  return ret;
}
//____________________________________________________________________
void
AliFMDMCCorrector::Fill(UShort_t d, Char_t r, TH2* esd, TH2* mc)
{
  // 
  // Fill comparison profiles
  // 
  // Parameters:
  //    d    Detector 
  //    r    Ring 
  //    esd  ESD histogram
  //    mc   MC histogram
  //
  if (!esd || !mc) return;
  TProfile2D* p = 0;
  switch (d) { 
  case 1:  p = fFMD1i;                                   break;
  case 2:  p = (r == 'I' || r == 'i' ? fFMD2i : fFMD2o); break;
  case 3:  p = (r == 'I' || r == 'i' ? fFMD3i : fFMD3o); break;
  }
  if (!p) return;

  for (Int_t iEta = 1; iEta <= esd->GetNbinsX(); iEta++) { 
    Double_t eta = esd->GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= esd->GetNbinsY(); iPhi++) { 
      Double_t phi  = esd->GetYaxis()->GetBinCenter(iPhi);
      Double_t mEsd = esd->GetBinContent(iEta,iPhi);
      Double_t mMc  = mc->GetBinContent(iEta,iPhi);
      
      p->Fill(eta, phi, (mMc > 0 ? mEsd / mMc : 0));
    }
  }
}

//____________________________________________________________________
Bool_t
AliFMDMCCorrector::CompareResults(AliForwardUtil::Histos& esd,
					  AliForwardUtil::Histos& mc)
{
  // 
  // Compare the result of analysing the ESD for 
  // the inclusive charged particle density to analysing 
  // MC truth 
  // 
  // Parameters:
  //    esd 
  //    mc 
  // 
  // Return:
  //   true 
  //
  Fill(1, 'I', esd.Get(1,'I'), mc.Get(1,'I'));
  Fill(2, 'I', esd.Get(2,'I'), mc.Get(2,'I'));
  Fill(2, 'O', esd.Get(2,'O'), mc.Get(2,'O'));
  Fill(3, 'I', esd.Get(3,'I'), mc.Get(3,'I'));
  Fill(3, 'O', esd.Get(3,'O'), mc.Get(3,'O'));

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDMCCorrector::DefineOutput(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  AliFMDCorrector::DefineOutput(dir);
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));

  fComps = new TList;
  fComps->SetName("esd_mc_comparison");
  d->Add(fComps);
}
//____________________________________________________________________
void
AliFMDMCCorrector::Print(Option_t* option) const
{
  // 
  // Print information
  // Parameters:
  //    option Not used 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  AliFMDCorrector::Print(option);
  std::cout << std::boolalpha
            << ind << " Use sec. map on MC:     " << fSecondaryForMC 
            << std::noboolalpha << std::endl;
}

//____________________________________________________________________
//
// EOF
//
	  


