// 
// This class calculates the inclusive charged particle density
// in each for the 5 FMD rings based on the MC truth.
//
// Input:
//   - AliMCEvent  MC truth event infromation
//
// Output:
//   - None
//
// Corrections used: 
//   - None
//
//
#include "AliFMDMCDensityCalculator.h"
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDStripIndex.h"
#include "AliMCEvent.h"
#include "AliESDFMD.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TProfile2D.h>

ClassImp(AliFMDMCDensityCalculator)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDMCDensityCalculator::~AliFMDMCDensityCalculator()
{
  // 
  // Destructor 
  //
  if (fComps)  fComps->Clear();
  if (fFMD1i)  delete fFMD1i;
  if (fFMD2i)  delete fFMD2i;
  if (fFMD2o)  delete fFMD2o;
  if (fFMD3i)  delete fFMD3i;
  if (fFMD3o)  delete fFMD3o;
  if (fFMD1iC) delete fFMD1iC;
  if (fFMD2iC) delete fFMD2iC;
  if (fFMD2oC) delete fFMD2oC;
  if (fFMD3iC) delete fFMD3iC;
  if (fFMD3oC) delete fFMD3oC;
}

//____________________________________________________________________
AliFMDMCDensityCalculator&
AliFMDMCDensityCalculator::operator=(const AliFMDMCDensityCalculator& o)
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
  AliFMDDensityCalculator::operator=(o);
  return *this;
}


//____________________________________________________________________
void
AliFMDMCDensityCalculator::Init(const TAxis& eAxis)
{
  // 
  // Initialize this object 
  // 
  // Parameters:
  //    etaAxis Eta axis to use 
  //
  AliFMDDensityCalculator::Init(eAxis);
  fFMD1i  = Make(1,'I',eAxis);
  fFMD2i  = Make(2,'I',eAxis);
  fFMD2o  = Make(2,'O',eAxis);
  fFMD3i  = Make(3,'I',eAxis);
  fFMD3o  = Make(3,'O',eAxis); 
  fFMD1iC = Make(1,'I');
  fFMD2iC = Make(2,'I');
  fFMD2oC = Make(2,'O');
  fFMD3iC = Make(3,'I');
  fFMD3oC = Make(3,'O');
 
  fComps->Add(fFMD1i);
  fComps->Add(fFMD2i);
  fComps->Add(fFMD2o);
  fComps->Add(fFMD3i);
  fComps->Add(fFMD3o);
  fComps->Add(fFMD1iC);
  fComps->Add(fFMD2iC);
  fComps->Add(fFMD2oC);
  fComps->Add(fFMD3iC);
  fComps->Add(fFMD3oC);
}
    

//____________________________________________________________________
Bool_t
AliFMDMCDensityCalculator::CalculateMC(const AliESDFMD&        fmd,
				       AliForwardUtil::Histos& hists)
{
  // 
  // Calculate the charged particle density from the MC track references. 
  // 
  // Parameters:
  //    fmd    Forward MC event
  //    hists  Histograms to fill
  //    vz     Interaction z coordinate @f$ v_z@f$
  //    vtxBin bin corresponding to @f$ v_z@f$
  // 
  // Return:
  //    true on success
  //
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = hists.Get(d,r);

      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  Float_t mult = fmd.Multiplicity(d,r,s,t);
	  
	  if (mult == 0 || mult > 20) continue;

	  Float_t phi = fmd.Phi(d,r,s,t) / 180 * TMath::Pi();
	  Float_t eta = fmd.Eta(d,r,s,t);

	  Float_t corr = Correction(d, r, s, t, 0, eta, false);
	  Float_t sig  = (corr <= 0 ? 0 : mult / corr);
	  h->Fill(eta,phi,sig);
	}
      }
    }
  }
  return kTRUE;
}

//____________________________________________________________________
TProfile2D*
AliFMDMCDensityCalculator::Make(UShort_t d, Char_t r, 
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
  ret->GetZaxis()->SetTitle("#LT incusive density ESD/MC#GT");
  ret->SetDirectory(0);
  return ret;
}
//____________________________________________________________________
TH2D*
AliFMDMCDensityCalculator::Make(UShort_t d, Char_t r) const
{
  // 
  // MAke comparison profiles
  // 
  // Parameters:
  //    d     Detector 
  //    r     Ring 
  // 
  // Return:
  //    Newly allocated profile object
  //
  TH2D* ret = new TH2D(Form("FMD%d%c_corr_mc_esd", d, r),
		       Form("ESD-MC correlation for FMD%d%c", d, r),
		       200, 0, 20, 200, 0, 20);
  ret->GetXaxis()->SetTitle("m_{incl} (MC)");
  ret->GetYaxis()->SetTitle("#Delta/#Delta_{mp} (ESD)");
  ret->GetZaxis()->SetTitle("Correlation");
  ret->SetDirectory(0);
  return ret;
}
//____________________________________________________________________
void
AliFMDMCDensityCalculator::Fill(UShort_t d, Char_t r, TH2* esd, TH2* mc)
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
  TH2D*       h = 0;
  switch (d) { 
  case 1:  
    p = fFMD1i;                                   
    h = fFMD1iC;
    break;
  case 2:  
    p = (r == 'I' || r == 'i' ? fFMD2i  : fFMD2o); 
    h = (r == 'I' || r == 'i' ? fFMD2iC : fFMD2oC); 
    break;
  case 3:  
    p = (r == 'I' || r == 'i' ? fFMD3i  : fFMD3o); 
    h = (r == 'I' || r == 'i' ? fFMD3iC : fFMD3oC); 
    break;
  }
  if (!p) return;

  for (Int_t iEta = 1; iEta <= esd->GetNbinsX(); iEta++) { 
    Double_t eta = esd->GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= esd->GetNbinsY(); iPhi++) { 
      Double_t phi  = esd->GetYaxis()->GetBinCenter(iPhi);
      Double_t mEsd = esd->GetBinContent(iEta,iPhi);
      Double_t mMc  = mc->GetBinContent(iEta,iPhi);
      
      p->Fill(eta, phi, (mMc > 0 ? mEsd / mMc : 0));
      h->Fill(mMc,mEsd);
    }
  }
}

//____________________________________________________________________
Bool_t
AliFMDMCDensityCalculator::CompareResults(AliForwardUtil::Histos& esd,
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
AliFMDMCDensityCalculator::DefineOutput(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  AliFMDDensityCalculator::DefineOutput(dir);
  TList* d = static_cast<TList*>(dir->FindObject(GetName()));

  fComps = new TList;
  fComps->SetName("esd_mc_comparison");
  d->Add(fComps);

}
//____________________________________________________________________
//
// EOF
//
	  


