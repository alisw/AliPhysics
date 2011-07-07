Color_t Color(UShort_t d, Char_t r ) const 
{ 
  return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	  + ((r == 'I' || r == 'i') ? 2 : -2));
}

void
RunCopyELossFit(UShort_t sys, UShort_t cms, Short_t field, Bool_t mc=false,
		const Char_t* path=0)
{
  RunCopyELossFit(sys == 1 ? "pp" : "PbPb", 
		  cms, 
		  field, 
		  mc,
		  path);
}

/** 
 * 
 * @param sys       Collision system 
 * @param cms       Center of mass energy per nucleon in GeV
 * @param field     Magnetic field 
 * 
 * @ingroup pwg2_forward_analysis_scripts
 */
void
RunCopyELossFit(const char* sys, UShort_t cms, Short_t field, bool mc=false,
		const char* path=0)
{
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
  gSystem->Load("libPWG2forward.so");

  AliFMDAnaParameters* p = AliFMDAnaParameters::Instance();
  p->SetRealData(!mc);
  p->SetEnergy(Float_t(cms));
  p->SetMagField(Float_t(field));
  p->SetCollisionSystem(sys);
  if (path) {
    p->SetBackgroundPath(path);
    p->SetEnergyPath(path);
    p->SetEventSelectionPath(path);
    p->SetSharingEfficiencyPath(path);
  }
  p->Init(true, AliFMDAnaParameters::kBackgroundCorrection|
	  AliFMDAnaParameters::kEnergyDistributions);
  
  Int_t    nVtx   = p->GetNvtxBins();
  Double_t minVtx = -p->GetVtxCutZ();
  Double_t maxVtx = -p->GetVtxCutZ();
  Int_t    nEta   = p->GetNetaBins();
  Double_t minEta = p->GetEtaMin();
  Double_t maxEta = p->GetEtaMax();
  
  TAxis vtxAxis(nVtx, minVtx, maxVtx);
  TAxis etaAxis(nEta, minEta, maxEta);
  AliFMDCorrELossFit* m = new AliFMDCorrELossFit;
  m->SetEtaAxis(etaAxis);
  
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
      
      for (UShort_t b = 1; b < nEta; b++) { 
	Double_t eta = etaAxis.GetBinCenter(b);
	Info("RunCopyELossFit", "FMD%d%c, bin %3d, eta %+8.4f", d, r, b, eta);
	if (eta < minEta || eta > maxEta) continue;
	TH1F* oldhist = p->GetEnergyDistribution(d,r,eta);
	if (!oldhist) {
	  Warning("RunCopyELossFit",
		  "Didn't find energy distribution for FMD%d%c, eta bin %3d", 
		  d, r, b);
	  continue;
	}
	TF1*  oldfunc = oldhist->GetFunction("FMDfitFunc");
	if (!oldfunc) {
	  // Warning("RunCopyELossFit",
	  // "Didn't find energy fit for FMD%d%c, eta bin %3d", 
	  // d, r, b);
	  continue;
	}
	Double_t chi2   = oldfunc->GetChisquare();
	UShort_t nu     = oldfunc->GetNDF();
	Double_t c      = oldfunc->GetParameter(0);
	Double_t delta  = oldfunc->GetParameter(1);
	Double_t xi     = oldfunc->GetParameter(2);
	Double_t a2 = 0, a3 = 0, ea2 = 0, ea3 = 0;
	if (oldfunc->GetNpar() >= 4) {
	  a2     = oldfunc->GetParameter(3);
	  ea2    = oldfunc->GetParError(3);
	}
	if (oldfunc->GetNpar() >= 5) { 
	  a3     = oldfunc->GetParameter(4);
	  ea3    = oldfunc->GetParError(4);
	}
	Int_t    n      = (a3 < 1e-5 ? (a2 < 1e-5 ? 1 : 2) : 3);
	Double_t ec     = oldfunc->GetParError(0);
	Double_t edelta = oldfunc->GetParError(1);
	Double_t exi    = oldfunc->GetParError(2);
	delta           = delta - xi * -0.22278298;
	edelta          = TMath::Sqrt(TMath::Power(edelta,2) - 
				      TMath::Power(0.22278298*exi,2));
	Info("RunCopyELossFit", "FMD%d%c etabin=%3d: n=%d\n" 
	     "  C=%f+/-%f, Delta=%f+/-%f, xi=%f+/-%f, a2=%f+/-%f, a3=%f+/-%f",
	     d, r, b, n, c, ec, delta, edelta, xi, exi, a2, ea2, a3, ea3);
	Double_t a[]  = { a2, a3, -9999 };
	Double_t ea[] = { ea2, ea3, -9999 };

	AliFMDCorrELossFit::ELossFit* fit = new 
	  AliFMDCorrELossFit::ELossFit(0, n, 
				       chi2,   nu, 
				       c,      ec, 
				       delta,  edelta, 
				       xi,     exi, 
				       0,      0, 
				       0,      0, 
				       a,      ea);
	fit->CalculateQuality();
	fit->fBin = b;
	fit->fDet = d;
	fit->fRing = r;
	fit->Print();
	
	Info("RunCopyELossFit", "Setting fit FMD%d%c etabin=%3d: %p", 
	     d, r, b, fit);
	m->SetFit(d, r, b, fit);
      }
    }
  }
  UShort_t isys = AliForwardUtil::ParseCollisionSystem(sys);
  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString fname(mgr.GetFileName(AliForwardCorrectionManager::kELossFits,
				isys, cms, field, mc));
  TFile* output = TFile::Open(fname.Data(), "RECREATE");
  if (!output) { 
    Warning("RunCopyELossFit", "Failed to open output file %s", fname.Data());
    return kFALSE;
  }
  m->Write(mgr.GetObjectName(AliForwardCorrectionManager::kELossFits));
  output->Write();
  output->Close();
  Info("RunCopyELossFit", 
       "File %s created.  It should be copied to %s and stored in SVN",
       fname.Data(), mgr.GetFileDir(AliForwardCorrectionManager::kELossFits));
  
}
//
// EOF
//
