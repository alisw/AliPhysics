Char_t *itssafilename  = "SPECTRA_ITSsa.root";
Char_t *itstpcfilename = "SPECTRA_ITSTPC_V2.root";
Char_t *tpctoffilename = "SPECTRA_TPCTOF.root";
Char_t *toffilename    = "SPECTRA_TOF.root";
Char_t *hydrofilename  = "HydroSpectra.root";

Char_t *itssaratiofilename  = "SPECTRA_ITSsa.root";
Char_t *itstpcratiofilename = "RATIOS_ITSTPC_V2.root";
Char_t *tpctofratiofilenameA = "RATIOSa_TPCTOF.root";
Char_t *tpctofratiofilenameB = "RATIOSb_TPCTOF.root";
Char_t *tofratiofilename    = "RATIOS_TOF.root";

const Int_t NptBins = 52;
Double_t ptBin[NptBins + 1] = {0.05, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

Double_t ptMin[5] = {0., 0., 0., 0., 0.};
Double_t ptMax[5] = {0., 0., 3.0, 3., 5.};

Char_t *centName[10] = {
  "0-5%",
  "5-10%",
  "10-20%",
  "20-30%",
  "30-40%",
  "40-50%",
  "50-60%",
  "60-70%",
  "70-80%",
  "80-90%"
};

Char_t *chargeName[2] = {"plus", "minus"};

Char_t *partChargeName[5][2] = {
  "", "",
  "", "",
  "#pi^{+}", "#pi^{-}",
  "K^{+}", "K^{-}",
  "p", "#bar{p}"
};

enum EPart_t {
  kPi, kPiPlus, kPiMinus,
  kKa, kKaPlus, kKaMinus,
  kPr, kPrPlus, kPrMinus
};

Char_t *ratioName[9] = {
  "pion", "pion_positive", "pion_negative",
  "kaon", "kaon_positive", "kaon_negative",
  "proton", "proton_positive", "proton_negative"
};

//______________________________________________________________

LoadLibraries()
{

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  //  gSystem->Load("libPWG2");
  gSystem->Load("libPWG2spectra");

}

//______________________________________________________________
//______________________________________________________________

Char_t *ITSsaPartName[5] = {"", "", "Pion", "Kaon", "Proton"};
Char_t *ITSsaChargeName[2] = {"Pos", "Neg"};
TH1D *
GetITSsaSpectrum(TFile *file, Int_t part, Int_t charge, Int_t cent, Bool_t cutSpectrum = kTRUE, Bool_t addSystematicError = kTRUE)
{
  /* pt limits for combined spectra */
  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0.1, 0.2, 0.3};
  Double_t ptMax[AliPID::kSPECIES] = {0., 0., 0.6, 0.5, 0.6};

  TList *list = (TList *)file->Get("output");
  TH1D *hin = (TH1D *)list->FindObject(Form("h_%s_%s_cen_%d", ITSsaPartName[part], ITSsaChargeName[charge], cent));
  if (!hin) return NULL;  

  /* get systematics */
  TFile *fsys = TFile::Open("SPECTRASYS_ITSsa.root");
  TH1 *hsys = fsys->Get(Form("hSystTot%s%s", ITSsaChargeName[charge], ITSsaPartName[part]));
			
  TH1D *h = new TH1D(Form("hITSsa_cent%d_%s_%s", cent, AliPID::ParticleName(part), chargeName[charge]), "ITSsa", NptBins, ptBin);
  Double_t pt, width, value, error, sys;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin[part] || pt > ptMax[part])) continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    /*** TEMP ADD SYS ***/
    if (addSystematicError) {
      sys = hsys->GetBinContent(bin) * value;
      error = TMath::Sqrt(error * error + sys * sys);
    }
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }

  h->SetTitle("ITSsa");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(1);
  h->SetFillStyle(0);
  h->SetFillColor(0);

  return h;
}

//______________________________________________________________

TH1D *
GetITSsaRatio(TFile *file, Int_t num, Int_t den, Int_t cent, Bool_t cutSpectrum = kTRUE, Bool_t addSystematicError = kTRUE)
{
  /* pt limits for combined spectra */
  //  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0.1, 0.2, 0.3};
  //  Double_t ptMax[AliPID::kSPECIES] = {0., 0., 0.6, 0.5, 0.6};

  TH1 *hnum, *hden;
  Double_t ptMin = 0., ptMax = 10.;
  switch (num) {
  case kPiMinus:
    ptMin = TMath::Min(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kPion, 1, cent, kFALSE, kFALSE);
    break;
  case kPiPlus:
    ptMin = TMath::Min(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kPion, 0, cent, kFALSE, kFALSE);
    break;
  case kPi:
    ptMin = TMath::Min(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kPion, 1, cent, kFALSE, kFALSE);
    hnum->Add(GetITSsaSpectrum(file, AliPID::kPion, 0, cent, kFALSE, kFALSE));
    break;
  case kKaMinus:
    ptMin = TMath::Min(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hnum = GetITSsaSpectrum(file, AliPID::kKaon, 1, cent, kFALSE, kFALSE);
    break;
  case kKaPlus:
    ptMin = TMath::Min(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hnum = GetITSsaSpectrum(file, AliPID::kKaon, 0, cent, kFALSE, kFALSE);
    break;
  case kKa:
    ptMin = TMath::Min(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hnum = GetITSsaSpectrum(file, AliPID::kKaon, 1, cent, kFALSE, kFALSE);
    hnum->Add(GetITSsaSpectrum(file, AliPID::kKaon, 0, cent, kFALSE, kFALSE));
    break;
  case kPrMinus:
    ptMin = TMath::Min(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kProton, 1, cent, kFALSE, kFALSE);
    break;
  case kPrPlus:
    ptMin = TMath::Min(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kProton, 0, cent, kFALSE, kFALSE);
    break;
  case kPr:
    ptMin = TMath::Min(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hnum = GetITSsaSpectrum(file, AliPID::kProton, 1, cent, kFALSE, kFALSE);
    hnum->Add(GetITSsaSpectrum(file, AliPID::kProton, 0, cent, kFALSE, kFALSE));
    break;
  }
  switch (den) {
  case kPiMinus:
    ptMin = TMath::Max(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kPion, 1, cent, kFALSE, kFALSE);
    break;
  case kPiPlus:
    ptMin = TMath::Max(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kPion, 0, cent, kFALSE, kFALSE);
    break;
  case kPi:
    ptMin = TMath::Max(ptMin, 0.1);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kPion, 1, cent, kFALSE, kFALSE);
    hden->Add(GetITSsaSpectrum(file, AliPID::kPion, 0, cent, kFALSE, kFALSE));
    break;
  case kKaMinus:
    ptMin = TMath::Max(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hden = GetITSsaSpectrum(file, AliPID::kKaon, 1, cent, kFALSE, kFALSE);
    break;
  case kKaPlus:
    ptMin = TMath::Max(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hden = GetITSsaSpectrum(file, AliPID::kKaon, 0, cent, kFALSE, kFALSE);
    break;
  case kKa:
    ptMin = TMath::Max(ptMin, 0.2);
    ptMax = TMath::Min(ptMax, 0.5);
    hden = GetITSsaSpectrum(file, AliPID::kKaon, 1, cent, kFALSE, kFALSE);
    hden->Add(GetITSsaSpectrum(file, AliPID::kKaon, 0, cent, kFALSE, kFALSE));
    break;
  case kPrMinus:
    ptMin = TMath::Max(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kProton, 1, cent, kFALSE, kFALSE);
    break;
  case kPrPlus:
    ptMin = TMath::Max(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kProton, 0, cent, kFALSE, kFALSE);
    break;
  case kPr:
    ptMin = TMath::Max(ptMin, 0.3);
    ptMax = TMath::Min(ptMax, 0.6);
    hden = GetITSsaSpectrum(file, AliPID::kProton, 1, cent, kFALSE, kFALSE);
    hden->Add(GetITSsaSpectrum(file, AliPID::kProton, 0, cent, kFALSE, kFALSE));
    break;
  }

  if (!hnum || !hden) return NULL;  

  Char_t sysname[1024];
  if (num == kPiMinus && den == kPiPlus)
    sprintf(sysname, "Pi_Pos2Neg");
  else if (num == kKaMinus && den == kKaPlus)
    sprintf(sysname, "K_Pos2Neg");
  else if (num == kPrMinus && den == kPrPlus)
    sprintf(sysname, "P_Pos2Neg");
  else if ((num == kKa || num == kKaPlus || num == kKaMinus)
	   && (den == kPi || den == kPiPlus || den == kPiMinus))
    sprintf(sysname, "K2Pi");
  else if ((num == kPr || num == kPrPlus || num == kPrMinus)
	   && (den == kPi || den == kPiPlus || den == kPiMinus))
    sprintf(sysname, "P2Pi");

  TH1D *hin = (TH1D *)hnum->Clone("hin");
  hin->Divide(hden);

  /* get systematics */
  TFile *fsys = TFile::Open("RATIOSYS_ITSsa.root");
  TH1 *hsys = fsys->Get(Form("hSystTot%s", sysname));
			
  TH1D *h = new TH1D(Form("hITSsa_cent%d_%s_%s", cent, ratioName[num], ratioName[den]), "ITSsa", NptBins, ptBin);
  Double_t pt, width, value, error, sys;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001) {
      //      printf("skipping %f because of sanity checks\n", pt);
      continue;
    }
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin || pt > ptMax)) {
      //      printf("skipping %f because of limits\n", pt);
      continue;
    }
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    /*** TEMP ADD SYS ***/
    if (addSystematicError) {
      sys = hsys->GetBinContent(bin) * value;
      error = TMath::Sqrt(error * error + sys * sys);
    }
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }

  h->SetTitle("ITSsa");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(1);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  return h;
}

//______________________________________________________________

Char_t *ITSTPCPartName[5] = {"", "", "Pion", "Kaon", "Proton"};
Char_t *ITSTPCChargeName[2] = {"Pos", "Neg"};
TH1D *
GetITSTPCSpectrum(TFile *file, Int_t part, Int_t charge, Int_t cent)
{
  TList *list = (TList *)file->Get("output");
  TH1D *hin = (TH1D *)list->FindObject(Form("h_%s_%s_cen_%d", ITSTPCPartName[part], ITSTPCChargeName[charge], cent + 1));
  if (!hin) return NULL;

  TH1D *h = new TH1D(Form("hITSTPC_cent%d_%s_%s", cent, AliPID::ParticleName(part), chargeName[charge]), "ITSTPC", NptBins, ptBin);
  Double_t pt, width, value, error;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }
  
#if 0
  /* add systematic error */
  Double_t sys;
  if (part == 2) sys = 0.5;
  else sys = 0.1;
  Double_t cont, conte;
  for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
    cont = h->GetBinContent(ipt + 1);
    conte = h->GetBinError(ipt + 1);
    conte = TMath::Sqrt(conte * conte + sys * sys * cont * cont);
    h->SetBinError(ipt + 1, conte);
  }
#endif
  
  h->SetTitle("ITSTPC");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(21);
  h->SetMarkerColor(2);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  return h;
}

//______________________________________________________________
//______________________________________________________________

Char_t *TPCTOFPartName[5] = {"", "", "pion", "kaon", "proton"};
Char_t *TPCTOFChargeName[2] = {"Pos", "Neg"};
TH1D *
GetTPCTOFSpectrum(TFile *file, Int_t part, Int_t charge, Int_t cent, Bool_t cutSpectrum = kTRUE)
{
  /* pt limits for combined spectra */
  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0., 0., 0.};
  Double_t ptMax[AliPID::kSPECIES] = {0., 0., 1.2, 1.2, 1.8};

  TH1D *hin = (TH1D *)file->Get(Form("%sFinal%s%d", TPCTOFPartName[part], TPCTOFChargeName[charge], cent));
  if (!hin) return NULL;

  TH1D *h = new TH1D(Form("hTPCTOF_cent%d_%s_%s", cent, AliPID::ParticleName(part), chargeName[charge]), "TPCTOF", NptBins, ptBin);
  Double_t pt, width, value, error;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin[part] || pt > ptMax[part])) continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }
  
  h->SetTitle("TPCTOF");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(22);
  h->SetMarkerColor(8);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  
  return h;
}

//______________________________________________________________

TH1D *
GetTPCTOFRatio(TFile *file, Int_t num, Int_t den, Int_t cent, Bool_t cutSpectrum = kTRUE)
{
  /* pt limits for combined spectra */
  Double_t ptMin_[9] = {
    0.0, 0.0, 0.0, 
    0., 0., 0., 
    0.5, 0.5, 0.5
  };
  Double_t ptMax_[9] = {
    1.2, 1.2, 1.2, 
    1.2, 1.2, 1.2, 
    1.8, 1.8, 1.8
  };

  Double_t ptMin = TMath::Max(ptMin_[num], ptMin_[den]);
  Double_t ptMax = TMath::Min(ptMax_[num], ptMax_[den]);

  Int_t part = 0, charge = 0;
  if (num == kPiMinus && den == kPiPlus) {
    part = AliPID::kPion;
    charge = 1;
  }
  else if (num == kKaMinus && den == kKaPlus) {
    part = AliPID::kKaon;
    charge = 1;
  }
  else if (num == kPrMinus && den == kPrPlus) {
    part = AliPID::kProton;
    charge = 1;
  }
  else if (num == kKaMinus && den == kPiMinus) {
    part = AliPID::kKaon;
    charge = 1;
  }
  else if (num == kKaPlus && den == kPiPlus) {
    part = AliPID::kKaon;
    charge = 0;
  }
  else if (num == kPrMinus && den == kPiMinus) {
    part = AliPID::kProton;
    charge = 1;
  }
  else if (num == kPrPlus && den == kPiPlus) {
    part = AliPID::kProton;
    charge = 0;
  }

  TH1D *hin = (TH1D *)file->Get(Form("%sFinal%s%d", TPCTOFPartName[part], TPCTOFChargeName[charge], cent));
  if (!hin) return NULL;

  TH1D *h = new TH1D(Form("hTPCTOF_cent%d_%s_%s", cent, ratioName[num], ratioName[den]), "TPCTOF", NptBins, ptBin);
  Double_t pt, width, value, error;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin || pt > ptMax)) continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }
  
  h->SetTitle("TPCTOF");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(22);
  h->SetMarkerColor(8);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  
  return h;
}

//______________________________________________________________
//______________________________________________________________

Char_t *TOFPartName[5] = {"", "", "pion", "kaon", "proton"};
Char_t *TOFChargeName[2] = {"positive", "negative"};
TH1D *
GetTOFSpectrum(TFile *file, Int_t part, Int_t charge, Int_t cent, Bool_t cutSpectrum = kTRUE)
{
  /* pt limits for combined spectra */
  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0.5, 0.45, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {0., 0., 3.0, 3.0, 4.5};

  TH1D *hin = (TH1D *)file->Get(Form("hFinal_cent%d_%s_%s", cent, TOFPartName[part], TOFChargeName[charge]));
  if (!hin) return NULL;

  /* get matching systematics */
  TFile *fsys = TFile::Open(Form("MATCHSYS_TOF_%s.root", TOFChargeName[charge]));
  TH1 *hsys = fsys->Get(Form("hErr%sMatch", ITSsaPartName[part]));
  TF1 *ffsys = new TF1("fsys", "[0] + [1] * x + [2] * TMath::Exp(-[3] * x)");
  ffsys->SetParameter(0, 0.02);
  ffsys->FixParameter(1, 0.);
  ffsys->SetParameter(2, 0.5);
  ffsys->SetParameter(3, 10.);
  hsys->Fit(ffsys, "W");
  ffsys->ReleaseParameter(1);
  hsys->Fit(ffsys, "W");
  hsys->Fit(ffsys, "W");
  hsys->Fit(ffsys, "W");
  hsys->Draw();
			
  TH1D *h = new TH1D(Form("hTOF_cent%d_%s_%s", cent, AliPID::ParticleName(part), chargeName[charge]), "TOF", NptBins, ptBin);
  Double_t pt, width, value, error, sys;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin[part] || pt > ptMax[part])) continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    /*** TEMP ADD SYS ***/
    sys = ffsys->Eval(pt) * value;
    error = TMath::Sqrt(error * error + sys * sys);
    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);

    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }
  
  h->SetTitle("TOF");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(23);
  h->SetMarkerColor(4);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  return h;
}

//______________________________________________________________

TH1D *
GetTOFRatio(TFile *file, Int_t num, Int_t den, Int_t cent, Bool_t cutSpectrum = kTRUE)
{
  /* pt limits for combined spectra */
  Double_t ptMin_[9] = {
    0.5, 0.5, 0.5, 
    0.45, 0.45, 0.45, 
    0.5, 0.5, 0.5
  };
  Double_t ptMax_[9] = {
    3.0, 3.0, 3.0, 
    3.0, 3.0, 3.0, 
    4.5, 4.5, 4.5
  };

  Double_t ptMin = TMath::Max(ptMin_[num], ptMin_[den]);
  Double_t ptMax = TMath::Min(ptMax_[num], ptMax_[den]);
  
  TH1D *hin = (TH1D *)file->Get(Form("hRatio_cent%d_%s_%s", cent, ratioName[num], ratioName[den]));
  if (!hin) return NULL;


#if 0
  /* get matching systematics */
  TFile *fsys = TFile::Open(Form("MATCHSYS_TOF_%s.root", TOFChargeName[charge]));
  TH1 *hsys = fsys->Get(Form("hErr%sMatch", ITSsaPartName[part]));
  TF1 *ffsys = new TF1("fsys", "[0] + [1] * x + [2] * TMath::Exp(-[3] * x)");
  ffsys->SetParameter(0, 0.02);
  ffsys->FixParameter(1, 0.);
  ffsys->SetParameter(2, 0.5);
  ffsys->SetParameter(3, 10.);
  hsys->Fit(ffsys, "W");
  ffsys->ReleaseParameter(1);
  hsys->Fit(ffsys, "W");
  hsys->Fit(ffsys, "W");
  hsys->Fit(ffsys, "W");
  hsys->Draw();
#endif
			
  TH1D *h = new TH1D(Form("hTOF_cent%d_%s_%s", cent, ratioName[num], ratioName[den]), "TOF", NptBins, ptBin);
  Double_t pt, width, value, error, sys;
  Int_t bin;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    /* get input bin */
    pt = h->GetBinCenter(ipt + 1);
    width = h->GetBinWidth(ipt + 1);
    bin = hin->FindBin(pt);
    /* sanity check */
    if (TMath::Abs(hin->GetBinCenter(bin) - pt) > 0.001 ||
	TMath::Abs(hin->GetBinWidth(bin) - width) > 0.001)
      continue;
    /* check pt limits */
    if (cutSpectrum && (pt < ptMin || pt > ptMax)) continue;
    /* copy bin */
    value = hin->GetBinContent(bin);
    error = hin->GetBinError(bin);
    /*** TEMP ADD SYS ***/
    //    sys = ffsys->Eval(pt) * value;
    //    error = TMath::Sqrt(error * error + sys * sys);
    //    h->SetBinContent(ipt + 1, value);
    //    h->SetBinError(ipt + 1, error);

    h->SetBinContent(ipt + 1, value);
    h->SetBinError(ipt + 1, error);
  }
  
  h->SetTitle("TOF");
  h->SetLineWidth(1);
  h->SetLineColor(1);
  h->SetMarkerStyle(23);
  h->SetMarkerColor(4);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  return h;
}

//______________________________________________________________

Char_t *HydroPartName[5] = {"", "", "Pion", "Kaon", "Proton"};
TGraph *
GetHydroSpectrum(TFile *file, Int_t part, Int_t charge, Int_t cent)
{
  TGraph *h = (TGraph *)file->Get(Form("%s_C%d", HydroPartName[part], cent));
  if (!h) return NULL;
  h->SetTitle("Hydro");
  h->SetLineWidth(2);
  h->SetLineColor(kYellow+1);
  h->SetMarkerStyle(24);
  h->SetMarkerColor(kYellow+1);
  h->SetFillStyle(0);
  h->SetFillColor(0);
  return h;
}

//______________________________________________________________

AddPointsToGraph(TH1D *h, TGraphErrors *g)
{
  Int_t ipoint;
  for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
    if (h->GetBinError(ipt + 1) <= 0. || h->GetBinContent(ipt + 1) <= 0.) continue;
    ipoint = g->GetN();
    g->SetPoint(ipoint, h->GetBinCenter(ipt + 1), h->GetBinContent(ipt + 1));
    g->SetPointError(ipoint, 0.5 * h->GetBinWidth(ipt + 1), h->GetBinError(ipt + 1));
  }
}

//______________________________________________________________

AddPointsToProfile(TH1D *h, TProfile *p)
{
  for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
    if (h->GetBinError(ipt + 1) <= 0. || h->GetBinContent(ipt + 1) <= 0.) continue;
    p->Fill(h->GetBinCenter(ipt + 1), h->GetBinContent(ipt + 1));
  }
}

//______________________________________________________________

CombineSpectra(TH1D *h, TObjArray *oa)
{

  /* combine with weighted mean */
  Double_t ptcen, w, sumw, mean, meane, cont, conte, syse, tote, minvalue, maxvalue, maxerror, maxmindiff;
  Int_t ptbin;
  TH1D *hin;
  for (Int_t ipt = 0; ipt < h->GetNbinsX(); ipt++) {
    mean = 0.;
    sumw = 0.;
    minvalue = kMaxInt;
    maxvalue = 0.;
    maxerror = 0.;
    ptcen = h->GetBinCenter(ipt + 1);
    TProfile prof(Form("prof_%d", ipt), "", 1, 0., 1.); /* to get RMS */
    for (Int_t ih = 0; ih < oa->GetEntries(); ih++) {
      hin = (TH1D *)oa->At(ih);
      ptbin = hin->FindBin(ptcen);
      if (ptbin <= 0.) continue;
      cont = hin->GetBinContent(ptbin);
      conte = hin->GetBinError(ptbin);
      if (cont == 0. || conte == 0.)
	continue;
      w = 1. / (conte * conte);
      mean += cont * w;
      sumw += w;
      if (cont < minvalue)
	minvalue = cont;
      if (cont > maxvalue)
	maxvalue = cont;
      if (conte > maxerror)
	maxerror = conte;
      prof.Fill(0., cont);
    }
    if (sumw == 0.) {
      mean = 0.;
      meane = 0.;	
      syse = 0.;
      tote = 0.;
    }
    else {
      mean /= sumw;
      meane = TMath::Sqrt(1. / sumw);
      syse = 0.5 * (maxvalue -  minvalue);
      tote = TMath::Sqrt(meane * meane + syse * syse);
      //      if (tote < maxerror)
      //	tote = maxerror;
    }

    //    printf("pt = %f, mean = %f, meane = %f, syse = %f, tote = %f\n", ptcen, mean, meane, syse, tote);
    
    //    syse = prof.GetBinError(1);
    h->SetBinContent(ipt + 1, mean);
    h->SetBinError(ipt + 1, tote);
    //    h->SetBinError(ipt + 1, TMath::Sqrt(meane * meane + syse * syse + mean * mean * 0.1 * 0.1)); /* add extra 10% systematics */
  }

}

//______________________________________________________________

AntiparticleParticleRatio(Int_t part, Int_t cent, Bool_t cutSpectra = kTRUE)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  
  TFile *itssafile = TFile::Open(itssafilename);
  //  TFile *itstpcfile = TFile::Open(itstpcfilename);
  TFile *tpctoffile = TFile::Open(tpctoffilename);
  TFile *toffile = TFile::Open(toffilename);
  //  TFile *hydrofile = TFile::Open(hydrofilename);

  TCanvas *cCanvas = new TCanvas("cCanvas");
  if (cent == -1) cCanvas->Divide(5, 2);
  TCanvas *cCanvasCombined = new TCanvas("cCanvasCombined");
  if (cent == -1) cCanvasCombined->Divide(5, 2);
  TCanvas *cCanvasRatio = new TCanvas("cCanvasRatio");
  if (cent == -1) cCanvasRatio->Divide(5, 2);
  TCanvas *cCanvasRatioFit = new TCanvas("cCanvasRatioFit");
  if (cent == -1) cCanvasRatioFit->Divide(5, 2);
  TH1D *hITSsa[2], *hITSTPC[2], *hTPCTOF[2], *hTOF[2];
  TGraphErrors *gCombined[2][10];
  TProfile *pCombined[2][10];
  TH1D *hCombined[2][10];
  TH1D *hRatio_ITSsa[10];
  TH1D *hRatio_ITSTPC[10];
  TH1D *hRatio_TPCTOF[10];
  TH1D *hRatio_TOF[10];
  TH1D *hRatio_COMB[10];
  for (Int_t icent = 0; icent < 10; icent++) {
    
    if (cent != -1 && icent != cent) continue;
    
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      gCombined[icharge][icent] = new TGraphErrors();
      pCombined[icharge][icent] = new TProfile(Form("pCombined_charge%d_cent%d", icharge, icent), "", NptBins, ptBin);
      hCombined[icharge][icent] = new TH1D(Form("hCombined_charge%d_cent%d", icharge, icent), "", NptBins, ptBin);
      TObjArray spectraArray;
      hITSsa[icharge] = GetITSsaSpectrum(itssafile, part, icharge, icent, cutSpectra);
      hITSTPC[icharge] = NULL;//GetITSTPCSpectrum(itstpcfile, part, icharge, icent, cutSpectra);
      hTPCTOF[icharge] = GetTPCTOFSpectrum(tpctoffile, part, icharge, icent, cutSpectra);
      hTOF[icharge] = GetTOFSpectrum(toffile, part, icharge, icent, cutSpectra);
      //    hHydro = GetHydroSpectrum(hydrofile, part, charge, icent);
      Double_t minimum = 0.;
      Double_t maximum = 0.;
      if (hITSsa[icharge]) {
	//	AddPointsToGraph(hITSsa[icharge], gCombined[icharge][icent]);
	//	AddPointsToProfile(hITSsa[icharge], pCombined[icharge][icent]);
	spectraArray.Add(hITSsa[icharge]);
	//	hITSsa[icharge]->DrawCopy("same");
	//	if (hITSsa[icharge]->GetMaximum() > maximum) maximum = hITSsa[icharge]->GetMaximum();
	//	if (hITSsa[icharge]->GetMinimum() < minimum) minimum = hITSsa[icharge]->GetMinimum();
      }
      if (hITSTPC[icharge]) {
	//	AddPointsToGraph(hITSTPC[icharge], gCombined[icharge][icent]);
	//	AddPointsToProfile(hITSTPC[icharge], pCombined[icharge][icent]);
	spectraArray.Add(hITSTPC[icharge]);
	//	hITSTPC[icharge]->DrawCopy("same");
	//	if (hITSTPC[icharge]->GetMaximum() > maximum) maximum = hITSTPC[icharge]->GetMaximum();
	//	if (hITSTPC[icharge]->GetMinimum() < minimum) minimum = hITSTPC[icharge]->GetMinimum();
      }
      if (hTPCTOF[icharge]) {
	//	AddPointsToGraph(hTPCTOF[icharge], gCombined[icharge][icent]);
	//	AddPointsToProfile(hTPCTOF[icharge], pCombined[icharge][icent]);
	spectraArray.Add(hTPCTOF[icharge]);
	//	hTPCTOF[icharge]->DrawCopy("same");
	//	if (hTPCTOF[icharge]->GetMaximum() > maximum) maximum = hTPCTOF[icharge]->GetMaximum();
	//	if (hTPCTOF[icharge]->GetMinimum() < minimum) minimum = hTPCTOF[icharge]->GetMinimum();
      }
      if (hTOF[icharge]) {
	//	AddPointsToGraph(hTOF[icharge], gCombined[icharge][icent]);
	//	AddPointsToProfile(hTOF[icharge], pCombined[icharge][icent]);
	spectraArray.Add(hTOF[icharge]);
	//	hTOF[icharge]->DrawCopy("same");
	//	if (hTOF[icharge]->GetMaximum() > maximum) maximum = hTOF[icharge]->GetMaximum();
	//	if (hTOF[icharge]->GetMinimum() < minimum) minimum = hTOF[icharge]->GetMinimum();
      }
      
      CombineSpectra(hCombined[icharge][icent], &spectraArray);

    }

    /* antiparticle/particle ratios */

    Char_t title[1024];
    sprintf(title, "%s/%s (%s);p_{T} (GeV/c);%s/%s;", partChargeName[part][1], partChargeName[part][0], centName[icent], partChargeName[part][1], partChargeName[part][0]);

    TH1D *hArea = new TH1D(Form("hArea_%d", icent), title, 100, ptMin[part], ptMax[part]);
    hArea->SetMaximum(1.2);
    hArea->SetMinimum(0.8);

    if (cent == -1)
      cCanvas->cd(icent + 1);
    else
      cCanvas->cd();

    hArea->DrawCopy();

    if (hITSsa[0] && hITSsa[1]) {
      hRatio_ITSsa[icent] = new TH1D(*hITSsa[1]);
      hRatio_ITSsa[icent]->Divide(hITSsa[0]);
      hRatio_ITSsa[icent]->DrawCopy("same");
    }

    if (hTPCTOF[0] && hTPCTOF[1]) {
      hRatio_TPCTOF[icent] = new TH1D(*hTPCTOF[1]);
      hRatio_TPCTOF[icent]->Divide(hTPCTOF[0]);
      hRatio_TPCTOF[icent]->DrawCopy("same");
    }

    if (hTOF[0] && hTOF[1]) {
      hRatio_TOF[icent] = new TH1D(*hTOF[1]);
      hRatio_TOF[icent]->Divide(hTOF[0]);
      hRatio_TOF[icent]->DrawCopy("same");
    }

    //    TLegend *legend = cCanvas->BuildLegend();
    //    legend->SetFillStyle(0);
    //    legend->SetFillColor(0);
    //    legend->DeleteEntry();

    if (cent == -1)
      cCanvasCombined->cd(icent + 1);
    else
      cCanvasCombined->cd();

    if (hCombined[0][icent] && hCombined[1][icent]) {
      hRatio_COMB[icent] = new TH1D(*hCombined[1][icent]);
      hRatio_COMB[icent]->Divide(hCombined[0][icent]);
      hRatio_COMB[icent]->SetMarkerStyle(20);
      hRatio_COMB[icent]->SetMarkerColor(2);
      hRatio_COMB[icent]->Fit("pol0", "");
      hRatio_COMB[icent]->SetMinimum(0.5);
      hRatio_COMB[icent]->SetMaximum(1.5);
      hRatio_COMB[icent]->SetTitle(title);
      //      hRatio_COMB[icent]->Draw("PHIST");
      //      pol0->SetRange(0., 5.);
      //      pol0->DrawCopy("same");
      //      hArea->DrawCopy();
      //      hRatio_COMB[icent]->Draw("same");
    
    }


  }

}

//______________________________________________________________

CompareSpectra(Int_t part, Int_t charge, Int_t cent = -1, Int_t ratio = kFALSE, Int_t fitfunc = -1, Bool_t cutSpectra = kTRUE) 
{

  gROOT->LoadMacro("HistoUtils.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  LoadLibraries();
  AliBWFunc bwf;
  bwf.SetVarType(AliBWFunc::kdNdpt);
  TF1 *fFitFunc = NULL;

  switch (fitfunc) {
  case 0:
    gROOT->LoadMacro("SpectraAnalysis.C");
    fFitFunc = STAR_BlastWave("fBW", AliPID::ParticleMass(part), 0.9, 0.1, 1.);
    fBW->SetParLimits(3, 0.5, 1.5);
    //    fBW->FixParameter(3, 1.);
    break;
  case 1:
    fFitFunc = bwf.GetLevi(AliPID::ParticleMass(part), AliPID::ParticleMass(part), 5., 1000.);
    break;
  case 2:
    fFitFunc = bwf.GetBoltzmann(AliPID::ParticleMass(part), AliPID::ParticleMass(part), 100.);
    break;
  case 3:
    fFitFunc = bwf.GetMTExp(AliPID::ParticleMass(part),AliPID::ParticleMass(part) , 100.);
    break;
  case 4:
    fFitFunc = bwf.GetPTExp(AliPID::ParticleMass(part), 100.);
    break;
  case 5:
    fFitFunc = bwf.GetBGBW(AliPID::ParticleMass(part), 0.5, 0.1, 1.e6);
    break;
  case 6:
    fFitFunc = new TF1("fpol9", "pol9", 0., 5.0);
    break;
  }
  if (fFitFunc) fFitFunc->SetLineWidth(2);

  TFile *itssafile = TFile::Open(ratio ? itssaratiofilename : itssafilename);
  //  TFile *itstpcfile = TFile::Open(itstpcfilename);
  if (part / 3 == charge / 3)
    Char_t *tpctofratiofilename = tpctofratiofilenameA;
  else
    Char_t *tpctofratiofilename = tpctofratiofilenameB;
  TFile *tpctoffile = TFile::Open(ratio ? tpctofratiofilename : tpctoffilename);
  TFile *toffile = TFile::Open(ratio ? tofratiofilename : toffilename);
  //  TFile *hydrofile = TFile::Open(hydrofilename);

  TCanvas *cCanvas = new TCanvas("cCanvas");
  if (cent == -1) cCanvas->Divide(5, 2);
  TCanvas *cCanvasCombined = new TCanvas("cCanvasCombined");
  TCanvas *cCanvasRatio = new TCanvas("cCanvasRatio");
  if (cent == -1) cCanvasRatio->Divide(5, 2);
  TCanvas *cCanvasRatioComb = new TCanvas("cCanvasRatioComb");
  if (cent == -1) cCanvasRatioComb->Divide(5, 2);
  TCanvas *cCanvasRatioFit = new TCanvas("cCanvasRatioFit");
  if (cent == -1) cCanvasRatioFit->Divide(5, 2);
  TPad *curpad = NULL;
  TH1D *hITSsa, *hITSTPC, *hTPCTOF, *hTOF;
  TGraph *hHydro;
  TGraphErrors *gCombined[10];
  TProfile *pCombined[10];
  TH1D *hCombined[10];
  TH1D *hRatio_ITSsa_ITSTPC[10];
  TH1D *hRatio_ITSsa_TPCTOF[10];
  TH1D *hRatio_ITSTPC_TPCTOF[10];
  TH1D *hRatio_ITSTPC_TOF[10];
  TH1D *hRatio_TPCTOF_TOF[10];
  TH1D *hRatio_ITSsa_TOF[10];
  for (Int_t icent = 0; icent < 10; icent++) {
    if (cent != -1 && icent != cent) continue;
    gCombined[icent] = new TGraphErrors();
    pCombined[icent] = new TProfile(Form("pCombined_cent%d", icent), "", NptBins, ptBin);
    hCombined[icent] = new TH1D(Form("hCombined_cent%d", icent), "", NptBins, ptBin);
    TObjArray spectraArray;
    hITSsa = ratio ? GetITSsaRatio(itssafile, part, charge, icent, cutSpectra): GetITSsaSpectrum(itssafile, part, charge, icent, cutSpectra);
    //    hITSTPC = GetITSTPCSpectrum(itstpcfile, part, charge, icent, cutSpectra);
    hTPCTOF = ratio ? GetTPCTOFRatio(tpctoffile, part, charge, icent, cutSpectra) : GetTPCTOFSpectrum(tpctoffile, part, charge, icent, cutSpectra);
    hTOF = ratio ? GetTOFRatio(toffile, part, charge, icent, cutSpectra) : GetTOFSpectrum(toffile, part, charge, icent, cutSpectra);
    //    hHydro = GetHydroSpectrum(hydrofile, part, charge, icent);
    if (cent == -1)
      curpad = (TPad *)cCanvas->cd(icent + 1);
    else
      curpad = (TPad *)cCanvas->cd();
    if (!ratio)
      TH1D *hArea = new TH1D(Form("hArea_%d", icent), Form("%s (%s);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{t}};", partChargeName[part][charge], centName[icent]), 100, 0., 5.);
    else
      TH1D *hArea = new TH1D(Form("hArea_%d", icent), Form("%s (%s);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{t}};", "generic ratio", centName[icent]), 100, 0., 5.);
    hArea->Draw();
    Double_t minimum = 0.001;
    Double_t maximum = 1000.;
    if (hITSsa) {
      AddPointsToGraph(hITSsa, gCombined[icent]);
      AddPointsToProfile(hITSsa, pCombined[icent]);
      spectraArray.Add(hITSsa);
      hITSsa->DrawCopy("same");
      if (hITSsa->GetMaximum() > maximum) maximum = hITSsa->GetMaximum();
      if (hITSsa->GetMinimum() < minimum) minimum = hITSsa->GetMinimum();
    }
    if (hITSTPC) {
      AddPointsToGraph(hITSTPC, gCombined[icent]);
      AddPointsToProfile(hITSTPC, pCombined[icent]);
      spectraArray.Add(hITSTPC);
      hITSTPC->DrawCopy("same");
      if (hITSTPC->GetMaximum() > maximum) maximum = hITSTPC->GetMaximum();
      if (hITSTPC->GetMinimum() < minimum) minimum = hITSTPC->GetMinimum();
    }
    if (hTPCTOF) {
      AddPointsToGraph(hTPCTOF, gCombined[icent]);
      AddPointsToProfile(hTPCTOF, pCombined[icent]);
      spectraArray.Add(hTPCTOF);
      hTPCTOF->DrawCopy("same");
      if (hTPCTOF->GetMaximum() > maximum) maximum = hTPCTOF->GetMaximum();
      if (hTPCTOF->GetMinimum() < minimum) minimum = hTPCTOF->GetMinimum();
    }
    if (hTOF) {
      AddPointsToGraph(hTOF, gCombined[icent]);
      AddPointsToProfile(hTOF, pCombined[icent]);
      spectraArray.Add(hTOF);
      hTOF->DrawCopy("same");
      if (hTOF->GetMaximum() > maximum) maximum = hTOF->GetMaximum();
      if (hTOF->GetMinimum() < minimum) minimum = hTOF->GetMinimum();
    }
    if (hHydro) {
      ;//hHydro->Draw("c,same");
    }
    TLegend *legend = curpad->BuildLegend();
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->DeleteEntry();

    hArea->SetMaximum(maximum * 1.1);
    hArea->SetMinimum(0.01);
    //    gPad->SetLogy();

    /*** RATIOS ***/

    /* switch canvas */
    if (cent == -1)
      curpad = (TPad *)cCanvasRatio->cd(icent + 1);
    else
      curpad = (TPad *)cCanvasRatio->cd();

   
    /* area histo */
    if (!ratio)
      TH1D *hAreaRatio = new TH1D(Form("hAreaRatio_%d", icent), Form("%s (%s);p_{T} (GeV/c);ratio;", partChargeName[part][charge], centName[icent]), 100, 0., 5.);
    else
      TH1D *hAreaRatio = new TH1D(Form("hAreaRatio_%d", icent), Form("%s (%s);p_{T} (GeV/c);ratio;", "generic ratio", centName[icent]), 100, 0., 5.);

    hAreaRatio->SetMaximum(1.5);
    hAreaRatio->SetMinimum(0.5);
    hAreaRatio->Draw();

    /* do ratios */
    if (hITSsa && hITSTPC) {
      hRatio_ITSsa_ITSTPC[icent] = new TH1D(*hITSsa);
      hRatio_ITSsa_ITSTPC[icent]->Divide(hITSTPC);
      hRatio_ITSsa_ITSTPC[icent]->SetNameTitle(Form("hRatio_ITSsa_ITSTPC_cent%d", icent), "ITSsa / ITSTPC");
      hRatio_ITSsa_ITSTPC[icent]->Draw("same");
    }
    if (hITSsa && hTPCTOF) {
      hRatio_ITSsa_TPCTOF[icent] = new TH1D(*hITSsa);
      hRatio_ITSsa_TPCTOF[icent]->Divide(hTPCTOF);
      hRatio_ITSsa_TPCTOF[icent]->SetNameTitle(Form("hRatio_ITSsa_TPCTOF_cent%d", icent), "ITSsa / TPCTOF");
      hRatio_ITSsa_TPCTOF[icent]->SetMarkerStyle(23);
      hRatio_ITSsa_TPCTOF[icent]->SetMarkerColor(4);
      hRatio_ITSsa_TPCTOF[icent]->Draw("same");
    }
    if (hITSTPC && hTPCTOF) {
      hRatio_ITSTPC_TPCTOF[icent] = new TH1D(*hITSTPC);
      hRatio_ITSTPC_TPCTOF[icent]->Divide(hTPCTOF);
      hRatio_ITSTPC_TPCTOF[icent]->SetNameTitle(Form("hRatio_ITSTPC_TPCTOF_cent%d", icent), "ITSTPC / TPCTOF");
      hRatio_ITSTPC_TPCTOF[icent]->Draw("same");
    }
    if (hTPCTOF && hTOF) {
      hRatio_TPCTOF_TOF[icent] = new TH1D(*hTPCTOF);
      hRatio_TPCTOF_TOF[icent]->Divide(hTOF);
      hRatio_TPCTOF_TOF[icent]->SetNameTitle(Form("hRatio_TPCTOF_TOF_cent%d", icent), "TPCTOF / TOF");
      hRatio_TPCTOF_TOF[icent]->Draw("same");
    }
    if (hITSsa && hTOF) {
      hRatio_ITSsa_TOF[icent] = new TH1D(*hITSsa);
      hRatio_ITSsa_TOF[icent]->Divide(hTOF);
      hRatio_ITSsa_TOF[icent]->SetNameTitle(Form("hRatio_ITSsa_TOF_cent%d", icent), "ITSsa / TOF");
      //      hRatio_ITSsa_TOF[icent]->SetMarkerStyle(25);
      //      hRatio_ITSsa_TOF[icent]->SetMarkerColor(2);
      hRatio_ITSsa_TOF[icent]->Draw("same");
    }

    /* legend */
    TLegend *legendRatio = curpad->BuildLegend();
    legendRatio->SetFillStyle(0);
    legendRatio->SetFillColor(0);
    legendRatio->DeleteEntry();

    CombineSpectra(hCombined[icent], &spectraArray);
    hCombined[icent]->SetFillStyle(0);
    hCombined[icent]->SetFillColor(kOrange + 1);
    hCombined[icent]->SetMarkerColor(kOrange+1);
    hCombined[icent]->SetMarkerStyle(24);
    hCombined[icent]->SetLineColor(kOrange+1);
    hCombined[icent]->SetLineWidth(2);
    hCombined[icent]->SetMarkerSize(0);
    
    //    hCombined[icent]->DrawCopy("same,E2");
    //    pCombined[icent]->DrawCopy("same");

    if (cent == -1)
      cCanvas->cd(icent + 1);
    else
      cCanvas->cd();
    //    hCombined[icent]->Draw("same, E2");

    cCanvasCombined->cd();
    if (cent == -1 && icent != 0)
      hCombined[icent]->Draw("E2,same");
    else
      hCombined[icent]->Draw("E2");

    //    cCanvasCombined->DrawClonePad();
    if (hITSsa) {
      hITSsa->DrawCopy("same");
    }
    if (hITSTPC) {
      hITSTPC->DrawCopy("same");
    }
    if (hTPCTOF) {
      hTPCTOF->DrawCopy("same");
    }
    if (hTOF) {
      hTOF->DrawCopy("same");
    }
    if (hHydro) {
      ;//hHydro->Draw("c,same");
    }

    if (cent == -1)
      cCanvasRatioComb->cd(icent + 1);
    else
      cCanvasRatioComb->cd();
    //    hCombined[icent]->Draw("same, E2");

    TH1 *hhr = HistoUtils_smartratio(hCombined[icent], hCombined[icent]);
    hhr->SetMaximum(1.25);
    hhr->SetMinimum(0.75);
    hhr->SetFillStyle(3001);
    hhr->SetTitle("combined error;p_{T} (GeV/c);ratio wrt. combined");
    hhr->Draw("e2");
    if (hITSsa) {
      hhr = HistoUtils_smartratio(hITSsa, hCombined[icent]);
      hhr->SetLineColor(1);
      hhr->SetLineWidth(2);
      hhr->Draw("e2,same");
    }
    if (hITSTPC) {
      hhr = HistoUtils_smartratio(hITSTPC, hCombined[icent]);
      hhr->SetLineColor(1);
      hhr->SetLineWidth(2);
      hhr->Draw("e2,same");
    }
    if (hTPCTOF) {
      hhr = HistoUtils_smartratio(hTPCTOF, hCombined[icent]);
      hhr->SetLineColor(8);
      hhr->SetLineWidth(2);
      hhr->Draw("e2,same");
    }
    if (hTOF) {
      hhr = HistoUtils_smartratio(hTOF, hCombined[icent]);
      hhr->SetLineColor(4);
      hhr->SetLineWidth(2);
      hhr->Draw("e2,same");
    }
    if (hHydro) {
      ;//hHydro->Draw("c,same");
    }


    if (!fFitFunc)
      continue;
    
    //    gCombined[icent]->Draw("p*");
    //    gCombined[icent]->Fit(fFitFunc, "0q", "", 0.5, 1.0);
    //    gCombined[icent]->Fit(fFitFunc, "0q", "", 0.2, 1.5);
    hCombined[icent]->Fit(fFitFunc, "0q", "", 0., 5.);
    fFitFunc->DrawCopy("same");
    printf("cent = %d, dNdy = %f +- %f\n", icent, fFitFunc->GetParameter(0), fFitFunc->GetParError(0));

    if (cent == -1)
      cCanvas->cd(icent + 1);
    else
      cCanvas->cd();
    fFitFunc->DrawCopy("same");

    if (cent == -1)
      cCanvasRatioFit->cd(icent + 1);
    else
      cCanvasRatioFit->cd();
    if (!ratio)
      TH1D *hAreaRatioFit = new TH1D(Form("hAreaRatioFit_%d", icent), Form("%s (%s);p_{T} (GeV/c);ratio wrt. fit;", partChargeName[part][charge], centName[icent]), 100, 0., 5.);
    else
      TH1D *hAreaRatioFit = new TH1D(Form("hAreaRatioFit_%d", icent), Form("%s (%s);p_{T} (GeV/c);ratio wrt. fit;", "generic ratio", centName[icent]), 100, 0., 5.);
    hAreaRatioFit->SetMaximum(1.5);
    hAreaRatioFit->SetMinimum(0.5);
    hAreaRatioFit->Draw();
    legend->Draw("same");

    if (hITSsa) {
      hITSsa->Divide(fFitFunc);
      hITSsa->DrawCopy("same");
    }
    if (hITSTPC) {
      hITSTPC->Divide(fFitFunc);
      hITSTPC->DrawCopy("same");
    }
    if (hTPCTOF) {
      hTPCTOF->Divide(fFitFunc);
      hTPCTOF->DrawCopy("same");
    }
    if (hTOF) {
      hTOF->Divide(fFitFunc);
      hTOF->DrawCopy("same");
    }

  }

}

//______________________________________________________________

ConvertSpectraNameITSsa(const Char_t *fileoutname)
{

  TFile *itssafile = TFile::Open(itssafilename);
  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  TH1D *hSpectrum;
  for (Int_t part = 2; part < AliPID::kSPECIES; part++) {
    for (Int_t charge = 0; charge < 2; charge++) {
      for (Int_t icent = 0; icent < 9; icent++) {
	hSpectrum = GetITSsaSpectrum(itssafile, part, charge, icent);
	if (!hSpectrum)
	  continue;
	hSpectrum->SetName(Form("cent%d_%s_%s", icent, AliPID::ParticleName(part), chargeName[charge]));
	fileout->cd();
	hSpectrum->Write();
      }
    }
  }
  fileout->Close();
}

//______________________________________________________________

ConvertSpectraNameITSTPC(const Char_t *fileoutname)
{

  TFile *itstpcfile = TFile::Open(itstpcfilename);
  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  TH1D *hSpectrum;
  for (Int_t part = 2; part < AliPID::kSPECIES; part++) {
    for (Int_t charge = 0; charge < 2; charge++) {
      for (Int_t icent = 0; icent < 9; icent++) {
	hSpectrum = GetITSTPCSpectrum(itstpcfile, part, charge, icent);
	if (!hSpectrum)
	  continue;
	hSpectrum->SetName(Form("cent%d_%s_%s", icent, AliPID::ParticleName(part), chargeName[charge]));
	fileout->cd();
	hSpectrum->Write();
      }
    }
  }
  fileout->Close();
}

//______________________________________________________________

ConvertSpectraNameTPCTOF(const Char_t *fileoutname)
{

  TFile *tpctoffile = TFile::Open(tpctoffilename);
  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  TH1D *hSpectrum;
  for (Int_t part = 2; part < AliPID::kSPECIES; part++) {
    for (Int_t charge = 0; charge < 2; charge++) {
      for (Int_t icent = 0; icent < 9; icent++) {
	hSpectrum = GetTPCTOFSpectrum(tpctoffile, part, charge, icent);
	if (!hSpectrum)
	  continue;
	hSpectrum->SetName(Form("cent%d_%s_%s", icent, AliPID::ParticleName(part), chargeName[charge]));
	fileout->cd();
	hSpectrum->Write();
      }
    }
  }
  fileout->Close();
}

//______________________________________________________________

ConvertSpectraNameTOF(const Char_t *fileoutname)
{

  TFile *toffile = TFile::Open(toffilename);
  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  TH1D *hSpectrum;
  for (Int_t part = 2; part < AliPID::kSPECIES; part++) {
    for (Int_t charge = 0; charge < 2; charge++) {
      for (Int_t icent = 0; icent < 9; icent++) {
	hSpectrum = GetTOFSpectrum(toffile, part, charge, icent);
	if (!hSpectrum)
	  continue;
	hSpectrum->SetName(Form("cent%d_%s_%s", icent, AliPID::ParticleName(part), chargeName[charge]));
	fileout->cd();
	hSpectrum->Write();
      }
    }
  }
  fileout->Close();
}

//______________________________________________________________

SummedSpectra(const Char_t *filename, const Char_t *fileoutname)
{

  TFile *filein = TFile::Open(filename);
  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  TH1D *hSpectrum, *hSummed;
  for (Int_t icent = 0; icent < 9; icent++) {
    for (Int_t ipart = 2; ipart < 5; ipart++) {
      for (Int_t icharge = 0; icharge < 2; icharge++) {
	hSpectrum = (TH1D *)filein->Get(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	if (!hSpectrum) continue;
	if (ipart == 2 && icharge == 0) {
	  hSummed = new TH1D(*hSpectrum);
	  hSummed->SetName(Form("cent%d_summed", icent));
	}
	else {
	  hSummed->Add(hSpectrum);
	}
      }
    }
    fileout->cd();
    if (hSummed)
      hSummed->Write();
  }
  fileout->Close();
}

//______________________________________________________________

CombinedSpectra(const Char_t *fileoutname, Int_t what = -1)
{

  TFile *itssafile = TFile::Open(itssafilename);
  //  TFile *itstpcfile = TFile::Open(itstpcfilename);
  TFile *tpctoffile = TFile::Open(tpctoffilename);
  TFile *toffile = TFile::Open(toffilename);

  TFile *fileout = TFile::Open(fileoutname, "RECREATE");

  TH1D *hITSsa, *hITSTPC, *hTPCTOF, *hTOF;
  TH1D *hCombined[10];
  for (Int_t part = 2; part < AliPID::kSPECIES; part++) {
    for (Int_t charge = 0; charge < 2; charge++) {
      for (Int_t icent = 0; icent < 10; icent++) {
	hCombined[icent] = new TH1D(Form("cent%d_%s_%s", icent, AliPID::ParticleName(part), chargeName[charge]), "", NptBins, ptBin);
	TObjArray spectraArray;
	hITSsa = GetITSsaSpectrum(itssafile, part, charge, icent);
	//hITSTPC = GetITSTPCSpectrum(itstpcfile, part, charge, icent);
	hTPCTOF = GetTPCTOFSpectrum(tpctoffile, part, charge, icent);
	hTOF = GetTOFSpectrum(toffile, part, charge, icent);
	if (hITSsa && (what == -1 || what == 0)) {
	  spectraArray.Add(hITSsa);
	}
	if (hITSTPC && (what == -1 || what == 1)) {
	  spectraArray.Add(hITSTPC);
	}
	if (hTPCTOF && (what == -1 || what == 2)) {
	  spectraArray.Add(hTPCTOF);
	}
	if (hTOF && (what == -1 || what == 3)) {
	  spectraArray.Add(hTOF);
	}
	
	CombineSpectra(hCombined[icent], &spectraArray);
	fileout->cd();
	hCombined[icent]->Write();
      }
    }
  }	  
  
  fileout->Close();

}

