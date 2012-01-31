/**************************************************************/
enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
const Char_t *chargeName[kNCharges] = {
  "positive",
  "negative",
};
/**************************************************************/
const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
/**************************************************************/
const Int_t NptBins = 46;
Double_t ptBin[NptBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
/**************************************************************/
Char_t *partLatex[AliPID::kSPECIES][2] = {
  "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
};
 /**************************************************************/

FinalSpectra(const Char_t *rawfilename = "TOF_rawSpectra.root", const Char_t *matchfilename = "TOF_matchingEfficiency.root", const Char_t *trackfilename = "TOF_trackingEfficiency.root", const Char_t *primaryfilename = "TOF_primaryFraction.root", const Char_t *electronfilename = "TOF_electronCorrection.root")
 {

   Int_t marker[2] = {20, 25};
   Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};

   TFile *rawfile = TFile::Open(rawfilename);
   TFile *matchfile = TFile::Open(matchfilename);
   TFile *trackfile = TFile::Open(trackfilename);
   TFile *primaryfile = TFile::Open(primaryfilename);
   TFile *electronfile = TFile::Open(electronfilename);

   TFile *fileout = TFile::Open("TOF_finalSpectra.root", "RECREATE");
   TH1D *hRaw, *hMatch, *hTrack, *hPrim, *hElectron;
   Char_t title[1024];
   for (Int_t icent = 0; icent < NcentralityBins; icent++)
     for (Int_t icharge = 0; icharge < kNCharges; icharge++)
       for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	 hRaw = (TH1D *)rawfile->Get(Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hRaw) {
	   printf("cannot find %s in %s\n", Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), rawfilename);
	   continue;
	 }
	 hMatch = (TH1D *)matchfile->Get(Form("hMatchEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hMatch) {
	   printf("cannot find %s in %s\n", Form("hMatchEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), matchfilename);
	   continue;
	 }
	 hTrack = (TH1D *)trackfile->Get(Form("hTrackingEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hTrack) {
	   printf("cannot find %s in %s\n", Form("hTrackingEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), trackfilename);
	   continue;
	 }
	 hPrim = (TH1D *)primaryfile->Get(Form("hPrimaryFrac_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (ipart != 3 && !hPrim) {
	   printf("cannot find %s in %s\n", Form("hPrimaryFrac_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), primaryfilename);
	   continue;
	 }
	 hElectron = (TH1D *)electronfile->Get("hElectronCorr_average");
	 if (!hElectron) {
	   printf("cannot find hElectronCorr_average in %s\n", electronfilename);
	   continue;
	 }
	 hRaw->Divide(hMatch);
	 hRaw->Divide(hTrack);
	 if (hPrim) hRaw->Multiply(hPrim);
	 if (ipart == 2) hRaw->Multiply(hElectron);
	 sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partLatex[ipart][icharge], centralityBin[icent], centralityBin[icent + 1]);
	 hRaw->SetName(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 hRaw->SetTitle(title);
	 hRaw->SetMarkerStyle(marker[icharge]);
	 hRaw->SetMarkerColor(color[ipart]);
	 fileout->cd();
	 hRaw->Write();
       }

   fileout->Close();

   FinalRatios();
 }


/**************************************************************/

FinalRatios(const Char_t *filename = "TOF_finalSpectra.root")
{

  printf("WARNING: skipping ratios\n");
  return;

  Int_t marker[2] = {20, 25};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};

  /* open data */
  TFile *filein = TFile::Open(filename);
  /* get spectra */
  TH1D *hSpectrum[NcentralityBins][AliPID::kSPECIES][kNCharges];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hSpectrum[icent][ipart][icharge] = (TH1D *)filein->Get(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
      }
    }
  }

  /* output */
  TFile *fileout = TFile::Open("TOF_finalRatios.root", "RECREATE");

  /* particle/anti-particle ratios */
  TH1D *hRatio;
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      if (!hSpectrum[icent][ipart][kNegative] ||
	  !hSpectrum[icent][ipart][kPositive]) continue;
      hRatio = new TH1D(*hSpectrum[icent][ipart][kNegative]);
      hRatio->Divide(hSpectrum[icent][ipart][kPositive]);
      hRatio->SetName(Form("hRatio_cent%d_%s_negative_%s_positive", icent, AliPID::ParticleName(ipart), AliPID::ParticleName(ipart)));
      hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partLatex[ipart][kNegative], partLatex[ipart][kPositive], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partLatex[ipart][kNegative], partLatex[ipart][kPositive]));
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerColor(color[ipart]);
      fileout->cd();
      hRatio->Write();
    }
  }

  /* kaon/pion ratios */
  TH1D *hSum1, *hSum2;
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      if (!hSpectrum[icent][ipart][kNegative] ||
	  !hSpectrum[icent][AliPID::kPion][icharge]) continue;
      hRatio = new TH1D(*hSpectrum[icent][AliPID::kKaon][icharge]);
      hRatio->Divide(hSpectrum[icent][AliPID::kPion][icharge]);
      hRatio->SetName(Form("hRatio_cent%d_kaon_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
      hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partLatex[AliPID::kKaon][icharge], partLatex[AliPID::kPion][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partLatex[AliPID::kKaon][icharge], partLatex[AliPID::kPion][icharge]));
      hRatio->SetMarkerStyle(marker[icharge]);
      hRatio->SetMarkerColor(color[AliPID::kKaon]);
      fileout->cd();
      hRatio->Write();
    }
    if (!hSpectrum[icent][AliPID::kKaon][kPositive] ||
	!hSpectrum[icent][AliPID::kKaon][kNegative] ||
	!hSpectrum[icent][AliPID::kPion][kPositive] ||
	!hSpectrum[icent][AliPID::kPion][kNegative]) continue;
    hSum1 = new TH1D(*hSpectrum[icent][AliPID::kKaon][kPositive]);
    hSum1->Add(hSpectrum[icent][AliPID::kKaon][kNegative]);
    hSum2 = new TH1D(*hSpectrum[icent][AliPID::kPion][kPositive]);
    hSum2->Add(hSpectrum[icent][AliPID::kPion][kNegative]);
    hRatio = new TH1D(*hSum1);
    hRatio->Divide(hSum2);
    hRatio->SetName(Form("hRatio_cent%d_kaon_pion", icent));
    hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (%d-%d%%);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partLatex[AliPID::kKaon][kPositive], partLatex[AliPID::kKaon][kNegative], partLatex[AliPID::kPion][kPositive], partLatex[AliPID::kPion][kNegative], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partLatex[AliPID::kKaon][kPositive], partLatex[AliPID::kKaon][kNegative], partLatex[AliPID::kPion][kPositive], partLatex[AliPID::kPion][kNegative]));
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(color[AliPID::kKaon]);
    fileout->cd();
    hRatio->Write();
  }
  
  /* proton/pion ratios */
  TH1D *hSum1, *hSum2;
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      if (!hSpectrum[icent][AliPID::kProton][icharge] ||
	  !hSpectrum[icent][AliPID::kPion][icharge]) continue;
      hRatio = new TH1D(*hSpectrum[icent][AliPID::kProton][icharge]);
      hRatio->Divide(hSpectrum[icent][AliPID::kPion][icharge]);
      hRatio->SetName(Form("hRatio_cent%d_proton_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
      hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partLatex[AliPID::kProton][icharge], partLatex[AliPID::kPion][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partLatex[AliPID::kProton][icharge], partLatex[AliPID::kPion][icharge]));
      hRatio->SetMarkerStyle(marker[icharge]);
      hRatio->SetMarkerColor(color[AliPID::kProton]);
      fileout->cd();
      hRatio->Write();
    }
    if (!hSpectrum[icent][AliPID::kProton][kPositive] ||
	!hSpectrum[icent][AliPID::kProton][kNegative] ||
	!hSpectrum[icent][AliPID::kPion][kPositive] ||
	!hSpectrum[icent][AliPID::kPion][kNegative]) continue;
    hSum1 = new TH1D(*hSpectrum[icent][AliPID::kProton][kPositive]);
    hSum1->Add(hSpectrum[icent][AliPID::kProton][kNegative]);
    hSum2 = new TH1D(*hSpectrum[icent][AliPID::kPion][kPositive]);
    hSum2->Add(hSpectrum[icent][AliPID::kPion][kNegative]);
    hRatio = new TH1D(*hSum1);
    hRatio->Divide(hSum2);
    hRatio->SetName(Form("hRatio_cent%d_proton_pion", icent));
    hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (%d-%d%%);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partLatex[AliPID::kProton][kPositive], partLatex[AliPID::kProton][kNegative], partLatex[AliPID::kPion][kPositive], partLatex[AliPID::kPion][kNegative], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partLatex[AliPID::kProton][kPositive], partLatex[AliPID::kProton][kNegative], partLatex[AliPID::kPion][kPositive], partLatex[AliPID::kPion][kNegative]));
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(color[AliPID::kProton]);
    fileout->cd();
    hRatio->Write();
  }
  

  fileout->Close();
}

