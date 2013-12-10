#include "CommonDefs.C"

/**************************************************************/
Char_t *partChargeName[AliPID::kSPECIES][2] = {
  "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
};
 /**************************************************************/

FinalSpectra()
{

  doFinalSpectra(1., "TOF_rawSpectra.root", "TOF_matchingEfficiency.root", "TOF_trackingEfficiency.root", "TOF_primaryFraction.root", "TOF_electronCorrection.root", "TOF_finalSpectra.root");
  doFinalSpectra(1., "TOF_rawSpectra.root", "TOF_matchingEfficiency.root", "TOF_trackingEfficiency.root", "TOF_primaryFraction.root", "TOF_electronCorrection.root", "TOF_finalSpectra_centCorr.root", kFALSE);

  doFinalRatios("TOF_finalSpectra.root", "TOF_finalRatios.root");

  FinalSpectra_mismatchCorrected();
}

FinalSpectra_mismatchCorrected()
{
  doFinalSpectra(1., "TOF_rawSpectra_mismatchCorrected.root", "TOF_matchingEfficiency_mismatchCorrected.root", "TOF_trackingEfficiency.root", "TOF_primaryFraction.root", "TOF_electronCorrection.root", "TOF_finalSpectra_mismatchCorrected.root");
  doFinalSpectra(1., "TOF_rawSpectra_mismatchCorrected.root", "TOF_matchingEfficiency_mismatchCorrected.root", "TOF_trackingEfficiency.root", "TOF_primaryFraction.root", "TOF_electronCorrection.root", "TOF_finalSpectra_mismatchCorrected_centCorr.root", kFALSE);

  doFinalRatios("TOF_finalSpectra_mismatchCorrected.root", "TOF_finalRatios_mismatchCorrected.root");
}

doFinalSpectra(Double_t scaleFact = 1., const Char_t *rawfilename = "TOF_rawSpectra.root", const Char_t *matchfilename = "TOF_matchingEfficiency.root", const Char_t *trackfilename = "TOF_trackingEfficiency.root", const Char_t *primaryfilename = "TOF_primaryFraction.root", const Char_t *electronfilename = "TOF_electronCorrection.root", const Char_t *outputfilename = "TOF_finalSpectra.root", Bool_t useMB = kTRUE)
 {

   Int_t marker[2] = {20, 25};
   Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};

   TFile *rawfile = TFile::Open(rawfilename);
   TFile *matchfile = TFile::Open(matchfilename);
   TFile *trackfile = TFile::Open(trackfilename);
   TFile *primaryfile = TFile::Open(primaryfilename);
   TFile *electronfile = TFile::Open(electronfilename);

   TFile *fileout = TFile::Open(outputfilename, "RECREATE");
   TH1D *hRaw, *hMatch, *hTrack, *hPrim, *hElectron;
   TH1D *hSpectrum[NcentralityBins+1][5][2];
   Char_t title[1024];
   for (Int_t icent = -1; icent < NcentralityBins; icent++)
     for (Int_t icharge = 0; icharge < kNCharges; icharge++)
       for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	 if (icent == -1)
	   hRaw = (TH1D *)rawfile->Get(Form("hRaw_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	 else
	   hRaw = (TH1D *)rawfile->Get(Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hRaw) {
	   printf("cannot find %s in %s\n", Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), rawfilename);
	   continue;
	 }
	 if (0)
	   hMatch = (TH1D *)matchfile->Get(Form("hMatchEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 else
	   hMatch = (TH1D *)matchfile->Get(Form("hMatchEff_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hMatch) {
	   printf("cannot find %s in %s\n", Form("hMatchEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), matchfilename);
	   continue;
	 }
	 if (0)
	   hTrack = (TH1D *)trackfile->Get(Form("hTrackingEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 else
	   hTrack = (TH1D *)trackfile->Get(Form("hTrackingEff_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	 if (!hTrack) {
	   printf("cannot find %s in %s\n", Form("hTrackingEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), trackfilename);
	   continue;
	 }
	 if (useMB || icent == -1) {
	   hPrim = (TH1D *)primaryfile->Get(Form("hPrimaryFrac_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	   if (ipart != 3 && !hPrim) {
	     printf("cannot find %s in %s\n", Form("hPrimaryFrac_MB_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), primaryfilename);
	     continue;
	   }
	 }
	 else {
	   hPrim = (TH1D *)primaryfile->Get(Form("hPrimaryFrac_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	   if (ipart != 3 && !hPrim) {
	     printf("cannot find %s in %s\n", Form("hPrimaryFrac_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]), primaryfilename);
	     continue;
	   }
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
	 if (icent == -1) 
	   sprintf(title, "%s (MB);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partChargeName[ipart][icharge]);
	 else
	   sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partChargeName[ipart][icharge], centralityBin[icent], centralityBin[icent + 1]);
	 if (icent == -1)
	   hRaw->SetName(Form("hFinal_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	 else
	   hRaw->SetName(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 hRaw->SetTitle(title);
	 hRaw->SetMarkerStyle(marker[icharge]);
	 hRaw->SetMarkerColor(color[ipart]);
	 fileout->cd();
	 hRaw->Scale(scaleFact);
	 hRaw->Write();
	 if (icent == -1)
	   hSpectrum[NcentralityBins][ipart][icharge] = hRaw;
	 else
	   hSpectrum[icent][ipart][icharge] = hRaw;
       }

   TH1 *hRatioMB;
   for (Int_t icent = 0; icent < NcentralityBins; icent++)
     for (Int_t icharge = 0; icharge < kNCharges; icharge++)
       for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	 hRatioMB = (TH1 *)hSpectrum[icent][ipart][icharge]->Clone(Form("hRatioMB_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	 hRatioMB->Divide(hRatioMB, hSpectrum[NcentralityBins][ipart][icharge], 1., 1., "");
	 fileout->cd();
	 hRatioMB->Write();
       }
   
   fileout->Close();

 }


/**************************************************************/

doFinalRatios(const Char_t *filename = "TOF_finalSpectra.root", const Char_t *outputfilename = "TOF_finalRatios.root")
{

  //  printf("WARNING: skipping ratios\n");
  //  return;

  Int_t marker[2] = {20, 25};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};

  /* open data */
  TFile *filein = TFile::Open(filename);
  /* get spectra */
  TH1D *hSpectrum[NcentralityBins+1][AliPID::kSPECIES][kNCharges];
  for (Int_t icent = 0; icent < NcentralityBins+1; icent++) {
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	if (icent == NcentralityBins)
	  hSpectrum[icent][ipart][icharge] = (TH1D *)filein->Get(Form("hFinal_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	else
	  hSpectrum[icent][ipart][icharge] = (TH1D *)filein->Get(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
      }
    }
  }

  /* output */
  TFile *fileout = TFile::Open(outputfilename, "RECREATE");

  /* particle/anti-particle ratios */
  TH1D *hRatio;
  TH1D *hPosNegRatio[NcentralityBins+1][5];
  for (Int_t icent = 0; icent < NcentralityBins+1; icent++) {
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      if (!hSpectrum[icent][ipart][kNegative] ||
	  !hSpectrum[icent][ipart][kPositive]) continue;
      hRatio = new TH1D(*hSpectrum[icent][ipart][kNegative]);
      hRatio->Divide(hSpectrum[icent][ipart][kPositive]);
      if (icent == NcentralityBins) {
	hRatio->SetName(Form("hRatio_MB_%s_negative_%s_positive", AliPID::ParticleName(ipart), AliPID::ParticleName(ipart)));
	hRatio->SetTitle(Form("%s/%s (MB);p_{T} (GeV/c);%s/%s;", partChargeName[ipart][kNegative], partChargeName[ipart][kPositive], partChargeName[ipart][kNegative], partChargeName[ipart][kPositive]));
      }
      else {
	hRatio->SetName(Form("hRatio_cent%d_%s_negative_%s_positive", icent, AliPID::ParticleName(ipart), AliPID::ParticleName(ipart)));
	hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partChargeName[ipart][kNegative], partChargeName[ipart][kPositive], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partChargeName[ipart][kNegative], partChargeName[ipart][kPositive]));
      }
      hRatio->SetMarkerStyle(20);
      hRatio->SetMarkerColor(color[ipart]);
      fileout->cd();
      hRatio->Write();
      hPosNegRatio[icent][ipart] = hRatio;
    }
  }
  
  TH1 *hRatioMB;
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      hRatioMB = (TH1 *)hPosNegRatio[icent][ipart]->Clone(Form("hDoubleRatioMB_cent%d_%s_negative_%s_positive", icent, AliPID::ParticleName(ipart), AliPID::ParticleName(ipart)));
      hRatioMB->Divide(hRatioMB, hPosNegRatio[NcentralityBins][ipart], 1., 1., "");
      fileout->cd();
      hRatioMB->Write();
    }
  
  /* kaon/pion ratios */
  TH1D *hSum1, *hSum2;
  TH1D *hKaToPi[NcentralityBins+1][3];
  for (Int_t icent = 0; icent < NcentralityBins+1; icent++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      if (!hSpectrum[icent][AliPID::kKaon][icharge] ||
	  !hSpectrum[icent][AliPID::kPion][icharge]) continue;
      hRatio = new TH1D(*hSpectrum[icent][AliPID::kKaon][icharge]);
      hRatio->Divide(hSpectrum[icent][AliPID::kPion][icharge]);
      if (icent == NcentralityBins) {
	hRatio->SetName(Form("hRatio_MB_kaon_%s_pion_%s", chargeName[icharge], chargeName[icharge]));
	hRatio->SetTitle(Form("%s/%s (MB);p_{T} (GeV/c);%s/%s;", partChargeName[AliPID::kKaon][icharge], partChargeName[AliPID::kPion][icharge], partChargeName[AliPID::kKaon][icharge], partChargeName[AliPID::kPion][icharge]));
      }
      else {
	hRatio->SetName(Form("hRatio_cent%d_kaon_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
	hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partChargeName[AliPID::kKaon][icharge], partChargeName[AliPID::kPion][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partChargeName[AliPID::kKaon][icharge], partChargeName[AliPID::kPion][icharge]));
      }
      hRatio->SetMarkerStyle(marker[icharge]);
      hRatio->SetMarkerColor(color[AliPID::kKaon]);
      fileout->cd();
      hRatio->Write();
      hKaToPi[icent][icharge] = hRatio;
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
    if (icent == NcentralityBins) {
      hRatio->SetName(Form("hRatio_MB_kaon_pion"));
      hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (MB);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partChargeName[AliPID::kKaon][kPositive], partChargeName[AliPID::kKaon][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative], partChargeName[AliPID::kKaon][kPositive], partChargeName[AliPID::kKaon][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative]));
    }
    else {
      hRatio->SetName(Form("hRatio_cent%d_kaon_pion", icent));
      hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (%d-%d%%);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partChargeName[AliPID::kKaon][kPositive], partChargeName[AliPID::kKaon][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partChargeName[AliPID::kKaon][kPositive], partChargeName[AliPID::kKaon][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative]));
    }
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(color[AliPID::kKaon]);
    fileout->cd();
    hRatio->Write();
    hKaToPi[icent][2] = hRatio;
  }

  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hRatioMB = (TH1 *)hKaToPi[icent][icharge]->Clone(Form("hDoubleRatioMB_cent%d_kaon_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
      hRatioMB->Divide(hRatioMB, hKaToPi[NcentralityBins][icharge], 1., 1., "");
      fileout->cd();
      hRatioMB->Write();
    }
    hRatioMB = (TH1 *)hKaToPi[icent][2]->Clone(Form("hDoubleRatioMB_cent%d_kaon_pion", icent));
    hRatioMB->Divide(hRatioMB, hKaToPi[NcentralityBins][2], 1., 1., "");
    fileout->cd();
    hRatioMB->Write();
  }

  
  /* proton/pion ratios */
  TH1D *hSum1, *hSum2;
  TH1D *hPrToPi[NcentralityBins+1][3];
  for (Int_t icent = 0; icent < NcentralityBins+1; icent++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      if (!hSpectrum[icent][AliPID::kProton][icharge] ||
	  !hSpectrum[icent][AliPID::kPion][icharge]) continue;
      hRatio = new TH1D(*hSpectrum[icent][AliPID::kProton][icharge]);
      hRatio->Divide(hSpectrum[icent][AliPID::kPion][icharge]);
      if (icent == NcentralityBins) {
	hRatio->SetName(Form("hRatio_MB_proton_%s_pion_%s", chargeName[icharge], chargeName[icharge]));
	hRatio->SetTitle(Form("%s/%s (MB);p_{T} (GeV/c);%s/%s;", partChargeName[AliPID::kProton][icharge], partChargeName[AliPID::kPion][icharge], partChargeName[AliPID::kProton][icharge], partChargeName[AliPID::kPion][icharge]));
      }
      else {
	hRatio->SetName(Form("hRatio_cent%d_proton_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
	hRatio->SetTitle(Form("%s/%s (%d-%d%%);p_{T} (GeV/c);%s/%s;", partChargeName[AliPID::kProton][icharge], partChargeName[AliPID::kPion][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partChargeName[AliPID::kProton][icharge], partChargeName[AliPID::kPion][icharge]));
      }
      hRatio->SetMarkerStyle(marker[icharge]);
      hRatio->SetMarkerColor(color[AliPID::kProton]);
      fileout->cd();
      hRatio->Write();
      hPrToPi[icent][icharge] = hRatio;
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
    if (icent == NcentralityBins) {
      hRatio->SetName(Form("hRatio_MB_proton_pion"));
      hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (MB);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partChargeName[AliPID::kProton][kPositive], partChargeName[AliPID::kProton][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative], partChargeName[AliPID::kProton][kPositive], partChargeName[AliPID::kProton][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative]));
    }
    else {
      hRatio->SetName(Form("hRatio_cent%d_proton_pion", icent));
      hRatio->SetTitle(Form("(%s+%s)/(%s+%s) (%d-%d%%);p_{T} (GeV/c);(%s+%s)/(%s+%s);", partChargeName[AliPID::kProton][kPositive], partChargeName[AliPID::kProton][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1], partChargeName[AliPID::kProton][kPositive], partChargeName[AliPID::kProton][kNegative], partChargeName[AliPID::kPion][kPositive], partChargeName[AliPID::kPion][kNegative]));
    }
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(color[AliPID::kProton]);
    fileout->cd();
    hRatio->Write();
    hPrToPi[icent][2] = hRatio;
  }
  
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      hRatioMB = (TH1 *)hPrToPi[icent][icharge]->Clone(Form("hDoubleRatioMB_cent%d_proton_%s_pion_%s", icent, chargeName[icharge], chargeName[icharge]));
      hRatioMB->Divide(hRatioMB, hPrToPi[NcentralityBins][icharge], 1., 1., "");
      fileout->cd();
      hRatioMB->Write();
    }
    hRatioMB = (TH1 *)hPrToPi[icent][2]->Clone(Form("hDoubleRatioMB_cent%d_proton_pion", icent));
    hRatioMB->Divide(hRatioMB, hPrToPi[NcentralityBins][2], 1., 1., "");
    fileout->cd();
    hRatioMB->Write();
  }


  fileout->Close();
}

