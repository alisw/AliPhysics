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

SystematicChecks()
{
  SystematicCheck();
  SystematicPlots_spectra();
}


SystematicPlots_spectra()
{


  enum EData_t {
    kuseTPCcrossedRows,
    kminTPCclusters_60,
    kminTPCclusters_80, 
    kscaleDCAxy_09,
    kscaleDCAxy_11,
    kscaleDCAxy_10x,
    ketaCut_75,
    ketaCut_85,
    kbkgFit_fixed_scaleSigma_09,
    kbkgFit_fixed_scaleSigma_11,
    kbkgFit_fixed_scaleTail_09,
    kbkgFit_fixed_scaleTail_11,
    kbkgFit_fixed_scaleSigma_09_scaleTail_11,
    kbkgFit_fixed_scaleSigma_11_scaleTail_09,
    ksignalFit_fixed_scaleSigma_09,
    ksignalFit_fixed_scaleSigma_11,
    ksignalFit_fixed_scaleTail_09,
    ksignalFit_fixed_scaleTail_11,
    ksignalFit_fixed_scaleSigma_09_scaleTail_11,
    ksignalFit_fixed_scaleSigma_11_scaleTail_09,
    kbkgFit_fixed,
    kbkgFit_free,
    ksignalFit_fixed,
    ksignalFit_free,
    ksignalbkgFit_free,
    kdefaultFit_fitElectrons,
    kdefaultFit_limitedRange,
    kfieldReversal
  };
  const Int_t ndata = 28;
  const Char_t *name[ndata] = {
    "useTPCcrossedRows",
    "minTPCclusters_60",
    "minTPCclusters_80", 
    "scaleDCAxy_09",
    "scaleDCAxy_11",
    "scaleDCAxy_10x",
    "etaCut_75",
    "etaCut_85",
    "bkgFit_fixed_scaleSigma_09",
    "bkgFit_fixed_scaleSigma_11",
    "bkgFit_fixed_scaleTail_09",
    "bkgFit_fixed_scaleTail_11",
    "bkgFit_fixed_scaleSigma_09_scaleTail_11",
    "bkgFit_fixed_scaleSigma_11_scaleTail_09",
    "signalFit_fixed_scaleSigma_09",
    "signalFit_fixed_scaleSigma_11",
    "signalFit_fixed_scaleTail_09",
    "signalFit_fixed_scaleTail_11",
    "signalFit_fixed_scaleSigma_09_scaleTail_11",
    "signalFit_fixed_scaleSigma_11_scaleTail_09",
    "bkgFit_fixed",
    "bkgFit_free",
    "signalFit_fixed",
    "signalFit_free",
    "signalbkgFit_free",
    "defaultFit_fitElectrons",
    "defaultFit_limitedRange",
    "fieldReversal"
  };
  Int_t marker[ndata] = {
    25, 20, 21,
    20, 21, 25,
    20, 21,
    20, 21, 20, 21, 20, 21,
    20, 21, 20, 21, 20, 21,
    20, 21, 20, 21, 25,
    20,
    20,
    20
  };
  Int_t color[ndata] = {
    4, 2, 8,
    2, 8, 4,
    2, 8,
    2, 2, 4, 4, 8, 8,
    2, 2, 4, 4, 8, 8,
    2, 8, 4,
    4,
    4,
    4
  };

  TFile *fin[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) 
    fin[idata] = TFile::Open(Form("SystematicCheck_finalSpectra_%s.root", name[idata]));

  TFile *fileout = TFile::Open("SystematicPlots_spectra.root", "RECREATE");
    
  TH1D *hin[ndata];
  TH1D *hArea = new TH1D("hArea", "", NptBins, ptBin);
  hArea->SetStats(kFALSE);
  Char_t title[1024];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	for (Int_t idata = 0; idata < ndata; idata++) {
	  if (!fin[idata] || !fin[idata]->IsOpen()) continue;
	  hin[idata] = (TH1D *)fin[idata]->Get(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	  if (!hin[idata]) {
	    printf("cannot find hFinal_cent%d_%s_%s in %s\n", icent, AliPID::ParticleName(ipart), chargeName[icharge], fin[idata]->GetName());
	    continue;
	  }
	  hin[idata]->SetMarkerStyle(marker[idata]);
	  hin[idata]->SetMarkerColor(color[idata]);
	}

	/* set common title */
	sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);relative variation wrt. standard;", partLatex[ipart][icharge], centralityBin[icent], centralityBin[icent + 1]);
	hArea->SetTitle(title);

	/* TPC quality cuts */
	hArea->SetMinimum(-0.05);
	hArea->SetMaximum(0.05);
	hArea->Draw();
	hin[kuseTPCcrossedRows]->Draw("same");
	hin[kminTPCclusters_60]->Draw("same");
	hin[kminTPCclusters_80]->Draw("same");
	TLegend *l = gPad->BuildLegend();
	l->DeleteEntry();
	l->SetFillColor(0);
	l->SetBorderSize(1);
	fileout->cd();
	gPad->Write(Form("cSpectraSys_TPCcuts_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

	/* DCA cuts */
	hArea->SetMinimum(-0.05);
	hArea->SetMaximum(0.05);
	hArea->Draw();
	hin[kscaleDCAxy_09]->Draw("same");
	hin[kscaleDCAxy_11]->Draw("same");
	hin[kscaleDCAxy_10x]->Draw("same");
	TLegend *l = gPad->BuildLegend();
	l->DeleteEntry();
	l->SetFillColor(0);
	l->SetBorderSize(1);
	fileout->cd();
	gPad->Write(Form("cSpectraSys_DCAcuts_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

	/* eta cut */
	hArea->SetMinimum(-0.05);
	hArea->SetMaximum(0.05);
	hArea->Draw();
	hin[ketaCut_75]->Draw("same");
	hin[ketaCut_85]->Draw("same");
	TLegend *l = gPad->BuildLegend();
	l->DeleteEntry();
	l->SetFillColor(0);
	l->SetBorderSize(1);
	fileout->cd();
	gPad->Write(Form("cSpectraSys_etaCut_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

	/* background fit */
	hArea->SetMinimum(-0.5);
	hArea->SetMaximum(0.5);
	hArea->Draw();
	hin[kbkgFit_fixed_scaleSigma_09]->Draw("same");
	hin[kbkgFit_fixed_scaleSigma_11]->Draw("same");
	hin[kbkgFit_fixed_scaleTail_09]->Draw("same");
	hin[kbkgFit_fixed_scaleTail_11]->Draw("same");
	hin[kbkgFit_fixed_scaleSigma_09_scaleTail_11]->Draw("same");
	hin[kbkgFit_fixed_scaleSigma_11_scaleTail_09]->Draw("same");
	TLegend *l = gPad->BuildLegend();
	l->DeleteEntry();
	l->SetFillColor(0);
	l->SetBorderSize(1);
	fileout->cd();
	gPad->Write(Form("cSpectraSys_bkgFit_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));

	/* signal fit */
	hArea->SetMinimum(-0.5);
	hArea->SetMaximum(0.5);
	hArea->Draw();
	hin[ksignalFit_fixed_scaleSigma_09]->Draw("same");
	hin[ksignalFit_fixed_scaleSigma_11]->Draw("same");
	hin[ksignalFit_fixed_scaleTail_09]->Draw("same");
	hin[ksignalFit_fixed_scaleTail_11]->Draw("same");
	hin[ksignalFit_fixed_scaleSigma_09_scaleTail_11]->Draw("same");
	hin[ksignalFit_fixed_scaleSigma_11_scaleTail_09]->Draw("same");
	TLegend *l = gPad->BuildLegend();
	l->DeleteEntry();
	l->SetFillColor(0);
	l->SetBorderSize(1);
	fileout->cd();
	gPad->Write(Form("cSpectraSys_signalFit_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
      }
    }
  }

  fileout->Close();
  for (Int_t idata = 0; idata < ndata; idata++) 
    if (fin[idata])
    fin[idata]->Close();
}

SystematicCheck()
{
  SystematicCheck("rawSpectra");
  SystematicCheck("matchingEfficiency");
  SystematicCheck("trackingEfficiency");
  SystematicCheck("primaryFraction");
  SystematicCheck("finalSpectra");
  SystematicCheck("finalRatios");
}

SystematicCheck(const Char_t *checkname)
{

  gROOT->LoadMacro("HistoUtils.C");

  const Int_t ndata = 30;
  const Char_t *name[ndata] = {
    "useTPCcrossedRows",
    "minTPCclusters_60",
    "minTPCclusters_80",
    "scaleDCAxy_09",
    "scaleDCAxy_11",
    "scaleDCAxy_10x",
    "DCAz_1",
    "etaCut_75",
    "etaCut_85",
    "bkgFit_fixed",
    "bkgFit_fixed_scaleSigma_09",
    "bkgFit_fixed_scaleSigma_09_scaleTail_11",
    "bkgFit_fixed_scaleSigma_11",
    "bkgFit_fixed_scaleSigma_11_scaleTail_09",
    "bkgFit_fixed_scaleTail_09",
    "bkgFit_fixed_scaleTail_11",
    "bkgFit_free",
    "mismatchCorrected",
    "defaultFit_fitElectrons",
    "defaultFit_limitedRange",
    "defaultFit_tightRange",
    "fieldReversal",
    "signalbkgFit_free",
    "signalFit_fixed_scaleSigma_09",
    "signalFit_fixed_scaleSigma_09_scaleTail_11",
    "signalFit_fixed_scaleSigma_11",
    "signalFit_fixed_scaleSigma_11_scaleTail_09",
    "signalFit_fixed_scaleTail_09",
    "signalFit_fixed_scaleTail_11",
    "signalFit_free"
  };


  const Char_t *title[ndata] = {
    "N_{crossed-rows} cut",
    "N_{clusters}^{min} = 60",
    "N_{clusters}^{min} = 80",
    "-10% DCA_{xy} cut",
    "+10% DCA_{xy} cut",
    "10x larger DCA_{xy} cut",
    "DCA_{z} < 1 cm",
    "|#eta| < 0.75",
    "|#eta| < 0.85",
    "background fit (fixed params)",
    "-10% #sigma background",
    "-10% #sigma +10% #tau background",
    "+10% #sigma background",
    "+10% #sigma -10% #tau background",
    "-10% #tau background",
    "+10% #tau background",
    "background fit (free params)",
    "with mismatch correction",
    "with electron background",
    "fit in limited range",
    "fit in tight range",
    "reversed magnetic field",
    "signal+background fit (free params)",
    "-10% #sigma signal",
    "-10% #sigma +10% #tau signal",
    "+10% #sigma signal",
    "+10% #sigma -10% #tau signal",
    "-10% #tau signal",
    "+10% #tau signal",
    "signal fit (free params)"
  };

#if 0
  Char_t filename1[1024], filename2[1024], ofilename[1024];
  for (Int_t idata = 1; idata < ndata; idata++) {
    sprintf(filename1, "standardPrimaryCuts_%s/TOF_%s.root", name[idata], checkname);
    sprintf(filename2, "standardPrimaryCuts/TOF_%s.root", checkname);
    sprintf(ofilename, "SystematicCheck_%s_%s.root", checkname, name[idata]);
    printf("%s %s %s\n", filename1, filename2, ofilename);
    HistoUtils_autosystematics(filename1, filename2, ofilename, Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[idata]));
  }
  return;
#endif
  
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[0], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[0]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[0]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[1], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[1]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[1]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[2], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[2]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[2]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[3], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[3]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[3]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[4], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[4]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[4]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[5], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[5]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[5]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[6], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[6]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[6]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[7], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[7]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[7]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[8], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[8]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[8]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[9], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[9]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[9]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[10], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[10]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[10]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[11], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[11]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[11]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[12], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[12]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[12]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[13], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[13]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[13]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[14], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[14]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[14]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[15], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[15]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[15]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[16], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[16]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[16]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[17], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[17]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[17]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[18], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[18]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[18]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[19], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[19]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[19]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[20], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[20]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[20]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[21], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[21]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[21]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[22], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[22]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[22]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[23], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[23]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[23]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[24], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[24]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[24]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[25], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[25]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[25]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[26], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[26]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[26]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[26], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[27]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[27]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[26], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[28]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[28]));
  HistoUtils_autosystematics(Form("standardPrimaryCuts_%s/TOF_%s.root", name[26], checkname), Form("standardPrimaryCuts/TOF_%s.root", checkname), Form("SystematicCheck_%s_%s.root", checkname, name[29]), Form("%s;p_{T} (GeV/c);relative variation wrt. default;", title[29]));

  
}

/**************************************************************/
/**************************************************************/

FinalSpectra_systematics_useTPCcrossedRows()
{
  
  TFile *fileout = TFile::Open("FinalSpectra_systematics_useTPCcrossedRows.root", "RECREATE");
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	FinalSpectra_systematics_useTPCcrossedRows(ipart, icharge, icent, fileout);
      }
  fileout->Close();
}

FinalSpectra_systematics_useTPCcrossedRows(Int_t ipart, Int_t icharge, Int_t icent, TFile *fileout = NULL)
{
  const Int_t ndata = 1;
  const Char_t *name[ndata] = {
    "standardPrimaryCuts_useTPCcrossedRows"
  };
  const Char_t *title[ndata] = {
    "TPC crossed-rows cut;p_{T} (GeV/c); ratio"
  };
  Int_t marker[ndata] = {22};
  Int_t color[ndata] = {2};
  
  TH1D *hArea = new TH1D("hArea", Form("%s (%d-%d%%);p_{T} (GeV/c);ratio;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]), NptBins, ptBin);
  hArea->SetMinimum(0.9);
  hArea->SetMaximum(1.1); 
  hArea->SetStats(kFALSE);
  hArea->Draw();
  
  TH1D *hr[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) {
    hr[idata] = FinalSpectra_systematics_ratio(name[idata], "standardPrimaryCuts", ipart, icharge, icent, name[idata], title[idata], marker[idata], color[idata]);
    hr[idata]->Draw("same");
  }
  
  gPad->SetGridy();
  TLegend *legend = gPad->BuildLegend();
  legend->DeleteEntry();
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  if (fileout) {
    fileout->cd();
    gPad->Write(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
  }
}

FinalSpectra_systematics_minTPCclusters()
{
  
  TFile *fileout = TFile::Open("FinalSpectra_systematics_minTPCclusters.root", "RECREATE");
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	FinalSpectra_systematics_minTPCclusters(ipart, icharge, icent, fileout);
      }
  fileout->Close();
}

FinalSpectra_systematics_minTPCclusters(Int_t ipart, Int_t icharge, Int_t icent, TFile *fileout = NULL)
{
  const Int_t ndata = 2;
  const Char_t *name[ndata] = {
    "standardPrimaryCuts_minTPCclusters_60",
    "standardPrimaryCuts_minTPCclusters_80"
  };
  const Char_t *title[ndata] = {
    "N_{TPC-cls}^{min} = 60;p_{T} (GeV/c); ratio",
    "N_{TPC-cls}^{min} = 80;p_{T} (GeV/c); ratio"
  };
  Int_t marker[ndata] = {22, 28};
  Int_t color[ndata] = {2, 8};

  TH1D *hArea = new TH1D("hArea", Form("%s (%d-%d%%);p_{T} (GeV/c);ratio;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]), NptBins, ptBin);
  hArea->SetMinimum(0.9);
  hArea->SetMaximum(1.1);
  hArea->SetStats(kFALSE);
  hArea->Draw();
  
  TH1D *hr[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) {
    hr[idata] = FinalSpectra_systematics_ratio(name[idata], "standardPrimaryCuts", ipart, icharge, icent, name[idata], title[idata], marker[idata], color[idata]);
    hr[idata]->Draw("same");
  }

  gPad->SetGridy();
  TLegend *legend = gPad->BuildLegend();
  legend->DeleteEntry();
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  if (fileout) {
    fileout->cd();
    gPad->Write(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
  }
}

FinalSpectra_systematics_scaleDCAxy()
{
  
  TFile *fileout = TFile::Open("FinalSpectra_systematics_scaleDCAxy.root", "RECREATE");
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	FinalSpectra_systematics_scaleDCAxy(ipart, icharge, icent, fileout);
      }
  fileout->Close();
}

FinalSpectra_systematics_scaleDCAxy(Int_t ipart, Int_t icharge, Int_t icent, TFile *fileout = NULL)
{
  const Int_t ndata = 2;
  const Char_t *name[ndata] = {
    "standardPrimaryCuts_scaleDCAxy_09",
    "standardPrimaryCuts_scaleDCAxy_11"
  };
  const Char_t *title[ndata] = {
    "-10% DCA_{xy} cut;p_{T} (GeV/c); ratio",
    "+10% DCA_{xy} cut;p_{T} (GeV/c); ratio"
  };
  Int_t marker[ndata] = {22, 28};
  Int_t color[ndata] = {2, 8};

  TH1D *hArea = new TH1D("hArea", Form("%s (%d-%d%%);p_{T} (GeV/c);ratio;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]), NptBins, ptBin);
  hArea->SetMinimum(0.9);
  hArea->SetMaximum(1.1);
  hArea->SetStats(kFALSE);
  hArea->Draw();
  
  TH1D *hr[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) {
    hr[idata] = FinalSpectra_systematics_ratio(name[idata], "standardPrimaryCuts", ipart, icharge, icent, name[idata], title[idata], marker[idata], color[idata]);
    hr[idata]->Draw("same");
  }

  gPad->SetGridy();
  TLegend *legend = gPad->BuildLegend();
  legend->DeleteEntry();
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  if (fileout) {
    fileout->cd();
    gPad->Write(Form("cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
  }
}

TH1D *
FinalSpectra_systematics_ratio(const Char_t *dirname1, const Char_t *dirname2, Int_t ipart, Int_t icharge, Int_t icent, const Char_t *name = "finalRatio", const Char_t *title = ";p_{T} (GeV/c);final yield ratio;", Int_t marker = 20, Int_t color = 2, Bool_t correlated = kFALSE)
{

  TH1D *hr = new TH1D("hr", "", NptBins, ptBin);

  /* open data */
  Char_t outfilename1[1024];
  sprintf(outfilename1, "%s/TOF_finalSpectra.root", dirname1);
  TFile *filein1 = TFile::Open(outfilename1);
  if (!filein1 || !filein1->IsOpen()) {
    printf("cannot open %s\n", outfilename1);
    return;
  }
  Char_t outfilename2[1024];
  sprintf(outfilename2, "%s/TOF_finalSpectra.root", dirname2);
  TFile *filein2 = TFile::Open(outfilename2);
  if (!filein2 || !filein2->IsOpen()) {
    printf("cannot open %s\n", outfilename2);
    return;
  }
  /* get data */
  TH1D *h1 = (TH1D *)filein1->Get(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
  if (!h1) {
    printf("cannot get hFinal_cent%d_%s_%s from %s\n", icent, AliPID::ParticleName(ipart), chargeName[icharge], outfilename1);
    return;
  }
  TH1D *h2 = (TH1D *)filein2->Get(Form("hFinal_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
  if (!h2) {
    printf("cannot get hFinal_cent%d_%s_%s from %s\n", icent, AliPID::ParticleName(ipart), chargeName[icharge], outfilename2);
    return;
  }
  /* ratio */
  if (correlated) hr->Divide(h1, h2, 1., 1., "B");
  else hr->Divide(h1, h2);
  hr->SetNameTitle(name, title);
  hr->SetMarkerStyle(marker);
  hr->SetMarkerColor(color);
  
  filein1->Close();
  filein2->Close();
  
  return hr;
}
