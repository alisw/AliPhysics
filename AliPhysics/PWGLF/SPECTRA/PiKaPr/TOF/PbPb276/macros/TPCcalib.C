Double_t tofReso = 85.;

/**************************************************************/
/*** HISTOS AND BINNING ***************************************/
/**************************************************************/

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
enum EHistoParam_t {
  kCentrality,
  kPtpc,
  kTPCNcls,
  kTPCdelta,
  kTPCgain,
  kTOFsignal,
  kNHistoParams
};
/**************************************************************/
const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
/**************************************************************/
const Int_t NpBins = 46;
Double_t pBin[NpBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
/**************************************************************/
const Int_t NclsBins = 200;
Double_t clsBin[NclsBins + 1];
Double_t clsMin = 0., clsMax = 200., clsStep = (clsMax - clsMin) / NclsBins;
/**************************************************************/
const Int_t NtpcdeltaBins = 200;
Double_t tpcdeltaBin[NtpcdeltaBins + 1];
Double_t tpcdeltaMin = -100., tpcdeltaMax = 100., tpcdeltaStep = (tpcdeltaMax - tpcdeltaMin) / NtpcdeltaBins;
/**************************************************************/
const Int_t NtpcgainBins = 200;
Double_t tpcgainBin[NtpcgainBins + 1];
Double_t tpcgainMin = 0., tpcgainMax = 2., tpcgainStep = (tpcgainMax - tpcgainMin) / NtpcgainBins;
/**************************************************************/
const Int_t NsigmaBins = 20;
Double_t sigmaBin[NsigmaBins + 1];
Double_t sigmaMin = -5., sigmaMax = 5., sigmaStep = (sigmaMax - sigmaMin) / NsigmaBins;
/**************************************************************/
Int_t NparamsBins[kNHistoParams] = {NcentralityBins, NpBins, NclsBins, NtpcdeltaBins, NtpcgainBins, NsigmaBins};
Double_t *paramBin[kNHistoParams] = {centralityBin, pBin, clsBin, tpcdeltaBin, tpcgainBin, sigmaBin};
/**************************************************************/

/**************************************************************/


TPCcalib(const Char_t *filename, Int_t evMax = kMaxInt, Int_t startEv = 0)
{
  
  /* include path for ACLic */
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  /* load libraries */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  /* build analysis task class */
  gROOT->LoadMacro("AliAnalysisParticle.cxx+g");
  gROOT->LoadMacro("AliAnalysisEvent.cxx+g");
  gROOT->LoadMacro("AliAnalysisTrack.cxx+g");

  /* open file, get tree and connect */
  TFile *filein = TFile::Open(filename);
  TTree *treein = (TTree *)filein->Get("aodTree");
  printf("got \"aodTree\": %d entries\n", treein->GetEntries());
  AliAnalysisEvent *analysisEvent = new AliAnalysisEvent();
  TClonesArray *analysisTrackArray = new TClonesArray("AliAnalysisTrack");
  AliAnalysisTrack *analysisTrack = NULL;
  treein->SetBranchAddress("AnalysisEvent", &analysisEvent);
  treein->SetBranchAddress("AnalysisTrack", &analysisTrackArray);

  /**************************************************************/
  /*** HISTOS ***************************************************/
  /**************************************************************/

  /* run-time binning */
  for (Int_t ibin = 0; ibin < NclsBins + 1; ibin++)
    clsBin[ibin] = clsMin + ibin * clsStep;
  for (Int_t ibin = 0; ibin < NtpcdeltaBins + 1; ibin++)
    tpcdeltaBin[ibin] = tpcdeltaMin + ibin * tpcdeltaStep;
  for (Int_t ibin = 0; ibin < NtpcgainBins + 1; ibin++)
    tpcgainBin[ibin] = tpcgainMin + ibin * tpcgainStep;
  for (Int_t ibin = 0; ibin < NsigmaBins + 1; ibin++)
    sigmaBin[ibin] = sigmaMin + ibin * sigmaStep;

  /* THnSparse */
  THnSparse *hTPCcalib[AliPID::kSPECIES][kNCharges];
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge< kNCharges; icharge++) {
      hTPCcalib[ipart][icharge] = new THnSparseF(Form("hTPCcalib_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", kNHistoParams, NparamsBins);
      for (Int_t iparam = 0; iparam < kNHistoParams; iparam++) {
	hTPCcalib[ipart][icharge]->SetBinEdges(iparam, paramBin[iparam]);
      }
    }
  }

  /**************************************************************/
  /**************************************************************/
  /**************************************************************/

  /* TOF PID response */
  AliTOFPIDResponse tofResponse;
  tofResponse.SetTimeResolution(tofReso);
  /* TPC PID response */
  AliTPCPIDResponse *tpcResponse = AliAnalysisTrack::GetTPCResponse();

  /* start stopwatch */
  TStopwatch timer;
  timer.Start();

  /* loop over events */
  Bool_t hastofpid;
  Int_t charge, index;
  UShort_t dedxN;
  Double_t ptpc, dedx, bethe, deltadedx, dedx_sigma, tpcsignal;
  Double_t p, time, time_sigma, timezero, timezero_sigma, tof, tof_sigma, texp, texp_sigma, deltat, deltat_sigma, tofsignal;
  Double_t tpctofsignal;
  Double_t param[kNHistoParams];
  for (Int_t iev = startEv; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 100 == 0) printf("iev = %d\n", iev);
    /* check vertex */
    if (!analysisEvent->AcceptVertex()) continue;
    /* check collision candidate */
    if (!analysisEvent->IsCollisionCandidate()) continue;
    /* check centrality quality */
    if (analysisEvent->GetCentralityQuality() != 0.) continue;

    /*** ACCEPTED EVENT ***/

    /* apply time-zero TOF correction */
    analysisEvent->ApplyTimeZeroTOFCorrection();

    /* get centrality */
    param[kCentrality] = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check accepted track */
      if (!analysisTrack->AcceptTrack()) continue;
      /* get charge */
      charge = analysisTrack->GetSign() > 0. ? kPositive : kNegative;

      /*** ACCEPTED TRACK ***/
      
      /* check TOF pid */
      if (!analysisTrack->HasTOFPID() || !analysisTrack->HasTPCPID()) continue;
      
      /*** ACCEPTED TRACK WITH TPC+TOF PID ***/
      
      /* apply expected time correction */
      analysisTrack->ApplyTOFExpectedTimeCorrection();
      
      /* get track info */
      p = analysisTrack->GetP();

      /* get TPC info */
      dedx = analysisTrack->GetTPCdEdx();
      dedxN = analysisTrack->GetTPCdEdxN();
      ptpc = analysisTrack->GetTPCmomentum();
      param[kPtpc] = ptpc;
      param[kTPCNcls] = dedxN;
      
      /* get TOF info */
      time = analysisTrack->GetTOFTime();
      time_sigma = tofReso;
      timezero = analysisEvent->GetTimeZeroTOF(p);
      timezero_sigma = analysisEvent->GetTimeZeroTOFSigma(p);
      tof = time - timezero;
      tof_sigma = TMath::Sqrt(time_sigma * time_sigma + timezero_sigma * timezero_sigma);

      /* loop over particle IDs */
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
	
	/* check rapidity */
	if (TMath::Abs(analysisTrack->GetY(AliPID::ParticleMass(ipart))) > 0.5) continue;

	/*** ACCEPTED TRACK WITHIN CORRECT RAPIDITY ***/
	
	/* TOF expected time */
	texp = analysisTrack->GetTOFExpTime(ipart);
	texp_sigma = analysisTrack->GetTOFExpTimeSigma(ipart);
	
	/* TOF signal */
	deltat = tof - texp;
	deltat_sigma = TMath::Sqrt(tof_sigma * tof_sigma + texp_sigma * texp_sigma);
	tofsignal = deltat / deltat_sigma;
	param[kTOFsignal] = tofsignal;

	/* check TOF PID with 3-sigma cut */
	if (TMath::Abs(tofsignal) > 3.) continue;

	/*** ACCEPTED TRACK WITH COMPATIBLE TOF PID ***/

	/* TPC signal */
	bethe = tpcResponse->GetExpectedSignal(ptpc, ipart);
	deltadedx = dedx - bethe;
	param[kTPCdelta] = deltadedx;
	param[kTPCgain] = dedx / bethe;
	
	/* fill histo */
	hTPCcalib[ipart][charge]->Fill(param);

      } /* end of loop over particle IDs */      
    } /* end of loop over tracks */
  } /* end of loop over events */
  
  /* start stopwatch */
  timer.Stop();
  timer.Print();
  
  /* output */
  TFile *fileout = TFile::Open(Form("TPCcalib.%s", filename), "RECREATE");
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hTPCcalib[ipart][icharge]->Write();
    }
  }
  fileout->Close();
  
}

/**************************************************************/

TPCcalib_delta(const Char_t *filename, Int_t ipart, Int_t icharge, Float_t tofsigmaMin = -3., Float_t tofsigmaMax = 3.)
{

  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0.5, 0.5, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {0.5, 0., 1.5, 1.5, 2.0};

  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");

  /* open data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hTPCcalib = (THnSparseF *)filein->Get(Form("hTPCcalib_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));

  /* set momentum range */
  hTPCcalib->GetAxis(kPtpc)->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);

  /* set TOF pid range */
  //  hTPCcalib->GetAxis(kTOFsignal)->SetRangeUser(tofsigmaMin + 0.001, tofsigmaMax - 0.001);

  /* minimum-bias projection */
  hTPCcalib->GetAxis(kCentrality)->SetRange(1, NcentralityBins);
  TH2 *hDeltaPtpc = hTPCcalib->Projection(kTPCdelta, kPtpc);
  TObjArray *oaTPCdelta_mb = HistoUtils_FitPeak(NULL, hDeltaPtpc, 5, 3., 3., 100, Form("hTPCdelta_mb_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  delete hDeltaPtpc;

  /* convert to bethe-blok */
  TH1D *hTPCdelta_mb = (TH1D *)oaTPCdelta_mb->At(1);
  TH1D *hTPCbethe_mb = TPCcalib_delta2bethe(hTPCdelta_mb, ipart, Form("hTPCbethe_mb_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  
  /* centrality projections */
  TObjArray *oaTPCdelta_cent[NcentralityBins];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {

    hTPCcalib->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
    hDeltaPtpc = hTPCcalib->Projection(kTPCdelta, kPtpc);
    oaTPCdelta_cent[icent] = HistoUtils_FitPeak(NULL, hDeltaPtpc, 5, 3., 3., 100, Form("hTPCdelta_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    delete hDeltaPtpc;
  }

}

/**************************************************************/

TF1 *
TPCcalib_BetheBlochAleph(Double_t min, Double_t max)
{

  TF1 *f = new TF1("fBetheBlochAleph", "[0] * AliExternalTrackParam::BetheBlochAleph(x, [1], [2], [3], [4], [5])", min, max);
  f->SetParameter(0, 50.);
  f->SetParameter(1, 0.0283086);
  f->SetParameter(2, 2.63394e+01);
  f->SetParameter(3, 5.04114e-11);
  f->SetParameter(4, 2.12543);
  f->SetParameter(5, 4.88663);
  return f;

};

/**************************************************************/

void
TPCcalib_fitBetheBloch(TGraphErrors *g, TF1 *f, Int_t param = -1)
{

  if (param >= 0) {
    for (Int_t ipar = 0; ipar < 6; ipar++)
      if (ipar != param)
	f->FixParameter(ipar, f->GetParameter(ipar));
      else
	f->ReleaseParameter(ipar);
  }
  else {
    for (Int_t ipar = 0; ipar < 6; ipar++)
      f->ReleaseParameter(ipar);
  }
  g->Fit(f);

}

/**************************************************************/

TPCcalib_bethe(const Char_t *filename, Int_t icent = -1, Float_t tofsigmaMin = -3., Float_t tofsigmaMax = 3.)
{

  Double_t ptMin[AliPID::kSPECIES] = {0., 0., 0.5, 0.5, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {0., 0., 1.5, 1.5, 2.0};

  /* set centrality name */
  Char_t centName[1024];
  if (icent < 0 || icent >= NcentralityBins)
    sprintf(centName, "cent0090");
  else
    sprintf(centName, "cent%02d%02d", , centralityBin[icent], centralityBin[icent + 1]);
  
  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");

  /* open data */
  TFile *filein = TFile::Open(filename);

  /* output */
  TFile *fileout = TFile::Open(Form("TPCcalib_bethe_%s.%s", centName, filename), "RECREATE");
  
  /* loop over particles and charges */
  TH1D *hTPCbethe[AliPID::kSPECIES][kNCharges];
  TGraphErrors *gTPCbethe = new TGraphErrors();
  gTPCbethe->SetName("gTPCbethe");
  Int_t nPoints = 0;
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {

      /* check momentum range */
      if (ptMax[ipart] == 0. || ptMax[ipart] < ptMin[ipart]) continue;
      
      /* get histo */
      THnSparseF *hTPCcalib = (THnSparseF *)filein->Get(Form("hTPCcalib_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
      
      /* set centrality range */
      if (icent < 0 || icent >= NcentralityBins)
	hTPCcalib->GetAxis(kCentrality)->SetRange(1, NcentralityBins);
      else
	hTPCcalib->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
      
      /* set TOF pid range */
      //  hTPCcalib->GetAxis(kTOFsignal)->SetRangeUser(tofsigmaMin + 0.001, tofsigmaMax - 0.001);
      
      /* set momentum range */
      hTPCcalib->GetAxis(kPtpc)->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
      
      /* minimum-bias projection */
      TH2 *hDeltaPtpc = hTPCcalib->Projection(kTPCdelta, kPtpc);
      TObjArray *oaTPCdelta = HistoUtils_FitPeak(NULL, hDeltaPtpc, 5, 3., 3., 100, Form("hTPCdelta_%s_%s_%s", centName, AliPID::ParticleName(ipart), chargeName[icharge]));
      delete hDeltaPtpc;
      
      /* convert to bethe-blok */
      TH1D *hTPCdelta = (TH1D *)oaTPCdelta->At(1);
      hTPCbethe[ipart][icharge] = TPCcalib_delta2bethe(hTPCdelta, ipart, Form("hTPCbethe_%s_%s_%s", centName, AliPID::ParticleName(ipart), chargeName[icharge]));
      
      /* delete array */
      delete oaTPCdelta;

      /* add points to graph */
      for (Int_t ibin = 0; ibin < hTPCbethe[ipart][icharge]->GetNbinsX(); ibin++) {
	gTPCbethe->SetPoint(nPoints, hTPCbethe[ipart][icharge]->GetBinCenter(ibin + 1), hTPCbethe[ipart][icharge]->GetBinContent(ibin + 1));
	gTPCbethe->SetPointError(nPoints, 0.5 * hTPCbethe[ipart][icharge]->GetBinWidth(ibin + 1), hTPCbethe[ipart][icharge]->GetBinError(ibin + 1));
	nPoints++;
      }
      
      /* write */
      fileout->cd();
      hTPCbethe[ipart][icharge]->Write();
    }
  }

  TF1 *fBethe = TPCcalib_BetheBlochAleph(0.1, 1.e3);
  for (Int_t i = 0; i < 50; i++)
    for (Int_t ii = 0; ii < 6; ii++)
      TPCcalib_fitBetheBloch(gTPCbethe, fBethe, ii);
  TPCcalib_fitBetheBloch(gTPCbethe, fBethe, -1);

  /* write */
  fileout->cd();
  gTPCbethe->Write("gTPCbethe");
  fBethe->Write("fBetheBlochAleph");

  /* close file */
  fileout->Close();
  
  return;

}

/**************************************************************/

TPCcalib_gain(const Char_t *filename, Int_t ipart, Int_t icharge, Float_t tofsigmaMin = -3., Float_t tofsigmaMax = 3.)
{

  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");

  /* open data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hTPCcalib = (THnSparseF *)filein->Get(Form("hTPCcalib_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));

  /* set TOF pid range */
  //  hTPCcalib->GetAxis(kTOFsignal)->SetRangeUser(tofsigmaMin + 0.001, tofsigmaMax - 0.001);

  /* minimum-bias projection */
  hTPCcalib->GetAxis(kCentrality)->SetRange(1, NcentralityBins);
  TH2 *hGainPtpc = hTPCcalib->Projection(kTPCgain, kPtpc);
  TObjArray *oaTPCgain_mb = HistoUtils_FitPeak(NULL, hGainPtpc, 5, 3., 3., 100, Form("hTPCgain_mb_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  delete hGainPtpc;
  TH1D *hmb = (TH1D *)oaTPCgain_mb->At(1);
  
  /* centrality projections */
  TObjArray *oaTPCgain_cent[NcentralityBins];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {

    hTPCcalib->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
    hGainPtpc = hTPCcalib->Projection(kTPCgain, kPtpc);
    oaTPCgain_cent[icent] = HistoUtils_FitPeak(NULL, hGainPtpc, 5, 3., 3., 100, Form("hTPCgain_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    TH1D *hcent = (TH1D *)oaTPCgain_cent[icent]->At(1);
    hcent->Divide(hmb);
    delete hGainPtpc;
  }

}

/**************************************************************/

TH1D *
TPCcalib_delta2bethe(TH1D *hDelta, Int_t ipart, const Char_t *name = "hTPCbethe")
{

  Int_t nbins = hDelta->GetNbinsX();
  Float_t *xbins = new Float_t[nbins + 1];
  for(Int_t ibin = 0; ibin <= nbins; ibin++)
    xbins[ibin] = hDelta->GetBinLowEdge(ibin + 1) / AliPID::ParticleMass(ipart);
  
  TH1D *hBethe = new TH1D(name, "", nbins, xbins);
  AliTPCPIDResponse tpcResponse;
  Double_t ptpc, delta, bethe, deltae;
  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    ptpc = hDelta->GetBinCenter(ibin + 1);
    delta = hDelta->GetBinContent(ibin + 1);
    deltae = hDelta->GetBinError(ibin + 1);
    bethe = tpcResponse.GetExpectedSignal(ptpc, ipart);
    hBethe->SetBinContent(ibin + 1, bethe + delta);
    hBethe->SetBinError(ibin + 1, deltae);
  }

  return hBethe;
}
