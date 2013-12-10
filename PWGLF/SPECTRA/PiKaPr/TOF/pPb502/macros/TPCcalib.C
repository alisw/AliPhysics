#include "CommonDefs.C"

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

  /* open enabled flag map */
  TH1F *hEnabledFlag = NULL;
  if (enabledChannelsFileName) {
    TFile *enabledfile = TFile::Open(enabledChannelsFileName);
    hEnabledFlag = (TH1F *)enabledfile->Get("hEnabledFlag");
  }

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

  /* Start stopwatch */
  TStopwatch timer;
  timer.Start();

  if (1) {
    TFile *calibfile = TFile::Open("TZEROcalibration.root");
    hCentrality_TZEROA_mean = (TH1 *)calibfile->Get("hCentrality_TZEROA_mean");
    hCentrality_TZEROA_sigma = (TH1 *)calibfile->Get("hCentrality_TZEROA_sigma");
    hCentrality_TZEROC_mean = (TH1 *)calibfile->Get("hCentrality_TZEROC_calib");
    hCentrality_TZEROC_sigma = (TH1 *)calibfile->Get("hCentrality_TZEROC_sigma");
    hCentrality_TOF_mean = (TH1 *)calibfile->Get("hCentrality_TOF_mean");
    hCentrality_TOF_TZEROA_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROA_mean");
    hCentrality_TOF_TZEROC_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROC_mean");
    hCentrality_TOF_TZEROTOF_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROTOF_mean");

    TFile *resofile = TFile::Open("TZEROresolution.root");
    hCentrality_TZEROA_reso = (TH1 *)resofile->Get("hTZEROA_reso");
    hCentrality_TZEROC_reso = (TH1 *)resofile->Get("hTZEROC_reso");
    hCentrality_TZEROAND_reso = (TH1 *)resofile->Get("hTZEROAND_reso");
  }
  Double_t TZEROA_mean;
  Double_t TZEROA_sigma;
  Double_t TZEROC_mean;
  Double_t TZEROC_sigma;
  Double_t TOF_mean;
  Double_t TOF_TZEROA_mean;
  Double_t TOF_TZEROC_mean;
  Double_t TOF_TZEROTOF_mean;
  Double_t TZEROA;
  Double_t TZEROA_reso;
  Bool_t hasTZEROA;
  Double_t TZEROC;
  Double_t TZEROC_reso;
  Bool_t hasTZEROC;
  Double_t TZEROAND;
  Double_t TZEROAND_reso;
  Bool_t hasTZEROAND;
  Double_t TZEROTOF;
  Double_t TZEROTOF_reso;
  Bool_t hasTZEROTOF;
  Double_t TZEROMEAN;
  Double_t TZEROMEAN_weight;
  Double_t TZEROBEST;
  Double_t TZEROBEST_reso;


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
    if (iev % 10000 == 0) printf("iev = %d\n", iev);
    /* check event */
    if (!analysisEvent->AcceptEvent(acceptEventType)) continue;

    /*** ACCEPTED EVENT ***/

    /* get centrality */
    param[kCentrality] = analysisEvent->GetCentralityPercentile(centralityEstimator);
    cent = analysisEvent->GetCentralityPercentile(centralityEstimator);


    /* TZERO corrections */
    Int_t icent;
    for (icent = 0; icent < NcentralityBins; icent++)
      if (cent < centralityBin[icent + 1])
	break;
    TZEROA_mean = hCentrality_TZEROA_mean ? hCentrality_TZEROA_mean->GetBinContent(icent + 1) : 0.;
    TZEROA_sigma = hCentrality_TZEROA_sigma ? hCentrality_TZEROA_sigma->GetBinContent(icent + 1) : 1000.;
    TZEROC_mean  = hCentrality_TZEROC_mean ? hCentrality_TZEROC_mean->GetBinContent(icent + 1) : 0.;
    TZEROC_sigma = hCentrality_TZEROC_sigma ? hCentrality_TZEROC_sigma->GetBinContent(icent + 1) : 1000.;

    TOF_mean = hCentrality_TOF_mean ? hCentrality_TOF_mean->GetBinContent(icent + 1) : 0.;
    TOF_TZEROA_mean = hCentrality_TOF_TZEROA_mean ? hCentrality_TOF_TZEROA_mean->GetBinContent(icent + 1) : 0.;
    TOF_TZEROC_mean = hCentrality_TOF_TZEROC_mean ? hCentrality_TOF_TZEROC_mean->GetBinContent(icent + 1) : 0.;
    TOF_TZEROTOF_mean = hCentrality_TOF_TZEROTOF_mean ? hCentrality_TOF_TZEROTOF_mean->GetBinContent(icent + 1) : 0.;

    TZEROA_reso = hCentrality_TZEROA_reso ? hCentrality_TZEROA_reso->GetBinContent(icent + 1) : 70.;
    TZEROC_reso = hCentrality_TZEROC_reso ? hCentrality_TZEROC_reso->GetBinContent(icent + 1) : 70.;
    TZEROAND_reso = hCentrality_TZEROAND_reso ? hCentrality_TZEROAND_reso->GetBinContent(icent + 1) : 50.;

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
      if (!analysisTrack->HasTOFPID(hEnabledFlag) || !analysisTrack->HasTPCPID()) continue;
      

      /* get TOF info */
      index = analysisTrack->GetTOFIndex();
      time = analysisTrack->GetTOFTime() - TOF_mean;
      time_sigma = tofReso;
      /* TZEROTOF */
      TZEROTOF = analysisEvent->GetTimeZeroTOF(analysisTrack->GetP());
      TZEROTOF_reso = analysisEvent->GetTimeZeroTOFSigma(analysisTrack->GetP());
      hasTZEROTOF = TZEROTOF_reso < 150.;
      if (hasTZEROTOF) {
	TZEROTOF += TOF_TZEROTOF_mean - TOF_mean;
	TZEROTOF_reso *= TZEROTOF_resoScaleFactor;
      }
      /* TZEROA */
      TZEROA = analysisEvent->GetTimeZeroT0(1) - TZEROA_shift;
      //      TZEROA_reso = TZEROA_resolution[icent];
      hasTZEROA = TMath::Abs(TZEROA - TZEROA_mean) < 3. * TZEROA_sigma;
      TZEROA += TOF_TZEROA_mean - TOF_mean;
      /* TZEROC */
      TZEROC = analysisEvent->GetTimeZeroT0(2) - TZEROC_shift;
      //      TZEROC_reso = TZEROC_resolution[icent];
      hasTZEROC = TMath::Abs(TZEROC - TZEROC_mean) < 3. * TZEROC_sigma;
      TZEROC += TOF_TZEROC_mean - TOF_mean;
      /* TZEROAND */
      TZEROAND = (TZEROA + TZEROC) * 0.5;
      //      TZEROAND_reso = TZEROAND_resolution[icent];
      hasTZEROAND = hasTZEROA && hasTZEROC;
      /* TZEROMEAN */
      TZEROMEAN = TZEROTOF / TZEROTOF_reso / TZEROTOF_reso;
      TZEROMEAN_weight = 1. / TZEROTOF_reso / TZEROTOF_reso;
      if (hasTZEROAND) {
	//	printf("TZEROAND\n");
	TZEROMEAN += TZEROAND / TZEROAND_reso / TZEROAND_reso;
	TZEROMEAN_weight = 1. / TZEROAND_reso / TZEROAND_reso;
      }
      else if (hasTZEROA) {
	//	printf("TZEROA\n");
	TZEROMEAN += TZEROA / TZEROA_reso / TZEROA_reso;
	TZEROMEAN_weight = 1. / TZEROA_reso / TZEROA_reso;
      }
      else if (hasTZEROC) {
	//	printf("TZEROC\n");
	TZEROMEAN += TZEROC / TZEROC_reso / TZEROC_reso;
	TZEROMEAN_weight = 1. / TZEROC_reso / TZEROC_reso;
      }
      timezero = TZEROMEAN / TZEROMEAN_weight;
      timezero_sigma = TMath::Sqrt(1. / TZEROMEAN_weight);
      /* TZEROBEST */
      TZEROBEST = TZEROTOF;
      TZEROBEST_reso = TZEROTOF_reso;
      if (hasTZEROAND && TZEROAND_reso < TZEROBEST_reso) {
	TZEROBEST = TZEROAND;
	TZEROBEST_reso = TZEROAND_reso;
      }
      else if (hasTZEROA && TZEROA_reso < TZEROBEST_reso) {
	TZEROBEST = TZEROA;
	TZEROBEST_reso = TZEROA_reso;
      }
      if (hasTZEROC && TZEROC_reso < TZEROBEST_reso) {
	TZEROBEST = TZEROC;
	TZEROBEST_reso = TZEROC_reso;
      }
      timezero = TZEROBEST;
      timezero_sigma = TZEROBEST_reso;

      /* DEBUG */
      //      timezero = 0.;//TZEROTOF;
      //      timezero_sigma = 203.854691;//TZEROTOF_reso;

      //      if (timezero == 0.)
      //	printf("%f %f\n", timezero, timezero_sigma);

      timezero_sigma *= scaletimezerosigma;

      if (resetTZERO) {
	timezero = 0.;
	timezero_sigma = timezero_spread;
      }


      tof = time - timezero;
      tof_sigma = TMath::Sqrt(time_sigma * time_sigma + timezero_sigma * timezero_sigma);
      
      /* TOF expected time */
      texp = analysisTrack->GetTOFExpTime(iipart);
      texp_sigma = analysisTrack->GetTOFExpTimeSigma(iipart) * scaletexpreso[iipart];
      
      /* TOF signal */
      deltat = tof - texp;
      deltat_sigma = TMath::Sqrt(tof_sigma * tof_sigma + texp_sigma * texp_sigma);
      tofsignal = deltat / deltat_sigma;
      


      /*** ACCEPTED TRACK WITH TPC+TOF PID ***/
      
      /* get track info */
      p = analysisTrack->GetP();

      /* get TPC info */
      dedx = analysisTrack->GetTPCdEdx();
      dedxN = analysisTrack->GetTPCdEdxN();
      ptpc = analysisTrack->GetTPCmomentum();
      param[kPtpc] = ptpc;
      param[kTPCNcls] = dedxN;
      
      /* get TOF info */
      time = analysisTrack->GetTOFTime() - TOF_mean;
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

#if 0
   1  p0           4.94077e+01           nan           nan  -4.26724e+00
   2  p1           2.83086e-02           nan           nan  -7.44772e+03
   3  p2           2.62767e+01           nan           nan   3.89671e-01
   4  p3          -1.00486e-03           nan           nan   7.91984e+04
   5  p4           2.43863e+00           nan           nan  -1.64563e-04
   6  p5           4.37961e+00           nan           nan           inf
#endif


  return f;

};

/**************************************************************/

void
TPCcalib_fitBetheBloch(TGraphErrors *g, TF1 *f, Int_t param = -1)
{

  TVirtualFitter::SetMaxIterations(1000000);
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
  g->Fit(f, "", "", 0.5, 5.);

}

/**************************************************************/

TPCcalib_bethe(const Char_t *filename, Int_t icent = -1, Float_t tofsigmaMin = -3., Float_t tofsigmaMax = 3.)
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


  AliAnalysisTrack t;
  
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
      hTPCcalib->GetAxis(kTOFsignal)->SetRangeUser(tofsigmaMin + 0.001, tofsigmaMax - 0.001);
      
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
  /* TPC PID response */
  AliTPCPIDResponse *tpcResponse = AliAnalysisTrack::GetTPCResponse();
  Double_t ptpc, delta, bethe, deltae;
  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    ptpc = hDelta->GetBinCenter(ibin + 1);
    delta = hDelta->GetBinContent(ibin + 1);
    deltae = hDelta->GetBinError(ibin + 1);
    bethe = tpcResponse->GetExpectedSignal(ptpc, ipart);
    hBethe->SetBinContent(ibin + 1, bethe + delta);
    hBethe->SetBinError(ibin + 1, deltae);
  }

  return hBethe;
}
