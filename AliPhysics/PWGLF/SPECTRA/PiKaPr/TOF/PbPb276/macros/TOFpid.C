#define DISPLAY 0

/* SIGNAL SHAPE */
Bool_t GAUSSIAN_SIGNAL = kFALSE;
Bool_t GAUSSIANTAIL_SIGNAL = kFALSE;
Bool_t GAUSSIANTAIL2_SIGNAL = kFALSE;
Bool_t GAUSSIANPLUSGAUSSIANTAIL_SIGNAL = kFALSE;
Bool_t GAUSSIANPLUSEXPONENTIAL_SIGNAL = kFALSE;
Bool_t EXPECTED_SIGNAL_TAIL = kTRUE;
Bool_t EXPECTED_SIGNAL_TEMPLATE = kTRUE;
Bool_t EXPECTED_SIGNAL_FIT = kFALSE;
/* SIGNAL PARAMETERS */
Bool_t FIX_SIGNAL_MEAN = kTRUE;
Bool_t FIX_SIGNAL_SIGMA = kTRUE;
Bool_t FIX_SIGNAL_TAIL = kTRUE;
Float_t SCALE_SIGNAL_SIGMA = 1.;
Float_t SCALE_SIGNAL_TAIL = 1.;
/* OTHER STUFF */
Char_t *SIGNAL_PARAM_FILE = NULL;//"signalParamFile.root";
Float_t DEFAULT_SIGNAL_MEAN = 0.;
Float_t MIN_SIGNAL_MEAN = -0.2;
Float_t MAX_SIGNAL_MEAN = 0.2.;
Float_t DEFAULT_SIGNAL_SIGMA = 1.;
Float_t MIN_SIGNAL_SIGMA = 0.8;
Float_t MAX_SIGNAL_SIGMA = 1.2;
Float_t DEFAULT_SIGNAL_TAIL = 1.;
Float_t MIN_SIGNAL_TAIL = 0.5;
Float_t MAX_SIGNAL_TAIL = 1.5;
/* BACKGROUND */
Bool_t EXPECTED_BACKGROUND_TAIL = kTRUE;
Bool_t EXPECTED_BACKGROUND_TEMPLATE = kTRUE;
Bool_t EXPECTED_BACKGROUND_FIT = kFALSE;
Bool_t GAUSSIAN_BACKGROUND = kFALSE;
Bool_t USE_ELECTRON_BACKGROUND = kTRUE;
/* BACKGROUND PARAMETERS */
Bool_t FIX_BACKGROUND_MEAN = kTRUE;
Bool_t FIX_BACKGROUND_SIGMA = kTRUE;
Bool_t FIX_BACKGROUND_TAIL = kTRUE;
Float_t SCALE_BACKGROUND_SIGMA = 1.;
Float_t SCALE_BACKGROUND_TAIL = 1.;
/* MISMATCH */
Bool_t NO_MISMATCH = kFALSE;
Bool_t EXPECTED_MISMATCH = kTRUE;
Bool_t EXPONENTIAL_MISMATCH = kFALSE;
Bool_t UNIFORM_MISMATCH = kFALSE;
Bool_t DOUBLEEXPONENTIAL_MISMATCH = kFALSE;
Bool_t UNIFORMPLUSEXPONENTIAL_MISMATCH = kFALSE;
/* AUTOMATIC CONSTRAINS */
Float_t FIT_ELECTRONS_PT_MIN = 0.;
Float_t FIT_ELECTRONS_PT_MAX = 0.;
Float_t FIT_MUONS_PT_MIN = 0.;
Float_t FIT_MUONS_PT_MAX = 0.;
Float_t FIT_PIONS_PT_MIN = 0.;
Float_t FIT_PIONS_PT_MAX = 5.;

Float_t CONSTRAINSIGNAL_LIMIT = 0.;
Float_t CONSTRAINTAIL_LIMIT = 0.;

enum EFitParams_t {
  kMean,
  kSigma,
  kTail,
  kTotalCounts,
  kIntegralCounts,
  kSignalCounts,
  kBkg1Counts,
  kBkg2Counts,
  kBkg3Counts,
  kBkg4Counts,
  kMismatchCounts,
  kNFitParams
};

/* fit output params name */
const Char_t *fitParamName[kNFitParams] = {
  "Mean",
  "Sigma",
  "Tail",
  "TotalCounts",
  "IntegralCounts",
  "SignalCounts",
  "Bkg1Counts",
  "Bkg2Counts",
  "Bkg3Counts",
  "Bkg4Counts",
  "MismatchCounts"
};

/* fit output params title */
const Char_t *fitParamTitle[kNFitParams] = {
  "Signal mean;p_{T} (GeV/c);#mu (ps)",
  "Signal sigma;p_{T} (GeV/c);#sigma (ps)",
  "Signal tail;p_{T} (GeV/c);#sigma_{tail} (ps)",
  "Total counts;p_{T} (GeV/c);counts",
  "Total counts within fit range;p_{T} (GeV/c);counts",
  "Signal counts;p_{T} (GeV/c);counts",
  "Background-1 counts;p_{T} (GeV/c);counts",
  "Background-2 counts;p_{T} (GeV/c);counts",
  "Background-3 counts;p_{T} (GeV/c);counts",
  "Background-4 counts;p_{T} (GeV/c);counts",
  "Mismatch counts within fit range;p_{T} (GeV/c);counts"
};

/**************************************************************/
/*** GENERATION OF TEMPLATE HISTOS ****************************/
/**************************************************************/

const Int_t NmismatchTrials = 1;
const Int_t NexpectedTrials = 1;

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
  kTOFsigma,
  kPt,
  kTPCsigma,
  kNHistoParams
};
/**************************************************************/
const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
/**************************************************************/
const Int_t NtofsigmaBins = 1750;
Double_t tofsigmaBin[NtofsigmaBins + 1];
Double_t tofsigmaMin = -100., tofsigmaMax = 250., tofsigmaStep = (tofsigmaMax - tofsigmaMin) / NtofsigmaBins;
/**************************************************************/
const Int_t NtofsignalBins = 2000;
Double_t tofsignalBin[NtofsignalBins + 1];
Double_t tofsignalMin = -48800., tofsignalMax = 48800., tofsignalStep = (tofsignalMax - tofsignalMin) / NtofsignalBins;
/**************************************************************/
const Int_t NptBins = 46;
Double_t ptBin[NptBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
/**************************************************************/
const Int_t NmtBins = 46;
Double_t mtBin[NmtBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
/**************************************************************/
const Int_t NtpcsigmaBins = 10;
Double_t tpcsigmaBin[NtpcsigmaBins + 1];
Double_t tpcsigmaMin = -5., tpcsigmaMax = 5., tpcsigmaStep = (tpcsigmaMax - tpcsigmaMin) / NtpcsigmaBins;
/**************************************************************/
Int_t NparamsBins[kNHistoParams] = {NcentralityBins, NtofsigmaBins, NptBins, NtpcsigmaBins};
Double_t *paramBin[kNHistoParams] = {centralityBin, tofsigmaBin, ptBin, tpcsigmaBin};
/**************************************************************/

/**************************************************************/
/**************************************************************/
/**************************************************************/

Double_t tofReso = 85.;
Double_t tofTail = 80.;
const Char_t *t0FillOnlineFileName = "~/ALICE.2011/ANALYSIS/TOFSpectraPbPb/t0fill/T0FillOnline.139465.extended.root";
Double_t t0Fill_offset = -1.26416e+04;
//const Char_t *t0FillOnlineFileName = "T0FillOnline.117116.extended.root";
//Double_t t0Fill_offset = -1.35543e+04;
const Char_t *enabledChannelsFileName = "~/ALICE.2011/ANALYSIS/TOFSpectraPbPb/enabledch/enabledChannels.139507.root";

/**************************************************************/

AliTOFGeometry tofGeo;
Float_t c = TMath::C() * 1.e2 / 1.e12; /* cm/ps */
Float_t c_1 = 1. / c;

Double_t
GenerateRandomHit(TH1F *hT0Fill, Double_t t0fill, Int_t index)
{
  
  Int_t det[5];
  Float_t length, timeexp, pos[3];
  
  /* compute length and expected time */
  tofGeo.GetVolumeIndices(index, det);
  tofGeo.GetPosPar(det, pos);
  length = 0.;
  for (Int_t i = 0; i < 3; i++) length += pos[i] * pos[i];
  length = TMath::Sqrt(length);
  timeexp = length * c_1;
  
  Double_t hittime = hT0Fill->GetRandom() - t0fill + timeexp;
  return hittime;

}

/**************************************************************/

TOFpid(const Char_t *filename, Int_t ipart, Int_t icharge, Int_t iipart, Bool_t electronCut = kFALSE, Bool_t cutOnTPC = kFALSE, Float_t tpcsignalMin = -2., Float_t tpcsignalMax = 2., Int_t evMax = kMaxInt, Int_t startEv = 0, Bool_t mcFlag = kFALSE)
{
  
  printf("****************************************\n");
  printf("RUNNING TOF PID:\n");
  printf("RAPIDITY-CUT:   %s\n", AliPID::ParticleName(ipart));
  printf("CHARGE:         %s\n", chargeName[icharge]);
  printf("PARTICLE-ID:    %s\n", AliPID::ParticleName(iipart));
  if (electronCut)
    printf("-> ELECTRON CUT REQUESTED\n");
  if (cutOnTPC)
    printf(" -> CUT-ON-TPC [%3.1f,%3.1f]\n", tpcsignalMin, tpcsignalMax);
  printf("****************************************\n");

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

  /* create TOF response with tail */
  gROOT->LoadMacro("~/ALICE.2011/ANALYSIS/TOFSpectraPbPb/macros/TOFsignal.C");
  TF1 *fTOFsignal = new TF1("fTOFsignal", TOFsignal, -2440., 2440., 4);
  fTOFsignal->SetParameter(0, 1.);
  fTOFsignal->SetParameter(1, 0.);
  fTOFsignal->SetParameter(2, tofReso);
  fTOFsignal->SetParameter(3, tofTail);
  
  /* open file, get tree and connect */
  TFile *filein = TFile::Open(filename);
  TTree *treein = (TTree *)filein->Get("aodTree");
  printf("got \"aodTree\": %d entries\n", treein->GetEntries());
  AliAnalysisEvent *analysisEvent = new AliAnalysisEvent();
  TClonesArray *analysisTrackArray = new TClonesArray("AliAnalysisTrack");
  AliAnalysisTrack *analysisTrack = NULL;
  treein->SetBranchAddress("AnalysisEvent", &analysisEvent);
  treein->SetBranchAddress("AnalysisTrack", &analysisTrackArray);

  /* open hT0fill for mismatch */
  TFile *filein_T0Fill = TFile::Open(t0FillOnlineFileName);
  TH1F *hT0Fill = (TH1F *)filein_T0Fill->Get("hT0Fill");
  Double_t t0fill = t0Fill_offset;

  /* open enabled flag map */
  TFile *enabledfile = TFile::Open(enabledChannelsFileName);
  TH1F *hEnabledFlag = (TH1F *)enabledfile->Get("hEnabledFlag");

  /**************************************************************/
  /*** HISTOS ***************************************************/
  /**************************************************************/

  /* run-time binning */
  for (Int_t ibin = 0; ibin < NtofsigmaBins + 1; ibin++)
    tofsigmaBin[ibin] = tofsigmaMin + ibin * tofsigmaStep;
  for (Int_t ibin = 0; ibin < NtofsignalBins + 1; ibin++)
    tofsignalBin[ibin] = tofsignalMin + ibin * tofsignalStep;
  for (Int_t ibin = 0; ibin < NtpcsigmaBins + 1; ibin++)
    tpcsigmaBin[ibin] = tpcsigmaMin + ibin * tpcsigmaStep;

  /* histos */
  TH1F *hAcceptedEvents = new TH1F(Form("hAcceptedEvents_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin);
  TH2I *hAcceptedTracks = new TH2I(Form("hAcceptedTracks_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin, NptBins, ptBin);

  TH3I *hTOFpid = new TH3I(Form("hTOFpid_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsigmaBins, tofsigmaBin);
  TH3I *hTOFmismatch = new TH3I(Form("hTOFmismatch_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsigmaBins, tofsigmaBin);
  TH3I *hTOFexpected[AliPID::kSPECIES];
  for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
    hTOFexpected[iiipart] = new TH3I(Form("hTOFexpected_%s_%s_%sID_%sBKG", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart), AliPID::ParticleName(iiipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsigmaBins, tofsigmaBin);
  }
  
  TH3I *hTOFpid_delta = new TH3I(Form("hTOFpid_delta_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsignalBins, tofsignalBin);
  TH3I *hTOFmismatch_delta = new TH3I(Form("hTOFmismatch_delta_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsignalBins, tofsignalBin);
  TH3I *hTOFexpected_delta[AliPID::kSPECIES];
  for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
    hTOFexpected_delta[iiipart] = new TH3I(Form("hTOFexpected_delta_%s_%s_%sID_%sBKG", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart), AliPID::ParticleName(iiipart)), "", NcentralityBins, centralityBin, NptBins, ptBin, NtofsignalBins, tofsignalBin);
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

  /* verbose */
  printf("***** RUNNING for %s %s *****\n", chargeName[icharge], AliPID::ParticleName(ipart));
  if (cutOnTPC) {
    printf("***** CUT-ON-TPC requested %3.1f-%3.1f *****\n", tpcsignalMin, tpcsignalMax);
  }

  /* loop over events */
  Bool_t hastofpid;
  Int_t charge, index;
  UShort_t dedxN;
  Double_t cent, p, pt, mt, tofsignal, tpcsignal;
  Double_t dedx, bethe, deltadedx, dedx_sigma, ptpc;
  Double_t time, time_sigma, timezero, timezero_sigma, tof, tof_sigma, texp, texp_sigma, deltat, deltat_sigma, tof_rnd, tof_th, signal_smear, timezero_smear, texp_smear;

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
    cent = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);

    /* fill histos */
    hAcceptedEvents->Fill(cent);

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check accepted track */
      if (!analysisTrack->AcceptTrack()) continue;
      /* check rapidity */
      if (TMath::Abs(analysisTrack->GetY(AliPID::ParticleMass(ipart))) > 0.5) continue;
      /* check charge */
      charge = analysisTrack->GetSign() > 0. ? kPositive : kNegative;
      if (charge != icharge) continue;

      /*** ACCEPTED TRACK ***/

      /* get track info */
      p = analysisTrack->GetP();
      pt = analysisTrack->GetPt();
      
      /* compute track mt */
      mt = TMath::Sqrt(pt * pt + AliPID::ParticleMass(ipart) * AliPID::ParticleMass(ipart));
	
      /* get TPC info */
      dedx = analysisTrack->GetTPCdEdx();
      dedxN = analysisTrack->GetTPCdEdxN();
      ptpc = analysisTrack->GetTPCmomentum();
      
      /* TPC signal */
      bethe = tpcResponse->GetExpectedSignal(ptpc, iipart);
      /* fix electron signal */
      if (iipart == AliPID::kElectron)
	bethe += 23.;
      deltadedx = dedx - bethe;
      dedx_sigma = 0.07 * bethe;
      tpcsignal = deltadedx / dedx_sigma;

      /* electronCut requested, remove electrons */
      if (electronCut) {
	/* TPC signal */
	bethe = tpcResponse->GetExpectedSignal(ptpc, AliPID::kElectron);
	/* fix electron signal */
	bethe += 23.;
	deltadedx = dedx - bethe;
	dedx_sigma = 0.07 * bethe;
	tpcsignal = deltadedx / dedx_sigma;
	if (TMath::Abs(tpcsignal) < 1.5) continue;
      }
      
      /* cut on TPC signal if requested */
      if (cutOnTPC && (tpcsignal < tpcsignalMin || tpcsignal > tpcsignalMax))
	continue;
      
      /* fill histos */
      hAcceptedTracks->Fill(cent, pt);
      
      /* set TOF pid flag */
      hastofpid = analysisTrack->HasTOFPID();
      /* check channel enabled */
      index = analysisTrack->GetTOFIndex();
      //      if (hEnabledFlag->GetBinContent(index + 1) == 0.) hastofpid = kFALSE;
      
      /* check TOF pid */
      if (!hastofpid)
	continue;
      
      /*** ACCEPTED TRACK WITH TOF PID ***/
      
      /* apply expected time correction */
      analysisTrack->ApplyTOFExpectedTimeCorrection();
      
      /* get TOF info */
      time = analysisTrack->GetTOFTime();
      time_sigma = tofReso;
      timezero = analysisEvent->GetTimeZeroTOF(p);
      timezero_sigma = analysisEvent->GetTimeZeroTOFSigma(p);
      tof = time - timezero;
      tof_sigma = TMath::Sqrt(time_sigma * time_sigma + timezero_sigma * timezero_sigma);
      
      /* TOF expected time */
      texp = analysisTrack->GetTOFExpTime(iipart);
      texp_sigma = analysisTrack->GetTOFExpTimeSigma(iipart);
      
      /* TOF signal */
      deltat = tof - texp;
      deltat_sigma = TMath::Sqrt(tof_sigma * tof_sigma + texp_sigma * texp_sigma);
      tofsignal = deltat / deltat_sigma;
      
      /* fill histo */
      hTOFpid->Fill(cent, pt, tofsignal);
      hTOFpid_delta->Fill(cent, p, deltat);

      /*** EXPECTED MISMATCH ***/
      
      /* loop to generated random hits */
      for (Int_t irnd = 0; irnd < NmismatchTrials; irnd++) {
	
	/* generate ramdom tof values */
	tof_rnd = GenerateRandomHit(hT0Fill, t0fill, index);

	/* TOF signal */
	deltat = tof_rnd - texp;
	tofsignal = deltat / deltat_sigma;
	
	/* fill histo */
	hTOFmismatch->Fill(cent, pt, tofsignal);
	hTOFmismatch_delta->Fill(cent, p, deltat);
	
      } /* end of loop over generated random hits */
	
      /*** EXPECTED SIGNALS ***/
      
      /* loop over other particles */
      for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
	
	/* generate expected tof value */
	tof_th = analysisTrack->GetTOFExpTime(iiipart);
	texp_sigma = analysisTrack->GetTOFExpTimeSigma(iiipart);
	 	
	/* loop over many trials */
	for (Int_t irnd = 0; irnd < NexpectedTrials; irnd++) {
	  
	  /* tof response smearing */
	  signal_smear = fTOFsignal->GetRandom();
	  /* timezero resolution smearing */
	  timezero_smear = gRandom->Gaus(0., timezero_sigma);
	  /* texp resolution smearing */
	  texp_smear = gRandom->Gaus(0., texp_sigma);
	  
	  /* deltat and sigma */
	  deltat = tof_th - texp + signal_smear + timezero_smear + texp_smear;
	  tofsignal = deltat / deltat_sigma;
	  
	  /* fill histo */
	  hTOFexpected[iiipart]->Fill(cent, pt, tofsignal);
	  hTOFexpected_delta[iiipart]->Fill(cent, p, deltat);
	    
	} /* end of loop over many trials */
      } /* end of loop over other particle */
    } /* end of loop over tracks */
  } /* end of loop over events */
  
  /* stop stopwatch */
  timer.Stop();
  timer.Print();
  
  /* output */
  TString outputstring = "TOFpid";
  if (electronCut)
    outputstring += "_electronCut";
  if (cutOnTPC)
    outputstring += Form("_cutOnTPC[%3.1f,%3.1f]", , tpcsignalMin, tpcsignalMax);
  outputstring += Form("_%s_%s_%sID.%s", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart), filename);
  TFile *fileout = TFile::Open(outputstring.Data(), "RECREATE");
  hAcceptedEvents->Write();
  hAcceptedTracks->Write();
  hTOFpid->Write();
  hTOFmismatch->Write();
  for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
    hTOFexpected[iiipart]->Write();
  hTOFpid_delta->Write();
  hTOFmismatch_delta->Write();
  for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
    hTOFexpected_delta[iiipart]->Write();
  
  fileout->Close();
  
}

//___________________________________________________________________________________

TOFspectra_defaultFit(const Char_t *filename)
{

Bool_t EXPECTED_SIGNAL_TEMPLATE = kTRUE;
Bool_t EXPECTED_SIGNAL_FIT = kFALSE;
Bool_t EXPECTED_BACKGROUND_TEMPLATE = kFALSE;
Bool_t EXPECTED_BACKGROUND_FIT = kTRUE;

}

//___________________________________________________________________________________

TOFspectra_defaultFit_fitElectrons(const Char_t *filename, Float_t electronLimit = 5.)
{


  TOFspectra(filename, electronLimit);
}

//___________________________________________________________________________________

TOFspectra_signalFit(Bool_t fixParams = kTRUE, Float_t scaleSigma = 1., Float_t scaleTail = 1.)
{

  EXPECTED_SIGNAL_TEMPLATE = kFALSE;
  EXPECTED_SIGNAL_FIT = kTRUE;
  FIX_SIGNAL_MEAN = fixParams;
  FIX_SIGNAL_SIGMA = fixParams;
  FIX_SIGNAL_TAIL = fixParams;
  SCALE_SIGNAL_SIGMA = scaleSigma;
  SCALE_SIGNAL_TAIL = scaleTail;
  
}

//___________________________________________________________________________________

TOFspectra_bkgFit(Bool_t fixParams = kTRUE, Float_t scaleSigma = 1., Float_t scaleTail = 1.)
{

  EXPECTED_BACKGROUND_TEMPLATE = kFALSE;
  EXPECTED_BACKGROUND_FIT = kTRUE;
  FIX_BACKGROUND_MEAN = fixParams;
  FIX_BACKGROUND_SIGMA = fixParams;
  FIX_BACKGROUND_TAIL = fixParams;
  SCALE_BACKGROUND_SIGMA = scaleSigma;
  SCALE_BACKGROUND_TAIL = scaleTail;
  
  TOFspectra(filename);

}

//___________________________________________________________________________________

void
TOFspectra(const Char_t *filename, Float_t electronLimit = 0.)
{

  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
	TOFspectrum(filename, ipart, icharge, ipart, icent, -1., -1., electronLimit);
  
}

//___________________________________________________________________________________

/* fit ranges */
Double_t fitPtMin[AliPID::kSPECIES] = {0.5, 0.5, 0.3, 0.4, 0.5};
Double_t fitPtMax[AliPID::kSPECIES] = {3.0, 3.0, 3.0, 3.0, 5.0};
Double_t fitSigmaMin[AliPID::kSPECIES] = {tofsigmaMin, tofsigmaMin, -25., -75., -65.};
Double_t fitSigmaMax[AliPID::kSPECIES] = {tofsigmaMax, tofsigmaMax, 225., 200., 100.};

void
TOFspectrum(const Char_t *filename, Int_t ipart, Int_t icharge, Int_t iipart, Int_t icent, Float_t ptMin = -1., Float_t ptMax = -1., Bool_t checkHistoFlag = kFALSE)
{

  printf("****************************************\n");
  printf("RUNNING SPECTRA FIT:\n");
  printf("RAPIDITY-CUT:   %s\n", AliPID::ParticleName(ipart));
  printf("CHARGE:         %s\n", chargeName[icharge]);
  printf("PARTICLE:       %s\n", AliPID::ParticleName(iipart));
  printf("CENTRALITY BIN: %d\n", icent);
  printf("****************************************\n");

  /* open data */
  TFile *filein = TFile::Open(filename);

  /* get number of events */
  TH1F *hAcceptedEvents = (TH1F *)filein->Get(Form("hAcceptedEvents_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
  if (!hAcceptedEvents) {
    printf("cannot find %s\n", Form("hAcceptedEvents_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
    return;
  }
  Double_t nevents;
  if (icent < 0 || icent >= NcentralityBins)
    nevents = hAcceptedEvents->Integral(1, NcentralityBins);
  else
    nevents = hAcceptedEvents->Integral(icent + 1, icent + 1);
  printf("N_EVENTS      : %d\n", nevents);
  printf("****************************************\n");
  
  /* get histos */
  TH2I *hAcceptedTracks = (TH2I *)filein->Get(Form("hAcceptedTracks_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
  if (!hAcceptedTracks) {
    printf("cannot find %s\n", Form("hAcceptedTracks_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
    //    return;
  }
  TH3I *hTOFpid = (TH3I *)filein->Get(Form("hTOFpid_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
  if (!hTOFpid) {
    printf("cannot find %s\n", Form("hTOFpid_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
    return;
  }
  TH3I *hTOFmismatch = (TH3I *)filein->Get(Form("hTOFmismatch_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
  if (!hTOFmismatch) {
    printf("cannot find %s\n", Form("hTOFpid_%s_%s_%sID", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)));
    return;
  }
  TH3I *hTOFexpected[AliPID::kSPECIES];
  for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
    hTOFexpected[iiipart] = (TH3I *)filein->Get(Form("hTOFexpected_%s_%s_%sID_%sBKG", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart), AliPID::ParticleName(iiipart)));
    if (!hTOFexpected[iiipart]) {
      printf("cannot find %s\n", Form("hTOFexpected_%s_%s_%sID_%sBKG", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart), AliPID::ParticleName(iiipart)));
      return;
    }
  }

  /* setup centrality range */
  if (icent < 0 || icent >= NcentralityBins) {
    printf("WARNING: undefined centrality -> using 00-90\% range\n");
    if (hAcceptedTracks) hAcceptedTracks->GetXaxis()->SetRange(1, NcentralityBins);
    hTOFpid->GetXaxis()->SetRange(1, NcentralityBins);
    hTOFmismatch->GetXaxis()->SetRange(1, NcentralityBins);
    for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
      hTOFexpected[iiipart]->GetXaxis()->SetRange(1, NcentralityBins);
  }
  else {
    printf("***** FITTING CENTRALITY-BIN [%02d, %02d] %% *****\n", centralityBin[icent], centralityBin[icent + 1]);
    if (hAcceptedTracks) hAcceptedTracks->GetXaxis()->SetRange(icent + 1, icent + 1);
    hTOFpid->GetXaxis()->SetRange(icent + 1, icent + 1);
    hTOFmismatch->GetXaxis()->SetRange(icent + 1, icent + 1);
    for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
      hTOFexpected[iiipart]->GetXaxis()->SetRange(icent + 1, icent + 1);
  }

  /* init flags */
  Bool_t requestedRange = kFALSE;
  Bool_t fitElectrons = kTRUE;
  Bool_t fitMuons = kTRUE;
  Bool_t fitPions = kTRUE;

  /* setup pt range if requested */
  if (ptMin > -0.001 && ptMax > -0.001 && ptMax > ptMin) {
    printf("***** FITTING PT-BIN [%f, %f] GeV/c *****\n", ptMin, ptMax);
    requestedRange = kTRUE;

    /* check electron-fit is allowed */
    fitElectrons = kTRUE;
    if ((ptMin + 0.001) < FIT_ELECTRONS_PT_MIN || (ptMax - 0.001) > FIT_ELECTRONS_PT_MAX)
      fitElectrons = kFALSE;
    /* check muon-fit is allowed */
    fitMuons = kTRUE;
    if ((ptMin + 0.001) < FIT_MUONS_PT_MIN || (ptMax - 0.001) > FIT_MUONS_PT_MAX)
      fitMuons = kFALSE;
    /* check pion-fit is allowed */
    fitPions = kTRUE;
    if ((ptMin + 0.001) < FIT_PIONS_PT_MIN || (ptMax - 0.001) > FIT_PIONS_PT_MAX)
      fitPions = kFALSE;

    
    hTOFpid->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hTOFmismatch->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
      hTOFexpected[iiipart]->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
  }

  /* output */
  Char_t outfilename[1024];
  if (icent < 0 || icent >= NcentralityBins)
    sprintf(outfilename, "TOFspectrum_cent0090_%s_%s_%sID.root", AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart));
  else {
    sprintf(outfilename, "TOFspectrum_cent%02d%02d_%s_%s_%sID.root", centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart));
  }
  TFile *fileout = TFile::Open(outfilename, "RECREATE");
  TDirectory *fitDir = fileout->mkdir("FitParams");
  /* canvas */
  TCanvas *canvas = new TCanvas("canvas");
  canvas->SetLogy();
  /* histo */
  TH1D *hFitParamHisto[kNFitParams];
  for (Int_t iparam = 0; iparam < kNFitParams; iparam++)
    hFitParamHisto[iparam] = new TH1D(Form("h%s", fitParamName[iparam]), fitParamTitle[iparam], NptBins, ptBin);

  /* loop over ptBins */
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    
    if (!requestedRange) {
      if ((ptBin[ipt] + 0.001) < fitPtMin[ipart] || (ptBin[ipt + 1] - 0.001) > fitPtMax[ipart]) continue;
      printf("***** FITTING PT-BIN [%f, %f] GeV/c *****\n", ptBin[ipt], ptBin[ipt + 1]);

      /* check electron-fit is allowed */
      fitElectrons = kTRUE;
      if ((ptBin[ipt] + 0.001) < FIT_ELECTRONS_PT_MIN || (ptBin[ipt + 1] - 0.001) > FIT_ELECTRONS_PT_MAX)
	fitElectrons = kFALSE;
      /* check muon-fit is allowed */
      fitMuons = kTRUE;
      if ((ptBin[ipt] + 0.001) < FIT_MUONS_PT_MIN || (ptBin[ipt + 1] - 0.001) > FIT_MUONS_PT_MAX)
	fitMuons = kFALSE;
      /* check pion-fit is allowed */
      fitPions = kTRUE;
      if ((ptBin[ipt] + 0.001) < FIT_PIONS_PT_MIN || (ptBin[ipt + 1] - 0.001) > FIT_PIONS_PT_MAX)
	fitPions = kFALSE;

      hTOFpid->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hTOFmismatch->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++)
	hTOFexpected[iiipart]->GetYaxis()->SetRange(ipt + 1, ipt + 1);
    }
    
    /* nsigma projections */
    TH1 *hSignal_py = hTOFpid->Project3D("z");
    TH1 *hMismatch_py = hTOFmismatch->Project3D("z");
    TH1 *hSignalExp_py[AliPID::kSPECIES];
    TH1 *hSignalExpTail_py[AliPID::kSPECIES];
    for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
      hSignalExp_py[iiipart] = hTOFexpected[iiipart]->Project3D("z");
      hSignalExpTail_py[iiipart] = hTOFexpected[iiipart]->Project3D("z");
    }

    /* prepare histos for the fitter */
    Int_t partbkg1[AliPID::kSPECIES] = {AliPID::kKaon, AliPID::kKaon, AliPID::kKaon, AliPID::kPion, AliPID::kPion};
    Int_t partbkg2[AliPID::kSPECIES] = {AliPID::kProton, AliPID::kProton, AliPID::kProton, AliPID::kProton, AliPID::kKaon};
    Int_t partbkg3[AliPID::kSPECIES] = {AliPID::kPion, AliPID::kPion, AliPID::kElectron, AliPID::kElectron, AliPID::kElectron};
    Int_t partbkg4[AliPID::kSPECIES] = {AliPID::kMuon, AliPID::kElectron, AliPID::kMuon, AliPID::kMuon, AliPID::kMuon};
    TH1 *hSigExp_py, *hBkgExp1_py, *hBkgExp2_py, *hBkgExp3_py, *hBkgExp4_py;
    hSigExp_py = EXPECTED_SIGNAL_TAIL ? hSignalExpTail_py[iipart] : hSignalExp_py[iipart];
    hBkgExp1_py = EXPECTED_BACKGROUND_TAIL ? hSignalExpTail_py[partbkg1[iipart]] : hSignalExp_py[partbkg1[iipart]];
    hBkgExp2_py = EXPECTED_BACKGROUND_TAIL ? hSignalExpTail_py[partbkg2[iipart]] : hSignalExp_py[partbkg2[iipart]];
    hBkgExp3_py = EXPECTED_BACKGROUND_TAIL ? hSignalExpTail_py[partbkg3[iipart]] : hSignalExp_py[partbkg3[iipart]];
    hBkgExp4_py = EXPECTED_BACKGROUND_TAIL ? hSignalExpTail_py[partbkg4[iipart]] : hSignalExp_py[partbkg4[iipart]];

    /* check histos if requested */
    if (checkHistoFlag) {
      TCanvas *cCheckHisto = new TCanvas("cCheckHisto");
      cCheckHisto->Divide(2, 3);
      cCheckHisto->cd(1);
      hSignal_py->Draw();
      cCheckHisto->cd(2);
      hSigExp_py->Draw();
      cCheckHisto->cd(3);
      hBkgExp1_py->Draw();
      cCheckHisto->cd(4);
      hBkgExp2_py->Draw();
      cCheckHisto->cd(5);
      hBkgExp3_py->Draw();
      cCheckHisto->cd(6);
      hMismatch_py->Draw();
      return;
    }

    Double_t rangeMin = fitSigmaMin[iipart], rangeMax = fitSigmaMax[iipart];
    Bool_t constrainSignal = kFALSE;
    Bool_t constrainBkg1 = kFALSE;
    Bool_t constrainBkg2 = kFALSE;
    Bool_t forceGaussianSignal = kFALSE;
    Bool_t fitBkg1 = kTRUE, fitBkg2 = kTRUE, fitBkg3 = kTRUE, fitBkg4 = kTRUE;

    /* check whether we can fit electrons */
    if (!fitElectrons) {
      printf("INHIBIT FIT ELECTRONS\n");
      if (partbkg1[iipart] == AliPID::kElectron) fitBkg1 = kFALSE;
      if (partbkg2[iipart] == AliPID::kElectron) fitBkg2 = kFALSE;
      if (partbkg3[iipart] == AliPID::kElectron) fitBkg3 = kFALSE;
      if (partbkg4[iipart] == AliPID::kElectron) fitBkg4 = kFALSE;
    }
    /* check whether we can fit muons */
    if (!fitMuons) {
      printf("INHIBIT FIT MUONS\n");
      if (partbkg1[iipart] == AliPID::kMuon) fitBkg1 = kFALSE;
      if (partbkg2[iipart] == AliPID::kMuon) fitBkg2 = kFALSE;
      if (partbkg3[iipart] == AliPID::kMuon) fitBkg3 = kFALSE;
      if (partbkg4[iipart] == AliPID::kMuon) fitBkg4 = kFALSE;
    }
    /* check whether we can fit pions */
    if (!fitPions) {
      printf("INHIBIT FIT PIONS\n");
      if (partbkg1[iipart] == AliPID::kPion) fitBkg1 = kFALSE;
      if (partbkg2[iipart] == AliPID::kPion) fitBkg2 = kFALSE;
      if (partbkg3[iipart] == AliPID::kPion) fitBkg3 = kFALSE;
      if (partbkg4[iipart] == AliPID::kPion) fitBkg4 = kFALSE;
    }


    /* fit */
    Double_t param[kNFitParams];
    Double_t param_err[kNFitParams];
    TOFpid_fit(hSignal_py, hSigExp_py, hBkgExp1_py, hBkgExp2_py, hBkgExp3_py, hBkgExp4_py, hMismatch_py, rangeMin, rangeMax, fitBkg1, fitBkg2, fitBkg3, fitBkg4, constrainSignal, constrainBkg1, constrainBkg2, forceGaussianSignal, param, param_err, canvas);

    /* check requested pt-range */
    if (requestedRange)
      return;

    /* write canvas */
    fitDir->cd();
    canvas->Write(Form("fitDisplay_ptBin_%3.2f_%3.2f", ptBin[ipt], ptBin[ipt + 1]));

    /* set histo */
    for (Int_t iparam = 0; iparam < kNFitParams; iparam++) {
      hFitParamHisto[iparam]->SetBinContent(ipt + 1, param[iparam]);
      hFitParamHisto[iparam]->SetBinError(ipt + 1, param_err[iparam]);
    }

    /* delete */
    delete hSignal_py;
    delete hMismatch_py;
    for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
      delete hSignalExp_py[iiipart];
      //      delete hSignalExpTail_py[iiipart];
    }
  }

  /* check requested pt-range */
  if (requestedRange)
    return;

  /*** POST-ANALYSIS ***/

  TDirectory *postDir = fileout->mkdir("PostAnalysis");
 
  /* compute overflows */
  TH1D *hOverflowCounts = new TH1D(*hFitParamHisto[kTotalCounts]);
  hOverflowCounts->SetNameTitle("hOverflowCounts", "Overflow counts: TotalCounts - IntegralCounts;p_{T} (GeV/c);counts");
  hOverflowCounts->Add(hFitParamHisto[kIntegralCounts], -1.);
  /* compute total mismatch counts */
  TH1D *hTotalMismatchCounts = new TH1D(*hFitParamHisto[kMismatchCounts]);
  hTotalMismatchCounts->SetNameTitle("hTotalMismatchCounts", "Total mismatch counts: MismatchCounts + OverflowCounts;p_{T} (GeV/c);counts");
  hTotalMismatchCounts->Add(hOverflowCounts);
  /* computed mismatch fraction */
  TH1D *hTotalMismatchFraction = new TH1D(*hTotalMismatchCounts);
  hTotalMismatchFraction->SetNameTitle("hTotalMismatchFraction", "Total mismatch fraction: TotalMismatchCounts / TotalCounts;p_{T} (GeV/c);");
  hTotalMismatchFraction->Divide(hFitParamHisto[kTotalCounts]);
  /* compute identified counts */
  TH1D *hIdentifiedCounts = new TH1D(*hFitParamHisto[kSignalCounts]);
  hIdentifiedCounts->SetNameTitle("hIdentifiedCounts", "Identified counts: SignalCounts + sum(BkgCounts);p_{T} (GeV/c);counts");
  hIdentifiedCounts->Add(hFitParamHisto[kBkg1Counts]);
  hIdentifiedCounts->Add(hFitParamHisto[kBkg2Counts]);
  hIdentifiedCounts->Add(hFitParamHisto[kBkg3Counts]);
  /* compute signal fraction */
  TH1D *hSignalFraction = new TH1D(*hFitParamHisto[kSignalCounts]);
  hSignalFraction->SetNameTitle("hSignalFraction", "Signal fraction: SignalCounts / IdentifiedCounts;p_{T} (GeV/c);");
  hSignalFraction->Divide(hSignalFraction, hIdentifiedCounts, 1., 1., "B");
  /* compute bkg1 fraction */
  TH1D *hBkg1Fraction = new TH1D(*hFitParamHisto[kBkg1Counts]);
  hBkg1Fraction->SetNameTitle("hBkg1Fraction", "Bkg1 fraction: Bkg1Counts / IdentifiedCounts;p_{T} (GeV/c);");
  hBkg1Fraction->Divide(hBkg1Fraction, hIdentifiedCounts, 1., 1., "B");
  /* compute bkg2 fraction */
  TH1D *hBkg2Fraction = new TH1D(*hFitParamHisto[kBkg2Counts]);
  hBkg2Fraction->SetNameTitle("hBkg2Fraction", "Bkg2 fraction: Bkg2Counts / IdentifiedCounts;p_{T} (GeV/c);");
  hBkg2Fraction->Divide(hBkg2Fraction, hIdentifiedCounts, 1., 1., "B");
  /* compute bkg3 fraction */
  TH1D *hBkg3Fraction = new TH1D(*hFitParamHisto[kBkg3Counts]);
  hBkg3Fraction->SetNameTitle("hBkg3Fraction", "Bkg3 fraction: Bkg3Counts / IdentifiedCounts;p_{T} (GeV/c);");
  hBkg3Fraction->Divide(hBkg3Fraction, hIdentifiedCounts, 1., 1., "B");
  /* compute bkg4 fraction */
  TH1D *hBkg4Fraction = new TH1D(*hFitParamHisto[kBkg4Counts]);
  hBkg4Fraction->SetNameTitle("hBkg4Fraction", "Bkg4 fraction: Bkg4Counts / IdentifiedCounts;p_{T} (GeV/c);");
  hBkg4Fraction->Divide(hBkg4Fraction, hIdentifiedCounts, 1., 1., "B");
  /* compute mismatch-correction counts */
  TH1D *hMismatchCorrectionCounts = new TH1D(*hTotalMismatchCounts);
  hMismatchCorrectionCounts->SetNameTitle("hMismatchCorrectionCounts", "Mismatch-correction counts: TotalMismatchCounts * SignalFraction;p_{T} (GeV/c);counts");
  hMismatchCorrectionCounts->Multiply(hSignalFraction);
  /* compute mismatch-corrected signal counts */
  TH1D *hMismatchCorrectedSignalCounts = new TH1D(*hFitParamHisto[kSignalCounts]);
  hMismatchCorrectedSignalCounts->SetNameTitle("hMismatchCorrectedSignalCounts", "Mismatch-corrected signal counts: SignalCounts + MismatchCorrectionCounts;p_{T} (GeV/c);counts");
  hMismatchCorrectedSignalCounts->Add(hMismatchCorrectionCounts);

  /* accepted tracks histo */
  if (hAcceptedTracks) {
    TH1D *hAcceptedTracks_py = hAcceptedTracks->ProjectionY();
    hAcceptedTracks_py->SetNameTitle("hAcceptedTracks", "Accepted tracks;p_{T} (GeV/c);");
    hAcceptedTracks_py->Sumw2();
  }

  /*** RAW SPECTRA ***/

  TDirectory *rawDir = fileout->mkdir("RawSpectra");

  /* compute normalized raw yield */
  TH1D *hNormalizedRawYield = new TH1D(*hFitParamHisto[kSignalCounts]);
  hNormalizedRawYield->SetNameTitle("hNormalizedRawYield", "Raw yield;p_{T} (GeV/c);#frac{1}{N_{ev}} #frac{d^{2}N}{dy dp_{T}}");
  TOFpid_normalize(hNormalizedRawYield, nevents);
  /* compute normalized mismatch-corrected raw yield */
  TH1D *hNormalizedMismatchCorrectedRawYield = new TH1D(*hMismatchCorrectedSignalCounts);
  hNormalizedMismatchCorrectedRawYield->SetNameTitle("hNormalizedMismatchCorrectedRawYield", "Mismatch-corrected raw yield;p_{T} (GeV/c);#frac{1}{N_{ev}} #frac{d^{2}N}{dy dp_{T}}");
  TOFpid_normalize(hNormalizedMismatchCorrectedRawYield, nevents);

  /*** OUTPUT ***/

  /* write fir params histos */
  fitDir->cd();
  for (Int_t iparam = 0; iparam < kNFitParams; iparam++)
    hFitParamHisto[iparam]->Write();
  /* write post-analysis histos */
  postDir->cd();
  hOverflowCounts->Write();
  hTotalMismatchCounts->Write();
  hTotalMismatchFraction->Write();
  hIdentifiedCounts->Write();
  hSignalFraction->Write();
  hBkg1Fraction->Write();
  hBkg2Fraction->Write();
  hBkg3Fraction->Write();
  hBkg4Fraction->Write();
  hMismatchCorrectionCounts->Write();
  hMismatchCorrectedSignalCounts->Write();
  if (hAcceptedTracks) hAcceptedTracks_py->Write();
  /* write raw spectra histos */
  rawDir->cd();
  hNormalizedRawYield->Write();
  hNormalizedMismatchCorrectedRawYield->Write();

  /* clean up */
  delete canvas;
  for (Int_t iparam = 0; iparam < kNFitParams; iparam++)
    delete hFitParamHisto[iparam];
  delete hOverflowCounts;
  delete hTotalMismatchCounts;
  delete hTotalMismatchFraction;
  delete hIdentifiedCounts;
  delete hSignalFraction;
  delete hBkg1Fraction;
  delete hBkg2Fraction;
  delete hBkg3Fraction;
  delete hBkg4Fraction;
  delete hMismatchCorrectionCounts;
  delete hMismatchCorrectedSignalCounts;
  if (hAcceptedTracks) {
    delete hAcceptedTracks_py;
  }
  delete hNormalizedRawYield;
  delete hNormalizedMismatchCorrectedRawYield;

  /* close file */
  fileout->Close();
}

//___________________________________________________________________________________

Float_t 
TOFpid_histomin(TH1 *h)
{

  for (Int_t ibin = 0; ibin < h->GetNbinsX(); ibin++)
    if (h->GetBinContent(ibin + 1) > 0.)
      return h->GetXaxis()->GetBinCenter(ibin + 1);
  return kMaxInt;
}

//___________________________________________________________________________________

Float_t 
TOFpid_histomax(TH1 *h)
{

  for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--)
    if (h->GetBinContent(ibin) > 0.)
      return h->GetXaxis()->GetBinCenter(ibin);
  return -kMaxInt;
}

//___________________________________________________________________________________

TOFpid_checkneg(TH1 *h)
{

  for (Int_t ibin = 0; ibin < h->GetNbinsX(); ibin++)
    if (h->GetBinContent(ibin + 1) <= 0.) {
      h->SetBinContent(ibin + 1, 1.e-300);
      //      h->SetBinError(ibin + 1, 0.1);
    }
}

//___________________________________________________________________________________

Double_t
TOFpid_fit(TH1 *hSignal, TH1 *hSigExp, TH1 *hBkgExp1, TH1 *hBkgExp2, TH1 *hBkgExp3, TH1 *hBkgExp4, TH1 *hMismatch, Double_t rangeMin, Double_t rangeMax, Bool_t fitBkg1, Bool_t fitBkg2, Bool_t fitBkg3, Bool_t fitBkg4, Bool_t constrainSignal, Bool_t constrainBkg1, Bool_t constrainBkg2, Bool_t forceGaussianSignal, Double_t *param = NULL, Double_t *param_err = NULL, TCanvas *canvas = NULL)
{

  /** ROOFIT ***/
  gSystem->Load("libRooFit");
  using namespace RooFit;
 /*** LOAD GAUSSIANTAIL CLASS ***/
  gSystem->Load("~/ALICE.2011/ANALYSIS/TOFSpectraPbPb/macros/RooFermiCutoff_cxx");
  gSystem->Load("~/ALICE.2011/ANALYSIS/TOFSpectraPbPb/macros/RooGaussianTail_cxx");

  /*** DEFINE FIT RANGE ***/

  printf("***** FIT RANGE DEFINITION *****\n");

  /* check mismatch histogram to define min/max fit range */
  //  rangeMin = TMath::Max(rangeMin, TOFpid_histomin(hMismatch));
  //  rangeMax = TMath::Min(rangeMax, TOFpid_histomax(hMismatch));
  /* fix zeroes */
  TOFpid_checkneg(hMismatch);
  
  /* define range */
  RooRealVar x("x", "n_{#sigma}", 0., rangeMin, rangeMax, "");
  printf("FIT RANGE DEFINED: %f -> %f\n", rangeMin, rangeMax);
  printf("********************************\n");

  /*** DEFINE HISTOGRAM DATA ***/
  
  /* define data to fit and background from input histogram */
  RooDataHist hdata("hdata", "hdata", x, hSignal);
  RooDataHist hsig("hsig", "hsig", x, hSigExp);
  RooDataHist hbkg1("hbkg1", "hbkg1", x, hBkgExp1);
  RooDataHist hbkg2("hbkg2", "hbkg2", x, hBkgExp2);
  RooDataHist hbkg3("hbkg3", "hbkg3", x, hBkgExp3);
  RooDataHist hbkg4("hbkg4", "hbkg4", x, hBkgExp4);
  RooDataHist hmismatch("hmismatch", "hmismatch", x, hMismatch);

  /*** DEFINE SIGNAL SHAPE ***/

  printf("***** SIGNAL SHAPE DEFINITION *****\n");

  /* variables */
  //  if (FIX_SIGNAL_MEAN) {
  //    RooRealVar mean("mean", "mean", DEFAULT_SIGNAL_MEAN, "");
  //    printf("FIXED SIGNAL_MEAN = %f\n", mean.getVal());
  //  }
  //  else {
  RooRealVar mean("mean", "mean", DEFAULT_SIGNAL_MEAN, MIN_SIGNAL_MEAN, MAX_SIGNAL_MEAN, "");
  //    printf("FREE SIGNAL_MEAN = %f [%f, %f]\n", mean.getVal(), MIN_SIGNAL_MEAN, MAX_SIGNAL_MEAN);
  //  }
  //  if (FIX_SIGNAL_SIGMA) {
  //    RooRealVar sigma("sigma", "sigma", DEFAULT_SIGNAL_SIGMA, "");
  //    printf("FIXED SIGNAL_SIGMA = %f\n", sigma.getVal());
  //  }
  //  else {
    RooRealVar sigma("sigma", "sigma", DEFAULT_SIGNAL_SIGMA, MIN_SIGNAL_SIGMA, MAX_SIGNAL_SIGMA), "");
//    printf("FREE SIGNAL_SIGMA = %f\n", sigma.getVal());
//  }
//  if (FIX_SIGNAL_TAIL) {
//    RooRealVar tail("tail", "tail", DEFAULT_SIGNAL_TAIL, "");
//    printf("FIXED SIGNAL_TAIL = %f\n", tail.getVal());
//  }
//  else {
    RooRealVar tail("tail", "tail", DEFAULT_SIGNAL_TAIL, MIN_SIGNAL_TAIL, MAX_SIGNAL_TAIL, "");
//    printf("FREE SIGNAL_TAIL = %f\n", tail.getVal());
//  }
  RooRealVar gaussianfrac("gaussianfrac", "gaussianfrac", 1., 0., 1., "");
  RooRealVar sigalpha("sigalpha", "sigalpha", 0., -10., 0.);

  /* shapes */
  if (GAUSSIAN_SIGNAL || forceGaussianSignal) {
    printf("USING GAUSSIAN_SIGNAL SHAPE\n");
    RooGaussian signal("signal", "signal", x, mean, sigma);
  }
  else if (GAUSSIANTAIL_SIGNAL) {
    printf("USING GAUSSIANTAIL_SIGNAL SHAPE\n");
    RooGaussianTail signal("signal", "signal", x, mean, sigma, tail);
  }
  else if (GAUSSIANTAIL2_SIGNAL) {
    printf("USING GAUSSIANTAIL2_SIGNAL SHAPE\n");
    RooGaussianTail signal("signal", "signal", x, mean, sigma, sigma);
  }
  else if (GAUSSIANPLUSGAUSSIANTAIL_SIGNAL) {
    printf("USING GAUSSIANPLUSGAUSSIANTAIL_SIGNAL SHAPE\n");
    RooGaussian gaussian("gaussian", "gaussian", x, mean, sigma);
    RooGaussianTail gaussiantail("gaussiantail", "gaussiantail", x, mean, sigma, tail);
    RooAddPdf signal("signal", "signal", RooArgList(gaussian, gaussiantail), gaussianfrac, kTRUE);

  }
  else if (GAUSSIANPLUSEXPONENTIAL_SIGNAL) {
    printf("USING GAUSSIANPLUSEXPONENTIAL_SIGNAL SHAPE\n");
    RooGaussian gaussian("gaussian", "gaussian", x, mean, sigma);
    RooExponential sigexpo("sigexpo", "sigexpo", x, sigalpha);
    RooRealVar sigcutoff("sigcutoff", "sigcutoff", 0.);
    RooRealVar sigpower("sigpower", "sigpower", 0.1);
    RooFermiCutoff sigfermi("sigfermi", "sigfermi", x, sigcutoff, sigpower);
    RooProdPdf exposignal("exposignal", "exposignal", sigfermi, sigexpo, -100.);
    RooRealVar gaussianfrac("gaussianfrac", "gaussianfrac", 0.9, 0.7, 1., "");
    RooAddPdf signal("signal", "signal", RooArgList(gaussian, exposignal), gaussianfrac, kTRUE);
  }
  else if (EXPECTED_SIGNAL_TEMPLATE) {
    printf("SHAPE OF SIGNAL FROM TEMPLATE HISTOGRAM\n");
    RooHistPdf signal("signal", "signal", x, hsig);
  }
  else if (EXPECTED_SIGNAL_FIT) {
    /* fitting bkg1 */
    TF1 *fGaus = (TF1 *)gROOT->GetFunction("gaus");
    hSigExp->Fit(fGaus, "0q");
    Double_t sig_mean = fGaus->GetParameter(1);
    Double_t sig_sigma = fGaus->GetParameter(2);
    mean.setVal(sig_mean);
    mean.setRange(sig_mean - 10., sig_mean + 10.);
    sigma.setVal(sig_sigma);
    sigma.setRange(sig_sigma * 0.5, sig_sigma * 1.5);
    tail.setVal(1.);
    tail.setRange(0.5, 5.);
    RooGaussianTail signal("signal", "signal", x, mean, sigma, tail);
    signal.fitTo(hsig, Range(sig_mean - 5. * sig_sigma, sig_mean + 10. * sig_sigma), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
#if DISPLAY
    TCanvas *cSig_fit = new TCanvas("cSig_fit");
    RooPlot *sig_frame = x.frame();
    hsig.plotOn(sig_frame);
    signal.plotOn(sig_frame, LineColor(kRed));
    sig_frame->Draw();
    cSig_fit->Update();
#endif
    printf("SIGNAL PARAMETERS AFTER FIT OF EXPECTED SIGNAL\n");
    printf("mean  = %f +- %f\n", mean.getVal(), mean.getError());
    printf("sigma = %f +- %f\n", sigma.getVal(), sigma.getError());
    printf("tail  = %f +- %f\n", tail.getVal(), tail.getError());
    /* scale parameters if requested */
    if (SCALE_SIGNAL_SIGMA != 1.) {
      printf("SCALE FITTED SIGNAL SIGMA BY %f\n", SCALE_SIGNAL_SIGMA);
      sigma.setVal(sigma.getVal() * SCALE_SIGNAL_SIGMA);
    }
    if (SCALE_SIGNAL_TAIL != 1.) {
      printf("SCALE FITTED SIGNAL TAIL BY %f\n", SCALE_SIGNAL_TAIL);
      tail.setVal(tail.getVal() * SCALE_SIGNAL_TAIL);
    }
    /* fix/release parameters if requested */
    if (FIX_SIGNAL_MEAN) {
      printf("SETTING FITTED SIGNAL MEAN AS CONSTANTS\n");
      mean.setConstant(kTRUE);
    }
    else {
      printf("SETTING FITTED SIGNAL MEAN AS FREE\n");
      //      mean.setRange(mean.getVal() - 0.25 * TMath::Abs(mean.getVal()), mean.getVal() + 0.25 * TMath::Abs(mean.getVal()));
      //      mean.setRange(mean.getVal() - 0.5 * TMath::Abs(mean.getVal()), mean.getVal() + 0.5 * TMath::Abs(mean.getVal()));
      mean.setRange(-0.1, 0.1);
    }
    if (FIX_SIGNAL_SIGMA) {
      printf("SETTING FITTED SIGNAL SIGMA AS CONSTANTS\n");
      sigma.setConstant(kTRUE);
    }
    else {
      printf("SETTING FITTED SIGNAL SIGMA AS FREE\n");
      sigma.setRange(sigma.getVal() * 0.75, sigma.getVal() * 1.25);
      //      sigma.setRange(sigma.getVal() - 0.5, sigma.getVal() + 0.5);
    }
    if (FIX_SIGNAL_TAIL) {
      printf("SETTING FITTED SIGNAL TAIL AS CONSTANTS\n");
      tail.setConstant(kTRUE);
    }
    else {
      printf("SETTING FITTED SIGNAL TAIL AS FREE\n");
      tail.setRange(0.75, 2.0);
      //      tail.setRange(tail.getVal() * 0.75, tail.getVal() * 1.25);
      //      tail.setRange(tail.getVal() - 0.5, tail.getVal() + 0.5);
    }
  }
  else {
    printf("SHAPE OF SIGNAL NOT DEFINE: using GAUSSIAN_SIGNAL\n");
    RooGaussian signal("signal", "signal", x, mean, sigma);
  }


#if 0
  if (constrainSignal) {
#if 0
  /* fix expected signal and constrain parameters if requested */
  signal.fitTo(hsig);
#if 0
  TCanvas *cConstrainSignal = new TCanvas("cConstrainSignal");
  RooPlot *xfr = x.frame();
  hsig.plotOn(xfr);
  signal.plotOn(xfr, LineColor(kRed));
  xfr->Draw();
  cConstrainSignal->Update();
#endif
  printf("SIGNAL PARAMETERS AFTER FIT OF EXPECTED SIGNAL\n");
  printf("mean  = %f +- %f\n", mean.getVal(), mean.getError());
  printf("sigma = %f +- %f\n", sigma.getVal(), sigma.getError());
  printf("tail  = %f +- %f\n", tail.getVal(), tail.getError());
  if (constrainSignal) {
    mean.setConstant(kTRUE);
    sigma.setConstant(kTRUE);
    tail.setConstant(kTRUE);
  printf("SIGNAL PARAMETERS CONSTRAINED AFTER FIT OF EXPECTED SIGNAL\n");
  }
#endif
  }
#endif

  if (constrainSignal) {
    mean.setConstant(kTRUE);
    sigma.setConstant(kTRUE);
    tail.setConstant(kTRUE);
    //    mean.setRange(-0.1, 0.1);
    //    sigma.setRange(0.95, 1.05);
    //    tail.setRange(0.95, 1.25);
  }

  printf("***********************************\n");

  /*** DEFINE IDENTIFIED BACKGROUND SHAPES ***/

  printf("***** IDENTIFIED BACKGROUND SHAPE DEFINITION *****\n");

  /* shapes */
if (EXPECTED_BACKGROUND_TEMPLATE) {
  printf("USING EXPECTED BACKGROUND TEMPLATE SHAPES FROM HISTOGRAMS\n");
    RooHistPdf bkg1("bkg1", "bkg1", x, hbkg1, 0);
    RooHistPdf bkg2("bkg2", "bkg2", x, hbkg2, 0);
      RooHistPdf bkg3("bkg3", "bkg3", x, hbkg3, 0);
      RooHistPdf bkg4("bkg4", "bkg4", x, hbkg4, 0);
  }
 else if (EXPECTED_BACKGROUND_FIT) {
    printf("USING EXPECTED BACKGROUND FITTED SHAPES FROM HISTOGRAMS\n");
    /* fitting bkg1 */
    TF1 *fGaus = (TF1 *)gROOT->GetFunction("gaus");
    hBkgExp1->Fit(fGaus, "0q");
    Double_t bkgexp1_mean = fGaus->GetParameter(1);
    Double_t bkgexp1_sigma = fGaus->GetParameter(2);
    RooRealVar mean_bkg1("mean_bkg1", "mean_bkg1", bkgexp1_mean, bkgexp1_mean - 10., bkgexp1_mean + 10., "");
    RooRealVar sigma_bkg1("sigma_bkg1", "sigma_bkg1", bkgexp1_sigma, bkgexp1_sigma * 0.5, bkgexp1_sigma * 1.5, "");
    RooRealVar tail_bkg1("tail_bkg1", "tail_bkg1", 1.0, 0.5, 5., "");
    RooGaussianTail bkg1("bkg1", "bkg1", x, mean_bkg1, sigma_bkg1, tail_bkg1);
    bkg1.fitTo(hbkg1, Range(bkgexp1_mean - 5. * bkgexp1_sigma, bkgexp1_mean + 10. * bkgexp1_sigma), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
#if DISPLAY
    TCanvas *cBkg1_fit = new TCanvas("cBkg1_fit");
    RooPlot *bkg1_frame = x.frame();
    hbkg1.plotOn(bkg1_frame);
    bkg1.plotOn(bkg1_frame, LineColor(kCyan+1));
    bkg1_frame->Draw();
    cBkg1_fit->Update();
#endif
    printf("BACKGROUND-1 PARAMETERS AFTER FIT OF EXPECTED BACKGROUND-1\n");
    printf("mean_bkg1  = %f +- %f\n", mean_bkg1.getVal(), mean_bkg1.getError());
    printf("sigma_bkg1 = %f +- %f\n", sigma_bkg1.getVal(), sigma_bkg1.getError());
    printf("tail_bkg1  = %f +- %f\n", tail_bkg1.getVal(), tail_bkg1.getError());
    /* scale parameters if requested */
    if (SCALE_BACKGROUND_SIGMA != 1.) {
      printf("SCALE FITTED BACKGROUND-1 SIGMA BY %f\n", SCALE_BACKGROUND_SIGMA);
      sigma_bkg1.setVal(sigma_bkg1.getVal() * SCALE_BACKGROUND_SIGMA);
    }
    if (SCALE_BACKGROUND_TAIL != 1.) {
      printf("SCALE FITTED BACKGROUND-1 TAIL BY %f\n", SCALE_BACKGROUND_TAIL);
      tail_bkg1.setVal(tail_bkg1.getVal() * SCALE_BACKGROUND_TAIL);
    }
    /* fix/release parameters if requested */
    if (FIX_BACKGROUND_MEAN) {
      printf("SETTING BACKGROUND-1 FITTED MEAN AS CONSTANTS\n");
      mean_bkg1.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-1 FITTED MEAN AS FREE\n");
      mean_bkg1.setRange(mean_bkg1.getVal() - 0.25 * TMath::Abs(mean_bkg1.getVal()), mean_bkg1.getVal() + 0.25 * TMath::Abs(mean_bkg1.getVal()));
    }
    if (FIX_BACKGROUND_SIGMA) {
      printf("SETTING BACKGROUND-1 FITTED SIGMA AS CONSTANTS\n");
      sigma_bkg1.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-1 FITTED SIGMA AS FREE\n");
      sigma_bkg1.setRange(sigma_bkg1.getVal() * 0.75, sigma_bkg1.getVal() * 1.25);
    }
    if (FIX_BACKGROUND_TAIL) {
      printf("SETTING BACKGROUND-1 FITTED TAIL AS CONSTANTS\n");
      tail_bkg1.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-1 FITTED TAIL AS FREE\n");
      tail_bkg1.setRange(tail_bkg1.getVal() * 0.75, tail_bkg1.getVal() * 1.25);
    }
    /* fitting bkg2 */
    TF1 *fGaus = (TF1 *)gROOT->GetFunction("gaus");
    hBkgExp2->Fit(fGaus, "0q");
    Double_t bkgexp2_mean = fGaus->GetParameter(1);
    Double_t bkgexp2_sigma = fGaus->GetParameter(2);
    RooRealVar mean_bkg2("mean_bkg2", "mean_bkg2", bkgexp2_mean, bkgexp2_mean - 10., bkgexp2_mean + 10., "");
    RooRealVar sigma_bkg2("sigma_bkg2", "sigma_bkg2", bkgexp2_sigma, bkgexp2_sigma * 0.5, bkgexp2_sigma * 1.5, "");
    RooRealVar tail_bkg2("tail_bkg2", "tail_bkg2", 1.0, 0.5, 5., "");
    RooGaussianTail bkg2("bkg2", "bkg2", x, mean_bkg2, sigma_bkg2, tail_bkg2);
    bkg2.fitTo(hbkg2, Range(bkgexp2_mean - 5. * bkgexp2_sigma, bkgexp2_mean + 10. * bkgexp2_sigma), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
#if DISPLAY
    TCanvas *cBkg2_fit = new TCanvas("cBkg2_fit");
    RooPlot *bkg2_frame = x.frame();
    hbkg2.plotOn(bkg2_frame);
    bkg2.plotOn(bkg2_frame,  LineColor(kGreen+1));
    bkg2_frame->Draw();
    cBkg2_fit->Update();
#endif
    printf("BACKGROUND-2 PARAMETERS AFTER FIT OF EXPECTED BACKGROUND-2\n");
    printf("mean_bkg2  = %f +- %f\n", mean_bkg2.getVal(), mean_bkg2.getError());
    printf("sigma_bkg2 = %f +- %f\n", sigma_bkg2.getVal(), sigma_bkg2.getError());
    printf("tail_bkg2  = %f +- %f\n", tail_bkg2.getVal(), tail_bkg2.getError());
    /* scale parameters if requested */
    if (SCALE_BACKGROUND_SIGMA != 1.) {
      printf("SCALE FITTED BACKGROUND-2 SIGMA BY %f\n", SCALE_BACKGROUND_SIGMA);
      sigma_bkg2.setVal(sigma_bkg2.getVal() * SCALE_BACKGROUND_SIGMA);
    }
    if (SCALE_BACKGROUND_TAIL != 1.) {
      printf("SCALE FITTED BACKGROUND-2 TAIL BY %f\n", SCALE_BACKGROUND_TAIL);
      tail_bkg2.setVal(tail_bkg2.getVal() * SCALE_BACKGROUND_TAIL);
    }
    /* fix/release parameters if requested */
    if (FIX_BACKGROUND_MEAN) {
      printf("SETTING BACKGROUND-2 FITTED MEAN AS CONSTANTS\n");
      mean_bkg2.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-2 FITTED MEAN AS FREE\n");
      mean_bkg2.setRange(mean_bkg2.getVal() - 0.25 * TMath::Abs(mean_bkg2.getVal()), mean_bkg2.getVal() + 0.25 * TMath::Abs(mean_bkg2.getVal()));
    }
    if (FIX_BACKGROUND_SIGMA) {
      printf("SETTING BACKGROUND-2 FITTED SIGMA AS CONSTANTS\n");
      sigma_bkg2.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-2 FITTED SIGMA AS FREE\n");
      sigma_bkg2.setRange(sigma_bkg2.getVal() * 0.75, sigma_bkg2.getVal() * 1.25);
    }
    if (FIX_BACKGROUND_TAIL) {
      printf("SETTING BACKGROUND-2 FITTED TAIL AS CONSTANTS\n");
      tail_bkg2.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-2 FITTED TAIL AS FREE\n");
      tail_bkg2.setRange(tail_bkg2.getVal() * 0.75, tail_bkg2.getVal() * 1.25);
    }
    /* fitting bkg3 */
    TF1 *fGaus = (TF1 *)gROOT->GetFunction("gaus");
    hBkgExp3->Fit(fGaus, "0q");
    Double_t bkgexp3_mean = fGaus->GetParameter(1);
    Double_t bkgexp3_sigma = fGaus->GetParameter(2);
    RooRealVar mean_bkg3("mean_bkg3", "mean_bkg3", bkgexp3_mean, bkgexp3_mean - 10., bkgexp3_mean + 10., "");
    RooRealVar sigma_bkg3("sigma_bkg3", "sigma_bkg3", bkgexp3_sigma, bkgexp3_sigma * 0.5, bkgexp3_sigma * 1.5, "");
    RooRealVar tail_bkg3("tail_bkg3", "tail_bkg3", 1., 0.5, 5., "");
    RooGaussianTail bkg3("bkg3", "bkg3", x, mean_bkg3, sigma_bkg3, tail_bkg3);
    bkg3.fitTo(hbkg3, Range(bkgexp3_mean - 5. * bkgexp3_sigma, bkgexp3_mean + 10. * bkgexp3_sigma), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
#if DISPLAY
    TCanvas *cBkg3_fit = new TCanvas("cBkg3_fit");
    RooPlot *bkg3_frame = x.frame();
    hbkg3.plotOn(bkg3_frame);
    bkg3.plotOn(bkg3_frame,  LineColor(kYellow+1));
    bkg3_frame->Draw();
    cBkg3_fit->Update();
#endif
    printf("BACKGROUND-3 PARAMETERS AFTER FIT OF EXPECTED BACKGROUND-3\n");
    printf("mean_bkg3  = %f +- %f\n", mean_bkg3.getVal(), mean_bkg3.getError());
    printf("sigma_bkg3 = %f +- %f\n", sigma_bkg3.getVal(), sigma_bkg3.getError());
    printf("tail_bkg3  = %f +- %f\n", tail_bkg3.getVal(), tail_bkg3.getError());
    /* scale parameters if requested */
    if (SCALE_BACKGROUND_SIGMA != 1.) {
      printf("SCALE FITTED BACKGROUND-3 SIGMA BY %f\n", SCALE_BACKGROUND_SIGMA);
      sigma_bkg3.setVal(sigma_bkg3.getVal() * SCALE_BACKGROUND_SIGMA);
    }
    if (SCALE_BACKGROUND_TAIL != 1.) {
      printf("SCALE FITTED BACKGROUND-3 TAIL BY %f\n", SCALE_BACKGROUND_TAIL);
      tail_bkg3.setVal(tail_bkg3.getVal() * SCALE_BACKGROUND_TAIL);
    }
    /* fix/release parameters if requested */
    if (FIX_BACKGROUND_MEAN) {
      printf("SETTING BACKGROUND-3 FITTED MEAN AS CONSTANTS\n");
      mean_bkg3.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-3 FITTED MEAN AS FREE\n");
      mean_bkg3.setRange(mean_bkg3.getVal() - 0.25 * TMath::Abs(mean_bkg3.getVal()), mean_bkg3.getVal() + 0.25 * TMath::Abs(mean_bkg3.getVal()));
    }
    if (FIX_BACKGROUND_SIGMA) {
      printf("SETTING BACKGROUND-3 FITTED SIGMA AS CONSTANTS\n");
      sigma_bkg3.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-3 FITTED SIGMA AS FREE\n");
      sigma_bkg3.setRange(sigma_bkg3.getVal() * 0.75, sigma_bkg3.getVal() * 1.25);
    }
    if (FIX_BACKGROUND_TAIL) {
      printf("SETTING BACKGROUND-3 FITTED TAIL AS CONSTANTS\n");
      tail_bkg3.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-3 FITTED TAIL AS FREE\n");
      tail_bkg3.setRange(tail_bkg3.getVal() * 0.75, tail_bkg3.getVal() * 1.25);
    }
    /* fitting bkg4 */
    TF1 *fGaus = (TF1 *)gROOT->GetFunction("gaus");
    hBkgExp4->Fit(fGaus, "0q");
    Double_t bkgexp4_mean = fGaus->GetParameter(1);
    Double_t bkgexp4_sigma = fGaus->GetParameter(2);
    RooRealVar mean_bkg4("mean_bkg4", "mean_bkg4", bkgexp4_mean, bkgexp4_mean - 10., bkgexp4_mean + 10., "");
    RooRealVar sigma_bkg4("sigma_bkg4", "sigma_bkg4", bkgexp4_sigma, bkgexp4_sigma * 0.5, bkgexp4_sigma * 1.5, "");
    RooRealVar tail_bkg4("tail_bkg4", "tail_bkg4", 1., 0.5, 5., "");
    RooGaussianTail bkg4("bkg4", "bkg4", x, mean_bkg4, sigma_bkg4, tail_bkg4);
    bkg4.fitTo(hbkg4, Range(bkgexp4_mean - 5. * bkgexp4_sigma, bkgexp4_mean + 10. * bkgexp4_sigma), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));
#if DISPLAY
    TCanvas *cBkg4_fit = new TCanvas("cBkg4_fit");
    RooPlot *bkg4_frame = x.frame();
    hbkg4.plotOn(bkg4_frame);
    bkg4.plotOn(bkg4_frame,  LineColor(kYellow+2));
    bkg4_frame->Draw();
    cBkg4_fit->Update();
#endif
    printf("BACKGROUND-4 PARAMETERS AFTER FIT OF EXPECTED BACKGROUND-4\n");
    printf("mean_bkg4  = %f +- %f\n", mean_bkg4.getVal(), mean_bkg4.getError());
    printf("sigma_bkg4 = %f +- %f\n", sigma_bkg4.getVal(), sigma_bkg4.getError());
    printf("tail_bkg4  = %f +- %f\n", tail_bkg4.getVal(), tail_bkg4.getError());
    /* scale parameters if requested */
    if (SCALE_BACKGROUND_SIGMA != 1.) {
      printf("SCALE FITTED BACKGROUND-4 SIGMA BY %f\n", SCALE_BACKGROUND_SIGMA);
      sigma_bkg4.setVal(sigma_bkg4.getVal() * SCALE_BACKGROUND_SIGMA);
    }
    if (SCALE_BACKGROUND_TAIL != 1.) {
      printf("SCALE FITTED BACKGROUND-4 TAIL BY %f\n", SCALE_BACKGROUND_TAIL);
      tail_bkg4.setVal(tail_bkg4.getVal() * SCALE_BACKGROUND_TAIL);
    }
    /* fix/release parameters if requested */
    if (FIX_BACKGROUND_MEAN) {
      printf("SETTING BACKGROUND-4 FITTED MEAN AS CONSTANTS\n");
      mean_bkg4.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-4 FITTED MEAN AS FREE\n");
      mean_bkg4.setRange(mean_bkg4.getVal() - 0.25 * TMath::Abs(mean_bkg4.getVal()), mean_bkg4.getVal() + 0.25 * TMath::Abs(mean_bkg4.getVal()));
    }
    if (FIX_BACKGROUND_SIGMA) {
      printf("SETTING BACKGROUND-4 FITTED SIGMA AS CONSTANTS\n");
      sigma_bkg4.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-4 FITTED SIGMA AS FREE\n");
      sigma_bkg4.setRange(sigma_bkg4.getVal() * 0.75, sigma_bkg4.getVal() * 1.25);
    }
    if (FIX_BACKGROUND_TAIL) {
      printf("SETTING BACKGROUND-4 FITTED TAIL AS CONSTANTS\n");
      tail_bkg4.setConstant(kTRUE);
    }
    else {
      printf("SETTING BACKGROUND-4 FITTED TAIL AS FREE\n");
      tail_bkg4.setRange(tail_bkg4.getVal() * 0.75, tail_bkg4.getVal() * 1.25);
    }
  }
  else if (GAUSSIAN_BACKGROUND) {
    printf("USING GAUSSIAN BACKGROUND SHAPES (not reccomended)\n"); 
    RooRealVar mean1("mean1", "mean1", hBkgExp1->GetMean(), hBkgExp1->GetMean() * 0.95, hBkgExp1->GetMean() * 1.05, "");
    RooRealVar sigma1("sigma1", "sigma1", hBkgExp1->GetRMS(), hBkgExp1->GetRMS() * 0.5, hBkgExp1->GetRMS() * 1.5, "");
    RooGaussian bkg1("bkg1", "bkg1", x, mean1, sigma1);
    
    RooRealVar mean2("mean2", "mean2", hBkgExp2->GetMean(), hBkgExp2->GetMean() * 0.95, hBkgExp2->GetMean() * 1.05, "");
    RooRealVar sigma2("sigma2", "sigma2", hBkgExp2->GetRMS(), hBkgExp2->GetRMS() * 0.5, hBkgExp2->GetRMS() * 1.5, "");
    RooGaussian bkg2("bkg2", "bkg2", x, mean2, sigma2);

    RooRealVar mean3("mean3", "mean3", hBkgExp3->GetMean(), hBkgExp3->GetMean() * 0.95, hBkgExp3->GetMean() * 1.05, "");
    RooRealVar sigma3("sigma3", "sigma3", hBkgExp3->GetRMS(), hBkgExp3->GetRMS() * 0.5, hBkgExp3->GetRMS() * 1.5, "");
    RooGaussian bkg3("bkg3", "bkg3", x, mean3, sigma3);
  }
  else {
    printf("SHAPE OF BACKGROUND NOT DEFINE: using EXPECTED_BACKGROUND\n");
    RooHistPdf bkg1("bkg1", "bkg1", x, hbkg1, 0);
    RooHistPdf bkg2("bkg2", "bkg2", x, hbkg2, 0);
    RooHistPdf bkg3("bkg3", "bkg3", x, hbkg3, 0);
    RooHistPdf bkg4("bkg4", "bkg4", x, hbkg4, 0);
  }
  printf("**************************************************\n");

  /*** DEFINE MISMATCH BACKGROUND SHAPE ***/

  printf("***** MISMATCH BACKGROUND SHAPE DEFINITION *****\n");
  
  /* variables and generic shapes */
  Double_t expectedCutoff = hBkgExp3->GetMean();
  Double_t expectedCutoffRMS = hBkgExp3->GetRMS();
  //    RooRealVar cutoff("cutoff", "cutoff", -30., -50., 0., "");
  RooRealVar cutoff("cutoff", "cutoff", expectedCutoff, expectedCutoff - 3. * expectedCutoffRMS, expectedCutoff, "");
  //  RooRealVar cutoff("cutoff", "cutoff", expectedCutoff, "");
  //  RooRealVar power("power", "power", 1., 0.5, 5.0, "");
  RooRealVar power("power", "power", 1., "");
  RooFermiCutoff fermi("fermi", "fermi", x, cutoff, power);
  RooRealVar alpha("alpha", "alpha", 0., -10., 0., "");
  RooExponential expo("expo", "expo", x, alpha);
  RooUniform uniform("uniform", "uniform", x);

  /* mismatch shape */
  if (EXPECTED_MISMATCH) {
    printf("USING EXPECTED MISMATCH SHAPE FROM HISTOGRAMS\n");
    RooHistPdf mismatch("mismatch", "mismatch", x, hmismatch, 0);
  }
  else if (UNIFORM_MISMATCH) {
    printf("USING UNIFORM BACKGROUND SHAPE WITH CUTOFF\n");
    RooProdPdf mismatch("mismatch", "mismatch", fermi, uniform, -100.);
  }
  else if (EXPONENTIAL_MISMATCH) {
    printf("USING EXPONENTIAL BACKGROUND SHAPE WITH CUTOFF\n");
    RooProdPdf mismatch("mismatch", "mismatch", fermi, expo, -100.);
  }
  else if (DOUBLEEXPONENTIAL_MISMATCH) {
    printf("USING DOUBLE-EXPONENTIAL BACKGROUND SHAPE WITH CUTOFF\n");
    RooRealVar alpha1("alpha1", "alpha1", 0., -10., 0., "");
    RooRealVar alpha2("alpha2", "alpha2", 0., -10., 0., "");
    RooRealVar frac("frac", "frac", 0.5, 0., 1., "");
    RooGenericPdf doubleexpo("doubleexpo", "TMath::Exp(alpha1 * x) + frac * TMath::Exp(alpha2 * x)", RooArgSet(x, alpha1, alpha2, frac));
    RooProdPdf mismatch("mismatch", "mismatch", fermi, doubleexpo, -100.);
  }
  else if (UNIFORMPLUSEXPONENTIAL_MISMATCH) {
    printf("USING UNIFORM-PLUS-EXPONENTIAL BACKGROUND SHAPE WITH CUTOFF\n");
    RooRealVar funiform("funiform", "funiform", 100., 0., 100000., "");
    RooRealVar fexpo("fexpo", "fexpo", 100., 0., 100000., "");
    RooAddPdf uniformplusexpo("uniformplusexpo", "uniformplusexpo", RooArgList(uniform, expo), RooArgList(funiform, fexpo), kFALSE);
    RooProdPdf mismatch("mismatch", "mismatch", fermi, uniformplusexpo, -100.);
  }
  else {
    printf("SHAPE OF MISMATCH NOT DEFINE: using EXPECTED_MISMATCH\n");
    RooHistPdf mismatch("mismatch", "mismatch", x, hmismatch, 0);
  }
  printf("************************************************\n");

  /*** DEFINE THE MODEL ***/

  printf("***** MODEL DEFINITION *****\n");

  /* variables */
  Double_t integral = hdata.sumEntries();
  RooRealVar nsignal("nsignal", "nsignal", 0.3 * integral, 0., integral);
  RooRealVar nbkg1("nbkg1", "nbkg1", 0.3 * integral, 0., integral);
  RooRealVar nbkg2("nbkg2", "nbkg2", 0.3 * integral, 0., integral);
  RooRealVar nbkg3("nbkg3", "nbkg3", 0.3 * integral, 0., integral);
  RooRealVar nbkg4("nbkg4", "nbkg4", 0.3 * integral, 0., integral);
  RooRealVar nmismatch("nmismatch", "nmismatch", 0.1 * integral, 0., integral);

if (!fitBkg1) {
  nbkg1.setVal(0.);
  nbkg1.setConstant(kTRUE);  
 }
if (!fitBkg2) {
  nbkg2.setVal(0.);
  nbkg2.setConstant(kTRUE);  
 }
if (!fitBkg3) {
  nbkg3.setVal(0.);
  nbkg3.setConstant(kTRUE);  
 }
if (!fitBkg4) {
  nbkg4.setVal(0.);
  nbkg4.setConstant(kTRUE);  
 }

  RooAddPdf model("model", "model p.d.f.", RooArgList(signal, bkg1, bkg2, bkg3, bkg4, mismatch), RooArgList(nsignal, nbkg1, nbkg2, nbkg3, nbkg4, nmismatch));

#if 0
  /* the model */
  if (USE_ELECTRON_BACKGROUND && fitElectrons && fitPions) {
    printf("USING ELECTRON BACKGROUND\n");
    if (NO_MISMATCH) {
      printf("NOT USING MISMATCH BACKGROUND\n");
      nmismatch.setVal(0.);
      RooAddPdf model("model", "model p.d.f.", RooArgList(signal, bkg1, bkg2, bkg3/*, bkg4*/), RooArgList(nsignal, nbkg1, nbkg2, nbkg3/*, nbkg4*/));
    }
    else {
      printf("USING MISMATCH BACKGROUND\n");
      RooAddPdf model("model", "model p.d.f.", RooArgList(signal, bkg1, bkg2, bkg3/*, bkg4*/, mismatch), RooArgList(nsignal, nbkg1, nbkg2, nbkg3/*, nbkg4*/, nmismatch));
    }
  }
  else if (!USE_ELECTRON_BACKGROUND || !fitElectrons) {
    printf("NOT USING ELECTRON BACKGROUND\n");
    nbkg3.setVal(0.);
    nbkg4.setVal(0.);
    if (NO_MISMATCH) {
      printf("NOT USING MISMATCH BACKGROUND\n");
      nmismatch.setVal(0.);
      RooAddPdf model("model", "model p.d.f.", RooArgList(signal, bkg1, bkg2), RooArgList(nsignal, nbkg1, nbkg2));
    }
    else {
      printf("USING MISMATCH BACKGROUND\n");
      RooAddPdf model("model", "model p.d.f.", RooArgList(signal, bkg1, bkg2, mismatch), RooArgList(nsignal, nbkg1, nbkg2, nmismatch));
    }
  }
  printf("****************************\n");
#endif




  /*** FIT ***/

  printf("***** FIT *****\n");

  printf("SIGNAL PARAMETERS BEFORE GLOBAL FIT\n");
  printf("mean  = %f +- %f\n", mean.getVal(), mean.getError());
  printf("sigma = %f +- %f\n", sigma.getVal(), sigma.getError());
  printf("tail  = %f +- %f\n", tail.getVal(), tail.getError());

  /* fit and draw */
  if (canvas) canvas->cd();
  //  model.chi2FitTo(hdata, Extended(kTRUE), Verbose(kFALSE), SumW2Error(kFALSE), Range(-40., 140.), Binned(kTRUE));
//  model.fitTo(hdata, Extended(kTRUE), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10), Range(-10., 10.));
model.fitTo(hdata, Range(rangeMin, rangeMax), Extended(kTRUE), SumW2Error(kFALSE), Verbose(kFALSE), PrintEvalErrors(10));

  printf("***************\n");

  /*** DRAW ***/
#if DISPLAY
  RooPlot *xframe = x.frame();
hdata.plotOn(xframe, XErrorSize(0), DrawOption("PZ"));
model.plotOn(xframe, LineWidth(2)/*, Precision(1.e-4)*/);
model.plotOn(xframe, Components(signal), LineWidth(2), LineColor(kRed)/*, Precision(1.e-4)*/);
  model.plotOn(xframe, Components(bkg1), LineWidth(2), LineStyle(kDashed), LineColor(kCyan+1));
  model.plotOn(xframe, Components(bkg2), LineWidth(2), LineStyle(kDashed), LineColor(kGreen+1));
  if (USE_ELECTRON_BACKGROUND) {
    model.plotOn(xframe, Components(bkg3), LineWidth(2), LineStyle(kDashed), LineColor(kYellow+1));
    model.plotOn(xframe, Components(bkg4), LineWidth(2), LineStyle(kDashed), LineColor(kYellow+2));
  }
  if (!NO_MISMATCH)
    model.plotOn(xframe, Components(mismatch), LineWidth(2), LineStyle(kDashed), LineColor(kGray+1));
  hSignal->SetFillColor(kYellow);
  hSignal->SetLineWidth(0);
  hSignal->SetFillStyle(0);
  hSignal->SetMinimum(0.1);
  hSignal->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
//  hSignal->Draw();
xframe->SetMinimum(0.1);
  xframe->Draw();
#endif
  if (canvas) canvas->Update();

  /*** COMPUTE CHI2 ***/
  Double_t datax, datapoint, datapoint_err, modelpoint;
  Double_t chi2 = 0.;
  Int_t ndf = 0;
  for (Int_t ibin = 0; ibin < hSignal->GetNbinsX(); ibin++) {
    datax = hSignal->GetBinCenter(ibin + 1);
    if (datax < rangeMin || datax > rangeMax) continue;
    datapoint = hSignal->GetBinContent(ibin + 1);
    datapoint_err = hSignal->GetBinError(ibin + 1);
    if (datapoint_err == 0.) continue;
    x.setVal(datax);
    modelpoint = model.getVal();
    chi2 += (datapoint - modelpoint) * (datapoint - modelpoint) / (datapoint_err * datapoint_err);
    ndf++;
  }


/*** PRINT FIT OUTPUT ***/
printf("***** CHI-SQUARE *****\n");
  printf("chi-square = %f\n", chi2);
  printf("NDF        = %d\n", ndf);
  printf("chi2/NDF   = %f\n", ndf > 0 ? chi2 / ndf : 0.);
  printf("***** SIGNAL-SHAPE INFO *****\n");
  printf("mean      = %f +- %f\n", mean.getVal(), mean.getError());
  printf("sigma     = %f +- %f\n", sigma.getVal(), sigma.getError());
  printf("tail      = %f +- %f\n", tail.getVal(), tail.getError());
  printf("pure gaus = %f +- %f\n", gaussianfrac.getVal(), gaussianfrac.getError());
printf("*****************************\n");
printf("***** COUNTS *****\n");
printf("total     = %f\n", hSignal->Integral(-1, -1));
  printf("integral  = %f\n", hdata.sumEntries());
  printf("nsignal   = %f +- %f\n", nsignal.getVal(), nsignal.getError());
  printf("nbkg1     = %f +- %f\n", nbkg1.getVal(), nbkg1.getError());
  printf("nbkg2     = %f +- %f\n", nbkg2.getVal(), nbkg2.getError());
  printf("nbkg3     = %f +- %f\n", nbkg3.getVal(), nbkg3.getError());
  printf("nbkg4     = %f +- %f\n", nbkg4.getVal(), nbkg4.getError());
  printf("nmismatch = %f +- %f\n", nmismatch.getVal(), nmismatch.getError());
  printf("******************\n");

  /*** OUTPUT FIT PARAMS ***/
  
  param[kMean] = mean.getVal();
  param_err[kMean] = mean.getError();
  param[kSigma] = sigma.getVal();
  param_err[kSigma] = sigma.getError();
  param[kTail] = tail.getVal();
  param_err[kTail] = tail.getError();
  param[kTotalCounts] = hSignal->GetEntries();
  param_err[kTotalCounts] = sqrt(hSignal->GetEntries());
  param[kIntegralCounts] = hdata.sumEntries();
  param_err[kIntegralCounts] = sqrt(hdata.sumEntries());
  param[kSignalCounts] = nsignal.getVal();
  param_err[kSignalCounts] = nsignal.getError();
  param[kBkg1Counts] = nbkg1.getVal();
  param_err[kBkg1Counts] = nbkg1.getError();
  param[kBkg2Counts] = nbkg2.getVal();
  param_err[kBkg2Counts] = nbkg2.getError();
  param[kBkg3Counts] = nbkg3.getVal();
  param_err[kBkg3Counts] = nbkg3.getError();
  param[kBkg4Counts] = nbkg4.getVal();
  param_err[kBkg4Counts] = nbkg4.getError();
  param[kMismatchCounts] = nmismatch.getVal();
  param_err[kMismatchCounts] = nmismatch.getError();

  return;
}

//___________________________________________________________________________________

TOFpid_efficiency(Int_t ipart, Int_t icharge, Int_t icent, Int_t useTPCcut = -1, Double_t minsignalFrac = 0., Int_t nTrials = 1000)
{

  /* prepare centrality name */
  Char_t centName[1024];
  if (icent < 0 || icent >= NcentralityBins)
    sprintf(centName, "cent0090");
  else
    sprintf(centName, "cent%02d%02d", centralityBin[icent], centralityBin[icent + 1]);

  /* fraction names */
  const Char_t *fractionName[AliPID::kSPECIES][AliPID::kSPECIES] = {
    "hSignalFraction", "hBkg4Fraction", "hBkg3Fraction", "hBkg1Fraction", "hBkg2Fraction",
    "hBkg4Fraction", "hSignalFraction", "hBkg3Fraction", "hBkg1Fraction", "hBkg2Fraction",
    "hBkg3Fraction", "hBkg4Fraction", "hSignalFraction", "hBkg1Fraction", "hBkg2Fraction",
    "hBkg3Fraction", "hBkg4Fraction", "hBkg1Fraction", "hSignalFraction", "hBkg2Fraction",
    "hBkg3Fraction", "hBkg4Fraction", "hBkg1Fraction", "hBkg2Fraction", "hSignalFraction"
  };

  enum ETPCcut_t {
    kCurrentDir,
    k11Cut,
    k12Cut,
    k21Cut,
    k22Cut,
    k23Cut,
    k32Cut,
    k33Cut,
    kNTPCcuts
  };
  const Char_t *tpccutdir[8] = {
    ".",
    "TOFpid_cutOnTPC[-1.0,1.0]",
    "TOFpid_cutOnTPC[1.0,2.0]",
    "TOFpid_cutOnTPC[-2.0,-1.0]",
    "TOFpid_cutOnTPC[-2.0,2.0]",
    "TOFpid_cutOnTPC[2.0,3.0]",
    "TOFpid_cutOnTPC[-3.0,-2.0]",
    "TOFpid_cutOnTPC[-3.0,3.0]"
  };

  /* get data */
  Char_t filename[1024];
  TH1D *hAccepted[AliPID::kSPECIES][kNTPCcuts];
  TH1D *hIdentified[AliPID::kSPECIES][kNTPCcuts];
  TH1D *hEfficiencyIn[AliPID::kSPECIES][kNTPCcuts];
  TH1D *hFraction[AliPID::kSPECIES][AliPID::kSPECIES][kNTPCcuts];
  for (Int_t itpccut = 0; itpccut < kNTPCcuts; itpccut++) {

    /* check whether we use this cut */
    if (useTPCcut >= 0 && itpccut != useTPCcut) 
      continue;

    for (Int_t iipart = 0; iipart < AliPID::kSPECIES; iipart++) {
      /* skip electrons and muons */
      if (iipart == AliPID::kMuon)
	continue;
      
      /* open file */
      sprintf(filename, "%s/TOFspectrum_%s_%s_%s_%sID.root", tpccutdir[itpccut], centName, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart));
      TFile *filein = TFile::Open(filename);
      if (!filein || !filein->IsOpen()) {
	printf("cannot open %s\n", filename);
	return;
      }
      
      /* get accepted tracks */
      hAccepted[iipart][itpccut] = (TH1D *)filein->Get("PostAnalysis/hAcceptedTracks");
      if (!hAccepted[iipart][itpccut]) {
	printf("cannot find PostAnalysis/hAcceptedTracks in %s\n", filename);
	return;
      }
      
      /* get identified tracks */
      hIdentified[iipart][itpccut] = (TH1D *)filein->Get("PostAnalysis/hIdentifiedCounts");
      if (!hIdentified[iipart][itpccut]) {
	printf("cannot find PostAnalysis/hIdentifiedCounts in %s\n", filename);
	return;
      }
      
      /* compute efficiency */
      hEfficiencyIn[iipart][itpccut] = new TH1D(*hIdentified[iipart][itpccut]);
      hEfficiencyIn[iipart][itpccut]->Divide(hEfficiencyIn[iipart][itpccut], hAccepted[iipart][itpccut], 1., 1., "B");
      
      /* get fractions */
      for (Int_t iiipart = 0; iiipart < AliPID::kSPECIES; iiipart++) {
	/* skip muons */
	if (iiipart == AliPID::kMuon)
	  continue;
	
	hFraction[iipart][iiipart][itpccut] = (TH1D *)filein->Get(Form("PostAnalysis/%s", fractionName[iipart][iiipart]));
	if (!hFraction[iipart][iiipart][itpccut]) {
	  printf("cannot find PostAnalysis/%s in %s\n", fractionName[iipart][iiipart], filename);
	  return;
	}
      }
    }
  }

  /* prepare output efficiency histos */
  TH1D *hEfficiencyOut[AliPID::kSPECIES];
  for (Int_t iipart = 0; iipart < AliPID::kSPECIES; iipart++) {
    hEfficiencyOut[iipart] = new TH1D(Form("hPIDEff_%d_%s_%s_%sID", icent, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(iipart)), "", NptBins, ptBin);
  }
 

  /*** MATRIX ***/

  const Int_t nPart = 4;
  Int_t partIndex[4] = {AliPID::kElectron, AliPID::kPion, AliPID::kKaon, AliPID::kProton};

  TMatrixD Meffin(nPart, 1);
  TMatrixD Meffout(nPart, 1);
  TMatrixD Mfrac(nPart, nPart);

  Double_t eff[4], effe[4], frac[4][4], frace[4][4], geneff, genfrac;
  Bool_t gotbestcut[4];

  TH1F *hEffTemp[4];
  for (Int_t iipart = 0; iipart < nPart; iipart++)
    hEffTemp[iipart] = new TH1F(Form("hEffTemp_%d", iipart), "", 20000, -10., 10.);

  /* loop over pt-bins */
  for (Int_t ibin = 0; ibin < NptBins; ibin++) {

    /* reset temp histos */
    for (Int_t iipart = 0; iipart < nPart; iipart++)
      hEffTemp[iipart]->Reset();

    /* get measured data */
    for (Int_t iipart = 0; iipart < nPart; iipart++) {

      /* select the best set of cuts */
      Double_t signalFrac, maxsignalFrac = minsignalFrac;
      Int_t bestcut = -1;
      gotbestcut[iipart] = kFALSE;
      for (Int_t itpccut = 0; itpccut < kNTPCcuts; itpccut++) {

	/* check whether we use this cut */
	if (useTPCcut >= 0 && itpccut != useTPCcut) 
	  continue;
	
	/* maximize the signal fraction */
	signalFrac = hFraction[partIndex[iipart]][partIndex[iipart]][itpccut]->GetBinContent(ibin + 1);
	if (signalFrac > maxsignalFrac) {
	  maxsignalFrac = signalFrac;
	  bestcut = itpccut;
	  gotbestcut[iipart] = kTRUE;
	}
      }
      if (bestcut < 0)
	continue;

      eff[iipart] = hEfficiencyIn[partIndex[iipart]][bestcut]->GetBinContent(ibin + 1);
      effe[iipart] = hEfficiencyIn[partIndex[iipart]][bestcut]->GetBinError(ibin + 1);
      for (Int_t iiipart = 0; iiipart < nPart; iiipart++) {
	frac[iipart][iiipart] = hFraction[partIndex[iipart]][partIndex[iiipart]][bestcut]->GetBinContent(ibin + 1);
	frace[iipart][iiipart] = hFraction[partIndex[iipart]][partIndex[iiipart]][bestcut]->GetBinError(ibin + 1);
      }
    }

    /* check best cuts */
    Bool_t skip = kFALSE;
    for (Int_t iipart = 0; iipart < nPart; iipart++)
      if (!gotbestcut[iipart])
	skip = kTRUE;
    if (skip) continue;
    
    /* loop over trials */
    for (Int_t itry = 0; itry < nTrials; itry++) {
      
      /* setup matrix */
      for (Int_t iipart = 0; iipart < nPart; iipart++) {
	geneff = gRandom->Gaus(eff[iipart], effe[iipart]);
	if (geneff < 0.) geneff = 0.;
	if (geneff > 1.) geneff = 1.;
	Meffin[iipart] = geneff != 0. ? 1. / geneff : 0.;
	for (Int_t iiipart = 0; iiipart < nPart; iiipart++) {
	  genfrac = gRandom->Gaus(frac[iipart][iiipart], frace[iipart][iiipart]);
	  if (genfrac < 0.) genfrac = 0.;
	  if (genfrac > 1.) genfrac = 1.;
	  Mfrac[iipart][iiipart] = genfrac;
	}
      }

      /* invert matrix */
      TDecompLU lu(Mfrac);
      TMatrixD Minv(nPart, nPart);
      if (!lu.Decompose()) continue;
      lu.Invert(Minv);
      Meffout = Minv * Meffin;
      
      /* fill histos */
      for (Int_t iipart = 0; iipart < nPart; iipart++) {
	if (Meffout.GetMatrixArray()[iipart] > 0.)
	  hEffTemp[iipart]->Fill(1. / Meffout.GetMatrixArray()[iipart]);
      }

    } /* end of loop over trials */

    
    /* get average and RMS */
    for (Int_t iipart = 0; iipart < nPart; iipart++) {
      hEfficiencyOut[partIndex[iipart]]->SetBinContent(ibin + 1, hEffTemp[iipart]->GetMean());
      hEfficiencyOut[partIndex[iipart]]->SetBinError(ibin + 1, hEffTemp[iipart]->GetRMS());
    }

  } /* end of loop over pt-bins */

  
  hEfficiencyOut[AliPID::kPion]->SetMarkerStyle(20);
  hEfficiencyOut[AliPID::kKaon]->SetMarkerStyle(20);
  hEfficiencyOut[AliPID::kProton]->SetMarkerStyle(20);
  hEfficiencyOut[AliPID::kElectron]->SetMarkerStyle(20);

  hEfficiencyOut[AliPID::kPion]->SetMarkerColor(4);
  hEfficiencyOut[AliPID::kKaon]->SetMarkerColor(8);
  hEfficiencyOut[AliPID::kProton]->SetMarkerColor(2);
  hEfficiencyOut[AliPID::kElectron]->SetMarkerColor(1);
 
  hEfficiencyOut[AliPID::kPion]->Draw();
  hEfficiencyOut[AliPID::kKaon]->Draw("same");
  hEfficiencyOut[AliPID::kProton]->Draw("same");
  hEfficiencyOut[AliPID::kElectron]->Draw("same");


  /* output */
  TFile *fileout = TFile::Open(Form("TOFpid_efficiency_cent%d_%s_%s.root", icent, AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");
  hEfficiencyOut[AliPID::kPion]->Write();
  hEfficiencyOut[AliPID::kKaon]->Write();
  hEfficiencyOut[AliPID::kProton]->Write();
  hEfficiencyOut[AliPID::kElectron]->Write();
  fileout->Close();

}

TOFpid_normalize(TH1D *h, Double_t nevents = 8.42446600000000000e+06)
{
  
  Double_t counts, counts_err;
  for (Int_t ibin = 0; ibin < h->GetNbinsX(); ibin++) {
    counts = h->GetBinContent(ibin + 1);
    counts_err = h->GetBinError(ibin + 1);
    counts /= h->GetBinWidth(ibin + 1);
    counts_err /= h->GetBinWidth(ibin + 1);
    h->SetBinContent(ibin + 1, counts);
    h->SetBinError(ibin + 1, counts_err);
  }
  
  h->Scale(1. / nevents);
  
}

TOFpid_normalizeAndwrite(const Char_t *fileoutname, const Char_t *corrstring = "")
{

  TFile *fpiplus = TFile::Open("TOFpid_spectrum_pion_plus.root");
  TFile *fpiminus = TFile::Open("TOFpid_spectrum_pion_minus.root");
  TFile *fkaplus = TFile::Open("TOFpid_spectrum_kaon_plus.root");
  TFile *fkaminus = TFile::Open("TOFpid_spectrum_kaon_minus.root");
  TFile *fprplus = TFile::Open("TOFpid_spectrum_proton_plus.root");
  TFile *fprminus = TFile::Open("TOFpid_spectrum_proton_minus.root");

  hpiplus = (TH1F *)fpiplus->Get(Form("hSignal%sCounts", corrstring));
  hpiminus = (TH1F *)fpiminus->Get(Form("hSignal%sCounts", corrstring));
  hkaplus = (TH1F *)fkaplus->Get(Form("hSignal%sCounts", corrstring));
  hkaminus = (TH1F *)fkaminus->Get(Form("hSignal%sCounts", corrstring));
  hprplus = (TH1F *)fprplus->Get(Form("hSignal%sCounts", corrstring));
  hprminus = (TH1F *)fprminus->Get(Form("hSignal%sCounts", corrstring));

  hpiplus->SetName("hpiplus");
  hpiminus->SetName("hpiminus");
  hkaplus->SetName("hkaplus");
  hkaminus->SetName("hkaminus");
  hprplus->SetName("hprplus");
  hprminus->SetName("hprminus");

  TOFpid_normalize(hpiplus);
  TOFpid_normalize(hpiminus);
  TOFpid_normalize(hkaplus);
  TOFpid_normalize(hkaminus);
  TOFpid_normalize(hprplus);
  TOFpid_normalize(hprminus);

  TFile *fileout = TFile::Open(fileoutname, "RECREATE");
  hpiplus->Write("hpiplus");
  hpiminus->Write("hpiminus");
  hkaplus->Write("hkaplus");
  hkaminus->Write("hkaminus");
  hprplus->Write("hprplus");
  hprminus->Write("hprminus");
  fileout->Close();

}

/**************************************************************/

TOFpid_rawSpectra(const Char_t *destdir = "default")
{

  Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };

  TFile *fileout = TFile::Open("TOF_rawSpectra.root", "RECREATE");
  TH1D *hRaw;
  Char_t title[1024];
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	hRaw = TOFpid_rawSpectra(destdir, ipart, icharge, icent);
	if (!hRaw) continue;
	hRaw->SetMarkerStyle(marker[icharge]);
	hRaw->SetMarkerColor(color[ipart]);
	hRaw->SetLineColor(1);
	hRaw->SetLineWidth(1);
	sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]);
	hRaw->SetTitle(title);
	hRaw->SetName(Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	fileout->cd();
	hRaw->Write();
      }

  fileout->Close();
}

TH1D *
TOFpid_rawSpectra(const Char_t *dirname, Int_t ipart, Int_t icharge, Int_t icent)
{

  /* open data */
  Char_t outfilename[1024];
  if (icent < 0 || icent >= NcentralityBins) {
    sprintf(outfilename, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", dirname, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  else {
    sprintf(outfilename, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", dirname, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  TFile *filein = TFile::Open(outfilename);
  if (!filein || !filein->IsOpen()) {
    printf("cannot open %s\n", outfilename);
    return NULL;
  }
  /* get data */
  TH1D *h = (TH1D *)filein->Get("RawSpectra/hNormalizedRawYield");
  if (!h) {
    printf("cannot get RawSpectra/hNormalizedRawYield from %s\n", outfilename);
    return NULL;
  }

  return h;

}

TOFpid_rawSpectra_mismatchCorrected(const Char_t *destdir = "default")
{

  Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };

  TFile *fileout = TFile::Open("TOF_rawSpectra_mismatchCorrected.root", "RECREATE");
  TH1D *hRaw;
  Char_t title[1024];
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	hRaw = TOFpid_rawSpectra_mismatchCorrected(destdir, ipart, icharge, icent);
	if (!hRaw) continue;
	hRaw->SetMarkerStyle(marker[icharge]);
	hRaw->SetMarkerColor(color[ipart]);
	hRaw->SetLineColor(1);
	hRaw->SetLineWidth(1);
	sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]);
	hRaw->SetTitle(title);
	hRaw->SetName(Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	fileout->cd();
	hRaw->Write();
      }

  fileout->Close();
}

TH1D *
TOFpid_rawSpectra_mismatchCorrected(const Char_t *dirname, Int_t ipart, Int_t icharge, Int_t icent)
{

  /* open data */
  Char_t outfilename[1024];
  if (icent < 0 || icent >= NcentralityBins) {
    sprintf(outfilename, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", dirname, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  else {
    sprintf(outfilename, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", dirname, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  TFile *filein = TFile::Open(outfilename);
  if (!filein || !filein->IsOpen()) {
    printf("cannot open %s\n", outfilename);
    return NULL;
  }
  /* get data */
  TH1D *h = (TH1D *)filein->Get("RawSpectra/hNormalizedMismatchCorrectedRawYield");
  if (!h) {
    printf("cannot get RawSpectra/hNormalizedRawYield from %s\n", outfilename);
    return NULL;
  }

  return h;

}

/**************************************************************/

TOFpid_electronCorrection()
{

    Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };

  
  TFile *fileout = TFile::Open("TOF_electronCorrection.root", "RECREATE");
  TH1D *hCorr;
  TProfile *pCorrAv = new TProfile("pCorrAv", "", NptBins, ptBin, "s");
  Char_t title[1024];
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < 3; ipart++) {
	hCorr = TOFpid_electronCorrection(ipart, icharge, icent);
	hCorr->SetMarkerStyle(marker[icharge]);
	hCorr->SetMarkerColor(color[ipart]);
	hCorr->SetLineColor(1);
	hCorr->SetLineWidth(1);
	sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);electron-subtracted pion fraction;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]);
	hCorr->SetTitle(title);
	hCorr->SetName(Form("hElectronCorr_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	fileout->cd();
	hCorr->Write();
	/* fill for average correction */
	for (Int_t ipt = 0; ipt < NptBins; ipt++) {
	  pCorrAv->Fill(hCorr->GetBinCenter(ipt + 1), hCorr->GetBinContent(ipt + 1));
	}
      }
  hCorr = new TH1D("hElectronCorr_average", ";p_{T} (GeV/c);electron-subtracted pion fraction;", NptBins, ptBin);
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    hCorr->SetBinContent(ipt + 1, pCorrAv->GetBinContent(ipt + 1));
    hCorr->SetBinError(ipt + 1, pCorrAv->GetBinError(ipt + 1));
  }
  hCorr->Write();
  fileout->Close();
}

TH1D *
TOFpid_electronCorrection(Int_t ipart, Int_t icharge, Int_t icent)
{

  Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };

  TH1D *hr = TOFpid_systematics_ratio("defaultFit_electronCut", "defaultFit", ipart, icharge, icent, "electronCorrection", "", marker[icharge], color[ipart], kTRUE);
  TF1 fOne("fOne", "1.", fitPtMin[ipart], fitPtMax[ipart]);
  hr->Add(&fOne, -1.);
  hr->Scale(1. / 0.866);
  hr->Add(&fOne, 1.);

  return hr;
}

/**************************************************************/

SystematicCheck(const Char_t *defaultname = "defaultFit")
{

  const Int_t ndata = 6;
  const Char_t *name[ndata] = {
    "signalFit_fixed_scaleSigma_09",
    "signalFit_fixed_scaleSigma_11",
    "signalFit_fixed_scaleTail_09",
    "signalFit_fixed_scaleTail_11",
    "signalFit_fixed_scaleSigma_09_scaleTail_11",
    "signalFit_fixed_scaleSigma_11_scaleTail_09"
  };
  for (Int_t idata = 0; idata < ndata; idata++)
    SystematicCheck(name[idata], defaultname);

}

SystematicCheck(const Char_t *checkname, const Char_t *defaultname = "defaultFit")
{

  gROOT->LoadMacro("HistoUtils.C");

  Char_t filename1[1024];
  Char_t filename2[1024];

  Int_t marker[2] = {20, 25};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };


  TFile *fileout = TFile::Open(Form("SystematicCheck_%s.root", checkname), "RECREATE");
  Char_t title[1024];
  TH1 *hd;
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {

	if (icent < 0 || icent >= NcentralityBins) {
	  sprintf(filename1, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", checkname, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
	  sprintf(filename2, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", defaultname, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
	}
	else {
	  sprintf(filename1, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", checkname, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
	  sprintf(filename2, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", defaultname, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
	}
	
	hd = HistoUtils_systematics(filename1, filename2, "RawSpectra/hNormalizedRawYield", NULL);
	hd->SetName(Form("hRaw_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);#frac{d^{2}N}{dy dp_{T}} (c/GeV);", partLatex[ipart][icharge], centralityBin[icent], centralityBin[icent + 1]);
	hd->SetTitle(title);
	hd->SetLineColor(1);
	hd->SetLineWidth(1);
	hd->SetMarkerStyle(marker[icharge]);
	hd->SetMarkerColor(color[ipart]);
	fileout->cd();
	hd->Write();
	delete hd;
      }
  fileout->Close();

}

/**************************************************************/

TOFpid_systematics()
{

  TFile *fileout = TFile::Open("TOFpid_systematics.root", "RECREATE");
  TH1D *hSys;
  for (Int_t icent = 0; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	hSys = TOFpid_systematics(ipart, icharge, icent);
	fileout->cd();
	hSys->Write(Form("hRawSys_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	delete hSys;
      }

  fileout->Close();
}

TH1D *
TOFpid_systematics(Int_t ipart, Int_t icharge, Int_t icent)
{
  
  TH1D *hSignalFit = TOFpid_systematics_signalFit(ipart, icharge, icent);
  TH1D *hBkgFit = TOFpid_systematics_bkgFit(ipart, icharge, icent);
  TH1D *hFitRange = TOFpid_systematics_fitRange(ipart, icharge, icent);

  TH1D *hSys = new TH1D("hSys", "", NptBins, ptBin);
  Double_t sigsys, bkgsys, rangesys, totsys = 0.;
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    sigsys = hSignalFit->GetBinError(ipt + 1);
    bkgsys = hBkgFit->GetBinError(ipt + 1);
    rangesys = hFitRange->GetBinError(ipt + 1);
    totsys = TMath::Sqrt(sigsys * sigsys + bkgsys * bkgsys + rangesys * rangesys);
    hSys->SetBinContent(ipt + 1, totsys);
  }

  hSys->Draw();

  delete hSignalFit;
  delete hBkgFit;
  delete hFitRange;

  return hSys;
}

TH1D *
TOFpid_systematics_fitRange(Int_t ipart, Int_t icharge, Int_t icent)
{
  
  TH1D *hArea = new TH1D("hArea", "", NptBins, ptBin);
  hArea->SetMinimum(0.8);
  hArea->SetMaximum(1.2);
  hArea->Draw();

  TH1D *hr = TOFpid_systematics_ratio("default_wideRange", "default", ipart, icharge, icent, "wide", "wide-range fit;p_{T} (GeV/c);ratio wrt. default;", 20, 4);
  hr->Draw("same");

  TH1D *hSys = new TH1D("hSys", "", NptBins, ptBin);
  hSys->SetFillStyle(0);
  hSys->SetMarkerSize(0);
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    Double_t val, max = 0.;
    if (hr->GetBinContent(ipt + 1) == 0.) continue;
    val = TMath::Abs(hr->GetBinContent(ipt + 1) - 1.);
    hSys->SetBinContent(ipt + 1, 1.);
    hSys->SetBinError(ipt + 1, val);
  }
  hSys->Draw("same, E2");
  
  delete hr;

  return hSys;
}

TH1D *
TOFpid_systematics_signalFit(Int_t ipart, Int_t icharge, Int_t icent)
{

  const Int_t ndata = 6;
  const Char_t *name[ndata] = {
    "signalFit_fixed_scaleSigma_09",
    "signalFit_fixed_scaleSigma_11",
    "signalFit_fixed_scaleTail_09",
    "signalFit_fixed_scaleTail_11",
    "signalFit_fixed_scaleSigma_09_scaleTail_11",
    "signalFit_fixed_scaleSigma_11_scaleTail_09"
  };
  const Char_t *title[ndata] = {
    "-10% #sigma;p_{T} (GeV/c); ratio",
    "+10% #sigma;p_{T} (GeV/c); ratio",
    "-10% #tau;p_{T} (GeV/c); ratio",
    "+10% #tau;p_{T} (GeV/c); ratio",
    "-10% #sigma, +10% #tau;p_{T} (GeV/c); ratio",
    "+10% #sigma, -10% #tau;p_{T} (GeV/c); ratio"
  };
  Int_t marker[ndata] = {22, 28, 22, 28, 22, 28};
  Int_t color[ndata] = {2, 2, 8, 8, 4, 4};

  TH1D *hArea = new TH1D("hArea", "", NptBins, ptBin);
  hArea->SetMinimum(0.8);
  hArea->SetMaximum(1.2);
  hArea->Draw();
  
  TH1D *hr[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) {
    hr[idata] = TOFpid_systematics_ratio(name[idata], "defaultFit", ipart, icharge, icent, name[idata], title[idata], marker[idata], color[idata]);
    hr[idata]->Draw("same");
  }

  TH1D *hSys = new TH1D("hSys", "", NptBins, ptBin);
  hSys->SetFillStyle(0);
  hSys->SetMarkerSize(0);
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    Double_t val, max = 0.;
    for (Int_t idata = 0; idata < ndata; idata++) {
      if (hr[idata]->GetBinContent(ipt + 1) == 0.) continue;
      val = TMath::Abs(hr[idata]->GetBinContent(ipt + 1) - 1.);
      if (val > 0.9) continue;
      if (val > max)
	max = val;
    }
    hSys->SetBinContent(ipt + 1, 1.);
    hSys->SetBinError(ipt + 1, max);
  }
  hSys->Draw("same, E2");

  //  delete hArea;
  //  for (Int_t idata = 0; idata < ndata; idata++)
  //    delete hr[idata];
  
  return hSys;
}

TH1D *
TOFpid_systematics_bkgFit(Int_t ipart, Int_t icharge, Int_t icent)
{

  const Int_t ndata = 6;
  const Char_t *name[ndata] = {
    "bkgFit_fixed_scaleSigma_09",
    "bkgFit_fixed_scaleSigma_11",
    "bkgFit_fixed_scaleTail_09",
    "bkgFit_fixed_scaleTail_11",
    "bkgFit_fixed_scaleSigma_09_scaleTail_11",
    "bkgFit_fixed_scaleSigma_11_scaleTail_09"
  };
  const Char_t *title[ndata] = {
    "-10% #sigma;p_{T} (GeV/c); ratio",
    "+10% #sigma;p_{T} (GeV/c); ratio",
    "-10% #tau;p_{T} (GeV/c); ratio",
    "+10% #tau;p_{T} (GeV/c); ratio",
    "-10% #sigma, +10% #tau;p_{T} (GeV/c); ratio",
    "+10% #sigma, -10% #tau;p_{T} (GeV/c); ratio"
  };
  Int_t marker[ndata] = {22, 28, 22, 28, 22, 28};
  Int_t color[ndata] = {2, 2, 8, 8, 4, 4};

  TH1D *hArea = new TH1D("hArea", "", NptBins, ptBin);
  hArea->SetMinimum(0.5);
  hArea->SetMaximum(1.5);
  hArea->Draw();
  
  TH1D *hr[ndata];
  for (Int_t idata = 0; idata < ndata; idata++) {
    hr[idata] = TOFpid_systematics_ratio(name[idata], "defaultFit", ipart, icharge, icent, name[idata], title[idata], marker[idata], color[idata]);
    hr[idata]->Draw("same");
  }

  TH1D *hSys = new TH1D("hSys", "", NptBins, ptBin);
  hSys->SetFillStyle(0);
  hSys->SetMarkerSize(0);
  for (Int_t ipt = 0; ipt < NptBins; ipt++) {
    Double_t val, max = 0.;
    for (Int_t idata = 0; idata < ndata; idata++) {
      if (hr[idata]->GetBinContent(ipt + 1) == 0.) continue;
      val = TMath::Abs(hr[idata]->GetBinContent(ipt + 1) - 1.);
      if (val > 0.9) continue;
      if (val > max)
	max = val;
    }
    hSys->SetBinContent(ipt + 1, 1.);
    hSys->SetBinError(ipt + 1, max);
  }
  hSys->Draw("same, E2");
  
  //  delete hArea;
  //  for (Int_t idata = 0; idata < ndata; idata++)
  //    delete hr[idata];

  return hSys;
}

TH1D *
TOFpid_systematics_ratio(const Char_t *dirname1, const Char_t *dirname2, Int_t ipart, Int_t icharge, Int_t icent, const Char_t *name = "rawRatio", const Char_t *title = ";p_{T} (GeV/c);raw yield ratio;", Int_t marker = 20, Int_t color = 2, Bool_t correlated = kFALSE)
{
  
  TH1D *hr = new TH1D("hr", "", NptBins, ptBin);

  /* open data */
  Char_t outfilename1[1024];
  Char_t outfilename2[1024];
  if (icent < 0 || icent >= NcentralityBins) {
    sprintf(outfilename1, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", dirname1, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
    sprintf(outfilename2, "%s/TOFspectrum_cent0090_%s_%s_%sID.root", dirname2, AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  else {
    sprintf(outfilename1, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", dirname1, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
    sprintf(outfilename2, "%s/TOFspectrum_cent%02d%02d_%s_%s_%sID.root", dirname2, centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge], AliPID::ParticleName(ipart));
  }
  TFile *filein1 = TFile::Open(outfilename1);
  if (!filein1 || !filein1->IsOpen()) {
    printf("cannot open %s\n", outfilename1);
    return;
  }
  TFile *filein2 = TFile::Open(outfilename2);
  if (!filein2 || !filein2->IsOpen()) {
    printf("cannot open %s\n", outfilename2);
    return;
  }
  /* get data */
  TH1D *h1 = (TH1D *)filein1->Get("RawSpectra/hNormalizedRawYield");
  if (!h1) {
    printf("cannot get RawSpectra/hNormalizedRawYield from %s\n", outfilename1);
    return;
  }
  TH1D *h2 = (TH1D *)filein2->Get("RawSpectra/hNormalizedRawYield");
  if (!h2) {
    printf("cannot get RawSpectra/hNormalizedRawYield from %s\n", outfilename2);
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

