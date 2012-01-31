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
  kPt,
  kTPCsignal,
  kTOFsignal,
  kTPCTOFsignal,
  kNHistoParams
};
/**************************************************************/
const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
/**************************************************************/
const Int_t NptBins = 46;
Double_t ptBin[NptBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
/**************************************************************/
const Int_t NsigmaBins = 200;
Double_t sigmaBin[NsigmaBins + 1];
Double_t sigmaMin = -10., sigmaMax = 10., sigmaStep = (sigmaMax - sigmaMin) / NsigmaBins;
/**************************************************************/
Int_t NparamsBins[kNHistoParams] = {NcentralityBins, NptBins, NsigmaBins, NsigmaBins, NsigmaBins};
Double_t *paramBin[kNHistoParams] = {centralityBin, ptBin, sigmaBin, sigmaBin, sigmaBin};
/**************************************************************/

/**************************************************************/


TPCTOFpid(const Char_t *filename, Int_t evMax = kMaxInt, Int_t startEv = 0)
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
  for (Int_t ibin = 0; ibin < NsigmaBins + 1; ibin++)
    sigmaBin[ibin] = sigmaMin + ibin * sigmaStep;

  /* THnSparse */
  THnSparse *hTPCTOFpid[AliPID::kSPECIES][kNCharges];
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge< kNCharges; icharge++) {
      hTPCTOFpid[ipart][icharge] = new THnSparseF(Form("hTPCTOFpid_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", kNHistoParams, NparamsBins);
      for (Int_t iparam = 0; iparam < kNHistoParams; iparam++) {
	hTPCTOFpid[ipart][icharge]->SetBinEdges(iparam, paramBin[iparam]);
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
      param[kPt] = analysisTrack->GetPt();

      /* get TPC info */
      dedx = analysisTrack->GetTPCdEdx();
      dedxN = analysisTrack->GetTPCdEdxN();
      ptpc = analysisTrack->GetTPCmomentum();
      
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
	
	/* TPC signal */
	bethe = tpcResponse->GetExpectedSignal(ptpc, ipart);
	deltadedx = dedx - bethe;
	dedx_sigma = tpcResponse->GetExpectedSigma(ptpc, dedxN, ipart);
	tpcsignal = deltadedx / dedx_sigma;
	param[kTPCsignal] = tpcsignal;
	
	/* TOF expected time */
	texp = analysisTrack->GetTOFExpTime(ipart);
	texp_sigma = analysisTrack->GetTOFExpTimeSigma(ipart);
	
	/* TOF signal */
	deltat = tof - texp;
	deltat_sigma = TMath::Sqrt(tof_sigma * tof_sigma + texp_sigma * texp_sigma);
	tofsignal = deltat / deltat_sigma;
	param[kTOFsignal] = tofsignal;

	/* TPC+TOF signal */
	tpctofsignal = TMath::Sqrt(tpcsignal * tpcsignal + tofsignal * tofsignal);
	param[kTPCTOFsignal] = tpctofsignal;

	/* fill histo */
	hTPCTOFpid[ipart][charge]->Fill(param);

      } /* end of loop over particle IDs */      
    } /* end of loop over tracks */
  } /* end of loop over events */
  
  /* start stopwatch */
  timer.Stop();
  timer.Print();
  
  /* output */
  TFile *fileout = TFile::Open(Form("TPCTOFpid.%s", filename), "RECREATE");
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hTPCTOFpid[ipart][icharge]->Write();
    }
  }
  fileout->Close();
  
}

