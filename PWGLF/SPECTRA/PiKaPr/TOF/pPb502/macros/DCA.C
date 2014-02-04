#include "CommonDefs.C"

enum EMCHisto_t {
  kPrimary,
  kWeakDecay,
  kMaterial,
  kNMCHistos
};
const Char_t *mchistoName[kNMCHistos] = {
  "primary",
  "weakdecay",
  "material"
};

DCAdata(const Char_t *filename, Int_t evMax = kMaxInt, Int_t startEv = 0)
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
  for (Int_t ibin = 0; ibin < NdcaBins + 1; ibin++)
    dcaBin[ibin] = dcaMin + ibin * dcaStep;
  
  /* histos */
  TH3I *hDCAopen[AliPID::kSPECIES][kNCharges];
  TH3I *hDCAcut[AliPID::kSPECIES][kNCharges];
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hDCAopen[ipart][icharge] = new TH3I(Form("hDCAopen_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin, NptBins, ptBin, NdcaBins, dcaBin);
      hDCAcut[ipart][icharge] = new TH3I(Form("hDCAcut_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin, NptBins, ptBin, NdcaBins, dcaBin);
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
  Int_t charge, index;
  UShort_t dedxN;
  Double_t cent, p, pt, mt, dca, tofsignal, tpcsignal, tpctofsignal;
  Double_t dedx, bethe, deltadedx, dedx_sigma, ptpc;
  Double_t time, time_sigma, timezero, timezero_sigma, tof, tof_sigma, texp, texp_sigma, deltat, deltat_sigma;

  /* open TZERO calibfile */
  TH1 *hCentrality_TZEROA_mean = NULL;
  TH1 *hCentrality_TZEROA_sigma = NULL;
  TH1 *hCentrality_TZEROC_mean = NULL;
  TH1 *hCentrality_TZEROC_sigma = NULL;
  TH1 *hCentrality_TOF_mean = NULL;
  TH1 *hCentrality_TOF_TZEROA_mean = NULL;
  TH1 *hCentrality_TOF_TZEROC_mean = NULL;
  TH1 *hCentrality_TOF_TZEROTOF_mean = NULL;
  TH1 *hCentrality_TZEROA_reso = NULL;
  TH1 *hCentrality_TZEROC_reso = NULL;
  TH1 *hCentrality_TZEROAND_reso = NULL;
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


  for (Int_t iev = startEv; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 100000 == 0) printf("iev = %d\n", iev);
    /* check event */
    if (!analysisEvent->AcceptEvent(acceptEventType)) continue;

    /*** ACCEPTED EVENT ***/

    /* apply time-zero TOF correction */
    analysisEvent->ApplyTimeZeroTOFCorrection();

    /* get centrality */
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
      /* check accepted track (no DCA cut) */
      if (!analysisTrack->AcceptTrack(kFALSE)) continue;
      /* get charge */
      charge = analysisTrack->GetSign() > 0. ? kPositive : kNegative;
      /* check TPC pid */
      if (!analysisTrack->HasTPCPID()) continue;
      /* check TOF pid */
      if (!analysisTrack->HasTOFPID(hEnabledFlag)) continue;

      /*** ACCEPTED TRACK WITH TPC+TOF PID ***/

      /* get track info */
      p = analysisTrack->GetP();
      pt = analysisTrack->GetPt();
      dca = analysisTrack->GetImpactParameter(0);

      /* get TPC info */
      dedx = analysisTrack->GetTPCdEdx();
      dedxN = analysisTrack->GetTPCdEdxN();
      ptpc = analysisTrack->GetTPCmomentum();
      
      /* apply expected time correction */
      //      analysisTrack->ApplyTOFExpectedTimeCorrection();
      
      /* get TOF info */
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

      tof = time - timezero;
      tof_sigma = TMath::Sqrt(time_sigma * time_sigma + timezero_sigma * timezero_sigma);
      
      /* loop over particle IDs */
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
	
	/* check rapidity */
	if (analysisTrack->GetY(AliPID::ParticleMass(ipart)) - rapidityShift > rapidityMaxCut || 
	    analysisTrack->GetY(AliPID::ParticleMass(ipart)) - rapidityShift < rapidityMinCut) continue;

	/*** ACCEPTED TRACK WITHIN CORRECT RAPIDTY ***/
	
	/* TPC signal */
	bethe = tpcResponse->GetExpectedSignal(ptpc, ipart);
	deltadedx = dedx - bethe;
	dedx_sigma = tpcResponse->GetExpectedSigma(ptpc, dedxN, ipart);
	tpcsignal = deltadedx / dedx_sigma;
	
	/* TOF expected time */
	texp = analysisTrack->GetTOFExpTime(ipart);
	texp_sigma = analysisTrack->GetTOFExpTimeSigma(ipart);
	
	/* TOF signal */
	deltat = tof - texp;
	deltat_sigma = TMath::Sqrt(tof_sigma * tof_sigma + texp_sigma * texp_sigma);
	tofsignal = deltat / deltat_sigma;
	
	/* TPC+TOF signal */
	tpctofsignal = TMath::Sqrt(tpcsignal * tpcsignal + tofsignal * tofsignal);

	/* check PID cuts */
	if (tpctofsignal > 2.) continue;

	/* fill histo */
	hDCAopen[ipart][charge]->Fill(cent, pt, dca);

	/* check accepted track (with DCA cut) */
	if (!analysisTrack->AcceptTrack(kTRUE)) continue;

	/* fill histo */
	hDCAcut[ipart][charge]->Fill(cent, pt, dca);
	
      } /* end of loop over particle IDs */
    } /* end of loop over tracks */
  } /* end of loop over events */
  
  /* stop stopwatch */
  timer.Stop();
  timer.Print();
  
  /* output */
  TFile *fileout = TFile::Open(Form("DCAdata.%s", filename), "RECREATE");
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hDCAopen[ipart][icharge]->Write();
      hDCAcut[ipart][icharge]->Write();
    }
  }

  fileout->Close();
  
}

/**************************************************************/

DCAmc(const Char_t *filename, Int_t evMax = kMaxInt, Int_t startEv = 0)
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
  for (Int_t ibin = 0; ibin < NdcaBins + 1; ibin++)
    dcaBin[ibin] = dcaMin + ibin * dcaStep;
  
  /* histos */
  TH3I *hDCAopen[kNMCHistos][AliPID::kSPECIES][kNCharges];
  TH3I *hDCAcut[kNMCHistos][AliPID::kSPECIES][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNMCHistos; ihisto++ ) {
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hDCAopen[ihisto][ipart][icharge] = new TH3I(Form("hDCAopen_%s_%s_%s", mchistoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin, NptBins, ptBin, NdcaBins, dcaBin);
	hDCAcut[ihisto][ipart][icharge] = new TH3I(Form("hDCAcut_%s_%s_%s", mchistoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin, NptBins, ptBin, NdcaBins, dcaBin);
      }
    }
  }
  
  /**************************************************************/
  /**************************************************************/
  /**************************************************************/

  /* start stopwatch */
  TStopwatch timer;
  timer.Start();

  /* loop over events */
  Int_t part, charge, type;
  Double_t cent, pt, mt, dca;

  for (Int_t iev = startEv; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 100000 == 0) printf("iev = %d\n", iev);
    /* check event */
    if (!analysisEvent->AcceptEvent(acceptEventType)) continue;

    /*** ACCEPTED EVENT ***/

    /* get centrality */
    cent = analysisEvent->GetCentralityPercentile(centralityEstimator);

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check defined PID */
      part = analysisTrack->GetMCPID();
      if (part < 0 || analysisTrack->GetSign() == 0.) continue;
      /* check accepted track (no DCA cut) */
      if (!analysisTrack->AcceptTrack(kFALSE)) continue;
      /* check rapidity */
      if (analysisTrack->GetY(AliPID::ParticleMass(part)) - rapidityShift > rapidityMaxCut || 
	  analysisTrack->GetY(AliPID::ParticleMass(part)) - rapidityShift < rapidityMinCut) continue;
      
      /*** ACCEPTED TRACK WITHIN CORRECT RAPIDTY ***/
      
      /* get track info */
      pt = analysisTrack->GetPt();
      dca = analysisTrack->GetImpactParameter(0);
      charge = analysisTrack->GetMCCharge() > 0. ? kPositive : kNegative;
      if (analysisTrack->IsMCPrimary())
	type = kPrimary;
      else if (analysisTrack->IsMCSecondaryWeakDecay())
	type = kWeakDecay;
      else if (analysisTrack->IsMCSecondaryMaterial())
	type = kMaterial;
      else
	continue;

      /* fill histo */
      hDCAopen[type][part][charge]->Fill(cent, pt, dca);
      
      /* check accepted track (with DCA cut) */
      if (!analysisTrack->AcceptTrack(kTRUE)) continue;
      
      /* fill histo */
      hDCAcut[type][part][charge]->Fill(cent, pt, dca);
      
    } /* end of loop over tracks */
  } /* end of loop over events */
  
  /* stop stopwatch */
  timer.Stop();
  timer.Print();
  
  /* output */
  TFile *fileout = TFile::Open(Form("DCAmc.%s", filename), "RECREATE");
  for (Int_t ihisto = 0; ihisto < kNMCHistos; ihisto++) {
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hDCAopen[ihisto][ipart][icharge]->Write();
	hDCAcut[ihisto][ipart][icharge]->Write();
      }
    }
  }

  fileout->Close();
  
}

/**************************************************************/

enum EFitParams_t {
  kIntegralCounts,
  kPrimaryCounts,
  kWeakDecayCounts,
  kMaterialCounts,
  kPrimaryIntegral,
  kWeakDecayIntegral,
  kMaterialIntegral,
  kDataCutIntegral,
  kPrimaryCutIntegral,
  kWeakDecayCutIntegral,
  kMaterialCutIntegral,
  kNFitParams
};

/* fit output params name */
const Char_t *fitParamName[kNFitParams] = {
  "IntegralCounts",
  "PrimaryCounts",
  "WeakDecayCounts",
  "MaterialCounts",
  "PrimaryIntegral",
  "WeakDecayIntegral",
  "MaterialIntegral",
  "DataCutIntegral",
  "PrimaryCutIntegral",
  "WeakDecayCutIntegral",
  "MaterialCutIntegral"
};

/* fit output params title */
const Char_t *fitParamTitle[kNFitParams] = {
  "Integral counts;p_{T} (GeV/c);",
  "Primary counts;p_{T} (GeV/c);",
  "Weak-decay counts;p_{T} (GeV/c);",
  "Material counts;p_{T} (GeV/c);",
  "Primary integral;p_{T} (GeV/c);",
  "Weak-decay integral;p_{T} (GeV/c);",
  "Material integral;p_{T} (GeV/c);",
  "Data cut integral;p_{T} (GeV/c);",
  "Primary cut integral;p_{T} (GeV/c);",
  "Weak-decay cut integral;p_{T} (GeV/c);",
  "Material cut integral;p_{T} (GeV/c);"
};

/* fit ranges */
Double_t fitPtMin[AliPID::kSPECIES] = {0.5, 0.5, 0.3, 0.4, 0.5};
Double_t fitPtMax[AliPID::kSPECIES] = {2.0, 2.0, 2.0, 2.0, 3.0};

/* rebin DCA */
Int_t rebindca = 10;

DCAdeltafeed(const Char_t *datafilename, const Char_t *mcfilename)
{
  for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
    DCAdeltafeed(datafilename, mcfilename, 2, icharge, -1);
    DCAdeltafeed(datafilename, mcfilename, 4, icharge, -1);
    for (Int_t icent = 0; icent < NcentralityBins; icent++) {
      DCAdeltafeed(datafilename, mcfilename, 2, icharge, icent);
      DCAdeltafeed(datafilename, mcfilename, 4, icharge, icent);
    }
  }
}

DCAdeltafeed(const Char_t *datafilename, const Char_t *mcfilename, Int_t ipart, Int_t icharge, Int_t icent, Float_t ptMin = -1., Float_t ptMax = -1., Bool_t checkHistoFlag = kFALSE)
{

  printf("****************************************\n");
  printf("RUNNING DCA FIT:\n");
  printf("RAPIDITY-CUT:   %s\n", AliPID::ParticleName(ipart));
  printf("CHARGE:         %s\n", chargeName[icharge]);
  printf("PARTICLE:       %s\n", AliPID::ParticleName(ipart));
  printf("CENTRALITY BIN: %d\n", icent);
  printf("****************************************\n");

  /* open data */
  TFile *datafilein = TFile::Open(datafilename);
  TFile *mcfilein = TFile::Open(mcfilename);

  /* get histos */
  TH3I *hDCAopen = (TH3I *)datafilein->Get(Form("hDCAopen_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAopen_primary = (TH3I *)mcfilein->Get(Form("hDCAopen_primary_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAopen_weakdecay = (TH3I *)mcfilein->Get(Form("hDCAopen_weakdecay_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAopen_material = (TH3I *)mcfilein->Get(Form("hDCAopen_material_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAcut = (TH3I *)datafilein->Get(Form("hDCAcut_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAcut_primary = (TH3I *)mcfilein->Get(Form("hDCAcut_primary_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAcut_weakdecay = (TH3I *)mcfilein->Get(Form("hDCAcut_weakdecay_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
  TH3I *hDCAcut_material = (TH3I *)mcfilein->Get(Form("hDCAcut_material_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));

  /* setup centrality range */
  if (icent < 0 || icent >= NcentralityBins) {
    printf("WARNING: undefined centrality -> using 00-90\% range\n");
    hDCAopen->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAopen_primary->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAopen_weakdecay->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAopen_material->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut_primary->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut_weakdecay->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut_material->GetXaxis()->SetRange(1, NcentralityBins);
  }
  else {
    printf("***** FITTING CENTRALITY-BIN [%02d, %02d] %% *****\n", centralityBin[icent], centralityBin[icent + 1]);
    hDCAopen->GetXaxis()->SetRange(icent + 1, icent + 1);
#if 0
    hDCAopen_primary->GetXaxis()->SetRange(icent + 1, icent + 1);
    hDCAopen_weakdecay->GetXaxis()->SetRange(icent + 1, icent + 1);
    hDCAopen_material->GetXaxis()->SetRange(icent + 1, icent + 1);
#else
    hDCAopen_primary->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAopen_weakdecay->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAopen_material->GetXaxis()->SetRange(1, NcentralityBins);
#endif
    hDCAcut->GetXaxis()->SetRange(icent + 1, icent + 1);
#if 0
    hDCAcut_primary->GetXaxis()->SetRange(icent + 1, icent + 1);
    hDCAcut_weakdecay->GetXaxis()->SetRange(icent + 1, icent + 1);
    hDCAcut_material->GetXaxis()->SetRange(icent + 1, icent + 1);
#else
    hDCAcut_primary->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut_weakdecay->GetXaxis()->SetRange(1, NcentralityBins);
    hDCAcut_material->GetXaxis()->SetRange(1, NcentralityBins);
#endif
  }

  /* setup pt range */
  Bool_t requestedRange = kFALSE;
  if (ptMin > -0.001 && ptMax > -0.001 && ptMax > ptMin) {
    printf("***** FITTING PT-BIN [%f, %f] GeV/c *****\n", ptMin, ptMax);
    requestedRange = kTRUE;
    hDCAopen->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAopen_primary->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAopen_weakdecay->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAopen_material->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAcut->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAcut_primary->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAcut_weakdecay->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
    hDCAcut_material->GetYaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
  }

  /* output */
  Char_t outfilename[1024];
  if (icent < 0 || icent >= NcentralityBins)
    sprintf(outfilename, "DCAdeltafeed_cent0090_%s_%s.root", AliPID::ParticleName(ipart), chargeName[icharge]);
  else {
    sprintf(outfilename, "DCAdeltafeed_cent%02d%02d_%s_%s.root", centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge]);
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
      hDCAopen->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAopen_primary->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAopen_weakdecay->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAopen_material->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAcut->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAcut_primary->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAcut_weakdecay->GetYaxis()->SetRange(ipt + 1, ipt + 1);
      hDCAcut_material->GetYaxis()->SetRange(ipt + 1, ipt + 1);
    }
    
    /* dca projections */
    TH1 *hDCAopen_py = hDCAopen->Project3D("z");
    TH1 *hDCAopen_primary_py = hDCAopen_primary->Project3D("z");
    TH1 *hDCAopen_weakdecay_py = hDCAopen_weakdecay->Project3D("z");
    TH1 *hDCAopen_material_py = hDCAopen_material->Project3D("z");
    TH1 *hDCAcut_py = hDCAcut->Project3D("z");
    TH1 *hDCAcut_primary_py = hDCAcut_primary->Project3D("z");
    TH1 *hDCAcut_weakdecay_py = hDCAcut_weakdecay->Project3D("z");
    TH1 *hDCAcut_material_py = hDCAcut_material->Project3D("z");

    /* rebin */
    hDCAopen_py->Rebin(rebindca);
    hDCAopen_primary_py->Rebin(rebindca);
    hDCAopen_weakdecay_py->Rebin(rebindca);
    hDCAopen_material_py->Rebin(rebindca);
    hDCAcut_py->Rebin(rebindca);
    hDCAcut_primary_py->Rebin(rebindca);
    hDCAcut_weakdecay_py->Rebin(rebindca);
    hDCAcut_material_py->Rebin(rebindca);

    /* check histos if requested */
    if (checkHistoFlag) {
      TCanvas *cCheckHisto = new TCanvas("cCheckHisto");
      cCheckHisto->Divide(2, 4);
      cCheckHisto->cd(1);
      hDCAopen_py->Draw();
      cCheckHisto->cd(2);
      hDCAcut_py->Draw();
      cCheckHisto->cd(3);
      hDCAopen_primary_py->Draw();
      cCheckHisto->cd(4);
      hDCAcut_primary_py->Draw();
      cCheckHisto->cd(5);
      hDCAopen_weakdecay_py->Draw();
      cCheckHisto->cd(6);
      hDCAcut_weakdecay_py->Draw();
      cCheckHisto->cd(7);
      hDCAopen_material_py->Draw();
      cCheckHisto->cd(8);
      hDCAcut_material_py->Draw();
      return;
    }

    Double_t param[kNFitParams];
    Double_t param_err[kNFitParams];

    param[kPrimaryIntegral] = hDCAopen_primary_py->Integral();
    param_err[kPrimaryIntegral] = TMath::Sqrt(hDCAopen_primary_py->Integral());
    param[kWeakDecayIntegral] = hDCAopen_weakdecay_py->Integral();
    param_err[kWeakDecayIntegral] = TMath::Sqrt(hDCAopen_weakdecay_py->Integral());
    param[kMaterialIntegral] = hDCAopen_material_py->Integral();
    param_err[kMaterialIntegral] = TMath::Sqrt(hDCAopen_material_py->Integral());
    
    param[kDataCutIntegral] = hDCAcut_py->Integral();
    param_err[kDataCutIntegral] = TMath::Sqrt(hDCAcut_py->Integral());
    param[kPrimaryCutIntegral] = hDCAcut_primary_py->Integral();
    param_err[kPrimaryCutIntegral] = TMath::Sqrt(hDCAcut_primary_py->Integral());
    param[kWeakDecayCutIntegral] = hDCAcut_weakdecay_py->Integral();
    param_err[kWeakDecayCutIntegral] = TMath::Sqrt(hDCAcut_weakdecay_py->Integral());
    param[kMaterialCutIntegral] = hDCAcut_material_py->Integral();
    param_err[kMaterialCutIntegral] = TMath::Sqrt(hDCAcut_material_py->Integral());

    /* fit */
    DCAfit(hDCAopen_py, hDCAopen_primary_py, hDCAopen_weakdecay_py, hDCAopen_material_py, param, param_err, canvas);

    
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
    delete hDCAopen_py;
    delete hDCAopen_primary_py;
    delete hDCAopen_weakdecay_py;
    delete hDCAopen_material_py;
    delete hDCAcut_py;
    delete hDCAcut_primary_py;
    delete hDCAcut_weakdecay_py;
    delete hDCAcut_material_py;

  }

  /* check requested pt-range */
  if (requestedRange)
    return;

  /*** POST-ANALYSIS ***/

  TDirectory *postDir = fileout->mkdir("PostAnalysis");
 
  /* compute fractions */
  TH1D *hPrimaryFraction = new TH1D(*hFitParamHisto[kPrimaryCounts]);
  hPrimaryFraction->Divide(hFitParamHisto[kPrimaryCounts], hFitParamHisto[kIntegralCounts], 1., 1., "B");
  TH1D *hWeakDecayFraction = new TH1D(*hFitParamHisto[kWeakDecayCounts]);
  hWeakDecayFraction->Divide(hFitParamHisto[kWeakDecayCounts], hFitParamHisto[kIntegralCounts], 1., 1., "B");
  TH1D *hMaterialFraction = new TH1D(*hFitParamHisto[kMaterialCounts]);
  hMaterialFraction->Divide(hFitParamHisto[kMaterialCounts], hFitParamHisto[kIntegralCounts], 1., 1., "B");

  /* compute scale factors */
  TH1D *hPrimaryScale = new TH1D(*hFitParamHisto[kPrimaryCounts]);
  hPrimaryScale->Divide(hFitParamHisto[kPrimaryCounts], hFitParamHisto[kPrimaryIntegral], 1., 1., "B");
  TH1D *hWeakDecayScale = new TH1D(*hFitParamHisto[kWeakDecayCounts]);
  hWeakDecayScale->Divide(hFitParamHisto[kWeakDecayCounts], hFitParamHisto[kWeakDecayIntegral], 1., 1., "B");
  TH1D *hMaterialScale = new TH1D(*hFitParamHisto[kMaterialCounts]);
  hMaterialScale->Divide(hFitParamHisto[kMaterialCounts], hFitParamHisto[kMaterialIntegral], 1., 1., "B");

  /* compute cut integrals */
  TH1D *hPrimaryCutIntegral = new TH1D(*hFitParamHisto[kPrimaryCutIntegral]);
  hPrimaryCutIntegral->Multiply(hPrimaryScale);
  TH1D *hWeakDecayCutIntegral = new TH1D(*hFitParamHisto[kWeakDecayCutIntegral]);
  hWeakDecayCutIntegral->Multiply(hWeakDecayScale);
  TH1D *hMaterialCutIntegral = new TH1D(*hFitParamHisto[kMaterialCutIntegral]);
  hMaterialCutIntegral->Multiply(hMaterialScale);

  /* compute cut fractions */
  TH1D *hPrimaryCutFraction = new TH1D(*hPrimaryCutIntegral);
  hPrimaryCutFraction->Divide(hPrimaryCutIntegral, hFitParamHisto[kDataCutIntegral], 1., 1., "B");
  TH1D *hWeakDecayCutFraction = new TH1D(*hWeakDecayCutIntegral);
  hWeakDecayCutFraction->Divide(hWeakDecayCutIntegral, hFitParamHisto[kDataCutIntegral], 1., 1., "B");
  TH1D *hMaterialCutFraction = new TH1D(*hMaterialCutIntegral);
  hMaterialCutFraction->Divide(hMaterialCutIntegral, hFitParamHisto[kDataCutIntegral], 1., 1., "B");


  /*** OUTPUT ***/

  /* write fir params histos */
  fitDir->cd();
  for (Int_t iparam = 0; iparam < kNFitParams; iparam++)
    hFitParamHisto[iparam]->Write();
  /* write post-analysis histos */
  postDir->cd();
  hPrimaryFraction->Write("hPrimaryFraction");
  hWeakDecayFraction->Write("hWeakDecayFraction");
  hMaterialFraction->Write("hMaterialFraction");
  hPrimaryCutFraction->Write("hPrimaryCutFraction");
  hWeakDecayCutFraction->Write("hWeakDecayCutFraction");
  hMaterialCutFraction->Write("hMaterialCutFraction");

  /* cleanup */
  delete hPrimaryFraction;
  delete hWeakDecayFraction;
  delete hMaterialFraction;
  delete hPrimaryScale;
  delete hWeakDecayScale;
  delete hMaterialScale;
  delete hPrimaryCutIntegral;
  delete hWeakDecayCutIntegral;
  delete hMaterialCutIntegral;
  delete hPrimaryCutFraction;
  delete hWeakDecayCutFraction;
  delete hMaterialCutFraction;
  delete canvas;
  for (Int_t iparam = 0; iparam < kNFitParams; iparam++)
    delete hFitParamHisto[iparam];

  /* close file */
  fileout->Close();


}

/**************************************************************/

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


Double_t
DCAfit(TH1 *hData, TH1 *hPrimary, TH1 *hWeakDecay, TH1 *hMaterial, Double_t *param = NULL, Double_t *param_err = NULL, TCanvas *canvas = NULL)
{

  /** ROOFIT ***/
  gSystem->Load("libRooFit");
  using namespace RooFit;

  /*** DEFINE FIT RANGE ***/

  printf("***** FIT RANGE DEFINITION *****\n");

  /* check material histogram to define min/max fit range */
  Double_t rangeMin = TMath::Max(dcaMin, TOFpid_histomin(hMaterial));
  Double_t rangeMax = TMath::Min(dcaMax, TOFpid_histomax(hMaterial));
  /* fix zeroes */
  TOFpid_checkneg(hMaterial);

  /* define range */
  RooRealVar x("x", "DCA_{xy}", 0., -5., 5., "");
  printf("FIT RANGE DEFINED: %f -> %f\n", dcaMin, dcaMax);
  printf("********************************\n");

  /*** DEFINE HISTOGRAM DATA ***/
  
  /* define data to fit and background from input histogram */
  RooDataHist hdata("hdata", "hdata", x, hData);
  RooDataHist hprimary("hprimary", "hprimary", x, hPrimary);
  RooDataHist hweakdecay("hweakdecay", "hweakdecay", x, hWeakDecay);
  RooDataHist hmaterial("hmaterial", "hmaterial", x, hMaterial);

  /*** DEFINE SHAPES ***/

  RooHistPdf primary("primary", "primary", x, hprimary);
  RooHistPdf weakdecay("weakdecay", "weakdecay", x, hweakdecay);
  RooHistPdf material("material", "material", x, hmaterial);

  /*** DEFINE MODEL ***/

  Double_t integral = hdata.sumEntries();
  RooRealVar nprimary("nprimary", "nprimary", 0.80 * integral, 0., integral);
  RooRealVar nweakdecay("nweakdecay", "nweakdecay", 0.1 * integral, 0., integral);
  RooRealVar nmaterial("nmaterial", "nmaterial", 0.1 * integral, 0., integral);

  RooAddPdf model("model", "model p.d.f.", RooArgList(primary, weakdecay, material), RooArgList(nprimary, nweakdecay, nmaterial));

  /*** FIT ***/

  if (canvas) canvas->cd();
  model.fitTo(hdata, Extended(kTRUE), SumW2Error(kFALSE));

  /*** DRAW ***/

  RooPlot *xframe = x.frame();
  hdata.plotOn(xframe, XErrorSize(0), DrawOption("PZ"));
  model.plotOn(xframe, LineWidth(2), Precision(1.e-4));
  model.plotOn(xframe, Components(primary), LineWidth(2), LineColor(kRed), Precision(1.e-4));
  model.plotOn(xframe, Components(weakdecay), LineWidth(2), LineStyle(kDashed), LineColor(kCyan+1));
  model.plotOn(xframe, Components(material), LineWidth(2), LineStyle(kDashed), LineColor(kGreen+1));
  xframe->SetMinimum(0.1);
  xframe->Draw();
  if (canvas) canvas->Update();

  printf("*****************************\n");
  printf("***** FRACTIONS *****\n");
  printf("primary     = %f +- %f\n", nprimary.getVal(), nprimary.getError());
  printf("weak-decay  = %f +- %f\n", nweakdecay.getVal(), nweakdecay.getError());
  printf("material    = %f +- %f\n", nmaterial.getVal(), nmaterial.getError());
  printf("******************\n");

  /*** OUTPUT FIT PARAMS ***/
  
  param[kIntegralCounts] = integral;
  param_err[kIntegralCounts] = sqrt(integral);
  param[kPrimaryCounts] = nprimary.getVal();
  param_err[kPrimaryCounts] = nprimary.getError();
  param[kWeakDecayCounts] = nweakdecay.getVal();
  param_err[kWeakDecayCounts] = nweakdecay.getError();
  param[kMaterialCounts] = nmaterial.getVal();
  param_err[kMaterialCounts] = nmaterial.getError();

  return;

}

enum EFitFunc_t {
  kExpo,
  kInverseExpo,
  kPowerLaw,
  kInversePowerLaw,
  kExpoPowerLaw
};

DCAdeltafeed_param()
{
#if 0
  DCAdeltafeed_param(2, 0, "WeakDecayCutFraction", kExpo);
  DCAdeltafeed_param(2, 1, "WeakDecayCutFraction", kExpo);
  DCAdeltafeed_param(4, 0, "WeakDecayCutFraction", kExpo);
  DCAdeltafeed_param(4, 1, "WeakDecayCutFraction", kExpo);
  DCAdeltafeed_param(2, 0, "MaterialCutFraction", kPowerLaw);
  DCAdeltafeed_param(2, 1, "MaterialCutFraction", kPowerLaw);
  DCAdeltafeed_param(4, 0, "MaterialCutFraction", kPowerLaw);
  DCAdeltafeed_param(4, 1, "MaterialCutFraction", kPowerLaw);
#endif
  DCAdeltafeed_param_bothcharges(2, 0, "PrimaryCutFraction", kInverseExpo);
  DCAdeltafeed_param_bothcharges(2, 1, "PrimaryCutFraction", kInverseExpo);
  DCAdeltafeed_param_bothcharges(4, 0, "PrimaryCutFraction", kInverseExpo);
  DCAdeltafeed_param_bothcharges(4, 1, "PrimaryCutFraction", kInverseExpo);
}

DCAdeltafeed_param_bothcharges(Int_t ipart, Int_t icharge, const Char_t *name, Int_t fitFunc)
{

  TVirtualFitter::SetMaxIterations(1000000);
  
  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");

  Char_t outfilename[1024];
  TH1D *hDeltaFeed_mb[kNCharges];
  TH1D *hDeltaFeed_cent[kNCharges][NcentralityBins];
  for (Int_t iicharge = 0; iicharge < kNCharges; iicharge++) {

    /* get minimum-bias results */
    sprintf(outfilename, "DCAdeltafeed_cent0090_%s_%s.root", AliPID::ParticleName(ipart), chargeName[iicharge]);
    TFile *filein = TFile::Open(outfilename);
    hDeltaFeed_mb[iicharge] = (TH1D *)filein->Get(Form("PostAnalysis/h%s", name));
    if (!hDeltaFeed_mb[iicharge]) {
      printf("cannot find PostAnalysis/h%s in %s\n", name, outfilename);
      return;
    }

    /* get centrality-binned results */
    for (Int_t icent = 0; icent < NcentralityBins; icent++) {
      sprintf(outfilename, "DCAdeltafeed_cent%02d%02d_%s_%s.root", centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[iicharge]);
      filein = TFile::Open(outfilename);
      if (!filein || !filein->IsOpen()) continue;
      hDeltaFeed_cent[iicharge][icent] = (TH1D *)filein->Get(Form("PostAnalysis/h%s", name));
      if (!hDeltaFeed_cent[iicharge][icent]) {
	printf("cannot find PostAnalysis/h%s in %s\n", name, outfilename);
	return;
      }
    }
  }

  /* define delta feed-down function */
  switch (fitFunc) {
  case kExpo:
    printf("expo function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expo", "[0] * ([1] + [2] * TMath::Exp([3] * x))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 0.01);
    fDeltaFeed->SetParLimits(1, 0., 0.1);
    fDeltaFeed->SetParameter(2, 0.1);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kInverseExpo:
    printf("inverse-expo function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expo", "[0] * ([1] + [2] * TMath::Exp([3] * x))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 1.);
    fDeltaFeed->SetParLimits(1, 0.9, 1.0);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParLimits(2, -1., 0.);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kPowerLaw:
    printf("powerlaw function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_powerlaw", "[0] * ([1] + [2] * TMath::Power(x, [3]))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 0.01);
    fDeltaFeed->SetParLimits(1, 0., 0.1);
    fDeltaFeed->SetParameter(2, 0.1);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kInversePowerLaw:
    printf("inverse-powerlaw function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_powerlaw", "[0] * ([1] + [2] * TMath::Power(x, [3]) + [4] * x)", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 1.);
    fDeltaFeed->SetParLimits(1, 0.9, 1.);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParameter(3, -1.);
    fDeltaFeed->SetParameter(4, 0.);
    fDeltaFeed->SetParLimits(4, 0., 1.)
    break;
  case kExpoPowerLaw:
    printf("expo+powerlaw function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expopowerlaw", "[0] * ([1] + [2] * TMath::Exp([3] * x) + [4] * TMath::Exp([5] * x))", 0.5, 5.);

    fDeltaFeed->SetParameter(1, 0.95);
    fDeltaFeed->SetParLimits(1, 0.9, 1.);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParLimits(2, -1000., 0.);
    fDeltaFeed->SetParameter(3, -1.);
    fDeltaFeed->SetParameter(4, -100.);
    fDeltaFeed->SetParLimits(4, -100., 0.);
    fDeltaFeed->SetParameter(5, -20.);
    fDeltaFeed->SetParLimits(5, -100., 10.);
    break;
  default:
    printf("fit function not defined\n");
    return;
    break;
  }
  fDeltaFeed->FixParameter(0, 1.);

  /* output */
  TFile *fileout = TFile::Open(Form("DCAdeltafeed_param_%s_%s_%s.root", name, AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");

  /* extra-material function */
  TF1 *fExtraMaterial = new TF1("fExtraMaterial", "1. + [0] * TMath::Exp([1] * x)", 0., 5.);
  fExtraMaterial->SetParameter(0, 0.1);
  fExtraMaterial->SetParameter(1, -0.5);

  /* loop over centrality bins */
  for (Int_t icent = -1; icent < NcentralityBins; icent++) {

    printf(" icent = %d\n", icent);
    
    /* fit material from proton-antiproton ratio and remove it */
    if (ipart == AliPID::kProton) {
      printf("fitting extra-material for positive protons\n");
      if (icent == -1) {
	TH1D *hr = new TH1D(*hDeltaFeed_mb[kPositive]);
	hr->Divide(hDeltaFeed_mb[kNegative]);
	hr->Fit(fExtraMaterial, "R", "", fitPtMin[ipart], fitPtMax[ipart]);
	hDeltaFeed_mb[kPositive]->Divide(fExtraMaterial);
      }
      else {
	TH1D *hr = new TH1D(*hDeltaFeed_cent[kPositive][icent]);
	hr->Divide(hDeltaFeed_cent[kNegative][icent]);
	hr->Fit(fExtraMaterial, "R", "", fitPtMin[ipart], fitPtMax[ipart]);
	hDeltaFeed_cent[kPositive][icent]->Divide(fExtraMaterial);
      }
    }
    
    /* build graph with both charges */
    TGraphErrors *gDeltaFeed = new TGraphErrors();
    Int_t npoints = 0;
    Double_t pt, pte, val, vale;
    for (Int_t iicharge = 0; iicharge < kNCharges; iicharge++) {
      for (Int_t ipt = 0; ipt < NptBins; ipt++) {
	if (icent == -1) {
	  if (hDeltaFeed_mb[iicharge]->GetBinError(ipt + 1) == 0.)
	    continue;
	  pt = hDeltaFeed_mb[iicharge]->GetBinCenter(ipt + 1);
	  pte = 0.5 * hDeltaFeed_mb[iicharge]->GetBinWidth(ipt + 1);
	  val = hDeltaFeed_mb[iicharge]->GetBinContent(ipt + 1);
	  vale = hDeltaFeed_mb[iicharge]->GetBinError(ipt + 1);
	  gDeltaFeed->SetPoint(npoints, pt, val);
	  gDeltaFeed->SetPointError(npoints, pte, vale);
	  npoints++;
	}
	else {
	  if (hDeltaFeed_cent[iicharge][icent]->GetBinError(ipt + 1) == 0.)
	    continue;
	  pt = hDeltaFeed_cent[iicharge][icent]->GetBinCenter(ipt + 1);
	  pte = 0.5 * hDeltaFeed_cent[iicharge][icent]->GetBinWidth(ipt + 1);
	  val = hDeltaFeed_cent[iicharge][icent]->GetBinContent(ipt + 1);
	  vale = hDeltaFeed_cent[iicharge][icent]->GetBinError(ipt + 1);
	  gDeltaFeed->SetPoint(npoints, pt, val);
	  gDeltaFeed->SetPointError(npoints, pte, vale);
	  npoints++;
	}
      }
    }

    /* fit graph */
    Int_t fitres = -1;
    fitres = gDeltaFeed->Fit(fDeltaFeed, "R", "", 0.5, 2.0);
    if (fitres != 0) {
      fDeltaFeed->SetParLimits(1, 0.9, 1.0);
      fDeltaFeed->SetParameter(2, -0.1);
      fDeltaFeed->SetParLimits(2, -1., 0.);
      fDeltaFeed->SetParameter(3, -1.);
    }
    fitres = -1;
    while(fitres != 0)
      fitres = gDeltaFeed->Fit(fDeltaFeed, "R", "", fitPtMin[ipart], fitPtMax[ipart]);
    
    gDeltaFeed->Draw("ap*");
    fDeltaFeed->Draw("same");
    gPad->Update();
    
    printf("done part = %d, charge = %d, cent = %d\n", ipart, icharge, icent);

    /* build fit profile */
    TProfile *pDeltaFeed = new TProfile("pDeltaFeed", "", NptBins, ptBin);
    HistoUtils_Function2Profile(fDeltaFeed, pDeltaFeed);

    /* put material back */
    if (ipart == AliPID::kProton && icharge == kPositive) {
      if (icent == -1) {
	hDeltaFeed_mb[kPositive]->Multiply(fExtraMaterial);
	pDeltaFeed->Multiply(fExtraMaterial);
      }
      else {
	hDeltaFeed_cent[kPositive][icent]->Multiply(fExtraMaterial);
	pDeltaFeed->Multiply(fExtraMaterial);
      }
    }
    
    /* write */
    fileout->cd();
    if (icent == -1) {
      hDeltaFeed_mb[icharge]->Write(Form("hDeltaFeed_MB"));
      fDeltaFeed->Write(Form("fDeltaFeed_MB"));
      pDeltaFeed->Write(Form("pDeltaFeed_MB"));
      fExtraMaterial->Write(Form("fExtraMaterial_MB"));
    }
    else {
      hDeltaFeed_cent[icharge][icent]->Write(Form("hDeltaFeed_cent%d", icent));
      fDeltaFeed->Write(Form("fDeltaFeed_cent%d", icent));
      pDeltaFeed->Write(Form("pDeltaFeed_cent%d", icent));
      fExtraMaterial->Write(Form("fExtraMaterial_cent%d", icent));
    }
  }
  
  fileout->Close();
}

DCAdeltafeed_param(Int_t ipart, Int_t icharge, const Char_t *name, Int_t fitFunc)
{

  TVirtualFitter::SetMaxIterations(1000000);

  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");
  
  /* get minimum-bias results */
  Char_t outfilename[1024];
  sprintf(outfilename, "DCAdeltafeed_cent0090_%s_%s.root", AliPID::ParticleName(ipart), chargeName[icharge]);
  TFile *filein = TFile::Open(outfilename);
  TH1D *hDeltaFeed_mb = (TH1D *)filein->Get(Form("PostAnalysis/h%s", name));
  if (!hDeltaFeed_mb) {
    printf("cannot find PostAnalysis/h%s in %s\n", name, outfilename);
    return;
  }

  /* get centrality-binned results */
  TH1D *hDeltaFeed_cent[NcentralityBins];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    sprintf(outfilename, "DCAdeltafeed_cent%02d%02d_%s_%s.root", centralityBin[icent], centralityBin[icent + 1], AliPID::ParticleName(ipart), chargeName[icharge]);
    filein = TFile::Open(outfilename);
    if (!filein || !filein->IsOpen()) continue;
    hDeltaFeed_cent[icent] = (TH1D *)filein->Get(Form("PostAnalysis/h%s", name));
    if (!hDeltaFeed_cent[icent]) {
      printf("cannot find PostAnalysis/h%s in %s\n", name, outfilename);
      return;
    }
  }
 
  /* define delta feed-down function */
  switch (fitFunc) {
  case kExpo:
    printf("expo function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expo", "[0] * ([1] + [2] * TMath::Exp([3] * x))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 0.01);
    fDeltaFeed->SetParLimits(1, 0., 0.1);
    fDeltaFeed->SetParameter(2, 0.1);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kInverseExpo:
    printf("inverse-expo function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expo", "[0] * ([1] + [2] * TMath::Exp([3] * x))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 1.);
    fDeltaFeed->SetParLimits(1, 0.9, 1.0);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kPowerLaw:
    printf("powerlaw function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_powerlaw", "[0] * ([1] + [2] * TMath::Power(x, [3]))", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 0.01);
    fDeltaFeed->SetParLimits(1, 0., 0.1);
    fDeltaFeed->SetParameter(2, 0.1);
    fDeltaFeed->SetParameter(3, -1.);
    break;
  case kInversePowerLaw:
    printf("inverse-powerlaw function selected\n");
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_powerlaw", "[0] * ([1] + [2] * TMath::Power(x, [3]) + [4] * x)", 0.5, 5.);
    fDeltaFeed->SetParameter(1, 1.);
    fDeltaFeed->SetParLimits(1, 0.9, 1.);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParameter(3, -1.);
    fDeltaFeed->SetParameter(4, 0.);
    fDeltaFeed->SetParLimits(4, 0., 1.)
    break;
  case kExpoPowerLaw:
    printf("expo+powerlaw function selected\n");
    //    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expopowerlaw", "[0] * ([1] + [2] * TMath::Exp([3] * x) + [4] * TMath::Exp([5] * x))", 0.5, 5.);
    TF1 *fDeltaFeed = new TF1("fDeltaFeed_expopowerlaw", "[0] * ([1] + [2] * TMath::Exp([3] * x) * TMath::Exp([5] / x))", 0.5, 5.);

    fDeltaFeed->SetParameter(1, 0.95);
    fDeltaFeed->SetParLimits(1, 0.9, 1.);
    fDeltaFeed->SetParameter(2, -0.1);
    fDeltaFeed->SetParLimits(2, -1000., 0.);
    fDeltaFeed->SetParameter(3, -1.);
    fDeltaFeed->SetParameter(4, -100.);
    fDeltaFeed->SetParLimits(4, -100., 0.);
    fDeltaFeed->SetParameter(5, -20.);
    fDeltaFeed->SetParLimits(5, -100., 10.);
    break;
  default:
    printf("fit function not defined\n");
    return;
    break;
  }
  fDeltaFeed->FixParameter(0, 1.);

  /* output */
  TFile *fileout = TFile::Open(Form("DCAdeltafeed_param_%s_%s_%s.root", name, AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");

  /* fit minimum-bias */
  TCanvas *cMinimumBias = new TCanvas("cMinimumBias");
  if (fitFunc == kExpoPowerLaw) {
    fDeltaFeed->FixParameter(4, 0.);
    hDeltaFeed_mb->Fit(fDeltaFeed, "IMRE", "", 1.0, 3.0);
    fDeltaFeed->FixParameter(1, fDeltaFeed->GetParameter(1));
    fDeltaFeed->FixParameter(2, fDeltaFeed->GetParameter(2));
    fDeltaFeed->FixParameter(3, fDeltaFeed->GetParameter(3));
    fDeltaFeed->ReleaseParameter(4);
    hDeltaFeed_mb->Fit(fDeltaFeed, "IMRE", "", 0.5, 3.0);
  }
  else {
    hDeltaFeed_mb->Fit(fDeltaFeed, "IMRE", "", 0.5, 3.0);
  }
  cMinimumBias->SetLogy();

  /* build fit profile */
  TProfile *pDeltaFeed = new TProfile("pDeltaFeed", "", NptBins, ptBin);
  HistoUtils_Function2Profile(fDeltaFeed, pDeltaFeed);

  /* save MB params */
  Double_t parMB[100];
  for (Int_t ipar = 0; ipar < fDeltaFeed->GetNpar(); ipar++)
    parMB[ipar] = fDeltaFeed->GetParameter(ipar);
  
  /* write */
  fileout->cd();
  hDeltaFeed_mb->Write("hDeltaFeed_mb");
  fDeltaFeed->Write("fDeltaFeed_mb");
  pDeltaFeed->Write("pDeltaFeed_mb");

  /* fix pt-depencence and release scale factor
     to fit centrality bins */
  TCanvas *cCentralityDependence = new TCanvas("cCentralityDependence");
  TH1D *hCentDep = new TH1D("hCentDep", "", NcentralityBins, centralityBin);
  //fDeltaFeed->ReleaseParameter(0);
  //  for (Int_t ipar = 1; ipar < fDeltaFeed->GetNpar(); ipar++)
  //    fDeltaFeed->FixParameter(ipar, fDeltaFeed->GetParameter(ipar));
  //fDeltaFeed->FixParameter(1, fDeltaFeed->GetParameter(1));
  TProfile *pDeltaFeed_cent[NcentralityBins]; = new TProfile("pDeltaFeed", "", NptBins, ptBin);
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    if (!hDeltaFeed_cent[icent]) continue;

    /* restore MB parameters */
    for (Int_t ipar = 0; ipar < fDeltaFeed->GetNpar(); ipar++)
      fDeltaFeed->SetParameter(ipar, parMB[ipar]);

    if (fitFunc == kExpoPowerLaw) {
      fDeltaFeed->FixParameter(4, 0.);
      hDeltaFeed_mb->Fit(fDeltaFeed, "IMRE", "", 1.0, 3.0);
      fDeltaFeed->FixParameter(1, fDeltaFeed->GetParameter(1));
      fDeltaFeed->FixParameter(2, fDeltaFeed->GetParameter(2));
      fDeltaFeed->FixParameter(3, fDeltaFeed->GetParameter(3));
      fDeltaFeed->ReleaseParameter(4);
      hDeltaFeed_cent[icent]->Fit(fDeltaFeed, "IMRE", "", 0.5, 3.0);
    }
    else {
      hDeltaFeed_cent[icent]->Fit(fDeltaFeed, "IMRE", "", 0.5, 3.0);
    }

    pDeltaFeed_cent[icent] = new TProfile(Form("pDeltaFeed_cent%d", icent), "", NptBins, ptBin);
    HistoUtils_Function2Profile(fDeltaFeed, pDeltaFeed_cent[icent]);
    hCentDep->SetBinContent(icent + 1, fDeltaFeed->GetParameter(0.));
    hCentDep->SetBinError(icent + 1, fDeltaFeed->GetParError(0.));

    /* write */
    fileout->cd();
    hDeltaFeed_cent[icent]->Write(Form("hDeltaFeed_cent%d", icent));
    fDeltaFeed->Write(Form("fDeltaFeed_cent%d", icent));
    pDeltaFeed_cent[icent]->Write(Form("pDeltaFeed_cent%d", icent));
    
  }

  /* fit centrality dependence */
  TF1 *fCentDep = (TF1 *)gROOT->GetFunction("pol3");
  //  hCentDep->Fit(fCentDep);
  
  /* build fit profile */
  TProfile *pCentDep = new TProfile("pCentDep", "", NcentralityBins, centralityBin);
  HistoUtils_Function2Profile(fCentDep, pCentDep);
  
  /* write */
  fileout->cd();
  hCentDep->Write("hCentDep");
  fCentDep->Write("fCentDep");
  pCentDep->Write("pCentDep");
  fileout->Close();

}

DCA_primaryFraction(const Char_t *destdir)
{

  Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"
  };

  TFile *fileout = TFile::Open("TOF_primaryFraction.root", "RECREATE");

  TProfile *pFrac;
  TH1D *hFrac = new TH1D("hFrac", "", NptBins, ptBin);
  Char_t title[1024];
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
    if (ipart == 3) continue;
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      TFile *filein = TFile::Open(Form("%s/DCAdeltafeed_param_PrimaryCutFraction_%s_%s.root", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
      for (Int_t icent = -1; icent < NcentralityBins; icent++) {
	if (icent == -1)
	  pFrac = (TProfile *)filein->Get(Form("pDeltaFeed_MB"));
	else
	  pFrac = (TProfile *)filein->Get(Form("pDeltaFeed_cent%d", icent));
	/* copy profile in TH1D */
	for (Int_t ipt = 0; ipt < NptBins; ipt++) {
	  hFrac->SetBinContent(ipt + 1, pFrac->GetBinContent(ipt + 1));
	  hFrac->SetBinError(ipt + 1, pFrac->GetBinError(ipt + 1));
	}
	if (icent == -1)
	  sprintf(title, "%s (MB);p_{T} (GeV/c);fraction of primary particle;", partLatex[ipart][icharge]);
	else
	  sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);fraction of primary particle;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]);
	hFrac->SetMarkerStyle(marker[icharge]);
	hFrac->SetMarkerColor(color[ipart]);
	hFrac->SetTitle(title);
	if (icent == -1)
	  hFrac->SetName(Form("hPrimaryFrac_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	else
	  hFrac->SetName(Form("hPrimaryFrac_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	fileout->cd();
	hFrac->Write();
      }
    }
  }

  fileout->Close();

}
