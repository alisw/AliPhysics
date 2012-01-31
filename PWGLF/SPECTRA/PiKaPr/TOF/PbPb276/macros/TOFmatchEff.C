enum EHisto_t {
  kAcceptedTracks,
  kMatchedTracks,
  kMismatchedTracks,
  kUncorrelatedTracks,
  kMatchedCorrelatedTracks,
  kMatchedGoodTracks,
  kNHistos
};
const Char_t *histoName[kNHistos] = {
  "hAcceptedTracks",
  "hMatchedTracks",
  "hMismatchedTracks",
  "hUncorrelatedTracks",
  "hMatchedCorrelatedTracks",
  "hMatchedGoodTracks"
};

enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
const Char_t *chargeName[kNCharges] = {
  "positive",
  "negative"
};

enum EParam_t {
  kCentrality,
  kPt,
  kEta,
  kPhi,
  kNParams
};

const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};

const Int_t NptBins = 46;
Double_t ptBin[NptBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

const Int_t NetaBins = 20;
Double_t etaMin = -1.;
Double_t etaMax = 1.;
Double_t etaStep = (etaMax - etaMin) / NetaBins;
Double_t etaBin[NetaBins + 1]; /* computed at run-time */

const Int_t NphiBins = 200;
Double_t phiMin = 0.;
Double_t phiMax = 2. * TMath::Pi();
Double_t phiStep = (phiMax - phiMin) / NphiBins;
Double_t phiBin[NphiBins + 1]; /* computed at run-time */

const Int_t NptsubBins = 5;
Double_t ptsubBin[NptsubBins + 1] = {0.2, 0.5, 1.0, 1.5, 3.0, 5.0};
Int_t ptsubBinMin[NptsubBins] = {0, 6, 16, 21, 36};
Int_t ptsubBinMax[NptsubBins] = {5, 15, 20, 35, 45};

Int_t multcentColor[10] = {
  kRed,
  kOrange+1,
  kOrange,
  kYellow,
  kYellow+1,
  kGreen,
  kGreen+1,
  kCyan+1,
  kBlue,
  kMagenta,
  //  kMagenta+1  
};

const Char_t *extendedChargeName[2] = {"positive", "negative"};

const Double_t kEpsilon = 0.001;

TOFmatchEff(const Char_t *filename, const Char_t *enabledfilename = NULL, Int_t evMax = kMaxInt)
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
  if (enabledfilename) {
    TFile *enabledfile = TFile::Open(enabledfilename);
    hEnabledFlag = (TH1F *)enabledfile->Get("hEnabledFlag");
  }

  /* binning */
  for (Int_t ieta = 0; ieta < NetaBins + 1; ieta++)
    etaBin[ieta] = etaMin + ieta * etaStep;
  for (Int_t iphi = 0; iphi < NphiBins + 1; iphi++)
    phiBin[iphi] = phiMin + iphi * phiStep;
  /* THnSparse */
  Int_t NparamBins[kNParams] = {NcentralityBins, NptBins, NetaBins, NphiBins};
  Double_t *paramBin[kNParams] = {centralityBin, ptBin, etaBin, phiBin};
  THnSparseF *hHisto[kNHistos][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hHisto[ihisto][icharge] = new THnSparseF(Form("hHisto_%s_%s", histoName[ihisto], chargeName[icharge]), "", kNParams, NparamBins);
      for (Int_t iparam = 0; iparam < kNParams; iparam++)
	hHisto[ihisto][icharge]->SetBinEdges(iparam, paramBin[iparam]);
    }
#if 0
  TH2F *hHistoPID[kNHistos][AliPID::kSPECIES][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hHistoPID[ihisto][ipart][icharge] = new TH2F(Form("hHistoPID_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin, NptBins, ptBin);
      }
#endif
  /* histos */
  TH2F *hFEAMap = new TH2F("hFEAMap", "", 72, 0., 18., 91, 0., 91.);

  /* start stopwatch */
  TStopwatch timer;
  timer.Start();
  
  /* loop over events */
  Double_t param[kNParams];
  Int_t charge;
  Int_t index, sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  AliTOFcalibHisto calib;
  calib.LoadCalibHisto();
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 1000 == 0) printf("iev = %d\n", iev);
    /* check vertex */
    if (!analysisEvent->AcceptVertex()) continue;
    /* check collision candidate */
    if (!analysisEvent->IsCollisionCandidate()) continue;
    /* check centrality quality */
    if (analysisEvent->GetCentralityQuality() != 0.) continue;

    /*** ACCEPTED EVENT ***/

    /* get centrality */
    param[kCentrality] = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);
    
    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check accepted track */
      if (!analysisTrack->AcceptTrack()) continue;

      /*** ACCEPTED TRACK ***/
      
      /* get track info */
      param[kPt] = analysisTrack->GetPt();
      param[kEta] = analysisTrack->GetEta();
      param[kPhi] = analysisTrack->GetPhi();
      charge = analysisTrack->GetSign() > 0. ? kPositive : kNegative;
      index = analysisTrack->GetTOFIndex();
      
      /* fill accepted tracks histos */
      hHisto[kAcceptedTracks][charge]->Fill(param);
#if 0
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
	if (TMath::Abs(analysisTrack->GetY(AliPID::ParticleMass(ipart))) > 0.5) continue;
	hHistoPID[kAcceptedTracks][ipart][charge]->Fill(param[kCentrality], param[kPt]);
      }
#endif

      /* check TOF PID */
      if (!analysisTrack->HasTOFPID()) continue;
      /* check channel enabled */
      //      if (hEnabledFlag && hEnabledFlag->GetBinContent(index + 1) == 0.) continue;

      /*** ACCEPTED TRACK WITH TOF SIGNAL ***/

      /* fill FEA map */
      sector = calib.GetCalibMap(AliTOFcalibHisto::kSector, index);
      sectorStrip = calib.GetCalibMap(AliTOFcalibHisto::kSectorStrip, index);
      padx = calib.GetCalibMap(AliTOFcalibHisto::kPadX, index);
      fea = padx / 12;
      hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
      hitmapy = sectorStrip;
      hFEAMap->Fill(hitmapx, hitmapy);
      
      /* fill matched tracks histos */
      hHisto[kMatchedTracks][charge]->Fill(param);
      if (!analysisTrack->IsMismatchMC()) {
	hHisto[kMatchedGoodTracks][charge]->Fill(param);
	hHisto[kMatchedCorrelatedTracks][charge]->Fill(param);
      }
      else {
	hHisto[kMismatchedTracks][charge]->Fill(param);
	if (!analysisTrack->IsUncorrelatedMismatchMC())
	  hHisto[kMatchedCorrelatedTracks][charge]->Fill(param);
	else
	  hHisto[kUncorrelatedTracks][charge]->Fill(param);
      }
#if 0
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
	if (TMath::Abs(analysisTrack->GetY(AliPID::ParticleMass(ipart))) > 0.5) continue;
	hHistoPID[kMatchedTracks][ipart][charge]->Fill(param[kCentrality], param[kPt]);
	if (!analysisTrack->IsMismatchMC()) {
	  hHistoPID[kMatchedGoodTracks][charge]->Fill(param);
	  hHistoPID[kMatchedCorrelatedTracks][charge]->Fill(param);
	}
	else {
	  hHistoPID[kMismatchedTracks][charge]->Fill(param);
	  if (!analysisTrack->IsUncorrelatedMismatchMC())
	    hHistoPID[kMatchedCorrelatedTracks][charge]->Fill(param);
	  else
	    hHistoPID[kUncorrelatedTracks][charge]->Fill(param);
	}
      }
#endif
    }
  }

  /* stop stopwatch */
  timer.Stop();
  timer.Print();

  TFile *fileout = TFile::Open(Form("TOFmatchEff.%s", filename), "RECREATE");
  hFEAMap->Write();
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hHisto[ihisto][icharge]->Write();
    }
#if 0
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hHistoPID[ihisto][ipart][icharge]->Write();
      }
#endif
  fileout->Close();
  
}

//_____________________________________________________________________________-

TOFmatchEff_efficiencyPt(const Char_t *filename)
{

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos][kNCharges];
  TH2F *hHistoPID[kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoPt_MB[kNHistos][kNCharges], *hHistoPt_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hHistoPIDPt_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hHistoPIDPt_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoAllPt_MB[kNHistos], *hHistoAllPt_centrality[NcentralityBins][kNHistos];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {

    /* INCLUSIVE */

    hHistoAllPt_MB[ihisto] = new TH1D(Form("hHistoAllPt_MB_%s", histoName[ihisto]), "", NptBins, ptBin);
    for (Int_t icent = 0; icent < NcentralityBins; icent++)
      hHistoAllPt_centrality[icent][ihisto] = new TH1D(Form("hHistoAllPt_centrality%d_%s", icent, histoName[ihisto]), "", NptBins, ptBin);

    /* SINGLE PARTICLE */

    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      /* get histo */
      hHisto[ihisto][icharge] = (THnSparseF *)filein->Get(Form("hHisto_%s_%s", histoName[ihisto], chargeName[icharge]));
      
      /* MB projection */
      hHistoPt_MB[ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPt);
      hHistoPt_MB[ihisto][icharge]->SetName(Form("hHistoPt_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hHistoPt_MB[ihisto][icharge]->Sumw2();
      hHistoAllPt_MB[ihisto]->Add(hHistoPt_MB[ihisto][icharge]);
	
      /* centrality projection */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {

	hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	hHistoPt_centrality[icent][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPt);
	hHistoPt_centrality[icent][ihisto][icharge]->SetName(Form("hHistoPt_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hHistoPt_centrality[icent][ihisto][icharge]->Sumw2();
	hHistoAllPt_centrality[icent][ihisto]->Add(hHistoPt_centrality[icent][ihisto][icharge]);
      }
      
#if 0
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {

	/* get histo */
	hHistoPID[ihisto][ipart][icharge] = (TH2F *)filein->Get(Form("hHistoPID_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));

	/* MB projection */
	hHistoPIDPt_MB[ihisto][ipart][icharge] = hHistoPID[ihisto][ipart][icharge]->ProjectionY("hpy");
	hHistoPIDPt_MB[ihisto][ipart][icharge]->SetName(Form("hHistoPIDPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hHistoPIDPt_MB[ihisto][ipart][icharge]->Sumw2();
	
	/* centrality projection */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	  hHistoPIDPt_centrality[icent][ihisto][ipart][icharge] = hHistoPID[ihisto][ipart][icharge]->ProjectionY("hpy", icent + 1, icent + 1);
	  hHistoPIDPt_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hHistoPIDPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hHistoPIDPt_centrality[icent][ihisto][ipart][icharge]->Sumw2();

	}
	
      }
#endif
      
    }
  }

  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyPt");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiencyPt_MB[kNHistos][kNCharges], *hEfficiencyPt_centrality[NcentralityBins][kNHistos][kNCharges], *hEfficiencyPt_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hEfficiencyPIDPt_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyPIDPt_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyPIDPt_ratioMB_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hEfficiencyAllPt_MB[kNHistos], *hEfficiencyAllPt_centrality[NcentralityBins][kNHistos], *hEfficiencyAllPt_ratioMB_centrality[NcentralityBins][kNHistos];

  TH1D *hFractionPt_MB[kNHistos][kNCharges], *hFractionPt_centrality[NcentralityBins][kNHistos][kNCharges], *hFractionPt_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    
    if (ihisto == kAcceptedTracks) continue;
    
    /* INCLUSIVE */
    
    /* MB efficiency */
    hEfficiencyAllPt_MB[ihisto] = new TH1D(*hHistoAllPt_MB[ihisto]);
    hEfficiencyAllPt_MB[ihisto]->SetName(Form("hEfficiencyAllPt_MB_%s", histoName[ihisto]));
    hEfficiencyAllPt_MB[ihisto]->SetLineWidth(2);
    hEfficiencyAllPt_MB[ihisto]->SetLineColor(1);
    hEfficiencyAllPt_MB[ihisto]->SetMarkerStyle(20);
    hEfficiencyAllPt_MB[ihisto]->SetMarkerColor(1);
    hEfficiencyAllPt_MB[ihisto]->Divide(hEfficiencyAllPt_MB[ihisto], hHistoAllPt_MB[kAcceptedTracks], 1, 1, "B");
    hEfficiencyAllPt_MB[ihisto]->Write();
    
    /* multiplicity/centrality efficiency */
    for (Int_t icent = 0; icent < NcentralityBins; icent++) {
      hEfficiencyAllPt_centrality[icent][ihisto] = new TH1D(*hHistoAllPt_centrality[icent][ihisto]);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetName(Form("hEfficiencyAllPt_centrality%d_%s", icent, histoName[ihisto]));
      hEfficiencyAllPt_centrality[icent][ihisto]->SetLineWidth(2);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetLineColor(multcentColor[icent]);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetMarkerStyle(20);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetMarkerColor(multcentColor[icent]);
      hEfficiencyAllPt_centrality[icent][ihisto]->Divide(hEfficiencyAllPt_centrality[icent][ihisto], hHistoAllPt_centrality[icent][kAcceptedTracks], 1, 1, "B");
      hEfficiencyAllPt_centrality[icent][ihisto]->Write();
      
      /* ratio wrt. MB */
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto] = new TH1D(*hEfficiencyAllPt_centrality[icent][ihisto]);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->SetName(Form("hEfficiencyAllPt_ratioMB_centrality%d_%s", icent, histoName[ihisto]));
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->SetLineWidth(2);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->SetLineColor(multcentColor[icent]);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->SetMarkerStyle(20);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->SetMarkerColor(multcentColor[icent]);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->Divide(hEfficiencyAllPt_MB[ihisto]);
      hEfficiencyAllPt_ratioMB_centrality[icent][ihisto]->Write();
    }
    
    /* SINGLE PARTICLE */
    
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      /* MB efficiency */
      hEfficiencyPt_MB[ihisto][icharge] = new TH1D(*hHistoPt_MB[ihisto][icharge]);
      hEfficiencyPt_MB[ihisto][icharge]->SetName(Form("hEfficiencyPt_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hEfficiencyPt_MB[ihisto][icharge]->SetLineWidth(2);
      hEfficiencyPt_MB[ihisto][icharge]->SetLineColor(1);
      hEfficiencyPt_MB[ihisto][icharge]->SetMarkerStyle(20);
      hEfficiencyPt_MB[ihisto][icharge]->SetMarkerColor(1);
      hEfficiencyPt_MB[ihisto][icharge]->Divide(hEfficiencyPt_MB[ihisto][icharge], hHistoPt_MB[kAcceptedTracks][icharge], 1, 1, "B");
      hEfficiencyPt_MB[ihisto][icharge]->Write();
      
      /* multiplicity/centrality efficiency */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	hEfficiencyPt_centrality[icent][ihisto][icharge] = new TH1D(*hHistoPt_centrality[icent][ihisto][icharge]);
	hEfficiencyPt_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyPt_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyPt_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyPt_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyPt_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyPt_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyPt_centrality[icent][ihisto][icharge]->Divide(hEfficiencyPt_centrality[icent][ihisto][icharge], hHistoPt_centrality[icent][kAcceptedTracks][icharge], 1, 1, "B");
	hEfficiencyPt_centrality[icent][ihisto][icharge]->Write();
	
	/* ratio wrt. MB */
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge] = new TH1D(*hEfficiencyPt_centrality[icent][ihisto][icharge]);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyPt_ratioMB_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->Divide(hEfficiencyPt_MB[ihisto][icharge]);
	hEfficiencyPt_ratioMB_centrality[icent][ihisto][icharge]->Write();
      }
     
#if 0
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {

	/* MB efficiency */
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge] = new TH1D(*hHistoPIDPt_MB[ihisto][ipart][icharge]);
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->SetName(Form("hEfficiencyPIDPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->SetLineWidth(2);
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->SetLineColor(1);
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->SetMarkerStyle(20);
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->SetMarkerColor(1);
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->Divide(hEfficiencyPIDPt_MB[ihisto][ipart][icharge], hHistoPIDPt_MB[kAcceptedTracks][ipart][icharge], 1, 1, "B");
	hEfficiencyPIDPt_MB[ihisto][ipart][icharge]->Write();
	
	/* multiplicity/centrality efficiency */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {

	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hHistoPIDPt_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyPIDPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);

	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge], hHistoPIDPt_centrality[icent][kAcceptedTracks][ipart][icharge], 1, 1, "B");
	  hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]->Write();

	  /* ratio wrt. MB */
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hEfficiencyPIDPt_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyPIDPt_ratioMB_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyPIDPt_MB[ihisto][ipart][icharge]);
	  hEfficiencyPIDPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->Write();
	}
	
	
      }
#endif
      
    }       
  }
  
  fileout->Close();
	     
}

//_____________________________________________________________________________-

TOFmatchEff_efficiencyPhi(const Char_t *filename)
{

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos][kNCharges];
  TH1D *hHistoPhi_MB[kNHistos][kNCharges], *hHistoPhi_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hHistoPhi_MB_pt[NptsubBins][kNHistos][kNCharges], *hHistoPhi_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      /* get histo */
      hHisto[ihisto][icharge] = (THnSparseF *)filein->Get(Form("hHisto_%s_%s", histoName[ihisto], chargeName[icharge]));
      
      /* MB projection */
      hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(0, 0);
      hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
      hHistoPhi_MB[ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPhi);
      hHistoPhi_MB[ihisto][icharge]->SetName(Form("hHistoPhi_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hHistoPhi_MB[ihisto][icharge]->Sumw2();
      /* pt bins */
      for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
	hHistoPhi_MB_pt[ipt][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPhi);
	hHistoPhi_MB_pt[ipt][ihisto][icharge]->SetName(Form("hHistoPhi_MB_pt%d_%s_%s", ipt, histoName[ihisto], chargeName[icharge]));
	hHistoPhi_MB_pt[ipt][ihisto][icharge]->Sumw2();
      }
      
      /* centrality projection */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(0, 0);
	hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	hHistoPhi_centrality[icent][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPhi);
	hHistoPhi_centrality[icent][ihisto][icharge]->SetName(Form("hHistoPhi_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hHistoPhi_centrality[icent][ihisto][icharge]->Sumw2();
	/* pt bins */
	for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	  hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	  hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	  hHistoPhi_centrality_pt[icent][ipt][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kPhi);
	  hHistoPhi_centrality_pt[icent][ipt][ihisto][icharge]->SetName(Form("hHistoPhi_centrality%d_pt%d_%s_%s", icent, ipt, histoName[ihisto], chargeName[icharge]));
	  hHistoPhi_centrality_pt[icent][ipt][ihisto][icharge]->Sumw2();
	}
      }
    }
  
  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyPhi");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiencyPhi_MB[kNHistos][kNCharges], *hEfficiencyPhi_centrality[NcentralityBins][kNHistos][kNCharges], *hEfficiencyPhi_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hEfficiencyPhi_MB_pt[NptsubBins][kNHistos][kNCharges], *hEfficiencyPhi_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges], *hEfficiencyPhi_ratioMB_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges];


  TH1D *hFractionPhi_MB[kNHistos][kNCharges], *hFractionPhi_centrality[NcentralityBins][kNHistos][kNCharges], *hFractionPhi_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      if (ihisto == kAcceptedTracks) continue;
      
      /* MB efficiency */
      hEfficiencyPhi_MB[ihisto][icharge] = new TH1D(*hHistoPhi_MB[ihisto][icharge]);
      hEfficiencyPhi_MB[ihisto][icharge]->SetName(Form("hEfficiencyPhi_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hEfficiencyPhi_MB[ihisto][icharge]->SetLineWidth(2);
      hEfficiencyPhi_MB[ihisto][icharge]->SetLineColor(1);
      hEfficiencyPhi_MB[ihisto][icharge]->SetMarkerStyle(20);
      hEfficiencyPhi_MB[ihisto][icharge]->SetMarkerColor(1);
      hEfficiencyPhi_MB[ihisto][icharge]->Divide(hEfficiencyPhi_MB[ihisto][icharge], hHistoPhi_MB[kAcceptedTracks][icharge], 1, 1, "B");
      hEfficiencyPhi_MB[ihisto][icharge]->Write();
      /* pt bins */
      for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge] = new TH1D(*hHistoPhi_MB_pt[ipt][ihisto][icharge]);
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->SetName(Form("hEfficiencyPhi_MB_pt%d_%s_%s", ipt, histoName[ihisto], chargeName[icharge]));
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->SetLineColor(1);
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->SetMarkerColor(1);
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->Divide(hEfficiencyPhi_MB_pt[ipt][ihisto][icharge], hHistoPhi_MB_pt[ipt][kAcceptedTracks][icharge], 1, 1, "B");
	hEfficiencyPhi_MB_pt[ipt][ihisto][icharge]->Write();
      }
      
      /* multiplicity/centrality efficiency */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	hEfficiencyPhi_centrality[icent][ihisto][icharge] = new TH1D(*hHistoPhi_centrality[icent][ihisto][icharge]);
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyPhi_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->Divide(hEfficiencyPhi_centrality[icent][ihisto][icharge], hHistoPhi_centrality[icent][kAcceptedTracks][icharge], 1, 1, "B");
	hEfficiencyPhi_centrality[icent][ihisto][icharge]->Write();
	
	/* ratio wrt. MB */
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge] = new TH1D(*hEfficiencyPhi_centrality[icent][ihisto][icharge]);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyPhi_ratioMB_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->Divide(hEfficiencyPhi_MB[ihisto][icharge]);
	hEfficiencyPhi_ratioMB_centrality[icent][ihisto][icharge]->Write();
      }
      
    }       
  }
  
  fileout->Close();
	     
}

//_____________________________________________________________________________-

TOFmatchEff_efficiencyEta(const Char_t *filename)
{

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos][kNCharges];
  TH1D *hHistoEta_MB[kNHistos][kNCharges], *hHistoEta_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hHistoEta_MB_pt[NptsubBins][kNHistos][kNCharges], *hHistoEta_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      /* get histo */
      hHisto[ihisto][icharge] = (THnSparseF *)filein->Get(Form("hHisto_%s_%s", histoName[ihisto], chargeName[icharge]));
      
      /* MB projection */
      hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(0, 0);
      hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
      hHistoEta_MB[ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kEta);
      hHistoEta_MB[ihisto][icharge]->SetName(Form("hHistoEta_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hHistoEta_MB[ihisto][icharge]->Sumw2();
      /* pt bins */
      for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
	hHistoEta_MB_pt[ipt][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kEta);
	hHistoEta_MB_pt[ipt][ihisto][icharge]->SetName(Form("hHistoEta_MB_pt%d_%s_%s", ipt, histoName[ihisto], chargeName[icharge]));
	hHistoEta_MB_pt[ipt][ihisto][icharge]->Sumw2();
      }
      
      /* centrality projection */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(0, 0);
	hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	hHistoEta_centrality[icent][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kEta);
	hHistoEta_centrality[icent][ihisto][icharge]->SetName(Form("hHistoEta_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hHistoEta_centrality[icent][ihisto][icharge]->Sumw2();
	/* pt bins */
	for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	  hHisto[ihisto][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	  hHisto[ihisto][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	  hHistoEta_centrality_pt[icent][ipt][ihisto][icharge] = hHisto[ihisto][icharge]->Projection(kEta);
	  hHistoEta_centrality_pt[icent][ipt][ihisto][icharge]->SetName(Form("hHistoEta_centrality%d_pt%d_%s_%s", icent, ipt, histoName[ihisto], chargeName[icharge]));
	  hHistoEta_centrality_pt[icent][ipt][ihisto][icharge]->Sumw2();
	}
      }
    }
  
  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyEta");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiencyEta_MB[kNHistos][kNCharges], *hEfficiencyEta_centrality[NcentralityBins][kNHistos][kNCharges], *hEfficiencyEta_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  TH1D *hEfficiencyEta_MB_pt[NptsubBins][kNHistos][kNCharges], *hEfficiencyEta_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges], *hEfficiencyEta_ratioMB_centrality_pt[NcentralityBins][NptsubBins][kNHistos][kNCharges];


  TH1D *hFractionEta_MB[kNHistos][kNCharges], *hFractionEta_centrality[NcentralityBins][kNHistos][kNCharges], *hFractionEta_ratioMB_centrality[NcentralityBins][kNHistos][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      
      if (ihisto == kAcceptedTracks) continue;
      
      /* MB efficiency */
      hEfficiencyEta_MB[ihisto][icharge] = new TH1D(*hHistoEta_MB[ihisto][icharge]);
      hEfficiencyEta_MB[ihisto][icharge]->SetName(Form("hEfficiencyEta_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
      hEfficiencyEta_MB[ihisto][icharge]->SetLineWidth(2);
      hEfficiencyEta_MB[ihisto][icharge]->SetLineColor(1);
      hEfficiencyEta_MB[ihisto][icharge]->SetMarkerStyle(20);
      hEfficiencyEta_MB[ihisto][icharge]->SetMarkerColor(1);
      hEfficiencyEta_MB[ihisto][icharge]->Divide(hEfficiencyEta_MB[ihisto][icharge], hHistoEta_MB[kAcceptedTracks][icharge], 1, 1, "B");
      hEfficiencyEta_MB[ihisto][icharge]->Write();
      /* pt bins */
      for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge] = new TH1D(*hHistoEta_MB_pt[ipt][ihisto][icharge]);
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->SetName(Form("hEfficiencyEta_MB_pt%d_%s_%s", ipt, histoName[ihisto], chargeName[icharge]));
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->SetLineColor(1);
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->SetMarkerColor(1);
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->Divide(hEfficiencyEta_MB_pt[ipt][ihisto][icharge], hHistoEta_MB_pt[ipt][kAcceptedTracks][icharge], 1, 1, "B");
	hEfficiencyEta_MB_pt[ipt][ihisto][icharge]->Write();
      }
      
      /* multiplicity/centrality efficiency */
      for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	hEfficiencyEta_centrality[icent][ihisto][icharge] = new TH1D(*hHistoEta_centrality[icent][ihisto][icharge]);
	hEfficiencyEta_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyEta_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyEta_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyEta_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyEta_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyEta_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyEta_centrality[icent][ihisto][icharge]->Divide(hEfficiencyEta_centrality[icent][ihisto][icharge], hHistoEta_centrality[icent][kAcceptedTracks][icharge], 1, 1, "B");
	hEfficiencyEta_centrality[icent][ihisto][icharge]->Write();
	
	/* ratio wrt. MB */
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge] = new TH1D(*hEfficiencyEta_centrality[icent][ihisto][icharge]);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->SetName(Form("hEfficiencyEta_ratioMB_centrality%d_%s_%s", icent, histoName[ihisto], chargeName[icharge]));
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->SetLineWidth(2);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->SetLineColor(multcentColor[icent]);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerStyle(20);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->SetMarkerColor(multcentColor[icent]);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->Divide(hEfficiencyEta_MB[ihisto][icharge]);
	hEfficiencyEta_ratioMB_centrality[icent][ihisto][icharge]->Write();
      }
      
    }       
  }
  
  fileout->Close();
	     
}

#if 0
//_____________________________________________________________________________-

TOFmatchMC_efficiencyCent(const Char_t *filename, const Char_t *particle = "", const Char_t *charge = "", const Char_t *trdmode = "")
{

  const Int_t npt = 4;
  Double_t pt[npt + 1] = {0., 0.5, 1.0, 1.5, 5.0};

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos];
  TH1D *hHisto_all[kNHistos], *hHisto_MB_all[kNHistos], *hHisto_pt[kNHistos][MAXMULTCENTBINS], *hHisto_MB_pt[kNHistos][MAXMULTCENTBINS];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    if (ihisto == kMatchedCorrelatedTracks) continue;
    /* get histo */
    hHisto[ihisto] = (THnSparseF *)filein->Get(histoName[ihisto]);
    /* set range limits */
    hHisto[ihisto]->GetAxis(kRapidity)->SetRangeUser(-0.5 + kEpsilon, 0.5 - kEpsilon);
    hHisto[ihisto]->GetAxis(kEta)->SetRangeUser(-0.8 + kEpsilon, 0.8 - kEpsilon);
    hHisto[ihisto]->GetAxis(kPt)->SetRangeUser(0. + kEpsilon, 5.0 - kEpsilon);
    /* select particle if requested */
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      if (TString(particle) == AliPID::ParticleName(ipart))
	hHisto[ihisto]->GetAxis(kParticle)->SetRange(ipart + 1, ipart + 1);
    /* select charge if requested */
    if (TString(charge) == "plus")
      hHisto[ihisto]->GetAxis(kCharge)->SetRange(1, 1);
    else if (TString(charge) == "minus")
      hHisto[ihisto]->GetAxis(kCharge)->SetRange(2, 2);
    /* select TRD mode if requested */
    if (TString(trdmode) == "trdout")
      hHisto[ihisto]->GetAxis(kTRDmode)->SetRange(1, 1);
    else if (TString(trdmode) == "notrdout")
      hHisto[ihisto]->GetAxis(kTRDmode)->SetRange(2, 2);
    /* all projection */
    hHisto_all[ihisto] = hHisto[ihisto]->Projection(kMultCent);
    hHisto_all[ihisto]->SetName(Form("%s_all", histoName[ihisto]));
    hHisto_all[ihisto]->Sumw2();
    /* MB all projection */
    hHisto_MB_all[ihisto] = new TH1D(*hHisto_all[ihisto]);
    hHisto_MB_all[ihisto]->SetName(Form("%s_MB_all", histoName[ihisto]));
    for (Int_t i = 0; i < hHisto_MB_all[ihisto]->GetNbinsX(); i++) {
      hHisto_MB_all[ihisto]->SetBinContent(i + 1, hHisto_all[ihisto]->Integral());
      hHisto_MB_all[ihisto]->SetBinError(i + 1, TMath::Sqrt(hHisto_all[ihisto]->Integral()));
    }
    /* pt projection */
    for (Int_t ibin = 0; ibin < npt; ibin++) {
      hHisto[ihisto]->GetAxis(kPt)->SetRangeUser(pt[ibin] + kEpsilon, pt[ibin + 1] - kEpsilon);
      hHisto_pt[ihisto][ibin] = hHisto[ihisto]->Projection(kMultCent);
      hHisto_pt[ihisto][ibin]->SetName(Form("%s_pt%3.1f-%3.1f", histoName[ihisto], pt[ibin], pt[ibin + 1]));
      hHisto_pt[ihisto][ibin]->Sumw2();
      /* MB pt projection */
      hHisto_MB_pt[ihisto][ibin] = new TH1D(*hHisto_pt[ihisto][ibin]);
      hHisto_MB_pt[ihisto][ibin]->SetName(Form("%s_MB_pt%3.1f-%3.1f", histoName[ihisto], pt[ibin], pt[ibin + 1]));
      for (Int_t i = 0; i < hHisto_MB_pt[ihisto][ibin]->GetNbinsX(); i++) {
      	hHisto_MB_pt[ihisto][ibin]->SetBinContent(i + 1, hHisto_pt[ihisto][ibin]->Integral());
	hHisto_MB_pt[ihisto][ibin]->SetBinError(i + 1, TMath::Sqrt(hHisto_pt[ihisto][ibin]->Integral()));
      }
    }
  }
  /*** matched correlated histos ***/
  /* all projection */
  hHisto_all[kMatchedCorrelatedTracks] = new TH1D(*hHisto_all[kMatchedTracks]);
  hHisto_all[kMatchedCorrelatedTracks]->Add(hHisto_all[kUncorrelatedTracks], -1.);
  hHisto_all[kMatchedCorrelatedTracks]->SetName(Form("%s_all", histoName[kMatchedCorrelatedTracks]));
  /* MB all projection */
  hHisto_MB_all[kMatchedCorrelatedTracks] = new TH1D(*hHisto_MB_all[kMatchedTracks]);
  hHisto_MB_all[kMatchedCorrelatedTracks]->Add(hHisto_MB_all[kUncorrelatedTracks], -1.);
  hHisto_MB_all[kMatchedCorrelatedTracks]->SetName(Form("%s_MB_all", histoName[kMatchedCorrelatedTracks]));
  /* pt projection */
  for (Int_t ibin = 0; ibin < npt; ibin++) {
    hHisto_pt[kMatchedCorrelatedTracks][ibin] = new TH1D(*hHisto_pt[kMatchedTracks][ibin]);
    hHisto_pt[kMatchedCorrelatedTracks][ibin]->Add(hHisto_pt[kUncorrelatedTracks][ibin], -1.);
    hHisto_pt[kMatchedCorrelatedTracks][ibin]->SetName(Form("%s_pt%3.1f-%3.1f", histoName[kMatchedCorrelatedTracks], pt[ibin], pt[ibin + 1]));
    /* MB pt projection */
    hHisto_MB_pt[kMatchedCorrelatedTracks][ibin] = new TH1D(*hHisto_MB_pt[kMatchedTracks][ibin]);
    hHisto_MB_pt[kMatchedCorrelatedTracks][ibin]->Add(hHisto_MB_pt[kUncorrelatedTracks][ibin], -1.);
    hHisto_MB_pt[kMatchedCorrelatedTracks][ibin]->SetName(Form("%s_MB_pt%3.1f-%3.1f", histoName[kMatchedCorrelatedTracks], pt[ibin], pt[ibin + 1]));
  }

  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyCent");
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    if (TString(particle) == AliPID::ParticleName(ipart))
      str.Insert(str.Length() - TString(".root").Length(), Form(".%s", AliPID::ParticleName(ipart)));
  if (TString(charge) == "plus")
    str.Insert(str.Length() - TString(".root").Length(), ".plus");
  else if (TString(charge) == "minus")
    str.Insert(str.Length() - TString(".root").Length(), ".minus");
  if (TString(trdmode) == "trdout")
    str.Insert(str.Length() - TString(".root").Length(), ".trdout");
  else if (TString(trdmode) == "notrdout")
    str.Insert(str.Length() - TString(".root").Length(), ".notrdout");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiency_all[kNHistos], *hEfficiency_MB_all[kNHistos], *hEfficiency_ratioMB_all[kNHistos], *hEfficiency_pt[kNHistos][MAXMULTCENTBINS], *hEfficiency_MB_pt[kNHistos][MAXMULTCENTBINS], *hEfficiency_ratioMB_pt[kNHistos][MAXMULTCENTBINS];


  //  TH1D *hEfficiency_MB[kNHistos], *hEfficiency_multcent[kNHistos][MAXMULTCENTBINS], *hEfficiency_ratioMB_multcent[kNHistos][MAXMULTCENTBINS];
  //  TH1D *hFraction_MB[kNHistos], *hFraction_multcent[kNHistos][MAXMULTCENTBINS], *hFraction_ratioMB_multcent[kNHistos][MAXMULTCENTBINS];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    if (ihisto == kAcceptedTracks) continue;
    /* all efficiency */
    hEfficiency_all[ihisto] = new TH1D(*hHisto_all[ihisto]);
    hEfficiency_all[ihisto]->SetName(Form("%s_efficiency_all", histoName[ihisto]));
    hEfficiency_all[ihisto]->SetLineWidth(2);
    hEfficiency_all[ihisto]->SetLineColor(1);
    hEfficiency_all[ihisto]->SetMarkerStyle(20);
    hEfficiency_all[ihisto]->SetMarkerColor(1);
    hEfficiency_all[ihisto]->Divide(hHisto_all[kAcceptedTracks]);
    hEfficiency_all[ihisto]->Write();
    /* MB all efficiency */
    hEfficiency_MB_all[ihisto] = new TH1D(*hHisto_MB_all[ihisto]);
    hEfficiency_MB_all[ihisto]->SetName(Form("%s_efficiency_MB_all", histoName[ihisto]));
    hEfficiency_MB_all[ihisto]->Divide(hHisto_MB_all[kAcceptedTracks]);
    hEfficiency_MB_all[ihisto]->Write();
    /* ratio wrt. MB */
    hEfficiency_ratioMB_all[ihisto] = new TH1D(*hEfficiency_all[ihisto]);
    hEfficiency_ratioMB_all[ihisto]->SetName(Form("%s_efficiency_ratioMB_all", histoName[ihisto]));
    hEfficiency_ratioMB_all[ihisto]->Divide(hEfficiency_MB_all[ihisto]);
    hEfficiency_ratioMB_all[ihisto]->Write();
    /* pt efficiency */
    for (Int_t ibin = 0; ibin < npt; ibin++) {
      hEfficiency_pt[ihisto][ibin] = new TH1D(*hHisto_pt[ihisto][ibin]);
      hEfficiency_pt[ihisto][ibin]->SetName(Form("%s_efficiency_pt%3.1f-%3.1f", histoName[ihisto], pt[ibin], pt[ibin + 1]));
      hEfficiency_pt[ihisto][ibin]->SetLineWidth(2);
      hEfficiency_pt[ihisto][ibin]->SetLineColor(multcentColor[ibin]);
      hEfficiency_pt[ihisto][ibin]->SetMarkerStyle(20);
      hEfficiency_pt[ihisto][ibin]->SetMarkerColor(multcentColor[ibin]);
      hEfficiency_pt[ihisto][ibin]->Divide(hHisto_pt[kAcceptedTracks][ibin]);
      hEfficiency_pt[ihisto][ibin]->Write();
      /* MB pt efficiency */
      hEfficiency_MB_pt[ihisto][ibin] = new TH1D(*hHisto_MB_pt[ihisto][ibin]);
      hEfficiency_MB_pt[ihisto][ibin]->SetName(Form("%s_efficiency_MB_pt%3.1f-%3.1f", histoName[ihisto], pt[ibin], pt[ibin + 1]));
      hEfficiency_MB_pt[ihisto][ibin]->Divide(hHisto_MB_pt[kAcceptedTracks][ibin]);
      hEfficiency_MB_pt[ihisto][ibin]->Write();
      /* ratio wrt. central */
      hEfficiency_ratioMB_pt[ihisto][ibin] = new TH1D(*hEfficiency_pt[ihisto][ibin]);
      hEfficiency_ratioMB_pt[ihisto][ibin]->SetName(Form("%s_efficiency_ratioMB_pt%3.1f-%3.1f", histoName[ihisto], pt[ibin], pt[ibin + 1]));
      hEfficiency_ratioMB_pt[ihisto][ibin]->Divide(hEfficiency_MB_pt[ihisto][ibin]);
      hEfficiency_ratioMB_pt[ihisto][ibin]->Write();
    }

#if 0

    if (ihisto == kAcceptedTracks || ihisto == kMatchedTracks) continue;
    /* MB fraction */
    hFraction_MB[ihisto] = new TH1D(*hHisto_MB[ihisto]);
    hFraction_MB[ihisto]->SetName(Form("%s_fraction_MB", histoName[ihisto]));
    hFraction_MB[ihisto]->SetLineWidth(2);
    hFraction_MB[ihisto]->SetLineColor(1);
    hFraction_MB[ihisto]->SetMarkerStyle(20);
    hFraction_MB[ihisto]->SetMarkerColor(1);
    hFraction_MB[ihisto]->Divide(hHisto_MB[kMatchedTracks]);
    hFraction_MB[ihisto]->Write();
    /* multiplicity/centrality fraction */
    for (Int_t ibin = 0; ibin < hHisto[kMatchedTracks]->GetAxis(kMultCent)->GetNbins(); ibin++) {
       hFraction_multcent[ihisto][ibin] = new TH1D(*hHisto_multcent[ihisto][ibin]);
       hFraction_multcent[ihisto][ibin]->SetName(Form("%s_fraction_multcent%d", histoName[ihisto], ibin));
       hFraction_multcent[ihisto][ibin]->SetLineWidth(2);
       hFraction_multcent[ihisto][ibin]->SetLineColor(multcentColor[ibin]);
       hFraction_multcent[ihisto][ibin]->SetMarkerStyle(20);
       hFraction_multcent[ihisto][ibin]->SetMarkerColor(multcentColor[ibin]);
       hFraction_multcent[ihisto][ibin]->Divide(hHisto_multcent[kMatchedTracks][ibin]);
       hFraction_multcent[ihisto][ibin]->Write();
       /* ratio wrt. MB */
       hFraction_ratioMB_multcent[ihisto][ibin] = new TH1D(*hFraction_multcent[ihisto][ibin]);
       hFraction_ratioMB_multcent[ihisto][ibin]->SetName(Form("%s_fraction_ratioMB_multcent%d", histoName[ihisto], ibin));
       hFraction_ratioMB_multcent[ihisto][ibin]->SetLineWidth(2);
       hFraction_ratioMB_multcent[ihisto][ibin]->SetLineColor(multcentColor[ibin]);
       hFraction_ratioMB_multcent[ihisto][ibin]->SetMarkerStyle(20);
       hFraction_ratioMB_multcent[ihisto][ibin]->SetMarkerColor(multcentColor[ibin]);
       hFraction_ratioMB_multcent[ihisto][ibin]->Divide(hFraction_MB[ihisto]);
       hFraction_ratioMB_multcent[ihisto][ibin]->Write();
    }

    #endif
 
  }
  fileout->Close();
  
}
#endif

//_____________________________________________________________________________-

TH1D *
TOFmatchMC_get(const Char_t *filename, const Char_t *name, Option_t *opt = "", Int_t color = 1, Int_t marker = 20)
{

  TFile *f = TFile::Open(filename);
  if (!f || !f->IsOpen())
    return NULL;
  TH1D *h = (TH1D *)f->Get(name);
  if (!h)
    return NULL;
  h->SetLineColor(color);
  h->SetLineWidth(2);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->Draw(opt);
  return h;
}

//_____________________________________________________________________________-

TH1D *
TOFmatchMC_divide(const Char_t *filename, const Char_t *name1, const Char_t *name2, Option_t *opt = "")
{

  TFile *f = TFile::Open(filename);
  if (!f || !f->IsOpen())
    return NULL;
  TH1D *h1 = (TH1D *)f->Get(name1);
  TH1D *h2 = (TH1D *)f->Get(name2);
  if (!h1 || !h2)
    return NULL;
  TH1D *hr = new TH1D(*h1);
  hr->Divide(h2);
  hr->Draw(opt);
  return hr;
}

//_____________________________________________________________________________-

TH1D *
TOFmatchMC_sub(const Char_t *filename, const Char_t *name1, const Char_t *name2, Option_t *opt = "", Int_t color = 1, Int_t marker = 20)
{

  TFile *f = TFile::Open(filename);
  if (!f || !f->IsOpen())
    return NULL;
  TH1D *h1 = (TH1D *)f->Get(name1);
  TH1D *h2 = (TH1D *)f->Get(name2);
  if (!h1 || !h2)
    return NULL;
  TH1D *hr = new TH1D(*h1);
  hr->Add(h2, -1.);
  hr->SetLineColor(color);
  hr->SetLineWidth(2);
  hr->SetMarkerColor(color);
  hr->SetMarkerStyle(marker);
  hr->Draw(opt);
  return hr;
}

//_____________________________________________________________________________-

TH1D *
TOFmatchMC_compare(const Char_t *filename1, const Char_t *filename2, const Char_t *name, Option_t *opt = "", Int_t color = 1, Int_t marker = 20)
{

  TFile *f1 = TFile::Open(filename1);
  TFile *f2 = TFile::Open(filename2);
  if (!f1 || !f2 || !f1->IsOpen() || !f2->IsOpen())
    return NULL;
  TH1D *h1 = (TH1D *)f1->Get(name);
  TH1D *h2 = (TH1D *)f2->Get(name);
  if (!h1 || !h2)
    return NULL;
  TH1D *hr = new TH1D(*h1);
  hr->Divide(h2);
  hr->SetLineColor(color);
  hr->SetLineWidth(2);
  hr->SetMarkerColor(color);
  hr->SetMarkerStyle(marker);
  hr->Draw(opt);
  return hr;
}

//_____________________________________________________________________________-

TH1D *
TOFmatchMC_comparesub(const Char_t *filename1, const Char_t *filename2, const Char_t *name1, const Char_t *name2, Option_t *opt = "", Int_t color = 1, Int_t marker = 20)
{

  TFile *f1 = TFile::Open(filename1);
  TFile *f2 = TFile::Open(filename2);
  if (!f1 || !f2 || !f1->IsOpen() || !f2->IsOpen())
    return NULL;
  TH1D *h11 = (TH1D *)f1->Get(name1);
  TH1D *h21 = (TH1D *)f2->Get(name1);
  TH1D *h12 = (TH1D *)f1->Get(name2);
  TH1D *h22 = (TH1D *)f2->Get(name2);
  if (!h11 || !h21 || !h12 || !h22)
    return NULL;
  TH1D *hs1 = new TH1D(*h11);
  hs1->Add(h12, -1.);
  TH1D *hs2 = new TH1D(*h21);
  hs2->Add(h22, -1.);
  TH1D *hr = new TH1D(*hs1);
  hr->Divide(hs2);
  hr->SetLineColor(color);
  hr->SetLineWidth(2);
  hr->SetMarkerColor(color);
  hr->SetMarkerStyle(marker);
  hr->Draw(opt);
  return hr;
}


TH1D *
TOFmatchEff_efficiencyPt_MB_plot(const Char_t *filename, Int_t ihisto = kMatchedTracks, Int_t icharge, Int_t marker = 20, Int_t color = 2, Option_t *opt = "")
{
  
  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  const Char_t *destdir = "matchingEfficiency_DATA";
  
  Double_t ptMin = 0.5;
  Double_t ptMax = 5.0;

  TF1 *fEff = new TF1("fEff", "[0] + [1] * x - [2] * TMath::Exp(-[3] * TMath::Power(x, [4]))", 0., 5.0);
  fEff->SetParameter(0, 0.5);
  fEff->SetParameter(1, 0.);
  fEff->SetParameter(2, 0.5);
  fEff->SetParameter(3, 1.);
  fEff->SetParameter(4, 2.);


  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyPt_MB_%s_%s", histoName[ihisto], chargeName[icharge]));
  hEfficiencyPt->Fit(fEff, "0", "", 0.5, 5.0);
  hEfficiencyPt->SetTitle(Form("%s tracks;p_{T} (GeV/c);acceptance #times efficiency;", extendedChargeName[icharge]));
  hEfficiencyPt->SetMinimum(0.4);
  hEfficiencyPt->SetMaximum(0.8);
  hEfficiencyPt->SetMarkerStyle(marker);
  hEfficiencyPt->SetMarkerColor(color);
  hEfficiencyPt->SetMarkerSize(1.5);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->SetLineColor(color);
  hEfficiencyPt->GetXaxis()->SetRangeUser(ptMin + 0.001, ptMax - 0.001);
  hEfficiencyPt->SetStats(kFALSE);
  cCanvas1->cd();
  hEfficiencyPt->Draw(opt);
  fEff->Draw("same");

  cCanvas1->SetGridx();
  cCanvas1->SetGridy();
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s.C", destdir, chargeName[icharge]));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s.root", destdir, chargeName[icharge]));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s.png", destdir, chargeName[icharge]));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s.eps", destdir, chargeName[icharge]));
  
  TH1D *hRatioPt = new TH1D(*hEfficiencyPt);
  hRatioPt->Divide(fEff);
  hRatioPt->SetTitle(Form("%s tracks;p_{T} (GeV/c);ratio wrt. fitted dependence;", extendedChargeName[icharge]));
  hRatioPt->SetMinimum(0.9);
  hRatioPt->SetMaximum(1.1);
  cCanvas2->cd();
  hRatioPt->Draw();

  cCanvas2->SetGridx();
  cCanvas2->SetGridy();
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s.C", destdir, chargeName[icharge]));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s.root", destdir, chargeName[icharge]));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s.png", destdir, chargeName[icharge]));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s.eps", destdir, chargeName[icharge]));
  

  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyPt;
}


//_____________________________________________________________________________-

TH1D *
TOFmatchEff_efficiencyPt_centrality_all_plot(const Char_t *filename, Int_t ihisto = kMatchedTracks, Int_t marker = 20, Int_t color = 1, Option_t *opt = "")
{
  
  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  const Char_t *destdir = "matchingEfficiency_DATA";
  
  Double_t ptMin = 0.5;
  Double_t ptMax = 5.0;

  TF1 *fMismatchFrac = new TF1("fMismatchFrac", "([0]+[1]*TMath::Exp(-[2]*TMath::Power(x,[3])))*[4]", 0., 5.0);
  fMismatchFrac->SetParameter(0, 0.0447133);
  fMismatchFrac->SetParameter(1, 0.179172);
  fMismatchFrac->SetParameter(2, 2.54333);
  fMismatchFrac->SetParameter(3, 1.16819);
  fMismatchFrac->SetParameter(4, 1.);
  
  TF1 *fEff = new TF1("fEff", "([0] + [1] * x - [2] * TMath::Exp(-[3] * TMath::Power(x, [4]))) * [5]", 0., 5.0);
  fEff->SetParameter(0, 0.5);
  fEff->SetParameter(1, 0.);
  fEff->SetParameter(2, 0.5);
  fEff->SetParameter(3, 1.);
  fEff->SetParameter(4, 2.);
  fEff->FixParameter(5, 1.);
  
  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyAllPt_MB_%s", histoName[ihisto]));
  hEfficiencyPt->Fit(fEff, "0", "", 0.5, 5.0);
  hEfficiencyPt->SetTitle("all particles;p_{T} (GeV/c);acceptance #times efficiency;");
  hEfficiencyPt->SetMinimum(0.4);
  hEfficiencyPt->SetMaximum(0.8);
  hEfficiencyPt->SetMarkerStyle(marker);
  hEfficiencyPt->SetMarkerColor(color);
  hEfficiencyPt->SetMarkerSize(1.5);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->SetLineColor(color);
  hEfficiencyPt->GetXaxis()->SetRangeUser(0.5 + 0.001, 5.0 - 0.001);
  hEfficiencyPt->SetStats(kFALSE);
  cCanvas1->cd();
  hEfficiencyPt->Draw(opt);
  fEff->Draw("same");

  cCanvas1->SetGridx();
  cCanvas1->SetGridy();
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.C", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.root", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.png", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.eps", destdir));

  TH1D *hRatioPt = new TH1D(*hEfficiencyPt);
  hRatioPt->Divide(fEff);
  hRatioPt->SetTitle("all particles;p_{T} (GeV/c);ratio wrt. fitted dependence;");
  hRatioPt->SetMinimum(0.9);
  hRatioPt->SetMaximum(1.1);
  cCanvas2->cd();
  hRatioPt->Draw();

  cCanvas2->SetGridx();
  cCanvas2->SetGridy();
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_all.C", destdir));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_all.root", destdir));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_all.png", destdir));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_all.eps", destdir));
  
  /* fix efficiency shape and release scale factor */
  fEff->FixParameter(0, fEff->GetParameter(0));
  fEff->FixParameter(1, fEff->GetParameter(1));
  fEff->FixParameter(2, fEff->GetParameter(2));
  fEff->FixParameter(3, fEff->GetParameter(3));
  fEff->FixParameter(4, fEff->GetParameter(4));
  fEff->ReleaseParameter(5);
  
  TH1D *hEfficiencyCent = new TH1D("hEfficiencyCent", "all particles;centrality percentile;acceptance x efficiency scale factor;", NcentralityBins, centralityBin);
  hEfficiencyCent->SetMinimum(0.95);
  hEfficiencyCent->SetMaximum(1.05);
  hEfficiencyCent->SetMarkerStyle(marker);
  hEfficiencyCent->SetMarkerColor(color);
  hEfficiencyCent->SetMarkerSize(1.5);
  hEfficiencyCent->SetLineWidth(2);
  hEfficiencyCent->SetLineColor(color);
  hEfficiencyCent->SetStats(kFALSE);
  
  
  TH1D *hEfficiencyPt_cent[NcentralityBins];
  TH1D *hRatioPt_cent[NcentralityBins];
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    
  
    hEfficiencyPt_cent[icent] = (TH1D *)filein->Get(Form("hEfficiencyAllPt_centrality%d_%s", icent, histoName[ihisto]));
    hEfficiencyPt_cent[icent]->Fit(fEff, "", "", 0.5, 5.0);
    
    hEfficiencyPt_cent[icent]->SetTitle(Form("all particles (%d-%d\%);p_{T} (GeV/c);acceptance #times efficiency;", (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hEfficiencyPt_cent[icent]->SetMinimum(0.2);
    hEfficiencyPt_cent[icent]->SetMaximum(0.8);
    hEfficiencyPt_cent[icent]->SetMarkerStyle(marker);
    hEfficiencyPt_cent[icent]->SetMarkerColor(color);
    hEfficiencyPt_cent[icent]->SetMarkerSize(1.5);
    hEfficiencyPt_cent[icent]->SetLineWidth(2);
    hEfficiencyPt_cent[icent]->SetLineColor(color);
    hEfficiencyPt_cent[icent]->GetXaxis()->SetRangeUser(0.5 + 0.001, 5.0 - 0.001);
    hEfficiencyPt_cent[icent]->SetStats(kFALSE);
    cCanvas1->cd();
    hEfficiencyPt_cent[icent]->Draw(opt);
    fEff->Draw("same");
    
    hEfficiencyCent->SetBinContent(icent + 1, fEff->GetParameter(5));
    hEfficiencyCent->SetBinError(icent + 1, fEff->GetParError(5));
    
    cCanvas1->SetGridx();
    cCanvas1->SetGridy();
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.C", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.root", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.png", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.eps", destdir, icent));
    
    hRatioPt_cent[icent] = new TH1D(*hEfficiencyPt_cent[icent]);
    hRatioPt_cent[icent]->Divide(fEff);
    hRatioPt_cent[icent]->SetTitle(Form("all particles (%d-%d\%);p_{T} (GeV/c);ratio wrt. fitted dependence;", (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hRatioPt_cent[icent]->SetMinimum(0.9);
    hRatioPt_cent[icent]->SetMaximum(1.1);
    cCanvas2->cd();
    hRatioPt_cent[icent]->Draw();

    cCanvas2->SetGridx();
    cCanvas2->SetGridy();
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.C", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.root", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.png", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.eps", destdir, icent));

    
  }

  TF1 *fEffCent = new TF1("fEffCent", "[0] - [1] * TMath::Exp(-[2] * TMath::Power(x, [3]))", 0., 90.);
  fEffCent->SetParameter(0, 1.02);
  fEffCent->SetParameter(1, 0.04);
  fEffCent->SetParameter(2, 0.001);
  fEffCent->SetParameter(3, 2.);
  //  hEfficiencyCent->Fit(fEffCent, "q0", "", 0., 90.);
  
  TCanvas *cCanvas3 = new TCanvas("cCanvas3");
  hEfficiencyCent->Draw();
  //  fEffCent->Draw("same");
  
  cCanvas3->SetGridx();
  cCanvas3->SetGridy();
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.C", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.root", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.png", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.eps", destdir));

  TCanvas *cCanvas4 = new TCanvas("cCanvas4");

  TH1D *hRatioCent = new TH1D(*hEfficiencyCent);
  hRatioCent->Divide(fEffCent);
  hRatioCent->SetTitle(Form("all particles;centrality percentile;ratio wrt. fitted dependence;"));
  hRatioCent->SetMinimum(0.95);
  hRatioCent->SetMaximum(1.05);
  cCanvas4->cd();
  hRatioCent->Draw();

  cCanvas4->SetGridx();
  cCanvas4->SetGridy();
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.C", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.root", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.png", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.eps", destdir));
  


  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyCent;
}


//_____________________________________________________________________________-

TH1D *
TOFmatchEff_efficiencyPt_centrality_all_plot_nomismatch(const Char_t *filename, Int_t ihisto = kMatchedTracks, Int_t marker = 20, Int_t color = 1, Option_t *opt = "")
{

  TF1 *fpol0 = (TF1 *)gROOT->GetFunction("pol0");
  fpol0->SetRange(0., 5.0);

  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  const Char_t *destdir = "matchingEfficiency_DATA";
  
  Double_t ptMin = 0.5;
  Double_t ptMax = 5.0;

  TF1 *fMismatchCorr = new TF1("fMismatchCorr", "1. - ([0]+[1]*TMath::Exp(-[2]*TMath::Power(x,[3])))*[4]", 0., 5.0);
  fMismatchCorr->SetParameter(0, 4.47133e-02);
  fMismatchCorr->SetParameter(1, 1.79172e-01);
  fMismatchCorr->SetParameter(2, 2.54333e+00);
  fMismatchCorr->SetParameter(3, 1.16819e+00);
  fMismatchCorr->SetParameter(4, 1.);

  TF1 *fMismatchScale = new TF1("fMismatchScale", "[0] + [1] * TMath::Exp(-[2] * x)", 0., 100.);
  fMismatchScale->SetParameter(0, -6.36877e-02);
  fMismatchScale->SetParameter(1, 1.74818e+00);
  fMismatchScale->SetParameter(2, 3.00818e-02);
  
  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyAllPt_MB_%s", histoName[ihisto]));
  if (ihisto != kMatchedCorrelatedTracks)
    hEfficiencyPt->Multiply(fMismatchCorr);
  hEfficiencyPt->SetTitle("all particles;p_{T} (GeV/c);acceptance #times efficiency;");
  hEfficiencyPt->SetMinimum(0.4);
  hEfficiencyPt->SetMaximum(0.8);
  hEfficiencyPt->SetMarkerStyle(marker);
  hEfficiencyPt->SetMarkerColor(color);
  hEfficiencyPt->SetMarkerSize(1.5);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->SetLineColor(color);
  hEfficiencyPt->GetXaxis()->SetRangeUser(0.5 + 0.001, 5.0 - 0.001);
  hEfficiencyPt->SetStats(kFALSE);
  cCanvas1->cd();
  hEfficiencyPt->Draw(opt);

  cCanvas1->SetGridx();
  cCanvas1->SetGridy();
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.C", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.root", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.png", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.eps", destdir));

  TH1D *hEfficiencyCent = new TH1D("hEfficiencyCent", "all particles;centrality percentile;acceptance x efficiency scale factor;", NcentralityBins, centralityBin);
  hEfficiencyCent->SetMinimum(0.95);
  hEfficiencyCent->SetMaximum(1.05);
  hEfficiencyCent->SetMarkerStyle(marker);
  hEfficiencyCent->SetMarkerColor(color);
  hEfficiencyCent->SetMarkerSize(1.5);
  hEfficiencyCent->SetLineWidth(2);
  hEfficiencyCent->SetLineColor(color);
  hEfficiencyCent->SetStats(kFALSE);
  
  
  TH1D *hEfficiencyPt_cent[NcentralityBins];
  TH1D *hRatioPt_cent[NcentralityBins];
  Double_t centmean;
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {
    
    centmean = 0.5 * (centralityBin[icent] + centralityBin[icent + 1]);
    hEfficiencyPt_cent[icent] = (TH1D *)filein->Get(Form("hEfficiencyAllPt_centrality%d_%s", icent, histoName[ihisto]));
    fMismatchCorr->SetParameter(4, fMismatchScale->Eval(centmean));
    if (ihisto != kMatchedCorrelatedTracks)
      hEfficiencyPt_cent[icent]->Multiply(fMismatchCorr);
    hEfficiencyPt_cent[icent]->Divide(hEfficiencyPt);
    hEfficiencyPt_cent[icent]->Fit(fpol0, "q0", "", 0.5, 5.0);
    hEfficiencyPt_cent[icent]->SetTitle(Form("all particles (%d-%d\%);p_{T} (GeV/c);acceptance #times efficiency;", (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hEfficiencyPt_cent[icent]->SetMinimum(0.2);
    hEfficiencyPt_cent[icent]->SetMaximum(0.8);
    hEfficiencyPt_cent[icent]->SetMarkerStyle(marker);
    hEfficiencyPt_cent[icent]->SetMarkerColor(color);
    hEfficiencyPt_cent[icent]->SetMarkerSize(1.5);
    hEfficiencyPt_cent[icent]->SetLineWidth(2);
    hEfficiencyPt_cent[icent]->SetLineColor(color);
    hEfficiencyPt_cent[icent]->GetXaxis()->SetRangeUser(0.5 + 0.001, 5.0 - 0.001);
    hEfficiencyPt_cent[icent]->SetStats(kFALSE);
    cCanvas1->cd();
    hEfficiencyPt_cent[icent]->Draw(opt);
    fpol0->Draw("same");
    
    hEfficiencyCent->SetBinContent(icent + 1, fpol0->GetParameter(0));
    hEfficiencyCent->SetBinError(icent + 1, fpol0->GetParError(0));
    
    cCanvas1->SetGridx();
    cCanvas1->SetGridy();
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.C", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.root", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.png", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.eps", destdir, icent));
    
    hRatioPt_cent[icent] = new TH1D(*hEfficiencyPt_cent[icent]);
    hRatioPt_cent[icent]->Divide(fpol0);
    hRatioPt_cent[icent]->SetTitle(Form("all particles (%d-%d\%);p_{T} (GeV/c);ratio wrt. fitted dependence;", (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hRatioPt_cent[icent]->SetMinimum(0.9);
    hRatioPt_cent[icent]->SetMaximum(1.1);
    cCanvas2->cd();
    hRatioPt_cent[icent]->Draw();

    cCanvas2->SetGridx();
    cCanvas2->SetGridy();
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.C", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.root", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.png", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.eps", destdir, icent));

    
  }

  TF1 *fEffCent = new TF1("fEffCent", "[0] - [1] * TMath::Exp(-[2] * TMath::Power(x, [3]))", 0., 90.);
  fEffCent->SetParameter(0, 1.02);
  fEffCent->SetParameter(1, 0.04);
  fEffCent->SetParameter(2, 0.001);
  fEffCent->SetParameter(3, 2.);
  //  hEfficiencyCent->Fit(fEffCent, "q0", "", 0., 90.);
  
  TCanvas *cCanvas3 = new TCanvas("cCanvas3");
  hEfficiencyCent->Draw();
  //  fEffCent->Draw("same");
  
  cCanvas3->SetGridx();
  cCanvas3->SetGridy();
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.C", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.root", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.png", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.eps", destdir));

  TCanvas *cCanvas4 = new TCanvas("cCanvas4");

  TH1D *hRatioCent = new TH1D(*hEfficiencyCent);
  hRatioCent->Divide(fEffCent);
  hRatioCent->SetTitle(Form("all particles;centrality percentile;ratio wrt. fitted dependence;"));
  hRatioCent->SetMinimum(0.95);
  hRatioCent->SetMaximum(1.05);
  cCanvas4->cd();
  hRatioCent->Draw();

  cCanvas4->SetGridx();
  cCanvas4->SetGridy();
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.C", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.root", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.png", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.eps", destdir));
  


  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyCent;
}


TH1D *
TOFmatchEff_rapidityCut(const Char_t *filename, Int_t ipart, Int_t icharge, Int_t color = 2, Int_t marker = 20, Option_t *opt = "")
{

  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiency_MB = (TH1D *)filein->Get(Form("hEfficiencyPt_MB_hMatchedTracks_%s", chargeName[icharge]));
  TH1D *hEfficiencyPID_MB = (TH1D *)filein->Get(Form("hEfficiencyPIDPt_MB_hMatchedTracks_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));

  TH1D *hRatio = new TH1D(*hEfficiencyPID_MB);
  hRatio->Divide(hEfficiency_MB);
  hRatio->SetLineColor(color);
  hRatio->SetLineWidth(2);
  hRatio->SetMarkerColor(color);
  hRatio->SetMarkerStyle(marker);
  hRatio->Draw(opt);
  return hRatio;
}

void
TOFmatchEff_checkFEAmap(const Char_t *filename)
{

  TFile *filein = TFile::Open(filename);
  TH2F *hFEAMap = (TH2F *)filein->Get("hFEAMap");
  TH2F *hEmptyFEA = new TH2F("hEmptyFEA", "", 72, 0., 18., 91, 0., 91.);

  Float_t nhits = hFEAMap->GetEntries();
  if (nhits <= 0) {
    printf("found not hits.\n");
    return;
  }
  printf("found %d hits\n", nhits);
  Float_t avhits = nhits / (hFEAMap->GetNbinsX() * hFEAMap->GetNbinsX());
  printf("on average %f hits/FEA\n", avhits);
  
  TH1F *hFEANormHit = new TH1F("hFEANormHit", "", 1000, 0., 3.);
  for (Int_t ibinx = 0; ibinx < hFEAMap->GetNbinsX(); ibinx++)
    for (Int_t ibiny = 0; ibiny < hFEAMap->GetNbinsY(); ibiny++) {
      if (hFEAMap->GetBinContent(ibinx + 1, ibiny + 1) == 0.) {
	hEmptyFEA->SetBinContent(ibinx + 1, ibiny + 1, 1.);
	continue;
      }
      hFEANormHit->Fill(hFEAMap->GetBinContent(ibinx + 1, ibiny + 1) / avhits);
    }

  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  cCanvas1->Divide(1, 2);
  cCanvas1->cd(1);
  hFEAMap->Draw("colz");
  cCanvas1->cd(2);
  hEmptyFEA->Draw("colz");

  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  hFEANormHit->Draw();


  TFile *fileout = TFile::Open(Form("TOFmatchEff_checkFEAmap.%s", filename), "RECREATE");
  hFEAMap->Write();
  hFEANormHit->Write();
  hEmptyFEA->Write();
  fileout->Close();
}
