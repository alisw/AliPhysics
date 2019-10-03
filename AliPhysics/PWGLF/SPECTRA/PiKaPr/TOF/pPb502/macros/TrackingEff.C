#include "CommonDefs.C"

const Char_t *destdir = "plots";

enum EHisto_t {
  kPrimaryTracks,
  kReconstructedTracks,
  kNHistos
};
const Char_t *histoName[kNHistos] = {
  "hPrimaryTracks",
  "hReconstructedTracks"
};

enum EParam_t {
  kCentrality,
  kPt,
  kEta,
  kPhi,
  kY,
  kNParams
};

TrackingEff(const Char_t *filename, Int_t evMax = kMaxInt)
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
  TClonesArray *analysisParticleArray = new TClonesArray("AliAnalysisParticle");
  AliAnalysisParticle *analysisParticle = NULL;
  treein->SetBranchAddress("AnalysisEvent", &analysisEvent);
  treein->SetBranchAddress("AnalysisTrack", &analysisTrackArray);
  treein->SetBranchAddress("AnalysisParticle", &analysisParticleArray);

  /* binning */
  for (Int_t iy = 0; iy < NyBins + 1; iy++)
    yBin[iy] = yMin + iy * yStep;
  for (Int_t ieta = 0; ieta < NetaBins + 1; ieta++)
    etaBin[ieta] = etaMin + ieta * etaStep;
  for (Int_t iphi = 0; iphi < NphiBins + 1; iphi++)
    phiBin[iphi] = phiMin + iphi * phiStep;
  /* THnSparse */
  Int_t NparamBins[kNParams] = {NcentralityBins, NptBins, NetaBins, NphiBins, NyBins};
  Double_t *paramBin[kNParams] = {centralityBin, ptBin, etaBin, phiBin, yBin};
  THnSparseF *hHisto[kNHistos][AliPID::kSPECIES][kNCharges];
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	hHisto[ihisto][ipart][icharge] = new THnSparseF(Form("hHisto_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), "", kNParams, NparamBins);
	for (Int_t iparam = 0; iparam < kNParams; iparam++)
	  hHisto[ihisto][ipart][icharge]->SetBinEdges(iparam, paramBin[iparam]);
      }

  /* selected particle array */
  AliAnalysisParticle *selParticle[1000000];
  for (Int_t idx = 0; idx < 1000000; idx++)
    selParticle[idx] = NULL;

  /* filled slot array */
  Int_t filledSlot[1000000];
  Int_t nFilledSlots = 0;


  /* loop over events */
  Double_t param[kNParams];
  Int_t part, charge;
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
    param[kCentrality] = analysisEvent->GetCentralityPercentile(centralityEstimator);

    /*** ACCEPTED EVENT ***/
    
    /* reset filled slots selected particle array */
    for (Int_t islot = 0; islot < nFilledSlots; islot++)
      selParticle[filledSlot[islot]] = NULL;
    /* reset filled slots */
    nFilledSlots = 0;

    /* loop over primaries */
    for (Int_t iprim = 0; iprim < analysisParticleArray->GetEntries(); iprim++) {
      /* get and check particle */
      analysisParticle = (AliAnalysisParticle *)analysisParticleArray->At(iprim);
      if (!analysisParticle || analysisParticle->GetPID() == -1 || analysisParticle->GetLabel() < 0 || analysisParticle->GetSign() == 0.) continue;

      /* check rapidity */
      if (analysisParticle->GetY() - rapidityShift > rapidityMaxCut ||
	  analysisParticle->GetY() - rapidityShift < rapidityMinCut)
	continue;

      /*** ACCEPTED PARTICLE ***/

      /* get particle info */
      part = analysisParticle->GetPID();
      param[kPt] = analysisParticle->GetPt();
      param[kEta] = analysisParticle->GetEta();
      param[kPhi] = analysisParticle->GetPhi();
      charge = analysisParticle->GetSign() > 0. ? kPositive : kNegative;

      /* fill histo and slot */
      hHisto[kPrimaryTracks][part][charge]->Fill(param);
      selParticle[analysisParticle->GetLabel()] = analysisParticle;
      filledSlot[nFilledSlots] = analysisParticle->GetLabel();
      nFilledSlots++;
      
    }

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get and check track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check accepted track */
      if (!analysisTrack->AcceptTrack()) continue;
      /* get corresponding particle */
      analysisParticle = selParticle[TMath::Abs(analysisTrack->GetLabel())];
      if (!analysisParticle) continue;

      /*** ACCEPTED PARTICLE WITH RECONSTRUCTED TRACK ***/
      
      /* get particle info */
      part = analysisParticle->GetPID();
      param[kPt] = analysisParticle->GetPt();
      param[kEta] = analysisParticle->GetEta();
      param[kPhi] = analysisParticle->GetPhi();
      charge = analysisParticle->GetSign() > 0. ? kPositive : kNegative;

      /* fill reconstructed histo */
      hHisto[kReconstructedTracks][part][charge]->Fill(param);
      
    }

  }

  TFile *fileout = TFile::Open(Form("TrackingEff.%s", filename), "RECREATE");
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++)
	hHisto[ihisto][ipart][icharge]->Write();
  fileout->Close();

  TString str = Form("TrackingEff.%s", filename);
  TrackingEff_efficiencyPt(str.Data());
  
}

//_____________________________________________________________________________-

TrackingEff_efficiencyPt(const Char_t *filename)
{

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoPt_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hHistoPt_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoAllPt_MB[kNHistos], *hHistoAllPt_centrality[NcentralityBins][kNHistos];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {

    /* INCLUSIVE */

    hHistoAllPt_MB[ihisto] = new TH1D(Form("hHistoAllPt_MB_%s", histoName[ihisto]), "", NptBins, ptBin);
    for (Int_t icent = 0; icent < NcentralityBins; icent++)
      hHistoAllPt_centrality[icent][ihisto] = new TH1D(Form("hHistoAllPt_centrality%d_%s", icent, histoName[ihisto]), "", NptBins, ptBin);

    /* SINGLE PARTICLE */

    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {

	/* get histo */
	hHisto[ihisto][ipart][icharge] = (THnSparseF *)filein->Get(Form("hHisto_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));

	/* MB projection */
	hHistoPt_MB[ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kPt);
	hHistoPt_MB[ihisto][ipart][icharge]->SetName(Form("hHistoPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hHistoPt_MB[ihisto][ipart][icharge]->Sumw2();
	hHistoAllPt_MB[ihisto]->Add(hHistoPt_MB[ihisto][ipart][icharge]);
	
	/* centrality projection */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	  hHisto[ihisto][ipart][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	  hHistoPt_centrality[icent][ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kPt);
	  hHistoPt_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hHistoPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hHistoPt_centrality[icent][ihisto][ipart][icharge]->Sumw2();
	  hHistoAllPt_centrality[icent][ihisto]->Add(hHistoPt_centrality[icent][ihisto][ipart][icharge]);
	}
      }
    }
  }    
  
  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyPt");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiencyPt_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyPt_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyPt_ratioMB_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hEfficiencyAllPt_MB[kNHistos], *hEfficiencyAllPt_centrality[NcentralityBins][kNHistos], *hEfficiencyAllPt_ratioMB_centrality[NcentralityBins][kNHistos];

  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {

    if (ihisto == kPrimaryTracks) continue;

    /* INCLUSIVE */

    /* MB efficiency */
    hEfficiencyAllPt_MB[ihisto] = new TH1D(*hHistoAllPt_MB[ihisto]);
    hEfficiencyAllPt_MB[ihisto]->SetName(Form("hEfficiencyAllPt_MB_%s", histoName[ihisto]));
    hEfficiencyAllPt_MB[ihisto]->SetLineWidth(2);
    hEfficiencyAllPt_MB[ihisto]->SetLineColor(1);
    hEfficiencyAllPt_MB[ihisto]->SetMarkerStyle(20);
    hEfficiencyAllPt_MB[ihisto]->SetMarkerColor(1);
    hEfficiencyAllPt_MB[ihisto]->Divide(hEfficiencyAllPt_MB[ihisto], hHistoAllPt_MB[kPrimaryTracks], 1, 1, "B");
    hEfficiencyAllPt_MB[ihisto]->Write();
    
    /* multiplicity/centrality efficiency */
    for (Int_t icent = 0; icent < NcentralityBins; icent++) {
      hEfficiencyAllPt_centrality[icent][ihisto] = new TH1D(*hHistoAllPt_centrality[icent][ihisto]);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetName(Form("hEfficiencyAllPt_centrality%d_%s", icent, histoName[ihisto]));
      hEfficiencyAllPt_centrality[icent][ihisto]->SetLineWidth(2);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetLineColor(multcentColor[icent]);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetMarkerStyle(20);
      hEfficiencyAllPt_centrality[icent][ihisto]->SetMarkerColor(multcentColor[icent]);
      hEfficiencyAllPt_centrality[icent][ihisto]->Divide(hEfficiencyAllPt_centrality[icent][ihisto], hHistoAllPt_centrality[icent][kPrimaryTracks], 1, 1, "B");
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

    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {

 	/* MB efficiency */
	hEfficiencyPt_MB[ihisto][ipart][icharge] = new TH1D(*hHistoPt_MB[ihisto][ipart][icharge]);
	hEfficiencyPt_MB[ihisto][ipart][icharge]->SetName(Form("hEfficiencyPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hEfficiencyPt_MB[ihisto][ipart][icharge]->SetLineWidth(2);
	hEfficiencyPt_MB[ihisto][ipart][icharge]->SetLineColor(1);
	hEfficiencyPt_MB[ihisto][ipart][icharge]->SetMarkerStyle(20);
	hEfficiencyPt_MB[ihisto][ipart][icharge]->SetMarkerColor(1);
	hEfficiencyPt_MB[ihisto][ipart][icharge]->Divide(hEfficiencyPt_MB[ihisto][ipart][icharge], hHistoPt_MB[kPrimaryTracks][ipart][icharge], 1, 1, "B");
	hEfficiencyPt_MB[ihisto][ipart][icharge]->Write();
	
	/* multiplicity/centrality efficiency */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hHistoPt_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyPt_centrality[icent][ihisto][ipart][icharge], hHistoPt_centrality[icent][kPrimaryTracks][ipart][icharge], 1, 1, "B");
	  hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]->Write();
	  
	  /* ratio wrt. MB */
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hEfficiencyPt_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyPt_ratioMB_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyPt_MB[ihisto][ipart][icharge]);
	  hEfficiencyPt_ratioMB_centrality[icent][ihisto][ipart][icharge]->Write();
	}
	
      }       
    }
  }
  
  fileout->Close();
	     
  TrackingEff_trackingEfficiency(str.Data());
}

//_____________________________________________________________________________-

TrackingEff_efficiencyEta(const Char_t *filename)
{

  /* get data */
  TFile *filein = TFile::Open(filename);
  THnSparseF *hHisto[kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoEta_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hHistoEta_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hHistoEta_MB_pt[NptsubBins][kNHistos][AliPID::kSPECIES][kNCharges], *hHistoEta_centrality_pt[NcentralityBins][NptsubBins][kNHistos][AliPID::kSPECIES][kNCharges];
  /* loop over histos */
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++)
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	
	/* get histo */
	hHisto[ihisto][ipart][icharge] = (THnSparseF *)filein->Get(Form("hHisto_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	
	/* MB projection */
	hHisto[ihisto][ipart][icharge]->GetAxis(kPt)->SetRange(0, 0);
	hHisto[ihisto][ipart][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
	hHistoEta_MB[ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kEta);
	hHistoEta_MB[ihisto][ipart][icharge]->SetName(Form("hHistoEta_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hHistoEta_MB[ihisto][ipart][icharge]->Sumw2();
	/* pt bins */
	for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	  hHisto[ihisto][ipart][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	  hHisto[ihisto][ipart][icharge]->GetAxis(kCentrality)->SetRange(0, 0);
	  hHistoEta_MB_pt[ipt][ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kEta);
	  hHistoEta_MB_pt[ipt][ihisto][ipart][icharge]->SetName(Form("hHistoEta_MB_pt%d_%s_%s_%s", ipt, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hHistoEta_MB_pt[ipt][ihisto][ipart][icharge]->Sumw2();
	}
	
	/* centrality projection */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	  hHisto[ihisto][ipart][icharge]->GetAxis(kPt)->SetRange(0, 0);
	  hHisto[ihisto][ipart][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	  hHistoEta_centrality[icent][ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kEta);
	  hHistoEta_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hHistoEta_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hHistoEta_centrality[icent][ihisto][ipart][icharge]->Sumw2();
	  /* pt bins */
	  for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	    hHisto[ihisto][ipart][icharge]->GetAxis(kPt)->SetRange(ptsubBinMin[ipt] + 1, ptsubBinMax[ipt] + 1);
	    hHisto[ihisto][ipart][icharge]->GetAxis(kCentrality)->SetRange(icent + 1, icent + 1);
	    hHistoEta_centrality_pt[icent][ipt][ihisto][ipart][icharge] = hHisto[ihisto][ipart][icharge]->Projection(kEta);
	    hHistoEta_centrality_pt[icent][ipt][ihisto][ipart][icharge]->SetName(Form("hHistoEta_centrality%d_pt%d_%s_%s_%s", icent, ipt, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	    hHistoEta_centrality_pt[icent][ipt][ihisto][ipart][icharge]->Sumw2();
	  }
	}
      }
  
  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".efficiencyEta");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");
  
  /* efficiencies/fractions and write */
  TH1D *hEfficiencyEta_MB[kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyEta_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyEta_ratioMB_centrality[NcentralityBins][kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hEfficiencyEta_MB_pt[NptsubBins][kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyEta_centrality_pt[NcentralityBins][NptsubBins][kNHistos][AliPID::kSPECIES][kNCharges], *hEfficiencyEta_ratioMB_centrality_pt[NcentralityBins][NptsubBins][kNHistos][AliPID::kSPECIES][kNCharges];
  
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	
	if (ihisto == kPrimaryTracks) continue;
	
	/* MB efficiency */
	hEfficiencyEta_MB[ihisto][ipart][icharge] = new TH1D(*hHistoEta_MB[ihisto][ipart][icharge]);
	hEfficiencyEta_MB[ihisto][ipart][icharge]->SetName(Form("hEfficiencyEta_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	hEfficiencyEta_MB[ihisto][ipart][icharge]->SetLineWidth(2);
	hEfficiencyEta_MB[ihisto][ipart][icharge]->SetLineColor(1);
	hEfficiencyEta_MB[ihisto][ipart][icharge]->SetMarkerStyle(20);
	hEfficiencyEta_MB[ihisto][ipart][icharge]->SetMarkerColor(1);
	hEfficiencyEta_MB[ihisto][ipart][icharge]->Divide(hEfficiencyEta_MB[ihisto][ipart][icharge], hHistoEta_MB[kPrimaryTracks][ipart][icharge], 1, 1, "B");
	hEfficiencyEta_MB[ihisto][ipart][icharge]->Write();
	/* pt bins */
	for (Int_t ipt = 0; ipt < NptsubBins; ipt++) {
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge] = new TH1D(*hHistoEta_MB_pt[ipt][ihisto][ipart][icharge]);
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->SetName(Form("hEfficiencyEta_MB_pt%d_%s_%s_%s", ipt, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->SetLineColor(1);
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->SetMarkerColor(1);
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->Divide(hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge], hHistoEta_MB_pt[ipt][kPrimaryTracks][ipart][icharge], 1, 1, "B");
	  hEfficiencyEta_MB_pt[ipt][ihisto][ipart][icharge]->Write();
	}
	
	/* multiplicity/centrality efficiency */
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hHistoEta_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyEta_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyEta_centrality[icent][ihisto][ipart][icharge], hHistoEta_centrality[icent][kPrimaryTracks][ipart][icharge], 1, 1, "B");
	  hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]->Write();
	  
	  /* ratio wrt. MB */
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge] = new TH1D(*hEfficiencyEta_centrality[icent][ihisto][ipart][icharge]);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetName(Form("hEfficiencyEta_ratioMB_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineWidth(2);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetLineColor(multcentColor[icent]);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerStyle(20);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->SetMarkerColor(multcentColor[icent]);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->Divide(hEfficiencyEta_MB[ihisto][ipart][icharge]);
	  hEfficiencyEta_ratioMB_centrality[icent][ihisto][ipart][icharge]->Write();
	}
	
      }       
  }
  
  fileout->Close();
  
}

//_____________________________________________________________________________-

TrackingEff_centralityDependence(const Char_t *filename)
{

  Double_t fitMin[AliPID::kSPECIES] = {0.5, 0.5, 0.5, 0.5, 0.5};
  Double_t fitMax[AliPID::kSPECIES] = {3.0, 3.0, 3.0, 3.0, 5.0};

  TF1 *pol0 = (TF1 *)gROOT->GetFunction("pol0");
  TF1 *pol1 = (TF1 *)gROOT->GetFunction("pol1");
  pol0->SetRange(0., 5.0);
  TFile *filein = TFile::Open(filename);

  /* output */
  TString str = filename;
  str.Insert(str.Length() - TString(".root").Length(), ".centralityDependence");
  TFile *fileout = TFile::Open(str.Data(), "RECREATE");

  TH1D *hEfficiencyPt_ratioMB;
  TH1D *hEfficiencyCentrality[kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hChi2Centrality[kNHistos][AliPID::kSPECIES][kNCharges];
  TH1D *hEfficiencyAllCentrality[kNHistos];
  TH1F *hHistoEffRatio = new TH1F("hHistoEffRatio", "", 100, 0.5, 1.5);
  for (Int_t ihisto = 0; ihisto < kNHistos; ihisto++) {
    if (ihisto == kPrimaryTracks) continue;

    /* INCLUSIVE */

    hEfficiencyAllCentrality[ihisto] = new TH1D(Form("hEfficiencyAllCentrality_ratioMB_%s", histoName[ihisto]), "", NcentralityBins, centralityBin);
    hEfficiencyAllCentrality[ihisto]->SetLineWidth(2);
    hEfficiencyAllCentrality[ihisto]->SetLineColor(1);
    hEfficiencyAllCentrality[ihisto]->SetMarkerStyle(20);
    hEfficiencyAllCentrality[ihisto]->SetMarkerColor(1);
    
    for (Int_t icent = 0; icent < NcentralityBins; icent++) {
      
      hEfficiencyPt_ratioMB = (TH1D *)filein->Get(Form("hEfficiencyAllPt_ratioMB_centrality%d_%s", icent, histoName[ihisto]));
      hEfficiencyPt_ratioMB->Fit(pol0, "q0", "", 0., 5.0);
      hEfficiencyAllCentrality[ihisto]->SetBinContent(icent + 1, pol0->GetParameter(0));
      hEfficiencyAllCentrality[ihisto]->SetBinError(icent + 1, pol0->GetParError(0));
      hEfficiencyPt_ratioMB->Add(pol0, -1.);
      hEfficiencyPt_ratioMB->Write();
      
    }
    hEfficiencyAllCentrality[ihisto]->Write();
    
    
    /* SINGLE PARTICLE */
    
    for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++)
      for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
	
	hEfficiencyCentrality[ihisto][ipart][icharge] = new TH1D(Form("hEfficiencyCentrality_ratioMB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), "", NcentralityBins, centralityBin);
	hEfficiencyCentrality[ihisto][ipart][icharge]->SetLineWidth(2);
	hEfficiencyCentrality[ihisto][ipart][icharge]->SetLineColor(particleColor[ipart]);
	hEfficiencyCentrality[ihisto][ipart][icharge]->SetMarkerStyle(chargeMarker[icharge]);
	hEfficiencyCentrality[ihisto][ipart][icharge]->SetMarkerColor(particleColor[ipart]);
	
	for (Int_t icent = 0; icent < NcentralityBins; icent++) {
	 
	  hEfficiencyPt_ratioMB = (TH1D *)filein->Get(Form("hEfficiencyPt_ratioMB_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));

	  hEfficiencyPt_ratioMB->Fit(pol0, "q0", "", fitMin[ipart], fitMax[ipart]);
	  hEfficiencyCentrality[ihisto][ipart][icharge]->SetBinContent(icent + 1, pol0->GetParameter(0));
	  hEfficiencyCentrality[ihisto][ipart][icharge]->SetBinError(icent + 1, pol0->GetParError(0));
	  pol0->SetRange(fitMin[ipart], fitMax[ipart]);
	  hEfficiencyPt_ratioMB->Add(pol0, -1.);
	  hEfficiencyPt_ratioMB->Write();

	}

	hEfficiencyCentrality[ihisto][ipart][icharge]->Write();
      }
  }
  
  fileout->Close();
}

//_____________________________________________________________________________-

TrackingEff_centralityDependenceFit(const Char_t *filename, Int_t ihisto = kReconstructedTracks)
{

  const Char_t *particleLabel[5][2] = {"", "", "", "", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};

  TF1 *fCentralityDependence = new TF1("fCentralityDependence", "[0] + [1] * (1. - TMath::Exp(-x / [2]))", 0., 90.);
  fCentralityDependence->SetParameter(0, 0.98);
  fCentralityDependence->SetParameter(1, 0.05);
  fCentralityDependence->SetParameter(2, 50.);

  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyAllCentrality;
  TH1D *hEfficiencyCentrality[AliPID::kSPECIES][kNCharges];
  hEfficiencyAllCentrality = (TH1D *)filein->Get(Form("hEfficiencyAllCentrality_ratioMB_%s", histoName[ihisto]));
  TCanvas *cFit = new TCanvas("cFit");
  hEfficiencyAllCentrality->Fit(fCentralityDependence);
  hEfficiencyAllCentrality->SetMaximum(1.05);
  hEfficiencyAllCentrality->SetMinimum(0.95);
  TCanvas *cRatios = new TCanvas("cRatios");
  cRatios->Divide(2, 3);
  for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < 2; icharge++) {
      cRatios->cd(icharge + 1 + 2 * (ipart - 2));
      cRatios->cd(icharge + 1 + 2 * (ipart - 2))->SetGridx();
      cRatios->cd(icharge + 1 + 2 * (ipart - 2))->SetGridy();
      hEfficiencyCentrality[ipart][icharge] = (TH1D *)filein->Get(Form("hEfficiencyCentrality_ratioMB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
      hEfficiencyCentrality[ipart][icharge]->Divide(fCentralityDependence);
      hEfficiencyCentrality[ipart][icharge]->SetMaximum(1.05);
      hEfficiencyCentrality[ipart][icharge]->SetMinimum(0.95);
      hEfficiencyCentrality[ipart][icharge]->SetTitle(Form("%s;centrality percentile;efficiency ratio wrt. fitted centrality dependence", particleLabel[ipart][icharge]));
      hEfficiencyCentrality[ipart][icharge]->SetStats(kFALSE);
      hEfficiencyCentrality[ipart][icharge]->Draw();
    }
  }
  

}


//_____________________________________________________________________________-

TrackingEff_efficiencyPt_MB_plot(const Char_t *filename)
{
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 2, 0, 20, 4);
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 2, 1, 25, 4);
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 3, 0, 20, 8);
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 3, 1, 25, 8);
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 4, 0, 20, 2);
  TrackingEff_efficiencyPt_MB_plot(filename, kReconstructedTracks, 4, 1, 25, 2);
}

TH1D *
TrackingEff_efficiencyPt_MB_plot(const Char_t *filename, Int_t ihisto = kReconstructedTracks, Int_t ipart, Int_t icharge, Int_t marker = 20, Int_t color = 2, Option_t *opt = "")
{
  
  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  
  Double_t ptMin[AliPID::kSPECIES] = {0.5, 0.5, 0.5, 0.5, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {3.0, 3.0, 3.0, 3.0, 5.0};

  TF1 *fEff = new TF1("fEff", "[0] + [1] * x - [2] * TMath::Exp(-[3] * TMath::Power(x, [4]))", 0.5, 5.0);
  fEff->SetParameter(0, 0.8);
  fEff->SetParameter(1, -0.01);
  fEff->SetParameter(2, 2.0);
  fEff->SetParameter(3, 1.);
  fEff->SetParameter(4, 2.);

  TFile *fileout = TFile::Open(Form("%s/efficiencyPt_MB_%s_%s.root", destdir, AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");

  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
  hEfficiencyPt->Fit(fEff, "0", "IME", 0.5, 5.0);
  hEfficiencyPt->SetTitle(Form("%s;p_{T} (GeV/c);acceptance #times efficiency;", partChargeName[ipart][icharge]));
  hEfficiencyPt->SetMinimum(0.4);
  hEfficiencyPt->SetMaximum(0.8);
  hEfficiencyPt->SetMarkerStyle(marker);
  hEfficiencyPt->SetMarkerColor(color);
  hEfficiencyPt->SetMarkerSize(1.5);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->SetLineColor(color);
  hEfficiencyPt->GetXaxis()->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
  hEfficiencyPt->SetStats(kFALSE);
  cCanvas1->cd();
  hEfficiencyPt->DrawCopy();
  fEff->DrawCopy("same");
  fileout->cd();
  hEfficiencyPt->Write("hEfficiencyPt");
  fEff->Write("fEff");
  
  cCanvas1->SetGridx();
  cCanvas1->SetGridy();
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s_%s.C", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s_%s.png", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_%s_%s.eps", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  
  TH1D *hRatioPt = new TH1D(*hEfficiencyPt);
  hRatioPt->Divide(fEff);
  hRatioPt->SetTitle(Form("%s;p_{T} (GeV/c);ratio wrt. fitted dependence;", partChargeName[ipart][icharge]));
  hRatioPt->SetMinimum(0.9);
  hRatioPt->SetMaximum(1.1);
  cCanvas2->cd();
  hRatioPt->DrawCopy();
  fileout->cd();
  hRatioPt->Write("hRatioPt");

  cCanvas2->SetGridx();
  cCanvas2->SetGridy();
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s_%s.C", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s_%s.png", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_%s_%s.eps", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  

  fileout->Close();

  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyPt;
}

//_____________________________________________________________________________-

TH1D *
TrackingEff_efficiencyPt_centrality_all_plot(const Char_t *filename, Int_t ihisto = kReconstructedTracks, Int_t marker = 20, Int_t color = 1, Option_t *opt = "")
{
  
  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  
  Double_t ptMin[AliPID::kSPECIES] = {0.5, 0.5, 0.5, 0.5, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {3.0, 3.0, 3.0, 3.0, 5.0};
  
  TF1 *fEff = new TF1("fEff", "[5] * ([0] + [1] * x - [2] * TMath::Exp(-[3] * TMath::Power(x, [4])))", 0.5, 5.0);
  fEff->SetParameter(0, 0.8);
  fEff->SetParameter(1, -0.01);
  fEff->SetParameter(2, 2.0);
  fEff->SetParameter(3, 1.);
  fEff->SetParameter(4, 2.);
  fEff->FixParameter(5, 1.);

  TFile *fileout = TFile::Open(Form("%s/efficiencyPt_MB_all.root", destdir), "RECREATE");
  
  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyAllPt_MB_%s", histoName[ihisto]));
  hEfficiencyPt->Fit(fEff, "0", "IME", 0.5, 5.0);
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
  hEfficiencyPt->DrawCopy();
  fEff->DrawCopy("same");
  fileout->cd();
  hEfficiencyPt->Write("hEfficiencyPt");
  fEff->Write("fEff");

  cCanvas1->SetGridx();
  cCanvas1->SetGridy();
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.C", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.png", destdir));
  cCanvas1->SaveAs(Form("%s/efficiencyPt_MB_all.eps", destdir));

  TH1D *hRatioPt = new TH1D(*hEfficiencyPt);
  hRatioPt->Divide(fEff);
  hRatioPt->SetTitle("all particles;p_{T} (GeV/c);ratio wrt. fitted dependence;");
  hRatioPt->SetMinimum(0.9);
  hRatioPt->SetMaximum(1.1);
  cCanvas2->cd();
  hRatioPt->DrawCopy();
  fileout->cd();
  hRatioPt->Write("hRatioPt");

  cCanvas2->SetGridx();
  cCanvas2->SetGridy();
  cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_MB_all.C", destdir));
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
    hEfficiencyPt_cent[icent]->Fit(fEff, "0", "IME", 1.0, 3.0);
    hEfficiencyPt_cent[icent]->Fit(fEff, "", "IME", 0.5, 5.0);
    
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
    hEfficiencyPt_cent[icent]->DrawCopy();
    fEff->DrawCopy("same");
    fileout->cd();
    hEfficiencyPt_cent[icent]->Write(Form("hEfficiencyPt_cent%d", icent));
    fEff->Write(Form("fEff_cent%d", icent));
    
    hEfficiencyCent->SetBinContent(icent + 1, fEff->GetParameter(5));
    hEfficiencyCent->SetBinError(icent + 1, fEff->GetParError(5));
    
    cCanvas1->SetGridx();
    cCanvas1->SetGridy();
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.C", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.png", destdir, icent));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_all.eps", destdir, icent));
    
    hRatioPt_cent[icent] = new TH1D(*hEfficiencyPt_cent[icent]);
    hRatioPt_cent[icent]->Divide(fEff);
    hRatioPt_cent[icent]->SetTitle(Form("all particles (%d-%d\%);p_{T} (GeV/c);ratio wrt. fitted dependence;", (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hRatioPt_cent[icent]->SetMinimum(0.9);
    hRatioPt_cent[icent]->SetMaximum(1.1);
    cCanvas2->cd();
    hRatioPt_cent[icent]->DrawCopy();
    fileout->cd();
    hRatioPt_cent[icent]->Write(Form("hRatio_cent%d", icent));

    cCanvas2->SetGridx();
    cCanvas2->SetGridy();
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.C", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.png", destdir, icent));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_all.eps", destdir, icent));

    
  }

  TF1 *fEffCent = new TF1("fEffCent", "[0] - [1] * TMath::Exp(-[2] * TMath::Power(x, [3]))", 0., 90.);
  fEffCent->SetParameter(0, 1.02);
  fEffCent->SetParameter(1, 0.04);
  fEffCent->SetParameter(2, 0.001);
  fEffCent->SetParameter(3, 2.);
  hEfficiencyCent->Fit(fEffCent, "q0", "IME", 0., 90.);
  
  TCanvas *cCanvas3 = new TCanvas("cCanvas3");
  hEfficiencyCent->DrawCopy();
  fEffCent->DrawCopy("same");
  fileout->cd();
  hEfficiencyCent->Write("hEfficiencyCent");
  fEffCent->Write("fEffCent");
  
  cCanvas3->SetGridx();
  cCanvas3->SetGridy();
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.C", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.png", destdir));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_all.eps", destdir));

  TCanvas *cCanvas4 = new TCanvas("cCanvas4");

  TH1D *hRatioCent = new TH1D(*hEfficiencyCent);
  hRatioCent->Divide(fEffCent);
  hRatioCent->SetTitle(Form("all particles;centrality percentile;ratio wrt. fitted dependence;"));
  hRatioCent->SetMinimum(0.95);
  hRatioCent->SetMaximum(1.05);
  cCanvas4->cd();
  hRatioCent->DrawCopy();
  fileout->cd();
  hRatioCent->Write("hRatioCent");

  cCanvas4->SetGridx();
  cCanvas4->SetGridy();
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.C", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.png", destdir));
  cCanvas4->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_ratioFit_all.eps", destdir));
  
  fileout->Close();

  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyCent;
}


//_____________________________________________________________________________-

TrackingEff_efficiencyPt_centrality_plot(const Char_t *filename)
{
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 2, 0, 20, 4);
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 2, 1, 25, 4);
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 3, 0, 20, 8);
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 3, 1, 25, 8);
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 4, 0, 20, 2);
  TrackingEff_efficiencyPt_centrality_plot(filename, kReconstructedTracks, 4, 1, 25, 2);
}

TH1D *
TrackingEff_efficiencyPt_centrality_plot(const Char_t *filename, Int_t ihisto = kReconstructedTracks, Int_t ipart, Int_t icharge, Int_t marker = 20, Int_t color = 2, Option_t *opt = "")
{

  TVirtualFitter::SetMaxIterations(1000000);

  /* load HistoUtils */
  gROOT->LoadMacro("HistoUtils.C");
  
  TCanvas *cCanvas1 = new TCanvas("cCanvas1");
  TCanvas *cCanvas2 = new TCanvas("cCanvas2");
  TCanvas *cCanvas5 = new TCanvas("cCanvas5");

  Double_t ptMin[AliPID::kSPECIES] = {0.5, 0.5, 0.5, 0.5, 0.5};
  Double_t ptMax[AliPID::kSPECIES] = {3.0, 3.0, 3.0, 3.0, 5.0};

  /* fit minimum-bias efficiency pt */

  TF1 *fEff = new TF1("fEff", "([0] * TMath::Exp(-TMath::Power([1] / x, [2])) + [3] * x) * TMath::Exp(TMath::Power([4] / x, [5]))", ptMin[ipart], ptMax[ipart]);
  fEff->SetParameter(0, 0.5);
  fEff->SetParameter(1, 0.1);
  fEff->SetParameter(2, 2.);
  fEff->SetParameter(3, 0.);
  fEff->SetParameter(4, 0.001);
  fEff->SetParameter(5, 1.);
  Int_t nPars = fEff->GetNpar();
  Int_t scalePar = 0.;
  TF1 *fEffCopy = new TF1(*fEff);

  TFile *fileout = TFile::Open(Form("%s/efficiencyPt_centrality_%s_%s.root", destdir, AliPID::ParticleName(ipart), chargeName[icharge]), "RECREATE");

  TFile *filein = TFile::Open(filename);
  TH1D *hEfficiencyPt = (TH1D *)filein->Get(Form("hEfficiencyPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
  hEfficiencyPt->Fit(fEff, "0", "IMRE", ptMin[ipart], ptMax[ipart]);

  /* build efficiency profile */
  for (Int_t ipar = 0; ipar < nPars; ipar++) {
    fEffCopy->SetParameter(ipar, fEff->GetParameter(ipar));
    fEffCopy->SetParError(ipar, fEff->GetParError(ipar));
  }
  TProfile *pEfficiencyPt = new TProfile(Form("pEfficiencyPt_MB_%s_%s_%s", histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), Form("%s;p_{T} (GeV/c);acceptance #times efficiency;", partChargeName[ipart][icharge]), NptBins, ptBin, "s");
  HistoUtils_Function2Profile(fEffCopy, pEfficiencyPt);

  /* draw */
  cCanvas1->cd();
  pEfficiencyPt->SetMarkerStyle(marker);
  pEfficiencyPt->SetMarkerColor(color);
  pEfficiencyPt->SetMarkerSize(0);
  pEfficiencyPt->SetLineWidth(2);
  pEfficiencyPt->SetLineColor(color);
  pEfficiencyPt->SetFillStyle(3001);
  pEfficiencyPt->SetFillColor(color);
  pEfficiencyPt->GetXaxis()->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
  pEfficiencyPt->SetStats(kFALSE);
  pEfficiencyPt->DrawCopy("E2");

  hEfficiencyPt->SetTitle(Form("%s;p_{T} (GeV/c);acceptance #times efficiency;", partChargeName[ipart][icharge]));
  hEfficiencyPt->SetMarkerStyle(marker);
  hEfficiencyPt->SetMarkerColor(color);
  hEfficiencyPt->SetMarkerSize(1.5);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->SetLineColor(color);
  hEfficiencyPt->GetXaxis()->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
  hEfficiencyPt->SetStats(kFALSE);
  hEfficiencyPt->DrawCopy("same");
  fEff->DrawCopy("same");
  cCanvas1->Update();

  /* write */
  fileout->cd();
  pEfficiencyPt->Write("pEfficiencyPt");
  hEfficiencyPt->Write("hEfficiencyPt");
  fEff->Write("fEff");

  TH1D *hEfficiencyCent = new TH1D("hEfficiencyCent", Form("%s;centrality percentile;acceptance x efficiency scale factor;", partChargeName[ipart][icharge]), NcentralityBins, centralityBin);
  hEfficiencyCent->SetMinimum(0.9);
  hEfficiencyCent->SetMaximum(1.1);
  hEfficiencyCent->SetMarkerStyle(marker);
  hEfficiencyCent->SetMarkerColor(color);
  hEfficiencyCent->SetMarkerSize(1.5);
  hEfficiencyCent->SetLineWidth(2);
  hEfficiencyCent->SetLineColor(color);
  hEfficiencyCent->SetStats(kFALSE);
  
  TProfile *pEfficiencyPt_cent[NcentralityBins];
  TH1D *hEfficiencyPt_cent[NcentralityBins];
  TH1D *hRatioPt_cent[NcentralityBins];

  /* fix efficiency shape and release scale factor */
  for (Int_t ipar = 0; ipar < nPars; ipar++)
    fEff->FixParameter(ipar, fEff->GetParameter(ipar));
  fEff->ReleaseParameter(scalePar);
    
  gStyle->SetOptStat(1100);
  for (Int_t icent = 0; icent < NcentralityBins; icent++) {

    hEfficiencyPt_cent[icent] = (TH1D *)filein->Get(Form("hEfficiencyPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]));
    hEfficiencyPt_cent[icent]->Fit(fEff, "", "IME", ptMin[ipart], ptMax[ipart]);
    
    /* build efficiency profile */
    fEffCopy->SetParameter(scalePar, fEff->GetParameter(scalePar));
    fEffCopy->SetParError(scalePar, fEff->GetParError(scalePar));
    pEfficiencyPt_cent[icent] = new TProfile(Form("pEfficiencyPt_centrality%d_%s_%s_%s", icent, histoName[ihisto], AliPID::ParticleName(ipart), chargeName[icharge]), Form("%s;p_{T} (GeV/c);acceptance #times efficiency;", partChargeName[ipart][icharge]), NptBins, ptBin, "s");
    HistoUtils_Function2Profile(fEffCopy, pEfficiencyPt_cent[icent]);
    
    /* draw */
    cCanvas1->cd();
    pEfficiencyPt_cent[icent]->SetMarkerStyle(marker);
    pEfficiencyPt_cent[icent]->SetMarkerColor(color);
    pEfficiencyPt_cent[icent]->SetMarkerSize(0);
    pEfficiencyPt_cent[icent]->SetLineWidth(2);
    pEfficiencyPt_cent[icent]->SetLineColor(color);
    pEfficiencyPt_cent[icent]->SetFillStyle(3001);
    pEfficiencyPt_cent[icent]->SetFillColor(color);
    pEfficiencyPt_cent[icent]->GetXaxis()->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
    pEfficiencyPt_cent[icent]->SetStats(kFALSE);
    pEfficiencyPt_cent[icent]->DrawCopy("E2");
    
    hEfficiencyPt_cent[icent]->SetTitle(Form("%s (%d-%d\%);p_{T} (GeV/c);acceptance #times efficiency;", partChargeName[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hEfficiencyPt_cent[icent]->SetMarkerStyle(marker);
    hEfficiencyPt_cent[icent]->SetMarkerColor(color);
    hEfficiencyPt_cent[icent]->SetMarkerSize(1.5);
    hEfficiencyPt_cent[icent]->SetLineWidth(2);
    hEfficiencyPt_cent[icent]->SetLineColor(color);
    hEfficiencyPt_cent[icent]->GetXaxis()->SetRangeUser(ptMin[ipart] + 0.001, ptMax[ipart] - 0.001);
    hEfficiencyPt_cent[icent]->SetStats(kFALSE);
    hEfficiencyPt_cent[icent]->DrawCopy("same");
    fEff->DrawCopy("same");
    cCanvas1->Update();

    fileout->cd();
    pEfficiencyPt_cent[icent]->Write(Form("pEfficiencyPt_cent%d", icent));
    hEfficiencyPt_cent[icent]->Write(Form("hEfficiencyPt_cent%d", icent));
    fEff->Write(Form("fEff_cent%d", icent));

    hEfficiencyCent->SetBinContent(icent + 1, fEff->GetParameter(5));
    hEfficiencyCent->SetBinError(icent + 1, fEff->GetParError(5));

    cCanvas1->SetGridx();
    cCanvas1->SetGridy();
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_%s_%s.C", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_%s_%s.png", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    cCanvas1->SaveAs(Form("%s/efficiencyPt_centrality%d_%s_%s.eps", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));

    hRatioPt_cent[icent] = new TH1D(*hEfficiencyPt_cent[icent]);
    hRatioPt_cent[icent]->Divide(fEff);
    hRatioPt_cent[icent]->SetTitle(Form("%s (%d-%d\%);p_{T} (GeV/c);ratio wrt. fitted dependence;", partChargeName[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]));
    hRatioPt_cent[icent]->SetMinimum(0.75);
    hRatioPt_cent[icent]->SetMaximum(1.25);
    cCanvas2->cd();
    hRatioPt_cent[icent]->DrawCopy();
    fileout->cd();
    hRatioPt_cent[icent]->Write(Form("hRatioPt_cent%d", icent));
    

    cCanvas2->SetGridx();
    cCanvas2->SetGridy();
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_%s_%s.C", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_%s_%s.png", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    cCanvas2->SaveAs(Form("%s/efficiencyPt_ratioFit_centrality%d_%s_%s.eps", destdir, icent, AliPID::ParticleName(ipart), chargeName[icharge]));
    
  }

  TCanvas *cCanvas3 = new TCanvas("cCanvas3");
  hEfficiencyCent->DrawCopy();
  fileout->cd();
  hEfficiencyCent->Write("hEfficiencyCent");
  
  cCanvas3->SetGridx();
  cCanvas3->SetGridy();
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_%s_%s.C", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_%s_%s.png", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));
  cCanvas3->SaveAs(Form("%s/efficiencyCentrality_scaleFactor_%s_%s.eps", destdir, AliPID::ParticleName(ipart), chargeName[icharge]));

  //  hEfficiencyPt->Add(fEff, -1.);
  return hEfficiencyCent;
}

TrackingEff_trackingEfficiency(const Char_t *filename, Bool_t useGFcorrection =kTRUE)
{

  Int_t marker[2] = {20, 21};
  Int_t color[AliPID::kSPECIES] = {1, 1, 4, 8, 2};
  Char_t *partLatex[AliPID::kSPECIES][2] = {
    "", "", "#pi^{+}", "K^{+}", "p",
    "", "", "#pi^{-}", "K^{-}", "#bar{p}"
  };

  TFile *filein = TFile::Open(filename);
  if (useGFcorrection)
    TFile *fileout = TFile::Open("TOF_trackingEfficiency.root", "RECREATE");
  else
    TFile *fileout = TFile::Open("TOF_trackingEfficiency_noGF.root", "RECREATE");
  TH1D *hEff;
  TF1 *fGF;
  Char_t title[1024];
  for (Int_t icent = -1; icent < NcentralityBins; icent++)
    for (Int_t icharge = 0; icharge < kNCharges; icharge++)
      for (Int_t ipart = 2; ipart < AliPID::kSPECIES; ipart++) {
	if (icent == -1)
	  hEff = (TH1D *)filein->Get(Form("hEfficiencyPt_MB_hReconstructedTracks_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	else
	  hEff = (TH1D *)filein->Get(Form("hEfficiencyPt_centrality%d_hReconstructedTracks_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	/* geant-fluka correction */
	fGF = TrackingEff_geantflukaCorrection(ipart, icharge);
	if (useGFcorrection)
	  hEff->Divide(fGF);
	if (icent == -1)
	  sprintf(title, "%s (MB);p_{T} (GeV/c);tracking efficiency;", partLatex[ipart][icharge]);
	else
	  sprintf(title, "%s (%d-%d%%);p_{T} (GeV/c);tracking efficiency;", partLatex[ipart][icharge], (Int_t)centralityBin[icent], (Int_t)centralityBin[icent + 1]);
	hEff->SetMarkerStyle(marker[icharge]);
	hEff->SetMarkerColor(color[ipart]);
	hEff->SetLineColor(1);
	hEff->SetLineWidth(1);
	hEff->SetTitle(title);
	if (icent == -1)
	  hEff->SetName(Form("hTrackingEff_MB_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]));
	else
	  hEff->SetName(Form("hTrackingEff_cent%d_%s_%s", icent, AliPID::ParticleName(ipart), chargeName[icharge]));
	fileout->cd();
	hEff->Write();
      }

  fileout->Close();
}

TF1 *
TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  //  else if (ipart == 4 && icharge == kNegative) {
  //    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  //  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t
TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t
TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t
TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}
