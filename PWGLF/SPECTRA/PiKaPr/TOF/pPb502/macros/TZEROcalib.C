#include "CommonDefs.C"

Int_t ncentbins = 5;
Int_t nvertexbins = 20;

TZEROcalib(const Char_t *filename = "data.root", const Char_t *calibfilename = NULL, Int_t evMax = kMaxInt)
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

  /* open calibfile */
  TH1 *hCentrality_TZEROA_mean = NULL;
  TH1 *hCentrality_TZEROA_sigma = NULL;
  TH1 *hCentrality_TZEROC_mean = NULL;
  TH1 *hCentrality_TZEROC_sigma = NULL;
  TH1 *hCentrality_TOF_mean = NULL;
  TH1 *hCentrality_TOF_TZEROA_mean = NULL;
  TH1 *hCentrality_TOF_TZEROC_mean = NULL;
  TH1 *hCentrality_TOF_TZEROTOF_mean = NULL;
  if (calibfilename) {
    TFile *calibfile = TFile::Open(calibfilename);
    hCentrality_TZEROA_mean = (TH1 *)calibfile->Get("hCentrality_TZEROA_mean");
    hCentrality_TZEROA_sigma = (TH1 *)calibfile->Get("hCentrality_TZEROA_sigma");
    hCentrality_TZEROC_mean = (TH1 *)calibfile->Get("hCentrality_TZEROC_calib");
    hCentrality_TZEROC_sigma = (TH1 *)calibfile->Get("hCentrality_TZEROC_sigma");
    hCentrality_TOF_mean = (TH1 *)calibfile->Get("hCentrality_TOF_mean");
    hCentrality_TOF_TZEROA_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROA_mean");
    hCentrality_TOF_TZEROC_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROC_mean");
    hCentrality_TOF_TZEROTOF_mean = (TH1 *)calibfile->Get("hCentrality_TOF_TZEROTOF_mean");
  }

  /* histos */
  TH1F *hCentrality = new TH1F("hCentrality", "", NcentralityBins, centralityBin);
  TH1F *hCentrality_AC = new TH1F("hCentrality_AC", "", NcentralityBins, centralityBin);
  TH1F *hCentrality_A = new TH1F("hCentrality_A", "", NcentralityBins, centralityBin);
  TH1F *hCentrality_C = new TH1F("hCentrality_C", "", NcentralityBins, centralityBin);
  TH1F *hCentrality_NONE = new TH1F("hCentrality_NONE", "", NcentralityBins, centralityBin);

  TH2F *hCentrality_TZEROA = new TH2F("hCentrality_TZEROA", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROC = new TH2F("hCentrality_TZEROC", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROOR = new TH2F("hCentrality_TZEROOR", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROAND = new TH2F("hCentrality_TZEROAND", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROTOF = new TH2F("hCentrality_TZEROTOF", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZERODIFF = new TH2F("hCentrality_TZERODIFF", "", NcentralityBins, centralityBin, 200, -2440., 2440.);

  TH2F *hCentrality_TZEROTOF_TZEROA = new TH2F("hCentrality_TZEROTOF_TZEROA", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROTOF_TZEROC = new TH2F("hCentrality_TZEROTOF_TZEROC", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROTOF_TZEROAND = new TH2F("hCentrality_TZEROTOF_TZEROAND", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TZEROA_TZEROC = new TH2F("hCentrality_TZEROA_TZEROC", "", NcentralityBins, centralityBin, 200, -2440., 2440.);

  TH2F *hVertex_TZEROA = new TH2F("hVertex_TZEROA", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TZEROC = new TH2F("hVertex_TZEROC", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TZEROOR = new TH2F("hVertex_TZEROOR", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TZEROAND = new TH2F("hVertex_TZEROAND", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TZEROTOF = new TH2F("hVertex_TZEROTOF", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TZERODIFF = new TH2F("hVertex_TZERODIFF", "", nvertexbins, -10., 10., 200, -2440., 2440.);


  TH2F *hCentrality_TOF = new TH2F("hCentrality_TOF", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TOF_TZEROA = new TH2F("hCentrality_TOF_TZEROA", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TOF_TZEROC = new TH2F("hCentrality_TOF_TZEROC", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TOF_TZEROAND = new TH2F("hCentrality_TOF_TZEROAND", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_TOF_TZEROTOF = new TH2F("hCentrality_TOF_TZEROTOF", "", NcentralityBins, centralityBin, 200, -2440., 2440.);
  TH2F *hCentrality_Resolution = new TH2F("hCentrality_Resolution", "", NcentralityBins, centralityBin, 20, 0., 100.);

  TH2F *hResolution_TZERODIFF = new TH2F("hResolution_TZERODIFF", "", 20, 0., 100., 200, -2440., 2440.);
  TH2F *hResolution_TOF_TZEROA = new TH2F("hResolution_TOF_TZEROA", "", 20, 0., 100., 200, -2440., 2440.);
  TH2F *hResolution_TOF_TZEROC = new TH2F("hResolution_TOF_TZEROC", "", 20, 0., 100., 200, -2440., 2440.);
  TH2F *hResolution_TOF_TZEROAND = new TH2F("hResolution_TOF_TZEROAND", "", 20, 0., 100., 200, -2440., 2440.);
  TH2F *hResolution_TOF_TZEROTOF = new TH2F("hResolution_TOF_TZEROTOF", "", 20, 0., 100., 200, -2440., 2440.);

  TH2F *hVertex_TOF = new TH2F("hVertex_TOF", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TOF_TZEROA = new TH2F("hVertex_TOF_TZEROA", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TOF_TZEROC = new TH2F("hVertex_TOF_TZEROC", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TOF_TZEROOR = new TH2F("hVertex_TOF_TZEROOR", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TOF_TZEROAND = new TH2F("hVertex_TOF_TZEROAND", "", nvertexbins, -10., 10., 200, -2440., 2440.);
  TH2F *hVertex_TOF_TZEROTOF = new TH2F("hVertex_TOF_TZEROTOF", "", nvertexbins, -10., 10., 200, -2440., 2440.);

  TH2F *hDiffDiff = new TH2F("hDiffDiff", "", 200, -2440., 2440., 200, -2440., 2440.);

  /* loop over events */
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 100000 == 0) printf("iev = %d\n", iev);
    /* check event */
    if (!analysisEvent->AcceptEvent(acceptEventType)) continue;
    
    /*** ACCEPTED EVENT ***/

    /* get vertex position */
    Double_t vertexz = analysisEvent->GetVertexZ();

    /* get centrality */
    Double_t cent = analysisEvent->GetCentralityPercentile(centralityEstimator);
    hCentrality->Fill(cent);

    Int_t icent;
    for (icent = 0; icent < NcentralityBins; icent++)
      if (cent < centralityBin[icent + 1])
	break;

    Double_t TZEROA_mean = hCentrality_TZEROA_mean ? hCentrality_TZEROA_mean->GetBinContent(icent + 1) : 0.;
    Double_t TZEROA_sigma = hCentrality_TZEROA_sigma ? hCentrality_TZEROA_sigma->GetBinContent(icent + 1) : 1000.;
    Double_t TZEROC_mean  = hCentrality_TZEROC_mean ? hCentrality_TZEROC_mean->GetBinContent(icent + 1) : 0.;
    Double_t TZEROC_sigma = hCentrality_TZEROC_sigma ? hCentrality_TZEROC_sigma->GetBinContent(icent + 1) : 1000.;

    Double_t TOF_mean = hCentrality_TOF_mean ? hCentrality_TOF_mean->GetBinContent(icent + 1) : 0.;
    Double_t TOF_TZEROA_mean = hCentrality_TOF_TZEROA_mean ? hCentrality_TOF_TZEROA_mean->GetBinContent(icent + 1) : 0.;
    Double_t TOF_TZEROC_mean = hCentrality_TOF_TZEROC_mean ? hCentrality_TOF_TZEROC_mean->GetBinContent(icent + 1) : 0.;
    Double_t TOF_TZEROTOF_mean = hCentrality_TOF_TZEROTOF_mean ? hCentrality_TOF_TZEROTOF_mean->GetBinContent(icent + 1) : 0.;
    
    
    /* TZERO */
    Double_t TZEROA = analysisEvent->GetTimeZeroT0(1) - TZEROA_shift;
    Bool_t hasTZEROA = TMath::Abs(TZEROA - TZEROA_mean) < 3. * TZEROA_sigma;
    Double_t TZEROC = analysisEvent->GetTimeZeroT0(2) - TZEROC_shift;
    Bool_t hasTZEROC = TMath::Abs(TZEROC - TZEROC_mean) < 3. * TZEROC_sigma;
    /* vertex correction */
    //    TZEROA += -TZEROvertexCorr * vertexz;
    //    TZEROC += TZEROvertexCorr * vertexz;
    /* alignment to TOF */
    TZEROA += TOF_TZEROA_mean - TOF_mean;
    TZEROC += TOF_TZEROC_mean - TOF_mean;
    /* TZEROAND */
    Double_t TZEROAND = (TZEROA + TZEROC) * 0.5;
    Bool_t hasTZEROAND = hasTZEROA && hasTZEROC;
    /* TZEROTOF */
    Double_t TZEROTOF = analysisEvent->GetTimeZeroTOF()[9];
    Bool_t hasTZEROTOF = analysisEvent->GetTimeZeroTOFSigma()[9] < 150.;
    TZEROTOF += TOF_TZEROTOF_mean - TOF_mean;

    if (hasTZEROA) {
      hCentrality_TZEROA->Fill(cent, TZEROA);
      hVertex_TZEROA->Fill(vertexz, TZEROA);
    }
    if (hasTZEROC) {
      hCentrality_TZEROC->Fill(cent, TZEROC);
      hVertex_TZEROC->Fill(vertexz, TZEROC);
    }
    if (hasTZEROAND) {
      hCentrality_TZEROAND->Fill(cent, TZEROAND);
      hVertex_TZEROAND->Fill(vertexz, TZEROAND);
    }
    if (hasTZEROTOF) {
      hCentrality_TZEROTOF->Fill(cent, TZEROTOF);
      hVertex_TZEROTOF->Fill(vertexz, TZEROTOF);
    }

    if (hasTZEROA && hasTZEROC && hasTZEROAND && hasTZEROTOF) {
      hCentrality_TZEROTOF_TZEROA->Fill(cent, TZEROTOF - TZEROA);
      hCentrality_TZEROTOF_TZEROC->Fill(cent, TZEROTOF - TZEROC);
      hCentrality_TZEROTOF_TZEROAND->Fill(cent, TZEROTOF - TZEROAND);
      hCentrality_TZEROA_TZEROC->Fill(cent, TZEROA - TZEROC);
    }

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      AliAnalysisTrack *analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check accepted track */
      if (!analysisTrack->AcceptTrack()) continue;
      /* check momentum */
      if (analysisTrack->GetP() < 0.9 || analysisTrack->GetP() > 1.1) continue;
      /* check TOF pid */
      if (!analysisTrack->HasTOFPID()) continue;
      
      TZEROTOF = analysisEvent->GetTimeZeroTOF(analysisTrack->GetP());
      TZEROTOF += TOF_TZEROTOF_mean - TOF_mean;
      hasTZEROTOF = analysisEvent->GetTimeZeroTOFSigma(analysisTrack->GetP()) < 150.;

      Double_t time = analysisTrack->GetTOFTime() - TOF_mean;
      Double_t texp = analysisTrack->GetTOFExpTime(2);
      Double_t deltat = time - texp;

      
      hCentrality_TOF->Fill(cent, deltat);
      hVertex_TOF->Fill(vertexz, deltat);
      if (hasTZEROA) {
	hCentrality_TOF_TZEROA->Fill(cent, deltat - TZEROA);
	hVertex_TOF_TZEROA->Fill(vertexz, deltat - TZEROA);
      }
      if (hasTZEROC) {
	hCentrality_TOF_TZEROC->Fill(cent, deltat - TZEROC);
	hVertex_TOF_TZEROC->Fill(vertexz, deltat - TZEROC);
      }
      if (hasTZEROAND) {
	hCentrality_TOF_TZEROAND->Fill(cent, deltat - TZEROAND);
	hVertex_TOF_TZEROAND->Fill(vertexz, deltat - TZEROAND);
      }
      if (hasTZEROTOF) {
	hCentrality_TOF_TZEROTOF->Fill(cent, deltat - TZEROTOF);
	hVertex_TOF_TZEROTOF->Fill(vertexz, deltat - TZEROTOF);
      }

    }

  } /* end of loop over events */
  
  /* output */
  if (!calibfilename) TFile *fileout = TFile::Open(Form("TZEROcalib.%s", filename), "RECREATE");
  else TFile *fileout = TFile::Open(Form("TZEROcheck.%s", filename), "RECREATE");

  hCentrality->Write();
  hCentrality_TZEROA->Write();
  hCentrality_TZEROC->Write();
  hCentrality_TZEROAND->Write();
  hCentrality_TZEROTOF->Write();
  hCentrality_TZERODIFF->Write();

  hVertex_TZEROA->Write();
  hVertex_TZEROC->Write();
  hVertex_TZEROAND->Write();
  hVertex_TZEROTOF->Write();
  hVertex_TZERODIFF->Write();

  hCentrality_TOF->Write();
  hCentrality_TOF_TZEROA->Write();
  hCentrality_TOF_TZEROC->Write();
  hCentrality_TOF_TZEROAND->Write();
  hCentrality_TOF_TZEROTOF->Write();
  hCentrality_Resolution->Write();

  hCentrality_TZEROTOF_TZEROA->Write();
  hCentrality_TZEROTOF_TZEROC->Write();
  hCentrality_TZEROTOF_TZEROAND->Write();
  hCentrality_TZEROA_TZEROC->Write();

  hVertex_TOF->Write();
  hVertex_TOF_TZEROA->Write();
  hVertex_TOF_TZEROC->Write();
  hVertex_TOF_TZEROAND->Write();
  hVertex_TOF_TZEROTOF->Write();

  hResolution_TZERODIFF->Write();
  hResolution_TOF_TZEROA->Write();
  hResolution_TOF_TZEROC->Write();
  hResolution_TOF_TZEROAND->Write();

  hDiffDiff->Write();

  fileout->Close();

  if (!calibfilename) {
    TZEROcalibration(Form("TZEROcalib.%s", filename));
    TZEROcalib(filename, "TZEROcalibration.root");
    TZEROresolution(Form("TZEROcheck.%s", filename));
  }

}

TZEROcalibration(const Char_t *filename = "TZEROcalib.data.root", Float_t min = 2., Float_t max = 1.)
{

  gROOT->LoadMacro("FitPeak.C");

  TFile *fin = TFile::Open(filename);
  TH2 *hCentrality_TOF = (TH2 *)fin->Get("hCentrality_TOF");
  TH2 *hCentrality_TZEROA = (TH2 *)fin->Get("hCentrality_TZEROA");
  TH2 *hCentrality_TZEROC = (TH2 *)fin->Get("hCentrality_TZEROC");

  TH2 *hCentrality_TOF_TZEROTOF = (TH2 *)fin->Get("hCentrality_TOF_TZEROTOF");
  TH2 *hCentrality_TOF_TZEROA = (TH2 *)fin->Get("hCentrality_TOF_TZEROA");
  TH2 *hCentrality_TOF_TZEROC = (TH2 *)fin->Get("hCentrality_TOF_TZEROC");
  TH2 *hCentrality_TOF_TZEROTOF = (TH2 *)fin->Get("hCentrality_TOF_TZEROTOF");

  TFile *fout = TFile::Open("TZEROcalibration.root", "RECREATE");
  FitPeak(hCentrality_TOF, 200., 2., 1., "hCentrality_TOF")->Write();
  FitPeak(hCentrality_TZEROA, 200., 2., 1., "hCentrality_TZEROA")->Write();
  FitPeak(hCentrality_TZEROC, 200., 2., 1., "hCentrality_TZEROC")->Write();

  FitPeak(hCentrality_TOF_TZEROA, 200., 2., 1., "hCentrality_TOF_TZEROA")->Write();
  FitPeak(hCentrality_TOF_TZEROC, 200., 2., 1., "hCentrality_TOF_TZEROC")->Write();
  FitPeak(hCentrality_TOF_TZEROTOF, 200., 2., 1., "hCentrality_TOF_TZEROTOF")->Write();

  fout->Close();
}

TZEROresolution(const Char_t *filename, Float_t min = 2., Float_t max = 2., Char_t *what = "Centrality")
{

  TFile *fin = TFile::Open(filename);
  hdiff = (TH2 *)fin->Get(Form("h%s_TZEROA_TZEROC", what));
  hA = (TH2 *)fin->Get(Form("h%s_TZEROTOF_TZEROA", what));
  hC = (TH2 *)fin->Get(Form("h%s_TZEROTOF_TZEROC", what));
  hAND = (TH2 *)fin->Get(Form("h%s_TZEROTOF_TZEROAND", what));

  gROOT->LoadMacro("FitPeak.C");

  
  TH1 *hdiff_sigma = (TH1 *)FitPeak(hdiff, 100., min, max, "sigma1")->At(1);
  
  TH1 *hA_sigma = (TH1 *)FitPeak(hA, 100., min, max, "sigma2")->At(1);
  
  TH1 *hC_sigma = (TH1 *)FitPeak(hC, 100., min, max, "sigma3")->At(1);
  
  TH1 *hAND_sigma = (TH1 *)FitPeak(hAND, 100., min, max, "sigma4")->At(1);

  TH1 *hTZEROTOF_reso = (TH1 *)hdiff_sigma->Clone("hTZEROTOF_reso");
  hTZEROTOF_reso->Reset();
  hTZEROTOF_reso->SetMarkerStyle(20);
  hTZEROTOF_reso->SetMarkerColor(2);
  hTZEROTOF_reso->SetLineColor(2);
  hTZEROTOF_reso->SetTitle("T0-TOF");
  TH1 *hTZEROA_reso = (TH1 *)hdiff_sigma->Clone("hTZEROA_reso");
  hTZEROA_reso->Reset();
  hTZEROA_reso->SetMarkerStyle(21);
  hTZEROA_reso->SetMarkerColor(8);
  hTZEROA_reso->SetLineColor(8);
  hTZEROA_reso->SetTitle("T0-TZEROA");
  TH1 *hTZEROC_reso = (TH1 *)hdiff_sigma->Clone("hTZEROC_reso");
  hTZEROC_reso->Reset();
  hTZEROC_reso->SetMarkerStyle(22);
  hTZEROC_reso->SetMarkerColor(4);
  hTZEROC_reso->SetLineColor(4);
  hTZEROC_reso->SetTitle("T0-TZEROC");
  TH1 *hTZEROAND_reso = (TH1 *)hdiff_sigma->Clone("hTZEROAND_reso");
  hTZEROAND_reso->Reset();
  hTZEROAND_reso->SetMarkerStyle(23);
  hTZEROAND_reso->SetMarkerColor(kYellow+1);
  hTZEROAND_reso->SetLineColor(kYellow+1);
  hTZEROAND_reso->SetTitle("T0-TZEROAND");

  for (Int_t i = 0; i < hTZEROTOF_reso->GetNbinsX(); i++) {

    Double_t diff = hdiff_sigma->GetBinContent(i+1);
    Double_t diff_e = hdiff_sigma->GetBinError(i+1);
    Double_t diff2 = diff*diff;
    Double_t diff2_e = 2.*diff*diff_e;
    Double_t a = hA_sigma->GetBinContent(i+1);
    Double_t a_e = hA_sigma->GetBinError(i+1);
    Double_t a2 = a*a;
    Double_t a2_e = 2.*a*a_e;
    Double_t c = hC_sigma->GetBinContent(i+1);
    Double_t c_e = hC_sigma->GetBinError(i+1);
    Double_t c2 = c*c;
    Double_t c2_e = 2.*c*c_e;
    Double_t ac = hAND_sigma->GetBinContent(i+1);
    Double_t ac_e = hAND_sigma->GetBinError(i+1);
    Double_t ac2 = ac*ac;
    Double_t ac2_e = 2.*ac*ac_e;

    if (diff == 0.) continue;
    if (a2 + c2 - diff2 < 0.) continue;

    printf("a=%f, c=%f, diff=%f, ac=%f\n", a, c, diff, ac);

    Double_t tofreso2 = (a2 + c2 - diff2) * 0.5;
    Double_t tofreso2_e = (a2_e + c2_e + diff2_e) * 0.5;
    Double_t tofreso = TMath::Sqrt(tofreso2);
    Double_t tofreso_e = tofreso2_e / (2. * tofreso);
    hTZEROTOF_reso->SetBinContent(i+1, tofreso);
    hTZEROTOF_reso->SetBinError(i+1, tofreso_e);

    if (a2 - tofreso*tofreso < 0.) continue;
    if (c2 - tofreso*tofreso < 0.) continue;
    if (ac2 - tofreso*tofreso < 0.) continue;

    Double_t tzeroareso2 = a2 - tofreso2;
    Double_t tzeroareso2_e = a2_e + tofreso2_e;
    Double_t tzeroareso = TMath::Sqrt(tzeroareso2);
    Double_t tzeroareso_e = tzeroareso2_e / (2. * tzeroareso);
    hTZEROA_reso->SetBinContent(i+1, tzeroareso);
    hTZEROA_reso->SetBinError(i+1, tzeroareso_e);

    Double_t tzerocreso2 = c2 - tofreso2;
    Double_t tzerocreso2_e = c2_e + tofreso2_e;
    Double_t tzerocreso = TMath::Sqrt(tzerocreso2);
    Double_t tzerocreso_e = tzerocreso2_e / (2. * tzerocreso);
    hTZEROC_reso->SetBinContent(i+1, tzerocreso);
    hTZEROC_reso->SetBinError(i+1, tzerocreso_e);

    Double_t tzeroacreso2 = ac2 - tofreso2;
    Double_t tzeroacreso2_e = ac2_e + tofreso2_e;
    Double_t tzeroacreso = TMath::Sqrt(tzeroacreso2);
    Double_t tzeroacreso_e = tzeroacreso2_e / (2. * tzeroacreso);
    hTZEROAND_reso->SetBinContent(i+1, tzeroacreso);
    hTZEROAND_reso->SetBinError(i+1, tzeroacreso_e);

    //   Double_t tzerocreso = TMath::Sqrt(c2 - tofreso*tofreso);
    //    hTZEROC_reso->SetBinContent(i+1, tzerocreso);

    //    Double_t tzeroacreso = TMath::Sqrt(ac2 - tofreso*tofreso);
    //    hTZEROAND_reso->SetBinContent(i+1, tzeroacreso);

    printf("TOF=%f, TZEROA=%f, TZEROC=%f, TZEROAND=%f\n", tofreso, tzeroareso, tzerocreso, tzeroacreso);
  }

  hTZEROTOF_reso->Draw();
  hTZEROA_reso->Draw("same");
  hTZEROC_reso->Draw("same");
  hTZEROAND_reso->Draw("same");

  TFile *fout = TFile::Open("TZEROresolution.root", "RECREATE");
  hTZEROTOF_reso->Write();
  hTZEROA_reso->Write();
  hTZEROC_reso->Write();
  hTZEROAND_reso->Write();
  fout->Close();
}
