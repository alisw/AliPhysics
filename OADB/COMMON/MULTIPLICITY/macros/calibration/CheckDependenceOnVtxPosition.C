#include "AliMultEstimator.h"
#include "AliMultSelectionCuts.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include <TString.h>
#include <TProfile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TFile.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TLatex.h>

void StyleSettings( TString format = ""){
  //gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptDate(0);   //show day and time
  gStyle->SetOptStat(0);  //show statistic
  gStyle->SetPalette(1,0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTextSize(0.5);
  gStyle->SetLabelSize(0.03,"xyz");
  gStyle->SetLabelOffset(0.002,"xyz");
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);

  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.09);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.13);


  TGaxis::SetMaxDigits(3);
  gErrorIgnoreLevel=kError;

  if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
}

//__________________________________________________________________________________________________________
void SetStyleTLatex( TLatex* text,
                     Size_t textSize,
                     Width_t lineWidth,
                     Color_t textColor = 1,
                     Font_t textFont = 42,
                     Bool_t kNDC = kTRUE,
                     Short_t align = 11
){
  if (kNDC) {text->SetNDC();}
  text->SetTextFont(textFont);
  text->SetTextColor(textColor);
  text->SetTextSize(textSize);
  text->SetLineWidth(lineWidth);
  text->SetTextAlign(align);
}


//________________________________________________________________
void DrawTProfileWithCorrespondingFunction(
                                            TProfile* prof,
                                            TF1* fit            = NULL,
                                            TString xAxisName   ="",
                                            TString yAxisName   ="",
                                            TString label       = "",
                                            TString outputName  = ""
                                          ){
    TCanvas* canvas = new TCanvas("canvas1", "canvas1", 600, 450);
    canvas->SetLeftMargin(0.08);
    canvas->SetRightMargin(0.015);
    canvas->SetTopMargin(0.03);
    canvas->SetBottomMargin(0.08);
    canvas->SetFillColor(0);


    prof->SetMarkerStyle(24);
    prof->SetMarkerSize(1);
    prof->SetLineWidth(3);
    prof->SetMarkerColor(kBlack);
    prof->SetLineColor(kBlack);
    prof->GetYaxis()->SetLabelFont(42);
    prof->GetXaxis()->SetLabelFont(42);
    prof->GetYaxis()->SetTitleFont(62);
    prof->GetXaxis()->SetTitleFont(62);
    prof->GetXaxis()->SetTitleOffset(0.9);
    prof->GetYaxis()->SetTitle(yAxisName.Data());
    prof->GetXaxis()->SetTitle(xAxisName.Data());
    prof->GetYaxis()->SetRangeUser(prof->GetMinimum(0)*0.85, prof->GetMaximum()*1.15);
    prof->SetTitle("");
    prof->Draw("pe");

    if (fit != NULL){
      fit->SetLineColor(kBlue+1);
      fit->SetLineStyle(7);
      fit->SetLineWidth(10);
      fit->Draw("same");
    }
    TLatex* tlabel = new TLatex(0.95, 0.93, "");
    SetStyleTLatex(tlabel, 0.035, 1, kBlack, 42, kTRUE, 31);
    tlabel->DrawLatex(0.95, 0.91, label.Data());
    TString variableName[5] = {"a","b","c","d","e"};
    if (fit != NULL){
      for (Int_t i = 0; i < fit->GetNpar() && i < 5; i++){
        tlabel->DrawLatex(0.95, 0.90-(i+1)*0.035, Form("%s = %3.4f #pm %3.5f",variableName[i].Data(), fit->GetParameter(i), fit->GetParError(i) ));
      }
      tlabel->DrawLatex(0.95, 0.90-(fit->GetNpar()+1)*0.035, Form("#chi^{2}/ndf = %3.4f",fit->GetChisquare()/ fit->GetNDF() ));
    }
    canvas->SaveAs(outputName.Data());
    delete canvas;
}

//________________________________________________________________
void CheckDependenceOnVtxPosition(
                                    TString lPeriodName     = "",
                                    TString nameInputFile   = "",
                                    TString nameOutputFile  = "",
                                    Bool_t isLowStat        = kFALSE
                                  ) {

    StyleSettings("pdf");
    // Function meant to generate calibration OADB
    //
    // --- input : nameInputFile, containing a TTree object
    Int_t collSys         = 0;

    TString collisionSystem = lPeriodName;
    if ( lPeriodName.Contains("LHC16q") || lPeriodName.Contains("LHC16t") || lPeriodName.Contains("LHC13b") || lPeriodName.Contains("LHC13c") || lPeriodName.Contains("LHC13d") || lPeriodName.Contains("LHC13e") ){
      collisionSystem = collisionSystem+", p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
      collSys         = 1;
    } else if (lPeriodName.Contains("LHC13f") ) {
      collisionSystem = collisionSystem+", Pb-p #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
      collSys         = 1;
    } else if (lPeriodName.Contains("LHC16r") ){
      collisionSystem = collisionSystem+", p-Pb #sqrt{#it{s}_{_{NN}}} = 8.16 TeV";
      collSys         = 1;
    } else if (lPeriodName.Contains("LHC16s") ) {
      collisionSystem = collisionSystem+", Pb-p #sqrt{#it{s}_{_{NN}}} = 8.16 TeV";
      collSys         = 1;
    } else if (lPeriodName.Contains("LHC17n") ) {
      collisionSystem = collisionSystem+", Xe-Xe #sqrt{#it{s}_{_{NN}}} = 5.44 TeV";
      collSys         = 2;
    } else if (lPeriodName.Contains("LHC16f_lowB") || lPeriodName.Contains("LHC17g") || lPeriodName.Contains("LHC18c")){
      collisionSystem = collisionSystem+", pp #sqrt{#it{s}} = 13 TeV, B = 0.2 T";
    } else if (lPeriodName.Contains("LHC15n") || lPeriodName.Contains("LHC17p") || lPeriodName.Contains("LHC17q")){
      collisionSystem = collisionSystem+", pp #sqrt{#it{s}} = 5 TeV";
    }

    AliMultSelectionCalibrator *lCalib = new AliMultSelectionCalibrator("lCalib");
    if ( lPeriodName.Contains("LHC16q") || lPeriodName.Contains("LHC16t") || lPeriodName.Contains("LHC13b") || lPeriodName.Contains("LHC13c")) {
      cout<<"Setting event selection criteria for p-Pb..."<<endl;
      lCalib->GetEventCuts()->SetVzCut(10.0);
      lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
      lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
      lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
      lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
      lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
      lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    } else if ( lPeriodName.Contains("LHC16r") || lPeriodName.Contains("LHC16s")  || lPeriodName.Contains("LHC13d") || lPeriodName.Contains("LHC13e") || lPeriodName.Contains("LHC13f") ) {
      cout<<"Setting event selection criteria for HI p-Pb and Pb-p ..."<<endl;
      lCalib->GetEventCuts()->SetVzCut(10.0);
      lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
      lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
      lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
      lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kTRUE);
      lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
      lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    } else if ( lPeriodName.Contains("LHC17n")  ){
      cout<<"Setting event selection criteria for Xe-Xe..."<<endl;
      lCalib->GetEventCuts()->SetVzCut(10.0);
      lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
      lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
      lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE );
      lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
      lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
      lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE );
    } else {
      lCalib->GetEventCuts()->SetVzCut(10.0);
      lCalib->GetEventCuts()->SetTriggerCut                (kTRUE);
      lCalib->GetEventCuts()->SetINELgtZEROCut             (kTRUE);
      lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kTRUE);
      lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kTRUE);
      lCalib->GetEventCuts()->SetVertexConsistencyCut      (kTRUE);
      lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
      lCalib->GetEventCuts()->SetIsNotIncompleteDAQ        (kTRUE);
    }


    TString nameOutputDir = Form("CalibrationQA/%s", lPeriodName.Data());
    gSystem->Exec("mkdir -p "+nameOutputDir);

    // Steps involved:
    //  (1) Set up basic I/O
    //  (2) Detect Runs From Input File
    //  (3) Create runwise histograms
    //  (4) Write runwise histograms to file and merge them
    //  (5) Fit histograms

    cout<<"=== START Reading tree ==="<<endl;
    cout<<" * Input File.....: "<<nameInputFile.Data()<<endl;
    cout<<endl;
    cout<<" * Event Selection Peformed: "<<endl;
    lCalib->GetEventCuts() -> Print();
    cout<<endl;

    // STEP 1: Basic I/O
    cout<<"(1) Opening File"<<endl;

    //Open File
    TFile *fInputFile = TFile::Open( nameInputFile.Data(), "READ");
    if(!fInputFile) {
      cout << "File "<<   nameInputFile.Data() << " not found! " <<   endl;
        return ;
    }
    //Locate TTree object
    TTree* fTree = (TTree*)fInputFile->FindObjectAny("fTreeEvent");
    if(!fTree) {
        cout << "fTreeEvent object not found!" << endl;
        return ;
    }

    //Event Selection Variables
    Bool_t fEvSel_IsNotPileupInMultBins     = kFALSE ;
    Bool_t fEvSel_Triggered                 = kFALSE ;
    Bool_t fEvSel_INELgtZERO                = kFALSE ;
    Bool_t fEvSel_PassesTrackletVsCluster   = kFALSE ;
    Bool_t fEvSel_HasNoInconsistentVertices = kFALSE ;
    Bool_t fEvSel_IsNotAsymmetricInVZERO    = kFALSE ;
    Bool_t fEvSel_IsNotIncompleteDAQ        = kFALSE ;
    Bool_t fEvSel_HasGoodVertex2016         = kFALSE ;
    Int_t fRunNumber;
    Float_t fAmplitude_V0A                  = 0.;
    Float_t fAmplitude_V0C                  = 0.;
    Float_t fAmplitude_V0AOnline            = 0.;
    Float_t fAmplitude_V0COnline            = 0.;
    Float_t fAmplitude_V0AEq                = 0.;
    Float_t fAmplitude_V0CEq                = 0.;
    Float_t fAmplitude_ADA                  = 0.;
    Float_t fAmplitude_ADC                  = 0.;
    Int_t fnSPDClusters0                    = 0;
    Int_t fnSPDClusters1                    = 0;
    Int_t fRefMultEta5                      = 0;
    Int_t fRefMultEta8                      = 0;
    Int_t fnTracklets                       = 0;
    Int_t fnSPDClusters                     = 0;
    Float_t fEvSel_VtxZ                     = 0.;
    Float_t fZnaTower                       = 0.;
    Float_t fZncTower                       = 0.;
    Int_t fZnaFired                         = 0;
    Int_t fZncFired                         = 0;

    //FIXME/CAUTION: non-zero if using tree without that branch
    Int_t fnContributors = 1000;
    UInt_t fEvSel_TriggerMask; //! save full info for checking later
    Bool_t fCheckTriggerType = kTRUE;
    AliVEvent::EOfflineTriggerTypes fTrigType = AliVEvent::kINT7;
    if (!lPeriodName.CompareTo("LHC17n"))
      fCheckTriggerType = kFALSE;

    //SetBranchAddresses for event Selection Variables
    //(multiplicity related will be done automatically!)
    fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTree->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTree->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTree->SetBranchAddress("fEvSel_TriggerMask",&fEvSel_TriggerMask);
    fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTree->SetBranchAddress("fRunNumber",&fRunNumber);
    fTree->SetBranchAddress("fnContributors", &fnContributors);
    fTree->SetBranchAddress("fEvSel_IsNotAsymmetricInVZERO", &fEvSel_IsNotAsymmetricInVZERO);
    fTree->SetBranchAddress("fEvSel_IsNotIncompleteDAQ", &fEvSel_IsNotIncompleteDAQ);
    fTree->SetBranchAddress("fEvSel_HasGoodVertex2016", &fEvSel_HasGoodVertex2016);
    fTree->SetBranchAddress("fEvSel_VtxZ", &fEvSel_VtxZ);

    fTree->SetBranchAddress("fAmplitude_V0A",&fAmplitude_V0A);
    fTree->SetBranchAddress("fAmplitude_V0C",&fAmplitude_V0C);
    fTree->SetBranchAddress("fAmplitude_OnlineV0A",&fAmplitude_V0AOnline);
    fTree->SetBranchAddress("fAmplitude_OnlineV0C",&fAmplitude_V0COnline);
    fTree->SetBranchAddress("fAmplitude_V0AEq",&fAmplitude_V0AEq);
    fTree->SetBranchAddress("fAmplitude_V0CEq",&fAmplitude_V0CEq);
    fTree->SetBranchAddress("fAmplitude_ADA",&fAmplitude_ADA);
    fTree->SetBranchAddress("fAmplitude_ADC",&fAmplitude_ADC);
    fTree->SetBranchAddress("fnSPDClusters0",&fnSPDClusters0);
    fTree->SetBranchAddress("fnSPDClusters1",&fnSPDClusters1);
    fTree->SetBranchAddress("fnSPDClusters",&fnSPDClusters);
    fTree->SetBranchAddress("fRefMultEta5",&fRefMultEta5);
    fTree->SetBranchAddress("fRefMultEta8",&fRefMultEta8);
    fTree->SetBranchAddress("fnTracklets",&fnTracklets);

    fTree->SetBranchAddress("fZnaFired",&fZnaFired);
    fTree->SetBranchAddress("fZnaTower",&fZnaTower);
    fTree->SetBranchAddress("fZncFired",&fZncFired);
    fTree->SetBranchAddress("fZncTower",&fZncTower);


    //============================================================
    // Read input
    //============================================================


    Int_t maxAmplitudeV0A     = 1000;
    Int_t maxAmplitudeADA     = 1000;
    Int_t maxSPDCl            = 500;
    Int_t maxRefMult          = 200;
    Int_t maxNTracklets       = 500;
    Int_t maxZN               = 2000;
    Int_t nBinsZvtx           = 120;
    if (!lPeriodName.CompareTo("LHC17n") || !lPeriodName.CompareTo("LHC18d2")){
      maxAmplitudeV0A     = 1000*30;
      maxSPDCl            = 500*10;
      maxRefMult          = 200*10;
      maxNTracklets       = 500*8;
      maxZN               = 2000*10;
    }
    if (isLowStat)
      nBinsZvtx               = 60;

    Long64_t lNEv = fTree->GetEntries();
    cout<<"(1) File opened, event count is "<<lNEv<<endl;

    cout<<"(2) Creating buffer, computing averages"<<endl;
    const int lMax = 1000;
    const int lMaxQuantiles = 10000;
    Int_t lRunNumbers[lMaxQuantiles];
    Long_t lRunStats[lMaxQuantiles];
    Int_t fNRunRanges = 0;

    for( Int_t ix=0; ix<1000;ix++){
        lRunNumbers[ix] = 0;
        lRunStats[ix] = 0;
    }

    auto h2V0MvsNTracklets      = new TH2F("h2V0MvsNTracklets","h2V0MvsNTracklets",1000,0,maxAmplitudeV0A,maxNTracklets,0,maxNTracklets);
    auto h2V0AvsVtxZ            = new TH2F("h2V0AvsVtxZ","h2V0AvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0CvsVtxZ            = new TH2F("h2V0CvsVtxZ","h2V0CvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0MvsVtxZ            = new TH2F("h2V0MvsVtxZ","h2V0MvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0AOnlinevsVtxZ      = new TH2F("h2V0AOnlinevsVtxZ","h2V0AOnlinevsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0COnlinevsVtxZ      = new TH2F("h2V0COnlinevsVtxZ","h2V0COnlinevsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0MOnlinevsVtxZ      = new TH2F("h2V0MOnlinevsVtxZ","h2V0MOnlinevsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2ADAvsVtxZ            = new TH2F("h2ADAvsVtxZ","h2ADAvsVtxZ",1000,0,maxAmplitudeADA,nBinsZvtx,-12,12);
    auto h2ADCvsVtxZ            = new TH2F("h2ADCvsVtxZ","h2ADCvsVtxZ",1000,0,maxAmplitudeADA,nBinsZvtx,-12,12);
    auto h2ADMvsVtxZ            = new TH2F("h2ADMvsVtxZ","h2ADMvsVtxZ",1000,0,maxAmplitudeADA,nBinsZvtx,-12,12);
    auto h2V0AEqvsVtxZ          = new TH2F("h2V0AEqvsVtxZ","h2V0AEqvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0CEqvsVtxZ          = new TH2F("h2V0CEqvsVtxZ","h2V0CEqvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2V0MEqvsVtxZ          = new TH2F("h2V0MEqvsVtxZ","h2V0MEqvsVtxZ",1000,0,maxAmplitudeV0A,nBinsZvtx,-12,12);
    auto h2SPDCl0vsVtxZ         = new TH2F("h2SPDCl0vsVtxZ","h2SPDCl0vsVtxZ",maxSPDCl,0,maxSPDCl,nBinsZvtx,-12,12);
    auto h2SPDCl1vsVtxZ         = new TH2F("h2SPDCl1vsVtxZ","h2SPDCl1vsVtxZ",maxSPDCl,0,maxSPDCl,nBinsZvtx,-12,12);
    auto h2SPDClvsVtxZ          = new TH2F("h2SPDClvsVtxZ","h2SPDClvsVtxZ",maxSPDCl,0,maxSPDCl,nBinsZvtx,-12,12);
    auto h2RefMultEta5vsVtxZ    = new TH2F("h2RefMultEta5vsVtxZ","h2RefMultEta5vsVtxZ",maxRefMult,0,maxRefMult,nBinsZvtx,-12,12);
    auto h2RefMultEta8vsVtxZ    = new TH2F("h2RefMultEta8vsVtxZ","h2RefMultEta8vsVtxZ",maxRefMult,0,maxRefMult,nBinsZvtx,-12,12);
    auto h2NTrackletsvsVtxZ     = new TH2F("h2NTrackletsvsVtxZ","h2NTrackletsvsVtxZ",maxNTracklets,0,maxNTracklets,nBinsZvtx,-12,12);
    auto h2ZNAvsVtxZ            = new TH2F("h2ZNAvsVtxZ","h2ZNAvsVtxZ",maxZN,-10,maxZN,60,-12,12);
    auto h2ZNCvsVtxZ            = new TH2F("h2ZNAvsVtxZ","h2ZNAvsVtxZ",maxZN,-10,maxZN,60,-12,12);

    auto hprofVtxZvsV0A         = new TProfile("hprofVtxZvsV0A","hprofVtxZvsV0A",nBinsZvtx,-12,12, 0,maxAmplitudeV0A);
    auto hprofVtxZvsV0C         = new TProfile("hprofVtxZvsV0C","hprofVtxZvsV0C",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0M         = new TProfile("hprofVtxZvsV0M","hprofVtxZvsV0M",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0AOnline   = new TProfile("hprofVtxZvsV0AOnline","hprofVtxZvsV0AOnline",nBinsZvtx,-12,12, 0,maxAmplitudeV0A);
    auto hprofVtxZvsV0COnline   = new TProfile("hprofVtxZvsV0COnline","hprofVtxZvsV0COnline",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0MOnline   = new TProfile("hprofVtxZvsV0MOnline","hprofVtxZvsV0MOnline",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0AEq       = new TProfile("hprofVtxZvsV0AEq","hprofVtxZvsV0AEq",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0CEq       = new TProfile("hprofVtxZvsV0CEq","hprofVtxZvsV0CEq",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsV0MEq       = new TProfile("hprofVtxZvsV0MEq","hprofVtxZvsV0MEq",nBinsZvtx,-12,12,0,maxAmplitudeV0A);
    auto hprofVtxZvsADA         = new TProfile("hprofVtxZvsADA","hprofVtxZvsADA",nBinsZvtx,-12,12, 0,maxAmplitudeADA);
    auto hprofVtxZvsADC         = new TProfile("hprofVtxZvsADC","hprofVtxZvsADC",nBinsZvtx,-12,12,0,maxAmplitudeADA);
    auto hprofVtxZvsADM         = new TProfile("hprofVtxZvsADM","hprofVtxZvsADM",nBinsZvtx,-12,12,0,maxAmplitudeADA);
    auto hprofVtxZvsSPDCl0      = new TProfile("hprofVtxZvsSPDCl0","hprofVtxZvsSPDCl0",nBinsZvtx,-12,12,0,maxSPDCl);
    auto hprofVtxZvsSPDCl1      = new TProfile("hprofVtxZvsSPDCl1","hprofVtxZvsSPDCl1",nBinsZvtx,-12,12,0,maxSPDCl);
    auto hprofVtxZvsSPDCl       = new TProfile("hprofVtxZvsSPDCl","hprofVtxZvsSPDCl",nBinsZvtx,-12,12,0,maxSPDCl);
    auto hprofVtxZvsRefMultEta5 = new TProfile("hprofVtxZvsRefMultEta5","hprofVtxZvsRefMultEta5",nBinsZvtx,-12,12,0,maxRefMult);
    auto hprofVtxZvsRefMultEta8 = new TProfile("hprofVtxZvsRefMultEta8","hprofVtxZvsRefMultEta8",nBinsZvtx,-12,12,0,maxRefMult);
    auto hprofVtxZvsNTracklets  = new TProfile("hprofVtxZvsNTracklets","hprofVtxZvsNTracklets",nBinsZvtx,-12,12,0,maxNTracklets);
    auto hprofVtxZvsZNA         = new TProfile("hprofVtxZvsZNA","hprofVtxZvsZNA",(Int_t)(nBinsZvtx/2),-12,12,-10,maxZN);
    auto hprofVtxZvsZNC         = new TProfile("hprofVtxZvsZNC","hprofVtxZvsZNC",(Int_t)(nBinsZvtx/2),-12,12,-10,maxZN);

    Int_t lNRuns = 0;
    Bool_t lNewRun = kTRUE;
    Int_t lThisRunIndex = -1;

    for(Long64_t iEv = 0; iEv<fTree->GetEntries(); iEv++) {
        if ( iEv % 100000 == 0 ) {
            Double_t complete = 100. * ( double ) ( iEv ) / ( double ) ( fTree->GetEntries() );
            cout << "Event # " << iEv << "/" << fTree->GetEntries() << " (" << complete << "%)" <<endl;;
        }
        fTree->GetEntry(iEv);
        //Perform Event selection
        Bool_t lSaveThisEvent = kTRUE; //let's be optimistic

        //Apply trigger mask (will only work if not kAny)
        Bool_t isSelected = 0;
        isSelected = fEvSel_TriggerMask & fTrigType;
        if(!isSelected && fCheckTriggerType) lSaveThisEvent = kFALSE;

        //Check Selections as they are in the fMultSelectionCuts Object
        if( lCalib->GetEventCuts()->GetTriggerCut()    && ! fEvSel_Triggered  ) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetINELgtZEROCut() && ! fEvSel_INELgtZERO ) lSaveThisEvent = kFALSE;
        if( TMath::Abs( fEvSel_VtxZ ) > lCalib->GetEventCuts()->GetVzCut()      ) lSaveThisEvent = kFALSE;
        //ADD ME HERE: Tracklets Vs Clusters Cut?
        if( lCalib->GetEventCuts()->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins    ) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster  ) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetNonZeroNContribs()          &&  fnContributors < 1 ) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetIsNotAsymmetricInVZERO()    && ! fEvSel_IsNotAsymmetricInVZERO) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetIsNotIncompleteDAQ()        && ! fEvSel_IsNotIncompleteDAQ) lSaveThisEvent = kFALSE;
        if( lCalib->GetEventCuts()->GetHasGoodVertex2016()         && ! fEvSel_HasGoodVertex2016) lSaveThisEvent = kFALSE;

        Int_t lIndex = -1;
        lNewRun = kTRUE;
        for(Int_t iRun=0; iRun<lNRuns; iRun++) {
            if( lRunNumbers[iRun] == fRunNumber ) {
                lNewRun = kFALSE;
                lIndex = iRun;
            }
        }
        if( lNewRun == kTRUE ) {
            cout<<"(Autodiscover) New Run Found: "<<fRunNumber<<", added as #"<<lNRuns<<" (so far: "<<lNRuns<<" runs)"<<endl;

              lRunNumbers[lNRuns] = fRunNumber;
              lIndex = lNRuns;
              lNRuns++;
              fNRunRanges++;
          }
        if ( lSaveThisEvent ) {
          h2V0AvsVtxZ->Fill(fAmplitude_V0A, fEvSel_VtxZ);
          h2V0CvsVtxZ->Fill(fAmplitude_V0C, fEvSel_VtxZ);
          h2V0MvsVtxZ->Fill(fAmplitude_V0A+fAmplitude_V0C, fEvSel_VtxZ);
          h2V0AOnlinevsVtxZ->Fill(fAmplitude_V0AOnline, fEvSel_VtxZ);
          h2V0COnlinevsVtxZ->Fill(fAmplitude_V0COnline, fEvSel_VtxZ);
          h2V0MOnlinevsVtxZ->Fill(fAmplitude_V0AOnline+fAmplitude_V0COnline, fEvSel_VtxZ);
          h2V0AEqvsVtxZ->Fill(fAmplitude_V0AEq, fEvSel_VtxZ);
          h2V0CEqvsVtxZ->Fill(fAmplitude_V0CEq, fEvSel_VtxZ);
          h2V0CEqvsVtxZ->Fill(fAmplitude_V0CEq+fAmplitude_V0AEq, fEvSel_VtxZ);
          h2ADAvsVtxZ->Fill(fAmplitude_ADA, fEvSel_VtxZ);
          h2ADCvsVtxZ->Fill(fAmplitude_ADC, fEvSel_VtxZ);
          h2ADMvsVtxZ->Fill(fAmplitude_ADA+fAmplitude_ADC, fEvSel_VtxZ);
          h2SPDCl0vsVtxZ->Fill(fnSPDClusters0, fEvSel_VtxZ);
          h2SPDCl1vsVtxZ->Fill(fnSPDClusters1, fEvSel_VtxZ);
          h2SPDClvsVtxZ->Fill(fnSPDClusters, fEvSel_VtxZ);
          h2RefMultEta5vsVtxZ->Fill(fRefMultEta5, fEvSel_VtxZ);
          h2RefMultEta8vsVtxZ->Fill(fRefMultEta8, fEvSel_VtxZ);
          h2NTrackletsvsVtxZ->Fill(fnTracklets, fEvSel_VtxZ);
          h2ZNAvsVtxZ->Fill((Bool_t)fZnaFired*fZnaTower+!((Bool_t)fZnaFired)*0, fEvSel_VtxZ);
          h2ZNCvsVtxZ->Fill((Bool_t)fZncFired*fZncTower+!((Bool_t)fZncFired)*0, fEvSel_VtxZ);

          h2V0MvsNTracklets->Fill(fAmplitude_V0A+fAmplitude_V0C, fnTracklets);

          hprofVtxZvsV0A->Fill(fEvSel_VtxZ, fAmplitude_V0A);
          hprofVtxZvsV0C->Fill(fEvSel_VtxZ, fAmplitude_V0C);
          hprofVtxZvsV0M->Fill(fEvSel_VtxZ, fAmplitude_V0A+fAmplitude_V0C);
          hprofVtxZvsV0AOnline->Fill(fEvSel_VtxZ, fAmplitude_V0AOnline);
          hprofVtxZvsV0COnline->Fill(fEvSel_VtxZ, fAmplitude_V0COnline);
          hprofVtxZvsV0MOnline->Fill(fEvSel_VtxZ, fAmplitude_V0AOnline+fAmplitude_V0COnline);
          hprofVtxZvsV0AEq->Fill(fEvSel_VtxZ, fAmplitude_V0AEq);
          hprofVtxZvsV0CEq->Fill(fEvSel_VtxZ, fAmplitude_V0CEq);
          hprofVtxZvsV0MEq->Fill(fEvSel_VtxZ, fAmplitude_V0CEq+fAmplitude_V0AEq);
          hprofVtxZvsADA->Fill(fEvSel_VtxZ, fAmplitude_ADA);
          hprofVtxZvsADC->Fill(fEvSel_VtxZ, fAmplitude_ADC);
          hprofVtxZvsADM->Fill(fEvSel_VtxZ, fAmplitude_ADA+fAmplitude_ADC);
          hprofVtxZvsSPDCl0->Fill(fEvSel_VtxZ, fnSPDClusters0);
          hprofVtxZvsSPDCl1->Fill(fEvSel_VtxZ, fnSPDClusters1);
          hprofVtxZvsSPDCl->Fill(fEvSel_VtxZ, fnSPDClusters);
          hprofVtxZvsRefMultEta5->Fill(fEvSel_VtxZ, fRefMultEta5);
          hprofVtxZvsRefMultEta8->Fill(fEvSel_VtxZ, fRefMultEta8);
          hprofVtxZvsNTracklets->Fill(fEvSel_VtxZ, fnTracklets);

          if (collSys > 0){
            hprofVtxZvsZNA->Fill(fEvSel_VtxZ, (Bool_t)fZnaFired*fZnaTower+!((Bool_t)fZnaFired)*0);
            hprofVtxZvsZNC->Fill(fEvSel_VtxZ, (Bool_t)fZncFired*fZncTower+!((Bool_t)fZncFired)*0);
          } else {
            hprofVtxZvsZNA->Fill(fEvSel_VtxZ, -(Bool_t)fZnaFired*fZnaTower+!((Bool_t)fZnaFired)*1e6);
            hprofVtxZvsZNC->Fill(fEvSel_VtxZ, -(Bool_t)fZncFired*fZncTower+!((Bool_t)fZncFired)*1e6);
          }
          //           sTree [ lIndex ] -> Fill();
        }

    }

    TF1* fitVtxZvsV0A = new TF1("fitVtxZvsV0A","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0A->Fit(fitVtxZvsV0A,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0C = new TF1("fitVtxZvsV0C","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0C->Fit(fitVtxZvsV0C,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0M = new TF1("fitVtxZvsV0M","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0M->Fit(fitVtxZvsV0M,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0AOnline = new TF1("fitVtxZvsV0AOnline","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0AOnline->Fit(fitVtxZvsV0AOnline,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0COnline = new TF1("fitVtxZvsV0COnline","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0COnline->Fit(fitVtxZvsV0COnline,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0MOnline = new TF1("fitVtxZvsV0MOnline","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0MOnline->Fit(fitVtxZvsV0MOnline,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0AEq = new TF1("fitVtxZvsV0AEq","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0AEq->Fit(fitVtxZvsV0AEq,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0CEq = new TF1("fitVtxZvsV0CEq","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0CEq->Fit(fitVtxZvsV0CEq,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsV0MEq = new TF1("fitVtxZvsV0MEq","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsV0MEq->Fit(fitVtxZvsV0MEq,"Q0EMRN","",-10,10);

    TF1* fitVtxZvsADA = new TF1("fitVtxZvsADA","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsADA->Fit(fitVtxZvsADA,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsADC = new TF1("fitVtxZvsADC","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsADC->Fit(fitVtxZvsADC,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsADM = new TF1("fitVtxZvsADM","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsADM->Fit(fitVtxZvsADM,"Q0EMRN","",-10,10);

    TF1* fitVtxZvsSPDCl0 = new TF1("fitVtxZvsSPDCl0","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsSPDCl0->Fit(fitVtxZvsSPDCl0,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsSPDCl1 = new TF1("fitVtxZvsSPDCl1","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsSPDCl1->Fit(fitVtxZvsSPDCl1,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsSPDCl = new TF1("fitVtxZvsSPDCl","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsSPDCl->Fit(fitVtxZvsSPDCl,"Q0EMRN","",-10,10);

    TF1* fitVtxZvsRefMultEta5 = new TF1("fitVtxZvsRefMultEta5","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsRefMultEta5->Fit(fitVtxZvsRefMultEta5,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsRefMultEta8 = new TF1("fitVtxZvsRefMultEta8","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsRefMultEta8->Fit(fitVtxZvsRefMultEta8,"Q0EMRN","",-10,10);

    TF1* fitVtxZvsNTracklets = new TF1("fitVtxZvsNTracklets","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsNTracklets->Fit(fitVtxZvsNTracklets,"Q0EMRN","",-10,10);

    TF1* fitVtxZvsZNA = new TF1("fitVtxZvsZNA","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsZNA->Fit(fitVtxZvsZNA,"Q0EMRN","",-10,10);
    TF1* fitVtxZvsZNC = new TF1("fitVtxZvsZNC","[4]+[3]*x+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",-10,10);
    hprofVtxZvsZNC->Fit(fitVtxZvsZNC,"Q0EMRN","",-10,10);

    //Write buffer to file
//     for(Int_t iRun=0; iRun<lNRuns; iRun++) sTree[iRun]->Write();
    TFile *fOutput = new TFile (nameOutputFile.Data(), "RECREATE");
      h2V0MvsNTracklets->Write();
      h2V0AvsVtxZ->Write();
      h2V0CvsVtxZ->Write();
      h2V0MvsVtxZ->Write();
      h2V0AOnlinevsVtxZ->Write();
      h2V0COnlinevsVtxZ->Write();
      h2V0MOnlinevsVtxZ->Write();
      h2V0AEqvsVtxZ->Write();
      h2V0CEqvsVtxZ->Write();
      h2V0MEqvsVtxZ->Write();
      h2ADAvsVtxZ->Write();
      h2ADCvsVtxZ->Write();
      h2ADMvsVtxZ->Write();
      h2SPDCl0vsVtxZ->Write();
      h2SPDCl1vsVtxZ->Write();
      h2SPDClvsVtxZ->Write();
      h2RefMultEta5vsVtxZ->Write();
      h2RefMultEta8vsVtxZ->Write();
      h2NTrackletsvsVtxZ->Write();
      h2ZNAvsVtxZ->Write();
      h2ZNCvsVtxZ->Write();

      hprofVtxZvsV0A->Write();
      hprofVtxZvsV0C->Write();
      hprofVtxZvsV0M->Write();
      hprofVtxZvsV0AOnline->Write();
      hprofVtxZvsV0COnline->Write();
      hprofVtxZvsV0MOnline->Write();
      hprofVtxZvsV0AEq->Write();
      hprofVtxZvsV0CEq->Write();
      hprofVtxZvsV0MEq->Write();
      hprofVtxZvsADA->Write();
      hprofVtxZvsADC->Write();
      hprofVtxZvsADM->Write();
      hprofVtxZvsSPDCl0->Write();
      hprofVtxZvsSPDCl1->Write();
      hprofVtxZvsSPDCl->Write();
      hprofVtxZvsRefMultEta5->Write();
      hprofVtxZvsRefMultEta8->Write();
      hprofVtxZvsNTracklets->Write();
      hprofVtxZvsZNA->Write();
      hprofVtxZvsZNC->Write();

      fitVtxZvsV0A->Write();
      fitVtxZvsV0C->Write();
      fitVtxZvsV0M->Write();
      fitVtxZvsV0AOnline->Write();
      fitVtxZvsV0COnline->Write();
      fitVtxZvsV0MOnline->Write();
      fitVtxZvsV0AEq->Write();
      fitVtxZvsV0CEq->Write();
      fitVtxZvsV0MEq->Write();
      fitVtxZvsADA->Write();
      fitVtxZvsADC->Write();
      fitVtxZvsADM->Write();
      fitVtxZvsSPDCl0->Write();
      fitVtxZvsSPDCl1->Write();
      fitVtxZvsSPDCl->Write();
      fitVtxZvsRefMultEta5->Write();
      fitVtxZvsRefMultEta8->Write();
      fitVtxZvsNTracklets->Write();
      fitVtxZvsZNA->Write();
      fitVtxZvsZNC->Write();

      fOutput->Write();
    fOutput->Close();

    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0A, fitVtxZvsV0A, "Z_{vtx} (cm)", "V0A signal (arb. units)", collisionSystem, Form("%s/V0AZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0C, fitVtxZvsV0C, "Z_{vtx} (cm)", "V0C signal (arb. units)", collisionSystem, Form("%s/V0CZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0M, fitVtxZvsV0M, "Z_{vtx} (cm)", "V0M signal (arb. units)", collisionSystem, Form("%s/V0MZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0AOnline, fitVtxZvsV0AOnline, "Z_{vtx} (cm)", "V0A signal online (arb. units)", collisionSystem, Form("%s/V0AOnlineZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0COnline, fitVtxZvsV0COnline, "Z_{vtx} (cm)", "V0C signal online (arb. units)", collisionSystem, Form("%s/V0COnlineZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0MOnline, fitVtxZvsV0MOnline, "Z_{vtx} (cm)", "V0M signal online (arb. units)", collisionSystem, Form("%s/V0MOnlineZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0AEq, fitVtxZvsV0AEq, "Z_{vtx} (cm)", "V0AEq signal (arb. units)", collisionSystem, Form("%s/V0AEqZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0CEq, fitVtxZvsV0CEq, "Z_{vtx} (cm)", "V0CEq signal (arb. units)", collisionSystem, Form("%s/V0CEqZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsV0MEq, fitVtxZvsV0MEq, "Z_{vtx} (cm)", "V0MEq signal (arb. units)", collisionSystem, Form("%s/V0MEqZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsADA, fitVtxZvsADA, "Z_{vtx} (cm)", "ADA signal (arb. units)", collisionSystem, Form("%s/ADAZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsADC, fitVtxZvsADC, "Z_{vtx} (cm)", "ADC signal (arb. units)", collisionSystem, Form("%s/ADCZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsADM, fitVtxZvsADM, "Z_{vtx} (cm)", "ADM signal (arb. units)", collisionSystem, Form("%s/ADMZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsSPDCl0, fitVtxZvsSPDCl0, "Z_{vtx} (cm)", "SPD cl. layer 0", collisionSystem, Form("%s/SPCCl0ZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsSPDCl1, fitVtxZvsSPDCl1, "Z_{vtx} (cm)", "SPD cl. layer 1", collisionSystem, Form("%s/SPCCl1ZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsSPDCl, fitVtxZvsSPDCl, "Z_{vtx} (cm)", "SPD cl.", collisionSystem, Form("%s/SPCClZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsRefMultEta5, fitVtxZvsRefMultEta5, "Z_{vtx} (cm)", "mult |#eta| < 0.5", collisionSystem, Form("%s/RefMult05ZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsRefMultEta8, fitVtxZvsRefMultEta8, "Z_{vtx} (cm)", "mult |#eta| < 0.8", collisionSystem, Form("%s/RefMult08ZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsNTracklets, fitVtxZvsNTracklets, "Z_{vtx} (cm)", "SPD tracklets", collisionSystem, Form("%s/TrackletsZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsZNA, fitVtxZvsZNA, "Z_{vtx} (cm)", "ZNA signal (arb. units)", collisionSystem, Form("%s/ZNAZVtzDep.pdf", nameOutputDir.Data()));
    DrawTProfileWithCorrespondingFunction( hprofVtxZvsZNC, fitVtxZvsZNC, "Z_{vtx} (cm)", "ZNC signal (arb. units)", collisionSystem, Form("%s/ZNCZVtzDep.pdf", nameOutputDir.Data()));

//     cout<<"(3) Inspect Runs and corresponding statistics: "<<endl;
//     for(Int_t iRun = 0; iRun<fNRunRanges; iRun++) {
//       cout<<" --- Run #"<<iRun<<", (#"<<lRunNumbers[iRun]<<")"<<endl; //, N(events) = "<<sTree[iRun]->GetEntries()<<endl;
//     }
//     cout<<endl;
}
