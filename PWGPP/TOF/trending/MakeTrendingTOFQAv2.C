/*
 *  fbellini@cern.ch - last update on 18/10/2016
 *  Macro to run the TOF QA trending by accessing the std QA output,
 *  to be mainly used with the automatic scripts to fill the QA repository.
 *  Launch with
 *  aliroot -l -b -q "MakeTrendingTOFQA.C(\"${fullpath}/QAresults.root\", ${run}, ...)
 *  The macro produces a file containing the tree of trending variables and the main plots.
 *  A feature that displays the plots in canvases must be enable when needed.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFcalib.h"
#include "AliTOFcalibHisto.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGrid.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNamed.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include <iostream>
using std::cout;
using std::endl;
#endif

Int_t MakeTrendingTOFQAv2(const TString qafilename, //full path of the QA output;
    Int_t runNumber,                                //run number
    const TString dirsuffix = "",                   //suffix for subdirectories
    const Bool_t isMC = kFALSE,                     //MC flag, to disable meaningless checks
    const Bool_t checkPIDqa = kTRUE,                //set to kTRUE to check PIDqa output for TOF
    const Bool_t fitSignalModel = kTRUE,            //set to kTRUE to perform fit with TOF signal model too for PID
    const Double_t RangeFitNsigmaPIDmin = -2.,      //set lower limit for fitting Nsigma_TOF_ID
    const Double_t RangeFitNsigmaPIDmax = 2.,       //set upper limit for fitting Nsigma_TOF_ID
    const Double_t RangeTrksForTOFResMin = 10.,     //set lower limit to the number of tracks requested to extract the TOF resolution
    const Double_t RangeTrksForTOFResMax = 100.,    //set upper limit to the number of tracks requested to extract the TOF resolution
    const TString ocdbStorage = "raw://",           //set the default ocdb storage
    const Bool_t drawAll = kFALSE,                  //enable display plots on canvas and save png
    const Bool_t saveHisto = kTRUE,                 //set to kTRUE to save histograms in root file
    const Bool_t savePng = kTRUE,                   //set to kTRUE to save histogram to png image
    const Bool_t isAutoTrend = kFALSE);             //set to kTRUE for automatic trending

///
/// Function to get the number of active TOF channels
Double_t GetGoodTOFChannelsRatio(Int_t run = -1, Bool_t saveMap = kFALSE, TString OCDBstorage = "raw://", Bool_t inEta08 = kFALSE);

///
/// Function to compute the efficiency with 2D histograms
void Compute2Deff(TH2F* num, TH2F* den, TH1F*& h, TString name, Bool_t x = kTRUE);

///
/// Function to compute the efficiency with errors
std::pair<Double_t, Double_t> ComputeEff(Double_t num, Double_t den, Double_t numE, Double_t denE);

///
/// Function to get a TDirectoryFile
Int_t GetTDir(TDirectoryFile*& d, TFile*& f, TString name);

///
/// Function to get a TList
Int_t GetList(TDirectoryFile* d, TList*& l, TString name, TString suffix);

///
/// Function to get a TH1F from a TList
TH1F* GetTH1F(TList* l, TString name);

///
/// Function to get a TH2F from a TList
TH2F* GetTH2F(TList* l, TString name);

///
/// Function to get a Profile from a TH2F
TProfile* MakeProfile(TH2F* h, TString name, Float_t minY, Float_t maxY, Bool_t xaxis = kTRUE);

///
///Function to setup the histogram style
void MakeUpHisto(TH1* histo, TString titleX = "", TString titleY = "", Int_t marker = 20, Color_t color = kBlue + 2, Int_t lineWidth = 1);
void MakeUpHisto(TH2* histo, TString titleX = "", TString titleY = "", TString titleZ = "", Int_t marker = 20, Color_t color = kBlue + 2, Int_t lineWidth = 1);

///
///Function to setup the Pad style
void SetupPad(TString opt);

///
///Function to get a standard legend
TLegend* GetLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2);

///
///Function to get a standard legend
TPaveText* GetPave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t col, TString txt);

///
///Function add a label in a canvas, indicating a missing Plot
void AddMissingLabel(const TString histoname);

///
///Function check if the histogram exists and draw it or add a label in a canvas, indicating a missing Plot

#define CheckAndWrite(obj)          \
  if (saveHisto) {                  \
    trendFile->cd();                \
    if (obj)                        \
      obj->Write();                 \
    else {                          \
      TNamed miss(#obj, "MISSING"); \
      miss.Write();                 \
    }                               \
  }

#define CheckAndPrint(c, name) \
  if (savePng)                 \
    c->Print(Form("%s/%i%s_%s.png", plotDir.Data(), runNumber, dirsuffix.Data(), name));

#define CheckAndDraw(h, prof, opt) \
  if (h) {                         \
    h->Draw(opt);                  \
    TProfile* PProf = prof;        \
    if (PProf)                     \
      PProf->Draw("same");         \
  } else                           \
    AddMissingLabel(#h);

/********************************************************************************/
TF1* TOFsignal(Double_t rangeMin, Double_t rangeMax);
TList* FitNSigma(TList* l, TString part, TF1* f, TF1* f2, const TString name = "hTOFpidSigma", const TString suffix = "_all");

/********************************************************************************/
Int_t MakeTrendingTOFQA(char* runlist,
    Int_t year = 2010,
    TString period = "LHC10c",
    TString pass = "cpass1_pass4",
    TString nameSuffix = "_barrel",
    Bool_t isMC = kFALSE,
    Int_t trainId = 0,
    Double_t RangeFitNsigmaPIDmin = -2.0,
    Double_t RangeFitNsigmaPIDmax = 2.0,
    Double_t RangeTrksForTOFResMin = 10.,  //set lower limit to the number of tracks requested to extract the TOF resolution
    Double_t RangeTrksForTOFResMax = 100., //set upper limit to the number of tracks requested to extract the TOF resolution
    Bool_t saveHisto = kTRUE,
    Bool_t checkPIDqa = kTRUE,
    Bool_t fitSignalModel = kTRUE,
    Bool_t drawAll = kTRUE,
    Bool_t IsOnGrid = kTRUE,
    TString dirsuffix = "")
{

  if (IsOnGrid)
    TGrid::Connect("alien://");

  Int_t filesCounter = 0;
  Int_t processedCounter = 0;
  if (!runlist) {
    printf("Invalid list of runs given as input: nothing done\n");
    return 1;
  }
  Int_t runNumber = -1;
  FILE* files = fopen(runlist, "r");
  if (files == NULL) {
    ::Error("MakeTrendingTOFQA", "Error opening file %s", runlist);
    return 1;
  }
  while (fscanf(files, "%d", &runNumber) == 1) {
    TString qafilename = "alien:///alice/";
    //get QAtrain output
    if (trainId == 0) {
      if (!isMC)
        qafilename += Form("data/%i/%s/000%d/%s/", year, period.Data(), runNumber, pass.Data());
      else
        qafilename += Form("sim/%i/%s/%d/", year, period.Data(), runNumber);
    } else {
      if (!isMC)
        qafilename += Form("data/%i/%s/000%d/ESDs/%s/QA%i/", year, period.Data(), runNumber, pass.Data(), trainId);
      else
        qafilename += Form("sim/%i/%s/%d/QA%i/", year, period.Data(), runNumber, trainId);
    }
    qafilename += "QAresults" + nameSuffix + ".root";

    Printf("============== Opening QA file (%s) for run %i ==============\n", qafilename.Data(), runNumber);

    //run post-analysis
    if (MakeTrendingTOFQAv2(qafilename, runNumber, dirsuffix, isMC, checkPIDqa, fitSignalModel, RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax, RangeTrksForTOFResMin, RangeTrksForTOFResMax, "raw://", drawAll, saveHisto, kTRUE) == 0) {
      processedCounter++;
    } else
      Printf("Post analysis not run on QA output %s i.e. input #%i", qafilename.Data(), filesCounter);
    filesCounter++;
  }
  Printf(":::: Successfully Processed %i/%i runs", processedCounter, filesCounter);
  return 0;
}

//---------------------------------------------------------------------------------
Int_t MakeTrendingTOFQAv2(const TString qafilename, //full path of the QA output;
    Int_t runNumber,                                //run number
    const TString dirsuffix,                        //suffix for subdirectories
    const Bool_t isMC,                              //MC flag, to disable meaningless checks
    const Bool_t checkPIDqa,                        //set to kTRUE to check PIDqa output for TOF
    const Bool_t fitSignalModel,                    //set to kTRUE to perform fit with TOF signal model too for PID
    const Double_t RangeFitNsigmaPIDmin,            //set lower limit for fitting Nsigma_TOF_ID
    const Double_t RangeFitNsigmaPIDmax,            //set upper limit for fitting Nsigma_TOF_ID
    const Double_t RangeTrksForTOFResMin,           //set lower limit to the number of tracks requested to extract the TOF resolution
    const Double_t RangeTrksForTOFResMax,           //set upper limit to the number of tracks requested to extract the TOF resolution
    const TString ocdbStorage,                      //set the default ocdb storage
    const Bool_t drawAll,                           //enable display plots on canvas and save png
    const Bool_t saveHisto,                         //set to kTRUE to save histograms in root file
    const Bool_t savePng,                           //set to kTRUE to save histogram to png image
    const Bool_t isAutoTrend)                       //set to kTRUE for automatic trending
{
  // macro to generate tree with TOF QA trending variables
  // access qa PWGPP output files
  if (!qafilename) {
    Printf("Error - Invalid input file");
    return 1;
  }

  /*set graphic style*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatBorderSize(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(10);

  TString treePostFileName = "trending.root";
  if (!isAutoTrend) {
    treePostFileName = Form("trending_%i.root", runNumber);
    if (!dirsuffix.IsNull())
      treePostFileName.ReplaceAll("trending_", Form("trending_%s_", dirsuffix.Data()));
  }

  TFile* fin = TFile::Open(qafilename, "READ");
  if (!fin || !fin->IsOpen()) {
    Printf("ERROR: QA output not found. Exiting...\n");
    return -1;
  } else
    Printf("INFO: QA output file %s open. \n", fin->GetName());

//access list directories
#define GetIt(Var, name)            \
  TDirectoryFile* Var = nullptr;    \
  if (GetTDir(Var, fin, name) != 0) \
    return -1;

  GetIt(tofQAdir, "TOF");
  GetIt(pidQAdir, "PIDqa");
#undef GetIt

//access histograms lists
#define GetIt(Var, name)                   \
  TList* Var = nullptr;                    \
  GetList(tofQAdir, Var, name, dirsuffix); \
  if (!Var)                                \
    ::Error("MakeTrendingTOFQAv2", "%s '%s' list for QA histograms is absent or not accessible\n", #Var, name);

  GetIt(generalList, "base_noPID");
  GetIt(timeZeroList, "timeZero_noPID");
  GetIt(pidList, "pid_noPID");
  TList* pidListT0 = 0x0;
  GetList(pidQAdir, pidListT0, "PIDqa", "");
  TList* tofPidListT0 = pidListT0 ? (TList*)pidListT0->FindObject("TOF") : 0x0;
  if (!pidQAdir || !tofPidListT0) {
    printf("WARNING: PIDqa histograms not available\n");
  }
  GetIt(trdList, "trd_noPID");
  GetIt(trgList, "trigger_noPID");
#undef GetIt

  //Check on the input
  if ((!generalList) && (!timeZeroList) && (!pidList)) {
    Printf("FATAL ERROR: no QA available \n");
    return -1;
  }

  Printf(":::: Getting post-analysis info for run %i", runNumber);
  TFile* trendFile = new TFile(treePostFileName, "recreate");
  if (!trendFile)
    Printf(":::: ERROR creating output file.");
  //
  Double_t avTime = -9999., peakTime = -9999., spreadTime = -9999., peakTimeErr = -9999., spreadTimeErr = -9999., negTimeRatio = -9999.;
  Double_t avRawTime = -9999., peakRawTime = -9999., spreadRawTime = -9999., peakRawTimeErr = -9999., spreadRawTimeErr = -9999.;
  Double_t avTot = -9999., peakTot = -9999., spreadTot = -9999., peakTotErr = -9999., spreadTotErr = -9999.;
  Double_t meanResTOF = -999., spreadResTOF = -999., meanResTOFerr = -999., spreadResTOFerr = -999.;
  Double_t orphansRatio = -9999., avL = -9999., negLratio = -9999.;
  Double_t matchEffIntegratedErr = -9999., matchEffIntegrated = -9999., matchEffLinFit1Gev = -9999., matchEffLinFit1GevErr = -9999.;
  Double_t avPiDiffTime = -9999., peakPiDiffTime = -9999., spreadPiDiffTime = -9999., peakPiDiffTimeErr = -9999., spreadPiDiffTimeErr = -9999.;
  Double_t avT0A = -9999., peakT0A = -9999., spreadT0A = -9999., peakT0AErr = -9999., spreadT0AErr = -9999.;
  Double_t avT0C = -9999., peakT0C = -9999., spreadT0C = -9999., peakT0CErr = -9999., spreadT0CErr = -9999.;
  Double_t avT0AC = -9999., peakT0AC = -9999., spreadT0AC = -9999., peakT0ACErr = -9999., spreadT0ACErr = -9999.;
  Double_t avT0res = -9999., peakT0res = -9999., spreadT0res = -9999., peakT0resErr = -9999., spreadT0resErr = -9999.;
  Double_t avT0fillRes = -9999., avT0T0Res = -9999.;
  Double_t StartTime_pBestT0 = 0.0, StartTime_pBestT0Err = 0.0, StartTime_pFillT0 = 0.0, StartTime_pFillT0Err = 0.0;
  Double_t StartTime_pTOFT0 = 0.0, StartTime_pTOFT0Err = 0.0, StartTime_pT0ACT0 = 0.0, StartTime_pT0ACT0Err = 0.0, StartTime_pT0AT0 = 0.0, StartTime_pT0AT0Err = 0.0, StartTime_pT0CT0 = 0.0, StartTime_pT0CT0Err = 0.0;
  Double_t StartTime_pBestT0_Res = 0.0, StartTime_pFillT0_Res = 0.0, StartTime_pTOFT0_Res = 0.0, StartTime_pT0ACT0_Res = 0.0, StartTime_pT0AT0_Res = 0.0, StartTime_pT0CT0_Res = 0.0;
  Float_t avMulti = 0;
  Float_t fractionEventsWHits = -9999.;
  /*number of good (HW ok && efficient && !noisy) TOF channels from OCDB*/
  Double_t goodChannelsRatio = 0.0;
  Double_t goodChannelsRatioInAcc = 0.0;

  TTree* ttree = new TTree("trending", "tree of trending variables"); // TTree to store the trending values for each run
  Char_t VarType = 'F';
#define SetBranch(var) ttree->Branch(#var, &var, Form("%s/%c", #var, VarType));
  ttree->Branch("run", &runNumber, "run/I"); //run number
  SetBranch(avMulti);                        //average number of hits/event on the TOF
  SetBranch(fractionEventsWHits);            //fraction of events with hits on the TOF
  VarType = 'D';
  SetBranch(goodChannelsRatio);      //ratio of good TOF channels
  SetBranch(goodChannelsRatioInAcc); //ratio of good TOF channels in |eta|<0.8
  SetBranch(avTime);                 //mean time
  SetBranch(peakTime);               //main peak time after fit
  SetBranch(spreadTime);             //spread of main peak of time after fit
  SetBranch(peakTimeErr);            //main peak time after fit error
  SetBranch(spreadTimeErr);          //spread of main peak of time after fit error
  SetBranch(negTimeRatio);           //negative time ratio
  SetBranch(avRawTime);              //mean raw time
  SetBranch(peakRawTime);            //mean peak of RAW TIME after fit
  SetBranch(spreadRawTime);          //spread of main peak of raw time after fit
  SetBranch(peakRawTimeErr);         //main peak raw  time after fit error
  SetBranch(spreadRawTimeErr);       //spread of  raw main peak of time after fit error
  SetBranch(avTot);                  //main peak tot
  SetBranch(peakTot);                // main peak of tot after fit
  SetBranch(spreadTot);              //spread of main peak of tot after fit
  SetBranch(peakTotErr);             // main peak of tot after fit
  SetBranch(spreadTotErr);           //spread of main peak of tot after fit
  SetBranch(meanResTOF);             //mean of gaussian fit to the t-texp distribution
  SetBranch(spreadResTOF);           //sigma of gaussian fit to the t-texp distribution
  SetBranch(meanResTOFerr);          //error on the mean of gaussian fit to the t-texp distribution
  SetBranch(spreadResTOFerr);        //error on the sigma of gaussian fit to the t-texp distribution
  SetBranch(orphansRatio);           //orphans ratio
  SetBranch(avL);                    //mean track length
  SetBranch(negLratio);              //ratio of tracks with track length <350 cm
  SetBranch(matchEffIntegrated);     //matching eff integrated in pt (1-10GeV/c)
  SetBranch(matchEffIntegratedErr);  //matching eff integrated in pt (1-10GeV/c)
  SetBranch(matchEffLinFit1Gev);     //matching eff fit param
  SetBranch(matchEffLinFit1GevErr);  ////matching eff fit param error
  SetBranch(avPiDiffTime);           //mean t-texp
  SetBranch(peakPiDiffTime);         //main peak t-texp after fit
  SetBranch(spreadPiDiffTime);       //spread of main peak t-texp after fit
  SetBranch(peakPiDiffTimeErr);      //main peak t-texp after fit error
  SetBranch(spreadPiDiffTimeErr);    //spread of main peak of t-texp after fit error
  SetBranch(avT0A);                  //main peak t0A
  SetBranch(peakT0A);                // main peak of t0A after fit
  SetBranch(spreadT0A);              //spread of main peak of t0A after fit
  SetBranch(peakT0AErr);             // main peak of t0A after fit
  SetBranch(spreadT0AErr);           //spread of main peak of t0A after fit
  SetBranch(avT0C);                  //main peak t0C
  SetBranch(peakT0C);                // main peak of t0C after fit
  SetBranch(spreadT0C);              //spread of main peak of t0C after fit
  SetBranch(peakT0CErr);             // main peak of t0C after fit
  SetBranch(spreadT0CErr);           //spread of main peak of t0C after fit
  SetBranch(avT0AC);                 //main peak t0AC
  SetBranch(peakT0AC);               // main peak of t0AC after fit
  SetBranch(spreadT0AC);             //spread of main peak of t0AC after fit
  SetBranch(peakT0ACErr);            // main peak of t0AC after fit
  SetBranch(spreadT0ACErr);          //spread of main peak of t0AC after fit
  SetBranch(avT0res);                //main peak t0AC
  SetBranch(peakT0res);              // main peak of t0AC after fit
  SetBranch(spreadT0res);            //spread of main peak of t0AC after fit
  SetBranch(peakT0resErr);           // main peak of t0AC after fit
  SetBranch(spreadT0resErr);         //spread of main peak of t0AC after fit
  SetBranch(avT0fillRes);            //t0 fill res
  SetBranch(avT0T0Res);              //t0 T0 res
  SetBranch(StartTime_pBestT0);      //T0Best
  SetBranch(StartTime_pBestT0Err);
  SetBranch(StartTime_pFillT0); //T0Fill
  SetBranch(StartTime_pFillT0Err);
  SetBranch(StartTime_pTOFT0);    //T0TOF
  SetBranch(StartTime_pTOFT0Err); //T0TOF
  SetBranch(StartTime_pT0ACT0);   //T0AC
  SetBranch(StartTime_pT0ACT0Err);
  SetBranch(StartTime_pT0AT0); //T0A
  SetBranch(StartTime_pT0AT0Err);
  SetBranch(StartTime_pT0CT0); //T0C
  SetBranch(StartTime_pT0CT0Err);
  SetBranch(StartTime_pBestT0_Res); //T0Best res
  SetBranch(StartTime_pFillT0_Res); //T0Fill res
  SetBranch(StartTime_pTOFT0_Res);  //T0TOF res
  SetBranch(StartTime_pT0ACT0_Res); //T0AC res
  SetBranch(StartTime_pT0AT0_Res);  //T0A res
  SetBranch(StartTime_pT0CT0_Res);  //T0C res
#undef SetBranch

  //save quantities for trending
  goodChannelsRatio = GetGoodTOFChannelsRatio(runNumber, kFALSE, ocdbStorage, kFALSE);
  goodChannelsRatioInAcc = GetGoodTOFChannelsRatio(runNumber, kFALSE, ocdbStorage, kTRUE);

  //--------------------------------- Multiplicity ----------------------------------//
  TH1F* hMulti = GetTH1F(generalList, "hTOFmulti_all");
  TH1F* hFractionEventsWhits = new TH1F("hFractionEventsWhits", "hFractionEventsWhits;fraction of events with hits (%)", 200, 0., 100.);
  Float_t fraction = 0.0;
  if (hMulti->GetEntries() > 0.0) {
    fraction = hMulti->GetBinContent(1) / hMulti->GetEntries();
    avMulti = hMulti->GetMean();
  } else {
    fraction = 0.0;
  }
  hFractionEventsWhits->Fill(fraction * 100.);
  CheckAndWrite(hMulti);

  //--------------------------------- T0F signal ----------------------------------//
#define FitForTime(Hname, fun, av, peak, spread, peakErr, spreadErr, col, mar) \
  TH1F* Hname = GetTH1F(generalList, Form("%s_all", #Hname));                  \
  Hname->SetDirectory(0);                                                      \
  if ((Hname) && (Hname->GetEntries() > 0)) {                                  \
    av = Hname->GetMean();                                                     \
    if (!isMC) {                                                               \
      Hname->Fit(fun, "RQ0", "", 200., 250.);                                  \
      if (Hname->GetFunction(fun)) {                                           \
        peak = (Hname->GetFunction(fun))->GetParameter(1);                     \
        spread = (Hname->GetFunction(fun))->GetParameter(2);                   \
        peakErr = (Hname->GetFunction(fun))->GetParError(1);                   \
        spreadErr = (Hname->GetFunction(fun))->GetParError(2);                 \
      }                                                                        \
    } else {                                                                   \
      printf("Reminder: Raw time not available in MC simulated data.");        \
    }                                                                          \
  }                                                                            \
  MakeUpHisto(Hname, "", "matched tracks", mar, col, 1);

  FitForTime(hRawTime, "landau", avRawTime, peakRawTime, spreadRawTime, peakRawTimeErr, spreadRawTimeErr, kGreen + 2, 21);
  FitForTime(hTime, "landau", avTime, peakTime, spreadTime, peakTimeErr, spreadTimeErr, kBlue + 2, 20);
  if (hTime && hTime->GetEntries() > 0)
    negTimeRatio = (hTime->Integral(1, 3) / hTime->Integral()) * 100.;
  //
  FitForTime(hTot, "gaus", avTot, peakTot, spreadTot, peakTotErr, spreadTotErr, kViolet - 3, 8);
#undef FitForTime

  if (hTot->GetEntries() > 1)
    orphansRatio = hTot->GetBinContent(1) / hTot->GetEntries();
  //
  TPaveText* tOrphans = GetPave(0.38, 0.63, 0.88, 0.7, kViolet - 3, Form("orphans/matched = %4.2f%%", orphansRatio * 100.));

  TH1F* hL = GetTH1F(generalList, "hMatchedL_all");
  if (hL->GetEntries() > 0) {
    avL = hL->GetMean();
    negLratio = hL->Integral(1, 750) / hL->GetEntries();
  }
  MakeUpHisto(hL, "", "matched tracks", 1, kBlue + 2, 1);
  TPaveText* tLength = GetPave(0.15, 0.83, 0.65, 0.87, kOrange - 3, Form("trk with L<350cm /matched = %4.2f%%", negLratio * 100.));

  CheckAndWrite(hRawTime);
  CheckAndWrite(hTime);
  CheckAndWrite(hTot);
  CheckAndWrite(hL);

  //--------------------------------- residuals -------------------------------------//
  TH2F* hDxPos4profile = GetTH2F(generalList, "hMatchedDxVsPt_all");
  TH2F* hTOFmatchedDzVsStrip = GetTH2F(generalList, "hMatchedDzVsStrip_all");
  CheckAndWrite(hDxPos4profile);
  CheckAndWrite(hTOFmatchedDzVsStrip);

  //--------------------------------- matching eff ----------------------------------//
  //matching as function of pT
  //matching as function of eta and phi out
  const Double_t maxPtEff2Fit = 10.0;
  const Double_t minPtEff2Fit = 1.0;

  TH1F* hMatchingVsPt = NULL;
  TH2F* hMatchingVsEtaPhiOut = NULL;
  TH1F* hMatchingVsEta = NULL;
  TH1F* hMatchingVsPhiOut = NULL;
  TH1F* hMatchingVsPhi = NULL;

  TH1F* hDenom = GetTH1F(generalList, "hPrimaryPt_all");
  if (hDenom) {
    hDenom->Sumw2();
    hMatchingVsPt = (TH1F*)(GetTH1F(generalList, "hMatchedPt_all"))->Clone("hMatchingVsPt");
    //set underflow bin to the matching efficiency integrated in 1-10 GeV/c
    Int_t imin = hMatchingVsPt->GetXaxis()->FindBin(minPtEff2Fit);
    Int_t imax = hMatchingVsPt->GetXaxis()->FindBin(maxPtEff2Fit);
    Double_t numErr, denErr;
    Double_t numN = hMatchingVsPt->IntegralAndError(imin, imax, numErr);
    Double_t denN = hDenom->IntegralAndError(imin, imax, denErr);
    std::pair<Double_t, Double_t> IntEff = ComputeEff(numN, denN, numErr, denErr);

    matchEffIntegrated = IntEff.first;
    matchEffIntegratedErr = IntEff.second;
    //get efficiency
    hMatchingVsPt->Sumw2();
    hMatchingVsPt->Divide(hMatchingVsPt, hDenom, 1., 1., "B");
    hMatchingVsPt->SetTitle("TOF matching efficiency as function of transverse momentum");
    hMatchingVsPt->GetYaxis()->SetRangeUser(0, 1.2);
  }

  if (hMatchingVsPt->GetEntries() > 0) {
    hMatchingVsPt->Fit("pol0", "", "", minPtEff2Fit, maxPtEff2Fit);
    hMatchingVsPt->Draw();
    if (hMatchingVsPt->GetFunction("pol0")) {
      matchEffLinFit1Gev = (hMatchingVsPt->GetFunction("pol0"))->GetParameter(0);
      matchEffLinFit1GevErr = (hMatchingVsPt->GetFunction("pol0"))->GetParError(0);
    }
  } else {
    printf("WARNING: matching efficiency plot has 0 entries. Skipped!\n");
  }
  MakeUpHisto(hMatchingVsPt, "#it{p}_{T} (GeV/#it{c})", "matching efficiency", 1, kBlue + 2, 2);

  TH2F* hDenom2D = GetTH2F(generalList, "hPrimaryEtaVsOutPhi_all");
  if (!hDenom2D) {
    //matching as function of eta
    hDenom = GetTH1F(generalList, "hPrimaryEta_all");
    if (hDenom) {
      hDenom->Sumw2();
      hMatchingVsEta = (TH1F*)(GetTH1F(generalList, "hMatchedEta_all"))->Clone("hMatchingVsEta");
      hMatchingVsEta->Sumw2();
      hMatchingVsEta->Divide(hMatchingVsEta, hDenom, 1., 1., "B");
    }
  } else {
    hMatchingVsEtaPhiOut = (TH2F*)(GetTH2F(generalList, "hMatchedEtaVsOutPhi_all"))->Clone("hMatchingVsEtaPhiOut");
    //
    Compute2Deff(hMatchingVsEtaPhiOut, hDenom2D, hMatchingVsEta, "hMatchingVsEta", kFALSE);
    Compute2Deff(hMatchingVsEtaPhiOut, hDenom2D, hMatchingVsPhiOut, "hMatchingVsPhiOut", kTRUE);
    //
    hMatchingVsEtaPhiOut->Sumw2();
    hMatchingVsEtaPhiOut->Divide(hMatchingVsEtaPhiOut, hDenom2D, 1., 1., "B");
    hMatchingVsEtaPhiOut->GetZaxis()->SetRangeUser(0., 1.);
    hMatchingVsEtaPhiOut->SetTitle("TOF matching efficiency as function of #eta and #phi_{TPC,out}");
    MakeUpHisto(hMatchingVsEtaPhiOut, "#phi_{TPC,out} (deg)", "#eta", "matching efficiency", 1, kBlue + 2, 2);
  }

  if (hMatchingVsEta) {
    hMatchingVsEta->SetTitle("TOF matching efficiency as function of #eta");
    hMatchingVsEta->GetYaxis()->SetRangeUser(0, 1.2);
    MakeUpHisto(hMatchingVsEta, "#eta", "matching efficiency", 1, kBlue + 2, 2);
  }

  if (hMatchingVsPhiOut) {
    hMatchingVsPhiOut->SetTitle("TOF matching efficiency as function of #phi_{TPC,out}");
    hMatchingVsPhiOut->GetYaxis()->SetRangeUser(0, 1.2);
    MakeUpHisto(hMatchingVsPhiOut, "#phi_{TPC,out}", "matching efficiency", 1, kBlue + 2, 2);
  }

  //matching as function of phi
  hDenom = GetTH1F(generalList, "hPrimaryPhi_all");
  if (hDenom) {
    hDenom->Sumw2();
    hMatchingVsPhi = (TH1F*)(GetTH1F(generalList, "hMatchedPhi_all"))->Clone("hMatchingVsPhi");
    hMatchingVsPhi->Sumw2();
    hMatchingVsPhi->Divide(hMatchingVsPhi, hDenom, 1., 1., "B");
    hMatchingVsPhi->SetTitle("TOF matching efficiency as function of phi");
    hMatchingVsPhi->GetYaxis()->SetRangeUser(0, 1.2);
  }
  MakeUpHisto(hMatchingVsPhi, "#phi (deg)", "matching efficiency", 1, kBlue + 2, 2);

  CheckAndWrite(hMatchingVsPt);
  CheckAndWrite(hMatchingVsEta);
  CheckAndWrite(hMatchingVsPhi);
  CheckAndWrite(hMatchingVsPhiOut);
  CheckAndWrite(hMatchingVsEtaPhiOut);

  //--------------------------------- t-texp ----------------------------------//
  TH2F* hBetaP = GetTH2F(pidList, "hMatchedBetaVsP_all");
  if (hBetaP)
    hBetaP->GetYaxis()->SetRangeUser(0., 1.2);

  TH1F* hMass = GetTH1F(pidList, "hMatchedMass_all");
  MakeUpHisto(hMass, "", "tracks", 1, kBlue + 2, 1);
  // hMass->SetFillColor(kAzure+10);
  // hMass->SetFillStyle(1001);
  hMass->Rebin(2);

  //pions
  TH1F* hPionDiff = GetTH1F(pidList, "hExpTimePi_all");
  if ((hPionDiff) && (hPionDiff->GetEntries() > 0)) {
    avPiDiffTime = hPionDiff->GetMean();
    hPionDiff->Fit("gaus", "", "", -1000., 500.);
    if (hPionDiff->GetFunction("gaus")) {
      peakPiDiffTime = (hPionDiff->GetFunction("gaus"))->GetParameter(1);
      spreadPiDiffTime = (hPionDiff->GetFunction("gaus"))->GetParameter(2);
      peakPiDiffTimeErr = (hPionDiff->GetFunction("gaus"))->GetParError(1);
      spreadPiDiffTimeErr = (hPionDiff->GetFunction("gaus"))->GetParError(2);
    }
  }
  TCanvas* cTimeCalib = new TCanvas("cTimeCalib", "cTimeCalib", 800, 600);
  gStyle->SetOptStat(10);
  gStyle->SetOptFit();
  cTimeCalib->cd();
  gPad->SetLogy();
  MakeUpHisto(hPionDiff, "", "", 1, kBlue + 1, 1);
  hPionDiff->Draw();
  TString plotDir(".");
  CheckAndPrint(cTimeCalib, "TOFtime");

  //Retrieve plots for t-texp-t0
  TH1F* hDiffTimeT0fillPion = GetTH1F(pidList, "hExpTimePiFillSub_all");
  TH2F* hDiffTimeT0TOFPion1GeV = GetTH2F(pidList, "hExpTimePiT0Sub1GeV_all");

  //Pion
  TH2F* hDiffTimePi = GetTH2F(pidList, "hExpTimePiVsP_all");
  hDiffTimePi->SetTitle("PIONS t-t_{exp,#pi} (from tracking) vs. P");

  //Kaon
  TH2F* hDiffTimeKa = GetTH2F(pidList, "hExpTimeKaVsP_all");
  hDiffTimeKa->SetTitle("KAONS t-t_{exp,K} (from tracking) vs. P");

  //Protons
  TH2F* hDiffTimePro = GetTH2F(pidList, "hExpTimeProVsP_all");
  hDiffTimePro->SetTitle("PROTONS t-t_{exp,p} (from tracking) vs. P");

  CheckAndWrite(hBetaP);
  CheckAndWrite(hMass);
  CheckAndWrite(hPionDiff);
  CheckAndWrite(hDiffTimeT0fillPion);
  CheckAndWrite(hDiffTimeT0TOFPion1GeV);
  CheckAndWrite(hDiffTimePi);
  CheckAndWrite(hDiffTimeKa);
  CheckAndWrite(hDiffTimePro);

  //---------------------------------TOF resolution plot ----------------------------------//
  Int_t min = hDiffTimeT0TOFPion1GeV->GetXaxis()->FindBin(RangeTrksForTOFResMin);
  Int_t max = hDiffTimeT0TOFPion1GeV->GetXaxis()->FindBin(RangeTrksForTOFResMax);
  TH1D* hResTOF = (TH1D*)hDiffTimeT0TOFPion1GeV->ProjectionY("hTOFres1GeV", min, max);
  hResTOF->Fit("gaus", "", "", -200., 100.);
  if (hResTOF->GetFunction("gaus")) {
    meanResTOF = (hResTOF->GetFunction("gaus"))->GetParameter(1);
    spreadResTOF = (hResTOF->GetFunction("gaus"))->GetParameter(2);
    meanResTOFerr = (hResTOF->GetFunction("gaus"))->GetParError(1);
    spreadResTOFerr = (hResTOF->GetFunction("gaus"))->GetParError(2);
  }

  TCanvas* cTOFresolution = new TCanvas("cTOFresolution", "TOFresolution", 1200, 500);
  gStyle->SetOptStat(10);
  gStyle->SetOptFit();
  cTOFresolution->Divide(2, 1);
  cTOFresolution->cd(1);
  gPad->SetLogz();
  hDiffTimeT0TOFPion1GeV->GetXaxis()->SetRangeUser(0., 50.);
  hDiffTimeT0TOFPion1GeV->Draw("colz");
  cTOFresolution->cd(2);
  hResTOF->GetXaxis()->SetRangeUser(-1000, 1000);
  MakeUpHisto(hResTOF, "", "", 1, kBlue + 1, 1);
  hResTOF->Draw();

  CheckAndPrint(cTOFresolution, "TOFresolution");
  CheckAndWrite(hResTOF);

  //--------------------------------- T0 vs multiplicity plots ----------------------------------//
  TList* ListOfOutput_T0 = new TList(); ///List of output for all plots related to T0

  TH2F* hT0TOFvsNtracks = GetTH2F(timeZeroList, "hT0TOFvsNtrk");
  TProfile* hT0TOFProfile = MakeProfile(hT0TOFvsNtracks, "hT0TOFProfile", -50., 50.);
  if (hT0TOFvsNtracks) { //Add result to output list
    ListOfOutput_T0->Add(hT0TOFvsNtracks);
    ListOfOutput_T0->Add(hT0TOFProfile);
  }

  TH2F* hT0ACvsNtracks = GetTH2F(timeZeroList, "hT0ACvsNtrk");
  TProfile* hT0ACProfile = MakeProfile(hT0ACvsNtracks, "hT0ACProfile", -50., 50.);
  if (hT0ACvsNtracks) { //Add result to output list
    ListOfOutput_T0->Add(hT0ACvsNtracks);
    ListOfOutput_T0->Add(hT0ACProfile);
  }

  TH2F* hT0AvsNtracks = GetTH2F(timeZeroList, "hT0AvsNtrk");
  TProfile* hT0AProfile = MakeProfile(hT0AvsNtracks, "hT0AProfile", -50., 50.);
  if (hT0AvsNtracks) { //Add result to output list
    ListOfOutput_T0->Add(hT0AvsNtracks);
    ListOfOutput_T0->Add(hT0AProfile);
  }

  TH2F* hT0CvsNtracks = GetTH2F(timeZeroList, "hT0CvsNtrk");
  TProfile* hT0CProfile = MakeProfile(hT0CvsNtracks, "hT0CProfile", -50., 50.);
  if (hT0CvsNtracks) { //Add result to output list
    ListOfOutput_T0->Add(hT0CvsNtracks);
    ListOfOutput_T0->Add(hT0CProfile);
  }

  TH2F* hStartTime = GetTH2F(timeZeroList, "hStartTime");
  TProfile* hStartTimeProfile = MakeProfile(hStartTime, "hStartTimeProfile", -600., 600., kFALSE);
  if (hStartTime) { //Add result to output list
    ListOfOutput_T0->Add(hStartTime);
    ListOfOutput_T0->Add(hStartTimeProfile);

    //Set Start Time information
#define GetT0Info(Label, Tz, TzErr)                                         \
  Tzbin = hStartTime->GetYaxis()->FindBin(Label);                           \
  if (Tzbin <= 0 || Tzbin > hStartTime->GetNbinsY())                        \
    ::Error("MakeTrendingTOFQAv2", "cannot find start time bin %s", Label); \
  Tz = hStartTimeProfile->GetBinContent(Tzbin);                             \
  TzErr = hStartTimeProfile->GetBinError(Tzbin);
    Int_t Tzbin = -1;
    GetT0Info("best_t0", StartTime_pBestT0, StartTime_pBestT0Err);
    GetT0Info("fill_t0", StartTime_pFillT0, StartTime_pFillT0Err);
    GetT0Info("tof_t0", StartTime_pTOFT0, StartTime_pTOFT0Err);
    GetT0Info("T0AC", StartTime_pT0ACT0, StartTime_pT0ACT0Err);
    GetT0Info("T0A", StartTime_pT0AT0, StartTime_pT0AT0Err);
    GetT0Info("T0C", StartTime_pT0CT0, StartTime_pT0CT0Err);
#undef GetT0Info
  }

  TH2F* hStartTimeRes = GetTH2F(timeZeroList, "hStartTimeRes");
  TProfile* hStartTimeResProfile = MakeProfile(hStartTimeRes, "hStartTimeResProfile", 0, 299., kFALSE);
  if (hStartTimeRes) { //Add result to output list
    ListOfOutput_T0->Add(hStartTimeRes);
    ListOfOutput_T0->Add(hStartTimeResProfile);

//Set Start Time Resolution information
#define GetT0Info(Label, Tz)                                                \
  Tzbin = hStartTimeRes->GetYaxis()->FindBin(Label);                        \
  if (Tzbin <= 0 || Tzbin > hStartTimeRes->GetNbinsY())                     \
    ::Error("MakeTrendingTOFQAv2", "cannot find start time bin %s", Label); \
  Tz = hStartTimeResProfile->GetBinContent(Tzbin);
    Int_t Tzbin = -1;
    GetT0Info("best_t0", StartTime_pBestT0_Res);
    GetT0Info("fill_t0", StartTime_pFillT0_Res);
    GetT0Info("tof_t0", StartTime_pTOFT0_Res);
    GetT0Info("T0AC", StartTime_pT0ACT0_Res);
    GetT0Info("T0A", StartTime_pT0AT0_Res);
    GetT0Info("T0C", StartTime_pT0CT0_Res);
#undef GetT0Info
  }

  // cout << "StartTimeBest  " << StartTime_pBestT0 << " ps" << endl;
  // cout << "StartTimeBestErr  " << StartTime_pBestT0Err << " ps" << endl;
  // cout << "Res StartTimeBest  " << StartTime_pBestT0_Res << " ps" << endl;
  // cout << "Res StartTimeBestErr  " << StartTime_pBestT0_ResErr << " ps" << endl;
  // cout << "StartTimeFill  " << StartTime_pFillT0 << " ps" << endl;
  // cout << "StartTimeFillErr  " << StartTime_pFillT0Err << " ps" << endl;
  // cout << "Res StartTimeFill  " << StartTime_pFillT0_Res << " ps" << endl;
  // cout << "Res StartTimeFillErr  " << StartTime_pFillT0_ResErr << " ps" << endl;
  // cout << "StartTimeTOF  " << StartTime_pTOFT0 << " ps" << endl;
  // cout << "Res StartTimeTOF  " << StartTime_pTOFT0_Res << " ps" << endl;
  // cout << "StartTimeT0AC  " << StartTime_pT0ACT0 << " ps" << endl;
  // cout << "Res StartTimeT0AC  " << StartTime_pT0ACT0_Res << " ps" << endl;
  // cout << "Res StartTimeT0ACErr  " << StartTime_pT0ACT0_ResErr << " ps" << endl;
  // cout << "StartTimeT0A  " << StartTime_pT0AT0 << " ps" << endl;
  // cout << "Res StartTimeT0A  " << StartTime_pT0AT0_Res << " ps" << endl;
  // cout << "Res StartTimeT0AErr  " << StartTime_pT0AT0_ResErr << " ps" << endl;
  // cout << "StartTimeT0C  " << StartTime_pT0CT0 << " ps" << endl;
  // cout << "Res StartTimeT0C  " << StartTime_pT0CT0_Res << " ps" << endl;

  //Drawing section
  TCanvas* cT0vsMultiplicity = new TCanvas("cT0vsMultiplicity", "T0TOF,T0C,T0A,TOAC vs N_TOF,", 1200, 800);
  cT0vsMultiplicity->Divide(2, 2);

  cT0vsMultiplicity->cd(1);
  SetupPad("SetLogx SetLogz SetGridx SetGridy");
  CheckAndDraw(hT0TOFvsNtracks, hT0TOFProfile, "COLZ");

  cT0vsMultiplicity->cd(2);
  SetupPad("SetLogx SetLogz SetGridx SetGridy");
  CheckAndDraw(hT0ACvsNtracks, hT0ACProfile, "COLZ");

  cT0vsMultiplicity->cd(3);
  SetupPad("SetLogx SetLogz SetGridx SetGridy");
  CheckAndDraw(hT0AvsNtracks, hT0AProfile, "COLZ");

  cT0vsMultiplicity->cd(4);
  SetupPad("SetLogx SetLogz SetGridx SetGridy");
  CheckAndDraw(hT0CvsNtracks, hT0CProfile, "COLZ");

  TCanvas* cStartTime = new TCanvas("cStartTime", "start time with different methods", 1200, 800);
  SetupPad("SetLogz SetGridx SetGridy");
  CheckAndDraw(hStartTime, nullptr, "COLZ");

  TCanvas* cStartTimeRes = new TCanvas("cStartTimeRes", "Resolution of start time methods", 1200, 800);
  SetupPad("SetLogz SetGridx SetGridy");
  CheckAndDraw(hStartTimeRes, nullptr, "COLZ");

  CheckAndPrint(cT0vsMultiplicity, "T0vsMultiplicity");
  CheckAndPrint(cStartTime, "StartTimeMethods");
  CheckAndPrint(cStartTimeRes, "StartTimeResolution");

  if (saveHisto) {
    trendFile->cd();
    ListOfOutput_T0->Write();
    delete ListOfOutput_T0;
  }

  TLine l11(0., 0., 5., 0.);
  l11.SetLineWidth(1);
  TLine l12(0., 1., 5., 1.);
  l12.SetLineWidth(1);
#define DrawPidPerformance(l)                          \
  SetupPad("SetLogx SetLogz SetGridx SetGridy");       \
  ((TH1*)l->At(0))->GetXaxis()->SetRangeUser(0.2, 5.); \
  l->At(0)->Draw("same COLZ");                         \
  ((TH1*)l->At(1))->Draw("same");                      \
  ((TH1*)l->At(2))->Draw("same");                      \
  if (fSignalModel) {                                  \
    ((TH1*)l->At(3))->Draw("same");                    \
    ((TH1*)l->At(4))->Draw("same");                    \
  }                                                    \
  l11.DrawClone("same");                               \
  l12.DrawClone("same");                               \
  if (saveHisto) {                                     \
    CheckAndWrite(((TH1*)l->At(0)));                   \
    CheckAndWrite(((TH1*)l->At(1)));                   \
    CheckAndWrite(((TH1*)l->At(2)));                   \
  }

  //--------------------------------- NSigma PID from TOF QA ----------------------------------//
  TF1* f = new TF1("f", "gaus", RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax); //Normal gaussian
  const Float_t ModelRangeFitNsigmaPIDmin = 2 * RangeFitNsigmaPIDmin;        //-1.0
  const Float_t ModelRangeFitNsigmaPIDmax = 2 * RangeFitNsigmaPIDmax;        //+2.0
  TF1* fSignalModel = nullptr;
  if (fitSignalModel)
    fSignalModel = TOFsignal(ModelRangeFitNsigmaPIDmin, ModelRangeFitNsigmaPIDmax); //Gaussian + tail

  //----- PION ------//
  TList* lparPi = FitNSigma(pidList, "Pi", f, fSignalModel);

  //----- KAON ------//
  TList* lparKa = FitNSigma(pidList, "Ka", f, fSignalModel);

  //----- PROTON ------//
  TList* lparPro = FitNSigma(pidList, "Pro", f, fSignalModel);

  TCanvas* cPidPerformance3 = new TCanvas("cPidPerformance3", "TOF PID performance - N_{#sigma}^{TOF}", 1200, 500);

  cPidPerformance3->Divide(3, 1);
  cPidPerformance3->cd(1);
  DrawPidPerformance(lparPi);

  cPidPerformance3->cd(2);
  DrawPidPerformance(lparKa);

  cPidPerformance3->cd(3);
  DrawPidPerformance(lparPro);

  cPidPerformance3->cd(1);
  TLegend* pidLeg = GetLegend(0.15, 0.76, 0.88, 0.88);
  pidLeg->SetNColumns(2);
  pidLeg->AddEntry(lparPi->At(1), "Mean", "lp");
  pidLeg->AddEntry(lparPi->At(2), Form("#sigma, Gaus fit (%2.1f,%2.1f)", RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax), "lp");
  if (fSignalModel) {
    pidLeg->AddEntry(lparPi->At(3), "Mean", "lp");
    pidLeg->AddEntry(lparPi->At(4), Form("#sigma, Gaus+Tail fit (%2.1f,%2.1f)", ModelRangeFitNsigmaPIDmin, ModelRangeFitNsigmaPIDmax), "lp");
  }
  pidLeg->Draw("same");

  CheckAndPrint(cPidPerformance3, "PID_sigmas");

  //--------------------------------- NSigma PID from PIDqa ----------------------------------//
  TH2F* hSigmaPiT0 = 0x0;
  TH2F* hSigmaKaT0 = 0x0;
  TH2F* hSigmaProT0 = 0x0;

  if (checkPIDqa && pidQAdir) {

    //------- PIONS ------//
    TList* lparPiT0 = FitNSigma(tofPidListT0, "pion", f, fSignalModel, "hNsigmaP_TOF_", "");

    //------- KAONS ------//
    TList* lparKaT0 = FitNSigma(tofPidListT0, "kaon", f, fSignalModel, "hNsigmaP_TOF_", "");

    //------- PROTONS ------//
    TList* lparProT0 = FitNSigma(tofPidListT0, "proton", f, fSignalModel, "hNsigmaP_TOF_", "");

    //Show in canvas
    TCanvas* cPidPerformance3T0 = new TCanvas("cPidPerformance3T0", "PID performance from PIDqa - N_{#sigma}^{TOF} with StartTime", 1200, 500);
    cPidPerformance3T0->Divide(3, 1);
    cPidPerformance3T0->cd(1);
    DrawPidPerformance(lparPiT0);

    cPidPerformance3T0->cd(2);
    DrawPidPerformance(lparKaT0);

    cPidPerformance3T0->cd(3);
    DrawPidPerformance(lparProT0);

    cPidPerformance3T0->cd(1);
    pidLeg->Draw("same");

    CheckAndPrint(cPidPerformance3T0, "PID_sigmasStartTime");
  }
#undef DrawPidPerformance

  //--------------------------------- T0 detector ----------------------------------//
#define FitWithGaus(Hname, av, peak, spread, peakErr, spreadErr, col) \
  TH1F* Hname = GetTH1F(timeZeroList, #Hname);                        \
  Hname->SetDirectory(0);                                             \
  if ((Hname) && (Hname->GetEntries() > 0)) {                         \
    av = Hname->GetMean();                                            \
    Hname->Fit("gaus", "RQ0", "", -1000., 1000.);                     \
    if (Hname->GetFunction("gaus")) {                                 \
      peak = (Hname->GetFunction("gaus"))->GetParameter(1);           \
      spread = (Hname->GetFunction("gaus"))->GetParameter(2);         \
      peakErr = (Hname->GetFunction("gaus"))->GetParError(1);         \
      spreadErr = (Hname->GetFunction("gaus"))->GetParError(2);       \
    }                                                                 \
  }                                                                   \
  MakeUpHisto(Hname, "", "events", 8, col, 2);                        \
  Hname->Rebin(2);

  FitWithGaus(hT0A, avT0A, peakT0A, spreadT0A, peakT0AErr, spreadT0AErr, kBlue);
  FitWithGaus(hT0C, avT0C, peakT0C, spreadT0C, peakT0CErr, spreadT0CErr, kGreen + 1);
  FitWithGaus(hT0AC, avT0AC, peakT0AC, spreadT0AC, peakT0ACErr, spreadT0ACErr, kRed + 1);

  TLegend* lT0 = GetLegend(0.7125881, 0.6052519, 0.979435, 0.7408306);
  lT0->SetTextSize(0.041);
  lT0->AddEntry(hT0AC, "T0 A&C", "L");
  lT0->AddEntry(hT0A, "T0 A", "L");
  lT0->AddEntry(hT0C, "T0 C", "L");

  FitWithGaus(hT0DetRes, avT0res, peakT0res, spreadT0res, peakT0resErr, spreadT0resErr, kMagenta + 1);

  TH1F* hT0fillRes = GetTH1F(timeZeroList, "hT0fillRes");
  if ((hT0fillRes) && (hT0fillRes->GetEntries() > 0)) {
    avT0fillRes = hT0fillRes->GetMean();
  } else
    avT0fillRes = StartTime_pFillT0_Res;
  TH1F* hT0T0Res = GetTH1F(timeZeroList, "hT0T0Res");
  if ((hT0T0Res) && (hT0T0Res->GetEntries() > 0)) {
    avT0T0Res = hT0T0Res->GetMean();
  } else
    avT0T0Res = StartTime_pT0ACT0_Res;

  if (saveHisto) {
    CheckAndWrite(hT0AC);
    CheckAndWrite(hT0A);
    CheckAndWrite(hT0C);
    CheckAndWrite(hT0DetRes);
  }

  //close input file
  fin->Close();

  if (drawAll) {
    TCanvas* cTrackProperties = new TCanvas("cTrackProperties", "summary of matched tracks properties", 1200, 500);
    cTrackProperties->Divide(3, 1);
    cTrackProperties->cd(1);
    gPad->SetLogy();
    hTime->Draw("");
    hRawTime->Draw("same");
    cTrackProperties->cd(2);
    gPad->SetLogy();
    hTot->Draw("");
    tOrphans->Draw();
    cTrackProperties->cd(3);
    gPad->SetLogy();
    hL->Draw("");
    tLength->Draw();
    //
    TCanvas* cPidPerformance = new TCanvas("cPidPerformance", "summary of pid performance", 900, 500);
    cPidPerformance->Divide(2, 1);
    cPidPerformance->cd(1);
    gPad->SetLogz();
    hBetaP->Draw("colz");
    cPidPerformance->cd(2);
    gPad->SetLogy();
    hMass->Draw("HIST");
    //
    TCanvas* cT0detector = new TCanvas("cT0detector", "T0 detector", 800, 600);
    cT0detector->Divide(2, 2);
    cT0detector->cd(1);
    gPad->SetGridx();
    hT0A->GetXaxis()->SetRangeUser(-1000, 1000);
    hT0A->Draw();
    cT0detector->cd(2);
    gPad->SetGridx();
    hT0C->GetXaxis()->SetRangeUser(-1000, 1000);
    hT0C->Draw();
    cT0detector->cd(3);
    gPad->SetGridx();
    hT0AC->GetXaxis()->SetRangeUser(-1000, 1000);
    hT0AC->Draw();
    hT0AC->SetTitle("timeZero measured by T0 detector");
    cT0detector->cd(4);
    hT0DetRes->Draw();
    //
    CheckAndPrint(cTrackProperties, "TrackProperties");
    CheckAndPrint(cPidPerformance, "PID");
    CheckAndPrint(cT0detector, "T0Detector");
  }

  TCanvas* cProfile = new TCanvas("cProfile", "cProfile", 50, 50, 750, 550);
  gPad->SetLogz();
  hTOFmatchedDzVsStrip->Draw("colz");
  TProfile* hDzProfile = MakeProfile(hTOFmatchedDzVsStrip, "hDzProfile", -3, 3);
  hDzProfile->Draw("same");

  TCanvas* cMatchingPerformance = new TCanvas("cMatchingPerformance", "summary of matching performance", 1200, 500);
  cMatchingPerformance->Divide(3, 1);
  cMatchingPerformance->cd(1);
  hMatchingVsPt->Draw();
  cMatchingPerformance->cd(2);
  hMatchingVsEta->Draw();
  cMatchingPerformance->cd(3);
  hMatchingVsPhi->Draw();

  TCanvas* cPidPerformance2 = new TCanvas("cPidPerformance2", "summary of pid performance - expected times", 1200, 500);
  cPidPerformance2->Divide(3, 1);
  cPidPerformance2->cd(1);
  gPad->SetLogz();
  hDiffTimePi->Draw("colz");
  cPidPerformance2->cd(2);
  gPad->SetLogz();
  hDiffTimeKa->Draw("colz");
  cPidPerformance2->cd(3);
  gPad->SetLogz();
  hDiffTimePro->Draw("colz");

  if (savePng) {
    cMatchingPerformance->Print(Form("%s/%i%s_MatchingEfficiency.png", plotDir.Data(), runNumber, dirsuffix.Data()));
    cProfile->Print(Form("%s/%i%s_ProfileDZvsStripNumber.png", plotDir.Data(), runNumber, dirsuffix.Data()));
    cPidPerformance2->Print(Form("%s/%i%s_PID_ExpTimes.png", plotDir.Data(), runNumber, dirsuffix.Data()));
  }
  //Fill tree, save to file histos and canvases and delete list of canvases
  ttree->Fill();
  printf("============== Saving trending quantities in tree for run %i ==============\n", runNumber);
  trendFile->cd();
  ttree->Write();
  trendFile->Close();
  return 0;
}

//----------------------------------------------------------
Double_t GetGoodTOFChannelsRatio(Int_t run, Bool_t saveMap, TString OCDBstorage, Bool_t inEta08)
{
  /*
   *    It retrieves from OCDB the number of good (= efficient && not noisy && HW ok) TOF channels.
   *    Optionally is saves the channel map
   */
  if (run <= 0) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): invalid run number. Please set a run number.\n");
    return 0.0;
  }

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(OCDBstorage.Data());
  cdb->SetRun(run);

  AliCDBEntry* cdbe = cdb->Get("TOF/Calib/Status");
  if (!cdbe) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): OCDB entry not available. Please, try again.\n");
    return 0.0;
  }

  AliTOFChannelOnlineStatusArray* array = (AliTOFChannelOnlineStatusArray*)cdbe->GetObject();
  TH2F* hOkMap = new TH2F("hOkMap", "Ok map (!noisy & !problematic & efficient);sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2F* hOkMapInAcceptance = new TH2F("hOkMapInAcceptance", Form("Good channels in |#eta|<0.8 - run %i;sector;strip", run), 72, 0., 18., 91, 0., 91.);
  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibHisto();
  AliTOFcalib calib;
  calib.Init();
  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  for (Int_t i = 0; i < array->GetSize(); i++) {
    sector = calibHisto.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calibHisto.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calibHisto.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    if (!(array->GetNoiseStatus(i) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) && (calib.IsChannelEnabled(i, kTRUE, kTRUE))) {
      hOkMap->Fill(hitmapx, hitmapy);
      /* 3 strips / side excluded from norm. factor to consider |eta|<0.8 */
      if (sectorStrip > 2 && sectorStrip < 89)
        hOkMapInAcceptance->Fill(hitmapx, hitmapy);
    }
  }
  Int_t nOk = (Int_t)hOkMap->GetEntries();
  Int_t nOkInAcc = (Int_t)hOkMapInAcceptance->GetEntries();
  Double_t ratioOk = nOk / 152928.;
  /* 96 channels * 6 strips * 18 modules excluded from norm. factor to consider |eta|<0.8 */
  Double_t ratioOkInAcc = nOkInAcc / (152928. - 6. * 18. * 96.);
  if (saveMap) {
    hOkMap->SaveAs(Form("run%i_OKChannelsMap.root", run));
    hOkMapInAcceptance->SaveAs(Form("run%i_OKChannelsInAccMap.root", run));
  }
  cout << "###    Run " << run << ": TOF channels ok = " << nOk << "/ total 152928 channels = " << ratioOk * 100. << "% of whole TOF" << endl;
  cout << "###    Run " << run << ": TOF channels in acc. ok = " << nOkInAcc << "/ total 142560 channels = " << ratioOkInAcc * 100. << "% of whole TOF" << endl;
  if (inEta08)
    return ratioOkInAcc;
  else
    return ratioOk;
}

//----------------------------------------------------------
void MakeUpHisto(TH1* histo, TString titleX, TString titleY, Int_t marker, Color_t color, Int_t lineWidth)
{
  if (!histo)
    return;
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(0.7);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetLineWidth(lineWidth);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0);
  if (!titleX.IsNull())
    histo->GetXaxis()->SetTitle(titleX.Data());
  if (!titleY.IsNull())
    histo->GetYaxis()->SetTitle(titleY.Data());
  histo->GetYaxis()->SetTitleOffset(1.35);
  histo->GetXaxis()->SetLabelSize(0.03);
  return;
}

//----------------------------------------------------------
void MakeUpHisto(TH2* histo, TString titleX, TString titleY, TString titleZ, Int_t marker, Color_t color, Int_t lineWidth)
{
  if (!histo)
    return;
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(0.7);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetLineWidth(lineWidth);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0);
  if (!titleX.IsNull())
    histo->GetXaxis()->SetTitle(titleX.Data());
  if (!titleY.IsNull())
    histo->GetYaxis()->SetTitle(titleY.Data());
  if (!titleZ.IsNull())
    histo->GetZaxis()->SetTitle(titleZ.Data());
  histo->GetYaxis()->SetTitleOffset(1.35);
  histo->GetXaxis()->SetLabelSize(0.03);
  return;
}

//----------------------------------------------------------
void AddMissingLabel(const TString histoname)
{
  TPaveText missing(0.3, 0.5, 0.7, 0.8, "NDC NB"); //For canvases missing plots
  missing.SetFillColor(0);
  missing.SetFillStyle(0);
  missing.SetLineColor(0);
  missing.SetLineWidth(0);
  missing.SetBorderSize(0);
  missing.AddText(Form("Plot%s%s Missing", histoname.IsNull() ? "" : " ", histoname.Data()));
  missing.Draw();
}

//----------------------------------------------------------
void Compute2Deff(TH2F* num, TH2F* den, TH1F*& h, TString name, Bool_t x)
{
  h = (TH1F*)(x ? num->ProjectionX(name) : num->ProjectionY(name));
  h->Sumw2();
  TH1F* hden = (TH1F*)(x ? den->ProjectionX("hDenominatorX") : den->ProjectionY("hDenominatorY"));
  h->Divide(h, hden, 1, 1, "B");
  delete hden;
}

//----------------------------------------------------------
std::pair<Double_t, Double_t> ComputeEff(Double_t num, Double_t den, Double_t numE, Double_t denE)
{
  std::pair<Double_t, Double_t> eff(-111, 0);
  if (den > 0) {
    //
    const Double_t ratio = num / den;
    Double_t ratioErr = 0;
    if (num > 0)
      ratioErr = TMath::Sqrt(TMath::Power(numE / num, 2.0) + TMath::Power(denE / den, 2.0));
    //
    ratioErr *= ratio;
    eff.first = ratio;
    eff.second = ratioErr;
  }
  Printf("Computed efficiency from %f (%f) / %f (%f) = %f (%f)", num, numE, den, denE, eff.first, eff.second);
  //
  return eff;
}

//----------------------------------------------------------
Int_t GetList(TDirectoryFile* d, TList*& l, TString name, TString suffix)
{
  if (!d) {
    ::Error("MakeTrendingTOFQAv2::GetList", "No input directory was given");
    return -1;
  }
  d->GetObject(Form("%s%s", name.Data(), suffix.Data()), l);
  if (!l) {
    ::Error("MakeTrendingTOFQAv2::GetList", "Cannot find TList %s in Dir. %s\n\tDir. Content:", name.Data(), d->GetName());
    d->ls();
    return -1;
  }
  return 0;
}

//----------------------------------------------------------
TH1F* GetTH1F(TList* l, TString name)
{
  if (!l) {
    ::Error("MakeTrendingTOFQAv2::GetTH1F", "No input TList was given");
    return 0x0;
  }
  TH1F* h = (TH1F*)l->FindObject(name);
  if (!h) {
    ::Error("MakeTrendingTOFQAv2::GetTH1F", "Cannot find TH1F %s in TList %s\n\tDir. Content:", name.Data(), l->GetName());
    l->ls();
    return 0x0;
  }
  return h;
}

//----------------------------------------------------------
TH2F* GetTH2F(TList* l, TString name)
{
  if (!l) {
    ::Error("MakeTrendingTOFQAv2::GetTH2F", "No input TList was given");
    return 0x0;
  }
  TH2F* h = (TH2F*)l->FindObject(name);
  if (!h) {
    ::Error("MakeTrendingTOFQAv2::GetTH2F", "Cannot find TH2F %s in TList %s\n\tDir. Content:", name.Data(), l->GetName());
    l->ls();
    return 0x0;
  }
  return h;
}

//----------------------------------------------------------
TProfile* MakeProfile(TH2F* h, TString name, Float_t min, Float_t max, Bool_t xaxis)
{
  if (!h)
    return nullptr;
  TAxis* axis = xaxis ? h->GetYaxis() : h->GetXaxis();
  TAxis* axislabel = xaxis ? h->GetXaxis() : h->GetYaxis();
  //
  const Int_t minbin = axis->FindBin(min);
  const Int_t maxbin = axis->FindBin(max);
  if (minbin < 1 || minbin > axis->GetNbins())
    ::Error("MakeTrendingTOFQAv2::MakeProfile", "Asking a profile larger than the histogram range: Min: %f outside [%f, %f]", min, axis->GetBinLowEdge(1), axis->GetBinLowEdge(axis->GetNbins()));
  if (maxbin < 1 || maxbin > axis->GetNbins())
    ::Error("MakeTrendingTOFQAv2::MakeProfile", "Asking a profile larger than the histogram range: Max: %f outside [%f, %f]", max, axis->GetBinLowEdge(1), axis->GetBinLowEdge(axis->GetNbins()));

  TProfile* p = nullptr;
  if (xaxis)
    p = h->ProfileX(name, minbin, maxbin);
  else
    p = h->ProfileY(name, minbin, maxbin);
  //
  p->GetYaxis()->SetTitle(axis->GetTitle());
  TString label = axislabel->GetBinLabel(1);
  if (!label.IsNull()) {
    for (Int_t i = 1; i <= axislabel->GetNbins(); i++) {
      p->GetXaxis()->SetBinLabel(i, axislabel->GetBinLabel(i));
    }
  }
  p->SetLineWidth(3);
  p->SetLineColor(1);
  p->SetFillStyle(0);
  return p;
}

//----------------------------------------------------------
void SetupPad(TString opt)
{
#define DoIt(what)         \
  if (opt.Contains(#what)) \
    gPad->what();
  DoIt(SetLogx);
  DoIt(SetLogy);
  DoIt(SetLogz);
  DoIt(SetGridx);
  DoIt(SetGridy);
#undef DoIt
}

//----------------------------------------------------------
Double_t TOFsignal(Double_t* x, Double_t* par)
{
  //Define function to fit TOF signal t-texp
  //as a gaussian + exponential tail
  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t tail = par[3];
  Double_t a = par[4];
  Double_t b = par[5];
  //Double_t c = par[6];
  if (x[0] <= (tail + mean))
    return norm * TMath::Gaus(x[0], mean, sigma) + a + b * x[0]; /*+ c * x[0] * x[0]*/
  else
    return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma)) + a + b * x[0]; /*+ c * x[0] * x[0]*/
}

//----------------------------------------------------------
TF1* TOFsignal(Double_t rangeMin, Double_t rangeMax)
{
  TF1* f = new TF1("fSignalModel", TOFsignal, rangeMin, rangeMax, 6);
  f->SetParNames("Norm", "Mean", "Sigma", "Tail", "Shift", "Slope" /*, "Square"*/);
  f->SetTitle("TOF Signal");
  f->SetParameter(0, 1.);
  f->SetParameter(1, 0.);
  f->SetParLimits(1, -2., 1.);
  f->SetParameter(2, 1.);
  f->SetParLimits(2, 0.5, 2.);
  f->SetParameter(3, 1.);
  f->SetParLimits(3, 0.5, 1.5);
  f->SetParameter(4, 1.);
  f->SetParLimits(4, 0., 1.e8);
  f->SetParameter(5, 0.);
  f->SetParLimits(5, -10., 10.);
  f->SetNpx(2000);
  f->SetLineColor(kRed + 1);
  return f;
}

//----------------------------------------------------------
TList* FitNSigma(TList* l, TString part, TF1* f, TF1* f2, const TString name, const TString suffix)
{
  TH2F* h = GetTH2F(l, Form("%s%s%s", name.Data(), part.Data(), suffix.Data()));
  h->GetYaxis()->SetRangeUser(-5., 5.);
  h->GetXaxis()->SetRangeUser(0.2, 10.);

  TList* lpars = new TList();
  lpars->SetOwner();
  lpars->Add(h);
  //fit with simple gaussian
  h->FitSlicesY(f);
  lpars->Add((TH1D*)gDirectory->Get(Form("%s_1", h->GetName()))->Clone(Form("%s_mean", h->GetName())));
  MakeUpHisto((TH1D*)lpars->Last(), "", "", 1, kBlack, 2);
  lpars->Add((TH1D*)gDirectory->Get(Form("%s_2", h->GetName()))->Clone(Form("%s_pull", h->GetName())));
  MakeUpHisto((TH1D*)lpars->Last(), "", "", 1, kRed, 2);

  //fit with signal model = gaussian + exponential tail
  if (f2) {
    f2->SetParameter(0, h->GetMaximum() * 0.6);
    f2->SetParLimits(0, 0., h->GetMaximum() * 1.2);
    f2->SetParameter(4, h->GetMaximum() * 0.25);
    f2->SetParLimits(4, 0., h->GetMaximum() * 0.5);
    h->FitSlicesY(f2, 0, -1, 0, "QR");
    Int_t i = 1;
    lpars->Add((TH1D*)gDirectory->Get(Form("%s_%i", h->GetName(), i))->Clone(Form("%s_%s", h->GetName(), f2->GetParName(i))));
    MakeUpHisto((TH1D*)lpars->Last(), "", "", 1, kBlue, 2);
    i = 2;
    lpars->Add((TH1D*)gDirectory->Get(Form("%s_%i", h->GetName(), i))->Clone(Form("%s_%s", h->GetName(), f2->GetParName(i))));
    MakeUpHisto((TH1D*)lpars->Last(), "", "", 1, kMagenta + 2, 2);
  }
  return lpars;
  // //***************************************************
  // //Save parameters obtained with the signal model fit
  // // ***************************************************
  // TF1* fSignalModel_bkg = new TF1("fSignalModel_bkg", "pol1", ModelRangeFitNsigmaPIDmin, ModelRangeFitNsigmaPIDmax);
  // fSignalModel_bkg->SetTitle("bkg");
  // fSignalModel_bkg->SetLineColor(kMagenta);
  // const Int_t ptcheck = 100;

  // for (Int_t jj = 0; jj < 3; jj++) {
  //   TCanvas* FitParameters = new TCanvas(Form("FitParameters%i", jj), Form("FitParameters%i", jj));
  //   FitParameters->Divide(2, 4);
  //   FitParameters->cd(1);
  //   f->Draw();
  //   fSignalModel->Draw("same");

  //   for (Int_t cc = 0; cc < npars; cc++) {
  //     FitParameters->cd(cc + 2);
  //     if (par[jj][cc]) {
  //       par[jj][cc]->DrawCopy();
  //       switch (jj) {
  //       case 0:
  //         if (cc == 1)
  //           hSigmaPi_mean->DrawCopy("same");
  //         if (cc == 2)
  //           hSigmaPi_pull->DrawCopy("same");
  //         break;
  //       case 1:
  //         if (cc == 1)
  //           hSigmaKa_mean->DrawCopy("same");
  //         if (cc == 2)
  //           hSigmaKa_pull->DrawCopy("same");
  //         break;
  //       case 2:
  //         if (cc == 1)
  //           hSigmaPro_mean->DrawCopy("same");
  //         if (cc == 2)
  //           hSigmaPro_pull->DrawCopy("same");
  //         break;
  //       default:
  //         break;
  //       }
  //     }
  //   }
  //   FitParameters->SaveAs(Form("TOFQA_FitParameters%i.pdf", jj));
  // }

  // TCanvas* fitcheck = new TCanvas("fitcheck", "fitcheck");
  // fitcheck->Divide(3, 1);
  // for (Int_t jj = 0; jj < 3; jj++) {
  //   fitcheck->cd(jj + 1);
  //   for (Int_t cc = 0; cc < npars; cc++) {
  //     fSignalModel->SetParameter(cc, par[jj][cc]->GetBinContent(ptcheck));
  //   }
  //   TH1D* proj = 0x0;
  //   TH1D* par1 = 0x0;
  //   TH1D* par2 = 0x0;

  //   switch (jj) {
  //   case 0:
  //     proj = (TH1D*)hSigmaPi->ProjectionY(Form("hSigmaPi_par%i", ptcheck), ptcheck, ptcheck);
  //     par1 = hSigmaPi_mean;
  //     par2 = hSigmaPi_pull;
  //     break;
  //   case 1:
  //     proj = (TH1D*)hSigmaKa->ProjectionY(Form("hSigmaKa_par%i", ptcheck), ptcheck, ptcheck);
  //     par1 = hSigmaKa_mean;
  //     par2 = hSigmaKa_pull;
  //     break;
  //   case 2:
  //     proj = (TH1D*)hSigmaPro->ProjectionY(Form("hSigmaPro_par%i", ptcheck), ptcheck, ptcheck);
  //     par1 = hSigmaPro_mean;
  //     par2 = hSigmaPro_pull;
  //     break;
  //   default:
  //     break;
  //   }

  //   f->SetParameter(0, par[jj][0]->GetBinContent(ptcheck));
  //   f->SetParameter(1, par1->GetBinContent(ptcheck));
  //   f->SetParameter(2, par2->GetBinContent(ptcheck));

  //   proj->DrawCopy();
  //   fSignalModel->DrawCopy("same");
  //   f->DrawCopy("same")->SetLineColor(kBlue);
  //   fSignalModel_bkg->SetParameters(fSignalModel->GetParameter(npars - 2), fSignalModel->GetParameter(npars - 1));
  //   fSignalModel_bkg->DrawCopy("same");
  // }
  // gPad->BuildLegend(0.72, 0.67, 0.88, 1);
  // fitcheck->SaveAs("TOFQA_fitcheck.pdf");
}

//----------------------------------------------------------
TLegend* GetLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.04);
  leg->SetShadowColor(0);
  return leg;
}

//----------------------------------------------------------
Int_t GetTDir(TDirectoryFile*& d, TFile*& f, TString name)
{
  f->GetObject(name, d);
  if (!d) {
    ::Error("MakeTrendingTOFQA::GetTDir", "TDirectoryFile %s not present in input file %s.\n", name.Data(), f->GetName());
    return -1;
  }
  return 0;
}

//----------------------------------------------------------
TPaveText* GetPave(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t col, TString txt)
{
  TPaveText* pav = new TPaveText(x1, y1, x2, y2, "NDC");
  pav->SetBorderSize(0);
  pav->SetTextSize(0.045);
  pav->SetFillColor(0); //white background
  pav->SetTextAlign(12);
  pav->SetTextColor(col);
  pav->AddText(txt);
  return pav;
}

#undef CheckAndWrite
#undef CheckAndPrint
#undef CheckAndDraw
