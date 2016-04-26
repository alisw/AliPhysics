/*
  fbellini@cern.ch - last update on 09/03/2016
  Macro to run the TOF QA trending by accessing the std QA output, 
  to be mainly used with the automatic scripts to fill the QA repository.
  Launch with 
  aliroot -l -b -q "MakeTrendingTOFQA.C(\"${fullpath}/QAresults.root\", ${run}, ...) 
  The macro produces a file containing the tree of trending variables and the main plots.
  A feature that displays the plots in canvases must be enable when needed.
*/
/*
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGrid.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPaveText.h"
#include "AliTOFcalib.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TProfile.h"
#include "AliTOFcalibHisto.h"
#include "AliTOFChannelOnlineStatusArray.h"
*/

///Functions with default parameters
Int_t MakeTrendingTOFQAv2(TString qafilename,                   //full path of the QA output;
			  Int_t runNumber,                      //run number
			  Bool_t isMC = kFALSE,                 //MC flag, to disable meaningless checks
			  Bool_t checkPIDqa = kTRUE,            //set to kTRUE to check PIDqa output for TOF
			  Double_t RangeFitNsigmaPIDmin = -2.,  //set lower limit for fitting Nsigma_TOF_ID
			  Double_t RangeFitNsigmaPIDmax = 2.,   //set upper limit for fitting Nsigma_TOF_ID
			  TString ocdbStorage = "raw://",       //set the default ocdb storage
			  Bool_t drawAll = kFALSE,              //enable display plots on canvas and save png
			  Bool_t saveHisto = kTRUE,             //set to kTRUE to save histograms in root file
			  Bool_t savePng = kTRUE );             //set to kTRUE to save histogram to png image

Double_t GetGoodTOFChannelsRatio(Int_t run = -1, Bool_t saveMap = kFALSE, TString OCDBstorage = "raw://", Bool_t inEta08 = kFALSE);

void MakeUpHisto(TH1* histo, TString titleY, Int_t marker = 20, Color_t color = kBlue+2);


Int_t MakeTrendingTOFQA(char * runlist,
			Int_t year = 2010,
			TString period = "LHC10c",
			TString pass = "cpass1_pass4",
			TString nameSuffix = "_barrel",
			Bool_t isMC = kFALSE,
			Int_t trainId = 0,
			Double_t RangeFitNsigmaPIDmin = -2.,
			Double_t RangeFitNsigmaPIDmax = 2.,
			Bool_t saveHisto = kTRUE,
			Bool_t checkPIDqa = kTRUE,
			Bool_t drawAll = kTRUE,
			Bool_t IsOnGrid = kTRUE)
{

  if (IsOnGrid) TGrid::Connect("alien://");

  Int_t filesCounter=0; 
  if (!runlist) {
    printf("Invalid list of runs given as input: nothing done\n");
    return 1;
  }	
  Int_t runNumber=-1;
  char infile[300]; 
  FILE * files = fopen(runlist, "r") ; 
  while (fscanf(files,"%d",&runNumber)==1 ){
    
    //get QAtrain output
    if (trainId==0){
      if (!isMC) sprintf(infile,"alien:///alice/data/%i/%s/000%d/%s/QAresults%s.root",year,period.Data(),runNumber,pass.Data(),nameSuffix.Data());
      else sprintf(infile,"alien:///alice/sim/%i/%s/%d/QAresults%s.root",year,period.Data(),runNumber,nameSuffix.Data());
    } else{
      if (!isMC) sprintf(infile,"alien:///alice/data/%i/%s/000%d/ESDs/%s/QA%i/QAresults%s.root",year,period.Data(),runNumber,pass.Data(),trainId,nameSuffix.Data());
      else sprintf(infile,"alien:///alice/sim/%i/%s/%d/QA%i/QAresults%s.root",year,period.Data(),runNumber,trainId,nameSuffix.Data());
    }
    
    Printf("============== Opening QA file(s) for run %i =======================\n",runNumber);
    
    //run post-analysis
    if (MakeTrendingTOFQAv2(infile, runNumber, isMC, checkPIDqa, RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax, "raw://", drawAll, saveHisto, kTRUE)==0){
      filesCounter++;
    } else Printf("Post analysis not run on QA output %s", infile);
  }
  Printf(":::: Processed %i runs", filesCounter);
  return 0;
}

//---------------------------------------------------------------------------------
Int_t MakeTrendingTOFQAv2(TString qafilename,             //full path of the QA output;
			  Int_t runNumber,                //run number
			  Bool_t isMC,                    //MC flag, to disable meaningless checks
			  Bool_t checkPIDqa,              //set to kTRUE to check PIDqa output for TOF
			  Double_t RangeFitNsigmaPIDmin,  //set lower limit for fitting Nsigma_TOF_ID
			  Double_t RangeFitNsigmaPIDmax,  //set upper limit for fitting Nsigma_TOF_ID		     
			  TString ocdbStorage,            //set the default ocdb storage
			  Bool_t drawAll ,                //enable display plots on canvas and save png
			  Bool_t saveHisto,               //set to kTRUE to save histograms in root file
			  Bool_t savePng)                 //set to kTRUE to save histogram to png image
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
  gStyle->SetTitleBorderSize(0)  ;
  gStyle->SetTitleFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite); 
  gStyle->SetStatBorderSize(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(10);

  TString treePostFileName=Form("trending_%i.root",runNumber);
  TFile * fin = TFile::Open(qafilename,"r");
  if (!fin) {
    Printf("ERROR: QA output not found. Exiting...\n");
    return -1;
  } else {
    Printf("INFO: QA output file %s open. \n",fin->GetName());
  }
  
  //access histograms lists
  char tofQAdirName[15]="TOF";
  char genListName[15]="base_noPID";
  char t0ListName[15]="timeZero_noPID";
  char pidListName[15]="pid_noPID"; 
  char trdListName[15]="trd_noPID";
  char trgListName[15]="trigger_noPID";
  char PIDqaListName[6]="PIDqa";
  TDirectoryFile * tofQAdir=(TDirectoryFile*)fin->Get(tofQAdirName);
  TDirectoryFile * pidQAdir=(TDirectoryFile*)fin->Get(PIDqaListName);
  if (!tofQAdir) {
    Printf("ERROR: TOF QA directory not present in input file.\n");
    return -1;
  }
  TList * generalList=(TList*)tofQAdir->Get(genListName);
  TList  *timeZeroList=(TList*)tofQAdir->Get(t0ListName);
  TList  *pidList=(TList*)tofQAdir->Get(pidListName);
  TList  *pidListT0=0x0;
  TList  *tofPidListT0=0x0;
  if (!pidQAdir) {
    printf("WARNING: PIDqa histograms not available\n");
  } else {
    pidListT0=(TList*)pidQAdir->Get(PIDqaListName);
    tofPidListT0=(TList*)pidListT0->FindObject("TOF");
  }
  TList  *trdList=(TList*)tofQAdir->Get(trdListName);
  TList  *trgList=(TList*)tofQAdir->Get(trgListName);
  
  if (!generalList) Printf("WARNING: general QA histograms absent or not accessible\n");
  if (!timeZeroList) Printf("WARNING: timeZero QA histograms absent or not accessible\n");
  if (!pidList) Printf("WARNING: PID QA histograms absent or not accessible\n");
  if (!trdList) Printf("WARNING: QA histograms for TRD checks absent or not accessible\n");
  if (!trgList) Printf("WARNING: QA histograms for trigger absent or not accessible\n");
  
  if ( (!generalList) && (!timeZeroList) && (!pidList) ){
    printf("ERROR: no QA available \n");
    return -1;
  }
  
  Printf(":::: Getting post-analysis info for run %i",runNumber);
  TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
  Double_t avTime=-9999., peakTime=-9999., spreadTime=-9999., peakTimeErr=-9999., spreadTimeErr=-9999., negTimeRatio=-9999.,
    avRawTime=-9999., peakRawTime=-9999., spreadRawTime=-9999., peakRawTimeErr=-9999., spreadRawTimeErr=-9999., avTot=-9999., peakTot=-9999.,spreadTot=-9999.,  peakTotErr=-9999.,spreadTotErr=-9999.,
    orphansRatio=-9999., avL=-9999., negLratio=-9999.,
    effPt1=-9999., effPt2=-9999., matchEffLinFit1Gev=-9999.,matchEffLinFit1GevErr=-9999.; 
  Double_t avPiDiffTime=-9999.,peakPiDiffTime=-9999., spreadPiDiffTime=-9999.,peakPiDiffTimeErr=-9999., spreadPiDiffTimeErr=-9999.;
  Double_t avT0A=-9999.,peakT0A=-9999., spreadT0A=-9999.,peakT0AErr=-9999., spreadT0AErr=-9999.;
  Double_t avT0C=-9999.,peakT0C=-9999., spreadT0C=-9999.,peakT0CErr=-9999., spreadT0CErr=-9999.;
  Double_t avT0AC=-9999.,peakT0AC=-9999., spreadT0AC=-9999.,peakT0ACErr=-9999., spreadT0ACErr=-9999.;
  Double_t avT0res=-9999.,peakT0res=-9999., spreadT0res=-9999.,peakT0resErr=-9999., spreadT0resErr=-9999.;
  Double_t avT0fillRes=-9999., avT0T0Res=-9999.;
  Float_t avMulti=0;
  Float_t fractionEventsWHits=-9999.;
  /*number of good (HW ok && efficient && !noisy) TOF channels from OCDB*/
  Double_t goodChannelRatio=0.0;
  Double_t goodChannelRatioInAcc=0.0;

  TTree * ttree=new TTree("trending","tree of trending variables");
  ttree->Branch("run",&runNumber,"run/I");
  ttree->Branch("avMulti",&avMulti,"avMulti/F");
  ttree->Branch("fractionEventsWHits",&fractionEventsWHits,"fractionEventsWHits/F");
  ttree->Branch("goodChannelsRatio",&goodChannelRatio,"goodChannelRatio/D");
  ttree->Branch("goodChannelsRatioInAcc",&goodChannelRatioInAcc,"goodChannelRatioInAcc/D");
  ttree->Branch("avTime",&avTime,"avTime/D"); //mean time
  ttree->Branch("peakTime",&peakTime,"peakTime/D"); //main peak time after fit
  ttree->Branch("spreadTime",&spreadTime,"spreadTime/D"); //spread of main peak of time after fit
  ttree->Branch("peakTimeErr",&peakTimeErr,"peakTimeErr/D"); //main peak time after fit error
  ttree->Branch("spreadTimeErr",&spreadTimeErr,"spreadTimeErr/D"); //spread of main peak of time after fit error
  ttree->Branch("negTimeRatio",&negTimeRatio,"negTimeRatio/D"); //negative time ratio
  
  ttree->Branch("avRawTime",&avRawTime,"avRawTime/D"); //mean raw time
  ttree->Branch("peakRawTime",&peakRawTime,"peakRawTime/D"); //mean peak of RAW TIME after fit
  ttree->Branch("spreadRawTime",&spreadRawTime,"spreadRawTime/D"); //spread of main peak of raw time after fit
  ttree->Branch("peakRawTimeErr",&peakRawTimeErr,"peakRawTimeErr/D"); //main peak raw  time after fit error
  ttree->Branch("spreadRawTimeErr",&spreadRawTimeErr,"spreadRawTimeErr/D"); //spread of  raw main peak of time after fit error
  
  ttree->Branch("avTot",&avTot,"avTot/D"); //main peak tot
  ttree->Branch("peakTot",&peakTot,"peakTot/D"); // main peak of tot after fit
  ttree->Branch("spreadTot",&spreadTot,"spreadTot/D"); //spread of main peak of tot after fit
  ttree->Branch("peakTotErr",&peakTotErr,"peakTotErr/D"); // main peak of tot after fit
  ttree->Branch("spreadTotErr",&spreadTotErr,"spreadTotErr/D"); //spread of main peak of tot after fit
  
  ttree->Branch("orphansRatio",&orphansRatio,"orphansRatio/D"); //orphans ratio

  ttree->Branch("avL",&avL,"avL/D"); //mean track length
  ttree->Branch("negLratio",&negLratio,"negLratio/D");//ratio of tracks with track length <350 cm
  ttree->Branch("effPt1",&effPt1,"effPt1/D");//matching eff at 1 GeV/c
  ttree->Branch("effPt2",&effPt2,"effPt2/D"); //matching eff at 2 GeV/c
  ttree->Branch("matchEffLinFit1Gev",&matchEffLinFit1Gev,"matchEffLinFit1Gev/D");//matching eff fit param 
  ttree->Branch("matchEffLinFit1GevErr",&matchEffLinFit1GevErr,"matchEffLinFit1GevErr/D");////matching eff fit param error
  
  ttree->Branch("avPiDiffTime",&avPiDiffTime,"avPiDiffTime/D"); //mean t-texp
  ttree->Branch("peakPiDiffTime",&peakPiDiffTime,"peakPiDiffTime/D"); //main peak t-texp after fit
  ttree->Branch("spreadPiDiffTime",&spreadPiDiffTime,"spreadPiDiffTime/D"); //spread of main peak t-texp after fit
  ttree->Branch("peakPiDiffTimeErr",&peakPiDiffTimeErr,"peakPiDiffTimeErr/D"); //main peak t-texp after fit error
  ttree->Branch("spreadPiDiffTimeErr",&spreadPiDiffTimeErr,"spreadPiDiffTimeErr/D"); //spread of main peak of t-texp after fit error

  ttree->Branch("avT0A",&avT0A,"avT0A/D"); //main peak t0A
  ttree->Branch("peakT0A",&peakT0A,"peakT0A/D"); // main peak of t0A after fit
  ttree->Branch("spreadT0A",&spreadT0A,"spreadTot/D"); //spread of main peak of t0A after fit
  ttree->Branch("peakT0AErr",&peakT0AErr,"peakT0AErr/D"); // main peak of t0A after fit
  ttree->Branch("spreadT0AErr",&spreadT0AErr,"spreadT0AErr/D"); //spread of main peak of t0A after fit

  ttree->Branch("avT0C",&avT0C,"avT0C/D"); //main peak t0C
  ttree->Branch("peakT0C",&peakT0C,"peakT0C/D"); // main peak of t0C after fit
  ttree->Branch("spreadT0C",&spreadT0C,"spreadT0C/D"); //spread of main peak of t0C after fit
  ttree->Branch("peakT0CErr",&peakT0CErr,"peakT0CErr/D"); // main peak of t0C after fit
  ttree->Branch("spreadT0CErr",&spreadT0CErr,"spreadT0CErr/D"); //spread of main peak of t0C after fit
 
  ttree->Branch("avT0AC",&avT0AC,"avT0AC/D"); //main peak t0AC
  ttree->Branch("peakT0AC",&peakT0AC,"peakT0AC/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0AC",&spreadT0AC,"spreadT0AC/D"); //spread of main peak of t0AC after fit
  ttree->Branch("peakT0ACErr",&peakT0ACErr,"peakT0ACErr/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0ACErr",&spreadT0ACErr,"spreadT0ACErr/D"); //spread of main peak of t0AC after fit

  ttree->Branch("avT0res",&avT0res,"avT0res/D"); //main peak t0AC
  ttree->Branch("peakT0res",&peakT0res,"peakT0res/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0res",&spreadT0res,"spreadT0res/D"); //spread of main peak of t0AC after fit
  ttree->Branch("peakT0resErr",&peakT0resErr,"peakT0resErr/D"); // main peak of t0AC after fit
  ttree->Branch("spreadT0resErr",&spreadT0resErr,"spreadT0resErr/D"); //spread of main peak of t0AC after fit
  ttree->Branch("avT0fillRes",&avT0fillRes,"avT0fillRes/D"); //t0 fill res
  ttree->Branch("avT0T0Res",&avT0T0Res,"avT0T0Res/D"); //t0 T0 res

  //save quantities for trending
  goodChannelRatio=(Double_t)GetGoodTOFChannelsRatio(runNumber, kFALSE, ocdbStorage, kFALSE);
  goodChannelRatioInAcc=(Double_t)GetGoodTOFChannelsRatio(runNumber, kFALSE, ocdbStorage, kTRUE);
	
  //--------------------------------- Multiplicity ----------------------------------//

  TH1F * hMulti = (TH1F*) generalList->FindObject("hTOFmulti_all");
  TH1F* hFractionEventsWhits = new TH1F("hFractionEventsWhits","hFractionEventsWhits;fraction of events with hits (%)",200,0.,100.);
  Float_t fraction=0.0;
  if (hMulti->GetEntries()>0.0) {
    fraction = ((Float_t) hMulti->GetBinContent(1))/((Float_t) hMulti->GetEntries());
    avMulti = hMulti->GetMean();
  } else { fraction=0.0; }
  hFractionEventsWhits->Fill(fraction*100.);
  
  //--------------------------------- T0F signal ----------------------------------//
  TH1F * hRawTime = (TH1F*)generalList->FindObject("hRawTime_all");
  if ((hRawTime)&&(hRawTime->GetEntries()>0)){
    avRawTime=hRawTime->GetMean();
    if (!isMC){
      hRawTime->Fit("landau","RQ0","",200.,250.);
      if (hRawTime->GetFunction("landau")) {
	peakRawTime=(hRawTime->GetFunction("landau"))->GetParameter(1);
	spreadRawTime=(hRawTime->GetFunction("landau"))->GetParameter(2);
	peakRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(1);
	spreadRawTimeErr=(hRawTime->GetFunction("landau"))->GetParError(2);
      }
    } else {
      printf("Reminder: Raw time not available in MC simulated data.");
    }
  }
  MakeUpHisto(hRawTime, "matched tracks", 21, kGreen+2);
  
  TH1F * hTime = (TH1F*)generalList->FindObject("hTime_all");
  if ((hTime)&&(hTime->GetEntries()>0)) {
    avTime=hTime->GetMean();
    hTime->Fit("landau","RQ0","",0.,50.);
    if (hTime->GetFunction("landau")) {
      peakTime=(hTime->GetFunction("landau"))->GetParameter(1);
      spreadTime=(hTime->GetFunction("landau"))->GetParameter(2);
      peakTimeErr=(hTime->GetFunction("landau"))->GetParError(1);
      spreadTimeErr=(hTime->GetFunction("landau"))->GetParError(2);
      negTimeRatio=((Double_t)hTime->Integral(1,3)*100.)/((Double_t)hTime->Integral());
    }
    MakeUpHisto(hTime, "matched tracks", 20, kBlue+2);
    
    TLegend *lTime = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
    lTime->SetTextSize(0.04281433);
    lTime->AddEntry(hRawTime, "raw","L");
    lTime->AddEntry(hTime, "ESD","L"); 
    lTime->SetFillColor(kWhite);
    lTime->SetShadowColor(0);
  }
      
  TH1F * hTot = (TH1F*)generalList->FindObject("hTot_all");
  if ((hTot)&&(hTot->GetEntries()>0)){
    avTot=hTot->GetMean();
    hTot->Fit("gaus","","",0.,50.);
    if (hTot->GetFunction("gaus")) {
      peakTot=(hTot->GetFunction("gaus"))->GetParameter(1);
      spreadTot=(hTot->GetFunction("gaus"))->GetParameter(2);
      peakTotErr=(hTot->GetFunction("gaus"))->GetParError(1);
      spreadTotErr=(hTot->GetFunction("gaus"))->GetParError(2);
    }
  }      
  MakeUpHisto(hTot, "matched tracks", 8, kViolet-3);
  
  char orphansTxt[200];
  if (hTot->GetEntries()>1){
    orphansRatio=((Float_t) hTot->GetBinContent(1))/((Float_t) hTot->GetEntries()) ;
  }
  sprintf(orphansTxt,"orphans/matched = %4.2f%%",orphansRatio*100.);
  TPaveText *tOrphans = new TPaveText(0.38,0.63,0.88,0.7, "NDC");
  tOrphans->SetBorderSize(0);
  tOrphans->SetTextSize(0.045);
  tOrphans->SetFillColor(0); //white background
  tOrphans->SetTextAlign(12);
  tOrphans->SetTextColor(kViolet-3);
  tOrphans->AddText(orphansTxt);
  
  TH1F * hL=(TH1F*)generalList->FindObject("hMatchedL_all");
  char negLengthTxt[200];
  if (hL->GetEntries()>0){
    avL=hL->GetMean();
    negLratio=(hL->Integral(1,750))/((Float_t) hL->GetEntries()) ;
  }
  MakeUpHisto(hL, "matched tracks", 1, kBlue+2);
  sprintf(negLengthTxt,"trk with L<350cm /matched = %4.2f%%", negLratio*100.);
  TPaveText *tLength = new TPaveText(0.15,0.83,0.65,0.87, "NDC");
  tLength->SetBorderSize(0);
  tLength->SetTextSize(0.04);
  tLength->SetFillColor(0); //white background
  tLength->SetTextAlign(11);
  tLength->SetTextColor(kOrange-3);
  tLength->AddText(negLengthTxt);

  //--------------------------------- residuals -------------------------------------//
  TH2F* hDxPos4profile = (TH2F*) generalList->FindObject("hMatchedDxVsPt_all");
  TH2F* hTOFmatchedDzVsStrip = (TH2F*)generalList->FindObject("hMatchedDzVsStrip_all");
 
  //--------------------------------- matching eff ----------------------------------//
  //matching as function of pT
  TH1F * hMatchingVsPt = new TH1F("hMatchingVsPt","Matching probability vs. Pt; Pt(GeV/c); matching probability", 50, 0., 5. );
  TH1F * hDenom = (TH1F*)generalList->FindObject("hPrimaryPt_all"); 
  if (hDenom) {  
    hDenom->Sumw2();
    hMatchingVsPt=(TH1F*) generalList->FindObject("hMatchedPt_all")->Clone(); 
    hMatchingVsPt->Sumw2();
    // hMatchingVsPt->Rebin(5);
    // hDenom->Rebin(5);
    hMatchingVsPt->Divide(hDenom);
    hMatchingVsPt->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsPt->SetTitle("TOF matching efficiency as function of transverse momentum");
    hMatchingVsPt->GetYaxis()->SetRangeUser(0,1.2); 
  }
  
  if (hMatchingVsPt->GetEntries()>0){
    hMatchingVsPt->Fit("pol0","","",1.0,10.);
    hMatchingVsPt->Draw();
    if (hMatchingVsPt->GetFunction("pol0")){
      matchEffLinFit1Gev=(hMatchingVsPt->GetFunction("pol0"))->GetParameter(0);
      matchEffLinFit1GevErr=(hMatchingVsPt->GetFunction("pol0"))->GetParError(0);	
      //printf("Matching efficiency fit param is %f +- %f\n",matchEffLinFit1Gev,matchEffLinFit1GevErr );
    }
  } else {
    printf("WARNING: matching efficiency plot has 0 entries. Skipped!\n");
  }
  MakeUpHisto(hMatchingVsPt, "efficiency", 1, kBlue+2);
  
  //matching as function of eta
  TH1F * hMatchingVsEta = new TH1F("hMatchingVsEta","Matching probability vs. #\Eta; #\Eta; matching probability", 20, -1., 1.);
  hDenom->Clear();
  hDenom=(TH1F*)generalList->FindObject("hPrimaryEta_all"); 
  if (hDenom) {  
    hDenom->Sumw2();
    hMatchingVsEta=(TH1F*) generalList->FindObject("hMatchedEta_all")->Clone(); 
    hMatchingVsEta->Sumw2();
    // hMatchingVsEta->Rebin(5);
    // hDenom->Rebin(5);
    hMatchingVsEta->Divide(hDenom);
    hMatchingVsEta->GetXaxis()->SetRangeUser(-1,1);
    hMatchingVsEta->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsEta->GetYaxis()->SetRangeUser(0,1.2);
    hMatchingVsEta->SetTitle("TOF matching efficiency as function of pseudorapidity");
  }
  MakeUpHisto(hMatchingVsEta, "efficiency", 1, kBlue+2);
  
  //matching as function of phi
  TH1F * hMatchingVsPhi = new TH1F("hMatchingVsPhi","Matching probability vs. Phi; Phi(rad); matching probability", 628, 0., 6.28);
  hDenom->Clear();
  hDenom=(TH1F*)generalList->FindObject("hPrimaryPhi_all");  
  if (hDenom) {  
    hDenom->Sumw2();
    hMatchingVsPhi=(TH1F*) generalList->FindObject("hMatchedPhi_all")->Clone(); 
    // hMatchingVsPhi->Rebin(2);
    // hDenom->Rebin(2);    
    hMatchingVsPhi->Sumw2();
    hMatchingVsPhi->Divide(hDenom);
    hMatchingVsPhi->GetYaxis()->SetTitle("matching efficiency");
    hMatchingVsPhi->SetTitle("TOF matching efficiency as function of phi");
    hMatchingVsPhi->GetYaxis()->SetRangeUser(0,1.2);
  }
  MakeUpHisto(hMatchingVsPhi, "efficiency", 1, kBlue+2);

  if (saveHisto) {
    trendFile->cd();
    hMulti->Write();
    hTime->Write();
    hRawTime->Write();
    hTot->Write();
    hL->Write();
    hDxPos4profile->Write();
    hTOFmatchedDzVsStrip->Write();
    hMatchingVsPt->Write();
    hMatchingVsEta->Write();
    hMatchingVsPhi->Write();
  }  
   
//--------------------------------- t-texp ----------------------------------//
  TH2F * hBetaP=(TH2F*)pidList->FindObject("hMatchedBetaVsP_all");
  if (hBetaP) hBetaP->GetYaxis()->SetRangeUser(0.,1.2);
  
  TH1F * hMass=(TH1F*)pidList->FindObject("hMatchedMass_all");
  MakeUpHisto(hMass, "tracks", 1, kBlue+2);
  // hMass->SetFillColor(kAzure+10);
  // hMass->SetFillStyle(1001);
  hMass->Rebin(2);
  
  //pions
  TH1F * hPionDiff=(TH1F*)pidList->FindObject("hExpTimePi_all"); 
  if ((hPionDiff)&&(hPionDiff->GetEntries()>0)) {
    avPiDiffTime=hPionDiff->GetMean();
    hPionDiff->Fit("gaus","","",-1000.,500.);
    if (hPionDiff->GetFunction("gaus")){
      peakPiDiffTime=(hPionDiff->GetFunction("gaus"))->GetParameter(1);
      spreadPiDiffTime=(hPionDiff->GetFunction("gaus"))->GetParameter(2);
      peakPiDiffTimeErr=(hPionDiff->GetFunction("gaus"))->GetParError(1);
      spreadPiDiffTimeErr=(hPionDiff->GetFunction("gaus"))->GetParError(2);
    }
  }
 
  TH1F * hDiffTimeT0fillPion=(TH1F*)pidList->FindObject("hExpTimePiFillSub_all"); 
  TH1F * hDiffTimeT0TOFPion1GeV=(TH1F*)pidList->FindObject("hExpTimePiT0Sub1GeV_all");  
  TH2F * hDiffTimePi=(TH2F*)pidList->FindObject("hExpTimePiVsP_all"); 
  hDiffTimePi->SetTitle("PIONS t-t_{exp,#pi} (from tracking) vs. P");

  //Kaon
  TH2F * hDiffTimeKa=(TH2F*)pidList->FindObject("hExpTimeKaVsP_all");  
  hDiffTimeKa->SetTitle("KAONS t-t_{exp,K} (from tracking) vs. P");

  //Protons
  TH2F * hDiffTimePro=(TH2F*)pidList->FindObject("hExpTimeProVsP_all"); 
  hDiffTimePro->SetTitle("PROTONS t-t_{exp,p} (from tracking) vs. P");
  
  if (saveHisto) {
    trendFile->cd();
    hBetaP->Write();
    hMass->Write();
    hPionDiff->Write();
    hDiffTimeT0fillPion->Write();
    hDiffTimeT0TOFPion1GeV->Write();
    hDiffTimePi->Write();
    hDiffTimeKa->Write();
    hDiffTimePro->Write();
  }

  //--------------------------------- T0 vs multiplicity plots ----------------------------------//
  TH2F * hT0TOFvsNtracks=(TH2F*)timeZeroList->FindObject("hT0TOFvsNtrk");
  Int_t binmin1 = hT0TOFvsNtracks->GetYaxis()->FindBin(-50.);
  Int_t binmax1 = hT0TOFvsNtracks->GetYaxis()->FindBin(50.);
  TProfile* hT0TOFProfile = (TProfile*)hT0TOFvsNtracks->ProfileX("hT0TOFProfile",binmin1,binmax1);
  hT0TOFProfile->SetLineWidth(3);
  hT0TOFProfile->SetLineColor(1);

  TH2F * hT0ACvsNtracks=(TH2F*)timeZeroList->FindObject("hT0ACvsNtrk");
  Int_t binmin2 = hT0ACvsNtracks->GetYaxis()->FindBin(-50.);
  Int_t binmax2 = hT0ACvsNtracks->GetYaxis()->FindBin(50.);
  TProfile* hT0ACProfile = (TProfile*)hT0ACvsNtracks->ProfileX("hT0ACProfile",binmin2,binmax2);
  hT0ACProfile->SetLineWidth(3);
  hT0ACProfile->SetLineColor(1);

  TH2F * hT0AvsNtracks=(TH2F*)timeZeroList->FindObject("hT0AvsNtrk");
  Int_t binmin3 = hT0AvsNtracks->GetYaxis()->FindBin(-50.);
  Int_t binmax3 = hT0AvsNtracks->GetYaxis()->FindBin(50.);
  TProfile* hT0AProfile = (TProfile*)hT0AvsNtracks->ProfileX("hT0AProfile",binmin3,binmax3);
  hT0AProfile->SetLineWidth(3);
  hT0AProfile->SetLineColor(1);

  TH2F * hT0CvsNtracks=(TH2F*)timeZeroList->FindObject("hT0CvsNtrk");
  Int_t binmin4 = hT0CvsNtracks->GetYaxis()->FindBin(-50.);
  Int_t binmax4 = hT0CvsNtracks->GetYaxis()->FindBin(50.);
  TProfile* hT0CProfile = (TProfile*)hT0CvsNtracks->ProfileX("hT0CProfile",binmin4,binmax4);
  hT0CProfile->SetLineWidth(3);
  hT0CProfile->SetLineColor(1);

  TH2F * hStartTime=(TH2F*)timeZeroList->FindObject("hStartTime");
  TH2F * hStartTimeRes=(TH2F*)timeZeroList->FindObject("hStartTimeRes");

  Int_t binminst = hStartTime->GetXaxis()->FindBin(-600.);
  Int_t binmaxst = hStartTime->GetXaxis()->FindBin(600.);

  TProfile* hStartTimeProfile = (TProfile*) hStartTime->ProfileY("hStartTimeProfile",binminst,binmaxst);
  hStartTimeProfile->SetFillStyle(0);
  hStartTimeProfile->SetLineWidth(3);
  hStartTimeProfile->SetLineColor(1);  

  TCanvas *cT0vsMultiplicty = new TCanvas("cT0vsMultiplicity","T0TOF,T0C,T0A,TOAC vs N_TOF,",1200,800);
  cT0vsMultiplicity->Divide(2,2);
  cT0vsMultiplicity->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hT0TOFvsNtracks->Draw("colz");
  hT0TOFProfile->Draw("same");

  cT0vsMultiplicity->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hT0ACvsNtracks->Draw("colz");
  hT0ACProfile->Draw("same");

  cT0vsMultiplicity->cd(3);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hT0AvsNtracks->Draw("colz");
  hT0AProfile->Draw("same");

  cT0vsMultiplicity->cd(4);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hT0CvsNtracks->Draw("colz");
  hT0CProfile->Draw("same");

  TCanvas *cStartTimeRes = new TCanvas("cStartTimeRes","Resolution of start time methods",1200,800);
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  hStartTimeRes->Draw("colz");

  TCanvas *cStartTime = new TCanvas("cStartTime","start time with different methods",1200,800);
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  hStartTime->Draw("colz");

  TString plotDir(".");

  if (savePng) {
    cT0vsMultiplicity->Print(Form("%s/%i_T0vsMultiplicity.png", plotDir.Data(), runNumber));
    cStartTime->Print(Form("%s/%i_StartTimeMethods.png", plotDir.Data(), runNumber));
    cStartTimeRes->Print(Form("%s/%i_StartTimeResolution.png", plotDir.Data(), runNumber));
  }
  if (saveHisto) {
    trendFile->cd();
    hT0TOFvsNtracks->Write();
    hT0TOFProfile->Write();
    hT0ACvsNtracks->Write();
    hT0ACProfile->Write();
    hT0CvsNtracks->Write();
    hT0CProfile->Write();
    hT0AvsNtracks->Write();
    hT0AProfile->Write();
    hStartTime->Write();
    hStartTimeRes->Write();
    hStartTimeProfile->Write();
  }

  //--------------------------------- NSigma PID from TOF QA ----------------------------------//
  TF1 * f = new TF1("f","gaus", RangeFitNsigmaPIDmin, RangeFitNsigmaPIDmax);
  TH2F * hSigmaPi=(TH2F*)pidList->FindObject("hTOFpidSigmaPi_all"); 
  hSigmaPi->GetYaxis()->SetRangeUser(-5.,5.);
  hSigmaPi->GetXaxis()->SetRangeUser(0.2,10.);
  hSigmaPi->FitSlicesY(f);
  TH1D * hSigmaPi_1 = (TH1D*)gDirectory->Get("hTOFpidSigmaPi_all_1");
  TH1D * hSigmaPi_2 = (TH1D*)gDirectory->Get("hTOFpidSigmaPi_all_2");
  //hSigmaPiT0->SetName("hSigmaPiT0");
  hSigmaPi_1->SetLineColor(1);
  hSigmaPi_1->SetLineWidth(2);
  hSigmaPi_2->SetLineColor(2);
  hSigmaPi_2->SetLineWidth(2);

  TH2F * hSigmaKa=(TH2F*)pidList->FindObject("hTOFpidSigmaKa_all"); 
  hSigmaKa->GetYaxis()->SetRangeUser(-5.,5.);
  hSigmaKa->GetXaxis()->SetRangeUser(0.2,10.);
  hSigmaKa->FitSlicesY();
  TH1D * hSigmaKa_1 = (TH1D*)gDirectory->Get("hTOFpidSigmaKa_all_1");
  TH1D * hSigmaKa_2 = (TH1D*)gDirectory->Get("hTOFpidSigmaKa_all_2");
  //hSigmaKaT0->SetName("hSigmaKaT0");
  hSigmaKa_1->SetLineColor(1);
  hSigmaKa_1->SetLineWidth(2);
  hSigmaKa_2->SetLineColor(2);
  hSigmaKa_2->SetLineWidth(2);

  TH2F * hSigmaPro=(TH2F*)pidList->FindObject("hTOFpidSigmaPro_all"); 
  hSigmaPro->GetYaxis()->SetRangeUser(-5.,5.);
  hSigmaPro->GetXaxis()->SetRangeUser(0.2,10.);
  hSigmaPro->FitSlicesY();
  TH1D * hSigmaPro_1 = (TH1D*)gDirectory->Get("hTOFpidSigmaPro_all_1");
  TH1D * hSigmaPro_2 = (TH1D*)gDirectory->Get("hTOFpidSigmaPro_all_2");
  //hSigmaProT0->SetName("hSigmaProT0");
  hSigmaPro_1->SetLineColor(1);
  hSigmaPro_1->SetLineWidth(2);
  hSigmaPro_2->SetLineColor(2);
  hSigmaPro_2->SetLineWidth(2);

  //TString plotDir(".");
 
  TCanvas *cPidPerformance3= new TCanvas("cPidPerformance3","summary of pid performance - sigmas",1200,500);
  cPidPerformance3->Divide(3,1);
  cPidPerformance3->cd(1);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hSigmaPi->GetXaxis()->SetRangeUser(0.2, 5.);
  hSigmaPi->Draw("colz");
  hSigmaPi_1->Draw("same");
  hSigmaPi_2->Draw("same");
  cPidPerformance3->cd(2);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hSigmaKa->GetXaxis()->SetRangeUser(0.2, 5.);
  hSigmaKa->Draw("colz");
  hSigmaKa_1->Draw("same");
  hSigmaKa_2->Draw("same");
  cPidPerformance3->cd(3);
  gPad->SetLogz();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  hSigmaPro->GetXaxis()->SetRangeUser(0.2, 5.);
  hSigmaPro->Draw("colz");
  hSigmaPro_1->Draw("same");
  hSigmaPro_2->Draw("same");

  if (savePng) cPidPerformance3->Print(Form("%s/%i_PID_sigmas.png", plotDir.Data(), runNumber));
  if (saveHisto){
    hSigmaPi->Write();
    hSigmaKa->Write();
    hSigmaPro->Write();
    hSigmaPi_1->Write();
    hSigmaPi_2->Write();
    hSigmaKa_1->Write();
    hSigmaKa_2->Write();
    hSigmaPro_1->Write();
    hSigmaPro_2->Write();
  }
  
 //--------------------------------- NSigma PID from PIDqa ----------------------------------//
  TH2F * hSigmaPiT0 = 0x0;
  TH2F * hSigmaKaT0 = 0x0;
  TH2F * hSigmaProT0 = 0x0;
  if(checkPIDqa && pidQAdir){
    //Pions pid
    hSigmaPiT0=(TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_pion");
    hSigmaPiT0->GetYaxis()->SetRangeUser(-5.,5.);
    hSigmaPiT0->GetXaxis()->SetRangeUser(0.2, 5.);
    hSigmaPiT0->FitSlicesY();
    
    TH1D * hSigmaPiT0_1 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_pion_1");
    TH1D * hSigmaPiT0_2 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_pion_2");
    //hSigmaPiT0->SetName("hSigmaPiT0");
    hSigmaPiT0_1->SetLineColor(1);
    hSigmaPiT0_1->SetLineWidth(2);
    hSigmaPiT0_2->SetLineColor(2);
    hSigmaPiT0_2->SetLineWidth(2);
    
    //Kaons pid
    hSigmaKaT0=(TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_kaon");
    hSigmaKaT0->GetYaxis()->SetRangeUser(-5.,5.);
    hSigmaKaT0->GetXaxis()->SetRangeUser(0.2, 5.);
    hSigmaKaT0->FitSlicesY();
    TH1D * hSigmaKaT0_1 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_kaon_1");
    TH1D * hSigmaKaT0_2 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_kaon_2");
    hSigmaKaT0_1->SetLineColor(1);
    hSigmaKaT0_1->SetLineWidth(2);
    hSigmaKaT0_2->SetLineColor(2);
    hSigmaKaT0_2->SetLineWidth(2);

    //protons pid
    hSigmaProT0 = (TH2F*)tofPidListT0->FindObject("hNsigmaP_TOF_proton");
    hSigmaProT0->GetYaxis()->SetRangeUser(-5.,5.);
    hSigmaProT0->GetXaxis()->SetRangeUser(0.2, 5.);
    hSigmaProT0->FitSlicesY();
    TH1D * hSigmaProT0_1 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_proton_1");
    TH1D * hSigmaProT0_2 = (TH1D*)gDirectory->Get("hNsigmaP_TOF_proton_2");
    hSigmaProT0_1->SetLineColor(1);
    hSigmaProT0_1->SetLineWidth(2);
    hSigmaProT0_2->SetLineColor(2);
    hSigmaProT0_2->SetLineWidth(2);

    TLine *l1=new TLine(0.,0.,5.,0.);
    TLine *l2=new TLine(0.,1.,5.,1.);
    TCanvas *cPidPerformance3T0 = new TCanvas("cPidPerformance3T0","PID performance from PIDqa - sigmas - with StartTime",1200,500);   
    cPidPerformance3T0->Divide(3,1);
    cPidPerformance3T0->cd(1);
    gPad->SetLogz();
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hSigmaPiT0->Draw("colz");
    hSigmaPiT0_1->Draw("same");
    hSigmaPiT0_2->Draw("same");
    l1->Draw("same");
    l2->Draw("same");
    cPidPerformance3T0->cd(2);
    gPad->SetLogz();
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hSigmaKaT0->Draw("colz");
    hSigmaKaT0_1->Draw("same");
    hSigmaKaT0_2->Draw("same");
    l1->Draw("same");
    l2->Draw("same");
    cPidPerformance3T0->cd(3);
    gPad->SetLogz();
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    hSigmaProT0->Draw("colz");
    hSigmaProT0_1->Draw("same");
    hSigmaProT0_2->Draw("same");
    l1->Draw("same");
    l2->Draw("same");
    if (savePng) cPidPerformance3T0->Print(Form("%s/%i_PID_sigmasStartTime.png", plotDir.Data(), runNumber));
    if (saveHisto) {
      trendFile->cd();
      hSigmaPiT0->Write();
      hSigmaPiT0_1->Write();
      hSigmaPiT0_2->Write();
      hSigmaKaT0->Write();
      hSigmaKaT0_1->Write();
      hSigmaKaT0_2->Write();
      hSigmaProT0->Write();
      hSigmaProT0_1->Write();
      hSigmaProT0_2->Write();
    }
  }
  

   //--------------------------------- T0 detector ----------------------------------//
   
  TH1F*hT0A=(TH1F*)timeZeroList->FindObject("hT0A");
  if ((hT0A)&&(hT0A->GetEntries()>0)) {
    avT0A = hT0A->GetMean();
    hT0A->Fit("gaus","RQ0", "", -1000., 1000.);
    if (hT0A->GetFunction("gaus")) {
      peakT0A=(hT0A->GetFunction("gaus"))->GetParameter(1);
      spreadT0A=(hT0A->GetFunction("gaus"))->GetParameter(2);
      peakT0AErr=(hT0A->GetFunction("gaus"))->GetParError(1);
      spreadT0AErr=(hT0A->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0A(gaus): mean = %f +- %f\n",peakT0A,peakT0AErr );
      // printf("Main peak T0A (gaus): spread = %f +- %f\n",spreadT0A,spreadT0AErr );
      //add integral of main peak over total
    }
  }
  MakeUpHisto(hT0A, "events", 8, kBlue);
  hT0A->Rebin(2);

  TH1F*hT0C=(TH1F*)timeZeroList->FindObject("hT0C");
  if ((hT0C)&&(hT0C->GetEntries()>0)) {
    avT0C=hT0C->GetMean();
    hT0C->Fit("gaus","RQ0","", -1000., 1000.);
    if (hT0C->GetFunction("gaus")) {
      peakT0C=(hT0C->GetFunction("gaus"))->GetParameter(1);
      spreadT0C=(hT0C->GetFunction("gaus"))->GetParameter(2);
      peakT0CErr=(hT0C->GetFunction("gaus"))->GetParError(1);
      spreadT0CErr=(hT0C->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0C(gaus): mean = %f +- %f\n",peakT0C,peakT0CErr );
      // printf("Main peak T0C (gaus): spread = %f +- %f\n",spreadT0C,spreadT0CErr );
      //add integral of main peak over total
    }
  }
  MakeUpHisto(hT0C, "events", 8, kGreen+1);
  hT0C->Rebin(2);
	
  TH1F*hT0AC=(TH1F*)timeZeroList->FindObject("hT0AC");
  if ((hT0AC)&&(hT0AC->GetEntries()>0)) {
    avT0AC=hT0AC->GetMean();
    hT0AC->Fit("gaus","RQ0", "",-1000., 1000.);
    if (hT0AC->GetFunction("gaus")) {
      peakT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(1);
      spreadT0AC=(hT0AC->GetFunction("gaus"))->GetParameter(2);
      peakT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(1);
      spreadT0ACErr=(hT0AC->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0AC(gaus): mean = %f +- %f\n",peakT0AC,peakT0ACErr );
      // printf("Main peak T0AC (gaus): spread = %f +- %f\n",spreadT0AC,spreadT0ACErr );	 
    }
  }
  MakeUpHisto(hT0AC, "events", 8, kRed+1);
  hT0AC->Rebin(2);
  
  TLegend *lT0 = new TLegend(0.7125881,0.6052519,0.979435,0.7408306,NULL,"brNDC");
  lT0->SetTextSize(0.041);
  lT0->AddEntry(hT0AC, "T0 A&C","L");
  lT0->AddEntry(hT0A, "T0 A","L"); 
  lT0->AddEntry(hT0C, "T0 C","L");
  lT0->SetFillColor(kWhite);
  lT0->SetShadowColor(0);
  
  TH1F*hT0res=(TH1F*)timeZeroList->FindObject("hT0DetRes");
  if ((hT0res)&&(hT0res->GetEntries()>0)) {
    avT0res=hT0res->GetMean();
    hT0res->Fit("gaus","Q0","");
    if (hT0res->GetFunction("gaus")) {
      peakT0res=(hT0res->GetFunction("gaus"))->GetParameter(1);
      spreadT0res=(hT0res->GetFunction("gaus"))->GetParameter(2);
      peakT0resErr=(hT0res->GetFunction("gaus"))->GetParError(1);
      spreadT0resErr=(hT0res->GetFunction("gaus"))->GetParError(2);	
      // printf("Main peak T0res(gaus): mean = %f +- %f\n",peakT0res,peakT0resErr );
      // printf("Main peak T0res (gaus): spread = %f +- %f\n",spreadT0res,spreadT0resErr );	 
    //add integral of main peak over total
    }
  }
  TH1F*hT0fillRes=(TH1F*)timeZeroList->FindObject("hT0fillRes");
  if ((hT0fillRes)&&(hT0fillRes->GetEntries()>0)) {
    avT0fillRes=hT0fillRes->GetMean();
    }
  TH1F*hT0T0Res=(TH1F*)timeZeroList->FindObject("hT0T0Res");
    if ((hT0T0Res)&&(hT0T0Res->GetEntries()>0)) {
      avT0T0Res=hT0T0Res->GetMean();
    }
  
  if (saveHisto) {
    trendFile->cd(); 
    hT0AC->Write();
    hT0A->Write();
    hT0C->Write();
    hT0res->Write();
    //hT0fillRes->Write();
    //hT0T0Res->Write();
  }	
  //Fill tree and save to file
  ttree->Fill();
  printf("==============  Saving trending quantities in tree for run %i ===============\n",runNumber);
  trendFile->cd();
  ttree->Write();
  trendFile->Close();
  
  //close input file
  fin->Close();

  TCanvas *cTrackProperties = 0x0;
  TCanvas *cProfile = 0x0;
  TCanvas *cMatchingPerformance = 0x0;
  TCanvas *cPidPerformance = 0x0;
  TCanvas *cPidPerformance2 = 0x0;
  TCanvas *cT0detector = 0x0;

  if (drawAll){
    cTrackProperties= new TCanvas("cTrackProperties","summary of matched tracks properties", 1200, 500);
    cTrackProperties->Divide(3,1);
    cTrackProperties->cd(1);
    gPad->SetLogy();
    hTime->Draw("");
    hRawTime ->Draw("same");
    cTrackProperties->cd(2);
    gPad->SetLogy();
    hTot->Draw("");
    tOrphans->Draw(); 
    cTrackProperties->cd(3);
    gPad->SetLogy(); 
    hL->Draw("");
    tLength->Draw(); 

    cPidPerformance= new TCanvas("cPidPerformance","summary of pid performance", 900,500);
    cPidPerformance->Divide(2,1);
    cPidPerformance->cd(1);
    gPad->SetLogz();
    hBetaP->Draw("colz");   
    cPidPerformance->cd(2);
    gPad->SetLogy();
    hMass->Draw("HIST");

    cT0detector= new TCanvas("cT0detector","T0 detector",800,600);
    cT0detector->Divide(2,1);
    cT0detector->cd(1);
    gPad->SetGridx();
    hT0AC->Draw("");
    hT0AC->SetTitle("timeZero measured by T0 detector");
    hT0A ->Draw("same");
    hT0C ->Draw("same");
    //lT0->Draw();  
    cT0detector->cd(2);
    hT0res->Draw();

    // TCanvas *cResiduals= new TCanvas("residuals","residuals", 900,500);
    // cResiduals->Divide(2,1);
    // cResiduals->cd(1);
    // gPad->SetLogz();
    // hDxPos4profile->GetYaxis()->SetRangeUser(-5.,5.);
    // hDxPos4profile->Draw("colz");
    // // profDxPos->SetLineColor(kRed);
    // // profDxPos->Draw("same");
    // cResiduals->cd(2);
    // gPad->SetLogz();
    // hDxNeg4profile->GetYaxis()->SetRangeUser(-5.,5.); 
    // hDxNeg4profile->Draw("colz");
    // // profDxNeg->SetLineColor(kBlue);
    // // profDxNeg->Draw("same"); 

    if (savePng) {
      cTrackProperties->Print(Form("%s/%i_TrackProperties.png",plotDir.Data(), runNumber));
      cPidPerformance->Print(Form("%s/%i_PID.png",plotDir.Data(), runNumber));
      cT0detector->Print(Form("%s/%i_T0Detector.png",plotDir.Data(), runNumber));
    }
  }
   
  cProfile = new TCanvas("cProfile","cProfile",50,50,750,550);
  cProfile->cd();
  gPad->SetLogz();
  hTOFmatchedDzVsStrip->Draw("colz");
  Int_t binmin = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(-3);
  Int_t binmax = hTOFmatchedDzVsStrip->GetYaxis()->FindBin(3);
  TProfile* hDzProfile = (TProfile*)hTOFmatchedDzVsStrip->ProfileX("hDzProfile",binmin, binmax);
  hDzProfile->SetLineWidth(3);
  hDzProfile->Draw("same");
    
  cMatchingPerformance= new TCanvas("cMatchingPerformance","summary of matching performance",1200,500);
  cMatchingPerformance->Divide(3,1);
  cMatchingPerformance->cd(1);
  hMatchingVsPt->Draw();
  cMatchingPerformance->cd(2);
  hMatchingVsEta->Draw();
  cMatchingPerformance->cd(3);
  hMatchingVsPhi->Draw();
    
  cPidPerformance2= new TCanvas("cPidPerformance2","summary of pid performance - expected times", 1200, 500);
  cPidPerformance2->Divide(3,1);
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
    cProfile->Print(Form("%s/%i_ProfileDZvsStripNumber.png",plotDir.Data(), runNumber));
    //cMatchingPerformance->Print(Form("%s/%i_Matching.png",plotDir.Data(), runNumber));
    cPidPerformance2->Print(Form("%s/%i_PID_ExpTimes.png",plotDir.Data(), runNumber));
  }
  return  0;
}


//----------------------------------------------------------
Double_t GetGoodTOFChannelsRatio(Int_t run, Bool_t saveMap , TString OCDBstorage , Bool_t inEta08)
{
  /*
    It retrieves from OCDB the number of good (= efficient && not noisy && HW ok) TOF channels.
    Optionally is saves the channel map
  */
  if (run<=0) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): invalid run number. Please set a run number.\n"); 
    return 0.0;
  }
  
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(OCDBstorage.Data());
  cdb->SetRun(run);
  
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/Status");
  if (!cdbe) {
    printf("MakeTrendingTOFqa.C - ERROR in CheckCalibStatus(): OCDB entry not available. Please, try again.\n");
    return 0.0;
  }  

  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  TH2F *hOkMap = new TH2F("hOkMap", "Ok map (!noisy & !problematic & efficient);sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2F *hOkMapInAcceptance = new TH2F("hOkMapInAcceptance",Form( "Good channels in |#eta|<0.8 - run %i;sector;strip",run), 72, 0., 18., 91, 0., 91.);
  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibHisto();
  AliTOFcalib calib;
  calib.Init();
  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  for (Int_t i = 0; i <  array->GetSize(); i++) {
    sector = calibHisto.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calibHisto.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calibHisto.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    if ( !(array->GetNoiseStatus(i) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad)   &&
	 (calib.IsChannelEnabled(i,kTRUE,kTRUE))) {
      hOkMap->Fill(hitmapx,hitmapy);
      /* 3 strips / side excluded from norm. factor to consider |eta|<0.8 */
      if (sectorStrip>2 && sectorStrip<89) hOkMapInAcceptance->Fill(hitmapx,hitmapy);
    }
  }
  Int_t nOk = (Int_t) hOkMap->GetEntries();
  Int_t nOkInAcc = (Int_t) hOkMapInAcceptance->GetEntries();
  Double_t ratioOk = nOk/152928.;
  /* 96 channels * 6 strips * 18 modules excluded from norm. factor to consider |eta|<0.8 */
  Double_t ratioOkInAcc = nOkInAcc/(152928.-6.*18.*96.); 
  if (saveMap) {
    hOkMap->SaveAs(Form("run%i_OKChannelsMap.root",run));
    hOkMapInAcceptance->SaveAs(Form("run%i_OKChannelsInAccMap.root",run));
  }
  cout << "###    Run " << run << ": TOF channels ok = " << nOk << "/ total 152928 channels = " << ratioOk*100. << "% of whole TOF" << endl;
  cout << "###    Run " << run << ": TOF channels in acc. ok = " << nOkInAcc << "/ total 142560 channels = " << ratioOkInAcc*100. << "% of whole TOF" << endl;
  if (inEta08) 
    return ratioOkInAcc;
  else
    return ratioOk;

}

//----------------------------------------------------------
void MakeUpHisto(TH1* histo, TString titleY, Int_t marker, Color_t color)
{
  if (!histo) return;
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(0.7);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetFillColor(kWhite);
  histo->SetFillStyle(0); 
  histo->GetYaxis()->SetTitle(titleY.Data());
  histo->GetYaxis()->SetTitleOffset(1.35);
  histo->GetXaxis()->SetLabelSize(0.03);
  
  return;
}
