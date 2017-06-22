#ifndef __CINT__
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliDAQ.h"
#endif
#include "map"
using namespace std;
// TODO read number of bits from AliVEvent?
#define NBITS 29
#define NMAXCLASSES 100

TString bitNames[NBITS] = {
    "kINT1",
    "kINT7",
    "kMUON",
    "kHighMultSPD",
    "kEMC1",
    "kINT5",
    "kINT7inMUON",
    "kMuonSingleHighPt7",
    "kMuonLikeLowPt7",
    "kMuonUnlikeLowPt7",
    "kEMC7",
    "kMuonSingleLowPt7",
    "kPHI1",
    "kPHI78",
    "kEMCEJE",
    "kEMCEGA",
    "kHighMultV0",
    "kSemiCentral",
    "kDG",
    "kZED",
    "kSPI78",
    "kINT8",
    "kMuonSingleLowPt8",
    "kMuonSingleHighPt8",
    "kMuonLikeLowPt8",
    "kMuonUnlikeLowPt8",
    "kMuonUnlikeLowPt0",
    "kUserDefined",
    "kTRD"
};

void SetHisto(TH1D* h, Bool_t setMinimumToZero=0);
void SetHisto(TH2D* h, Bool_t setMinimumToZero=0);
void AddFillSeparationLines(TH1* h, map<Int_t,Int_t> &fills);
void AddPeriodSeparationLines(TH1* h,  map<Int_t,TString> &periods);

//void periodLevelQA(TString inputFileName ="/afs/cern.ch/work/a/aliqaevs/www/data/2012/LHC12h/pass1/trending.root"){
void periodLevelQA(TString inputFileName ="trending.root"){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(1);

  TFile* f = new TFile(inputFileName.Data());
  TTree* t = (TTree*) f->Get("trending");
  TObjArray* classes = new TObjArray();
  TObjString* activeDetectors = new TObjString();
  TObjString* partition = new TObjString();
  TObjString* lhcState = new TObjString();
  TObjString* lhcPeriod = new TObjString();
  Int_t run                = 0;
  Int_t fill               = 0;
  Double_t run_duration    = 0;
  Int_t nBCsPerOrbit       = 0;
  Double_t refCounts       = 0;
  Double_t mu              = 0;
  Double_t lumi_seen       = 0;
  Double_t interactionRate = 0;
  ULong64_t class_l0b[NMAXCLASSES]         = {0};
  ULong64_t class_l0a[NMAXCLASSES]         = {0};
  ULong64_t class_l1b[NMAXCLASSES]         = {0};
  ULong64_t class_l1a[NMAXCLASSES]         = {0};
  ULong64_t class_l2b[NMAXCLASSES]         = {0};
  ULong64_t class_l2a[NMAXCLASSES]         = {0};
  Double_t  class_lifetime[NMAXCLASSES]    = {0};
  Double_t  class_lumi[NMAXCLASSES]        = {0};
  ULong64_t alias_recorded[NBITS]          = {0};
  ULong64_t alias_reconstructed[NBITS]     = {0};
  ULong64_t alias_accepted[NBITS]          = {0};
  ULong64_t alias_acc_step1[NBITS]         = {0};
  ULong64_t alias_acc_step2[NBITS]         = {0};
  ULong64_t alias_acc_step3[NBITS]         = {0};
  ULong64_t alias_acc_step4[NBITS]         = {0};
  ULong64_t alias_acc_step5[NBITS]         = {0};
  ULong64_t alias_acc_step6[NBITS]         = {0};
  ULong64_t alias_acc_step7[NBITS]         = {0};
  ULong64_t alias_acc_step8[NBITS]         = {0};
  ULong64_t alias_acc_step9[NBITS]         = {0};
  Double_t alias_l0b_rate[NBITS]           = {0};
  Double_t alias_lifetime[NBITS]           = {0};
  Double_t alias_lumi_recorded[NBITS]      = {0};
  Double_t alias_lumi_reconstructed[NBITS] = {0};
  Double_t alias_lumi_accepted[NBITS]      = {0};
  Int_t timeStart = 0;
  Int_t timeEnd = 0;
  Double_t meanV0MOn = 0;
  Double_t meanV0MOf = 0;
  Double_t meanOFO = 0;
  Double_t meanTKL = 0;
  Double_t meanErrV0MOn = 0;
  Double_t meanErrV0MOf = 0;
  Double_t meanErrOFO = 0;
  Double_t meanErrTKL = 0;
  Double_t meanV0MOnHM = 0;
  Double_t meanV0MOfHM = 0;
  Double_t meanOFOHM = 0;
  Double_t meanTKLHM = 0;
  Double_t meanErrV0MOnHM = 0;
  Double_t meanErrV0MOfHM = 0;
  Double_t meanErrOFOHM = 0;
  Double_t meanErrTKLHM = 0;
  Double_t thresholdV0M = 0;
  Double_t aV0MOnVsOfVal = 0;
  Double_t bV0MOnVsOfVal = 0;
  Double_t aV0MOnVsOfErr = 0;
  Double_t bV0MOnVsOfErr = 0;
  Double_t aSPDOnVsOfVal = 0;
  Double_t bSPDOnVsOfVal = 0;
  Double_t aSPDOnVsOfErr = 0;
  Double_t bSPDOnVsOfErr = 0;
  Double_t aV0MOnVsOfOADB = 0;
  Double_t bV0MOnVsOfOADB = 0;
  Double_t aSPDOnVsOfOADB = 0;
  Double_t bSPDOnVsOfOADB = 0;

  TH2F* hHistStat = new TH2F();
  
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("fill",&fill);
  t->SetBranchAddress("bcs",&nBCsPerOrbit);
  t->SetBranchAddress("run_duration",&run_duration);
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("interactionRate",&interactionRate);
  t->SetBranchAddress("refCounts",&refCounts);
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("class_l0b",&class_l0b);
  t->SetBranchAddress("class_l0a",&class_l0a);
  t->SetBranchAddress("class_l1b",&class_l1b);
  t->SetBranchAddress("class_l1a",&class_l1a);
  t->SetBranchAddress("class_l2b",&class_l2b);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->SetBranchAddress("class_lifetime",&class_lifetime);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("alias_recorded",&alias_recorded);
  t->SetBranchAddress("alias_reconstructed",&alias_reconstructed);
  t->SetBranchAddress("alias_accepted",&alias_accepted);
  t->SetBranchAddress("alias_acc_step1",&alias_acc_step1);
  t->SetBranchAddress("alias_acc_step2",&alias_acc_step2);
  t->SetBranchAddress("alias_acc_step3",&alias_acc_step3);
  t->SetBranchAddress("alias_acc_step4",&alias_acc_step4);
  t->SetBranchAddress("alias_acc_step5",&alias_acc_step5);
  t->SetBranchAddress("alias_acc_step6",&alias_acc_step6);
  t->SetBranchAddress("alias_acc_step7",&alias_acc_step7);
  t->SetBranchAddress("alias_acc_step8",&alias_acc_step8);
  t->SetBranchAddress("alias_acc_step9",&alias_acc_step9);
  t->SetBranchAddress("alias_l0b_rate",&alias_lifetime);
  t->SetBranchAddress("alias_lifetime",&alias_lifetime);
  t->SetBranchAddress("alias_lumi_recorded",&alias_lumi_recorded);
  t->SetBranchAddress("alias_lumi_reconstructed",&alias_lumi_reconstructed);
  t->SetBranchAddress("alias_lumi_accepted",&alias_lumi_accepted);
  t->SetBranchAddress("activeDetectors",&activeDetectors);
  t->SetBranchAddress("partition",&partition);
  t->SetBranchAddress("lhcState",&lhcState);
  t->SetBranchAddress("lhcPeriod",&lhcPeriod);
  t->SetBranchAddress("timeStart",&timeStart);
  t->SetBranchAddress("timeEnd",&timeEnd);
  t->SetBranchAddress("meanV0MOn",&meanV0MOn);
  t->SetBranchAddress("meanV0MOf",&meanV0MOf);
  t->SetBranchAddress("meanOFO",&meanOFO);
  t->SetBranchAddress("meanTKL",&meanTKL);
  t->SetBranchAddress("meanErrV0MOn",&meanErrV0MOn);
  t->SetBranchAddress("meanErrV0MOf",&meanErrV0MOf);
  t->SetBranchAddress("meanErrOFO",&meanErrOFO);
  t->SetBranchAddress("meanErrTKL",&meanErrTKL);
  t->SetBranchAddress("meanV0MOnHM",&meanV0MOnHM);
  t->SetBranchAddress("meanV0MOfHM",&meanV0MOfHM);
  t->SetBranchAddress("meanOFOHM",&meanOFOHM);
  t->SetBranchAddress("meanTKLHM",&meanTKLHM);
  t->SetBranchAddress("meanErrV0MOnHM",&meanErrV0MOnHM);
  t->SetBranchAddress("meanErrV0MOfHM",&meanErrV0MOfHM);
  t->SetBranchAddress("meanErrOFOHM",&meanErrOFOHM);
  t->SetBranchAddress("meanErrTKLHM",&meanErrTKLHM);
  t->SetBranchAddress("thresholdV0M",&thresholdV0M);
  t->SetBranchAddress("aV0MOnVsOfVal",&aV0MOnVsOfVal);
  t->SetBranchAddress("bV0MOnVsOfVal",&bV0MOnVsOfVal);
  t->SetBranchAddress("aV0MOnVsOfErr",&aV0MOnVsOfErr);
  t->SetBranchAddress("bV0MOnVsOfErr",&bV0MOnVsOfErr);
  t->SetBranchAddress("aSPDOnVsOfVal",&aSPDOnVsOfVal);
  t->SetBranchAddress("bSPDOnVsOfVal",&bSPDOnVsOfVal);
  t->SetBranchAddress("aSPDOnVsOfErr",&aSPDOnVsOfErr);
  t->SetBranchAddress("bSPDOnVsOfErr",&bSPDOnVsOfErr);
  t->SetBranchAddress("aV0MOnVsOfOADB",&aV0MOnVsOfOADB);
  t->SetBranchAddress("bV0MOnVsOfOADB",&bV0MOnVsOfOADB);
  t->SetBranchAddress("aSPDOnVsOfOADB",&aSPDOnVsOfOADB);
  t->SetBranchAddress("bSPDOnVsOfOADB",&bSPDOnVsOfOADB);

  t->SetBranchAddress("hHistStat",&hHistStat);

  Int_t nRuns = t->GetEntries();
  const Int_t nDetectors=19;
  TH2D* hActiveDetectors    = new TH2D("hActiveDetectors","Active detectors",nRuns,0,nRuns,nDetectors,0,nDetectors);
  TH1D* hInteractionRate    = new TH1D("hInteractionRate","INEL interaction rate [Hz]",nRuns,0,nRuns);
  TH1D* hMu                 = new TH1D("hMu","Average number of INEL collisions per BC",nRuns,0,nRuns);
  TH1D* hBCs                = new TH1D("hBCs","Number of colliding bunches",nRuns,0,nRuns);
  TH1D* hDuration           = new TH1D("hDuration","Duration, s",nRuns,0,nRuns);
  TH1D* hLumiSeen           = new TH1D("hLumiSeen","Luminosity seen, nb-1",nRuns,0,nRuns);
  TH2D* hClassL0BvsRun      = new TH2D("hClassL0BVsRun","Class L0B vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassL2AvsRun      = new TH2D("hClassL2AVsRun","Class L2A vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassLifetimeVsRun = new TH2D("hClassLifetimeVsRun","Lifetime class-by-class vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassLumiVsRun     = new TH2D("hClassLumiVsRun","Luminosity class-by-class vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hRecorded           = new TH2D("hRecorded","Recorded",nRuns,0,nRuns,1,0,1);  
  TH2D* hReconstructed      = new TH2D("hReconstructed","Reconstructed",nRuns,0,nRuns,1,0,1);  
  TH2D* hAccepted           = new TH2D("hAccepted","Accepted",nRuns,0,nRuns,1,0,1);  
  TH2D* hAccStep1           = new TH2D("hAccStep1","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep2           = new TH2D("hAccStep2","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep3           = new TH2D("hAccStep3","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep4           = new TH2D("hAccStep4","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep5           = new TH2D("hAccStep5","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep6           = new TH2D("hAccStep6","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep7           = new TH2D("hAccStep7","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep8           = new TH2D("hAccStep8","",nRuns,0,nRuns,1,0,1);
  TH2D* hAccStep9           = new TH2D("hAccStep9","",nRuns,0,nRuns,1,0,1);
  TH2D* hLumiRecorded       = new TH2D("hLumiRecorded","Lumi recorded",nRuns,0,nRuns,1,0,1);  
  TH2D* hLumiReconstructed  = new TH2D("hLumiReconstructed","Lumi reconstructed",nRuns,0,nRuns,1,0,1);  
  TH2D* hLumiAccepted       = new TH2D("hLumiAccepted","Lumi accepted",nRuns,0,nRuns,1,0,1);  
  TH1D* hMeanV0MOn          = new TH1D("hMeanV0MOn","Mean V0M online in kINT7",nRuns,0,nRuns);
  TH1D* hMeanV0MOf          = new TH1D("hMeanV0MOf","Mean V0M offline in kINT7",nRuns,0,nRuns);
  TH1D* hMeanOFO            = new TH1D("hMeanOFO","Mean fired outer FO chips in kINT7",nRuns,0,nRuns);
  TH1D* hMeanTKL            = new TH1D("hMeanTKL","Mean number of tracklets in kINT7",nRuns,0,nRuns);
  TH1D* hMeanV0MOnHM        = new TH1D("hMeanV0MOnHM","Mean V0M online in kHighMultV0",nRuns,0,nRuns);
  TH1D* hMeanV0MOfHM        = new TH1D("hMeanV0MOfHM","Mean V0M offline in kHighMultV0",nRuns,0,nRuns);
  TH1D* hMeanOFOHM          = new TH1D("hMeanOFOHM","Mean fired outer FO chips in kHighMultV0",nRuns,0,nRuns);
  TH1D* hMeanTKLHM          = new TH1D("hMeanTKLHM","Mean number of tracklets in kHighMultV0",nRuns,0,nRuns);
  TH1D* hThresholdV0M       = new TH1D("hThresholdV0M","V0M threshold",nRuns,0,nRuns);
  TH1D* hThresholdV0Mnorm   = new TH1D("hThresholdV0Mnorm","V0M threshold / <V0M>",nRuns,0,nRuns);
  TH1D* hV0MOnVsOfA         = new TH1D("hV0MOnVsOfA","V0MOnVsOf cut parameter A",nRuns,0,nRuns);
  TH1D* hV0MOnVsOfB         = new TH1D("hV0MOnVsOfB","V0MOnVsOf cut parameter B",nRuns,0,nRuns);
  TH1D* hSPDOnVsOfA         = new TH1D("hSPDOnVsOfA","SPDOnVsOf cut parameter A",nRuns,0,nRuns);
  TH1D* hSPDOnVsOfB         = new TH1D("hSPDOnVsOfB","SPDOnVsOf cut parameter B",nRuns,0,nRuns);
  TH1D* hV0MOnVsOfA_OADB    = new TH1D("hV0MOnVsOfA_OADB","V0MOnVsOf cut parameter A in OADB",nRuns,0,nRuns);
  TH1D* hV0MOnVsOfB_OADB    = new TH1D("hV0MOnVsOfB_OADB","V0MOnVsOf cut parameter B in OADB",nRuns,0,nRuns);
  TH1D* hSPDOnVsOfA_OADB    = new TH1D("hSPDOnVsOfA_OADB","SPDOnVsOf cut parameter A in OADB",nRuns,0,nRuns);
  TH1D* hSPDOnVsOfB_OADB    = new TH1D("hSPDOnVsOfB_OADB","SPDOnVsOf cut parameter B in OADB",nRuns,0,nRuns);
  
  map<Int_t,Int_t> fills;
  map<Int_t,TString> periods;

  TString detName[nDetectors]={"ACORDE","PMD","FMD","HMPID","CPV","PHOS","EMCAL","MUONTRK","MUONTRG","T0","VZERO","AD","ZDC","ITSSPD","ITSSDD","ITSSSD","TPC","TRD","TOF"};
  for (Int_t iDet=0;iDet<nDetectors;iDet++) hActiveDetectors->GetYaxis()->SetBinLabel(iDet+1,detName[iDet].Data());
  hClassL0BvsRun     ->SetBit(TH1::kCanRebin);
  hClassL2AvsRun     ->SetBit(TH1::kCanRebin);
  hClassLifetimeVsRun->SetBit(TH1::kCanRebin);
  hClassLumiVsRun    ->SetBit(TH1::kCanRebin);
  hRecorded          ->SetBit(TH1::kCanRebin);
  hReconstructed     ->SetBit(TH1::kCanRebin);
  hAccepted          ->SetBit(TH1::kCanRebin);
  hLumiRecorded      ->SetBit(TH1::kCanRebin);
  hLumiReconstructed ->SetBit(TH1::kCanRebin);
  hLumiAccepted      ->SetBit(TH1::kCanRebin);

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("local:///cvmfs/alice.cern.ch/calibration/data/2017/OCDB");
//  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("raw://");

  for (Int_t r=0;r<nRuns;r++){
    t->GetEntry(r);
//    if (run==253660) continue;
//    if (run>=255042 && run<=255076) continue;
//    if (TMath::Abs(aV0MOnVsOfVal)<1e-5) continue;
//    if (TMath::Abs(bV0MOnVsOfVal)<1e-5) continue;
//    if (TMath::Abs(aSPDOnVsOfVal)<1e-5) continue;
//    if (TMath::Abs(bSPDOnVsOfVal)<1e-5) continue;    

//    if (!partition->String().Contains("PHYSICS_1")) continue;
//    if (!lhcState->String().Contains("STABLE")) continue;
//    if (!lhcPeriod->String().Contains("LHC15o")) continue;
    Double_t thr = 0;
    if (1) {
      man->SetRun(run);
      AliCDBEntry* triggerEntry = man->Get("VZERO/Trigger/Data");
      AliVZEROTriggerData* trigData = triggerEntry ? (AliVZEROTriggerData*) triggerEntry->GetObject() : 0;
      thr = trigData ? trigData->GetCentralityV0AThrLow() : 0;
    }
    char* srun = Form("%i",run);
    hInteractionRate->Fill(srun,interactionRate);
    hMu             ->Fill(srun,mu);
    hBCs            ->Fill(srun,nBCsPerOrbit);
    hDuration       ->Fill(srun,run_duration);
    hLumiSeen       ->Fill(srun,lumi_seen);
    hMeanV0MOn      ->Fill(srun,meanV0MOn);
    hMeanV0MOf      ->Fill(srun,meanV0MOf);
    hMeanOFO        ->Fill(srun,meanOFO);
    hMeanTKL        ->Fill(srun,meanTKL);
    hMeanV0MOn->SetBinError(hMeanV0MOn->GetXaxis()->FindBin(srun),meanErrV0MOn);
    hMeanV0MOf->SetBinError(hMeanV0MOf->GetXaxis()->FindBin(srun),meanErrV0MOf);
    hMeanOFO  ->SetBinError(hMeanOFO->GetXaxis()->FindBin(srun)  ,meanErrOFO);
    hMeanTKL  ->SetBinError(hMeanTKL->GetXaxis()->FindBin(srun)  ,meanErrTKL);
    hMeanV0MOnHM   ->Fill(srun,meanV0MOnHM);
    hMeanV0MOfHM   ->Fill(srun,meanV0MOfHM);
    hMeanOFOHM     ->Fill(srun,meanOFOHM);
    hMeanTKLHM     ->Fill(srun,meanTKLHM);
    hMeanV0MOnHM->SetBinError(hMeanV0MOnHM->GetXaxis()->FindBin(srun),meanErrV0MOnHM);
    hMeanV0MOfHM->SetBinError(hMeanV0MOfHM->GetXaxis()->FindBin(srun),meanErrV0MOfHM);
    hMeanOFOHM  ->SetBinError(hMeanOFOHM->GetXaxis()->FindBin(srun)  ,meanErrOFOHM);
    hMeanTKLHM  ->SetBinError(hMeanTKLHM->GetXaxis()->FindBin(srun)  ,meanErrTKLHM);
    hV0MOnVsOfA    ->Fill(srun,aV0MOnVsOfVal);
    hV0MOnVsOfB    ->Fill(srun,bV0MOnVsOfVal);
    hSPDOnVsOfA    ->Fill(srun,aSPDOnVsOfVal);
    hSPDOnVsOfB    ->Fill(srun,bSPDOnVsOfVal);
    hV0MOnVsOfA->SetBinError(hV0MOnVsOfA->GetXaxis()->FindBin(srun),aV0MOnVsOfErr);
    hV0MOnVsOfB->SetBinError(hV0MOnVsOfB->GetXaxis()->FindBin(srun),bV0MOnVsOfErr);
    hSPDOnVsOfA->SetBinError(hSPDOnVsOfA->GetXaxis()->FindBin(srun),aSPDOnVsOfErr);
    hSPDOnVsOfB->SetBinError(hSPDOnVsOfB->GetXaxis()->FindBin(srun),bSPDOnVsOfErr);
    hV0MOnVsOfA_OADB->Fill(srun,aV0MOnVsOfOADB);
    hV0MOnVsOfB_OADB->Fill(srun,bV0MOnVsOfOADB);
    hSPDOnVsOfA_OADB->Fill(srun,aSPDOnVsOfOADB);
    hSPDOnVsOfB_OADB->Fill(srun,bSPDOnVsOfOADB);
    
    hThresholdV0M->Fill(srun,thr);
    hThresholdV0Mnorm->Fill(srun,meanV0MOn>1e-4 ? thr/meanV0MOn : 0);
    hThresholdV0Mnorm->SetBinError(hThresholdV0Mnorm->GetXaxis()->FindBin(srun),meanV0MOn>1e-4 ? thr/meanV0MOn*meanErrV0MOn/meanV0MOn : 0);
    printf("threshold %f\n",thr);
    fills[run]=fill;
    periods[run]=lhcPeriod->String();

    for (Int_t iDet=0;iDet<nDetectors;iDet++) {
     hActiveDetectors->Fill(srun,detName[iDet].Data(),activeDetectors->String().Contains(detName[iDet].Data()));
    }
    for (Int_t i=0;i<classes->GetEntriesFast();i++){
      TString className = classes->At(i)->GetName();
      if      (className.Contains("-A-"))         continue;
      else if (className.Contains("-C-"))         continue;
      else if (className.Contains("-E-"))         continue;
      else if (className.Contains("-AC-"))        continue;
      else if (className.Contains("-ACE-"))       continue;
      else if (className.Contains("-GA-"))        continue;
      else if (className.Contains("-GC-"))        continue;
      else if (className.Contains("1A-ABCE-"))    continue;
      else if (className.Contains("1C-ABCE-"))    continue;
      else if (className.Contains("C0LSR-ABCE-")) continue;
      hClassL0BvsRun     ->Fill(srun,className.Data(),Double_t(class_l0b[i]));
      hClassL2AvsRun     ->Fill(srun,className.Data(),Double_t(class_l2a[i]));
      hClassLifetimeVsRun->Fill(srun,className.Data(),class_lifetime[i]);
      hClassLumiVsRun    ->Fill(srun,className.Data(),class_lumi[i]);
    }
    for (Int_t ibit=0;ibit<NBITS;ibit++) {
      if (alias_reconstructed[ibit]<1) continue;
      hRecorded         ->Fill(srun,bitNames[ibit],alias_recorded[ibit]);
      hReconstructed    ->Fill(srun,bitNames[ibit],alias_reconstructed[ibit]);
      hAccepted         ->Fill(srun,bitNames[ibit],alias_accepted[ibit]);
      hAccStep1         ->Fill(srun,bitNames[ibit],alias_acc_step1[ibit]);
      hAccStep2         ->Fill(srun,bitNames[ibit],alias_acc_step2[ibit]);
      hAccStep3         ->Fill(srun,bitNames[ibit],alias_acc_step3[ibit]);
      hAccStep4         ->Fill(srun,bitNames[ibit],alias_acc_step4[ibit]);
      hAccStep5         ->Fill(srun,bitNames[ibit],alias_acc_step5[ibit]);
      hAccStep6         ->Fill(srun,bitNames[ibit],alias_acc_step6[ibit]);
      hAccStep7         ->Fill(srun,bitNames[ibit],alias_acc_step7[ibit]);
      hAccStep8         ->Fill(srun,bitNames[ibit],alias_acc_step8[ibit]);
      hAccStep9         ->Fill(srun,bitNames[ibit],alias_acc_step9[ibit]);
      hLumiRecorded     ->Fill(srun,bitNames[ibit],alias_lumi_recorded[ibit]);
      hLumiReconstructed->Fill(srun,bitNames[ibit],alias_lumi_reconstructed[ibit]);
      hLumiAccepted     ->Fill(srun,bitNames[ibit],alias_lumi_accepted[ibit]);
    }
  }
  hInteractionRate   ->LabelsDeflate("X");
  hMu                ->LabelsDeflate("X");
  hBCs               ->LabelsDeflate("X");
  hDuration          ->LabelsDeflate("X");
  hLumiSeen          ->LabelsDeflate("X");
  hActiveDetectors   ->LabelsDeflate("X");
  hMeanV0MOn         ->LabelsDeflate("X");
  hMeanV0MOf         ->LabelsDeflate("X");
  hMeanOFO           ->LabelsDeflate("X");
  hMeanTKL           ->LabelsDeflate("X");
  hMeanV0MOnHM       ->LabelsDeflate("X");
  hMeanV0MOfHM       ->LabelsDeflate("X");
  hMeanOFOHM         ->LabelsDeflate("X");
  hMeanTKLHM         ->LabelsDeflate("X");
  hV0MOnVsOfA        ->LabelsDeflate("X");
  hV0MOnVsOfB        ->LabelsDeflate("X");
  hSPDOnVsOfA        ->LabelsDeflate("X");
  hSPDOnVsOfB        ->LabelsDeflate("X");
  hV0MOnVsOfA_OADB   ->LabelsDeflate("X");
  hV0MOnVsOfB_OADB   ->LabelsDeflate("X");
  hSPDOnVsOfA_OADB   ->LabelsDeflate("X");
  hSPDOnVsOfB_OADB   ->LabelsDeflate("X");
  hThresholdV0M      ->LabelsDeflate("X");
  hThresholdV0Mnorm  ->LabelsDeflate("X");
  hClassL0BvsRun     ->LabelsDeflate("XY");
  hClassL2AvsRun     ->LabelsDeflate("XY");
  hClassLifetimeVsRun->LabelsDeflate("XY");
  hClassLumiVsRun    ->LabelsDeflate("XY");
  hRecorded          ->LabelsDeflate("XY");
  hReconstructed     ->LabelsDeflate("XY");
  hAccepted          ->LabelsDeflate("XY");
  hAccStep1          ->LabelsDeflate("XY");
  hAccStep2          ->LabelsDeflate("XY");
  hAccStep3          ->LabelsDeflate("XY");
  hAccStep4          ->LabelsDeflate("XY");
  hAccStep5          ->LabelsDeflate("XY");
  hAccStep6          ->LabelsDeflate("XY");
  hAccStep7          ->LabelsDeflate("XY");
  hAccStep8          ->LabelsDeflate("XY");
  hAccStep9          ->LabelsDeflate("XY");
  hLumiRecorded      ->LabelsDeflate("XY");
  hLumiReconstructed ->LabelsDeflate("XY");
  hLumiAccepted      ->LabelsDeflate("XY");
  
  SetHisto(hActiveDetectors,1);
  SetHisto(hInteractionRate,1);
  SetHisto(hMu);
  SetHisto(hBCs,1);
  SetHisto(hDuration,1);
  SetHisto(hLumiSeen,1);
  SetHisto(hMeanV0MOn);
  SetHisto(hMeanV0MOf);
  SetHisto(hMeanOFO);
  SetHisto(hMeanTKL);
  SetHisto(hMeanV0MOnHM);
  SetHisto(hMeanV0MOfHM);
  SetHisto(hMeanOFOHM);
  SetHisto(hMeanTKLHM);
  SetHisto(hV0MOnVsOfA);
  SetHisto(hV0MOnVsOfB);
  SetHisto(hSPDOnVsOfA);
  SetHisto(hSPDOnVsOfB);
  SetHisto(hV0MOnVsOfA_OADB);
  SetHisto(hV0MOnVsOfB_OADB);
  SetHisto(hSPDOnVsOfA_OADB);
  SetHisto(hSPDOnVsOfB_OADB);
  SetHisto(hThresholdV0M);
  SetHisto(hThresholdV0Mnorm);
  SetHisto(hClassL0BvsRun);
  SetHisto(hClassL2AvsRun);
  SetHisto(hClassLifetimeVsRun);
  SetHisto(hClassLumiVsRun);
  SetHisto(hRecorded);
  SetHisto(hReconstructed);
  SetHisto(hAccepted);
  SetHisto(hAccStep1);
  SetHisto(hAccStep2);
  SetHisto(hAccStep3);
  SetHisto(hAccStep4);
  SetHisto(hAccStep5);
  SetHisto(hAccStep6);
  SetHisto(hAccStep7);
  SetHisto(hAccStep8);
  SetHisto(hAccStep9);
  SetHisto(hLumiRecorded);
  SetHisto(hLumiReconstructed);
  SetHisto(hLumiAccepted);

  TCanvas* cActiveDetectors = new TCanvas("active_detectors","Active Detectors",1800,500);
  cActiveDetectors->SetMargin(0.05,0.01,0.18,0.06);
  hActiveDetectors->GetYaxis()->SetLabelOffset(0.001);
  hActiveDetectors->SetMaximum(2);
  hActiveDetectors->Draw("col");
  AddFillSeparationLines(hActiveDetectors,fills);
  AddPeriodSeparationLines(hActiveDetectors,periods);
  gPad->Print("detectors.png");
  gPad->Print("global_properties.pdf(");

  TCanvas* cInteractionRate = new TCanvas("cInteractionRate","Interaction Rate",1800,500);
  cInteractionRate->SetMargin(0.05,0.01,0.18,0.06);
  hInteractionRate->Draw();
  AddFillSeparationLines(hInteractionRate,fills);
  AddPeriodSeparationLines(hInteractionRate,periods);
  gPad->Print("rate.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMu = new TCanvas("mu","mu",1800,500);
  cMu->SetMargin(0.05,0.01,0.18,0.06);
  hMu->Draw("h");
  AddFillSeparationLines(hMu,fills);
  AddPeriodSeparationLines(hMu,periods);
  gPad->Print("mu.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cBCs = new TCanvas("bcs","bcs",1800,500);
  cBCs->SetMargin(0.05,0.01,0.18,0.06);
  hBCs->Draw("h");
  AddFillSeparationLines(hBCs,fills);
  AddPeriodSeparationLines(hBCs,periods);
  gPad->Print("bcs.png");
  gPad->Print("global_properties.pdf");
  
  TCanvas* cDuration = new TCanvas("duration","duration",1800,500);
  cDuration->SetMargin(0.05,0.01,0.18,0.06);
  hDuration->SetTitle(Form("Duration in seconds: total= %.0f s = %.0f h",hDuration->Integral(),hDuration->Integral()/3600));
  hDuration->SetFillColor(kBlue);
  hDuration->Draw("h");
  AddFillSeparationLines(hDuration,fills);
  AddPeriodSeparationLines(hDuration,periods);
  gPad->Print("duration.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cLumiSeen = new TCanvas("lumiseen","lumi seen",1800,500);
  cLumiSeen->SetMargin(0.05,0.01,0.18,0.06);
  hLumiSeen->SetTitle(Form("Luminosity seen [1/ub]: total= %.3f",hLumiSeen->Integral()));
  hLumiSeen->SetFillColor(kBlue);
  hLumiSeen->Draw("h");
  AddFillSeparationLines(hLumiSeen,fills);
  AddPeriodSeparationLines(hLumiSeen,periods);
  gPad->Print("lumi_seen.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanV0MOn = new TCanvas("meanV0MOn","meanV0MOn",1800,500);
  cMeanV0MOn->SetMargin(0.05,0.01,0.18,0.06);
  hMeanV0MOn->SetFillColor(0);
  hMeanV0MOn->Draw("");
  AddFillSeparationLines(hMeanV0MOn,fills);
  AddPeriodSeparationLines(hMeanV0MOn,periods);
  gPad->Print("meanV0MOn.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanV0MOf = new TCanvas("meanV0MOf","meanV0MOf",1800,500);
  cMeanV0MOf->SetMargin(0.05,0.01,0.18,0.06);
  hMeanV0MOf->Draw("");
  AddFillSeparationLines(hMeanV0MOf,fills);
  AddPeriodSeparationLines(hMeanV0MOf,periods);
  gPad->Print("meanV0MOf.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanOFO = new TCanvas("meanOFO","meanOFO",1800,500);
  cMeanOFO->SetMargin(0.05,0.01,0.18,0.06);
  hMeanOFO->Draw("");
  AddFillSeparationLines(hMeanOFO,fills);
  AddPeriodSeparationLines(hMeanOFO,periods);
  gPad->Print("meanOFO.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanTKL = new TCanvas("meanTKL","meanTKL",1800,500);
  cMeanTKL->SetMargin(0.05,0.01,0.18,0.06);
  hMeanTKL->Draw("");
  AddFillSeparationLines(hMeanTKL,fills);
  AddPeriodSeparationLines(hMeanTKL,periods);
  gPad->Print("meanTKL.png");
  gPad->Print("global_properties.pdf");
  
  TCanvas* cMeanV0MOnHM = new TCanvas("meanV0MOnHM","meanV0MOnHM",1800,500);
  cMeanV0MOnHM->SetMargin(0.05,0.01,0.18,0.06);
  hMeanV0MOnHM->SetFillColor(0);
  hMeanV0MOnHM->Draw("");
  AddFillSeparationLines(hMeanV0MOnHM,fills);
  AddPeriodSeparationLines(hMeanV0MOnHM,periods);
  gPad->Print("meanV0MOnHM.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanV0MOfHM = new TCanvas("meanV0MOfHM","meanV0MOfHM",1800,500);
  cMeanV0MOfHM->SetMargin(0.05,0.01,0.18,0.06);
  hMeanV0MOfHM->Draw("");
  AddFillSeparationLines(hMeanV0MOfHM,fills);
  AddPeriodSeparationLines(hMeanV0MOfHM,periods);
  gPad->Print("meanV0MOfHM.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanOFOHM = new TCanvas("meanOFOHM","meanOFOHM",1800,500);
  cMeanOFOHM->SetMargin(0.05,0.01,0.18,0.06);
  hMeanOFOHM->Draw("");
  AddFillSeparationLines(hMeanOFOHM,fills);
  AddPeriodSeparationLines(hMeanOFOHM,periods);
  gPad->Print("meanOFOHM.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cMeanTKLHM = new TCanvas("meanTKLHM","meanTKLHM",1800,500);
  cMeanTKLHM->SetMargin(0.05,0.01,0.18,0.06);
  hMeanTKLHM->Draw("");
  AddFillSeparationLines(hMeanTKLHM,fills);
  AddPeriodSeparationLines(hMeanTKLHM,periods);
  gPad->Print("meanTKLHM.png");
  gPad->Print("global_properties.pdf");
  
  TCanvas* cThresholdV0M = new TCanvas("thresholdV0M","thresholdV0M",1800,500);
  cThresholdV0M->SetMargin(0.05,0.01,0.18,0.06);
  hThresholdV0M->Draw("");
  AddFillSeparationLines(hThresholdV0M,fills);
  AddPeriodSeparationLines(hThresholdV0M,periods);
  gPad->Print("thresholdV0M.png");
  gPad->Print("global_properties.pdf");

  TCanvas* cThresholdV0Mnorm = new TCanvas("thresholdV0Mnorm","thresholdV0Mnorm",1800,500);
  cThresholdV0Mnorm->SetMargin(0.05,0.01,0.18,0.06);
  hThresholdV0Mnorm->Draw("e");
  AddFillSeparationLines(hThresholdV0Mnorm,fills);
  AddPeriodSeparationLines(hThresholdV0Mnorm,periods);
  gPad->Print("thresholdV0Mnorm.png");
  
  gPad->Print("global_properties.pdf)");
  
  TCanvas* cV0MOnVsOfA = new TCanvas("aV0MOnVsOf","aV0MOnVsOf",1800,500);
  cV0MOnVsOfA->SetMargin(0.05,0.01,0.18,0.06);
  hV0MOnVsOfA_OADB->SetLineColor(kRed);
  Int_t maxBin1,minBin1,maxBin2,minBin2;
  Double_t ymin,ymax;
  maxBin1 = hV0MOnVsOfA_OADB->GetMaximumBin();
  minBin1 = hV0MOnVsOfA_OADB->GetMinimumBin();
  maxBin2 = hV0MOnVsOfA->GetMaximumBin();
  minBin2 = hV0MOnVsOfA->GetMinimumBin();
  ymin = TMath::Min(hV0MOnVsOfA_OADB->GetBinContent(minBin1),hV0MOnVsOfA->GetBinContent(minBin2));
  ymax = TMath::Max(hV0MOnVsOfA_OADB->GetBinContent(maxBin1),hV0MOnVsOfA->GetBinContent(maxBin2));
  hV0MOnVsOfA_OADB->SetMinimum(ymin-(ymax-ymin)*0.1);
  hV0MOnVsOfA_OADB->SetMaximum(ymax+(ymax-ymin)*0.1);
  hV0MOnVsOfA_OADB->Draw("");
  hV0MOnVsOfA->Draw("e same");
  AddFillSeparationLines(hV0MOnVsOfA,fills);
  AddPeriodSeparationLines(hV0MOnVsOfA,periods);
  gPad->Print("V0MOnVsOfA.png");

  TCanvas* cV0MOnVsOfB = new TCanvas("V0MOnVsOfB","V0MOnVsOfB",1800,500);
  cV0MOnVsOfB->SetMargin(0.05,0.01,0.18,0.06);
  hV0MOnVsOfB_OADB->SetLineColor(kRed);
  maxBin1 = hV0MOnVsOfB_OADB->GetMaximumBin();
  minBin1 = hV0MOnVsOfB_OADB->GetMinimumBin();
  maxBin2 = hV0MOnVsOfB->GetMaximumBin();
  minBin2 = hV0MOnVsOfB->GetMinimumBin();
  ymin = TMath::Min(hV0MOnVsOfB_OADB->GetBinContent(minBin1),hV0MOnVsOfB->GetBinContent(minBin2));
  ymax = TMath::Max(hV0MOnVsOfB_OADB->GetBinContent(maxBin1),hV0MOnVsOfB->GetBinContent(maxBin2));
  hV0MOnVsOfB_OADB->SetMinimum(ymin-(ymax-ymin)*0.1);
  hV0MOnVsOfB_OADB->SetMaximum(ymax+(ymax-ymin)*0.1);
  hV0MOnVsOfB_OADB->Draw("");
  hV0MOnVsOfB->Draw("e same");
  AddFillSeparationLines(hV0MOnVsOfB,fills);
  AddPeriodSeparationLines(hV0MOnVsOfB,periods);
  gPad->Print("V0MOnVsOfB.png");

  TCanvas* cSPDOnVsOfA = new TCanvas("aSPDOnVsOf","aSPDOnVsOf",1800,500);
  cSPDOnVsOfA->SetMargin(0.05,0.01,0.18,0.06);
  hSPDOnVsOfA_OADB->SetLineColor(kRed);
  maxBin1 = hSPDOnVsOfA_OADB->GetMaximumBin();
  minBin1 = hSPDOnVsOfA_OADB->GetMinimumBin();
  maxBin2 = hSPDOnVsOfA->GetMaximumBin();
  minBin2 = hSPDOnVsOfA->GetMinimumBin();
  ymin = TMath::Min(hSPDOnVsOfA_OADB->GetBinContent(minBin1),hSPDOnVsOfA->GetBinContent(minBin2));
  ymax = TMath::Max(hSPDOnVsOfA_OADB->GetBinContent(maxBin1),hSPDOnVsOfA->GetBinContent(maxBin2));
  hSPDOnVsOfA_OADB->SetMinimum(ymin-(ymax-ymin)*0.1);
  hSPDOnVsOfA_OADB->SetMaximum(ymax+(ymax-ymin)*0.1);

  hSPDOnVsOfA_OADB->Draw("");
  hSPDOnVsOfA->Draw("e same");
  AddFillSeparationLines(hSPDOnVsOfA,fills);
  AddPeriodSeparationLines(hSPDOnVsOfA,periods);
  gPad->Print("SPDOnVsOfA.png");

  TCanvas* cSPDOnVsOfB = new TCanvas("SPDOnVsOfB","SPDOnVsOfB",1800,500);
  cSPDOnVsOfB->SetMargin(0.05,0.01,0.18,0.06);
  hSPDOnVsOfB_OADB->SetLineColor(kRed);
  maxBin1 = hSPDOnVsOfB_OADB->GetMaximumBin();
  minBin1 = hSPDOnVsOfB_OADB->GetMinimumBin();
  maxBin2 = hSPDOnVsOfB->GetMaximumBin();
  minBin2 = hSPDOnVsOfB->GetMinimumBin();
  ymin = TMath::Min(hSPDOnVsOfB_OADB->GetBinContent(minBin1),hSPDOnVsOfB->GetBinContent(minBin2));
  ymax = TMath::Max(hSPDOnVsOfB_OADB->GetBinContent(maxBin1),hSPDOnVsOfB->GetBinContent(maxBin2));
  hSPDOnVsOfB_OADB->SetMinimum(ymin-(ymax-ymin)*0.1);
  hSPDOnVsOfB_OADB->SetMaximum(ymax+(ymax-ymin)*0.1);

  hSPDOnVsOfB_OADB->Draw("");
  hSPDOnVsOfB->Draw("e same");
  AddFillSeparationLines(hSPDOnVsOfB,fills);
  AddPeriodSeparationLines(hSPDOnVsOfB,periods);
  gPad->Print("SPDOnVsOfB.png");

  
  TFile* fglobal = new TFile("global_properties.root","recreate");
  hActiveDetectors->Write();
  hInteractionRate->Write();
  hMu->Write();
  hBCs->Write();
  hDuration->Write();
  hLumiSeen->Write();
  hMeanTKL->Write();
  hMeanTKLHM->Write();
  hMeanV0MOf->Write();
  hMeanV0MOn->Write();
  fglobal->Close();

  TFile* fclassL0B = new TFile("class_L0B_counts.root","recreate");
  TCanvas* cClassL0B = new TCanvas("cClassL0B","Class L0B vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  hClassL0BvsRun->Draw("col");
  gPad->Print("class_L0B_counts.pdf(");
  hClassL0BvsRun->Write();
  TCanvas* cL0B = new TCanvas("cL0B","L0B vs run",1800,500);
  cL0B->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassL0BvsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassL0BvsRun->ProjectionX(Form("hClassL0BvsRun_%02i",i),i,i);
    h->SetTitle(Form("%s L0B counts: %.0f",hClassL0BvsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
    SetHisto(h);
    h->Draw();
    AddFillSeparationLines(h,fills);
    AddPeriodSeparationLines(h,periods);
    gPad->Print("class_L0B_counts.pdf");
    h->Write();
  }
  gPad->Print("class_L0B_counts.pdf]");
  fclassL0B->Close();
  
  TFile* fclassL2A = new TFile("class_L2A_counts.root","recreate");
  TCanvas* cClassL2A = new TCanvas("cClassL2A","Class L2A vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  hClassL2AvsRun->Draw("col");
  gPad->Print("class_L2A_counts.pdf(");
  hClassL2AvsRun->Write();
  TCanvas* cL2A = new TCanvas("cL2A","L2A vs run",1800,500);
  cL2A->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassL2AvsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassL2AvsRun->ProjectionX(Form("hClassL2AvsRun_%02i",i),i,i);
    h->SetTitle(Form("%s L2A counts: %.0f",hClassL2AvsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
    SetHisto(h);
    h->Draw();
    AddFillSeparationLines(h,fills);
    AddPeriodSeparationLines(h,periods);
    gPad->Print("class_L2A_counts.pdf");
    h->Write();
  }
  gPad->Print("class_L2A_counts.pdf]");
  fclassL2A->Close();
  
  TFile* fclassLifetime = new TFile("class_lifetime.root","recreate");
  TCanvas* cClassLifetime = new TCanvas("cClassLifetime","Lifetime class-by-class vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  hClassLifetimeVsRun->Draw("col");
  gPad->Print("class_lifetime.pdf(");
  hClassLifetimeVsRun->Write();
  TCanvas* cLifetime = new TCanvas("cLifetime","Lifetime vs run",1800,500);
  gPad->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassLifetimeVsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassLifetimeVsRun->ProjectionX(Form("hClassLifetimeVsRun_%02i",i),i,i);
    h->SetTitle(Form("%s lifetime",hClassLifetimeVsRun->GetYaxis()->GetBinLabel(i)));
    SetHisto(h);
    h->Draw("");
    AddFillSeparationLines(h,fills);
    AddPeriodSeparationLines(h,periods);
    gPad->Print("class_lifetime.pdf");
    h->Write();
  }
  gPad->Print("class_lifetime.pdf]");
  fclassLifetime->Close();
  
  TFile* fclassLumi = new TFile("class_lumi.root","recreate");
  TCanvas* cClassLumi = new TCanvas("cClassLumi","Luminosity class-by-class vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  hClassLumiVsRun->Draw("col");
  gPad->Print("class_lumi.pdf(");

  TCanvas* clumi = new TCanvas("clumi","lumi vs run",1800,500);
  gPad->SetMargin(0.05,0.01,0.18,0.06);

  for (Int_t i=1;i<=hClassLumiVsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassLumiVsRun->ProjectionX(Form("hClassLumiVsRun_%s",hClassLumiVsRun->GetYaxis()->GetBinLabel(i)),i,i);
    h->SetTitle(Form("%s luminosity [ub-1]: %.0f",hClassLumiVsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
    SetHisto(h);
    h->Draw("");
    AddFillSeparationLines(h,fills);
    AddPeriodSeparationLines(h,periods);
    gPad->Print("class_lumi.pdf");
    h->Write();
  }
  gPad->Print("class_lumi.pdf]");
  fclassLumi->Close();

  TCanvas* dummy = new TCanvas("dummy","dummy",1800,500);
  gPad->SetMargin(0.05,0.01,0.18,0.06);
  gPad->Print("alias_event_statistics.pdf[");
  gPad->Print("alias_lumi_statistics.pdf[");
  gPad->Print("accepted_fraction.pdf[");
  gPad->Print("rejected_fraction.pdf[");

  TCanvas* cCounts           = new TCanvas("c_alias_counts"   ,"c_alias_counts"   ,1800,500);
  TCanvas* cLumi             = new TCanvas("c_lumi"           ,"c_lumi"           ,1800,500);
  TCanvas* cAcceptedFraction = new TCanvas("accepted_fraction","accepted fraction",1800,500);
  TCanvas* cRejectedFraction = new TCanvas("rejected_fraction","rejected fraction",1800,500);
  cCounts->SetMargin(0.05,0.01,0.18,0.06);
  cLumi->SetMargin(0.05,0.01,0.18,0.06);
  cAcceptedFraction->SetMargin(0.05,0.01,0.18,0.06);
  cRejectedFraction->SetMargin(0.05,0.01,0.18,0.06);
  cRejectedFraction->SetLogy();

  TFile* fstat = new TFile("alias_statistics.root","recreate");
  for (Int_t i=1;i<=hRecorded->GetNbinsY();i++) {
    char* bitName = hRecorded->GetYaxis()->GetBinLabel(i);
    printf("bit=%i %s\n",i,bitName);
    TH1D* hRecorded1D = hRecorded->ProjectionX(Form("hRecorded_%s",bitName),i,i);
    TH1D* hReconstructed1D = hReconstructed->ProjectionX(Form("hReconstructed_%s",bitName),i,i);
    TH1D* hAccepted1D  = hAccepted->ProjectionX(Form("hAccepted_%s",bitName),i,i);
    TH1D* hAccStep1_1D = hAccStep1->ProjectionX(Form("hAccStep1_%s",bitName),i,i);
    TH1D* hAccStep2_1D = hAccStep2->ProjectionX(Form("hAccStep2_%s",bitName),i,i);
    TH1D* hAccStep3_1D = hAccStep3->ProjectionX(Form("hAccStep3_%s",bitName),i,i);
    TH1D* hAccStep4_1D = hAccStep4->ProjectionX(Form("hAccStep4_%s",bitName),i,i);
    TH1D* hAccStep5_1D = hAccStep5->ProjectionX(Form("hAccStep5_%s",bitName),i,i);
    TH1D* hAccStep6_1D = hAccStep6->ProjectionX(Form("hAccStep6_%s",bitName),i,i);
    TH1D* hAccStep7_1D = hAccStep7->ProjectionX(Form("hAccStep7_%s",bitName),i,i);
    TH1D* hAccStep8_1D = hAccStep8->ProjectionX(Form("hAccStep8_%s",bitName),i,i);
    TH1D* hAccStep9_1D = hAccStep9->ProjectionX(Form("hAccStep9_%s",bitName),i,i);
    
    TH1D* hRejected1D = hReconstructed->ProjectionX(Form("hRejected1D_%s",bitName),i,i);
    TH1D* hLumiRecorded1D = hLumiRecorded->ProjectionX(Form("hLumiRecorded_%s",bitName),i,i);
    TH1D* hLumiReconstructed1D = hLumiReconstructed->ProjectionX(Form("hLumiReconstructed_%s",bitName),i,i);
    TH1D* hLumiAccepted1D = hLumiAccepted->ProjectionX(Form("hLumiAccepted_%s",bitName),i,i);
    hRejected1D->Add(hAccepted1D,-1);
    
    if (hReconstructed1D->Integral()<1) continue;
    SetHisto(hRecorded1D);
    SetHisto(hReconstructed1D);
    SetHisto(hAccepted1D);
    SetHisto(hAccStep1_1D);
    SetHisto(hAccStep2_1D);
    SetHisto(hAccStep3_1D);
    SetHisto(hAccStep4_1D);
    SetHisto(hAccStep5_1D);
    SetHisto(hAccStep6_1D);
    SetHisto(hAccStep7_1D);
    SetHisto(hAccStep8_1D);
    SetHisto(hAccStep9_1D);
    SetHisto(hRejected1D);
    SetHisto(hLumiRecorded1D);
    SetHisto(hLumiReconstructed1D);
    SetHisto(hLumiAccepted1D);
    hRecorded1D->Sumw2();
    hReconstructed1D->Sumw2();
    hAccepted1D->Sumw2();
//    hRejected1D->Sumw2();
    hRecorded1D->SetLineColor(kRed+2);
    hRecorded1D->SetFillColor(kRed+2);
    hReconstructed1D->SetLineColor(kBlue);
    hReconstructed1D->SetFillColor(0);
    hAccepted1D->SetLineColor(kGreen+2);
    hAccepted1D->SetFillColor(kGreen+2);
    hLumiRecorded1D->SetLineColor(kRed+2);
    hLumiRecorded1D->SetFillColor(kRed+2);
    hLumiReconstructed1D->SetLineColor(kBlue);
    hLumiReconstructed1D->SetFillColor(kBlue);
    hLumiAccepted1D->SetLineColor(kGreen+2);
    hLumiAccepted1D->SetFillColor(kGreen+2);

    TH1D* hAcceptedFraction = (TH1D*) hReconstructed1D->Clone(Form("hAcceptedFraction_%s",bitName));
    TH1D* hAccStep1Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep1Fraction_%s",bitName));
    TH1D* hAccStep2Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep2Fraction_%s",bitName));
    TH1D* hAccStep3Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep3Fraction_%s",bitName));
    TH1D* hAccStep4Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep4Fraction_%s",bitName));
    TH1D* hAccStep5Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep5Fraction_%s",bitName));
    TH1D* hAccStep6Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep61raction_%s",bitName));
    TH1D* hAccStep7Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep71raction_%s",bitName));
    TH1D* hAccStep8Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep81raction_%s",bitName));
    TH1D* hAccStep9Fraction = (TH1D*) hReconstructed1D->Clone(Form("hAccStep91raction_%s",bitName));
    hAcceptedFraction->SetTitle(Form("Accepted fraction: %s",bitName));
    hAcceptedFraction->Divide(hAccepted1D,hReconstructed1D,1,1,"B");
    hAccStep1Fraction->Divide(hAccStep1_1D,hReconstructed1D,1,1,"B");
    hAccStep2Fraction->Divide(hAccStep2_1D,hReconstructed1D,1,1,"B");
    hAccStep3Fraction->Divide(hAccStep3_1D,hReconstructed1D,1,1,"B");
    hAccStep4Fraction->Divide(hAccStep4_1D,hReconstructed1D,1,1,"B");
    hAccStep5Fraction->Divide(hAccStep5_1D,hReconstructed1D,1,1,"B");
    hAccStep6Fraction->Divide(hAccStep6_1D,hReconstructed1D,1,1,"B");
    hAccStep7Fraction->Divide(hAccStep7_1D,hReconstructed1D,1,1,"B");
    hAccStep8Fraction->Divide(hAccStep8_1D,hReconstructed1D,1,1,"B");
    hAccStep9Fraction->Divide(hAccStep9_1D,hReconstructed1D,1,1,"B");
    hAcceptedFraction->SetFillColor(0);
    hAcceptedFraction->SetLineWidth(2);
    TH1D* hRejectedFraction = (TH1D*) hReconstructed1D->Clone(Form("hRejectedFraction_%s",bitName));
    hRejectedFraction->SetTitle(Form("Rejected fraction: %s",bitName));
    hRejectedFraction->Divide(hRejected1D,hReconstructed1D,1,1,"B");
    hRejectedFraction->SetFillColor(0);
    hRejectedFraction->SetLineWidth(2);
    hAccStep1Fraction->SetLineColor(kBlue);
    hAccStep2Fraction->SetLineColor(kRed);
    hAccStep3Fraction->SetLineColor(kMagenta);
    hAccStep4Fraction->SetLineColor(kBlue+2);
    hAccStep5Fraction->SetLineColor(kRed+2);
    hAccStep6Fraction->SetLineColor(kMagenta+2);
    hAccStep7Fraction->SetLineColor(kGreen+2);
    hAccStep8Fraction->SetLineColor(kYellow+1);
    hAccStep9Fraction->SetLineColor(kGreen+1);
    cAcceptedFraction->cd();
    hAcceptedFraction->SetTitle(Form("%s: average accepted fraction = %.3f",bitName,hAccepted1D->Integral()/hReconstructed1D->Integral()));
    Float_t min[10];
    min[0] = hAcceptedFraction->GetBinContent(hAcceptedFraction->GetMinimumBin());
    min[1] = hAccStep1Fraction->GetBinContent(hAccStep1Fraction->GetMinimumBin());
    min[2] = hAccStep2Fraction->GetBinContent(hAccStep2Fraction->GetMinimumBin());
    min[3] = hAccStep3Fraction->GetBinContent(hAccStep3Fraction->GetMinimumBin());
    min[4] = hAccStep4Fraction->GetBinContent(hAccStep4Fraction->GetMinimumBin());
    min[5] = hAccStep5Fraction->GetBinContent(hAccStep5Fraction->GetMinimumBin());
    min[6] = hAccStep6Fraction->GetBinContent(hAccStep6Fraction->GetMinimumBin());
    min[7] = hAccStep7Fraction->GetBinContent(hAccStep7Fraction->GetMinimumBin());
    min[8] = hAccStep8Fraction->GetBinContent(hAccStep8Fraction->GetMinimumBin());
    min[9] = hAccStep9Fraction->GetBinContent(hAccStep9Fraction->GetMinimumBin());
    Float_t max[10];
    max[0] = hAcceptedFraction->GetBinContent(hAcceptedFraction->GetMaximumBin());
    max[1] = hAccStep1Fraction->GetBinContent(hAccStep1Fraction->GetMaximumBin());
    max[2] = hAccStep2Fraction->GetBinContent(hAccStep2Fraction->GetMaximumBin());
    max[3] = hAccStep3Fraction->GetBinContent(hAccStep3Fraction->GetMaximumBin());
    max[4] = hAccStep4Fraction->GetBinContent(hAccStep4Fraction->GetMaximumBin());
    max[5] = hAccStep5Fraction->GetBinContent(hAccStep5Fraction->GetMaximumBin());
    max[6] = hAccStep6Fraction->GetBinContent(hAccStep6Fraction->GetMaximumBin());
    max[7] = hAccStep7Fraction->GetBinContent(hAccStep7Fraction->GetMaximumBin());
    max[8] = hAccStep8Fraction->GetBinContent(hAccStep8Fraction->GetMaximumBin());
    max[9] = hAccStep9Fraction->GetBinContent(hAccStep9Fraction->GetMaximumBin());
    
    Float_t elmax = TMath::MaxElement(10,max);
    Float_t elmin = TMath::MinElement(10,min);
    
    hAcceptedFraction->SetMinimum(elmin-0.1*(elmax-elmin));
    hAcceptedFraction->SetMaximum(elmax+0.1*(elmax-elmin));
    hAccStep1Fraction->SetMinimum(elmin-0.1*(elmax-elmin));
    hAccStep1Fraction->SetMaximum(1.0);
    hAccStep1Fraction->SetTitle(hAcceptedFraction->GetTitle());
    // hAcceptedFraction->Draw();
    hAccStep1Fraction->Draw();
    hAccStep2Fraction->Draw("same");
    hAccStep3Fraction->Draw("same");
    hAccStep4Fraction->Draw("same");
    hAccStep5Fraction->Draw("same");
    hAccStep6Fraction->Draw("same");
    hAccStep7Fraction->Draw("same");
    hAccStep8Fraction->Draw("same");
    hAccStep9Fraction->Draw("same");
    AddFillSeparationLines(hAcceptedFraction,fills);
    AddPeriodSeparationLines(hAcceptedFraction,periods);

    TLegend* legAcc = new TLegend(0.08,0.22,0.25,0.6);
//    legAcc->AddEntry(hAcceptedFraction,"accepted");
    legAcc->AddEntry(hAccStep1Fraction,"V0A & V0C");
    legAcc->AddEntry(hAccStep2Fraction,"+ no ZDCBG");
    legAcc->AddEntry(hAccStep3Fraction,"+ no ClsVsTklBG");
    legAcc->AddEntry(hAccStep4Fraction,"+ no V0C012vsTklBG");
    legAcc->AddEntry(hAccStep5Fraction,"+ no V0MOnVsOfPileup");
    legAcc->AddEntry(hAccStep6Fraction,"+ no SPDOnVsOfPileup");
    legAcc->AddEntry(hAccStep7Fraction,"+ no SPDVtxPileup");
    legAcc->AddEntry(hAccStep8Fraction,"+ no V0PFPileup");
    legAcc->AddEntry(hAccStep9Fraction,"+ no V0Casym");
    legAcc->Draw();

//    alias_acc_step1[ibit]     = Int_t(hHistStat->GetBinContent(3,j));
//    alias_acc_step2[ibit]     = Int_t(hHistStat->GetBinContent(4,j));
//    alias_acc_step3[ibit]     = Int_t(hHistStat->GetBinContent(5,j));
//    alias_acc_step4[ibit]     = Int_t(hHistStat->GetBinContent(6,j));
//    alias_acc_step5[ibit]     = Int_t(hHistStat->GetBinContent(7,j));
//    alias_acc_step6[ibit]     = Int_t(hHistStat->GetBinContent(8,j));
//    alias_acc_step7[ibit]     = Int_t(hHistStat->GetBinContent(9,j));
//    fHistStat->Fill("all",trigger,all);
//    fHistStat->Fill("accepted",trigger,accepted);
//    fHistStat->Fill("V0A & V0C",trigger,v0and);
//    fHistStat->Fill("+ !SPDClsVsTrkBG",trigger,plusNoSPDClsVsTrkBG);
//    fHistStat->Fill("+ !V0C012vsTklBG",trigger,plusNoV0C012vsTklBG);
//    fHistStat->Fill("+ !V0MOnVsOfPileup",trigger,plusNoV0MOnVsOfPileup);
//    fHistStat->Fill("+ !SPDOnVsOfPileup",trigger,plusNoSPDOnVsOfPileup);
//    fHistStat->Fill("+ !SPDVtxPileup",trigger,plusNoSPDVtxPileup);
//    fHistStat->Fill("+ !V0PFPileup",trigger,plusNoV0PFPileup);
//    fHistStat->Fill("+ !V0Casym",trigger,plusNoV0Casym);

    gPad->Print("accepted_fraction.pdf");
    gPad->Print(Form("accepted_fraction_%s.png",bitName));
    
    cRejectedFraction->cd();
    hRejectedFraction->SetMaximum(1);
    hRejectedFraction->SetTitle(Form("%s: average rejected fraction = %.3f",bitName,hRejected1D->Integral()/hReconstructed1D->Integral()));
    hRejectedFraction->Draw();
    AddFillSeparationLines(hRejectedFraction,fills);
    AddPeriodSeparationLines(hRejectedFraction,periods);
    gPad->Print("rejected_fraction.pdf");

    hRecorded1D->Write();
    hReconstructed1D->Write();
    hAccepted1D->Write();
    hAcceptedFraction->Write();
    hAccStep1Fraction->Write();
    hAccStep2Fraction->Write();
    hAccStep3Fraction->Write();
    hAccStep4Fraction->Write();
    hAccStep5Fraction->Write();
    hAccStep6Fraction->Write();
    hAccStep7Fraction->Write();
    hAccStep8Fraction->Write();
    hAccStep9Fraction->Write();
    hRejectedFraction->Write();

    
    if (hRecorded1D->Integral()<1) continue;
    
    cCounts->cd();
    hRecorded1D->SetTitle(Form("%s trigger counts: recorded=%0.f, total=%.0f, accepted=%.0f",bitName,hRecorded1D->Integral(),hReconstructed1D->Integral(),hAccepted1D->Integral()));
    hReconstructed1D->SetTitle(Form("%s trigger counts: recorded=%0.f, total=%.0f, accepted=%.0f",bitName,hRecorded1D->Integral(),hReconstructed1D->Integral(),hAccepted1D->Integral()));
    hAccepted1D->SetTitle(Form("%s trigger counts: recorded=%0.f, total=%.0f, accepted=%.0f",bitName,hRecorded1D->Integral(),hReconstructed1D->Integral(),hAccepted1D->Integral()));
    hRecorded1D->Draw("h");
    hReconstructed1D->Draw("h same");
    hAccepted1D->Draw("h same");
    AddFillSeparationLines(hAccepted1D,fills);
    AddPeriodSeparationLines(hAccepted1D,periods);

    gPad->RedrawAxis();
    gPad->Print("alias_event_statistics.pdf");

    cLumi->cd();
    hLumiRecorded1D->SetTitle(Form("%s luminosity [ub-1]: recorded=%.0g, total=%.0g, accepted=%.0g",bitName,hLumiRecorded1D->Integral(),hLumiReconstructed1D->Integral(),hLumiAccepted1D->Integral()));
    hLumiReconstructed1D->SetTitle(Form("%s luminosity [ub-1]: recorded=%.0g, total=%.0g, accepted=%.0g",bitName,hLumiRecorded1D->Integral(),hLumiReconstructed1D->Integral(),hLumiAccepted1D->Integral()));
    hLumiAccepted1D->SetTitle(Form("%s luminosity [ub-1]: recorded=%.0g, total=%.0g, accepted=%.0g",bitName,hLumiRecorded1D->Integral(),hLumiReconstructed1D->Integral(),hLumiAccepted1D->Integral()));
    hLumiRecorded1D->Draw("h");
    hLumiReconstructed1D->Draw("h same");
    hLumiAccepted1D->Draw("h same");
    AddFillSeparationLines(hLumiAccepted1D,fills);
    AddPeriodSeparationLines(hLumiAccepted1D,periods);
    gPad->RedrawAxis();
    gPad->Print("alias_lumi_statistics.pdf");
  }
  dummy->Print("alias_event_statistics.pdf]");
  dummy->Print("alias_lumi_statistics.pdf]");
  dummy->Print("accepted_fraction.pdf]");
  dummy->Print("rejected_fraction.pdf]");
  fstat->Close();
}

void SetHisto(TH1D* h, Bool_t setMinimumToZero){
  h->SetTitleFont(43);
  h->SetTitleSize(25);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetLabelSize(25);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTitleOffset(0.6);
  h->GetYaxis()->SetDecimals(1);
  h->SetLineColor(kBlue);
  h->SetFillColor(0);
  h->LabelsOption("av");
  h->SetLineWidth(2);
  if (setMinimumToZero) h->SetMinimum(0);
}

void SetHisto(TH2D* h, Bool_t setMinimumToZero){
  h->SetTitleFont(43);
  h->SetTitleSize(15);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetZaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetLabelSize(15);
  h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetDecimals(1);
  h->GetXaxis()->LabelsOption("av");
//  h->GetYaxis()->LabelsOption("a");
  h->SetLineWidth(2);
  h->SetLineColor(kBlue);
  if (setMinimumToZero) h->SetMinimum(0);
}

void AddFillSeparationLines(TH1* h, map<Int_t,Int_t> &fills){
//  return;
  gPad->Update();
  Double_t ymin = gPad->GetUymin();
  Double_t ymax = gPad->GetUymax();
  TLine * fillSeparationLine = new TLine();
  fillSeparationLine->SetLineColor(kRed);
  fillSeparationLine->SetLineWidth(1);
  for(Int_t iBin = 1; iBin < h->GetXaxis()->GetNbins(); iBin++) {
    UInt_t runOld = atoi(h->GetXaxis()->GetBinLabel(iBin));
    UInt_t runNew = atoi(h->GetXaxis()->GetBinLabel(iBin + 1));
    if (fills[runOld]==fills[runNew]) continue;
    fillSeparationLine->DrawLine(iBin,ymin,iBin,ymax);
  }
}

void AddPeriodSeparationLines(TH1* h,  map<Int_t,TString> &periods){
  gPad->Update();
  Double_t ymin = gPad->GetUymin();
  Double_t ymax = gPad->GetUymax();
  TLine * periodSeparationLine = new TLine();
  periodSeparationLine->SetLineColor(kMagenta);
  periodSeparationLine->SetLineWidth(3);
  for(Int_t iBin = 1; iBin < h->GetXaxis()->GetNbins(); iBin++) {
    UInt_t runOld = atoi(h->GetXaxis()->GetBinLabel(iBin));
    UInt_t runNew = atoi(h->GetXaxis()->GetBinLabel(iBin + 1));
    if (periods[runOld].EqualTo(periods[runNew])) continue;
    periodSeparationLine->DrawLine(iBin,ymin,iBin,ymax);
  }
}

