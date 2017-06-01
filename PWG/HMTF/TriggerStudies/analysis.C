#include "TPad.h"
#include "TFile.h"
#include "TStyle.h"
#include "TList.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TObjString.h"
#include "TBits.h"
#include "TLine.h"
#include "TChain.h"
#ifndef __CINT__
#include "AliTriggerBCMask.h"
#include "AliAnalysisTriggerScalers.h"
#endif
#include "TMath.h"

//enum triggerAliases {kINT=0,kSH1,kSH2,kVHM,kVHMNOPF,kVHMSPD1,kSH2NOPF,kSH2SPD1,kSH1NOPF,NTRIGGERS};
enum triggerAliases {kINT=0,kSH1,kSH2,kVHM,NTRIGGERS};

Bool_t IsOutOfBunchPileup(TBits* fIR1, TBits* fIR2, Int_t fRunNumber, UShort_t fBC);

void analysis(){
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetHistLineWidth(2);
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetHistLineColor(kBlue);
  
  TFile* fMean = new TFile("mean.root");
  TH1D* hMeanOfV0M = (TH1D*) fMean->Get("hMeanOfflineV0M");
  TH1D* hMeanOnV0M = (TH1D*) fMean->Get("hMeanOnlineV0M");
  TH1D* hMeanOfTKL = (TH1D*) fMean->Get("hMeanTKL");
  TH1D* hMeanOfCL1 = (TH1D*) fMean->Get("hMeanCL1");
  TH1D* hMeanOnFOR = (TH1D*) fMean->Get("hMeanOnlineFO");
  TH1D* hMeanOfFOR = (TH1D*) fMean->Get("hMeanOfflineFO");
  TChain* fTree = new TChain("events");
  fTree->Add("/alice/cernbox/hm/lhc15f/AnalysisResults.*.root");
  fTree->Add("/alice/cernbox/hm/lhc15h/AnalysisResults.*.root");
  fTree->Add("/alice/cernbox/hm/lhc15i/AnalysisResults.*.root");
  fTree->Add("/alice/cernbox/hm/lhc15j/AnalysisResults.*.root");
  fTree->Add("/alice/cernbox/hm/lhc15l/AnalysisResults.*.root");
//  fTree->Add("/alice/cernbox/hm/lhc15l/AnalysisResults.239319.root");
//  fTree->Add("/alice/cernbox/hm/lhc15l/AnalysisResults.241393.root");
//  fTree->Add("/alice/cernbox/hm/lhc15l/AnalysisResults.239319.root");
//  fTree->Add("/alice/cernbox/hm/lhc15l/AnalysisResults.*.root");
  
  TObjString* fClassesFired = new TObjString();
  Bool_t fIsIncomplete;
  Int_t fRunNumber;
  UShort_t fBC;
  Float_t fV0ATime;
  Float_t fV0CTime;
  Bool_t fBBFlag[64];
  Bool_t fBGFlag[64];
  Bool_t fIsFriend=0;
  ULong64_t fBBFlagPF[21];
  ULong64_t fBGFlagPF[21];
  Float_t fMTotV0A;
  Float_t fMTotV0C;
  Float_t fMRingV0A[4];
  Float_t fMRingV0C[4];
  UShort_t fTriggerChargeA;
  UShort_t fTriggerChargeC;
  TBits* fIR1 = new TBits();
  TBits* fIR2 = new TBits();
  TBits* fIRT0A = new TBits();
  TBits* fIRT0C = new TBits();
  TBits* fIRTVX = new TBits();
  UInt_t fNofTracklets;
  Int_t fNofITSClusters[6];
  TBits* fFOmap = new TBits();
  TBits* fFiredChipMap = new TBits();
  UInt_t fInputsL0;
  UInt_t fVertexContributors;
  UInt_t fPileupContributors;
  Bool_t fIsPileupSPD;
  Int_t fV0ADecision;
  Int_t fV0CDecision;
  Float_t fT0A[5];
  Float_t fT0C[5];
  Float_t fTVX[5];

  fTree->SetBranchAddress("fIsIncomplete",&fIsIncomplete);
  fTree->SetBranchAddress("fRunNumber",&fRunNumber);
  fTree->SetBranchAddress("fBC",&fBC);
  fTree->SetBranchAddress("fClassesFired",&fClassesFired);
  fTree->SetBranchAddress("fInputsL0",&fInputsL0);
  fTree->SetBranchAddress("fVertexContributors",&fVertexContributors);
  fTree->SetBranchAddress("fPileupContributors",&fPileupContributors);
  fTree->SetBranchAddress("fIsPileupSPD",&fIsPileupSPD);
  fTree->SetBranchAddress("fV0ATime",&fV0ATime);
  fTree->SetBranchAddress("fV0CTime",&fV0CTime);
  fTree->SetBranchAddress("fBBFlag",&fBBFlag);
  fTree->SetBranchAddress("fBGFlag",&fBGFlag);
  fTree->SetBranchAddress("fIsFriend",&fIsFriend);
  fTree->SetBranchAddress("fBBFlagPF",&fBBFlagPF);
  fTree->SetBranchAddress("fBGFlagPF",&fBGFlagPF);
  fTree->SetBranchAddress("fMTotV0A",&fMTotV0A);
  fTree->SetBranchAddress("fMTotV0C",&fMTotV0C);
  fTree->SetBranchAddress("fMRingV0A",&fMRingV0A);
  fTree->SetBranchAddress("fMRingV0C",&fMRingV0C);
  fTree->SetBranchAddress("fTriggerChargeA",&fTriggerChargeA);
  fTree->SetBranchAddress("fTriggerChargeC",&fTriggerChargeC);
  fTree->SetBranchAddress("fIR1",&fIR1);
  fTree->SetBranchAddress("fIR2",&fIR2);
  fTree->SetBranchAddress("fNofTracklets",&fNofTracklets);
  fTree->SetBranchAddress("fNofITSClusters",&fNofITSClusters);
  fTree->SetBranchAddress("fFOmap",&fFOmap);
  fTree->SetBranchAddress("fFiredChipMap",&fFiredChipMap);
  fTree->SetBranchAddress("fV0ADecision",&fV0ADecision);
  fTree->SetBranchAddress("fV0CDecision",&fV0CDecision);
  fTree->SetBranchAddress("fT0A",&fT0A);
  fTree->SetBranchAddress("fT0C",&fT0C);
  fTree->SetBranchAddress("fTVX",&fTVX);
  
  TH3D* hEventCounter = new TH3D("hEventCounter","",1,0,1,1,0,1,1,0,1);
  
  TH3D* hPileupBcRun = new TH3D("hPileupBcRun",";pileup;bc;",40,-20,20,4,0,4,1,0,1);

  TH3D* hClsVsTklall = new TH3D("hClsVsTklall",";;tracklets;clusters",NTRIGGERS,0,NTRIGGERS,400,0,400,400,0,4000);
  TH3D* hClsVsTklbgd = new TH3D("hClsVsTklbgd",";;tracklets;clusters",NTRIGGERS,0,NTRIGGERS,400,0,400,400,0,4000);
  TH3D* hClsVsTklsig = new TH3D("hClsVsTklsig",";;tracklets;clusters",NTRIGGERS,0,NTRIGGERS,400,0,400,400,0,4000);
  TH3D* hClsVsTklsel = new TH3D("hClsVsTklsel",";;tracklets;clusters",NTRIGGERS,0,NTRIGGERS,400,0,400,400,0,4000);
  
  TH3D* hBBAall = new TH3D("hBBAall",";BBA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBBCall = new TH3D("hBBCall",";BBC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGAall = new TH3D("hBGAall",";BGA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGCall = new TH3D("hBGCall",";BGC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBBAbgd = new TH3D("hBBAbgd",";BBA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBBCbgd = new TH3D("hBBCbgd",";BBC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGAbgd = new TH3D("hBGAbgd",";BGA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGCbgd = new TH3D("hBGCbgd",";BGC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBBAsig = new TH3D("hBBAsig",";BBA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBBCsig = new TH3D("hBBCsig",";BBC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGAsig = new TH3D("hBGAsig",";BGA flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  TH3D* hBGCsig = new TH3D("hBGCsig",";BGC flags;",NTRIGGERS,0,NTRIGGERS,33,0,33,1,0,1);
  
  TH3D* hV0ATime   = new TH3D("hV0ATime"  ,";;V0A time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hV0CTime   = new TH3D("hV0CTime  ",";;V0C time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hV0ATimeBB = new TH3D("hV0ATimeBB",";;V0A time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hV0CTimeBB = new TH3D("hV0CTimeBB",";;V0C time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hV0ATimeBG = new TH3D("hV0ATimeBG",";;V0A time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hV0CTimeBG = new TH3D("hV0CTimeBG",";;V0C time, ns;",NTRIGGERS,0,NTRIGGERS,180,-45,45,1,0,1);
  TH3D* hOnV0M     = new TH3D("hOnV0M","",NTRIGGERS,0,NTRIGGERS, 220,0,11000,1,0,1);
  TH3D* hOfV0M     = new TH3D("hOfV0M","",NTRIGGERS,0,NTRIGGERS, 200,0, 2000,1,0,1);
  TH3D* hOfTKL     = new TH3D("hOfTKL","",NTRIGGERS,0,NTRIGGERS, 200,0,  200,1,0,1);
  TH3D* hOfCL1     = new TH3D("hOfCL1","",NTRIGGERS,0,NTRIGGERS, 200,0,  400,1,0,1);
  TH3D* hOnFOR     = new TH3D("hOnFOR","",NTRIGGERS,0,NTRIGGERS, 200,0,  200,1,0,1);
  TH3D* hOfFOR     = new TH3D("hOfFOR","",NTRIGGERS,0,NTRIGGERS, 300,0,  300,1,0,1);
  TH3D* hOnV0Mnorm = new TH3D("hOnV0Mnorm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  TH3D* hOfV0Mnorm = new TH3D("hOfV0Mnorm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  TH3D* hOfTKLnorm = new TH3D("hOfTKLnorm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  TH3D* hOfCL1norm = new TH3D("hOfCL1norm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  TH3D* hOnFORnorm = new TH3D("hOnFORnorm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  TH3D* hOfFORnorm = new TH3D("hOfFORnorm","",NTRIGGERS,0,NTRIGGERS, 150,0,  15,1,0,1);
  
  TH3D* hOnFORvsOfFOR_all      = new TH3D("hOnFORvsOfFOR_all"     ,"",200,0,200,300,0,300,1,0,1);
  TH3D* hOnFORvsOfFOR_sel_mod0 = new TH3D("hOnFORvsOfFOR_sel_mod0","",200,0,200,300,0,300,1,0,1);
  TH3D* hOnFORvsOfFOR_sel_mod1 = new TH3D("hOnFORvsOfFOR_sel_mod1","",200,0,200,300,0,300,1,0,1);
  TH3D* hOnFORvsOfFOR_sel_mod2 = new TH3D("hOnFORvsOfFOR_sel_mod2","",200,0,200,300,0,300,1,0,1);
  TH3D* hOnFORvsOfFOR_sel_mod3 = new TH3D("hOnFORvsOfFOR_sel_mod3","",200,0,200,300,0,300,1,0,1);

  TH3D* hOnV0MvsOfV0M_all      = new TH3D("hOnV0MvsOfV0M_all"      ,"",200,0,2000,220,0,11000,1,0,1);
  TH3D* hOnV0MvsOfV0M_sel_mod0 = new TH3D("hOnV0MvsOfV0M_sel_mod0 ","",200,0,2000,220,0,11000,1,0,1);
  TH3D* hOnV0MvsOfV0M_sel_mod1 = new TH3D("hOnV0MvsOfV0M_sel_mod1 ","",200,0,2000,220,0,11000,1,0,1);
  TH3D* hOnV0MvsOfV0M_sel_mod2 = new TH3D("hOnV0MvsOfV0M_sel_mod2 ","",200,0,2000,220,0,11000,1,0,1);
  TH3D* hOnV0MvsOfV0M_sel_mod3 = new TH3D("hOnV0MvsOfV0M_sel_mod3 ","",200,0,2000,220,0,11000,1,0,1);

  TH2D* hV0C012vsV0C3 = new TH2D("hV0C012vsV0C3","",800,0,800,300,0,300);
 
  TH1D* h0TVXall = new TH1D("h0TVXall","",1,0,1);
  TH1D* h0TVXbad = new TH1D("h0TVXbad","",1,0,1);
  TH1D* h0VIRall = new TH1D("h0VIRall","",1,0,1);
  TH1D* h0VIRbad = new TH1D("h0VIRbad","",1,0,1);
  TH1D* h0TVXbcs = new TH1D("h0TVXbcs","",3564,0,3564);
  

  
  printf("Entries: %i\n",fTree->GetEntries());
  Int_t nSPD1  = 0;
  Int_t nWrong = 0;
  Bool_t tr[NTRIGGERS];
  for (Int_t ev=0;ev<fTree->GetEntries();ev++){
//  for (Int_t ev=0;ev<100000;ev++){
    fTree->GetEntry(ev);
    if (ev%100000==0) printf("Event=%10i Run=%i\n",ev, fRunNumber);
    
    const char* srun = Form("%i",fRunNumber);

    if (fClassesFired->String().Contains("C0TVX-B")){
      if (fIsIncomplete || fInputsL0==0) continue;
      if (!(fInputsL0 & 1<<2)) printf("0TVX input not fired");
      h0TVXall->Fill(srun,1);
      h0TVXbad->Fill(srun,0);
      if (!fIR2->TestBitNumber(90)) {
        h0TVXbad->Fill(srun,1);
        h0TVXbcs->Fill(fBC);
        printf("0TVX:");
        for (Int_t i=-90;i<=90;i++) printf("%i",fIR2->TestBitNumber(90+i));
        printf("\n");
      }
    }

    if (fClassesFired->String().Contains("CINT7-B")){
      if (fIsIncomplete || fInputsL0==0) continue;
      if (!(fInputsL0 & 1<<3)) continue;
      h0VIRall->Fill(srun,1);
      h0VIRbad->Fill(srun,0);
      if (!fIR1->TestBitNumber(90)) {
        printf("0VIR:");
        h0VIRbad->Fill(srun,1);
        for (Int_t i=-90;i<=90;i++) printf("%i",fIR1->TestBitNumber(90+i));
        printf("\n");
      }
    }
//    continue;
    Bool_t mb      = fClassesFired->String().Contains("CINT7-B");
    Bool_t sh1     = fClassesFired->String().Contains("CVHMSH1-B");
    Bool_t sh2     = fClassesFired->String().Contains("CVHMSH2-B");
    Bool_t vhm     = fClassesFired->String().Contains("CVHMV0M-B");
    Bool_t vhmnopf = fClassesFired->String().Contains("CVHMV0M-B-NOPF");
    Bool_t vhmspd1 = fClassesFired->String().Contains("CVHMV0M-B-SPD1");
    Bool_t sh2nopf = fClassesFired->String().Contains("CVHMSH2-B-NOPF");
    Bool_t sh2spd1 = fClassesFired->String().Contains("CVHMSH2-B-SPD1");
    Bool_t sh1nopf = fClassesFired->String().Contains("CVHMSH1-B-NOPF");
    
    
    
    tr[kINT]     = mb;
    tr[kSH1]     = sh1;
    tr[kSH2]     = sh2;
    tr[kVHM]     = vhm;
    
    if (mb)      hEventCounter->Fill("all","CINT7-B-NOPF"  ,srun,1.);
    if (vhmnopf) hEventCounter->Fill("all","CVHMV0M-B-NOPF",srun,1.);
    if (vhmspd1) hEventCounter->Fill("all","CVHMV0M-B-SPD1",srun,1.);
    if (sh2nopf) hEventCounter->Fill("all","CVHMSH2-B-NOPF",srun,1.);
    if (sh2spd1) hEventCounter->Fill("all","CVHMSH2-B-SPD1",srun,1.);
    if (sh1nopf) hEventCounter->Fill("all","CVHMSH1-B-NOPF",srun,1.);
    
    if (fIsIncomplete || fInputsL0==0) continue;

    Bool_t VBA = fInputsL0 & 1<< 0;
    Bool_t VBC = fInputsL0 & 1<< 1;
    Bool_t V0M = fInputsL0 & 1<< 9;
    Bool_t VHM = fInputsL0 & 1<<10;
    Bool_t SH2 = fInputsL0 & 1<<13;
    if (mb  && (!VBA || !VBC)) { printf("run=%i fInputsL0=%i mb  mismatch: VBA=%i VBC=%i\n",fRunNumber,fInputsL0,VBA,VBC); continue; }
    if (sh2 && (!VHM || !SH2)) { printf("run=%i fInputsL0=%i sh2 mismatch: VHM=%i SH2=%i\n",fRunNumber,fInputsL0,VHM,SH2); continue; }
    if (vhm && (!VHM || !V0M)) { printf("run=%i fInputsL0=%i vhm mismatch: VHM=%i V0M=%i\n",fRunNumber,fInputsL0,VHM,V0M); continue; }
    
    if (mb)      hEventCounter->Fill("complete","CINT7-B-NOPF"  ,srun,1.);
    if (vhmnopf) hEventCounter->Fill("complete","CVHMV0M-B-NOPF",srun,1.);
    if (vhmspd1) hEventCounter->Fill("complete","CVHMV0M-B-SPD1",srun,1.);
    if (sh2nopf) hEventCounter->Fill("complete","CVHMSH2-B-NOPF",srun,1.);
    if (sh2spd1) hEventCounter->Fill("complete","CVHMSH2-B-SPD1",srun,1.);
    if (sh1nopf) hEventCounter->Fill("complete","CVHMSH1-B-NOPF",srun,1.);

    Double_t meanOfV0M = hMeanOfV0M->GetBinContent(hMeanOfV0M->GetXaxis()->FindBin(srun));
    Double_t meanOnV0M = hMeanOnV0M->GetBinContent(hMeanOnV0M->GetXaxis()->FindBin(srun));
    Double_t meanOfFOR = hMeanOfFOR->GetBinContent(hMeanOfFOR->GetXaxis()->FindBin(srun));
    Double_t meanOnFOR = hMeanOnFOR->GetBinContent(hMeanOnFOR->GetXaxis()->FindBin(srun));
    Double_t meanOfTKL = hMeanOfTKL->GetBinContent(hMeanOfTKL->GetXaxis()->FindBin(srun));
    Double_t meanOfCL1 = hMeanOfCL1->GetBinContent(hMeanOfCL1->GetXaxis()->FindBin(srun));
    
    Int_t nBBA=0;      for (Int_t b=0;b<32;b++) nBBA+=fBBFlag[b+32];
    Int_t nBBC=0;      for (Int_t b=0;b<32;b++) nBBC+=fBBFlag[b];
    Int_t nBGA=0;      for (Int_t b=0;b<32;b++) nBGA+=fBGFlag[b+32];
    Int_t nBGC=0;      for (Int_t b=0;b<32;b++) nBGC+=fBGFlag[b];
    
    Bool_t bgdV0     = (fV0ADecision!=1 || fV0CDecision!=1);
    Bool_t goodV0    = fV0ADecision==1 && fV0CDecision==1;
    Int_t nClusters  = fNofITSClusters[0]+fNofITSClusters[1];
    Int_t nTracklets = fNofTracklets;
    Double_t V0C012  = fMRingV0C[0]+fMRingV0C[1]+fMRingV0C[2];
    Double_t V0C3    = fMRingV0C[3];
    Double_t onFOR   = fFOmap->CountBits(400);
    Double_t ofFOR   = fFiredChipMap->CountBits(400);
    Double_t ofCL1   = fNofITSClusters[1];
    Double_t ofTKL   = fNofTracklets;
    Double_t onV0M   = fTriggerChargeA+fTriggerChargeC;
    Double_t ofV0M   = fMTotV0A+fMTotV0C;

    Bool_t bgdSPD  = (nClusters<85 ? 0 : nClusters>85.+(2000.-85.)/(150.-7.)*(nTracklets-7.));
    Bool_t goodSPD = nClusters < 65+4*nTracklets;
    Bool_t isV0Casym = (V0C012 < 160.) || (V0C3 > 12.*TMath::Power(0.01*(V0C012-160.), 1.7)); 
    Bool_t isOutOfBunchPileup = IsOutOfBunchPileup(fIR1, fIR2, fRunNumber, fBC);
    Bool_t isEventSelected = goodV0 && !isOutOfBunchPileup && !fIsPileupSPD && goodSPD && isV0Casym;
    Bool_t isOutlierFOR = onFOR<200./300.*(ofFOR-50);
    Bool_t isOutlierV0M = onV0M<5000./800.*(ofV0M-200);

    if (goodV0){
      if (mb)      hEventCounter->Fill("good_V0","CINT7-B-NOPF",srun,1.);
      if (vhmnopf) hEventCounter->Fill("good_V0","CVHMV0M-B-NOPF",srun,1.);
      if (vhmspd1) hEventCounter->Fill("good_V0","CVHMV0M-B-SPD1",srun,1.);
      if (sh2nopf) hEventCounter->Fill("good_V0","CVHMSH2-B-NOPF",srun,1.);
      if (sh2spd1) hEventCounter->Fill("good_V0","CVHMSH2-B-SPD1",srun,1.);
      if (sh1nopf) hEventCounter->Fill("good_V0","CVHMSH1-B-NOPF",srun,1.);
      if (goodSPD) {
        if (mb)      hEventCounter->Fill("good_SPD","CINT7-B-NOPF",srun,1.);
        if (vhmnopf) hEventCounter->Fill("good_SPD","CVHMV0M-B-NOPF",srun,1.);
        if (vhmspd1) hEventCounter->Fill("good_SPD","CVHMV0M-B-SPD1",srun,1.);
        if (sh2nopf) hEventCounter->Fill("good_SPD","CVHMSH2-B-NOPF",srun,1.);
        if (sh2spd1) hEventCounter->Fill("good_SPD","CVHMSH2-B-SPD1",srun,1.);
        if (sh1nopf) hEventCounter->Fill("good_SPD","CVHMSH1-B-NOPF",srun,1.);
        if (!fIsPileupSPD) {
          if (mb)      hEventCounter->Fill("no_same_bunch_pileup","CINT7-B-NOPF",srun,1.);
          if (vhmnopf) hEventCounter->Fill("no_same_bunch_pileup","CVHMV0M-B-NOPF",srun,1.);
          if (vhmspd1) hEventCounter->Fill("no_same_bunch_pileup","CVHMV0M-B-SPD1",srun,1.);
          if (sh2nopf) hEventCounter->Fill("no_same_bunch_pileup","CVHMSH2-B-NOPF",srun,1.);
          if (sh2spd1) hEventCounter->Fill("no_same_bunch_pileup","CVHMSH2-B-SPD1",srun,1.);
          if (sh1nopf) hEventCounter->Fill("no_same_bunch_pileup","CVHMSH1-B-NOPF",srun,1.);
          if (!isOutOfBunchPileup) {
            if (mb)      hEventCounter->Fill("no_out_of_bunch_pileup","CINT7-B-NOPF",srun,1.);
            if (vhmnopf) hEventCounter->Fill("no_out_of_bunch_pileup","CVHMV0M-B-NOPF",srun,1.);
            if (vhmspd1) hEventCounter->Fill("no_out_of_bunch_pileup","CVHMV0M-B-SPD1",srun,1.);
            if (sh2nopf) hEventCounter->Fill("no_out_of_bunch_pileup","CVHMSH2-B-NOPF",srun,1.);
            if (sh2spd1) hEventCounter->Fill("no_out_of_bunch_pileup","CVHMSH2-B-SPD1",srun,1.);
            if (sh1nopf) hEventCounter->Fill("no_out_of_bunch_pileup","CVHMSH1-B-NOPF",srun,1.);
            if (isV0Casym) {
              if (mb)      hEventCounter->Fill("selected","CINT7-B-NOPF",srun,1.);
              if (vhmnopf) hEventCounter->Fill("selected","CVHMV0M-B-NOPF",srun,1.);
              if (vhmspd1) hEventCounter->Fill("selected","CVHMV0M-B-SPD1",srun,1.);
              if (sh2nopf) hEventCounter->Fill("selected","CVHMSH2-B-NOPF",srun,1.);
              if (sh2spd1) hEventCounter->Fill("selected","CVHMSH2-B-SPD1",srun,1.);
              if (sh1nopf) hEventCounter->Fill("selected","CVHMSH1-B-NOPF",srun,1.);
            }
          }
        }
      }
    }

    if (goodV0 && !isOutOfBunchPileup && !fIsPileupSPD && goodSPD) hV0C012vsV0C3->Fill(V0C012,V0C3);
    
    hPileupBcRun->Fill(0.,0.,srun,0.);
    if (mb & isOutlierFOR) for (Int_t i=-19;i<=20;i++) if (fIR2->TestBitNumber(90+i)) hPileupBcRun->Fill(i-1,fBC%4,srun,1);
////    if (sh2) for (Int_t i=-19;i<=20;i++) if (fIR1->TestBitNumber(90+i)) hPileupBcRun->Fill(i-1,fBC%4,srun,1);


    Double_t onV0Mnorm = onV0M/meanOnV0M;
    Double_t ofV0Mnorm = ofV0M/meanOfV0M;
    Double_t onFORnorm = onFOR/meanOnFOR;
    Double_t ofFORnorm = ofFOR/meanOfFOR;
    Double_t ofTKLnorm = ofTKL/meanOfTKL;
    Double_t ofCL1norm = ofCL1/meanOfCL1;
    
    if (ofFORnorm>=15) ofFOR=14.9; 
    if (onFORnorm>=15) onFOR=14.9;
    if (ofTKLnorm>=15) ofTKL=14.9;
    if (ofCL1norm>=15) ofCL1=14.9;
    if (onV0Mnorm>=15) onV0M=14.9;
    if (ofV0Mnorm>=15) ofV0M=14.9;
    if (ofFOR>=300)     ofFOR=299.9; 
    if (onFOR>=200)     onFOR=199.9;
    if (ofTKL>=200)     ofTKL=199.9;
    if (ofCL1>=400)     ofCL1=399.9;
    if (onV0M>=11000)   onV0M=10999.9;
    if (ofV0M>=2000)    ofV0M=1999.9;
    if (nTracklets>=400) nTracklets=399;
    if (nClusters>=4000) nClusters=3999;
    
    for (Int_t t=0;t<NTRIGGERS;t++){
//    for (Int_t t=0;t<1;t++){
      if (!tr[t]) continue;

      hClsVsTklall->Fill(t,nTracklets,nClusters,1);
      hBBAall->Fill(t,nBBA,srun,1);
      hBBCall->Fill(t,nBBC,srun,1);
      hBGAall->Fill(t,nBGA,srun,1);
      hBGCall->Fill(t,nBGC,srun,1);
      
      if (fV0ATime>=-45 && fV0ATime<=45) hV0ATime->Fill(t,fV0ATime,srun,1);
      if (fV0CTime>=-45 && fV0CTime<=45) hV0CTime->Fill(t,fV0CTime,srun,1);
      
      if (fV0ADecision==1) hV0ATimeBB->Fill(t,fV0ATime,srun,1);
      if (fV0CDecision==1) hV0CTimeBB->Fill(t,fV0CTime,srun,1);
      if (fV0ADecision==2) hV0ATimeBG->Fill(t,fV0ATime,srun,1);
      if (fV0CDecision==2) hV0CTimeBG->Fill(t,fV0CTime,srun,1);

      if (goodV0 && !isOutOfBunchPileup) hClsVsTklsel->Fill(t,nTracklets,nClusters,1);
      
      if (bgdSPD) {
        hClsVsTklbgd->Fill(t,nTracklets,nClusters,1);
        hBBAbgd->Fill(t,nBBA,srun,1);
        hBBCbgd->Fill(t,nBBC,srun,1);
        hBGAbgd->Fill(t,nBGA,srun,1);
        hBGCbgd->Fill(t,nBGC,srun,1);
      }
      if (t==0) {
        hOnFORvsOfFOR_all->Fill(ofFOR,onFOR,srun,1);
        hOnV0MvsOfV0M_all->Fill(ofV0M,onV0M,srun,1);
      }
      
      if (!isEventSelected) continue;
      
      hClsVsTklsig->Fill(t,nTracklets,nClusters,1);
      
      hBBAsig->Fill(t,nBBA,srun,1);
      hBBCsig->Fill(t,nBBC,srun,1);
      hBGAsig->Fill(t,nBGA,srun,1);
      hBGCsig->Fill(t,nBGC,srun,1);

      hOnV0M->Fill(t,onV0M,srun,1);
      hOfV0M->Fill(t,ofV0M,srun,1);
      hOfTKL->Fill(t,ofTKL,srun,1);
      hOfCL1->Fill(t,ofCL1,srun,1);
      hOnFOR->Fill(t,onFOR,srun,1);
      hOfFOR->Fill(t,ofFOR,srun,1);
      
      hOnV0Mnorm->Fill(t,onV0Mnorm,srun,1);
      hOfV0Mnorm->Fill(t,ofV0Mnorm,srun,1);
      hOfTKLnorm->Fill(t,ofTKLnorm,srun,1);
      hOfCL1norm->Fill(t,ofCL1norm,srun,1);
      hOnFORnorm->Fill(t,onFORnorm,srun,1);
      hOfFORnorm->Fill(t,ofFORnorm,srun,1);
      if (t==0) {
      if (fBC%4==0) hOnFORvsOfFOR_sel_mod0->Fill(ofFOR,onFOR,srun,1);
      if (fBC%4==1) hOnFORvsOfFOR_sel_mod1->Fill(ofFOR,onFOR,srun,1);
      if (fBC%4==2) hOnFORvsOfFOR_sel_mod2->Fill(ofFOR,onFOR,srun,1);
      if (fBC%4==3) hOnFORvsOfFOR_sel_mod3->Fill(ofFOR,onFOR,srun,1);
      
      if (fBC%4==0) hOnV0MvsOfV0M_sel_mod0->Fill(ofV0M,onV0M,srun,1);
      if (fBC%4==1) hOnV0MvsOfV0M_sel_mod1->Fill(ofV0M,onV0M,srun,1);
      if (fBC%4==2) hOnV0MvsOfV0M_sel_mod2->Fill(ofV0M,onV0M,srun,1);
      if (fBC%4==3) hOnV0MvsOfV0M_sel_mod3->Fill(ofV0M,onV0M,srun,1);
      }
    }
    
    

  }
  
  TFile* fout = new TFile("results.root","recreate");
//  hEventCounter->Write();
//  hClsVsTklAll->Write();
//  hClsVsTkl->Write();
//  hOnVsOffSPD1->Write();
//  hOnVsOffSPD2->Write();
//  hOnVsOffSPD1vsRun->Write();
//  hOnVsOffSPD2vsRun->Write();
//  hOnVsOffV0M1->Write();
//  hOnVsOffV0M2->Write();
//  hPileup->Write();
//  hBBA_SMB_all->Write();
//  hBBC_SMB_all->Write();
//  hBGA_SMB_all->Write();
//  hBGC_SMB_all->Write();
//  hBBA_SMB_sig->Write();
//  hBBC_SMB_sig->Write();
//  hBGA_SMB_sig->Write();
//  hBGC_SMB_sig->Write();
//  hBBA_SMB_bgd->Write();
//  hBBC_SMB_bgd->Write();
//  hBGA_SMB_bgd->Write();
//  hBGC_SMB_bgd->Write();
//  hBBA_VMB_all->Write();
//  hBBC_VMB_all->Write();
//  hBGA_VMB_all->Write();
//  hBGC_VMB_all->Write();
//  hBBA_VMB_sig->Write();
//  hBBC_VMB_sig->Write();
//  hBGA_VMB_sig->Write();
//  hBGC_VMB_sig->Write();
//  hBBA_VMB_bgd->Write();
//  hBBC_VMB_bgd->Write();
//  hBGA_VMB_bgd->Write();
//  hBGC_VMB_bgd->Write();
//  hBBA_SHM_all->Write();
//  hBBC_SHM_all->Write();
//  hBGA_SHM_all->Write();
//  hBGC_SHM_all->Write();
//  hBBA_VHM_all->Write();
//  hBBC_VHM_all->Write();
//  hBGA_VHM_all->Write();
//  hBGC_VHM_all->Write();
//  hBBA_SHM_bgd->Write();
//  hBBC_SHM_bgd->Write();
//  hBGA_SHM_bgd->Write();
//  hBGC_SHM_bgd->Write();
//  hBBA_VHM_bgd->Write();
//  hBBC_VHM_bgd->Write();
//  hBGA_VHM_bgd->Write();
//  hBGC_VHM_bgd->Write();
//  hBBA_SHM_sig->Write();
//  hBBC_SHM_sig->Write();
//  hBGA_SHM_sig->Write();
//  hBGC_SHM_sig->Write();
//  hBBA_VHM_sig->Write();
//  hBBC_VHM_sig->Write();
//  hBGA_VHM_sig->Write();
//  hBGC_VHM_sig->Write();
//  hOnlineV0M->Write();
//  hOnlineSPD->Write();
//  hOnlineVHM->Write();
//  hOnlineSHM->Write();
//  hOfflineV0M->Write();
//  hOfflineVHM->Write();
//  hOfflineV0Mcut->Write();
//  hOfflineVHMcut->Write();
//  hOfflineSPD->Write();
//  hOfflineSHM->Write();
//  hV0ATime->Write();
//  hV0CTime->Write();
//  hV0ATimeNoOutOfBunch->Write();
//  hV0CTimeNoOutOfBunch->Write();
//  hV0ATimeBB->Write();

  
  hEventCounter->Write();
  hPileupBcRun->Write();
  hClsVsTklall->Write();
  hClsVsTklbgd->Write();
  hClsVsTklsig->Write();
  hClsVsTklsel->Write();
  hBBAall->Write();
  hBBCall->Write();
  hBGAall->Write();
  hBGCall->Write();
  hBBAbgd->Write();
  hBBCbgd->Write();
  hBGAbgd->Write();
  hBGCbgd->Write();
  hBBAsig->Write();
  hBBCsig->Write();
  hBGAsig->Write();
  hBGCsig->Write();
  hV0ATime->Write();
  hV0CTime->Write();
  hV0ATimeBB->Write();
  hV0CTimeBB->Write();
  hV0ATimeBG->Write();
  hV0CTimeBG->Write();
  hOnV0M->Write();
  hOfV0M->Write();
  hOfTKL->Write();
  hOfCL1->Write();
  hOnFOR->Write();
  hOfFOR->Write();
  hOnV0Mnorm->Write();
  hOfV0Mnorm->Write();
  hOfTKLnorm->Write();
  hOfCL1norm->Write();
  hOnFORnorm->Write();
  hOfFORnorm->Write();
  hOnFORvsOfFOR_all->Write();
  hOnFORvsOfFOR_sel_mod0->Write();
  hOnFORvsOfFOR_sel_mod1->Write();
  hOnFORvsOfFOR_sel_mod2->Write();
  hOnFORvsOfFOR_sel_mod3->Write();
  hOnV0MvsOfV0M_all->Write();
  hOnV0MvsOfV0M_sel_mod0 ->Write();
  hOnV0MvsOfV0M_sel_mod1->Write();
  hOnV0MvsOfV0M_sel_mod2->Write();
  hOnV0MvsOfV0M_sel_mod3->Write();
  hV0C012vsV0C3->Write();
  h0TVXbcs->Write();
  h0TVXall->Write();
  h0TVXbad->Write();
  h0VIRall->Write();
  h0VIRbad->Write();
  fout->Close();
}


Bool_t IsOutOfBunchPileup(TBits* fIR1, TBits* fIR2, Int_t fRunNumber, UShort_t fBC) {
  Bool_t isOutOfBunchPileup = 0;
  if ((fRunNumber>=225000 && fRunNumber<=225719) ||
      (fRunNumber>=226062 && fRunNumber<=226500) ||
      (fRunNumber>=233912 && fRunNumber<=234050) 
      ) isOutOfBunchPileup=0; // isolated bunches
  else {
    Int_t pf1 = -4-fBC%4; // standard D-1;D0;D+1
    Int_t pf2 = +7-fBC%4; // standard D-1;D0;D+1
    if (fRunNumber<=238621) { pf1+=4; pf2+=4; } // D0;D+1;D+2: shifted settings in LHC15f-j
    if (fRunNumber==240404) { pf1+=4; pf2+=4; } // test in LHC15k with D0;D+1;D+2
    if (fRunNumber==240411) { pf1+=0; pf2-=4; } // test in LHC15k with D-1;D0
    if (fRunNumber==240427) { pf1-=4; pf2-=4; } // test in LHC15k with D-2;D-1;D0
    
    // temporary staff to be adjusted with friend info
    Int_t pf2maxForT0 = pf2;
    Int_t ir1skip = 0;
    if (fRunNumber<=228948) { pf2maxForT0=8; }  // temporary: 0UBA|0UBC IR2 with afterpulses starting from +9
    if (fRunNumber<=238187) { ir1skip=2;     }  // temporary: for 0VIR(>=2) or VBA&VBC or VBA|VBC
    
    for (Int_t i=pf1;i<=pf2;i++) {
      if (i==0) continue;
      if (i<=pf2maxForT0) isOutOfBunchPileup|=fIR2->TestBitNumber(90+i); // T0-based clean-up
      if (i>0 && i<=ir1skip) continue; // skip next 2 for old IR definitions
      isOutOfBunchPileup|=fIR1->TestBitNumber(90+i); // V0-based clean-up
    }
  }
  return isOutOfBunchPileup;
}
