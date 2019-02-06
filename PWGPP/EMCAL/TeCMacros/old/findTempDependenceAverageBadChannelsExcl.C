// Fill and save histograms with info from the DA object
#if !defined( __CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TTree.h>
#include "/opt/alice/alice/v4-18-Test/EMCAL/AliCaloCalibSignal.h"
#include "/opt/alice/alice/v4-18-Test/STEER/AliCDBEntry.h"

#endif

// global EMCal numbers
const int kNSM = 20; 
const int kNSectors = 2;
const int kNSides = 2;
const int kNPer = 15;

// these other numbers are per SuperModule
const int kNCol = 48;
const int kNRow = 24;
const int kNStrip = kNCol / 2;

// amplitude cut, in ADCs
const int kAmpCut = -1; // -1 => include everything for now..
 
// which gains to plot?
const int towerGain = 1;
const int stripGain = 0; 
const int fileNum=629;  // Number of files in the "hist2" directory
const int fileNumTemp=908;

void findTempDependenceAverageBadChannelsExcl(const int runno=219496,
	     const int debug = 0)

{

 char idref[100];
 char idln[100];
 char name1[30];
 char name2[30];
 char name3[30];
 char name4[30];
 Int_t ntowers;

 TH1F *hRatiosForTowers[kNSM][fileNum];
 TH1F *herrorsForTowers[kNSM][fileNum];
 TH2F *hTowerTempDep[kNSM][1152];
 TProfile *hTowerTempDepProf[kNSM][1152];
 TH2F *hTowerErrorTempDep[kNSM][1152];

 // 15 Feb 2018: add bad channels removal

 TH2I *fBad[kNSM][kNPer];
 TH1I *fBadLn[kNSM][kNPer];
 for (int iSM=0; iSM<kNSM; iSM++) {
   for (int iPer=0; iPer<kNPer; iPer++) {
     sprintf(idref, "fBad%d_%d", iSM, iPer);
     sprintf(idln, "fBadLn%d_%d", iSM, iPer);
     fBad[iSM][iPer]=new TH2I (idref, " ", kNCol, 0, kNCol, kNRow, 0., kNRow);
     fBadLn[iSM][iPer]=new TH1I (idln, " ", 1152, -0.5, 1151.5);
   }
 }
 
 Int_t fRunNumber[15]={225305, 225705, 235709, 236892, 238896, 239396, 240610, 241001, 241521, 244340, 244531, 244911, 245496, 246089, 246390};
 TFile *fil1 = TFile::Open("EMCALBadChannels.root");
 AliOADBContainer *contBC = new AliOADBContainer("");
 contBC->InitFromFile("EMCALBadChannels.root","AliEMCALBadChannels");

 for (Int_t l=0; l<15; l++) {
   
   //Printf("");
   //Printf("Period = %d", l);
   
    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRunNumber[l]);
    if (!arrayBC)
      {
        Printf("No external hot channel set for run number: %d", fRunNumber[l]);
      }
    
    for (Int_t i=0; i<10; ++i)
      {
        fBad[i][l] = (TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
	
        if (!fBad[i][l])
          {
            Printf("Can not get EMCALBadChannelMap_Mod%d",i);
            continue;
          }
      }
 }
 
 for (Int_t lp=0; lp<15; lp++) {
   for (Int_t iSM=0; iSM<kNSM; iSM++) {
     ntowers=-1;
     for (Int_t ibinY=1; ibinY<=kNRow; ibinY++) {
       for (Int_t ibinX=1; ibinX<=kNCol; ibinX++) {
         ntowers=ntowers+1;
	 fBadLn[iSM][lp]->SetBinContent(ntowers+1, fBad[iSM][lp]->GetBinContent(ibinX, ibinY));
       }       
     }
   }
 }
 
 // end 15 Feb 2018

 for (Int_t iSM=0; iSM<kNSM; iSM++) {
   for (Int_t i=0; i<1152; i++) {
     sprintf(name1, "hTowerErrorTempDep_%d_%d", iSM, i);
     hTowerErrorTempDep[iSM][i] = new TH2F(name1, " ", 120, 21.0, 27.0, 200, 0., 1.0);
     sprintf(name3, "hTowerTempDep_%d_%d", iSM, i);
     hTowerTempDep[iSM][i] = new TH2F(name3, " ", 120, 21.0, 27.0, 1500, 0.5, 2.0);
     sprintf(name2, "hTowerTempDepProf_%d_%d", iSM, i);

     if ((iSM==6)||(iSM==8)||(iSM==9)||(iSM==7)) hTowerTempDepProf[iSM][i] = new TProfile(name2, " ", 120, 21.0, 27.0, 0.5, 2.0, "s");
     if ((iSM==1)||(iSM==2)||(iSM==4)||(iSM==3)||(iSM==0)||(iSM==5)) hTowerTempDepProf[iSM][i] = new TProfile(name2, " ", 120, 21.0, 27.0, 0.5, 2.0, "s");
     // if (iSM==3) hTowerTempDepProf[iSM][i] = new TProfile(name2, " ", 120, 21.0, 27.0, 0.7, 1.4, "s");
     // if ((iSM==6)||(iSM==7)) hTowerTempDepProf[iSM][i] = new TProfile(name2, " ", 60, 22.0, 25.0, 0.7, 1.4, "s");
     // if (iSM==5) hTowerTempDepProf[iSM][i] = new TProfile(name2, " ", 60, 25.0, 28.0, 0.7, 1.4, "s");
   }
 }

 TFile *fil1 = TFile::Open("LED/factorsFolder/factorsNew.root");
 TH1F *dataRunNumberHisto=(TH1F*)fil1->Get("dataRunNumberHisto");
   for (int iSM=0; iSM<kNSM; iSM++) {
     for (Int_t ig=0; ig<fileNum; ig++) {
       sprintf(name2, "hRatios_%d_%d", iSM, ig);
       hRatiosForTowers[iSM][ig]=(TH1F*)fil1->Get(name2);
       sprintf(name3, "herrorOfRatios_%d_%d", iSM, ig);
       herrorsForTowers[iSM][ig]=(TH1F*)fil1->Get(name3);
     }
   }

 TFile *fil2 = TFile::Open("TemperatureFiles/tempDependenceAllRuns.root");
 TH1F *runNumberHisto=(TH1F*)fil2->Get("hRunNumber");
 TH1F *tempHistoSM0=(TH1F*)fil2->Get("hTemperatureSM_0");
 TH1F *tempHistoSM1=(TH1F*)fil2->Get("hTemperatureSM_1");
 TH1F *tempHistoSM2=(TH1F*)fil2->Get("hTemperatureSM_2");
 TH1F *tempHistoSM3=(TH1F*)fil2->Get("hTemperatureSM_3");
 TH1F *tempHistoSM4=(TH1F*)fil2->Get("hTemperatureSM_4");
 TH1F *tempHistoSM5=(TH1F*)fil2->Get("hTemperatureSM_5");
 TH1F *tempHistoSM6=(TH1F*)fil2->Get("hTemperatureSM_6");
 TH1F *tempHistoSM7=(TH1F*)fil2->Get("hTemperatureSM_7");
 TH1F *tempHistoSM8=(TH1F*)fil2->Get("hTemperatureSM_8");
 TH1F *tempHistoSM9=(TH1F*)fil2->Get("hTemperatureSM_9");

 Int_t tempbin;
 for (int iSM=0; iSM<kNSM; iSM++) {
   for (Int_t k=0; k<1152; k++) {
     for (Int_t m=0; m<fileNum; m++) {
       //if (m<91) continue;
       for (Int_t j=0; j<fileNumTemp; j++) {
	 if(runNumberHisto->GetBinContent(j+1)==dataRunNumberHisto->GetBinContent(m+1)) {

	   Int_t iperiod=0;
	   for (Int_t ip=0; ip<kNPer; ip++) {
	     if ((runNumberHisto->GetBinContent(j+1)>=fRunNumber[ip])&&(runNumberHisto->GetBinContent(j+1)<fRunNumber[ip+1]))
	       { iperiod=ip;  break;}
	   }

	   
	   if(fBadLn[iSM][iperiod]->GetBinContent(k+1)==0) { // make it in 1D!!!


	   
	   if (hRatiosForTowers[iSM][m]->GetBinContent(k+1)!=0.) {
	     if (iSM==0) hTowerTempDep[iSM][k]->Fill(tempHistoSM0->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==0) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM0->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==0) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM0->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==1) hTowerTempDep[iSM][k]->Fill(tempHistoSM1->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==1) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM1->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==1) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM1->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==2) hTowerTempDep[iSM][k]->Fill(tempHistoSM2->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==2) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM2->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==2) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM2->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==3) hTowerTempDep[iSM][k]->Fill(tempHistoSM3->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==3) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM3->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==3) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM3->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==4) hTowerTempDep[iSM][k]->Fill(tempHistoSM4->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==4) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM4->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==4) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM4->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==5) hTowerTempDep[iSM][k]->Fill(tempHistoSM5->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==5) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM5->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==5) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM5->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==6) hTowerTempDep[iSM][k]->Fill(tempHistoSM6->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==6) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM6->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==6) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM6->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==7) hTowerTempDep[iSM][k]->Fill(tempHistoSM7->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==7) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM7->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==7) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM7->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==8) hTowerTempDep[iSM][k]->Fill(tempHistoSM8->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==8) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM8->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==8) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM8->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==9) hTowerTempDep[iSM][k]->Fill(tempHistoSM9->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==9) hTowerErrorTempDep[iSM][k]->Fill(tempHistoSM9->GetBinContent(j+1), herrorsForTowers[iSM][m]->GetBinContent(k+1));
	     if (iSM==9) hTowerTempDepProf[iSM][k]->Fill(tempHistoSM9->GetBinContent(j+1), hRatiosForTowers[iSM][m]->GetBinContent(k+1));

	   }
	 }
	   
	   break;
	 }
       }
     }

   }
 }

 TFile* outputFile = TFile::Open("temperatureDependenceOfTheRatioBadChannelsExcl.root", "recreate");
 if (!outputFile || !outputFile->IsOpen()) {
   Printf("Cannot Open OUTPUT file");
   return;
 }
 for (int iSM=0; iSM<kNSM; iSM++) {
   for (Int_t k=0; k<1152; k++) {
     hTowerTempDep[iSM][k]->Write();
     hTowerTempDepProf[iSM][k]->Write();
     hTowerErrorTempDep[iSM][k]->Write();
   }
 }

}

