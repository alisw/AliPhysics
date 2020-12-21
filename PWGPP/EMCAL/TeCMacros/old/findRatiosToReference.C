// Make the ratio plots between current run number and the reference run
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

// these other numbers are per SuperModule
const int kNCol = 48;
const int kNRow = 24;
const int kNStrip = kNCol / 2;

// amplitude cut, in ADCs
const int kAmpCut = -1; // -1 => include everything for now..
 
// which gains to plot?
const int towerGain = 1;
const int stripGain = 0;
//old ref 243522
void findRatiosToReference(const char *listfile = "histosFilelist.txt", const int fileNum=591, const int refRun=230701, 
	     const int debug = 0)
 
{

 char * sideStr[] = {"A", "C"};
 char filename[100];
 TH2F *h1[kNSM];
 TH2F *hAmp[kNSM];
 TH2F *herrorCurrentAmp[kNSM];
 TH2F *herrorRefLED[kNSM];
 TH1F *hLEDMonError[kNSM];
 TH1F *hLEDMonRefRunError[kNSM];
 TH2F *herrorRefAmp[kNSM];
 TH2F *hRefAmp[kNSM];
 TH1F *hRefLED[kNSM];
 TH2F *hCurAmp[kNSM];
 TH1F *hCurLED[kNSM];
 char id[100];
 char id2[100];
 char idref[100];
 char title[100]; 
 char name1[30];
 char name2[30];
 char name3[30];
 char refrunfile[30];

 TH1F *hRatios[kNSM][fileNum];
 TH1F *herrorOfRatios[kNSM][fileNum];

 Int_t nfiles=-1;

 for (Int_t iSM=0; iSM<kNSM; iSM++) {
   for (Int_t k=0; k<fileNum; k++) {
     sprintf(name2, "hRatios_%d_%d", iSM, k);
     hRatios[iSM][k] = new TH1F(name2, " ", 1152, -0.5, 1151.5);
     sprintf(name3, "herrorOfRatios_%d_%d", iSM, k);
     herrorOfRatios[iSM][k] = new TH1F(name3, " ", 1152, -0.5, 1151.5);
   }
 }

 TH1F *dataRunNumberHisto=new TH1F("dataRunNumberHisto", "  " , fileNum, -0.5, fileNum-0.5);
 
 sprintf(refrunfile, "hist5/LED_000%dNew.root", refRun);
 TFile *fil1 = TFile::Open(refrunfile);
 for (int iSM=0; iSM<kNSM; iSM++) {
   sprintf(idref, "hAmpOverMon%02d%d", iSM, towerGain);
   sprintf(name1, "hsmLEDAmpRmsOverMean%d", iSM);
   sprintf(name3, "hsmAmpRmsOverMean%d", iSM);
   sprintf(name2, "hCount%02d%d", iSM, towerGain);
   hAmp[iSM]=(TH2F*)fil1->Get(idref);
   hLEDMonRefRunError[iSM]=(TH1F*)fil1->Get(name1);
   herrorRefAmp[iSM]=(TH2F*)fil1->Get(name3);
   sprintf(id, "h%02d%d", iSM, towerGain);
   hRefAmp[iSM]=(TH2F*)fil1->Get(id);
   sprintf(id2, "hStrip%02d%d", iSM, towerGain);
   hRefLED[iSM]=(TH1F*)fil1->Get(id2);
 }

  ifstream fin(listfile);
  Int_t runno = 0;
  char fileRoot[100];

  while ( fin.good() ) {
    fin >> runno >> fileRoot;
    if ( fin.good() ) {
      cout << " fileRoot " << fileRoot << " runno " << runno << endl;

   TFile *f = TFile::Open(fileRoot);
   if (!f || !f->IsOpen()) {
     continue;
   }

   nfiles=nfiles+1;
   dataRunNumberHisto->Fill(nfiles, runno);
   Printf("1");
   
   for (int iSM=0; iSM<kNSM; iSM++) {
     sprintf(idref, "hAmpOverMon%02d%d", iSM, towerGain);
     sprintf(name1, "hsmLEDAmpRmsOverMean%d", iSM);
     sprintf(name3, "hsmAmpRmsOverMean%d", iSM);
     sprintf(name2, "hCount%02d%d", iSM, towerGain);
     h1[iSM]=(TH2F*)f->Get(idref);
     hLEDMonError[iSM]=(TH1F*)f->Get(name1);
     herrorCurrentAmp[iSM]=(TH2F*)f->Get(name2);
     sprintf(id, "h%02d%d", iSM, towerGain);
     hCurAmp[iSM]=(TH2F*)f->Get(id);
     sprintf(id2, "hStrip%02d%d", iSM, towerGain);
     hCurLED[iSM]=(TH1F*)f->Get(id2);
   }

   Double_t errorOfFactor=0.0;
   Double_t errorNum, errorDen;
   Double_t factorT;
   Int_t ntowers;
   for (int iSM=0; iSM<kNSM; iSM++) {
     ntowers=-1;
     for (int ibinY=1; ibinY<=kNRow; ibinY++) {
       for (int ibinX=1; ibinX<=kNCol; ibinX++) {
	 ntowers=ntowers+1;

 	 //if (h1[iSM]->GetBinContent(ibinX, ibinY)>0.0) 	 errorOfFactor=TMath::Sqrt(((herrorRef[iSM]->GetBinContent(ibinX, ibinY)/h1[iSM]->GetBinContent(ibinX, ibinY))*(herrorRef[iSM]->GetBinContent(ibinX, ibinY)/h1[iSM]->GetBinContent(ibinX, ibinY)))+((hAmp[iSM]->GetBinContent(ibinX, ibinY)*herrorCurrent[iSM]->GetBinContent(ibinX, ibinY))/(h1[iSM]->GetBinContent(ibinX, ibinY)*h1[iSM]->GetBinContent(ibinX, ibinY)))*((hAmp[iSM]->GetBinContent(ibinX, ibinY)*herrorCurrent[iSM]->GetBinContent(ibinX, ibinY))/(h1[iSM]->GetBinContent(ibinX, ibinY)*h1[iSM]->GetBinContent(ibinX, ibinY))));

	 errorNum=(hAmp[iSM]->GetBinContent(ibinX, ibinY))*TMath::Sqrt((herrorRefAmp[iSM]->GetBinContent(ibinX, ibinY)/hRefAmp[iSM]->GetBinContent(ibinX, ibinY))*(herrorRefAmp[iSM]->GetBinContent(ibinX, ibinY)/hRefAmp[iSM]->GetBinContent(ibinX, ibinY))+(hLEDMonRefRunError[iSM]->GetBinContent(ibinX, ibinY)/hRefLED[iSM]->GetBinContent(ibinX, ibinY))*(hLEDMonRefRunError[iSM]->GetBinContent(ibinX, ibinY)/hRefLED[iSM]->GetBinContent(ibinX, ibinY)));
	 
	 errorDen=(h1[iSM]->GetBinContent(ibinX, ibinY))*TMath::Sqrt((herrorCurrentAmp[iSM]->GetBinContent(ibinX, ibinY)/hCurAmp[iSM]->GetBinContent(ibinX, ibinY))*(herrorCurrentAmp[iSM]->GetBinContent(ibinX, ibinY)/hCurAmp[iSM]->GetBinContent(ibinX, ibinY))+ (hLEDMonError[iSM]->GetBinContent(ibinX, ibinY)/hCurLED[iSM]->GetBinContent(ibinX, ibinY))*(hLEDMonError[iSM]->GetBinContent(ibinX, ibinY)/hCurLED[iSM]->GetBinContent(ibinX, ibinY)));

	 if (h1[iSM]->GetBinContent(ibinX, ibinY)!=0.0) factorT=hAmp[iSM]->GetBinContent(ibinX, ibinY)/h1[iSM]->GetBinContent(ibinX, ibinY);
	 if ((hAmp[iSM]->GetBinContent(ibinX, ibinY)!=0.0)&&(h1[iSM]->GetBinContent(ibinX, ibinY)!=0.0))
	 errorOfFactor=factorT*TMath::Sqrt((errorNum/hAmp[iSM]->GetBinContent(ibinX, ibinY))*(errorNum/hAmp[iSM]->GetBinContent(ibinX, ibinY))+(errorDen/h1[iSM]->GetBinContent(ibinX, ibinY))*(errorDen/h1[iSM]->GetBinContent(ibinX, ibinY)));


	 if (iSM==0)  Printf("factor = %f, error = %f", factorT, errorOfFactor);    
 	 if (h1[iSM]->GetBinContent(ibinX, ibinY)!=0.0) {
 	     hRatios[iSM][nfiles]->Fill(ntowers, hAmp[iSM]->GetBinContent(ibinX, ibinY)/h1[iSM]->GetBinContent(ibinX, ibinY));
 	     if (hAmp[iSM]->GetBinContent(ibinX, ibinY)!=0.0) herrorOfRatios[iSM][nfiles]->Fill(ntowers, errorOfFactor);
 	   }
	 else
	   if ((hAmp[iSM]->GetBinContent(ibinX, ibinY)==0.0)||(h1[iSM]->GetBinContent(ibinX, ibinY)==0.0))
	     { 
	       hRatios[iSM][nfiles]->Fill(ntowers, 0.0);
	       herrorOfRatios[iSM][nfiles]->Fill(ntowers, 0.0); 
	     }
	 
       }
     }

   }
    
   f->Close();
    }
 }

 TFile* outputFile = TFile::Open("factorsFolder/factorsNewAllSMs.root", "recreate");
 if (!outputFile || !outputFile->IsOpen()) {
   Printf("Cannot Open OUTPUT file");
   return;
 }

 dataRunNumberHisto->Write();
   for (int iSM=0; iSM<kNSM; iSM++) {
     hAmp[iSM]->Write();
     for (Int_t ig=0; ig<fileNum; ig++) {
       hRatios[iSM][ig]->Write();
       herrorOfRatios[iSM][ig]->Write();
     }
   }
}

