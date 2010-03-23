#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGrid.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliPDG.h"
#include "STEER/AliESDEvent.h"
#include "STEER/AliESDZDC.h"

#endif
void CheckAlienZDCESD(Int_t year=2010, const char* period="10a", 
	Int_t nRun=0, Int_t recoPass=1, Int_t nMaxFiles=1,
	Bool_t plot=kTRUE, Bool_t esdWordCut=kFALSE)
{
  
 if(nRun==0){
   printf("\n\n YOU MUST PROVIDE A RUN NUMBER!!! \n\n");
   return;
 }
  
 // Histogram definition
 // ----------------------------------------------------------------
/* TH2F *centroidZNsideC = new TH2F("centroidZNsideC","Impact point over ZNC",100,-5.,5.,100,-5.,5.);
  centroidZNsideC->SetXTitle("X_{ZNC} (cm)");
  centroidZNsideC->SetYTitle("Y_{ZNC} (cm)");
 TH2F * centroidZNsideA = new TH2F("centroidZNsideA","Impact point over ZNA",100,-5.,5.,100,-5.,5.);
  centroidZNsideA->SetXTitle("X_{ZNA} (cm)");
  centroidZNsideA->SetYTitle("Y_{ZNA} (cm)");
*/
 TH1F * enZNC = new TH1F("enZNC", "ZNC signal",100,0.,2000.);
  enZNC->SetXTitle("E (GeV)");
 TH1F * enZPC = new TH1F("enZPC", "ZPC signal",100,0.,2000.);
  enZPC->SetXTitle("E (GeV)");
 TH1F * enZNA = new TH1F("enZNA", "ZNA signal",100,0.,2000.);
  enZNA->SetXTitle("E (GeV)");
 TH1F * enZPA = new TH1F("enZPA", "ZPA signal",100,0.,2000.);
  enZPA->SetXTitle("E (GeV)");
 TH1D * enZEM1 = new TH1D("enZEM1", "Energy in ZEM1",100,0.,2000.);
  enZEM1->SetXTitle("E (GeV)");
 TH1D * enZEM2 = new TH1D("enZEM2", "Energy in ZEM2",100,0.,2000.);
  enZEM2->SetXTitle("E (GeV)");
 // ----------------------------------------------------------------
 TH1D * hZNCTow[5]; TH1D * hZPCTow[5]; 
 TH1D * hZNATow[5]; TH1D * hZPATow[5]; 
 char nomehistznc[30], nomehistzpc[30], nomehistzna[30], nomehistzpa[30];
 for(Int_t i=0; i<5; i++){
   sprintf(nomehistznc,"ZNC-pm%d",i);
   hZNCTow[i] = new TH1D(nomehistznc, nomehistznc, 100, 0.,1000.);
   sprintf(nomehistzpc,"ZPC-pm%d",i);
   hZPCTow[i] = new TH1D(nomehistzpc, nomehistzpc, 100, 0.,1000.);
   sprintf(nomehistzna,"ZNA-pm%d",i);
   hZNATow[i] = new TH1D(nomehistzna, nomehistzna, 100, 0.,1000.);
   sprintf(nomehistzpa,"ZPA-pm%d",i);
   hZPATow[i] = new TH1D(nomehistzpa, nomehistzpa, 100, 0.,1000.);
 }
 //
/* TH1D *hSumQZNC = new TH1D("hSumQZNC", "hSumQZNC", 100, 0., 1000.);
 TH1D *hSumQZPC = new TH1D("hSumQZPC", "hSumQZPC", 100, 0., 1000.);
 TH1D *hSumQZNA = new TH1D("hSumQZNA", "hSumQZNA", 100, 0., 1000.);
 TH1D *hSumQZPA = new TH1D("hSumQZPA", "hSumQZPA", 100, 0., 1000.);
*/
 //
 TH1F *hESDword = new TH1F("hESDword","hESDword",6,0.5,6.5);
 hESDword->SetXTitle("ZDC trigger pattern"); 

 TGrid::Connect("alien:",0,0,"t");
 gSystem->Exec(Form("gbbox find \"/alice/data/%d/LHC%s/000%d/ESDs/pass%d\" \"AliESDs.root\" > ESDFiles.txt",
       year, period, nRun, recoPass));
 FILE* listruns=fopen("ESDFiles.txt","r");
 
 char esdFileName[200], filnamalien[200];
 char yperiod, dirESD;
 // 
 Int_t nAnalyzedFiles=0;
 Int_t nevPhys=0, nevZNC=0, nevZPC=0, nevZNA=0, nevZPA=0, nevZEM1=0, nevZEM2=0;
 
 while(!feof(listruns)){
  
 if(nAnalyzedFiles!=nMaxFiles){
  
  int st = fscanf(listruns,"%s\n",esdFileName);    
  Char_t directory[100];
  sprintf(directory,"/alice/data/%d",year);
  if(!strstr(esdFileName,directory)) continue;
  sscanf(esdFileName,"/alice/data/%d/LHC%s/000%d/ESDs/pass%d/%s/AliESDs.root",&year,&yperiod,&nRun,&recoPass,&dirESD);
  sprintf(filnamalien,"alien://%s",esdFileName);
  printf("\n Opening file: %s\n",filnamalien);
 
  // open the ESD file
  TFile* esdFile = TFile::Open(filnamalien);
  if(!esdFile) {
    Error("CheckZDCESD", "opening ESD file %s failed",filnamalien);
    return;
  }
  
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if(!tree) {
    Error("CheckZDCESD", "No ESD tree found");
    return;
  }
  
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("AliESDHeader*", 1);
  tree->SetBranchStatus("AliESDRun*", 1);
  tree->SetBranchStatus("AliESDZDC*", 1);
  tree->SetBranchStatus("PrimaryVertex*", 1);
  tree->SetBranchStatus("SPDVertex*", 1);
  tree->SetBranchStatus("AliESDVZERO*", 1);

  esd->ReadFromTree(tree);
  
  Int_t nevents = (Int_t)tree->GetEntries();
  printf("\n    No. of events in ESD tree = %d\n", nevents);
  for(Int_t iEvent=0; iEvent<nevents; iEvent++){
    // get the event summary data
    tree->GetEvent(iEvent);
    
    //printf("    ev. type %d\n",esd->GetEventType());
    //
    if(esd->GetEventType() == 7){
      nevPhys++;
      if(!esd) {
        Error("CheckESD", "no ESD object found for event %d", iEvent);
        return;
      }
    
      Double_t sumQznc=0., sumQzpc=0., sumQzna=0., sumQzpa=0.;
    
      AliESDZDC *esdZDC = esd->GetESDZDC();
      //Double_t centrZNC={-999.,-999.}, centrZNA={-999.,-999.};
      //esdZDC->GetZNCentroidInpp(centrZNC, centrZNA);
      //Short_t npart = esdZDC->GetZDCParticipants();
      Double_t energyZNC = esdZDC->GetZDCN1Energy();
      Double_t energyZPC = esdZDC->GetZDCP1Energy();
      Double_t energyZNA = esdZDC->GetZDCN2Energy();
      Double_t energyZPA = esdZDC->GetZDCP2Energy();
      Double_t energyZEM1 = esdZDC->GetZDCEMEnergy(0);
      Double_t energyZEM2 = esdZDC->GetZDCEMEnergy(1);
      const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
      const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
      const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
      const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
      UInt_t iWord = esdZDC->GetESDQuality();
      //
      if((iWord & 0x00000001) == 0x00000001){
        nevZNA++;   
	hESDword->Fill(1.);
      }
      if((iWord & 0x00000002) == 0x00000002){
        nevZPA++;   
	hESDword->Fill(2.);
      }
      if((iWord & 0x00000004) == 0x00000004){
        nevZEM1++;  
	hESDword->Fill(3.);
      }
      if((iWord & 0x00000008) == 0x00000008){
        nevZEM2++;  
	hESDword->Fill(4.);
      }
      if((iWord & 0x00000010) == 0x00000010){
        nevZNC++;   
	hESDword->Fill(5.);
      }
      if((iWord & 0x00000020) == 0x00000020){
        nevZPC++;   
	hESDword->Fill(6.);
      }
      //if(centrZNC[0]!=-999. && centrZNC[1]!=-999) centroidZNsideC->Fill(centrZNC[0], centrZNC[1]);
      //if(centrZNA[0]!=-999. && centrZNA[1]!=-999) centroidZNsideA->Fill(centrZNA[0], centrZNA[1]);
      enZNC->Fill(energyZNC);
      enZPC->Fill(energyZPC);
      enZNA->Fill(energyZNA);
      enZPA->Fill(energyZPA);
      enZEM1->Fill(energyZEM1);
      enZEM2->Fill(energyZEM2);
      //
      for(Int_t jj=0; jj<5; jj++){
        if(esdWordCut){
          if((iWord & 0x00000010) == 0x00000010) hZNCTow[jj]->Fill(towZNC[jj]);
          if((iWord & 0x00000020) == 0x00000020) hZPCTow[jj]->Fill(towZPC[jj]);
          if((iWord & 0x00000001) == 0x00000001) hZNATow[jj]->Fill(towZNA[jj]);
          if((iWord & 0x00000002) == 0x00000002) hZPATow[jj]->Fill(towZPA[jj]);
        }
	else{
	  hZNCTow[jj]->Fill(towZNC[jj]);
          hZPCTow[jj]->Fill(towZPC[jj]);
	  hZNATow[jj]->Fill(towZNA[jj]);
          hZPATow[jj]->Fill(towZPA[jj]);
        }
        //
      	if(jj!=0){
          if(esdWordCut){
            if((iWord & 0x00000010) == 0x00000010) sumQznc += towZNC[jj];
	    if((iWord & 0x00000020) == 0x00000020) sumQzpc += towZPC[jj];
            if((iWord & 0x00000001) == 0x00000001) sumQzna += towZNA[jj];
	    if((iWord & 0x00000002) == 0x00000002) sumQzpa += towZPA[jj];
          }
	  else{
	    sumQznc += towZNC[jj];
	    sumQzpc += towZPC[jj];
	    sumQzna += towZNA[jj];
	    sumQzpa += towZPA[jj];
	  }
        }
      }
      //
/*      if(esdWordCut){
        if((iWord & 0x00000010) == 0x00000010) {
	  hSumQZNC->Fill(sumQznc);
	  //
	  //if(centrZNA[0]!=-999. && centrZNA[1]!=-999) centroidZNsideC->Fill(centrZNA[0], centrZNA[1]);
        }
	//
        if((iWord & 0x00000020) == 0x00000020) {
	  hSumQZPC->Fill(sumQzpc);
        }
	
      }
      else{
        hSumQZNC->Fill(sumQznc);
        //
        hSumQZPC->Fill(sumQzpc);
      }
      //
      hSumQZNA->Fill(sumQzna);
      hSumQZPA->Fill(sumQzpa);
*/
    }
    
  }
  
  nAnalyzedFiles++;
  esdFile->Close();
  }//if(nAnalyzedFiles<=nMaxFiles)
  else{
   printf("\t %d files analyzed\n\n",nMaxFiles);
   break;
  }
 } // while closing
  
 printf("    No. of events over threshold: ZNA: %d  ZPA: %d  ZEM1: %d "
        " ZEM2: %d  ZNC: %d  ZPC: %d\n\n", 
        nevZNA, nevZPA, nevZEM1, nevZEM2, nevZNC, nevZPC);
	
 if(plot){  
  //***********************************************************
  // #### ROOT initialization
  gROOT->Reset();
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptTitle(1);
  if(esdWordCut) gStyle->SetOptStat(1111111);
  else gStyle->SetOptStat(1111);
  gStyle->SetOptFit(0);
  gStyle->SetTitleTextColor(4);
  gStyle->SetStatTextColor(4);
  gStyle->SetStatX(0.92);
  gStyle->SetStatY(0.92);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.04);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.16); 
   

  //-------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","ZDCs + ZEMs signals",400,0,500,800);
  c1->Divide(2,3);
  c1->cd(1);
  gPad->SetLogy(1);
  enZNC->Draw("");
  enZNC->SetLineColor(kBlue);
  enZNC->SetFillColor(kBlue);
  c1->cd(2);
  gPad->SetLogy(1);
  enZPC->Draw("");
  enZPC->SetLineColor(kBlue+3);
  enZPC->SetFillColor(kBlue+3);
  c1->cd(3);
  gPad->SetLogy(1);
  enZEM1->SetLineColor(kRed);
  enZEM1->SetFillColor(kRed);
  enZEM1->Draw("");
  c1->cd(4);
  gPad->SetLogy(1);
  enZEM2->SetLineColor(kRed);
  enZEM2->SetFillColor(kRed);
  enZEM2->Draw("");
  c1->cd(5);
  gPad->SetLogy(1);
  enZNA->Draw("");
  enZNA->SetLineColor(kRed);
  enZNA->SetFillColor(kRed);
  c1->cd(6);
  gPad->SetLogy(1);
  enZPA->Draw("");
  enZPA->SetLineColor(kRed+1);
  enZPA->SetFillColor(kRed+1);  
  
  //-------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","Side C ZDCs",0,0,800,400);
  c3->Divide(5,2);
  c3->cd(1);
  gPad->SetLogy(1);
  hZNCTow[0]->SetLineColor(kBlue);
  hZNCTow[0]->SetFillColor(kBlue);
  hZNCTow[0]->Draw("");
  c3->cd(2);
  gPad->SetLogy(1);
  hZNCTow[1]->SetLineColor(kBlue);
  hZNCTow[1]->SetFillColor(kBlue);
  hZNCTow[1]->Draw("");
  c3->cd(3);
  gPad->SetLogy(1);
  hZNCTow[2]->SetLineColor(kBlue);
  hZNCTow[2]->SetFillColor(kBlue);
  hZNCTow[2]->Draw("");
  c3->cd(4);
  gPad->SetLogy(1);
  hZNCTow[3]->SetLineColor(kBlue);
  hZNCTow[3]->SetFillColor(kBlue);
  hZNCTow[3]->Draw("");
  c3->cd(5);
  gPad->SetLogy(1);
  hZNCTow[4]->SetLineColor(kBlue);
  hZNCTow[4]->SetFillColor(kBlue);
  hZNCTow[4]->Draw("");
  //
  c3->cd(6);
  gPad->SetLogy(1);
  hZPCTow[0]->SetLineColor(kBlue+3);
  hZPCTow[0]->SetFillColor(kBlue+3);
  hZPCTow[0]->Draw("");
  c3->cd(7);
  gPad->SetLogy(1);
  hZPCTow[1]->SetLineColor(kBlue+3);
  hZPCTow[1]->SetFillColor(kBlue+3);
  hZPCTow[1]->Draw("");
  c3->cd(8);
  gPad->SetLogy(1);
  hZPCTow[2]->SetLineColor(kBlue+3);
  hZPCTow[2]->SetFillColor(kBlue+3);
  hZPCTow[2]->Draw("");
  c3->cd(9);
  gPad->SetLogy(1);
  hZPCTow[3]->SetLineColor(kBlue+3);
  hZPCTow[3]->SetFillColor(kBlue+3);
  hZPCTow[3]->Draw("");
  c3->cd(10);
  gPad->SetLogy(1);
  hZPCTow[4]->SetLineColor(kBlue+3);
  hZPCTow[4]->SetFillColor(kBlue+3);
  hZPCTow[4]->Draw("");
  
  
  //-------------------------------------------------
  TCanvas *c32 = new TCanvas("c32","side A ZDCs",700,0,800,400);
  c32->Divide(5,2);
  c32->cd(1);
  gPad->SetLogy(1);
  hZNATow[0]->SetLineColor(kRed);
  hZNATow[0]->SetFillColor(kRed);
  hZNATow[0]->Draw("");
  c32->cd(2);
  gPad->SetLogy(1);
  hZNATow[1]->SetLineColor(kRed);
  hZNATow[1]->SetFillColor(kRed);
  hZNATow[1]->Draw("");
  c32->cd(3);
  gPad->SetLogy(1);
  hZNATow[2]->SetLineColor(kRed);
  hZNATow[2]->SetFillColor(kRed);
  hZNATow[2]->Draw("");
  c32->cd(4);
  gPad->SetLogy(1);
  hZNATow[3]->SetLineColor(kRed);
  hZNATow[3]->SetFillColor(kRed);
  hZNATow[3]->Draw("");
  c32->cd(5);
  gPad->SetLogy(1);
  hZNATow[4]->SetLineColor(kRed);
  hZNATow[4]->SetFillColor(kRed);
  hZNATow[4]->Draw("");
  //
  c32->cd(6);
  gPad->SetLogy(1);
  hZPATow[0]->SetLineColor(kRed+1);
  hZPATow[0]->SetFillColor(kRed+1);
  hZPATow[0]->Draw("");
  c32->cd(7);
  gPad->SetLogy(1);
  hZPATow[1]->SetLineColor(kRed+1);
  hZPATow[1]->SetFillColor(kRed+1);
  hZPATow[1]->Draw("");
  c32->cd(8);
  gPad->SetLogy(1);
  hZPATow[2]->SetLineColor(kRed+1);
  hZPATow[2]->SetFillColor(kRed+1);
  hZPATow[2]->Draw("");
  c32->cd(9);
  gPad->SetLogy(1);
  hZPATow[3]->SetLineColor(kRed+1);
  hZPATow[3]->SetFillColor(kRed+1);
  hZPATow[3]->Draw("");
  c32->cd(10);
  gPad->SetLogy(1);
  hZPATow[4]->SetLineColor(kRed+1);
  hZPATow[4]->SetFillColor(kRed+1);
  hZPATow[4]->Draw("");
 }
   
  TFile * fileout = new TFile("ESDhistos.root","recreate");
  fileout->cd();
  //centroidZNsideC->Write();
  //centroidZNsideA->Write();
  enZNC->Write();
  enZNA->Write();
  enZPC->Write();
  enZPA->Write();
  enZEM1->Write();
  enZEM2->Write();
  for(Int_t jj=0; jj<5; jj++){
      hZNCTow[jj]->Write();
      hZPCTow[jj]->Write();
      hZNATow[jj]->Write();
      hZPATow[jj]->Write();
  }
  /*hSumQZNC->Write();
  hSumQZPC->Write();
  hSumQZNA->Write();
  hSumQZPA->Write();*/
  //
  hESDword->Write();
  //
  fileout->Close();

}
