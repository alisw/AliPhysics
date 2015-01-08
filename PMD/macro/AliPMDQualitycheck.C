//********************************************************************
//     Macro that checks PMD data Quality 
//                                         using the ESD classes
//     Modified by P.V.K.S.BABA (J.U.)/Ajay(IOP) for doing Quality check.
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"

#include "AliESD.h"
#endif
double Xcon[96] =
  { 75.133, 54.204, 53.254, 32.326, 31.376,10.447,
    75.133, 54.204, 53.254, 32.326, 31.376,10.447,
    75.133, 54.204, 53.254, 32.326, 31.376,10.447,
    75.133, 54.204, 53.254, 32.326, 31.376,10.447,
    -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
    -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
    -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
    -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
    9.167, -32.543, -33.493, -75.133,
    9.167, -32.543, -33.493, -75.133,
    9.167, -32.543, -33.493, -75.133,
    9.167, -32.543, -33.493, -75.133,
    9.167, -32.543, -33.493, -75.133,
    9.167, -32.543, -33.493, -75.133,
    -9.167, 32.543, 33.493, 75.133,
    -9.167, 32.543, 33.493, 75.133,
    -9.167, 32.543, 33.493, 75.133,
    -9.167, 32.543, 33.493, 75.133,
    -9.167, 32.543, 33.493, 75.133,
    -9.167, 32.543, 33.493, 75.133,};

double Ycon[96] =
  {86.475,  86.475,  86.475, 86.475,  86.475,  86.475,
   38.225,  38.225,  38.225, 38.225,  38.225,  38.225,
   37.325,  37.325,  37.325, 37.325,  37.325,  37.325,
   -10.925, -10.925, -10.925, -10.925, -10.925, -10.925,
   -86.475, -86.475, -86.475, -86.475, -86.475, -86.475,
   -38.225,  -38.225,  -38.225, -38.225, -38.225, -38.225,
   -37.325,  -37.325,  -37.325, -37.325,  -37.325,  -37.325
   10.925, 10.925, 10.925, 10.925, 10.925, 10.925,

   86.475,  86.475, 86.475,  86.475,
   62.225,  62.225, 62.225,  62.225,
   61.325,  61.325, 61.325,  61.325,
   37.075, 37.075, 37.075, 37.075
   36.175,  36.175, 36.175,  36.175,
   11.925, 11.925, 11.925 , 11.925,
   -86.475,  -86.475, -86.475,  -86.475,
   -62.225,  -62.225, -62.225,  -62.225,
   -61.325,  -61.325, -61.325,  -61.325,
   -37.075,  -37.075, -37.075,  -37.075
   -36.175,  -36.175, -36.175,  -36.175,
   -11.925, -11.925,  -11.925 , -11.925 };


extern TStyle *gStyle;

// Default is that it will take 1000 Events sample. Want to process less give (Number) .
Int_t AliPMDQualitycheck(Int_t nevt=1000) { 
AliCDBManager::Instance()->SetRun(0);
  TStopwatch timer;
  char PRINTSMN ;

//  printf("Whether You want plotof All SMN , Give y or n ? \n " ) ;
//  scanf("%1s" , &PRINTSMN );
 
  TH2F *hP1 = new TH2F("hP1","XY of Clusters",100,-100.,100.,100,-100.,100.);
  TH1F *hC2 = new TH1F("hC2","CPV  PHI",200,-1,9);
  TH1F *hP2 = new TH1F("hP2","PRE  PHI",200,-1,9);
  TH1F *hC3 = new TH1F("hC3","CPV  Clus",30,0.,500.);
  TH1F *hP3 = new TH1F("hP3","PRE  N-gammalike",20,0.,500.);
  TH1F *hP4 = new TH1F("hP4","PRE  EDEP",30,0.,1000.);
  TH1F *hC5 = new TH1F("hC5","CPV  n-cell",20,0.,100.);
  TH1F *hP5 = new TH1F("hP5","PMD  n-cell",20,0.,100.);
  TH2F *hCP0 = new TH2F("hCP0","PRE CLUS Quad.1 vs 2",150,0.,300.,150,0.,300.);
  TH2F *hCP1 = new TH2F("hCP1","PRE CLUS Quad.3 vs 4",150,0.,300.,150,0.,300.);
  TH2F *hCP2 = new TH2F("hCP2","PRE EDEP Quad.3 vs 4",50,0.,300.,50,0.,300.);
  TH2F *hCP3 = new TH2F("hCP3","PRE EDEP vs Tot Clus ",10,0.,1000.,10,0.,300.);
  TH2F *hCP4 = new TH2F("hCP4","PRE Clus vs CPV Clus ",150,0.,200.,150,0.,200.);

  TH2F *hSM1 = new TH2F("hSM1","PRE Cluster XY",200,-100,100,200,-100,100);
  TH2F *hSM2 = new TH2F("hSM2","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM3 = new TH2F("hSM3","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM4 = new TH2F("hSM4","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM5 = new TH2F("hSM5","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM6 = new TH2F("hSM6","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM7 = new TH2F("hSM7","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM8 = new TH2F("hSM8","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM9 = new TH2F("hSM9","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM10 = new TH2F("hSM10","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM11 = new TH2F("hSM11","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM12 = new TH2F("hSM12","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM13 = new TH2F("hSM13","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM14 = new TH2F("hSM14","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM15 = new TH2F("hSM15","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM16 = new TH2F("hSM16","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM17 = new TH2F("hSM17","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM18 = new TH2F("hSM18","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM19 = new TH2F("hSM19","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM20 = new TH2F("hSM20","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM21 = new TH2F("hSM21","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM22 = new TH2F("hSM22","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM23 = new TH2F("hSM23","",999,-100.0,100.0,999,-100.0,100.0);
  TH2F *hSM24 = new TH2F("hSM24","",999,-100.0,100.0,999,-100.0,100.0);
  TH1F *hSM = new TH1F("hSM","Plot of all 24 Super Modules",24,0,24);

// Star processing.
  
  AliPMDUtility *cc = new AliPMDUtility(); 
 
  TFile *ef=TFile::Open("AliESDtree.root");
  if (!ef || !ef->IsOpen()) {cerr<<"Cant open AliESDtree.root !\n"; return 1;}
  AliESD* event = new AliESD;
  TTree* tree = (TTree*) ef->Get("esdTree");
  if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
  tree->SetBranchAddress("ESD", &event);
  Int_t SMN;
  Int_t n=0;
  Float_t TotCPVClus     ;
  Float_t TotPREClus     ;
  Float_t TotPREEdep  ;
  Float_t TotCPVCell     ;
  Float_t TotPRECell     ;
  Float_t PREcluQUAD[4]  ;
  Float_t CPVcluQUAD[4]  ;
  Float_t PREadcQUAD[4]  ;
  Float_t CPVADCQUAD[4]  ;
  Float_t PREcelQUAD[4]  ;
  Float_t CPVcelQUAD[4]  ;
 
 
  //******* The loop over events
  while (tree->GetEvent(n)) {
    cout<<endl<<"Processing event number : "<<n++<<endl;
    
    
    Int_t CPVhits = 0 ;
    Int_t PMDclus = 0 ;
    Int_t npmdcl=event->GetNumberOfPmdTracks();
    //cout<<"Number of PMD tracks : "<<npmdcl<<endl; 
    
    //****** The loop over PMD clusters

    for (Int_t kk=0; kk<4;kk++) {
        CPVcluQUAD[kk] = 0.0 ;
        PREcluQUAD[kk] = 0.0 ;
        CPVcelQUAD[kk] = 0.0 ;
        PREcelQUAD[kk] = 0.0 ;
        PREadcQUAD[kk] = 0.0 ;
     }
    while (npmdcl--) {
      AliESDPmdTrack *pmdtr = event->GetPmdTrack(npmdcl);
      
      Int_t   det   = pmdtr->GetDetector();
      Float_t clsX  = pmdtr->GetClusterX();
      Float_t clsY  = pmdtr->GetClusterY();
      Float_t clsZ  = pmdtr->GetClusterZ();
      Float_t ncell = pmdtr->GetClusterCells();
      Float_t adc   = pmdtr->GetClusterADC();
      Float_t pid   = pmdtr->GetClusterPID();

//    cout << " "  << det << " " << ncell << " " << adc << " " << pid << endl ;      
      cc->SetXYZ(clsX,clsY,clsZ);
      //cc->SetPxPyPz(clsX,clsY,clsZ);
      cc->CalculateEta();
      cc->CalculatePhi();
      Float_t eta = cc->GetEta();
      Float_t phi = cc->GetPhi();
 
// Calculating S.Module Number from Cluster .
            CalculateSMN(clsX, clsY,SMN);
      if( det == 1)
	{
          if(SMN >= 0 && SMN <= 5) {
                 ++CPVcluQUAD[0] ;
           CPVADCQUAD[0] =+ adc  ;
           CPVcelQUAD[0] =+ ncell ;
          }
          if(SMN >= 6 && SMN <=11) {
                 ++CPVcluQUAD[1] ;
           CPVADCQUAD[1] =+ adc  ;
           CPVcelQUAD[1] =+ ncell ;
          }
          if(SMN >=12 && SMN <=17) {
                 ++CPVcluQUAD[2] ;
           CPVADCQUAD[2] =+ adc  ;
           CPVcelQUAD[2] =+ ncell ;
          }
          if(SMN >=18 && SMN <=23) {
                 ++CPVcluQUAD[3] ;
           CPVADCQUAD[3] =+ adc  ;
           CPVcelQUAD[3] =+ ncell ;
          }

	  if(eta >= 2.3 && eta <= 3.5)
	  {
	    hC2->Fill(phi);
	  }
	}
      if( det == 0)
	{
          if(SMN >= 0 && SMN <= 5) { 
                 ++PREcluQUAD[0] ;
           PREadcQUAD[0] =+ adc  ;    
           PREcelQUAD[0] =+ ncell ;    
          }    
          if(SMN >= 6 && SMN <=11) { 
                 ++PREcluQUAD[1] ;
           PREadcQUAD[1] =+ adc  ;    
           PREcelQUAD[1] =+ ncell ;    
          }    
          if(SMN >=12 && SMN <=17) { 
                 ++PREcluQUAD[2] ;
           PREadcQUAD[2] =+ adc  ;    
           PREcelQUAD[2] =+ ncell ;    
          }    
          if(SMN >=18 && SMN <=23) { 
                 ++PREcluQUAD[3] ;
           PREadcQUAD[3] =+ adc  ;    
           PREcelQUAD[3] =+ ncell ;    
          }    
          if ( n <= 100 ) { 
                         hSM->Fill(SMN);
            if(SMN == 0)hSM1->Fill(-clsX,clsY);
            if(SMN == 0)hSM1->Fill(-clsX,clsY);
            if(SMN == 1)hSM2->Fill(-clsX,clsY);
            if(SMN == 2)hSM3->Fill(-clsX,clsY);
            if(SMN == 3)hSM4->Fill(-clsX,clsY);
            if(SMN == 4)hSM5->Fill(-clsX,clsY);
            if(SMN == 5)hSM6->Fill(-clsX,clsY);
            if(SMN == 6)hSM7->Fill(-clsX,clsY);
            if(SMN == 7)hSM8->Fill(-clsX,clsY);
            if(SMN == 8)hSM9->Fill(-clsX,clsY);
            if(SMN == 9)hSM10->Fill(-clsX,clsY);
            if(SMN ==10)hSM11->Fill(-clsX,clsY);
            if(SMN ==11)hSM12->Fill(-clsX,clsY);
            if(SMN ==12)hSM13->Fill(-clsX,clsY);
            if(SMN ==13)hSM14->Fill(-clsX,clsY);
            if(SMN ==14)hSM15->Fill(-clsX,clsY);
            if(SMN ==15)hSM16->Fill(-clsX,clsY);
            if(SMN ==16)hSM17->Fill(-clsX,clsY);
            if(SMN ==17)hSM18->Fill(-clsX,clsY);
            if(SMN ==18)hSM19->Fill(-clsX,clsY);
            if(SMN ==19)hSM20->Fill(-clsX,clsY);
            if(SMN ==20)hSM21->Fill(-clsX,clsY);
            if(SMN ==21)hSM22->Fill(-clsX,clsY);
            if(SMN ==22)hSM23->Fill(-clsX,clsY);
            if(SMN ==23)hSM24->Fill(-clsX,clsY);
          }     
	  if(eta >= 2.3 && eta <= 3.5)
	  {
	    hP2->Fill(phi);
	  }
	    hP1->Fill(clsX,clsY);
	}
    }
    for (Int_t k=0; k<4;k++) {
        TotCPVClus =+ CPVcluQUAD [k] ;
        TotPREClus =+ PREcluQUAD [k] ;
        TotCPVCell =+ CPVcelQUAD [k] ;
        TotPRECell =+ PREcelQUAD [k] ;
        TotPREEdep =+ PREadcQUAD [k] ;     
    }
        Float_t TotCPVPREClus = TotPREClus + TotCPVClus ;
        Float_t TotCPVPRECell = TotPRECell + TotCPVCell ;
	  if(eta >= 2.3 && eta <= 3.5)
	  {
	    hC3->Fill(TotCPVClus);
	    hP3->Fill(TotPREClus);
	    hP4->Fill(TotPREEdep);
	    hP5->Fill(TotPRECell);
           hCP0->Fill(PREcluQUAD[0],PREcluQUAD[1]);
           hCP1->Fill(PREcluQUAD[2],PREcluQUAD[3]);
           hCP2->Fill(PREadcQUAD[2],PREadcQUAD[3]);
           hCP3->Fill(TotPREEdep,TotCPVPREClus);
           hCP4->Fill(TotPREClus,TotCPVClus);
	  }
           TotCPVClus == 0.0; 
           TotPREClus == 0.0; 
           TotCPVCell == 0.0; 
           TotPRECell == 0.0; 
           TotPREEdep == 0.0; 
    if(n >= nevt)break;    
  }

  gStyle->SetOptStat(110000);
  gStyle->SetOptFit(1);

//  if ( PRINTSMN == 'y' ) {
  TCanvas *cP1 = new TCanvas("cP1","Cluster XY", 10,10, 600, 600);
  cP1->Range(-100, -100,100 ,100 );
  hSM1->SetMarkerColor(2);
  hSM1->Draw();
  hSM1->GetXaxis()->SetTitle("Cluster X");
  hSM1->GetYaxis()->SetTitle("Cluster Y");
  hSM2->SetMarkerColor(2);
  hSM2->Draw("same");
  hSM3->SetMarkerColor(2);
  hSM3->Draw("same");
  hSM4->SetMarkerColor(2);
  hSM4->Draw("same");
  hSM5->SetMarkerColor(2);
  hSM5->Draw("same");
  hSM6->SetMarkerColor(2);
  hSM6->Draw("same");
  hSM7->SetMarkerColor(4);
  hSM7->Draw("same");
  hSM8->SetMarkerColor(4);
  hSM8->Draw("same");
  hSM9->SetMarkerColor(4);
  hSM9->Draw("same");
  hSM10->SetMarkerColor(4);
  hSM10->Draw("same");
  hSM11->SetMarkerColor(4);
  hSM11->Draw("same");
  hSM12->SetMarkerColor(4);
  hSM12->Draw("same");
  hSM13->SetMarkerColor(6);
  hSM13->Draw("same");
  hSM14->SetMarkerColor(6);
  hSM14->Draw("same");
  hSM15->SetMarkerColor(6);
  hSM15->Draw("same");
  hSM16->SetMarkerColor(6);
  hSM16->Draw("same");
  hSM17->SetMarkerColor(6);
  hSM17->Draw("same");
  hSM18->SetMarkerColor(6);
  hSM18->Draw("same");
  hSM19->SetMarkerColor(8);
  hSM19->Draw("same");
  hSM20->SetMarkerColor(8);
  hSM20->Draw("same");
  hSM21->SetMarkerColor(8);
  hSM21->Draw("same");
  hSM22->SetMarkerColor(8);
  hSM22->Draw("same");
  hSM23->SetMarkerColor(8);
  hSM23->Draw("same");
  hSM24->SetMarkerColor(8);
  hSM24->Draw("same");

  Int_t linsav = gStyle->GetLineWidth();
  DrawPMDBoundary();
  DrawPMDBoundarySM1();
  DrawPMDBoundarySM2();
  DrawPMDBoundarySM3();
  DrawPMDBoundarySM4();
  cP1->Print("ClusterXY.gif");
  

  TCanvas *cP2 = new TCanvas("PHI CPV / PMD","",10,10, 600,600);
  cP2->Divide(1,2);
  cP2->cd(1);
  cP2->SetFillColor(0);
  hC2->SetLineColor(4);
  hC2->Draw();
  cP2->cd(2);
  hP2->SetLineColor(2);
  hP2->Draw();
  cP2->Print("CPVPREphi.gif");

// }

  TCanvas *cP2 = new TCanvas("cP2","",10,10,600,600);
  cP2->cd();
  hSM->SetFillColor(2);
  hSM->Draw();
  cP2->Print("AllSMN.gif");



  TCanvas *cP3 = new TCanvas("CPV Clus  PRE Clus Correlations", " ",10,10, 600,600);
  cP3->Divide(2,2);
  cP3->cd(1);
  hCP0->SetMarkerColor(9);
  hCP0->Draw();
  cP3->cd(2);
  hCP1->SetMarkerColor(6);
  hCP1->Draw();
  cP3->cd(3);
  hP3->SetLineColor(2);
  hP3->Draw();
  cP3->cd(4);
  hCP4->SetMarkerColor(3);
  hCP4->Draw();
  cP3->Print("CPVPREClus.gif");

  TCanvas *cP6 = new TCanvas("CPV Clus PRE Adc"," ",10,10, 600,600);
  cP6->Divide(1,2);
  cP6->cd(1);
  hC3->SetLineColor(4);
  hC3->Draw();
  cP6->cd(2);
  hP4->SetLineColor(2);
  hP4->Draw();
  cP6->Print("CPVPREAdc.gif");

  timer.Stop();
  timer.Print();
  return 0;
  }

void CalculateSMN( Float_t clsX, Float_t clsY, Int_t & smn)
{
  //smn = 0;

  //---------------------------------------------------------------------
  if((clsX <= Xcon[0]) && (clsX >= Xcon[1]) &&
    (clsY <= Ycon[0]) && (clsY >= Ycon[6])) smn = 0 ;

  else if((clsX <=Xcon[2]) && (clsX >= Xcon[3]) &&
  (clsY <= Ycon[1]) && (clsY >= Ycon[7]))smn = 1 ;

  else if((clsX <=Xcon[4]) && (clsX >= Xcon[5]) &&
  (clsY <= Ycon[3]) && (clsY >= Ycon[8]))smn = 2 ;

  else if((clsX <= Xcon[0]) && (clsX >= Xcon[1]) &&
  (clsY <= Ycon[12]) && (clsY >= Ycon[18])) smn = 3 ;

  else if((clsX <=Xcon[2]) && (clsX >= Xcon[3]) &&
  (clsY <= Ycon[12]) && (clsY >= Ycon[18]))smn = 4 ;

  else if((clsX <=Xcon[4]) && (clsX >= Xcon[5]) &&
  (clsY <= Ycon[12]) && (clsY >= Ycon[18]))smn = 5 ;
  //------------------------------------------------------------------
  else if((clsX >= Xcon[24]) && (clsX <= Xcon[25]) &&
  (clsY >= Ycon[24]) && (clsY <= Ycon[30])) smn = 6 ;

  else if((clsX >=Xcon[26]) && (clsX <= Xcon[27]) &&
  (clsY >= Ycon[25]) && (clsY <= Ycon[31]))smn = 7 ;

  else if((clsX >=Xcon[28]) && (clsX <= Xcon[29]) &&
  (clsY >= Ycon[26]) && (clsY <= Ycon[32]))smn = 8 ;

  else if((clsX >= Xcon[24]) && (clsX <= Xcon[25]) &&
  (clsY >= Ycon[36]) && (clsY <= Ycon[42])) smn = 9 ;

  else if((clsX >=Xcon[26]) && (clsX <= Xcon[27]) &&
  (clsY >= Ycon[36]) && (clsY <= Ycon[42]))smn = 10;

  else if((clsX >=Xcon[28]) && (clsX <= Xcon[29]) &&
  (clsY >= Ycon[36]) && (clsY <= Ycon[42]))smn = 11;
  //------------------------------------------------------------------
  else if((clsX <= Xcon[48]) && (clsX >= Xcon[49]) &&
  (clsY <= Ycon[48]) && (clsY >= Ycon[52])) smn = 12 ;

  else if((clsX <=Xcon[50]) && (clsX >= Xcon[51]) &&
  (clsY <= Ycon[48]) && (clsY >= Ycon[52]))smn = 13 ;

  else if((clsX <=Xcon[48]) && (clsX >= Xcon[49]) &&
  (clsY <= Ycon[56]) && (clsY >= Ycon[60]))smn = 14 ;

  else if((clsX <=Xcon[50]) && (clsX >= Xcon[51]) &&
  (clsY <= Ycon[56]) && (clsY >= Ycon[60]))smn = 15 ;

  else if((clsX <=Xcon[48]) && (clsX >= Xcon[49]) &&
  (clsY <= Ycon[64]) && (clsY >= Ycon[68]))smn = 16 ;

  else if((clsX <=Xcon[50]) && (clsX >= Xcon[51]) &&
  (clsY <= Ycon[64]) && (clsY >= Ycon[68]))smn = 17 ;
  //--------------------------------------------------------------
  else if((clsX >= Xcon[72]) && (clsX <= Xcon[73]) &&
  (clsY >= Ycon[72]) && (clsY <= Ycon[76])) smn = 18 ;

  else if((clsX >=Xcon[74]) && (clsX <= Xcon[75]) &&
  (clsY >= Ycon[72]) && (clsY <= Ycon[76]))smn = 19 ;

  else if((clsX >=Xcon[72]) && (clsX <= Xcon[73]) &&
  (clsY >= Ycon[80]) && (clsY <= Ycon[84]))smn = 20 ;

  else if((clsX >=Xcon[74]) && (clsX <= Xcon[75]) &&
  (clsY >= Ycon[80]) && (clsY <= Ycon[84]))smn = 21;

  else if((clsX >= Xcon[72]) && (clsX <= Xcon[73]) &&
  (clsY >= Ycon[88]) && (clsY <= Ycon[92])) smn = 22 ;

  else if((clsX >=Xcon[74]) && (clsX <= Xcon[75]) &&
  (clsY >= Ycon[88]) && (clsY <= Ycon[92]))smn = 23 ;
  else smn = 111;

 }

void DrawPMDBoundary()
{
  TH2F *h = new TH2F("h","",200,-100,100,200,-100,100);
  gStyle->SetLineWidth(2);
  gStyle->SetLineColor(2);
  TLine * l;
  l = new TLine(75.1333, 86.475, -75.1333, 86.475); l->Draw("same");
  l = new TLine(-75.1333, 86.470,-75.1333, -86.475); l->Draw("same");
  l = new TLine(-75.1333, -86.475,75.1333, -86.475); l->Draw("same");
  l = new TLine(75.1333, -86.475,75.1333, 86.475); l->Draw("same");
}
void DrawPMDBoundarySM1()
{
  TH2F *hsm1 = new TH2F("hsm1","",200,-100,100,200,-100,100);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(4);
  TLine * l;
  l = new TLine(-75.1333, 86.475, -10.447,  86.475); l->Draw("same");
  l = new TLine(-10.447, 86.475, -10.446, -10.925); l->Draw("same");
  l = new TLine(-10.446, -10.925, -75.1333,-10.925); l->Draw("same");
  l = new TLine(-75.1333,-10.925, -75.1333, 86.475); l->Draw("same");
}
void DrawPMDBoundarySM2()
{
  TH2F *hsm2 = new TH2F("hsm2","",200,-100,100,200,-100,100);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(4);
  TLine * l;
  l = new TLine(75.1333, -86.475, 10.446,  -86.475); l->Draw("same");
  l = new TLine(10.446,  -86.475, 10.446,  10.925); l->Draw("same");
  l = new TLine(10.446,   10.925, 75.1333, 10.925); l->Draw("same");
  l = new TLine(75.1333,  10.925, 75.1333, -86.475); l->Draw("same");
}

void DrawPMDBoundarySM3()
{
  TH2F *hsm3 = new TH2F("hsm3","",200,-100,100,200,-100,100);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  TLine * l;
  l = new TLine(  -9.167, 86.475, 75.1333, 86.475); l->Draw("same");
  l = new TLine(75.1333,86.475, 75.1333, 11.925); l->Draw("same");
  l = new TLine(75.1333,11.925,   -9.167,  11.925); l->Draw("same");
  l = new TLine(  -9.167, 11.925,   -9.167,  86.475); l->Draw("same");
}
void DrawPMDBoundarySM4()
{
  TH2F *hsm4 = new TH2F("hsm4","",200,-100,100,200,-100,100);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  TLine * l;
  l = new TLine(9.167, -86.475, -75.1333,-86.475); l->Draw("same");
  l = new TLine(-75.1333,-86.475, -75.1333,-11.925); l->Draw("same");
  l = new TLine(-75.1333,-11.925, 9.167, -11.925); l->Draw("same");
  l = new TLine(9.167, -11.925, 9.167, -86.475); l->Draw("same");
}

