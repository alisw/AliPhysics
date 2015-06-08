#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <TRandom.h>
#include <TSystem.h>
#include <TDatime.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>

#include <AliCDBManager.h>
#include <AliTRDcalibDB.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>
#include <AliGeometry.h>
#include <AliPID.h>

#include "../TRD/AliTRDgeometry.h"
#include "../TRD/Cal/AliTRDCalROC.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalDet.h"

#endif

Int_t GetPlane(Int_t d);
Int_t GetSector(Int_t d);
Int_t GetChamber(Int_t d);
void ReadCDD(Bool_t residual = kFALSE, Int_t ndet = 319);

void ReadCDD(Bool_t residual, Int_t ndet){


  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
 

  //definition number of detector et histos
  Int_t Ndet = AliTRDgeometry::kNdet;
  TH1F *Vdriftdet   = new TH1F("Vdriftdet","Vdrift det",Ndet,0,Ndet);
  TH1F *Gaindet   = new TH1F("Gaindet","Gain det",Ndet,0,Ndet);
  TH1F *Vdriftdets   = new TH1F("Vdriftdets","Vdrift",20,1.0,2.0);
  TH1F *Gaindets   = new TH1F("Gaindets","Gain",15,0.0,2.0);


  TH1F *T0pad   = new TH1F("T0pad","T0pad",900,-1.0,1.0);
  TH1F *Vdriftpad   = new TH1F("Vdriftpad","Vdriftpad",100,0.9,2.1);
  TH1F *Gainpad   = new TH1F("Gainpad","Gainpad",150,0.2,1.8);

  //////////////////
  //Datebase used 
  //////////////////
  AliCDBManager *man = AliCDBManager::Instance();
 

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("TRD/Calib/LocalGainFactor","local://.");
  man->SetSpecificStorage("TRD/Calib/LocalT0","local://.");
  man->SetSpecificStorage("TRD/Calib/LocalVdrift","local://.");
  man->SetSpecificStorage("TRD/Calib/ChamberT0","local://.");
  if(! residual){
    man->SetSpecificStorage("TRD/Calib/ChamberVdrift","local://.");
    man->SetSpecificStorage("TRD/Calib/ChamberGainFactor","local://.");
  }
  man->SetRun(0);
 
  
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDgeometry *geo = new AliTRDgeometry();

   
  ////////////////////////////////////////////
  // drift velocity and gain for detectors
  ///////////////////////////////////////////

  //drift velocity par detector
  for(Int_t k = 0; k < Ndet; k++){
  Vdriftdet->Fill(k+0.5,cal->GetVdriftAverage(k));
  Vdriftdets->Fill(cal->GetVdriftAverage(k));
  //cout << "Vdrift average: " << cal->GetVdriftAverage(k) << endl;
  }

  TCanvas *c1 = new TCanvas("c1","Calibchamberdrift",50,50,600,800);
  c1->Divide(1,2);
  c1->cd(1);
  Vdriftdet->SetStats(0);
  Vdriftdet->SetXTitle("Det Number");
  Vdriftdet->SetYTitle("Mean Drift Velocity [cm/#mu s]");
  Vdriftdet->Draw();

  //gain par detector
  for(Int_t k = 0; k < Ndet; k++){
  Gaindet->Fill(k+0.5,cal->GetGainFactorAverage(k));
  Gaindets->Fill(cal->GetGainFactorAverage(k));
  //cout << "Vdrift average: " << cal->GetVdriftAverage(k) << endl;
  }

  //TCanvas *c1g = new TCanvas("c1g","Calibchambergain",50,50,600,800);
  c1->cd(2);
  Gaindet->SetStats(0);
  Gaindet->SetXTitle("Det Number");
  Gaindet->SetYTitle("Mean Gain");
  Gaindet->Draw();


  TCanvas *c1u = new TCanvas("c1u","Calibchamberdriftu",50,50,600,800);
  c1u->Divide(1,2);
  c1u->cd(1);
  Vdriftdets->SetStats(0);
  Vdriftdets->SetXTitle("Mean Drift Velocity [cm/#mu s]");
  Vdriftdets->SetYTitle("Number of detectors");
  Vdriftdets->Draw();

  c1u->cd(2);
  Gaindets->SetStats(0);
  Gaindets->SetXTitle("Mean Gain");
  Gaindets->SetYTitle("Number of detectors");
  Gaindets->Draw();

  /////////////////////////////////
  // Maps of one detector choosen
  /////////////////////////////////

  Int_t rowMaxd = geo->GetRowMax(GetPlane(ndet),GetChamber(ndet),GetSector(ndet));
  Int_t colMaxd = geo->GetColMax(GetPlane(ndet));
    
  TH2F *hZYvd     = new TH2F("hZYvd"   ,"Y vs Z", rowMaxd,0,rowMaxd,colMaxd,0,colMaxd);
  TH2F *hZYg     = new TH2F("hZYg"   ,"Y vs Z", rowMaxd,0,rowMaxd,colMaxd,0,colMaxd);
  TH2F *hZYT0     = new TH2F("hZYT0"   ,"Y vs Z", rowMaxd,0,rowMaxd,colMaxd,0,colMaxd);
 
  
  for (Int_t  col = 0;  col <  colMaxd;  col++) {
    for (Int_t  row = 0;  row <  rowMaxd;  row++) {
      hZYvd->Fill(row,col,cal->GetVdrift(ndet,col,row));
      hZYg->Fill(row,col,cal->GetGainFactor(ndet,col,row));
      hZYT0->Fill(row,col,cal->GetT0(ndet,col,row));
    }//boucle row
  }//boucle col

  ////////////////////////////////
  // coefficients per pad
  ///////////////////////////////

  for(Int_t isector = 0; isector < 18; isector++){
    for(Int_t iplane =0; iplane < 6; iplane++){
      for(Int_t ichamb =0; ichamb < 5; ichamb++){
	
	Int_t det = AliTRDgeometry::GetDetector(iplane,ichamb,isector);
	Int_t rowMax = geo->GetRowMax(iplane,ichamb,isector);
	Int_t colMax = geo->GetColMax(iplane);
	
	for (Int_t  col = 0;  col <  colMax;  col++) {
	  for (Int_t  row = 0;  row <  rowMax;  row++) {
	    Gainpad->Fill(cal->GetGainFactor(det,col,row));
	    T0pad->Fill(cal->GetT0(det,col,row));
	    Vdriftpad->Fill(cal->GetVdrift(det,col,row));
	  }//boucle row
	}//boucle col
	
      }//loop chamber
    }//loop plane  
  }

  /////////////////////////////
  // Plot the choosen detector
  /////////////////////////////
    
  TCanvas *c2vd = new TCanvas("c2vd","choosen detector",50,50,600,800);
  c2vd->Divide(3,1);
  c2vd->cd(1);
  hZYvd->SetXTitle("row Number");
  hZYvd->SetYTitle("col Number");
  hZYvd->SetZTitle("Drift velocity");
  hZYvd->SetStats(0);
  hZYvd->Draw("LEGO");

  c2vd->cd(2);
  hZYT0->SetXTitle("row Number");
  hZYT0->SetYTitle("col Number");
  hZYT0->SetZTitle("T0 [time bin]");
  hZYT0->SetStats(0);
  hZYT0->Draw("SURF1");
  
  c2vd->cd(3);
  hZYg->SetXTitle("row Number");
  hZYg->SetYTitle("col Number");
  hZYg->SetZTitle("gain factor");
  hZYg->SetStats(0);
  hZYg->Draw("SURF1");

  /////////////////////////////
  // Plot the pad info
  /////////////////////////////

  TCanvas *c5vd = new TCanvas("c5vd","PadInfo",50,50,600,800);
  c5vd->Divide(3,1);
  c5vd->cd(1);
  Vdriftpad->SetXTitle("v_{d} [cm/#mus]");
  Vdriftpad->SetYTitle("number of pads");
  Vdriftpad->SetStats(0);
  Vdriftpad->Draw();

  c5vd->cd(2);
  T0pad->SetXTitle("T0 [time bin]");
  T0pad->SetYTitle("number of pads");
  T0pad->SetStats(0);
  T0pad->Draw();
  

  c5vd->cd(3);
  Gainpad->SetXTitle("gain factor");
  Gainpad->SetYTitle("number of pads");
  Gainpad->SetStats(0);
  Gainpad->Draw();


}
//____________________________________________________________________________________________
Int_t GetPlane(Int_t d) 
{
  //
  // Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t GetChamber(Int_t d) 
{
  //
  // Reconstruct the chamber number from the detector number
  //
  Int_t fgkNplan = 6;

  return ((Int_t) (d % 30) / fgkNplan);

}

//_____________________________________________________________________________
Int_t GetSector(Int_t d) 
{
  //
  // Reconstruct the sector number from the detector number
  //
  Int_t fg = 30;

  return ((Int_t) (d / fg));

}
