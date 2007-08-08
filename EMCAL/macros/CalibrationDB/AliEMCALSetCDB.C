// Script to create calibration parameters and store them into CDB
// Thre sets of calibration parameters can be created:
// 1) equal parameters
// 2) randomly distributed parameters for decalibrated detector silumations
// 3) gausian  distributed parameters for decalibrated detector silumations
// 
// Modified from PHOS script for EMCAL by Gustavo Conesa
// May 3, 2007: Modified by Aleksei Pavlinov

//.x $ALICE_ROOT/EMCAL/macros/CalibrationDB/AliEMCALSetCDB.C

#if !defined(__CINT__)
#include <TControlBar.h>
#include <TString.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>

#include "AliRun.h"
#include "AliEMCALCalibData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALSetCDB()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
  menu->AddButton("Help to run EMCAL CDB","Help()",
		  "Explains how to use EMCAL CDS menus");
  menu->AddButton("Equal CC","SetCC(0)",
		  "Set equal calibration coefficients");
  menu->AddButton("Random De-calibration","SetCC(1)",
		  "Set random decalibration calibration coefficients");
  menu->AddButton("Gaussian De-calibration","SetCC(2)",
		  "Set gausian decalibration calibration coefficients");
  menu->AddButton("Read equal CC","GetCC(0)",
		  "Read initial equal calibration coefficients");
  menu->AddButton("Read random CC","GetCC(1)",
		  "Read random decalibration calibration coefficients");
  menu->AddButton("Read gaussian CC","GetCC(2)",
		  "Read gausian decalibration calibration coefficients");
  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\nSet calibration parameters and write them into ALICE CDB. Press button \"Equal CC\" to create equal pedestals and gain factors. Press button \"Decalibrate\" to create random pedestals and gain factors to imitate decalibrated detector\n";
  printf(string);
}

//------------------------------------------------------------------------
void SetCC(Int_t flag=0)
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: all calibration coefficients random (decalibration)
  //   flag=2: all calibration coefficients have Gaussian random distribution (decalibration)
  // Author: Boris Polishchuk (Boris.Polichtchouk at cern.ch)

  TString DBFolder;
  Int_t firstRun   =  0; // What is this
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
    firstRun  =  0;
    lastRun   =  0;
    objFormat = "EMCAL initial gain factors and pedestals";
  }
  else if (flag == 1) {
    DBFolder  ="local://DeCalibDB";
    firstRun  =  0;
    lastRun   = 12;
    objFormat = "EMCAL random pedestals and ADC gain factors (12x48x24)";
  }
  else if (flag == 2) {
    DBFolder  ="local://DeCalibDB"; // create directory DeCalibDB in current directory
    firstRun  =  0;
    lastRun   = 12; // Why 12 ?
    objFormat = "EMCAL random pedestals and gausian ADC gain factors (12x48x24)";
  }
  
  AliEMCALCalibData *calibda=new AliEMCALCalibData("EMCAL");
  
  Float_t fADCpedestal = 0.009;
  Float_t fADCchannel  = 0.0153;  // 250 GeV / (16*1024)
  //  Float_t fADCchannel  = 0.00305;
  Float_t rDecalibration  = 0.1 * fADCchannel; // 10% decalibration - just a guess
  Float_t cc=0, ped;

  TRandom rn;
  Int_t nSMod  = 12;
  Int_t nCol   = 48;
  Int_t nRow   = 24;
  Int_t nRow2  = 12; //Modules 11 and 12 are half modules

  for(Int_t supermodule=0; supermodule < nSMod; supermodule++) {
    for(Int_t column=0; column< nCol; column++) {
      if(supermodule >= 10)
 	nRow = nRow2;
      for(Int_t row=0; row< nRow; row++) {
        cc  = fADCchannel;
        ped = fADCpedestal;
	if (flag == 1) {
	  // Decalibration:
	  // Spread calibration coefficients uniformly with
	  // Cmax/Cmin = 5, (Cmax-Cmin)/2 = 0.0015
	  // and pedestals 0.005 +-10%
	  //	  fADCchannel  = rn.Uniform(0.00075,0.00375);
	  cc  = rn.Uniform(0.00140,0.00160);
	  ped = 0;
	} else if (flag == 2) { // Gaussian
	  cc  = rn.Gaus(fADCchannel, rDecalibration);
	  ped = rn.Uniform(0.0045,0.0055);
	}
	calibda->SetADCchannel (supermodule,column,row, cc);
	calibda->SetADCpedestal(supermodule,column,row, ped);
 	cout<<"Set SM: "<<supermodule<<" col "<<column<<" row "<<row
 	    <<" cc "<< cc <<" ped "<<ped<<endl;
      }
    }
  }

  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Aleksei Pavlinov");
  
  AliCDBId id("EMCAL/Calib/Data",firstRun,lastRun); // create in EMCAL/Calib/Data DBFolder 

  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage(DBFolder.Data());
  loc->Put(calibda, id, &md);

}

//------------------------------------------------------------------------
void GetCC(Int_t flag=0)
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: all calibration coefficients random (decalibration)
  // Author: Yuri.Kharlov at cern.ch

  TString DBFolder;
  int drawKey=1;

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
  }
  else if (flag == 1 || flag == 2) {
    //    DBFolder  ="local://DeCalibDB"; // Get DB in current folder
    DBFolder  ="local:///data/r22b/ALICE/PROD/CALIBRATION_May_2007/PI0/DeCalibDB"; // absolute name - Aug 31, 2007
    TString HOST(gSystem->Getenv("HOST"));
    if(HOST.Contains("pc")) { // pdsf; May 31, 2007
      DBFolder  ="local:///eliza5/alice/pavlinov/PROD/CALIBRATION_May_2007/PI0/10GEV/DECALIB/DeCalibDB";
    }
  }

  AliEMCALCalibData* clb  = (AliEMCALCalibData*)
    (AliCDBManager::Instance()
     ->GetStorage(DBFolder.Data())
     ->Get("EMCAL/Calib/Data",
	   gAlice->GetRunNumber())->GetObject());

  static const Int_t nSMod = 12;
  static const Int_t nCol  = 48;
  Int_t nRow  = 24;
  Int_t nRow2 = 12; //Modules 11 and 12 are half modules
  Int_t nCC   = 0;

  TH2F *hPed[nSMod], *hGain[nSMod];
  TH1F *hCCSum = new TH1F("hCCSum"," CC summary (in MeV) ", 200, 0.0, 20.);
  TH1F *hPedSum = new TH1F("hPedSum"," pedestal summary (in MeV) ", 100, 4., 6.);

  TCanvas *cPed=0, *cGain=0, *cPed2=0, *cGain2=0;
  if(drawKey>1) {
    cPed   = new TCanvas("cPed" ,"Pedestals Mod 0-5"   , 10,10,400,800);
    cGain  = new TCanvas("cGain","Gain factors Mod 0-5", 410,10,400,800);
    cPed2  = new TCanvas("cPed2","Pedestals SMod 6-11", 10,10,400,800);
    cGain2 = new TCanvas("cGain2","Gain factors SMod 6-11", 410,10,400,800);
    cPed   ->Divide(2,3);
    cGain  ->Divide(2,3);
    cPed2  ->Divide(2,3);
    cGain2 ->Divide(2,3);
  }
  TCanvas *cSum   = new TCanvas("cSum" ,"summary"   , 10,10,600,800);
  cSum->Divide(1,2); 

  cout<<endl;
  for (Int_t supermodule=0; supermodule<nSMod; supermodule++) {

    if(supermodule >= 10)
      nRow = nRow2;

    TString namePed="hPed";
    namePed+=supermodule;
    TString titlePed="Pedestals in supermodule ";
    titlePed+=supermodule;
    hPed[supermodule] = new TH2F(namePed.Data(),titlePed.Data(),
			    nCol,1.,1.*nCol,nRow,1.,1.*nRow);

    TString nameGain="hGain";
    nameGain+=supermodule;
    TString titleGain="Gain factors in supermodule ";
    titleGain+=supermodule;
    hGain[supermodule] = new TH2F(nameGain.Data(),titleGain.Data(),
				    nCol,1.,1.*nCol,nRow,1.,1.*nRow);
    for (Int_t column=0; column<nCol; column++) {
      for (Int_t row=0; row<nRow; row++) {
	Float_t ped  = clb->GetADCpedestal(supermodule,column,row);
	Float_t gain = clb->GetADCchannel (supermodule,column,row);
	//cout<<"Get SM: "<<supermodule<<" col "<<column<<" row "<<row
	//<<" chan "<<gain<<endl;
	hPed[supermodule] ->SetBinContent(column+1,row+1,ped*1.e+3);  // in mev
	hGain[supermodule]->SetBinContent(column+1,row+1,gain*1.e+3); // in mev

        hPedSum->Fill(ped*1.e+3);
        hCCSum->Fill(gain*1.e+3); 

        nCC++;
      }
    }
    cout<<" Fill cc for SM "<< supermodule << " nCC "<< nCC << endl;
    
    if(drawKey>1) {
      if(supermodule < 6){
        cPed ->cd(supermodule+1);
        hPed[supermodule]->Draw("lego2");
        cGain->cd(supermodule+1);
        hGain[supermodule]->Draw("lego2");
      }
      else{
        cPed2 ->cd(supermodule-5);
        hPed[supermodule]->Draw("lego2");
        cGain2->cd(supermodule-5);
        hGain[supermodule]->Draw("lego2");
      }
    }

  }
  cout << " Get "<<nCC<<" calibration coeffs"<<endl;

  cSum->cd(1);
  gStyle->SetOptFit(111);
  hCCSum->Fit("gaus");
  hCCSum->SetLineWidth(2);
  hCCSum->GetFunction("gaus")->SetLineColor(2);
  
  cSum->cd(2);
  hPedSum->Draw();
  hPedSum->SetLineWidth(2);

  cSum->Update();

  /*
  cPed   ->Print("pedestals_SM_0_6.eps");
  cGain  ->Print("gains_SM_0_5.eps");
  cPed2  ->Print("pedestals_SM_6_11.eps");
  cGain2 ->Print("gains_SM_6_11.eps");
  */
}
