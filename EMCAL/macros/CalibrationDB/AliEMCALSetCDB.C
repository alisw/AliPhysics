
// Script to create calibration parameters and store them into CDB
// Two sets of calibration parameters can be created:
// 1) equal parameters
// 2) randomly distributed parameters for decalibrated detector silumations
// Modified from PHOS script for EMCAL by Gustavo Conesa

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TCanvas.h"

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
  menu->AddButton("Decalibrate","SetCC(1)",
		  "Set random decalibration calibration coefficients");
  menu->AddButton("Read equal CC","GetCC(0)",
		  "Read initial equal calibration coefficients");
  menu->AddButton("Read random CC","GetCC(1)",
		  "Read random decalibration calibration coefficients");
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
  // Author: Boris Polishchuk (Boris.Polichtchouk at cern.ch)

  TString DBFolder;
  Int_t firstRun   =  0;
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
  
  AliEMCALCalibData *calibda=new AliEMCALCalibData("EMCAL");
  
  Float_t fADCpedestal = 0.0;
  Float_t fADCchannel  = 0.0153;  // 250 GeV / (16*1024)

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
	if (flag == 1) {
	  // Decalibration:
	  // Spread calibration coefficients uniformly around central value 
	  fADCchannel  = rn.Uniform(0.00140,0.00160);
	  fADCpedestal = 0;
	}
	calibda->SetADCchannel (supermodule,column,row,fADCchannel);
	calibda->SetADCpedestal(supermodule,column,row,fADCpedestal);
 	cout<<"Set SM: "<<supermodule<<" col "<<column<<" row "<<row
 	    <<" chan "<<fADCchannel<<" ped fADCpedestal"<<endl;
      }
    }
  }

  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Boris Polichtchouk");
  
  AliCDBId id("EMCAL/Calib/Data",firstRun,lastRun);

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

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
  }
  else if (flag == 1) {
    DBFolder  ="local://DeCalibDB";
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

  TH2F *hPed[nSMod], *hGain[nSMod];
  TCanvas *cPed   = new TCanvas("cPed" ,"Pedestals Mod 1-6"   , 10,10,400,800);
  TCanvas *cGain  = new TCanvas("cGain","Gain factors Mod 1-6", 410,10,400,800);
  TCanvas *cPed2  = new TCanvas("cPed2","Pedestals SMod 7-12", 10,10,400,800);
  TCanvas *cGain2 = new TCanvas("cGain2","Gain factors SMod 7-12", 410,10,400,800);
  cPed   ->Divide(2,3);
  cGain  ->Divide(2,3);
  cPed2  ->Divide(2,3);
  cGain2 ->Divide(2,3);
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
	//  <<" chan "<<gain<<endl;
	hPed[supermodule] ->SetBinContent(column,row,ped);
	hGain[supermodule]->SetBinContent(column,row,gain);
      }
    }
    if(supermodule < 7){
      cPed ->cd(supermodule);
      hPed[supermodule] ->Draw("lego2");
      cGain->cd(supermodule);
      hGain[supermodule]->Draw("lego2");
    }
    else{
      cPed2 ->cd(supermodule-6);
      hPed[supermodule] ->Draw("lego2");
      cGain2->cd(supermodule-6);
      hGain[supermodule]->Draw("lego2");
    }

  }
  cPed   ->Print("pedestals_SM_1_6.eps");
  cGain  ->Print("gains_SM_1-6.eps");
  cPed2  ->Print("pedestals_SM_7-12.eps");
  cGain2 ->Print("gains_SM_7-12.eps");
}
