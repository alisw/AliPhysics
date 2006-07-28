/* $Id$ */

// Script to create calibration parameters and store them into CDB
// Two sets of calibration parameters can be created:
// 1) equal parameters
// 2) randomly distributed parameters for decalibrated detector silumations

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliRun.h"
#include "AliPHOSCalibData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

static const Int_t nMod =  5;
static const Int_t nCol = 56;
static const Int_t nRow = 64;

void AliPHOSSetCDB()
{
  TControlBar *menu = new TControlBar("vertical","PHOS CDB");
  menu->AddButton("Help to run PHOS CDB","Help()",
		  "Explains how to use PHOS CDS menus");
  menu->AddButton("Equal CC","SetCC(0)",
		  "Set equal calibration coefficients");
  menu->AddButton("Decalibrate","SetCC(1)",
		  "Set random decalibration calibration coefficients");
  menu->AddButton("Set calibration equal to invers decalibration coefficients","SetCC(2)",
		  "Set calibration coefficients inverse to decalibration ones");
  menu->AddButton("Residual calibration","SetCC(3)",
		  "Set residual decalibration calibration coefficients");

  menu->AddButton("Read equal CC","GetCC(0)",
		  "Read initial equal calibration coefficients");
  menu->AddButton("Read random CC","GetCC(1)",
		  "Read random decalibration calibration coefficients");
  menu->AddButton("Read inverse CC","GetCC(2)",
		  "Read calibration coefficients inverse to decalibration ones");
  menu->AddButton("Read residual CC","GetCC(3)",
		  "Read residial calibration coefficients");
  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\nSet calibration parameters and write them into ALICE CDB.
Press button \"Equal CC\" to create equal pedestals and gain factors.
Press button \"Decalibrate\" to create random pedestals and gain factors to imitate decalibrated detector\n";
  printf(string);
}

//------------------------------------------------------------------------
void SetCC(Int_t flag=0)
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: decalibration coefficients
  //   flag=2: calibration coefficients equal to inverse decalibration ones
  // Author: Boris Polishchuk (Boris.Polichtchouk at cern.ch)

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    =  0;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  AliPHOSCalibData* cdb = 0;

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
    firstRun  =  0;
    lastRun   =  0;
    objFormat = "PHOS initial pedestals and ADC gain factors (5x64x56)";
    cdb = new AliPHOSCalibData();
    cdb->CreateNew();
  }

  else if (flag == 1) {
    DBFolder  ="local://DeCalibDB";
    firstRun  =  100;
    lastRun   =  100;
    objFormat = "PHOS decalibration pedestals and ADC gain factors (5x64x56)";
 
    cdb = new AliPHOSCalibData();    
    cdb->RandomEmc();
    cdb->RandomCpv();
  }
  
  else if (flag == 2) {
    // First read decalibration DB
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    AliCDBManager::Instance()->SetSpecificStorage("PHOS","local://DeCalibDB");
    AliPHOSCalibData* cdbDecalib = new AliPHOSCalibData(100);

    DBFolder  ="local://InverseCalibDB";
    firstRun  = 200;
    lastRun   = 200;
    objFormat = "PHOS calibration parameters equal to inverse decalibration ones (5x64x56)";
    cdb = new AliPHOSCalibData();    

    // Read EMC decalibration parameters and put inverse values to a new artificial CDB

    for (Int_t module=1; module<=nMod; module++) {
      for (Int_t column=1; column<=nCol; column++) {
	for (Int_t row=1; row<=nRow; row++) {
	  Float_t valueGain = cdbDecalib->GetADCchannelEmc (module,column,row);
	  Float_t valuePed  = cdbDecalib->GetADCpedestalEmc(module,column,row);
	  cdb->SetADCchannelEmc(module,column,row,1./valueGain);
	  cdb->SetADCpedestalEmc(module,column,row,valuePed);
	}
      }
    }

    // Read CPV decalibration parameters and put inverse values to a new artificial CDB

    for (Int_t module=1; module<=nMod; module++) {
      for (Int_t column=1; column<=nCol*2; column++) {
	for (Int_t row=1; row<=nRow; row++) {
	  Float_t valueGain = cdbDecalib->GetADCchannelCpv (module,column,row);
	  Float_t valuePed  = cdbDecalib->GetADCpedestalCpv(module,column,row);
	  cdb->SetADCchannelCpv(module,column,row,1./valueGain);
	  cdb->SetADCpedestalCpv(module,column,row,valuePed);
	}
      }
    }
  }
  else if (flag == 3) {
    DBFolder  ="local://ResidualCalibDB";
    firstRun  =  0;
    lastRun   =  999999;
    objFormat = "PHOS residual ADC gain factors (5x64x56)";
    
    cdb = new AliPHOSCalibData();    
    cdb->RandomEmc(0.98,1.02);
    cdb->RandomCpv(0.00115,0.00125);
  }
  
  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Boris Polichtchouk");
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS",DBFolder.Data());

  cdb->WriteEmc(firstRun,lastRun,&md);
  cdb->WriteCpv(firstRun,lastRun,&md);

}

//------------------------------------------------------------------------
void GetCC(Int_t flag=0)
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: decalibration coefficients
  //   flag=2: calibration coefficients equal to inverse decalibration ones
  // Author: Yuri.Kharlov at cern.ch

  TString DBFolder;
  Int_t runNumber;

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
    runNumber = 0;
  }
  else if (flag == 1) {
    DBFolder  ="local://DeCalibDB";
    runNumber = 100;
  }
  else if (flag == 2) {
    DBFolder  ="local://InverseCalibDB";
    runNumber = 200;
  }
  else if (flag == 2) {
    DBFolder  ="local://ResidualCalibDB";
    runNumber = 0;
  }

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS",DBFolder.Data());

  AliPHOSCalibData* clb  = new AliPHOSCalibData(runNumber);

  TH2::AddDirectory(kFALSE);

  TH2F* hPed[5];
  TH2F* hGain[5];

  TCanvas *cPed  = new TCanvas("cPed" ,"PHOS EMC Pedestals"   , 10,10,400,800);
  TCanvas *cGain = new TCanvas("cGain","PHOS EMC Gain factors",410,10,400,800);
  cPed ->Divide(1,5);
  cGain->Divide(1,5);

  for (Int_t module=1; module<=nMod; module++) {
    TString namePed="hPed";
    namePed+=module;
    TString titlePed="Pedestals in module ";
    titlePed+=module;
    hPed[module-1] = new TH2F(namePed.Data(),titlePed.Data(),
			    nCol,1.,1.*nCol,nRow,1.,1.*nRow);

    TString nameGain="hGain";
    nameGain+=module;
    TString titleGain="Gain factors in module ";
    titleGain+=module;
    hGain[module-1] = new TH2F(nameGain.Data(),titleGain.Data(),
			    nCol,1.,1.*nCol,nRow,1.,1.*nRow);

    for (Int_t column=1; column<=nCol; column++) {
      for (Int_t row=1; row<=nRow; row++) {
	Float_t ped  = clb->GetADCpedestalEmc(module,column,row);
	Float_t gain = clb->GetADCchannelEmc (module,column,row);
	hPed[module-1]->SetBinContent(column,row,ped);
	hGain[module-1]->SetBinContent(column,row,gain);
      }
    }
    cPed ->cd(module);
    hPed[module-1]->Draw("lego2");
    cGain->cd(module);
    hGain[module-1]->Draw("lego2");
  }
  cPed ->Print("pedestals.eps");
  cGain->Print("gains.eps");
}
