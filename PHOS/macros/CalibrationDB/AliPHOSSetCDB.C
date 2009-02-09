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
		  "Explains how to use PHOS CDB menus");
  menu->AddButton("Create ideal calibration","SetCC(0)",
		  "Set equal CC");
  menu->AddButton("Create full decalibration","SetCC(1)",
		  "Set random decalibration CC");
  menu->AddButton("Create residual decalibration","SetCC(2)",
		  "Set residual decalibration calibration coefficients");

  menu->AddButton("Read ideal calibration","GetCC(0)",
		  "Read equal calibration coefficients");
  menu->AddButton("Read full decalibration","GetCC(1)",
		  "Read random decalibration calibration coefficients");
  menu->AddButton("Read residual calibration","GetCC(2)",
		  "Read residial calibration coefficients");
  menu->AddButton("Exit","gApplication->Terminate(0)","Quit aliroot session");
  menu->Show();
}

//-----------------------------------------------------------------------
void Help()
{
  TString string="\nSet calibration parameters and write them into ALICE OCDB:";
  string += "\n\tPress button \"Equal CC\" to create equal calibration coefficients;";
  string += "\n\tPress button \"Full decalibration\" to create \n\t random calibration coefficients with 20\% spread \n\t to imitate fully decalibrated detector;";
  string += "\n\tPress button \"Residual decalibration\" to create \n\t random calibration coefficients with +-2\% spread\n\t to imitate decalibrated detector after initial calibration.\n";
  printf("%s",string.Data());
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
    // Ideal calibration with all channels at nominal value 0.005
    DBFolder  ="local://InitCalibDB";
    firstRun  =  0;
    lastRun   =  999999;
    objFormat = "PHOS ideal pedestals and ADC gain factors (5x64x56)";
    cdb = new AliPHOSCalibData();
    cdb->CreateNew();
  }

  else if (flag == 1) {
    // Full decalibration is +-10% of the nominal value
    DBFolder  ="local://FullDecalibDB";
    firstRun  =  0;
    lastRun   =  999999;
    objFormat = "PHOS fully decalibrated calibration coefficients (5x64x56)";
 
    cdb = new AliPHOSCalibData();    
    cdb->RandomEmc(0.045,0.055);
    cdb->RandomCpv(0.0008,0.0016);
  }
  
  else if (flag == 2) {
    // Residual decalibration is +-1% of the nominal value
    DBFolder  ="local://ResidualCalibDB";
    firstRun  =  0;
    lastRun   =  999999;
    objFormat = "PHOS residual calibration coefficients (5x64x56)";
    
    cdb = new AliPHOSCalibData();    
    cdb->RandomEmc(0.00495,0.00505);
    cdb->RandomCpv(0.00115,0.00125);
  }
  
  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Boris Polichtchouk");
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if(gSystem->Getenv("STORAGE")){
    cout << "Setting specific storage" << endl;
    AliCDBManager::Instance()->SetSpecificStorage("PHOS/*",DBFolder.Data());
  }

  cdb->WriteEmc(firstRun,lastRun,&md);
  cdb->WriteCpv(firstRun,lastRun,&md);
  cdb->WriteEmcBadChannelsMap(firstRun,lastRun,&md);

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

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TString DBFolder;
  Int_t runNumber;

  if      (flag == 0) {
    DBFolder  ="local://InitCalibDB";
    runNumber = 0;
  }
  else if (flag == 1) {
    DBFolder  ="local://FullDecalibDB";
    runNumber = 0;
  }
  else if (flag == 2) {
    DBFolder  ="local://ResidualCalibDB";
    runNumber = 0;
  }

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if(gSystem->Getenv("STORAGE")){
    cout << "Setting specific storage" << endl;
    AliCDBManager::Instance()->SetSpecificStorage("PHOS/*",DBFolder.Data());
  }

  AliPHOSCalibData* clb  = new AliPHOSCalibData(runNumber);

  TH2::AddDirectory(kFALSE);

  TH2F* hGain[5];

  TCanvas *cGain = new TCanvas("cGain","PHOS EMC Gain factors",10,10,700,500);
  cGain->Divide(3,2);

  for (Int_t module=1; module<=nMod; module++) {
    TString nameGain="hGain";
    nameGain+=module;
    TString titleGain="Gain factors in module ";
    titleGain+=module;
    hGain[module-1] = new TH2F(nameGain.Data(),titleGain.Data(),
			    nCol,1.,1.*nCol,nRow,1.,1.*nRow);

    for (Int_t column=1; column<=nCol; column++) {
      for (Int_t row=1; row<=nRow; row++) {
        Float_t gain = clb->GetADCchannelEmc (module,column,row);
	hGain[module-1]->SetBinContent(column,row,gain);
      }
    }
    cGain->cd(module);
    hGain[module-1]->SetMinimum(0);
    hGain[module-1]->Draw("colz");
  }
  cGain->Print("gains.eps");
}
