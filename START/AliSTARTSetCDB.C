/* $Id$ */

// Script to create calibration parameters and store them into CDB
// Two sets of calibration parameters can be created:
// 1) equal parameters
// 2) randomly distributed parameters for decalibrated detector silumations

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "AliRun.h"
#include "AliSTARTCalibData.h"
#include "AliSTARTAlignData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliSTARTSetCDB()
{
  TControlBar *menu = new TControlBar("vertical","START CDB");
  menu->AddButton("Set Calib","SetCC()",
		  "Set calibration coefficients");
  menu->AddButton("Set Align","SetAC()",
		  "Set alignment coefficients");
  menu->AddButton("Read calibration CC","GetCC()",
		  "Read calibration  coefficients");
  menu->AddButton("Read alignment CC","GetAC()",
		  "Read face detector position ");
  menu->Show();
}


//------------------------------------------------------------------------
void SetAC()
{
  // Writing alignment coefficients into the Condition DB
  // Arguments:
  
  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Align";
  firstRun  =  0;
  lastRun   =  10;
  objFormat = "START array Z positions";

  
  AliSTARTAlignData *alignda=new AliSTARTAlignData("START");
  alignda-> SetZposition (67.9,373);
  alignda->Print();
  
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
   AliCDBManager::Instance()->SetSpecificStorage("START",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="START/Align/Positions";


  //  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("START");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
 if(storage) {
   AliCDBId id(fPath.Data(),firstRun,lastRun);

   storage->Put(alignda, id, &md);
 }
}
//------------------------------------------------------------------------
void SetCC()
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Calib";
  firstRun  =  0;
  lastRun   =  10;
  objFormat = "START initial gain factors, time delay, slewnig";

  AliSTARTCalibData *calibda=new AliSTARTCalibData("START");
  
  Float_t fGain = 1;
  Float_t fTimeDelay  = 200;
  
  TRandom rn;
  
  for(Int_t ipmt=0; ipmt<24; ipmt++) {
    calibda->SetGain (fGain,ipmt);
    calibda->SetTimeDelayCFD(fTimeDelay,ipmt);
    calibda->SetTimeDelayLED(fTimeDelay,ipmt);
    calibda->SetWalk(ipmt,"data/re.root");
    calibda->SetSlewingLED(ipmt,"data/CFD-LED.txt");
    calibda->SetSlewingRec(ipmt,"data/CFD-LED.txt");
    Double_t value=calibda->GetSlewingLED(ipmt,300);
    Double_t rec= calibda->GetSlewingRec(ipmt, value);
    cout<<" in "<<value<<" out "<<rec<<endl;
  }
  calibda->Print();
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");

  //  AliCDBManager::Instance()->SetSpecificStorage("START",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="START/Calib/Gain_TimeDelay_Slewing_Walk";


  // AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("START");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun,lastRun);
    storage->Put(calibda, id, &md);
  }
}

//------------------------------------------------------------------------
void GetCC()
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  
  TString DBFolder;
  
  DBFolder  ="local://Calib";
  Int_t nRun=gAlice->GetRunNumber();
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *stor2 = man->GetStorage("local://Calib");
  AliCDBEntry *entry;
  entry = stor2->Get("START/Calib/Gain_TimeDelay_Slewing_Walk");
   
  AliSTARTCalibData *clb = (AliSTARTCalibData*)entry->GetObject();
  clb->Print();
}
//------------------------------------------------------------------------
void GetAC()
{
  // Read align coefficients into the Calibration DB
  // Arguments:
  
  TString DBFolder;
  
  DBFolder  ="local://Align";
  Int_t nRun=gAlice->GetRunNumber();
  
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBStorage *stor2 = man->GetStorage("local://Align");
  AliCDBEntry *entry;
  entry = stor2->Get("START/Align/Positions", nRun);
   
  AliSTARTAlignData *aln = (AliSTARTAlignData*)entry->GetObject();
  aln->Print();
}
