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
#include "AliT0CalibData.h"
//#include "AliT0AlignData.h"
#include "AliT0Align.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliT0SetCDB()
{
  TControlBar *menu = new TControlBar("vertical","T0 CDB");
  menu->AddButton("Set time delay","SetTimeDelay()",
		  "Set time delay");
  menu->AddButton("Set walk","SetWalk()",
		  "Set slewing coorection");
  menu->AddButton("Set Align","SetAC()",
		  "Set alignment coefficients");
  menu->AddButton("Set LookUpTable","SetLookUp()",
                  "Set LookUp table");
  menu->AddButton("Read time delay","GetTimeDelay()",
		  "Read time delay");
  menu->AddButton("Read walk","GetWalk()",
		  "Read amplitude-time correction");
  menu->AddButton("Read alignment CC","GetAC()",
		  "Read face detector position ");
  menu->AddButton("Read Lookup","GetLookUp()",
		  "Read Lookup table ");
  menu->Show();
}


//------------------------------------------------------------------------
void SetAC()
{
  // Writing alignment coefficients into the Condition DB
  // Arguments:
  
  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 99999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Align";
  firstRun  =  0;
  lastRun   =  99999;
  objFormat = "T0 array  positions";

  AliT0Align *al = new AliT0Align(1,835615);
  al->Run();

   /* 
  AliT0AlignData *alignda=new AliT0AlignData("T0");
  alignda-> SetZposition (67.9,373);
  alignda->Print();
  
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
   AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Align/Positions";


  //  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("T0");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
   AliCDBId id(fPath.Data(),firstRun,lastRun);

   storage->Put(alignda, id, &md);
 }
   */
}
//------------------------------------------------------------------------
void SetTimeDelay()
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 99999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Calib";
  firstRun  =  0;
  lastRun   =  999999;
  objFormat = "T0 initial time delay";

  AliT0CalibTimeEq *calibda=new AliT0CalibTimeEq("T0");
  
   Float_t fTimeDelay  = 1000;
 for(Int_t ipmt=0; ipmt<24; ipmt++) {
   calibda->SetTimeEq(ipmt,fTimeDelay+ipmt*100);
 
  }
 
  calibda->Print();
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  //  AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/TimeDelay";


  // AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("T0");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun,lastRun);
    storage->Put(calibda, id, &md);
  }
}

//------------------------------------------------------------------------
void SetWalk()
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 999999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Calib";
  firstRun  =  0;
  lastRun   =  999999;
  objFormat = "T0 initial slewnig correction";

  AliT0CalibWalk *calibda=new AliT0CalibWalk("T0");
  
   
  TRandom rn;
  
  for(Int_t ipmt=0; ipmt<24; ipmt++) {
    calibda->SetWalk(ipmt);
    calibda->SetAmpLEDRec(ipmt);
  }
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  //  AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/Slewing_Walk";


  // AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("T0");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun,lastRun);
    storage->Put(calibda, id, &md);
  }
}

//------------------------------------------------------------------------
void GetTimeDelay()
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  
  TString DBFolder;
  
  DBFolder  ="local://Calib";
  Int_t nRun=gAlice->GetRunNumber();
  
  AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
  AliCDBEntry* entry = stor->Get("T0/Calib/TimeDelay",0);
  
  AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entry->GetObject();
  clb->Print();
  for (Int_t i=0; i<24; i++) {
   cout<<clb->GetTimeEq(i)<<endl;
   //  cout<<" equalizing CFD "<<(clb->GetTimeDelayCFD(i)-clb->GetTimeDelayCFD(0))<<endl;
 
  }
  
}
//------------------------------------------------------------------------
void GetWalk()
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  
  TString DBFolder;
  
  DBFolder  ="local://Calib";
  Int_t nRun=gAlice->GetRunNumber();
  
   AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
   AliCDBEntry* entry = stor->Get("T0/Calib/Slewing_Walk",0);
   
   AliT0CalibWalk *clb = (AliT0CalibWalk*)entry->GetObject();
   Int_t ipmt=0;
   //  cin>>" enter channel number">>ipmt;
   TGraph *gr = clb->GetAmpLEDRec(ipmt); 
   gr->Draw("AP");
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
  entry = stor2->Get("T0/Align/Positions", nRun);
   
  AliT0AlignData *aln = (AliT0AlignData*)entry->GetObject();
  aln->Print();
}
//------------------------------------------------------------------------
void SetLookUp()
{
  // Writing Lookup table into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Calib";
  firstRun  =  0;
  lastRun   =  10;
  objFormat = "T0 Lookup Table";

  AliT0CalibData *calibda=new AliT0CalibData("T0");

//  calibda->ReadAsciiLookup("lookUpTable.txt");
  calibda->ReadAsciiLookup("lookUpTable_tanay.txt");

  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  //  AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());

  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/LookUp_Table";


  // AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("T0");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun,lastRun);
    storage->Put(calibda, id, &md);
  }
}
//------------------------------------------------------------------------
void GetLookUp()
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:

  TString DBFolder;

  //  DBFolder  ="local://Calib";
  //   Int_t nRun=gAlice->GetRunNumber();
  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBStorage *stor = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
  // cout<<" GetLookUp :: "<<stor<<endl;
  AliCDBEntry *entry;
  //entry = stor->Get("T0/Calib/LookUp_Table",2,0,0);
  entry = stor->Get("T0/Calib/LookUp_Table",1);
  //cout<<"entry="<<entry<<endl;
   cout<<" AliT0CalibData ::GetLookUp :: "<<entry<<endl;
  AliT0CalibData *clb = (AliT0CalibData*)entry->GetObject();
  cout<<" AliT0CalibData *clb "<<clb <<endl;
  //cout<<"clb->a="<<clb->GetA()<<endl;
  //  clb->Dump();
  for (Int_t i=0; i<20; i++) 
    clb->PrintLookupNames("all",i);

}
