/* $Id: AliT0SetCDB.C 22796 2007-12-06 11:32:28Z alla $ */

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
  menu->AddButton("Set LookUpTable","SetLookUp()",
                  "Set LookUp table");
  menu->AddButton("Read time delay","GetTimeDelay()",
		  "Read time delay");
  menu->AddButton("Read walk","GetWalk()",
		  "Read amplitude-time correction");
  menu->AddButton("Read Lookup","GetLookUp()",
		  "Read Lookup table ");
  menu->Show();
}


//------------------------------------------------------------------------
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
  calibda-> ComputeOnlineParams("t0treeDA08000025762005.10.root");
  /*
   Float_t fTimeDelay  = 1000;
 for(Int_t ipmt=0; ipmt<24; ipmt++) {
   calibda->SetTimeEq(ipmt,fTimeDelay+ipmt*100);
 
  }
  */
  calibda->Print();
  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //AliCDBManager::Instance()->SetDefaultStorage("local:///home/alla/alice/testFeb08/OCDB/");

  //  AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/TimeDelay";
  //TString fPath="TimeDelay";
  cout<<fPath.Data()<<endl;

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
  lastRun   =  9999999;
  objFormat = "T0 initial slewnig correction";

  AliT0CalibWalk *calibda=new AliT0CalibWalk("T0");
  //  cout<<"AliT0CalibWalk "<< calibda<<endl;
  //  calibda->Dump();
  const char *filename="t0tree08000025765005.10.root";
 calibda->MakeWalkCorrGraph(filename);

 
  /*  
  TRandom rn;
  
  for(Int_t ipmt=0; ipmt<24; ipmt++) {
    calibda->SetWalk(ipmt);
    calibda->SetAmpLEDRec(ipmt);
  }
  */


  //Store calibration data into database
 AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 //AliCDBManager::Instance()->SetDefaultStorage("local:///home/alla/alice/testFeb08/OCDB/");

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
       // AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///home/alla/alice/testFeb08/OCDB/");
  AliCDBEntry* entry = stor->Get("T0/Calib/TimeDelay",25068);
  
  AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entry->GetObject();
  //  clb->Print();
  for (Int_t i=0; i<24; i++) {
    cout<<i<<"  "<<clb->GetTimeEq(i)<<endl;
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
  //  AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///home/alla/alice/testFeb08/OCDB/");
 
       AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
      AliCDBEntry* entry = stor->Get("T0/Calib/Slewing_Walk",25081);

      //  AliT0Parameters* param = AliT0Parameters::Instance();
      //param->Init();
 
   
   AliT0CalibWalk *clb = (AliT0CalibWalk*)entry->GetObject();
   //  cin>>" enter channel number">>ipmt;
   
   TCanvas *c1 = new TCanvas("c1", "CFD C side",0,48,1280,951);
   c1->Divide(4,3);
   for (Int_t i=0; i<12; i++) {
     c1->cd(i+1);
     // TGraph* fu = param ->GetWalk(ipmt);
     //     TGraph *gr = clb->GetWalk(i);
       TGraph *gr = clb->GetAmpLEDRec(i); 
     //     cout<<"   "<<gr<<endl; 
     if(gr) {
       gr->SetMarkerStyle(7);
       gr->Draw("AP");
     }
   }
   TCanvas *c2 = new TCanvas("c2", "CFD A side",0,48,1280,951);
   c2->Divide(4,3);
   for (Int_t i=12; i<24; i++) {
     c2->cd(i+1-12);
     //     if(i==15) continue;
     // TGraph *gr = clb->GetWalk(i); 
     //  cout<<i<<" "<<gr->GetN()<<" "<<endl;
      TGraph *gr = clb->GetAmpLEDRec(i); 
      //  TGraph* fu = param ->GetWalk(pmt+12);
      cout<<i<<" "<<gr->GetN()<<" "<<endl;
     
     gr->SetMarkerStyle(7);
     //     gr->SetMarkerSize(20);
     gr->Draw("AP");
   }
    
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
