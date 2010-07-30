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


void AliT0SetCDBcosmic()
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
  menu->AddButton("SetLatency","setLat()",
		  "SetLatency ");
  menu->AddButton("ReadLatency","readLat()",
		  "print Latency ");
  menu->Show();
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
void SetTimeDelay()
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:

  firstRun  =  125001;
  lastRun   =  125800;
  Int_t beamPeriod =  1;
  char*   objFormat = "T0 initial time delay";

  DBFolder  ="local://Calib";

  //    Int_t shift[24] = {0, 0,0,0,0,0,0,0,0,0,0,0, 
  //	     0, 0,0,0,0,0,0,0,0,0,0,0};
  //   Int_t shift[24] = {0, };
  //run 114786
  //  Int_t shift[24] = {-17, 0, 14, 48, 36, 35, 16 , -23, , 72, 50, 44,-15, 
  //  		     0, 13, -32, -19, 32, -2, 20, 43, 20, 78, 43, 118};    
  //run 115318
  //  Int_t shift[24] = {-17, 0, 14, 50, 36, 34, 14, -22, 72, 50, 42, -15,
  //  		     0, 13, -32, -19, 32, -2, 20, 43, 20, 78, 43, 118};    
  //run 116562
  //    Int_t shift[24] = {-17, 0, 14, 50, 36, 34, 14, -22, 72, 50, 42, -15,
  //   		     0, 13, -32, -19, 32, -2, 20, 43, 40, 78, 43, 118};    
  // run 117118
  //  Int_t shift[24] = {0, 0, 14, 50, 37, 35, 16, -22, 72, 52, 42, -14,
  //    		     0, 14, -32, -19, 31, -1, 20, 43, 37, 79, 42, 119};    
  // run 117112
  //  Int_t shift[24] = {0, 0, 14, 5, 37, 34, 15, -22, 72, 50, 42, -15,
  // 		     0, 14, -32, -19, 31, -1, 20, 43, 37, 79, 42, 119};    
   // run 118000
  //  Int_t shift[24] = { 17,    0,  14,  47,  34 , 33, 18,  -21,   72,  51,  42,  -15, 
  //		       -15,  0,   -47,  -36,  22,   -15,  9, 26,   1, 60, 23, 106};    
  // run119163
  //  Int_t shift[24] = {18, 0, 16, 49, 36, 34, 22, -22, 71, 51, 43, -14,
  //	      -15, 0, -48, -37, 19, -17, 9, 30, 4, 62,  24, 103};

  // run 120076
  //    Int_t shift[24] = {19,0,16,51,36,37,24,-15,74,54,43,-13,
  //     -18,0,-42,-34,16,-16,8,27,4,57,24,104};    

  // run 120244
  //  Int_t shift[24] = {14,0,13,48,34,35,22,-20,73,52,40,-14,
  // 		     0, 16,-27,-16, 32,0,20,45,19,72,40,121};    
  // run 120824
  // Int_t shift[24] = {34,0,11,51,37,35,21,-20,70,47,40,-12,
  //		     0, 16,-27,-16, 32,0,20,45,19,72,40,119};    
  //124187

  //  Int_t shift[24] = {0, 0, 15, 50, 36, 35, 23, -20, 83, 53, 44, -11,
  //		     0 , 20 , -26 , -8 , 36 , 4 , 21 , 47 , 20 , 78 , 52 , 123};
  //124702
  //  Int_t shift[24] = {20, 0 ,15, 50, 37, 36, 23, -20, 77, 54, 44, -12, 
  //  		     0, 19, -25, -13, 35, 1, 23, 46, 20, 77, 50, 122};

  //125085
  //  Int_t shift[24] = {0, 0, 15, 51, 36, 34, 26, -20, 81, 53, 44, -12, 
  //	       0 , 20 , -26 , -7 , 36 , 4 , 22 , 48 , 18 , 79 , 50 , 121};
//

//  125097
  Float_t  shift[24]={8,0, 16, 50, 36, 34, 25, -20, 78, 54, 42, -11, 
		      0, 21, -21, -7, 38, 5, 27, 50, 20, 78, 53, 124};

  //125295
  //  Float_t  shift[24]={16, 0,  2, 2, 30, 3, -5, -8, 5, -4, 16, -5, 
  //			 0, 10, 10, 12, 13, -24, 15, 26, -2, 17, 10, -30};
  //125842
  //   Float_t shift[24]={16, 0, 1, 3,  32,   4, -2, -8, 5, -3, 15, -5, 
  //		      0, 10, 9, 12, 14, -23, 16, 26, 0, 19, 9, -31};



   //126407
  //  Float_t shift[24]={16, 0,   1, 4, 33, 4, 0, -7, 6, -3, 15, -4,
  //		     0 , 10, 10, 13, 13, -23, 16, 28, 0, 19, 11, -29};
//
   //for ( Int_t indexfile=filestart; indexfile < filestop+1;indexfile++ ) 

  AliT0CalibTimeEq *calibda=new AliT0CalibTimeEq("T0");
  //  calibda-> ComputeOnlineParams("t0treeDA08000025762005.10.root");
  
   for(Int_t ipmt=0; ipmt<24; ipmt++) {
   calibda->SetTimeEq(ipmt,shift[ipmt]);
 
  }
  
  calibda->Print();
  //Store calibration data into database
  // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetDefaultStorage("local:///scratch/alla/alice/Jun10/TestCDB/");
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
    //  AliCDBId id(fPath.Data(),firstRun,AliCDBRunRange::Infinity() );
     AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(calibda, id, &md);
  }
}

//------------------------------------------------------------------------
void SetWalk()
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   = 1000;
  Int_t lastRun    = 999999999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://Calib";
  objFormat = "T0 initial slewnig correction";

  AliT0CalibWalk *calibda=new AliT0CalibWalk("T0");
  //  cout<<"AliT0CalibWalk "<< calibda<<endl;
  //  calibda->Dump();
  const char *filename="amphist616.root";
 calibda->MakeWalkCorrGraph(filename);

 
  //Store calibTestCDB/T0/Calib/Slewing_Walk/ration data into database
 // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 AliCDBManager::Instance()->SetDefaultStorage("local:///home/alla/alice/Mar10/TestCDB/");

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
  // Int_t nRun=gAlice->GetRunNumber();
  
  //     AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
     //    AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://");
  AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///scratch/alla/alice/Jun10/TestCDB/");
  AliCDBEntry* entry = stor->Get("T0/Calib/TimeDelay",125095);
  
  AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entry->GetObject();
  //  clb->Print();
  for (Int_t i=0; i<24; i++) {
    cout<<"  "<<clb->GetTimeEq(i)<<" ";
   //  cout<<" equalizing CFD "<<(clb->GetTimeDelayCFD(i)-clb->GetTimeDelayCFD(0))<<endl;
  }
  
    cout<<endl;
}
//------------------------------------------------------------------------
void GetWalk()
{
  // Read calibration coefficients into the Calibration DB
  // Arguments:
  
  TString DBFolder;
  
  DBFolder  ="local://Calib";
  Int_t runNumber=127001;
  // Int_t nRun=gAlice->GetRunNumber();
  AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///scratch/alla/alice/Jun10/TestCDB/");
  // AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///home/alla/alice/May10/TestCDB/");
 
  //  AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
   AliCDBEntry* entry = stor->Get("T0/Calib/Slewing_Walk",runNumber);

      //  AliT0Parameters* param = AliT0Parameters::Instance();
      //param->Init();
 
   
   AliT0CalibWalk *clb = (AliT0CalibWalk*)entry->GetObject();
   //  cin>>" enter channel number">>ipmt;
   TString buf3;
   TCanvas *c1 = new TCanvas("c1", "LED-CFD C side",0,48,1280,951);
   c1->Divide(4,3);
   for (Int_t i=0; i<12; i++) {
     c1->cd(i+1);
     TGraph *gr = clb->GetAmpLED(i); 
     if(gr) {
       gr->GetXaxis()->SetTitle("led-cfd");
       gr->GetYaxis()->SetTitle("MIPs");
       gr->SetMarkerStyle(20);
       gr->Draw("AP");
     }
   }
   buf3 = Form("ampplots/ampLEDminCFD_C_%i.png",runNumber);
  c1->Print(buf3);

   
   TCanvas *c2 = new TCanvas("c2", "LED-CFD  A side",0,48,1280,951);
   c2->Divide(4,3);
   for (Int_t i=12; i<24; i++) {
     c2->cd(i+1-12);
     TGraph *gr = clb->GetAmpLED(i); 
     gr->GetXaxis()->SetTitle("led-cfd");
     gr->GetYaxis()->SetTitle("MIPs");
     gr->SetMarkerStyle(20);
     gr->Draw("AP");
   }
  buf3 = Form("ampplots/ampLEDminCFD_A_%i.png",runNumber);

   c2->Print(buf3);
   
   TCanvas *c3 = new TCanvas("c3", "QTC C side",0,48,1280,951);
   c3->Divide(4,3);
   for (Int_t i=0; i<12; i++) {
     c3->cd(i+1);
     TGraph *gr = clb->GetQTC(i);
     if(gr) {
      gr->SetTitle(Form("PMT%i",i));
      gr->GetXaxis()->SetTitle("qtc");
       gr->GetYaxis()->SetTitle("MIPs");
       gr->SetMarkerStyle(20);
       gr->Draw("AP");
     }
   }
   buf3 = Form("ampplots/ampQTC_C_%i.png",runNumber);

   c3->Print(buf3);
   
   TCanvas *c4 = new TCanvas("c4", "QTC  A side",0,48,1280,951);
   c4->Divide(4,3);
   for (Int_t i=12; i<24; i++) {
     c4->cd(i+1-12);
     TGraph *gr = clb->GetQTC(i); 
     //   TGraph *gr = clb->GetWalk(i); 
     
     //  TGraph *gr = clb->GetAmpLEDRec(i);              
     gr->SetTitle(Form("PMT%i",i));
     gr->GetXaxis()->SetTitle("qtc");
     gr->GetYaxis()->SetTitle("MIPs");
     gr->SetMarkerStyle(20);
     //     gr->SetMarkerSize(20);
     gr->Draw("AP");
   }
   buf3 = Form("ampplots/ampQTC_A_%i.png",runNumber);

   c4->Print(buf3);
   
   TCanvas *c5 = new TCanvas("c5", "walk LED-CFD C side",0,48,1280,951);
   c5->Divide(4,3);
   for (Int_t i=0; i<12; i++) {
     c5->cd(i+1);
     TGraph *gr = clb->GetAmpLEDRec(i); 
          
     if(gr) {
     gr->SetTitle(Form("PMT%i",i));
       gr->GetXaxis()->SetTitle("led-cfd");
       gr->GetYaxis()->SetTitle("walk");
       gr->SetMarkerStyle(20);
       gr->Draw("AP");
       if(i==0) gr->Print();
     }
   }

   buf3 = Form("ampplots/walkLEDminCFD_C_%i.png",runNumber);
   c5->Print(buf3);

   TCanvas *c6 = new TCanvas("c6", "walk LED-CFD  A side",0,48,1280,951);
   c6->Divide(4,3);
   for (Int_t i=12; i<24; i++) {
     c6->cd(i+1-12);
     TGraph *gr = clb->GetAmpLEDRec(i);              
     gr->SetTitle(Form("PMT%i",i));
     gr->GetXaxis()->SetTitle("led-cfd");
     gr->GetYaxis()->SetTitle("walk");
     gr->SetMarkerStyle(20);
     gr->Draw("AP");
     
   }

   buf3 = Form("ampplots/walkLEDminCFD_A_%i.png",runNumber);
   c6->Print(buf3);
   
   TCanvas *c7 = new TCanvas("c7", "walk QTC C side",0,48,1280,951);
   c7->Divide(4,3);
   for (Int_t i=0; i<12; i++) {
     c7->cd(i+1);
     TGraph *gr = clb->GetWalk(i); 
          
     if(gr) {
     gr->SetTitle(Form("PMT%i",i));
       gr->GetXaxis()->SetTitle("qtc");
       gr->GetYaxis()->SetTitle("walk");
       gr->SetMarkerStyle(20);
       gr->Draw("AP");
       if(i==0) gr->Print();
     }
   }
   buf3 = Form("ampplots/walkQTC_C_%i.png",runNumber);
   c7->Print(buf3);

   TCanvas *c8 = new TCanvas("c8", "walk QTC  A side",0,48,1280,951);
   c8->Divide(4,3);
   for (Int_t i=12; i<24; i++) {
     c8->cd(i+1-12);
    TGraph *gr = clb->GetWalk(i);              
     gr->SetTitle(Form("PMT%i",i));
      gr->GetXaxis()->SetTitle("qtc");
     gr->GetYaxis()->SetTitle("walk");
     gr->SetMarkerStyle(20);
     gr->Draw("AP");
   }
   buf3 = Form("ampplots/walkQTC_A_%i.png",runNumber);
   c8->Print(buf3);
 


   
}
//------------------------------------------------------------------------
void SetLookUp()
{
  // Writing Lookup table into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 9999999;
  Int_t beamPeriod =  1;
  char* objFormat = "T0 Lookup Table";

  AliT0CalibData *calibda=new AliT0CalibData("T0");

//  calibda->ReadAsciiLookup("lookUpTable.txt");
  calibda->ReadAsciiLookup("lookUpTable_tanay.txt");

  //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local:///home/alla/alice/testOct09/TestCDB");
    
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
//--------------------------------------------------------
void setLat()
{
  // Arguments:
  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Alla");
  metaData.SetComment("Latency");
  //Store calibration data into database
  
  AliT0CalibLatency *calibda=new AliT0CalibLatency("T0");
  
  calibda->SetLatencyHPTDC(9000);

  //124702
  //  calibda-> SetLatencyL1(8.91358e+03);
  //  calibda-> SetLatencyL1A(8.91352e+03);
  //  calibda-> SetLatencyL1C(8.91361e+03);

 //125097 
  calibda-> SetLatencyL1 (8.91406e+03)  ;
  calibda-> SetLatencyL1A( 8.91401e+03);
  calibda-> SetLatencyL1C (8.91412e+03) ;

 //run 125295
 ///   calibda-> SetLatencyL1(8.914520e+03) ;
 //   calibda-> SetLatencyL1A( 8.914860e+03) ;
 //   calibda-> SetLatencyL1C(8.914180e+03);
  //125842
  //   calibda-> SetLatencyL1(8.91306e+03);
  //   calibda-> SetLatencyL1A (8.91338e+03);
  //  calibda-> SetLatencyL1C (8.91274e+03);
     //126407
  // calibda-> SetLatencyL1 (8.91345e+03);
  ///  calibda->SetLatencyL1A (8.91378e+03);
  //  calibda->SetLatencyL1C (8.91311e+03);

 
  Int_t beamPeriod =  1;
  char*   objFormat = "T0 initial time delay";
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");

 TString fPath="T0/Calib/Latency";
  
  AliCDBManager::Instance()->SetDefaultStorage("local:///scratch/alla/alice/Jun10/TestCDB");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    //    AliCDBId id(fPath.Data(), 126401 ,AliCDBRunRange::Infinity());
     AliCDBId id(fPath.Data(), 125001 , 125200);
    storage->Put(calibda, id, &metaData);
  calibda->Print();

  }
}
//--------------------------------------------------------
void readLat()
{
  // Arguments:
  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Alla");
  metaData.SetComment("Latency");
  //Store calibration data into database
    AliCDBStorage *stor =AliCDBManager::Instance()->GetStorage("local:///scratch/alla/alice/Jun10/TestCDB/");
  AliCDBEntry* entry = stor->Get("T0/Calib/Latency",124401);

  AliT0CalibLatency *calibda=(AliT0CalibLatency*)entry->GetObject();
  calibda->Print();

}
