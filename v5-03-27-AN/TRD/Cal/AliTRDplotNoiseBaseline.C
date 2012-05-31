#if !defined( __CINT__) || defined(__MAKECINT__)


#include <Riostream.h>
#include <TSystem.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGrid.h>


#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"


#include "../TRD/AliTRDarrayF.h"
#include "../TRD/AliTRDCalibPadStatus.h"
#include "../TRD/Cal/AliTRDCalPadStatus.h"
#include "../TRD/Cal/AliTRDCalDet.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalROC.h"
#include "../TRD/AliTRDcalibDB.h"


#endif


//void PlotNoiseBaseline(Int_t run, Int_t sm, Int_t det, const char * pathdatabase="local:///d/alice12/bailhache/TestShuttle/database/", const char * pathreferencefile="local:///d/alice12/bailhache/TestShuttle/reference")
//void PlotNoiseBaseline(Int_t run=34529, Int_t sm=0, Int_t det=0, const char * pathdatabase="alien://Folder=/alice/data/2008/LHC08b/OCDB/", const char * pathreferencedatabase="alien://Folder=/alice/data/2008/LHC08b/Reference/")
//void PlotNoiseBaseline(Int_t run=1, Int_t sm=0, Int_t det=0, const char * pathdatabase="local:///d/alice12/bailhache/AliAnalysisTask/v4-13-Head/SHUTTLE/TestShuttle/TestCDB/", const char * pathreferencedatabase="local:///d/alice12/bailhache/AliAnalysisTask/v4-13-Head/SHUTTLE/TestShuttle/TestReference/")
void AliTRDplotNoiseBaseline(Int_t run=34529, Int_t sm=0, Int_t det=0, const char * pathdatabase="alien://Folder=/alice/data/2008/LHC08b/OCDB/", const char * pathreferencedatabase="alien://Folder=/alice/data/2008/LHC08b/Reference/")
{

  //TGrid::Connect("alien://",0,0,"t");

  AliCDBManager *CDB = AliCDBManager::Instance();
  CDB->SetDefaultStorage("local://$ALICE_ROOT/OCDB"); 
  CDB->SetSpecificStorage("TRD/Calib/PadNoise",pathdatabase);
  CDB->SetSpecificStorage("TRD/Calib/DetNoise",pathdatabase);
  CDB->SetSpecificStorage("TRD/Calib/PadStatus",pathdatabase);
  CDB->SetRun(run);

  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();

  //const AliTRDCalDet *u = cal->GetNoiseDet();

  AliTRDCalDet *u = new AliTRDCalDet("u","u");
  for(Int_t k = 0; k < 540; k++){
    u->SetValue(k,10.0);
  }

  //Style
  //************************
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  //Build the Cal Pad
  //********************************
  Int_t smi = sm*30;
  AliTRDCalPad *ki = new AliTRDCalPad("testnoise","testnoise");
  for(Int_t k = 0; k < 540; k++){
    ki->SetCalROC(k,(AliTRDCalROC *) cal->GetNoiseROC(k));   
  }

  // padstatus 2D
  Int_t smn = (Int_t) det/30;
  if((smn==0) || (smn==1) || (smn==2) || (smn==9) || (smn==10) || (smn==11)) smn = 1;
  if((smn==3) || (smn==4) || (smn==5) || (smn==12) || (smn==13) || (smn==14)) smn = 2;
  if((smn==6) || (smn==7) || (smn==8) || (smn==15) || (smn==16) || (smn==17)) smn = 3;
  TString name("TRD/DAQData/PadStatus");
  name += smn;
  //name += 3;
  AliCDBEntry *entrypadstatus = AliCDBManager::Instance()->Get("TRD/Calib/PadStatus",run);
  if(!entrypadstatus) return;
  AliTRDCalPadStatus *lo = (AliTRDCalPadStatus *)entrypadstatus->GetObject();
  AliCDBEntry *entryo = AliCDBManager::Instance()->GetStorage(pathreferencedatabase)->Get(name, run);
  if(!entryo) return;
  AliTRDCalibPadStatus *calpad = (AliTRDCalibPadStatus *) entryo->GetObject();
  if(!calpad) return;
 

  // Plot
  //***********
 
  
  // noise 2D
  TCanvas *cnoise = new TCanvas((const char*)"noise1",(const char*)"noise1",50,50,600,800);
  cnoise->Divide(3,2);
  cnoise->cd(1);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,0,u,0,0.0,3.5,-1))->Draw("colz");
  cnoise->cd(2);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,1,u,0,0.0,3.5,-1))->Draw("colz");
  cnoise->cd(3);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,2,u,0,0.0,3.5,-1))->Draw("colz");
  cnoise->cd(4);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,3,u,0,0.0,3.5,-1))->Draw("colz");
  cnoise->cd(5);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,4,u,0,0.0,3.5,-1))->Draw("colz");
  cnoise->cd(6);
  ((TH2F *)ki->MakeHisto2DSmPl(sm,5,u,0,0.0,3.5,-1))->Draw("colz");
  

  // Pad Status 
  TCanvas *cpadstatus = new TCanvas((const char*)"padstatus",(const char*)"padstatus",50,50,600,800);
  cpadstatus->Divide(3,2);
  cpadstatus->cd(1);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,0))->Draw("colz");
  cpadstatus->cd(2);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,1))->Draw("colz");
  cpadstatus->cd(3);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,2))->Draw("colz");
  cpadstatus->cd(4);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,3))->Draw("colz");
  cpadstatus->cd(5);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,4))->Draw("colz");
  cpadstatus->cd(6);
  ((TH2F *)lo->MakeHisto2DSmPl(sm,5))->Draw("colz");

    
 
  // reference data 

  TCanvas *cpoui = new TCanvas((const char*)"cpoui",(const char*)"cpoui",50,50,600,800);
  cpoui->cd();
  ((TH2F *)calpad->GetHisto(det))->Draw("lego");


  AliTRDCalROC *ouip = calpad->GetCalRocMean(det);
  TCanvas *cpouilo = new TCanvas((const char*)"cpouilo",(const char*)"cpouilo",50,50,600,800);
  cpouilo->Divide(2,1);
  cpouilo->cd(1);
  ((TH1F *)ouip->MakeHisto1D(8.5,10.5,-1,10.0))->Draw();
  //((TH1F *)ouip->MakeHisto1D(0.85,1.05,-1))->Draw();
  cpouilo->cd(2);
  ((TH2F *)ouip->MakeHisto2D(8.5,10.5,-1,10.0))->Draw("colz");
  //((TH2F *)ouip->MakeHisto2D(0.85,1.05,-1))->Draw("colz");

  AliTRDCalROC *ouiphy = calpad->GetCalRocRMS(det);
  TCanvas *cpouiloh = new TCanvas((const char*)"cpouiloh",(const char*)"cpouiloh",50,50,600,800);
  cpouiloh->Divide(2,1);
  cpouiloh->cd(1);
  ((TH1F *)ouiphy->MakeHisto1D(0.1,4.5,-1,10.0))->Draw();
  //((TH1F *)ouiphy->MakeHisto1D(0.01,0.45,-1))->Draw();
  cpouiloh->cd(2);
  ((TH2F *)ouiphy->MakeHisto2D(0.1,4.5,-1,10.0))->Draw("colz");
  //((TH2F *)ouiphy->MakeHisto2D(0.01,0.45,-1))->Draw("colz");
  
    

 
}
