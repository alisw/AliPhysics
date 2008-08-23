// Script to create reconstruction parameters and store them into CDB
// Author: Yuri Kharlov

/* $Id$ */

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"

#include "AliEMCALRecParam.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALSetRecParamCDB()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
  menu->AddButton("Help to run EMCAL CDB","Help()",
		  "Explains how to use EMCAL CDS menus");
  menu->AddButton("Set RecParam","SetRecParam()",
		  "Set clusterization parameters");
  menu->AddButton("Get RecParam","GetRecParam()",
		  "Get clusterization parameters");
  menu->AddButton("Exit","gApplication->Terminate(0)","Quit aliroot session");
  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\n\nSet reconstruction parameters and write them into ALICE OCDB. \nPress button \"Set RecParam\" to create an object AliEMCALRecParam \nand store it to OCDB. \nPress button \"Get RecParam\" to read reconstruction parameters from OCDB \nand print then to stdout.\n\n";
  printf(string);
}

//------------------------------------------------------------------------
void SetRecParam()
{
  // Create an object AliEMCALRecParam and store it to OCDB

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    =  999999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  DBFolder  ="local://LocalCDB";
  objFormat = "EMCAL reconstruction parameters";

  // Create reconstruction parameter object and set parameter values
  
  AliEMCALRecParam *recParamDB = new AliEMCALRecParam();
  //Clusterization
  recParamDB->SetClusteringThreshold(0.5);
  recParamDB->SetW0(4.5);
  recParamDB->SetMinECut(0.45);
  recParamDB->SetUnfold(kFALSE);
  recParamDB->SetLocMaxCut(0.03);

  //Track matching
  recParamDB->SetTrkCutX(6.0);
  recParamDB->SetTrkCutY(6.0);
  recParamDB->SetTrkCutZ(6.0);
  recParamDB->SetTrkCutR(10.0);
  recParamDB->SetTrkCutAlphaMin(-50.0);
  recParamDB->SetTrkCutAlphaMax( 50.0);
  recParamDB->SetTrkCutAngle(10000.0);      // i.e. exclude this for the moment

  //PID

  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      
      recParamDB->SetGamma(i,j,0.);
      recParamDB->SetHadron(i,j,0.);
      recParamDB->SetPiZero5to10(i,j, 0.);
      recParamDB->SetPiZero10to60(i,j,0.);
    }
  } 

  recParamDB->SetGamma(0,0,  0.038022);
  recParamDB->SetGamma(0,1, -0.0001883);
  recParamDB->SetGamma(0,2,  5.449e-06);
  
  recParamDB->SetGamma(1,0,  0.207313);
  recParamDB->SetGamma(1,1, -0.000978);
  recParamDB->SetGamma(1,2,  0.00001634);
  
  recParamDB->SetGamma(2,0,  0.043364);
  recParamDB->SetGamma(2,1, -0.0002048);
  recParamDB->SetGamma(2,2,  8.661e-06);
  recParamDB->SetGamma(2,3, -1.353e-07);
  
  recParamDB->SetGamma(3,0,  0.265004);
  recParamDB->SetGamma(3,1,  0.061298);
  recParamDB->SetGamma(3,2, -0.003203);
  recParamDB->SetGamma(3,3,  4.73e-05);
  
  recParamDB->SetGamma(4,0,  0.243579);
  recParamDB->SetGamma(4,1, -1.614e-05);
  
  recParamDB->SetGamma(5,0,  0.002942);
  recParamDB->SetGamma(5,1, -3.976e-05);
  
  recParamDB->SetHadron(0,0,  0.011945 / 3.);
  recParamDB->SetHadron(0,1,  0.000386 / 3.);
  recParamDB->SetHadron(0,2, -0.000014 / 3.);
  recParamDB->SetHadron(0,3,  1.336e-07 / 3.);
  
  recParamDB->SetHadron(1,0,  0.496544);
  recParamDB->SetHadron(1,1, -0.003226);
  recParamDB->SetHadron(1,2,  0.00001678);
  
  recParamDB->SetHadron(2,0,  0.144838);
  recParamDB->SetHadron(2,1, -0.002954);
  recParamDB->SetHadron(2,2,  0.00008754);
  recParamDB->SetHadron(2,3, -7.587e-07);
  
  recParamDB->SetHadron(3,0,  1.264461 / 7.);
  recParamDB->SetHadron(3,1,  0.002097 / 7.);
  
  recParamDB->SetHadron(4,0,  0.261950);
  recParamDB->SetHadron(4,1, -0.001078);
  recParamDB->SetHadron(4,2,  0.00003237);
  recParamDB->SetHadron(4,3, -3.241e-07);
  recParamDB->SetHadron(4,4,  0.);
  recParamDB->SetHadron(4,5,  0.);
  recParamDB->SetHadron(5,0,  0.010317);
  recParamDB->SetHadron(5,1,  0.);
  recParamDB->SetHadron(5,2,  0.);
  recParamDB->SetHadron(5,3,  0.);
  recParamDB->SetHadron(5,4,  0.);
  recParamDB->SetHadron(5,5,  0.);
  
  recParamDB->SetPiZero5to10(0,0, 0.009138);
  recParamDB->SetPiZero5to10(0,1, 0.0006377);
  
  recParamDB->SetPiZero5to10(1,0, 0.08);
  
  recParamDB->SetPiZero5to10(2,0, -0.061119);
  recParamDB->SetPiZero5to10(2,1,  0.019013);
  
  recParamDB->SetPiZero5to10(3,0,  0.2);
  
  recParamDB->SetPiZero5to10(4,0,  0.252044);
  recParamDB->SetPiZero5to10(4,1, -0.002315);
  
  recParamDB->SetPiZero5to10(5,0,  0.002942);
  recParamDB->SetPiZero5to10(5,1, -3.976e-05);
  
  recParamDB->SetPiZero10to60(0,0,  0.009138);
  recParamDB->SetPiZero10to60(0,1,  0.0006377);
  
  recParamDB->SetPiZero10to60(1,0,  1.272837);
  recParamDB->SetPiZero10to60(1,1, -0.069708);
  recParamDB->SetPiZero10to60(1,2,  0.001568);
  recParamDB->SetPiZero10to60(1,3, -1.162e-05);
  
  recParamDB->SetPiZero10to60(2,0,  0.139703);
  recParamDB->SetPiZero10to60(2,1,  0.003687);
  recParamDB->SetPiZero10to60(2,2, -0.000568);
  recParamDB->SetPiZero10to60(2,3,  1.498e-05);
  recParamDB->SetPiZero10to60(2,4, -1.174e-07);
  
  recParamDB->SetPiZero10to60(3,0, -0.826367);
  recParamDB->SetPiZero10to60(3,1,  0.096951);
  recParamDB->SetPiZero10to60(3,2, -0.002215);
  recParamDB->SetPiZero10to60(3,3,  2.523e-05);
  
  recParamDB->SetPiZero10to60(4,0,  0.249890);
  recParamDB->SetPiZero10to60(4,1, -0.000063);
  
  recParamDB->SetPiZero10to60(5,0,  0.002942);
  recParamDB->SetPiZero10to60(5,1, -3.976e-05);

  // raw signal fitting
  recParamDB->SetHighLowGainFactor(16.);
  recParamDB->SetOrderParameter(2);
  recParamDB->SetTau(2.35);
  recParamDB->SetNoiseThreshold(3);
  recParamDB->SetNPedSamples(5);
  
  // Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Yuri Kharlov");
  
  AliCDBId id("EMCAL/Config/RecParam",firstRun,lastRun);

  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage(DBFolder.Data());
  loc->Put(recParamDB, id, &md);
  recParamDB->Print();
}

//------------------------------------------------------------------------
void GetRecParam()
{
  // Read reconstruction parameters from OCDB

  TString DBFolder;

  //  DBFolder  ="local://LocalCDB";
  DBFolder  ="local://$ALICE_ROOT";
  Int_t runNumber = 0;

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("EMCAL/*",DBFolder.Data());

  AliCDBEntry* cdbEntry = AliCDBManager::Instance()->Get("EMCAL/Config/RecParam/",runNumber);
  if (cdbEntry == 0) {
    cerr << "No CDBEntry found at path "<<DBFolder.Data()<<"/"<<"EMCAL/Config/RecParam/"<<endl;
    return;
  }

  AliEMCALRecParam* recParam  = (AliEMCALRecParam*)cdbEntry->GetObject();
  if (recParam != 0) recParam->Print("");
  else
    cerr << "GetRecParam(): no AliEMCALRecParam object is found in OCDB" << endl;

}
