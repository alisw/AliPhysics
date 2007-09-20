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
  recParamDB->SetClusteringThreshold(0.5);
  recParamDB->SetW0(4.5);
  recParamDB->SetMinECut(0.45);

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

  DBFolder  ="local://LocalCDB";
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
