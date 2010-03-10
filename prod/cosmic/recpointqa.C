#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliQAv1.h"
#include "AliQAManager.h"
#include "AliESDEvent.h"
#include "AliRecoParam.h"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

const char * ClassName() 
{
  return "recpointqa" ; 
}

void recpointqa() 
{
  Int_t run = atoi(gSystem->Getenv("RUNNUM")) ; 
  TString detectors(gSystem->Getenv("DET")) ; 
  
  printf("Set the CDB storage location...\n");
  AliCDBManager * man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://"); //set alien right path for CDB
  man->SetRun(run); 

  AliLog::SetGlobalDebugLevel(0) ;  

  if ( ! AliGeomManager::GetGeometry() ) AliGeomManager::LoadGeometry();
  AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QAref");

  AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kQAMODE) ;     
  qam->SetEventSpecie(qam->GetEventSpecieFromESD()) ;
  qam->SetWriteExpert();    
  TString detectorsW = qam->Run(detectors,AliQAv1::kRECPOINTS);

  // The summary 
  printf("\n\n********** Summary for run %d **********\n", run) ; 
  printf("     detectors present in the run        : %s\n", detectors.Data()) ; 
  printf("     detectors present in the run with QA: %s\n", detectorsW.Data()) ; 
  printf("     %d events have been processed \n", qam->GetCurrentEvent());
}
