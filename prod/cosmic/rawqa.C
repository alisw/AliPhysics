#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliDAQ.h"
#include "AliGRPObject.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAManager.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliGeomManager.h"
#include "AliRecoParam.h"

TString ClassName() { return "rawqa" ; } 

//________________________________qa______________________________________
void rawqa(Char_t * filename, Int_t run, AliRecoParam::EventSpecie_t es=AliRecoParam::kDefault) 
{	
//  TGrid * grid = TGrid::Connect("alien://") ; 
//  TString filename ; 
//  if (grid) {
//    filename = Form("alien:///alice/data/2009/LHC09c/000085034/raw/09000085034023.40.root") ;  
//  } else {
//    filename = "raw.root" ;       
//  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
//  man->SetDefaultStorage("local://cdb2");
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliCDBEntry *  entry    = man->Get("GRP/GRP/Data");  
  if (!entry) 
    return ;
  AliGRPObject * fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if (!fGRPData) {
    printf("ERROR: No GRP entry found in OCDB") ; 
    return ; 
  }
  Int_t activeDetectors = fGRPData->GetDetectorMask();
  const char * detNameOff[AliDAQ::kNDetectors] = { 
    "ITS",
    "ITS",
    "ITS",
    "TPC",
    "TRD",
    "TOF",
    "HMPID",
    "PHOS",
    "CPV",
    "PMD",
    "MUON",
    "MUON",
    "FMD",
    "T0",
    "VZERO", // Name to be changed to V0 ?
    "ZDC",
    "ACORDE",
    "TRG",
    "EMCAL",
    "DAQ_TEST",
    "HLT"
  } ; 
  TString detectors  = ""; 
  TString detectorsW ; 
  for(Int_t iDet = 0; iDet < (AliDAQ::kNDetectors-1); iDet++) {
    if ((activeDetectors >> iDet) & 0x1) {
      if (!detectors.Contains(detNameOff[iDet])) {
        detectors +=detNameOff[iDet] ; 
        detectors += " " ;  
      }
    }
  }
  AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QAref") ;

  AliLog::SetGlobalDebugLevel(0) ; 
	
  AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ; 
  qam->SetEventSpecie(AliRecoParam::kCosmic) ; 
  AliQAv1::Instance()->SetEventSpecie(es) ; 
//  TString detectorsW = ""; 
  UShort_t eventsProcessed = 0 ; 
  UShort_t filesProcessed  = 1 ; 
  AliGeomManager::LoadGeometry();
  printf("INFO: Proccessing detectors %s from file %s\n", detectors.Data(), filename) ;
  if ( !detectors.IsNull() ) {
	  qam->SetMaxEvents(-1) ; 
    qam->SetTasks(Form("%d", AliQAv1::kRAWS)); 		
    AliRawReader * rawReader = new AliRawReaderRoot(Form("alien://%s", filename));
    detectorsW = qam->Run(detectors, rawReader) ;
    qam->Reset() ;
  } else {
    printf("ERROR: No valid detectors found") ; 
  } 
  eventsProcessed += qam->GetCurrentEvent() ; 
	//qam->Merge(run) ; 
	
	// The summary 
  printf("\n\n********** Summary for run %d **********", run) ; 
  printf("     data file                           : %s\n", filename);
  printf("     detectors present in the run        : %s\n", detectors.Data()) ; 
  printf("     detectors present in the run with QA: %s\n", detectorsW.Data()) ; 
  printf("     number of files/events processed    : %d/%d\n", filesProcessed, eventsProcessed) ; 
  //const char * qaFileName = AliQAv1::GetQAResultFileName() ;
  //  TFile * qaResult = TFile::Open(qaFileName) ; 
  TFile * qaResult = TFile::Open("QA.root") ; 
  if ( qaResult ) {
    AliQAv1 * qa = dynamic_cast<AliQAv1 *>(qaResult->Get(AliQAv1::GetQAName())) ; 
    if ( qa) {
      for (Int_t index = 0 ; index < AliQAv1::kNDET ; index++)
        if (detectorsW.Contains(AliQAv1::GetDetName(index))) 
          qa->Show((AliQAv1::GetDetIndex(AliQAv1::GetDetName(index)))) ;
    } else {
      printf("ERROR: %s not found in %s !", AliQAv1::GetQAName(), AliQAv1::GetQAResultFileName()) ; 
    }
  } else {
    printf("ERROR: %s has not been produced !", AliQAv1::GetQAResultFileName()) ; 
  }
}
