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
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAManager.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
//#include "AliTRDrawStreamBase.h"
#include "AliGeomManager.h"

TString ClassName() { return "rawqa" ; } 

//________________________________qa______________________________________
void rawqa(Char_t * filename, Int_t run) 
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
  
  // set the location of reference data 
  //AliQAv1::SetQARefStorage(Form("%s%s/", AliQAv1::GetQARefDefaultStorage(), year)) ;  
  AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QAref") ;
	
  AliLog::SetGlobalDebugLevel(0) ; 
	
  	
  Bool_t detIn[AliDAQ::kNDetectors] = {kFALSE} ;
  const char * detNameOff[AliDAQ::kNDetectors] = {"ITS", "ITS", "ITS", "TPC", "TRD", "TOF", "HMPID", "PHOS", "PHOS", 
          "PMD", "MUON", "MUON", "FMD", "T0", "VZERO", "ZDC", "ACORDE", "TRG", 
          "EMCAL", "DAQ_TEST", "HLT"} ; 
	
  AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ; 
  qam->SetEventSpecie(AliRecoParam::kCosmic) ; 
  AliQAv1::Instance()->SetEventSpecie(AliRecoParam::kCosmic) ; 
  TString detectors  = ""; 
  TString detectorsW = ""; 
  UShort_t eventsProcessed = 0 ; 
  UShort_t filesProcessed  = 1 ; 
  man->SetRun(run);
  AliGeomManager::LoadGeometry();
  printf("INFO: Proccessing file %s\n", filename) ;
  // check which detectors are present 
  TString alienName = Form("alien://%s", filename) ;
  AliRawReader * rawReader = new AliRawReaderRoot(alienName.Data());
  //AliTRDrawStreamBase::SetRawStreamVersion("TB");
  while ( rawReader->NextEvent() ) {
    man->SetRun(rawReader->GetRunNumber());
    AliLog::Flush();
    UChar_t * data ; 
    while (rawReader->ReadNextData(data)) {
      Int_t detID = rawReader->GetDetectorID();
      if (detID < 0 || detID >= AliDAQ::kNDetectors) {
        printf("INFO: Wrong detector ID! Skipping payload...\n");
        continue;
      }
      detIn[detID] = kTRUE ; 
    }
    for (Int_t detID = 0; detID < AliDAQ::kNDetectors ; detID++) {
      if (detIn[detID]) {
        if ( ! detectors.Contains(detNameOff[detID]) ) {
          detectors.Append(detNameOff[detID]) ;
          detectors.Append(" ") ;
        }
      }
    }
    if ( !detectors.IsNull() )
      break ; 
  }
  if ( !detectors.IsNull() ) {
	  qam->SetMaxEvents(-1) ; 
    qam->SetTasks(Form("%d", AliQAv1::kRAWS)); 			
    detectorsW = qam->Run(detectors, rawReader) ;
    qam->Reset() ;
  } else {
       	  printf("ERROR: No valid detectors found") ; 
  } 
  delete rawReader ;
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
