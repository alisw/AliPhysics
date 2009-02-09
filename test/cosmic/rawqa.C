//#include <iostream>
//#include <fstream>
//
//#include <TAlienCollection.h>
//#include <TFile.h>
//#include <TGrid.h>
//#include <TGridResult.h>
//#include <TMath.h>
//#include <TROOT.h>
//#include <TString.h>
//#include <TSystem.h>
//
//#include "AliCDBManager.h"
//#include "AliDAQ.h"
//#include "AliLog.h"
//#include "AliQA.h"
//#include "AliQADataMakerSteer.h"
//#include "AliRawReader.h"
//#include "AliRawReaderRoot.h"
//#include "AliTRDrawStreamBase.h"
//#include "AliGeomManager.h"

TString ClassName() { return "rawqa" ; } 

//________________________________qa______________________________________
void rawqa(const char * filename) 
{	
  // retrieve evironment variables
  const char * year = gSystem->Getenv("YEAR") ; 
  const TString baseDir(gSystem->Getenv("BASEDIR")) ;
  const Int_t runNumber = atoi(gSystem->Getenv("RUNNUM")) ; 

  // build the default storage of OCDB 
  TString sfilename(filename) ; 
  sfilename = sfilename.Strip() ; 
  sfilename = sfilename.Strip(TString::kLeading) ;
  sfilename.Prepend("alien://") ; 

  baseDir.Append("/") ; 
  TString temp(filename) ;
  temp = temp.ReplaceAll(baseDir, "") ; 
  temp = temp.Strip(TString::kLeading, '/') ; 
  baseDir.Append(temp(0, temp.Index("/"))) ;  
  char * kDefaultOCDBStorage = Form("alien://folder=%s/OCDB/", baseDir.Data()) ; 

  // set the location of reference data 
  //AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), year)) ;  
  AliQA::SetQARefStorage("local://$ALICE_ROOT/OCDB") ;
	
  AliLog::SetGlobalDebugLevel(0) ; 
	
  // connect to the grid 
  TGrid * grid = 0x0 ; 
  grid = TGrid::Connect("alien://") ; 
	
  Bool_t detIn[AliDAQ::kNDetectors] = {kFALSE} ;
  char * detNameOff[AliDAQ::kNDetectors] = {"ITS", "ITS", "ITS", "TPC", "TRD", "TOF", "HMPID", "PHOS", "PHOS", 
          "PMD", "MUON", "MUON", "FMD", "T0", "VZERO", "ZDC", "ACORDE", "TRG", 
          "EMCAL", "DAQ_TEST", "HLT"} ; 
	
  AliQADataMakerSteer qas("rec") ; 
  TString detectors  = ""; 
  TString detectorsW = ""; 
  UShort_t eventsProcessed = 0 ; 
  UShort_t filesProcessed  = 1 ; 
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(kDefaultOCDBStorage) ;  
  man->SetRun(runNumber) ; 
  AliGeomManager::LoadGeometry();
  printf("INFO: Proccessing file %s\n", filename) ;
  // check which detectors are present 
  AliRawReader * rawReader = new AliRawReaderRoot(sfilename);
  AliTRDrawStreamBase::SetRawStreamVersion("TB");
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
	  //qas->SetMaxEvents(1) ; 
    qas.SetTasks(Form("%d", AliQA::kRAWS)); 			
    detectorsW = qas.Run(detectors, rawReader) ;
    qas.Reset() ;
  } else {
       	  printf("ERROR: No valid detectors found") ; 
  } 
  delete rawReader ;
  eventsProcessed += qas.GetCurrentEvent() ; 
	//qas.Merge(runNumber) ; 
	
	// The summary 
  printf("\n\n********** Summary for run %d **********", runNumber) ; 
  printf("     data file                           : %s\n", filename);
  printf("     detectors present in the run        : %s\n", detectors.Data()) ; 
  printf("     detectors present in the run with QA: %s\n", detectorsW.Data()) ; 
  printf("     number of files/events processed    : %d/%d\n", filesProcessed, eventsProcessed) ; 
  TFile * qaResult = TFile::Open(AliQA::GetQAResultFileName()) ; 
  if ( qaResult ) {
    AliQA * qa = dynamic_cast<AliQA *>(qaResult->Get(AliQA::GetQAName())) ; 
    if ( qa) {
      for (Int_t index = 0 ; index < AliQA::kNDET ; index++)
        if (detectorsW.Contains(AliQA::GetDetName(index))) 
          qa->Show((AliQA::GetDetIndex(AliQA::GetDetName(index)))) ;
    } else {
      printf("ERROR: %s not found in %s !", AliQA::GetQAName(), AliQA::GetQAResultFileName()) ; 
    }
  } else {
    printf("ERROR: %s has not been produced !", AliQA::GetQAResultFileName()) ; 
  }
}
