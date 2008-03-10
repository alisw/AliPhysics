#include <iostream>
#include <fstream>

#include <TAlienCollection.h>
#include <TFile.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQADataMakerSteer.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliDAQ.h"

TString ClassName() { return "rawqa" ; } 

//________________________________qa______________________________________
void rawqa(const Int_t runNumber, Int_t maxFiles = 10, const char* year = "08") 
{	
	const char * kDefaultOCDBStorage = Form("alien://folder=/alice/data/20%s/LHC%sa/OCDB/", year, year) ; 
	
	UInt_t maxEvents = 99999 ;
	if ( maxFiles < 0 ) {
		maxEvents = TMath::Abs(maxFiles) ; 
		maxFiles = 99 ; 
	}
	AliLog::SetGlobalDebugLevel(0) ; 
	// connect to the grid 
	TGrid * grid = 0x0 ; 
       	grid = TGrid::Connect("alien://") ; 
		
	Bool_t detIn[AliDAQ::kNDetectors] = {kFALSE} ;
	char * detNameOff[AliDAQ::kNDetectors] = {"ITS", "ITS", "ITS", "TPC", "TRD", "TOF", "HMPID", "PHOS", "PHOS", "PMD", "MUON", "MUON", "FMD", "T0", "VZERO", "ZDC", "ACORDE", "TRG", "EMCAL", "DAQ_TEST", "HLT"} ; 
	// make the file name pattern year and run number
	TString pattern;
	pattern.Form("%9d",runNumber);
	pattern.ReplaceAll(" ", "0") ; 
	pattern.Prepend(year);
	pattern.Append("*0.root");

	// find the files associated to this run
	TGridResult * result = 0x0 ; 
	Bool_t local = kFALSE ; 
	if (grid) { // get the list of files from AliEn directly 
	  TString baseDir; 
	  baseDir.Form("/alice/data/20%s/",year);
	  result = grid->Query(baseDir, pattern) ;  
	} else {
	  TString collectionFile(pattern) ; 
	  collectionFile.Append(".xml") ; 
	  if ( gSystem->AccessPathName(collectionFile) == 0 ) { // get the list of files from an a-priori created collection file
	    TAlienCollection collection(collectionFile.Data(), maxFiles) ; 
	    result = collection.GetGridResult("", 0, 0); 
	  } else { // get the list of files from the local current directory 
	    local = kTRUE ; 
	    char line[100] ; 
	    sprintf(line, ".! ls %s*.root > tempo.txt", pattern.Data()) ; 
	    gROOT->ProcessLine(line) ;
	  }
	}
	AliLog::Flush();
	ifstream in ; 
	if (local) 
		in.open("tempo.txt", ifstream::in) ; 

	AliCDBManager* man = AliCDBManager::Instance();
	man->SetDefaultStorage(kDefaultOCDBStorage) ;  
	AliQA::SetQARefStorage("alien://folder=/alice/QA/2008/") ; 
	man->SetSpecificStorage("*", AliQA::GetQARefStorage());
	AliQADataMakerSteer qas ; 
	TString detectors  = ""; 
	TString detectorsW = ""; 
	UShort_t file = 0 ; 
	UShort_t filesProcessed = 0 ; 
	for ( file = 0 ; file < maxFiles ; file++) {
		TString fileName ; 
		if ( local) {
			in >> fileName ;
		} 
		else 
			fileName = result->GetKey(file, "turl");
	       	if ( fileName == "" )  
		        break ;
		if ( fileName.Contains("tag") )
			continue; 
                filesProcessed++ ;
		char input[200] ; 
		if (local) 
			sprintf(input, "%s", fileName.Data()) ; 
		else 
			sprintf(input, "%s", result->GetKey(file, "turl")); 
		AliInfo(Form("Proccessing file # %d --> %s", file, input)) ;
		AliLog::Flush();
		// check which detectors are present 
		AliRawReader * rawReader = new AliRawReaderRoot(input);
		while ( rawReader->NextEvent() ) {
			man->SetRun(rawReader->GetRunNumber());
			AliLog::Flush();
			UChar_t * data ; 
			while (rawReader->ReadNextData(data)) {
				Int_t detID = rawReader->GetDetectorID();
				if (detID < 0 || detID >= AliDAQ::kNDetectors) {
					AliError("Wrong detector ID! Skipping payload...");
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
		AliInfo(Form("Current = %d Max = %d", qas.GetCurrentEvent(), maxEvents)) ; 
		AliLog::Flush();
		if ( qas.GetCurrentEvent() > maxEvents) 
			break ;
		// TEMPORARY REMOVAL OF TRD!!!
		detectors.ReplaceAll("TRD", "") ;
		// TEMPORARY REMOVAL OF TRD!!!
		if ( !detectors.IsNull() ) {
			qas.SetMaxEvents(maxEvents) ; 	
			detectorsW = qas.Run(detectors, rawReader) ;
			qas.Reset() ;
		} else {
			AliError("No valid detectors found") ; 
		} 
		delete rawReader ;
	}
	AliLog::Flush();
	qas.Merge(runNumber) ; 
	
	AliLog::Flush();
	// The summary 
	AliInfo(Form("\n\n********** Summary for run %d **********", runNumber)) ; 
	printf("     detectors present in the run        : %s\n", detectors.Data()) ; 
	printf("     detectors present in the run with QA: %s\n", detectorsW.Data()) ; 
	printf("     number of files processed           : %d\n", filesProcessed) ; 
	TFile * qaResult = TFile::Open(AliQA::GetQAResultFileName()) ; 
	AliQA * qa = dynamic_cast<AliQA *>(qaResult->Get("QA")) ; 
	for (Int_t index = 0 ; index < AliQA::kNDET ; index++)
		if (detectorsW.Contains(AliQA::GetDetName(AliQA::DETECTORINDEX(index)))) 
			qa->ShowStatus(AliQA::DETECTORINDEX(index)) ;
}
