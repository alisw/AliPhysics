// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Sebastian Bablok <Sebastian.Bablok@ift.uib.no>        *
 *                  Kenneth Aamodt                                        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * @file   AliHLTPreprocessor.cxx
 * @author Kenneth Aamodt, Sebastian Bablok
 * @date   2007-12-06
 * @brief  Implementation of the HLT preprocessor (used by the Offline Shuttle) 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPreprocessor.h"

//#include <AliCDBMetaData.h>
//#include <AliCDBEntry.h>

#include <AliCDBMetaData.h>


#include <TObjString.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>


ClassImp(AliHLTPreprocessor)

const Int_t AliHLTPreprocessor::fgkHuffmanTablesNum = 6;

const char* AliHLTPreprocessor::fgkHLTPreproc = "HLT";

const char* AliHLTPreprocessor::fgkHuffmanFileBase = "huffmanData_";

const char* AliHLTPreprocessor::fgkHuffmanFileDetector = "TPC_";	// at the moment only one

const char* AliHLTPreprocessor::fgkTempHistoFileName = "HLTTemperatureHistograms.root";

AliHLTPreprocessor::AliHLTPreprocessor(AliShuttleInterface* shuttle) 
  :
  AliPreprocessor(fgkHLTPreproc, shuttle),
  fRun(0),
  fStartTime(0),
  fEndTime(0)
{
// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


AliHLTPreprocessor::~AliHLTPreprocessor() {
// see header file for function documentation
}

void AliHLTPreprocessor::Initialize(Int_t run, UInt_t startTime, 
			UInt_t endTime) {
// see header file for function documentation
	fRun = run;
	fStartTime = startTime;
	fEndTime = endTime;

	TString msg("Preprocessor for HLT initialized for run: ");
	msg += run;
//	Log(msg.Data());
}


UInt_t AliHLTPreprocessor::Process(TMap* dcsAliasMap) {
// see header file for function documentation
	UInt_t retVal = 0;
//	const char* localFileName = 0;

	if (!GetHLTStatus()) {
		return 0;
	}

	// get Huffman tables
	for (Int_t i = 0; i < fgkHuffmanTablesNum; i++) {
		TString runNumberString;
		runNumberString.Form("%08d", fRun);
		TString filename(fgkHuffmanFileBase);
		filename += fgkHuffmanFileDetector;
		filename += runNumberString;
		filename += "_0x23000";
		filename += i;
		filename += "0";
		filename += i;
		filename += ".root";

		//spec 0x23000Y0Y -> huffmanData_<detector>_<runnumber>_<specification>.root
		TList* HLTlist = GetFileSources(kHLT, filename.Data());	
		if (!HLTlist) {
			Log("Missing list for the HLT");
		    continue;
		}

		if (HLTlist->GetSize() != 1) {
			Log(Form("Problem on the size of the list: %d (HLT)", 
						HLTlist->GetSize()));
			continue;
  		}

		TObjString* location = (TObjString*) HLTlist->At(0);
		if (location == 0) {
			Log("Error in location HLT list.");
			continue;
		}
		TString localFileName = GetFile(kHLT, filename.Data(), 
					location->String().Data()); 
		
/*		
		TFile localFile(localFileName);
		
		AliCDBMetaData meta("Jennifer Wagner");
		TString name("huffmanData_");
		name += kDetector;
		name += "Patch_";
		name += i;

		if (!(Store("CalibTPC", name.Data(), (TObject*) &localFile, &meta, 0, kTRUE))) {
*/
		if (!(StoreReferenceFile(localFileName.Data(), filename.Data()))) {
       		TString msg("Storing of object '");
			msg += filename;
			msg += "' to Reference Storage failed!";
        	Log(msg.Data());
			retVal = 1; 
			// I think this is then really an error and should return an error code
		}
	}
	
	// get Temp Histogram map
	TList* HLTlist = GetFileSources(kHLT, fgkTempHistoFileName);
	if (!HLTlist) {
    	Log("Missing list for the HLT");
		return 0;
	}

	if (HLTlist->GetSize() != 1) {
		Log(Form("Problem on the size of the list: %d (HLT)", HLTlist->GetSize()));
		return 0;
	}

	TObjString* location = (TObjString*) HLTlist->At(0);
	if (location == 0) {
		Log("Error in location HLT list.");
		return 0;
	}
	TString localFileName = GetFile(kHLT, fgkTempHistoFileName, 
				location->String().Data());
/*
	TFile localFile(localFileName);
	AliCDBMetaData meta("Sebastian Bablok");

	if (!(Store("Calib", kTempHistoFileName, (TObject*) &localFile, &meta, 0, kTRUE))) {
*/
	if (!(StoreReferenceFile(localFileName.Data(), fgkTempHistoFileName))) {
		TString msg("Storing of object '");
		msg += fgkTempHistoFileName;
		msg += "' to Reference Storage failed!";
		Log(msg.Data());
		retVal = 1;
		// I think this is then really an error and should return an error code
	}

	return retVal;
}


Bool_t AliHLTPreprocessor::ProcessDCS() {
// see header file for function documentation
    return kFALSE;
}

