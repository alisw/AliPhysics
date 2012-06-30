// $Id: AliHLTSamplePreprocessor.cxx 23039 2007-12-13 20:53:02Z richterm $

/**************************************************************************
 * This file is property of and copyright by the                          * 
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

/// @file   AliHLTSamplePreprocessor.cxx
/// @author Kenneth Aamodt, Sebastian Bablok
/// @date   2007-12-06
/// @brief  HLT Preprocessor plugin for the AliHLTComp library
/// 

#include "AliHLTSamplePreprocessor.h"
#include "AliPreprocessor.h"

#include <AliCDBMetaData.h>
#include <TObjString.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>

ClassImp(AliHLTSamplePreprocessor)

AliHLTSamplePreprocessor::AliHLTSamplePreprocessor()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTSamplePreprocessor::fgkTempHistoFileName = "HLTTemperatureHistograms.root";

AliHLTSamplePreprocessor::~AliHLTSamplePreprocessor()
{
  // see header file for function documentation
}

void AliHLTSamplePreprocessor::Initialize(Int_t /*run*/, UInt_t /*startTime*/, 
					  UInt_t /*endTime*/)
{
  // see header file for function documentation
}


UInt_t AliHLTSamplePreprocessor::Process(TMap* /*dcsAliasMap*/)
{
  // see header file for function documentation
  UInt_t retVal = 0;
  if (GetTempHisto() != 0) {
    // unable to fetch the temperature histogram from HLT
    retVal = 1;
    // but if set to 1, then also file from GetHuffmanTables won't be saved !!!
  }

  return retVal;
}

UInt_t AliHLTSamplePreprocessor::GetTempHisto()
{
  // see header file for function documentation

	UInt_t retVal = 0;
    // get Temp Histogram map
    TList* hltList = GetFileSources(AliPreprocessor::kHLT, fgkTempHistoFileName);
    if (!hltList) {
        Log("Missing list for the HLT");
   
        return 1;
    }

	if (hltList->GetSize() == 0) {
		Log("No Temperature histogram produced inside the HLT by a DA for this run.");
		return retVal; 
		// return no error -> DA might not have run, but other file shall be saved.
	} else if (hltList->GetSize() > 1) {
        Log(Form("Problem on the size of the list: %d (HLT)", hltList->GetSize()));
        return 0; // might have to be changed, when there will be more than one histogram file
    }

    TObjString* location = (TObjString*) hltList->At(0);
    if (location == 0) {
        Log("Error in location HLT list.");
        return 0;
    }
    TString localFileName = GetFile(AliPreprocessor::kHLT, fgkTempHistoFileName,
                location->String().Data());
	if (!(localFileName.Length() > 0)) {
		Log("Local file name for Temperature Histogram has zero length.");
		return 1;
	}
	
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
    }

	return retVal;
}
