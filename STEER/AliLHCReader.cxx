/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class to read the file coming from DCS containing the information        //
//  from LHC. Everything is stored in a TMap, where:                         //
//  Key   --> DP name, as passed by LHC                                      // 
//  value --> TObjArray of AliDCSArray objects                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <time.h>

#include <TObjArray.h>
#include <TObjString.h>
#include <TObject.h>
#include <TString.h>
#include <TMap.h>
#include <TSystem.h>

#include "AliDCSArray.h"
#include "AliLHCReader.h"
#include "AliLog.h"

//--------------------------------------------------------------------------
AliLHCReader::AliLHCReader():
	TObject(),
	fStartTime(0),
	fEndTime(0)
{
	// default ctor 
}

//--------------------------------------------------------------------------
AliLHCReader::~AliLHCReader()
{
	//
	// dtor 
	//
}

//--------------------------------------------------------------------------
TMap* AliLHCReader::ReadLHCDP(TString filename)
{
	//
	// reading the file with the inputs
	//

       	if( gSystem->AccessPathName( filename.Data() ) ) {
		AliError(Form( "file (%s) not found", filename.Data() ) );
		return NULL;
	}

	ifstream *file = new ifstream ( filename.Data() );
	if (!*file) {
		AliError(Form("Error opening file (%s) !",filename.Data()));
		file->close();
		delete file;
		return NULL;
	}
	TMap* mapLHC = new TMap();
	mapLHC->SetOwner(1);
	TString strLine;
	while(strLine.ReadLine(*file)){
		// tokenize the line with tabs
		TObjArray* tokens = strLine.Tokenize(" \t");
		Int_t ntokens = tokens->GetEntriesFast();
		if (ntokens < 3){  
			// requiring at least the DP name, the format, and the number of entries
			delete tokens;
			continue;
		}
		TObjString* lhcDPname = (TObjString*)tokens->At(0);
		TString lhcDPtype = ((TObjString*)tokens->At(1))->String();
		AliDebug(2,Form("lhcDPtype = %s",lhcDPtype.Data()));
		TObjArray* typeTokens = lhcDPtype.Tokenize(":");
		if (typeTokens->GetEntriesFast() < 2 ){  
			// requiring the the type and the number of elements for each measurement
			AliError(Form("The format does not match the expected one, skipping the current line for DP = %s", lhcDPtype.Data()));
			delete typeTokens;
			continue;
		}
		TString type = ((TObjString*)typeTokens->At(0))->String();
		AliDebug(2,Form("type = %s",type.Data()));
      		Int_t nelements = (((TObjString*)typeTokens->At(1))->String()).Atoi();
		AliDebug(2,Form("nelements = %i",nelements));
		Int_t nentries = (((TObjString*)tokens->At(2))->String()).Atoi();
		AliDebug(2,Form("nentries = %i",nentries));
		Int_t nValuesPerEntry = nelements+1;
		Int_t nfixed = 3; // n. of fixed entries
		TObjArray* array = new TObjArray();
		array->SetOwner(1);
		for (Int_t ientry=0; ientry< nentries; ientry ++){
			Int_t indextime = nfixed+nValuesPerEntry*ientry+nelements;
			TString strTimestamp = ((TObjString*)tokens->At(indextime))->String();
			TObjArray* timeTokens = strTimestamp.Tokenize(".");
			if (timeTokens->GetEntriesFast() < 2 ){  
				// requiring both the seconds and the nseconds for the timestamp
				AliError(Form("The timestamp format does not match the expected one, skipping entry %d for DP = %s", ientry, lhcDPtype.Data()));
				continue;
			}
			time_t seconds = time_t((((TObjString*)timeTokens->At(0))->String()).Atoi());
			Int_t nseconds = Int_t((((TObjString*)timeTokens->At(1))->String()).Atoi());
			TTimeStamp* timestamp = new TTimeStamp(seconds, nseconds);
			AliDebug(2,Form("Timestamp in unix time = %d (s) = %d (ns)",timestamp->GetSec(),timestamp->GetNanoSec()));
			if (fStartTime!=0 && fEndTime!=0 && (fStartTime > (UInt_t)timestamp->GetSec() || fEndTime < (UInt_t)timestamp->GetSec())){
				// error in case the measurement is not within the data taking time interval
				AliError(Form("Timestamp for entry %d of DP %s not in [%d,%d]", ientry, lhcDPtype.Data(),fStartTime,fEndTime));
				continue;
			}
			if (type == "i"){
				Int_t* value = new Int_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = (((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atoi();
					AliDebug(2,Form("Value at index %d = %d",nfixed+ielement+ientry*nValuesPerEntry,value[ielement]));
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			}
			else if (type == "f"){
				Float_t* value = new Float_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = (((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atof();
					AliDebug(2,Form("Value at index %d = %f",nfixed+ielement+ientry*nValuesPerEntry,value[ielement]));
				} 
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			} 
			else if (type == "s"){
				TString* value = new TString[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = ((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String();
					AliDebug(2,Form("Value at index %d = %s",nfixed+ielement+ientry*nValuesPerEntry,value[ielement].Data()));
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			}
			else{
				AliError("Non-expected type");
				return NULL;
			} 
		}
		mapLHC->Add(lhcDPname,array);			
	}
	return mapLHC;
}

				






