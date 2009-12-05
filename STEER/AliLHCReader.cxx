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
#include <TError.h>

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
	//	mapLHC->SetOwner(1);
	TString strLine;
	Int_t lhcEntries;
	Int_t nBlocks = 0;
	//	TObjArray** array;
	Int_t nline =0;
	while(strLine.ReadLine(*file)){
		nline++;
		AliDebug(3,Form("line n. = %d",nline));
		// tokenize the line with tabs
		//if (strLine.BeginsWith("=")) continue;
		//if (!strLine.CompareTo("END_OF_BLOCK")) {
		if (strLine.Contains("END_OF_BLOCK")) {
			AliInfo("END_OF_BLOCK");
			nBlocks++;
			continue;
		}
		TObjArray* tokens = strLine.Tokenize("\t");
		Int_t ntokens = tokens->GetEntriesFast();
		AliDebug(3,Form("Number of tokens = %d",ntokens));
		if (ntokens == 2 && !(((TObjString*)tokens->At(0))->String()).CompareTo("LHC_ENTRIES")){
			lhcEntries = (((TObjString*)tokens->At(1))->String()).Atoi();
			AliInfo(Form("LHC entries = %d",lhcEntries));
			AliDebug(3,Form("LHC entries = %d",lhcEntries));
			continue;
		}
		if (ntokens == 1 && !(((TObjString*)tokens->At(0))->String()).CompareTo("END_OF_DATA")){
			AliInfo("End of file reached");
			break;
		}
		if (ntokens < 4){  
			AliInfo(Form("Wrong number of tokens --> # tokens = %d at line %d",ntokens,nline));
			// requiring at least the index of the DP, the DP name, the format, and the number of entries
			delete tokens;
			continue;
		}
		Int_t lhcDPindex = (((TObjString*)tokens->At(0))->String()).Atoi();
		AliDebug(2,Form("lhcDPindex = %d",lhcDPindex));
		TObjString* lhcDPname = (TObjString*)tokens->At(1);
		TString lhcDPtype = ((TObjString*)tokens->At(2))->String();
		AliInfo(Form("lhcDPname = %s",(lhcDPname->String()).Data()));
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
		Int_t nentries = (((TObjString*)tokens->At(3))->String()).Atoi();
		AliDebug(2,Form("nentries = %i",nentries));
		Int_t nValuesPerEntry = nelements+1;
		Int_t nfixed = 4; // n. of fixed entries
		TObjArray* array;
		if (mapLHC->GetValue(lhcDPname)==0x0){
			array = new TObjArray();
			array->SetOwner(1);
		}
		else{
			TPair* p = mapLHC->RemoveEntry(lhcDPname);
			array = (TObjArray*)p->Value();
			AliDebug(2,Form("entry found! --> %p",array));
		}
					
		for (Int_t ientry=0; ientry< nentries; ientry ++){
			Int_t indextime = nfixed+nValuesPerEntry*ientry+nelements;
			TString strTimestamp = ((TObjString*)tokens->At(indextime))->String();
			//			TObjArray* timeTokens = strTimestamp.Tokenize(".");
			//if (timeTokens->GetEntriesFast() < 2 ){  
			//	// requiring both the seconds and the nseconds for the timestamp
			//		AliError(Form("The timestamp format does not match the expected one, skipping entry %d for DP = %s", ientry, lhcDPtype.Data()));
			//	continue;
			//}
			//time_t seconds = time_t((((TObjString*)timeTokens->At(0))->String()).Atoi());
			//Int_t nseconds = Int_t((((TObjString*)timeTokens->At(1))->String()).Atoi());
			//TTimeStamp* timestamp = new TTimeStamp(seconds, (Int_t)(nseconds*1E8));
			Double_t timestamp = strTimestamp.Atof();
			AliDebug(2,Form("Timestamp in unix time = %f (s)",timestamp));
			if (fStartTime!=0 && fEndTime!=0 && (fStartTime > timestamp || fEndTime < timestamp)){
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
			else if (type == "b"){
				Bool_t* value = new Bool_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = Bool_t((((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atoi());
					AliDebug(2,Form("Value at index %d = %d",nfixed+ielement+ientry*nValuesPerEntry,Int_t(value[ielement])));
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			}
			else if (type == "f"){ // the floats should be considered as doubles
				Double_t* value = new Double_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					TString tempstr = (TString)(((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String());
					value[ielement] = (((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atof();
					AliDebug(2,Form("Value at index %d = %f from string %s",nfixed+ielement+ientry*nValuesPerEntry,value[ielement],tempstr.Data()));
				} 
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			} 
			else if (type == "s"){
				TObjArray* value = new TObjArray();
				for (Int_t ielement=0; ielement<nelements; ielement++){
					TObjString* strobj = ((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry));
					AliDebug(2,Form("Value at index %d = %s",nfixed+ielement+ientry*nValuesPerEntry,(strobj->String()).Data()));
					value->Add(strobj);
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
			}
			else{
				AliError(Form("Non-expected type %s",type.Data()));
				return NULL;
			} 
		}
		mapLHC->Add(lhcDPname,array);			
	}
	//mapLHC->Print();
	return mapLHC;
}

				






