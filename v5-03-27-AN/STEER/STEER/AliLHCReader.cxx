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
  gSystem->ExpandPathName(filename);
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
	mapLHC->SetOwnerKeyValue();
	TString strLine;
	Int_t lhcEntries;
	Int_t nBlocks = 0;
	Int_t nline =0;
	while(strLine.ReadLine(*file)){
		nline++;
		AliDebug(3,Form("line n. = %d",nline));
		// tokenize the line with tabs
		//if (strLine.BeginsWith("=")) continue;
		//if (!strLine.CompareTo("END_OF_BLOCK")) {
		if (strLine.Contains("END_OF_BLOCK")) {
			AliDebug(2,"END_OF_BLOCK");
			nBlocks++;
			continue;
		}
		TObjArray* tokens = strLine.Tokenize("\t");
		Int_t ntokens = tokens->GetEntriesFast();
		//
		if ( strLine.Contains("LHC_ENTRIES ") ) {
		  // RS: special treatment for "LHC_ENTRIES N" record, which is space (and not tab) separated)
		  delete tokens;
		  tokens = strLine.Tokenize(" ");
		  ntokens = tokens->GetEntriesFast();
		}
		//
		AliDebug(3,Form("Number of tokens = %d",ntokens));
		if (ntokens == 2 && !(((TObjString*)tokens->At(0))->String()).CompareTo("LHC_ENTRIES")){
			lhcEntries = (((TObjString*)tokens->At(1))->String()).Atoi();
			AliInfo(Form("LHC entries = %d",lhcEntries));
			AliDebug(3,Form("LHC entries = %d",lhcEntries));
			continue;
		}
		if (ntokens == 1 && !(((TObjString*)tokens->At(0))->String()).CompareTo("END_OF_DATA")){
			AliDebug(2,"End of file reached");
			delete tokens;
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
		AliDebug(2,Form("lhcDPname = %s",(lhcDPname->String()).Data()));
		AliDebug(2,Form("lhcDPtype = %s",lhcDPtype.Data()));
		TObjArray* typeTokens = lhcDPtype.Tokenize(":");
		if (typeTokens->GetEntriesFast() < 2 ){  
			// requiring the the type and the number of elements for each measurement
			AliError(Form("The format does not match the expected one, skipping the current line for DP = %s", lhcDPtype.Data()));
			delete typeTokens;
			delete tokens;
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
			mapLHC->Add(new TObjString(lhcDPname->String()),array);			
		}
		else{
			array = (TObjArray*)mapLHC->GetValue(lhcDPname);
			AliDebug(2,Form("entry found! --> %p",array));
		}
					
		for (Int_t ientry=0; ientry< nentries; ientry ++){
			Int_t indextime = nfixed+nValuesPerEntry*ientry+nelements;
			TString strTimestamp = ((TObjString*)tokens->At(indextime))->String();
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
				delete[] value;
			}
			else if (type == "b"){
				Bool_t* value = new Bool_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = Bool_t((((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atoi());
					AliDebug(2,Form("Value at index %d = %d",nfixed+ielement+ientry*nValuesPerEntry,Int_t(value[ielement])));
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
				delete[] value;
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
				delete[] value;
			} 
			else if (type == "s"){
				TObjArray* value = new TObjArray();
				value->SetOwner(1);
				for (Int_t ielement=0; ielement<nelements; ielement++){
				  TObjString* strobj = (new TObjString(((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()));
					AliDebug(2,Form("Value at index %d = %s",nfixed+ielement+ientry*nValuesPerEntry,(strobj->String()).Data()));
					value->Add(strobj);
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
				delete value;
			}
			else{
				AliError(Form("Non-expected type %s",type.Data()));
				delete typeTokens;
				delete tokens;	
				file->close();
				delete file;	
				return NULL;
			} 
		}
		delete typeTokens;
		delete tokens;
	}
	file->close();
	delete file;
	return mapLHC;
}

//--------------------------------------------------------------------------
TObjArray* AliLHCReader::ReadSingleLHCDP(TString filename, TString alias)
{
	//
	// reading the file with the inputs for the selected alias
	// returning the TObjArray containing the information only for the current alias
	//
  gSystem->ExpandPathName(filename);
       	if( gSystem->AccessPathName( filename.Data() ) ) {
		AliError(Form( "file (%s) not found", filename.Data() ) );
		return NULL;
	}

	TString selection = gSystem->GetFromPipe(Form("grep -P '^\\d+\\s+%s+\\s' %s",alias.Data(), filename.Data()));

	if (selection.Length() == 0) {
		AliError(Form("Alias %s not fouond in LHC Data file, returning a null pointer",alias.Data()));
		return NULL;
	}

	Int_t nline =0;

	TObjArray* tokenslines = selection.Tokenize("\n");
	Int_t ntokenslines = tokenslines->GetEntriesFast();
	AliDebug(3,Form("Number of tokenslines = %d",ntokenslines));

	TObjArray* array = new TObjArray(); // array to be returned
	array->SetOwner(1);

	for (Int_t iline=0; iline<ntokenslines; iline++){
		TString strLine = ((TObjString*)tokenslines->At(iline))->String();
		AliDebug(4,Form("***************** line = %s\n",strLine.Data()));
		TObjArray* tokens = strLine.Tokenize("\t");
		Int_t ntokens = tokens->GetEntriesFast();
		AliDebug(3,Form("Number of tokens = %d",ntokens));
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
		AliDebug(2,Form("lhcDPname = %s",(lhcDPname->String()).Data()));
		AliDebug(2,Form("lhcDPtype = %s",lhcDPtype.Data()));
		TObjArray* typeTokens = lhcDPtype.Tokenize(":");
		if (typeTokens->GetEntriesFast() < 2 ){  
			// requiring the the type and the number of elements for each measurement
			AliError(Form("The format does not match the expected one, skipping the current line for DP = %s", lhcDPtype.Data()));
			delete typeTokens;
			delete tokens;
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
		for (Int_t ientry=0; ientry< nentries; ientry ++){
			Int_t indextime = nfixed+nValuesPerEntry*ientry+nelements;
			TString strTimestamp = ((TObjString*)tokens->At(indextime))->String();
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
				delete[] value;
			}
			else if (type == "b"){
				Bool_t* value = new Bool_t[nelements];
				for (Int_t ielement=0; ielement<nelements; ielement++){
					value[ielement] = Bool_t((((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()).Atoi());
					AliDebug(2,Form("Value at index %d = %d",nfixed+ielement+ientry*nValuesPerEntry,Int_t(value[ielement])));
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
				delete[] value;
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
				delete[] value;
			} 
			else if (type == "s"){
				TObjArray* value = new TObjArray();
				value->SetOwner(1);
				for (Int_t ielement=0; ielement<nelements; ielement++){
				  TObjString* strobj = (new TObjString(((TObjString*)tokens->At(nfixed+ielement+ientry*nValuesPerEntry))->String()));
					AliDebug(2,Form("Value at index %d = %s",nfixed+ielement+ientry*nValuesPerEntry,(strobj->String()).Data()));
					value->Add(strobj);
				}
				AliDCSArray* dcs = new AliDCSArray(nelements,value,timestamp);
				array->Add(dcs);
				delete value;
			}
			else{
				AliError(Form("Non-expected type %s",type.Data()));
				delete typeTokens;
				delete tokens;	
				delete tokenslines;
				return NULL;
			} 
		}
		delete typeTokens;
		delete tokens;
	}
	delete tokenslines;
	return array;
}







