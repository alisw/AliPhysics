// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no         *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTGlobalPreprocessor.cxx
//  @author Matthias Richter
//  @date   2010-08-20
//  @brief  HLT Preprocessor plugin for global HLT
// 

#include "AliHLTGlobalPreprocessor.h"
#include "AliPreprocessor.h"
#include "AliHLTMisc.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliLog.h"
#include <cassert>
#include <cerrno>

#include <TObjArray.h>
#include <TStreamerInfo.h>
//#include <AliDCSValue.h>
//#include <TTimeStamp.h>


ClassImp(AliHLTGlobalPreprocessor)

AliHLTGlobalPreprocessor::AliHLTGlobalPreprocessor()
  : AliHLTModulePreprocessor()
{
  // constructor
  // HLT preprocessor for global HLT objects
  //
  // Produced OCDB objects:
  // - HLT/Calib/Streamerinfo
  //   The streamer info object is produced by the ROOTSchemaEvolutionComponent
  //   See ProcessStreamerInfo() for details.
}

const char* AliHLTGlobalPreprocessor::fgkStreamerInfoAlias="StreamerInfo";
const char* AliHLTGlobalPreprocessor::fgkStreamerInfoName="StreamerInfo";
const char* AliHLTGlobalPreprocessor::fgkStreamerInfoType="Calib";

AliHLTGlobalPreprocessor::~AliHLTGlobalPreprocessor()
{
  // destructor
}


void AliHLTGlobalPreprocessor::Initialize(Int_t /*run*/, UInt_t /*startTime*/, 
					  UInt_t /*endTime*/)
{
  // initializes AliHLTGlobalPreprocessor

}

UInt_t AliHLTGlobalPreprocessor::Process(TMap* /*dcsAliasMap*/)
{
  Int_t returnValue = ProcessStreamerInfo();
  if (returnValue < 0) {
	AliInfo(Form("Processing for %s failed with return code %d", fgkStreamerInfoAlias, returnValue));
  }
  return 0; // return success
}

Int_t AliHLTGlobalPreprocessor::GetModuleNumber()
{
  // get module number of this preprocessor, corresponds to the position
  // in the detector bit field, or 0 if no corresponding detector existing
  Int_t modulenumber = 0;
  //modulenumber = AliHLTModulePreprocessor::DetectorBitMask("GRP");
  return modulenumber;
}

Int_t AliHLTGlobalPreprocessor::ProcessStreamerInfo() {
	// get file sources
	TList* list = GetFileSources(AliPreprocessor::kHLT, fgkStreamerInfoAlias);
	if ((!list) || (list->GetEntries() == 0)) {
		AliInfo(Form("No sources for %s found",fgkStreamerInfoAlias));
		return -1; // no sources
	}
	bool bStore = false;
    TObjArray* clone = NULL;
	// get existing object or create new one
    AliCDBEntry* entry = GetFromOCDB(fgkStreamerInfoType, fgkStreamerInfoName);
    if (entry && entry->GetObject()) {
    	TObject* cloneObj = entry->GetObject()->Clone();
    	if (cloneObj) clone = dynamic_cast<TObjArray*>(cloneObj);
    } else {
    	clone = new TObjArray();
    	bStore = true;
    }
    if (!clone) {
    	AliError(Form("Could not clone %s, %s", fgkStreamerInfoType, fgkStreamerInfoName));
    	return -2; // no clone
    }
	// loop over all sources
	TObjLink *lnk = list->FirstLink();
	while (lnk) {
		TObject* obj = lnk->GetObject();
		TObjString* objStr = dynamic_cast<TObjString*>(obj);
		if (!objStr) {
			AliError(Form("GetFileSources returned TList with no TObjString entry?! %s", obj->ClassName()));
			// continue with next list entry
			lnk = lnk->Next();
			continue;
		}
	  	TString fileName = GetFile(AliPreprocessor::kHLT,fgkStreamerInfoAlias ,objStr->GetString().Data());
	  	if (fileName.Length() == 0) {
	    		AliError(Form("Could not get %d-%s-%s", AliPreprocessor::kHLT,fgkStreamerInfoAlias ,objStr->GetString().Data()));
	  	}
	  	// time to process the file...
	  	TFile* f = new TFile(fileName.Data(), "READ");
	  	if (!f || !f->IsOpen()) {
	  		AliError(Form("Could not open %s", objStr->GetString().Data()));
	  		return -3;
	  	}
	  	// loop over objects and create new TObjArrary to feed into merger...
	  	TObjArray* streamerinfos = new TObjArray(100);
	  	TList* keys = f->GetListOfKeys();
	  	TObjLink *lnkFile = keys->FirstLink();
	  	while (lnkFile) {
	  		// processing
	  		TObject* streamerobj = f->Get(lnkFile->GetObject()->GetName());
	  		TStreamerInfo* streamer=dynamic_cast<TStreamerInfo*>(streamerobj);
	  		if (!streamer) {
	  			AliError(Form("StreamerInfo object has wrong class type %s, expecting TStreamerInfo", streamerobj->ClassName()));
	  		} else {
	  			streamerinfos->Add(streamer);
	  		}
	  		lnkFile = lnkFile->Next();
	  	}
	  	if (streamerinfos->GetEntriesFast()!=0) {
	  		bStore |= AliHLTMisc::Instance().MergeStreamerInfo(clone, streamerinfos, 1)>0;
	  	}
	  	delete streamerinfos;
	  	f->Close();
	  	delete f;
	  	lnk = lnk->Next();
	}
	// store if necessary
	if (bStore) {
	    AliCDBMetaData* metaData=entry?entry->GetMetaData():NULL;
	    AliCDBMetaData* newMetaData=NULL;
	    if (!metaData) {
	    	newMetaData=new AliCDBMetaData;
	    	if (newMetaData) {
	    		metaData=newMetaData;
	    		metaData->SetBeamPeriod(0);
	    		metaData->SetResponsible("ALICE HLT alice-hlt-core@cern.ch");
	    		metaData->SetComment("Streamer info for HLTOUT payload");
	    		//metaData->SetAliRootVersion(ALIROOT_VERSION);
	    	} else {
	    		return -ENOMEM;
	    	}
	    }
	     Store(fgkStreamerInfoType, fgkStreamerInfoName, clone, metaData, 0, kTRUE);
	     if (newMetaData) {
	     	delete newMetaData;
	        newMetaData=NULL;
	        metaData=NULL;
	     }
	} else {
		if (entry) {
	      AliInfo(Form("StreamerInfo object in OCDB is already up-to-date, skipping new object"));
	      //entry->PrintId();
		}
	}
    delete clone;
	return 0;
}


