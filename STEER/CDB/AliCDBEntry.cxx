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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBEntry						   //
//  container for an object, it identity (AliCDBId)  		   //
//  and its metaData (AliCDBMetaData) 				   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBEntry.h"
#include "AliLog.h"

ClassImp(AliCDBEntry)

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry():
fObject(NULL),
fId(),
fMetaData(NULL), 
fIsOwner(kFALSE){
// default constructor

}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBId& id, 
			AliCDBMetaData* metaData, Bool_t owner):
fObject(object), 
fId(id), 
fMetaData(metaData), 
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBPath& path, 
			const AliCDBRunRange& runRange,
                        AliCDBMetaData* metaData,Bool_t owner):
fObject(object), 
fId(path, runRange, -1, -1), 
fMetaData(metaData),
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBPath& path, 
			const AliCDBRunRange& runRange,
			Int_t version, AliCDBMetaData* metaData, Bool_t owner):
fObject(object), 
fId(path, runRange, version, -1), 
fMetaData(metaData),
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBPath& path, 
			const AliCDBRunRange& runRange,
			Int_t version, Int_t subVersion, 
			AliCDBMetaData* metaData, Bool_t owner):
fObject(object),
fId(path, runRange, version, subVersion), 
fMetaData(metaData), 
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}


//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBPath& path, 
			Int_t firstRun, Int_t lastRun, 
			AliCDBMetaData* metaData, Bool_t owner):
fObject(object),
fId(path, firstRun, lastRun, -1, -1), 
fMetaData(metaData), 
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry(TObject* object, const AliCDBPath& path, 
			Int_t firstRun, Int_t lastRun,
			Int_t version, AliCDBMetaData* metaData,
			Bool_t owner):
fObject(object),
fId(path, firstRun, lastRun, version, -1),
fMetaData(metaData), 
fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::AliCDBEntry( TObject* object, const AliCDBPath& path, 
			Int_t firstRun, Int_t lastRun,
			Int_t version, Int_t subVersion,
			AliCDBMetaData* metaData, Bool_t owner):
fObject(object),
fId(path, firstRun, lastRun, version, subVersion),
fMetaData(metaData), fIsOwner(owner){
// constructor
    fMetaData->SetObjectClassName(fObject->ClassName());
}

//_____________________________________________________________________________
AliCDBEntry::~AliCDBEntry() {
// destructor

	if (fIsOwner) {
		if (fObject) {
			delete fObject;
		}

		if (fMetaData) {
			delete fMetaData;
		}
	}
}

//_____________________________________________________________________________
void AliCDBEntry::PrintId() const {
 
	AliInfo(Form("%s",fId.ToString().Data()));

}
