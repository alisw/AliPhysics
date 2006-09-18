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
//  class AliCDBMetaData					   //
//  Set of data describing the object  				   //
//  but not used to identify the object 			   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include <TObjString.h>

ClassImp(AliCDBMetaData)

//_____________________________________________________________________________
AliCDBMetaData::AliCDBMetaData() :
TObject(),
fObjectClassName(""),
fResponsible(""),
fBeamPeriod(0),
fAliRootVersion(""),
fComment(""),
fProperties()	
{
// default constructor

	fProperties.SetOwner(1);
}

//_____________________________________________________________________________
AliCDBMetaData::AliCDBMetaData(const char *responsible, UInt_t beamPeriod,
				const char* alirootVersion, const char* comment) :
TObject(),
fObjectClassName(""),
fResponsible(responsible),
fBeamPeriod(beamPeriod),
fAliRootVersion(alirootVersion),
fComment(comment),
fProperties()	
{
// constructor

	fProperties.SetOwner(1);
}

//_____________________________________________________________________________
AliCDBMetaData::~AliCDBMetaData() {
// destructor

}

//_____________________________________________________________________________
void AliCDBMetaData::SetProperty(const char* property, TObject* object) {
// add something to the list of properties

	fProperties.Add(new TObjString(property), object);
}

//_____________________________________________________________________________
TObject* AliCDBMetaData::GetProperty(const char* property) const {
// get a property specified by its name (property)

	return fProperties.GetValue(property);
}

//_____________________________________________________________________________
Bool_t AliCDBMetaData::RemoveProperty(const char* property) {
// removes a property

	TObjString objStrProperty(property);
	TObjString* aKey = (TObjString*) fProperties.Remove(&objStrProperty);	

	if (aKey) {
		delete aKey;
		return kTRUE;
	} else {
		return kFALSE;
	}
}

//_____________________________________________________________________________
void AliCDBMetaData::PrintMetaData() {
// print the object's metaData

	TString message;
	if(fObjectClassName != "")
		message += Form("\tObject's class name:	%s\n", fObjectClassName.Data());
	if(fResponsible != "")
		message += Form("\tResponsible:		%s\n", fResponsible.Data());
	if(fBeamPeriod != 0)
		message += Form("\tBeam period:		%d\n", fBeamPeriod);
	if(fAliRootVersion != "")
		message += Form("\tAliRoot version:	%s\n", fAliRootVersion.Data());
	if(fComment != "")
		message += Form("\tComment:		%s\n", fComment.Data());
	if(fProperties.GetEntries() > 0){
		message += "\tProperties key names:";

		TIter iter(fProperties.GetTable());
		TPair* aPair;
		while ((aPair = (TPair*) iter.Next())) {
			message += Form("\t\t%s\n", ((TObjString* ) aPair->Key())->String().Data());
		}
	}
	message += '\n';
	AliInfo(Form("**** Object's MetaData parameters **** \n%s", message.Data()));
}
