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

//
//  This class is a container for the data queried from DAQ's logbook and logbook_shuttle tables.
//  It holds the run number, the start time and end time of the run,
//  and the array of the detectors' status (Unprocessed, Inactive, Failed, Done)

#include "AliShuttleLogbookEntry.h"
#include "AliLog.h"
#include "TTimeStamp.h"
#include <TString.h>
#include <TObjString.h>
#include <TMap.h>

ClassImp(AliShuttleLogbookEntry)

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry() :
TObject(),
fRun(-1),
fRunParameters(0)
{
  // default constructor

	const UInt_t nDet = AliShuttleInterface::NDetectors();
	memset(fDetectorStatus, kUnprocessed, nDet*sizeof(Status));
	fRunParameters.SetOwner(1);
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry(Int_t run, Status* status) :
TObject(),
fRun(run),
fRunParameters(0)
{
// default constructor

	const UInt_t nDet = AliShuttleInterface::NDetectors();
	memset(fDetectorStatus, kUnprocessed, nDet*sizeof(Status));
	if(status) SetDetectorStatus(status);
	fRunParameters.SetOwner(1);
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::~AliShuttleLogbookEntry() {
// destructor

}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry(const AliShuttleLogbookEntry &c) :
TObject(),
fRun(c.fRun),
fRunParameters(0)
{
  // copy constructor

	SetDetectorStatus(c.GetDetectorStatus());
	fRunParameters.SetOwner(1);
	TIter iter(c.fRunParameters.GetTable());
	TPair* aPair = 0;
	while((aPair = dynamic_cast<TPair*>(iter.Next()))){
		TObjString* aKey= dynamic_cast<TObjString*>(aPair->Key());
		TObjString* aValue= dynamic_cast<TObjString*>(aPair->Value());
		fRunParameters.Add(aKey->Clone(), aValue->Clone());
	}

}

//______________________________________________________________________________________________
AliShuttleLogbookEntry &AliShuttleLogbookEntry::operator=(const AliShuttleLogbookEntry &c)
{
  // assigment operator

	if (this != &c)
		((AliShuttleLogbookEntry &) c).Copy(*this);
	return *this;
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::Copy(TObject& c) const
{
  // copy function

	AliShuttleLogbookEntry& target = (AliShuttleLogbookEntry &) c;

	target.fRun = fRun;
	target.fRunParameters.SetOwner(1);
	TIter iter(fRunParameters.GetTable());
	TPair* aPair = 0;
	while((aPair = dynamic_cast<TPair*>(iter.Next()))){
		TObjString* aKey= dynamic_cast<TObjString*>(aPair->Key());
		TObjString* aValue= dynamic_cast<TObjString*>(aPair->Value());
		target.fRunParameters.Add(aKey->Clone(), aValue->Clone());
	}

	target.SetDetectorStatus(GetDetectorStatus());
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::Status AliShuttleLogbookEntry::GetDetectorStatus(const char* detName) const
{
// get detector status from detector code

	return GetDetectorStatus(AliShuttleInterface::GetDetPos(detName));
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::Status AliShuttleLogbookEntry::GetDetectorStatus(Int_t detPos) const
{
// get detector status from detector code

	if(detPos < 0 || detPos >= (Int_t) AliShuttleInterface::NDetectors()) {
		AliError(Form("Invalid parameter: %d", detPos));
		return kUnprocessed;
	}
	return fDetectorStatus[detPos];
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(const char* detName, Status status)
{
// set detector status from detector code

	Int_t detPos = AliShuttleInterface::GetDetPos(detName);
	if(detPos<0) return;
	SetDetectorStatus(detPos, status);
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(const char* detName, const char* statusName)
{
// set detector status from detector code

	Int_t detPos = AliShuttleInterface::GetDetPos(detName);
	if(detPos<0) return;
	SetDetectorStatus(detPos, statusName);
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(Status* status)
{
// set detector status from detector code

	for(UInt_t i=0; i < AliShuttleInterface::NDetectors(); i++){
		fDetectorStatus[i] = status[i];
	}
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(UInt_t detPos, Status status)
{
// set detector status from detector code

	if(detPos >= AliShuttleInterface::NDetectors()) {
		AliError(Form("Shuttle has only %d subdetectors!", AliShuttleInterface::NDetectors()));
		return;
	}
	fDetectorStatus[detPos] = status;
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(UInt_t detPos, const char* statusName)
{
// set detector status from detector code

	if(detPos >= AliShuttleInterface::NDetectors()) {
		AliError(Form("Shuttle has only %d subdetectors!", AliShuttleInterface::NDetectors()));
		return;
	}
	TString statusString(statusName);
	if(statusString.Contains("UNPROCESSED", TString::kIgnoreCase)){
		SetDetectorStatus(detPos, kUnprocessed);
	} else if (statusString.Contains("INACTIVE", TString::kIgnoreCase)) {
		SetDetectorStatus(detPos, kInactive);
	} else if (statusString.Contains("FAILED", TString::kIgnoreCase)) {
		SetDetectorStatus(detPos, kFailed);
	} else if (statusString.Contains("DONE", TString::kIgnoreCase)) {
		SetDetectorStatus(detPos, kDone);
	} else {
		AliError(Form("Invalid status name: %s", statusName));
	}
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::IsDone() const{
// return TRUE if all subdetectors are in status DONE, FAILED or INACTIVE

	for(UInt_t i=0; i < AliShuttleInterface::NDetectors(); i++){
		if(fDetectorStatus[i] == kUnprocessed) return kFALSE;
	}
	return kTRUE;
}

//______________________________________________________________________________________________
const char* AliShuttleLogbookEntry::GetDetectorStatusName(Status status)
{
  // returns a name (string) of the detector status

      switch (status){
	    case kUnprocessed: return "UNPROCESSED";
	    case kInactive: return "INACTIVE";
	    case kFailed: return "FAILED";
	    case kDone: return "DONE";
     }
     return 0;

}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::Print(Option_t* option) const
{
  // print current shuttle logbook entry

	TString message = "\n*** Run parameters ***\n";
	TTimeStamp startTimeStamp(GetStartTime());
	TTimeStamp endTimeStamp(GetEndTime());
	message += Form("\tRun \t\t%d\n", fRun);
	message += Form("\tStarted \t%s\n", startTimeStamp.AsString("s"));
	message += Form("\tFinished \t%s\n", endTimeStamp.AsString("s"));
	message += "\n*** Detector status ***\n";

	for(UInt_t i=0; i < AliShuttleInterface::NDetectors(); i++)
		message += Form("\t%2d - %s: %s\n", i, AliShuttleInterface::GetDetName(i),
					GetDetectorStatusName(fDetectorStatus[i]));

	AliInfo(Form("option: %s",option));
	TString optionStr(option);
	if(optionStr=="all"){
		message += "\nPrinting full list of run parameters\n";
		message += "\tParameter                      Value\n";
		TIter iter(fRunParameters.GetTable());
		TPair* aPair = 0;
		while((aPair = dynamic_cast<TPair*>(iter.Next()))){
			TObjString* aKey= dynamic_cast<TObjString*>(aPair->Key());
			TObjString* aValue= dynamic_cast<TObjString*>(aPair->Value());
			TString keyStr=aKey->GetName();
			if(keyStr != "log"){
				message += Form("\t%s ", aKey->GetName());
				if(keyStr.Length()<30) message.Append(' ', 30-keyStr.Length());
				message += Form("%s\n", aValue->GetName());
			} else {
				message += "\tlog                            ...\n";
			}
		}
	}

	AliInfo(Form("%s",message.Data()));
}
//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetRunParameter(const char* key, const char* value){
// set a run parameter (read from the DAQ logbook)

	TObjString* keyObj = new TObjString(key);
	if (fRunParameters.Contains(key)) {
		AliWarning(Form("Parameter %s already existing and it will be replaced.", key));
		delete fRunParameters.Remove(keyObj);

	}
	fRunParameters.Add(keyObj, new TObjString(value));
	AliDebug(2, Form("Number of parameters: %d", fRunParameters.GetEntries()));
}
//______________________________________________________________________________________________
const char* AliShuttleLogbookEntry::GetRunParameter(const char* key) const{
// get a run parameter

	TObjString* value = dynamic_cast<TObjString*> (fRunParameters.GetValue(key));
	if(!value) {
		AliError(Form("No such parameter: %s", key));
		return 0;
	}
	return value->GetName();
}
