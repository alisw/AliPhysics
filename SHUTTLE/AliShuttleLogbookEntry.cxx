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

// TODO test only!
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TObjArray.h>

ClassImp(AliShuttleLogbookEntry)

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry() :
TObject(),
fRun(-1),
fStartTime(0),
fEndTime(0),
//fDetectorStatus(0),
fServer(0)
{
  // default constructor

	const UInt_t nDet = AliShuttle::NDetectors();
//	fDetectorStatus = new Status[nDet];
	memset(fDetectorStatus, kUnprocessed, nDet*sizeof(Status));
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry(Int_t run, UInt_t startTime, UInt_t endTime, Status* status) :
TObject(),
fRun(run),
fStartTime(startTime),
fEndTime(endTime),
//fDetectorStatus(0),
fServer(0)
{
  // default constructor

	const UInt_t nDet = AliShuttle::NDetectors();
//	fDetectorStatus = new Status[nDet];
	memset(fDetectorStatus, kUnprocessed, nDet*sizeof(Status));
	if(status) SetDetectorStatus(status);
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::~AliShuttleLogbookEntry() {
// destructor

	if(fServer){
		if(fServer->IsConnected()) fServer->Close();
		delete fServer;
	}
//	if(fDetectorStatus) delete[] fDetectorStatus; fDetectorStatus=0;
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::AliShuttleLogbookEntry(const AliShuttleLogbookEntry &c) :
TObject(),
fRun(c.fRun),
fStartTime(c.fStartTime),
fEndTime(c.fEndTime),
fServer(0)
{
  // copy constructor

  SetDetectorStatus(c.GetDetectorStatus());
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
	target.fStartTime = fStartTime;
	target.fEndTime = fEndTime;

	target.SetDetectorStatus(GetDetectorStatus());
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::Status AliShuttleLogbookEntry::GetDetectorStatus(const char* detCode) const
{
// get detector status from detector code

	return GetDetectorStatus(AliShuttle::GetDetPos(detCode));
}

//______________________________________________________________________________________________
AliShuttleLogbookEntry::Status AliShuttleLogbookEntry::GetDetectorStatus(Int_t detPos) const
{
// get detector status from detector code

	if(detPos < 0 || detPos >= (Int_t) AliShuttle::NDetectors()) {
		AliError(Form("Invalid parameter: %d", detPos));
		return kUnprocessed;
	}
	return fDetectorStatus[detPos];
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(const char* detCode, Status status)
{
// set detector status from detector code

	Int_t detPos = AliShuttle::GetDetPos(detCode);
	if(detPos<0) return;
	SetDetectorStatus(detPos, status);
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(Status* status)
{
// set detector status from detector code

	for(UInt_t i=0; i < AliShuttle::NDetectors(); i++){
		fDetectorStatus[i] = status[i];
	}
}

//______________________________________________________________________________________________
void AliShuttleLogbookEntry::SetDetectorStatus(UInt_t detPos, Status status)
{
// set detector status from detector code

	if(detPos >= AliShuttle::NDetectors()) {
		AliError(Form("Shuttle has only %d subdetectors!", AliShuttle::NDetectors()));
		return;
	}
	fDetectorStatus[detPos] = status;
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::IsDone() const{
// return TRUE if all subdetectors are in status DONE, FAILED or INACTIVE

	for(UInt_t i=0; i < AliShuttle::NDetectors(); i++){
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
void AliShuttleLogbookEntry::Print(Option_t* /*option*/) const
{
  // print current shuttle logbook entry

	TString message = "\n*** Run parameters ***\n";
	TTimeStamp startTimeStamp(fStartTime);
	TTimeStamp endTimeStamp(fEndTime);
	message += Form("\tRun \t\t%d\n", fRun);
	message += Form("\tStarted \t%s\n", startTimeStamp.AsString("s"));
	message += Form("\tFinished \t%s\n", endTimeStamp.AsString("s"));
	message += "\n*** Detector status ***\n";

	for(UInt_t i=0; i < AliShuttle::NDetectors(); i++)
		message += Form("\t%2d - %s: %s\n", i, AliShuttle::GetDetCode(i),
					GetDetectorStatusName(fDetectorStatus[i]));

	AliInfo(Form("%s",message.Data()));
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::Connect(){
// Connect to MySQL Server of the DAQ logbook

	// check connection: if already connected return
	if(fServer && fServer->IsConnected()) return kTRUE;

	fServer = TSQLServer::Connect("mysql://pcald30.cern.ch","offline","alice");

	if (!fServer || !fServer->IsConnected()) {
		AliError("Can't establish connection to DAQ log book DB!");
		if(fServer) delete fServer;
		return kFALSE;
	}

	// Get table
	TSQLResult* aResult=0;
	aResult = fServer->GetTables("REFSYSLOG");
	delete aResult;
	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::QueryShuttleLogbook(Int_t runNumber)
{
// Query DAQ's Shuttle logbook and fills detector status array

	Int_t run;
	if(runNumber < 0) {
		run = GetRun();
	} else{
		run = runNumber;
	}

	// check connection, in case connect
	if(!Connect()) return kFALSE;

	TString sqlQuery;
	sqlQuery = Form("select * from logbook_shuttle where run = %d", run);

	TSQLResult* aResult = fServer->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	// TODO Check field count!
	if (aResult->GetFieldCount() != 24) {
		AliError("Invalid SQL result field number!");
		delete aResult;
		return kFALSE;
	}

	TSQLRow* aRow;
	while ((aRow = aResult->Next())) {
		TString runString(aRow->GetField(0), aRow->GetFieldLength(0));
		Int_t run = runString.Atoi();

		Status detStatus[24];

		// loop on detectors
		for(UInt_t ii = 0; ii < 24; ii++){
			TString detCode(aResult->GetFieldName(ii));
			Int_t detPos = AliShuttle::GetDetPos(detCode.Data());
			if(detPos < 0) continue;
			TString statusString(aRow->GetField(ii), aRow->GetFieldLength(ii));
			if(statusString == "UNPROCESSED"){
				detStatus[detPos] = AliShuttleLogbookEntry::kUnprocessed;
			} else if (statusString == "INACTIVE") {
				detStatus[detPos] = AliShuttleLogbookEntry::kInactive;
			} else if (statusString == "FAILED") {
				detStatus[detPos] = AliShuttleLogbookEntry::kFailed;
			} else if (statusString == "DONE") {
				detStatus[detPos] = AliShuttleLogbookEntry::kDone;
			}
		}

		SetRun(run);
		SetDetectorStatus(detStatus);
		delete aRow;
	}

	Print("");

	delete aResult;
	return kTRUE;
}
//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::UpdateShuttleLogbook()
{
  // Update Shuttle logbook table - TEST ONLY, USE WITH CARE!


	if(!Connect()) return kFALSE;

	TString sqlQuery("update logbook_shuttle set ");

	for(UInt_t i=0; i < AliShuttle::NDetectors(); i++){
		sqlQuery += Form("%s=\"%s\"", AliShuttle::GetDetCode(i), GetDetectorStatusName(fDetectorStatus[i]));
		if(i < AliShuttle::NDetectors()-1) sqlQuery += ", ";
	}

	sqlQuery += Form(" where run=%d;",GetRun());

	AliInfo(Form("sqlQuery: %s", sqlQuery.Data()));

	TSQLResult* aResult;
	aResult = fServer->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	delete aResult;

	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::UpdateShuttleLogbook(const char* detCode, Status status)
{
  // Update Shuttle logbook table - TEST ONLY, USE WITH CARE!


	if(AliShuttle::GetDetPos(detCode) < 0) return kFALSE;
	SetDetectorStatus(detCode, status);
	if(!Connect()) return kFALSE;

	TString sqlQuery("update logbook_shuttle set ");


	sqlQuery += Form("%s=\"%s\" ", detCode, GetDetectorStatusName(status));

	sqlQuery += Form("where run=%d;",GetRun());

	AliInfo(Form("sqlQuery: %s", sqlQuery.Data()));

	TSQLResult* aResult;
	aResult = fServer->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	delete aResult;

	return kTRUE;
}

//______________________________________________________________________________________________
Bool_t AliShuttleLogbookEntry::InsertNewRun(Int_t runNumber)
{
  // Update Shuttle logbook table - TEST ONLY, USE WITH CARE!

	if(runNumber<=0 && GetRun()<=0) return kFALSE;
	if(runNumber>0) SetRun(runNumber);
	if(!Connect()) return kFALSE;

	TString sqlQuery = Form("insert into logbook_shuttle (run) values (%d);", GetRun());

	AliInfo(Form("sqlQuery: %s", sqlQuery.Data()));

	TSQLResult* aResult;
	aResult = fServer->Query(sqlQuery);
	if (!aResult) {
		AliError(Form("Can't execute query <%s>!", sqlQuery.Data()));
		return kFALSE;
	}

	delete aResult;

	UpdateShuttleLogbook();

	return kTRUE;
}
