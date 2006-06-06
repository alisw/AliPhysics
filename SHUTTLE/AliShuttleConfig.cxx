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

/*
$Log$
Revision 1.7  2006/05/12 09:07:16  colla
12/05/06
New configuration complete

Revision 1.2  2006/03/07 07:52:34  hristov
New version (B.Yordanov)

Revision 1.4  2005/11/19 14:20:31  byordano
logbook config added to AliShuttleConfig

Revision 1.3  2005/11/17 19:24:25  byordano
TList changed to TObjArray in AliShuttleConfig

Revision 1.2  2005/11/17 14:43:23  byordano
import to local CVS

Revision 1.1.1.1  2005/10/28 07:33:58  hristov
Initial import as subdirectory in AliRoot

Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.3  2005/08/30 09:13:02  byordano
some docs added

*/


//
// This class keeps the AliShuttle configuration.
// It reads the configuration for LDAP server.
// For every child entry in basedn which has schema type 'shuttleConfig'
// it creates a detector configuration. This configuration includes:
// DCS server host and port and the set of aliases for which data from
// will be retrieved (used by AliShuttle).
//


#include "AliShuttleConfig.h"

#include "AliLog.h"

#include <TSystem.h>
#include <TObjString.h>
#include <TLDAPResult.h>
#include <TLDAPEntry.h>
#include <TLDAPAttribute.h>

AliShuttleConfig::ConfigHolder::ConfigHolder(const TLDAPEntry* entry):
	fIsValid(kFALSE)
{
	TLDAPAttribute* anAttribute;
	
	anAttribute = entry->GetAttribute("det");
	if (!anAttribute) {
		AliError("Invalid configuration! Can't get detector name.");
		return;
	}
	fDetector = anAttribute->GetValue();
	if (!fDetector.Length()) {
		AliError("Detector name can't be an empty string!")
		return;
	}

	anAttribute = entry->GetAttribute("DCSHost");
	if (!anAttribute) {
		AliError("Invalid configuration! Can't get DCSHost.");
		return;
	}
	fDCSHost = anAttribute->GetValue();
	if (!fDCSHost.Length()) {
		AliError("Host can't be an empty string!")
		return;
	}

	anAttribute = entry->GetAttribute("DCSPort");
        if (!anAttribute) {
		AliError("Invalid configuration! Can't get DCSPort.");
		return;
        }
	TString portStr = anAttribute->GetValue();
	if (!portStr.Length()) {
		AliError("port can't be an empty string!")
		return;
	}
	fDCSPort = portStr.Atoi();

	anAttribute = entry->GetAttribute("DCSAlias");
        if (!anAttribute) {
		AliError("Invalid configuration! Can't get alias attribute.");
		return;
	}
	const char* anAlias;
	while ((anAlias	= anAttribute->GetValue())) {
		fDCSAliases.AddLast(new TObjString(anAlias));
	}

	anAttribute = entry->GetAttribute("DAQFileIDs");
        if (!anAttribute) {
		AliError("Invalid configuration! Can't get DAQFileIDs attribute.");
		return;
	}
	const char* aFileID;
	while ((aFileID	= anAttribute->GetValue())) {
		fDAQFileIDs.AddLast(new TObjString(aFileID));
	}
		
	fIsValid = kTRUE;
}

AliShuttleConfig::ConfigHolder::~ConfigHolder() {
	fDCSAliases.Delete();
}

ClassImp(AliShuttleConfig)

AliShuttleConfig::AliShuttleConfig(const char* host, Int_t port, 
	const char* binddn, const char* password, const char* basedn):
	fIsValid(kFALSE),
	fProcessAll(kFALSE)
{
	//
	// host: ldap server host
	// port: ldap server port
	// binddn: binddn used for ldap binding (simple bind is used!).
	// password: password for binddn
	// basedn: this is basedn whose childeren entries which have
	// (objectClass=shuttleConfig) will be used as detector configurations.
	//

	TLDAPServer aServer(host, port, binddn, password, 3);
	
	if (!aServer.IsConnected()) {
		AliError(Form("Can't connect to ldap server %s:%d", 
				host, port));
		return;
	}

	// reads configuration for the shuttle running on this machine
	
	fShuttleInstanceHost = gSystem->HostName();
	TString queryFilter = "(ShuttleHost=";
	queryFilter += fShuttleInstanceHost;
	queryFilter += ")";	
	
	TLDAPResult* aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL,
			queryFilter.Data());

	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}

	if (aResult->GetCount() == 0) {
		AliError(Form("No Shuttle instance for host = %s!",fShuttleInstanceHost.Data()));
		AliError(Form("All detectors will be processed."));
		fProcessAll=kTRUE;
	}

	if (aResult->GetCount() > 1) {
		AliError(Form("More than one Shuttle instance for host %s!",fShuttleInstanceHost.Data()));
		return;
	}
	
	TLDAPEntry* anEntry;
	TLDAPAttribute* anAttribute;

	if(!fProcessAll){
		anEntry = aResult->GetNext();
		anAttribute = anEntry->GetAttribute("detectors");
		const char *detName;
		while((detName = anAttribute->GetValue())){
			TObjString *objDet= new TObjString(detName);
			fProcessedDetectors.Add(objDet);
		}
	}

	aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL,
			"(objectClass=AliShuttleDetector)");
	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}

	while ((anEntry = aResult->GetNext())) {
		ConfigHolder* aHolder = new ConfigHolder(anEntry);
		delete anEntry;

		if (!aHolder->IsValid()) {
			AliError("This entry is going to be skipped!");
			delete aHolder;

			continue;
		}

		TObjString* detStr = new TObjString(aHolder->GetDetector());
		fDetectorMap.Add(detStr, aHolder);
		fDetectorList.AddLast(detStr);
	}	
	
	delete aResult;


	aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL,
			"(objectClass=AliShuttleGlobalConfig)");
	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}

	if (aResult->GetCount() == 0) {
		AliError("Can't find DAQ logbook configuration!");
		return;
	}

	if (aResult->GetCount() > 1) {
		AliError("More than one DAQ logbook configuration found!");
		return;
	}

	anEntry = aResult->GetNext();
	
	anAttribute = anEntry->GetAttribute("DAQLogbookHost");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookHost attribute!");
		return;
	}
	fDAQLogBookHost = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("DAQLogbookUser");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookUser attribute!");
		return;
	}
	fDAQLogBookUser = anAttribute->GetValue();
	
	anAttribute = anEntry->GetAttribute("DAQLogbookPassword");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookPassword attribute!");
		return;
	}
	fDAQLogBookPassword = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("DAQFileSystemHost");
	if (!anAttribute) {
		AliError("Can't find DAQFileSystemHost attribute!");
		return;
	}
	fDAQFSHost = anAttribute->GetValue();

	delete anEntry;
	delete aResult;
	
	fIsValid = kTRUE;
}

AliShuttleConfig::~AliShuttleConfig() {
	fDetectorMap.DeleteAll();
}

const TObjArray* AliShuttleConfig::GetDetectors() const {
	//
	// returns collection of TObjString which contains the name
	// of every detector which is in the configuration.
	//

	return &fDetectorList;
}

Bool_t AliShuttleConfig::HasDetector(const char* detector) const {
	//
	// checks for paricular detector in the configuration.
	//
	return fDetectorMap.GetValue(detector) != NULL;
}

const char* AliShuttleConfig::GetDCSHost(const char* detector) const {
	//
	// returns DCS server host used by particular detector
	//
	
	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
	if (!aHolder) {
		AliError(Form("There isn't configuration for detector: %s",
			detector));
		return NULL;
	}

	return aHolder->GetDCSHost();
}

Int_t AliShuttleConfig::GetDCSPort(const char* detector) const {
	//
        // returns DCS server port used by particular detector
        //


	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return 0;
        }

	return aHolder->GetDCSPort();
}

const TObjArray* AliShuttleConfig::GetDCSAliases(const char* detector) const {
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector
	//

	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetDCSAliases();
}

const TObjArray* AliShuttleConfig::GetDAQFileIDs(const char* detector) const {
	//
	// returns collection of TObjString which represents the set of DAQ file IDs
	// which used for data retrieval for particular detector
	//

	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetDAQFileIDs();
}

Bool_t AliShuttleConfig::HostProcessDetector(const char* detector) const {
	// return TRUE if detector is handled by host or if fProcessAll is TRUE

	if(fProcessAll) return kTRUE;
	TIter iter(&fProcessedDetectors);
	TObjString* detName;
	while((detName = (TObjString*) iter.Next())){
		if(detName->String() == detector) return kTRUE;
	}
	return kFALSE;
}

void AliShuttleConfig::Print(Option_t* /*option*/) const {
	
	TString result;
	result += '\n';

	result += "Shuttle running on host: ";
	result += fShuttleInstanceHost;
	result += '\n';
	result += "Detectors handled by this host: ";
	TIter it(&fProcessedDetectors);
	TObjString* aDet;
	while ((aDet = (TObjString*) it.Next())) {
		result += aDet->String();
		result += ' ';
	}
	result += '\n';
	if(fProcessAll) result += "ProcessAll is ON";

	result += '\n';
	result += '\n';

	result += "DAQ LogBook Host: ";
	result += fDAQLogBookHost;
	result += '\n';
	result += "DAQ LogBook User: ";
	result += fDAQLogBookUser;
	result += '\n';
	result += "DAQ LogBook Password: ";
	result.Append('*', fDAQLogBookPassword.Length());
	result += '\n';
	result += '\n';
	result += "DAQ File System Host: ";
	result += fDAQFSHost;
	result += '\n';

	TIter iter(fDetectorMap.GetTable());
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		ConfigHolder* aHolder = (ConfigHolder*) aPair->Value();
		result += '\n';
		result += " Detector: ";
		result += aHolder->GetDetector();
		result += '\n';	
		result += " DCS Host: ";
		result += aHolder->GetDCSHost();
		result += '\n';
		result += " DCS Port: ";
		result += aHolder->GetDCSPort();
		result += '\n';

		result += " DCS Aliases: ";
		const TObjArray* aliases = aHolder->GetDCSAliases(); 	
		TIter it(aliases);
		TObjString* anAlias;	
		while ((anAlias = (TObjString*) it.Next())) {
			result += anAlias->String();
			result += ' ';
		}	
		
		result += '\n';

		result += " DAQ File IDs: ";
		const TObjArray* fileIDs = aHolder->GetDAQFileIDs(); 	
		TIter it2(fileIDs);		
		TObjString* aFileID;	
		while ((aFileID = (TObjString*) it2.Next())) {
			result += aFileID->String();
			result += ' ';
		}	
		result += '\n';
	}

	AliInfo(result);
}
