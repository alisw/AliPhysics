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

#include <TObjString.h>
#include <TLDAPResult.h>
#include <TLDAPEntry.h>
#include <TLDAPAttribute.h>

AliShuttleConfig::ConfigHolder::ConfigHolder(const TLDAPEntry* entry):
	fIsValid(kFALSE)
{
	TLDAPAttribute* anAttribute;
	
	anAttribute = entry->GetAttribute("dt");
	if (!anAttribute) {
		AliError("Invalid configuration! Can't get detector name.");
		return;
	}
	fDetector = anAttribute->GetValue();
	if (!fDetector.Length()) {
		AliError("Detector name can't be an empty string!")
		return;
	}

	anAttribute = entry->GetAttribute("ipHost");
	if (!anAttribute) {
		AliError("Invalid configuration! Can't get ipHost.");
		return;
	}
	fHost = anAttribute->GetValue();
	if (!fHost.Length()) {
		AliError("Host can't be an empty string!")
		return;
	}

	anAttribute = entry->GetAttribute("ipServicePort");
        if (!anAttribute) {
		AliError("Invalid configuration! Can't get ipServicePort.");
		return;
        }
	TString portStr = anAttribute->GetValue();
	if (!portStr.Length()) {
		AliError("ipServicePort can't be an empty string!")
		return;
	}
	fPort = portStr.Atoi();

	anAttribute = entry->GetAttribute("alias");
        if (!anAttribute) {
		AliError("Invalid configuration! Can't get alias attribute.");
		return;
	}
	const char* anAlias;
	while ((anAlias	= anAttribute->GetValue())) {
		fAliases.Add(new TObjString(anAlias));
	}
		
	fIsValid = kTRUE;
}

AliShuttleConfig::ConfigHolder::~ConfigHolder() {
	fAliases.Delete();
}

ClassImp(AliShuttleConfig)

AliShuttleConfig::AliShuttleConfig(const char* host, Int_t port, 
	const char* binddn, const char* password, const char* basedn):
	fIsValid(kFALSE)
{
	//
	// host: ldap server host
	// port: ldap server port
	// binddn: binddn used for ldap binding (simple bind is used!).
	// password: password for binddn
	// basedn: this is basedn whose childeren entries which have 
	// (objectClass=shuttleConfig) will be used as detector configurations.
	//

	TLDAPServer aServer(host, port, binddn, password);
	
	if (!aServer.IsConnected()) {
		AliError(Form("Can't connect to ldap server %s:%d", 
				host, port));
		return;
	}

	TLDAPResult* aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL,
			"(objectClass=shuttleConfig)");
	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}
	
	TLDAPEntry* anEntry;
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
		fDetectorList.Add(detStr);
	}	
	
	delete aResult;

	fIsValid = kTRUE;
}

AliShuttleConfig::~AliShuttleConfig() {
	fDetectorMap.DeleteAll();
}

const TList* AliShuttleConfig::GetDetectors() const {
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

const char* AliShuttleConfig::GetHost(const char* detector) const {
	//
	// returns DCS server host used by particular detector
	//
	
	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
	if (!aHolder) {
		AliError(Form("There isn't configuration for detector: %s",
			detector));
		return NULL;
	}

	return aHolder->GetHost();
}

Int_t AliShuttleConfig::GetPort(const char* detector) const {
	//
        // returns DCS server port used by particular detector
        //


	ConfigHolder* aHolder = (ConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return 0;
        }

	return aHolder->GetPort();
}

const TList* AliShuttleConfig::GetAliases(const char* detector) const {
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

	return aHolder->GetAliases();
}

void AliShuttleConfig::Print(Option_t* /*option*/) const {
	
	TString result;
	result += '\n';
	
	TIter iter(fDetectorMap.GetTable());
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		ConfigHolder* aHolder = (ConfigHolder*) aPair->Value();
		result += '\n';
		result += " Detector: ";
		result += aHolder->GetDetector();
		result += '\n';	
		result += " Host: ";
		result += aHolder->GetHost();
		result += '\n';
		result += " Port: ";
		result += aHolder->GetPort();
		result += '\n';

		result += " Aliases: ";
		const TList* aliases = aHolder->GetAliases(); 	
		TIter it(aliases);		
		TObjString* anAlias;	
		while ((anAlias = (TObjString*) it.Next())) {
			result += anAlias->String();
			result += ' ';
		}	
		
		result += '\n';
	}

	AliInfo(result);
}
