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
Revision 1.9  2006/10/02 16:38:39  jgrosseo
update (alberto):
fixed memory leaks
storing of objects that failed to be stored to the grid before
interfacing of shuttle status table in daq system

Revision 1.8  2006/08/15 10:50:00  jgrosseo
effc++ corrections (alberto)

Revision 1.7  2006/07/20 09:54:40  jgrosseo
introducing status management: The processing per subdetector is divided into several steps,
after each step the status is stored on disk. If the system crashes in any of the steps the Shuttle
can keep track of the number of failures and skips further processing after a certain threshold is
exceeded. These thresholds can be configured in LDAP.

Revision 1.6  2006/07/19 10:09:55  jgrosseo
new configuration, accesst to DAQ FES (Alberto)

Revision 1.5  2006/07/10 13:01:41  jgrosseo
enhanced storing of last sucessfully processed run (alberto)

Revision 1.4  2006/06/12 09:11:16  jgrosseo
coding conventions (Alberto)

Revision 1.3  2006/06/06 14:26:40  jgrosseo
o) removed files that were moved to STEER
o) shuttle updated to follow the new interface (Alberto)

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
#include "AliShuttleInterface.h"

#include "AliLog.h"

#include <TSystem.h>
#include <TObjString.h>
#include <TLDAPResult.h>
#include <TLDAPEntry.h>
#include <TLDAPAttribute.h>


AliShuttleConfig::AliShuttleConfigHolder::AliShuttleConfigHolder(const TLDAPEntry* entry):
fDetector(""),
fDCSHost(""),
fDCSPort(0),
fIsValid(kFALSE),
fSkipDCSQuery(kFALSE)
{
// constructor of the shuttle configuration holder

	TLDAPAttribute* anAttribute;
	fDCSAliases = new TObjArray();
	fDCSAliases->SetOwner(1);

	anAttribute = entry->GetAttribute("det"); // MUST
	fDetector = anAttribute->GetValue();

	anAttribute = entry->GetAttribute("DCSHost"); // MAY
	if (!anAttribute) {
		AliWarning(
			Form("%s has not DCS host entry - Shuttle will skip DCS data query!",
				fDetector.Data()));
		fIsValid = kTRUE;
		fSkipDCSQuery = kTRUE;
		return;
	}

	fDCSHost = anAttribute->GetValue();

	anAttribute = entry->GetAttribute("DCSPort"); // MAY
        if (!anAttribute) {
		AliError(Form("Invalid configuration! %s has DCS Host but no port number!",
				fDetector.Data()));
		return;
        }
	TString portStr = anAttribute->GetValue();
	fDCSPort = portStr.Atoi();

	anAttribute = entry->GetAttribute("DCSAlias"); // MAY
        if (!anAttribute) {
		AliError(Form("Invalid configuration! %s has DCS host settings but no DCSAlias entries!",
				fDetector.Data()));
		return;
	}

	const char* anAlias;
	while ((anAlias	= anAttribute->GetValue())) {
		fDCSAliases->AddLast(new TObjString(anAlias));
	}

	fIsValid = kTRUE;


}

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleConfigHolder::~AliShuttleConfigHolder()
{
// destructor of the shuttle configuration holder

	delete fDCSAliases;
}

ClassImp(AliShuttleConfig)

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleConfig(const char* host, Int_t port,
	const char* binddn, const char* password, const char* basedn):
	fIsValid(kFALSE),
	fDAQlbHost(""), fDAQlbUser(""), fDAQlbPass(""),
	fMaxRetries(0), fPPTimeOut(0), fDetectorMap(), fDetectorList(),
	fShuttleInstanceHost(""), fProcessedDetectors(), fProcessAll(kFALSE)
{
	//
	// host: ldap server host
	// port: ldap server port
	// binddn: binddn used for ldap binding (simple bind is used!).
	// password: password for binddn
	// basedn: this is basedn whose childeren entries which have
	// (objectClass=shuttleConfig) will be used as detector configurations.
	//

	fDetectorMap.SetOwner();
	fDetectorList.SetOwner(0); //fDetectorList and fDetectorMap share the same object!
	fProcessedDetectors.SetOwner();

	TLDAPServer aServer(host, port, binddn, password, 3);

	if (!aServer.IsConnected()) {
		AliError(Form("Can't connect to ldap server %s:%d",
				host, port));
		return;
	}

	// reads configuration for the shuttle running on this machine

	fShuttleInstanceHost = gSystem->HostName();
	TString queryFilter = Form("(ShuttleHost=%s)", fShuttleInstanceHost.Data());

	TLDAPResult* aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL, queryFilter.Data());

	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}

	if (aResult->GetCount() == 0) {
		AliError(Form("No Shuttle instance for host = %s!",
					fShuttleInstanceHost.Data()));
		AliError(Form("All detectors will be processed."));
		fProcessAll=kTRUE;
	}

	if (aResult->GetCount() > 1) {
		AliError(Form("More than one Shuttle instance for host %s!",
					fShuttleInstanceHost.Data()));
		delete aResult;
		return;
	}

	TLDAPEntry* anEntry = 0;
	TLDAPAttribute* anAttribute = 0;

	if(!fProcessAll){
		anEntry = aResult->GetNext();
		anAttribute = anEntry->GetAttribute("detectors");
		const char *detName;
		while((detName = anAttribute->GetValue())){
			TObjString *objDet= new TObjString(detName);
			fProcessedDetectors.Add(objDet);
		}
	}

	delete anEntry; delete aResult;

	// Detector configuration (DCS Archive DB settings)

	aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL, "(objectClass=AliShuttleDetector)");
	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s", basedn));
		return;
	}


	while ((anEntry = aResult->GetNext())) {
		AliShuttleConfigHolder* aHolder = new AliShuttleConfigHolder(anEntry);
		delete anEntry;

		if (!aHolder->IsValid()) {
			AliError("Detector configuration error!");
			delete aHolder;
			delete aResult;
			return;
		}

		TObjString* detStr = new TObjString(aHolder->GetDetector());
		fDetectorMap.Add(detStr, aHolder);
		fDetectorList.AddLast(detStr);
	}

	delete aResult;

	// Global configuration (DAQ logbook)

	aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL,
			"(objectClass=AliShuttleGlobalConfig)");
	if (!aResult) {
		AliError(Form("Can't find configuration with base DN: %s",
				basedn));
		return;
	}

	if (aResult->GetCount() == 0) {
		AliError("Can't find DAQ logbook configuration!");
		delete aResult;
		return;
	}

	if (aResult->GetCount() > 1) {
		AliError("More than one DAQ logbook configuration found!");
		delete aResult;
		return;
	}

	anEntry = aResult->GetNext();

	anAttribute = anEntry->GetAttribute("DAQLogbookHost");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookHost attribute!");
		delete anEntry; delete aResult;
		return;
	}
	fDAQlbHost = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("DAQLogbookUser");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookUser attribute!");
		delete aResult; delete anEntry;
		return;
	}
	fDAQlbUser = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("DAQLogbookPassword");
	if (!anAttribute) {
		AliError("Can't find DAQLogbookPassword attribute!");
		delete aResult; delete anEntry;
		return;
	}
	fDAQlbPass = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("MaxRetries");
	if (!anAttribute) {
		AliError("Can't find MaxRetries attribute!");
		delete aResult; delete anEntry;
		return;
	}
	TString tmpStr = anAttribute->GetValue();
	fMaxRetries = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("PPTimeOut");
	if (!anAttribute) {
		AliError("Can't find PPTimeOut attribute!");
		delete aResult; delete anEntry;
		return;
	}
	tmpStr = anAttribute->GetValue();
	fPPTimeOut = tmpStr.Atoi();

	delete aResult; delete anEntry;

	// FES configuration (FES logbook and hosts)

	for(int iSys=0;iSys<3;iSys++){
		queryFilter = Form("(system=%s)", AliShuttleInterface::fkSystemNames[iSys]);
		aResult = aServer.Search(basedn, LDAP_SCOPE_ONELEVEL, queryFilter.Data());
		if (!aResult) {
			AliError(Form("Can't find configuration for system: %s",
					AliShuttleInterface::fkSystemNames[iSys]));
			return;
		}

		if (aResult->GetCount() != 1 ) {
			AliError("Error in FES configuration!");
			delete aResult;
			return;
		}

		anEntry = aResult->GetNext();

		anAttribute = anEntry->GetAttribute("LogbookHost");
		if (!anAttribute) {
			AliError(Form ("Can't find LogbookHost attribute for %s!!",
						AliShuttleInterface::fkSystemNames[iSys]));
			delete aResult; delete anEntry;
			return;
		}
		fFESlbHost[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("LogbookUser");
		if (!anAttribute) {
			AliError(Form ("Can't find LogbookUser attribute for %s!!",
						AliShuttleInterface::fkSystemNames[iSys]));
			delete aResult; delete anEntry;
			return;
		}
		fFESlbUser[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("LogbookPassword");
		if (!anAttribute) {
			AliError(Form ("Can't find LogbookPassword attribute for %s!!",
						AliShuttleInterface::fkSystemNames[iSys]));
			delete aResult; delete anEntry;
			return;
		}
		fFESlbPass[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("FSHost");
		if (!anAttribute) {
			AliError(Form ("Can't find FSHost attribute for %s!!",
						AliShuttleInterface::fkSystemNames[iSys]));
			delete aResult; delete anEntry;
			return;
		}
		fFESHost[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("FSUser");
		if (!anAttribute) {
			AliError(Form ("Can't find FSUser attribute for %s!!",
						AliShuttleInterface::fkSystemNames[iSys]));
			delete aResult; delete anEntry;
			return;
		}
		fFESUser[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("FSPassword");
		if (anAttribute) fFESPass[iSys] = anAttribute->GetValue();

		delete aResult; delete anEntry;
	}

	fIsValid = kTRUE;
}

//______________________________________________________________________________________________
AliShuttleConfig::~AliShuttleConfig()
{
// destructor

	fDetectorMap.DeleteAll();
	fDetectorList.Clear();
	fProcessedDetectors.Delete();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetDetectors() const
{
	//
	// returns collection of TObjString which contains the name
	// of every detector which is in the configuration.
	//

	return &fDetectorList;
}

//______________________________________________________________________________________________
Bool_t AliShuttleConfig::HasDetector(const char* detector) const
{
	//
	// checks for paricular detector in the configuration.
	//
	return fDetectorMap.GetValue(detector) != NULL;
}

//______________________________________________________________________________________________
const char* AliShuttleConfig::GetDCSHost(const char* detector) const
{
	//
	// returns DCS server host used by particular detector
	//
	
	AliShuttleConfigHolder* aHolder = (AliShuttleConfigHolder*) fDetectorMap.GetValue(detector);
	if (!aHolder) {
		AliError(Form("There isn't configuration for detector: %s",
			detector));
		return NULL;
	}

	return aHolder->GetDCSHost();
}

//______________________________________________________________________________________________
Int_t AliShuttleConfig::GetDCSPort(const char* detector) const
{
	//
        // returns DCS server port used by particular detector
        //


	AliShuttleConfigHolder* aHolder = (AliShuttleConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return 0;
        }

	return aHolder->GetDCSPort();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetDCSAliases(const char* detector) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector
	//

	AliShuttleConfigHolder* aHolder = (AliShuttleConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetDCSAliases();
}

//______________________________________________________________________________________________
Bool_t AliShuttleConfig::HostProcessDetector(const char* detector) const
{
	// return TRUE if detector is handled by host or if fProcessAll is TRUE

	if(fProcessAll) return kTRUE;
	TIter iter(&fProcessedDetectors);
	TObjString* detName;
	while((detName = (TObjString*) iter.Next())){
		if(detName->String() == detector) return kTRUE;
	}
	return kFALSE;
}

//______________________________________________________________________________________________
void AliShuttleConfig::Print(Option_t* /*option*/) const
{
// print configuration
	
	TString result;
	result += '\n';

	result += Form("\nShuttle running on %s \n\n", fShuttleInstanceHost.Data());

	if(fProcessAll) {
		result += Form("All detectors will be processed! \n\n");
	} else {
		result += "Detectors processed by this host: ";
		TIter it(&fProcessedDetectors);
		TObjString* aDet;
		while ((aDet = (TObjString*) it.Next())) {
			result += Form("%s ", aDet->String().Data());
		}
		result += "\n\n";
	}

	result += Form("PP time out = %d - Max total retries = %d\n\n", fPPTimeOut, fMaxRetries);

	result += Form("DAQ Logbook Configuration \n \tHost: %s - User: %s - ",
		fDAQlbHost.Data(), fDAQlbUser.Data());

	result += "Password: ";
	result.Append('*', fDAQlbPass.Length());
	result += "\n\n";

	for(int iSys=0;iSys<3;iSys++){
		result += Form("FES Configuration for %s system\n", AliShuttleInterface::fkSystemNames[iSys]);
		result += Form("\tLogbook host: \t%s - \tUser: %s\n",
						fFESlbHost[iSys].Data(), fFESlbUser[iSys].Data());
		// result += Form("Logbook Password:",fFESlbPass[iSys].Data());
		result += Form("\tFES host: \t%s - \tUser: %s\n\n", fFESHost[iSys].Data(), fFESUser[iSys].Data());
		// result += Form("FES Password:",fFESPass[iSys].Data());
	}

	TIter iter(fDetectorMap.GetTable());
	TPair* aPair;
	while ((aPair = (TPair*) iter.Next())) {
		AliShuttleConfigHolder* aHolder = (AliShuttleConfigHolder*) aPair->Value();
		if(aHolder->SkipDCSQuery()) continue;
		result += Form("DCS archive DB configuration for *** %s *** \n", aHolder->GetDetector());
		result += Form("\tAmanda server: %s:%d \n", aHolder->GetDCSHost(), aHolder->GetDCSPort());

		result += "\tDCS Aliases: ";
		const TObjArray* aliases = aHolder->GetDCSAliases();
		TIter it(aliases);
		TObjString* anAlias;
		while ((anAlias = (TObjString*) it.Next())) {
			result += Form("%s ", anAlias->String().Data());
		}

		result += "\n\n";

	}

	if(!fIsValid) result += "\n\n********** !!!!! Configuration is INVALID !!!!! **********\n";

	AliInfo(result);
}
