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
Revision 1.27  2007/12/17 03:23:32  jgrosseo
several bugfixes
added "empty preprocessor" as placeholder for Acorde in FDR

Revision 1.26  2007/12/07 19:14:36  acolla
in AliShuttleTrigger:

Added automatic collection of new runs on a regular time basis (settable from the configuration)

in AliShuttleConfig: new members

- triggerWait: time to wait for DIM trigger (s) before starting automatic collection of new runs
- mode: run mode (test, prod) -> used to build log folder (logs or logs_PROD)

in AliShuttle:

- logs now stored in logs/#RUN/DET_#RUN.log

Revision 1.25  2007/11/26 16:58:37  acolla
Monalisa configuration added: host and table name

Revision 1.24  2007/10/24 10:44:08  acolla

debug AliInfo removed

Revision 1.23  2007/09/28 15:27:40  acolla

AliDCSClient "multiSplit" option added in the DCS configuration
in AliDCSMessage: variable MAX_BODY_SIZE set to 500000

Revision 1.22  2007/09/27 16:53:13  acolla
Detectors can have more than one AMANDA server. SHUTTLE queries the servers sequentially,
merges the dcs aliases/DPs in one TMap and sends it to the preprocessor.

Revision 1.21  2007/04/27 07:06:48  jgrosseo
GetFileSources returns empty list in case of no files, but successful query
No mails sent in testmode

Revision 1.20  2007/04/04 10:33:36  jgrosseo
1) Storing of files to the Grid is now done _after_ your preprocessors succeeded. This is transparent, which means that you can still use the same functions (Store, StoreReferenceData) to store files to the Grid. However, the Shuttle first stores them locally and transfers them after the preprocessor finished. The return code of these two functions has changed from UInt_t to Bool_t which gives you the success of the storing.
In case of an error with the Grid, the Shuttle will retry the storing later, the preprocessor does not need to be run again.

2) The meaning of the return code of the preprocessor has changed. 0 is now success and any other value means failure. This value is stored in the log and you can use it to keep details about the error condition.

3) New function StoreReferenceFile to _directly_ store a file (without opening it) to the reference storage.

4) The memory usage of the preprocessor is monitored. If it exceeds 2 GB it is terminated.

5) New function AliPreprocessor::ProcessDCS(). If you do not need to have DCS data in all cases, you can skip the processing by implemting this function and returning kFALSE under certain conditions. E.g. if there is a certain run type.
If you always need DCS data (like before), you do not need to implement it.

6) The run type has been added to the monitoring page

Revision 1.19  2007/02/28 10:41:56  acolla
Run type field added in SHUTTLE framework. Run type is read from "run type" logbook and retrieved by
AliPreprocessor::GetRunType() function.
Added some ldap definition files.

Revision 1.18  2007/01/23 19:20:03  acolla
Removed old ldif files, added TOF, MCH ldif files. Added some options in
AliShuttleConfig::Print. Added in Ali Shuttle: SetShuttleTempDir and
SetShuttleLogDir

Revision 1.17  2007/01/18 11:17:47  jgrosseo
changing spaces to tabs ;-)

Revision 1.16  2007/01/18 11:10:35  jgrosseo
adding the possibility of defining DCS alias and data points with patterns
first pattern introduced: [N..M] to add all names between the two digits, this works also recursively.

Revision 1.15  2007/01/15 18:27:11  acolla
implementation of sending mail to subdetector expert in case the preprocessor fails.
shuttle.schema updated with expert's email entry

Revision 1.13  2006/12/07 08:51:26  jgrosseo
update (alberto):
table, db names in ldap configuration
added GRP preprocessor
DCS data can also be retrieved by data point

Revision 1.12  2006/11/16 16:16:48  jgrosseo
introducing strict run ordering flag
removed giving preprocessor name to preprocessor, they have to know their name themselves ;-)

Revision 1.11  2006/11/06 14:23:04  jgrosseo
major update (Alberto)
o) reading of run parameters from the logbook
o) online offline naming conversion
o) standalone DCSclient package

Revision 1.10  2006/10/20 15:22:59  jgrosseo
o) Adding time out to the execution of the preprocessors: The Shuttle forks and the parent process monitors the child
o) Merging Collect, CollectAll, CollectNew function
o) Removing implementation of empty copy constructors (declaration still there!)

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

#include <Riostream.h>
#include "AliShuttleConfig.h"
#include "AliShuttleInterface.h"

#include "AliLog.h"

#include <TSystem.h>
#include <TObjString.h>
#include <TLDAPResult.h>
#include <TLDAPEntry.h>
#include <TLDAPAttribute.h>
#include <TKey.h>


AliShuttleConfig::AliShuttleDCSConfigHolder::AliShuttleDCSConfigHolder(const TLDAPEntry* entry):
fDCSHost(""),
fDCSPort(0),
fMultiSplit(100),
fDCSAliases(0),
fDCSDataPoints(0),
fDCSAliasesComp(0),
fDCSDataPointsComp(0),
fIsValid(kFALSE)
{
// constructor of the shuttle DCS configuration holder

	TLDAPAttribute* anAttribute;
	fDCSAliases = new TObjArray();
	fDCSAliases->SetOwner(1);
	fDCSDataPoints = new TObjArray();
	fDCSDataPoints->SetOwner(1);
	fDCSAliasesComp = new TObjArray();
	fDCSAliasesComp->SetOwner(1);
	fDCSDataPointsComp = new TObjArray();
	fDCSDataPointsComp->SetOwner(1);

	
	anAttribute = entry->GetAttribute("dcsHost"); 
	if (!anAttribute)
	{
		AliError("Unexpected: no DCS host!");
		return;
	}

	fDCSHost = anAttribute->GetValue();

	anAttribute = entry->GetAttribute("dcsPort");
        if (!anAttribute)
	{
		AliError("Unexpected: no DCS port!");
		return;
        }
	TString portStr = anAttribute->GetValue();
	fDCSPort = portStr.Atoi();

	anAttribute = entry->GetAttribute("multiSplit"); // MAY
        if (anAttribute)
	{
		TString multiSplitStr = anAttribute->GetValue();
		fMultiSplit = multiSplitStr.Atoi();
		if(fMultiSplit == 0) {
			AliError("MultiSplit must be a positive integer!");
			return;
		}
		
	}
	
	anAttribute = entry->GetAttribute("dcsAlias"); // MAY
        if (anAttribute)
	{
		const char* anAlias;
		while ((anAlias	= anAttribute->GetValue()))
		{
			fDCSAliasesComp->AddLast(new TObjString(anAlias));
			ExpandAndAdd(fDCSAliases, anAlias);
		}
	}

	anAttribute = entry->GetAttribute("dcsDP"); // MAY
        if (anAttribute)
	{
		const char* aDataPoint;
		while ((aDataPoint = anAttribute->GetValue()))
		{
		fDCSDataPointsComp->AddLast(new TObjString(aDataPoint));
		ExpandAndAdd(fDCSDataPoints, aDataPoint);
		}
	}
	
	fIsValid = kTRUE;
}
//______________________________________________________________________________________________
void AliShuttleConfig::AliShuttleDCSConfigHolder::ExpandAndAdd(TObjArray* target, const char* entry)
{
	//
	// adds <entry> to <target> applying expanding of the name
	// [N..M] creates M-N+1 names with the corresponding digits
	//

	TString entryStr(entry);

	Int_t begin = entryStr.Index("[");
	Int_t end = entryStr.Index("]");
	if (begin != -1 && end != -1 && end > begin)
	{
		TString before(entryStr(0, begin));
		TString after(entryStr(end+1, entryStr.Length()));

		AliDebug(2, Form("Found [] pattern. Splitted input string %s %s", before.Data(), after.Data()));

		Int_t dotdot = entryStr.Index("..");

		TString nStr(entryStr(begin+1, dotdot-begin-1));
		TString mStr(entryStr(dotdot+2, end-dotdot-2));

		AliDebug(2, Form("Found [N..M] pattern. %s %s", nStr.Data(), mStr.Data()));

		if (nStr.IsDigit() && mStr.IsDigit())
		{
			Int_t n = nStr.Atoi();
			Int_t m = mStr.Atoi();

			Int_t nDigits = nStr.Length();
			TString formatStr;
			formatStr.Form("%%s%%0%dd%%s", nDigits);

			AliDebug(2, Form("Format string is %s", formatStr.Data()));

			for (Int_t current = n; current<=m; ++current)
			{
				TString newEntry;
				newEntry.Form(formatStr.Data(), before.Data(), current, after.Data());

				AliDebug(2, Form("Calling recursive with %s", newEntry.Data()));

				// and go recursive
				ExpandAndAdd(target, newEntry);
			}

			// return here because we processed the entries already recursively.
			return;
		}
	}

	AliDebug(2, Form("Adding name %s", entry));
	target->AddLast(new TObjString(entry));
}

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleDCSConfigHolder::~AliShuttleDCSConfigHolder()
{
// destructor of the shuttle configuration holder

	delete fDCSAliases;
	delete fDCSDataPoints;
	delete fDCSAliasesComp;
	delete fDCSDataPointsComp;	
}

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleDetConfigHolder::AliShuttleDetConfigHolder(const TLDAPEntry* entry):
fDetector(""),
fDCSConfig(),
fResponsibles(0),
fIsValid(kFALSE),
fSkipDCSQuery(kFALSE),
fStrictRunOrder(kFALSE)
{
// constructor of the shuttle configuration holder

	TLDAPAttribute* anAttribute;
	
	fResponsibles = new TObjArray();
	fResponsibles->SetOwner(1);
	fDCSConfig = new TObjArray();
	fDCSConfig->SetOwner(1);

	anAttribute = entry->GetAttribute("det"); // MUST
        if (!anAttribute)
	{
		AliError(Form("No \"det\" attribute!"));
		return;
        }
	fDetector = anAttribute->GetValue();

	anAttribute = entry->GetAttribute("strictRunOrder"); // MAY
        if (!anAttribute)
	{
		AliWarning(Form("%s did not set strictRunOrder flag - the default is FALSE",
				fDetector.Data()));
        } else {
		TString strictRunStr = anAttribute->GetValue();
		if (!(strictRunStr == "0" || strictRunStr == "1"))
		{
			AliError("strictRunOrder flag must be 0 or 1!");
			return;
		}
		fStrictRunOrder = (Bool_t) strictRunStr.Atoi();
	}

	anAttribute = entry->GetAttribute("responsible"); // MAY
        if (!anAttribute)
	{
		AliDebug(2, "Warning! No \"responsible\" attribute!");
        }
	else
	{
		const char* aResponsible;
		while ((aResponsible = anAttribute->GetValue()))
		{
			fResponsibles->AddLast(new TObjString(aResponsible));
		}
	}

	fIsValid = kTRUE;
}

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleDetConfigHolder::~AliShuttleDetConfigHolder()
{
// destructor of the shuttle configuration holder

	delete fResponsibles;
	delete fDCSConfig;
}

//______________________________________________________________________________________________
const char* AliShuttleConfig::AliShuttleDetConfigHolder::GetDCSHost(Int_t iServ) const
{
	//
	// returns DCS server host 
	//
	
	if (iServ < 0 || iServ >= GetNServers()) return 0;
	
	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetDCSHost();
}

//______________________________________________________________________________________________
Int_t AliShuttleConfig::AliShuttleDetConfigHolder::GetDCSPort(Int_t iServ) const
{
	//
	// returns DCS server port 
	//
	
	if (iServ < 0 || iServ >= GetNServers()) return 0;
	
	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetDCSPort();
}

//______________________________________________________________________________________________
Int_t AliShuttleConfig::AliShuttleDetConfigHolder::GetMultiSplit(Int_t iServ) const
{
	//
	// returns DCS "multi split" value
	//
	
	if (iServ < 0 || iServ >= GetNServers()) return 0;
	
	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetMultiSplit();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::AliShuttleDetConfigHolder::GetDCSAliases(Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector
	//
	
	if (iServ < 0 || iServ >= GetNServers()) return 0;
	
	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetDCSAliases();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::AliShuttleDetConfigHolder::GetDCSDataPoints(Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector
	//

	if (iServ < 0 || iServ >= GetNServers()) return 0;

	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetDCSDataPoints();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::AliShuttleDetConfigHolder::GetCompactDCSAliases(Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector (Compact style)
	//

	if (iServ < 0 || iServ >= GetNServers()) return 0;

	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetCompactDCSAliases();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::AliShuttleDetConfigHolder::GetCompactDCSDataPoints(Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector (Compact style)
	//

	if (iServ < 0 || iServ >= GetNServers()) return 0;

	AliShuttleDCSConfigHolder *aHolder = dynamic_cast<AliShuttleDCSConfigHolder *> 
							(fDCSConfig->At(iServ));
	
	return aHolder->GetCompactDCSDataPoints();
}

//______________________________________________________________________________________________
void AliShuttleConfig::AliShuttleDetConfigHolder::AddDCSConfig(AliShuttleDCSConfigHolder* holder)
{
	//
	// adds a DCS configuration set in the array of DCS configurations
	// 
	
	if(!holder) return;
	fDCSConfig->AddLast(holder);
}

ClassImp(AliShuttleConfig)

//______________________________________________________________________________________________
AliShuttleConfig::AliShuttleConfig(const char* host, Int_t port,
	const char* binddn, const char* password, const char* basedn):
	fConfigHost(host), 
	fAlienPath(""), 
	fDAQlbHost(""), 
	fDAQlbPort(), 
	fDAQlbUser(""), 
	fDAQlbPass(""),
	fDAQlbDB(""), 
	fDAQlbTable(""), 
	fShuttlelbTable(""), 
	fRunTypelbTable(""),
	fPasswdFilePath(""),
	fTerminateFilePath(""), 
	fMaxRetries(0),
	fPPTimeOut(0), 
	fDCSTimeOut(0), 
	fDCSRetries(0), 
	fDCSQueryOffset(0),
	fDCSDelay(0),
	fPPMaxMem(0), 
	fMonitorHost(""), 
	fMonitorTable(""), 
	fTriggerWait(3600),
	fShuttleFileSystem("/"),
	fFreeDiskWarningThreshold(20),
	fFreeDiskFatalThreshold(10),
	fRunMode(kTest),
	fDetectorMap(), 
	fDetectorList(),
	fAdmin(),
	fShuttleInstanceHost(""), 
	fProcessedDetectors(), 
	fKeepDCSMap(kFALSE),
	fKeepTempFolder(kFALSE),
	fSendMail(kFALSE),
	fProcessAll(kFALSE), 
	fIsValid(kFALSE)

{
	//
	// host: ldap server host
	// port: ldap server port
	// binddn: binddn used for ldap binding (simple bind is used!).
	// password: password for binddn
	// basedn: this is basedn whose childeren entries which have
	//

	fDetectorMap.SetOwner(1);
	fDetectorList.SetOwner(0); //fDetectorList and fDetectorMap share the same object!
	fProcessedDetectors.SetOwner();
	
	for (int i=0; i<6; i++)
	{
		fAdmin[i] = new TObjArray(); 
		fAdmin[i]->SetOwner();
	}

	TLDAPServer aServer(host, port, binddn, password, 3);

	if (!aServer.IsConnected()) {
		AliError(Form("Can't connect to ldap server %s:%d",
				host, port));
		return;
	}

	// reads configuration for the shuttle running on this machine
	
	TLDAPResult* aResult = 0;
	TLDAPEntry* anEntry = 0;
	
	TList dcsList;
	dcsList.SetOwner(1);
	TList detList;
	detList.SetOwner(1);
	TList globalList;
	globalList.SetOwner(1);
	TList sysList;
	sysList.SetOwner(1);
	TList hostList;
	hostList.SetOwner(1);
	
	aResult = aServer.Search(basedn, LDAP_SCOPE_SUBTREE, 0, 0);
	
	if (!aResult) 
	{
		AliError(Form("Can't find configuration with base DN: %s !", basedn));
		return;
	}
	
	while ((anEntry = aResult->GetNext())) 
	{
		TString dn = anEntry->GetDn();
		
		if (dn.BeginsWith("dcsHost=")) 
		{
			dcsList.Add(anEntry);
		} 
		else if (dn.BeginsWith("det="))
		{
			detList.Add(anEntry);
		}
		else if (dn.BeginsWith("name=globalConfig"))
		{
			globalList.Add(anEntry);
		}
		else if (dn.BeginsWith("system="))
		{
			sysList.Add(anEntry);
		}
		else if (dn.BeginsWith("shuttleHost="))
		{
			hostList.Add(anEntry);
		}
		else 
		{
			delete anEntry;
		}
	
	}
	delete aResult;
	
	Int_t result=0;
	
	result += SetGlobalConfig(&globalList);
	result += SetSysConfig(&sysList);
	result += SetPasswords();
	result += SetDetConfig(&detList,&dcsList);
	result += SetHostConfig(&hostList);
	
	if(result) 
	{
		AliError("Configuration is INVALID!");
	}
	else 
	{
		fIsValid = kTRUE;
	}
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
Int_t AliShuttleConfig::GetNServers(const char* detector) const
{
	//
	// returns number of DCS servers for detector
	//
	
	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
	if (!aHolder) {
		AliError(Form("There isn't configuration for detector: %s",
			detector));
		return 0;
	}

	return aHolder->GetNServers();
}


//______________________________________________________________________________________________
const char* AliShuttleConfig::GetDCSHost(const char* detector, Int_t iServ) const
{
	//
	// returns i-th DCS server host used by particular detector
	//
	
	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
	if (!aHolder) {
		AliError(Form("There isn't configuration for detector: %s",
			detector));
		return NULL;
	}

	return aHolder->GetDCSHost(iServ);
}

//______________________________________________________________________________________________
Int_t AliShuttleConfig::GetDCSPort(const char* detector, Int_t iServ) const
{
	//
        // returns i-th DCS server port used by particular detector
        //


	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return 0;
        }

	return aHolder->GetDCSPort(iServ);
}

//______________________________________________________________________________________________
Int_t AliShuttleConfig::GetMultiSplit(const char* detector, Int_t iServ) const
{
	//
        // returns i-th DCS "multi split" value
        //


	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return 0;
        }

	return aHolder->GetMultiSplit(iServ);
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetDCSAliases(const char* detector, Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the i-th set of aliases
	// which used for data retrieval for particular detector
	//

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetDCSAliases(iServ);
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetDCSDataPoints(const char* detector, Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector
	//

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetDCSDataPoints(iServ);
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetCompactDCSAliases(const char* detector, Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the i-th set of aliases
	// which used for data retrieval for particular detector (Compact style)
	//

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetCompactDCSAliases(iServ);
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetCompactDCSDataPoints(const char* detector, Int_t iServ) const
{
	//
	// returns collection of TObjString which represents the set of aliases
	// which used for data retrieval for particular detector (Compact style)
	//

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetCompactDCSDataPoints(iServ);
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetResponsibles(const char* detector) const
{
	//
	// returns collection of TObjString which represents the list of mail addresses
	// of the detector's responsible(s)
	//

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder) {
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return NULL;
        }

	return aHolder->GetResponsibles();
}

//______________________________________________________________________________________________
const TObjArray* AliShuttleConfig::GetAdmins(Int_t sys) const
{
	//
	// returns collection of TObjString which represents the list of mail addresses
	// of the system's administrator(s) (DAQ, DCS, HLT, Global, Amanda, DQM)
	//

	if (sys < 0 || sys > 5) return 0;
	return fAdmin[sys];
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
Bool_t AliShuttleConfig::StrictRunOrder(const char* detector) const
{
	// return TRUE if detector wants strict run ordering of stored data

	AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) fDetectorMap.GetValue(detector);
        if (!aHolder)
	{
                AliError(Form("There isn't configuration for detector: %s",
                        detector));
                return kTRUE;
        }

	return aHolder->StrictRunOrder();
}

//______________________________________________________________________________________________
UInt_t AliShuttleConfig::SetGlobalConfig(TList* list)
{
	// Set the global configuration (DAQ Logbook + preprocessor monitoring settings)


	TLDAPEntry* anEntry = 0;
	TLDAPAttribute* anAttribute = 0;
	
	if (list->GetEntries() == 0) 
	{
		AliError("Global configuration not found!");
		return 1;
	} 
	else if (list->GetEntries() > 1)
	{
		AliError("More than one global configuration found!");
		return 2;
	}
	
	anEntry = dynamic_cast<TLDAPEntry*> (list->At(0));
	
	if (!anEntry)
	{
		AliError("Unexpected! Global list does not contain a TLDAPEntry");
		return 3;
	} 
	
	
	anAttribute = anEntry->GetAttribute("AlienPath");
	if (!anAttribute) {
		AliError("Can't find AlienPath attribute!");
		return 4;
	}
	fAlienPath = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("daqLbHost");
	if (!anAttribute) {
		AliError("Can't find daqLbHost attribute!");
		return 4;
	}
	fDAQlbHost = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("daqLbPort"); // MAY
	if (anAttribute)
	{
		fDAQlbPort = ((TString) anAttribute->GetValue()).Atoi();
	} else {
		fDAQlbPort = 3306; // mysql
	}

	anAttribute = anEntry->GetAttribute("daqLbUser");
	if (!anAttribute) {
		AliError("Can't find daqLbUser attribute!");
		return 4;
	}
	fDAQlbUser = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("daqLbDB");
	if (!anAttribute) {
		AliError("Can't find daqLbDB attribute!");
		return 4;
	}
	fDAQlbDB = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("daqLbTable");
	if (!anAttribute) {
		AliError("Can't find daqLbTable attribute!");
		return 4;
	}
	fDAQlbTable = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("shuttleLbTable");
	if (!anAttribute) {
		AliError("Can't find shuttleLbTable attribute!");
		return 4;
	}
	fShuttlelbTable = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("runTypeLbTable");
	if (!anAttribute) {
		AliError("Can't find runTypeLbTable attribute!");
		return 4;
	}
	fRunTypelbTable = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("ppmaxRetries");
	if (!anAttribute) {
		AliError("Can't find ppmaxRetries attribute!");
		return 4;
	}
	TString tmpStr = anAttribute->GetValue();
	fMaxRetries = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("terminateFilePath");
	if (anAttribute) {
		fTerminateFilePath = anAttribute->GetValue();
	}

	anAttribute = anEntry->GetAttribute("ppTimeOut");
	if (!anAttribute) {
		AliError("Can't find ppTimeOut attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fPPTimeOut = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("dcsTimeOut");
	if (!anAttribute) {
		AliError("Can't find dcsTimeOut attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fDCSTimeOut = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("nDCSretries");
	if (!anAttribute) {
		AliError("Can't find dcsTimeOut attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fDCSRetries = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("DCSQueryOffset");
	if (!anAttribute) {
		AliError("Can't find DCSQueryOffset attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fDCSQueryOffset = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("DCSDelay");
	if (!anAttribute) {
		AliError("Can't find DCSDelay attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fDCSDelay = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("ppMaxMem");
	if (!anAttribute) {
		AliError("Can't find ppMaxMem attribute!");
		return 4;
	}
	tmpStr = anAttribute->GetValue();
	fPPMaxMem = tmpStr.Atoi();
	
	anAttribute = anEntry->GetAttribute("monitorHost");
	if (!anAttribute) {
		AliError("Can't find monitorHost attribute!");
		return 4;
	}
	fMonitorHost = anAttribute->GetValue();
	
	anAttribute = anEntry->GetAttribute("monitorTable");
	if (!anAttribute) {
		AliError("Can't find monitorTable attribute!");
		return 4;
	}
	fMonitorTable = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("triggerWait"); // MAY
	if (!anAttribute) {
		AliWarning(Form("triggerWait not set! default = %d", fTriggerWait));
	}
	tmpStr = anAttribute->GetValue();
	fTriggerWait = tmpStr.Atoi();

        anAttribute = anEntry->GetAttribute("ShuttleFileSystem"); 
	if (!anAttribute) {
		AliWarning(Form("ShuttleFileSystem not set! default = %s", fShuttleFileSystem.Data()));
	}
	fShuttleFileSystem = anAttribute->GetValue();

	anAttribute = anEntry->GetAttribute("FreeDiskWarningThreshold"); // MAY
	if (!anAttribute) {
		AliWarning(Form("FreeDiskWarningThreshold not set! default = %d", fFreeDiskWarningThreshold));
	}
	tmpStr = anAttribute->GetValue();
	fFreeDiskWarningThreshold = tmpStr.Atoi();

	anAttribute = anEntry->GetAttribute("FreeDiskFatalThreshold"); // MAY
	if (!anAttribute) {
		AliWarning(Form("FreeDiskFatalThreshold not set! default = %d", fFreeDiskFatalThreshold));
	}
	tmpStr = anAttribute->GetValue();
	fFreeDiskFatalThreshold = tmpStr.Atoi();		

	anAttribute = anEntry->GetAttribute("mode");
	if (!anAttribute) {
		AliWarning("Run mode not set! Running in test mode.");
	} else {
	  tmpStr = anAttribute->GetValue();
	  if (tmpStr == "test")
	  {
	    fRunMode = kTest;
	  } else if (tmpStr == "prod") {
	    fRunMode = kProd;
	  } else {
	    AliWarning(Form("Not a valid run mode: %s", tmpStr.Data()));		
	    AliWarning("Valid run modes are \"test\" and \"prod\". Running in test mode.");
	  }
	}

	anAttribute = anEntry->GetAttribute("keepDCSMap"); // MAY
        if (!anAttribute)
	{
		AliWarning("keepDCSMap flag not set - default is FALSE");
        } else {
		TString keepDCSMapStr = anAttribute->GetValue();
		if (!(keepDCSMapStr == "0" || keepDCSMapStr == "1"))
		{
			AliError("keepDCSMap flag must be 0 or 1!");
			return 4;
		}
		fKeepDCSMap = (Bool_t) keepDCSMapStr.Atoi();
	}
	
	anAttribute = anEntry->GetAttribute("keepTempFolder"); // MAY
        if (!anAttribute)
	{
		AliWarning("keepTempFolder flag not set - default is FALSE");
        } else {
		TString keepTempFolderStr = anAttribute->GetValue();
		if (!(keepTempFolderStr == "0" || keepTempFolderStr == "1"))
		{
			AliError("keepTempFolder flag must be 0 or 1!");
			return 4;
		}
		fKeepTempFolder = (Bool_t) keepTempFolderStr.Atoi();
	}
	

	anAttribute = anEntry->GetAttribute("shuttleAdmin"); // MAY
        if (!anAttribute)
	{
		AliDebug(2, "Warning! No \"shuttleAdmin\" attribute!");
        }
	else
	{
		const char* anAdmin;
		while ((anAdmin = anAttribute->GetValue()))
		{
			fAdmin[kGlobal]->AddLast(new TObjString(anAdmin));
		}
	}

	anAttribute = anEntry->GetAttribute("amandaAdmin"); // MAY
        if (!anAttribute)
	{
		AliDebug(2, "Warning! No \"amandaAdmin\" attribute!");
        }
	else
	{
		const char* anAdmin;
		while ((anAdmin = anAttribute->GetValue()))
		{
			fAdmin[kAmanda]->AddLast(new TObjString(anAdmin));
		}
	}
	
	anAttribute = anEntry->GetAttribute("sendMail"); // MAY
        if (!anAttribute)
	{
		AliWarning("sendMail flag not set - default is FALSE");
        } else {
		TString sendMailStr = anAttribute->GetValue();
		if (!(sendMailStr == "0" || sendMailStr == "1"))
		{
			AliError("sendMail flag must be 0 or 1!");
			return 4;
		}
		fSendMail = (Bool_t) sendMailStr.Atoi();
	}
						
	anAttribute = anEntry->GetAttribute("passwdFilePath");
	if (!anAttribute) {
		AliError("Can't find Passwords File Path attribute!");
		return 4;
	}
	fPasswdFilePath = anAttribute->GetValue();

	return 0;
}

//______________________________________________________________________________________________
UInt_t AliShuttleConfig::SetSysConfig(TList* list)
{
	// Set the online FXS configuration (DAQ + DCS + HLT + DQM)


	TLDAPEntry* anEntry = 0;
	TLDAPAttribute* anAttribute = 0;
	
	if (list->GetEntries() != 4) 
	{
		AliError(Form("Wrong number of online systems found: %d !", list->GetEntries()));
		return 1;
	} 

	TIter iter(list);
	Int_t count = 0;
	SystemCode iSys=kDAQ;
	
	while ((anEntry = dynamic_cast<TLDAPEntry*> (iter.Next())))
	{
		anAttribute = anEntry->GetAttribute("system");
		TString sysName = anAttribute->GetValue();
		
		if (sysName == "DAQ") 
		{
			iSys = kDAQ;
			count += 1;
		}
		else if (sysName == "DCS")
		{
			iSys = kDCS;
			count += 10;
		}
		else if (sysName == "HLT")
		{
			iSys = kHLT;
			count += 100;
		}
		else if (sysName == "DQM")
		{
			iSys = kDQM; 
			count += 1000;
		}
		
		anAttribute = anEntry->GetAttribute("dbHost");
		if (!anAttribute) {
			AliError(Form ("Can't find dbHost attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}
		fFXSdbHost[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("dbPort"); // MAY
		if (anAttribute)
		{
			fFXSdbPort[iSys] = ((TString) anAttribute->GetValue()).Atoi();
		} else {
			fFXSdbPort[iSys] = 3306; // mysql
		}

		anAttribute = anEntry->GetAttribute("dbUser");
		if (!anAttribute) {
			AliError(Form ("Can't find dbUser attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}
		fFXSdbUser[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("dbName");
		if (!anAttribute) {
			AliError(Form ("Can't find dbName attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}

		fFXSdbName[iSys] = anAttribute->GetValue();
		anAttribute = anEntry->GetAttribute("dbTable");
		if (!anAttribute) {
			AliError(Form ("Can't find dbTable attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}
		fFXSdbTable[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("fxsHost");
		if (!anAttribute) {
			AliError(Form ("Can't find fxsHost attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}
		fFXSHost[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("fxsPort"); // MAY
		if (anAttribute)
		{
			fFXSPort[iSys] = ((TString) anAttribute->GetValue()).Atoi();
		} else {
			fFXSPort[iSys] = 22; // scp port number
		}

		anAttribute = anEntry->GetAttribute("fxsUser");
		if (!anAttribute) {
			AliError(Form ("Can't find fxsUser attribute for %s!!",
						AliShuttleInterface::GetSystemName(iSys)));
			return 5;
		}
		fFXSUser[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("fxsPasswd");
		if (anAttribute) fFXSPass[iSys] = anAttribute->GetValue();

		anAttribute = anEntry->GetAttribute("fxsBaseFolder");
		if (anAttribute) fFXSBaseFolder[iSys] = anAttribute->GetValue();
	
		anAttribute = anEntry->GetAttribute("fxsAdmin"); // MAY
        	if (!anAttribute)
		{
			AliDebug(2, "Warning! No \"fxsAdmin\" attribute!");
        	}
		else
		{
			const char* anAdmin;
			while ((anAdmin = anAttribute->GetValue()))
			{
				fAdmin[iSys]->AddLast(new TObjString(anAdmin));
			}
		}
	
	}
	
	if(count != 1111) {
		AliError(Form("Wrong system configuration! (code = %d)", count));
		return 6;
	}
	
	return 0;
}

//______________________________________________________________________________________________
UInt_t AliShuttleConfig::SetPasswords(){
	
	AliInfo("Setting Passwords");

	// Retrieving Passwords for DAQ lb, DAQ/DCS/HLT/DQM FXS

	ifstream *inputfile = new ifstream(fPasswdFilePath.Data());
	if (!*inputfile) {
		AliError(Form("Error opening file %s !", fPasswdFilePath.Data()));
		inputfile->close();
		delete inputfile;
		return 1;
	}

	TString line;
	Int_t nPwd=0;
	Int_t nPwdFake=0;
	while (line.ReadLine(*inputfile)) {
		TObjArray *tokens = line.Tokenize(" \t");
		TString system = ((TObjString *)tokens->At(0))->String(); 
		TString password = ((TObjString *)tokens->At(1))->String();
		if (system.Contains("DAQ_LB")){
			fDAQlbPass=password;
			nPwd++;
			AliDebug(3,Form("DAQ_LB: Password %s for %s found", password.Data(), system.Data()));
		}
		else if (system.Contains("DAQ_DB")){
			fFXSdbPass[0]=password;
			nPwd++;
			AliDebug(3,Form("DAQ_DB: Password %s for %s found", password.Data(), system.Data()));
		}
		else if (system.Contains("DCS_DB")){
			fFXSdbPass[1]=password;
			nPwd++;
			AliDebug(3,Form("DCS_DB: Password %s for %s found", password.Data(), system.Data()));
		}
		else if (system.Contains("HLT_DB")){
			fFXSdbPass[2]=password;
			nPwd++;
			AliDebug(3,Form("HLT_DB: Password %s for %s found", password.Data(), system.Data()));
		}
		else if (system.Contains("DQM_DB")){
			fFXSdbPass[3]=password;
			nPwd++;
			AliDebug(3,Form("DQM_DB: Password %s for %s found", password.Data(), system.Data()));
		}
		else {
			nPwdFake++;
			AliDebug(3,Form("%i fake line(s) found in file %s", nPwdFake, fPasswdFilePath.Data()));
			continue;
		}
		delete tokens;
	}

	inputfile->close();
	delete inputfile;

	if (nPwd!=5){
		AliError(Form("Wrong file for DAQ Logbook password found %s (some passwors missing), please Check!", fPasswdFilePath.Data()));
		return 2;
	}

	return 0;

}
//______________________________________________________________________________________________
UInt_t AliShuttleConfig::SetDetConfig(TList* detList, TList* dcsList)
{
	// Set the detector configuration (general settings + DCS amanda server and alias/DP lists)

	TLDAPEntry* anEntry = 0;
	
	TIter iter(detList);
	while ((anEntry = dynamic_cast<TLDAPEntry*> (iter.Next())))
	{
		
		AliShuttleDetConfigHolder* detHolder = new AliShuttleDetConfigHolder(anEntry);

		if (!detHolder->IsValid()) {
			AliError("Detector configuration error!");
			delete detHolder;
			return 7;
		}

		TObjString* detStr = new TObjString(detHolder->GetDetector());
		
		// Look for DCS Configuration
		TIter dcsIter(dcsList);
		TLDAPEntry *dcsEntry = 0;
		while ((dcsEntry = dynamic_cast<TLDAPEntry*> (dcsIter.Next())))
		{
			TString dn = dcsEntry->GetDn();
			if(dn.Contains(detStr->GetName())) {
				AliDebug(2, Form("Found DCS configuration: dn = %s",dn.Data()));
				AliShuttleDCSConfigHolder* dcsHolder = new AliShuttleDCSConfigHolder(dcsEntry);
				if (!dcsHolder->IsValid()) {
					AliError("DCS configuration error!");
					delete detHolder;
					delete dcsHolder;
					return 7;
				}
				detHolder->AddDCSConfig(dcsHolder);
			}
		}
		
		
		fDetectorMap.Add(detStr, detHolder);
		fDetectorList.AddLast(detStr);
	}
	
	return 0;
}

//______________________________________________________________________________________________
UInt_t AliShuttleConfig::SetHostConfig(TList* list)
{
	// Set the Shuttle machines configuration (which detectors processes each machine)
	
	TLDAPEntry* anEntry = 0;
	TLDAPAttribute* anAttribute = 0;
	
	fShuttleInstanceHost = gSystem->HostName();
	
	TIter iter(list);
	while ((anEntry = dynamic_cast<TLDAPEntry*> (iter.Next())))
	{
	
		TString dn(anEntry->GetDn());
		if (!dn.Contains(Form("shuttleHost=%s", fShuttleInstanceHost.Data()))) continue;
		
		if (!fProcessAll)
		{
			anAttribute = anEntry->GetAttribute("detectors");
			const char *detName;
			while((detName = anAttribute->GetValue())){
				TObjString *objDet= new TObjString(detName);
				fProcessedDetectors.Add(objDet);
			}
		}
	}	
	
	return 0;
}


//______________________________________________________________________________________________
void AliShuttleConfig::Print(Option_t* option) const
{
// print configuration
// options : "": print configuration for all detectors, aliases and DPs in compacted format
// 	     "uncompact": print configuration for all detectors, aliases and DPs in uncompacted format
// 	     "DET": print configuration for DET, aliases and DPs in compacted format
//	     "DET, uncompact": print configuration for DET, aliases and DPs in uncompacted format

	TString result;
	result += '\n';
	
	TString mode = "test";
	if (fRunMode == kProd) mode = "production";

	result += "########################################################################\n";
	result += Form(" Shuttle configuration from %s - Run Mode: <%s> \n", 
					fConfigHost.Data(), mode.Data());
	result += "########################################################################\n";
	result += Form("\nShuttle running on %s \n", fShuttleInstanceHost.Data());

	if(fProcessAll) {
		result += Form("All detectors will be processed! \n");
	} else {
		result += "Detectors processed by this host: ";
		TIter it(&fProcessedDetectors);
		TObjString* aDet;
		while ((aDet = (TObjString*) it.Next())) {
			result += Form("%s ", aDet->String().Data());
		}
		result += "\n";
	}

	result += Form("PP time out = %d - DCS time out = %d - max PP mem size = %d KB - max retries = %d "
		       "- DIM trigger waiting timeout = %d\n", 
				fPPTimeOut, fDCSTimeOut, fPPMaxMem, fMaxRetries, fTriggerWait);
	result += Form("FLAGS: keepDCSMap = %d - keepTempFolder = %d - SendMail = %d \n", 
				fKeepDCSMap, fKeepTempFolder, fSendMail);
	const TObjArray* shuttleAdmins = GetAdmins(kGlobal);
	if (shuttleAdmins->GetEntries() != 0)
	{
		result += "SHUTTLE administrator(s): ";
		TIter it(shuttleAdmins);
		TObjString* anAdmin;
		while ((anAdmin = (TObjString*) it.Next()))
		{
			result += Form("%s ", anAdmin->String().Data());
		}
		result += "\n";
	}
	const TObjArray* amandaAdmins = GetAdmins(kAmanda);
	if (amandaAdmins->GetEntries() != 0)
	{
		result += "Amanda server administrator(s): ";
		TIter it(amandaAdmins);
		TObjString* anAdmin;
		while ((anAdmin = (TObjString*) it.Next()))
		{
			result += Form("%s ", anAdmin->String().Data());
		}
		result += "\n\n";
	}
	result += "------------------------------------------------------\n";

	result += Form("Logbook Configuration \n\n \tHost: %s:%d; \tUser: %s; ",
		fDAQlbHost.Data(), fDAQlbPort, fDAQlbUser.Data());

	result += Form("\tDB: %s; \tTables: %s, %s, %s\n",
		fDAQlbDB.Data(), fDAQlbTable.Data(), fShuttlelbTable.Data(), fRunTypelbTable.Data());

	result += Form("Terminate file path: %s\n", fTerminateFilePath.Data());

	result += "\n\n";
	
	result += "------------------------------------------------------\n";
	result += "FXS configuration\n\n";

	for(int iSys=0;iSys<4;iSys++){
		result += Form("*** %s ***\n", AliShuttleInterface::GetSystemName(iSys));
		result += Form("\tDB  host: %s:%d; \tUser: %s; \tName: %s; \tTable: %s\n",
						fFXSdbHost[iSys].Data(), fFXSdbPort[iSys], fFXSdbUser[iSys].Data(),
						fFXSdbName[iSys].Data(), fFXSdbTable[iSys].Data());
		result += Form("\tFXS host: %s:%d; \tUser: %s\n", fFXSHost[iSys].Data(), fFXSPort[iSys],
						fFXSUser[iSys].Data());
		const TObjArray* fxsAdmins = GetAdmins(iSys);
		if (fxsAdmins->GetEntries() != 0)
		{
			result += "\tAdministrator(s): ";
			TIter it(fxsAdmins);
			TObjString* anAdmin;
			while ((anAdmin = (TObjString*) it.Next()))
			{
				result += Form("%s ", anAdmin->String().Data());
			}
			result += "\n\n";
		}
	}

	result += "------------------------------------------------------\n";
	result += "MonaLisa configuration\n\n";
	
	result += Form("\tHost: %s; \tTable name: %s",
		fMonitorHost.Data(), fMonitorTable.Data());
		
	result += "\n\n";
		
	TString optStr(option);

	result += "------------------------------------------------------\n";
	result += "Detector-specific configuration\n\n";
	
	TIter iter(fDetectorMap.GetTable());
	TPair* aPair = 0;
	
	while ((aPair = (TPair*) iter.Next())) {
		AliShuttleDetConfigHolder* aHolder = (AliShuttleDetConfigHolder*) aPair->Value();
		if (optStr != "" && !optStr.Contains(aHolder->GetDetector()) && 
				optStr.CompareTo("uncompact",TString::kIgnoreCase) != 0 )
				continue;
		
		result += Form("*** %s *** \n", aHolder->GetDetector());

		const TObjArray* responsibles = aHolder->GetResponsibles();
		if (responsibles->GetEntries() != 0)
		{
			result += "\tDetector responsible(s): ";
			TIter it(responsibles);
			TObjString* aResponsible;
			while ((aResponsible = (TObjString*) it.Next()))
			{
				result += Form("%s ", aResponsible->String().Data());
			}
			result += "\n";
		}

		result += Form("\tStrict run ordering: %s \n\n", aHolder->StrictRunOrder() ? "YES" : "NO");
		
		const TObjArray* dcsConfig = aHolder->GetDCSConfig();
		
		AliShuttleDCSConfigHolder* dcsHolder = 0;
		TIter dcsIter(dcsConfig);
		Int_t count=0;
		while ((dcsHolder = dynamic_cast<AliShuttleDCSConfigHolder*> (dcsIter.Next())))
		{
			result += Form("\tAmanda server [%d]: %s:%d - MultiSplit = %d\n", count,
				dcsHolder->GetDCSHost(), dcsHolder->GetDCSPort(), dcsHolder->GetMultiSplit());

			const TObjArray* aliases = 0;
			if (optStr.Contains("uncompact",TString::kIgnoreCase))
			{
				aliases = dcsHolder->GetDCSAliases();
			} else {
				aliases = dcsHolder->GetCompactDCSAliases();
			}

			if (aliases->GetEntries() != 0)
			{
				result += Form("\tDCS Aliases [%d]: ", count);
				TIter it(aliases);
				TObjString* anAlias;
				while ((anAlias = (TObjString*) it.Next()))
				{
					result += Form("%s ", anAlias->String().Data());
				}
				result += "\n";
			}

			const TObjArray* dataPoints = 0;
			if (optStr.Contains("uncompact",TString::kIgnoreCase))
			{
				dataPoints = dcsHolder->GetDCSDataPoints();
			} else {
				dataPoints = dcsHolder->GetCompactDCSDataPoints();
			}
			if (dataPoints->GetEntries() != 0)
			{
				result += Form("\tDCS Data Points [%d]: ", count);
				TIter it(dataPoints);
				TObjString* aDataPoint;
				while ((aDataPoint = (TObjString*) it.Next())) {
					result += Form("%s ", aDataPoint->String().Data());
				}
				result += "\n";
			}
			count++;
			result += "\n";
		}
	}
	if(!fIsValid) result += "\n\n********** !!!!! Configuration is INVALID !!!!! **********\n";

	AliInfo(result);
}
