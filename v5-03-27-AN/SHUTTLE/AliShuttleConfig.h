#ifndef ALI_SHUTTLE_CONFIG_H
#define ALI_SHUTTLE_CONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class keeps the AliShuttle configuration.
// It reads the configuration for LDAP server.
// For more info see AliShuttleConfig.cxx
//

#include <TObject.h>
#include <TString.h>
#include <TObjArray.h>
#include <TMap.h>
#include <TLDAPServer.h>

class AliShuttleConfig: public TObject {
public:
	enum RunMode {kTest=0, kProd};
        enum SystemCode { kDAQ = 0, kDCS, kHLT, kDQM, kGlobal, kAmanda};
	
	AliShuttleConfig(const char* host, Int_t port = LDAP_PORT,
			const char* binddn = 0, const char* password = 0,
			const char* basedn = "o=alice,dc=cern,dc=ch");
	virtual ~AliShuttleConfig();

	Bool_t IsValid() const {return fIsValid;};

	const char* GetConfigHost() const {return fConfigHost.Data();}

	const char* GetAlienPath() const {return fAlienPath.Data();}

	const char* GetDAQlbHost() const {return fDAQlbHost.Data();}
	UInt_t 	    GetDAQlbPort() const {return fDAQlbPort;}
	const char* GetDAQlbUser() const {return fDAQlbUser.Data();}
	const char* GetDAQlbPass() const {return fDAQlbPass.Data();}
	const char* GetDAQlbDB() const {return fDAQlbDB.Data();}
	const char* GetDAQlbTable() const {return fDAQlbTable.Data();}
	const char* GetShuttlelbTable() const {return fShuttlelbTable.Data();}
	const char* GetRunTypelbTable() const {return fRunTypelbTable.Data();}

	const char* GetFXSHost(Int_t sys) const {return fFXSHost[sys].Data();}
	UInt_t 	    GetFXSPort(Int_t sys) const {return fFXSPort[sys];}
	const char* GetFXSUser(Int_t sys) const {return fFXSUser[sys].Data();}
	const char* GetFXSPass(Int_t sys) const {return fFXSPass[sys].Data();}
	const char* GetFXSBaseFolder(Int_t sys) const {return fFXSBaseFolder[sys].Data();}

	const char* GetFXSdbHost(Int_t sys) const {return fFXSdbHost[sys].Data();}
	UInt_t 	    GetFXSdbPort(Int_t sys) const {return fFXSdbPort[sys];}
	const char* GetFXSdbUser(Int_t sys) const {return fFXSdbUser[sys].Data();}
	const char* GetFXSdbPass(Int_t sys) const {return fFXSdbPass[sys].Data();}
	const char* GetFXSdbName(Int_t sys) const {return fFXSdbName[sys].Data();}
	const char* GetFXSdbTable(Int_t sys) const {return fFXSdbTable[sys].Data();}

	const TObjArray* GetAdmins(Int_t sys) const;
	
	Bool_t SendMail() const { return fSendMail; }
	
	Int_t GetMaxRetries() const { return fMaxRetries; }

	Int_t GetPPTimeOut() const { return fPPTimeOut; }
	Int_t GetDCSTimeOut() const { return fDCSTimeOut; }
	Int_t GetDCSRetries() const { return fDCSRetries; }
	UInt_t GetDCSQueryOffset() const { return fDCSQueryOffset; }
	UInt_t GetDCSDelay() const { return fDCSDelay; }
	Int_t GetPPMaxMem() const { return fPPMaxMem; }

	const char* GetTerminateFilePath() const { return fTerminateFilePath; }

	Bool_t KeepDCSMap() const { return fKeepDCSMap; }
	Bool_t KeepTempFolder() const { return fKeepTempFolder; }
	
		
	const char* GetMonitorHost() const {return fMonitorHost.Data();}
	const char* GetMonitorTable() const {return fMonitorTable.Data();}

	Int_t GetTriggerWait() const {return fTriggerWait;}
	const char* GetShuttleFileSystem() const {return fShuttleFileSystem.Data();}
	Int_t GetFreeDiskWarningThreshold() const {return fFreeDiskWarningThreshold;}
	Int_t GetFreeDiskFatalThreshold() const {return fFreeDiskFatalThreshold;}
	
	RunMode GetRunMode() const {return fRunMode;}
	
	const TObjArray* GetDetectors() const;

	Bool_t HasDetector(const char* detector) const;
	
	Int_t GetNServers(const char* detector) const;
	
	const char* GetDCSHost(const char* detector, Int_t iServ) const;
	Int_t GetDCSPort(const char* detector, Int_t iServ) const;
	
	Int_t GetMultiSplit(const char* detector, Int_t iServ) const;
	
	const TObjArray* GetDCSAliases(const char* detector, Int_t iServ) const;
	const TObjArray* GetDCSDataPoints(const char* detector, Int_t iServ) const;
	const TObjArray* GetCompactDCSAliases(const char* detector, Int_t iServ) const;
	const TObjArray* GetCompactDCSDataPoints(const char* detector, Int_t iServ) const;
	
	const TObjArray* GetResponsibles(const char* detector) const;
	Bool_t StrictRunOrder(const char* detector) const;

	void SetProcessAll(Bool_t flag=kTRUE) {fProcessAll=flag;}
	Bool_t ProcessAll() const {return fProcessAll;}

	Bool_t HostProcessDetector(const char* detector) const;

	virtual void Print(Option_t* option = NULL) const;

private:

	class AliShuttleDCSConfigHolder: public TObject {
	public:
		AliShuttleDCSConfigHolder(const TLDAPEntry* entry);
		~AliShuttleDCSConfigHolder();
		
		const char* GetDCSHost() const {return fDCSHost;}
		Int_t GetDCSPort() const {return fDCSPort;}
		Int_t GetMultiSplit() const {return fMultiSplit;}
		
		const TObjArray* GetDCSAliases() const {return fDCSAliases;}
		const TObjArray* GetDCSDataPoints() const {return fDCSDataPoints;}
		const TObjArray* GetCompactDCSAliases() const {return fDCSAliasesComp;}
		const TObjArray* GetCompactDCSDataPoints() const {return fDCSDataPointsComp;}
		
		void ExpandAndAdd(TObjArray* target, const char* entry);
		Bool_t IsValid() const {return fIsValid;}

	private:
		AliShuttleDCSConfigHolder(const AliShuttleDCSConfigHolder& other);  // not implemented
		AliShuttleDCSConfigHolder& operator= (const AliShuttleDCSConfigHolder& other);  // not implemented

		TString fDCSHost; 	// Host name of the DCS server
		Int_t 	fDCSPort; 	// port of the DCS server
		Int_t   fMultiSplit;     // N. of DPs to be queried in each query chunk
		TObjArray* fDCSAliases; // List of DCS aliases to be retrieved
		TObjArray* fDCSDataPoints;     // List of DCS data points to be retrieved
		TObjArray* fDCSAliasesComp;    // Compact list of DCS aliases to be printed
		TObjArray* fDCSDataPointsComp; // Compact list of DCS data points to be printed
		Bool_t  fIsValid;  	       // flag for the validity of the configuration

		ClassDef(AliShuttleDCSConfigHolder, 0);
	};

	class AliShuttleDetConfigHolder: public TObject {
	public:
		AliShuttleDetConfigHolder(const TLDAPEntry* entry);
		~AliShuttleDetConfigHolder();

		const char* GetDetector() const {return fDetector.Data();}
		const TObjArray* GetDCSConfig() const {return fDCSConfig;}
		void AddDCSConfig(AliShuttleDCSConfigHolder* holder);

		Int_t GetNServers() const {return fDCSConfig ? fDCSConfig->GetEntries() : 0;}
		
		const char* GetDCSHost(Int_t iServ) const;
		Int_t GetDCSPort(Int_t iServ) const;
		Int_t GetMultiSplit(Int_t iServ) const;
		
		const TObjArray* GetDCSAliases(Int_t iServ) const;
		const TObjArray* GetDCSDataPoints(Int_t iServ) const;
		const TObjArray* GetCompactDCSAliases(Int_t iServ) const;
		const TObjArray* GetCompactDCSDataPoints(Int_t iServ) const;
		
		const TObjArray* GetResponsibles() const {return fResponsibles;}		
		Bool_t IsValid() const {return fIsValid;}
		Bool_t SkipDCSQuery() const {return fSkipDCSQuery;}
		Bool_t StrictRunOrder() const {return fStrictRunOrder;}

	private:
		AliShuttleDetConfigHolder(const AliShuttleDetConfigHolder& other);  // not implemented
		AliShuttleDetConfigHolder& operator= (const AliShuttleDetConfigHolder& other);  // not implemented

		TString fDetector;  	// Detector name
		TObjArray* fDCSConfig;  // Array of DCS configuration objects (AliShuttleDCSConfigHolder)
		TObjArray* fResponsibles; // List of email addresses of the detector's responsible(s)
		Bool_t fIsValid;  	// flag for the validity of the configuration
		Bool_t fSkipDCSQuery; 	// flag - if TRUE (-> DCS config empty) skip DCS archive data query
		Bool_t fStrictRunOrder; // flag - if TRUE connect data in a strict run ordering

		ClassDef(AliShuttleDetConfigHolder, 0);
	};


	AliShuttleConfig(const AliShuttleConfig& other);
	AliShuttleConfig& operator= (const AliShuttleConfig& other);
	UInt_t SetGlobalConfig(TList* globalList);
	UInt_t SetSysConfig(TList* sysList);
	UInt_t SetPasswords();
	UInt_t SetDetConfig(TList* detList, TList* dcsList);
	UInt_t SetHostConfig(TList* hostList);
	
	TString fConfigHost;  		// Host of the Shuttle configuration LDAP server

	TString fAlienPath;		// Alien Path for CDB/Reference objs in the Alien catalogue

	TString fDAQlbHost;		// Host of the DAQ logbook MySQL Server
	UInt_t  fDAQlbPort;		// port of the DAQ logbook MySQL Server
	TString fDAQlbUser;  		// username of the DAQ logbook MySQL Server
	TString fDAQlbPass; 		// password of the DAQ logbook MySQL Server
	TString fDAQlbDB; 		// DB name of the DAQ logbook MySQL Server
	TString fDAQlbTable; 		// Table name of the DAQ logbook 
	TString fShuttlelbTable; 	// Table name of the Shuttle logbook
	TString fRunTypelbTable; 	// Table name of the Run typr logbook
	TString fPasswdFilePath; 	// Path for the local file where the passwords are stored

	TString fFXSHost[4]; 		// Host of the [DAQ, DCS, HLT, DQM] File eXchange Server
	UInt_t  fFXSPort[4]; 		// Port of the [DAQ, DCS, HLT, DQM] File eXchange Server
	TString fFXSUser[4]; 		// username of the [DAQ, DCS, HLT, DQM] File eXchange Server
	TString fFXSPass[4]; 		// password of the [DAQ, DCS, HLT, DQM] File eXchange Server
	TString fFXSBaseFolder[4]; 	// base folder of the [DAQ, DCS, HLT, DQM] File eXchange Server

	TString fFXSdbHost[4];		// Host of the [DAQ, DCS, HLT, DQM] FXS database
	UInt_t  fFXSdbPort[4];		// Port of the [DAQ, DCS, HLT, DQM] FXS database
	TString fFXSdbUser[4];  	// username of the [DAQ, DCS, HLT, DQM] FXS database
	TString fFXSdbPass[4]; 		// password of the [DAQ, DCS, HLT, DQM] FXS database
	TString fFXSdbName[4]; 		// name of the [DAQ, DCS, HLT, DQM] FXS database
	TString fFXSdbTable[4]; 	// Table name of the [DAQ, DCS, HLT, DQM] FXS database

	TString	fTerminateFilePath;     // File for inter-process communcation to terminate the Shuttle

	Int_t fMaxRetries;        	// number of retries of a failed preprocessor
	
	Int_t fPPTimeOut;         	// timeout until a preprocessor is killed
	Int_t fDCSTimeOut;         	// timeout until the query to DCS is terminated
	Int_t fDCSRetries;         	// number of retries until the query to DCS connection is terminated
	UInt_t fDCSQueryOffset;         // offset for the DCS DPs query, to be subracted from the StartTime and
	                                // added to the EndTime
	UInt_t fDCSDelay;         	// delay to make the query to DCS
	Int_t fPPMaxMem;          	// maximum allowed memory until a preprocessor is killed
	
	TString fMonitorHost;           // host of the MonaLisa monitoring server
	TString fMonitorTable;          // Monalisa's SHUTTLE table name
	
	Int_t fTriggerWait;           	// time to wait for DIM trigger before starting new collection 
	TString fShuttleFileSystem;     // path of the Shuttle file system
	Int_t fFreeDiskWarningThreshold;   // threshold for free space in the Shuttle file system to send a mail to the responsibles
	Int_t fFreeDiskFatalThreshold;     // threshold for free space in the Shuttle file system to send a mail to the responsibles and terminate the Shuttle
	RunMode  fRunMode;          	// Working mode (0=test; 1=prod)

	TMap fDetectorMap; 		// Map of the detector-by-detector configuration
	TObjArray fDetectorList; 	// List of detectors with valid configuration
	
	TObjArray *fAdmin[6]; 		// Array of system administrators' email addresses (DAQ, DCS, HLT, DQM, Global, Amanda)

	TString fShuttleInstanceHost; 	// Instance of the SHUTTLE
	TObjArray fProcessedDetectors; 	// list of the detector to be processed by this machine
	
	Bool_t fKeepDCSMap;		// Flag to keep DCS map also in case of success 
	Bool_t fKeepTempFolder; 	// Flag to keep temp folder also in case of success

	Bool_t fSendMail; 		// Send mail flag
	
	Bool_t fProcessAll; 		// flag indicating that all detectors will be processed
	Bool_t fIsValid;  		// flag for the validity of the configuration

	ClassDef(AliShuttleConfig, 0);
};

#endif

