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
	AliShuttleConfig(const char* host, Int_t port = LDAP_PORT,
			const char* binddn = 0, const char* password = 0,
			const char* basedn = "o=alice,dc=cern,dc=ch");
	virtual ~AliShuttleConfig();

	Bool_t IsValid() const {return fIsValid;};

	const char* GetDAQLogBookHost() const {return fDAQLogBookHost.Data();}
	const char* GetDAQLogBookUser() const {return fDAQLogBookUser.Data();}
	const char* GetDAQLogBookPassword() const {return fDAQLogBookPassword.Data();}

	const char* GetDAQFSHost() const {return fDAQFSHost.Data();}

	const TObjArray* GetDetectors() const;

	Bool_t HasDetector(const char* detector) const;
	const char* GetDCSHost(const char* detector) const;
	Int_t GetDCSPort(const char* detector) const;
	const TObjArray* GetDCSAliases(const char* detector) const;
	const TObjArray* GetDAQFileIDs(const char* detector) const;

	void SetProcessAll(Bool_t flag=kTRUE) {fProcessAll=flag;}
	Bool_t ProcessAll() const {return fProcessAll;}

	Bool_t HostProcessDetector(const char* detector) const;

	virtual void Print(Option_t* option = NULL) const;

private:

	class AliShuttleConfigHolder: public TObject {
	public:
		AliShuttleConfigHolder(const TLDAPEntry* entry);
		~AliShuttleConfigHolder();

		const char* GetDetector() const {return fDetector.Data();};
		const char* GetDCSHost() const {return fDCSHost.Data();};
		Int_t GetDCSPort() const {return fDCSPort;};
		const TObjArray* GetDCSAliases() const {return &fDCSAliases;};
		const TObjArray* GetDAQFileIDs() const {return &fDAQFileIDs;};

		Bool_t IsValid() const {return fIsValid;};

	private:
		TString fDetector;  	// Detector name
		TString fDCSHost; 	// Host name of the DCS server
		Int_t 	fDCSPort; 	// port of the DCS server
		TObjArray fDCSAliases; 	// List of DCS aliases to be retrieved
		TObjArray fDAQFileIDs; 	// list of IDs of the files to be retrived from DAQ
		Bool_t fIsValid;  	// flag for the validity of the configuration


		ClassDef(AliShuttleConfigHolder, 0);
	};


	Bool_t fIsValid;  		// flag for the validity of the configuration

	TString fDAQLogBookHost;	// Host of the DAQ logbook MySQL Server
	TString fDAQLogBookUser;  	// username of the DAQ logbook MySQL Server
	TString fDAQLogBookPassword; 	// password of the DAQ logbook MySQL Server

	TString fDAQFSHost; 		// Host of the DAQ file system

	TMap fDetectorMap; 		// Map of the detector-by-detector configuration
	TObjArray fDetectorList; 	// List of detectors with valid configuration

	TString fShuttleInstanceHost; 	// Instance of the SHUTTLE
	TObjArray fProcessedDetectors; 	// list of the detector to be processed by this machine
	Bool_t fProcessAll; 		// flag indicating that all detectors will be processed

	ClassDef(AliShuttleConfig, 0);
};

#endif

