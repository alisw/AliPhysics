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

	const char* GetDAQlbHost() const {return fDAQlbHost.Data();}
	const char* GetDAQlbUser() const {return fDAQlbUser.Data();}
	const char* GetDAQlbPass() const {return fDAQlbPass.Data();}

	const char* GetFESHost(Int_t system) const {return fFESHost[system].Data();}
	const char* GetFESUser(Int_t system) const {return fFESUser[system].Data();}
	const char* GetFESPass(Int_t system) const {return fFESPass[system].Data();}

	const char* GetFESlbHost(Int_t system) const {return fFESlbHost[system].Data();}
	const char* GetFESlbUser(Int_t system) const {return fFESlbUser[system].Data();}
	const char* GetFESlbPass(Int_t system) const {return fFESlbPass[system].Data();}

	Int_t GetMaxRetries() const { return fMaxRetries; }

	Int_t GetPPTimeOut() const { return fPPTimeOut; }

  const TObjArray* GetDetectors() const;

	Bool_t HasDetector(const char* detector) const;
	const char* GetDCSHost(const char* detector) const;
	Int_t GetDCSPort(const char* detector) const;
	const TObjArray* GetDCSAliases(const char* detector) const;
	Bool_t StrictRunOrder(const char* detector) const;

	void SetProcessAll(Bool_t flag=kTRUE) {fProcessAll=flag;}
	Bool_t ProcessAll() const {return fProcessAll;}

	Bool_t HostProcessDetector(const char* detector) const;

	virtual void Print(Option_t* option = NULL) const;

private:

	class AliShuttleConfigHolder: public TObject {
	public:
		AliShuttleConfigHolder(const TLDAPEntry* entry);
		~AliShuttleConfigHolder();

		const char* GetDetector() const {return fDetector.Data();}
		const char* GetDCSHost() const {return fDCSHost.Data();}
		Int_t GetDCSPort() const {return fDCSPort;}
		const TObjArray* GetDCSAliases() const {return fDCSAliases;}

		Bool_t IsValid() const {return fIsValid;}
		Bool_t SkipDCSQuery() const {return fSkipDCSQuery;}
		Bool_t StrictRunOrder() const {return fStrictRunOrder;}

	private:
		TString fDetector;  	// Detector name
		TString fDCSHost; 	// Host name of the DCS server
		Int_t 	fDCSPort; 	// port of the DCS server
		TObjArray* fDCSAliases; // List of DCS aliases to be retrieved
		Bool_t fIsValid;  	// flag for the validity of the configuration
		Bool_t fSkipDCSQuery; 	// flag - if TRUE (-> DCS config empty) skip DCS archive data query
		Bool_t fStrictRunOrder; // flag - if TRUE connect data in a strict run ordering


		ClassDef(AliShuttleConfigHolder, 0);
	};


	Bool_t fIsValid;  		//! flag for the validity of the configuration

	TString fDAQlbHost;		//! Host of the DAQ logbook MySQL Server
	TString fDAQlbUser;  		//! username of the DAQ logbook MySQL Server
	TString fDAQlbPass; 		//! password of the DAQ logbook MySQL Server

	TString fFESHost[3]; 		//! Host of the [DAQ, DCS, HLT] File Exchange Server
	TString fFESUser[3]; 		//! username of the [DAQ, DCS, HLT] File Exchange Server
	TString fFESPass[3]; 		//! password of the [DAQ, DCS, HLT] File Exchange Server

	TString fFESlbHost[3];		//! Host of the [DAQ, DCS, HLT] FES logbook
	TString fFESlbUser[3];  	//! username of the [DAQ, DCS, HLT] FES logbook
	TString fFESlbPass[3]; 		//! password of the [DAQ, DCS, HLT] FES logbook

	Int_t fMaxRetries;        // number of retries of a failed preprocessor

	Int_t fPPTimeOut;         // timeout until a preprocessor is canceled

	TMap fDetectorMap; 		//! Map of the detector-by-detector configuration
	TObjArray fDetectorList; 	//! List of detectors with valid configuration

	TString fShuttleInstanceHost; 	//! Instance of the SHUTTLE
	TObjArray fProcessedDetectors; 	//! list of the detector to be processed by this machine
	Bool_t fProcessAll; 		//! flag indicating that all detectors will be processed

	ClassDef(AliShuttleConfig, 0);
};

#endif

