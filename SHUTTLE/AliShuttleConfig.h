#ifndef ALI_SHUTTLE_CONFIG_H
#define ALI_SHUTTLE_CONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class keeps the AliShuttle configuration.
// It reads the configuration for LDAP server.
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
	Bool_t ProcessAll() {return fProcessAll;}

	Bool_t HostProcessDetector(const char* detector) const;

	virtual void Print(Option_t* option = NULL) const;

private:

	class ConfigHolder: public TObject {
		TString fDetector;
		TString fDCSHost;
		Int_t 	fDCSPort;
		TObjArray fDCSAliases;
		TObjArray fDAQFileIDs;
		Bool_t fIsValid;

	public:
		ConfigHolder(const TLDAPEntry* entry);
		~ConfigHolder();

		const char* GetDetector() const {return fDetector.Data();};
		const char* GetDCSHost() const {return fDCSHost.Data();};
		Int_t GetDCSPort() const {return fDCSPort;};
		const TObjArray* GetDCSAliases() const {return &fDCSAliases;};
		const TObjArray* GetDAQFileIDs() const {return &fDAQFileIDs;};

		Bool_t IsValid() const {return fIsValid;};

		ClassDef(ConfigHolder, 0);
	};


	Bool_t fIsValid;

	TString fDAQLogBookHost;
	TString fDAQLogBookUser;
	TString fDAQLogBookPassword;

	TString fDAQFSHost;

	TMap fDetectorMap;
	TObjArray fDetectorList;

	TString fShuttleInstanceHost;
	TObjArray fProcessedDetectors;
	Bool_t fProcessAll;

	ClassDef(AliShuttleConfig, 0);
};

#endif

