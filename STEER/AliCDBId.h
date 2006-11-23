#ifndef ALI_CDB_ID_H
#define ALI_CDB_ID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBId						   //
//  Identity of an object stored into a database:  		   //
//  path, run validity range, version, subVersion 		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliCDBPath.h"
#include "AliCDBRunRange.h"

#include <TObject.h>

class AliCDBId: public TObject {

public:

	AliCDBId();

	AliCDBId(const AliCDBId& other);

	AliCDBId(const AliCDBPath& path, const AliCDBRunRange& runRange,
		Int_t version = -1, Int_t subVersion = -1);

	AliCDBId(const AliCDBPath& path, 
		Int_t firstRun , Int_t lastRun , Int_t verison = -1, 
		Int_t subVersion = -1); 

	virtual ~AliCDBId();

	const AliCDBPath& 	GetAliCDBPath() const {return fPath;};
	const TString& 		GetPath() const {return fPath.GetPath();};
	const TString& 		GetLevel0() const {return fPath.GetLevel0();};
	const TString& 		GetLevel1() const {return fPath.GetLevel1();};
	const TString& 		GetLevel2() const {return fPath.GetLevel2();};
	Bool_t 			IsWildcard() const {return fPath.IsWildcard();};

	void 	SetPath(const char* path) {fPath.SetPath(path);};

	const 		AliCDBRunRange& GetAliCDBRunRange() const {return fRunRange;};
	AliCDBRunRange& GetAliCDBRunRange() {return fRunRange;};
	Int_t 		GetFirstRun() const {return fRunRange.GetFirstRun();};
	Int_t 		GetLastRun() const {return fRunRange.GetLastRun();};	
	void 		SetFirstRun(Int_t firstRun) {fRunRange.SetFirstRun(firstRun);};
	void 		SetLastRun(Int_t lastRun) {fRunRange.SetLastRun(lastRun);};
	void 		SetRunRange(Int_t firstRun, Int_t lastRun) 
			{fRunRange.SetRunRange(firstRun, lastRun);};


	Bool_t 	IsAnyRange() const {return fRunRange.IsAnyRange();};


	Int_t 	GetVersion() const {return fVersion;};
	Int_t 	GetSubVersion() const {return fSubVersion;};
	void 	SetVersion(Int_t version) {fVersion = version;};	
	void 	SetSubVersion(Int_t subVersion) {fSubVersion = subVersion;};

	const TString& 	GetLastStorage() const {return fLastStorage;};
	void 		SetLastStorage(TString& lastStorage){fLastStorage = lastStorage;};

	Bool_t IsValid() const; 
	Bool_t IsSpecified() const {return !(IsWildcard() || IsAnyRange());};

	Bool_t HasVersion() const {return fVersion >= 0;};
	Bool_t HasSubVersion() const {return fSubVersion >= 0;};

	Bool_t Comprises(const AliCDBId& other) const
		{return fPath.Comprises(other.fPath)
		         && fRunRange.Comprises(other.fRunRange);};

	virtual Bool_t IsEqual(const TObject *obj) const;

	TString ToString() const;

private:

	AliCDBPath fPath;		// path	
	AliCDBRunRange fRunRange;	// run range
	Int_t fVersion;			// version
	Int_t fSubVersion;		// subversion
	TString fLastStorage;		// previous storage place (new, grid, local, dump)

	ClassDef(AliCDBId, 1);
};

#endif
