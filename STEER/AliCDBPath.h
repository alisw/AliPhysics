#ifndef ALI_CDB_PATH_H
#define ALI_CDB_PATH_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBPath						   //
//  Path string identifying the object:  			   //
//  "level0/level1/level2" 					   //
//  (was: "Detector/DBType/DetSpecType") 		           //
//  (example: "ZDC/Calib/Pedestals") 		         	   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

class AliCDBPath: public TObject {

public:
	
	AliCDBPath();

	AliCDBPath(const AliCDBPath& other);

	AliCDBPath(const char* level0, const char* level1, 
			const char* level2);

	AliCDBPath(const char* path);

	AliCDBPath(const TString& path);

	virtual ~AliCDBPath();


	const TString& GetPath() const {return fPath;}
	void  SetPath(const char* path) {fPath=path; InitPath();}

	const char* GetLevel(Int_t i) const;

	Bool_t IsValid() const {return fIsValid;}

	Bool_t IsWildcard() const {return fIsWildcard;}

	Bool_t Level0Comprises(const TString& str) const;
	Bool_t Level1Comprises(const TString& str) const;
	Bool_t Level2Comprises(const TString& str) const;

	Bool_t Comprises(const AliCDBPath& other) const;

private:

	Bool_t IsWord(const TString& str);

	void InitPath();

	void Init();

	TString fPath;		// detector pathname (Detector/DBType/SpecType)
	TString fLevel0;	// level0 name (ex. detector: ZDC, TPC...)
	TString fLevel1;	// level1 name (ex. DB type, Calib, Align)
	TString fLevel2;	// level2 name (ex. DetSpecType, pedestals, gain...)

	Bool_t fIsValid;	// validity flag
	Bool_t fIsWildcard;	// wildcard flag
	
	ClassDef(AliCDBPath, 1); 
};

#endif
