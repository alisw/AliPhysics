#ifndef ALI_META_DATA_H
#define ALI_META_DATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBMetaData					   //
//  Set of data describing the object  				   //
//  but not used to identify the object 			   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMap.h>

class AliCDBMetaData: public TObject {
	
public:
	AliCDBMetaData();
	AliCDBMetaData(const char *responsible, UInt_t beamPeriod=0, const char* alirootVersion="", const char* comment="");
	virtual ~AliCDBMetaData();

	void 		SetObjectClassName(const char* name) 
				{fObjectClassName = name;};
	const char* 	GetObjectClassName() const 
				{return fObjectClassName.Data();};
	
        void 		SetResponsible(const char* yourName) 
				{fResponsible = yourName;};
	const char* 	GetResponsible() const 
				{return fResponsible.Data();};

	void 		SetBeamPeriod(UInt_t period) 
				{fBeamPeriod = period;};
	UInt_t 		GetBeamPeriod() const 
				{return fBeamPeriod;};

	void 		SetAliRootVersion(const char* version)
				{fAliRootVersion = version;};
	const char* 	GetAliRootVersion() const 
				{return fAliRootVersion.Data();};

	void 		SetComment(const char* comment) 
				{fComment = comment;};
	const char* 	GetComment() const 
				{return fComment.Data();};

	void 		SetProperty(const char* property, TObject* object);
	TObject* 	GetProperty(const char* property) const;
	Bool_t 		RemoveProperty(const char* property);
	
	void PrintMetaData();
	
private:

	TString fObjectClassName; 	// object's class name
	TString fResponsible; 		// object's responsible person
	UInt_t  fBeamPeriod;		// beam period
	TString fAliRootVersion;	// AliRoot version
	TString fComment;		// extra comments
	//TList fCalibRuns;		
	
	TMap fProperties;		// list of object specific properties

	ClassDef(AliCDBMetaData, 1);
};

#endif
