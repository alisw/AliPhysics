// Author: Mihai Niculescu 2012

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, all rights reserved. 					 *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          										 *
 * full copyright notice.                                                																						 *
 **************************************************************************/
 
#ifndef AliEveApplication_H
#define AliEveApplication_H

#include <TRint.h>

class AliEveApplication : public TRint
{
public:
	AliEveApplication(const char* appClassName, Int_t* argc, char** argv, void* options = 0, Int_t numOptions = 0, Bool_t noLogo = kFALSE);
	virtual ~AliEveApplication();
	
	void Init(); // Initialize AliEve & Rint Environment
	
private:
	AliEveApplication(const AliEveApplication&);               // not implemented
	AliEveApplication& operator=(const AliEveApplication&);    // not implemented
	
public:

	ClassDef(AliEveApplication, 0); // AliEve application
};
#endif
