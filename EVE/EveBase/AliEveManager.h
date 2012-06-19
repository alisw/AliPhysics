// Author: Mihai Niculescu 2012

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, all rights reserved. 					 *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          										 *
 * full copyright notice.                                                																						 *
 **************************************************************************/
 
#ifndef AliEveManager_H
#define AliEveManager_H

#include <TEveManager.h>

class AliEveManager : public TEveManager
{
public:
	AliEveManager(UInt_t w, UInt_t h, Bool_t map_window=kTRUE, Option_t* opt="FI");
	~AliEveManager();

	static AliEveManager* Create(Bool_t map_window=kTRUE, Option_t* opt="FIV");
	
	void CloseEveWindow();
	void Terminate();
protected:
	void Init();

public:

	ClassDef(AliEveManager, 0); // Eve application manager
};
#endif
