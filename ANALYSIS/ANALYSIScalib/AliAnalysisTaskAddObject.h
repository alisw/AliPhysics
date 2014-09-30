#ifndef ALIANALYSISTASKADDOBJECT_H
#define ALIANALYSISTASKADDOBJECT_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

//*************************************************************************
// Class AliAnalysisTaskAddObject
// Test Task to add an object to the new ESDfriends file 
//*************************************************************************

class TH1D;

#include "AliAnalysisTask.h"

class AliESDInputHandler;
class AliESDEvent;
class AliESDfriend;

class AliAnalysisTaskAddObject : public AliAnalysisTask
{
 public:

	AliAnalysisTaskAddObject();
	AliAnalysisTaskAddObject(const char *name);
	virtual ~AliAnalysisTaskAddObject();
	// Implementation of interface methods
	virtual void CreateOutputObjects();
	virtual void Exec(Option_t *option);
	virtual void Terminate(Option_t *option);
	virtual void ConnectInputData(Option_t *option = "");
		
 private:
	
	AliAnalysisTaskAddObject(const AliAnalysisTaskAddObject &);
	AliAnalysisTaskAddObject& operator=(const AliAnalysisTaskAddObject&);
	
	AliESDEvent  *fESDInput;        // ESD input object
	AliESDfriend *fESDfriendInput;  // ESD input friend object
	AliESDInputHandler *fESDhandler;     // Pointer to ESD input handler
	TH1D* fh; // histogram
	
 ClassDef(AliAnalysisTaskAddObject,1); // AliAnalysisTask to create an extra object
};

#endif

