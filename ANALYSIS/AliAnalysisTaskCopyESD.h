#ifndef ALIANALYSISTASKCOPYESD_H
#define ALIANALYSISTASKCOPYESD_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTaskFilter.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"

class AliAnalysisTaskCopyESD : public AliAnalysisTaskFilter
{
 public:
	AliAnalysisTaskCopyESD();
	AliAnalysisTaskCopyESD(const char* name);
	virtual ~AliAnalysisTaskCopyESD() {;}
	// Implementation of interface methods
	virtual void   UserCreateOutputObjects();
	virtual void   Init();
	virtual void   LocalInit() {Init();}
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *option);
	
	
 private:
	AliAnalysisTaskCopyESD(const AliAnalysisTaskCopyESD&);
	AliAnalysisTaskCopyESD& operator=(const AliAnalysisTaskCopyESD&);
	
	AliESDEvent* fESDEvent;
	AliESDfriend* fESDfriend;

	
	ClassDef(AliAnalysisTaskCopyESD, 1);
};

#endif
