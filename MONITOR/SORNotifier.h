/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/
#ifndef SORNOTIFIER_H
#define SORNOTIFIER_H

//______________________________________________________________________________
//
// ECS Start-of-Run notifier
//
// This class "listens" to the SOR coming from the ECS.
//

// DIM
#include <dic.hxx>

class AliOnlineRecoTrigger;

class SORNotifier: public DimUpdatedInfo
{
public:
        SORNotifier(AliOnlineRecoTrigger* trigger): 
	  DimUpdatedInfo("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS", -1), fRun(-1), fTrigger(trigger) {}

	void infoHandler();

	void errorHandler(int severity, int code, char *msg);

	int GetRun() const {return fRun;}
private:
	SORNotifier(const SORNotifier& other);
	SORNotifier& operator = (const SORNotifier& other);

	int fRun;

	AliOnlineRecoTrigger* fTrigger;
};

#endif
