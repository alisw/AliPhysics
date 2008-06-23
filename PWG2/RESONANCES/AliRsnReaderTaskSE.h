/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//----------------------------------------------------------------------------------
//  Class AliRsnReaderTaskSE
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#ifndef AliRsnReaderTaskSE_H
#define AliRsnReaderTaskSE_H

#include "AliAnalysisTaskSE.h"
//#include "AliRsnReader.h"
//#include "AliRsnPID.h"

class AliESDEvent;
class AliRsnPID;
class AliRsnReader;

class AliRsnReaderTaskSE : public AliAnalysisTaskSE
{
public:

    AliRsnReaderTaskSE();
	AliRsnReaderTaskSE(const char *name);
	virtual ~AliRsnReaderTaskSE() {Clear();}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

	// setters
	void SetReader(AliRsnReader *reader) {fReader = reader;}
	void SetPID(AliRsnPID *pid) {fPID = pid;}

	// getters
	AliRsnReader* GetReader() {return fReader;}
	AliRsnPID*    GetPID() {return fPID;}

private:

    AliRsnReaderTaskSE(const AliRsnReaderTaskSE&);
	AliRsnReaderTaskSE& operator=(const AliRsnReaderTaskSE&);

	AliRsnReader *fReader;     // read manager
	AliRsnPID    *fPID;        // particle identification manager
	AliRsnEvent  *fRsnEvent;   // output events in the AliRsnEvent format

	ClassDef(AliRsnReaderTaskSE, 0); // implementation of RsnReader as AnalysisTaskSE
};

#endif
