/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//----------------------------------------------------------------------------------
//  Class AliRsnReaderTask
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// original author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for Analysis Framework
// by    : R. Vernet                          (email: renaud.vernet@cern.ch)
//----------------------------------------------------------------------------------

#ifndef ALIRSNREADERTASK_H
#define ALIRSNREADERTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliRsnReader.h"
#include "AliRsnPID.h"

class TTree;
class TBranch;
class AliESDEvent;

class AliRsnReaderTask : public AliAnalysisTaskSE
{
public:

    enum ESource {
        kESD = 0,
        kAOD,
        kMC
    };

    AliRsnReaderTask(ESource source = kESD);
	AliRsnReaderTask(const char *name, ESource source = kESD);
	virtual ~AliRsnReaderTask() {Clear();}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
	
	// setters
	void SetReader(AliRsnReader *reader) {fReader = reader;}
	void SetPID(AliRsnPID *pid) {fPID = pid;}
	void SetSource(ESource source) {fSource = source;}
	
	// getters
	AliRsnReader* GetReader() {return fReader;}
	AliRsnPID*    GetPID() {return fPID;}

private:

    AliRsnReaderTask(const AliRsnReaderTask&);
	AliRsnReaderTask& operator=(const AliRsnReaderTask&);
	
	ESource fSource;   // source of data
	
	AliRsnReader* fReader;     // read manager
	AliRsnPID*    fPID;        // particle identification manager
	TClonesArray* fRsnEvents;  // output events in the AliRsnEvent format
	
	ClassDef(AliRsnReaderTask, 0); // implementation of RsnReader as AnalysisTaskSE
};

#endif
