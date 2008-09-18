//
// Class AliRsnReaderTaskSE
//
// An AnalysisTask object to convert any kind of source event type (ESD/AOD/MC)
// into the RSN internal format (AliRsnEvent).
// The output of this task is a TTree with converted events, which is saved in a file
// and can then be processed as many times as desired, to build invariant mass spectra.
// ---
// original author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
// adapted for Analysis Framework by: R. Vernet (renaud.vernet@cern.ch)
//

#ifndef AliRsnReaderTaskSE_H
#define AliRsnReaderTaskSE_H

#include "AliRsnAnalysisTaskSEBase.h"

class AliRsnPID;
class AliESDEvent;
class AliRsnReader;

class AliRsnReaderTaskSE : public AliRsnAnalysisTaskSEBase
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

//     void          SetReader(AliRsnReader *reader) {fReader = reader;}
//     void          SetPID(AliRsnPID *pid) {fPID = pid;}
//     AliRsnReader* GetReader() {return fReader;}
//     AliRsnPID*    GetPID() {return fPID;}
    AliRsnEvent*  GetCurrentEvent() {return fRsnEvent;}

private:

    AliRsnReaderTaskSE(const AliRsnReaderTaskSE &copy) :
      AliRsnAnalysisTaskSEBase(copy),fRsnEvent(0x0) { /*nothing*/ }
    AliRsnReaderTaskSE& operator=(const AliRsnReaderTaskSE&)
      { /*nothing*/ return (*this); }

    AliRsnEvent  *fRsnEvent;   // output events in the AliRsnEvent format

    ClassDef(AliRsnReaderTaskSE, 1); // implementation of RsnReader as AnalysisTaskSE
};

#endif
