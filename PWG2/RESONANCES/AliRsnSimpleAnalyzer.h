/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//=========================================================================
// Class AliRsnSimpleAnalyzer
//
// Implementation of the event processing which returns histograms of
// invariant mass for resonances and backgrounds evaluation.
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//=========================================================================

#ifndef AliRsnSimpleAnalyzer_H
#define AliRsnSimpleAnalyzer_H

#include <TNamed.h>

class TObjArray;
class AliRsnSimpleFunction;
class AliRsnEventBuffer;

class AliRsnSimpleAnalyzer : public TNamed
{
  public:

    AliRsnSimpleAnalyzer(Int_t bufferSize = 1000);
    virtual ~AliRsnSimpleAnalyzer() {Clear();}

    // setters
    void    SetBufferSize(Int_t value) {fBufferSize = value;}
    void    SetMixMultiplicityCut(Int_t cut) {fMixMultCut = cut;}
    void    SetMixVzCut(Double_t cut) {fMixVzCut = cut;}
    void    SetNMix(Int_t value) {fNMix = value;}

    // getters
    TObjArray*   GetSingle() {return fSingle;}
    TObjArray*   GetMix()    {return fMix;}

    // working routines
    void         Init();
    virtual void Clear(Option_t *option = "C");
    void         Add(AliRsnSimpleFunction *pair);
    void         Process(AliRsnEvent *event);

  private:

    AliRsnSimpleAnalyzer(const AliRsnSimpleAnalyzer &copy) :
        TNamed(copy),fBufferSize(1000),
        fMixMultCut(10),fMixVzCut(0.5),fNMix(10),
        fSingle(0x0),fMix(0x0),fBuffer(0x0) { }
    AliRsnSimpleAnalyzer& operator=(const AliRsnSimpleAnalyzer & /*copy*/) { return (*this); }
    void ProcessEvents(TObjArray *pairs, AliRsnEvent *event1, AliRsnEvent *event2 = 0x0);

    Int_t              fBufferSize;  // size of buffer
    Int_t              fMixMultCut;  // multiplicity cut for event mixing
    Double_t           fMixVzCut;    // difference in VZ cut for event mixing
    Int_t              fNMix;        // number of events for mixing

    TObjArray         *fSingle;      // collection of particle pairs to be read for 1-event analysis
    TObjArray         *fMix;         // collection of particle pairs to be read for event-mixing
    AliRsnEventBuffer *fBuffer;      // buffer for event mixing

    // Rsn analyzer implementation
    ClassDef(AliRsnSimpleAnalyzer,1)
};

#endif
