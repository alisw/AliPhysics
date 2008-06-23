/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//=========================================================================
// Class AliRsnAnalyzerSimple
//
// Implementation of the event processing which returns histograms of
// invariant mass for resonances and backgrounds evaluation.
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//=========================================================================

#ifndef AliRsnAnalyzerSimple_H
#define AliRsnAnalyzerSimple_H

#include "AliRsnPID.h"

class TObjArray;
class AliRsnPairSimple;
class AliRsnEventBuffer;

class AliRsnAnalyzerSimple : public TObject
{
public:

    AliRsnAnalyzerSimple(Int_t bufferSize = 1000);
    virtual ~AliRsnAnalyzerSimple() {Clear();}

    // setters
    void    SetBufferSize(Int_t value) {fBufferSize = value;}
    void    SetMixMultiplicityCut(Int_t cut) {fMixMultCut = cut;}
    void    SetMixVzCut(Double_t cut) {fMixVzCut = cut;}
    void    SetNMix(Int_t value) {fNMix = value;}

    // getters
    TObjArray* GetPairs() {return fPairs;}
    TObjArray* GetMixPairs() {return fMixPairs;}

    // working routines
    void         Init();
    virtual void Clear(Option_t *option = "C");
    void         AddPair(AliRsnPairSimple *pair);
    Stat_t       Process(AliRsnEvent *event);

private:

    AliRsnAnalyzerSimple(const AliRsnAnalyzerSimple &copy) :
      TObject(copy),fBufferSize(1000),
      fMixMultCut(10),fMixVzCut(0.5),fNMix(10),
      fPairs(0x0),fMixPairs(0x0),fBuffer(0x0) { }
    AliRsnAnalyzerSimple& operator=(const AliRsnAnalyzerSimple & /*copy*/) { return (*this); }

    Int_t         fBufferSize;       // size of buffer
    Int_t         fMixMultCut;       // multiplicity cut for event mixing
    Double_t      fMixVzCut;         // difference in VZ cut for event mixing
    Int_t         fNMix;             // number of events for mixing

    TObjArray         *fPairs;       //! collection of particle pairs to be read for 1-event analysis
    TObjArray         *fMixPairs;    //! collection of particle pairs to be read for event-mixing
    AliRsnEventBuffer *fBuffer;      //! buffer for event mixing

    // Rsn analyzer implementation
    ClassDef(AliRsnAnalyzerSimple,1)
};

#endif
