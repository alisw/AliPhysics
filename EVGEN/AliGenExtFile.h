#ifndef ALIGENEXTFILE_H
#define ALIGENEXTFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Event generator that can read events from a files.
// The reading is performed by a realisation of AliGenReader specific to the file format.
// Author: andreas.morsch@cern.ch

///
#include <functional>
///

#include "AliGenMC.h"
class AliGenReader;


class TTree;

class AliGenExtFile : public AliGenMC
{
 public:
    AliGenExtFile();
    AliGenExtFile(Int_t npart);
     virtual ~AliGenExtFile();
    // Initialise
    virtual void Init();
    // generate event
    virtual void Generate();
    void SetReader(AliGenReader* reader) {fReader = reader;}
    void SetStartEvent(Int_t startEvent) {fStartEvent = startEvent;}
    AliGenReader* Reader() const {return fReader;}

    ///
    virtual void SetMultiplicityTrigger(Double_t multCut) { fSetMultTrig = kTRUE; fMultCut = multCut; }               // Enable base multiplicity trigger and set the cut
    virtual void SetPtTrigger(Double_t ptCut) { fSetPtTrig = kTRUE; fPtCut = ptCut; }                                 // Enable base pT trigger and set the cut
    virtual void SetUserTrigger(std::function<Bool_t(AliStack*)> fn) { fSetUserTrig = kTRUE; fUserTrigger = fn; }     // Enable user trigger
    ///

 protected:
    void CdEventFile();
    const Text_t     *fFileName;      //! File to read from
    AliGenReader     *fReader;        //! Reader to read the file
    Int_t  fStartEvent; //! Start event number

    ///
    // Base Multiplicity and pT triggers
    Int_t  fLimitDiscardedEvents;                                                             // Limit on the events discarded consecutively by the trigger
    Bool_t fSetMultTrig;                                                                      // TRUE if the multiplicity trigger is set
    Bool_t fSetPtTrig;                                                                        // TRUE if the pT trigger is set
    Bool_t fSetUserTrig;

    Bool_t MultiplicityTrigger(AliStack *stack);                                              // Base multiplicity trigger method
    Bool_t PtTrigger(AliStack *stack);                                                        // Base pT trigger method

    Int_t fMultCut;
    Double_t fPtCut;

    // Custom trigger
    std::function<Bool_t(AliStack*)> fUserTrigger;
    ///

 private:
    AliGenExtFile(const AliGenExtFile &ext);
    AliGenExtFile & operator=(const AliGenExtFile & rhs);

  ClassDef(AliGenExtFile,2) //Generate particles from external file
};
#endif
