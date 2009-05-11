#ifndef ALIGENPYTHIAEVENTHEADER_H
#define ALIGENPYTHIAEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenPythiaEventHeader : public AliGenEventHeader
{
 public:
    AliGenPythiaEventHeader();
    AliGenPythiaEventHeader(const char* name);
    virtual ~AliGenPythiaEventHeader() {}
    // Getters
    Int_t    ProcessType()  {return fProcessType;}
    // Setters
    void     SetProcessType(Int_t type)  {fProcessType = type;}
    Int_t    Trials() {return fTrials;}
    void     SetTrials(Int_t trials) {fTrials = trials;}
    void     AddJet(Float_t px, Float_t py, Float_t pz, Float_t e);
    void     AddUQJet(Float_t px, Float_t py, Float_t pz, Float_t e);
    Int_t    NTriggerJets() {return fNJets;}
    Int_t    NUQTriggerJets() {return fNUQJets;}
    void     TriggerJet(Int_t i, Float_t p[4]);
    void     UQJet(Int_t i, Float_t p[4]);
    Double_t GetXJet() {return fXJet;}
    Double_t GetYJet() {return fYJet;}
    Double_t GetInMediumLength() {return fInMediumLength;}
    void     SetXYJet(Double_t x, Double_t y);
    void     SetInMediumLength(Double_t l) {fInMediumLength = l;}
    void     SetZQuench(Double_t z[4]);
    void     GetZQuench(Double_t z[4]);
    void     SetPtHard(Float_t pthard) {fPtHard = pthard;}
    Float_t  GetPtHard() {return fPtHard;}    
	
	    
protected:
    Int_t    fProcessType;               // PYTHIA process id for this event 
    Int_t    fTrials;                    // Number of trials to fulfill trigger condition
    Int_t    fNJets;                     // Number of triggered jets
    Int_t    fNUQJets;                   // Number of unquenched
    Double_t fXJet;                      // Jet production point (x)
    Double_t fYJet;                      // Jet production point (y)
    Double_t fInMediumLength;            // In medium length
    Float_t  fJets[4][10];               // Trigger jets
    Float_t  fUQJets[4][10];             // Unquenched trigger jets
    Double_t fZquench[4];                // Quenching fraction
    Float_t  fPtHard;                    // pT hard
    ClassDef(AliGenPythiaEventHeader,5)  // Event header for Pythia event
};
	
	

#endif
