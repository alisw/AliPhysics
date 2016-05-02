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
    Int_t    ProcessType() const {return fProcessType;}
    void     SetNMPI(Int_t mpi) {fNMPI = mpi;}
    void     SetNSuperpositions(Int_t sp) { fNSuperpositions = sp; }
    // Setters
    void     SetProcessType(Int_t type)  {fProcessType = type;}
    Int_t    Trials() const {return fTrials;}
    void     SetTrials(Int_t trials) {fTrials = trials;}
    void     AddJet(Float_t px, Float_t py, Float_t pz, Float_t e);
    void     AddUQJet(Float_t px, Float_t py, Float_t pz, Float_t e);
    Int_t    NTriggerJets() const {return fNJets;}
    Int_t    NUQTriggerJets() const {return fNUQJets;}
    void     TriggerJet(Int_t i, Float_t p[4]) const;
    void     UQJet(Int_t i, Float_t p[4]) const;
    Double_t GetXJet() const {return fXJet;}
    Double_t GetYJet() const {return fYJet;}
    Double_t GetInMediumLength() const  {return fInMediumLength;}
    Double_t GetImpactParameter() const {return fImpactParameter;}
    void     SetXYJet(Double_t x, Double_t y);
    void     SetImpactParameter(Double_t b) {fImpactParameter = b;}
    void     SetInMe(Double_t l) {fInMediumLength = l;}
    void     SetZQuench(Double_t z[4]);
    void     GetZQuench(Double_t z[4]) const;
    void     SetPtHard(Float_t pthard) {fPtHard = pthard;}
    Float_t  GetPtHard() const {return fPtHard;}    
    void     SetXsection(Float_t xsec) {fXsection = xsec;}
    Float_t  GetXsection() const {return fXsection;}
    Int_t    GetNMPI() const {return fNMPI;}
    Int_t    GetNSuperpositions() const { return fNSuperpositions; }
	    
protected:
    Int_t    fProcessType;               // PYTHIA process id for this event 
    Int_t    fTrials;                    // Number of trials to fulfill trigger condition
    Int_t    fNJets;                     // Number of triggered jets
    Int_t    fNUQJets;                   // Number of unquenched jets
    Int_t    fNMPI;                      // numbers of MPI
    Int_t    fNSuperpositions;           // numbers of superimposed events
    Double_t fXJet;                      // Jet production point (x)
    Double_t fYJet;                      // Jet production point (y)
    Double_t fInMediumLength;            // In medium length
    Double_t fImpactParameter;           // Impact parameter for Q-Pythia
    Float_t  fJets[4][10];               // Trigger jets
    Float_t  fUQJets[4][10];             // Unquenched trigger jets
    Double_t fZquench[4];                // Quenching fraction
    Float_t  fPtHard;                    // pT hard
    Float_t  fXsection;                  // Cross-section

    ClassDef(AliGenPythiaEventHeader,9)  // Event header for Pythia event
};
	
	

#endif
