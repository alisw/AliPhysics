#ifndef ALIGENHERWIG_H
#define ALIGENHERWIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Generator using HERWIG as an external generator
// The main HERWIG options are accessable for the user through this interface.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include "AliGenMC.h"
#include <TString.h>
#include <TArrayI.h>
#include <AliRndm.h>
#include <AliStructFuncType.h>

class THerwig6;
class TArrayI;
class TParticle;
class TClonesArray;


class AliGenHerwig : public AliGenMC

{
    enum {kNoTrigger, kHardProcesses, kDirectPhotons};

 public:
    AliGenHerwig();
    AliGenHerwig(Int_t npart);
    virtual ~AliGenHerwig();
    virtual void    Generate();
    virtual void    Init();
    virtual void    InitJimmy();
    // set centre of mass energy
    virtual void    SetBeamMomenta(Float_t p1=7000., Float_t p2=7000.)
	{fMomentum1 = p1; fMomentum2 = p2;}
    virtual void    SetProcess(Int_t proc)            {fProcess = proc;}
    virtual void    KeepFullEvent();
    virtual void    SetDecaysOff(Int_t flag=1)        {fDecaysOff = flag;}
    virtual void    SetTrigger(Int_t flag=kNoTrigger) {fTrigger   = flag;}
    virtual void    SetFlavor(Int_t flag=0)           {fFlavor    = flag;}
    virtual void    SetSelectAll(Int_t flag=0)        {fSelectAll = flag;}
    virtual void    SetStrucFunc(StrucFunc_t func = kCTEQ5L)
      {fStrucFunc = func;}
    virtual void    SetPtHardMin(Double_t pt) {fPtHardMin=pt;}
    virtual void    SetPtRMS(Double_t pt) {fPtRMS=pt;}
    virtual void    SetMaxPr(Int_t i) {fMaxPr=i;}
    virtual void    SetMaxErrors(Int_t i) {fMaxErrors=i;}
    virtual void    FinishRun();
    virtual void    FinishRunJimmy();
    virtual void    SetEnSoft(Double_t e) {fEnSoft=e;}

    virtual void    SetHardProcessFile(TString filename) {fFileName=TString(filename);};
    virtual void    SetEventListRange(Int_t eventFirst=-1, Int_t eventLast=-1);

 protected:
    Bool_t SelectFlavor(Int_t pid);

 protected:
    TString     fAutPDF;         // PDF group
    Int_t       fModPDF;         // PDF set
    StrucFunc_t fStrucFunc;      //Structure Function
    Int_t       fKeep;           // Flag to keep full event information
    Int_t       fDecaysOff;      // Flag to turn off decays of pi0, K_s, D, Lambda, sigma
    Int_t       fTrigger;        // Trigger type
    Int_t       fSelectAll;      // Flag to write the full event
    Int_t       fFlavor;         // Selected particle flavor 4: charm+beauty 5: beauty
    Float_t     fEnergyCMS;      // Centre of mass energy
    Float_t     fMomentum1;      // Momentum of projectile
    Float_t     fMomentum2;      // Momentum of target
    Float_t     fKineBias;       // Bias from kinematic selection
    Int_t       fTrials;         // Number of trials
    Float_t     fXsection;       // Cross-section
    THerwig6    *fHerwig;        // Herwig
    Int_t       fProcess;        // Process number
    Double_t    fPtHardMin;      // lower pT-hard cut
    Double_t    fPtRMS;          // intrinsic pt of incoming hadrons
    Int_t       fMaxPr;          // maximum number of events to print out
    Int_t       fMaxErrors;      // maximum number of errors allowed
    Double_t    fEnSoft;          // change on soft energy distribution
    Int_t       fEv1Pr;          // first event to be printed
    Int_t       fEv2Pr;          // last event to be printed
    TString     fFileName;       //!Name of file to read from hard scattering

 private:
    AliGenHerwig(const AliGenHerwig &Herwig);
    AliGenHerwig &  operator=(const AliGenHerwig & rhs);

    // check if particle is selected as parent particle
    Bool_t ParentSelected(Int_t ip);
    // check if particle is selected as child particle
    Bool_t ChildSelected(Int_t ip);
    // adjust the weight from kinematic cuts
    void   AdjustWeights();
    // check seleted daughters
    Bool_t DaughtersSelection(TParticle* iparticle, TClonesArray* particles);
    // check if stable
    Bool_t Stable(TParticle*  particle);

    void InitPDF();

    ClassDef(AliGenHerwig,1) // AliGenerator interface to Herwig
};
#endif





