#ifndef ALIGENPYTHIA_H
#define ALIGENPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliGenMC.h"
#include "AliPythia.h"
#include "TArrayF.h"

class AliPythia;
class TParticle;

class AliGenPythia : public AliGenMC
{
 public:
    AliGenPythia();
    AliGenPythia(Int_t npart);
    AliGenPythia(const AliGenPythia &Pythia);
    virtual ~AliGenPythia();
    virtual void    Generate();
    virtual void    Init();
    virtual void    SetEventListRange(Int_t eventFirst=-1, Int_t eventLast=-1);
    // select process type
    virtual void    SetProcess(Process_t proc = kPyCharm) {fProcess = proc;}
    // select structure function
    virtual void    SetStrucFunc(StrucFunc_t func = kGRV_HO) {fStrucFunc = func;}
    // select pt of hard scattering 
    virtual void    SetPtHard(Float_t ptmin = 0, Float_t ptmax = 1.e10)
	{fPtHardMin = ptmin; fPtHardMax = ptmax; }
    // set centre of mass energy
    virtual void    SetEnergyCMS(Float_t energy = 5500) {fEnergyCMS = energy;}
    // treat protons as inside nuclei
    virtual void    SetNuclei(Int_t a1, Int_t a2);
    // get cross section of process
    virtual Float_t GetXsection() {return fXsection;}      
    virtual void    FinishRun();
    
    // Assignment Operator
    AliGenPythia & operator=(const AliGenPythia & rhs);
 private:
    Int_t  GenerateMB();
    virtual void    MakeHeader();    
 protected:
    TClonesArray* fParticles;     //Particle  List
    
    Process_t   fProcess;         //Process type
    StrucFunc_t fStrucFunc;       //Structure Function
    Float_t     fEnergyCMS;       //Centre of mass energy
    Float_t     fKineBias;        //!Bias from kinematic selection
    Int_t       fTrials;          //!Number of trials
    Int_t       fFlavorSelect;    //Heavy Flavor Selection
    Float_t     fXsection;        //Cross-section
    AliPythia   *fPythia;         //!Pythia 
    Float_t     fPtHardMin;       //lower pT-hard cut 
    Float_t     fPtHardMax;       //higher pT-hard cut
    Int_t       fNucA1;           //mass number nucleus side 1
    Int_t       fNucA2;           //mass number nucleus side 2
    Bool_t      fFullEvent;       //!Write Full event if true
    AliDecayer  *fDecayer;        //!Pointer to the decayer instance
    Int_t       fDebugEventFirst; //!First event to debug
    Int_t       fDebugEventLast;  //!Last  event to debug
    TArrayF     fEventVertex;     //!The current event vertex
    
 private:
    // adjust the weight from kinematic cuts
    void   AdjustWeights();

    ClassDef(AliGenPythia,1) // AliGenerator interface to Pythia
};
#endif





