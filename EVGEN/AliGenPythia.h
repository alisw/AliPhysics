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

  typedef enum {kFlavorSelection, kParentSelection} StackFillOpt_t;
  typedef enum {kCountAll, kCountParents, kCountTrackables} CountMode_t;

    AliGenPythia();
    AliGenPythia(Int_t npart);
    AliGenPythia(const AliGenPythia &Pythia);
    virtual ~AliGenPythia();
    virtual void    Generate();
    virtual void    Init();
    // set a cut on the Z coord. of the primary vertex (cm)
    //
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
    virtual void    SetJetEtaRange(Float_t etamin = -20., Float_t etamax = 20.)
	{fEtaMinJet = etamin; fEtaMaxJet = etamax;}
    virtual void    SetJetPhiRange(Float_t phimin = -180., Float_t phimax = 180.)
	{fPhiMinJet = TMath::Pi()*phimin/180.; fPhiMaxJet = TMath::Pi()*phimax/180.;}
    virtual void    SetGammaEtaRange(Float_t etamin = -20., Float_t etamax = 20.)
	{fEtaMinGamma = etamin; fEtaMaxGamma = etamax;}
    virtual void    SetGammaPhiRange(Float_t phimin = -180., Float_t phimax = 180.)
	{fPhiMinGamma = TMath::Pi()*phimin/180.; fPhiMaxGamma = TMath::Pi()*phimax/180.;}
    // Set option for feed down from higher family
    virtual void SetFeedDownHigherFamily(Bool_t opt) {
      fFeedDownOpt = opt;
    }
    // Set option for selecting particles kept in stack according to flavor
    // or to parent selection
    virtual void SetStackFillOpt(StackFillOpt_t opt) {
      fStackFillOpt = opt;
    }
    // Set fragmentation option
    virtual void SetFragmentation(const Bool_t opt) {
      fFragmentation = opt;
    }
    // Set counting mode
    virtual void SetCountMode(const CountMode_t mode) {
      fCountMode = mode;
    }
    
    // get cross section of process
    virtual Float_t GetXsection() {return fXsection;}      
    virtual void    FinishRun();
    Bool_t CheckTrigger(TParticle* jet1, TParticle* jet2);
    
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
    Float_t     fEtaMinJet;      // Minimum eta of triggered Jet
    Float_t     fEtaMaxJet;      // Maximum eta of triggered Jet
    Float_t     fPhiMinJet;      // Minimum phi of triggered Jet
    Float_t     fPhiMaxJet;      // Maximum phi of triggered Jet

    Float_t     fEtaMinGamma;    // Minimum eta of triggered gamma
    Float_t     fEtaMaxGamma;    // Maximum eta of triggered gamma
    Float_t     fPhiMinGamma;    // Minimum phi of triggered gamma
    Float_t     fPhiMaxGamma;    // Maximum phi of triggered gamma

    StackFillOpt_t fStackFillOpt; // Stack filling with all particles with
                                  // that flavour or only with selected
                                  // parents and their decays
    Bool_t fFeedDownOpt;          // Option to set feed down from higher
                                  // quark families (e.g. b->c)
    Bool_t fFragmentation;        // Option to activate fragmentation by Pythia
    //
    // Options for counting when the event will be finished.
    // fCountMode = kCountAll         --> All particles that end up in the
    //                                    stack are counted
    // fCountMode = kCountParents     --> Only selected parents are counted
    // fCountMode = kCountTrackabless --> Only particles flagged for tracking
    //                                     are counted
    CountMode_t fCountMode;

 private:
    // adjust the weight from kinematic cuts
    void   AdjustWeights();

    ClassDef(AliGenPythia,2) // AliGenerator interface to Pythia
};
#endif





