#ifndef ALIGENPYTHIA_H
#define ALIGENPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Generator using the TPythia interface (via AliPythia)
// to generate pp collisions.
// Using SetNuclei() also nuclear modifications to the structure functions
// can be taken into account. This makes, of course, only sense for the
// generation of the products of hard processes (heavy flavor, jets ...)
//
// andreas.morsch@cern.ch
//

#include "AliGenMC.h"
#include "AliPythia.h"

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
    virtual void    SetStrucFunc(StrucFunc_t func = kGRVHO) {fStrucFunc = func;}
    // select pt of hard scattering 
    virtual void    SetPtHard(Float_t ptmin = 0, Float_t ptmax = 1.e10)
	{fPtHardMin = ptmin; fPtHardMax = ptmax; }
    virtual void    SetYHard(Float_t ymin = -1.e10, Float_t ymax = 1.e10)
	{fYHardMin = ymin; fYHardMax = ymax; }
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
    virtual Float_t GetXsection() const {return fXsection;}
    // Getters
    virtual Process_t    GetProcess() {return fProcess;}
    virtual StrucFunc_t  GetStrucFunc() {return fStrucFunc;}
    virtual void         GetPtHard(Float_t& ptmin, Float_t& ptmax)
	{ptmin = fPtHardMin; ptmax = fPtHardMax = ptmax;}
    virtual Float_t      GetEnergyCMS() {return fEnergyCMS;}
    virtual void         GetNuclei(Int_t&  a1, Int_t& a2)
	{a1 = fNucA1; a2 = fNucA2;}
    virtual void         GetJetEtaRange(Float_t& etamin, Float_t& etamax)
	{etamin = fEtaMinJet; etamax = fEtaMaxJet;}
    virtual void         GetJetPhiRange(Float_t& phimin, Float_t& phimax)
	{phimin = fPhiMinJet*180./TMath::Pi(); phimax = fPhiMaxJet*180/TMath::Pi();}
    virtual void         GetGammaEtaRange(Float_t& etamin, Float_t& etamax)
	{etamin = fEtaMinGamma; etamax = fEtaMaxGamma;}
    virtual void         GetGammaPhiRange(Float_t& phimin, Float_t& phimax)
	{phimin = fPhiMinGamma*180./TMath::Pi(); phimax = fPhiMaxGamma*180./TMath::Pi();}
    //
    virtual void FinishRun();
    Bool_t CheckTrigger(TParticle* jet1, TParticle* jet2) const;
    
    // Assignment Operator
    AliGenPythia & operator=(const AliGenPythia & rhs);
 protected:
    // adjust the weight from kinematic cuts
    void   AdjustWeights();
    Int_t  GenerateMB();
    void   MakeHeader() const;    

    TClonesArray* fParticles;     //Particle  List
    
    Process_t   fProcess;         //Process type
    StrucFunc_t fStrucFunc;       //Structure Function
    Float_t     fEnergyCMS;       //Centre of mass energy
    Float_t     fKineBias;        //!Bias from kinematic selection
    Int_t       fTrials;          //!Number of trials for current event
    Int_t       fTrialsRun;       //!Number of trials for run
    Float_t     fQ;               //Mean Q
    Float_t     fX1;              //Mean x1
    Float_t     fX2;              //Mean x2
    Int_t       fNev;             //Number of events 
    Int_t       fFlavorSelect;    //Heavy Flavor Selection
    Float_t     fXsection;        //Cross-section
    AliPythia   *fPythia;         //!Pythia 
    Float_t     fPtHardMin;       //lower pT-hard cut 
    Float_t     fPtHardMax;       //higher pT-hard cut
    Float_t     fYHardMin;        //lower  y-hard cut 
    Float_t     fYHardMax;        //higher y-hard cut
    Int_t       fNucA1;           //mass number nucleus side 1
    Int_t       fNucA2;           //mass number nucleus side 2
    Bool_t      fFullEvent;       //!Write Full event if true
    AliDecayer  *fDecayer;        //!Pointer to the decayer instance
    Int_t       fDebugEventFirst; //!First event to debug
    Int_t       fDebugEventLast;  //!Last  event to debug
    Float_t     fEtaMinJet;       //Minimum eta of triggered Jet
    Float_t     fEtaMaxJet;       //Maximum eta of triggered Jet
    Float_t     fPhiMinJet;       //Minimum phi of triggered Jet
    Float_t     fPhiMaxJet;       //Maximum phi of triggered Jet

    Float_t     fEtaMinGamma;     // Minimum eta of triggered gamma
    Float_t     fEtaMaxGamma;     // Maximum eta of triggered gamma
    Float_t     fPhiMinGamma;     // Minimum phi of triggered gamma
    Float_t     fPhiMaxGamma;     // Maximum phi of triggered gamma

    StackFillOpt_t fStackFillOpt; // Stack filling with all particles with
                                  // that flavour or only with selected
                                  // parents and their decays
    Bool_t fFeedDownOpt;          // Option to set feed down from higher
                                  // quark families (e.g. b->c)
    Bool_t fFragmentation;        // Option to activate fragmentation by Pythia
    //

    CountMode_t fCountMode;        // Options for counting when the event will be finished.
    // fCountMode = kCountAll         --> All particles that end up in the
    //                                    stack are counted
    // fCountMode = kCountParents     --> Only selected parents are counted
    // fCountMode = kCountTrackabless --> Only particles flagged for tracking
    //                                     are counted
    //
    ClassDef(AliGenPythia,3) // AliGenerator interface to Pythia
};
#endif





