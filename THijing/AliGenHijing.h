#ifndef ALIGENHIJING_H
#define ALIGENHIJING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Generator using HIJING as an external generator
// The main HIJING options are accessable for the user through this interface.
// andreas.morsch@cern.ch

#include "AliGenMC.h"
#include "AliGenHijingEventHeader.h"
#include <TString.h>

class THijing;
class TParticle;
class TClonesArray;
class TGraph;
class TF1;


class AliGenHijing : public AliGenMC
{
    enum {kNoTrigger, kHardProcesses, kDirectPhotons};

 public:
    AliGenHijing();
    AliGenHijing(Int_t npart);
    virtual ~AliGenHijing();
    virtual void    Generate();
    virtual void    Init();
    virtual void    SetSeed(UInt_t seed);
    // set centre of mass energy
    virtual void    SetEnergyCMS(Float_t energy=5500) {fEnergyCMS=energy;}
    virtual void    SetReferenceFrame(TString frame="CMS")
	{fFrame=frame;}

    virtual void    SetImpactParameterRange(Float_t bmin = 0, Float_t bmax = 15.)
	{fMinImpactParam=bmin; fMaxImpactParam=bmax;}
    virtual void    KeepFullEvent();
    virtual void    SetJetQuenching(Int_t flag=1)     {fQuench     = flag;}
    virtual void    SetShadowing(Int_t flag=1)        {fShadowing  = flag;}
    virtual void    SetDecaysOff(Int_t flag=1)        {fDecaysOff  = flag;}
    virtual void    SetTrigger(Int_t flag=kNoTrigger) {fTrigger    = flag;}
    virtual void    SetFlavor(Int_t flag=0)           {fFlavor     = flag;}
    virtual void    SetEvaluate(Int_t flag=0)         {fEvaluate   = flag;}
    virtual void    SetSelectAll(Int_t flag=0)        {fSelectAll  = flag;}
    virtual void    SetRadiation(Int_t flag=3)        {fRadiation  = flag;}    
    virtual void    SetSpectators(Int_t spects=1)     {fSpectators = spects;}
    virtual void    SetPtHardMin(Float_t ptmin)       {fPtHardMin  = ptmin;}
    virtual void    SetPtHardMax(Float_t ptmax)       {fPtHardMax  = ptmax;}
    virtual void    SetPtJet(Float_t ptmin)           {fPtMinJet   = ptmin;}
    virtual void    SetSimpleJets(Int_t flag=0)       {fSimpleJet  = flag;}
    virtual void    SetNoGammas(Int_t flag=0)         {fNoGammas   = flag;}
	    
    virtual void    SetJetEtaRange(Float_t etamin = -20., Float_t etamax = 20.)
	{fEtaMinJet = etamin; fEtaMaxJet = etamax;}
    virtual void    SetJetPhiRange(Float_t phimin = -180., Float_t phimax = 180.)
	{fPhiMinJet = TMath::Pi()*phimin/180.; fPhiMaxJet = TMath::Pi()*phimax/180.;}
    virtual void    SetBoostLHC(Int_t flag = 0)       {fLHC        = flag;}
    virtual void    SetRandomPz(Bool_t flag = 0)      {fRandomPz   = flag;}
    virtual void    SwitchOffHeavyQuarks(Bool_t flag = kTRUE) {fNoHeavyQuarks = flag;}
    virtual void    SetSigmaNN(Float_t val)           {fSigmaNN    = val;}    
    virtual void    SetNoElas(Bool_t b)               {fNoElas     = b; }	  
    virtual void    SetDataDrivenSpectators()	      {fDataFragmentation = kTRUE;}

// Getters
    virtual TString GetReferenceFrame()  const {return fFrame;}
    virtual void    GetImpactParameterRange(Float_t& bmin, Float_t& bmax) const
	{bmin = fMinImpactParam; bmax = fMaxImpactParam;}
    virtual Int_t   GetJetQuenching()                    const {return fQuench;}
    virtual Int_t   GetShadowing()                       const {return fShadowing;}
    virtual Int_t   GetTrigger()                         const {return fTrigger;}
    virtual Int_t   GetFlavor()                          const {return fFlavor;}
    virtual Int_t   GetRadiation()                       const {return fRadiation;}    
    virtual Int_t   GetSpectators()                      const {return fSpectators;}
    virtual Float_t GetPtHardMin()                       const {return fPtHardMin;}
    virtual Float_t GetPtHardMax()                       const {return fPtHardMax;}
    virtual Float_t GetPtJet()                           const {return fPtMinJet;}
    virtual void    GetJetEtaRange(Float_t& etamin, Float_t& etamax)      const 
	{etamin = fEtaMinJet; etamax = fEtaMaxJet;}
    virtual void    GetJetPhiRange(Float_t& phimin, Float_t& phimax)      const
	{phimin = fPhiMinJet*180./TMath::Pi(); phimax = fPhiMaxJet*180./TMath::Pi();}
     THijing       *GetTHijing()                         const {return fHijing;}
    virtual Float_t GetSigmaNN()                         const {return fSigmaNN;}
    virtual Bool_t  HasDataDrivenSpectators()	         const {return fDataFragmentation;}
    virtual Int_t   GetFreeProjSpecn()	         const {return  fFreeProjSpecn;}
    virtual Int_t   GetFreeProjSpecp()	         const {return  fFreeProjSpecp;}
    virtual Int_t   GetFreeTargSpecn()	         const {return  fFreeTargSpecn;}
    virtual Int_t   GetFreeTargSpecp()	         const {return  fFreeTargSpecp;}

// Physics Routines
    virtual Bool_t  ProvidesCollisionGeometry() const {return kTRUE;}
    virtual void    EvaluateCrossSections();
    virtual TGraph* CrossSection()     {return fDsigmaDb;}
    virtual TGraph* BinaryCollisions() {return fDnDb;}
    virtual Bool_t  CheckTrigger();
//
 protected:
    Bool_t SelectFlavor(Int_t pid);
    void   MakeHeader();
 protected:
    TString     fFrame;          // Reference frame 
    Float_t     fMinImpactParam; // minimum impact parameter
    Float_t     fMaxImpactParam; // maximum impact parameter	
    Int_t       fKeep;           // Flag to keep full event information
    Int_t       fQuench;         // Flag to switch on jet quenching
    Int_t       fShadowing;      // Flag to switch on nuclear effects on parton distribution function
    Int_t       fDecaysOff;      // Flag to turn off decays of pi0, K_s, D, Lambda, sigma
    Int_t       fTrigger;        // Trigger type
    Int_t       fEvaluate;       // Evaluate total and partial cross-sections
    Int_t       fSelectAll;      // Flag to write the full event
    Int_t       fFlavor;         // Selected particle flavor 4: charm+beauty 5: beauty
    Float_t     fKineBias;       // Bias from kinematic selection
    Int_t       fTrials;         // Number of trials
    Float_t     fXsection;       // Cross-section
    THijing    *fHijing;         //!Hijing
    Float_t     fPtHardMin;      // lower pT-hard cut 
    Float_t     fPtHardMax;      // higher pT-hard cut
    Int_t       fSpectators;     // put spectators on stack
    TGraph*     fDsigmaDb;       // dSigma/db for the system
    TGraph*     fDnDb;           // dNBinaryCollisions/db
    Float_t     fPtMinJet;       // Minimum Pt of triggered Jet
    Float_t     fEtaMinJet;      // Minimum eta of triggered Jet
    Float_t     fEtaMaxJet;      // Maximum eta of triggered Jet
    Float_t     fPhiMinJet;      // At least one of triggered Jets must be in this
    Float_t     fPhiMaxJet;      // phi range
    Int_t       fRadiation;      // Flag to switch on/off initial and final state radiation
    Int_t       fSimpleJet;      // Flag to produce simple tiggered jet topology
    Int_t       fNoGammas;       // Don't write gammas if flag "on"
// ZDC proposal (by Chiara) to store num. of SPECTATORS protons and neutrons
    Int_t 	fProjectileSpecn;// Num. of spectator neutrons from projectile nucleus
    Int_t 	fProjectileSpecp;// Num. of spectator protons from projectile nucleus
    Int_t 	fTargetSpecn;	 // Num. of spectator neutrons from target nucleus
    Int_t 	fTargetSpecp;	 // Num. of spectator protons from target nucleus
    Int_t       fLHC;            // Assume LHC as lab frame
    Bool_t      fRandomPz;       // Randomise sign of pz  event by event
    Bool_t      fNoHeavyQuarks;  // If true no heavy quarks are produced
    AliGenHijingEventHeader     fHeader; // MC Header
    Float_t     fSigmaNN;        // If not -1 set sigmaNN (HIPR1) 
    Bool_t      fNoElas;         // If true switch off elastic scattering
    // Added by Chiara to consider fragments production from data
    Bool_t      fDataFragmentation; // Spectators w. data driven correction
    Int_t 	fFreeProjSpecn;// Num. of spectator neutrons from projectile nucleus
    Int_t 	fFreeProjSpecp;// Num. of spectator protons from projectile nucleus
    Int_t 	fFreeTargSpecn;    // Num. of spectator neutrons from target nucleus
    Int_t 	fFreeTargSpecp;    // Num. of spectator protons from target nucleus
    TF1		*fFragmNeutrons;   // data driven correction for nuclear fragment formation
    TF1		*fFragmProtons;   // data driven correction for nuclear fragment formation

 private:
    AliGenHijing(const AliGenHijing &Hijing);
    AliGenHijing &  operator=(const AliGenHijing & rhs);

    // adjust the weight from kinematic cuts
    void   AdjustWeights();
    // check seleted daughters
    Bool_t DaughtersSelection(const TParticle* iparticle);
    // check if stable
    Bool_t Stable(const TParticle*  particle) const;
    // calculate no. of free spectators from data driven parametrization
    Int_t FreeSpectatorsn(Float_t b, Int_t nSpecn);
    Int_t FreeSpectatorsp(Float_t b, Int_t nSpecp);
    
    ClassDef(AliGenHijing, 10) // AliGenerator interface to Hijing
};
#endif
