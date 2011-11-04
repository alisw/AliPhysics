#ifndef ALIGENAMPT_H
#define ALIGENAMPT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenMC.h"
#include <TString.h>

class TAmpt;
class TParticle;
class TClonesArray;
class TGraph;
class AliGenHijingEventHeader;

class AliGenAmpt : public AliGenMC
{
  enum {kNoTrigger, kHardProcesses, kDirectPhotons};

  public:
    AliGenAmpt();
    AliGenAmpt(Int_t npart);
    virtual ~AliGenAmpt();

    virtual void    SetDecay(Bool_t b)                {fDecay = b;}
    virtual TAmpt  *Ampt() { return fAmpt; }
    virtual void    Generate();
    virtual void    Init();
    virtual void    SetEnergyCMS(Float_t energy=5500) {fEnergyCMS=energy;}
    virtual void    SetReferenceFrame(TString frame="CMS")
	{fFrame=frame;}
    virtual void    SetImpactParameterRange(Float_t bmin = 0, Float_t bmax = 15.)
	{fMinImpactParam=bmin; fMaxImpactParam=bmax;}
    virtual void    KeepFullEvent()                   {fKeep      =1;     }
    virtual void    SetJetQuenching(Int_t flag=1)     {fQuench     = flag;}
    virtual void    SetShadowing(Int_t flag=1)        {fShadowing  = flag;}
    virtual void    SetDecaysOff(Int_t flag=1)        {fDecaysOff  = flag;}
    virtual void    SetTrigger(Int_t flag=kNoTrigger) {fTrigger    = flag;}
    virtual void    SetFlavor(Int_t flag=0)           {fFlavor     = flag;}
    virtual void    SetEvaluate(Int_t flag=0)         {fEvaluate   = flag;}
    virtual void    SetSelectAll(Int_t flag=0)        {fSelectAll  = flag;}
    virtual void    SetRadiation(Int_t flag=3)        {fRadiation  = flag;}    
    virtual void    SetSpectators(Int_t spects=1)     {fSpectators = spects;}
    virtual void    SetDecayer(AliDecayer *decayer)   {fDecayer = decayer;}
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
    virtual void    SetIsoft(Int_t i)                         {fIsoft = i;  }
    virtual void    SetNtMax(Int_t max)                       {fNtMax = max;}
    virtual void    SetIpop(Int_t pop)                        {fIpop  = pop;}
    virtual void    SetXmu(Float_t m)                         {fXmu   = m;  }
    virtual void    SetAlpha(Float_t alpha)                   {fAlpha = alpha;            }
    virtual void    SetStringFrag(Float_t a, Float_t b)       {fStringA = a; fStringB = b;}
	    
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
    const TClonesArray *GetParticles() const { return &fParticles; }

    // Physics Routines
    virtual Bool_t  ProvidesCollisionGeometry() const {return kTRUE;}
    virtual void    EvaluateCrossSections();
    virtual TGraph* CrossSection()     {return fDsigmaDb;}
    virtual TGraph* BinaryCollisions() {return fDnDb;}
    virtual Bool_t  CheckTrigger();
    virtual Bool_t  IsThisAKnownParticle(TParticle *thisGuy);

  protected:
    Bool_t      SelectFlavor(Int_t pid);
    void        MakeHeader();

    AliDecayer              *fDecayer;
    TString                  fFrame;           // Reference frame 
    Float_t                  fMinImpactParam;  // minimum impact parameter
    Float_t                  fMaxImpactParam;  // maximum impact parameter	
    Int_t                    fKeep;            // Flag to keep full event information
    Int_t                    fQuench;          // Flag to switch on jet quenching
    Int_t                    fShadowing;       // Flag to switch on nuclear effects on parton distribution function
    Int_t                    fDecaysOff;       // Flag to turn off decays of pi0, K_s, D, Lambda, sigma
    Int_t                    fTrigger;         // Trigger type
    Int_t                    fEvaluate;        // Evaluate total and partial cross-sections
    Int_t                    fSelectAll;       // Flag to write the full event
    Int_t                    fFlavor;          // Selected particle flavor 4: charm+beauty 5: beauty
    Float_t                  fKineBias;        // Bias from kinematic selection
    Int_t                    fTrials;          // Number of trials
    Float_t                  fXsection;        // Cross-section
    TAmpt                   *fAmpt;            //!Ampt
    Float_t                  fPtHardMin;       // lower pT-hard cut 
    Float_t                  fPtHardMax;       // higher pT-hard cut
    Int_t                    fSpectators;      // put spectators on stack
    TGraph*                  fDsigmaDb;        // dSigma/db for the system
    TGraph*                  fDnDb;            // dNBinaryCollisions/db
    Float_t                  fPtMinJet;        // Minimum Pt of triggered Jet
    Float_t                  fEtaMinJet;       // Minimum eta of triggered Jet
    Float_t                  fEtaMaxJet;       // Maximum eta of triggered Jet
    Float_t                  fPhiMinJet;       // At least one of triggered Jets must be in this
    Float_t                  fPhiMaxJet;       // phi range
    Int_t                    fRadiation;       // Flag to switch on/off initial and final state radiation
    Int_t                    fSimpleJet;       // Flag to produce simple tiggered jet topology
    Int_t                    fNoGammas;        // Don't write gammas if flag "on"
    Int_t 	             fProjectileSpecn; // Num. of spectator neutrons from projectile nucleus
    Int_t 	             fProjectileSpecp; // Num. of spectator protons from projectile nucleus
    Int_t 	             fTargetSpecn;     // Num. of spectator neutrons from target nucleus
    Int_t 	             fTargetSpecp;     // Num. of spectator protons from target nucleus
    Int_t                    fLHC;             // Assume LHC as lab frame
    Bool_t                   fRandomPz;        // Randomise sign of pz event by event
    Bool_t                   fNoHeavyQuarks;   // If true no heavy quarks are produced
    Int_t                    fIsoft;           // ISOFT (D=1): select Default AMPT or String Melting
    Int_t                    fNtMax;           // NTMAX: number of timesteps (D=150)
    Int_t                    fIpop;            // (D=1,yes;0,no) flag for popcorn mechanism(netbaryon stopping)
    Float_t                  fXmu;             // parton screening mass in fm^(-1) (D=3.2264d0)
    Float_t                  fAlpha;           // alpha running (fixed) coupling
    Float_t                  fStringA;         // string frag parameter A
    Float_t                  fStringB;         // string frag parameter B
    Float_t                  fEventTime;       // The event time
    AliGenHijingEventHeader *fHeader;          // header
    Bool_t                   fDecay;           // decay "long-lived" particles

  private:
    AliGenAmpt(const AliGenAmpt &Ampt);
    AliGenAmpt &  operator=(const AliGenAmpt & rhs);

    // adjust the weight from kinematic cuts
    void   AdjustWeights();
    // check seleted daughters
    Bool_t DaughtersSelection(TParticle* iparticle);
    // check if stable
    Bool_t Stable(TParticle*  particle) const;

    ClassDef(AliGenAmpt, 4) // AliGenerator interface to Ampt
};
#endif
