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
class AliGenHerwigEventHeader;
class TArrayI;
class TParticle;
class TClonesArray;


class AliGenHerwig : public AliGenMC

{
    enum {kNoTrigger, kHardProcesses, kDirectPhotons};
    enum {kHeJets = 1500, kHeDirectGamma = 1800};
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
    virtual void    SetPtHardMax(Double_t pt) {fPtHardMax=pt;}
    virtual void    SetWeightPower(Double_t ptweight) {fWeightPower=ptweight;}
    virtual void    SetPtRMS(Double_t pt) {fPtRMS=pt;}
    virtual void    SetMaxPr(Int_t i) {fMaxPr=i;}
    virtual void    SetMaxErrors(Int_t i) {fMaxErrors=i;}
    virtual void    SetSeed(UInt_t seed) {GetRandom()->SetSeed(seed);}
    virtual void    FinishRun();
    virtual void    FinishRunJimmy();
    virtual void    SetEnSoft(Double_t e) {fEnSoft=e;}

    virtual void    SetHardProcessFile(TString filename) {fFileName=TString(filename);};
    virtual Bool_t CheckParton(const TParticle* parton1, const TParticle* parton2);

    virtual void         GetPartonEtaRange(Float_t& etamin, Float_t& etamax) const
	{etamin = fEtaMinParton; etamax = fEtaMaxParton;}
    virtual void         GetPartonPhiRange(Float_t& phimin, Float_t& phimax) const
	{phimin = fPhiMinParton*180./TMath::Pi(); phimax = fPhiMaxParton*180/TMath::Pi();}
    virtual void         GetGammaEtaRange(Float_t& etamin, Float_t& etamax) const
	{etamin = fEtaMinGamma; etamax = fEtaMaxGamma;}
    virtual void         GetGammaPhiRange(Float_t& phimin, Float_t& phimax) const
	{phimin = fPhiMinGamma*180./TMath::Pi(); phimax = fPhiMaxGamma*180./TMath::Pi();}

    virtual void    SetPartonEtaRange(Float_t etamin = -20., Float_t etamax = 20.)
	{fEtaMinParton = etamin; fEtaMaxParton = etamax;}
    // Phi range for jet trigger
    virtual void    SetPartonPhiRange(Float_t phimin = 0., Float_t phimax = 360.)
	{fPhiMinParton = TMath::Pi()*phimin/180.; fPhiMaxParton = TMath::Pi()*phimax/180.;}
    // Eta range for gamma trigger 
    virtual void    SetGammaEtaRange(Float_t etamin = -20., Float_t etamax = 20.)
	{fEtaMinGamma = etamin; fEtaMaxGamma = etamax;}
    // Phi range for gamma trigger
    virtual void    SetGammaPhiRange(Float_t phimin = 0., Float_t phimax = 360.)
	{fPhiMinGamma = TMath::Pi()*phimin/180.; fPhiMaxGamma = TMath::Pi()*phimax/180.;}

 protected:
    Bool_t SelectFlavor(Int_t pid) const;
    void   MakeHeader();    
 protected:
    TString     fAutPDF;         // PDF group
    Int_t       fModPDF;         // PDF set
    StrucFunc_t fStrucFunc;      // Structure Function
    Int_t       fKeep;           // Flag to keep full event information
    Int_t       fDecaysOff;      // Flag to turn off decays of pi0, K_s, D, Lambda, sigma
    Int_t       fTrigger;        // Trigger type
    Int_t       fSelectAll;      // Flag to write the full event
    Int_t       fFlavor;         // Selected particle flavor 4: charm+beauty 5: beauty
    Float_t     fMomentum1;      // Momentum of projectile
    Float_t     fMomentum2;      // Momentum of target
    Float_t     fKineBias;       // Bias from kinematic selection
    Int_t       fTrials;         // Number of trials
    Float_t     fXsection;       // Cross-section
    THerwig6    *fHerwig;        // Herwig
    Int_t       fProcess;        // Process number
    Double_t    fPtHardMin;      // higher pT-hard cut
    Double_t    fPtHardMax;      // lower pT-hard cut
    Double_t    fPtHardGen;      // generated pT-hard 
    Double_t    fWeightPower;   //power of the event weight
    Double_t    fPtRMS;          // intrinsic pt of incoming hadrons
    Int_t       fMaxPr;          // maximum number of events to print out
    Int_t       fMaxErrors;      // maximum number of errors allowed
    Double_t    fEnSoft;          // change on soft energy distribution
    TString     fFileName;       //!Name of file to read from hard scattering
    Float_t     fEtaMinParton;         //Minimum eta of parton shower
    Float_t     fEtaMaxParton;         //Maximum eta of parton shower
    Float_t     fPhiMinParton;         //Minimum phi of parton shower
    Float_t     fPhiMaxParton;         //Maximum phi of parton shower
    Float_t     fEtaMinGamma;       // Minimum eta of triggered gamma
    Float_t     fEtaMaxGamma;       // Maximum eta of triggered gamma
    Float_t     fPhiMinGamma;       // Minimum phi of triggered gamma
    Float_t     fPhiMaxGamma;       // Maximum phi of triggered gamma
    AliGenHerwigEventHeader* fHeader;  //! Event header
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
    Bool_t DaughtersSelection(const TParticle* iparticle, const TClonesArray* particles);
    // check if stable
    Bool_t Stable(const TParticle*  particle) const;

    void InitPDF();

    ClassDef(AliGenHerwig,2) // AliGenerator interface to Herwig
};
#endif





