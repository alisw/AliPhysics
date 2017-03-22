#ifndef ALIPYTHIA_H
#define ALIPYTHIA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPythia6.h>
#include <AliRndm.h>
#include <AliStructFuncType.h>
#include "PythiaProcesses.h"
#include "AliOmegaDalitz.h"
#include "AliDecayerExodus.h"
class AliFastGlauber;
class AliQuenchingWeights;

class AliPythia : public TPythia6, public AliRndm
{

 public:
    virtual ~AliPythia(){;}
    // convert to compressed code and print result (for debugging only)
    virtual Int_t CheckedLuComp(Int_t kf);
    // Pythia initialisation for selected processes
    virtual void ProcInit
	(Process_t process, Float_t energy, StrucFunc_t strucfunc, Int_t itune = -1);
    // treat protons as inside nuclei
    virtual void  SetNuclei(Int_t a1, Int_t a2, Int_t pdf);
    // Set colliding nuclei ("p","n",...)
    virtual void  SetCollisionSystem(TString projectile, TString target) { fProjectile = projectile; fTarget = target; }
    // Print particle properties
    virtual void PrintParticles();
    virtual void ResetDecayTable();
    virtual void SetDecayTable();
    virtual void SetWeightPower(Double_t pow); // use weighted cross sections
    virtual void Pyevnw();
    virtual void Pycell(Int_t& nclus);
    virtual void Pyclus(Int_t& nclus);
    virtual void Pyshow(Int_t ip1, Int_t ip2, Double_t qmax);
    virtual void Pyshowq(Int_t ip1, Int_t ip2, Double_t qmax);
    virtual void Pyrobo(Int_t imi, Int_t ima, Double_t the, Double_t phi, Double_t bex, Double_t bey, Double_t bez);
    virtual void Pytune(Int_t itune);
    virtual void Py2ent(Int_t idx, Int_t pdg1, Int_t pdg2, Double_t p);
    virtual void InitQuenching(Float_t bmin, Float_t bmax, Float_t k, Int_t iECMethod, Float_t zmax = 0.97, Int_t ngmax = 30);
    virtual void SetPyquenParameters(Double_t t0, Double_t tau0, Int_t nf, Int_t iengl, Int_t iangl);
    virtual void Pyquen(Double_t a, Int_t ibf, Double_t b);
    virtual void Qpygin0();
    virtual void GetQuenchingParameters(Double_t& xp, Double_t& yp, Double_t z[4]);
    virtual Int_t GetNMPI();
    // return instance of the singleton
    static  AliPythia* Instance();
    virtual void Quench();
    void DalitzDecays();
    // Dalitz and resonance decays from EXODUS
    void PizeroDalitz();
    void EtaDalitz();
    void RhoDirect();
    void OmegaDalitz();
    void OmegaDirect();
    void EtaprimeDalitz();
    void PhiDalitz();
    void PhiDirect();
    void JPsiDirect();
    // Assignment Operator
    AliPythia & operator=(const AliPythia & rhs);
    void Copy(TObject&) const;
 protected:
    Process_t             fProcess;           // Process type
    Float_t               fEcms;              // Centre of mass energy
    StrucFunc_t           fStrucFunc;         // Structure function
    TString               fProjectile;        // Projectile
    TString               fTarget;            // Target
    Int_t                 fDefMDCY[501];      //  ! Default decay switches per particle
    Int_t                 fDefMDME[2001];     //  ! Default decay switches per mode
    Double_t              fZQuench[4];        //  ! Quenching fractions for this even
    Double_t              fXJet;              //  ! Jet production point X
    Double_t              fYJet;              //  ! Jet production point Y
    Int_t                 fNGmax;             //    Maximum number of radiated gluons in quenching
    Float_t               fZmax;              //    Maximum energy loss in quenching
    AliFastGlauber*       fGlauber;           //  ! The Glauber model
    AliQuenchingWeights*  fQuenchingWeights;  //  ! The Quenching Weights model
    Int_t                 fItune;             //  ! Pythia tune 
    AliOmegaDalitz        fOmegaDalitz;       //  ! omega dalitz decayer
    AliDecayerExodus      fExodus;            // ! EXODUS decayer
    static AliPythia*     fgAliPythia;        // Pointer to single instance
 private: 
    AliPythia();
    AliPythia(const AliPythia& pythia);
    void ConfigHeavyFlavor();
    void AtlasTuning();
    void AtlasTuningMC09();
    ClassDef(AliPythia,1) //ALICE UI to PYTHIA
};

#endif





