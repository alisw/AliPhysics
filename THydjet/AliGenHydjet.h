#ifndef ALIGENHYDJET_H
#define ALIGENHYDJET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Generator using Hydjet as an external generator
// The main Hydjet options are accessable for the user through this interface.
// rafael.diaz.valdes@cern.ch


#include "AliGenMC.h"
#include <TString.h>

class THydjet;
class TParticle;
class TClonesArray;

class AliGenHydjet : public AliGenMC
{

 public:
    AliGenHydjet(Int_t npart=0);
    virtual ~AliGenHydjet();
    virtual void    Generate();
    virtual void    Init();
    // set initial beam parameters
    virtual void    SetEnergyCMS(Float_t energy=5500.) {fEnergyCMS=energy;}
    virtual void    SetReferenceFrame(TString frame="CMS") {fFrame=frame;}
    virtual void    SetProjectileWeigth(Int_t a=207){fAtomicWeigth=a;}
    virtual void    SetCentralityType(Int_t ifb=0){fIfbtype = ifb;}
    virtual void    SetFixedImpactParameter(Float_t bfix=0) {fFixImpactParam = bfix;}
    virtual void    SetImpactParameterRange(Float_t bmin = 0., Float_t bmax = 1.)
                                     {fMinImpactParam=bmin; fMaxImpactParam=bmax;}
    virtual void    SetMeanSoftMultiplicity(Int_t mult=20000){fSoftMult=mult;}
    // set hydro parameters
    virtual void    SetJetProduction(Int_t nhsel = 2){fJetProd=nhsel;}
    virtual void    SetMaxLongitudinalFlow(Float_t yflow = 5.){fYflow=yflow;}
    virtual void    SetMaxTransverseFlow(Float_t tflow = 1.){fTflow=tflow;}
    virtual void    SetSoftMultFraction(Float_t softfract = 1.){fSoftFract=softfract;}
    // set input PYTHIA parameters
    virtual void    SetDiJetProd(Int_t flag=1){fDijetProd=flag;}
    virtual void    SetMinPtHard(Float_t ptmin=10.){fMinPtHard=ptmin;}
    virtual void    SetStructFunction(Int_t mstp =7){fStructFunction=mstp;} // CTEQ5M
    virtual void    SetMultipleInteractions(Int_t flag=0){fMultipleInt=flag;}

    virtual void    SetFlavor(Int_t flag=0)           {fFlavor     = flag;}
    virtual void    SetSelectAll(Int_t flag=0)        {fSelectAll  = flag;}
    virtual void    AddHeader(AliGenEventHeader* header);

// Getters
   // virtual TString GetReferenceFrame()  const {return fFrame;}

//
    AliGenHydjet &  operator=(const AliGenHydjet & rhs);
 protected:
    Bool_t SelectFlavor(Int_t pid);
    void   MakeHeader();
 protected:
    //initial parameters
    TString     fFrame;          // Reference frame
    Float_t     fAtomicWeigth;   // Projectile-Target atomic weight
    Int_t       fIfbtype;        // centrality type
    Float_t     fFixImpactParam; // fixed impact parameter
    Float_t     fMinImpactParam; // minimum impact parameter
    Float_t     fMaxImpactParam; // maximum impact parameter
    Int_t       fSoftMult;       // mean soft multiplicity
    //hydro parameters
    Int_t       fJetProd;        // flag Jet production (nhsel)
    Float_t     fYflow;          // max longitudinal flow
    Float_t     fTflow;          // max transverse flow
    Float_t     fSoftFract;      // Soft multiplicity fraction
    // PYTHIA parameters
    Int_t       fDijetProd;      // flag dijet production
    Float_t     fMinPtHard;      // min pt hard
    Int_t       fStructFunction; // Structure Function (default CTEQ5M)
    Int_t       fMultipleInt;    // flag multiple interaction

    THydjet    *fHydjet;         //!Hydjet

    Int_t       fSelectAll;      // Flag to write the full event
    Int_t       fFlavor;         // Selected particle flavor 4: charm+beauty 5: beauty

 private:
    AliGenHydjet(const AliGenHydjet &Hijing);
    void Copy(TObject &rhs) const;
    // check if stable
    Bool_t Stable(TParticle*  particle) const;

    ClassDef(AliGenHydjet, 1) // AliGenerator interface to Hydjet
};
#endif





