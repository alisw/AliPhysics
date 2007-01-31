#ifndef ALIGENPARAM_H
#define ALIGENPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class to generate particles from using parametrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
// (For example AliGenMUONlib)
//
// andreas.morsch@cern.ch
//

#include "AliGenMC.h"

class AliPythia;
class TParticle;
class AliGenLib;
class TF1;

typedef enum { kAnalog, kNonAnalog} Weighting_t;
//-------------------------------------------------------------
class AliGenParam : public AliGenMC
{
 public:
    AliGenParam();
    AliGenParam(Int_t npart, AliGenLib * Library, Int_t param, char* tname = 0);
    AliGenParam(Int_t npart, Int_t param, const char* tname = 0, const char* name = 0);
    AliGenParam(Int_t npart, Int_t param,
		Double_t (*PtPara)(Double_t*, Double_t*),
		Double_t (*YPara )(Double_t*, Double_t*),
		Int_t    (*IpPara)(TRandom*)           );
     
    virtual ~AliGenParam();
    virtual void Generate();
    virtual void Init();
    // select particle type
    virtual void SetParam(Int_t param) {fParam = param;}
    // force decay type
    virtual void SetWeighting(Weighting_t flag = kAnalog) {fAnalog = flag;}	
    virtual void SetDeltaPt(Float_t delta=0.01) {fDeltaPt = delta;}
    virtual void SetDecayer(AliDecayer* decayer) {fDecayer = decayer;}
    virtual void Draw(const char * opt);
    TF1 *  GetPt() { return fPtPara;}
    TF1 *  GetY() {return fYPara;}
    Float_t GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax, Float_t phiMin, Float_t phiMax);

 protected:
    Double_t (*fPtParaFunc)(Double_t*, Double_t*); //! Pointer to Pt parametrisation function
    Double_t (*fYParaFunc )(Double_t*, Double_t*); //! Pointer to Y parametrisation function
    Int_t    (*fIpParaFunc )(TRandom*);    //! Pointer to particle type parametrisation function
    TF1* fPtPara;              // Transverse momentum parameterisation
    TF1* fYPara;               // Rapidity parameterisation
    Int_t       fParam;        // Parameterisation type 
    Float_t     fdNdy0;        // central multiplicity per event
    Float_t     fYWgt;         // Y-weight
    Float_t     fPtWgt;        // Pt-weight
    Float_t     fBias;         // Biasing factor
    Int_t       fTrials;       // Number of trials
    Float_t     fDeltaPt;      // pT sampling in steps of fDeltaPt
    AliDecayer  *fDecayer;     // ! Pointer to pythia object for decays

 private:
    AliGenParam(const AliGenParam &Param);
    AliGenParam & operator=(const AliGenParam & rhs);

    ClassDef(AliGenParam,1) // Generator using parameterised pt- and y-distribution
};
#endif










