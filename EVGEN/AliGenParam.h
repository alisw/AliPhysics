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
    AliGenParam(Int_t npart, const AliGenLib * Library, Int_t param,   const char*  tname = 0);
    AliGenParam(Int_t npart, Int_t param, const char* tname = 0, const char*  name  = 0);
    AliGenParam(Int_t npart, Int_t param,
		Double_t (*PtPara)(const Double_t*, const Double_t*),
		Double_t (*YPara )(const Double_t*, const Double_t*),
	      Double_t (*V2Para)(const Double_t*, const Double_t*),
		Int_t    (*IpPara)(TRandom*)           );
     
    virtual ~AliGenParam();
    virtual void GenerateN(Int_t ntimes);
    virtual void Generate();
    virtual void Init();
    // select particle type
    virtual void SetParam(Int_t param) {fParam = param;}
    //Setting the flag for Background transportation while using SetForceDecay()
    void SetSelectAll(Bool_t selectall) {fSelectAll = selectall;}
    // force decay type
    virtual void SetWeighting(Weighting_t flag = kAnalog) {fAnalog = flag;}	
    virtual void SetDeltaPt(Float_t delta=0.01) {fDeltaPt = delta;}
    virtual void SetDecayer(AliDecayer* decayer) {fDecayer = decayer;}
    virtual void SetForceGammaConversion(Bool_t force=kTRUE) {fForceConv = force;}
    virtual void Draw(const char * opt);
    TF1 *  GetPt() { return fPtPara;}
    TF1 *  GetY() {return fYPara;}
    Float_t GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax, Float_t phiMin, Float_t phiMax);

    static double ScreenVar(double Z, double e0, double eps){ return 136/pow(Z,0.333333)*e0/eps/(1-eps); }
    static double ScreenFunc1(double d);
    static double ScreenFunc2(double d);
    static double AuxScreenFunc1(double d, double Fz){ return 3*ScreenFunc1(d)-ScreenFunc2(d)-Fz; }
    static double AuxScreenFunc2(double d, double Fz){ return 3*0.5*ScreenFunc1(d)-0.5*ScreenFunc2(d)-Fz; }
    static TVector3 OrthogonalVector(TVector3 &inVec);
    double EnergyFraction(double Z, double E);
    double PolarAngle(double E);
    Int_t ForceGammaConversion(TClonesArray *particles, Int_t nPart);
  
protected:
    Double_t (*fPtParaFunc)(const Double_t*, const Double_t*); //! Pointer to Pt parametrisation function
    Double_t (*fYParaFunc )(const Double_t*, const Double_t*); //! Pointer to Y parametrisation function
    Int_t    (*fIpParaFunc )(TRandom*);    //! Pointer to particle type parametrisation function
  Double_t (*fV2ParaFunc )(const Double_t*, const Double_t*);//! Pointer to V2 parametrisation function
    TF1* fPtPara;              // Transverse momentum parameterisation
    TF1* fYPara;               // Rapidity parameterisation
  TF1*        fV2Para;       // v2 parametrization
  TF1*        fdNdPhi;       // Phi distribution depending on v2
    Int_t       fParam;        // Parameterisation type 
    Float_t     fdNdy0;        // central multiplicity per event
    Float_t     fYWgt;         // Y-weight
    Float_t     fPtWgt;        // Pt-weight
    Float_t     fBias;         // Biasing factor
    Int_t       fTrials;       // Number of trials
    Float_t     fDeltaPt;      // pT sampling in steps of fDeltaPt
    Bool_t      fSelectAll;    // Flag for transportation of Background while using SetForceDecay()
    AliDecayer  *fDecayer;     // ! Pointer to pythia object for decays
    Bool_t      fForceConv;    //

private:
    AliGenParam(const AliGenParam &Param);
    AliGenParam & operator=(const AliGenParam & rhs);

    ClassDef(AliGenParam, 2) // Generator using parameterised pt- and y-distribution
};
#endif










