#ifndef ALIGENPARAM_H
#define ALIGENPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"
#include "TF1.h"
#include "TArrayI.h"
#include "AliPythia.h"
#include "GenTypeDefs.h"
//-------------------------------------------------------------
class AliGenParam : public AliGenerator
{
 public:
    AliGenParam();
    AliGenParam(Int_t npart, Param_t param);
    AliGenParam(Int_t npart, Param_t param,
		Double_t (*PtPara)(Double_t*, Double_t*),
		Double_t (*YPara )(Double_t*, Double_t*),
		Int_t    (*IpPara)()                      ); 
    virtual ~AliGenParam();
    virtual void Generate();
    virtual void Init();
    // select particle type
    virtual void SetParam(Param_t param=jpsi_p) {fParam=param;}
    // force decay type
    virtual void SetForceDecay(Decay_t decay=dimuon) {fForceDecay=decay;}
    virtual void SetWeighting(Weighting_t flag=analog) {fAnalog=flag;}	
    virtual void SetCutOnChild(Int_t flag=0) {fCutOnChild=flag;}
    virtual void SetChildMomentumRange(Float_t pmin=0, Float_t pmax=1.e10)
	{fChildPMin = pmin; fChildPMax = pmax;}
    virtual void SetChildPtRange(Float_t ptmin=0, Float_t ptmax=20.)
	{fChildPtMin = ptmin; fChildPtMax = ptmax;}
    virtual void SetChildPhiRange(Float_t phimin=-180., Float_t phimax=180)
	{fChildPhiMin = TMath::Pi()*phimin/180;
	fChildPhiMax = TMath::Pi()*phimax/180;}
    virtual void SetChildThetaRange(Float_t thetamin=0, Float_t thetamax=180)
	{fChildThetaMin = TMath::Pi()*thetamin/180;
	fChildThetaMax = TMath::Pi()*thetamax/180;}
    virtual void SetDeltaPt(Float_t delta=0.01) {fDeltaPt=delta;}
    
	    
 protected:
    Double_t (*fPtParaFunc)(Double_t*, Double_t*); //! Pointer to Pt parametrisation function
    Double_t (*fYParaFunc )(Double_t*, Double_t*); //! Pointer to Y parametrisation function
    Int_t    (*fIpParaFunc )();    //! Pointer to particle type parametrisation function
    TF1* fPtPara;              // Transverse momentum parameterisation
    TF1* fYPara;               // Rapidity parameterisation
    Param_t     fParam;        // Parameterisation type 
    Float_t     fdNdy0;        // central multiplicity per event
    Float_t     fYWgt;         // Y-weight
    Float_t     fPtWgt;        // Pt-weight
    Weighting_t fAnalog;       // Flag for anolog or pt-weighted generation
    Float_t     fBias;         // Biasing factor
    Int_t       fTrials;       // Number of trials
    Decay_t     fForceDecay;   // Decay channel forced
    Int_t       fCutOnChild;   // Cuts on decay products (children)  are enabled/disabled
    Float_t     fChildPtMin;   // Children minimum pT
    Float_t     fChildPtMax;   // Children maximum pT
    Float_t     fChildPMin;    // Children minimum p
    Float_t     fChildPMax;    // Children maximum p
    Float_t     fChildPhiMin;  // Children minimum phi
    Float_t     fChildPhiMax;  // Children maximum phi
    Float_t     fChildThetaMin;// Children minimum theta
    Float_t     fChildThetaMax;// Children maximum theta
    Float_t     fDeltaPt;      // pT sampling in steps of fDeltaPt
    TArrayI     fChildSelect;  // Children to be selected from decay products
    AliPythia   *fPythia;      // Pointer to pythia object for decays
 private:
    // check if particle is selected as child
    Bool_t ChildSelected(Int_t ip);
    // all kinematic selection goes here
    Bool_t KinematicSelection(TParticle *particle);

  ClassDef(AliGenParam,1) // Generator using parameterised pt- and y-distribution
};
#endif










