#ifndef ALIGENERATOR_H
#define ALIGENERATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////
//                                                       //
//  Class to generate the particles for the MC           //
//  The base class is empty                              //
//                                                       //
///////////////////////////////////////////////////////////

#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TGenerator.h"

typedef enum { none, perEvent, perTrack } VertexSmear_t;

class AliGenerator : public TNamed
{

protected:
    static TGenerator* fgMCEvGen;
    Float_t     fThetaMin;     //Minimum theta of generation in radians
    Float_t     fThetaMax;     //Maximum theta of generation in radians
    Float_t     fPhiMin;       //Minimum phi of generation in radians
    Float_t     fPhiMax;       //Maximum phi of generation in radians
    Float_t     fPMin;         //Minimum momentum of generation in GeV/c
    Float_t     fPMax;         //Minimum momentum of generation in GeV/c
    Float_t     fPtMin;        //Minimum transverse momentum
    Float_t     fPtMax;        //Maximum transverse momentum
    Float_t     fYMin;         //Minimum rapidity
    Float_t     fYMax;         //Maximum rapidity
    TArrayF     fVMin;         //Minimum Decaylength
    TArrayF     fVMax;         //Minimum Decaylength    
    Int_t       fNpart;        //Maximum number of particles per event
    Float_t     fParentWeight; //Parent Weight
    Float_t     fChildWeight;  //ChildWeight
    Int_t       fTrackit;      // Track the generated final state particle if 1
    Int_t       fAnalog;       //Flaf for anolog or pt-weighted generation
   //
    VertexSmear_t     fVertexSmear; //Vertex Smearing mode
    Int_t       fTrackIt;    // if 1 Track final state particles 
    TArrayF     fOrigin;     //Origin of event
    TArrayF     fOsigma;     //Sigma of the Origin of event

    enum {kThetaRange=1, kVertexRange=2, kPhiRange=4, kPtRange=8,
	  kYRange=32, kMomentumRange=16};

 public:
    AliGenerator();
    AliGenerator(Int_t npart);
    virtual ~AliGenerator();
    virtual void Init();
    virtual void SetOrigin(Float_t ox, Float_t oy, Float_t oz)
	{fOrigin[0]=ox;fOrigin[1]=oy;fOrigin[2]=oz;}
    virtual void SetOrigin(const TLorentzVector &o)
	{fOrigin[0]=o[0];fOrigin[1]=o[1];fOrigin[2]=o[2];}
    virtual void SetSigma(Float_t sx, Float_t sy, Float_t sz)
	{fOsigma[0]=sx;fOsigma[1]=sy;fOsigma[2]=sz;}
    virtual void SetMomentumRange(Float_t pmin=0, Float_t pmax=1.e10)
	{fPMin = pmin; fPMax = pmax; SetBit(kMomentumRange);}
    virtual void SetPtRange(Float_t ptmin=0, Float_t ptmax=1.e10)
	{fPtMin = ptmin; fPtMax = ptmax; SetBit(kPtRange);}
    virtual void SetPhiRange(Float_t phimin=-180., Float_t phimax=180)
	{fPhiMin = TMath::Pi()*phimin/180;
	fPhiMax = TMath::Pi()*phimax/180; SetBit(kPhiRange);}
    virtual void SetYRange(Float_t ymin=-100, Float_t ymax=100)
	{fYMin=ymin; fYMax=ymax; SetBit(kYRange);}
    virtual void SetVRange(Float_t vxmin, Float_t vxmax,
			   Float_t vymin, Float_t vymax,
			   Float_t vzmin, Float_t vzmax)
	{
	    fVMin[0]=vxmin; fVMin[1]=vymin; fVMin[2]=vzmin;
	    fVMax[0]=vxmax; fVMax[1]=vymax; fVMax[2]=vzmax;
	    SetBit(kVertexRange);
	}
    virtual void SetNumberParticles(Int_t npart=100) {fNpart=npart;}
    virtual Int_t NumberParticles() {return fNpart;}
    virtual void SetThetaRange(Float_t thetamin=0, Float_t thetamax=180)
	{fThetaMin = TMath::Pi()*thetamin/180;
	fThetaMax = TMath::Pi()*thetamax/180; SetBit(kThetaRange);}
    virtual void Generate()=0;
    virtual void SetParentWeight(Float_t wgt) {fParentWeight=wgt;}
    virtual void SetChildWeight(Float_t wgt)  {fChildWeight=wgt;}    
    virtual void SetAnalog(Int_t flag=1) {fAnalog=flag;}	
    virtual void SetVertexSmear(VertexSmear_t smear) {fVertexSmear = smear;}
    virtual void SetTrackingFlag(Int_t flag=1) {fTrackIt=flag;}
 	    
    virtual void SetMC(TGenerator *theMC) 
	{if (!fgMCEvGen) fgMCEvGen =theMC;}

  // Getters

    virtual void GetOrigin(Float_t &ox, Float_t &oy, Float_t &oz)
	{ox=fOrigin[0];oy=fOrigin[1];oz=fOrigin[2];}
    virtual void GetOrigin(TLorentzVector &o)
	{o[0]=fOrigin[0];o[1]=fOrigin[1];o[2]=fOrigin[2];o[3]=0;}

    ClassDef(AliGenerator,1)
};

#endif
