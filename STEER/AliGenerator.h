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

class TGenerator;

#include "TLorentzVector.h"
#include "TArrayF.h"
#include "AliRndm.h"

typedef enum { kNoSmear, kPerEvent, kPerTrack } VertexSmear_t;
typedef enum { kExternal, kInternal}            VertexSource_t;

class AliGenerator : public TNamed, public AliRndm
{

 public:
    AliGenerator();
    AliGenerator(Int_t npart);
    AliGenerator(const AliGenerator &gen);
    virtual ~AliGenerator();
    virtual void Init();
    void Copy(AliGenerator &gen) const;
    virtual void SetOrigin(Float_t ox, Float_t oy, Float_t oz);
    virtual void SetOrigin(const TLorentzVector &o);
    virtual void SetSigma(Float_t sx, Float_t sy, Float_t sz);
    virtual void SetMomentumRange(Float_t pmin=0, Float_t pmax=1.e10);
    virtual void SetPtRange(Float_t ptmin=0, Float_t ptmax=100);
    virtual void SetPhiRange(Float_t phimin=-180., Float_t phimax=180);
    virtual void SetYRange(Float_t ymin=-100, Float_t ymax=100);
    virtual void SetVRange(Float_t vxmin, Float_t vxmax,
			   Float_t vymin, Float_t vymax,
			   Float_t vzmin, Float_t vzmax);
    virtual void SetNumberParticles(Int_t npart=100) {fNpart=npart;}
    virtual Int_t NumberParticles() const {return fNpart;}
    virtual void SetThetaRange(Float_t thetamin=0, Float_t thetamax=180);
    virtual void Generate()=0;
    virtual void SetParentWeight(Float_t wgt) {fParentWeight=wgt;}
    virtual void SetChildWeight(Float_t wgt)  {fChildWeight=wgt;}    
    virtual void SetAnalog(Int_t flag=1) {fAnalog=flag;}	
    virtual void SetVertexSmear(VertexSmear_t smear) {fVertexSmear = smear;}
    virtual void SetVertexSource(VertexSource_t smear) {fVertexSource = kInternal;}    
    virtual void SetTrackingFlag(Int_t flag=1) {fTrackIt=flag;}
    void Vertex();
    void VertexExternal();
    virtual void VertexInternal();
    virtual void FinishRun(){;}
    
    virtual void SetMC(TGenerator *theMC) 
	{if (!fgMCEvGen) fgMCEvGen =theMC;}

    AliGenerator & operator=(const AliGenerator &gen);

  // Getters

    virtual void GetOrigin(Float_t &ox, Float_t &oy, Float_t &oz) const
	{ox=fOrigin.At(0);oy=fOrigin.At(1);oz=fOrigin.At(2);}
    virtual void GetOrigin(TLorentzVector &o) const
	{o[0]=fOrigin.At(0);o[1]=fOrigin.At(1);o[2]=fOrigin.At(2);o[3]=0;}

protected:
    static TGenerator* fgMCEvGen; // Pointer to the generator
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
    Int_t       fAnalog;       //Flaf for anolog or pt-weighted generation
   //
    VertexSmear_t     fVertexSmear;  //Vertex Smearing mode
    VertexSource_t    fVertexSource; //Vertex source (internal/external)    
    Int_t       fTrackIt;    // if 1 Track final state particles 
    TArrayF     fOrigin;     // Origin of event
    TArrayF     fOsigma;     // Sigma of the Origin of even
    TArrayF     fVertex;     //! Vertex of current event
    
    enum {kThetaRange=1, kVertexRange=2, kPhiRange=4, kPtRange=8,
	  kYRange=32, kMomentumRange=16};

    ClassDef(AliGenerator,1) // Base class for event generators
};

#endif
