#ifndef ROOT_TMevSim
#define ROOT_TMevSim

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMevSim                                                              //
//                                                                      //
// This class implements an interface to the MevSim event generator.    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGenerator.h"
#include "MevSimCommon.h"
#include "TMevSimPartTypeParams.h"
#include "TObjArray.h"

class TMevSim : public TGenerator {

protected:

  Int_t fNEvents;
  Int_t fModelType;
  Int_t fReacPlaneCntrl;
  Float_t fPsiRMean, fPsiRStDev;
  Float_t fMultFacMean, fMultFacStDev;
  Float_t fPtCutMin, fPtCutMax;
  Float_t fEtaCutMin, fEtaCutMax;
  Float_t fPhiCutMin, fPhiCutMax;
  Float_t fNStDevMult, fNStDevTemp, fNStDevSigma, fNStDevExpVel, fNStdDevPSIr, fNStDevVn, fNStDevMultFac;
  Int_t fNIntegPts;
  Int_t fNScanPts;
  Int_t firand;  
  TClonesArray *fParticleTypeParameters;
   
// Copied from AliGeant
  enum {kMaxParticles = 35};

  Int_t fNPDGCodes;               // Number of PDG codes known by G3

  Int_t fPDGCode[kMaxParticles];  // Translation table of PDG codes

 public:   
   // Constructors and destructors
   
   TMevSim(Int_t nEvents = 1, Int_t modelType=1, Int_t reacPlaneCntrl=4,
	   Float_t psiRMean=0.0, Float_t psiRStDev=0.0, Float_t multFacMean=1.0, Float_t multFacStDev=0.05,
	   Float_t ptCutMin = 0.01, Float_t ptCutMax = 3.0, Float_t etaCutMin=(-4.5), Float_t etaCutMax = 4.5, 
	   Float_t phiCutMin=0.0, Float_t phiCutMax=360.0, Int_t irand=87266);

   TMevSim(TMevSim& mevsim);                    // copy constructor

   virtual            ~TMevSim();
   
   // Assignment operator
  
   virtual TMevSim& operator=(TMevSim& mevsim);

   // Mandatory TGenerator functions
   
   virtual void        Initialize();

   virtual void        GenerateEvent();

   virtual Int_t       ImportParticles(TClonesArray *particles, Option_t *option="");

   //Parameters for the generation:
   
   virtual void        SetNEvents(Int_t nEvents );
   virtual Int_t       GetNEvents() const;

   virtual Int_t       GetNPidTypes() const;
 
   virtual void        SetModelType(Int_t modelType);
   virtual Int_t       GetModelType() const;

   virtual void        SetReacPlaneCntrl(Int_t reacPlaneCntrl);
   virtual Int_t       GetReacPlaneCntrl() const;

   virtual void        SetPsiRParams(Float_t psiRMean,  Float_t psiRStDev);
   virtual Float_t     GetPsiRMean() const;
    virtual Float_t     GetPsiRStDev() const;

   virtual void        SetMultFacParams(Float_t multFacMean,  Float_t multFacStDev);
   virtual Float_t     GetMultFacMean() const;
   virtual Float_t     GetMultFacStDev() const;

   // Pt and geometry cut

   virtual void        SetPtCutRange(Float_t ptCutMin,  Float_t ptCutMax);
   virtual Float_t     GetPtCutMin() const;
   virtual Float_t     GetPtCutMax() const;

   virtual void        SetEtaCutRange(Float_t etaCutMin,  Float_t etaCutMax);
   virtual Float_t     GetEtaCutMin() const;
   virtual Float_t     GetEtaCutMax() const;

   virtual void        SetPhiCutRange(Float_t phiCutMin,  Float_t phiCutMax);
   virtual Float_t     GetPhiCutMin() const;
   virtual Float_t     GetPhiCutMax() const;

   // StDev

   virtual void        SetNStDevMult(Float_t nStDevMult);
   virtual Float_t       GetNStDevMult() const;

   virtual void        SetNStDevTemp(Float_t nStDevTemp);
   virtual Float_t       GetNStDevTemp() const;

   virtual void        SetNStDevSigma(Float_t nStDevSigma);
   virtual Float_t       GetNStDevSigma() const;

   virtual void        SetNStDevExpVel(Float_t nStDevExpVel);
   virtual Float_t       GetNStDevExpVel() const;

   virtual void        SetNStDevPSIr(Float_t nStDevPSIr);
   virtual Float_t       GetNStDevPSIr() const;

   virtual void        SetNStDevVn(Float_t nStDevVn);
   virtual Float_t       GetNStDevVn() const;

   virtual void        SetNStDevMultFac(Float_t nStDevMultFac);
   virtual Float_t       GetNStDevMultFac() const;

   // grid

   virtual void        SetNIntegPts(Int_t nIntegPts);
   virtual Int_t       GetNintegPts() const;

   virtual void        SetNScanPts(Int_t nScanPts);
   virtual Int_t       GetNScanPts() const;

   // adding particles

   virtual void        AddPartTypeParams(TMevSimPartTypeParams *params);
   virtual void        SetPartTypeParams(Int_t index, TMevSimPartTypeParams *params);
   virtual void        GetPartTypeParamsByIndex(Int_t index, TMevSimPartTypeParams *params);
   virtual void        GetPartTypeParamsByGPid(Int_t gpid, TMevSimPartTypeParams *params);

   // conversion   

   virtual Int_t       PDGFromId(Int_t gpid) const;
   virtual void        DefineParticles();
   

   ClassDef(TMevSim,1)  //Interface to MevSim Event Generator
     
};

#endif






