#ifndef ALIGENZDC_H
#define ALIGENZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////
//                                                //  
// Test pc generator for ZDC (taking into account //
// Fermi smearing, beam divergence and crossing)  //
//                                                //
////////////////////////////////////////////////////


#include <TMath.h>
 
#include "AliGenerator.h"

 
class AliGenZDC : public AliGenerator {

public:
  AliGenZDC();
  AliGenZDC(Int_t npart);
  virtual      ~AliGenZDC() {}
  virtual void Init();
  virtual void Generate();
  
  // Parameters that could be set for generation
  virtual void SetParticle(Int_t ipart) {fIpart=ipart;};
  virtual void SetMomentum(Float_t ptot) {fPMin=ptot; fPMax=ptot;};
  virtual void SetDirection(Float_t zpsrp, Float_t cosx, Float_t cosy, Float_t cosz)
                {fPseudoRapidity=zpsrp; fCosx=cosx; fCosy=cosy; fCosz=cosz;};
  virtual void SetFermi(Int_t Fflag) {fFermiflag=Fflag;};
  virtual void SetDiv(Float_t bmdiv, Float_t bmcra, Int_t iflcr) 
                {fBeamDiv=bmdiv; fBeamCrossAngle=bmcra; fBeamCrossPlane=iflcr;};
  
  // Getters 
  Double_t GetFermi2p(Int_t key) {return fProbintp[key];}
  Double_t GetFermi2n(Int_t key) {return fProbintn[key];}
  Float_t  GetInMomentum(Int_t key) {return fPInit[key];};
  Float_t  GetBoostMomentum(Int_t key) {return fBoostP[key];};
  Float_t  GetDivMomentum(Int_t key) {return fDivP[key];};
  Float_t  GetTrackMomentum(Int_t key) {return fPTrack[key];};
  
  // Fermi smearing, beam divergence and crossing angle       	       
  virtual void FermiTwoGaussian(Double_t A, Float_t Z, Double_t* pp, 
                Double_t* probintp, Double_t* probintn);
  virtual void ExtractFermi(Int_t id, Double_t* pp, Double_t* probintp, 
                Double_t* probintn, Double_t* pFermi);
  virtual void BeamDivCross(Int_t icross, Float_t divergence, Float_t crossangle, 
                Int_t crossplane, Double_t* pLab);
  virtual void AddAngle(Double_t theta1, Double_t phi1, Double_t theta2,
  	        Double_t phi2, Double_t* angle);
 
protected:
  Int_t    fIpart;              // Particle to generate
  Float_t  fCosx;               // Cos x of particle
  Float_t  fCosy;               // Cos y of particle
  Float_t  fCosz;               // Cos z of particle
  Float_t  fPseudoRapidity;     // Pseudo Rapidity of particle
  Int_t    fFermiflag;          // Fermi momentum flag
  Float_t  fBeamDiv;            // Beam divergence
  Float_t  fBeamCrossAngle;     // Beam crossing angle
  Int_t    fBeamCrossPlane;     // Beam crossing plane
  Double_t fProbintp[201];      // for protons
  Double_t fProbintn[201];      // for neutrons
  Double_t fPp[201];            // for protons
  Double_t fP[3];               // temporary momentum
  Float_t  fPInit[3];           // initial momentum 
  Float_t  fBoostP[3];          // boosted momentum
  Float_t  fDivP[3];            // divergence
  Float_t  fPTrack[3];          // track momentum
  
   ClassDef(AliGenZDC,1)  // Generator for AliZDC class
};

#endif
