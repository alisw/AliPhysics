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
  void Init();
  void Generate();
  
  // Fermi smearing, beam divergence and crossing angle       	       
  void FermiTwoGaussian(Float_t A, Double_t *pp, 
        Double_t *probintp, Double_t *probintn);
  void ExtractFermi(Int_t id, Double_t *pp, Double_t *probintp, 
        Double_t *probintn, Double_t *pFermi);
  void BeamDivCross(Int_t icross, Float_t divergence, Float_t crossangle, 
        Int_t crossplane, Double_t *pLab);
  void AddAngle(Double_t theta1, Double_t phi1, Double_t theta2,
  	        Double_t phi2, Double_t *angle);
 
  
  // Parameters that could be set for generation
  void SetParticle(Int_t ipart) {fIpart=ipart;};
  void SetMomentum(Float_t ptot) {fPMin=ptot; fPMax=ptot;};
  void SetDirection(Float_t zpsrp, Float_t cosx, Float_t cosy, Float_t cosz)
          {fPseudoRapidity=zpsrp; fCosx=cosx; fCosy=cosy; fCosz=cosz;};
  void SetFermi(Int_t Fflag) {fFermiflag=Fflag;};
  void SetDiv(Float_t bmdiv, Float_t bmcra, Int_t iflcr) 
          {fBeamDiv=bmdiv; fBeamCrossAngle=bmcra; fBeamCrossPlane=iflcr;};
  void SetDebug(Int_t idebu) {fDebugOpt = idebu;};
  
  // Getters 
  Double_t GetFermi2p(Int_t key) {return fProbintp[key];}
  Double_t GetFermi2n(Int_t key) {return fProbintn[key];}

protected:
  Int_t    fIpart;              // Particle to be generated
  Float_t  fCosx;               // Director cos of the track - x direction
  Float_t  fCosy;               // Director cos of the track - y direction 
  Float_t  fCosz;               // Director cos of the track - z direction
  Float_t  fPseudoRapidity;     // Pseudorapidity (!=0 -> eta of the particle)
                                // (=0 -> director cos of the track)
  Int_t    fFermiflag;          // Fermi momentum flag (=1 -> Fermi smearing)
  Float_t  fBeamDiv;            // Beam divergence (angle in rad)
  Float_t  fBeamCrossAngle;     // Beam crossing angle (angle in rad)
  Int_t    fBeamCrossPlane;     // Beam crossing plane 
                                // (=1 -> horizontal, =2 -> vertical plane)
  Double_t fProbintp[201];      // Protons momentum distribution due to Fermi 
  Double_t fProbintn[201];      // Neutrons momentum distribution due to Fermi 
  Double_t fPp[201];            // Spectator momenta
  Int_t    fDebugOpt;		// Option for debugging [0->No debug, 1->Screen
  				//  prints, 2->ASCII data file]
  
 private:
  AliGenZDC(const AliGenZDC & gen);
  AliGenZDC & operator=(const AliGenZDC & gen);

   ClassDef(AliGenZDC,1)  	// Generator for AliZDC class
};

#endif
