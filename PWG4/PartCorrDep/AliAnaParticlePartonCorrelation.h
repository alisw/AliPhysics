#ifndef ALIANAPARTICLEPARTONCORRELATION_H
#define ALIANAPARTICLEPARTONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of particle-parton correlation
// Particle (for example direct gamma) must be found in a previous analysis 
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT ---
class TH2F ;

// --- ANALYSIS ---
#include "AliAnaPartCorrBaseClass.h"

class AliAnaParticlePartonCorrelation : public AliAnaPartCorrBaseClass {
       
public:
  
  AliAnaParticlePartonCorrelation() ; // default ctor
  virtual ~AliAnaParticlePartonCorrelation() {;} //virtual dtor              
  
  TList * GetCreateOutputObjects();
  
  void InitParameters();
    
  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 
  
  void Print(const Option_t * opt) const;
  
private:
  
  TH2F * fhDeltaEtaNearParton; //! Difference of parton eta and prompt trigger particle eta
  TH2F * fhDeltaPhiNearParton; //! Difference of parton phi and prompt trigger particle phi
  TH2F * fhDeltaPtNearParton;  //! Difference of parton pT and prompt trigger particle pT
  TH2F * fhPtRatNearParton;    //! Ratio of parton pT and prompt trigger particle pT
  
  TH2F * fhDeltaEtaAwayParton; //! Difference of parton eta and prompt trigger particle eta
  TH2F * fhDeltaPhiAwayParton; //! Difference of parton phi and prompt trigger particle phi
  TH2F * fhDeltaPtAwayParton;  //! Difference of parton pT and prompt trigger particle pT
  TH2F * fhPtRatAwayParton;    //! Ratio of parton pT and prompt trigger particle pT
  
  AliAnaParticlePartonCorrelation(const AliAnaParticlePartonCorrelation & g) ;               // cpy ctor
  AliAnaParticlePartonCorrelation & operator = (const AliAnaParticlePartonCorrelation & g) ; // cpy assignment
  
  ClassDef(AliAnaParticlePartonCorrelation,1)
  
} ;


#endif //ALIANAPARTICLEPARTONCORRELATION_H



