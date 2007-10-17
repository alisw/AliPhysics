#ifndef ALIANAGAMMAHADRON_H
#define ALIANAGAMMAHADRON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.4  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.3.4.2  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of gamma - hadron correlations
//*-- Author: Gustavo Conesa (INFN-LNF)

#include "AliAnaGammaCorrelation.h"

class AliAnaGammaHadron : public AliAnaGammaCorrelation {

public: 
  
  AliAnaGammaHadron() ; // default ctor
  AliAnaGammaHadron(const AliAnaGammaHadron & g) ; // cpy ctor
  AliAnaGammaHadron & operator = (const AliAnaGammaHadron & g) ;//cpy assignment
  virtual ~AliAnaGammaHadron() ; //virtual dtor

  TList * GetCreateOutputObjects();

  void InitParameters();

  void Print(const Option_t * opt) const;
 
  void MakeGammaCorrelation(TParticle *pGamma, TClonesArray * plCTS,   TClonesArray * plNe) ;
  void MakeGammaChargedCorrelation(TParticle *pGamma, TClonesArray * pl) ;
  void MakeGammaNeutralCorrelation(TParticle *pGamma, TClonesArray * pl)  ;

  private:
  
  //Histograms
  TH2F * fhPhiCharged  ; 
  TH2F * fhPhiNeutral   ; 
  TH2F * fhEtaCharged  ; 
  TH2F * fhEtaNeutral   ; 
  TH2F * fhDeltaPhiGammaCharged  ;  
  TH2F * fhDeltaPhiGammaNeutral   ; 
  TH2F * fhDeltaEtaGammaCharged  ; 
  TH2F * fhDeltaEtaGammaNeutral  ; 
  TH2F * fhDeltaPhiChargedPt  ;

  TH2F * fhCorrelationGammaNeutral  ; //Neutral hadron correlation histogram 
  TH2F * fhCorrelationGammaCharged  ; //Charged hadron correlation histogram
  
  ClassDef(AliAnaGammaHadron,1)
} ;
 

#endif //ALIANAGAMMAHADRON_H



