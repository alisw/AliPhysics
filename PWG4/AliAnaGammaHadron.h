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
  TH2F * fhPhiCharged  ; //Phi distribution of charged particles
  TH2F * fhPhiNeutral   ;  //Phi distribution of neutral particles
  TH2F * fhEtaCharged  ; //Eta distribution of charged particles
  TH2F * fhEtaNeutral   ; //Eta distribution of neutral particles
  TH2F * fhDeltaPhiGammaCharged  ;  //Difference of charged particle phi and prompt gamma phi as function of gamma pT
  TH2F * fhDeltaPhiGammaNeutral   ;  //Difference of neutral particle phi and prompt gamma phi as function of gamma pT
  TH2F * fhDeltaEtaGammaCharged  ;  //Difference of charged particle eta and prompt gamma eta as function of gamma pT
  TH2F * fhDeltaEtaGammaNeutral  ;  //Difference of neutral particle eta and prompt gamma eta as function of gamma pT
  TH2F * fhDeltaPhiChargedPt  ;  //Difference of charged particle eta and prompt gamma eta as function of charged pT

  TH2F * fhCorrelationGammaNeutral  ; //Neutral hadron correlation histogram 
  TH2F * fhCorrelationGammaCharged  ; //Charged hadron correlation histogram
  
  ClassDef(AliAnaGammaHadron,1)
} ;
 

#endif //ALIANAGAMMAHADRON_H



