#ifndef ALIANAGAMMAJETFINDER_H
#define ALIANAGAMMAJETFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of gamma-jet (standard finder) correlation 
//-- Author: Gustavo Conesa (INFN-LNF)

#include "AliAnaGammaCorrelation.h"
     
class AliAnaGammaJetFinder : public AliAnaGammaCorrelation {
       
  public: 
       
       AliAnaGammaJetFinder() ; // default ctor
       AliAnaGammaJetFinder(const AliAnaGammaJetFinder & g) ; // cpy ctor
       AliAnaGammaJetFinder & operator = (const AliAnaGammaJetFinder & g) ;//cpy assignment
       virtual ~AliAnaGammaJetFinder() ; //virtual dtor
              
       TList * GetCreateOutputObjects();
       
       void InitParameters();
       
       void Print(const Option_t * opt) const;
       
       void MakeGammaCorrelation( TParticle * pGamma, TClonesArray *pl, TClonesArray *)  ;
       
  private:
       
       TH2F * fhDeltaEtaJet;
       TH2F * fhDeltaPhiJet;
       TH2F * fhDeltaPtJet;
       TH2F * fhPtRatJet;
       
       ClassDef(AliAnaGammaJetFinder,1)
 } ;


#endif //ALIANAGAMMAJETFINDER_H



