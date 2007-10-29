#ifndef ALIANAGAMMAPARTON_H
#define ALIANAGAMMAPARTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.3  2007/09/26 11:07:19  schutz
 * Update classes for the new analysis framwork
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of gamma-parton correlation 
//-- Author: Gustavo Conesa (INFN-LNF)

#include "AliAnaGammaCorrelation.h"
     
class AliAnaGammaParton : public AliAnaGammaCorrelation {
       
  public: 
       
       AliAnaGammaParton() ; // default ctor
       AliAnaGammaParton(const AliAnaGammaParton & g) ; // cpy ctor
       AliAnaGammaParton & operator = (const AliAnaGammaParton & g) ;//cpy assignment
       virtual ~AliAnaGammaParton() ; //virtual dtor
              
       TList * GetCreateOutputObjects();
       
       void InitParameters();
       
       void Print(const Option_t * opt) const;
       
       void MakeGammaCorrelation( TParticle * pGamma, TClonesArray *pl, TClonesArray *)  ;
       
  private:
       
       TH2F * fhDeltaEtaParton; //Difference of parton eta and prompt gamma eta
       TH2F * fhDeltaPhiParton; //Difference of parton phi and prompt gamma phi
       TH2F * fhDeltaPtParton; //Difference of parton pT and prompt gamma pT
       TH2F * fhPtRatParton; //Ratio of parton pT and prompt gamma pT
       
       ClassDef(AliAnaGammaParton,1)
 } ;


#endif //ALIANAGAMMAPARTON_H



