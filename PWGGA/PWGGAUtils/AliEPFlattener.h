#ifndef ALIEPFLATTENER_H
#define ALIEPFLATTENER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class to perform flattening of the event plane distribution
//       It stores necessary parameterizations and apply when requested
//
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---
class TH2 ;
#include "TNamed.h"

class AliEPFlattener : public TNamed {

 public:
  
  AliEPFlattener() ;
  AliEPFlattener(const char * name) ; //To separate different runs use names
  AliEPFlattener(const AliEPFlattener & fl) ; 
  virtual ~AliEPFlattener() ;
  AliEPFlattener & operator = (const AliEPFlattener & flat);

public:

  Double_t MakeFlat(Double_t oldPhi,Double_t centrality)const ; //Apply (centrality-dependent) flattening to oldPhi
  void SetParameterization(TH2 * h) ;  //Set Parameterization to use (see code for the meaning of parameters

private:
  Int_t fNCentrBins ; // Number of centrality bins
  Int_t fNHarmonics ; // Number of harmonics used in parameterization
  Int_t fNparam ;     // Total number of parameters (fNCentrBins*fNHarmonics)
  Double32_t *fParam ;  //[fNparam][-1.,1.,16] array of flattening parameters

  ClassDef(AliEPFlattener,1) 

} ;

#endif //  ALIEPFLATTENER_H
