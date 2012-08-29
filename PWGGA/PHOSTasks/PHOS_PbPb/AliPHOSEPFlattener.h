#ifndef ALIPHOSEPFLATTENER_H
#define ALIPHOSEPFLATTENER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPHOSDigit.h 56876 2012-06-05 20:48:55Z fca $ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.34  2006/04/22 10:30:17  hristov
 * Add fEnergy to AliPHOSDigit and operate with EMC amplitude in energy units (Yu.Kharlov)
 *
 * Revision 1.33  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Class to perform flattening of the event plane distribution
//       It stores necessary parameterizations and apply when requested
//
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---
class TH1D ;
#include "TNamed.h"

class AliPHOSEPFlattener : public TNamed {

 public:
  
  AliPHOSEPFlattener() ;
  AliPHOSEPFlattener(const char * name) ; //To separate different runs use names
  AliPHOSEPFlattener(const AliPHOSEPFlattener & fl) ; 
  virtual ~AliPHOSEPFlattener() ;
  AliPHOSEPFlattener & operator = (const AliPHOSEPFlattener & flat);

public:

  Double_t MakeFlat(Double_t oldPhi,Double_t centrality)const ; //Apply (centrality-dependent) flattening to oldPhi
  void SetParameterization(TH1D * h) ;  //Set Parameterization to use (see code for the meaning of parameters

private:
  UShort_t fNCentrBins ; //Number of centrality bins
  UShort_t fNHarmonics ; //Number of harmonics used in parameterization
  Double32_t * fParam ;  //[-1,1.,16] [fNCentrBins*fNHarmonics]

  ClassDef(AliPHOSEPFlattener,1) 

} ;

#endif //  ALIPHOSEPFLATTENER_H
