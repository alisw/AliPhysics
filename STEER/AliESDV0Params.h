#ifndef ALIESDV0PARAMS_H
#define ALIESDV0PARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class - parameterization 
//          This class is part of the Event Summary Data set of classes
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"

class AliESDV0Params : public TObject{
  friend class AliESDv0;
 public:
  AliESDV0Params();
 private:
  Double_t  fPSigmaOffsetD0;        // sigma offset DCA
  Double_t  fPSigmaOffsetAP0;       // sigma offset AP
  // effective sigma DCA params    
  Double_t  fPSigmaMaxDE;           // maximal allowed sigma DCA
  Double_t  fPSigmaOffsetDE;        // offset sigma DCA
  Double_t  fPSigmaCoefDE;          // sigma coefiecient 
  Double_t  fPSigmaRminDE;          // max radius  - with momentum dependence 
  // effective sigma PA params
  Double_t  fPSigmaBase0APE;        // base sigma PA
  Double_t  fPSigmaMaxAPE;          // maximal sigma PA
  Double_t  fPSigmaR0APE;           // radial dependent part   - coeficient
  Double_t  fPSigmaR1APE;           // radial dependent part   - offset 
  Double_t  fPSigmaP0APE;           // momentum dependent part - coeficient
  Double_t  fPSigmaP1APE;           // momentum dependent part - offset
  // minimax parameters
  Double_t fPMinFractionAP0;        // minimal allowed fraction of effective params - PA
  Double_t fPMaxFractionAP0;        // maximal allowed fraction of effective params - PA
  Double_t fPMinAP0;                // minimal minimax - PA sigma 
  //
  Double_t fPMinFractionD0;         // minimal allowed fraction of effective params - DCA
  Double_t fPMaxFractionD0;         // maximal allowed fraction of effective params - DCA
  Double_t fPMinD0;                 // minimal minimax - DCA sigma
  //
  ClassDef(AliESDV0Params,1)      // ESD V0 vertex - error and likelihood parameterization constant
};



#endif
