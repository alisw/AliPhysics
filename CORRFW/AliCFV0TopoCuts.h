/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// Class to cut on V0 topology
//   -> support for :
//                    DCA between V0 daughters
//                    V0 daughter impact parameters wrt primary vertex
//                    cosine of V0 pointing angle
//
///////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
///////////////////////////////////////////////////////////////////////////

#ifndef ALICFV0TOPOCUTS_H
#define ALICFV0TOPOCUTS_H

#include "AliCFCutBase.h"

class TObject;
class AliVEvent;

class AliCFV0TopoCuts : public AliCFCutBase 
{
 public :
  AliCFV0TopoCuts () ;
  AliCFV0TopoCuts (const Char_t* name, const Char_t* title) ;
  AliCFV0TopoCuts (const AliCFV0TopoCuts& c) ;
  AliCFV0TopoCuts& operator=(const AliCFV0TopoCuts& c) ;
  virtual ~AliCFV0TopoCuts() { } ;
  Bool_t IsSelected(TObject* v0) ;
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  void   SetEvtInfo(TObject* evt) {fEvent = (AliVEvent*)evt;}
  void   SetMaxDcaDaughters (Double32_t dca)  {fMaxDcaDaughters = dca;}
  void   SetMinDcaNeg       (Double32_t dca)  {fMinDcaNeg = dca;}
  void   SetMinDcaPos       (Double32_t dca)  {fMinDcaPos = dca;}
  void   SetMinCosPointAngle(Double32_t cos)  {fMinCosP   = cos;}
  
 private :
  Double32_t   fMaxDcaDaughters ; // max. dca between V0 daughters
  Double32_t   fMinDcaNeg ;       // min impact parameter (aka dca to prim. vertex) of neg. daughter
  Double32_t   fMinDcaPos ;       // min impact parameter (aka dca to prim. vertex) of pos. daughter
  Double32_t   fMinCosP ;         // min cosine of pointing angle
  AliVEvent*   fEvent;            // pointer to current event (needed for cuts related to PV position)
  
  ClassDef(AliCFV0TopoCuts,0);
};

#endif
