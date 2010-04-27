#ifndef ALITRDEVENTCUTS_H
#define ALITRDEVENTCUTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Event cut class for the TRD Performance Train                          //
//                                                                        //
// author                                                                 //
// Markus Fasel <m.fasel@gsi.de>                                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TObjArray;
class AliESDEvent;
class AliTRDeventCuts : public TNamed
{
public:
  AliTRDeventCuts();
  AliTRDeventCuts(const Char_t *name);
  ~AliTRDeventCuts();

  Bool_t IsSelected(AliESDEvent *event, Bool_t col=kTRUE);

  void AddTrigger(const Char_t *name);
  void SetVertexN(Int_t n);
  void SetVertexZ(Double_t z);

protected:
  Bool_t CheckTrigger(const Char_t *name);

private:
  AliTRDeventCuts(const AliTRDeventCuts &);
  AliTRDeventCuts &operator=(const AliTRDeventCuts &);

  TObjArray *fTriggerNames; // Container for Trigger names
  Int_t     fVertexN;       // Min number of contributors to Vertex
  Double_t  fVertexZ;       // Max Abs(z) of the reconstructed Vertex

  ClassDef(AliTRDeventCuts, 1)   //Event Cut class for TRD PWG1 Train 
};
#endif
