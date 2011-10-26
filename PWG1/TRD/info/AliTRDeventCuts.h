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
  AliTRDeventCuts(const AliTRDeventCuts &ref);
  ~AliTRDeventCuts();

  void    AddTrigger(const Char_t *name);
  Bool_t  CheckTrigger(const Char_t *name);
  Bool_t  IsSelected(AliESDEvent *event, Bool_t col=kTRUE);

  void    Print(Option_t *opt="") const;

  void    SetVertexN(Int_t n) { fVertexN = n; };
  void    SetVertexZ(Double_t z) { fVertexZ = z; };
  void    SetBunchSelection(Int_t n, Int_t bunches[]);

private:
  AliTRDeventCuts &operator=(const AliTRDeventCuts &);

  TObjArray *fTriggerNames; // Container for Trigger names
  Int_t     *fBunches;      // List of bunches accepted for analysis
  Int_t     fVertexN;       // Min number of contributors to Vertex
  Double_t  fVertexZ;       // Max Abs(z) of the reconstructed Vertex

  ClassDef(AliTRDeventCuts, 1)   //Event Cut class for TRD PWG1 Train 
};
#endif
