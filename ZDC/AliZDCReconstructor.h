#ifndef ALIZDCRECONSTRUCTOR_H
#define ALIZDCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"

class TF1;
class AliLoader;


class AliZDCReconstructor: public AliReconstructor {
public:
  AliZDCReconstructor();
  virtual ~AliZDCReconstructor();

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         Reconstruct(AliRunLoader* runLoader, 
                                   AliRawReader* rawReader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;

private:
  AliZDCReconstructor(const AliZDCReconstructor& reconstructor);
  AliZDCReconstructor& operator = (const AliZDCReconstructor& reconstructor);

  void                 ReconstructEvent(AliLoader* loader, Int_t znraw,
                                        Int_t zpraw, Int_t zemraw) const;

  TF1*   fZNCen;     //! Nspectator n true vs. EZN
  TF1*   fZNPer;     //! Nspectator n true vs. EZN
  TF1*   fZPCen;     //! Nspectator p true vs. EZP
  TF1*   fZPPer;     //! Nspectator p true vs. EZP
  TF1*   fZDCCen;    //! Nspectators true vs. EZDC
  TF1*   fZDCPer;    //! Nspectators true vs. EZDC
  TF1*   fbCen;      //! b vs. EZDC
  TF1*   fbPer;      //! b vs. EZDC
  TF1*   fZEMn;      //! Nspectators n from ZEM energy
  TF1*   fZEMp;      //! Nspectators p from ZEM energy
  TF1*   fZEMsp;     //! Nspectators from ZEM energy
  TF1*   fZEMb;      //! b from ZEM energy

  ClassDef(AliZDCReconstructor, 0)   // class for the ZDC reconstruction
};

#endif
