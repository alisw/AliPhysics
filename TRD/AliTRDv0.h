#ifndef ALITRDV0_H
#define ALITRDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Manager and hits classes for set: TRD version 0                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
 
#include "AliTRD.h"

class AliTRDsim;

class AliTRDv0 : public AliTRD {

 public:

  AliTRDv0();
  AliTRDv0(const char *name, const char *title);
  virtual ~AliTRDv0();

  virtual Int_t    IsVersion() const           { return 0;    };
  virtual void     Init();

  virtual void     CreateGeometry();
  virtual void     CreateMaterials();

  virtual void     StepManager();
 
  virtual void     SetHits()                   { fHitsOn = 1; };
          void     SetTR(Bool_t )              { };

          Bool_t   GetTR() const               { return 0;    };

 protected:

          Int_t    fHitsOn;     //  Used to switch hits on

  ClassDef(AliTRDv0,2)          //  Transition Radiation Detector version 0 (fast simulator)

};

#endif
