#ifndef ALITRDV2_H
#define ALITRDV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set: TRD version 2   //
////////////////////////////////////////////////////////

#include "AliTRDv1.h"

//_____________________________________________________________________________
class AliTRDv2 : public AliTRDv1 {

 public:

  AliTRDv2();
  AliTRDv2(const char *name, const char *title);
  AliTRDv2(const AliTRDv2 &trd);
  virtual ~AliTRDv2();
  AliTRDv2 &operator=(const AliTRDv2 &trd);

  virtual void       Copy(TObject &trd) const;
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual Int_t      IsVersion() const          { return 2; };

 protected:

 private:
   
  ClassDef(AliTRDv2,1) // Transition Radiation Detector version 2 (slow simulator,detailed geometry)

};

#endif
