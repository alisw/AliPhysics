#ifndef ALIPHOSEVENTSELECTION_CXX
#define ALIPHOSEVENTSELECTION_CXX

 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class for selection of PHOS events, create a sub-class inheriting from
// this class to customize your event selection logback using
// AliPHOSClusterSelectionLogbackTask
//
// Uses IsEqual() to identify scope of Event selection!!! (Map from
// EventSelection to Event cluster array)
//
// Authors : Henrik Qvigstad
// Date    : 16.01.2014
/* $Id$ */


class AliVEvent;
class AliESDCaloEvent;
class AliAODCaloEvent;

#include "TObject.h"

class AliPHOSEventSelection : public TObject {
 public:
  AliPHOSEventSelection();
  virtual ~AliPHOSEventSelection();

  virtual Bool_t  IsEqual(const TObject* obj) const; // recommended: reimplement
  virtual ULong_t Hash() const;


 protected:
  AliPHOSEventSelection(const AliPHOSEventSelection&); // not implemented
  AliPHOSEventSelection& operator=(const AliPHOSEventSelection&); // not implemented

  ClassDef(AliPHOSEventSelection, 1);
};

#endif
