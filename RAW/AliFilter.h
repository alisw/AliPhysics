#ifndef ALIFILTER_H
#define ALIFILTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// base class for high level filter algorithms
// Derived classes should implement a default constructor and
// the virtual method Filter
//

#include <TObject.h>

class AliRawEvent;
class AliESD;


class AliFilter: public TObject {
public:
  virtual Bool_t       Filter(AliRawEvent* event, AliESD* esd) = 0;

private:
  ClassDef(AliFilter, 0)   // base class for high level filter algorithms
};

#endif
