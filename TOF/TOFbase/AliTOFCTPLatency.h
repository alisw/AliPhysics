#ifndef ALITOFCTPLATENCY_H
#define ALITOFCTPLATENCY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the CTPLatency object to be stored
// * in OCDB in order to apply CTPLatency correction during 
// * reconstruction. 
// *
// *
// *

#include "TObject.h"

class AliTOFCTPLatency :
public TObject
{

 public:

  AliTOFCTPLatency(); // default constructor
  virtual ~AliTOFCTPLatency(); // default destructor
  AliTOFCTPLatency(const AliTOFCTPLatency &source); // copy constructor
  AliTOFCTPLatency &operator=(const AliTOFCTPLatency &source); // operator=
  Float_t GetCTPLatency() const {return fCTPLatency;}; // getter
  void SetCTPLatency(Float_t value) {fCTPLatency = value;}; // setter

 private:

  Float_t fCTPLatency; // CTP latency (ps)

  ClassDef(AliTOFCTPLatency, 1);
};

#endif /* ALITOFCTPLATENCY_H */
