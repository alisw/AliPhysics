#ifndef ALIUNICORANALGLOBAL_H
#define ALIUNICORANALGLOBAL_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// event global variable analyzer
//=============================================================================

#include "AliUnicorAnal.h"
class AliUnicorEvent;

//=============================================================================
class AliUnicorAnalGlobal : public AliUnicorAnal {
   
 public:
  AliUnicorAnalGlobal(const char *nam="global"); // constructor
  virtual ~AliUnicorAnalGlobal(){}               // destructor
  void Process(AliUnicorEvent *ev) const;        // fill histograms

  ClassDef(AliUnicorAnalGlobal,1)
};
//=============================================================================
#endif
