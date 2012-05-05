//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionD0.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIOND0_H
#define ALIDXHFEPARTICLESELECTIOND0_H

#include "AliDxHFEParticleSelection.h"

/**
 * @class AliDxHFEParticleSelectionD0
 * D0 selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionD0 : public AliDxHFEParticleSelection {
  public:
  /// constructor
  AliDxHFEParticleSelectionD0(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionD0();

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual bool IsSelected(AliVParticle* p);

 protected:

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionD0(const AliDxHFEParticleSelectionD0&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionD0& operator=(const AliDxHFEParticleSelectionD0&);

  ClassDef(AliDxHFEParticleSelectionD0, 1);
};

#endif
