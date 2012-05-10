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

  /// overloaded from AliDxHFEParticleSelection: init the control objects
  virtual int InitControlObjects();

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual bool IsSelected(AliVParticle* p);

 protected:
  /// overloaded from AliDxHFEParticleSelection: histogram particle properties
  virtual int HistogramParticleProperties(AliVParticle* p, bool selected=true);

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionD0(const AliDxHFEParticleSelectionD0&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionD0& operator=(const AliDxHFEParticleSelectionD0&);

  THnSparse* fD0Properties; //! the particle properties of selected particles
  // TODO: at the moment the dimensions of the different THnSparse objects are different
  // needs to be consolidated
  // TODO: one might need particle properties of all and/or at different cut stages

  ClassDef(AliDxHFEParticleSelectionD0, 1);
};

#endif
