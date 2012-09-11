//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionMCD0.h
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  D0 MC selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIONMCD0_H
#define ALIDXHFEPARTICLESELECTIONMCD0_H

#include "AliDxHFEParticleSelectionD0.h"
#include "AliDxHFEToolsMC.h"

/**
 * @class AliDxHFEParticleSelectionMCD0
 * Monte Carlo D0 selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionMCD0 : public AliDxHFEParticleSelectionD0 {
  public:
  /// constructor
  AliDxHFEParticleSelectionMCD0(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionMCD0();

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual int IsSelected(AliVParticle* p, const AliVEvent *pEvent=NULL);

  virtual THnSparse* DefineTHnSparse() const;
  // TODO: function can be renamed to better describe what it's doing
  virtual int DefineParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;

  /// check MC criteria
  int CheckMC(AliVParticle* p, const AliVEvent* pEvent);

  /// clear internal memory
  virtual void Clear(const char* option="");

 protected:

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionMCD0(const AliDxHFEParticleSelectionMCD0&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionMCD0& operator=(const AliDxHFEParticleSelectionMCD0&);

  AliDxHFEToolsMC fMCTools;  // MC selction tools
  int fResultMC;             // Result on MC check
  int fOriginMother;         // Holds info on the original mother particle
  static const char* fgkTrackControlBinNames[]; //! bin labels for track control histogram


  ClassDef(AliDxHFEParticleSelectionMCD0, 1);
};

#endif
