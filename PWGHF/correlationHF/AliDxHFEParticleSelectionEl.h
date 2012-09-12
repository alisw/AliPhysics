//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionEl.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  Electron selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIONEL_H
#define ALIDXHFEPARTICLESELECTIONEL_H

#include "AliDxHFEParticleSelection.h"
class AliPID;
class AliPIDResponse;
class AliHFEcuts;
class AliHFEvarManager;
class AliHFEpid;
class AliHFEpidBase;
class AliHFEtools;
class AliVEvent;
class AliCFManager;
class TH1;

/**
 * @class AliDxHFEParticleSelectionEl
 * Electron selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionEl : public AliDxHFEParticleSelection {
  public:
  /// constructor
  AliDxHFEParticleSelectionEl(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionEl();

  enum {
    kCutHFE = 0,
    kCutPID = 1,
    kNCuts
  };

  enum {
    kRecKineITSTPC=0,
    kRecPrim,
    kHFEcutsITS,
    kHFEcutsTOF,
    kHFEcutsTPC,
    kPID,
    kSelected,
    kNCutLabels
  };

  ///overloaded from AliDxHFEParticleSelection: Init
  virtual int Init();

  /// overloaded from AliDxHFEParticleSelection: init the control objects
  virtual int InitControlObjects();
  virtual THnSparse* DefineTHnSparse() const;

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual int IsSelected(AliVParticle* p, const AliVEvent* pEvent);

  virtual int HistogramParticleProperties(AliVParticle* p, int selected);

  // TODO: function can be renamed to better describe what it's doing
  virtual int DefineParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;

  /// set cuts object: a type cast check is implemented in the method
  virtual void SetCuts(TObject* /*cuts*/, int /*level*/=0);

 protected:

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionEl(const AliDxHFEParticleSelectionEl&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionEl& operator=(const AliDxHFEParticleSelectionEl&);

  /// check cut of specified step, e.g.
  bool ProcessCutStep(Int_t cutStep, AliVParticle *track);

  AliHFEpid*    fPID;                //! the PID object
  THnSparse*    fElectronProperties; // the particle properties of selected particles
  TH1*          fWhichCut;           // effective cut for a rejected particle
  AliHFEcuts*   fCuts;               //! Cuts for HF electrons
  AliCFManager* fCFM;                //! Correction Framework Manager

  static const char* fgkTrackControlBinNames[]; //! bin labels for track control histogram
  static const char* fgkCutBinNames[]; //! bin labels for cuts histogram

  ClassDef(AliDxHFEParticleSelectionEl, 2);
};

#endif
