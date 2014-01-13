//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionMCEl.h
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  El MC selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIONMCEL_H
#define ALIDXHFEPARTICLESELECTIONMCEL_H

#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEToolsMC.h"

class THnSparse;
class TList;

/**
 * @class AliDxHFEParticleSelectionMCEl
 * Monte Carlo electron selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionMCEl : public AliDxHFEParticleSelectionEl {
  public:
  /// constructor
  AliDxHFEParticleSelectionMCEl(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionMCEl();

  //Overload Init function from electron selection class
  virtual int Init();

  // Setting up control objects: overloaded from AliDxHFEParticleSelection
  virtual THnSparse* DefineTHnSparse();

  virtual int FillParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;
  virtual AliVParticle* CreateParticle(AliVParticle* track);


  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual int IsSelected(AliVParticle* p, const AliVEvent *pEvent=NULL);

  /// check MC criteria
  int CheckMC(AliVParticle* p, const AliVEvent* pEvent);

  /// Return the result on check on initial quark
  int GetOriginMother() const {return fOriginMother;}

  /// clear internal memory
  virtual void Clear(const char* option="");

  enum{
    kAllPassingSelection=-1,
    kHadron=0,
    kNonHFE=1,
    kHFE=2,
    kOnlyc=3,
    kOnlyb=4,
  };

  enum{
    kMCFirst=0,
    kMCLast=1,
    kMCOnly=2
  };

  enum{
    kPi0HF=0,
    kPi0=1,
    kEtaHF=2,
    kEta=3,
    kGammaHF=4,
    kGamma=5,
    kJPsi=6,
    kOtherHF=7,
    kOther=8,
    kNrBackground=9
  };


 protected:

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionMCEl(const AliDxHFEParticleSelectionMCEl&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionMCEl& operator=(const AliDxHFEParticleSelectionMCEl&);

  int ParseArguments(const char* arguments);
  /// TODO: check if the label definitions can be used from the ToolsMC
  static const char* fgkPDGMotherBinLabels[];
  static const char* fgkPDGBinLabels[];
                    
  AliDxHFEToolsMC fMCTools;            // MC selection tools
  TList*          fHistoList;          // list of histograms
  int fOriginMother;                   //  Holds the origin motherquark (process)
  int fResultMC;                       // Holds information on check on MC
  Bool_t fUseKine;                     // Run over MC stack
  vector<int> fMotherPDGs;             // list off mothers counted as background
  Bool_t fUseMCReco;                   // Run over all MC reconstructed tracks
  Int_t fSelectionStep;                // Where to stop track selection
  Bool_t          fStoreCutStepInfo;   // Whether to store which cut step is the effective one
  Int_t fElSelection;                  // If there should be a selection on the electrons
  Bool_t fStoreOnlyMCElectrons;        // Store only MC truth electrons, discard hadrons
  Int_t fMCInfo;                       // What to check and in which order
  Bool_t fRemoveEfromD0;               // Removing electrons decaying from a D0
  Bool_t fRemoveSecondary;             // whether or not to remove secondary tracks
  Bool_t fUseGenerator;                // whether or not to check generator for particle
  Int_t fGenerator;                    // which generator was used to produce particle


  ClassDef(AliDxHFEParticleSelectionMCEl, 4);
};

#endif
