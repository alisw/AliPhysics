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

class TH1;

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
  virtual int InitControlObjects();

  virtual THnSparse* DefineTHnSparse();
  virtual int FillParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;
  virtual AliVParticle* CreateParticle(AliVParticle* track);

  // Get the D0 efficiency
  virtual double GetD0Eff(AliVParticle* tr);

  /// check MC criteria
  int CheckMC(AliVParticle* p, const AliVEvent* pEvent);

  /// Flag to run over MC "stack". Not used at the moment
  void SetUseKine(bool kine){fUseKine=kine;}

  /// clear internal memory
  virtual void Clear(const char* option="");

  enum{
    kMCFirst=0,
    kMCLast=1,
    kMCOnly=2
  };

 protected:
  virtual int HistogramParticleProperties(AliVParticle* p, int selected=1);

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelectionMCD0(const AliDxHFEParticleSelectionMCD0&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionMCD0& operator=(const AliDxHFEParticleSelectionMCD0&);

  int ParseArguments(const char* arguments);

  AliDxHFEToolsMC fMCTools;  // MC selction tools
  TH1* fPDGnotMCD0;          // holds PDG of not MC truth D0s
  int fResultMC;             // Result on MC check
  int fOriginMother;         // Holds info on the original mother particle
  bool fUseKine;             // Whether to run over MC particles (true) or Reco (false)
  THnSparse* fD0PropertiesKine; //the particle properties of selected particles
  Bool_t fStoreOnlyMCD0;     // Store only MC truth D0s, discard the rest
  Int_t fMCInfo;             // What to check and in which order
  Bool_t fRequireD0toKpi;    // Only take D0 to Kpi
  Short_t fSystem;           // Which system pp/PbPb
  Short_t fUseCentrality;    // Using centrality or not
  ClassDef(AliDxHFEParticleSelectionMCD0, 4);
};

#endif
