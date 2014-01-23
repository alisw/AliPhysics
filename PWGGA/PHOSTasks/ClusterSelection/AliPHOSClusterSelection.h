#ifndef ALIPHOSCLUSTERSELECTION_CXX
#define ALIPHOSCLUSTERSELECTION_CXX

 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class for selection of PHOS clusters
// Authors : Henrik Qvigstad
// Date    : 16.01.2014
/* $Id$ */


class AliVCluster;
class AliESDCaloCluster;
class AliAODCaloCluster;
class AliVEvent;

#include "TObject.h"

class AliPHOSClusterSelection : TObject {
 public:
  AliPHOSClusterSelection();
  virtual ~AliPHOSClusterSelection();
  
  // Selection functions
  virtual Bool_t IsSelected(AliVCluster* cluster) const;

  virtual Bool_t IsSelectedCPV(AliVCluster* cluster) const;
  virtual Bool_t IsSelectedUnfolded(AliVCluster* cluster) const;
  virtual Bool_t IsSelectedDisp(AliVCluster* cluster) const;
  virtual Bool_t IsSelectedDispCore(AliVCluster* cluster) const;
  virtual Bool_t IsSelectedTOF(AliVCluster* cluster) const;

  // Configuration Functions
  AliPHOSClusterSelection* SetMinChargedParticleTrackDistance(Float_t distance);
  AliPHOSClusterSelection* SetNotUnfolded(Bool_t notUnfolded);
  AliPHOSClusterSelection* SetMaxDispR2(Float_t maxR2);
  AliPHOSClusterSelection* SetMaxDispCoreR2(Float_t maxR2);
  AliPHOSClusterSelection* SetMaxTOF(Float_t maxTOF);
  
  AliPHOSClusterSelection* SetMinSelection();

  virtual TString ToString() const;
  static Float_t SetMinChargedParticleTrackDistance(const TString string);
  
 protected:
  AliPHOSClusterSelection(const AliPHOSClusterSelection&); // not implemented
  AliPHOSClusterSelection& operator=(const AliPHOSClusterSelection&); // not implemented
  
  // Selection Parameters
  Float_t fMinChargedParticleTrackDistance; // CPV, Charged Particle Veto
  Bool_t fNotUnfolded; // if true, rejects Unfolded Clusters
  Float_t fMaxDispR2; // dispersion cut
  Float_t fMaxDispCoreR2; // dispersion cut of core cells
  Float_t fMaxTOF; // TOF cut

  AliVEvent* GetCurrentEvent() const;
  
  ClassDef(AliPHOSClusterSelection, 1);
};

#endif
