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
  virtual Bool_t IsSelectedTOF(AliVCluster* cluster) const;

  // Configuration Functions
  AliPHOSClusterSelection* SetMinChargedParticleTrackDistance(Float_t distance);
  AliPHOSClusterSelection* SetNotUnfolded(Bool_t notUnfolded);
  AliPHOSClusterSelection* SetMaxDispR2(Float_t maxR2);
  AliPHOSClusterSelection* SetCoreRadius(Float_t radius);
  AliPHOSClusterSelection* SetMaxTOF(Float_t maxTOF);
  
  AliPHOSClusterSelection* SetMinSelection();

  virtual TString ToString() const;
  static Float_t  GetMinChargedParticleTrackDistance(const TString string);
  static Bool_t   GetUnfolded(const TString string);
  static Float_t  GetMaxDispR2(const TString string);
  static Float_t  GetCoreRadius(const TString string);
  static Float_t  GetMaxTOF(const TString string);


 protected:
  AliPHOSClusterSelection(const AliPHOSClusterSelection&); // not implemented
  AliPHOSClusterSelection& operator=(const AliPHOSClusterSelection&); // not implemented
  
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge, Double_t mf);//Copied from Pi0Flow task
  void  AliPHOSClusterSelection::EvalCoreLambdas(AliVCluster * clu, AliVCaloCells * cells,Double_t &m02, Double_t &m20);//Copied from Pi0Flow task
  Bool_t AliAnalysisTaskPi0Flow::TestLambda(Double_t pt,Double_t l1,Double_t l2);//Copied from Pi0FlowTask
  // Selection Parameters
  Float_t fMinChargedParticleTrackDistance; // CPV, Charged Particle Veto
  Bool_t fNotUnfolded; // if true, rejects Unfolded Clusters
  Float_t fMaxDispR2; // dispersion cut
  Float_t fCoreRadius; // Core radius
  Float_t fMaxTOF; // TOF cut

  AliVEvent* GetCurrentEvent() const;

  ClassDef(AliPHOSClusterSelection, 1);
};

#endif

