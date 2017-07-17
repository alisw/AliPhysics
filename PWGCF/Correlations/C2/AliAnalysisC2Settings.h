/// \class AliAnalysisC2Settings
/// \brief Encapsulation of settings and their setters
///
/// \author Christian Bourjau <cbourjau@cern.ch>, University of Copenhagen, Denmark
/// \date Apr 04, 2016

#ifndef AliAnalysisC2Settings_cxx
#define AliAnalysisC2Settings_cxx

#include "TObject.h"
#include "TString.h"

class AliAnalysisC2Settings : public TObject {
  typedef std::vector< Double_t > edgeContainer;

 public:
  AliAnalysisC2Settings();

  // Types of data this analysis can access
  enum {kMCTRUTH, kRECON};
  // The type of data this task is accessing
  Int_t fDataType;
  // Vector holding the bin edges along eta in the backward region
  edgeContainer fEtaEdgesBwd;
  // Vector holding the bin edges along eta in the central region
  edgeContainer fEtaEdgesIts;
  // Vector holding the bin edges along eta in the forward region
  edgeContainer fEtaEdgesFwd;
  // Lower bound phi acceptance (should always be 0)
  Double_t fPhiAcceptanceLowEdge;
  // Upper bound phi acceptance (should always be 2pi)
  Double_t fPhiAcceptanceUpEdge;
  // Number of bins used along phi
  // fPhiBins must be divisable by 2, but not by 4; so that we can later shift it by pi/2 (2pi is total int.)
  // The idea is to have the deltaPhi histogram with a bin centered arround 0
  Int_t fNPhiBins;
  // Vector holding the bin edges for pT
  // Every histogram binning with a pT axis uses this setting
  edgeContainer fPtBinEdges;
  // Vector holding the bin edges in "centrality" (ie. multiplicity)
  edgeContainer fMultBinEdges;
  // Lower edge of the Zvtx acceptance region in cm
  Double_t fZVtxAcceptanceLowEdge;
  // Upper edge of the Z_vtx acceptance region in cm
  Double_t fZVtxAcceptanceUpEdge;
  // Number of bins used along Z_vtx
  Int_t fNZvtxBins;
  // Max tangential distance of closest approach to primary vertex in order to be still valid
  // Double_t fMaxDcaTang;
  // Max longitudinal distance of closest approach to primary vertex in order to be still valid
  // Double_t fMaxDcaLong;
  // Multiplicity estimator used for this analysis
  TString fMultEstimator;
  // Options for MultEstimator. This would be better down in a subclass, but CINT chokes on that...
  // Use these predefined estimator strings to set the estimator for this analysis like:
  //    fMultEstimator = fMultEstimatorV0M;
  TString fMultEstimatorV0M;
  TString fMultEstimatorRefMult08;
  TString fMultEstimatorValidTracks;

  // Require this trigger in the event selection. This is needed to
  // mitigate the missing physics selection in the ITS stand alone
  // runs (eg. 15j, at the time of writing)
  // TString fTrigger_str;
  // Predefined strings of triggers. Use these to set the trigger for this analysis like:
  //     fTrigger = fTriggerCint7;
  // TString fTriggerCint7;
  // TString fTriggerVhmV0M;

  // If availible, trigger selection should be done throught the physics selection trigger mass
  UInt_t fTriggerMask;
  // Use FMD in forward region; requires AliAODForwardMult object
  Bool_t fUseFMD;
  /// Use V0 amplitudes in the forward region
  Bool_t fUseV0;
  /// Use SPD Clusters in the central region
  Bool_t fUseSPDclusters;
  /// Use SPD tracklets in the central region
  Bool_t fUseSPDtracklets;

private:
  ClassDef(AliAnalysisC2Settings, 1);
};
#endif
