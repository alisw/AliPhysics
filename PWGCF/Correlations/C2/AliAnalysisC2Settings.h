/// \class AliAnalysisC2Settings
/// \brief Encapsulation of settings and their setters
///
/// \author Christian Bourjau <cbourjau@cern.ch>, University of Copenhagen, Denmark
/// \date Apr 04, 2016

#ifndef AliAnalysisC2Settings_cxx
#define AliAnalysisC2Settings_cxx

#include "TMath.h"
#include "TObject.h"
#include <vector>

typedef std::vector< Double_t > edgeContainer;

class AliAnalysisC2Settings : public TObject {
 public:
  AliAnalysisC2Settings();

  // Types of data this analysis can access
  static const Int_t kMCTRUTH = 0;
  static const Int_t kRECON = 1;
  // The type of data this task is accessing
  Int_t fDataType;
  // The lower bound of the eta acceptance region
  Double_t fEtaAcceptanceLowEdge;
  // The upper bound of the eta acceptance region
  Double_t fEtaAcceptanceUpEdge;
  // Number of bins used along eta
  Int_t fNEtaBins;
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
  Double_t fMaxDcaTang;
  // Max longitudinal distance of closest approach to primary vertex in order to be still valid
  Double_t fMaxDcaLong;

 private:
  ClassDef(AliAnalysisC2Settings, 1);
};
#endif
