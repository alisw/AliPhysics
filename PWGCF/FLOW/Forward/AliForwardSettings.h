#ifndef AliForwardSettings_cxx
#define AliForwardSettings_cxx
/**
 * @file AliForwardSettings.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "TObject.h"
#include "TString.h"
#include "TH3.h"
#include "TH2.h"

class AliForwardSettings : public TObject {
  typedef std::vector< Double_t > edgeContainer;

 public:
  AliForwardSettings();

  // Types of data this analysis can access
  enum {kMCTRUTH, kRECON};
  // The type of data this task is accessing
  Int_t fDataType;

  // Lower bound phi acceptance (should always be 0)
  Double_t fPhiAcceptanceLowEdge;
  // Upper bound phi acceptance (should always be 2pi)
  Double_t fPhiAcceptanceUpEdge;
  // Lower bound eta acceptance (should always be -6.0)
  Double_t fEtaLowEdge;
  // Upper bound eta acceptance (should always be 6.0)
  Double_t fEtaUpEdge;
  // Number of bins used along phi
  // fPhiBins must be divisable by 2, but not by 4; so that we can later shift it by pi/2 (2pi is total int.)
  // The idea is to have the deltaPhi histogram with a bin centered arround 0
  Int_t fNPhiBins;
  // Lower edge of the Zvtx acceptance region in cm
  Double_t fZVtxAcceptanceLowEdge;
  // Upper edge of the Z_vtx acceptance region in cm
  Double_t fZVtxAcceptanceUpEdge;
  // Number of bins used along Z_vtx
  Int_t fNZvtxBins;

  // type of analysis
  TString qctype;

  Int_t fnoSamples;
  Int_t fNRefEtaBins; // eta bins in reference histograms
  Int_t fNDiffEtaBins; // eta bins in differential histograms
  Int_t fCentBins; // bins in centrality

  UShort_t fFlowFlags;     //  Flow flags, e.g., eta-gap, sat. vtx.

  TH3F* nuacentral;
  TH3F* nuaforward;

  bool doNUA;

  Double_t gap;
  Double_t maxpt;
  Double_t minpt;
  Bool_t mc;
  Bool_t esd;

  Int_t tracktype;
  UShort_t nua_mode;
  Bool_t useTPC;
  Bool_t useSPD;
  Bool_t use_primaries;
  Bool_t use_primaries_cen;
  Bool_t use_primaries_fwd;
  Bool_t etagap;
  Bool_t makeFakeHoles;
  Int_t fnoClusters;
  Double_t fCutChargedDCAxyMax;
  Double_t fCutChargedDCAzMax;
  TString centrality_estimator;
  // return true if good event

  // flags used for method of cumulant
  enum EFlowFlags {
    kStdQC   = 0x0001, // Standard QC{2} and QC{4} calculations
    kEtaGap  = 0x0002, // QC{2} w/ an eta-gap
    k3Cor    = 0x0004, // 3 correlator method for QC{2} w/ an eta-gap
    kSymEta  = 0x0008, // Symmetrize ref flow in std. QC{2} and QC{4} around eta = 0
    kSatVtx  = 0x0010, // Do satellite vertex input (currently not implemented)
    kNUAcorr = 0x0020, // Apply full NUA correction
    kFMD     = 0x0040, // Use FMD for forward flow
    kVZERO   = 0x0080, // Use VZERO for forward flow
    kSPD     = 0x0100, // SPD object flag
    kMC      = 0x0200, // MC object flag
    kTracks  = 0x1000, // Use tracks for reference flow
    kTPC     = 0x3000, // Use TPC tracks
  };

  // flags used for method of cumulant
  enum {
    kNormal   = 0x0001, // Standard QC{2} and QC{4} calculations
    kFill  = 0x0002, // QC{2} w/ an eta-gap
    kInterpolate    = 0x0004, // 3 correlator method for QC{2} w/ an eta-gap
  };

  enum {
    kTPCOnly = 128, // TPC only tracks
    kHybrid = 768, // TPC only tracks
    kGlobal = 32, // Global tracks
    kGlobalLoose = 64, // Global tracks
    kGlobalComb = 96,
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };
  // definition of different variables to save
  enum {
    kWA = 1,           // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
    kWA2,           // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
    kWB,               // multiplicity for all particles in subevent B (note subevent B can NOT be the entire event)
    kW2,               // <w2>
    k3pWeight,         // M(M-1)(M-1) or (mp*M-2mq)(M-1)
    kW4,               // <w4>
    kW4Four,           // <w4*four>
    kW2Two,            // <w2*two>
    kCosphi1A,         // <cos(phi1)> for subevent A
    kSinphi1A,         // <sin(phi1)> for subevent A
    kCosphi1B,         // <cos(phi1)> for subevent B
    kSinphi1B,         // <sin(phi1)> for subevent B
    kCosphi1phi2p,     // <cos(phi1+phi2)>
    kCosphi1phi2m,     // <cos(phi1-phi2)>
    kSinphi1phi2p,     // <sin(phi1+phi2)>
    kCosphi1phi2phi3m, // <cos(phi1-phi2-phi3)>
    kSinphi1phi2phi3m, // <sin(phi1-phi2-phi3)>
    kCosphi1phi2phi3p, // <cos(phi1+phi2-phi3)>
    kSinphi1phi2phi3p,  // <sin(phi1+phi2-phi3)>
  };

  // definition of different variables to save
  enum {
    kN2 = 1,
    kD2
  };

private:
  ClassDef(AliForwardSettings, 2);
};
#endif
