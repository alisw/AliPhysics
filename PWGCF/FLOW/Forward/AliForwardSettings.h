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
  //~AliForwardSettings();                                       // destructor
  //AliForwardSettings(const AliForwardSettings &);             // copy constructor
  //AliForwardSettings & operator=(const AliForwardSettings &); // assignment
  // Lower bound phi acceptance (should always be 0)
  Double_t fPhiAcceptanceLowEdge;
  // Upper bound phi acceptance (should always be 2pi)
  Double_t fPhiAcceptanceUpEdge;
  // Number of bins used along phi
  // fPhiBins must be divisable by 2, but not by 4; so that we can later shift it by pi/2 (2pi is total int.)
  // The idea is to have the deltaPhi histogram with a bin centered arround 0
  Int_t fNPhiBins;

  // Lower bound eta acceptance (should always be -6.0)
  Double_t fEtaLowEdge;
  // Upper bound eta acceptance (should always be 6.0)
  Double_t fEtaUpEdge;
  // Lower edge of the Zvtx acceptance region in cm
  Double_t fZVtxAcceptanceLowEdge;
  // Upper edge of the Z_vtx acceptance region in cm
  Double_t fZVtxAcceptanceUpEdge;
  // Number of bins used along Z_vtx
  Int_t fNZvtxBins;


  Int_t fnoSamples;
  Int_t fNRefEtaBins; // eta bins in reference histograms
  Int_t fNDiffEtaBins; // eta bins in differential histograms
  Int_t fCentBins; // bins in centrality
  Int_t fCentUpEdge; // up edge in centrality

  TH3F* nuacentral;
  TH3F* nuaforward;
  TH3F* seccorr_fwd;
  TH3F* seccorr_cent;
  TH3F* nuehist;
  bool doNUA;
  bool doNUE;

  Double_t gap;
  Double_t minpt;
  Double_t maxpt;
  Bool_t mc;
  Bool_t esd;

  Int_t tracktype;
  UShort_t nua_mode;
  UShort_t ref_mode;
  Bool_t useTPC;
  Bool_t useSPD;
  Bool_t useITS;
  Bool_t use_primaries_cen;
  Bool_t use_primaries_fwd;
  Bool_t use_primaries_fwdref;
  Bool_t useEventcuts;
  TString centrality_estimator;
  Bool_t etagap;
  Bool_t makeFakeHoles;
  Int_t fnoClusters;
  Double_t fCutChargedDCAxyMax;
  Double_t fCutChargedDCAzMax;
  Bool_t doPt;
  Bool_t stdQC;
  Bool_t sec_corr;
  Bool_t a5;
  TString fileName;
  Int_t fMaxConsequtiveStrips;
  Bool_t standard_only;
  Double_t fmdlowcut;
  Double_t fmdhighcut;
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
    kSPDref   = 0x0001, // Standard QC{2} and QC{4} calculations
    kITSref  = 0x0002, // QC{2} w/ an eta-gap
    kTPCref    = 0x0004, // 3 correlator method for QC{2} w/ an eta-gap
    kFMDref = 0x0008
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
    kWA2,              // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
    kWB,               // multiplicity for all particles in subevent B (note subevent B can NOT be the entire event)
    k3pWeight,         // M(M-1)(M-1) or (mp*M-2mq)(M-1)
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
    kSinphi1phi2phi3p, // <sin(phi1+phi2-phi3)>
  };


  // definition of different variables to save
  // enum {
  //   kWA = 1,           // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
  //   kWA2,           // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
  //   kWB,               // multiplicity for all particles in subevent B (note subevent B can NOT be the entire event)
  //   kW2,               // <w2>
  //   k3pWeight,         // M(M-1)(M-1) or (mp*M-2mq)(M-1)
  //   kW2Two,            // <w2*two>
  //   kW4,               // <w4>
  //   kW4Four,           // <w4*four>
  //   kCosphi1A,         // <cos(phi1)> for subevent A
  //   kSinphi1A,         // <sin(phi1)> for subevent A
  //   kCosphi1B,         // <cos(phi1)> for subevent B
  //   kSinphi1B,         // <sin(phi1)> for subevent B
  //   kCosphi1phi2p,     // <cos(phi1+phi2)>
  //   kCosphi1phi2m,     // <cos(phi1-phi2)>
  //   kSinphi1phi2p,     // <sin(phi1+phi2)>
  //   kCosphi1phi2phi3m, // <cos(phi1-phi2-phi3)>
  //   kSinphi1phi2phi3m, // <sin(phi1-phi2-phi3)>
  //   kCosphi1phi2phi3p, // <cos(phi1+phi2-phi3)>
  //   kSinphi1phi2phi3p,  // <sin(phi1+phi2-phi3)>
  // };

  // enum {
  //   kW2A =1,               // <w2>
  //   kW2B,               // <w2>
  //   kW2TwoA,            // <w2*two>
  //   kW2TwoB,            // <w2*two>
  //   kW4A,               // <w4>
  //   kW4B,               // <w4>
  //   kW4FourA,           // <w4*four>
  //   kW4FourB,           // <w4*four>
  //   kW4FourTwoA,
  //   kW4FourTwoB,
  //   kW4ThreeTwoA,
  //   kW4ThreeTwoB
  // };
  Int_t kCountBin = 0;
  Int_t kMBin = 1;
  Int_t kMeanBin = 2;

  Int_t dW2A         = 1; // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
  Int_t dW2TwoA      = 2; // <w2*two>
  Int_t dW2B         = 3; // multiplicity for all particles in subevent B (note subevent B can NOT be the entire event)
  Int_t dW2TwoB      = 4; // <w2*two>  Int_t kW4          = 3; // <w4>
  Int_t dW4          = 5;
  Int_t dW4Four      = 6;

  Int_t dW4FourTwo   = 1;
  Int_t dW4ThreeTwo  = 2;
  Int_t dW4_mixed    = 3;
  Int_t dWTwoTwoN    = 4; // Numerator of R_{n,n; 2}
  Int_t dWTwoTwoD    = 5; // Denominator of R_{n,n; 2}

  Int_t rW2          = 1; // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
  Int_t rW2Two       = 2; // <w2*two>
  Int_t rW4          = 3;
  Int_t rW4Four      = 4;


  Int_t  kW2 =1;               // <w2>
  Int_t  kW2Two=2;             // <w2*two>
  Int_t  kW4=3;               // <w4>
  Int_t  kW4Four=4;           // <w4*four>
  // Int_t  kW4FourTwo=5;
  // Int_t  kW4ThreeTwo=6;
  Int_t track_sample;
  // enum {
  //   kW2 =1,               // <w2>
  //   kW2Two,             // <w2*two>
  //   kW4,               // <w4>
  //   kW4Four,           // <w4*four>
  //   kW4FourTwo,
  //   kW4ThreeTwo
  // };
  // definition of different variables to save
  enum {
    kN2 = 1,
    kD2
  };

  Int_t nua_runnumber;
  TH3F* correct_nua_mc;
  Int_t run_list;
  
  Bool_t second_analysis;
  Bool_t SC_analysis;
  Bool_t decorr_analysis;
  Bool_t normal_analysis;

  Int_t runnumber;

private:
  ClassDef(AliForwardSettings, 1);
};
#endif
