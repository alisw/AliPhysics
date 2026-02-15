#ifndef AliAnalysisTaskFM_moments_H
#define AliAnalysisTaskFM_moments_H

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliEventPoolManager.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TNtuple.h"

const Int_t _numM = 52, _maxPtBin = 5, _numEnvs = 6;

class TH1F;
class TH2F;
class TNtuple;
class TString;
class TList;
class AliAODEvent;
class AliMCEvent;
class AliPIDResponse;
class AliEventPoolManager;

class AliAnalysisTaskFM_moments : public AliAnalysisTaskSE {
public:
AliAnalysisTaskFM_moments();
AliAnalysisTaskFM_moments(const char *name);
virtual ~AliAnalysisTaskFM_moments();
virtual void UserCreateOutputObjects();
virtual void UserExec(Option_t *option);
virtual void Terminate(Option_t *option);

void SetTwoTrackCuts(Double_t deta, Double_t dphi, Bool_t twotrack, Bool_t QA) { fdeta = deta; fdphi = dphi; flag2Track = twotrack; _FLAG_DETADPHI_ = QA; }
void SetSharingFraction(Double_t fshfr, Bool_t sharity) { fSharedFraction = fshfr; flagSharity = sharity; }
void SetVtxCut(Double_t fVx, Double_t fVy, Double_t fVz) { _vxMax = fVx; _vyMax = fVy; _vzMax = fVz; }
void SetIsMC(Bool_t isMC, TString gen) { _FLAG_MC_ = isMC; fGenName = gen; }
void SetMixedEvent(Bool_t mixed) { _FLAG_MIXED_ = mixed; }
void SetCentLim(Int_t f_minCent, Int_t f_maxCent) { _minCent = f_minCent; _maxCent = f_maxCent; }
void SetEta(Double_t f_minEta, Double_t f_maxEta) { _minEta = f_minEta; _maxEta = f_maxEta; }
void Setfbit(Int_t fb) { filterBit = fb; }
void SetYear(Int_t year) { _YEAR_ = year; }
void SetPSbinning(TArrayI MB, TArrayI NB, Int_t Mm) { _mBin2 = MB; _mNBin2 = NB; _maxPhaseSpaceBins = Mm; }
void SetPtArray(TArrayD f_ptArray, Int_t nbins) { _ptArray = f_ptArray; _numPtBin = nbins; }
void SetRejectElectrons(Bool_t fRejectElectron) { flagRejEls = fRejectElectron; }
void SetRejectPileup(Bool_t RejectPileup) { _FLAG_PILEUP_ = RejectPileup; }
void SetBfield(Int_t bf) { fBfield = bf; }
void SetDCAXYRangeMax(Int_t dcaxy) { _dcaXYMax = dcaxy; }
void SetDCAZRangeMax(Double_t dcaz) { _dcaZMax = dcaz; }
void SetITSClusterCut(Double_t itscut) { fITSCls = itscut; }
void SetTPCClusterCut(Double_t tpccut) { fTPCCls = tpccut; }
void SetnTPCrossedrows(Double_t tpcrowcut) { fTPCRows = tpcrowcut; }
void SetSelfAffine(Bool_t selfaffine) { flagSelfAff = selfaffine; }
void SetSharedCls(Double_t fSharedCls, Double_t fSharedrows, Double_t fFindable, Bool_t fSharedMean) { _sharedClsMax = fSharedCls; _sharedRowsMax = fSharedrows; _findableClsMin = fFindable; flagShClPro = fSharedMean; }

protected:
  AliEventCuts fEventCuts;
  void FillTrackInfo();
  void FillMCTrackInfo();
  void SetEventMixing(Bool_t mixed);
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
  Float_t CalculateDPhiStar(Float_t phi1, Float_t eta1, Float_t pt1, Int_t charge1, Float_t phi2, Float_t eta2, Float_t pt2, Int_t charge2,Float_t bSign);
  Float_t SharedClusterFraction(TBits &, TBits &, TBits &, TBits &);
  void GetPtBin(Double_t);
  void CalculateNFMs(Bool_t _isGen);
  Bool_t GetDCA(AliAODTrack *track, Double_t dca[2]);
  void DataPosting();

private:

// Event related (Data and MC)
AliAODEvent *fAOD;    // AOD object to access event-level data.
AliMCEvent *fMCEvent; // MC event object to access Monte Carlo truth information.
AliEventPool* _pool; // Event pool object to access event-level data.
Bool_t _FLAG_MC_;        // Flag to indicate whether the analysis is running on MC data.
AliEventPoolManager *fEvPoolMgr; // Event pool manager object to access event-level data.
Bool_t _FLAG_PILEUP_;    // Flag to enable or disable pile-up rejection.
Int_t _YEAR_;          // Year of data taking (e.g., 2010, 2015, 2018).
TString fGenName;     // Name of the generator used for MC events (e.g., "Hijing").
Int_t counter;        // General-purpose counter for internal use.
Double_t fdeta;       // Maximum allowed difference in pseudorapidity (Δη) for two-track cuts.
Double_t fdphi;       // Maximum allowed difference in azimuthal angle (Δφ) for two-track cuts.
Double_t fSharedFraction; // Maximum allowed fraction of shared clusters between tracks.
Bool_t flag2Track;    // Flag to enable or disable two-track cuts.
Bool_t flagSharity;   // Flag to enable or disable shared cluster analysis.
Bool_t flagRejEls;    // Flag to enable or disable electron rejection.
Bool_t flagSelfAff;   // Flag to enable or disable self-affine analysis.
Double_t fITSCls;     // ITS cluster cut value (minimum number of ITS clusters required).
Double_t fTPCCls;     // TPC cluster cut value (minimum number of TPC clusters required).
Double_t fTPCRows;    // Minimum number of TPC crossed rows required.
Bool_t flagShClPro;   // Flag to enable or disable shared cluster profiling.
Int_t fBfield;        // Magnetic field polarity: -1 for negative, 1 for positive, 0 for no cut.
Bool_t _FLAG_DETADPHI_;  // Flag to enable or disable QA for two-track cuts.
Bool_t _FLAG_MIXED_;    // Flag to enable or disable mixed event analysis.
Int_t _numPtBin;        // Number of transverse momentum (pT) bins for the analysis.

TList *hList[_numEnvs];  // Output lists for histograms and QA plots for different environments.
TList *fNtupleList[_maxPtBin][_numEnvs]; // Output lists for TNtuples for reconstructed data or MC.
TList *fNtupleListGen[_maxPtBin];     // Output lists for TNtuples for generated MC data.

AliPIDResponse *fPIDResponse; // PID response object for particle identification.
AliAODMCHeader *mcHeader;     // MC header object to access event-level MC information.

Int_t filterBit;              // Filter bit for track selection (e.g., 128 for TPC-only tracks).

Double_t _vxMax;              // Maximum allowed x-coordinate of the primary vertex.
Double_t _vyMax;              // Maximum allowed y-coordinate of the primary vertex.
Double_t _vzMax;              // Maximum allowed z-coordinate of the primary vertex.
Double_t _minCent;             // Minimum centrality percentile for event selection.
Double_t _maxCent;             // Maximum centrality percentile for event selection.
Double_t _minEta;              // Minimum pseudorapidity for track selection.
Double_t _maxEta;              // Maximum pseudorapidity for track selection.
Double_t _sharedClsMax;       // Maximum allowed number of shared clusters.
Double_t _sharedRowsMax;      // Maximum allowed number of shared rows.
Double_t _findableClsMin;     // Minimum fraction of findable clusters required.
Int_t _dcaXYMax;              // Maximum allowed DCA in the xy-plane.
Double_t _dcaZMax;            // Maximum allowed DCA in the z-direction.
Int_t foundPart[_maxPtBin][_numEnvs]; // Number of found particles in each pt bin and environment
Int_t _maxPhaseSpaceBins;                   // Maximum number of phase space bins (e.g., 82 or 123).
TArrayD _ptArray;              // Array of pT bin edges.
TArrayI _mBin2, _mNBin2;       // Arrays for phase space binning (e.g., for self-affine analysis).

std::vector<Bool_t> _ptBinFlags;    // vector of Boolean for pt bin

std::vector<TBits> clusmap;   // vector of TBits for cluster map
std::vector<TBits> sharedmap; // vector of TBits for shared map

Double_t c_field; //! variable for magnetic field
Double_t c_percentile; //! variable for centrality percentile

// Histograms for event level
TH1F *hEventCounter;    //! Histogram to track events etc
TH1F *hTrackCounter;    //! Histogram to track tracks etc
TH1F *hTrackCounterGen; //! Histogram to track gen tracks etc
TH1F *hEnvsCounter;     //! Histogram to track different envs etc
TH1F *hQACent[2];   //! Histogram for centrality distribution raw and sel
TH1F *hQAVx, *hQAVy;        //! Histograms for Vx and Vy
TH1F *hQAVz[2];     //! Histograms for Vz raw and sel

// Histograms for DCA distributions
TH1F *hDCAxy[_numEnvs];       //! DCAxy distribution for different environments
TH1F *hDCAz[_numEnvs];        //! DCAz distribution for different environments
TH1F *hDCAxy_imp[_numEnvs];   //! DCAxy calculated with impact parameter
TH1F *hDCAz_imp[_numEnvs];    //! DCAz calculated with impact parameter
TH2F *hDCAxypT[_numEnvs];     //! DCAxy vs pT
TH2F *hDCAzpT[_numEnvs];      //! DCAz vs pT

// Histograms for ITS and TPC cluster distributions
TH1F *hnITScls[_numEnvs];     //! ITS cluster distribution
TH1F *hnITScls2[_numEnvs];    //! ITS cluster distribution (alternative)
TH1F *hnTPCcls[_numEnvs];     //! TPC cluster distribution
TH1F *hnTPCcls2[_numEnvs];    //! TPC cluster distribution (alternative)
TH1F *hnTPCcrossedrows[_numEnvs]; //! TPC crossed rows distribution
TH1F *hmissingCls[_numEnvs]; //! Missing clusters distribution
TH1F *htpcTgl[_numEnvs]; //! TPC track tangent distribution

// Histograms for shared cluster distributions
TH1F *hNShCls[_numEnvs];      //! TPC shared cluster distribution
TH1F *hNShClsNcls[_numEnvs];  //! TPC shared cluster fraction (sharedcls/ncls)
TH2F *hNShClsFra[_numEnvs];  //! TPC shared cluster fraction (sharedcls/ncls vs sharedcls/ncrrows)
TH1F *hNShClsNcrows[_numEnvs]; //! TPC shared cluster fraction (sharedcls/ncrrows)
TH2F *hNFoundClsFra[_numEnvs]; //! TPC found cluster fraction (sharedcls/ncls vs ncrrows/findablecls)
TH2F *hNShClsVsPt[_numEnvs];  //! TPC shared cluster fraction vs pT

// Histograms for findable cluster distributions
TH1F *hNFcls[_numEnvs];       //! TPC findable cluster distribution
TH1F *hNFindClsFra[_numEnvs]; //! TPC findable cluster fraction (fncrrows/findablecls)
TH2F *hNFindClsVsPt[_numEnvs]; //! TPC findable cluster fraction vs pT

// Histograms for TPC and ITS signals
TH1F *hTPCsignal[_numEnvs];   //! TPC signal distribution
TH1F *hTPCsignalN[_numEnvs]; //! TPC signal distribution (alternative)
TH2F *hTPCsignalvsPt[_numEnvs]; //! TPC signal vs pT
TH2F *hTPCsignalvsPtN[_numEnvs]; //! TPC signal vs pT (alternative)
TH2F *hTPCsignalvsPtTuned[_numEnvs]; //! TPC signal vs pT (tuned)
TH2F *hTPCsignalvsPtot[_numEnvs]; //! TPC signal vs Ptot
TH1F *hTPCsignalTuned[_numEnvs]; //! TPC signal distribution (tuned)
TH1F *hITSsignal[_numEnvs];   //! ITS signal distribution

// Histograms for chi2 distributions
TH1F *hChi2TPC[_numEnvs];     //! TPC chi2 distribution
TH1F *hChi2ITS[_numEnvs];     //! ITS chi2 distribution
TH1F *hChi2perclTPC[_numEnvs]; //! TPC chi2 per cluster distribution

// Histograms for kinematic distributions
TH1F *hPtDis[_numEnvs];       //! Pt distribution
TH1F *hPtotDis[_numEnvs];     //! Ptot distribution
TH1F *hEtaDis[_numEnvs];      //! Eta distribution
TH1F *hPhiDis[_numEnvs];      //! Phi distribution
TH1F *hPtDisGen;           //! Pt distribution for MC generated tracks
TH1F *hEtaDisGen;          //! Eta distribution for MC generated tracks
TH1F *hPhiDisGen;          //! Phi distribution for MC generated tracks

// Histograms for golden chi2
TH1F *hGoldenChi2[_numEnvs];  //! Golden Chi2 distribution

// Histograms for particle counts
TH1F *hNumPions[_numEnvs];    //! Number of Pions
TH2F *hNumPionsCent[_numEnvs]; //! Number of Pions vs Centrality
TH1F *hNumKaons[_numEnvs];    //! Number of Kaons
TH2F *hNumKaonsCent[_numEnvs]; //! Number of Kaons vs Centrality
TH1F *hNumProtons[_numEnvs];  //! Number of Protons
TH2F *hNumProtonsCent[_numEnvs]; //! Number of Protons vs Centrality

// Histograms for pt, eta, phi, and multiplicity bins
TH1F *hPtBin[_maxPtBin][_numEnvs];       //! Pt distribution for each pt bin and environment
TH1F *hEtaBin[_maxPtBin][_numEnvs];      //! Eta distribution for each pt bin and environment
TH1F *hPhiBin[_maxPtBin][_numEnvs];      //! Phi distribution for each pt bin and environment
TH1F *hMultBin[_maxPtBin][_numEnvs];     //! Multiplicity distribution for each pt bin and environment

// Histograms for MC generated tracks
TH1F *hPtBinGen[_maxPtBin];    //! Pt distribution for MC generated tracks
TH1F *hEtaBinGen[_maxPtBin];   //! Eta distribution for MC generated tracks
TH1F *hPhiBinGen[_maxPtBin];   //! Phi distribution for MC generated tracks
TH1F *hMultBinGen[_maxPtBin];  //! Multiplicity distribution for MC generated tracks

// 2D eta-phi histograms for phase space bins
TH2F *hEtaPhiBin[_maxPtBin][_numM][_numEnvs];       //! Eta-phi distribution for each pt bin, phase space bin, and environment
TH2F *hEtaPhiBinGen[_maxPtBin][_numM][1];    //! Eta-phi distribution for MC generated tracks, [1] is added as a dummy to avoid compilation in NFM calculation

// TNtuples for storing M2, average bin content, and factorial moments
TNtuple *fntpMBin[_maxPtBin][_numM][_numEnvs];      //! TNtuple for factorial moments for each pt bin, phase space bin, and environment
TNtuple *fntpMBinGen[_maxPtBin][_numM];   //! TNtuple for factorial moments for MC generated tracks

// Histograms for dEta and dPhi
TH1F *hdEta[_numEnvs];                     //! Histogram for dEta in different environments
TH1F *hdPhi[_numEnvs];                     //! Histogram for dPhi in different environments

// Histograms for dEta-dPhi for each pt bin
TH2F *hDEDPraw[_maxPtBin][_numEnvs];         //! Raw dEta-dPhi for each pt bin
TH2F *hDEDPSel[_maxPtBin][_numEnvs];         //! Selected dEta-dPhi for each pt bin

// Histograms for dEta-dPhi for same charges
TH2F *hDEDPrawChSame[_maxPtBin][2][3][_numEnvs]; //! dEta-dPhi for same charges (pt bins, range, charge)

// Histograms for dEta-dPhi for different charges
TH2F *hDEDPrawChDiff[_maxPtBin][2][_numEnvs];    //! dEta-dPhi for different charges (pt bins, range)

// Histograms for dEta-dPhi for pt-ordered tracks
TH2F *hDEDPrawPtOrd[_maxPtBin][4][3][2][_numEnvs]; //! dEta-dPhi for pt-ordered tracks (pt bins, order, charge, range)

AliAnalysisTaskFM_moments(const AliAnalysisTaskFM_moments &); // not implemented
AliAnalysisTaskFM_moments & operator=(const AliAnalysisTaskFM_moments &); // not implemented
ClassDef(AliAnalysisTaskFM_moments, 1);
};

#endif
