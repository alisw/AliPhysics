#ifndef ALI2PCORRELATIONS_H
#define ALI2PCORRELATIONS_H

/// \file Ali2PCorrelations.h
/// \brief Class for collecting data for two-particle correlation functions
///

#include <TNamed.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>

/// \class Ali2PCorrelations
/// \brief Encapsulates the needed process and data for building
/// two-particle correlation functions
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date May, 2021

class TList;
class TH1;
class TH2;
class TH3;
class THn;
class AliVEvent;
class AliVTrack;
class AliVParticle;

class Ali2PCorrelations : public TNamed {
public:
  /// \enum TrackPairs
  /// \brief The track combinations hadled by the class
  enum TrackPairs {
    kOO = 0,     ///< one-one pairs
    kOT,         ///< one-two pairs
    kTO,         ///< two-one pairs
    kTT,         ///< two-two pairs
    nTrackPairs  ///< the number of track pairs
  };

  static const char *TrackPairsNames[nTrackPairs];

                              Ali2PCorrelations();
                              Ali2PCorrelations(const char *name);
  virtual                    ~Ali2PCorrelations();
                              /// \brief Use only singles, do not produce pairs data
                              /// \param so flag to activate deactivate
  void                        SetSinglesOnly(Bool_t so) { fSinglesOnly = so; }
                              /// \brief Use calibration weights
                              /// \param uw flag to activate deactivate the use of weights
  void                        SetUseWeights(Bool_t uw) { fUseWeights = uw; }
                              /// \brief Use simulated tracks
                              /// \param us flag to activate deactivate the use of simulated tracks
  void                        SetUseSimulation(Bool_t us) { fUseSimulation = us; }
                              /// \brief Get if simulation is being used
                              /// \return kTRUE if track simulation is being used kFALSE otherwise
  Bool_t                      GetUseSimulation() { return fUseSimulation; }
                              /// \brief Set if both tracks are of the same charge sign
                              /// \param ss kTRUE if both tracks have the same charge kFALSE otherwise

  Bool_t                      ConfigureBinning(const char *configstring);
  void                        ConfigureResonances(const char *confstring);
  TString                     GetBinningConfigurationString() const;
  TString                     GetResonancesConfigurationString() const;


  Bool_t                      SetWeigths(const TH3F *h3_1, const TH3F *h3_2);
  Bool_t                      SetEfficiencyCorrection(const TH1F *h1_1, const TH1F *h1_2);
  Bool_t                      SetPairEfficiencyCorrection(const THn *h11, const THn *h12, const THn *h22, const THn *h21);
  Bool_t                      SetSimultationPdfs(const TObjArray *pluspdf, const TObjArray *minuspdf);
                              /// \brief set the number of simulated events per passed event
                              /// \param nevents the number of events to simulate per passed event
  void                        SetSimEventsPerEvent(Int_t nevents) { fSimEventsPerEvent = nevents; }
                              /// \brief get the number of events being simulated per passed event
                              /// \return the number of events being simulated per passed event
  Int_t                       GetSimEventsPerEvent() { return fSimEventsPerEvent; }

  void                        Initialize();
                              /// Get the histograms list
                              /// \return the histograms list
  TList                      *GetHistogramsList() { return fOutput; }
  Bool_t                      StartEvent(Float_t centrality, Float_t vertexZ);
  Bool_t                      ProcessTrack(Int_t trkId, AliVTrack *trk);
  Bool_t                      ProcessTrack(Int_t trkId, AliVParticle *part);
  Bool_t                      ProcessTrack(Int_t trkId, Int_t charge, Float_t pT, Float_t eta, Float_t phi);
  void                        ProcessEventData();
  void                        FinalizeProcess();

private:
  void                        ProcessLikeSignPairs(Int_t bank);
  void                        FlagConversionsAndResonances();
  void                        ProcessUnlikeSignPairs();

private:
  TList                      *fOutput;                      //!<! Output histograms list

  Float_t                     fVertexZ;                     ///< the current event vertex \f$z\f$ coordinate
  Int_t                       fIxVertexZ;                   ///< bin index for the current event vertex \f$z\f$ coordinate
  Float_t                     fCentrality;                  ///< current event centrality in percentage

  Bool_t                      fSinglesOnly;                 ///< true if only collect singles data
  Bool_t                      fUseWeights;                  ///< kTRUE if correction weights must be utilized
  Bool_t                      fUseSimulation;               ///< kTRUE if particle production from stored profiles must be used

  /* the arrays with tracks 1 and 2 information */
  Int_t                       fArraySize;                   ///< the size of the array to collect accepted track information
  Int_t                       fNoOfTracks1;                 ///< the number of stored track 1 tracks
  Int_t                      *fId_1;                        //!<! the array of track 1 Ids
  Int_t                      *fCharge_1;                    //!<! the array of track 1 charge
  Int_t                      *fIxEta_1;                     //!<! the array of track 1 eta bin index
  Int_t                      *fIxPhi_1;                     //!<! the array of track 1 phi bin index
  Int_t                      *fIxPt_1;                      //!<! the array of track 1 pT bin index
  UInt_t                     *fFlags_1;                     //!<! the array of track 1 flags
  Float_t                    *fPt_1;                        //!<! the array of track 1 \f$p_T\f$
  Float_t                    *fEta_1;                       //!<! the array of track 1 \f$\eta\f$
  Float_t                    *fPhi_1;                       //!<! the array of track 1 \f$\varphi\f$
  Float_t                    *fCorrection_1;                //!<! the array of the correction to apply to track 1
  Int_t                       fNoOfTracks2;                 ///< the number of stored track 2 tracks
  Int_t                      *fId_2;                        //!<! the array of track 2 Ids
  Int_t                      *fCharge_2;                    //!<! the array of track 2 charge
  Int_t                      *fIxEta_2;                     //!<! the array of track 2 eta bin index
  Int_t                      *fIxPhi_2;                     //!<! the array of track 2 phi bin index
  Int_t                      *fIxPt_2;                      //!<! the array of track 2 pT bin index
  UInt_t                     *fFlags_2;                     //!<! the array of track 2 flags
  Float_t                    *fPt_2;                        //!<! the array of track 2 \f$p_T\f$
  Float_t                    *fEta_2;                       //!<! the array of track 2 \f$\eta\f$
  Float_t                    *fPhi_2;                       //!<! the array of track 2 \f$\varphi\f$
  Float_t                    *fCorrection_2;                //!<! the array of the correction to apply to track 2

  Float_t                    *fCorrectionWeights_1;         //!<! structure with the track 1 correction weights
  Float_t                    *fCorrectionWeights_2;         //!<! structure with the track 2 correction weights
  Float_t                    *fEfficiencyCorrection_1;      //!<! structure with the track 1 efficiency correction weights
  Float_t                    *fEfficiencyCorrection_2;      //!<! structure with the track 2 efficiency correction weights
  const THn                  *fPairsEfficiency_PP;          //!<! pair efficiency correction for the plus plus pair
  const THn                  *fPairsEfficiency_PM;          //!<! pair efficiency correction for the plus minus pair
  const THn                  *fPairsEfficiency_MM;          //!<! pair efficiency correction for the minus minus pair
  const THn                  *fPairsEfficiency_MP;          //!<! pair efficiency correction for the minus plus pair
  const TObjArray            *fPositiveTrackPdfs;           //!<! the positive tracks density distributions
  const TObjArray            *fNegativeTrackPdfs;           //!<! the negative tracks density distributions
  TH3F                       *fPositiveTrackCurrentPdf;     //!<! current event positive tracks density distribution
  TH3F                       *fNegativeTrackCurrentPdf;     //!<! current event negative tracks density distribution
  Int_t                       fSimEventsPerEvent;           ///< the number of simulated events produced per real event

  Int_t                       fNBins_vertexZ;               ///< the \f$z\f$ vertex component number of bins
  Double_t                    fMin_vertexZ;                 ///< the minimum \f$z\f$ vertex component value
  Double_t                    fMax_vertexZ;                 ///< the maximum \f$z\f$ vertex component value
  Double_t                    fWidth_vertexZ;               ///< the \f$z\f$ vertex component bin width

  Double_t                    fNBinsPhiShift;               ///< the number of bins the phi origin is shifted

  Int_t                       fNBins_pt_1;                  ///< the track 1 \f$p_T\f$ number of bins
  Double_t                    fMin_pt_1;                    ///< the track 1 minimum \f$p_T\f$ value
  Double_t                    fMax_pt_1;                    ///< the track 1 maximum \f$p_T\f$ value
  Double_t                    fWidth_pt_1;                  ///< the track 1 \f$p_T\f$ bin width
  Int_t                       fNBins_phi_1;                 ///< the track 1 \f$\phi\f$ number of bins
  Double_t                    fMin_phi_1;                   ///< the track 1 minimum \f$\phi\f$ value
  Double_t                    fMax_phi_1;                   ///< the track 1 maximum \f$\phi\f$ value
  Double_t                    fWidth_phi_1;                 ///< the track 1 \f$\phi\f$ bin width
  Int_t                       fNBins_eta_1;                 ///< the track 1 \f$\eta\f$ number of bins
  Double_t                    fMin_eta_1;                   ///< the track 1 minimum \f$\eta\f$ value
  Double_t                    fMax_eta_1;                   ///< the track 1 maximum \f$\eta\f$ value
  Double_t                    fWidth_eta_1;                 ///< the track 1 \f$\eta\f$ bin width
  Int_t                       fNBins_etaPhi_1;              ///< the track 1 combined \f$\eta, \phi\f$ number of bins
  Int_t                       fNBins_etaPhiPt_1;            ///< the track 1 combined \f$\eta, \phi, p_T\f$ number of bins
  Int_t                       fNBins_zEtaPhiPt_1;           ///< the combined event \f$z\f$ vertex component and track 1 \f$\eta, \phi, p_T\f$ number of bins

  Int_t                       fNBins_pt_2;                  ///< the track 2 \f$p_T\f$ number of bins
  Double_t                    fMin_pt_2;                    ///< the track 2 minimum \f$p_T\f$ value
  Double_t                    fMax_pt_2;                    ///< the track 2 maximum \f$p_T\f$ value
  Double_t                    fWidth_pt_2;                  ///< the track 2 \f$p_T\f$ bin width
  Int_t                       fNBins_phi_2;                 ///< the track 2 \f$\phi\f$ number of bins
  Double_t                    fMin_phi_2;                   ///< the track 2 minimum \f$\phi\f$ value
  Double_t                    fMax_phi_2;                   ///< the track 2 maximum \f$\phi\f$ value
  Double_t                    fWidth_phi_2;                 ///< the track 2 \f$\phi\f$ bin width
  Int_t                       fNBins_eta_2;                 ///< the track 2 \f$\eta\f$ number of bins
  Double_t                    fMin_eta_2;                   ///< the track 2 minimum \f$\eta\f$ value
  Double_t                    fMax_eta_2;                   ///< the track 2 maximum \f$\eta\f$ value
  Double_t                    fWidth_eta_2;                 ///< the track 2 \f$\eta\f$ bin width
  Int_t                       fNBins_etaPhi_2;              ///< the track 2 combined \f$\eta, \phi\f$ number of bins
  Int_t                       fNBins_etaPhiPt_2;            ///< the track 2 combined \f$\eta, \phi, p_T\f$ number of bins
  Int_t                       fNBins_zEtaPhiPt_2;           ///< the combined event \f$z\f$ vertex component and track 2 \f$\eta, \phi, p_T\f$ number of bins

  Int_t                       fNBins_deltaphi;              ///< the pair \f$\Delta\phi\f$ number of bins
  Double_t                    fMin_deltaphi;                ///< the pair minimum \f$\Delta\phi\f$ value
  Double_t                    fMax_deltaphi;                ///< the pair maximum \f$\Delta\phi\f$ value
  Double_t                    fWidth_deltaphi;              ///< the pair \f$\Delta\phi\f$ bin width
  Int_t                       fNBins_deltaeta;              ///< the pair \f$\Delta\eta\f$ number of bins
  Double_t                    fMin_deltaeta;                ///< the pair minimum \f$\Delta\eta\f$ value
  Double_t                    fMax_deltaeta;                ///< the pair maximum \f$\Delta\eta\f$ value
  Double_t                    fWidth_deltaeta;              ///< the pair \f$\Delta\eta\f$ bin width

  Double_t                    fN1_1;                        ///< weighted number of track 1 tracks for current event
  Double_t                    fN1_2;                        ///< weighted number of track 2 tracks for current event
  Double_t                    fNnw1_1;                      ///< not weighted number of track 1 tracks for current event
  Double_t                    fNnw1_2;                      ///< not weighted number of track 2 tracks for current event
  Double_t                    fSum1Pt_1;                    ///< accumulated sum of weighted track 1 \f$p_T\f$ for current event
  Double_t                    fSum1Pt_2;                    ///< accumulated sum of weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum1Ptnw_1;                  ///< accumulated sum of not weighted track 1 \f$p_T\f$ for current event
  Double_t                    fSum1Ptnw_2;                  ///< accumulated sum of not weighted track 2 \f$p_T\f$ for current event

  /* histograms */
  TH1F                       *fhN1_1_vsPt;                    //!<! track 1 weighted single particle distribution vs \f$p_T\f$
  TH2F                       *fhN1_1_vsEtaPhi;                //!<! track 1 weighted single particle distribution vs \f$\eta,\;\phi\f$
  TH2F                       *fhSum1Pt_1_vsEtaPhi;            //!<! track 1 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  TH3F                       *fhN1_1_vsZEtaPhiPt;             //!<! track 1 single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
  TH1F                       *fhN1_2_vsPt;                    //!<! track 2 weighted single particle distribution vs \f$p_T\f$
  TH2F                       *fhN1_2_vsEtaPhi;                //!<! track 2 weighted single particle distribution vs \f$\eta,\;\phi\f$
  TH2F                       *fhSum1Pt_2_vsEtaPhi;            //!<! track 2 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  TH3F                       *fhN1_2_vsZEtaPhiPt;             //!<! track 2 single particle distribution vs \f$\mbox{vtx}_z,\;\eta,\;\phi,\;p_T\f$
  TH2F                       *fhN2_12_vsPtPt[4];              //!<! track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  TH2F                       *fhN2_12_vsDEtaDPhi[4];          //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
  TH2F                       *fhSum2PtPt_12_vsDEtaDPhi[4];    //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
  TH2F                       *fhSum2DptDpt_12_vsDEtaDPhi[4];  //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations

  /* versus centrality  profiles */
  TProfile                   *fhN1_1_vsC;                   //!<! track 1 weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Pt_1_vsC;               //!<! track 1 accumulated sum of weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1nw_1_vsC;                 //!<! track 1 un-weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Ptnw_1_vsC;             //!<! track 1 accumulated sum of un-weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1_2_vsC;                   //!<! track 2 weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Pt_2_vsC;               //!<! track 2 accumulated sum of weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1nw_2_vsC;                 //!<! track 2 un-weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Ptnw_2_vsC;             //!<! track 2 accumulated sum of un-weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN2_12_vsC[4];               //!<! weighted accumulated two particle distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2PtPt_12_vsC[4];         //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2DptDpt_12_vsC[4];       //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhN2nw_12_vsC[4];             //!<! un-weighted accumulated two particle distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2PtPtnw_12_vsC[4];       //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2DptDptnw_12_vsC[4];     //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations

  static Int_t                fgkNoOfResonances;            ///< the number of resonances conversions to consider
  static Double_t             fgkMass[16];                  ///< the masses of resonances conversions to consider
  static Double_t             fgkChildMass[2][16];          ///< the masses of the resonances / conversions products
  static Double_t             fgkMassThreshold[16];         ///< the resonance / conversion mass threshold modulus
  Int_t                       fThresholdMult[16];           ///< the threshold multiplier, in 1/4 modulus units (i.e, four is one modulus, zero disable it)
  TH2F                       *fhResonanceRoughMasses;       ///< the resonance approximate invariant mass histogram
  TH2F                       *fhResonanceMasses;            ///< the resonance invariant mass histogram
  TH2F                       *fhDiscardedResonanceMasses;   ///< the discarded resonance invariant mass histogram


private:
  /// Copy constructor
  /// Not allowed. Forced private.
  Ali2PCorrelations(const Ali2PCorrelations&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  Ali2PCorrelations& operator=(const Ali2PCorrelations&);

  /// \cond CLASSIMP
  ClassDef(Ali2PCorrelations,1);
  /// \endcond
};

extern inline Float_t GetSquaredInvMass(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
extern inline Float_t GetSquaredInvMassCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
extern inline Float_t checkIfResonance(Int_t ires, Bool_t fill, Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2);

#endif // ALI2PCORRELATIONS_H
