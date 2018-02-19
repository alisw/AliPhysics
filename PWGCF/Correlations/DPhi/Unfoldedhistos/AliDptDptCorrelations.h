#ifndef ALIDPTDPTCORRELATIONS_H
#define ALIDPTDPTCORRELATIONS_H

/// \file AliDptDptCorrelations.h
/// \brief Class for collecting data for \f$(\Delta p_T, \Delta p_T)\f$ correlation functions
///

#include <TNamed.h>

/// \class AliDptDptCorrelations
/// \brief Encapsulates the needed process and data for building
/// \f$(\Delta p_T, \Delta p_T)\f$ correlation functions
///
/// Based on the work of P.Pujahari &amp C. Pruneau at Wayne State University
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Dec 04, 2016

class TList;
class TH1;
class TH2;
class TH3;
class THn;
class AliVEvent;
class AliVTrack;
class TParticle;
class AliAODMCParticle;

class AliDptDptCorrelations : public TNamed {
public:
                              AliDptDptCorrelations();
                              AliDptDptCorrelations(const char *name);
  virtual                    ~AliDptDptCorrelations();
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
  void                        SetSameSign(Bool_t ss) { fSameSign = ss; }
                              /// \brief Produce data for all track combinations
                              /// \param all kTRUE if all track combinations have to be produced kFALSE otherwise
  void                        SetAllCombinations(Bool_t all) { fAllCombinations = all; }
                              /// \brief The requested charge for first track
                              /// \param q charge, greater than one if positive negative otherwise
  void                        SetRequestedCharge_1(Int_t q) { fRequestedCharge_1 = q; }
                              /// \brief The requested charge for second track
                              /// \param q charge, greater than one if positive negative otherwise
  void                        SetRequestedCharge_2(Int_t q) { fRequestedCharge_2 = q; }

  void                        ConfigureBinning(const char *configstring);
  TString                     GetBinningConfigurationString() const;


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
  Bool_t                      ProcessTrack(Int_t trkId, TParticle *part);
  Bool_t                      ProcessTrack(Int_t trkId, AliAODMCParticle *part);
  Bool_t                      ProcessTrack(Int_t trkId, Int_t charge, Float_t pT, Float_t eta, Float_t phi);
  void                        ProcessEventData();
  void                        FinalizeProcess();

private:
  void                        ProcessLikeSignPairs(Int_t bank);
  void                        ProcessNotHalfSymmLikeSignPairs(Int_t bank);
  void                        ProcessUnlikeSignPairs();
  void                        ProcessNotHalfSymmUnlikeSignPairs();

  void                        fillHistoWithArray(TH1 * h, double * array, int size);
  void                        fillHistoWithArray(TH2 * h, double * array, int size1, int size2);
  void                        fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3);
  void                        fillHistoWithArray(TH1 * h, float * array, int size);
  void                        fillHistoWithArray(TH2 * h, float * array, int size1, int size2);
  void                        fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3);

private:
  TList                      *fOutput;                      //!<! Output histograms list

  Float_t                     fVertexZ;                     ///< the current event vertex \f$z\f$ coordinate
  Int_t                       fIxVertexZ;                   ///< bin index for the current event vertex \f$z\f$ coordinate
  Float_t                     fCentrality;                  ///< current event centrality in percentage

  Bool_t                      fHalfSymmetrize;              ///< kTRUE if half symmetrizing for memory layout reduction
  Bool_t                      fSinglesOnly;                 ///< kTRUE if not pair calculations, just single particle calculations
  Bool_t                      fUseWeights;                  ///< kTRUE if correction weights must be utilized
  Bool_t                      fUseSimulation;               ///< kTRUE if particle production from stored profiles must be used
  Bool_t                      fSameSign;                    ///< kTRUE if the same charge particles are utilized to build correlations
  Bool_t                      fAllCombinations;             ///< kTRUE if all track pair combinations are utilized to build correlations
  Int_t                       fRequestedCharge_1;           ///< requested charge sign for the first particle
  Int_t                       fRequestedCharge_2;           ///< requested charge sign for the second particle

  /* the arrays with tracks 1 and 2 information */
  Int_t                       fArraySize;                   ///< the size of the array to collect accepted track information
  Int_t                       fNoOfTracks1;                 ///< the number of stored track 1 tracks
  Int_t                      *fId_1;                        //!<! the array of track 1 Ids
  Int_t                      *fCharge_1;                    //!<! the array of track 1 charge
  Int_t                      *fIxEtaPhi_1;                  //!<! the array of track 1 combined eta phi bin index
  Int_t                      *fIxPt_1;                      //!<! the array of track 1 pT bin index
  Float_t                    *fPt_1;                        //!<! the array of track 1 pT
  Float_t                    *fCorrection_1;                //!<! the array of the correction to apply to track 1
  Int_t                       fNoOfTracks2;                 ///< the number of stored track 2 tracks
  Int_t                      *fId_2;                        //!<! the array of track 2 Ids
  Int_t                      *fCharge_2;                    //!<! the array of track 2 charge
  Int_t                      *fIxEtaPhi_2;                  //!<! the array of track 2 combined eta phi bin index
  Int_t                      *fIxPt_2;                      //!<! the array of track 2 pT bin index
  Float_t                    *fPt_2;                        //!<! the array of track 2 pT
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

  Int_t                       fNBins_etaPhi_12;             ///< the track 1 and 2 combined \f$\eta, \phi\f$ number of bins

  Double_t                    fN1_1;                        ///< weighted number of track 1 tracks for current event
  Double_t                    fN1_2;                        ///< weighted number of track 2 tracks for current event
  Double_t                    fN2_12;                       ///< weighted number of track1 track 2 pairs for current event
  Double_t                    fNnw1_1;                      ///< not weighted number of track 1 tracks for current event
  Double_t                    fNnw1_2;                      ///< not weighted number of track 2 tracks for current event
  Double_t                    fNnw2_12;                     ///< not weighted number of track1 track 2 pairs for current event
  Double_t                    fSum1Pt_1;                    ///< accumulated sum of weighted track 1 \f$p_T\f$ for current event
  Double_t                    fSum1Pt_2;                    ///< accumulated sum of weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum2PtPt_12;                 ///< accumulated sum of weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current event
  Double_t                    fSum1Ptnw_1;                  ///< accumulated sum of not weighted track 1 \f$p_T\f$ for current event
  Double_t                    fSum1Ptnw_2;                  ///< accumulated sum of not weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum2PtPtnw_12;               ///< accumulated sum of not weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current event
  Double_t                    fSum2NPt_12;                  ///< accumulated sum of weighted number of track 1 tracks times weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum2PtN_12;                  ///< accumulated sum of weighted track 1 \f$p_T\f$ times weighted number of track 2 tracks for current event
  Double_t                    fSum2NPtnw_12;                ///< accumulated sum of not weighted number of track 1 tracks times not weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum2PtNnw_12;                ///< accumulated sum of not weighted track 1 \f$p_T\f$ times not weighted number of track  tracks for current event

  Double_t                   *fN1_1_vsPt;                   //!<! storage for track 1 weighted single particle distribution vs \f$p_T\f$
  Double_t                   *fN1_1_vsEtaPhi;               //!<! storage for track 1 weighted single particle distribution vs \f$\eta,\;\phi\f$
  Double_t                   *fSum1Pt_1_vsEtaPhi;           //!<! storage for track 1 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  Float_t                    *fN1_1_vsZEtaPhiPt;            //!<! storage for track 1 single particle distribution vs \f$\mbox{vtx}_z,\;\eta,\;\phi,\;p_T\f$
  Double_t                   *fN1_2_vsPt;                   //!<! storage for track 2 weighted single particle distribution vs \f$p_T\f$
  Double_t                   *fN1_2_vsEtaPhi;               //!<! storage for track 2 weighted single particle distribution vs \f$\eta,\;\phi\f$
  Double_t                   *fSum1Pt_2_vsEtaPhi;           //!<! storage for track 2 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  Float_t                    *fN1_2_vsZEtaPhiPt;            //!<! storage for track 2 single particle distribution vs \f$\mbox{vtx}_z,\;\eta,\;\phi,\;p_T\f$
  Double_t                   *fN2_12_vsPtPt;                //!<! storage for track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  Float_t                    *fN2_12_vsEtaPhi;              //!<! storage for track 1 and 2 weighted two particle distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2PtPt_12_vsEtaPhi;        //!<! storage for track 1 and 2 weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2PtN_12_vsEtaPhi;         //!<! storage for track 1 and 2 weighted accumulated \f${p_T}_1 n_2\f$ distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2NPt_12_vsEtaPhi;         //!<! storage for track 1 and 2 weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$

  /* histograms */
  TH1F                       *fhN1_1_vsPt;                  //!<! track 1 weighted single particle distribution vs \f$p_T\f$
  TH2F                       *fhN1_1_vsEtaPhi;              //!<! track 1 weighted single particle distribution vs \f$\eta,\;\phi\f$
  TH2F                       *fhSum1Pt_1_vsEtaPhi;          //!<! track 1 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  TH3F                       *fhN1_1_vsZEtaPhiPt;           //!<! track 1 single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$
  TH1F                       *fhN1_2_vsPt;                  //!<! track 2 weighted single particle distribution vs \f$p_T\f$
  TH2F                       *fhN1_2_vsEtaPhi;              //!<! track 2 weighted single particle distribution vs \f$\eta,\;\phi\f$
  TH2F                       *fhSum1Pt_2_vsEtaPhi;          //!<! track 2 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  TH3F                       *fhN1_2_vsZEtaPhiPt;           //!<! track 2 single particle distribution vs \f$\mbox{vtx}_z,\;\eta,\;\phi,\;p_T\f$
  TH2F                       *fhN2_12_vsPtPt;               //!<! track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  TH1F                       *fhN2_12_vsEtaPhi;             //!<! track 1 and 2 weighted two particle distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2PtPt_12_vsEtaPhi;       //!<! track 1 and 2 weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2PtN_12_vsEtaPhi;        //!<! track 1 and 2 weighted accumulated \f${p_T}_1 n_2\f$ distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2NPt_12_vsEtaPhi;        //!<! track 1 and 2 weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  /* versus centrality  profiles */
  TProfile                   *fhN1_1_vsC;                   //!<! track 1 weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Pt_1_vsC;               //!<! track 1 accumulated sum of weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1nw_1_vsC;                 //!<! track 1 un-weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Ptnw_1_vsC;              //!<! track 1 accumulated sum of un-weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1_2_vsC;                   //!<! track 2 weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Pt_2_vsC;               //!<! track 2 accumulated sum of weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN1nw_2_vsC;                 //!<! track 2 un-weighted single particle distribution vs event centrality
  TProfile                   *fhSum1Ptnw_2_vsC;             //!<! track 2 accumulated sum of un-weighted \f$p_T\f$ vs event centrality
  TProfile                   *fhN2_12_vsC;                  //!<! track 1 and 2 weighted two particle distribution vs event centrality
  TProfile                   *fhSum2PtPt_12_vsC;            //!<! track 1 and 2 weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality
  TProfile                   *fhSum2PtN_12_vsC;             //!<! track 1 and 2 weighted accumulated \f${p_T}_1 n_2\f$ distribution vs event centrality
  TProfile                   *fhSum2NPt_12_vsC;             //!<! track 1 and 2 weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs event centrality
  TProfile                   *fhN2nw_12_vsC;                //!<! track 1 and 2 un-weighted two particle distribution vs event centrality
  TProfile                   *fhSum2PtPtnw_12_vsC;          //!<! track 1 and 2 un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality
  TProfile                   *fhSum2PtNnw_12_vsC;           //!<! track 1 and 2 un-weighted accumulated \f${p_T}_1 n_2\f$ distribution vs event centrality
  TProfile                   *fhSum2NPtnw_12_vsC;           //!<! track 1 and 2 un-weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs event centrality


private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliDptDptCorrelations(const AliDptDptCorrelations&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliDptDptCorrelations& operator=(const AliDptDptCorrelations&);

  /// \cond CLASSIMP
  ClassDef(AliDptDptCorrelations,2);
  /// \endcond
};

#endif // ALIDPTDPTCORRELATIONS_H
