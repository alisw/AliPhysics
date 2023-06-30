#ifndef ALIDPTDPTCORRELATIONS_H
#define ALIDPTDPTCORRELATIONS_H

/// \file AliDptDptCorrelations.h
/// \brief Class for collecting unfolded data for two particle correlation functions
///

#include "AliTwoParticleCorrelationsBase.h"

/// \class AliDptDptCorrelations
/// \brief Encapsulates the needed process and data for building
/// two-particle correlation functions
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
class AliVParticle;

class AliDptDptCorrelations : public AliTwoParticleCorrelationsBase {
public:
                              AliDptDptCorrelations();
                              AliDptDptCorrelations(const char *name);
  virtual                    ~AliDptDptCorrelations();
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

  virtual Bool_t              ConfigureBinning(const char *configstring);
  virtual TString             GetBinningConfigurationString() const;

  void                        Initialize();
                              /// Get the histograms list
                              /// \return the histograms list
  Bool_t                      StartEvent(Float_t centrality, Float_t vertexZ);
  bool ProcessTrack(int pid, float pT, float eta, float phi);
  void                        ProcessEventData();
  void FinalizeProcess();

 private:
  enum kLSUS {
    kLS,
    kUS,
    kALL
  };
  enum kSYMM {
    kHALF,
    kFULL
  };

  template <kLSUS kind, kSYMM half>
  void ProcessPairs();

  void fillHistoWithArray(TH1* h, double* array, int size);
  void fillHistoWithArray(TH2* h, double* array, int size1, int size2);
  void fillHistoWithArray(TH3* h, double* array, int size1, int size2, int size3);
  void fillHistoWithArray(TH1* h, float* array, int size);
  void fillHistoWithArray(TH2* h, float* array, int size1, int size2);
  void fillHistoWithArray(TH3* h, float* array, int size1, int size2, int size3);

 private:
  Bool_t                      fHalfSymmetrize;              ///< kTRUE if half symmetrizing for memory layout reduction
  Bool_t                      fSameSign;                    ///< kTRUE if the same charge particles are utilized to build correlations
  Bool_t                      fAllCombinations;             ///< kTRUE if all track pair combinations are utilized to build correlations
  Int_t                       fRequestedCharge_1;           ///< requested charge sign for the first particle
  Int_t                       fRequestedCharge_2;           ///< requested charge sign for the second particle

  /* the arrays with tracks 1 and 2 information */
  int* fIxEtaPhi;     //!<! the array of track combined eta phi bin index
  int* fIxPt;         //!<! the array of track pT bin index
  float* fCorrection; //!<! the array of the correction to apply to track

  Int_t fNBins_etaPhi_12; ///< the track 1 and 2 combined \f$\eta, \phi\f$ number of bins

  Double_t                    fN2_12;                       ///< weighted number of track1 track 2 pairs for current event
  Double_t                    fNnw2_12;                     ///< not weighted number of track1 track 2 pairs for current event
  Double_t                    fSum2PtPt_12;                 ///< accumulated sum of weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current event
  Double_t                    fSum2PtPtnw_12;               ///< accumulated sum of not weighted track 1 track 2 \f${p_T}_1 {p_T}_2\f$ for current event
  Double_t                    fSum2NPt_12;                  ///< accumulated sum of weighted number of track 1 tracks times weighted track 2 \f$p_T\f$ for current event
  Double_t                    fSum2PtN_12;                  ///< accumulated sum of weighted track 1 \f$p_T\f$ times weighted number of track 2 tracks for current event
  Double_t                    fSum2NPtnw_12;                ///< accumulated sum of not weighted number of track 1 tracks times not weighted track 2 \f$p_T\f$ for current event
  Double_t fSum2PtNnw_12;                                   ///< accumulated sum of not weighted track 1 \f$p_T\f$ times not weighted number of track  tracks for current event

  std::vector<double*> fN1_vsPt;                            //!<! storage for track 1 and 2 weighted single particle distribution vs \f$p_T\f$
  std::vector<double*> fN1_vsEtaPhi;                        //!<! storage for track 1 and 2 weighted single particle distribution vs \f$\eta,\;\phi\f$
  std::vector<double*> fSum1Pt_vsEtaPhi;                    //!<! storage for track 1 and 2 accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$
  std::vector<float*> fN1_vsZEtaPhiPt;                      //!<! storage for track 1 and 2 single particle distribution vs \f$\mbox{vtx}_z,\;\eta,\;\phi,\;p_T\f$
  Double_t                   *fN2_12_vsPtPt;                //!<! storage for track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  Float_t                    *fN2_12_vsEtaPhi;              //!<! storage for track 1 and 2 weighted two particle distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2PtPt_12_vsEtaPhi;        //!<! storage for track 1 and 2 weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2PtN_12_vsEtaPhi;         //!<! storage for track 1 and 2 weighted accumulated \f${p_T}_1 n_2\f$ distribution vs \f$\eta,\;\phi\f$
  Float_t                    *fSum2NPt_12_vsEtaPhi;         //!<! storage for track 1 and 2 weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$

  /* histograms */
  TH2F                       *fhN2_12_vsPtPt;               //!<! track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  TH1F                       *fhN2_12_vsEtaPhi;             //!<! track 1 and 2 weighted two particle distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2PtPt_12_vsEtaPhi;       //!<! track 1 and 2 weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2PtN_12_vsEtaPhi;        //!<! track 1 and 2 weighted accumulated \f${p_T}_1 n_2\f$ distribution vs \f$\eta,\;\phi\f$
  TH1F                       *fhSum2NPt_12_vsEtaPhi;        //!<! track 1 and 2 weighted accumulated \f$n_1 {p_T}_2\f$ distribution vs \f$\eta,\;\phi\f$
  /* versus centrality  profiles */
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
  ClassDef(AliDptDptCorrelations, 5);
  /// \endcond
};

#endif // ALIDPTDPTCORRELATIONS_H
