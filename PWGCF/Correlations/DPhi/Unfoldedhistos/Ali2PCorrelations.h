#ifndef ALI2PCORRELATIONS_H
#define ALI2PCORRELATIONS_H

/// \file Ali2PCorrelations.h
/// \brief Class for collecting differential data for two-particle correlation functions
///

#include "AliTwoParticleCorrelationsBase.h"

/// \class Ali2PCorrelations
/// \brief Encapsulates the needed process and data for building
/// two-particle correlation functions
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date May, 2021

class AliVEvent;
class AliVTrack;
class AliVParticle;

class Ali2PCorrelations : public AliTwoParticleCorrelationsBase {
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

  Bool_t                      ConfigureBinning(const char *configstring);
  TString                     GetBinningConfigurationString() const;

  Bool_t                      SetPtAvg(const TH2 *h2_1, const TH2 *h2_2);

  void                        Initialize();
                              /// Get the histograms list
                              /// \return the histograms list
  Bool_t                      StartEvent(Float_t centrality, Float_t vertexZ);
  Bool_t                      ProcessTrack(Int_t trkId, Int_t charge, Float_t pT, Float_t eta, Float_t phi);
  void                        ProcessEventData();
  void                        FinalizeProcess();

private:
  void                        ProcessLikeSignPairs(Int_t bank);
  void                        ProcessUnlikeSignPairs();

private:
  /* the average pT histograms */
  const TH2                  *fhPositivePtAverage;          //!<! corrected avg pT for positive tracks in eta phi
  const TH2                  *fhNegativePtAverage;          //!<! corrected avg pT for negative tracks in eta phi

  /* the arrays with tracks 1 and 2 information */
  Int_t                      *fId_1;                        //!<! the array of track 1 Ids
  Int_t                      *fCharge_1;                    //!<! the array of track 1 charge
  Int_t                      *fIxEta_1;                     //!<! the array of track 1 eta bin index
  Int_t                      *fIxPhi_1;                     //!<! the array of track 1 phi bin index
  Int_t                      *fIxPt_1;                      //!<! the array of track 1 pT bin index
  float                      *fAvgPt_1;                     //!<! the avg pT associated to the bin of track 1
  Float_t                    *fCorrection_1;                //!<! the array of the correction to apply to track 1
  Int_t                      *fId_2;                        //!<! the array of track 2 Ids
  Int_t                      *fCharge_2;                    //!<! the array of track 2 charge
  Int_t                      *fIxEta_2;                     //!<! the array of track 2 eta bin index
  Int_t                      *fIxPhi_2;                     //!<! the array of track 2 phi bin index
  Int_t                      *fIxPt_2;                      //!<! the array of track 2 pT bin index
  float                      *fAvgPt_2;                     //!<! the avg pT associated to the bin of track 2
  Float_t                    *fCorrection_2;                //!<! the array of the correction to apply to track 2

  Int_t                       fNBins_deltaphi;              ///< the pair \f$\Delta\phi\f$ number of bins
  Double_t                    fMin_deltaphi;                ///< the pair minimum \f$\Delta\phi\f$ value
  Double_t                    fMax_deltaphi;                ///< the pair maximum \f$\Delta\phi\f$ value
  Double_t                    fWidth_deltaphi;              ///< the pair \f$\Delta\phi\f$ bin width
  Int_t                       fNBins_deltaeta;              ///< the pair \f$\Delta\eta\f$ number of bins
  Double_t                    fMin_deltaeta;                ///< the pair minimum \f$\Delta\eta\f$ value
  Double_t                    fMax_deltaeta;                ///< the pair maximum \f$\Delta\eta\f$ value
  Double_t                    fWidth_deltaeta;              ///< the pair \f$\Delta\eta\f$ bin width

  /* histograms */
  TH2F                       *fhN2_12_vsPtPt[4];              //!<! track 1 and 2 weighted two particle distribution vs \f${p_T}_1, {p_T}_2\f$
  TH2F                       *fhN2_12_vsDEtaDPhi[4];          //!<! two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
  TH2F                       *fhSum2PtPt_12_vsDEtaDPhi[4];    //!<! two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations
  TH2F                       *fhSum2DptDpt_12_vsDEtaDPhi[4];  //!<! two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ 1-1,1-2,2-1,2-2, combinations

  /* versus centrality  profiles */
  TProfile                   *fhN2_12_vsC[4];               //!<! weighted accumulated two particle distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2PtPt_12_vsC[4];         //!<! weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2DptDpt_12_vsC[4];       //!<! weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhN2nw_12_vsC[4];             //!<! un-weighted accumulated two particle distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2PtPtnw_12_vsC[4];       //!<! un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations
  TProfile                   *fhSum2DptDptnw_12_vsC[4];     //!<! un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality 1-1,1-2,2-1,2-2, combinations

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

#endif // ALI2PCORRELATIONS_H
