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
                              Ali2PCorrelations();
                              Ali2PCorrelations(const char *name);
  virtual                    ~Ali2PCorrelations();

  Bool_t                      ConfigureBinning(const char *configstring);
  TString GetBinningConfigurationString() const;

  bool SetPtAvg(std::vector<const TH2*>);

  void                        Initialize();
                              /// Get the histograms list
                              /// \return the histograms list
  Bool_t                      StartEvent(Float_t centrality, Float_t vertexZ);
  bool ProcessTrack(int pid, float pT, float eta, float phi);
  void                        ProcessEventData();
  void FinalizeProcess();

 private:
  template<bool chargedhadrons>
  void ProcessPairs();

 private:
  /* the average pT histograms */
  std::vector<const TH2*> fhPtAverage; //!<! corrected avg pT for each of the species

  /* the arrays with track information */
  int* fIxEta;        //!<! the array of track eta bin index
  int* fIxPhi;        //!<! the array of track phi bin index
  int* fIxPt;         //!<! the array of track pT bin index
  float* fAvgPt;      //!<! the array of the avg pT associated to the bin of track
  float* fCorrection; //!<! the array of the correction to apply to track

  Int_t fNBins_deltaphi;    ///< the pair \f$\Delta\phi\f$ number of bins
  Double_t fMin_deltaphi;   ///< the pair minimum \f$\Delta\phi\f$ value
  Double_t fMax_deltaphi;   ///< the pair maximum \f$\Delta\phi\f$ value
  Double_t fWidth_deltaphi; ///< the pair \f$\Delta\phi\f$ bin width
  Int_t fNBins_deltaeta;    ///< the pair \f$\Delta\eta\f$ number of bins
  Double_t fMin_deltaeta;   ///< the pair minimum \f$\Delta\eta\f$ value
  Double_t fMax_deltaeta;   ///< the pair maximum \f$\Delta\eta\f$ value
  Double_t fWidth_deltaeta; ///< the pair \f$\Delta\eta\f$ bin width

  /* histograms */
  std::vector<std::vector<TH2F*>> fhN2_12_vsPtPt;                ///< two-particle distribution vs \f${p_T}_1, {p_T}_2\f$
  std::vector<std::vector<TH2F*>> fhN2_12_vsDEtaDPhi;            ///< two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$
  std::vector<std::vector<TH2F*>> fhN2_12_vsDEtaDPhi_na;         ///< two-particle distribution vs \f$\Delta\eta,\;\Delta\phi\f$, no aliasing
  std::vector<std::vector<TH2F*>> fhSum2PtPt_12_vsDEtaDPhi;      ///< two-particle  \f$\sum {p_T}_1 {p_T}_2\f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
  std::vector<std::vector<TH2F*>> fhSum2DptDpt_12_vsDEtaDPhi;    ///< two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$
  std::vector<std::vector<TH2F*>> fhSum2DptDpt_12_vsDEtaDPhi_na; ///< two-particle  \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$, no aliasing

  /* versus centrality/multiplicity  profiles */
  std::vector<std::vector<TProfile*>> fhN2_12_vsC;           ///< weighted accumulated two particle distribution vs event centrality/multiplicity
  std::vector<std::vector<TProfile*>> fhSum2PtPt_12_vsC;     ///< weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
  std::vector<std::vector<TProfile*>> fhSum2DptDpt_12_vsC;   ///< weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs event centrality/multiplicity
  std::vector<std::vector<TProfile*>> fhN2nw_12_vsC;         ///< un-weighted accumulated two particle distribution vs event centrality/multiplicity
  std::vector<std::vector<TProfile*>> fhSum2PtPtnw_12_vsC;   ///< un-weighted accumulated \f${p_T}_1 {p_T}_2\f$ distribution vs event centrality/multiplicity
  std::vector<std::vector<TProfile *>>
      fhSum2DptDptnw_12_vsC; ///< un-weighted accumulated \f$\sum ({p_T}_1- <{p_T}_1>) ({p_T}_2 - <{p_T}_2>) \f$ distribution vs \f$\Delta\eta,\;\Delta\phi\f$ distribution vs event centrality/multiplicity

  /* have the acumulators allocated from the beginning */
  std::vector<std::vector<double>> fN2_12;
  std::vector<std::vector<double>> fSum2PtPt_12;
  std::vector<std::vector<double>> fSum2DptDpt_12;
  std::vector<std::vector<double>> fNnw2_12;
  std::vector<std::vector<double>> fSum2PtPtnw_12;
  std::vector<std::vector<double>> fSum2DptDptnw_12;

  private:
  /// Copy constructor
  /// Not allowed. Forced private.
  Ali2PCorrelations(const Ali2PCorrelations&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  Ali2PCorrelations& operator=(const Ali2PCorrelations&);

  /// \cond CLASSIMP
  ClassDef(Ali2PCorrelations, 6);
  /// \endcond
};

#endif // ALI2PCORRELATIONS_H
