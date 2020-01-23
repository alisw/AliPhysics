#ifndef ALICSPAIRANALYSIS_H
#define ALICSPAIRANALYSIS_H

/// \file AliCSPairAnalysis.h
/// \brief Class for analyzing particle pairs
///

#include <TNamed.h>
#include <THn.h>

/// \class AliCSPairAnalysis
/// \brief Encapsulates the needed process and data for analyzing
/// particle pairs
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 30, 2017

class TH3F;
class THn;
class TList;
class TObjArray;
class AliVParticle;
class AliAODMCParticle;
class AliVEvent;

#include "AliVTrack.h"

class AliCSPairAnalysis : public TNamed {
public:
                              AliCSPairAnalysis();
                              AliCSPairAnalysis(const char *name);
  virtual                    ~AliCSPairAnalysis();

  Bool_t                      ConfigureBinning(const char *configstring);
  TString                     GetBinningConfigurationString() const;

  void                        Initialize();
                              /// Stores the singles efficiency histograms reference for further use
                              /// \param h3_1 the efficiency of the first track
                              /// \param h3_2 the efficiency of the second track
                              /// \return kTRUE always
  Bool_t                      SetSinglesEfficiency(const TH3F *h3_1, const TH3F *h3_2)
                              { fPlusEfficiency = h3_1; fMinusEfficiency = h3_2; return kTRUE; }
                              /// Stores the pairs efficiency histograms reference for further use
                              /// \param h11 the efficiency for the first first pair
                              /// \param h12 the efficiency for the first second pair
                              /// \param h22 the efficiency for the second second pair
                              /// \param h21 the efficiency for the second first pair
                              /// \return kTRUE always
Bool_t                        SetPairEfficiency(const THn *h11, const THn *h12, const THn *h22, const THn *h21)
                              { fPairEfficiency_PP = h11; fPairEfficiency_PM = h12; fPairEfficiency_MM = h22; fPairEfficiency_MP = h21; return kTRUE; }
                              /// Get the histograms list
                              /// \return the histograms list
  TList                      *GetHistogramsList() { return fOutput; }
  Bool_t                      StartEvent(Float_t vertexZ);
  Bool_t                      ProcessTrack(Int_t, AliVTrack *trk);
  Bool_t                      ProcessTrack(Int_t, AliVParticle *par);
  void                        ProcessEventData();
  void                        FinalizeProcess();

private:
  static const Int_t          kgHistosDimension;            ///< the number of dimensions of the multidimensional histograms
  TList                      *fOutput;                      //!<! Output histograms list

  Int_t                       fNBins_vertexZ;               ///< the \f$z\f$ vertex component number of bins
  Double_t                    fMin_vertexZ;                 ///< the minimum \f$z\f$ vertex component value
  Double_t                    fMax_vertexZ;                 ///< the maximum \f$z\f$ vertex component value
  Double_t                    fWidth_vertexZ;               ///< the \f$z\f$ vertex component bin width

  Double_t                    fNBinsPhiShift;               ///< the number of bins the phi origin is shifted

  Int_t                       fNBins_pt;                    ///< the \f$p_T\f$ number of bins
  Double_t                    fMin_pt;                      ///< the minimum \f$p_T\f$ value
  Double_t                    fMax_pt;                      ///< the maximum \f$p_T\f$ value
  Double_t                    fWidth_pt;                    ///< the \f$p_T\f$ bin width
  Int_t                       fNBins_phi;                   ///< the \f$\phi\f$ number of bins
  Double_t                    fMin_phi;                     ///< the minimum \f$\phi\f$ value
  Double_t                    fMax_phi;                     ///< the maximum \f$\phi\f$ value
  Double_t                    fWidth_phi;                   ///< the \f$\phi\f$ bin width
  Int_t                       fNBins_eta;                   ///< the \f$\eta\f$ number of bins
  Double_t                    fMin_eta;                     ///< the minimum \f$\eta\f$ value
  Double_t                    fMax_eta;                     ///< the maximum \f$\eta\f$ value
  Double_t                    fWidth_eta;                   ///< the \f$\eta\f$ bin width

  /* the arrays with plus and minus tracks */
  TObjArray                  *fPlusTracksArray;             ///< the array of positive tracks
  TObjArray                  *fMinusTracksArray;            ///< the array of negative tracks

  /* the singles efficiency structures */
  const TH3F                 *fPlusEfficiency;              ///< the positive track efficiency
  const TH3F                 *fMinusEfficiency;             ///< the negative track efficiency
  Float_t                  ***fEventPlusEfficiency;         //!<! the processed event efficiency for the positive track, according to event zvtx
  Float_t                  ***fEventMinusEfficiency;        //!<! the processed event efficiency for the negative track, according to event zvtx

  /* the pairs efficiencies */
  const THn                  *fPairEfficiency_PP;           ///< the positive positive tracks pair efficiency
  const THn                  *fPairEfficiency_PM;           ///< the positive negative tracks pair efficiency
  const THn                  *fPairEfficiency_MM;           ///< the negative negative tracks pair efficiency
  const THn                  *fPairEfficiency_MP;           ///< the negative positive tracks pair efficiency

  /* histograms */
  TH1F                       *fhNplus;                      //!<! positive particle multiplicity
  TH1F                       *fhNminus;                     //!<! negative particle multiplicity
  THnF                       *fhPPDeltaEtaDeltaPhi;         //!<! positive-positive delta eta delta phi separation vs pT, pT
  THnF                       *fhPMDeltaEtaDeltaPhi;         //!<! positive-negative delta eta delta phi separation vs pT, pT
  THnF                       *fhMPDeltaEtaDeltaPhi;         //!<! positive-negative delta eta delta phi separation vs pT, pT
  THnF                       *fhMMDeltaEtaDeltaPhi;         //!<! negative-negative delta eta delta phi separation vs pT, pT

private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSPairAnalysis(const AliCSPairAnalysis&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSPairAnalysis& operator=(const AliCSPairAnalysis&);

  /// \cond CLASSIMP
  ClassDef(AliCSPairAnalysis,1);
  /// \endcond
};

#endif // ALICSPAIRANALYSIS_H
