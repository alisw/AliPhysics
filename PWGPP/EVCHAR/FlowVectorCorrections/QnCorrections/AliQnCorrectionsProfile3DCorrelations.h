#ifndef ALIQNCORRECTIONS_PROFILE3DCORRELATIONS_H
#define ALIQNCORRECTIONS_PROFILE3DCORRELATIONS_H

/// \file AliQnCorrectionsProfile3DCorrelations.h
/// \brief Three detector correlation components based set of profiles with harmonic support for the Q vector correction framework

#include "AliQnCorrectionsHistogramBase.h"

class AliQnCorrectionsQnVector;

/// \class AliQnCorrectionsProfile3DCorrelations
/// \brief Base class for three detectors correlation components based set of profiles with harmonic support
///
/// Provides profile histograms for storing component, XX, XY, YX, YY, based
/// information for a set of predefined harmonics defined at creation
/// time and for a set of three different subdetectors generically identified
/// as A, B and C. The user can select the harmonic addressing procedure so that
/// it will possible to ask for just one harmonic support and assign
/// to it any desired number.
///
/// As any histogram derived from AliQnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor together with the required number of
/// harmonics and an optional harmonic expected numbering scheme.
/// Of course,  the base name and base title for the different
/// histograms has also to be provided.
///
/// The harmonic map passed should contain an ordered array with
/// as many items as requested harmonics that provides the external
/// number to be used for request the corresponding harmonic.
/// Requesting five harmonics without maps is equivalent to pass
/// {1,2,3,4,5} as map. Requesting just support for the harmonic
/// four will require a map {4}.
///
/// Externally the harmonic number is addressed as usual. An additional
/// harmonic multiplier field allows to handle mxn harmonics. n is always
/// the external harmonic required internally it si handled as well as
/// n but all the information manipulated is really associated to mxn.
/// Only in the histograms name it appears the proper mxn harmonic to
/// not confuse the external user which browse the histograms.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 19, 2016
class AliQnCorrectionsProfile3DCorrelations : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsProfile3DCorrelations();
  AliQnCorrectionsProfile3DCorrelations(
      const char *name,
      const char *title,
      const char *nameA,
      const char *nameB,
      const char *nameC,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsProfile3DCorrelations();

  Bool_t CreateCorrelationComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t nHarmonicMultiplier = 1, Int_t *harmonicMap = NULL);
  virtual Bool_t AttachHistograms(TList *histogramList);
  /// wrong call for this class invoke base class behavior
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup)
  { return AliQnCorrectionsHistogramBase::AttachHistograms(histogramList, bUsedChannel, nChannelGroup); }

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behavior
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer, nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetXXBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXXBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinError(const char *comb, Int_t harmonic, Long64_t bin);

  virtual void Fill(const AliQnCorrectionsQnVector *QnA,
      const AliQnCorrectionsQnVector *QnB,
      const AliQnCorrectionsQnVector *QnC,
      const Float_t *variableContainer);

private:
  THnF ***fXXValues;            //!<! XX component histogram for each requested harmonic
  THnF ***fXYValues;            //!<! XY component histogram for each requested harmonic
  THnF ***fYXValues;            //!<! YX component histogram for each requested harmonic
  THnF ***fYYValues;            //!<! YY component histogram for each requested harmonic
  THnI  *fEntries;             //!<! Cumulates the number on each of the event classes
  TString fNameA;               ///< the name of the A detector
  TString fNameB;               ///< the name of the B detector
  TString fNameC;               ///< the name of the C detector
  Int_t fHarmonicMultiplier;    ///< the multiplier for the harmonic number
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsProfile3DCorrelations, 1);
  /// \endcond
};

#endif
