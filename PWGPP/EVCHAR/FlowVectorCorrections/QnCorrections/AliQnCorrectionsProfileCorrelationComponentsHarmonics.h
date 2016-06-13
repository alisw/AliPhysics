#ifndef ALIQNCORRECTIONS_PROFILECORRCOMPHARM_H
#define ALIQNCORRECTIONS_PROFILECORRCOMPHARM_H

/// \file AliQnCorrectionsProfileCorrelationComponentsHarmonics.h
/// \brief Correlation components based set of profiles with harmonic support for the Q vector correction framework

#include "AliQnCorrectionsHistogramBase.h"

/// \class AliQnCorrectionsProfileCorrelationComponentsHarmonics
/// \brief Base class for the correlation components based set of profiles
///
/// Provides profile histograms for storing component, XX, XY, YX, YY, based
/// information for a set of predefined harmonics defined at creation
/// time. The user can select the harmonic addressing procedure so that
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
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 19, 2016
class AliQnCorrectionsProfileCorrelationComponentsHarmonics : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsProfileCorrelationComponentsHarmonics();
  AliQnCorrectionsProfileCorrelationComponentsHarmonics(
      const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsProfileCorrelationComponentsHarmonics();

  Bool_t CreateCorrelationComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t *harmonicMap = NULL);
  virtual Bool_t AttachHistograms(TList *histogramList);
  /// wrong call for this class invoke base class behavior
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup)
  { return AliQnCorrectionsHistogramBase::AttachHistograms(histogramList, bUsedChannel, nChannelGroup); }

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behavior
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer, nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetXXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinError(Int_t harmonic, Long64_t bin);

  virtual void FillXX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillXY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillYX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillYY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);

private:
  THnF **fXXValues;            //!<! XX component histogram for each requested harmonic
  THnF **fXYValues;            //!<! XY component histogram for each requested harmonic
  THnF **fYXValues;            //!<! YX component histogram for each requested harmonic
  THnF **fYYValues;            //!<! YY component histogram for each requested harmonic
  UInt_t fXXharmonicFillMask;  //!<! keeps track of harmonic XX component filled values
  UInt_t fXYharmonicFillMask;  //!<! keeps track of harmonic XY component filled values
  UInt_t fYXharmonicFillMask;  //!<! keeps track of harmonic YX component filled values
  UInt_t fYYharmonicFillMask;  //!<! keeps track of harmonic YY component filled values
  UInt_t fFullFilled;          //!<! mask for the fully filled condition
  THnI  *fEntries;             //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsProfileCorrelationComponentsHarmonics, 1);
  /// \endcond
};

#endif
