#ifndef ALIQNCORRECTIONS_PROFILECOMP_H
#define ALIQNCORRECTIONS_PROFILECOMP_H

/// \file AliQnCorrectionsProfileComponents.h
/// \brief Component based set of profiles for the Q vector correction framework

#include "AliQnCorrectionsHistogramBase.h"

/// \class AliQnCorrectionsProfileComponents
/// \brief Base class for the components based set of profiles
///
/// Provides profile histograms for storing component, X, Y, based
/// information for a set of predefined harmonics defined at creation
/// time. The user can select the harmonic addressing procedure so that
/// it will possible to ask for just one harmonic support and assign
/// to it any desired number.
///
/// As any histogram derived from AliQnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor. Of course,  the base name and base
/// title for the different histograms has to be also provided.
/// At creation time the required number of harmonics and an optional
/// expected harmonic numbering scheme has to be passed.
///
/// The harmonic map passed should contain an ordered array with
/// as many items as requested harmonics that provides the external
/// number to be used for request the corresponding harmonic.
/// Requesting five harmonics without maps is equivalent to pass
/// {1,2,3,4,5} as map. Requesting just support for the harmonic
/// four will require a map {4}.
///
/// When you fill the histograms care is taken for not increasing
/// the number of entries until all components for the whole set of
/// harmonics have been filled. If you try to fill twice a harmonic
/// component before the whole set is filled you will get an execution
/// error because you are doing something that shall be corrected
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 15, 2016
class AliQnCorrectionsProfileComponents : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsProfileComponents();
  AliQnCorrectionsProfileComponents(const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsProfileComponents();

  Bool_t CreateComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t *harmonicMap = NULL);
  virtual Bool_t AttachHistograms(TList *histogramList);
  /// wrong call for this class invoke base class behavior
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup)
  { return AliQnCorrectionsHistogramBase::AttachHistograms(histogramList, bUsedChannel, nChannelGroup); }

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behavior
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer, nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYBinError(Int_t harmonic, Long64_t bin);

  virtual void FillX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);

private:
  THnF **fXValues;            //!<! X component histogram for each requested harmonic
  THnF **fYValues;            //!<! Y component histogram for each requested harmonic
  UInt_t fXharmonicFillMask;  //!<! keeps track of harmonic X component filled values
  UInt_t fYharmonicFillMask;  //!<! keeps track of harmonic Y component filled values
  UInt_t fFullFilled;         //!<! mask for the fully filled condition
  THnI  *fEntries;            //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsProfileComponents, 1);
  /// \endcond
};

#endif
