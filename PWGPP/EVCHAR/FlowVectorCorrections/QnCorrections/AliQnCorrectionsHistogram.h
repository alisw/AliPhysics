#ifndef ALIQNCORRECTIONS_HISTOGRAM_H
#define ALIQNCORRECTIONS_HISTOGRAM_H

/// \file AliQnCorrectionsHistogram.h
/// \brief Single multidimensional histograms

#include "AliQnCorrectionsHistogramBase.h"

/// \class AliQnCorrectionsHistogramChannelized
/// \brief Single histogram class for the Q vector correction histograms
///
/// Encapsulates a multi dimensional histogram. Each dimension
/// corresponds to one of the event classes variables so,
/// the number of dimensions matches the number of variables within
/// the set passed in the constructor.
///
/// The involved histograms can be created on the fly when needed,
/// and included in a provided list. They are not destroyed because
/// the are not own by the class but by the involved list.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jun 16, 2016
class AliQnCorrectionsHistogram : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsHistogram();
  AliQnCorrectionsHistogram(const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs);
  virtual ~AliQnCorrectionsHistogram();

  Bool_t CreateHistogram(TList *histogramList);

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behaviour
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer,nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetBinContent(Long64_t bin);
  virtual Float_t GetBinError(Long64_t bin);

  virtual void Fill(const Float_t *variableContainer, Float_t weight);
  /// wrong call for this class invoke base class behavior
  virtual void Fill(const Float_t *variableContainer, Int_t nChannel, Float_t weight)
  { AliQnCorrectionsHistogramBase::Fill(variableContainer, nChannel, weight); }
private:
  THnF *fValues;              //!<! Cumulates values for each of the event classes

  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsHistogram, 1);
  /// \endcond
};

#endif
