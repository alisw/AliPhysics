#ifndef ALIQNCORRECTIONS_PROFILECORRCOMP_H
#define ALIQNCORRECTIONS_PROFILECORRCOMP_H

/// \file AliQnCorrectionsProfileCorrelationComponents.h
/// \brief Correlation components based set of profiles for the Q vector correction framework

#include "AliQnCorrectionsHistogramBase.h"

/// \class AliQnCorrectionsProfileCorrelationComponents
/// \brief Base class for the correlation components based set of profiles
///
/// Provides profile histograms for storing component, XX, XY, YX, YY, based
/// information.
///
/// As any histogram derived from AliQnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor.
/// Of course,  the base name and base title for the different
/// histograms has also to be provided.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date May 08, 2016
class AliQnCorrectionsProfileCorrelationComponents : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsProfileCorrelationComponents();
  AliQnCorrectionsProfileCorrelationComponents(
      const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsProfileCorrelationComponents();

  Bool_t CreateCorrelationComponentsProfileHistograms(TList *histogramList);
  virtual Bool_t AttachHistograms(TList *histogramList);
  /// wrong call for this class invoke base class behavior
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup)
  { return AliQnCorrectionsHistogramBase::AttachHistograms(histogramList, bUsedChannel, nChannelGroup); }

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behavior
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer, nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetXXBinContent(Long64_t bin);
  virtual Float_t GetXYBinContent(Long64_t bin);
  virtual Float_t GetYXBinContent(Long64_t bin);
  virtual Float_t GetYYBinContent(Long64_t bin);
  virtual Float_t GetXXBinError(Long64_t bin);
  virtual Float_t GetXYBinError(Long64_t bin);
  virtual Float_t GetYXBinError(Long64_t bin);
  virtual Float_t GetYYBinError(Long64_t bin);

  virtual void FillXX(const Float_t *variableContainer, Float_t weight);
  virtual void FillXY(const Float_t *variableContainer, Float_t weight);
  virtual void FillYX(const Float_t *variableContainer, Float_t weight);
  virtual void FillYY(const Float_t *variableContainer, Float_t weight);

private:
  THnF *fXXValues;            //!<! XX component histogram
  THnF *fXYValues;            //!<! XY component histogram
  THnF *fYXValues;            //!<! YX component histogram
  THnF *fYYValues;            //!<! YY component histogram
  UInt_t fXXXYYXYYFillMask;   //!<! keeps track of component filled values
  UInt_t fFullFilled;          //!<! mask for the fully filled condition
  THnI  *fEntries;             //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsProfileCorrelationComponents, 1);
  /// \endcond
};

#endif
