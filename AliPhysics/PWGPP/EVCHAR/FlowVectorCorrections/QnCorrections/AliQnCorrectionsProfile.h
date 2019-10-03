#ifndef ALIQNCORRECTIONS_PROFILE_H
#define ALIQNCORRECTIONS_PROFILE_H

/// \file AliQnCorrectionsProfile.h
/// \brief Single multidimensional profile class for the Q vector correction framework

#include "AliQnCorrectionsHistogramBase.h"

/// \class AliQnCorrectionsProfile
/// \brief Single profile class for the Q vector correction histograms
///
/// Encapsulates a multi dimensional profile. Each dimension
/// corresponds to one of the event classes variables so,
/// the number of dimensions matches the number of variables within
/// the set passed in the constructor.
///
/// The involved histograms can be created on the fly when needed,
/// and included in a provided list, or attached to existing ones
/// from a given list. In any case they are not destroyed because
/// the are not own by the class but by the involved list.
///
/// GetBinContent (once the intended bin is obtained by mean
/// of GetBin) returns in the profile way
/// \f[
///    \frac{\Sigma \mbox{fValues(bin)}}{\mbox{fEntries(bin)}}
/// \f]
/// while GetBinError returns the standard deviation of the values
/// in the interested bin
/// \f[
///    \sqrt{\frac{\Sigma \mbox{fValues}^2\mbox{(bin)}}{\mbox{fEntries(bin)}}
///          - \left(\frac{\Sigma \mbox{fValues(bin)}}{\mbox{fEntries(bin)}}\right)^2}
/// \f]
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 11, 2016
class AliQnCorrectionsProfile : public AliQnCorrectionsHistogramBase {
public:
  AliQnCorrectionsProfile();
  AliQnCorrectionsProfile(const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsProfile();

  Bool_t CreateProfileHistograms(TList *histogramList);
  virtual Bool_t AttachHistograms(TList *histogramList);
  /// wrong call for this class invoke base class behavior
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup)
  { return AliQnCorrectionsHistogramBase::AttachHistograms(histogramList, bUsedChannel, nChannelGroup); }

  virtual Long64_t GetBin(const Float_t *variableContainer);
  /// wrong call for this class invoke base class behavior
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsHistogramBase::GetBin(variableContainer, nChannel); }
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetBinContent(Long64_t bin);
  virtual Float_t GetBinError(Long64_t bin);

  virtual void Fill(const Float_t *variableContainer, Float_t weight);
  /// wrong call for this class invoke base class behavior
  virtual void Fill(const Float_t *variableContainer, Int_t nChannel, Float_t weight)
  { AliQnCorrectionsHistogramBase::Fill(variableContainer, nChannel, weight); }
private:
  THnF *fValues;   //!<! Cumulates values for each of the event classes
  THnI *fEntries;  //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsProfile, 1);
  /// \endcond
};

#endif
