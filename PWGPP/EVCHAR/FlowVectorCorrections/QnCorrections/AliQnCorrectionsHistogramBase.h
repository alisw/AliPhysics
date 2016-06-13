#ifndef ALIQNCORRECTIONS_HISTOGRAMSBASE_H
#define ALIQNCORRECTIONS_HISTOGRAMSBASE_H

/// \file AliQnCorrectionsHistogramBase.h
/// \brief Multidimensional profile histograms base class for the Q vector correction framework

#include "THn.h"

class AliQnCorrectionsEventClassVariablesSet;

/// \class AliQnCorrectionsHistogramBase
/// \brief Base class for the Q vector correction histograms
///
/// Basically stores the set of variables that identify
/// the different event classes the involved histograms
/// are storing information about. It also stores (in its
/// parent) the base name and base title for the different
/// histograms it will encapsulate.
///
/// The passed at construction option parameter is the option for the
/// the computation of the  errors in the descendant profiles. Possible
/// values for the options are:
///
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
///
/// The encapsulated bin axes values provide an efficient
/// runtime storage for computing bin numbers.
///
/// Provides the interface for the whole set of histogram
/// classes providing error information that helps debugging.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 11, 2016
class AliQnCorrectionsHistogramBase : public TNamed {
protected:
  /// \typedef QnCorrectionHistogramErrorMode
  /// \brief The type of bin errors supported by the framework histograms.
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  typedef enum {
    kERRORMEAN = 0,                 ///< the bin errors are the standard error on the mean
    kERRORSPREAD                    ///< the bin errors are the standard deviation
  } QnCorrectionHistogramErrorMode;
public:
  AliQnCorrectionsHistogramBase();
  AliQnCorrectionsHistogramBase(const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs,
      Option_t *option="");
  virtual ~AliQnCorrectionsHistogramBase();

  virtual Bool_t AttachHistograms(TList *histogramList);
  virtual Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup);


  virtual Long64_t GetBin(const Float_t *variableContainer);
  virtual Long64_t GetBin(const Float_t *variableContainer, Int_t nChannel);
  /// Check the validity of the content of the passed bin
  /// Pure virtual function
  /// \param bin the bin to check its content validity
  /// \return kTRUE if the content is valid kFALSE otherwise
  virtual Bool_t BinContentValidated(Long64_t bin) = 0;

  virtual Float_t GetBinContent(Long64_t bin);
  virtual Float_t GetXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXXBinContent(Long64_t bin);
  virtual Float_t GetXYBinContent(Long64_t bin);
  virtual Float_t GetYXBinContent(Long64_t bin);
  virtual Float_t GetYYBinContent(Long64_t bin);
  virtual Float_t GetXXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinContent(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinContent(Int_t harmonic, Long64_t bin);

  virtual Float_t GetBinError(Long64_t bin);
  virtual Float_t GetXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXXBinError(Long64_t bin);
  virtual Float_t GetXYBinError(Long64_t bin);
  virtual Float_t GetYXBinError(Long64_t bin);
  virtual Float_t GetYYBinError(Long64_t bin);
  virtual Float_t GetXXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinError(Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinError(Int_t harmonic, Long64_t bin);

  virtual void Fill(const Float_t *variableContainer, Float_t weight);
  virtual void Fill(const Float_t *variableContainer, Int_t nChannel, Float_t weight);
  virtual void FillX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillXX(const Float_t *variableContainer, Float_t weight);
  virtual void FillXY(const Float_t *variableContainer, Float_t weight);
  virtual void FillYX(const Float_t *variableContainer, Float_t weight);
  virtual void FillYY(const Float_t *variableContainer, Float_t weight);
  virtual void FillXX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillXY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillYX(Int_t harmonic, const Float_t *variableContainer, Float_t weight);
  virtual void FillYY(Int_t harmonic, const Float_t *variableContainer, Float_t weight);

protected:
  void FillBinAxesValues(const Float_t *variableContainer, Int_t chgrpId = -1);
  THnF* DivideTHnF(THnF* values, THnI* entries);
  void CopyTHnF(THnF *hDest, THnF *hSource, Int_t *binsArray);
  void CopyTHnFDimension(THnF *hDest, THnF *hSource, Int_t *binsArray, Int_t dimension);

  AliQnCorrectionsEventClassVariablesSet fEventClassVariables;  //!<! The variables set that determines the event classes
  Double_t *fBinAxesValues;                                  //!<! Runtime place holder for computing bin number
  QnCorrectionHistogramErrorMode fErrorMode;                 //!<! The error type for the current instance
  /// \cond CLASSIMP
  ClassDef(AliQnCorrectionsHistogramBase, 1);
  /// \endcond
  static const char *szChannelAxisTitle;                 ///< The title for the channel extra axis
  static const char *szGroupAxisTitle;                   ///< The title for the channel group extra axis
  static const char *szGroupHistoPrefix;                 ///< The prefix for the name of the group histograms
  static const char *szEntriesHistoSuffix;               ///< The suffix for the name of the entries histograms
  static const char *szXComponentSuffix;                 ///< The suffix for the name of X component histograms
  static const char *szYComponentSuffix;                 ///< The suffix for the name of Y component histograms
  static const char *szXXCorrelationComponentSuffix;     ///< The suffix for the name of XX correlation component histograms
  static const char *szXYCorrelationComponentSuffix;     ///< The suffix for the name of XY correlation component histograms
  static const char *szYXCorrelationComponentSuffix;     ///< The suffix for the name of YX correlation component histograms
  static const char *szYYCorrelationComponentSuffix;     ///< The suffix for the name of YY correlation component histograms
  static const Int_t nMaxHarmonicNumberSupported;        ///< The maximum external harmonic number the framework support
  static const UInt_t harmonicNumberMask[];              ///< Mask for each external harmonic number
  static const UInt_t correlationXXmask;                 ///< Maks for XX correlation component
  static const UInt_t correlationXYmask;                 ///< Maks for XY correlation component
  static const UInt_t correlationYXmask;                 ///< Maks for YX correlation component
  static const UInt_t correlationYYmask;                 ///< Maks for YY correlation component
  static const Int_t nMinNoOfEntriesValidated;           ///< The minimum number of entries for validating a bin content
};

/// Fills the axes values for the current passed variable container
///
/// Core of the GetBin members. Stores the current values of the involved
/// variables in the internal place holder. Space is prepared for potential
/// channel or group id.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param chgrpId additional optional channel or group Id
inline void AliQnCorrectionsHistogramBase::FillBinAxesValues(const Float_t *variableContainer, Int_t chgrpId) {
  for (Int_t var = 0; var < fEventClassVariables.GetEntriesFast(); var++) {
    fBinAxesValues[var] = variableContainer[fEventClassVariables.At(var)->GetVariableId()];
  }
  fBinAxesValues[fEventClassVariables.GetEntriesFast()] = chgrpId;
}


#endif
