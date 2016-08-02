#ifndef ALIQNCORRECTIONS_DETECTORCONFCHANNEL_H
#define ALIQNCORRECTIONS_DETECTORCONFCHANNEL_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsDetectorConfigurationChannels.h
/// \brief Channel detector configuration class for Q vector correction framework
///

#include "AliQnCorrectionsCorrectionsSetOnInputData.h"
#include "AliQnCorrectionsDataVectorChannelized.h"
#include "AliQnCorrectionsDetectorConfigurationBase.h"

class AliQnCorrectionsProfileComponents;

/// \class AliQnCorrectionsDetectorConfigurationChannels
/// \brief Channel detector configuration within Q vector correction framework
///
/// A channel detector within the Q vector correction framework is defined
/// as one for which its data vectors involve azimuthal angles and channels
/// susceptible of weighting and / or grouping and / or calibration, etc.
///
/// According to that, the proper channelized data vector is used and an extra
/// Q vector builder is incorporated.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class AliQnCorrectionsDetectorConfigurationChannels :
    public AliQnCorrectionsDetectorConfigurationBase {
public:
  friend class AliQnCorrectionsCorrectionStepBase;
  friend class AliQnCorrectionsDetector;
  AliQnCorrectionsDetectorConfigurationChannels();
  AliQnCorrectionsDetectorConfigurationChannels(const char *name,
      AliQnCorrectionsEventClassVariablesSet *eventClassesVariables,
      Int_t nNoOfChannels,
      Int_t nNoOfHarmonics,
      Int_t *harmonicMap = NULL);
  virtual ~AliQnCorrectionsDetectorConfigurationChannels();

  /// Gets the number of channels
  /// \return the number of channels of the associated detector
  Int_t GetNoOfChannels() { return fNoOfChannels; }
  /// Gets the used channels mask
  /// \return the used channels mask
  const Bool_t *GetUsedChannelsMask() const { return fUsedChannel; }
  /// Gets the channels groups
  /// \return the group associated to each channel
  const Int_t *GetChannelsGroups() const { return fChannelGroup; }
  /// Gets the hard coded group weights
  /// \return the groups hard coded weights
  const Float_t *GetHardCodedGroupWeights() const { return fHardCodedGroupWeights; }
  /// Get if the detector configuration is own by a tracking detector
  /// \return FALSE, this is a hit / channel detector configuration
  virtual Bool_t GetIsTrackingDetector() const
  { return kFALSE; }


  void SetChannelsScheme(Bool_t *bUsedChannel, Int_t *nChannelGroup, Float_t *hardCodedGroupWeights = NULL);

  /* QA section */
  /// Sets the variable id used for centrality in QA histograms.
  /// It must be one of the event class variables.
  /// \param id id for the variable used for centrality
  void SetQACentralityVar(Int_t id) { fQACentralityVarId = id; }
  /// Sets the characteristics of the multiplicity axis in QA histograms
  /// \param nbins the number of bins
  /// \param min minimum multiplicity value
  /// \param max maximum multiplicity value
  void SetQAMultiplicityAxis(Int_t nbins, Float_t min, Float_t max)
  { fQAnBinsMultiplicity = nbins; fQAMultiplicityMin = min; fQAMultiplicityMax = max; }

  virtual void AttachCorrectionsManager(AliQnCorrectionsManager *manager);

  virtual void CreateSupportDataStructures();
  virtual Bool_t CreateSupportHistograms(TList *list);
  virtual Bool_t CreateQAHistograms(TList *list);
  virtual Bool_t CreateNveQAHistograms(TList *list);

  /// Activate the processing for the passed harmonic
  /// \param harmonic the desired harmonic number to activate
  virtual void ActivateHarmonic(Int_t harmonic)
  { AliQnCorrectionsDetectorConfigurationBase::ActivateHarmonic(harmonic); fRawQnVector.ActivateHarmonic(harmonic); }
  virtual Bool_t AttachCorrectionInputs(TList *list);
  virtual void AfterInputsAttachActions();
  virtual Bool_t ProcessCorrections(const Float_t *variableContainer);
  virtual Bool_t ProcessDataCollection(const Float_t *variableContainer);

  virtual void AddCorrectionOnInputData(AliQnCorrectionsCorrectionOnInputData *correctionOnInputData);

  virtual Bool_t AddDataVector(const Float_t *variableContainer, Double_t phi, Double_t weight, Int_t channelId);

  virtual void BuildQnVector();
  void BuildRawQnVector();
  virtual void IncludeQnVectors(TList *list);
  virtual void FillOverallInputCorrectionStepList(TList *list) const;
  virtual void FillOverallQnVectorCorrectionStepList(TList *list) const;
  virtual void ReportOnCorrections(TList *steps, TList *calib, TList *apply) const;

  /// Checks if the current content of the variable bank applies to
  /// the detector configuration for the passed channel.
  /// \param variableContainer pointer to the variable content bank
  /// \param nChannel the interested external channel number
  /// \return kTRUE if the current content applies to the configuration
  virtual Bool_t IsSelected(const Float_t *variableContainer, Int_t nChannel)
    { return ((fUsedChannel[nChannel]) ? ((fCuts != NULL) ? fCuts->IsSelected(variableContainer) : kTRUE) : kFALSE); }
  /// wrong call for this class invoke base class behavior
  virtual Bool_t IsSelected(const Float_t *variableContainer)
  { return AliQnCorrectionsDetectorConfigurationBase::IsSelected(variableContainer); }

  virtual void ClearConfiguration();

private:
  static const char *szRawQnVectorName;   ///< the name of the raw Qn vector from raw data without input data corrections
  AliQnCorrectionsQnVector fRawQnVector;     ///< Q vector from input data before pre-processing
  Int_t fNoOfChannels;                    ///< The number of channels associated
  /// array, which of the detector channels is used for this configuration
  Bool_t *fUsedChannel;                   //[fNoOfChannels]
  /// array, mapping external to internal channel id. TODO: this has to go to a more generic histogram support
  Int_t *fChannelMap;                     //[fNoOfChannels]
  /// array, the group to which the channel pertains
  Int_t *fChannelGroup;                   //[fNoOfChannels]
  /// array, group hard coded weight
  Float_t *fHardCodedGroupWeights;         //[fNoOfChannels]
  AliQnCorrectionsCorrectionsSetOnInputData fInputDataCorrections; ///< set of corrections to apply on input data vectors

  /* QA section */
  void FillQAHistograms(const Float_t *variableContainer);
  static const char *szQAMultiplicityHistoName; ///< QA multiplicity histograms name
  static const char *szQAQnAverageHistogramName; ///< name and title for plain Qn vector components average QA histograms
  Int_t fQACentralityVarId;   ///< the id of the variable used for centrality in QA histograms
  Int_t fQAnBinsMultiplicity; ///< number of bins for multiplicity in QA histograms
  Float_t fQAMultiplicityMin; ///< minimum multiplicity value
  Float_t fQAMultiplicityMax; ///< maximum multiplicity value
  TH3F *fQAMultiplicityBefore3D; //!<! 3D channel multiplicity histogram before input equalization
  TH3F *fQAMultiplicityAfter3D;  //!<! 3D channel multiplicity histogram after input equalization
  AliQnCorrectionsProfileComponents *fQAQnAverageHistogram; //!<! the plain average Qn components QA histogram

private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliQnCorrectionsDetectorConfigurationChannels(const AliQnCorrectionsDetectorConfigurationChannels &);
  /// Assignment operator
  /// Not allowed. Forced private.
  AliQnCorrectionsDetectorConfigurationChannels& operator= (const AliQnCorrectionsDetectorConfigurationChannels &);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsDetectorConfigurationChannels, 2);
/// \endcond
};

/// New data vector for the detector configuration.
/// A check is made to match the channel Id with the ones assigned
/// to the detector configuration and then an additional one to
/// see if the current variable bank content passes
/// the associated cuts. If so, the data vector is stored.
/// \param variableContainer pointer to the variable content bank
/// \param phi azimuthal angle
/// \param weight the weight of the data vector. Ignored for track detector configurations.
/// \param channelId the channel Id that originates the data vector. Ignored for track detector configurations.
/// \return kTRUE if the data vector was accepted and stored
inline Bool_t AliQnCorrectionsDetectorConfigurationChannels::AddDataVector(
    const Float_t *variableContainer, Double_t phi, Double_t weight, Int_t channelId) {
  if (IsSelected(variableContainer, channelId)) {
    /// add the data vector to the bank
    new (fDataVectorBank->ConstructedAt(fDataVectorBank->GetEntriesFast()))
      AliQnCorrectionsDataVectorChannelized(channelId, phi, weight);
    return kTRUE;
  }
  return kFALSE;
}

/// Builds raw Qn vector before Q vector corrections and before input
/// data corrections but considering the chosen calibration method.
/// This is a channelized configuration so this Q vector will NOT be
/// the one to be used for subsequent Q vector corrections.
inline void AliQnCorrectionsDetectorConfigurationChannels::BuildRawQnVector() {
  fTempQnVector.Reset();

  for(Int_t ixData = 0; ixData < fDataVectorBank->GetEntriesFast(); ixData++){
    AliQnCorrectionsDataVectorChannelized *dataVector = static_cast<AliQnCorrectionsDataVectorChannelized *>(fDataVectorBank->At(ixData));
    fTempQnVector.Add(dataVector->Phi(), dataVector->Weight());
  }
  fTempQnVector.CheckQuality();
  fTempQnVector.Normalize(fQnNormalizationMethod);
  fRawQnVector.Set(&fTempQnVector, kFALSE);
}

/// Builds Qn vector before Q vector corrections and after input corrections
/// and considering the chosen calibration method.
/// The built Q vector is the one to be used for
/// subsequent Q vector corrections.
inline void AliQnCorrectionsDetectorConfigurationChannels::BuildQnVector() {
  fTempQnVector.Reset();
  fTempQ2nVector.Reset();

  for(Int_t ixData = 0; ixData < fDataVectorBank->GetEntriesFast(); ixData++){
    AliQnCorrectionsDataVectorChannelized *dataVector = static_cast<AliQnCorrectionsDataVectorChannelized *>(fDataVectorBank->At(ixData));
    fTempQnVector.Add(dataVector->Phi(), dataVector->EqualizedWeight());
    fTempQ2nVector.Add(dataVector->Phi(), dataVector->EqualizedWeight());
  }
  fTempQnVector.CheckQuality();
  fTempQ2nVector.CheckQuality();
  fTempQnVector.Normalize(fQnNormalizationMethod);
  fTempQ2nVector.Normalize(fQnNormalizationMethod);
  fPlainQnVector.Set(&fTempQnVector, kFALSE);
  fPlainQ2nVector.Set(&fTempQ2nVector, kFALSE);
  fCorrectedQnVector.Set(&fTempQnVector, kFALSE);
  fCorrectedQ2nVector.Set(&fTempQ2nVector, kFALSE);
}


/// Ask for processing corrections for the involved detector configuration
///
/// The request is transmitted to the incoming data correction steps
/// and then to Q vector correction steps.
/// The first not applied correction step breaks the loop and kFALSE is returned
/// \return kTRUE if all correction steps were applied
inline Bool_t AliQnCorrectionsDetectorConfigurationChannels::ProcessCorrections(const Float_t *variableContainer) {

  /* first we build the raw Q vector with the chosen calibration */
  BuildRawQnVector();

  /* then we transfer the request to the input data correction steps */
  for (Int_t ixCorrection = 0; ixCorrection < fInputDataCorrections.GetEntries(); ixCorrection++) {
    if (fInputDataCorrections.At(ixCorrection)->ProcessCorrections(variableContainer))
      continue;
    else
      return kFALSE;
  }

  /* input corrections were applied so let's build the Q vector with the chosen calibration */
  BuildQnVector();

  /* now let's propagate it to Q vector corrections */
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    if (fQnVectorCorrections.At(ixCorrection)->ProcessCorrections(variableContainer))
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

/// Ask for processing corrections data collection for the involved detector configuration
///
/// The request is transmitted to the incoming data correction steps
/// and then to Q vector correction steps.
/// The first not applied correction step should break the loop after collecting the data and kFALSE is returned
/// \return kTRUE if all correction steps were applied
inline Bool_t AliQnCorrectionsDetectorConfigurationChannels::ProcessDataCollection(const Float_t *variableContainer) {

  /* we transfer the request to the input data correction steps */
  for (Int_t ixCorrection = 0; ixCorrection < fInputDataCorrections.GetEntries(); ixCorrection++) {
    if (fInputDataCorrections.At(ixCorrection)->ProcessDataCollection(variableContainer))
      continue;
    else
      return kFALSE;
  }

  /* check whether QA histograms must be filled */
  FillQAHistograms(variableContainer);

  /* now let's propagate it to Q vector corrections */
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    if (fQnVectorCorrections.At(ixCorrection)->ProcessDataCollection(variableContainer))
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

/// Clean the configuration to accept a new event
///
/// Transfers the order to the Q vector correction steps then
/// to the input data correction steps and finally
/// cleans the own Q vector and the input data vector bank
/// for accepting the next event.
inline void AliQnCorrectionsDetectorConfigurationChannels::ClearConfiguration() {
  /* transfer the order to the Q vector corrections */
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    fQnVectorCorrections.At(ixCorrection)->ClearCorrectionStep();
  }
  /* transfer the order to the data vector corrections */
  for (Int_t ixCorrection = 0; ixCorrection < fInputDataCorrections.GetEntries(); ixCorrection++) {
    fInputDataCorrections.At(ixCorrection)->ClearCorrectionStep();
  }
  /* clean the raw Q vector */
  fRawQnVector.Reset();
  /* clean the own Q vector */
  fPlainQnVector.Reset();
  fPlainQ2nVector.Reset();
  fCorrectedQnVector.Reset();
  fCorrectedQ2nVector.Reset();
  /* and now clear the the input data bank */
  fDataVectorBank->Clear("C");
}

#endif // ALIQNCORRECTIONS_DETECTORCONFCHANNEL_H
