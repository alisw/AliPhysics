/// \file AliQnCorrectionsProfileChannelized.cxx
/// \brief Implementation of multidimensional channelized profile class

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileChannelized.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsProfileChannelized);
/// \endcond

/// Default constructor
AliQnCorrectionsProfileChannelized::AliQnCorrectionsProfileChannelized() :
    AliQnCorrectionsHistogramBase() {

  fValues = NULL;
  fEntries = NULL;
  fUsedChannel = NULL;
  fChannelGroup = NULL;
  fNoOfChannels = 0;
  fActualNoOfChannels = 0;
  fChannelMap = NULL;
}

/// Normal constructor
///
/// Stores the set of variables that identify the
/// different event classes passing them to its parent
/// and prepares the object for actual histogram
/// creation or attachment
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param ecvs the event classes variables set
/// \param nNoOfChannels the number of channels associated
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
AliQnCorrectionsProfileChannelized::AliQnCorrectionsProfileChannelized(const char *name,
    const char *title,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Int_t nNoOfChannels,
    Option_t *option) : AliQnCorrectionsHistogramBase(name, title, ecvs, option) {

  fValues = NULL;
  fEntries = NULL;
  fUsedChannel = NULL;
  fChannelGroup = NULL;
  fNoOfChannels = nNoOfChannels;
  fActualNoOfChannels = 0;
  fChannelMap = NULL;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsProfileChannelized::~AliQnCorrectionsProfileChannelized() {

  if (fUsedChannel != NULL) delete [] fUsedChannel;
  if (fChannelGroup != NULL) delete [] fChannelGroup;
  if (fChannelMap != NULL) delete [] fChannelMap;
}



/// Creates the support histograms for the profile function
///
/// Based in the event classes variables set in the parent class
/// and the channel information passed as parameters
/// the values and entries multidimensional histograms are
/// created.
///
/// Both histograms are added to the passed histogram list
///
/// The actual number of channels is stored and a mask from
/// external channel number to histogram channel number. If
/// bUsedChannel is NULL all channels
/// within fNoOfChannels are assigned to this profile.
/// If nChannelGroup is NULL all channels assigned to this
/// profile are allocated to the same group.
/// \param histogramList list where the histograms have to be added
/// \param bUsedChannel array of booleans one per each channel
/// \param nChannelGroup array of group number for each channel
/// \return true if properly created
Bool_t AliQnCorrectionsProfileChannelized::CreateProfileHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup) {
  /* let's build the histograms names and titles */
  TString histoName = GetName();
  TString histoTitle = GetTitle();
  TString entriesHistoName = GetName(); entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle(); entriesHistoTitle += szEntriesHistoSuffix;

  /* we open space for channel variable as well */
  Int_t nVariables = fEventClassVariables.GetEntriesFast();
  Double_t *minvals = new Double_t[nVariables+1];
  Double_t *maxvals = new Double_t[nVariables+1];
  Int_t *nbins = new Int_t[nVariables+1];

  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);

  /* lets consider now the channel information */
  fUsedChannel = new Bool_t[fNoOfChannels];
  fChannelGroup = new Int_t[fNoOfChannels];
  fChannelMap = new Int_t[fNoOfChannels];

  fActualNoOfChannels = 0;
  for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
    if (bUsedChannel != NULL) {
      fUsedChannel[ixChannel] = bUsedChannel[ixChannel];
    }
    else {
      fUsedChannel[ixChannel] = kTRUE;
    }
    if (nChannelGroup != NULL) {
      fChannelGroup[ixChannel] = nChannelGroup[ixChannel];
    }
    else {
      fChannelGroup[ixChannel] = 0;
    }

    if (fUsedChannel[ixChannel]) {
      fChannelMap[ixChannel] = fActualNoOfChannels;
      fActualNoOfChannels++;
    }
  }

  /* There will be a wrong external view of the channel number especially */
  /* manifested when there are holes in the channel assignment */
  /* so, lets complete the dimension information */
  /* WARNING: be aware that ROOT does not keep label information when projecting THn */
  minvals[nVariables] = -0.5;
  maxvals[nVariables] = -0.5 + fActualNoOfChannels;
  nbins[nVariables] = fActualNoOfChannels;

  /* create the values and entries multidimensional histograms */
  fValues = new THnF((const char *) histoName, (const char *) histoTitle,nVariables+1,nbins,minvals,maxvals);
  fEntries = new THnI((const char *) entriesHistoName, (const char *) entriesHistoTitle,nVariables+1,nbins,minvals,maxvals);

  /* now let's set the proper binning and label on each axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fValues->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fEntries->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fValues->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
  }

  /* and now the channel axis */
  fValues->GetAxis(nVariables)->SetTitle(szChannelAxisTitle);
  fEntries->GetAxis(nVariables)->SetTitle(szChannelAxisTitle);

  /* and now set the proper channel labels if needed */
  if (fActualNoOfChannels != fNoOfChannels) {
    for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
      if (fUsedChannel[ixChannel]) {
        fValues->GetAxis(nVariables)->SetBinLabel(fChannelMap[ixChannel]+1, Form("%d", ixChannel));
        fEntries->GetAxis(nVariables)->SetBinLabel(fChannelMap[ixChannel]+1, Form("%d", ixChannel));
      }
    }
  }

  fValues->Sumw2();

  histogramList->Add(fValues);
  histogramList->Add(fEntries);

  delete [] minvals;
  delete [] maxvals;
  delete [] nbins;

  return kTRUE;
}


/// Get the bin number for the current variable content and passed channel
///
/// The bin number identifies the event class the current
/// variable content points to under the passed channel.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsProfileChannelized::GetBin(const Float_t *variableContainer, Int_t nChannel) {

  FillBinAxesValues(variableContainer, fChannelMap[nChannel]);
  /* store the channel number */
  return fEntries->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// If the number of entries is lower
/// than the minimum number of entries to validate it
/// the bin content is not considered valid and kFALSE is returned,
/// otherwise kTRUE is returned
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t AliQnCorrectionsProfileChannelized::BinContentValidated(Long64_t bin) {
  Int_t nEntries = Int_t(fEntries->GetBinContent(bin));

  if (nEntries < fMinNoOfEntriesToValidate) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
}

/// Get the bin content for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// is requested. If the bin content is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileChannelized::GetBinContent(Long64_t bin) {

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content error
Float_t AliQnCorrectionsProfileChannelized::GetBinError(Long64_t bin) {

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fValues->GetBinContent(bin);
    Float_t error2 = fValues->GetBinError2(bin);

    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
    case kERRORMEAN:
      /* standard error on the mean of the bin values */
      return serror / TMath::Sqrt(nEntries);
      break;
    case kERRORSPREAD:
      /* standard deviation of the bin values */
      return serror;
      break;
    default:
      return 0.0;
    }
  }
}

/// Fills the histogram
///
/// The involved bin is computed according to the current variables
/// content and the passed external channel number. The bin is then
/// increased by the given weight and the entries also increased properly.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \param weight the increment in the bin content
void AliQnCorrectionsProfileChannelized::Fill(const Float_t *variableContainer, Int_t nChannel, Float_t weight) {
  /* keep the total entries in fValues updated */
  Double_t nEntries = fValues->GetEntries();

  FillBinAxesValues(variableContainer, fChannelMap[nChannel]);
  /* and now update the bin */
  fValues->Fill(fBinAxesValues, weight);
  fValues->SetEntries(nEntries + 1);
  fEntries->Fill(fBinAxesValues, 1.0);
}

