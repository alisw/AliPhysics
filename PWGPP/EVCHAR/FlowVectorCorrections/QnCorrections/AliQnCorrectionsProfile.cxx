/// \file AliQnCorrectionsProfile.cxx
/// \brief Implementation of the single multidimensional profile class 

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfile.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsProfile);
/// \endcond

/// Default constructor
AliQnCorrectionsProfile::AliQnCorrectionsProfile(): AliQnCorrectionsHistogramBase() {
  fValues = NULL;
  fEntries = NULL;
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
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
AliQnCorrectionsProfile::AliQnCorrectionsProfile(const char *name,
    const char *title,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Option_t *option):
    AliQnCorrectionsHistogramBase(name, title, ecvs, option) {
  fValues = NULL;
  fEntries = NULL;
}

/// Default destructor
///
/// Does nothing because none of the members are own at destruction time
AliQnCorrectionsProfile::~AliQnCorrectionsProfile() {

}

/// Creates the support histograms for the profile function
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// Both histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \return true if properly created
Bool_t AliQnCorrectionsProfile::CreateProfileHistograms(TList *histogramList) {
  /* let's build the histograms names and titles */
  TString histoName = GetName();
  TString histoTitle = GetTitle();
  TString entriesHistoName = GetName(); entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle(); entriesHistoTitle += szEntriesHistoSuffix;

  Int_t nVariables = fEventClassVariables.GetEntriesFast();

  Double_t *minvals = new Double_t[nVariables];
  Double_t *maxvals = new Double_t[nVariables];
  Int_t *nbins = new Int_t[nVariables];

  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);

  /* create the values and entries multidimensional histograms */
  fValues = new THnF((const char *) histoName, (const char *) histoTitle,nVariables,nbins,minvals,maxvals);
  fEntries = new THnI((const char *) entriesHistoName, (const char *) entriesHistoTitle,nVariables,nbins,minvals,maxvals);

  /* now let's set the proper binning and label on each axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fValues->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fEntries->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fValues->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
  }

  fValues->Sumw2();

  histogramList->Add(fValues);
  histogramList->Add(fEntries);

  delete [] minvals;
  delete [] maxvals;
  delete [] nbins;

  return kTRUE;
}

/// Attaches existing histograms as the support histograms for the profile function
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t AliQnCorrectionsProfile::AttachHistograms(TList *histogramList) {
  /* let's build the histograms names */
  TString histoName = GetName();
  TString entriesHistoName = GetName(); entriesHistoName += szEntriesHistoSuffix;

  /* initialize. Remember we don't own the histograms */
  fEntries = NULL;
  fValues = NULL;

  fEntries = (THnI *) histogramList->FindObject((const char*) entriesHistoName);
  if (fEntries != NULL && fEntries->GetEntries() != 0) {
    fValues = (THnF *) histogramList->FindObject((const char *)histoName);
    if (fValues == NULL)
      return kFALSE;
  }
  else
    return kFALSE;

  return kTRUE;
}

/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsProfile::GetBin(const Float_t *variableContainer) {
  FillBinAxesValues(variableContainer);
  return fEntries->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// If the number of entries is lower
/// than the minimum number of entries to validate it
/// the bin content is not considered valid and kFALSE is returned,
/// otherwise kTRUE is returned
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t AliQnCorrectionsProfile::BinContentValidated(Long64_t bin) {
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
Float_t AliQnCorrectionsProfile::GetBinContent(Long64_t bin) {

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
Float_t AliQnCorrectionsProfile::GetBinError(Long64_t bin) {
  Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
  Float_t values = fValues->GetBinContent(bin);
  Float_t error2 = fValues->GetBinError2(bin);

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
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
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// \param variableContainer the current variables conten addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsProfile::Fill(const Float_t *variableContainer, Float_t weight) {
  /* keep the total entries in fValues updated */
  Double_t nEntries = fValues->GetEntries();

  FillBinAxesValues(variableContainer);
  fValues->Fill(fBinAxesValues, weight);
  fValues->SetEntries(nEntries + 1);
  fEntries->Fill(fBinAxesValues, 1.0);
}

