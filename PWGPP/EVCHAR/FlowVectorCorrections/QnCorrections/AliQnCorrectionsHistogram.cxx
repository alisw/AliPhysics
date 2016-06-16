/// \file AliQnCorrectionsHistogram.cxx
/// \brief Implementation of the single multidimensional histograms

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsHistogram.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsHistogram);
/// \endcond

/// Default constructor
AliQnCorrectionsHistogram::AliQnCorrectionsHistogram() :
    AliQnCorrectionsHistogramBase() {
  fValues = NULL;
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
AliQnCorrectionsHistogram::AliQnCorrectionsHistogram(const char *name,
      const char *title,
      AliQnCorrectionsEventClassVariablesSet &ecvs) :
          AliQnCorrectionsHistogramBase(name, title, ecvs) {
  fValues = NULL;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsHistogram::~AliQnCorrectionsHistogram() {

}


/// Creates the support histogram for the histogram function
///
/// Based in the event classes variables set in the parent class
/// the values multidimensional histogram is created.
///
/// The histogram is added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \return true if properly created
Bool_t AliQnCorrectionsHistogram::CreateHistogram(TList *histogramList) {
  /* let's build the histograms names and titles */
  TString histoName = GetName();
  TString histoTitle = GetTitle();
  TString entriesHistoName = GetName(); entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle(); entriesHistoTitle += szEntriesHistoSuffix;

  /* we open space for channel for event class variables */
  Int_t nVariables = fEventClassVariables.GetEntriesFast();
  Double_t *minvals = new Double_t[nVariables];
  Double_t *maxvals = new Double_t[nVariables];
  Int_t *nbins = new Int_t[nVariables];

  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);

  /* create the values multidimensional histogram */
  fValues = new THnF((const char *) histoName, (const char *) histoTitle,nVariables,nbins,minvals,maxvals);

  /* now let's set the proper binning and label on each axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fValues->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fValues->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
  }

  fValues->Sumw2();

  histogramList->Add(fValues);

  delete [] minvals;
  delete [] maxvals;
  delete [] nbins;

  return kTRUE;
}

/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to under the passed channel.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsHistogram::GetBin(const Float_t *variableContainer) {

  FillBinAxesValues(variableContainer);
  /* store the channel number */
  return fValues->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// This kind of histograms cannot validate the bin content so, it
/// is always valid.
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t AliQnCorrectionsHistogram::BinContentValidated(Long64_t) {

  return kTRUE;
}

/// Get the bin content for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// is requested.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogram::GetBinContent(Long64_t bin) {

  return fValues->GetBinContent(bin);
}

/// Get the bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// error is requested.
///
/// \param bin the interested bin number
/// \return the bin number content error
Float_t AliQnCorrectionsHistogram::GetBinError(Long64_t bin) {

  return fValues->GetBinError(bin);
}

/// Fills the histogram
///
/// The involved bin is computed according to the current variables
/// content and the passed external channel number. The bin is then
/// increased by the given weight.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogram::Fill(const Float_t *variableContainer, Float_t weight) {
  /* keep the total entries in fValues updated */
  Double_t nEntries = fValues->GetEntries();

  FillBinAxesValues(variableContainer);
  /* and now update the bin */
  fValues->Fill(fBinAxesValues, weight);
  fValues->SetEntries(nEntries + 1);
}


