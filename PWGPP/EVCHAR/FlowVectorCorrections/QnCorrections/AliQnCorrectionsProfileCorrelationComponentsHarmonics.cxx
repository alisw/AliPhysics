/// \file AliQnCorrectionsProfileCorrelationComponentsHarmonics.cxx
/// \brief Implementation of the multidimensional correlation components based set of profiles with harmonic support

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileCorrelationComponentsHarmonics.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsProfileCorrelationComponentsHarmonics);
/// \endcond

/// Default constructor
AliQnCorrectionsProfileCorrelationComponentsHarmonics::AliQnCorrectionsProfileCorrelationComponentsHarmonics() :
    AliQnCorrectionsHistogramBase() {

  fXXValues = NULL;
  fXYValues = NULL;
  fYXValues = NULL;
  fYYValues = NULL;
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
  fFullFilled = 0x0000;
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
AliQnCorrectionsProfileCorrelationComponentsHarmonics::AliQnCorrectionsProfileCorrelationComponentsHarmonics(const char *name,
    const char *title,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Option_t *option) :
        AliQnCorrectionsHistogramBase(name, title, ecvs, option) {

  fXXValues = NULL;
  fXYValues = NULL;
  fYXValues = NULL;
  fYYValues = NULL;
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
  fFullFilled = 0x0000;
  fEntries = NULL;
}

/// Default destructor
///
/// Returns the only taken memory, the harmonic histograms storage,
/// the own histograms and other members are not own at destruction time
AliQnCorrectionsProfileCorrelationComponentsHarmonics::~AliQnCorrectionsProfileCorrelationComponentsHarmonics() {

  if (fXXValues != NULL)
    delete [] fXXValues;
  if (fXYValues != NULL)
    delete [] fXYValues;
  if (fYXValues != NULL)
    delete [] fYXValues;
  if (fYYValues != NULL)
    delete [] fYYValues;
}

/// Creates the XX, XY, YX, YY correlation components support histograms
/// for the profile function
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// For each harmonic number fout values histograms are created, XX,
/// XY, YX and YY. The histograms are organized to support external harmonic
/// number. By default the external harmonic number is always considered to
/// start by one. If no map is passed as parameter the external harmonic
/// numbers are considered as: 1, 2, ..., nNoOfHarmonic.
/// If the user wants a different assignment he has to provide an
/// ordered map, for instance: four harmonics with external harmonic numbers
/// 2, 4, 6 and 8 will require nNoOfHarmonics = 4 and harmonicMap = [2, 4, 6, 8].
/// The fully filled condition is computed and stored
///
/// The whole set of histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \param nNoOfHarmonics the desired number of harmonics
/// \param harmonicMap ordered array with the external number of the harmonics
/// \return true if properly created
Bool_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::CreateCorrelationComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t *harmonicMap) {
  /* let's build the histograms names and titles */
  TString histoXXName = GetName(); histoXXName += szXXCorrelationComponentSuffix;
  TString histoXYName = GetName(); histoXYName += szXYCorrelationComponentSuffix;
  TString histoYXName = GetName(); histoYXName += szYXCorrelationComponentSuffix;
  TString histoYYName = GetName(); histoYYName += szYYCorrelationComponentSuffix;
  TString histoXXTitle = GetTitle(); histoXXTitle += szXXCorrelationComponentSuffix;
  TString histoXYTitle = GetTitle(); histoXYTitle += szXYCorrelationComponentSuffix;
  TString histoYXTitle = GetTitle(); histoYXTitle += szYXCorrelationComponentSuffix;
  TString histoYYTitle = GetTitle(); histoYYTitle += szYYCorrelationComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle += szXXCorrelationComponentSuffix;
  entriesHistoTitle += szXYCorrelationComponentSuffix;
  entriesHistoTitle += szYXCorrelationComponentSuffix;
  entriesHistoTitle += szYYCorrelationComponentSuffix;
  entriesHistoTitle += szEntriesHistoSuffix;

  /* check whether within the supported harmonic range */
  Int_t nHigherHarmonic = nNoOfHarmonics;
  if (harmonicMap != NULL) {
    nHigherHarmonic = harmonicMap[nNoOfHarmonics - 1];
  }
  if (nMaxHarmonicNumberSupported < nHigherHarmonic) {
    AliFatal(Form("You requested support for harmonic %d but the highest harmonic supported by the framework is currently %d",
        nHigherHarmonic, nMaxHarmonicNumberSupported));
  }

  /* let's support the external harmonic number map */
  /* external harmonic number will always start from one */
  Int_t nNumberOfSlots = 1;
  if (harmonicMap != NULL) {
    /* the highest harmonic number within the map if any */
    nNumberOfSlots += harmonicMap[nNoOfHarmonics - 1];
  }
  else {
    nNumberOfSlots += nNoOfHarmonics;
  }

  /* now allocate the slots for the values histograms */
  fXXValues = new THnF *[nNumberOfSlots];
  fXYValues = new THnF *[nNumberOfSlots];
  fYXValues = new THnF *[nNumberOfSlots];
  fYYValues = new THnF *[nNumberOfSlots];
  /* and initiallize them */
  for (Int_t i = 0; i < nNumberOfSlots; i++) {
    fXXValues[i] = NULL;
    fXYValues[i] = NULL;
    fYXValues[i] = NULL;
    fYYValues[i] = NULL;
  }

  /* now prepare the construction of the histograms */
  Int_t nVariables = fEventClassVariables.GetEntriesFast();

  Double_t *minvals = new Double_t[nVariables];
  Double_t *maxvals = new Double_t[nVariables];
  Int_t *nbins = new Int_t[nVariables];
  TString sVariableLabels = "";

  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);

  /* create the values multidimensional histograms for each harmonic */
  Int_t currentHarmonic = 0;
  for (Int_t i = 0; i < nNoOfHarmonics; i++) {
    if (harmonicMap != NULL) {
      currentHarmonic = harmonicMap[i];
    }
    else {
      currentHarmonic++;
    }
    fXXValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoXXName, currentHarmonic),
        Form("%s h%d", (const char *) histoXXTitle, currentHarmonic),
        nVariables,nbins,minvals,maxvals);
    fXYValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoXYName, currentHarmonic),
        Form("%s h%d", (const char *) histoXYTitle, currentHarmonic),
        nVariables,nbins,minvals,maxvals);
    fYXValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoYXName, currentHarmonic),
        Form("%s h%d", (const char *) histoYXTitle, currentHarmonic),
        nVariables,nbins,minvals,maxvals);
    fYYValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoYYName, currentHarmonic),
        Form("%s h%d", (const char *) histoYYTitle, currentHarmonic),
        nVariables,nbins,minvals,maxvals);

    /* now let's set the proper binning and label on each axis */
    for (Int_t var = 0; var < nVariables; var++) {
      fXXValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fXXValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fXXValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fXXValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fXYValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fXYValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fXYValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fXYValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fYXValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fYXValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fYXValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fYXValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fYYValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fYYValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
      fYYValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      fYYValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
    }

    /* ask for square sum accumulation */
    fXXValues[currentHarmonic]->Sumw2();
    fXYValues[currentHarmonic]->Sumw2();
    fYXValues[currentHarmonic]->Sumw2();
    fYYValues[currentHarmonic]->Sumw2();

    /* and finally add the histograms to the list */
    histogramList->Add(fXXValues[currentHarmonic]);
    histogramList->Add(fXYValues[currentHarmonic]);
    histogramList->Add(fYXValues[currentHarmonic]);
    histogramList->Add(fYYValues[currentHarmonic]);

    /* and update the fully filled condition */
    fFullFilled |= harmonicNumberMask[currentHarmonic];
  }

  /* create the entries multidimensional histogram */
  fEntries = new THnI((const char *) entriesHistoName, (const char *) entriesHistoTitle,nVariables,nbins,minvals,maxvals);

  /* now let's set the proper binning and label on each entries histogram axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fEntries->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
  }

  /* and finally add the entries histogram to the list */
  histogramList->Add(fEntries);

  delete [] minvals;
  delete [] maxvals;
  delete [] nbins;

  return kTRUE;
}

/// Attaches existing histograms as the support histograms for XX, XY, YX, YY
/// correlation component of the profile function for different harmonics
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// The harmonic map is inferred from the found histograms within the list
/// that match the naming scheme.
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::AttachHistograms(TList *histogramList) {
  /* let's build the histograms names */
  TString histoXXName = GetName(); histoXXName += szXXCorrelationComponentSuffix;
  TString histoXYName = GetName(); histoXYName += szXYCorrelationComponentSuffix;
  TString histoYXName = GetName(); histoYXName += szYXCorrelationComponentSuffix;
  TString histoYYName = GetName(); histoYYName += szYYCorrelationComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;

  /* initialize. Remember we don't own the histograms */
  fEntries = NULL;
  if (fXXValues != NULL) {
    delete [] fXXValues;
    fXXValues = NULL;
  }
  if (fXYValues != NULL) {
    delete [] fXYValues;
    fXYValues = NULL;
  }
  if (fYXValues != NULL) {
    delete [] fYXValues;
    fYXValues = NULL;
  }
  if (fYYValues != NULL) {
    delete [] fYYValues;
    fYYValues = NULL;
  }
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
  fFullFilled = 0x0000;

  fEntries = (THnI *) histogramList->FindObject((const char*) entriesHistoName);
  if (fEntries != NULL && fEntries->GetEntries() != 0) {
    /* allocate enough space for the supported harmonic numbers */
    fXXValues = new THnF *[nMaxHarmonicNumberSupported + 1];
    fXYValues = new THnF *[nMaxHarmonicNumberSupported + 1];
    fYXValues = new THnF *[nMaxHarmonicNumberSupported + 1];
    fYYValues = new THnF *[nMaxHarmonicNumberSupported + 1];

    /* search the multidimensional histograms for each harmonic */
    Int_t currentHarmonic = 0;
    for (Int_t i = 0; i < nMaxHarmonicNumberSupported; i++) {
      currentHarmonic++;

      fXXValues[currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoXXName, currentHarmonic));
      fXYValues[currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoXYName, currentHarmonic));
      fYXValues[currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoYXName, currentHarmonic));
      fYYValues[currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoYYName, currentHarmonic));

      /* and update the fully filled condition whether applicable */
      if ((fXXValues[currentHarmonic]  != NULL) && (fXYValues[currentHarmonic] != NULL)
          && (fYXValues[currentHarmonic] != NULL) && (fYYValues[currentHarmonic] != NULL))
      fFullFilled |= harmonicNumberMask[currentHarmonic];
    }
  }
  else
    return kFALSE;

  /* check that we actually got something */
  if (fFullFilled != 0x0000)
    return kTRUE;
  else
    return kFALSE;
}

/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetBin(const Float_t *variableContainer) {
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
Bool_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::BinContentValidated(Long64_t bin) {
  Int_t nEntries = Int_t(fEntries->GetBinContent(bin));

  if (nEntries < fMinNoOfEntriesToValidate) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
}

/// Get the XX correlation component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetXXBinContent(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fXXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXXValues[harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XY correlation component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetXYBinContent(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fXYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXYValues[harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YX correlation component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetYXBinContent(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fYXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYXValues[harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YY correlation component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetYYBinContent(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fYYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYYValues[harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XX correlation component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetXXBinError(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fXXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXXValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fXXValues[harmonic]->GetBinError2(bin);

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

/// Get the XY correlation component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetXYBinError(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fXYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXYValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fXYValues[harmonic]->GetBinError2(bin);

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

/// Get the YX correlation component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetYXBinError(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fYXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYXValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fYXValues[harmonic]->GetBinError2(bin);

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

/// Get the YY correlation component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfileCorrelationComponentsHarmonics::GetYYBinError(Int_t harmonic, Long64_t bin) {

  /* sanity check */
  if (fYYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYYValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fYYValues[harmonic]->GetBinError2(bin);

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

/// Fills the XX correlation component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsProfileCorrelationComponentsHarmonics::FillXX(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  /* first the sanity checks */
  if (fXXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
  }

  if (fXXharmonicFillMask & harmonicNumberMask[harmonic]) {
    AliFatal(Form("Filling twice the harmonic %d before entries update in histogram %s.\n" \
        "   This means you probably have not updated the other components for this harmonic. FIX IT, PLEASE.", harmonic, GetName()));
  }

  /* now it's safe to continue */

  /* keep total entries in fValues updated */
  Double_t nEntries = fXXValues[harmonic]->GetEntries();

  FillBinAxesValues(variableContainer);
  fXXValues[harmonic]->Fill(fBinAxesValues, weight);
  fXXValues[harmonic]->SetEntries(nEntries + 1);

  /* update harmonic fill mask */
  fXXharmonicFillMask |= harmonicNumberMask[harmonic];

  /* now check if time for updating entries histogram */
  if (fXXharmonicFillMask != fFullFilled) return;
  if (fXYharmonicFillMask != fFullFilled) return;
  if (fYXharmonicFillMask != fFullFilled) return;
  if (fYYharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
}

/// Fills the XY correlation component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsProfileCorrelationComponentsHarmonics::FillXY(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  /* first the sanity checks */
  if (fXYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
  }

  if (fXYharmonicFillMask & harmonicNumberMask[harmonic]) {
    AliFatal(Form("Filling twice the harmonic %d before entries update in histogram %s.\n" \
        "   This means you probably have not updated the other components for this harmonic. FIX IT, PLEASE.", harmonic, GetName()));
  }

  /* now it's safe to continue */

  /* keep total entries in fValues updated */
  Double_t nEntries = fXYValues[harmonic]->GetEntries();

  FillBinAxesValues(variableContainer);
  fXYValues[harmonic]->Fill(fBinAxesValues, weight);
  fXYValues[harmonic]->SetEntries(nEntries + 1);

  /* update harmonic fill mask */
  fXYharmonicFillMask |= harmonicNumberMask[harmonic];

  /* now check if time for updating entries histogram */
  if (fXXharmonicFillMask != fFullFilled) return;
  if (fXYharmonicFillMask != fFullFilled) return;
  if (fYXharmonicFillMask != fFullFilled) return;
  if (fYYharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
}

/// Fills the YX correlation component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsProfileCorrelationComponentsHarmonics::FillYX(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  /* first the sanity checks */
  if (fYXValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
  }

  if (fYXharmonicFillMask & harmonicNumberMask[harmonic]) {
    AliFatal(Form("Filling twice the harmonic %d before entries update in histogram %s.\n" \
        "   This means you probably have not updated the other components for this harmonic. FIX IT, PLEASE.", harmonic, GetName()));
  }

  /* now it's safe to continue */

  /* keep total entries in fValues updated */
  Double_t nEntries = fYXValues[harmonic]->GetEntries();

  FillBinAxesValues(variableContainer);
  fYXValues[harmonic]->Fill(fBinAxesValues, weight);
  fYXValues[harmonic]->SetEntries(nEntries + 1);

  /* update harmonic fill mask */
  fYXharmonicFillMask |= harmonicNumberMask[harmonic];

  /* now check if time for updating entries histogram */
  if (fXXharmonicFillMask != fFullFilled) return;
  if (fXYharmonicFillMask != fFullFilled) return;
  if (fYXharmonicFillMask != fFullFilled) return;
  if (fYYharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
}

/// Fills the YY correlation component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsProfileCorrelationComponentsHarmonics::FillYY(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  /* first the sanity checks */
  if (fYYValues[harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d in correlation component histogram %s. FIX IT, PLEASE.", harmonic, GetName()));
  }

  if (fYYharmonicFillMask & harmonicNumberMask[harmonic]) {
    AliFatal(Form("Filling twice the harmonic %d before entries update in histogram %s.\n" \
        "   This means you probably have not updated the other components for this harmonic. FIX IT, PLEASE.", harmonic, GetName()));
  }

  /* now it's safe to continue */

  /* keep total entries in fValues updated */
  Double_t nEntries = fYYValues[harmonic]->GetEntries();

  FillBinAxesValues(variableContainer);
  fYYValues[harmonic]->Fill(fBinAxesValues, weight);
  fYYValues[harmonic]->SetEntries(nEntries + 1);

  /* update harmonic fill mask */
  fYYharmonicFillMask |= harmonicNumberMask[harmonic];

  /* now check if time for updating entries histogram */
  if (fXXharmonicFillMask != fFullFilled) return;
  if (fXYharmonicFillMask != fFullFilled) return;
  if (fYXharmonicFillMask != fFullFilled) return;
  if (fYYharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXXharmonicFillMask = 0x0000;
  fXYharmonicFillMask = 0x0000;
  fYXharmonicFillMask = 0x0000;
  fYYharmonicFillMask = 0x0000;
}

