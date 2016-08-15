/// \file AliQnCorrectionsProfile3DCorrelations.cxx
/// \brief Implementation of the multidimensional correlation components based set of profiles
/// for three Qn vectors with harmonic support

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliLog.h"

#include "AliQnCorrectionsProfile3DCorrelations.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsProfile3DCorrelations);
/// \endcond

///< the number of Qn supported
#define CORRELATIONSNOOFQNVECTORS 3

/// Default constructor
AliQnCorrectionsProfile3DCorrelations::AliQnCorrectionsProfile3DCorrelations() :
    AliQnCorrectionsHistogramBase(), fNameA(""), fNameB(""), fNameC("") {

  fXXValues = NULL;
  fXYValues = NULL;
  fYXValues = NULL;
  fYYValues = NULL;
  fEntries = NULL;
  fHarmonicMultiplier = 1;
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
/// \param nameA A detector name
/// \param nameB B detector name
/// \param nameC C detector name
/// \param ecvs the event classes variables set
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
AliQnCorrectionsProfile3DCorrelations::AliQnCorrectionsProfile3DCorrelations(const char *name,
    const char *title,
    const char *nameA,
    const char *nameB,
    const char *nameC,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Option_t *option) :
        AliQnCorrectionsHistogramBase(name, title, ecvs, option), fNameA(nameA), fNameB(nameB), fNameC(nameC) {

  fXXValues = NULL;
  fXYValues = NULL;
  fYXValues = NULL;
  fYYValues = NULL;
  fEntries = NULL;
  fHarmonicMultiplier = 1;
}

/// Default destructor
///
/// Returns the only taken memory, the harmonic histograms storage,
/// the own histograms and other members are not own at destruction time
AliQnCorrectionsProfile3DCorrelations::~AliQnCorrectionsProfile3DCorrelations() {

  if (fXXValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXXValues[ixComb] != NULL)
        delete [] fXXValues[ixComb];
    }
    delete [] fXXValues;
  }
  if (fXYValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXYValues[ixComb] != NULL)
        delete [] fXYValues[ixComb];
    }
    delete [] fXYValues;
  }
  if (fYXValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYXValues[ixComb] != NULL)
        delete [] fYXValues[ixComb];
    }
    delete [] fYXValues;
  }
  if (fYYValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYYValues[ixComb] != NULL)
        delete [] fYYValues[ixComb];
    }
    delete [] fYYValues;
  }
}

/// Creates the XX, XY, YX, YY correlation components support histograms
/// for the profile function for each combination of Qn vectors
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// For each Qn vector comibination, three,
/// and for each harmonic number four values histograms are created, XX,
/// XY, YX and YY. The histograms are organized to support external harmonic
/// number. By default the external harmonic number is always considered to
/// start by one. If no map is passed as parameter the external harmonic
/// numbers are considered as: 1, 2, ..., nNoOfHarmonic.
/// If the user wants a different assignment he has to provide an
/// ordered map, for instance: four harmonics with external harmonic numbers
/// 2, 4, 6 and 8 will require nNoOfHarmonics = 4 and harmonicMap = [2, 4, 6, 8].
/// The fully filled condition is computed and stored
///
/// The potential situation where the Qn vector has an harmonic multiplier
/// is properly supported
///
/// The whole set of histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \param nNoOfHarmonics the desired number of harmonics
/// \param nHarmonicMultiplier the multiplier for the harmonic number
/// \param harmonicMap ordered array with the external number of the harmonics
/// \return true if properly created
Bool_t AliQnCorrectionsProfile3DCorrelations::CreateCorrelationComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t nHarmonicMultiplier, Int_t *harmonicMap) {

  /* store the potential harmonic multiplier */
  fHarmonicMultiplier = nHarmonicMultiplier;
  /* for now on, everything is handled as if the multiplier were m=1 so, we only consider n */

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

  /* now allocate the slots for the values histograms for each Qn vector correlation combination */
  fXXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fXYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fYXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fYYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    fXXValues[ixComb] = new THnF *[nNumberOfSlots];
    fXYValues[ixComb] = new THnF *[nNumberOfSlots];
    fYXValues[ixComb] = new THnF *[nNumberOfSlots];
    fYYValues[ixComb] = new THnF *[nNumberOfSlots];
    /* and initiallize them */
    for (Int_t i = 0; i < nNumberOfSlots; i++) {
      fXXValues[ixComb][i] = NULL;
      fXYValues[ixComb][i] = NULL;
      fYXValues[ixComb][i] = NULL;
      fYYValues[ixComb][i] = NULL;
    }
  }

  /* now prepare the construction of the histograms */
  Int_t nVariables = fEventClassVariables.GetEntriesFast();

  Double_t *minvals = new Double_t[nVariables];
  Double_t *maxvals = new Double_t[nVariables];
  Int_t *nbins = new Int_t[nVariables];
  TString sVariableLabels = "";

  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);

  /* create the values multidimensional histograms for each Qn vector correlation combination and harmonic */
  const char *combNames[CORRELATIONSNOOFQNVECTORS] = {fNameA.Data(),fNameB.Data(),fNameC.Data() };
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    Int_t currentHarmonic = 0;
    for (Int_t i = 0; i < nNoOfHarmonics; i++) {
      if (harmonicMap != NULL) {
        currentHarmonic = harmonicMap[i];
      }
      else {
        currentHarmonic++;
      }
      /* let's build the histograms names and titles */
      TString BaseName = Form("%s %sx%s", GetName(), combNames[ixComb], combNames[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]);
      TString BaseTitle = Form("%s %sx%s", GetTitle(), combNames[ixComb], combNames[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]);
      TString histoXXName = BaseName; histoXXName += szXXCorrelationComponentSuffix;
      TString histoXYName = BaseName; histoXYName += szXYCorrelationComponentSuffix;
      TString histoYXName = BaseName; histoYXName += szYXCorrelationComponentSuffix;
      TString histoYYName = BaseName; histoYYName += szYYCorrelationComponentSuffix;
      TString histoXXTitle = BaseTitle; histoXXTitle += szXXCorrelationComponentSuffix;
      TString histoXYTitle = BaseTitle; histoXYTitle += szXYCorrelationComponentSuffix;
      TString histoYXTitle = BaseTitle; histoYXTitle += szYXCorrelationComponentSuffix;
      TString histoYYTitle = BaseTitle; histoYYTitle += szYYCorrelationComponentSuffix;

      fXXValues[ixComb][currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoXXName, currentHarmonic * fHarmonicMultiplier),
          Form("%s h%d", (const char *) histoXXTitle, currentHarmonic),
          nVariables,nbins,minvals,maxvals);
      fXYValues[ixComb][currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoXYName, currentHarmonic * fHarmonicMultiplier),
          Form("%s h%d", (const char *) histoXYTitle, currentHarmonic),
          nVariables,nbins,minvals,maxvals);
      fYXValues[ixComb][currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoYXName, currentHarmonic * fHarmonicMultiplier),
          Form("%s h%d", (const char *) histoYXTitle, currentHarmonic),
          nVariables,nbins,minvals,maxvals);
      fYYValues[ixComb][currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoYYName, currentHarmonic * fHarmonicMultiplier),
          Form("%s h%d", (const char *) histoYYTitle, currentHarmonic),
          nVariables,nbins,minvals,maxvals);

      /* now let's set the proper binning and label on each axis */
      for (Int_t var = 0; var < nVariables; var++) {
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      }

      /* ask for square sum accumulation */
      fXXValues[ixComb][currentHarmonic]->Sumw2();
      fXYValues[ixComb][currentHarmonic]->Sumw2();
      fYXValues[ixComb][currentHarmonic]->Sumw2();
      fYYValues[ixComb][currentHarmonic]->Sumw2();

      /* and finally add the histograms to the list */
      histogramList->Add(fXXValues[ixComb][currentHarmonic]);
      histogramList->Add(fXYValues[ixComb][currentHarmonic]);
      histogramList->Add(fYXValues[ixComb][currentHarmonic]);
      histogramList->Add(fYYValues[ixComb][currentHarmonic]);
    }
  }
  /* now the entries histogram name and title */
  TString entriesHistoName = GetName();
  entriesHistoName = entriesHistoName + fNameA + fNameB + fNameC;
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle = entriesHistoTitle + fNameA + fNameB + fNameC;
  entriesHistoTitle += szXXCorrelationComponentSuffix;
  entriesHistoTitle += szXYCorrelationComponentSuffix;
  entriesHistoTitle += szYXCorrelationComponentSuffix;
  entriesHistoTitle += szYYCorrelationComponentSuffix;
  entriesHistoTitle += szEntriesHistoSuffix;

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
/// and for the different Qn vector correlation combinations.
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// The harmonic map is inferred from the found histograms within the list
/// that match the naming scheme.
///
/// The potential situation where the Qn vectors have an harmonic multiplier
/// is properly supported
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t AliQnCorrectionsProfile3DCorrelations::AttachHistograms(TList *histogramList) {
  /* initialize. Remember we don't own the histograms */
  AliInfo("");
  fEntries = NULL;
  if (fXXValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXXValues[ixComb] != NULL)
        delete [] fXXValues[ixComb];
    }
    delete [] fXXValues;
  }
  if (fXYValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXYValues[ixComb] != NULL)
        delete [] fXYValues[ixComb];
    }
    delete [] fXYValues;
  }
  if (fYXValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYXValues[ixComb] != NULL)
        delete [] fYXValues[ixComb];
    }
    delete [] fYXValues;
  }
  if (fYYValues != NULL) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYYValues[ixComb] != NULL)
        delete [] fYYValues[ixComb];
    }
    delete [] fYYValues;
  }

  /* let's build the entries histogram name */
  TString entriesHistoName = GetName();
  entriesHistoName = entriesHistoName + fNameA + fNameB + fNameC;
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;

  UInt_t harmonicFilledMask = 0x0000;
  fEntries = (THnI *) histogramList->FindObject((const char*) entriesHistoName);
  if (fEntries != NULL && fEntries->GetEntries() != 0) {
    /* allocate enough space for the supported harmonic numbers */
    fXXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fXYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fYXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fYYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      fXXValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fXYValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fYXValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fYYValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
    }
    /* search the multidimensional histograms for each harmonic and Qn vecto correlation combination */
    const char *combNames[CORRELATIONSNOOFQNVECTORS] = {fNameA.Data(),fNameB.Data(),fNameC.Data() };
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      Int_t currentHarmonic = 0;
      /* let's build the histograms names */
      TString BaseName = Form("%s %sx%s", GetName(), combNames[ixComb], combNames[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]);
      TString BaseTitle = Form("%s %sx%s", GetTitle(), combNames[ixComb], combNames[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]);
      TString histoXXName = BaseName; histoXXName += szXXCorrelationComponentSuffix;
      TString histoXYName = BaseName; histoXYName += szXYCorrelationComponentSuffix;
      TString histoYXName = BaseName; histoYXName += szYXCorrelationComponentSuffix;
      TString histoYYName = BaseName; histoYYName += szYYCorrelationComponentSuffix;
      for (Int_t i = 0; i < nMaxHarmonicNumberSupported; i++) {
        currentHarmonic++;

        fXXValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoXXName, currentHarmonic * fHarmonicMultiplier));
        fXYValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoXYName, currentHarmonic * fHarmonicMultiplier));
        fYXValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoYXName, currentHarmonic * fHarmonicMultiplier));
        fYYValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoYYName, currentHarmonic * fHarmonicMultiplier));

        /* update the correcto condition */
        if ((fXXValues[ixComb][currentHarmonic]  != NULL) && (fXYValues[ixComb][currentHarmonic] != NULL)
            && (fYXValues[ixComb][currentHarmonic] != NULL) && (fYYValues[ixComb][currentHarmonic] != NULL))
          harmonicFilledMask |= harmonicNumberMask[currentHarmonic];
      }
    }
  }
  else {
    AliInfo(Form("Calibration histogram %s NOT FOUND", (const char*) entriesHistoName));
    return kFALSE;
  }

  /* check that we actually got something */
  if (harmonicFilledMask != 0x0000)
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
Long64_t AliQnCorrectionsProfile3DCorrelations::GetBin(const Float_t *variableContainer) {
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
Bool_t AliQnCorrectionsProfile3DCorrelations::BinContentValidated(Long64_t bin) {
  Int_t nEntries = Int_t(fEntries->GetBinContent(bin));

  if (nEntries < fMinNoOfEntriesToValidate) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
}

/// Get the XX correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfile3DCorrelations::GetXXBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fXXValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXXValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XY correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfile3DCorrelations::GetXYBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fXYValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXYValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YX correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfile3DCorrelations::GetYXBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fYXValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYXValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YY correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfile3DCorrelations::GetYYBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fYYValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s  in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYYValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XX correlation component bin content error for the passed bin number
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfile3DCorrelations::GetXXBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fXXValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXXValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fXXValues[ixComb][harmonic]->GetBinError2(bin);

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
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfile3DCorrelations::GetXYBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fXYValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXYValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fXYValues[ixComb][harmonic]->GetBinError2(bin);

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
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfile3DCorrelations::GetYXBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fYXValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYXValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fYXValues[ixComb][harmonic]->GetBinError2(bin);

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
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsProfile3DCorrelations::GetYYBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    AliFatal(Form("Accessing non existing Qn vector correlation combination %s. FIX IT, PLEASE.", comb));

  /* sanity check */
  if (fYYValues[ixComb][harmonic] == NULL) {
    AliFatal(Form("Accessing non allocated harmonic %d of Qn vector combination %s in correlation component histogram %s. FIX IT, PLEASE.", harmonic, comb, GetName()));
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  }
  else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYYValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fYYValues[ixComb][harmonic]->GetBinError2(bin);

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

/// Fills the correlation component for the different Qn vector correlation combinations
/// and for all handled harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the corresponding values.
/// The entries count is updated accordingly.
///
/// It is considered that the three Qn vectors have the same harmonic
/// structure including the harmonic multiplier. If this is not the case
/// and that situation should be supported this member must be modified.
/// \param QnA A Qn vector
/// \param QnB B Qn vector
/// \param QnC C Qn vector
/// \param variableContainer the current variables content addressed by var Id
void AliQnCorrectionsProfile3DCorrelations::Fill(const AliQnCorrectionsQnVector *QnA,
    const AliQnCorrectionsQnVector *QnB,
    const AliQnCorrectionsQnVector *QnC,
    const Float_t *variableContainer) {

  /* first the sanity checks */
  if (!((QnA->IsGoodQuality()) && (QnB->IsGoodQuality()) && (QnC->IsGoodQuality()))) return;
  if ((QnA->GetHarmonicMultiplier() != QnB->GetHarmonicMultiplier()) || (QnA->GetHarmonicMultiplier() != QnC->GetHarmonicMultiplier())) {
    AliFatal("Your are accessing here with Qn vectors with different harmonic multipliers. FIX IT, PLEASE.");
    return;
  }

  /* let's get the axis information */
  FillBinAxesValues(variableContainer);

  /* consider all combinations */
  const AliQnCorrectionsQnVector *combQn[CORRELATIONSNOOFQNVECTORS] = {QnA,QnB,QnC};
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    /* and all harmonics */
    Int_t nCurrentHarmonic = QnA->GetFirstHarmonic();
    while (nCurrentHarmonic != -1) {
      /* first the sanity checks */
      if (fXXValues[ixComb][nCurrentHarmonic] == NULL) {
        AliFatal(Form("Non allocated harmonic %d in 3D correlation component histogram %s. FIX IT, PLEASE.", nCurrentHarmonic, GetName()));
      }

      /* keep total entries in fValues updated */
      Double_t nXXEntries = fXXValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nXYEntries = fXYValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nYXEntries = fYXValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nYYEntries = fYYValues[ixComb][nCurrentHarmonic]->GetEntries();

      fXXValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues, combQn[ixComb]->Qx(nCurrentHarmonic) * combQn[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]->Qx(nCurrentHarmonic));
      fXYValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues, combQn[ixComb]->Qx(nCurrentHarmonic) * combQn[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]->Qy(nCurrentHarmonic));
      fYXValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues, combQn[ixComb]->Qy(nCurrentHarmonic) * combQn[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]->Qx(nCurrentHarmonic));
      fYYValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues, combQn[ixComb]->Qy(nCurrentHarmonic) * combQn[(ixComb+1)%CORRELATIONSNOOFQNVECTORS]->Qy(nCurrentHarmonic));

      fXXValues[ixComb][nCurrentHarmonic]->SetEntries(nXXEntries + 1);
      fXYValues[ixComb][nCurrentHarmonic]->SetEntries(nXYEntries + 1);
      fYXValues[ixComb][nCurrentHarmonic]->SetEntries(nYXEntries + 1);
      fYYValues[ixComb][nCurrentHarmonic]->SetEntries(nYYEntries + 1);

      nCurrentHarmonic = QnA->GetNextHarmonic(nCurrentHarmonic);
    }
  }

  /* update the profile entries */
  fEntries->Fill(fBinAxesValues, 1.0);
}
