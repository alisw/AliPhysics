/// \file AliQnCorrectionsHistogramBase.cxx
/// \brief Implementation of the multidimensional profile base class

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsHistogramBase.h"
#include "AliLog.h"

const char *AliQnCorrectionsHistogramBase::szChannelAxisTitle = "Channel number";
const char *AliQnCorrectionsHistogramBase::szGroupAxisTitle = "Channels group";
const char *AliQnCorrectionsHistogramBase::szGroupHistoPrefix = "Group";
const char *AliQnCorrectionsHistogramBase::szEntriesHistoSuffix = "_entries";
const char *AliQnCorrectionsHistogramBase::szXComponentSuffix = "X";
const char *AliQnCorrectionsHistogramBase::szYComponentSuffix = "Y";
const char *AliQnCorrectionsHistogramBase::szXXCorrelationComponentSuffix = "XX";
const char *AliQnCorrectionsHistogramBase::szXYCorrelationComponentSuffix = "XY";
const char *AliQnCorrectionsHistogramBase::szYXCorrelationComponentSuffix = "YX";
const char *AliQnCorrectionsHistogramBase::szYYCorrelationComponentSuffix = "YY";
const Int_t AliQnCorrectionsHistogramBase::nMaxHarmonicNumberSupported = 15;
const UInt_t AliQnCorrectionsHistogramBase::harmonicNumberMask[] =
  {0x0000,0x0002,0x0004,0x0008,0x0010,0x0020,0x0040,0x0080,
   0x0100,0x0200,0x0400,0x0800,0x1000,0x2000,0x4000,0x8000};
const UInt_t AliQnCorrectionsHistogramBase::correlationXXmask = 0x0001;
const UInt_t AliQnCorrectionsHistogramBase::correlationXYmask = 0x0002;
const UInt_t AliQnCorrectionsHistogramBase::correlationYXmask = 0x0004;
const UInt_t AliQnCorrectionsHistogramBase::correlationYYmask = 0x0008;
const Int_t AliQnCorrectionsHistogramBase::nDefaultMinNoOfEntriesValidated = 2;

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsHistogramBase);
/// \endcond

/// Default constructor
AliQnCorrectionsHistogramBase::AliQnCorrectionsHistogramBase() :
  TNamed(),
  fEventClassVariables(),
  fBinAxesValues(NULL) {

  fErrorMode = kERRORMEAN;
  fMinNoOfEntriesToValidate = nDefaultMinNoOfEntriesValidated;
}

/// Default destructor
///
/// restores the taken memory for the bin axes values bank
AliQnCorrectionsHistogramBase::~AliQnCorrectionsHistogramBase() {
  if (fBinAxesValues != NULL)
    delete [] fBinAxesValues;
}

/// Normal constructor
///
/// Basically stores the set of variables that identify
/// the different event classes the involved histograms
/// are storing information about
///
/// This is the base class. For simplicity and consistency
/// we leave open an extra variable storage that will be
/// used by the channelized histograms.
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param ecvs the event classes variables set
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
AliQnCorrectionsHistogramBase::AliQnCorrectionsHistogramBase(const char *name,
    const char *title,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Option_t *option) :
  TNamed(name, title),
  fEventClassVariables(ecvs),
  fBinAxesValues(NULL) {

  /* one place more for storing the channel number by inherited classes */
  fBinAxesValues = new Double_t[fEventClassVariables.GetEntries() + 1];

  TString opt = option;
  opt.ToLower();
  fErrorMode = kERRORMEAN;
  if (opt.Contains("s")) fErrorMode = kERRORSPREAD;
  fMinNoOfEntriesToValidate = nDefaultMinNoOfEntriesValidated;
}

/// Attaches existing histograms as the supporting histograms
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
///
/// \param histogramList list where the histograms have to be located
/// \return the associated bin to the current variables content
Bool_t AliQnCorrectionsHistogramBase::AttachHistograms(TList *histogramList) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a channelized profile without passing channels information. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::AttachHistograms()"));
  return kFALSE;
}

/// Attaches existing histograms as the support histograms for the profile function
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
///
/// \param histogramList list where the histograms have to be located
/// \param bUsedChannel array of booleans one per each channel
/// \param nChannelGroup array of group number for each channel
/// \return true if properly attached else false
Bool_t AliQnCorrectionsHistogramBase::AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a non channelized profile but passing channels information. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::AttachHistograms()"));
  return kFALSE;
}


/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsHistogramBase::GetBin(const Float_t *variableContainer) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a channelized profile without passing the channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetBin()"));
  return -1;
}

/// Get the bin number for the current variable content and channel number
///
/// The bin number identifies the event class the current
/// variable content points to and the passed channel
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsHistogramBase::GetBin(const Float_t *variableContainer, Int_t nChannel) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a non channelized profile passing a channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetBin()"));
  return -1;
}

/// Get the bin content for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetBinContent(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetXBinContent or GetYBinContent, or GetXXBinContent ... GetYYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetBinContent()"));
  return 0.0;
}

/// Get the X component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetXBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXXBinContent ... GetYYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXBinContent()"));
  return 0.0;
}

/// Get the Y component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetYBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXXBinContent ... GetYYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYBinContent()"));
  return 0.0;
}

/// Get the bin XX component content for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetXXBinContent(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent, GetYBinContent or GetXXBinContent(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXXBinContent()"));
  return 0.0;
}

/// Get the bin XY component content for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetXYBinContent(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent, GetYBinContent or GetXYBinContent(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXYBinContent()"));
  return 0.0;
}

/// Get the bin YX component content for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetYXBinContent(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent, GetYBinContent or GetYXBinContent(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYXBinContent()"));
  return 0.0;
}

/// Get the bin YY component content for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetYYBinContent(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent, GetYBinContent or GetYYBinContent(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYYBinContent()"));
  return 0.0;
}

/// Get the bin XX component content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetXXBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent or GetYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXXBinContent()"));
  return 0.0;
}

/// Get the bin XY component content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetXYBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent or GetYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXYBinContent()"));
  return 0.0;
}

/// Get the bin YX component content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetYXBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent or GetYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYXBinContent()"));
  return 0.0;
}

/// Get the bin YY component content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsHistogramBase::GetYYBinContent(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinContent, or GetXBinContent or GetYBinContent, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYYBinContent()"));
  return 0.0;
}

/// Get the bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin number content error
Float_t AliQnCorrectionsHistogramBase::GetBinError(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetXBinError or GetYBinError, or GetXXBinError ... GetYYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetBinError()"));
  return 0.0;
}

/// Get the X component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetXBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXXBinError ... GetYYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXBinError()"));
  return 0.0;
}

/// Get the Y component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetYBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXXBinError ... GetYYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYBinError()"));
  return 0.0;
}

/// Get the XX component bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetXXBinError(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, GetXBinError, GetYBinError or GetXXBinError(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXXBinError()"));
  return 0.0;
}

/// Get the XY component bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetXYBinError(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, GetXBinError, GetYBinError or GetXYBinError(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXYBinError()"));
  return 0.0;
}

/// Get the YX component bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetYXBinError(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, GetXBinError, GetYBinError or GetYXBinError(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYXBinError()"));
  return 0.0;
}

/// Get the YY component bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetYYBinError(Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, GetXBinError, GetYBinError or GetYYBinError(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYYBinError()"));
  return 0.0;
}

/// Get the XX component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetXXBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXBinError or GetYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXXBinError()"));
  return 0.0;
}

/// Get the XY component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetXYBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXBinError or GetYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetXYBinError()"));
  return 0.0;
}

/// Get the YX component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetYXBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXBinError or GetYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYXBinError()"));
  return 0.0;
}

/// Get the YY component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t AliQnCorrectionsHistogramBase::GetYYBinError(Int_t harmonic, Long64_t bin) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   GetBinError, or GetXBinError or GetYBinError, or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::GetYYBinError()"));
  return 0.0;
}

/// Fills the histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::Fill(const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   FillX or FillY, or FillXX ... FillYY or you have instantiated a base class or you are using\n" \
      "a channelized profile without passing a channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::Fill()"));
}

/// Fills the histogram
///
/// The involved bin is computed according to the current variables
/// content and passed channel number. The bin is then increased by the given
/// weight and the entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::Fill(const Float_t *variableContainer, Int_t nChannel, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   FillX or FillY, or FillXX ... FillYY or you have instantiated a base class or you are using\n" \
      "a non channelized profile passing a channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::Fill()"));
}

/// Fills the X component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillX(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill or FillXX ... FillYY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillX()"));
}

/// Fills the Y component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillY(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill or FillXX ... FillYY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillY()"));
}

/// Fills the XX component histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillXX(const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX, FillY or FillXX(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillXX()"));
}

/// Fills the XY component histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillXY(const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX, FillY  or FillXY(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillXY()"));
}

/// Fills the YX component histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillYX(const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX, FillY  or FillYX(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillYX()"));
}

/// Fills the YY component histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillYY(const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX, FillY  or FillXX(harmonic), or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillYY()"));
}

/// Fills the XX component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillXX(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX or FillY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillXX()"));
}

/// Fills the XY component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillXY(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX or FillY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillXY()"));
}

/// Fills the YX component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillYX(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX or FillY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillYX()"));
}

/// Fills the YY component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight and the
/// entries also increased properly.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void AliQnCorrectionsHistogramBase::FillYY(Int_t harmonic, const Float_t *variableContainer, Float_t weight) {
  AliFatal(Form("You have reached base member %s. This means either you should have used\n" \
      "   Fill, FillX or FillY or you have instantiated a base class. FIX IT, PLEASE.",
      "AliQnCorrectionsHistogramBase::FillYY()"));
}

/// Divide two THn histograms
///
/// Creates a value / error multidimensional histogram from
/// a values and entries multidimensional histograms.
/// The validation histogram is filled according to entries threshold value.
/// \param hValues the values multidimensional histogram
/// \param hEntries the entries multidimensional histogram
/// \param hValid optional multidimensional histogram where validation information is stored
/// \return the values / error multidimensional histogram
THnF* AliQnCorrectionsHistogramBase::DivideTHnF(THnF *hValues, THnI *hEntries, THnC *hValid) {

  THnF *hResult =  (THnF*) THn::CreateHn(hValues->GetName(), hValues->GetTitle(), hValues);

  Double_t value, error2;
  Int_t nEntries;
  Bool_t bErrorMessage = kFALSE;
  Int_t nNotValidatedBins = 0;

  for (Long64_t bin = 0; bin < hResult->GetNbins(); bin++){
    value = hValues->GetBinContent(bin);
    nEntries = Int_t(hEntries->GetBinContent(bin));
    error2 = hValues->GetBinError2(bin);

    if (nEntries < fMinNoOfEntriesToValidate) {
      /* bin content not validated */
      hResult->SetBinContent(bin, 0.0);
      hResult->SetBinError(bin, 0.0);
      if (hValid != NULL) hValid->SetBinContent(bin, 0.0);
      if (value != 0.0) {
        bErrorMessage = kTRUE;
        nNotValidatedBins++;
      }
    }
    else {
      Double_t average = value / nEntries;
      Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
      switch (fErrorMode) {
      case kERRORMEAN:
        /* standard error on the mean of the bin values */
        hResult->SetBinContent(bin, average);
        hResult->SetBinError(bin, serror / TMath::Sqrt(nEntries));
        break;
      case kERRORSPREAD:
        /* standard deviation of the bin values */
        hResult->SetBinContent(bin, average);
        hResult->SetBinError(bin, serror);
        break;
      }
      if (hValid != NULL) hValid->SetBinContent(bin, 1.0);
    }
    hResult->SetEntries(hValues->GetEntries());
  }
  if (bErrorMessage) {
    AliError(Form("There are %d bins whose bin content were not validated! histogram: %s.\n" \
        "   Minimum number of entries to validate: %d.",
        nNotValidatedBins,
        hValues->GetName(),
        fMinNoOfEntriesToValidate));
  }
  return hResult;
}

/// Starts the copy of two THnF histograms.
///
/// Source should not have channel/group structure (axis)
/// while dest should. The channel/group involved in dest
/// should be stored in binsArray[nVariables].
/// \param hDest the histogram that will receive the copy
/// \param hSource the histogram to copy
/// \param binsArray the array to build the bin numbers on each dimension
void AliQnCorrectionsHistogramBase::CopyTHnF(THnF *hDest, THnF *hSource, Int_t *binsArray) {

  CopyTHnFDimension(hDest, hSource, binsArray, 0);
}

/// Process a dimension in the copy of two THnF histograms process.
///
/// Source should not have channel/group structure (axis)
/// while dest should. The channel/group involved in dest
/// should be stored in binsArray[nVariables].
/// It is called in a recursive way until the number of
/// event class variables (dimensions in the source) is
/// exhausted then, the corresponding bins are copied.
/// \param hDest the histogram that will receive the copy
/// \param hSource the histogram to copy
/// \param binsArray the array to build the bin numbers on each dimension
/// \param dimension the current dimension being handled
void AliQnCorrectionsHistogramBase::CopyTHnFDimension(THnF *hDest, THnF *hSource, Int_t *binsArray, Int_t dimension) {

  /* are all variables settled */
  if (dimension < fEventClassVariables.GetEntriesFast()) {
    /* no then, scan this dimension and move to the next one */
    for (Long64_t bin = 0; bin < hSource->GetAxis(dimension)->GetNbins(); bin++) {
      binsArray[dimension] = bin + 1;
      CopyTHnFDimension(hDest, hSource, binsArray, dimension +1);
    }
  }
  else {
    /* all variables have a new bin configuration, let's do the actual copy */
    Double_t value = hSource->GetBinContent(binsArray);
    Double_t error = hSource->GetBinError(binsArray);
    hDest->SetBinContent(binsArray, value);
    hDest->SetBinError(binsArray,error);
  }
}

