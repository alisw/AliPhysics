/// AliEmcalEmbeddingQA
///
/// Performs QA and counting for embedded events which have passed
/// the data event selection

#include <TString.h>
#include <TList.h>

#include <AliLog.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliEmcalEmbeddingQA.h"

/// \cond CLASSIMP
ClassImp(AliEmcalEmbeddingQA);
/// \endcond

/**
 * Constructor
 */
AliEmcalEmbeddingQA::AliEmcalEmbeddingQA():
  fInitialized(false),
  fHistManager("AliEmcalEmbeddingQA")
{
}

/**
 * Initialize the embedding QA histograms. Should only be called once
 *
 * @return true if the histograms were successfully initialized (this is the value of fInitialized)
 */
bool AliEmcalEmbeddingQA::Initialize()
{
  if (fInitialized) {
    AliError("QA histograms are already initialized!");
    return fInitialized;
  }

  // Create histograms
  // Determine and define properties to create the histograms
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (!embeddingHelper) {
    AliError("The embedding helper is not available!");
    return fInitialized;
  }
  TString histName;
  TString histTitle;
  // Extract the value from embedding helper for convenience
  int nPtHardBins = embeddingHelper->GetNPtHardBins();

  // Cross section
  histName = "fHistXsectionAfterSel";
  histTitle = "Pythia Cross Section After Event Selection;p_{T} hard bin; XSection";
  fHistManager.CreateTProfile(histName, histTitle, nPtHardBins + 1, -1, nPtHardBins);

  // Trials
  histName = "fHistTrialsAfterSel";
  histTitle = "Number of Pythia Trials After Event Selecdtion;p_{T} hard bin;Trials";
  fHistManager.CreateTH1(histName, histTitle, nPtHardBins + 1, -1, nPtHardBins);

  // Pt hard spectra
  histName = "fHistPtHardAfterSel";
  histTitle = "p_{T} Hard Spectra After Event Selection;p_{T} hard;Counts";
  fHistManager.CreateTH1(histName, histTitle, 500, 0, 1000);

  fInitialized = true;
  return fInitialized;
}

/**
 * Add the QA plots to output list.
 *
 * @param list List to which the histograms should be added
 */
bool AliEmcalEmbeddingQA::AddQAPlotsToList(TList * list)
{
  if (!list) {
    AliError("Output list is null!");
    return false;
  }

  if (!fInitialized) {
    AliWarning("Embedding QA histograms are not initialized! Attempting to initialize...");
    bool res = Initialize();
    if (!fInitialized) {
      AliFatal("Unable to initialize the embedding QA hists!");
      // No need to return here, as it will crash anyway
    }
  }

  // Add all histograms to output list
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    list->Add(obj);
  }

  return true;
}

/**
 * Record the embedded event properties for accepted events. It should only be called if
 * the data event is selected. The unselected data is available in the main embedding task.
 */
void AliEmcalEmbeddingQA::RecordEmbeddedEventProperties()
{
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  double ptHardBin = embeddingHelper->GetPtHardBin();

  // Fill trials, xsec, pt hard
  fHistManager.FillTH1("fHistTrialsAfterSel", ptHardBin, embeddingHelper->GetPythiaTrials());
  fHistManager.FillProfile("fHistXsectionAfterSel", ptHardBin, embeddingHelper->GetPythiaXSection());
  fHistManager.FillTH1("fHistPtHardAfterSel", embeddingHelper->GetPythiaPtHard());
}

