// AliEmcalCorrectionCellEnergyCompression
//

#include <TGrid.h>
#include <TFile.h>

#include "AliEmcalCorrectionCellEnergyCompression.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellEnergyCompression);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellEnergyCompression> AliEmcalCorrectionCellEnergyCompression::reg("AliEmcalCorrectionCellEnergyCompression");

/**
 * Default constructor
 */
AliEmcalCorrectionCellEnergyCompression::AliEmcalCorrectionCellEnergyCompression() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellEnergyCompression"),
  fCellEnergyCompressionLoss(NULL),
  fCellTimeCompressionLoss(NULL),
  fDoEnergyCompression(false),
  fEnergy(0.),
  fTruncationMaxE(250),
  fEnergyCompression(0.0153),
  fCompressionBias(0.),
  fDoTimeCompression(false),
  fTime(0.),
  fTimeCompression(0.1),
  fTimeShift(0.5)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellEnergyCompression::~AliEmcalCorrectionCellEnergyCompression()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellEnergyCompression::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  GetProperty("DoEnergyCompression", fDoEnergyCompression);
  GetProperty("TruncationMaxE", fTruncationMaxE);
  GetProperty("EnergyCompression", fEnergyCompression);
  GetProperty("CompressionBias", fCompressionBias);

  GetProperty("DoTimeCompression", fDoTimeCompression);
  GetProperty("TimeCompression", fTimeCompression);
  GetProperty("TimeShift", fTimeShift);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellEnergyCompression::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  // Create my user objects.
  if (fCreateHisto){
    fCellEnergyCompressionLoss = new TH1F("CellEnergyCompressionLoss","CellEnergyCompressionLoss;E_{before} - E_{after} (MeV)",400, -20, 20);
    fOutput->Add(fCellEnergyCompressionLoss);
    fCellTimeCompressionLoss = new TH1F("CellTimeCompressionLoss","CellTimeCompressionLoss;t_{before} - t_{after} (ns)",400, -20, 20);
    fOutput->Add(fCellTimeCompressionLoss);

    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionCellEnergyCompression::ExecOnce()
{
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellEnergyCompression::Run()
{
  AliEmcalCorrectionComponent::Run();

  Short_t  absId  =-1;
  Double_t ecell = 0;
  Double_t tcell = 0;
  Double_t efrac = 0;
  Int_t  mclabel = -1;

  // Loop over all EMCal cells
  for (Int_t iCell = 0; iCell < fCaloCells->GetNumberOfCells(); iCell++){

    // Get cell
    Bool_t getCellResult = fCaloCells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);
    if (!getCellResult) {
      AliWarning(TString::Format("Could not get cell %i from cell collection %s", iCell, fCaloCells->GetName()));
    }

    // Get high gain attribute in addition to cell
    // NOTE: GetCellHighGain() uses the cell position, not cell index, and thus should _NOT_ be used!
    Bool_t cellHighGain = fCaloCells->GetHighGain(iCell);

    if(fDoEnergyCompression){
      Double_t ecellNoComp = ecell;
      SetEnergy(ecell);
      ecell = GetEnergy();
      if(fCreateHisto){
        fCellEnergyCompressionLoss->Fill(1000*(ecellNoComp - ecell));
      }
    }
    if(fDoTimeCompression){
      Double_t tcellNoComp = tcell;
      SetTime(tcell);
      tcell = GetTime();
      if(fCreateHisto){
        fCellTimeCompressionLoss->Fill(1e9*(tcellNoComp - tcell));
      }
    }

    if (ecell > 0.) {
      fCaloCells->SetCell(iCell, absId, ecell, tcell, mclabel, efrac, cellHighGain);
    }

  }

  return kTRUE;
}

/**
 * Set Energy. For data compression, energy is stored in int16_t (Short_t) (as in O2 DataFormats/EMCAL/...Cell.cxx)
 */
void AliEmcalCorrectionCellEnergyCompression::SetEnergy(Double_t e)
{
  double truncatedEnergy = e;
  if (truncatedEnergy < 0.) {
    truncatedEnergy = 0.;
  } else if (truncatedEnergy > fTruncationMaxE) {
    truncatedEnergy = fTruncationMaxE;
  }

  fEnergy = static_cast<Short_t>((truncatedEnergy / fEnergyCompression) + fCompressionBias);
}

/**
 * Get Energy. Reconvert to double
 */
Double_t AliEmcalCorrectionCellEnergyCompression::GetEnergy() const
{
  return fEnergy * fEnergyCompression;
}

/**
 * Set Time. For data compression, energy is stored in int16_t (Short_t) (as in O2 DataFormats/EMCAL/...Cell.cxx)
 */
void AliEmcalCorrectionCellEnergyCompression::SetTime(Double_t t)
{
  fTime = static_cast<Short_t>((t + fTimeShift) / fTimeCompression);
}

/**
* Get Time. Reconvert to double
*/
Double_t AliEmcalCorrectionCellEnergyCompression::GetTime() const
{
  return fTime * fTimeCompression;
}
