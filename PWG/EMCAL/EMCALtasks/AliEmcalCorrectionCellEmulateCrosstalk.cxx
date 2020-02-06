// AliEmcalCorrectionCellEmulateCrosstalk
//

#include <map>
#include <vector>
#include <string>

#include <TObjArray.h>
#include <TFile.h>

#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include <AliAODEvent.h>

#include "AliEmcalCorrectionCellEmulateCrosstalk.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellEmulateCrosstalk);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellEmulateCrosstalk> AliEmcalCorrectionCellEmulateCrosstalk::reg("AliEmcalCorrectionCellEmulateCrosstalk");

/**
 * Default constructor
 */
AliEmcalCorrectionCellEmulateCrosstalk::AliEmcalCorrectionCellEmulateCrosstalk() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellEmulateCrosstalk"),
  fCellEnergyDistBefore(0),
  fCellEnergyDistAfter(0),
  fTCardCorrClusEnerConserv(kFALSE),
  fRandom(0),
  fRandomizeTCard(kTRUE),
  fTCardCorrMinAmp(0.01),
  fTCardCorrMinInduced(0), 
  fTCardCorrMaxInducedELeak(0), 
  fTCardCorrMaxInduced(100),
  fPrintOnce(kFALSE),
  fAODCellsTmp(0x0)
{
  for(Int_t i = 0; i < fgkNsm;    i++)
  {
    fTCardCorrInduceEnerProb   [i] = 0;
    fTCardCorrInduceEnerFracMax[i] = 100;
    fTCardCorrInduceEnerFracMin[i] =-100;
    
    for(Int_t j = 0; j < 4 ; j++)
    {
      fTCardCorrInduceEner         [j][i] =  0 ;
      fTCardCorrInduceEnerFrac     [j][i] =  0 ;
      fTCardCorrInduceEnerFracP1   [j][i] =  0 ;
      fTCardCorrInduceEnerFracWidth[j][i] =  0 ;
    }
  }
  
  ResetArrays();
}

/**
 * Destructor
 */
AliEmcalCorrectionCellEmulateCrosstalk::~AliEmcalCorrectionCellEmulateCrosstalk()
{
  fAODCellsTmp->DeleteContainer();
  delete fAODCellsTmp;
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellEmulateCrosstalk::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  AliWarning("Init EMCAL crosstalk emulation");
  
  GetProperty("conservEnergy", fTCardCorrClusEnerConserv);
  GetProperty("randomizeTCardInducedEnergy", fRandomizeTCard);
  GetProperty("inducedTCardMinimumCellEnergy", fTCardCorrMinAmp);
  GetProperty("inducedTCardMaximum", fTCardCorrMaxInduced);
  GetProperty("inducedTCardMaximumELeak", fTCardCorrMaxInducedELeak);
  GetProperty("inducedTCardMinimum", fTCardCorrMinInduced);

  // Handle array initialization
  // For the format of these values, see the "default" yaml configuration file
  const std::map <std::string, Float_t (*)[fgkNsm]> properties2D = {
    {"inducedEnergyLossConstant", fTCardCorrInduceEner},
    {"inducedEnergyLossFraction", fTCardCorrInduceEnerFrac},
    {"inducedEnergyLossFractionP1", fTCardCorrInduceEnerFracP1},
    {"inducedEnergyLossFractionWidth", fTCardCorrInduceEnerFracWidth}
  };
  const std::map <std::string, Float_t *> properties1D = {
    {"inducedEnergyLossMinimumFraction", fTCardCorrInduceEnerFracMin},
    {"inducedEnergyLossMaximumFraction", fTCardCorrInduceEnerFracMax},
    {"inducedEnergyLossProbability", fTCardCorrInduceEnerProb}
  };
  RetrieveAndSetProperties(properties2D);
  RetrieveAndSetProperties(properties1D);
  
  if (!fRecoUtils) {
    fRecoUtils  = new AliEMCALRecoUtils;
  }

  GetProperty("printConfiguration", fPrintOnce);
  if (fPrintOnce) {
    PrintTCardParam();
  }

  return kTRUE;
}

/**
 * Handles assigning properties to 2D array class members.
 *
 * @param[in] val Array class member from map
 * @param[in] property Values extracted from the YAML config to be assigned to val
 * @param[in] iSM Super-module number (0 indexed)
 * @param[in] name Name of the property which is being assigned
 */
void AliEmcalCorrectionCellEmulateCrosstalk::SetProperty(Float_t val[][fgkNsm], std::vector<double> & property, unsigned int iSM, const std::string & name)
{
  for (unsigned int iProperty = 0; iProperty < property.size(); iProperty++) {
    val[iProperty][iSM] = property.at(iProperty);
  }
}

/**
 * Handles assigning properties to 1D array class members.
 *
 * @param[in] val Array class member from map
 * @param[in] property Values extracted from the YAML config to be assigned to val. It is expected to be contain one value.
 * @param[in] iSM Super-module number (0 indexed)
 * @param[in] name Name of the property which is being assigned
 */
void AliEmcalCorrectionCellEmulateCrosstalk::SetProperty(Float_t val[fgkNsm], std::vector<double> & property, unsigned int iSM, const std::string & name)
{
  if (property.size() != 1) {
    AliErrorStream() << "Trying to set single value for property " << name << ", but given " << property.size() << " values. Please check your configuration!\n";
  }
  else {
    val[iSM] = property.at(0);
  }
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellEmulateCrosstalk::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_{cell} (GeV)",7000,0,70);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_{cell} (GeV)",7000,0,70);
    fOutput->Add(fCellEnergyDistAfter);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellEmulateCrosstalk::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  if (!fEventManager.InputEvent()) {
    AliError("Event ptr = 0, returning");
    return kFALSE;
  }
 
  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  if (fCaloCells->GetNumberOfCells()<=0)
  {
    AliDebug(2, Form("Number of EMCAL cells = %d, returning", fCaloCells->GetNumberOfCells()));
    return kFALSE;
  }

  if(fCreateHisto)
    FillCellQA(fCellEnergyDistBefore); // "before" QA
  
  // CELL CROSSTALK EMULATION -------------------------------------------------------
  
  // Compute the induced cell energies by T-Card correlation emulation, ONLY MC
  MakeCellTCardCorrelation();
  
  // Add to existing cells the found induced energies in MakeCellTCardCorrelation() if new signal is larger than 10 MeV.
  AddInducedEnergiesToExistingCells();
  
  // Add new cells with found induced energies in MakeCellTCardCorrelation() if new signal is larger than 10 MeV.
  AddInducedEnergiesToNewCells();
  
  // -------------------------------------------------------
  
  if(fCreateHisto)
    FillCellQA(fCellEnergyDistAfter); // "after" QA
  
  ResetArrays();
  
  return kTRUE;
}

/**
 * Recover each cell amplitude and absId and induce energy
 * in cells in cross of the same T-Card
 */
void AliEmcalCorrectionCellEmulateCrosstalk::MakeCellTCardCorrelation()
{
  Int_t    id     = -1;
  Float_t  amp    = -1;
  
  if (fPrintOnce) {
    PrintTCardParam();
    fPrintOnce = 0;
  }
  
  Int_t nCells = fCaloCells->GetNumberOfCells();
  
  // Loop on all cells with signal
  for (Int_t icell = 0; icell < nCells; icell++)
  {
    id  = fCaloCells->GetCellNumber(icell);
    amp = fCaloCells->GetAmplitude (icell); // fCaloCells->GetCellAmplitude(id);

    if ( amp <= fTCardCorrMinAmp ) continue ;

    //
    // First get the SM, col-row of this tower
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
    fGeom->GetCellIndex(id,imod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);
    
    //
    // Determine randomly if we want to create a correlation for this cell,
    // depending the SM number of the cell
    if ( fTCardCorrInduceEnerProb[imod] < 1 )
    {  
      Float_t rand = fRandom.Uniform(0, 1);
      
      if ( rand > fTCardCorrInduceEnerProb[imod] ) continue;
    }
    
    AliDebug(1,Form("Reference cell absId %d, iEta %d, iPhi %d, amp %2.3f",id,ieta,iphi,amp));

    //
    // Get the absId of the cells in the cross and same T-Card
    Int_t absIDup = -1;
    Int_t absIDdo = -1;
    Int_t absIDlr  = -1;
    Int_t absIDuplr = -1;
    Int_t absIDdolr = -1;
    
    Int_t absIDup2 = -1;
    Int_t absIDup2lr = -1;
    Int_t absIDdo2 = -1;
    Int_t absIDdo2lr = -1;
    
    // Only 2 columns in the T-Card, +1 for even and -1 for odd with respect reference cell
    Int_t colShift = 0;
    if (  (ieta%2) && ieta <= AliEMCALGeoParams::fgkEMCALCols-1 ) colShift = -1;
    if ( !(ieta%2) && ieta >= 0 )                                 colShift = +1;
    
    absIDlr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+colShift);
    
    // Check up / down cells from reference cell not out of SM and in same T-Card
    if (  iphi < AliEMCALGeoParams::fgkEMCALRows-1 )
    {
      absIDup   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
      absIDuplr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta+colShift);
    }
    
    if (  iphi > 0 )
    {
      absIDdo   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
      absIDdolr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta+colShift);
    }
    
    // Check 2 up / 2 down cells from reference cell not out of SM
    if (  iphi < AliEMCALGeoParams::fgkEMCALRows-2 )
    {
      absIDup2   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+2, ieta);
      absIDup2lr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+2, ieta+colShift);
    }
    
    if (  iphi > 1 )
    {
      absIDdo2   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-2, ieta);
      absIDdo2lr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-2, ieta+colShift);
    }
    
    // In same T-Card?
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi+1)/8) ) { absIDup  = -1 ; absIDuplr  = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi-1)/8) ) { absIDdo  = -1 ; absIDdolr  = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi+2)/8) ) { absIDup2 = -1 ; absIDup2lr = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi-2)/8) ) { absIDdo2 = -1 ; absIDdo2lr = -1 ; }
    
    // Calculate induced energy to T-Card cells
    //
    AliDebug(1,Form("cell up %d:"  ,absIDup));
    CalculateInducedEnergyInTCardCell(absIDup   , id, imod, amp, 0);    
    AliDebug(1,Form("cell down %d:",absIDdo));
    CalculateInducedEnergyInTCardCell(absIDdo   , id, imod, amp, 0);
    
    AliDebug(1,Form("cell up left-right %d:"  ,absIDuplr));
    CalculateInducedEnergyInTCardCell(absIDuplr , id, imod, amp, 1);    
    AliDebug(1,Form("cell down left-right %d:",absIDdolr));
    CalculateInducedEnergyInTCardCell(absIDdolr , id, imod, amp, 1);
    
    AliDebug(1,Form("cell left-right %d:",absIDlr));
    CalculateInducedEnergyInTCardCell(absIDlr   , id, imod, amp, 2);
    
    AliDebug(1,Form("cell up 2nd row %d:"  ,absIDup2));
    CalculateInducedEnergyInTCardCell(absIDup2  , id, imod, amp, 3);    
    AliDebug(1,Form("cell down 2nd row %d:",absIDdo2));
    CalculateInducedEnergyInTCardCell(absIDdo2  , id, imod, amp, 3);
    
    AliDebug(1,Form("cell up left-right 2nd row %d:"  ,absIDup2lr));
    CalculateInducedEnergyInTCardCell(absIDup2lr, id, imod, amp, 3);    
    AliDebug(1,Form("cell down left-right 2nd row %d:",absIDdo2lr));
    CalculateInducedEnergyInTCardCell(absIDdo2lr, id, imod, amp, 3);
    
  } // cell loop
  
}

//_______________________________________________________
/// Calculate the induced energy in a cell belonging to the
/// same T-Card as the reference cell.
/// Used in MakeCellTCardCorrelation()
/// \param absId Id number of cell in same T-Card as reference cell
/// \param absIdRef Id number of reference cell
/// \param sm Supermodule number of cell 
/// \param ampRef Amplitude of the reference cell
/// \param cellCase Type of cell with respect reference cell 0: up or down, 1: up or down on the diagonal, 2: left or right, 3: 2nd row up/down both left/right
//_______________________________________________________
void AliEmcalCorrectionCellEmulateCrosstalk::CalculateInducedEnergyInTCardCell
(Int_t absId, Int_t absIdRef, Int_t sm, Float_t ampRef, Int_t cellCase) 
{
  // Check that the cell exists
  if ( !AcceptCell(absId) ) return ; 
  
  // Get the fraction
  Float_t frac = fTCardCorrInduceEnerFrac[cellCase][sm] + ampRef * fTCardCorrInduceEnerFracP1[cellCase][sm];
  
  // Use an absolute minimum and maximum fraction if calculated one is out of range
  if ( frac < fTCardCorrInduceEnerFracMin[sm] ) frac = fTCardCorrInduceEnerFracMin[sm];
  if ( frac > fTCardCorrInduceEnerFracMax[sm] ) frac = fTCardCorrInduceEnerFracMax[sm];   
  
  AliDebug(1,Form("\t fraction %2.3f",frac));
  
  // Randomize the induced fraction, if requested
  if ( fRandomizeTCard )
  {
    frac = fRandom.Gaus(frac, fTCardCorrInduceEnerFracWidth[cellCase][sm]);
    
    AliDebug(1,Form("\t randomized fraction %2.3f",frac));
  }
  
  // If too small or negative, do nothing else
  if ( frac < 0.0001 ) return;
  
  // Calculate induced energy
  Float_t inducedE = fTCardCorrInduceEner[cellCase][sm] + ampRef * frac;
  
  // Check if we induce too much energy, in such case use a constant value
  if ( fTCardCorrMaxInduced < inducedE ) inducedE = fTCardCorrMaxInduced;
  
  AliDebug(1,Form("\t induced E %2.3f",inducedE));
  
  // Add the induced energy, check if cell existed
  // Check that the induced+amp is large enough to avoid extra linearity effects
  // typically of the order of the clusterization cell energy cut
  // But if it is below 1 ADC, typically 10 MeV, also do it, to match Beam test linearity
  Float_t amp = fCaloCells->GetCellAmplitude(absId) ;
  if ( (amp+inducedE) > fTCardCorrMinInduced || inducedE < fTCardCorrMaxInducedELeak )
  {
    fTCardCorrCellsEner[absId] += inducedE;
    
    // If original energy of cell was null, create new one 
    if ( amp < 0.01 ) fTCardCorrCellsNew[absId] = kTRUE;
  }
  else return ;
  
  // Subtract the added energy to main cell, if energy conservation is requested
  if ( fTCardCorrClusEnerConserv )
  fTCardCorrCellsEner[absIdRef] -= inducedE;
}


/**
 * Add to existing cells the found induced energies in MakeCellTCardCorrelation() if new signal is larger than 10 MeV.
 * Need to destroy/create the default cells list and do a copy from the old to the new via a temporal arrat fAODCellsTmp
 * Not too nice or fast, but it works.
 */
void AliEmcalCorrectionCellEmulateCrosstalk::AddInducedEnergiesToExistingCells()
{
  Short_t   absId   = -1;
  Double_t  amp     = -1;
  Double_t  time    = -1;
  Int_t     mclabel = -1;
  Double_t  efrac   = 0.;
  
  // Add the induced energy to the cells and copy them into a new temporal container
  // used in AddInducedEnergiesToNewCells() to refill the default cells list fCaloCells
  // Create the data member only once. Done here, not sure where to do this properly in the framework.
  //
  if ( !fAODCellsTmp )
    fAODCellsTmp = new AliAODCaloCells("tmpCaloCellsAOD","tmpCaloCellsAOD",AliAODCaloCells::kEMCALCell);
  
  Int_t nCells = fCaloCells->GetNumberOfCells();
  fAODCellsTmp->CreateContainer(nCells);

  for (Int_t icell = 0; icell < nCells; icell++)
  {
    // Get cell
    fCaloCells->GetCell(icell, absId, amp, time, mclabel, efrac);
        
    amp+=fTCardCorrCellsEner[absId];

    // Set new amplitude in new temporal container
    fAODCellsTmp->SetCell(icell, absId, amp, time, mclabel, efrac);
  }
}
  
/**
 * Add news cells if the found induced energies in MakeCellTCardCorrelation() is larger than 10 MeV.
 */
void AliEmcalCorrectionCellEmulateCrosstalk::AddInducedEnergiesToNewCells()
{
  // Count how many new cells
  //
  Int_t nCells    = fCaloCells->GetNumberOfCells();
  Int_t nCellsNew = 0;
  for(Int_t j = 0; j < fgkNEMCalCells; j++)
  {
    // Newly created?
    if ( !fTCardCorrCellsNew[j] ) continue;
    
    // Accept only if at least 10 MeV
    if (  fTCardCorrCellsEner[j] < 0.01 ) continue;
    
    nCellsNew++;
  }
    
  // Delete default cells container and create new one with larger arrays
  // to contain new cells.
  //
  fCaloCells->DeleteContainer();
  fCaloCells->CreateContainer(nCells+nCellsNew);
  
  // Recover the original input from fAODCellsTmp filled in AddInducedEnergiesToExistingCells()
  //
  Short_t   absId      = -1;
  Float_t   amp        = -1;
  Double_t  time       =  0;
  Int_t     mclabel    = -1;
  Double_t  efrac      = -1;
  Bool_t    highgain   =  1;
  Int_t     cellNumber = -1;
  
  for(Int_t j = 0; j < nCells; j++)
  {
    absId      = fAODCellsTmp->GetCellNumber(j);
    amp        = fAODCellsTmp->GetAmplitude(j);
    time       = fAODCellsTmp->GetTime(j);
    mclabel    = fAODCellsTmp->GetMCLabel(j);
    efrac      = fAODCellsTmp->GetEFraction(j);
    highgain   = fAODCellsTmp->GetHighGain(j);
    fCaloCells->SetCell(j, absId, amp, time, mclabel, efrac, highgain);
  }
  
  // Add the new cells
  //
  for(Int_t j = 0; j < fgkNEMCalCells; j++)
  {
    // Newly created?
    if ( !fTCardCorrCellsNew[j] ) continue;
    
    // Accept only if at least 10 MeV
    if (  fTCardCorrCellsEner[j] < 0.01 ) continue;
    
    // Add new cell
    cellNumber = nCells;
    absId      = j;
    amp        = fTCardCorrCellsEner[j];
    time       = 615.*1e-9;
    mclabel    = -1;
    efrac      = 0.;

    Int_t ok = fCaloCells->SetCell(cellNumber, absId, amp, time, mclabel, efrac,1);
    
    if ( !ok ) AliError("Induced new cell could not be added!");
      
    nCells++;
  }
  
  if ( nCellsNew > 0 ) fCaloCells->Sort();  
}

/**
 * Reset arrays containing information for all possible cells.
 */
void AliEmcalCorrectionCellEmulateCrosstalk::ResetArrays()
{
  for(Int_t j = 0; j < fgkNEMCalCells; j++)
  {
    fTCardCorrCellsEner[j] = 0.;
    fTCardCorrCellsNew [j] = kFALSE;
  }

  if ( fAODCellsTmp ) fAODCellsTmp->DeleteContainer();
}

/**
 * Reject cell if acceptance criteria not passed:
 *   * correct cell number
 *
 * \param absID: absolute cell ID number
 *
 * \return bool quality of cell, exists or not
 */
Bool_t AliEmcalCorrectionCellEmulateCrosstalk::AcceptCell(Int_t absID)
{
  if ( absID < 0 || absID >= 24*48*fGeom->GetNumberOfSuperModules() )
    return kFALSE;
  
  Int_t imod = -1,iTower = -1, iIphi = -1, iIeta = -1;
  if (!fGeom->GetCellIndex(absID,imod,iTower,iIphi,iIeta))
    return kFALSE;

  return kTRUE;
}

/**
 * Print parameters for T-Card correlation emulation.
 */
void AliEmcalCorrectionCellEmulateCrosstalk::PrintTCardParam()
{
  printf("T-Card emulation activated, energy conservation <%d>, randomize E <%d>, induced energy parameters:\n",
         fTCardCorrClusEnerConserv,fRandomizeTCard);
  
  printf("T-Card emulation super-modules fraction: Min cell E %2.1f MeV; "
         "induced Min E %2.1f MeV; Max at low E %2.1f MeV; Max E %2.2f GeV\n",
         fTCardCorrMinAmp         *1000, fTCardCorrMinInduced*1000,
         fTCardCorrMaxInducedELeak*1000, fTCardCorrMaxInduced              );

  for(Int_t ism = 0; ism < fgkNsm; ism++)
  {
    printf("\t sm %d, fraction %2.3f, E frac abs min %2.3e max %2.3e \n",
        ism, fTCardCorrInduceEnerProb[ism],fTCardCorrInduceEnerFracMin[ism],fTCardCorrInduceEnerFracMax[ism]);

    for(Int_t icell = 0; icell < 4; icell++)
    {
      printf("\t \t cell type %d, c %2.4e, p0 %2.4e, p1 %2.4e, sigma %2.4e \n",
          icell,fTCardCorrInduceEner[icell][ism],fTCardCorrInduceEnerFrac[icell][ism],
          fTCardCorrInduceEnerFracP1[icell][ism],fTCardCorrInduceEnerFracWidth[icell][ism]);
    }
  }
}
