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
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_{cell} (GeV)",1000,0,10);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_{cell} (GeV)",1000,0,10);
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
    
    //
    // Check if they are not declared bad or exist
    Bool_t okup   = AcceptCell(absIDup   );
    Bool_t okdo   = AcceptCell(absIDdo   );
    Bool_t oklr   = AcceptCell(absIDlr   );
    Bool_t okuplr = AcceptCell(absIDuplr );
    Bool_t okdolr = AcceptCell(absIDdolr );
    Bool_t okup2  = AcceptCell(absIDup2  );
    Bool_t okdo2  = AcceptCell(absIDdo2  );
    Bool_t okup2lr= AcceptCell(absIDup2lr);
    Bool_t okdo2lr= AcceptCell(absIDdo2lr);
    
    AliDebug(1,Form("Same T-Card cells:\n \t up %d (%d), down %d (%d), left-right %d (%d), up-lr %d (%d), down-lr %d (%d)\n"
                    "\t up2 %d (%d), down2 %d (%d), up2-lr %d (%d), down2-lr %d (%d)",
                    absIDup ,okup ,absIDdo ,okdo ,absIDlr,oklr,absIDuplr ,okuplr ,absIDdolr ,okdolr ,
                    absIDup2,okup2,absIDdo2,okdo2,             absIDup2lr,okup2lr,absIDdo2lr,okdo2lr));

    //
    // Generate some energy for the nearby cells in same TCard , depending on this cell energy
    // Check if originally the tower had no or little energy, in which case tag it as new
    Float_t fracupdown     = fTCardCorrInduceEnerFrac[0][imod]+amp*fTCardCorrInduceEnerFracP1[0][imod];
    Float_t fracupdownleri = fTCardCorrInduceEnerFrac[1][imod]+amp*fTCardCorrInduceEnerFracP1[1][imod];
    Float_t fracleri       = fTCardCorrInduceEnerFrac[2][imod]+amp*fTCardCorrInduceEnerFracP1[2][imod];
    Float_t frac2nd        = fTCardCorrInduceEnerFrac[3][imod]+amp*fTCardCorrInduceEnerFracP1[3][imod];
    
    AliDebug(1,Form("Fraction for SM %d (min %2.3f, max %2.3f):\n"
                    "\t up-down   : c %2.3e, p1 %2.3e, p2 %2.4e, sig %2.3e, fraction %2.3f\n"
                    "\t up-down-lr: c %2.3e, p1 %2.3e, p2 %2.4e, sig %2.3e, fraction %2.3f\n"
                    "\t left-right: c %2.3e, p1 %2.3e, p2 %2.4e, sig %2.3e, fraction %2.3f\n"
                    "\t 2nd row   : c %2.3e, p1 %2.3e, p2 %2.4e, sig %2.3e, fraction %2.3f",
                    imod, fTCardCorrInduceEnerFracMin[imod], fTCardCorrInduceEnerFracMax[imod],
                    fTCardCorrInduceEner[0][imod],fTCardCorrInduceEnerFrac[0][imod],fTCardCorrInduceEnerFracP1[0][imod],fTCardCorrInduceEnerFracWidth[0][imod],fracupdown,
                    fTCardCorrInduceEner[1][imod],fTCardCorrInduceEnerFrac[1][imod],fTCardCorrInduceEnerFracP1[1][imod],fTCardCorrInduceEnerFracWidth[1][imod],fracupdownleri,
                    fTCardCorrInduceEner[2][imod],fTCardCorrInduceEnerFrac[2][imod],fTCardCorrInduceEnerFracP1[2][imod],fTCardCorrInduceEnerFracWidth[2][imod],fracleri,
                    fTCardCorrInduceEner[3][imod],fTCardCorrInduceEnerFrac[3][imod],fTCardCorrInduceEnerFracP1[3][imod],fTCardCorrInduceEnerFracWidth[3][imod],frac2nd));
    
    if( fracupdown     < fTCardCorrInduceEnerFracMin[imod] ) fracupdown     = fTCardCorrInduceEnerFracMin[imod];
    if( fracupdown     > fTCardCorrInduceEnerFracMax[imod] ) fracupdown     = fTCardCorrInduceEnerFracMax[imod];
    if( fracupdownleri < fTCardCorrInduceEnerFracMin[imod] ) fracupdownleri = fTCardCorrInduceEnerFracMin[imod];
    if( fracupdownleri > fTCardCorrInduceEnerFracMax[imod] ) fracupdownleri = fTCardCorrInduceEnerFracMax[imod];
    if( fracleri       < fTCardCorrInduceEnerFracMin[imod] ) fracleri       = fTCardCorrInduceEnerFracMin[imod];
    if( fracleri       > fTCardCorrInduceEnerFracMax[imod] ) fracleri       = fTCardCorrInduceEnerFracMax[imod];
    if( frac2nd        < fTCardCorrInduceEnerFracMin[imod] ) frac2nd        = fTCardCorrInduceEnerFracMin[imod];
    if( frac2nd        > fTCardCorrInduceEnerFracMax[imod] ) frac2nd        = fTCardCorrInduceEnerFracMax[imod];
    
    // Randomize the induced fraction, if requested
    if(fRandomizeTCard)
    {
      fracupdown     = fRandom.Gaus(fracupdown    ,fTCardCorrInduceEnerFracWidth[0][imod]);
      fracupdownleri = fRandom.Gaus(fracupdownleri,fTCardCorrInduceEnerFracWidth[1][imod]);
      fracleri       = fRandom.Gaus(fracleri      ,fTCardCorrInduceEnerFracWidth[2][imod]);
      frac2nd        = fRandom.Gaus(frac2nd       ,fTCardCorrInduceEnerFracWidth[3][imod]);
      
      AliDebug(1,Form("Randomized fraction: up-down %2.3f; up-down-left-right %2.3f; left-right %2.3f; 2nd row %2.3f",
                      fracupdown,fracupdownleri,fracleri,frac2nd));
    }
    
    // Calculate induced energy
    Float_t indEupdown     = fTCardCorrInduceEner[0][imod]+amp*fracupdown;
    Float_t indEupdownleri = fTCardCorrInduceEner[1][imod]+amp*fracupdownleri;
    Float_t indEleri       = fTCardCorrInduceEner[2][imod]+amp*fracleri;
    Float_t indE2nd        = fTCardCorrInduceEner[3][imod]+amp*frac2nd;
    
    AliDebug(1,Form("Induced energy: up-down %2.3f; up-down-left-right %2.3f; left-right %2.3f; 2nd row %2.3f",
                    indEupdown,indEupdownleri,indEleri,indE2nd));    
    
    // Check if we induce too much energy, in such case use a constant value
    if ( fTCardCorrMaxInduced < indE2nd        ) indE2nd        = fTCardCorrMaxInduced;
    if ( fTCardCorrMaxInduced < indEupdownleri ) indEupdownleri = fTCardCorrMaxInduced;
    if ( fTCardCorrMaxInduced < indEupdown     ) indEupdown     = fTCardCorrMaxInduced;
    if ( fTCardCorrMaxInduced < indEleri       ) indEleri       = fTCardCorrMaxInduced;
    
    AliDebug(1,Form("Induced energy, saturated?: up-down %2.3f; up-down-left-right %2.3f; left-right %2.3f; 2nd row %2.3f",
                    indEupdown,indEupdownleri,indEleri,indE2nd));

    //
    // Add the induced energy, check if cell existed
    if ( okup )
    {
      fTCardCorrCellsEner[absIDup] += indEupdown;
      
      if ( fCaloCells->GetCellAmplitude(absIDup) < 0.01 ) fTCardCorrCellsNew[absIDup] = kTRUE;
    }
    
    if ( okdo )
    {
      fTCardCorrCellsEner[absIDdo] += indEupdown;
      
      if ( fCaloCells->GetCellAmplitude(absIDdo) < 0.01 ) fTCardCorrCellsNew[absIDdo] = kTRUE;
    }
    
    if ( oklr )
    {
      fTCardCorrCellsEner[absIDlr] += indEleri;
      
      if ( fCaloCells->GetCellAmplitude(absIDlr) < 0.01 ) fTCardCorrCellsNew[absIDlr]  = kTRUE;
    }
    
    if ( okuplr )
    {
      fTCardCorrCellsEner[absIDuplr] += indEupdownleri;
      
      if ( fCaloCells->GetCellAmplitude(absIDuplr ) < 0.01 ) fTCardCorrCellsNew[absIDuplr]  = kTRUE;
    }
    
    if ( okdolr )
    {
      fTCardCorrCellsEner[absIDdolr] += indEupdownleri;
      
      if ( fCaloCells->GetCellAmplitude(absIDdolr ) < 0.01 ) fTCardCorrCellsNew[absIDdolr]  = kTRUE;
    }
    
    if ( okup2 )
    {
      fTCardCorrCellsEner[absIDup2] += indE2nd;
      
      if ( fCaloCells->GetCellAmplitude(absIDup2) < 0.01 ) fTCardCorrCellsNew[absIDup2] = kTRUE;
    }
    
    if ( okup2lr )
    {
      fTCardCorrCellsEner[absIDup2lr] += indE2nd;
      
      if ( fCaloCells->GetCellAmplitude(absIDup2lr) < 0.01 ) fTCardCorrCellsNew[absIDup2lr] = kTRUE;
    }
    
    if ( okdo2 )
    {
      fTCardCorrCellsEner[absIDdo2] += indE2nd;
      
      if ( fCaloCells->GetCellAmplitude(absIDdo2) < 0.01 ) fTCardCorrCellsNew[absIDdo2] = kTRUE;
    }
    
    if ( okdo2lr )
    {
      fTCardCorrCellsEner[absIDdo2lr] += indE2nd;
      
      if ( fCaloCells->GetCellAmplitude(absIDdo2lr) < 0.01 ) fTCardCorrCellsNew[absIDdo2lr] = kTRUE;
    }
    
    //
    // Subtract the added energy to main cell, if energy conservation is requested
    if ( fTCardCorrClusEnerConserv )
    {
      if ( oklr    ) fTCardCorrCellsEner[id] -= indEleri;
      if ( okuplr  ) fTCardCorrCellsEner[id] -= indEupdownleri;
      if ( okdolr  ) fTCardCorrCellsEner[id] -= indEupdownleri;
      if ( okup    ) fTCardCorrCellsEner[id] -= indEupdown;
      if ( okdo    ) fTCardCorrCellsEner[id] -= indEupdown;
      if ( okup2   ) fTCardCorrCellsEner[id] -= indE2nd;
      if ( okup2lr ) fTCardCorrCellsEner[id] -= indE2nd;
      if ( okdo2   ) fTCardCorrCellsEner[id] -= indE2nd;
      if ( okdo2lr ) fTCardCorrCellsEner[id] -= indE2nd;
    } // conserve energy
    
  } // cell loop
  
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
  printf("T-Card emulation super-modules fraction: Min cell E %2.2f Max induced E %2.2f\n",
        fTCardCorrMinAmp,fTCardCorrMaxInduced);

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
