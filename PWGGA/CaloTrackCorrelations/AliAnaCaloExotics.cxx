/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- ROOT system ---
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TH3F.h>
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaCaloExotics.h"
#include "AliCaloTrackReader.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliMCAnalysisUtils.h"
#include "TCustomBinning.h"

// --- Detectors --- 
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaCaloExotics) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters.
/// Init histogram arrays to 0.
//__________________________________________
AliAnaCaloExotics::AliAnaCaloExotics() :
AliAnaCaloTrackCorrBaseClass(),  

fCellAmpMin(),                         fEMinForExo(0),                         
fExoCut(0),                            fNCellHighCut(0),
fTimeCutMin(-10000),                   fTimeCutMax(10000),
fHighEnergyCutSM(0),                   fHighNCellsCutSM(0), 
fLowEnergyCutSM3(0),                   fLowNCellsCutSM3(0),
fEventMaxNumberOfStrips(0),
fLowEnergyCutSM3Strip(0),              fLowNCellsCutSM3Strip(0),
fCellEnMax(0),                         fConstantTimeShift(0),                 
fClusterMomentum(),

fLED20(0),                             fLED12(0),   
fLED20Time(0),                         fLED12Time(0),
fAcceptEvent(1),                       fEventNStripActive(0),

fFillCellHisto(1),                     fFillAllCellEventParamHisto(1),
fFill1CellHisto(0),    
fFillStripHisto(1),
fFillMatchingHisto(0),                 fFillSameDiffFracHisto(0),
fFillExoEnMinCut(0),                   fFillAllCellSameTCardHisto(0),
fFillPerSMHisto(0),                    fFillClusterColRowHisto(0),
fFillOpenTimeHisto(1),                 fFillExo50ns(0),
fFillClusterHistoAfterEventCut(1),     fFillTrackHistoAfterEventCut(0),

// Histograms
fhNClusterPerEventNCellHigh20(0),      fhNClusterPerEventNCellHigh12(0),        
fhNClusterPerEventExotic(0),           
fhNClusterPerEventExotic1Cell(0),      fhNClusterPerEventExoticNCell(0),

fh2NClusterPerEventNCellHigh20(0),     fh2NClusterPerEventNCellHigh12(0),        
fh2NClusterPerEventExotic(0),           
fh2NClusterPerEventExotic1Cell(0),     fh2NClusterPerEventExoticNCell(0),
fh2NClusterPerEventExoticAmpMax(0),

fhExoticityEClus(0),                   
fhExoticityEClusPerSM(0),              fhExoticityEMaxCell(0),
fhExoticityEClusTrackMatch(0),         fhExoticity1Cell(0),

fhNCellsPerCluster(0),                 fhNCellsPerClusterW(0),                 fhNCellsPerClusterTimeDiff(0),  
fhNCellsPerClusterEMaxCell(0),         fhNCellsPerClusterWEMaxCell(0),   
fhNCellsPerClusterOpenTime(0),         fhNCellsPerClusterWOpenTime(0),
fhNCellsPerClusterEMaxCellOpenTime(0), fhNCellsPerClusterWEMaxCellOpenTime(0),  
fhNCellsPerClusterExo(0),              
fhNCellsPerClusterPerSM(0),            fhNCellsPerClusterWPerSM(0),              
fhNCellsPerClusterTrackMatch(0),       fhNCellsPerClusterExoTrackMatch(0),    
fhNCellsPerClusterM02(0),

fhEtaPhiGridExoEnCut(0),                              
fhEtaPhiGridEnExoCut(0),               fhEtaPhiGridEn1Cell(0),
fhEtaPhiGridEnHighNCells(0),           fhEtaPhiGridNCellEnCut(0),

fhTimeEnergyExo(0),                    fhTimeEnergy1Cell(0),                
fhTimeDiffClusCellExo(0),              fhTimeDiffClusCellDiffTCardExo(0),      fhTimeDiffClusCellSameTCardExo(0),
fhTimeDiffWClusCellExo(0),             fhTimeDiffAmpClusCellExo(0), 
fhTimeEnergyM02(0),                    fhTimeDiffClusCellM02(0),  
fhTimeEnergyNCells(0),                 fhTimeEnergyNCellsW(0),                 
fhTimeNCellCut(0),                

fhM02EnergyExo(0),                     fhM20EnergyExoM02MinCut(0),                                                 

// Other ncell
fhNCellsPerClusterSame(0),             fhNCellsPerClusterDiff(0),
fhNCellsPerClusterSame5(0),            fhNCellsPerClusterDiff5(0),
fhNCellsPerClusterSameW (0),           fhNCellsPerClusterDiffW (0),  
fhNCellsPerClusterSameDiff(0),         fhNCellsPerClusterSameDiffW(0),         fhNCellsPerClusterSameDiffTimeDiff(0),       
fhNCellsPerClusterSameFrac(0),         fhNCellsPerClusterSameFracExo(0),
fhNCellsPerClusterSameFracW(0),        fhNCellsPerClusterSameFracWExo(0),

// Other Exoticity definitions
fhExoSame(0),                          fhExoDiff(0), 
fhExoSame5(0),                         fhExoDiff5(0),

fhFracEnDiffSame(0),                   fhFracNCellDiffSame(0),                 fhFracEnNCellDiffSame(0), 
fhFracEnDiffSameW(0),                  fhFracNCellDiffSameW(0),                fhFracEnNCellDiffSameW(0), 
fhFracEnDiffSame5(0),                  fhFracNCellDiffSame5(0),                fhFracEnNCellDiffSame5(0),

fhFracEnDiffSameExo(0),                fhFracNCellDiffSameExo(0),              fhFracEnNCellDiffSameExo(0), 
fhFracEnDiffSameWExo(0),               fhFracNCellDiffSameWExo(0),             fhFracEnNCellDiffSameWExo(0), 
fhFracEnDiffSame5Exo(0),               fhFracNCellDiffSame5Exo(0),             fhFracEnNCellDiffSame5Exo(0),

fhFracEnDiffSameEnCut(0),              fhFracNCellDiffSameEnCut(0),            fhFracEnNCellDiffSameEnCut(0), 
fhFracEnDiffSameWEnCut(0),             fhFracNCellDiffSameWEnCut(0),           fhFracEnNCellDiffSameWEnCut(0), 
fhFracEnDiffSame5EnCut(0),             fhFracNCellDiffSame5EnCut(0),           fhFracEnNCellDiffSame5EnCut(0),

fhCellEnSameExo(0),                    fhCellEnDiffExo(0),
fhCellEnNCellWOpenTime(0),             fhCellEnNCellW(0),
fhCellEnNCellWEMaxOpenTime(0),         fhCellEnNCellWEMax(0),
fhCellTimeDiffNCellWOpenTime(0),       fhCellTimeDiffNCellW(0),
fhCellTimeDiffNCellWEMaxOpenTime(0),   fhCellTimeDiffNCellWEMax(0),

fhCellMaxClusterEnOpenTime(0),         fhCellMaxClusterEn(0),
fhCellMaxClusterEnRatioOpenTime(0),    fhCellMaxClusterEnRatio(0),
fhCellMaxClusterEnRatioNCellWOpenTime(0), fhCellMaxClusterEnRatioNCellW(0),
fhCellMaxClusterEnRatioExo(0),

// Exoticity for different E min cut
fhExoticityWEClus(0),                  fhNCellsPerClusterExoW(0),             
fhTimeEnergyExoW(0),                   fhM02EnergyExoW(0),  
fhExoticityECellMinCut(0),             

fhNCellsPerClusterMinEnCut(0),
fhNCellsPerClusterDiffMinEnCut(0),     fhNCellsPerClusterSameMinEnCut(0),

// All cells in same T-Card
fhExoticityEClusAllSameTCard(0),       fhM02EnergyAllSameTCard(0),
fhNCellsPerClusterAllSameTCard(0),     fhEtaPhiGridExoEnCutSameFracCut(0),
fhExoticityEClusAllSameTCardW(0),      fhM02EnergyAllSameTCardW(0),
fhNCellsPerClusterAllSameTCardW(0),    fhEtaPhiGridExoEnCutSameFracCutW(0),
fhExoticityEClusAllSameTCardMinEnCut(0), fhM02EnergyAllSameTCardMinEnCut(0),
fhNCellsPerClusterAllSameTCardMinEnCut(0),

// 50 ns cut
          
fhExoticity50nsEClus(0),               fhExoticity50ns1Cell(0),
fhExoticity50nsEClusAllSameTCard(0),   
fhNCellsPerClusterExo50ns(0),          fhM02EnergyExo50ns(0),

// Track matching vs exoticity
fhTrackMatchedDEtaNegExo(0),           fhTrackMatchedDPhiNegExo(0),            fhTrackMatchedDEtaDPhiNegExo(0),
fhTrackMatchedDEtaPosExo(0),           fhTrackMatchedDPhiPosExo(0),            fhTrackMatchedDEtaDPhiPosExo(0),
fhEOverPExo(0),

// Track matching of 1 cell clusters
fhTrackMatchedDEtaNeg1Cell(0),         fhTrackMatchedDPhiNeg1Cell(0),          fhTrackMatchedDEtaDPhiNeg1Cell(0),
fhTrackMatchedDEtaPos1Cell(0),         fhTrackMatchedDPhiPos1Cell(0),          fhTrackMatchedDEtaDPhiPos1Cell(0),
fhEOverP1Cell(0),

// Cells
fhCellExoAmp(0),                                               
fhCellExoAmpTime(0),                    fhCellExoGrid(0),    
fhCellGridTimeHighNCell20(0),           fhCellGridTimeHighNCell12(0),
fhNStripsPerEvent(0),                   fhNStripsPerEventPerSM(0),
fhNStripsPerEventAccept(0),             fhNStripsPerEventAcceptPerSM(0),
fhSM3NCellsSumEnSuspiciousEvents(0),
fhSumEnCellsPerSMEventSuspicious(0),    fhNCellsPerSMEventSuspicious(0),
fhNStripsPerEventSuspicious(0),         fhNStripsPerEventSuspiciousPerSM(0)
{        
  AddToHistogramsName("AnaCaloExotic_");
  
  // Init main cuts
  //
  fCellAmpMin = 0.5;
  fEMinForExo = 10.0;
  fExoCut     = 0.97;
  fNCellHighCut = 20;
  
  fEnergyBins [0] =   7; fEnergyBins [1] =  10;  fEnergyBins [2] =  16;
  fEnergyBins [3] =  22; fEnergyBins [4] =  30;  fEnergyBins [5] =  50;
  fEnergyBins [6] =  75; fEnergyBins [7] = 100;  fEnergyBins [8] = 125;
  fEnergyBins [9] = 150; fEnergyBins[10] = 175;  fEnergyBins[11] = 200;
  fEnergyBins [12]= 250; fEnergyBins[13] = 300;  
  
  fCellEnMins[0] = 0.5;  fCellEnMins[1] = 1.0; fCellEnMins[2] = 2.0;
  fCellEnMax = 15;
  
  fHighEnergyCutSM = 500.; fHighNCellsCutSM = 100;
  fLowEnergyCutSM3 = 20. ; fLowNCellsCutSM3 = 20;
  
  fEventMaxNumberOfStrips = 0; 
  fHighEnergyCutStrip[0] = 80; fHighEnergyCutStrip[1] = 55; 
  fHighNCellsCutStrip[0] = 24; fHighNCellsCutStrip[1] = 15;
  fLowEnergyCutSM3Strip  = 100; // open
  fLowNCellsCutSM3Strip  = 100; // open
 
  // Init to zero
  //
  for(Int_t i = 0; i < fgkNEBins; i++) 
  {
    fhM02ExoNCells        [i] = 0;
    fhM02ExoNCellsNotAllSameTCard[i] = 0;
    fhM02Exo50nsNCells    [i] = 0;
    fhClusterColRowExo [0][i] = 0;
    fhClusterColRowExo [1][i] = 0;  
//  fhClusterColRowExoW[0][i] = 0;
//  fhClusterColRowExoW[1][i] = 0;  
//  fhClusterColRow    [0][i] = 0;
//  fhClusterColRow    [1][i] = 0;
    fhClusterColRowPerSMHighNCell[i] = 0;
    fhNCellsSameDiffExo   [i] = 0;  
    fhEnSameDiffExo       [i] = 0; 
    fhEnNCellsSameDiffExo [i] = 0;  
  }
  
  for(Int_t i = 0; i < 20; i++)  
  {
    fhNCellsPerClusterExoPerSM[i] = 0;
    fEventNStripActiveSM      [i] = 0;
  }
  
  for(Int_t iemin = 0; iemin < fgkNCellEnMinBins; iemin++)
  {
    fhCellGridTime              [iemin] = 0;
    fhSumEnCells                [iemin] = 0;                     
    fhNCells                    [iemin] = 0;
    fhAverSumEnCells            [iemin] = 0;
    fhSumEnCellsNHigh20         [iemin] = 0;                     
    fhNCellsNHigh20             [iemin] = 0;
    fhAverSumEnCellsNHigh20     [iemin] = 0;
    fhSumEnCellsAcceptEvent     [iemin] = 0;                     
    fhNCellsAcceptEvent         [iemin] = 0;
    fhAverSumEnCellsAcceptEvent [iemin] = 0;
    fhSumEnCellsPerSM           [iemin] = 0;                     
    fhNCellsPerSM               [iemin] = 0;
    fhAverSumEnCellsPerSM       [iemin] = 0;
    fhSumEnCellsPerSMNHigh20    [iemin] = 0;                     
    fhNCellsPerSMNHigh20        [iemin] = 0;
    fhAverSumEnCellsPerSMNHigh20[iemin] = 0;
    fhSumEnCellsPerSMAcceptEvent     [iemin] = 0;                     
    fhNCellsPerSMAcceptEvent         [iemin] = 0;
    fhAverSumEnCellsPerSMAcceptEvent [iemin] = 0;
    
    fhSumEnCellsAcceptEventStrip[iemin] = 0;                     
    fhNCellsAcceptEventStrip    [iemin] = 0;
    fhSumEnCellsPerSMAcceptEventStrip[iemin] = 0;                     
    fhNCellsPerSMAcceptEventStrip[iemin] = 0;
    
    fhSumEnCellsAcceptEventBoth[iemin] = 0;                     
    fhNCellsAcceptEventBoth    [iemin] = 0;
    fhSumEnCellsPerSMAcceptEventBoth[iemin] = 0;                     
    fhNCellsPerSMAcceptEventBoth[iemin] = 0;
    
    fhSumEnCellsPerStrip        [iemin] = 0;                     
    fhNCellsPerStrip            [iemin] = 0;
    fhSumEnCellsPerStripPerSM   [iemin] = 0;                     
    fhNCellsPerStripPerSM       [iemin] = 0;
    
    fhSumEnCellsPerStripNHigh20     [iemin] = 0;                     
    fhNCellsPerStripNHigh20         [iemin] = 0;
    fhSumEnCellsPerStripPerSMNHigh20[iemin] = 0;                     
    fhNCellsPerStripPerSMNHigh20    [iemin] = 0;
    
    fhSumEnCellsPerStripAcceptEventStrip     [iemin] = 0;                     
    fhNCellsPerStripAcceptEventStrip         [iemin] = 0;
    fhSumEnCellsPerStripPerSMAcceptEventStrip[iemin] = 0;                     
    fhNCellsPerStripPerSMAcceptEventStrip    [iemin] = 0;
    
    fhSumEnCellsPerStripAcceptEventBoth     [iemin] = 0;                     
    fhNCellsPerStripAcceptEventBoth         [iemin] = 0;
    fhSumEnCellsPerStripPerSMAcceptEventBoth[iemin] = 0;                     
    fhNCellsPerStripPerSMAcceptEventBoth    [iemin] = 0;
    
    fhNSumEnCells                [iemin] = 0;                     
    fhNSumEnCellsAcceptEvent     [iemin] = 0;                     
    fhNSumEnCellsAcceptEventStrip[iemin] = 0;                     
    fhNSumEnCellsAcceptEventBoth[iemin] = 0;                     

    fhNSumEnCellsPerSM                [iemin] = 0;                     
    fhNSumEnCellsPerSMAcceptEvent     [iemin] = 0;                     
    fhNSumEnCellsPerSMAcceptEventStrip[iemin] = 0;                    
    fhNSumEnCellsPerSMAcceptEventBoth[iemin] = 0;                    
    
    fhNSumEnCellsPerStrip                [iemin] = 0;                     
    fhNSumEnCellsPerStripAcceptEvent     [iemin] = 0;                     
    fhNSumEnCellsPerStripAcceptEventStrip[iemin] = 0;                     
    fhNSumEnCellsPerStripAcceptEventBoth[iemin] = 0;                     
    
    fhNSumEnCellsPerStripPerSM                [iemin] = 0;                     
    fhNSumEnCellsPerStripPerSMAcceptEvent     [iemin] = 0;                     
    fhNSumEnCellsPerStripPerSMAcceptEventStrip[iemin] = 0;   
    fhNSumEnCellsPerStripPerSMAcceptEventBoth[iemin] = 0;   
    
    if ( iemin == 3 ) continue;
    
    fhFracNCells                [iemin] = 0;
    fhFracSumEnCells            [iemin] = 0;
    fhFracNCellsNHigh20         [iemin] = 0;
    fhFracSumEnCellsNHigh20     [iemin] = 0;   
    fhFracNCellsAcceptEvent     [iemin] = 0;
    fhFracSumEnCellsAcceptEvent [iemin] = 0;  
    fhFracNCellsPerSM           [iemin] = 0;
    fhFracSumEnCellsPerSM       [iemin] = 0;
    fhFracNCellsPerSMNHigh20    [iemin] = 0;
    fhFracSumEnCellsPerSMNHigh20[iemin] = 0;  
    fhFracNCellsPerSMAcceptEvent     [iemin] = 0;
    fhFracSumEnCellsPerSMAcceptEvent [iemin] = 0;  
  }
  
  for(Int_t i = 0; i < 4; i++)
  { 
    for(Int_t j = 0; j < 4; j++)
    {
      fhCellEnDiffColRowDiff      [i][j] = 0;
      fhCellEnSameColRowDiff      [i][j] = 0;
      fhCellEnDiffColRowDiffExoCut[i][j] = 0;
      fhCellEnSameColRowDiffExoCut[i][j] = 0;
      
      fhTimeDiffClusSameCellColRowDiff      [i][j] = 0; 
      fhTimeDiffClusDiffCellColRowDiff      [i][j] = 0;
      fhTimeDiffClusSameCellColRowDiffExoCut[i][j] = 0;
      fhTimeDiffClusDiffCellColRowDiffExoCut[i][j] = 0; 
    }
  }

  for (Int_t ieta = 0; ieta < 24; ieta++)
  {
    for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
    {
      for (Int_t ism = 0; ism < 20; ism++)
      {   
        fEnCellsStrip[icut][ism][ieta] = 0.; 
        fnCellsStrip [icut][ism][ieta] = 0 ; 
      }
    }
  }
  
  for (Int_t icase = 0; icase < 4; icase++)
  {
    fhEventCutClusterEnergyTime    [icase] = 0;
    fhEventCutClusterEnergy        [icase] = 0;
    fhEventCutClusterEnergyECellMax[icase] = 0;
    fhEventCutClusterEtaPhiGrid    [icase] = 0; 
    fhEventCutBunchCrossing        [icase] = 0;
    fhEventCutVertexX              [icase] = 0;
    fhEventCutVertexY              [icase] = 0;
    fhEventCutVertexZ              [icase] = 0;
    fhEventCutTrackPt              [icase] = 0;
    fhEventCutTrackEtaPhi          [icase] = 0;
  }
}

//
//____________________________________________________________
/// Fill histograms related to cells only.
/// \param cells: cells info list container
//____________________________________________________________
void AliAnaCaloExotics::CellHistograms(AliVCaloCells *cells)
{  
  Int_t ncells = cells->GetNumberOfCells();
  if( ncells <= 0 ) return;
  
  AliDebug(1,Form("%s cell entries %d", GetCalorimeterString().Data(), ncells));

  Int_t    icol   = -1, icolAbs = -1;
  Int_t    irow   = -1, irowAbs = -1;
  Int_t    iRCU   = -1;
  Float_t  amp    = 0.;
  Double_t time   = 0.;
  Int_t    id     = -1;
  //Bool_t   highG  = kFALSE;
  Float_t  exoticity = -1000;
  Float_t  exoticity50ns = -1000;
  
  Int_t    bc     = GetReader()->GetInputEvent()->GetBunchCrossNumber();
    
  Int_t   nCells[fgkNCellEnMinBins] ;
  Float_t eCells[fgkNCellEnMinBins] ;
  
  Int_t   nCellsPerSM[fgkNCellEnMinBins][20] ;
  Float_t eCellsPerSM[fgkNCellEnMinBins][20] ;  
 
  for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
  {
    nCells[icut] = 0 ;
    eCells[icut] = 0 ;
    for(Int_t ism = 0; ism < 20; ism++)
    {
      nCellsPerSM[icut][ism] = 0 ;
      eCellsPerSM[icut][ism] = 0 ;
    }
  }
  
  for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++)
  {
    if ( cells->GetCellNumber(iCell) < 0 || cells->GetAmplitude(iCell) < fCellAmpMin ) continue; 
    
    AliDebug(2,Form("Cell : amp %f, absId %d", cells->GetAmplitude(iCell), cells->GetCellNumber(iCell)));
    
    Int_t nModule = GetModuleNumberCellIndexesAbsCaloMap(cells->GetCellNumber(iCell),GetCalorimeter(), 
                                                         icol   , irow, iRCU,
                                                         icolAbs, irowAbs    );
    
    AliDebug(2,Form("\t module %d, column %d (%d), row %d (%d)", nModule,icolAbs,icol,irowAbs,irow));
    
    amp     = cells->GetAmplitude(iCell);
    
    time    = cells->GetTime(iCell);
    time *= 1.0e9;
    time-=fConstantTimeShift;
    
    id      = cells->GetCellNumber(iCell);
    
    //highG   = cells->GetCellHighGain(id);
    
    exoticity =  1-GetCaloUtils()->GetECross(id,cells,bc)/amp;
    
    GetReader()->GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(50); 
    exoticity50ns =  1-GetCaloUtils()->GetECross(id,cells,bc)/amp;
    GetReader()->GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(1e6); // Back to default 

    // Fill histograms
    
    fhCellExoAmp    ->Fill(amp       , exoticity, GetEventWeight());
    fhCellExoAmpTime->Fill(amp , time, exoticity, GetEventWeight());
    
    if ( fFillExo50ns )    
      fhCellExo50nsAmp->Fill(amp, exoticity50ns, GetEventWeight());

    if ( amp > fEMinForExo )
      fhCellExoGrid->Fill(icolAbs, irowAbs, exoticity, GetEventWeight());

    // If there is at least one cluster with nCellW>20 or 12, check the cell distribution
    if ( fLED20 )
    {
      fhCellGridTimeHighNCell20->Fill(icolAbs, irowAbs, fLED20Time, GetEventWeight());
    }
  
    if ( fLED12 )
    {
      fhCellGridTimeHighNCell12->Fill(icolAbs, irowAbs, fLED12Time, GetEventWeight());
    }
    
    for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
    {
      if ( amp > fCellEnMins[icut] && amp < fCellEnMax )
      {
        fhCellGridTime[icut]->Fill(icolAbs, irowAbs, time, GetEventWeight());
        
        nCells[icut]++;
        eCells[icut]+=amp;
        nCellsPerSM[icut][nModule]++;
        eCellsPerSM[icut][nModule]+=amp;  
      }
    }
  } // Cell loop

//  if ( nCellsPerSM[0][3] < 2 || eCellsPerSM[0][3] < 1 )
//    printf("SM3 nCells %d, sum E cells %2.2f \n",nCellsPerSM[0][3],eCellsPerSM[0][3] );
  
  // LED event rejection trial
  // use the total N cells or E sum for E cell > 0.5
  // Low activity in SM3 and very large activiy on any of the other SM
  //
  fAcceptEvent = kTRUE;
  if ( nCellsPerSM[0][3] <= fLowNCellsCutSM3 || eCellsPerSM[0][3] <= fLowEnergyCutSM3 )
  {
    for(Int_t ism = 0; ism < 20; ism++)
    {
      if ( ism == 3 ) continue;
      
      if ( nCellsPerSM[0][ism] >= fHighNCellsCutSM ) fAcceptEvent = kFALSE;
      if ( eCellsPerSM[0][ism] >= fHighEnergyCutSM ) fAcceptEvent = kFALSE;
    }
   
    if ( !fAcceptEvent )
    {
      printf("Reject event: ");
      for(Int_t jsm = 0; jsm < 20; jsm++)
      {
        if ( nCellsPerSM[0][jsm] > 0 ) 
          printf("\t SM%d: ncells %d; sum E %3.1f \n",jsm,nCellsPerSM[0][jsm],eCellsPerSM[0][jsm]);
        
      }
    }
  }
  else
  {
    Bool_t suspiciousEvent = kFALSE;

    for(Int_t ism = 0; ism < 20; ism++)
    {
      if ( ism == 3 ) continue;
      
      if ( nCellsPerSM[0][ism] >= fHighNCellsCutSM ) suspiciousEvent = kTRUE;
      if ( eCellsPerSM[0][ism] >= fHighEnergyCutSM ) suspiciousEvent = kTRUE;
    }
  
    if ( suspiciousEvent ) 
    {
      fhSM3NCellsSumEnSuspiciousEvents->Fill(nCellsPerSM[0][3], eCellsPerSM[0][3], GetEventWeight());  
      if ( fFillStripHisto ) 
        fhNStripsPerEventSuspicious->Fill(fEventNStripActive, GetEventWeight());
      for(Int_t ism = 0; ism < 20; ism++)
      {
        fhSumEnCellsPerSMEventSuspicious->Fill(eCellsPerSM[0][ism], ism, GetEventWeight());
        fhNCellsPerSMEventSuspicious    ->Fill(nCellsPerSM[0][ism], ism, GetEventWeight());
        if ( fFillStripHisto ) 
          fhNStripsPerEventSuspiciousPerSM->Fill(fEventNStripActiveSM[ism], ism, GetEventWeight());
      }
    }
  }
  
  if ( fAcceptEvent && fFillStripHisto )
  {
    fhNStripsPerEventAccept->Fill(fEventNStripActive, GetEventWeight());
    for (Int_t ism = 0; ism < 20; ism++)
      fhNStripsPerEventAcceptPerSM->Fill(fEventNStripActiveSM[ism], ism, GetEventWeight());
    
   
    for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
    { 
      for (Int_t ieta = 0; ieta < 24; ieta++)
      {
        for (Int_t ism = 0; ism < 20; ism++)
        {
          fhSumEnCellsPerStripEventAccept[icut]->Fill(fEnCellsStrip[icut][ism][ieta], GetEventWeight());
          fhNCellsPerStripEventAccept    [icut]->Fill(fnCellsStrip [icut][ism][ieta], GetEventWeight());

          fhSumEnCellsPerStripPerSMEventAccept[icut]->Fill(fEnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
          fhNCellsPerStripPerSMEventAccept    [icut]->Fill(fnCellsStrip [icut][ism][ieta], ism, GetEventWeight());
          
          fhNSumEnCellsPerStripAcceptEvent     [icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], GetEventWeight());
          fhNSumEnCellsPerStripPerSMAcceptEvent[icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
          
          if ( fEventNStripActive <= fEventMaxNumberOfStrips )
          {
            fhSumEnCellsPerStripAcceptEventBoth     [icut]->Fill(fEnCellsStrip[icut][ism][ieta], GetEventWeight());
            fhNCellsPerStripAcceptEventBoth         [icut]->Fill(fnCellsStrip [icut][ism][ieta], GetEventWeight());
            
            fhSumEnCellsPerStripPerSMAcceptEventBoth[icut]->Fill(fEnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
            fhNCellsPerStripPerSMAcceptEventBoth    [icut]->Fill(fnCellsStrip [icut][ism][ieta], ism, GetEventWeight());
            
            fhNSumEnCellsPerStripAcceptEventBoth     [icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], GetEventWeight());
            fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
          }
        } // ism
      } // ieta
    } // cut
  }
  
  // LED Event rejection
  //
  
  if ( fFillAllCellEventParamHisto  < 1 ) return;
  
  for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
  {
    Float_t averECells = 0;
    if ( nCells[icut] > 0 ) averECells = eCells[icut] / nCells[icut];
    
    fhSumEnCells [icut]->Fill(eCells[icut], GetEventWeight());
    fhNCells     [icut]->Fill(nCells[icut], GetEventWeight());
    fhNSumEnCells[icut]->Fill(eCells[icut],nCells[icut], GetEventWeight());
    
    if ( fFillAllCellEventParamHisto  > 1 )
      fhAverSumEnCells[icut]->Fill(averECells  , GetEventWeight());

    if ( fLED20 ) 
    {
      fhSumEnCellsNHigh20[icut]->Fill(eCells[icut], GetEventWeight());
      fhNCellsNHigh20    [icut]->Fill(nCells[icut], GetEventWeight());
      if ( fFillAllCellEventParamHisto  > 1 )
        fhAverSumEnCellsNHigh20[icut]->Fill(averECells  , GetEventWeight());
    }
    
    if ( fAcceptEvent )
    {
      fhSumEnCellsAcceptEvent [icut]->Fill(eCells[icut], GetEventWeight());
      fhNCellsAcceptEvent     [icut]->Fill(nCells[icut], GetEventWeight());
      fhNSumEnCellsAcceptEvent[icut]->Fill(eCells[icut],nCells[icut], GetEventWeight());

      if ( fFillAllCellEventParamHisto  > 1 )
        fhAverSumEnCellsAcceptEvent[icut]->Fill(averECells, GetEventWeight());
    }
 
    if ( fEventNStripActive <= fEventMaxNumberOfStrips )
    {
      fhSumEnCellsAcceptEventStrip [icut]->Fill(eCells[icut], GetEventWeight());
      fhNCellsAcceptEventStrip     [icut]->Fill(nCells[icut], GetEventWeight());
      fhNSumEnCellsAcceptEventStrip[icut]->Fill(eCells[icut],nCells[icut], GetEventWeight());
      if ( fAcceptEvent )
      {
        fhSumEnCellsAcceptEventBoth [icut]->Fill(eCells[icut], GetEventWeight());
        fhNCellsAcceptEventBoth     [icut]->Fill(nCells[icut], GetEventWeight());
        fhNSumEnCellsAcceptEventBoth[icut]->Fill(eCells[icut],nCells[icut], GetEventWeight());
      }
    }
    
    for(Int_t ism = 0; ism < 20; ism++)
    {
      averECells = 0;
      if ( nCellsPerSM[icut][ism] > 0 ) 
        averECells = eCellsPerSM[icut][ism] / nCellsPerSM[icut][ism];
      
      fhSumEnCellsPerSM [icut]->Fill(eCellsPerSM[icut][ism], ism, GetEventWeight());
      fhNCellsPerSM     [icut]->Fill(nCellsPerSM[icut][ism], ism, GetEventWeight());
      fhNSumEnCellsPerSM[icut]->Fill(eCellsPerSM[icut][ism],nCellsPerSM[icut][ism], ism, GetEventWeight());

      if ( fFillAllCellEventParamHisto  > 1 )
        fhAverSumEnCellsPerSM[icut]->Fill(averECells, ism, GetEventWeight());
      
      if ( fLED20 ) 
      {
        fhSumEnCellsPerSMNHigh20[icut]->Fill(eCellsPerSM[icut][ism], ism, GetEventWeight());
        fhNCellsPerSMNHigh20    [icut]->Fill(nCellsPerSM[icut][ism], ism, GetEventWeight());
        if ( fFillAllCellEventParamHisto  > 1 )
          fhAverSumEnCellsPerSMNHigh20[icut]->Fill(averECells, ism, GetEventWeight());
      }
      
      if ( fAcceptEvent )
      {
        fhSumEnCellsPerSMAcceptEvent [icut]->Fill(eCellsPerSM[icut][ism], ism, GetEventWeight());
        fhNCellsPerSMAcceptEvent     [icut]->Fill(nCellsPerSM[icut][ism], ism, GetEventWeight());
        fhNSumEnCellsPerSMAcceptEvent[icut]->Fill(eCellsPerSM[icut][ism],nCellsPerSM[icut][ism], ism, GetEventWeight());

        if ( fFillAllCellEventParamHisto  > 1 )
          fhAverSumEnCellsPerSMAcceptEvent[icut]->Fill(averECells, ism, GetEventWeight());
      }
      
      if ( fEventNStripActive <= fEventMaxNumberOfStrips )
      {
        fhSumEnCellsPerSMAcceptEventStrip [icut]->Fill(eCellsPerSM[icut][ism], ism, GetEventWeight());
        fhNCellsPerSMAcceptEventStrip     [icut]->Fill(nCellsPerSM[icut][ism], ism, GetEventWeight());
        fhNSumEnCellsPerSMAcceptEventStrip[icut]->Fill(eCellsPerSM[icut][ism],nCellsPerSM[icut][ism], ism, GetEventWeight());
        
        if ( fAcceptEvent )
        {
          fhSumEnCellsPerSMAcceptEventBoth [icut]->Fill(eCellsPerSM[icut][ism], ism, GetEventWeight());
          fhNCellsPerSMAcceptEventBoth     [icut]->Fill(nCellsPerSM[icut][ism], ism, GetEventWeight());
          fhNSumEnCellsPerSMAcceptEventBoth[icut]->Fill(eCellsPerSM[icut][ism],nCellsPerSM[icut][ism], ism, GetEventWeight());
        }
      }
    } // Per SM
    
    if ( fFillAllCellEventParamHisto < 2 ) continue; 
    
    if ( icut  == 0 ) continue; 
    
    Float_t frNCells = 0;
    if ( nCells[0] > 0 ) frNCells = (1.*nCells[icut]) / nCells[0];
    
    Float_t frEnCells = 0;
    if ( eCells[0] > 0 ) frEnCells = eCells[icut] / eCells[0];

    fhFracNCells    [icut-1]->Fill(frNCells , GetEventWeight());
    fhFracSumEnCells[icut-1]->Fill(frEnCells, GetEventWeight());
    
    if ( fLED20 ) 
    {
      fhFracNCellsNHigh20    [icut-1]->Fill(frNCells , GetEventWeight());
      fhFracSumEnCellsNHigh20[icut-1]->Fill(frEnCells, GetEventWeight());
    }
    
    if ( fAcceptEvent )
    {
      fhFracNCellsAcceptEvent    [icut-1]->Fill(frNCells , GetEventWeight());
      fhFracSumEnCellsAcceptEvent[icut-1]->Fill(frEnCells, GetEventWeight());
    }
    
    for(Int_t ism = 0; ism < 20; ism++)
    {
      frNCells = 0;
      if ( nCellsPerSM[0][ism] > 0 ) 
        frNCells = (1.*nCellsPerSM[icut][ism]) / nCellsPerSM[0][ism];
      
      Float_t frEnCells = 0;
      if ( eCellsPerSM[0][ism] > 0 ) 
        frEnCells = eCellsPerSM[icut][ism] / eCellsPerSM[0][ism];
      
      fhFracNCellsPerSM    [icut-1]->Fill(frNCells , ism, GetEventWeight());
      fhFracSumEnCellsPerSM[icut-1]->Fill(frEnCells, ism, GetEventWeight());
      
      if ( fLED20 ) 
      {
        fhFracNCellsPerSMNHigh20    [icut-1]->Fill(frNCells , ism, GetEventWeight());
        fhFracSumEnCellsPerSMNHigh20[icut-1]->Fill(frEnCells, ism, GetEventWeight());
      }
      
      if ( fAcceptEvent )
      {
        fhFracNCellsPerSMAcceptEvent    [icut-1]->Fill(frNCells , ism, GetEventWeight());
        fhFracSumEnCellsPerSMAcceptEvent[icut-1]->Fill(frEnCells, ism, GetEventWeight());
      }
    } // Per SM
  }
}

//
//____________________________________________________________
/// Fill histograms related to strips (2x24, eta,phi) only.
/// \param cells: cells info list container
//____________________________________________________________
void AliAnaCaloExotics::StripHistograms(AliVCaloCells *cells)
{  
  Int_t ncells = cells->GetNumberOfCells();
  if( ncells <= 0 ) return;
  
  fEventNStripActive = 0;

  Float_t amp1   = 0., amp2   = 0. ;
  Int_t   absId1 = -1, absId2 = -1 ;
  
  for (Int_t ism = 0; ism < 20; ism++)
  {
    fEventNStripActiveSM[ism] = 0;
  
    for (Int_t ieta = 0; ieta < 48; ieta=ieta+2)
    {
      for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
      {
        fEnCellsStrip[icut][ism][ieta/2] = 0.; 
        fnCellsStrip [icut][ism][ieta/2] = 0 ; 
      }
      
      for (Int_t iphi = 0; iphi < 24; iphi++)
      {
        absId1 = GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(ism, iphi, ieta);
        if ( absId1 < 0 || absId1 > 17664 ) continue;
        
        absId2 = GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(ism, iphi, ieta+1);
        if ( absId2 < 0 || absId2 > 17664 ) continue;   
        
        amp1 = cells->GetCellAmplitude(absId1);
        amp2 = cells->GetCellAmplitude(absId2);
//        printf("ieta %d iphi %d: a) id %d amp %2.2f; b) id %d amp %2.2f\n",
//               ieta, iphi,absId1,amp1,absId2,amp2);

        for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
        {
          if ( amp1 > fCellEnMins[icut] && amp1 < fCellEnMax )
          {            
            fnCellsStrip [icut][ism][ieta/2]++;
            fEnCellsStrip[icut][ism][ieta/2]+=amp1;
//            if( icut == 0 && (fnCellsStrip [icut][ism][ieta/2] > 0 || fEnCellsStrip[icut][ism][ieta/2] > 0))
//              printf("Task: cell1 %d amp %f; SM %d, strip %d, n cells %d, sum E %f \n",
//                     absId1, amp1,
//                     ism,ieta,
//                     fnCellsStrip [0][ism][ieta/2],
//                     fEnCellsStrip[0][ism][ieta/2] 
//                     );
          }
          
          if ( amp2 > fCellEnMins[icut] && amp2 < fCellEnMax )
          {            
            fnCellsStrip [icut][ism][ieta/2]++;
            fEnCellsStrip[icut][ism][ieta/2]+=amp2;
//            if( icut == 0 && (fnCellsStrip [icut][ism][ieta/2] > 0 || fEnCellsStrip[icut][ism][ieta/2] > 0))
//              printf("Task: cell2 %d amp %f; SM %d, strip %d, n cells %d , sum E %f \n",
//                     absId2, amp2,
//                     ism,ieta,
//                     fnCellsStrip [0][ism][ieta/2],
//                     fEnCellsStrip[0][ism][ieta/2] 
//                     );
          }
        } //  e cut
        
      }// iphi
      
      // Fill histograms per Strip
      for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
      {   
        fhSumEnCellsPerStrip     [icut]->Fill(fEnCellsStrip[icut][ism][ieta/2], GetEventWeight());
        fhNCellsPerStrip         [icut]->Fill(fnCellsStrip [icut][ism][ieta/2], GetEventWeight());
        
        fhSumEnCellsPerStripPerSM[icut]->Fill(fEnCellsStrip[icut][ism][ieta/2], ism, GetEventWeight());
        fhNCellsPerStripPerSM    [icut]->Fill(fnCellsStrip [icut][ism][ieta/2], ism, GetEventWeight());
       
        fhNSumEnCellsPerStrip     [icut]->Fill(fEnCellsStrip[icut][ism][ieta/2],fnCellsStrip[icut][ism][ieta/2], GetEventWeight());
        fhNSumEnCellsPerStripPerSM[icut]->Fill(fEnCellsStrip[icut][ism][ieta/2],fnCellsStrip[icut][ism][ieta/2], ism, GetEventWeight());
      
        if ( fLED20 ) 
        {
          fhSumEnCellsPerStripNHigh20     [icut]->Fill(fEnCellsStrip[icut][ism][ieta/2], GetEventWeight());
          fhNCellsPerStripNHigh20         [icut]->Fill(fnCellsStrip [icut][ism][ieta/2], GetEventWeight());
        
          fhSumEnCellsPerStripPerSMNHigh20[icut]->Fill(fEnCellsStrip[icut][ism][ieta/2], ism, GetEventWeight());
          fhNCellsPerStripPerSMNHigh20    [icut]->Fill(fnCellsStrip [icut][ism][ieta/2], ism, GetEventWeight());
        }
      }
    } // ieta
  }// sm    
  
  // Count per event over event cut
  // Low activity on SM3 for emin = 0.5
  Bool_t bSM3StripsLowActivity = kTRUE;
  for (Int_t ieta = 0; ieta < 24; ieta++)
  {
    if ( fEnCellsStrip[0][3][ieta] > fLowEnergyCutSM3Strip || 
         fnCellsStrip [0][3][ieta] > fLowNCellsCutSM3Strip   ) 
      bSM3StripsLowActivity = kFALSE;
  }
  
  if ( bSM3StripsLowActivity )
  {
    
    Int_t   maxNCells = 20;
    Float_t maxECells = 45;
    for (Int_t ism = 0; ism < 20; ism++)
    {
      if ( ism == 3 ) continue ;
     
      maxNCells = fHighNCellsCutStrip[0];
      maxECells = fHighEnergyCutStrip[0];
      if ( ism == 10 || ism == 11 || 
           ism == 18 || ism == 19   ) 
      {
        maxNCells = fHighNCellsCutStrip[1];
        maxECells = fHighEnergyCutStrip[1];
      }
      
      for (Int_t ieta = 0; ieta < 24; ieta++)
      {
//        if(fnCellsStrip[0][ism][ieta] > 0 || fEnCellsStrip[0][ism][ieta] > 0 )
//          printf("Task: SM %d, strip %d, n cells %d < %d , sum E %f < %f \n",
//               ism,ieta,
//               fnCellsStrip [0][ism][ieta],maxNCells,
//               fEnCellsStrip[0][ism][ieta],maxECells 
//               );
        if( fEnCellsStrip[0][ism][ieta] >= maxECells || 
            fnCellsStrip [0][ism][ieta] >= maxNCells   )
        {
          fEventNStripActive++;
          fEventNStripActiveSM[ism]++;          
        }
      } // ieta
    } // ism
  } // bSM03
  
  fhNStripsPerEvent->Fill(fEventNStripActive, GetEventWeight());
  for (Int_t ism = 0; ism < 20; ism++)
    fhNStripsPerEventPerSM->Fill(fEventNStripActiveSM[ism], ism, GetEventWeight());
  
  //printf("Task: Number of strips %d\n",fEventNStripActive);
  
  // Strip activity if cut on number of active strips
  if ( fEventNStripActive <= fEventMaxNumberOfStrips )
  {
    for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
    { 
      for (Int_t ieta = 0; ieta < 24; ieta++)
      {
        for (Int_t ism = 0; ism < 20; ism++)
        {
          fhSumEnCellsPerStripAcceptEventStrip     [icut]->Fill(fEnCellsStrip[icut][ism][ieta], GetEventWeight());
          fhNCellsPerStripAcceptEventStrip         [icut]->Fill(fnCellsStrip [icut][ism][ieta], GetEventWeight());
          
          fhSumEnCellsPerStripPerSMAcceptEventStrip[icut]->Fill(fEnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
          fhNCellsPerStripPerSMAcceptEventStrip    [icut]->Fill(fnCellsStrip [icut][ism][ieta], ism, GetEventWeight());
          
          fhNSumEnCellsPerStripAcceptEventStrip     [icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], GetEventWeight());
          fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut]->Fill(fEnCellsStrip[icut][ism][ieta],fnCellsStrip[icut][ism][ieta], ism, GetEventWeight());
        } // icut
      } // ieta
    } // sm
    
  } // strip event cut
  
}

//____________________________________________________________________________
/// Fill clusters related histograms, execute here the loop of clusters
/// apply basic selection cuts (track matching, goondess, exoticity, timing)
/// and the call to the different methods
/// filling different type of histograms:
/// * Basic cluster histograms for good or bad clusters
/// * Exotic cluster histograms
/// * Cells in cluster
/// * Invariant mass
/// * Matched clusters histograms
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCaloExotics::ClusterHistograms(const TObjArray *caloClusters,
                                                AliVCaloCells* cells)
{
  fLED20 = kFALSE;
  fLED12 = kFALSE;   
  fLED20Time = -1e6;
  fLED12Time = -1e6;
  
  Int_t  nCaloClusters = caloClusters->GetEntriesFast() ;
  Int_t  bc            = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  // Get vertex, not needed.
  Double_t v[3] = {0,0,0};
  //GetReader()->GetVertex(v);
  Int_t nNCellW20    = 0;     
  Int_t nNCellW12    = 0;        
  Int_t nExotic      = 0;
  Int_t nExotic1Cell = 0;
  Int_t nExoticNCell = 0;
  
  Float_t highestAmpMaxExo = 0;
  
  AliDebug(1,Form("In %s there are %d clusters", GetCalorimeterString().Data(), nCaloClusters));
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++)
  {
    AliDebug(1,Form("Cluster: %d/%d, data %d",iclus+1,nCaloClusters,GetReader()->GetDataType()));
    
    AliVCluster* clus =  (AliVCluster*) caloClusters->At(iclus);
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    Float_t en  = fClusterMomentum.E();
    Float_t eta = fClusterMomentum.Eta();
    Float_t phi = GetPhi(fClusterMomentum.Phi());
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if ( IsFiducialCutOn() ) 
      in =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum.Eta(),fClusterMomentum.Phi(),GetCalorimeter()) ;
    
    if ( !in )
    {
      AliDebug(1,Form("Remove cluster with phi %2.2f and eta %2.2f",phi*TMath::RadToDeg(),eta));
      continue;
    }
    
    Int_t    nCaloCellsPerCluster = clus->GetNCells();
    //Int_t  nLabel = clus->GetNLabels();
    //Int_t *labels = clus->GetLabels();
    
    // Cluster mathed with track?
    //
    Bool_t     matched  = 0;
    AliVTrack *track    = 0;
    Double_t   eOverP   = 0;
    Double_t   tmom     = 0;
    Bool_t     positive = kFALSE;
    if ( fFillMatchingHisto )
    {
      matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(), GetReader()->GetInputEvent());
      track   = GetCaloUtils()->GetMatchedTrack(clus, GetReader()->GetInputEvent());
      
      if ( track )
      {
        tmom   = track->P();
        eOverP = en/tmom;
        
        positive = (track->Charge()>0);
      }
    }
    
    // Residuals
    Float_t deta  = clus->GetTrackDz();
    Float_t dphi  = clus->GetTrackDx();  
    
    // Get the 1-fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    Int_t icolMax  = -1, icolMaxAbs = -1;
    Int_t irowMax  = -1, irowMaxAbs = -1;
    Int_t iRCU     = -1;
    Int_t nSM  = GetModuleNumberCellIndexesAbsCaloMap(absIdMax,GetCalorimeter(), 
                                                      icolMax   , irowMax, iRCU,
                                                      icolMaxAbs, irowMaxAbs    );
        
    // Get time of max cell/cluster
    Double_t tmax  = cells->GetCellTime(absIdMax);
    tmax*=1.e9;
    tmax-=fConstantTimeShift;
    
    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    Float_t enOrg  = ampMax/(1-maxCellFraction); // cluster energy without non linearity correction

    Float_t exoticity  = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc)/ampMax;

    Float_t exoticity50ns = -1;
    if ( fFillExo50ns )
    {
      GetReader()->GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(50); 
      exoticity50ns =  1-GetCaloUtils()->GetECross(absIdMax,cells,bc)/ampMax;
      GetReader()->GetCaloUtils()->GetEMCALRecoUtils()->SetExoticCellDiffTimeCut(1e6); // Back to default 
    }
    
    Float_t m02  = clus->GetM02();
    Float_t m20  = clus->GetM20();
    
    Int_t ebin = -1;
    for(Int_t i = 0; i < fgkNEBins-1; i++) 
    {
      if(en >= fEnergyBins[i] && en < fEnergyBins[i+1] ) ebin = i;
    }
    
    // Cell cluster loop
    //
    Int_t   nCellSame   = 0, nCellSame5  = 0; 
    Int_t   nCellDiff   = 0, nCellDiff5  = 0;
    Int_t   nCellW      = 0, nCellSameW  = 0,  nCellDiffW = 0;
    Int_t   nCellTimeDiff = 0, nCellSameTimeDiff = 0, nCellDiffTimeDiff = 0;
    Int_t   rowDiff     = -100, colDiff = -100;
    Int_t   rowDiffAbs  = -100, colDiffAbs = -100;
    Float_t enSame      = 0,    enDiff  = 0;
    Float_t enSame5     = 0,    enDiff5 = 0;
    Float_t enSameW     = 0,    enDiffW = 0;
    Float_t enSameTimeDiff = 0, enDiffTimeDiff = 0;

    const Int_t nMinEnCut = 18;
    Int_t  nCellSameMinEn[nMinEnCut],  nCellDiffMinEn[nMinEnCut];
    for(Int_t imin = 0; imin < 20; imin++)
    {
      nCellSameMinEn[imin] = 0;
      nCellDiffMinEn[imin] = 0;
    }
    
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId     = clus->GetCellsAbsId()[ipos];  
      
      Float_t amp       = cells->GetCellAmplitude(absId);
      
      if ( absId == absIdMax || amp < 0.1 ) continue;
      
      Bool_t  sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
      colDiffAbs = TMath::Abs(colDiff);
      rowDiffAbs = TMath::Abs(rowDiff);
      
      if ( sameTCard ) 
      { 
        nCellSame ++; 
        enSame += amp; 
        if ( rowDiffAbs <= 1 && colDiffAbs <= 1 ) 
        {
          enSame5 += amp;
          nCellSame5++;
        }        
      }
      else             
      { 
        nCellDiff ++; 
        enDiff += amp; 
        if ( rowDiffAbs <= 1 && colDiffAbs <= 1 ) 
        {
          enDiff5 += amp;
          nCellDiff5++;
        }
      }
      
      for(Int_t imin = 0; imin < nMinEnCut; imin++)
      {
        Float_t enCellMin  = 0.2 + imin*0.1;
        if ( amp > enCellMin )
        {
          if ( sameTCard ) nCellSameMinEn[imin]++;
          else             nCellDiffMinEn[imin]++;
        }
      }
      
      Float_t weight    = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(amp, en);
      
      if( weight > 0.01 ) 
      {
        nCellW++;
        if ( sameTCard ) { nCellSameW++; enSameW+=amp; }
        else             { nCellDiffW++; enDiffW+=amp; }
      }
      
      Double_t time  = cells->GetCellTime(absId);
      time*=1.e9;
      time-=fConstantTimeShift;
      
      Float_t tdiff = tmax-time;
      if( tdiff < 50 ) 
      {
        nCellTimeDiff++;
        if ( sameTCard ) { nCellSameTimeDiff++; enSameTimeDiff+=amp; }
        else             { nCellDiffTimeDiff++; enDiffTimeDiff+=amp; }
      }
    } // cells in cluster loop
    
    Float_t fracEnDiffSame  = 0, fracNCellDiffSame  = 0, fracEnNCellDiffSame  = 0;
    Float_t fracEnDiffSameW = 0, fracNCellDiffSameW = 0, fracEnNCellDiffSameW = 0;
    Float_t fracEnDiffSame5 = 0, fracNCellDiffSame5 = 0, fracEnNCellDiffSame5 = 0;
    
    if ( fFillSameDiffFracHisto )
    {
      if ( enSame    > 0 ) fracEnDiffSame      = enDiff / enSame;
      if ( nCellSame > 0 ) fracNCellDiffSame   = nCellDiff*1. / nCellSame*1.;
      if ( nCellDiff > 0 && nCellSame > 0 && enSame > 0)
        fracEnNCellDiffSame = (enDiff/nCellDiff*1.) / (enSame/nCellSame*1.);
      
      if ( enSameW    > 0 ) fracEnDiffSameW      = enDiffW / enSameW;
      if ( nCellSameW > 0 ) fracNCellDiffSameW   = nCellDiffW*1. / nCellSameW*1.;
      if ( nCellDiffW > 0 && nCellSameW > 0 && enSameW > 0)
        fracEnNCellDiffSameW = (enDiffW/nCellDiffW) / (enSameW/nCellSameW);
      
      if ( enSame5    > 0 ) fracEnDiffSame5      = enDiff5 / enSame5;
      if ( nCellSame5 > 0 ) fracNCellDiffSame5   = nCellDiff5*1. / nCellSame5*1.;
      if ( nCellDiff5 > 0 && nCellSame5 > 0 && enSame5 > 0)
        fracEnNCellDiffSame5 = (enDiff5/nCellDiff5*1.) / (enSame5/nCellSame5*1.);
    }

    AliDebug(1,Form("cluster: E %2.3f, F+ %2.3f, eta %2.3f, phi %2.3f, col %d, row %d, ncells %d,"
                    "match %d; cell max: id %d, en %2.3f, time %2.3f, m02 %2.2f",
                    en,exoticity,eta,phi*TMath::RadToDeg(), icolMaxAbs, irowMaxAbs, nCaloCellsPerCluster,
                    matched, absIdMax,ampMax,tmax,m02));  

    // Counters
    if ( nCellW > 20 ) { nNCellW20++; fLED20 = kTRUE; fLED20Time = tmax; }  
    if ( nCellW > 12 ) { nNCellW12++; fLED12 = kTRUE; fLED12Time = tmax; }
    if ( en > 5 )
    {
      if ( exoticity > fExoCut ) nExotic++;
      if ( nCaloCellsPerCluster == 1 ) nExotic1Cell++;
      if ( exoticity > fExoCut && nCaloCellsPerCluster > 1 ) nExoticNCell++;
    }
    
    if ( exoticity > fExoCut && ampMax > highestAmpMaxExo ) 
    {
      //printf("Exotic en %f, amp %f\n",en,ampMax);
      highestAmpMaxExo = ampMax;
    }
  
    //
    // Fill histograms related to single cluster 
    //
    
    // Cluster/Cell max Time
    //
    fhTimeEnergyNCells ->Fill(en, tmax, nCaloCellsPerCluster, GetEventWeight());
    fhTimeEnergyNCellsW->Fill(en, tmax, nCellW              , GetEventWeight());
    
    if ( nCellW > fNCellHighCut )
      fhTimeNCellCut->Fill(tmax, GetEventWeight());

    if ( nCaloCellsPerCluster > 1 )
    {
      fhTimeEnergyExo  ->Fill(en, tmax, exoticity , GetEventWeight());
      fhTimeEnergyM02  ->Fill(en, tmax, m02       , GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhTimeEnergy1Cell->Fill(en, tmax,            GetEventWeight());

//    if(maxCellFraction < 0.2) printf("frac %2.2f, en cell %2.2f, en cluster %2.2f, refrac %2.2f\n",
//                                     maxCellFraction,ampMax,en,ampMax/en);
    
    if ( fFillOpenTimeHisto )
    {
      if ( nCaloCellsPerCluster > 1 )
      {
        //if(en > 10) printf("Amp %f, en %f, enOrg %f, frac %f\n",ampMax,en,enOrg,1-maxCellFraction);
        fhCellMaxClusterEnOpenTime           ->Fill(enOrg, ampMax                   , GetEventWeight());
        fhCellMaxClusterEnRatioOpenTime      ->Fill(en   , 1-maxCellFraction        , GetEventWeight());
        fhCellMaxClusterEnRatioNCellWOpenTime->Fill(en   , 1-maxCellFraction, nCellW, GetEventWeight());
      }
      
      fhNCellsPerClusterOpenTime             ->Fill(en    , nCaloCellsPerCluster       , GetEventWeight());
      fhNCellsPerClusterWOpenTime            ->Fill(en    , nCellW                     , GetEventWeight());
      fhNCellsPerClusterEMaxCellOpenTime     ->Fill(ampMax, nCaloCellsPerCluster       , GetEventWeight());
      fhNCellsPerClusterWEMaxCellOpenTime    ->Fill(ampMax, nCellW                     , GetEventWeight());
      
      for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
      {
        Int_t   absId     = clus->GetCellsAbsId()[ipos];  
        
        Float_t amp       = cells->GetCellAmplitude(absId);
        
        if ( absId == absIdMax || amp < 0.1 ) continue;
        
        fhCellEnNCellWOpenTime    ->Fill(en    , amp, nCellW, GetEventWeight());
        fhCellEnNCellWEMaxOpenTime->Fill(ampMax, amp, nCellW, GetEventWeight());
        
        Double_t time  = cells->GetCellTime(absId);
        time*=1.e9;
        time-=fConstantTimeShift;
        
        Float_t tdiff = tmax-time;
        
        fhCellTimeDiffNCellWOpenTime    ->Fill(en    , tdiff, nCellW, GetEventWeight());
        fhCellTimeDiffNCellWEMaxOpenTime->Fill(ampMax, tdiff, nCellW, GetEventWeight());
      }
    }
    
    //----------------------------------------    
    // Apply time cut 
    //----------------------------------------    
    //
    if ( tmax < fTimeCutMin || tmax> fTimeCutMax ) continue;
    
    // N cells per cluster
    //
    fhNCellsPerCluster             ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
    fhNCellsPerClusterW            ->Fill(en, nCellW                         , GetEventWeight());
    fhNCellsPerClusterTimeDiff     ->Fill(en, nCellTimeDiff                  , GetEventWeight());
    fhNCellsPerClusterExo          ->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    if ( fFillExo50ns )
      fhNCellsPerClusterExo50ns    ->Fill(en, nCaloCellsPerCluster, exoticity50ns, GetEventWeight());

    fhNCellsPerClusterM02          ->Fill(en, nCaloCellsPerCluster, m02      , GetEventWeight());
    
    fhNCellsPerClusterEMaxCell     ->Fill(ampMax, nCaloCellsPerCluster       , GetEventWeight());
    fhNCellsPerClusterWEMaxCell    ->Fill(ampMax, nCellW                     , GetEventWeight());
    
    if ( fFillPerSMHisto ) 
    {
      fhNCellsPerClusterPerSM        ->Fill(en, nCaloCellsPerCluster, nSM      , GetEventWeight());
      fhNCellsPerClusterWPerSM       ->Fill(en, nCellW              , nSM      , GetEventWeight());
      fhNCellsPerClusterExoPerSM[nSM]->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    }
    
    if ( matched && fFillMatchingHisto  && track )
    {
      fhNCellsPerClusterTrackMatch   ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
      fhNCellsPerClusterExoTrackMatch->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    }
    
    // Fill histograms for clusters with 1 cell
    //
    if ( nCaloCellsPerCluster == 1 )
    {
      if ( fFill1CellHisto )
      {
        fhExoticity1Cell->Fill(en, exoticity, GetEventWeight());
        if ( fFillExo50ns )
          fhExoticity50ns1Cell->Fill(en, exoticity50ns, GetEventWeight());
        
        fhEtaPhiGridEn1Cell->Fill(icolMaxAbs, irowMaxAbs, en, GetEventWeight());
        
        if ( matched && fFillMatchingHisto && track )
        {
          fhEOverP1Cell->Fill(en, eOverP, GetEventWeight());
          
          if ( positive )
          {
            fhTrackMatchedDEtaPos1Cell->Fill(en, deta, GetEventWeight());
            fhTrackMatchedDPhiPos1Cell->Fill(en, dphi, GetEventWeight());
            
            if ( en > fEMinForExo ) 
            {
              fhTrackMatchedDEtaDPhiPos1Cell->Fill(deta, dphi, GetEventWeight());
            }
          }
          else 
          {
            fhTrackMatchedDEtaNeg1Cell->Fill(en, deta, GetEventWeight());
            fhTrackMatchedDPhiNeg1Cell->Fill(en, dphi, GetEventWeight());
            
            if(en > fEMinForExo) 
            {
              fhTrackMatchedDEtaDPhiNeg1Cell->Fill(deta, dphi, GetEventWeight());
            }
          }
        } // matched 1 cell histograms
      } // 1 cell histograms
      
      continue;
    }

    //----------------------------------------    
    // Now only clusters with more than 1 cell
    //----------------------------------------
    
    // Exoticity
    //
    fhExoticityEClus   ->Fill(en    , exoticity , GetEventWeight());
    fhExoticityEMaxCell->Fill(ampMax, exoticity , GetEventWeight());
    
    if ( fFillPerSMHisto ) 
      fhExoticityEClusPerSM ->Fill(en, exoticity, nSM, GetEventWeight());
    
    if ( matched && fFillMatchingHisto && track )
      fhExoticityEClusTrackMatch->Fill(en, exoticity, GetEventWeight());      
    
    if ( fFillExo50ns )
      fhExoticity50nsEClus->Fill(en, exoticity50ns, GetEventWeight());
    
    // Calculate exoticity for different E min cuts
    //
    if ( fFillExoEnMinCut )
    {
      // Recalculate exoticity with different cell minimum energy thresholds
      for(Int_t imin = 0; imin < 41; imin++)
      {
        Float_t enCellMin  = 0.1 + imin*0.05;
        Float_t exoMod = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc,enCellMin)/ampMax;
        //        if(en > 5 && exoticity > 0.97 && imin == 0)
        //        {
        //          printf("en %f, imin %d, emin %f exo org %f; exoW %f, exo recalc %f\n",
        //                 en, imin, enCellMin,exoticity,exoticityW, exoMod);
        //        }
        
        fhExoticityECellMinCut->Fill(en, exoMod, enCellMin, GetEventWeight());
      }
      
      // Calculate exoticity with cut on w
      //
      Float_t exoticityW = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc,0,kTRUE,en)/ampMax;
      
      fhExoticityWEClus     ->Fill(en, exoticityW      , GetEventWeight());
      fhTimeEnergyExoW      ->Fill(en, tmax, exoticityW, GetEventWeight());
      fhM02EnergyExoW       ->Fill(en, m02 , exoticityW, GetEventWeight());
      fhNCellsPerClusterExoW->Fill(en, nCaloCellsPerCluster, exoticityW, GetEventWeight());
    }
    
    // Acceptance
    //
    if ( en > fEMinForExo )
    {
      fhEtaPhiGridExoEnCut  ->Fill(icolMaxAbs, irowMaxAbs, exoticity, GetEventWeight());
      fhEtaPhiGridNCellEnCut->Fill(icolMaxAbs, irowMaxAbs, nCaloCellsPerCluster, GetEventWeight());
    }
    
    if ( exoticity > fExoCut )
      fhEtaPhiGridEnExoCut  ->Fill(icolMaxAbs, irowMaxAbs, en       , GetEventWeight());
   
    if ( nCellW > fNCellHighCut )
      fhEtaPhiGridEnHighNCells->Fill(icolMaxAbs, irowMaxAbs, en       , GetEventWeight());

    // Cluster energy vs cell max energy
    //
    fhCellMaxClusterEn           ->Fill(enOrg, ampMax                      , GetEventWeight());
    fhCellMaxClusterEnRatio      ->Fill(en   , 1-maxCellFraction           , GetEventWeight());
    fhCellMaxClusterEnRatioNCellW->Fill(en   , 1-maxCellFraction, nCellW   , GetEventWeight());
    fhCellMaxClusterEnRatioExo   ->Fill(en   , 1-maxCellFraction, exoticity, GetEventWeight());
    
    // Cell cluster loop
    //
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId     = clus->GetCellsAbsId()[ipos];  
     
      Float_t amp       = cells->GetCellAmplitude(absId);

      if ( absId == absIdMax || amp < 0.1 ) continue;

      fhCellEnNCellW    ->Fill(en    , amp, nCellW, GetEventWeight());
      fhCellEnNCellWEMax->Fill(ampMax, amp, nCellW, GetEventWeight());
      
      Float_t weight    = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(amp, en);
      //
      //      if( weight > 0.01 ) 
      //      {
      //        if ( ebin >= 0 && ebin < fgkNEBins-1 )
      //        {
      //          fhClusterColRowExoW[icolMax%2][ebin]->Fill(colDiff, rowDiff, exoticity, GetEventWeight());
      //        }
      //      }
      
      Double_t time  = cells->GetCellTime(absId);
      time*=1.e9;
      time-=fConstantTimeShift;
      
      Float_t tdiff = tmax-time; 
      
      fhTimeDiffClusCellExo->Fill(en , tdiff, exoticity, GetEventWeight());
      fhTimeDiffClusCellM02->Fill(en , tdiff, m02      , GetEventWeight());  
      
      if ( weight > 0.01 )
        fhTimeDiffWClusCellExo->Fill(en , tdiff, exoticity, GetEventWeight());
      
      fhCellTimeDiffNCellW    ->Fill(en    , tdiff, nCellW, GetEventWeight());
      fhCellTimeDiffNCellWEMax->Fill(ampMax, tdiff, nCellW, GetEventWeight());
      
      if ( en > fEMinForExo )
        fhTimeDiffAmpClusCellExo->Fill(amp, tdiff, exoticity, GetEventWeight());
      
      Bool_t  sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
      colDiffAbs = TMath::Abs(colDiff);
      rowDiffAbs = TMath::Abs(rowDiff);
      
      if ( fFillClusterColRowHisto && ebin >= 0 && ebin < fgkNEBins-1 )
      {
        fhClusterColRowExo[icolMax%2][ebin]->Fill(colDiff, rowDiff, exoticity, GetEventWeight());
        //if(exoticity < fExoCut)
        //  fhClusterColRow   [icolMax%2][ebin]->Fill(colDiff, rowDiff, GetEventWeight());
        if ( fFillPerSMHisto && nCellW > fNCellHighCut )
          fhClusterColRowPerSMHighNCell[ebin]->Fill(colDiff, rowDiff, nSM, GetEventWeight());
      }
            
      if ( sameTCard ) 
      { 
        fhCellEnSameExo->Fill(en, amp, exoticity, GetEventWeight());
        fhTimeDiffClusCellSameTCardExo->Fill(en , tdiff, exoticity, GetEventWeight());

        if( colDiffAbs < 4 && rowDiffAbs < 4 ) 
        {
          if(colDiffAbs > 1) printf("Same colDiff %d rowDiff %d\n",colDiffAbs,rowDiffAbs);
         
          if ( exoticity > fExoCut )
          {
            fhCellEnSameColRowDiffExoCut[colDiffAbs][rowDiffAbs]->Fill(en, amp, GetEventWeight());
            fhTimeDiffClusSameCellColRowDiffExoCut[colDiffAbs][rowDiffAbs]->Fill(en, tdiff, GetEventWeight());
          }
          else
          {
            fhCellEnSameColRowDiff[colDiffAbs][rowDiffAbs]->Fill(en, amp, GetEventWeight());
            fhTimeDiffClusSameCellColRowDiff[colDiffAbs][rowDiffAbs]->Fill(en, tdiff, GetEventWeight());
          }
        }

      }
      else             
      { 
        fhCellEnDiffExo->Fill(en, amp, exoticity, GetEventWeight());
        fhTimeDiffClusCellDiffTCardExo->Fill(en , tdiff, exoticity, GetEventWeight());

        if( colDiffAbs < 4 && rowDiffAbs < 4 ) 
        {
          if ( exoticity > fExoCut )
          {
            fhCellEnDiffColRowDiffExoCut[colDiffAbs][rowDiffAbs]->Fill(en, amp, GetEventWeight());
            fhTimeDiffClusDiffCellColRowDiffExoCut[colDiffAbs][rowDiffAbs]->Fill(en, tdiff, GetEventWeight());
          }
          else
          {
            fhCellEnDiffColRowDiff[colDiffAbs][rowDiffAbs]->Fill(en, amp, GetEventWeight());
            fhTimeDiffClusDiffCellColRowDiff[colDiffAbs][rowDiffAbs]->Fill(en, tdiff, GetEventWeight());
          }
        }
      } // diff T-Card     
      
    } // Fill cell-cluster histogram loop
      
//    if ( en > 10 ) 
//    {
//      printf("total %d, same %d, diff %d, same5 %d, diff3 %d, same5+diff3 %d\n",
//             nCaloCellsPerCluster,nCellSame,nCellDiff,nCellSame5,nCellDiff5,nCellSame5+nCellDiff5);
//      
//      printf("E %2.2f, Esame %2.2f, Ediff %2.2f, Esame5 %2.2f, Ediff3 %2.2f\n",
//             en,enSame,enDiff,enSame5,enDiff5);
//    }
    
    fhNCellsPerClusterSameDiff ->Fill(en, nCellSame , nCellDiff , GetEventWeight());
    fhNCellsPerClusterSameDiffW->Fill(en, nCellSameW, nCellDiffW, GetEventWeight());

    Float_t frac  = (1.*nCellSame) /(nCaloCellsPerCluster-1.);
    
    fhNCellsPerClusterSameFrac    ->Fill(en, frac ,            GetEventWeight());
    fhNCellsPerClusterSameFracExo ->Fill(en, frac , exoticity, GetEventWeight());
    
    if ( fFillAllCellSameTCardHisto )
    {
      Float_t fracW = 0; 
      if ( nCellW > 0 ) fracW = (1.*nCellSameW)/(nCellW);
      
      fhNCellsPerClusterSameFracW   ->Fill(en, fracW,            GetEventWeight());
      fhNCellsPerClusterSameFracWExo->Fill(en, fracW, exoticity, GetEventWeight());
      
      if ( nCellDiff == 0 )
      {
        fhNCellsPerClusterAllSameTCard->Fill(en, nCaloCellsPerCluster, GetEventWeight());
        fhExoticityEClusAllSameTCard  ->Fill(en, exoticity           , GetEventWeight());
        if ( fFillExo50ns )
          fhExoticity50nsEClusAllSameTCard  ->Fill(en, exoticity50ns , GetEventWeight());

        fhM02EnergyAllSameTCard       ->Fill(en, m02                 , GetEventWeight());
        if ( en > fEMinForExo  )
          fhEtaPhiGridExoEnCutSameFracCut->Fill(icolMaxAbs, irowMaxAbs, exoticity, GetEventWeight());
      }
      else
      {
        if ( ebin >= 0 && ebin < fgkNEBins-1 )
          fhM02ExoNCellsNotAllSameTCard[ebin]->Fill(m02, exoticity, nCaloCellsPerCluster, GetEventWeight()); ;
      }
      
      if ( nCellDiffW == 0 )
      {
        fhNCellsPerClusterAllSameTCardW->Fill(en, nCaloCellsPerCluster, GetEventWeight());
        fhExoticityEClusAllSameTCardW  ->Fill(en, exoticity           , GetEventWeight());
        fhM02EnergyAllSameTCardW       ->Fill(en, m02                 , GetEventWeight());
        if ( en > fEMinForExo  )
          fhEtaPhiGridExoEnCutSameFracCutW->Fill(icolMaxAbs, irowMaxAbs, exoticity, GetEventWeight());
      }
      
      for(Int_t imin = 0; imin < nMinEnCut; imin++)
      {
        Float_t enCellMin  = 0.2 + imin*0.1;
        fhNCellsPerClusterMinEnCut    ->Fill(en, nCellDiffMinEn[imin]+nCellSameMinEn[imin], enCellMin, GetEventWeight());
        fhNCellsPerClusterDiffMinEnCut->Fill(en, nCellDiffMinEn[imin], enCellMin, GetEventWeight());
        fhNCellsPerClusterSameMinEnCut->Fill(en, nCellSameMinEn[imin], enCellMin, GetEventWeight());

        if ( nCellDiffMinEn[imin] == 0 )
        {
          fhNCellsPerClusterAllSameTCardMinEnCut->Fill(en, nCaloCellsPerCluster, enCellMin, GetEventWeight());
          fhExoticityEClusAllSameTCardMinEnCut  ->Fill(en, exoticity           , enCellMin, GetEventWeight());
          fhM02EnergyAllSameTCardMinEnCut       ->Fill(en, m02                 , enCellMin, GetEventWeight());
        }
      }
    } // fFillAllCellSameTCardHisto 
    
    fhNCellsPerClusterSame    ->Fill(en, nCellSame  , GetEventWeight());
    fhNCellsPerClusterDiff    ->Fill(en, nCellDiff  , GetEventWeight());
    if ( ebin >= 0 && ebin < fgkNEBins-1 )
    {
      fhNCellsSameDiffExo[ebin]->Fill(nCellSame, nCellDiff, exoticity, GetEventWeight());
      fhEnSameDiffExo    [ebin]->Fill(enSame   , enDiff   , exoticity, GetEventWeight());
    }
    
    if( fFillSameDiffFracHisto )
    {
      fhNCellsPerClusterSame5   ->Fill(en, nCellSame5 , GetEventWeight());
      fhNCellsPerClusterDiff5   ->Fill(en, nCellDiff5 , GetEventWeight());
      fhNCellsPerClusterSameW   ->Fill(en, nCellSameW , GetEventWeight());
      fhNCellsPerClusterDiffW   ->Fill(en, nCellDiffW , GetEventWeight());   
      
      fhNCellsPerClusterSameTimeDiff    ->Fill(en, nCellSameTimeDiff, GetEventWeight());
      fhNCellsPerClusterDiffTimeDiff    ->Fill(en, nCellDiffTimeDiff, GetEventWeight());  
      fhNCellsPerClusterSameDiffTimeDiff->Fill(en, nCellSameTimeDiff, nCellDiffTimeDiff, GetEventWeight());

      fhFracEnDiffSame      ->Fill(en, fracEnDiffSame      , GetEventWeight()); 
      fhFracNCellDiffSame   ->Fill(en, fracNCellDiffSame   , GetEventWeight());
      fhFracEnNCellDiffSame ->Fill(en, fracEnNCellDiffSame , GetEventWeight());
      
      fhFracEnDiffSameW     ->Fill(en, fracEnDiffSameW     , GetEventWeight()); 
      fhFracNCellDiffSameW  ->Fill(en, fracNCellDiffSameW  , GetEventWeight());
      fhFracEnNCellDiffSameW->Fill(en, fracEnNCellDiffSameW, GetEventWeight());
      
      fhFracEnDiffSame5     ->Fill(en, fracEnDiffSame5     , GetEventWeight()); 
      fhFracNCellDiffSame5  ->Fill(en, fracNCellDiffSame5  , GetEventWeight());
      fhFracEnNCellDiffSame5->Fill(en, fracEnNCellDiffSame5, GetEventWeight());
 
      fhFracEnDiffSameExo      ->Fill(en, fracEnDiffSame      , exoticity, GetEventWeight()); 
      fhFracNCellDiffSameExo   ->Fill(en, fracNCellDiffSame   , exoticity, GetEventWeight());
      fhFracEnNCellDiffSameExo ->Fill(en, fracEnNCellDiffSame , exoticity, GetEventWeight());
      
      fhFracEnDiffSameWExo     ->Fill(en, fracEnDiffSameW     , exoticity, GetEventWeight()); 
      fhFracNCellDiffSameWExo  ->Fill(en, fracNCellDiffSameW  , exoticity, GetEventWeight());
      fhFracEnNCellDiffSameWExo->Fill(en, fracEnNCellDiffSameW, exoticity, GetEventWeight());
      
      fhFracEnDiffSame5Exo     ->Fill(en, fracEnDiffSame5     , exoticity, GetEventWeight()); 
      fhFracNCellDiffSame5Exo  ->Fill(en, fracNCellDiffSame5  , exoticity, GetEventWeight());
      fhFracEnNCellDiffSame5Exo->Fill(en, fracEnNCellDiffSame5, exoticity, GetEventWeight());
      
      if ( en > fEMinForExo )
      {
        fhFracEnDiffSameEnCut      ->Fill(fracEnDiffSame      , GetEventWeight()); 
        fhFracNCellDiffSameEnCut   ->Fill(fracNCellDiffSame   , GetEventWeight());
        fhFracEnNCellDiffSameEnCut ->Fill(fracEnNCellDiffSame , GetEventWeight());
        
        fhFracEnDiffSameWEnCut     ->Fill(fracEnDiffSameW     , GetEventWeight()); 
        fhFracNCellDiffSameWEnCut  ->Fill(fracNCellDiffSameW  , GetEventWeight());
        fhFracEnNCellDiffSameWEnCut->Fill(fracEnNCellDiffSameW, GetEventWeight());
        
        fhFracEnDiffSame5EnCut     ->Fill(fracEnDiffSame5     , GetEventWeight()); 
        fhFracNCellDiffSame5EnCut  ->Fill(fracNCellDiffSame5  , GetEventWeight());
        fhFracEnNCellDiffSame5EnCut->Fill(fracEnNCellDiffSame5, GetEventWeight());
      }
      
      Float_t rEnSame  = 1-enSame /ampMax; 
      Float_t rEnDiff  = 1-enDiff /ampMax; 
      Float_t rEnSame5 = 1-enSame5/ampMax; 
      Float_t rEnDiff5 = 1-enDiff5/ampMax;
      
      fhExoSame ->Fill(en, rEnSame , GetEventWeight());
      fhExoDiff ->Fill(en, rEnDiff , GetEventWeight());
      fhExoSame5->Fill(en, rEnSame5, GetEventWeight());
      fhExoDiff5->Fill(en, rEnDiff5, GetEventWeight());
      
      if ( ebin >= 0 && ebin < fgkNEBins-1 )
      {
        if ( nCellSame > 0 && nCellDiff > 0 )
          fhEnNCellsSameDiffExo[ebin]->Fill(enSame/nCellSame, enDiff/nCellDiff, exoticity, GetEventWeight());
      }
    }
    
    // Shower shape
    //
    fhM02EnergyExo ->Fill(en, m02, exoticity , GetEventWeight());
    
    if ( m02 > 0.1 )
      fhM20EnergyExoM02MinCut->Fill(en, m20, exoticity, GetEventWeight());
    
    if ( ebin >= 0 && ebin < fgkNEBins-1 )
      fhM02ExoNCells[ebin]->Fill(m02, exoticity, nCaloCellsPerCluster, GetEventWeight()); 
    
    if ( fFillExo50ns )
    {
      fhM02EnergyExo50ns ->Fill(en, m02, exoticity50ns, GetEventWeight());
      if ( ebin >= 0 && ebin < fgkNEBins-1 )
        fhM02Exo50nsNCells[ebin]->Fill(m02, exoticity50ns, nCaloCellsPerCluster, GetEventWeight()); 
    }
    
    // Track matching
    //
    if ( matched && fFillMatchingHisto && track )
    {
      fhEOverPExo->Fill(en, eOverP, exoticity, GetEventWeight());
      
      if ( positive )
      {
        fhTrackMatchedDEtaPosExo->Fill(en, deta, exoticity, GetEventWeight());
        fhTrackMatchedDPhiPosExo->Fill(en, dphi, exoticity, GetEventWeight());
        
        if ( en > fEMinForExo ) 
        {
          fhTrackMatchedDEtaDPhiPosExo->Fill(deta, dphi, exoticity, GetEventWeight());
        }
      }
      else 
      {
        fhTrackMatchedDEtaNegExo->Fill(en, deta, exoticity, GetEventWeight());
        fhTrackMatchedDPhiNegExo->Fill(en, dphi, exoticity, GetEventWeight());
        
        if ( en > fEMinForExo ) 
        {
          fhTrackMatchedDEtaDPhiNegExo->Fill(deta, dphi, exoticity, GetEventWeight());
        }
      }
  
    } // matched
    
  } // Cluster loop
  
  // Counters per event
  fhNClusterPerEventNCellHigh20->Fill(nNCellW20   , GetEventWeight());      
  fhNClusterPerEventNCellHigh12->Fill(nNCellW12   , GetEventWeight());        
  fhNClusterPerEventExotic     ->Fill(nExotic     , GetEventWeight());
  fhNClusterPerEventExotic1Cell->Fill(nExotic1Cell, GetEventWeight());
  fhNClusterPerEventExoticNCell->Fill(nExoticNCell, GetEventWeight());
  
  fh2NClusterPerEventNCellHigh20->Fill(nCaloClusters, nNCellW20   , GetEventWeight());      
  fh2NClusterPerEventNCellHigh12->Fill(nCaloClusters, nNCellW12   , GetEventWeight());        
  fh2NClusterPerEventExotic     ->Fill(nCaloClusters, nExotic     , GetEventWeight());
  fh2NClusterPerEventExotic1Cell->Fill(nCaloClusters, nExotic1Cell, GetEventWeight());
  fh2NClusterPerEventExoticNCell->Fill(nCaloClusters, nExoticNCell, GetEventWeight());
  
  //if ( highestAmpMaxExo > 0 ) printf("Exotic final amp %f\n",highestAmpMaxExo);

  fh2NClusterPerEventExoticAmpMax->Fill(highestAmpMaxExo, nCaloClusters, GetEventWeight());
  
}

//____________________________________________________________________________
/// Fill basic clusters related histograms 
/// for different Event activity selections
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCaloExotics::ClusterHistogramsAfterEventCut
(const TObjArray *caloClusters, AliVCaloCells* cells)
{
  Int_t  nCaloClusters = caloClusters->GetEntriesFast() ;
  Int_t  bc            = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  // Get vertex, not needed.
  Double_t v[3] = {0,0,0};
    
  AliDebug(1,Form("In %s there are %d clusters", GetCalorimeterString().Data(), nCaloClusters));
  
  //if ( bc < 0 || bc > 3563 ) printf("Unexpected BC: %d\n", bc);
  
  fhEventCutBunchCrossing[0]->Fill(bc, GetEventWeight());
  fhEventCutVertexX      [0]->Fill(GetVertex(0)[0], GetEventWeight());
  fhEventCutVertexY      [0]->Fill(GetVertex(0)[1], GetEventWeight());
  fhEventCutVertexZ      [0]->Fill(GetVertex(0)[2], GetEventWeight());
  
  if ( fAcceptEvent )
  {
    fhEventCutBunchCrossing[1]->Fill(bc, GetEventWeight());
    fhEventCutVertexX      [1]->Fill(GetVertex(0)[0], GetEventWeight());
    fhEventCutVertexY      [1]->Fill(GetVertex(0)[1], GetEventWeight());
    fhEventCutVertexZ      [1]->Fill(GetVertex(0)[2], GetEventWeight());
  }
  
  if ( fEventNStripActive <= fEventMaxNumberOfStrips )
  {
    fhEventCutBunchCrossing[2]->Fill(bc, GetEventWeight());
    fhEventCutVertexX      [2]->Fill(GetVertex(0)[0], GetEventWeight());
    fhEventCutVertexY      [2]->Fill(GetVertex(0)[1], GetEventWeight());
    fhEventCutVertexZ      [2]->Fill(GetVertex(0)[2], GetEventWeight());
    
    if ( fAcceptEvent )
    {
      fhEventCutBunchCrossing[3]->Fill(bc, GetEventWeight());
      fhEventCutVertexX      [3]->Fill(GetVertex(0)[0], GetEventWeight());
      fhEventCutVertexY      [3]->Fill(GetVertex(0)[1], GetEventWeight());
      fhEventCutVertexZ      [3]->Fill(GetVertex(0)[2], GetEventWeight());
    }
  }
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++)
  {
    AliDebug(1,Form("Cluster: %d/%d, data %d",iclus+1,nCaloClusters,GetReader()->GetDataType()));
    
    AliVCluster* clus =  (AliVCluster*) caloClusters->At(iclus);
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    Float_t en  = fClusterMomentum.E();
    
    Int_t    nCaloCellsPerCluster = clus->GetNCells();
    //Int_t  nLabel = clus->GetNLabels();
    //Int_t *labels = clus->GetLabels();
    
    // Get the 1-fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    Int_t icolMax  = -1, icolMaxAbs = -1;
    Int_t irowMax  = -1, irowMaxAbs = -1;
    Int_t iRCU     = -1;
    //Int_t nSM  =
    GetModuleNumberCellIndexesAbsCaloMap(absIdMax,GetCalorimeter(), 
                                         icolMax   , irowMax, iRCU,
                                         icolMaxAbs, irowMaxAbs    );
    
    // Get time of max cell/cluster
    Double_t tmax  = cells->GetCellTime(absIdMax);
    tmax*=1.e9;
    tmax-=fConstantTimeShift;
    
    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    Float_t enOrg  = ampMax/(1-maxCellFraction); // cluster energy without non linearity correction
    
    Float_t exoticity  = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc)/ampMax;
    
    Float_t m02  = clus->GetM02();
    
//    if(fAcceptEvent)                                 
//      if ( fEventNStripActive <= fEventMaxNumberOfStrips )
//        if ( fAcceptEvent && fEventNStripActive <= fEventMaxNumberOfStrips )
//
    
    // Take only non exotic clusters
    if ( exoticity > 0.95 || m02 < 0.15 || nCaloCellsPerCluster < 3 ) continue;
    
    fhEventCutClusterEnergyTime[0]->Fill(en, tmax, GetEventWeight());
    
    if ( fAcceptEvent )
      fhEventCutClusterEnergyTime[1]->Fill(en, tmax, GetEventWeight());
    
    if ( fEventNStripActive <= fEventMaxNumberOfStrips )
    {
      fhEventCutClusterEnergyTime[2]->Fill(en, tmax, GetEventWeight());
      
      if ( fAcceptEvent )
        fhEventCutClusterEnergyTime[3]->Fill(en, tmax, GetEventWeight());
    }
    
    // Take clusters within good time window
    if ( tmax < fTimeCutMin || tmax> fTimeCutMax ) continue;

    fhEventCutClusterEnergy[0]->Fill(en, GetEventWeight());
    fhEventCutClusterEnergyECellMax[0]->Fill(enOrg, ampMax, GetEventWeight());
   
    if ( en > fEMinForExo ) 
      fhEventCutClusterEtaPhiGrid[0]->Fill(icolMaxAbs, irowMaxAbs, GetEventWeight());
   
    if ( fAcceptEvent )
    {
      fhEventCutClusterEnergy[1]->Fill(en, GetEventWeight());
      fhEventCutClusterEnergyECellMax[1]->Fill(enOrg, ampMax, GetEventWeight());
     
      if ( en > fEMinForExo ) 
        fhEventCutClusterEtaPhiGrid[1]->Fill(icolMaxAbs, irowMaxAbs, GetEventWeight());
    }
   
    if ( fEventNStripActive <= fEventMaxNumberOfStrips )
    {
      fhEventCutClusterEnergy[2]->Fill(en, GetEventWeight());
      fhEventCutClusterEnergyECellMax[2]->Fill(enOrg, ampMax, GetEventWeight());
      
      if ( en > fEMinForExo ) 
        fhEventCutClusterEtaPhiGrid[2]->Fill(icolMaxAbs, irowMaxAbs, GetEventWeight());
      
      if ( fAcceptEvent )
      {
        fhEventCutClusterEnergy[3]->Fill(en, GetEventWeight());
        fhEventCutClusterEnergyECellMax[3]->Fill(enOrg, ampMax, GetEventWeight());
        
        if ( en > fEMinForExo ) 
          fhEventCutClusterEtaPhiGrid[3]->Fill(icolMaxAbs, irowMaxAbs, GetEventWeight());
      } // fAcceptEvent
    } // good strips
    
  } // cluster loop
  
}


//____________________________________________________________________________
/// Fill basic track related histograms 
/// for different Event activity selections
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCaloExotics::TrackHistogramsAfterEventCut()
{
  Int_t  nTracks =  GetCTSTracks()->GetEntriesFast() ;  
  
  AliDebug(1,Form("There are %d clusters",  nTracks));
  
  // Loop over tracks
  for(Int_t itra = 0; itra < nTracks; itra++)
  {
    AliDebug(1,Form("Track: %d/%d",itra+1,nTracks));
    
    AliVTrack* track =  (AliVTrack*) GetCTSTracks()->At(itra);
    
    Float_t pt  = track->Pt();
    Float_t eta = track->Eta();
    Float_t phi = track->Phi ();
    
    fhEventCutTrackPt[0]->Fill(pt, GetEventWeight());
    
    if ( fAcceptEvent )
      fhEventCutTrackPt[1]->Fill(pt, GetEventWeight());
    
    if ( fEventNStripActive <= fEventMaxNumberOfStrips )
      fhEventCutTrackPt[2]->Fill(pt, GetEventWeight());
    
    if ( fEventNStripActive <= fEventMaxNumberOfStrips  && fAcceptEvent )
      fhEventCutTrackPt[3]->Fill(pt, GetEventWeight());
    
    if ( pt > 0.5 )
    {
      fhEventCutTrackEtaPhi[0]->Fill(eta, phi, GetEventWeight());
      
      if ( fAcceptEvent )
        fhEventCutTrackEtaPhi[1]->Fill(eta, phi, GetEventWeight());
      
      if ( fEventNStripActive <= fEventMaxNumberOfStrips )
        fhEventCutTrackEtaPhi[2]->Fill(eta, phi, GetEventWeight());
      
      if ( fEventNStripActive <= fEventMaxNumberOfStrips  && fAcceptEvent )
        fhEventCutTrackEtaPhi[3]->Fill(eta, phi, GetEventWeight());
    }
    
  } // track loop
  
}


//_________________________________________________
/// Save parameters used for analysis in a string.
//_________________________________________________
TObjString * AliAnaCaloExotics::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaCaloExotics ---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"Cell Amplitude > %2.1f GeV;",fCellAmpMin) ;
  parList+=onePar ;

  snprintf(onePar,buffersize,"E min Exo > %2.1f GeV;",fEMinForExo) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"F+ > %0.2f GeV;",fExoCut) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"NcellsW > %d;",fNCellHighCut) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"SM: nCell > %d, Sum E > %2.0f;",fHighNCellsCutSM,fHighEnergyCutSM) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"SM3: nCell < %d, Sum E < %2.0f;",fLowNCellsCutSM3,fLowEnergyCutSM3) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"Strip: nCell > %d-%d, Sum E > %2.0f-%2.0f; Event N Strips < %d;",
           fHighNCellsCutStrip[0], fHighNCellsCutStrip[1],
           fHighEnergyCutStrip[0], fHighEnergyCutStrip[1], fEventMaxNumberOfStrips) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"SM3 strips: nCell < %d, Sum E < %2.0f;",fLowNCellsCutSM3Strip,fLowEnergyCutSM3Strip) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"%2.0f < time < %2.0f ns;",fTimeCutMin,fTimeCutMax) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"time shift = %2.0f;",fConstantTimeShift) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"fill cell histo: %d;",fFillCellHisto) ;
  parList+=onePar ;

  snprintf(onePar,buffersize,"fill cluster with 1 cell histo: %d;",fFill1CellHisto) ;
  parList+=onePar ;

  snprintf(onePar,buffersize,"fill cells in event histo: %d;",fFillAllCellEventParamHisto) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"fill strip histo: %d;",fFillStripHisto) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"fill cluster track-matching histo: %d;",fFillMatchingHisto) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"fill cluster exoticity 50ns histo: %d;",fFillExo50ns) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"fill histo after event cuts: clusters %d; tracks %d",
           fFillClusterHistoAfterEventCut, fFillTrackHistoAfterEventCut) ;
  parList+=onePar ;
  
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
	
  return new TObjString(parList) ;
}

//___________________________________________________
/// Create histograms to be saved in output file and
/// store them in the output container.
//___________________________________________________
TList * AliAnaCaloExotics::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ExoticHistos") ; 
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  // Histogram binning and ranges
  // 

  // Energy bins
  Int_t nptbins = GetHistogramRanges()->GetHistoPtBins();           
  Float_t ptmax = GetHistogramRanges()->GetHistoPtMax();           
  Float_t ptmin = GetHistogramRanges()->GetHistoPtMin();
  
  TCustomBinning eBinning;
  eBinning.SetMinimum(0);
  eBinning.AddStep(15,0.5);                   // 30
  if(ptmax > 15)  eBinning.AddStep( 30, 1.0); // 15
  if(ptmax > 30)  eBinning.AddStep( 60, 2.5); // 12
  if(ptmax > 60)  eBinning.AddStep(100, 5.0); // 8 
  if(ptmax > 100) eBinning.AddStep(200,10.0); // 10
  if(ptmax > 200) eBinning.AddStep(300,20.0); // 5
  TArrayD eBinsArray;
  eBinning.CreateBinEdges(eBinsArray);
 
  TCustomBinning e2Binning;
  e2Binning.SetMinimum(0.1);
  e2Binning.AddStep(  2,0.10); // 19
  e2Binning.AddStep(  5,0.25); // 12
  if(ptmax >  5) e2Binning.AddStep( 10, 0.50); // 10
  if(ptmax > 10) e2Binning.AddStep( 25, 1.00); // 15
  if(ptmax > 25) e2Binning.AddStep( 50, 2.50); //  8
  if(ptmax > 50) e2Binning.AddStep(100, 5.00); // 10
  if(ptmax > 100)e2Binning.AddStep(200,10.00); // 10
  if(ptmax > 200)e2Binning.AddStep(300,25.00); // 4
  TArrayD e2BinsArray;
  e2Binning.CreateBinEdges(e2BinsArray);
  
  TCustomBinning eminBinning;
  eminBinning.SetMinimum(0.15);
  eminBinning.AddStep(2.05,0.1); 
  TArrayD eminBinsArray;
  eminBinning.CreateBinEdges(eminBinsArray);
  
  // Exoticity
//Int_t nexobins  = 201; Float_t exomin  = -1 ; Float_t exomax  = 1.01;

  TCustomBinning fBinning;
  fBinning.SetMinimum(-1.0);
  fBinning.AddStep(0.000,0.2000); // 5
  fBinning.AddStep(0.500,0.1000); // 5
  fBinning.AddStep(0.700,0.0500); // 4
  fBinning.AddStep(0.800,0.0200); // 5
  fBinning.AddStep(0.850,0.0100); // 5
  fBinning.AddStep(0.900,0.0050); // 20 
  fBinning.AddStep(1.002,0.0020); // 51 
  TArrayD fBinsArray;
  fBinning.CreateBinEdges(fBinsArray);
  
  // N cells in cluster
//Int_t nceclbins  = GetHistogramRanges()->GetHistoNClusterCellBins(); 
  Int_t   nceclmax = GetHistogramRanges()->GetHistoNClusterCellMax(); 
  Int_t   nceclmin = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  
  TCustomBinning nBinning;
  nBinning.SetMinimum(nceclmin-0.5);
  nBinning.AddStep(nceclmax+0.5,1);
  TArrayD nBinsArray;
  nBinning.CreateBinEdges(nBinsArray);

  TCustomBinning n2Binning;
  n2Binning.SetMinimum(nceclmin-0.5);
  n2Binning.AddStep(2*nceclmax+0.5,1);
  TArrayD n2BinsArray;
  n2Binning.CreateBinEdges(n2BinsArray);
  
  TCustomBinning nsameBinning;
  nsameBinning.SetMinimum(-0.5);
  nsameBinning.AddStep(16.5,1);
  TArrayD nsameBinsArray;
  nsameBinning.CreateBinEdges(nsameBinsArray);
  
  // Super-module bins
  Int_t totalSM = fLastModule-fFirstModule+1;

  TCustomBinning smBinning;
  smBinning.SetMinimum(fFirstModule-0.5);
  smBinning.AddStep(fLastModule+0.5,1); 
  TArrayD smBinsArray;
  smBinning.CreateBinEdges(smBinsArray);
  
  // Cell column-row histograms, see base class for data members setting
//Int_t   ncolcell    = fNMaxColsFull+2;
  Float_t colcellmin  = -1.5;
  Float_t colcellmax  = fNMaxColsFull+0.5;
  
//Int_t   nrowcell    = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  Float_t rowcellmin  = fNMaxRowsFullMin-1.5;
  Float_t rowcellmax  = fNMaxRowsFullMax+0.5;
  
  TCustomBinning rowBinning;
  rowBinning.SetMinimum(rowcellmin-1.5);
  rowBinning.AddStep(rowcellmax+0.5,1); 
  TArrayD rowBinsArray;
  rowBinning.CreateBinEdges(rowBinsArray);
  
  TCustomBinning colBinning;
  colBinning.SetMinimum(colcellmin-1.5);
  colBinning.AddStep(colcellmax+0.5,1);   
  TArrayD colBinsArray;
  colBinning.CreateBinEdges(colBinsArray);
  
  // Shower shape
//  Int_t ssbins  = GetHistogramRanges()->GetHistoShowerShapeBins();  
//  Float_t ssmax = GetHistogramRanges()->GetHistoShowerShapeMax();  
//  Float_t ssmin = GetHistogramRanges()->GetHistoShowerShapeMin();
  
  TCustomBinning ssBinning;
  ssBinning.SetMinimum(-0.005);
  ssBinning.AddStep(0.505,0.005); //101 
  ssBinning.AddStep(1.005,0.025); // 20
  ssBinning.AddStep(3.05,0.1);    // 20
  TArrayD ssBinsArray;
  ssBinning.CreateBinEdges(ssBinsArray);

  // Fractions
  TCustomBinning fracBinning;
  fracBinning.SetMinimum(-0.1);
  fracBinning.AddStep(1.1,0.01);   
  TArrayD fracBinsArray;
  fracBinning.CreateBinEdges(fracBinsArray);

  TCustomBinning frac2Binning;
  frac2Binning.SetMinimum(-0.1);
  frac2Binning.AddStep(2.1,0.01);   
  TArrayD frac2BinsArray;
  frac2Binning.CreateBinEdges(frac2BinsArray);
  
  // Cluster time bins
  Int_t ntimebins = GetHistogramRanges()->GetHistoTimeBins();         
  Float_t timemax = GetHistogramRanges()->GetHistoTimeMax();         
  Float_t timemin = GetHistogramRanges()->GetHistoTimeMin();  
  
  TCustomBinning tBinning;
  tBinning.SetMinimum(timemin);
  tBinning.AddStep(timemax,(timemax-timemin)/ntimebins);   
  TArrayD tBinsArray;
  tBinning.CreateBinEdges(tBinsArray);

  TCustomBinning t2Binning;
  t2Binning.SetMinimum(-610);
  t2Binning.AddStep(610,20);   
  TArrayD t2BinsArray;
  t2Binning.CreateBinEdges(t2BinsArray);
  
  // Difference in time within cluster
  Int_t tdbins  = GetHistogramRanges()->GetHistoDiffTimeBins() ;    
  Float_t tdmax = GetHistogramRanges()->GetHistoDiffTimeMax();     
  Float_t tdmin = GetHistogramRanges()->GetHistoDiffTimeMin();
  
  TCustomBinning tdBinning;
  tdBinning.SetMinimum(tdmin);
  tdBinning.AddStep(tdmax,(tdmax-tdmin)/tdbins);   
  TArrayD tdBinsArray;
  tdBinning.CreateBinEdges(tdBinsArray);
  
  // TM residuals
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  TCustomBinning retaBinning;
  retaBinning.SetMinimum(resetamin);
  retaBinning.AddStep(resetamax,(resetamax-resetamin)/nresetabins);   
  TArrayD retaBinsArray;
  retaBinning.CreateBinEdges(retaBinsArray);
 
  TCustomBinning rphiBinning;
  rphiBinning.SetMinimum(resphimin);
  rphiBinning.AddStep(resphimax,(resphimax-resphimin)/nresphibins);   
  TArrayD rphiBinsArray;
  rphiBinning.CreateBinEdges(rphiBinsArray);
     
  // E over p bins
  Int_t nPoverEbins  = GetHistogramRanges()->GetHistoEOverPBins();       
  Float_t eOverPmax  = GetHistogramRanges()->GetHistoEOverPMax();       
  Float_t eOverPmin  = GetHistogramRanges()->GetHistoEOverPMin();
  
  TCustomBinning eopBinning;
  eopBinning.SetMinimum(eOverPmin);
  eopBinning.AddStep(resphimax,(eOverPmax-eOverPmin)/nPoverEbins);   
  TArrayD eopBinsArray;
  eopBinning.CreateBinEdges(eopBinsArray);
  
  // Cluster size
  TCustomBinning sizeBinning;
  sizeBinning.SetMinimum(-8.5);
  sizeBinning.AddStep(8.5,1);   
  TArrayD sizeBinsArray;
  sizeBinning.CreateBinEdges(sizeBinsArray);
  
  //
  // Init histograms
  //
  
  // Event counters
  fhNClusterPerEventNCellHigh20 = new TH1F 
  ("hNClusterPerEventNCellHigh20",
   "#it{n}_{cluster} per event with #it{n}_{cell}^{#it{w}} > 20",
   1000,0,1000); 
  fhNClusterPerEventNCellHigh20->SetXTitle("#it{n}_{cluster}^{#it{n}_{cell}^{#it{w}} > 20}");
  fhNClusterPerEventNCellHigh20->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventNCellHigh20);
 
  fhNClusterPerEventNCellHigh12 = new TH1F 
  ("hNClusterPerEventNCellHigh12",
   "#it{n}_{cluster} per event with #it{n}_{cell}^{#it{w}} > 12",
   1000,0,1000); 
  fhNClusterPerEventNCellHigh12->SetXTitle("#it{n}_{cluster}^{#it{n}_{cell}^{#it{w}} > 12}");
  fhNClusterPerEventNCellHigh12->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventNCellHigh12);
  
  fhNClusterPerEventExotic = new TH1F 
  ("hNClusterPerEventExotic",
   "#it{n}_{cluster} per event with #it{E} > 5 GeV, #it{F}_{+} > 0.97",
   200,0,200); 
  fhNClusterPerEventExotic->SetXTitle("#it{n}_{cluster}^{#it{F}_{+} > 0.97}");
  fhNClusterPerEventExotic->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventExotic);
  
  fhNClusterPerEventExotic1Cell = new TH1F 
  ("hNClusterPerEventExotic1Cell",
   "#it{n}_{cluster} per event with #it{E} > 5 GeV, #it{n}_{cells} = 1",
   200,0,200); 
  fhNClusterPerEventExotic1Cell->SetXTitle("#it{n}_{cluster}^{#it{n}_{cells} = 1}");
  fhNClusterPerEventExotic1Cell->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventExotic1Cell);

  fhNClusterPerEventExoticNCell = new TH1F 
  ("hNClusterPerEventExoticNCell",
   "#it{n}_{cluster} per event with #it{E} > 5 GeV, #it{n}_{cells} > 1, #it{F}_{+} > 0.97",
   200,0,200); 
  fhNClusterPerEventExoticNCell->SetXTitle("#it{n}_{cluster}^{#it{F}_{+} > 0.97, #it{n}_{cells} > 1}");
  fhNClusterPerEventExoticNCell->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventExoticNCell);
  
  //////////
  
  fh2NClusterPerEventNCellHigh20 = new TH2F 
  ("h2NClusterPerEventNCellHigh20",
   "#it{n}_{cluster} per event with #it{n}_{cell}^{#it{w}} > 20",
   500,0,500,500,0,500); 
  fh2NClusterPerEventNCellHigh20->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventNCellHigh20->SetYTitle("#it{n}_{cluster}^{#it{n}_{cell}^{#it{w}} > 20}");
  fh2NClusterPerEventNCellHigh20->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventNCellHigh20);
  
  fh2NClusterPerEventNCellHigh12 = new TH2F 
  ("h2NClusterPerEventNCellHigh12",
   "#it{n}_{cluster} per event, total vs with #it{n}_{cell}^{#it{w}} > 12",
   500,0,500,500,0,500); 
  fh2NClusterPerEventNCellHigh12->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventNCellHigh12->SetYTitle("#it{n}_{cluster}^{#it{n}_{cell}^{#it{w}} > 12}");
  fh2NClusterPerEventNCellHigh12->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventNCellHigh12);
  
  fh2NClusterPerEventExotic = new TH2F 
  ("h2NClusterPerEventExotic",
   "#it{n}_{cluster} per event, total vs with #it{E} > 5 GeV, #it{F}_{+} > 0.97",
   200,0,200,200,0,200); 
  fh2NClusterPerEventExotic->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventExotic->SetYTitle("#it{n}_{cluster}^{#it{F}_{+} > 0.97}");
  fh2NClusterPerEventExotic->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventExotic);
  
  fh2NClusterPerEventExotic1Cell = new TH2F 
  ("h2NClusterPerEventExotic1Cell",
   "#it{n}_{cluster} per event, total vs with #it{E} > 5 GeV, #it{n}_{cells} = 1",
   200,0,200,200,0,200); 
  fh2NClusterPerEventExotic1Cell->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventExotic1Cell->SetYTitle("#it{n}_{cluster}^{#it{n}_{cells} = 1}");
  fh2NClusterPerEventExotic1Cell->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventExotic1Cell);
  
  fh2NClusterPerEventExoticNCell = new TH2F 
  ("h2NClusterPerEventExoticNCell",
   "#it{n}_{cluster} per event, total vs with #it{E} > 5 GeV, #it{n}_{cells} > 1, #it{F}_{+} > 0.97",
   200,0,200,200,0,200); 
  fh2NClusterPerEventExoticNCell->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventExoticNCell->SetYTitle("#it{n}_{cluster}^{#it{F}_{+} > 0.97, #it{n}_{cells} > 1}");
  fh2NClusterPerEventExoticNCell->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventExoticNCell);
  
  fh2NClusterPerEventExoticAmpMax  = new TH2F 
  ("h2NClusterPerEventExoticAmpMax",
   "#it{n}_{cluster} per event vs highest energy cluster cell with #it{F}_{+} > 0.97",
   nptbins,ptmin,ptmax,100,0,100); 
  fh2NClusterPerEventExoticAmpMax->SetXTitle("#it{E}_{cell}^{#it{F}_{+} > 0.97}");
  fh2NClusterPerEventExoticAmpMax->SetYTitle("#it{n}_{cluster}");
  fh2NClusterPerEventExoticAmpMax->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventExoticAmpMax);
  
  // Cluster Exoticity 2D
  //
  fhExoticityEClus = new TH2F 
  ("hExoticityEClus","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhExoticityEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClus->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClus);    
  
  fhExoticityEMaxCell = new TH2F 
  ("hExoticityEMaxCell","cell #it{F}_{+} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhExoticityEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhExoticityEMaxCell->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEMaxCell);    
  
  if ( fFillPerSMHisto )
  {
    fhExoticityEClusPerSM = new TH3F 
    ("hExoticityEClusPerSM","cell #it{F}_{+} vs #it{E}_{cluster}, vs SM, #it{n}_{cluster}^{cell} > 1",
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,  fBinsArray.GetArray(), 
     smBinsArray.GetSize() - 1, smBinsArray.GetArray());
    //nptbins,ptmin,ptmax, nexobins,exomin,exomax, totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhExoticityEClusPerSM->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoticityEClusPerSM->SetYTitle("#it{F}_{+}");
    fhExoticityEClusPerSM->SetZTitle("SM");
    outputContainer->Add(fhExoticityEClusPerSM);   
  }
  
  if ( fFill1CellHisto )
  {
    fhExoticity1Cell = new TH2F 
    ("hExoticity1Cell","cell #it{F}_{+} vs #it{E}, #it{n}_{cluster}^{cell} = 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhExoticity1Cell->SetXTitle("#it{E} (GeV) ");
    fhExoticity1Cell->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticity1Cell);    
  }
  
  if ( fFillExo50ns )
  {
    fhExoticity50nsEClus = new TH2F 
    ("hExoticity50nsEClus","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhExoticity50nsEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoticity50nsEClus->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticity50nsEClus);    
    
    if ( fFill1CellHisto )
    {
      fhExoticity50ns1Cell = new TH2F 
      ("hExoticity50ns1Cell","cell #it{F}_{+} vs #it{E}, #it{n}_{cluster}^{cell} = 1",
       //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
       eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
      fhExoticity50ns1Cell->SetXTitle("#it{E} (GeV) ");
      fhExoticity50ns1Cell->SetYTitle("#it{F}_{+}");
      outputContainer->Add(fhExoticity50ns1Cell);    
    }
  }
  
  // N cells per cluster
  //
  if ( fFillOpenTimeHisto )
  {
    fhNCellsPerClusterOpenTime  = new TH2F 
    ("hNCellsPerClusterOpenTime","#it{n}_{cells} vs #it{E}_{cluster}, no time cut",
     //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
      eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
     n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
    fhNCellsPerClusterOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterOpenTime->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterOpenTime);
    
    fhNCellsPerClusterWOpenTime  = new TH2F 
    ("hNCellsPerClusterWOpenTime","#it{n}_{cells} with w > 0 vs #it{E}_{cluster}, no time cut",
     //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
      eBinsArray.GetSize() - 1,   eBinsArray.GetArray(),  
     n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
    fhNCellsPerClusterWOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterWOpenTime->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterWOpenTime);
    
    fhNCellsPerClusterEMaxCellOpenTime  = new TH2F 
    ("hNCellsPerClusterEMaxCellOpenTime","#it{n}_{cells} vs #it{E}_{cell}^{max}, no time cut",
     //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
      eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
     n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
    fhNCellsPerClusterEMaxCellOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV)");
    fhNCellsPerClusterEMaxCellOpenTime->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterEMaxCellOpenTime);
    
    fhNCellsPerClusterWEMaxCellOpenTime  = new TH2F 
    ("hNCellsPerClusterWEMaxCellOpenTime","#it{n}_{cells} with w > 0 vs #it{E}_{cell}^{max}, no time cut",
     //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
      eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
     n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
    fhNCellsPerClusterWEMaxCellOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV)");
    fhNCellsPerClusterWEMaxCellOpenTime->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterWEMaxCellOpenTime);
  }
  
  fhNCellsPerCluster  = new TH2F 
  ("hNCellsPerCluster","#it{n}_{cells} vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
    eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
   n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
  fhNCellsPerCluster->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerCluster->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerCluster);

  fhNCellsPerClusterW  = new TH2F 
  ("hNCellsPerClusterW","#it{n}_{cells} with w > 0 vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
    eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
   n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
  fhNCellsPerClusterW->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterW->SetYTitle("#it{n}^{#it{w}}_{cells}");
  outputContainer->Add(fhNCellsPerClusterW);

  fhNCellsPerClusterTimeDiff  = new TH2F 
  ("hNCellsPerClusterTimeDiff","#it{n}_{cells} with #Delta #it{t} < 50 vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
    eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
   n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
  fhNCellsPerClusterTimeDiff->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterTimeDiff->SetYTitle("#it{n}^{#Delta #it{t} < 50}_{cells}");
  outputContainer->Add(fhNCellsPerClusterTimeDiff);
  
  fhNCellsPerClusterEMaxCell  = new TH2F 
  ("hNCellsPerClusterEMaxCell","#it{n}_{cells} vs #it{E}_{cell}^{max}",
   //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
    eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
   n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
  fhNCellsPerClusterEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV)");
  fhNCellsPerClusterEMaxCell->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerClusterEMaxCell);

  fhNCellsPerClusterWEMaxCell  = new TH2F 
  ("hNCellsPerClusterWEMaxCell","#it{n}_{cells} with w > 0 vs #it{E}_{cell}^{max}",
   //nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
   eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
   n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
  fhNCellsPerClusterWEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV)");
  fhNCellsPerClusterWEMaxCell->SetYTitle("#it{n}_{cells}^{#it{w}}");
  outputContainer->Add(fhNCellsPerClusterWEMaxCell);
  
  fhNCellsPerClusterExo  = new TH3F 
  ("hNCellsPerClusterExo","# cells per cluster vs #it{E}_{cluster} vs #it{F}_{+}",
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1, fBinsArray.GetArray());
   //nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
  fhNCellsPerClusterExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterExo->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterExo);
 
  if ( fFillExo50ns )
  {
    fhNCellsPerClusterExo50ns  = new TH3F 
    ("hNCellsPerClusterExo50ns","# cells per cluster vs #it{E}_{cluster} vs #it{F}_{+}",
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    //nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
    fhNCellsPerClusterExo50ns->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterExo50ns->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterExo50ns->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsPerClusterExo50ns);
  }
  
  if ( fFillPerSMHisto )
  {
    fhNCellsPerClusterPerSM  = new TH3F 
    ("hNCellsPerClusterPerSM","# cells per cluster vs #it{E}_{cluster} vs #it{F}_{+}",
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,  nBinsArray.GetArray(), 
     smBinsArray.GetSize() - 1, smBinsArray.GetArray());
    //nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhNCellsPerClusterPerSM->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterPerSM->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterPerSM->SetZTitle("SM");
    outputContainer->Add(fhNCellsPerClusterPerSM);
    
    fhNCellsPerClusterWPerSM  = new TH3F 
    ("hNCellsPerClusterWPerSM","# cells per cluster vs #it{E}_{cluster} vs #it{F}_{+}",
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,  nBinsArray.GetArray(), 
     smBinsArray.GetSize() - 1, smBinsArray.GetArray());
    //nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhNCellsPerClusterWPerSM->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterWPerSM->SetYTitle("#it{n}_{cells}^{#it{w}}");
    fhNCellsPerClusterWPerSM->SetZTitle("SM");
    outputContainer->Add(fhNCellsPerClusterWPerSM);
    
    for(Int_t imod = 0; imod < totalSM; imod++)
    {
      if(imod < fFirstModule || imod > fLastModule) continue;
      
      fhNCellsPerClusterExoPerSM[imod]  = new TH3F 
      (Form("hNCellsPerClusterExo_SM%d",imod),
       Form("# cells per cluster vs #it{E}_{cluster} vs exoticity, SM%d",imod),
       eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
       nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1, fBinsArray.GetArray());
      //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
      fhNCellsPerClusterExoPerSM[imod]->SetXTitle("#it{E}_{cluster} (GeV)");
      fhNCellsPerClusterExoPerSM[imod]->SetYTitle("#it{n}_{cells}");
      fhNCellsPerClusterExoPerSM[imod]->SetZTitle("#it{F}_{+}");
      outputContainer->Add(fhNCellsPerClusterExoPerSM[imod]);
    }
  }
  
  if ( fFillMatchingHisto )
  {
    fhExoticityEClusTrackMatch = new TH2F 
    ("hExoticityEClusTrackMatch","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1, track matched",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhExoticityEClusTrackMatch->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoticityEClusTrackMatch->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticityEClusTrackMatch);    
    
    fhNCellsPerClusterTrackMatch  = new TH2F 
    ("hNCellsPerClusterTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
     //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
     n2BinsArray.GetSize() - 1,  n2BinsArray.GetArray());
    fhNCellsPerClusterTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterTrackMatch->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterTrackMatch);
    
    fhNCellsPerClusterExoTrackMatch  = new TH3F 
    ("hNCellsPerClusterExoTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
     //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax, nexobinsS,exominS,exomaxS); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhNCellsPerClusterExoTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterExoTrackMatch->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterExoTrackMatch->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsPerClusterExoTrackMatch);
  }
  
  fhNCellsPerClusterM02  = new TH3F 
  ("hNCellsPerClusterM02","# cells per cluster vs #it{E}_{cluster} vs  #sigma^{2}_{long}",
   //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,100,0,0.5); 
    eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
   ssBinsArray.GetSize() - 1, ssBinsArray.GetArray());
  fhNCellsPerClusterM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterM02->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhNCellsPerClusterM02);
  
  // Different n cells in cluster depending T-Card 
  //
  fhNCellsPerClusterSameDiff  = new TH3F 
  ("hNCellsPerClusterSameDiff","#it{n}_{cells} in same vs different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, 17,0,17,nceclbins,nceclmin,nceclmax); 
       eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
   nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
       nBinsArray.GetSize() - 1,     nBinsArray.GetArray());
  fhNCellsPerClusterSameDiff->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameDiff->SetYTitle("#it{n}_{cells, same T-Card}");
  fhNCellsPerClusterSameDiff->SetZTitle("#it{n}_{cells, diff T-Card}");
  outputContainer->Add(fhNCellsPerClusterSameDiff);
  
  fhNCellsPerClusterSame  = new TH2F 
  ("hNCellsPerClusterSame","# cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, 17,0,17); 
   eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
   nsameBinsArray.GetSize() - 1,  nsameBinsArray.GetArray());
  fhNCellsPerClusterSame->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSame->SetYTitle("#it{n}_{cells, same T-Card}");
  outputContainer->Add(fhNCellsPerClusterSame);
  
  fhNCellsPerClusterDiff  = new TH2F 
  ("hNCellsPerClusterDiff","# cells per cluster in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhNCellsPerClusterDiff->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterDiff->SetYTitle("#it{n}_{cells, diff T-Card}");
  outputContainer->Add(fhNCellsPerClusterDiff);
  
  for(Int_t i = 0; i < fgkNEBins-1; i++) 
  {
    fhNCellsSameDiffExo[i] = new TH3F 
    (Form("hNCellsSameDiffExo_Ebin%d",i),
     Form("#it{n}_{cells-same} vs #it{n}_{cells-diff}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
     //17,0,17,nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
     nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,     nBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhNCellsSameDiffExo[i]->SetXTitle("#it{n}_{cells}^{same}");
    fhNCellsSameDiffExo[i]->SetYTitle("#it{n}_{cells}^{diff}");
    fhNCellsSameDiffExo[i]->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsSameDiffExo[i]);
    
    fhEnSameDiffExo[i] = new TH3F 
    (Form("hEnSameDiffExo_Ebin%d",i),
     Form("#Sigma #it{E}_{same}^{cells} vs #Sigma #it{E}_{diff}^{cells}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
     //200, 0, 20, 200, 0, 20, nexobinsS,exominS,exomaxS); 
     e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
     e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhEnSameDiffExo[i]->SetXTitle("#Sigma #it{E}_{same}^{cells} (GeV)");
    fhEnSameDiffExo[i]->SetYTitle("#Sigma #it{E}_{diff}^{cells} (GeV)");
    fhEnSameDiffExo[i]->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEnSameDiffExo[i]);
  } 
  
  // Only cells with weight
  fhNCellsPerClusterSameDiffW  = new TH3F 
  ("hNCellsPerClusterSameDiffW","#it{n}^{#it{w}}_{cells} in same vs different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, 17,0,17,nceclbins,nceclmin,nceclmax); 
   eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
   nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
   nBinsArray.GetSize() - 1,     nBinsArray.GetArray());
  fhNCellsPerClusterSameDiffW->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameDiffW->SetYTitle("#it{n}_{cells, same T-Card}^{#it{w}}");
  fhNCellsPerClusterSameDiffW->SetZTitle("#it{n}_{cells, diff T-Card}^{#it{w}}");
  outputContainer->Add(fhNCellsPerClusterSameDiffW);
  
  fhNCellsPerClusterSameFrac  = new TH2F 
  ("hNCellsPerClusterSameFrac","Fraction of # cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, 101,-0.005,1.005); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray()); 
  fhNCellsPerClusterSameFrac->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameFrac->SetYTitle("#it{n}_{cells, same T-Card} / (#it{n}_{cells}-1)");
  outputContainer->Add(fhNCellsPerClusterSameFrac);
  
  fhNCellsPerClusterSameFracExo  = new TH3F 
  ("hNCellsPerClusterSameFracExo","Fraction of # cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax, 101,-0.005,1.005, nexobins,exomin,exomax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
  fhNCellsPerClusterSameFracExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameFracExo->SetYTitle("#it{n}_{cells, same T-Card} / (#it{n}_{cells}-1)");
  fhNCellsPerClusterSameFracExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterSameFracExo); 
  
  if ( fFillSameDiffFracHisto )
  {
    fhNCellsPerClusterSameFracW  = new TH2F 
    ("hNCellsPerClusterSameFracW","Fraction of #it{n}^{#it{w}}_{cells} in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, 101,-0.005,1.005); 
     eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     fracBinsArray.GetSize() - 1, fracBinsArray.GetArray()); 
    fhNCellsPerClusterSameFracW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameFracW->SetYTitle("#it{n}^{#it{w}}_{cells, same T-Card} / #it{n}^{#it{w}}_{cells}");
    outputContainer->Add(fhNCellsPerClusterSameFracW);
    
    fhNCellsPerClusterSameFracWExo  = new TH3F 
    ("hNCellsPerClusterSameFracWExo","Fraction of #it{n}^{#it{w}}_{cells} in same T-Card as max #it{E} cell vs #it{E}_{cluster} vs #it{F}_{+}",
     //nptbins,ptmin,ptmax, 101,-0.005,1.005, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhNCellsPerClusterSameFracWExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameFracWExo->SetYTitle("#it{n}^{#it{w}}_{cells, same T-Card} / #it{n}^{#it{w}}_{cells}");
    fhNCellsPerClusterSameFracWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsPerClusterSameFracWExo); 
    
    fhNCellsPerClusterSame5  = new TH2F 
    ("hNCellsPerClusterSame5","# cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
     nptbins,ptmin,ptmax, 7,0,7); 
    fhNCellsPerClusterSame5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSame5->SetYTitle("#it{n}_{cells, same T-Card}");
    outputContainer->Add(fhNCellsPerClusterSame5);
    
    fhNCellsPerClusterDiff5  = new TH2F 
    ("hNCellsPerClusterDiff5","# cells per cluster in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
     nptbins,ptmin,ptmax, 7,0,7); 
    fhNCellsPerClusterDiff5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterDiff5->SetYTitle("#it{n}_{cells, diff T-Card}");
    outputContainer->Add(fhNCellsPerClusterDiff5);
    
    fhNCellsPerClusterSameW   = new TH2F 
    ("hNCellsPerClusterSameW","# cells per cluster with #it{w} in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, 17,0,17); 
        eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
    nsameBinsArray.GetSize() - 1,  nsameBinsArray.GetArray());
    fhNCellsPerClusterSameW ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameW ->SetYTitle("#it{n}_{cells, same T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterSameW );
    
    fhNCellsPerClusterDiffW   = new TH2F 
    ("hNCellsPerClusterDiffW","# cells per cluster with #it{w} in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhNCellsPerClusterDiffW ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterDiffW ->SetYTitle("#it{n}_{cells, diff T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterDiffW );
  
    fhNCellsPerClusterSameTimeDiff   = new TH2F 
    ("hNCellsPerClusterSameTimeDiff","# cells per cluster with #Delta #it{t}_{diff}<50 ns in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, 17,0,17); 
         eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
     nsameBinsArray.GetSize() - 1,  nsameBinsArray.GetArray());
    fhNCellsPerClusterSameTimeDiff ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameTimeDiff ->SetYTitle("#it{n}_{cells, same T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterSameTimeDiff );
    
    fhNCellsPerClusterDiffTimeDiff   = new TH2F 
    ("hNCellsPerClusterDiffTimeDiff","# cells per cluster with #Delta #it{t}_{diff}<50 in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhNCellsPerClusterDiffTimeDiff ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterDiffTimeDiff ->SetYTitle("#it{n}_{cells, diff T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterDiffTimeDiff );
    
    fhNCellsPerClusterSameDiffTimeDiff  = new TH3F 
    ("hNCellsPerClusterSameDiffTimeDiff","#it{n}^{#it{w}}_{cells} in same vs different T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, 17,0,17,nceclbins,nceclmin,nceclmax); 
         eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
     nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,     nBinsArray.GetArray());
    fhNCellsPerClusterSameDiffTimeDiff->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameDiffTimeDiff->SetYTitle("#it{n}_{cells, same T-Card}^{#Delta #it{t}<50}");
    fhNCellsPerClusterSameDiffTimeDiff->SetZTitle("#it{n}_{cells, diff T-Card}^{#Delta #it{t}<50}");
    outputContainer->Add(fhNCellsPerClusterSameDiffTimeDiff);
    
    // Cluster Exoticity other definitions
    //
    
    fhExoSame = new TH2F 
    ("hExoSame","cell #it{F}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoSame->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoSame->SetYTitle("#it{F}_{same}");
    outputContainer->Add(fhExoSame);    
    
    fhExoDiff = new TH2F 
    ("hExoDiff","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoDiff->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoDiff->SetYTitle("#it{F}_{diff}");
    outputContainer->Add(fhExoDiff);   
    
    fhExoSame5 = new TH2F 
    ("hExoSame5","cell #it{F}_{same-5} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoSame5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoSame5->SetYTitle("#it{F}_{same-5}");
    outputContainer->Add(fhExoSame5);    
    
    fhExoDiff5 = new TH2F 
    ("hExoDiff5","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoDiff5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoDiff5->SetYTitle("#it{F}_{diff-5}");
    outputContainer->Add(fhExoDiff5);   
    
    //
    fhFracEnDiffSame    = new TH2F 
    ("hFracEnDiffSame","cell #Sigma #it{E}_{diff}/#Sigma #it{E}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSame->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSame->SetYTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
    outputContainer->Add(fhFracEnDiffSame);  
    
    fhFracNCellDiffSame = new TH2F 
    ("hFracNCellDiffSame","cell #it{n}_{diff}/#it{n}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSame->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSame->SetYTitle("#it{n}_{diff}/#it{n}_{same}");
    outputContainer->Add(fhFracNCellDiffSame); 
    
    fhFracEnNCellDiffSame = new TH2F 
    ("hFracEnNCellDiffSame","cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSame->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSame->SetYTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
    outputContainer->Add(fhFracEnNCellDiffSame); 
    
    fhFracEnDiffSameW    = new TH2F 
    ("hFracEnDiffSameW","cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSameW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSameW->SetYTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
    outputContainer->Add(fhFracEnDiffSameW);  
    
    fhFracNCellDiffSameW = new TH2F 
    ("hFracNCellDiffSameW","cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSameW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSameW->SetYTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
    outputContainer->Add(fhFracNCellDiffSameW); 
    
    fhFracEnNCellDiffSameW = new TH2F 
    ("hFracEnNCellDiffSameW","cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSameW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSameW->SetYTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
    outputContainer->Add(fhFracEnNCellDiffSameW); 
    
    fhFracEnDiffSame5    = new TH2F 
    ("hFracEnDiffSame5","cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSame5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSame5->SetYTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
    outputContainer->Add(fhFracEnDiffSame5);  
    
    fhFracNCellDiffSame5 = new TH2F 
    ("hFracNCellDiffSame5","cell #it{n}_{diff}^{next}/#it{n}_{same}^{next} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSame5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSame5->SetYTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
    outputContainer->Add(fhFracNCellDiffSame5); 
    
    fhFracEnNCellDiffSame5 = new TH2F 
    ("hFracEnNCellDiffSame5","cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSame5->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSame5->SetYTitle("(#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next})");
    outputContainer->Add(fhFracEnNCellDiffSame5); 
    
    //
    fhFracEnDiffSameExo    = new TH3F 
    ("hFracEnDiffSameExo","cell #Sigma #it{E}_{diff}/#Sigma #it{E}_{same} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSameExo->SetYTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
    fhFracEnDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSameExo);  
    
    fhFracNCellDiffSameExo = new TH3F 
    ("hFracNCellDiffSameExo","cell #it{n}_{diff}/#it{n}_{same} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSameExo->SetYTitle("#it{n}_{diff}/#it{n}_{same}");
    fhFracNCellDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSameExo); 
    
    fhFracEnNCellDiffSameExo = new TH3F 
    ("hFracEnNCellDiffSameExo","cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    
    fhFracEnNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSameExo->SetYTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
    fhFracEnNCellDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnNCellDiffSameExo); 
    
    fhFracEnDiffSameWExo    = new TH3F 
    ("hFracEnDiffSameWExo","cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSameWExo->SetYTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
    fhFracEnDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSameWExo);  
    
    fhFracNCellDiffSameWExo = new TH3F 
    ("hFracNCellDiffSameWExo","cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSameWExo->SetYTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
    fhFracNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSameWExo); 
    
    fhFracEnNCellDiffSameWExo = new TH3F 
    ("hFracEnNCellDiffSameWExo","cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSameWExo->SetYTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
    fhFracEnNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnNCellDiffSameWExo); 
    
    fhFracEnDiffSame5Exo    = new TH3F 
    ("hFracEnDiffSame5Exo","cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray()); 
    fhFracEnDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnDiffSame5Exo->SetYTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
    fhFracEnDiffSame5Exo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSame5Exo);  
    
    fhFracNCellDiffSame5Exo = new TH3F 
    ("hFracNCellDiffSame5Exo","cell #it{n}_{diff}^{next}/#it{n}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracNCellDiffSame5Exo->SetYTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
    fhFracNCellDiffSame5Exo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSame5Exo); 
    
    fhFracEnNCellDiffSame5Exo = new TH3F 
    ("hFracEnNCellDiffSame5Exo","cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray()); 
    fhFracEnNCellDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhFracEnNCellDiffSame5Exo->SetYTitle("(#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next})");
    fhFracEnNCellDiffSame5Exo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnNCellDiffSame5Exo); 
    
    //
    fhFracEnDiffSameEnCut    = new TH1F 
    ("hFracEnDiffSameEnCut",
     Form("cell #Sigma #it{E}_{diff}/#Sigma #it{E}_{same}, #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracEnDiffSameEnCut->SetXTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
    outputContainer->Add(fhFracEnDiffSameEnCut);  
    
    fhFracNCellDiffSameEnCut = new TH1F 
    ("hFracNCellDiffSameEnCut",
     Form("cell #it{n}_{diff}/#it{n}_{same}, #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10); 
    fhFracNCellDiffSameEnCut->SetXTitle("#it{n}_{diff}/#it{n}_{same}");
    outputContainer->Add(fhFracNCellDiffSameEnCut); 
    
    fhFracEnNCellDiffSameEnCut = new TH1F 
    ("hFracEnNCellDiffSameEncut",
     Form("cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}), #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracEnNCellDiffSameEnCut->SetXTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
    outputContainer->Add(fhFracEnNCellDiffSameEnCut); 
    
    fhFracEnDiffSameWEnCut    = new TH1F 
    ("hFracEnDiffSameWEnCut",
     Form("cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracEnDiffSameWEnCut->SetXTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
    outputContainer->Add(fhFracEnDiffSameWEnCut);  
    
    fhFracNCellDiffSameWEnCut = new TH1F 
    ("hFracNCellDiffSameWEnCut",
     Form("cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracNCellDiffSameWEnCut->SetXTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
    outputContainer->Add(fhFracNCellDiffSameWEnCut); 
    
    fhFracEnNCellDiffSameWEnCut = new TH1F 
    ("hFracEnNCellDiffSameWEnCut",
     Form("cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})#it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10); 
    fhFracEnNCellDiffSameWEnCut->SetXTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
    outputContainer->Add(fhFracEnNCellDiffSameWEnCut); 
    
    fhFracEnDiffSame5EnCut    = new TH1F 
    ("hFracEnDiffSame5EnCut",
     Form("cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}, #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracEnDiffSame5EnCut->SetXTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
    outputContainer->Add(fhFracEnDiffSame5EnCut);  
    
    fhFracNCellDiffSame5EnCut = new TH1F 
    ("hFracNCellDiffSame5EnCut",
     Form("cell #it{n}_{diff}^{next}/#it{n}_{same}^{next}, #it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracNCellDiffSame5EnCut->SetXTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
    outputContainer->Add(fhFracNCellDiffSame5EnCut); 
    
    fhFracEnNCellDiffSame5EnCut = new TH1F 
    ("hFracEnNCellDiffSame5EnCut",
     Form("cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next})#it{E}_{cluster} > %2.1f, #it{n}_{cluster}^{cell} > 1",fEMinForExo),
     1000,0,10);
    fhFracEnNCellDiffSame5EnCut->SetXTitle("(#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next})");
    outputContainer->Add(fhFracEnNCellDiffSame5EnCut); 
    
    for(Int_t i = 0; i < fgkNEBins-1; i++) 
    {
      fhEnNCellsSameDiffExo[i] = new TH3F 
      (Form("hEnNCellsSameDiffExo_Ebin%d",i),
       Form("#Sigma #it{E}_{same}^{cells}/#it{n}_{cells}^{same} vs #Sigma #it{E}_{diff}^{cells}/#it{n}_{cells}^{diff}, %2.1f < #it{E} < %2.1f GeV",
            fEnergyBins[i],fEnergyBins[i+1]),
       //100, 0, 10, 100, 0, 10, nexobinsS,exominS,exomaxS); 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
      fhEnNCellsSameDiffExo[i]->SetXTitle("#Sigma #it{E}_{same}^{cells}/#it{n}_{cells}^{same}  (GeV)");
      fhEnNCellsSameDiffExo[i]->SetYTitle("#Sigma #it{E}_{diff}^{cells}/#it{n}_{cells}^{diff}  (GeV)");
      fhEnNCellsSameDiffExo[i]->SetZTitle("#it{F}_{+}");
      outputContainer->Add(fhEnNCellsSameDiffExo[i]);
    }
  }

  fhCellEnSameExo = new TH3F 
  ("hCellEnSameExo","#it{E}_{cluster} vs #it{E}_{cell}^{same} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhCellEnSameExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellEnSameExo->SetYTitle("#it{E}_{cell}^{same} (GeV)");
  fhCellEnSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnSameExo);    
  
  fhCellEnDiffExo = new TH3F 
  ("hCellEnDiffExo","#it{E}_{cluster} vs #it{E}_{cell}^{diff} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhCellEnDiffExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellEnDiffExo->SetYTitle("#it{E}_{cell}^{diff} (GeV)");
  fhCellEnDiffExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnDiffExo);     

  for(Int_t icoldiff = 0; icoldiff < 4; icoldiff++)
  {
    for(Int_t irowdiff = 0; irowdiff < 4; irowdiff++)
    {
      if ( irowdiff == 0 && icoldiff == 0 ) continue;
      
      fhCellEnDiffColRowDiff[icoldiff][irowdiff] = new TH2F 
      (Form("hCellEnDiff_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{diff}, #Delta col=%d - #Delta row=%d, #it{F}_{+}<%0.2f",icoldiff,irowdiff,fExoCut),
       //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray());
      fhCellEnDiffColRowDiff[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhCellEnDiffColRowDiff[icoldiff][irowdiff] ->SetYTitle("#it{E}_{cell}^{diff} (GeV)");
      outputContainer->Add(fhCellEnDiffColRowDiff[icoldiff][irowdiff] );     
      
      fhCellEnDiffColRowDiffExoCut[icoldiff][irowdiff] = new TH2F 
      (Form("hCellEnDiffExo_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{diff}, #Delta col=%d - #Delta row=%d, #it{F}_{+}>%0.2f",icoldiff,irowdiff,fExoCut),
       //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray());
      fhCellEnDiffColRowDiffExoCut[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhCellEnDiffColRowDiffExoCut[icoldiff][irowdiff] ->SetYTitle("#it{E}_{cell}^{diff} (GeV)");
      outputContainer->Add(fhCellEnDiffColRowDiffExoCut[icoldiff][irowdiff] );  
      
      fhTimeDiffClusDiffCellColRowDiff[icoldiff][irowdiff] = new TH2F 
      (Form("hTimeDiffClusDiffCell_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{diff}, #Delta col=%d - #Delta row=%d, #it{F}_{+}<%0.2f",icoldiff,irowdiff,fExoCut),
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       tdBinsArray.GetSize() - 1, tdBinsArray.GetArray());
      fhTimeDiffClusDiffCellColRowDiff[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhTimeDiffClusDiffCellColRowDiff[icoldiff][irowdiff] ->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
      outputContainer->Add(fhTimeDiffClusDiffCellColRowDiff[icoldiff][irowdiff] );     
      
      fhTimeDiffClusDiffCellColRowDiffExoCut[icoldiff][irowdiff] = new TH2F 
      (Form("hTimeDiffClusDiffCellExo_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #Delta #it{t}^{max-sec}, #Delta col=%d - #Delta row=%d, #it{F}_{+}>%0.2f",icoldiff,irowdiff,fExoCut),
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       tdBinsArray.GetSize() - 1, tdBinsArray.GetArray());
      fhTimeDiffClusDiffCellColRowDiffExoCut[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhTimeDiffClusDiffCellColRowDiffExoCut[icoldiff][irowdiff] ->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
      outputContainer->Add(fhTimeDiffClusDiffCellColRowDiffExoCut[icoldiff][irowdiff] );  
      
      if ( icoldiff > 1 ) continue;
      
      fhCellEnSameColRowDiff[icoldiff][irowdiff]  = new TH2F 
      (Form("hCellEnSame_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{same}, #Delta col=%d - #Delta row=%d, #it{F}_{+}<%0.2f",icoldiff,irowdiff,fExoCut),
       //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray());
      fhCellEnSameColRowDiff[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhCellEnSameColRowDiff[icoldiff][irowdiff] ->SetYTitle("#it{E}_{cell}^{same} (GeV)");
      outputContainer->Add(fhCellEnSameColRowDiff[icoldiff][irowdiff] );    
      
      fhCellEnSameColRowDiffExoCut[icoldiff][irowdiff]  = new TH2F 
      (Form("hCellEnSameExo_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{same}, #Delta col=%d - #Delta row=%d, #it{F}_{+}>%0.2f",icoldiff,irowdiff,fExoCut),
       //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
        eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       e2BinsArray.GetSize() - 1, e2BinsArray.GetArray());
      fhCellEnSameColRowDiffExoCut[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhCellEnSameColRowDiffExoCut[icoldiff][irowdiff] ->SetYTitle("#it{E}_{cell}^{same} (GeV)");
      outputContainer->Add(fhCellEnSameColRowDiffExoCut[icoldiff][irowdiff] ); 
      
      fhTimeDiffClusSameCellColRowDiff[icoldiff][irowdiff] = new TH2F 
      (Form("hTimeDiffClusSameCell_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #it{E}_{cell}^{diff}, #Delta col=%d - #Delta row=%d, #it{F}_{+}<%0.2f",icoldiff,irowdiff,fExoCut),
       eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       tdBinsArray.GetSize() - 1, tdBinsArray.GetArray());
      fhTimeDiffClusSameCellColRowDiff[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhTimeDiffClusSameCellColRowDiff[icoldiff][irowdiff] ->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
      outputContainer->Add(fhTimeDiffClusSameCellColRowDiff[icoldiff][irowdiff] );     
      
      fhTimeDiffClusSameCellColRowDiffExoCut[icoldiff][irowdiff] = new TH2F 
      (Form("hTimeDiffClusSameCellExo_DiffCol%d_DiffRow%d",icoldiff,irowdiff),
       Form("#it{E}_{cluster} vs #Delta #it{t}^{max-sec}, #Delta col=%d - #Delta row=%d, #it{F}_{+}>%0.2f",icoldiff,irowdiff,fExoCut),
       eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
       tdBinsArray.GetSize() - 1, tdBinsArray.GetArray());
      fhTimeDiffClusSameCellColRowDiffExoCut[icoldiff][irowdiff] ->SetXTitle("#it{E}_{cluster} (GeV)");
      fhTimeDiffClusSameCellColRowDiffExoCut[icoldiff][irowdiff] ->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
      outputContainer->Add(fhTimeDiffClusSameCellColRowDiffExoCut[icoldiff][irowdiff] );  
    }
  }
  
  if ( fFillOpenTimeHisto )
  {
    fhCellEnNCellWOpenTime = new TH3F 
    ("hCellEnNCellWOpenTime","#it{E}_{cluster} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}, no time cut",
     //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhCellEnNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhCellEnNCellWOpenTime->SetYTitle("#it{E}_{cell} (GeV)");
    fhCellEnNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
    outputContainer->Add(fhCellEnNCellWOpenTime);  
    
    fhCellEnNCellWEMaxOpenTime = new TH3F 
    ("hCellEnNCellWEMaxOpenTime","#it{E}_{cell}^{max} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}, no time cut",
     //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhCellEnNCellWEMaxOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV)");
    fhCellEnNCellWEMaxOpenTime->SetYTitle("#it{E}_{cell} (GeV)");
    fhCellEnNCellWEMaxOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
    outputContainer->Add(fhCellEnNCellWEMaxOpenTime);  
  }
  
  fhCellEnNCellW = new TH3F 
  ("hCellEnNCellW","#it{E}_{cluster} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellW->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellEnNCellW->SetYTitle("#it{E}_{cell} (GeV)");
  fhCellEnNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellW);  
 
  fhCellEnNCellWEMax = new TH3F 
  ("hCellEnNCellWEMax","#it{E}_{cell}^{max} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellWEMax->SetXTitle("#it{E}_{cell}^{max} (GeV)");
  fhCellEnNCellWEMax->SetYTitle("#it{E}_{cell} (GeV)");
  fhCellEnNCellWEMax->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellWEMax);  

  if ( fFillOpenTimeHisto )
  {
    fhCellTimeDiffNCellWOpenTime = new TH3F 
    ("hCellTimeDiffNCellWOpenTime","#it{E}_{cluster} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}, no time cut",
     //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhCellTimeDiffNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhCellTimeDiffNCellWOpenTime->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
    fhCellTimeDiffNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
    outputContainer->Add(fhCellTimeDiffNCellWOpenTime);  
    
    fhCellTimeDiffNCellWEMaxOpenTime = new TH3F 
    ("hCellTimeDiffNCellWEMaxOpenTime","#it{E}_{cell}^{max} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}, no time cut",
     //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhCellTimeDiffNCellWEMaxOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV)");
    fhCellTimeDiffNCellWEMaxOpenTime->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
    fhCellTimeDiffNCellWEMaxOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
    outputContainer->Add(fhCellTimeDiffNCellWEMaxOpenTime);  
  }
  
  fhCellTimeDiffNCellW = new TH3F 
  ("hCellTimeDiffNCellW","#it{E}_{cluster} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellW->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellTimeDiffNCellW->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
  fhCellTimeDiffNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellW);  
  
  fhCellTimeDiffNCellWEMax = new TH3F 
  ("hCellTimeDiffNCellWEMax","#it{E}_{cell}^{max} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellWEMax->SetXTitle("#it{E}_{cell}^{max} (GeV)");
  fhCellTimeDiffNCellWEMax->SetYTitle("#Delta #it{t}^{max-sec} (ns)");
  fhCellTimeDiffNCellWEMax->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellWEMax);  
  
  fhCellMaxClusterEn = new TH2F 
  ("hCellMaxClusterEn","#it{E}_{cluster}^{org} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray());
  fhCellMaxClusterEn->SetXTitle("#it{E}_{cluster}^{org} (GeV)");
  fhCellMaxClusterEn->SetYTitle("#it{E}_{cell}^{max} (GeV)");
  outputContainer->Add(fhCellMaxClusterEn);  

  fhCellMaxClusterEnRatio = new TH2F 
  ("hCellMaxClusterEnRatio","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 202,0,1.01); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray());
  fhCellMaxClusterEnRatio->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellMaxClusterEnRatio->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  outputContainer->Add(fhCellMaxClusterEnRatio);  
  
  if ( fFillOpenTimeHisto )
  {
    fhCellMaxClusterEnOpenTime = new TH2F 
    ("hCellMaxClusterEnOpenTime","#it{E}_{cluster}^{org} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1, no time cut",
     //nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray());
    fhCellMaxClusterEnOpenTime->SetXTitle("#it{E}_{cluster}^{org} (GeV)");
    fhCellMaxClusterEnOpenTime->SetYTitle("#it{E}_{cell}^{max} (GeV)");
    outputContainer->Add(fhCellMaxClusterEnOpenTime);    
    
    fhCellMaxClusterEnRatioOpenTime = new TH2F 
    ("hCellMaxClusterEnRatioOpenTime","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org}, #it{n}_{cluster}^{cell} > 1, no time cut",
     //nptbins,ptmin,ptmax, 202,0,1.01); 
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     fracBinsArray.GetSize() - 1, fracBinsArray.GetArray());
    fhCellMaxClusterEnRatioOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhCellMaxClusterEnRatioOpenTime->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
    outputContainer->Add(fhCellMaxClusterEnRatioOpenTime);  
    
    fhCellMaxClusterEnRatioNCellWOpenTime = new TH3F 
    ("hCellMaxClusterEnRatioNCellWOpenTime","#it{E}_{cluster}vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{n}_{cell}^{#it{w}}, #it{n}_{cluster}^{cell} > 1, no time cut",
     //nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
        nBinsArray.GetSize() - 1,    nBinsArray.GetArray());
    fhCellMaxClusterEnRatioNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV)");
    fhCellMaxClusterEnRatioNCellWOpenTime->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
    fhCellMaxClusterEnRatioNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
    outputContainer->Add(fhCellMaxClusterEnRatioNCellWOpenTime);  
  }
  
  fhCellMaxClusterEnRatioNCellW = new TH3F 
  ("hCellMaxClusterEnRatioNCellW","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{n}_{cell}^{#it{w}}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,    nBinsArray.GetArray());
  fhCellMaxClusterEnRatioNCellW->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellMaxClusterEnRatioNCellW->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  fhCellMaxClusterEnRatioNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellMaxClusterEnRatioNCellW);  
 
  fhCellMaxClusterEnRatioExo = new TH3F 
  ("hCellMaxClusterEnRatioExo","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 101,0,1.01, nexobins,exomin,exomax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
  fhCellMaxClusterEnRatioExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellMaxClusterEnRatioExo->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  fhCellMaxClusterEnRatioExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellMaxClusterEnRatioExo);  
  
  // Cluster acceptance
  //
  fhEtaPhiGridExoEnCut  = new TH3F 
  ("hEtaPhiGridExoEnCut",
   Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1",fEMinForExo),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
  fhEtaPhiGridExoEnCut->SetXTitle("column-#eta");
  fhEtaPhiGridExoEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridExoEnCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhEtaPhiGridExoEnCut);
  
  fhEtaPhiGridEnExoCut  = new TH3F 
  ("hEtaPhiGridEnExoCut", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{F}_{+} > %2.2f",fExoCut),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
  fhEtaPhiGridEnExoCut->SetXTitle("column-#eta");
  fhEtaPhiGridEnExoCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridEnExoCut->SetZTitle("#it{E} (GeV)");
  outputContainer->Add(fhEtaPhiGridEnExoCut);
  
  fhEtaPhiGridEnHighNCells  = new TH3F 
  ("hEtaPhiGridEnHighNCells", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}^{#it{w}} > %d",fNCellHighCut),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
  fhEtaPhiGridEnHighNCells->SetXTitle("column-#eta");
  fhEtaPhiGridEnHighNCells->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridEnHighNCells->SetZTitle("#it{E} (GeV)");
  outputContainer->Add(fhEtaPhiGridEnHighNCells);
 
  fhEtaPhiGridNCellEnCut  = new TH3F 
  ("hEtaPhiGridNCellEnCut",
   Form("colum (#eta) vs row (#varphi) vs #it{n}_{cells}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1",fEMinForExo),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nceclbins,nceclmin,nceclmax); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,   nBinsArray.GetArray());
  fhEtaPhiGridNCellEnCut->SetXTitle("column-#eta");
  fhEtaPhiGridNCellEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridNCellEnCut->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhEtaPhiGridNCellEnCut);
  
  if ( fFill1CellHisto )
  {
    fhEtaPhiGridEn1Cell  = new TH3F 
    ("hEtaPhiGridEn1Cell","colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}=1",
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins,ptmin,ptmax); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
       eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
    fhEtaPhiGridEn1Cell->SetXTitle("column-#eta");
    fhEtaPhiGridEn1Cell->SetYTitle("row-#varphi (rad)");
    fhEtaPhiGridEn1Cell->SetZTitle("#it{E} (GeV)");
    outputContainer->Add(fhEtaPhiGridEn1Cell);
  }
  
  // Timing and energy
  fhTimeEnergyExo  = new TH3F 
  ("hTimeEnergyExo","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}>1",
   //nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1, fBinsArray.GetArray());
  fhTimeEnergyExo->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyExo->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeEnergyExo);
 
  if ( fFill1CellHisto )
  {
    fhTimeEnergy1Cell  = new TH2F 
    ("hTimeEnergy1Cell","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}=1",
     //nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     tBinsArray.GetSize() - 1, tBinsArray.GetArray());
    fhTimeEnergy1Cell->SetXTitle("#it{E} (GeV)");
    fhTimeEnergy1Cell->SetYTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhTimeEnergy1Cell);
  }
  
  fhTimeDiffClusCellExo  = new TH3F 
  ("hTimeDiffClusCellExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1",
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhTimeDiffClusCellExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellExo);

  fhTimeDiffClusCellDiffTCardExo  = new TH3F 
  ("hTimeDiffClusCellDiffTCardExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1, diff T-Card",
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhTimeDiffClusCellDiffTCardExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellDiffTCardExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellDiffTCardExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellDiffTCardExo);

  fhTimeDiffClusCellSameTCardExo  = new TH3F 
  ("hTimeDiffClusCellSameTCardExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1, same T-Card",
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhTimeDiffClusCellSameTCardExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellSameTCardExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellSameTCardExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellSameTCardExo);
  
  fhTimeDiffWClusCellExo  = new TH3F 
  ("hTimeDiffWClusCellExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} for cells with w>0 vs #it{F}_{+}, #it{n}_{cells}>1",
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhTimeDiffWClusCellExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffWClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffWClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffWClusCellExo);
  
  fhTimeDiffAmpClusCellExo  = new TH3F 
  ("hTimeDiffAmpClusCellExo",
   Form("#it{E}_{cell i} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1, #it{E}_{cluster}>%2.1f GeV",fEMinForExo),
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhTimeDiffAmpClusCellExo->SetXTitle("#it{E}_{cell i} (GeV)");
  fhTimeDiffAmpClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffAmpClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffAmpClusCellExo);
  
  fhTimeEnergyM02  = new TH3F 
  ("hTimeEnergyM02","#it{E}_{cluster} vs #it{t}_{cluster} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   //nptbins,ptmin,ptmax, ntimebins,timemin,timemax, 100,0,0.5); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
    tBinsArray.GetSize() - 1,  tBinsArray.GetArray(), 
   ssBinsArray.GetSize() - 1, ssBinsArray.GetArray());
  fhTimeEnergyM02->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyM02->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeEnergyM02);
  
  fhTimeDiffClusCellM02  = new TH3F 
  ("hTimeDiffClusCellM02","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   //nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, 100,0,0.5); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
   ssBinsArray.GetSize() - 1, ssBinsArray.GetArray());
  fhTimeDiffClusCellM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellM02->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeDiffClusCellM02);
  
  fhTimeEnergyNCells  = new TH3F 
  ("hTimeEnergyNCells","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{n}_{cells}",
   //nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nceclbins,nceclmin,nceclmax); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
   nBinsArray.GetSize() - 1, nBinsArray.GetArray());
  fhTimeEnergyNCells->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyNCells->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyNCells->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhTimeEnergyNCells);

  fhTimeEnergyNCellsW  = new TH3F 
  ("hTimeEnergyNCellsW","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{n}_{cells} with #it{w} > 0",
   //nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nceclbins,nceclmin,nceclmax); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
   nBinsArray.GetSize() - 1, nBinsArray.GetArray());
  fhTimeEnergyNCellsW->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyNCellsW->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyNCellsW->SetZTitle("#it{n}_{cells}^{#it{w}}");
  outputContainer->Add(fhTimeEnergyNCellsW);

  fhTimeNCellCut  = new TH1F 
  ("hTimeNCellCut",Form("#it{t}_{cluster} for #it{n}_{cells}^{#it{w}} > %d",fNCellHighCut),
   4000,-1000,1000); 
  fhTimeNCellCut->SetXTitle("#it{t}_{cluster} (ns)");
  outputContainer->Add(fhTimeNCellCut);
  
  // Shower shape
  //
  fhM02EnergyExo  = new TH3F 
  ("hM02EnergyExo","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhM02EnergyExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyExo->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM02EnergyExo); 
  
  if ( fFillExo50ns )
  {
    fhM02EnergyExo50ns  = new TH3F 
    ("hM02EnergyExo50ns","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
     //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nexobinsS,exominS,exomaxS); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhM02EnergyExo50ns->SetXTitle("#it{E}_{cluster} (GeV)");
    fhM02EnergyExo50ns->SetYTitle("#sigma^{2}_{long}");
    fhM02EnergyExo50ns->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhM02EnergyExo50ns); 
  }
  
  fhM20EnergyExoM02MinCut  = new TH3F 
  ("hM20EnergyExoM02MinCut","#sigma^{2}_{short} vs #it{E}_{cluster} vs #it{F}_{+}, #sigma^{2}_{long} > 0.1",
   //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/2, nexobinsS,exominS,exomaxS); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhM20EnergyExoM02MinCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM20EnergyExoM02MinCut->SetYTitle("#sigma^{2}_{short}");
  fhM20EnergyExoM02MinCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM20EnergyExoM02MinCut);  
  
  for(Int_t i = 0; i < fgkNEBins-1; i++) 
  {
    fhM02ExoNCells[i] = new TH3F 
    (Form("hM02ExoNCells_Ebin%d",i),
     Form("#sigma^{2}_{long} vs #it{F}_{+} vs #it{n}_{cells}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
     //100,0,0.5,nexobinsS,exominS,exomaxS,nceclbins,nceclmin,nceclmax); 
     ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(),
      fBinsArray.GetSize() - 1,  fBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhM02ExoNCells[i]->SetXTitle("#sigma^{2}_{long}");
    fhM02ExoNCells[i]->SetYTitle("#it{F}_{+}");
    fhM02ExoNCells[i]->SetZTitle("#it{n}_{cells}");
    outputContainer->Add(fhM02ExoNCells[i]); 
    
    if ( fFillExo50ns )
    {
      fhM02Exo50nsNCells[i] = new TH3F 
      (Form("hM02Exo50nsNCells_Ebin%d",i),
       Form("#sigma^{2}_{long} vs #it{F}_{+} vs #it{n}_{cells}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
       //100,0,0.5,nexobinsS,exominS,exomaxS,nceclbins,nceclmin,nceclmax); 
       ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(),
        fBinsArray.GetSize() - 1,  fBinsArray.GetArray(), 
        nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
      fhM02Exo50nsNCells[i]->SetXTitle("#sigma^{2}_{long}");
      fhM02Exo50nsNCells[i]->SetYTitle("#it{F}_{+}");
      fhM02Exo50nsNCells[i]->SetZTitle("#it{n}_{cells}");
      outputContainer->Add(fhM02Exo50nsNCells[i]); 
    }
  }
  
  if (  fFillClusterColRowHisto )
  {
    if ( fFillPerSMHisto )
    {
      for(Int_t i = 0; i < fgkNEBins-1; i++)
      {
        fhClusterColRowPerSMHighNCell[i] = new TH3F 
        (Form("hClusterColRowPerSMHighNCell_Ebin%d",i),
         Form("column vs row vs SM, %2.1f < #it{E} < %2.1f GeV, #it{n}_{cells}^{#it{w}} > %d",
              fEnergyBins[i],fEnergyBins[i+1],fNCellHighCut),
         //17,-8.5,8.5,17,-8.5,8.5,totalSM,fFirstModule-0.5,fLastModule+0.5); 
         sizeBinsArray.GetSize() - 1, sizeBinsArray.GetArray(), 
         sizeBinsArray.GetSize() - 1, sizeBinsArray.GetArray(), 
           smBinsArray.GetSize() - 1,   smBinsArray.GetArray());
        fhClusterColRowPerSMHighNCell[i]->SetXTitle("column");
        fhClusterColRowPerSMHighNCell[i]->SetYTitle("row");
        fhClusterColRowPerSMHighNCell[i]->SetZTitle("SM");
        outputContainer->Add(fhClusterColRowPerSMHighNCell[i]);
      }
    }
    
    for(Int_t i = 0; i < fgkNEBins-1; i++)
    {
      for(Int_t j = 0; j < 2; j++) 
      {
        fhClusterColRowExo[j][i] = new TH3F 
        (Form("hClusterColRowExo_Ebin%d_Col%d",i,j),
         Form("column vs row vs #it{F}_{+}, %2.1f < #it{E} < %2.1f GeV, column %d",fEnergyBins[i],fEnergyBins[i+1],j),
         //17,-8.5,8.5,17,-8.5,8.5,nexobinsS,exominS,exomaxS); 
         sizeBinsArray.GetSize() - 1, sizeBinsArray.GetArray(), 
         sizeBinsArray.GetSize() - 1, sizeBinsArray.GetArray(), 
            fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
        fhClusterColRowExo[j][i]->SetXTitle("column");
        fhClusterColRowExo[j][i]->SetYTitle("row");
        fhClusterColRowExo[j][i]->SetZTitle("#it{F}_{+}");
        outputContainer->Add(fhClusterColRowExo[j][i]); 
        
        //      fhClusterColRow[j][i] = new TH2F 
        //      (Form("hClusterColRow_Ebin%d_Col%d",i,j),
        //       Form("column vs row, #it{F}_{+}<%2.2f, %2.1f < #it{E} < %2.1f GeV, column %d",fExoCut, fEnergyBins[i],fEnergyBins[i+1],j),
        //       17,-8.5,8.5,17,-8.5,8.5); 
        //      fhClusterColRow[j][i]->SetXTitle("column");
        //      fhClusterColRow[j][i]->SetYTitle("row");
        //      outputContainer->Add(fhClusterColRow[j][i]); 
        
        //      fhClusterColRowExoW [j][i] = new TH3F 
        //      (Form("hClusterColRowExoW_Ebin%d_Col%d",i,j),
        //       Form("column vs row vs #it{F}_{+}, %2.1f < #it{E} < %2.1f GeV, #it{w} > 0, column %d",fEnergyBins[i],fEnergyBins[i+1],j),
        //       17,-8.5,8.5,17,-8.5,8.5,nexobinsS,exominS,exomaxS); 
        //      fhClusterColRowExoW [j][i]->SetXTitle("column");
        //      fhClusterColRowExoW [j][i]->SetYTitle("row");
        //      fhClusterColRowExoW [j][i]->SetZTitle("#it{F}_{+}");
        //      outputContainer->Add(fhClusterColRowExoW [j][i]); 
      }
    }
  }
  
  if ( fFillExoEnMinCut )
  {
    fhExoticityECellMinCut = new TH3F 
    ("hExoticityECellMinCut","cell #it{F}_{+} for different #it{E}_{cell}^{min} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray(), 
     eminBinsArray.GetSize() - 1, eminBinsArray.GetArray());
    //nptbins,ptmin,ptmax, nexobins,exomin,exomax, 40,0.075,2.075); 
    fhExoticityECellMinCut->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoticityECellMinCut->SetYTitle("#it{F}_{+}");
    fhExoticityECellMinCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
    outputContainer->Add(fhExoticityECellMinCut); 
    
    fhExoticityWEClus = new TH2F 
    ("hExoticityWEClus","cell #it{F}_{+}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhExoticityWEClus->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoticityWEClus->SetYTitle("#it{F}_{+}^{#it{w}}");
    outputContainer->Add(fhExoticityWEClus);  
    
    fhNCellsPerClusterExoW  = new TH3F 
    ("hNCellsPerClusterExoW","# cells per cluster vs #it{E}_{cluster} vs #it{F}_{+}^{#it{w}}",
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1, nBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    //nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
    fhNCellsPerClusterExoW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterExoW->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterExoW->SetZTitle("#it{F}_{+}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterExoW);
    
    fhTimeEnergyExoW  = new TH3F 
    ("hTimeEnergyExoW","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}^{#it{w}}, #it{n}_{cells}>1",
     //nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhTimeEnergyExoW->SetXTitle("#it{E} (GeV)");
    fhTimeEnergyExoW->SetYTitle("#it{t}_{cluster} (ns)");
    fhTimeEnergyExoW->SetZTitle("#it{F}_{+}^{#it{w}}");
    outputContainer->Add(fhTimeEnergyExoW);
    
    fhM02EnergyExoW  = new TH3F 
    ("hM02EnergyExoW","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}^{#it{w}}",
     //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nexobinsS,exominS,exomaxS); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
    fhM02EnergyExoW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhM02EnergyExoW->SetYTitle("#sigma^{2}_{long}");
    fhM02EnergyExoW->SetZTitle("#it{F}_{+}^{#it{w}}");
    outputContainer->Add(fhM02EnergyExoW); 
  }
  
  fhNCellsPerClusterMinEnCut  = new TH3F 
  ("hNCellsPerClusterMinEnCut","# cells per cluster vs #it{E}_{cluster} vs #it{E}_{cell}^{min}",
      eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,     nBinsArray.GetArray(), 
   eminBinsArray.GetSize() - 1,  eminBinsArray.GetArray());
  fhNCellsPerClusterMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterMinEnCut->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
  outputContainer->Add(fhNCellsPerClusterMinEnCut);
  
  fhNCellsPerClusterSameMinEnCut  = new TH3F 
  ("hNCellsPerClusterSameMinEnCut","# cells per cluster in same T-Card vs #it{E}_{cluster} vs #it{E}_{cell}^{min}",
       eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
   nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
    eminBinsArray.GetSize() - 1,  eminBinsArray.GetArray());
  fhNCellsPerClusterSameMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameMinEnCut->SetYTitle("#it{n}_{cells}^{same}");
  fhNCellsPerClusterSameMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
  outputContainer->Add(fhNCellsPerClusterSameMinEnCut);
 
  fhNCellsPerClusterDiffMinEnCut  = new TH3F 
  ("hNCellsPerClusterDiffMinEnCut","# cells per cluster in diff. T-Card vs #it{E}_{cluster} vs #it{E}_{cell}^{min}",
      eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,     nBinsArray.GetArray(), 
   eminBinsArray.GetSize() - 1,  eminBinsArray.GetArray());
  fhNCellsPerClusterDiffMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterDiffMinEnCut->SetYTitle("#it{n}_{cells}^{diff}");
  fhNCellsPerClusterDiffMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
  outputContainer->Add(fhNCellsPerClusterDiffMinEnCut);
  
  if ( fFillAllCellSameTCardHisto )
  {
    for(Int_t i = 0; i < fgkNEBins-1; i++) 
    {
      fhM02ExoNCellsNotAllSameTCard[i] = new TH3F 
      (Form("hM02ExoNCellsNotAllSameTCard_Ebin%d",i),
       Form("#sigma^{2}_{long} vs #it{F}_{+} vs #it{n}_{cells}, %2.1f < #it{E} < %2.1f GeV, #it{n}_{cells-diff}^{#it{w}} > 0",fEnergyBins[i],fEnergyBins[i+1]),
       //100,0,0.5,nexobinsS,exominS,exomaxS,nceclbins,nceclmin,nceclmax); 
       ssBinsArray.GetSize() - 1, ssBinsArray.GetArray(),
        fBinsArray.GetSize() - 1,  fBinsArray.GetArray(), 
        nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
      fhM02ExoNCellsNotAllSameTCard[i]->SetXTitle("#sigma^{2}_{long}");
      fhM02ExoNCellsNotAllSameTCard[i]->SetYTitle("#it{F}_{+}");
      fhM02ExoNCellsNotAllSameTCard[i]->SetZTitle("#it{n}_{cells}");
      outputContainer->Add(fhM02ExoNCellsNotAllSameTCard[i]); 
    }
    
    fhExoticityEClusAllSameTCard = new TH2F 
    ("hExoticityEClusAllSameTCard","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoticityEClusAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoticityEClusAllSameTCard->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticityEClusAllSameTCard);    
    
    if ( fFillExo50ns )
    {
      fhExoticity50nsEClusAllSameTCard = new TH2F 
      ("hExoticity50nsEClusAllSameTCard","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}",
       //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
       eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1, fBinsArray.GetArray());
      fhExoticity50nsEClusAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
      fhExoticity50nsEClusAllSameTCard->SetYTitle("#it{F}_{+}");
      outputContainer->Add(fhExoticity50nsEClusAllSameTCard);    
    }
    
    fhNCellsPerClusterAllSameTCard  = new TH2F 
    ("hNCellsPerClusterAllSameTCard","# cells per cluster vs #it{E}_{cluster}, #it{n}_{cells}=#it{n}_{cells-same}",
     //nptbins,ptmin,ptmax, 17,0,17); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray());
    fhNCellsPerClusterAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterAllSameTCard->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterAllSameTCard);
    
    fhM02EnergyAllSameTCard  = new TH2F 
    ("hM02EnergyNCellAllSameTCard","#sigma^{2}_{long} vs #it{E}_{cluster}, #it{n}_{cells} = #it{n}_{cells-same}",
     //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     ssBinsArray.GetSize() - 1, ssBinsArray.GetArray());
    fhM02EnergyAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
    fhM02EnergyAllSameTCard->SetYTitle("#sigma^{2}_{long}");
    outputContainer->Add(fhM02EnergyAllSameTCard); 
    
    fhEtaPhiGridExoEnCutSameFracCut = new TH3F 
    ("hEtaPhiGridExoEnCutSameFracCut",
     Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1, #it{n}_{cells-diff}^{#it{w}} = 0",fEMinForExo),
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
    fhEtaPhiGridExoEnCutSameFracCut->SetXTitle("column-#eta");
    fhEtaPhiGridExoEnCutSameFracCut->SetYTitle("row-#varphi (rad)");
    fhEtaPhiGridExoEnCutSameFracCut->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEtaPhiGridExoEnCutSameFracCut);
    
    //
    
    fhExoticityEClusAllSameTCardW = new TH2F 
    ("hExoticityEClusAllSameTCardW","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}, #it{n}_{cells-diff}^{#it{w}} = 0",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoticityEClusAllSameTCardW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoticityEClusAllSameTCardW->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticityEClusAllSameTCardW);    
    
    fhNCellsPerClusterAllSameTCardW  = new TH2F 
    ("hNCellsPerClusterAllSameTCardW","# cells per cluster vs #it{E}_{cluster}, #it{n}_{cells}=#it{n}_{cells-same}, #it{n}_{cells-diff}^{#it{w}} = 0",
     //nptbins,ptmin,ptmax, 17,0,17); 
         eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
         nBinsArray.GetSize() - 1, nBinsArray.GetArray());
    fhNCellsPerClusterAllSameTCardW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterAllSameTCardW->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterAllSameTCardW);
    
    fhM02EnergyAllSameTCardW  = new TH2F 
    ("hM02EnergyNCellAllSameTCardW","#sigma^{2}_{long} vs #it{E}_{cluster}, #it{n}_{cells} = #it{n}_{cells-same}, #it{n}_{cells-diff}^{#it{w}} = 0",
     //nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     ssBinsArray.GetSize() - 1, ssBinsArray.GetArray());
    fhM02EnergyAllSameTCardW->SetXTitle("#it{E}_{cluster} (GeV)");
    fhM02EnergyAllSameTCardW->SetYTitle("#sigma^{2}_{long}");
    outputContainer->Add(fhM02EnergyAllSameTCardW); 
    
    fhEtaPhiGridExoEnCutSameFracCutW = new TH3F 
    ("hEtaPhiGridExoEnCutSameFracCutW",
     Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1, #it{n}_{cells-diff}^{#it{w}} = 0",fEMinForExo),
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
    fhEtaPhiGridExoEnCutSameFracCutW->SetXTitle("column-#eta");
    fhEtaPhiGridExoEnCutSameFracCutW->SetYTitle("row-#varphi (rad)");
    fhEtaPhiGridExoEnCutSameFracCutW->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEtaPhiGridExoEnCutSameFracCutW);

    //
    
    fhExoticityEClusAllSameTCardMinEnCut = new TH3F 
    ("hExoticityEClusAllSameTCardMinEnCut","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}, #it{n}_{cells-diff}^{E cut} = 0",
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray(),
     eminBinsArray.GetSize() - 1, eminBinsArray.GetArray());
    fhExoticityEClusAllSameTCardMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
    fhExoticityEClusAllSameTCardMinEnCut->SetYTitle("#it{F}_{+}");
    fhExoticityEClusAllSameTCardMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
    outputContainer->Add(fhExoticityEClusAllSameTCardMinEnCut);    
    
    fhNCellsPerClusterAllSameTCardMinEnCut  = new TH3F 
    ("hNCellsPerClusterAllSameTCardMinEnCut","# cells per cluster vs #it{E}_{cluster}, #it{n}_{cells}=#it{n}_{cells-same}, #it{n}_{cells-diff}^{E cut} = 0",
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
         nBinsArray.GetSize() - 1,     nBinsArray.GetArray(),
      eminBinsArray.GetSize() - 1,  eminBinsArray.GetArray());
    fhNCellsPerClusterAllSameTCardMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterAllSameTCardMinEnCut->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterAllSameTCardMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
    outputContainer->Add(fhNCellsPerClusterAllSameTCardMinEnCut);
    
    fhM02EnergyAllSameTCardMinEnCut  = new TH3F 
    ("hM02EnergyNCellAllSameTCardMinEnCut","#sigma^{2}_{long} vs #it{E}_{cluster}, #it{n}_{cells} = #it{n}_{cells-same}, #it{n}_{cells-diff}^{E cut} = 0",
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
       ssBinsArray.GetSize() - 1,   ssBinsArray.GetArray(),
     eminBinsArray.GetSize() - 1, eminBinsArray.GetArray());
    fhM02EnergyAllSameTCardMinEnCut->SetXTitle("#it{E}_{cluster} (GeV)");
    fhM02EnergyAllSameTCardMinEnCut->SetYTitle("#sigma^{2}_{long}");
    fhM02EnergyAllSameTCardMinEnCut->SetZTitle("#it{E}_{cell}^{min} (GeV)");
    outputContainer->Add(fhM02EnergyAllSameTCardMinEnCut); 
  }
  
  // Track matching
  //
  if ( fFillMatchingHisto )
  {
    fhTrackMatchedDEtaNegExo  = new TH3F
    ("hTrackMatchedDEtaNegExo","d#eta of cluster-negative track vs cluster energy vs #it{F}_{+}",
     //nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     retaBinsArray.GetSize() - 1, retaBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDEtaNegExo->SetYTitle("d#eta");
    fhTrackMatchedDEtaNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDEtaNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDPhiNegExo  = new TH3F
    ("hTrackMatchedDPhiNegExo","d#varphi of cluster-negative track vs cluster energy vs #it{F}_{+}",
     //nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     rphiBinsArray.GetSize() - 1, rphiBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDPhiNegExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaDPhiNegExo  = new TH3F
    ("hTrackMatchedDEtaDPhiNegExo",
     Form("d#eta vs d#varphi of cluster- negative track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
     //nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
     retaBinsArray.GetSize() - 1, retaBinsArray.GetArray(), 
     rphiBinsArray.GetSize() - 1, rphiBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDEtaDPhiNegExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiNegExo->SetXTitle("d#eta");
    fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaPosExo  = new TH3F
    ("hTrackMatchedDEtaPosExo","d#eta of cluster-positive track vs cluster energy vs #it{F}_{+}",
     //nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     retaBinsArray.GetSize() - 1, retaBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDEtaPosExo->SetYTitle("d#eta");
    fhTrackMatchedDEtaPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDEtaPosExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDPhiPosExo  = new TH3F
    ("hTrackMatchedDPhiPosExo","d#varphi of cluster-positive track vs cluster energy vs #it{F}_{+}",
     //nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
        eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
     rphiBinsArray.GetSize() - 1, rphiBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDPhiPosExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaDPhiPosExo  = new TH3F
    ("hTrackMatchedDEtaDPhiPosExo",
     Form("d#eta vs d#varphi of cluster-positive track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
     //nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
     retaBinsArray.GetSize() - 1, retaBinsArray.GetArray(), 
     rphiBinsArray.GetSize() - 1, rphiBinsArray.GetArray(), 
        fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
    fhTrackMatchedDEtaDPhiPosExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiPosExo->SetXTitle("d#eta");
    fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhEOverPExo = new TH3F
    ("hEOverPExo",
    "Track matches #it{E}/#it{p} vs #it{F}_{+}",
     //nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax, nexobinsS,exominS,exomaxS);
       eBinsArray.GetSize() - 1,   eBinsArray.GetArray(), 
     eopBinsArray.GetSize() - 1, eopBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
    fhEOverPExo->SetYTitle("#it{E}/#it{p}");
    fhEOverPExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhEOverPExo->SetZTitle("#it{F}_{+}");
    
    outputContainer->Add(fhTrackMatchedDEtaNegExo) ;
    outputContainer->Add(fhTrackMatchedDPhiNegExo) ;
    outputContainer->Add(fhTrackMatchedDEtaPosExo) ;
    outputContainer->Add(fhTrackMatchedDPhiPosExo) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiNegExo) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiPosExo) ;
    outputContainer->Add(fhEOverPExo);
    
    
    if ( fFill1CellHisto )
    {
      fhTrackMatchedDEtaNeg1Cell  = new TH2F
      ("hTrackMatchedDEtaNeg1Cell",
       "d#eta of cluster-negative track vs cluster energy, #it{n}_{cell}=1",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNeg1Cell->SetYTitle("d#eta");
      fhTrackMatchedDEtaNeg1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNeg1Cell  = new TH2F
      ("hTrackMatchedDPhiNeg1Cell",
       "d#varphi of cluster-negative track vs cluster energy, #it{n}_{cell}=1",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNeg1Cell->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDPhiNeg1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiNeg1Cell  = new TH2F
      ("hTrackMatchedDEtaDPhiNeg1Cell",
       Form("d#eta vs d#varphi of cluster-negative track, E > %2.2f, #it{n}_{cell}=1",fEMinForExo),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiNeg1Cell->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDEtaDPhiNeg1Cell->SetXTitle("d#eta");
      
      fhTrackMatchedDEtaPos1Cell  = new TH2F
      ("hTrackMatchedDEtaPos1Cell",
       "d#eta of cluster-positive track vs cluster energy, #it{n}_{cell}=1",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaPos1Cell->SetYTitle("d#eta");
      fhTrackMatchedDEtaPos1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiPos1Cell  = new TH2F
      ("hTrackMatchedDPhiPos1Cell",
       "d#varphi of cluster-positive track vs cluster energy, #it{n}_{cell}=1",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiPos1Cell->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDPhiPos1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiPos1Cell  = new TH2F
      ("hTrackMatchedDEtaDPhiPos1Cell",
       Form("d#eta vs d#varphi of cluster-positive track, E > %2.2f, #it{n}_{cell}=1",fEMinForExo),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiPos1Cell->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDEtaDPhiPos1Cell->SetXTitle("d#eta");
      
      fhEOverP1Cell = new TH2F
      ("hEOverP1Cell",
       "Track matches #it{E}/#it{p}, #it{n}_{cell}=1",
       nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
      fhEOverP1Cell->SetYTitle("#it{E}/#it{p}");
      fhEOverP1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNeg1Cell) ;
      outputContainer->Add(fhTrackMatchedDPhiNeg1Cell) ;
      outputContainer->Add(fhTrackMatchedDEtaPos1Cell) ;
      outputContainer->Add(fhTrackMatchedDPhiPos1Cell) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhiNeg1Cell) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhiPos1Cell) ;
      outputContainer->Add(fhEOverP1Cell);
    }
  }
  
  // Calorimeter cells
  //
  if ( fFillCellHisto )
  {
    fhCellExoAmp     = new TH2F 
    ("hCellExoAmp","cell #it{F}_{+} vs #it{E}_{cell}",
     //nptbins,ptmin,ptmax/2, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhCellExoAmp->SetXTitle("#it{E}_{cell} (GeV)");
    fhCellExoAmp->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoAmp);    
    
    if ( fFillExo50ns )
    {
      fhCellExo50nsAmp     = new TH2F 
      ("hCellExo50nsAmp","cell #it{F}_{+} vs #it{E}_{cell}",
       //nptbins,ptmin,ptmax/2, nexobins,exomin,exomax); 
       eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1, fBinsArray.GetArray());
      fhCellExo50nsAmp->SetXTitle("#it{E}_{cell} (GeV)");
      fhCellExo50nsAmp->SetYTitle("#it{F}_{+}");
      outputContainer->Add(fhCellExo50nsAmp);    
    }
    
    fhCellExoAmpTime = new TH3F 
    ("hCellExoAmpTime","Cell #it{F}_{+} vs #it{E}_{cell} vs time",
     //nptbins,ptmin,ptmax/2, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhCellExoAmpTime->SetXTitle("#it{E}_{cell} (GeV)");
    fhCellExoAmpTime->SetYTitle("#it{t}_{cell} (ns)");
    fhCellExoAmpTime->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoAmpTime);    
    
    fhCellExoGrid    = new TH3F 
    ("hCellExoGrid",
     Form("Cell hits row-column vs #it{F}_{+} for #it{E}_{cell} > %2.1f",fEMinForExo), 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, nexobinsS,exominS,exomaxS); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
    fhCellExoGrid->SetYTitle("row (phi direction)");
    fhCellExoGrid->SetXTitle("column (eta direction)");
    fhCellExoGrid->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoGrid);
    
    fhCellGridTimeHighNCell20    = new TH3F 
    ("hCellGridTimeHighNCell20",
     "Cell hits row-column vs cluster time for #it{n}_{cell}^{#it{w}} > 20", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTimeHighNCell20->SetYTitle("row (phi direction)");
    fhCellGridTimeHighNCell20->SetXTitle("column (eta direction)");
    fhCellGridTimeHighNCell20->SetZTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhCellGridTimeHighNCell20);
    
    fhCellGridTimeHighNCell12    = new TH3F 
    ("hCellGridTimeHighNCell12",
     "Cell hits row-column vs cluster time for #it{n}_{cell}^{#it{w}} > 12", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTimeHighNCell12->SetYTitle("row (phi direction)");
    fhCellGridTimeHighNCell12->SetXTitle("column (eta direction)");
    fhCellGridTimeHighNCell12->SetZTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhCellGridTimeHighNCell12);
    
    for(Int_t icut = 0; icut < fgkNCellEnMinBins; icut++)
    {
      fhCellGridTime[icut]    = new TH3F 
      (Form("hCellGridTime_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("Cell hits row-column vs cell time for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
       rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
        t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
      fhCellGridTime[icut]->SetYTitle("row (phi direction)");
      fhCellGridTime[icut]->SetXTitle("column (eta direction)");
      fhCellGridTime[icut]->SetZTitle("#it{t}_{cell} (ns)");
      outputContainer->Add(fhCellGridTime[icut]);
      
      // Per Strip
      //
      if ( fFillStripHisto )
      {
        // Count cells in strip
        fhNCellsPerStrip[icut] = new TH1F 
        (Form("hNCellsPerStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5);
        fhNCellsPerStrip[icut]->SetYTitle("Counts per event");
        fhNCellsPerStrip[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStrip[icut]);   
        
        fhNCellsPerStripAcceptEventStrip[icut] = new TH1F 
        (Form("hNCellsPerStripAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5);
        fhNCellsPerStripAcceptEventStrip[icut]->SetYTitle("Counts per event");
        fhNCellsPerStripAcceptEventStrip[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripAcceptEventStrip[icut]);           
  
        fhNCellsPerStripAcceptEventBoth[icut] = new TH1F 
        (Form("hNCellsPerStripAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5);
        fhNCellsPerStripAcceptEventBoth[icut]->SetYTitle("Counts per event");
        fhNCellsPerStripAcceptEventBoth[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripAcceptEventBoth[icut]);         
        
        fhNCellsPerStripEventAccept[icut] = new TH1F 
        (Form("hNCellsPerStripEventAccept_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5);
        fhNCellsPerStripEventAccept[icut]->SetYTitle("Counts per event");
        fhNCellsPerStripEventAccept[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripEventAccept[icut]);        
        
        fhNCellsPerStripNHigh20[icut] = new TH1F 
        (Form("hNCellsPerStripNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5);
        fhNCellsPerStripNHigh20[icut]->SetYTitle("Counts per event");
        fhNCellsPerStripNHigh20[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripNHigh20[icut]);   
        
        // Per SM
        fhNCellsPerStripPerSM[icut] = new TH2F 
        (Form("hNCellsPerStripPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNCellsPerStripPerSM[icut]->SetZTitle("Counts per event");
        fhNCellsPerStripPerSM[icut]->SetYTitle("SM");
        fhNCellsPerStripPerSM[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripPerSM[icut]);   
        
        fhNCellsPerStripPerSMAcceptEventStrip[icut] = new TH2F 
        (Form("hNCellsPerStripPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNCellsPerStripPerSMAcceptEventStrip[icut]->SetZTitle("Counts per event");
        fhNCellsPerStripPerSMAcceptEventStrip[icut]->SetYTitle("SM");
        fhNCellsPerStripPerSMAcceptEventStrip[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripPerSMAcceptEventStrip[icut]);   
  
        fhNCellsPerStripPerSMAcceptEventBoth[icut] = new TH2F 
        (Form("hNCellsPerStripPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNCellsPerStripPerSMAcceptEventBoth[icut]->SetZTitle("Counts per event");
        fhNCellsPerStripPerSMAcceptEventBoth[icut]->SetYTitle("SM");
        fhNCellsPerStripPerSMAcceptEventBoth[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripPerSMAcceptEventBoth[icut]);   
        
        fhNCellsPerStripPerSMEventAccept[icut] = new TH2F 
        (Form("hNCellsPerStripPerSMEventAccept_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNCellsPerStripPerSMEventAccept[icut]->SetZTitle("Counts per event");
        fhNCellsPerStripPerSMEventAccept[icut]->SetYTitle("SM");
        fhNCellsPerStripPerSMEventAccept[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripPerSMEventAccept[icut]);  
        
        fhNCellsPerStripPerSMNHigh20[icut] = new TH2F 
        (Form("hNCellsPerStripPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax),
         48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNCellsPerStripPerSMNHigh20[icut]->SetZTitle("Counts per event");
        fhNCellsPerStripPerSMNHigh20[icut]->SetYTitle("SM");
        fhNCellsPerStripPerSMNHigh20[icut]->SetXTitle("#it{n}_{cells}^{strip}");
        outputContainer->Add(fhNCellsPerStripPerSMNHigh20[icut]);   
        
        // Sum strip cells energy
        //
        fhSumEnCellsPerStrip[icut] = new TH1F 
        (Form("hSumEnCellsPerStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500);
        fhSumEnCellsPerStrip[icut]->SetYTitle("Counts per event");
        fhSumEnCellsPerStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStrip[icut]);   
        
        fhSumEnCellsPerStripAcceptEventStrip[icut] = new TH1F 
        (Form("hSumEnCellsPerStripAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500);
        fhSumEnCellsPerStripAcceptEventStrip[icut]->SetYTitle("Counts per event");
        fhSumEnCellsPerStripAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripAcceptEventStrip[icut]);   
        
        fhSumEnCellsPerStripAcceptEventBoth[icut] = new TH1F 
        (Form("hSumEnCellsPerStripAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500);
        fhSumEnCellsPerStripAcceptEventBoth[icut]->SetYTitle("Counts per event");
        fhSumEnCellsPerStripAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripAcceptEventBoth[icut]);   
        
        fhSumEnCellsPerStripEventAccept[icut] = new TH1F 
        (Form("hSumEnCellsPerStripEventAccept_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500);
        fhSumEnCellsPerStripEventAccept[icut]->SetYTitle("Counts per event");
        fhSumEnCellsPerStripEventAccept[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripEventAccept[icut]);
        
        fhSumEnCellsPerStripNHigh20[icut] = new TH1F 
        (Form("hSumEnCellsPerStripNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV,at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500);
        fhSumEnCellsPerStripNHigh20[icut]->SetYTitle("Counts per event");
        fhSumEnCellsPerStripNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripNHigh20[icut]);   
        
        // Per Strip per SM
        fhSumEnCellsPerStripPerSM[icut] = new TH2F 
        (Form("hSumEnCellsPerStripPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhSumEnCellsPerStripPerSM[icut]->SetZTitle("Counts per event");
        fhSumEnCellsPerStripPerSM[icut]->SetYTitle("SM");
        fhSumEnCellsPerStripPerSM[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripPerSM[icut]);   
        
        fhSumEnCellsPerStripPerSMAcceptEventStrip[icut] = new TH2F 
        (Form("hSumEnCellsPerStripPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetZTitle("Counts per event");
        fhSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetYTitle("SM");
        fhSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripPerSMAcceptEventStrip[icut]);

        fhSumEnCellsPerStripPerSMAcceptEventBoth[icut] = new TH2F 
        (Form("hSumEnCellsPerStripPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetZTitle("Counts per event");
        fhSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetYTitle("SM");
        fhSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripPerSMAcceptEventBoth[icut]);
        
        fhSumEnCellsPerStripPerSMEventAccept[icut] = new TH2F 
        (Form("hSumEnCellsPerStripPerSMEventAccept_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhSumEnCellsPerStripPerSMEventAccept[icut]->SetZTitle("Counts per event");
        fhSumEnCellsPerStripPerSMEventAccept[icut]->SetYTitle("SM");
        fhSumEnCellsPerStripPerSMEventAccept[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripPerSMEventAccept[icut]);
        
        fhSumEnCellsPerStripPerSMNHigh20[icut] = new TH2F 
        (Form("hSumEnCellsPerStripPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV,at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax), 
         1500,0,1500, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhSumEnCellsPerStripPerSMNHigh20[icut]->SetZTitle("Counts per event");
        fhSumEnCellsPerStripPerSMNHigh20[icut]->SetYTitle("SM");
        fhSumEnCellsPerStripPerSMNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhSumEnCellsPerStripPerSMNHigh20[icut]);   
        
        // N vs Sum strip cells energy
        //
        fhNSumEnCellsPerStrip[icut] = new TH2F 
        (Form("hNSumEnCellsPerStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5);
        fhNSumEnCellsPerStrip[icut]->SetZTitle("Counts per event");
        fhNSumEnCellsPerStrip[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStrip[icut]);   
        
        fhNSumEnCellsPerStripAcceptEvent[icut] = new TH2F 
        (Form("hNSumEnCellsPerStripAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5);
        fhNSumEnCellsPerStripAcceptEvent[icut]->SetZTitle("Counts per event");
        fhNSumEnCellsPerStripAcceptEvent[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripAcceptEvent[icut]);   
        
        fhNSumEnCellsPerStripAcceptEventStrip[icut] = new TH2F 
        (Form("hNSumEnCellsPerStripAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5);
        fhNSumEnCellsPerStripAcceptEventStrip[icut]->SetZTitle("Counts per event");
        fhNSumEnCellsPerStripAcceptEventStrip[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripAcceptEventStrip[icut]); 
        
        fhNSumEnCellsPerStripAcceptEventBoth[icut] = new TH2F 
        (Form("hNSumEnCellsPerStripAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5);
        fhNSumEnCellsPerStripAcceptEventBoth[icut]->SetZTitle("Counts per event");
        fhNSumEnCellsPerStripAcceptEventBoth[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripAcceptEventBoth[icut]); 
        
        fhNSumEnCellsPerStripPerSM[icut] = new TH3F 
        (Form("hNSumEnCellsPerStripPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNSumEnCellsPerStripPerSM[icut]->SetZTitle("SM");
        fhNSumEnCellsPerStripPerSM[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripPerSM[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripPerSM[icut]);   
        
        fhNSumEnCellsPerStripPerSMAcceptEvent[icut] = new TH3F 
        (Form("hNSumEnCellsPerStripPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNSumEnCellsPerStripPerSMAcceptEvent[icut]->SetZTitle("SM");
        fhNSumEnCellsPerStripPerSMAcceptEvent[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripPerSMAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripPerSMAcceptEvent[icut]);   
        
        fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut] = new TH3F 
        (Form("hNSumEnCellsPerStripPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetZTitle("SM");
        fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripPerSMAcceptEventStrip[icut]);   
        
        fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut] = new TH3F 
        (Form("hNSumEnCellsPerStripPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
         Form("#Sigma #it{E}_{#it{i}} vs #it{n}_{cells} in strip for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
         150,0,1500, 48,-0.5,47.5, totalSM, fFirstModule-0.5, fLastModule+0.5);
        fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetZTitle("SM");
        fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetYTitle("#it{n}_{cells}^{strip}");
        fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}}^{strip} (GeV)");
        outputContainer->Add(fhNSumEnCellsPerStripPerSMAcceptEventBoth[icut]);   
      }
      
      if ( fFillAllCellEventParamHisto < 1 ) continue;      

      // N cells in event
      //
      fhNCells[icut] = new TH1F 
      (Form("hNCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax),
       17664,0,17664);
      fhNCells[icut]->SetYTitle("Counts per event");
      fhNCells[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCells[icut]);  
      
      fhNCellsNHigh20[icut] = new TH1F 
      (Form("hNCellsNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax),
       17664,0,17664);
      fhNCellsNHigh20[icut]->SetYTitle("Counts per event");
      fhNCellsNHigh20[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsNHigh20[icut]);   
      
      fhNCellsAcceptEvent[icut] = new TH1F 
      (Form("hNCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax),
       17664,0,17664);
      fhNCellsAcceptEvent[icut]->SetYTitle("Counts per event");
      fhNCellsAcceptEvent[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsAcceptEvent[icut]);   
   
      fhNCellsAcceptEventStrip[icut] = new TH1F 
      (Form("hNCellsAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax),
       17664,0,17664);
      fhNCellsAcceptEventStrip[icut]->SetYTitle("Counts per event");
      fhNCellsAcceptEventStrip[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsAcceptEventStrip[icut]);   
      
      fhNCellsAcceptEventBoth[icut] = new TH1F 
      (Form("hNCellsAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax),
       17664,0,17664);
      fhNCellsAcceptEventBoth[icut]->SetYTitle("Counts per event");
      fhNCellsAcceptEventBoth[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsAcceptEventBoth[icut]);   
      
      // Per SM
      fhNCellsPerSM[icut] = new TH2F 
      (Form("hNCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, per SM",fCellEnMins[icut],fCellEnMax),
       1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNCellsPerSM[icut]->SetZTitle("Counts per event");
      fhNCellsPerSM[icut]->SetYTitle("SM");
      fhNCellsPerSM[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerSM[icut]);   
     
      fhNCellsPerSMNHigh20[icut] = new TH2F 
      (Form("hNCellsPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20, per SM",fCellEnMins[icut],fCellEnMax),
       1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNCellsPerSMNHigh20[icut]->SetZTitle("Counts per event");
      fhNCellsPerSMNHigh20[icut]->SetYTitle("SM");
      fhNCellsPerSMNHigh20[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerSMNHigh20[icut]);   
      
      fhNCellsPerSMAcceptEvent[icut] = new TH2F 
      (Form("hNCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3, per SM",fCellEnMins[icut],fCellEnMax),
       1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNCellsPerSMAcceptEvent[icut]->SetZTitle("Counts per event");
      fhNCellsPerSMAcceptEvent[icut]->SetYTitle("SM");
      fhNCellsPerSMAcceptEvent[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerSMAcceptEvent[icut]);   
 
      fhNCellsPerSMAcceptEventStrip[icut] = new TH2F 
      (Form("hNCellsPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3, per SM",fCellEnMins[icut],fCellEnMax),
       1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNCellsPerSMAcceptEventStrip[icut]->SetZTitle("Counts per event");
      fhNCellsPerSMAcceptEventStrip[icut]->SetYTitle("SM");
      fhNCellsPerSMAcceptEventStrip[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerSMAcceptEventStrip[icut]);   
      
      fhNCellsPerSMAcceptEventBoth[icut] = new TH2F 
      (Form("hNCellsPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3, per SM",fCellEnMins[icut],fCellEnMax),
       1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNCellsPerSMAcceptEventBoth[icut]->SetZTitle("Counts per event");
      fhNCellsPerSMAcceptEventBoth[icut]->SetYTitle("SM");
      fhNCellsPerSMAcceptEventBoth[icut]->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerSMAcceptEventBoth[icut]);   
      
      // Sum of cells energy in event
      fhSumEnCells[icut] = new TH1F 
      (Form("hSumEnCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhSumEnCells[icut]->SetYTitle("Counts per event");
      fhSumEnCells[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCells[icut]);   
      
      fhSumEnCellsNHigh20[icut] = new TH1F 
      (Form("hSumEnCellsNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhSumEnCellsNHigh20[icut]->SetYTitle("Counts per event");
      fhSumEnCellsNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsNHigh20[icut]);   
      
      fhSumEnCellsAcceptEvent[icut] = new TH1F 
      (Form("hSumEnCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhSumEnCellsAcceptEvent[icut]->SetYTitle("Counts per event");
      fhSumEnCellsAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsAcceptEvent[icut]);   
 
      fhSumEnCellsAcceptEventStrip[icut] = new TH1F 
      (Form("hSumEnCellsAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhSumEnCellsAcceptEventStrip[icut]->SetYTitle("Counts per event");
      fhSumEnCellsAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsAcceptEventStrip[icut]);   
      
      fhSumEnCellsAcceptEventBoth[icut] = new TH1F 
      (Form("hSumEnCellsAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhSumEnCellsAcceptEventBoth[icut]->SetYTitle("Counts per event");
      fhSumEnCellsAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsAcceptEventBoth[icut]);   
      
      // Sum En vs N cells
      
      fhNSumEnCells[icut] = new TH2F 
      (Form("hNSumEnCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 1104,0,17664);
      fhNSumEnCells[icut]->SetZTitle("Counts per event");
      fhNSumEnCells[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCells[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCells[icut]);   
 
      fhNSumEnCellsAcceptEvent[icut] = new TH2F 
      (Form("hNSumEnCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 1104,0,17664);
      fhNSumEnCellsAcceptEvent[icut]->SetZTitle("Counts per event");
      fhNSumEnCellsAcceptEvent[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsAcceptEvent[icut]);   
      
      fhNSumEnCellsAcceptEventStrip[icut] = new TH2F 
      (Form("hNSumEnCellsAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 1104,0,17664);
      fhNSumEnCellsAcceptEventStrip[icut]->SetZTitle("Counts per event");
      fhNSumEnCellsAcceptEventStrip[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsAcceptEventStrip[icut]);   
      
      fhNSumEnCellsAcceptEventBoth[icut] = new TH2F 
      (Form("hNSumEnCellsAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 1104,0,17664);
      fhNSumEnCellsAcceptEventBoth[icut]->SetZTitle("Counts per event");
      fhNSumEnCellsAcceptEventBoth[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsAcceptEventBoth[icut]); 
      
      // Per SM
      fhSumEnCellsPerSM[icut] = new TH2F 
      (Form("hSumEnCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhSumEnCellsPerSM[icut]->SetZTitle("Counts per event");
      fhSumEnCellsPerSM[icut]->SetYTitle("SM");
      fhSumEnCellsPerSM[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsPerSM[icut]);   
      
      fhSumEnCellsPerSMNHigh20[icut] = new TH2F 
      (Form("hSumEnCellsPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhSumEnCellsPerSMNHigh20[icut]->SetZTitle("Counts per event");
      fhSumEnCellsPerSMNHigh20[icut]->SetYTitle("SM");
      fhSumEnCellsPerSMNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsPerSMNHigh20[icut]);   
  
      fhSumEnCellsPerSMAcceptEvent[icut] = new TH2F 
      (Form("hSumEnCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, low Activity SM3, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhSumEnCellsPerSMAcceptEvent[icut]->SetZTitle("Counts per event");
      fhSumEnCellsPerSMAcceptEvent[icut]->SetYTitle("SM");
      fhSumEnCellsPerSMAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsPerSMAcceptEvent[icut]);   
   
      fhSumEnCellsPerSMAcceptEventStrip[icut] = new TH2F 
      (Form("hSumEnCellsPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, low Activity SM3, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhSumEnCellsPerSMAcceptEventStrip[icut]->SetZTitle("Counts per event");
      fhSumEnCellsPerSMAcceptEventStrip[icut]->SetYTitle("SM");
      fhSumEnCellsPerSMAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsPerSMAcceptEventStrip[icut]);   
      
      fhSumEnCellsPerSMAcceptEventBoth[icut] = new TH2F 
      (Form("hSumEnCellsPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, low Activity SM3, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhSumEnCellsPerSMAcceptEventBoth[icut]->SetZTitle("Counts per event");
      fhSumEnCellsPerSMAcceptEventBoth[icut]->SetYTitle("SM");
      fhSumEnCellsPerSMAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhSumEnCellsPerSMAcceptEventBoth[icut]);   
      
      // Sum En vs N cells
      
      fhNSumEnCellsPerSM[icut] = new TH3F 
      (Form("hNSumEnCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 144,0,1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNSumEnCellsPerSM[icut]->SetZTitle("SM");
      fhNSumEnCellsPerSM[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsPerSM[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsPerSM[icut]);   
      
      fhNSumEnCellsPerSMAcceptEvent[icut] = new TH3F 
      (Form("hNSumEnCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 144,0,1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNSumEnCellsPerSMAcceptEvent[icut]->SetZTitle("SM");
      fhNSumEnCellsPerSMAcceptEvent[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsPerSMAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsPerSMAcceptEvent[icut]);   
      
      fhNSumEnCellsPerSMAcceptEventStrip[icut] = new TH3F 
      (Form("hNSumEnCellsPerSMAcceptEventStrip_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 144,0,1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNSumEnCellsPerSMAcceptEventStrip[icut]->SetZTitle("SM");
      fhNSumEnCellsPerSMAcceptEventStrip[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsPerSMAcceptEventStrip[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsPerSMAcceptEventStrip[icut]);   
      
      fhNSumEnCellsPerSMAcceptEventBoth[icut] = new TH3F 
      (Form("hNSumEnCellsPerSMAcceptEventBoth_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       200,0,2000, 144,0,1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhNSumEnCellsPerSMAcceptEventBoth[icut]->SetZTitle("SM");
      fhNSumEnCellsPerSMAcceptEventBoth[icut]->SetYTitle("#it{n}_{cells}");
      fhNSumEnCellsPerSMAcceptEventBoth[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
      outputContainer->Add(fhNSumEnCellsPerSMAcceptEventBoth[icut]);   
      
      if ( fFillAllCellEventParamHisto < 2 ) continue ;

      // Average energy
      //
      fhAverSumEnCells[icut] = new TH1F 
      (Form("hAverSumEnCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhAverSumEnCells[icut]->SetYTitle("Counts per event");
      fhAverSumEnCells[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCells[icut]);  
      
      fhAverSumEnCellsNHigh20[icut] = new TH1F 
      (Form("hAverSumEnCellsNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhAverSumEnCellsNHigh20[icut]->SetYTitle("Counts per event");
      fhAverSumEnCellsNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCellsNHigh20[icut]);  
      
      fhAverSumEnCellsAcceptEvent[icut] = new TH1F 
      (Form("hAverSumEnCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000);
      fhAverSumEnCellsAcceptEvent[icut]->SetYTitle("Counts per event");
      fhAverSumEnCellsAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCellsAcceptEvent[icut]);  
      
      fhAverSumEnCellsPerSM[icut] = new TH2F 
      (Form("hAverSumEnCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhAverSumEnCellsPerSM[icut]->SetZTitle("Counts per event");
      fhAverSumEnCellsPerSM[icut]->SetYTitle("SM");
      fhAverSumEnCellsPerSM[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCellsPerSM[icut]);  
      
      fhAverSumEnCellsPerSMNHigh20[icut] = new TH2F 
      (Form("hAverSumEnCellsPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhAverSumEnCellsPerSMNHigh20[icut]->SetZTitle("Counts per event");
      fhAverSumEnCellsPerSMNHigh20[icut]->SetYTitle("SM");
      fhAverSumEnCellsPerSMNHigh20[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCellsPerSMNHigh20[icut]);  
      
      fhAverSumEnCellsPerSMAcceptEvent[icut] = new TH2F 
      (Form("hAverSumEnCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, reject high activity events in any SM but SM3, per SM",fCellEnMins[icut],fCellEnMax), 
       10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhAverSumEnCellsPerSMAcceptEvent[icut]->SetZTitle("Counts per event");
      fhAverSumEnCellsPerSMAcceptEvent[icut]->SetYTitle("SM");
      fhAverSumEnCellsPerSMAcceptEvent[icut]->SetXTitle("#Sigma #it{E}_{#it{i}} / #it{n}_{cells} (GeV)");
      outputContainer->Add(fhAverSumEnCellsPerSMAcceptEvent[icut]);  
      
      if ( icut == 0 ) continue ;

      // Fraction of n cells
      //
      fhFracNCells[icut-1] = new TH1F 
      (Form("hFracNCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracNCells[icut-1]->SetYTitle("Counts per event");
      fhFracNCells[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                           fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCells[icut-1]);   
      
      fhFracNCellsNHigh20[icut-1] = new TH1F 
      (Form("hFracNCellsNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracNCellsNHigh20[icut-1]->SetYTitle("Counts per event");
      fhFracNCellsNHigh20[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                                  fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCellsNHigh20[icut-1]); 
      
      fhFracNCellsAcceptEvent[icut-1] = new TH1F 
      (Form("hFracNCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}, reject high activity events in any SM but SM3", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracNCellsAcceptEvent[icut-1]->SetYTitle("Counts per event");
      fhFracNCellsAcceptEvent[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                                      fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCellsAcceptEvent[icut-1]);  
      
      // Per SM

      fhFracNCellsPerSM[icut-1] = new TH2F 
      (Form("hFracNCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracNCellsPerSM[icut-1]->SetZTitle("Counts per event");
      fhFracNCellsPerSM[icut-1]->SetYTitle("SM");
      fhFracNCellsPerSM[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                                fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCellsPerSM[icut-1]);   
      
      fhFracNCellsPerSMNHigh20[icut-1] = new TH2F 
      (Form("hFracNCellsPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracNCellsPerSMNHigh20[icut-1]->SetZTitle("Counts per event");
      fhFracNCellsPerSMNHigh20[icut-1]->SetYTitle("SM");
      fhFracNCellsPerSMNHigh20[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                                       fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCellsPerSMNHigh20[icut-1]);    
      
      fhFracNCellsPerSMAcceptEvent[icut-1] = new TH2F 
      (Form("hFracNCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}, reject high activity events in any SM but SM3, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracNCellsPerSMAcceptEvent[icut-1]->SetZTitle("Counts per event");
      fhFracNCellsPerSMAcceptEvent[icut-1]->SetYTitle("SM");
      fhFracNCellsPerSMAcceptEvent[icut-1]->SetXTitle(Form("#it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f} / #it{n}_{cells}^{#it{E}_{#it{i}}>%1.1f}", 
                                                           fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracNCellsPerSMAcceptEvent[icut-1]);   
      
      // Fraction of summed energy
      //
      fhFracSumEnCells[icut-1] = new TH1F 
      (Form("hFracSumEnCells_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracSumEnCells[icut-1]->SetYTitle("Counts per event");
      fhFracSumEnCells[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
                                               fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracSumEnCells[icut-1]);   
      
      fhFracSumEnCellsNHigh20[icut-1] = new TH1F 
      (Form("hFracSumEnCellsNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20", fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracSumEnCellsNHigh20[icut-1]->SetYTitle("Counts per event");
      fhFracSumEnCellsNHigh20[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
                                                      fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracSumEnCellsNHigh20[icut-1]);   
      
      fhFracSumEnCellsAcceptEvent[icut-1] = new TH1F 
      (Form("hFracSumEnCellsAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}, reject high activity events in any SM but SM3", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005);
      fhFracSumEnCellsAcceptEvent[icut-1]->SetYTitle("Counts per event");
      fhFracSumEnCellsAcceptEvent[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}",fCellEnMins[icut]));
      outputContainer->Add(fhFracSumEnCellsAcceptEvent[icut-1]);   
      
      // Per SM
      
      fhFracSumEnCellsPerSM[icut-1] = new TH2F 
      (Form("hFracSumEnCellsPerSM_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracSumEnCellsPerSM[icut-1]->SetZTitle("Counts per event");
      fhFracSumEnCellsPerSM[icut-1]->SetYTitle("SM");
      fhFracSumEnCellsPerSM[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
                                                    fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracSumEnCellsPerSM[icut-1]);   
      
      fhFracSumEnCellsPerSMNHigh20[icut-1] = new TH2F 
      (Form("hFracSumEnCellsPerSMNHigh20_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}, at least 1 cluster with #it{n}_{cells}^{#it{w}}>20, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracSumEnCellsPerSMNHigh20[icut-1]->SetZTitle("Counts per event");
      fhFracSumEnCellsPerSMNHigh20[icut-1]->SetYTitle("SM");
      fhFracSumEnCellsPerSMNHigh20[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
                                                           fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracSumEnCellsPerSMNHigh20[icut-1]);   
      
      fhFracSumEnCellsPerSMAcceptEvent[icut-1] = new TH2F 
      (Form("hFracSumEnCellsPerSMAcceptEvent_EnMin%1.1f_Max%2.0f",fCellEnMins[icut],fCellEnMax),
       Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}, reject high activity events in any SM but SM3, per SM", 
            fCellEnMins[icut], fCellEnMins[0]), 
       201,0,1.005, totalSM, fFirstModule-0.5, fLastModule+0.5);
      fhFracSumEnCellsPerSMAcceptEvent[icut-1]->SetZTitle("Counts per event");
      fhFracSumEnCellsPerSMAcceptEvent[icut-1]->SetYTitle("SM");
      fhFracSumEnCellsPerSMAcceptEvent[icut-1]->SetXTitle(Form("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f} / #Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>%1.1f}", 
                                                          fCellEnMins[icut], fCellEnMins[0]));
      outputContainer->Add(fhFracSumEnCellsPerSMAcceptEvent[icut-1]);   
      
    }
    
    fhSM3NCellsSumEnSuspiciousEvents = new TH2F 
    (Form("hSM3NCellsSumEnSuspiciousEvents_EnMin%1.1f_Max%2.0f",fCellEnMins[0],fCellEnMax),
     Form("After event rejection, SM3 activity when high energy events in other SM, for %1.1f < #it{E}_{cell} < %2.0f GeV",fCellEnMins[0],fCellEnMax),
     1152,0,1152,100,0,100);
    fhSM3NCellsSumEnSuspiciousEvents->SetZTitle("Counts");
    fhSM3NCellsSumEnSuspiciousEvents->SetYTitle("#it{n}_{cells}^{#it{E}_{#it{i}}>0.5}");
    fhSM3NCellsSumEnSuspiciousEvents->SetXTitle("#Sigma #it{E}_{#it{i}}^{#it{E}_{#it{i}}>0.5}");
    outputContainer->Add(fhSM3NCellsSumEnSuspiciousEvents);   
    
    fhSumEnCellsPerSMEventSuspicious = new TH2F 
    (Form("hSumEnCellsPerSMEventSuspicious_EnMin%1.1f_Max%2.0f",fCellEnMins[0],fCellEnMax),
     Form("#Sigma #it{E}_{#it{i}} for %1.1f < #it{E}_{cell} < %2.0f GeV, per SM, suspicious",fCellEnMins[0],fCellEnMax), 
     10000,0,10000, totalSM, fFirstModule-0.5, fLastModule+0.5);
    fhSumEnCellsPerSMEventSuspicious->SetZTitle("Counts per event");
    fhSumEnCellsPerSMEventSuspicious->SetYTitle("SM");
    fhSumEnCellsPerSMEventSuspicious->SetXTitle("#Sigma #it{E}_{#it{i}} (GeV)");
    outputContainer->Add(fhSumEnCellsPerSMEventSuspicious);
    
    fhNCellsPerSMEventSuspicious = new TH2F 
    (Form("hNCellsPerSMEventSuspicious_EnMin%1.1f_Max%2.0f",fCellEnMins[0],fCellEnMax),
     Form("#it{n}_{cells} for %1.1f < #it{E}_{cell} < %2.0f GeV, per SM",fCellEnMins[0],fCellEnMax),
     1152, 0, 1152, totalSM, fFirstModule-0.5, fLastModule+0.5);
    fhNCellsPerSMEventSuspicious->SetZTitle("Counts per event");
    fhNCellsPerSMEventSuspicious->SetYTitle("SM");
    fhNCellsPerSMEventSuspicious->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerSMEventSuspicious);   
    
    if ( fFillStripHisto )
    {
      fhNStripsPerEventSuspicious = new TH1F 
      ("hNStripsPerEventSuspicious","#it{n}_{strips} active, Low SM3 activity and high energy strips in other SM",450,0,450);
      fhNStripsPerEventSuspicious->SetYTitle("Counts");
      fhNStripsPerEventSuspicious->SetXTitle("#it{n}_{strips}^{high activity}");
      outputContainer->Add(fhNStripsPerEventSuspicious);
      
      fhNStripsPerEventSuspiciousPerSM = new TH2F 
      ("hNStripsPerEventSuspiciousPerSM","#it{n}_{strips} active per SM, Low SM3 activity and high energy strips in other SM",
       24,0,24,totalSM,fFirstModule-0.5, fLastModule+0.5);
      fhNStripsPerEventSuspiciousPerSM->SetZTitle("Counts");
      fhNStripsPerEventSuspiciousPerSM->SetYTitle("SM");
      fhNStripsPerEventSuspiciousPerSM->SetXTitle("#it{n}_{strips}^{high activity}");
      outputContainer->Add(fhNStripsPerEventSuspiciousPerSM);
    }
  } 
  
  if ( fFillStripHisto )
  {
    fhNStripsPerEvent = new TH1F 
    ("hNStripsPerEvent","#it{n}_{strips} active, Low SM3 activity and high energy strips in other SM",450,0,450);
    fhNStripsPerEvent->SetYTitle("Counts");
    fhNStripsPerEvent->SetXTitle("#it{n}_{strips}^{high activity}");
    outputContainer->Add(fhNStripsPerEvent);
    
    fhNStripsPerEventPerSM = new TH2F 
    ("hNStripsPerEventPerSM","#it{n}_{strips} active per SM, Low SM3 activity and high energy strips in other SM",
     24,0,24,totalSM,fFirstModule-0.5, fLastModule+0.5);
    fhNStripsPerEventPerSM->SetZTitle("Counts");
    fhNStripsPerEventPerSM->SetYTitle("SM");
    fhNStripsPerEventPerSM->SetXTitle("#it{n}_{strips}^{high activity}");
    outputContainer->Add(fhNStripsPerEventPerSM);
    
    if ( fFillCellHisto && fFillAllCellEventParamHisto > 0 )
    {
      fhNStripsPerEventAccept = new TH1F 
      ("hNStripsPerEventAccept","#it{n}_{strips} active, Low SM3 activity and high energy strips in other SM, event accepted",450,0,450);
      fhNStripsPerEventAccept->SetYTitle("Counts");
      fhNStripsPerEventAccept->SetXTitle("#it{n}_{strips}^{high activity}");
      outputContainer->Add(fhNStripsPerEventAccept);
      
      fhNStripsPerEventAcceptPerSM= new TH2F 
      ("hNStripsPerEventAcceptPerSM","#it{n}_{strips} active per SM, Low SM3 activity and high energy strips in other SM, event accepted",
       24,0,24,totalSM,fFirstModule-0.5, fLastModule+0.5);
      fhNStripsPerEventAcceptPerSM->SetZTitle("Counts");
      fhNStripsPerEventAcceptPerSM->SetYTitle("SM");
      fhNStripsPerEventAcceptPerSM->SetXTitle("#it{n}_{strips}^{high activity}");
      outputContainer->Add(fhNStripsPerEventAcceptPerSM);
    }
  }
  
  TString evtSelName[] = {"Open","PassLEDSM","PassLEDStrip","PassLEDBoth"};
  if ( fFillClusterHistoAfterEventCut )
  {
    for (Int_t icase = 0; icase < 4; icase++)
    {
      fhEventCutClusterEnergyTime[icase]  = new TH2F 
      (Form("hEventCutClusterEnergyTime_%s",evtSelName[icase].Data()),
       Form("#it{E}_{cluster} vs #it{t}_{cluster} after exotic cuts, %s",evtSelName[icase].Data()),
       eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
       tBinsArray.GetSize() - 1, tBinsArray.GetArray()); 
      fhEventCutClusterEnergyTime[icase]->SetXTitle("#it{E} (GeV)");
      fhEventCutClusterEnergyTime[icase]->SetYTitle("#it{t}_{cluster} (ns)");
      outputContainer->Add(fhEventCutClusterEnergyTime[icase]);
      
      fhEventCutClusterEnergy[icase]  = new TH1F 
      (Form("hEventCutClusterEnergy_%s",evtSelName[icase].Data()),
       Form("#it{E}_{cluster} after time and exotic cuts, %s",evtSelName[icase].Data()),
       eBinsArray.GetSize() - 1, eBinsArray.GetArray()) ;
      fhEventCutClusterEnergy[icase]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEventCutClusterEnergy[icase]);
      
      fhEventCutClusterEnergyECellMax[icase]  = new TH2F 
      (Form("hEventCutClusterEnergyECellMax_%s",evtSelName[icase].Data()),
       Form("#it{E}_{cluster} vs #it{E}^{max}_{cell} after time and exotic cuts, %s",evtSelName[icase].Data()),
       eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
       eBinsArray.GetSize() - 1, eBinsArray.GetArray()); 
      fhEventCutClusterEnergyECellMax[icase]->SetXTitle("#it{E}^{org}_{cluster} (GeV)");
      fhEventCutClusterEnergyECellMax[icase]->SetYTitle("#it{E}^{max}_{cell} (GeV)");
      outputContainer->Add(fhEventCutClusterEnergyECellMax[icase]);
      
      fhEventCutClusterEtaPhiGrid[icase] = new TH2F 
      (Form("hEventCutClusterEtaPhiGrid_%s",evtSelName[icase].Data()),
       Form("colum (#eta) vs row (#varphi) after time and exotic cuts, #it{E}_{cluster}> %2.1f, %s",
            fEMinForExo,evtSelName[icase].Data()),
       colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
       rowBinsArray.GetSize() - 1, rowBinsArray.GetArray());
      fhEventCutClusterEtaPhiGrid[icase]->SetXTitle("column-#eta");
      fhEventCutClusterEtaPhiGrid[icase]->SetYTitle("row-#varphi (rad)");
      outputContainer->Add(fhEventCutClusterEtaPhiGrid[icase]);
      
      fhEventCutBunchCrossing[icase] = new TH1I 
      (Form("hEventCutBunchCrossing_%s"        ,evtSelName[icase].Data()),
       Form("Selected event bunch crossing, %s",evtSelName[icase].Data()),
       3564, -0.5, 3563.5) ;
      fhEventCutBunchCrossing[icase]->SetXTitle("Bunch Crossing");
      outputContainer->Add(fhEventCutBunchCrossing[icase]);
      
      fhEventCutVertexX[icase] = new TH1F 
      (Form("hEventCutVertexX_%s"   ,evtSelName[icase].Data()),
       Form("Selected event vx, %s",evtSelName[icase].Data()),
       200,-1,1) ;
      fhEventCutVertexX[icase]->SetXTitle("#it{v_{x}} (cm)");
      outputContainer->Add(fhEventCutVertexX[icase]);
      
      fhEventCutVertexY[icase] = new TH1F 
      (Form("hEventCutVertexY_%s"   ,evtSelName[icase].Data()),
       Form("Selected event vy, %s",evtSelName[icase].Data()),
       200,-1,1) ;
      fhEventCutVertexY[icase]->SetXTitle("#it{v_{y}} (cm)");
      outputContainer->Add(fhEventCutVertexY[icase]);
      
      fhEventCutVertexZ[icase] = new TH1F 
      (Form("hEventCutVertexZ_%s"   ,evtSelName[icase].Data()),
       Form("Selected event vz, %s",evtSelName[icase].Data()),
       200,-50,50) ;
      fhEventCutVertexZ[icase]->SetXTitle("#it{v_{z}} (cm)");
      outputContainer->Add(fhEventCutVertexZ[icase]);
    }
  }
  
  if ( fFillTrackHistoAfterEventCut )
   {
 
     
     for (Int_t icase = 0; icase < 4; icase++)
     {
       fhEventCutTrackPt[icase]  = new TH1F 
       (Form("hEventCutTrackPt_%s",evtSelName[icase].Data()),
        Form("Track #it{p}_{T}, %s",evtSelName[icase].Data()),
        200,0,20); 
       fhEventCutTrackPt[icase]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       outputContainer->Add(fhEventCutTrackPt[icase]);
       
       fhEventCutTrackEtaPhi[icase]  = new TH2F 
       (Form("hEventCutTrackEtaPhi_%s",evtSelName[icase].Data()),
        Form("Track #it{p}_{T} > 0.5,  #eta vs #varphi, %s",evtSelName[icase].Data()),
        90,-0.9,0.9,90,0,TMath::TwoPi()); 
       fhEventCutTrackEtaPhi[icase]->SetYTitle("#eta");
       fhEventCutTrackEtaPhi[icase]->SetZTitle("#varphi");
       outputContainer->Add(fhEventCutTrackEtaPhi[icase]);
     }
   }
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}


//______________________________
/// Check if the calorimeter setting is ok, if not abort.
//______________________________
void AliAnaCaloExotics::Init()
{
  if(GetCalorimeter() != kPHOS && GetCalorimeter() !=kEMCAL)
    AliFatal(Form("Wrong calorimeter name <%s>", GetCalorimeterString().Data()));
  
  if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
    AliFatal("Analysis of reconstructed data, MC reader not aplicable");
}


//_________________________________________________________
/// Print some relevant parameters set for the analysis.
//_________________________________________________________
void AliAnaCaloExotics::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",GetCalorimeterString().Data());
  printf("Min Amplitude : %2.1f GeV/c\n", fCellAmpMin) ;
  printf("Min Energy for exotic : %2.1f GeV/c\n", fEMinForExo) ;
  printf("Exoticity cut: %0.2f \n", fExoCut) ;
  printf("NCell cut: %d \n", fNCellHighCut) ;
  printf("Time range: [%2.2f,%2.2f] ns\n",fTimeCutMin,fTimeCutMax);
  
  printf("SM : nCell >= %d - Sum E >= %2.0f\n",fHighNCellsCutSM,fHighEnergyCutSM) ;
  printf("SM3: nCell <= %d - Sum E <= %2.0f\n",fLowNCellsCutSM3,fLowEnergyCutSM3) ;
  printf("Strip: nCell > %d-%d - Sum E > %2.0f-%2.0f; Event N Strips <= %d\n",
         fHighNCellsCutStrip[0], fHighNCellsCutStrip[1],
         fHighEnergyCutStrip[0], fHighEnergyCutStrip[1], fEventMaxNumberOfStrips) ;
  printf("SM3 strips: nCell <= %d - Sum E <= %2.0f\n",fLowNCellsCutSM3Strip,fLowEnergyCutSM3Strip) ;

  printf("Min Cell Energy cut: ");
  for(Int_t i = 0; i < fgkNCellEnMinBins; i++) printf("%d) E %1.2f; ", i,fCellEnMins[i] );
  printf("\n");
  
  printf("Fill cell histo: %d\n"          , fFillCellHisto) ;
  printf("Fill all cell event histo: %d\n", fFillAllCellEventParamHisto) ;
  printf("Fill strip histo : %d\n"        , fFillStripHisto) ;
  printf("Fill 1 cell cluster histo: %d\n", fFill1CellHisto) ;
  printf("Fill Matching histo: %d\n"      , fFillMatchingHisto) ;
  printf("Fill Exoticity 50ns cut: %d\n"  , fFillExo50ns) ;
  printf("Fill histo after event cut: cluster %d - track %d\n", 
         fFillClusterHistoAfterEventCut, fFillTrackHistoAfterEventCut) ;
}


//_____________________________________________________
/// Main task method, call all the methods filling QA histograms.
//_____________________________________________________
void  AliAnaCaloExotics::MakeAnalysisFillHistograms()
{
  AliDebug(1,"Start");

  // Get List with CaloClusters , calo Cells, init min amplitude
  TObjArray     * caloClusters = NULL;
  AliVCaloCells * cells        = 0x0;
  
  if      (GetCalorimeter() == kPHOS)
  {
    caloClusters = GetPHOSClusters();
    cells        = GetPHOSCells();
  }
  else if (GetCalorimeter() == kEMCAL)
  {
    caloClusters = GetEMCALClusters();
    cells        = GetEMCALCells();
  }
  
  if( !caloClusters || !cells )
  {
    AliWarning(Form("AliAnaCaloExotics::MakeAnalysisFillHistograms() - No CaloClusters or CaloCells available"));
    return; 
  }

  if(caloClusters->GetEntriesFast() == 0) return ;
  
  //printf("Exotic Task: N cells %d, N clusters %d \n",cells->GetNumberOfCells(),caloClusters->GetEntriesFast());
  
  // Clusters
  ClusterHistograms(caloClusters,cells);
  
  // Strips of cells 2x24
  if ( fFillStripHisto )  StripHistograms(cells);
  
  // Cells  
  if ( fFillCellHisto )  CellHistograms(cells);
  
  if ( fFillClusterHistoAfterEventCut )
    ClusterHistogramsAfterEventCut(caloClusters,cells);
  
  if ( fFillTrackHistoAfterEventCut )
    TrackHistogramsAfterEventCut();
  
  AliDebug(1,"End");
}

