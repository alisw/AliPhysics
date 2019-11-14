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
fLED20(0),                             fLED12(0),   
fLED20Time(0),                         fLED12Time(0),
fFillCellHisto(0),                     fFill1CellHisto(0),                    
fFillMatchingHisto(0),                 fFillSameDiffFracHisto(0),
fConstantTimeShift(0),                 fClusterMomentum(),         

// Histograms
fhNClusterPerEventNCellHigh20(0),      fhNClusterPerEventNCellHigh12(0),        
fhNClusterPerEventExotic(0),           
fhNClusterPerEventExotic1Cell(0),      fhNClusterPerEventExoticNCell(0),

fh2NClusterPerEventNCellHigh20(0),     fh2NClusterPerEventNCellHigh12(0),        
fh2NClusterPerEventExotic(0),           
fh2NClusterPerEventExotic1Cell(0),     fh2NClusterPerEventExoticNCell(0),
fh2NClusterPerEventExoticAmpMax(0),

fhExoticityECellMinCut(0),             
fhExoticityEClus(0),                   fhExoticityWEClus(0),
fhExoticityEClusAllSameTCard(0),
fhExoticityEClusPerSM(0),              fhExoticityEMaxCell(0),
fhExoticityEClusTrackMatch(0),         fhExoticity1Cell(0),

fhNCellsPerCluster(0),                 fhNCellsPerClusterW(0),  
fhNCellsPerClusterEMaxCell(0),         fhNCellsPerClusterWEMaxCell(0),   
fhNCellsPerClusterOpenTime(0),         fhNCellsPerClusterWOpenTime(0),
fhNCellsPerClusterEMaxCellOpenTime(0), fhNCellsPerClusterWEMaxCellOpenTime(0),  
fhNCellsPerClusterAllSameTCard(0),     
fhNCellsPerClusterExo(0),              fhNCellsPerClusterExoW(0),
fhNCellsPerClusterPerSM(0),            fhNCellsPerClusterWPerSM(0),              
fhNCellsPerClusterTrackMatch(0),       fhNCellsPerClusterExoTrackMatch(0),    
fhNCellsPerClusterM02(0),

fhEtaPhiGridExoEnCut(0),               fhEtaPhiGridExoEnCutSameFracCut(0),               
fhEtaPhiGridEnExoCut(0),               fhEtaPhiGridEn1Cell(0),
fhEtaPhiGridEnHighNCells(0),           fhEtaPhiGridNCellEnCut(0),

fhTimeEnergyExo(0),                    fhTimeEnergyExoW(0),                    
fhTimeEnergy1Cell(0),                
fhTimeDiffClusCellExo(0),              fhTimeDiffAmpClusCellExo(0),   
fhTimeEnergyM02(0),                    fhTimeDiffClusCellM02(0),  
fhTimeEnergyNCells(0),                 fhTimeEnergyNCellsW(0),                 fhTimeNCellCut(0),                

fhM02EnergyAllSameTCard(0),      
fhM02EnergyExo(0),                     fhM02EnergyExoW(0),           
fhM20EnergyExoM02MinCut(0),                                                 

// Other ncell
fhNCellsPerClusterSame(0),             fhNCellsPerClusterDiff(0),
fhNCellsPerClusterSame5(0),            fhNCellsPerClusterDiff5(0),
fhNCellsPerClusterSameW (0),           fhNCellsPerClusterDiffW (0),  
fhNCellsPerClusterSameDiff(0),         
fhNCellsPerClusterSameFrac(0),         fhNCellsPerClusterSameFracExo(0),

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
fhCellGridTime05(0),                    fhCellGridTime1(0),
fhCellGridTime2(0),                     fhCellGridTime5(0),
fhCellTimeNCell20(0),                   fhCellTimeNCell12(0),

fhSumEnCells05(0),                      fhSumEnCells1(0), 
fhSumEnCells2(0),                       fhSumEnCells5(0),
fhNCells05(0),                          fhNCells1(0), 
fhNCells2(0),                           fhNCells5(0),
fhAverSumEnCells05(0),                  fhAverSumEnCells1(0), 
fhAverSumEnCells2(0),                   fhAverSumEnCells5(0),
fhFracNCells1(0),                       fhFracNCells2(0),
fhFracNCells5(0),
fhFracSumEnCells1(0),                   fhFracSumEnCells2(0),                   
fhFracSumEnCells5(0),
fhSumEnCells05NHigh20(0),               fhSumEnCells1NHigh20(0), 
fhSumEnCells2NHigh20(0),                fhSumEnCells5NHigh20(0),
fhNCells05NHigh20(0),                   fhNCells1NHigh20(0), 
fhNCells2NHigh20(0),                    fhNCells5NHigh20(0),
fhAverSumEnCells05NHigh20(0),           fhAverSumEnCells1NHigh20(0), 
fhAverSumEnCells2NHigh20(0),            fhAverSumEnCells5NHigh20(0),
fhFracNCells1NHigh20(0),                fhFracNCells2NHigh20(0),
fhFracNCells5NHigh20(0),
fhFracSumEnCells1NHigh20(0),            fhFracSumEnCells2NHigh20(0),                   
fhFracSumEnCells5NHigh20(0)
{        
  AddToHistogramsName("AnaCaloExotic_");
  
  fCellAmpMin = 0.5;
  fEMinForExo = 10.0;
  fExoCut     = 0.97;
  fNCellHighCut = 20;
  
  fEnergyBins [0] =   7; fEnergyBins [1] =  10;  fEnergyBins [2] =  16;
  fEnergyBins [3] =  22; fEnergyBins [4] =  30;  fEnergyBins [5] =  50;
  fEnergyBins [6] =  75; fEnergyBins [7] = 100;  fEnergyBins [8] = 125;
  fEnergyBins [9] = 150; fEnergyBins[10] = 175;  fEnergyBins[11] = 200;
  
  for(Int_t i = 0; i < fgkNEBins; i++) 
  {
    fhM02ExoNCells        [i] = 0;
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
  
//  fCellMinEnBins[0] = 0.10; fCellMinEnBins[1] = 0.25; fCellMinEnBins[2] = 0.50;
//  fCellMinEnBins[3] = 0.75; fCellMinEnBins[4] = 1.00; fCellMinEnBins[5] = 1.25;
//
//  for(Int_t i = 0; i < fgkNCellMinEnBins; i++) 
//  {
//    fhExoticityECellMinCut[i] = 0;
//  }
  
  for(Int_t i = 0; i < 20; i++)  fhNCellsPerClusterExoPerSM[i] = 0;
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
  
  Int_t    bc     = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  Int_t nCells05 = 0;
  Int_t nCells1  = 0;
  Int_t nCells2  = 0;
  Int_t nCells5  = 0;

  Float_t eCells05 = 0;
  Float_t eCells1  = 0;
  Float_t eCells2  = 0;
  Float_t eCells5  = 0;
  
  for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++)
  {
    if ( cells->GetCellNumber(iCell) < 0 ||  cells->GetAmplitude(iCell) < fCellAmpMin ) continue; // CPV 
    
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
    
    // Fill histograms
    
    fhCellExoAmp    ->Fill(amp       , exoticity, GetEventWeight());
    fhCellExoAmpTime->Fill(amp , time, exoticity, GetEventWeight());
    
    if ( amp > fEMinForExo )
      fhCellExoGrid->Fill(icolAbs, irowAbs, exoticity, GetEventWeight());

    // If there is at least one cluster with nCellW>20 or 12, check the cell distribution
    if ( fLED20 )
    {
      fhCellGridTimeHighNCell20->Fill(icolAbs, irowAbs, fLED20Time, GetEventWeight());
      fhCellTimeNCell20->Fill(fLED20Time,time);
    }
  
    if ( fLED12 )
    {
      fhCellGridTimeHighNCell12->Fill(icolAbs, irowAbs, fLED12Time, GetEventWeight());
      fhCellTimeNCell12->Fill(fLED12Time,time);
    }
    
    if ( amp > 0.5 )
    {
      fhCellGridTime05->Fill(icolAbs, irowAbs, time, GetEventWeight());
      nCells05++;
      eCells05+=amp;
    }
    
    if ( amp > 1 )
    {
      fhCellGridTime1->Fill(icolAbs, irowAbs, time, GetEventWeight());
      nCells1++;
      eCells1+=amp;
    }
    
    if ( amp > 2 )
    {
      fhCellGridTime2->Fill(icolAbs, irowAbs, time, GetEventWeight());
      nCells2++;
      eCells2+=amp;
    }
    
    if ( amp > 5 )
    {
      fhCellGridTime5->Fill(icolAbs, irowAbs, time, GetEventWeight());
      nCells5++;
      eCells5+=amp;
    }
  } // Cell loop
  
  Float_t averECells05 = 0;
  Float_t averECells1  = 0;
  Float_t averECells2  = 0;
  Float_t averECells5  = 0;
  
  if ( nCells05 > 0 ) averECells05 = eCells05 / nCells05;
  if ( nCells1  > 0 ) averECells1  = eCells1  / nCells1 ;
  if ( nCells2  > 0 ) averECells2  = eCells2  / nCells2 ;
  if ( nCells5  > 0 ) averECells5  = eCells5  / nCells5 ;
  
  Float_t frNCells1 = 0;
  Float_t frNCells2 = 0;
  Float_t frNCells5 = 0;
  if ( nCells05 > 0 )
  {
    frNCells1 = (1.*nCells1) / nCells05;
    frNCells2 = (1.*nCells2) / nCells05;
    frNCells5 = (1.*nCells5) / nCells05;
  }
  
  Float_t frEn1 = 0;
  Float_t frEn2 = 0;
  Float_t frEn5 = 0;
  if ( eCells05 > 0 )
  {
    frEn1 = eCells1 / eCells05;
    frEn2 = eCells2 / eCells05;
    frEn5 = eCells5 / eCells05;
  }
  
  fhSumEnCells05->Fill(eCells05, GetEventWeight());
  fhSumEnCells1 ->Fill(eCells1 , GetEventWeight());
  fhSumEnCells2 ->Fill(eCells2 , GetEventWeight());
  fhSumEnCells5 ->Fill(eCells5 , GetEventWeight());
  
  fhNCells05    ->Fill(nCells05, GetEventWeight());
  fhNCells1     ->Fill(nCells1 , GetEventWeight());
  fhNCells2     ->Fill(nCells2 , GetEventWeight());
  fhNCells5     ->Fill(nCells5 , GetEventWeight());
  
  fhAverSumEnCells05->Fill(averECells05, GetEventWeight());
  fhAverSumEnCells1 ->Fill(averECells1 , GetEventWeight());
  fhAverSumEnCells2 ->Fill(averECells2 , GetEventWeight());
  fhAverSumEnCells5 ->Fill(averECells5 , GetEventWeight());
  
  fhFracNCells1->Fill(frNCells1, GetEventWeight());
  fhFracNCells2->Fill(frNCells2, GetEventWeight());
  fhFracNCells5->Fill(frNCells5, GetEventWeight());
  
  fhFracSumEnCells1 ->Fill(frEn1, GetEventWeight());
  fhFracSumEnCells2 ->Fill(frEn2, GetEventWeight());
  fhFracSumEnCells5 ->Fill(frEn5, GetEventWeight());

  if ( fLED20 )  
  {
    fhSumEnCells05NHigh20->Fill(eCells05, GetEventWeight());
    fhSumEnCells1NHigh20 ->Fill(eCells1 , GetEventWeight());
    fhSumEnCells2NHigh20 ->Fill(eCells2 , GetEventWeight());
    fhSumEnCells5NHigh20 ->Fill(eCells5 , GetEventWeight());
    
    fhNCells05NHigh20    ->Fill(nCells05, GetEventWeight());
    fhNCells1NHigh20     ->Fill(nCells1 , GetEventWeight());
    fhNCells2NHigh20     ->Fill(nCells2 , GetEventWeight());
    fhNCells5NHigh20     ->Fill(nCells5 , GetEventWeight());
    
    fhAverSumEnCells05NHigh20->Fill(averECells05, GetEventWeight());
    fhAverSumEnCells1NHigh20 ->Fill(averECells1 , GetEventWeight());
    fhAverSumEnCells2NHigh20 ->Fill(averECells2 , GetEventWeight());
    fhAverSumEnCells5NHigh20 ->Fill(averECells5 , GetEventWeight());
    
    fhFracNCells1NHigh20->Fill(frNCells1, GetEventWeight());
    fhFracNCells2NHigh20->Fill(frNCells2, GetEventWeight());
    fhFracNCells5NHigh20->Fill(frNCells5, GetEventWeight());
    
    fhFracSumEnCells1NHigh20 ->Fill(frEn1, GetEventWeight());
    fhFracSumEnCells2NHigh20 ->Fill(frEn2, GetEventWeight());
    fhFracSumEnCells5NHigh20 ->Fill(frEn5, GetEventWeight()); 
  }
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
    Bool_t matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(), GetReader()->GetInputEvent());
     
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
    Float_t exoticityW = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc,0,kTRUE,en)/ampMax;

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
    Int_t   nCellW      = 0, nCellSameW  = 0,  nCellDiffW  = 0;
    Int_t   rowDiff     = -100, colDiff = -100;
    Float_t enSame      = 0,    enDiff  = 0;
    Float_t enSame5     = 0,    enDiff5 = 0;
    Float_t enSameW     = 0,    enDiffW = 0;

    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId     = clus->GetCellsAbsId()[ipos];  
      
      Float_t amp       = cells->GetCellAmplitude(absId);
      
      if ( absId == absIdMax || amp < 0.1 ) continue;
      
      Bool_t  sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
      
      if ( sameTCard ) 
      { 
        nCellSame ++; 
        enSame += amp; 
        if ( TMath::Abs(rowDiff) <= 1 && TMath::Abs(colDiff) <= 1 ) 
        {
          enSame5 += amp;
          nCellSame5++;
        }        
      }
      else             
      { 
        nCellDiff ++; 
        enDiff += amp; 
        if ( TMath::Abs(rowDiff) <= 1 && TMath::Abs(colDiff) <= 1 ) 
        {
          enDiff5 += amp;
          nCellDiff5++;
        }
      }
      
      Float_t weight    = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(amp, en);
      
      if( weight > 0.01 ) 
      {
        nCellW++;
        if ( sameTCard ) { nCellSameW++; enSameW+=amp; }
          else           { nCellDiffW++; enDiffW+=amp; }
      }
    } // cells in cluster loop
    
    Float_t fracEnDiffSame       = 0, fracNCellDiffSame    = 0, fracEnNCellDiffSame  = 0;
    Float_t fracEnDiffSameW      = 0, fracNCellDiffSameW   = 0, fracEnNCellDiffSameW = 0;
    Float_t fracEnDiffSame5      = 0, fracNCellDiffSame5   = 0, fracEnNCellDiffSame5 = 0;
    
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
      if ( exoticity > 0.97 ) nExotic++;
      if ( nCaloCellsPerCluster == 1 ) nExotic1Cell++;
      if ( exoticity > 0.97 && nCaloCellsPerCluster > 1 ) nExoticNCell++;
    }
    
    if ( exoticity > 0.97 && ampMax > highestAmpMaxExo ) 
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
      fhTimeEnergyExoW ->Fill(en, tmax, exoticityW, GetEventWeight());
      fhTimeEnergyM02  ->Fill(en, tmax, m02       , GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhTimeEnergy1Cell->Fill(en, tmax,            GetEventWeight());
    
//    if(maxCellFraction < 0.2) printf("frac %2.2f, en cell %2.2f, en cluster %2.2f, refrac %2.2f\n",
//                                     maxCellFraction,ampMax,en,ampMax/en);
    
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
    
    // Apply time cut for other parameters
    //
    if ( tmax < fTimeCutMin || tmax> fTimeCutMax ) continue;

    // Exoticity
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      fhExoticityEClus   ->Fill(en    , exoticity , GetEventWeight());
      fhExoticityWEClus  ->Fill(en    , exoticityW, GetEventWeight());
      fhExoticityEMaxCell->Fill(ampMax, exoticity , GetEventWeight());
      
      fhExoticityEClusPerSM ->Fill(en, exoticity, nSM, GetEventWeight());
      
      if ( matched && fFillMatchingHisto )
        fhExoticityEClusTrackMatch->Fill(en, exoticity, GetEventWeight());
      
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
      
    }
    else if ( fFill1CellHisto )
      fhExoticity1Cell->Fill(en, exoticity, GetEventWeight());
    
    // N cells per cluster
    //
    fhNCellsPerCluster               ->Fill(en, nCaloCellsPerCluster            , GetEventWeight());
    fhNCellsPerClusterW              ->Fill(en, nCellW                          , GetEventWeight());
    fhNCellsPerClusterExo            ->Fill(en, nCaloCellsPerCluster, exoticity , GetEventWeight());
    fhNCellsPerClusterExoW           ->Fill(en, nCaloCellsPerCluster, exoticityW, GetEventWeight());
    fhNCellsPerClusterM02            ->Fill(en, nCaloCellsPerCluster, m02       , GetEventWeight());

    fhNCellsPerClusterEMaxCell       ->Fill(ampMax, nCaloCellsPerCluster       , GetEventWeight());
    fhNCellsPerClusterWEMaxCell      ->Fill(ampMax, nCellW                     , GetEventWeight());
    
    fhNCellsPerClusterPerSM          ->Fill(en, nCaloCellsPerCluster, nSM      , GetEventWeight());
    fhNCellsPerClusterWPerSM         ->Fill(en, nCellW              , nSM      , GetEventWeight());
    fhNCellsPerClusterExoPerSM[nSM]  ->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
 
    if ( matched && fFillMatchingHisto )
    {
      fhNCellsPerClusterTrackMatch   ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
      fhNCellsPerClusterExoTrackMatch->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    }
    
    // Acceptance
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      if ( en > fEMinForExo )
      {
        fhEtaPhiGridExoEnCut  ->Fill(icolMaxAbs, irowMaxAbs, exoticity, GetEventWeight());
        fhEtaPhiGridNCellEnCut->Fill(icolMaxAbs, irowMaxAbs, nCaloCellsPerCluster, GetEventWeight());
      }
      
      if ( exoticity > fExoCut )
        fhEtaPhiGridEnExoCut  ->Fill(icolMaxAbs, irowMaxAbs, en       , GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhEtaPhiGridEn1Cell     ->Fill(icolMaxAbs, irowMaxAbs, en       , GetEventWeight());
    
    if ( nCellW > fNCellHighCut )
      fhEtaPhiGridEnHighNCells->Fill(icolMaxAbs, irowMaxAbs, en       , GetEventWeight());

    if ( nCaloCellsPerCluster > 1 )
    {
      fhCellMaxClusterEn           ->Fill(enOrg, ampMax                       , GetEventWeight());
      fhCellMaxClusterEnRatio      ->Fill(en    , 1-maxCellFraction           , GetEventWeight());
      fhCellMaxClusterEnRatioNCellW->Fill(en    , 1-maxCellFraction, nCellW   , GetEventWeight());
      fhCellMaxClusterEnRatioExo   ->Fill(en    , 1-maxCellFraction, exoticity, GetEventWeight());
    }
    
    // Cell cluster loop
    //
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId     = clus->GetCellsAbsId()[ipos];  
     
      Float_t amp       = cells->GetCellAmplitude(absId);

      if ( absId == absIdMax || amp < 0.1 ) continue;

      Bool_t  sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
      
      if ( ebin >= 0 && ebin < fgkNEBins-1 )
      {
        fhClusterColRowExo[icolMax%2][ebin]->Fill(colDiff, rowDiff, exoticity, GetEventWeight());
        //if(exoticity < fExoCut)
        //  fhClusterColRow   [icolMax%2][ebin]->Fill(colDiff, rowDiff, GetEventWeight());
        if ( nCellW > fNCellHighCut )
          fhClusterColRowPerSMHighNCell[ebin]->Fill(colDiff, rowDiff, nSM, GetEventWeight());
      }
            
      if ( sameTCard ) 
      { 
        fhCellEnSameExo->Fill(en, amp, exoticity, GetEventWeight());
      }
      else             
      { 
        fhCellEnDiffExo->Fill(en, amp, exoticity, GetEventWeight());
      }

      fhCellEnNCellW    ->Fill(en    , amp, nCellW, GetEventWeight());
      fhCellEnNCellWEMax->Fill(ampMax, amp, nCellW, GetEventWeight());

//      Float_t weight    = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(amp, en);
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
      
      fhCellTimeDiffNCellW    ->Fill(en    , tdiff, nCellW, GetEventWeight());
      fhCellTimeDiffNCellWEMax->Fill(ampMax, tdiff, nCellW, GetEventWeight());

      
      if ( en > fEMinForExo )
        fhTimeDiffAmpClusCellExo->Fill(amp, tdiff, exoticity, GetEventWeight());
      
    } // Fill cell-cluster histogram loop
      
//    if ( en > 10 ) 
//    {
//      printf("total %d, same %d, diff %d, same5 %d, diff3 %d, same5+diff3 %d\n",
//             nCaloCellsPerCluster,nCellSame,nCellDiff,nCellSame5,nCellDiff5,nCellSame5+nCellDiff5);
//      
//      printf("E %2.2f, Esame %2.2f, Ediff %2.2f, Esame5 %2.2f, Ediff3 %2.2f\n",
//             en,enSame,enDiff,enSame5,enDiff5);
//    }
    
    fhNCellsPerClusterSameDiff->Fill(en, nCellSame, nCellDiff, GetEventWeight());

    if( fFillSameDiffFracHisto )
    {
      fhNCellsPerClusterSame    ->Fill(en, nCellSame  , GetEventWeight());
      fhNCellsPerClusterDiff    ->Fill(en, nCellDiff  , GetEventWeight());
      fhNCellsPerClusterSame5   ->Fill(en, nCellSame5 , GetEventWeight());
      fhNCellsPerClusterDiff5   ->Fill(en, nCellDiff5 , GetEventWeight());
      fhNCellsPerClusterSameW   ->Fill(en, nCellSameW , GetEventWeight());
      fhNCellsPerClusterDiffW   ->Fill(en, nCellDiffW , GetEventWeight());    
    }
    
    if ( nCaloCellsPerCluster > 1 && fFillSameDiffFracHisto )
    {
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
      
      Float_t frac = (1.*nCellSame)/(nCaloCellsPerCluster-1.);
//      if ( frac > 0.9 ) 
//        printf("E %2.2f, n cells %d, n in T-Card / n cluster = %1.2f\n",
//               en, nCaloCellsPerCluster, frac);
      fhNCellsPerClusterSameFrac   ->Fill(en, frac,            GetEventWeight());
      fhNCellsPerClusterSameFracExo->Fill(en, frac, exoticity, GetEventWeight());
      
      if ( frac > 0.99 )
      {
        fhNCellsPerClusterAllSameTCard->Fill(en, nCaloCellsPerCluster, GetEventWeight());
        fhExoticityEClusAllSameTCard  ->Fill(en, exoticity           , GetEventWeight());
        fhM02EnergyAllSameTCard       ->Fill(en, m02                 , GetEventWeight());
        if ( en > fEMinForExo && fFillSameDiffFracHisto )
          fhEtaPhiGridExoEnCutSameFracCut->Fill(icolMaxAbs, irowMaxAbs, exoticity, GetEventWeight());
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
        fhNCellsSameDiffExo[ebin]->Fill(nCellSame, nCellDiff, exoticity, GetEventWeight());
        fhEnSameDiffExo    [ebin]->Fill(enSame   , enDiff   , exoticity, GetEventWeight());
        if ( nCellSame > 0 && nCellDiff > 0 )
          fhEnNCellsSameDiffExo[ebin]->Fill(enSame/nCellSame, enDiff/nCellDiff, exoticity, GetEventWeight());
      }
    }
    
    // Shower shape
    //
    if ( nCaloCellsPerCluster > 1 )
    {      
      fhM02EnergyExo ->Fill(en, m02, exoticity , GetEventWeight());
      fhM02EnergyExoW->Fill(en, m02, exoticityW, GetEventWeight());
            
      if ( m02 > 0.1 )
        fhM20EnergyExoM02MinCut->Fill(en, m20, exoticity, GetEventWeight());
      
      if ( ebin >= 0 && ebin < fgkNEBins-1 )
        fhM02ExoNCells[ebin]->Fill(m20, exoticity, nCaloCellsPerCluster, GetEventWeight()); ;
    }

    // Track matching
    //
    if ( matched && fFillMatchingHisto )
    {
      AliVTrack *track = GetCaloUtils()->GetMatchedTrack(clus, GetReader()->GetInputEvent());
      
      if(!track) continue;
      
      Double_t tmom   = track->P();
      Double_t eOverP = en/tmom;

      Bool_t positive = kFALSE;
      if(track) positive = (track->Charge()>0);
      
      // Residuals
      Float_t deta  = clus->GetTrackDz();
      Float_t dphi  = clus->GetTrackDx();
      
      if ( nCaloCellsPerCluster > 1 )
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
      
      } // more than 1 cell
      else if ( fFill1CellHisto )
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
      } // 1 cell clusters
      
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
  
  snprintf(onePar,buffersize,"%2.0f < time < %2.0f ns;",fTimeCutMin,fTimeCutMax) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"time shift = %2.0f;",fConstantTimeShift) ;
  parList+=onePar ;
  
  snprintf(onePar,buffersize,"fill cell histo: %d;",fFillCellHisto) ;
  parList+=onePar ;

  snprintf(onePar,buffersize,"fill cluster with 1 cell histo: %d;",fFill1CellHisto) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"fill cluster track-matching histo: %d;",fFillMatchingHisto) ;
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
  eBinning.AddStep(15,0.5);                 // 30
  if(ptmax > 15)  eBinning.AddStep(30,1);   // 15
  if(ptmax > 30)  eBinning.AddStep(60,2.5); // 12
  if(ptmax > 60)  eBinning.AddStep(100,5);  // 8 
  if(ptmax > 100) eBinning.AddStep(200,10); // 10
  if(ptmax > 200) eBinning.AddStep(300,20); // 5
  TArrayD eBinsArray;
  eBinning.CreateBinEdges(eBinsArray);
 
  TCustomBinning e2Binning;
  e2Binning.SetMinimum(0.1);
  e2Binning.AddStep(2,0.1);  // 19
  e2Binning.AddStep(10,0.2); // 40
  e2Binning.AddStep(20,0.5); // 20
  e2Binning.AddStep(30,1.0); // 10
  e2Binning.AddStep(50,2.0); // 20
  e2Binning.AddStep(100,5.); // 25
  TArrayD e2BinsArray;
  e2Binning.CreateBinEdges(e2BinsArray);
  
  TCustomBinning eminBinning;
  eminBinning.SetMinimum(0.075);
  eminBinning.AddStep(2.075,0.05); 
  TArrayD eminBinsArray;
  eminBinning.CreateBinEdges(eminBinsArray);
  
  // Exoticity
//Int_t nexobins  = 201; Float_t exomin  = -1 ; Float_t exomax  = 1.01;

  TCustomBinning fBinning;
  fBinning.SetMinimum(-1.0);
  fBinning.AddStep(0.0 ,0.100); // 10
  fBinning.AddStep(0.50,0.050); // 25
  fBinning.AddStep(0.90,0.020); // 20 
  fBinning.AddStep(1.01,0.002); // 50 
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
   "#it{n}_{cluster} per event with #it{n}_{cell}^{w} > 20",
   1000,0,1000); 
  fhNClusterPerEventNCellHigh20->SetXTitle("#it{n}_{cluster}^{#it{n}_{cell}^{w} > 20}");
  fhNClusterPerEventNCellHigh20->SetYTitle("Counts per event");
  outputContainer->Add(fhNClusterPerEventNCellHigh20);
 
  fhNClusterPerEventNCellHigh12 = new TH1F 
  ("hNClusterPerEventNCellHigh12",
   "#it{n}_{cluster} per event with #it{n}_{cell}^{w} > 12",
   1000,0,1000); 
  fhNClusterPerEventNCellHigh12->SetXTitle("#it{n}_{cluster}^{#it{n}_{cell}^{w} > 12}");
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
   "#it{n}_{cluster} per event with #it{n}_{cell}^{w} > 20",
   500,0,500,500,0,500); 
  fh2NClusterPerEventNCellHigh20->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventNCellHigh20->SetYTitle("#it{n}_{cluster}^{#it{n}_{cell}^{w} > 20}");
  fh2NClusterPerEventNCellHigh20->SetZTitle("Counts per event");
  outputContainer->Add(fh2NClusterPerEventNCellHigh20);
  
  fh2NClusterPerEventNCellHigh12 = new TH2F 
  ("h2NClusterPerEventNCellHigh12",
   "#it{n}_{cluster} per event, total vs with #it{n}_{cell}^{w} > 12",
   500,0,500,500,0,500); 
  fh2NClusterPerEventNCellHigh12->SetXTitle("#it{n}_{cluster}");
  fh2NClusterPerEventNCellHigh12->SetYTitle("#it{n}_{cluster}^{#it{n}_{cell}^{w} > 12}");
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
  
  fhExoticityEClus = new TH2F 
  ("hExoticityEClus","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhExoticityEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClus->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClus);    

  fhExoticityWEClus = new TH2F 
  ("hExoticityWEClus","cell #it{F}_{+}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhExoticityWEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityWEClus->SetYTitle("#it{F}_{+}^{#it{w}}");
  outputContainer->Add(fhExoticityWEClus);    
  
  fhExoticityEMaxCell = new TH2F 
  ("hExoticityEMaxCell","cell #it{F}_{+} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
   eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhExoticityEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhExoticityEMaxCell->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEMaxCell);    
  
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
  
  // N cells per cluster
  //
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
  
  
  // Different n cells definitions
  //
  fhNCellsPerClusterSameDiff  = new TH3F 
  ("hNCellsPerClusterSameDiff","# cells per cluster in same vs different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   //nptbins,ptmin,ptmax, 17,0,17,nceclbins,nceclmin,nceclmax); 
       eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
   nsameBinsArray.GetSize() - 1, nsameBinsArray.GetArray(), 
       nBinsArray.GetSize() - 1,     nBinsArray.GetArray());
  fhNCellsPerClusterSameDiff->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameDiff->SetYTitle("#it{n}_{cells, same T-Card}");
  fhNCellsPerClusterSameDiff->SetZTitle("#it{n}_{cells, diff T-Card}");
  outputContainer->Add(fhNCellsPerClusterSameDiff);
  
  if ( fFillSameDiffFracHisto )
  {
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
    ("hNCellsPerClusterSameW ","# cells per cluster with #it{w} in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, 17,0,17); 
        eBinsArray.GetSize() - 1,      eBinsArray.GetArray(), 
    nsameBinsArray.GetSize() - 1,  nsameBinsArray.GetArray());
    fhNCellsPerClusterSameW ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterSameW ->SetYTitle("#it{n}_{cells, same T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterSameW );
    
    fhNCellsPerClusterDiffW   = new TH2F 
    ("hNCellsPerClusterDiffW ","# cells per cluster with #it{w} in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
     //nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
     eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
     nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
    fhNCellsPerClusterDiffW ->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterDiffW ->SetYTitle("#it{n}_{cells, diff T-Card}^{#it{w}}");
    outputContainer->Add(fhNCellsPerClusterDiffW );
    
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
    
    // Cluster Exoticity other definitions
    //
    
    fhExoSame = new TH2F 
    ("hExoSame","cell #it{F}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoSame->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoSame->SetYTitle("#it{F}_{same}");
    outputContainer->Add(fhExoSame);    
    
    fhExoDiff = new TH2F 
    ("hExoDiff","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoDiff->SetYTitle("#it{F}_{diff}");
    outputContainer->Add(fhExoDiff);   
    
    fhExoSame5 = new TH2F 
    ("hExoSame5","cell #it{F}_{same-5} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoSame5->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoSame5->SetYTitle("#it{F}_{same-5}");
    outputContainer->Add(fhExoSame5);    
    
    fhExoDiff5 = new TH2F 
    ("hExoDiff5","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoDiff5->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoDiff5->SetYTitle("#it{F}_{diff-5}");
    outputContainer->Add(fhExoDiff5);   
    
    //
    fhFracEnDiffSame    = new TH2F 
    ("hFracEnDiffSame","cell #Sigma #it{E}_{diff}/#Sigma #it{E}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSame->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSame->SetYTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
    outputContainer->Add(fhFracEnDiffSame);  
    
    fhFracNCellDiffSame = new TH2F 
    ("hFracNCellDiffSame","cell #it{n}_{diff}/#it{n}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSame->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSame->SetYTitle("#it{n}_{diff}/#it{n}_{same}");
    outputContainer->Add(fhFracNCellDiffSame); 
    
    fhFracEnNCellDiffSame = new TH2F 
    ("hFracEnNCellDiffSame","cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSame->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnNCellDiffSame->SetYTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
    outputContainer->Add(fhFracEnNCellDiffSame); 
    
    fhFracEnDiffSameW    = new TH2F 
    ("hFracEnDiffSameW","cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSameW->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSameW->SetYTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
    outputContainer->Add(fhFracEnDiffSameW);  
    
    fhFracNCellDiffSameW = new TH2F 
    ("hFracNCellDiffSameW","cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSameW->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSameW->SetYTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
    outputContainer->Add(fhFracNCellDiffSameW); 
    
    fhFracEnNCellDiffSameW = new TH2F 
    ("hFracEnNCellDiffSameW","cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSameW->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnNCellDiffSameW->SetYTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
    outputContainer->Add(fhFracEnNCellDiffSameW); 
    
    fhFracEnDiffSame5    = new TH2F 
    ("hFracEnDiffSame5","cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnDiffSame5->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSame5->SetYTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
    outputContainer->Add(fhFracEnDiffSame5);  
    
    fhFracNCellDiffSame5 = new TH2F 
    ("hFracNCellDiffSame5","cell #it{n}_{diff}^{next}/#it{n}_{same}^{next} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracNCellDiffSame5->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSame5->SetYTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
    outputContainer->Add(fhFracNCellDiffSame5); 
    
    fhFracEnNCellDiffSame5 = new TH2F 
    ("hFracEnNCellDiffSame5","cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next}) vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
     nptbins,ptmin,ptmax, 200,0,2); 
    fhFracEnNCellDiffSame5->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnNCellDiffSame5->SetYTitle("(#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next})");
    outputContainer->Add(fhFracEnNCellDiffSame5); 
    
    //
    fhFracEnDiffSameExo    = new TH3F 
    ("hFracEnDiffSameExo","cell #Sigma #it{E}_{diff}/#Sigma #it{E}_{same} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSameExo->SetYTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
    fhFracEnDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSameExo);  
    
    fhFracNCellDiffSameExo = new TH3F 
    ("hFracNCellDiffSameExo","cell #it{n}_{diff}/#it{n}_{same} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSameExo->SetYTitle("#it{n}_{diff}/#it{n}_{same}");
    fhFracNCellDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSameExo); 
    
    fhFracEnNCellDiffSameExo = new TH3F 
    ("hFracEnNCellDiffSameExo","cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    
    fhFracEnNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnNCellDiffSameExo->SetYTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
    fhFracEnNCellDiffSameExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnNCellDiffSameExo); 
    
    fhFracEnDiffSameWExo    = new TH3F 
    ("hFracEnDiffSameWExo","cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSameWExo->SetYTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
    fhFracEnDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSameWExo);  
    
    fhFracNCellDiffSameWExo = new TH3F 
    ("hFracNCellDiffSameWExo","cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSameWExo->SetYTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
    fhFracNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSameWExo); 
    
    fhFracEnNCellDiffSameWExo = new TH3F 
    ("hFracEnNCellDiffSameWExo","cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());
    fhFracEnNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnNCellDiffSameWExo->SetYTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
    fhFracEnNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnNCellDiffSameWExo); 
    
    fhFracEnDiffSame5Exo    = new TH3F 
    ("hFracEnDiffSame5Exo","cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray()); 
    fhFracEnDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracEnDiffSame5Exo->SetYTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
    fhFracEnDiffSame5Exo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracEnDiffSame5Exo);  
    
    fhFracNCellDiffSame5Exo = new TH3F 
    ("hFracNCellDiffSame5Exo","cell #it{n}_{diff}^{next}/#it{n}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray());  
    fhFracNCellDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhFracNCellDiffSame5Exo->SetYTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
    fhFracNCellDiffSame5Exo->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhFracNCellDiffSame5Exo); 
    
    fhFracEnNCellDiffSame5Exo = new TH3F 
    ("hFracEnNCellDiffSame5Exo","cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
     //nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
         eBinsArray.GetSize() - 1,     eBinsArray.GetArray(), 
     frac2BinsArray.GetSize() - 1, frac2BinsArray.GetArray(), 
         fBinsArray.GetSize() - 1,     fBinsArray.GetArray()); 
    fhFracEnNCellDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV) ");
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
    
    fhExoticityEClusAllSameTCard = new TH2F 
    ("hExoticityEClusAllSameTCard","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}",
     //nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhExoticityEClusAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoticityEClusAllSameTCard->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticityEClusAllSameTCard);    
    
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
  }

  fhCellEnSameExo = new TH3F 
  ("hCellEnSameExo","#it{E}_{cluster} vs #it{E}_{cell}^{same} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhCellEnSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnSameExo->SetYTitle("#it{E}_{cell}^{same} (GeV) ");
  fhCellEnSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnSameExo);    
  
  fhCellEnDiffExo = new TH3F 
  ("hCellEnDiffExo","#it{E}_{cluster} vs #it{E}_{cell}^{diff} vs #it{F}_{+}",
   //nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    fBinsArray.GetSize() - 1,  fBinsArray.GetArray());
  fhCellEnDiffExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnDiffExo->SetYTitle("#it{E}_{cell}^{diff} (GeV) ");
  fhCellEnDiffExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnDiffExo);     

  fhCellEnNCellWOpenTime = new TH3F 
  ("hCellEnNCellWOpenTime","#it{E}_{cluster} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}, no time cut",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnNCellWOpenTime->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellWOpenTime);  

  fhCellEnNCellWEMaxOpenTime = new TH3F 
  ("hCellEnNCellWEMaxOpenTime","#it{E}_{cell}^{max} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}, no time cut",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellWEMaxOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhCellEnNCellWEMaxOpenTime->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellWEMaxOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellWEMaxOpenTime);  
  
  fhCellEnNCellW = new TH3F 
  ("hCellEnNCellW","#it{E}_{cluster} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellW->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnNCellW->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellW);  
 
  fhCellEnNCellWEMax = new TH3F 
  ("hCellEnNCellWEMax","#it{E}_{cell}^{max} vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   e2BinsArray.GetSize() - 1, e2BinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellEnNCellWEMax->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhCellEnNCellWEMax->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellWEMax->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellWEMax);  

  fhCellTimeDiffNCellWOpenTime = new TH3F 
  ("hCellTimeDiffNCellWOpenTime","#it{E}_{cluster} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}, no time cut",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellTimeDiffNCellWOpenTime->SetYTitle("#Delta #it{t}^{max-sec} (ns) ");
  fhCellTimeDiffNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellWOpenTime);  
  
  fhCellTimeDiffNCellWEMaxOpenTime = new TH3F 
  ("hCellTimeDiffNCellWEMaxOpenTime","#it{E}_{cell}^{max} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}, no time cut",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellWEMaxOpenTime->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhCellTimeDiffNCellWEMaxOpenTime->SetYTitle("#Delta #it{t}^{max-sec} (ns) ");
  fhCellTimeDiffNCellWEMaxOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellWEMaxOpenTime);  
  
  fhCellTimeDiffNCellW = new TH3F 
  ("hCellTimeDiffNCellW","#it{E}_{cluster} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellW->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellTimeDiffNCellW->SetYTitle("#Delta #it{t}^{max-sec} (ns) ");
  fhCellTimeDiffNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellW);  
  
  fhCellTimeDiffNCellWEMax = new TH3F 
  ("hCellTimeDiffNCellWEMax","#it{E}_{cell}^{max} vs #Delta #it{t}^{max-sec} vs #it{n}_{cell}^{#it{w}}",
   //nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
    eBinsArray.GetSize() - 1,  eBinsArray.GetArray(), 
   tdBinsArray.GetSize() - 1, tdBinsArray.GetArray(), 
    nBinsArray.GetSize() - 1,  nBinsArray.GetArray());
  fhCellTimeDiffNCellWEMax->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhCellTimeDiffNCellWEMax->SetYTitle("#Delta #it{t}^{max-sec} (ns) ");
  fhCellTimeDiffNCellWEMax->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellTimeDiffNCellWEMax);  
  
  fhCellMaxClusterEn = new TH2F 
  ("hCellMaxClusterEn","#it{E}_{cluster}^{org} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray());
  fhCellMaxClusterEn->SetXTitle("#it{E}_{cluster}^{org} (GeV) ");
  fhCellMaxClusterEn->SetYTitle("#it{E}_{cell}^{max} (GeV) ");
  outputContainer->Add(fhCellMaxClusterEn);  

  fhCellMaxClusterEnOpenTime = new TH2F 
  ("hCellMaxClusterEnOpenTime","#it{E}_{cluster}^{org} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1, no time cut",
   //nptbins,ptmin,ptmax, nptbins,ptmin,ptmax); 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
   eBinsArray.GetSize() - 1, eBinsArray.GetArray());
  fhCellMaxClusterEnOpenTime->SetXTitle("#it{E}_{cluster}^{org} (GeV) ");
  fhCellMaxClusterEnOpenTime->SetYTitle("#it{E}_{cell}^{max} (GeV) ");
  outputContainer->Add(fhCellMaxClusterEnOpenTime);    
  
  fhCellMaxClusterEnRatio = new TH2F 
  ("hCellMaxClusterEnRatio","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 202,0,1.01); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray());
  fhCellMaxClusterEnRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatio->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  outputContainer->Add(fhCellMaxClusterEnRatio);  
 
  fhCellMaxClusterEnRatioOpenTime = new TH2F 
  ("hCellMaxClusterEnRatioOpenTime","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org}, #it{n}_{cluster}^{cell} > 1, no time cut",
   //nptbins,ptmin,ptmax, 202,0,1.01); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray());
  fhCellMaxClusterEnRatioOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioOpenTime->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  outputContainer->Add(fhCellMaxClusterEnRatioOpenTime);  
  
  fhCellMaxClusterEnRatioNCellWOpenTime = new TH3F 
  ("hCellMaxClusterEnRatioNCellWOpenTime","#it{E}_{cluster}vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{n}_{cell}^{#it{w}}, #it{n}_{cluster}^{cell} > 1, no time cut",
   //nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,    nBinsArray.GetArray());
  fhCellMaxClusterEnRatioNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellWOpenTime->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  fhCellMaxClusterEnRatioNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellMaxClusterEnRatioNCellWOpenTime);  
  
  fhCellMaxClusterEnRatioNCellW = new TH3F 
  ("hCellMaxClusterEnRatioNCellW","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{n}_{cell}^{#it{w}}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      nBinsArray.GetSize() - 1,    nBinsArray.GetArray());
  fhCellMaxClusterEnRatioNCellW->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellW->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}^{org}");
  fhCellMaxClusterEnRatioNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellMaxClusterEnRatioNCellW);  
 
  fhCellMaxClusterEnRatioExo = new TH3F 
  ("hCellMaxClusterEnRatioExo","#it{E}_{cluster} vs #it{E}_{cell}^{max}/#it{E}_{cluster}^{org} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   //nptbins,ptmin,ptmax, 101,0,1.01, nexobins,exomin,exomax); 
      eBinsArray.GetSize() - 1,    eBinsArray.GetArray(), 
   fracBinsArray.GetSize() - 1, fracBinsArray.GetArray(), 
      fBinsArray.GetSize() - 1,    fBinsArray.GetArray());
  fhCellMaxClusterEnRatioExo->SetXTitle("#it{E}_{cluster} (GeV) ");
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
  fhEtaPhiGridExoEnCut->SetXTitle("column-#eta ");
  fhEtaPhiGridExoEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridExoEnCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhEtaPhiGridExoEnCut);
  
  if ( fFillSameDiffFracHisto )
  {
    fhEtaPhiGridExoEnCutSameFracCut = new TH3F 
    ("hEtaPhiGridExoEnCutSameFracCut",
     Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1, #it{n}_{cells}^{same}/(#it{n}_{cells}-1)=1",fEMinForExo),
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
       fBinsArray.GetSize() - 1,   fBinsArray.GetArray());
    fhEtaPhiGridExoEnCutSameFracCut->SetXTitle("column-#eta ");
    fhEtaPhiGridExoEnCutSameFracCut->SetYTitle("row-#varphi (rad)");
    fhEtaPhiGridExoEnCutSameFracCut->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEtaPhiGridExoEnCutSameFracCut);
  }
  
  fhEtaPhiGridEnExoCut  = new TH3F 
  ("hEtaPhiGridEnExoCut", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{F}_{+} > %2.2f",fExoCut),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
  fhEtaPhiGridEnExoCut->SetXTitle("column-#eta ");
  fhEtaPhiGridEnExoCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridEnExoCut->SetZTitle("#it{E} (GeV)");
  outputContainer->Add(fhEtaPhiGridEnExoCut);
  
  fhEtaPhiGridEnHighNCells  = new TH3F 
  ("hEtaPhiGridEnHighNCells", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}^{#it{w}} > %d",fNCellHighCut),
   //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
   colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
   rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
     eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
  fhEtaPhiGridEnHighNCells->SetXTitle("column-#eta ");
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
  fhEtaPhiGridNCellEnCut->SetXTitle("column-#eta ");
  fhEtaPhiGridNCellEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridNCellEnCut->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhEtaPhiGridNCellEnCut);
  
  if ( fFill1CellHisto )
  {
    fhEtaPhiGridEn1Cell  = new TH3F 
    ("hEtaPhiGridEn1Cell", "colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}=1",
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins,ptmin,ptmax); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
       eBinsArray.GetSize() - 1,   eBinsArray.GetArray());
    fhEtaPhiGridEn1Cell->SetXTitle("column-#eta ");
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
  
  fhTimeDiffAmpClusCellExo  = new TH3F 
  ("hTimeDiffAmpClusCellExo",
   Form("#it{E}_{cell i} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1"),
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
  }
  
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
  
  // Track matching
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
    fhCellExoAmp->SetXTitle("#it{E}_{cell} (GeV) ");
    fhCellExoAmp->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoAmp);    
    
    fhCellExoAmpTime = new TH3F 
    ("hCellExoAmpTime","Cell #it{F}_{+} vs #it{E}_{cell} vs time",
     //nptbins,ptmin,ptmax/2, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
     eBinsArray.GetSize() - 1, eBinsArray.GetArray(), 
     tBinsArray.GetSize() - 1, tBinsArray.GetArray(), 
     fBinsArray.GetSize() - 1, fBinsArray.GetArray());
    fhCellExoAmpTime->SetXTitle("#it{E}_{cell} (GeV) ");
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
     "Cell hits row-column vs cluster time for #it{n}_{cell}^{w} > 20", 
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
     "Cell hits row-column vs cluster time for #it{n}_{cell}^{w} > 12", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTimeHighNCell12->SetYTitle("row (phi direction)");
    fhCellGridTimeHighNCell12->SetXTitle("column (eta direction)");
    fhCellGridTimeHighNCell12->SetZTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhCellGridTimeHighNCell12);
    
    fhCellGridTime05    = new TH3F 
    ("hCellGridTime05",
     "Cell hits row-column vs cell time for #it{E}_{cell} > 0.5 GeV", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTime05->SetYTitle("row (phi direction)");
    fhCellGridTime05->SetXTitle("column (eta direction)");
    fhCellGridTime05->SetZTitle("#it{t}_{cell} (ns)");
    outputContainer->Add(fhCellGridTime05);

    fhCellGridTime1    = new TH3F 
    ("hCellGridTime1",
     "Cell hits row-column vs cell time for #it{E}_{cell} > 1 GeV", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTime1->SetYTitle("row (phi direction)");
    fhCellGridTime1->SetXTitle("column (eta direction)");
    fhCellGridTime1->SetZTitle("#it{t}_{cell} (ns)");
    outputContainer->Add(fhCellGridTime1);

    fhCellGridTime2    = new TH3F 
    ("hCellGridTime2",
     "Cell hits row-column vs cell time for #it{E}_{cell} > 2 GeV", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTime2->SetYTitle("row (phi direction)");
    fhCellGridTime2->SetXTitle("column (eta direction)");
    fhCellGridTime2->SetZTitle("#it{t}_{cell} (ns)");
    outputContainer->Add(fhCellGridTime2);    
 
    fhCellGridTime5    = new TH3F 
    ("hCellGridTime5",
     "Cell hits row-column vs cell time for #it{E}_{cell} > 5 GeV", 
     //ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, 61,-610,610); 
     colBinsArray.GetSize() - 1, colBinsArray.GetArray(), 
     rowBinsArray.GetSize() - 1, rowBinsArray.GetArray(), 
      t2BinsArray.GetSize() - 1,  t2BinsArray.GetArray());
    fhCellGridTime5->SetYTitle("row (phi direction)");
    fhCellGridTime5->SetXTitle("column (eta direction)");
    fhCellGridTime5->SetZTitle("#it{t}_{cell} (ns)");
    outputContainer->Add(fhCellGridTime5);    
    
    fhCellTimeNCell20 = new TH2F 
    ("hCellTimeNCell20",
     "All cell time vs cluster time for cluster #it{n}_{cell}^{w} > 20", 
     tBinsArray.GetSize() - 1,   tBinsArray.GetArray(),
     tBinsArray.GetSize() - 1,   tBinsArray.GetArray());
    fhCellTimeNCell20->SetYTitle("#it{t}_{cell} (ns)");
    fhCellTimeNCell20->SetXTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhCellTimeNCell20);
    
    fhCellTimeNCell12 = new TH2F 
    ("hCellTimeNCell12",
     "All cell time vs cluster time for cluster #it{n}_{cell}^{w} > 12", 
     tBinsArray.GetSize() - 1,   tBinsArray.GetArray(),
     tBinsArray.GetSize() - 1,   tBinsArray.GetArray());
    fhCellTimeNCell12->SetYTitle("#it{t}_{cell} (ns)");
    fhCellTimeNCell12->SetXTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhCellTimeNCell12);
    
    fhSumEnCells05 = new TH1F 
    ("hSumEnCells05","#Sigma #it{E} for #it{E}_{i} > 0.5 GeV", 1000,0,1000);
    fhSumEnCells05->SetYTitle("Counts per event");
    fhSumEnCells05->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells05);                     

    fhSumEnCells1 = new TH1F 
    ("hSumEnCells1","#Sigma #it{E} for #it{E}_{i} > 1 GeV", 1000,0,1000);
    fhSumEnCells1->SetYTitle("Counts per event");
    fhSumEnCells1->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells1);         
  
    fhSumEnCells2 = new TH1F 
    ("hSumEnCells2","#Sigma #it{E} for #it{E}_{i} > 2 GeV", 1000,0,1000);
    fhSumEnCells2->SetYTitle("Counts per event");
    fhSumEnCells2->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells2);    
    
    fhSumEnCells5 = new TH1F 
    ("hSumEnCells5","#Sigma #it{E} for #it{E}_{i} > 5 GeV", 1000,0,1000);
    fhSumEnCells5->SetYTitle("Counts per event");
    fhSumEnCells5->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells5);    


    fhNCells05 = new TH1F 
    ("hNCells05","#it{n}_{cells}for #it{E}_{i} > 0.5 GeV", 17000,0,17000);
    fhNCells05->SetYTitle("Counts per event");
    fhNCells05->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells05);                     
    
    fhNCells1 = new TH1F 
    ("hNCells1","#it{n}_{cells} for #it{E}_{i} > 1 GeV", 17000,0,17000);
    fhNCells1->SetYTitle("Counts per event");
    fhNCells1->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells1);         
    
    fhNCells2 = new TH1F 
    ("hNCells2","#it{n}_{cells} for #it{E}_{i} > 2 GeV", 17000,0,17000);
    fhNCells2->SetYTitle("Counts per event");
    fhNCells2->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells2);    
    
    fhNCells5 = new TH1F 
    ("hNCells5","#it{n}_{cells} for #it{E}_{i} > 5 GeV", 17000,0,17000);
    fhNCells5->SetYTitle("Counts per event");
    fhNCells5->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells5);        
    
    fhAverSumEnCells05 = new TH1F 
    ("hAverSumEnCells05","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 0.5 GeV", 1000,0,1000);
    fhAverSumEnCells05->SetYTitle("Counts per event");
    fhAverSumEnCells05->SetXTitle("#Sigma #it{E}/#it{n}_{cells} (GeV)");
    outputContainer->Add(fhAverSumEnCells05);                     
    
    fhAverSumEnCells1 = new TH1F 
    ("hAverSumEnCells1","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 1 GeV", 1000,0,1000);
    fhAverSumEnCells1->SetYTitle("Counts per event");
    fhAverSumEnCells1->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhAverSumEnCells1);         
    
    fhAverSumEnCells2 = new TH1F 
    ("hAverSumEnCells2","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 2 GeV", 1000,0,1000);
    fhAverSumEnCells2->SetYTitle("Counts per event");
    fhAverSumEnCells2->SetXTitle("#Sigma #it{E}/#it{n}_{cells} (GeV)");
    outputContainer->Add(fhAverSumEnCells2);    
    
    fhAverSumEnCells5 = new TH1F 
    ("hAverSumEnCells5","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 5 GeV", 1000,0,1000);
    fhAverSumEnCells5->SetYTitle("Counts per event");
    fhAverSumEnCells5->SetXTitle("#Sigma #it{E}/#it{n}_{cells}");
    outputContainer->Add(fhAverSumEnCells5);    
    
    fhFracSumEnCells1 = new TH1F 
    ("hFracSumEnCells1","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 1 GeV", 201,0,1.005);
    fhFracSumEnCells1->SetYTitle("Counts per event");
    fhFracSumEnCells1->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells1);         
    
    fhFracSumEnCells2 = new TH1F 
    ("hFracSumEnCells2","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 2 GeV", 201,0,1.005);
    fhFracSumEnCells2->SetYTitle("Counts per event");
    fhFracSumEnCells2->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells2);    
    
    fhFracSumEnCells5 = new TH1F 
    ("hFracSumEnCells5","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 5 GeV", 201,0,1.005);
    fhFracSumEnCells5->SetYTitle("Counts per event");
    fhFracSumEnCells5->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells5);    

    fhFracNCells1 = new TH1F 
    ("hFracNCells1","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 1 GeV", 201,0,1.005);
    fhFracNCells1->SetYTitle("Counts per event");
    fhFracNCells1->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells1);         
    
    fhFracNCells2 = new TH1F 
    ("hFracNCells2","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 2 GeV", 201,0,1.005);
    fhFracNCells2->SetYTitle("Counts per event");
    fhFracNCells2->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells2);    
    
    fhFracNCells5 = new TH1F 
    ("hFracNCells5","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 5 GeV", 201,0,1.005);
    fhFracNCells5->SetYTitle("Counts per event");
    fhFracNCells5->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells5);       
    
    ///////////////////
    
    fhSumEnCells05NHigh20 = new TH1F 
    ("hSumEnCells05NHigh20","#Sigma #it{E} for #it{E}_{i} > 0.5 GeV", 1000,0,1000);
    fhSumEnCells05NHigh20->SetYTitle("Counts per event");
    fhSumEnCells05NHigh20->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells05NHigh20);                     
    
    fhSumEnCells1NHigh20 = new TH1F 
    ("hSumEnCells1NHigh20","#Sigma #it{E} for #it{E}_{i} > 1 GeV", 1000,0,1000);
    fhSumEnCells1NHigh20->SetYTitle("Counts per event");
    fhSumEnCells1NHigh20->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells1NHigh20);         
    
    fhSumEnCells2NHigh20 = new TH1F 
    ("hSumEnCells2NHigh20","#Sigma #it{E} for #it{E}_{i} > 2 GeV", 1000,0,1000);
    fhSumEnCells2NHigh20->SetYTitle("Counts per event");
    fhSumEnCells2NHigh20->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells2NHigh20);    
    
    fhSumEnCells5NHigh20 = new TH1F 
    ("hSumEnCells5NHigh20","#Sigma #it{E} for #it{E}_{i} > 5 GeV", 1000,0,1000);
    fhSumEnCells5NHigh20->SetYTitle("Counts per event");
    fhSumEnCells5NHigh20->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhSumEnCells5NHigh20);    
    
    
    fhNCells05NHigh20 = new TH1F 
    ("hNCells05NHigh20","#it{n}_{cells}for #it{E}_{i} > 0.5 GeV", 17000,0,17000);
    fhNCells05NHigh20->SetYTitle("Counts per event");
    fhNCells05NHigh20->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells05NHigh20);                     
    
    fhNCells1NHigh20 = new TH1F 
    ("hNCells1NHigh20","#it{n}_{cells} for #it{E}_{i} > 1 GeV", 17000,0,17000);
    fhNCells1NHigh20->SetYTitle("Counts per event");
    fhNCells1NHigh20->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells1NHigh20);         
    
    fhNCells2NHigh20 = new TH1F 
    ("hNCells2NHigh20","#it{n}_{cells} for #it{E}_{i} > 2 GeV", 17000,0,17000);
    fhNCells2NHigh20->SetYTitle("Counts per event");
    fhNCells2NHigh20->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells2NHigh20);    
    
    fhNCells5NHigh20 = new TH1F 
    ("hNCells5NHigh20","#it{n}_{cells} for #it{E}_{i} > 5 GeV", 17000,0,17000);
    fhNCells5NHigh20->SetYTitle("Counts per event");
    fhNCells5NHigh20->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells5NHigh20);        
    
    fhAverSumEnCells05NHigh20 = new TH1F 
    ("hAverSumEnCells05NHigh20","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 0.5 GeV", 1000,0,1000);
    fhAverSumEnCells05NHigh20->SetYTitle("Counts per event");
    fhAverSumEnCells05NHigh20->SetXTitle("#Sigma #it{E}/#it{n}_{cells} (GeV)");
    outputContainer->Add(fhAverSumEnCells05NHigh20);                     
    
    fhAverSumEnCells1NHigh20 = new TH1F 
    ("hAverSumEnCells1NHigh20","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 1 GeV", 1000,0,1000);
    fhAverSumEnCells1NHigh20->SetYTitle("Counts per event");
    fhAverSumEnCells1NHigh20->SetXTitle("#Sigma #it{E} (GeV)");
    outputContainer->Add(fhAverSumEnCells1NHigh20);         
    
    fhAverSumEnCells2NHigh20 = new TH1F 
    ("hAverSumEnCells2NHigh20","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 2 GeV", 1000,0,1000);
    fhAverSumEnCells2NHigh20->SetYTitle("Counts per event");
    fhAverSumEnCells2NHigh20->SetXTitle("#Sigma #it{E}/#it{n}_{cells} (GeV)");
    outputContainer->Add(fhAverSumEnCells2NHigh20);    
    
    fhAverSumEnCells5NHigh20 = new TH1F 
    ("hAverSumEnCells5NHigh20","#Sigma #it{E}/#it{n}_{cells} for #it{E}_{i} > 5 GeV", 1000,0,1000);
    fhAverSumEnCells5NHigh20->SetYTitle("Counts per event");
    fhAverSumEnCells5NHigh20->SetXTitle("#Sigma #it{E}/#it{n}_{cells}");
    outputContainer->Add(fhAverSumEnCells5NHigh20);    
    
    fhFracSumEnCells1NHigh20 = new TH1F 
    ("hFracSumEnCells1NHigh20","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 1 GeV", 201,0,1.005);
    fhFracSumEnCells1NHigh20->SetYTitle("Counts per event");
    fhFracSumEnCells1NHigh20->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells1NHigh20);         
    
    fhFracSumEnCells2NHigh20 = new TH1F 
    ("hFracSumEnCells2NHigh20","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 2 GeV", 201,0,1.005);
    fhFracSumEnCells2NHigh20->SetYTitle("Counts per event");
    fhFracSumEnCells2NHigh20->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells2NHigh20);    
    
    fhFracSumEnCells5NHigh20 = new TH1F 
    ("hFracSumEnCells5NHigh20","#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5} for #it{E}_{i} > 5 GeV", 201,0,1.005);
    fhFracSumEnCells5NHigh20->SetYTitle("Counts per event");
    fhFracSumEnCells5NHigh20->SetXTitle("#Sigma #it{E}/#Sigma #it{E}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracSumEnCells5NHigh20);    
    
    fhFracNCells1NHigh20 = new TH1F 
    ("hFracNCells1NHigh20","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 1 GeV", 201,0,1.005);
    fhFracNCells1NHigh20->SetYTitle("Counts per event");
    fhFracNCells1NHigh20->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells1NHigh20);         
    
    fhFracNCells2NHigh20 = new TH1F 
    ("hFracNCells2NHigh20","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 2 GeV", 201,0,1.005);
    fhFracNCells2NHigh20->SetYTitle("Counts per event");
    fhFracNCells2NHigh20->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells2NHigh20);    
    
    fhFracNCells5NHigh20 = new TH1F 
    ("hFracNCells5NHigh20","#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5} for #it{E}_{i} > 5 GeV", 201,0,1.005);
    fhFracNCells5NHigh20->SetYTitle("Counts per event");
    fhFracNCells5NHigh20->SetXTitle("#it{n}_{cells}/#it{n}_{cells}_{#it{E}_{i}>0.5}");
    outputContainer->Add(fhFracNCells5NHigh20);       
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
  printf("Exoticity cut: %2.1f \n", fExoCut) ;
  printf("NCell cut: %d \n", fNCellHighCut) ;
  printf("Time range: [%2.2f,%2.2f] ns\n",fTimeCutMin,fTimeCutMax);
  printf("Fill cell histo : %d GeV/c\n", fFillCellHisto) ;
  printf("Fill 1 cell cluster histo : %d GeV/c\n", fFill1CellHisto) ;
  printf("Fill Matching histo : %d GeV/c\n", fFillMatchingHisto) ;
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
  
  // Cells  
  if ( fFillCellHisto )  CellHistograms(cells);
  
  AliDebug(1,"End");
}

