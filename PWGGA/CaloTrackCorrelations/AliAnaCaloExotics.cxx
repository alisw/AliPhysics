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
fFillCellHisto(0),                     fFill1CellHisto(0),                    
fFillMatchingHisto(0),                 fConstantTimeShift(0),                 
fClusterMomentum(),         

// Histograms
fhExoticityEClus(0),                   fhExoticityEClusAllSameTCard(0),
fhExoticityEClusPerSM(0),              fhExoticityEMaxCell(0),
fhExoticityEClusTrackMatch(0),         fhExoticity1Cell(0),

fhNCellsPerCluster(0),                 fhNCellsPerClusterAllSameTCard(0),      
fhNCellsPerClusterPerSM(0),            fhNCellsPerClusterWPerSM(0),              
fhNCellsPerClusterExo(0),              
fhNCellsPerClusterTrackMatch(0),       fhNCellsPerClusterExoTrackMatch(0),    
fhNCellsPerClusterM02(0),

fhEtaPhiGridExoEnCut(0),               fhEtaPhiGridExoEnCutSameFracCut(0),               
fhEtaPhiGridEnExoCut(0),               fhEtaPhiGridEn1Cell(0),
fhEtaPhiGridEnHighNCells(0),           fhEtaPhiGridNCellEnCut(0),

fhTimeEnergyExo(0),                    fhTimeEnergy1Cell(0),                
fhTimeDiffClusCellExo(0),              fhTimeDiffAmpClusCellExo(0),   
fhTimeEnergyM02(0),                    fhTimeDiffClusCellM02(0),  
fhTimeEnergyNCells(0),                 fhTimeEnergyNCellsW(0),                 fhTimeNCellCut(0),                

fhM02EnergyNCell(0),                   fhM02EnergyAllSameTCard(0),      
fhM02EnergyExo(0),                     fhM02EnergyExoZoomIn(0),               
fhM20EnergyExoM02MinCut(0),                                                 

// Other ncell
fhNCellsPerClusterW (0),               
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
fhCellExoAmpTime(0),                    fhCellExoGrid(0)                     
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
  }
  
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

  } // Cell loop
  
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
  Int_t  nCaloClusters = caloClusters->GetEntriesFast() ;
  Int_t  bc            = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  // Get vertex, not needed.
  Double_t v[3] = {0,0,0};
  //GetReader()->GetVertex(v);
  
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
     
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
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
    
    Float_t exoticity = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc)/ampMax;
    
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

    AliDebug(1,Form("cluster: E %2.3f, F+ %2.3f, eta %2.3f, phi %2.3f, col %d, row %d, ncells %d,"
                    "match %d; cell max: id %d, en %2.3f, time %2.3f, m02 %2.2f",
                    en,exoticity,eta,phi*TMath::RadToDeg(), icolMaxAbs, irowMaxAbs, nCaloCellsPerCluster,
                    matched, absIdMax,ampMax,tmax,m02));  

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
      fhTimeEnergyExo  ->Fill(en, tmax, exoticity, GetEventWeight());
      fhTimeEnergyM02  ->Fill(en, tmax, m02      , GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhTimeEnergy1Cell->Fill(en, tmax,            GetEventWeight());
    
//    if(maxCellFraction < 0.2) printf("frac %2.2f, en cell %2.2f, en cluster %2.2f, refrac %2.2f\n",
//                                     maxCellFraction,ampMax,en,ampMax/en);
    fhCellMaxClusterEnRatioNCellWOpenTime->Fill(en, 1-maxCellFraction, nCellW, GetEventWeight());
    
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId     = clus->GetCellsAbsId()[ipos];  
      
      Float_t amp       = cells->GetCellAmplitude(absId);
      
      if ( absId == absIdMax || amp < 0.1 ) continue;
      
      fhCellEnNCellWOpenTime->Fill(en, amp, nCellW, GetEventWeight());
    }
    
    // Apply time cut for other parameters
    //
    if ( tmax < fTimeCutMin || tmax> fTimeCutMax ) continue;

    // Exoticity
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      fhExoticityEClus   ->Fill(en    , exoticity, GetEventWeight());
      fhExoticityEMaxCell->Fill(ampMax, exoticity, GetEventWeight());
      
      fhExoticityEClusPerSM ->Fill(en, exoticity, nSM, GetEventWeight());
      
      if ( matched && fFillMatchingHisto )
        fhExoticityEClusTrackMatch->Fill(en, exoticity, GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhExoticity1Cell->Fill(en, exoticity, GetEventWeight());
    
    // N cells per cluster
    //
    fhNCellsPerCluster               ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
    fhNCellsPerClusterExo            ->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    fhNCellsPerClusterM02            ->Fill(en, nCaloCellsPerCluster, m02      , GetEventWeight());
    
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

    fhCellMaxClusterEnRatioNCellW->Fill(en, 1-maxCellFraction, nCellW   , GetEventWeight());
    fhCellMaxClusterEnRatioExo   ->Fill(en, 1-maxCellFraction, exoticity, GetEventWeight());
    
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

      fhCellEnNCellW->Fill(en, amp, nCellW, GetEventWeight());

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
    
    fhNCellsPerClusterW       ->Fill(en, nCellW     , GetEventWeight());
    fhNCellsPerClusterSame    ->Fill(en, nCellSame  , GetEventWeight());
    fhNCellsPerClusterDiff    ->Fill(en, nCellDiff  , GetEventWeight());
    fhNCellsPerClusterSame5   ->Fill(en, nCellSame5 , GetEventWeight());
    fhNCellsPerClusterDiff5   ->Fill(en, nCellDiff5 , GetEventWeight());
    fhNCellsPerClusterSameW   ->Fill(en, nCellSameW , GetEventWeight());
    fhNCellsPerClusterDiffW   ->Fill(en, nCellDiffW , GetEventWeight());    
    fhNCellsPerClusterSameDiff->Fill(en, nCellSame, nCellDiff, GetEventWeight());
    
    if ( nCaloCellsPerCluster > 1 )
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
 
      
      fhFracEnDiffSameExo      ->Fill(en, fracEnDiffSame     , exoticity, GetEventWeight()); 
      fhFracNCellDiffSameExo   ->Fill(en, fracNCellDiffSame  , exoticity, GetEventWeight());
      fhFracEnNCellDiffSameExo ->Fill(en, fracEnNCellDiffSame, exoticity, GetEventWeight());
      
      fhFracEnDiffSameWExo     ->Fill(en, fracEnDiffSameW     , GetEventWeight()); 
      fhFracNCellDiffSameWExo  ->Fill(en, fracNCellDiffSameW  , GetEventWeight());
      fhFracEnNCellDiffSameWExo->Fill(en, fracEnNCellDiffSameW, GetEventWeight());
      
      fhFracEnDiffSame5Exo     ->Fill(en, fracEnDiffSame5     , GetEventWeight()); 
      fhFracNCellDiffSame5Exo  ->Fill(en, fracNCellDiffSame5  , GetEventWeight());
      fhFracEnNCellDiffSame5Exo->Fill(en, fracEnNCellDiffSame5, GetEventWeight());
      
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
        if ( en > fEMinForExo )
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
      }
    }
    
    // Shower shape
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      fhM02EnergyNCell->Fill(en, m02, nCaloCellsPerCluster, GetEventWeight());
      
      fhM02EnergyExo->Fill(en, m02, exoticity, GetEventWeight());
      
      fhM02EnergyExoZoomIn->Fill(en, m02, exoticity, GetEventWeight());
      
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
 
  snprintf(onePar,buffersize,"fill 1 cell histo: %d;",fFill1CellHisto) ;
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
    
  Int_t totalSM = fLastModule-fFirstModule+1;

  // Histogram binning and ranges
  // 
  Int_t nptbins       = GetHistogramRanges()->GetHistoPtBins(); 	        
  Float_t ptmax       = GetHistogramRanges()->GetHistoPtMax();           
  Float_t ptmin       = GetHistogramRanges()->GetHistoPtMin();
  
//  Int_t nphibins      = GetHistogramRanges()->GetHistoPhiBins();           
//  Float_t phimax      = GetHistogramRanges()->GetHistoPhiMax();          
//  Float_t phimin      = GetHistogramRanges()->GetHistoPhiMin();
//  
//  Int_t netabins      = GetHistogramRanges()->GetHistoEtaBins();          
//  Float_t etamax      = GetHistogramRanges()->GetHistoEtaMax();          
//  Float_t etamin      = GetHistogramRanges()->GetHistoEtaMin();  
  
  Int_t ntimebins     = GetHistogramRanges()->GetHistoTimeBins();         
  Float_t timemax     = GetHistogramRanges()->GetHistoTimeMax();         
  Float_t timemin     = GetHistogramRanges()->GetHistoTimeMin();       
 
  Int_t nceclbins     = GetHistogramRanges()->GetHistoNClusterCellBins(); 
  Int_t   nceclmax    = GetHistogramRanges()->GetHistoNClusterCellMax(); 
  Int_t   nceclmin    = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  
  Int_t ssbins        = GetHistogramRanges()->GetHistoShowerShapeBins();  
  Float_t ssmax       = GetHistogramRanges()->GetHistoShowerShapeMax();  
  Float_t ssmin       = GetHistogramRanges()->GetHistoShowerShapeMin();
  
  Int_t tdbins        = GetHistogramRanges()->GetHistoDiffTimeBins() ;    
  Float_t tdmax       = GetHistogramRanges()->GetHistoDiffTimeMax();     
  Float_t tdmin       = GetHistogramRanges()->GetHistoDiffTimeMin();
  
  // TM residuals
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  Int_t nPoverEbins   = GetHistogramRanges()->GetHistoEOverPBins();       
  Float_t eOverPmax   = GetHistogramRanges()->GetHistoEOverPMax();       
  Float_t eOverPmin   = GetHistogramRanges()->GetHistoEOverPMin();
  
  // Exoticity
  Int_t nexobins  = 201; Float_t exomin  = -1 ; Float_t exomax  = 1.01;
  // For TH3, reduced binning
  Int_t nexobinsS = 82 ; Float_t exominS = 0.7975; Float_t exomaxS = 1.0025;
    
  // Cell column-row histograms, see base class for data members setting
  //fNMaxColsFull+2,-1.5,fNMaxColsFull+0.5, fNMaxRowsFull+2,-1.5,fNMaxRowsFull+0.5
  Int_t   ncolcell    = fNMaxColsFull+2;
  Float_t colcellmin  = -1.5;
  Float_t colcellmax  = fNMaxColsFull+0.5;
  
  Int_t   nrowcell    = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  Float_t rowcellmin  = fNMaxRowsFullMin-1.5;
  Float_t rowcellmax  = fNMaxRowsFullMax+0.5;
  
  //
  // Init histograms
  //
  
  // Cluster Exoticity 2D
  //
  fhExoticityEClus = new TH2F 
  ("hExoticityEClus","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClus->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClus);    

  fhExoticityEClusAllSameTCard = new TH2F 
  ("hExoticityEClusAllSameTCard","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cell} > 1, #it{n}_{cells} = #it{n}_{cells-same}",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEClusAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClusAllSameTCard->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClusAllSameTCard);    
  
  fhExoticityEMaxCell = new TH2F 
  ("hExoticityEMaxCell","cell #it{F}_{+} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhExoticityEMaxCell->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEMaxCell);    
  
  fhExoticityEClusPerSM = new TH3F 
  ("hExoticityEClusPerSM","cell #it{F}_{+} vs #it{E}_{cluster}, vs SM, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax, totalSM,fFirstModule-0.5,fLastModule+0.5); 
  fhExoticityEClusPerSM->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClusPerSM->SetYTitle("#it{F}_{+}");
  fhExoticityEClusPerSM->SetZTitle("SM");
  outputContainer->Add(fhExoticityEClusPerSM);   
  
  if ( fFill1CellHisto )
  {
    fhExoticity1Cell = new TH2F 
    ("hExoticity1Cell","cell #it{F}_{+} vs #it{E}, #it{n}_{cluster}^{cell} = 1",
     nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
    fhExoticity1Cell->SetXTitle("#it{E} (GeV) ");
    fhExoticity1Cell->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticity1Cell);    
  }
  
  // N cells per cluster
  fhNCellsPerCluster  = new TH2F 
  ("hNCellsPerCluster","# cells per cluster vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, nceclbins*2,nceclmin,nceclmax*2); 
  fhNCellsPerCluster->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerCluster->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerCluster);

  fhNCellsPerClusterAllSameTCard  = new TH2F 
  ("hNCellsPerClusterAllSameTCard","# cells per cluster vs #it{E}_{cluster}, #it{n}_{cells}=#it{n}_{cells-same}",
   nptbins,ptmin,ptmax, 17,0,17); 
  fhNCellsPerClusterAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterAllSameTCard->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerClusterAllSameTCard);
  
  fhNCellsPerClusterExo  = new TH3F 
  ("hNCellsPerClusterExo","# cells per cluster vs #it{E}_{cluster} vs exoticity",
   nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
  fhNCellsPerClusterExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterExo->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterExo);
  
  fhNCellsPerClusterPerSM  = new TH3F 
  ("hNCellsPerClusterPerSM","# cells per cluster vs #it{E}_{cluster} vs exoticity",
   nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
  fhNCellsPerClusterPerSM->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterPerSM->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterPerSM->SetZTitle("SM");
  outputContainer->Add(fhNCellsPerClusterPerSM);

  fhNCellsPerClusterWPerSM  = new TH3F 
  ("hNCellsPerClusterWPerSM","# cells per cluster vs #it{E}_{cluster} vs exoticity",
   nptbins/2,ptmin,ptmax, nceclbins,nceclmin,nceclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
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
     nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
    fhNCellsPerClusterExoPerSM[imod]->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterExoPerSM[imod]->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterExoPerSM[imod]->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsPerClusterExoPerSM[imod]);
  }
  
  if ( fFillMatchingHisto )
  {
    fhExoticityEClusTrackMatch = new TH2F 
    ("hExoticityEClusTrackMatch","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1, track matched",
     nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
    fhExoticityEClusTrackMatch->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhExoticityEClusTrackMatch->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticityEClusTrackMatch);    
    
    fhNCellsPerClusterTrackMatch  = new TH2F 
    ("hNCellsPerClusterTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
     nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerClusterTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterTrackMatch->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterTrackMatch);
    
    fhNCellsPerClusterExoTrackMatch  = new TH3F 
    ("hNCellsPerClusterExoTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
     nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax, nexobinsS,exominS,exomaxS); 
    fhNCellsPerClusterExoTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNCellsPerClusterExoTrackMatch->SetYTitle("#it{n}_{cells}");
    fhNCellsPerClusterExoTrackMatch->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsPerClusterExoTrackMatch);
  }
  
  fhNCellsPerClusterM02  = new TH3F 
  ("hNCellsPerClusterM02","# cells per cluster vs #it{E}_{cluster} vs  #sigma^{2}_{long}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,100,0,0.5); 
  fhNCellsPerClusterM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterM02->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhNCellsPerClusterM02);
  
  
  // Different n cells definitions
  //
  fhNCellsPerClusterW   = new TH2F 
  ("hNCellsPerClusterW ","# cells per cluster with #it{w} > 0 vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerClusterW ->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterW ->SetYTitle("#it{n}_{cells}^{#it{w}}");
  outputContainer->Add(fhNCellsPerClusterW );
  
  fhNCellsPerClusterSame  = new TH2F 
  ("hNCellsPerClusterSame","# cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, 17,0,17); 
  fhNCellsPerClusterSame->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSame->SetYTitle("#it{n}_{cells, same T-Card}");
  outputContainer->Add(fhNCellsPerClusterSame);
  
  fhNCellsPerClusterDiff  = new TH2F 
  ("hNCellsPerClusterDiff","# cells per cluster in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
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
   nptbins,ptmin,ptmax, 17,0,17); 
  fhNCellsPerClusterSameW ->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameW ->SetYTitle("#it{n}_{cells, same T-Card}^{#it{w}}");
  outputContainer->Add(fhNCellsPerClusterSameW );
  
  fhNCellsPerClusterDiffW   = new TH2F 
  ("hNCellsPerClusterDiffW ","# cells per cluster with #it{w} in different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerClusterDiffW ->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterDiffW ->SetYTitle("#it{n}_{cells, diff T-Card}^{#it{w}}");
  outputContainer->Add(fhNCellsPerClusterDiffW );
  
  fhNCellsPerClusterSameDiff  = new TH3F 
  ("hNCellsPerClusterSameDiff","# cells per cluster in same vs different T-Card as max #it{E} cell vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, 17,0,17,nceclbins,nceclmin,nceclmax); 
  fhNCellsPerClusterSameDiff->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameDiff->SetYTitle("#it{n}_{cells, same T-Card}");
  fhNCellsPerClusterSameDiff->SetZTitle("#it{n}_{cells, diff T-Card}");
  outputContainer->Add(fhNCellsPerClusterSameDiff);
  
  fhNCellsPerClusterSameFrac  = new TH2F 
  ("hNCellsPerClusterSameFrac","Fraction of # cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, 101,-0.005,1.005); 
  fhNCellsPerClusterSameFrac->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameFrac->SetYTitle("#it{n}_{cells, same T-Card} / (#it{n}_{cells}-1)");
  outputContainer->Add(fhNCellsPerClusterSameFrac);
  
  fhNCellsPerClusterSameFracExo  = new TH3F 
  ("hNCellsPerClusterSameFracExo","Fraction of # cells per cluster in same T-Card as max #it{E} cell vs #it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax, 101,-0.005,1.005, nexobins,exomin,exomax); 
  fhNCellsPerClusterSameFracExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterSameFracExo->SetYTitle("#it{n}_{cells, same T-Card} / (#it{n}_{cells}-1)");
  fhNCellsPerClusterSameFracExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterSameFracExo); 
  
  // Cluster Exoticity other definitions
  //
  fhExoSame = new TH2F 
  ("hExoSame","cell #it{F}_{same} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoSame->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoSame->SetYTitle("#it{F}_{same}");
  outputContainer->Add(fhExoSame);    
  
  fhExoDiff = new TH2F 
  ("hExoDiff","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoDiff->SetYTitle("#it{F}_{diff}");
  outputContainer->Add(fhExoDiff);   
 
  fhExoSame5 = new TH2F 
  ("hExoSame5","cell #it{F}_{same-5} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoSame5->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoSame5->SetYTitle("#it{F}_{same-5}");
  outputContainer->Add(fhExoSame5);    
  
  fhExoDiff5 = new TH2F 
  ("hExoDiff5","cell #it{F}_{diff} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
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
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracEnDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracEnDiffSameExo->SetYTitle("#Sigma #it{E}_{diff}/#Sigma #it{E}_{same}");
  fhFracEnDiffSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracEnDiffSameExo);  
  
  fhFracNCellDiffSameExo = new TH3F 
  ("hFracNCellDiffSameExo","cell #it{n}_{diff}/#it{n}_{same} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracNCellDiffSameExo->SetYTitle("#it{n}_{diff}/#it{n}_{same}");
  fhFracNCellDiffSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracNCellDiffSameExo); 
  
  fhFracEnNCellDiffSameExo = new TH3F 
  ("hFracEnNCellDiffSameExo","cell (#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracEnNCellDiffSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracEnNCellDiffSameExo->SetYTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
  fhFracEnNCellDiffSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracEnNCellDiffSameExo); 
  
  fhFracEnDiffSameWExo    = new TH3F 
  ("hFracEnDiffSameWExo","cell #Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracEnDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracEnDiffSameWExo->SetYTitle("#Sigma #it{E}_{diff}^{#it{w}}/#Sigma #it{E}_{same}^{#it{w}}");
  fhFracEnDiffSameWExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracEnDiffSameWExo);  
  
  fhFracNCellDiffSameWExo = new TH3F 
  ("hFracNCellDiffSameWExo","cell #it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracNCellDiffSameWExo->SetYTitle("#it{n}_{diff}^{#it{w}}/#it{n}_{same}^{#it{w}}");
  fhFracNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracNCellDiffSameWExo); 
  
  fhFracEnNCellDiffSameWExo = new TH3F 
  ("hFracEnNCellDiffSameWExo","cell (#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracEnNCellDiffSameWExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracEnNCellDiffSameWExo->SetYTitle("(#Sigma #it{E}_{diff}^{#it{w}}/#it{n}_{diff}^{#it{w}})/(#Sigma #it{E}_{same}^{#it{w}}/#it{n}_{same}^{#it{w}})");
  fhFracEnNCellDiffSameWExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracEnNCellDiffSameWExo); 
  
  fhFracEnDiffSame5Exo    = new TH3F 
  ("hFracEnDiffSame5Exo","cell #Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracEnDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracEnDiffSame5Exo->SetYTitle("#Sigma #it{E}_{diff}^{next}/#Sigma #it{E}_{same}^{next}");
  fhFracEnDiffSame5Exo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracEnDiffSame5Exo);  
  
  fhFracNCellDiffSame5Exo = new TH3F 
  ("hFracNCellDiffSame5Exo","cell #it{n}_{diff}^{next}/#it{n}_{same}^{next} vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 101,0,1.01,nexobins,exomin,exomax); 
  fhFracNCellDiffSame5Exo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhFracNCellDiffSame5Exo->SetYTitle("#it{n}_{diff}^{next}/#it{n}_{same}^{next}");
  fhFracNCellDiffSame5Exo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhFracNCellDiffSame5Exo); 
  
  fhFracEnNCellDiffSame5Exo = new TH3F 
  ("hFracEnNCellDiffSame5Exo","cell (#Sigma #it{E}_{diff}^{next}/#it{n}_{diff}^{next})/(#Sigma #it{E}_{same}^{next}/#it{n}_{same}^{next}) vs #it{E}_{cluster} vs #it{F}_{+}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, 200/2,0,2,nexobins,exomin,exomax); 
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
  fhFracEnNCellDiffSame->SetXTitle("(#Sigma #it{E}_{diff}/#it{n}_{diff})/(#Sigma #it{E}_{same}/#it{n}_{same})");
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
     17,0,17,nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
    fhNCellsSameDiffExo[i]->SetXTitle("#it{n}_{cells-same}");
    fhNCellsSameDiffExo[i]->SetYTitle("#it{n}_{cells-diff}");
    fhNCellsSameDiffExo[i]->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhNCellsSameDiffExo[i]);
    
    fhEnSameDiffExo[i] = new TH3F 
    (Form("hEnSameDiffExo_Ebin%d",i),
     Form("#Sigma #it{E}_{same}^{cells} vs #Sigma #it{E}_{diff}^{cells}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
     200, 0, 20, 200, 0, 20, nexobinsS,exominS,exomaxS); 
    fhEnSameDiffExo[i]->SetXTitle("#Sigma #it{E}_{same}^{cells} (GeV)");
    fhEnSameDiffExo[i]->SetYTitle("#Sigma #it{E}_{diff}^{cells} (GeV)");
    fhEnSameDiffExo[i]->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEnSameDiffExo[i]);
  }
  
  fhCellEnSameExo = new TH3F 
  ("hCellEnSameExo","#it{E}_{cluster} (GeV) vs #it{E}_{cell}^{same} vs #it{F}_{+}",
   nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
  fhCellEnSameExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnSameExo->SetYTitle("#it{E}_{cell}^{same} (GeV) ");
  fhCellEnSameExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnSameExo);    
  
  fhCellEnDiffExo = new TH3F 
  ("hCellEnDiffExo","#it{E}_{cluster} (GeV) vs #it{E}_{cell}^{diff} vs #it{F}_{+}",
   nptbins,ptmin,ptmax, 200,0,20, nexobins,exomin,exomax); 
  fhCellEnDiffExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnDiffExo->SetYTitle("#it{E}_{cell}^{diff} (GeV) ");
  fhCellEnDiffExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellEnDiffExo);     

  fhCellEnNCellWOpenTime = new TH3F 
  ("hCellEnNCellWOpenTime","#it{E}_{cluster} (GeV) vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}, no time cut",
   nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
  fhCellEnNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnNCellWOpenTime->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellWOpenTime);  

  fhCellEnNCellW = new TH3F 
  ("hCellEnNCellW","#it{E}_{cluster} (GeV) vs #it{E}_{cell} vs #it{n}_{cell}^{#it{w}}",
   nptbins,ptmin,ptmax, 200,0,20, nceclbins,nceclmin,nceclmax); 
  fhCellEnNCellW->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellEnNCellW->SetYTitle("#it{E}_{cell} (GeV) ");
  fhCellEnNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellEnNCellW);  

  fhCellMaxClusterEnRatioNCellWOpenTime = new TH3F 
  ("fhCellMaxClusterEnRatioNCellWOpenTime","#it{E}_{cluster} (GeV) vs #it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{n}_{cell}^{#it{w}}, no time cut",
   nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
  fhCellMaxClusterEnRatioNCellWOpenTime->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellWOpenTime->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellWOpenTime->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellMaxClusterEnRatioNCellWOpenTime);  
  
  fhCellMaxClusterEnRatioNCellW = new TH3F 
  ("fhCellMaxClusterEnRatioNCellW","#it{E}_{cluster} (GeV) vs #it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{n}_{cell}^{#it{w}}",
   nptbins,ptmin,ptmax, 101,0,1.01, nceclbins,nceclmin,nceclmax); 
  fhCellMaxClusterEnRatioNCellW->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellW->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioNCellW->SetZTitle("#it{n}_{cell}^{#it{w}}");
  outputContainer->Add(fhCellMaxClusterEnRatioNCellW);  
 
  fhCellMaxClusterEnRatioExo = new TH3F 
  ("fhCellMaxClusterEnRatioExo","#it{E}_{cluster} (GeV) vs #it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax, 101,0,1.01, nexobins,exomin,exomax); 
  fhCellMaxClusterEnRatioExo->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioExo->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster} (GeV) ");
  fhCellMaxClusterEnRatioExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellMaxClusterEnRatioExo);  
  
  // Cluster acceptance
  //
  fhEtaPhiGridExoEnCut  = new TH3F 
  ("hEtaPhiGridExoEnCut",
   Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1",fEMinForExo),
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
  fhEtaPhiGridExoEnCut->SetXTitle("column-#eta ");
  fhEtaPhiGridExoEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridExoEnCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhEtaPhiGridExoEnCut);
  
  fhEtaPhiGridExoEnCutSameFracCut = new TH3F 
  ("hEtaPhiGridExoEnCutSameFracCut",
   Form("colum (#eta) vs row (#varphi) vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1, #it{n}_{cells}^{same}/(#it{n}_{cells}-1)=1",fEMinForExo),
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nexobinsS,exominS,exomaxS); 
  fhEtaPhiGridExoEnCutSameFracCut->SetXTitle("column-#eta ");
  fhEtaPhiGridExoEnCutSameFracCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridExoEnCutSameFracCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhEtaPhiGridExoEnCutSameFracCut);
  
  fhEtaPhiGridEnExoCut  = new TH3F 
  ("hEtaPhiGridEnExoCut", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{F}_{+} > %2.2f",fExoCut),
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
  fhEtaPhiGridEnExoCut->SetXTitle("column-#eta ");
  fhEtaPhiGridEnExoCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridEnExoCut->SetZTitle("#it{E} (GeV)");
  outputContainer->Add(fhEtaPhiGridEnExoCut);
  
  fhEtaPhiGridEnHighNCells  = new TH3F 
  ("hEtaPhiGridEnHighNCells", Form("colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}^{#it{w}} > %d",fNCellHighCut),
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins/2,ptmin,ptmax); 
  fhEtaPhiGridEnHighNCells->SetXTitle("column-#eta ");
  fhEtaPhiGridEnHighNCells->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridEnHighNCells->SetZTitle("#it{E} (GeV)");
  outputContainer->Add(fhEtaPhiGridEnHighNCells);
 
  fhEtaPhiGridNCellEnCut  = new TH3F 
  ("hEtaPhiGridNCellEnCut",
   Form("colum (#eta) vs row (#varphi) vs #it{n}_{cells}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1",fEMinForExo),
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nceclbins,nceclmin,nceclmax); 
  fhEtaPhiGridNCellEnCut->SetXTitle("column-#eta ");
  fhEtaPhiGridNCellEnCut->SetYTitle("row-#varphi (rad)");
  fhEtaPhiGridNCellEnCut->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhEtaPhiGridNCellEnCut);
  
  if ( fFill1CellHisto )
  {
    fhEtaPhiGridEn1Cell  = new TH3F 
    ("hEtaPhiGridEn1Cell", "colum (#eta) vs row (#varphi) vs #it{E}, #it{n}_{cells}=1",
     ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,nptbins,ptmin,ptmax); 
    fhEtaPhiGridEn1Cell->SetXTitle("column-#eta ");
    fhEtaPhiGridEn1Cell->SetYTitle("row-#varphi (rad)");
    fhEtaPhiGridEn1Cell->SetZTitle("#it{E} (GeV)");
    outputContainer->Add(fhEtaPhiGridEn1Cell);
  }
  
  // Timing and energy
  fhTimeEnergyExo  = new TH3F 
  ("hTimeEnergyExo","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
  fhTimeEnergyExo->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyExo->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeEnergyExo);
  
  if ( fFill1CellHisto )
  {
    fhTimeEnergy1Cell  = new TH2F 
    ("hTimeEnergy1Cell","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}=1",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhTimeEnergy1Cell->SetXTitle("#it{E} (GeV)");
    fhTimeEnergy1Cell->SetYTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhTimeEnergy1Cell);
  }
  
  fhTimeDiffClusCellExo  = new TH3F 
  ("hTimeDiffClusCellExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
  fhTimeDiffClusCellExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellExo);
  
  fhTimeDiffAmpClusCellExo  = new TH3F 
  ("hTimeDiffAmpClusCellExo",
   Form("#it{E}_{cell i} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1"),
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
  fhTimeDiffAmpClusCellExo->SetXTitle("#it{E}_{cell i} (GeV)");
  fhTimeDiffAmpClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffAmpClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffAmpClusCellExo);
  
  fhTimeEnergyM02  = new TH3F 
  ("hTimeEnergyM02","#it{E}_{cluster} vs #it{t}_{cluster} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, 100,0,0.5); 
  fhTimeEnergyM02->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyM02->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeEnergyM02);
  
  fhTimeDiffClusCellM02  = new TH3F 
  ("hTimeDiffClusCellM02","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, 100,0,0.5); 
  fhTimeDiffClusCellM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellM02->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeDiffClusCellM02);
  
  fhTimeEnergyNCells  = new TH3F 
  ("hTimeEnergyNCells","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{n}_{cells}",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nceclbins,nceclmin,nceclmax); 
  fhTimeEnergyNCells->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyNCells->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyNCells->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhTimeEnergyNCells);

  fhTimeEnergyNCellsW  = new TH3F 
  ("hTimeEnergyNCellsW","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{n}_{cells} with #it{w} > 0",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nceclbins,nceclmin,nceclmax); 
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
  fhM02EnergyNCell  = new TH3F 
  ("hM02EnergyNCell","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{n}_{cells}",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nceclbins,nceclmin,nceclmax); 
  fhM02EnergyNCell->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyNCell->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyNCell->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhM02EnergyNCell); 
  
  fhM02EnergyAllSameTCard  = new TH2F 
  ("hM02EnergyNCellAllSameTCard","#sigma^{2}_{long} vs #it{E}_{cluster}, #it{n}_{cells} = #it{n}_{cells-same}",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhM02EnergyAllSameTCard->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyAllSameTCard->SetYTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhM02EnergyAllSameTCard); 
  
  fhM02EnergyExo  = new TH3F 
  ("hM02EnergyExo","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nexobinsS,exominS,exomaxS); 
  fhM02EnergyExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyExo->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM02EnergyExo); 

  fhM02EnergyExoZoomIn  = new TH3F 
  ("hM02EnergyExoZoomIn","#sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax,100,0,0.5,nexobinsS,exominS,exomaxS); 
  fhM02EnergyExoZoomIn->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyExoZoomIn->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyExoZoomIn->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM02EnergyExoZoomIn); 
  
  fhM20EnergyExoM02MinCut  = new TH3F 
  ("hM20EnergyExoM02MinCut","#sigma^{2}_{short} vs #it{E}_{cluster} vs #it{F}_{+}, #sigma^{2}_{long} > 0.1",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/2, nexobinsS,exominS,exomaxS); 
  fhM20EnergyExoM02MinCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM20EnergyExoM02MinCut->SetYTitle("#sigma^{2}_{short}");
  fhM20EnergyExoM02MinCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM20EnergyExoM02MinCut);  
  
  for(Int_t i = 0; i < fgkNEBins-1; i++) 
  {
    fhM02ExoNCells[i] = new TH3F 
    (Form("hM02ExoNCells_Ebin%d",i),
     Form("#sigma^{2}_{long} vs #it{F}_{+} vs #it{n}_{cells}, %2.1f < #it{E} < %2.1f GeV",fEnergyBins[i],fEnergyBins[i+1]),
     100,0,0.5,nexobinsS,exominS,exomaxS,nceclbins,nceclmin,nceclmax); 
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
     17,-8.5,8.5,17,-8.5,8.5,totalSM,fFirstModule-0.5,fLastModule+0.5); 
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
       17,-8.5,8.5,17,-8.5,8.5,nexobinsS,exominS,exomaxS); 
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
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDEtaNegExo->SetYTitle("d#eta");
    fhTrackMatchedDEtaNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDEtaNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDPhiNegExo  = new TH3F
    ("hTrackMatchedDPhiNegExo","d#varphi of cluster-negative track vs cluster energy vs #it{F}_{+}",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDPhiNegExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaDPhiNegExo  = new TH3F
    ("hTrackMatchedDEtaDPhiNegExo",
     Form("d#eta vs d#varphi of cluster- negative track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDEtaDPhiNegExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiNegExo->SetXTitle("d#eta");
    fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaPosExo  = new TH3F
    ("hTrackMatchedDEtaPosExo","d#eta of cluster-positive track vs cluster energy vs #it{F}_{+}",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDEtaPosExo->SetYTitle("d#eta");
    fhTrackMatchedDEtaPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDEtaPosExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDPhiPosExo  = new TH3F
    ("hTrackMatchedDPhiPosExo","d#varphi of cluster-positive track vs cluster energy vs #it{F}_{+}",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDPhiPosExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
    fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhTrackMatchedDEtaDPhiPosExo  = new TH3F
    ("hTrackMatchedDEtaDPhiPosExo",
     Form("d#eta vs d#varphi of cluster-positive track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
    fhTrackMatchedDEtaDPhiPosExo->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiPosExo->SetXTitle("d#eta");
    fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");
    
    fhEOverPExo = new TH3F
    ("hEOverPExo",
     "Track matches #it{E}/#it{p} vs #it{F}_{+}",
     nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax, nexobinsS,exominS,exomaxS);
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
     nptbins,ptmin,ptmax/2, nexobins,exomin,exomax); 
    fhCellExoAmp->SetXTitle("#it{E}_{cell} (GeV) ");
    fhCellExoAmp->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoAmp);    
    
    fhCellExoAmpTime = new TH3F 
    ("hCellExoAmpTime","Cell #it{F}_{+} vs #it{E}_{cell} vs time",
     nptbins,ptmin,ptmax/2, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
    fhCellExoAmpTime->SetXTitle("#it{E}_{cell} (GeV) ");
    fhCellExoAmpTime->SetYTitle("#it{t}_{cell} (ns)");
    fhCellExoAmpTime->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoAmpTime);    
    
    fhCellExoGrid    = new TH3F 
    ("hCellExoGrid",
     Form("Cell hits row-column vs #it{F}_{+} for #it{E}_{cell} > %2.1f",fEMinForExo), 
     ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, nexobinsS,exominS,exomaxS); 
    fhCellExoGrid->SetYTitle("row (phi direction)");
    fhCellExoGrid->SetXTitle("column (eta direction)");
    fhCellExoGrid->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoGrid);
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

