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
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TH3F.h>
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaCalorimeterQA.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliMCAnalysisUtils.h"

// --- Detectors --- 
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaCalorimeterQA) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters.
//__________________________________________
AliAnaCalorimeterQA::AliAnaCalorimeterQA() :
AliAnaCaloTrackCorrBaseClass(),  

// Switches
fFillAllCellTimeHisto(kTRUE),
fFillAllPosHisto(kFALSE),              fFillAllPosHisto2(kFALSE),
fFillAllTH3(kFALSE),                   
fFillAllTMHisto(kTRUE),                
fFillAllPi0Histo(kTRUE),               fFillInvMassOpenAngle(kFALSE),               
fFillPi0PairDiffTime(kFALSE),          fFillInvMassInEMCALWithPHOSDCalAcc(kFALSE), 
fCorrelate(kTRUE),                     fStudyBadClusters(kFALSE),               
fStudyClustersAsymmetry(kFALSE),       fStudyExotic(kFALSE),
fStudyWeight(kFALSE),

// Parameters and cuts
fNModules(12),                         fNRCU(2),
fNMaxCols(48),                         fNMaxRows(24),  
fTimeCutMin(-10000),                   fTimeCutMax(10000),
fCellAmpMin(0),                        
fEMCALCellAmpMin(0),                   fPHOSCellAmpMin(0),                    
fEMCALClusterNCellMin(0),              fPHOSClusterNCellMin(0),                    

// Invariant mass
fInvMassMinECut(0),                    fInvMassMaxECut(0),                    
fInvMassMinM02Cut(0),                  fInvMassMaxM02Cut(0),                    
fInvMassMaxOpenAngle(0),               fInvMassMaxTimeDifference(0),

// Exotic
fExoNECrossCuts(0),                    fExoECrossCuts(),
fExoNDTimeCuts(0),                     fExoDTimeCuts(),    

fClusterMomentum(),                    fClusterMomentum2(),
fPrimaryMomentum(),
fConstantTimeShift(0),
// Histograms
fhE(0),                                fhPt(0),                                
fhPhi(0),                              fhEta(0),                               
fhEtaPhi(0),                           fhEtaPhiE(0),
fhECharged(0),                         fhPtCharged(0),             
fhPhiCharged(0),                       fhEtaCharged(0),                        
fhEtaPhiCharged(0),                    fhEtaPhiECharged(0), 

// Invariant mass
fhIM(0),                               fhIMSame(0),             fhIMDiff(0),                              
fhIMDCAL(0),                           fhIMDCALSame(0),         fhIMDCALDiff(0),
fhIMDCALPHOS(0),                       fhIMDCALPHOSSame(0), 
fhIMEMCALPHOS(0),                      fhIMEMCALPHOSSame(0), 
fhAsym(0), 
fhOpAngle(0),                          fhIMvsOpAngle(0),
fhNCellsPerCluster(0),                 fhNCellsPerClusterNoCut(0),             fhNClusters(0),

// Timing
fhClusterTimeEnergy(0),                fhCellTimeSpreadRespectToCellMax(0),  
fhCellIdCellLargeTimeSpread(0),        fhClusterPairDiffTimeE(0),              fhClusterPairDiffTimeESameMod(0),
fhClusterMaxCellCloseCellRatio(0),     fhClusterMaxCellCloseCellDiff(0), 
fhClusterMaxCellDiff(0),               fhClusterMaxCellDiffNoCut(0), 
//fhClusterMaxCellDiffAverageTime(0),    fhClusterMaxCellDiffWeightedTime(0),    
fhClusterMaxCellECross(0),
fhLambda0(0),                          fhLambda1(0),                          // fhDispersion(0),

// bad clusters
fhBadClusterEnergy(0),                 fhBadClusterTimeEnergy(0),              fhBadClusterEtaPhi(0),            
fhBadClusterPairDiffTimeE(0),          fhBadCellTimeSpreadRespectToCellMax(0), 
fhBadClusterMaxCellCloseCellRatio(0),  fhBadClusterMaxCellCloseCellDiff(0),    fhBadClusterMaxCellDiff(0),
fhBadClusterMaxCellDiffAverageTime(0), fhBadClusterMaxCellDiffWeightedTime(0),
fhBadClusterMaxCellECross(0),
fhBadClusterLambda0(0),                fhBadClusterLambda1(0),
fhBadClusterDeltaIEtaDeltaIPhiE0(0),   fhBadClusterDeltaIEtaDeltaIPhiE2(0),          
fhBadClusterDeltaIEtaDeltaIPhiE6(0),   fhBadClusterDeltaIA(0), 

// Position
fhRNCells(0),                          fhXNCells(0),               
fhYNCells(0),                          fhZNCells(0),
fhRE(0),                               fhXE(0),                    
fhYE(0),                               fhZE(0),    
fhXYZ(0),
fhRCellE(0),                           fhXCellE(0),                
fhYCellE(0),                           fhZCellE(0),
fhXYZCell(0),
fhDeltaCellClusterRNCells(0),          fhDeltaCellClusterXNCells(0),
fhDeltaCellClusterYNCells(0),          fhDeltaCellClusterZNCells(0),
fhDeltaCellClusterRE(0),               fhDeltaCellClusterXE(0),     
fhDeltaCellClusterYE(0),               fhDeltaCellClusterZE(0),

// Cells
fhNCells(0),                           fhNCellsCutAmpMin(0),
fhAmplitude(0),                        fhAmpId(0),                             
fhEtaPhiAmpCell(0),                    fhEtaPhiCell(0),
fhTime(0),                             //fhTimeVz(0),
fhTimeId(0),                           fhTimeAmp(0),
fhAmpIdLowGain(0),                     fhTimeIdLowGain(0),                     fhTimeAmpLowGain(0),

fhCellECross(0),

fhEMCALPHOSCorrNClusters(0),           fhEMCALPHOSCorrEClusters(0),     
fhEMCALPHOSCorrNCells(0),              fhEMCALPHOSCorrECells(0),
fhEMCALDCALCorrNClusters(0),           fhEMCALDCALCorrEClusters(0),     
fhEMCALDCALCorrNCells(0),              fhEMCALDCALCorrECells(0),
fhDCALPHOSCorrNClusters(0),            fhDCALPHOSCorrEClusters(0),     
fhDCALPHOSCorrNCells(0),               fhDCALPHOSCorrECells(0),
fhCaloV0SCorrNClusters(0),             fhCaloV0SCorrEClusters(0),              
fhCaloV0SCorrNCells(0),                fhCaloV0SCorrECells(0),
fhCaloV0MCorrNClusters(0),             fhCaloV0MCorrEClusters(0),  
fhCaloV0MCorrNCells(0),                fhCaloV0MCorrECells(0),
fhCaloTrackMCorrNClusters(0),          fhCaloTrackMCorrEClusters(0), 
fhCaloTrackMCorrNCells(0),             fhCaloTrackMCorrECells(0),
fhCaloCenNClusters(0),                 fhCaloCenEClusters(0),
fhCaloCenNCells(0),                    fhCaloCenECells(0),
fhCaloEvPNClusters(0),                 fhCaloEvPEClusters(0),
fhCaloEvPNCells(0),                    fhCaloEvPECells(0),

// Super-Module dependent histograms
fhEMod(0),                             fhAmpMod(0),                            
fhEWeirdMod(0),                        fhAmpWeirdMod(0),                            
fhTimeMod(0),  
fhNClustersMod(0),                     fhNCellsMod(0),
fhNCellsSumAmpPerMod(0),               fhNClustersSumEnergyPerMod(0),
fhNCellsPerClusterMod(0),              fhNCellsPerClusterModNoCut(0), 
fhNCellsPerClusterWeirdMod(0),         fhNCellsPerClusterWeirdModNoCut(0), 

fhGridCells(0),                        fhGridCellsE(0),                        fhGridCellsTime(0),
fhGridCellsLowGain(0),                 fhGridCellsELowGain(0),                 fhGridCellsTimeLowGain(0),
fhTimeAmpPerRCU(0),                    fhIMMod(0),

// Weight studies
fhECellClusterRatio(0),                fhECellClusterLogRatio(0),                 
fhEMaxCellClusterRatio(0),             fhEMaxCellClusterLogRatio(0),                
fhECellTotalRatio(0),                  fhECellTotalLogRatio(0),
fhECellTotalRatioMod(0),               fhECellTotalLogRatioMod(0),

fhExoL0ECross(0),                      fhExoL1ECross(0),

// MC and reco
fhRecoMCE(),                           fhRecoMCPhi(),                          fhRecoMCEta(), 
fhRecoMCDeltaE(),                      fhRecoMCRatioE(),                      
fhRecoMCDeltaPhi(),                    fhRecoMCDeltaEta(),               

// MC only
fhGenMCE(),                            fhGenMCPt(),                            fhGenMCEtaPhi(),
fhGenMCAccE(),                         fhGenMCAccPt(),                         fhGenMCAccEtaPhi(),

// Matched MC
fhEMVxyz(0),                           fhEMR(0),                   
fhHaVxyz(0),                           fhHaR(0),
fh1EOverP(0),                          fh2dR(0),                   
fh2EledEdx(0),                         fh2MatchdEdx(0),
fh1EOverPR02(0),                       fh1EleEOverP(0),
fhMCEle1EOverP(0),                     fhMCEle1dR(0),                          fhMCEle2MatchdEdx(0),
fhMCChHad1EOverP(0),                   fhMCChHad1dR(0),                        fhMCChHad2MatchdEdx(0),
fhMCNeutral1EOverP(0),                 fhMCNeutral1dR(0),                      fhMCNeutral2MatchdEdx(0),        
fhMCEle1EOverPR02(0),                  fhMCChHad1EOverPR02(0),                 fhMCNeutral1EOverPR02(0),
fhMCEle1EleEOverP(0),                  fhMCChHad1EleEOverP(0),                 fhMCNeutral1EleEOverP(0),
fhTrackMatchedDEtaNeg(0),              fhTrackMatchedDPhiNeg(0),               fhTrackMatchedDEtaDPhiNeg(0),
fhTrackMatchedDEtaPos(0),              fhTrackMatchedDPhiPos(0),               fhTrackMatchedDEtaDPhiPos(0),
fhTrackMatchedDEtaNegMod(0),           fhTrackMatchedDPhiNegMod(0),             
fhTrackMatchedDEtaPosMod(0),           fhTrackMatchedDPhiPosMod(0)
{
  // Weight studies
  for(Int_t i =0; i < 12; i++)
  {
    for(Int_t j = 0; j < 4; j++)
    {    
      for(Int_t k = 0; k < 3; k++)
      {
        fhLambda0ForW0AndCellCuts    [i][j][k] = 0;
//      fhLambda1ForW0AndCellCuts    [i][j][k] = 0;
        fhLambda0ForW0AndCellCutsEta0[i][j][k] = 0;
      }
    }
    for(Int_t j = 0; j < 5; j++)
    {
      fhLambda0ForW0MC[i][j] = 0;
//    fhLambda1ForW0MC[i][j] = 0;
    }
  }
  
  // Cluster size
  fhDeltaIEtaDeltaIPhiE0[0] = 0 ;         fhDeltaIEtaDeltaIPhiE2[0] = 0;          fhDeltaIEtaDeltaIPhiE6[0] = 0; 
  fhDeltaIEtaDeltaIPhiE0[1] = 0 ;         fhDeltaIEtaDeltaIPhiE2[1] = 0;          fhDeltaIEtaDeltaIPhiE6[1] = 0; 
  fhDeltaIA[0]              = 0 ;         fhDeltaIAL0[0]            = 0;          fhDeltaIAL1[0]            = 0;
  fhDeltaIA[1]              = 0 ;         fhDeltaIAL0[1]            = 0;          fhDeltaIAL1[1]            = 0;                         
  fhDeltaIANCells[0]        = 0 ;         fhDeltaIANCells[1]        = 0;
  fhDeltaIAMC[0]            = 0 ;         fhDeltaIAMC[1]            = 0;
  fhDeltaIAMC[2]            = 0 ;         fhDeltaIAMC[3]            = 0;
  
  // Exotic
  for (Int_t ie = 0; ie < 10 ; ie++) 
  {
    fhExoDTime[ie] = 0;
    for (Int_t idt = 0; idt < 5 ; idt++) 
    {
      fhExoNCell    [ie][idt] = 0;
      fhExoL0       [ie][idt] = 0;
      fhExoL1       [ie][idt] = 0;
      fhExoECross   [ie][idt] = 0;
      fhExoTime     [ie][idt] = 0;
      fhExoL0NCell  [ie][idt] = 0;
      fhExoL1NCell  [ie][idt] = 0;
    } 
  }
  
  // MC
  
  for(Int_t i = 0; i < 7; i++)
  {
    fhRecoMCE[i][0]         = 0; fhRecoMCE[i][1]        = 0;  
    fhRecoMCPhi[i][0]       = 0; fhRecoMCPhi[i][1]      = 0;  
    fhRecoMCEta[i][0]       = 0; fhRecoMCEta[i][1]      = 0;  
    fhRecoMCDeltaE[i][0]    = 0; fhRecoMCDeltaE[i][1]   = 0;  
    fhRecoMCRatioE[i][0]    = 0; fhRecoMCRatioE[i][1]   = 0;  
    fhRecoMCDeltaPhi[i][0]  = 0; fhRecoMCDeltaPhi[i][1] = 0;  
    fhRecoMCDeltaEta[i][0]  = 0; fhRecoMCDeltaEta[i][1] = 0;
  }
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhGenMCE[i]         = 0;
    fhGenMCPt[i]        = 0;
    fhGenMCEtaPhi[i]    = 0;
    fhGenMCAccE[i]      = 0;
    fhGenMCAccPt[i]     = 0;
    fhGenMCAccEtaPhi[i] = 0;
  }
  
  InitParameters();
}

//______________________________________________________________________________________________________________________
/// Fill histograms with parameters of clusters declared bad when calling the method.
/// \param clus: calorimeter cluster pointer.
/// \param caloClusters: full list of clusters, used for comparison of bad and good clusters timing.
/// \param cells: list of cells, needed to get the cell with highest energy in cluster.
/// \param absIdMax: id number of cell with highest energy in the cluster.
/// \param maxCellFraction: ratio E_cell_max/ E_cluster.
/// \param eCrossFrac: exoticity fraction.
/// \param tmax: time of highest energy cell in cluster.
///
//______________________________________________________________________________________________________________________
void AliAnaCalorimeterQA::BadClusterHistograms(AliVCluster* clus, const TObjArray *caloClusters, AliVCaloCells * cells,
                                               Int_t absIdMax, Double_t maxCellFraction, Float_t eCrossFrac,
                                               Double_t tmax)
{
  //  printf("AliAnaCalorimeterQA::BadClusterHistograms() - Event %d - Calorimeter %s \n \t  E %f, n cells %d, max cell absId %d, maxCellFrac %f\n",
  //         GetReader()->GetEventNumber(), GetCalorimeterString().Data(), 
  //         clus->E(),clus->GetNCells(),absIdMax,maxCellFraction);
      
  Double_t tof    = clus->GetTOF()*1.e9;
  if(tof > 400) tof-=fConstantTimeShift;
  
  Float_t  energy = clus->E();
  
  fhBadClusterEnergy       ->Fill(energy,                  GetEventWeight());
  fhBadClusterTimeEnergy   ->Fill(energy, tof            , GetEventWeight());
  fhBadClusterMaxCellDiff  ->Fill(energy, maxCellFraction, GetEventWeight());
  fhBadClusterMaxCellECross->Fill(energy, eCrossFrac     , GetEventWeight());

  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi += TMath::TwoPi();
  
  if(energy > 0.5) fhBadClusterEtaPhi->Fill(fClusterMomentum.Eta(), phi, GetEventWeight());

  fhBadClusterLambda0->Fill(energy, clus->GetM02());
  fhBadClusterLambda1->Fill(energy, clus->GetM20());
  
  if(fStudyClustersAsymmetry) ClusterAsymmetryHistograms(clus,absIdMax,kFALSE);
  
  // Clusters in event time difference bad minus good
  
  for(Int_t iclus2 = 0; iclus2 < caloClusters->GetEntriesFast(); iclus2++ )
  {
    AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(iclus2);
    
    if(clus->GetID() == clus2->GetID()) continue;
    
    Float_t maxCellFraction2 = 0.;
    Int_t absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(cells, clus2,maxCellFraction2);
    
    if(IsGoodCluster(absIdMax2, clus->GetM02(), clus->GetNCells(), cells) &&  clus2->GetM02() > 0.1 )
    {
      Double_t tof2   = clus2->GetTOF()*1.e9;     
      if(tof2>400) tof2-=fConstantTimeShift;
      
      fhBadClusterPairDiffTimeE  ->Fill(clus->E(), (tof-tof2), GetEventWeight());
    }
  } // loop
  
  // Max cell compared to other cells in cluster
  if(fFillAllCellTimeHisto) 
  {
    // Get some time averages
    Double_t timeAverages[2] = {0.,0.};
    CalculateAverageTime(clus, cells, timeAverages);

    fhBadClusterMaxCellDiffAverageTime ->Fill(clus->E(), tmax-timeAverages[0], GetEventWeight());
    fhBadClusterMaxCellDiffWeightedTime->Fill(clus->E(), tmax-timeAverages[1], GetEventWeight());
  }           
  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) 
  {
    Int_t absId  = clus->GetCellsAbsId()[ipos]; 
    if(absId!=absIdMax && cells->GetCellAmplitude(absIdMax) > 0.01)
    {
      Float_t frac = cells->GetCellAmplitude(absId)/cells->GetCellAmplitude(absIdMax);
      
      fhBadClusterMaxCellCloseCellRatio->Fill(clus->E(), frac, GetEventWeight());
      fhBadClusterMaxCellCloseCellDiff ->Fill(clus->E(), cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId), GetEventWeight());
      
      if(fFillAllCellTimeHisto) 
      {
        Double_t time  = cells->GetCellTime(absId);
        GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
        
        Float_t diff = (tmax-(time*1e9-fConstantTimeShift));
        fhBadCellTimeSpreadRespectToCellMax->Fill(clus->E(), diff, GetEventWeight());
        
      } 
    } // Not max
  } // loop
}

//______________________________________________________________________
/// Calculate time averages and weights in clusters.
/// \param clus: calorimeter cluster pointer.
/// \param cells: list of cells, needed to get the cell with highest energy in cluster.
/// \param timeAverages: return average weighted time in cluster.
//______________________________________________________________________
void AliAnaCalorimeterQA::CalculateAverageTime(AliVCluster *clus,
                                               AliVCaloCells* cells,  
                                               Double_t timeAverages[2])
{
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0, amp = 0;
//Int_t    absIdMax = -1;
    
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) 
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    // Recalibrate cell energy if needed
    amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,GetCalorimeter(), id);
    
    energy    += amp;
    
    if(amp> ampMax) 
    {
      ampMax   = amp;
//    absIdMax = id;
    }
  } // energy loop       
  
  // Calculate average time of cells in cluster and weighted average
  Double_t aTime  = 0; 
  Double_t wTime  = 0;
  Float_t  wTot   = 0;
  Double_t time   = 0;
  Int_t    id     =-1;
  Double_t w      = 0;
  Int_t    ncells = clus->GetNCells();
  
  for (Int_t ipos = 0; ipos < ncells; ipos++) 
  {
    id   = clus ->GetCellsAbsId()[ipos];
    amp  = cells->GetCellAmplitude(id);
    time = cells->GetCellTime(id);
    
    // Recalibrate energy and time
    GetCaloUtils()->RecalibrateCellAmplitude(amp , GetCalorimeter(), id);    
    GetCaloUtils()->RecalibrateCellTime     (time, GetCalorimeter(), id, GetReader()->GetInputEvent()->GetBunchCrossNumber());

    w      = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(cells->GetCellAmplitude(id),energy);
    aTime += time*1e9;
    wTime += time*1e9 * w;
    wTot  += w;
  }        
  
  if(ncells > 0) aTime /= ncells;
  else           aTime  = 0;
  
  if(wTot   > 0) wTime /= wTot;
  else           wTime  = 0;

  timeAverages[0] = aTime;        
  timeAverages[1] = wTime;
}

//____________________________________________________________
/// Fill histograms related to cells only.
//____________________________________________________________
void AliAnaCalorimeterQA::CellHistograms(AliVCaloCells *cells)
{
  Int_t ncells = cells->GetNumberOfCells();
  if( ncells > 0 ) fhNCells->Fill(ncells, GetEventWeight()) ; // Not ok for PHOS with CPV

  Int_t   ncellsCut = 0;
  Float_t ecellsCut = 0;
  
  AliDebug(1,Form("%s cell entries %d", GetCalorimeterString().Data(), ncells));
  
  // Init arrays and used variables
  Int_t   *nCellsInModule = new Int_t  [fNModules];
  Float_t *eCellsInModule = new Float_t[fNModules];
  
  for(Int_t imod = 0; imod < fNModules; imod++ )
  {
    nCellsInModule[imod] = 0 ;
    eCellsInModule[imod] = 0.;
  }
  
  Int_t    icol   = -1, icolAbs = -1;
  Int_t    irow   = -1, irowAbs = -1;
  Int_t    iRCU   = -1;
  Float_t  amp    = 0.;
  Double_t time   = 0.;
  Int_t    id     = -1;
  Bool_t   highG  = kFALSE;
  Float_t  recalF = 1.;  
  Int_t    bc     = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++)
  {
    if ( cells->GetCellNumber(iCell) < 0 ) continue; // CPV 
    
    AliDebug(2,Form("Cell : amp %f, absId %d", cells->GetAmplitude(iCell), cells->GetCellNumber(iCell)));
   
    Int_t nModule = GetModuleNumberCellIndexesAbsCaloMap(cells->GetCellNumber(iCell),GetCalorimeter(), 
                                                         icol   , irow, iRCU,
                                                         icolAbs, irowAbs    );
    
    AliDebug(2,Form("\t module %d, column %d (%d), row %d (%d)", nModule,icolAbs,icol,irowAbs,irow));
    
    if(nModule < fNModules) 
    {	
      //Check if the cell is a bad channel
      if(GetCaloUtils()->IsBadChannelsRemovalSwitchedOn())
      {
        if(GetCalorimeter()==kEMCAL)
        {
          if(GetCaloUtils()->GetEMCALChannelStatus(nModule,icol,irow)) continue;
        }
        else 
        {
          if(GetCaloUtils()->GetPHOSChannelStatus(nModule,icol,irow) ) continue;
        }
      } // use bad channel map
      
      amp     = cells->GetAmplitude(iCell)*recalF;
      time    = cells->GetTime(iCell);
      id      = cells->GetCellNumber(iCell);
      highG   = cells->GetCellHighGain(id);
      if(IsDataMC()) highG = kTRUE; // MC does not distinguish High and Low, put them all in high
      
      // Amplitude recalibration if set
      GetCaloUtils()->RecalibrateCellAmplitude(amp,  GetCalorimeter(), id);
      
      // Time recalibration if set
      GetCaloUtils()->RecalibrateCellTime     (time, GetCalorimeter(), id, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
      
      // Transform time to ns
      time *= 1.0e9;
      time-=fConstantTimeShift;

      if(time < fTimeCutMin || time > fTimeCutMax)
      {
        AliDebug(1,Form("Remove cell with Time %f",time));
        continue;
      }
      
      // Remove exotic cells, defined only for EMCAL
      if(GetCalorimeter()==kEMCAL && 
         GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(id, cells, bc)) continue;
      
      fhAmplitude->Fill(amp,          GetEventWeight());
      fhAmpId    ->Fill(amp, id     , GetEventWeight());
      fhAmpMod   ->Fill(amp, nModule, GetEventWeight());
      fhAmpWeirdMod->Fill(amp, nModule, GetEventWeight());
      if(!highG) fhAmpIdLowGain->Fill(amp, id, GetEventWeight());
        
      // E cross for exotic cells
      if(amp > 0.05)
      {
        fhCellECross->Fill(amp, 1-GetECross(id,cells)/amp, GetEventWeight());
        ecellsCut+=amp ;
        if(fStudyWeight) eCellsInModule[nModule]+=amp ;
      }
      
      if ( amp > fCellAmpMin )
      {
        ncellsCut++    ;
        nCellsInModule[nModule]++    ;

        if(!fStudyWeight) eCellsInModule[nModule]+=amp ;
        
        fhGridCells ->Fill(icolAbs, irowAbs, GetEventWeight());
        fhGridCellsE->Fill(icolAbs, irowAbs, amp             );
        
        if(!highG)
        {
          fhGridCellsLowGain ->Fill(icolAbs, irowAbs, GetEventWeight());
          fhGridCellsELowGain->Fill(icolAbs, irowAbs, amp             );
        }
        
        if(fFillAllCellTimeHisto)
        {
          //printf("%s: time %g\n",GetCalorimeterString().Data(), time);
          
//          Double_t v[3] = {0,0,0}; //vertex ;
//          GetReader()->GetVertex(v);          
//          if(amp > 0.5) fhTimeVz   ->Fill(TMath::Abs(v[2]), time, GetEventWeight());
          
          fhTime   ->Fill(time,       GetEventWeight());
          fhTimeId ->Fill(time, id  , GetEventWeight());
          fhTimeAmp->Fill(amp , time, GetEventWeight());
            
          fhGridCellsTime->Fill(icolAbs, irowAbs, time);
          if(!highG) fhGridCellsTimeLowGain->Fill(icolAbs, irowAbs, time);
            
          fhTimeMod->Fill(time, nModule, GetEventWeight());
          fhTimeAmpPerRCU[nModule*fNRCU+iRCU]->Fill(amp, time, GetEventWeight());
          
          if(!highG)
          {
            fhTimeIdLowGain ->Fill(time, id  , GetEventWeight());
            fhTimeAmpLowGain->Fill(amp , time, GetEventWeight());
          }
        }
      }
      
      // Get Eta-Phi position of Cell
      if(fFillAllPosHisto)
      {
        if ( GetCalorimeter() == kEMCAL && GetCaloUtils()->IsEMCALGeoMatrixSet() )
        {
          Float_t celleta = 0.;
          Float_t cellphi = 0.;
          GetEMCALGeometry()->EtaPhiFromIndex(id, celleta, cellphi); 

          if ( cellphi < 0 ) cellphi+=TMath::TwoPi();
          
          if(fFillAllTH3)
            fhEtaPhiAmpCell->Fill(celleta, cellphi, amp, GetEventWeight());
          else
            fhEtaPhiCell   ->Fill(celleta, cellphi,      GetEventWeight());
          
          Double_t cellpos[] = {0, 0, 0};
          GetEMCALGeometry()->GetGlobal(id, cellpos);
            
          fhXCellE->Fill(cellpos[0], amp, GetEventWeight())  ;
          fhYCellE->Fill(cellpos[1], amp, GetEventWeight())  ;
          fhZCellE->Fill(cellpos[2], amp, GetEventWeight())  ;
            
          Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]);//+cellpos[2]*cellpos[2]);
          fhRCellE ->Fill(rcell, amp, GetEventWeight())  ;
            
          fhXYZCell->Fill(cellpos[0], cellpos[1], cellpos[2], GetEventWeight())  ;
        } // EMCAL Cells
        else if ( GetCalorimeter() == kPHOS && GetCaloUtils()->IsPHOSGeoMatrixSet() )
        {
          TVector3 xyz;
          Int_t relId[4], module;
          Float_t xCell, zCell;
          
          GetPHOSGeometry()->AbsToRelNumbering(id,relId);
          module = relId[0];
          GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
          GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
          
          Float_t rcell = TMath::Sqrt(xyz.X()*xyz.X()+xyz.Y()*xyz.Y());
          
          fhXCellE ->Fill(xyz.X(), amp, GetEventWeight())  ;
          fhYCellE ->Fill(xyz.Y(), amp, GetEventWeight())  ;
          fhZCellE ->Fill(xyz.Z(), amp, GetEventWeight())  ;
          fhRCellE ->Fill(rcell  , amp, GetEventWeight())  ;
            
          fhXYZCell->Fill(xyz.X(), xyz.Y(), xyz.Z(), GetEventWeight())  ;
        } // PHOS cells
      } // Fill cell position histograms
      
    } // N modules
  } // Cell loop
  
  // Fill the cells after the cut on min amplitude and bad/exotic channels
  if( ncellsCut > 0 ) fhNCellsCutAmpMin->Fill(ncellsCut, GetEventWeight()) ;
  
  // Number of cells per module
  for(Int_t imod = 0; imod < fNModules; imod++ )
  {
    AliDebug(1,Form("Module %d, calo %s, N cells %d, sum Amp %f", imod, GetCalorimeterString().Data(), nCellsInModule[imod], eCellsInModule[imod]));

    fhNCellsMod     ->Fill(nCellsInModule[imod], imod, GetEventWeight()) ;
    fhSumCellsAmpMod->Fill(eCellsInModule[imod], imod, GetEventWeight()) ;
    
    fhNCellsSumAmpPerMod[imod]->Fill(eCellsInModule[imod], nCellsInModule[imod], GetEventWeight());
  }
  
  // Check energy distribution in calorimeter for selected cells
  if(fStudyWeight)
  {
    for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++)
    {
      if ( cells->GetCellNumber(iCell) < 0 ) continue; // CPV 

      AliDebug(2,Form("Cell : amp %f, absId %d", cells->GetAmplitude(iCell), cells->GetCellNumber(iCell)));
      
      Int_t nModule = GetModuleNumberCellIndexes(cells->GetCellNumber(iCell),GetCalorimeter(), icol, irow, iRCU);
      
      AliDebug(2,Form("\t module %d, column %d, row %d", nModule,icol,irow));
      
      if(nModule < fNModules)
      {
        //Check if the cell is a bad channel
        if(GetCaloUtils()->IsBadChannelsRemovalSwitchedOn())
        {
          if(GetCalorimeter()==kEMCAL)
          {
            if(GetCaloUtils()->GetEMCALChannelStatus(nModule, icol, irow)) continue;
          }
          else
          {
            if(GetCaloUtils()->GetPHOSChannelStatus (nModule, icol, irow) ) continue;
          }
        } // use bad channel map
        
        amp     = cells->GetAmplitude(iCell)*recalF;
        time    = cells->GetTime(iCell);
        id      = cells->GetCellNumber(iCell);
        
        // Amplitude recalibration if set
        GetCaloUtils()->RecalibrateCellAmplitude(amp,  GetCalorimeter(), id);
        
        // Time recalibration if set
        GetCaloUtils()->RecalibrateCellTime     (time, GetCalorimeter(), id,
                                                 GetReader()->GetInputEvent()->GetBunchCrossNumber());
        
        // Transform time to ns
        time *= 1.0e9;
        time -= fConstantTimeShift;
        
        if(time < fTimeCutMin || time > fTimeCutMax)
        {
          AliDebug(1,Form("Remove cell with Time %f",time));
          continue;
        }
        
        // Remove exotic cells, defined only for EMCAL
        if(GetCalorimeter()==kEMCAL &&
           GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(id, cells, bc)) continue;
        
        // E cross for exotic cells
        if(amp > 0.05)
        {
          if(ecellsCut > 0)
          {
            Float_t ratio    = amp/ecellsCut;
            fhECellTotalRatio    ->Fill(ecellsCut,           ratio , GetEventWeight());
            fhECellTotalLogRatio ->Fill(ecellsCut,TMath::Log(ratio), GetEventWeight());
          }
          
          if(eCellsInModule[nModule] > 0)
          {
            Float_t ratioMod = amp/eCellsInModule[nModule];
            fhECellTotalRatioMod   [nModule]->Fill(eCellsInModule[nModule],           ratioMod , GetEventWeight());
            fhECellTotalLogRatioMod[nModule]->Fill(eCellsInModule[nModule],TMath::Log(ratioMod), GetEventWeight());
          }
        } // amp > 0.5
      } // nMod > 0 < Max
    } // cell loop
  } // weight studies
  
  delete [] nCellsInModule;
  delete [] eCellsInModule;
}

//__________________________________________________________________________
/// Fill histograms releated to cluster cell position.
//__________________________________________________________________________
void AliAnaCalorimeterQA::CellInClusterPositionHistograms(AliVCluster* clus)
{
  Int_t nCaloCellsPerCluster = clus->GetNCells();
    
  UShort_t * indexList = clus->GetCellsAbsId();
  
  Float_t pos[3];
  clus->GetPosition(pos);
  
  Float_t clEnergy = clus->E();
  
  // Loop on cluster cells
  for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++)
  {
    //	printf("Index %d\n",ipos);
    Int_t absId  = indexList[ipos]; 
    
    //Get position of cell compare to cluster
    
    if ( GetCalorimeter() == kEMCAL && GetCaloUtils()->IsEMCALGeoMatrixSet() )
    {
      Double_t cellpos[] = {0, 0, 0};
      GetEMCALGeometry()->GetGlobal(absId, cellpos);
      
      fhDeltaCellClusterXNCells->Fill(pos[0]-cellpos[0], nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterYNCells->Fill(pos[1]-cellpos[1], nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterZNCells->Fill(pos[2]-cellpos[2], nCaloCellsPerCluster, GetEventWeight()) ;
      
      fhDeltaCellClusterXE->Fill(pos[0]-cellpos[0], clEnergy, GetEventWeight())  ;
      fhDeltaCellClusterYE->Fill(pos[1]-cellpos[1], clEnergy, GetEventWeight())  ;
      fhDeltaCellClusterZE->Fill(pos[2]-cellpos[2], clEnergy, GetEventWeight())  ;
      
      Float_t r     = TMath::Sqrt(pos[0]    *pos[0]     + pos[1]    * pos[1]    );
      Float_t rcell = TMath::Sqrt(cellpos[0]*cellpos[0] + cellpos[1]* cellpos[1]);
      
      fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterRE     ->Fill(r-rcell, clEnergy            , GetEventWeight())  ;
    } // EMCAL and its matrices are available
    else if ( GetCalorimeter() == kPHOS && GetCaloUtils()->IsPHOSGeoMatrixSet() )
    {
      TVector3 xyz;
      Int_t relId[4], module;
      Float_t xCell, zCell;
      
      GetPHOSGeometry()->AbsToRelNumbering(absId,relId);
      module = relId[0];
      GetPHOSGeometry()->RelPosInModule(relId,xCell,zCell);
      GetPHOSGeometry()->Local2Global(module,xCell,zCell,xyz);
      
      fhDeltaCellClusterXNCells->Fill(pos[0]-xyz.X(), nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterYNCells->Fill(pos[1]-xyz.Y(), nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterZNCells->Fill(pos[2]-xyz.Z(), nCaloCellsPerCluster, GetEventWeight()) ;
      
      fhDeltaCellClusterXE->Fill(pos[0]-xyz.X(), clEnergy, GetEventWeight())  ;
      fhDeltaCellClusterYE->Fill(pos[1]-xyz.Y(), clEnergy, GetEventWeight())  ;
      fhDeltaCellClusterZE->Fill(pos[2]-xyz.Z(), clEnergy, GetEventWeight())  ;
      
      Float_t r     = TMath::Sqrt(pos[0]  * pos[0]  + pos[1]  * pos[1] );
      Float_t rcell = TMath::Sqrt(xyz.X() * xyz.X() + xyz.Y() * xyz.Y());
      
      fhDeltaCellClusterRNCells->Fill(r-rcell, nCaloCellsPerCluster, GetEventWeight()) ;
      fhDeltaCellClusterRE     ->Fill(r-rcell, clEnergy            , GetEventWeight()) ;
    } // PHOS and its matrices are available
  } // cluster cell loop
}

//_____________________________________________________________________________________
// Study the shape of the cluster in cell units terms.
//_____________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterAsymmetryHistograms(AliVCluster* clus, Int_t absIdMax,
                                                     Bool_t goodCluster)
{
  // No use to study clusters with less than 4 cells
  if( clus->GetNCells() <= 3 ) return;
  
  Int_t dIeta = 0;
  Int_t dIphi = 0;
  
  Int_t ietaMax=-1; Int_t iphiMax = 0; Int_t rcuMax = 0;
  Int_t smMax = GetModuleNumberCellIndexes(absIdMax,GetCalorimeter(), ietaMax, iphiMax, rcuMax);
  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t absId = clus->GetCellsAbsId()[ipos];
    
    Int_t ieta=-1; Int_t iphi = 0; Int_t rcu = 0;
    Int_t sm = GetModuleNumberCellIndexes(absId,GetCalorimeter(), ieta, iphi, rcu);
    
    if(dIphi < TMath::Abs(iphi-iphiMax)) dIphi = TMath::Abs(iphi-iphiMax);
    
    if(smMax==sm)
    {
      if(dIeta < TMath::Abs(ieta-ietaMax)) dIeta = TMath::Abs(ieta-ietaMax);
    }
    else
    {
      Int_t ietaShift    = ieta;
      Int_t ietaMaxShift = ietaMax;
      if (ieta > ietaMax)  ietaMaxShift+=48;
      else                 ietaShift   +=48;
      if(dIeta < TMath::Abs(ietaShift-ietaMaxShift)) dIeta = TMath::Abs(ietaShift-ietaMaxShift);
    }
  }// Fill cell-cluster histogram loop
  
  Float_t dIA = 1.*(dIphi-dIeta)/(dIeta+dIphi);

  if(goodCluster)
  {
    // Was cluster matched?
    Bool_t matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(),GetReader()->GetInputEvent());
    
    if     (clus->E() < 2 ) fhDeltaIEtaDeltaIPhiE0[matched]->Fill(dIeta, dIphi, GetEventWeight());
    else if(clus->E() < 6 ) fhDeltaIEtaDeltaIPhiE2[matched]->Fill(dIeta, dIphi, GetEventWeight());
    else                    fhDeltaIEtaDeltaIPhiE6[matched]->Fill(dIeta, dIphi, GetEventWeight());
    
    fhDeltaIA[matched]->Fill(clus->E(), dIA, GetEventWeight());
    
    if(clus->E() > 0.5)
    {
      fhDeltaIAL0    [matched]->Fill(clus->GetM02()   , dIA, GetEventWeight());
      fhDeltaIAL1    [matched]->Fill(clus->GetM20()   , dIA, GetEventWeight());
      fhDeltaIANCells[matched]->Fill(clus->GetNCells(), dIA, GetEventWeight());
    }
    
    // Origin of  clusters
    Int_t  nLabel = clus->GetNLabels();
    Int_t* labels = clus->GetLabels();
    
    if(IsDataMC())
    {
      Int_t tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabel, GetReader(),GetCalorimeter());
      if(   GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) && 
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0)    &&
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)    &&
           !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
        fhDeltaIAMC[0]->Fill(clus->E(), dIA, GetEventWeight()); // Pure Photon
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron) &&
               !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)  ){
        fhDeltaIAMC[1]->Fill(clus->E(), dIA, GetEventWeight()); // Pure electron
      }
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) && 
                GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)  ){
        fhDeltaIAMC[2]->Fill(clus->E(), dIA, GetEventWeight()); // Converted cluster
      }
      else if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)){ 
        fhDeltaIAMC[3]->Fill(clus->E(), dIA, GetEventWeight()); // Hadrons
      }
      
    }  // MC
  } // good cluster
  else
  {
    if     (clus->E() < 2 ) fhBadClusterDeltaIEtaDeltaIPhiE0->Fill(dIeta, dIphi, GetEventWeight());
    else if(clus->E() < 6 ) fhBadClusterDeltaIEtaDeltaIPhiE2->Fill(dIeta, dIphi, GetEventWeight());
    else                    fhBadClusterDeltaIEtaDeltaIPhiE6->Fill(dIeta, dIphi, GetEventWeight());
    
    fhBadClusterDeltaIA->Fill(clus->E(), dIA, GetEventWeight());
  }
}

//__________________________________________________________________________________________________________________
/// Fill calorimeter cluster related histograms. Only good clusters, excluded bad previously.
/// \param clus: calorimeter cluster pointer.
/// \param caloClusters: full list of clusters, used for comparison of bad and good clusters timing.
/// \param cells: list of cells, needed to get the cell with highest energy in cluster.
/// \param absIdMax: id number of cell with highest energy in the cluster.
/// \param maxCellFraction: ratio E_cell_max/ E_cluster.
/// \param eCrossFrac: exoticity fraction.
/// \param tmax: time of highest energy cell in cluster.
//__________________________________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterHistograms(AliVCluster* clus, const TObjArray *caloClusters, AliVCaloCells * cells,
                                            Int_t absIdMax, Double_t maxCellFraction, Float_t eCrossFrac,
                                            Double_t tmax)
{
  Double_t tof = clus->GetTOF()*1.e9;
  if(tof>400) tof-=fConstantTimeShift;
  
  fhLambda0             ->Fill(clus->E(), clus->GetM02()       , GetEventWeight());
  fhLambda1             ->Fill(clus->E(), clus->GetM20()       , GetEventWeight());
  //  fhDispersion          ->Fill(clus->E(), clus->GetDispersion(), GetEventWeight());
  
  fhClusterMaxCellDiff  ->Fill(clus->E(), maxCellFraction, GetEventWeight());
  fhClusterMaxCellECross->Fill(clus->E(), eCrossFrac     , GetEventWeight());
  fhClusterTimeEnergy   ->Fill(clus->E(), tof            , GetEventWeight());
  
  if(fStudyClustersAsymmetry) ClusterAsymmetryHistograms(clus,absIdMax,kTRUE);
  
  Int_t    nModule = GetModuleNumber(clus);
  Int_t    nCaloCellsPerCluster = clus->GetNCells();
  
  // Clusters in event time difference
  for(Int_t iclus2 = 0; iclus2 < caloClusters->GetEntriesFast(); iclus2++ )
  {
    AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(iclus2);
    
    if( clus->GetID() == clus2->GetID() ) continue;
    
    if( clus->GetM02() > 0.01 && clus2->GetM02() > 0.01 )
    {
      Int_t    nModule2 = GetModuleNumber(clus2);

      Double_t tof2   = clus2->GetTOF()*1.e9;    
      if(tof2>400) tof2-=fConstantTimeShift;

      fhClusterPairDiffTimeE  ->Fill(clus->E(), tof-tof2, GetEventWeight());
      
      if ( nModule2 == nModule )
        fhClusterPairDiffTimeESameMod->Fill(clus->E(), tof-tof2, GetEventWeight());
    }
  }        
  
  if(nCaloCellsPerCluster > 1)
  {
    // Check time of cells respect to max energy cell
    
//    if(fFillAllCellTimeHisto) 
//    {
//      // Get some time averages
//      Double_t timeAverages[2] = {0.,0.};
//      CalculateAverageTime(clus, cells, timeAverages);
//
//      fhClusterMaxCellDiffAverageTime      ->Fill(clus->E(), tmax-timeAverages[0], GetEventWeight());
//      fhClusterMaxCellDiffWeightedTime     ->Fill(clus->E(), tmax-timeAverages[1], GetEventWeight());
//    }
    
    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t absId  = clus->GetCellsAbsId()[ipos];             
      if( absId == absIdMax || cells->GetCellAmplitude(absIdMax) < 0.01 ) continue;
      
      Float_t frac = cells->GetCellAmplitude(absId)/cells->GetCellAmplitude(absIdMax);            
      fhClusterMaxCellCloseCellRatio->Fill(clus->E(), frac, GetEventWeight());
      fhClusterMaxCellCloseCellDiff ->Fill(clus->E(), cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId), GetEventWeight());
      
      if(fFillAllCellTimeHisto) 
      {
        Double_t time  = cells->GetCellTime(absId);
        GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
        
        Float_t diff = (tmax-(time*1.0e9-fConstantTimeShift));
        fhCellTimeSpreadRespectToCellMax->Fill(clus->E(), diff, GetEventWeight());
        if(TMath::Abs(TMath::Abs(diff) > 100) && clus->E() > 1 ) fhCellIdCellLargeTimeSpread->Fill(absId, GetEventWeight());
      }
      
    } // Fill cell-cluster histogram loop
    
  } // Check time and energy of cells respect to max energy cell if cluster of more than 1 cell
    
  Float_t e   = fClusterMomentum.E();
  Float_t pt  = fClusterMomentum.Pt();
  Float_t eta = fClusterMomentum.Eta();
  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  AliDebug(1,Form("cluster: E %2.3f, pT %2.3f, eta %2.3f, phi %2.3f",e,pt,eta,phi*TMath::RadToDeg()));
  
  fhE     ->Fill(e, GetEventWeight());
  if(nModule >=0 && nModule < fNModules)
  {
    fhEMod     ->Fill(e, nModule, GetEventWeight());
    fhEWeirdMod->Fill(e, nModule, GetEventWeight()); // different binning
  }
  
  fhPt     ->Fill(pt , GetEventWeight());
  fhPhi    ->Fill(phi, GetEventWeight());
  fhEta    ->Fill(eta, GetEventWeight());
  
  if ( fFillAllTH3 )  fhEtaPhiE->Fill(eta, phi, e, GetEventWeight());
  else if ( e > 0.5 ) fhEtaPhi ->Fill(eta, phi,    GetEventWeight());

  // Cells per cluster
  fhNCellsPerCluster->Fill(e, nCaloCellsPerCluster, GetEventWeight());

  if(e > 100) fhNCellsPerClusterWeirdMod->Fill(nCaloCellsPerCluster, nModule, GetEventWeight());
  
  // Position
  if(fFillAllPosHisto2)
  {
    Float_t pos[3] ;     
    clus->GetPosition(pos);
    
    fhXE     ->Fill(pos[0], e, GetEventWeight());
    fhYE     ->Fill(pos[1], e, GetEventWeight());
    fhZE     ->Fill(pos[2], e, GetEventWeight());
    
    fhXYZ    ->Fill(pos[0], pos[1], pos[2], GetEventWeight());
    
    fhXNCells->Fill(pos[0], nCaloCellsPerCluster, GetEventWeight());
    fhYNCells->Fill(pos[1], nCaloCellsPerCluster, GetEventWeight());
    fhZNCells->Fill(pos[2], nCaloCellsPerCluster, GetEventWeight());
      
    Float_t rxyz = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]);//+pos[2]*pos[2]);
      
    fhRE     ->Fill(rxyz, e                   , GetEventWeight());
    fhRNCells->Fill(rxyz, nCaloCellsPerCluster, GetEventWeight());
  }
  
  if( nModule >= 0 && nModule < fNModules ) fhNCellsPerClusterMod[nModule]->Fill(e, nCaloCellsPerCluster, GetEventWeight());
}

//____________________________________________________________________________
/// Fill clusters related histograms, execute here the loop of clusters
/// apply basic selection cuts (track matching, goondess, exoticity, timing)
/// and the call to the different methods
/// filling different type of histograms:
/// * Basic cluster histograms for good or bad clusters
/// * Excotic cluster histograms
/// * Cells in cluster
/// * Invariant mass
/// * Matched clusters histograms
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCalorimeterQA::ClusterLoopHistograms(const TObjArray *caloClusters,
                                                AliVCaloCells* cells)
{
  Int_t  nLabel                = 0  ;
  Int_t *labels                = 0x0;
  Int_t  nCaloClusters         = caloClusters->GetEntriesFast() ;
  Int_t  nCaloClustersAccepted = 0  ;
  Int_t  nCaloCellsPerCluster  = 0  ;
  Bool_t matched               = kFALSE;
  Int_t  nModule               =-1  ;
  
  // Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
//GetReader()->GetVertex(v);

  Int_t   *nClustersInModule  = new Int_t  [fNModules];
  Float_t *energyInModule     = new Float_t[fNModules];
  for(Int_t imod = 0; imod < fNModules; imod++ )
  {
    nClustersInModule[imod] = 0;
    energyInModule   [imod] = 0;
  }
  
  AliDebug(1,Form("In %s there are %d clusters", GetCalorimeterString().Data(), nCaloClusters));
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++)
  {
    AliDebug(1,Form("Cluster: %d/%d, data %d",iclus+1,nCaloClusters,GetReader()->GetDataType()));
    
    AliVCluster* clus =  (AliVCluster*) caloClusters->At(iclus);
    
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    // Cut on time of clusters
    Double_t tof = clus->GetTOF()*1.e9;
    if(tof>400) tof-=fConstantTimeShift;

    if( tof < fTimeCutMin || tof > fTimeCutMax )
    { 
      AliDebug(1,Form("Remove cluster with TOF %f",tof));
      continue;
    }    
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if(IsFiducialCutOn()) in =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum.Eta(),fClusterMomentum.Phi(),GetCalorimeter()) ;
    if(!in) continue;
    
    // MC labels
    nLabel = clus->GetNLabels();
    labels = clus->GetLabels();
    
    // SuperModule number of cluster
    nModule = GetModuleNumber(clus);
    
    // Cells per cluster
    nCaloCellsPerCluster = clus->GetNCells();
    
    // Cluster mathed with track?
    matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(), GetReader()->GetInputEvent());
    
    // Get time of max cell
    Double_t tmax  = cells->GetCellTime(absIdMax);
    GetCaloUtils()->RecalibrateCellTime(tmax, GetCalorimeter(), absIdMax,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    tmax*=1.e9;
    tmax-=fConstantTimeShift;
    //
    // Fill histograms related to single cluster 
    //

    // Fill some histograms before applying the exotic cell / bad map cut
    if(fStudyBadClusters) 
    { 
      fhNCellsPerClusterNoCut->Fill(clus->E(), nCaloCellsPerCluster, GetEventWeight());
      
      if(nModule >=0 && nModule < fNModules)
      {
        fhNCellsPerClusterModNoCut[nModule]->Fill(clus->E(), nCaloCellsPerCluster, GetEventWeight());
        if(clus->E() > 100) fhNCellsPerClusterWeirdModNoCut->Fill(nCaloCellsPerCluster, nModule, GetEventWeight());
      }
      
      fhClusterMaxCellDiffNoCut->Fill(clus->E(), maxCellFraction, GetEventWeight());
    }
    
    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    GetCaloUtils()->RecalibrateCellAmplitude(ampMax,GetCalorimeter(), absIdMax);
    
    if(fStudyExotic) ExoticHistograms(absIdMax, ampMax, clus, cells);
    
    // Check bad clusters if requested and rejection was not on
    Bool_t goodCluster = IsGoodCluster(absIdMax, clus->GetM02(), nCaloCellsPerCluster, cells);
              
    Float_t eCrossFrac = 0;
    if(ampMax > 0.01) eCrossFrac = 1-GetECross(absIdMax,cells)/ampMax;
    
    AliDebug(1,Form("Accept cluster? %d",goodCluster));
          
    if(!goodCluster) 
    {
      if ( fStudyBadClusters ) BadClusterHistograms(clus, caloClusters, 
                                                    cells, absIdMax, 
                                                    maxCellFraction, 
                                                    eCrossFrac, tmax);
      continue;
    }
        

    ClusterHistograms(clus, caloClusters, cells, absIdMax, 
                      maxCellFraction, eCrossFrac, tmax);
    
    nCaloClustersAccepted++;
    nModule = GetModuleNumber(clus);
    if(nModule >=0 && nModule < fNModules && fClusterMomentum.E() > 2*fCellAmpMin)
    {
     nClustersInModule[nModule]++;
     if(clus->E() > 0.5)
      energyInModule   [nModule] += clus->E();
    }
    
    // Cluster weights
    if(fStudyWeight) WeightHistograms(clus, cells);
    
    // Cells in cluster position
    if(fFillAllPosHisto) CellInClusterPositionHistograms(clus);
    
    // Fill histograms related to single cluster, mc vs data
    Int_t  mcOK = kFALSE;
    Int_t  pdg  = -1;
    if(IsDataMC() && nLabel > 0 && labels) 
      mcOK = ClusterMCHistograms(matched, labels, nLabel, pdg);

    // Matched clusters with tracks, also do some MC comparison, needs input from ClusterMCHistograms
    if( matched &&  fFillAllTMHisto)
      ClusterMatchedWithTrackHistograms(clus,mcOK,pdg);
    
    // Invariant mass
    // Try to reduce background with a mild shower shape cut and no more than 1 maxima 
    // in cluster and remove low energy clusters
        
    if (   fFillAllPi0Histo 
        && nCaloClusters > 1 
        && nCaloCellsPerCluster > 1 
        && GetCaloUtils()->GetNumberOfLocalMaxima(clus,cells) == 1 
        && clus->GetM02() < fInvMassMaxM02Cut 
        && clus->GetM02() > fInvMassMinM02Cut 
        && clus->E() > fInvMassMinECut 
        && clus->E() < fInvMassMaxECut 
        )
      InvariantMassHistograms(iclus, nModule, caloClusters,cells);
    
  } // Cluster loop
  
  // Number of clusters histograms
  if(nCaloClustersAccepted > 0) fhNClusters->Fill(nCaloClustersAccepted, GetEventWeight());
  
  // Number of clusters per module
  for(Int_t imod = 0; imod < fNModules; imod++ )
  { 
    AliDebug(1,Form("Module %d calo %s clusters %d, sum E %f", imod, GetCalorimeterString().Data(), nClustersInModule[imod], energyInModule[imod]));

    fhNClustersMod        ->Fill(nClustersInModule[imod], imod, GetEventWeight());
    fhSumClustersEnergyMod->Fill(energyInModule   [imod], imod, GetEventWeight());
    
    fhNClustersSumEnergyPerMod[imod]->Fill(energyInModule[imod], nClustersInModule[imod], GetEventWeight());
  }
  
  delete [] nClustersInModule;
  delete [] energyInModule;
}

//__________________________________________________________________________________
/// Fill histograms depending on the MC origin information.
/// Only possible for simulations.
//__________________________________________________________________________________
Bool_t AliAnaCalorimeterQA::ClusterMCHistograms(Bool_t matched,const Int_t * labels,
                                                Int_t nLabels, Int_t & pdg )
{
  if(!labels || nLabels<=0)
  {
    AliWarning(Form("Strange, labels array %p, n labels %d", labels,nLabels));
    return kFALSE;
  }
  
  AliDebug(1,Form("Primaries: nlabels %d",nLabels));
  
  Float_t e   = fClusterMomentum.E();
  Float_t eta = fClusterMomentum.Eta();
  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  AliAODMCParticle * aodprimary  = 0x0;
  TParticle * primary = 0x0;
  
  // Play with the MC stack if available
  Int_t label = labels[0];
  
  if(label < 0) 
  {
    AliDebug(1,Form(" *** bad label ***:  label %d", label));
    return kFALSE;
  }
  
  Int_t pdg0    =-1; Int_t status = -1; Int_t iMother = -1; Int_t iParent = -1;
  Float_t vxMC  = 0; Float_t vyMC  = 0;
  Float_t eMC   = 0; //Float_t ptMC= 0;
  Float_t phiMC = 0; Float_t etaMC = 0;
  Int_t charge  = 0;
  
  // Check the origin.
  Int_t tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabels, GetReader(),GetCalorimeter());
  
  if     ( GetReader()->ReadStack() && 
          !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown))
  { // If MC stack and known tag
    
    if( label >= GetMCStack()->GetNtrack()) 
    {
      AliDebug(1,Form("*** large label ***:  label %d, n tracks %d", label, GetMCStack()->GetNtrack()));
      return kFALSE;
    }
    
    primary = GetMCStack()->Particle(label);
    iMother = label;
    pdg0    = TMath::Abs(primary->GetPdgCode());
    pdg     = pdg0;
    status  = primary->GetStatusCode();
    vxMC    = primary->Vx();
    vyMC    = primary->Vy();
    iParent = primary->GetFirstMother();
    
    AliDebug(1,"Cluster most contributing mother:");
    AliDebug(1,Form("\t Mother label %d, pdg %d, %s, status %d, parent %d",iMother, pdg0, primary->GetName(),status, iParent));
    
    
    // Get final particle, no conversion products
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion))
    {
      // Get the parent
      primary = GetMCStack()->Particle(iParent);
      pdg = TMath::Abs(primary->GetPdgCode());
      
      AliDebug(2,"Converted cluster!. Find before conversion:");

      while((pdg == 22 || pdg == 11) && status != 1)
      {
        Int_t iMotherOrg = iMother;
        iMother = iParent;
        primary = GetMCStack()->Particle(iMother);
        status  = primary->GetStatusCode();
        pdg     = TMath::Abs(primary->GetPdgCode());
        iParent = primary->GetFirstMother();
        
        // If gone too back and non stable, assign the decay photon/electron
        // there are other possible decays, ignore them for the moment
        if(pdg==111 || pdg==221)
        {
          primary = GetMCStack()->Particle(iMotherOrg);
          break;
        }
        
        if( iParent < 0 )
        {
          iParent = iMother;
          break;
        }
        
        AliDebug(2,Form("\t pdg %d, index %d, %s, status %d",pdg, iMother,  primary->GetName(),status));
      }	

      AliDebug(1,"Converted Cluster mother before conversion:");
      AliDebug(1,Form("\t Mother label %d, pdg %d, %s, status %d, parent %d",iMother, pdg, primary->GetName(), status, iParent));
      
    }
    
    // Overlapped pi0 (or eta, there will be very few), get the meson
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
       GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta))
    {
      AliDebug(2,"Overlapped Meson decay!, Find it:");

      while(pdg != 111 && pdg != 221)
      {     
        //printf("iMother %d, pdg %d, iParent %d, pdg %d\n",iMother,pdg,iParent,GetMCStack()->Particle(iParent)->GetPdgCode());
        iMother = iParent;
        primary = GetMCStack()->Particle(iMother);
        status  = primary->GetStatusCode();
        pdg     = TMath::Abs(primary->GetPdgCode());
        iParent = primary->GetFirstMother();

        if( iParent < 0 ) break;
        
        AliDebug(2,Form("\t pdg %d, %s, index %d",pdg,  primary->GetName(),iMother));
        
        if(iMother==-1) 
        {
          AliWarning("Tagged as Overlapped photon but meson not found, why?");
          //break;
        }
      }
        
      AliDebug(2,Form("Overlapped %s decay, label %d",primary->GetName(),iMother));
    }
    
    eMC    = primary->Energy();
    //ptMC   = primary->Pt();
    phiMC  = primary->Phi();
    etaMC  = primary->Eta();
    pdg    = TMath::Abs(primary->GetPdgCode());
    charge = (Int_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
  }
  else if( GetReader()->ReadAODMCParticles() && 
          !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown))
  {
    // If MC AOD and known tag
    // Get the list of MC particles
    if(!GetReader()->GetAODMCParticles())
      AliFatal("MCParticles not available!");
    
    aodprimary = (AliAODMCParticle*) (GetReader()->GetAODMCParticles())->At(label);
    iMother = label;
    pdg0    = TMath::Abs(aodprimary->GetPdgCode());
    pdg     = pdg0;
    status  = aodprimary->IsPrimary();
    vxMC    = aodprimary->Xv();
    vyMC    = aodprimary->Yv();
    iParent = aodprimary->GetMother();
    
    AliDebug(1,"Cluster most contributing mother:");
    AliDebug(1,Form("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d, parent %d",
                    iMother, pdg0, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary(), iParent));
    
    //Get final particle, no conversion products
    if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion) )
    {
      AliDebug(2,"Converted cluster!. Find before conversion:");
      
      // Get the parent
      aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles())->At(iParent);
      pdg = TMath::Abs(aodprimary->GetPdgCode());
        
      while ((pdg == 22 || pdg == 11) && !aodprimary->IsPhysicalPrimary()) 
      {
        Int_t iMotherOrg = iMother;
        iMother    = iParent;
        aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles())->At(iMother);
        status     = aodprimary->IsPrimary();
        iParent    = aodprimary->GetMother();
        pdg        = TMath::Abs(aodprimary->GetPdgCode());

        // If gone too back and non stable, assign the decay photon/electron
        // there are other possible decays, ignore them for the moment
        if( pdg == 111 || pdg == 221 )
        {
          aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles())->At(iMotherOrg);
          break;
        }        
        
        if( iParent < 0 )
        {
          iParent = iMother;
          break;
        }
        
        AliDebug(2,Form("\t pdg %d, index %d, Primary? %d, Physical Primary? %d",
                        pdg, iMother, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary()));
      }	
      
      AliDebug(1,"Converted Cluster mother before conversion:");
      AliDebug(1,Form("\t Mother label %d, pdg %d, parent %d, Primary? %d, Physical Primary? %d",
                      iMother, pdg, iParent, aodprimary->IsPrimary(), aodprimary->IsPhysicalPrimary()));
    
    }
    
    //Overlapped pi0 (or eta, there will be very few), get the meson
    if(GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) || 
       GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta))
    {
      AliDebug(2,Form("Overlapped Meson decay!, Find it: PDG %d, mom %d",pdg, iMother));
  
      while(pdg != 111 && pdg != 221)
      {
        iMother    = iParent;
        aodprimary = (AliAODMCParticle*)(GetReader()->GetAODMCParticles())->At(iMother);
        status     = aodprimary->IsPrimary();
        iParent    = aodprimary->GetMother();
        pdg        = TMath::Abs(aodprimary->GetPdgCode());

        if( iParent < 0 ) break;
        
        AliDebug(2,Form("\t pdg %d, index %d",pdg, iMother));
        
        if(iMother==-1) 
        {
          AliWarning("Tagged as Overlapped photon but meson not found, why?");
          //break;
        }
      }	
      
      AliDebug(2,Form("Overlapped %s decay, label %d",aodprimary->GetName(),iMother));
    }	
    
    status = aodprimary->IsPrimary();
    eMC    = aodprimary->E();
    //ptMC   = aodprimary->Pt();
    phiMC  = aodprimary->Phi();
    etaMC  = aodprimary->Eta();
    pdg    = TMath::Abs(aodprimary->GetPdgCode());
    charge = aodprimary->Charge();
  }
  
  //Float_t vz = primary->Vz();
  Float_t rVMC = TMath::Sqrt(vxMC*vxMC + vyMC*vyMC);
  if( ( pdg == 22 || TMath::Abs(pdg) == 11 ) && status != 1 )
  {
    fhEMVxyz   ->Fill(vxMC, vyMC, GetEventWeight());//,vz);
    fhEMR      ->Fill(e   , rVMC, GetEventWeight());
  }
  
  //printf("reco e %f, pt %f, phi %f, eta %f \n", e, pt, phi, eta);
  //printf("prim e %f, pt %f, phi %f, eta %f \n", eMC,ptMC,phiMC ,etaMC );
  //printf("vertex: vx %f, vy %f, vz %f, r %f \n", vxMC, vyMC, vz, r);
  
  // Overlapped pi0 (or eta, there will be very few)
  Int_t mcIndex = -1;
  if     ( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0     ) )
  {
    mcIndex = kmcPi0;
  } // Overlapped pizero decay
  else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta     ) )
  {
    mcIndex = kmcEta;
  } // Overlapped eta decay
  else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton  ) )
  {
    if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion))
       mcIndex = kmcPhotonConv ;
    else
       mcIndex = kmcPhoton ;
  } // photon
  else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron) )
  {
    mcIndex = kmcElectron;
    fhEMVxyz->Fill(vxMC, vyMC, GetEventWeight());//,vz);
    fhEMR   ->Fill(e   , rVMC, GetEventWeight());
  }
  else if(charge == 0)
  {
    mcIndex = kmcNeHadron;
    fhHaVxyz->Fill(vxMC, vyMC, GetEventWeight());//,vz);
    fhHaR   ->Fill(e   , rVMC, GetEventWeight());
  }
  else if(charge!=0)
  {
    mcIndex = kmcChHadron;
    fhHaVxyz->Fill(vxMC, vyMC, GetEventWeight());//,vz);
    fhHaR   ->Fill(e   , rVMC, GetEventWeight());
  }

  //printf("mc index %d\n",mcIndex);
  
  if( mcIndex >= 0  && mcIndex < 7  && e > 0.5 && eMC > 0.5)
  {
    fhRecoMCE       [mcIndex][(matched)]->Fill(e  , eMC      , GetEventWeight());
    fhRecoMCEta     [mcIndex][(matched)]->Fill(eta, etaMC    , GetEventWeight());
    fhRecoMCPhi     [mcIndex][(matched)]->Fill(phi, phiMC    , GetEventWeight());
    fhRecoMCRatioE  [mcIndex][(matched)]->Fill(e  , e/eMC    , GetEventWeight());
    fhRecoMCDeltaE  [mcIndex][(matched)]->Fill(e  , eMC-e    , GetEventWeight());
    fhRecoMCDeltaPhi[mcIndex][(matched)]->Fill(e  , phiMC-phi, GetEventWeight());
    fhRecoMCDeltaEta[mcIndex][(matched)]->Fill(e  , etaMC-eta, GetEventWeight());
  }
  
  if( primary || aodprimary ) return kTRUE ;
  else                        return kFALSE;
}

//_________________________________________________________________________________________________________
/// Histograms for clusters matched with tracks.
/// \param clus: cluster pointer.
/// \param okPrimary: primary was found before.
/// \param pdg: primary pdg number found earlier.
//_________________________________________________________________________________________________________
void AliAnaCalorimeterQA::ClusterMatchedWithTrackHistograms(AliVCluster *clus, Bool_t okPrimary, Int_t pdg)
{
  Float_t e   = fClusterMomentum.E();
  Float_t pt  = fClusterMomentum.Pt();
  Float_t eta = fClusterMomentum.Eta();
  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi +=TMath::TwoPi();

  fhECharged   ->Fill(e  , GetEventWeight());
  fhPtCharged  ->Fill(pt , GetEventWeight());
  fhPhiCharged ->Fill(phi, GetEventWeight());
  fhEtaCharged ->Fill(eta, GetEventWeight());
    
  if ( fFillAllTH3 )  fhEtaPhiECharged->Fill(eta, phi, e, GetEventWeight());
  else if ( e > 0.5 ) fhEtaPhiCharged ->Fill(eta, phi,    GetEventWeight());
  
  // Study the track and matched cluster if track exists.
    
  AliVTrack *track = GetCaloUtils()->GetMatchedTrack(clus, GetReader()->GetInputEvent());
  
  if(!track) return ;
    
  Double_t tpt   = track->Pt();
  Double_t tmom  = track->P();
  Double_t dedx  = track->GetTPCsignal();
  Int_t    nITS  = track->GetNcls(0);
  Int_t    nTPC  = track->GetNcls(1);
  Bool_t positive = kFALSE;
  if(track) positive = (track->Charge()>0);
  
  // Residuals
  Float_t deta  = clus->GetTrackDz();
  Float_t dphi  = clus->GetTrackDx();
  Double_t  dR  = TMath::Sqrt(dphi*dphi + deta*deta);
  Int_t nModule = GetModuleNumber(clus);

  if( TMath::Abs(dphi) < 999 )
  {
    if(positive)
    {
      fhTrackMatchedDEtaPos->Fill(e, deta, GetEventWeight());
      fhTrackMatchedDPhiPos->Fill(e, dphi, GetEventWeight());
        
      if(e > 0.5) 
      {
        fhTrackMatchedDEtaPosMod ->Fill(deta, nModule, GetEventWeight());
        fhTrackMatchedDPhiPosMod ->Fill(dphi, nModule, GetEventWeight());
        fhTrackMatchedDEtaDPhiPos->Fill(deta, dphi   , GetEventWeight());
      }
    }
    else 
    {
      fhTrackMatchedDEtaNeg->Fill(e, deta, GetEventWeight());
      fhTrackMatchedDPhiNeg->Fill(e, dphi, GetEventWeight());
        
      if(e > 0.5) 
      {
        fhTrackMatchedDEtaNegMod ->Fill(deta, nModule, GetEventWeight());
        fhTrackMatchedDPhiNegMod ->Fill(dphi, nModule, GetEventWeight());
        fhTrackMatchedDEtaDPhiNeg->Fill(deta, dphi   , GetEventWeight());
      }
    }
  }
  
  Double_t eOverP = e/tmom;
  fh1EOverP->Fill(tpt, eOverP, GetEventWeight());
  if(e > 0.5 && tpt > 0.5)  fh1EOverPMod->Fill(eOverP, nModule, GetEventWeight());
    
  if(dR < 0.02)
  {
    fh1EOverPR02->Fill(tpt, eOverP, GetEventWeight());
    if(e > 0.5 && tpt > 0.5) fh1EOverPR02Mod->Fill(eOverP, nModule, GetEventWeight());
    
    if(dedx > 60 && dedx < 100) 
    {
      fh1EleEOverP->Fill(tpt, eOverP, GetEventWeight());
      if(e > 0.5 && tpt > 0.5)  fh1EleEOverPMod->Fill(eOverP, nModule, GetEventWeight());
    }
  }
  
  fh2dR->Fill(e, dR, GetEventWeight());
  fh2MatchdEdx->Fill(tmom, dedx, GetEventWeight());
  
  if(e > 0.5 && tmom > 0.5) 
  {
    fh2dRMod->Fill(dR, nModule, GetEventWeight());
    fh2MatchdEdxMod->Fill(dedx, nModule, GetEventWeight());
  }

  if(dR < 0.02 && eOverP > 0.6 && eOverP < 1.2
     && clus->GetNCells() > 1 && nITS > 3 && nTPC > 20) 
  {
    fh2EledEdx->Fill(tmom, dedx, GetEventWeight());
    if(e > 0.5 && tmom > 0.5) fh2EledEdxMod->Fill(dedx, nModule, GetEventWeight());
  }
  
  if(IsDataMC() && okPrimary)
  { 
    Double_t  charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    
    if(TMath::Abs(pdg) == 11)
    {
      fhMCEle1EOverP   ->Fill(tpt , eOverP, GetEventWeight());
      fhMCEle1dR       ->Fill(dR  ,         GetEventWeight());
      fhMCEle2MatchdEdx->Fill(tmom,  dedx , GetEventWeight());
        
      if(dR < 0.02)
      {
        fhMCEle1EOverPR02->Fill(tpt, eOverP, GetEventWeight());
        if(dedx > 60 && dedx < 100) fhMCEle1EleEOverP->Fill(tpt, eOverP, GetEventWeight());
      }
    }
    else if(charge!=0)
    {
      fhMCChHad1EOverP   ->Fill(tpt , eOverP, GetEventWeight());
      fhMCChHad1dR       ->Fill(dR  ,         GetEventWeight());
      fhMCChHad2MatchdEdx->Fill(tmom, dedx  , GetEventWeight());
        
      if(dR < 0.02)
      {
        fhMCChHad1EOverPR02->Fill(tpt, eOverP, GetEventWeight());
        if(dedx > 60 && dedx < 100) fhMCChHad1EleEOverP->Fill(tpt, eOverP, GetEventWeight());
      }
    }
    else if(charge == 0)
    {
      fhMCNeutral1EOverP   ->Fill(tpt , eOverP, GetEventWeight());
      fhMCNeutral1dR       ->Fill(dR  ,         GetEventWeight());
      fhMCNeutral2MatchdEdx->Fill(tmom, dedx  , GetEventWeight());
        
      if(dR < 0.02)
      {
        fhMCNeutral1EOverPR02->Fill(tpt, eOverP, GetEventWeight());
        if(dedx > 60 && dedx < 100) fhMCNeutral1EleEOverP->Fill(tpt, eOverP, GetEventWeight());
      }
    }
  } // DataMC
}

//___________________________________
/// Correlate information from PHOS and EMCal
/// and DCal and with V0 and track multiplicity
//___________________________________
void AliAnaCalorimeterQA::Correlate()
{
  // Clusters arrays
  TObjArray * caloClustersEMCAL = GetEMCALClusters();
  TObjArray * caloClustersPHOS  = GetPHOSClusters();
  
  if(!caloClustersEMCAL || !caloClustersPHOS)
  {
    AliDebug(1,Form("PHOS (%p) or EMCAL (%p) clusters array not available, do not correlate",caloClustersPHOS,caloClustersEMCAL));
    return ;
  }
  
  // Cells arrays
  AliVCaloCells * cellsEMCAL = GetEMCALCells();
  AliVCaloCells * cellsPHOS  = GetPHOSCells();
  
  if(!cellsEMCAL || !cellsPHOS)
  {
    AliDebug(1,Form("PHOS (%p) or EMCAL (%p) cells array ot available, do not correlate",cellsPHOS,cellsEMCAL));
    return ;
  }

  // Clusters parameters
  Int_t nclEMCAL = 0;
  Int_t nclDCAL  = 0;
  Int_t nclPHOS  = 0;

  Float_t sumClusterEnergyEMCAL = 0;
  Float_t sumClusterEnergyDCAL  = 0;
  Float_t sumClusterEnergyPHOS  = 0;
  
  Int_t iclus = 0;
  Float_t energy = 0;
  AliVCluster* cluster = 0;
  for(iclus = 0 ; iclus <  caloClustersEMCAL->GetEntriesFast() ; iclus++) 
  {
    cluster = (AliVCluster*)caloClustersEMCAL->At(iclus);
    Float_t energy = cluster->E();
    
    if( energy < 0.5 ) continue;
    
    if(cluster->GetCellsAbsId()[0] < 12288)
    {
      nclEMCAL++;
      sumClusterEnergyEMCAL += energy;
    }
    else
    {
      nclDCAL++;
      sumClusterEnergyDCAL  += energy;
    }
  }
  
  for(iclus = 0 ; iclus <  caloClustersPHOS ->GetEntriesFast(); iclus++) 
  {
    cluster = (AliVCluster*) caloClustersPHOS->At(iclus);

    energy = cluster->E();
    
    if( energy < 0.5 ) continue;
    
    nclPHOS++;
    sumClusterEnergyPHOS += energy;
  }
  
  // Cells parameters
  Int_t ncellsEMCAL = 0 ;
  Int_t ncellsDCAL  = 0 ;
  Int_t ncellsPHOS  = 0;
  
  Float_t sumCellEnergyEMCAL = 0;
  Float_t sumCellEnergyDCAL  = 0;
  Float_t sumCellEnergyPHOS  = 0;
  Int_t icell = 0;
  for(icell = 0 ; icell < cellsEMCAL->GetNumberOfCells()  ; icell++) 
  {
    Float_t  amp = cellsEMCAL->GetAmplitude(icell);
    Int_t cellId = cellsEMCAL->GetCellNumber(icell);
    
    if (amp < fEMCALCellAmpMin) continue;
    
    if( cellId < 12288) 
    {
      ncellsEMCAL++;
      sumCellEnergyEMCAL += amp;
    }
    else
    {
      ncellsDCAL++;
      sumCellEnergyDCAL  += amp;
    }
  }
  
  for(icell = 0 ; icell <  cellsPHOS->GetNumberOfCells(); icell++) 
  {
    Float_t  amp = cellsPHOS->GetAmplitude(icell);
    Int_t cellId = cellsPHOS->GetCellNumber(icell);

    if ( cellId < 0 ) continue ; // CPV
    
    if ( amp < fPHOSCellAmpMin ) continue;
    
    ncellsPHOS++;
    sumCellEnergyPHOS += amp;
  }
  
  // Fill Histograms
  
  fhEMCALPHOSCorrNClusters->Fill(nclEMCAL             , nclPHOS             , GetEventWeight());
  fhEMCALPHOSCorrEClusters->Fill(sumClusterEnergyEMCAL, sumClusterEnergyPHOS, GetEventWeight());
  fhEMCALPHOSCorrNCells   ->Fill(ncellsEMCAL          , ncellsPHOS          , GetEventWeight());
  fhEMCALPHOSCorrECells   ->Fill(sumCellEnergyEMCAL   , sumCellEnergyPHOS   , GetEventWeight());
  
  fhEMCALDCALCorrNClusters->Fill(nclEMCAL             , nclDCAL             , GetEventWeight());
  fhEMCALDCALCorrEClusters->Fill(sumClusterEnergyEMCAL, sumClusterEnergyDCAL, GetEventWeight());
  fhEMCALDCALCorrNCells   ->Fill(ncellsEMCAL          , ncellsDCAL          , GetEventWeight());
  fhEMCALDCALCorrECells   ->Fill(sumCellEnergyEMCAL   , sumCellEnergyDCAL   , GetEventWeight());
  
  fhDCALPHOSCorrNClusters ->Fill(nclDCAL              , nclPHOS             , GetEventWeight());
  fhDCALPHOSCorrEClusters ->Fill(sumClusterEnergyDCAL , sumClusterEnergyPHOS, GetEventWeight());
  fhDCALPHOSCorrNCells    ->Fill(ncellsDCAL           , ncellsPHOS          , GetEventWeight());
  fhDCALPHOSCorrECells    ->Fill(sumCellEnergyDCAL    , sumCellEnergyPHOS   , GetEventWeight());
  
  Int_t   v0S = GetV0Signal(0)       + GetV0Signal(1);
  Int_t   v0M = GetV0Multiplicity(0) + GetV0Multiplicity(1);
  Int_t   trM = GetTrackMultiplicity();
  Float_t cen = GetEventCentrality();
  Float_t ep  = GetEventPlaneAngle();
  
  Int_t   ncl              = nclPHOS;
  Float_t sumClusterEnergy = sumClusterEnergyPHOS;
  Int_t   ncells           = ncellsPHOS;
  Float_t sumCellEnergy    = sumCellEnergyPHOS;
  
  if ( GetCalorimeter() == kEMCAL )
  {
    ncl              = nclEMCAL              + nclDCAL;
    sumClusterEnergy = sumClusterEnergyEMCAL + sumClusterEnergyDCAL;
    ncells           = ncellsEMCAL           + ncellsDCAL;
    sumCellEnergy    = sumCellEnergyEMCAL    + sumCellEnergyDCAL;
  }
  
  fhCaloV0MCorrNClusters   ->Fill(v0M, ncl             , GetEventWeight());
  fhCaloV0MCorrEClusters   ->Fill(v0M, sumClusterEnergy, GetEventWeight());
  fhCaloV0MCorrNCells      ->Fill(v0M, ncells          , GetEventWeight());
  fhCaloV0MCorrECells      ->Fill(v0M, sumCellEnergy   , GetEventWeight());
  
  fhCaloV0SCorrNClusters   ->Fill(v0S, ncl             , GetEventWeight());
  fhCaloV0SCorrEClusters   ->Fill(v0S, sumClusterEnergy, GetEventWeight());
  fhCaloV0SCorrNCells      ->Fill(v0S, ncells          , GetEventWeight());
  fhCaloV0SCorrECells      ->Fill(v0S, sumCellEnergy   , GetEventWeight());
  
  fhCaloTrackMCorrNClusters->Fill(trM, ncl             , GetEventWeight());
  fhCaloTrackMCorrEClusters->Fill(trM, sumClusterEnergy, GetEventWeight());
  fhCaloTrackMCorrNCells   ->Fill(trM, ncells          , GetEventWeight());
  fhCaloTrackMCorrECells   ->Fill(trM, sumCellEnergy   , GetEventWeight());
  
  fhCaloCenNClusters       ->Fill(cen, ncl             , GetEventWeight());
  fhCaloCenEClusters       ->Fill(cen, sumClusterEnergy, GetEventWeight());
  fhCaloCenNCells          ->Fill(cen, ncells          , GetEventWeight());
  fhCaloCenECells          ->Fill(cen, sumCellEnergy   , GetEventWeight());
  
  fhCaloEvPNClusters       ->Fill(ep , ncl             , GetEventWeight());
  fhCaloEvPEClusters       ->Fill(ep , sumClusterEnergy, GetEventWeight());
  fhCaloEvPNCells          ->Fill(ep , ncells          , GetEventWeight());
  fhCaloEvPECells          ->Fill(ep , sumCellEnergy   , GetEventWeight());
  
  AliDebug(1,"Correlate():");
  AliDebug(1,Form("\t EMCAL: N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f",
                  ncellsEMCAL,nclEMCAL, sumCellEnergyEMCAL,sumClusterEnergyEMCAL));  
  AliDebug(1,Form("\t DCAL : N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f",
                  ncellsDCAL,nclDCAL, sumCellEnergyDCAL,sumClusterEnergyDCAL));
  AliDebug(1,Form("\t PHOS : N cells %d, N clusters  %d, summed E cells %f, summed E clusters %f",
                  ncellsPHOS,nclPHOS,sumCellEnergyPHOS,sumClusterEnergyPHOS));
  AliDebug(1,Form("\t V0 : Signal %d, Multiplicity  %d, Track Multiplicity %d", v0S,v0M,trM));
  AliDebug(1,Form("\t centrality : %f, Event plane angle %f", cen,ep));
}

//_________________________________________________
/// Save parameters used for analysis in a string.
//_________________________________________________
TObjString * AliAnaCalorimeterQA::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaCalorimeterQA ---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Time Cut : %2.2f < T < %2.2f ns;",fTimeCutMin, fTimeCutMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Cell Amplitude: PHOS > %2.2f GeV, EMCAL > %2.2f GeV;",fPHOSCellAmpMin, fEMCALCellAmpMin) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"N cells per cluster: PHOS >= %d, EMCAL >= %d;",fPHOSClusterNCellMin, fEMCALClusterNCellMin) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Inv. Mass   %2.1f < E_cl < %2.1f GeV;",fInvMassMinECut, fInvMassMaxECut) ;
  parList+=onePar ;    
  snprintf(onePar,buffersize,"Inv. Mass   %2.1f < M02 < %2.1f GeV;",fInvMassMinM02Cut, fInvMassMaxM02Cut) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"Cluster pair opening angle < %2.1f rad;",fInvMassMaxOpenAngle) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"Cluster pair time difference < %2.1f rad;",fInvMassMaxTimeDifference) ;
  parList+=onePar ;

  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
	
  return new TObjString(parList) ;
}

//_________________________________________________________________________________
/// Fill histograms with exoticity parameters.
//_________________________________________________________________________________
void AliAnaCalorimeterQA::ExoticHistograms(Int_t absIdMax, Float_t ampMax,
                                           AliVCluster *clus, AliVCaloCells* cells)
{
  if(ampMax < 0.01) 
  {
    AliDebug(1,Form("Low amplitude energy %f",ampMax));
    return;
  }
    
  Float_t  l0   = clus->GetM02();
  Float_t  l1   = clus->GetM20();
  Float_t  en   = clus->E();
  Int_t    nc   = clus->GetNCells();  
  Double_t tmax = clus->GetTOF()*1.e9-fConstantTimeShift; // recalibrated elsewhere
  
  Float_t eCrossFrac = 1-GetECross(absIdMax,cells, 10000000)/ampMax;

  if(en > 5) 
  {
    fhExoL0ECross->Fill(eCrossFrac, l0, GetEventWeight());
    fhExoL1ECross->Fill(eCrossFrac, l1, GetEventWeight());
  }
  
  for(Int_t ie = 0; ie < fExoNECrossCuts; ie++)
  {    
    for(Int_t idt = 0; idt < fExoNDTimeCuts; idt++)
    {  
      eCrossFrac = 1-GetECross(absIdMax,cells, fExoDTimeCuts[idt])/ampMax;
      
      if(eCrossFrac > fExoECrossCuts[ie])
      {
        // Exotic
        fhExoL0    [ie][idt]->Fill(en,l0  , GetEventWeight());
        fhExoL1    [ie][idt]->Fill(en,l1  , GetEventWeight());
        fhExoTime  [ie][idt]->Fill(en,tmax, GetEventWeight());
        
        if(en > 5) 
        {
          fhExoL0NCell[ie][idt]->Fill(nc, l0, GetEventWeight());
          fhExoL1NCell[ie][idt]->Fill(nc, l1, GetEventWeight());
        } 
        
        // Diff time, do for one cut in time
        if(idt == 0)
        {
          for (Int_t icell = 0; icell < clus->GetNCells(); icell++) 
          {
            Int_t absId  = clus->GetCellsAbsId()[icell]; 
            Double_t time  = cells->GetCellTime(absId);
            GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
            time *= 1e9;
            time -= fConstantTimeShift;
            
            Float_t diff = tmax-time;
            fhExoDTime[ie]->Fill(en, diff, GetEventWeight());
          }
        }
      }
      else
      {
        fhExoECross[ie][idt]->Fill(en, eCrossFrac, GetEventWeight());
        fhExoNCell [ie][idt]->Fill(en, nc        , GetEventWeight());
      }
    } // D time cut loop
  } // e cross cut loop
}

//___________________________________________________
/// Create histograms to be saved in output file and
/// store them in the output container.
//___________________________________________________
TList * AliAnaCalorimeterQA::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("QAHistos") ; 
  
  // Init the number of modules, set in the class AliCalorimeterUtils
    
  fNModules = GetCaloUtils()->GetNumberOfSuperModulesUsed();
  if(GetCalorimeter()==kPHOS && fNModules > 4) fNModules = 4;
  
  // Histogram binning and ranges
    
  Int_t nptbins     = GetHistogramRanges()->GetHistoPtBins(); 	        Float_t ptmax     = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin     = GetHistogramRanges()->GetHistoPtMin();
  Int_t nfineptbins = GetHistogramRanges()->GetHistoFinePtBins(); 	    Float_t ptfinemax = GetHistogramRanges()->GetHistoFinePtMax();       Float_t ptfinemin = GetHistogramRanges()->GetHistoFinePtMin();
  Int_t nphibins    = GetHistogramRanges()->GetHistoPhiBins();     	    Float_t phimax    = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin    = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins    = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax    = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin    = GetHistogramRanges()->GetHistoEtaMin();	
  Int_t nmassbins   = GetHistogramRanges()->GetHistoMassBins();         Float_t massmax   = GetHistogramRanges()->GetHistoMassMax(); 	       Float_t massmin   = GetHistogramRanges()->GetHistoMassMin();
  Int_t nasymbins   = GetHistogramRanges()->GetHistoAsymmetryBins();    Float_t asymmax   = GetHistogramRanges()->GetHistoAsymmetryMax();    Float_t asymmin   = GetHistogramRanges()->GetHistoAsymmetryMin();
  Int_t nPoverEbins = GetHistogramRanges()->GetHistoPOverEBins();       Float_t eOverPmax = GetHistogramRanges()->GetHistoPOverEMax();       Float_t eOverPmin = GetHistogramRanges()->GetHistoPOverEMin();
  Int_t ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();         Float_t dedxmax   = GetHistogramRanges()->GetHistodEdxMax();         Float_t dedxmin   = GetHistogramRanges()->GetHistodEdxMin();
  Int_t ndRbins     = GetHistogramRanges()->GetHistodRBins();           Float_t dRmax     = GetHistogramRanges()->GetHistodRMax();           Float_t dRmin     = GetHistogramRanges()->GetHistodRMin();
  Int_t ntimebins   = GetHistogramRanges()->GetHistoTimeBins();         Float_t timemax   = GetHistogramRanges()->GetHistoTimeMax();         Float_t timemin   = GetHistogramRanges()->GetHistoTimeMin();       
  Int_t nclbins     = GetHistogramRanges()->GetHistoNClustersBins();    Int_t   nclmax    = GetHistogramRanges()->GetHistoNClustersMax();    Int_t   nclmin    = GetHistogramRanges()->GetHistoNClustersMin(); 
  Int_t ncebins     = GetHistogramRanges()->GetHistoNCellsBins();       Int_t   ncemax    = GetHistogramRanges()->GetHistoNCellsMax();       Int_t   ncemin    = GetHistogramRanges()->GetHistoNCellsMin(); 
  Int_t nceclbins   = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nceclmax  = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nceclmin  = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  Int_t nvdistbins  = GetHistogramRanges()->GetHistoVertexDistBins();   Float_t vdistmax  = GetHistogramRanges()->GetHistoVertexDistMax();   Float_t vdistmin  = GetHistogramRanges()->GetHistoVertexDistMin();
  Int_t rbins       = GetHistogramRanges()->GetHistoRBins();            Float_t rmax      = GetHistogramRanges()->GetHistoRMax();            Float_t rmin      = GetHistogramRanges()->GetHistoRMin(); 
  Int_t xbins       = GetHistogramRanges()->GetHistoXBins();            Float_t xmax      = GetHistogramRanges()->GetHistoXMax();            Float_t xmin      = GetHistogramRanges()->GetHistoXMin(); 
  Int_t ybins       = GetHistogramRanges()->GetHistoYBins();            Float_t ymax      = GetHistogramRanges()->GetHistoYMax();            Float_t ymin      = GetHistogramRanges()->GetHistoYMin(); 
  Int_t zbins       = GetHistogramRanges()->GetHistoZBins();            Float_t zmax      = GetHistogramRanges()->GetHistoZMax();            Float_t zmin      = GetHistogramRanges()->GetHistoZMin(); 
  Int_t ssbins      = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax     = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin     = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t tdbins      = GetHistogramRanges()->GetHistoDiffTimeBins() ;    Float_t tdmax     = GetHistogramRanges()->GetHistoDiffTimeMax();     Float_t tdmin     = GetHistogramRanges()->GetHistoDiffTimeMin();
  
  Int_t nv0sbins    = GetHistogramRanges()->GetHistoV0SignalBins();          Int_t nv0smax = GetHistogramRanges()->GetHistoV0SignalMax();          Int_t nv0smin = GetHistogramRanges()->GetHistoV0SignalMin(); 
  Int_t nv0mbins    = GetHistogramRanges()->GetHistoV0MultiplicityBins();    Int_t nv0mmax = GetHistogramRanges()->GetHistoV0MultiplicityMax();    Int_t nv0mmin = GetHistogramRanges()->GetHistoV0MultiplicityMin(); 
  Int_t ntrmbins    = GetHistogramRanges()->GetHistoTrackMultiplicityBins(); Int_t ntrmmax = GetHistogramRanges()->GetHistoTrackMultiplicityMax(); Int_t ntrmmin = GetHistogramRanges()->GetHistoTrackMultiplicityMin(); 
  
  // Calorimeter grid
    
  // EMCAL
  fNMaxCols = 48;
  fNMaxRows = 24;
  fNRCU     = 2 ;
  // PHOS
  if(GetCalorimeter()==kPHOS)
  {
    fNMaxCols = 56;
    fNMaxRows = 64;
    fNRCU     = 4 ;
  }
  
    
  // Init histograms
    
  fhE  = new TH1F ("hE","#it{E} reconstructed clusters ", nptbins*5,ptmin,ptmax*5);  
  fhE->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhE);
  
  fhPt  = new TH1F ("hPt","#it{p}_{T} reconstructed clusters", nptbins,ptmin,ptmax);
  fhPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPt);
  
  fhPhi  = new TH1F ("hPhi","#phi reconstructed clusters ",nphibins,phimin,phimax);
  fhPhi->SetXTitle("#phi (rad)");
  outputContainer->Add(fhPhi);
  
  fhEta  = new TH1F ("hEta","#eta reconstructed clusters ",netabins,etamin,etamax);
  fhEta->SetXTitle("#eta ");
  outputContainer->Add(fhEta);
  
  if(fFillAllTH3)
  {
    fhEtaPhiE  = new TH3F ("hEtaPhiE","#eta vs #phi vs energy, reconstructed clusters",
                           netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
    fhEtaPhiE->SetXTitle("#eta ");
    fhEtaPhiE->SetYTitle("#phi (rad)");
    fhEtaPhiE->SetZTitle("#it{E} (GeV) ");
    outputContainer->Add(fhEtaPhiE);
  }
  else
  {
    fhEtaPhi  = new TH2F ("hEtaPhi","#eta vs #phi for #it{E} > 0.5 GeV, reconstructed clusters",
                           netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi->SetXTitle("#eta ");
    fhEtaPhi->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhi);
  }
  
  fhClusterTimeEnergy  = new TH2F ("hClusterTimeEnergy","energy vs TOF, reconstructed clusters",
                                   nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
  fhClusterTimeEnergy->SetXTitle("#it{E} (GeV) ");
  fhClusterTimeEnergy->SetYTitle("TOF (ns)");
  outputContainer->Add(fhClusterTimeEnergy);
  
  fhClusterPairDiffTimeE = new TH2F("hClusterPairDiffTimeE","cluster pair time difference vs E, only good clusters",
                                    nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhClusterPairDiffTimeE->SetXTitle("#it{E}_{cluster} (GeV)");
  fhClusterPairDiffTimeE->SetYTitle("#Delta #it{t} (ns)");
  outputContainer->Add(fhClusterPairDiffTimeE);  

  fhClusterPairDiffTimeESameMod = new TH2F("hClusterPairDiffTimeESameMod","cluster pair time difference vs E, only good clusters",
                                    nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
  fhClusterPairDiffTimeESameMod->SetXTitle("#it{E}_{cluster} (GeV)");
  fhClusterPairDiffTimeESameMod->SetYTitle("#Delta #it{t} (ns)");
  outputContainer->Add(fhClusterPairDiffTimeESameMod);  
  
  fhLambda0  = new TH2F ("hLambda0","shower shape, #lambda^{2}_{0} vs E",
                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhLambda0->SetXTitle("#it{E}_{cluster}");
  fhLambda0->SetYTitle("#lambda^{2}_{0}");
  outputContainer->Add(fhLambda0); 
  
  fhLambda1  = new TH2F ("hLambda1","shower shape, #lambda^{2}_{1} vs E for bad cluster ",
                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
  fhLambda1->SetXTitle("#it{E}_{cluster}");
  fhLambda1->SetYTitle("#lambda^{2}_{1}");
  outputContainer->Add(fhLambda1); 
  
//  fhDispersion  = new TH2F ("hDispersion","shower shape, Dispersion^{2} vs E for bad cluster ",
//                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
//  fhDispersion->SetXTitle("#it{E}_{cluster}");
//  fhDispersion->SetYTitle("Dispersion");
//  outputContainer->Add(fhDispersion);       
  
  fhClusterMaxCellCloseCellRatio  = new TH2F ("hClusterMaxCellCloseCellRatio","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                              nptbins,ptmin,ptmax, 100,0,1.); 
  fhClusterMaxCellCloseCellRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhClusterMaxCellCloseCellRatio->SetYTitle("#it{E}_{cell i}/#it{E}_{cell max}");
  outputContainer->Add(fhClusterMaxCellCloseCellRatio);
  
  fhClusterMaxCellCloseCellDiff  = new TH2F ("hClusterMaxCellCloseCellDiff","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                             nptbins,ptmin,ptmax, 500,0,100.); 
  fhClusterMaxCellCloseCellDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhClusterMaxCellCloseCellDiff->SetYTitle("#it{E}_{cell max}-#it{E}_{cell i} (GeV)");
  outputContainer->Add(fhClusterMaxCellCloseCellDiff);
  
  fhClusterMaxCellDiff  = new TH2F ("hClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                    nptbins,ptmin,ptmax, 500,0,1.); 
  fhClusterMaxCellDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhClusterMaxCellDiff->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
  outputContainer->Add(fhClusterMaxCellDiff);  
  
  if(fStudyBadClusters)
  {
    fhClusterMaxCellDiffNoCut  = new TH2F ("hClusterMaxCellDiffNoCut","energy vs difference of cluster energy - max cell energy / cluster energy",
                                           nptbins,ptmin,ptmax, 500,0,1.); 
    fhClusterMaxCellDiffNoCut->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhClusterMaxCellDiffNoCut->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
    outputContainer->Add(fhClusterMaxCellDiffNoCut);  
  }
  
  fhClusterMaxCellECross  = new TH2F ("hClusterMaxCellECross","1 - Energy in cross around max energy cell / max energy cell vs cluster energy, good clusters",
                                       nptbins,ptmin,ptmax, 400,-1,1.); 
  fhClusterMaxCellECross->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhClusterMaxCellECross->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell max}");
  outputContainer->Add(fhClusterMaxCellECross);    
  
  if(fStudyBadClusters)
  {
    fhNCellsPerClusterNoCut  = new TH2F ("hNCellsPerClusterNoCut","# cells per cluster vs energy, no bad clusters cut",
                                         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerClusterNoCut->SetXTitle("#it{E} (GeV)");
    fhNCellsPerClusterNoCut->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterNoCut);
  }
  
  fhNCellsPerCluster  = new TH2F ("hNCellsPerCluster","# cells per cluster vs energy",nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerCluster->SetXTitle("#it{E} (GeV)");
  fhNCellsPerCluster->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerCluster);
  
  if(fStudyBadClusters)
  {
    fhNCellsPerClusterWeirdModNoCut  = new TH2F ("hNCellsPerClusterWeirdNoCutMod","# cells per cluster vs energy, E > 100, no bad clusters cut, per SM",
                                                 nceclbins,nceclmin,nceclmax,fNModules,0,fNModules); 
    fhNCellsPerClusterWeirdModNoCut->SetYTitle("SM number");
    fhNCellsPerClusterWeirdModNoCut->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterWeirdModNoCut);
  }
  
  fhNCellsPerClusterWeirdMod  = new TH2F ("hNCellsPerClusterWeirdMod","# cells per cluster, E > 100 GeV, per SM", 
                                          nceclbins*2,nceclmin,nceclmax*2,fNModules,0,fNModules); 
  fhNCellsPerClusterWeirdMod->SetYTitle("SM number");
  fhNCellsPerClusterWeirdMod->SetXTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerClusterWeirdMod);
  
  fhNClusters  = new TH1F ("hNClusters","# clusters", nclbins,nclmin,nclmax); 
  fhNClusters->SetXTitle("#it{n}_{clusters}");
  outputContainer->Add(fhNClusters);

  if(fStudyBadClusters)
  {
    fhBadClusterEnergy  = new TH1F ("hBadClusterEnergy","Bad cluster energy", nptbins,ptmin,ptmax); 
    fhBadClusterEnergy->SetXTitle("#it{E}_{cluster} (GeV) ");
    outputContainer->Add(fhBadClusterEnergy);

    fhBadClusterEtaPhi  = new TH2F ("hBadClusterEtaPhi","Bad cluster, #eta vs #phi, #it{E} > 0.5 GeV", 
                                    netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhBadClusterEtaPhi->SetXTitle("#eta ");
    fhBadClusterEtaPhi->SetXTitle("#phi (rad) ");
    outputContainer->Add(fhBadClusterEtaPhi);
    
    fhBadClusterLambda0  = new TH2F ("hBadClusterLambda0","Bad cluster,shower shape, #lambda^{2}_{0} vs E",
                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhBadClusterLambda0->SetXTitle("#it{E}_{cluster}");
    fhBadClusterLambda0->SetYTitle("#lambda^{2}_{0}");
    outputContainer->Add(fhBadClusterLambda0); 
    
    fhBadClusterLambda1  = new TH2F ("hBadClusterLambda1","Bad cluster,shower shape, #lambda^{2}_{1} vs E for bad cluster ",
                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhBadClusterLambda1->SetXTitle("#it{E}_{cluster}");
    fhBadClusterLambda1->SetYTitle("#lambda^{2}_{1}");
    outputContainer->Add(fhBadClusterLambda1); 
    
    fhBadClusterMaxCellCloseCellRatio  = new TH2F ("hBadClusterMaxCellCloseCellRatio","energy vs ratio of max cell / neighbour cell constributing cell, reconstructed bad clusters",
                                                   nptbins,ptmin,ptmax, 100,0,1.); 
    fhBadClusterMaxCellCloseCellRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhBadClusterMaxCellCloseCellRatio->SetYTitle("ratio");
    outputContainer->Add(fhBadClusterMaxCellCloseCellRatio);
    
    fhBadClusterMaxCellCloseCellDiff  = new TH2F ("hBadClusterMaxCellCloseCellDiff","energy vs ratio of max cell - neighbour cell constributing cell, reconstructed bad clusters",
                                                  nptbins,ptmin,ptmax, 500,0,100); 
    fhBadClusterMaxCellCloseCellDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhBadClusterMaxCellCloseCellDiff->SetYTitle("#it{E}_{cell max} - #it{E}_{cell i} (GeV)");
    outputContainer->Add(fhBadClusterMaxCellCloseCellDiff);    
    
    fhBadClusterMaxCellDiff  = new TH2F ("hBadClusterMaxCellDiff","energy vs difference of cluster energy - max cell energy / cluster energy for bad clusters",
                                         nptbins,ptmin,ptmax, 500,0,1.); 
    fhBadClusterMaxCellDiff->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhBadClusterMaxCellDiff->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max}) / #it{E}_{cluster}");
    outputContainer->Add(fhBadClusterMaxCellDiff);
    
    fhBadClusterTimeEnergy  = new TH2F ("hBadClusterTimeEnergy","energy vs TOF of reconstructed bad clusters",
                                        nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhBadClusterTimeEnergy->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhBadClusterTimeEnergy->SetYTitle("#it{t} (ns)");
    outputContainer->Add(fhBadClusterTimeEnergy);    
    
    fhBadClusterPairDiffTimeE = new TH2F("hBadClusterPairDiffTimeE","cluster pair time difference (bad - good) vs E from bad cluster",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
    fhBadClusterPairDiffTimeE->SetXTitle("#it{E}_{bad cluster} (GeV)");
    fhBadClusterPairDiffTimeE->SetYTitle("#Delta #it{t} (ns)");
    outputContainer->Add(fhBadClusterPairDiffTimeE);    
    
    fhBadClusterMaxCellECross  = new TH2F ("hBadClusterMaxCellECross","1 - #it{E}_{+} around max energy cell / max energy cell vs cluster energy, bad clusters",
                                        nptbins,ptmin,ptmax, 400,-1,1.); 
    fhBadClusterMaxCellECross->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhBadClusterMaxCellECross->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell max}");
    outputContainer->Add(fhBadClusterMaxCellECross);        
    
    if(fFillAllCellTimeHisto) 
    {
      fhBadCellTimeSpreadRespectToCellMax = new TH2F ("hBadCellTimeSpreadRespectToCellMax","#it{t}_{cell max}-#it{t}_{cell i} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
      fhBadCellTimeSpreadRespectToCellMax->SetXTitle("#it{E} (GeV)");
      fhBadCellTimeSpreadRespectToCellMax->SetYTitle("#Delta #it{t}_{cell max - i} (ns)");
      outputContainer->Add(fhBadCellTimeSpreadRespectToCellMax);
      
      fhBadClusterMaxCellDiffAverageTime = new TH2F ("hBadClusterMaxCellDiffAverageTime","#it{t}_{cell max}-#it{t}_{average} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
      fhBadClusterMaxCellDiffAverageTime->SetXTitle("#it{E} (GeV)");
      fhBadClusterMaxCellDiffAverageTime->SetYTitle("#Delta #it{t}_{cell max - average} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffAverageTime);
            
      fhBadClusterMaxCellDiffWeightedTime = new TH2F ("hBadClusterMaxCellDiffWeightedTime","#it{t}_{cell max}-#it{t}_{weighted} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
      fhBadClusterMaxCellDiffWeightedTime->SetXTitle("#it{E} (GeV)");
      fhBadClusterMaxCellDiffWeightedTime->SetYTitle("#Delta #it{t}_{cell max - weighted} (ns)");
      outputContainer->Add(fhBadClusterMaxCellDiffWeightedTime);
      
    }  
    
  }
  
  if(fStudyExotic)
  {
    fhExoL0ECross  = new TH2F("hExoL0_ECross",
                               "#lambda^{2}_{0} vs 1-#it{E}_{+}/#it{E}_{max} for E > 5 GeV",
                               400,0,1,ssbins,ssmin,ssmax); 
    fhExoL0ECross ->SetXTitle("1-#it{E}_{+}/#it{E}_{cell max}");
    fhExoL0ECross ->SetYTitle("#lambda^{2}_{0}");
    outputContainer->Add(fhExoL0ECross) ;     

    fhExoL1ECross  = new TH2F("hExoL1_ECross",
                              "#lambda^{2}_{1} vs 1-#it{E}_{+}/#it{E}_{max} for E > 5 GeV",
                              400,0,1,ssbins,ssmin,ssmax); 
    fhExoL1ECross ->SetXTitle("1-#it{E}_{+}/#it{E}_{cell max}");
    fhExoL1ECross ->SetYTitle("#lambda^{2}_{1}");
    outputContainer->Add(fhExoL1ECross) ;  
    
    for(Int_t ie = 0; ie <fExoNECrossCuts; ie++)
    {  
      
      fhExoDTime[ie]  = new TH2F(Form("hExoDTime_ECross%d",ie),
                                 Form("#Delta time = t_{max}-t_{cells} vs #it{E}_{cluster} for exotic, 1-#it{E}_{+}/#it{E}_{max} < %2.2f",fExoECrossCuts[ie]),
                                 nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
      fhExoDTime[ie] ->SetYTitle("#Delta #it{t} (ns)");
      fhExoDTime[ie] ->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhExoDTime[ie]) ; 
      
      for(Int_t idt = 0; idt < fExoNDTimeCuts; idt++)
      {        
        fhExoNCell[ie][idt]  = new TH2F(Form("hExoNCell_ECross%d_DT%d",ie,idt),
                                     Form("N cells per cluster vs E cluster, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t < %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                     nptbins,ptmin,ptmax,nceclbins,nceclmin,nceclmax); 
        fhExoNCell[ie][idt] ->SetYTitle("#it{n}_cells");
        fhExoNCell[ie][idt] ->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhExoNCell[ie][idt]) ; 
        
        fhExoL0   [ie][idt]  = new TH2F(Form("hExoL0_ECross%d_DT%d",ie,idt),
                                     Form("#lambda^{2}_{0} vs E cluster for exotic, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t = %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhExoL0   [ie][idt] ->SetYTitle("#lambda^{2}_{0}");
        fhExoL0   [ie][idt] ->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhExoL0[ie][idt]) ; 

        fhExoL1   [ie][idt]  = new TH2F(Form("hExoL1_ECross%d_DT%d",ie,idt),
                                        Form("#lambda^{2}_{1} vs E cluster for exotic, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t = %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhExoL1   [ie][idt] ->SetYTitle("#lambda^{2}_{1}");
        fhExoL1   [ie][idt] ->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhExoL1[ie][idt]) ; 
        
        fhExoECross[ie][idt]  = new TH2F(Form("hExoECross_ECross%d_DT%d",ie,idt),
                                      Form("#it{E} cross for cells vs E cell, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t < %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                      nptbins,ptmin,ptmax,400,0,1); 
        fhExoECross[ie][idt] ->SetYTitle("1-#it{E}_{+}/#it{E}_{cell max}");
        fhExoECross[ie][idt] ->SetXTitle("#it{E}_{cell} (GeV)");
        outputContainer->Add(fhExoECross[ie][idt]) ; 
        
        fhExoTime  [ie][idt]  = new TH2F(Form("hExoTime_ECross%d_DT%d",ie,idt),
                                        Form("Time of cluster (max cell) vs E cluster for exotic, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t = %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                        nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhExoTime  [ie][idt] ->SetYTitle("#it{t}_{max} (ns)");
        fhExoTime  [ie][idt] ->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhExoTime[ie][idt]) ; 

        fhExoL0NCell[ie][idt]  = new TH2F(Form("hExoL0_NCell%d_DT%d",ie,idt),
                                          Form("#lambda^{2}_{0} vs N cells per clusters for E > 5 GeV, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t = %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                          nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhExoL0NCell[ie][idt] ->SetYTitle("#it{n}_{cells}");
        fhExoL0NCell[ie][idt] ->SetXTitle("#lambda^{2}_{0}");
        outputContainer->Add(fhExoL0NCell[ie][idt]) ;  
        
        fhExoL1NCell[ie][idt]  = new TH2F(Form("hExoL1_NCell%d_DT%d",ie,idt),
                                          Form("#lambda^{2}_{1} vs N cells per clusters for E > 5 GeV, 1-#it{E}_{+}/#it{E}_{max} < %2.2f, #Delta t = %2.0f",fExoECrossCuts[ie],fExoDTimeCuts[idt]),
                                          nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhExoL1NCell[ie][idt] ->SetYTitle("#it{n}_{cells}");
        fhExoL1NCell[ie][idt] ->SetXTitle("#lambda^{2}_{1}");
        outputContainer->Add(fhExoL1NCell[ie][idt]) ;  
      }
    } 
  }
  
  // Cluster size in terms of cells
    
  if(fStudyClustersAsymmetry)
  {
    fhDeltaIEtaDeltaIPhiE0[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE0"," Cluster size in columns vs rows for E < 2 GeV, #it{n}_{cells} > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE0[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE0[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE0[0]); 
    
    fhDeltaIEtaDeltaIPhiE2[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE2"," Cluster size in columns vs rows for 2 <E < 6 GeV, #it{n}_{cells} > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE2[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE2[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE2[0]); 
    
    fhDeltaIEtaDeltaIPhiE6[0]  = new TH2F ("hDeltaIEtaDeltaIPhiE6"," Cluster size in columns vs rows for E > 6 GeV, #it{n}_{cells} > 3",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE6[0]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE6[0]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE6[0]); 
    
    fhDeltaIA[0]  = new TH2F ("hDeltaIA"," Cluster *asymmetry* in cell units vs E",
                              nptbins,ptmin,ptmax,21,-1.05,1.05); 
    fhDeltaIA[0]->SetXTitle("#it{E}_{cluster}");
    fhDeltaIA[0]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIA[0]); 
    
    fhDeltaIAL0[0]  = new TH2F ("hDeltaIAL0"," Cluster *asymmetry* in cell units vs #lambda^{2}_{0}",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL0[0]->SetXTitle("#lambda^{2}_{0}");
    fhDeltaIAL0[0]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL0[0]); 
    
    fhDeltaIAL1[0]  = new TH2F ("hDeltaIAL1"," Cluster *asymmetry* in cell units vs #lambda^{2}_{1}",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL1[0]->SetXTitle("#lambda^{2}_{1}");
    fhDeltaIAL1[0]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL1[0]); 
    
    fhDeltaIANCells[0]  = new TH2F ("hDeltaIANCells"," Cluster *asymmetry* in cell units vs N cells in cluster",
                                    nceclbins,nceclmin,nceclmax,21,-1.05,1.05); 
    fhDeltaIANCells[0]->SetXTitle("#it{n}_{cell in cluster}");
    fhDeltaIANCells[0]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIANCells[0]); 
    
    
    fhDeltaIEtaDeltaIPhiE0[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE0Charged"," Cluster size in columns vs rows for E < 2 GeV, #it{n}_{cells} > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE0[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE0[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE0[1]); 
    
    fhDeltaIEtaDeltaIPhiE2[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE2Charged"," Cluster size in columns vs rows for 2 <E < 6 GeV, #it{n}_{cells} > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE2[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE2[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE2[1]); 
    
    fhDeltaIEtaDeltaIPhiE6[1]  = new TH2F ("hDeltaIEtaDeltaIPhiE6Charged"," Cluster size in columns vs rows for E > 6 GeV, #it{n}_{cells} > 3, matched with track",
                                           50,0,50,50,0,50); 
    fhDeltaIEtaDeltaIPhiE6[1]->SetXTitle("#Delta Column");
    fhDeltaIEtaDeltaIPhiE6[1]->SetYTitle("#Delta Row");
    outputContainer->Add(fhDeltaIEtaDeltaIPhiE6[1]); 
    
    fhDeltaIA[1]  = new TH2F ("hDeltaIACharged"," Cluster *asymmetry* in cell units vs E, matched with track",
                              nptbins,ptmin,ptmax,21,-1.05,1.05); 
    fhDeltaIA[1]->SetXTitle("#it{E}_{cluster}");
    fhDeltaIA[1]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIA[1]); 
    
    fhDeltaIAL0[1]  = new TH2F ("hDeltaIAL0Charged"," Cluster *asymmetry* in cell units vs #lambda^{2}_{0}, matched with track",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL0[1]->SetXTitle("#lambda^{2}_{0}");
    fhDeltaIAL0[1]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL0[1]); 
    
    fhDeltaIAL1[1]  = new TH2F ("hDeltaIAL1Charged"," Cluster *asymmetry* in cell units vs #lambda^{2}_{1}, matched with track",
                                ssbins,ssmin,ssmax,21,-1.05,1.05); 
    fhDeltaIAL1[1]->SetXTitle("#lambda^{2}_{1}");
    fhDeltaIAL1[1]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIAL1[1]); 
    
    fhDeltaIANCells[1]  = new TH2F ("hDeltaIANCellsCharged"," Cluster *asymmetry* in cell units vs N cells in cluster, matched with track",
                                    nceclbins,nceclmin,nceclmax,21,-1.05,1.05); 
    fhDeltaIANCells[1]->SetXTitle("#it{n}_{cell in cluster}");
    fhDeltaIANCells[1]->SetYTitle("#it{A}_{cell in cluster}");
    outputContainer->Add(fhDeltaIANCells[1]); 
    
    if(IsDataMC()){
      TString particle[]={"Photon","Electron","Conversion","Hadron"};
      for (Int_t iPart = 0; iPart < 4; iPart++) {
        
        fhDeltaIAMC[iPart]  = new TH2F (Form("hDeltaIA_MC%s",particle[iPart].Data()),Form(" Cluster *asymmetry* in cell units vs E, from %s",particle[iPart].Data()),
                                        nptbins,ptmin,ptmax,21,-1.05,1.05); 
        fhDeltaIAMC[iPart]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIAMC[iPart]->SetYTitle("#it{A}_{cell in cluster}");
        outputContainer->Add(fhDeltaIAMC[iPart]);     
      }
    }
    
    if(fStudyBadClusters)
    {
      fhBadClusterDeltaIEtaDeltaIPhiE0  = new TH2F ("hBadClusterDeltaIEtaDeltaIPhiE0"," Cluster size in columns vs rows for E < 2 GeV, #it{n}_{cells} > 3",
                                                    50,0,50,50,0,50); 
      fhBadClusterDeltaIEtaDeltaIPhiE0->SetXTitle("#Delta Column");
      fhBadClusterDeltaIEtaDeltaIPhiE0->SetYTitle("#Delta Row");
      outputContainer->Add(fhBadClusterDeltaIEtaDeltaIPhiE0); 
      
      fhBadClusterDeltaIEtaDeltaIPhiE2  = new TH2F ("hBadClusterDeltaIEtaDeltaIPhiE2"," Cluster size in columns vs rows for 2 <E < 6 GeV, #it{n}_{cells} > 3",
                                                    50,0,50,50,0,50); 
      fhBadClusterDeltaIEtaDeltaIPhiE2->SetXTitle("#Delta Column");
      fhBadClusterDeltaIEtaDeltaIPhiE2->SetYTitle("#Delta Row");
      outputContainer->Add(fhBadClusterDeltaIEtaDeltaIPhiE2); 
      
      fhBadClusterDeltaIEtaDeltaIPhiE6  = new TH2F ("hBadClusterDeltaIEtaDeltaIPhiE6"," Cluster size in columns vs rows for E > 6 GeV, #it{n}_{cells} > 3",
                                                    50,0,50,50,0,50); 
      fhBadClusterDeltaIEtaDeltaIPhiE6->SetXTitle("#Delta Column");
      fhBadClusterDeltaIEtaDeltaIPhiE6->SetYTitle("#Delta Row");
      outputContainer->Add(fhBadClusterDeltaIEtaDeltaIPhiE6); 
      
      fhBadClusterDeltaIA  = new TH2F ("hBadClusterDeltaIA"," Cluster *asymmetry* in cell units vs E",
                                       nptbins,ptmin,ptmax,21,-1.05,1.05); 
      fhBadClusterDeltaIA->SetXTitle("#it{E}_{cluster}");
      fhBadClusterDeltaIA->SetYTitle("#it{A}_{cell in cluster}");
      outputContainer->Add(fhBadClusterDeltaIA); 
    }
  }
  
  if(fStudyWeight)
  {
    fhECellClusterRatio  = new TH2F ("hECellClusterRatio"," cell energy / cluster energy vs cluster energy",
                                     nptbins,ptmin,ptmax, 100,0,1.); 
    fhECellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhECellClusterRatio->SetYTitle("#it{E}_{cell i}/#it{E}_{cluster}");
    outputContainer->Add(fhECellClusterRatio);
    
    fhECellClusterLogRatio  = new TH2F ("hECellClusterLogRatio"," Log(cell energy / cluster energy) vs cluster energy",
                                        nptbins,ptmin,ptmax, 100,-10,0); 
    fhECellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhECellClusterLogRatio->SetYTitle("Log(#it{E}_{cell i}/#it{E}_{cluster})");
    outputContainer->Add(fhECellClusterLogRatio);
    
    fhEMaxCellClusterRatio  = new TH2F ("hEMaxCellClusterRatio"," max cell energy / cluster energy vs cluster energy",
                                        nptbins,ptmin,ptmax, 100,0,1.); 
    fhEMaxCellClusterRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhEMaxCellClusterRatio->SetYTitle("#it{E}_{max cell}/#it{E}_{cluster}");
    outputContainer->Add(fhEMaxCellClusterRatio);
    
    fhEMaxCellClusterLogRatio  = new TH2F ("hEMaxCellClusterLogRatio"," Log(max cell energy / cluster energy) vs cluster energy",
                                           nptbins,ptmin,ptmax, 100,-10,0); 
    fhEMaxCellClusterLogRatio->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhEMaxCellClusterLogRatio->SetYTitle("Log (#it{E}_{max cell}/#it{E}_{cluster})");
    outputContainer->Add(fhEMaxCellClusterLogRatio);
    
    fhECellTotalRatio  = new TH2F ("hECellTotalRatio"," cell energy / sum all energy vs all energy",
                                     nptbins*2,ptmin,ptmax*2, 100,0,1.);
    fhECellTotalRatio->SetXTitle("#it{E}_{total} (GeV) ");
    fhECellTotalRatio->SetYTitle("#it{E}_{cell i}/#it{E}_{total}");
    outputContainer->Add(fhECellTotalRatio);
    
    fhECellTotalLogRatio  = new TH2F ("hECellTotalLogRatio"," Log(cell energy / sum all energy) vs all energy",
                                        nptbins*2,ptmin,ptmax*2, 100,-10,0);
    fhECellTotalLogRatio->SetXTitle("#it{E}_{total} (GeV) ");
    fhECellTotalLogRatio->SetYTitle("Log(#it{E}_{cell i}/#it{E}_{total})");
    outputContainer->Add(fhECellTotalLogRatio);
    
    fhECellTotalRatioMod    = new TH2F*[fNModules];
    fhECellTotalLogRatioMod = new TH2F*[fNModules];
    
    for(Int_t imod = 0; imod < fNModules; imod++)
    {
      fhECellTotalRatioMod[imod]  = new TH2F (Form("hECellTotalRatio_Mod%d",imod),
                                              Form("#cell energy / sum all energy vs all energy in Module %d",imod),
                                              nptbins*2,ptmin,ptmax*2, 100,0,1.);
      fhECellTotalRatioMod[imod]->SetXTitle("#it{E} (GeV)");
      fhECellTotalRatioMod[imod]->SetYTitle("#it{n}_{cells}");
      outputContainer->Add(fhECellTotalRatioMod[imod]);
      
      fhECellTotalLogRatioMod[imod]  = new TH2F (Form("hECellTotalLogRatio_Mod%d",imod),
                                              Form("Log(cell energy / sum all energy) vs all energy in Module %d",imod),
                                              nptbins*2,ptmin,ptmax*2, 100,-10,0);
      fhECellTotalLogRatioMod[imod]->SetXTitle("#it{E} (GeV)");
      fhECellTotalLogRatioMod[imod]->SetYTitle("#it{n}_{cells}");
      outputContainer->Add(fhECellTotalLogRatioMod[imod]);
    }
    
    // To be done properly with setters and data members ...
    Float_t cellEmin [] = {0.05,0.1,0.15,0.2};
    Float_t cellTmin [] = {50.,100.,1000000.};

    for(Int_t iw = 0; iw < 12; iw++)
    {
      Float_t w0 = 4+0.05*iw; // 3+0.25*iw;
      for(Int_t iEmin = 0; iEmin < 4; iEmin++)
      {    
        for(Int_t iTmin = 0; iTmin < 3; iTmin++)
        {
          fhLambda0ForW0AndCellCuts[iw][iEmin][iTmin]  = new TH2F (Form("hLambda0ForW0%d_CellEMin%d_TimeMax%d",iw,iEmin,iTmin),
                                                            Form("#lambda^{2}_{0} vs E, w0=%1.2f, cell E>%2.2f MeV, |t|<%2.0f ns",
                                                                 w0, cellEmin[iEmin], cellTmin[iTmin]),
                                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhLambda0ForW0AndCellCuts[iw][iEmin][iTmin]->SetXTitle("#it{E}_{cluster}");
          fhLambda0ForW0AndCellCuts[iw][iEmin][iTmin]->SetYTitle("#lambda^{2}_{0}");
          outputContainer->Add(fhLambda0ForW0AndCellCuts[iw][iEmin][iTmin]); 

          
//          fhLambda1ForW0AndCellCuts[iw][iEmin][iTmin]  = new TH2F (Form("hLambda1ForW0%d_CellEMin%d_TimeMax%d",iw,iEmin,iTmin),
//                                                            Form("#lambda^{2}_{1} vs E, w0=%1.2f, cell E>%2.2f MeV, |t|<%2.0f ns"",
//                                                                 w0, cellEmin[iEmin], cellTmin[iTmin]),
//                                                            nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
//          fhLambda1ForW0AndCellCuts[iw][iEmin][iTmin]->SetXTitle("#it{E}_{cluster}");
//          fhLambda1ForW0AndCellCuts[iw][iEmin][iTmin]->SetYTitle("#lambda^{2}_{1}");
//          outputContainer->Add(fhLambda1ForW0AndCellCuts[iw][iEmin][iTmin]); 

          fhLambda0ForW0AndCellCutsEta0[iw][iEmin][iTmin]  = new TH2F (Form("hLambda0ForW0%d_CellEMin%d_TimeMax%d_Eta0",iw,iEmin,iTmin),
                                                                Form("#lambda^{2}_{0} vs E, w0=%1.2f, cell E>%2.2f MeV, |t|<%2.0f ns, |#eta| < 0.15",
                                                                     w0, cellEmin[iEmin], cellTmin[iTmin]),
                                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhLambda0ForW0AndCellCutsEta0[iw][iEmin][iTmin]->SetXTitle("#it{E}_{cluster}");
          fhLambda0ForW0AndCellCutsEta0[iw][iEmin][iTmin]->SetYTitle("#lambda^{2}_{0}");
          outputContainer->Add(fhLambda0ForW0AndCellCutsEta0[iw][iEmin][iTmin]); 
        }
      }
      
      if(IsDataMC())
      {
        TString mcnames[] = {"Photon", "Electron","Conversion","Pi0","Hadron"};
        for(Int_t imc = 0; imc < 5; imc++)
        {
          fhLambda0ForW0MC[iw][imc]  = new TH2F (Form("hLambda0ForW0%d_MC%s",iw,mcnames[imc].Data()),
                                                 Form("shower shape, #lambda^{2}_{0} vs E, w0 = %1.1f, for MC %s",w0,mcnames[imc].Data()),
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          fhLambda0ForW0MC[iw][imc]->SetXTitle("#it{E}_{cluster}");
          fhLambda0ForW0MC[iw][imc]->SetYTitle("#lambda^{2}_{0}");
          outputContainer->Add(fhLambda0ForW0MC[iw][imc]); 
          
//          fhLambda1ForW0MC[iw][imc]  = new TH2F (Form("hLambda1ForW0%d_MC%s",iw,mcnames[imc].Data()),
//                                                 Form("shower shape, #lambda^{2}_{1} vs E, w0 = %1.1f, for MC %s",w0,mcnames[imc].Data()),
//                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
//          fhLambda1ForW0MC[iw][imc]->SetXTitle("#it{E}_{cluster}");
//          fhLambda1ForW0MC[iw][imc]->SetYTitle("#lambda^{2}_{1}");
//          outputContainer->Add(fhLambda1ForW0MC[iw][imc]); 
        }
      }
    }     
  }
  
  // Track Matching
  if(fFillAllTMHisto)
  {
    Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
    Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
    Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
    Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
    Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
    Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
    
    fhTrackMatchedDEtaNeg  = new TH2F("hTrackMatchedDEtaNeg","d#eta of cluster-negative track vs cluster energy",
                                   nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaNeg->SetYTitle("d#eta");
    fhTrackMatchedDEtaNeg->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDPhiNeg  = new TH2F("hTrackMatchedDPhiNeg","d#phi of cluster-negative track vs cluster energy",
                                   nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiNeg->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiNeg  = new TH2F("hTrackMatchedDEtaDPhiNeg","d#eta vs d#phi of cluster- negative track vs cluster energy",
                                       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiNeg->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiNeg->SetXTitle("d#eta");
    
    fhTrackMatchedDEtaPos  = new TH2F("hTrackMatchedDEtaPos","d#eta of cluster-positive track vs cluster energy",
                                      nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaPos->SetYTitle("d#eta");
    fhTrackMatchedDEtaPos->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDPhiPos  = new TH2F("hTrackMatchedDPhiPos","d#phi of cluster-positive track vs cluster energy",
                                      nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhiPos->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiPos  = new TH2F("hTrackMatchedDEtaDPhiPos","d#eta vs d#phi of cluster-positive track vs cluster energy",
                                          nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiPos->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhiPos->SetXTitle("d#eta");

    
    fhTrackMatchedDEtaNegMod  = new TH2F("hTrackMatchedDEtaNegPerModule","d#eta of cluster-negative track vs module, E > 0.5 GeV",
                                         nresetabins,resetamin,resetamax,fNModules,0,fNModules);
    fhTrackMatchedDEtaNegMod->SetXTitle("d#eta");
    fhTrackMatchedDEtaNegMod->SetYTitle("Module");
    
    fhTrackMatchedDPhiNegMod  = new TH2F("hTrackMatchedDPhiNegPerModule","d#phi of cluster-negative track vs module, E > 0.5 GeV",
                                         nresetabins,resetamin,resetamax,fNModules,0,fNModules);
    fhTrackMatchedDPhiNegMod->SetXTitle("d#phi (rad)");
    fhTrackMatchedDPhiNegMod->SetYTitle("Module");
    
    fhTrackMatchedDEtaPosMod  = new TH2F("hTrackMatchedDEtaPosPerModule","d#eta of cluster-positive track vs module, E > 0.5 GeV",
                                         nresetabins,resetamin,resetamax,fNModules,0,fNModules);
    fhTrackMatchedDEtaPosMod->SetXTitle("d#eta");
    fhTrackMatchedDEtaPosMod->SetYTitle("Module");
    
    fhTrackMatchedDPhiPosMod  = new TH2F("hTrackMatchedDPhiPosPerModule","d#phi of cluster-positive track vs module, E > 0.5 GeV",
                                         nresetabins,resetamin,resetamax,fNModules,0,fNModules);
    fhTrackMatchedDPhiPosMod->SetXTitle("d#phi (rad)");
    fhTrackMatchedDPhiPosMod->SetYTitle("Module");
  
    outputContainer->Add(fhTrackMatchedDEtaNeg) ;
    outputContainer->Add(fhTrackMatchedDPhiNeg) ;
    outputContainer->Add(fhTrackMatchedDEtaPos) ;
    outputContainer->Add(fhTrackMatchedDPhiPos) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiNeg) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiPos) ;
    
    outputContainer->Add(fhTrackMatchedDEtaNegMod) ;
    outputContainer->Add(fhTrackMatchedDPhiNegMod) ;
    outputContainer->Add(fhTrackMatchedDEtaPosMod) ;
    outputContainer->Add(fhTrackMatchedDPhiPosMod) ;
    
    fhECharged  = new TH1F ("hECharged","#it{E} reconstructed clusters, matched with track", nptbins,ptmin,ptmax);
    fhECharged->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhECharged);
    
    fhPtCharged  = new TH1F ("hPtCharged","#it{p}_{T} reconstructed clusters, matched with track", nptbins,ptmin,ptmax);
    fhPtCharged->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtCharged);
    
    fhPhiCharged  = new TH1F ("hPhiCharged","#phi reconstructed clusters, matched with track",nphibins,phimin,phimax);
    fhPhiCharged->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiCharged);
    
    fhEtaCharged  = new TH1F ("hEtaCharged","#eta reconstructed clusters, matched with track",netabins,etamin,etamax);
    fhEtaCharged->SetXTitle("#eta ");
    outputContainer->Add(fhEtaCharged);
    
    if(fFillAllTH3)
    {
      fhEtaPhiECharged  = new TH3F ("hEtaPhiECharged","#eta vs #phi, reconstructed clusters, matched with track",
                                    netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
      fhEtaPhiECharged->SetXTitle("#eta ");
      fhEtaPhiECharged->SetYTitle("#phi ");
      fhEtaPhiECharged->SetZTitle("#it{E} (GeV) ");
      outputContainer->Add(fhEtaPhiECharged);	
    }
    else
    {
      fhEtaPhiCharged  = new TH2F ("hEtaPhiCharged","#eta vs #phi for #it{E} > 0.5 GeV, reconstructed clusters, with matched track",
                            netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhiCharged->SetXTitle("#eta ");
      fhEtaPhiCharged->SetYTitle("#phi (rad)");
      outputContainer->Add(fhEtaPhiCharged);
    }
    
    fh1EOverP = new TH2F("h1EOverP","TRACK matches #it{E}/#it{p}",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fh1EOverP->SetYTitle("#it{E}/#it{p}");
    fh1EOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fh1EOverP);
    
    fh2dR = new TH2F("h2dR","TRACK matches #Delta #it{R}",nptbins,ptmin,ptmax,ndRbins,dRmin,dRmax);
    fh2dR->SetYTitle("#Delta #it{R} (rad)");
    fh2dR->SetXTitle("#it{E} cluster (GeV)");
    outputContainer->Add(fh2dR) ;
    
    fh2MatchdEdx = new TH2F("h2MatchdEdx","#it{dE/dx} vs. #it{p} for all matches",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2MatchdEdx->SetXTitle("p (GeV/#it{c})");
    fh2MatchdEdx->SetYTitle("<#it{dE/dx}>");
    outputContainer->Add(fh2MatchdEdx);
    
    fh2EledEdx = new TH2F("h2EledEdx","#it{dE/dx} vs. #it{p} for electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fh2EledEdx->SetXTitle("p (GeV/#it{c})");
    fh2EledEdx->SetYTitle("<#it{dE/dx}>");
    outputContainer->Add(fh2EledEdx) ;
    
    fh1EOverPR02 = new TH2F("h1EOverPR02","TRACK matches #it{E}/#it{p}, all",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fh1EOverPR02->SetYTitle("#it{E}/#it{p}");
    fh1EOverPR02->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fh1EOverPR02);	

    fh1EleEOverP = new TH2F("h1EleEOverP","Electron candidates #it{E}/#it{p} (60<#it{dE/dx}<100)",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fh1EleEOverP->SetYTitle("#it{E}/#it{p}");
    fh1EleEOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fh1EleEOverP);
    
    
    // Projections per SM
    
    fh1EOverPMod = new TH2F("h1EOverP_PerModule","TRACK matches #it{E}/#it{p}, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",nPoverEbins,eOverPmin,eOverPmax,fNModules,0,fNModules);
    fh1EOverPMod->SetXTitle("#it{E}/#it{p}");
    fh1EOverPMod->SetYTitle("Module");
    outputContainer->Add(fh1EOverPMod);
    
    fh2dRMod = new TH2F("h2dR_PerModule","TRACK matches #Delta #it{R}, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",ndRbins,dRmin,dRmax,fNModules,0,fNModules);
    fh2dRMod->SetXTitle("#Delta #it{R} (rad)");
    fh2dRMod->SetYTitle("Module");
    outputContainer->Add(fh2dRMod) ;
    
    fh2MatchdEdxMod = new TH2F("h2MatchdEdx_PerModule","#it{dE/dx} vs. #it{p} for all matches, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",ndedxbins,dedxmin,dedxmax,fNModules,0,fNModules);
    fh2MatchdEdxMod->SetYTitle("Module");
    fh2MatchdEdxMod->SetXTitle("<#it{dE/dx}>");
    outputContainer->Add(fh2MatchdEdxMod);
    
    fh2EledEdxMod = new TH2F("h2EledEdx_PerModule","#it{dE/dx} vs. #it{p} for electrons, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",ndedxbins,dedxmin,dedxmax,fNModules,0,fNModules);
    fh2EledEdxMod->SetYTitle("Module");
    fh2EledEdxMod->SetXTitle("<#it{dE/dx}>");
    outputContainer->Add(fh2EledEdxMod) ;
    
    fh1EOverPR02Mod = new TH2F("h1EOverPR02_PerModule","TRACK matches #it{E}/#it{p}, all, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",nPoverEbins,eOverPmin,eOverPmax,fNModules,0,fNModules);
    fh1EOverPR02Mod->SetXTitle("#it{E}/#it{p}");
    fh1EOverPR02Mod->SetYTitle("Module");
    outputContainer->Add(fh1EOverPR02Mod);	
    
    fh1EleEOverPMod = new TH2F("h1EleEOverP_PerModule","Electron candidates #it{E}/#it{p} (60<#it{dE/dx}<100), #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",nPoverEbins,eOverPmin,eOverPmax,fNModules,0,fNModules);
    fh1EleEOverPMod->SetXTitle("#it{E}/#it{p}");
    fh1EleEOverPMod->SetYTitle("Module");
    outputContainer->Add(fh1EleEOverPMod);
  }
  
  if(fFillAllPi0Histo)
  {
    fhIM  = new TH2F ("hIM","Cluster pairs (EMCAL or PHOS) Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1",
                      nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIM->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
    fhIM->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
    outputContainer->Add(fhIM);

    fhIMDiff  = new TH2F ("hIMDiff","Cluster pairs (EMCAL or PHOS) Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, different SM",
                          nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIMDiff->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
    fhIMDiff->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
    outputContainer->Add(fhIMDiff);

    fhIMSame  = new TH2F ("hIMSame","Cluster pairs (EMCAL or PHOS) Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same SM",
                          nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
    fhIMSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
    fhIMSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
    outputContainer->Add(fhIMSame);
    
    if(fNModules > 12 && (GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL))
    {
      fhIMDCAL  = new TH2F ("hIMDCAL","Cluster pairs in DCAL Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1",
                            nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMDCAL->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMDCAL->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMDCAL);

      fhIMDCALDiff  = new TH2F ("hIMDCALDiff","Cluster pairs in DCAL Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, different SM",
                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMDCALDiff->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMDCALDiff->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMDCALDiff);

      fhIMDCALSame  = new TH2F ("hIMDCALSame","Cluster pairs in DCAL Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same SM",
                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMDCALSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMDCALSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMDCALSame);
      
      fhIMDCALPHOS  = new TH2F ("hIMDCALPHOS","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1",
                                nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMDCALPHOS->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMDCALPHOS->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMDCALPHOS);

      fhIMDCALPHOSSame  = new TH2F ("hIMDCALPHOSSame","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same #phi sector",
                                    nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMDCALPHOSSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMDCALPHOSSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMDCALPHOSSame);
    }

    if(fFillInvMassInEMCALWithPHOSDCalAcc)
    {
      fhIMEMCALPHOS  = new TH2F ("hIMEMCALPHOS","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1",
                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMEMCALPHOS->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMEMCALPHOS->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMEMCALPHOS);
      
      fhIMEMCALPHOSSame  = new TH2F ("hIMEMCALPHOSSame","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same #phi sector",
                                     nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMEMCALPHOSSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMEMCALPHOSSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMEMCALPHOSSame);    
    }
    
    if(fFillPi0PairDiffTime)
    {
      fhClusterPairDiffTimeEPi0Mass = new TH2F("hClusterPairDiffTimeEPi0Mass","cluster pair time difference vs E, only good clusters",
                                               nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
      fhClusterPairDiffTimeEPi0Mass->SetXTitle("#it{E}_{cluster} (GeV)");
      fhClusterPairDiffTimeEPi0Mass->SetYTitle("#Delta #it{t} (ns)");
      outputContainer->Add(fhClusterPairDiffTimeEPi0Mass);  
      
      fhClusterPairDiffTimeEPi0MassSame = new TH2F("hClusterPairDiffTimeEPi0MassSame","cluster pair time difference vs E, only good clusters",
                                                   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
      fhClusterPairDiffTimeEPi0MassSame->SetXTitle("#it{E}_{cluster} (GeV)");
      fhClusterPairDiffTimeEPi0MassSame->SetYTitle("#Delta #it{t} (ns)");
      outputContainer->Add(fhClusterPairDiffTimeEPi0MassSame);  
      
      if(fNModules > 12 && (GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL))
      {
        fhClusterPairDiffTimeEPi0MassDCal = new TH2F("hClusterPairDiffTimeEPi0MassDCal","cluster pair time difference vs E, only good clusters",
                                                     nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
        fhClusterPairDiffTimeEPi0MassDCal->SetXTitle("#it{E}_{cluster} (GeV)");
        fhClusterPairDiffTimeEPi0MassDCal->SetYTitle("#Delta #it{t} (ns)");
        outputContainer->Add(fhClusterPairDiffTimeEPi0MassDCal);  
        
        fhClusterPairDiffTimeEPi0MassDCalSame = new TH2F("hClusterPairDiffTimeEPi0MassDCalSame","cluster pair time difference vs E, only good clusters",
                                                         nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
        fhClusterPairDiffTimeEPi0MassDCalSame->SetXTitle("#it{E}_{cluster} (GeV)");
        fhClusterPairDiffTimeEPi0MassDCalSame->SetYTitle("#Delta #it{t} (ns)");
        outputContainer->Add(fhClusterPairDiffTimeEPi0MassDCalSame);  
      }
    }
    
    fhAsym  = new TH2F ("hAssym","Cluster pairs Asymmetry vs reconstructed pair energy",nptbins,ptmin,ptmax,nasymbins,asymmin,asymmax); 
    fhAsym->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
    fhAsym->SetYTitle("#it{Asymmetry}");
    outputContainer->Add(fhAsym);	
    
    if(fFillInvMassOpenAngle)
    {
      fhOpAngle  = new TH2F ("hOpeningAngle","Cluster pairs opening angle vs reconstructed pair #it{p}_{T}, ncell > 1",
                             nptbins,ptmin,ptmax, 180,0,TMath::Pi()); 
      fhOpAngle->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhOpAngle->SetYTitle("Opening angle (degrees)");
      outputContainer->Add(fhOpAngle);
      
      Int_t nBinOpAngle = TMath::Nint(fInvMassMaxOpenAngle*TMath::RadToDeg()); // bin of 1 degree size
      fhIMvsOpAngle  = new TH2F ("hIMvsOpAngle","Cluster pairs Invariant mass vs reconstructed pair opening angle, ncell > 1",
                                 nBinOpAngle,0.,fInvMassMaxOpenAngle,nmassbins,massmin,massmax); 
      fhIMvsOpAngle->SetXTitle("Opening angle (degrees)");
      fhIMvsOpAngle->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMvsOpAngle);  
    }
  }
  
  if(fFillAllPosHisto2)
  {
    fhXYZ  = new TH3F ("hXYZ","Cluster: #it{x} vs #it{y} vs #it{z}",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax);
    fhXYZ->SetXTitle("#it{x} (cm)");
    fhXYZ->SetYTitle("#it{y} (cm)");
    fhXYZ->SetZTitle("#it{z} (cm) ");
    outputContainer->Add(fhXYZ);  
    
    fhXE  = new TH2F ("hXE","Cluster X position vs cluster energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXE->SetXTitle("#it{x} (cm)");
    fhXE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhXE);
    
    fhYE  = new TH2F ("hYE","Cluster Y position vs cluster energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYE->SetXTitle("#it{y} (cm)");
    fhYE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhYE);
    
    fhZE  = new TH2F ("hZE","Cluster Z position vs cluster energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZE->SetXTitle("#it{z} (cm)");
    fhZE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhZE);    
    
    fhRE  = new TH2F ("hRE","Cluster R position vs cluster energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhRE);
    
    fhXNCells  = new TH2F ("hXNCells","Cluster X position vs N Cells per Cluster",xbins,xmin,xmax,nceclbins,nceclmin,nceclmax); 
    fhXNCells->SetXTitle("#it{x} (cm)");
    fhXNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhXNCells);
    
    fhYNCells  = new TH2F ("hYNCells","Cluster Y position vs N Cells per Cluster",ybins,ymin,ymax,nceclbins,nceclmin,nceclmax); 
    fhYNCells->SetXTitle("#it{y} (cm)");
    fhYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhYNCells);

    fhZNCells  = new TH2F ("hZNCells","Cluster Z position vs N Cells per Cluster",zbins,zmin,zmax,nceclbins,nceclmin,nceclmax); 
    fhZNCells->SetXTitle("#it{z} (cm)");
    fhZNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhZNCells);
        
    fhRNCells  = new TH2F ("hRNCells","Cluster R position vs N Cells per Cluster",rbins,rmin,rmax,nceclbins,nceclmin,nceclmax); 
    fhRNCells->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhRNCells);
  }
  
  if(fFillAllPosHisto)
  {
    fhRCellE  = new TH2F ("hRCellE","Cell R position vs cell energy",rbins,rmin,rmax,nptbins,ptmin,ptmax); 
    fhRCellE->SetXTitle("r = #sqrt{x^{2}+y^{2}} (cm)");
    fhRCellE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhRCellE);
    
    fhXCellE  = new TH2F ("hXCellE","Cell X position vs cell energy",xbins,xmin,xmax,nptbins,ptmin,ptmax); 
    fhXCellE->SetXTitle("#it{x} (cm)");
    fhXCellE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhXCellE);
    
    fhYCellE  = new TH2F ("hYCellE","Cell Y position vs cell energy",ybins,ymin,ymax,nptbins,ptmin,ptmax); 
    fhYCellE->SetXTitle("#it{y} (cm)");
    fhYCellE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhYCellE);
    
    fhZCellE  = new TH2F ("hZCellE","Cell Z position vs cell energy",zbins,zmin,zmax,nptbins,ptmin,ptmax); 
    fhZCellE->SetXTitle("#it{z} (cm)");
    fhZCellE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhZCellE);
    
    fhXYZCell  = new TH3F ("hXYZCell","Cell : #it{x} vs #it{y} vs #it{z}",xbins,xmin,xmax,ybins,ymin,ymax,zbins,zmin,zmax);
    fhXYZCell->SetXTitle("#it{x} (cm)");
    fhXYZCell->SetYTitle("#it{y} (cm)");
    fhXYZCell->SetZTitle("#it{z} (cm)");
    outputContainer->Add(fhXYZCell);
    
    Float_t dx = TMath::Abs(xmin)+TMath::Abs(xmax);
    Float_t dy = TMath::Abs(ymin)+TMath::Abs(ymax);
    Float_t dz = TMath::Abs(zmin)+TMath::Abs(zmax);
    Float_t dr = TMath::Abs(rmin)+TMath::Abs(rmax);
    
    fhDeltaCellClusterRNCells  = new TH2F ("hDeltaCellClusterRNCells","Cluster-Cell R position vs N Cells per Cluster",rbins*2,-dr,dr,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterRNCells->SetXTitle("#it{r} = #sqrt{x^{2}+y^{2}}, #it{r}_{clus}-#it{r}_{cell}  (cm)");
    fhDeltaCellClusterRNCells->SetYTitle("#it{n}_{cells per cluster}");
    outputContainer->Add(fhDeltaCellClusterRNCells);
    
    fhDeltaCellClusterXNCells  = new TH2F ("hDeltaCellClusterXNCells","Cluster-Cell X position vs N Cells per Cluster",xbins*2,-dx,dx,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterXNCells->SetXTitle("#it{x}_{clus}-#it{x}_{cell} (cm)");
    fhDeltaCellClusterXNCells->SetYTitle("#it{n}_{cells per cluster}");
    outputContainer->Add(fhDeltaCellClusterXNCells);
    
    fhDeltaCellClusterYNCells  = new TH2F ("hDeltaCellClusterYNCells","Cluster-Cell Y position vs N Cells per Cluster",ybins*2,-dy,dy,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterYNCells->SetXTitle("#it{y}_{clus}-#it{y}_{cell} (cm)");
    fhDeltaCellClusterYNCells->SetYTitle("N cells per cluster");
    outputContainer->Add(fhDeltaCellClusterYNCells);
    
    fhDeltaCellClusterZNCells  = new TH2F ("hDeltaCellClusterZNCells","Cluster-Cell Z position vs N Cells per Cluster",zbins*2,-dz,dz,nceclbins,nceclmin,nceclmax); 
    fhDeltaCellClusterZNCells->SetXTitle("#it{z}_{clus}-#it{z}_{cell} (cm)");
    fhDeltaCellClusterZNCells->SetYTitle("#it{n}_{cells per cluster}");
    outputContainer->Add(fhDeltaCellClusterZNCells);
    
    fhDeltaCellClusterRE  = new TH2F ("hDeltaCellClusterRE","Cluster-Cell R position vs cluster energy",rbins*2,-dr,dr,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterRE->SetXTitle("#it{r} = #sqrt{x^{2}+y^{2}}, #it{r}_{clus}-#it{r}_{cell} (cm)");
    fhDeltaCellClusterRE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhDeltaCellClusterRE);		
    
    fhDeltaCellClusterXE  = new TH2F ("hDeltaCellClusterXE","Cluster-Cell X position vs cluster energy",xbins*2,-dx,dx,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterXE->SetXTitle("#it{x}_{clus}-#it{x}_{cell} (cm)");
    fhDeltaCellClusterXE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhDeltaCellClusterXE);
    
    fhDeltaCellClusterYE  = new TH2F ("hDeltaCellClusterYE","Cluster-Cell Y position vs cluster energy",ybins*2,-dy,dy,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterYE->SetXTitle("#it{y}_{clus}-#it{y}_{cell} (cm)");
    fhDeltaCellClusterYE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhDeltaCellClusterYE);
    
    fhDeltaCellClusterZE  = new TH2F ("hDeltaCellClusterZE","Cluster-Cell Z position vs cluster energy",zbins*2,-dz,dz,nptbins,ptmin,ptmax); 
    fhDeltaCellClusterZE->SetXTitle("#it{z}_{clus}-#it{z}_{cell} (cm)");
    fhDeltaCellClusterZE->SetYTitle("#it{E} (GeV)");
    outputContainer->Add(fhDeltaCellClusterZE);
    
    if(fFillAllTH3)
    {
      fhEtaPhiAmpCell  = new TH3F ("hEtaPhiAmpCell","Cell #eta vs cell #phi vs cell energy",
                                   netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
      fhEtaPhiAmpCell->SetXTitle("#eta ");
      fhEtaPhiAmpCell->SetYTitle("#phi (rad)");
      fhEtaPhiAmpCell->SetZTitle("#it{E} (GeV) ");
      outputContainer->Add(fhEtaPhiAmpCell);
    }
    else
    {
      fhEtaPhiCell  = new TH2F ("hEtaPhiCell","Cell #eta vs cell #phi vs cell energy",
                                netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhiCell->SetXTitle("#eta ");
      fhEtaPhiCell->SetYTitle("#phi (rad)");
      outputContainer->Add(fhEtaPhiCell);
    }
  }
  
  // Calorimeter cells
    
  fhNCells  = new TH1F ("hNCells","# cells", ncebins,ncemin+0.5,ncemax); 
  fhNCells->SetXTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCells);

  fhNCellsCutAmpMin  = new TH1F ("hNCellsCutAmpMin",Form("# cells amp > %1.2f-%1.2f",fEMCALCellAmpMin,fPHOSCellAmpMin), ncebins,ncemin+0.5,ncemax);
  fhNCellsCutAmpMin->SetXTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsCutAmpMin);
  
  fhAmplitude  = new TH1F ("hAmplitude","#it{E}_{cell}", nptbins*2,ptmin,ptmax);
  fhAmplitude->SetXTitle("#it{E}_{cell} (GeV)");
  outputContainer->Add(fhAmplitude);
  
  fhAmpId  = new TH2F ("hAmpId","#it{E}_{cell}", nfineptbins,ptfinemin,ptfinemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules); 
  fhAmpId->SetXTitle("#it{E}_{cell} (GeV)");
  outputContainer->Add(fhAmpId);

  fhAmpIdLowGain  = new TH2F ("hAmpIdLG","Low gain: #it{E}_{cell}", nfineptbins,15,ptfinemax+15,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules);
  fhAmpIdLowGain->SetXTitle("#it{E}_{cell} (GeV)");
  outputContainer->Add(fhAmpIdLowGain);
  
  if(fFillAllCellTimeHisto)
  {
    fhCellTimeSpreadRespectToCellMax = new TH2F ("hCellTimeSpreadRespectToCellMax","t_{cell max}-t_{cell i} per cluster", nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
    fhCellTimeSpreadRespectToCellMax->SetXTitle("#it{E} (GeV)");
    fhCellTimeSpreadRespectToCellMax->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
    outputContainer->Add(fhCellTimeSpreadRespectToCellMax);
    
//    fhClusterMaxCellDiffAverageTime = new TH2F ("hClusterMaxCellDiffAverageTime","t_{cell max}-t_{average} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
//    fhClusterMaxCellDiffAverageTime->SetXTitle("#it{E} (GeV)");
//    fhClusterMaxCellDiffAverageTime->SetYTitle("#Delta #it{t}_{cell max - average} (ns)");
//    outputContainer->Add(fhClusterMaxCellDiffAverageTime);
//        
//    fhClusterMaxCellDiffWeightedTime = new TH2F ("hClusterMaxCellDiffWeightedTime","t_{cell max}-t_{weighted} per cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax); 
//    fhClusterMaxCellDiffWeightedTime->SetXTitle("#it{E} (GeV)");
//    fhClusterMaxCellDiffWeightedTime->SetYTitle("#Delta #it{t}_{cell max - weighted} (ns)");
//    outputContainer->Add(fhClusterMaxCellDiffWeightedTime);
    
    fhCellIdCellLargeTimeSpread= new TH1F ("hCellIdCellLargeTimeSpread","Cells with time 100 ns larger than cell max in cluster ", 
                                           fNMaxCols*fNMaxRows*fNModules,0,fNMaxCols*fNMaxRows*fNModules); 
    fhCellIdCellLargeTimeSpread->SetXTitle("Absolute Cell Id");
    outputContainer->Add(fhCellIdCellLargeTimeSpread);
    
    fhTime  = new TH1F ("hTime","#it{t}_{cell}",ntimebins,timemin,timemax); 
    fhTime->SetXTitle("#it{t}_{cell} (ns)");
    outputContainer->Add(fhTime);
    
//    fhTimeVz  = new TH2F ("hTimeVz","#it{t}_{cell} vs vertex, amplitude > 0.5 GeV",100, 0, 50,ntimebins,timemin,timemax); 
//    fhTimeVz->SetXTitle("|v_{z}| (cm)");
//    fhTimeVz->SetYTitle("#it{t}_{cell} (ns)");
//    outputContainer->Add(fhTimeVz);
    
    fhTimeId  = new TH2F ("hTimeId","#it{t}_{cell} vs Absolute Id",
                          ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules); 
    fhTimeId->SetXTitle("#it{t}_{cell} (ns)");
    fhTimeId->SetYTitle("Cell Absolute Id");
    outputContainer->Add(fhTimeId);
    
    fhTimeAmp  = new TH2F ("hTimeAmp","#it{t}_{cell} vs #it{E}_{cell}",nptbins*2,ptmin,ptmax,ntimebins,timemin,timemax); 
    fhTimeAmp->SetYTitle("#it{t}_{cell} (ns)");
    fhTimeAmp->SetXTitle("#it{E}_{cell} (GeV)");
    outputContainer->Add(fhTimeAmp);
    
    fhTimeIdLowGain  = new TH2F ("hTimeIdLG","Low gain: #it{t}_{cell} vs Absolute Id",
                          ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules);
    fhTimeIdLowGain->SetXTitle("#it{t}_{cell} (ns)");
    fhTimeIdLowGain->SetYTitle("Cell Absolute Id");
    outputContainer->Add(fhTimeIdLowGain);
    
    fhTimeAmpLowGain  = new TH2F ("hTimeAmpLG","Low gain: #it{t}_{cell} vs #it{E}_{cell}",nptbins*2,ptmin,ptmax,ntimebins,timemin,timemax);
    fhTimeAmpLowGain->SetYTitle("#it{t}_{cell} (ns)");
    fhTimeAmpLowGain->SetXTitle("#it{E}_{cell} (GeV)");
    outputContainer->Add(fhTimeAmpLowGain);
  }
  
  fhCellECross  = new TH2F ("hCellECross","1 - Energy in cross around cell /  cell energy",
                            nptbins,ptmin,ptmax, 400,-1,1.); 
  fhCellECross->SetXTitle("#it{E}_{cell} (GeV) ");
  fhCellECross->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell}");
  outputContainer->Add(fhCellECross);    
  
  
  if(fCorrelate)
  {
    // PHOS vs EMCAL
    fhEMCALPHOSCorrNClusters  = new TH2F ("hEMCALPHOSCorrNClusters","# clusters in EMCAL vs PHOS", nclbins,nclmin,nclmax,nclbins,nclmin,nclmax); 
    fhEMCALPHOSCorrNClusters->SetXTitle("number of clusters in EMCAL");
    fhEMCALPHOSCorrNClusters->SetYTitle("number of clusters in PHOS");
    outputContainer->Add(fhEMCALPHOSCorrNClusters);
    
    fhEMCALPHOSCorrEClusters  = new TH2F ("hEMCALPHOSCorrEClusters","summed energy of clusters in EMCAL vs PHOS", nptbins,ptmin,ptmax*2,nptbins,ptmin,ptmax*2);
    fhEMCALPHOSCorrEClusters->SetXTitle("#Sigma #it{E} of clusters in EMCAL (GeV)");
    fhEMCALPHOSCorrEClusters->SetYTitle("#Sigma #it{E} of clusters in PHOS (GeV)");
    outputContainer->Add(fhEMCALPHOSCorrEClusters);
    
    fhEMCALPHOSCorrNCells  = new TH2F ("hEMCALPHOSCorrNCells","# Cells in EMCAL vs PHOS", ncebins,ncemin,ncemax, ncebins,ncemin,ncemax); 
    fhEMCALPHOSCorrNCells->SetXTitle("number of Cells in EMCAL");
    fhEMCALPHOSCorrNCells->SetYTitle("number of Cells in PHOS");
    outputContainer->Add(fhEMCALPHOSCorrNCells);
    
    fhEMCALPHOSCorrECells  = new TH2F ("hEMCALPHOSCorrECells","summed energy of Cells in EMCAL vs PHOS", nptbins*2,ptmin,ptmax*4,nptbins*2,ptmin,ptmax*4);
    fhEMCALPHOSCorrECells->SetXTitle("#Sigma #it{E} of Cells in EMCAL (GeV)");
    fhEMCALPHOSCorrECells->SetYTitle("#Sigma #it{E} of Cells in PHOS (GeV)");
    outputContainer->Add(fhEMCALPHOSCorrECells);

    // DCal vs EMCAL
    fhEMCALDCALCorrNClusters  = new TH2F ("hEMCALDCALCorrNClusters","# clusters in EMCAL vs DCAL", nclbins,nclmin,nclmax,nclbins,nclmin,nclmax); 
    fhEMCALDCALCorrNClusters->SetXTitle("number of clusters in EMCAL");
    fhEMCALDCALCorrNClusters->SetYTitle("number of clusters in DCAL");
    outputContainer->Add(fhEMCALDCALCorrNClusters);
    
    fhEMCALDCALCorrEClusters  = new TH2F ("hEMCALDCALCorrEClusters","summed energy of clusters in EMCAL vs DCAL", nptbins,ptmin,ptmax*2,nptbins,ptmin,ptmax*2);
    fhEMCALDCALCorrEClusters->SetXTitle("#Sigma #it{E} of clusters in EMCAL (GeV)");
    fhEMCALDCALCorrEClusters->SetYTitle("#Sigma #it{E} of clusters in DCAL (GeV)");
    outputContainer->Add(fhEMCALDCALCorrEClusters);
    
    fhEMCALDCALCorrNCells  = new TH2F ("hEMCALDCALCorrNCells","# Cells in EMCAL vs DCAL", ncebins,ncemin,ncemax, ncebins,ncemin,ncemax); 
    fhEMCALDCALCorrNCells->SetXTitle("number of Cells in EMCAL");
    fhEMCALDCALCorrNCells->SetYTitle("number of Cells in DCAL");
    outputContainer->Add(fhEMCALDCALCorrNCells);
    
    fhEMCALDCALCorrECells  = new TH2F ("hEMCALDCALCorrECells","summed energy of Cells in EMCAL vs DCAL", nptbins*2,ptmin,ptmax*4,nptbins*2,ptmin,ptmax*4);
    fhEMCALDCALCorrECells->SetXTitle("#Sigma #it{E} of Cells in EMCAL (GeV)");
    fhEMCALDCALCorrECells->SetYTitle("#Sigma #it{E} of Cells in DCAL (GeV)");
    outputContainer->Add(fhEMCALDCALCorrECells);

    
    // DCAL vs PHOS
    fhDCALPHOSCorrNClusters  = new TH2F ("hDCALPHOSCorrNClusters","# clusters in DCAL vs PHOS", nclbins,nclmin,nclmax,nclbins,nclmin,nclmax); 
    fhDCALPHOSCorrNClusters->SetXTitle("number of clusters in DCAL");
    fhDCALPHOSCorrNClusters->SetYTitle("number of clusters in PHOS");
    outputContainer->Add(fhDCALPHOSCorrNClusters);
    
    fhDCALPHOSCorrEClusters  = new TH2F ("hDCALPHOSCorrEClusters","summed energy of clusters in DCAL vs PHOS", nptbins,ptmin,ptmax*2,nptbins,ptmin,ptmax*2);
    fhDCALPHOSCorrEClusters->SetXTitle("#Sigma #it{E} of clusters in DCAL (GeV)");
    fhDCALPHOSCorrEClusters->SetYTitle("#Sigma #it{E} of clusters in PHOS (GeV)");
    outputContainer->Add(fhDCALPHOSCorrEClusters);
    
    fhDCALPHOSCorrNCells  = new TH2F ("hDCALPHOSCorrNCells","# Cells in DCAL vs PHOS", ncebins,ncemin,ncemax, ncebins,ncemin,ncemax); 
    fhDCALPHOSCorrNCells->SetXTitle("number of Cells in DCAL");
    fhDCALPHOSCorrNCells->SetYTitle("number of Cells in PHOS");
    outputContainer->Add(fhDCALPHOSCorrNCells);
    
    fhDCALPHOSCorrECells  = new TH2F ("hDCALPHOSCorrECells","summed energy of Cells in DCAL vs PHOS", nptbins*2,ptmin,ptmax*4,nptbins*2,ptmin,ptmax*4);
    fhDCALPHOSCorrECells->SetXTitle("#Sigma #it{E} of Cells in DCAL (GeV)");
    fhDCALPHOSCorrECells->SetYTitle("#Sigma #it{E} of Cells in PHOS (GeV)");
    outputContainer->Add(fhDCALPHOSCorrECells);
    
    // Calorimeter vs. V0 signal
      
    fhCaloV0SCorrNClusters  = new TH2F ("hCaloV0SNClusters",Form("# clusters in %s vs V0 signal",GetCalorimeterString().Data()), nv0sbins,nv0smin,nv0smax,nclbins,nclmin,nclmax); 
    fhCaloV0SCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrNClusters->SetYTitle(Form("number of clusters in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0SCorrNClusters);
    
    fhCaloV0SCorrEClusters  = new TH2F ("hCaloV0SEClusters",Form("summed energy of clusters in %s vs V0 signal",GetCalorimeterString().Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax*2);
    fhCaloV0SCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0SCorrEClusters->SetYTitle(Form("#Sigma #it{E} of clusters in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0SCorrEClusters);
    
    fhCaloV0SCorrNCells  = new TH2F ("hCaloV0SNCells",Form("# Cells in %s vs V0 signal",GetCalorimeterString().Data()), nv0sbins,nv0smin,nv0smax, ncebins,ncemin,ncemax); 
    fhCaloV0SCorrNCells->SetXTitle("V0 signal");
    fhCaloV0SCorrNCells->SetYTitle(Form("number of Cells in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0SCorrNCells);
    
    fhCaloV0SCorrECells  = new TH2F ("hCaloV0SECells",Form("summed energy of Cells in %s vs V0 signal",GetCalorimeterString().Data()), nv0sbins,nv0smin,nv0smax,nptbins,ptmin,ptmax*2);
    fhCaloV0SCorrECells->SetXTitle("V0 signal");
    fhCaloV0SCorrECells->SetYTitle(Form("#Sigma #it{E} of Cells in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0SCorrECells);    
    
    // Calorimeter vs V0 multiplicity
      
    fhCaloV0MCorrNClusters  = new TH2F ("hCaloV0MNClusters",Form("# clusters in %s vs V0 signal",GetCalorimeterString().Data()), nv0mbins,nv0mmin,nv0mmax,nclbins,nclmin,nclmax); 
    fhCaloV0MCorrNClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrNClusters->SetYTitle(Form("number of clusters in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0MCorrNClusters);
    
    fhCaloV0MCorrEClusters  = new TH2F ("hCaloV0MEClusters",Form("summed energy of clusters in %s vs V0 signal",GetCalorimeterString().Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax*2);
    fhCaloV0MCorrEClusters->SetXTitle("V0 signal");
    fhCaloV0MCorrEClusters->SetYTitle(Form("#Sigma #it{E} of clusters in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0MCorrEClusters);
    
    fhCaloV0MCorrNCells  = new TH2F ("hCaloV0MNCells",Form("# Cells in %s vs V0 signal",GetCalorimeterString().Data()), nv0mbins,nv0mmin,nv0mmax, ncebins,ncemin,ncemax); 
    fhCaloV0MCorrNCells->SetXTitle("V0 signal");
    fhCaloV0MCorrNCells->SetYTitle(Form("number of Cells in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0MCorrNCells);
    
    fhCaloV0MCorrECells  = new TH2F ("hCaloV0MECells",Form("summed energy of Cells in %s vs V0 signal",GetCalorimeterString().Data()), nv0mbins,nv0mmin,nv0mmax,nptbins,ptmin,ptmax*2);
    fhCaloV0MCorrECells->SetXTitle("V0 signal");
    fhCaloV0MCorrECells->SetYTitle(Form("#Sigma #it{E} of Cells in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloV0MCorrECells);    
    
    //Calorimeter VS Track multiplicity
    fhCaloTrackMCorrNClusters  = new TH2F ("hCaloTrackMNClusters",Form("# clusters in %s vs # tracks",GetCalorimeterString().Data()), ntrmbins,ntrmmin,ntrmmax,nclbins,nclmin,nclmax); 
    fhCaloTrackMCorrNClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrNClusters->SetYTitle(Form("number of clusters in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloTrackMCorrNClusters);
    
    fhCaloTrackMCorrEClusters  = new TH2F ("hCaloTrackMEClusters",Form("summed energy of clusters in %s vs # tracks",GetCalorimeterString().Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax*2);
    fhCaloTrackMCorrEClusters->SetXTitle("# tracks");
    fhCaloTrackMCorrEClusters->SetYTitle(Form("#Sigma #it{E} of clusters in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloTrackMCorrEClusters);
    
    fhCaloTrackMCorrNCells  = new TH2F ("hCaloTrackMNCells",Form("# Cells in %s vs # tracks",GetCalorimeterString().Data()), ntrmbins,ntrmmin,ntrmmax, ncebins,ncemin,ncemax); 
    fhCaloTrackMCorrNCells->SetXTitle("# tracks");
    fhCaloTrackMCorrNCells->SetYTitle(Form("number of Cells in %s",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloTrackMCorrNCells);
    
    fhCaloTrackMCorrECells  = new TH2F ("hCaloTrackMECells",Form("summed energy of Cells in %s vs # tracks",GetCalorimeterString().Data()), ntrmbins,ntrmmin,ntrmmax,nptbins,ptmin,ptmax*2);
    fhCaloTrackMCorrECells->SetXTitle("# tracks");
    fhCaloTrackMCorrECells->SetYTitle(Form("#Sigma #it{E} of Cells in %s (GeV)",GetCalorimeterString().Data()));
    outputContainer->Add(fhCaloTrackMCorrECells);    
    
    fhCaloCenNClusters  = new TH2F ("hCaloCenNClusters","# clusters in calorimeter vs centrality",100,0,100,nclbins,nclmin,nclmax);
    fhCaloCenNClusters->SetYTitle("number of clusters in calorimeter");
    fhCaloCenNClusters->SetXTitle("Centrality");
    outputContainer->Add(fhCaloCenNClusters);
    
    fhCaloCenEClusters  = new TH2F ("hCaloCenEClusters","summed energy of clusters in calorimeter vs centrality",100,0,100,nptbins,ptmin,ptmax*2);
    fhCaloCenEClusters->SetYTitle("#Sigma #it{E} of clusters in calorimeter (GeV)");
    fhCaloCenEClusters->SetXTitle("Centrality");
    outputContainer->Add(fhCaloCenEClusters);
    
    fhCaloCenNCells  = new TH2F ("hCaloCenNCells","# Cells in calorimeter vs centrality",100,0,100,ncebins,ncemin,ncemax);
    fhCaloCenNCells->SetYTitle("number of Cells in calorimeter");
    fhCaloCenNCells->SetXTitle("Centrality");
    outputContainer->Add(fhCaloCenNCells);
    
    fhCaloCenECells  = new TH2F ("hCaloCenECells","summed energy of Cells in calorimeter vs centrality",100,0,100,nptbins*2,ptmin,ptmax*4);
    fhCaloCenECells->SetYTitle("#Sigma #it{E} of Cells in calorimeter (GeV)");
    fhCaloCenECells->SetXTitle("Centrality");
    outputContainer->Add(fhCaloCenECells);

    fhCaloEvPNClusters  = new TH2F ("hCaloEvPNClusters","# clusters in calorimeter vs event plane angle",100,0,TMath::Pi(),nclbins,nclmin,nclmax);
    fhCaloEvPNClusters->SetYTitle("number of clusters in calorimeter");
    fhCaloEvPNClusters->SetXTitle("Event plane angle (rad)");
    outputContainer->Add(fhCaloEvPNClusters);
    
    fhCaloEvPEClusters  = new TH2F ("hCaloEvPEClusters","summed energy of clusters in calorimeter vs  event plane angle",100,0,TMath::Pi(),nptbins,ptmin,ptmax*2);
    fhCaloEvPEClusters->SetYTitle("#Sigma #it{E} of clusters in calorimeter (GeV)");
    fhCaloEvPEClusters->SetXTitle("Event plane angle (rad)");
    outputContainer->Add(fhCaloEvPEClusters);
    
    fhCaloEvPNCells  = new TH2F ("hCaloEvPNCells","# Cells in calorimeter vs  event plane angle",100,0,TMath::Pi(),ncebins,ncemin,ncemax);
    fhCaloEvPNCells->SetYTitle("number of Cells in calorimeter");
    fhCaloEvPNCells->SetXTitle("Event plane angle (rad)");
    outputContainer->Add(fhCaloEvPNCells);
    
    fhCaloEvPECells  = new TH2F ("hCaloEvPECells","summed energy of Cells in calorimeter vs  event plane angle",100,0,TMath::Pi(),nptbins*2,ptmin,ptmax*4);
    fhCaloEvPECells->SetYTitle("#Sigma #it{E} of Cells in calorimeter (GeV)");
    fhCaloEvPECells->SetXTitle("Event plane angle (rad)");
    outputContainer->Add(fhCaloEvPECells);
  } // correlate calorimeters
  
  // Module histograms
  
  fhEMod  = new TH2F ("hE_Mod","Cluster reconstructed Energy in each present Module",nptbins,ptmin,ptmax,fNModules,0,fNModules); 
  fhEMod->SetXTitle("#it{E} (GeV)");
  fhEMod->SetYTitle("Module");
  outputContainer->Add(fhEMod);
  
  fhAmpMod  = new TH2F ("hAmp_Mod","Cell energy in each present Module",nptbins,ptmin,ptmax,fNModules,0,fNModules); 
  fhAmpMod->SetXTitle("#it{E} (GeV)");
  fhAmpMod->SetYTitle("Module");
  outputContainer->Add(fhAmpMod);
  
  fhEWeirdMod  = new TH2F ("hEWeird_Mod","Cluster reconstructed Energy in each present Module, ridiculously large E",200,0,10000,fNModules,0,fNModules); 
  fhEWeirdMod->SetXTitle("#it{E} (GeV)");
  fhEWeirdMod->SetYTitle("Module");
  outputContainer->Add(fhEWeirdMod);
  
  fhAmpWeirdMod  = new TH2F ("hAmpWeird_Mod","Cell energy in each present Module, ridiculously large E",200,0,10000,fNModules,0,fNModules); 
  fhAmpWeirdMod->SetXTitle("#it{E} (GeV)");
  fhAmpWeirdMod->SetYTitle("Module");
  outputContainer->Add(fhAmpWeirdMod);
  
  if(fFillAllCellTimeHisto) 
  {
    fhTimeMod  = new TH2F ("hTime_Mod","Cell time in each present Module",ntimebins,timemin,timemax,fNModules,0,fNModules); 
    fhTimeMod->SetXTitle("t (ns)");
    fhTimeMod->SetYTitle("Module");
    outputContainer->Add(fhTimeMod);
  }
  
  fhNClustersMod  = new TH2F ("hNClusters_Mod","# clusters vs Module", nclbins,nclmin+0.5,nclmax,fNModules,0,fNModules); 
  fhNClustersMod->SetXTitle("number of clusters");
  fhNClustersMod->SetYTitle("Module");
  outputContainer->Add(fhNClustersMod);
  
  fhNCellsMod  = new TH2F ("hNCells_Mod","# cells vs Module", ncebins,ncemin+0.5,ncemax,fNModules,0,fNModules); 
  fhNCellsMod->SetXTitle("#it{n}_{cells}");
  fhNCellsMod->SetYTitle("Module");
  outputContainer->Add(fhNCellsMod);

  fhSumClustersEnergyMod  = new TH2F ("hSumClustersEnergy_Mod","# clusters vs Module", 1000, 0, 2000,fNModules,0,fNModules); 
  fhSumClustersEnergyMod->SetXTitle("#Sigma_{clusters} #it{E} (GeV)");
  fhSumClustersEnergyMod->SetYTitle("Module");
  outputContainer->Add(fhSumClustersEnergyMod);
  
  fhSumCellsAmpMod  = new TH2F ("hSumCellsAmp_Mod","# cells vs Module", 1000, 0, 2000,fNModules,0,fNModules); 
  fhSumCellsAmpMod->SetXTitle("#Sigma_{cells} #it{Amp} (GeV)");
  fhSumCellsAmpMod->SetYTitle("Module");
  outputContainer->Add(fhSumCellsAmpMod);
  
  Int_t colmaxs = fNMaxCols;
  Int_t rowmaxs = fNMaxRows;
  if(GetCalorimeter()==kEMCAL)
  {
    colmaxs=2*fNMaxCols;
    rowmaxs=Int_t(fNModules/2)*fNMaxRows;
  }
  else
  {
    rowmaxs=fNModules*fNMaxRows;
  }
  
  fhGridCells  = new TH2F ("hGridCells",Form("Entries in grid of cells"), 
                           colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5); 
  fhGridCells->SetYTitle("row (phi direction)");
  fhGridCells->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCells);
  
  fhGridCellsE  = new TH2F ("hGridCellsE","Accumulated energy in grid of cells", 
                            colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5); 
  fhGridCellsE->SetYTitle("row (phi direction)");
  fhGridCellsE->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCellsE);
  
  fhGridCellsLowGain  = new TH2F ("hGridCellsLG",Form("Low gain: Entries in grid of cells"),
                           colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5);
  fhGridCellsLowGain->SetYTitle("row (phi direction)");
  fhGridCellsLowGain->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCellsLowGain);
  
  fhGridCellsELowGain  = new TH2F ("hGridCellsELG","Low gain: Accumulated energy in grid of cells",
                            colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5);
  fhGridCellsELowGain->SetYTitle("row (phi direction)");
  fhGridCellsELowGain->SetXTitle("column (eta direction)");
  outputContainer->Add(fhGridCellsELowGain);
  
  if(fFillAllCellTimeHisto)
  { 
    fhGridCellsTime  = new TH2F ("hGridCellsTime","Accumulated time in grid of cells",
                                 colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5);
    fhGridCellsTime->SetYTitle("row (phi direction)");
    fhGridCellsTime->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsTime);
    
    fhGridCellsTimeLowGain  = new TH2F ("hGridCellsTimeLG","Low gain: Accumulated time in grid of cells",
                                        colmaxs+2,-1.5,colmaxs+0.5, rowmaxs+2,-1.5,rowmaxs+0.5);
    fhGridCellsTimeLowGain->SetYTitle("row (phi direction)");
    fhGridCellsTimeLowGain->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsTimeLowGain);
  }
  
  fhNCellsSumAmpPerMod       = new TH2F*[fNModules];  
  
  fhNClustersSumEnergyPerMod = new TH2F*[fNModules];
  
  fhIMMod                    = new TH2F*[fNModules];
  
  fhNCellsPerClusterMod      = new TH2F*[fNModules];
  
  if(fStudyBadClusters) 
    fhNCellsPerClusterModNoCut = new TH2F*[fNModules];
  
  if(fFillAllCellTimeHisto) 
    fhTimeAmpPerRCU = new TH2F*[fNModules*fNRCU];       
  
  for(Int_t imod = 0; imod < fNModules; imod++)
  {
    fhNCellsSumAmpPerMod[imod] = new TH2F (Form("hNCellsSumAmp_Mod%d",imod),
                                             Form("# cells in SM vs sum of cells energy in Module %d",imod), 
                                             nptbins,ptmin,ptmax*4, ncebins,ncemin,ncemax); 
    fhNCellsSumAmpPerMod[imod]->SetXTitle("#Sigma #it{Amplitude} (GeV)");
    fhNCellsSumAmpPerMod[imod]->SetYTitle("#Sigma #it{n}_{cells}");
    outputContainer->Add(fhNCellsSumAmpPerMod[imod]);

    fhNClustersSumEnergyPerMod[imod] = new TH2F (Form("hNClustersSumEnergy_Mod%d",imod),
                                           Form("# clusters in SM vs sum of clusters energy in Module %d",imod), 
                                           nptbins,ptmin,ptmax*4, nclbins,nclmin,nclmax); 
    fhNClustersSumEnergyPerMod[imod]->SetXTitle("#Sigma #it{E} (GeV)");
    fhNClustersSumEnergyPerMod[imod]->SetYTitle("#Sigma #it{n}_{clusters}");
    outputContainer->Add(fhNClustersSumEnergyPerMod[imod]);
    
    fhNCellsPerClusterMod[imod]  = new TH2F (Form("hNCellsPerCluster_Mod%d",imod),
                                             Form("# cells per cluster vs cluster energy in Module %d",imod), 
                                             nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerClusterMod[imod]->SetXTitle("#it{E} (GeV)");
    fhNCellsPerClusterMod[imod]->SetYTitle("#it{n}_{cells per cluster}");
    outputContainer->Add(fhNCellsPerClusterMod[imod]);
    
    if(fStudyBadClusters)
    {
      fhNCellsPerClusterModNoCut[imod]  = new TH2F (Form("hNCellsPerClusterNoCut_Mod%d",imod),
                                                    Form("# cells per cluster vs cluster energy in Module %d, no cut",imod), 
                                                    nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
      fhNCellsPerClusterModNoCut[imod]->SetXTitle("#it{E} (GeV)");
      fhNCellsPerClusterModNoCut[imod]->SetYTitle("#it{n}_{cells per cluster}");
      outputContainer->Add(fhNCellsPerClusterModNoCut[imod]);
    }
    
    if(fFillAllCellTimeHisto) 
    {
      for(Int_t ircu = 0; ircu < fNRCU; ircu++)
      {
        fhTimeAmpPerRCU[imod*fNRCU+ircu]  = new TH2F (Form("hTimeAmp_Mod%d_RCU%d",imod,ircu),
                                                      Form("#it{E}_{cell} vs #it{t}_{cell} in Module %d, RCU %d ",imod,ircu), 
                                                      nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetXTitle("#it{E} (GeV)");
        fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetYTitle("#it{t} (ns)");
        outputContainer->Add(fhTimeAmpPerRCU[imod*fNRCU+ircu]);
      }
    }
    
    if(fFillAllPi0Histo)
    {
      fhIMMod[imod]  = new TH2F (Form("hIM_Mod%d",imod),
                                 Form("Cluster pairs Invariant mass vs reconstructed pair energy in Module %d, n cell > 1",imod),
                                 nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
      fhIMMod[imod]->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
      fhIMMod[imod]->SetYTitle("#it{M}_{cluster pairs} (GeV/#it{c}^{2})");
      outputContainer->Add(fhIMMod[imod]);
    }
  }
  
  // Monte Carlo Histograms
  
  TString particleName[] = {
    "Photon",        "Pi0",         "Eta",
    "Electron",      "PhotonConv",
    "NeutralHadron", "ChargedHadron"      };

  if(IsDataMC())
  {
    for(Int_t iPart = 0; iPart < 7; iPart++)
    {
      for(Int_t iCh = 0; iCh < 2; iCh++)
      {
        fhRecoMCRatioE[iPart][iCh]  = new TH2F (Form("hRecoMCRatioE_%s_Match%d",particleName[iPart].Data(),iCh),
                                                Form("Reconstructed/Generated E, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                                nptbins, ptmin, ptmax, 200,0,2); 
        fhRecoMCRatioE[iPart][iCh]->SetYTitle("#it{E}_{reconstructed}/#it{E}_{generated}");
        fhRecoMCRatioE[iPart][iCh]->SetXTitle("#it{E}_{reconstructed} (GeV)");
        outputContainer->Add(fhRecoMCRatioE[iPart][iCh]);
        
        
        fhRecoMCDeltaE[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaE_%s_Match%d",particleName[iPart].Data(),iCh),
                                                Form("Generated - Reconstructed E, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                                nptbins, ptmin, ptmax, nptbins*2,-ptmax,ptmax); 
        fhRecoMCDeltaE[iPart][iCh]->SetYTitle("#Delta #it{E} (GeV)");
        fhRecoMCDeltaE[iPart][iCh]->SetXTitle("#it{E}_{reconstructed} (GeV)");
        outputContainer->Add(fhRecoMCDeltaE[iPart][iCh]);
        
        fhRecoMCDeltaPhi[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaPhi_%s_Match%d",particleName[iPart].Data(),iCh),
                                                  Form("Generated - Reconstructed #phi, %s, Matched %d",particleName[iPart].Data(),iCh),
                                                  nptbins, ptmin, ptmax, nphibins*2,-phimax,phimax); 
        fhRecoMCDeltaPhi[iPart][iCh]->SetYTitle("#Delta #phi (rad)");
        fhRecoMCDeltaPhi[iPart][iCh]->SetXTitle("#it{E}_{reconstructed} (GeV)");
        outputContainer->Add(fhRecoMCDeltaPhi[iPart][iCh]);
        
        fhRecoMCDeltaEta[iPart][iCh]  = new TH2F (Form("hRecoMCDeltaEta_%s_Match%d",particleName[iPart].Data(),iCh),
                                                  Form("Generated - Reconstructed #eta, %s, Matched %d",particleName[iPart].Data(),iCh),
                                                  nptbins, ptmin, ptmax,netabins*2,-etamax,etamax); 
        fhRecoMCDeltaEta[iPart][iCh]->SetYTitle("#Delta #eta ");
        fhRecoMCDeltaEta[iPart][iCh]->SetXTitle("#it{E}_{reconstructed} (GeV)");
        outputContainer->Add(fhRecoMCDeltaEta[iPart][iCh]);
        
        fhRecoMCE[iPart][iCh]  = new TH2F (Form("hRecoMCE_%s_Match%d",particleName[iPart].Data(),iCh),
                                           Form("#it{E} distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh),
                                           nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhRecoMCE[iPart][iCh]->SetXTitle("#it{E}_{rec} (GeV)");
        fhRecoMCE[iPart][iCh]->SetYTitle("#it{E}_{gen} (GeV)");
        outputContainer->Add(fhRecoMCE[iPart][iCh]);	  
        
        fhRecoMCPhi[iPart][iCh]  = new TH2F (Form("hRecoMCPhi_%s_Match%d",particleName[iPart].Data(),iCh),
                                             Form("#phi distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh),
                                             nphibins,phimin,phimax, nphibins,phimin,phimax); 
        fhRecoMCPhi[iPart][iCh]->SetXTitle("#phi_{reconstructed} (rad)");
        fhRecoMCPhi[iPart][iCh]->SetYTitle("#phi_{generated} (rad)");
        outputContainer->Add(fhRecoMCPhi[iPart][iCh]);
        
        fhRecoMCEta[iPart][iCh]  = new TH2F (Form("hRecoMCEta_%s_Match%d",particleName[iPart].Data(),iCh),
                                             Form("#eta distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                             netabins,etamin,etamax,netabins,etamin,etamax); 
        fhRecoMCEta[iPart][iCh]->SetXTitle("#eta_{reconstructed} ");
        fhRecoMCEta[iPart][iCh]->SetYTitle("#eta_{generated} ");
        outputContainer->Add(fhRecoMCEta[iPart][iCh]);
      }
    }  
    
    // Pure MC
      
    for(Int_t iPart = 0; iPart < 4; iPart++)
    {
      fhGenMCE [iPart]     = new TH1F(Form("hGenMCE_%s",particleName[iPart].Data()) ,
                                      Form("#it{E} of generated %s",particleName[iPart].Data()),
                                      nptbins,ptmin,ptmax);
      
      fhGenMCPt[iPart]     = new TH1F(Form("hGenMCPt_%s",particleName[iPart].Data()) ,
                                     Form("#it{p}_{T} of generated %s",particleName[iPart].Data()),
                                     nptbins,ptmin,ptmax);

      fhGenMCEtaPhi[iPart] = new TH2F(Form("hGenMCEtaPhi_%s",particleName[iPart].Data()),
                                      Form("Y vs #phi of generated %s",particleName[iPart].Data()),
                                      200,-1,1,360,0,TMath::TwoPi());
    	
      fhGenMCE [iPart]    ->SetXTitle("#it{E} (GeV)");
      fhGenMCPt[iPart]    ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhGenMCEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCEtaPhi[iPart]->SetYTitle("#phi (rad)");

      outputContainer->Add(fhGenMCE     [iPart]);
      outputContainer->Add(fhGenMCPt    [iPart]);
      outputContainer->Add(fhGenMCEtaPhi[iPart]);
      
      
      fhGenMCAccE [iPart]     = new TH1F(Form("hGenMCAccE_%s",particleName[iPart].Data()) ,
                                        Form("#it{E} of generated %s",particleName[iPart].Data()),
                                        nptbins,ptmin,ptmax);
      fhGenMCAccPt[iPart]     = new TH1F(Form("hGenMCAccPt_%s",particleName[iPart].Data()) ,
                                        Form("#it{p}_{T} of generated %s",particleName[iPart].Data()),
                                        nptbins,ptmin,ptmax);
      fhGenMCAccEtaPhi[iPart] = new TH2F(Form("hGenMCAccEtaPhi_%s",particleName[iPart].Data()),
                                         Form("Y vs #phi of generated %s",particleName[iPart].Data()),
                                         netabins,etamin,etamax,nphibins,phimin,phimax);
    	
      fhGenMCAccE [iPart]    ->SetXTitle("#it{E} (GeV)");
      fhGenMCAccPt[iPart]    ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhGenMCAccEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCAccEtaPhi[iPart]->SetYTitle("#phi (rad)");
      
      outputContainer->Add(fhGenMCAccE     [iPart]);
      outputContainer->Add(fhGenMCAccPt    [iPart]);
      outputContainer->Add(fhGenMCAccEtaPhi[iPart]);
      
    }    
    
    // Vertex of generated particles
    
    fhEMVxyz  = new TH2F ("hEMVxyz","Production vertex of reconstructed ElectroMagnetic particles",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhEMVxyz->SetXTitle("#it{v}_{x}");
    fhEMVxyz->SetYTitle("#it{v}_{y}");
    //fhEMVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhEMVxyz);
    
    fhHaVxyz  = new TH2F ("hHaVxyz","Production vertex of reconstructed hadrons",nvdistbins,vdistmin,vdistmax,nvdistbins,vdistmin,vdistmax);//,100,0,500); 
    fhHaVxyz->SetXTitle("#it{v}_{x}");
    fhHaVxyz->SetYTitle("#it{v}_{y}");
    //fhHaVxyz->SetZTitle("v_{z}");
    outputContainer->Add(fhHaVxyz);
    
    fhEMR  = new TH2F ("hEMR","Distance to production vertex of reconstructed ElectroMagnetic particles vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhEMR->SetXTitle("#it{E} (GeV)");
    fhEMR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhEMR);
    
    fhHaR  = new TH2F ("hHaR","Distance to production vertex of reconstructed Hadrons vs E rec",nptbins,ptmin,ptmax,nvdistbins,vdistmin,vdistmax); 
    fhHaR->SetXTitle("#it{E} (GeV)");
    fhHaR->SetYTitle("TMath::Sqrt(v_{x}^{2}+v_{y}^{2})");
    outputContainer->Add(fhHaR);
    
    // Track Matching
    
    fhMCEle1EOverP = new TH2F("hMCEle1EOverP","TRACK matches #it{E}/#it{p}, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCEle1EOverP->SetYTitle("#it{E}/#it{p}");
    fhMCEle1EOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCEle1EOverP);
    
    fhMCEle1dR = new TH1F("hMCEle1dR","TRACK matches dR, MC electrons",ndRbins,dRmin,dRmax);
    fhMCEle1dR->SetXTitle("#Delta #it{R} (rad)");
    outputContainer->Add(fhMCEle1dR) ;
    
    fhMCEle2MatchdEdx = new TH2F("hMCEle2MatchdEdx","#it{dE/dx} vs. #it{p} for all matches, MC electrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCEle2MatchdEdx->SetXTitle("#it{p} (GeV/#it{c})");
    fhMCEle2MatchdEdx->SetYTitle("<#it{dE/dx}>");
    outputContainer->Add(fhMCEle2MatchdEdx);
    
    fhMCChHad1EOverP = new TH2F("hMCChHad1EOverP","TRACK matches #it{E}/#it{p}, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCChHad1EOverP->SetYTitle("#it{E}/#it{p}");
    fhMCChHad1EOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCChHad1EOverP);
    
    fhMCChHad1dR = new TH1F("hMCChHad1dR","TRACK matches dR, MC charged hadrons",ndRbins,dRmin,dRmax);
    fhMCChHad1dR->SetXTitle("#Delta R (rad)");
    outputContainer->Add(fhMCChHad1dR) ;
    
    fhMCChHad2MatchdEdx = new TH2F("hMCChHad2MatchdEdx","#it{dE/dx} vs. #it{p} for all matches, MC charged hadrons",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCChHad2MatchdEdx->SetXTitle("#it{p} (GeV/#it{c})");
    fhMCChHad2MatchdEdx->SetYTitle("#it{dE/dx}>");
    outputContainer->Add(fhMCChHad2MatchdEdx);
    
    fhMCNeutral1EOverP = new TH2F("hMCNeutral1EOverP","TRACK matches #it{E}/#it{p}, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCNeutral1EOverP->SetYTitle("#it{E}/#it{p}");
    fhMCNeutral1EOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCNeutral1EOverP);
    
    fhMCNeutral1dR = new TH1F("hMCNeutral1dR","TRACK matches dR, MC neutrals",ndRbins,dRmin,dRmax);
    fhMCNeutral1dR->SetXTitle("#Delta #it{R} (rad)");
    outputContainer->Add(fhMCNeutral1dR) ;
    
    fhMCNeutral2MatchdEdx = new TH2F("hMCNeutral2MatchdEdx","#it{dE/dx} vs. #it{p} for all matches, MC neutrals",nptbins,ptmin,ptmax,ndedxbins,dedxmin,dedxmax);
    fhMCNeutral2MatchdEdx->SetXTitle("#it{p} (GeV/#it{c})");
    fhMCNeutral2MatchdEdx->SetYTitle("#it{dE/dx}>");
    outputContainer->Add(fhMCNeutral2MatchdEdx);
    
    fhMCEle1EOverPR02 = new TH2F("hMCEle1EOverPR02","TRACK matches #it{E}/#it{p}, MC electrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCEle1EOverPR02->SetYTitle("#it{E}/#it{p}");
    fhMCEle1EOverPR02->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCEle1EOverPR02);
    
    fhMCChHad1EOverPR02 = new TH2F("hMCChHad1EOverPR02","TRACK matches #it{E}/#it{p}, MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCChHad1EOverPR02->SetYTitle("#it{E}/#it{p}");
    fhMCChHad1EOverPR02->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCChHad1EOverPR02);
    
    fhMCNeutral1EOverPR02 = new TH2F("hMCNeutral1EOverPR02","TRACK matches #it{E}/#it{p}, MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCNeutral1EOverPR02->SetYTitle("#it{E}/#it{p}");
    fhMCNeutral1EOverPR02->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCNeutral1EOverPR02);
    
    fhMCEle1EleEOverP = new TH2F("hMCEle1EleEOverP","Electron candidates #it{E}/#it{p} (60<dEdx<100), MC electrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCEle1EleEOverP->SetYTitle("#it{E}/#it{p}");
    fhMCEle1EleEOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCEle1EleEOverP);

    fhMCChHad1EleEOverP = new TH2F("hMCEle1EleEOverP","Electron candidates #it{E}/#it{p} (60<dEdx<100), MC charged hadrons",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCChHad1EleEOverP->SetYTitle("#it{E}/#it{p}");
    fhMCChHad1EleEOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCChHad1EleEOverP);

    fhMCNeutral1EleEOverP = new TH2F("hMCNeutral1EleEOverP","Electron candidates #it{E}/#it{p} (60<dEdx<100), MC neutrals",nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhMCNeutral1EleEOverP->SetYTitle("#it{E}/#it{p}");
    fhMCNeutral1EleEOverP->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhMCNeutral1EleEOverP);
  }
  
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}

//______________________________________________________________________________________
/// Get energy in cross axis around maximum cell. For both calorimeters.
/// For EMCal, same procedure as in AliCalorimeterUtils
//______________________________________________________________________________________
Float_t AliAnaCalorimeterQA::GetECross(Int_t absID, AliVCaloCells* cells, Float_t dtcut)
{
  Int_t icol =-1, irow=-1,iRCU = -1;   
  Int_t imod = GetModuleNumberCellIndexes(absID, GetCalorimeter(), icol, irow, iRCU);
    
  if ( GetCalorimeter() == kEMCAL)
  {
    // Get close cells index, energy and time, not in corners
    
    Int_t absID1 = -1;
    Int_t absID2 = -1;
    
    if( irow < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = GetCaloUtils()->GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(imod, irow+1, icol);
    if( irow > 0 )                                absID2 = GetCaloUtils()->GetEMCALGeometry()->GetAbsCellIdFromCellIndexes(imod, irow-1, icol);
    
    // In case of cell in eta = 0 border, depending on SM shift the cross cell index
    Int_t absID3 = -1;
    Int_t absID4 = -1;
    
    if     ( icol == AliEMCALGeoParams::fgkEMCALCols - 1 && !(imod%2) )
    {
      absID3 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod+1, irow, 0);
      absID4 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod  , irow, icol-1); 
    }
    else if( icol == 0 && imod%2 )
    {
      absID3 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod  , irow, icol+1);
      absID4 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod-1, irow, AliEMCALGeoParams::fgkEMCALCols-1); 
    }
    else
    {
      if( icol < AliEMCALGeoParams::fgkEMCALCols-1 )
        absID3 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod, irow, icol+1);
      if( icol > 0 )    
        absID4 = GetCaloUtils()->GetEMCALGeometry()-> GetAbsCellIdFromCellIndexes(imod, irow, icol-1);
    }
    
    // Recalibrate cell energy if needed
    //Float_t  ecell = cells->GetCellAmplitude(absID);
    //GetCaloUtils()->RecalibrateCellAmplitude(ecell,GetCalorimeter(), absID);
    Double_t tcell = cells->GetCellTime(absID);
    GetCaloUtils()->RecalibrateCellTime(tcell, GetCalorimeter(), absID,GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    
    Float_t  ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
    Double_t tcell1  = 0, tcell2  = 0, tcell3  = 0, tcell4  = 0;
    
    if(absID1 > 0 )
    {
      ecell1 = cells->GetCellAmplitude(absID1);
      GetCaloUtils()->RecalibrateCellAmplitude(ecell1, GetCalorimeter(), absID1);
      tcell1 = cells->GetCellTime(absID1);
      GetCaloUtils()->RecalibrateCellTime     (tcell1, GetCalorimeter(), absID1,GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    }
      
    if(absID2 > 0 )
    {
      ecell2 = cells->GetCellAmplitude(absID2);
      GetCaloUtils()->RecalibrateCellAmplitude(ecell2, GetCalorimeter(), absID2);
      tcell2 = cells->GetCellTime(absID2);
      GetCaloUtils()->RecalibrateCellTime     (tcell2, GetCalorimeter(), absID2, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    }
      
    if(absID3 > 0 )
    {
      ecell3 = cells->GetCellAmplitude(absID3);
      GetCaloUtils()->RecalibrateCellAmplitude(ecell3, GetCalorimeter(), absID3);
      tcell3 = cells->GetCellTime(absID3);
      GetCaloUtils()->RecalibrateCellTime     (tcell3, GetCalorimeter(), absID3, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    }
      
    if(absID4 > 0 )
    {
      ecell4 = cells->GetCellAmplitude(absID4);
      GetCaloUtils()->RecalibrateCellAmplitude(ecell4, GetCalorimeter(), absID4);
      tcell4 = cells->GetCellTime(absID4);
      GetCaloUtils()->RecalibrateCellTime     (tcell4, GetCalorimeter(), absID4, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    }
        
    if(TMath::Abs(tcell-tcell1)*1.e9 > dtcut) ecell1 = 0 ;
    if(TMath::Abs(tcell-tcell2)*1.e9 > dtcut) ecell2 = 0 ;
    if(TMath::Abs(tcell-tcell3)*1.e9 > dtcut) ecell3 = 0 ;
    if(TMath::Abs(tcell-tcell4)*1.e9 > dtcut) ecell4 = 0 ;
    
    return ecell1+ecell2+ecell3+ecell4;
  }
  else // PHOS
  { 
    Int_t absId1 = -1, absId2 = -1, absId3 = -1, absId4 = -1;
    
    Int_t relId1[] = { imod+1, 0, irow+1, icol   };
    Int_t relId2[] = { imod+1, 0, irow-1, icol   };
    Int_t relId3[] = { imod+1, 0, irow  , icol+1 };
    Int_t relId4[] = { imod+1, 0, irow  , icol-1 };
    
    GetCaloUtils()->GetPHOSGeometry()->RelToAbsNumbering(relId1, absId1);
    GetCaloUtils()->GetPHOSGeometry()->RelToAbsNumbering(relId2, absId2);
    GetCaloUtils()->GetPHOSGeometry()->RelToAbsNumbering(relId3, absId3);
    GetCaloUtils()->GetPHOSGeometry()->RelToAbsNumbering(relId4, absId4);
    
    Float_t  ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
    
    if(absId1 > 0 ) ecell1 = cells->GetCellAmplitude(absId1);
    if(absId2 > 0 ) ecell2 = cells->GetCellAmplitude(absId2);
    if(absId3 > 0 ) ecell3 = cells->GetCellAmplitude(absId3);
    if(absId4 > 0 ) ecell4 = cells->GetCellAmplitude(absId4);
    
    return ecell1+ecell2+ecell3+ecell4;
  }
}

//___________________________________________________________________________________________________________
/// Fill Invariant mass histograms.
/// \param iclus: main loop cluster index
/// \param nModule: (super)module number of iclus.
/// \param caloClusters: array with the clusters.
/// \param cells: array with cells.
//___________________________________________________________________________________________________________
void AliAnaCalorimeterQA::InvariantMassHistograms(Int_t iclus,  Int_t nModule, const TObjArray* caloClusters,
                                                  AliVCaloCells * cells) 
{
  AliDebug(1,"Start");
  
  //Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
  //GetReader()->GetVertex(v);
  
  Int_t nModule2      = -1;
  Int_t nCaloClusters = caloClusters->GetEntriesFast();
  
  Float_t phi1 = fClusterMomentum.Phi();
  if(phi1 < 0) phi1 += TMath::TwoPi();    
  
  Double_t tof1 =  ((AliVCluster*) caloClusters->At(iclus))->GetTOF()*1.e9;
  if(tof1>400) tof1-=fConstantTimeShift;

  for(Int_t jclus = iclus + 1 ; jclus < nCaloClusters ; jclus++) 
  {
    AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(jclus);

    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus2, maxCellFraction);

    Double_t tof2 =  clus2->GetTOF()*1.e9;
    if(tof2>400) tof2-=fConstantTimeShift;

    Double_t diffTof = tof1-tof2;
    
    // Try to reduce background with a mild shower shape cut and no more 
    // than 1 local maximum in cluster and remove low energy clusters

    if(  !IsGoodCluster(absIdMax, clus2->GetM02(), clus2->GetNCells(), cells) 
       || GetCaloUtils()->GetNumberOfLocalMaxima(clus2,cells) > 1 
       || clus2->GetM02() > fInvMassMaxM02Cut
       || clus2->GetM02() < fInvMassMinM02Cut
       || clus2->E() < fInvMassMinECut 
       || clus2->E() > fInvMassMaxECut  
       || TMath::Abs(diffTof) > fInvMassMaxTimeDifference
       ) continue;
    
    // Get cluster kinematics
    clus2->GetMomentum(fClusterMomentum2,v);
    
    // Check only certain regions
    Bool_t in2 = kTRUE;
    if(IsFiducialCutOn()) in2 =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum2.Eta(),fClusterMomentum2.Phi(),GetCalorimeter()) ;
    if(!in2) continue;	

    Float_t  pairPt = (fClusterMomentum+fClusterMomentum2).Pt();
    
    // Opening angle cut, avoid combination of DCal and EMCal clusters
    Double_t angle  = fClusterMomentum.Angle(fClusterMomentum2.Vect());
     
    if ( fFillInvMassOpenAngle ) fhOpAngle->Fill(pairPt, angle, GetEventWeight()) ;
    
    if( angle > fInvMassMaxOpenAngle ) continue;
        
    // Get module of cluster
    nModule2 = GetModuleNumber(clus2);
    
    // Fill histograms
    Float_t mass   = (fClusterMomentum+fClusterMomentum2).M ();
    Float_t asym   = TMath::Abs( fClusterMomentum.E() - fClusterMomentum2.E() ) /
                               ( fClusterMomentum.E() + fClusterMomentum2.E() );

    // All modules
    // Combine EMCal or PHOS clusters
        
    Float_t phi2 = fClusterMomentum2.Phi();
    if(phi2 < 0) phi2 += TMath::TwoPi();
    
    AliDebug(1,Form("Selected pair: pT %f, mass %f, asym %f, angle %f, diffTof %f, SM1 %d, SM2 %d\n",pairPt,mass,asym,angle,diffTof,nModule,nModule2));
            
    Bool_t inPi0Window = kFALSE;
    if(mass < 0.18 && mass > 0.1)  inPi0Window = kTRUE ;
    
    if ( nModule < 12 && nModule2 < 12 ) 
    {
      fhIM    ->Fill(pairPt, mass, GetEventWeight());
            
      if( fFillPi0PairDiffTime && inPi0Window )
        fhClusterPairDiffTimeEPi0Mass->Fill(pairPt,  diffTof, GetEventWeight());
      
      if ( nModule == nModule2 ) 
      {
        fhIMSame->Fill(pairPt, mass, GetEventWeight());
                
        if( fFillPi0PairDiffTime && inPi0Window )
          fhClusterPairDiffTimeEPi0MassSame->Fill(pairPt,  diffTof, GetEventWeight());
      }
      else        
      {
        fhIMDiff->Fill(pairPt, mass, GetEventWeight());
      }
    }
    // Combine DCal clusters
    else if ( ( GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL ) &&
            nModule > 11 && nModule2 > 11  && fNModules > 12 ) 
    {
      fhIMDCAL->Fill(pairPt, mass, GetEventWeight());
      
      if( fFillPi0PairDiffTime && inPi0Window )
        fhClusterPairDiffTimeEPi0MassDCal->Fill(pairPt,  diffTof, GetEventWeight());
      
      if ( nModule == nModule2 )
      {
        fhIMDCALSame->Fill(pairPt, mass, GetEventWeight());
        
        if( fFillPi0PairDiffTime && inPi0Window )
          fhClusterPairDiffTimeEPi0MassDCalSame->Fill(pairPt,  diffTof, GetEventWeight());
      }
      else         
      {
        fhIMDCALDiff->Fill(pairPt, mass, GetEventWeight());
      }
    } 
    
    if ( fFillInvMassOpenAngle ) fhIMvsOpAngle->Fill(mass, angle, GetEventWeight());

    // Single module
    if(nModule == nModule2 && nModule >= 0 && nModule < fNModules)
      fhIMMod[nModule]->Fill(pairPt, mass, GetEventWeight());
    
    // Asymetry histograms
    fhAsym->Fill(pairPt, asym, GetEventWeight());
  } // 2nd cluster loop
  
  //
  // Combine PHOS and DCal
  //
  
  if( ( GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL ) &&
     fNModules > 12 && nModule > 11)
  {
    AliDebug(1,"Check DCal-PHOS pairs\n");

    Int_t sector1 = -1;
    Int_t sector2 = -1;
    
    if(phi1 >= 260*TMath::DegToRad() && phi1 < 280) sector1 = 0;
    if(phi1 >= 280*TMath::DegToRad() && phi1 < 300) sector1 = 1;
    if(phi1 >= 300*TMath::DegToRad() && phi1 < 320) sector1 = 2;
    
    for(Int_t jclus = 0 ; jclus < GetPHOSClusters()->GetEntriesFast() ; jclus++) 
    {
      AliVCluster* clus2 =  (AliVCluster*) GetPHOSClusters()->At(jclus);
            
      Float_t maxCellFraction = 0.;
      Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus2, maxCellFraction);
      
      // Try to reduce background, remove low energy clusters
      if(  !IsGoodCluster(absIdMax, clus2->GetM02(), clus2->GetNCells(), cells) 
         || clus2->E() < fInvMassMinECut 
         || clus2->E() > fInvMassMaxECut  
         ) continue;
      
      // Get cluster kinematics
      clus2->GetMomentum(fClusterMomentum2,v);
      
      // Fill histograms
      
      Float_t mass   = (fClusterMomentum+fClusterMomentum2).M ();
      Float_t pairPt = (fClusterMomentum+fClusterMomentum2).Pt();
      //Float_t asym   = TMath::Abs( fClusterMomentum.E() - fClusterMomentum2.E() ) /
      //( fClusterMomentum.E() + fClusterMomentum2.E() );
     
      fhIMDCALPHOS->Fill(pairPt, mass, GetEventWeight());
      
      Float_t phiPHOS = fClusterMomentum2.Phi();
      if(phiPHOS < 0) phiPHOS += TMath::TwoPi();
      
      if(phiPHOS >= 260*TMath::DegToRad() && phiPHOS < 280) sector2 = 0;
      if(phiPHOS >= 280*TMath::DegToRad() && phiPHOS < 300) sector2 = 1;
      if(phiPHOS >= 300*TMath::DegToRad() && phiPHOS < 320) sector2 = 2;
      
      if(sector1 == sector2) fhIMDCALPHOSSame->Fill(pairPt, mass, GetEventWeight());
      
    } // PHOS cluster loop
  } // DCal-PHOS combination


   // Select EMCal clusters in DCal eta acceptance
   // Cross check combination of DCal-PHOS pairs
  if( fFillInvMassInEMCALWithPHOSDCalAcc && TMath::Abs(fClusterMomentum.Eta() > 0.22))
  {
    AliDebug(1,"Check EMCAL(DCal)-EMCAL(PHOS) pairs\n");

    Int_t sector1 = -1;
    Int_t sector2 = -1;
    
    if(phi1 >=  80*TMath::DegToRad() && phi1 < 100) sector1 = 0;
    if(phi1 >= 100*TMath::DegToRad() && phi1 < 120) sector1 = 1;
    if(phi1 >= 120*TMath::DegToRad() && phi1 < 140) sector1 = 2;
    if(phi1 >= 140*TMath::DegToRad() && phi1 < 160) sector1 = 3;
    if(phi1 >= 160*TMath::DegToRad() && phi1 < 180) sector1 = 4;
    if(phi1 >= 180*TMath::DegToRad() && phi1 < 190) sector1 = 5;
    
    for(Int_t jclus = 0 ; jclus < caloClusters->GetEntriesFast() ; jclus++) 
    {
      AliVCluster* clus2 =  (AliVCluster*) caloClusters->At(jclus);
      
      Float_t maxCellFraction = 0.;
      Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus2, maxCellFraction);
      
      Double_t tof2 =  clus2->GetTOF()*1.e9;
      if(tof2>400) tof2-=fConstantTimeShift;

      Double_t diffTof = TMath::Abs(tof1-tof2);
      
      // Try to reduce background with a mild shower shape cut and no more 
      // than 1 local maximum in cluster and remove low energy clusters
      if(  !IsGoodCluster(absIdMax, clus2->GetM02(), clus2->GetNCells(), cells) 
         || GetCaloUtils()->GetNumberOfLocalMaxima(clus2,cells) > 1 
         || clus2->GetM02() > fInvMassMaxM02Cut 
         || clus2->GetM02() < fInvMassMinM02Cut
         || clus2->E() < fInvMassMinECut 
         || clus2->E() > fInvMassMaxECut  
         || TMath::Abs(diffTof) > fInvMassMaxTimeDifference
         ) continue;
      
      // Get cluster kinematics
      clus2->GetMomentum(fClusterMomentum2,v);
      
      // Select EMCal clusters in PHOS acceptance
      if(TMath::Abs(fClusterMomentum2.Eta() > 0.13)) continue ; 
      
      // Fill histograms
      
      Float_t mass   = (fClusterMomentum+fClusterMomentum2).M ();
      Float_t pairPt = (fClusterMomentum+fClusterMomentum2).Pt();
      //Float_t asym   = TMath::Abs( fClusterMomentum.E() - fClusterMomentum2.E() ) /
      //( fClusterMomentum.E() + fClusterMomentum2.E() );
      
      fhIMEMCALPHOS->Fill(pairPt, mass, GetEventWeight());
      
      Float_t phiPHOS = fClusterMomentum2.Phi();
      if(phiPHOS < 0) phiPHOS += TMath::TwoPi();
      
      if(phiPHOS >=  80*TMath::DegToRad() && phiPHOS < 100) sector2 = 0;
      if(phiPHOS >= 100*TMath::DegToRad() && phiPHOS < 120) sector2 = 1;
      if(phiPHOS >= 120*TMath::DegToRad() && phiPHOS < 140) sector2 = 2;
      if(phiPHOS >= 140*TMath::DegToRad() && phiPHOS < 160) sector2 = 3;
      if(phiPHOS >= 160*TMath::DegToRad() && phiPHOS < 180) sector2 = 4;
      if(phiPHOS >= 180*TMath::DegToRad() && phiPHOS < 190) sector2 = 5;
      
      if(sector1 == sector2) fhIMEMCALPHOSSame->Fill(pairPt, mass, GetEventWeight());
      
    } // PHOS cluster loop
  } // EMCal(DCal)-EMCAL(PHOS) combination

  AliDebug(1,"End");
}

//______________________________
// Check if the calorimeter setting is ok, if not abort.
//______________________________
void AliAnaCalorimeterQA::Init()
{
  if(GetCalorimeter() != kPHOS && GetCalorimeter() !=kEMCAL)
    AliFatal(Form("Wrong calorimeter name <%s>", GetCalorimeterString().Data()));
  
  //if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
  //  AliFatal("Analysis of reconstructed data, MC reader not aplicable");
}

//________________________________________
/// Initialize the parameters of the analysis.
//________________________________________
void AliAnaCalorimeterQA::InitParameters()
{ 
  AddToHistogramsName("AnaCaloQA_");
  
  fNModules        = 22; // set maximum to maximum number of EMCAL modules
  fNRCU            = 2;  // set maximum number of RCU in EMCAL per SM
  
  fTimeCutMin      = -9999999;
  fTimeCutMax      =  9999999;
  
  fEMCALCellAmpMin = 0.2; // 200 MeV
  fPHOSCellAmpMin  = 0.2; // 200 MeV
  fCellAmpMin      = 0.2; // 200 MeV

  fEMCALClusterNCellMin = 2; // at least 2
  fPHOSClusterNCellMin  = 3; // at least 3
  
  fInvMassMinECut  = 0.5; // 500 MeV
  fInvMassMaxECut  = 10; 

  fInvMassMinM02Cut = 0.1; 
  fInvMassMaxM02Cut = 0.4; 
  
  fInvMassMaxTimeDifference = 200; // ns, quite open 
  
  fInvMassMaxOpenAngle = 100*TMath::DegToRad(); // 100 degrees
  
  // Exotic studies
  fExoNECrossCuts  = 10 ;
  fExoNDTimeCuts   = 4  ;
  
  fExoDTimeCuts [0] = 1.e4 ; fExoDTimeCuts [1] = 50.0 ; fExoDTimeCuts [2] = 25.0 ; fExoDTimeCuts [3] = 10.0 ;
  fExoECrossCuts[0] = 0.80 ; fExoECrossCuts[1] = 0.85 ; fExoECrossCuts[2] = 0.90 ; fExoECrossCuts[3] = 0.92 ; fExoECrossCuts[4] = 0.94 ;
  fExoECrossCuts[5] = 0.95 ; fExoECrossCuts[6] = 0.96 ; fExoECrossCuts[7] = 0.97 ; fExoECrossCuts[8] = 0.98 ; fExoECrossCuts[9] = 0.99 ;
}

//_____________________________________________________________________________
/// Identify cluster as exotic or not.
//_____________________________________________________________________________
Bool_t AliAnaCalorimeterQA::IsGoodCluster(Int_t absIdMax, Float_t m02, 
                                          Int_t nCellsPerCluster, AliVCaloCells* cells)
{
  //if(!fStudyBadClusters) return kTRUE;
    
  if(GetCalorimeter() == kEMCAL) 
  {
 
    if(  m02 < 0.05 || nCellsPerCluster < fEMCALClusterNCellMin ) return kFALSE ; // mild shower shape cut for exotics

    if(!GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster())
    {
      return  !(GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(absIdMax,cells,(GetReader()->GetInputEvent())->GetBunchCrossNumber()));
    }
    else
    {
      return kTRUE;
    }
  }
  else // PHOS
  {
    if(  m02 < 0.05 || nCellsPerCluster < fPHOSClusterNCellMin ) return kFALSE ; // mild shower shape cut for exotics

    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    GetCaloUtils()->RecalibrateCellAmplitude(ampMax, GetCalorimeter(), absIdMax);
    
    if(ampMax < 0.01) return kFALSE;
    
    if(1-GetECross(absIdMax,cells)/ampMax > 0.95) return kFALSE;
    else                                          return kTRUE;
  }
}

//_________________________________________________________
/// Print some relevant parameters set for the analysis.
//_________________________________________________________
void AliAnaCalorimeterQA::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",GetCalorimeterString().Data());
  printf("Time Cut: %3.1f < TOF  < %3.1f\n"   , fTimeCutMin, fTimeCutMax);
  printf("EMCAL Min Amplitude : %2.1f GeV/c\n", fEMCALCellAmpMin) ;
  printf("PHOS  Min Amplitude : %2.1f GeV/c\n", fPHOSCellAmpMin) ;
  printf("EMCAL Min n cells   : %d\n"         , fEMCALClusterNCellMin) ;
  printf("PHOS  Min n cells   : %d\n"         , fPHOSClusterNCellMin) ;
  
  printf("Inv. Mass %2.1f < E_clus < %2.1f GeV/c\n"  , fInvMassMinECut  , fInvMassMaxECut  ) ;
  printf("Inv. Mass %2.1f < M02_clus < %2.1f GeV/c\n", fInvMassMinM02Cut, fInvMassMaxM02Cut) ;
  printf("Inv. Mass open angle : %2.1f deg\n"        , fInvMassMaxOpenAngle*TMath::RadToDeg()) ;
  printf("Inv. Mass time difference: %2.1f ns\n"     , fInvMassMaxTimeDifference) ;
}

//_____________________________________________________
/// Main task method, call all the methods filling QA histograms.
//_____________________________________________________
void  AliAnaCalorimeterQA::MakeAnalysisFillHistograms()
{
  AliDebug(1,"Start");

  //Print("");
  
  // Play with the MC stack if available
  if(IsDataMC()) MCHistograms();
  
  // Correlate Calorimeters and V0 and track Multiplicity
  if(fCorrelate)	Correlate();

  // Get List with CaloClusters , calo Cells, init min amplitude
  TObjArray     * caloClusters = NULL;
  AliVCaloCells * cells        = 0x0;
  if      (GetCalorimeter() == kPHOS)
  {
    fCellAmpMin  = fPHOSCellAmpMin;
    caloClusters = GetPHOSClusters();
    cells        = GetPHOSCells();
  }
  else if (GetCalorimeter() == kEMCAL)
  {
    fCellAmpMin  = fEMCALCellAmpMin;
    caloClusters = GetEMCALClusters();
    cells        = GetEMCALCells();
  }
  else
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END", GetCalorimeterString().Data()));
  
  if( !caloClusters || !cells )
  {
    AliFatal(Form("AliAnaCalorimeterQA::MakeAnalysisFillHistograms() - No CaloClusters or CaloCells available"));
    return; // trick coverity
  }
  
  if(caloClusters->GetEntriesFast() == 0) return ;
  
  //printf("QA: N cells %d, N clusters %d \n",cells->GetNumberOfCells(),caloClusters->GetEntriesFast());
  
  // Clusters
  ClusterLoopHistograms(caloClusters,cells);
  
  // Cells  
  CellHistograms(cells);
  
  AliDebug(1,"End");
}

//______________________________________
/// Access the generated particle kinematics
/// and fill the corresponding histograms.
/// Not dependent on cluster/cell.
//______________________________________
void AliAnaCalorimeterQA::MCHistograms()
{
  Int_t    pdg    =  0 ;
  Int_t    status =  0 ;
  Int_t    nprim  =  0 ;
  
  TParticle        * primStack = 0;
  AliAODMCParticle * primAOD   = 0;
  
  // Get the MC arrays and do some checks before filling MC histograms
    
  // Get the ESD MC particles container
  AliStack * stack = 0;
  if( GetReader()->ReadStack() )
  {
    stack = GetMCStack();
    if(!stack ) return;
    nprim = stack->GetNtrack();
  }
  
  // Get the AOD MC particles container
  TClonesArray * mcparticles = 0;
  if( GetReader()->ReadAODMCParticles() )
  {
    mcparticles = GetReader()->GetAODMCParticles();
    if( !mcparticles ) return;
    nprim = mcparticles->GetEntriesFast();
  }
  
  //printf("N primaries %d\n",nprim);
  for(Int_t i=0 ; i < nprim; i++)
  {
    if(GetReader()->AcceptOnlyHIJINGLabels() && !GetReader()->IsHIJINGLabel(i)) continue ;
    
    // Get the generated particles, check that it is primary (not parton, resonance)
    // and get its momentum. Different way to recover from ESD or AOD
    if(GetReader()->ReadStack())
    {
      primStack = stack->Particle(i) ;
      if(!primStack)
      {
        AliWarning("ESD primaries pointer not available!!");
        continue;
      }
      
      pdg    = primStack->GetPdgCode();
      status = primStack->GetStatusCode();
      
      //printf("Input: i %d, %s, pdg %d, status %d \n",i, primStack->GetName(), pdg, status);
      
      if ( status > 11 ) continue; //Working for PYTHIA and simple generators, check for HERWIG, HIJING?
      
      // Protection against floating point exception
      if ( primStack->Energy() == TMath::Abs(primStack->Pz()) ||
          (primStack->Energy() - primStack->Pz()) < 1e-3      ||
          (primStack->Energy() + primStack->Pz()) < 0           )  continue ; 
      
      //printf("Take : i %d, %s, pdg %d, status %d \n",i, primStack->GetName(), pdg, status);
      
      // Photon kinematics
      primStack->Momentum(fPrimaryMomentum);
    }
    else
    {
      primAOD = (AliAODMCParticle *) mcparticles->At(i);
      if(!primAOD)
      {
        AliWarning("AOD primaries pointer not available!!");
        continue;
      }
      
      pdg    = primAOD->GetPdgCode();
      status = primAOD->GetStatus();
      
      //printf("Input: i %d, %s, pdg %d, status %d \n",i, primAOD->GetName(), pdg, status);

      if (!primAOD->IsPrimary()) continue; //accept all which is not MC transport generated. Don't know how to avoid partons
      
      if ( status > 11 ) continue; //Working for PYTHIA and simple generators, check for HERWIG
   
      //Protection against floating point exception
      if ( primAOD->E() == TMath::Abs(primAOD->Pz()) || 
          (primAOD->E() - primAOD->Pz()) < 1e-3      || 
          (primAOD->E() + primAOD->Pz()) < 0           )  continue ; 

      //printf("Take : i %d, %s, pdg %d, status %d \n",i, primAOD->GetName(), pdg, status);
      
      // Kinematics
      fPrimaryMomentum.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());
    }

    Float_t eMC    = fPrimaryMomentum.E();
    if(eMC < 0.2) continue;
    
    Float_t ptMC   = fPrimaryMomentum.E();
    
    Float_t etaMC  = fPrimaryMomentum.Eta();
      
    // Only particles in |eta| < 1
    if (TMath::Abs(etaMC) > 1) continue;
    
    Float_t phiMC  = fPrimaryMomentum.Phi();
    if(phiMC < 0)
      phiMC  += TMath::TwoPi();
    
    Int_t mcIndex = -1;
    if      (pdg==22)  mcIndex = kmcPhoton;
    else if (pdg==111) mcIndex = kmcPi0;
    else if (pdg==221) mcIndex = kmcEta;
    else if (TMath::Abs(pdg)==11) mcIndex = kmcElectron;
    
    if( mcIndex >=0 )
    {
      fhGenMCE [mcIndex]->Fill( eMC, GetEventWeight());
      fhGenMCPt[mcIndex]->Fill(ptMC, GetEventWeight());
      if(eMC > 0.5) fhGenMCEtaPhi[mcIndex]->Fill(etaMC, phiMC, GetEventWeight());
      
      Bool_t inacceptance = kTRUE;
      // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
      if( IsFiducialCutOn() && !GetFiducialCut()->IsInFiducialCut(fPrimaryMomentum.Eta(),fPrimaryMomentum.Phi(),GetCalorimeter()) )
        inacceptance = kFALSE ;
      
      if(IsRealCaloAcceptanceOn()) // defined on base class
      {
        if(GetReader()->ReadStack()          &&
           !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), primStack)) inacceptance = kFALSE ;
        if(GetReader()->ReadAODMCParticles() &&
           !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), primAOD  )) inacceptance = kFALSE ;
      }
      
      if(!inacceptance) continue;
      
      fhGenMCAccE [mcIndex]->Fill( eMC, GetEventWeight());
      fhGenMCAccPt[mcIndex]->Fill(ptMC, GetEventWeight());
      if(eMC > 0.5) fhGenMCAccEtaPhi[mcIndex]->Fill(etaMC, phiMC, GetEventWeight());
    }
  }
}

//_________________________________________________________________________________
/// Check cluster weights, check the effect of different w0 parameter on shower shape
/// Check effect of time and energy cuts at cell level on the shower shape
//_________________________________________________________________________________
void AliAnaCalorimeterQA::WeightHistograms(AliVCluster *clus, AliVCaloCells* cells)
{
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  Float_t  ampMax = 0;
  Int_t    bc     = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  Float_t  enOrg  = clus->E();

  // Do study when there are enough cells in cluster
  if(clus->GetNCells() < 3) return ;
  
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) 
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp, GetCalorimeter(), id);
    
    energy    += amp;
    
    if ( amp > ampMax )
      ampMax = amp;
    
  } // energy loop       
  
  if ( energy <=0 )
  {
    AliWarning(Form("Wrong calculated energy %f",energy));
    return;
  }
  
  // Remove non linearity correction
  clus->SetE(energy);
  
  fhEMaxCellClusterRatio   ->Fill(energy, ampMax/energy            , GetEventWeight());
  fhEMaxCellClusterLogRatio->Fill(energy, TMath::Log(ampMax/energy), GetEventWeight());
  
  // Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++) 
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    // Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp, GetCalorimeter(), id);
    
    fhECellClusterRatio   ->Fill(energy, amp/energy            , GetEventWeight());
    fhECellClusterLogRatio->Fill(energy, TMath::Log(amp/energy), GetEventWeight());
  }        
  
  //
  // Recalculate shower shape for different W0 and cell cuts
  //
  if ( GetCalorimeter() == kEMCAL )
  {
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    Float_t w0org = GetCaloUtils()->GetEMCALRecoUtils()->GetW0();
    
    Int_t tagMC = -1;
    if(IsDataMC() && clus->GetNLabels() > 0)
    {
      Int_t tag = GetMCAnalysisUtils()->CheckOrigin(clus->GetLabels(),clus->GetNLabels(), GetReader(),GetCalorimeter());
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)   &&
         !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0)      &&
         !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta)      &&
         !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
        tagMC = 0;
      } // Pure Photon
      else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCElectron) &&
              !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
        tagMC = 1;
      } // Electron
      else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCConversion)        ){
        tagMC = 2;
      } // Conversion
      else if( GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0) ){
        tagMC = 3;
      }// Pi0
      else if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEta) &&
              !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton) ){
        tagMC = 4;
      } // Hadron
    } // Is MC
    
    // To be done properly with setters and data members ... 
    Float_t cellEmin [] = {0.05,0.1,0.15,0.2}; 
    Float_t cellTmin [] = {50.,100.,1000000.}; 

    for(Int_t iw = 0; iw < 12; iw++)
    {
      for(Int_t iEmin = 0; iEmin < 4; iEmin++)
      {
        for(Int_t iTmin = 0; iTmin < 3; iTmin++)
        {
          GetCaloUtils()->GetEMCALRecoUtils()->SetW0(4+iw*0.05);
          
          Float_t newEnergy = 0; 
          Float_t l0   = 0, l1   = 0;
          Float_t disp = 0, dEta = 0, dPhi    = 0;
          Float_t sEta = 0, sPhi = 0, sEtaPhi = 0;
          
          GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParametersWithCellCuts(GetEMCALGeometry(),  
                                                                                                   GetReader()->GetEMCALCells(), clus, 
                                                                                                   cellEmin[iEmin], cellTmin[iTmin], bc, 
                                                                                                   newEnergy, l0, l1, disp, dEta, dPhi, 
                                                                                                   sEta, sPhi, sEtaPhi);

          

          fhLambda0ForW0AndCellCuts[iw][iEmin][iTmin]->Fill(enOrg, l0, GetEventWeight());
//        fhLambda1ForW0AndCellCuts[iw][iEmin][iTmin]->Fill(enOrg, l0, GetEventWeight());
          
          if(TMath::Abs(fClusterMomentum.Eta()) < 0.15)
            fhLambda0ForW0AndCellCutsEta0[iw][iEmin][iTmin]->Fill(enOrg, l0, GetEventWeight());
          
        } // E cell loop
      } // Time loop
      
      if(IsDataMC() && tagMC >= 0)
      {
        fhLambda0ForW0MC[iw][tagMC]->Fill(energy, clus->GetM02(), GetEventWeight());
//      fhLambda1ForW0MC[iw][tagMC]->Fill(energy, clus->GetM20(), GetEventWeight());
      }
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    
    GetCaloUtils()->GetEMCALRecoUtils()->SetW0(w0org);
    
  } // EMCAL
  
  clus->SetE(enOrg);
}


/// Comment out, done in AliEMCALRecoUtils, except few details, keep for the moment in this manner
////___________________________________________________________________________________________________________________
///// Calculates new center of gravity in the local EMCAL-module coordinates
///// and tranfers into global ALICE coordinates.
///// Calculates Dispersion and main axis.
////___________________________________________________________________________________________________________________
//void AliAnaCalorimeterQA::RecalculateClusterShowerShapeParametersWithCellCut(const AliEMCALGeometry * geom,
//                                                                             AliVCaloCells* cells, AliVCluster * cluster,
//                                                                             Float_t eCellMin, Float_t & newEnergy,
//                                                                             Float_t & l0,   Float_t & l1,
//                                                                             Float_t & disp, Float_t & dEta, Float_t & dPhi,
//                                                                             Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi)
//{
//  if(!cluster)
//  {
//    AliWarning("Cluster pointer null!");
//    return;
//  }
//  
//  Double_t eCell       = 0.;
//  Float_t  fraction    = 1.;
//  Float_t  recalFactor = 1.;
//  
//  Int_t    iSupMod = -1;
//  Int_t    iTower  = -1;
//  Int_t    iIphi   = -1;
//  Int_t    iIeta   = -1;
//  Int_t    iphi    = -1;
//  Int_t    ieta    = -1;
//  Double_t etai    = -1.;
//  Double_t phii    = -1.;
//  
//  Int_t    nstat   = 0 ;
//  Float_t  wtot    = 0.;
//  Double_t w       = 0.;
//  Double_t etaMean = 0.;
//  Double_t phiMean = 0.;
//  
//  newEnergy = 0;
//  
//  Bool_t  shared = GetCaloUtils()-> IsClusterSharedByTwoSuperModules(geom,cluster);
//  
//  Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(cluster, cells);
//  
////  Float_t simuTotWeight = 0;
////  if(GetCaloUtils()->IsMCECellClusFracCorrectionOn())
////  {
////    simuTotWeight =  GetCaloUtils()->RecalibrateClusterEnergyWeightCell(cluster, cells,energy);
////    simuTotWeight/= energy;
////  }
//  
//  // Loop on cells, get weighted parameters
//  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
//  {
//    // Get from the absid the supermodule, tower and eta/phi numbers
//    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
//    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
//    
//    // Get the cell energy, if recalibration is on, apply factors
//    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
//    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
//    
//    if(GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
//    {
//      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
//    }
//    
//    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
//    
//    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
//    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
//    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
//    
//    if(energy > 0 && eCell > eCellMin)
//    {
//      if(GetCaloUtils()->IsMCECellClusFracCorrectionOn())
//        eCell*=GetCaloUtils()->GetMCECellClusFracCorrection(eCell,energy);// /simuTotWeight;
//      
//      newEnergy+=eCell;
//      
//      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);
//      
//      // correct weight, ONLY in simulation
//      //w *= (fWSimu[0] - fWSimu[1] * w );
//      
//      etai=(Double_t)ieta;
//      phii=(Double_t)iphi;
//      
//      if(w > 0.0)
//      {
//        wtot += w ;
//        nstat++;
//        //Shower shape
//        sEta     += w * etai * etai ;
//        etaMean  += w * etai ;
//        sPhi     += w * phii * phii ;
//        phiMean  += w * phii ;
//        sEtaPhi  += w * etai * phii ;
//      }
//    }
//    else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f", eCell, energy));
//  }//cell loop
//  
//  // Normalize to the weight
//  if (wtot > 0)
//  {
//    etaMean /= wtot ;
//    phiMean /= wtot ;
//  }
//  //else
//  //  AliError(Form("Wrong weight %f", wtot));
//  
//  // Calculate dispersion
//  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
//  {
//    // Get from the absid the supermodule, tower and eta/phi numbers
//    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
//    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
//    
//    // Get the cell energy, if recalibration is on, apply factors
//    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
//    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
//    if (GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
//    {
//      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
//    }
//    
//    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
//    
//    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
//    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
//    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
//    
//    if(energy > 0 && eCell > eCellMin)
//    {
//      if(GetCaloUtils()->IsMCECellClusFracCorrectionOn())
//        eCell*=GetCaloUtils()->GetMCECellClusFracCorrection(eCell,energy); // /simuTotWeight;
//      
//      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);
//      
//      //correct weight, ONLY in simulation
//      //w *= (fWSimu[0] - fWSimu[1] * w );
//      
//      etai=(Double_t)ieta;
//      phii=(Double_t)iphi;
//      if(w > 0.0)
//      {
//        disp +=  w *((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
//        dEta +=  w * (etai-etaMean)*(etai-etaMean) ;
//        dPhi +=  w * (phii-phiMean)*(phii-phiMean) ;
//      }
//    }
//    else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f", eCell, energy));
//  } // cell loop
//  
//  // Normalize to the weigth and set shower shape parameters
//  if (wtot > 0 && nstat > 1)
//  {
//    disp    /= wtot ;
//    dEta    /= wtot ;
//    dPhi    /= wtot ;
//    sEta    /= wtot ;
//    sPhi    /= wtot ;
//    sEtaPhi /= wtot ;
//    
//    sEta    -= etaMean * etaMean ;
//    sPhi    -= phiMean * phiMean ;
//    sEtaPhi -= etaMean * phiMean ;
//    
//    l0 = (0.5 * (sEta + sPhi) + TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
//    l1 = (0.5 * (sEta + sPhi) - TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
//  }
//  else
//  {
//    l0   = 0. ;
//    l1   = 0. ;
//    dEta = 0. ; dPhi = 0. ; disp    = 0. ;
//    sEta = 0. ; sPhi = 0. ; sEtaPhi = 0. ;
//  }
//}



