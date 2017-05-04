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
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliMCAnalysisUtils.h"
#include "TCustomBinning.h"

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
fFillAllTMHisto(kTRUE),                fFillClusterMaxCellHisto(kFALSE),   
fFillAllPi0Histo(kTRUE),               fFillInvMassOpenAngle(kFALSE),               
fFillPi0PairDiffTime(kFALSE),          fFillInvMassInEMCALWithPHOSDCalAcc(kFALSE), 
fFillEBinAcceptanceHisto(kFALSE),      fFillAllClusterHistograms(kTRUE),
fFillAllCellHistograms(kTRUE),         fFillAllCellAbsIdHistograms(kTRUE),
fCorrelate(kTRUE),                     fStudyBadClusters(kFALSE),               
fStudyClustersAsymmetry(kFALSE),       fStudyExotic(kFALSE),
fStudyWeight(kFALSE),                  fStudyTCardCorrelation(kFALSE), 
fStudyM02Dependence (kFALSE),

// Parameters and cuts
fTimeCutMin(-10000),                   fTimeCutMax(10000),
fCellAmpMin(0),                        
fEMCALCellAmpMin(0),                   fPHOSCellAmpMin(0),     
fEMCALClusterM02Min(0),
fEMCALClusterNCellMin(0),              fPHOSClusterNCellMin(0),  
fNEBinCuts(0),

// Invariant mass
fInvMassMinECut(0),                    fInvMassMaxECut(0),                    
fInvMassMinM02Cut(0),                  fInvMassMaxM02Cut(0),                    
fInvMassMaxOpenAngle(0),               fInvMassMaxTimeDifference(0),

//// Exotic
//fExoNECrossCuts(0),                    fExoECrossCuts(),
//fExoNDTimeCuts(0),                     fExoDTimeCuts(),    

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
fhNCellsPerCluster(0),                 fhNCellsPerClusterNoCut(0),             
fhNClusters(0),

// Timing
fhClusterTimeEnergy(0),                fhCellTimeSpreadRespectToCellMax(0),  
fhCellIdCellLargeTimeSpread(0),        fhClusterPairDiffTimeE(0),              fhClusterPairDiffTimeESameMod(0),
fhClusterMaxCellCloseCellRatio(0),     fhClusterMaxCellCloseCellDiff(0), 
fhClusterMaxCellDiff(0),               fhClusterMaxCellDiffNoCut(0), 
//fhClusterMaxCellDiffAverageTime(0),    fhClusterMaxCellDiffWeightedTime(0),    
fhClusterMaxCellECross(0),
// fhDispersion(0),
fhLambda0(0),                          fhLambda1(0),                           fhNLocMax(0),  

fhEnergyTMEtaResidual1Cell(0),         fhEnergyTMPhiResidual1Cell(0),
fhColRowExoticHighE1CellPosTime(0),    fhColRowExoticHighE1CellNegTime(0),     fhColRowExoticHighE1CellNulTime(0),
fhEnergyTMEtaResidualExotic(0),        fhEnergyTMPhiResidualExotic(0),
fhColRowExoticHighEPosTime(0),         fhColRowExoticHighENegTime(0),          fhColRowExoticHighENulTime(0),
fhColRowHighEPosTime(0),               fhColRowHighENegTime(0),                fhColRowHighENulTime(0),
fhEnergyTMEtaResidualTCardCorrNoSelection1Cell(0),  fhEnergyTMPhiResidualTCardCorrNoSelection1Cell(0),
fhEnergyTMEtaResidualTCardCorrNoSelectionExotic(0), fhEnergyTMPhiResidualTCardCorrNoSelectionExotic(0),


// bad clusters
fhBadClusterEnergy(0),                 fhBadClusterTimeEnergy(0),              fhBadClusterEtaPhi(0),            
fhBadClusterPairDiffTimeE(0),          fhBadCellTimeSpreadRespectToCellMax(0), 
fhBadClusterMaxCellCloseCellRatio(0),  fhBadClusterMaxCellCloseCellDiff(0),    fhBadClusterMaxCellDiff(0),
//fhBadClusterMaxCellDiffAverageTime(0), fhBadClusterMaxCellDiffWeightedTime(0),
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
fhTimeId(0),                           fhTimeL1UnCorrId(0),                    fhTimeAmp(0),
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

//fhExoL0ECross(0),                      fhExoL1ECross(0),

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
fhTrackMatchedDEtaPosMod(0),           fhTrackMatchedDPhiPosMod(0),
fhClusterTimeEnergyM02(0),             fhCellTimeSpreadRespectToCellMaxM02(0), 
fhClusterMaxCellCloseCellRatioM02(0),  fhClusterMaxCellCloseCellDiffM02(0),  
fhClusterMaxCellDiffM02(0),            fhClusterMaxCellECrossM02(0),           fhNCellsPerClusterM02(0) 
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
  
//  // Exotic
//  for (Int_t ie = 0; ie < 10 ; ie++) 
//  {
//    fhExoDTime[ie] = 0;
//    for (Int_t idt = 0; idt < 5 ; idt++) 
//    {
//      fhExoNCell    [ie][idt] = 0;
//      fhExoL0       [ie][idt] = 0;
//      fhExoL1       [ie][idt] = 0;
//      fhExoECross   [ie][idt] = 0;
//      fhExoTime     [ie][idt] = 0;
//      fhExoL0NCell  [ie][idt] = 0;
//      fhExoL1NCell  [ie][idt] = 0;
//    } 
//  }
  
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
  
  for(Int_t i = 0; i < 14; i++)
  {
    fhEBinClusterEtaPhi[i] = 0 ;
    fhEBinClusterColRow[i] = 0 ;
    fhEBinCellColRow   [i] = 0 ; 
  }
  
  for(Int_t bc = 0; bc < 4; bc++) fhTimePerSMPerBC[bc] = 0 ;
  
  // TCard correl studies
  for(Int_t tm = 0; tm < 2;  tm++)
  { 
    fhEnergyTime1Cell       [tm] = 0;
    fhColRowExoticLowE1Cell [tm] = 0;
    fhColRowExoticHighE1Cell[tm] = 0;
    
    fhEnergyTimeExotic [tm] = 0;
    fhColRowExoticLowE [tm] = 0 ;       
    fhColRowExoticHighE[tm] = 0 ;       
    fhColRowExotic2ndCellDiffLowE [tm] = 0 ;
    fhColRowExotic2ndCellDiffHighE[tm] = 0 ;
    fhColRowExotic2ndCellSameLowE [tm] = 0 ;
    fhColRowExotic2ndCellSameHighE[tm] = 0 ;

    fhEnergyTimeTCardCorrNoSelection1Cell [tm] = 0;
    fhEnergyTimeTCardCorrNoSelectionExotic[tm] = 0;
    
    fhColRowTCardCorrNoSelectionExoticLowE [tm] = 0 ;       
    fhColRowTCardCorrNoSelectionExoticHighE[tm] = 0 ;       

    fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE [tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE [tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[tm] = 0 ;

    fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE [tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE [tm] = 0 ;
    fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[tm] = 0 ;
    
    fhColRowTCardCorrNoSelectionLowE [tm] = 0 ;       
    fhColRowTCardCorrNoSelectionHighE[tm] = 0 ;  
    
    fhLambda0TCardCorrNoSelection[tm] = 0 ;  
    fhLambda1TCardCorrNoSelection[tm] = 0 ;      
    fhLambda0NLM1TCardCorrNoSelection[tm] = 0 ;  
    fhLambda1NLM1TCardCorrNoSelection[tm] = 0 ;      
    fhLambda0NLM2TCardCorrNoSelection[tm] = 0 ;  
    fhLambda1NLM2TCardCorrNoSelection[tm] = 0 ;  
    fhLambdaRTCardCorrNoSelection[tm] = 0 ;  
    fhNLocMaxTCardCorrNoSelection[tm] = 0 ;  
    
    fhEMaxRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhEMaxRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhEMaxRatNLM3TCardCorrNoSelection[tm] = 0 ;     
    fhE2ndRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndRatNLM3TCardCorrNoSelection[tm] = 0 ;      
    fhE2ndEMaxRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndEMaxRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndEMaxRatNLM3TCardCorrNoSelection[tm] = 0 ;

    fhE2ndSameRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndSameRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndSameRatNLM3TCardCorrNoSelection[tm] = 0 ;      
    fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[tm] = 0 ;
    
    fhECellClusRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhECellClusRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhECellClusRatNLM3TCardCorrNoSelection[tm] = 0 ;         
    fhLogECellNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhLogECellNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhLogECellNLM3TCardCorrNoSelection[tm] = 0 ;     
    fhECellWeightNLM1TCardCorrNoSelection [tm] = 0 ;  
    fhECellWeightNLM2TCardCorrNoSelection [tm] = 0 ;  
    fhECellWeightNLM3TCardCorrNoSelection [tm] = 0 ;      

    fhECellSameClusRatNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhECellSameClusRatNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhECellSameClusRatNLM3TCardCorrNoSelection[tm] = 0 ;         
    fhLogECellSameNLM1TCardCorrNoSelection[tm] = 0 ;  
    fhLogECellSameNLM2TCardCorrNoSelection[tm] = 0 ;  
    fhLogECellSameNLM3TCardCorrNoSelection[tm] = 0 ;     
    fhECellSameWeightNLM1TCardCorrNoSelection [tm] = 0 ;  
    fhECellSameWeightNLM2TCardCorrNoSelection [tm] = 0 ;  
    fhECellSameWeightNLM3TCardCorrNoSelection [tm] = 0 ;      

    
    fhNCellsTCardCorrNoSelection [tm] = 0 ;        
    fhNCellsTCardCorrWithWeightNoSelection     [tm] = 0 ;  
    fhNCellsTCardCorrRatioWithWeightNoSelection[tm] = 0 ; 
    fhExoticTCardCorrNoSelection [tm] = 0 ;        
    fhNCellsTCardSameAndDiffFraction      [tm] = 0;
    fhNCellsTCardSameAndDiffFractionExotic[tm] = 0;

    fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[tm] = 0;
    fhSameRowDiffColAndTCardCellsTimeDiffClusterE  [tm] = 0;
    fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[tm] = 0;
    fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE  [tm] = 0;    
    fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[tm] = 0;
    fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo  [tm] = 0;
    fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[tm] = 0;
    fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo  [tm] = 0;

    for(Int_t i = 0; i < 6; i++)
    {
      for(Int_t j = 0; j < 6; j++)
      {
        fhLambda0TCardCorrelNCell[i][j][tm] = 0 ;  
        fhLambda1TCardCorrelNCell[i][j][tm] = 0 ;          
        fhLambda0NLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhLambda1NLM1TCardCorrelNCell[i][j][tm] = 0 ;          
        fhLambda0NLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhLambda1NLM2TCardCorrelNCell[i][j][tm] = 0 ;  
//      fhLambdaRTCardCorrelNCell[i][j][tm] = 0 ;  
        fhNLocMaxTCardCorrelNCell[i][j][tm] = 0 ;  
 
        fhEMaxRatNLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhEMaxRatNLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhEMaxRatNLM3TCardCorrelNCell[i][j][tm] = 0 ;          
        fhE2ndRatNLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhE2ndRatNLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhE2ndRatNLM3TCardCorrelNCell[i][j][tm] = 0 ;          
        fhE2ndEMaxRatNLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhE2ndEMaxRatNLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhE2ndEMaxRatNLM3TCardCorrelNCell[i][j][tm] = 0 ;  
        
        fhECellClusRatNLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhECellClusRatNLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhECellClusRatNLM3TCardCorrelNCell[i][j][tm] = 0 ;              
        fhLogECellNLM1TCardCorrelNCell[i][j][tm] = 0 ;  
        fhLogECellNLM2TCardCorrelNCell[i][j][tm] = 0 ;  
        fhLogECellNLM3TCardCorrelNCell[i][j][tm] = 0 ;          
        fhECellWeightNLM1TCardCorrelNCell [i][j][tm] = 0 ;  
        fhECellWeightNLM2TCardCorrelNCell [i][j][tm] = 0 ;  
        fhECellWeightNLM3TCardCorrelNCell [i][j][tm] = 0 ;     
        
        fhMassEClusTCardCorrelNCell[i][j][tm] = 0 ;  
//      fhMassEPairTCardCorrelNCell[i][j][tm] = 0 ;  
        fhExoticTCardCorrelNCell   [i][j][tm] = 0 ;  
        fhTimeDiffTCardCorrelNCell [i][j][tm] = 0 ;  
        fhTimeDiffExoTCardCorrelNCell[i][j][tm] = 0 ;  
        fhColRowTCardCorrelNCellLowE [i][j][tm] = 0 ; 
        fhColRowTCardCorrelNCellHighE[i][j][tm] = 0 ; 
      }
      
//      fhLambda0TCardCorrelN[i][tm] = 0 ;  
//      fhNCellsTCardCorrelN [i][tm] = 0 ;   
//      fhExoticTCardCorrelN [i][tm] = 0 ;   
//      fhColRowTCardCorrelNLowE [i][tm] = 0 ; 
//      fhColRowTCardCorrelNHighE[i][tm] = 0 ; 
//      
//      fhLambda0TCardCorrelNAllSameTCard[i][tm] = 0 ;  
//      fhNCellsTCardCorrelNAllSameTCard [i][tm] = 0 ;   
//      fhExoticTCardCorrelNAllSameTCard [i][tm] = 0 ;   
//      fhColRowTCardCorrelNAllSameTCardLowE [i][tm] = 0 ; 
//      fhColRowTCardCorrelNAllSameTCardHighE[i][tm] = 0 ;
//      
//      fhLambda0TCardCorrelNExotic[i][tm] = 0 ;  
//      fhNCellsTCardCorrelNExotic [i][tm] = 0 ; 
//      fhColRowTCardCorrelNLowEExotic [i][tm] = 0 ; 
//      fhColRowTCardCorrelNHighEExotic[i][tm] = 0 ; 
//      
//      fhLambda0TCardCorrelNAllSameTCardExotic[i][tm] = 0 ;  
//      fhNCellsTCardCorrelNAllSameTCardExotic [i][tm] = 0 ; 
//      fhColRowTCardCorrelNAllSameTCardLowEExotic [i][tm] = 0 ; 
//      fhColRowTCardCorrelNAllSameTCardHighEExotic[i][tm] = 0 ; 
//      fhLambda0TCardCorrelNearRow[i][tm] = 0 ; 
//      fhNCellsTCardCorrelNearRow [i][tm] = 0 ; 
    }
    
//    for(Int_t i = 0; i < 4; i++)
//    {
//      fhLambda0TCardCorrel2ndMax[i][tm] = 0 ; 
//      fhNCellsTCardCorrel2ndMax [i][tm] = 0 ; 
//      fhLambda0TCardCorrelExotic[i][tm] = 0 ;  
//      fhNCellsTCardCorrelExotic [i][tm] = 0 ; 
//    }
    
    for(Int_t i = 0; i < 14; i++)
    {
      fhLambda0Exoticity[i][tm] = 0;
      fhLambda1Exoticity[i][tm] = 0;
//    fhLambdaRExoticity[i][tm] = 0;
      fhNCellsExoticity [i][tm] = 0;
      fhTimeExoticity   [i][tm] = 0;
      fhLambda0Lambda1  [i][tm] = 0;

//    fhLambda0ExoticityAllSameTCard[i][tm] = 0;
//    fhLambda1ExoticityAllSameTCard[i][tm] = 0;
//    fhLambdaRExoticityAllSameTCard[i][tm] = 0;
//    fhNCellsExoticityAllSameTCard [i][tm] = 0;
//    fhLambda0Lambda1AllSameTCard  [i][tm] = 0;

      fhNCellsTCardSameAndDiff      [i][tm] = 0;
      fhNCellsTCardSameAndDiffExotic[i][tm] = 0;
    }
    
    for(Int_t j = 0; j < 6; j++)
    { 
      for(Int_t k = 0; k < 6; k++)
      {
        fhLambda0ExoticityPerNCell[j][k][tm] = 0;
        fhLambda1ExoticityPerNCell[j][k][tm] = 0;
//      fhLambdaRExoticityPerNCell[j][k][tm] = 0;
      }
    }

//    for(Int_t i = 0; i < 7; i++)
//    {
//      fhLambda0TCardCorrel[i][tm] = 0 ; 
//      fhNCellsTCardCorrel [i][tm] = 0 ; 
//      fhExoticTCardCorrel [i][tm] = 0 ; 
//      
//      fhLambda0TCardCorrelOtherTCard[i][tm] = 0 ; 
//      fhNCellsTCardCorrelOtherTCard [i][tm] = 0 ; 
//      fhExoticTCardCorrelOtherTCard [i][tm] = 0 ; 
//      fhColRowTCardCorrelOtherTCardLowE [i][tm] = 0 ; 
//      fhColRowTCardCorrelOtherTCardHighE[i][tm] = 0 ; 
//    }
    
    for(Int_t i = 0; i < 12; i++)
    {
      fhTCardCorrECellMaxDiff[i][tm]=0;                 
      fhTCardCorrEClusterDiff[i][tm]=0;                  
//    fhTCardCorrECellMaxRat [i][tm]=0;                 
//    fhTCardCorrEClusterRat [i][tm]=0;                  
      fhTCardCorrTCellMaxDiff[i][tm]=0;                 
      
      fhTCardCorrECellMaxDiffExo[i][tm]=0;               
      fhTCardCorrEClusterDiffExo[i][tm]=0;               
//    fhTCardCorrECellMaxRatExo [i][tm]=0;             
//    fhTCardCorrEClusterRatExo [i][tm]=0;            
      fhTCardCorrTCellMaxDiffExo[i][tm]=0;         
    }
  }
  
  for(Int_t i = 0; i < fNEBinCuts; i++)
  {
    fhTMPhiResidualExoticity[i]  = 0;
    fhTMEtaResidualExoticity[i]  = 0;
    fhTMPhiResidualExoticityLooseCut[i]  = 0;
    fhTMEtaResidualExoticityLooseCut[i]  = 0;
//  fhTMPhiResidualExoticityAllSameTCard[i]  = 0;
//  fhTMEtaResidualExoticityAllSameTCard[i]  = 0;
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
  
  if(fFillClusterMaxCellHisto)
  {
    fhBadClusterMaxCellDiff  ->Fill(energy, maxCellFraction, GetEventWeight());
    fhBadClusterMaxCellECross->Fill(energy, eCrossFrac     , GetEventWeight());
  }
  
  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi += TMath::TwoPi();
  
  if(energy > 0.5) fhBadClusterEtaPhi->Fill(fClusterMomentum.Eta(), phi, GetEventWeight());

  fhBadClusterLambda0->Fill(energy, clus->GetM02());
  fhBadClusterLambda1->Fill(energy, clus->GetM20());
  
  if(fStudyClustersAsymmetry) ClusterAsymmetryHistograms(clus,absIdMax,kFALSE);
  
  if(fFillPi0PairDiffTime)
  {
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
  }
  
//  // Max cell compared to other cells in cluster
//  if(fFillAllCellTimeHisto && fFillClusterMaxCellHisto) 
//  {
//    // Get some time averages
//    Double_t timeAverages[2] = {0.,0.};
//    CalculateAverageTime(clus, cells, timeAverages);
//
//    fhBadClusterMaxCellDiffAverageTime ->Fill(clus->E(), tmax-timeAverages[0], GetEventWeight());
//    fhBadClusterMaxCellDiffWeightedTime->Fill(clus->E(), tmax-timeAverages[1], GetEventWeight());
//  }           
  
  if(fFillClusterMaxCellHisto)
  {
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
  if(!fFillAllCellHistograms) return;
  
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
  Double_t timeL1UnCorr= 0.;
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
    
    if ( nModule < fFirstModule || nModule > fLastModule ) 
    {
      AliDebug(1,Form("Cell module out of range %d",nModule));
      continue ;
    }
    
    AliDebug(2,Form("\t module %d, column %d (%d), row %d (%d)", nModule,icolAbs,icol,irowAbs,irow));
    
    // Check if the cell is a bad channel
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
    timeL1UnCorr=time;
    // Correction of L1 phase if set
    GetCaloUtils()->RecalibrateCellTimeL1Phase(time,0 , nModule, GetReader()->GetInputEvent()->GetBunchCrossNumber());
    // Transform time to ns
    time *= 1.0e9;
    time-=fConstantTimeShift;
    timeL1UnCorr *= 1.0e9;
    timeL1UnCorr-=fConstantTimeShift;
    
    if(time < fTimeCutMin || time > fTimeCutMax)
    {
      AliDebug(1,Form("Remove cell with Time %f",time));
      continue;
    }
    
    // Remove exotic cells, defined only for EMCAL
    if(GetCalorimeter()==kEMCAL && 
       GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(id, cells, bc)) continue;
    
    Double_t binWidthCorrection=1;
    if(amp>=10)binWidthCorrection=1.0/4;
    if(amp>=20)binWidthCorrection=1.0/10;
    fhAmplitude  ->Fill(amp,          GetEventWeight());
    fhAmpMod     ->Fill(amp, nModule, GetEventWeight());
    fhAmpWeirdMod->Fill(amp, nModule, GetEventWeight());
    
    if(fFillAllCellAbsIdHistograms)
    {
      fhAmpId->Fill(amp, id     , GetEventWeight()*binWidthCorrection);
      
      if(!highG) fhAmpIdLowGain->Fill(amp, id, GetEventWeight());
    }
    
    if(fFillEBinAcceptanceHisto)
    {
      for(Int_t ie = 0; ie < fNEBinCuts; ie++)
      {
        if( amp >= fEBinCuts[ie] && amp < fEBinCuts[ie+1] )
        {            
          fhEBinCellColRow[ie]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
        }
      }
    }
    
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
        
        fhTime    ->Fill(time,       GetEventWeight());
        fhTimeAmp ->Fill(amp , time, GetEventWeight());
        
        if(fFillAllCellAbsIdHistograms) 
        {
          fhTimeId ->Fill(time, id  , GetEventWeight());
          
          if(GetCaloUtils()->IsL1PhaseInTimeRecalibrationOn()==1)
            fhTimeL1UnCorrId ->Fill(timeL1UnCorr, id  , GetEventWeight());
        }
        
        Int_t bc = (GetReader()->GetInputEvent()->GetBunchCrossNumber())%4;
        fhTimePerSMPerBC[bc]->Fill(time, nModule, GetEventWeight());
        
        fhGridCellsTime->Fill(icolAbs, irowAbs, time);
        if(!highG) fhGridCellsTimeLowGain->Fill(icolAbs, irowAbs, time);
        
        fhTimeMod->Fill(time, nModule, GetEventWeight());
        
        if(fFillAllCellAbsIdHistograms)
          fhTimeAmpPerRCU[nModule*fNRCU+iRCU]->Fill(amp, time, GetEventWeight());
        
        if(!highG)
        {
          if(fFillAllCellAbsIdHistograms) 
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
    
    
  } // Cell loop
  
  // Fill the cells after the cut on min amplitude and bad/exotic channels
  if( ncellsCut > 0 ) fhNCellsCutAmpMin->Fill(ncellsCut, GetEventWeight()) ;
  
  // Number of cells per module
  for(Int_t imod = 0; imod < fNModules; imod++ )
  {
    if ( imod < fFirstModule || imod > fLastModule ) continue ;

    AliDebug(1,Form("Module %d, calo %s, N cells %d, sum Amp %f", imod, GetCalorimeterString().Data(), nCellsInModule[imod], eCellsInModule[imod]));

    fhNCellsMod     ->Fill(nCellsInModule[imod], imod, GetEventWeight()) ;
    fhSumCellsAmpMod->Fill(eCellsInModule[imod], imod, GetEventWeight()) ;
    
    if(fFillAllCellAbsIdHistograms)
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

//___________________________________________________
/// Check effect of T-Card channels correlation on shower shape/energy
/// for EMCal
///
/// \param clus: cluster pointer
/// \param cells: list with all cells
/// \param matched: bool with loose track-matching info
/// \param absIdMax: id of highest energy cell in cluster
/// \param exoticity: cross energy fraction
///
//___________________________________________________
void AliAnaCalorimeterQA::ChannelCorrelationInTCard(AliVCluster* clus, AliVCaloCells* cells, 
                                                    Bool_t matched, Int_t absIdMax, Float_t exoticity) 
{
  // Get the col and row of the leading cluster cell

  Int_t icol = -1, irow = -1, iRCU = -1, icolAbs = -1, irowAbs = -1;
  GetModuleNumberCellIndexesAbsCaloMap(absIdMax,GetCalorimeter(), icol, irow, iRCU, icolAbs, irowAbs);

  Float_t energy = clus->E();
  Int_t   ebin = -1;
  for(Int_t ie = 0; ie < fNEBinCuts; ie++)
  {
    if( energy >= fEBinCuts[ie] && energy < fEBinCuts[ie+1] ) ebin = ie;
  }
  
  Int_t   ncells = clus->GetNCells();
  Double_t time  = clus->GetTOF()*1.e9;
  time-=fConstantTimeShift;
  Float_t deta   = clus->GetTrackDz();
  Float_t dphi   = clus->GetTrackDx();
  
  //if(clus->GetNCells()==1 && clus->E() > 4)printf("TCard E %f, NCells %d\n",energy,clus->GetNCells());

  if(fStudyExotic)
  {
    if ( ncells == 1 )
    {
      if(energy >= 5 && energy < 8)
        fhColRowExoticLowE1Cell [matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
      else if(energy >= 8)
      {
        fhColRowExoticHighE1Cell[matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
        
        if     ( time >  5) fhColRowExoticHighE1CellPosTime->Fill(icolAbs,irowAbs,GetEventWeight());
        else if( time < -5) fhColRowExoticHighE1CellNegTime->Fill(icolAbs,irowAbs,GetEventWeight());
        else                fhColRowExoticHighE1CellNulTime->Fill(icolAbs,irowAbs,GetEventWeight());
      }
      
      fhEnergyTime1Cell[matched]->Fill(energy,time,GetEventWeight());
      fhEnergyTMEtaResidual1Cell->Fill(energy,deta,GetEventWeight());
      fhEnergyTMPhiResidual1Cell->Fill(energy,dphi,GetEventWeight());
    }
    else if(exoticity > 0.97)
    {
      if(energy >= 5 && energy < 8)
        fhColRowExoticLowE [matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
      else if(energy >= 8)
      {
        fhColRowExoticHighE[matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
        
        if     ( time >  5) fhColRowExoticHighEPosTime->Fill(icolAbs,irowAbs,GetEventWeight());
        else if( time < -5) fhColRowExoticHighENegTime->Fill(icolAbs,irowAbs,GetEventWeight());
        else                fhColRowExoticHighENulTime->Fill(icolAbs,irowAbs,GetEventWeight());
      }
      
      fhEnergyTimeExotic[matched]->Fill(energy,time,GetEventWeight());
      fhEnergyTMEtaResidualExotic->Fill(energy,deta,GetEventWeight());
      fhEnergyTMPhiResidualExotic->Fill(energy,dphi,GetEventWeight());
      
      for (Int_t ipos = 0; ipos < ncells; ipos++) 
      {
        Int_t absId  = clus->GetCellsAbsId()[ipos];   
        
        Float_t  eCell = cells->GetCellAmplitude(absId);        
        GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
              
        // consider cells with enough energy weight and not the reference one
        Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
        
        if( absId == absIdMax || weight < 0.01 ) continue;
        
        Int_t rowDiff = -100, colDiff = -100;
        Bool_t sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
        
        // Get the col and row of the secondary cluster cell
        Int_t icol2 = -1, irow2 = -1, iRCU2 = -1, icolAbs2 = -1, irowAbs2 = -1;
        GetModuleNumberCellIndexesAbsCaloMap(absId,GetCalorimeter(), icol2, irow2, iRCU2, icolAbs2, irowAbs2);
        
        if ( !sameTCard )
        {
          if(energy >= 5 && energy < 8)
            fhColRowExotic2ndCellDiffLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
          else if(energy >= 8)
            fhColRowExotic2ndCellDiffHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        }
        else
        {
          if(energy >= 5 && energy < 8)
            fhColRowExotic2ndCellSameLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
          else if(energy >= 8)
            fhColRowExotic2ndCellSameHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        }
      }
      
    }
    else if ( energy > 8 )
    {
      if     ( time >  5) fhColRowHighEPosTime->Fill(icolAbs,irowAbs,GetEventWeight());
      else if( time < -5) fhColRowHighENegTime->Fill(icolAbs,irowAbs,GetEventWeight());
      else                fhColRowHighENulTime->Fill(icolAbs,irowAbs,GetEventWeight());
    }
    
    if(ebin > -1)
    {
      fhTMPhiResidualExoticityLooseCut[ebin]->Fill(exoticity,dphi);
      fhTMEtaResidualExoticityLooseCut[ebin]->Fill(exoticity,deta);
    }
  }
  else if ( energy > 8  && ncells > 1 && exoticity < 0.97)
  {
    if     ( time >  5) fhColRowHighEPosTime->Fill(icolAbs,irowAbs,GetEventWeight());
    else if( time < -5) fhColRowHighENegTime->Fill(icolAbs,irowAbs,GetEventWeight());
    else                fhColRowHighENulTime->Fill(icolAbs,irowAbs,GetEventWeight());
  }
  
  // Clean the sample
  
  // away from dead region
  if ( clus->GetDistanceToBadChannel() < 5 ) return ;  
  
  // in center of SM
  Int_t etaRegion = -1, phiRegion = -1;
  GetCaloUtils()->GetEMCALSubregion(clus,cells,etaRegion,phiRegion);
  // Region 0: center of SM ~0.18<|eta|<0.55
  if ( etaRegion !=0 ) return ;
  
  if(fStudyExotic)
  {
    if ( ncells == 1 )
    {
      fhEnergyTimeTCardCorrNoSelection1Cell[matched]->Fill(energy,time,GetEventWeight());
      fhEnergyTMEtaResidualTCardCorrNoSelection1Cell->Fill(energy,deta,GetEventWeight());
      fhEnergyTMPhiResidualTCardCorrNoSelection1Cell->Fill(energy,dphi,GetEventWeight());
    }
    else if(exoticity > 0.97)
    {
      fhEnergyTimeTCardCorrNoSelectionExotic[matched]->Fill(energy,time,GetEventWeight());
      fhEnergyTMEtaResidualTCardCorrNoSelectionExotic->Fill(energy,deta,GetEventWeight());
      fhEnergyTMPhiResidualTCardCorrNoSelectionExotic->Fill(energy,dphi,GetEventWeight());
    }
  }
  
  Float_t m02   = clus->GetM02();
  Float_t m20   = clus->GetM20();
  Float_t lamR  = 0;
  if ( m02 > 0.001 ) lamR = m20/m02;
  
  Int_t   absIdList[ncells]; 
  Float_t maxEList [ncells];
  Int_t nlm  = GetCaloUtils()->GetNumberOfLocalMaxima(clus, cells, absIdList, maxEList) ; 
//Int_t nlm  = GetCaloUtils()->GetNumberOfLocalMaxima(clus,cells);

  //
  // Correlation to max
  //
  Int_t nCellWithWeight = 1;
  Bool_t nearRow = kFALSE;
  Bool_t nearCol = kFALSE;
  Int_t nCorr    = 0;
  Int_t nCorrNo  = 0;
  Int_t sameCol  = 0;
  Int_t other    = 0;
  Int_t sameRow  = 0;
  Float_t  eCellMax = cells->GetCellAmplitude(absIdMax);  
  Double_t tCellMax = cells->GetCellTime(absIdMax);      
  //printf("Org E %2.2f, t %2.2f\n",eCellMax,tCellMax*1e9);
  GetCaloUtils()->RecalibrateCellAmplitude(eCellMax, GetCalorimeter(), absIdMax);
  GetCaloUtils()->RecalibrateCellTime(tCellMax, GetCalorimeter(), absIdMax, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
  //printf("New E %2.2f, t %2.2f\n",eCellMax,tCellMax*1e9);
  
  tCellMax *= 1.0e9;
  tCellMax-=fConstantTimeShift;
  
  // correlation not max cells
//  Int_t nCorr2   = 0;
//  Int_t sameCol2 = 0;
//  Int_t other2   = 0;
//  Int_t sameRow2 = 0;
  
  // Get second highest energy cell
//  Int_t absId2ndMax = -1;
  Float_t emax2nd = 0;
  Bool_t sameTCard2ndMax = kFALSE;
//Int_t  rowDiff2 = -100;
//Int_t  colDiff2 = -100;
   
  Float_t  eCellSameRowSameTCardNearCol = 0.;
  Float_t  eCellSameRowDiffTCardNearCol = 0.;
  Double_t tCellSameRowSameTCardNearCol = 0.;
  Double_t tCellSameRowDiffTCardNearCol = 0.;
  
  //printf("Cluster E %2.2f, ncells %d, absIdMax %d, eCell Max %2.2f\n", energy, ncells, absIdMax, cells->GetCellAmplitude(absIdMax));
  
  //
  // Loop on the cluster cells, define correlations
  //
  for (Int_t ipos = 0; ipos < ncells; ipos++) 
  {
    Int_t absId  = clus->GetCellsAbsId()[ipos];   
    
    Float_t  eCell = cells->GetCellAmplitude(absId);
    Double_t tCell = cells->GetCellTime(absId);      

    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
    GetCaloUtils()->RecalibrateCellTime(tCell, GetCalorimeter(), absId, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    tCell *= 1.0e9;
    tCell-=fConstantTimeShift;
    
    // consider cells with enough energy weight and not the reference one
    Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
    
    if( absId == absIdMax || weight < 0.01 ) continue;

    if     (nlm==1)  
    {
      fhECellClusRatNLM1TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM1TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM1TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
    }
    else if(nlm==2)  
    {
      fhECellClusRatNLM2TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM2TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM2TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
    }
    else         
    {
      fhECellClusRatNLM3TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM3TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM3TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
    }
    
    Int_t rowDiff = -100, colDiff = -100;
    Bool_t sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);
    
    if(sameTCard)
    {
      if     (nlm==1)  
      {
        fhECellSameClusRatNLM1TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
        fhLogECellSameNLM1TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
        fhECellSameWeightNLM1TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
      }
      else if(nlm==2)  
      {
        fhECellSameClusRatNLM2TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
        fhLogECellSameNLM2TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
        fhECellSameWeightNLM2TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
      }
      else         
      {
        fhECellSameClusRatNLM3TCardCorrNoSelection[matched]->Fill(energy, eCell/energy, GetEventWeight());
        fhLogECellSameNLM3TCardCorrNoSelection    [matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
        fhECellSameWeightNLM3TCardCorrNoSelection [matched]->Fill(energy, weight, GetEventWeight());
      }
    }
    
    //if(eCellMax < eCell) printf("Check: E max %f (id %d), E sec %f (id %d)\n",eCellMax,absIdMax, eCell,absId);
    
    nCellWithWeight++;
            
    //printf("\t cell %d, absId %d, E %2.2f, w %2.2f, tcard %d\n", ipos, absId, eCell, weight, sameTCard);

    Int_t indexType = -1;
    if ( sameTCard ) 
    {
      nCorr++;
      
      if(TMath::Abs(rowDiff) == 1) nearRow = kTRUE;
      if(TMath::Abs(colDiff) == 1) nearCol = kTRUE;
      
      if      ( rowDiff == 0  && colDiff != 0 ) 
      {
        if ( nearCol ) indexType = 6;         
        else           indexType = 7;

        sameRow++; 
        /*printf("\t \t E %2.2f, Same row, diff row %d, col %d\n",eCell,rowDiff,colDiff);*/
      }
      else if ( rowDiff != 0  && colDiff == 0 ) 
      {
        if ( nearRow ) indexType = 8;         
        else           indexType = 9;
        
        sameCol++; 
        /*printf("\t \t E %2.2f, Same col, diff row %d, col %d\n",eCell,rowDiff,colDiff);*/
      }
      else                                      
      {
        if ( nearRow && nearCol) indexType = 10;
        else                     indexType = 11;

        other++; 
        /*printf("\t \t E %2.2f, Diff row/col, diff row %d, col %d\n",eCell,rowDiff,colDiff);*/
      }
    }
    else
    {
      nCorrNo++;
      
      if      ( rowDiff == 0  && colDiff != 0 ) 
      {
        if ( nearCol ) indexType = 0;         
        else           indexType = 1;

      }
      else if ( rowDiff != 0  && colDiff == 0 ) 
      {
        if ( nearRow ) indexType = 2;
        else           indexType = 3;
      }
      else                                      
      {
        if ( nearCol && nearRow ) indexType = 4;
        else                      indexType = 5;
      }
   }  
    
    if ( rowDiff == 0 && TMath::Abs(colDiff) == 1 )
    {
      if(sameTCard) 
      {
        eCellSameRowSameTCardNearCol = eCell;
        tCellSameRowSameTCardNearCol = tCell;
      }
      else
      {
        eCellSameRowDiffTCardNearCol = eCell;
        tCellSameRowDiffTCardNearCol = tCell;
      }
    }
    
    if( indexType >=0 )
    {
      Float_t eCellDiff = eCellMax - eCell;
      Float_t eClusDiff = energy   - eCell;      
//    Float_t eCellRat  = eCell / eCellMax;
//    Float_t eClusRat  = eCell / energy  ;
      Float_t tCellDiff = tCellMax - tCell;

      fhTCardCorrECellMaxDiff[indexType][matched]->Fill(energy, eCellDiff, GetEventWeight());
      fhTCardCorrEClusterDiff[indexType][matched]->Fill(energy, eClusDiff, GetEventWeight());
//    fhTCardCorrECellMaxRat [indexType][matched]->Fill(energy, eCellRat , GetEventWeight());
//    fhTCardCorrEClusterRat [indexType][matched]->Fill(energy, eClusRat , GetEventWeight());
      fhTCardCorrTCellMaxDiff[indexType][matched]->Fill(energy, tCellDiff, GetEventWeight());
      
      if ( fStudyExotic && exoticity > 0.97 )
      {
        fhTCardCorrECellMaxDiffExo[indexType][matched]->Fill(energy, eCellDiff, GetEventWeight());
        fhTCardCorrEClusterDiffExo[indexType][matched]->Fill(energy, eClusDiff, GetEventWeight());
//      fhTCardCorrECellMaxRatExo [indexType][matched]->Fill(energy, eCellRat , GetEventWeight());
//      fhTCardCorrEClusterRatExo [indexType][matched]->Fill(energy, eClusRat , GetEventWeight());
        fhTCardCorrTCellMaxDiffExo[indexType][matched]->Fill(energy, tCellDiff, GetEventWeight());
      }      
    }
    
    if ( fStudyExotic && exoticity > 0.97 )
    {  
      // Get the col and row of the secondary cluster cell
      Int_t icol2 = -1, irow2 = -1, iRCU2 = -1, icolAbs2 = -1, irowAbs2 = -1;
      GetModuleNumberCellIndexesAbsCaloMap(absId,GetCalorimeter(), icol2, irow2, iRCU2, icolAbs2, irowAbs2);
      
      if ( !sameTCard )
      {
        if(energy >= 5 && energy < 8)
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        else if(energy >= 8)
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
      }
      else
      {
        if(energy >= 5 && energy < 8)
          fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        else if(energy >= 8)
          fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
      }
    }
    
    if ( eCell > emax2nd ) 
    {
      emax2nd     = eCell;
//    absId2ndMax = absId;
      if(sameTCard)
      {
        sameTCard2ndMax = kTRUE;
//        rowDiff2 = rowDiff;
//        colDiff2 = colDiff;
      }
      else
      {
        sameTCard2ndMax = kFALSE;
//        rowDiff2 = -100;
//        colDiff2 = -100;
      }
    }
    
//    //
//    // Other TCard correlations
//    //
//    if ( sameTCard ) continue;
//    
//    for (Int_t ipos2 = 0; ipos2 < ncells; ipos2++) 
//    {
//      Int_t absId2  = clus->GetCellsAbsId()[ipos2];   
//      
//      eCell = cells->GetCellAmplitude(absId2);      
//      // consider cells with enough energy weight and not the reference one
//      weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
//      
//      if( absId2 == absIdMax || absId2 == absId || weight < 0.01 ) continue;
//            
//      rowDiff = -100, colDiff = -100;
//      Bool_t sameTCard2 = GetCaloUtils()->IsAbsIDsFromTCard(absId,absId2,rowDiff,colDiff);
//      
//      if(sameTCard2)
//      {
//        nCorr2++;
//        if      ( rowDiff == 0  && colDiff != 0 ) sameRow2++; 
//        else if ( rowDiff != 0  && colDiff == 0 ) sameCol2++; 
//        else                                      other2++;
//      }
//      //printf("\t cell %d, absId %d, E %2.2f, w %2.2f, tcard %d\n", ipos, absId, eCell, weight, sameTCard);
//    } // second cluster cell lopp for secondary TCard correlations
  } // cluster cell loop
  
  Float_t ratioNcells = nCellWithWeight/(ncells*1.);
  fhNCellsTCardCorrRatioWithWeightNoSelection[matched]->Fill(energy, ratioNcells, GetEventWeight());
  //printf("E %2.2f, ncells %d, nCellWithWeight %d, ratio %2.2f\n",energy,ncells,nCellWithWeight,ratioNcells);
  
  // If only one relevant cell, it makes no sense to continue
  if ( nCellWithWeight <= 1 ) return;

  // It should not happen, unless very exotic clusters
  if ( m02 < 0.001 )
  {
    printf("AliAnaCalorimeterQA: M02 %f, M20 %f, E %2.3f, ncell %d, n with weight %d; max cell E %2.3f\n",
           m02,m20,energy,ncells,nCellWithWeight,eCellMax);
  }
  
  //printf("\t Same col %d, same row %d, diff other %d\n",sameCol,sameRow,other);
  //printf("\t Second cell: E %2.2f, absId %d, correl %d, rowDiff %d, rowCol %d\n",emax,absId2ndMax,sameTCard2,rowDiff2, colDiff2);
   
  //
  // Fill histograms for different cell correlation criteria
  //
  if(energy >= 5 && energy < 8)
    fhColRowTCardCorrNoSelectionLowE [matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
  else if(energy >= 8)
    fhColRowTCardCorrNoSelectionHighE[matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
  
  if ( fStudyExotic && exoticity > 0.97 )
  {
    if(energy >= 5 && energy < 8)
      fhColRowTCardCorrNoSelectionExoticLowE [matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
    else if(energy >= 8)
      fhColRowTCardCorrNoSelectionExoticHighE[matched]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
  }

  fhLambda0TCardCorrNoSelection[matched]->Fill(energy, m02, GetEventWeight());
  fhLambda1TCardCorrNoSelection[matched]->Fill(energy, m20, GetEventWeight());
  
  if     ( nlm == 1 )
  {
    fhLambda0NLM1TCardCorrNoSelection[matched]->Fill(energy, m02, GetEventWeight());
    fhLambda1NLM1TCardCorrNoSelection[matched]->Fill(energy, m20, GetEventWeight()); 
  }
  else if( nlm == 2 )
  {
    fhLambda0NLM2TCardCorrNoSelection[matched]->Fill(energy, m02, GetEventWeight());
    fhLambda1NLM2TCardCorrNoSelection[matched]->Fill(energy, m20, GetEventWeight()); 
  }
  
  fhLambdaRTCardCorrNoSelection[matched]->Fill(energy,lamR, GetEventWeight());
  fhNLocMaxTCardCorrNoSelection[matched]->Fill(energy, nlm, GetEventWeight());
  fhExoticTCardCorrNoSelection [matched]->Fill(energy, exoticity, GetEventWeight());
  
  if     (nlm==1)  
  {
    fhEMaxRatNLM1TCardCorrNoSelection    [matched]->Fill(energy, eCellMax/energy , GetEventWeight());
    fhE2ndRatNLM1TCardCorrNoSelection    [matched]->Fill(energy, emax2nd/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM1TCardCorrNoSelection[matched]->Fill(energy, emax2nd/eCellMax, GetEventWeight());
  }
  else if(nlm==2)  
  {
    fhEMaxRatNLM2TCardCorrNoSelection    [matched]->Fill(energy, eCellMax/energy     , GetEventWeight());
    fhE2ndRatNLM2TCardCorrNoSelection    [matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM2TCardCorrNoSelection[matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
  }
  else         
  {
    fhEMaxRatNLM3TCardCorrNoSelection    [matched]->Fill(energy, eCellMax/energy     , GetEventWeight());
    fhE2ndRatNLM3TCardCorrNoSelection    [matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM3TCardCorrNoSelection[matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
  }

  if(sameTCard2ndMax)
  {
    if     (nlm==1)  
    {
      fhE2ndSameRatNLM1TCardCorrNoSelection    [matched]->Fill(energy, emax2nd/energy  , GetEventWeight());
      fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[matched]->Fill(energy, emax2nd/eCellMax, GetEventWeight());
    }
    else if(nlm==2)  
    {
      fhE2ndSameRatNLM2TCardCorrNoSelection    [matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
      fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
    }
    else         
    {
      fhE2ndSameRatNLM3TCardCorrNoSelection    [matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
      fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
    }
  }
  
  fhNCellsTCardCorrNoSelection          [matched]->Fill(energy, ncells, GetEventWeight());
  fhNCellsTCardCorrWithWeightNoSelection[matched]->Fill(energy, nCellWithWeight, GetEventWeight());

  if(eCellSameRowSameTCardNearCol > 0 && eCellSameRowDiffTCardNearCol > 0)
  { 
    Float_t eDiff = eCellSameRowSameTCardNearCol - eCellSameRowDiffTCardNearCol ;
    Float_t tDiff = tCellSameRowSameTCardNearCol - tCellSameRowDiffTCardNearCol ;

    fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[matched]->Fill(energy  , eDiff, GetEventWeight());
    fhSameRowDiffColAndTCardCellsTimeDiffClusterE  [matched]->Fill(energy  , tDiff, GetEventWeight());
    fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[matched]->Fill(eCellMax, eDiff, GetEventWeight());
    fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE  [matched]->Fill(eCellMax, tDiff, GetEventWeight());
    
    if ( fStudyExotic && exoticity > 0.97 )
    {
      fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[matched]->Fill(energy  , eDiff, GetEventWeight());
      fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo  [matched]->Fill(energy  , tDiff, GetEventWeight());
      fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[matched]->Fill(eCellMax, eDiff, GetEventWeight());
      fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo  [matched]->Fill(eCellMax, tDiff, GetEventWeight());
    }
  }
    
  ///////////////////////////
  Int_t nCorrInd = nCorr;
  if(nCorr > 4) nCorrInd = 5; 
  
  Int_t nCorrNoInd = nCorrNo;
  if(nCorrNoInd > 4) nCorrNoInd = 5; 
  
//  fhLambda0TCardCorrelN[nCorrInd][matched]->Fill(energy, m02, GetEventWeight());
//  fhNCellsTCardCorrelN [nCorrInd][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//  fhExoticTCardCorrelN [nCorrInd][matched]->Fill(energy, exoticity, GetEventWeight());
//  
//  if ( energy >= 2 && energy < 8 ) 
//    fhColRowTCardCorrelNLowE [nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//  if ( energy >= 8 ) 
//    fhColRowTCardCorrelNHighE[nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());

  if     ( nlm == 1 ) 
  {
    fhLambda0NLM1TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m02, GetEventWeight());
    fhLambda1NLM1TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m20, GetEventWeight());
  }
  else if( nlm == 2 ) 
  {
    fhLambda0NLM2TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m02, GetEventWeight());
    fhLambda1NLM2TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m20, GetEventWeight());
  }
  
  fhLambda0TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m02, GetEventWeight());
  fhLambda1TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, m20, GetEventWeight());
//fhLambdaRTCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy,lamR, GetEventWeight());
  fhNLocMaxTCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, nlm, GetEventWeight());
  
  if(fStudyExotic)
    fhExoticTCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, exoticity, GetEventWeight());
  
  if     (nlm==1) 
  {
    fhEMaxRatNLM1TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, eCellMax/energy , GetEventWeight());
    fhE2ndRatNLM1TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, emax2nd/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM1TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, emax2nd/eCellMax, GetEventWeight());
  }
  else if(nlm==2)  
  {
    fhEMaxRatNLM2TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, eCellMax/energy, GetEventWeight());
    fhE2ndRatNLM2TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM2TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
  }
  else    
  {
    fhEMaxRatNLM3TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, eCellMax/energy, GetEventWeight());
    fhE2ndRatNLM3TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, maxEList[1]/energy  , GetEventWeight());
    fhE2ndEMaxRatNLM3TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, maxEList[1]/eCellMax, GetEventWeight());
  }
  
  // Time diff in cluster, depending nCells
  for (Int_t ipos = 0; ipos < ncells; ipos++) 
  {
    Int_t absId  = clus->GetCellsAbsId()[ipos];   
    
    Float_t  eCell = cells->GetCellAmplitude(absId);
    Double_t tCell = cells->GetCellTime(absId);      
    
    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
    GetCaloUtils()->RecalibrateCellTime(tCell, GetCalorimeter(), absId, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    tCell *= 1.0e9;
    tCell-=fConstantTimeShift;
    
    // consider cells with enough energy weight and not the reference one
    Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
    
    if( absId == absIdMax || weight < 0.01 ) continue;
    
    Float_t tDiffMaxSecondary = tCellMax - tCell; 
    //printf("Time max %f, second %f, diff %f\n",tCellMax,tCell,tDiffMaxSecondary);
    fhTimeDiffTCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, tDiffMaxSecondary, GetEventWeight());
    if ( fStudyExotic && exoticity > 0.97 ) 
    {
      fhTimeDiffExoTCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, tDiffMaxSecondary, GetEventWeight());
      // Get the col and row of the secondary cluster cell
      Int_t icol2 = -1, irow2 = -1, iRCU2 = -1, icolAbs2 = -1, irowAbs2 = -1;
      GetModuleNumberCellIndexesAbsCaloMap(absId,GetCalorimeter(), icol2, irow2, iRCU2, icolAbs2, irowAbs2);
      
      Int_t rowDiff = -100, colDiff = -100;
      Bool_t sameTCard = GetCaloUtils()->IsAbsIDsFromTCard(absIdMax,absId,rowDiff,colDiff);

      if ( !sameTCard )
      {
        if ( energy >= 5 && energy < 8 && nCorr == 0 )
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        else if ( energy >= 8 && nCorr == 0 )
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
      }
      else
      {
        if ( energy >= 5 && energy < 8 && nCorrNo == 0 )
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE [matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
        else if ( energy >= 8 && nCorrNo == 0 )
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[matched]->Fill(icolAbs2,irowAbs2,GetEventWeight()) ;
      }
    } // exotic
    
    if     (nlm==1)  
    {
      fhECellClusRatNLM1TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM1TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM1TCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, weight, GetEventWeight());
    }
    else if(nlm==2)  
    {
      fhECellClusRatNLM2TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM2TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM2TCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, weight, GetEventWeight());
    }
    else         
    {
      fhECellClusRatNLM3TCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, eCell/energy, GetEventWeight());
      fhLogECellNLM3TCardCorrelNCell    [nCorrInd][nCorrNoInd][matched]->Fill(energy, TMath::Log(eCell), GetEventWeight());
      fhECellWeightNLM3TCardCorrelNCell [nCorrInd][nCorrNoInd][matched]->Fill(energy, weight, GetEventWeight());
    }
  } //cell loop
  
  // Invariant mass for clusters looking like photons, depending number of cells
  if(m02 > fInvMassMinM02Cut && m02 < fInvMassMaxM02Cut)
  {
    for(Int_t jclus = 0 ; jclus < GetEMCALClusters()->GetEntriesFast() ; jclus++) 
    {
      AliVCluster* clus2 =  (AliVCluster*) GetEMCALClusters()->At(jclus);
      
      Float_t maxCellFraction = 0.;
      Int_t absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(cells, clus2, maxCellFraction);
      
      Double_t tof2 =  clus2->GetTOF()*1.e9;
      if(tof2>400) tof2-=fConstantTimeShift;
      
      Double_t diffTof = tCellMax-tof2;
      
      // Try to reduce background with a mild shower shape cut and no more 
      // than 1 local maximum in cluster and remove low energy clusters
      
      if(   absIdMax == absIdMax2   
         || !IsGoodCluster(absIdMax2, clus2->GetM02(), clus2->GetNCells(), cells) 
         || GetCaloUtils()->GetNumberOfLocalMaxima(clus2,cells) > 1 
         || clus2->GetM02() > fInvMassMaxM02Cut
         || clus2->GetM02() < fInvMassMinM02Cut
         || clus2->E() < fInvMassMinECut 
         || clus2->E() > fInvMassMaxECut  
         || TMath::Abs(diffTof) > fInvMassMaxTimeDifference
         ) continue;
      
      // Get cluster kinematics
      Double_t v[3] = {0,0,0}; //vertex ;
      clus2->GetMomentum(fClusterMomentum2,v);
      
      // Check only certain regions
      Bool_t in2 = kTRUE;
      if(IsFiducialCutOn()) in2 =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum2.Eta(),fClusterMomentum2.Phi(),GetCalorimeter()) ;
      if(!in2) continue;	
      
      //Float_t  pairE = (fClusterMomentum+fClusterMomentum2).E();
      
      // Opening angle cut, avoid combination of DCal and EMCal clusters
      Double_t angle  = fClusterMomentum.Angle(fClusterMomentum2.Vect());
      
      if( angle > fInvMassMaxOpenAngle ) continue;
      
      // Fill histograms
      Float_t mass   = (fClusterMomentum+fClusterMomentum2).M ();
      fhMassEClusTCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(energy, mass, GetEventWeight());
      //fhMassEPairTCardCorrelNCell[nCorrInd][nCorrNoInd][matched]->Fill(pairE , mass, GetEventWeight());
    }
  }
  
  if ( energy >= 5 && energy < 8) 
    fhColRowTCardCorrelNCellLowE [nCorrInd][nCorrNoInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
  else if ( energy >= 8 ) 
    fhColRowTCardCorrelNCellHighE[nCorrInd][nCorrNoInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
  
//  if(nCorrNo == 0)
//  {
//    fhLambda0TCardCorrelNAllSameTCard[nCorrInd][matched]->Fill(energy, m02, GetEventWeight());
//    fhNCellsTCardCorrelNAllSameTCard [nCorrInd][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    fhExoticTCardCorrelNAllSameTCard [nCorrInd][matched]->Fill(energy, exoticity, GetEventWeight());
//    if ( energy >= 2 && energy < 8 ) 
//      fhColRowTCardCorrelNAllSameTCardLowE[nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//    else if ( energy >= 8 ) 
//      fhColRowTCardCorrelNAllSameTCardHighE[nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//  }

  ///////////////////////////

  if ( fStudyExotic && exoticity > 0.97 )
  {
    if ( energy >= 5 && energy < 8) 
      fhColRowTCardCorrelNCellExoticLowE [nCorrInd][nCorrNoInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
      else if ( energy >= 8 ) 
        fhColRowTCardCorrelNCellExoticHighE[nCorrInd][nCorrNoInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());

//    fhLambda0TCardCorrelNExotic[nCorrInd][matched]->Fill(energy, m02, GetEventWeight());
//    fhNCellsTCardCorrelNExotic [nCorrInd][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    
//    if ( energy > 5 && energy <= 8 ) 
//      fhColRowTCardCorrelNLowEExotic [nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//    if ( energy > 8 ) 
//      fhColRowTCardCorrelNHighEExotic[nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//    
//    if(nCorrNo == 0)
//    {
//      fhLambda0TCardCorrelNAllSameTCardExotic[nCorrInd][matched]->Fill(energy, m02, GetEventWeight());
//      fhNCellsTCardCorrelNAllSameTCardExotic [nCorrInd][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//      
//      if ( energy > 2 && energy <=8 ) 
//        fhColRowTCardCorrelNAllSameTCardLowEExotic [nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//      if ( energy > 8 ) 
//        fhColRowTCardCorrelNAllSameTCardHighEExotic[nCorrInd][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//    }
//    
//    Int_t indexExo = -1;
//    if      (!nearRow &&  nearCol ) indexExo = 0;
//    else if (!nearRow && !nearCol ) indexExo = 1;
//    else if ( nearRow &&  nearCol ) indexExo = 2;
//    else if ( nearRow && !nearCol ) indexExo = 3;
//    
//    if(indexExo >= 0)
//    {
//      fhLambda0TCardCorrelExotic[indexExo][matched]->Fill(energy, m02, GetEventWeight());
//      fhNCellsTCardCorrelExotic [indexExo][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    }
  }
  
 
  if( ebin > -1 )
  {
    fhLambda0Lambda1        [ebin][matched]->Fill(m20, m02, GetEventWeight());
    fhNCellsTCardSameAndDiff[ebin][matched]->Fill(nCorrNo, nCorr, GetEventWeight());
    //  if(nCorrNo == 0)
    //      fhLambda0Lambda1AllSameTCard  [ebin][matched]->Fill(m20, m02, GetEventWeight());
    
    if(fStudyExotic)
    {
      fhLambda0Exoticity[ebin][matched]->Fill(exoticity, m02, GetEventWeight());
      fhLambda1Exoticity[ebin][matched]->Fill(exoticity, m20, GetEventWeight());
    //fhLambdaRExoticity[ebin][matched]->Fill(exoticity,lamR, GetEventWeight());
      fhNCellsExoticity [ebin][matched]->Fill(exoticity, nCellWithWeight, GetEventWeight());
      fhTimeExoticity   [ebin][matched]->Fill(exoticity, tCellMax, GetEventWeight());
      
      if(energy > 8)
      {
        fhLambda0ExoticityPerNCell[nCorrInd][nCorrNoInd][matched]->Fill(exoticity, m02, GetEventWeight());
        fhLambda1ExoticityPerNCell[nCorrInd][nCorrNoInd][matched]->Fill(exoticity, m20, GetEventWeight());
      //fhLambdaRExoticityPerNCell[nCorrInd][nCorrNoInd][matched]->Fill(exoticity,lamR, GetEventWeight());
      }
      
      //      if(nCorrNo == 0)
      //      {
      //        fhLambda0ExoticityAllSameTCard[ebin][matched]->Fill(exoticity, m02, GetEventWeight());
      //        fhLambda1ExoticityAllSameTCard[ebin][matched]->Fill(exoticity, m20, GetEventWeight());
      //        fhLambdaRExoticityAllSameTCard[ebin][matched]->Fill(exoticity,lamR, GetEventWeight());
      //        fhNCellsExoticityAllSameTCard [ebin][matched]->Fill(exoticity, nCellWithWeight, GetEventWeight());
      //      }
      
      if ( exoticity > 0.97 )
        fhNCellsTCardSameAndDiffExotic[ebin][matched]->Fill(nCorrNo, nCorr, GetEventWeight());
      
      // Track matching residuals
      fhTMPhiResidualExoticity[ebin]->Fill(exoticity,dphi);
      fhTMEtaResidualExoticity[ebin]->Fill(exoticity,deta);
      
      //    if(nCorrNo==0)
      //    {
      //      fhTMPhiResidualExoticityAllSameTCard[ebin]->Fill(exoticity,dphi);
      //      fhTMEtaResidualExoticityAllSameTCard[ebin]->Fill(exoticity,deta);
      //    }
    }
  }

  Float_t nCellRat = nCorr*1. / ((nCorr+nCorrNo)*1.); 
  fhNCellsTCardSameAndDiffFraction[matched]->Fill(energy, nCellRat, GetEventWeight());
  if ( fStudyExotic && exoticity > 0.97 )
    fhNCellsTCardSameAndDiffFractionExotic[matched]->Fill(energy, nCellRat, GetEventWeight());
  
//  if(nCorr > 0)
//  {
//    Int_t index = -1;
//    if      (!sameRow &&  sameCol && !other ) index = 0;
//    else if (!sameRow && !sameCol &&  other ) index = 1;
//    else if (!sameRow &&  sameCol &&  other ) index = 2;
//    else if ( sameRow &&  sameCol && !other ) index = 3;
//    else if ( sameRow && !sameCol &&  other ) index = 4;
//    else if ( sameRow &&  sameCol &&  other ) index = 5;
//    else if ( sameRow && !sameCol && !other ) index = 6;
//    else printf("case not considered: sameRow %d, sameCol %d, other %d, nearRow %d\n",sameRow,sameCol,other,nearRow); 
//    
//    if(index >= 0)
//    {
//      fhLambda0TCardCorrel[index][matched]->Fill(energy, m02, GetEventWeight());
//      fhNCellsTCardCorrel [index][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//      fhExoticTCardCorrel [index][matched]->Fill(energy, exoticity, GetEventWeight());
//    }
//
//    // Comment out, no special effect observed    
//    Int_t indexNR = -1; 
//    if ( nearRow )
//    {
//      if      (!sameRow &&  sameCol && !other ) indexNR = 0;
//      else if (!sameRow && !sameCol &&  other ) indexNR = 1;
//      else if (!sameRow &&  sameCol &&  other ) indexNR = 2;
//      else if ( sameRow &&  sameCol && !other ) indexNR = 3;
//      else if ( sameRow && !sameCol &&  other ) indexNR = 4;
//      else if ( sameRow &&  sameCol &&  other ) indexNR = 5;
//      else printf("\t near row case not considered!: sameRow %d, sameCol %d, other %d\n",sameRow,sameCol,other); 
//    }
//    
//    if ( indexNR >= 0 )
//    {
//      fhLambda0TCardCorrelNearRow[indexNR][matched]->Fill(energy, m02, GetEventWeight());
//      fhNCellsTCardCorrelNearRow [indexNR][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    }
//    
//    if ( sameTCard2ndMax )
//    {
//      Int_t index2nd = -1;
//      if     ( TMath::Abs(rowDiff2) == 1 && TMath::Abs(colDiff2) != 1 ) index2nd = 0;
//      else if( TMath::Abs(rowDiff2) != 1 && TMath::Abs(colDiff2) == 1 ) index2nd = 1;
//      else if( TMath::Abs(rowDiff2) == 1 && TMath::Abs(colDiff2) == 1 ) index2nd = 2;
//      else                                                              index2nd = 3;
//      
//      fhLambda0TCardCorrel2ndMax[index2nd][matched]->Fill(energy, m02, GetEventWeight());
//      fhNCellsTCardCorrel2ndMax [index2nd][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    }
//  }
//  
//  Int_t indexOtherTCard = -1;
//  if      ( nCorr == 0 && nCorr2 == 0 ) indexOtherTCard = 6;
//  else if ( nCorr == 0 && nCorr2  > 0 )
//  {
//    if      (  sameRow2 && !sameCol2 && !other2) indexOtherTCard = 0; 
//    else if ( !sameRow2 &&  sameCol2 && !other2) indexOtherTCard = 1;
//    else                                         indexOtherTCard = 2;
//  }  
//  else if ( nCorr  > 0 && nCorr2  > 0 )
//  {
//    if      (  sameRow2 && !sameCol2 && !other2) indexOtherTCard = 3; 
//    else if ( !sameRow2 &&  sameCol2 && !other2) indexOtherTCard = 4;
//    else                                         indexOtherTCard = 5;
//  }
//  
//  if ( indexOtherTCard >= 0 )
//  {
//    fhLambda0TCardCorrelOtherTCard[indexOtherTCard][matched]->Fill(energy, m02, GetEventWeight());
//    fhNCellsTCardCorrelOtherTCard [indexOtherTCard][matched]->Fill(energy, nCellWithWeight, GetEventWeight());
//    fhExoticTCardCorrelOtherTCard [indexOtherTCard][matched]->Fill(energy, exoticity, GetEventWeight());
//    if ( energy >= 2 && energy < 8 ) 
//      fhColRowTCardCorrelOtherTCardLowE[indexOtherTCard][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//    else if ( energy >= 8 ) 
//      fhColRowTCardCorrelOtherTCardHighE[indexOtherTCard][matched]->Fill(icolAbs, irowAbs, GetEventWeight());
//  } 
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
  fhClusterTimeEnergy   ->Fill(clus->E(), tof            , GetEventWeight());

  fhLambda0             ->Fill(clus->E(), clus->GetM02()       , GetEventWeight());
  fhLambda1             ->Fill(clus->E(), clus->GetM20()       , GetEventWeight());
//fhDispersion          ->Fill(clus->E(), clus->GetDispersion(), GetEventWeight());
  fhNLocMax             ->Fill(clus->E(), GetCaloUtils()->GetNumberOfLocalMaxima(clus,cells), GetEventWeight());

  if(fFillClusterMaxCellHisto)
  {
    fhClusterMaxCellDiff  ->Fill(clus->E(), maxCellFraction, GetEventWeight());
    fhClusterMaxCellECross->Fill(clus->E(), eCrossFrac     , GetEventWeight());
  }
  
  if(fStudyM02Dependence)
  {
    fhClusterMaxCellDiffM02  ->Fill(clus->E(), maxCellFraction, clus->GetM02(), GetEventWeight());
    fhClusterMaxCellECrossM02->Fill(clus->E(), eCrossFrac     , clus->GetM02(), GetEventWeight());
    fhClusterTimeEnergyM02   ->Fill(clus->E(), tof            , clus->GetM02(), GetEventWeight());
  }
  
  if(fStudyClustersAsymmetry) ClusterAsymmetryHistograms(clus,absIdMax,kTRUE);
  
  Int_t    nModule = GetModuleNumber(clus);
  Int_t    nCaloCellsPerCluster = clus->GetNCells();
  
  // Clusters in event time difference
  if(fFillPi0PairDiffTime)
  {
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
    } // loop 
  } // fill cluster pair time diff        
  
  if(nCaloCellsPerCluster > 1 && 
     (fFillClusterMaxCellHisto || fStudyM02Dependence || fFillAllCellTimeHisto))
  {
    // Check time of cells respect to max energy cell
    
//    if(fFillAllCellTimeHisto && fFillClusterMaxCellHisto) 
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
      
      if(fFillClusterMaxCellHisto)
      {
        fhClusterMaxCellCloseCellRatio->Fill(clus->E(), frac, GetEventWeight());
        fhClusterMaxCellCloseCellDiff ->Fill(clus->E(), cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId), GetEventWeight());
      }
      
      if(fStudyM02Dependence)
      {
        fhClusterMaxCellCloseCellRatioM02->Fill(clus->E(), frac, clus->GetM02(), GetEventWeight());
        fhClusterMaxCellCloseCellDiffM02 ->Fill(clus->E(), cells->GetCellAmplitude(absIdMax)-cells->GetCellAmplitude(absId), clus->GetM02(),GetEventWeight());
      }
      
      if(fFillAllCellTimeHisto) 
      {
        Double_t time  = cells->GetCellTime(absId);
        GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
        
        Float_t diff = (tmax-(time*1.0e9-fConstantTimeShift));
        
        if(fFillClusterMaxCellHisto)
          fhCellTimeSpreadRespectToCellMax->Fill(clus->E(), diff, GetEventWeight());
        
        if(fStudyM02Dependence) 
          fhCellTimeSpreadRespectToCellMaxM02->Fill(clus->E(), diff, clus->GetM02(), GetEventWeight());
        
        if( fFillAllCellAbsIdHistograms && fFillAllCellHistograms && 
           TMath::Abs(TMath::Abs(diff) > 100) && clus->E() > 1 ) 
          fhCellIdCellLargeTimeSpread->Fill(absId, GetEventWeight());
      }
      
    } // Fill cell-cluster histogram loop
  } // Check time and energy of cells respect to max energy cell if cluster of more than 1 cell
  
  Float_t e   = fClusterMomentum.E();
  Float_t pt  = fClusterMomentum.Pt();
  Float_t eta = fClusterMomentum.Eta();
  Float_t phi = GetPhi(fClusterMomentum.Phi());
  
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
  
  if(fFillEBinAcceptanceHisto)
  {
    Int_t icol = -1, irow = -1, iRCU = -1, icolAbs = -1, irowAbs = -1;
    GetModuleNumberCellIndexesAbsCaloMap(absIdMax,GetCalorimeter(), icol, irow, iRCU, icolAbs, irowAbs);
    
    for(Int_t ie = 0; ie < fNEBinCuts; ie++)
    {
      if( e >= fEBinCuts[ie] && e < fEBinCuts[ie+1] )
      {
        fhEBinClusterEtaPhi[ie]->Fill(eta,phi,GetEventWeight()) ;
        
        fhEBinClusterColRow[ie]->Fill(icolAbs,irowAbs,GetEventWeight()) ;
      }
    }
  }
  else if ( fFillAllTH3 ) fhEtaPhiE->Fill(eta, phi, e, GetEventWeight());
  else if ( e > 0.5     ) fhEtaPhi ->Fill(eta, phi,    GetEventWeight());
  
  // Cells per cluster
  fhNCellsPerCluster->Fill(e, nCaloCellsPerCluster, GetEventWeight());
  
  if(fStudyM02Dependence)
    fhNCellsPerClusterM02->Fill(e, nCaloCellsPerCluster, clus->GetM02(),GetEventWeight());

  if(e > 100) 
    fhNCellsPerClusterWeirdMod->Fill(nCaloCellsPerCluster, nModule, GetEventWeight());
    
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
/// * Exotic cluster histograms
/// * Cells in cluster
/// * Invariant mass
/// * Matched clusters histograms
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCalorimeterQA::ClusterLoopHistograms(const TObjArray *caloClusters,
                                                AliVCaloCells* cells)
{
  if(!fFillAllClusterHistograms) return;

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
        
    // SuperModule number of cluster
    nModule = GetModuleNumber(clus);
    if ( nModule < fFirstModule || nModule > fLastModule ) 
    {
      AliDebug(1,Form("Cluster module out of range %d",nModule));
      continue ;
    }
    
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    // Cut on time of clusters
    Double_t tof = clus->GetTOF()*1.e9;
    if(tof>400) tof-=fConstantTimeShift;

    if( tof < fTimeCutMin || tof > fTimeCutMax )
    { 
      AliDebug(1,Form("Remove cluster with TOF %2.2f",tof));
      continue;
    }    
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if(IsFiducialCutOn()) in =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum.Eta(),fClusterMomentum.Phi(),GetCalorimeter()) ;
    if(!in)
    {
      AliDebug(1,Form("Remove cluster with phi %2.2f and eta %2.2f",
                      GetPhi(fClusterMomentum.Phi())*TMath::RadToDeg(),fClusterMomentum.Eta()));
      continue;
    }
    
    // MC labels
    nLabel = clus->GetNLabels();
    labels = clus->GetLabels();
    
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
      
      if(fFillClusterMaxCellHisto)
        fhClusterMaxCellDiffNoCut->Fill(clus->E(), maxCellFraction, GetEventWeight());
    }
    
    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    GetCaloUtils()->RecalibrateCellAmplitude(ampMax,GetCalorimeter(), absIdMax);
    
    //if(fStudyExotic) ExoticHistograms(absIdMax, ampMax, clus, cells);

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

    //
    ClusterHistograms(clus, caloClusters, cells, absIdMax, 
                      maxCellFraction, eCrossFrac, tmax);
    
    nCaloClustersAccepted++;
    
    //
    if(fStudyTCardCorrelation) ChannelCorrelationInTCard(clus, cells, matched, absIdMax, eCrossFrac);
    
    //
    if(nModule >=0 && nModule < fNModules && fClusterMomentum.E() > 2*fCellAmpMin)
    {
     nClustersInModule[nModule]++;
     if(clus->E() > 0.5)
      energyInModule  [nModule] += clus->E();
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
    if ( imod < fFirstModule || imod > fLastModule ) continue ;

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
    
  // Play with the MC stack if available
  Int_t label = labels[0];
  
  if(label < 0) 
  {
    AliDebug(1,Form(" *** bad label ***:  label %d", label));
    return kFALSE;
  }
  
  if( label >= GetMC()->GetNumberOfTracks()) 
  {
    AliDebug(1,Form("*** large label ***:  label %d, n tracks %d", label, GetMC()->GetNumberOfTracks()));
    return kFALSE;
  }

  Float_t e   = fClusterMomentum.E();
  Float_t eta = fClusterMomentum.Eta();
  Float_t phi = fClusterMomentum.Phi();
  if(phi < 0) phi +=TMath::TwoPi();
  
  AliAODMCParticle * aodprimary  = 0x0;
  TParticle * primary = 0x0;

  Int_t pdg0    =-1; Int_t status = -1; Int_t iMother = -1; Int_t iParent = -1;
  Float_t vxMC  = 0; Float_t vyMC  = 0;
  Float_t eMC   = 0; //Float_t ptMC= 0;
  Float_t phiMC = 0; Float_t etaMC = 0;
  Int_t charge  = 0;
  
  // Check the origin.
  Int_t tag = GetMCAnalysisUtils()->CheckOrigin(labels,nLabels, GetReader(),GetCalorimeter());
  
  if     ( GetReader()->ReadStack() && 
          !GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCUnknown))
  { 
    primary = GetMC()->Particle(label);
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
      primary = GetMC()->Particle(iParent);
      pdg = TMath::Abs(primary->GetPdgCode());
      
      AliDebug(2,"Converted cluster!. Find before conversion:");

      while((pdg == 22 || pdg == 11) && status != 1)
      {
        Int_t iMotherOrg = iMother;
        iMother = iParent;
        primary = GetMC()->Particle(iMother);
        status  = primary->GetStatusCode();
        pdg     = TMath::Abs(primary->GetPdgCode());
        iParent = primary->GetFirstMother();
        
        // If gone too back and non stable, assign the decay photon/electron
        // there are other possible decays, ignore them for the moment
        if(pdg==111 || pdg==221)
        {
          primary = GetMC()->Particle(iMotherOrg);
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
        //printf("iMother %d, pdg %d, iParent %d, pdg %d\n",iMother,pdg,iParent,GetMC()->Particle(iParent)->GetPdgCode());
        iMother = iParent;
        primary = GetMC()->Particle(iMother);
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
    aodprimary = (AliAODMCParticle*) GetMC()->GetTrack(label);
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
      aodprimary = (AliAODMCParticle*) GetMC()->GetTrack(iParent);
      pdg = TMath::Abs(aodprimary->GetPdgCode());
        
      while ((pdg == 22 || pdg == 11) && !aodprimary->IsPhysicalPrimary()) 
      {
        Int_t iMotherOrg = iMother;
        iMother    = iParent;
        aodprimary = (AliAODMCParticle*) GetMC()->GetTrack(iMother);
        status     = aodprimary->IsPrimary();
        iParent    = aodprimary->GetMother();
        pdg        = TMath::Abs(aodprimary->GetPdgCode());

        // If gone too back and non stable, assign the decay photon/electron
        // there are other possible decays, ignore them for the moment
        if( pdg == 111 || pdg == 221 )
        {
          aodprimary = (AliAODMCParticle*) GetMC()->GetTrack(iMotherOrg);
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
        aodprimary = (AliAODMCParticle*) GetMC()->GetTrack(iMother);
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
  snprintf(onePar,buffersize,"Cluster M02: EMCAL > %2.2f ;",fEMCALClusterM02Min) ;
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

//___________________________________________________
/// Create histograms to be saved in output file and
/// store them in the output container.
//___________________________________________________
TList * AliAnaCalorimeterQA::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("QAHistos") ; 
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  Int_t totalSM = fLastModule-fFirstModule+1;

  //printf("N SM %d, first SM %d, last SM %d, total %d\n",fNModules,fFirstModule,fLastModule, totalSM);
  
  // Histogram binning and ranges
  // 
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
  
  // TM residuals
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  // Cell column-row histograms, see base class for data members setting
  //fNMaxColsFull+2,-1.5,fNMaxColsFull+0.5, fNMaxRowsFull+2,-1.5,fNMaxRowsFull+0.5
  Int_t   ncolcell   = fNMaxColsFull+2;
  Float_t colcellmin = -1.5;
  Float_t colcellmax = fNMaxColsFull+0.5;
  
  Int_t   nrowcell   = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  Float_t rowcellmin = fNMaxRowsFullMin-1.5;
  Float_t rowcellmax = fNMaxRowsFullMax+0.5;
  
  //
  // Init histograms
  //
  
  // Calorimeter cluster histograms
  //
  if(fFillAllClusterHistograms)
  {
    fhE  = new TH1F ("hE","#it{E} reconstructed clusters ", nptbins*5,ptmin,ptmax*5);  
    fhE->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhE);
    
    fhPt  = new TH1F ("hPt","#it{p}_{T} reconstructed clusters", nptbins,ptmin,ptmax);
    fhPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPt);
    
    fhNClusters  = new TH1F ("hNClusters","# clusters", nclbins,nclmin,nclmax); 
    fhNClusters->SetXTitle("#it{n}_{clusters}");
    outputContainer->Add(fhNClusters);
    
    fhNCellsPerCluster  = new TH2F ("hNCellsPerCluster","# cells per cluster vs energy",nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
    fhNCellsPerCluster->SetXTitle("#it{E} (GeV)");
    fhNCellsPerCluster->SetYTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerCluster);
    
    fhNCellsPerClusterWeirdMod  = new TH2F ("hNCellsPerClusterWeirdMod","# cells per cluster, E > 100 GeV, per SM", 
                                            nceclbins*2,nceclmin,nceclmax*2,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhNCellsPerClusterWeirdMod->SetYTitle("SM number");
    fhNCellsPerClusterWeirdMod->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsPerClusterWeirdMod);
    
    // Acceptance plots
    //
    fhPhi  = new TH1F ("hPhi","#varphi reconstructed clusters ",nphibins,phimin,phimax);
    fhPhi->SetXTitle("#varphi (rad)");
    outputContainer->Add(fhPhi);
    
    fhEta  = new TH1F ("hEta","#eta reconstructed clusters ",netabins,etamin,etamax);
    fhEta->SetXTitle("#eta ");
    outputContainer->Add(fhEta);
    
    if(fFillEBinAcceptanceHisto)
    {
      for(Int_t ie=0; ie<fNEBinCuts; ie++)
      {
        fhEBinClusterEtaPhi[ie] = new TH2F
        (Form("hEBin%d_Cluster_EtaPhi",ie),
         Form("#eta vs #varphi, cluster, %2.2f<#it{p}_{T}<%2.2f GeV/#it{c}",fEBinCuts[ie],fEBinCuts[ie+1]),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEBinClusterEtaPhi[ie]->SetYTitle("#varphi (rad)");
        fhEBinClusterEtaPhi[ie]->SetXTitle("#eta");
        outputContainer->Add(fhEBinClusterEtaPhi[ie]) ;
        
        fhEBinClusterColRow[ie] = new TH2F
        (Form("hEBin%d_Cluster_ColRow",ie),
         Form("column vs row, cluster max E cell, %2.2f<#it{p}_{T}<%2.2f GeV/#it{c}",fEBinCuts[ie],fEBinCuts[ie+1]),
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhEBinClusterColRow[ie]->SetYTitle("row");
        fhEBinClusterColRow[ie]->SetXTitle("column");
        outputContainer->Add(fhEBinClusterColRow[ie]) ;
        
        if(fFillAllCellHistograms)
        {
          fhEBinCellColRow[ie] = new TH2F
          (Form("hEBin%d_Cell_ColRow",ie),
           Form("column vs row, cell, %2.2f<#it{p}_{T}<%2.2f GeV/#it{c}",fEBinCuts[ie],fEBinCuts[ie+1]),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhEBinCellColRow[ie]->SetYTitle("row");
          fhEBinCellColRow[ie]->SetXTitle("column");
          outputContainer->Add(fhEBinCellColRow[ie]) ;
        }
      }
    }
    else
    {
      if(fFillAllTH3)
      {
        fhEtaPhiE  = new TH3F ("hEtaPhiE","#eta vs #varphi vs energy, reconstructed clusters",
                               netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
        fhEtaPhiE->SetXTitle("#eta ");
        fhEtaPhiE->SetYTitle("#varphi (rad)");
        fhEtaPhiE->SetZTitle("#it{E} (GeV) ");
        outputContainer->Add(fhEtaPhiE);
      }
      else 
      {
        fhEtaPhi  = new TH2F ("hEtaPhi","#eta vs #varphi for #it{E} > 0.5 GeV, reconstructed clusters",
                              netabins,etamin,etamax,nphibins,phimin,phimax); 
        fhEtaPhi->SetXTitle("#eta ");
        fhEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhi);
      }
    }
    
    fhClusterTimeEnergy  = new TH2F ("hClusterTimeEnergy","energy vs TOF, reconstructed clusters",
                                     nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhClusterTimeEnergy->SetXTitle("#it{E} (GeV) ");
    fhClusterTimeEnergy->SetYTitle("TOF (ns)");
    outputContainer->Add(fhClusterTimeEnergy);
    
    if(fFillPi0PairDiffTime)
    {
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
    }
    
    // Shower shape
    //
    //fhDispersion  = new TH2F ("hDispersion","shower shape, Dispersion^{2} vs E for bad cluster ",
    //                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    //fhDispersion->SetXTitle("#it{E}_{cluster} (GeV)");
    //fhDispersion->SetYTitle("Dispersion");
    //outputContainer->Add(fhDispersion);      
    
    fhLambda0  = new TH2F ("hLambda0","shower shape, #lambda^{2}_{0} vs E",
                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLambda0->SetXTitle("#it{E}_{cluster} (GeV)");
    fhLambda0->SetYTitle("#lambda^{2}_{0}");
    outputContainer->Add(fhLambda0); 
    
    fhLambda1  = new TH2F ("hLambda1","shower shape, #lambda^{2}_{1} vs E",
                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLambda1->SetXTitle("#it{E}_{cluster} (GeV)");
    fhLambda1->SetYTitle("#lambda^{2}_{1}");
    outputContainer->Add(fhLambda1); 
    
    fhNLocMax  = new TH2F ("hNLocMax","#it{n}_{LM}vs E",
                           nptbins,ptmin,ptmax,10,0,10); 
    fhNLocMax->SetXTitle("#it{E}_{cluster} (GeV)");
    fhNLocMax->SetYTitle("#it{n}_{LM}");
    outputContainer->Add(fhNLocMax); 
    
    if(fStudyTCardCorrelation)
    {
      fhColRowHighEPosTime = new TH2F
      ("hColRowHighEPosTime",
       "column vs row, exo < 0.97, E > 8 GeV, t > 5 ns",
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhColRowHighEPosTime->SetYTitle("row");
      fhColRowHighEPosTime->SetXTitle("column");
      outputContainer->Add(fhColRowHighEPosTime) ;
      
      fhColRowHighENegTime = new TH2F
      ("hColRowHighENegTime",
       "column vs row, exo < 0.97, E > 8 GeV, t < -5 ns",
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhColRowHighENegTime->SetYTitle("row");
      fhColRowHighENegTime->SetXTitle("column");
      outputContainer->Add(fhColRowHighENegTime) ;
      
      fhColRowHighENulTime = new TH2F
      ("hColRowHighENulTime",
       "column vs row, exo < 0.97, E > 8 GeV, -5 < t < 5 ns",
       ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhColRowHighENulTime->SetYTitle("row");
      fhColRowHighENulTime->SetXTitle("column");
      outputContainer->Add(fhColRowHighENulTime) ;
      
      if(fStudyExotic)
      {
        fhColRowExoticHighE1CellPosTime = new TH2F
        ("hColRowExoticHighE1CellPosTime",
         "column vs row, 1 cell, E > 8 GeV, t > 5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighE1CellPosTime->SetYTitle("row");
        fhColRowExoticHighE1CellPosTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighE1CellPosTime) ;
        
        fhColRowExoticHighEPosTime = new TH2F
        ("hColRowExoticHighEPosTime",
         "column vs row, exo > 0.97, E > 8 GeV, t > 5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighEPosTime->SetYTitle("row");
        fhColRowExoticHighEPosTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighEPosTime) ;
        
        fhColRowExoticHighE1CellNegTime = new TH2F
        ("hColRowExoticHighE1CellNegTime",
         "column vs row, 1 cell, E > 8 GeV, t < -5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighE1CellNegTime->SetYTitle("row");
        fhColRowExoticHighE1CellNegTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighE1CellNegTime) ;
        
        fhColRowExoticHighENegTime = new TH2F
        ("hColRowExoticHighENegTime",
         "column vs row, exo > 0.97, E > 8 GeV, t < -5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighENegTime->SetYTitle("row");
        fhColRowExoticHighENegTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighENegTime) ;
        
        fhColRowExoticHighE1CellNulTime = new TH2F
        ("hColRowExoticHighE1CellNulTime",
         "column vs row, 1 cell, E > 8 GeV, -5 < t < 5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighE1CellNulTime->SetYTitle("row");
        fhColRowExoticHighE1CellNulTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighE1CellNulTime) ;
        
        fhColRowExoticHighENulTime = new TH2F
        ("hColRowExoticHighENulTime",
         "column vs row, exo > 0.97, E > 8 GeV, -5 < t < 5 ns",
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowExoticHighENulTime->SetYTitle("row");
        fhColRowExoticHighENulTime->SetXTitle("column");
        outputContainer->Add(fhColRowExoticHighENulTime) ;
      }
      
      TString add[] = {"","TrackMatched"};
      for(Int_t tm = 0; tm < 2; tm++)
      {      
        fhColRowTCardCorrNoSelectionLowE[tm] = new TH2F
        (Form("hColRowTCardCorrNoSelectionLowE%s",add[tm].Data()),
         Form("column vs row, max E cell for TCard correlation selected clusters, 5 < E < 8  GeV %s",add[tm].Data()),
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowTCardCorrNoSelectionLowE[tm]->SetYTitle("row");
        fhColRowTCardCorrNoSelectionLowE[tm]->SetXTitle("column");
        outputContainer->Add(fhColRowTCardCorrNoSelectionLowE[tm]) ;
        
        fhColRowTCardCorrNoSelectionHighE[tm] = new TH2F
        (Form("hColRowTCardCorrNoSelectionHighE%s",add[tm].Data()),
         Form("column vs row, max E cell for TCard correlation selected clusters, E > 8 GeV %s",add[tm].Data()),
         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        fhColRowTCardCorrNoSelectionHighE[tm]->SetYTitle("row");
        fhColRowTCardCorrNoSelectionHighE[tm]->SetXTitle("column");
        outputContainer->Add(fhColRowTCardCorrNoSelectionHighE[tm]) ;
        
        //
        
        fhNCellsTCardCorrNoSelection[tm]  = new TH2F 
        (Form("hNCellsTCardCorrNoSelection%s",add[tm].Data()),
         Form("# custer # cells vs #it{E} %s",add[tm].Data()),
         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        fhNCellsTCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhNCellsTCardCorrNoSelection[tm]->SetYTitle("#it{n}_{cells}");
        outputContainer->Add(fhNCellsTCardCorrNoSelection[tm]);
        
        fhNCellsTCardCorrWithWeightNoSelection[tm]  = new TH2F 
        (Form("hNCellsTCardCorrWithWeightNoSelection%s",add[tm].Data()),
         Form("custer # cells vs #it{E}, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        fhNCellsTCardCorrWithWeightNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhNCellsTCardCorrWithWeightNoSelection[tm]->SetYTitle("#it{n}_{cells}");
        outputContainer->Add(fhNCellsTCardCorrWithWeightNoSelection[tm]);
        
        fhNCellsTCardCorrRatioWithWeightNoSelection[tm]  = new TH2F 
        (Form("hNCellsTCardCorrRatioWithWeightNoSelection%s",add[tm].Data()),
         Form("custer # cells vs #it{E}, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax, 100,0,1); 
        fhNCellsTCardCorrRatioWithWeightNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhNCellsTCardCorrRatioWithWeightNoSelection[tm]->SetYTitle("#it{n}^{w>0.01}_{cells} / #it{n}_{cells}");
        outputContainer->Add(fhNCellsTCardCorrRatioWithWeightNoSelection[tm]);
        
        fhLambda0TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda0TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{0} vs #it{E} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda0TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda0TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{0}");
        outputContainer->Add(fhLambda0TCardCorrNoSelection[tm]); 
        
        fhLambda1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{1} vs #it{E} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda1TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda1TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{1}");
        outputContainer->Add(fhLambda1TCardCorrNoSelection[tm]); 
        
        fhLambda0NLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda0NLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{0} vs #it{E}, nlm=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda0NLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda0NLM1TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{0}");
        outputContainer->Add(fhLambda0NLM1TCardCorrNoSelection[tm]); 
        
        fhLambda1NLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda1NLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{1} vs #it{E}, nlm=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda1NLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda1NLM1TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{1}");
        outputContainer->Add(fhLambda1NLM1TCardCorrNoSelection[tm]); 
        
        fhLambda0NLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda0NLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{0} vs #it{E}, nlm=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda0NLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda0NLM2TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{0}");
        outputContainer->Add(fhLambda0NLM2TCardCorrNoSelection[tm]); 
        
        fhLambda1NLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambda1NLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{2}_{1} vs #it{E}, nlm=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLambda1NLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambda1NLM2TCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{1}");
        outputContainer->Add(fhLambda1NLM2TCardCorrNoSelection[tm]); 
        
        
        fhLambdaRTCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLambdaRTCardCorrNoSelection%s",add[tm].Data()),
         Form("#lambda^{1}_{0}/#lambda^{2}_{0} vs #it{E} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,110,0,1.1); 
        fhLambdaRTCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhLambdaRTCardCorrNoSelection[tm]->SetYTitle("#lambda^{2}_{1}/#lambda^{2}_{0}");
        outputContainer->Add(fhLambdaRTCardCorrNoSelection[tm]); 
        
        fhNLocMaxTCardCorrNoSelection[tm]  = new TH2F 
        (Form("hNLocMaxTCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{n}_{LM} vs E %s",add[tm].Data()),
         nptbins,ptmin,ptmax,10,0,10); 
        fhNLocMaxTCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhNLocMaxTCardCorrNoSelection[tm]->SetYTitle("#it{n}_{LM}");
        outputContainer->Add(fhNLocMaxTCardCorrNoSelection[tm]); 
        
        fhEMaxRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hEMaxRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs E, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhEMaxRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhEMaxRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
        outputContainer->Add(fhEMaxRatNLM1TCardCorrNoSelection[tm]); 
        
        fhEMaxRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hEMaxRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs E, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhEMaxRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhEMaxRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
        outputContainer->Add(fhEMaxRatNLM2TCardCorrNoSelection[tm]); 
        
        fhEMaxRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hEMaxRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs E, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhEMaxRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhEMaxRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
        outputContainer->Add(fhEMaxRatNLM3TCardCorrNoSelection[tm]); 
        
        fhE2ndRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd max}/#it{E}_{cluster} vs E, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndRatNLM1TCardCorrNoSelection[tm]); 
        
        fhE2ndRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs E, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndRatNLM2TCardCorrNoSelection[tm]); 
        
        fhE2ndRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs E, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndRatNLM3TCardCorrNoSelection[tm]);
        
        
        fhE2ndEMaxRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndEMaxRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max} vs E, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndEMaxRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndEMaxRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndEMaxRatNLM1TCardCorrNoSelection[tm]); 
        
        fhE2ndEMaxRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndEMaxRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max} vs E, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndEMaxRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndEMaxRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndEMaxRatNLM2TCardCorrNoSelection[tm]); 
        
        fhE2ndEMaxRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndEMaxRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max} vs E, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndEMaxRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndEMaxRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndEMaxRatNLM3TCardCorrNoSelection[tm]);
        
        /////
        
        fhE2ndSameRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd max}/#it{E}_{cluster} vs E, 2nd in same TCard as leading, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndSameRatNLM1TCardCorrNoSelection[tm]); 
        
        fhE2ndSameRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs E, 2nd in same TCard as leading, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndSameRatNLM2TCardCorrNoSelection[tm]); 
        
        fhE2ndSameRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs E, 2nd in same TCard as leading, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
        outputContainer->Add(fhE2ndSameRatNLM3TCardCorrNoSelection[tm]);
        
        
        fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameEMaxRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max} vs E, 2nd in same TCard as leading, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndSameEMaxRatNLM1TCardCorrNoSelection[tm]); 
        
        fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameEMaxRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max} vs E, 2nd in same TCard as leading, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndSameEMaxRatNLM2TCardCorrNoSelection[tm]); 
        
        fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hE2ndSameEMaxRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max} vs E, 2nd in same TCard as leading, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
        outputContainer->Add(fhE2ndSameEMaxRatNLM3TCardCorrNoSelection[tm]);
        
        /////
        
        fhECellClusRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellClusRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E_{cluster}, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellClusRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellClusRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellClusRatNLM1TCardCorrNoSelection[tm]); 
        
        fhECellClusRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellClusRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellClusRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellClusRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellClusRatNLM2TCardCorrNoSelection[tm]); 
        
        fhECellClusRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellClusRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellClusRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellClusRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellClusRatNLM3TCardCorrNoSelection[tm]);
        
        fhLogECellNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E_{cluster}, #it{n}_{LM}=1, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellNLM1TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellNLM1TCardCorrNoSelection[tm]); 
        
        fhLogECellNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E, #it{n}_{LM}=2, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellNLM2TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellNLM2TCardCorrNoSelection[tm]); 
        
        fhLogECellNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E, #it{n}_{LM}>2, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellNLM3TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellNLM3TCardCorrNoSelection[tm]);
        
        
        fhECellWeightNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellWeightNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs E_{cluster}, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellWeightNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellWeightNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellWeightNLM1TCardCorrNoSelection[tm]); 
        
        fhECellWeightNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellWeightNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs E, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellWeightNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellWeightNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellWeightNLM2TCardCorrNoSelection[tm]); 
        
        fhECellWeightNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellWeightNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))vs E, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellWeightNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellWeightNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellWeightNLM3TCardCorrNoSelection[tm]);
        
        
        fhECellSameClusRatNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameClusRatNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E_{cluster}, cell from same T-Card as leading, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellSameClusRatNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameClusRatNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellSameClusRatNLM1TCardCorrNoSelection[tm]); 
        
        fhECellSameClusRatNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameClusRatNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E, cell from same T-Card as leading, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellSameClusRatNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameClusRatNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellSameClusRatNLM2TCardCorrNoSelection[tm]); 
        
        fhECellSameClusRatNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameClusRatNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{E}_{cell}/#it{E}_{cluster} vs E, cell from same T-Card as leading, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhECellSameClusRatNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameClusRatNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
        outputContainer->Add(fhECellSameClusRatNLM3TCardCorrNoSelection[tm]);
        
        fhLogECellSameNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellSameNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E_{cluster}, cell from same T-Card as leading, #it{n}_{LM}=1, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellSameNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellSameNLM1TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellSameNLM1TCardCorrNoSelection[tm]); 
        
        fhLogECellSameNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellSameNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E, #it{n}_{LM}=2, cell from same T-Card as leading, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellSameNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellSameNLM2TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellSameNLM2TCardCorrNoSelection[tm]); 
        
        fhLogECellSameNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hLogECellSameNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("log(#it{E}_{cell}) vs E, #it{n}_{LM}>2, cell from same T-Card as leading, w > 0.01 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,150,-3,3); 
        fhLogECellSameNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhLogECellSameNLM3TCardCorrNoSelection[tm]->SetYTitle("log(#it{E}_{cell})");
        outputContainer->Add(fhLogECellSameNLM3TCardCorrNoSelection[tm]);
        
        
        fhECellSameWeightNLM1TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameWeightNLM1TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs E_{cluster}, cell from same T-Card as leading, #it{n}_{LM}=1 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellSameWeightNLM1TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameWeightNLM1TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellSameWeightNLM1TCardCorrNoSelection[tm]); 
        
        fhECellSameWeightNLM2TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameWeightNLM2TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs E,  cell from same T-Card as leading, #it{n}_{LM}=2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellSameWeightNLM2TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameWeightNLM2TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellSameWeightNLM2TCardCorrNoSelection[tm]); 
        
        fhECellSameWeightNLM3TCardCorrNoSelection[tm]  = new TH2F 
        (Form("hECellSameWeightNLM3TCardCorrNoSelection%s",add[tm].Data()),
         Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))vs E, cell from same T-Card as leading, #it{n}_{LM}>2 %s",add[tm].Data()),
         nptbins,ptmin,ptmax,90,0,4.5); 
        fhECellSameWeightNLM3TCardCorrNoSelection[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhECellSameWeightNLM3TCardCorrNoSelection[tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
        outputContainer->Add(fhECellSameWeightNLM3TCardCorrNoSelection[tm]);
        
        fhExoticTCardCorrNoSelection[tm]  = new TH2F 
        (Form("hExoticTCardCorrNoSelection%s",add[tm].Data()),
         Form("exoticity vs #it{E} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,200,-1,1); 
        fhExoticTCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
        fhExoticTCardCorrNoSelection[tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
        outputContainer->Add(fhExoticTCardCorrNoSelection[tm]); 
        
        fhNCellsTCardSameAndDiffFraction[tm]  = new TH2F 
        (Form("hNCellsTCardSameAndDiffFraction%s",add[tm].Data()),
         Form("#it{n}_{cells} same TCard vs diff TCard fraction, w > 0.01, %s",add[tm].Data()),
         nptbins,ptmin,ptmax,100,0,1); 
        fhNCellsTCardSameAndDiffFraction[tm]->SetXTitle("#it{E} (GeV)");
        fhNCellsTCardSameAndDiffFraction[tm]->SetYTitle("#it{n}_{cells} - same TCard / #it{n}_{cells} - total");
        outputContainer->Add(fhNCellsTCardSameAndDiffFraction[tm]);             
        
        fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[tm] = new TH2F 
        (Form("hSameRowDiffColAndTCardCellsEnergyDiffClusterE%s",add[tm].Data()),
         Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} vs #it{E}_{cluster} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,200,-10,10); 
        fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[tm]->SetYTitle("#it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} (GeV)");
        outputContainer->Add(fhSameRowDiffColAndTCardCellsEnergyDiffClusterE[tm]); 
        
        fhSameRowDiffColAndTCardCellsTimeDiffClusterE[tm] = new TH2F 
        (Form("hSameRowDiffColAndTCardCellsTimeDiffClusterE%s",add[tm].Data()),
         Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} vs #it{E}_{cluster} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,200,-100,100); 
        fhSameRowDiffColAndTCardCellsTimeDiffClusterE[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
        fhSameRowDiffColAndTCardCellsTimeDiffClusterE[tm]->SetYTitle("#it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} (ns)");
        outputContainer->Add(fhSameRowDiffColAndTCardCellsTimeDiffClusterE[tm]); 
        
        fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[tm] = new TH2F 
        (Form("hSameRowDiffColAndTCardCellsEnergyDiffCellMaxE%s",add[tm].Data()),
         Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} vs #it{E}_{cell max} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,200,-10,10); 
        fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[tm]->SetXTitle("#it{E}_{cell max} (GeV)");
        fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[tm]->SetYTitle("#it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} (GeV)");
        outputContainer->Add(fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxE[tm]); 
        
        fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm] = new TH2F 
        (Form("hSameRowDiffColAndTCardCellsTimeDiffCellMaxE%s",add[tm].Data()),
         Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} vs #it{E}_{cell max} %s",add[tm].Data()),
         nptbins,ptmin,ptmax,200,-100,100); 
        fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]->SetXTitle("#it{E}_{cell max} (GeV)");
        fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]->SetYTitle("#it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} (ns)");
        outputContainer->Add(fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]); 
        
        if(fStudyExotic)
        {
          fhEnergyTime1Cell[tm]  = new TH2F 
          (Form("hEnergyTime1Cell%s",add[tm].Data()),
           Form("#it{t} vs #it{E}, 1 cells cluster %s",add[tm].Data()),
           nptbins,ptmin,ptmax,300,-150,150); 
          fhEnergyTime1Cell[tm]->SetXTitle("#it{E} (GeV)");
          fhEnergyTime1Cell[tm]->SetYTitle("#it{t} (ns)");
          outputContainer->Add(fhEnergyTime1Cell[tm]); 
          
          fhEnergyTimeExotic[tm]  = new TH2F 
          (Form("hEnergyTimeExotic%s",add[tm].Data()),
           Form("#it{t} vs #it{E},  exo > 0.97, %s",add[tm].Data()),
           nptbins,ptmin,ptmax,300,-150,150); 
          fhEnergyTimeExotic[tm]->SetXTitle("#it{E} (GeV)");
          fhEnergyTimeExotic[tm]->SetYTitle("#it{t} (ns)");
          outputContainer->Add(fhEnergyTimeExotic[tm]); 
          
          fhEnergyTimeTCardCorrNoSelection1Cell[tm]  = new TH2F 
          (Form("hEnergyTimeTCardCorrNoSelection1Cell%s",add[tm].Data()),
           Form("#it{t} vs #it{E}, 1 cells cluster %s",add[tm].Data()),
           nptbins,ptmin,ptmax,300,-150,150); 
          fhEnergyTimeTCardCorrNoSelection1Cell[tm]->SetXTitle("#it{E} (GeV)");
          fhEnergyTimeTCardCorrNoSelection1Cell[tm]->SetYTitle("#it{t} (ns)");
          outputContainer->Add(fhEnergyTimeTCardCorrNoSelection1Cell[tm]); 
          
          fhEnergyTimeTCardCorrNoSelectionExotic[tm]  = new TH2F 
          (Form("hEnergyTimeTCardCorrNoSelectionExotic%s",add[tm].Data()),
           Form("#it{t} vs #it{E},  exo > 0.97, %s",add[tm].Data()),
           nptbins,ptmin,ptmax,300,-150,150); 
          fhEnergyTimeTCardCorrNoSelectionExotic[tm]->SetXTitle("#it{E} (GeV)");
          fhEnergyTimeTCardCorrNoSelectionExotic[tm]->SetYTitle("#it{t} (ns)");
          outputContainer->Add(fhEnergyTimeTCardCorrNoSelectionExotic[tm]); 
          
          fhColRowExoticLowE1Cell[tm] = new TH2F
          (Form("hColRowExoticLowE1Cell%s",add[tm].Data()),
           Form("column vs row, 1 cell, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExoticLowE1Cell[tm]->SetYTitle("row");
          fhColRowExoticLowE1Cell[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExoticLowE1Cell[tm]) ;
          
          fhColRowExoticHighE1Cell[tm] = new TH2F
          (Form("hColRowExoticHighE1Cell%s",add[tm].Data()),
           Form("column vs row, 1 cell, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExoticHighE1Cell[tm]->SetYTitle("row");
          fhColRowExoticHighE1Cell[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExoticHighE1Cell[tm]) ;
          
          fhColRowExoticLowE[tm] = new TH2F
          (Form("hColRowExoticLowE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExoticLowE[tm]->SetYTitle("row");
          fhColRowExoticLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExoticLowE[tm]) ;
          
          fhColRowExoticHighE[tm] = new TH2F
          (Form("hColRowExoticHighE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExoticHighE[tm]->SetYTitle("row");
          fhColRowExoticHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExoticHighE[tm]) ;
          
          fhColRowExotic2ndCellDiffLowE[tm] = new TH2F
          (Form("hColRowExotic2ndCellDiffLowE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExotic2ndCellDiffLowE[tm]->SetYTitle("row");
          fhColRowExotic2ndCellDiffLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExotic2ndCellDiffLowE[tm]) ;
          
          fhColRowExotic2ndCellDiffHighE[tm] = new TH2F
          (Form("hColRowExotic2ndCellDiffHighE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExotic2ndCellDiffHighE[tm]->SetYTitle("row");
          fhColRowExotic2ndCellDiffHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExotic2ndCellDiffHighE[tm]) ;
          
          fhColRowExotic2ndCellSameLowE[tm] = new TH2F
          (Form("hColRowExotic2ndCellSameLowE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExotic2ndCellSameLowE[tm]->SetYTitle("row");
          fhColRowExotic2ndCellSameLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExotic2ndCellSameLowE[tm]) ;
          
          fhColRowExotic2ndCellSameHighE[tm] = new TH2F
          (Form("hColRowExotic2ndCellSameHighE%s",add[tm].Data()),
           Form("column vs row, max E cell, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowExotic2ndCellSameHighE[tm]->SetYTitle("row");
          fhColRowExotic2ndCellSameHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowExotic2ndCellSameHighE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExoticLowE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExoticLowE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExoticLowE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExoticLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExoticLowE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExoticHighE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExoticHighE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExoticHighE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExoticHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExoticHighE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellDiffLowE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellDiffLowE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellDiffHighE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellDiffHighE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellSameLowE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellSameLowE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellSameHighE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellSameHighE[tm]) ;        
          
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameLowE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellDiffNoSameHighE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, 5 < E < 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffLowE[tm]) ;
          
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[tm] = new TH2F
          (Form("hColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE%s",add[tm].Data()),
           Form("column vs row, max E cell for TCard correlation selected clusters, exo > 0.97, E > 8 GeV %s",add[tm].Data()),
           ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[tm]->SetYTitle("row");
          fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[tm]->SetXTitle("column");
          outputContainer->Add(fhColRowTCardCorrNoSelectionExotic2ndCellSameNoDiffHighE[tm]) ;
          
          fhNCellsTCardSameAndDiffFractionExotic[tm]  = new TH2F 
          (Form("hNCellsTCardSameAndDiffFraction_Exotic%s",add[tm].Data()),
           Form("#it{n}_{cells} same TCard vs diff TCard fraction, w > 0.01, exo > 0.97 %s",add[tm].Data()),
           nptbins,ptmin,ptmax,100,0,1);  
          fhNCellsTCardSameAndDiffFractionExotic[tm]->SetXTitle("#it{E} (GeV)");
          fhNCellsTCardSameAndDiffFractionExotic[tm]->SetYTitle("#it{n}_{cells} - same TCard / #it{n}_{cells} - total");
          outputContainer->Add(fhNCellsTCardSameAndDiffFractionExotic[tm]);   
          
          fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[tm] = new TH2F 
          (Form("hSameRowDiffColAndTCardCellsEnergyDiffClusterEExo%s",add[tm].Data()),
           Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} vs #it{E}_{cluster}, exo > 0.97 %s",add[tm].Data()),
           nptbins,ptmin,ptmax,200,-10,10); 
          fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
          fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[tm]->SetYTitle("#it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} (GeV)");
          outputContainer->Add(fhSameRowDiffColAndTCardCellsEnergyDiffClusterEExo[tm]); 
          
          fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo[tm] = new TH2F 
          (Form("hSameRowDiffColAndTCardCellsTimeDiffClusterEExo%s",add[tm].Data()),
           Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} vs #it{E}_{cluster}, exo > 0.97 %s",add[tm].Data()),
           nptbins,ptmin,ptmax,200,-100,100); 
          fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo[tm]->SetXTitle("#it{E}_{cluster} (GeV)");
          fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo[tm]->SetYTitle("#it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} (ns)");
          outputContainer->Add(fhSameRowDiffColAndTCardCellsTimeDiffClusterEExo[tm]); 
          
          fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[tm] = new TH2F 
          (Form("hSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo%s",add[tm].Data()),
           Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} vs #it{E}_{cell max}, exo > 0.97 %s",add[tm].Data()),
           nptbins,ptmin,ptmax,200,-10,10); 
          fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[tm]->SetXTitle("#it{E}_{cell max} (GeV)");
          fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[tm]->SetYTitle("#it{E}_{cell}^{same TCard}-#it{E}_{cell}^{diff TCard} (GeV)");
          outputContainer->Add(fhSameRowDiffColAndTCardCellsEnergyDiffCellMaxEExo[tm]); 
          
          fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo[tm] = new TH2F 
          (Form("hSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo%s",add[tm].Data()),
           Form("#Delta row = 0, |#Delta col = 1|, with respect to leading cell, #it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} vs #it{E}_{cell max}, exo > 0.97 %s",add[tm].Data()),
           nptbins,ptmin,ptmax,200,-100,100); 
          fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo[tm]->SetXTitle("#it{E}_{cell max} (GeV)");
          fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo[tm]->SetYTitle("#it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} (ns)");
          outputContainer->Add(fhSameRowDiffColAndTCardCellsTimeDiffCellMaxEExo[tm]); 
        }
        
        for(Int_t i = 0; i < 6; i++)
        { 
          for(Int_t j = 0; j < 6; j++)
          {
            fhLambda0TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda0TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{0} vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda0TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda0TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{0}");
            outputContainer->Add(fhLambda0TCardCorrelNCell[i][j][tm]); 
            
            fhLambda1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{1} vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda1TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{1}");
            outputContainer->Add(fhLambda1TCardCorrelNCell[i][j][tm]); 
            
            fhLambda0NLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda0NLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{0} vs #it{E}, nlm=1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda0NLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda0NLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{0}");
            outputContainer->Add(fhLambda0NLM1TCardCorrelNCell[i][j][tm]); 
            
            fhLambda1NLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda1NLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{1} vs #it{E}, nlm=1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda1NLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda1NLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{1}");
            outputContainer->Add(fhLambda1NLM1TCardCorrelNCell[i][j][tm]); 
            
            fhLambda0NLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda0NLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{0} vs #it{E}, nlm=2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda0NLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda0NLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{0}");
            outputContainer->Add(fhLambda0NLM2TCardCorrelNCell[i][j][tm]); 
            
            fhLambda1NLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLambda1NLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#lambda^{2}_{1} vs #it{E}, nlm=2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
            fhLambda1NLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhLambda1NLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{1}");
            outputContainer->Add(fhLambda1NLM2TCardCorrelNCell[i][j][tm]); 
            
            
            //          fhLambdaRTCardCorrelNCell[i][j][tm]  = new TH2F 
            //          (Form("hLambdaRTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
            //           Form("#lambda^{2}_{1}/#lambda^{2}_{0} vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
            //           nptbins,ptmin,ptmax,110,0,1.1); 
            //          fhLambdaRTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            //          fhLambdaRTCardCorrelNCell[i][j][tm]->SetYTitle("#lambda^{2}_{1}/#lambda^{2}_{0}");
            //          outputContainer->Add(fhLambdaRTCardCorrelNCell[i][j][tm]); 
            
            fhNLocMaxTCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hNLocMaxTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{n}_{LM} vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,10,0,10); 
            fhNLocMaxTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhNLocMaxTCardCorrelNCell[i][j][tm]->SetYTitle("#it{n}_{LM}");
            outputContainer->Add(fhNLocMaxTCardCorrelNCell[i][j][tm]); 
            
            fhEMaxRatNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hEMaxRatNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhEMaxRatNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhEMaxRatNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
            outputContainer->Add(fhEMaxRatNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhEMaxRatNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hEMaxRatNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhEMaxRatNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhEMaxRatNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
            outputContainer->Add(fhEMaxRatNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhEMaxRatNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hEMaxRatNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhEMaxRatNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhEMaxRatNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{max}/#it{E}_{cluster}");
            outputContainer->Add(fhEMaxRatNLM3TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndRatNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndRatNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndRatNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndRatNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cluster}");
            outputContainer->Add(fhE2ndRatNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndRatNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndRatNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndRatNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndRatNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
            outputContainer->Add(fhE2ndRatNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndRatNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndRatNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndRatNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndRatNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster}");
            outputContainer->Add(fhE2ndRatNLM3TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndEMaxRatNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndEMaxRatNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max} vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndEMaxRatNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndEMaxRatNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd max}/#it{E}_{cell}^{max}");
            outputContainer->Add(fhE2ndEMaxRatNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndEMaxRatNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndEMaxRatNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max} vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndEMaxRatNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndEMaxRatNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
            outputContainer->Add(fhE2ndEMaxRatNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhE2ndEMaxRatNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hE2ndEMaxRatNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}^{2nd loc max}/#it{E}_{cluster} vs #it{E}_{cell}^{max}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhE2ndEMaxRatNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhE2ndEMaxRatNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}^{2nd loc max}/#it{E}_{cell}^{max}");
            outputContainer->Add(fhE2ndEMaxRatNLM3TCardCorrelNCell[i][j][tm]); 
            
            
            fhECellClusRatNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellClusRatNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhECellClusRatNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhECellClusRatNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
            outputContainer->Add(fhECellClusRatNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhECellClusRatNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellClusRatNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhECellClusRatNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhECellClusRatNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
            outputContainer->Add(fhECellClusRatNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhECellClusRatNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellClusRatNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{E}_{cell}/#it{E}_{cluster} vs #it{E}_{cluster}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,100,0,1); 
            fhECellClusRatNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhECellClusRatNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
            outputContainer->Add(fhECellClusRatNLM3TCardCorrelNCell[i][j][tm]); 
            
            fhLogECellNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLogECellNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("log(#it{E}_{cell}) vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,150,-3,3); 
            fhLogECellNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhLogECellNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("log(#it{E}_{cell})");
            outputContainer->Add(fhLogECellNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhLogECellNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLogECellNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("log(#it{E}_{cell}) vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,150,-3,3); 
            fhLogECellNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhLogECellNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("log(#it{E}_{cell})");
            outputContainer->Add(fhLogECellNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhLogECellNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hLogECellNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("log(#it{E}_{cell}) vs #it{E}_{cluster}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,150,-3,3); 
            fhLogECellNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhLogECellNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("log(#it{E}_{cell})");
            outputContainer->Add(fhLogECellNLM3TCardCorrelNCell[i][j][tm]); 
            
            fhECellWeightNLM1TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellWeightNLM1TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs #it{E}_{cluster}, #it{n}_{LM} = 1, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,90,0,4.5); 
            fhECellWeightNLM1TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhECellWeightNLM1TCardCorrelNCell[i][j][tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
            outputContainer->Add(fhECellWeightNLM1TCardCorrelNCell[i][j][tm]); 
            
            fhECellWeightNLM2TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellWeightNLM2TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs #it{E}_{cluster}, #it{n}_{LM} = 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,90,0,4.5); 
            fhECellWeightNLM2TCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhECellWeightNLM2TCardCorrelNCell[i][j][tm]->SetYTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
            outputContainer->Add(fhECellWeightNLM2TCardCorrelNCell[i][j][tm]); 
            
            fhECellWeightNLM3TCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hECellWeightNLM3TCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster})) vs #it{E}_{cluster}, #it{n}_{LM} > 2, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,90,0,4.5); 
            fhECellWeightNLM3TCardCorrelNCell[i][j][tm]->SetXTitle("#it{w}=Max(4,5+log(#it{E}_{cell}/#it{E}_{cluster}))");
            fhECellWeightNLM3TCardCorrelNCell[i][j][tm]->SetYTitle("Log. weight");
            outputContainer->Add(fhECellWeightNLM3TCardCorrelNCell[i][j][tm]); 
            
            fhMassEClusTCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hMassEClusTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{M}_{#gamma #gamma} vs #it{E}_{cluster}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
            fhMassEClusTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{cluster} (GeV)");
            fhMassEClusTCardCorrelNCell[i][j][tm]->SetYTitle("#it{M}_{#gamma #gamma}");
            outputContainer->Add(fhMassEClusTCardCorrelNCell[i][j][tm]);         
            
            //          fhMassEPairTCardCorrelNCell[i][j][tm]  = new TH2F 
            //          (Form("hMassEPairTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
            //           Form("#it{M}_{#gamma #gamma} vs #it{E}_{pair}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
            //           nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
            //          fhMassEPairTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E}_{pair} (GeV)");
            //          fhMassEPairTCardCorrelNCell[i][j][tm]->SetYTitle("#it{M}_{#gamma #gamma}");
            //          outputContainer->Add(fhMassEPairTCardCorrelNCell[i][j][tm]); 
            
            fhTimeDiffTCardCorrelNCell[i][j][tm]  = new TH2F 
            (Form("hTimeDiffTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("#it{t}_{cell}^{max}-#it{t}_{cell}^{other} vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             nptbins,ptmin,ptmax,300,-150,150); 
            fhTimeDiffTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
            fhTimeDiffTCardCorrelNCell[i][j][tm]->SetYTitle("#it{t}_{cell}^{max}-#it{t}_{cell}^{other}");
            outputContainer->Add(fhTimeDiffTCardCorrelNCell[i][j][tm]); 
            
            fhColRowTCardCorrelNCellLowE[i][j][tm] = new TH2F
            (Form("hColRowTCardCorrelNCellLowE_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("column vs row, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
            fhColRowTCardCorrelNCellLowE[i][j][tm]->SetYTitle("row");
            fhColRowTCardCorrelNCellLowE[i][j][tm]->SetXTitle("column");
            outputContainer->Add(fhColRowTCardCorrelNCellLowE[i][j][tm]) ;
            
            fhColRowTCardCorrelNCellHighE[i][j][tm] = new TH2F
            (Form("hColRowTCardCorrelNCellHighE_Same%d_Diff%d%s",i,j,add[tm].Data()),
             Form("column vs row,N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
             ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
            fhColRowTCardCorrelNCellHighE[i][j][tm]->SetYTitle("row");
            fhColRowTCardCorrelNCellHighE[i][j][tm]->SetXTitle("column");
            outputContainer->Add(fhColRowTCardCorrelNCellHighE[i][j][tm]) ;
            
            if(fStudyExotic)
            {
              fhExoticTCardCorrelNCell[i][j][tm]  = new TH2F 
              (Form("hExoticTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
               Form("exoticity vs #it{E}, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
               nptbins,ptmin,ptmax,200,-1,1); 
              fhExoticTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
              fhExoticTCardCorrelNCell[i][j][tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
              outputContainer->Add(fhExoticTCardCorrelNCell[i][j][tm]); 
              
              fhTimeDiffExoTCardCorrelNCell[i][j][tm]  = new TH2F 
              (Form("hTimeDiffExoTCardCorrelNCell_Same%d_Diff%d%s",i,j,add[tm].Data()),
               Form("#it{t}_{cell}^{max}-#it{t}_{cell}^{other} vs #it{E}, N cells with  w > 0.01, exoticity > 0.97, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
               nptbins,ptmin,ptmax,300,-150,150); 
              fhTimeDiffExoTCardCorrelNCell[i][j][tm]->SetXTitle("#it{E} (GeV)");
              fhTimeDiffExoTCardCorrelNCell[i][j][tm]->SetYTitle("#it{t}_{cell}^{max}-#it{t}_{cell}^{other}");
              outputContainer->Add(fhTimeDiffExoTCardCorrelNCell[i][j][tm]); 
              
              fhColRowTCardCorrelNCellExoticLowE[i][j][tm] = new TH2F
              (Form("hColRowTCardCorrelNCellExoticLowE_Same%d_Diff%d%s",i,j,add[tm].Data()),
               Form("column vs row, N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
               ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
              fhColRowTCardCorrelNCellExoticLowE[i][j][tm]->SetYTitle("row");
              fhColRowTCardCorrelNCellExoticLowE[i][j][tm]->SetXTitle("column");
              outputContainer->Add(fhColRowTCardCorrelNCellExoticLowE[i][j][tm]) ;
              
              fhColRowTCardCorrelNCellExoticHighE[i][j][tm] = new TH2F
              (Form("hColRowTCardCorrelNCellExoticHighE_Same%d_Diff%d%s",i,j,add[tm].Data()),
               Form("column vs row,N cells with  w > 0.01, TCard same = %d, diff =%d %s",i,j,add[tm].Data()),
               ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
              fhColRowTCardCorrelNCellExoticHighE[i][j][tm]->SetYTitle("row");
              fhColRowTCardCorrelNCellExoticHighE[i][j][tm]->SetXTitle("column");
              outputContainer->Add(fhColRowTCardCorrelNCellExoticHighE[i][j][tm]) ;
            }
          }
          
          ///////////
          ///////////
          
          //        fhLambda0TCardCorrelN[i][tm]  = new TH2F 
          //        (Form("hLambda0TCardCorrelN_Case%d%s",i,add[tm].Data()),
          //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, N corr = %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          //        fhLambda0TCardCorrelN[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhLambda0TCardCorrelN[i][tm]->SetYTitle("#lambda^{2}_{0}");
          //        outputContainer->Add(fhLambda0TCardCorrelN[i][tm]); 
          //        
          //        fhNCellsTCardCorrelN[i][tm]  = new TH2F 
          //        (Form("hNCellsTCardCorrelN_Case%d%s",i,add[tm].Data()),
          //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, N corr = %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
          //        fhNCellsTCardCorrelN[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhNCellsTCardCorrelN[i][tm]->SetYTitle("#it{n}_{cells}");
          //        outputContainer->Add(fhNCellsTCardCorrelN[i][tm]);
          //        
          //        fhExoticTCardCorrelN[i][tm]  = new TH2F 
          //        (Form("hExoticTCardCorrelN_Case%d%s",i,add[tm].Data()),
          //         Form("exoticity vs #it{E}, max E cell correl with TCard cell, N corr = %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,200,-1,1); 
          //        fhExoticTCardCorrelN[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhExoticTCardCorrelN[i][tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          //        outputContainer->Add(fhExoticTCardCorrelN[i][tm]); 
          //        
          //        fhColRowTCardCorrelNLowE[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNLowE_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, E > 2 GeV, N corr = %d %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNLowE[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNLowE[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNLowE[i][tm]) ;
          //        
          //        fhColRowTCardCorrelNHighE[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNHighE_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, E > 8 GeV, N corr = %d %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNHighE[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNHighE[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNHighE[i][tm]) ;
          //      
          //        ////////
          //        ////////
          //                
          //        fhLambda0TCardCorrelNExotic[i][tm]  = new TH2F 
          //        (Form("hLambda0TCardCorrelN_Exotic_Case%d%s",i,add[tm].Data()),
          //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, exo > 0.97, N corr = %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          //        fhLambda0TCardCorrelNExotic[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhLambda0TCardCorrelNExotic[i][tm]->SetYTitle("#lambda^{2}_{0}");
          //        outputContainer->Add(fhLambda0TCardCorrelNExotic[i][tm]); 
          //        
          //        fhNCellsTCardCorrelNExotic[i][tm]  = new TH2F 
          //        (Form("hNCellsTCardCorrelN_Exotic_Case%d%s",i,add[tm].Data()),
          //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, exo > 0.97, N corr = %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
          //        fhNCellsTCardCorrelNExotic[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhNCellsTCardCorrelNExotic[i][tm]->SetYTitle("#it{n}_{cells}");
          //        outputContainer->Add(fhNCellsTCardCorrelNExotic[i][tm]);
          //        
          //        fhColRowTCardCorrelNLowEExotic[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNLowEExotic_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, exo > 0.97, E > 2 GeV, N corr = %d %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNLowEExotic[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNLowEExotic[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNLowEExotic[i][tm]) ;
          //        
          //        fhColRowTCardCorrelNHighEExotic[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNHighEExotic_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, exo > 0.97, E > 8 GeV, N corr = %d %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNHighEExotic[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNHighEExotic[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNHighEExotic[i][tm]) ;
          //
          //        ///////////
          //        ///////////
          //        
          //        fhLambda0TCardCorrelNAllSameTCard[i][tm]  = new TH2F 
          //        (Form("hLambda0TCardCorrelNAllSameTCard_Case%d%s",i,add[tm].Data()),
          //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          //        fhLambda0TCardCorrelNAllSameTCard[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhLambda0TCardCorrelNAllSameTCard[i][tm]->SetYTitle("#lambda^{2}_{0}");
          //        outputContainer->Add(fhLambda0TCardCorrelNAllSameTCard[i][tm]); 
          //        
          //        fhNCellsTCardCorrelNAllSameTCard[i][tm]  = new TH2F 
          //        (Form("hNCellsTCardCorrelNAllSameTCard_Case%d%s",i,add[tm].Data()),
          //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
          //        fhNCellsTCardCorrelNAllSameTCard[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhNCellsTCardCorrelNAllSameTCard[i][tm]->SetYTitle("#it{n}_{cells}");
          //        outputContainer->Add(fhNCellsTCardCorrelNAllSameTCard[i][tm]);
          //        
          //        fhExoticTCardCorrelNAllSameTCard[i][tm]  = new TH2F 
          //        (Form("hExoticTCardCorrelNAllSameTCard_Case%d%s",i,add[tm].Data()),
          //         Form("exoticity vs #it{E}, max E cell correl with TCard cell, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,200,-1,1); 
          //        fhExoticTCardCorrelNAllSameTCard[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhExoticTCardCorrelNAllSameTCard[i][tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          //        outputContainer->Add(fhExoticTCardCorrelNAllSameTCard[i][tm]); 
          //                
          //        fhColRowTCardCorrelNAllSameTCardLowE[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNAllSameTCardLowE_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, E > 2 GeV, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNAllSameTCardLowE[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNAllSameTCardLowE[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNAllSameTCardLowE[i][tm]) ;
          //        
          //        fhColRowTCardCorrelNAllSameTCardHighE[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNAllSameTCardHighE_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, E > 8 GeV, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNAllSameTCardHighE[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNAllSameTCardHighE[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNAllSameTCardHighE[i][tm]) ;
          //        
          //        ////////
          //        
          //        fhLambda0TCardCorrelNAllSameTCardExotic[i][tm]  = new TH2F 
          //        (Form("hLambda0TCardCorrelNAllSameTCard_Exotic_Case%d%s",i,add[tm].Data()),
          //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, exo > 0.97, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
          //        fhLambda0TCardCorrelNAllSameTCardExotic[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhLambda0TCardCorrelNAllSameTCardExotic[i][tm]->SetYTitle("#lambda^{2}_{0}");
          //        outputContainer->Add(fhLambda0TCardCorrelNAllSameTCardExotic[i][tm]); 
          //        
          //        fhNCellsTCardCorrelNAllSameTCardExotic[i][tm]  = new TH2F 
          //        (Form("hNCellsTCardCorrelNAllSameTCard_Exotic_Case%d%s",i,add[tm].Data()),
          //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, exo > 0.97, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
          //        fhNCellsTCardCorrelNAllSameTCardExotic[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhNCellsTCardCorrelNAllSameTCardExotic[i][tm]->SetYTitle("#it{n}_{cells}");
          //        outputContainer->Add(fhNCellsTCardCorrelNAllSameTCardExotic[i][tm]);
          //        
          //        fhColRowTCardCorrelNAllSameTCardLowEExotic[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNAllSameTCardLowEExotic_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, exo > 0.97, E > 2 GeV, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNAllSameTCardLowEExotic[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNAllSameTCardLowEExotic[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNAllSameTCardLowEExotic[i][tm]) ;
          //        
          //        fhColRowTCardCorrelNAllSameTCardHighEExotic[i][tm] = new TH2F
          //        (Form("hColRowTCardCorrelNAllSameTCardHighEExotic_Case%d%s",i,add[tm].Data()),
          //         Form("column vs row, max E cell correl with TCard cell, exo > 0.97, E > 8 GeV, N corr = %d, no other TCard cells %s",i,add[tm].Data()),
          //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
          //        fhColRowTCardCorrelNAllSameTCardHighEExotic[i][tm]->SetYTitle("row");
          //        fhColRowTCardCorrelNAllSameTCardHighEExotic[i][tm]->SetXTitle("column");
          //        outputContainer->Add(fhColRowTCardCorrelNAllSameTCardHighEExotic[i][tm]) ;
        }
        
        //      for(Int_t i = 0; i < 7; i++)
        //      {
        //        fhLambda0TCardCorrel[i][tm]  = new TH2F 
        //        (Form("hLambda0TCardCorrel_Case%d%s",i,add[tm].Data()),
        //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        //        fhLambda0TCardCorrel[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhLambda0TCardCorrel[i][tm]->SetYTitle("#lambda^{2}_{0}");
        //        outputContainer->Add(fhLambda0TCardCorrel[i][tm]); 
        //        
        //        fhNCellsTCardCorrel[i][tm]  = new TH2F 
        //        (Form("hNCellsTCardCorrel_Case%d%s",i,add[tm].Data()),
        //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        //        fhNCellsTCardCorrel[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhNCellsTCardCorrel[i][tm]->SetYTitle("#it{n}_{cells}");
        //        outputContainer->Add(fhNCellsTCardCorrel[i][tm]);
        //        
        //        fhExoticTCardCorrel[i][tm]  = new TH2F 
        //        (Form("hExoticTCardCorrel_Case%d%s",i,add[tm].Data()),
        //         Form("exoticity vs #it{E}, max E cell correl with TCard cell, N corr = %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,200,-1,1); 
        //        fhExoticTCardCorrel[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhExoticTCardCorrel[i][tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
        //        outputContainer->Add(fhExoticTCardCorrel[i][tm]); 
        //      }
        //      
        //      for(Int_t i = 0; i < 4; i++)
        //      {
        //        fhLambda0TCardCorrelExotic[i][tm]  = new TH2F 
        //        (Form("hLambda0TCardCorrel_Exotic_Case%d%s",i,add[tm].Data()),
        //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, exo>0.97, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        //        fhLambda0TCardCorrelExotic[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhLambda0TCardCorrelExotic[i][tm]->SetYTitle("#lambda^{2}_{0}");
        //        outputContainer->Add(fhLambda0TCardCorrelExotic[i][tm]); 
        //        
        //        fhNCellsTCardCorrelExotic[i][tm]  = new TH2F 
        //        (Form("hNCellsTCardCorrel_Exotic_Case%d%s",i,add[tm].Data()),
        //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, exot > 0.97,case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        //        fhNCellsTCardCorrelExotic[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhNCellsTCardCorrelExotic[i][tm]->SetYTitle("#it{n}_{cells}");
        //        outputContainer->Add(fhNCellsTCardCorrelExotic[i][tm]);
        //      }
        
        
        for(Int_t i = 0; i < fNEBinCuts; i++)
        {
          if(fStudyExotic)
          {
            fhLambda0Exoticity[i][tm]  = new TH2F 
            (Form("hLambda0Exoticity_EBin%d%s",i,add[tm].Data()),
             Form("#lambda^{2}_{0} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
             200,-1,1,ssbins,ssmin,ssmax); 
            fhLambda0Exoticity[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            fhLambda0Exoticity[i][tm]->SetYTitle("#lambda^{2}_{0}");
            outputContainer->Add(fhLambda0Exoticity[i][tm]);    
            
            fhLambda1Exoticity[i][tm]  = new TH2F 
            (Form("hLambda1Exoticity_EBin%d%s",i,add[tm].Data()),
             Form("#lambda^{2}_{1} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
             200,-1,1,ssbins,ssmin,ssmax); 
            fhLambda1Exoticity[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            fhLambda1Exoticity[i][tm]->SetYTitle("#lambda^{2}_{1}");
            outputContainer->Add(fhLambda1Exoticity[i][tm]);    
            
            //        fhLambdaRExoticity[i][tm]  = new TH2F 
            //        (Form("hLambdaRExoticity_EBin%d%s",i,add[tm].Data()),
            //         Form("#lambda^{2}_{1}/#lambda^{2}_{0} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
            //         200,-1,1,110,0,1.1); 
            //        fhLambdaRExoticity[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            //        fhLambdaRExoticity[i][tm]->SetYTitle("#lambda^{2}_{1}/#lambda^{2}_{0}");
            //        outputContainer->Add(fhLambdaRExoticity[i][tm]);    
            
            fhNCellsExoticity[i][tm]  = new TH2F 
            (Form("hNCellsExoticity_EBin%d%s",i,add[tm].Data()),
             Form("#it{n}_{cells} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
             200,-1,1,nceclbins,nceclmin,nceclmax); 
            fhNCellsExoticity[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            fhNCellsExoticity[i][tm]->SetYTitle("#it{n}_{cells}");
            outputContainer->Add(fhNCellsExoticity[i][tm]); 
            
            fhTimeExoticity[i][tm]  = new TH2F 
            (Form("hTimeExoticity_EBin%d%s",i,add[tm].Data()),
             Form("#it{t} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
             200,-1,1,100,-25,25); 
            fhTimeExoticity[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            fhTimeExoticity[i][tm]->SetYTitle("#it{t} (ns)");
            outputContainer->Add(fhTimeExoticity[i][tm]); 
            
            fhNCellsTCardSameAndDiffExotic[i][tm]  = new TH2F 
            (Form("hNCellsTCardSameAndDiff_Exotic_EBin%d%s",i,add[tm].Data()),
             Form("#it{n}_{cells} same TCard vs diff TCard, w > 0.01, %2.2f<#it{E}<%2.2f GeV, exo > 0.97 %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
             nceclbins,nceclmin,nceclmax,nceclbins,nceclmin,nceclmax); 
            fhNCellsTCardSameAndDiffExotic[i][tm]->SetXTitle("#it{n}_{cells} - diff TCard");
            fhNCellsTCardSameAndDiffExotic[i][tm]->SetYTitle("#it{n}_{cells} - same TCard");
            outputContainer->Add(fhNCellsTCardSameAndDiffExotic[i][tm]);            
            
            //        fhLambda0ExoticityAllSameTCard[i][tm]  = new TH2F 
            //        (Form("hLambda0ExoticityAllSameTCard_EBin%d%s",i,add[tm].Data()),
            //         Form("#lambda^{2}_{0} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
            //         200,-1,1,ssbins,ssmin,ssmax); 
            //        fhLambda0ExoticityAllSameTCard[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            //        fhLambda0ExoticityAllSameTCard[i][tm]->SetYTitle("#lambda^{2}_{0}");
            //        outputContainer->Add(fhLambda0ExoticityAllSameTCard[i][tm]);    
            //
            //        fhLambda1ExoticityAllSameTCard[i][tm]  = new TH2F 
            //        (Form("hLambda1ExoticityAllSameTCard_EBin%d%s",i,add[tm].Data()),
            //         Form("#lambda^{2}_{1} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
            //         200,-1,1,ssbins,ssmin,ssmax); 
            //        fhLambda1ExoticityAllSameTCard[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            //        fhLambda1ExoticityAllSameTCard[i][tm]->SetYTitle("#lambda^{2}_{1}");
            //        outputContainer->Add(fhLambda1ExoticityAllSameTCard[i][tm]);    
            //        
            //        fhLambdaRExoticityAllSameTCard[i][tm]  = new TH2F 
            //        (Form("hLambdaRExoticityAllSameTCard_EBin%d%s",i,add[tm].Data()),
            //         Form("#lambda^{2}_{1}/#lambda^{2}_{0} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
            //         200,-1,1,110,0,1.1); 
            //        fhLambdaRExoticityAllSameTCard[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            //        fhLambdaRExoticityAllSameTCard[i][tm]->SetYTitle("#lambda^{2}_{1}/#lambda^{2}_{0}");
            //        outputContainer->Add(fhLambdaRExoticityAllSameTCard[i][tm]);    
            //        
            //        fhNCellsExoticityAllSameTCard[i][tm]  = new TH2F 
            //        (Form("hNCellsExoticityAllSameTCard_EBin%d%s",i,add[tm].Data()),
            //         Form("#it{n}_{cells} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
            //         200,-1,1,nceclbins,nceclmin,nceclmax); 
            //        fhNCellsExoticityAllSameTCard[i][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
            //        fhNCellsExoticityAllSameTCard[i][tm]->SetYTitle("#it{n}_{cells}");
            //        outputContainer->Add(fhNCellsExoticityAllSameTCard[i][tm]); 
            //        
          }
          
          fhNCellsTCardSameAndDiff[i][tm]  = new TH2F 
          (Form("hNCellsTCardSameAndDiff_EBin%d%s",i,add[tm].Data()),
           Form("#it{n}_{cells} same TCard vs diff TCard, w > 0.01, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
           nceclbins,nceclmin,nceclmax,nceclbins,nceclmin,nceclmax); 
          fhNCellsTCardSameAndDiff[i][tm]->SetXTitle("#it{n}_{cells} - diff TCard");
          fhNCellsTCardSameAndDiff[i][tm]->SetYTitle("#it{n}_{cells} - same TCard");
          outputContainer->Add(fhNCellsTCardSameAndDiff[i][tm]);       
          
          fhLambda0Lambda1[i][tm]  = new TH2F 
          (Form("hLambda0Lambda1_EBin%d%s",i,add[tm].Data()),
           Form("#lambda^{2}_{0} vs #lambda^{2}_{1}, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
           ssbins,ssmin,ssmax,ssbins,ssmin,ssmax); 
          fhLambda0Lambda1[i][tm]->SetXTitle("#lambda^{2}_{1}");
          fhLambda0Lambda1[i][tm]->SetYTitle("#lambda^{2}_{0}");
          outputContainer->Add(fhLambda0Lambda1[i][tm]);   
          
          //        fhLambda0Lambda1AllSameTCard[i][tm]  = new TH2F 
          //        (Form("hLambda0Lambda1AllSameTCard_EBin%d%s",i,add[tm].Data()),
          //         Form("#lambda^{2}_{0} vs #lambda^{2}_{1}, , all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV %s",fEBinCuts[i],fEBinCuts[i+1],add[tm].Data()),
          //         ssbins,ssmin,ssmax,ssbins,ssmin,ssmax); 
          //        fhLambda0Lambda1AllSameTCard[i][tm]->SetXTitle("#lambda^{2}_{1}");
          //        fhLambda0Lambda1AllSameTCard[i][tm]->SetYTitle("#lambda^{2}_{0}");
          //        outputContainer->Add(fhLambda0Lambda1AllSameTCard[i][tm]);   
        }
        
        if(fStudyExotic)
        {
          for(Int_t j = 0; j < 6; j++)
          {
            for(Int_t k = 0; k < 6; k++)
            {
              fhLambda0ExoticityPerNCell[j][k][tm]  = new TH2F 
              (Form("hLambda0Exoticity_NCell_Same%d_Diff%d%s",j,k,add[tm].Data()),
               Form("#lambda^{2}_{0} vs #it{exoticity}, #it{n}_{cell} TCard same = %d, diff =%d, #it{E}>8 GeV %s",j,k,add[tm].Data()),
               200,-1,1,ssbins,ssmin,ssmax); 
              fhLambda0ExoticityPerNCell[j][k][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
              fhLambda0ExoticityPerNCell[j][k][tm]->SetYTitle("#lambda^{2}_{0}");
              outputContainer->Add(fhLambda0ExoticityPerNCell[j][k][tm]);      
              
              fhLambda1ExoticityPerNCell[j][k][tm]  = new TH2F 
              (Form("hLambda1Exoticity_NCell_Same%d_Diff%d%s",j,k,add[tm].Data()),
               Form("#lambda^{2}_{1} vs #it{exoticity}, #it{n}_{cell} TCard same = %d, diff =%d, #it{E}>8 GeV %s",j,k,add[tm].Data()),
               200,-1,1,ssbins,ssmin,ssmax); 
              fhLambda1ExoticityPerNCell[j][k][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
              fhLambda1ExoticityPerNCell[j][k][tm]->SetYTitle("#lambda^{2}_{1}");
              outputContainer->Add(fhLambda1ExoticityPerNCell[j][k][tm]);  
              
              //          fhLambdaRExoticityPerNCell[j][k][tm]  = new TH2F 
              //          (Form("hLambdaRExoticity_NCell_Same%d_Diff%d%s",j,k,add[tm].Data()),
              //           Form("#lambda^{2}_{1}/#lambda^{2}_{0} vs #it{exoticity}, #it{n}_{cell} TCard same = %d, diff =%d, #it{E}>8 GeV %s",j,k,add[tm].Data()),
              //           200,-1,1,110,0,1.1); 
              //          fhLambdaRExoticityPerNCell[j][k][tm]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
              //          fhLambdaRExoticityPerNCell[j][k][tm]->SetYTitle("#lambda^{2}_{1}/#lambda^{2}_{0}");
              //          outputContainer->Add(fhLambdaRExoticityPerNCell[j][k][tm]);  
            }
          }
        }
        
        //      for(Int_t i = 0; i < 6; i++)
        //      {
        //        fhLambda0TCardCorrelNearRow[i][tm]  = new TH2F 
        //        (Form("hLambda0TCardCorrelNearRow_Case%d%s",i,add[tm].Data()),
        //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with TCard cell, one TCard cell is 1 row away, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        //        fhLambda0TCardCorrelNearRow[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhLambda0TCardCorrelNearRow[i][tm]->SetYTitle("#lambda^{2}_{0}");
        //        outputContainer->Add(fhLambda0TCardCorrelNearRow[i][tm]); 
        //        
        //        fhNCellsTCardCorrelNearRow[i][tm]  = new TH2F 
        //        (Form("hNCellsTCardCorrelNearRow_Case%d%s",i,add[tm].Data()),
        //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with TCard cell, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        //        fhNCellsTCardCorrelNearRow[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhNCellsTCardCorrelNearRow[i][tm]->SetYTitle("#it{n}_{cells}");
        //        outputContainer->Add(fhNCellsTCardCorrelNearRow[i][tm]);
        //      }
        //      
        //      for(Int_t i = 0; i < 4; i++)
        //      {
        //        fhLambda0TCardCorrel2ndMax[i][tm]  = new TH2F 
        //        (Form("hLambda0TCardCorrel2ndMax_Case%d%s",i,add[tm].Data()),
        //         Form("#lambda^{2}_{0} vs #it{E}, max E cell correl with 2nd max TCard cell, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        //        fhLambda0TCardCorrel2ndMax[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhLambda0TCardCorrel2ndMax[i][tm]->SetYTitle("#lambda^{2}_{0}");
        //        outputContainer->Add(fhLambda0TCardCorrel2ndMax[i][tm]); 
        //        
        //        fhNCellsTCardCorrel2ndMax[i][tm]  = new TH2F 
        //        (Form("hNCellsTCardCorrel2ndMax_Case%d%s",i,add[tm].Data()),
        //         Form("custer # cells vs #it{E}, w > 0.01, max E cell correl with 2nd max TCard cell, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        //        fhNCellsTCardCorrel2ndMax[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhNCellsTCardCorrel2ndMax[i][tm]->SetYTitle("#it{n}_{cells}");
        //        outputContainer->Add(fhNCellsTCardCorrel2ndMax[i][tm]);
        //      }
        //      
        //      for(Int_t i = 0; i < 7; i++)
        //      {
        //        fhLambda0TCardCorrelOtherTCard[i][tm]  = new TH2F 
        //        (Form("hLambda0TCardCorrelOtherTCard_Case%d%s",i,add[tm].Data()),
        //         Form("#lambda^{2}_{0} vs #it{E}, correlation of cells in different TCards, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        //        fhLambda0TCardCorrelOtherTCard[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhLambda0TCardCorrelOtherTCard[i][tm]->SetYTitle("#lambda^{2}_{0}");
        //        outputContainer->Add(fhLambda0TCardCorrelOtherTCard[i][tm]); 
        //        
        //        fhNCellsTCardCorrelOtherTCard[i][tm]  = new TH2F 
        //        (Form("hNCellsTCardCorrelOtherTCard_Case%d%s",i,add[tm].Data()),
        //         Form("custer # cells vs #it{E}, w > 0.01, correlation of cells in different TCards, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
        //        fhNCellsTCardCorrelOtherTCard[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhNCellsTCardCorrelOtherTCard[i][tm]->SetYTitle("#it{n}_{cells}");
        //        outputContainer->Add(fhNCellsTCardCorrelOtherTCard[i][tm]);
        //        
        //        fhExoticTCardCorrelOtherTCard[i][tm]  = new TH2F 
        //        (Form("hExoticTCardCorrelOtherTCard_Case%d%s",i,add[tm].Data()),
        //         Form("exoticity vs #it{E}, w > 0.01, correlation of cells in different TCards, case %d %s",i,add[tm].Data()),
        //         nptbins,ptmin,ptmax,200,-1,1); 
        //        fhExoticTCardCorrelOtherTCard[i][tm]->SetXTitle("#it{E} (GeV)");
        //        fhExoticTCardCorrelOtherTCard[i][tm]->SetYTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
        //        outputContainer->Add(fhExoticTCardCorrelOtherTCard[i][tm]); 
        //                
        //        fhColRowTCardCorrelOtherTCardLowE[i][tm] = new TH2F
        //        (Form("hColRowTCardCorrelOtherTCardLowE_Case%d%s",i,add[tm].Data()),
        //         Form("column vs row for different 2 TCard correlation cases, E > 2 GeV, case %d %s",i,add[tm].Data()),
        //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        //        fhColRowTCardCorrelOtherTCardLowE[i][tm]->SetYTitle("row");
        //        fhColRowTCardCorrelOtherTCardLowE[i][tm]->SetXTitle("column");
        //        outputContainer->Add(fhColRowTCardCorrelOtherTCardLowE[i][tm]) ;
        //
        //        fhColRowTCardCorrelOtherTCardHighE[i][tm] = new TH2F
        //        (Form("hColRowTCardCorrelOtherTCardHighE_Case%d%s",i,add[tm].Data()),
        //         Form("column vs row for different 2 TCard correlation cases, E > 8 GeV, case %d %s",i,add[tm].Data()),
        //         ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
        //        fhColRowTCardCorrelOtherTCardHighE[i][tm]->SetYTitle("row");
        //        fhColRowTCardCorrelOtherTCardHighE[i][tm]->SetXTitle("column");
        //        outputContainer->Add(fhColRowTCardCorrelOtherTCardHighE[i][tm]) ;
        //      }
        
        for(Int_t i = 0; i < 12; i++)
        {
          fhTCardCorrECellMaxDiff[i][tm]  = new TH2F 
          (Form("hTCardCorrECellMaxDiff_Case%d%s",i,add[tm].Data()),
           Form("#it{E}_{cell}^{max}-#it{E}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, case %d %s",i,add[tm].Data()),
           nptbins,ptmin,ptmax,210,-1,20); 
          fhTCardCorrECellMaxDiff[i][tm]->SetXTitle("#it{E} (GeV)");
          fhTCardCorrECellMaxDiff[i][tm]->SetYTitle("#it{E}_{cell}^{max}-#it{E}_{cell} (GeV)");
          outputContainer->Add(fhTCardCorrECellMaxDiff[i][tm]); 
          
          fhTCardCorrEClusterDiff[i][tm]  = new TH2F 
          (Form("hTCardCorrEClusterDiff_Case%d%s",i,add[tm].Data()),
           Form("#it{E}_{cluster}-#it{E}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, case %d %s",i,add[tm].Data()),
           nptbins,ptmin,ptmax,210,-1,20); 
          fhTCardCorrEClusterDiff[i][tm]->SetXTitle("#it{E} (GeV)");
          fhTCardCorrEClusterDiff[i][tm]->SetYTitle("#it{E}_{cluster}-#it{E}_{cell} (GeV)");
          outputContainer->Add(fhTCardCorrEClusterDiff[i][tm]); 
          
          //        fhTCardCorrECellMaxRat[i][tm]  = new TH2F 
          //        (Form("hTCardCorrECellMaxRat_Case%d%s",i,add[tm].Data()),
          //         Form("#it{E}_{cell}/#it{E}_{cell}^{max} vs #it{E}_{cluster}, for (un)correlated cells in TCard, case %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,110,0,1.1); 
          //        fhTCardCorrECellMaxRat[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhTCardCorrECellMaxRat[i][tm]->SetYTitle("#it{E}_{cell}/#it{E}^{max}_{cell}");
          //        outputContainer->Add(fhTCardCorrECellMaxRat[i][tm]); 
          //        
          //        fhTCardCorrEClusterRat[i][tm]  = new TH2F 
          //        (Form("hTCardCorrEClusterRat_Case%d%s",i,add[tm].Data()),
          //         Form("#it{E}_{cell}/#it{E}_{cluster} vs #it{E}_{cluster}, for (un)correlated cells in TCard, case %d %s",i,add[tm].Data()),
          //         nptbins,ptmin,ptmax,110,0,1.1); 
          //        fhTCardCorrEClusterRat[i][tm]->SetXTitle("#it{E} (GeV)");
          //        fhTCardCorrEClusterRat[i][tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
          //        outputContainer->Add(fhTCardCorrEClusterRat[i][tm]); 
          
          fhTCardCorrTCellMaxDiff[i][tm]  = new TH2F 
          (Form("hTCardCorrTCellMaxDiff_Case%d%s",i,add[tm].Data()),
           Form("#it{t}_{cell}^{max}-#it{t}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, case %d %s",i,add[tm].Data()),
           nptbins,ptmin,ptmax,1000,-100,100); 
          fhTCardCorrTCellMaxDiff[i][tm]->SetXTitle("#it{E} (GeV)");
          fhTCardCorrTCellMaxDiff[i][tm]->SetYTitle("#it{t}_{cell}^{max}-#it{t}_{cell} (ns)");
          outputContainer->Add(fhTCardCorrTCellMaxDiff[i][tm]); 
          
          if(fStudyExotic)
          {
            fhTCardCorrECellMaxDiffExo[i][tm]  = new TH2F 
            (Form("hTCardCorrECellMaxDiffExo_Case%d%s",i,add[tm].Data()),
             Form("#it{E}_{cell}^{max}-#it{E}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, exoticity > 0.97, case %d %s",i,add[tm].Data()),
             nptbins,ptmin,ptmax,210,-1,20); 
            fhTCardCorrECellMaxDiffExo[i][tm]->SetXTitle("#it{E} (GeV)");
            fhTCardCorrECellMaxDiffExo[i][tm]->SetYTitle("#it{E}_{cell}^{max}-#it{E}_{cell} (GeV)");
            outputContainer->Add(fhTCardCorrECellMaxDiffExo[i][tm]); 
            
            fhTCardCorrEClusterDiffExo[i][tm]  = new TH2F 
            (Form("hTCardCorrEClusterDiffExo_Case%d%s",i,add[tm].Data()),
             Form("#it{E}_{cluster}-#it{E}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, exoticity > 0.97, case %d %s",i,add[tm].Data()),
             nptbins,ptmin,ptmax,210,-1,20); 
            fhTCardCorrEClusterDiffExo[i][tm]->SetXTitle("#it{E} (GeV)");
            fhTCardCorrEClusterDiffExo[i][tm]->SetYTitle("#it{E}_{cluster}-#it{E}_{cell} (GeV)");
            outputContainer->Add(fhTCardCorrEClusterDiffExo[i][tm]); 
            
            //        fhTCardCorrECellMaxRatExo[i][tm]  = new TH2F 
            //        (Form("hTCardCorrECellMaxRatExo_Case%d%s",i,add[tm].Data()),
            //         Form("#it{E}_{cell}/#it{E}_{cell}^{max} vs #it{E}_{cluster}, for (un)correlated cells in TCard, exoticity > 0.97, case %d %s",i,add[tm].Data()),
            //         nptbins,ptmin,ptmax,110,0,1.1); 
            //        fhTCardCorrECellMaxRatExo[i][tm]->SetXTitle("#it{E} (GeV)");
            //        fhTCardCorrECellMaxRatExo[i][tm]->SetYTitle("#it{E}_{cell}/#it{E}^{max}_{cell}");
            //        outputContainer->Add(fhTCardCorrECellMaxRatExo[i][tm]); 
            //        
            //        fhTCardCorrEClusterRatExo[i][tm]  = new TH2F 
            //        (Form("hTCardCorrEClusterRatExo_Case%d%s",i,add[tm].Data()),
            //         Form("#it{E}_{cell}/#it{E}_{cluster} vs #it{E}_{cluster}, for (un)correlated cells in TCard, exoticity > 0.97, case %d %s",i,add[tm].Data()),
            //         nptbins,ptmin,ptmax,110,0,1.1); 
            //        fhTCardCorrEClusterRatExo[i][tm]->SetXTitle("#it{E} (GeV)");
            //        fhTCardCorrEClusterRatExo[i][tm]->SetYTitle("#it{E}_{cell}/#it{E}_{cluster}");
            //        outputContainer->Add(fhTCardCorrEClusterRatExo[i][tm]); 
            
            fhTCardCorrTCellMaxDiffExo[i][tm]  = new TH2F 
            (Form("hTCardCorrTCellMaxDiffExo_Case%d%s",i,add[tm].Data()),
             Form("#it{t}_{cell}^{max}-#it{t}_{cell} vs #it{E}_{cluster}, for (un)correlated cells in TCard, exoticity > 0.97, case %d %s",i,add[tm].Data()),
             nptbins,ptmin,ptmax,1000,-100,100); 
            fhTCardCorrTCellMaxDiffExo[i][tm]->SetXTitle("#it{E} (GeV)");
            fhTCardCorrTCellMaxDiffExo[i][tm]->SetYTitle("#it{t}_{cell}^{max}-#it{t}_{cell} (ns)");
            outputContainer->Add(fhTCardCorrTCellMaxDiffExo[i][tm]); 
          }
        }
      } // neutral or charged
      
      if(fStudyExotic)
      {
        fhEnergyTMEtaResidual1Cell  = new TH2F("hEnergyTMEtaResidual1Cell","#Delta #eta_{cluster-track} vs #it{E}, n cell = 1",
                                               nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
        fhEnergyTMEtaResidual1Cell->SetXTitle("#it{E} (GeV)");
        fhEnergyTMEtaResidual1Cell->SetYTitle("#Delta #eta_{cluster-track}");
        outputContainer->Add(fhEnergyTMEtaResidual1Cell);    
        
        fhEnergyTMPhiResidual1Cell  = new TH2F("hEnergyTMPhiResidual1Cell","#Delta #varphi_{cluster-track} vs #it{E}, n cell = 1",
                                               nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
        fhEnergyTMPhiResidual1Cell->SetXTitle("#it{E} (GeV)");
        fhEnergyTMPhiResidual1Cell->SetYTitle("#Delta #varphi_{cluster-track}");
        outputContainer->Add(fhEnergyTMPhiResidual1Cell);   
        
        fhEnergyTMEtaResidualExotic  = new TH2F("hEnergyTMEtaResidualExotic","#Delta #eta_{cluster-track} vs #it{E}, exo > 0.97",
                                                nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
        fhEnergyTMEtaResidualExotic->SetXTitle("#it{E} (GeV)");
        fhEnergyTMEtaResidualExotic->SetYTitle("#Delta #eta_{cluster-track}");
        outputContainer->Add(fhEnergyTMEtaResidualExotic);    
        
        fhEnergyTMPhiResidualExotic  = new TH2F("hEnergyTMPhiResidualExotic","#Delta #varphi_{cluster-track} vs #it{E}, exo > 0.97",
                                                nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
        fhEnergyTMPhiResidualExotic->SetXTitle("#it{E} (GeV)");
        fhEnergyTMPhiResidualExotic->SetYTitle("#Delta #varphi_{cluster-track}");
        outputContainer->Add(fhEnergyTMPhiResidualExotic);   
        
        fhEnergyTMEtaResidualTCardCorrNoSelection1Cell  = new TH2F("hEnergyTMEtaResidualTCardCorrNoSelection1Cell","#Delta #eta_{cluster-track} vs #it{E}, n cell = 1",
                                                                   nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
        fhEnergyTMEtaResidualTCardCorrNoSelection1Cell->SetXTitle("#it{E} (GeV)");
        fhEnergyTMEtaResidualTCardCorrNoSelection1Cell->SetYTitle("#Delta #eta_{cluster-track}");
        outputContainer->Add(fhEnergyTMEtaResidualTCardCorrNoSelection1Cell);    
        
        fhEnergyTMPhiResidualTCardCorrNoSelection1Cell  = new TH2F("hEnergyTMPhiResidualTCardCorrNoSelection1Cell","#Delta #varphi_{cluster-track} vs #it{E}, n cell = 1",
                                                                   nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
        fhEnergyTMPhiResidualTCardCorrNoSelection1Cell->SetXTitle("#it{E} (GeV)");
        fhEnergyTMPhiResidualTCardCorrNoSelection1Cell->SetYTitle("#Delta #varphi_{cluster-track}");
        outputContainer->Add(fhEnergyTMPhiResidualTCardCorrNoSelection1Cell);   
        
        fhEnergyTMEtaResidualTCardCorrNoSelectionExotic  = new TH2F("hEnergyTMEtaResidualTCardCorrNoSelectionExotic","#Delta #eta_{cluster-track} vs #it{E}, exo > 0.97",
                                                                    nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
        fhEnergyTMEtaResidualTCardCorrNoSelectionExotic->SetXTitle("#it{E} (GeV)");
        fhEnergyTMEtaResidualTCardCorrNoSelectionExotic->SetYTitle("#Delta #eta_{cluster-track}");
        outputContainer->Add(fhEnergyTMEtaResidualTCardCorrNoSelectionExotic);    
        
        fhEnergyTMPhiResidualTCardCorrNoSelectionExotic  = new TH2F("hEnergyTMPhiResidualTCardCorrNoSelectionExotic","#Delta #varphi_{cluster-track} vs #it{E}, exo > 0.97",
                                                                    nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
        fhEnergyTMPhiResidualTCardCorrNoSelectionExotic->SetXTitle("#it{E} (GeV)");
        fhEnergyTMPhiResidualTCardCorrNoSelectionExotic->SetYTitle("#Delta #varphi_{cluster-track}");
        outputContainer->Add(fhEnergyTMPhiResidualTCardCorrNoSelectionExotic);   
        
        for(Int_t i = 0; i < fNEBinCuts; i++)
        {
          fhTMPhiResidualExoticity[i]  = new TH2F 
          (Form("hTMPhiResidual_EBin%d",i),
           Form("#Delta #varphi_{cluster-track} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
           200,-1,1,nresphibins,resphimin,resphimax); 
          fhTMPhiResidualExoticity[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          fhTMPhiResidualExoticity[i]->SetYTitle("#Delta #varphi_{cluster-track}");
          outputContainer->Add(fhTMPhiResidualExoticity[i]);    
          
          fhTMEtaResidualExoticity[i]  = new TH2F 
          (Form("hTMEtaResidual_EBin%d",i),
           Form("#Delta #eta_{cluster-track} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
           200,-1,1,nresetabins,resetamin,resetamax); 
          fhTMEtaResidualExoticity[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          fhTMEtaResidualExoticity[i]->SetYTitle("#Delta #eta_{cluster-track}");
          outputContainer->Add(fhTMEtaResidualExoticity[i]);   
          
          fhTMPhiResidualExoticityLooseCut[i]  = new TH2F 
          (Form("hTMPhiResidual_LooseCut_EBin%d",i),
           Form("#Delta #varphi_{cluster-track} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
           200,-1,1,nresphibins,resphimin,resphimax); 
          fhTMPhiResidualExoticityLooseCut[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          fhTMPhiResidualExoticityLooseCut[i]->SetYTitle("#Delta #varphi_{cluster-track}");
          outputContainer->Add(fhTMPhiResidualExoticityLooseCut[i]);    
          
          fhTMEtaResidualExoticityLooseCut[i]  = new TH2F 
          (Form("hTMEtaResidual_LooseCut_EBin%d",i),
           Form("#Delta #eta_{cluster-track} vs #it{exoticity}, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
           200,-1,1,nresetabins,resetamin,resetamax); 
          fhTMEtaResidualExoticityLooseCut[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          fhTMEtaResidualExoticityLooseCut[i]->SetYTitle("#Delta #eta_{cluster-track}");
          outputContainer->Add(fhTMEtaResidualExoticityLooseCut[i]);   
          
          //      fhTMPhiResidualExoticityAllSameTCard[i]  = new TH2F 
          //      (Form("hTMPhiResidualAllSameTCard_EBin%d",i),
          //       Form("#Delta #varphi_{cluster-track} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
          //       200,-1,1,nresphibins,resphimin,resphimax); 
          //      fhTMPhiResidualExoticityAllSameTCard[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          //      fhTMPhiResidualExoticityAllSameTCard[i]->SetYTitle("#Delta #varphi_{cluster-track}");
          //      outputContainer->Add(fhTMPhiResidualExoticityAllSameTCard[i]);    
          //      
          //      fhTMEtaResidualExoticityAllSameTCard[i]  = new TH2F 
          //      (Form("hTMEtaResidualAllSameTCard_EBin%d",i),
          //       Form("#Delta #eta_{cluster-track} vs #it{exoticity}, all cells same TCard as leading, %2.2f<#it{E}<%2.2f GeV",fEBinCuts[i],fEBinCuts[i+1]),
          //       200,-1,1,nresetabins,resetamin,resetamax); 
          //      fhTMEtaResidualExoticityAllSameTCard[i]->SetXTitle("#it{F}_{+}=1-#it{E}_{+}/#it{E}_{lead cell}");
          //      fhTMEtaResidualExoticityAllSameTCard[i]->SetYTitle("#Delta #eta_{cluster-track}");
          //      outputContainer->Add(fhTMEtaResidualExoticityAllSameTCard[i]);    
        }
      }
    } // TCard correlation studies
    
    if(fFillClusterMaxCellHisto)
    {
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
      
      fhClusterMaxCellECross  = new TH2F ("hClusterMaxCellECross","1 - Energy in cross around max energy cell / max energy cell vs cluster energy, good clusters",
                                          nptbins,ptmin,ptmax, 400,-1,1.); 
      fhClusterMaxCellECross->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellECross->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell max}");
      outputContainer->Add(fhClusterMaxCellECross);    
    }
    
    if(fStudyBadClusters)
    {
      fhNCellsPerClusterNoCut  = new TH2F ("hNCellsPerClusterNoCut","# cells per cluster vs energy, no bad clusters cut",
                                           nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
      fhNCellsPerClusterNoCut->SetXTitle("#it{E} (GeV)");
      fhNCellsPerClusterNoCut->SetYTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerClusterNoCut);
      
      fhNCellsPerClusterWeirdModNoCut  = new TH2F ("hNCellsPerClusterWeirdNoCutMod","# cells per cluster vs energy, E > 100, no bad clusters cut, per SM",
                                                   nceclbins,nceclmin,nceclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
      fhNCellsPerClusterWeirdModNoCut->SetYTitle("SM number");
      fhNCellsPerClusterWeirdModNoCut->SetXTitle("#it{n}_{cells}");
      outputContainer->Add(fhNCellsPerClusterWeirdModNoCut);
      
      fhBadClusterEnergy  = new TH1F ("hBadClusterEnergy","Bad cluster energy", nptbins,ptmin,ptmax); 
      fhBadClusterEnergy->SetXTitle("#it{E}_{cluster} (GeV) ");
      outputContainer->Add(fhBadClusterEnergy);
      
      fhBadClusterEtaPhi  = new TH2F ("hBadClusterEtaPhi","Bad cluster, #eta vs #varphi, #it{E} > 0.5 GeV", 
                                      netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhBadClusterEtaPhi->SetXTitle("#eta ");
      fhBadClusterEtaPhi->SetXTitle("#varphi (rad) ");
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
      
      fhBadClusterTimeEnergy  = new TH2F ("hBadClusterTimeEnergy","energy vs TOF of reconstructed bad clusters",
                                          nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
      fhBadClusterTimeEnergy->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhBadClusterTimeEnergy->SetYTitle("#it{t} (ns)");
      outputContainer->Add(fhBadClusterTimeEnergy);   
      
      if(fFillPi0PairDiffTime)
      {
        fhBadClusterPairDiffTimeE = new TH2F("hBadClusterPairDiffTimeE","cluster pair time difference (bad - good) vs E from bad cluster",nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
        fhBadClusterPairDiffTimeE->SetXTitle("#it{E}_{bad cluster} (GeV)");
        fhBadClusterPairDiffTimeE->SetYTitle("#Delta #it{t} (ns)");
        outputContainer->Add(fhBadClusterPairDiffTimeE);    
      }
      
      if( fFillClusterMaxCellHisto)
      {
        fhClusterMaxCellDiffNoCut  = new TH2F ("hClusterMaxCellDiffNoCut","energy vs difference of cluster energy - max cell energy / cluster energy",
                                               nptbins,ptmin,ptmax, 500,0,1.); 
        fhClusterMaxCellDiffNoCut->SetXTitle("#it{E}_{cluster} (GeV) ");
        fhClusterMaxCellDiffNoCut->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
        outputContainer->Add(fhClusterMaxCellDiffNoCut);  
        
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
          
          //      fhBadClusterMaxCellDiffAverageTime = new TH2F ("hBadClusterMaxCellDiffAverageTime","#it{t}_{cell max}-#it{t}_{average} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
          //      fhBadClusterMaxCellDiffAverageTime->SetXTitle("#it{E} (GeV)");
          //      fhBadClusterMaxCellDiffAverageTime->SetYTitle("#Delta #it{t}_{cell max - average} (ns)");
          //      outputContainer->Add(fhBadClusterMaxCellDiffAverageTime);
          //            
          //      fhBadClusterMaxCellDiffWeightedTime = new TH2F ("hBadClusterMaxCellDiffWeightedTime","#it{t}_{cell max}-#it{t}_{weighted} from bad cluster", nptbins,ptmin,ptmax, tdbins,tdmin,tdmax);
          //      fhBadClusterMaxCellDiffWeightedTime->SetXTitle("#it{E} (GeV)");
          //      fhBadClusterMaxCellDiffWeightedTime->SetYTitle("#Delta #it{t}_{cell max - weighted} (ns)");
          //      outputContainer->Add(fhBadClusterMaxCellDiffWeightedTime);      
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
        if(imod < fFirstModule || imod > fLastModule) continue;

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
      fhTrackMatchedDEtaNeg  = new TH2F("hTrackMatchedDEtaNeg","d#eta of cluster-negative track vs cluster energy",
                                        nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNeg->SetYTitle("d#eta");
      fhTrackMatchedDEtaNeg->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNeg  = new TH2F("hTrackMatchedDPhiNeg","d#varphi of cluster-negative track vs cluster energy",
                                        nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNeg->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDPhiNeg->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiNeg  = new TH2F("hTrackMatchedDEtaDPhiNeg","d#eta vs d#varphi of cluster- negative track vs cluster energy",
                                            nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiNeg->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDEtaDPhiNeg->SetXTitle("d#eta");
      
      fhTrackMatchedDEtaPos  = new TH2F("hTrackMatchedDEtaPos","d#eta of cluster-positive track vs cluster energy",
                                        nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaPos->SetYTitle("d#eta");
      fhTrackMatchedDEtaPos->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiPos  = new TH2F("hTrackMatchedDPhiPos","d#varphi of cluster-positive track vs cluster energy",
                                        nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiPos->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDPhiPos->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiPos  = new TH2F("hTrackMatchedDEtaDPhiPos","d#eta vs d#varphi of cluster-positive track vs cluster energy",
                                            nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiPos->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDEtaDPhiPos->SetXTitle("d#eta");
      
      
      fhTrackMatchedDEtaNegMod  = new TH2F("hTrackMatchedDEtaNegPerModule","d#eta of cluster-negative track vs module, E > 0.5 GeV",
                                           nresetabins,resetamin,resetamax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fhTrackMatchedDEtaNegMod->SetXTitle("d#eta");
      fhTrackMatchedDEtaNegMod->SetYTitle("Module");
      
      fhTrackMatchedDPhiNegMod  = new TH2F("hTrackMatchedDPhiNegPerModule","d#varphi of cluster-negative track vs module, E > 0.5 GeV",
                                           nresetabins,resetamin,resetamax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fhTrackMatchedDPhiNegMod->SetXTitle("d#varphi (rad)");
      fhTrackMatchedDPhiNegMod->SetYTitle("Module");
      
      fhTrackMatchedDEtaPosMod  = new TH2F("hTrackMatchedDEtaPosPerModule","d#eta of cluster-positive track vs module, E > 0.5 GeV",
                                           nresetabins,resetamin,resetamax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fhTrackMatchedDEtaPosMod->SetXTitle("d#eta");
      fhTrackMatchedDEtaPosMod->SetYTitle("Module");
      
      fhTrackMatchedDPhiPosMod  = new TH2F("hTrackMatchedDPhiPosPerModule","d#varphi of cluster-positive track vs module, E > 0.5 GeV",
                                           nresetabins,resetamin,resetamax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fhTrackMatchedDPhiPosMod->SetXTitle("d#varphi (rad)");
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
      
      fhPhiCharged  = new TH1F ("hPhiCharged","#varphi reconstructed clusters, matched with track",nphibins,phimin,phimax);
      fhPhiCharged->SetXTitle("#varphi (rad)");
      outputContainer->Add(fhPhiCharged);
      
      fhEtaCharged  = new TH1F ("hEtaCharged","#eta reconstructed clusters, matched with track",netabins,etamin,etamax);
      fhEtaCharged->SetXTitle("#eta ");
      outputContainer->Add(fhEtaCharged);
      
      if(fFillAllTH3)
      {
        fhEtaPhiECharged  = new TH3F ("hEtaPhiECharged","#eta vs #varphi, reconstructed clusters, matched with track",
                                      netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
        fhEtaPhiECharged->SetXTitle("#eta ");
        fhEtaPhiECharged->SetYTitle("#varphi ");
        fhEtaPhiECharged->SetZTitle("#it{E} (GeV) ");
        outputContainer->Add(fhEtaPhiECharged);	
      }
      else
      {
        fhEtaPhiCharged  = new TH2F ("hEtaPhiCharged","#eta vs #varphi for #it{E} > 0.5 GeV, reconstructed clusters, with matched track",
                                     netabins,etamin,etamax,nphibins,phimin,phimax); 
        fhEtaPhiCharged->SetXTitle("#eta ");
        fhEtaPhiCharged->SetYTitle("#varphi (rad)");
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
      
      fh1EOverPMod = new TH2F
      ("h1EOverP_PerModule","TRACK matches #it{E}/#it{p}, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       nPoverEbins,eOverPmin,eOverPmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh1EOverPMod->SetXTitle("#it{E}/#it{p}");
      fh1EOverPMod->SetYTitle("Module");
      outputContainer->Add(fh1EOverPMod);
      
      fh2dRMod = new TH2F
      ("h2dR_PerModule","TRACK matches #Delta #it{R}, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       ndRbins,dRmin,dRmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh2dRMod->SetXTitle("#Delta #it{R} (rad)");
      fh2dRMod->SetYTitle("Module");
      outputContainer->Add(fh2dRMod) ;
      
      fh2MatchdEdxMod = new TH2F
      ("h2MatchdEdx_PerModule","#it{dE/dx} vs. #it{p} for all matches, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       ndedxbins,dedxmin,dedxmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh2MatchdEdxMod->SetYTitle("Module");
      fh2MatchdEdxMod->SetXTitle("<#it{dE/dx}>");
      outputContainer->Add(fh2MatchdEdxMod);
      
      fh2EledEdxMod = new TH2F
      ("h2EledEdx_PerModule","#it{dE/dx} vs. #it{p} for electrons, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       ndedxbins,dedxmin,dedxmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh2EledEdxMod->SetYTitle("Module");
      fh2EledEdxMod->SetXTitle("<#it{dE/dx}>");
      outputContainer->Add(fh2EledEdxMod) ;
      
      fh1EOverPR02Mod = new TH2F
      ("h1EOverPR02_PerModule","TRACK matches #it{E}/#it{p}, all, #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       nPoverEbins,eOverPmin,eOverPmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh1EOverPR02Mod->SetXTitle("#it{E}/#it{p}");
      fh1EOverPR02Mod->SetYTitle("Module");
      outputContainer->Add(fh1EOverPR02Mod);	
      
      fh1EleEOverPMod = new TH2F
      ("h1EleEOverP_PerModule","Electron candidates #it{E}/#it{p} (60<#it{dE/dx}<100), #it{E}_{cl}&#it{p}_{tr}>0.5 Gev/#it{c}",
       nPoverEbins,eOverPmin,eOverPmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
      fh1EleEOverPMod->SetXTitle("#it{E}/#it{p}");
      fh1EleEOverPMod->SetYTitle("Module");
      outputContainer->Add(fh1EleEOverPMod);
    }
    
    if(fFillAllPi0Histo)
    {
      if ( fFirstModule < 12 )
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
        
        if(fFillInvMassInEMCALWithPHOSDCalAcc)
        {
          fhIMEMCALPHOS  = new TH2F ("hIMEMCALPHOS","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1",
                                     nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
          fhIMEMCALPHOS->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
          fhIMEMCALPHOS->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
          outputContainer->Add(fhIMEMCALPHOS);
          
          fhIMEMCALPHOSSame  = new TH2F ("hIMEMCALPHOSSame","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same #varphi sector",
                                         nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
          fhIMEMCALPHOSSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
          fhIMEMCALPHOSSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
          outputContainer->Add(fhIMEMCALPHOSSame);    
        }
      }
      
      if ( fNModules > 12 && (GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL) && fLastModule > 11 )
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
        
        fhIMDCALPHOSSame  = new TH2F ("hIMDCALPHOSSame","Cluster pairs in DCAL-PHOS Invariant mass vs reconstructed pair #it{p}_{T}, ncell > 1, same #varphi sector",
                                      nptbins,ptmin,ptmax,nmassbins,massmin,massmax); 
        fhIMDCALPHOSSame->SetXTitle("#it{p}_{T, cluster pairs} (GeV) ");
        fhIMDCALPHOSSame->SetYTitle("M_{cluster pairs} (GeV/#it{c}^{2})");
        outputContainer->Add(fhIMDCALPHOSSame);
      }
      
      if(fFillPi0PairDiffTime)
      {
        if ( fFirstModule < 12 )
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
        }
        
        if ( fNModules > 12 && (GetCalorimeter() == kEMCAL || GetCalorimeter() == kDCAL)  && fLastModule > 11 )
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
        fhEtaPhiAmpCell  = new TH3F ("hEtaPhiAmpCell","Cell #eta vs cell #varphi vs cell energy",
                                     netabins,etamin,etamax,nphibins,phimin,phimax,nptbins,ptmin,ptmax); 
        fhEtaPhiAmpCell->SetXTitle("#eta ");
        fhEtaPhiAmpCell->SetYTitle("#varphi (rad)");
        fhEtaPhiAmpCell->SetZTitle("#it{E} (GeV) ");
        outputContainer->Add(fhEtaPhiAmpCell);
      }
      else
      {
        fhEtaPhiCell  = new TH2F ("hEtaPhiCell","Cell #eta vs cell #varphi vs cell energy",
                                  netabins,etamin,etamax,nphibins,phimin,phimax); 
        fhEtaPhiCell->SetXTitle("#eta ");
        fhEtaPhiCell->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiCell);
      }
    }
    
    if(fFillAllCellTimeHisto && fFillClusterMaxCellHisto)
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
    }
    
    // Module histograms
    
    fhEMod  = new TH2F ("hE_Mod","Cluster reconstructed Energy in each present Module",nptbins,ptmin,ptmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhEMod->SetXTitle("#it{E} (GeV)");
    fhEMod->SetYTitle("Module");
    outputContainer->Add(fhEMod);
    
    fhEWeirdMod  = new TH2F ("hEWeird_Mod","Cluster reconstructed Energy in each present Module, ridiculously large E",200,0,10000,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhEWeirdMod->SetXTitle("#it{E} (GeV)");
    fhEWeirdMod->SetYTitle("Module");
    outputContainer->Add(fhEWeirdMod);
    
    fhNClustersMod  = new TH2F ("hNClusters_Mod","# clusters vs Module", nclbins,nclmin+0.5,nclmax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhNClustersMod->SetXTitle("number of clusters");
    fhNClustersMod->SetYTitle("Module");
    outputContainer->Add(fhNClustersMod);
    
    fhSumClustersEnergyMod  = new TH2F ("hSumClustersEnergy_Mod","# clusters vs Module", 1000, 0, 2000,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhSumClustersEnergyMod->SetXTitle("#Sigma_{clusters} #it{E} (GeV)");
    fhSumClustersEnergyMod->SetYTitle("Module");
    outputContainer->Add(fhSumClustersEnergyMod);
    
    fhNClustersSumEnergyPerMod = new TH2F*[fNModules];
    
    fhIMMod                    = new TH2F*[fNModules];
    
    fhNCellsPerClusterMod      = new TH2F*[fNModules];
    
    if(fStudyBadClusters) 
      fhNCellsPerClusterModNoCut = new TH2F*[fNModules];
    
    
    for(Int_t imod = 0; imod < fNModules; imod++)
    {
      if(imod < fFirstModule || imod > fLastModule) continue;
      
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
    
    if(fStudyM02Dependence)
    {
      fhClusterTimeEnergyM02  = new TH3F ("hClusterTimeEnergyM02","energy vs TOF vs M02, reconstructed clusters",
                                          11,5.5,16.5,45,-25,20,100,0.,1); 
      fhClusterTimeEnergyM02->SetXTitle("#it{E} (GeV) ");
      fhClusterTimeEnergyM02->SetYTitle("TOF (ns)");
      fhClusterTimeEnergyM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterTimeEnergyM02);    
      
      if(fFillAllCellTimeHisto)
      {
        fhCellTimeSpreadRespectToCellMaxM02 = new TH3F ("hCellTimeSpreadRespectToCellMaxM02","t_{cell max}-t_{cell i} vs E vs M02", 
                                                        11,5.5,16.5,100,-100,100,100,0.,1); 
        fhCellTimeSpreadRespectToCellMaxM02->SetXTitle("#it{E} (GeV)");
        fhCellTimeSpreadRespectToCellMaxM02->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
        fhCellTimeSpreadRespectToCellMaxM02->SetZTitle("#lambda_{0}^{2}");
        outputContainer->Add(fhCellTimeSpreadRespectToCellMaxM02);
      }
      
      fhClusterMaxCellCloseCellRatioM02  = new TH3F ("hClusterMaxCellCloseCellRatioM02","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                                     11,5.5,16.5, 100,0,1.,100,0.,1); 
      fhClusterMaxCellCloseCellRatioM02->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellCloseCellRatioM02->SetYTitle("#it{E}_{cell i}/#it{E}_{cell max}");
      fhClusterMaxCellCloseCellRatioM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellCloseCellRatioM02);
      
      fhClusterMaxCellCloseCellDiffM02  = new TH3F ("hClusterMaxCellCloseCellDiffM02","energy vs ratio of max cell / neighbour cell, reconstructed clusters",
                                                    11,5.5,16.5, 100,0,100,100,0.,1); 
      fhClusterMaxCellCloseCellDiffM02->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellCloseCellDiffM02->SetYTitle("#it{E}_{cell max}-#it{E}_{cell i} (GeV)");
      fhClusterMaxCellCloseCellDiffM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellCloseCellDiffM02);
      
      fhClusterMaxCellDiffM02  = new TH3F ("hClusterMaxCellDiffM02","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                           11,5.5,16.5, 100,0,1.,100,0.,1); 
      fhClusterMaxCellDiffM02->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellDiffM02->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
      fhClusterMaxCellDiffM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellDiffM02);  
      
      fhClusterMaxCellECrossM02  = new TH3F ("hClusterMaxCellECrossM02",
                                             "1 - Energy in cross around max energy cell / max energy cell vs cluster energy, good clusters",
                                             11,5.5,16.5, 100,0.0,1.,100,0.,1); 
      fhClusterMaxCellECrossM02->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellECrossM02->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell max}");
      fhClusterMaxCellECrossM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellECrossM02);
      
      fhNCellsPerClusterM02  = new TH3F ("hNCellsPerClusterM02","# cells per cluster vs energy",
                                         11,5.5,16.5, 35,0.5,35.5,100,0.,1); 
      fhNCellsPerClusterM02->SetXTitle("#it{E} (GeV)");
      fhNCellsPerClusterM02->SetYTitle("#it{n}_{cells}");
      fhNCellsPerClusterM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsPerClusterM02);    
    }
  }
  
  // Calorimeter cells
  //
  if(fFillAllCellHistograms)
  {
    fhNCells  = new TH1F ("hNCells","# cells", ncebins,ncemin+0.5,ncemax); 
    fhNCells->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCells);
    
    fhNCellsCutAmpMin  = new TH1F ("hNCellsCutAmpMin",Form("# cells amp > %1.2f-%1.2f",fEMCALCellAmpMin,fPHOSCellAmpMin), ncebins,ncemin+0.5,ncemax);
    fhNCellsCutAmpMin->SetXTitle("#it{n}_{cells}");
    outputContainer->Add(fhNCellsCutAmpMin);
    
    fhAmplitude  = new TH1F ("hAmplitude","#it{E}_{cell}", nptbins,ptmin,ptmax/2);
    fhAmplitude->SetXTitle("#it{E}_{cell} (GeV)");
    outputContainer->Add(fhAmplitude);
    
    fhAmpMod  = new TH2F ("hAmp_Mod","Cell energy in each present Module",nptbins,ptmin,ptmax/2,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhAmpMod->SetXTitle("#it{E} (GeV)");
    fhAmpMod->SetYTitle("Module");
    outputContainer->Add(fhAmpMod);
    
    fhAmpWeirdMod  = new TH2F ("hAmpWeird_Mod","Cell energy in each present Module, ridiculously large E",200,0,10000,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhAmpWeirdMod->SetXTitle("#it{E} (GeV)");
    fhAmpWeirdMod->SetYTitle("Module");
    outputContainer->Add(fhAmpWeirdMod);
    
    if(fFillAllCellTimeHisto) 
    {
      fhTimeMod  = new TH2F ("hTime_Mod","Cell time in each present Module",ntimebins,timemin,timemax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
      fhTimeMod->SetXTitle("t (ns)");
      fhTimeMod->SetYTitle("Module");
      outputContainer->Add(fhTimeMod);
    }
    
    if(fFillAllCellAbsIdHistograms)
    {
      //..Create the fhAmpId TH2D with increasing binwidth
      //..0-10 GeV (0.05), 10-20 GeV (0.2), 20-30 GeV (0.5)
      Double_t binWidth=(ptfinemax-ptfinemin)/nfineptbins;
      TCustomBinning xBinning;
      xBinning.SetMinimum(ptfinemin);
      xBinning.AddStep(ptfinemax,binWidth);     //..first entries of the array are the set ranges and bins
      xBinning.AddStep(ptfinemax*2,binWidth*4); //..expand the previously defined range by 2 but increase the bin width
      xBinning.AddStep(ptfinemax*4,binWidth*10);//..expand the previously defined range by 4 but increase the bin width
      
      TCustomBinning yBinning;
      yBinning.SetMinimum(0);
      yBinning.AddStep(fNMaxRows*fNMaxCols*fNModules,1); //..add cells with binwidth 1
      
      TArrayD xbinsArray;
      xBinning.CreateBinEdges(xbinsArray);
      TArrayD ybinsArray;
      yBinning.CreateBinEdges(ybinsArray);
      
      fhAmpId  = new TH2F ("hAmpId","#it{E}_{cell}", xbinsArray.GetSize() - 1, xbinsArray.GetArray(), ybinsArray.GetSize() - 1, ybinsArray.GetArray());
      fhAmpId->SetXTitle("#it{E}_{cell} (GeV)");
      outputContainer->Add(fhAmpId);
      
      fhAmpIdLowGain  = new TH2F ("hAmpIdLG","Low gain: #it{E}_{cell}", nfineptbins,15,ptfinemax+15,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules);
      fhAmpIdLowGain->SetXTitle("#it{E}_{cell} (GeV)");
      outputContainer->Add(fhAmpIdLowGain);
    }
    
    if(fFillAllCellTimeHisto)
    {
      fhTime  = new TH1F ("hTime","#it{t}_{cell}",ntimebins,timemin,timemax); 
      fhTime->SetXTitle("#it{t}_{cell} (ns)");
      outputContainer->Add(fhTime);
      
      //    fhTimeVz  = new TH2F ("hTimeVz","#it{t}_{cell} vs vertex, amplitude > 0.5 GeV",100, 0, 50,ntimebins,timemin,timemax); 
      //    fhTimeVz->SetXTitle("|v_{z}| (cm)");
      //    fhTimeVz->SetYTitle("#it{t}_{cell} (ns)");
      //    outputContainer->Add(fhTimeVz);
      
      fhTimeAmp  = new TH2F ("hTimeAmp","#it{t}_{cell} vs #it{E}_{cell}",nptbins,ptmin,ptmax/2,ntimebins,timemin,timemax); 
      fhTimeAmp->SetYTitle("#it{t}_{cell} (ns)");
      fhTimeAmp->SetXTitle("#it{E}_{cell} (GeV)");
      outputContainer->Add(fhTimeAmp);
      
      fhTimeAmpLowGain  = new TH2F ("hTimeAmpLG","Low gain: #it{t}_{cell} vs #it{E}_{cell}",nptbins,ptmin,ptmax/2,ntimebins,timemin,timemax);
      fhTimeAmpLowGain->SetYTitle("#it{t}_{cell} (ns)");
      fhTimeAmpLowGain->SetXTitle("#it{E}_{cell} (GeV)");
      outputContainer->Add(fhTimeAmpLowGain);
      
      if(fFillAllCellAbsIdHistograms)
      {
        fhCellIdCellLargeTimeSpread= new TH1F ("hCellIdCellLargeTimeSpread","Cells with time 100 ns larger than cell max in cluster ", 
                                               fNMaxCols*fNMaxRows*fNModules,0,fNMaxCols*fNMaxRows*fNModules); 
        fhCellIdCellLargeTimeSpread->SetXTitle("Absolute Cell Id");
        outputContainer->Add(fhCellIdCellLargeTimeSpread);
        
        fhTimeId  = new TH2F ("hTimeId","#it{t}_{cell} vs Absolute Id",
                              ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules); 
        fhTimeId->SetXTitle("#it{t}_{cell} (ns)");
        fhTimeId->SetYTitle("Cell Absolute Id");
        outputContainer->Add(fhTimeId);
        
        fhTimeIdLowGain  = new TH2F ("hTimeIdLG","Low gain: #it{t}_{cell} vs Absolute Id",
                                     ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules);
        fhTimeIdLowGain->SetXTitle("#it{t}_{cell} (ns)");
        fhTimeIdLowGain->SetYTitle("Cell Absolute Id");
        outputContainer->Add(fhTimeIdLowGain);
        
        if(GetCaloUtils()->IsL1PhaseInTimeRecalibrationOn()==1)
        {
          fhTimeL1UnCorrId  = new TH2F ("hTimeL1UnCorrId","#it{t}_{cell} vs Absolute Id",
                                        ntimebins,timemin,timemax,fNMaxRows*fNMaxCols*fNModules,0,fNMaxRows*fNMaxCols*fNModules);
          fhTimeL1UnCorrId->SetXTitle("#it{t}_{cell} (ns)");
          fhTimeL1UnCorrId->SetYTitle("Cell Absolute Id");
          outputContainer->Add(fhTimeL1UnCorrId);
        }
      }
      
      for(Int_t bc = 0; bc < 4; bc++)
      {
        fhTimePerSMPerBC[bc]  = new TH2F (Form("hTimePerSM_BC%d",bc),
                                          Form("#it{t}_{cell} vs super-module, for BC/4=%d",bc),
                                          ntimebins,timemin,timemax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
        fhTimePerSMPerBC[bc]->SetXTitle("#it{t}_{cell} (ns)");
        fhTimePerSMPerBC[bc]->SetYTitle("Module");
        outputContainer->Add(fhTimePerSMPerBC[bc]);
      }
    }
    
    fhCellECross  = new TH2F ("hCellECross","1 - Energy in cross around cell /  cell energy",
                              nptbins,ptmin,ptmax/2, 400,-1,1.); 
    fhCellECross->SetXTitle("#it{E}_{cell} (GeV) ");
    fhCellECross->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell}");
    outputContainer->Add(fhCellECross);    
    
    fhNCellsMod  = new TH2F ("hNCells_Mod","# cells vs Module", ncebins,ncemin+0.5,ncemax,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhNCellsMod->SetXTitle("#it{n}_{cells}");
    fhNCellsMod->SetYTitle("Module");
    outputContainer->Add(fhNCellsMod);
    
    fhSumCellsAmpMod  = new TH2F ("hSumCellsAmp_Mod","# cells vs Module", 1000, 0, 2000,totalSM,fFirstModule-0.5,fLastModule+0.5); 
    fhSumCellsAmpMod->SetXTitle("#Sigma_{cells} #it{Amp} (GeV)");
    fhSumCellsAmpMod->SetYTitle("Module");
    outputContainer->Add(fhSumCellsAmpMod);
    
    fhGridCells  = new TH2F ("hGridCells",Form("Entries in grid of cells"), 
                             ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax); 
    fhGridCells->SetYTitle("row (phi direction)");
    fhGridCells->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCells);
    
    fhGridCellsE  = new TH2F ("hGridCellsE","Accumulated energy in grid of cells", 
                              ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax); 
    fhGridCellsE->SetYTitle("row (phi direction)");
    fhGridCellsE->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsE);
    
    fhGridCellsLowGain  = new TH2F ("hGridCellsLG",Form("Low gain: Entries in grid of cells"),
                                    ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
    fhGridCellsLowGain->SetYTitle("row (phi direction)");
    fhGridCellsLowGain->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsLowGain);
    
    fhGridCellsELowGain  = new TH2F ("hGridCellsELG","Low gain: Accumulated energy in grid of cells",
                                     ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
    fhGridCellsELowGain->SetYTitle("row (phi direction)");
    fhGridCellsELowGain->SetXTitle("column (eta direction)");
    outputContainer->Add(fhGridCellsELowGain);
    
    if(fFillAllCellTimeHisto)
    { 
      fhGridCellsTime  = new TH2F ("hGridCellsTime","Accumulated time in grid of cells",
                                   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhGridCellsTime->SetYTitle("row (phi direction)");
      fhGridCellsTime->SetXTitle("column (eta direction)");
      outputContainer->Add(fhGridCellsTime);
      
      fhGridCellsTimeLowGain  = new TH2F ("hGridCellsTimeLG","Low gain: Accumulated time in grid of cells",
                                          ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax);
      fhGridCellsTimeLowGain->SetYTitle("row (phi direction)");
      fhGridCellsTimeLowGain->SetXTitle("column (eta direction)");
      outputContainer->Add(fhGridCellsTimeLowGain);
    }
    
    if(fFillAllCellAbsIdHistograms)
    {
      fhNCellsSumAmpPerMod       = new TH2F*[fNModules];  
      
      for(Int_t imod = 0; imod < fNModules; imod++)
      {
        if(imod < fFirstModule || imod > fLastModule) continue;
        
        fhNCellsSumAmpPerMod[imod] = new TH2F (Form("hNCellsSumAmp_Mod%d",imod),
                                               Form("# cells in SM vs sum of cells energy in Module %d",imod), 
                                               nptbins,ptmin,ptmax*4, ncebins,ncemin,ncemax); 
        fhNCellsSumAmpPerMod[imod]->SetXTitle("#Sigma #it{Amplitude} (GeV)");
        fhNCellsSumAmpPerMod[imod]->SetYTitle("#Sigma #it{n}_{cells}");
        outputContainer->Add(fhNCellsSumAmpPerMod[imod]);}
    } 
    
    if(fFillAllCellTimeHisto && fFillAllCellAbsIdHistograms) 
    {
      fhTimeAmpPerRCU = new TH2F*[fNModules*fNRCU];  
      
      for(Int_t imod = 0; imod < fNModules; imod++)
      {
        for(Int_t ircu = 0; ircu < fNRCU; ircu++)
        {
          if( ircu ==1 && 
             (imod == 10 || imod== 11 || imod == 18 || imod == 19) 
             ) continue;
          
          fhTimeAmpPerRCU[imod*fNRCU+ircu]  = new TH2F (Form("hTimeAmp_Mod%d_RCU%d",imod,ircu),
                                                        Form("#it{E}_{cell} vs #it{t}_{cell} in Module %d, RCU %d ",imod,ircu), 
                                                        nptbins,ptmin,ptmax/2,ntimebins,timemin,timemax); 
          fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetXTitle("#it{E} (GeV)");
          fhTimeAmpPerRCU[imod*fNRCU+ircu]->SetYTitle("#it{t} (ns)");
          outputContainer->Add(fhTimeAmpPerRCU[imod*fNRCU+ircu]);
          
        }
      }
    }
  }
  
  // Detectors correlation
  //
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
  
  // Monte Carlo Histograms
  //
  if(IsDataMC())
  {    
    TString particleName[] = {
      "Photon",        "Pi0",         "Eta",
      "Electron",      "PhotonConv",
      "NeutralHadron", "ChargedHadron"      };
    
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
                                      Form("Y vs #varphi of generated %s",particleName[iPart].Data()),
                                      200,-1,1,360,0,TMath::TwoPi());
      
      fhGenMCE [iPart]    ->SetXTitle("#it{E} (GeV)");
      fhGenMCPt[iPart]    ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhGenMCEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCEtaPhi[iPart]->SetYTitle("#varphi (rad)");
      
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
                                         Form("Y vs #varphi of generated %s",particleName[iPart].Data()),
                                         netabins,etamin,etamax,nphibins,phimin,phimax);
      
      fhGenMCAccE [iPart]    ->SetXTitle("#it{E} (GeV)");
      fhGenMCAccPt[iPart]    ->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhGenMCAccEtaPhi[iPart]->SetXTitle("#eta");
      fhGenMCAccEtaPhi[iPart]->SetYTitle("#varphi (rad)");
      
      outputContainer->Add(fhGenMCAccE     [iPart]);
      outputContainer->Add(fhGenMCAccPt    [iPart]);
      outputContainer->Add(fhGenMCAccEtaPhi[iPart]);
      
    }    
    
    if(fFillAllClusterHistograms)
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
                                                    Form("Generated - Reconstructed #varphi, %s, Matched %d",particleName[iPart].Data(),iCh),
                                                    nptbins, ptmin, ptmax, nphibins*2,-phimax,phimax); 
          fhRecoMCDeltaPhi[iPart][iCh]->SetYTitle("#Delta #varphi (rad)");
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
                                               Form("#varphi distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh),
                                               nphibins,phimin,phimax, nphibins,phimin,phimax); 
          fhRecoMCPhi[iPart][iCh]->SetXTitle("#varphi_{reconstructed} (rad)");
          fhRecoMCPhi[iPart][iCh]->SetYTitle("#varphi_{generated} (rad)");
          outputContainer->Add(fhRecoMCPhi[iPart][iCh]);
          
          fhRecoMCEta[iPart][iCh]  = new TH2F (Form("hRecoMCEta_%s_Match%d",particleName[iPart].Data(),iCh),
                                               Form("#eta distribution, reconstructed vs generated, %s, Matched %d",particleName[iPart].Data(),iCh), 
                                               netabins,etamin,etamax,netabins,etamin,etamax); 
          fhRecoMCEta[iPart][iCh]->SetXTitle("#eta_{reconstructed} ");
          fhRecoMCEta[iPart][iCh]->SetYTitle("#eta_{generated} ");
          outputContainer->Add(fhRecoMCEta[iPart][iCh]);
        }
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
    
  fTimeCutMin      = -9999999;
  fTimeCutMax      =  9999999;
  
  fEMCALCellAmpMin = 0.2; // 200 MeV
  fPHOSCellAmpMin  = 0.2; // 200 MeV
  fCellAmpMin      = 0.2; // 200 MeV

  fEMCALClusterM02Min = 0.05; 
  
  fEMCALClusterNCellMin = 2; // at least 2
  fPHOSClusterNCellMin  = 3; // at least 3
  
  fInvMassMinECut  = 0.5; // 500 MeV
  fInvMassMaxECut  = 10; 

  fInvMassMinM02Cut = 0.1; 
  fInvMassMaxM02Cut = 0.4; 
  
  fInvMassMaxTimeDifference = 200; // ns, quite open 
  
  fInvMassMaxOpenAngle = 100*TMath::DegToRad(); // 100 degrees
  
//  // Exotic studies
//  fExoNECrossCuts  = 10 ;
//  fExoNDTimeCuts   = 4  ;
//  
//  fExoDTimeCuts [0] = 1.e4 ; fExoDTimeCuts [1] = 50.0 ; fExoDTimeCuts [2] = 25.0 ; fExoDTimeCuts [3] = 10.0 ;
//  fExoECrossCuts[0] = 0.80 ; fExoECrossCuts[1] = 0.85 ; fExoECrossCuts[2] = 0.90 ; fExoECrossCuts[3] = 0.92 ; fExoECrossCuts[4] = 0.94 ;
//  fExoECrossCuts[5] = 0.95 ; fExoECrossCuts[6] = 0.96 ; fExoECrossCuts[7] = 0.97 ; fExoECrossCuts[8] = 0.98 ; fExoECrossCuts[9] = 0.99 ;

//  fNEBinCuts = 14;
//  fEBinCuts[0] = 0.;  fEBinCuts[1] = 0.3;  fEBinCuts[2] = 0.5;
//  fEBinCuts[3] = 1.;  fEBinCuts[4] = 2. ;  fEBinCuts[5] = 3. ;
//  fEBinCuts[6] = 4.;  fEBinCuts[7] = 6. ;  fEBinCuts[8] = 8. ;
//  fEBinCuts[9] = 10.; fEBinCuts[10]= 12.;  fEBinCuts[11]= 16.;
//  fEBinCuts[12]= 20.; fEBinCuts[13]= 50.;  fEBinCuts[14]= 100.;
//  for(Int_t i = fNEBinCuts; i < 15; i++) fEBinCuts[i] = 1000.;

  fNEBinCuts = 7;
  fEBinCuts[0] = 2. ; fEBinCuts[1] = 4. ;  fEBinCuts[2] = 6. ;
  fEBinCuts[3] = 8. ; fEBinCuts[4] = 10.;  fEBinCuts[5] = 12.;
  fEBinCuts[6] = 16.; fEBinCuts[7] = 20.;  
  for(Int_t i = fNEBinCuts+1; i < 15; i++) fEBinCuts[i] = 1000.;
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
 
    if(  m02 < fEMCALClusterM02Min || nCellsPerCluster < fEMCALClusterNCellMin ) return kFALSE ; // mild shower shape cut for exotics

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
  printf("EMCAL Min M02   : %2.2f\n"          , fEMCALClusterM02Min) ;
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
  if ( !GetMC() ) return;
 
  Int_t    pdg    =  0 ;
  Int_t    status =  0 ;
  Int_t    nprim  = GetMC()->GetNumberOfTracks();
  
  TParticle        * primStack = 0;
  AliAODMCParticle * primAOD   = 0;
  
  //printf("N primaries %d\n",nprim);
  for(Int_t i=0 ; i < nprim; i++)
  {
    if ( !GetReader()->AcceptParticleMCLabel( i ) ) continue ;
    
    // Get the generated particles, check that it is primary (not parton, resonance)
    // and get its momentum. Different way to recover from ESD or AOD
    if(GetReader()->ReadStack())
    {
      primStack = GetMC()->Particle(i) ;
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
      primAOD = (AliAODMCParticle *) GetMC()->GetTrack(i);
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



