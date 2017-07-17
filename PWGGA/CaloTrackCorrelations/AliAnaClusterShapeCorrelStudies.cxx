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
#include <TH3F.h>
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaClusterShapeCorrelStudies.h"
#include "AliCaloTrackReader.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVEvent.h"

// --- Detectors --- 
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaClusterShapeCorrelStudies) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters.
/// Init all the histogram arrays to 0.
//__________________________________________
AliAnaClusterShapeCorrelStudies::AliAnaClusterShapeCorrelStudies() :
AliAnaCaloTrackCorrBaseClass(),  

// Switches
fStudyShape(kFALSE),                   fStudyShapeParam(kFALSE),
fStudyWeight(kFALSE),               
fStudyTCardCorrelation(kFALSE),        fStudyExotic(kFALSE),

// Parameters and cuts
fM02Min(0),                            fNCellMin(0),            
fMinDistToBad(0),                      fNEBinCuts(0),
fEMinShape(0),                         fEMaxShape(100),
fNCellMinShape(-1),

fdEdXMinEle(0),                        fdEdXMaxEle(0),
fdEdXMinHad(0),                        fdEdXMaxHad(0),

// Invariant mass
fInvMassMinECut(0),                    fInvMassMaxECut(0),                    
fInvMassMinM02Cut(0),                  fInvMassMaxM02Cut(0),                    
fInvMassMaxOpenAngle(0),               fInvMassMaxTimeDifference(0),

fConstantTimeShift(0),

fClusterMomentum(),                    fClusterMomentum2(),                    
fCaloCellList(NULL),                   fCaloClusList(NULL),

// Histograms

// TCard correl and shape and exoticity
fhEnergyTMEtaResidual1Cell(0),         fhEnergyTMPhiResidual1Cell(0),
fhColRowExoticHighE1CellPosTime(0),    fhColRowExoticHighE1CellNegTime(0),     fhColRowExoticHighE1CellNulTime(0),
fhEnergyTMEtaResidualExotic(0),        fhEnergyTMPhiResidualExotic(0),
fhColRowExoticHighEPosTime(0),         fhColRowExoticHighENegTime(0),          fhColRowExoticHighENulTime(0),
fhColRowHighEPosTime(0),               fhColRowHighENegTime(0),                fhColRowHighENulTime(0),
fhEnergyTMEtaResidualTCardCorrNoSelection1Cell(0),  fhEnergyTMPhiResidualTCardCorrNoSelection1Cell(0),
fhEnergyTMEtaResidualTCardCorrNoSelectionExotic(0), fhEnergyTMPhiResidualTCardCorrNoSelectionExotic(0),

// Shape studies
//fhCellTimeSpreadRespectToCellMaxM02(0), 
//fhClusterMaxCellCloseCellDiffM02(0),  
fhClusterMaxCellCloseCellRatioM02(0),  fhClusterMaxCellECrossM02(0),
fhInvMassNCellSM(0),                   fhInvMassNCellSMSame(0),
fhColRowM02(0),                        fhColRowM02NCellCut(0),
fhEMaxCellTimeM02SM(0),                fhEMaxCellTimeNCellSM(0),               fhESecCellTimeNCellSM(0),

// Weight studies
fhECellClusterRatio(0),                fhECellClusterLogRatio(0),                 
fhEMaxCellClusterRatio(0),             fhEMaxCellClusterLogRatio(0),                
fhECellTotalRatio(0),                  fhECellTotalLogRatio(0),
fhECellTotalRatioMod(0),               fhECellTotalLogRatioMod(0)
{
  for(Int_t i=0; i < 3; i++)
  {
    // Shower shape dependence
    fhClusterTimeEnergyM02 [i] = 0;  
    fhClusterMaxCellDiffM02[i] = 0;
    fhNCellsPerClusterM02  [i] = 0;
    fhNCellsPerClusterM20  [i] = 0;
    
    fhNCellsPerClusterMEta    [i] = 0;
    fhNCellsPerClusterMPhi    [i] = 0;
    fhNCellsPerClusterMEtaPhi [i] = 0;
    fhNCellsPerClusterMEtaPhiA[i] = 0;
    
    fhSMNCell              [i] = 0;
    fhSMNCellM02           [i] = 0;
    fhSMM02                [i] = 0;
    fhSMM02NoCut           [i] = 0;
    fhColM02               [i] = 0;
    fhRowM02               [i] = 0;

    fhOriginE              [i] = 0;
    fhOriginM02            [i] = 0;
    
    // Cluster asymmetry
    fhDeltaIEtaDeltaIPhi   [i] = 0;
    fhDeltaIA              [i] = 0;
    fhDeltaIAM02           [i] = 0;         
    fhDeltaIAM20           [i] = 0;
    fhDeltaIANCells        [i] = 0;
    fhDeltaIAOrigin        [i] = 0;
    
    fhDeltaIEtaDeltaIPhiTot[i] = 0;
    fhDeltaIATot           [i] = 0;
    fhDeltaIATotM02        [i] = 0;  
    fhDeltaIATotM20        [i] = 0;
    fhDeltaIATotNCells     [i] = 0;
    fhDeltaIATotOrigin     [i] = 0;    
  }
  
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
  
  for(Int_t i = 0; i < 20; i++)
  {
    fhNCellsPerClusterM02M20PerSM  [i] = 0;
    fhESecCellEMaxCellM02NCellPerSM[i] = 0;           
    fhESecCellEClusterM02NCellPerSM[i] = 0;
    fhESecCellLogM02NCellPerSM     [i] = 0;   
    
    fhEMaxCellEClusterM02NCellPerSM[i] = 0;            
    fhEMaxCellLogM02NCellPerSM     [i] = 0;

    fhEMaxESecCellNCellLowM02PerSM [i] = 0;   
    fhEMaxECrossNCellLowM02PerSM   [i] = 0;  
    fhEMaxESecCellNCellHighM02PerSM[i] = 0;   
    fhEMaxECrossNCellHighM02PerSM  [i] = 0;  
    
    fhColRowFromCellMaxLowM02PerSM[i][0] = 0;
    fhColRowFromCellMaxLowM02PerSM[i][1] = 0;
    
    fhColRowFromCellMaxHighM02PerSM[i][0] = 0;
    fhColRowFromCellMaxHighM02PerSM[i][1] = 0;
    
    for(Int_t j = 0; j < 3; j++)
    {
      fhColRowFromCellMaxEMaxSecDiffLowM02PerSM [i][0][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffLowM02PerSM [i][1][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][0][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][1][j] = 0;
 
      fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM [i][0][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM [i][1][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][0][j] = 0;
      fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][1][j] = 0;
      
      fhColRowFromCellMaxECellClusterRatLowM02PerSM [i][0][j] = 0;
      fhColRowFromCellMaxECellClusterRatLowM02PerSM [i][1][j] = 0;
      fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][0][j] = 0;
      fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][1][j] = 0;
    }
  }

  fhColRowFromCellMaxLowM02[0] = 0;
  fhColRowFromCellMaxLowM02[1] = 0;
  
  fhColRowFromCellMaxHighM02[0] = 0;
  fhColRowFromCellMaxHighM02[1] = 0;
    
  InitParameters();
}


//___________________________________________________
/// Check effect of T-Card channels correlation on shower shape/energy
/// for EMCal
///
/// \param clus: cluster pointer
/// \param matched: bool with loose track-matching info
/// \param absIdMax: id of highest energy cell in cluster
/// \param exoticity: cross energy fraction
///
//___________________________________________________
void AliAnaClusterShapeCorrelStudies::ChannelCorrelationInTCard
(AliVCluster* clus, Bool_t matched,Int_t absIdMax, Float_t exoticity) 
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
  if(time>400) time-=fConstantTimeShift;
  
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
        
        Float_t  eCell = fCaloCellList->GetCellAmplitude(absId);        
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
    
  // in center of SM
  Int_t etaRegion = -1, phiRegion = -1;
  GetCaloUtils()->GetEMCALSubregion(clus,fCaloCellList,etaRegion,phiRegion);
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
  Int_t nlm  = GetCaloUtils()->GetNumberOfLocalMaxima(clus, fCaloCellList, absIdList, maxEList) ; 
//Int_t nlm  = GetCaloUtils()->GetNumberOfLocalMaxima(clus,fCaloCellList);

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
  Float_t  eCellMax = fCaloCellList->GetCellAmplitude(absIdMax);  
  Double_t tCellMax = fCaloCellList->GetCellTime(absIdMax);      
  //printf("Org E %2.2f, t %2.2f\n",eCellMax,tCellMax*1e9);
  GetCaloUtils()->RecalibrateCellAmplitude(eCellMax, GetCalorimeter(), absIdMax);
  GetCaloUtils()->RecalibrateCellTime(tCellMax, GetCalorimeter(), absIdMax, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
  //printf("New E %2.2f, t %2.2f\n",eCellMax,tCellMax*1e9);
  
  tCellMax *= 1.0e9;
  if(tCellMax>400) tCellMax-=fConstantTimeShift;
  
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
  
  //printf("Cluster E %2.2f, ncells %d, absIdMax %d, eCell Max %2.2f\n", energy, ncells, absIdMax, fCaloCellList->GetCellAmplitude(absIdMax));
  
  //
  // Loop on the cluster cells, define correlations
  //
  for (Int_t ipos = 0; ipos < ncells; ipos++) 
  {
    Int_t absId  = clus->GetCellsAbsId()[ipos];   
    
    Float_t  eCell = fCaloCellList->GetCellAmplitude(absId);
    Double_t tCell = fCaloCellList->GetCellTime(absId);      

    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
    GetCaloUtils()->RecalibrateCellTime(tCell, GetCalorimeter(), absId, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    tCell *= 1.0e9;
    if(tCell>400) tCell-=fConstantTimeShift;
    
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
//      eCell = fCaloCellList->GetCellAmplitude(absId2);      
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
    printf("AliAnaClusterShapeCorrelStudies: M02 %f, M20 %f, E %2.3f, ncell %d, n with weight %d; max cell E %2.3f\n",
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

  fhTimeTCardCorrNoSelection   [matched]->Fill(energy,time, GetEventWeight());
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
    
    Float_t  eCell = fCaloCellList->GetCellAmplitude(absId);
    Double_t tCell = fCaloCellList->GetCellTime(absId);      
    
    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
    GetCaloUtils()->RecalibrateCellTime(tCell, GetCalorimeter(), absId, GetReader()->GetInputEvent()->GetBunchCrossNumber());    
    tCell *= 1.0e9;
    if(tCell>400) tCell-=fConstantTimeShift;
    
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
    for(Int_t jclus = 0 ; jclus < fCaloClusList->GetEntriesFast() ; jclus++) 
    {
      AliVCluster* clus2 =  (AliVCluster*) fCaloClusList->At(jclus);
      
      Float_t maxCellFraction = 0.;
      Int_t absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(fCaloCellList, clus2, maxCellFraction);
      
      Double_t tof2 =  clus2->GetTOF()*1.e9;
      if(tof2>400) tof2-=fConstantTimeShift;
      
      Double_t diffTof = tCellMax-tof2;
      
      // Try to reduce background with a mild shower shape cut and no more 
      // than 1 local maximum in cluster and remove low energy clusters
      
      if(   absIdMax == absIdMax2   
         || !IsGoodCluster(absIdMax2, clus2->GetM02(), clus2->GetNCells()) 
         || GetCaloUtils()->GetNumberOfLocalMaxima(clus2,fCaloCellList) > 1 
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

//__________________________________________________________________________________________________________________
/// Fill mostly TH3 histograms with cluster energy vs shower shape or cluster asymmetry vs another parameter.
///
/// \param clus: calorimeter cluster pointer.
/// \param absIdMax: id number of cell with highest energy in the cluster.
/// \param maxFrac: ratio E_cell_max/ E_cluster.
/// \param eCrossFrac: exoticity fraction.
/// \param tmax: time of highest energy cell in cluster.
/// \param matchedPID: 0-neutral, 1-matched to electron track, 2-matched to hadronic track.
/// \param mcIndex: index of origin tagging as photon, pion, etc.
//__________________________________________________________________________________________________________________
void AliAnaClusterShapeCorrelStudies::ClusterShapeHistograms
(AliVCluster* clus , Int_t   absIdMax, Double_t maxFrac  , 
 Float_t eCrossFrac, Float_t eCellMax, Double_t tmax     ,
 Int_t   matchedPID, Int_t    mcIndex)
{
  // By definition a cluster has at least 1 cell, 
  // and shape only makes sense with at least 2
  // in case fNCellMin was open, check again the size
  Int_t nCaloCellsPerCluster = clus->GetNCells();
  if ( nCaloCellsPerCluster < 2 ) return; 
  
  Float_t energy = clus->E();
  Float_t m02    = clus->GetM02();
  Float_t m20    = clus->GetM20();
  Int_t   nCell  = 0;

  Int_t   dIeta    = 0;
  Int_t   dIphi    = 0;
  Int_t   dIetaNeg = 0;
  Int_t   dIphiNeg = 0;
  Int_t   dIetaPos = 0;
  Int_t   dIphiPos = 0;
    
  // Loop on cells in cluster to get number of cells with significant energy
  for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
  {
    Int_t   absId = clus ->GetCellsAbsId()[ipos];
    Float_t eCell = fCaloCellList->GetCellAmplitude(absId) ;
    
    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);

    if( absId == absIdMax || eCell < 0.01 ) continue;
    
    Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
    
    if( weight < 0.01 ) continue;
    
    nCell++;
  }

  if ( nCell < 1 ) return; 

  //
  // Cluster location
  //
  Int_t ietaMax=-1, iphiMax = 0, rcuMax = 0, icolAbs = -1, irowAbs = -1;  
  Int_t smMax = GetModuleNumberCellIndexesAbsCaloMap(absIdMax, GetCalorimeter(), 
                                                     ietaMax, iphiMax, rcuMax, icolAbs, irowAbs);

  if ( matchedPID == 0 && energy > fEMinShape && energy < fEMaxShape )
  {
    fhColRowM02->Fill(icolAbs,irowAbs,m02,GetEventWeight()) ;
    
    if ( nCell > fNCellMinShape )
      fhColRowM02NCellCut->Fill(icolAbs,irowAbs,m02,GetEventWeight()) ;
  }
  //
  
  //
  // Clean the sample with very strict cut on acceptance, select only
  // in center of SM
  //
  Int_t etaRegion = -1, phiRegion = -1;
  GetCaloUtils()->GetEMCALSubregion(clus,fCaloCellList,etaRegion,phiRegion);
  // Region 0: center of SM ~0.18<|eta|<0.55
  if ( etaRegion !=0 ) return ;

  // Loop on cells in cluster to get cell cluster asymmetry and 
  // other correlation parameters
  for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
  {  
    Int_t   absId = clus ->GetCellsAbsId()[ipos];
    Float_t eCell = fCaloCellList->GetCellAmplitude(absId) ;
    
    GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
    
    if( absId == absIdMax || eCell < 0.01 ) continue;
    
    Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
    
    if( weight < 0.01 ) continue;

    Float_t fracCell = eCell/eCellMax;            
    fhClusterMaxCellCloseCellRatioM02->Fill(energy, fracCell, m02, GetEventWeight());
    
    Float_t fracClus = (energy-eCell)/energy;
    
    Float_t logECell = TMath::Log(eCell);
    
    //Float_t ampDiff = fCaloCellList->GetCellAmplitude(absIdMax)-fCaloCellList->GetCellAmplitude(absId);
    //fhClusterMaxCellCloseCellDiffM02 ->Fill(energy,ampDiff, m02,GetEventWeight());
    
    //Double_t time  = fCaloCellList->GetCellTime(absId);
    //GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
    //
    //Float_t tdiff = (tmax-(time*1.0e9-fConstantTimeShift));
    //fhCellTimeSpreadRespectToCellMaxM02->Fill(energy, tdiff, m02, GetEventWeight());   
    
    if ( energy > fEMinShape && energy < fEMaxShape  && matchedPID == 0 ) 
    {
      fhESecCellEMaxCellM02NCellPerSM[smMax]->Fill(fracCell, nCell, m02, GetEventWeight());
      fhESecCellEClusterM02NCellPerSM[smMax]->Fill(fracClus, nCell, m02, GetEventWeight());
      fhESecCellLogM02NCellPerSM     [smMax]->Fill(logECell, nCell, m02, GetEventWeight());
    
      if ( m02 > 0.1 && m02 < 0.3 )
      {
        fhEMaxESecCellNCellLowM02PerSM [smMax]->Fill(eCellMax, eCell     , nCell, GetEventWeight());
        fhEMaxECrossNCellLowM02PerSM   [smMax]->Fill(eCellMax, eCrossFrac, nCell, GetEventWeight());
        
        Double_t time  = fCaloCellList->GetCellTime(absId);
        GetCaloUtils()->RecalibrateCellTime(time, GetCalorimeter(), absId,GetReader()->GetInputEvent()->GetBunchCrossNumber());
        time*=1.e9;
        if(time>400) time-=fConstantTimeShift;
        
        fhESecCellTimeNCellSM->Fill(time, smMax, nCell, GetEventWeight());
      }
      else if  ( m02 > 0.5 && m02 < 2 )
      {
        fhEMaxESecCellNCellHighM02PerSM[smMax]->Fill(eCellMax, eCell     , nCell, GetEventWeight());
        fhEMaxECrossNCellHighM02PerSM  [smMax]->Fill(eCellMax, eCrossFrac, nCell, GetEventWeight());
      }
    }
    
    //if(fStudyShapeParam)
    {
      //////
      // Cluster asymmetry in cell units
      //////
      Int_t ieta=-1; Int_t iphi = 0; Int_t rcu = 0;
      Int_t sm = GetModuleNumberCellIndexes(absId,GetCalorimeter(), ieta, iphi, rcu);
      
      if(dIphi < TMath::Abs(iphi-iphiMax)) dIphi = TMath::Abs(iphi-iphiMax);
      if(iphi-iphiMax < 0 && dIphiNeg > iphi-iphiMax) dIphiNeg = iphi-iphiMax;
      if(iphi-iphiMax > 0 && dIphiPos < iphi-iphiMax) dIphiPos = iphi-iphiMax;
      
      if(smMax==sm)
      {
        if(dIeta < TMath::Abs(ieta-ietaMax)) dIeta = TMath::Abs(ieta-ietaMax);
        if(ieta-ietaMax < 0 && dIetaNeg > ieta-ietaMax) dIetaNeg = ieta-ietaMax;
        if(ieta-ietaMax > 0 && dIetaPos < ieta-ietaMax) dIetaPos = ieta-ietaMax;
        
        if ( energy > fEMinShape && energy < fEMaxShape  && matchedPID == 0 ) 
        {
          Int_t nCellBin = 2;
          if      ( nCell == 2 || nCell == 3 ) nCellBin = 0;
          else if ( nCell == 4 || nCell == 5 ) nCellBin = 1;
          if ( m02 > 0.1 && m02 < 0.3 )
          {
            fhColRowFromCellMaxLowM02PerSM[smMax][ietaMax%2]->Fill(ietaMax-ieta, iphiMax-iphi, nCell, GetEventWeight());
            fhColRowFromCellMaxLowM02            [ietaMax%2]->Fill(ietaMax-ieta, iphiMax-iphi,        GetEventWeight());
            
            fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi,(eCellMax-eCell)/eCellMax, GetEventWeight());
            fhColRowFromCellMaxEMaxSecDiffLowM02PerSM    [smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi, eCellMax-eCell, GetEventWeight());
            fhColRowFromCellMaxECellClusterRatLowM02PerSM[smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi, eCell/energy  , GetEventWeight());
          }
          else if (m02 > 0.5 && m02 < 2)
          {
            fhColRowFromCellMaxHighM02PerSM[smMax][ietaMax%2]->Fill(ietaMax-ieta, iphiMax-iphi, nCell, GetEventWeight());
            fhColRowFromCellMaxHighM02            [ietaMax%2]->Fill(ietaMax-ieta, iphiMax-iphi,        GetEventWeight());
            
            fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi,(eCellMax-eCell)/eCellMax, GetEventWeight());
            fhColRowFromCellMaxEMaxSecDiffHighM02PerSM    [smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi, eCellMax-eCell, GetEventWeight());
            fhColRowFromCellMaxECellClusterRatHighM02PerSM[smMax][ietaMax%2][nCellBin]->Fill(ietaMax-ieta, iphiMax-iphi, eCell/energy  , GetEventWeight());
          }
        }        

      }
      else
      {
        Int_t ietaShift    = ieta;
        Int_t ietaMaxShift = ietaMax;
        
        if (ieta > ietaMax)  ietaMaxShift+=48;
        else                 ietaShift   +=48;
        
        if(dIeta < TMath::Abs(ietaShift-ietaMaxShift)) dIeta = TMath::Abs(ietaShift-ietaMaxShift);
        if(ietaShift-ietaMaxShift < 0 && dIetaNeg > ietaShift-ietaMaxShift) dIetaNeg = ietaShift-ietaMaxShift;
        if(ietaShift-ietaMaxShift > 0 && dIetaPos < ietaShift-ietaMaxShift) dIetaPos = ietaShift-ietaMaxShift;
      }
      
    }
  } // Fill cell-cluster histogram loop
  
  // Col-Row histogram, fill emax/ecluster ratio for highest energy cell.
  if ( energy > fEMinShape && energy < fEMaxShape  && matchedPID == 0 ) 
  {
    fhNCellsPerClusterM02M20PerSM[smMax]->Fill(m20, nCell, m02, GetEventWeight());

    Int_t nCellBin = 2;
    if      ( nCell == 2 || nCell == 3 ) nCellBin = 0;
    else if ( nCell == 4 || nCell == 5 ) nCellBin = 1;
    if ( m02 > 0.1 && m02 < 0.3 )
      fhColRowFromCellMaxECellClusterRatLowM02PerSM [smMax][ietaMax%2][nCellBin]->Fill(0., 0., eCellMax/energy, GetEventWeight());
    else if (m02 > 0.5 && m02 < 2)
      fhColRowFromCellMaxECellClusterRatHighM02PerSM[smMax][ietaMax%2][nCellBin]->Fill(0., 0., eCellMax/energy, GetEventWeight());
  }        
  
  if(fStudyExotic)
    fhClusterMaxCellECrossM02->Fill(energy, eCrossFrac, m02, GetEventWeight());
  
  //
  // Fill histograms only for PID
  //  
  fhClusterMaxCellDiffM02[matchedPID]->Fill(energy, maxFrac, m02, GetEventWeight());
  fhClusterTimeEnergyM02 [matchedPID]->Fill(energy, tmax   , m02, GetEventWeight());
  fhNCellsPerClusterM02  [matchedPID]->Fill(energy, nCell  , m02, GetEventWeight());
  fhNCellsPerClusterM20  [matchedPID]->Fill(energy, nCell  , m20, GetEventWeight());
  
  fhSMNCell              [matchedPID]->Fill(energy, smMax  , nCell, GetEventWeight());
  fhSMM02NoCut           [matchedPID]->Fill(energy, smMax  , m02  , GetEventWeight());
  
  if ( energy > fEMinShape && energy < fEMaxShape ) 
  {
    fhSMNCellM02[matchedPID]->Fill(smMax , nCell, m02, GetEventWeight());
    
    if ( matchedPID == 0 )
    {
      if ( m02 > 0.1 && m02 < 0.3 ) 
        fhEMaxCellTimeNCellSM->Fill(tmax, smMax, nCell, GetEventWeight());

      fhEMaxCellTimeM02SM->Fill(tmax, smMax, m02, GetEventWeight());

      fhEMaxCellEClusterM02NCellPerSM[smMax]->Fill(maxFrac, nCell, m02, GetEventWeight());
      fhEMaxCellLogM02NCellPerSM     [smMax]->Fill(TMath::Log(eCellMax), nCell, m02, GetEventWeight());
    }
  } // energy bin
  
  if ( nCell > fNCellMinShape ) // it makes sense only for significant size histograms
  { 
    fhSMM02 [matchedPID]->Fill(energy, smMax  , m02, GetEventWeight());
    fhColM02[matchedPID]->Fill(energy, ietaMax, m02, GetEventWeight());
    fhRowM02[matchedPID]->Fill(energy, iphiMax, m02, GetEventWeight());
  }

  // Different shower shape parameters
  if ( fStudyShapeParam )
  {
    // cluster asymmetry
    Float_t dIA    = 1.*(dIphi-dIeta)/(dIeta+dIphi);
    Float_t dIATot = 1.*((dIphiPos-dIphiNeg)-(dIetaPos-dIetaNeg))/((dIetaPos-dIetaNeg)+(dIphiPos-dIphiNeg));
    
    //  Int_t dIphiMin = TMath::Abs(dIphiNeg);
    //  Int_t dIetaMin = TMath::Abs(dIetaNeg);
    //  if(dIphiMin > dIphiPos) dIphiMin = dIphiPos ; 
    //  if(dIetaMin > dIetaPos) dIetaMin = dIetaPos ; 
    //  Float_t dIAMin = 0;
    //  if(dIphiMin > 0 && dIetaMin > 0) dIAMin = 1.*(dIphiMin-dIetaMin)/(dIetaMin+dIphiMin);
    
    AliDebug(1,Form("E %2.2f, nCell %d, dPhi %d, dEta %d, dIA %2.2f, match %d",energy,nCell, dIphi,dIeta,dIA,matchedPID));
    
    //  if(nCell > 5)
    //  {
    //    printf("E %2.2f, nCell %d, dPhi %d, dEta %d, dIA %2.2f, match %d\n",
    //           energy,nCell, dIphi,dIeta,dIA,matchedPID);
    //    printf("\t dPhiNeg %d, dPhiPos %d, dEtaNeg %d, dEtaPos %d, dIATot %2.2f, dIAMin %2.2f\n",
    //           dIphiNeg,dIphiPos,dIetaNeg,dIetaPos,dIATot,dIAMin);
    //  }

    fhDeltaIANCells        [matchedPID]->Fill(energy, nCell  , dIA   , GetEventWeight());
    fhDeltaIATotNCells     [matchedPID]->Fill(energy, nCell  , dIATot, GetEventWeight());
    
    if ( nCell > fNCellMinShape ) // it makes sense only for significant size histograms
    {
      fhDeltaIEtaDeltaIPhi[matchedPID]->Fill(energy, dIeta, dIphi, GetEventWeight());    
      fhDeltaIA           [matchedPID]->Fill(energy, dIA         , GetEventWeight());
      fhDeltaIAM02        [matchedPID]->Fill(energy, m02  , dIA  , GetEventWeight());
      fhDeltaIAM20        [matchedPID]->Fill(energy, m20  , dIA  , GetEventWeight());
      
      fhDeltaIEtaDeltaIPhiTot[matchedPID]->Fill(energy,dIetaPos-dIetaNeg, dIphiPos-dIphiNeg, GetEventWeight());    
      fhDeltaIATot           [matchedPID]->Fill(energy, dIATot       , GetEventWeight());
      fhDeltaIATotM02        [matchedPID]->Fill(energy, m02  , dIATot, GetEventWeight());
      fhDeltaIATotM20        [matchedPID]->Fill(energy, m20  , dIATot, GetEventWeight());
    }
    
    if ( IsDataMC() && mcIndex > -1 && mcIndex < 10 && nCell > fNCellMinShape )
    {
      fhDeltaIAOrigin   [matchedPID]->Fill(energy, mcIndex, dIA   , GetEventWeight()); 
      fhDeltaIATotOrigin[matchedPID]->Fill(energy, mcIndex, dIATot, GetEventWeight()); 
    }
    
    Float_t l0   = 0., l1   = 0.;
    Float_t dispp= 0., dEta = 0., dPhi    = 0.;
    Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
    if ( GetCalorimeter() == kEMCAL )
    {
      GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), fCaloCellList, clus,
                                                                                   l0, l1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);
      
      Float_t sEtaPhiA = -1000.;
      if(sEta+sPhi>0.0001) sEtaPhiA = (sPhi-sEta)/(sEta+sPhi);
      
      AliDebug(2,Form("Recalculate shower shape org: m02 %2.2f, m20 %2.2f, disp %2.2f;"
                      " new: m02 %2.2f, m20 %2.2f, disp %2.2f; "
                      "mEta %2.2f, mPhi %2.2f, mEtaPhi %2.2f, A_EtaPhi %2.2f; dEta %2.2f dPhi %2.2f",
                      m02,m20,clus->GetDispersion(),l0,l1,dispp,sEta,sPhi,sEtaPhi,sEtaPhiA,dEta,dPhi));
      
      fhNCellsPerClusterMEta    [matchedPID]->Fill(energy, nCell, sEta    , GetEventWeight());
      fhNCellsPerClusterMPhi    [matchedPID]->Fill(energy, nCell, sPhi    , GetEventWeight());
      fhNCellsPerClusterMEtaPhi [matchedPID]->Fill(energy, nCell, sEtaPhi , GetEventWeight());
      fhNCellsPerClusterMEtaPhiA[matchedPID]->Fill(energy, nCell, sEtaPhiA, GetEventWeight());
    }
  }
  
  // Invariant mass for clusters looking like photons, depending number of cells
  if ( matchedPID == 0 && energy > fEMinShape && energy < fEMaxShape &&
       m02 > fInvMassMinM02Cut && m02 < fInvMassMaxM02Cut )
  {
    for(Int_t jclus = 0 ; jclus < fCaloClusList->GetEntriesFast() ; jclus++) 
    {
      AliVCluster* clus2 =  (AliVCluster*) fCaloClusList->At(jclus);
      
      Float_t maxCellFraction = 0.;
      Int_t absIdMax2 = GetCaloUtils()->GetMaxEnergyCell(fCaloCellList, clus2, maxCellFraction);
      
      Double_t tof2 =  clus2->GetTOF()*1.e9;
      if(tof2>400) tof2-=fConstantTimeShift;
      
      Double_t diffTof = tmax-tof2;
      
      // Try to reduce background with a mild shower shape cut and no more 
      // than 1 local maximum in cluster and remove low energy clusters
      
      if(   absIdMax == absIdMax2   
         || !IsGoodCluster(absIdMax2, clus2->GetM02(), clus2->GetNCells()) 
         || GetCaloUtils()->GetNumberOfLocalMaxima(clus2,fCaloCellList) > 1 
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
      fhInvMassNCellSM->Fill(mass, nCell, smMax, GetEventWeight());
      
      Int_t smMax2 = GetModuleNumber(clus2);
      if(smMax == smMax2) fhInvMassNCellSMSame->Fill(mass, nCell, smMax, GetEventWeight());
    }
  }
  
  // Check the origin.
  if ( IsDataMC() && mcIndex > -1 && mcIndex < 10)
  {
    fhOriginE  [matchedPID]->Fill(energy, mcIndex,      GetEventWeight());
    fhOriginM02[matchedPID]->Fill(energy, mcIndex, m02, GetEventWeight());    
  } // MC
}

//____________________________________________________________________________
/// Check if cluster is matched to a track and do a track ID
/// either electrons or hadrons. 
/// filling different type of histograms:
///
/// \param clus: AliVCluster pointer
/// \param matchedPID: 0-neutral, 1-electron, 2-hadron
//____________________________________________________________________________
void  AliAnaClusterShapeCorrelStudies::ClusterMatchedToTrackPID
(AliVCluster *clus, Int_t & matchedPID)
{
  
  AliVTrack *track = GetCaloUtils()->GetMatchedTrack(clus, GetReader()->GetInputEvent());
  
  if(!track) 
  {
    matchedPID = -1;
    return ;
  }
  
//  Double_t tpt   = track->Pt();
//  Double_t tmom  = track->P();
  Double_t dedx  = track->GetTPCsignal();
//  Int_t    nITS  = track->GetNcls(0);
//  Int_t    nTPC  = track->GetNcls(1);
//  Bool_t positive = kFALSE;
//  if(track) positive = (track->Charge()>0);
  
  // Residuals
//  Float_t deta  = clus->GetTrackDz();
//  Float_t dphi  = clus->GetTrackDx();
//  Double_t  dR  = TMath::Sqrt(dphi*dphi + deta*deta);
//  Int_t nModule = GetModuleNumber(clus);

  // Electron or else?
  
  // Init at least once
  if(fdEdXMinEle == 0 || fdEdXMaxEle == 0 || fdEdXMinHad == 0 || fdEdXMaxHad == 0) 
    InitdEdXParameters();
  
  if      ( dedx >= fdEdXMinEle && dedx < fdEdXMaxEle ) matchedPID = 1; 
  else if ( dedx >= fdEdXMinHad && dedx < fdEdXMaxHad ) matchedPID = 2;  
  else
  {
    AliDebug(1,Form("dEdX out of range %2.2f",dedx));
    matchedPID = -1;
  }
}

//____________________________________________________________________________
/// Fill clusters related histograms, execute here the loop of clusters
/// apply basic selection cuts (track matching, gooness, exoticity, timing)
/// and the call to the different methods
/// filling different type of histograms:
/// * Cluster Asymmetry
/// * Dependence of m02 and energy on different parameters
/// * Effect of different weight in clusters
/// * Cluster shape vs ncells, in same or different T-Card 
//____________________________________________________________________________
void AliAnaClusterShapeCorrelStudies::ClusterLoopHistograms()
{  
  Int_t  nCaloClusters         = fCaloClusList->GetEntriesFast() ;
  Int_t  nCaloCellsPerCluster  = 0  ;
  Bool_t matched               = kFALSE;
  Int_t  nModule               =-1  ;
  
  // Get vertex for photon momentum calculation and event selection
  Double_t v[3] = {0,0,0}; //vertex ;
                           //GetReader()->GetVertex(v);
  
  AliDebug(1,Form("In %s there are %d clusters", GetCalorimeterString().Data(), nCaloClusters));
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++)
  {
    AliDebug(1,Form("Cluster: %d/%d, data %d",iclus+1,nCaloClusters,GetReader()->GetDataType()));
    
    AliVCluster* clus =  (AliVCluster*) fCaloClusList->At(iclus);
        
    // away from dead region
    if ( clus->GetDistanceToBadChannel() < fMinDistToBad ) return ;  

    // SuperModule number of cluster
    nModule = GetModuleNumber(clus);
    if ( nModule < fFirstModule || nModule > fLastModule ) 
    {
      AliDebug(1,Form("Cluster module out of range %d",nModule));
      continue ;
    }
    
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
    //
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(fCaloCellList, clus, maxCellFraction);
    
    // Cut on time of clusters
    Double_t tof = clus->GetTOF()*1.e9;
    if(tof>400) tof-=fConstantTimeShift;
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    Float_t e   = fClusterMomentum.E();
    Float_t pt  = fClusterMomentum.Pt();
    Float_t eta = fClusterMomentum.Eta();
    Float_t phi = GetPhi(fClusterMomentum.Phi());
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if(IsFiducialCutOn()) 
      in = GetFiducialCut()->IsInFiducialCut(fClusterMomentum.Eta(),fClusterMomentum.Phi(),GetCalorimeter()) ;
    
    if(!in)
    {
      AliDebug(1,Form("Remove cluster with phi %2.2f and eta %2.2f", phi*TMath::RadToDeg(), eta));
      continue;
    }
    
    AliDebug(1,Form("cluster: E %2.3f, pT %2.3f, eta %2.3f, phi %2.3f",e,pt,eta,phi*TMath::RadToDeg()));
    
    // Select the cluster
    //
    nCaloCellsPerCluster = clus->GetNCells();

    Bool_t goodCluster = IsGoodCluster(absIdMax, clus->GetM02(), nCaloCellsPerCluster);
    
    AliDebug(1,Form("Accept cluster? %d",goodCluster));
    
    if(!goodCluster) continue;

    // MC origin finding
    //
    Int_t mcTag   =  0;
    Int_t mcIndex = -1;
    if ( IsDataMC() && fStudyShape )
    {
      mcTag = GetMCAnalysisUtils()->CheckOrigin(clus->GetLabels(), clus->GetNLabels(), GetMC());
      
      if      ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPi0        ) ||
                GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCEta        ) ) mcIndex = 0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPi0Decay   ) ||
                GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCEtaDecay   ) ||
                GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCOtherDecay ) ) mcIndex = 1;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPhoton     ) ) mcIndex = 2;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCElectron   ) ) mcIndex = 3;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPion       ) ) mcIndex = 4;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCKaon       ) ) mcIndex = 5;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCProton     ) ) mcIndex = 6;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCAntiProton ) ) mcIndex = 7;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCNeutron    ) ) mcIndex = 8;
      else if ( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCAntiNeutron) ) mcIndex = 9;
    }

    // Cluster mathed with track? and what kind?
    //
    matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(), GetReader()->GetInputEvent());
 
    Int_t matchedPID = 0;
    if ( matched && fStudyShape )
      ClusterMatchedToTrackPID(clus, matchedPID);
    
    // Get amp and time of max cell, recalibrate and calculate things
    //
    Int_t    bc     = (GetReader()->GetInputEvent())->GetBunchCrossNumber();
    Double_t tmax   = fCaloCellList->GetCellTime(absIdMax);
    Float_t  ampMax = fCaloCellList->GetCellAmplitude(absIdMax);

    GetCaloUtils()->RecalibrateCellTime(tmax, GetCalorimeter(), absIdMax, bc);
    tmax*=1.e9;
    if(tmax>400) tmax-=fConstantTimeShift;
    
    GetCaloUtils()->RecalibrateCellAmplitude(ampMax, GetCalorimeter(), absIdMax);
    
    Float_t eCrossFrac = 0;
    if ( ampMax > 0.01 ) 
      eCrossFrac = 1-GetCaloUtils()->GetECross(absIdMax,fCaloCellList,bc)/ampMax;
    
    
    // Call analysis method filling histograms
    //
  
    //
    if ( fStudyShape  && matchedPID >= 0 && matchedPID < 3 )
      ClusterShapeHistograms(clus, absIdMax, maxCellFraction, eCrossFrac, ampMax, tmax, matchedPID, mcIndex);
    
    //
    if ( fStudyTCardCorrelation ) 
      ChannelCorrelationInTCard(clus, matched, absIdMax, eCrossFrac);
    
    // 
    if ( fStudyWeight ) 
      WeightHistograms(clus, mcTag);
    
  } // Cluster loop
  
}

//_________________________________________________
/// Save parameters used for analysis in a string.
//_________________________________________________
TObjString * AliAnaClusterShapeCorrelStudies::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaClusterShapeCorrelStudies ---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Cluster M02: > %2.2f ; n cells > %d;  dist to bad>%2.1f;",fM02Min, fNCellMin, fMinDistToBad) ;
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
TList * AliAnaClusterShapeCorrelStudies::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ClusterShapeStudies") ; 
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  //Int_t totalSM = fLastModule-fFirstModule+1;
  //printf("N SM %d, first SM %d, last SM %d, total %d\n",fNModules,fFirstModule,fLastModule, totalSM);

  
  // MC origin
  TString mcParticleStringLabel[] = {"Merged #gamma#gamma","Decay #gamma","Direct #gamma","e^{#pm}","#pi^{#pm}","k^{#pm}","p","#bar{p}","n","#bar{n}"};
  //TString mcParticleStringTitle[] = {"MergedPhoton","DecayPhoton","DirectPhoton","Electron","Pion","Kaon","Proton","AntiProton","Neutron","AntiNeutron"};

  // track-match PID matching
  TString matchCase[] = {"Neutral","Electron","Hadron"};
  
  // Histogram binning and ranges
  // 
  Int_t nptbins     = GetHistogramRanges()->GetHistoPtBins(); 	        Float_t ptmax     = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin     = GetHistogramRanges()->GetHistoPtMin();
  Int_t nmassbins   = GetHistogramRanges()->GetHistoMassBins();         Float_t massmax   = GetHistogramRanges()->GetHistoMassMax(); 	       Float_t massmin   = GetHistogramRanges()->GetHistoMassMin();
 
  Int_t ntimebins   = GetHistogramRanges()->GetHistoTimeBins();         Float_t timemax   = GetHistogramRanges()->GetHistoTimeMax();         Float_t timemin   = GetHistogramRanges()->GetHistoTimeMin();       
  Int_t tdbins      = GetHistogramRanges()->GetHistoDiffTimeBins() ;    Float_t tdmax     = GetHistogramRanges()->GetHistoDiffTimeMax();     Float_t tdmin     = GetHistogramRanges()->GetHistoDiffTimeMin();

  Int_t nceclbins   = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nceclmax  = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nceclmin  = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  Int_t ssbins      = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax     = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin     = GetHistogramRanges()->GetHistoShowerShapeMin();
  
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
  
  // E bins in TH3
  Int_t nEbins = 16  ;
  Float_t minE =  2.5;
  Float_t maxE = 18.5;
  
  // shower shape bins in TH3
  Int_t nShShBins = 200;
  Float_t minShSh = 0.;
  Float_t maxShSh = 2.;
  
  // Asymmetry bins
  Int_t asyBins  = 21;
  Float_t asyMax = 1.05;
  Float_t asyMin = -1*asyMax;
  
  // n cell bins for TH3
  Int_t cellBins  = 15;
  Float_t cellMax = 15;
  Float_t cellMin = 0;
  
  //
  // Init histograms
  //
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
      
      fhTimeTCardCorrNoSelection[tm]  = new TH2F 
      (Form("hTimeTCardCorrNoSelection%s",add[tm].Data()),
       Form("#it{time} vs #it{E} %s",add[tm].Data()),
       nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
      fhTimeTCardCorrNoSelection[tm]->SetXTitle("#it{E} (GeV)");
      fhTimeTCardCorrNoSelection[tm]->SetYTitle("#it{time} (ns)");
      outputContainer->Add(fhTimeTCardCorrNoSelection[tm]); 
      
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
       nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
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
       nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
      fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]->SetXTitle("#it{E}_{cell max} (GeV)");
      fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]->SetYTitle("#it{t}_{cell}^{same TCard}-#it{t}_{cell}^{diff TCard} (ns)");
      outputContainer->Add(fhSameRowDiffColAndTCardCellsTimeDiffCellMaxE[tm]); 
      
      if(fStudyExotic)
      {
        fhEnergyTime1Cell[tm]  = new TH2F 
        (Form("hEnergyTime1Cell%s",add[tm].Data()),
         Form("#it{t} vs #it{E}, 1 cells cluster %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhEnergyTime1Cell[tm]->SetXTitle("#it{E} (GeV)");
        fhEnergyTime1Cell[tm]->SetYTitle("#it{t} (ns)");
        outputContainer->Add(fhEnergyTime1Cell[tm]); 
        
        fhEnergyTimeExotic[tm]  = new TH2F 
        (Form("hEnergyTimeExotic%s",add[tm].Data()),
         Form("#it{t} vs #it{E},  exo > 0.97, %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhEnergyTimeExotic[tm]->SetXTitle("#it{E} (GeV)");
        fhEnergyTimeExotic[tm]->SetYTitle("#it{t} (ns)");
        outputContainer->Add(fhEnergyTimeExotic[tm]); 
        
        fhEnergyTimeTCardCorrNoSelection1Cell[tm]  = new TH2F 
        (Form("hEnergyTimeTCardCorrNoSelection1Cell%s",add[tm].Data()),
         Form("#it{t} vs #it{E}, 1 cells cluster %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
        fhEnergyTimeTCardCorrNoSelection1Cell[tm]->SetXTitle("#it{E} (GeV)");
        fhEnergyTimeTCardCorrNoSelection1Cell[tm]->SetYTitle("#it{t} (ns)");
        outputContainer->Add(fhEnergyTimeTCardCorrNoSelection1Cell[tm]); 
        
        fhEnergyTimeTCardCorrNoSelectionExotic[tm]  = new TH2F 
        (Form("hEnergyTimeTCardCorrNoSelectionExotic%s",add[tm].Data()),
         Form("#it{t} vs #it{E},  exo > 0.97, %s",add[tm].Data()),
         nptbins,ptmin,ptmax,ntimebins,timemin,timemax); 
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
         nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
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
         nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
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
           nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
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
             nptbins,ptmin,ptmax,tdbins,tdmin,tdmax); 
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
  
  // Cluster size in terms of cells and shape TH3
  if(fStudyShape)
  {
    fhColRowM02 = new TH3F
    (Form("hColRowM02"),
     Form("column vs row vs M02, %2.2f < #it{E} < %2.2f GeV for Neutral",fEMinShape,fEMaxShape),
     ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,40,0.,2.);
    fhColRowM02->SetYTitle("row");
    fhColRowM02->SetXTitle("column");
    fhColRowM02->SetZTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhColRowM02) ;
    
    fhColRowM02NCellCut = new TH3F
    (Form("hColRowM02NCellCut"),
     Form("column vs row vs M02, %2.2f<#it{E}<%2.2f GeV #it{n}_{cells}^{w>0.01} > %d for Neutral",
          fEMinShape,fEMaxShape,fNCellMinShape),
     ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax,40,0.,2.);
    fhColRowM02NCellCut->SetYTitle("row");
    fhColRowM02NCellCut->SetXTitle("column");
    fhColRowM02NCellCut->SetZTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhColRowM02NCellCut) ;

    fhInvMassNCellSM  = new TH3F 
    ("hInvMassNCellSM",
     Form("%2.2f<#it{E}_{1}<%2.2f GeV, %2.2f<#it{E}_{2}<%2.2f GeV, %2.2f<#lambda^{2}_{0}<%2.2f"
          "#it{M}_{#gamma #gamma} vs #it{n}_{1, cells}^{w>0.01} vs SM number trig cluster",
          fEMinShape,fEMaxShape,fInvMassMinECut,fInvMassMaxECut,fInvMassMinM02Cut,fInvMassMaxM02Cut),
     nmassbins,massmin,massmax,cellBins,cellMin,cellMax,fNModules,-0.5,fNModules-0.5); 
    fhInvMassNCellSM->SetZTitle("SM number");
    fhInvMassNCellSM->SetYTitle("#it{n}_{cells}^{w>0.01}");
    fhInvMassNCellSM->SetXTitle("#it{M}_{#gamma #gamma}");
    outputContainer->Add(fhInvMassNCellSM);         
    
    fhInvMassNCellSMSame  = new TH3F 
    ("hInvMassNCellSMSame",
     Form("%2.2f<#it{E}_{1}<%2.2f GeV, %2.2f<#it{E}_{2}<%2.2f GeV, %2.2f<#lambda^{2}_{0}<%2.2f"
          "#it{M}_{#gamma #gamma} vs #it{n}_{1, cells}^{w>0.01} vs SM number both cluster",
          fEMinShape,fEMaxShape,fInvMassMinECut,fInvMassMaxECut,fInvMassMinM02Cut,fInvMassMaxM02Cut),
     nmassbins,massmin,massmax,cellBins,cellMin,cellMax,fNModules,-0.5,fNModules-0.5); 
    fhInvMassNCellSMSame->SetZTitle("SM number");
    fhInvMassNCellSMSame->SetYTitle("#it{n}_{cells}^{w>0.01}");
    fhInvMassNCellSMSame->SetXTitle("#it{M}_{#gamma #gamma}");
    outputContainer->Add(fhInvMassNCellSMSame);       
    
    //
        
    fhEMaxCellTimeM02SM  = new TH3F 
    ("hEMaxCellTimeM02SM",
     Form("time vs #lambda_{0}^{2} vs SM number, "
          "%2.2f<#it{E}<%2.2f GeV, #it{n}_{cells}^{w>0.01}>%d",fEMinShape,fEMaxShape,fNCellMinShape),
     45,-25.5,20.5, fNModules,-0.5,fNModules-0.5,nShShBins,minShSh,maxShSh); 
    fhEMaxCellTimeM02SM->SetZTitle("#lambda_{0}^{2}");
    fhEMaxCellTimeM02SM->SetYTitle("SM number");
    fhEMaxCellTimeM02SM->SetXTitle("time (ns)");
    outputContainer->Add(fhEMaxCellTimeM02SM);  
    
    fhEMaxCellTimeNCellSM  = new TH3F 
    ("hEMaxCellTimeNCellSM",
     Form("time vs #it{n}_{cells}^{w>0.01} vs SM number, "
          "%2.2f<#it{E}<%2.2f GeV, 0.1<#it{sigma}_{long}<0.3",fEMinShape,fEMaxShape),
     45,-25.5,20.5, fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax); 
    fhEMaxCellTimeNCellSM->SetZTitle("#it{n}_{cells}^{w>0.01}");
    fhEMaxCellTimeNCellSM->SetYTitle("SM number");
    fhEMaxCellTimeNCellSM->SetXTitle("time (ns)");
    outputContainer->Add(fhEMaxCellTimeNCellSM);  

    fhESecCellTimeNCellSM  = new TH3F 
    ("hESecCellTimeNCellSM",
     Form("secondary cell time vs #it{n}_{cells}^{w>0.01} vs SM number, "
          "%2.2f<#it{E}<%2.2f GeV, 0.1<#it{sigma}_{long}<0.3",fEMinShape,fEMaxShape),
     50,-100,100, fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax); 
    fhESecCellTimeNCellSM->SetZTitle("#it{n}_{cells}^{w>0.01}");
    fhESecCellTimeNCellSM->SetYTitle("SM number");
    fhESecCellTimeNCellSM->SetXTitle("time (ns)");
    outputContainer->Add(fhESecCellTimeNCellSM);  
    
    //
    
    Int_t   nbinsdeltacells = 19 ;
    Float_t   mindeltacells =-9.5;
    Float_t   maxdeltacells = 9.5;
    for(Int_t col = 0; col < 2; col++)
    {
      fhColRowFromCellMaxLowM02[col]  = new TH2F 
      (Form("hColRowFromCellMaxLowM02_Col%d",col),
       Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{n}_{cells}^{w>0.01} for 0.1 < #lambda_{0}^{2} < 0.3, "
            "%2.2f<#it{E}<%2.2f GeV, colum %d",fEMinShape,fEMaxShape,col),
       nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells); 
      fhColRowFromCellMaxLowM02[col]->SetXTitle("#Delta column_{max-secondary}");
      fhColRowFromCellMaxLowM02[col]->SetYTitle("#Delta row_{max-secondary}");
      outputContainer->Add(fhColRowFromCellMaxLowM02[col]);   
            
      fhColRowFromCellMaxHighM02[col]  = new TH2F 
      (Form("hColRowFromCellMaxHighM02_Col%d",col),
       Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{n}_{cells}^{w>0.01} for 0.5 < #lambda_{0}^{2} < 2, "
            "%2.2f<#it{E}<%2.2f GeV, colum %d",fEMinShape,fEMaxShape,col),
       nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells); 
      fhColRowFromCellMaxHighM02[col]->SetXTitle("#Delta column_{max-secondary}");
      fhColRowFromCellMaxHighM02[col]->SetYTitle("#Delta row_{max-secondary}");
      outputContainer->Add(fhColRowFromCellMaxHighM02[col]);        
    }//odd/pair col

    //
    
    for(Int_t i = 0; i < fNModules; i++)
    {
      
      fhNCellsPerClusterM02M20PerSM [i] = new TH3F
      (Form("hNCellsPerClusterM02M20_SM%d",i),
       Form(" vs #lambda_{0}^{2} vs #it{n}_{cells}^{w>0.01} vs #lambda_{1}^{2}, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       nShShBins/2,minShSh,maxShSh/2, cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhNCellsPerClusterM02M20PerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhNCellsPerClusterM02M20PerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhNCellsPerClusterM02M20PerSM[i]->SetXTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhNCellsPerClusterM02M20PerSM[i]);  

      fhEMaxCellEClusterM02NCellPerSM[i]  = new TH3F 
      (Form("hEMaxCellEClusterM02NCell_SM%d",i),
       Form("(#it{E}_{cluster} - #it{E}_{max cell})/#it{E}_{cluster} vs #lambda_{0}^{2} vs #it{n}_{cells}^{w>0.01}, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       50,0,1., cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhEMaxCellEClusterM02NCellPerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhEMaxCellEClusterM02NCellPerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxCellEClusterM02NCellPerSM[i]->SetXTitle("(#it{E}_{cluster} - #it{E}_{max cell})/ #it{E}_{cluster}");
      outputContainer->Add(fhEMaxCellEClusterM02NCellPerSM[i]);  
      
      fhEMaxCellLogM02NCellPerSM[i]  = new TH3F 
      (Form("hEMaxCellLogM02NCell_SM%d",i),
       Form("log(#it{E}_{max cell}) vs #lambda_{0}^{2} vs SM number, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       150,-3,3, cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhEMaxCellLogM02NCellPerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhEMaxCellLogM02NCellPerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxCellLogM02NCellPerSM[i]->SetXTitle("log(#it{E}_{max cell})");
      outputContainer->Add(fhEMaxCellLogM02NCellPerSM[i]);  
      
      fhESecCellEMaxCellM02NCellPerSM[i]  = new TH3F 
      (Form("hESecCellEMaxCellM02NCell_SM%d",i),
       Form("#it{E}_{cell}/#it{E}_{cell max} vs #lambda_{0}^{2} vs #it{n}_{cells}^{w>0.01}, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       50,0,1., cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhESecCellEMaxCellM02NCellPerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhESecCellEMaxCellM02NCellPerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhESecCellEMaxCellM02NCellPerSM[i]->SetXTitle("#it{E}_{cell}/#it{E}_{cell max}");
      outputContainer->Add(fhESecCellEMaxCellM02NCellPerSM[i]);  
      
      fhESecCellEClusterM02NCellPerSM[i]  = new TH3F 
      (Form("hESecCellEClusterM02NCell_SM%d",i),
       Form("(#it{E}_{cluster} - #it{E}_{cell})/#it{E}_{cluster} vs #lambda_{0}^{2} vs #it{n}_{cells}^{w>0.01}, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       50,0,1., cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhESecCellEClusterM02NCellPerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhESecCellEClusterM02NCellPerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhESecCellEClusterM02NCellPerSM[i]->SetXTitle("(#it{E}_{cluster} - #it{E}_{cell})/ #it{E}_{cluster}");
      outputContainer->Add(fhESecCellEClusterM02NCellPerSM[i]);  
      
      fhESecCellLogM02NCellPerSM[i]  = new TH3F 
      (Form("hESecCellLogM02NCell_SM%d",i),
       Form("log(#it{E}_{cell}) vs #lambda_{0}^{2} vs #it{n}_{cells}^{w>0.01}, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,fNCellMinShape),
       150,-3,3, cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhESecCellLogM02NCellPerSM[i]->SetZTitle("#lambda_{0}^{2}");
      fhESecCellLogM02NCellPerSM[i]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhESecCellLogM02NCellPerSM[i]->SetXTitle("log(#it{E}_{cell})");
      outputContainer->Add(fhESecCellLogM02NCellPerSM[i]);  
      
      fhEMaxESecCellNCellLowM02PerSM[i]  = new TH3F 
      (Form("hEMaxESecCellNCellLowM02_SM%d",i),
       Form("#it{E}_{cell}^{max} vs #it{E}_{cell}^{secondary} vs #it{n}_{cells}^{w>0.01} vs 0.1 < #lambda_{0}^{2} < 0.3, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       24,0,12, 48,0,12,cellBins,cellMin,cellMax); 
      fhEMaxESecCellNCellLowM02PerSM[i]->SetZTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxESecCellNCellLowM02PerSM[i]->SetYTitle("#it{E}_{cell}^{secondary} (GeV)");
      fhEMaxESecCellNCellLowM02PerSM[i]->SetXTitle("#it{E}_{cell}^{max} (GeV)");
      outputContainer->Add(fhEMaxESecCellNCellLowM02PerSM[i]);      
      
      fhEMaxECrossNCellLowM02PerSM[i]  = new TH3F 
      (Form("hEMaxECrossNCellLowM02_SM%d",i),
       Form("#it{E}_{cell}^{max} vs exoticity vs #it{n}_{cells}^{w>0.01} vs 0.1 < #lambda_{0}^{2} < 0.3, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       24,0,12, 40,0.6,1.,cellBins,cellMin,cellMax); 
      fhEMaxECrossNCellLowM02PerSM[i]->SetZTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxECrossNCellLowM02PerSM[i]->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell}^{max}");
      fhEMaxECrossNCellLowM02PerSM[i]->SetXTitle("#it{E}_{cell}^{max} (GeV)");
      outputContainer->Add(fhEMaxECrossNCellLowM02PerSM[i]); 
      
      fhEMaxESecCellNCellHighM02PerSM[i]  = new TH3F 
      (Form("hEMaxESecCellNCellHighM02_SM%d",i),
       Form("#it{E}_{cell}^{max} vs #it{E}_{cell}^{secondary} vs #it{n}_{cells}^{w>0.01} vs 0.5 < #lambda_{0}^{2} < 2, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       24,0,12, 48,0,12,cellBins,cellMin,cellMax); 
      fhEMaxESecCellNCellHighM02PerSM[i]->SetZTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxESecCellNCellHighM02PerSM[i]->SetYTitle("#it{E}_{cell}^{secondary} (GeV)");
      fhEMaxESecCellNCellHighM02PerSM[i]->SetXTitle("#it{E}_{cell}^{max} (GeV)");
      outputContainer->Add(fhEMaxESecCellNCellHighM02PerSM[i]);      
      
      fhEMaxECrossNCellHighM02PerSM[i]  = new TH3F 
      (Form("hEMaxECrossNCellHighM02_SM%d",i),
       Form("#it{E}_{cell}^{max} vs exoticity vs #it{n}_{cells}^{w>0.01} for 0.5 < #lambda_{0}^{2} < 2, "
            "%2.2f<#it{E}<%2.2f GeV, SM=%d",fEMinShape,fEMaxShape,i),
       24,0,12, 40,0.6,1.,cellBins,cellMin,cellMax); 
      fhEMaxECrossNCellHighM02PerSM[i]->SetZTitle("#it{n}_{cells}^{w>0.01}");
      fhEMaxECrossNCellHighM02PerSM[i]->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell}^{max}");
      fhEMaxECrossNCellHighM02PerSM[i]->SetXTitle("#it{E}_{cell}^{max} (GeV)");
      outputContainer->Add(fhEMaxECrossNCellHighM02PerSM[i]);    

      Int_t   nbinsdeltacells = 19 ;
      Float_t   mindeltacells =-9.5;
      Float_t   maxdeltacells = 9.5;
      for(Int_t col = 0; col < 2; col++)
      {
        fhColRowFromCellMaxLowM02PerSM[i][col]  = new TH3F 
        (Form("hColRowFromCellMaxLowM02_SM%d_Col%d",i,col),
         Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{n}_{cells}^{w>0.01} for 0.1 < #lambda_{0}^{2} < 0.3, "
              "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d",fEMinShape,fEMaxShape,i,col),
         nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,cellBins,cellMin,cellMax); 
        fhColRowFromCellMaxLowM02PerSM[i][col]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhColRowFromCellMaxLowM02PerSM[i][col]->SetXTitle("#Delta column_{max-secondary}");
        fhColRowFromCellMaxLowM02PerSM[i][col]->SetYTitle("#Delta row_{max-secondary}");
        outputContainer->Add(fhColRowFromCellMaxLowM02PerSM[i][col]);   
        
        fhColRowFromCellMaxHighM02PerSM[i][col]  = new TH3F 
        (Form("hColRowFromCellMaxHighM02_SM%d_Col%d",i,col),
         Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{n}_{cells}^{w>0.01} for 0.5 < #lambda_{0}^{2} < 2, "
              "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d",fEMinShape,fEMaxShape,i,col),
         nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,cellBins,cellMin,cellMax); 
        fhColRowFromCellMaxHighM02PerSM[i][col]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhColRowFromCellMaxHighM02PerSM[i][col]->SetXTitle("#Delta column_{max-secondary}");
        fhColRowFromCellMaxHighM02PerSM[i][col]->SetYTitle("#Delta row_{max-secondary}");
        outputContainer->Add(fhColRowFromCellMaxHighM02PerSM[i][col]);   
        
        for(Int_t j = 0; j < 3; j++)
        {
          fhColRowFromCellMaxEMaxSecDiffLowM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxEMaxSecDiffLowM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #Delta #it{E}_{max-secondary} for 0.1 < #lambda_{0}^{2} < 0.3, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,120,0,12); 
          fhColRowFromCellMaxEMaxSecDiffLowM02PerSM[i][col][j]->SetZTitle("#Delta #it{E}_{max-secondary} (GeV)");
          fhColRowFromCellMaxEMaxSecDiffLowM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxEMaxSecDiffLowM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxEMaxSecDiffLowM02PerSM[i][col][j]);   
          
          fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxEMaxSecDiffHighM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #Delta #it{E}_{max-secondary} for 0.5 < #lambda_{0}^{2} < 2, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,120,0,12); 
          fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][col][j]->SetZTitle("#Delta #it{E}_{max-secondary} (GeV)");
          fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxEMaxSecDiffHighM02PerSM[i][col][j]);   

          fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxEMaxSecDiffFracLowM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #Delta #it{E}_{max-secondary}/#it{E}_{max} for 0.1 < #lambda_{0}^{2} < 0.3, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,40,0,1); 
          fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[i][col][j]->SetZTitle("#Delta #it{E}_{max-secondary}/#it{E}_{max}");
          fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxEMaxSecDiffFracLowM02PerSM[i][col][j]);   
          
          fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxEMaxSecDiffFracHighM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #Delta #it{E}_{max-secondary}/#it{E}_{max} for 0.5 < #lambda_{0}^{2} < 2, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,40,0,1); 
          fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][col][j]->SetZTitle("#Delta #it{E}_{max-secondary}/#it{E}_{max}");
          fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxEMaxSecDiffFracHighM02PerSM[i][col][j]);   
          
          fhColRowFromCellMaxECellClusterRatLowM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxECellClusterRatLowM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{E}_{cell}/#it{E}_{cluster} for 0.1 < #lambda_{0}^{2} < 0.3, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,50,0,1); 
          fhColRowFromCellMaxECellClusterRatLowM02PerSM[i][col][j]->SetZTitle("#it{E}_{cell}/#it{E}_{cluster}");
          fhColRowFromCellMaxECellClusterRatLowM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxECellClusterRatLowM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxECellClusterRatLowM02PerSM[i][col][j]);   
          
          fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][col][j]  = new TH3F 
          (Form("hColRowFromCellMaxECellClusterRatHighM02_SM%d_Col%d_NCellBin%d",i,col,j),
           Form("cell col_{max}-#col_{secondary} vs cell #row_{max}-#row_{secondary} vs #it{E}_{cell}/#it{E}_{cluster} for 0.5 < #lambda_{0}^{2} < 2, "
                "%2.2f<#it{E}<%2.2f GeV, SM=%d, colum %d, #it{n}_{cells}^{w>0.01} bin %d",fEMinShape,fEMaxShape,i,col,j),
           nbinsdeltacells,mindeltacells,maxdeltacells,nbinsdeltacells,mindeltacells,maxdeltacells,50,0,1); 
          fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][col][j]->SetZTitle("#it{E}_{cell}/#it{E}_{cluster}");
          fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][col][j]->SetXTitle("#Delta column_{max-secondary}");
          fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][col][j]->SetYTitle("#Delta row_{max-secondary}");
          outputContainer->Add(fhColRowFromCellMaxECellClusterRatHighM02PerSM[i][col][j]);   

        } // 3 n cell bins
      }//odd/pair col
  
    } // SM
    
    for(Int_t imatch = 0; imatch < 3; imatch++)
    {      
      if ( fStudyShapeParam )
      {
        fhDeltaIEtaDeltaIPhi[imatch]  = new TH3F 
        (Form("hDeltaIEtaDeltaIPhi_%s",matchCase[imatch].Data()),
         Form("Cluster max size  with respect main cell in columns vs rows vs E for %s, #it{n}_{cells}^{w>0.01}>%d",matchCase[imatch].Data(),fNCellMinShape),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,cellBins,cellMin,cellMax); 
        fhDeltaIEtaDeltaIPhi[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIEtaDeltaIPhi[imatch]->SetYTitle("#Delta Column");
        fhDeltaIEtaDeltaIPhi[imatch]->SetZTitle("#Delta Row");
        outputContainer->Add(fhDeltaIEtaDeltaIPhi[imatch]); 
        
        fhDeltaIEtaDeltaIPhiTot[imatch]  = new TH3F 
        (Form("hDeltaIEtaDeltaIPhiTot_%s",matchCase[imatch].Data()),
         Form("Cluster size in columns vs rows, minus main cell, vs E for %s, #it{n}_{cells}^{w>0.01}>%d",matchCase[imatch].Data(),fNCellMinShape),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,cellBins,cellMin,cellMax); 
        fhDeltaIEtaDeltaIPhiTot[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIEtaDeltaIPhiTot[imatch]->SetYTitle("#Delta Column");
        fhDeltaIEtaDeltaIPhiTot[imatch]->SetZTitle("#Delta Row");
        outputContainer->Add(fhDeltaIEtaDeltaIPhiTot[imatch]); 
        
        fhDeltaIA[imatch]  = new TH2F
        (Form("hDeltaIA_%s",matchCase[imatch].Data()),
         Form("Cluster *asymmetry* in cell units vs E for %s",matchCase[imatch].Data()),
         nptbins,ptmin,ptmax,asyBins,asyMin,asyMax); 
        fhDeltaIA[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIA[imatch]->SetYTitle("#it{A}_{cell in cluster}");
        outputContainer->Add(fhDeltaIA[imatch]); 
        
        fhDeltaIATot[imatch]  = new TH2F
        (Form("hDeltaIATot_%s",matchCase[imatch].Data()),
         Form("Cluster *total asymmetry* in cell units vs E for %s",matchCase[imatch].Data()),
         nptbins,ptmin,ptmax,asyBins,asyMin,asyMax); 
        fhDeltaIATot[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIATot[imatch]->SetYTitle("#it{A}_{cell in cluster}^{total}");
        outputContainer->Add(fhDeltaIATot[imatch]); 
        
        fhDeltaIAM02[imatch]  = new TH3F 
        (Form("hDeltaIAM02_%s",matchCase[imatch].Data()),
         Form("Cluster *asymmetry* in cell units vs #lambda^{2}_{0} for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,nShShBins,minShSh,maxShSh,asyBins,asyMin,asyMax); 
        fhDeltaIAM02[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIAM02[imatch]->SetYTitle("#lambda^{2}_{0}");
        fhDeltaIAM02[imatch]->SetZTitle("#it{A}_{cell in cluster}");
        outputContainer->Add(fhDeltaIAM02[imatch]); 
        
        fhDeltaIATotM02[imatch]  = new TH3F 
        (Form("hDeltaIATotM02_%s",matchCase[imatch].Data()),
         Form("Cluster *total asymmetry* in cell units vs #lambda^{2}_{0} for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,nShShBins,minShSh,maxShSh,asyBins,asyMin,asyMax); 
        fhDeltaIATotM02[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIATotM02[imatch]->SetYTitle("#lambda^{2}_{0}");
        fhDeltaIATotM02[imatch]->SetZTitle("#it{A}_{cell in cluster}^{total}");
        outputContainer->Add(fhDeltaIATotM02[imatch]); 
        
        fhDeltaIAM20[imatch]  = new TH3F 
        (Form("hDeltaIAM20_%s",matchCase[imatch].Data()),
         Form("Cluster *asymmetry* in cell units vs #lambda^{2}_{1} for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,(Int_t)nShShBins/1.5,minShSh,(Int_t)maxShSh/1.5,asyBins,asyMin,asyMax); 
        fhDeltaIAM20[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIAM20[imatch]->SetYTitle("#lambda^{2}_{1}");
        fhDeltaIAM20[imatch]->SetZTitle("#it{A}_{cell in cluster}");
        outputContainer->Add(fhDeltaIAM20[imatch]); 
        
        fhDeltaIATotM20[imatch]  = new TH3F 
        (Form("hDeltaIATotM20_%s",matchCase[imatch].Data()),
         Form("Cluster *total asymmetry* in cell units vs #lambda^{2}_{1} for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,(Int_t)nShShBins/1.5,minShSh,(Int_t)maxShSh/1.5,asyBins,asyMin,asyMax); 
        fhDeltaIATotM20[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIATotM20[imatch]->SetYTitle("#lambda^{2}_{1}");
        fhDeltaIATotM20[imatch]->SetZTitle("#it{A}_{cell in cluster}^{total}");
        outputContainer->Add(fhDeltaIATotM20[imatch]); 
        
        fhDeltaIANCells[imatch]  = new TH3F 
        (Form("hDeltaIANCells_%s",matchCase[imatch].Data()),
         Form("Cluster *asymmetry* in cell units vs N cells in cluster for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,asyBins,asyMin,asyMax); 
        fhDeltaIANCells[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIANCells[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhDeltaIANCells[imatch]->SetZTitle("#it{A}_{cell in cluster}");
        outputContainer->Add(fhDeltaIANCells[imatch]); 
        
        fhDeltaIATotNCells[imatch]  = new TH3F 
        (Form("hDeltaIATotNCells_%s",matchCase[imatch].Data()),
         Form("Cluster *total asymmetry* in cell units vs N cells in cluster for %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,asyBins,asyMin,asyMax); 
        fhDeltaIATotNCells[imatch]->SetXTitle("#it{E}_{cluster}");
        fhDeltaIATotNCells[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhDeltaIATotNCells[imatch]->SetZTitle("#it{A}_{cell in cluster}^{total}");
        outputContainer->Add(fhDeltaIATotNCells[imatch]); 
        
        if ( IsDataMC() )
        {
          fhDeltaIAOrigin[imatch]  = new TH3F 
          (Form("hDeltaIAOrigin_%s",matchCase[imatch].Data()),
           Form("Cluster *asymmetry* in cell units vs E vs primary ID, for ID %s",matchCase[imatch].Data()),
           nEbins,minE,maxE,10,-0.5,9.5,asyBins,asyMin,asyMax); 
          fhDeltaIAOrigin[imatch]->SetXTitle("#it{E}_{cluster}");
          fhDeltaIAOrigin[imatch]->SetYTitle("particle");
          fhDeltaIAOrigin[imatch]->SetZTitle("#it{A}_{cell in cluster}");
          for(Int_t ilabel = 1; ilabel <=10; ilabel++)
          fhDeltaIAOrigin[imatch]->GetYaxis()->SetBinLabel(ilabel,mcParticleStringLabel[ilabel-1]);
          outputContainer->Add(fhDeltaIAOrigin[imatch]);    
          
          fhDeltaIATotOrigin[imatch]  = new TH3F 
          (Form("hDeltaIATotOrigin_%s",matchCase[imatch].Data()),
           Form("Cluster *total asymmetry* in cell units vs E vs primary ID, for ID %s",matchCase[imatch].Data()),
           nEbins,minE,maxE,10,-0.5,9.5,asyBins,asyMin,asyMax); 
          fhDeltaIATotOrigin[imatch]->SetXTitle("#it{E}_{cluster}");
          fhDeltaIATotOrigin[imatch]->SetYTitle("particle");
          fhDeltaIATotOrigin[imatch]->SetZTitle("#it{A}_{cell in cluster}");
          for(Int_t ilabel = 1; ilabel <=10; ilabel++)
          fhDeltaIATotOrigin[imatch]->GetYaxis()->SetBinLabel(ilabel,mcParticleStringLabel[ilabel-1]);
          outputContainer->Add(fhDeltaIATotOrigin[imatch]);    
        }
      } // Cluster asymmetry
      
      fhClusterMaxCellDiffM02[imatch]  = new TH3F 
      (Form("hClusterMaxCellDiffM02_%s",matchCase[imatch].Data()),
       Form("#it{E}_{cluster} vs (#it{E}_{cluster} - #it{E}_{cell max})/#it{E}_{cluster} vs #lambda_{0}^{2} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE, 20,0,1.,nShShBins,minShSh,maxShSh); 
      fhClusterMaxCellDiffM02[imatch]->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellDiffM02[imatch]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
      fhClusterMaxCellDiffM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellDiffM02[imatch]);  
      
      fhClusterTimeEnergyM02[imatch]  = new TH3F 
      (Form("hClusterTimeEnergyM02_%s",matchCase[imatch].Data()),
       Form("#it{E} vs TOF vs #lambda_{0}^{2} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE,45,-25.5,20.5,nShShBins,minShSh,maxShSh); 
      fhClusterTimeEnergyM02[imatch]->SetXTitle("#it{E} (GeV) ");
      fhClusterTimeEnergyM02[imatch]->SetYTitle("TOF (ns)");
      fhClusterTimeEnergyM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterTimeEnergyM02[imatch]);    
      
      fhNCellsPerClusterM02[imatch]  = new TH3F 
      (Form("hNCellsPerClusterM02_%s",matchCase[imatch].Data()),
       Form("#it{E} vs #it{n}_{cells} vs #lambda_{0}^{2} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE,cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhNCellsPerClusterM02[imatch]->SetXTitle("#it{E} (GeV)");
      fhNCellsPerClusterM02[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhNCellsPerClusterM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsPerClusterM02[imatch]); 

      fhNCellsPerClusterM20[imatch]  = new TH3F 
      (Form("hNCellsPerClusterM20_%s",matchCase[imatch].Data()),
       Form("#it{E} vs #it{n}_{cells} vs #lambda_{1}^{2} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE,cellBins,cellMin,cellMax,(Int_t)nShShBins/1.5,minShSh,(Int_t)maxShSh/1.5); 
      fhNCellsPerClusterM20[imatch]->SetXTitle("#it{E} (GeV)");
      fhNCellsPerClusterM20[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      fhNCellsPerClusterM20[imatch]->SetZTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhNCellsPerClusterM20[imatch]); 
      
      fhSMM02[imatch]  = new TH3F 
      (Form("hSMM02_%s",matchCase[imatch].Data()),
       Form("#it{E} vs SM number vs #lambda_{0}^{2}, #it{n}_{cells}^{w>0.01}>%d, for ID %s",fNCellMinShape,matchCase[imatch].Data()),
       nEbins,minE,maxE,fNModules,-0.5,fNModules-0.5,nShShBins,minShSh,maxShSh); 
      fhSMM02[imatch]->SetXTitle("#it{E} (GeV)");
      fhSMM02[imatch]->SetYTitle("SM number");
      fhSMM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhSMM02[imatch]); 
      
      fhSMM02NoCut[imatch]  = new TH3F 
      (Form("hSMM02NoCut_%s",matchCase[imatch].Data()),
       Form("#it{E} vs SM number vs #lambda_{0}^{2} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE,fNModules,-0.5,fNModules-0.5,nShShBins,minShSh,maxShSh); 
      fhSMM02NoCut[imatch]->SetXTitle("#it{E} (GeV)");
      fhSMM02NoCut[imatch]->SetYTitle("SM number");
      fhSMM02NoCut[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhSMM02NoCut[imatch]); 
      
      fhSMNCell[imatch]  = new TH3F 
      (Form("hSMNCell_%s",matchCase[imatch].Data()),
       Form("#it{E} vs SM number vs  #it{n}_{cells}^{w>0.01} for ID %s",matchCase[imatch].Data()),
       nEbins,minE,maxE,fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax); 
      fhSMNCell[imatch]->SetXTitle("#it{E} (GeV)");
      fhSMNCell[imatch]->SetYTitle("SM number");
      fhSMNCell[imatch]->SetZTitle("#it{n}_{cells}^{w>0.01}");
      outputContainer->Add(fhSMNCell[imatch]); 
      
      fhSMNCellM02[imatch]  = new TH3F 
      (Form("hSMNCellM02_%s",matchCase[imatch].Data()),
       Form("SM number vs #it{n}_{cells}^{w>0.01} vs #lambda_{0}^{2}, "
            "%2.2f<#it{E}<%2.2f GeV, for ID %s",fEMinShape,fEMaxShape,matchCase[imatch].Data()),
       fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
      fhSMNCellM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      fhSMNCellM02[imatch]->SetXTitle("SM number");
      fhSMNCellM02[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
      outputContainer->Add(fhSMNCellM02[imatch]); 
      
      fhColM02[imatch]  = new TH3F 
      (Form("hColM02_%s",matchCase[imatch].Data()),
       Form("#it{E} vs column number vs #lambda_{0}^{2}, #it{n}_{cells}^{w>0.01}>%d, for ID %s",fNCellMinShape,matchCase[imatch].Data()),
       nEbins,minE,maxE,48,-0.5,47.5,nShShBins,minShSh,maxShSh); 
      fhColM02[imatch]->SetXTitle("#it{E} (GeV)");
      fhColM02[imatch]->SetYTitle("column number");
      fhColM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhColM02[imatch]); 

      fhRowM02[imatch]  = new TH3F 
      (Form("hRowM02_%s",matchCase[imatch].Data()),
       Form("#it{E} vs row number vs #lambda_{0}^{2}, #it{n}_{cells}^{w>0.01}>%d, for ID %s",fNCellMinShape,matchCase[imatch].Data()),
       nEbins,minE,maxE,24,-0.5,23.5,nShShBins,minShSh,maxShSh); 
      fhRowM02[imatch]->SetXTitle("#it{E} (GeV)");
      fhRowM02[imatch]->SetYTitle("row number");
      fhRowM02[imatch]->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhRowM02[imatch]); 
      
      if ( fStudyShapeParam && GetCalorimeter() == kEMCAL )
      {
        fhNCellsPerClusterMEta[imatch]  = new TH3F 
        (Form("hNCellsPerClusterMEta_%s",matchCase[imatch].Data()),
         Form("#it{E} vs #it{n}_{cells} vs #sigma_{#eta}^{2} for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
        fhNCellsPerClusterMEta[imatch]->SetXTitle("#it{E} (GeV)");
        fhNCellsPerClusterMEta[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhNCellsPerClusterMEta[imatch]->SetZTitle("#sigma_{#eta}^{2}");
        outputContainer->Add(fhNCellsPerClusterMEta[imatch]); 
        
        fhNCellsPerClusterMPhi[imatch]  = new TH3F 
        (Form("hNCellsPerClusterMPhi_%s",matchCase[imatch].Data()),
         Form("#it{E} vs #it{n}_{cells} vs #sigma_{#varphi}^{2} for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,nShShBins,minShSh,maxShSh); 
        fhNCellsPerClusterMPhi[imatch]->SetXTitle("#it{E} (GeV)");
        fhNCellsPerClusterMPhi[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhNCellsPerClusterMPhi[imatch]->SetZTitle("#sigma_{#varphi}^{2}");
        outputContainer->Add(fhNCellsPerClusterMPhi[imatch]); 
        
        fhNCellsPerClusterMEtaPhi[imatch]  = new TH3F 
        (Form("hNCellsPerClusterMEtaPhi_%s",matchCase[imatch].Data()),
         Form("#it{E} vs #it{n}_{cells} vs #sigma_{#eta#varphi}^{2} for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,nShShBins,-1*maxShSh,maxShSh); 
        fhNCellsPerClusterMEtaPhi[imatch]->SetXTitle("#it{E} (GeV)");
        fhNCellsPerClusterMEtaPhi[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhNCellsPerClusterMEtaPhi[imatch]->SetZTitle("#sigma_{#eta#varphi}^{2}");
        outputContainer->Add(fhNCellsPerClusterMEtaPhi[imatch]); 
        
        fhNCellsPerClusterMEtaPhiA[imatch]  = new TH3F 
        (Form("hNCellsPerClusterMEtaPhiA_%s",matchCase[imatch].Data()),
         Form("#it{E} vs #it{n}_{cells} vs (#sigma_{#varphi}^{2}-#sigma_{#eta}^{2})/(#sigma_{#varphi}^{2}+#sigma_{#eta}^{2}) for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,cellBins,cellMin,cellMax,nShShBins,-1*maxShSh,maxShSh); 
        fhNCellsPerClusterMEtaPhiA[imatch]->SetXTitle("#it{E} (GeV)");
        fhNCellsPerClusterMEtaPhiA[imatch]->SetYTitle("#it{n}_{cells}^{w>0.01}");
        fhNCellsPerClusterMEtaPhiA[imatch]->SetZTitle("(#sigma_{#varphi}^{2}-#sigma_{#eta}^{2})/(#sigma_{#varphi}^{2}+#sigma_{#eta}^{2})");
        outputContainer->Add(fhNCellsPerClusterMEtaPhiA[imatch]); 
      }
      
      if ( IsDataMC() )
      {
        fhOriginE[imatch]  = new TH2F 
        (Form("hOrigin_%s",matchCase[imatch].Data()),
         Form("#it{E} vs origin for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE, 10,-0.5,9.5); 
        fhOriginE[imatch]->SetXTitle("#it{E} (GeV)");
        fhOriginE[imatch]->SetYTitle("particle");
        for(Int_t ilabel = 1; ilabel <=10; ilabel++)
          fhOriginE[imatch]->GetYaxis()->SetBinLabel(ilabel ,mcParticleStringLabel[ilabel-1]);
        outputContainer->Add(fhOriginE[imatch]); 
        
        fhOriginM02[imatch]  = new TH3F 
        (Form("hOriginM02_%s",matchCase[imatch].Data()),
         Form("#it{E} vs origin vs #lambda_{0}^{2} for ID %s",matchCase[imatch].Data()),
         nEbins,minE,maxE,10,-0.5,9.5,nShShBins,minShSh,maxShSh); 
        fhOriginM02[imatch]->SetXTitle("#it{E} (GeV)");
        fhOriginM02[imatch]->SetYTitle("particle");
        fhOriginM02[imatch]->SetZTitle("#lambda_{0}^{2}");
        for(Int_t ilabel = 1; ilabel <=10; ilabel++)
          fhOriginM02[imatch]->GetYaxis()->SetBinLabel(ilabel,mcParticleStringLabel[ilabel-1]);
        outputContainer->Add(fhOriginM02[imatch]); 
      } // MC
    } // match loop
    
    //    fhCellTimeSpreadRespectToCellMaxM02 = new TH3F 
    //    ("hCellTimeSpreadRespectToCellMaxM02",
    //     "#it{E} vs t_{cell max}-t_{cell i} vs #lambda_{0}^{2}", 
    //     nEbins,minE,maxE,100,-100,100,nShShBins,minShSh,maxShSh); 
    //    fhCellTimeSpreadRespectToCellMaxM02->SetXTitle("#it{E} (GeV)");
    //    fhCellTimeSpreadRespectToCellMaxM02->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
    //    fhCellTimeSpreadRespectToCellMaxM02->SetZTitle("#lambda_{0}^{2}");
    //    outputContainer->Add(fhCellTimeSpreadRespectToCellMaxM02);
    
    fhClusterMaxCellCloseCellRatioM02  = new TH3F 
    ("hClusterMaxCellCloseCellRatioM02","#it{E} vs #it{E}_{cell-i}/#it{E}_{cell max} vs #lambda_{0}^{2}",
     nEbins,minE,maxE, 20,0,1.,nShShBins,minShSh,maxShSh); 
    fhClusterMaxCellCloseCellRatioM02->SetXTitle("#it{E}_{cluster} (GeV) ");
    fhClusterMaxCellCloseCellRatioM02->SetYTitle("#it{E}_{cell i}/#it{E}_{cell max}");
    fhClusterMaxCellCloseCellRatioM02->SetZTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhClusterMaxCellCloseCellRatioM02);
    
    //    fhClusterMaxCellCloseCellDiffM02  = new TH3F 
    //    ("hClusterMaxCellCloseCellDiffM02",
    //     "#it{E} vs #it{E}_{cell max}-#it{E}_{cell i} vs #lambda_{0}^{2}",
    //     nEbins,minE,maxE, 40,0,20,nShShBins,minShSh,maxShSh); 
    //    fhClusterMaxCellCloseCellDiffM02->SetXTitle("#it{E}_{cluster} (GeV) ");
    //    fhClusterMaxCellCloseCellDiffM02->SetYTitle("#it{E}_{cell max}-#it{E}_{cell i} (GeV)");
    //    fhClusterMaxCellCloseCellDiffM02->SetZTitle("#lambda_{0}^{2}");
    //    outputContainer->Add(fhClusterMaxCellCloseCellDiffM02);   
    
    if(fStudyExotic)
    {
      fhClusterMaxCellECrossM02  = new TH3F 
      ("hClusterMaxCellECrossM02",
       "#it{E} vs exoticity vs #lambda_{0}^{2}",
       nEbins,minE,maxE, 40,0.6,1.,nShShBins,minShSh,maxShSh); 
      fhClusterMaxCellECrossM02->SetXTitle("#it{E}_{cluster} (GeV) ");
      fhClusterMaxCellECrossM02->SetYTitle("1- #it{E}_{cross}/#it{E}_{cell max}");
      fhClusterMaxCellECrossM02->SetZTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhClusterMaxCellECrossM02);
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

  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}

//______________________________
/// Check if the calorimeter setting is ok, if not abort.
//______________________________
void AliAnaClusterShapeCorrelStudies::Init()
{
  if(GetCalorimeter() != kPHOS && GetCalorimeter() !=kEMCAL)
    AliFatal(Form("Wrong calorimeter name <%s>", GetCalorimeterString().Data()));
  
  //if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
  //  AliFatal("Analysis of reconstructed data, MC reader not aplicable");
}

//________________________________________
/// Initialize the parameters of the analysis.
//________________________________________
void AliAnaClusterShapeCorrelStudies::InitParameters()
{ 
  AddToHistogramsName("AnaClusShapeCorr_");
   
  fM02Min = 0.05; 
  
  fNCellMin = 3; // at least 3
 
  fMinDistToBad = 5;
  
  fInvMassMinECut  = 0.5; // 500 MeV
  fInvMassMaxECut  = 12; 

  fInvMassMinM02Cut = 0.1; 
  fInvMassMaxM02Cut = 0.4; 
  
  fInvMassMaxTimeDifference = 70; // ns
  
  fInvMassMaxOpenAngle = 100*TMath::DegToRad(); // 100 degrees
  
  fNEBinCuts = 7;
  fEBinCuts[0] = 2. ; fEBinCuts[1] = 4. ;  fEBinCuts[2] = 6. ;
  fEBinCuts[3] = 8. ; fEBinCuts[4] = 10.;  fEBinCuts[5] = 12.;
  fEBinCuts[6] = 16.; fEBinCuts[7] = 20.;  
  for(Int_t i = fNEBinCuts+1; i < 15; i++) fEBinCuts[i] = 1000.;
  
  fEMinShape     =  8;
  fEMaxShape     = 12;
  fNCellMinShape =  4;
}

//________________________________________
/// Initialize the dedx range depending on run number or MC
//________________________________________
void AliAnaClusterShapeCorrelStudies::InitdEdXParameters()
{
  if(IsDataMC())
  {
    fdEdXMinEle = 80;
    fdEdXMaxEle =100;
    fdEdXMinHad = 45;
    fdEdXMaxHad = 78;
    
    AliInfo(Form("dEdX cuts init for MC: %d<dedx_ele<%d; %d<dedx_had<%d",
                 fdEdXMinEle,fdEdXMaxEle,fdEdXMinHad,fdEdXMaxHad));
  } // MC
  else
  {
    Int_t  runNumber = GetReader()->GetInputEvent()->GetRunNumber();

    // data
    if      ( runNumber < 146861 ) // LHC11a, LHC10 check for both
    { 
      fdEdXMinEle = 72;
      fdEdXMaxEle = 90;
      fdEdXMinHad = 40;
      fdEdXMaxHad = 70;
    }
    else if ( runNumber < 156000 ) // LHC11c, LHC11b, check for b
    {
      fdEdXMinEle = 60;
      fdEdXMaxEle = 72;
      fdEdXMinHad = 35;
      fdEdXMaxHad = 56;
    }
    else if ( runNumber < 165000 ) // LHC11d and other, check for other
    {
      fdEdXMinEle = 74;
      fdEdXMaxEle = 90;
      fdEdXMinHad = 40;
      fdEdXMaxHad = 70;
    }
    else    // LHC12cdfhi and beyond, check for >=LHC13!
    {
      fdEdXMinEle = 78;
      fdEdXMaxEle = 95;
      fdEdXMinHad = 40;
      fdEdXMaxHad = 74;
    }
    
    AliInfo(Form("dEdX cuts init for run %d: %d<dedx_ele<%d; %d<dedx_had<%d",
                 runNumber,fdEdXMinEle,fdEdXMaxEle,fdEdXMinHad,fdEdXMaxHad));
  } // data
}

//_____________________________________________________________________________
/// Identify cluster as exotic or not.
/// \return true if good
///
/// \param absIdMax: absolute ID of main cell in clusrter 
/// \param m02: shower shape main axis of cluster
/// \param nCellsPerCluster: number of cells in cluster
//_____________________________________________________________________________
Bool_t AliAnaClusterShapeCorrelStudies::IsGoodCluster(Int_t absIdMax, Float_t m02, 
                                                      Int_t nCellsPerCluster)
{    
  if(  m02 < fM02Min || nCellsPerCluster < fNCellMin ) return kFALSE ; // mild shower shape cut for exotics
  
  Int_t bc = (GetReader()->GetInputEvent())->GetBunchCrossNumber();

  if(GetCalorimeter() == kEMCAL) 
  {
    if(!GetCaloUtils()->GetEMCALRecoUtils()->IsRejectExoticCluster())
      return  !(GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(absIdMax,fCaloCellList,bc));
    else
      return kTRUE;
  }
  else // PHOS
  {
    Float_t ampMax = fCaloCellList->GetCellAmplitude(absIdMax);
    GetCaloUtils()->RecalibrateCellAmplitude(ampMax, GetCalorimeter(), absIdMax);
    
    if(ampMax < 0.01) return kFALSE;
    if(1-GetCaloUtils()->GetECross(absIdMax,fCaloCellList,bc)/ampMax > 0.95) 
      return kFALSE;
    else                                          
      return kTRUE;
  }
}

//_________________________________________________________
/// Print some relevant parameters set for the analysis.
//_________________________________________________________
void AliAnaClusterShapeCorrelStudies::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",GetCalorimeterString().Data());
  printf("Min n cells    : %d\n"   , fNCellMin) ;
  printf("Min dist to bad: %2.1f\n", fMinDistToBad) ;
  printf("Min M02        : %1.2f\n", fM02Min) ;
  
  printf("Inv. Mass %2.1f < E_clus < %2.1f GeV/c\n"  , fInvMassMinECut  , fInvMassMaxECut  ) ;
  printf("Inv. Mass %2.1f < M02_clus < %2.1f GeV/c\n", fInvMassMinM02Cut, fInvMassMaxM02Cut) ;
  printf("Inv. Mass open angle : %2.1f deg\n"        , fInvMassMaxOpenAngle*TMath::RadToDeg()) ;
  printf("Inv. Mass time difference: %2.1f ns\n"     , fInvMassMaxTimeDifference) ;
}

//_____________________________________________________
/// Main task method.
/// It just recovers the list of clusters or cells depending
/// on the calorimeter and passes it to the method doing the
/// loop on the clusters and doing the different analysis.
//_____________________________________________________
void  AliAnaClusterShapeCorrelStudies::MakeAnalysisFillHistograms()
{
  AliDebug(1,"Start");

  //Print("");
  
  // Get List with CaloClusters , calo Cells, init min amplitude
  if      (GetCalorimeter() == kPHOS)
  {
    fCaloClusList = GetPHOSClusters();
    fCaloCellList = GetPHOSCells();
  }
  else if (GetCalorimeter() == kEMCAL)
  {
    fCaloClusList = GetEMCALClusters();
    fCaloCellList = GetEMCALCells();
  }
  else
    AliFatal(Form("AliAnaClusterShapeCorrelStudies::MakeAnalysisFillHistograms() - Wrong calorimeter name <%s>, END", GetCalorimeterString().Data()));
  
  if( !fCaloClusList || !fCaloCellList )
  {
    AliFatal(Form("AliAnaClusterShapeCorrelStudies::MakeAnalysisFillHistograms() - No CaloClusters or CaloCells available"));
    return; // trick coverity
  }
  
  if(fCaloClusList->GetEntriesFast() == 0) return ;
  
  AliDebug(1,Form("N cells %d, N clusters %d \n",fCaloCellList->GetNumberOfCells(),fCaloClusList->GetEntriesFast()));
  
  // Clusters
  ClusterLoopHistograms();
    
  AliDebug(1,"End");
}

//_________________________________________________________________________________
/// Check cluster weights, check the effect of different w0 parameter on shower shape
/// Check effect of time and energy cuts at cell level on the shower shape
///
/// \param clus: cluster pointer
/// \param mcTag: mc origin map
//_________________________________________________________________________________
void AliAnaClusterShapeCorrelStudies::WeightHistograms(AliVCluster *clus, Int_t mcTag)
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
    Float_t amp = fCaloCellList->GetCellAmplitude(id);
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
    Float_t amp = fCaloCellList->GetCellAmplitude(id);
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
    
    Int_t mcIndex = -1;
    if(IsDataMC() && clus->GetNLabels() > 0)
    {      
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPhoton)   &&
         !GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPi0)      &&
         !GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCEta)      &&
         !GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCConversion)        ){
        mcIndex = 0;
      } // Pure Photon
      else if( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCElectron) &&
              !GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCConversion)        ){
        mcIndex = 1;
      } // Electron
      else if( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCConversion)        ){
        mcIndex = 2;
      } // Conversion
      else if( GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPi0) ){
        mcIndex = 3;
      }// Pi0
      else if(!GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCEta) &&
              !GetMCAnalysisUtils()->CheckTagBit(mcTag, AliMCAnalysisUtils::kMCPhoton) ){
        mcIndex = 4;
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
      
      if(IsDataMC() && mcIndex >= 0)
      {
        fhLambda0ForW0MC[iw][mcIndex]->Fill(energy, clus->GetM02(), GetEventWeight());
//      fhLambda1ForW0MC[iw][mcIndex]->Fill(energy, clus->GetM20(), GetEventWeight());
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



