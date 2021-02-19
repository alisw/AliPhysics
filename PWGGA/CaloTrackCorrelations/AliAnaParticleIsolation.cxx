/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TClass.h>
#include <TH2F.h>
#include "TDatabasePDG.h"
#include <TCustomBinning.h>

// --- AliRoot/Analysis system ---
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"

// --- CaloTrackCorr classes ---
#include "AliAnaParticleIsolation.h"
#include "AliCaloTrackReader.h"
#include "AliIsolationCut.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliCaloTrackParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"

// --- Detectors ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

/// \cond CLASSIMP
ClassImp(AliAnaParticleIsolation) ;
/// \endcond

//__________________________________________________
/// Default constructor. Initialize parameters.
//__________________________________________________
AliAnaParticleIsolation::AliAnaParticleIsolation() :
AliAnaCaloTrackCorrBaseClass(),
fIsoDetector(-1),                 fIsoDetectorString(""),
fFillTMHisto(0),                  fFillSSHisto(1),      
fFillPerSMHistograms(0),          fFillPerTCardIndexHistograms(0),         fTCardIndex(-1),                
fFillEMCALRegionHistograms(0),   
fFillOverlapHistograms(0),                        
fStudyTracksInCone(0),            fStudyMCConversionRadius(0),             fFillTrackOriginHistograms(0),
fFillTaggedDecayHistograms(0),    fNDecayBits(0),
fDecayBits(),                     fFillNLMHistograms(0),                   
fFillOnlyTH3Histo(0),             fFillIsolationControlHistograms(0),
fLeadingOnly(0),                  fCheckLeadingWithNeutralClusters(0),
fSelectPrimariesInCone(0),        fMakePrimaryPi0DecayStudy(0),
fMinCellsAngleOverlap(0),         fNumberMCParticleCases(fgkNmcTypes),

fStudyPtCutInCone(0),             fNPtCutsInCone(0),                        
fMinPtCutInCone(),                fMaxPtCutInCone(),
fStudyEtaCutInCone(0),            fNEtaCutsInCone(0),                       fEtaCutInCone(),
fStudyRCutInCone(0),              fNRCutsInCone(0),                         fRCutInCone(),
fStudyNCellsCut(0),               fNNCellsInCandidate(0),                   fNCellsInCandidate(),      
fNCellsWithWeight(0),             fTrigSupMod(-1),
fStudyExoticTrigger(0),           fNExoCutInCandidate(0),                   fExoCutInCandidate(),
fMomentum(),                      fMomIso(),
fMomDaugh1(),                     fMomDaugh2(),
fTrackVector(),                   fProdVertex(),
fCluster(0),                      fClustersArr(0),                          fCaloCells(0),                
fIsExoticTrigger(0),              fClusterExoticity(1),
// Histograms

fhPtInConeExoTrigger(0),          fhPtClusterInConeExoTrigger(0),           fhPtTrackInConeExoTrigger(0),
fhPtTrackInConeOtherBCPileUpSPD(0), fhPtTrackInConeVtxBC0(0),
fhPtTrackInConeBC0PileUpSPD(0),
fhPerpConeSumPtTOFBC0(0),         fhPtInPerpConeTOFBC0(0),
fhEtaPhiInPerpConeTOFBC0(0),
fhPtM02SumPtCone(0),             fhPtM02SumPtConeCharged(0),
fhPtM02SumPtConeCent(0x0),       fhPtM02SumPtConeChargedCent(0x0),          
fhPtM02SumPtConeCentMC(0x0),     fhPtM02SumPtConeChargedCentMC(0x0),          
fhConeSumPtExoTrigger(0),        fhConeSumPtClusterExoTrigger(0),            fhConeSumPtTrackExoTrigger(0),                      

fhPtPrimMCPi0DecayPairOutOfCone(0),
fhPtPrimMCPi0DecayPairOutOfAcceptance(0),
fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap(0),
fhPtPrimMCPi0DecayPairAcceptInConeLowPt(0),
fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap(0),
fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE(0),
fhPtPrimMCPi0DecayPairNoOverlap(0),
fhPtPrimMCPi0DecayIsoPairOutOfCone(0),
fhPtPrimMCPi0DecayIsoPairOutOfAcceptance(0),
fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap(0),
fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt(0),
fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap(0),
fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE(0),
fhPtPrimMCPi0DecayIsoPairNoOverlap(0),
fhPtPrimMCPi0Overlap(0),                    fhPtPrimMCPi0IsoOverlap(0),
fhPtPrimMCEtaDecayPairOutOfCone(0),
fhPtPrimMCEtaDecayPairOutOfAcceptance(0),
fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap(0),
fhPtPrimMCEtaDecayPairAcceptInConeLowPt(0),
fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap(0),
fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE(0),
fhPtPrimMCEtaDecayPairNoOverlap(0),
fhPtPrimMCEtaDecayIsoPairOutOfCone(0),
fhPtPrimMCEtaDecayIsoPairOutOfAcceptance(0),
fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap(0),
fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt(0),
fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap(0),
fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE(0),
fhPtPrimMCEtaDecayIsoPairNoOverlap(0),
fhPtPrimMCEtaOverlap(0),                    fhPtPrimMCEtaIsoOverlap(0), 

// PileUp
fhTimeENoCut(0),                  fhTimeESPD(0),                  fhTimeESPDMulti(0),
fhTimeNPileUpVertSPD(0),          fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0),  fhTimePileUpMainVertexZDiamond(0),

fhPtClusterInConePerRCut(0),         fhPtTrackInConePerRCut(0),           
fhConeSumPtClusterPerRCut(0),        fhConeSumPtTrackPerRCut(0),

fhConeNClusterPerMinPtCut(0),        
fhEtaBandConeNClusterPerMinPtCut(0),   fhConeNClusterSubEtaBandPerMinPtCut(0), 
fhPhiBandConeNClusterPerMinPtCut(0),   fhConeNClusterSubPhiBandPerMinPtCut(0), 

fhConeNTrackPerMinPtCut(0),         
fhPerpConeNTrackPerMinPtCut(0),      fhConeNTrackSubPerpPerMinPtCut(0), 
fhEtaBandConeNTrackPerMinPtCut(0),   fhConeNTrackSubEtaBandPerMinPtCut(0), 
fhPhiBandConeNTrackPerMinPtCut(0),   fhConeNTrackSubPhiBandPerMinPtCut(0), 

fhConeNClusterPerMinPtCutCent(0),    
fhEtaBandConeNClusterPerMinPtCutCent(0),fhConeNClusterSubEtaBandPerMinPtCutCent(0),     
fhPhiBandConeNClusterPerMinPtCutCent(0),fhConeNClusterSubPhiBandPerMinPtCutCent(0),   

fhConeNTrackPerMinPtCutCent(0),     
fhPerpConeNTrackPerMinPtCutCent(0),  fhConeNTrackSubPerpPerMinPtCutCent(0),     
fhEtaBandConeNTrackPerMinPtCutCent(0),fhConeNTrackSubEtaBandPerMinPtCutCent(0),     
fhPhiBandConeNTrackPerMinPtCutCent(0),fhConeNTrackSubPhiBandPerMinPtCutCent(0),   

fhConeSumPtClusterPerMinPtCut(0),    
fhEtaBandConeSumPtClusterPerMinPtCut(0),fhConeSumPtClusterSubEtaBandPerMinPtCut(0),
fhPhiBandConeSumPtClusterPerMinPtCut(0),fhConeSumPtClusterSubPhiBandPerMinPtCut(0),

fhConeSumPtTrackPerMinPtCut(0),     
fhPerpConeSumPtTrackPerMinPtCut(0),  fhConeSumPtTrackSubPerpPerMinPtCut(0),
fhEtaBandConeSumPtTrackPerMinPtCut(0),fhConeSumPtTrackSubEtaBandPerMinPtCut(0),
fhPhiBandConeSumPtTrackPerMinPtCut(0),fhConeSumPtTrackSubPhiBandPerMinPtCut(0),

fhConeSumPtClusterPerMinPtCutCent(0),
fhEtaBandConeSumPtClusterPerMinPtCutCent(0),fhConeSumPtClusterSubEtaBandPerMinPtCutCent(0), 
fhPhiBandConeSumPtClusterPerMinPtCutCent(0),fhConeSumPtClusterSubPhiBandPerMinPtCutCent(0), 

fhConeSumPtTrackPerMinPtCutCent(0), 
fhPerpConeSumPtTrackPerMinPtCutCent(0), fhConeSumPtTrackSubPerpPerMinPtCutCent(0), 
fhEtaBandConeSumPtTrackPerMinPtCutCent(0),fhConeSumPtTrackSubEtaBandPerMinPtCutCent(0), 
fhPhiBandConeSumPtTrackPerMinPtCutCent(0),fhConeSumPtTrackSubPhiBandPerMinPtCutCent(0), 

fhConeNSubEtaBandPerMinPtCut(0),         fhConeNSubPhiBandPerMinPtCut(0),
fhConeSumPtSubEtaBandPerMinPtCut(0),     fhConeSumPtSubPhiBandPerMinPtCut(0),
fhConeNSubEtaBandPerMinPtCutCent(0),     fhConeNSubPhiBandPerMinPtCutCent(0), 
fhConeSumPtSubEtaBandPerMinPtCutCent(0), fhConeSumPtSubPhiBandPerMinPtCutCent(0),

fhConeSumPtClusterPerMaxPtCut(0),    fhConeSumPtTrackPerMaxPtCut(0),     
fhConeSumPtTrackPerEtaCut(0),        

fhPtClusterInConePerNCellCut(0),     fhPtTrackInConePerNCellCut(0),      
fhConeSumPtClusterPerNCellCut(0),    fhConeSumPtTrackPerNCellCut(0),       
fhPtClusterInConePerExoCut(0),       fhPtTrackInConePerExoCut(0),         
fhConeSumPtClusterPerExoCut(0),      fhConeSumPtTrackPerExoCut(0),              

fhConeSumPtTrackTOFBC0(0), fhConeSumPtTrackTOFBCN(0), fhConeSumPtTrackTOFNo(0),
fhPtTrackInConeTOFBC0 (0), fhPtTrackInConeTOFBCN (0), fhPtTrackInConeTOFNo (0),
fhPhiTrackInCone(0), fhEtaTrackInCone(0), fhEtaPhiTrackInCone(0),
fhPhiTrackInConeTOFBC0(0), fhPhiTrackInConeTOFBCN(0), fhPhiTrackInConeTOFNo(0),
fhEtaTrackInConeTOFBC0(0), fhEtaTrackInConeTOFBCN(0), fhEtaTrackInConeTOFNo(0),
fhEtaPhiTrackInConeTOFBC0(0), fhEtaPhiTrackInConeTOFBCN(0), fhEtaPhiTrackInConeTOFNo(0),
fhTrackTOFInCone(0),       fhTrackTOFInConeBC0(0),    fhTrackTOFInConeExoTrigger(0),

fhConeSumPtTrackITSRefitOnSPDOn(0),   fhConeSumPtTrackITSRefitOnSPDOff(0),    fhConeSumPtTrackITSRefitOffSPDOff(0),
fhPtTrackInConeITSRefitOnSPDOn(0),    fhPtTrackInConeITSRefitOnSPDOff(0) ,    fhPtTrackInConeITSRefitOffSPDOff(0),
fhPhiTrackInConeITSRefitOnSPDOn(0),   fhPhiTrackInConeITSRefitOnSPDOff(0),    fhPhiTrackInConeITSRefitOffSPDOff(0),
fhEtaTrackInConeITSRefitOnSPDOn(0),   fhEtaTrackInConeITSRefitOnSPDOff(0),    fhEtaTrackInConeITSRefitOffSPDOff(0),
fhEtaPhiTrackInConeITSRefitOnSPDOn(0),fhEtaPhiTrackInConeITSRefitOnSPDOff(0), fhEtaPhiTrackInConeITSRefitOffSPDOff(0),
fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn(0), fhPtTrackInConeTOFBC0ITSRefitOnSPDOn(0), 
fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn(0), fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn(0),fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn(0),
fhPerpConeSumPtITSRefitOnSPDOn (0),       fhPtInPerpConeITSRefitOnSPDOn(0),        fhEtaPhiInPerpConeITSRefitOnSPDOn(0),
fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn (0), fhPtInPerpConeTOFBC0ITSRefitOnSPDOn (0), fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn(0)
{
  InitParameters();
  
  for(Int_t ibit =0; ibit < AliNeutralMesonSelection::fgkMaxNDecayBits; ibit++)
  {
    for(Int_t iso =0; iso < 2; iso++)
    {
      fhPtDecay       [iso][ibit] = 0;
      fhEtaPhiDecay   [iso][ibit] = 0;
      fhPtLambda0Decay[iso][ibit] = 0;
      for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
        fhPtDecayMC[iso][ibit][imc]    = 0;
    }
  }
  
  for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
  {
    fhPtM02SumPtConeMC[imc] = 0;
    fhPtM02SumPtConeChargedMC[imc] = 0;
   
    for(Int_t i = 0; i < 2 ; i++)
    {
      for(Int_t ishsh = 0; ishsh < 2 ; ishsh++)
      {
        fhPtMC[imc][i][ishsh] = 0;
      }
      
      fhConeSumPtM02CutMC          [imc][i] = 0;
      fhPtEtaPhiMC                 [imc][i] = 0;
      fhPtLambda0MC                [imc][i] = 0;
      fhPtLambda0MCConv            [imc][i] = 0;      
      fhPtLambda0MCWith1Overlap    [imc][i] = 0;
      fhPtLambda0MCConvWith1Overlap[imc][i] = 0;      
      fhPtLambda0MCWithNoOverlap    [imc][i] = 0;
      fhPtLambda0MCConvWithNoOverlap[imc][i] = 0;
      fhPtNOverlap                 [imc][i] = 0;
      fhPtNOverlapConv             [imc][i] = 0;
      fhTrackMatchedDEtaMC         [imc][i] = 0 ;
      fhTrackMatchedDPhiMC         [imc][i] = 0 ;
      fhTrackMatchedDEtaDPhiMC     [imc][i] = 0 ;
    }
  }

  for(Int_t imc = 0; imc < 4; imc++)
  {
    fhPtTrackInConeMCPrimary       [imc] = 0;
    fhPtTrackInConeMCSecondary     [imc] = 0;  
    fhPtTrackInConeMCPrimaryGener  [imc] = 0;
    fhPtTrackInConeMCSecondaryGener[imc] = 0;  
  }
  
  for(Int_t i = 0; i < 2 ; i++)
  {
    fhConeSumPtM02Cut[i] = 0;
    fhConeSumPtCentM02Cut[i] = 0;
    for(Int_t ishsh = 0; ishsh < 2 ; ishsh++)
    {
      fhPt          [i][ishsh] = 0 ; 
      fhPtEtaPhi    [i][ishsh] = 0 ;
      fhPtCentrality[i][ishsh] = 0 ;                 
      fhPtEventPlane[i][ishsh] = 0 ;      
      fhPtNLocMax   [i][ishsh] = 0 ;
    }
    
    fhPtExoTrigger[i] = 0 ;
    
    fhTrackMatchedDEta[i] = 0 ;             fhTrackMatchedDPhi[i] = 0 ;   fhTrackMatchedDEtaDPhi  [i] = 0 ;
    fhdEdx            [i] = 0 ;             fhEOverP          [i] = 0 ;   fhTrackMatchedMCParticle[i] = 0 ;
    fhPtLambda0       [i] = 0 ;             fhPtLambda0TRD    [i] = 0 ;   fhPtLambda0Cent         [i] = 0 ;   
    
    // Number of local maxima in cluster
    fhPtLambda0LocMax1[i] = 0 ;              fhPtLambda1LocMax1[i] = 0 ;
    fhPtLambda0LocMax2[i] = 0 ;              fhPtLambda1LocMax2[i] = 0 ;
    fhPtLambda0LocMaxN[i] = 0 ;              fhPtLambda1LocMaxN[i] = 0 ;
    
    fhMCConversionVertex[i] = fhMCConversionVertexTRD[i] = 0 ;
    for(Int_t iR = 0; iR < 6; iR++) fhMCConversionLambda0RcutTRD[iR][i] = 0;            
    
    for(Int_t ieta = 0; ieta < 4; ieta++) 
    {  
      for(Int_t iphi = 0; iphi < 3; iphi++) 
      {
//        fhLam0EMCALRegion   [i][ieta][iphi] = 0;
//        fhLam0EMCALRegionTRD[i][ieta][iphi] = 0;
//        
//        for(Int_t ireg = 0; ireg < 6; ireg++) 
//        {
//          fhLam0EMCALRegionMCConvRcut   [i][ieta][iphi][ireg] = 0;
//          fhLam0EMCALRegionTRDMCConvRcut[i][ieta][iphi][ireg] = 0;
//        }
        
        for(Int_t ism =0; ism < 20; ism++)
        {
          fhLam0EMCALRegionPerSM            [i][ieta][iphi][ism] = 0; 
          fhConeSumPtTrackEMCALRegionPerSM     [ieta][iphi][ism] = 0; 
          fhConeSumPtClusterEMCALRegionPerSM   [ieta][iphi][ism] = 0; 
        }
      }
    }
    for(Int_t j = 0; j < 7; j++) fhEtaPhiLam0BinPtBin[i][j] = 0 ;
  }
  
  // Acceptance
  for(Int_t i = 0; i < fgkNmcPrimTypes; i++)
  {
    fhConeSumPtPrimMC          [i] = 0;
    fhConeSumPtChargedPrimMC   [i] = 0;
    fhConeSumPtCenPrimMC       [i] = 0;
    fhConeSumPtCenChargedPrimMC[i] = 0;

    fhConeSumPtUESubPrimMC          [i] = 0;
    fhConeSumPtUESubChargedPrimMC   [i] = 0;
    fhConeSumPtUESubCenPrimMC       [i] = 0;
    fhConeSumPtUESubCenChargedPrimMC[i] = 0;

    fhConeSumPtPrimMCEmb          [i] = 0;
    fhConeSumPtChargedPrimMCEmb   [i] = 0;
    fhConeSumPtCenPrimMCEmb       [i] = 0;
    fhConeSumPtCenChargedPrimMCEmb[i] = 0;

    fhConeSumPtUESubPrimMCEmb          [i] = 0;
    fhConeSumPtUESubChargedPrimMCEmb   [i] = 0;
    fhConeSumPtUESubCenPrimMCEmb       [i] = 0;
    fhConeSumPtUESubCenChargedPrimMCEmb[i] = 0;

    fhConeSumPtNeutralChargedRatioPrimMC        [i] = 0;
    fhConeSumPtPerpConeNeutralChargedRatioPrimMC[i] = 0;
    
    fhPtPrimMCiso[i] = 0;
    fhPtPrimMC   [i] = 0;
    fhEtaPrimMC  [i] = 0;
    fhPhiPrimMC  [i] = 0;
  }
  
  // Pile-Up
  
  for(Int_t i = 0 ; i < 7 ; i++)
  {
    fhPtInConePileUp[i] = 0 ;
    fhPtPileUp   [i][0] = 0 ;
    fhPtPileUp   [i][1] = 0 ;
  }
  
  for(Int_t i = 0 ; i < 3 ; i++)
  {
    fhPtTrackInConeDCA    [i] = 0 ;
    fhPtTrackInPerpConeDCA[i] = 0 ;
  }
  
   for(Int_t i = 0 ; i < 4 ; i++)
   {
     fhPtClusterInConePerNCellPerSM [i]=0;
     fhPtTrackInConePerNCellPerSM   [i]=0;
     fhConeSumPtClusterPerNCellPerSM[i]=0;
     fhConeSumPtTrackPerNCellPerSM  [i]=0;
   }
  
  fhPtPerSM[0] = 0;
  fhPtPerSM[1] = 0;
  for(Int_t ism =0; ism < 20; ism++)
  {
    for(Int_t iso =0; iso < 2; iso++)
    {
      fhPtLambda0PerSM     [iso][ism] = 0;
      fhPtLambda0PerSMNCellCut[iso][ism] = 0;
      fhPtNCellPerSM       [iso][ism] = 0;
      fhPtNCellLowM02PerSM [iso][ism] = 0; 
      fhPtNCellHighM02PerSM[iso][ism] = 0;  
    }
    
    fhConeSumPtPerSM       [ism] = 0;
    fhConeSumPtClusterPerSM[ism] = 0;
    fhConeSumPtTrackPerSM  [ism] = 0;
    
    fhPtInConePerSM        [ism] = 0;
    fhPtClusterInConePerSM [ism] = 0;
    fhPtTrackInConePerSM   [ism] = 0;
  }
  
  fhPtPerTCardIndex[0] = 0;
  fhPtPerTCardIndex[1] = 0;
  for(Int_t itc =0; itc < 16; itc++)
  {
    fhPtLambda0PerTCardIndex    [0][itc] = 0;
    fhPtLambda0PerTCardIndex    [1][itc] = 0;
    
    fhConeSumPtPerTCardIndex       [itc] = 0;
    fhConeSumPtClusterPerTCardIndex[itc] = 0;
    fhConeSumPtTrackPerTCardIndex  [itc] = 0;
    
    fhPtInConePerTCardIndex        [itc] = 0;
    fhPtClusterInConePerTCardIndex [itc] = 0;
    fhPtTrackInConePerTCardIndex   [itc] = 0;
  }
  
  for(Int_t icut = 0; icut < 20; icut++) 
  {
    fConeptsumTrackPerMinCut       [icut] = 0;
    fConeNTrackPerMinCut           [icut] = 0;
    
    fConeptsumPerpTrackPerMinCut   [icut] = 0;
    fConeNPerpTrackPerMinCut       [icut] = 0;
    fConeptsumTrackSubPerpPerMinCut[icut] = 0;
    fConeNTrackSubPerpPerMinCut    [icut] = 0;
    
    fConeptsumEtaBandTrackPerMinCut   [icut] = 0;
    fConeNEtaBandTrackPerMinCut       [icut] = 0;
    fConeptsumTrackSubEtaBandPerMinCut[icut] = 0;
    fConeNTrackSubEtaBandPerMinCut    [icut] = 0;
    
    fConeptsumPhiBandTrackPerMinCut   [icut] = 0;
    fConeNPhiBandTrackPerMinCut       [icut] = 0;
    fConeptsumTrackSubPhiBandPerMinCut[icut] = 0;
    fConeNTrackSubPhiBandPerMinCut    [icut] = 0;
    
    fConeptsumClusterPerMinCut       [icut] = 0;
    fConeNClusterPerMinCut           [icut] = 0;
    
    fConeptsumEtaBandClusterPerMinCut   [icut] = 0;
    fConeNEtaBandClusterPerMinCut       [icut] = 0;
    fConeptsumClusterSubEtaBandPerMinCut[icut] = 0;
    fConeNClusterSubEtaBandPerMinCut    [icut] = 0;
    
    fConeptsumPhiBandClusterPerMinCut   [icut] = 0;
    fConeNPhiBandClusterPerMinCut       [icut] = 0;
    fConeptsumClusterSubPhiBandPerMinCut[icut] = 0;
    fConeNClusterSubPhiBandPerMinCut    [icut] = 0;
  }
}

//_____________________________________________________________________________
/// Fill some histograms to understand pile-up.
/// Remember to relax time cuts in the reader for time related histograms, filled only for isolated clusters
//_____________________________________________________________________________
void AliAnaParticleIsolation::FillPileUpHistograms(AliCaloTrackParticleCorrelation* pCandidate)
{ 
  Bool_t  isolated   = pCandidate->IsIsolated();
  Float_t energy     = pCandidate->E();
  Float_t pt         = pCandidate->Pt();
  Float_t weightTrig = pCandidate->GetWeight();
  Float_t time       = pCandidate->GetTime();

  if(GetReader()->IsPileUpFromSPD())               fhPtPileUp[0][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromEMCal())             fhPtPileUp[1][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtPileUp[2][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtPileUp[3][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtPileUp[4][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtPileUp[5][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtPileUp[6][isolated]->Fill(pt, GetEventWeight()*weightTrig) ; 
  
  // Fill histo now only for isolated clusters
  if ( !isolated ) return;

  AliVEvent * event = GetReader()->GetInputEvent();
  
  fhTimeENoCut->Fill(energy, time, GetEventWeight());
  if(GetReader()->IsPileUpFromSPD())     fhTimeESPD     ->Fill(energy, time, GetEventWeight());
  if(event->IsPileupFromSPDInMultBins()) fhTimeESPDMulti->Fill(energy, time, GetEventWeight());
  
  if(energy < 8) return; // Fill time figures for high energy clusters not too close to trigger threshold
  
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
  
  // N pile up vertices
  Int_t nVerticesSPD    = -1;
  Int_t nVerticesTracks = -1;
  
  if      (esdEv)
  {
    nVerticesSPD    = esdEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = esdEv->GetNumberOfPileupVerticesTracks();
    
  }//ESD
  else if (aodEv)
  {
    nVerticesSPD    = aodEv->GetNumberOfPileupVerticesSPD();
    nVerticesTracks = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  fhTimeNPileUpVertSPD  ->Fill(time, nVerticesSPD   , GetEventWeight());
  fhTimeNPileUpVertTrack->Fill(time, nVerticesTracks, GetEventWeight());
  
  //printf("Is SPD %d, Is SPD Multi %d, n spd %d, n track %d\n",
  //       GetReader()->IsPileUpFromSPD(),event->IsPileupFromSPDInMultBins(),nVerticesSPD,nVerticesTracks);
  
  Int_t ncont = -1;
  Float_t z1 = -1, z2 = -1;
  Float_t diamZ = -1;
  for(Int_t iVert=0; iVert<nVerticesSPD;iVert++)
  {
    if      (esdEv)
    {
      const AliESDVertex* pv=esdEv->GetPileupVertexSPD(iVert);
      ncont=pv->GetNContributors();
      z1 = esdEv->GetPrimaryVertexSPD()->GetZ();
      z2 = pv->GetZ();
      diamZ = esdEv->GetDiamondZ();
    }//ESD
    else if (aodEv)
    {
      AliAODVertex *pv=aodEv->GetVertex(iVert);
      if(pv->GetType()!=AliAODVertex::kPileupSPD) continue;
      ncont=pv->GetNContributors();
      z1=aodEv->GetPrimaryVertexSPD()->GetZ();
      z2=pv->GetZ();
      diamZ = aodEv->GetDiamondZ();
    }// AOD
    
    Double_t distZ  = TMath::Abs(z2-z1);
    diamZ  = TMath::Abs(z2-diamZ);
    
    fhTimeNPileUpVertContributors  ->Fill(time, ncont, GetEventWeight());
    fhTimePileUpMainVertexZDistance->Fill(time, distZ, GetEventWeight());
    fhTimePileUpMainVertexZDiamond ->Fill(time, diamZ, GetEventWeight());
    
  }// loop
}

//______________________________________________________________
/// Fill shower Shape control histograms and other like decay tagged clusters
///  and NLM dependent histograms
//_____________________________________________________________
void AliAnaParticleIsolation::FillShowerShapeControlHistograms
(AliCaloTrackParticleCorrelation  *pCandidate, Int_t mcIndex, Int_t noverlaps)
{
  if( !fFillSSHisto && !fFillTaggedDecayHistograms) return;
  
  Bool_t  isolated   = pCandidate->IsIsolated();
  Float_t m02        = pCandidate->GetM02() ;
  Float_t pt         = pCandidate->Pt();
  Float_t eta        = pCandidate->Eta();
  Float_t phi        = GetPhi(pCandidate->Phi());
  Int_t   iSM        = pCandidate->GetSModNumber();
  Int_t   nMaxima    = pCandidate->GetNLM();
  Int_t   mcTag      = pCandidate->GetTag() ;
  Float_t weightTrig = pCandidate->GetWeight();

  // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
  //
  if ( fFillTaggedDecayHistograms )
  {
    Int_t decayTag = pCandidate->DecayTag();
    if(decayTag < 0) decayTag = 0;
    
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      if(!GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit])) continue;
      
      if ( fFillSSHisto ) fhPtLambda0Decay[isolated][ibit]->Fill(pt, m02, GetEventWeight()*weightTrig);
      
      // In case it was not done on the trigger selection task
      // apply here a shower shape cut to select photons
      if ( m02 > fM02Narrow[1] || m02 < fM02Narrow[0] ) continue;
      
      fhPtDecay    [isolated][ibit]->Fill(pt,       GetEventWeight()*weightTrig);
      fhEtaPhiDecay[isolated][ibit]->Fill(eta, phi, GetEventWeight()*weightTrig);
      
      if ( IsDataMC() && mcIndex < fNumberMCParticleCases && !GetReader()->AreMCPromptPhotonsSelected() )
      {
        fhPtDecayMC[isolated][ibit][mcIndex]->Fill(pt, GetEventWeight()*weightTrig);
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
          fhPtDecayMC[isolated][ibit][kmcPhoton]->Fill(pt, GetEventWeight()*weightTrig);
        
        if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) && 
            fNumberMCParticleCases > kmcPi0DecayLostPair)
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtDecayMC[isolated][ibit][kmcPi0DecayLostPair]->Fill(pt, GetEventWeight()*weightTrig);
          else if( mcIndex == kmcEtaDecay ) fhPtDecayMC[isolated][ibit][kmcEtaDecayLostPair]->Fill(pt, GetEventWeight()*weightTrig);
        }
      }
    } // bit loop
  } // decay histograms
  
  //
  // Shower shape dependent histograms
  //
  if ( fFillSSHisto )
  {
    fhPtLambda0[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);

    if ( IsHighMultiplicityAnalysisOn() )
      fhPtLambda0Cent[isolated]->Fill(pt, m02, GetEventCentrality(), GetEventWeight()*weightTrig);

    if ( fFillPerSMHistograms )
    {
      fhPtPerSM       [isolated]     ->Fill(pt,iSM, GetEventWeight()*weightTrig);
      fhPtLambda0PerSM[isolated][iSM]->Fill(pt,m02, GetEventWeight()*weightTrig);
      
      if ( fStudyNCellsCut )
      {
        if ( fNCellsWithWeight > 4 ) 
          fhPtLambda0PerSMNCellCut[isolated][iSM]->Fill(pt, m02, GetEventWeight()*weightTrig);
        
        fhPtNCellPerSM[isolated][iSM]->Fill(pt, fNCellsWithWeight, GetEventWeight()*weightTrig);
        
        if      ( m02 > 0.1 && m02 <= 0.3 ) fhPtNCellLowM02PerSM [isolated][iSM]->Fill(pt, fNCellsWithWeight, GetEventWeight()*weightTrig);
        else if ( m02 > 0.5 && m02 <= 2   ) fhPtNCellHighM02PerSM[isolated][iSM]->Fill(pt, fNCellsWithWeight, GetEventWeight()*weightTrig);
      }
    }

    if ( fFillPerTCardIndexHistograms )
    {
      fhPtPerTCardIndex       [isolated]             ->Fill(pt,fTCardIndex, GetEventWeight()*weightTrig);
      fhPtLambda0PerTCardIndex[isolated][fTCardIndex]->Fill(pt,m02, GetEventWeight()*weightTrig);
    }
    
    //
    // MC
    //
    if ( IsDataMC() && mcIndex < fNumberMCParticleCases && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
      {
        fhPtLambda0MC[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        if ( fNCellsWithWeight > 4 &&  fStudyNCellsCut )
          fhPtLambda0MCNCellCut[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      }
      
      if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) &&
           fNumberMCParticleCases > kmcPi0DecayLostPair )
      {
        if     ( mcIndex == kmcPi0Decay ) 
          fhPtLambda0MC[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        else if( mcIndex == kmcEtaDecay ) 
          fhPtLambda0MC[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);

        if ( fNCellsWithWeight > 4 && fStudyNCellsCut )
        {
          if     ( mcIndex == kmcPi0Decay ) 
            fhPtLambda0MCNCellCut[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          else if( mcIndex == kmcEtaDecay ) 
            fhPtLambda0MCNCellCut[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig); 
        }
      } // lost
      
      fhPtLambda0MC[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      if ( fNCellsWithWeight > 4 && fStudyNCellsCut )
        fhPtLambda0MCNCellCut[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);

      if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
      {
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtLambda0MCConv[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) &&
            fNumberMCParticleCases > kmcPi0DecayLostPair )
        {
          if     ( mcIndex == kmcPi0Decay ) 
            fhPtLambda0MCConv[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          else if( mcIndex == kmcEtaDecay ) 
            fhPtLambda0MCConv[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        }
        
        fhPtLambda0MCConv[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      } // Conversion
      
      //
      // Check overlaps
      //
      if ( fFillOverlapHistograms && noverlaps >= 0 )
      {
        if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtNOverlap[kmcPhoton][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
        
        if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) &&
             fNumberMCParticleCases > kmcPi0DecayLostPair )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtNOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
          else if( mcIndex == kmcEtaDecay ) fhPtNOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
        }
        
        fhPtNOverlap[mcIndex][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
        {
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtNOverlapConv[kmcPhoton][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) 
              fhPtNOverlapConv[kmcPi0DecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
            else if( mcIndex == kmcEtaDecay ) 
              fhPtNOverlapConv[kmcEtaDecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
          }
          
          fhPtNOverlapConv[mcIndex][isolated]->Fill(pt, noverlaps, GetEventWeight()*weightTrig);
        } // Conversion
        
        if ( noverlaps == 1 )
        {
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtLambda0MCWith1Overlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) &&
               fNumberMCParticleCases > kmcPi0DecayLostPair )
          {
            if      ( mcIndex == kmcPi0Decay ) 
              fhPtLambda0MCWith1Overlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            else if ( mcIndex == kmcEtaDecay ) 
              fhPtLambda0MCWith1Overlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          }
          
          fhPtLambda0MCWith1Overlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
          {
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtLambda0MCConvWith1Overlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
            {
              if      ( mcIndex == kmcPi0Decay ) 
                fhPtLambda0MCConvWith1Overlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
              else if ( mcIndex == kmcEtaDecay ) 
                fhPtLambda0MCConvWith1Overlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            }
            
            fhPtLambda0MCConvWith1Overlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          } // Conversion
        } // At least 1 overlap
        else if (noverlaps == 0 ) // No overlap
        {
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtLambda0MCWithNoOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) &&
              fNumberMCParticleCases > kmcPi0DecayLostPair )
          {
            if     ( mcIndex == kmcPi0Decay ) 
              fhPtLambda0MCWithNoOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            else if( mcIndex == kmcEtaDecay ) 
              fhPtLambda0MCWithNoOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          }
          
          fhPtLambda0MCWithNoOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
          {
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtLambda0MCConvWithNoOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
            {
              if     ( mcIndex == kmcPi0Decay ) 
                fhPtLambda0MCConvWithNoOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
              else if( mcIndex == kmcEtaDecay ) 
                fhPtLambda0MCConvWithNoOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
            }
            
            fhPtLambda0MCConvWithNoOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
          } // Conversion
        } // more than 1 overlap
      }
    } // MC
    
    if ( GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 &&
         GetModuleNumber(pCandidate) >= GetFirstSMCoveredByTRD()  )
    {
      fhPtLambda0TRD[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
    }
    
    if ( fFillNLMHistograms )
    {
      if     (nMaxima==1)
      {
        fhPtLambda0LocMax1[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        fhPtLambda1LocMax1[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      }
      else if(nMaxima==2)
      {
        fhPtLambda0LocMax2[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        fhPtLambda1LocMax2[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      }
      else
      {
        fhPtLambda0LocMaxN[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
        fhPtLambda1LocMaxN[isolated]->Fill(pt, m02, GetEventWeight()*weightTrig);
      }
    }
  } // SS histo fill
  
}

//_____________________________________________________________________________________________________________________
/// Fill Track matching and Shower Shape control histograms.
//_____________________________________________________________________________________________________________________
void AliAnaParticleIsolation::FillTrackMatchingControlHistograms
 (AliCaloTrackParticleCorrelation  *pCandidate, Int_t mcIndex)
{
  // Track matching dependent histograms
  if ( !fFillTMHisto || !fCluster ) return;
  
  Bool_t  isolated   = pCandidate->IsIsolated();
  Float_t pt         = pCandidate->Pt();
  Float_t energy     = pCandidate->E();
  Int_t   mcTag      = pCandidate->GetTag() ;
  Float_t weightTrig = pCandidate->GetWeight();

  Float_t dZ  = fCluster->GetTrackDz();
  Float_t dR  = fCluster->GetTrackDx();
  
  //    if(fCluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
  //    {
  //      dR = 2000., dZ = 2000.;
  //      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(fCluster->GetID(),dZ,dR);
  //    }
  
  //printf("ParticleIsolation: dPhi %f, dEta %f\n",dR,dZ);
  if ( fhTrackMatchedDEta[isolated] && TMath::Abs(dR) < 999 )
  {
    fhTrackMatchedDEta[isolated]->Fill(pt, dZ, GetEventWeight()*weightTrig);
    fhTrackMatchedDPhi[isolated]->Fill(pt, dR, GetEventWeight()*weightTrig);
    if ( pt > 0.5 ) 
      fhTrackMatchedDEtaDPhi[isolated]->Fill(dZ, dR, GetEventWeight()*weightTrig);
   
    if ( IsDataMC() && mcIndex < fNumberMCParticleCases && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      fhTrackMatchedDEtaMC[mcIndex][isolated]->Fill(pt, dZ, GetEventWeight()*weightTrig);
      fhTrackMatchedDPhiMC[mcIndex][isolated]->Fill(pt, dR, GetEventWeight()*weightTrig);
      if ( pt > 0.5 ) 
        fhTrackMatchedDEtaDPhiMC[mcIndex][isolated]->Fill(dZ, dR, GetEventWeight()*weightTrig);
    }
  }
  
  // Check dEdx and E/p of matched clusters
  
  if ( TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05 )
  {
    AliVTrack *track = GetCaloUtils()->GetMatchedTrack(fCluster, GetReader()->GetInputEvent());
    
    if(track)
    {
      Float_t dEdx = track->GetTPCsignal();
      fhdEdx[isolated]->Fill(pt, dEdx, GetEventWeight()*weightTrig);
      
      Float_t eOverp = energy/track->P();
      fhEOverP[isolated]->Fill(pt,  eOverp, GetEventWeight()*weightTrig);
    }
    //else
    //  printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
    
    if ( IsDataMC() && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      if ( !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion)  )
      {
        if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                  GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 2.5, GetEventWeight()*weightTrig);
        else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 0.5, GetEventWeight()*weightTrig);
        else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 1.5, GetEventWeight()*weightTrig);
        else                                                                                   
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 3.5, GetEventWeight()*weightTrig);
      }
      else
      {
        if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                  GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 6.5, GetEventWeight()*weightTrig);
        else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 4.5, GetEventWeight()*weightTrig);
        else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) 
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 5.5, GetEventWeight()*weightTrig);
        else                                                                                   
          fhTrackMatchedMCParticle[isolated]->Fill(pt, 7.5, GetEventWeight()*weightTrig);
      }
    }  // MC
  } // match window
}

//______________________________________________________
/// Save parameters used for analysis.
//______________________________________________________
TObjString *  AliAnaParticleIsolation::GetAnalysisCuts()
{
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar, buffersize,"--- AliAnaParticleIsolation ---:") ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"Cand. Detector: %s;",fIsoDetectorString.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillTMHisto=%d;",fFillTMHisto) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillSSHisto=%d;",fFillSSHisto) ;
  parList+=onePar ;
  
  // Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  // Get parameters set in IC class.
  parList += GetIsolationCut()->GetICParametersList() ;
  
  return new TObjString(parList) ;
}

//________________________________________________________
/// Create histograms to be saved in output file and
/// store them in outputContainer.
//________________________________________________________
TList *  AliAnaParticleIsolation::GetCreateOutputObjects()
{  
  TList * outputContainer = new TList() ;
  outputContainer->SetName("IsolatedParticleHistos") ;
  
  //Float_t ptthre    = GetIsolationCut()->GetPtThreshold();
  //Float_t ptsumthre = GetIsolationCut()->GetSumPtThreshold();
  //Float_t ptfrac    = GetIsolationCut()->GetPtFraction();
  Float_t r         = GetIsolationCut()->GetConeSize();
  Int_t   method    = GetIsolationCut()->GetICMethod() ;
  Int_t   particle  = GetIsolationCut()->GetParticleTypeInCone() ;
  
  // Histogram ranges and bins
  //
  Int_t   nptbins  = GetHistogramRanges()->GetHistoPtBins();
  Int_t   nphibins = GetHistogramRanges()->GetHistoPhiBins();
  Int_t   netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax    = GetHistogramRanges()->GetHistoPtMax();
  Float_t phimax   = GetHistogramRanges()->GetHistoPhiMax();
  Float_t etamax   = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin    = GetHistogramRanges()->GetHistoPtMin();
  Float_t phimin   = GetHistogramRanges()->GetHistoPhiMin();
  Float_t etamin   = GetHistogramRanges()->GetHistoEtaMin();
  Int_t   ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();
  Float_t ssmax    = GetHistogramRanges()->GetHistoShowerShapeMax();
  Float_t ssmin    = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t   ntimebins= GetHistogramRanges()->GetHistoTimeBins();
  Float_t timemax  = GetHistogramRanges()->GetHistoTimeMax();
  Float_t timemin  = GetHistogramRanges()->GetHistoTimeMin();
  
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  Int_t   ndedxbins   = GetHistogramRanges()->GetHistodEdxBins();
  Float_t dedxmax     = GetHistogramRanges()->GetHistodEdxMax();
  Float_t dedxmin     = GetHistogramRanges()->GetHistodEdxMin();
  Int_t   nPoverEbins = GetHistogramRanges()->GetHistoEOverPBins();
  Float_t pOverEmax   = GetHistogramRanges()->GetHistoEOverPMax();
  Float_t pOverEmin   = GetHistogramRanges()->GetHistoEOverPMin();
  
  Int_t   nptsumbins    = GetHistogramRanges()->GetHistoNPtSumBins();
  Float_t ptsummax      = GetHistogramRanges()->GetHistoPtSumMax();
  Float_t ptsummin      = GetHistogramRanges()->GetHistoPtSumMin();
  Int_t   nptinconebins = GetHistogramRanges()->GetHistoNPtInConeBins();
  Float_t ptinconemax   = GetHistogramRanges()->GetHistoPtInConeMax();
  Float_t ptinconemin   = GetHistogramRanges()->GetHistoPtInConeMin();
  
  Int_t   ncenbin  = GetHistogramRanges()->GetHistoCentralityBins()  ;
  Float_t cenmin   = GetHistogramRanges()->GetHistoCentralityMin()   ;
  Float_t cenmax   = GetHistogramRanges()->GetHistoCentralityMax()   ;
  
  Int_t   nmultbin = GetHistogramRanges()->GetHistoTrackMultiplicityBins();
  Int_t   multmax  = GetHistogramRanges()->GetHistoTrackMultiplicityMax ();
  Int_t   multmin  = GetHistogramRanges()->GetHistoTrackMultiplicityMin ();
 
  Int_t   nclbin   = GetHistogramRanges()->GetHistoNClustersBins();
  Int_t   nclmax   = GetHistogramRanges()->GetHistoNClustersMax ();
  Int_t   nclmin   = GetHistogramRanges()->GetHistoNClustersMin ();
  
  Int_t cellBins   = GetHistogramRanges()->GetHistoNClusterCellBins(); 
  Int_t cellMax    = GetHistogramRanges()->GetHistoNClusterCellMax(); 
  Int_t cellMin    = GetHistogramRanges()->GetHistoNClusterCellMin();
  
  Int_t   nnovbin  = GetHistogramRanges()->GetHistoNoverlapBins();
  Float_t novmin   = GetHistogramRanges()->GetHistoNoverlapMin() ;
  Float_t novmax   = GetHistogramRanges()->GetHistoNoverlapMax() ;
  
  // TH3 histograms
  InitHistoRangeArrays();

  TArrayD ptBinsArray     = GetHistogramRanges()->GetHistoPtArr();
  TArrayD ptWideBinsArray = GetHistogramRanges()->GetHistoWidePtArr();
  TArrayD cenBinsArray    = GetHistogramRanges()->GetHistoCentralityArr();
  TArrayD etaBinsArray    = GetHistogramRanges()->GetHistoEtaArr();
  TArrayD phiBinsArray    = GetHistogramRanges()->GetHistoPhiArr();
  TArrayD ssBinsArray     = GetHistogramRanges()->GetHistoShowerShapeArr();
  TArrayD nceBinsArray    = GetHistogramRanges()->GetHistoNClusterCellArr();
  TArrayD nclBinsArray    = GetHistogramRanges()->GetHistoNClustersArr();
  TArrayD mulBinsArray    = GetHistogramRanges()->GetHistoTrackMultiplicityArr();
  TArrayD sum0BinsArray   = GetHistogramRanges()->GetHistoPtSumArr();
  TArrayD sueBinsArray    = GetHistogramRanges()->GetHistoPtSumSubArr();
  TArrayD ratBinsArray    = GetHistogramRanges()->GetHistoRatioArr();
  TArrayD sumBinsArray;
  if ( method >= AliIsolationCut::kSumBkgSubIC ) 
    sumBinsArray = GetHistogramRanges()->GetHistoPtSumSubArr();
  else 
    sumBinsArray = GetHistogramRanges()->GetHistoPtSumArr();

  TCustomBinning  minPtBinning;
  minPtBinning.SetMinimum(-0.5);
  minPtBinning.AddStep(fNPtCutsInCone, 1); 
  TArrayD minPtBinsArray;
  minPtBinning.CreateBinEdges(minPtBinsArray);
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  // Strings for histograms names, titles
  //
  TString sThreshold[] = {"",""};
  TString ueType = "";
  if      ( method == AliIsolationCut::kSumPtIC ||  
            method >= AliIsolationCut::kSumBkgSubIC )
  {
    sThreshold[0] = Form(", %2.2f < #Sigma #it{p}_{T}^{in cone} < %2.2f GeV/#it{c}",
                         GetIsolationCut()->GetSumPtThreshold()+GetIsolationCut()->GetSumPtThresholdGap(),
                         GetIsolationCut()->GetSumPtThresholdMax());
    if ( GetIsolationCut()->GetSumPtThresholdMax() > 200 )
      sThreshold[0] = Form(", #Sigma #it{p}_{T}^{in cone} > %2.2f GeV/#it{c}",
                           GetIsolationCut()->GetSumPtThreshold()+GetIsolationCut()->GetSumPtThresholdGap());
    
    sThreshold[1] = Form(", #Sigma #it{p}_{T}^{in cone} < %2.2f GeV/#it{c}",
                         GetIsolationCut()->GetSumPtThreshold());
    
    if      ( method == AliIsolationCut::kSumBkgSubIC ) 
      ueType = ", UE #perp cones" ; 
    else if ( method == AliIsolationCut::kSumBkgSubEtaBandIC ) 
      ueType = ", UE #eta band"   ;  
    else if ( method == AliIsolationCut::kSumBkgSubPhiBandIC ) 
      ueType = ", UE #varphi band";  
    
    sThreshold[0]+=ueType;
    sThreshold[1]+=ueType;
  }
  else if ( method == AliIsolationCut::kPtThresIC)
  {
    sThreshold[0] = Form(", %2.2f < #it{p}_{T}^{th} < %2.2f GeV/#it{c}",
                         GetIsolationCut()->GetPtThreshold(),GetIsolationCut()->GetPtThresholdMax());
    if(GetIsolationCut()->GetSumPtThreshold() > 200)
      sThreshold[0] = Form(", #it{p}_{T}^{th} = %2.2f GeV/#it{c}",
                           GetIsolationCut()->GetPtThreshold());
    
    sThreshold[1] = sThreshold[0];
  }
  else if ( method == AliIsolationCut::kPtFracIC)
  {
    sThreshold[0] = Form(", #Sigma #it{p}_{T}^{in cone}/#it{p}_{T}^{trig} = %2.2f" ,
                         GetIsolationCut()->GetPtFraction());
    sThreshold[1] = sThreshold[0];
  }
  
  TString sParticle = ", x^{  0,#pm}";
  if      ( particle == AliIsolationCut::kOnlyNeutral )  sParticle = ", x^{ 0}";
  else if ( particle == AliIsolationCut::kOnlyCharged )  sParticle = ", x^{  #pm}";
  
  TString parTitle[2];
  parTitle[0] = Form("#it{R} = %2.2f%s%s",GetIsolationCut()->GetConeSize(),sThreshold[0].Data(),sParticle.Data());
  parTitle[1] = Form("#it{R} = %2.2f%s%s",GetIsolationCut()->GetConeSize(),sThreshold[1].Data(),sParticle.Data());
  TString parTitleR   = Form("#it{R} = %2.2f%s%s"        ,GetIsolationCut()->GetConeSize(),sParticle.Data(),ueType.Data());
  TString parTitleRCh = Form("#it{R} = %2.2f, x^{ #pm}%s",GetIsolationCut()->GetConeSize(),ueType.Data());
  TString parTitleRNe = Form("#it{R} = %2.2f, x^{ 0}%s"  ,GetIsolationCut()->GetConeSize(),ueType.Data());
  
  if( GetIsolationCut()->GetMinDistToTrigger() > 0 )
  {
    parTitle[0]  = Form("%2.2f<#it{R}<%2.2f%s%s"     ,
                        GetIsolationCut()->GetMinDistToTrigger(),GetIsolationCut()->GetConeSize(),sThreshold[0].Data(),sParticle.Data());
    parTitle[1]  = Form("%2.2f<#it{R}<%2.2f%s%s"     ,
                        GetIsolationCut()->GetMinDistToTrigger(),GetIsolationCut()->GetConeSize(),sThreshold[1].Data(),sParticle.Data());
    
    parTitleR   = Form("%2.2f<#it{R}<%2.2f%s"       ,
                       GetIsolationCut()->GetMinDistToTrigger(),GetIsolationCut()->GetConeSize(),sParticle.Data());
    parTitleRCh = Form("%2.2f<#it{R}<%2.2f, x^{#pm}",
                       GetIsolationCut()->GetMinDistToTrigger(),GetIsolationCut()->GetConeSize());
    parTitleRNe = Form("%2.2f<#it{R}<%2.2f, x^{0}"  ,
                       GetIsolationCut()->GetMinDistToTrigger(),GetIsolationCut()->GetConeSize());
  }
  
  TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
  
  // MC histograms title and name
  TString mcPartType[] = { "#gamma"                 , "#gamma_{prompt}"    , "#gamma_{fragmentation}",
                           "#pi^{0} (merged #gamma)", "#gamma_{#pi decay}" , "#gamma_{#pi decay} lost companion",
                           "#eta (merged #gamma)"   , "#gamma_{#eta decay}", "#gamma_{#eta decay} lost companion",
                           "#gamma_{other decay}"   , "e^{#pm}"            , "hadrons?"} ;
  
  TString mcPartName[] = { "Photon","PhotonPrompt","PhotonFrag",
                           "Pi0"   ,"Pi0Decay"    ,"Pi0DecayLostPair",
                           "Eta"   ,"EtaDecay"    ,"EtaDecayLostPair",
                           "OtherDecay","Electron","Hadron"} ;
  
  // Primary MC histograms title and name
  TString pptype[] = { "#gamma"         , "#gamma_{#pi decay}"    , "#gamma_{#eta decay}", "#gamma_{other decay}",
                       "#gamma_{prompt}", "#gamma_{fragmentation}", "#gamma_{ISR}"       ,
                       "#pi^{0}"        , "#eta"} ;
  
  TString ppname[] = { "Photon"      , "PhotonPi0Decay","PhotonEtaDecay", "PhotonOtherDecay",
                       "PhotonPrompt", "PhotonFrag"    , "PhotonISR"    ,
                       "Pi0"         , "Eta"} ;
  
  TString isoName [] = {"NoIso","Iso"};
  TString isoTitle[] = {"Not isolated"  ,"Isolated"};
  
  TString m02Name [] = {"Wide","Narrow"};
  TString m02Title[] = { Form(", %2.2f < #sigma_{long}^{2} < %2.2f",fM02Wide  [0],fM02Wide  [1]), 
                         Form(", %2.2f < #sigma_{long}^{2} < %2.2f",fM02Narrow[0],fM02Narrow[1]) };
  
  Int_t nShSh = 2;
  if( !fFillSSHisto ) 
  {
    nShSh       = 1;
    m02Name [0] = ""; 
    m02Name [1] = ""; 
    m02Title[0] = ""; 
    m02Title[1] = ""; 
  }
  
  // Reference histograms
  
  for(Int_t iso = 0; iso < 2; iso++)
  {
    for(Int_t ishsh = 0; ishsh < nShSh; ishsh++)
    {      
      if ( fFillOnlyTH3Histo ) continue;
      
      fhPt[iso][ishsh]  = new TH1F
      (Form("hPt%s%s",isoName[iso].Data(),m02Name[ishsh].Data()),
       Form("%s%s, %s",
            isoTitle[iso].Data(), m02Title[ishsh].Data(), parTitle[iso].Data()),
       nptbins,ptmin,ptmax);
      fhPt[iso][ishsh]->SetYTitle("#it{counts}");
      fhPt[iso][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPt[iso][ishsh]) ;
      
      fhPtEtaPhi[iso][ishsh]  = new TH3F
      (Form("hPtEtaPhi%s%s",isoName[iso].Data(), m02Name[ishsh].Data()),
       Form("%s%s, %s",
            isoTitle[iso].Data(), m02Title[ishsh].Data(), parTitle[iso].Data()),
       ptWideBinsArray.GetSize() - 1, ptWideBinsArray.GetArray(),
          etaBinsArray.GetSize() - 1,    etaBinsArray.GetArray(),      
          phiBinsArray.GetSize() - 1,    phiBinsArray.GetArray());
      fhPtEtaPhi[iso][ishsh]->SetYTitle("#eta");
      fhPtEtaPhi[iso][ishsh]->SetZTitle("#varphi (rad)");
      fhPtEtaPhi[iso][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtEtaPhi[iso][ishsh]) ;
      
      if ( IsDataMC() && !GetReader()->AreMCPromptPhotonsSelected())
      {
        // For histograms in arrays, index in the array, corresponding to any particle origin
        
        for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
        {
          if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
          {
            if ( mcPartName[imc].Contains("Prompt") ) continue;
            if ( mcPartName[imc].Contains("Frag"  ) ) continue;
            if ( mcPartName[imc].Contains("ISR"   ) ) continue;
          }
          
          fhPtMC[imc][iso][ishsh]  = new TH1F
          (Form("hPt%s%sMC%s",isoName[iso].Data(),m02Name[ishsh].Data(),mcPartName[imc].Data()),
           Form("%s%s, %s %s",
                isoTitle[iso].Data(), m02Title[ishsh].Data(), mcPartType[imc].Data(), parTitle[iso].Data()),
           nptbins,ptmin,ptmax);
          fhPtMC[imc][iso][ishsh]->SetYTitle("#it{counts}");
          fhPtMC[imc][iso][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtMC[imc][iso][ishsh]) ;
          
          if ( ishsh > 0 ) continue;
          
          fhPtEtaPhiMC[imc][iso]  = new TH3F
          (Form("hPtEtaPhi%sMC%s",isoName[iso].Data(), mcPartName[imc].Data()),
           Form("%s %s, %s",
                isoTitle[iso].Data(), mcPartType[imc].Data(), parTitle[iso].Data()),
           ptWideBinsArray.GetSize() - 1, ptWideBinsArray.GetArray(),
              etaBinsArray.GetSize() - 1,    etaBinsArray.GetArray(),      
              phiBinsArray.GetSize() - 1,    phiBinsArray.GetArray());
          fhPtEtaPhiMC[imc][iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhPtEtaPhiMC[imc][iso]->SetYTitle("#eta");
          fhPtEtaPhiMC[imc][iso]->SetZTitle("#varphi (rad)");
          outputContainer->Add(fhPtEtaPhiMC[imc][iso]) ;
        }
      }
    }
    
    if ( fStudyExoticTrigger )
    {
      fhPtExoTrigger[iso]  = new TH1F
      (Form("hPt%sExoTrigger",isoName[iso].Data()),
       Form("%s %s, #it{F}_{+}>0.97",isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax);
      fhPtExoTrigger[iso]->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtExoTrigger[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtExoTrigger[iso]) ;
    }
    
    // Histograms for tagged candidates as decay
    if ( fFillTaggedDecayHistograms )
    {
      for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
      {        
        fhPtDecay[iso][ibit]  =
        new TH1F(Form("hPtDecay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
                 Form("Number of %s leading pi0 decay particles vs #it{p}_{T}, bit %d, %s",
                      isoTitle[iso].Data(),fDecayBits[ibit],parTitle[iso].Data()),
                 nptbins,ptmin,ptmax);
        fhPtDecay[iso][ibit]->SetYTitle("#it{counts}");
        fhPtDecay[iso][ibit]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtDecay[iso][ibit]) ;
        
        fhEtaPhiDecay[iso][ibit]  =
        new TH2F(Form("hEtaPhiDecay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
                 Form("Number of %s leading Pi0 decay particles #eta vs #varphi, bit %d, %s",
                      isoTitle[iso].Data(),fDecayBits[ibit],parTitle[iso].Data()),
                 netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiDecay[iso][ibit]->SetXTitle("#eta");
        fhEtaPhiDecay[iso][ibit]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiDecay[iso][ibit]) ;
        
        if ( fFillSSHisto )
        {
          fhPtLambda0Decay[iso][ibit]  = new TH2F
          (Form("hPtLambda0Decay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
           Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}, decay bit %d, %s",
                isoTitle[iso].Data(), fDecayBits[ibit], parTitle[iso].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0Decay[iso][ibit]->SetYTitle("#sigma_{long}^{2}");
          fhPtLambda0Decay[iso][ibit]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0Decay[iso][ibit]) ;
        }
        
        if ( IsDataMC() && !GetReader()->AreMCPromptPhotonsSelected() )
        {
          for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
          {
            if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
            {
              if ( mcPartName[imc].Contains("Prompt") ) continue;
              if ( mcPartName[imc].Contains("Frag"  ) ) continue;
              if ( mcPartName[imc].Contains("ISR"   ) ) continue;
            }
            
            fhPtDecayMC[iso][ibit][imc]  =
            new TH1F(Form("hPtDecay%s_bit%d_MC%s",isoName[iso].Data(),fDecayBits[ibit],mcPartName[imc].Data()),
                     Form("#it{p}_{T} of %s, decay bit %d,  %s, %s",
                          isoTitle[iso].Data(),fDecayBits[ibit],mcPartType[imc].Data(),parTitle[iso].Data()),
                     nptbins,ptmin,ptmax);
            fhPtDecayMC[iso][ibit][imc]->SetYTitle("#it{counts}");
            fhPtDecayMC[iso][ibit][imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtDecayMC[iso][ibit][imc]) ;
          }// MC particle loop
        }// MC
      } // bit loop
    }// decay
  } // iso / no iso
  
  
  if ( method == AliIsolationCut::kSumPtIC || 
       method >= AliIsolationCut::kSumBkgSubIC )
  {
    if ( !IsHighMultiplicityAnalysisOn() )
    {
      fhPtM02SumPtCone = new TH3F
      ("hPtM02SumPtCone",Form("%s",parTitleR.Data()),
       ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
       sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
      fhPtM02SumPtCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtM02SumPtCone->SetYTitle("#sigma_{long}^{2}");
      fhPtM02SumPtCone->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
      outputContainer->Add(fhPtM02SumPtCone) ;
    }
    else 
    {
      //printf("*** N centrality bins %d\n",GetNCentrBin());
      fhPtM02SumPtConeCent = new TH3F*[GetNCentrBin()] ;
      for(Int_t icent = 0; icent < GetNCentrBin(); icent++)
      {
        fhPtM02SumPtConeCent[icent] = new TH3F
        (Form("hPtM02SumPtCone_Cent%d",icent),
         Form("%s, centrality bin %d",parTitleR.Data(), icent),
         ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
        fhPtM02SumPtConeCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtM02SumPtConeCent[icent]->SetYTitle("#sigma_{long}^{2}");
        fhPtM02SumPtConeCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhPtM02SumPtConeCent[icent]) ;
      }
    }
    
    if ( particle == AliIsolationCut::kNeutralAndCharged )
    {
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhPtM02SumPtConeCharged = new TH3F
        ("hPtM02SumPtConeCharged",Form("%s",parTitleRCh.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
          ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
        fhPtM02SumPtConeCharged->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtM02SumPtConeCharged->SetYTitle("#sigma_{long}^{2}");
        fhPtM02SumPtConeCharged->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhPtM02SumPtConeCharged) ;
      }
      else
      {
        //printf("*** N centrality bins %d\n",GetNCentrBin());
        fhPtM02SumPtConeChargedCent = new TH3F*[GetNCentrBin()] ;
        for(Int_t icent = 0; icent < GetNCentrBin(); icent++)
        {
          fhPtM02SumPtConeChargedCent[icent] = new TH3F
          (Form("hPtM02SumPtConeCharged_Cent%d",icent),
           Form("%s, centrality bin %d",parTitleRCh.Data(), icent),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
            ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
          fhPtM02SumPtConeChargedCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhPtM02SumPtConeChargedCent[icent]->SetYTitle("#sigma_{long}^{2}");
          fhPtM02SumPtConeChargedCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
          outputContainer->Add(fhPtM02SumPtConeChargedCent[icent]) ;
        } // cen loop
      } // high mult
    }
    
    if ( IsDataMC() && !IsHighMultiplicityAnalysisOn() && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
      {
        if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
        {
          if ( mcPartName[imc].Contains("Prompt") ) continue;
          if ( mcPartName[imc].Contains("Frag"  ) ) continue;
          if ( mcPartName[imc].Contains("ISR"   ) ) continue;
        }
        
        fhPtM02SumPtConeMC[imc] = new TH3F
        (Form("hPtM02SumPtCone_MC%s",mcPartName[imc].Data()),
         Form("%s, MC %s", parTitleR.Data(), mcPartType[imc].Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
          ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
        fhPtM02SumPtConeMC[imc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtM02SumPtConeMC[imc]->SetYTitle("#sigma_{long}^{2}");
        fhPtM02SumPtConeMC[imc]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhPtM02SumPtConeMC[imc]) ;
        
        if ( particle == AliIsolationCut::kNeutralAndCharged )
        {
          fhPtM02SumPtConeChargedMC[imc] = new TH3F
          (Form("hPtM02SumPtConeCharged_MC%s",mcPartName[imc].Data()),
           Form("%s, MC %s", parTitleRCh.Data(), mcPartType[imc].Data()),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
            ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
          fhPtM02SumPtConeChargedMC[imc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhPtM02SumPtConeChargedMC[imc]->SetYTitle("#sigma_{long}^{2}");
          fhPtM02SumPtConeChargedMC[imc]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
          outputContainer->Add(fhPtM02SumPtConeChargedMC[imc]) ;
        }
      } // MC particle loop 
    } // MC, not high mult
    
    // Centrality dependent MC
    if ( IsDataMC() && IsHighMultiplicityAnalysisOn() && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      //printf("INIT, ncent %d, mc cases %d\n",GetNCentrBin(),fNumberMCParticleCases);
      fhPtM02SumPtConeCentMC = new TH3F*[GetNCentrBin()*fNumberMCParticleCases] ;
      if ( particle == AliIsolationCut::kNeutralAndCharged )
        fhPtM02SumPtConeChargedCentMC = new TH3F*[GetNCentrBin()*fNumberMCParticleCases] ;
      
      for(Int_t icent = 0; icent < GetNCentrBin(); icent++)
      {
        for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
        {
          Int_t index = icent*fNumberMCParticleCases+imc;
          
          if ( GetReader()->AreMCFragmentationPhotonsRejected() )
          {
            if ( mcPartName[imc].Contains("Frag"  ) ) continue;
            if ( mcPartName[imc].Contains("Prompt") ) continue;
            if ( mcPartName[imc].Contains("ISR"   ) ) continue;
          }
          
          //printf("\t icent %d, imc %d, index %d\n",icent,imc,index);

          fhPtM02SumPtConeCentMC[index] = new TH3F
          (Form("hPtM02SumPtCone_Cent%d_MC%s",icent,mcPartName[imc].Data()),
           Form("%s, MC %s, cent bin %d", parTitleR.Data(), mcPartType[imc].Data(),icent),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
            ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
          fhPtM02SumPtConeCentMC[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhPtM02SumPtConeCentMC[index]->SetYTitle("#sigma_{long}^{2}");
          fhPtM02SumPtConeCentMC[index]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
          outputContainer->Add(fhPtM02SumPtConeCentMC[index]) ;
          
          if ( particle == AliIsolationCut::kNeutralAndCharged )
          {
            fhPtM02SumPtConeChargedCentMC[index] = new TH3F
            (Form("hPtM02SumPtConeCharged_Cent%d_MC%s",icent,mcPartName[imc].Data()),
             Form("%s, MC %s, cent bin %d", parTitleRCh.Data(), mcPartType[imc].Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
              ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray()); 
            fhPtM02SumPtConeChargedCentMC[index]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhPtM02SumPtConeChargedCentMC[index]->SetYTitle("#sigma_{long}^{2}");
            fhPtM02SumPtConeChargedCentMC[index]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhPtM02SumPtConeChargedCentMC[index]) ;
          }
        } // MC particle loop
      } // centrality
    } // MC, high mult
    
    for(Int_t ishsh = 0; ishsh < nShSh; ishsh++)
    {
      if ( fFillOnlyTH3Histo ) continue;
      
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeSumPtM02Cut[ishsh] = new TH2F
        (Form("hConeSumPtM02%s",m02Name[ishsh].Data()),
         Form("%s%s",parTitleR.Data(),m02Title[ishsh].Data()),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax); 
        fhConeSumPtM02Cut[ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtM02Cut[ishsh]->SetYTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtM02Cut[ishsh]) ;
      }
      else
      {
        fhConeSumPtCentM02Cut[ishsh] = new TH3F
        (Form("hConeSumPtCentM02%s",m02Name[ishsh].Data()),
         Form("%s%s",parTitleR.Data(),m02Title[ishsh].Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());         
        fhConeSumPtCentM02Cut[ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCentM02Cut[ishsh]->SetYTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        fhConeSumPtCentM02Cut[ishsh]->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtCentM02Cut[ishsh]) ;
      }
      
      if ( IsDataMC() && !IsHighMultiplicityAnalysisOn() && !GetReader()->AreMCPromptPhotonsSelected())
      {
        for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
        {
          if ( GetReader()->AreMCFragmentationPhotonsRejected() )
          {
            if ( mcPartName[imc].Contains("Frag"  ) ) continue;
            if ( mcPartName[imc].Contains("Prompt") ) continue;
            if ( mcPartName[imc].Contains("ISR"   ) ) continue;
          }

          fhConeSumPtM02CutMC[imc][ishsh] = new TH2F
          (Form("hConeSumPtM02%s_MC%s",m02Name[ishsh].Data(),mcPartName[imc].Data()),
           Form("%s%s, MC %s",parTitleR.Data(), m02Title[ishsh].Data(), mcPartType[imc].Data()),
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax); 
          fhConeSumPtM02CutMC[imc][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtM02CutMC[imc][ishsh]->SetYTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtM02CutMC[imc][ishsh]) ;
        }
      } // MC
    } // Fill sum pt vs pT for shower shape cut
  } // Fill TH3 and pT iso shower shape cut
  
  for(Int_t iso = 0; iso < 2; iso++)
  {
    if ( fFillOnlyTH3Histo ) continue;
    
    for(Int_t ishsh = 0; ishsh < nShSh; ishsh++)
    {
      if ( fFillNLMHistograms )
      {
        fhPtNLocMax[iso][ishsh]  = new TH2F
        (Form("hPtNLocMax%s%s",isoName[iso].Data(),m02Name[ishsh].Data()),
         Form("%s%s, %s", isoTitle[iso].Data(), m02Title[ishsh].Data(), parTitle[iso].Data()),
         nptbins,ptmin,ptmax,10,0,10);
        fhPtNLocMax[iso][ishsh]->SetYTitle("#it{n}_{LM}");
        fhPtNLocMax[iso][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtNLocMax[iso][ishsh]) ;
      }
      
      if ( IsHighMultiplicityAnalysisOn() )
      {
        fhPtCentrality[iso][ishsh]  = new TH2F
        (Form("hPtCentrality%s%s",isoName[iso].Data(),m02Name[ishsh].Data()),
         Form("%s%s, %s",parTitle[iso].Data(), m02Title[ishsh].Data(), isoTitle[iso].Data()),
         nptbins,ptmin,ptmax, ncenbin,cenmin,cenmax);
        fhPtCentrality[iso][ishsh]->SetYTitle("centrality");
        fhPtCentrality[iso][ishsh]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
        outputContainer->Add(fhPtCentrality[iso][ishsh]) ;
        
        fhPtEventPlane[iso][ishsh]  = new TH2F
        (Form("hPtEventPlane%s%s",isoName[iso].Data(),m02Name[ishsh].Data()),
         Form("%s%s, %s",parTitle[iso].Data(), m02Title[ishsh].Data(), isoTitle[iso].Data()),
         nptbins,ptmin,ptmax, 100,0,TMath::Pi());
        fhPtEventPlane[iso][ishsh]->SetYTitle("Event plane angle (rad)");
        fhPtEventPlane[iso][ishsh]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtEventPlane[iso][ishsh]) ;
      }
    }
  }
  
  if ( fFillPerSMHistograms )
  {
    for(Int_t ism = 0; ism < fNModules; ism++)
    {
      if ( ism < fFirstModule || ism > fLastModule ) continue;
      
      fhConeSumPtPerSM[ism]  = new TH2F
      (Form("hConePtSum_SM%d",ism),
       Form("Cone track+cluster #Sigma #it{p}_{T} for %s, SM %d",parTitleR.Data(),ism),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtPerSM[ism]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtPerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtPerSM[ism]) ;
      
      fhPtInConePerSM[ism]  = new TH2F
      (Form("hPtInCone_SM%d",ism),
       Form("Cone track+cluster #it{p}_{T} for %s, SM %d", parTitleR.Data(),ism),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConePerSM[ism]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConePerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInConePerSM[ism]) ;
    }
  }
  
  if ( fFillPerTCardIndexHistograms )
  {
    for(Int_t itc = 0; itc < 16; itc++)
    {        
      fhConeSumPtPerTCardIndex[itc]  = new TH2F
      (Form("hConePtSum_TC%d",itc),
       Form("#it{p}_{T}#Sigma #it{p}_{T} for %s, TC index %d",parTitleR.Data(),itc),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtPerTCardIndex[itc]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtPerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtPerTCardIndex[itc]) ;
      
      fhPtInConePerTCardIndex[itc]  = new TH2F
      (Form("hPtInCone_TC%d",itc),
       Form("Cone track+cluster #it{p}_{T} for %s, TC index %d", parTitleR.Data(),itc),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConePerTCardIndex[itc]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConePerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInConePerTCardIndex[itc]) ;
    }
  }
  
  if ( fStudyExoticTrigger )
  {
    fhConeSumPtExoTrigger  = new TH2F
    ("hConePtSumExoTrigger",
     Form("Cone #Sigma #it{p}_{T} for %s, exo trigger",parTitleR.Data()),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPtExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtExoTrigger) ;
    
    fhPtInConeExoTrigger  = new TH2F
    ("hPtInConeExoTrigger",
     Form("Cone #it{p}_{T} for %s, exotic trigger",parTitleR.Data()),
     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInConeExoTrigger) ;
  }
  
  // Cluster only histograms
  if ( GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyCharged )
  {      
    if ( fFillPerSMHistograms )
    {
      for(Int_t ism = 0; ism < fNModules; ism++)
      {
        if ( ism < fFirstModule || ism > fLastModule ) continue;
        
        fhConeSumPtClusterPerSM[ism]  = new TH2F
        (Form("hConePtSumCluster_SM%d",ism),
         Form("Clusters cone #Sigma #it{p}_{T} for %s, SM %d",parTitleR.Data(),ism),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerSM[ism]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterPerSM[ism]) ;
        
        fhPtClusterInConePerSM[ism]  = new TH2F
        (Form("hPtClusterInCone_SM%d",ism),
         Form("Clusters cone #it{p}_{T} for %s, SM %d", parTitleR.Data(),ism),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtClusterInConePerSM[ism]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtClusterInConePerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtClusterInConePerSM[ism]) ;
      }
    }
    
    if ( fFillPerTCardIndexHistograms )
    {
      for(Int_t itc = 0; itc < 16; itc++)
      {
        fhConeSumPtClusterPerTCardIndex[itc]  = new TH2F
        (Form("hConePtSumCluster_TC%d",itc),
         Form("Cluster cone #Sigma #it{p}_{T} %s, TC index %d",parTitleR.Data(),itc),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerTCardIndex[itc]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterPerTCardIndex[itc]) ;
        
        fhPtClusterInConePerTCardIndex[itc]  = new TH2F
        (Form("hPtClusterInCone_TC%d",itc),
         Form("Clusters cone #it{p}_{T} for %s, TC index %d", parTitleR.Data(),itc),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtClusterInConePerTCardIndex[itc]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtClusterInConePerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtClusterInConePerTCardIndex[itc]) ;
      }
    }
    
    if ( fStudyPtCutInCone )
    {
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeNClusterPerMinPtCut = new TH2F
        ("hConeNClusterPerMinPtCut",
         Form("Clusters, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5, nclbin,nclmin,nclmax);
        fhConeNClusterPerMinPtCut->SetYTitle("#it{N}^{cluster}");
        fhConeNClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNClusterPerMinPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterPerMinPtCut) ;
        
        fhConeSumPtClusterPerMinPtCut = new TH2F
        ("hConePtSumClusterPerMinPtCut",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtClusterPerMinPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMinPtCut) ;
        
        // Eta band
        
        fhEtaBandConeNClusterPerMinPtCut = new TH2F
        ("hEtaBandConeNClusterPerMinPtCut",
         Form("Different track #it{p}_{T} cuts in #eta-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5, nclbin,nclmin,nclmax);
        fhEtaBandConeNClusterPerMinPtCut->SetYTitle("#it{N}^{cluster}");
        fhEtaBandConeNClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhEtaBandConeNClusterPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeNClusterPerMinPtCut) ;
        
        fhEtaBandConeSumPtClusterPerMinPtCut = new TH2F
        ("hEtaBandConePtSumClusterPerMinPtCut",
         Form("Cluster #Sigma #it{p}_{T}, different #it{p}_{T} cuts in #eta-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhEtaBandConeSumPtClusterPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhEtaBandConeSumPtClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhEtaBandConeSumPtClusterPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeSumPtClusterPerMinPtCut) ;
        
        fhConeNClusterSubEtaBandPerMinPtCut = new TH2F
        ("hConeNClusterSubEtaBandPerMinPtCut",
         Form("N clusters - UE in #eta-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
        fhConeNClusterSubEtaBandPerMinPtCut->SetYTitle("#it{N}^{cluster}-#it{N}^{cluster}_{UE #eta-band}");
        fhConeNClusterSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNClusterSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterSubEtaBandPerMinPtCut) ;
        
        fhConeSumPtClusterSubEtaBandPerMinPtCut = new TH2F
        ("hConePtSumClusterSubEtaBandPerMinPtCut",
         Form("Cluster cone #Sigma #it{p}_{T} - UE in #eta-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
        fhConeSumPtClusterSubEtaBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
        fhConeSumPtClusterSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtClusterSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterSubEtaBandPerMinPtCut) ;
        
        // Phi band
        
        fhPhiBandConeNClusterPerMinPtCut = new TH2F
        ("hPhiBandConeNClusterPerMinPtCut",
         Form("Different cluster #it{p}_{T} cuts in #varphi-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5, nclbin,nclmin,nclmax);
        fhPhiBandConeNClusterPerMinPtCut->SetYTitle("#it{N}^{cluster}");
        fhPhiBandConeNClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPhiBandConeNClusterPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeNClusterPerMinPtCut) ;
        
        fhPhiBandConeSumPtClusterPerMinPtCut = new TH2F
        ("hPhiBandConePtSumClusterPerMinPtCut",
         Form("Cluster #Sigma #it{p}_{T}-UE, different #it{p}_{T} cuts in #varphi-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPhiBandConeSumPtClusterPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPhiBandConeSumPtClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPhiBandConeSumPtClusterPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeSumPtClusterPerMinPtCut) ;
        
        fhConeNClusterSubPhiBandPerMinPtCut = new TH2F
        ("hConeNClusterSubPhiBandPerMinPtCut",
         Form("N clusters - UE in #varphi-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
        fhConeNClusterSubPhiBandPerMinPtCut->SetYTitle("#it{N}^{cluster}-#it{N}^{cluster}_{UE #varphi-band}");
        fhConeNClusterSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNClusterSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterSubPhiBandPerMinPtCut) ;
        
        fhConeSumPtClusterSubPhiBandPerMinPtCut = new TH2F
        ("hConePtSumClusterSubPhiBandPerMinPtCut",
         Form("Cluster cone #Sigma #it{p}_{T} - UE in #varphi-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
        fhConeSumPtClusterSubPhiBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
        fhConeSumPtClusterSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtClusterSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterSubPhiBandPerMinPtCut) ;
        
        fhConeSumPtClusterPerMaxPtCut = new TH2F
        ("hConePtSumClusterPerMaxPtCut",
         Form("Cluster #Sigma #it{p}_{T}, different max #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMaxPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMaxPtCut->SetXTitle("#it{p}_{T, max} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMaxPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMaxPtCut) ;
      }
      else
      {
        fhConeNClusterPerMinPtCutCent = new TH3F
        ("hConeNClusterPerMinPtCutCent",
         Form("Clusters, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           nclBinsArray.GetSize() - 1,   nclBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray());         
        fhConeNClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNClusterPerMinPtCutCent->SetYTitle("#it{N}^{cluster}");
        fhConeNClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterPerMinPtCutCent) ;
        
        fhConeSumPtClusterPerMinPtCutCent = new TH3F
        ("hConePtSumClusterPerMinPtCutCent",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,  sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtClusterPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMinPtCutCent) ;
        
        // Eta band
        
        fhEtaBandConeNClusterPerMinPtCutCent = new TH3F
        ("hEtaBandConeNClusterPerMinPtCutCent",
         Form("Clusters, different min #it{p}_{T} cut in #eta band for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           nclBinsArray.GetSize() - 1,   nclBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray());           
        fhEtaBandConeNClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhEtaBandConeNClusterPerMinPtCutCent->SetYTitle("#it{N}^{cluster}");
        fhEtaBandConeNClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhEtaBandConeNClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeNClusterPerMinPtCutCent) ;
        
        fhEtaBandConeSumPtClusterPerMinPtCutCent = new TH3F
        ("hEtaBandConePtSumClusterPerMinPtCutCent",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cut in #eta band for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,  sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhEtaBandConeSumPtClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhEtaBandConeSumPtClusterPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhEtaBandConeSumPtClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhEtaBandConeSumPtClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeSumPtClusterPerMinPtCutCent) ;
        
        fhConeNClusterSubEtaBandPerMinPtCutCent = new TH3F
        ("hConeNClusterSubEtaBandPerMinPtCutCent",
         Form("N Clusters - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nclbin*2,-nclmax,nclmax,GetNCentrBin(),0,100);
        fhConeNClusterSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNClusterSubEtaBandPerMinPtCutCent->SetYTitle("#it{N}^{cluster}-#it{N}^{cluster}_{UE #eta-band}");
        fhConeNClusterSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNClusterSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterSubEtaBandPerMinPtCutCent) ;
        
        fhConeSumPtClusterSubEtaBandPerMinPtCutCent = new TH3F
        ("hConePtSumClusterSubEtaBandPerMinPtCutCent",
         Form("Cluster #Sigma #it{p}_{T} - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtClusterSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtClusterSubEtaBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
        fhConeSumPtClusterSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterSubEtaBandPerMinPtCutCent) ;
        
        // Phi band
        
        fhPhiBandConeNClusterPerMinPtCutCent = new TH3F
        ("hPhiBandConeNClusterPerMinPtCutCent",
         Form("Clusters, different min #it{p}_{T} cut in #varphi band cone for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           nclBinsArray.GetSize() - 1,   nclBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray());   
        fhPhiBandConeNClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPhiBandConeNClusterPerMinPtCutCent->SetYTitle("#it{N}^{cluster}");
        fhPhiBandConeNClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPhiBandConeNClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeNClusterPerMinPtCutCent) ;
        
        fhPhiBandConeSumPtClusterPerMinPtCutCent = new TH3F
        ("hPhiBandConePtSumClusterPerMinPtCutCent",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cut in #varphi band for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,   sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhPhiBandConeSumPtClusterPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPhiBandConeSumPtClusterPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPhiBandConeSumPtClusterPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPhiBandConeSumPtClusterPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeSumPtClusterPerMinPtCutCent) ;
        
        fhConeNClusterSubPhiBandPerMinPtCutCent = new TH3F
        ("hConeNClusterSubPhiBandPerMinPtCutCent",
         Form("N Clusters - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nclbin*2,-nclmax,nclmax,GetNCentrBin(),0,100);
        fhConeNClusterSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNClusterSubPhiBandPerMinPtCutCent->SetYTitle("#it{N}^{cluster}-#it{N}^{cluster}_{UE #varphi-band}");
        fhConeNClusterSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNClusterSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterSubPhiBandPerMinPtCutCent) ;
        
        fhConeSumPtClusterSubPhiBandPerMinPtCutCent = new TH3F
        ("hConePtSumClusterSubPhiBandPerMinPtCutCent",
         Form("Cluster #Sigma #it{p}_{T} - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtClusterSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtClusterSubPhiBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
        fhConeSumPtClusterSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterSubPhiBandPerMinPtCutCent) ;
      }
    }
    
    if ( fStudyRCutInCone )
    {
      fhConeSumPtClusterPerRCut = new TH2F
      ("hConePtSumClusterPerRCut","Cluster #Sigma #it{p}_{T}, different #it{R} cuts",
       fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClusterPerRCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterPerRCut->SetXTitle("#it{R}");
      for(Int_t i = 1; i <= fNRCutsInCone; i++)
        fhConeSumPtClusterPerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtClusterPerRCut) ;
      
      fhPtClusterInConePerRCut = new TH2F
      ("hPtClusterInConePerRCut","Cluster #it{p}_{T}, different #it{R} cuts",
       fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtClusterInConePerRCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
      fhPtClusterInConePerRCut->SetXTitle("#it{R}");
      for(Int_t i = 1; i <= fNRCutsInCone; i++)
        fhPtClusterInConePerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
      outputContainer->Add(fhPtClusterInConePerRCut) ;
    }
    
    if ( fStudyNCellsCut )
    {
      fhConeSumPtClusterPerNCellCut = new TH2F
      ("hConePtSumClusterPerNCellCut","Cluster #Sigma #it{p}_{T}, different #it{N}_{cell} cuts",
       fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClusterPerNCellCut->SetYTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
      fhConeSumPtClusterPerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
      for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
        fhConeSumPtClusterPerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
      outputContainer->Add(fhConeSumPtClusterPerNCellCut) ;
      
      fhPtClusterInConePerNCellCut = new TH2F
      ("hPtClusterInConePerNCellCut","Cluster #it{p}_{T}, different #it{N}_{cell} cuts",
       fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtClusterInConePerNCellCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
      fhPtClusterInConePerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
      for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
        fhPtClusterInConePerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
      outputContainer->Add(fhPtClusterInConePerNCellCut) ;
      
      for(Int_t ishsh = 0; ishsh < 4; ishsh++)
      {
        fhPtClusterInConePerNCellPerSM[ishsh] = new TH3F
        (Form("hPtClusterInConePerNCellPerSM_ShSh%d",ishsh),
         Form("Cluster #it{p}_{T} vs #it{n}_{cell} vs SM , 8 < #it{p}_{T}^{trig} < 12 GeV, sh. shape bin %d",ishsh),
         200,0,20,fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax);
        fhPtClusterInConePerNCellPerSM[ishsh]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhPtClusterInConePerNCellPerSM[ishsh]->SetYTitle("SM number");
        fhPtClusterInConePerNCellPerSM[ishsh]->SetXTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        outputContainer->Add(fhPtClusterInConePerNCellPerSM[ishsh]) ;
        
        fhConeSumPtClusterPerNCellPerSM[ishsh] = new TH3F
        (Form("hConeSumPtClusterPerNCellPerSM_ShSh%d",ishsh),
         Form("Cluster #Sigma #it{p}_{T} in cone vs #it{n}_{cell} vs SM , 8 < #it{p}_{T}^{trig} < 12 GeV, sh. shape bin %d",ishsh),
         200,0,50,fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax);
        fhConeSumPtClusterPerNCellPerSM[ishsh]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhConeSumPtClusterPerNCellPerSM[ishsh]->SetYTitle("SM number");
        fhConeSumPtClusterPerNCellPerSM[ishsh]->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterPerNCellPerSM[ishsh]) ;
      }
    }
    
    if ( fStudyExoticTrigger )
    {  
      fhConeSumPtClusterExoTrigger  = new TH2F
      ("hConePtSumClusterExoTrigger",
       Form("Cluster #Sigma #it{p}_{T} in isolation cone for %s, exo trigger",parTitleR.Data()),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClusterExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterExoTrigger) ;
      
      fhConeSumPtClusterPerExoCut = new TH2F
      ("hConePtSumClusterPerExoCut","Cluster #Sigma #it{p}_{T}, different exoticity cuts",
       fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClusterPerExoCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterPerExoCut->SetXTitle("exoticity");
      for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
        fhConeSumPtClusterPerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
      outputContainer->Add(fhConeSumPtClusterPerExoCut) ;
      
      fhPtClusterInConePerExoCut = new TH2F
      ("hPtClusterInConePerExoCut","Cluster #it{p}_{T}, different exoticity cuts",
       fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtClusterInConePerExoCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
      fhPtClusterInConePerExoCut->SetXTitle("exoticity");
      for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
        fhPtClusterInConePerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
      outputContainer->Add(fhPtClusterInConePerExoCut) ;
      
      fhPtClusterInConeExoTrigger  = new TH2F
      ("hPtClusterInConeExoTrigger",
       Form("Cluster cone #it{p}_{T} for %s, exotic trigger",parTitleR.Data()),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtClusterInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtClusterInConeExoTrigger) ;
    }
  }
  
  // Track only histograms
  if ( GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral )
  {
    if ( fFillPerSMHistograms )
    {
      for(Int_t ism = 0; ism < fNModules; ism++)
      {
        if ( ism < fFirstModule || ism > fLastModule ) continue;
        
        fhPtTrackInConePerSM[ism]  = new TH2F
        (Form("hPtTrackInCone_SM%d",ism),
         Form("Track cone #it{p}_{T} for %s, SM %d", parTitleR.Data(),ism),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConePerSM[ism]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtTrackInConePerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConePerSM[ism]) ;
        
        fhConeSumPtTrackPerSM[ism]  = new TH2F
        (Form("hConePtSumTrack_SM%d",ism),
         Form("Track cone #Sigma #it{p}_{T} for %s, SM %d", parTitleR.Data(),ism),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerSM[ism]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerSM[ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtTrackPerSM[ism]) ;
      }
    } 
    
    if ( fFillPerTCardIndexHistograms )
    {
      for(Int_t itc = 0; itc < 16; itc++)
      {            
        fhPtTrackInConePerTCardIndex[itc]  = new TH2F
        (Form("hPtTrackInCone_TC%d",itc),
         Form("Tracks cone #it{p}_{T} for %s, TC index %d", parTitleR.Data(),itc),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConePerTCardIndex[itc]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtTrackInConePerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConePerTCardIndex[itc]) ;
        
        fhConeSumPtTrackPerTCardIndex[itc]  = new TH2F
        (Form("hConePtSumTrack_TC%d",itc),
         Form("Track cone #Sigma #it{p}_{T} for %s, TC index %d",parTitleR.Data(),itc),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerTCardIndex[itc]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerTCardIndex[itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtTrackPerTCardIndex[itc]) ;
      }
    }
    
    if ( IsDataMC() && fFillTrackOriginHistograms && IsGeneratedParticlesAnalysisOn()  )
    {
      TString mcChPartName[] = {"Pion","Kaon","Proton","Other"};
      for(Int_t imc = 0; imc < 4; imc++)
      {
        if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
        {
          if ( mcPartName[imc].Contains("Prompt") ) continue;
          if ( mcPartName[imc].Contains("Frag"  ) ) continue;
          if ( mcPartName[imc].Contains("ISR"   ) ) continue;
        }
        
        fhPtTrackInConeMCPrimary[imc]  = new TH2F
        (Form("hPtTrackInCone_Primary_%s",mcChPartName[imc].Data()),
         Form("Track cone reco #it{p}_{T} for %s, primary MC %s",parTitleR.Data(),mcChPartName[imc].Data()),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeMCPrimary[imc]->SetYTitle("#it{p}_{T in cone}^{reco} (GeV/#it{c})");
        fhPtTrackInConeMCPrimary[imc]->SetXTitle("#it{p}_{T}^{reco} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeMCPrimary[imc]) ;
        
        fhPtTrackInConeMCSecondary[imc]  = new TH2F
        (Form("hPtTrackInCone_Secondary_%s",mcChPartName[imc].Data()),
         Form("Track cone reco #it{p}_{T} for %s, secondary MC %s",parTitleR.Data(),mcChPartName[imc].Data()),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeMCSecondary[imc]->SetYTitle("#it{p}_{T in cone}^{reco} (GeV/#it{c})");
        fhPtTrackInConeMCSecondary[imc]->SetXTitle("#it{p}_{T}^{reco} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeMCSecondary[imc]) ;
        
        fhPtTrackInConeMCPrimaryGener[imc]  = new TH2F
        (Form("hPtTrackInCone_Gener_Primary_%s",mcChPartName[imc].Data()),
         Form("Track cone gener #it{p}_{T} for %s, primary MC %s",parTitleR.Data(),mcChPartName[imc].Data()),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeMCPrimaryGener[imc]->SetYTitle("#it{p}_{T in cone}^{gener} (GeV/#it{c})");
        fhPtTrackInConeMCPrimaryGener[imc]->SetXTitle("#it{p}_{T}^{gener} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeMCPrimaryGener[imc]) ;
        
        fhPtTrackInConeMCSecondaryGener[imc]  = new TH2F
        (Form("hPtTrackInCone_Gener_Secondary_%s",mcChPartName[imc].Data()),
         Form("Track cone gener #it{p}_{T} for %s, secondary MC %s",parTitleR.Data(),mcChPartName[imc].Data()),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeMCSecondaryGener[imc]->SetYTitle("#it{p}_{T in cone}^{gener} (GeV/#it{c})");
        fhPtTrackInConeMCSecondaryGener[imc]->SetXTitle("#it{p}_{T}^{gener} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeMCSecondaryGener[imc]) ;
      }
    }
    
    if ( fStudyPtCutInCone )
    {
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeNTrackPerMinPtCut = new TH2F
        ("hConeNTrackPerMinPtCut",
         Form("N tracks, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
        fhConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackPerMinPtCut) ;
        
        fhConeSumPtTrackPerMinPtCut = new TH2F
        ("hConePtSumTrackPerMinPtCut",
         Form("Track cone #Sigma #it{p}_{T}, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMinPtCut) ;
        
        // Perpendicular cones
        
        fhPerpConeNTrackPerMinPtCut = new TH2F
        ("hPerpConeNTrackPerMinPtCut",
         Form("Different track #it{p}_{T} cuts in #perp cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhPerpConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
        fhPerpConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPerpConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPerpConeNTrackPerMinPtCut) ;
        
        fhPerpConeSumPtTrackPerMinPtCut = new TH2F
        ("hPerpConePtSumTrackPerMinPtCut",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in #perp cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPerpConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPerpConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPerpConeSumPtTrackPerMinPtCut) ;
        
        fhConeNTrackSubPerpPerMinPtCut = new TH2F
        ("hConeNTrackSubPerpPerMinPtCut",
         Form("N tracks - UE in #perp cones, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
        fhConeNTrackSubPerpPerMinPtCut->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #perp}");
        fhConeNTrackSubPerpPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNTrackSubPerpPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubPerpPerMinPtCut) ;
        
        fhConeSumPtTrackSubPerpPerMinPtCut = new TH2F
        ("hConePtSumTrackSubPerpPerMinPtCut",
         Form("Track cone #Sigma #it{p}_{T}- UE in #perp cones, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
        fhConeSumPtTrackSubPerpPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #perp} (GeV/#it{c})");
        fhConeSumPtTrackSubPerpPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtTrackSubPerpPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubPerpPerMinPtCut) ;
        
        // Eta band
        
        fhEtaBandConeNTrackPerMinPtCut = new TH2F
        ("hEtaBandConeNTrackPerMinPtCut",
         Form("Different track #it{p}_{T} cuts in #eta-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhEtaBandConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
        fhEtaBandConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhEtaBandConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeNTrackPerMinPtCut) ;
        
        fhEtaBandConeSumPtTrackPerMinPtCut = new TH2F
        ("hEtaBandConePtSumTrackPerMinPtCut",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in #eta-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhEtaBandConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhEtaBandConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhEtaBandConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeSumPtTrackPerMinPtCut) ;
        
        fhConeNTrackSubEtaBandPerMinPtCut = new TH2F
        ("hConeNTrackSubEtaBandPerMinPtCut",
         Form("N tracks - UE in #eta-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
        fhConeNTrackSubEtaBandPerMinPtCut->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #eta-band}");
        fhConeNTrackSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNTrackSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubEtaBandPerMinPtCut) ;
        
        fhConeSumPtTrackSubEtaBandPerMinPtCut = new TH2F
        ("hConePtSumTrackSubEtaBandPerMinPtCut",
         Form("Track cone #Sigma #it{p}_{T} - UE in #eta-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
        fhConeSumPtTrackSubEtaBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
        fhConeSumPtTrackSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtTrackSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubEtaBandPerMinPtCut) ;
        
        // Phi band
        
        fhPhiBandConeNTrackPerMinPtCut = new TH2F
        ("hPhiBandConeNTrackPerMinPtCut",
         Form("Different track #it{p}_{T} cuts in #varphi-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhPhiBandConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
        fhPhiBandConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPhiBandConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeNTrackPerMinPtCut) ;
        
        fhPhiBandConeSumPtTrackPerMinPtCut = new TH2F
        ("hPhiBandConePtSumTrackPerMinPtCut",
         Form("Track #Sigma #it{p}_{T}-UE, different #it{p}_{T} cuts in #varphi-band for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPhiBandConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPhiBandConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhPhiBandConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeSumPtTrackPerMinPtCut) ;
        
        fhConeNTrackSubPhiBandPerMinPtCut = new TH2F
        ("hConeNTrackSubPhiBandPerMinPtCut",
         Form("N tracks - UE in #varphi-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
        fhConeNTrackSubPhiBandPerMinPtCut->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #varphi-band}");
        fhConeNTrackSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNTrackSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubPhiBandPerMinPtCut) ;
        
        fhConeSumPtTrackSubPhiBandPerMinPtCut = new TH2F
        ("hConePtSumTrackSubPhiBandPerMinPtCut",
         Form("Track cone #Sigma #it{p}_{T} - UE in #varphi-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
        fhConeSumPtTrackSubPhiBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
        fhConeSumPtTrackSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtTrackSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubPhiBandPerMinPtCut) ;
        
        fhConeSumPtTrackPerMaxPtCut = new TH2F
        ("hConePtSumTrackPerMaxPtCut",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMaxPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMaxPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMaxPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMaxPtCut) ;
      }      
      else
      {
        fhConeNTrackPerMinPtCutCent = new TH3F
        ("hConeNTrackPerMinPtCutCent",
         Form("Tracks, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           mulBinsArray.GetSize() - 1,   mulBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray());   
        fhConeNTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNTrackPerMinPtCutCent->SetYTitle("#it{N}^{track}");
        fhConeNTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackPerMinPtCutCent) ;
        
        fhConeSumPtTrackPerMinPtCutCent = new TH3F
        ("hConePtSumTrackPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T}, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,  sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtTrackPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMinPtCutCent) ;
        
        // Perpendicular cones
        
        fhPerpConeNTrackPerMinPtCutCent = new TH3F
        ("hPerpConeNTrackPerMinPtCutCent",
         Form("Tracks, different min #it{p}_{T} cut in #perp cone for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           mulBinsArray.GetSize() - 1,   mulBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhPerpConeNTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPerpConeNTrackPerMinPtCutCent->SetYTitle("#it{N}^{track}");
        fhPerpConeNTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPerpConeNTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPerpConeNTrackPerMinPtCutCent) ;
        
        fhPerpConeSumPtTrackPerMinPtCutCent = new TH3F
        ("hPerpConePtSumTrackPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T}, different min #it{p}_{T} cut in #perp cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,  sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhPerpConeSumPtTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPerpConeSumPtTrackPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPtTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPerpConeSumPtTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPerpConeSumPtTrackPerMinPtCutCent) ;
        
        fhConeNTrackSubPerpPerMinPtCutCent = new TH3F
        ("hConeNTrackSubPerpPerMinPtCutCent",
         Form("N Tracks - UE in #perp cones, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin*2,-multmax,multmax,GetNCentrBin(),0,100);
        fhConeNTrackSubPerpPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNTrackSubPerpPerMinPtCutCent->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #perp}");
        fhConeNTrackSubPerpPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackSubPerpPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubPerpPerMinPtCutCent) ;
        
        fhConeSumPtTrackSubPerpPerMinPtCutCent = new TH3F
        ("hConePtSumTrackSubPerpPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T} - UE in #perp cones, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtTrackSubPerpPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtTrackSubPerpPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #perp} (GeV/#it{c})");
        fhConeSumPtTrackSubPerpPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackSubPerpPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubPerpPerMinPtCutCent) ;
        
        // Eta band
        
        fhEtaBandConeNTrackPerMinPtCutCent = new TH3F
        ("hEtaBandConeNTrackPerMinPtCutCent",
         Form("Tracks, different min #it{p}_{T} cut in #eta band for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin,multmin,multmax,GetNCentrBin(),0,100);
        fhEtaBandConeNTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhEtaBandConeNTrackPerMinPtCutCent->SetYTitle("#it{N}^{track}");
        fhEtaBandConeNTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhEtaBandConeNTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeNTrackPerMinPtCutCent) ;
        
        fhEtaBandConeSumPtTrackPerMinPtCutCent = new TH3F
        ("hEtaBandConePtSumTrackPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T}, different min #it{p}_{T} cut in #eta band for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,   sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,    cenBinsArray.GetArray()); 
        fhEtaBandConeSumPtTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhEtaBandConeSumPtTrackPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhEtaBandConeSumPtTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhEtaBandConeSumPtTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhEtaBandConeSumPtTrackPerMinPtCutCent) ;
        
        fhConeNTrackSubEtaBandPerMinPtCutCent = new TH3F
        ("hConeNTrackSubEtaBandPerMinPtCutCent",
         Form("N Tracks - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin*2,-multmax,multmax,GetNCentrBin(),0,100);
        fhConeNTrackSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNTrackSubEtaBandPerMinPtCutCent->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #eta-band}");
        fhConeNTrackSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubEtaBandPerMinPtCutCent) ;
        
        fhConeSumPtTrackSubEtaBandPerMinPtCutCent = new TH3F
        ("hConePtSumTrackSubEtaBandPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T} - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtTrackSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtTrackSubEtaBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
        fhConeSumPtTrackSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubEtaBandPerMinPtCutCent) ;
        
        // Phi band
        
        fhPhiBandConeNTrackPerMinPtCutCent = new TH3F
        ("hPhiBandConeNTrackPerMinPtCutCent",
         Form("Tracks, different min #it{p}_{T} cut in #varphi band for %s",parTitleR.Data()),
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           mulBinsArray.GetSize() - 1,   mulBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray());
        fhPhiBandConeNTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPhiBandConeNTrackPerMinPtCutCent->SetYTitle("#it{N}^{track}");
        fhPhiBandConeNTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPhiBandConeNTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeNTrackPerMinPtCutCent) ;
        
        fhPhiBandConeSumPtTrackPerMinPtCutCent = new TH3F
        ("hPhiBandConePtSumTrackPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T}, different min #it{p}_{T} cut in #varphi band for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,  sum0BinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhPhiBandConeSumPtTrackPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhPhiBandConeSumPtTrackPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPhiBandConeSumPtTrackPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhPhiBandConeSumPtTrackPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhPhiBandConeSumPtTrackPerMinPtCutCent) ;
        
        fhConeNTrackSubPhiBandPerMinPtCutCent = new TH3F
        ("hConeNTrackSubPhiBandPerMinPtCutCent",
         Form("N Tracks - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin*2,-multmax,multmax,GetNCentrBin(),0,100);
        fhConeNTrackSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeNTrackSubPhiBandPerMinPtCutCent->SetYTitle("#it{N}^{track}-#it{N}^{track}_{UE #varphi-band}");
        fhConeNTrackSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackSubPhiBandPerMinPtCutCent) ;
        
        fhConeSumPtTrackSubPhiBandPerMinPtCutCent = new TH3F
        ("hConePtSumTrackSubPhiBandPerMinPtCutCent",
         Form("Track #Sigma #it{p}_{T} - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
         //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
         minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
        fhConeSumPtTrackSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
        fhConeSumPtTrackSubPhiBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
        fhConeSumPtTrackSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackSubPhiBandPerMinPtCutCent) ;
      }
    }
    
    if ( fStudyEtaCutInCone )
    {
      fhConeSumPtTrackPerEtaCut = new TH2F
      ("hConePtSumTrackPerEtaCut",
       Form("Track cone #Sigma #it{p}_{T}, different #eta cuts for %s",parTitleR.Data()),
       fNEtaCutsInCone,0.5,fNEtaCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackPerEtaCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackPerEtaCut->SetXTitle("#eta_{max}");
      for(Int_t i = 1; i <= fNEtaCutsInCone; i++)
        fhConeSumPtTrackPerEtaCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fEtaCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtTrackPerEtaCut) ;
    }
    
    if ( fStudyRCutInCone )
    {
      fhConeSumPtTrackPerRCut = new TH2F
      ("hConePtSumTrackPerRCut","Track #Sigma #it{p}_{T}, different #it{R} cuts",
       fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackPerRCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackPerRCut->SetXTitle("#it{R}");
      for(Int_t i = 1; i <= fNRCutsInCone; i++)
        fhConeSumPtTrackPerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtTrackPerRCut) ;
      
      fhPtTrackInConePerRCut = new TH2F
      ("hPtTrackInConePerRCut","Track #it{p}_{T}, different #it{R} cuts",
       fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtTrackInConePerRCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
      fhPtTrackInConePerRCut->SetXTitle("#it{R}");
      for(Int_t i = 1; i <= fNRCutsInCone; i++)
        fhPtTrackInConePerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
      outputContainer->Add(fhPtTrackInConePerRCut) ;
    }
    
    if ( fStudyNCellsCut )
    {
      fhConeSumPtTrackPerNCellCut = new TH2F
      ("hConePtSumTrackPerNCellCut","Track #Sigma #it{p}_{T}, different #it{N}_{cell} cuts",
       fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackPerNCellCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackPerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
      for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
        fhConeSumPtTrackPerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
      outputContainer->Add(fhConeSumPtTrackPerNCellCut) ;
      
      fhPtTrackInConePerNCellCut = new TH2F
      ("hPtTrackInConePerNCellCut","Track #it{p}_{T}, different #it{N}_{cell} cuts",
       fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtTrackInConePerNCellCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
      fhPtTrackInConePerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
      for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
        fhPtTrackInConePerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
      outputContainer->Add(fhPtTrackInConePerNCellCut) ;
      
      for(Int_t ishsh = 0; ishsh < 4; ishsh++)
      {
        fhPtTrackInConePerNCellPerSM[ishsh] = new TH3F
        (Form("hPtTrackInConePerNCellPerSM_ShSh%d",ishsh),
         Form("Track #it{p}_{T} vs #it{n}_{cell} vs SM , 8 < #it{p}_{T}^{trig} < 12 GeV, sh. shape bin %d",ishsh),
         200,0,20,fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax);
        fhPtTrackInConePerNCellPerSM[ishsh]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhPtTrackInConePerNCellPerSM[ishsh]->SetYTitle("SM number");
        fhPtTrackInConePerNCellPerSM[ishsh]->SetXTitle("#it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConePerNCellPerSM[ishsh]) ;
        
        fhConeSumPtTrackPerNCellPerSM[ishsh] = new TH3F
        (Form("hConeSumPtTrackPerNCellPerSM_ShSh%d",ishsh),
         Form("Track #Sigma #it{p}_{T} in cone vs #it{n}_{cell} vs SM , 8 < #it{p}_{T}^{trig} < 12 GeV, sh. shape bin %d",ishsh),
         200,0,50,fNModules,-0.5,fNModules-0.5,cellBins,cellMin,cellMax);
        fhConeSumPtTrackPerNCellPerSM[ishsh]->SetZTitle("#it{n}_{cells}^{w>0.01}");
        fhConeSumPtTrackPerNCellPerSM[ishsh]->SetYTitle("SM number");
        fhConeSumPtTrackPerNCellPerSM[ishsh]->SetXTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtTrackPerNCellPerSM[ishsh]) ;
      }
    }
    
    if ( fStudyExoticTrigger )
    { 
      fhPtTrackInConeExoTrigger  = new TH2F
      ("hPtTrackInConeExoTrigger",
       Form("Track cone #it{p}_{T} for %s, exotic trigger",parTitleR.Data()),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeExoTrigger) ;
      
      fhConeSumPtTrackPerExoCut = new TH2F
      ("hConePtSumTrackPerExoCut","Track #Sigma #it{p}_{T}, different exoticity cuts",
       fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackPerExoCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackPerExoCut->SetXTitle("exoticity");
      for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
        fhConeSumPtTrackPerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
      outputContainer->Add(fhConeSumPtTrackPerExoCut) ;
      
      fhPtTrackInConePerExoCut = new TH2F
      ("hPtTrackInConePerExoCut","Track #it{p}_{T}, different exoticity cuts",
       fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
      fhPtTrackInConePerExoCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
      fhPtTrackInConePerExoCut->SetXTitle("exoticity");
      for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
        fhPtTrackInConePerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
      outputContainer->Add(fhPtTrackInConePerExoCut) ;
      
      fhConeSumPtTrackExoTrigger  = new TH2F
      ("hConePtSumTrackExoTrigger",
       Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, exo trigger",r),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackExoTrigger) ;
    }
    
    if ( fStudyTracksInCone )
    {
      Int_t ntofbins = 1000;
      Int_t mintof = -500;
      Int_t maxtof =  500;
      
      fhTrackTOFInCone  = new TH2F 
      ("hTrackTOFInCone","TOF signal vs track #it{p}_{T}", 
       nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
      fhTrackTOFInCone->SetYTitle("TOF signal (ns)");
      fhTrackTOFInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhTrackTOFInCone);
      
      if ( fStudyExoticTrigger )
      {
        fhTrackTOFInConeExoTrigger  = new TH2F 
        ("hTrackTOFInConeExoTrigger","TOF signal vs track #it{p}_{T}, exoticity > 0.97", 
         nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
        fhTrackTOFInConeExoTrigger->SetYTitle("TOF signal (ns)");
        fhTrackTOFInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhTrackTOFInConeExoTrigger);
      }
      
      fhTrackTOFInConeBC0  = new TH2F 
      ("hTrackTOFInConeBC0","TOF signal vs track #it{p}_{T}, BC=0", 
       nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
      fhTrackTOFInConeBC0->SetYTitle("TOF signal (ns)");
      fhTrackTOFInConeBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhTrackTOFInConeBC0);
      
      fhPtTrackInConeVtxBC0  = new TH2F("hPtTrackInConeVtxBC0",
                                        Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, TOF from BC==0",r),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeVtxBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeVtxBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeVtxBC0) ;
      
      fhEtaPhiTrackInCone = new TH2F("hEtaPhiTrackInCone",
                                     Form("#eta vs #varphi of Tracks in cone for %s",parTitleR.Data()),
                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInCone->SetXTitle("#eta");
      fhEtaPhiTrackInCone->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInCone) ;
      
      fhEtaTrackInCone = new TH2F("hEtaTrackInCone",
                                  Form("#eta vs #it{p}_{T} of Tracks in cone for %s",parTitleR.Data()),
                                  nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInCone->SetYTitle("#eta");
      fhEtaTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInCone) ;
      
      fhPhiTrackInCone = new TH2F("hPhiTrackInCone",
                                  Form("#varphi vs #it{p}_{T} of Tracks in cone for %s",parTitleR.Data()),
                                  nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInCone->SetYTitle("#varphi (rad)");
      fhPhiTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInCone) ;
      
      //
      // Different track cuts:
      //
      // TOF info
      //
      fhConeSumPtTrackTOFBC0  = new TH2F("hConePtSumTrackTOFBC0",
                                         Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track TOF BC=0",r),
                                         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackTOFBC0->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackTOFBC0->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackTOFBC0) ;
      
      fhConeSumPtTrackTOFBCN  = new TH2F("hConePtSumTrackTOFBCN",
                                         Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track TOF BC!=0",r),
                                         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackTOFBCN->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackTOFBCN->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackTOFBCN) ;
      
      fhConeSumPtTrackTOFNo  = new TH2F("hConePtSumTrackTOFNo",
                                        Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track no TOF",r),
                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackTOFNo->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackTOFNo->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackTOFNo) ;
      
      fhPtTrackInConeTOFBC0  = new TH2F("hPtTrackInConeTOFBC0",
                                        Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, TOF from BC=0",r),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeTOFBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeTOFBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeTOFBC0) ;
      
      fhPtTrackInConeTOFBCN  = new TH2F("hPtTrackInConeTOFBCN",
                                        Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0",r),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeTOFBCN->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeTOFBCN->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeTOFBCN) ;
      
      fhPtTrackInConeTOFNo  = new TH2F("hPtTrackInConeTOFNo",
                                       Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, no TOF",r),
                                       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeTOFNo->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeTOFNo->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeTOFNo) ;
      
      
      fhEtaPhiTrackInConeTOFBC0 = new TH2F("hEtaPhiTrackInConeTOFBC0",
                                           Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, TOF BC=0",r),
                                           netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeTOFBC0->SetXTitle("#eta");
      fhEtaPhiTrackInConeTOFBC0->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeTOFBC0) ;  
      
      fhEtaPhiTrackInConeTOFBCN = new TH2F("hEtaPhiTrackInConeTOFBCN",
                                           Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, TOF BC!=0",r),
                                           netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeTOFBCN->SetXTitle("#eta");
      fhEtaPhiTrackInConeTOFBCN->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeTOFBCN) ;
      
      fhEtaPhiTrackInConeTOFNo = new TH2F("hEtaPhiTrackInConeTOFNo",
                                          Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, no TOF",r),
                                          netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeTOFNo->SetXTitle("#eta");
      fhEtaPhiTrackInConeTOFNo->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeTOFNo) ;
      
      fhEtaTrackInConeTOFBC0 = new TH2F("hEtaTrackInConeTOFBC0",
                                        Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC=0",r),
                                        nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeTOFBC0->SetYTitle("#eta");
      fhEtaTrackInConeTOFBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeTOFBC0) ;  
      
      fhEtaTrackInConeTOFBCN = new TH2F("hEtaTrackInConeTOFBCN",
                                        Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC!=0",r),
                                        nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeTOFBCN->SetYTitle("#eta");
      fhEtaTrackInConeTOFBCN->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeTOFBCN) ;
      
      fhEtaTrackInConeTOFNo = new TH2F("hEtaTrackInConeTOFNo",
                                       Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, no TOF",r),
                                       nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeTOFNo->SetYTitle("#eta");
      fhEtaTrackInConeTOFNo->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeTOFNo) ;
      
      fhPhiTrackInConeTOFBC0 = new TH2F("hPhiTrackInConeTOFBC0",
                                        Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC=0",r),
                                        nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeTOFBC0->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeTOFBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeTOFBC0) ;  
      
      fhPhiTrackInConeTOFBCN = new TH2F("hPhiTrackInConeTOFBCN",
                                        Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC!=0",r),
                                        nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeTOFBCN->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeTOFBCN->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeTOFBCN) ;
      
      fhPhiTrackInConeTOFNo = new TH2F("hPhiTrackInConeTOFNo",
                                       Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, no TOF",r),
                                       nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeTOFNo->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeTOFNo->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeTOFNo) ;
      
      //
      // ITS info
      //
      fhConeSumPtTrackITSRefitOnSPDOn  = new TH2F("hConePtSumTrackITSRefitOnSPDOn",
                                                  Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track ITS Refit SPD On",r),
                                                  nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackITSRefitOnSPDOn->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackITSRefitOnSPDOn->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackITSRefitOnSPDOn) ;
      
      fhConeSumPtTrackITSRefitOnSPDOff  = new TH2F("hConePtSumTrackITSRefitOnSPDOff",
                                                   Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track ITS Refit SPD Off",r),
                                                   nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackITSRefitOnSPDOff->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackITSRefitOnSPDOff->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackITSRefitOnSPDOff) ;
      
      fhConeSumPtTrackITSRefitOffSPDOff  = new TH2F("hConePtSumTrackITSRefitOffSPDOff",
                                                    Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track no ITS Refit SPD Off",r),
                                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackITSRefitOffSPDOff->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackITSRefitOffSPDOff->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackITSRefitOffSPDOff) ;
      
      fhPtTrackInConeITSRefitOnSPDOn  = new TH2F("hPtTrackInConeITSRefitOnSPDOn",
                                                 Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, TOF from BC=0",r),
                                                 nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeITSRefitOnSPDOn->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeITSRefitOnSPDOn) ;
      
      fhPtTrackInConeITSRefitOnSPDOff  = new TH2F("hPtTrackInConeITSRefitOnSPDOff",
                                                  Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0",r),
                                                  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeITSRefitOnSPDOff->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeITSRefitOnSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeITSRefitOnSPDOff) ;
      
      fhPtTrackInConeITSRefitOffSPDOff  = new TH2F("hPtTrackInConeITSRefitOffSPDOff",
                                                   Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, no ITS Refit SPD Off",r),
                                                   nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeITSRefitOffSPDOff->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeITSRefitOffSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeITSRefitOffSPDOff) ;
      
      
      fhEtaPhiTrackInConeITSRefitOnSPDOn = new TH2F("hEtaPhiTrackInConeITSRefitOnSPDOn",
                                                    Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD On",r),
                                                    netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeITSRefitOnSPDOn->SetXTitle("#eta");
      fhEtaPhiTrackInConeITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeITSRefitOnSPDOn) ;  
      
      fhEtaPhiTrackInConeITSRefitOnSPDOff = new TH2F("hEtaPhiTrackInConeITSRefitOnSPDOff",
                                                     Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD Off",r),
                                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeITSRefitOnSPDOff->SetXTitle("#eta");
      fhEtaPhiTrackInConeITSRefitOnSPDOff->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeITSRefitOnSPDOff) ;
      
      fhEtaPhiTrackInConeITSRefitOffSPDOff = new TH2F("hEtaPhiTrackInConeITSRefitOffSPDOff",
                                                      Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, no ITS Refit SPD Off",r),
                                                      netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeITSRefitOffSPDOff->SetXTitle("#eta");
      fhEtaPhiTrackInConeITSRefitOffSPDOff->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeITSRefitOffSPDOff) ;
      
      fhEtaTrackInConeITSRefitOnSPDOn = new TH2F("hEtaTrackInConeITSRefitOnSPDOn",
                                                 Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD On",r),
                                                 nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeITSRefitOnSPDOn->SetYTitle("#eta");
      fhEtaTrackInConeITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeITSRefitOnSPDOn) ;  
      
      fhEtaTrackInConeITSRefitOnSPDOff = new TH2F("hEtaTrackInConeITSRefitOnSPDOff",
                                                  Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD Off",r),
                                                  nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeITSRefitOnSPDOff->SetYTitle("#eta");
      fhEtaTrackInConeITSRefitOnSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeITSRefitOnSPDOff) ;
      
      fhEtaTrackInConeITSRefitOffSPDOff = new TH2F("hEtaTrackInConeITSRefitOffSPDOff",
                                                   Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, no ITS Refit SPD Off",r),
                                                   nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeITSRefitOffSPDOff->SetYTitle("#eta");
      fhEtaTrackInConeITSRefitOffSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeITSRefitOffSPDOff) ;
      
      fhPhiTrackInConeITSRefitOnSPDOn = new TH2F("hPhiTrackInConeITSRefitOnSPDOn",
                                                 Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD On",r),
                                                 nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeITSRefitOnSPDOn) ;  
      
      fhPhiTrackInConeITSRefitOnSPDOff = new TH2F("hPhiTrackInConeITSRefitOnSPDOff",
                                                  Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, ITS Refit SPD Off",r),
                                                  nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeITSRefitOnSPDOff->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeITSRefitOnSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeITSRefitOnSPDOff) ;
      
      fhPhiTrackInConeITSRefitOffSPDOff = new TH2F("hPhiTrackInConeITSRefitOffSPDOff",
                                                   Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, no ITS Refit SPD Off",r),
                                                   nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeITSRefitOffSPDOff->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeITSRefitOffSPDOff->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeITSRefitOffSPDOff) ;
      
      //
      // TOF and ITS info
      //
      fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn  = new TH2F("hConePtSumTrackTOFBC0ITSRefitOnSPDOn",
                                                        Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, track TOF BC=0, track ITS Refit SPD On",r),
                                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn) ;
      
      fhPtTrackInConeTOFBC0ITSRefitOnSPDOn  = new TH2F("hPtTrackInConeTOFBC0ITSRefitOnSPDOn",
                                                       Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, TOF from BC=0, track ITS Refit SPD On",r),
                                                       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeTOFBC0ITSRefitOnSPDOn) ;
      
      fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn = new TH2F("hEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn",
                                                          Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f, TOF BC=0, track ITS Refit SPD On",r),
                                                          netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#eta");
      fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn) ;  
      
      fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn = new TH2F("hEtaTrackInConeTOFBC0ITSRefitOnSPDOn",
                                                       Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC=0, track ITS Refit SPD On",r),
                                                       nptbins,ptmin,ptmax,netabins,-1,1);
      fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#eta");
      fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn) ;  
      
      fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn = new TH2F("hPhiTrackInConeTOFBC0ITSRefitOnSPDOn",
                                                       Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f, TOF BC=0, track ITS Refit SPD On",r),
                                                       nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn) ;  
      
      //
      // DCA
      //
      TString dcaName[] = {"xy","z","Cons"} ;
      Int_t ndcabins = 400;
      Int_t mindca = -2;
      Int_t maxdca =  2;
      
      for(Int_t i = 0 ; i < 3 ; i++)
      {
        fhPtTrackInConeDCA[i]  = new TH2F(Form("hPtTrackInConeDCA%s",dcaName[i].Data()),
                                          Form("Track DCA%s vs #it{p}_{T}^{track} in cone for trigger #it{p}_{T} >10 GeV/#it{c}",dcaName[i].Data()),
                                          nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtTrackInConeDCA[i]->SetXTitle("#it{p}_{T}^{} (GeV/#it{c})");
        fhPtTrackInConeDCA[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtTrackInConeDCA[i]);
        
        fhPtTrackInPerpConeDCA[i]  = new TH2F(Form("hPtTrackInPerpConeDCA%s",dcaName[i].Data()),
                                              Form("Track DCA%s vs #it{p}_{T}^{track} in perpendicular cone for trigger #it{p}_{T} >10 GeV/#it{c}",dcaName[i].Data()),
                                              nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtTrackInPerpConeDCA[i]->SetXTitle("#it{p}_{T}^{} (GeV/#it{c})");
        fhPtTrackInPerpConeDCA[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtTrackInPerpConeDCA[i]);
      }
      
      // TOF info
      fhPerpConeSumPtTOFBC0  = new TH2F("hPerpConePtSumTOFBC0",
                                        Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, TOF BC=0",r),
                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPerpConeSumPtTOFBC0->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhPerpConeSumPtTOFBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPerpConeSumPtTOFBC0) ;
      
      fhPtInPerpConeTOFBC0  = new TH2F("hPtInPerpConeTOFBC0",
                                       Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, TOF BC=0",r),
                                       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInPerpConeTOFBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInPerpConeTOFBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInPerpConeTOFBC0) ;
      
      fhEtaPhiInPerpConeTOFBC0= new TH2F("hEtaPhiInPerpConeTOFBC0",
                                         Form("#eta vs #varphi of all Tracks, TOF BC=0"),
                                         netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiInPerpConeTOFBC0->SetXTitle("#eta");
      fhEtaPhiInPerpConeTOFBC0->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiInPerpConeTOFBC0) ;
      
      // ITS info
      fhPerpConeSumPtITSRefitOnSPDOn  = new TH2F("hPerpConePtSumITSRefitOnSPDOn",
                                                 Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, ITS Refit, SPD On",r),
                                                 nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPerpConeSumPtITSRefitOnSPDOn->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhPerpConeSumPtITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPerpConeSumPtITSRefitOnSPDOn) ;
      
      fhPtInPerpConeITSRefitOnSPDOn  = new TH2F("hPtInPerpConeITSRefitOnSPDOn",
                                                Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, ITS Refit, SPD On",r),
                                                nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInPerpConeITSRefitOnSPDOn->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInPerpConeITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInPerpConeITSRefitOnSPDOn) ;
      
      fhEtaPhiInPerpConeITSRefitOnSPDOn= new TH2F("hEtaPhiInPerpConeITSRefitOnSPDOn",
                                                  Form("#eta vs #varphi of all Tracks, ITS Refit, SPD On"),
                                                  netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiInPerpConeITSRefitOnSPDOn->SetXTitle("#eta");
      fhEtaPhiInPerpConeITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiInPerpConeITSRefitOnSPDOn) ;
      
      
      // TOF and ITS info
      fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn  = new TH2F("hPerpConePtSumTOFBC0ITSRefitOnSPDOn",
                                                       Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, TOF BC=0, ITS refit, SPD on",r),
                                                       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn) ;
      
      fhPtInPerpConeTOFBC0ITSRefitOnSPDOn  = new TH2F("hPtInPerpConeTOFBC0ITSRefitOnSPDOn",
                                                      Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f, TOF BC=0, ITS refit, SPD on",r),
                                                      nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInPerpConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInPerpConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInPerpConeTOFBC0ITSRefitOnSPDOn) ;
      
      fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn = new TH2F("hEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn",
                                                         Form("#eta vs #varphi of all Tracks, TOF BC=0, ITS refit, SPD on"),
                                                         netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn->SetXTitle("#eta");
      fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn) ;
      
    }
  }
  
  
  if ( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged &&
      fStudyPtCutInCone )
  {
    if ( !IsHighMultiplicityAnalysisOn() )
    {
      // Eta band  
      
      fhConeNSubEtaBandPerMinPtCut = new TH2F
      ("hConeNSubEtaBandPerMinPtCut",
       Form("N - UE in #eta-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
       fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
      fhConeNSubEtaBandPerMinPtCut->SetYTitle("#it{N}-#it{N}_{UE #eta-band}");
      fhConeNSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
      fhConeNSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeNSubEtaBandPerMinPtCut) ;
      
      fhConeSumPtSubEtaBandPerMinPtCut = new TH2F
      ("hConePtSumSubEtaBandPerMinPtCut",
       Form("#Sigma #it{p}_{T} - UE in #eta-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
       fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
      fhConeSumPtSubEtaBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
      fhConeSumPtSubEtaBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
      fhConeSumPtSubEtaBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtSubEtaBandPerMinPtCut) ;
      
      // Phi band
      
      fhConeNSubPhiBandPerMinPtCut = new TH2F
      ("hConeNSubPhiBandPerMinPtCut",
       Form("N - UE in #varphi-band, different #it{p}_{T} cuts in isolation cone for %s",parTitleR.Data()),
       fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin*2,-multmax,multmax);
      fhConeNSubPhiBandPerMinPtCut->SetYTitle("#it{N}-#it{N}_{UE #varphi-band}");
      fhConeNSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
      fhConeNSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeNSubPhiBandPerMinPtCut) ;
      
      fhConeSumPtSubPhiBandPerMinPtCut = new TH2F
      ("hConePtSumSubPhiBandPerMinPtCut",
       Form(" cone #Sigma #it{p}_{T} - UE in #varphi-band, different #it{p}_{T} cut for %s",parTitleR.Data()),
       fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins*2,-ptsummax,ptsummax);
      fhConeSumPtSubPhiBandPerMinPtCut->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
      fhConeSumPtSubPhiBandPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
      fhConeSumPtSubPhiBandPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtSubPhiBandPerMinPtCut) ;
    }
    else 
    {
      // Eta band
      
      fhConeNSubEtaBandPerMinPtCutCent = new TH3F
      ("hConeNSubEtaBandPerMinPtCutCent",
       Form("N s - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
       fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin*2,-multmax,multmax,GetNCentrBin(),0,100);
      fhConeNSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
      fhConeNSubEtaBandPerMinPtCutCent->SetYTitle("#it{N}^-#it{N}_{UE #eta-band}");
      fhConeNSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeNSubEtaBandPerMinPtCutCent) ;
      
      fhConeSumPtSubEtaBandPerMinPtCutCent = new TH3F
      ("hConePtSumSubEtaBandPerMinPtCutCent",
       Form(" #Sigma #it{p}_{T} - UE in #eta band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
       //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
       minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
         cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
      fhConeSumPtSubEtaBandPerMinPtCutCent->SetZTitle("Centrality (%)");
      fhConeSumPtSubEtaBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #eta-band} (GeV/#it{c})");
      fhConeSumPtSubEtaBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtSubEtaBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtSubEtaBandPerMinPtCutCent) ;
      
      // Phi band
      
      fhConeNSubPhiBandPerMinPtCutCent = new TH3F
      ("hConeNSubPhiBandPerMinPtCutCent",
       Form("N - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
       fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin*2,-multmax,multmax,GetNCentrBin(),0,100);
      fhConeNSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
      fhConeNSubPhiBandPerMinPtCutCent->SetYTitle("#it{N}-#it{N}_{UE #varphi-band}");
      fhConeNSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeNSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeNSubPhiBandPerMinPtCutCent) ;
      
      fhConeSumPtSubPhiBandPerMinPtCutCent = new TH3F
      ("hConePtSumSubPhiBandPerMinPtCutCent",
       Form("#Sigma #it{p}_{T} - UE in #varphi band, different min #it{p}_{T} cut in cone for %s",parTitleR.Data()),
       //fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax,GetNCentrBin(),0,100);
       minPtBinsArray.GetSize() - 1, minPtBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1,   sueBinsArray.GetArray(),
         cenBinsArray.GetSize() - 1,   cenBinsArray.GetArray()); 
      fhConeSumPtSubPhiBandPerMinPtCutCent->SetZTitle("Centrality (%)");
      fhConeSumPtSubPhiBandPerMinPtCutCent->SetYTitle("#Sigma #it{p}_{T}-#Sigma #it{p}_{T, UE #varphi-band} (GeV/#it{c})");
      fhConeSumPtSubPhiBandPerMinPtCutCent->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
      for(Int_t i = 1; i <= fNPtCutsInCone; i++)
        fhConeSumPtSubPhiBandPerMinPtCutCent->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
      outputContainer->Add(fhConeSumPtSubPhiBandPerMinPtCutCent) ;
    }
    
  } // Neutral and charged and fStudyPtInCone
  
  TString region[] = {"ITS","TPC","TRD","TOF","Top EMCal","In EMCal"}; // conversion regions
  
  for(Int_t iso = 0; iso < 2; iso++)
  {
    if ( fFillTMHisto )
    {
      fhTrackMatchedDEta[iso]  = new TH2F
      (Form("hTrackMatchedDEta%s",isoName[iso].Data()),
       Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEta[iso]->SetYTitle("d#eta");
      fhTrackMatchedDEta[iso]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhi[iso]  = new TH2F
      (Form("hTrackMatchedDPhi%s",isoName[iso].Data()),
       Form("%s - d#varphi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhi[iso]->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDPhi[iso]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhi[iso]  = new TH2F
      (Form("hTrackMatchedDEtaDPhi%s",isoName[iso].Data()),
       Form("%s - d#eta vs d#varphi of cluster-track, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhi[iso]->SetYTitle("d#varphi (rad)");
      fhTrackMatchedDEtaDPhi[iso]->SetXTitle("d#eta");
      
      outputContainer->Add(fhTrackMatchedDEta[iso]) ;
      outputContainer->Add(fhTrackMatchedDPhi[iso]) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhi[iso]) ;
      
      fhdEdx[iso]  = new TH2F
      (Form("hdEdx%s",isoName[iso].Data()),
       Form("%s - Matched track <d#it{E}/d#it{x}> vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhdEdx[iso]->SetXTitle("#it{E} (GeV)");
      fhdEdx[iso]->SetYTitle("<d#it{E}/d#it{x}>");
      outputContainer->Add(fhdEdx[iso]);
      
      fhEOverP[iso]  = new TH2F
      (Form("hEOverP%s",isoName[iso].Data()),
       Form("%s - Matched track #it{E}/#it{p} vs cluster, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhEOverP[iso]->SetXTitle("#it{E} (GeV)");
      fhEOverP[iso]->SetYTitle("#it{E}/#it{p}");
      outputContainer->Add(fhEOverP[iso]);
      
      if ( IsDataMC() && !GetReader()->AreMCPromptPhotonsSelected() )
      {
        for(int imc = 0; imc < fNumberMCParticleCases; imc++)
        {
          if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
          {
            if ( mcPartName[imc].Contains("Prompt") ) continue;
            if ( mcPartName[imc].Contains("Frag"  ) ) continue;
            if ( mcPartName[imc].Contains("ISR"   ) ) continue;
          }
          
          fhTrackMatchedDEtaMC[imc][iso] = new TH2F
          (Form("hTrackMatchedDEta%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
           Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
           nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
          fhTrackMatchedDEtaMC[imc][iso]->SetYTitle("d#eta");
          fhTrackMatchedDEtaMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
          
          fhTrackMatchedDPhiMC[imc][iso] = new TH2F
          (Form("hTrackMatchedDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
           Form("%s - d#varphi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
           nptbins,ptmin,ptmax,nresetabins,resphimin,resphimax);
          fhTrackMatchedDPhiMC[imc][iso]->SetYTitle("d#varphi (rad)");
          fhTrackMatchedDPhiMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
          
          fhTrackMatchedDEtaDPhiMC[imc][iso]  = new TH2F
          (Form("hTrackMatchedDEtaDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
           Form("%s - d#eta vs d#varphi of cluster-track, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
           nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
          fhTrackMatchedDEtaDPhiMC[imc][iso]->SetYTitle("d#varphi (rad)");
          fhTrackMatchedDEtaDPhiMC[imc][iso]->SetXTitle("d#eta");
          
          outputContainer->Add(fhTrackMatchedDEtaMC[imc][iso]) ;
          outputContainer->Add(fhTrackMatchedDPhiMC[imc][iso]) ;
          outputContainer->Add(fhTrackMatchedDEtaDPhiMC[imc][iso]);
        }
        
        fhTrackMatchedMCParticle[iso]  = new TH2F
        (Form("hTrackMatchedMCParticle%s",isoName[iso].Data()),
         Form("%s - Origin of particle vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,8,0,8);
        fhTrackMatchedMCParticle[iso]->SetXTitle("#it{E} (GeV)");
        //fhTrackMatchedMCParticle[iso]->SetYTitle("Particle type");
        
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(1 ,"Photon");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(2 ,"Electron");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(4 ,"Rest");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
        fhTrackMatchedMCParticle[iso]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
        
        outputContainer->Add(fhTrackMatchedMCParticle[iso]);
      }
    }
    
    if ( fFillSSHisto )
    {
      fhPtLambda0[iso]  = new TH2F
      (Form("hPtLambda0%s",isoName[iso].Data()),
       Form("%s, %s",
            isoTitle[iso].Data(), parTitle[iso].Data()),
       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhPtLambda0[iso]->SetYTitle("#sigma_{long}^{2}");
      fhPtLambda0[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtLambda0[iso]) ;
      
      if ( IsHighMultiplicityAnalysisOn() )
      {
        fhPtLambda0Cent[iso]  = new TH3F
        (Form("hPtLambda0Cent%s",isoName[iso].Data()),
         Form("%s, %s",
              isoTitle[iso].Data(), parTitle[iso].Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
          ssBinsArray.GetSize() - 1,  ssBinsArray.GetArray(),      
         cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
        fhPtLambda0Cent[iso]->SetYTitle("#sigma_{long}^{2}");
        fhPtLambda0Cent[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtLambda0Cent[iso]->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPtLambda0Cent[iso]) ;
      }
      
      if ( fFillPerSMHistograms )
      {
        Int_t totalSM = fLastModule-fFirstModule+1;
        
        fhPtPerSM[iso] = new TH2F
        (Form("hPtPerSM_%s",isoName[iso].Data()),
         Form("%s candidate #it{p}_{T} and SM number, %s",isoTitle[iso].Data(), parTitle[iso].Data()),
         nptbins,ptmin,ptmax,totalSM,fFirstModule-0.5,fLastModule+0.5);
        fhPtPerSM[iso]->SetYTitle("SuperModule ");
        fhPtPerSM[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPerSM[iso]) ;
        
        for(Int_t ism = 0; ism < fNModules; ism++)
        {
          if ( ism < fFirstModule || ism > fLastModule ) continue;
          
          fhPtLambda0PerSM[iso][ism]  = new TH2F
          (Form("hPtLambda0_SM%d_%s",ism,isoName[iso].Data()),
           Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}, SM %d, %s",isoTitle[iso].Data(), ism, parTitle[iso].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0PerSM[iso][ism]->SetYTitle("#sigma_{long}^{2}");
          fhPtLambda0PerSM[iso][ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0PerSM[iso][ism]) ;
          
          if ( fStudyNCellsCut )
          {
            fhPtLambda0PerSMNCellCut[iso][ism]  = new TH2F
            (Form("hPtLambda0_NCellCut_SM%d_%s",ism,isoName[iso].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}, n_{cell}^{w>0} > 4, SM %d, %s",isoTitle[iso].Data(), ism, parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0PerSMNCellCut[iso][ism]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0PerSMNCellCut[iso][ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtLambda0PerSMNCellCut[iso][ism]) ;
            
            fhPtNCellPerSM[iso][ism]  = new TH2F
            (Form("hPtNCell_SM%d_%s",ism,isoName[iso].Data()),
             Form("%s cluster : #it{p}_{T} vs n_{cell}^{w>0}, SM %d, %s",isoTitle[iso].Data(), ism, parTitle[iso].Data()),
             nptbins,ptmin,ptmax,30,-0.5,29.5);
            fhPtNCellPerSM[iso][ism]->SetYTitle("n_{cell}^{w>0}");
            fhPtNCellPerSM[iso][ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtNCellPerSM[iso][ism]) ;
            
            fhPtNCellLowM02PerSM[iso][ism]  = new TH2F
            (Form("hPtNCell_LowM02_SM%d_%s",ism,isoName[iso].Data()),
             Form("%s cluster : #it{p}_{T} vs n_{cell}^{w>0}, 0.1<#sigma_{long}^{2}<0.3, SM %d, %s",isoTitle[iso].Data(), ism, parTitle[iso].Data()),
             nptbins,ptmin,ptmax,30,-0.5,29.5);
            fhPtNCellLowM02PerSM[iso][ism]->SetYTitle("n_{cell}^{w>0}");
            fhPtNCellLowM02PerSM[iso][ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtNCellLowM02PerSM[iso][ism]) ;
            
            fhPtNCellHighM02PerSM[iso][ism]  = new TH2F
            (Form("hPtNCell_HighM02_SM%d_%s",ism,isoName[iso].Data()),
             Form("%s cluster : #it{p}_{T} vs n_{cell}^{w>0}, 0.5<#sigma_{long}^{2}<2, SM %d, %s",isoTitle[iso].Data(), ism, parTitle[iso].Data()),
             nptbins,ptmin,ptmax,30,-0.5,29.5);
            fhPtNCellHighM02PerSM[iso][ism]->SetYTitle("n_{cell}^{w>0}");
            fhPtNCellHighM02PerSM[iso][ism]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtNCellHighM02PerSM[iso][ism]) ;
          }
        }
      }
      
      if ( fFillPerTCardIndexHistograms )
      {          
        fhPtPerTCardIndex[iso] = new TH2F
        (Form("hPtPerTCardIndex_%s",isoName[iso].Data()),
         Form("%s candidate #it{p}_{T} and super-module number, %s",isoTitle[iso].Data(), parTitle[iso].Data()),
         nptbins,ptmin,ptmax,16,-0.5,15.5);
        fhPtPerTCardIndex[iso]->SetYTitle("T-Card cell max index");
        fhPtPerTCardIndex[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPerTCardIndex[iso]) ;
        
        for(Int_t itc = 0; itc < 16; itc++)
        {            
          fhPtLambda0PerTCardIndex[iso][itc]  = new TH2F
          (Form("hPtLambda0_TC%d_%s",itc,isoName[iso].Data()),
           Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}, TC index %d, %s",isoTitle[iso].Data(), itc, parTitle[iso].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0PerTCardIndex[iso][itc]->SetYTitle("#sigma_{long}^{2}");
          fhPtLambda0PerTCardIndex[iso][itc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0PerTCardIndex[iso][itc]) ;
        }
      }
      
      if ( IsDataMC() && !GetReader()->AreMCPromptPhotonsSelected() )
      {
        for(Int_t imc = 0; imc < fNumberMCParticleCases; imc++)
        {
          if ( GetReader()->AreMCFragmentationPhotonsRejected() )
          {
            if ( mcPartName[imc].Contains("Frag"  ) ) continue;
            if ( mcPartName[imc].Contains("Prompt") ) continue;
            if ( mcPartName[imc].Contains("ISR"   ) ) continue;
          }

          fhPtLambda0MC[imc][iso]  = new TH2F
          (Form("hPtLambda0%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
           Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s",
                isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0MC[imc][iso]->SetYTitle("#sigma_{long}^{2}");
          fhPtLambda0MC[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
          outputContainer->Add( fhPtLambda0MC[imc][iso]) ;
          
          fhPtLambda0MCConv[imc][iso]  = new TH2F
          (Form("hPtLambda0%s_MC%sConv",isoName[iso].Data(),mcPartName[imc].Data()),
           Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, from conversion",
                isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0MCConv[imc][iso]->SetYTitle("#sigma_{long}^{2}");
          fhPtLambda0MCConv[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
          outputContainer->Add( fhPtLambda0MCConv[imc][iso]) ;
          
          if ( fStudyNCellsCut )
          {
            fhPtLambda0MCNCellCut[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%sNCellCut",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}, n_{cell}^{w>0.1} > 4: %s %s",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCNCellCut[imc][iso]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0MCNCellCut[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCNCellCut[imc][iso]) ;
          }
          
          if ( fFillOverlapHistograms )
          {
            fhPtLambda0MCWith1Overlap[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, 1 overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCWith1Overlap[imc][iso]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0MCWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCWith1Overlap[imc][iso]) ;
            
            fhPtLambda0MCConvWith1Overlap[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, from conversion, 1 overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCConvWith1Overlap[imc][iso]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0MCConvWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCConvWith1Overlap[imc][iso]) ;
            
            fhPtLambda0MCWithNoOverlap[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%s_NoOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, no overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCWithNoOverlap[imc][iso]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0MCWithNoOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCWithNoOverlap[imc][iso]) ;
            
            fhPtLambda0MCConvWithNoOverlap[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%sConv_NoOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, from conversion, no overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCConvWithNoOverlap[imc][iso]->SetYTitle("#sigma_{long}^{2}");
            fhPtLambda0MCConvWithNoOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCConvWithNoOverlap[imc][iso]) ;
            
            fhPtNOverlap[imc][iso]  = new TH2F
            (Form("hPtNOverlaps%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, 1 overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,nnovbin,novmin,novmax);
            fhPtNOverlap[imc][iso]->SetYTitle("#it{N} overlaps");
            fhPtNOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtNOverlap[imc][iso]) ;
            
            fhPtNOverlapConv[imc][iso]  = new TH2F
            (Form("hPtNOverlaps%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s cluster : #it{p}_{T} vs #sigma_{long}^{2}: %s %s, from conversion, 1 overlap",
                  isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle[iso].Data()),
             nptbins,ptmin,ptmax,nnovbin,novmin,novmax);
            fhPtNOverlapConv[imc][iso]->SetYTitle("#it{N} overlaps");
            fhPtNOverlapConv[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtNOverlapConv[imc][iso]) ;
          }
          
        }
      }
      
      if ( fIsoDetector==kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 )
      {
        fhPtLambda0TRD[iso]  = new TH2F
        (Form("hPtLambda0TRD%s",isoName[iso].Data()),
         Form("%s cluster: #it{p}_{T} vs #sigma_{long}^{2}, SM behind TRD, %s",
              isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0TRD[iso]->SetYTitle("#sigma_{long}^{2}");
        fhPtLambda0TRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLambda0TRD[iso]) ;
      }
      
      if ( fFillNLMHistograms )
      {
        fhPtLambda0LocMax1[iso]  = new TH2F
        (Form("hPtLambda0LocMax1%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}=1, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0LocMax1[iso]->SetYTitle("#sigma_{long}^{2}");
        fhPtLambda0LocMax1[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda0LocMax1[iso]) ;
        
        fhPtLambda1LocMax1[iso]  = new TH2F
        (Form("hPtLambda1LocMax1%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}=1, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda1LocMax1[iso]->SetYTitle("#sigma_{short}^{2}");
        fhPtLambda1LocMax1[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda1LocMax1[iso]) ;
        
        fhPtLambda0LocMax2[iso]  = new TH2F
        (Form("hPtLambda0LocMax2%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}=2, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0LocMax2[iso]->SetYTitle("#sigma_{long}^{2}");
        fhPtLambda0LocMax2[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda0LocMax2[iso]) ;
        
        fhPtLambda1LocMax2[iso]  = new TH2F
        (Form("hPtLambda1LocMax2%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}=2, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda1LocMax2[iso]->SetYTitle("#sigma_{short}^{2}");
        fhPtLambda1LocMax2[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda1LocMax2[iso]) ;
        
        fhPtLambda0LocMaxN[iso]  = new TH2F
        ( Form("hPtLambda0LocMaxN%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}>2, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0LocMaxN[iso]->SetYTitle("#sigma_{long}^{2}");
        fhPtLambda0LocMaxN[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda0LocMaxN[iso]) ;
        
        fhPtLambda1LocMaxN[iso]  = new TH2F
        (Form("hPtLambda1LocMaxN%s",isoName[iso].Data()),
         Form("%s cluster (#eta) pairs: #it{n}_{LM}>2, %s",isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda1LocMaxN[iso]->SetYTitle("#sigma_{short}^{2}");
        fhPtLambda1LocMaxN[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c}");
        outputContainer->Add(fhPtLambda1LocMaxN[iso]) ;
      } // NLM
    } // SS histo
    
    if ( GetCalorimeter() == kEMCAL && fFillEMCALRegionHistograms )
    {
      for(Int_t ieta = 0; ieta < 4; ieta++) 
      {  
        for(Int_t iphi = 0; iphi < 3; iphi++) 
        {
          for(Int_t ism = 0; ism < fNModules; ism++) 
          {
            if ( ism < fFirstModule || ism > fLastModule ) continue;
            
            fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism] = 
            new TH2F(Form("hLam0_%s_eta%d_phi%d_sm%d",isoName[iso].Data(),ieta,iphi,ism),
                     Form("%s, #it{p}_{T} vs #sigma_{long}^{2}, sm %d, region #eta %d, #varphi %d",
                          isoTitle[iso].Data(),ism,ieta,iphi),
                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]->SetYTitle("#sigma_{long}^{2}");
            fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]->SetXTitle("#it{p}_{T} (GeV)");
            outputContainer->Add(fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]) ;
            
            if(iso==0)
            {
              fhConeSumPtTrackEMCALRegionPerSM[ieta][iphi][ism]  = new TH2F
              (Form("hConePtSumTrack_eta%d_phi%d_sm%d",ieta,iphi,ism),
               Form("Track #Sigma #it{p}_{T}, sm %d, region #eta %d, #varphi %d",ism,ieta,iphi),
               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
              fhConeSumPtTrackEMCALRegionPerSM[ieta][iphi][ism]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
              fhConeSumPtTrackEMCALRegionPerSM[ieta][iphi][ism]->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
              outputContainer->Add(fhConeSumPtTrackEMCALRegionPerSM[ieta][iphi][ism]) ;
              
              fhConeSumPtClusterEMCALRegionPerSM[ieta][iphi][ism]  = new TH2F
              (Form("hConePtSumCluster_eta%d_phi%d_sm%d",ieta,iphi,ism),
               Form("Track #Sigma #it{p}_{T}, sm %d, region #eta %d, #varphi %d",ism,ieta,iphi),
               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
              fhConeSumPtClusterEMCALRegionPerSM[ieta][iphi][ism]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
              fhConeSumPtClusterEMCALRegionPerSM[ieta][iphi][ism]->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
              outputContainer->Add(fhConeSumPtClusterEMCALRegionPerSM[ieta][iphi][ism]) ;
            }
          } // ism
        } // iphi 
      } // ieta
      
      Float_t ptLimit[] = {2,3,4,5,6,8,10,12};
      for(Int_t ipt = 0; ipt < 7; ipt++)
      {
        fhEtaPhiLam0BinPtBin[iso][ipt]  = new TH2F
        (Form("hEtaPhiLam0BinPtBin%d%s",ipt,isoName[iso].Data()),
         Form("%s, #eta vs #varphi in #it{p}_{T}=[%2.1f,%2.1f] GeV/#it{c} and #lambda^{2}_{0}=[0.3,0.4]",
              isoTitle[iso].Data(),ptLimit[ipt],ptLimit[ipt+1]),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiLam0BinPtBin[iso][ipt]->SetYTitle("#varphi (rad)");
        fhEtaPhiLam0BinPtBin[iso][ipt]->SetXTitle("#eta");
        outputContainer->Add(fhEtaPhiLam0BinPtBin[iso][ipt]) ;
      }
    } // regions in EMCal
    
    if (IsDataMC() && fStudyMCConversionRadius && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      fhMCConversionVertex[iso] = new TH2F
      (Form("hMCPhotonConversionVertex%s",isoName[iso].Data()),
       Form("%s, cluster from converted photon, #it{p}_{T} vs vertex distance, %s",
            isoTitle[iso].Data(),parTitle[iso].Data()),
       nptbins,ptmin,ptmax,500,0,500);
      fhMCConversionVertex[iso]->SetYTitle("#it{R} (cm)");
      fhMCConversionVertex[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCConversionVertex[iso]) ;
      
      if ( GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 )
      {
        fhMCConversionVertexTRD[iso] = new TH2F
        (Form("hMCPhotonConversionVertexTRD%s",isoName[iso].Data()),
         Form("%s, cluster from converted photon, #it{p}_{T} vs vertex distance, %s, SM covered by TRD",
              isoTitle[iso].Data(),parTitle[iso].Data()),
         nptbins,ptmin,ptmax,500,0,500);
        fhMCConversionVertexTRD[iso]->SetYTitle("#it{R} (cm)");
        fhMCConversionVertexTRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCConversionVertexTRD[iso]) ;
      }
      
      for(Int_t iR = 0; iR < 6; iR++)
      {
        fhMCConversionLambda0Rcut[iR][iso] = new TH2F
        (Form("hMCPhotonConversionLambda0%s_R%d",isoName[iso].Data(),iR),
         Form("%s, cluster from converted photon, #it{p}_{T} vs #sigma_{long}^{2}, conversion in %s",
              isoTitle[iso].Data(),region[iR].Data()),
         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCConversionLambda0Rcut[iR][iso]->SetYTitle("#sigma_{long}^{2}");
        fhMCConversionLambda0Rcut[iR][iso]->SetXTitle("#it{p}_{T} (GeV)");
        outputContainer->Add(fhMCConversionLambda0Rcut[iR][iso]) ;
        
        if(GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0)
        {
          fhMCConversionLambda0RcutTRD[iR][iso] = new TH2F
          (Form("hMCPhotonConversionLambda0TRD%s_R%d",isoName[iso].Data(),iR),
           Form("%s, cluster from converted photon, #it{p}_{T} vs #sigma_{long}^{2}, conversion in %s, SM covered by TRD",
                isoTitle[iso].Data(),region[iR].Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCConversionLambda0RcutTRD[iR][iso]->SetYTitle("#sigma_{long}^{2}");
          fhMCConversionLambda0RcutTRD[iR][iso]->SetXTitle("#it{p}_{T} (GeV)");
          outputContainer->Add(fhMCConversionLambda0RcutTRD[iR][iso]) ;
        }
      }
    } // conversion radius
    
  } // control histograms for isolated and non isolated objects
  
  if ( IsPileUpAnalysisOn() )
  {
    if ( fStudyTracksInCone )
    {
      fhPtTrackInConeOtherBCPileUpSPD  = new TH2F
      ("hPtTrackInConeOtherBCPileUpSPD",
       Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0, pile-up from SPD",r),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeOtherBCPileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeOtherBCPileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeOtherBCPileUpSPD) ;
      
      fhPtTrackInConeBC0PileUpSPD  = new TH2F
      ("hPtTrackInConeBC0PileUpSPD",
       Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0, pile-up from SPD",r),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeBC0PileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeBC0PileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeBC0PileUpSPD) ;
    }
    
    for (Int_t i = 0; i < 7 ; i++)
    {
      fhPtInConePileUp[i]  = new TH2F
      (Form("hPtInConePileUp%s",pileUpName[i].Data()),
       Form("#it{p}_{T} in isolation cone for #it{R} =  %2.2f, from pile-up (%s)",r,pileUpName[i].Data()),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConePileUp[i]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConePileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInConePileUp[i]) ;
    }
  }
  
  // Generator level primary particles
  if ( IsDataMC() && IsGeneratedParticlesAnalysisOn()  )
  {
    // For histograms in arrays, index in the array, corresponding to any particle origin
    
    for(Int_t i = 0; i < fgkNmcPrimTypes; i++)
    {
      if ( GetReader()->AreMCPromptPhotonsSelected() && !ppname[i].Contains("Prompt") ) continue;
      
      if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
      {
        if ( ppname[i].Contains("Prompt") ) continue;
        if ( ppname[i].Contains("Frag"  ) ) continue;
        if ( ppname[i].Contains("ISR"   ) ) continue;
      }

      if ( !fFillOnlyTH3Histo )
      {
        fhPtPrimMC[i]  = new TH1F
        (Form("hPtPrim_MC%s",ppname[i].Data()),
         Form("primary photon  %s, %s",pptype[i].Data(),parTitle[1].Data()),
         nptbins,ptmin,ptmax);
        fhPtPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMC[i]) ;
        
        fhPtPrimMCiso[i]  = new TH1F
        (Form("hPtPrim_MCiso%s",ppname[i].Data()),
         Form("primary isolated photon %s, %s",pptype[i].Data(),parTitle[1].Data()),
         nptbins,ptmin,ptmax);
        fhPtPrimMCiso[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCiso[i]) ;
        
        fhEtaPrimMC[i]  = new TH2F
        (Form("hEtaPrim_MC%s",ppname[i].Data()),
         Form("primary photon %s : #eta vs #it{p}_{T}, %s",pptype[i].Data(),parTitle[1].Data()),
         nptbins,ptmin,ptmax,200,-2,2);
        fhEtaPrimMC[i]->SetYTitle("#eta");
        fhEtaPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaPrimMC[i]) ;
        
        fhPhiPrimMC[i]  = new TH2F
        (Form("hPhiPrim_MC%s",ppname[i].Data()),
         Form("primary photon %s : #varphi vs #it{p}_{T}, %s",pptype[i].Data(),parTitle[1].Data()),
         nptbins,ptmin,ptmax,200,0.,TMath::TwoPi());
        fhPhiPrimMC[i]->SetYTitle("#varphi (rad)");
        fhPhiPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPhiPrimMC[i]) ;
      }
      
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeSumPtPrimMC[i] = new TH2F
        (Form("hConeSumPtPrim_MC%s",ppname[i].Data()),
         Form("primary #gamma %s: %s",pptype[i].Data(),parTitleR.Data()),
         //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
           ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sum0BinsArray.GetSize() - 1,sum0BinsArray.GetArray());
        fhConeSumPtPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPrimMC[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPrimMC[i]) ;

        if ( particle == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtChargedPrimMC[i] = new TH2F
          (Form("hConeSumPtChargedPrim_MC%s",ppname[i].Data()),
           Form("primary #gamma %s: %s",pptype[i].Data(),parTitleRCh.Data()),
           //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
             ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
           sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray());
          fhConeSumPtChargedPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtChargedPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtChargedPrimMC[i]) ;
        }

        if (  method >= AliIsolationCut::kSumBkgSubIC )
        {
          fhConeSumPtUESubPrimMC[i] = new TH2F
          (Form("hConeSumPtUESubPrim_MC%s",ppname[i].Data()),
           Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleR.Data()),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
          fhConeSumPtUESubPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtUESubPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}-UE (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtUESubPrimMC[i]) ;

          if ( particle == AliIsolationCut::kNeutralAndCharged )
          {
            fhConeSumPtUESubChargedPrimMC[i] = new TH2F
            (Form("hConeSumPtUESubChargedPrim_MC%s",ppname[i].Data()),
             Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleRCh.Data()),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
            fhConeSumPtUESubChargedPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhConeSumPtUESubChargedPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{ch}-UE (GeV/#it{c})");
            outputContainer->Add(fhConeSumPtUESubChargedPrimMC[i]) ;
          }
        }

        fhConeSumPtNeutralChargedRatioPrimMC[i] = new TH2F
        (Form("hConeSumPtNeutralChargedRatioPrim_MC%s",ppname[i].Data()),
         Form("primary #gamma %s: %s",pptype[i].Data(),parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ratBinsArray.GetSize() - 1, ratBinsArray.GetArray());
        fhConeSumPtNeutralChargedRatioPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtNeutralChargedRatioPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{neutral} / #Sigma #it{p}_{T}^{charged}");
        outputContainer->Add(fhConeSumPtNeutralChargedRatioPrimMC[i]) ;

        fhConeSumPtPerpConeNeutralChargedRatioPrimMC[i] = new TH2F
        (Form("hConeSumPtPerpConeNeutralChargedRatioPrim_MC%s",ppname[i].Data()),
         Form("primary #gamma %s: %s, #perp cone",pptype[i].Data(),parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ratBinsArray.GetSize() - 1, ratBinsArray.GetArray());
        fhConeSumPtPerpConeNeutralChargedRatioPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPerpConeNeutralChargedRatioPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{neutral, #perp cone} / #Sigma #it{p}_{T}^{charged, #perp cone}");
        outputContainer->Add(fhConeSumPtPerpConeNeutralChargedRatioPrimMC[i]) ;

       if ( IsEmbedingAnalysisOn() )
       {
         fhConeSumPtPrimMCEmb[i] = new TH2F
         (Form("hConeSumPtPrimEmb_MC%s",ppname[i].Data()),
          Form("primary #gamma %s: %s",pptype[i].Data(),parTitleR.Data()),
          //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
          sum0BinsArray.GetSize() - 1,sum0BinsArray.GetArray());
         fhConeSumPtPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
         fhConeSumPtPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
         outputContainer->Add(fhConeSumPtPrimMCEmb[i]) ;

         if ( particle == AliIsolationCut::kNeutralAndCharged )
         {
           fhConeSumPtChargedPrimMCEmb[i] = new TH2F
           (Form("hConeSumPtChargedPrimEmb_MC%s",ppname[i].Data()),
            Form("primary #gamma %s: %s",pptype[i].Data(),parTitleRCh.Data()),
            //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
              ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
            sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray());
           fhConeSumPtChargedPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
           fhConeSumPtChargedPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
           outputContainer->Add(fhConeSumPtChargedPrimMCEmb[i]) ;
         }

         if (  method >= AliIsolationCut::kSumBkgSubIC )
         {
           fhConeSumPtUESubPrimMCEmb[i] = new TH2F
           (Form("hConeSumPtUESubPrimEmb_MC%s",ppname[i].Data()),
            Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleR.Data()),
             ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
            sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
           fhConeSumPtUESubPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
           fhConeSumPtUESubPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T}-UE (GeV/#it{c})");
           outputContainer->Add(fhConeSumPtUESubPrimMCEmb[i]) ;

           if ( particle == AliIsolationCut::kNeutralAndCharged )
           {
             fhConeSumPtUESubChargedPrimMCEmb[i] = new TH2F
             (Form("hConeSumPtUESubChargedPrimEmb_MC%s",ppname[i].Data()),
              Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleRCh.Data()),
               ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
              sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
             fhConeSumPtUESubChargedPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
             fhConeSumPtUESubChargedPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T}^{ch}-UE (GeV/#it{c})");
             outputContainer->Add(fhConeSumPtUESubChargedPrimMCEmb[i]) ;
           }
         }
       } // embedded

      } // low mult
      else
      {
        fhConeSumPtCenPrimMC[i] = new TH3F
        (Form("hConeSumPtCenPrim_MC%s",ppname[i].Data()),
         Form("primary #gamma  %s: %s",pptype[i].Data(),parTitleR.Data()),
         //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
           ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
         sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray(),
          cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
        fhConeSumPtCenPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCenPrimMC[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCenPrimMC[i]->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtCenPrimMC[i]) ;

        if ( particle == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtCenChargedPrimMC[i] = new TH3F
          (Form("hConeSumPtCenChargedPrim_MC%s",ppname[i].Data()),
           Form("primary #gamma  %s: %s",pptype[i].Data(),parTitleRCh.Data()),
           //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
             ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
           sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray(),
            cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
          fhConeSumPtCenChargedPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtCenChargedPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
          fhConeSumPtCenChargedPrimMC[i]->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtCenChargedPrimMC[i]) ;
        }

        if (  method >= AliIsolationCut::kSumBkgSubIC )
        {
          fhConeSumPtUESubCenPrimMC[i] = new TH3F
          (Form("hConeSumPtUESubCenPrim_MC%s",ppname[i].Data()),
           Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleR.Data()),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
          fhConeSumPtUESubCenPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtUESubCenPrimMC[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtUESubCenPrimMC[i]->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtUESubCenPrimMC[i]) ;

          if ( particle == AliIsolationCut::kNeutralAndCharged )
          {
            fhConeSumPtUESubCenChargedPrimMC[i] = new TH3F
            (Form("hConeSumPtUESubCenChargedPrim_MC%s",ppname[i].Data()),
             Form("primary #gamma  %s: %s, UE subtracted",pptype[i].Data(),parTitleRCh.Data()),
             //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
             cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
            fhConeSumPtUESubCenChargedPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhConeSumPtUESubCenChargedPrimMC[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
            fhConeSumPtUESubCenChargedPrimMC[i]->SetZTitle("Centrality (%)");
            outputContainer->Add(fhConeSumPtUESubCenChargedPrimMC[i]) ;
          }
        }

        if ( IsEmbedingAnalysisOn() )
        {
          fhConeSumPtCenPrimMCEmb[i] = new TH3F
          (Form("hConeSumPtCenPrimEmb_MC%s",ppname[i].Data()),
           Form("primary #gamma  %s: %s",pptype[i].Data(),parTitleR.Data()),
           //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
             ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
           sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray(),
            cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
          fhConeSumPtCenPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtCenPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtCenPrimMCEmb[i]->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtCenPrimMCEmb[i]) ;

          if ( particle == AliIsolationCut::kNeutralAndCharged )
          {
            fhConeSumPtCenChargedPrimMCEmb[i] = new TH3F
            (Form("hConeSumPtCenChargedPrimEmb_MC%s",ppname[i].Data()),
             Form("primary #gamma  %s: %s",pptype[i].Data(),parTitleRCh.Data()),
             //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
               ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
             sum0BinsArray.GetSize() - 1, sum0BinsArray.GetArray(),
              cenBinsArray.GetSize() - 1,  cenBinsArray.GetArray());
            fhConeSumPtCenChargedPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhConeSumPtCenChargedPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
            fhConeSumPtCenChargedPrimMCEmb[i]->SetZTitle("Centrality (%)");
            outputContainer->Add(fhConeSumPtCenChargedPrimMCEmb[i]) ;
          }

          if (  method >= AliIsolationCut::kSumBkgSubIC )
          {
            fhConeSumPtUESubCenPrimMCEmb[i] = new TH3F
            (Form("hConeSumPtUESubCenPrimEmb_MC%s",ppname[i].Data()),
             Form("primary #gamma %s: %s, UE subtracted",pptype[i].Data(),parTitleR.Data()),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
             cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
            fhConeSumPtUESubCenPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhConeSumPtUESubCenPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
            fhConeSumPtUESubCenPrimMCEmb[i]->SetZTitle("Centrality (%)");
            outputContainer->Add(fhConeSumPtUESubCenPrimMCEmb[i]) ;

            if ( particle == AliIsolationCut::kNeutralAndCharged )
            {
              fhConeSumPtUESubCenChargedPrimMCEmb[i] = new TH3F
              (Form("hConeSumPtUESubCenChargedPrimEmb_MC%s",ppname[i].Data()),
               Form("primary #gamma  %s: %s, UE subtracted",pptype[i].Data(),parTitleRCh.Data()),
               //nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
                ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
               cenBinsArray.GetSize() - 1, cenBinsArray.GetArray());
              fhConeSumPtUESubCenChargedPrimMCEmb[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              fhConeSumPtUESubCenChargedPrimMCEmb[i]->SetYTitle("#Sigma #it{p}_{T}^{ch} (GeV/#it{c})");
              fhConeSumPtUESubCenChargedPrimMCEmb[i]->SetZTitle("Centrality (%)");
              outputContainer->Add(fhConeSumPtUESubCenChargedPrimMCEmb[i]) ;
            }
          }
        } // embedded

      } // high mult
    }
    
    if ( fMakePrimaryPi0DecayStudy )
    {
      fhPtPrimMCPi0DecayPairAcceptInConeLowPt  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPt",
       Form("primary photon  %s : #it{p}_{T}, pair in cone, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPt) ;
      
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPt",
       Form("isolated primary photon %s, pair in cone : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt) ;
      
      fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPtNoOverlap",
       Form("primary photon  %s, no overlap, pair in cone : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap) ;
      
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPtNoOverlap",
       Form("isolated primary photon  %s, pair in cone,no overlap : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap) ;
      
      fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPtNoOverlapCaloE",
       Form("primary photon  %s, no overlap, pair in cone, E > calo min: #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE) ;
      
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPtNoOverlapCaloE",
       Form("isolated primary photon  %s, pair in cone,no overlap, E > calo min: #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE) ;
      
      
      fhPtPrimMCPi0DecayPairNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairNoOverlap",
       Form("primary photon  %s, no overlap: #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairNoOverlap) ;
      
      fhPtPrimMCPi0DecayIsoPairNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairNoOverlap",
       Form("isolated primary photon  %s, no overlap: #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairNoOverlap) ;
      
      fhPtPrimMCPi0DecayPairOutOfCone  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairOutOfCone",
       Form("primary photon %s : #it{p}_{T}, pair out of cone, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfCone) ;
      
      fhPtPrimMCPi0DecayIsoPairOutOfCone  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairOutOfCone",
       Form("isolated primary photon %s, pair out of cone : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfCone) ;
      
      fhPtPrimMCPi0DecayPairOutOfAcceptance  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairOutOfAcceptance",
       Form("primary photon %s : #it{p}_{T}, pair out of acceptance, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfAcceptance) ;
      
      fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonPi0DecayPairOutOfAcceptanceNoOverlap",
       Form("primary photon %s : #it{p}_{T}, pair out of acceptance, no overlap, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap) ;
      
      fhPtPrimMCPi0DecayIsoPairOutOfAcceptance  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairOutOfAcceptance",
       Form("isolated primary photon %s, pair out of acceptance : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfAcceptance) ;
      
      fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonPi0DecayPairOutOfAcceptanceNoOverlap",
       Form("isolated primary photon %s, pair out of acceptance, no overlap : #it{p}_{T}, %s",
            pptype[kmcPrimPi0Decay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap) ;
      
      fhPtPrimMCPi0Overlap  = new TH1F
      ("hPtPrim_MCPi0Overlap",
       Form("primary %s, overlap: #it{p}_{T}, %s",
            pptype[kmcPrimPi0].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0Overlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0Overlap) ;
      
      fhPtPrimMCPi0IsoOverlap  = new TH1F
      ("hPtPrim_MCisoPi0Overlap",
       Form("primary %s, overlap: #it{p}_{T}, %s",
            pptype[kmcPrimPi0].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCPi0IsoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCPi0IsoOverlap) ;
      
      
      fhPtPrimMCEtaDecayPairAcceptInConeLowPt  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPt",
       Form("primary photon  %s : #it{p}_{T}, pair in cone, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPt) ;
      
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPt",
       Form("isolated primary photon %s, pair in cone : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt) ;
      
      fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPtNoOverlap",
       Form("primary photon  %s, no overlap, pair in cone : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap) ;
      
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPtNoOverlap",
       Form("isolated primary photon  %s, pair in cone,no overlap : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap) ;
      
      fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPtNoOverlapCaloE",
       Form("primary photon  %s, no overlap, pair in cone, E > calo min: #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE) ;
      
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPtNoOverlapCaloE",
       Form("isolated primary photon  %s, pair in cone,no overlap, E > calo min: #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE) ;
      
      
      fhPtPrimMCEtaDecayPairNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairNoOverlap",
       Form("primary photon  %s, no overlap: #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairNoOverlap) ;
      
      fhPtPrimMCEtaDecayIsoPairNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairNoOverlap",
       Form("isolated primary photon  %s, no overlap: #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairNoOverlap) ;
      
      fhPtPrimMCEtaDecayPairOutOfCone  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairOutOfCone",
       Form("primary photon %s : #it{p}_{T}, pair out of cone, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfCone) ;
      
      fhPtPrimMCEtaDecayIsoPairOutOfCone  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairOutOfCone",
       Form("isolated primary photon %s, pair out of cone : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfCone) ;
      
      fhPtPrimMCEtaDecayPairOutOfAcceptance  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairOutOfAcceptance",
       Form("primary photon %s : #it{p}_{T}, pair out of acceptance, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfAcceptance) ;
      
      fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap  = new TH1F
      ("hPtPrim_MCPhotonEtaDecayPairOutOfAcceptanceNoOverlap",
       Form("primary photon %s : #it{p}_{T}, pair out of acceptance, no overlap, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap) ;
      
      fhPtPrimMCEtaDecayIsoPairOutOfAcceptance  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairOutOfAcceptance",
       Form("isolated primary photon %s, pair out of acceptance : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfAcceptance) ;
      
      fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap  = new TH1F
      ("hPtPrim_MCisoPhotonEtaDecayPairOutOfAcceptanceNoOverlap",
       Form("isolated primary photon %s, pair out of acceptance, no overlap : #it{p}_{T}, %s",
            pptype[kmcPrimEtaDecay].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap) ;
      
      fhPtPrimMCEtaOverlap  = new TH1F
      ("hPtPrim_MCEtaOverlap",
       Form("primary %s, overlap: #it{p}_{T}, %s",
            pptype[kmcPrimEta].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaOverlap) ;
      
      fhPtPrimMCEtaIsoOverlap  = new TH1F
      ("hPtPrim_MCisoEtaOverlap",
       Form("primary %s, overlap: #it{p}_{T}, %s",
            pptype[kmcPrimEta].Data(),parTitle[1].Data()),
       nptbins,ptmin,ptmax);
      fhPtPrimMCEtaIsoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCEtaIsoOverlap) ;
    }
    
  }//Histos with MC primary particles generator level
  
  if ( IsPileUpAnalysisOn() )
  {
    for (Int_t i = 0; i < 7 ; i++)
    {
      for(Int_t iso = 0; iso < 2; iso++)
      {
        fhPtPileUp[i][iso]  = new TH1F
        (Form("hPt%sPileUp%s",isoName[iso].Data(),pileUpName[i].Data()),
         Form("Number of %s particles vs #it{p}_{T}, %s, pile-up event by %s",
              isoTitle[iso].Data(),parTitle[iso].Data(),pileUpName[i].Data()),
         nptbins,ptmin,ptmax);
        fhPtPileUp[i][iso]->SetYTitle("d#it{N} / #it{p}_{T}");
        fhPtPileUp[i][iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPileUp[i][iso]) ;
      }
    }
    
    fhTimeENoCut  = new TH2F 
    ("hTimeE_NoCut","time of cluster vs E of clusters, no cut", 
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeENoCut->SetXTitle("#it{E} (GeV)");
    fhTimeENoCut->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeENoCut);
    
    fhTimeESPD  = new TH2F 
    ("hTimeE_SPD","time of cluster vs E of clusters, SPD cut", 
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPD->SetXTitle("#it{E} (GeV)");
    fhTimeESPD->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPD);
    
    fhTimeESPDMulti  = new TH2F 
    ("hTimeE_SPDMulti","time of cluster vs E of clusters, SPD multi cut", 
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPDMulti->SetXTitle("#it{E} (GeV)");
    fhTimeESPDMulti->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPDMulti);
    
    fhTimeNPileUpVertSPD  = new TH2F 
    ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", 
     ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F 
    ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex",
     ntimebins,timemin,timemax, 50,0,50 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhTimeNPileUpVertContributors  = new TH2F 
    ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", 
     ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);
    
    fhTimePileUpMainVertexZDistance  = new TH2F 
    ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",
     ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDistance->SetYTitle("distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);
    
    fhTimePileUpMainVertexZDiamond  = new TH2F 
    ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",
     ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);
  }
  
  // Get histograms created and filled in AliIsolationCut and put them in the output list
  //
  if ( fFillIsolationControlHistograms )
  {
    GetIsolationCut()->SetHistogramRanges(GetHistogramRanges()); // important, not initialized inside AliIsolationCut
    GetIsolationCut()->SetNCentrBins(GetNCentrBin()); 
    TList * isoHistos = GetIsolationCut()->GetCreateOutputObjects() ;
    
    for(Int_t iIso = 0; iIso < isoHistos->GetEntries(); iIso++) 
    { 
      TH1 * obj = (TH1*) isoHistos->At(iIso);
      obj->SetName(Form("IsoCut_%s",obj->GetName()));
      outputContainer->Add(obj) ;
    }
    
    delete isoHistos;
  }
  
  // Return full list
  //
  return outputContainer ;
}

//____________________________________________________
/// MC histogram index depending on origin of candidate.
//____________________________________________________
Int_t AliAnaParticleIsolation::GetMCIndex(Int_t mcTag)
{
  if ( !IsDataMC() ) return -1;
  
  if     (GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPrompt))
  {
    return kmcPrompt;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCFragmentation) ||
          GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCISR))
  {
    return kmcFragment;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0))
  {
    return kmcPi0;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta))
  {
    return kmcEta;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0Decay))
  {
    return kmcPi0Decay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEtaDecay))
  {
    return kmcEtaDecay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
  {
    return kmcOtherDecay;
  }
  else if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))
  {
    return kmcElectron;
  }
  else // anything else
  {
    // careful can contain also other decays, to be checked.
    return kmcHadron;
  }
}

//__________________________________
/// Do some checks and init stuff.
//__________________________________
void AliAnaParticleIsolation::Init()
{
  if(!GetReader()->IsCTSSwitchedOn() && GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    AliFatal("STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!");
}

//____________________________________________
/// Initialize the parameters of the analysis
/// with default values.
//____________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  SetInputAODName("CaloTrackParticle");
  SetAODObjArrayName("IsolationCone");
  AddToHistogramsName("AnaIsolation_");
  
  fIsoDetectorString = "EMCAL" ;
  fIsoDetector       = kEMCAL  ;
    
  fMinCellsAngleOverlap = 3.;
  
  fLeadingOnly = kTRUE;
  fCheckLeadingWithNeutralClusters = kTRUE;
  
  fNDecayBits = 1;
  fDecayBits[0] = AliNeutralMesonSelection::kPi0;
  fDecayBits[1] = AliNeutralMesonSelection::kEta;
  fDecayBits[2] = AliNeutralMesonSelection::kPi0RightSide;
  fDecayBits[3] = AliNeutralMesonSelection::kEtaRightSide;
  fDecayBits[4] = AliNeutralMesonSelection::kEtaLeftSide;
  fDecayBits[5] = AliNeutralMesonSelection::kEtaBothSides;
  fDecayBits[6] = AliNeutralMesonSelection::kPi0LeftSide ; // Leave it last since likely not used
  fDecayBits[7] = AliNeutralMesonSelection::kPi0BothSides; // Leave it last since likely not used

  fM02Narrow[0] = 0.1;
  fM02Narrow[1] = 0.3;
  
  fM02Wide  [0] = 0.4;
  fM02Wide  [1] = 2.0;
  
  fNPtCutsInCone = 20;
  for(Int_t i = 0; i < 20; i++) fMinPtCutInCone[i] = 0.2+i*0.1;
  for(Int_t i = 0; i < 20; i++) fMaxPtCutInCone[i] = 3+i*1;

  fNEtaCutsInCone = 10;
  for(Int_t i = 0; i < 10; i++) fEtaCutInCone[i] = 0.8-i*0.05;

  fNRCutsInCone = 10;
  for(Int_t i = 0; i < 10; i++) fRCutInCone[i] = (i+1)*0.05;

  fNNCellsInCandidate = 20;
  for(Int_t i = 0; i < 20; i++) fNCellsInCandidate[i] = i+1;

  fNExoCutInCandidate = 20;
  for(Int_t i = 0; i < 20; i++) fExoCutInCandidate[i] = 0.61+0.02*i;
  
  fFillIsolationControlHistograms = kTRUE;
}

//_________________________________________________________________________________________
/// Check if the what of the selected isolation candidates is leading particle
/// in the same hemisphere comparing with all the candidates, all the tracks or all the clusters.
//_________________________________________________________________________________________
Bool_t AliAnaParticleIsolation::IsTriggerTheNearSideEventLeadingParticle(Int_t & idLeading)
{
  Double_t ptTrig      = GetMinPt();
  Double_t phiTrig     = 0 ;
  Int_t index          =-1 ;
  AliCaloTrackParticleCorrelation* pLeading = 0;
  
  // Loop on stored AOD particles, find leading trigger on the selected list, with at least min pT.
  
  for(Int_t iaod = 0; iaod < GetInputAODBranch()->GetEntriesFast() ; iaod++)
  {
    AliCaloTrackParticleCorrelation* particle =  (AliCaloTrackParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    particle->SetLeadingParticle(kFALSE); // set it later
    
    // Vertex cut in case of mixing
    if(GetMixedEvent())
    {
      Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
      if(check ==  0) continue;
      if(check == -1) return kFALSE; // not sure if it is correct.
    }
    
    //check if it is low pt trigger particle
    if ( particle->Pt() < GetIsolationCut()->GetPtThreshold() ||
         particle->Pt() < GetIsolationCut()->GetSumPtThreshold() )
    {
      continue ; //trigger should not come from underlying event
    }

    // find the leading particles with highest momentum
    if (particle->Pt() > ptTrig)
    {
      ptTrig   = particle->Pt() ;
      phiTrig  = particle->Phi();
      index    = iaod     ;
      pLeading = particle ;
    }
  }// finish search of leading trigger particle on the AOD branch.
  
  if(index < 0) return kFALSE;
  
  idLeading = index;
  
  //printf("AOD leading pT %2.2f, ID %d\n", pLeading->Pt(),pLeading->GetCaloLabel(0));
  
  if(phiTrig < 0 ) phiTrig += TMath::TwoPi();
  
  // Compare if it is the leading of all tracks
  
  for(Int_t ipr = 0;ipr < GetCTSTracks()->GetEntriesFast() ; ipr ++ )
  {
    AliVTrack * track = (AliVTrack *) (GetCTSTracks()->At(ipr)) ;
    
    // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
    // do not count the candidate or the daughters of the candidate
    // in the isolation conte
    if ( pLeading->GetDetectorTag() == AliFiducialCut::kCTS ) // make sure conversions are tagged as kCTS!!!
    {
      Int_t  trackID   = GetReader()->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
      Bool_t contained = kFALSE;
      
      for(Int_t i = 0; i < 4; i++) 
      {
        if( trackID == pLeading->GetTrackLabel(i) ) contained = kTRUE;
      }
      
      if ( contained ) continue ;
    }
    
    fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
    Float_t pt   = fTrackVector.Pt();
    Float_t phi  = fTrackVector.Phi() ;
    if(phi < 0) phi+=TMath::TwoPi();
    
    //skip this event if near side associated particle pt larger than trigger
    
    Float_t deltaPhi = phiTrig-phi;
    //
    // Calculate deltaPhi shift so that for the particles on the opposite side
    // it is defined between 90 and 270 degrees
    // Shift [-360,-90]  to [0, 270]
    // and [270,360] to [-90,0]
    if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();

    if(pt > ptTrig && deltaPhi < TMath::PiOver2())  return kFALSE;
  
  }// track loop
  
  // Compare if it is leading of all calorimeter clusters
  
  if(fCheckLeadingWithNeutralClusters)
  {
    // Select the calorimeter cluster list
    TObjArray * nePl = 0x0;
    if      (pLeading->GetDetectorTag() == kPHOS )
      nePl = GetPHOSClusters();
    else
      nePl = GetEMCALClusters();
    
    if(!nePl) return kTRUE; // Do the selection just with the tracks if no calorimeter is available.
    
    for(Int_t ipr = 0;ipr < nePl->GetEntriesFast() ; ipr ++ )
    {
      AliVCluster * cluster = (AliVCluster *) (nePl->At(ipr)) ;
      
      if(cluster->GetID() == pLeading->GetCaloLabel(0) || cluster->GetID() == pLeading->GetCaloLabel(1) ) continue ;
      
      cluster->GetMomentum(fMomentum,GetVertex(0));
      
      Float_t pt   = fMomentum.Pt();
      Float_t phi  = fMomentum.Phi() ;
      if(phi < 0) phi+=TMath::TwoPi();
      
      if(IsTrackMatched(cluster,GetReader()->GetInputEvent())) continue ; // avoid charged clusters, already covered by tracks, or cluster merging with track.
      
      // skip this event if near side associated particle pt larger than trigger
      // not really needed for calorimeter, unless DCal is included
     
      Float_t deltaPhi = phiTrig-phi;
      if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
      if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();

      if(pt > ptTrig && deltaPhi < TMath::PiOver2()) return kFALSE ;

    }// cluster loop
  } // check neutral clusters
  
  idLeading = index ;
  pLeading->SetLeadingParticle(kTRUE);
  
  AliDebug(1,Form("Particle AOD with index %d is leading with pT %2.2f",idLeading, pLeading->Pt()));
  
  return kTRUE;
}

//__________________________________________________
/// Do analysis and fill aods.
/// Search for the isolated photon in GetCalorimeter() with GetMinPt() < pt < GetMaxPt()
/// and if the particle is leading in the near side (if requested).
//__________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillAOD()
{
  if ( !GetInputAODBranch() )
    AliFatal(Form("No input particles in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
  
  if ( strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliCaloTrackParticleCorrelation") )
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliCaloTrackParticleCorrelation> and not <%s>",
                  GetInputAODBranch()->GetClass()->GetName()));
  
  Int_t n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0, coneptlead = 0;
  
  //Loop on AOD branch, filled previously in AliAnaPhoton, find leading particle to do isolation only with it
  Int_t idLeading = -1 ;
  Int_t iaod0 = 0;
  Int_t naod  = GetInputAODBranch()->GetEntriesFast();
  
  AliDebug(1,Form("Input aod branch entries %d", naod));
  
  if ( IsLeadingOnlyOn() )
  {
    Bool_t leading = IsTriggerTheNearSideEventLeadingParticle(idLeading);
    if ( !leading )
    {
      AliDebug(1,"Not leading; End fill AODs");
      return;
    }
    iaod0 = idLeading  ; // first entry in particle loop
    naod  = idLeading+1; // last entry in particle loop
  }
  
  // Check isolation of list of candidate particles or leading particle
  
  for(Int_t iaod = iaod0; iaod < naod; iaod++)
  {
    AliCaloTrackParticleCorrelation * aodinput =  (AliCaloTrackParticleCorrelation*) (GetInputAODBranch()->At(iaod));

    // Check isolation only of clusters in fiducial region
    
    if ( IsFiducialCutOn() )
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(aodinput->Eta(), aodinput->Phi(), aodinput->GetDetectorTag()) ;
      if(! in ) continue ;
    }
    
    // If too small or too large pt, skip
    Float_t pt = aodinput->Pt();
    if ( pt < GetMinPt() || pt > GetMaxPt() ) continue ;

    // Check if it is low pt trigger particle
    if ( pt < GetIsolationCut()->GetPtThreshold() ||  pt < GetIsolationCut()->GetSumPtThreshold() )
    {
      continue ; //trigger should not come from underlying event
    }
    
    //After cuts, study isolation
    n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0; coneptlead = 0;
    GetIsolationCut()->MakeIsolationCut(aodinput, GetReader(),
                                        kTRUE, kFALSE, GetAODObjArrayName(), 
                                        0x0, 0x0,
                                        GetCalorimeter(), GetCaloPID(),
                                        n, nfrac, coneptsum, coneptlead, isolated,
                                        GetEventWeight()*aodinput->GetWeight(),
                                        GetEventCentrality(), GetEventCentralityBin());

    aodinput->SetIsolated(isolated);
    
    AliDebug(1,Form("Particle isolated? %i; if so with index %d",isolated,iaod));
  } // particle isolation loop
}

//_________________________________________________________
/// Do analysis and fill histograms.
//_________________________________________________________
void  AliAnaParticleIsolation::MakeAnalysisFillHistograms()
{
  // In case of simulated data, fill acceptance histograms
  // Not ready for multiple case analysis.

  if ( IsDataMC() && IsGeneratedParticlesAnalysisOn() )
    FillAcceptanceHistograms();
  
  //Loop on stored AOD
  Int_t naod = GetInputAODBranch()->GetEntriesFast();

  AliDebug(1,Form("Histo aod branch entries %d", naod));
  
  Int_t   method = GetIsolationCut()->GetICMethod() ;
  Float_t sumThr = GetIsolationCut()->GetSumPtThreshold() ;
  Float_t sumGap = GetIsolationCut()->GetSumPtThresholdGap() ;
  Int_t partInCone = GetIsolationCut()->GetParticleTypeInCone() ;
  Int_t   icent  = GetEventCentralityBin();

  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliCaloTrackParticleCorrelation* aod =  (AliCaloTrackParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if ( IsLeadingOnlyOn() && !aod->IsLeadingParticle() ) continue; // Try to isolate only leading cluster or track
    
    // Check isolation only of clusters in fiducial region
    if ( IsFiducialCutOn() )
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(aod->Eta(),aod->Phi(),aod->GetDetectorTag()) ;
      if(! in ) continue ;
    }
    
    // If too small or too large pt, skip
    //
    Float_t pt         = aod->Pt();
    if ( pt < GetMinPt() || pt > GetMaxPt() ) continue ;
    
    Int_t mcTag        = aod->GetTag() ;
    Int_t mcIndex      = GetMCIndex(mcTag);
    Float_t weightTrig = aod->GetWeight();
    
    // Get isolation candidate kine and decision parameters
    //
    Bool_t  isolated = aod->IsIsolated();
    Float_t energy   = aod->E();
    Float_t phi      = GetPhi(aod->Phi());
    Float_t eta      = aod->Eta();
    Float_t m02      = aod->GetM02();
    Int_t   iSM      = aod->GetSModNumber();
    
    AliDebug(1,Form("pt %1.1f, eta %1.1f, phi %1.1f, Isolated %d, weight %2.3f",
                    pt, eta, phi, isolated, weightTrig));

    Float_t coneptLeadCluster = aod->GetNeutralLeadPtInCone();
    Float_t coneptsumCluster  = aod->GetNeutralPtSumInCone();
    Float_t coneptLeadTrack   = aod->GetChargedLeadPtInCone();
    Float_t coneptsumTrack    = aod->GetChargedPtSumInCone();
    Float_t coneptsum         = coneptsumTrack + coneptsumCluster;

    if ( partInCone == AliIsolationCut::kOnlyCharged )
      coneptsum  = coneptsumTrack;
    if ( partInCone == AliIsolationCut::kOnlyNeutral)
      coneptsum  = coneptsumCluster;
    
    Float_t coneptLead = coneptLeadTrack;
    if ( partInCone == AliIsolationCut::kNeutralAndCharged )
    {
      if ( coneptLeadCluster > coneptLeadTrack )
        coneptLead = coneptLeadCluster;
    }
    else if ( partInCone == AliIsolationCut::kOnlyNeutral )
    {
      coneptLead = coneptLeadCluster;
    }
    
    AliDebug(1,Form("Particle %d Energy Sum in Isolation Cone %2.2f, Leading pT in cone %2.2f",
                    iaod, coneptsum, coneptLead));
     
    Bool_t narrow       = kFALSE;
    Bool_t inM02Windows = kTRUE;
    
    if      ( m02 > fM02Narrow[0] && m02 < fM02Narrow[1] ) narrow = kTRUE;
    else if ( m02 > fM02Wide  [0] && m02 < fM02Wide  [1] ) narrow = kFALSE; 
    else inM02Windows = kFALSE; // skip clusters out of both ranges

    if ( method == AliIsolationCut::kSumPtIC || 
         method >= AliIsolationCut::kSumBkgSubIC )
    {
      if ( !IsHighMultiplicityAnalysisOn() )
        fhPtM02SumPtCone           ->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
      else if ( icent >= 0 && GetNCentrBin() > 0 && icent < GetNCentrBin() )
        fhPtM02SumPtConeCent[icent]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
      
      if ( partInCone == AliIsolationCut::kNeutralAndCharged )
      {
        if ( !IsHighMultiplicityAnalysisOn() )
          fhPtM02SumPtConeCharged           ->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
        else if ( icent >= 0 && GetNCentrBin() > 0 && icent < GetNCentrBin() )
          fhPtM02SumPtConeChargedCent[icent]->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
      }

      if ( inM02Windows && !fFillOnlyTH3Histo )
      {
        if ( !IsHighMultiplicityAnalysisOn() )
          fhConeSumPtM02Cut    [narrow]->Fill(pt, coneptsum, GetEventWeight()*weightTrig);
        else 
          fhConeSumPtCentM02Cut[narrow]->Fill(pt, coneptsum, GetEventCentrality(), GetEventWeight()*weightTrig);
      }

      if ( IsDataMC() && mcIndex < fNumberMCParticleCases && !GetReader()->AreMCPromptPhotonsSelected() )
      {
        if ( !IsHighMultiplicityAnalysisOn() )
        {
          fhPtM02SumPtConeMC[mcIndex]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtM02SumPtConeMC[kmcPhoton]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
          
          if ( partInCone == AliIsolationCut::kNeutralAndCharged )
          {
            fhPtM02SumPtConeChargedMC[mcIndex]->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
            
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtM02SumPtConeChargedMC[kmcPhoton]->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
          }
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) && 
              fNumberMCParticleCases > kmcPi0DecayLostPair )
          {
            if      ( mcIndex == kmcPi0Decay )
              fhPtM02SumPtConeMC[kmcPi0DecayLostPair]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
            else if ( mcIndex == kmcEtaDecay )
              fhPtM02SumPtConeMC[kmcEtaDecayLostPair]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
          }

          if ( inM02Windows  && !fFillOnlyTH3Histo )
          {
            fhConeSumPtM02CutMC[mcIndex][narrow]->Fill(pt, coneptsum, GetEventWeight()*weightTrig);

            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhConeSumPtM02CutMC[kmcPhoton][narrow]->Fill(pt, coneptsum, GetEventWeight()*weightTrig);
            
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) && 
                fNumberMCParticleCases > kmcPi0DecayLostPair )
            {
              if      ( mcIndex == kmcPi0Decay )
                fhConeSumPtM02CutMC[kmcPi0DecayLostPair][narrow]->Fill(pt, coneptsum, GetEventWeight()*weightTrig);
              else if ( mcIndex == kmcEtaDecay )
                fhConeSumPtM02CutMC[kmcEtaDecayLostPair][narrow]->Fill(pt, coneptsum, GetEventWeight()*weightTrig);
            }
          } // in m02 window          
        } // MC, not high mult
        else if ( icent < GetNCentrBin() && icent >=0  && GetNCentrBin() > 0 )
        {
          Int_t index   = icent*fNumberMCParticleCases+mcIndex;
          Int_t indexPh = icent*fNumberMCParticleCases+kmcPhoton;
          fhPtM02SumPtConeCentMC[index]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
          
          if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtM02SumPtConeCentMC[indexPh]->Fill(pt, m02, coneptsum, GetEventWeight()*weightTrig);
          
          if ( partInCone == AliIsolationCut::kNeutralAndCharged )
          {
            fhPtM02SumPtConeChargedCentMC[index]->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
            
            if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtM02SumPtConeChargedCentMC[indexPh]->Fill(pt, m02, coneptsumTrack, GetEventWeight()*weightTrig);
          }
        } // MC, high mult
      }
    } // TH3 histo and TH2 pT iso

    //---------------------------------------------------------------
    // Recover original cluster if requested, needed for some studies
    //---------------------------------------------------------------
    
    if ( fFillOverlapHistograms     || fFillTMHisto        || fFillPerTCardIndexHistograms ||
         fFillEMCALRegionHistograms || fStudyExoticTrigger || fStudyNCellsCut)
    {
      Int_t iclus = -1;
      fCluster = 0;
      fIsExoticTrigger = kFALSE;
      fClusterExoticity = 1;
      
      if      ( GetCalorimeter() == kEMCAL ) { fClustersArr = GetEMCALClusters(); fCaloCells = GetEMCALCells() ; }
      else if ( GetCalorimeter() == kPHOS  ) { fClustersArr = GetPHOSClusters (); fCaloCells = GetPHOSCells () ; }
      
      if ( fClustersArr )
      {
        Int_t  clusterID = aod->GetCaloLabel(0) ;
        
        if ( clusterID < 0 )
          AliWarning(Form("ID of cluster = %d, not possible!", clusterID));
        else
        {
          fCluster = FindCluster(fClustersArr,clusterID,iclus);
          
          // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
          Float_t maxCellFraction = 0.;
          Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(GetEMCALCells(),fCluster,maxCellFraction);
          fTrigSupMod = GetModuleNumber(fCluster);

          if ( fStudyExoticTrigger )
          {
            Int_t bc = GetReader()->GetInputEvent()->GetBunchCrossNumber();
            
            Float_t  eCellMax = fCaloCells->GetCellAmplitude(absIdMax);  
            //Double_t tCellMax = fCaloCells->GetCellTime(absIdMax);      
            
            GetCaloUtils()->RecalibrateCellAmplitude(eCellMax, GetCalorimeter(), absIdMax);
            //GetCaloUtils()->RecalibrateCellTime     (tCellMax, GetCalorimeter(), absIdMax, bc);    
            
            //fClusterExoticity = 1-GetCaloUtils()->GetEMCALRecoUtils()->GetECross(absIdMax,tCellMax,fCaloCells,bc)/eCellMax; // EMCAL
            fClusterExoticity = 1-GetCaloUtils()->GetECross(absIdMax,fCaloCells,bc)/eCellMax; // PHOS and EMCal
            
            if ( fClusterExoticity > 0.97 ) fIsExoticTrigger = kTRUE ;
            //fIsExoticTrigger = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(absIdMax,GetEMCALCells(),bc);
            //if ( fIsExoticTrigger ) 
            //printf("Isolation: IsExotic? %d, E %f, ncells %d, exoticity %f\n", 
            //       fIsExoticTrigger, aod->E(), fCluster->GetNCells(), exoticity);
          }
          
          if ( fStudyNCellsCut )
          {
            // Init for this trigger cluster
            fNCellsWithWeight = 0;
            
            // Loop on cells in cluster to get number of cells with significant energy
            for (Int_t ipos = 0; ipos < fCluster->GetNCells(); ipos++) 
            {
              Int_t   absId = fCluster  ->GetCellsAbsId()[ipos];
              Float_t eCell = fCaloCells->GetCellAmplitude(absId) ;
              
              GetCaloUtils()->RecalibrateCellAmplitude(eCell, GetCalorimeter(), absId);
              
              if( absId == absIdMax || eCell < 0.01 ) continue;
              
              Float_t weight = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell, energy);
              
              if( weight < 0.01 ) continue;
              
              fNCellsWithWeight++;
            } // loop
          } // study n cells
          
          if ( fFillPerTCardIndexHistograms )
          {
            Int_t ietaMax=-1, iphiMax = 0, rcuMax = 0;  
            GetModuleNumberCellIndexes(absIdMax, GetCalorimeter(),ietaMax, iphiMax, rcuMax);
            
            Int_t rowTCard = Int_t(iphiMax%8);
            Int_t colTCard = Int_t(ietaMax%2);
            fTCardIndex    = rowTCard+8*colTCard; 
            //printf("Iso: row TCard %d, col TCard %d, iTCard %d\n",rowTCard,colTCard,fTCardIndex);
          }
          
        } // cluster found
        
      } // cluster array found
    } // recover cluster
    
    // Check number of overlaps in cluster
    //
    Int_t noverlaps = -1;
    if ( fFillOverlapHistograms && fCluster )
    {
      const UInt_t nlabels = fCluster->GetNLabels();
      Int_t overpdg[nlabels];
      Int_t overlab[nlabels];
      noverlaps = GetMCAnalysisUtils()->GetNOverlaps(fCluster->GetLabels(), nlabels, mcTag, -1,
                                                     GetMC(), overpdg, overlab);
    }
    
    //---------------------------------------------------------------
    // Fill pt/sum pT distribution of particles in cone 
    // for different selection cases: per SM, exotics, TCard, ...
    // Similar to the methods in AliIsolationCut, with many more additional histograms
    //---------------------------------------------------------------
  
    if ( fStudyExoticTrigger || fFillPerSMHistograms || fFillPerTCardIndexHistograms  || 
         fStudyPtCutInCone   || fStudyRCutInCone     || fStudyNCellsCut               ||  
         fStudyTracksInCone  || IsPileUpAnalysisOn() || IsHighMultiplicityAnalysisOn()   )    
    {
      StudyTracksInCone   (aod);
      StudyClustersInCone (aod);
      
      if ( fStudyExoticTrigger && fIsExoticTrigger )
        fhConeSumPtExoTrigger  ->Fill(pt, coneptsum, GetEventWeight()*weightTrig);
      
      if ( fFillPerSMHistograms )     
        fhConeSumPtPerSM[aod->GetSModNumber()]->Fill(pt,coneptsum, GetEventWeight()*weightTrig);
      
      if ( fFillPerTCardIndexHistograms )     
        fhConeSumPtPerTCardIndex[fTCardIndex]->Fill(pt,coneptsum, GetEventWeight()*weightTrig);
    }
    
    //---------------------------------------------------------------
    // Check tracks or clusters in UE bands or perpendicular cones to trigger
    //---------------------------------------------------------------
    
    if ( fStudyTracksInCone || fStudyPtCutInCone )
      StudyTracksUEInCone(aod);
    
     if ( fStudyPtCutInCone )
      StudyClustersUEInCone(aod);

    //---------------------------------------------------------------    
    //---------------------------------------------------------------
    // Do not fill in non isolated histogram clusters in isolation gap
    //---------------------------------------------------------------
    //---------------------------------------------------------------

    if ( !isolated && method >= AliIsolationCut::kSumPtIC && coneptsum < sumThr+sumGap ) 
    {
      AliDebug(1,Form("Candidate not isolated but in gap: pt %1.1f, eta %1.1f, phi %1.1f, cone sum pT %2.2f",
                      pt, eta, phi, coneptsum));
      continue;
    }
    
    //---------------------------------------------------------------
    // Fill Shower shape  histograms
    //---------------------------------------------------------------
    
    FillShowerShapeControlHistograms(aod, mcIndex, noverlaps);
    
    //---------------------------------------------------------------
    // Fill track matching histograms
    //---------------------------------------------------------------
    
    FillTrackMatchingControlHistograms(aod, mcIndex);
    
    //---------------------------------------------------------------
    // EMCAL SM regions
    //---------------------------------------------------------------
    
    if ( fFillEMCALRegionHistograms && GetCalorimeter() == kEMCAL ) 
      StudyEMCALRegions(pt, phi, eta, m02, 
                        coneptsumTrack, coneptsumCluster, isolated, iSM);
    
    //---------------------------------------------------------------
    // Conversion radius in MC
    //---------------------------------------------------------------
    
    if ( IsDataMC() && fStudyMCConversionRadius ) 
      StudyMCConversionRadius(pt, iSM, isolated, m02, mcTag, aod->GetLabel());
    
    //---------------------------------------------------------------
    // Fill histograms to undertand pile-up before other cuts applied
    // Remember to relax time cuts in the reader for time related histograms
    //---------------------------------------------------------------
    if ( IsPileUpAnalysisOn() )
      FillPileUpHistograms(aod);//aod->GetCaloLabel(0));
    
    //---------------------------------------------------------------
    // Isolated/ Non isolated, photon/bkg histograms
    //---------------------------------------------------------------
   
    if ( fStudyExoticTrigger && fIsExoticTrigger )
    {
      fhPtExoTrigger[isolated]->Fill(pt, GetEventWeight()*weightTrig);        
    }
    
    // Fill depending shower shape narrow or wide or iso or not iso
    // On non photon analysis, iso or non iso filled
    //
    if ( !inM02Windows || fFillOnlyTH3Histo ) continue; // it is on the wide or narrow window if those are selected, see above
    
    fhPt[isolated][narrow]->Fill(pt    , GetEventWeight()*weightTrig);
    
    fhPtEtaPhi[isolated][narrow]->Fill(pt, eta, phi, GetEventWeight()*weightTrig);
    
    if ( IsDataMC() && mcIndex < fNumberMCParticleCases && !GetReader()->AreMCPromptPhotonsSelected() )
    {
      // For histograms in arrays, index in the array, corresponding to any particle origin
      if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
      {
        fhPtMC      [kmcPhoton][isolated][narrow]->Fill(  pt, GetEventWeight()*weightTrig);
        fhPtEtaPhiMC[kmcPhoton][isolated]->Fill(pt, eta, phi, GetEventWeight()*weightTrig);
      }
      
      if ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) && 
          fNumberMCParticleCases > kmcPi0DecayLostPair )
      {
        if      ( mcIndex == kmcPi0Decay )
        {
          fhPtMC      [kmcPi0DecayLostPair][isolated][narrow]->Fill(  pt, GetEventWeight()*weightTrig);
          fhPtEtaPhiMC[kmcPi0DecayLostPair][isolated]->Fill(pt, eta, phi, GetEventWeight()*weightTrig);
        }
        else if ( mcIndex == kmcEtaDecay )
        {
          fhPtMC      [kmcEtaDecayLostPair][isolated][narrow]->Fill(  pt, GetEventWeight()*weightTrig);
          fhPtEtaPhiMC[kmcEtaDecayLostPair][isolated]->Fill(pt, eta, phi, GetEventWeight()*weightTrig);
        }
      }
      
      fhPtMC      [mcIndex][isolated][narrow]->Fill(  pt, GetEventWeight()*weightTrig);
      fhPtEtaPhiMC[mcIndex][isolated]->Fill(pt, eta, phi, GetEventWeight()*weightTrig);
    }//Histograms with MC
    
    if ( fFillNLMHistograms )
      fhPtNLocMax[isolated][narrow] ->Fill(pt, aod->GetNLM(), GetEventWeight()*weightTrig) ;
    
    if ( IsHighMultiplicityAnalysisOn() )
    {
      fhPtCentrality[isolated][narrow]->Fill(pt, GetEventCentrality(), GetEventWeight()*weightTrig) ;
      fhPtEventPlane[isolated][narrow]->Fill(pt, GetEventPlaneAngle(), GetEventWeight()*weightTrig) ;
    }
  }// aod loop
}

//______________________________________________________
/// Fill primary generated particle acceptance histograms
/// if MC data is available. Only when particle is in the calorimeter.
/// Rethink if CTS is used.
//______________________________________________________
void AliAnaParticleIsolation::FillAcceptanceHistograms()
{
  if ( !GetMC() ) return;

  AliDebug(1,"Start");
  
  Double_t photonY   = -100 ;
  Double_t photonE   = -1 ;
  Double_t photonPt  = -1 ;
  Double_t photonPhi =  100 ;
  Double_t photonEta = -1 ;
  
  Int_t    pdg       =  0 ;
  Int_t    status    =  0 ;
  Int_t    tag       =  0 ;
  Int_t    mcIndex   =  0 ;
  Int_t    nprim     = GetMC()->GetNumberOfTracks();
  
  Bool_t   ok        = kFALSE;
  Int_t    momLabel  = -1;
  
  AliVParticle * primary = 0;
  
  TString genName = "";
  
  // Calorimeter cluster merging angle
  // angle smaller than 3 cells  6 cm (0.014) in EMCal, 2.2 cm in PHOS (0.014*(2.2/6))
  Float_t overlapAngle = 0;
  Float_t minECalo     = 0;
  if      (GetCalorimeter()==kEMCAL)
  {
    overlapAngle = fMinCellsAngleOverlap*0.014  ;
    minECalo = GetReader()->GetEMCALEMin();
  }
  else if (GetCalorimeter()==kPHOS )
  {
    overlapAngle = fMinCellsAngleOverlap*0.00382;
    minECalo = GetReader()->GetPHOSEMin();
  }

  Int_t firstParticle = 0;
  
  // Loop only over likely final particles not partons
  if ( GetReader()->GetGenPythiaEventHeader() ) 
    firstParticle = GetMCAnalysisUtils()->GetPythiaMaxPartParent();
  
  Int_t partInConeType = GetIsolationCut()->GetParticleTypeInCone() ;
  Int_t isoMethod      = GetIsolationCut()->GetICMethod();
  Float_t coneSize     = GetIsolationCut()->GetConeSize();
  Float_t coneSizeGap  = GetIsolationCut()->GetConeSizeBandGap();

  Float_t centrality   = GetEventCentrality();
  
  for(Int_t i = firstParticle ; i < nprim; i++)
  {
    if ( !GetReader()->AcceptParticleMCLabel( i ) ) continue ;
    
    primary = GetMC()->GetTrack(i) ;
    if(!primary)
    {
      AliWarning("primaries pointer not available!!");
      continue;
    }
    
    pdg    = primary->PdgCode();
    status = primary->MCStatusCode();
    
    // Protection against floating point exception
    if ( primary->E() == TMath::Abs(primary->Pz()) || 
        (primary->E() - primary->Pz()) < 1e-3      ||
        (primary->E() + primary->Pz()) < 0           )  continue ; 
    
    //printf("i %d, %s %d  %s %d \n",i, GetMC()->Particle(i)->GetName(), GetMC()->Particle(i)->GetPdgCode(),
    //       primary->GetName(), primary->GetPdgCode());
    
    // Select only photons in the final state
    if(pdg != 22  && pdg!=111 && pdg !=221) continue ;
    
    // Consider only final state particles, but this depends on generator,
    // status 1 is the usual one, in case of not being ok, leave the possibility
    // to not consider this.
    if( pdg == 22 && status != 1 &&
       GetMCAnalysisUtils()->GetMCGenerator() != AliMCAnalysisUtils::kBoxLike ) continue ;
    
    // Not interested in large rapidity, cut those
    photonY = 0.5*TMath::Log((primary->E()+primary->Pz())/(primary->E()-primary->Pz())) ;
    if ( TMath::Abs(photonY) > 1.0 ) continue;
    
    // Photon kinematics
    primary->Momentum(fMomentum);
    
    // If too small or too large pt, skip, same cut as for data analysis
    photonPt  = fMomentum.Pt () ;
    
    if ( photonPt < GetMinPt() || photonPt > GetMaxPt() ) continue ;
    
    // Avoid too large pT particle in pT hard binned productions
    if ( GetReader()->GetGenPythiaEventHeader() ) 
    {
      Float_t pTHard = (GetReader()->GetGenPythiaEventHeader())->GetPtHard();
      
      if ( pTHard / photonPt < 0.5 ) 
      {
        AliInfo(Form("Reject primary particle pdg %d, pT %f larger than pT hard %f, ratio %f",
                     pdg, photonPt, pTHard, pTHard/photonPt));
        //GetMCAnalysisUtils()->PrintAncestry(GetMC(),i);
        continue;
      }
    }
    
    photonE   = fMomentum.E  () ;
    photonEta = fMomentum.Eta() ;
    photonPhi = fMomentum.Phi() ;
    if ( photonPhi < 0 ) photonPhi+=TMath::TwoPi();
    
    // Check if photons hit the Calorimeter acceptance
    if ( IsRealCaloAcceptanceOn() && fIsoDetector!=kCTS ) // defined on base class
    {
      if ( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primary) ) continue ;
    }
    
    // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
    if ( !GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),fIsoDetector) ) continue ;
    
    // Get tag of this particle photon from fragmentation, decay, prompt ...
    // Set the origin of the photon.
    tag = GetMCAnalysisUtils()->CheckOrigin(i, GetMC(),
                                            GetReader()->GetNameOfMCEventHederGeneratorToAccept(),
                                            photonE); // Not used, should be cluster
    
    if ( GetReader()->AreMCPromptPhotonsSelected() && 
        !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) ) continue;

    if ( GetReader()->AreMCFragmentationPhotonsRejected() ) 
    {
      if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) ) continue;
      if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) ) continue;
      if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) ) continue;
    }
    
    if ( pdg == 22 && !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
    {
      // A conversion photon from a hadron, skip this kind of photon
      // printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!\n ");
      // GetMCAnalysisUtils()->PrintMCTag(tag);
      
      continue;
    }
    
    // Particle ID and pT dependent Weight
    Int_t   index    = GetReader()->GetCocktailGeneratorAndIndex(i, genName);
    Float_t weightPt = GetParticlePtWeight(photonPt, pdg, genName, index) ; 
    //
    
    // Check the origin of the photon or if it is a pi0, assing a tag
    Int_t pi0d1Label = -1, pi0d2Label = -1;
    Bool_t overlapPi0 = kTRUE;
    if( pdg == 111 || pdg == 221 )
    {
      if( pdg == 111 ) mcIndex = kmcPrimPi0;
      else             mcIndex = kmcPrimEta;
      
      //printf("check meson pdg %d\n",pdg);
      
      // Get the labels of the decay particles, remove them from isolation cone
      // Get also the opening angle and check if decays likely overlap
      Bool_t okpi0 = kFALSE;
      Int_t ndaugh = GetMCAnalysisUtils()->GetNDaughters(i,GetMC(), okpi0);
      // printf("OK pi0 %d, ndaugh %d\n",okpi0,ndaugh);
      Int_t d1Pdg = 0, d1Status = 0; Bool_t ok1 = kFALSE;
      Int_t d2Pdg = 0, d2Status = 0; Bool_t ok2 = kFALSE;
      
      if ( ndaugh > 0 ) fMomDaugh1 = GetMCAnalysisUtils()->GetDaughter(0,i,GetMC(),d1Pdg, d1Status,ok1, pi0d1Label,fProdVertex);
      if ( ndaugh > 1 ) fMomDaugh2 = GetMCAnalysisUtils()->GetDaughter(1,i,GetMC(),d2Pdg, d2Status,ok2, pi0d2Label,fProdVertex);
      
      //printf("pi0 daug %d: a) %d, b) %d\n", ndaugh,pi0d1Label,pi0d2Label);
      //if ( ndaugh !=2 ) printf("PDG: %d, %d\n",d1Pdg,d2Pdg);
      
      // Select decays in 2 photons
      if ( ndaugh!=2 || (d2Pdg != d1Pdg && d1Pdg!=22) ) okpi0 = kFALSE;
      
      // Check overlap of decays
      if ( okpi0 && fMakePrimaryPi0DecayStudy )
      {
        Float_t d12Angle = fMomDaugh1.Angle(fMomDaugh2.Vect());
        if(d12Angle > overlapAngle) overlapPi0 = kFALSE;
        //printf("  -- d12 angle %2.3f, angle limit %2.3f, overlap %d\n",d12Angle,overlapAngle,overlapPi0);
      }
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) )
    {
      mcIndex = kmcPrimPrompt;
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) )
    {
      mcIndex = kmcPrimFrag ;
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) )
    {
      mcIndex = kmcPrimISR;
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) )
    {
      mcIndex   = kmcPrimPi0Decay;
      fMomentum = GetMCAnalysisUtils()->GetMotherWithPDG(i, 111, GetMC(),ok, momLabel);        
      weightPt  = GetParticlePtWeight(fMomentum.Pt(), 111, genName, index) ; 
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) )
    {
      mcIndex   = kmcPrimEtaDecay;
      fMomentum = GetMCAnalysisUtils()->GetMotherWithPDG(i, 221, GetMC(),ok, momLabel);        
      weightPt  = GetParticlePtWeight(fMomentum.Pt(), 221, genName, index) ; 
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) )
    {
      mcIndex = kmcPrimOtherDecay;
    }
    else
    {
      // Other decay but from non final state particle
      mcIndex = kmcPrimOtherDecay;
    }//Other origin
    
    //printf("mcIndex %d\n",mcIndex);
    
    // ////////////////////ISO MC/////////////////////////
    Float_t dR              = 0. ;
    Float_t sumPtInCone     = 0. ;
    Float_t sumPtInConeCh   = 0. ;
    Float_t sumPtInConeNe   = 0. ;

    Float_t etaBandPtSumCh  = 0. ;
    Float_t phiBandPtSumCh  = 0. ;
    Float_t perpConePtSumCh = 0. ;
    Float_t perpConePtSumNe = 0. ;

    Float_t etaBandPtSumNe  = 0. ;
    Float_t phiBandPtSumNe  = 0. ;
    
    AliVParticle * mcisop = 0;
    
    Int_t partInConeStatus = -1, partInConeMother = -1;
    Double_t partInConePt = 0, partInConeE = 0, partInConeEta = 0, partInConePhi = 0;
    Int_t partInConeCharge = 0, npart = 0, partInConePDG = 0;
    Bool_t physPrimary = kFALSE;
    
    for(Int_t ip = firstParticle; ip < nprim ; ip++)
    {
      if ( ip==i ) continue;
      
      if( (pdg==111 || pdg==221) && ( ip == pi0d1Label || ip == pi0d2Label ) )
      {
        //printf("Do not count pi0 decays in cone when isolating pi0 \n");
        continue;
      }
      
      mcisop = GetMC()->GetTrack(ip);
      if( !mcisop ) continue;
      
      partInConeStatus = mcisop->MCStatusCode();
      
      // Consider only final state particles, but this depends on generator,
      // status 1 is the usual one, in case of not being ok, leave the possibility
      // to not consider this.
      if( partInConeStatus  != 1 &&
         GetMCAnalysisUtils()->GetMCGenerator()!= AliMCAnalysisUtils::kBoxLike ) continue ;
      
      // Protection against floating point exception
      if ( mcisop->E() == TMath::Abs(mcisop->Pz())    ||
          (mcisop->E() - mcisop->Pz()) < 1e-3         ||
           mcisop->Pt() < 0.01  || mcisop->Eta() > 10 ||
          (mcisop->E() + mcisop->Pz()) < 0               )  continue ; 
      
      partInConeMother = mcisop->GetMother();
      partInConePt     = mcisop->Pt();
      partInConeE      = mcisop->E();
      partInConeEta    = mcisop->Eta();
      partInConePhi    = GetPhi(mcisop->Phi());
      partInConePDG    = mcisop->PdgCode();
      //partInConeCharge = TMath::Abs((Int_t) TDatabasePDG::Instance()->GetParticle(partInConePDG)->Charge());
      partInConeCharge = TMath::Abs(mcisop->Charge());
      physPrimary      = mcisop->IsPhysicalPrimary(); 
      mcisop->Momentum(fMomIso);
      
      if ( partInConeMother == i ) continue;

      // Distance to photon candidate
      //
      dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);

      // Angles between trigger and track
      Double_t dEta = photonEta - partInConeEta;
      Double_t dPhi = photonPhi - partInConePhi;

      // Shift phi angle when trigger is close to 0 or 360
      if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
      if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();

      //
      // Apply acceptance and energy/pt cut for particles in cone
      // select only charged final particles and photons
      if ( fSelectPrimariesInCone )
      {
        if ( partInConeCharge > 0 ) // charged pT cut and acceptance
        {
          if ( partInConeType == AliIsolationCut::kOnlyNeutral ) continue;

          if ( partInConePt < GetReader()->GetCTSPtMin () ) continue;
          
          if ( !GetReader()->GetFiducialCut()->IsInFiducialCut(fMomIso.Eta(),fMomIso.Phi(),kCTS) ) continue ;
          
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // Fill the histograms at +-45 degrees in phi
          // from trigger particle, perpedicular to trigger axis in phi
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if ( dR > coneSize && physPrimary &&
               isoMethod == AliIsolationCut::kSumBkgSubIC  )
          {
            Double_t dPhiPlu = dPhi + TMath::PiOver2();
            Double_t dPhiMin = dPhi - TMath::PiOver2();

            Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
            Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;

            Bool_t fillPerp = kFALSE;
            if ( TMath::Sqrt(argPlu) < coneSize ) fillPerp = kTRUE ;
            if ( TMath::Sqrt(argMin) < coneSize ) fillPerp = kTRUE ;

            if ( fillPerp )
              perpConePtSumCh += partInConePt;
          }
        }
        else  // photons E cut and acceptance, discard hadrons (neutron, k0L, ...)
        {
          if ( partInConeType == AliIsolationCut::kOnlyCharged ) continue;

          if ( partInConePDG != 22 ) continue; // just photons
          
          if ( partInConeE  <= minECalo ) continue;
          
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // Fill the histograms at +-45 degrees in phi
          // from trigger particle, perpedicular to trigger axis in phi
          // Apply only TPC acceptance cut for comparison with perpendicular charged
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if ( dR > coneSize && physPrimary &&
               isoMethod == AliIsolationCut::kSumBkgSubIC &&
               GetReader()->GetFiducialCut()->IsInFiducialCut(fMomIso.Eta(),fMomIso.Phi(),kCTS) )
          {
            Double_t dPhiPlu = dPhi + TMath::PiOver2();
            Double_t dPhiMin = dPhi - TMath::PiOver2();

            Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
            Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;

            Bool_t fillPerp = kFALSE;
            if ( TMath::Sqrt(argPlu) < coneSize ) fillPerp = kTRUE ;
            if ( TMath::Sqrt(argMin) < coneSize ) fillPerp = kTRUE ;

            if ( fillPerp )
              perpConePtSumNe += partInConePt;
          }
          
          if ( !GetReader()->GetFiducialCut()->IsInFiducialCut(fMomIso.Eta(),fMomIso.Phi(),GetCalorimeter()) ) continue ;
          
          if ( IsRealCaloAcceptanceOn() ) // defined on base class
          {
            if ( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), mcisop) ) continue ;
          }
        }
      }
      
//      printf("Input status %d, E %f, pt %f, pz %f, pdg %d, gamma eta %f, phi %f, par in cone eta %f phi %f\n",
//             partInConeStatus, partInConeE, partInConePt, mcisop->Pz(), partInConePDG, 
//             photonEta, photonPhi, partInConeEta, partInConePhi);
      
      // UE estimation
      //
      if ( physPrimary && dR > coneSize+coneSizeGap )
      {
        if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
        {
          // Phi band
          //
          Bool_t takeIt = kTRUE;

          // Look only half TPC with respect candidate, avoid opposite side jet
          if ( TMath::Abs(dPhi) > TMath::PiOver2() )  takeIt = kFALSE;

          // Within eta cone size
          if ( TMath::Abs(dEta) < coneSize+coneSizeGap && takeIt )
          {
            if ( partInConeCharge > 0 )
              phiBandPtSumCh += partInConePt;
            else
              phiBandPtSumNe += partInConePt;
          } // phi band

          // Eta band
          //
          //takeIt = kTRUE;

          // Within phi cone size
          if ( TMath::Abs(dPhi) < coneSize+coneSizeGap )//&& takeIt )
          {
            if ( partInConeCharge > 0 )
              etaBandPtSumCh += partInConePt;
            else
              etaBandPtSumNe += partInConePt;
          } // eta band

        } // out of cone

      } // Physical primary, out of cone

      // In cone
      //
      if ( dR < GetIsolationCut()->GetMinDistToTrigger() )
        continue;
      
      if ( dR > GetIsolationCut()->GetConeSize() )
        continue;
      
      if ( fFillTrackOriginHistograms && 
           partInConeCharge > 0 &&  TMath::Abs(partInConePDG) != 11 ) // exclude electrons and neutrals
      {
        Int_t mcChTag = 3;
        if      ( TMath::Abs(partInConePDG) == 211  )  mcChTag = 0;
        else if ( TMath::Abs(partInConePDG) == 321  )  mcChTag = 1; 
        else if ( TMath::Abs(partInConePDG) == 2212 )  mcChTag = 2; 
        
        if ( physPrimary )
          fhPtTrackInConeMCPrimaryGener  [mcChTag]->Fill(photonPt, partInConePt, GetEventWeight()*weightPt);
        else
          fhPtTrackInConeMCSecondaryGener[mcChTag]->Fill(photonPt, partInConePt, GetEventWeight()*weightPt);
        
        //printf("Selected particles pdg %d, status %d\n", partInConePDG, partInConeStatus);
      }
      
      if ( !physPrimary ) continue ; 
      
      if ( partInConeCharge > 0 )
      {
        sumPtInConeCh += partInConePt;
        //printf("Charged %d pT %f\n",ip,partInConePt);
      }
      else
        sumPtInConeNe += partInConePt;
      
      if ( partInConePt > GetIsolationCut()->GetPtThreshold() &&
           partInConePt < GetIsolationCut()->GetPtThresholdMax() ) npart++;
      
//      printf("Selected: status %d, E %f, pt %f, pz %f, pdg %d, eta %f, phi %f, par in cone eta %f phi %f\n",
//             partInConeStatus, partInConeE, partInConePt, mcisop->Pz(), partInConePDG, 
//             photonEta, photonPhi, partInConeEta, partInConePhi);      
    }
    
//    printf("FillAcc::Generated Method %d cone %0.2f: sumPt Ch %2.1f, Ne %2.1f; "
//           "eta band Ch %2.1f, Ne %2.1f; phi band Ch %2.1f, Ne %2.1f; perp Ch %2.1f \n ",
//           isoMethod,coneSize,sumPtInConeCh,sumPtInConeNe,
//           etaBandPtSumCh,etaBandPtSumNe,phiBandPtSumCh,phiBandPtSumNe,perpConePtSumCh);
    ///////END ISO MC/////////////////////////
    
    // ////////////////////ISO EMBED DATA/////////////////////////
    Float_t sumPtInConeEmb     = 0. ;
    Float_t sumPtInConeChEmb   = 0. ;
    Float_t sumPtInConeNeEmb   = 0. ;

    Float_t etaBandPtSumChEmb  = 0. ;
    Float_t phiBandPtSumChEmb  = 0. ;
    Float_t perpConePtSumChEmb = 0. ;

    Float_t etaBandPtSumNeEmb  = 0. ;
    Float_t phiBandPtSumNeEmb  = 0. ;
    Float_t gapCandidate = 0.0143*2.5;
    if ( GetIsolationCut()->GetMinDistToTrigger() > 0 )
      gapCandidate = GetIsolationCut()->GetMinDistToTrigger();

    if ( IsEmbedingAnalysisOn() )
    {
      // Clusters
      // Recover  clusters
      TObjArray * plNe = 0x0;
      if      (GetCalorimeter()==kPHOS )
        plNe = GetReader()->GetPHOSClusters();
      else if (GetCalorimeter()==kEMCAL)
        plNe = GetReader()->GetEMCALClusters();

      // Get vertex for cluster momentum calculation
      Double_t vertex[] = {0,0,0} ; //vertex ;
      if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
        GetReader()->GetVertex(vertex);

      for(Int_t icalo=0; icalo < plNe->GetEntriesFast(); icalo++)
      {
        AliVCluster * calo = dynamic_cast<AliVCluster *>(plNe->At(icalo)) ;

        // const Int_t mcLabel = calo->GetLabel();
        const Int_t nlabels = calo->GetNLabels();

        if ( nlabels != 0  ) continue; // Accept only purely MC clusters

        // Reject track-matched clusters
        if ( GetIsolationCut()->IsTrackMatchedClusterRejectionInConeOn()  )
        {
          Bool_t bRes = kFALSE, bEoP = kFALSE;
          Bool_t matched = GetCaloPID()->IsTrackMatched(calo, GetReader()->GetCaloUtils(),
                                                        GetReader()->GetInputEvent(),
                                                        bEoP,bRes);
          if ( partInConeType == AliIsolationCut::kNeutralAndCharged && matched ) continue ;
        }

        calo->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line

        partInConePt     = fMomentum.Pt();
        partInConeE      = fMomentum.E();
        partInConeEta    = fMomentum.Eta();
        partInConePhi    = GetPhi(fMomentum.Phi());

        // Distance to photon candidate
        //
        dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);

        // Angles between trigger and cluster
        Double_t dEta = photonEta - partInConeEta;
        Double_t dPhi = photonPhi - partInConePhi;

        // Shift phi angle when trigger is close to 0 or 360
        if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
        if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();

        // UE estimation
        //
        if (  dR > coneSize+coneSizeGap )
        {
          if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
          {
            // Phi band
            //
            Bool_t takeIt = kTRUE;

            // Look only half TPC with respect candidate, avoid opposite side jet
            if ( TMath::Abs(dPhi) > TMath::PiOver2() )  takeIt = kFALSE;

            // Within eta cone size
            if ( TMath::Abs(dEta) < coneSize+coneSizeGap && takeIt )
            {
              phiBandPtSumNeEmb += partInConePt;
            } // phi band

            // Eta band
            //
            //takeIt = kTRUE;

            // Within phi cone size
            if ( TMath::Abs(dPhi) < coneSize+coneSizeGap )//  && takeIt )
            {
              etaBandPtSumNeEmb += partInConePt;
            } // eta band

          } // out of cone

        } // Physical primary, out of cone

        // In cone
        //
        // Avoid clusters too close to generated photon
        if ( dR < gapCandidate )
          continue;

        if ( dR > GetIsolationCut()->GetConeSize() )
          continue;

        sumPtInConeNeEmb += partInConePt;
      }

      // Tracks
      TObjArray * plCh = GetReader()->GetCTSTracks();
      for(Int_t itrack=0; itrack < plCh->GetEntriesFast(); itrack++)
      {
        AliVTrack* track = dynamic_cast<AliVTrack*>(plCh->At(itrack));

        Int_t label  = track->GetLabel();

        if ( label >= 0  ) continue; // Accept only purely MC clusters

        partInConePt     = track->Pt();
        partInConeEta    = track->Eta();
        partInConePhi    = GetPhi(track->Phi());

        // Distance to photon candidate
        //
        dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);

        // Angles between trigger and track
        Double_t dEta = photonEta - partInConeEta;
        Double_t dPhi = photonPhi - partInConePhi;

        // Shift phi angle when trigger is close to 0 or 360
        if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
        if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();

        // UE estimation
        //
        if (  dR > coneSize+coneSizeGap )
        {
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // Fill the histograms at +-45 degrees in phi
          // from trigger particle, perpedicular to trigger axis in phi
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if ( isoMethod == AliIsolationCut::kSumBkgSubIC  )
          {
            Double_t dPhiPlu = dPhi + TMath::PiOver2();
            Double_t dPhiMin = dPhi - TMath::PiOver2();

            Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
            Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;

            Bool_t fillPerp = kFALSE;
            if ( TMath::Sqrt(argPlu) < coneSize ) fillPerp = kTRUE ;
            if ( TMath::Sqrt(argMin) < coneSize ) fillPerp = kTRUE ;

            if ( fillPerp )
              perpConePtSumChEmb += partInConePt;
          }

          // BANDS

          if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
          {
            // Phi band
            //
            Bool_t takeIt = kTRUE;

            // Look only half TPC with respect candidate, avoid opposite side jet
            if ( TMath::Abs(dPhi) > TMath::PiOver2() )  takeIt = kFALSE;

            // Within eta cone size
            if ( TMath::Abs(dEta) < coneSize+coneSizeGap  &&  takeIt )
            {
              phiBandPtSumChEmb += partInConePt;
            } // phi band

            // Eta band
            //
            //takeIt = kTRUE;

            // Within phi cone size
            if ( TMath::Abs(dPhi) < coneSize+coneSizeGap ) // && takeIt )
            {
              etaBandPtSumChEmb += partInConePt;
            } // eta band

          } // out of cone

        } // Physical primary, out of cone

        // In cone
        //
        if ( dR < gapCandidate )
          continue;

        if ( dR > GetIsolationCut()->GetConeSize() )
          continue;

        sumPtInConeChEmb += partInConePt;
        //printf("Track %d pT %f\n",itrack, partInConePt);

      }

//      printf("FillAcc::Embedded Method %d cone %0.2f: "
//             "sumPt Ch %2.1f, Ne %2.1f; eta band Ch %2.1f, Ne %2.1f; "
//             "phi band Ch %2.1f, Ne %2.1f; perp Ch %2.1f \n ",
//             isoMethod,coneSize,sumPtInConeChEmb,sumPtInConeNeEmb,
//             etaBandPtSumChEmb,etaBandPtSumNeEmb,phiBandPtSumChEmb,
//             phiBandPtSumNeEmb,perpConePtSumChEmb);

      // Add generated particles to embedded
      sumPtInConeChEmb+=sumPtInConeCh;
      sumPtInConeNeEmb+=sumPtInConeNe;

      etaBandPtSumChEmb+=etaBandPtSumCh;
      phiBandPtSumChEmb+=phiBandPtSumCh;
      etaBandPtSumNeEmb+=etaBandPtSumNe;
      phiBandPtSumNeEmb+=phiBandPtSumNe;

      perpConePtSumChEmb+=perpConePtSumCh;
    }
    ///////END ISO EMBED DATA/////////////////////////

    // Get detectors acceptance
    // Do it once, needed for Band UE estimation and excess area determination
    //if ( fMakeConeExcessCorr || fICMethod >= kSumBkgSubEtaBandIC )  // comment not needed in all cases, but it does not take much, in case UE study done in AliAnaParticleIsolation
    GetIsolationCut()->GetDetectorAngleLimits(GetReader(),GetCalorimeter());

    // Calculate how much of the cone got out the detectors acceptance
    //
    Float_t excessAreaChEta = 1;
    Float_t excessAreaNeEta = 1;
    Float_t excessAreaNePhi = 1;
    Float_t excessChEta     = 0;
    Float_t excessNeEta     = 0;
    Float_t excessNePhi     = 0;
    if ( GetIsolationCut()->IsConeExcessCorrected() )
    {
      GetIsolationCut()->CalculateExcessAreaFractionForChargedAndNeutral(photonEta  , photonPhi,
                                                                         excessChEta, excessAreaChEta,
                                                                         excessNeEta, excessAreaNeEta,
                                                                         excessNePhi, excessAreaNePhi);
      sumPtInConeCh *= excessAreaChEta ;
      sumPtInConeNe *= (excessAreaNeEta*excessAreaNePhi) ;

      sumPtInConeChEmb *= excessAreaChEta ;
      sumPtInConeNeEmb *= (excessAreaNeEta*excessAreaNePhi) ;
    }

    //printf("Centrality %d\n",GetEventCentrality());
    //printf("SumPt in cone ne %2.2f, ch %2.2f;\n \t emb ne %2.2f, ch %2.2f\n",
    //       sumPtInConeNe, sumPtInConeCh, sumPtInConeNeEmb, sumPtInConeChEmb);

    // Total energy in cone
    //
    sumPtInCone = sumPtInConeNe + sumPtInConeCh  ;
    sumPtInConeEmb = sumPtInConeNeEmb + sumPtInConeChEmb  ;

    // Normalize estimated UE to cone size
    //
    Double_t coneptsumUESub   = sumPtInCone;
    Double_t coneptsumUESubNe = sumPtInConeNe;
    Double_t coneptsumUESubCh = sumPtInConeCh;

    Double_t coneptsumBkgCh    = 0.;
    Double_t coneptsumBkgNe    = 0.;

    //Double_t coneptsumBkgChRaw = 0.;
    //Double_t coneptsumBkgNeRaw = 0.;

    Double_t coneptsumUESubEmb   = sumPtInConeEmb;
    Double_t coneptsumUESubNeEmb = sumPtInConeNeEmb;
    Double_t coneptsumUESubChEmb = sumPtInConeChEmb;

    Double_t coneptsumBkgChEmb    = 0.;
    Double_t coneptsumBkgNeEmb    = 0.;

    //printf("SumPt in UE perp %2.2f; phi ne %2.2f; ch %2.2f; eta ne %2.2f, ch %2.2f\n"
    //       "\t Emb: perp %2.2f; phi ne %2.2f; ch %2.2f; eta ne %2.2f, ch %2.2f\n", 
    //       perpConePtSumCh, etaBandPtSumNe, etaBandPtSumCh, phiBandPtSumNe, phiBandPtSumCh,
    //       perpConePtSumChEmb, etaBandPtSumNeEmb, etaBandPtSumChEmb, phiBandPtSumNeEmb, phiBandPtSumChEmb
    //       );

    if ( isoMethod == AliIsolationCut::kSumBkgSubIC )
    {
      coneptsumBkgCh = perpConePtSumCh/2.;
      coneptsumBkgNe = perpConePtSumNe/2.;

      coneptsumBkgChEmb = perpConePtSumChEmb/2.;
      //printf("centrality %f, neutral/charged %f\n",centrality,GetNeutralOverChargedRatio(centrality));
    }
    else if ( isoMethod >= AliIsolationCut::kSumBkgSubEtaBandIC ) // eta or phi band
    {
      //printf("UE band\n");
      Float_t  etaBandPtSumChNorm = 0;
      Float_t  phiBandPtSumChNorm = 0;
      Float_t  etaBandPtSumNeNorm = 0;
      Float_t  phiBandPtSumNeNorm = 0;
      Float_t  etaBandPtSumChNormEmb = 0;
      Float_t  phiBandPtSumChNormEmb = 0;
      Float_t  etaBandPtSumNeNormEmb = 0;
      Float_t  phiBandPtSumNeNormEmb = 0;
      //Float_t  coneptsumChSub     = 0 ;
      //Float_t  coneptsumNeSub     = 0 ;

      // Normalize background to cone area

      // Check clusters band
      //
      // If DCal or PHOS trigger clusters, it does not make sense,
      // just use track based estimation later
      Bool_t checkClustersBand = kTRUE;

      if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
      {
        if      ( GetCalorimeterString() != "EMCAL" )
          checkClustersBand = kFALSE;
        // if declared as emcal because analyze all calo together
        // but within phos/dcal phi acceptance
        else if ( photonPhi > 4.35 && photonPhi < 5.71 )
          checkClustersBand = kFALSE;

        //        if ( !checkClustersBand ) printf("Cluster band not checked! det %d phi %f\n",
        //                                         pCandidate->GetDetectorTag(),photonPhi);
      }

      if ( partInConeType !=  AliIsolationCut::kOnlyCharged && checkClustersBand )
      {
        GetIsolationCut()->CalculateUEBandClusterNormalization
        (photonEta         , photonPhi         ,
         excessNeEta       , excessNePhi       ,
         excessAreaNeEta   , excessAreaNePhi   ,
         etaBandPtSumNe    , phiBandPtSumCh    ,
         etaBandPtSumNeNorm, phiBandPtSumChNorm);

        if ( IsEmbedingAnalysisOn () )
        {
          GetIsolationCut()->CalculateUEBandClusterNormalization
          (photonEta         , photonPhi         ,
           excessNeEta       , excessNePhi       ,
           excessAreaNeEta   , excessAreaNePhi   ,
           etaBandPtSumNeEmb    , phiBandPtSumChEmb    ,
           etaBandPtSumNeNormEmb, phiBandPtSumChNormEmb);
        }

        if      ( isoMethod == AliIsolationCut::kSumBkgSubEtaBandIC )
        {
          //coneptsumBkgNeRaw  = etaBandPtSumNe;
          coneptsumBkgNe     = etaBandPtSumNeNorm;
          coneptsumBkgNeEmb  = etaBandPtSumNeNormEmb;
        }
        else  if( isoMethod == AliIsolationCut::kSumBkgSubPhiBandIC )
        {
          //coneptsumBkgNeRaw  = phiBandPtSumNe;
          coneptsumBkgNe     = phiBandPtSumNeNorm;
          coneptsumBkgNeEmb  = phiBandPtSumNeNormEmb;
        }

        //coneptsumNeSub = sumPtInConeNe - coneptsumBkgNe;

        //        printf("Neutral: sumpT %2.2f, \n \t phi: sum pT %2.2f, sumpT norm %2.2f;\n"
        //               "\t eta: sum pT %2.2f, sumpT norm %2.2f;\n"
        //               "\t subtracted:  %2.2f\n",
        //               coneptsumNe,
        //               phiBandPtSumNe, phiBandPtSumNeNorm,
        //               etaBandPtSumNe, etaBandPtSumNeNorm,
        //               coneptsumNeSub);
      } // Neutrals in cone

      //printf("Pass cluster\n");

      if ( partInConeType != AliIsolationCut::kOnlyNeutral )
      {
        GetIsolationCut()->CalculateUEBandTrackNormalization
        (photonEta         ,
         excessChEta       , excessAreaChEta,
         etaBandPtSumCh    , phiBandPtSumCh,
         etaBandPtSumChNorm, phiBandPtSumChNorm);

        if ( IsEmbedingAnalysisOn() )
        {
          GetIsolationCut()->CalculateUEBandTrackNormalization
          (photonEta         ,
           excessChEta       , excessAreaChEta,
           etaBandPtSumChEmb    , phiBandPtSumChEmb,
           etaBandPtSumChNormEmb, phiBandPtSumChNormEmb);
        }

        if      ( isoMethod == AliIsolationCut::kSumBkgSubEtaBandIC )
        {
          //coneptsumBkgChRaw = etaBandPtSumCh;
          coneptsumBkgCh    = etaBandPtSumChNorm;
          coneptsumBkgChEmb = etaBandPtSumChNormEmb;
        }
        else  if( isoMethod == AliIsolationCut::kSumBkgSubPhiBandIC )
        {
          //coneptsumBkgChRaw = phiBandPtSumCh;
          coneptsumBkgCh    = phiBandPtSumChNorm;
          coneptsumBkgChEmb = phiBandPtSumChNormEmb;
        }

        //coneptsumChSub  = sumPtInConeCh - coneptsumBkgCh;

        //        printf("Charged: sumpT %2.2f, \n \t phi: sum pT %2.2f, sumpT norm %2.2f;\n"
        //               "\t eta: sum pT %2.2f, sumpT norm %2.2f;\n"
        //               "\t subtracted: %2.2f\n",
        //               coneptsumCh,
        //               phiBandPtSumCh, phiBandPtSumChNorm,
        //               etaBandPtSumCh, etaBandPtSumChNorm,
        //               coneptsumChSub);

      } // Charged in cone

      // In case of DCal and PHOS, estimate the bkg in clusters from tracks
      //
      if ( partInConeType == AliIsolationCut::kNeutralAndCharged && !checkClustersBand )
      {
        coneptsumBkgNe = coneptsumBkgCh*GetIsolationCut()->GetNeutralOverChargedRatio(centrality);

        if ( IsEmbedingAnalysisOn() )
          coneptsumBkgNeEmb = coneptsumBkgChEmb*GetIsolationCut()->GetNeutralOverChargedRatio(centrality);

        //coneptsumNeSub = sumPtInConeNe - coneptsumBkgNe;
      }

      //     // Uncomment if perp cone is hacked to be calculated
      //      printf("UE BKG per method:\n \t Perp Cone %2.2f (cl %2.2f, tr %2.2f)\n \t Phi Band %2.2f (cl %2.2f, tr %2.2f)\n \t Eta Band %2.2f (cl %2.2f, tr %2.2f)\n",
      //             perpPtSumCh+perpPtSumCh*fNeutralOverChargedRatio, perpPtSumCh*fNeutralOverChargedRatio,perpPtSumCh,
      //             phiBandPtSumChNorm+phiBandPtSumClusterNorm,phiBandPtSumClusterNorm,phiBandPtSumChNorm,
      //             etaBandPtSumChNorm+etaBandPtSumClusterNorm,etaBandPtSumClusterNorm,etaBandPtSumChNorm);
    } // UE subtraction by Eta band

    coneptsumUESubNe -= coneptsumBkgNe;
    coneptsumUESubCh -= coneptsumBkgCh;

    if ( IsEmbedingAnalysisOn() )
    {
      coneptsumUESubNeEmb -= coneptsumBkgNeEmb;
      coneptsumUESubChEmb -= coneptsumBkgChEmb;
    }

    // Calculated in case of fICMethod == kSumBkgSubEtaBandIC,
    // reset excess areas to 1 if not used in final result.
    if ( !GetIsolationCut()->IsConeExcessCorrected() )
    {
      excessAreaChEta = 1;
      excessAreaNeEta = 1;
      excessAreaNePhi = 1;
    }

    if      ( partInConeType == AliIsolationCut::kNeutralAndCharged )
      coneptsumUESub  = coneptsumUESubNe * excessAreaNePhi * excessAreaNeEta + coneptsumUESubCh * excessAreaChEta;
    else if ( partInConeType ==AliIsolationCut::kOnlyCharged        )
      coneptsumUESub  = coneptsumUESubCh * excessAreaChEta;
    else if ( partInConePhi == AliIsolationCut::kOnlyNeutral        )
      coneptsumUESub  = coneptsumUESubNe * excessAreaNePhi*excessAreaNeEta;

    if ( IsEmbedingAnalysisOn() )
    {
      if      ( partInConeType == AliIsolationCut::kNeutralAndCharged )
        coneptsumUESubEmb  = coneptsumUESubNeEmb * excessAreaNePhi * excessAreaNeEta + coneptsumUESubCh * excessAreaChEta;
      else if ( partInConeType ==AliIsolationCut::kOnlyCharged        )
        coneptsumUESubEmb  = coneptsumUESubChEmb * excessAreaChEta;
      else if ( partInConePhi == AliIsolationCut::kOnlyNeutral        )
        coneptsumUESubEmb  = coneptsumUESubNeEmb * excessAreaNePhi*excessAreaNeEta;
    }

    // Fill the histograms

    if ( !IsHighMultiplicityAnalysisOn() )
    {
      fhConeSumPtPrimMC[mcIndex]      ->Fill(photonPt, sumPtInCone, GetEventWeight()*weightPt) ;
      if ( !GetReader()->AreMCPromptPhotonsSelected() )
        fhConeSumPtPrimMC[kmcPrimPhoton]->Fill(photonPt, sumPtInCone, GetEventWeight()*weightPt) ;

      if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
      {
        fhConeSumPtChargedPrimMC[mcIndex]      ->Fill(photonPt, sumPtInConeCh, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtChargedPrimMC[kmcPrimPhoton]->Fill(photonPt, sumPtInConeCh, GetEventWeight()*weightPt) ;
      }

      if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
      {
        fhConeSumPtUESubPrimMC[mcIndex]      ->Fill(photonPt, coneptsumUESub, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtUESubPrimMC[kmcPrimPhoton]->Fill(photonPt, coneptsumUESub, GetEventWeight()*weightPt) ;

        if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtUESubChargedPrimMC[mcIndex]      ->Fill(photonPt, coneptsumUESubCh, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtUESubChargedPrimMC[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubCh, GetEventWeight()*weightPt) ;
        }
      }

      if ( perpConePtSumCh > 0 && perpConePtSumNe > 0 )
      {
        fhConeSumPtPerpConeNeutralChargedRatioPrimMC[mcIndex]      ->Fill(photonPt, perpConePtSumNe/perpConePtSumCh, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtPerpConeNeutralChargedRatioPrimMC[kmcPrimPhoton]->Fill(photonPt, perpConePtSumNe/perpConePtSumCh, GetEventWeight()*weightPt) ;
      }

      if ( sumPtInConeCh > 0 && sumPtInConeNe > 0 )
      {
        fhConeSumPtNeutralChargedRatioPrimMC[mcIndex]      ->Fill(photonPt, sumPtInConeNe/sumPtInConeCh, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtNeutralChargedRatioPrimMC[kmcPrimPhoton]->Fill(photonPt, sumPtInConeNe/sumPtInConeCh, GetEventWeight()*weightPt) ;
      }

      if ( IsEmbedingAnalysisOn() )
      {
        fhConeSumPtPrimMCEmb[mcIndex]      ->Fill(photonPt, sumPtInConeEmb, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, sumPtInConeEmb, GetEventWeight()*weightPt) ;

        if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtChargedPrimMCEmb[mcIndex]      ->Fill(photonPt, sumPtInConeChEmb, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtChargedPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, sumPtInConeChEmb, GetEventWeight()*weightPt) ;
        }

        if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
        {
          fhConeSumPtUESubPrimMCEmb[mcIndex]      ->Fill(photonPt, coneptsumUESubEmb, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtUESubPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubEmb, GetEventWeight()*weightPt) ;

          if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
          {
            fhConeSumPtUESubChargedPrimMCEmb[mcIndex]      ->Fill(photonPt, coneptsumUESubChEmb, GetEventWeight()*weightPt) ;
            if ( !GetReader()->AreMCPromptPhotonsSelected() )
              fhConeSumPtUESubChargedPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubChEmb, GetEventWeight()*weightPt) ;
          }
        }
      }
    }
    else
    {
      fhConeSumPtCenPrimMC [mcIndex]      ->Fill(photonPt, sumPtInCone, centrality, GetEventWeight()*weightPt) ;
      if ( !GetReader()->AreMCPromptPhotonsSelected() )
        fhConeSumPtCenPrimMC [kmcPrimPhoton]->Fill(photonPt, sumPtInCone, centrality, GetEventWeight()*weightPt) ;

      if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
      {
        fhConeSumPtCenChargedPrimMC [mcIndex]      ->Fill(photonPt, sumPtInConeCh, centrality, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtCenChargedPrimMC [kmcPrimPhoton]->Fill(photonPt, sumPtInConeCh, centrality, GetEventWeight()*weightPt) ;
      }

      if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
      {
        fhConeSumPtUESubCenPrimMC[mcIndex]      ->Fill(photonPt, coneptsumUESub, centrality, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtUESubCenPrimMC[kmcPrimPhoton]->Fill(photonPt, coneptsumUESub, centrality, GetEventWeight()*weightPt) ;

        if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtUESubCenChargedPrimMC[mcIndex]      ->Fill(photonPt, coneptsumUESubCh, centrality, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtUESubCenChargedPrimMC[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubCh, centrality, GetEventWeight()*weightPt) ;
        }
      }

      if ( IsEmbedingAnalysisOn() )
      {
        fhConeSumPtCenPrimMCEmb [mcIndex]      ->Fill(photonPt, sumPtInConeEmb, centrality, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhConeSumPtCenPrimMCEmb [kmcPrimPhoton]->Fill(photonPt, sumPtInConeEmb, centrality, GetEventWeight()*weightPt) ;

        if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
        {
          fhConeSumPtCenChargedPrimMCEmb [mcIndex]      ->Fill(photonPt, sumPtInConeChEmb, centrality, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtCenChargedPrimMCEmb [kmcPrimPhoton]->Fill(photonPt, sumPtInConeChEmb, centrality, GetEventWeight()*weightPt) ;
        }

        if ( isoMethod >= AliIsolationCut::kSumBkgSubIC )
        {
          fhConeSumPtUESubCenPrimMCEmb[mcIndex]      ->Fill(photonPt, coneptsumUESubEmb, centrality, GetEventWeight()*weightPt) ;
          if ( !GetReader()->AreMCPromptPhotonsSelected() )
            fhConeSumPtUESubCenPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubEmb, centrality, GetEventWeight()*weightPt) ;

          if ( partInConeType == AliIsolationCut::kNeutralAndCharged )
          {
            fhConeSumPtUESubCenChargedPrimMCEmb[mcIndex]      ->Fill(photonPt, coneptsumUESubChEmb, centrality, GetEventWeight()*weightPt) ;
            if ( !GetReader()->AreMCPromptPhotonsSelected() )
              fhConeSumPtUESubCenChargedPrimMCEmb[kmcPrimPhoton]->Fill(photonPt, coneptsumUESubChEmb, centrality, GetEventWeight()*weightPt) ;
          }
        }
      }
    }
    
    if ( !fFillOnlyTH3Histo )
    {
      if ( !GetReader()->AreMCPromptPhotonsSelected() )
      {
        fhPtPrimMC [kmcPrimPhoton]->Fill(photonPt,            GetEventWeight()*weightPt) ;
        fhEtaPrimMC[kmcPrimPhoton]->Fill(photonPt, photonEta, GetEventWeight()*weightPt) ;
        fhPhiPrimMC[kmcPrimPhoton]->Fill(photonPt, photonPhi, GetEventWeight()*weightPt) ;
      }
      
      fhPtPrimMC [mcIndex]->Fill(photonPt,            GetEventWeight()*weightPt) ;
      fhEtaPrimMC[mcIndex]->Fill(photonPt, photonEta, GetEventWeight()*weightPt) ;
      fhPhiPrimMC[mcIndex]->Fill(photonPt, photonPhi, GetEventWeight()*weightPt) ;
    }

    // In case the photon is a decay from pi0 or eta,
    // study how the decay kinematics affects the isolation
    Int_t  ndaugh   = -1;
    Bool_t okpi0    =  0, ok1     =  0, ok2     =  0;
    Int_t  pi0label = -1, d1Label = -1, d2Label = -1;
    Bool_t d2Acc   = kTRUE, overlap = kTRUE;
    Int_t  d2AbsId = -1;
    Float_t dRdaugh2 = 0, d12Angle = 0;
    
    if ( fMakePrimaryPi0DecayStudy )
    {
      if( (mcIndex == kmcPrimPi0Decay || mcIndex == kmcPrimEtaDecay ) )
      {
        if(mcIndex == kmcPrimPi0Decay) GetMCAnalysisUtils()->GetMotherWithPDG(i,111,GetMC(),okpi0, pi0label);
        else                           GetMCAnalysisUtils()->GetMotherWithPDG(i,221,GetMC(),okpi0, pi0label);
        
        if(okpi0)
        {
          ndaugh = GetMCAnalysisUtils()->GetNDaughters(pi0label,GetMC(), okpi0);
          if(ndaugh==2)
          {
            Int_t d1Pdg = 0, d1Status = 0;
            fMomDaugh1 = GetMCAnalysisUtils()->GetDaughter(0,pi0label,GetMC(),d1Pdg, d1Status,ok1, d1Label,fProdVertex);
            Int_t d2Pdg = 0, d2Status = 0;
            fMomDaugh2 = GetMCAnalysisUtils()->GetDaughter(1,pi0label,GetMC(),d2Pdg, d2Status,ok2, d2Label,fProdVertex);
            if(d2Pdg != d1Pdg && d1Pdg!=22) okpi0 = kFALSE;
            
            // Check the momentum and location of second daughter
            if(okpi0)
            {
              // assign current trigger to first daughter
              if(d1Label!=i)
              {
                Int_t tmpLabel = d2Label;
                d2Label = d1Label;
                d1Label = tmpLabel;
                fMomentum  = fMomDaugh2;
                fMomDaugh2 = fMomDaugh1;
                fMomDaugh1 = fMomentum;
              }
              
              // Check if photons hit the Calorimeter acceptance
              if(IsRealCaloAcceptanceOn() && fIsoDetector!=kCTS) // defined on base class
                d2Acc = GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector,fMomDaugh2.Eta(),
                                                                            fMomDaugh2.Theta(),fMomDaugh2.Phi(),d2AbsId) ;
              
              //printf("D2  (eta %2.2f,phi %2.2f)in real calo %d, with absId %d\n",
              //       daugh2mom.Eta(), daugh2mom.Phi()*TMath::RadToDeg(),d2Acc,d2AbsId);
              
              // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
              if(d2Acc) d2Acc = GetReader()->GetFiducialCut()->IsInFiducialCut(fMomDaugh2.Eta(),fMomDaugh2.Phi(),fIsoDetector);
              //printf("D2 fidcut %d\n",d2Acc);
              
              Float_t phiDaugh2 = fMomDaugh2.Phi();
              if(phiDaugh2 < 0) phiDaugh2+=TMath::TwoPi();
              dRdaugh2 = GetIsolationCut()->Radius(photonEta, photonPhi, fMomDaugh2.Eta(),phiDaugh2);
              
              // Opening angle, check if pairs will likely overlap
              d12Angle = fMomDaugh1.Angle(fMomDaugh2.Vect());
              if(d12Angle > overlapAngle) overlap = kFALSE;
              
            }
          }
        }
        
        //printf("Check mother of label %d: mom label %d, okmom %d ndaugh %d, daugh label1 %d, label2 %d, ok1 %d, ok2 %d, R %2.3f, opening angle %2.3f, overlap %d\n",
        //       i, pi0label,okpi0,ndaugh,d1Label,d2Label,ok1,ok2, dRdaugh2, d12Angle, overlap);
        
        if(mcIndex == kmcPrimPi0Decay)
        {
          // Second decay out of cone
          if ( dRdaugh2 > GetIsolationCut()->GetConeSize() || 
               dRdaugh2 < GetIsolationCut()->GetMinDistToTrigger() )
            fhPtPrimMCPi0DecayPairOutOfCone->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCPi0DecayPairOutOfAcceptance->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap) fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          }
          
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCPi0DecayPairNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay pt smaller than threshold
          if ( d2Acc && 
               dRdaugh2 < GetIsolationCut()->GetConeSize() && 
               dRdaugh2 > GetIsolationCut()->GetMinDistToTrigger() &&
               fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCPi0DecayPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap)
            {
              fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight()*weightPt);
            }
          }
        } // pi0 decay
        else // eta decay
        {
          // Second decay out of cone
          if ( dRdaugh2 > GetIsolationCut()->GetConeSize() || 
               dRdaugh2 < GetIsolationCut()->GetMinDistToTrigger() )
            fhPtPrimMCEtaDecayPairOutOfCone->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCEtaDecayPairOutOfAcceptance->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap) fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          }
          
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCEtaDecayPairNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay pt smaller than threshold
          if ( d2Acc && 
               dRdaugh2 < GetIsolationCut()->GetConeSize() && 
               dRdaugh2 > GetIsolationCut()->GetMinDistToTrigger() &&
               fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCEtaDecayPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap)
            {
              fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight()*weightPt);
            }
          }
        } // eta decay
        
      } // eta or pi0 decay
      
      if(overlapPi0)
      {
        if( mcIndex == kmcPrimPi0) fhPtPrimMCPi0Overlap->Fill(photonPt, GetEventWeight()*weightPt);
        if( mcIndex == kmcPrimEta) fhPtPrimMCEtaOverlap->Fill(photonPt, GetEventWeight()*weightPt);
      }
    }
    
    // Check if candidate is isolated.
    //
    // Only for isolation based on sum pT. 
    // In case of UE subtraction it needs revision!!!, 
    // implementing same approach as in reconstructed clusters.
    
    Bool_t isolated  = kFALSE;
    if ( ( isoMethod == AliIsolationCut::kSumPtIC || 
           isoMethod >= AliIsolationCut::kSumBkgSubIC ) &&
         ( sumPtInCone < GetIsolationCut()->GetSumPtThreshold() || 
           sumPtInCone > GetIsolationCut()->GetSumPtThresholdMax() ) )
      isolated = kTRUE;
    
    //printf("Primary in cone n Part %d sum pT in cone %2.2f\n",npart, sumPtInCone);
    
    if ( GetIsolationCut()->GetICMethod() == AliIsolationCut::kPtThresIC && npart == 0 )
      isolated = kTRUE;
    
    if ( isolated )
    {
      if ( !fFillOnlyTH3Histo )
      {
        fhPtPrimMCiso [mcIndex]      ->Fill(photonPt, GetEventWeight()*weightPt) ;
        if ( !GetReader()->AreMCPromptPhotonsSelected() )
          fhPtPrimMCiso [kmcPrimPhoton]->Fill(photonPt, GetEventWeight()*weightPt) ;
      }
      
      if ( fMakePrimaryPi0DecayStudy )
      {
        if ( mcIndex == kmcPrimPi0Decay )
        {
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCPi0DecayIsoPairNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of cone
          if ( dRdaugh2 > GetIsolationCut()->GetConeSize() || 
               dRdaugh2 < GetIsolationCut()->GetMinDistToTrigger() )
            fhPtPrimMCPi0DecayIsoPairOutOfCone->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCPi0DecayIsoPairOutOfAcceptance->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap) fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          }
          
          // Second decay pt smaller than threshold
          if ( d2Acc && 
               dRdaugh2 < GetIsolationCut()->GetConeSize() &&
               dRdaugh2 > GetIsolationCut()->GetMinDistToTrigger() &&
               fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold() )
          {
            fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap)
            {
              fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight()*weightPt);
            }
          }
        }// pi0 decay
        else if( mcIndex == kmcPrimEtaDecay )
        {
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCEtaDecayIsoPairNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of cone
          if ( dRdaugh2 > GetIsolationCut()->GetConeSize() || 
               dRdaugh2 < GetIsolationCut()->GetMinDistToTrigger() )
            fhPtPrimMCEtaDecayIsoPairOutOfCone->Fill(photonPt, GetEventWeight()*weightPt);
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCEtaDecayIsoPairOutOfAcceptance->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap) fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          }
          
          // Second decay pt smaller than threshold
          if ( d2Acc && 
               dRdaugh2 < GetIsolationCut()->GetConeSize() &&
               dRdaugh2 > GetIsolationCut()->GetMinDistToTrigger() &&
               fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold() )
          {
            fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight()*weightPt);
            if(!overlap)
            {
              fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight()*weightPt);
            }
          }
        }// eta decay
        
        if(overlapPi0)
        {
          if( mcIndex == kmcPrimPi0 ) fhPtPrimMCPi0IsoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
          if( mcIndex == kmcPrimEta ) fhPtPrimMCEtaIsoOverlap->Fill(photonPt, GetEventWeight()*weightPt);
        }
      }
    } // isolated
    
  }//loop on primaries
  
  AliDebug(1,"End");
  
}

//_____________________________________________________________
/// Print some relevant parameters set for the analysis.
//_____________________________________________________________
void AliAnaParticleIsolation::Print(const Option_t * opt) const
{
  if (! opt )
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Calorimeter for isolation = %s (%d) \n",  GetCalorimeterString().Data(),fIsoDetector) ;
  printf("Detector for candidate isolation = %s \n", fIsoDetectorString.Data()) ;
  printf("Select leading cand. %d, within neutrals %d\n",fLeadingOnly,fCheckLeadingWithNeutralClusters);
  printf("Fill tagged decays histo %d, n decay bits %d\n",
         fFillTaggedDecayHistograms,fNDecayBits);
  
  printf("Active histogram filled: TM %d, SS %d, per SM %d, per TCard %d,\n"
         "EMCal region %d, NLM %d, TH3 only %d, AliIsolationCut histo %d\n",
         fFillTMHisto, fFillSSHisto, fFillPerSMHistograms, fFillPerTCardIndexHistograms, 
         fFillEMCALRegionHistograms, fFillNLMHistograms, fFillOnlyTH3Histo, fFillIsolationControlHistograms);
  
  printf("Studies: Tracks in cone %d, Conversion radius %d;\n"
         " pt in cone cuts %d (n=%d); eta cuts %d (n=%d);\n"
         " r cuts %d (n=%d); n cell cuts %d (n=%d); exotic cuts %d (n=%d)\n ",
         fStudyTracksInCone,fStudyMCConversionRadius, 
         fStudyPtCutInCone,fNPtCutsInCone,
         fStudyEtaCutInCone,fNEtaCutsInCone,
         fStudyRCutInCone,fNRCutsInCone,
         fStudyNCellsCut,fNNCellsInCandidate,
         fStudyExoticTrigger,fNExoCutInCandidate);
  
  if ( IsDataMC() ) 
    printf("MC analysis: select prim %d, prim pi0 decay %d;\n histo active: Overlaps %d, primary overlap angle %2.2f\n",
           fSelectPrimariesInCone,fMakePrimaryPi0DecayStudy,
           fFillOverlapHistograms,fMinCellsAngleOverlap);
  
  printf("    \n") ;
} 

//_____________________________________________________________
/// Set the detector for the analysis.
//_____________________________________________________________
void AliAnaParticleIsolation::SetTriggerDetector(TString det)
{
  fIsoDetectorString = det;
  
  if     (det=="EMCAL") fIsoDetector = kEMCAL;
  else if(det=="PHOS" ) fIsoDetector = kPHOS;
  else if(det=="CTS")   fIsoDetector = kCTS;
  else if(det=="DCAL")  fIsoDetector = kDCAL;
  else if(det.Contains("DCAL") && det.Contains("PHOS")) fIsoDetector = kDCALPHOS;
  else AliFatal(Form("Detector < %s > not known!", det.Data()));
}

//_________________________________________________________
/// Set the detector for the analysis.
//_________________________________________________________
void AliAnaParticleIsolation::SetTriggerDetector(Int_t det)
{
  fIsoDetector = det;
  
  if     (det==kEMCAL)    fIsoDetectorString = "EMCAL";
  else if(det==kPHOS )    fIsoDetectorString = "PHOS";
  else if(det==kCTS)      fIsoDetectorString = "CTS";
  else if(det==kDCAL)     fIsoDetectorString = "DCAL";
  else if(det==kDCALPHOS) fIsoDetectorString = "DCAL_PHOS";
  else AliFatal(Form("Detector < %d > not known!", det));
}

//_________________________________________________________
/// Check shower shape and cluster/trac sum pT in cone 
/// in different EMCal regions. Also check the cluster eta-phi
/// distribution in shower shape tail.
///
/// \param pt: cluster pT
/// \param phi: cluster phi
/// \param eta: cluster eta
/// \param m02: shower shape long
/// \param coneptsumTrack: sum of tracks pT in cone
/// \param coneptsumCluster: sum of clusters pT in cone
/// \param isolated: bool
/// \param iSM: super module number
/// 
//_________________________________________________________
void AliAnaParticleIsolation::StudyEMCALRegions
(Float_t pt,  Float_t phi, Float_t eta, Float_t m02, 
 Float_t coneptsumTrack, Float_t coneptsumCluster, 
 Bool_t isolated, Int_t iSM) 
{  
  if ( !fCluster ) return;
  
  // Get the predefined regions
  Int_t etaRegion = -1, phiRegion = -1;
  GetCaloUtils()->GetEMCALSubregion(fCluster,GetReader()->GetEMCALCells(),etaRegion,phiRegion);
  
  if ( etaRegion >= 0 && etaRegion < 4 && phiRegion >=0 && phiRegion < 3 ) 
  {
    fhLam0EMCALRegionPerSM            [isolated][etaRegion][phiRegion][iSM]->Fill(pt, m02, GetEventWeight());
    fhConeSumPtTrackEMCALRegionPerSM            [etaRegion][phiRegion][iSM]->Fill(pt, coneptsumTrack  , GetEventWeight());
    fhConeSumPtClusterEMCALRegionPerSM          [etaRegion][phiRegion][iSM]->Fill(pt, coneptsumCluster, GetEventWeight());
  }
  
  // Check eta-phi distribution for shower shape tail region
  if ( m02 >=0.3 && m02 <= 0.4 )
  {
    Float_t ptLimit[] = {2,3,4,5,6,8,10,12};
    Int_t ptbin = -1;
    for(Int_t ipt = 0; ipt < 7; ipt++)  
    {
      if( pt >= ptLimit[ipt] && pt < ptLimit[ipt+1]  )
      {
        ptbin = ipt;
        break;
      }
    }
    
    if ( ptbin >= 0 ) 
      fhEtaPhiLam0BinPtBin[isolated][ptbin]->Fill(eta, phi, GetEventWeight());
  }
}

//_________________________________________________________
/// Check shower shape for different conversion radius
///
/// \param pt: cluster pT
/// \param isolated: bool
/// \param iSM: supermodule number
/// \param m02: shower shape long
/// \param mcTag: mc particle type tag
/// \param label: cluster mc label
/// 
//_________________________________________________________
void AliAnaParticleIsolation::StudyMCConversionRadius
(Float_t pt, Bool_t isolated, Int_t iSM, 
 Float_t m02, Int_t mcTag, Int_t label) 
{
  // Do it for converted photon clusters
  Bool_t convPhoton = kFALSE;
  if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
      GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)     &&
     !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)        &&
     !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)           )
    convPhoton = kTRUE;
  
  if ( !convPhoton ) return;
  
  // Check that it has daughters
  Int_t pdg  = 0, status  = 0, momLabel  = -1;
  Int_t pdgD = 0, statusD = 0, daugLabel = -1;
  Bool_t ok = kFALSE, okD = kFALSE;
  
  //fPrimaryMom = 
  GetMCAnalysisUtils()->GetMother(label,GetMC(), pdg, status, ok, momLabel);     
  //fMomentum = 
  GetMCAnalysisUtils()->GetDaughter(0,momLabel,GetMC(),pdgD, statusD, okD, daugLabel, fProdVertex);
  
  if(!okD) return;
  
  Float_t prodR = TMath::Sqrt(fProdVertex.X()*fProdVertex.X()+fProdVertex.Y()*fProdVertex.Y());
  
  //printf("Iso %d, Conversion: mom pdg %d (stat %d), 1st daugher %d (stat %d), mom label %d, org label %d, daugh label %d, prodR %f\n", isolated, pdg,status, pdgD, statusD, 
  //       momLabel, aod->GetLabel(),daugLabel,prodR);
  
  fhMCConversionVertex[isolated]->Fill(pt,prodR,GetEventWeight());
  if(GetCalorimeter() == kEMCAL && GetFirstSMCoveredByTRD() >= 0 && iSM >= GetFirstSMCoveredByTRD() )
    fhMCConversionVertexTRD[isolated]->Fill(pt,prodR,GetEventWeight());
  
  Int_t convR = -1;
  if      ( prodR < 75.  ) convR = 0;
  else if ( prodR < 275. ) convR = 1;
  else if ( prodR < 375. ) convR = 2;
  else if ( prodR < 400. ) convR = 3;
  else if ( prodR < 430. ) convR = 4;
  else                     convR = 5;
  
  if ( convR >= 0 )
  {
    fhMCConversionLambda0Rcut[convR][isolated]->Fill(pt,m02,GetEventWeight());
    if ( GetCalorimeter() == kEMCAL && GetFirstSMCoveredByTRD() >= 0 && iSM >= GetFirstSMCoveredByTRD() )
      fhMCConversionLambda0Rcut[convR][isolated]->Fill(pt,m02,GetEventWeight());
    
    //
    // EMCAL SM regions
    //
    //            if ( GetCalorimeter() == kEMCAL && fFillEMCALRegionHistograms )
    //            {              
    //              Int_t etaRegion = -1, phiRegion = -1;
    //              
    //              if ( cluster ) GetCaloUtils()->GetEMCALSubregion(cluster,GetReader()->GetEMCALCells(),etaRegion,phiRegion);
    //              
    //              if( etaRegion >= 0 && etaRegion < 4 && phiRegion >=0 && phiRegion < 3 ) 
    //              {
    //                fhLam0EMCALRegionMCConvRcut[isolated][etaRegion][phiRegion][convR]->Fill(pt,m02, GetEventWeight());
    //                
    //                if ( GetFirstSMCoveredByTRD() >= 0 && iSM >= GetFirstSMCoveredByTRD()  )
    //                  fhLam0EMCALRegionTRDMCConvRcut[isolated][etaRegion][phiRegion][convR]->Fill(pt, m02, GetEventWeight());
    //                
    //              } // region found
    //            } // check region
    
  } // convR > 0
} 

//______________________________________________________________________________________________________________
/// Get the cluster pT or sum of pT in isolation cone.
//______________________________________________________________________________________________________________
void AliAnaParticleIsolation::StudyClustersInCone(AliCaloTrackParticleCorrelation * aodParticle)
{  
  if ( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t coneptLeadCluster = aodParticle->GetNeutralLeadPtInCone();
  Float_t coneptsumCluster  = aodParticle->GetNeutralPtSumInCone();
  Float_t ptTrig            = aodParticle->Pt();
  Float_t weightTrig        = aodParticle->GetWeight();
  
  // Recover reference arrays with clusters and tracks
  TObjArray * refclusters = aodParticle->GetObjArray(GetAODObjArrayName()+"Clusters");
  if ( !refclusters )
  {
    if ( fStudyExoticTrigger && fIsExoticTrigger )
      fhConeSumPtClusterExoTrigger->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
    
    if ( fFillPerSMHistograms )     
      fhConeSumPtClusterPerSM[aodParticle->GetSModNumber()]->Fill(ptTrig,0., GetEventWeight()*weightTrig);
    
    if ( fFillPerTCardIndexHistograms )     
      fhConeSumPtClusterPerTCardIndex[fTCardIndex]->Fill(ptTrig,0., GetEventWeight()*weightTrig);
    
    if ( coneptLeadCluster > 0  || coneptsumCluster > 0 ) 
      AliError(Form("No ref tracks!!! sum %f, lead %f",coneptsumCluster,coneptLeadCluster));
    
    return ;
  }
  
  // Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t ptcone = 0;
  Float_t coneptsumClusterPerMaxCut[20];
  Float_t coneptsumClusterPerRCut  [10];
  
  Float_t coneptsumClusterPerNCellCut[20];
  Float_t coneptsumClusterPerExoCut  [20];
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++)
    {
      fConeNClusterPerMinCut    [icut] = 0;
      fConeptsumClusterPerMinCut[icut] = 0;
      coneptsumClusterPerMaxCut [icut] = 0;
    }
  }
  
  if ( fStudyRCutInCone )
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      coneptsumClusterPerRCut[icut] = 0;
    }
  }
  
  Int_t ishsh = -1;
  if ( fStudyNCellsCut )
  {
    Float_t m02 = aodParticle->GetM02();
    if      ( m02 > 0.1 && m02 <= 0.3 ) ishsh = 0;
    else if ( m02 > 0.3 && m02 <= 0.4 ) ishsh = 1;  
    else if ( m02 > 0.4 && m02 <= 1.0 ) ishsh = 2;  
    else if ( m02 > 1.0 && m02 <= 3.0 ) ishsh = 3;  
    
    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      coneptsumClusterPerNCellCut[icut] = 0;
    }
  }
  
  if ( fStudyExoticTrigger )
  {
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      coneptsumClusterPerExoCut[icut] = 0;
    }
  }
  
  for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
  {
    AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
    
    calo->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line
    
    ptcone = fMomentum.Pt();
    
    if ( fFillPerSMHistograms ) 
    {
      fhPtInConePerSM       [aodParticle->GetSModNumber()]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      fhPtClusterInConePerSM[aodParticle->GetSModNumber()]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
    }
    
    if ( fFillPerTCardIndexHistograms ) 
    {
      fhPtInConePerTCardIndex       [fTCardIndex]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      fhPtClusterInConePerTCardIndex[fTCardIndex]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
    }
    
    if ( fStudyExoticTrigger && fIsExoticTrigger )
    {
      fhPtInConeExoTrigger        ->Fill(ptTrig , ptcone, GetEventWeight()*weightTrig);
      fhPtClusterInConeExoTrigger ->Fill(ptTrig , ptcone, GetEventWeight()*weightTrig);
    }
    
    if ( IsPileUpAnalysisOn() )
    {
      if(GetReader()->IsPileUpFromSPD())               fhPtInConePileUp[0]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig, ptcone, GetEventWeight()*weightTrig);
    }
    
    if ( fStudyPtCutInCone )
    {
      for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
      {
        if ( ptcone > fMinPtCutInCone[icut] ) 
        {
          fConeptsumClusterPerMinCut[icut]+=ptcone;
          fConeNClusterPerMinCut    [icut]++;
        }
        
        if ( ptcone < fMaxPtCutInCone[icut] ) coneptsumClusterPerMaxCut[icut]+=ptcone;
      }
    }
    
    if ( fStudyRCutInCone )
    {
      Float_t distance = GetIsolationCut()->Radius(aodParticle->Eta(), aodParticle->Phi(), fMomentum.Eta(), GetPhi(fMomentum.Phi()));
      for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
      {
        if ( distance < fRCutInCone[icut] ) 
        {
          coneptsumClusterPerRCut[icut]+=ptcone;
          fhPtClusterInConePerRCut->Fill(icut+1, ptcone, GetEventWeight()*weightTrig);
        }
      }
    }
    
    if ( fStudyNCellsCut )
    {
      if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
        fhPtClusterInConePerNCellPerSM[ishsh]->Fill(ptcone, fTrigSupMod, fNCellsWithWeight);
      
      for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
      {
        if ( fNCellsWithWeight >= fNCellsInCandidate[icut] ) 
        {
          coneptsumClusterPerNCellCut[icut]+=ptcone;
          fhPtClusterInConePerNCellCut->Fill(icut+1, ptcone, GetEventWeight()*weightTrig);
        }
      }
    }
    
    if ( fStudyExoticTrigger )
    {
      for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
      {
        if ( fClusterExoticity < fExoCutInCandidate[icut] ) 
        {
          coneptsumClusterPerExoCut[icut]+=ptcone;
          fhPtClusterInConePerExoCut->Fill(icut+1, ptcone, GetEventWeight()*weightTrig);
        }
      }
    }
  }
  
  if ( fStudyExoticTrigger && fIsExoticTrigger )
    fhConeSumPtClusterExoTrigger  ->Fill(ptTrig, coneptsumCluster  , GetEventWeight()*weightTrig);
  
  if ( fFillPerSMHistograms )     
    fhConeSumPtClusterPerSM[aodParticle->GetSModNumber()]->Fill(ptTrig,coneptsumCluster, GetEventWeight()*weightTrig);
  
  if ( fFillPerTCardIndexHistograms )     
    fhConeSumPtClusterPerTCardIndex[fTCardIndex]->Fill(ptTrig,coneptsumCluster, GetEventWeight()*weightTrig);
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeNClusterPerMinPtCut    ->Fill(icut, fConeNClusterPerMinCut    [icut], GetEventWeight()*weightTrig);
        fhConeSumPtClusterPerMinPtCut->Fill(icut, fConeptsumClusterPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeSumPtClusterPerMaxPtCut->Fill(icut, coneptsumClusterPerMaxCut [icut], GetEventWeight()*weightTrig);
      }
      else
      {
        fhConeNClusterPerMinPtCutCent    ->Fill(icut, fConeNClusterPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeSumPtClusterPerMinPtCutCent->Fill(icut, fConeptsumClusterPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
      }
    }
  }
  
  if ( fStudyRCutInCone )
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      fhConeSumPtClusterPerRCut->Fill(icut, coneptsumClusterPerRCut[icut], GetEventWeight()*weightTrig);
    }
  }
  
  if ( fStudyNCellsCut )
  {     
    if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
      fhConeSumPtClusterPerNCellPerSM[ishsh]->Fill(coneptsumCluster, fTrigSupMod, fNCellsWithWeight);
    
    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      fhConeSumPtClusterPerNCellCut->Fill(icut, coneptsumClusterPerNCellCut[icut], GetEventWeight()*weightTrig);
    }
  }
  
  if ( fStudyExoticTrigger )
  { 
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      fhConeSumPtClusterPerExoCut->Fill(icut, coneptsumClusterPerExoCut[icut], GetEventWeight()*weightTrig);
    }
  }  
}

//________________________________________________________________________________________________
/// Get the cluster pT or sum of pT in eta-phi bands.
/// Fill additional histograms not done in AliIsolationCut
/// Right now only for the case fStudyPtInCone, but other possibilities to be considered
//________________________________________________________________________________________________
void AliAnaParticleIsolation::StudyClustersUEInCone(AliCaloTrackParticleCorrelation * pCandidate)
{
  if( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kOnlyCharged ) return ;
  
  // Select the Calorimeter of the photon
  TObjArray * caloList = 0x0;
  if      (GetCalorimeter() == kPHOS )
    caloList = GetPHOSClusters();
  else if (GetCalorimeter() == kEMCAL)
    caloList    = GetEMCALClusters();
  
  if ( !caloList )
  {
    AliWarning(Form("TObjArray with %s clusters is NULL!",GetCalorimeterString().Data()));
    return;
  }
  
  // Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t conesize    = GetIsolationCut()->GetConeSize();
  Float_t conesizegap = GetIsolationCut()->GetConeSizeBandGap();
  Float_t emcEtaSize  = GetIsolationCut()->GetEMCEtaSize();
  Float_t emcPhiSize  = GetIsolationCut()->GetEMCPhiSize();
  Float_t distToTrig  = GetIsolationCut()->GetMinDistToTrigger();

  // Bands normalization
  //
  Float_t coneArea = conesize*conesize*TMath::Pi(); // Area = pi R^2, isolation cone area
  Float_t coneAreaGap = (conesize+conesizegap)*(conesize+conesizegap)*TMath::Pi(); // Area = pi R^2, isolation cone+gap area

  if ( distToTrig > 0 )
  {
    coneArea    -= distToTrig*distToTrig*TMath::Pi();
    coneAreaGap -= distToTrig*distToTrig*TMath::Pi();
  }
  
  // Area of band, rectangle minus isolation region
  Float_t etaBandArea = 2*(conesize+conesizegap) * emcEtaSize - coneAreaGap;
  Float_t phiBandArea = 2*(conesize+conesizegap) * emcPhiSize - coneAreaGap;
  //printf("Area band eta %f phi %f\n",etaBandArea,phiBandArea);
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fConeptsumEtaBandClusterPerMinCut   [icut] = 0;
      fConeNEtaBandClusterPerMinCut       [icut] = 0;
      fConeptsumClusterSubEtaBandPerMinCut[icut] = 0;
      fConeNClusterSubEtaBandPerMinCut    [icut] = 0;
      
      fConeptsumPhiBandClusterPerMinCut   [icut] = 0;
      fConeNPhiBandClusterPerMinCut       [icut] = 0;
      fConeptsumClusterSubPhiBandPerMinCut[icut] = 0;
      fConeNClusterSubPhiBandPerMinCut    [icut] = 0;
    }
  }
  
  //Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  Float_t weightTrig= pCandidate->GetWeight();
  
  Float_t phiBandPtSumCluster = 0;
  Float_t etaBandPtSumCluster = 0;
  
  for(Int_t icalo=0; icalo < caloList->GetEntriesFast(); icalo++)
   {
     AliVCluster* calo = (AliVCluster *) caloList->At(icalo);
     
    if(!calo)
    {
      AliWarning("Cluster not available?");
      continue;
    }
    
     // Do not count the candidate (photon or pi0) or the daughters of the candidate
     if ( calo->GetID() == pCandidate->GetCaloLabel(0) ||
          calo->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;
     
     // Skip matched clusters with tracks in case of neutral+charged analysis
     if ( GetIsolationCut()->IsTrackMatchedClusterRejectionInConeOn() )
     {
       if( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kNeutralAndCharged &&
           IsTrackMatched(calo, GetReader()->GetInputEvent()) ) continue ;
     }
     
     calo->GetMomentum(fMomentum,vertex) ; //Assume that come from vertex in straight line
     
     Float_t ptCluster  = fMomentum.Pt();
     Float_t etaCluster = fMomentum.Eta();
     Float_t phiCluster = GetPhi(fMomentum.Phi());

     // Angles between trigger and track
     Double_t dEta = etaTrig - etaCluster;
     Double_t dPhi = phiTrig - phiCluster;

     // Shift phi angle when trigger is close to 0 or 360
     if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
     if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();
     
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Eta or phi band with respect the trigger cone
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, etaCluster, phiCluster);

     if ( rad > conesize+conesizegap )
    {
      // Phi band
      //
      Bool_t takeIt = kTRUE;

      // Look only half TPC with respect candidate, avoid opposite side jet 
      if ( TMath::Abs(phiTrig-phiCluster) > TMath::PiOver2() )  takeIt = kFALSE;
      
      // Within eta cone size
      if ( TMath::Abs(dEta) < conesize+conesizegap &&  takeIt ) 
      {
        phiBandPtSumCluster += ptCluster;
        //printf("cluster %d, pT %f in phi band\n",icalo,ptCluster);

        if ( fStudyPtCutInCone )
        {
          for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
          {
            if ( ptCluster > fMinPtCutInCone[icut] ) 
            {
              fConeptsumPhiBandClusterPerMinCut[icut]+=ptCluster;
              fConeNPhiBandClusterPerMinCut    [icut]++;
            }          
          }
        } // min pt
      } // phi band
      
      // Eta band
      //
      //takeIt = kTRUE;
      
      // Within phi cone size
      if ( TMath::Abs(dPhi) < conesize+conesizegap )// && takeIt )
      {
        etaBandPtSumCluster += ptCluster;
        //printf("cluster %d, pT %f in phi band\n",icalo,ptCluster);

        if ( fStudyPtCutInCone )
        {
          for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
          {
            if ( ptCluster > fMinPtCutInCone[icut] ) 
            {
              fConeptsumEtaBandClusterPerMinCut[icut]+=ptCluster;
              fConeNEtaBandClusterPerMinCut    [icut]++;
            }          
          }
        } // min pt
       
      } // eta band
      
    } // out of cone
   
  } // track loop
  
  if ( fStudyPtCutInCone )
  {    
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      // Normalize to cone area
      // For simplicity do not correct for excess of detector acceptance
      //
//      if (icut == 0 )printf("cut %d before N eta %f phi %f, sum eta %f phi %f\n",icut,
//             fConeptsumEtaBandClusterPerMinCut[icut], fConeptsumPhiBandClusterPerMinCut[icut],
//             fConeNEtaBandClusterPerMinCut    [icut], fConeNPhiBandClusterPerMinCut    [icut]);
      
      if ( etaBandArea > 0 )
      {
        fConeptsumEtaBandClusterPerMinCut[icut] *= ( coneArea / etaBandArea );
        fConeNEtaBandClusterPerMinCut    [icut] *= ( coneArea / etaBandArea );
      }
      
      if ( phiBandArea > 0 )
      {
        fConeptsumPhiBandClusterPerMinCut[icut] *= ( coneArea / phiBandArea );
        fConeNPhiBandClusterPerMinCut    [icut] *= ( coneArea / phiBandArea );
      }
//      
//     if (icut == 0) printf("cut %d after N eta %f phi %f, sum eta %f phi %f\n",icut,
//                fConeptsumEtaBandClusterPerMinCut[icut], fConeptsumPhiBandClusterPerMinCut[icut],
//                fConeNEtaBandClusterPerMinCut    [icut], fConeNPhiBandClusterPerMinCut    [icut]);
      
      // Subtract to content in trigger cone (it must be filled before on another method StudyClustersInCone())
      //
      fConeptsumClusterSubEtaBandPerMinCut[icut] = fConeptsumClusterPerMinCut[icut] - fConeptsumEtaBandClusterPerMinCut[icut];
      fConeNClusterSubEtaBandPerMinCut    [icut] = fConeNClusterPerMinCut    [icut] - fConeNEtaBandClusterPerMinCut    [icut];
      
      fConeptsumClusterSubPhiBandPerMinCut[icut] = fConeptsumClusterPerMinCut[icut] - fConeptsumPhiBandClusterPerMinCut[icut];
      fConeNClusterSubPhiBandPerMinCut    [icut] = fConeNClusterPerMinCut    [icut] - fConeNPhiBandClusterPerMinCut    [icut];
      
      // Fill histograms
      //      
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhEtaBandConeSumPtClusterPerMinPtCut->Fill(icut+1, fConeptsumEtaBandClusterPerMinCut[icut], GetEventWeight()*weightTrig);
        fhEtaBandConeNClusterPerMinPtCut    ->Fill(icut+1, fConeNEtaBandClusterPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhConeSumPtClusterSubEtaBandPerMinPtCut->Fill(icut+1, fConeptsumClusterSubEtaBandPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeNClusterSubEtaBandPerMinPtCut    ->Fill(icut+1, fConeNClusterSubEtaBandPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhPhiBandConeSumPtClusterPerMinPtCut->Fill(icut+1, fConeptsumPhiBandClusterPerMinCut[icut], GetEventWeight()*weightTrig);
        fhPhiBandConeNClusterPerMinPtCut    ->Fill(icut+1, fConeNPhiBandClusterPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhConeSumPtClusterSubPhiBandPerMinPtCut->Fill(icut+1, fConeptsumClusterSubPhiBandPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeNClusterSubPhiBandPerMinPtCut    ->Fill(icut+1, fConeNClusterSubPhiBandPerMinCut    [icut], GetEventWeight()*weightTrig);
      }
      else
      {
        fhEtaBandConeNClusterPerMinPtCutCent    ->Fill(icut, fConeNEtaBandClusterPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhEtaBandConeSumPtClusterPerMinPtCutCent->Fill(icut, fConeptsumEtaBandClusterPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhConeSumPtClusterSubEtaBandPerMinPtCutCent->Fill(icut+1, fConeptsumClusterSubEtaBandPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeNClusterSubEtaBandPerMinPtCutCent    ->Fill(icut+1, fConeNClusterSubEtaBandPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhPhiBandConeNClusterPerMinPtCutCent    ->Fill(icut, fConeNPhiBandClusterPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhPhiBandConeSumPtClusterPerMinPtCutCent->Fill(icut, fConeptsumPhiBandClusterPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhConeSumPtClusterSubPhiBandPerMinPtCutCent->Fill(icut+1, fConeptsumClusterSubPhiBandPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeNClusterSubPhiBandPerMinPtCutCent    ->Fill(icut+1, fConeNClusterSubPhiBandPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
      }
      
      // Sum charged and neutral pT
      // Do here ONLY if StudyClustersUEInCone is called after StudyTracksUEUnCone!!!
      if ( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kNeutralAndCharged )
      {
        if ( !IsHighMultiplicityAnalysisOn() )
        {
          fhConeSumPtSubEtaBandPerMinPtCut->Fill(icut+1, 
                                                 fConeptsumClusterSubEtaBandPerMinCut[icut]+fConeptsumTrackSubEtaBandPerMinCut[icut],
                                                 GetEventWeight()*weightTrig);
          fhConeNSubEtaBandPerMinPtCut    ->Fill(icut+1, 
                                                 fConeNClusterSubEtaBandPerMinCut    [icut]+fConeNTrackSubEtaBandPerMinCut    [icut],
                                                 GetEventWeight()*weightTrig);
          
          fhConeSumPtSubPhiBandPerMinPtCut->Fill(icut+1, 
                                                 fConeptsumClusterSubPhiBandPerMinCut[icut]+fConeptsumTrackSubPhiBandPerMinCut[icut], 
                                                 GetEventWeight()*weightTrig);
          fhConeNSubPhiBandPerMinPtCut    ->Fill(icut+1, 
                                                 fConeNClusterSubPhiBandPerMinCut    [icut]+fConeNTrackSubPhiBandPerMinCut    [icut], 
                                                 GetEventWeight()*weightTrig);
        }
        else
        {
          fhConeSumPtSubEtaBandPerMinPtCutCent->Fill(icut+1, 
                                                 fConeptsumClusterSubEtaBandPerMinCut[icut]+fConeptsumTrackSubEtaBandPerMinCut[icut],
                                                 GetEventCentrality(), GetEventWeight()*weightTrig);
          fhConeNSubEtaBandPerMinPtCutCent    ->Fill(icut+1, 
                                                 fConeNClusterSubEtaBandPerMinCut    [icut]+fConeNTrackSubEtaBandPerMinCut    [icut],
                                                 GetEventCentrality(), GetEventWeight()*weightTrig);
          
          fhConeSumPtSubPhiBandPerMinPtCutCent->Fill(icut+1, 
                                                 fConeptsumClusterSubPhiBandPerMinCut[icut]+fConeptsumTrackSubPhiBandPerMinCut[icut], 
                                                 GetEventCentrality(), GetEventWeight()*weightTrig);
          fhConeNSubPhiBandPerMinPtCutCent    ->Fill(icut+1, 
                                                 fConeNClusterSubPhiBandPerMinCut    [icut]+fConeNTrackSubPhiBandPerMinCut    [icut], 
                                                 GetEventCentrality(), GetEventWeight()*weightTrig);
        }
      } // neutral and charged in cone, sum both
    }
  } // fStudyPtInCone
  
}


//___________________________________________________________________________________________________________
/// Get the track pT or sum of pT in isolation cone.
//___________________________________________________________________________________________________________
void AliAnaParticleIsolation::StudyTracksInCone(AliCaloTrackParticleCorrelation * aodParticle)
{  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  Float_t coneptLeadTrack = aodParticle->GetChargedLeadPtInCone();
  Float_t coneptsumTrack  = aodParticle->GetChargedPtSumInCone();
  Float_t ptTrig          = aodParticle->Pt();
  Float_t weightTrig      = aodParticle->GetWeight();
  
  // Recover reference arrays with clusters and tracks
  TObjArray * reftracks   = aodParticle->GetObjArray(GetAODObjArrayName()+"Tracks");
  if ( !reftracks )
  {    
    if ( fStudyExoticTrigger && fIsExoticTrigger )
      fhConeSumPtTrackExoTrigger->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
    
    if ( fFillPerSMHistograms )     
      fhConeSumPtTrackPerSM[aodParticle->GetSModNumber()]->Fill(ptTrig,0., GetEventWeight()*weightTrig);
    
    if ( fFillPerTCardIndexHistograms )     
      fhConeSumPtTrackPerTCardIndex[fTCardIndex]->Fill(ptTrig,0., GetEventWeight()*weightTrig);
    
    if ( fStudyTracksInCone )
    {
      fhConeSumPtTrackTOFNo ->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackTOFBC0->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackTOFBCN->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackITSRefitOnSPDOn  ->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackITSRefitOffSPDOff->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackITSRefitOnSPDOff ->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
      fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, 0., GetEventWeight()*weightTrig);
    }
    
    if ( coneptLeadTrack > 0  || coneptsumTrack > 0 ) 
      AliError(Form("No ref tracks!!! sum %f, lead %f",coneptsumTrack,coneptLeadTrack));
    
    return ;
  }
  
  Double_t bz = GetReader()->GetInputEvent()->GetMagneticField();
  
  Float_t pTtrack  = 0;
  Float_t phitrack = 0;
  Float_t etatrack = 0;

  Float_t coneptsumTrackPerMaxCut[20];
  Float_t coneptsumTrackPerEtaCut[10];
  Float_t coneptsumTrackPerNCellCut[20];
  Float_t coneptsumTrackPerExoCut[20];
  Float_t coneptsumTrackPerRCut  [10];
  Float_t coneptsumTrackTOFBC0 = 0;
  Float_t coneptsumTrackTOFBCN = 0;
  Float_t coneptsumTrackTOFNo  = 0;
  Float_t coneptsumTrackITSRefitOnSPDOn   = 0;
  Float_t coneptsumTrackITSRefitOnSPDOff  = 0;
  Float_t coneptsumTrackITSRefitOffSPDOff = 0;
  Float_t coneptsumTrackTOFBC0ITSRefitOnSPDOn = 0;
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fConeptsumTrackPerMinCut[icut] = 0;
      coneptsumTrackPerMaxCut [icut] = 0;
      fConeNTrackPerMinCut    [icut] = 0;
    }
  }
  
  if ( fStudyEtaCutInCone )
  {
    for(Int_t icut = 0; icut < fNEtaCutsInCone; icut++) 
    {
      coneptsumTrackPerEtaCut[icut] = 0;
    }
  }
  
  if ( fStudyRCutInCone )
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      coneptsumTrackPerRCut[icut] = 0;
    }
  }
  
  Int_t ishsh = -1;
  if ( fStudyNCellsCut )
  {
    Float_t m02 = aodParticle->GetM02();
    if      ( m02 > 0.1 && m02 <= 0.3 ) ishsh = 0;
    else if ( m02 > 0.3 && m02 <= 0.4 ) ishsh = 1;  
    else if ( m02 > 0.4 && m02 <= 1.0 ) ishsh = 2;  
    else if ( m02 > 1.0 && m02 <= 3.0 ) ishsh = 3;  
    
    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      coneptsumTrackPerNCellCut[icut] = 0;
    }
  }
  
  if ( fStudyExoticTrigger )
  { 
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      coneptsumTrackPerExoCut[icut] = 0;
    }
  }
  
  for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
  {
    AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
    
    pTtrack  = track->Pt();
    
    if(fFillPerSMHistograms)   
    {
      fhPtInConePerSM     [aodParticle->GetSModNumber()]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      fhPtTrackInConePerSM[aodParticle->GetSModNumber()]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
    }
    
    if(fFillPerTCardIndexHistograms)   
    {
      fhPtInConePerTCardIndex     [fTCardIndex]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      fhPtTrackInConePerTCardIndex[fTCardIndex]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
    }
    
    if ( IsDataMC() && GetMC() )
    {
      Int_t trackLabel = TMath::Abs(track->GetLabel());
      
      AliVParticle * mcpart = GetMC()->GetTrack(trackLabel);
      if( !mcpart ) continue;
      
      Int_t  partInConeCharge = TMath::Abs(mcpart->Charge());
      Int_t  partInConePDG    = mcpart->PdgCode();
      Bool_t physPrimary      = mcpart->IsPhysicalPrimary();
      
      if (  fFillTrackOriginHistograms && 
            partInConeCharge > 0 &&  TMath::Abs(partInConePDG) != 11 ) // exclude electrons and neutrals
      {
        Int_t mcChTag = 3;
        if      ( TMath::Abs(partInConePDG) == 211  )  mcChTag = 0;
        else if ( TMath::Abs(partInConePDG) == 321  )  mcChTag = 1; 
        else if ( TMath::Abs(partInConePDG) == 2212 )  mcChTag = 2; 
        
        if ( physPrimary )
          fhPtTrackInConeMCPrimary  [mcChTag]->Fill(ptTrig , pTtrack, GetEventWeight()*weightTrig);
        else
          fhPtTrackInConeMCSecondary[mcChTag]->Fill(ptTrig , pTtrack, GetEventWeight()*weightTrig);
      }
    }
    
    if ( fStudyExoticTrigger && fIsExoticTrigger )
    {
      fhPtInConeExoTrigger      ->Fill(ptTrig , pTtrack, GetEventWeight()*weightTrig);
      fhPtTrackInConeExoTrigger ->Fill(ptTrig , pTtrack, GetEventWeight()*weightTrig);
    }
    
    if(fStudyPtCutInCone)
    {
      for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
      {
        if ( pTtrack > fMinPtCutInCone[icut] ) 
        {
          fConeptsumTrackPerMinCut[icut]+=pTtrack;
          fConeNTrackPerMinCut    [icut]++;
        }
        
        if ( pTtrack < fMaxPtCutInCone[icut] ) coneptsumTrackPerMaxCut[icut]+=pTtrack;
      }
    }
    
    if(fStudyEtaCutInCone)
    {
      for(Int_t icut = 0; icut < fNEtaCutsInCone; icut++) 
      {
        if ( TMath::Abs(track->Eta()) < fEtaCutInCone[icut] ) coneptsumTrackPerEtaCut[icut]+=pTtrack;
      }
    }
    
    if(fStudyRCutInCone)
    {
      Float_t distance = GetIsolationCut()->Radius(aodParticle->Eta(), aodParticle->Phi(), track->Eta(), track->Phi());
      for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
      {
        if ( distance < fRCutInCone[icut] ) 
        {
          coneptsumTrackPerRCut[icut]+=pTtrack;
          fhPtTrackInConePerRCut->Fill(icut+1, pTtrack, GetEventWeight()*weightTrig);
        }
      }
    }
    
    if(fStudyNCellsCut)
    {
      for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
      {
        if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
          fhPtTrackInConePerNCellPerSM[ishsh]->Fill(pTtrack, fTrigSupMod, fNCellsWithWeight);
        
        if ( fNCellsWithWeight >= fNCellsInCandidate[icut] ) 
        {
          coneptsumTrackPerNCellCut[icut]+=pTtrack;
          fhPtTrackInConePerNCellCut->Fill(icut+1, pTtrack, GetEventWeight()*weightTrig);
        }
      }
    }
    
    if(fStudyExoticTrigger)
    { 
      for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
      {
        if ( fClusterExoticity < fExoCutInCandidate[icut] ) 
        {
          coneptsumTrackPerExoCut[icut]+=pTtrack;
          fhPtTrackInConePerExoCut->Fill(icut+1, pTtrack, GetEventWeight()*weightTrig);
        }
      }
    }
    
    Bool_t okTOF = kFALSE ;
    Int_t trackBC = 0;
    if(fStudyTracksInCone)
    {
      phitrack = track->Phi();
      etatrack = track->Eta();
      
      fhPhiTrackInCone->Fill(pTtrack, phitrack, GetEventWeight()*weightTrig);
      fhEtaTrackInCone->Fill(pTtrack, etatrack, GetEventWeight()*weightTrig);
      fhEtaPhiTrackInCone->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
      
      // TOF
      ULong_t status = track->GetStatus();
      okTOF   = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
      trackBC = track->GetTOFBunchCrossing(bz);
      Double32_t tof = track->GetTOFsignal()*1e-3;    
      
      Int_t vtxBC = GetReader()->GetVertexBC();
      if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA) fhPtTrackInConeVtxBC0->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      
      if(okTOF)
      {
        fhTrackTOFInCone->Fill(pTtrack,tof,GetEventWeight()*weightTrig);
        
        if(fStudyExoticTrigger && fIsExoticTrigger)
          fhTrackTOFInConeExoTrigger->Fill(pTtrack,tof,GetEventWeight()*weightTrig);
        
        if(trackBC == 0) 
        {
          fhPtTrackInConeTOFBC0 ->Fill(ptTrig , pTtrack , GetEventWeight()*weightTrig);
          fhPhiTrackInConeTOFBC0->Fill(pTtrack, phitrack, GetEventWeight()*weightTrig);
          fhEtaTrackInConeTOFBC0->Fill(pTtrack, etatrack, GetEventWeight()*weightTrig);
          fhEtaPhiTrackInConeTOFBC0->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
          fhTrackTOFInConeBC0   ->Fill(pTtrack, tof     , GetEventWeight()*weightTrig);
          coneptsumTrackTOFBC0 += pTtrack;
        }
        else 
        {
          fhPtTrackInConeTOFBCN ->Fill(ptTrig , pTtrack , GetEventWeight()*weightTrig);
          fhPhiTrackInConeTOFBCN->Fill(pTtrack, phitrack, GetEventWeight()*weightTrig);
          fhEtaTrackInConeTOFBCN->Fill(pTtrack, etatrack, GetEventWeight()*weightTrig);
          fhEtaPhiTrackInConeTOFBCN->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
          coneptsumTrackTOFBCN += pTtrack;
        }
      }
      else
      {
        fhPtTrackInConeTOFNo ->Fill(ptTrig, pTtrack , GetEventWeight()*weightTrig);
        fhPhiTrackInConeTOFNo->Fill(ptTrig, phitrack, GetEventWeight()*weightTrig);
        fhEtaTrackInConeTOFNo->Fill(ptTrig, etatrack, GetEventWeight()*weightTrig);
        fhEtaPhiTrackInConeTOFNo->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
        coneptsumTrackTOFNo  += pTtrack;
      }
      
      
      Bool_t bITSRefit    = (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
      Bool_t bConstrained = (!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1));
      //printf("Track %d, pt %2.2f, eta %2.2f, phi %2.2f, SPDRefit %d, refit %d\n",
      //       itrack, pTtrack, etatrack, phitrack, bConstrained, bITSRefit);
      
      if(bConstrained)
      {
        if(bITSRefit)
        {
          coneptsumTrackITSRefitOnSPDOff  += pTtrack;
          fhPtTrackInConeITSRefitOnSPDOff ->Fill(ptTrig, pTtrack , GetEventWeight()*weightTrig);
          fhPhiTrackInConeITSRefitOnSPDOff->Fill(ptTrig, phitrack, GetEventWeight()*weightTrig);
          fhEtaTrackInConeITSRefitOnSPDOff->Fill(ptTrig, etatrack, GetEventWeight()*weightTrig);
          fhEtaPhiTrackInConeITSRefitOnSPDOff->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
        }
        else
        {
          coneptsumTrackITSRefitOffSPDOff += pTtrack;
          fhPtTrackInConeITSRefitOffSPDOff ->Fill(ptTrig, pTtrack , GetEventWeight()*weightTrig);
          fhPhiTrackInConeITSRefitOffSPDOff->Fill(ptTrig, phitrack, GetEventWeight()*weightTrig);
          fhEtaTrackInConeITSRefitOffSPDOff->Fill(ptTrig, etatrack, GetEventWeight()*weightTrig);
          fhEtaPhiTrackInConeITSRefitOffSPDOff->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
        }
      }
      else
      {
        coneptsumTrackITSRefitOnSPDOn   += pTtrack;
        fhPtTrackInConeITSRefitOnSPDOn ->Fill(ptTrig, pTtrack , GetEventWeight()*weightTrig);
        fhPhiTrackInConeITSRefitOnSPDOn->Fill(ptTrig, phitrack, GetEventWeight()*weightTrig);
        fhEtaTrackInConeITSRefitOnSPDOn->Fill(ptTrig, etatrack, GetEventWeight()*weightTrig);
        fhEtaPhiTrackInConeITSRefitOnSPDOn->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
      }
      
      if(okTOF && trackBC == 0 && !bConstrained)
      {
        fhPtTrackInConeTOFBC0ITSRefitOnSPDOn ->Fill(ptTrig , pTtrack , GetEventWeight()*weightTrig);
        fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(pTtrack, phitrack, GetEventWeight()*weightTrig);
        fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(pTtrack, etatrack, GetEventWeight()*weightTrig);
        fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(etatrack, phitrack, GetEventWeight()*weightTrig);
        coneptsumTrackTOFBC0ITSRefitOnSPDOn += pTtrack;
      }
      
      // DCA
      //
      if(ptTrig > 10)
      {
        Double_t dca[2]   = {1e6,1e6};
        Double_t covar[3] = {1e6,1e6,1e6};
        
        Double_t dcaCons  = -999;
        if ( GetReader()->GetDataType() == AliCaloTrackReader::kAOD )
        {
          AliAODTrack * aodTrack = dynamic_cast<AliAODTrack*>(track);
          dcaCons = aodTrack->DCA();
        }
        
        track->PropagateToDCA(GetReader()->GetInputEvent()->GetPrimaryVertex(),bz,100.,dca,covar);
        
        if(dcaCons == -999)
        {
          fhPtTrackInConeDCA[0]->Fill(pTtrack,  dca[0], GetEventWeight()*weightTrig);
          fhPtTrackInConeDCA[1]->Fill(pTtrack,  dca[1], GetEventWeight()*weightTrig);
        }
        else
        {
          fhPtTrackInConeDCA[2]->Fill(pTtrack, dcaCons, GetEventWeight()*weightTrig);
        }
      } // trigger pt cut
    }
    
    if(IsPileUpAnalysisOn())
    {
      if(GetReader()->IsPileUpFromSPD())             
      {  
        fhPtInConePileUp[0]            ->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
        if(fStudyTracksInCone)
        {
          if(okTOF && trackBC!=0 )                         fhPtTrackInConeOtherBCPileUpSPD->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
          if(okTOF && trackBC==0 )                         fhPtTrackInConeBC0PileUpSPD    ->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig); 
        }
      }
      
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig, pTtrack, GetEventWeight()*weightTrig);
    }
  }
  
  if ( fFillPerSMHistograms )     
    fhConeSumPtTrackPerSM[aodParticle->GetSModNumber()]->Fill(ptTrig, coneptsumTrack, GetEventWeight()*weightTrig);
  
  if ( fFillPerTCardIndexHistograms )     
    fhConeSumPtTrackPerTCardIndex[fTCardIndex]->Fill(ptTrig, coneptsumTrack, GetEventWeight()*weightTrig);
  
  if ( fStudyExoticTrigger && fIsExoticTrigger )
    fhConeSumPtTrackExoTrigger  ->Fill(ptTrig, coneptsumTrack  , GetEventWeight()*weightTrig);
  
  if ( fStudyTracksInCone )
  {
    fhConeSumPtTrackTOFBC0->Fill(ptTrig, coneptsumTrackTOFBC0, GetEventWeight()*weightTrig);
    fhConeSumPtTrackTOFBCN->Fill(ptTrig, coneptsumTrackTOFBCN, GetEventWeight()*weightTrig);
    fhConeSumPtTrackTOFNo ->Fill(ptTrig, coneptsumTrackTOFNo , GetEventWeight()*weightTrig);
    
    fhConeSumPtTrackITSRefitOnSPDOn  ->Fill(ptTrig, coneptsumTrackITSRefitOnSPDOn  , GetEventWeight()*weightTrig);
    fhConeSumPtTrackITSRefitOffSPDOff->Fill(ptTrig, coneptsumTrackITSRefitOffSPDOff, GetEventWeight()*weightTrig);
    fhConeSumPtTrackITSRefitOnSPDOff ->Fill(ptTrig, coneptsumTrackITSRefitOnSPDOff , GetEventWeight()*weightTrig);
    
    fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, coneptsumTrackTOFBC0ITSRefitOnSPDOn, GetEventWeight()*weightTrig);
  }
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhConeSumPtTrackPerMinPtCut->Fill(icut+1, fConeptsumTrackPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeSumPtTrackPerMaxPtCut->Fill(icut+1, coneptsumTrackPerMaxCut [icut], GetEventWeight()*weightTrig);
        fhConeNTrackPerMinPtCut    ->Fill(icut+1, fConeNTrackPerMinCut    [icut], GetEventWeight()*weightTrig);
      }
      else
      {
        fhConeNTrackPerMinPtCutCent    ->Fill(icut, fConeNTrackPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeSumPtTrackPerMinPtCutCent->Fill(icut, fConeptsumTrackPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
      }
    }
  }
  
  if ( fStudyEtaCutInCone )
  {
    for(Int_t icut = 0; icut < fNEtaCutsInCone; icut++) 
    {
      fhConeSumPtTrackPerEtaCut ->Fill(icut+1, coneptsumTrackPerEtaCut[icut], GetEventWeight()*weightTrig);
    }
  }
  
  if ( fStudyRCutInCone )
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      fhConeSumPtTrackPerRCut ->Fill(icut+1, coneptsumTrackPerRCut[icut], GetEventWeight()*weightTrig);
    }
  }
  
  if ( fStudyNCellsCut )
  {
    if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
      fhConeSumPtTrackPerNCellPerSM[ishsh]->Fill(coneptsumTrack, fTrigSupMod, fNCellsWithWeight);
    
    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      fhConeSumPtTrackPerNCellCut->Fill(icut, coneptsumTrackPerNCellCut[icut], GetEventWeight()*weightTrig);
    }
  }
  
  if ( fStudyExoticTrigger )
  { 
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      fhConeSumPtTrackPerExoCut->Fill(icut, coneptsumTrackPerExoCut[icut], GetEventWeight()*weightTrig);
    }
  }
}

//________________________________________________________________________________________________
/// Get the track pT or sum of pT in a cone at 45 degrees from trigger or in eta-phi bands.
/// Fill additional histograms not done in AliIsolationCut
//________________________________________________________________________________________________
void AliAnaParticleIsolation::StudyTracksUEInCone(AliCaloTrackParticleCorrelation * pCandidate)
{
  if( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kOnlyNeutral ) return ;
  
  Float_t conesize    = GetIsolationCut()->GetConeSize();
  Float_t conesizegap = GetIsolationCut()->GetConeSizeBandGap();
  Float_t tpcEtaSize  = GetIsolationCut()->GetTPCEtaSize();
  Float_t tpcPhiSize  = GetIsolationCut()->GetTPCPhiSize();
  Float_t distToTrig  = GetIsolationCut()->GetMinDistToTrigger();

  // Bands normalization
  //
  Float_t coneArea = conesize*conesize*TMath::Pi(); // Area = pi R^2, isolation cone area
  Float_t coneAreaGap = (conesize+conesizegap)*(conesize+conesizegap)*TMath::Pi(); // Area = pi R^2, isolation cone area

  if ( distToTrig > 0 )
  {
    coneArea    -= distToTrig*distToTrig*TMath::Pi();
    coneAreaGap -= distToTrig*distToTrig*TMath::Pi();
  }
  
  // Area of band, rectangle minus isolation region
  Float_t etaBandArea = 2*(conesize+conesizegap) * tpcEtaSize - coneAreaGap;
  Float_t phiBandArea = 2*(conesize+conesizegap) * tpcPhiSize - coneAreaGap;
  //printf("Area band eta %f phi %f\n",etaBandArea,phiBandArea);

  //Double_t sumptPerp = 0. ;
  Double_t sumptPerpBC0 = 0. ;
  Double_t sumptPerpITSSPD = 0. ;
  Double_t sumptPerpBC0ITSSPD = 0.;
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fConeptsumPerpTrackPerMinCut   [icut] = 0;
      fConeNPerpTrackPerMinCut       [icut] = 0;
      fConeptsumTrackSubPerpPerMinCut[icut] = 0;
      fConeNTrackSubPerpPerMinCut    [icut] = 0;
      
      fConeptsumEtaBandTrackPerMinCut   [icut] = 0;
      fConeNEtaBandTrackPerMinCut       [icut] = 0;
      fConeptsumTrackSubEtaBandPerMinCut[icut] = 0;
      fConeNTrackSubEtaBandPerMinCut    [icut] = 0;
      
      fConeptsumPhiBandTrackPerMinCut   [icut] = 0;
      fConeNPhiBandTrackPerMinCut       [icut] = 0;
      fConeptsumTrackSubPhiBandPerMinCut[icut] = 0;
      fConeNTrackSubPhiBandPerMinCut    [icut] = 0;
    }
  }
  
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  Float_t weightTrig= pCandidate->GetWeight();

  Double_t bz = GetReader()->GetInputEvent()->GetMagneticField();
  
  TObjArray * trackList = GetCTSTracks() ;
  for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
  {
    AliVTrack* track = (AliVTrack *) trackList->At(itrack);
    
    if(!track)
    {
      AliWarning("Track not available?");
      continue;
    }
    
    // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
    // do not count the candidate or the daughters of the candidate
    // in the isolation conte
    if ( pCandidate->GetDetectorTag() == kCTS ) // make sure conversions are tagged as kCTS!!!
    {
      Int_t  trackID   = GetReader()->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
      Bool_t contained = kFALSE;
      
      for(Int_t i = 0; i < 4; i++) 
      {
        if( trackID == pCandidate->GetTrackLabel(i) ) contained = kTRUE;
      }
      
      if ( contained ) continue ;
    }
    
    Float_t ptTrack = track->Pt();
    Float_t etaTrack= track->Eta();
    Float_t phiTrack= track->Phi();
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Eta or phi band with respect the trigger cone
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, etaTrack, phiTrack);

    // Angles between trigger and track
    Double_t dEta = etaTrig - etaTrack;
    Double_t dPhi = phiTrig - phiTrack;

    // Shift phi angle when trigger is close to 0 or 360
    if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
    if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();
    
    if ( rad > conesize+conesizegap )
    {
      // Phi band
      //
      Bool_t takeIt = kTRUE;

      // Look only half TPC with respect candidate, avoid opposite side jet 
      if ( TMath::Abs(dPhi) > TMath::PiOver2() )  takeIt = kFALSE;
      
      // Within eta cone size
      if ( TMath::Abs(dEta) < conesize+conesizegap &&  takeIt ) 
      {
        //phiBandPtSumTrack += ptTrack;
        
        if ( fStudyPtCutInCone )
        {
          for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
          {
            if ( ptTrack > fMinPtCutInCone[icut] ) 
            {
              fConeptsumPhiBandTrackPerMinCut[icut]+=ptTrack;
              fConeNPhiBandTrackPerMinCut    [icut]++;
            }          
          }
        } // min pt
      } // phi band
      
      // Eta band
      //
      //takeIt = kTRUE;
      
      // Within phi cone size
      if ( TMath::Abs(dPhi) < conesize+conesizegap ) //&& takeIt )
      {
        //etaBandPtSumTrack += ptTrack;
        if ( fStudyPtCutInCone )
        {
          for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
          {
            if ( ptTrack > fMinPtCutInCone[icut] ) 
            {
              fConeptsumEtaBandTrackPerMinCut[icut]+=ptTrack;
              fConeNEtaBandTrackPerMinCut    [icut]++;
            }          
          }
        } // min pt
       
      } // eta band
      
    } // out of cone
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fill the histograms at +-45 degrees in phi from trigger particle, 
    // perpedicular to trigger axis in phi
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    Double_t dPhiPlu = dPhi + TMath::PiOver2();
    Double_t dPhiMin = dPhi - TMath::PiOver2();
    
    Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
    Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;
    
    Bool_t fillPerp = kFALSE;
    if ( TMath::Sqrt(argPlu) < conesize ) fillPerp = kTRUE ;
    if ( TMath::Sqrt(argMin) < conesize ) fillPerp = kTRUE ;
    
    if ( fillPerp ) 
    {
      //sumptPerp+=track->Pt();
      
      if ( fStudyPtCutInCone )
      {
        for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
        {
          if ( ptTrack > fMinPtCutInCone[icut] ) 
          {
            fConeptsumPerpTrackPerMinCut[icut]+=ptTrack;
            fConeNPerpTrackPerMinCut    [icut]++;
          }          
        }
      }
      
      if ( fStudyTracksInCone )
      {
        ULong_t status = track->GetStatus();
        Bool_t okTOF = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
        Int_t trackBC = track->GetTOFBunchCrossing(bz);
        //Double32_t tof = track->GetTOFsignal()*1e-3;    
        
        if ( okTOF && trackBC == 0 )
        {
          fhPtInPerpConeTOFBC0    ->Fill(  ptTrig, ptTrack, GetEventWeight()*weightTrig);
          fhEtaPhiInPerpConeTOFBC0->Fill(etaTrack,phiTrack, GetEventWeight()*weightTrig);
          
          sumptPerpBC0+=ptTrack;
        }
        
        Bool_t bConstrained = (!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1));
        //Bool_t bITSRefit    = (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
        if(!bConstrained) 
        {
          fhPtInPerpConeITSRefitOnSPDOn    ->Fill(ptTrig  , ptTrack, GetEventWeight()*weightTrig);
          fhEtaPhiInPerpConeITSRefitOnSPDOn->Fill(etaTrack,phiTrack, GetEventWeight()*weightTrig);
          
          sumptPerpITSSPD+=ptTrack;
        }
        
        if(okTOF && trackBC == 0 && !bConstrained)
        {
          fhPtInPerpConeTOFBC0ITSRefitOnSPDOn    ->Fill( ptTrig,  ptTrack, GetEventWeight()*weightTrig);
          fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn->Fill(etaTrack,phiTrack, GetEventWeight()*weightTrig);
          
          sumptPerpBC0ITSSPD+=ptTrack;
        }
        
        if ( ptTrig > 10 )
        {
          Double_t dca[2]   = {1e6,1e6};
          Double_t covar[3] = {1e6,1e6,1e6};
          
          Double_t dcaCons  = -999;
          if ( GetReader()->GetDataType() == AliCaloTrackReader::kAOD )
          {
            AliAODTrack * aodTrack = dynamic_cast<AliAODTrack*>(track);
            dcaCons = aodTrack->DCA();
          }
          
          track->PropagateToDCA(GetReader()->GetInputEvent()->GetPrimaryVertex(),bz,100.,dca,covar);
          
          if(dcaCons == -999)
          {
            fhPtTrackInPerpConeDCA[0]->Fill(ptTrack,  dca[0], GetEventWeight()*weightTrig);
            fhPtTrackInPerpConeDCA[1]->Fill(ptTrack,  dca[1], GetEventWeight()*weightTrig);
          }
          else
          {
            fhPtTrackInPerpConeDCA[2]->Fill(ptTrack, dcaCons, GetEventWeight()*weightTrig);
          }
        } // trigger pt cut for DCA
        
      } // study tracks in cone
    } // r in cone
  } // track loop
  
  if ( fStudyPtCutInCone )
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      // Normalize to cone area
      // For simplicity do not correct for excess of detector acceptance
      //
      
      // Perp cone normalization
      fConeptsumPerpTrackPerMinCut[icut]/=2.;
      fConeNPerpTrackPerMinCut    [icut]/=2.;
      
//      printf("cut %d before N eta %f phi %f, sum eta %f phi %f\n",icut,
//             fConeptsumEtaBandTrackPerMinCut[icut], fConeptsumPhiBandTrackPerMinCut[icut],
//             fConeNEtaBandTrackPerMinCut    [icut], fConeNPhiBandTrackPerMinCut    [icut]);
      
      if ( etaBandArea > 0 )
      {
        fConeptsumEtaBandTrackPerMinCut[icut] *= ( coneArea / etaBandArea );
        fConeNEtaBandTrackPerMinCut    [icut] *= ( coneArea / etaBandArea );
      }
      
      if ( phiBandArea > 0 )
      {
        fConeptsumPhiBandTrackPerMinCut[icut] *= ( coneArea / phiBandArea );
        fConeNPhiBandTrackPerMinCut    [icut] *= ( coneArea / phiBandArea );
      }
      
//      printf("cut %d after N eta %f phi %f, sum eta %f phi %f\n",icut,
//                fConeptsumEtaBandTrackPerMinCut[icut], fConeptsumPhiBandTrackPerMinCut[icut],
//                fConeNEtaBandTrackPerMinCut    [icut], fConeNPhiBandTrackPerMinCut    [icut]);
      
      // Subtract to content in trigger cone (it must be filled before on another method StudyTracksInCone())
      //
      fConeptsumTrackSubPerpPerMinCut   [icut] = fConeptsumTrackPerMinCut[icut] - fConeptsumPerpTrackPerMinCut   [icut];
      fConeNTrackSubPerpPerMinCut       [icut] = fConeNTrackPerMinCut    [icut] - fConeNPerpTrackPerMinCut       [icut];
      
      fConeptsumTrackSubEtaBandPerMinCut[icut] = fConeptsumTrackPerMinCut[icut] - fConeptsumEtaBandTrackPerMinCut[icut];
      fConeNTrackSubEtaBandPerMinCut    [icut] = fConeNTrackPerMinCut    [icut] - fConeNEtaBandTrackPerMinCut    [icut];
      
      fConeptsumTrackSubPhiBandPerMinCut[icut] = fConeptsumTrackPerMinCut[icut] - fConeptsumPhiBandTrackPerMinCut[icut];
      fConeNTrackSubPhiBandPerMinCut    [icut] = fConeNTrackPerMinCut    [icut] - fConeNPhiBandTrackPerMinCut    [icut];
      
      // Fill histograms
      //
      if ( !IsHighMultiplicityAnalysisOn() )
      {
        fhPerpConeSumPtTrackPerMinPtCut->Fill(icut+1, fConeptsumPerpTrackPerMinCut[icut], GetEventWeight()*weightTrig);
        fhPerpConeNTrackPerMinPtCut    ->Fill(icut+1, fConeNPerpTrackPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubPerpPerMinPtCut->Fill(icut+1, fConeptsumTrackSubPerpPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeNTrackSubPerpPerMinPtCut    ->Fill(icut+1, fConeNTrackSubPerpPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhEtaBandConeSumPtTrackPerMinPtCut->Fill(icut+1, fConeptsumEtaBandTrackPerMinCut[icut], GetEventWeight()*weightTrig);
        fhEtaBandConeNTrackPerMinPtCut    ->Fill(icut+1, fConeNEtaBandTrackPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubEtaBandPerMinPtCut->Fill(icut+1, fConeptsumTrackSubEtaBandPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeNTrackSubEtaBandPerMinPtCut    ->Fill(icut+1, fConeNTrackSubEtaBandPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhPhiBandConeSumPtTrackPerMinPtCut->Fill(icut+1, fConeptsumPhiBandTrackPerMinCut[icut], GetEventWeight()*weightTrig);
        fhPhiBandConeNTrackPerMinPtCut    ->Fill(icut+1, fConeNPhiBandTrackPerMinCut    [icut], GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubPhiBandPerMinPtCut->Fill(icut+1, fConeptsumTrackSubPhiBandPerMinCut[icut], GetEventWeight()*weightTrig);
        fhConeNTrackSubPhiBandPerMinPtCut    ->Fill(icut+1, fConeNTrackSubPhiBandPerMinCut    [icut], GetEventWeight()*weightTrig);
      }
      else 
      {
        fhPerpConeNTrackPerMinPtCutCent    ->Fill(icut, fConeNPerpTrackPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhPerpConeSumPtTrackPerMinPtCutCent->Fill(icut, fConeptsumPerpTrackPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubPerpPerMinPtCutCent->Fill(icut+1, fConeptsumTrackSubPerpPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeNTrackSubPerpPerMinPtCutCent    ->Fill(icut+1, fConeNTrackSubPerpPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhEtaBandConeNTrackPerMinPtCutCent    ->Fill(icut, fConeNEtaBandTrackPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhEtaBandConeSumPtTrackPerMinPtCutCent->Fill(icut, fConeptsumEtaBandTrackPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubEtaBandPerMinPtCutCent->Fill(icut+1, fConeptsumTrackSubEtaBandPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeNTrackSubEtaBandPerMinPtCutCent    ->Fill(icut+1, fConeNTrackSubEtaBandPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhPhiBandConeNTrackPerMinPtCutCent    ->Fill(icut, fConeNPhiBandTrackPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhPhiBandConeSumPtTrackPerMinPtCutCent->Fill(icut, fConeptsumPhiBandTrackPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        
        fhConeSumPtTrackSubPhiBandPerMinPtCutCent->Fill(icut+1, fConeptsumTrackSubPhiBandPerMinCut[icut], GetEventCentrality(), GetEventWeight()*weightTrig);
        fhConeNTrackSubPhiBandPerMinPtCutCent    ->Fill(icut+1, fConeNTrackSubPhiBandPerMinCut    [icut], GetEventCentrality(), GetEventWeight()*weightTrig);
      }
    }
  }
  
  if ( fStudyTracksInCone ) 
  {
    fhPerpConeSumPtTOFBC0         ->Fill(ptTrig, sumptPerpBC0   , GetEventWeight()*weightTrig);
    fhPerpConeSumPtITSRefitOnSPDOn->Fill(ptTrig, sumptPerpITSSPD, GetEventWeight()*weightTrig);
    fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, sumptPerpBC0ITSSPD, GetEventWeight()*weightTrig);
  }
  
}
