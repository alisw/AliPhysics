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

// --- Analysis system ---
#include "AliAnaParticleIsolation.h"
#include "AliCaloTrackReader.h"
#include "AliMCEvent.h"
#include "AliIsolationCut.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliVParticle.h"
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliMCAnalysisUtils.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

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
fReMakeIC(0),                     fMakeSeveralIC(0),
fFillTMHisto(0),                  fFillSSHisto(1),                          
fFillEMCALRegionHistograms(0),    fFillUEBandSubtractHistograms(1), 
fFillCellHistograms(0),
fFillOverlapHistograms(0),                        
fStudyTracksInCone(0),            fStudyMCConversionRadius(0),
fFillTaggedDecayHistograms(0),    fNDecayBits(0),
fDecayBits(),                     fDecayTagsM02Cut(0),
fFillNLMHistograms(0),
fLeadingOnly(0),                  fCheckLeadingWithNeutralClusters(0),
fSelectPrimariesInCone(0),        fMakePrimaryPi0DecayStudy(0),
fFillBackgroundBinHistograms(0),  fNBkgBin(0),
fFillPtTrigBinHistograms(0),      fNPtTrigBin(0),
fMinCellsAngleOverlap(0),
// Several IC
fNCones(0),                       fNPtThresFrac(0),
fConeSizes(),                     fPtThresholds(),
fPtFractions(),                   fSumPtThresholds(),
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
fhEIso(0),                        fhPtIso(0),
fhPtCentralityIso(0),             fhPtEventPlaneIso(0),
fhPtNLocMaxIso(0),
fhPhiIso(0),                      fhEtaIso(0),                              fhEtaPhiIso(0),
fhEtaPhiNoIso(0),
fhENoIso(0),                      fhPtNoIso(0),                             fhPtNLocMaxNoIso(0),
fhEIsoExoTrigger(0),              fhENoIsoExoTrigger(0),
fhPtIsoExoTrigger(0),             fhPtNoIsoExoTrigger(0),

fhPtInCone(0),
fhPtClusterInCone(0),             fhPtCellInCone(0),                        fhPtTrackInCone(0),
fhPtInConeExoTrigger(0),          fhPtClusterInConeExoTrigger(0),           fhPtTrackInConeExoTrigger(0),
fhPtTrackInConeOtherBCPileUpSPD(0), fhPtTrackInConeVtxBC0(0),
fhPtTrackInConeBC0PileUpSPD(0),
fhPtInConePileUp(),               fhPtInConeCent(0),

fhPerpConeSumPt(0),               fhPerpConeSumPtTOFBC0(0),               
fhPtInPerpCone(0),                fhPtInPerpConeTOFBC0(0),
fhEtaPhiInConeCluster(0),         fhEtaPhiCluster(0),
fhEtaPhiInConeTrack(0),           fhEtaPhiTrack(0),
fhEtaPhiInPerpCone(0),            fhEtaPhiInPerpConeTOFBC0(0),
fhEtaBandClusterEtaPhi(0),        fhPhiBandClusterEtaPhi(0),
fhEtaBandTrackEtaPhi(0),          fhPhiBandTrackEtaPhi(0),
fhEtaBandClusterPt(0),            fhPhiBandClusterPt(0),
fhEtaBandTrackPt(0),              fhPhiBandTrackPt(0),
fhEtaBandCell(0),                 fhPhiBandCell(0),
fhConePtLead(0),                  fhConePtLeadCluster(0),                   fhConePtLeadTrack(0),
fhConePtLeadClustervsTrack(0),    fhConePtLeadClusterTrackFrac(0),
fhConeSumPt(0),                   //fhPtLambda0Eiso(0),                       
fhConeSumPtCellTrack(0),
fhConeSumPtCell(0),               fhConeSumPtCluster(0),                    fhConeSumPtTrack(0),
fhConeSumPtExoTrigger(0),         fhConeSumPtClusterExoTrigger(0),          fhConeSumPtTrackExoTrigger(0),                      

fhConeSumPtEtaBandUECluster(0),             fhConeSumPtPhiBandUECluster(0),
fhConeSumPtEtaBandUETrack(0),               fhConeSumPtPhiBandUETrack(0),
fhConeSumPtEtaBandUECell(0),                fhConeSumPtPhiBandUECell(0),
fhConeSumPtTrigEtaPhi(0),
fhConeSumPtCellTrackTrigEtaPhi(0),
fhConeSumPtEtaBandUEClusterTrigEtaPhi(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhi(0),
fhConeSumPtEtaBandUETrackTrigEtaPhi(0),     fhConeSumPtPhiBandUETrackTrigEtaPhi(0),
fhConeSumPtEtaBandUECellTrigEtaPhi(0),      fhConeSumPtPhiBandUECellTrigEtaPhi(0),
fhConeSumPtEtaUESub(0),                     fhConeSumPtPhiUESub(0),
fhConeSumPtEtaUESubTrigEtaPhi(0),           fhConeSumPtPhiUESubTrigEtaPhi(0),
fhConeSumPtEtaUESubTrackCell(0),            fhConeSumPtPhiUESubTrackCell(0),
fhConeSumPtEtaUESubTrackCellTrigEtaPhi(0),  fhConeSumPtPhiUESubTrackCellTrigEtaPhi(0),
fhConeSumPtEtaUENormCluster(0),             fhConeSumPtPhiUENormCluster(0),
fhConeSumPtEtaUESubCluster(0),              fhConeSumPtPhiUESubCluster(0),
fhConeSumPtEtaUESubClusterTrigEtaPhi(0),    fhConeSumPtPhiUESubClusterTrigEtaPhi(0),
fhConeSumPtEtaUESubCell(0),                 fhConeSumPtPhiUESubCell(0),
fhConeSumPtEtaUESubCellTrigEtaPhi(0),       fhConeSumPtPhiUESubCellTrigEtaPhi(0),
fhConeSumPtEtaUENormTrack(0),               fhConeSumPtPhiUENormTrack(0),
fhConeSumPtEtaUESubTrack(0),                fhConeSumPtPhiUESubTrack(0),
fhConeSumPtEtaUESubTrackTrigEtaPhi(0),      fhConeSumPtPhiUESubTrackTrigEtaPhi(0),
fhFractionTrackOutConeEta(0),               fhFractionTrackOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConeEta(0),             fhFractionClusterOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConePhi(0),             fhFractionClusterOutConePhiTrigEtaPhi(0),
fhFractionCellOutConeEta(0),                fhFractionCellOutConeEtaTrigEtaPhi(0),
fhFractionCellOutConePhi(0),                fhFractionCellOutConePhiTrigEtaPhi(0),
fhConeSumPtClustervsTrack(0),               fhConeSumPtClusterTrackFrac(0),
fhConeSumPtEtaUESubClustervsTrack(0),       fhConeSumPtPhiUESubClustervsTrack(0),
fhConeSumPtCellvsTrack(0),
fhConeSumPtEtaUESubCellvsTrack(0),          fhConeSumPtPhiUESubCellvsTrack(0),
fhEtaBandClustervsTrack(0),                 fhPhiBandClustervsTrack(0),
fhEtaBandNormClustervsTrack(0),             fhPhiBandNormClustervsTrack(0),
fhEtaBandCellvsTrack(0),                    fhPhiBandCellvsTrack(0),
fhEtaBandNormCellvsTrack(0),                fhPhiBandNormCellvsTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiTrack(0),    fhConeSumPtSubNormvsConeSumPtTotPhiTrack(0),
fhConeSumPtSubvsConeSumPtTotEtaTrack(0),    fhConeSumPtSubNormvsConeSumPtTotEtaTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiCluster(0),  fhConeSumPtSubNormvsConeSumPtTotPhiCluster(0),
fhConeSumPtSubvsConeSumPtTotEtaCluster(0),  fhConeSumPtSubNormvsConeSumPtTotEtaCluster(0),
fhConeSumPtSubvsConeSumPtTotPhiCell(0),     fhConeSumPtSubNormvsConeSumPtTotPhiCell(0),
fhConeSumPtSubvsConeSumPtTotEtaCell(0),     fhConeSumPtSubNormvsConeSumPtTotEtaCell(0),
fhConeSumPtVSUETracksEtaBand(0),            fhConeSumPtVSUETracksPhiBand(0),
fhConeSumPtVSUEClusterEtaBand(0),           fhConeSumPtVSUEClusterPhiBand(0),
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
fhPtLeadConeBin(0),                         fhSumPtConeBin(0),
fhSumPtConeAfterEtaBandUESubBin(0),  
fhPtLeadConeBinMC(0),                       
fhSumPtConeAfterEtaBandUESubBinMC(0),
fhSumPtConeBinMC(0),
fhPtLeadConeBinDecay(0),                    fhSumPtConeBinDecay(0),
fhPtLeadConeBinLambda0(0),                  fhSumPtConeBinLambda0(0),
fhSumPtConeAfterEtaBandUESubBinLambda0(0),  fhPtLeadConeBinLambda0MC(0), 
fhSumPtConeAfterEtaBandUESubBinLambda0MC(0),fhSumPtConeBinLambda0MC(0),
fhPtTrigBinPtLeadCone(0),                   fhPtTrigBinSumPtCone(0),

fhPtTrigBinSumPtTrackCone(0),               fhPtTrigBinSumPtClusterCone(0),
fhPtTrigBinPtLeadConeMC(0),                 fhPtTrigBinSumPtConeMC(0),
fhPtTrigBinSumPtTrackConeMC(0),             fhPtTrigBinSumPtClusterConeMC(0),
fhPtTrigBinPtLeadConeDecay(0),              fhPtTrigBinSumPtConeDecay(0),
fhPtTrigBinSumPtTrackConeDecay(0),          fhPtTrigBinSumPtClusterConeDecay(0),
fhPtTrigBinLambda0vsPtLeadCone(0),          fhPtTrigBinLambda0vsSumPtCone(0),
fhPtTrigBinLambda0vsSumPtTrackCone(0),      fhPtTrigBinLambda0vsSumPtClusterCone(0),
fhPtTrigBinLambda0vsPtLeadConeMC(0),        fhPtTrigBinLambda0vsSumPtConeMC(0),
fhPtTrigBinLambda0vsSumPtTrackConeMC(0),    fhPtTrigBinLambda0vsSumPtClusterConeMC(0),
fhPtTrigBinLambda0vsSumPtConeMCNoOverlap(0),fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap(0), 
fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap(0), fhPtTrigBinLambda0vsSumPtConeMC1Overlap(0), 
fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap(0),fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap(0), 


// Number of local maxima in cluster
fhNLocMax(),
fhELambda0LocMax1(),              fhELambda1LocMax1(),
fhELambda0LocMax2(),              fhELambda1LocMax2(),
fhELambda0LocMaxN(),              fhELambda1LocMaxN(),
// PileUp
fhEIsoPileUp(),                   fhPtIsoPileUp(),
fhENoIsoPileUp(),                 fhPtNoIsoPileUp(),
fhTimeENoCut(0),                  fhTimeESPD(0),                  fhTimeESPDMulti(0),
fhTimeNPileUpVertSPD(0),          fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0),  fhTimePileUpMainVertexZDiamond(0),

fhPtClusterInConePerRCut(0),         fhPtClusterInConePerRCutLargePtTrig(0),
fhPtTrackInConePerRCut(0),           fhPtTrackInConePerRCutLargePtTrig(0),
fhConeSumPtClusterPerRCut(0),        fhConeSumPtClusterPerRCutLargePtTrig(0),
fhConeSumPtTrackPerRCut(0),          fhConeSumPtTrackPerRCutLargePtTrig(0),
fhConeNClusterPerMinPtCut(0),        fhConeNClusterPerMinPtCutLargePtTrig(0),
fhConeNTrackPerMinPtCut(0),          fhConeNTrackPerMinPtCutLargePtTrig(0),
fhPerpConeNTrackPerMinPtCut(0),      fhPerpConeNTrackPerMinPtCutLargePtTrig(0),
fhConeSumPtClusterPerMinPtCut(0),    fhConeSumPtClusterPerMinPtCutLargePtTrig(0),
fhConeSumPtTrackPerMinPtCut(0),      fhConeSumPtTrackPerMinPtCutLargePtTrig(0),
fhPerpConeSumPtTrackPerMinPtCut(0),  fhPerpConeSumPtTrackPerMinPtCutLargePtTrig(0),
fhConeSumPtClusterPerMaxPtCut(0),    fhConeSumPtClusterPerMaxPtCutLargePtTrig(0),
fhConeSumPtTrackPerMaxPtCut(0),      fhConeSumPtTrackPerMaxPtCutLargePtTrig(0),
fhConeSumPtTrackPerEtaCut(0),        fhConeSumPtTrackPerEtaCutLargePtTrig(0),

fhPtClusterInConePerNCellCut(0),     fhPtClusterInConePerNCellCutLargePtTrig(0),  
fhPtTrackInConePerNCellCut(0),       fhPtTrackInConePerNCellCutLargePtTrig(0),    
fhConeSumPtClusterPerNCellCut(0),    fhConeSumPtClusterPerNCellCutLargePtTrig(0), 
fhConeSumPtTrackPerNCellCut(0),      fhConeSumPtTrackPerNCellCutLargePtTrig(0), 
fhPtClusterInConePerExoCut(0),       fhPtClusterInConePerExoCutLargePtTrig(0), 
fhPtTrackInConePerExoCut(0),         fhPtTrackInConePerExoCutLargePtTrig(0), 
fhConeSumPtClusterPerExoCut(0),      fhConeSumPtClusterPerExoCutLargePtTrig(0), 
fhConeSumPtTrackPerExoCut(0),        fhConeSumPtTrackPerExoCutLargePtTrig(0),      

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
  
  for(Int_t i = 0; i < 5 ; i++)
  {
    fConeSizes[i]      = 0 ;
    
    for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
      fhSumPtLeadingPtMC[imc][i] = 0 ;
    
    for(Int_t j = 0; j < 5 ; j++)
    {
      fhPtThresIsolated             [i][j] = 0 ;
      fhPtFracIsolated              [i][j] = 0 ;
      fhSumPtIsolated               [i][j] = 0 ;
      
      fhEtaPhiPtThresIso            [i][j] = 0 ;
      fhEtaPhiPtThresDecayIso       [i][j] = 0 ;
      fhPtPtThresDecayIso           [i][j] = 0 ;
      
      fhEtaPhiPtFracIso             [i][j] = 0 ;
      fhEtaPhiPtFracDecayIso        [i][j] = 0 ;
      fhPtPtFracDecayIso            [i][j] = 0 ;
      fhPtPtSumDecayIso             [i][j] = 0 ;
      fhPtSumDensityIso             [i][j] = 0 ;
      fhPtSumDensityDecayIso        [i][j] = 0 ;
      fhEtaPhiSumDensityIso         [i][j] = 0 ;
      fhEtaPhiSumDensityDecayIso    [i][j] = 0 ;
      fhPtFracPtSumIso              [i][j] = 0 ;
      fhPtFracPtSumDecayIso         [i][j] = 0 ;
      fhEtaPhiFracPtSumIso          [i][j] = 0 ;
      fhEtaPhiFracPtSumDecayIso     [i][j] = 0 ;
      
      for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
      {
        fhPtThresIsolatedMC[imc][i][j] = 0 ;
        fhPtFracIsolatedMC [imc][i][j] = 0 ;
        fhSumPtIsolatedMC  [imc][i][j] = 0 ;
        
      }
    }
  }
    
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
  
  for(Int_t i = 0; i < 5 ; i++)
  {
    fPtFractions    [i] = 0 ;
    fPtThresholds   [i] = 0 ;
    fSumPtThresholds[i] = 0 ;
    
    fhSumPtLeadingPt    [i] = 0 ;
    fhPtLeadingPt       [i] = 0 ;
    fhPerpSumPtLeadingPt[i] = 0 ;
    fhPerpPtLeadingPt   [i] = 0 ;
  }
  
  for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
  {
    fhPtNoIsoMC  [imc]    = 0;
    fhPtIsoMC    [imc]    = 0;
    fhPhiIsoMC   [imc]    = 0;
    fhEtaIsoMC   [imc]    = 0;
    
    for(Int_t i = 0; i < 2 ; i++)
    {
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
    fhTrackMatchedDEta[i] = 0 ;             fhTrackMatchedDPhi[i] = 0 ;   fhTrackMatchedDEtaDPhi  [i] = 0 ;
    fhdEdx            [i] = 0 ;             fhEOverP          [i] = 0 ;   fhTrackMatchedMCParticle[i] = 0 ;
    fhELambda0        [i] = 0 ;             fhPtLambda0       [i] = 0 ;   //fhELambda1        [i] = 0 ;
    fhELambda0TRD     [i] = 0 ;             fhPtLambda0TRD    [i] = 0 ;   //fhELambda1TRD     [i] = 0 ;
    
    // Number of local maxima in cluster
    fhNLocMax        [i] = 0 ;
    fhELambda0LocMax1[i] = 0 ;              fhELambda1LocMax1[i] = 0 ;
    fhELambda0LocMax2[i] = 0 ;              fhELambda1LocMax2[i] = 0 ;
    fhELambda0LocMaxN[i] = 0 ;              fhELambda1LocMaxN[i] = 0 ;
    
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
    fhPtPrimMCiso[i] = 0;
    fhEPrimMC    [i] = 0;
    fhPtPrimMC   [i] = 0;
    fhEtaPrimMC  [i] = 0;
    fhPhiPrimMC  [i] = 0;
  }
  
  // Pile-Up
  
  for(Int_t i = 0 ; i < 7 ; i++)
  {
    fhPtInConePileUp[i] = 0 ;
    fhEIsoPileUp    [i] = 0 ;
    fhPtIsoPileUp   [i] = 0 ;
    fhENoIsoPileUp  [i] = 0 ;
    fhPtNoIsoPileUp [i] = 0 ;
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
}

//_______________________________________________________________________________________________
/// Get the clusters pT or sum of pT in phi/eta bands or at 45 degrees from trigger.
//_______________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                  Float_t & etaBandPtSum, Float_t & phiBandPtSum)
{
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t conesize   = GetIsolationCut()->GetConeSize();
  
  // Select the Calorimeter
  TObjArray * pl = 0x0;
  if      (GetCalorimeter() == kPHOS )
    pl    = GetPHOSClusters();
  else if (GetCalorimeter() == kEMCAL)
    pl    = GetEMCALClusters();
  
  if(!pl) return ;
  
  // Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
  for(Int_t icluster=0; icluster < pl->GetEntriesFast(); icluster++)
  {
    AliVCluster* cluster = (AliVCluster *) pl->At(icluster);
    
    if ( !cluster )
    {
      AliWarning("Cluster not available?");
      continue;
    }
    
    // Do not count the candidate (photon or pi0) or the daughters of the candidate
    if(cluster->GetID() == pCandidate->GetCaloLabel(0) ||
       cluster->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;
    
    // Remove matched clusters to tracks if Neutral and Track info is used
    if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged &&
       IsTrackMatched(cluster,GetReader()->GetInputEvent())) continue ;
    
    cluster->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line
    
    // Exclude particles in cone
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, fMomentum.Eta(), fMomentum.Phi());
    
    // Histograms of eta and phi for all clusters
    if(fFillUEBandSubtractHistograms > 1)
      fhEtaPhiCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
  
    if(rad < conesize)
    {
    	// Histograms for all clusters in cone
      if(fFillUEBandSubtractHistograms > 1)
        fhEtaPhiInConeCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
     
      continue ;
    }
      
    // Fill histogram for UE in phi band in EMCal acceptance
    if(fMomentum.Eta() > (etaTrig-conesize) && fMomentum.Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=fMomentum.Pt();
      if(fFillUEBandSubtractHistograms > 1) 
        fhPhiBandClusterEtaPhi->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
      fhPhiBandClusterPt    ->Fill(ptTrig, fMomentum.Pt (), GetEventWeight());      
    }
    
    // Fill histogram for UE in eta band in EMCal acceptance
    if(fMomentum.Phi() > (phiTrig-conesize) && fMomentum.Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=fMomentum.Pt();
      if(fFillUEBandSubtractHistograms > 1)
        fhEtaBandClusterEtaPhi->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
      fhEtaBandClusterPt ->Fill(ptTrig, fMomentum.Pt(), GetEventWeight());
    }
  }
  
  fhConeSumPtEtaBandUECluster->Fill(ptTrig, etaBandPtSum, GetEventWeight());
  fhConeSumPtPhiBandUECluster->Fill(ptTrig, phiBandPtSum, GetEventWeight());
  
  if(fFillUEBandSubtractHistograms > 1)
  {
    fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSum *GetEventWeight()); // Check
    fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSum *GetEventWeight()); // Check
  }
}

//________________________________________________________________________________________________
/// Get the cells amplitude or sum of amplitude in phi/eta bands or at 45 degrees from trigger.
//________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloCellUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                      Float_t & etaBandPtSumCells, Float_t & phiBandPtSumCells)
{
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  Float_t conesize = GetIsolationCut()->GetConeSize();
  
  Float_t phiTrig = pCandidate->Phi();
  if(phiTrig<0) phiTrig += TMath::TwoPi();
    
  Float_t etaTrig = pCandidate->Eta();
  
  if(pCandidate->GetDetectorTag()==kEMCAL)
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    Int_t absId = -999;
    
    if (eGeom->GetAbsCellIdFromEtaPhi(etaTrig,phiTrig,absId))
    {
      if(!eGeom->CheckAbsCellId(absId)) return ;
      
      // Get absolute (col,row) of trigger particle
      Int_t nSupMod = eGeom->GetSuperModuleNumber(absId);
      Int_t nModule = -1;
      Int_t imEta=-1, imPhi=-1;
      Int_t ieta =-1, iphi =-1;
      
      if (eGeom->GetCellIndex(absId,nSupMod,nModule,imPhi,imEta))
      {
        eGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,imPhi,imEta,iphi,ieta);
        
        Int_t colTrig = ieta;
        if (nSupMod % 2) colTrig = AliEMCALGeoParams::fgkEMCALCols + ieta ;
          
        Int_t rowTrig = iphi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
        
        Int_t sqrSize = int(conesize/0.0143);
        
        AliVCaloCells * cells = GetEMCALCells();
        
        Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ; // 24*(16/3) 5 full-size Sectors (2 SM) + 1 one-third Sector (2 SM)
        Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols;
        //  printf("nTotalRows %i, nTotalCols %i\n",nTotalRows,nTotalCols);
        // Loop on cells in eta band
        
        Int_t irowmin = rowTrig-sqrSize;
        if ( irowmin < 0 ) irowmin = 0 ;
          
        Int_t irowmax = rowTrig+sqrSize;
        if ( irowmax > AliEMCALGeoParams::fgkEMCALRows ) irowmax = AliEMCALGeoParams::fgkEMCALRows;
        
        for(Int_t irow = irowmin; irow <irowmax; irow++)
        {
          for(Int_t icol = 0; icol < nTotalCols; icol++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if(inSector==5) continue;
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;
            
            // Exclude cells in cone
            if(TMath::Abs(icol-colTrig) < sqrSize || TMath::Abs(irow-rowTrig) < sqrSize){
              continue ;
            }
              
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId)) continue;
            
            etaBandPtSumCells += cells->GetCellAmplitude(iabsId);
            
            fhEtaBandCell->Fill(colTrig, rowTrig, GetEventWeight());
            
            //          printf("ETA inSupMod %i,irowLoc %i,icolLoc %i, iabsId %i, etaBandPtSumCells %f\n",nSupMod,irowLoc,icolLoc,iabsId,etaBandPtSumCells);
          }
        }
          
        Int_t icolmin = colTrig-sqrSize;
        if ( icolmin < 0 ) icolmin = 0;
        
        Int_t icolmax = colTrig+sqrSize;
        if ( icolmax > AliEMCALGeoParams::fgkEMCALCols ) icolmax = AliEMCALGeoParams::fgkEMCALCols;
	      
        // Loop on cells in phi band
        for(Int_t icol = icolmin; icol < icolmax; icol++)
        {
          for(Int_t irow = 0; irow < nTotalRows; irow++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if ( inSector == 5 ) continue ;
            
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            //    printf("icol %i, irow %i, inSector %i\n",icol,irow ,inSector);
              
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              //	printf("icol < AliEMCALGeoParams::fgkEMCALCols %i\n",AliEMCALGeoParams::fgkEMCALCols );
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              //      printf("icol > AliEMCALGeoParams::fgkEMCALCols -1 %i\n",AliEMCALGeoParams::fgkEMCALCols -1 );
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;   // Stesso problema di sopra //
            
            // Exclude cells in cone
            if(TMath::Abs(icol-colTrig) < sqrSize) {
              //printf("TMath::Abs(icol-colTrig) %i < sqrSize %i\n",TMath::Abs(icol-colTrig) ,sqrSize);continue ;
            }
            if(TMath::Abs(irow-rowTrig) < sqrSize) {
              //printf("TMath::Abs(irow-rowTrig) %i < sqrSize %i\n",TMath::Abs(irow-rowTrig) ,sqrSize);continue ;
            }
            
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId))
            {
              AliWarning(Form("!eGeom->CheckAbsCellId(iabsId=%i) inSupMod %i irowLoc %i icolLoc %i",iabsId,inSupMod, irowLoc, icolLoc));
              continue;
            }
              
            phiBandPtSumCells += cells->GetCellAmplitude(iabsId);
            
            fhPhiBandCell->Fill(colTrig, rowTrig, GetEventWeight());
            //printf("inSupMod %i,irowLoc %i,icolLoc %i, iabsId %i, phiBandPtSumCells %f\n",nSupMod,irowLoc,icolLoc,iabsId,phiBandPtSumCells);
          }
        }
      }
    }
  }
  
  Float_t ptTrig = pCandidate->Pt();
  
  fhConeSumPtEtaBandUECell          ->Fill(ptTrig ,          etaBandPtSumCells, GetEventWeight());
  fhConeSumPtPhiBandUECell          ->Fill(ptTrig ,          phiBandPtSumCells, GetEventWeight());
  fhConeSumPtEtaBandUECellTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSumCells *GetEventWeight()); // Check
  fhConeSumPtPhiBandUECellTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSumCells *GetEventWeight()); // Check
}

//________________________________________________________________________________________________
/// Get the track pT or sum of pT in phi/eta bands or at 45 degrees from trigger.
//________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateTrackUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                                   Float_t & etaBandPtSum, Float_t & phiBandPtSum)
{
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  Float_t conesize   = GetIsolationCut()->GetConeSize();
  
  Double_t sumptPerp = 0. ;
  Double_t sumptPerpBC0 = 0. ;
  Double_t sumptPerpITSSPD = 0. ;
  Double_t sumptPerpBC0ITSSPD = 0.;
  
  Float_t coneptsumPerpTrackPerMinCut[20];
  Float_t coneNPerpTrackPerMinCut    [20];
  
  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      coneptsumPerpTrackPerMinCut[icut] = 0;
      coneNPerpTrackPerMinCut    [icut] = 0;
    }
  }
  
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
  Double_t bz = GetReader()->GetInputEvent()->GetMagneticField();
  
  TObjArray * trackList   = GetCTSTracks() ;
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
    
    // Histogram of eta:phi for all tracks
    if(fFillUEBandSubtractHistograms > 1)
      fhEtaPhiTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
    
    //exclude particles in cone
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, track->Eta(), track->Phi());
    if(rad < conesize)
    {
      // Histogram of eta:phi for all tracks in cone
      if(fFillUEBandSubtractHistograms > 1)
        fhEtaPhiInConeTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
      continue ;
    }
    
    // Fill histogram for UE in phi band
    if(track->Eta() > (etaTrig-conesize) && track->Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=track->Pt();
      if(fFillUEBandSubtractHistograms > 1) 
        fhPhiBandTrackEtaPhi->Fill(track->Eta(), track->Phi(), GetEventWeight());
      fhPhiBandTrackPt->Fill(ptTrig      , track->Pt (), GetEventWeight());
    }
    
    // Fill histogram for UE in eta band in EMCal acceptance
    if(track->Phi() > (phiTrig-conesize) && track->Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=track->Pt();
      if(fFillUEBandSubtractHistograms > 1)
        fhEtaBandTrackEtaPhi->Fill(track->Eta(), track->Phi(), GetEventWeight());
      fhEtaBandTrackPt    ->Fill(ptTrig      , track->Pt (), GetEventWeight());
    }
     
    // Fill the histograms at +-45 degrees in phi from trigger particle, perpedicular to trigger axis in phi
    //++++++++
    Double_t dPhi   = phiTrig - track->Phi() + TMath::PiOver2();
    Double_t dEta   = etaTrig - track->Eta();
    Double_t arg    = dPhi*dPhi + dEta*dEta;
 
    if(TMath::Sqrt(arg) < conesize)
    {
      sumptPerp+=track->Pt();

      fhPtInPerpCone->Fill(ptTrig, track->Pt(), GetEventWeight());
      
      if(fStudyTracksInCone)
      {
        fhEtaPhiInPerpCone->Fill(track->Eta(),track->Phi(), GetEventWeight());

        ULong_t status = track->GetStatus();
        Bool_t okTOF = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
        Int_t trackBC = track->GetTOFBunchCrossing(bz);
        //Double32_t tof = track->GetTOFsignal()*1e-3;    
        
        if(okTOF && trackBC == 0)
        {
          fhPtInPerpConeTOFBC0->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeTOFBC0->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpBC0+=track->Pt();
        }
        
        Bool_t bConstrained = (!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1));
        //Bool_t bITSRefit    = (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
        if(!bConstrained) 
        {
          fhPtInPerpConeITSRefitOnSPDOn->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeITSRefitOnSPDOn->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpITSSPD+=track->Pt();
        }
        
        if(okTOF && trackBC == 0 && !bConstrained)
        {
          fhPtInPerpConeTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpBC0ITSSPD+=track->Pt();
        }
      }
    }
    
    //----------
    dPhi = phiTrig - track->Phi() - TMath::PiOver2();
    arg  = dPhi*dPhi + dEta*dEta;
                                  
    if(TMath::Sqrt(arg) < conesize)
    {
      sumptPerp+=track->Pt();

      if(fStudyPtCutInCone)
      {
        for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
        {
          if ( track->Pt() > fMinPtCutInCone[icut] ) 
          {
            coneptsumPerpTrackPerMinCut[icut]+=track->Pt();
            coneNPerpTrackPerMinCut    [icut]++;
          }          
        }
      }
      
      fhPtInPerpCone->Fill(ptTrig, track->Pt(), GetEventWeight());
            
      if(fStudyTracksInCone)
      {
        fhEtaPhiInPerpCone->Fill(track->Eta(),track->Phi(), GetEventWeight());

        ULong_t status = track->GetStatus();
        Bool_t okTOF = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
        Int_t trackBC = track->GetTOFBunchCrossing(bz);
        //Double32_t tof = track->GetTOFsignal()*1e-3;    
        
        if(okTOF && trackBC == 0)
        {
          fhPtInPerpConeTOFBC0->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeTOFBC0->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpBC0+=track->Pt();
        }

        Bool_t bConstrained = (!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1));
        //Bool_t bITSRefit    = (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
        if(!bConstrained) 
        {
          fhPtInPerpConeITSRefitOnSPDOn->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeITSRefitOnSPDOn->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpITSSPD+=track->Pt();
        }
        
        if(okTOF && trackBC == 0 && !bConstrained)
        {
          fhPtInPerpConeTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, track->Pt(), GetEventWeight());
          fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn->Fill(track->Eta(),track->Phi(), GetEventWeight());
          
          sumptPerpBC0ITSSPD+=track->Pt();
        }
        
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
            fhPtTrackInPerpConeDCA[0]->Fill(track->Pt(),  dca[0], GetEventWeight());
            fhPtTrackInPerpConeDCA[1]->Fill(track->Pt(),  dca[1], GetEventWeight());
          }
          else
          {
            fhPtTrackInPerpConeDCA[2]->Fill(track->Pt(), dcaCons, GetEventWeight());
          }
        } // trigger pt cut for DCA
        
      } // study tracks in cone
    } // r in cone
  } // track loop
  
  fhPerpConeSumPt           ->Fill(ptTrig, sumptPerp   , GetEventWeight());
  fhConeSumPtEtaBandUETrack ->Fill(ptTrig, etaBandPtSum, GetEventWeight());
  fhConeSumPtPhiBandUETrack ->Fill(ptTrig, phiBandPtSum, GetEventWeight());
  
  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fhPerpConeSumPtTrackPerMinPtCut->Fill(icut+1, coneptsumPerpTrackPerMinCut[icut], GetEventWeight());
      fhPerpConeNTrackPerMinPtCut    ->Fill(icut+1, coneNPerpTrackPerMinCut    [icut], GetEventWeight());
      
      if ( ptTrig > 10 ) 
      {
        fhPerpConeSumPtTrackPerMinPtCutLargePtTrig->Fill(icut+1, coneptsumPerpTrackPerMinCut[icut], GetEventWeight());   
        fhPerpConeNTrackPerMinPtCutLargePtTrig    ->Fill(icut+1, coneNPerpTrackPerMinCut    [icut], GetEventWeight());   
      }
    }
  }

  
  if(fFillUEBandSubtractHistograms > 1)
  {
    fhConeSumPtEtaBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSum *GetEventWeight()); // check
    fhConeSumPtPhiBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSum *GetEventWeight()); // check
  }
  
  if(fStudyTracksInCone) 
  {
    fhPerpConeSumPtTOFBC0         ->Fill(ptTrig, sumptPerpBC0   , GetEventWeight());
    fhPerpConeSumPtITSRefitOnSPDOn->Fill(ptTrig, sumptPerpITSSPD, GetEventWeight());
    fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, sumptPerpBC0ITSSPD, GetEventWeight());
  }

}

//_____________________________________________________________________________________________________________________________________
/// Normalize phi/eta band per area unit
//_____________________________________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate, Float_t coneptsumCluster,
                                                                  Float_t coneptsumCell,  Float_t coneptsumTrack,
                                                                  Float_t &sumEtaUESub ,  Float_t &sumPhiUESub, Int_t mcIndex)
{
  Float_t etaUEptsumTrack   = 0 ;
  Float_t phiUEptsumTrack   = 0 ;
  Float_t etaUEptsumCluster = 0 ;
  Float_t phiUEptsumCluster = 0 ;
  Float_t etaUEptsumCell    = 0 ;
  Float_t phiUEptsumCell    = 0 ;
  
  Int_t   partTypeInCone    = GetIsolationCut()->GetParticleTypeInCone();
  
  // Do the normalization
  
  Float_t conesize  = GetIsolationCut()->GetConeSize();
  Float_t coneA     = conesize*conesize*TMath::Pi(); // A = pi R^2, isolation cone area
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  Float_t pt        = pCandidate->Pt();
  Float_t m02       = pCandidate->GetM02();
  Int_t  mcTag      = pCandidate->GetTag() ;
  
  // ------ //
  // Tracks //
  // ------ //
  Float_t phiUEptsumTrackNorm  = 0 ;
  Float_t etaUEptsumTrackNorm  = 0 ;
  Float_t coneptsumTrackSubPhi = 0 ;
  Float_t coneptsumTrackSubEta = 0 ;
  Float_t coneptsumTrackSubPhiNorm = 0 ;
  Float_t coneptsumTrackSubEtaNorm = 0 ;
  
  if ( partTypeInCone != AliIsolationCut::kOnlyNeutral )
  {
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateTrackUEBand   (pCandidate,etaUEptsumTrack  ,phiUEptsumTrack  );// rajouter ici l'histo eta phi
    
    // Fill histograms
    if(fFillUEBandSubtractHistograms > 1)
    {
      fhConeSumPtVSUETracksEtaBand->Fill(coneptsumTrack, etaUEptsumTrack, GetEventWeight());
      fhConeSumPtVSUETracksPhiBand->Fill(coneptsumTrack, phiUEptsumTrack, GetEventWeight());
    }
    
    Float_t correctConeSumTrack    = 1;
    Float_t correctConeSumTrackPhi = 1;
    
    GetIsolationCut()->CalculateUEBandTrackNormalization(GetReader(),etaTrig, phiTrig,
                                                         phiUEptsumTrack,etaUEptsumTrack,
                                                         phiUEptsumTrackNorm,etaUEptsumTrackNorm,
                                                         correctConeSumTrack,correctConeSumTrackPhi);
    
    coneptsumTrackSubPhi = coneptsumTrack - phiUEptsumTrackNorm;
    coneptsumTrackSubEta = coneptsumTrack - etaUEptsumTrackNorm;
    
    fhConeSumPtPhiUESubTrack  ->Fill(ptTrig, coneptsumTrackSubPhi, GetEventWeight());
    fhConeSumPtEtaUESubTrack  ->Fill(ptTrig, coneptsumTrackSubEta, GetEventWeight());

    fhConeSumPtPhiUENormTrack ->Fill(ptTrig, phiUEptsumTrackNorm , GetEventWeight());
    fhConeSumPtEtaUENormTrack ->Fill(ptTrig, etaUEptsumTrackNorm , GetEventWeight());

    if(coneptsumTrack > 0)
    {
      coneptsumTrackSubPhiNorm = coneptsumTrackSubPhi/coneptsumTrack;
      coneptsumTrackSubEtaNorm = coneptsumTrackSubEta/coneptsumTrack;
    }
    
    if(fFillUEBandSubtractHistograms > 1)
    {       
      fhConeSumPtPhiUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubPhi *GetEventWeight()); // check
      fhConeSumPtEtaUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubEta *GetEventWeight()); // check
      
      fhFractionTrackOutConeEta          ->Fill(ptTrig ,         correctConeSumTrack-1, GetEventWeight());
      fhFractionTrackOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig,correctConeSumTrack-1 *GetEventWeight()); // check
      
      fhConeSumPtSubvsConeSumPtTotPhiTrack    ->Fill(coneptsumTrack, coneptsumTrackSubPhi    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotPhiTrack->Fill(coneptsumTrack, coneptsumTrackSubPhiNorm, GetEventWeight());
      fhConeSumPtSubvsConeSumPtTotEtaTrack    ->Fill(coneptsumTrack, coneptsumTrackSubEta    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotEtaTrack->Fill(coneptsumTrack, coneptsumTrackSubEtaNorm, GetEventWeight());
    }
  }
  
  // ------------------------ //
  // EMCal Clusters and cells //
  // ------------------------ //
  Float_t phiUEptsumClusterNorm  = 0 ;
  Float_t etaUEptsumClusterNorm  = 0 ;
  Float_t coneptsumClusterSubPhi = 0 ;
  Float_t coneptsumClusterSubEta = 0 ;
  Float_t coneptsumClusterSubPhiNorm = 0 ;
  Float_t coneptsumClusterSubEtaNorm = 0 ;
  Float_t phiUEptsumCellNorm     = 0 ;
  Float_t etaUEptsumCellNorm     = 0 ;
  Float_t coneptsumCellSubPhi    = 0 ;
  Float_t coneptsumCellSubEta    = 0 ;
  Float_t coneptsumCellSubPhiNorm = 0 ;
  Float_t coneptsumCellSubEtaNorm = 0 ;
  
  if ( partTypeInCone != AliIsolationCut::kOnlyCharged )
  {
    // -------------- //
    // EMCal clusters //
    // -------------- //
    
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateCaloUEBand    (pCandidate,etaUEptsumCluster,phiUEptsumCluster);// rajouter ici l'histo eta phi

    // Fill histograms
    if(fFillUEBandSubtractHistograms > 1)
    {
      fhConeSumPtVSUEClusterEtaBand->Fill(coneptsumCluster, etaUEptsumCluster, GetEventWeight());
      fhConeSumPtVSUEClusterPhiBand->Fill(coneptsumCluster, phiUEptsumCluster, GetEventWeight());
    }
    
    Float_t correctConeSumClusterEta = 1;
    Float_t correctConeSumClusterPhi = 1;
    
    GetIsolationCut()->CalculateUEBandClusterNormalization(GetReader(),etaTrig, phiTrig,
                                                           phiUEptsumCluster,etaUEptsumCluster,
                                                           phiUEptsumClusterNorm,etaUEptsumClusterNorm,
                                                           correctConeSumClusterEta,correctConeSumClusterPhi);
    
    // In case that cone is out of eta and phi side, we are over correcting, not too often with the current cuts ...
    // Comment if not used
    //  Float_t coneBadCellsCoeff   =1;
    //  Float_t etaBandBadCellsCoeff=1;
    //  Float_t phiBandBadCellsCoeff=1;
    //  GetIsolationCut()->GetCoeffNormBadCell(pCandidate,   GetReader(),coneBadCellsCoeff,etaBandBadCellsCoeff,phiBandBadCellsCoeff) ;
    
    //coneptsumCluster=coneptsumCluster*coneBadCellsCoeff*correctConeSumClusterEta*correctConeSumClusterPhi;
    
    coneptsumClusterSubPhi = coneptsumCluster - phiUEptsumClusterNorm;
    coneptsumClusterSubEta = coneptsumCluster - etaUEptsumClusterNorm;
        
    fhConeSumPtPhiUENormCluster          ->Fill(ptTrig ,          phiUEptsumClusterNorm , GetEventWeight());
    fhConeSumPtPhiUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubPhi, GetEventWeight());
    fhConeSumPtEtaUENormCluster          ->Fill(ptTrig ,          etaUEptsumClusterNorm , GetEventWeight());
    fhConeSumPtEtaUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubEta, GetEventWeight());
    
    if(coneptsumCluster!=0)
    {
      coneptsumClusterSubPhiNorm = coneptsumClusterSubPhi/coneptsumCluster;
      coneptsumClusterSubEtaNorm = coneptsumClusterSubEta/coneptsumCluster;
    }
    
    if(fFillUEBandSubtractHistograms > 1)
    {
      fhConeSumPtPhiUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubPhi *GetEventWeight()); // check
      fhConeSumPtEtaUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubEta *GetEventWeight()); // check
      
      fhFractionClusterOutConeEta          ->Fill(ptTrig ,          correctConeSumClusterEta-1, GetEventWeight());
      fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterEta-1 *GetEventWeight()); // check
      fhFractionClusterOutConePhi          ->Fill(ptTrig ,          correctConeSumClusterPhi-1, GetEventWeight());
      fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterPhi-1 *GetEventWeight()); // check
      
      fhConeSumPtSubvsConeSumPtTotPhiCluster    ->Fill(coneptsumCluster,coneptsumClusterSubPhi    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotPhiCluster->Fill(coneptsumCluster,coneptsumClusterSubPhiNorm, GetEventWeight());
      fhConeSumPtSubvsConeSumPtTotEtaCluster    ->Fill(coneptsumCluster,coneptsumClusterSubEta    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotEtaCluster->Fill(coneptsumCluster,coneptsumClusterSubEtaNorm, GetEventWeight());
    }

    // ----------- //
    // EMCal Cells //
    // ----------- //
    
    if ( fFillCellHistograms )
    {
      // Sum the pT in the phi or eta band for clusters or tracks
      CalculateCaloCellUEBand(pCandidate,etaUEptsumCell   ,phiUEptsumCell   );
      
      // Move to AliIsolationCut the calculation not the histograms??
      
      //Careful here if EMCal limits changed .. 2010 (4 SM) to 2011-12 (10 SM), for the moment consider 100 deg in phi
      Float_t emcEtaSize = 0.7*2; // TO FIX
      Float_t emcPhiSize = TMath::DegToRad()*100.; // TO FIX
      
      if(((2*conesize*emcPhiSize)-coneA)!=0)phiUEptsumCellNorm = phiUEptsumCell*(coneA / ((2*conesize*emcPhiSize)-coneA));
      if(((2*conesize*emcEtaSize)-coneA)!=0)etaUEptsumCellNorm = etaUEptsumCell*(coneA / ((2*conesize*emcEtaSize)-coneA));
      
      // Need to correct coneptsumCluster by the fraction of the cone out of the calorimeter cut acceptance!
      
      Float_t correctConeSumCellEta = 1;
      if(TMath::Abs(etaTrig)+conesize > emcEtaSize/2.)
      {
        Float_t excess = TMath::Abs(etaTrig) + conesize - emcEtaSize/2.;
        correctConeSumCellEta = GetIsolationCut()->CalculateExcessAreaFraction(excess);
        //printf("Excess EMC-Eta %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f\n",excess,coneA, excessA, angle*TMath::RadToDeg(), correctConeSumClusterEta);
        // Need to correct phi band surface if part of the cone falls out of track cut acceptance!
        if(((2*(conesize-excess)*emcPhiSize)-(coneA-correctConeSumCellEta))!=0)phiUEptsumCellNorm = phiUEptsumCell*(coneA / ((2*(conesize-excess)*emcPhiSize)-(coneA-correctConeSumCellEta)));
      }
      
      Float_t correctConeSumCellPhi = 1;
      //printf("EMCPhiTrig %2.2f, conesize %2.2f, sum %2.2f, rest %2.2f \n",phiTrig*TMath::RadToDeg(),conesize*TMath::RadToDeg(),(phiTrig+conesize)*TMath::RadToDeg(),(phiTrig-conesize)*TMath::RadToDeg() );
      if((phiTrig+conesize > 180*TMath::DegToRad()) ||
         (phiTrig-conesize <  80*TMath::DegToRad()))
      {
        Float_t excess = 0;
        if( phiTrig+conesize > 180*TMath::DegToRad() ) excess = conesize + phiTrig - 180*TMath::DegToRad() ;
        else                                           excess = conesize - phiTrig +  80*TMath::DegToRad() ;
        
        correctConeSumCellPhi = GetIsolationCut()->CalculateExcessAreaFraction(excess);
        //printf("Excess EMC-Phi %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f\n",excess,coneA, excessA, angle*TMath::RadToDeg(), correctConeSumClusterPhi);
        
        // Need to correct eta band surface if part of the cone falls out of track cut acceptance!
        if(((2*(conesize-excess)*emcEtaSize)-(coneA-correctConeSumCellPhi))!=0)etaUEptsumCellNorm = etaUEptsumCell*(coneA / ((2*(conesize-excess)*emcEtaSize)-(coneA-correctConeSumCellPhi)));
      }
      
      // In case that cone is out of eta and phi side, we are over correcting, not too often with the current cuts ...
      coneptsumCellSubPhi = coneptsumCell*correctConeSumCellEta*correctConeSumCellPhi - phiUEptsumCellNorm;
      coneptsumCellSubEta = coneptsumCell*correctConeSumCellEta*correctConeSumCellPhi - etaUEptsumCellNorm;
      
      fhConeSumPtPhiUESubCell           ->Fill(ptTrig ,          coneptsumCellSubPhi, GetEventWeight());
      fhConeSumPtPhiUESubCellTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumCellSubPhi *GetEventWeight()); // check
      fhConeSumPtEtaUESubCell           ->Fill(ptTrig ,          coneptsumCellSubEta, GetEventWeight());
      fhConeSumPtEtaUESubCellTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumCellSubEta *GetEventWeight()); // check
      
      fhFractionCellOutConeEta          ->Fill(ptTrig ,          correctConeSumCellEta-1, GetEventWeight());
      fhFractionCellOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumCellEta-1 *GetEventWeight()); // check
      fhFractionCellOutConePhi          ->Fill(ptTrig ,          correctConeSumCellPhi-1, GetEventWeight());
      fhFractionCellOutConePhiTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumCellPhi-1 *GetEventWeight()); // check
      if ( coneptsumCell > 0.01 )
      {
        coneptsumCellSubPhiNorm = coneptsumCellSubPhi/coneptsumCell;
        coneptsumCellSubEtaNorm = coneptsumCellSubEta/coneptsumCell;
      }
      
      fhConeSumPtSubvsConeSumPtTotPhiCell    ->Fill(coneptsumCell,coneptsumCellSubPhi    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotPhiCell->Fill(coneptsumCell,coneptsumCellSubPhiNorm, GetEventWeight());
      fhConeSumPtSubvsConeSumPtTotEtaCell    ->Fill(coneptsumCell,coneptsumCellSubEta    , GetEventWeight());
      fhConeSumPtSubNormvsConeSumPtTotEtaCell->Fill(coneptsumCell,coneptsumCellSubEtaNorm, GetEventWeight());
    }
  }
  
  sumPhiUESub = coneptsumClusterSubPhi + coneptsumTrackSubPhi;
  sumEtaUESub = coneptsumClusterSubEta + coneptsumTrackSubEta;
  
  if ( partTypeInCone == AliIsolationCut::kNeutralAndCharged )
  {
    // --------------------------- //
    // Tracks and clusters in cone //
    // --------------------------- //
    fhConeSumPtPhiUESub ->Fill(ptTrig,  sumPhiUESub, GetEventWeight());
    fhConeSumPtEtaUESub ->Fill(ptTrig,  sumEtaUESub, GetEventWeight());
    
    if(fFillUEBandSubtractHistograms > 1)
    {
      fhConeSumPtPhiUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumPhiUESub *GetEventWeight()); // check
      fhConeSumPtEtaUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumEtaUESub *GetEventWeight()); // check
      
      fhEtaBandClustervsTrack    ->Fill(etaUEptsumCluster    ,etaUEptsumTrack    , GetEventWeight());
      fhPhiBandClustervsTrack    ->Fill(phiUEptsumCluster    ,phiUEptsumTrack    , GetEventWeight());
      fhEtaBandNormClustervsTrack->Fill(etaUEptsumClusterNorm,etaUEptsumTrackNorm, GetEventWeight());
      fhPhiBandNormClustervsTrack->Fill(phiUEptsumClusterNorm,phiUEptsumTrackNorm, GetEventWeight());
      
      fhConeSumPtEtaUESubClustervsTrack->Fill(coneptsumClusterSubEta, coneptsumTrackSubEta, GetEventWeight());
      fhConeSumPtPhiUESubClustervsTrack->Fill(coneptsumClusterSubPhi, coneptsumTrackSubPhi, GetEventWeight());
      
      
      // Get the sum of pt in cone after UE subtraction
      // assign a bin to the candidate, depending on this quantity
      // see the shower shape in those bins.
      if(fFillBackgroundBinHistograms)
      {
        // Get the background bin for this cone and trigger
        Int_t ptsumAfterEtaBandUESubBin  = -1;
        
        AliDebug(1,Form("pT cand: %2.2f, In cone pT: Sum %2.2f, n bins %d",pt,coneptsumTrack+coneptsumCluster,fNBkgBin));
        
        for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
        {
          if( sumEtaUESub >= fBkgBinLimit[ibin] && sumEtaUESub < fBkgBinLimit[ibin+1]) ptsumAfterEtaBandUESubBin = ibin;
        }
        
        // Fill the histograms per bin of pt sum after UE subtraction in eta band
        
        if ( ptsumAfterEtaBandUESubBin  >= 0 )
        {
          AliDebug(1,Form("\t Sum bin %d [%2.2f,%2.2f]" , ptsumAfterEtaBandUESubBin ,fBkgBinLimit[ptsumAfterEtaBandUESubBin] ,fBkgBinLimit[ptsumAfterEtaBandUESubBin +1]));
          
          fhSumPtConeAfterEtaBandUESubBin[ptsumAfterEtaBandUESubBin]->Fill(pt, GetEventWeight());
          
          if(fFillSSHisto) fhSumPtConeAfterEtaBandUESubBinLambda0[ptsumAfterEtaBandUESubBin]->Fill(pt, m02, GetEventWeight());
        }
        
        if(IsDataMC())
        {
          Int_t  ptsumAfterEtaBandSubBinMC =  ptsumAfterEtaBandUESubBin+mcIndex*fNBkgBin;
          
          if( ptsumAfterEtaBandUESubBin  >=0 )
          {
            fhSumPtConeAfterEtaBandUESubBinMC [ ptsumAfterEtaBandSubBinMC]->Fill(pt, GetEventWeight());
            if(fFillSSHisto)  fhSumPtConeAfterEtaBandUESubBinLambda0MC [ ptsumAfterEtaBandSubBinMC]->Fill(pt, m02, GetEventWeight());
          }
          
          if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
          {
            ptsumAfterEtaBandSubBinMC  =  ptsumAfterEtaBandUESubBin+kmcPhoton*fNBkgBin;
            
            if( ptsumAfterEtaBandUESubBin  >=0 )
            {
              fhSumPtConeAfterEtaBandUESubBinMC [ ptsumAfterEtaBandSubBinMC]->Fill(pt, GetEventWeight());
              if(fFillSSHisto)  fhSumPtConeAfterEtaBandUESubBinLambda0MC [ ptsumAfterEtaBandSubBinMC]->Fill(pt, m02, GetEventWeight());
            }
          }
          
          // Check if decay and if pair is lost
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay )
            {
              ptsumAfterEtaBandSubBinMC  =  ptsumAfterEtaBandUESubBin+kmcPi0DecayLostPair*fNBkgBin;
            }
            else if(mcIndex == kmcEtaDecay)
            {
              ptsumAfterEtaBandSubBinMC  =  ptsumAfterEtaBandUESubBin+kmcEtaDecayLostPair*fNBkgBin;
            }
            else
              AliFatal(Form("Lost decay Bit assigned to bad case, mcIndex %d",mcIndex));
            
            if( ptsumAfterEtaBandUESubBin  >=0 )
            {
              fhSumPtConeAfterEtaBandUESubBinMC [ ptsumAfterEtaBandSubBinMC]->Fill(pt);
              if(fFillSSHisto)  fhSumPtConeBinLambda0MC [ ptsumAfterEtaBandSubBinMC]->Fill(pt, m02, GetEventWeight());
            }
            
          } // check decays with lost pairs
          
        } // MC data
      } // background dependent bins
      
    }
    
    // ------------------------ //
    // Tracks and cells in cone //
    // ------------------------ //
    
    if(fFillCellHistograms)
    {
      Double_t sumPhiUESubTrackCell = coneptsumCellSubPhi + coneptsumTrackSubPhi;
      Double_t sumEtaUESubTrackCell = coneptsumCellSubEta + coneptsumTrackSubEta;
      
      fhConeSumPtPhiUESubTrackCell          ->Fill(ptTrig ,          sumPhiUESubTrackCell, GetEventWeight());
      fhConeSumPtPhiUESubTrackCellTrigEtaPhi->Fill(etaTrig, phiTrig, sumPhiUESubTrackCell *GetEventWeight()); // check
      fhConeSumPtEtaUESubTrackCell          ->Fill(ptTrig ,          sumEtaUESubTrackCell, GetEventWeight());
      fhConeSumPtEtaUESubTrackCellTrigEtaPhi->Fill(etaTrig, phiTrig, sumEtaUESubTrackCell *GetEventWeight()); // check
      
      fhEtaBandCellvsTrack    ->Fill(etaUEptsumCell    , etaUEptsumTrack    , GetEventWeight());
      fhPhiBandCellvsTrack    ->Fill(phiUEptsumCell    , phiUEptsumTrack    , GetEventWeight());
      fhEtaBandNormCellvsTrack->Fill(etaUEptsumCellNorm, etaUEptsumTrackNorm, GetEventWeight());
      fhPhiBandNormCellvsTrack->Fill(phiUEptsumCellNorm, phiUEptsumTrackNorm, GetEventWeight());
      
      fhConeSumPtEtaUESubCellvsTrack->Fill(coneptsumCellSubEta, coneptsumTrackSubEta, GetEventWeight());
      fhConeSumPtPhiUESubCellvsTrack->Fill(coneptsumCellSubPhi, coneptsumTrackSubPhi, GetEventWeight());
    }
  }
}

//______________________________________________________________________________________________________________
/// Get the cluster pT or sum of pT in isolation cone.
//______________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                        Float_t & coneptsumCluster, Float_t & coneptLeadCluster)
{
  coneptLeadCluster = 0;
  coneptsumCluster  = 0;
  
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;

  Float_t ptTrig = aodParticle->Pt();
  
  // Recover reference arrays with clusters and tracks
  TObjArray * refclusters = aodParticle->GetObjArray(GetAODObjArrayName()+"Clusters");
  if(!refclusters)
  {
    fhConeSumPtCluster ->Fill(ptTrig, 0., GetEventWeight());
    fhConePtLeadCluster->Fill(ptTrig, 0., GetEventWeight());
    if(fStudyExoticTrigger && fIsExoticTrigger)
      fhConeSumPtClusterExoTrigger->Fill(ptTrig, 0., GetEventWeight());

    if(coneptLeadCluster > 0  || coneptsumCluster > 0) 
      AliError(Form("No ref tracks!!! sum %f, lead %f",coneptsumCluster,coneptLeadCluster));
    
    return ;
  }
  
  // Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t ptcone = 0;
  Float_t coneNClusterPerMinCut    [20];
  Float_t coneptsumClusterPerMinCut[20];
  Float_t coneptsumClusterPerMaxCut[20];
  Float_t coneptsumClusterPerRCut  [10];
  
  Float_t coneptsumClusterPerNCellCut[20];
  Float_t coneptsumClusterPerExoCut  [20];

  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++)
    {
      coneNClusterPerMinCut    [icut] = 0;
      coneptsumClusterPerMinCut[icut] = 0;
      coneptsumClusterPerMaxCut[icut] = 0;
    }
  }
  
  if(fStudyRCutInCone)
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      coneptsumClusterPerRCut[icut] = 0;
    }
  }
  
  Int_t ishsh = -1;
  if(fStudyNCellsCut)
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
  
  if(fStudyExoticTrigger)
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
    
    fhPtInCone       ->Fill(ptTrig, ptcone, GetEventWeight());
    fhPtClusterInCone->Fill(ptTrig, ptcone, GetEventWeight());
    
    if(fStudyExoticTrigger && fIsExoticTrigger)
    {
      fhPtInConeExoTrigger        ->Fill(ptTrig , ptcone, GetEventWeight());
      fhPtClusterInConeExoTrigger ->Fill(ptTrig , ptcone, GetEventWeight());
    }
    
    if(IsPileUpAnalysisOn())
    {
      if(GetReader()->IsPileUpFromSPD())               fhPtInConePileUp[0]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig, ptcone, GetEventWeight());
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig, ptcone, GetEventWeight());
    }
    
    if(IsHighMultiplicityAnalysisOn()) fhPtInConeCent->Fill(GetEventCentrality(), ptcone, GetEventWeight());
    
    coneptsumCluster+=ptcone;
    if(ptcone > coneptLeadCluster) coneptLeadCluster = ptcone;
    
    if(fStudyPtCutInCone)
    {
      for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
      {
        if ( ptcone > fMinPtCutInCone[icut] ) 
        {
          coneptsumClusterPerMinCut[icut]+=ptcone;
          coneNClusterPerMinCut    [icut]++;
        }
        
        if ( ptcone < fMaxPtCutInCone[icut] ) coneptsumClusterPerMaxCut[icut]+=ptcone;
      }
    }
    
    if(fStudyRCutInCone)
    {
      Float_t distance = GetIsolationCut()->Radius(aodParticle->Eta(), aodParticle->Phi(), fMomentum.Eta(), GetPhi(fMomentum.Phi()));
      for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
      {
        if ( distance < fRCutInCone[icut] ) 
        {
          coneptsumClusterPerRCut[icut]+=ptcone;
          fhPtClusterInConePerRCut->Fill(icut+1, ptcone, GetEventWeight());
          if(ptTrig > 10) fhPtClusterInConePerRCutLargePtTrig->Fill(icut+1, ptcone, GetEventWeight());
        }
      }
    }
    
    if(fStudyNCellsCut)
    {
      if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
        fhPtClusterInConePerNCellPerSM[ishsh]->Fill(ptcone, fTrigSupMod, fNCellsWithWeight);
      
      for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
      {
        if ( fNCellsWithWeight >= fNCellsInCandidate[icut] ) 
        {
          coneptsumClusterPerNCellCut[icut]+=ptcone;
          fhPtClusterInConePerNCellCut->Fill(icut+1, ptcone, GetEventWeight());
          if(ptTrig > 10) fhPtClusterInConePerNCellCutLargePtTrig->Fill(icut+1, ptcone, GetEventWeight());
        }
      }
    }
    
    if(fStudyExoticTrigger)
    {
      for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
      {
        if ( fClusterExoticity < fExoCutInCandidate[icut] ) 
        {
          coneptsumClusterPerExoCut[icut]+=ptcone;
          fhPtClusterInConePerExoCut->Fill(icut+1, ptcone, GetEventWeight());
          if(ptTrig > 10) fhPtClusterInConePerExoCutLargePtTrig->Fill(icut+1, ptcone, GetEventWeight());
        }
      }
    }
  }
  
  fhConeSumPtCluster ->Fill(ptTrig, coneptsumCluster , GetEventWeight());
  fhConePtLeadCluster->Fill(ptTrig, coneptLeadCluster, GetEventWeight());
  if(fStudyExoticTrigger && fIsExoticTrigger)
    fhConeSumPtClusterExoTrigger  ->Fill(ptTrig, coneptsumCluster  , GetEventWeight());
  
  aodParticle->SetNeutralLeadPtInCone(coneptLeadCluster);
  aodParticle->SetNeutralPtSumInCone(coneptsumCluster);
  
  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fhConeNClusterPerMinPtCut    ->Fill(icut, coneNClusterPerMinCut    [icut], GetEventWeight());
      fhConeSumPtClusterPerMinPtCut->Fill(icut, coneptsumClusterPerMinCut[icut], GetEventWeight());
      fhConeSumPtClusterPerMaxPtCut->Fill(icut, coneptsumClusterPerMaxCut[icut], GetEventWeight());

      if ( ptTrig > 10 ) 
      {
        fhConeNClusterPerMinPtCutLargePtTrig    ->Fill(icut, coneNClusterPerMinCut    [icut], GetEventWeight());
        fhConeSumPtClusterPerMinPtCutLargePtTrig->Fill(icut, coneptsumClusterPerMinCut[icut], GetEventWeight());
        fhConeSumPtClusterPerMaxPtCutLargePtTrig->Fill(icut, coneptsumClusterPerMaxCut[icut], GetEventWeight());
      }
    }
  }

  if(fStudyRCutInCone)
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      fhConeSumPtClusterPerRCut->Fill(icut, coneptsumClusterPerRCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtClusterPerRCutLargePtTrig->Fill(icut, coneptsumClusterPerRCut[icut], GetEventWeight());      
    }
  }

  if(fStudyNCellsCut)
  {     
    if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
      fhConeSumPtClusterPerNCellPerSM[ishsh]->Fill(coneptsumCluster, fTrigSupMod, fNCellsWithWeight);

    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      fhConeSumPtClusterPerNCellCut->Fill(icut, coneptsumClusterPerNCellCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtClusterPerNCellCutLargePtTrig->Fill(icut, coneptsumClusterPerNCellCut[icut], GetEventWeight());      
    }
  }
  
  if(fStudyExoticTrigger)
  { 
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      fhConeSumPtClusterPerExoCut->Fill(icut, coneptsumClusterPerExoCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtClusterPerExoCutLargePtTrig->Fill(icut, coneptsumClusterPerExoCut[icut], GetEventWeight());      
    }
  }  
}

//______________________________________________________________________________________________________
/// Get the cell amplityde or sum of amplitudes in isolation cone.
/// Missing: Remove signal cells in cone in case the trigger is a cluster.
//______________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloCellSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                            Float_t & coneptsumCell)
{
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
  if( !fFillCellHistograms ) return ;
  
  Float_t conesize = GetIsolationCut()->GetConeSize();
  
  Float_t  ptTrig  = aodParticle->Pt();
  Float_t  phiTrig = aodParticle->Phi();
  if(phiTrig<0) phiTrig += TMath::TwoPi();
  Float_t  etaTrig = aodParticle->Eta();
  
  if(aodParticle->GetDetectorTag()==kEMCAL)
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    Int_t absId = -999;
    
    if (eGeom->GetAbsCellIdFromEtaPhi(etaTrig,phiTrig,absId))
    {
      if(!eGeom->CheckAbsCellId(absId)) return ;
      
      // Get absolute (col,row) of trigger particle
      Int_t nSupMod = eGeom->GetSuperModuleNumber(absId);
      Int_t nModule = -1;
      Int_t imEta=-1, imPhi=-1;
      Int_t ieta =-1, iphi =-1;
      
      if (eGeom->GetCellIndex(absId,nSupMod,nModule,imPhi,imEta))
      {
        Int_t iEta=-1, iPhi=-1;
        eGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,imPhi,imEta,iphi,ieta);
        
        Int_t colTrig = iEta;
        if (nSupMod % 2) colTrig = AliEMCALGeoParams::fgkEMCALCols + iEta ;
        Int_t rowTrig = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
        
        Int_t sqrSize = int(conesize/0.0143);
        
        AliVCaloCells * cells = GetEMCALCells();
        
        // Loop on cells in cone
        for(Int_t irow = rowTrig-sqrSize; irow < rowTrig+sqrSize; irow++)
        {
          for(Int_t icol = colTrig-sqrSize; icol < colTrig+sqrSize; icol++)
          {
            Int_t inSector = int(irow/AliEMCALGeoParams::fgkEMCALRows);
            if(inSector==5) continue;
            
            Int_t inSupMod = -1;
            Int_t icolLoc  = -1;
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              inSupMod = 2*inSector + 1;
              icolLoc  = icol;
            }
            else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1)
            {
              inSupMod = 2*inSector;
              icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
            }
            
            Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*inSector ;
            
            Int_t iabsId = eGeom->GetAbsCellIdFromCellIndexes(inSupMod,irowLoc,icolLoc);
            if(!eGeom->CheckAbsCellId(iabsId)) continue;
            
            fhPtCellInCone->Fill(ptTrig, cells->GetCellAmplitude(iabsId), GetEventWeight());
            coneptsumCell += cells->GetCellAmplitude(iabsId);
          }
        }
      }
    }
  }
  
  fhConeSumPtCell->Fill(ptTrig, coneptsumCell, GetEventWeight());
}

//___________________________________________________________________________________________________________
/// Get the track pT or sum of pT in isolation cone.
//___________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateTrackSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                         Float_t & coneptsumTrack, Float_t & coneptLeadTrack)
{
  coneptLeadTrack = 0;
  coneptsumTrack  = 0;
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  Float_t  ptTrig = aodParticle->Pt();
  
  // Recover reference arrays with clusters and tracks
  TObjArray * reftracks   = aodParticle->GetObjArray(GetAODObjArrayName()+"Tracks");
  if(!reftracks)
  {
    fhConeSumPtTrack      ->Fill(ptTrig, 0., GetEventWeight());
    if(fStudyExoticTrigger && fIsExoticTrigger)
      fhConeSumPtTrackExoTrigger->Fill(ptTrig, 0., GetEventWeight());
    
    if(fStudyTracksInCone)
    {
      fhConeSumPtTrackTOFNo ->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackTOFBC0->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackTOFBCN->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackITSRefitOnSPDOn  ->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackITSRefitOffSPDOff->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackITSRefitOnSPDOff ->Fill(ptTrig, 0., GetEventWeight());
      fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, 0., GetEventWeight());
    }
    
    fhConePtLeadTrack     ->Fill(ptTrig, 0., GetEventWeight());
    
    if(coneptLeadTrack > 0  || coneptsumTrack > 0) 
      AliError(Form("No ref tracks!!! sum %f, lead %f",coneptsumTrack,coneptLeadTrack));
    return ;
  }
  
  Double_t bz = GetReader()->GetInputEvent()->GetMagneticField();
  
  Float_t pTtrack  = 0;
  Float_t phitrack = 0;
  Float_t etatrack = 0;
  Float_t coneNTrackPerMinCut    [20];
  Float_t coneptsumTrackPerMinCut[20];
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

  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      coneptsumTrackPerMinCut[icut] = 0;
      coneptsumTrackPerMaxCut[icut] = 0;
      coneNTrackPerMinCut    [icut] = 0;
    }
  }
  
  if(fStudyEtaCutInCone)
  {
    for(Int_t icut = 0; icut < fNEtaCutsInCone; icut++) 
    {
      coneptsumTrackPerEtaCut[icut] = 0;
    }
  }

  if(fStudyRCutInCone)
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      coneptsumTrackPerRCut[icut] = 0;
    }
  }

  Int_t ishsh = -1;
  if(fStudyNCellsCut)
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
  
  if(fStudyExoticTrigger)
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
    
    fhPtInCone      ->Fill(ptTrig , pTtrack, GetEventWeight());
    fhPtTrackInCone ->Fill(ptTrig , pTtrack, GetEventWeight());
    
    if( IsDataMC() && GetMC() )
    {
      Int_t trackLabel = TMath::Abs(track->GetLabel());
      
      AliVParticle * mcpart = GetMC()->GetTrack(trackLabel);
      if( !mcpart ) continue;
      
      Int_t  partInConeCharge = TMath::Abs(mcpart->Charge());
      Int_t  partInConePDG    = mcpart->PdgCode();
      Bool_t physPrimary      = mcpart->IsPhysicalPrimary();
      
      if ( partInConeCharge > 0 &&  TMath::Abs(partInConePDG) != 11 ) // exclude electrons and neutrals
      {
        Int_t mcChTag = 3;
        if      ( TMath::Abs(partInConePDG) == 211  )  mcChTag = 0;
        else if ( TMath::Abs(partInConePDG) == 321  )  mcChTag = 1; 
        else if ( TMath::Abs(partInConePDG) == 2212 )  mcChTag = 2; 
        if(physPrimary)
          fhPtTrackInConeMCPrimary  [mcChTag]->Fill(ptTrig , pTtrack, GetEventWeight());
        else
          fhPtTrackInConeMCSecondary[mcChTag]->Fill(ptTrig , pTtrack, GetEventWeight());
      }
    }
    
    if(fStudyExoticTrigger && fIsExoticTrigger)
    {
      fhPtInConeExoTrigger      ->Fill(ptTrig , pTtrack, GetEventWeight());
      fhPtTrackInConeExoTrigger ->Fill(ptTrig , pTtrack, GetEventWeight());
    }
    
    if(IsHighMultiplicityAnalysisOn()) fhPtInConeCent->Fill(GetEventCentrality(), pTtrack, GetEventWeight());
    
    coneptsumTrack+=pTtrack;
    if(pTtrack > coneptLeadTrack) coneptLeadTrack = pTtrack;
    
    if(fStudyPtCutInCone)
    {
      for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
      {
        if ( pTtrack > fMinPtCutInCone[icut] ) 
        {
          coneptsumTrackPerMinCut[icut]+=pTtrack;
          coneNTrackPerMinCut    [icut]++;
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
          fhPtTrackInConePerRCut->Fill(icut+1, pTtrack, GetEventWeight());
          if(ptTrig > 10) fhPtTrackInConePerRCutLargePtTrig->Fill(icut+1, pTtrack, GetEventWeight());
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
          fhPtTrackInConePerNCellCut->Fill(icut+1, pTtrack, GetEventWeight());
          if(ptTrig > 10) fhPtTrackInConePerNCellCutLargePtTrig->Fill(icut+1, pTtrack, GetEventWeight());
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
          fhPtTrackInConePerExoCut->Fill(icut+1, pTtrack, GetEventWeight());
          if(ptTrig > 10) fhPtTrackInConePerExoCutLargePtTrig->Fill(icut+1, pTtrack, GetEventWeight());
        }
      }
    }
    
    Bool_t okTOF = kFALSE ;
    Int_t trackBC = 0;
    if(fStudyTracksInCone)
    {
      phitrack = track->Phi();
      etatrack = track->Eta();

      fhPhiTrackInCone->Fill(pTtrack, phitrack, GetEventWeight());
      fhEtaTrackInCone->Fill(pTtrack, etatrack, GetEventWeight());
      fhEtaPhiTrackInCone->Fill(etatrack, phitrack, GetEventWeight());
      
      // TOF
      ULong_t status = track->GetStatus();
      okTOF   = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
      trackBC = track->GetTOFBunchCrossing(bz);
      Double32_t tof = track->GetTOFsignal()*1e-3;    
      
      Int_t vtxBC = GetReader()->GetVertexBC();
      if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA) fhPtTrackInConeVtxBC0->Fill(ptTrig, pTtrack, GetEventWeight());
      
      if(okTOF)
      {
        fhTrackTOFInCone->Fill(pTtrack,tof,GetEventWeight());
        
        if(fStudyExoticTrigger && fIsExoticTrigger)
          fhTrackTOFInConeExoTrigger->Fill(pTtrack,tof,GetEventWeight());
        
        if(trackBC == 0) 
        {
          fhPtTrackInConeTOFBC0 ->Fill(ptTrig , pTtrack , GetEventWeight());
          fhPhiTrackInConeTOFBC0->Fill(pTtrack, phitrack, GetEventWeight());
          fhEtaTrackInConeTOFBC0->Fill(pTtrack, etatrack, GetEventWeight());
          fhEtaPhiTrackInConeTOFBC0->Fill(etatrack, phitrack, GetEventWeight());
          fhTrackTOFInConeBC0   ->Fill(pTtrack, tof     , GetEventWeight());
          coneptsumTrackTOFBC0 += pTtrack;
        }
        else 
        {
          fhPtTrackInConeTOFBCN ->Fill(ptTrig , pTtrack , GetEventWeight());
          fhPhiTrackInConeTOFBCN->Fill(pTtrack, phitrack, GetEventWeight());
          fhEtaTrackInConeTOFBCN->Fill(pTtrack, etatrack, GetEventWeight());
          fhEtaPhiTrackInConeTOFBCN->Fill(etatrack, phitrack, GetEventWeight());
          coneptsumTrackTOFBCN += pTtrack;
        }
      }
      else
      {
        fhPtTrackInConeTOFNo ->Fill(ptTrig, pTtrack , GetEventWeight());
        fhPhiTrackInConeTOFNo->Fill(ptTrig, phitrack, GetEventWeight());
        fhEtaTrackInConeTOFNo->Fill(ptTrig, etatrack, GetEventWeight());
        fhEtaPhiTrackInConeTOFNo->Fill(etatrack, phitrack, GetEventWeight());
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
          fhPtTrackInConeITSRefitOnSPDOff ->Fill(ptTrig, pTtrack , GetEventWeight());
          fhPhiTrackInConeITSRefitOnSPDOff->Fill(ptTrig, phitrack, GetEventWeight());
          fhEtaTrackInConeITSRefitOnSPDOff->Fill(ptTrig, etatrack, GetEventWeight());
          fhEtaPhiTrackInConeITSRefitOnSPDOff->Fill(etatrack, phitrack, GetEventWeight());
        }
        else
        {
          coneptsumTrackITSRefitOffSPDOff += pTtrack;
          fhPtTrackInConeITSRefitOffSPDOff ->Fill(ptTrig, pTtrack , GetEventWeight());
          fhPhiTrackInConeITSRefitOffSPDOff->Fill(ptTrig, phitrack, GetEventWeight());
          fhEtaTrackInConeITSRefitOffSPDOff->Fill(ptTrig, etatrack, GetEventWeight());
          fhEtaPhiTrackInConeITSRefitOffSPDOff->Fill(etatrack, phitrack, GetEventWeight());
        }
      }
      else
      {
        coneptsumTrackITSRefitOnSPDOn   += pTtrack;
        fhPtTrackInConeITSRefitOnSPDOn ->Fill(ptTrig, pTtrack , GetEventWeight());
        fhPhiTrackInConeITSRefitOnSPDOn->Fill(ptTrig, phitrack, GetEventWeight());
        fhEtaTrackInConeITSRefitOnSPDOn->Fill(ptTrig, etatrack, GetEventWeight());
        fhEtaPhiTrackInConeITSRefitOnSPDOn->Fill(etatrack, phitrack, GetEventWeight());
      }
      
      if(okTOF && trackBC == 0 && !bConstrained)
      {
        fhPtTrackInConeTOFBC0ITSRefitOnSPDOn ->Fill(ptTrig , pTtrack , GetEventWeight());
        fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(pTtrack, phitrack, GetEventWeight());
        fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(pTtrack, etatrack, GetEventWeight());
        fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn->Fill(etatrack, phitrack, GetEventWeight());
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
          fhPtTrackInConeDCA[0]->Fill(pTtrack,  dca[0], GetEventWeight());
          fhPtTrackInConeDCA[1]->Fill(pTtrack,  dca[1], GetEventWeight());
        }
        else
        {
          fhPtTrackInConeDCA[2]->Fill(pTtrack, dcaCons, GetEventWeight());
        }
      } // trigger pt cut
    }
    
    if(IsPileUpAnalysisOn())
    {
      if(GetReader()->IsPileUpFromSPD())             
      {  
        fhPtInConePileUp[0]            ->Fill(ptTrig, pTtrack, GetEventWeight());
        if(fStudyTracksInCone)
        {
          if(okTOF && trackBC!=0 )                         fhPtTrackInConeOtherBCPileUpSPD->Fill(ptTrig, pTtrack, GetEventWeight());
          if(okTOF && trackBC==0 )                         fhPtTrackInConeBC0PileUpSPD    ->Fill(ptTrig, pTtrack, GetEventWeight()); 
        }
      }
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig, pTtrack, GetEventWeight());
    }
  }

  fhConeSumPtTrack ->Fill(ptTrig, coneptsumTrack , GetEventWeight());
  fhConePtLeadTrack->Fill(ptTrig, coneptLeadTrack, GetEventWeight());
  if(fStudyExoticTrigger && fIsExoticTrigger)
    fhConeSumPtTrackExoTrigger  ->Fill(ptTrig, coneptsumTrack  , GetEventWeight());
  
  if(fStudyTracksInCone)
  {
    fhConeSumPtTrackTOFBC0->Fill(ptTrig, coneptsumTrackTOFBC0, GetEventWeight());
    fhConeSumPtTrackTOFBCN->Fill(ptTrig, coneptsumTrackTOFBCN, GetEventWeight());
    fhConeSumPtTrackTOFNo ->Fill(ptTrig, coneptsumTrackTOFNo , GetEventWeight());

    fhConeSumPtTrackITSRefitOnSPDOn  ->Fill(ptTrig, coneptsumTrackITSRefitOnSPDOn  , GetEventWeight());
    fhConeSumPtTrackITSRefitOffSPDOff->Fill(ptTrig, coneptsumTrackITSRefitOffSPDOff, GetEventWeight());
    fhConeSumPtTrackITSRefitOnSPDOff ->Fill(ptTrig, coneptsumTrackITSRefitOnSPDOff , GetEventWeight());
    
    fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn->Fill(ptTrig, coneptsumTrackTOFBC0ITSRefitOnSPDOn, GetEventWeight());
  }
  
  aodParticle->SetChargedLeadPtInCone(coneptLeadTrack);
  aodParticle->SetChargedPtSumInCone(coneptsumTrack);
  
  if(fStudyPtCutInCone)
  {
    for(Int_t icut = 0; icut < fNPtCutsInCone; icut++) 
    {
      fhConeSumPtTrackPerMinPtCut->Fill(icut+1, coneptsumTrackPerMinCut[icut], GetEventWeight());
      fhConeSumPtTrackPerMaxPtCut->Fill(icut+1, coneptsumTrackPerMaxCut[icut], GetEventWeight());
      fhConeNTrackPerMinPtCut    ->Fill(icut+1, coneNTrackPerMinCut    [icut], GetEventWeight());

      if ( ptTrig > 10 ) 
      {
        fhConeSumPtTrackPerMinPtCutLargePtTrig->Fill(icut+1, coneptsumTrackPerMinCut[icut], GetEventWeight());   
        fhConeSumPtTrackPerMaxPtCutLargePtTrig->Fill(icut+1, coneptsumTrackPerMaxCut[icut], GetEventWeight());
        fhConeNTrackPerMinPtCutLargePtTrig    ->Fill(icut+1, coneNTrackPerMinCut    [icut], GetEventWeight());   
      }
    }
  }
  
  if(fStudyEtaCutInCone)
  {
    for(Int_t icut = 0; icut < fNEtaCutsInCone; icut++) 
    {
      fhConeSumPtTrackPerEtaCut ->Fill(icut+1, coneptsumTrackPerEtaCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtTrackPerEtaCutLargePtTrig->Fill(icut+1, coneptsumTrackPerEtaCut[icut], GetEventWeight());
    }
  }
  
  if(fStudyRCutInCone)
  {
    for(Int_t icut = 0; icut < fNRCutsInCone; icut++) 
    {
      fhConeSumPtTrackPerRCut ->Fill(icut+1, coneptsumTrackPerRCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtTrackPerRCutLargePtTrig->Fill(icut+1, coneptsumTrackPerRCut[icut], GetEventWeight());
    }
  }
  
  if(fStudyNCellsCut)
  {
    if ( ptTrig > 8 && ptTrig < 12 && ishsh >=0 )
      fhConeSumPtTrackPerNCellPerSM[ishsh]->Fill(coneptsumTrack, fTrigSupMod, fNCellsWithWeight);

    for(Int_t icut = 0; icut < fNNCellsInCandidate; icut++) 
    {
      fhConeSumPtTrackPerNCellCut->Fill(icut, coneptsumTrackPerNCellCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtTrackPerNCellCutLargePtTrig->Fill(icut, coneptsumTrackPerNCellCut[icut], GetEventWeight());      
    }
  }
  
  if(fStudyExoticTrigger)
  { 
    for(Int_t icut = 0; icut < fNExoCutInCandidate; icut++) 
    {
      fhConeSumPtTrackPerExoCut->Fill(icut, coneptsumTrackPerExoCut[icut], GetEventWeight());
      if ( ptTrig > 10 ) fhConeSumPtTrackPerExoCutLargePtTrig->Fill(icut, coneptsumTrackPerExoCut[icut], GetEventWeight());      
    }
  }
}

//_____________________________________________________________________________
/// Fill some histograms to understand pile-up.
//_____________________________________________________________________________
void AliAnaParticleIsolation::FillPileUpHistograms(Float_t energy, Float_t time)
{  
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

//_____________________________________________________________________________________________________________________
/// Fill Track matching and Shower Shape control histograms.
//_____________________________________________________________________________________________________________________
void AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms
(AliAODPWG4ParticleCorrelation  *pCandidate,
 Float_t coneptsum,  Float_t coneptsumTrack, Float_t coneptsumClust, 
 Float_t coneleadpt, Int_t mcIndex)
{
  if(!fFillTMHisto && !fFillSSHisto && !fFillBackgroundBinHistograms && !fFillTaggedDecayHistograms) return;
  
  Int_t  nMaxima   = pCandidate->GetNLM();
  Int_t  mcTag     = pCandidate->GetTag() ;
  Bool_t isolated  = pCandidate->IsIsolated();
  
  Float_t m02    = pCandidate->GetM02() ;
  Float_t energy = pCandidate->E();
  Float_t pt     = pCandidate->Pt();
  Float_t eta    = pCandidate->Eta();
  Float_t phi    = pCandidate->Phi();
  if(phi<0) phi+= TMath::TwoPi();
  
  // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
  //
  if(fFillTaggedDecayHistograms)
  {
    Int_t decayTag = pCandidate->DecayTag();
    if(decayTag < 0) decayTag = 0;
    
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      if(!GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit])) continue;
      
      if(fFillSSHisto) fhPtLambda0Decay[isolated][ibit]->Fill(pt, m02, GetEventWeight());
      
      // In case it was not done on the trigger selection task
      // apply here a shower shape cut to select photons
      if( m02 > fDecayTagsM02Cut ) continue;
      
      fhPtDecay    [isolated][ibit]->Fill(pt,       GetEventWeight());
      fhEtaPhiDecay[isolated][ibit]->Fill(eta, phi, GetEventWeight());
      
      if(IsDataMC())
      {
        fhPtDecayMC[isolated][ibit][mcIndex]->Fill(pt, GetEventWeight());
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
          fhPtDecayMC[isolated][ibit][kmcPhoton]->Fill(pt, GetEventWeight());
        
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtDecayMC[isolated][ibit][kmcPi0DecayLostPair]->Fill(pt, GetEventWeight());
          else if( mcIndex == kmcEtaDecay ) fhPtDecayMC[isolated][ibit][kmcEtaDecayLostPair]->Fill(pt, GetEventWeight());
        }
      }
    } // bit loop
  } // decay histograms
  
  // Get the max pt leading in cone or the sum of pt in cone
  // assign a bin to the candidate, depending on both quantities
  // see the shower shape in those bins.
  if(fFillBackgroundBinHistograms)
  {
    // Get the background bin for this cone and trigger
    Int_t ptsumBin  = -1;
    Int_t leadptBin = -1;
    
    AliDebug(1,Form("pT cand: %2.2f, In cone pT: Sum %2.2f, Lead %2.2f, n bins %d",pt,coneptsum,coneleadpt,fNBkgBin));
    
    for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
    {
      if( coneptsum  >= fBkgBinLimit[ibin] && coneptsum  < fBkgBinLimit[ibin+1]) ptsumBin  = ibin;
      if( coneleadpt >= fBkgBinLimit[ibin] && coneleadpt < fBkgBinLimit[ibin+1]) leadptBin = ibin;
    }
    
    // Fill the histograms per bin of pt lead or pt sum
    
    if ( leadptBin >= 0 )
    {
      AliDebug(1,Form("\t Lead bin %d [%2.2f,%2.2f]", leadptBin,fBkgBinLimit[leadptBin],fBkgBinLimit[leadptBin+1]));
      
      fhPtLeadConeBin[leadptBin]->Fill(pt, GetEventWeight());
      
      if(fFillSSHisto)
        fhPtLeadConeBinLambda0[leadptBin]->Fill(pt, m02, GetEventWeight());
      
      if ( leadptBin == 0 )
        AliDebug(1,Form("No track/clusters in isolation cone: cand pt %2.2f GeV/c, track multiplicity %d, N clusters %d",
                        pt, GetTrackMultiplicity(),GetEMCALClusters()->GetEntriesFast()));
    }
    
    if ( ptsumBin  >= 0 )
    {
      AliDebug(1,Form("\t Sum bin %d [%2.2f,%2.2f]" , ptsumBin ,fBkgBinLimit[ptsumBin] ,fBkgBinLimit[ptsumBin +1]));
      
      fhSumPtConeBin[ptsumBin]->Fill(pt, GetEventWeight());
      
      if(fFillSSHisto) fhSumPtConeBinLambda0[ptsumBin]->Fill(pt, m02, GetEventWeight());
    }
    
    // Check if it was a decay
    if( fFillTaggedDecayHistograms && m02 < fDecayTagsM02Cut )
    {
      Int_t decayTag = pCandidate->DecayTag();
      if(decayTag < 0) decayTag = 0;
      
      for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
      {
        if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
        {
          Int_t leadptBinDecay = leadptBin+ibit*fNBkgBin;
          Int_t  ptsumBinDecay =  ptsumBin+ibit*fNBkgBin;
          if( leadptBin >=0 ) fhPtLeadConeBinDecay[leadptBinDecay]->Fill(pt, GetEventWeight());
          if( ptsumBin  >=0 ) fhSumPtConeBinDecay [ ptsumBinDecay]->Fill(pt, GetEventWeight());
        }
      }
    }
    
    if(IsDataMC())
    {
      Int_t leadptBinMC = leadptBin+mcIndex*fNBkgBin;
      Int_t  ptsumBinMC =  ptsumBin+mcIndex*fNBkgBin;
      
      if( leadptBin >=0 )
      {
        fhPtLeadConeBinMC[leadptBinMC]->Fill(pt, GetEventWeight());
        if(fFillSSHisto) fhPtLeadConeBinLambda0MC[leadptBinMC]->Fill(pt, m02, GetEventWeight());
      }
      
      if( ptsumBin  >=0 )
      {
        fhSumPtConeBinMC [ ptsumBinMC]->Fill(pt, GetEventWeight());
        if(fFillSSHisto)  fhSumPtConeBinLambda0MC [ ptsumBinMC]->Fill(pt, m02, GetEventWeight());
      }
      
      if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
      {
        leadptBinMC = leadptBin+kmcPhoton*fNBkgBin;
        ptsumBinMC  =  ptsumBin+kmcPhoton*fNBkgBin;
        if( leadptBin >=0 )
        {
          fhPtLeadConeBinMC[leadptBinMC]->Fill(pt, GetEventWeight());
          if(fFillSSHisto) fhPtLeadConeBinLambda0MC[leadptBinMC]->Fill(pt, m02, GetEventWeight());
        }
        
        if( ptsumBin  >=0 )
        {
          fhSumPtConeBinMC [ ptsumBinMC]->Fill(pt, GetEventWeight());
          if(fFillSSHisto)  fhSumPtConeBinLambda0MC [ ptsumBinMC]->Fill(pt, m02, GetEventWeight());
        }
      }
      
      // Check if decay and if pair is lost
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
      {
        if     ( mcIndex == kmcPi0Decay )
        {
          leadptBinMC = leadptBin+kmcPi0DecayLostPair*fNBkgBin;
          ptsumBinMC  =  ptsumBin+kmcPi0DecayLostPair*fNBkgBin;
        }
        else if(mcIndex == kmcEtaDecay)
        {
          leadptBinMC = leadptBin+kmcEtaDecayLostPair*fNBkgBin;
          ptsumBinMC  =  ptsumBin+kmcEtaDecayLostPair*fNBkgBin;
        }
        else
          AliFatal(Form("Lost decay Bit assigned to bad case, mcIndex %d",mcIndex));
        
        if( leadptBin >=0 )
        {
          fhPtLeadConeBinMC[leadptBinMC]->Fill(pt, GetEventWeight());
          if(fFillSSHisto) fhPtLeadConeBinLambda0MC[leadptBinMC]->Fill(pt, m02, GetEventWeight());
        }
        
        if( ptsumBin  >=0 )
        {
          fhSumPtConeBinMC [ ptsumBinMC]->Fill(pt);
          if(fFillSSHisto)  fhSumPtConeBinLambda0MC [ ptsumBinMC]->Fill(pt, m02, GetEventWeight());
        }
        
      } // check decays with lost pairs
      
    } // MC data
  } // background dependent bins

  // Fill histograms on selected pt bins of the trigger particle
  Int_t ptTrigBin  = -1;
  if(fFillPtTrigBinHistograms)
  {
    for(Int_t ibin = 0; ibin < fNPtTrigBin; ibin++)
    {
      if( pt  >= fPtTrigBinLimit[ibin] && pt < fPtTrigBinLimit[ibin+1]) ptTrigBin  = ibin;
    }
    
    // Fill the histograms if the bin is found.
    if ( ptTrigBin >= 0 )
    {
      AliDebug(1,Form("Trigger pT %f, bin %d [%2.2f,%2.2f]",pt,ptTrigBin,fPtTrigBinLimit[ptTrigBin],fPtTrigBinLimit[ptTrigBin+1]));
      
      fhPtTrigBinPtLeadCone[ptTrigBin]->Fill(coneleadpt, GetEventWeight());
      fhPtTrigBinSumPtCone [ptTrigBin]->Fill(coneptsum , GetEventWeight());
      
      if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
      {
        fhPtTrigBinSumPtTrackCone  [ptTrigBin]->Fill(coneptsumTrack, GetEventWeight());
        fhPtTrigBinSumPtClusterCone[ptTrigBin]->Fill(coneptsumClust, GetEventWeight());
      }
      
      if(fFillSSHisto)
      {
        fhPtTrigBinLambda0vsPtLeadCone[ptTrigBin]->Fill(coneleadpt, m02, GetEventWeight());
        fhPtTrigBinLambda0vsSumPtCone [ptTrigBin]->Fill(coneptsum , m02, GetEventWeight());
        
        if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
        {
          fhPtTrigBinLambda0vsSumPtTrackCone  [ptTrigBin]->Fill(coneptsumTrack, m02, GetEventWeight());
          fhPtTrigBinLambda0vsSumPtClusterCone[ptTrigBin]->Fill(coneptsumClust, m02, GetEventWeight());
        }
      }
      
      // Check if it was a decay
      if( fFillTaggedDecayHistograms && m02 < fDecayTagsM02Cut )
      {
        Int_t decayTag = pCandidate->DecayTag();
        if(decayTag < 0) decayTag = 0;
        
        for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
          {
            Int_t binDecay = ptTrigBin+ibit*fNPtTrigBin;
            if( binDecay > 0 )
            {
              fhPtTrigBinPtLeadConeDecay[binDecay]->Fill(coneleadpt, GetEventWeight());
              fhPtTrigBinSumPtConeDecay [binDecay]->Fill(coneptsum , GetEventWeight());
              
              if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
              {
                fhPtTrigBinSumPtTrackConeDecay  [binDecay]->Fill(coneptsumTrack, GetEventWeight());
                fhPtTrigBinSumPtClusterConeDecay[binDecay]->Fill(coneptsumClust, GetEventWeight());
              }
            }
          }
        }
      } // decay
      
      if( IsDataMC() )
      {
        Int_t ptTrigBinMC = ptTrigBin+mcIndex*fNPtTrigBin;
        
        fhPtTrigBinPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt    , GetEventWeight());
        fhPtTrigBinSumPtConeMC [ptTrigBinMC]->Fill(coneptsum     , GetEventWeight());
       
        if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
        {
          fhPtTrigBinSumPtTrackConeMC  [ptTrigBinMC]->Fill(coneptsumTrack, GetEventWeight());
          fhPtTrigBinSumPtClusterConeMC[ptTrigBinMC]->Fill(coneptsumClust, GetEventWeight());
        }
      } // MC
      
    } // proper pT bin found 
  } // pT trigger bins
  
  //
  // Shower shape dependent histograms
  //
  if(fFillSSHisto)
  {
    //fhPtLambda0Eiso->Fill(pt, m02, coneptsum);
    
    fhELambda0 [isolated]->Fill(energy, m02, GetEventWeight());
    fhPtLambda0[isolated]->Fill(pt,     m02, GetEventWeight());
    //fhELambda1 [isolated]->Fill(energy, m20);
    
    //
    // MC
    //
    if(IsDataMC())
    {
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
        fhPtLambda0MC[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
      
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
      {
        if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MC[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
        else if( mcIndex == kmcEtaDecay ) fhPtLambda0MC[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
      }
      
      fhPtLambda0MC[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
      
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
      {
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtLambda0MCConv[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCConv[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
          else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCConv[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
        }
        
        fhPtLambda0MCConv[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
      } // Conversion
      
      //
      // Check overlaps
      //
      Int_t noverlaps = 0;
      if ( fFillOverlapHistograms && fCluster )
      {
        const UInt_t nlabels = fCluster->GetNLabels();
        Int_t overpdg[nlabels];
        Int_t overlab[nlabels];
        noverlaps = GetMCAnalysisUtils()->GetNOverlaps(fCluster->GetLabels(), nlabels, mcTag, -1,
                                                       GetMC(), overpdg, overlab);
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtNOverlap[kmcPhoton][isolated]->Fill(pt, noverlaps, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtNOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight());
          else if( mcIndex == kmcEtaDecay ) fhPtNOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight());
        }
        
        fhPtNOverlap[mcIndex][isolated]->Fill(pt, noverlaps, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
        {
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtNOverlapConv[kmcPhoton][isolated]->Fill(pt, noverlaps, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhPtNOverlapConv[kmcPi0DecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight());
            else if( mcIndex == kmcEtaDecay ) fhPtNOverlapConv[kmcEtaDecayLostPair][isolated]->Fill(pt, noverlaps, GetEventWeight());
          }
          
          fhPtNOverlapConv[mcIndex][isolated]->Fill(pt, noverlaps, GetEventWeight());
        } // Conversion
        
        
        if ( noverlaps == 1 )
        {
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtLambda0MCWith1Overlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCWith1Overlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
            else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCWith1Overlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
          }
          
          fhPtLambda0MCWith1Overlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
          {
            if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtLambda0MCConvWith1Overlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
            
            if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
            {
              if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCConvWith1Overlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
              else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCConvWith1Overlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
            }
            
            fhPtLambda0MCConvWith1Overlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
          } // Conversion
        } // At least 1 overlap
        else if (noverlaps == 0 ) // No overlap
        {
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtLambda0MCWithNoOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCWithNoOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
            else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCWithNoOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
          }
          
          fhPtLambda0MCWithNoOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
          {
            if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              fhPtLambda0MCConvWithNoOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
            
            if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
            {
              if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCConvWithNoOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
              else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCConvWithNoOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
            }
            
            fhPtLambda0MCConvWithNoOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
          } // Conversion
        } // more than 1 overlap
      }
      
      //
      // Sum in cone and shower shape in different pT bins of candidate and overlap condition
      //
      if(fFillPtTrigBinHistograms && fFillSSHisto)
      {
        Int_t ptTrigBinMC        = ptTrigBin+mcIndex  *fNPtTrigBin;
        Int_t ptTrigBinMCPhoton  = ptTrigBin+kmcPhoton*fNPtTrigBin;
        Int_t ptTrigBinMCPi0Lost = ptTrigBin+kmcPi0DecayLostPair*fNPtTrigBin;
        Int_t ptTrigBinMCEtaLost = ptTrigBin+kmcEtaDecayLostPair*fNPtTrigBin;
        
        if ( ptTrigBin >= 0 )
        {
          fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , m02, GetEventWeight());
          fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, m02, GetEventWeight());

          if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
          {
            fhPtTrigBinLambda0vsSumPtTrackConeMC  [ptTrigBinMC]->Fill(coneptsumTrack, m02, GetEventWeight());
            fhPtTrigBinLambda0vsSumPtClusterConeMC[ptTrigBinMC]->Fill(coneptsumClust, m02, GetEventWeight());
          }
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          {
            fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMCPhoton]->Fill(coneptsum , m02, GetEventWeight());
            fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMCPhoton]->Fill(coneleadpt, m02, GetEventWeight());

            if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinLambda0vsSumPtTrackConeMC  [ptTrigBinMCPhoton]->Fill(coneptsumTrack, m02, GetEventWeight());
              fhPtTrigBinLambda0vsSumPtClusterConeMC[ptTrigBinMCPhoton]->Fill(coneptsumClust, m02, GetEventWeight());
            }
          }
          
          if( mcIndex == kmcPi0Decay )
          {
            fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMCPi0Lost]->Fill(coneptsum , m02, GetEventWeight());
            fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMCPi0Lost]->Fill(coneleadpt, m02, GetEventWeight());

            if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinLambda0vsSumPtTrackConeMC  [ptTrigBinMCPi0Lost]->Fill(coneptsumTrack, m02, GetEventWeight());
              fhPtTrigBinLambda0vsSumPtClusterConeMC[ptTrigBinMCPi0Lost]->Fill(coneptsumClust, m02, GetEventWeight());
            }
          }
          
          if( mcIndex == kmcEtaDecay )
          {
            fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMCEtaLost]->Fill(coneptsum , m02, GetEventWeight());
            fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMCEtaLost]->Fill(coneleadpt, m02, GetEventWeight());

            if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinLambda0vsSumPtTrackConeMC  [ptTrigBinMCEtaLost]->Fill(coneptsumTrack, m02, GetEventWeight());
              fhPtTrigBinLambda0vsSumPtClusterConeMC[ptTrigBinMCEtaLost]->Fill(coneptsumClust, m02, GetEventWeight());
            }
          }
          
          if(fFillOverlapHistograms)
          {
            if ( noverlaps == 0 )
            {
              fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[ptTrigBinMC]->Fill(coneptsum, m02, GetEventWeight());
              
              if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
              {
                fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap  [ptTrigBinMC]->Fill(coneptsumTrack, m02, GetEventWeight());
                fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[ptTrigBinMC]->Fill(coneptsumClust, m02, GetEventWeight());
              }
              
              if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              {
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[ptTrigBinMCPhoton]->Fill(coneptsum, m02, GetEventWeight());
                
                if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                {
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap  [ptTrigBinMCPhoton]->Fill(coneptsumTrack, m02, GetEventWeight());
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[ptTrigBinMCPhoton]->Fill(coneptsumClust, m02, GetEventWeight());
                }
              }
              
              if( mcIndex == kmcPi0Decay )
              {
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[ptTrigBinMCPi0Lost]->Fill(coneptsum, m02, GetEventWeight());
                
                if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                {
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap  [ptTrigBinMCPi0Lost]->Fill(coneptsumTrack, m02, GetEventWeight());
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[ptTrigBinMCPi0Lost]->Fill(coneptsumClust, m02, GetEventWeight());
                }
              }
              
              if( mcIndex == kmcEtaDecay )
              {
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[ptTrigBinMCEtaLost]->Fill(coneptsum, m02, GetEventWeight());
                
                if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                {
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap  [ptTrigBinMCEtaLost]->Fill(coneptsumTrack, m02, GetEventWeight());
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[ptTrigBinMCEtaLost]->Fill(coneptsumClust, m02, GetEventWeight());
                }
              }
              
              
            } // nover = 0
            else if ( noverlaps == 1 )
            {
              fhPtTrigBinLambda0vsSumPtConeMC1Overlap[ptTrigBinMC]->Fill(coneptsum, m02, GetEventWeight());
              if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
              {
                fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap  [ptTrigBinMC]->Fill(coneptsumTrack, m02, GetEventWeight());
                fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[ptTrigBinMC]->Fill(coneptsumClust, m02, GetEventWeight());
              }
              
              if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
              {
                fhPtTrigBinLambda0vsSumPtConeMC1Overlap[ptTrigBinMCPhoton]->Fill(coneptsum, m02, GetEventWeight());
                
                if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                {
                  fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap  [ptTrigBinMCPhoton]->Fill(coneptsumTrack, m02, GetEventWeight());
                  fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[ptTrigBinMCPhoton]->Fill(coneptsumClust, m02, GetEventWeight());
                }
                
                if( mcIndex == kmcPi0Decay )
                {
                  fhPtTrigBinLambda0vsSumPtConeMC1Overlap[ptTrigBinMCPi0Lost]->Fill(coneptsum, m02, GetEventWeight());
                  
                  if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                  {
                    fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap  [ptTrigBinMCPi0Lost]->Fill(coneptsumTrack, m02, GetEventWeight());
                    fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[ptTrigBinMCPi0Lost]->Fill(coneptsumClust, m02, GetEventWeight());
                  }
                }
                
                if( mcIndex == kmcEtaDecay )
                {
                  fhPtTrigBinLambda0vsSumPtConeMC1Overlap[ptTrigBinMCEtaLost]->Fill(coneptsum, m02, GetEventWeight());
                  
                  if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                  {
                    fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap  [ptTrigBinMCEtaLost]->Fill(coneptsumTrack, m02, GetEventWeight());
                    fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[ptTrigBinMCEtaLost]->Fill(coneptsumClust, m02, GetEventWeight());
                  }
                }
              }
            } // nover = 1
          } // fill overlaps
          
        } // pt bin exists
      } // fFillPtTrigBinHistograms && fFillSSHisto
      
    } // MC
    
    if(GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0 &&
       GetModuleNumber(pCandidate) >= GetFirstSMCoveredByTRD()  )
    {
      fhELambda0TRD [isolated]->Fill(energy, m02, GetEventWeight());
      fhPtLambda0TRD[isolated]->Fill(pt    , m02, GetEventWeight());
      //fhELambda1TRD [isolated]->Fill(energy, m20 );
    }
    
    if(fFillNLMHistograms)
    {
      fhNLocMax[isolated]->Fill(energy, nMaxima, GetEventWeight());
      if     (nMaxima==1)
      {
        fhELambda0LocMax1[isolated]->Fill(energy, m02, GetEventWeight());
        fhELambda1LocMax1[isolated]->Fill(energy, m02, GetEventWeight());
      }
      else if(nMaxima==2)
      {
        fhELambda0LocMax2[isolated]->Fill(energy, m02, GetEventWeight());
        fhELambda1LocMax2[isolated]->Fill(energy, m02, GetEventWeight());
      }
      else
      {
        fhELambda0LocMaxN[isolated]->Fill(energy, m02, GetEventWeight());
        fhELambda1LocMaxN[isolated]->Fill(energy, m02, GetEventWeight());
      }
    }
  } // SS histo fill
  
  // Track matching dependent histograms
  if ( fFillTMHisto && fCluster )
  {
    Float_t dZ  = fCluster->GetTrackDz();
    Float_t dR  = fCluster->GetTrackDx();
    
//    if(fCluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
//    {
//      dR = 2000., dZ = 2000.;
//      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(fCluster->GetID(),dZ,dR);
//    }
    
    //printf("ParticleIsolation: dPhi %f, dEta %f\n",dR,dZ);
    if(fhTrackMatchedDEta[isolated] && TMath::Abs(dR) < 999)
    {
      fhTrackMatchedDEta[isolated]->Fill(energy, dZ, GetEventWeight());
      fhTrackMatchedDPhi[isolated]->Fill(energy, dR, GetEventWeight());
      if(energy > 0.5) fhTrackMatchedDEtaDPhi[isolated]->Fill(dZ, dR, GetEventWeight());
      if(IsDataMC())
      {
        fhTrackMatchedDEtaMC[mcIndex][isolated]->Fill(energy, dZ, GetEventWeight());
        fhTrackMatchedDPhiMC[mcIndex][isolated]->Fill(energy, dR, GetEventWeight());
        if(energy > 0.5) fhTrackMatchedDEtaDPhiMC[mcIndex][isolated]->Fill(dZ, dR, GetEventWeight());
      }
    }
    
    // Check dEdx and E/p of matched clusters
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {
      AliVTrack *track = GetCaloUtils()->GetMatchedTrack(fCluster, GetReader()->GetInputEvent());
      
      if(track)
      {
        Float_t dEdx = track->GetTPCsignal();
        fhdEdx[isolated]->Fill(fCluster->E(), dEdx, GetEventWeight());
        
        Float_t eOverp = fCluster->E()/track->P();
        fhEOverP[isolated]->Fill(fCluster->E(),  eOverp, GetEventWeight());
      }
      //else
      //  printf("AliAnaParticleIsolation::FillTrackMatchingShowerShapeHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
      
      
      if(IsDataMC())
      {
        if ( !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion)  )
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                     GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 2.5, GetEventWeight());
          else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 0.5, GetEventWeight());
          else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 1.5, GetEventWeight());
          else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 3.5, GetEventWeight());
          
        }
        else
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)      ||
                     GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 6.5, GetEventWeight());
          else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 4.5, GetEventWeight());
          else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[isolated]->Fill(energy, 5.5, GetEventWeight());
          else                                                                                   fhTrackMatchedMCParticle[isolated]->Fill(energy, 7.5, GetEventWeight());
        }
      }  // MC
    } // match window
  }// TM histos fill
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
  snprintf(onePar, buffersize,"Isolation Cand. Detector: %s;",fIsoDetectorString.Data()) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"Fill histo UE %d",fFillUEBandSubtractHistograms) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fReMakeIC =%d (Flag for reisolation during histogram filling);",fReMakeIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fMakeSeveralIC=%d (Flag for isolation with several cuts at the same time );",fMakeSeveralIC) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillTMHisto=%d (Flag for track matching histograms);",fFillTMHisto) ;
  parList+=onePar ;
  snprintf(onePar, buffersize,"fFillSSHisto=%d (Flag for shower shape histograms);",fFillSSHisto) ;
  parList+=onePar ;
  
  if(fMakeSeveralIC)
  {
    snprintf(onePar, buffersize,"fNCones =%d (Number of cone sizes);",fNCones) ;
    parList+=onePar ;
    snprintf(onePar, buffersize,"fNPtThresFrac=%d (Flag for isolation with several cuts at the same time);",fNPtThresFrac) ;
    parList+=onePar ;
    
    for(Int_t icone = 0; icone < fNCones ; icone++)
    {
      snprintf(onePar, buffersize,"fConeSizes[%d]=%1.2f (isolation cone size);",icone, fConeSizes[icone]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtThresholds[%d]=%1.2f (isolation pt threshold);",ipt, fPtThresholds[ipt]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fPtFractions[%d]=%1.2f (isolation pt fraction threshold);",ipt, fPtFractions[ipt]) ;
      parList+=onePar ;
    }
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      snprintf(onePar, buffersize,"fSumPtThresholds[%d]=%1.2f (isolation sum pt threshold);",ipt, fSumPtThresholds[ipt]) ;
      parList+=onePar ;
    }
  }
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in IC class.
  if(!fMakeSeveralIC)parList += GetIsolationCut()->GetICParametersList() ;
  
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
  Int_t   nPoverEbins = GetHistogramRanges()->GetHistoPOverEBins();
  Float_t pOverEmax   = GetHistogramRanges()->GetHistoPOverEMax();
  Float_t pOverEmin   = GetHistogramRanges()->GetHistoPOverEMin();
  
  Int_t   nptsumbins    = GetHistogramRanges()->GetHistoNPtSumBins();
  Float_t ptsummax      = GetHistogramRanges()->GetHistoPtSumMax();
  Float_t ptsummin      = GetHistogramRanges()->GetHistoPtSumMin();
  Int_t   nptinconebins = GetHistogramRanges()->GetHistoNPtInConeBins();
  Float_t ptinconemax   = GetHistogramRanges()->GetHistoPtInConeMax();
  Float_t ptinconemin   = GetHistogramRanges()->GetHistoPtInConeMin();
  
  Int_t   nmultbin = GetHistogramRanges()->GetHistoTrackMultiplicityBins();
  Int_t   multmax  = GetHistogramRanges()->GetHistoTrackMultiplicityMax ();
  Int_t   multmin  = GetHistogramRanges()->GetHistoTrackMultiplicityMin ();
  
  // n cell bins for TH3
  Int_t cellBins  = 15;
  Float_t cellMax = 15;
  Float_t cellMin = 0;
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
  
  //Float_t ptthre    = GetIsolationCut()->GetPtThreshold();
  //Float_t ptsumthre = GetIsolationCut()->GetSumPtThreshold();
  //Float_t ptfrac    = GetIsolationCut()->GetPtFraction();
  Float_t r         = GetIsolationCut()->GetConeSize();
  Int_t   method    = GetIsolationCut()->GetICMethod() ;
  Int_t   particle  = GetIsolationCut()->GetParticleTypeInCone() ;
  
  TString sThreshold = "";
  if      ( method == AliIsolationCut::kSumPtIC )
  {
    sThreshold = Form(", %2.2f < #Sigma #it{p}_{T}^{in cone} < %2.2f GeV/#it{c}",
                      GetIsolationCut()->GetSumPtThreshold(), GetIsolationCut()->GetSumPtThresholdMax());
    if(GetIsolationCut()->GetSumPtThresholdMax() > 200)
      sThreshold = Form(", #Sigma #it{p}_{T}^{in cone} = %2.2f GeV/#it{c}",
                        GetIsolationCut()->GetSumPtThreshold());
  }
  else if ( method == AliIsolationCut::kPtThresIC)
  {
    sThreshold = Form(", %2.2f < #it{p}_{T}^{th} < %2.2f GeV/#it{c}",
                      GetIsolationCut()->GetPtThreshold(),GetIsolationCut()->GetPtThresholdMax());
    if(GetIsolationCut()->GetSumPtThreshold() > 200)
      sThreshold = Form(", #it{p}_{T}^{th} = %2.2f GeV/#it{c}",
                        GetIsolationCut()->GetPtThreshold());
  }
  else if ( method == AliIsolationCut::kPtFracIC)
    sThreshold = Form(", #Sigma #it{p}_{T}^{in cone}/#it{p}_{T}^{trig} = %2.2f" ,
                      GetIsolationCut()->GetPtFraction());
  
  TString sParticle = ", x^{0,#pm}";
  if      ( particle == AliIsolationCut::kOnlyNeutral )  sParticle = ", x^{0}";
  else if ( particle == AliIsolationCut::kOnlyCharged )  sParticle = ", x^{#pm}";
  
  TString parTitle    = Form("#it{R} = %2.2f%s%s"     ,GetIsolationCut()->GetConeSize(),sThreshold.Data(),sParticle.Data());
  TString parTitleR   = Form("#it{R} = %2.2f%s"       ,GetIsolationCut()->GetConeSize(),sParticle.Data());
  TString parTitleRCh = Form("#it{R} = %2.2f, x^{#pm}",GetIsolationCut()->GetConeSize());
  TString parTitleRNe = Form("#it{R} = %2.2f, x^{0}"  ,GetIsolationCut()->GetConeSize());
  
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
  
  // Not Isolated histograms, reference histograms
  
  fhENoIso  = new TH1F("hENoIso",
                       Form("Number of not isolated leading particles vs #it{p}_{T}, %s",parTitle.Data()),
                       nptbins,ptmin,ptmax);
  fhENoIso->SetYTitle("#it{counts}");
  fhENoIso->SetXTitle("E (GeV/#it{c})");
  outputContainer->Add(fhENoIso) ;
  
  fhPtNoIso  = new TH1F("hPtNoIso",
                        Form("Number of not isolated leading particles vs #it{p}_{T}, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
  fhPtNoIso->SetYTitle("#it{counts}");
  fhPtNoIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNoIso) ;
  
  if(fStudyExoticTrigger)
  {
    fhENoIsoExoTrigger   = new TH1F
    ("hENoIsoExoTrigger",
     Form("Number of not isolated particles vs E, %s, exot>0.97",parTitle.Data()),
     nptbins,ptmin,ptmax);
    fhENoIsoExoTrigger->SetYTitle("d#it{N} / d#it{E}");
    fhENoIsoExoTrigger->SetXTitle("#it{E} (GeV/#it{c})");
    outputContainer->Add(fhENoIsoExoTrigger) ;
    
    fhPtNoIsoExoTrigger  = new TH1F
    ("hPtNoIsoExoTrigger",
     Form("Number of not isolated particles vs #it{p}_{T}, %s, exot>0.97",parTitle.Data()),
     nptbins,ptmin,ptmax);
    fhPtNoIsoExoTrigger->SetYTitle("d#it{N} / #it{p}_{T}");
    fhPtNoIsoExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNoIsoExoTrigger) ;
  }
  
  fhEtaPhiNoIso  = new TH2F("hEtaPhiNoIso",
                            Form("Number of not isolated leading particles #eta vs #varphi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiNoIso->SetXTitle("#eta");
  fhEtaPhiNoIso->SetYTitle("#varphi (rad)");
  outputContainer->Add(fhEtaPhiNoIso) ;
  
  if(IsDataMC())
  {
    // For histograms in arrays, index in the array, corresponding to any particle origin
    
    for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
    {
      fhPtNoIsoMC[imc]  = new TH1F(Form("hPtNoIsoMC%s",mcPartName[imc].Data()),
                                   Form("#it{p}_{T} of NOT isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax);
      fhPtNoIsoMC[imc]->SetYTitle("#it{counts}");
      fhPtNoIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtNoIsoMC[imc]) ;
      
      fhPtIsoMC[imc]  = new TH1F(Form("hPtMC%s",mcPartName[imc].Data()),
                                 Form("#it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                 nptbins,ptmin,ptmax);
      fhPtIsoMC[imc]->SetYTitle("#it{counts}");
      fhPtIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtIsoMC[imc]) ;
      
      fhPhiIsoMC[imc]  = new TH2F(Form("hPhiMC%s",mcPartName[imc].Data()),
                                  Form("#varphi vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                  nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhPhiIsoMC[imc]->SetYTitle("#varphi (rad)");
      fhPhiIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPhiIsoMC[imc]) ;
      
      fhEtaIsoMC[imc]  = new TH2F(Form("hEtaMC%s",mcPartName[imc].Data()),
                                  Form("#eta vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                  nptbins,ptmin,ptmax,netabins,etamin,etamax);
      fhEtaIsoMC[imc]->SetYTitle("#eta");
      fhEtaIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhEtaIsoMC[imc]) ;
    }
  }
  
  // Histograms for tagged candidates as decay
  if(fFillTaggedDecayHistograms)
  {
    TString isoName [] = {"NoIso","Iso"};
    TString isoTitle[] = {"Not isolated"  ,"isolated"};
    
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      for(Int_t iso = 0; iso < 2; iso++)
      {
        if(fMakeSeveralIC && iso) continue;
        fhPtDecay[iso][ibit]  =
        new TH1F(Form("hPtDecay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
                 Form("Number of %s leading pi0 decay particles vs #it{p}_{T}, bit %d, %s",isoTitle[iso].Data(),fDecayBits[ibit],parTitle.Data()),
                 nptbins,ptmin,ptmax);
        fhPtDecay[iso][ibit]->SetYTitle("#it{counts}");
        fhPtDecay[iso][ibit]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtDecay[iso][ibit]) ;
        
        fhEtaPhiDecay[iso][ibit]  =
        new TH2F(Form("hEtaPhiDecay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
                 Form("Number of %s leading Pi0 decay particles #eta vs #varphi, bit %d, %s",isoTitle[iso].Data(),fDecayBits[ibit],parTitle.Data()),
                 netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiDecay[iso][ibit]->SetXTitle("#eta");
        fhEtaPhiDecay[iso][ibit]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiDecay[iso][ibit]) ;
        
        if(fFillSSHisto)
        {
          fhPtLambda0Decay[iso][ibit]  = new TH2F
          (Form("hPtLambda0Decay%s_bit%d",isoName[iso].Data(),fDecayBits[ibit]),
           Form("%s cluster : #it{p}_{T} vs #lambda_{0}, decay bit %d, %s",isoTitle[iso].Data(), fDecayBits[ibit], parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0Decay[iso][ibit]->SetYTitle("#lambda_{0}^{2}");
          fhPtLambda0Decay[iso][ibit]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0Decay[iso][ibit]) ;
        }
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
          {
            fhPtDecayMC[iso][ibit][imc]  =
            new TH1F(Form("hPtDecay%s_bit%d_MC%s",isoName[iso].Data(),fDecayBits[ibit],mcPartName[imc].Data()),
                     Form("#it{p}_{T} of %s, decay bit %d,  %s, %s",isoTitle[iso].Data(),fDecayBits[ibit],mcPartType[imc].Data(),parTitle.Data()),
                     nptbins,ptmin,ptmax);
            fhPtDecayMC[iso][ibit][imc]->SetYTitle("#it{counts}");
            fhPtDecayMC[iso][ibit][imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtDecayMC[iso][ibit][imc]) ;
          }// MC particle loop
        }// MC
      } // bit loop
    } //iso loop
  }// decay
  
  if(!fMakeSeveralIC)
  {
    TString isoName [] = {"NoIso","Iso"};
    TString isoTitle[] = {"Not isolated"  ,"isolated"};
    
    fhEIso   = new TH1F("hE",
                        Form("Number of isolated particles vs E, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
    fhEIso->SetYTitle("d#it{N} / d#it{E}");
    fhEIso->SetXTitle("#it{E} (GeV/#it{c})");
    outputContainer->Add(fhEIso) ;
    
    fhPtIso  = new TH1F("hPt",
                        Form("Number of isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                        nptbins,ptmin,ptmax);
    fhPtIso->SetYTitle("d#it{N} / #it{p}_{T}");
    fhPtIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtIso) ;
    
    if(fStudyExoticTrigger)
    {
      fhEIsoExoTrigger   = new TH1F
      ("hEIsoExoTrigger",
       Form("Number of isolated particles vs E, %s, exot>0.97",parTitle.Data()),
       nptbins,ptmin,ptmax);
      fhEIsoExoTrigger->SetYTitle("d#it{N} / d#it{E}");
      fhEIsoExoTrigger->SetXTitle("#it{E} (GeV/#it{c})");
      outputContainer->Add(fhEIsoExoTrigger) ;
      
      fhPtIsoExoTrigger  = new TH1F
      ("hPtIsoExoTrigger",
       Form("Number of isolated particles vs #it{p}_{T}, %s, exot>0.97",parTitle.Data()),
       nptbins,ptmin,ptmax);
      fhPtIsoExoTrigger->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtIsoExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtIsoExoTrigger) ;
    }
    
    fhPhiIso  = new TH2F("hPhi",
                         Form("Number of isolated particles vs #varphi, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhPhiIso->SetYTitle("#varphi (rad)");
    fhPhiIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPhiIso) ;
    
    fhEtaIso  = new TH2F("hEta",
                         Form("Number of isolated particles vs #eta, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhEtaIso) ;
    
    fhEtaPhiIso  = new TH2F("hEtaPhiIso",
                            Form("Number of isolated particles #eta vs #varphi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiIso->SetXTitle("#eta");
    fhEtaPhiIso->SetYTitle("#varphi (rad)");
    outputContainer->Add(fhEtaPhiIso) ;
    
    if(IsHighMultiplicityAnalysisOn())
    {
      fhPtCentralityIso  = new TH2F("hPtCentrality",
                                    Form("centrality vs #it{p}_{T} for isolated particles, %s",parTitle.Data()),
                                    nptbins,ptmin,ptmax, 100,0,100);
      fhPtCentralityIso->SetYTitle("centrality");
      fhPtCentralityIso->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPtCentralityIso) ;
      
      fhPtEventPlaneIso  = new TH2F("hPtEventPlane",
                                    Form("event plane angle vs #it{p}_{T} for isolated particles, %s",parTitle.Data()),
                                    nptbins,ptmin,ptmax, 100,0,TMath::Pi());
      fhPtEventPlaneIso->SetYTitle("Event plane angle (rad)");
      fhPtEventPlaneIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtEventPlaneIso) ;
    }
    
    if(fFillNLMHistograms)
    {
      fhPtNLocMaxIso  = new TH2F("hPtNLocMax",
                                 Form("Number of isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                                 nptbins,ptmin,ptmax,10,0,10);
      fhPtNLocMaxIso->SetYTitle("#it{NLM}");
      fhPtNLocMaxIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      
      fhPtNLocMaxNoIso  = new TH2F("hPtNLocMaxNoIso",
                                   Form("Number of not isolated particles vs #it{p}_{T}, %s",parTitle.Data()),
                                   nptbins,ptmin,ptmax,10,0,10);
      fhPtNLocMaxNoIso->SetYTitle("#it{NLM}");
      fhPtNLocMaxNoIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtNLocMaxNoIso) ;
    }

    fhConePtLead  = new TH2F("hConePtLead",
                            Form("Track or Cluster  leading #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                            nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhConePtLead->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
    fhConePtLead->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConePtLead) ;
    
    fhConeSumPt  = new TH2F
    ("hConePtSum",
     Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",r),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPt->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPt) ;
    
    if(fStudyExoticTrigger)
    {
      fhConeSumPtExoTrigger  = new TH2F
      ("hConePtSumExoTrigger",
       Form("#Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f, exo trigger",r),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtExoTrigger) ;
    }
    
//    fhPtLambda0Eiso = new TH3F
//    (Form("hPtLambda0Eiso"),
//     Form("ABCD 3D Matrix: #it{p}_{T} vs #lambda_{0}^{2} vs E_{T}^{iso}, %s",parTitle.Data()),70,0.,70.,200,0.,2.,110,-10.,100.);
//    fhPtLambda0Eiso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
//    fhPtLambda0Eiso->SetYTitle("#lambda_{0}^{2}");
//    fhPtLambda0Eiso->SetZTitle("E_{T}^{iso} (GeV/#it{c})");
//    outputContainer->Add(fhPtLambda0Eiso) ;
    
    fhConeSumPtTrigEtaPhi  = new TH2F("hConePtSumTrigEtaPhi",
                                      Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} in isolation cone for %s",parTitleR.Data()),
                                      netabins,etamin,etamax,nphibins,phimin,phimax);
    fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
    fhConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
    outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
    
    fhPtInCone  = new TH2F("hPtInCone",
                           Form("#it{p}_{T} of clusters and tracks in isolation cone for %s",parTitleR.Data()),
                           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInCone) ;
    
    if(fStudyExoticTrigger)
    {
      fhPtInConeExoTrigger  = new TH2F("hPtInConeExoTrigger",
                                       Form("#it{p}_{T}  in isolation cone for %s, exotic trigger",parTitleR.Data()),
                                       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtInConeExoTrigger) ;
    }
    
    if(fFillBackgroundBinHistograms)
    {
      fhPtLeadConeBin              = new TH1F*[fNBkgBin];
      fhSumPtConeBin               = new TH1F*[fNBkgBin];
      if(fFillUEBandSubtractHistograms) fhSumPtConeAfterEtaBandUESubBin = new TH1F*[fNBkgBin];
      
      if(fFillSSHisto)
      {
        fhPtLeadConeBinLambda0     = new TH2F*[fNBkgBin];
        fhSumPtConeBinLambda0      = new TH2F*[fNBkgBin];
        if(fFillUEBandSubtractHistograms) fhSumPtConeAfterEtaBandUESubBinLambda0 = new TH2F*[fNBkgBin];
      }
      
      if(fFillTaggedDecayHistograms)
      {
        fhPtLeadConeBinDecay       = new TH1F*[fNBkgBin*fNDecayBits];
        fhSumPtConeBinDecay        = new TH1F*[fNBkgBin*fNDecayBits];
      }
      
      if(IsDataMC())
      {
        fhPtLeadConeBinMC          = new TH1F*[fNBkgBin*fgkNmcTypes];
        fhSumPtConeBinMC           = new TH1F*[fNBkgBin*fgkNmcTypes];
        if(fFillUEBandSubtractHistograms) fhSumPtConeAfterEtaBandUESubBinMC = new TH1F*[fNBkgBin*fgkNmcTypes];
        
        if(fFillSSHisto)
        {
          fhPtLeadConeBinLambda0MC = new TH2F*[fNBkgBin*fgkNmcTypes];
          fhSumPtConeBinLambda0MC  = new TH2F*[fNBkgBin*fgkNmcTypes];
          if(fFillUEBandSubtractHistograms) fhSumPtConeAfterEtaBandUESubBinLambda0MC = new TH2F*[fNBkgBin*fgkNmcTypes];
        }
      }
      
      for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
      {        
        fhPtLeadConeBin[ibin]  = new TH1F
        (Form("hPtLeadCone_Bin%d",ibin),
         Form("cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
              fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
        fhPtLeadConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtLeadConeBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLeadConeBin[ibin]) ;
        
        fhSumPtConeBin[ibin]  = new TH1F
        (Form("hSumPtCone_Bin%d",ibin),
         Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
              fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
        fhSumPtConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhSumPtConeBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhSumPtConeBin[ibin]) ;
        
        if(fFillUEBandSubtractHistograms)
        {
          fhSumPtConeAfterEtaBandUESubBin[ibin]  = new TH1F
          (Form("hSumPtConeAfterEtaBandUESub_Bin%d",ibin),
           Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
          fhSumPtConeAfterEtaBandUESubBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhSumPtConeAfterEtaBandUESubBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSumPtConeAfterEtaBandUESubBin[ibin]) ;
        }
        
        if(fFillTaggedDecayHistograms)
        {
          for(Int_t idecay = 0; idecay < fNDecayBits; idecay++)
          {
            Int_t bindecay = ibin+idecay*fNBkgBin;

            fhPtLeadConeBinDecay[bindecay]  = new TH1F
            (Form("hPtLeadCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
            fhPtLeadConeBinDecay[bindecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtLeadConeBinDecay[bindecay]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtLeadConeBinDecay[bindecay]) ;
            
            fhSumPtConeBinDecay[bindecay]  = new TH1F
            (Form("hSumPtCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c},  %s",
                  fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
            fhSumPtConeBinDecay[bindecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhSumPtConeBinDecay[bindecay]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhSumPtConeBinDecay[bindecay]) ;
          }
        }

        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
          {
            Int_t binmc = ibin+imc*fNBkgBin;
            fhPtLeadConeBinMC[binmc]  = new TH1F
            (Form("hPtLeadCone_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
             Form("in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, MC %s, %s",
                  fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax);
            fhPtLeadConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtLeadConeBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtLeadConeBinMC[binmc]) ;
            
            fhSumPtConeBinMC[binmc]  = new TH1F
            (Form("hSumPtCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
             Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                  fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax);
            fhSumPtConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhSumPtConeBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhSumPtConeBinMC[binmc]) ;
            
            if(fFillUEBandSubtractHistograms)
            {
              fhSumPtConeAfterEtaBandUESubBinMC[binmc]  = new TH1F
              (Form("hSumPtConeAfterEtaBandUESub_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                    fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax);
              fhSumPtConeAfterEtaBandUESubBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
              fhSumPtConeAfterEtaBandUESubBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              outputContainer->Add(fhSumPtConeAfterEtaBandUESubBinMC[binmc]) ;
            }
          } // MC particle loop
        }
        
        if(fFillSSHisto)
        {
          fhPtLeadConeBinLambda0[ibin]  = new TH2F
          (Form("hPtLeadConeLambda0_Bin%d",ibin),
           Form("#lambda_{0}, in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLeadConeBinLambda0[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtLeadConeBinLambda0[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLeadConeBinLambda0[ibin]) ;
          
          fhSumPtConeBinLambda0[ibin]  = new TH2F
          (Form("hSumPtConeLambda0_Bin%d",ibin),
           Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhSumPtConeBinLambda0[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhSumPtConeBinLambda0[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSumPtConeBinLambda0[ibin]) ;
          
          if(fFillUEBandSubtractHistograms)
          {
            fhSumPtConeAfterEtaBandUESubBinLambda0[ibin]  = new TH2F
            (Form("hSumPtConeAfterEtaBandUESubLambda0_Bin%d",ibin),
             Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
                  fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhSumPtConeAfterEtaBandUESubBinLambda0[ibin]->SetYTitle("#lambda_{0}^{2}");
            fhSumPtConeAfterEtaBandUESubBinLambda0[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhSumPtConeAfterEtaBandUESubBinLambda0[ibin]) ;
          }
          
          if(IsDataMC())
          {
            for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
            {
              Int_t binmc = ibin+imc*fNBkgBin;
              fhPtLeadConeBinLambda0MC[binmc]  = new TH2F
              (Form("hPtLeadConeLambda0_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
               Form("#lambda_{0}, in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, MC %s, %s",
                    fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLeadConeBinLambda0MC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtLeadConeBinLambda0MC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              outputContainer->Add(fhPtLeadConeBinLambda0MC[binmc]) ;
              
              fhSumPtConeBinLambda0MC[binmc]  = new TH2F
              (Form("hSumPtConeLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                    fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhSumPtConeBinLambda0MC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhSumPtConeBinLambda0MC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              outputContainer->Add(fhSumPtConeBinLambda0MC[binmc]) ;
              
              if(fFillUEBandSubtractHistograms)
              {
                fhSumPtConeAfterEtaBandUESubBinLambda0MC[binmc]  = new TH2F
                (Form("hSumPtConeAfterEtaBandUESubLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
                 Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                      fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
                fhSumPtConeAfterEtaBandUESubBinLambda0MC[binmc]->SetYTitle("#lambda_{0}^{2}");
                fhSumPtConeAfterEtaBandUESubBinLambda0MC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
                outputContainer->Add(fhSumPtConeAfterEtaBandUESubBinLambda0MC[binmc]) ;
              }
            } // MC particle loop
          }
        } // shower shape on

      } // pt in cone bin loop
    } // bkg cone pt bin histograms

    if(fFillPtTrigBinHistograms)
    {
      fhPtTrigBinPtLeadCone = new TH1F*[fNPtTrigBin];
      fhPtTrigBinSumPtCone  = new TH1F*[fNPtTrigBin];
  
      fhPtTrigBinPtLeadConeDecay = new TH1F*[fNPtTrigBin*fNDecayBits];
      fhPtTrigBinSumPtConeDecay  = new TH1F*[fNPtTrigBin*fNDecayBits];
     
      if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
      {
        fhPtTrigBinSumPtTrackCone   = new TH1F*[fNPtTrigBin];
        fhPtTrigBinSumPtClusterCone = new TH1F*[fNPtTrigBin];
        
        fhPtTrigBinSumPtTrackConeDecay   = new TH1F*[fNPtTrigBin*fNDecayBits];
        fhPtTrigBinSumPtClusterConeDecay = new TH1F*[fNPtTrigBin*fNDecayBits];
      } 
      
      if(IsDataMC())
      {
        fhPtTrigBinPtLeadConeMC = new TH1F*[fNPtTrigBin*fgkNmcTypes];
        fhPtTrigBinSumPtConeMC  = new TH1F*[fNPtTrigBin*fgkNmcTypes];
        if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
        {
          fhPtTrigBinSumPtTrackConeMC   = new TH1F*[fNPtTrigBin*fgkNmcTypes];
          fhPtTrigBinSumPtClusterConeMC = new TH1F*[fNPtTrigBin*fgkNmcTypes];
        }
      }
      
      if(fFillSSHisto)
      {
        fhPtTrigBinLambda0vsPtLeadCone = new TH2F*[fNPtTrigBin];
        fhPtTrigBinLambda0vsSumPtCone  = new TH2F*[fNPtTrigBin];
        if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
        {
          fhPtTrigBinLambda0vsSumPtTrackCone   = new TH2F*[fNPtTrigBin];
          fhPtTrigBinLambda0vsSumPtClusterCone = new TH2F*[fNPtTrigBin];
        }
        
        if(IsDataMC())
        {
          fhPtTrigBinLambda0vsPtLeadConeMC = new TH2F*[fNPtTrigBin*fgkNmcTypes];
          fhPtTrigBinLambda0vsSumPtConeMC  = new TH2F*[fNPtTrigBin*fgkNmcTypes];
          if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
          {
            fhPtTrigBinLambda0vsSumPtTrackConeMC   = new TH2F*[fNPtTrigBin*fgkNmcTypes];
            fhPtTrigBinLambda0vsSumPtClusterConeMC = new TH2F*[fNPtTrigBin*fgkNmcTypes];
          }
          
          if(fFillOverlapHistograms)
          {
            fhPtTrigBinLambda0vsSumPtConeMCNoOverlap  = new TH2F*[fNPtTrigBin*fgkNmcTypes];
            fhPtTrigBinLambda0vsSumPtConeMC1Overlap   = new TH2F*[fNPtTrigBin*fgkNmcTypes];
            if ( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap   = new TH2F*[fNPtTrigBin*fgkNmcTypes];
              fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap = new TH2F*[fNPtTrigBin*fgkNmcTypes];
              fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap    = new TH2F*[fNPtTrigBin*fgkNmcTypes];
              fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap  = new TH2F*[fNPtTrigBin*fgkNmcTypes];
            }
          }
        }
      }
      
      for(Int_t ibin = 0; ibin < fNPtTrigBin; ibin++)
      {
        fhPtTrigBinPtLeadCone[ibin]  = new TH1F
        (Form("hPtTrigBin_PtLeadCone_Bin%d",ibin),
         Form("#it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
              fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
        fhPtTrigBinPtLeadCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtTrigBinPtLeadCone[ibin]->SetXTitle("#it{p}_{T}^{in cone} (GeV/#it{c})");
        outputContainer->Add(fhPtTrigBinPtLeadCone[ibin]) ;
        
        fhPtTrigBinSumPtCone[ibin]  = new TH1F
        (Form("hPtTrigBin_SumPtCone_Bin%d",ibin),
         Form("#Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
              fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptsumbins,ptsummin,ptsummax);
        fhPtTrigBinSumPtCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtTrigBinSumPtCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
        outputContainer->Add(fhPtTrigBinSumPtCone[ibin]) ;
        
        if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
        {
          fhPtTrigBinSumPtTrackCone[ibin]  = new TH1F
          (Form("hPtTrigBin_SumPtTrackCone_Bin%d",ibin),
           Form("#Sigma #it{p}_{T}^{in cone}_{track} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRCh.Data()),nptsumbins,ptsummin,ptsummax);
          fhPtTrigBinSumPtTrackCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhPtTrigBinSumPtTrackCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinSumPtTrackCone[ibin]) ;
          
          fhPtTrigBinSumPtClusterCone[ibin]  = new TH1F
          (Form("hPtTrigBin_SumPtClusterCone_Bin%d",ibin),
           Form("#Sigma #it{p}_{T}^{in cone}_{cluster} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRNe.Data()),nptsumbins,ptsummin,ptsummax);
          fhPtTrigBinSumPtClusterCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
          fhPtTrigBinSumPtClusterCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinSumPtClusterCone[ibin]) ;
        }
        
        if(fFillTaggedDecayHistograms)
        {
          for(Int_t idecay = 0; idecay < fNDecayBits; idecay++)
          {
            Int_t binDecay = ibin+idecay*fNPtTrigBin;
            
            fhPtTrigBinPtLeadConeDecay[binDecay]  = new TH1F
            (Form("hPtTrigBin_PtLeadCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax);
            fhPtTrigBinPtLeadConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinPtLeadConeDecay[binDecay]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinPtLeadConeDecay[binDecay]) ;
            
            fhPtTrigBinSumPtConeDecay[binDecay]  = new TH1F
            (Form("hPtTrigBin_SumPtCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, #Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptsumbins,ptsummin,ptsummax);
            fhPtTrigBinSumPtConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinSumPtConeDecay[binDecay]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinSumPtConeDecay[binDecay]) ;
            
            if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinSumPtTrackConeDecay[binDecay]  = new TH1F
              (Form("hPtTrigBin_SumPtTrackCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
               Form("Decay bit %d, #Sigma #it{p}_{T}^{in cone}_{track} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                    fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRCh.Data()),nptsumbins,ptsummin,ptsummax);
              fhPtTrigBinSumPtTrackConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
              fhPtTrigBinSumPtTrackConeDecay[binDecay]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinSumPtTrackConeDecay[binDecay]) ;
              
              fhPtTrigBinSumPtClusterConeDecay[binDecay]  = new TH1F
              (Form("hPtTrigBin_SumPtClusterCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
               Form("Decay bit %d, #Sigma #it{p}_{T}^{in cone}_{cluster} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                    fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRNe.Data()),nptsumbins,ptsummin,ptsummax);
              fhPtTrigBinSumPtClusterConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
              fhPtTrigBinSumPtClusterConeDecay[binDecay]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinSumPtClusterConeDecay[binDecay]) ;
            }
          }
        }
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
          {
            Int_t binmc = ibin+imc*fNPtTrigBin;
            fhPtTrigBinPtLeadConeMC[binmc]  = new TH1F
            (Form("hPtTrigBin_PtLeadCone_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
             Form("#it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, MC %s, %s",
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptbins,ptmin,ptmax);
            fhPtTrigBinPtLeadConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinPtLeadConeMC[binmc]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinPtLeadConeMC[binmc]) ;
            
            fhPtTrigBinSumPtConeMC[binmc]  = new TH1F
            (Form("hPtTrigBin_SumPtCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
             Form("#Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),nptsumbins,ptsummin,ptsummax);
            fhPtTrigBinSumPtConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinSumPtConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinSumPtConeMC[binmc]) ;
            
            if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
            {
              fhPtTrigBinSumPtTrackConeMC[binmc]  = new TH1F
              (Form("hPtTrigBin_SumPtTrackCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#Sigma #it{p}_{T}^{in cone}_{track}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRCh.Data()),nptsumbins,ptsummin,ptsummax);
              fhPtTrigBinSumPtTrackConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
              fhPtTrigBinSumPtTrackConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinSumPtTrackConeMC[binmc]) ;
              
              fhPtTrigBinSumPtClusterConeMC[binmc]  = new TH1F
              (Form("hPtTrigBin_SumPtClusterCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#Sigma #it{p}_{T}^{in cone}_{cluster}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRNe.Data()),nptsumbins,ptsummin,ptsummax);
              fhPtTrigBinSumPtClusterConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
              fhPtTrigBinSumPtClusterConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinSumPtClusterConeMC[binmc]) ;
            }
          } // MC particle loop
        } // MC
        
        if(fFillSSHisto)
        {
          fhPtTrigBinLambda0vsPtLeadCone[ibin]  = new TH2F
          (Form("hPtTrigBin_PtLeadConeVSLambda0_Bin%d",ibin),
           Form("#lambda_{0} vs #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtTrigBinLambda0vsPtLeadCone[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtTrigBinLambda0vsPtLeadCone[ibin]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinLambda0vsPtLeadCone[ibin]) ;
          
          fhPtTrigBinLambda0vsSumPtCone[ibin]  = new TH2F
          (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d",ibin),
           Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleR.Data()),nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
          fhPtTrigBinLambda0vsSumPtCone[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtTrigBinLambda0vsSumPtCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinLambda0vsSumPtCone[ibin]) ;
          
          if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
          {
            fhPtTrigBinLambda0vsSumPtTrackCone[ibin]  = new TH2F
            (Form("hPtTrigBin_SumPtTrackConeVSLambda0_Bin%d",ibin),
             Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{track} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRCh.Data()),nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
            fhPtTrigBinLambda0vsSumPtTrackCone[ibin]->SetYTitle("#lambda_{0}^{2}");
            fhPtTrigBinLambda0vsSumPtTrackCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinLambda0vsSumPtTrackCone[ibin]) ;         
            
            fhPtTrigBinLambda0vsSumPtClusterCone[ibin]  = new TH2F
            (Form("hPtTrigBin_SumPtClusterConeVSLambda0_Bin%d",ibin),
             Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{cluster} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitleRNe.Data()),nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
            fhPtTrigBinLambda0vsSumPtClusterCone[ibin]->SetYTitle("#lambda_{0}^{2}");
            fhPtTrigBinLambda0vsSumPtClusterCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinLambda0vsSumPtClusterCone[ibin]) ;
          }
          
          if(IsDataMC())
          {
            for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
            {              
              Int_t binmc = ibin+imc*fNPtTrigBin;
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]  = new TH2F
              (Form("hPtTrigBin_PtLeadConeVSLambda0_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
               Form("#lambda_{0} vs #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),
               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinLambda0vsPtLeadConeMC[binmc]) ;
              
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]  = new TH2F
              (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),
               nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinLambda0vsSumPtConeMC[binmc]) ;
              
              if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
              {
                fhPtTrigBinLambda0vsSumPtTrackConeMC[binmc]  = new TH2F
                (Form("hPtTrigBin_SumPtTrackConeVSLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
                 Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{track}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                      fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRCh.Data()),
                 nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                fhPtTrigBinLambda0vsSumPtTrackConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
                fhPtTrigBinLambda0vsSumPtTrackConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
                outputContainer->Add(fhPtTrigBinLambda0vsSumPtTrackConeMC[binmc]) ;
                
                fhPtTrigBinLambda0vsSumPtClusterConeMC[binmc]  = new TH2F
                (Form("hPtTrigBin_SumPtClusterConeVSLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
                 Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{cluster}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                      fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRNe.Data()),
                 nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                fhPtTrigBinLambda0vsSumPtClusterConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
                fhPtTrigBinLambda0vsSumPtClusterConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
                outputContainer->Add(fhPtTrigBinLambda0vsSumPtClusterConeMC[binmc]) ;
              }
              
              if(fFillOverlapHistograms)
              {  
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[binmc]  = new TH2F
                (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d_MC_NoOverlap%s",ibin,mcPartName[imc].Data()),
                 Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, No Overlaps %s, %s",
                      fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),
                 nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
                outputContainer->Add(fhPtTrigBinLambda0vsSumPtConeMCNoOverlap[binmc]) ;

                fhPtTrigBinLambda0vsSumPtConeMC1Overlap[binmc]  = new TH2F
                (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d_MC_1Overlap%s",ibin,mcPartName[imc].Data()),
                 Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, 1 Overlap %s, %s",
                      fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleR.Data()),
                 nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                fhPtTrigBinLambda0vsSumPtConeMC1Overlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                fhPtTrigBinLambda0vsSumPtConeMC1Overlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
                outputContainer->Add(fhPtTrigBinLambda0vsSumPtConeMC1Overlap[binmc]) ;

                if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
                {
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap[binmc]  = new TH2F
                  (Form("hPtTrigBin_SumPtTrackConeVSLambda0_Bin%d_MC_NoOverlap%s",ibin,mcPartName[imc].Data()),
                   Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{track}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, No Overlaps %s, %s",
                        fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRCh.Data()),
                   nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
                  outputContainer->Add(fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap[binmc]) ;
                  
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[binmc]  = new TH2F
                  (Form("hPtTrigBin_SumPtClusterConeVSLambda0_Bin%d_MC_NoOverlap%s",ibin,mcPartName[imc].Data()),
                   Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{cluster}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, No Overlaps %s, %s",
                        fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRNe.Data()),
                   nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
                  outputContainer->Add(fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap[binmc]) ;                
                  
                  fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap[binmc]  = new TH2F
                  (Form("hPtTrigBin_SumPtTrackConeVSLambda0_Bin%d_MC_1Overlap%s",ibin,mcPartName[imc].Data()),
                   Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{track}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, 1 Overlap %s, %s",
                        fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRCh.Data()),
                   nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                  fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                  fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{track} (GeV/#it{c})");
                  outputContainer->Add(fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap[binmc]) ;
                  
                  fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[binmc]  = new TH2F
                  (Form("hPtTrigBin_SumPtClusterConeVSLambda0_Bin%d_MC_1Overlap%s",ibin,mcPartName[imc].Data()),
                   Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}_{cluster}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC, 1 Overlap %s, %s",
                        fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitleRNe.Data()),
                   nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
                  fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[binmc]->SetYTitle("#lambda_{0}^{2}");
                  fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone}_{cluster} (GeV/#it{c})");
                  outputContainer->Add(fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap[binmc]) ;

                }               
              } // Overlap histograms
            } // MC particle loop
          } // MC
        } // SS histo
      } // pt trig bin loop
    } // pt trig bin histograms
    
    if(IsHighMultiplicityAnalysisOn())
    {
      fhPtInConeCent  = new TH2F("hPtInConeCent",
                                 Form("#it{p}_{T} in isolation cone for %s",parTitleR.Data()),
                                 100,0,100,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConeCent->SetXTitle("centrality");
      outputContainer->Add(fhPtInConeCent) ;
    }
    
    // Cluster only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyCharged)
    {
      fhConeSumPtCluster  = new TH2F
      ("hConePtSumCluster",
       Form("Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",r),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtCluster) ;
      
      if(fStudyExoticTrigger)
      {
        fhConeSumPtClusterExoTrigger  = new TH2F
        ("hConePtSumClusterExoTrigger",
         Form("Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f, exo trigger",r),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterExoTrigger) ;
      }
      
      if(fStudyPtCutInCone)
      {
        fhConeNClusterPerMinPtCut = new TH2F
        ("hConeNClusterPerMinPtCut",
         Form("N clusters, different min #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin,multmin,multmax);
        fhConeNClusterPerMinPtCut->SetYTitle("#it{N}^{cluster}");
        fhConeNClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNClusterPerMinPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterPerMinPtCut) ;
        
        fhConeNClusterPerMinPtCutLargePtTrig = new TH2F
        ("hConeNClusterPerMinPtCutLargePtTrig",
         Form("N cluster, different min #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nmultbin,multmin,multmax);
        fhConeNClusterPerMinPtCutLargePtTrig->SetYTitle("#it{N}^{cluster}");
        fhConeNClusterPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNClusterPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNClusterPerMinPtCutLargePtTrig) ;

        fhConeSumPtClusterPerMinPtCut = new TH2F
        ("hConePtSumClusterPerMinPtCut",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMinPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMinPtCut) ;
        
        fhConeSumPtClusterPerMinPtCutLargePtTrig = new TH2F
        ("hConePtSumClusterPerMinPtCutLargePtTrig",
         Form("Cluster #Sigma #it{p}_{T}, different min #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMinPtCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMinPtCutLargePtTrig) ;

        
        fhConeSumPtClusterPerMaxPtCut = new TH2F
        ("hConePtSumClusterPerMaxPtCut",
         Form("Cluster #Sigma #it{p}_{T}, different max #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMaxPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMaxPtCut->SetXTitle("#it{p}_{T, max} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMaxPtCut->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMaxPtCut) ;
        
        fhConeSumPtClusterPerMaxPtCutLargePtTrig = new TH2F
        ("hConePtSumClusterPerMaxPtCutLargePtTrig",
         Form("Cluster #Sigma #it{p}_{T}, different max #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,-0.5,fNPtCutsInCone-0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerMaxPtCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerMaxPtCutLargePtTrig->SetXTitle("#it{p}_{T, max} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtClusterPerMaxPtCutLargePtTrig->GetXaxis()->SetBinLabel(i ,Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerMaxPtCutLargePtTrig) ;
      }
      
      if(fStudyRCutInCone)
      {
        fhConeSumPtClusterPerRCut = new TH2F
        ("hConePtSumClusterPerRCut","Cluster #Sigma #it{p}_{T}, different #it{R} cuts",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerRCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerRCut->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhConeSumPtClusterPerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerRCut) ;
        
        fhConeSumPtClusterPerRCutLargePtTrig = new TH2F
        ("hConePtSumClusterPerRCutLargePtTrig","Cluster #Sigma #it{p}_{T}, different #it{R} cuts, #it{p}_{T}^{trig} > 10 GeV",
        fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerRCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerRCutLargePtTrig->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhConeSumPtClusterPerRCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerRCutLargePtTrig) ;
        
        fhPtClusterInConePerRCut = new TH2F
        ("hPtClusterInConePerRCut","Cluster #it{p}_{T}, different #it{R} cuts",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerRCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerRCut->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhPtClusterInConePerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhPtClusterInConePerRCut) ;
        
        fhPtClusterInConePerRCutLargePtTrig = new TH2F
        ("hPtClusterInConePerRCutLargePtTrig","Cluster #it{p}_{T}, different #it{R} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerRCutLargePtTrig->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerRCutLargePtTrig->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhPtClusterInConePerRCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhPtClusterInConePerRCutLargePtTrig) ;
      }
      
      if(fStudyNCellsCut)
      {
        fhConeSumPtClusterPerNCellCut = new TH2F
        ("hConePtSumClusterPerNCellCut","Cluster #Sigma #it{p}_{T}, different #it{N}_{cell} cuts",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerNCellCut->SetYTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtClusterPerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhConeSumPtClusterPerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerNCellCut) ;
        
        fhConeSumPtClusterPerNCellCutLargePtTrig = new TH2F
        ("hConePtSumClusterPerNCellCutLargePtTrig","Cluster #Sigma #it{p}_{T}, different #it{N}_{cell} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerNCellCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtClusterPerNCellCutLargePtTrig->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhConeSumPtClusterPerNCellCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerNCellCutLargePtTrig) ;
        
        fhPtClusterInConePerNCellCut = new TH2F
        ("hPtClusterInConePerNCellCut","Cluster #it{p}_{T}, different #it{N}_{cell} cuts",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerNCellCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhPtClusterInConePerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhPtClusterInConePerNCellCut) ;
        
        fhPtClusterInConePerNCellCutLargePtTrig = new TH2F
        ("hPtClusterInConePerNCellCutLargePtTrig","Cluster #it{p}_{T}, different #it{N}_{cell} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerNCellCutLargePtTrig->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerNCellCutLargePtTrig->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhPtClusterInConePerNCellCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhPtClusterInConePerNCellCutLargePtTrig) ;

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
      
      if(fStudyExoticTrigger)
      {  
        fhConeSumPtClusterPerExoCut = new TH2F
        ("hConePtSumClusterPerExoCut","Cluster #Sigma #it{p}_{T}, different exoticity cuts",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerExoCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerExoCut->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhConeSumPtClusterPerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerExoCut) ;
        
        fhConeSumPtClusterPerExoCutLargePtTrig = new TH2F
        ("hConePtSumClusterPerExoCutLargePtTrig","Cluster #Sigma #it{p}_{T}, different exoticity cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClusterPerExoCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterPerExoCutLargePtTrig->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhConeSumPtClusterPerExoCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtClusterPerExoCutLargePtTrig) ;
        
        fhPtClusterInConePerExoCut = new TH2F
        ("hPtClusterInConePerExoCut","Cluster #it{p}_{T}, different exoticity cuts",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerExoCut->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerExoCut->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhPtClusterInConePerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhPtClusterInConePerExoCut) ;
        
        fhPtClusterInConePerExoCutLargePtTrig = new TH2F
        ("hPtClusterInConePerExoCutLargePtTrig","Cluster #it{p}_{T}, different exoticity cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtClusterInConePerExoCutLargePtTrig->SetYTitle("#it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPtClusterInConePerExoCutLargePtTrig->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhPtClusterInConePerExoCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhPtClusterInConePerExoCutLargePtTrig) ;
      }
      
      fhConePtLeadCluster  = new TH2F("hConeLeadPtCluster",
                                    Form("Cluster leading in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhConePtLeadCluster->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLeadCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadCluster) ;

      
      if(fFillCellHistograms)
      {
        fhConeSumPtCell  = new TH2F("hConePtSumCell",
                                    Form("Cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUECluster  = new TH2F("hConePtSumEtaBandUECluster",
                                                "#Sigma cluster #it{p}_{T} in UE Eta Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUECluster) ;
        
        fhConeSumPtPhiBandUECluster  = new TH2F("hConePtSumPhiBandUECluster",
                                                "#Sigma cluster #it{p}_{T} UE Phi Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtPhiBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPhiBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUECluster) ;
       
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhConeSumPtEtaBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumEtaBandUEClusterTrigEtaPhi",
                                                            "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE Eta Band",
                                                            netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhi) ;
          
          fhConeSumPtPhiBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumPhiBandUEClusterTrigEtaPhi",
                                                            "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} UE Phi Band",
                                                            netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhi) ;
          
          fhEtaBandClusterEtaPhi  = new TH2F("hEtaBandClusterEtaPhi",
                                             Form("#eta vs #varphi of clusters in #eta band isolation cone for #it{R} =  %2.2f",r),
                                             netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaBandClusterEtaPhi->SetXTitle("#eta");
          fhEtaBandClusterEtaPhi->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaBandClusterEtaPhi) ;
          
          fhPhiBandClusterEtaPhi  = new TH2F("hPhiBandClusterEtaPhi",
                                             Form("#eta vs #varphi of clusters in #varphi band isolation cone for #it{R} =  %2.2f",r),
                                             netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhPhiBandClusterEtaPhi->SetXTitle("#eta");
          fhPhiBandClusterEtaPhi->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhPhiBandClusterEtaPhi) ;
          
          fhEtaPhiInConeCluster= new TH2F("hEtaPhiInConeCluster",
                                          Form("#eta vs #varphi of clusters in cone for #it{R} =  %2.2f",r),
                                          netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaPhiInConeCluster->SetXTitle("#eta");
          fhEtaPhiInConeCluster->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiInConeCluster) ;
          
          fhEtaPhiCluster= new TH2F("hEtaPhiCluster",
                                    Form("#eta vs #varphi of all clusters"),
                                    netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaPhiCluster->SetXTitle("#eta");
          fhEtaPhiCluster->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiCluster) ;
        }
          
        if(fFillCellHistograms)
        {
          fhConeSumPtEtaBandUECell  = new TH2F("hConePtSumEtaBandUECell",
                                               "#Sigma cell #it{p}_{T} in UE Eta Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtEtaBandUECell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaBandUECell) ;
          
          fhConeSumPtPhiBandUECell  = new TH2F("hConePtSumPhiBandUECell",
                                               "#Sigma cell #it{p}_{T} UE Phi Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPhiBandUECell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiBandUECell) ;
          
          fhConeSumPtEtaBandUECellTrigEtaPhi  = new TH2F("hConePtSumEtaBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #varphi, #Sigma cell #it{p}_{T} in UE Eta Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaBandUECellTrigEtaPhi) ;
          
          fhConeSumPtPhiBandUECellTrigEtaPhi  = new TH2F("hConePtSumPhiBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #varphi, #Sigma cell #it{p}_{T} UE Phi Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiBandUECellTrigEtaPhi) ;
        }
        
        fhEtaBandClusterPt  = new TH2F("hEtaBandClusterPt",
                                     Form("#it{p}_{T} of clusters in #eta band isolation cone for #it{R} =  %2.2f",r),
                                     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhEtaBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhEtaBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandClusterPt) ;
        
        fhPhiBandClusterPt  = new TH2F("hPhiBandClusterPt",
                                     Form("#it{p}_{T} of clusters in #varphi band isolation cone for #it{R} =  %2.2f",r),
                                     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhPhiBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhPhiBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandClusterPt) ;        
      }
      
      fhPtClusterInCone  = new TH2F("hPtClusterInCone",
                                    Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtClusterInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtClusterInCone) ;

      if(fStudyExoticTrigger)
      {
        fhPtClusterInConeExoTrigger  = new TH2F("hPtClusterInConeExoTrigger",
                                      Form("#it{p}_{T} of clusters in isolation cone for #it{R} = %2.2f, exotic trigger",r),
                                      nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtClusterInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtClusterInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtClusterInConeExoTrigger) ;
      }
      
      if(fFillCellHistograms)
      {
        fhPtCellInCone  = new TH2F("hPtCellInCone",
                                   Form("#it{p}_{T} of cells in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,1000,0,50);
        fhPtCellInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtCellInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtCellInCone) ;
        
        fhEtaBandCell  = new TH2F("hEtaBandCell",
                                  Form("#col vs #row of cells in #eta band isolation cone for #it{R} =  %2.2f",r),
                                  96,0,95,128,0,127);
        fhEtaBandCell->SetXTitle("#col");
        fhEtaBandCell->SetYTitle("#row");
        outputContainer->Add(fhEtaBandCell) ;
        
        fhPhiBandCell  = new TH2F("hPhiBandCell",
                                  Form("#col vs #row of cells in #varphi band isolation cone for #it{R} =  %2.2f",r),
                                  96,0,95,128,0,127);
        fhPhiBandCell->SetXTitle("#col");
        fhPhiBandCell->SetYTitle("#row");
        outputContainer->Add(fhPhiBandCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUENormCluster  = new TH2F("hConeSumPtEtaUENormCluster",
                                              Form("Clusters #Sigma #it{p}_{T} in normalized #eta band, #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUENormCluster->SetYTitle("#Sigma #it{p}_{T}^{#eta-band}_{norm} (GeV/#it{c})");
        fhConeSumPtEtaUENormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUENormCluster) ;
        
        fhConeSumPtPhiUENormCluster  = new TH2F("hConeSumPtPhiUENormCluster",
                                              Form("Clusters #Sigma #it{p}_{T} in normalized #varphi band, #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUENormCluster->SetYTitle("#Sigma #it{p}_{T}^{#varphi-band}_{norm} (GeV/#it{c})");
        fhConeSumPtPhiUENormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUENormCluster) ;
        
        fhConeSumPtEtaUESubCluster  = new TH2F("hConeSumPtEtaUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUESubCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubCluster) ;
        
        fhConeSumPtPhiUESubCluster  = new TH2F("hConeSumPtPhiUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUESubCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPhiUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubCluster) ;
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhConeSumPtEtaUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubClusterTrigEtaPhi",
                                                           Form("Trigger #eta vs #varphi, Clusters #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                           netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubClusterTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubClusterTrigEtaPhi",
                                                           Form("Trigger #eta vs #varphi, Clusters #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                           netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubClusterTrigEtaPhi) ;
        }
        
        if(fFillCellHistograms)
        {
          fhConeSumPtEtaUESubCell  = new TH2F("hConeSumPtEtaUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtEtaUESubCell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubCell) ;
          
          fhConeSumPtPhiUESubCell  = new TH2F("hConeSumPtPhiUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtPhiUESubCell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubCell) ;
          
          fhConeSumPtEtaUESubCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #varphi, Cells #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #varphi, Cells #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubCellTrigEtaPhi) ;
        }
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhFractionClusterOutConeEta  = new TH2F("hFractionClusterOutConeEta",
                                                  Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance",r),
                                                  nptbins,ptmin,ptmax,100,0,1);
          fhFractionClusterOutConeEta->SetYTitle("#it{fraction}");
          fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionClusterOutConeEta) ;
          
          fhFractionClusterOutConeEtaTrigEtaPhi  = new TH2F("hFractionClusterOutConeEtaTrigEtaPhi",
                                                            Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance, in trigger #eta-#varphi ",r),
                                                            netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;
          
          fhFractionClusterOutConePhi  = new TH2F("hFractionClusterOutConePhi",
                                                  Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #varphi acceptance",r),
                                                  nptbins,ptmin,ptmax,100,0,1);
          fhFractionClusterOutConePhi->SetYTitle("#it{fraction}");
          fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionClusterOutConePhi) ;
          
          fhFractionClusterOutConePhiTrigEtaPhi  = new TH2F("hFractionClusterOutConePhiTrigEtaPhi",
                                                            Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #varphi acceptance, in trigger #eta-#varphi ",r),
                                                            netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;
          
          
          fhConeSumPtSubvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCluster",
                                                            Form("#Sigma #it{p}_{T} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCluster);
          
          fhConeSumPtSubNormvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCluster",
                                                                Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                                nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCluster);
          
          fhConeSumPtSubvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCluster",
                                                            Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCluster);
          
          fhConeSumPtSubNormvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCluster",
                                                                Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                                nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCluster);
          
          fhConeSumPtVSUEClusterEtaBand  = new TH2F("hConeSumPtVSUEClusterEtaBand",
                                                    Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #eta band for cluster (before normalization), R=%2.2f",r),
                                                    nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
          fhConeSumPtVSUEClusterEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUEClusterEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtVSUEClusterEtaBand);
          
          fhConeSumPtVSUEClusterPhiBand  = new TH2F("hConeSumPtVSUEClusterPhiBand",
                                                    Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #varphi band for cluster (before normalization), R=%2.2f",r),
                                                    nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
          fhConeSumPtVSUEClusterPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUEClusterPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtVSUEClusterPhiBand);
        }
        
        if(fFillCellHistograms)
        {
          fhFractionCellOutConeEta  = new TH2F("hFractionCellOutConeEta",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConeEta->SetYTitle("#it{fraction}");
          fhFractionCellOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConeEta) ;
          
          fhFractionCellOutConeEtaTrigEtaPhi  = new TH2F("hFractionCellOutConeEtaTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance, in trigger #eta-#varphi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConeEtaTrigEtaPhi) ;
          
          fhFractionCellOutConePhi  = new TH2F("hFractionCellOutConePhi",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #varphi acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConePhi->SetYTitle("#it{fraction}");
          fhFractionCellOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConePhi) ;
          
          fhFractionCellOutConePhiTrigEtaPhi  = new TH2F("hFractionCellOutConePhiTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #varphi acceptance, in trigger #eta-#varphi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConePhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConePhiTrigEtaPhi) ;
          
          
          fhConeSumPtSubvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubNormvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCell);
          
          fhConeSumPtSubNormvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCell);
        }
      }
    }
    
    // Track only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    {
      fhConePtLeadTrack  = new TH2F
      ("hConeLeadPtTrack",
       Form("Track leading in isolation cone for #it{R} = %2.2f",r),
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhConePtLeadTrack->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLeadTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadTrack) ;
      
      fhPtTrackInCone  = new TH2F
      ("hPtTrackInCone",
       Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f",r),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInCone) ;

      if(IsDataMC())
      {
        TString mcChPartName[] = {"Pion","Kaon","Proton","Other"};
        for(Int_t imc = 0; imc < 4; imc++)
        {
          fhPtTrackInConeMCPrimary[imc]  = new TH2F
          (Form("hPtTrackInCone_Primary_%s",mcChPartName[imc].Data()),
           Form("reconstructed #it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, primary MC %s",r,mcChPartName[imc].Data()),
           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
          fhPtTrackInConeMCPrimary[imc]->SetYTitle("#it{p}_{T in cone}^{reco} (GeV/#it{c})");
          fhPtTrackInConeMCPrimary[imc]->SetXTitle("#it{p}_{T}^{reco} (GeV/#it{c})");
          outputContainer->Add(fhPtTrackInConeMCPrimary[imc]) ;
          
          fhPtTrackInConeMCSecondary[imc]  = new TH2F
          (Form("hPtTrackInCone_Secondary_%s",mcChPartName[imc].Data()),
           Form("reconstructed #it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, secondary MC %s",r,mcChPartName[imc].Data()),
           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
          fhPtTrackInConeMCSecondary[imc]->SetYTitle("#it{p}_{T in cone}^{reco} (GeV/#it{c})");
          fhPtTrackInConeMCSecondary[imc]->SetXTitle("#it{p}_{T}^{reco} (GeV/#it{c})");
          outputContainer->Add(fhPtTrackInConeMCSecondary[imc]) ;
          
          fhPtTrackInConeMCPrimaryGener[imc]  = new TH2F
          (Form("hPtTrackInCone_Gener_Primary_%s",mcChPartName[imc].Data()),
           Form("generated #it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, primary MC %s",r,mcChPartName[imc].Data()),
           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
          fhPtTrackInConeMCPrimaryGener[imc]->SetYTitle("#it{p}_{T in cone}^{gener} (GeV/#it{c})");
          fhPtTrackInConeMCPrimaryGener[imc]->SetXTitle("#it{p}_{T}^{gener} (GeV/#it{c})");
          outputContainer->Add(fhPtTrackInConeMCPrimaryGener[imc]) ;
          
          fhPtTrackInConeMCSecondaryGener[imc]  = new TH2F
          (Form("hPtTrackInCone_Gener_Secondary_%s",mcChPartName[imc].Data()),
           Form("generated #it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, secondary MC %s",r,mcChPartName[imc].Data()),
           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
          fhPtTrackInConeMCSecondaryGener[imc]->SetYTitle("#it{p}_{T in cone}^{gener} (GeV/#it{c})");
          fhPtTrackInConeMCSecondaryGener[imc]->SetXTitle("#it{p}_{T}^{gener} (GeV/#it{c})");
          outputContainer->Add(fhPtTrackInConeMCSecondaryGener[imc]) ;
        }
        
      }
      
      if(fStudyExoticTrigger)
      {
        fhPtTrackInConeExoTrigger  = new TH2F("hPtTrackInConeExoTrigger",
                                                Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f, exotic trigger",r),
                                                nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeExoTrigger->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtTrackInConeExoTrigger->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeExoTrigger) ;
      }
      
      if(fStudyPtCutInCone)
      {
        fhConeNTrackPerMinPtCut = new TH2F
        ("hConeNTrackPerMinPtCut",
         Form("N tracks, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
        fhConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackPerMinPtCut) ;
        
        fhConeNTrackPerMinPtCutLargePtTrig = new TH2F
        ("hConeNTrackPerMinPtCutLargePtTrig",
         Form("N tracks, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
        fhConeNTrackPerMinPtCutLargePtTrig->SetYTitle("#it{N}^{track}");
        fhConeNTrackPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeNTrackPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeNTrackPerMinPtCutLargePtTrig) ;

        fhConeSumPtTrackPerMinPtCut = new TH2F
        ("hConePtSumTrackPerMinPtCut",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMinPtCut) ;
        
        fhConeSumPtTrackPerMinPtCutLargePtTrig = new TH2F
        ("hConePtSumTrackPerMinPtCutLargePtTrig",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMinPtCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMinPtCutLargePtTrig) ;

        if(fFillUEBandSubtractHistograms)
        {
          fhPerpConeNTrackPerMinPtCut = new TH2F
          ("hPerpConeNTrackPerMinPtCut",
           Form("N tracks, different #it{p}_{T} cuts in perpendicular cone for #it{R} = %2.2f",r),
           fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
          fhPerpConeNTrackPerMinPtCut->SetYTitle("#it{N}^{track}");
          fhPerpConeNTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
          for(Int_t i = 1; i <= fNPtCutsInCone; i++)
            fhPerpConeNTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
          outputContainer->Add(fhPerpConeNTrackPerMinPtCut) ;
          
          fhPerpConeNTrackPerMinPtCutLargePtTrig = new TH2F
          ("hPerpConeNTrackPerMinPtCutLargePtTrig",
           Form("N tracks, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
           fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nmultbin,multmin,multmax);
          fhPerpConeNTrackPerMinPtCutLargePtTrig->SetYTitle("#it{N}^{track}");
          fhPerpConeNTrackPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
          for(Int_t i = 1; i <= fNPtCutsInCone; i++)
            fhPerpConeNTrackPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
          outputContainer->Add(fhPerpConeNTrackPerMinPtCutLargePtTrig) ;
          
          fhPerpConeSumPtTrackPerMinPtCut = new TH2F
          ("hPerpConePtSumTrackPerMinPtCut",
           Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in perpendicular cone for #it{R} = %2.2f",r),
           fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
          fhPerpConeSumPtTrackPerMinPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPerpConeSumPtTrackPerMinPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
          for(Int_t i = 1; i <= fNPtCutsInCone; i++)
            fhPerpConeSumPtTrackPerMinPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
          outputContainer->Add(fhPerpConeSumPtTrackPerMinPtCut) ;
          
          fhPerpConeSumPtTrackPerMinPtCutLargePtTrig = new TH2F
          ("hPerpConePtSumTrackPerMinPtCutLargePtTrig",
           Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in perpendicular cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
           fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
          fhPerpConeSumPtTrackPerMinPtCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPerpConeSumPtTrackPerMinPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
          for(Int_t i = 1; i <= fNPtCutsInCone; i++)
            fhPerpConeSumPtTrackPerMinPtCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMinPtCutInCone[i-1]));
          outputContainer->Add(fhPerpConeSumPtTrackPerMinPtCutLargePtTrig) ;
        }
        
        fhConeSumPtTrackPerMaxPtCut = new TH2F
        ("hConePtSumTrackPerMaxPtCut",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMaxPtCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMaxPtCut->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMaxPtCut->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMaxPtCut) ;
        
        fhConeSumPtTrackPerMaxPtCutLargePtTrig = new TH2F
        ("hConePtSumTrackPerMaxPtCutLargePtTrig",
         Form("Track #Sigma #it{p}_{T}, different #it{p}_{T} cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
         fNPtCutsInCone,0.5,fNPtCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerMaxPtCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerMaxPtCutLargePtTrig->SetXTitle("#it{p}_{T, min} (GeV/#it{c})");
        for(Int_t i = 1; i <= fNPtCutsInCone; i++)
          fhConeSumPtTrackPerMaxPtCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.1f",fMaxPtCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerMaxPtCutLargePtTrig) ;
      }
      
      if(fStudyEtaCutInCone)
      {
        fhConeSumPtTrackPerEtaCut = new TH2F("hConePtSumTrackPerEtaCut",
                                             Form("Track #Sigma #it{p}_{T}, different #eta cuts in isolation cone for #it{R} = %2.2f",r),
                                             fNEtaCutsInCone,0.5,fNEtaCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerEtaCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerEtaCut->SetXTitle("#eta_{max}");
          for(Int_t i = 1; i <= fNEtaCutsInCone; i++)
        fhConeSumPtTrackPerEtaCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fEtaCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerEtaCut) ;
        
        fhConeSumPtTrackPerEtaCutLargePtTrig = new TH2F("hConePtSumTrackPerEtaCutLargePtTrig",
                                                        Form("Track #Sigma #it{p}_{T}, different #eta cuts in isolation cone for #it{R} = %2.2f, #it{p}_{T}^{trig} > 10 GeV",r),
                                                        fNEtaCutsInCone,0.5,fNEtaCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerEtaCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerEtaCutLargePtTrig->SetXTitle("#eta_{max}");
        for(Int_t i = 1; i <= fNEtaCutsInCone; i++)
          fhConeSumPtTrackPerEtaCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fEtaCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerEtaCutLargePtTrig) ;
      }
      
      if(fStudyRCutInCone)
      {
        fhConeSumPtTrackPerRCut = new TH2F
        ("hConePtSumTrackPerRCut","Track #Sigma #it{p}_{T}, different #it{R} cuts",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerRCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerRCut->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhConeSumPtTrackPerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerRCut) ;
        
        fhConeSumPtTrackPerRCutLargePtTrig = new TH2F
        ("hConePtSumTrackPerRCutLargePtTrig","Track #Sigma #it{p}_{T}, different #it{R} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerRCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerRCutLargePtTrig->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhConeSumPtTrackPerRCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerRCutLargePtTrig) ;
        
        fhPtTrackInConePerRCut = new TH2F
        ("hPtTrackInConePerRCut","Track #it{p}_{T}, different #it{R} cuts",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerRCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerRCut->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhPtTrackInConePerRCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhPtTrackInConePerRCut) ;
        
        fhPtTrackInConePerRCutLargePtTrig = new TH2F
        ("hPtTrackInConePerRCutLargePtTrig","Track #it{p}_{T}, different #it{R} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNRCutsInCone,0.5,fNRCutsInCone+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerRCutLargePtTrig->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerRCutLargePtTrig->SetXTitle("#it{R}");
        for(Int_t i = 1; i <= fNRCutsInCone; i++)
          fhPtTrackInConePerRCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fRCutInCone[i-1]));
        outputContainer->Add(fhPtTrackInConePerRCutLargePtTrig) ;
      }
      
      if(fStudyNCellsCut)
      {
        fhConeSumPtTrackPerNCellCut = new TH2F
        ("hConePtSumTrackPerNCellCut","Track #Sigma #it{p}_{T}, different #it{N}_{cell} cuts",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerNCellCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhConeSumPtTrackPerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerNCellCut) ;
        
        fhConeSumPtTrackPerNCellCutLargePtTrig = new TH2F
        ("hConePtSumTrackPerNCellCutLargePtTrig","Track #Sigma #it{p}_{T}, different #it{N}_{cell} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerNCellCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerNCellCutLargePtTrig->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhConeSumPtTrackPerNCellCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerNCellCutLargePtTrig) ;
        
        fhPtTrackInConePerNCellCut = new TH2F
        ("hPtTrackInConePerNCellCut","Track #it{p}_{T}, different #it{N}_{cell} cuts",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerNCellCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerNCellCut->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhPtTrackInConePerNCellCut->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhPtTrackInConePerNCellCut) ;
        
        fhPtTrackInConePerNCellCutLargePtTrig = new TH2F
        ("hPtTrackInConePerNCellCutLargePtTrig","Track #it{p}_{T}, different #it{N}_{cell} cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNNCellsInCandidate,0.5,fNNCellsInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerNCellCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerNCellCutLargePtTrig->SetXTitle("#it{N}_{cell}^{min}");
        for(Int_t i = 1; i <= fNNCellsInCandidate; i++)
          fhPtTrackInConePerNCellCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%d",fNCellsInCandidate[i-1]));
        outputContainer->Add(fhPtTrackInConePerNCellCutLargePtTrig) ;
        
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
      
      if(fStudyExoticTrigger)
      { 
        fhConeSumPtTrackPerExoCut = new TH2F
        ("hConePtSumTrackPerExoCut","Track #Sigma #it{p}_{T}, different exoticity cuts",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerExoCut->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerExoCut->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhConeSumPtTrackPerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerExoCut) ;
        
        fhConeSumPtTrackPerExoCutLargePtTrig = new TH2F
        ("hConePtSumTrackPerExoCutLargePtTrig","Track #Sigma #it{p}_{T}, different exoticity cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackPerExoCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackPerExoCutLargePtTrig->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhConeSumPtTrackPerExoCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhConeSumPtTrackPerExoCutLargePtTrig) ;
        
        fhPtTrackInConePerExoCut = new TH2F
        ("hPtTrackInConePerExoCut","Track #it{p}_{T}, different exoticity cuts",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerExoCut->SetYTitle("#it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerExoCut->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhPtTrackInConePerExoCut->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhPtTrackInConePerExoCut) ;
        
        fhPtTrackInConePerExoCutLargePtTrig = new TH2F
        ("hPtTrackInConePerExoCutLargePtTrig","Track #it{p}_{T}, different exoticity cuts, #it{p}_{T}^{trig} > 10 GeV",
         fNExoCutInCandidate,0.5,fNExoCutInCandidate+0.5,nptsumbins,ptsummin,ptsummax);
        fhPtTrackInConePerExoCutLargePtTrig->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhPtTrackInConePerExoCutLargePtTrig->SetXTitle("exoticity");
        for(Int_t i = 1; i <= fNExoCutInCandidate; i++)
          fhPtTrackInConePerExoCutLargePtTrig->GetXaxis()->SetBinLabel(i, Form("%2.2f",fExoCutInCandidate[i-1]));
        outputContainer->Add(fhPtTrackInConePerExoCutLargePtTrig) ;
      }

      
      fhConeSumPtTrack  = new TH2F
      ("hConePtSumTrack",
       Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrack) ;
      
      if(fStudyExoticTrigger)
      {
        fhConeSumPtTrackExoTrigger  = new TH2F
        ("hConePtSumTrackExoTrigger",
         Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f, exo trigger",r),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtTrackExoTrigger->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackExoTrigger->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtTrackExoTrigger) ;
      }
      
      if(fStudyTracksInCone)
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

        if(fStudyExoticTrigger)
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
                                       Form("#eta vs #varphi of Tracks in cone for #it{R} = %2.2f",r),
                                       netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiTrackInCone->SetXTitle("#eta");
        fhEtaPhiTrackInCone->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiTrackInCone) ;
        
        fhEtaTrackInCone = new TH2F("hEtaTrackInCone",
                                    Form("#eta vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f",r),
                                    nptbins,ptmin,ptmax,netabins,-1,1);
        fhEtaTrackInCone->SetYTitle("#eta");
        fhEtaTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaTrackInCone) ;

        fhPhiTrackInCone = new TH2F("hPhiTrackInCone",
                                    Form("#varphi vs #it{p}_{T} of Tracks in cone for #it{R} = %2.2f",r),
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
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUETrack  = new TH2F("hConePtSumEtaBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Eta Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUETrack) ;
        
        fhConeSumPtPhiBandUETrack  = new TH2F("hConePtSumPhiBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Phi Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtPhiBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPhiBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUETrack) ;
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhConeSumPtEtaBandUETrackTrigEtaPhi  = new TH2F("hConePtSumEtaBandUETrackTrigEtaPhi",
                                                          "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE Eta Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhi) ;
          
          fhConeSumPtPhiBandUETrackTrigEtaPhi  = new TH2F("hConePtSumPhiBandUETrackTrigEtaPhi",
                                                          "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE Phi Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhi) ;
          
          fhEtaBandTrackEtaPhi  = new TH2F("hEtaBandTrackEtaPhi",
                                           Form("#eta vs #varphi of tracks in #eta band isolation cone for #it{R} =  %2.2f",r),
                                           netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaBandTrackEtaPhi->SetXTitle("#eta");
          fhEtaBandTrackEtaPhi->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaBandTrackEtaPhi) ;
          
          fhPhiBandTrackEtaPhi  = new TH2F("hPhiBandTrackEtaPhi",
                                           Form("#eta vs #varphi of tracks in #varphi band isolation cone for #it{R} =  %2.2f",r),
                                           netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhPhiBandTrackEtaPhi->SetXTitle("#eta");
          fhPhiBandTrackEtaPhi->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhPhiBandTrackEtaPhi) ;
        }
        
        fhEtaBandTrackPt  = new TH2F("hEtaBandTrackPt",
                                   Form("#it{p}_{T} of tracks in #eta band isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhEtaBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhEtaBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandTrackPt) ;
        
        fhPhiBandTrackPt  = new TH2F("hPhiBandTrackPt",
                                   Form("#eta vs #varphi of tracks in #varphi band isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhPhiBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhPhiBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandTrackPt) ;
        
        fhConeSumPtEtaUENormTrack  = new TH2F("hConeSumPtEtaUENormTrack",
                                             Form("Tracks #Sigma #it{p}_{T} in normalized #eta band, #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUENormTrack->SetYTitle("#Sigma #it{p}_{T}^{#eta-band}_{norm} (GeV/#it{c})");
        fhConeSumPtEtaUENormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUENormTrack) ;
        
        fhConeSumPtPhiUENormTrack  = new TH2F("hConeSumPtPhiUENormTrack",
                                             Form("Tracks #Sigma #it{p}_{T} in normalized #varphi band, #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUENormTrack->SetYTitle("#Sigma #it{p}_{T}^{#varphi-band}_{norm} (GeV/#it{c})");
        fhConeSumPtPhiUENormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUENormTrack) ;

        
        fhConeSumPtEtaUESubTrack  = new TH2F("hConeSumPtEtaUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUESubTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubTrack) ;
        
        fhConeSumPtPhiUESubTrack  = new TH2F("hConeSumPtPhiUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUESubTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPhiUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubTrack) ;
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhConeSumPtEtaUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackTrigEtaPhi",
                                                         Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubTrackTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackTrigEtaPhi",
                                                         Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubTrackTrigEtaPhi) ;
          
          fhFractionTrackOutConeEta  = new TH2F("hFractionTrackOutConeEta",
                                                Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance",r),
                                                nptbins,ptmin,ptmax,100,0,1);
          fhFractionTrackOutConeEta->SetYTitle("#it{fraction}");
          fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionTrackOutConeEta) ;
          
          fhFractionTrackOutConeEtaTrigEtaPhi  = new TH2F("hFractionTrackOutConeEtaTrigEtaPhi",
                                                          Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance, in trigger #eta-#varphi ",r),
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
          
          fhConeSumPtSubvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubvsConeSumPtTotPhiTrack",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiTrack);
          
          fhConeSumPtSubNormvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiTrack",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiTrack);
          
          fhConeSumPtSubvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubvsConeSumPtTotEtaTrack",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaTrack);
          
          fhConeSumPtSubNormvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaTrack",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaTrack);
        }
        
        // UE in perpendicular cone
        fhPerpConeSumPt  = new TH2F("hPerpConePtSum",
                                    Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPerpConeSumPt->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeSumPt) ;
        
        fhPtInPerpCone  = new TH2F("hPtInPerpCone",
                                   Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInPerpCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtInPerpCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInPerpCone) ;
        
        if(fStudyTracksInCone)
        {
          fhEtaPhiInPerpCone= new TH2F("hEtaPhiInPerpCone",
                                       Form("#eta vs #varphi of all Tracks"),
                                       netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaPhiInPerpCone->SetXTitle("#eta");
          fhEtaPhiInPerpCone->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiInPerpCone) ;

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
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhEtaPhiTrack= new TH2F("hEtaPhiTrack",
                                  Form("#eta vs #varphi of all Tracks"),
                                  netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaPhiTrack->SetXTitle("#eta");
          fhEtaPhiTrack->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiTrack) ;
          
          fhEtaPhiInConeTrack= new TH2F("hEtaPhiInConeTrack",
                                        Form("#eta vs #varphi of Tracks in cone for #it{R} =  %2.2f",r),
                                        netabins,-1,1,nphibins,0,TMath::TwoPi());
          fhEtaPhiInConeTrack->SetXTitle("#eta");
          fhEtaPhiInConeTrack->SetYTitle("#varphi");
          outputContainer->Add(fhEtaPhiInConeTrack) ;
          
          fhConeSumPtVSUETracksEtaBand  = new TH2F("hConeSumPtVSUETracksEtaBand",
                                                   Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #eta band for tracks (before normalization), R=%2.2f",r),
                                                   nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
          fhConeSumPtVSUETracksEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUETracksEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtVSUETracksEtaBand);
          
          fhConeSumPtVSUETracksPhiBand  = new TH2F("hConeSumPtVSUETracksPhiBand",
                                                   Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #varphi band for tracks (before normalization), R=%2.2f",r),
                                                   nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
          fhConeSumPtVSUETracksPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUETracksPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtVSUETracksPhiBand);
        }
      }
    }
    
    if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged )
    {
      fhConeSumPtClustervsTrack   = new TH2F("hConePtSumClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
      fhConeSumPtClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClustervsTrack) ;

      fhConeSumPtClusterTrackFrac   = new TH2F("hConePtSumClusterTrackFraction",
                                             Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track} in isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,200,0,5);
      fhConeSumPtClusterTrackFrac->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
      fhConeSumPtClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterTrackFrac) ;

      
      fhConePtLeadClustervsTrack   = new TH2F("hConePtLeadClustervsTrack",
                                             Form("Track vs Cluster lead #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhConePtLeadClustervsTrack->SetXTitle("#it{p}^{leading cluster}_{T} (GeV/#it{c})");
      fhConePtLeadClustervsTrack->SetYTitle("#it{p}^{leading track}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadClustervsTrack) ;
      
      fhConePtLeadClusterTrackFrac   = new TH2F("hConePtLeadClusterTrackFraction",
                                               Form(" #it{p}^{leading cluster}_{T}/#it{p}^{leading track}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,200,0,5);
      fhConePtLeadClusterTrackFrac->SetYTitle("#it{p}^{leading cluster}_{T}/ #it{p}^{leading track}_{T}");
      fhConePtLeadClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadClusterTrackFrac) ;

      
      if(fFillCellHistograms)
      {
        fhConeSumPtCellvsTrack   = new TH2F("hConePtSumCellvsTrack",
                                            Form("Track vs cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
        fhConeSumPtCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCellvsTrack) ;
        
        fhConeSumPtCellTrack = new TH2F("hConePtSumCellTrack",
                                        Form("Track and Cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCellTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCellTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCellTrack) ;
        
        fhConeSumPtCellTrackTrigEtaPhi  = new TH2F("hConePtSumCellTrackTrigEtaPhi",
                                                   Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                                   netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtCellTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtCellTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtCellTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtCellTrackTrigEtaPhi) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUESub  = new TH2F("hConeSumPtEtaUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESub) ;
        
        fhConeSumPtPhiUESub  = new TH2F("hConeSumPtPhiUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtPhiUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESub) ;
        
        if(fFillUEBandSubtractHistograms > 1)
        {
          fhConeSumPtEtaUESubTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrigEtaPhi",
                                                    Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                    netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}(GeV/#it{c})");
          fhConeSumPtEtaUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrigEtaPhi",
                                                    Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                    netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}(GeV/#it{c})");
          fhConeSumPtPhiUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubTrigEtaPhi) ;
          
          fhConeSumPtEtaUESubClustervsTrack   = new TH2F("hConePtSumEtaUESubClustervsTrack",
                                                         Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #eta band in isolation cone for #it{R} =  %2.2f",r),
                                                         1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
          
          fhConeSumPtPhiUESubClustervsTrack   = new TH2F("hConePhiUESubPtSumClustervsTrack",
                                                         Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                                         1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
          
          fhEtaBandClustervsTrack   = new TH2F("hEtaBandClustervsTrack",
                                               Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",r),
                                               nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhEtaBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhEtaBandClustervsTrack) ;
          
          fhPhiBandClustervsTrack   = new TH2F("hPhiBandClustervsTrack",
                                               Form("Track vs Cluster #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                               nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
          fhPhiBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhPhiBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhPhiBandClustervsTrack) ;
          
          fhEtaBandNormClustervsTrack   = new TH2F("hEtaBandNormClustervsTrack",
                                                   Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",r),
                                                   nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhEtaBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhEtaBandNormClustervsTrack) ;
          
          fhPhiBandNormClustervsTrack   = new TH2F("hPhiBandNormClustervsTrack",
                                                   Form("Track vs Cluster #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                                   nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhPhiBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhPhiBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhPhiBandNormClustervsTrack) ;
        }
        
        if(fFillCellHistograms)
        {
          fhConeSumPtEtaUESubCellvsTrack   = new TH2F("hConePtSumEtaUESubCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub #eta band in isolation cone for #it{R} =  %2.2f",r),
                                                      1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtEtaUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhConeSumPtEtaUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubCellvsTrack) ;
          
          fhConeSumPtPhiUESubCellvsTrack   = new TH2F("hConePhiUESubPtSumCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                                      1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtPhiUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhConeSumPtPhiUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubCellvsTrack) ;
          
          fhEtaBandCellvsTrack   = new TH2F("hEtaBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhEtaBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhEtaBandCellvsTrack) ;
          
          fhPhiBandCellvsTrack   = new TH2F("hPhiBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
          fhPhiBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhPhiBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhPhiBandCellvsTrack) ;
          
          fhEtaBandNormCellvsTrack   = new TH2F("hEtaBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhEtaBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhEtaBandNormCellvsTrack) ;
          
          fhPhiBandNormCellvsTrack   = new TH2F("hPhiBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhPhiBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T}^{cell} (GeV/#it{c})");
          fhPhiBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhPhiBandNormCellvsTrack) ;
          
          fhConeSumPtEtaUESubTrackCell  = new TH2F("hConeSumPtEtaUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtEtaUESubTrackCell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCell) ;
          
          fhConeSumPtPhiUESubTrackCell  = new TH2F("hConeSumPtPhiUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
          fhConeSumPtPhiUESubTrackCell->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubTrackCell) ;
          
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubTrackCellTrigEtaPhi) ;
        }
      }
    }
    
    TString region[] = {"ITS","TPC","TRD","TOF","Top EMCal","In EMCal"}; // conversion regions

    for(Int_t iso = 0; iso < 2; iso++)
    {
      if(fFillTMHisto)
      {
        fhTrackMatchedDEta[iso]  = new TH2F
        (Form("hTrackMatchedDEta%s",isoName[iso].Data()),
         Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEta[iso]->SetYTitle("d#eta");
        fhTrackMatchedDEta[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDPhi%s",isoName[iso].Data()),
         Form("%s - d#varphi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhi[iso]->SetYTitle("d#varphi (rad)");
        fhTrackMatchedDPhi[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDEtaDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDEtaDPhi%s",isoName[iso].Data()),
         Form("%s - d#eta vs d#varphi of cluster-track, %s",isoTitle[iso].Data(),parTitle.Data()),
         nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDEtaDPhi[iso]->SetYTitle("d#varphi (rad)");
        fhTrackMatchedDEtaDPhi[iso]->SetXTitle("d#eta");
          
        outputContainer->Add(fhTrackMatchedDEta[iso]) ;
        outputContainer->Add(fhTrackMatchedDPhi[iso]) ;
        outputContainer->Add(fhTrackMatchedDEtaDPhi[iso]) ;
         
        if(IsDataMC())
        {
          for(int imc = 0; imc < fgkNmcTypes; imc++)
          {
            fhTrackMatchedDEtaMC[imc][iso] = new TH2F
            (Form("hTrackMatchedDEta%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
             nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
            fhTrackMatchedDEtaMC[imc][iso]->SetYTitle("d#eta");
            fhTrackMatchedDEtaMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
            
            fhTrackMatchedDPhiMC[imc][iso] = new TH2F
            (Form("hTrackMatchedDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
             Form("%s - d#varphi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
             nptbins,ptmin,ptmax,nresetabins,resphimin,resphimax);
            fhTrackMatchedDPhiMC[imc][iso]->SetYTitle("d#varphi (rad)");
            fhTrackMatchedDPhiMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
              
              fhTrackMatchedDEtaDPhiMC[imc][iso]  = new TH2F
              (Form("hTrackMatchedDEtaDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s - d#eta vs d#varphi of cluster-track, %s",isoTitle[iso].Data(),parTitle.Data()),
               nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
              fhTrackMatchedDEtaDPhiMC[imc][iso]->SetYTitle("d#varphi (rad)");
              fhTrackMatchedDEtaDPhiMC[imc][iso]->SetXTitle("d#eta");
              
              outputContainer->Add(fhTrackMatchedDEtaMC[imc][iso]) ;
              outputContainer->Add(fhTrackMatchedDPhiMC[imc][iso]) ;
              outputContainer->Add(fhTrackMatchedDEtaDPhiMC[imc][iso]);
          }
          
        }
        
        fhdEdx[iso]  = new TH2F
        (Form("hdEdx%s",isoName[iso].Data()),
         Form("%s - Matched track <d#it{E}/d#it{x}> vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
        fhdEdx[iso]->SetXTitle("#it{E} (GeV)");
        fhdEdx[iso]->SetYTitle("<d#it{E}/d#it{x}>");
        outputContainer->Add(fhdEdx[iso]);
        
        fhEOverP[iso]  = new TH2F
        (Form("hEOverP%s",isoName[iso].Data()),
         Form("%s - Matched track #it{E}/#it{p} vs cluster, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
        fhEOverP[iso]->SetXTitle("#it{E} (GeV)");
        fhEOverP[iso]->SetYTitle("#it{E}/#it{p}");
        outputContainer->Add(fhEOverP[iso]);
        
        if(IsDataMC())
        {
          fhTrackMatchedMCParticle[iso]  = new TH2F
          (Form("hTrackMatchedMCParticle%s",isoName[iso].Data()),
           Form("%s - Origin of particle vs cluster #it{E}, %s",isoTitle[iso].Data(),parTitle.Data()),
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
      
      if(fFillSSHisto)
      {
        fhELambda0[iso]  = new TH2F
        (Form("hELambda0%s",isoName[iso].Data()),
         Form("%s cluster : #it{E} vs #lambda_{0}, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhELambda0[iso]->SetYTitle("#lambda_{0}^{2}");
        fhELambda0[iso]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhELambda0[iso]) ;
        
//        fhELambda1[iso]  = new TH2F
//        (Form("hELambda1%s",isoName[iso].Data()),
//         Form("%s cluster: #it{E} vs #lambda_{1}, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//        fhELambda1[iso]->SetYTitle("#lambda_{1}^{2}");
//        fhELambda1[iso]->SetXTitle("#it{E} (GeV)");
//        outputContainer->Add(fhELambda1[iso]) ;
        
        fhPtLambda0[iso]  = new TH2F
        (Form("hPtLambda0%s",isoName[iso].Data()),
         Form("%s cluster : #it{p}_{T} vs #lambda_{0}, %s",isoTitle[iso].Data(), parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhPtLambda0[iso]->SetYTitle("#lambda_{0}^{2}");
        fhPtLambda0[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLambda0[iso]) ;
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
          {
            fhPtLambda0MC[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
                                                Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MC[imc][iso]->SetYTitle("#lambda_{0}^{2}");
            fhPtLambda0MC[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MC[imc][iso]) ;

            fhPtLambda0MCConv[imc][iso]  = new TH2F
            (Form("hPtLambda0%s_MC%sConv",isoName[iso].Data(),mcPartName[imc].Data()),
                                                Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCConv[imc][iso]->SetYTitle("#lambda_{0}^{2}");
            fhPtLambda0MCConv[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCConv[imc][iso]) ;
            
            if(fFillOverlapHistograms)
            {
              fhPtLambda0MCWith1Overlap[imc][iso]  = new TH2F
              (Form("hPtLambda0%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, 1 overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCWith1Overlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCWith1Overlap[imc][iso]) ;
              
              fhPtLambda0MCConvWith1Overlap[imc][iso]  = new TH2F
              (Form("hPtLambda0%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, 1 overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCConvWith1Overlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCConvWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCConvWith1Overlap[imc][iso]) ;
              
              fhPtLambda0MCWithNoOverlap[imc][iso]  = new TH2F
              (Form("hPtLambda0%s_MC%s_NoOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, no overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCWithNoOverlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCWithNoOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCWithNoOverlap[imc][iso]) ;
              
              fhPtLambda0MCConvWithNoOverlap[imc][iso]  = new TH2F
              (Form("hPtLambda0%s_MC%sConv_NoOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, no overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCConvWithNoOverlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCConvWithNoOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCConvWithNoOverlap[imc][iso]) ;
              
              fhPtNOverlap[imc][iso]  = new TH2F
              (Form("hPtNOverlaps%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, 1 overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,10,0,10);
              fhPtNOverlap[imc][iso]->SetYTitle("#it{N} overlaps");
              fhPtNOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtNOverlap[imc][iso]) ;
              
              fhPtNOverlapConv[imc][iso]  = new TH2F
              (Form("hPtNOverlaps%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, 1 overlap",
                    isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
               nptbins,ptmin,ptmax,10,0,10);
              fhPtNOverlapConv[imc][iso]->SetYTitle("#it{N} overlaps");
              fhPtNOverlapConv[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtNOverlapConv[imc][iso]) ;
            }
          
          }
        }
        
        if(fIsoDetector==kEMCAL &&  GetFirstSMCoveredByTRD() >= 0)
        {
          fhPtLambda0TRD[iso]  = new TH2F
          (Form("hPtLambda0TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{p}_{T} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhPtLambda0TRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0TRD[iso]) ;
          
          fhELambda0TRD[iso]  = new TH2F
          (Form("hELambda0TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{E} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0TRD[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0TRD[iso]) ;
          
//          fhELambda1TRD[iso]  = new TH2F
//          (Form("hELambda1TRD%s",isoName[iso].Data()),
//           Form("%s cluster: #it{E} vs #lambda_{1}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//          fhELambda1TRD[iso]->SetYTitle("#lambda_{1}^{2}");
//          fhELambda1TRD[iso]->SetXTitle("#it{E} (GeV)");
//          outputContainer->Add(fhELambda1TRD[iso]) ;
        }
        
        if(fFillNLMHistograms)
        {
          fhNLocMax[iso] = new TH2F
          (Form("hNLocMax%s",isoName[iso].Data()),
           Form("%s - Number of local maxima in cluster, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,10,0,10);
          fhNLocMax[iso]->SetYTitle("#it{NLM}");
          fhNLocMax[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhNLocMax[iso]) ;
          
          fhELambda0LocMax1[iso]  = new TH2F
          (Form("hELambda0LocMax1%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax1[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax1[iso]) ;
          
          fhELambda1LocMax1[iso]  = new TH2F
          (Form("hELambda1LocMax1%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax1[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax1[iso]) ;
          
          fhELambda0LocMax2[iso]  = new TH2F
          (Form("hELambda0LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax2[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax2[iso]) ;
          
          fhELambda1LocMax2[iso]  = new TH2F
          (Form("hELambda1LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax2[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax2[iso]) ;
          
          fhELambda0LocMaxN[iso]  = new TH2F
          ( Form("hELambda0LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMaxN[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMaxN[iso]) ;
          
          fhELambda1LocMaxN[iso]  = new TH2F
          (Form("hELambda1LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),
           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMaxN[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMaxN[iso]) ;
        } // NLM
      } // SS histo
      
      if(GetCalorimeter() == kEMCAL && fFillEMCALRegionHistograms)
      {
        for(Int_t ieta = 0; ieta < 4; ieta++) 
        {  
          for(Int_t iphi = 0; iphi < 3; iphi++) 
          {
//            fhLam0EMCALRegion[iso][ieta][iphi] = 
//            new TH2F(Form("hLam0_%s_eta%d_phi%d",isoName[iso].Data(),ieta,iphi),
//                     Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, region #eta %d, #varphi %d",
//                           isoTitle[iso].Data(),ieta,iphi),
//                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//            fhLam0EMCALRegion[iso][ieta][iphi]->SetYTitle("#lambda_{0}^{2}");
//            fhLam0EMCALRegion[iso][ieta][iphi]->SetXTitle("#it{p}_{T} (GeV)");
//            outputContainer->Add(fhLam0EMCALRegion[iso][ieta][iphi]) ;
//            
//            if(GetFirstSMCoveredByTRD() >= 0)
//            {
//              fhLam0EMCALRegionTRD[iso][ieta][iphi] = 
//              new TH2F(Form("hLam0TRD_%s_eta%d_phi%d",isoName[iso].Data(),ieta,iphi),
//                       Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, region #eta %d, #varphi %d, SM covered by TRD",
//                             isoTitle[iso].Data(),ieta,iphi),
//                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//              fhLam0EMCALRegionTRD[iso][ieta][iphi]->SetYTitle("#lambda_{0}^{2}");
//              fhLam0EMCALRegionTRD[iso][ieta][iphi]->SetXTitle("#it{p}_{T} (GeV)");
//              outputContainer->Add(fhLam0EMCALRegionTRD[iso][ieta][iphi]) ;
//            } // TRD
            
            for(Int_t ism = 0; ism < fNModules; ism++) 
            {
              if(ism < fFirstModule || ism> fLastModule) continue;

              fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism] = 
              new TH2F(Form("hLam0_%s_eta%d_phi%d_sm%d",isoName[iso].Data(),ieta,iphi,ism),
                       Form("%s, #it{p}_{T} vs #lambda_{0}^{2}, sm %d, region #eta %d, #varphi %d",
                            isoTitle[iso].Data(),ism,ieta,iphi),
                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]->SetYTitle("#lambda_{0}^{2}");
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
      
      if(IsDataMC() && fStudyMCConversionRadius)
      {
        fhMCConversionVertex[iso] = new TH2F(Form("hMCPhotonConversionVertex%s",isoName[iso].Data()),
                                             Form("%s, cluster from converted photon, #it{p}_{T} vs vertex distance, %s",
                                                  isoTitle[iso].Data(),parTitle.Data()),
                                             nptbins,ptmin,ptmax,500,0,500);
        fhMCConversionVertex[iso]->SetYTitle("#it{R} (cm)");
        fhMCConversionVertex[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhMCConversionVertex[iso]) ;
        
        if(GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0)
        {
          fhMCConversionVertexTRD[iso] = new TH2F(Form("hMCPhotonConversionVertexTRD%s",isoName[iso].Data()),
                                                  Form("%s, cluster from converted photon, #it{p}_{T} vs vertex distance, %s, SM covered by TRD",
                                                       isoTitle[iso].Data(),parTitle.Data()),
                                                  nptbins,ptmin,ptmax,500,0,500);
          fhMCConversionVertexTRD[iso]->SetYTitle("#it{R} (cm)");
          fhMCConversionVertexTRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhMCConversionVertexTRD[iso]) ;
        }
        
        for(Int_t iR = 0; iR < 6; iR++)
        {
          fhMCConversionLambda0Rcut[iR][iso] = new TH2F(Form("hMCPhotonConversionLambda0%s_R%d",isoName[iso].Data(),iR),
                                                        Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s",
                                                             isoTitle[iso].Data(),region[iR].Data()),
                                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhMCConversionLambda0Rcut[iR][iso]->SetYTitle("#lambda_{0}^{2}");
          fhMCConversionLambda0Rcut[iR][iso]->SetXTitle("#it{p}_{T} (GeV)");
          outputContainer->Add(fhMCConversionLambda0Rcut[iR][iso]) ;
          
          if(GetCalorimeter() == kEMCAL &&  GetFirstSMCoveredByTRD() >= 0)
          {
            fhMCConversionLambda0RcutTRD[iR][iso] = new TH2F(Form("hMCPhotonConversionLambda0TRD%s_R%d",isoName[iso].Data(),iR),
                                                             Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s, SM covered by TRD",
                                                                  isoTitle[iso].Data(),region[iR].Data()),
                                                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhMCConversionLambda0RcutTRD[iR][iso]->SetYTitle("#lambda_{0}^{2}");
            fhMCConversionLambda0RcutTRD[iR][iso]->SetXTitle("#it{p}_{T} (GeV)");
            outputContainer->Add(fhMCConversionLambda0RcutTRD[iR][iso]) ;
          }
        }
          
//          if(GetCalorimeter() == kEMCAL && fFillEMCALRegionHistograms)
//          {
//            for(Int_t ieta = 0; ieta < 4; ieta++) 
//            {  
//              for(Int_t iphi = 0; iphi < 3; iphi++) 
//              {
//                for(Int_t iReg = 0; iReg < 6; iReg++) 
//                {
//                  fhLam0EMCALRegionMCConvRcut[iso][ieta][iphi][iReg] = 
//                  new TH2F(Form("hMCPhotonConversionLambda0%s_R%d_eta%d_phi%d",isoName[iso].Data(),iReg,ieta,iphi),
//                           Form("%s,cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s, region #eta %d, #varphi %d",
//                                isoTitle[iso].Data(),region[iReg].Data(),ieta,iphi),
//                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//                  fhLam0EMCALRegionMCConvRcut[iso][ieta][iphi][iReg]->SetYTitle("#lambda_{0}^{2}");
//                  fhLam0EMCALRegionMCConvRcut[iso][ieta][iphi][iReg]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
//                  outputContainer->Add(fhLam0EMCALRegionMCConvRcut[iso][ieta][iphi][iReg]) ;
//                  
//                  if(GetFirstSMCoveredByTRD() >= 0)
//                  {
//                    fhLam0EMCALRegionTRDMCConvRcut[iso][ieta][iphi][iReg] = 
//                    new TH2F(Form("hMCPhotonConversionLambda0TRD%s_R%d_eta%d_phi%d",isoName[iso].Data(),iReg,ieta,iphi),
//                             Form("%s,cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s, region #eta %d, #varphi %d, SM covered by TRD",
//                                  isoTitle[iso].Data(),region[iReg].Data(),ieta,iphi),
//                             nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//                    fhLam0EMCALRegionTRDMCConvRcut[iso][ieta][iphi][iReg]->SetYTitle("#lambda_{0}^{2}");
//                    fhLam0EMCALRegionTRDMCConvRcut[iso][ieta][iphi][iReg]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
//                    outputContainer->Add(fhLam0EMCALRegionTRDMCConvRcut[iso][ieta][iphi][iReg]) ;
//                  } // TRD
//                  
//                } // iR
//              } // iphi 
//            } // ieta
//          } // regions in EMCal
      } // conversion radius
      
    } // control histograms for isolated and non isolated objects
    
    if(IsPileUpAnalysisOn())
    {
      if(fStudyTracksInCone)
      {
        fhPtTrackInConeOtherBCPileUpSPD  = new TH2F("hPtTrackInConeOtherBCPileUpSPD",
                                                    Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0, pile-up from SPD",r),
                                                    nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeOtherBCPileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtTrackInConeOtherBCPileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeOtherBCPileUpSPD) ;
        
        fhPtTrackInConeBC0PileUpSPD  = new TH2F("hPtTrackInConeBC0PileUpSPD",
                                                Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0, pile-up from SPD",r),
                                                nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtTrackInConeBC0PileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtTrackInConeBC0PileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtTrackInConeBC0PileUpSPD) ;
      }
      
      for (Int_t i = 0; i < 7 ; i++)
      {
        fhPtInConePileUp[i]  = new TH2F(Form("hPtInConePileUp%s",pileUpName[i].Data()),
                                        Form("#it{p}_{T} in isolation cone for #it{R} =  %2.2f, from pile-up (%s)",r,pileUpName[i].Data()),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInConePileUp[i]->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtInConePileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInConePileUp[i]) ;
      }
    }
    
    if(IsDataMC())
    {
      // For histograms in arrays, index in the array, corresponding to any particle origin
      
      for(Int_t i = 0; i < fgkNmcPrimTypes; i++)
      {
        fhEPrimMC[i]  = new TH1F(Form("hEPrim_MC%s",ppname[i].Data()),
                                 Form("primary photon  %s : #it{E}, %s",pptype[i].Data(),parTitle.Data()),
                                 nptbins,ptmin,ptmax);
        fhEPrimMC[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEPrimMC[i]) ;

        fhPtPrimMC[i]  = new TH1F(Form("hPtPrim_MC%s",ppname[i].Data()),
                                 Form("primary photon  %s : #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                 nptbins,ptmin,ptmax);
        fhPtPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMC[i]) ;

        fhPtPrimMCiso[i]  = new TH1F(Form("hPtPrim_MCiso%s",ppname[i].Data()),
                                     Form("primary isolated photon %s : #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                     nptbins,ptmin,ptmax);
        fhPtPrimMCiso[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCiso[i]) ;
        
        fhEtaPrimMC[i]  = new TH2F(Form("hEtaPrim_MC%s",ppname[i].Data()),
                                   Form("primary photon %s : #eta vs #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax,200,-2,2);
        fhEtaPrimMC[i]->SetYTitle("#eta");
        fhEtaPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaPrimMC[i]) ;
        
        fhPhiPrimMC[i]  = new TH2F(Form("hPhiPrim_MC%s",ppname[i].Data()),
                                   Form("primary photon %s : #varphi vs #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax,200,0.,TMath::TwoPi());
        fhPhiPrimMC[i]->SetYTitle("#varphi (rad)");
        fhPhiPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPhiPrimMC[i]) ;
      }
      
      if(fMakePrimaryPi0DecayStudy)
      {
        fhPtPrimMCPi0DecayPairAcceptInConeLowPt  = new TH1F("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPt",
                                                            Form("primary photon  %s : #it{p}_{T}, pair in cone, %s",pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                            nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPt) ;
        
        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPt",
                                                               Form("isolated primary photon %s, pair in cone : #it{p}_{T}, %s",
                                                                    pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                               nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt) ;
        
        fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap  = new TH1F("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPtNoOverlap",
                                                                     Form("primary photon  %s, no overlap, pair in cone : #it{p}_{T}, %s",
                                                                          pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                                     nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap) ;
        
        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPtNoOverlap",
                                                                        Form("isolated primary photon  %s, pair in cone,no overlap : #it{p}_{T}, %s",
                                                                             pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                                        nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap) ;

        fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F("hPtPrim_MCPhotonPi0DecayPairAcceptInConeLowPtNoOverlapCaloE",
                                                                     Form("primary photon  %s, no overlap, pair in cone, E > calo min: #it{p}_{T}, %s",
                                                                          pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                                     nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE) ;

        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairAcceptInConeLowPtNoOverlapCaloE",
                                                                        Form("isolated primary photon  %s, pair in cone,no overlap, E > calo min: #it{p}_{T}, %s",
                                                                             pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                                        nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE) ;

        
        fhPtPrimMCPi0DecayPairNoOverlap  = new TH1F("hPtPrim_MCPhotonPi0DecayPairNoOverlap",
                                                                     Form("primary photon  %s, no overlap: #it{p}_{T}, %s",
                                                                          pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                                     nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairNoOverlap) ;

        fhPtPrimMCPi0DecayIsoPairNoOverlap  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairNoOverlap",
                                                    Form("isolated primary photon  %s, no overlap: #it{p}_{T}, %s",
                                                         pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                    nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairNoOverlap) ;
        
        fhPtPrimMCPi0DecayPairOutOfCone  = new TH1F("hPtPrim_MCPhotonPi0DecayPairOutOfCone",
                                                    Form("primary photon %s : #it{p}_{T}, pair out of cone, %s",pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                    nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfCone) ;
        
        fhPtPrimMCPi0DecayIsoPairOutOfCone  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairOutOfCone",
                                                       Form("isolated primary photon %s, pair out of cone : #it{p}_{T}, %s",
                                                            pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                       nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfCone) ;
        
        fhPtPrimMCPi0DecayPairOutOfAcceptance  = new TH1F("hPtPrim_MCPhotonPi0DecayPairOutOfAcceptance",
                                                          Form("primary photon %s : #it{p}_{T}, pair out of acceptance, %s",pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                          nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfAcceptance) ;
        
        fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap  = new TH1F("hPtPrim_MCPhotonPi0DecayPairOutOfAcceptanceNoOverlap",
                                                          Form("primary photon %s : #it{p}_{T}, pair out of acceptance, no overlap, %s",pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                          nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap) ;
        
        fhPtPrimMCPi0DecayIsoPairOutOfAcceptance  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairOutOfAcceptance",
                                                             Form("isolated primary photon %s, pair out of acceptance : #it{p}_{T}, %s",
                                                                  pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                             nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfAcceptance) ;
 
        fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap  = new TH1F("hPtPrim_MCisoPhotonPi0DecayPairOutOfAcceptanceNoOverlap",
                                                             Form("isolated primary photon %s, pair out of acceptance, no overlap : #it{p}_{T}, %s",
                                                                  pptype[kmcPrimPi0Decay].Data(),parTitle.Data()),
                                                             nptbins,ptmin,ptmax);
        fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap) ;
        
        fhPtPrimMCPi0Overlap  = new TH1F("hPtPrim_MCPi0Overlap",
                                         Form("primary %s, overlap: #it{p}_{T}, %s",
                                              pptype[kmcPrimPi0].Data(),parTitle.Data()),
                                         nptbins,ptmin,ptmax);
        fhPtPrimMCPi0Overlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0Overlap) ;
        
        fhPtPrimMCPi0IsoOverlap  = new TH1F("hPtPrim_MCisoPi0Overlap",
                                         Form("primary %s, overlap: #it{p}_{T}, %s",
                                              pptype[kmcPrimPi0].Data(),parTitle.Data()),
                                         nptbins,ptmin,ptmax);
        fhPtPrimMCPi0IsoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCPi0IsoOverlap) ;

        
        fhPtPrimMCEtaDecayPairAcceptInConeLowPt  = new TH1F("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPt",
                                                            Form("primary photon  %s : #it{p}_{T}, pair in cone, %s",pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                            nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPt) ;
        
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPt",
                                                               Form("isolated primary photon %s, pair in cone : #it{p}_{T}, %s",
                                                                    pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                               nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt) ;
        
        fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap  = new TH1F("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPtNoOverlap",
                                                                     Form("primary photon  %s, no overlap, pair in cone : #it{p}_{T}, %s",
                                                                          pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                     nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap) ;
        
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPtNoOverlap",
                                                                        Form("isolated primary photon  %s, pair in cone,no overlap : #it{p}_{T}, %s",
                                                                             pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                        nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap) ;
        
        fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F("hPtPrim_MCPhotonEtaDecayPairAcceptInConeLowPtNoOverlapCaloE",
                                                                          Form("primary photon  %s, no overlap, pair in cone, E > calo min: #it{p}_{T}, %s",
                                                                               pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                          nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE) ;
        
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairAcceptInConeLowPtNoOverlapCaloE",
                                                                             Form("isolated primary photon  %s, pair in cone,no overlap, E > calo min: #it{p}_{T}, %s",
                                                                                  pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                             nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE) ;
        
        
        fhPtPrimMCEtaDecayPairNoOverlap  = new TH1F("hPtPrim_MCPhotonEtaDecayPairNoOverlap",
                                                    Form("primary photon  %s, no overlap: #it{p}_{T}, %s",
                                                         pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                    nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairNoOverlap) ;
        
        fhPtPrimMCEtaDecayIsoPairNoOverlap  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairNoOverlap",
                                                       Form("isolated primary photon  %s, no overlap: #it{p}_{T}, %s",
                                                            pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                       nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairNoOverlap) ;
        
        fhPtPrimMCEtaDecayPairOutOfCone  = new TH1F("hPtPrim_MCPhotonEtaDecayPairOutOfCone",
                                                    Form("primary photon %s : #it{p}_{T}, pair out of cone, %s",pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                    nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfCone) ;
        
        fhPtPrimMCEtaDecayIsoPairOutOfCone  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairOutOfCone",
                                                       Form("isolated primary photon %s, pair out of cone : #it{p}_{T}, %s",
                                                            pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                       nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairOutOfCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfCone) ;
        
        fhPtPrimMCEtaDecayPairOutOfAcceptance  = new TH1F("hPtPrim_MCPhotonEtaDecayPairOutOfAcceptance",
                                                          Form("primary photon %s : #it{p}_{T}, pair out of acceptance, %s",pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                          nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfAcceptance) ;
        
        fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap  = new TH1F("hPtPrim_MCPhotonEtaDecayPairOutOfAcceptanceNoOverlap",
                                                                   Form("primary photon %s : #it{p}_{T}, pair out of acceptance, no overlap, %s",pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                   nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap) ;
        
        fhPtPrimMCEtaDecayIsoPairOutOfAcceptance  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairOutOfAcceptance",
                                                             Form("isolated primary photon %s, pair out of acceptance : #it{p}_{T}, %s",
                                                                  pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                             nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairOutOfAcceptance->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfAcceptance) ;
        
        fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap  = new TH1F("hPtPrim_MCisoPhotonEtaDecayPairOutOfAcceptanceNoOverlap",
                                                                      Form("isolated primary photon %s, pair out of acceptance, no overlap : #it{p}_{T}, %s",
                                                                           pptype[kmcPrimEtaDecay].Data(),parTitle.Data()),
                                                                      nptbins,ptmin,ptmax);
        fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap) ;
        
        fhPtPrimMCEtaOverlap  = new TH1F("hPtPrim_MCEtaOverlap",
                                         Form("primary %s, overlap: #it{p}_{T}, %s",
                                              pptype[kmcPrimEta].Data(),parTitle.Data()),
                                         nptbins,ptmin,ptmax);
        fhPtPrimMCEtaOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaOverlap) ;
        
        fhPtPrimMCEtaIsoOverlap  = new TH1F("hPtPrim_MCisoEtaOverlap",
                                            Form("primary %s, overlap: #it{p}_{T}, %s",
                                                 pptype[kmcPrimEta].Data(),parTitle.Data()),
                                            nptbins,ptmin,ptmax);
        fhPtPrimMCEtaIsoOverlap->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtPrimMCEtaIsoOverlap) ;

      }
      
    }//Histos with MC
    
  }
  
  if(fMakeSeveralIC)
  {
    const Int_t buffersize = 255;
    char name[buffersize];
    char title[buffersize];
    for(Int_t icone = 0; icone<fNCones; icone++)
    {
      // sum pt in cone vs. pt leading
      snprintf(name, buffersize,"hSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}#it{p}_{T} (GeV/#it{c})");//#Sigma #it{p}_{T}
      fhSumPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhSumPtLeadingPt[icone]) ;
      
      // pt in cone vs. pt leading
      snprintf(name, buffersize,"hPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtLeadingPt[icone] ->SetYTitle("#it{p}_{T}^{cone} (GeV/#it{c})");
      fhPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPtLeadingPt[icone]) ;
      
      // sum pt in cone vs. pt leading in the forward region (for background subtraction studies)
      snprintf(name, buffersize,"hPerpSumPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPerpSumPtLeadingPt[icone]  = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhPerpSumPtLeadingPt[icone] ->SetYTitle("#sum_{cone}#it{p}_{T} (GeV/#it{c})");//#Sigma #it{p}_{T}
      fhPerpSumPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPerpSumPtLeadingPt[icone]) ;
      
      // pt in cone vs. pt leading in the forward region (for background subtraction studies)
      snprintf(name, buffersize,"hPerpPtLeadingPt_Cone_%d",icone);
      snprintf(title, buffersize,"#it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSizes[icone]);
      fhPerpPtLeadingPt[icone]  = new TH2F(name, title,  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPerpPtLeadingPt[icone] ->SetYTitle("#it{p}_{T}^{cone} (GeV/#it{c})");
      fhPerpPtLeadingPt[icone] ->SetXTitle("#it{p}_{T}^{leading} (GeV/#it{c})");
      outputContainer->Add(fhPerpPtLeadingPt[icone]) ;
      
      if(IsDataMC())
      {
        for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
        {
          snprintf(name , buffersize,"hSumPtLeadingPt_MC%s_Cone_%d",mcPartName[imc].Data(),icone);
          snprintf(title, buffersize,"Candidate %s #it{p}_{T} vs cone #Sigma #it{p}_{T} for #it{R}=%2.2f",mcPartType[imc].Data(),fConeSizes[icone]);
          fhSumPtLeadingPtMC[imc][icone] = new TH2F(name, title,nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhSumPtLeadingPtMC[imc][icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhSumPtLeadingPtMC[imc][icone]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
          outputContainer->Add(fhSumPtLeadingPtMC[imc][icone]) ;
        }
      }//Histos with MC
      
      for(Int_t ipt = 0; ipt<fNPtThresFrac;ipt++)
      {
        snprintf(name, buffersize,"hPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
        fhPtThresIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtThresIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtThresIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        fhPtFracIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtFracIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hSumPt_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhSumPtIsolated[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
        // fhSumPtIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhSumPtIsolated[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhSumPtIsolated[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for density in #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhPtSumDensityIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPtSumDensityIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hPtFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #it{p}_{T} distribution for PtFracPtSum in #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhPtFracPtSumIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
        //fhPtSumIsolated[icone][ipt]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPtFracPtSumIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtFracPtSumIso[icone][ipt]) ;
        
        // eta:phi
        snprintf(name, buffersize,"hEtaPhiPtThres_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
        fhEtaPhiPtThresIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtThresIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtThresIso[icone][ipt]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiPtThresIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiPtFracIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtFracIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtFracIso[icone][ipt]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiPtFracIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtSumIso[icone][ipt]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiPtSumIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#varphi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiSumDensityIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiSumDensityIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiSumDensityIso[icone][ipt]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#varphi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiFracPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiFracPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiFracPtSumIso[icone][ipt]->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiFracPtSumIso[icone][ipt]) ;
        
        if(fFillTaggedDecayHistograms)
        {
          // pt decays isolated
          snprintf(name, buffersize,"hPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
          fhPtPtThresDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtPtThresDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtThresDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhPtPtFracDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);
          fhPtPtFracDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtFracDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhPtPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtPtSumDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtPtSumDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for density in #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhPtSumDensityDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtSumDensityDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtSumDensityDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hPtFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #it{p}_{T} distribution for PtFracPtSum in #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhPtFracPtSumDecayIso[icone][ipt]  = new TH1F(name, title,nptbins,ptmin,ptmax);//,nptsumbins,ptsummin,ptsummax);
          //  fhPtPtSumDecayIso[icone]->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhPtFracPtSumDecayIso[icone][ipt]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtFracPtSumDecayIso[icone][ipt]) ;
          
          // eta:phi decays
          snprintf(name, buffersize,"hEtaPhiPtThres_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
          fhEtaPhiPtThresDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiPtThresDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiPtFracDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiPtFracDecayIso[icone][ipt]) ;
          
          
          snprintf(name, buffersize,"hEtaPhiPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#varphi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiPtSumDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#varphi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiSumDensityDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiSumDensityDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#varphi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetYTitle("#varphi (rad)");
          outputContainer->Add(fhEtaPhiFracPtSumDecayIso[icone][ipt]) ;
          
        }
        
        if(IsDataMC())
        {
          for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
          {
            snprintf(name , buffersize,"hPtThreshMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #it{p}_{T}^{th}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fPtThresholds[ipt]);
            fhPtThresIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhPtThresIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhPtThresIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtThresIsolatedMC[imc][icone][ipt]) ;
            
            
            snprintf(name , buffersize,"hPtFracMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #Sigma #it{p}_{T}^{in cone}/#it{p}_{T}^{trig}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fPtFractions[ipt]);
            fhPtFracIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhPtFracIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhPtFracIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhPtFracIsolatedMC[imc][icone][ipt]) ;
            
            snprintf(name , buffersize,"hSumPtMC%s_Cone_%d_Pt%d",mcPartName[imc].Data(),icone,ipt);
            snprintf(title, buffersize,"Isolated %s #it{p}_{T} for #it{R}=%2.2f and #Sigma #it{p}_{T}^{in cone}=%2.2f",
                     mcPartType[imc].Data(),fConeSizes[icone], fSumPtThresholds[ipt]);
            fhSumPtIsolatedMC[imc][icone][ipt] = new TH1F(name, title,nptbins,ptmin,ptmax);
            fhSumPtIsolatedMC[imc][icone][ipt]->SetYTitle("#it{counts}");
            fhSumPtIsolatedMC[imc][icone][ipt]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add(fhSumPtIsolatedMC[imc][icone][ipt]) ;
          }
        }//Histos with MC
      }//icone loop
    }//ipt loop
  }
  
  if(IsPileUpAnalysisOn())
  {
    for (Int_t i = 0; i < 7 ; i++)
    {
      fhEIsoPileUp[i]   = new TH1F(Form("hEPileUp%s",pileUpName[i].Data()),
                                   Form("Number of isolated particles vs E, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhEIsoPileUp[i]->SetYTitle("d#it{N} / d#it{E}");
      fhEIsoPileUp[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEIsoPileUp[i]) ;
      
      fhPtIsoPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                                   Form("Number of isolated particles vs #it{p}_{T}, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhPtIsoPileUp[i]->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtIsoPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtIsoPileUp[i]) ;
      
      fhENoIsoPileUp[i]   = new TH1F(Form("hENoIsoPileUp%s",pileUpName[i].Data()),
                                     Form("Number of not isolated particles vs E, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                     nptbins,ptmin,ptmax);
      fhENoIsoPileUp[i]->SetYTitle("d#it{N} / dE");
      fhENoIsoPileUp[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhENoIsoPileUp[i]) ;
      
      fhPtNoIsoPileUp[i]  = new TH1F(Form("hPtNoIsoPileUp%s",pileUpName[i].Data()),
                                     Form("Number of not isolated particles vs #it{p}_{T}, %s, pile-up event by %s",parTitle.Data(),pileUpName[i].Data()),
                                     nptbins,ptmin,ptmax);
      fhPtNoIsoPileUp[i]->SetYTitle("d#it{N} / #it{p}_{T}");
      fhPtNoIsoPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtNoIsoPileUp[i]) ;
    }
    
    fhTimeENoCut  = new TH2F ("hTimeE_NoCut","time of cluster vs E of clusters, no cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeENoCut->SetXTitle("#it{E} (GeV)");
    fhTimeENoCut->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeENoCut);
    
    fhTimeESPD  = new TH2F ("hTimeE_SPD","time of cluster vs E of clusters, SPD cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPD->SetXTitle("#it{E} (GeV)");
    fhTimeESPD->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPD);
    
    fhTimeESPDMulti  = new TH2F ("hTimeE_SPDMulti","time of cluster vs E of clusters, SPD multi cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeESPDMulti->SetXTitle("#it{E} (GeV)");
    fhTimeESPDMulti->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeESPDMulti);
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 50,0,50 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
    fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
    fhTimeNPileUpVertContributors->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertContributors);
    
    fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDistance->SetYTitle("distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDistance->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDistance);
    
    fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",ntimebins,timemin,timemax,100,0,50);
    fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance #it{z} (cm) ");
    fhTimePileUpMainVertexZDiamond->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePileUpMainVertexZDiamond);
  }
    
  return outputContainer ;
}

//____________________________________________________
/// MC histogram index depending on origin of candidate.
//____________________________________________________
Int_t AliAnaParticleIsolation::GetMCIndex(Int_t mcTag)
{
  if(!IsDataMC()) return -1;
  
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
  // In case of several cone and thresholds analysis, open the cuts for the filling of the
  // track and cluster reference arrays in cone when done in the MakeAnalysisFillAOD().
  // The different cones, thresholds are tested for this list of tracks, clusters.
  if(fMakeSeveralIC)
  {
    AliInfo("Open default isolation cuts for multiple Isolation analysis");
    GetIsolationCut()->SetPtThreshold(100);
    GetIsolationCut()->SetPtFraction(100);
    GetIsolationCut()->SetConeSize(1);
  }
  
  if(!GetReader()->IsCTSSwitchedOn() && GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    AliFatal("STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!");
}

//____________________________________________
/// Initialize the parameters of the analysis
/// with default values.
//____________________________________________
void AliAnaParticleIsolation::InitParameters()
{
  SetInputAODName("PWG4Particle");
  SetAODObjArrayName("IsolationCone");
  AddToHistogramsName("AnaIsolation_");
  
  fIsoDetectorString = "EMCAL" ;
  fIsoDetector       = kEMCAL  ;
  
  fReMakeIC = kFALSE ;
  fMakeSeveralIC = kFALSE ;
  
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

  fDecayTagsM02Cut = 0.27;
  
  fNBkgBin = 11;
  fBkgBinLimit[ 0] = 00.0; fBkgBinLimit[ 1] = 00.2; fBkgBinLimit[ 2] = 00.3; fBkgBinLimit[ 3] = 00.4; fBkgBinLimit[ 4] = 00.5;
  fBkgBinLimit[ 5] = 01.0; fBkgBinLimit[ 6] = 01.5; fBkgBinLimit[ 7] = 02.0; fBkgBinLimit[ 8] = 03.0; fBkgBinLimit[ 9] = 05.0;
  fBkgBinLimit[10] = 10.0; fBkgBinLimit[11] = 100.;
  for(Int_t ibin = fNBkgBin+1; ibin < 20; ibin++) fBkgBinLimit[ibin] = 00.0;

  fNPtTrigBin = 7;
  fPtTrigBinLimit[ 0] =  8; fPtTrigBinLimit[ 1] = 10; fPtTrigBinLimit[ 2] = 12; fPtTrigBinLimit[ 3] = 16; fPtTrigBinLimit[ 4] = 20;
  fPtTrigBinLimit[ 5] = 25; fPtTrigBinLimit[ 6] = 30; fPtTrigBinLimit[ 7] = 50;
  for(Int_t ibin = fNPtTrigBin+1; ibin < 20; ibin++) fPtTrigBinLimit[ibin] = 10000.0;
  
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

  
  //----------- Several IC-----------------
  fNCones             = 5 ;
  fNPtThresFrac       = 5 ;
  fConeSizes      [0] = 0.1;     fConeSizes      [1] = 0.2;   fConeSizes      [2] = 0.3; fConeSizes      [3] = 0.4;  fConeSizes      [4] = 0.5;
  fPtThresholds   [0] = 1.;      fPtThresholds   [1] = 2.;    fPtThresholds   [2] = 3.;  fPtThresholds   [3] = 4.;   fPtThresholds   [4] = 5.;
  fPtFractions    [0] = 0.05;    fPtFractions    [1] = 0.075; fPtFractions    [2] = 0.1; fPtFractions    [3] = 1.25; fPtFractions    [4] = 1.5;
  fSumPtThresholds[0] = 1.;      fSumPtThresholds[1] = 2.;    fSumPtThresholds[2] = 3.;  fSumPtThresholds[3] = 4.;   fSumPtThresholds[4] = 5.;
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
  AliAODPWG4ParticleCorrelation* pLeading = 0;
  
  // Loop on stored AOD particles, find leading trigger on the selected list, with at least min pT.
  
  for(Int_t iaod = 0; iaod < GetInputAODBranch()->GetEntriesFast() ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* particle =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    particle->SetLeadingParticle(kFALSE); // set it later
    
    // Vertex cut in case of mixing
    if(GetMixedEvent())
    {
      Int_t check = CheckMixedEventVertex(particle->GetCaloLabel(0), particle->GetTrackLabel(0));
      if(check ==  0) continue;
      if(check == -1) return kFALSE; // not sure if it is correct.
    }
    
    //check if it is low pt trigger particle
    if((particle->Pt() < GetIsolationCut()->GetPtThreshold() ||
        particle->Pt() < GetIsolationCut()->GetSumPtThreshold()) &&
       !fMakeSeveralIC)
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
  if(!GetInputAODBranch())
    AliFatal(Form("No input particles in AOD with name branch < %s >, STOP",GetInputAODName().Data()));
  
  if(strcmp(GetInputAODBranch()->GetClass()->GetName(), "AliAODPWG4ParticleCorrelation"))
    AliFatal(Form("Wrong type of AOD object, change AOD class name in input AOD: It should be <AliAODPWG4ParticleCorrelation> and not <%s>",
                  GetInputAODBranch()->GetClass()->GetName()));
  
  Int_t n = 0, nfrac = 0;
  Bool_t  isolated  = kFALSE ;
  Float_t coneptsum = 0, coneptlead = 0;
  TObjArray * pl    = 0x0; ;
  
  //Select the calorimeter for candidate isolation with neutral particles
  if      (GetCalorimeter() == kPHOS )
    pl = GetPHOSClusters();
  else if (GetCalorimeter() == kEMCAL)
    pl = GetEMCALClusters();
  
  //Loop on AOD branch, filled previously in AliAnaPhoton, find leading particle to do isolation only with it
  Int_t idLeading = -1 ;
  Int_t iaod0 = 0;
  Int_t naod  = GetInputAODBranch()->GetEntriesFast();
  
  AliDebug(1,Form("Input aod branch entries %d", naod));
  
  if(IsLeadingOnlyOn())
  {
    Bool_t leading = IsTriggerTheNearSideEventLeadingParticle(idLeading);
    if(!leading)
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
    AliAODPWG4ParticleCorrelation * aodinput =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));

    // Check isolation only of clusters in fiducial region
    
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(aodinput->Eta(), aodinput->Phi(), aodinput->GetDetectorTag()) ;
      if(! in ) continue ;
    }
    
    //If too small or too large pt, skip
    Float_t pt = aodinput->Pt();
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ;

    
    //check if it is low pt trigger particle
    if( ( pt < GetIsolationCut()->GetPtThreshold() ||  pt < GetIsolationCut()->GetSumPtThreshold() ) &&
       !fMakeSeveralIC)
    {
      continue ; //trigger should not come from underlying event
    }
    
    //After cuts, study isolation
    n=0; nfrac = 0; isolated = kFALSE; coneptsum = 0; coneptlead = 0;
    GetIsolationCut()->MakeIsolationCut(GetCTSTracks(),pl,
                                        GetReader(), GetCaloPID(),
                                        kTRUE, aodinput, GetAODObjArrayName(),
                                        n,nfrac,coneptsum,coneptlead,isolated);
    
    if(!fMakeSeveralIC) aodinput->SetIsolated(isolated);
    
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
  if(IsDataMC() && !fMakeSeveralIC) FillAcceptanceHistograms();
  
  //Loop on stored AOD
  Int_t naod = GetInputAODBranch()->GetEntriesFast();
  
  AliDebug(1,Form("Histo aod branch entries %d", naod));
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4ParticleCorrelation* aod =  (AliAODPWG4ParticleCorrelation*) (GetInputAODBranch()->At(iaod));
    
    if(IsLeadingOnlyOn() && !aod->IsLeadingParticle()) continue; // Try to isolate only leading cluster or track
    
    // Check isolation only of clusters in fiducial region
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(aod->Eta(),aod->Phi(),aod->GetDetectorTag()) ;
      if(! in ) continue ;
    }
    
    Float_t pt         = aod->Pt();
    
    //If too small or too large pt, skip
    if(pt < GetMinPt() || pt > GetMaxPt() ) continue ;
    
    Int_t mcTag        = aod->GetTag() ;
    Int_t mcIndex      = GetMCIndex(mcTag);
    
    // --- In case of redoing isolation from delta AOD ----
    // Not standard case, not used since its implementation
    if(fMakeSeveralIC)
    {
      //Analysis of multiple IC at same time
      MakeSeveralICAnalysis(aod,mcIndex);
      continue;
    }
    
    // --- In case of redoing isolation multiple cuts ----
    
    if(fReMakeIC)
    {
      //In case a more strict IC is needed in the produced AOD
      Bool_t  isolated = kFALSE;
      Int_t   n = 0, nfrac = 0;
      Float_t coneptsum = 0, coneptlead = 0;
      
      // Recover reference arrays with clusters and tracks
      TObjArray * refclusters = aod->GetObjArray(GetAODObjArrayName()+"Clusters");
      TObjArray * reftracks   = aod->GetObjArray(GetAODObjArrayName()+"Tracks");
      
      GetIsolationCut()->MakeIsolationCut(reftracks,   refclusters,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, aod, "",
                                          n,nfrac,coneptsum,coneptlead,isolated);
    }
    
    Bool_t  isolated   = aod->IsIsolated();
    Float_t energy     = aod->E();
    Float_t phi        = GetPhi(aod->Phi());
    Float_t eta        = aod->Eta();
    Float_t m02        = aod->GetM02();
    Int_t   iSM        = aod->GetSModNumber();
    
    AliDebug(1,Form("pt %1.1f, eta %1.1f, phi %1.1f, Isolated %d",pt, eta, phi, isolated));
       
    //---------------------------------------------------------------
    // Recover original cluster if requested, needed for some studies
    //---------------------------------------------------------------
    if ( fFillOverlapHistograms     || fFillTMHisto || 
         fFillEMCALRegionHistograms || fStudyExoticTrigger || fStudyNCellsCut)
    {
      Int_t iclus = -1;
      fCluster = 0;
      fIsExoticTrigger = kFALSE;
      fClusterExoticity = 1;
      
      if     (GetCalorimeter() == kEMCAL) { fClustersArr = GetEMCALClusters(); fCaloCells = GetEMCALCells() ; }
      else if(GetCalorimeter() == kPHOS ) { fClustersArr = GetPHOSClusters (); fCaloCells = GetPHOSCells () ; }
      
      if(fClustersArr)
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
            
            if(fClusterExoticity > 0.97) fIsExoticTrigger = kTRUE ;
            //fIsExoticTrigger = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(absIdMax,GetEMCALCells(),bc);
            //if ( fIsExoticTrigger ) 
            //printf("Isolation: IsExotic? %d, E %f, ncells %d, exoticity %f\n", 
            //       fIsExoticTrigger, aod->E(), fCluster->GetNCells(), exoticity);
          }
          
          if(fStudyNCellsCut)
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
        } // cluster found
        
      } // cluster array found
    }
      
    //---------------------------------------------------------------
    // Fill pt/sum pT distribution of particles in cone or in UE band
    //---------------------------------------------------------------
    
    Float_t coneptLeadCluster= 0;
    Float_t coneptLeadTrack  = 0;
    Float_t coneptsumCluster = 0;
    Float_t coneptsumTrack   = 0;
    Float_t coneptsumCell    = 0;
    Float_t coneptsumSubEtaBand = 0;
    Float_t coneptsumSubPhiBand = 0;
    
    CalculateTrackSignalInCone   (aod,coneptsumTrack  , coneptLeadTrack  );
    CalculateCaloSignalInCone    (aod,coneptsumCluster, coneptLeadCluster);
    CalculateCaloCellSignalInCone(aod,coneptsumCell);
    
    if(GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kNeutralAndCharged)
    {
      fhConeSumPtClustervsTrack ->Fill(coneptsumCluster, coneptsumTrack , GetEventWeight());
      fhConePtLeadClustervsTrack->Fill(coneptLeadCluster,coneptLeadTrack, GetEventWeight());

      if(coneptsumTrack  > 0) fhConeSumPtClusterTrackFrac ->Fill(pt, coneptsumCluster /coneptsumTrack , GetEventWeight());
      if(coneptLeadTrack > 0) fhConePtLeadClusterTrackFrac->Fill(pt, coneptLeadCluster/coneptLeadTrack, GetEventWeight());

      if(fFillCellHistograms)
      {
        fhConeSumPtCellvsTrack        ->Fill(coneptsumCell, coneptsumTrack,               GetEventWeight());
        fhConeSumPtCellTrack          ->Fill(pt,            coneptsumTrack+coneptsumCell, GetEventWeight());
        fhConeSumPtCellTrackTrigEtaPhi->Fill(eta,      phi, coneptsumTrack+coneptsumCell *GetEventWeight()); // check
      }
    }
    
    fhConeSumPt              ->Fill(pt,       coneptsumTrack+coneptsumCluster, GetEventWeight());
    fhConeSumPtTrigEtaPhi    ->Fill(eta, phi, coneptsumTrack+coneptsumCluster *GetEventWeight()); // check
    if(fStudyExoticTrigger && fIsExoticTrigger)
      fhConeSumPtExoTrigger  ->Fill(pt, coneptsumTrack+coneptsumCluster, GetEventWeight());
    
    Float_t coneptLead = coneptLeadTrack;
    if(coneptLeadCluster > coneptLeadTrack) coneptLead = coneptLeadCluster;
    fhConePtLead->Fill(pt, coneptLead, GetEventWeight());
    
    AliDebug(1,Form("Particle %d Energy Sum in Isolation Cone %2.2f, Leading pT in cone %2.2f",
             iaod, coneptsumTrack+coneptsumCluster, coneptLead));
    
    //---------------------------------------------------------------
    // Normalize phi/eta band per area unit
    //---------------------------------------------------------------
    if(fFillUEBandSubtractHistograms)
      CalculateNormalizeUEBandPerUnitArea(aod, coneptsumCluster, coneptsumCell, coneptsumTrack, coneptsumSubEtaBand, coneptsumSubPhiBand, mcIndex) ;
        
    //---------------------------------------------------------------
    // EMCAL SM regions
    //---------------------------------------------------------------
    if(fFillEMCALRegionHistograms && GetCalorimeter() == kEMCAL) 
      StudyEMCALRegions(pt, phi, eta, m02, 
                        coneptsumTrack, coneptsumCluster, isolated, iSM);
        
    //---------------------------------------------------------------
    // Conversion radius in MC
    //---------------------------------------------------------------
    if(IsDataMC() && fStudyMCConversionRadius) 
      StudyMCConversionRadius(pt, iSM, isolated, m02, mcTag, aod->GetLabel());

    //---------------------------------------------------------------
    // Fill Shower shape and track matching histograms
    //---------------------------------------------------------------
    
    FillTrackMatchingShowerShapeControlHistograms(aod, coneptsumTrack+coneptsumCluster, 
                                                  coneptsumTrack, coneptsumCluster, coneptLead, mcIndex);
    
    //---------------------------------------------------------------
    // Isolated/ Non isolated histograms
    //---------------------------------------------------------------
    
    if(isolated)
    {
      AliDebug(1,Form("Particle %d ISOLATED: fill histograms", iaod));
      
      fhEIso      ->Fill(energy,      GetEventWeight());
      fhPtIso     ->Fill(pt    ,      GetEventWeight());
      fhPhiIso    ->Fill(pt    , phi, GetEventWeight());
      fhEtaIso    ->Fill(pt    , eta, GetEventWeight());
      fhEtaPhiIso ->Fill(eta   , phi, GetEventWeight());
      
      if(fStudyExoticTrigger && fIsExoticTrigger)
      {
        fhEIsoExoTrigger ->Fill(energy,  GetEventWeight());
        fhPtIsoExoTrigger->Fill(pt,      GetEventWeight());        
      }
      
      if(IsDataMC())
      {
        // For histograms in arrays, index in the array, corresponding to any particle origin
        if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
        {
          fhPtIsoMC [kmcPhoton]->Fill(pt,      GetEventWeight());
          fhPhiIsoMC[kmcPhoton]->Fill(pt, phi, GetEventWeight());
          fhEtaIsoMC[kmcPhoton]->Fill(pt, eta, GetEventWeight());
        }
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay )
          {
            fhPtIsoMC [kmcPi0DecayLostPair]->Fill(pt,      GetEventWeight());
            fhPhiIsoMC[kmcPi0DecayLostPair]->Fill(pt, phi, GetEventWeight());
            fhEtaIsoMC[kmcPi0DecayLostPair]->Fill(pt, eta, GetEventWeight());
          }
          else if( mcIndex == kmcEtaDecay )
          {
            fhPtIsoMC [kmcEtaDecayLostPair]->Fill(pt,      GetEventWeight());
            fhPhiIsoMC[kmcEtaDecayLostPair]->Fill(pt, phi, GetEventWeight());
            fhEtaIsoMC[kmcEtaDecayLostPair]->Fill(pt, eta, GetEventWeight());
          }
        }

        fhPtIsoMC [mcIndex]->Fill(pt,      GetEventWeight());
        fhPhiIsoMC[mcIndex]->Fill(pt, phi, GetEventWeight());
        fhEtaIsoMC[mcIndex]->Fill(pt, eta, GetEventWeight());
      }//Histograms with MC
      
      if(fFillNLMHistograms)
        fhPtNLocMaxIso ->Fill(pt, aod->GetNLM(), GetEventWeight()) ;
      
      if(IsHighMultiplicityAnalysisOn())
      {
        fhPtCentralityIso ->Fill(pt, GetEventCentrality(), GetEventWeight()) ;
        fhPtEventPlaneIso ->Fill(pt, GetEventPlaneAngle(), GetEventWeight()) ;
      }

      if(IsPileUpAnalysisOn())
      {
        if(GetReader()->IsPileUpFromSPD())               { fhEIsoPileUp[0] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[0]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromEMCal())             { fhEIsoPileUp[1] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[1]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDOrEMCal())        { fhEIsoPileUp[2] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[2]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDAndEMCal())       { fhEIsoPileUp[3] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[3]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())    { fhEIsoPileUp[4] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[4]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())    { fhEIsoPileUp[5] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[5]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) { fhEIsoPileUp[6] ->Fill(energy, GetEventWeight()) ; fhPtIsoPileUp[6]->Fill(pt, GetEventWeight()) ; }
        
        // Fill histograms to undertand pile-up before other cuts applied
        // Remember to relax time cuts in the reader
        FillPileUpHistograms(energy, aod->GetTime());//aod->GetCaloLabel(0));
      }

    }//Isolated histograms
    else // NON isolated
    {
      AliDebug(1,Form("Particle %d NOT ISOLATED, fill histograms", iaod));
      
      fhENoIso        ->Fill(energy,  GetEventWeight());
      fhPtNoIso       ->Fill(pt,      GetEventWeight());
      fhEtaPhiNoIso   ->Fill(eta,phi, GetEventWeight());
      
      if(fStudyExoticTrigger && fIsExoticTrigger)
      {
        fhENoIsoExoTrigger ->Fill(energy,  GetEventWeight());
        fhPtNoIsoExoTrigger->Fill(pt,      GetEventWeight());        
      }
      
      if(IsDataMC())
      {
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtNoIsoMC[kmcPhoton]->Fill(pt, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost)  )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtNoIsoMC[kmcPi0DecayLostPair]->Fill(pt, GetEventWeight());
          else if( mcIndex == kmcEtaDecay ) fhPtNoIsoMC[kmcEtaDecayLostPair]->Fill(pt, GetEventWeight());
        }
        
        fhPtNoIsoMC[mcIndex]->Fill(pt, GetEventWeight());
      }
      
      if(fFillNLMHistograms)
        fhPtNLocMaxNoIso ->Fill(pt, aod->GetNLM(), GetEventWeight());
      
      if(IsPileUpAnalysisOn())
      {
        if(GetReader()->IsPileUpFromSPD())                { fhENoIsoPileUp[0] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[0]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromEMCal())              { fhENoIsoPileUp[1] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[1]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDOrEMCal())         { fhENoIsoPileUp[2] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[2]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDAndEMCal())        { fhENoIsoPileUp[3] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[3]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())     { fhENoIsoPileUp[4] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[4]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())     { fhENoIsoPileUp[5] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[5]->Fill(pt, GetEventWeight()) ; }
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())  { fhENoIsoPileUp[6] ->Fill(energy, GetEventWeight()) ; fhPtNoIsoPileUp[6]->Fill(pt, GetEventWeight()) ; }
      }
    } // non iso
    
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
  
  //Double_t photonY   = -100 ;
  Double_t photonE   = -1 ;
  Double_t photonPt  = -1 ;
  Double_t photonPhi =  100 ;
  Double_t photonEta = -1 ;
  
  Int_t    pdg       =  0 ;
  Int_t    status    =  0 ;
  Int_t    tag       =  0 ;
  Int_t    mcIndex   =  0 ;
  Int_t    nprim     = GetMC()->GetNumberOfTracks();
  
  AliVParticle * primary = 0;
  
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

  for(Int_t i=0 ; i < nprim; i++)
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
      
      //photonY   = 0.5*TMath::Log((primary->E()+primary->Pz())/(primary->Energy()-prim->Pz())) ;
      
      //Photon kinematics
      primary->Momentum(fMomentum);
        
    // Select only photons in the final state
    if(pdg != 22  && pdg!=111 && pdg !=221) continue ;
    
    // Consider only final state particles, but this depends on generator,
    // status 1 is the usual one, in case of not being ok, leave the possibility
    // to not consider this.
    if( pdg == 22 && status != 1 &&
        GetMCAnalysisUtils()->GetMCGenerator() != AliMCAnalysisUtils::kBoxLike ) continue ;
    
    // If too small or too large pt, skip, same cut as for data analysis
    photonPt  = fMomentum.Pt () ;
    
    if(photonPt < GetMinPt() || photonPt > GetMaxPt() ) continue ;
    
    photonE   = fMomentum.E  () ;
    photonEta = fMomentum.Eta() ;
    photonPhi = fMomentum.Phi() ;
    
    if(photonPhi < 0) photonPhi+=TMath::TwoPi();
    
    // Check if photons hit the Calorimeter acceptance
    if(IsRealCaloAcceptanceOn() && fIsoDetector!=kCTS) // defined on base class
    {
      if ( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primary) ) continue ;
    }
    
    // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
    if(!GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),fIsoDetector)) continue ;
    
    // Get tag of this particle photon from fragmentation, decay, prompt ...
    // Set the origin of the photon.
    tag = GetMCAnalysisUtils()->CheckOrigin(i, GetMC());
    
    if(pdg == 22 && !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
    {
      // A conversion photon from a hadron, skip this kind of photon
      // printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!\n ");
      // GetMCAnalysisUtils()->PrintMCTag(tag);
      
      continue;
    }
    
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
      if( ndaugh!=2 || (d2Pdg != d1Pdg && d1Pdg!=22)) okpi0 = kFALSE;
      
      // Check overlap of decays
      if( okpi0 && fMakePrimaryPi0DecayStudy)
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
      mcIndex = kmcPrimPi0Decay;
    }
    else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) )
    {
      mcIndex = kmcPrimEtaDecay;
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
    Double_t sumPtInCone = 0; Double_t dR=0. ;
    
    AliVParticle * mcisop = 0;
    
    Int_t partInConeStatus = -1, partInConeMother = -1;
    Double_t partInConePt = 0, partInConeE = 0, partInConeEta = 0, partInConePhi = 0;
    Int_t partInConeCharge = 0, npart = 0, partInConePDG = 0;
    Bool_t physPrimary = kFALSE;
    
    for(Int_t ip = 0; ip < nprim ; ip++)
    {
      if(ip==i) continue;
      
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
      if( partInConeStatus  > 1 &&
         GetMCAnalysisUtils()->GetMCGenerator()!= AliMCAnalysisUtils::kBoxLike ) continue ;
      
      // Protection against floating point exception
      if ( mcisop->E() == TMath::Abs(mcisop->Pz())    ||
          (mcisop->E() - mcisop->Pz()) < 1e-3         ||
          mcisop->Pt() < 0.01  || mcisop->Eta() > 10  || 
          (mcisop->E() + mcisop->Pz()) < 0               )  continue ; 
      
      partInConeMother = mcisop->GetMother();
      partInConePt     = mcisop->Pt();
      partInConeE      = mcisop->E();
      partInConeEta    = mcisop->Eta();
      partInConePhi    = mcisop->Phi();
      partInConePDG    = mcisop->PdgCode();
      //partInConeCharge = TMath::Abs((Int_t) TDatabasePDG::Instance()->GetParticle(partInConePDG)->Charge());
      partInConeCharge = TMath::Abs(mcisop->Charge());
      physPrimary      = mcisop->IsPhysicalPrimary(); 
      mcisop->Momentum(fMomIso);
      
      if( partInConeMother == i ) continue;
      
      //
      // Apply acceptance and energy/pt cut for particles in cone
      if(fSelectPrimariesInCone)
      {
        if( partInConeCharge > 0) // charged pT cut and acceptance
        {
          if( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kOnlyNeutral ) continue;
          
          if( partInConePt < GetReader()->GetCTSPtMin () ) continue;
          
          if(!GetReader()->GetFiducialCut()->IsInFiducialCut(fMomIso.Eta(),fMomIso.Phi(),kCTS)) continue ;
        }
        else // neutrals E cut and acceptance
        {
          if( GetIsolationCut()->GetParticleTypeInCone() == AliIsolationCut::kOnlyCharged ) continue;
          
          if( partInConeE  <= minECalo ) continue;
          
          if(!GetReader()->GetFiducialCut()->IsInFiducialCut(fMomIso.Eta(),fMomIso.Phi(),GetCalorimeter())) continue ;
          
          if(IsRealCaloAcceptanceOn()) // defined on base class
          {
            if ( !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), mcisop) ) continue ;
          }
        }
      }
      
      //printf("status %d, E %f, pt %f, pz %f, pdg %d, gamma eta %f, phi %f, par in cone eta %f phi %f\n",
      //       partInConeStatus, partInConeE, partInConePt, mcisopAOD->Pz(), partInConePDG, photonEta, photonPhi, partInConeEta, partInConePhi);
      
      dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);
      
      if(dR > GetIsolationCut()->GetConeSize())
        continue;
      
      if ( partInConeCharge > 0 &&  TMath::Abs(partInConePDG) != 11 ) // exclude electrons and neutrals
      {
        Int_t mcChTag = 3;
        if      ( TMath::Abs(partInConePDG) == 211  )  mcChTag = 0;
        else if ( TMath::Abs(partInConePDG) == 321  )  mcChTag = 1; 
        else if ( TMath::Abs(partInConePDG) == 2212 )  mcChTag = 2; 
        
        if(physPrimary)
          fhPtTrackInConeMCPrimaryGener  [mcChTag]->Fill(photonPt, partInConePt, GetEventWeight());
        else
          fhPtTrackInConeMCSecondaryGener[mcChTag]->Fill(photonPt, partInConePt, GetEventWeight());

        //printf("Selected particles pdg %d, status %d\n", partInConePDG, partInConeStatus);
      }
      

      if ( !physPrimary ) continue ; 
          
      sumPtInCone += partInConePt;
      if(partInConePt > GetIsolationCut()->GetPtThreshold() &&
         partInConePt < GetIsolationCut()->GetPtThresholdMax()) npart++;
    }
    
    ///////END ISO MC/////////////////////////
    
    // Fill the histograms, only those in the defined calorimeter acceptance
    
    fhEtaPrimMC[kmcPrimPhoton]->Fill(photonPt, photonEta, GetEventWeight()) ;
    fhPhiPrimMC[kmcPrimPhoton]->Fill(photonPt, photonPhi, GetEventWeight()) ;
    fhEPrimMC  [kmcPrimPhoton]->Fill(photonE ,            GetEventWeight()) ;
    fhPtPrimMC [kmcPrimPhoton]->Fill(photonPt,            GetEventWeight()) ;
    
    fhEtaPrimMC[mcIndex]->Fill(photonPt, photonEta, GetEventWeight()) ;
    fhPhiPrimMC[mcIndex]->Fill(photonPt, photonPhi, GetEventWeight()) ;
    fhEPrimMC  [mcIndex]->Fill(photonE ,            GetEventWeight()) ;
    fhPtPrimMC [mcIndex]->Fill(photonPt,            GetEventWeight()) ;
    
    // In case the photon is a decay from pi0 or eta,
    // study how the decay kinematics affects the isolation
    Int_t  ndaugh   = -1;
    Bool_t okpi0    =  0, ok1     =  0, ok2     =  0;
    Int_t  pi0label = -1, d1Label = -1, d2Label = -1;
    Bool_t d2Acc   = kTRUE, overlap = kTRUE;
    Int_t  d2AbsId = -1;
    Float_t dRdaugh2 = 0, d12Angle = 0;
    
    if(fMakePrimaryPi0DecayStudy)
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
          if(dRdaugh2 > GetIsolationCut()->GetConeSize())
            fhPtPrimMCPi0DecayPairOutOfCone->Fill(photonPt, GetEventWeight());
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCPi0DecayPairOutOfAcceptance->Fill(photonPt, GetEventWeight());
            if(!overlap) fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight());
          }
          
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCPi0DecayPairNoOverlap->Fill(photonPt, GetEventWeight());
          
          // Second decay pt smaller than threshold
          if(d2Acc && dRdaugh2 < GetIsolationCut()->GetConeSize() &&
             fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCPi0DecayPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight());
            if(!overlap)
            {
              fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight());
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight());
            }
          }
        } // pi0 decay
        else // eta decay
        {
          // Second decay out of cone
          if(dRdaugh2 > GetIsolationCut()->GetConeSize())
            fhPtPrimMCEtaDecayPairOutOfCone->Fill(photonPt, GetEventWeight());
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCEtaDecayPairOutOfAcceptance->Fill(photonPt, GetEventWeight());
            if(!overlap) fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight());
          }
          
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCEtaDecayPairNoOverlap->Fill(photonPt, GetEventWeight());
          
          // Second decay pt smaller than threshold
          if(d2Acc && dRdaugh2 < GetIsolationCut()->GetConeSize() &&
             fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCEtaDecayPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight());
            if(!overlap)
            {
              fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight());
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight());
            }
          }
        } // eta decay
        
      } // eta or pi0 decay

      if(overlapPi0)
      {
        if( mcIndex == kmcPrimPi0) fhPtPrimMCPi0Overlap->Fill(photonPt, GetEventWeight());
        if( mcIndex == kmcPrimEta) fhPtPrimMCEtaOverlap->Fill(photonPt, GetEventWeight());
      }
    }
    
    // Isolated?
    Bool_t isolated = kFALSE;
    if(GetIsolationCut()->GetICMethod() == AliIsolationCut::kSumPtIC   &&
       (sumPtInCone < GetIsolationCut()->GetSumPtThreshold() ||
        sumPtInCone > GetIsolationCut()->GetSumPtThresholdMax()))
      isolated = kTRUE;
    
    if(GetIsolationCut()->GetICMethod() == AliIsolationCut::kPtThresIC &&
       npart == 0)
      isolated = kTRUE;
    
    if(isolated)
    {
      fhPtPrimMCiso [mcIndex]      ->Fill(photonPt, GetEventWeight()) ;
      fhPtPrimMCiso [kmcPrimPhoton]->Fill(photonPt, GetEventWeight()) ;
      
      if(fMakePrimaryPi0DecayStudy)
      {
        if( mcIndex == kmcPrimPi0Decay )
        {
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCPi0DecayIsoPairNoOverlap->Fill(photonPt, GetEventWeight());
          
          // Second decay out of cone
          if(dRdaugh2 > GetIsolationCut()->GetConeSize())
            fhPtPrimMCPi0DecayIsoPairOutOfCone->Fill(photonPt, GetEventWeight());
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCPi0DecayIsoPairOutOfAcceptance->Fill(photonPt, GetEventWeight());
            if(!overlap) fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight());
          }
          
          // Second decay pt smaller than threshold
          if(d2Acc && dRdaugh2 < GetIsolationCut()->GetConeSize() &&
             fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight());
            if(!overlap)
            {
              fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight());
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight());
            }
          }
        }// pi0 decay
        else if( mcIndex == kmcPrimEtaDecay )
        {
          // Not Overlapped decay
          if(!overlap) fhPtPrimMCEtaDecayIsoPairNoOverlap->Fill(photonPt, GetEventWeight());
          
          // Second decay out of cone
          if(dRdaugh2 > GetIsolationCut()->GetConeSize())
            fhPtPrimMCEtaDecayIsoPairOutOfCone->Fill(photonPt, GetEventWeight());
          
          // Second decay out of acceptance
          if(!ok2 || !d2Acc || fMomDaugh2.E() <= minECalo)
          {
            fhPtPrimMCEtaDecayIsoPairOutOfAcceptance->Fill(photonPt, GetEventWeight());
            if(!overlap) fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap->Fill(photonPt, GetEventWeight());
          }
          
          // Second decay pt smaller than threshold
          if(d2Acc && dRdaugh2 < GetIsolationCut()->GetConeSize() &&
             fMomDaugh2.E() < GetIsolationCut()->GetPtThreshold())
          {
            fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt->Fill(photonPt, GetEventWeight());
            if(!overlap)
            {
              fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap->Fill(photonPt, GetEventWeight());
              if(fMomDaugh2.E() > minECalo) fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE->Fill(photonPt, GetEventWeight());
            }
          }
        }// eta decay
        
        if(overlapPi0)
        {
          if( mcIndex == kmcPrimPi0 ) fhPtPrimMCPi0IsoOverlap->Fill(photonPt, GetEventWeight());
          if( mcIndex == kmcPrimEta ) fhPtPrimMCEtaIsoOverlap->Fill(photonPt, GetEventWeight());
        }
      }
    } // isolated

  }//loop on primaries
  
  AliDebug(1,"End");
  
}

//_____________________________________________________________________________________
/// Isolation Cut Analysis for both methods and different pt cuts and cones.
//_____________________________________________________________________________________
void  AliAnaParticleIsolation::MakeSeveralICAnalysis(AliAODPWG4ParticleCorrelation* ph,
                                                     Int_t mcIndex)
{
  Float_t ptC   = ph->Pt();
  Float_t etaC  = ph->Eta();
  Float_t phiC  = ph->Phi();
  if(phiC<0) phiC += TMath::TwoPi();
  Int_t   tag   = ph->GetTag();

  Int_t   decayTag = 0;
  if(fFillTaggedDecayHistograms)
  {
    decayTag = ph->DecayTag();
    if(decayTag < 0) decayTag = 0; 
  }

  AliDebug(1,Form("Isolate pT %2.2f, decay tag %d",ptC, decayTag));
  
  //Keep original setting used when filling AODs, reset at end of analysis
  Float_t ptthresorg = GetIsolationCut()->GetPtThreshold();
  Float_t ptfracorg  = GetIsolationCut()->GetPtFraction();
  Float_t ptsumcorg  = GetIsolationCut()->GetSumPtThreshold();
  Float_t rorg       = GetIsolationCut()->GetConeSize();
  
  Float_t coneptsum = 0, coneptlead = 0;
  Int_t   n    [10][10];//[fNCones][fNPtThresFrac];
  Int_t   nfrac[10][10];//[fNCones][fNPtThresFrac];
  Bool_t  isolated  = kFALSE;
  
  // Fill hist with all particles before isolation criteria
  fhENoIso     ->Fill(ph->E(),   GetEventWeight());
  fhPtNoIso    ->Fill(ptC,       GetEventWeight());
  fhEtaPhiNoIso->Fill(etaC,phiC, GetEventWeight());
  
  if(IsDataMC())
  {
    if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
      fhPtNoIsoMC[kmcPhoton]->Fill(ptC, GetEventWeight());
    
    if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
    {
      if     ( mcIndex==kmcPi0Decay ) fhPtNoIsoMC[kmcPi0DecayLostPair]->Fill(ptC, GetEventWeight());
      else if( mcIndex==kmcEtaDecay ) fhPtNoIsoMC[kmcEtaDecayLostPair]->Fill(ptC, GetEventWeight());
    }
    
    fhPtNoIsoMC[mcIndex]->Fill(ptC, GetEventWeight());
  }
  
  // Candidates tagged as decay in another analysis (AliAnaPi0EbE)
  if(fFillTaggedDecayHistograms && decayTag > 0)
  {
    for(Int_t ibit = 0; ibit < fNDecayBits; ibit++)
    {
      if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[ibit]))
      {
        fhPtDecay    [0][ibit]->Fill(ptC),      GetEventWeight();
        fhEtaPhiDecay[0][ibit]->Fill(etaC,phiC, GetEventWeight());
        
        if(IsDataMC())
        {
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhPtDecayMC[0][ibit][kmcPhoton]->Fill(ptC, GetEventWeight());

          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
           if      (mcIndex==kmcPi0Decay) fhPtDecayMC[0][ibit][kmcPi0DecayLostPair]->Fill(ptC, GetEventWeight());
           else if (mcIndex==kmcEtaDecay) fhPtDecayMC[0][ibit][kmcEtaDecayLostPair]->Fill(ptC, GetEventWeight());
          }
          
          fhPtDecayMC[0][ibit][mcIndex]->Fill(ptC, GetEventWeight());
        }
      } // bit ok
    } // bit loop
  } // decay histograms
  
  // Get vertex for photon momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  // Loop on cone sizes
  for(Int_t icone = 0; icone<fNCones; icone++)
  {
    // Recover reference arrays with clusters and tracks
    TObjArray * refclusters = ph->GetObjArray(GetAODObjArrayName()+"Clusters");
    TObjArray * reftracks   = ph->GetObjArray(GetAODObjArrayName()+"Tracks");
    
    //If too small or too large pt, skip
    if(ptC < GetMinPt() || ptC > GetMaxPt() ) continue ;
    
    //In case a more strict IC is needed in the produced AOD
    
    isolated = kFALSE; coneptsum = 0; coneptlead = 0;
    
    GetIsolationCut()->SetSumPtThreshold(100);
    GetIsolationCut()->SetPtThreshold(100);
    GetIsolationCut()->SetPtFraction(100);
    GetIsolationCut()->SetConeSize(fConeSizes[icone]);
    
    // Retreive pt tracks to fill histo vs. pt leading
    //Fill pt distribution of particles in cone
    //fhPtLeadingPt(),fhPerpSumPtLeadingPt(),fhPerpPtLeadingPt(),
    
    // Tracks in perpendicular cones
    Double_t sumptPerp = 0. ;
    TObjArray * trackList   = GetCTSTracks() ;
    for(Int_t itrack=0; itrack < trackList->GetEntriesFast(); itrack++)
    {
      AliVTrack* track = (AliVTrack *) trackList->At(itrack);
      //fill the histograms at forward range
      if(!track)
      {
        AliDebug(1,"Track not available?");
        continue;
      }
      
      Double_t dPhi = phiC - track->Phi() + TMath::PiOver2();
      Double_t dEta = etaC - track->Eta();
      Double_t arg  = dPhi*dPhi + dEta*dEta;
      Double_t pTrack = TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py());
        
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhPerpPtLeadingPt[icone]->Fill(ptC, pTrack, GetEventWeight());
        sumptPerp+=track->Pt();
      }
      
      dPhi = phiC - track->Phi() - TMath::PiOver2();
      arg  = dPhi*dPhi + dEta*dEta;
      if(TMath::Sqrt(arg) < fConeSizes[icone])
      {
        fhPerpPtLeadingPt[icone]->Fill(ptC, pTrack, GetEventWeight());
        sumptPerp+=track->Pt();
      }
    }
    
    fhPerpSumPtLeadingPt[icone]->Fill(ptC, sumptPerp, GetEventWeight());
    
    // Tracks in isolation cone, pT distribution and sum
    if(reftracks && GetIsolationCut()->GetParticleTypeInCone()!= AliIsolationCut::kOnlyNeutral)
    {
      for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
      {
        AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
        
        Float_t rad = GetIsolationCut()->Radius(etaC, phiC, track->Eta(), track->Phi());
        
        if(rad > fConeSizes[icone]) continue ;
        
        fhPtLeadingPt[icone]->Fill(ptC, track->Pt(), GetEventWeight());
        coneptsum += track->Pt();
      }
    }
    
    // Clusters in isolation cone, pT distribution and sum
    if(refclusters && GetIsolationCut()->GetParticleTypeInCone()!= AliIsolationCut::kOnlyCharged)
    {
      for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
      {
        AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
        
        calo->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line
        
        Float_t rad = GetIsolationCut()->Radius(etaC, phiC, fMomentum.Eta(), fMomentum.Phi());
        
        if(rad > fConeSizes[icone]) continue ;
        
        fhPtLeadingPt[icone]->Fill(ptC, fMomentum.Pt(), GetEventWeight());
        coneptsum += fMomentum.Pt();
      }
    }
    
    fhSumPtLeadingPt[icone]->Fill(ptC, coneptsum, GetEventWeight());
    
    if(IsDataMC())
    {
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
        fhSumPtLeadingPtMC[kmcPhoton][icone]->Fill(ptC, coneptsum, GetEventWeight()) ;
      
      if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
      {
        if      ( mcIndex==kmcPi0Decay ) fhSumPtLeadingPtMC[kmcPi0DecayLostPair][icone]->Fill(ptC, coneptsum, GetEventWeight()) ;
        else if ( mcIndex==kmcEtaDecay ) fhSumPtLeadingPtMC[kmcEtaDecayLostPair][icone]->Fill(ptC, coneptsum, GetEventWeight()) ;
      }
      
      fhSumPtLeadingPtMC[mcIndex][icone]->Fill(ptC, coneptsum, GetEventWeight()) ;
    }
    
    ///////////////////
    
    //Loop on pt thresholds
    for(Int_t ipt = 0; ipt < fNPtThresFrac ; ipt++)
    {
      n    [icone][ipt]=0;
      nfrac[icone][ipt]=0;
      GetIsolationCut()->SetPtThreshold(fPtThresholds[ipt]);
      GetIsolationCut()->SetPtFraction(fPtFractions[ipt]) ;
      GetIsolationCut()->SetSumPtThreshold(fSumPtThresholds[ipt]);
      
      GetIsolationCut()->MakeIsolationCut(reftracks, refclusters,
                                          GetReader(), GetCaloPID(),
                                          kFALSE, ph, "",
                                          n[icone][ipt],nfrac[icone][ipt],
                                          coneptsum, coneptlead, isolated);
      
      // Normal pT threshold cut
      
      AliDebug(1,Form("Cone size %1.1f, ptThres  %1.1f, sumptThresh  %1.1f",fConeSizes[icone],fPtThresholds[ipt],fSumPtThresholds[ipt]));
      AliDebug(1,Form("\t n %d, nfrac %d, coneptsum %2.2f",n[icone][ipt],nfrac[icone][ipt],coneptsum));
      AliDebug(1,Form("pt %1.1f, eta %1.1f, phi %1.1f",ptC, etaC, phiC));
      
      if(n[icone][ipt] == 0)
      {
        AliDebug(1,"Filling pt threshold loop");
        
        fhPtThresIsolated [icone][ipt]->Fill(ptC,        GetEventWeight());
        fhEtaPhiPtThresIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtThresDecayIso    [icone][ipt]->Fill(ptC,        GetEventWeight());
            fhEtaPhiPtThresDecayIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtThresIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     (mcIndex == kmcPi0Decay) fhPtThresIsolatedMC[kmcPi0DecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
            else if(mcIndex == kmcEtaDecay) fhPtThresIsolatedMC[kmcEtaDecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          }
          
          fhPtThresIsolatedMC[mcIndex][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          
        }
      }
      
      // pt in cone fraction
      if(nfrac[icone][ipt] == 0)
      {
        AliDebug(1,"Filling frac loop");
        
        fhPtFracIsolated [icone][ipt]->Fill(ptC ,       GetEventWeight());
        fhEtaPhiPtFracIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtFracDecayIso    [icone][ipt]->Fill(ptC ,       GetEventWeight());
            fhEtaPhiPtFracDecayIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhPtFracIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhPtFracIsolatedMC[kmcPi0DecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
            else if( mcIndex == kmcEtaDecay ) fhPtFracIsolatedMC[kmcEtaDecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          }
          
          fhPtFracIsolatedMC[mcIndex][icone][ipt]->Fill(ptC, GetEventWeight()) ;
        }
      }
      
      AliDebug(1,Form("Checking IC method : %i",GetIsolationCut()->GetICMethod()));
      
      //Pt threshold on pt cand/ sum in cone histograms
      if(coneptsum < fSumPtThresholds[ipt])
      {
        AliDebug(1,"Filling sum loop");
        
        fhSumPtIsolated [icone][ipt]->Fill(ptC,        GetEventWeight()) ;
        fhEtaPhiPtSumIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight()) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtPtSumDecayIso[icone][ipt]->Fill(ptC,            GetEventWeight());
            fhEtaPhiPtSumDecayIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight()) ;
          }
        }
        
        if(IsDataMC())
        {
          if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
            fhSumPtIsolatedMC[kmcPhoton][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhSumPtIsolatedMC[kmcPi0DecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
            else if( mcIndex == kmcEtaDecay ) fhSumPtIsolatedMC[kmcEtaDecayLostPair][icone][ipt]->Fill(ptC, GetEventWeight()) ;
          }
          
          fhSumPtIsolatedMC[mcIndex][icone][ipt]->Fill(ptC, GetEventWeight()) ;
        }
      }
      
      // pt sum pt frac method
      //    if( ((fPtFractions[ipt]*ptC < fSumPtThresholds[ipt]) && (coneptsum < fSumPtThresholds[ipt])) || ((fPtFractions[ipt]*ptC > fSumPtThresholds[ipt]) && (coneptsum < fPtFractions[ipt]*ptC)) )
      
      if(coneptsum < fPtFractions[ipt]*ptC)
      {
        AliDebug(1,"Filling PtFrac PtSum loop");
        
        fhPtFracPtSumIso    [icone][ipt]->Fill(ptC,        GetEventWeight()) ;
        fhEtaPhiFracPtSumIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight()) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtFracPtSumDecayIso    [icone][ipt]->Fill(ptC,        GetEventWeight());
            fhEtaPhiFracPtSumDecayIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
          }
        }
      }
      
      // density method
      Float_t cellDensity = GetIsolationCut()->GetCellDensity( ph, GetReader());
      if(coneptsum < fSumPtThresholds[ipt]*cellDensity)
      {
        AliDebug(1,"Filling density loop");
        
        fhPtSumDensityIso    [icone][ipt]->Fill(ptC ,       GetEventWeight()) ;
        fhEtaPhiSumDensityIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight()) ;
        
        if( fFillTaggedDecayHistograms && decayTag > 0 && fNDecayBits > 0)
        {
          if(GetNeutralMesonSelection()->CheckDecayBit(decayTag,fDecayBits[0]))
          {
            fhPtSumDensityDecayIso    [icone][ipt]->Fill(ptC ,       GetEventWeight());
            fhEtaPhiSumDensityDecayIso[icone][ipt]->Fill(etaC, phiC, GetEventWeight());
          }
        }
      }
    }//pt thresh loop
    
    
  }//cone size loop
  
  // Reset original parameters for AOD analysis
  GetIsolationCut()->SetPtThreshold(ptthresorg);
  GetIsolationCut()->SetPtFraction(ptfracorg);
  GetIsolationCut()->SetSumPtThreshold(ptsumcorg);
  GetIsolationCut()->SetConeSize(rorg);
}

//_____________________________________________________________
/// Print some relevant parameters set for the analysis.
//_____________________________________________________________
void AliAnaParticleIsolation::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("ReMake Isolation          = %d \n",  fReMakeIC) ;
  printf("Make Several Isolation    = %d \n",  fMakeSeveralIC) ;
  printf("Calorimeter for isolation = %s \n",  GetCalorimeterString().Data()) ;
  printf("Detector for candidate isolation = %s \n", fIsoDetectorString.Data()) ;
  printf("Subtract UE from cone sum pT histo fill %d \n",fFillUEBandSubtractHistograms) ;
  
  if(fMakeSeveralIC)
  {
    printf("N Cone Sizes       =     %d\n", fNCones) ;
    printf("Cone Sizes          =    \n") ;
    for(Int_t i = 0; i < fNCones; i++)
      printf("  %1.2f;",  fConeSizes[i]) ;
    printf("    \n") ;
    
    printf("N pT thresholds/fractions = %d\n", fNPtThresFrac) ;
    printf(" pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtThresholds[i]) ;
    
    printf("    \n") ;
    
    printf(" pT fractions         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fPtFractions[i]) ;
    
    printf("    \n") ;
    
    printf("sum pT thresholds         =    \n") ;
    for(Int_t i = 0; i < fNPtThresFrac; i++)
      printf("   %2.2f;",  fSumPtThresholds[i]) ;
  }
  
  printf("    \n") ;
} 

//_____________________________________________________________
/// Set the detector for the analysis.
//_____________________________________________________________
void AliAnaParticleIsolation::SetTriggerDetector(TString & det)
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

