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
#include <TClass.h>
#include <TH2F.h>
#include "TParticle.h"
#include "TDatabasePDG.h"

// --- Analysis system ---
#include "AliAnaParticleIsolation.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliIsolationCut.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODMCParticle.h"
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
fFillTMHisto(0),                  fFillSSHisto(1),                          fFillEMCALRegionSSHistograms(0),
fFillUEBandSubtractHistograms(1), fFillCellHistograms(0),
fFillOverlapHistograms(0),
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
fMomentum(),                      fMomIso(),
fMomDaugh1(),                     fMomDaugh2(),
fTrackVector(),                   fProdVertex(),
fCluster(0),                      fClustersArr(0),
// Histograms
fhEIso(0),                        fhPtIso(0),
fhPtCentralityIso(0),             fhPtEventPlaneIso(0),
fhPtNLocMaxIso(0),
fhPhiIso(0),                      fhEtaIso(0),                              fhEtaPhiIso(0),
fhEtaPhiNoIso(0),
fhENoIso(0),                      fhPtNoIso(0),                             fhPtNLocMaxNoIso(0),
fhPtInCone(0),
fhPtClusterInCone(0),             fhPtCellInCone(0),                        fhPtTrackInCone(0),
fhPtTrackInConeOtherBC(0),        fhPtTrackInConeOtherBCPileUpSPD(0),
fhPtTrackInConeBC0(0),            fhPtTrackInConeVtxBC0(0),
fhPtTrackInConeBC0PileUpSPD(0),
fhPtInConePileUp(),               fhPtInConeCent(0),
fhPerpConeSumPt(0),               fhPtInPerpCone(0),
fhEtaPhiInConeCluster(0),         fhEtaPhiCluster(0),
fhEtaPhiInConeTrack(0),           fhEtaPhiTrack(0),
fhEtaBandCluster(0),              fhPhiBandCluster(0),
fhEtaBandTrack(0),                fhPhiBandTrack(0),
fhEtaBandCell(0),                 fhPhiBandCell(0),
fhConePtLead(0),                  fhConePtLeadCluster(0),                   fhConePtLeadTrack(0),
fhConePtLeadClustervsTrack(0),    fhConePtLeadClusterTrackFrac(0),
fhConeSumPt(0),                   fhConeSumPtCellTrack(0),
fhConeSumPtCell(0),               fhConeSumPtCluster(0),                    fhConeSumPtTrack(0),
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
fhConeSumPtEtaUESubCluster(0),              fhConeSumPtPhiUESubCluster(0),
fhConeSumPtEtaUESubClusterTrigEtaPhi(0),    fhConeSumPtPhiUESubClusterTrigEtaPhi(0),
fhConeSumPtEtaUESubCell(0),                 fhConeSumPtPhiUESubCell(0),
fhConeSumPtEtaUESubCellTrigEtaPhi(0),       fhConeSumPtPhiUESubCellTrigEtaPhi(0),
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
fhPtLeadConeBinMC(0),                       fhSumPtConeBinMC(0),
fhPtLeadConeBinDecay(0),                    fhSumPtConeBinDecay(0),
fhPtLeadConeBinLambda0(0),                  fhSumPtConeBinLambda0(0),
fhPtLeadConeBinLambda0MC(0),                fhSumPtConeBinLambda0MC(0),
fhPtTrigBinPtLeadCone(0),                   fhPtTrigBinSumPtCone(0),
fhPtTrigBinPtLeadConeMC(0),                 fhPtTrigBinSumPtConeMC(0),
fhPtTrigBinPtLeadConeDecay(0),              fhPtTrigBinSumPtConeDecay(0),
fhPtTrigBinLambda0vsPtLeadCone(0),          fhPtTrigBinLambda0vsSumPtCone(0),
fhPtTrigBinLambda0vsPtLeadConeMC(0),        fhPtTrigBinLambda0vsSumPtConeMC(0),
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
fhTimePileUpMainVertexZDistance(0), fhTimePileUpMainVertexZDiamond(0)
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
      fhPtLambda0MCWithNOverlap    [imc][i] = 0;
      fhPtLambda0MCConvWithNOverlap[imc][i] = 0;
      fhPtNOverlap                 [imc][i] = 0;
      fhPtNOverlapConv             [imc][i] = 0;
      fhTrackMatchedDEtaMC         [imc][i] = 0 ;
      fhTrackMatchedDPhiMC         [imc][i] = 0 ;
      fhTrackMatchedDEtaDPhiMC     [imc][i] = 0 ;
    }
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
          fhLam0EMCALRegionPerSM[i][ieta][iphi][ism] = 0; 
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
    fhEtaPhiCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
    if(rad < conesize)
    {
    	// Histogrqmw for all clusters in cone
      fhEtaPhiInConeCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
      continue ;
    }
      
    // Fill histogram for UE in phi band in EMCal acceptance
    if(fMomentum.Eta() > (etaTrig-conesize) && fMomentum.Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=fMomentum.Pt();
      fhPhiBandCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
      
    }
    
    // Fill histogram for UE in eta band in EMCal acceptance
    if(fMomentum.Phi() > (phiTrig-conesize) && fMomentum.Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=fMomentum.Pt();
      fhEtaBandCluster->Fill(fMomentum.Eta(), fMomentum.Phi(), GetEventWeight());
    }
  }
  
  fhConeSumPtEtaBandUECluster          ->Fill(ptTrig  ,         etaBandPtSum, GetEventWeight());
  fhConeSumPtPhiBandUECluster          ->Fill(ptTrig  ,         phiBandPtSum, GetEventWeight());
  fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSum *GetEventWeight()); // Check
  fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSum *GetEventWeight()); // Check
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
  
  Double_t sumptPerp= 0. ;
  Float_t ptTrig    = pCandidate->Pt() ;
  Float_t phiTrig   = pCandidate->Phi();
  Float_t etaTrig   = pCandidate->Eta();
  
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
      Int_t  trackID   = TMath::Abs(track->GetID()) ;
      Bool_t contained = kFALSE;
      
      for(Int_t i = 0; i < 4; i++) 
      {
        if( trackID == pCandidate->GetTrackLabel(i) ) contained = kTRUE;
      }
      
      if ( contained ) continue ;
    }
    
    // Histogram of eta:phi for all tracks
    fhEtaPhiTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
    
    //exclude particles in cone
    Float_t rad = GetIsolationCut()->Radius(etaTrig, phiTrig, track->Eta(), track->Phi());
    if(rad < conesize)
    {
      // Histogram of eta:phi for all tracks in cone
      fhEtaPhiInConeTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
      continue ;
    }
    
    // Fill histogram for UE in phi band
    if(track->Eta() > (etaTrig-conesize) && track->Eta()  < (etaTrig+conesize))
    {
      phiBandPtSum+=track->Pt();
      fhPhiBandTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
    }
    
    // Fill histogram for UE in eta band in EMCal acceptance
    if(track->Phi() > (phiTrig-conesize) && track->Phi() < (phiTrig+conesize))
    {
      etaBandPtSum+=track->Pt();
      fhEtaBandTrack->Fill(track->Eta(), track->Phi(), GetEventWeight());
    }
     
    // Fill the histograms at +-45 degrees in phi from trigger particle, perpedicular to trigger axis in phi
    Double_t pTrack = TMath::Sqrt(track->Px()*track->Px()+track->Py()*track->Py());
    Double_t dPhi   = phiTrig - track->Phi() + TMath::PiOver2();
    Double_t dEta   = etaTrig - track->Eta();
    Double_t arg    = dPhi*dPhi + dEta*dEta;
 
    if(TMath::Sqrt(arg) < conesize)
    {
      fhPtInPerpCone->Fill(ptTrig, pTrack, GetEventWeight());
      sumptPerp+=track->Pt();
    }
    
    dPhi = phiTrig - track->Phi() - TMath::PiOver2();
    arg  = dPhi*dPhi + dEta*dEta;
                                  
    if(TMath::Sqrt(arg) < conesize)
    {
      fhPtInPerpCone->Fill(ptTrig, pTrack, GetEventWeight());
      sumptPerp+=track->Pt();
    }
  }
  
  fhPerpConeSumPt                    ->Fill(ptTrig ,          sumptPerp   , GetEventWeight());
  fhConeSumPtEtaBandUETrack          ->Fill(ptTrig ,          etaBandPtSum, GetEventWeight());
  fhConeSumPtPhiBandUETrack          ->Fill(ptTrig ,          phiBandPtSum, GetEventWeight());
  fhConeSumPtEtaBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSum *GetEventWeight()); // check
  fhConeSumPtPhiBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSum *GetEventWeight()); // check
}

//_____________________________________________________________________________________________________________________________________
/// Normalize phi/eta band per area unit
//_____________________________________________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate, Float_t coneptsumCluster,
                                                                  Float_t coneptsumCell,          Float_t coneptsumTrack,
                                                                  Float_t &etaBandptsumTrackNorm, Float_t &etaBandptsumClusterNorm)
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
  
  
  // ------ //
  // Tracks //
  // ------ //
  Float_t phiUEptsumTrackNorm  = 0 ;
  Float_t etaUEptsumTrackNorm  = 0 ;
  Float_t coneptsumTrackSubPhi = 0 ;
  Float_t coneptsumTrackSubEta = 0 ;
  Float_t coneptsumTrackSubPhiNorm = 0 ;
  Float_t coneptsumTrackSubEtaNorm = 0 ;
  etaBandptsumTrackNorm = 0 ;
  
  if ( partTypeInCone != AliIsolationCut::kOnlyNeutral )
  {
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateTrackUEBand   (pCandidate,etaUEptsumTrack  ,phiUEptsumTrack  );// rajouter ici l'histo eta phi
    
    // Fill histograms
    fhConeSumPtVSUETracksEtaBand->Fill(coneptsumTrack, etaUEptsumTrack, GetEventWeight());
    fhConeSumPtVSUETracksPhiBand->Fill(coneptsumTrack, phiUEptsumTrack, GetEventWeight());
    
    
    Float_t correctConeSumTrack    = 1;
    Float_t correctConeSumTrackPhi = 1;
    
    GetIsolationCut()->CalculateUEBandTrackNormalization(GetReader(),etaTrig, phiTrig,
                                                         phiUEptsumTrack,etaUEptsumTrack,
                                                         phiUEptsumTrackNorm,etaUEptsumTrackNorm,
                                                         correctConeSumTrack,correctConeSumTrackPhi);
    
    coneptsumTrackSubPhi = coneptsumTrack - phiUEptsumTrackNorm;
    coneptsumTrackSubEta = coneptsumTrack - etaUEptsumTrackNorm;
    
    etaBandptsumTrackNorm = etaUEptsumTrackNorm;
    
    fhConeSumPtPhiUESubTrack           ->Fill(ptTrig ,          coneptsumTrackSubPhi, GetEventWeight());
    fhConeSumPtPhiUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubPhi *GetEventWeight()); // check
    fhConeSumPtEtaUESubTrack           ->Fill(ptTrig ,          coneptsumTrackSubEta, GetEventWeight());
    fhConeSumPtEtaUESubTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrackSubEta *GetEventWeight()); // check
    
    fhFractionTrackOutConeEta          ->Fill(ptTrig ,         correctConeSumTrack-1, GetEventWeight());
    fhFractionTrackOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig,correctConeSumTrack-1 *GetEventWeight()); // check
    
    if(coneptsumTrack > 0)
    {
      coneptsumTrackSubPhiNorm = coneptsumTrackSubPhi/coneptsumTrack;
      coneptsumTrackSubEtaNorm = coneptsumTrackSubEta/coneptsumTrack;
    }
    
    fhConeSumPtSubvsConeSumPtTotPhiTrack    ->Fill(coneptsumTrack, coneptsumTrackSubPhi    , GetEventWeight());
    fhConeSumPtSubNormvsConeSumPtTotPhiTrack->Fill(coneptsumTrack, coneptsumTrackSubPhiNorm, GetEventWeight());
    fhConeSumPtSubvsConeSumPtTotEtaTrack    ->Fill(coneptsumTrack, coneptsumTrackSubEta    , GetEventWeight());
    fhConeSumPtSubNormvsConeSumPtTotEtaTrack->Fill(coneptsumTrack, coneptsumTrackSubEtaNorm, GetEventWeight());
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
  etaBandptsumClusterNorm = 0;
  
  if ( partTypeInCone != AliIsolationCut::kOnlyCharged )
  {
    // -------------- //
    // EMCal clusters //
    // -------------- //
    
    // Sum the pT in the phi or eta band for clusters or tracks
    CalculateCaloUEBand    (pCandidate,etaUEptsumCluster,phiUEptsumCluster);// rajouter ici l'histo eta phi
    
    // Fill histograms
    fhConeSumPtVSUEClusterEtaBand->Fill(coneptsumCluster, etaUEptsumCluster, GetEventWeight());
    fhConeSumPtVSUEClusterPhiBand->Fill(coneptsumCluster, phiUEptsumCluster, GetEventWeight());
    
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
    
    etaBandptsumClusterNorm = etaUEptsumClusterNorm;
    
    fhConeSumPtPhiUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubPhi, GetEventWeight());
    fhConeSumPtPhiUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubPhi *GetEventWeight()); // check
    fhConeSumPtEtaUESubCluster           ->Fill(ptTrig ,          coneptsumClusterSubEta, GetEventWeight());
    fhConeSumPtEtaUESubClusterTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumClusterSubEta *GetEventWeight()); // check
    
    fhFractionClusterOutConeEta          ->Fill(ptTrig ,          correctConeSumClusterEta-1, GetEventWeight());
    fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterEta-1 *GetEventWeight()); // check
    fhFractionClusterOutConePhi          ->Fill(ptTrig ,          correctConeSumClusterPhi-1, GetEventWeight());
    fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaTrig, phiTrig, correctConeSumClusterPhi-1 *GetEventWeight()); // check
    
    if(coneptsumCluster!=0)
    {
      coneptsumClusterSubPhiNorm = coneptsumClusterSubPhi/coneptsumCluster;
      coneptsumClusterSubEtaNorm = coneptsumClusterSubEta/coneptsumCluster;
    }
    
    fhConeSumPtSubvsConeSumPtTotPhiCluster    ->Fill(coneptsumCluster,coneptsumClusterSubPhi    , GetEventWeight());
    fhConeSumPtSubNormvsConeSumPtTotPhiCluster->Fill(coneptsumCluster,coneptsumClusterSubPhiNorm, GetEventWeight());
    fhConeSumPtSubvsConeSumPtTotEtaCluster    ->Fill(coneptsumCluster,coneptsumClusterSubEta    , GetEventWeight());
    fhConeSumPtSubNormvsConeSumPtTotEtaCluster->Fill(coneptsumCluster,coneptsumClusterSubEtaNorm, GetEventWeight());
    
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
  
  if ( partTypeInCone == AliIsolationCut::kNeutralAndCharged )
  {
    // --------------------------- //
    // Tracks and clusters in cone //
    // --------------------------- //
    
    Double_t sumPhiUESub = coneptsumClusterSubPhi + coneptsumTrackSubPhi;
    Double_t sumEtaUESub = coneptsumClusterSubEta + coneptsumTrackSubEta;
    
    fhConeSumPtPhiUESub          ->Fill(ptTrig ,          sumPhiUESub, GetEventWeight());
    fhConeSumPtPhiUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumPhiUESub *GetEventWeight()); // check
    fhConeSumPtEtaUESub          ->Fill(ptTrig ,          sumEtaUESub, GetEventWeight());
    fhConeSumPtEtaUESubTrigEtaPhi->Fill(etaTrig, phiTrig, sumEtaUESub *GetEventWeight()); // check
    
    fhEtaBandClustervsTrack    ->Fill(etaUEptsumCluster    ,etaUEptsumTrack    , GetEventWeight());
    fhPhiBandClustervsTrack    ->Fill(phiUEptsumCluster    ,phiUEptsumTrack    , GetEventWeight());
    fhEtaBandNormClustervsTrack->Fill(etaUEptsumClusterNorm,etaUEptsumTrackNorm, GetEventWeight());
    fhPhiBandNormClustervsTrack->Fill(phiUEptsumClusterNorm,phiUEptsumTrackNorm, GetEventWeight());
    
    fhConeSumPtEtaUESubClustervsTrack->Fill(coneptsumClusterSubEta, coneptsumTrackSubEta, GetEventWeight());
    fhConeSumPtPhiUESubClustervsTrack->Fill(coneptsumClusterSubPhi, coneptsumTrackSubPhi, GetEventWeight());
    
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
  
  // Recover reference arrays with clusters and tracks
  TObjArray * refclusters = aodParticle->GetObjArray(GetAODObjArrayName()+"Clusters");
  if(!refclusters) return ;
  
  Float_t ptTrig = aodParticle->Pt();
  
  // Get vertex for cluster momentum calculation
  Double_t vertex[] = {0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    GetReader()->GetVertex(vertex);
  
  Float_t ptcone        = 0;
  Float_t distToTrigger = 0;
  for(Int_t icalo=0; icalo < refclusters->GetEntriesFast(); icalo++)
  {
    AliVCluster* calo = (AliVCluster *) refclusters->At(icalo);
    calo->GetMomentum(fMomentum,vertex) ;//Assume that come from vertex in straight line
    
    ptcone = fMomentum.Pt();
    
    fhPtInCone       ->Fill(ptTrig, ptcone, GetEventWeight());
    fhPtClusterInCone->Fill(ptTrig, ptcone, GetEventWeight());
    
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
  }
  
  fhConeSumPtCluster ->Fill(ptTrig, coneptsumCluster , GetEventWeight());
  fhConePtLeadCluster->Fill(ptTrig, coneptLeadCluster, GetEventWeight());
  
  aodParticle->SetNeutralLeadPtInCone(coneptLeadCluster);
  aodParticle->SetNeutralPtSumInCone(coneptsumCluster);
}

//______________________________________________________________________________________________________
/// Get the cell amplityde or sum of amplitudes in isolation cone.
/// Missing: Remove signal cells in cone in case the trigger is a cluster.
//______________________________________________________________________________________________________
void AliAnaParticleIsolation::CalculateCaloCellSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle,
                                                            Float_t & coneptsumCell)
{
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyCharged ) return ;
  
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
  if( GetIsolationCut()->GetParticleTypeInCone()==AliIsolationCut::kOnlyNeutral ) return ;
  
  // Recover reference arrays with clusters and tracks
  TObjArray * reftracks   = aodParticle->GetObjArray(GetAODObjArrayName()+"Tracks");
  if(!reftracks) return ;
  
  Float_t  ptTrig = aodParticle->Pt();
  Double_t bz     = GetReader()->GetInputEvent()->GetMagneticField();
  
  Float_t pTtrack       = 0;
  Float_t distToTrigger = 0;
  for(Int_t itrack=0; itrack < reftracks->GetEntriesFast(); itrack++)
  {
    AliVTrack* track = (AliVTrack *) reftracks->At(itrack);
    
    pTtrack       = track->Pt();
    
    fhPtInCone     ->Fill(ptTrig, pTtrack, GetEventWeight());
    fhPtTrackInCone->Fill(ptTrig, pTtrack, GetEventWeight());
    
    if(IsPileUpAnalysisOn())
    {
      ULong_t status = track->GetStatus();
      Bool_t okTOF = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
      //Double32_t tof = track->GetTOFsignal()*1e-3;
      Int_t trackBC = track->GetTOFBunchCrossing(bz);
      
      if     ( okTOF && trackBC!=0 ) fhPtTrackInConeOtherBC->Fill(ptTrig,pTtrack, GetEventWeight());
      else if( okTOF && trackBC==0 ) fhPtTrackInConeBC0    ->Fill(ptTrig,pTtrack, GetEventWeight());
      
      Int_t vtxBC = GetReader()->GetVertexBC();
      if(vtxBC == 0 || vtxBC==AliVTrack::kTOFBCNA) fhPtTrackInConeVtxBC0->Fill(ptTrig, pTtrack, GetEventWeight());
      
      if(GetReader()->IsPileUpFromSPD())             {   fhPtInConePileUp[0]            ->Fill(ptTrig, pTtrack, GetEventWeight());
        if(okTOF && trackBC!=0 )                         fhPtTrackInConeOtherBCPileUpSPD->Fill(ptTrig, pTtrack, GetEventWeight());
        if(okTOF && trackBC==0 )                         fhPtTrackInConeBC0PileUpSPD    ->Fill(ptTrig, pTtrack, GetEventWeight()); }
      if(GetReader()->IsPileUpFromEMCal())             fhPtInConePileUp[1]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtInConePileUp[2]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtInConePileUp[3]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtInConePileUp[4]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtInConePileUp[5]->Fill(ptTrig, pTtrack, GetEventWeight());
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtInConePileUp[6]->Fill(ptTrig, pTtrack, GetEventWeight());
    }
    
    if(IsHighMultiplicityAnalysisOn()) fhPtInConeCent->Fill(GetEventCentrality(), pTtrack, GetEventWeight());
    
    coneptsumTrack+=pTtrack;
    if(pTtrack > coneptLeadTrack) coneptLeadTrack = pTtrack;
  }

  fhConeSumPtTrack ->Fill(ptTrig, coneptsumTrack , GetEventWeight());
  fhConePtLeadTrack->Fill(ptTrig, coneptLeadTrack, GetEventWeight());
  
  aodParticle->SetChargedLeadPtInCone(coneptLeadTrack);
  aodParticle->SetChargedPtSumInCone(coneptsumTrack);
}

//_____________________________________________________________________________
/// Fill some histograms to understand pile-up.
//_____________________________________________________________________________
void AliAnaParticleIsolation::FillPileUpHistograms(Float_t energy, Float_t time)//Int_t clusterID)
{
//  if ( clusterID < 0 )
//  {
//    AliWarning(Form("ID of cluster = %d, not possible!", clusterID));
//    return;
//  }
  
//  Int_t iclus = -1;
//  TObjArray* clusters = 0x0;
//  if     (GetCalorimeter() == kEMCAL) clusters = GetEMCALClusters();
//  else if(GetCalorimeter() == kPHOS ) clusters = GetPHOSClusters();
//  
//  Float_t energy = 0;
//  Float_t time   = -1000;
//  
//  if(clusters)
//  {
//    AliVCluster *cluster = FindCluster(clusters,clusterID,iclus);
//    energy = cluster->E();
//    time   = cluster->GetTOF()*1e9;
//  }
  
  //printf("E %f, time %f\n",energy,time);
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
void AliAnaParticleIsolation::FillTrackMatchingShowerShapeControlHistograms(AliAODPWG4ParticleCorrelation  *pCandidate,
                                                                            Float_t coneptsum, Float_t coneleadpt,
                                                                            Int_t mcIndex)
{
  if(!fFillTMHisto && !fFillSSHisto && !fFillBackgroundBinHistograms && !fFillTaggedDecayHistograms) return;
  
  Int_t  clusterID = pCandidate->GetCaloLabel(0) ;
  Int_t  nMaxima   = pCandidate->GetNLM();
  Int_t  mcTag     = pCandidate->GetTag() ;
  Bool_t isolated  = pCandidate->IsIsolated();

  if ( clusterID < 0 )
  {
    AliWarning(Form("ID of cluster = %d, not possible!", clusterID));
    return;
  }

  Float_t m02    = pCandidate->GetM02() ;
  Float_t energy = pCandidate->E();
  Float_t pt     = pCandidate->Pt();
  Float_t eta    = pCandidate->Eta();
  Float_t phi    = pCandidate->Phi();
  if(phi<0) phi+= TMath::TwoPi();
  
  // Recover original cluster if requested
  //
  if(fFillOverlapHistograms || fFillTMHisto)
  {
    Int_t iclus = -1;
    if     (GetCalorimeter() == kEMCAL) fClustersArr = GetEMCALClusters();
    else if(GetCalorimeter() == kPHOS ) fClustersArr = GetPHOSClusters();
    
    if(!fClustersArr) return;
    
    fCluster = FindCluster(fClustersArr,clusterID,iclus);
  }

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
  
  if(fFillPtTrigBinHistograms)
  {
    // Get the background bin for this cone and trigger
    Int_t ptTrigBin  = -1;
    
    for(Int_t ibin = 0; ibin < fNPtTrigBin; ibin++)
    {
      if( pt  >= fPtTrigBinLimit[ibin] && coneptsum  < fPtTrigBinLimit[ibin+1]) ptTrigBin  = ibin;
    }
    
    // Fill the histograms per pT candidate bin of pt lead or pt sum
    
    if ( ptTrigBin >= 0 )
    {
      AliDebug(1,Form("Trigger pT %f, bin %d [%2.2f,%2.2f]",pt,ptTrigBin,fPtTrigBinLimit[ptTrigBin],fPtTrigBinLimit[ptTrigBin+1]));

      fhPtTrigBinPtLeadCone[ptTrigBin]->Fill(coneleadpt, GetEventWeight());
      fhPtTrigBinSumPtCone [ptTrigBin]->Fill(coneptsum , GetEventWeight());
        
      if(fFillSSHisto)
      {
        fhPtTrigBinLambda0vsPtLeadCone[ptTrigBin]->Fill(coneleadpt, m02, GetEventWeight());
        fhPtTrigBinLambda0vsSumPtCone [ptTrigBin]->Fill(coneptsum , m02, GetEventWeight());
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
          }
        }
      }
    }
    
    if( IsDataMC() )
    {
      Int_t ptTrigBinMC = ptTrigBin+mcIndex*fNPtTrigBin;
      
      if( ptTrigBin >=0 )
      {
        fhPtTrigBinPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, GetEventWeight());
        fhPtTrigBinSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , GetEventWeight());
          
        if(fFillSSHisto)
        {
          fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, m02, GetEventWeight());
          fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , m02, GetEventWeight());
        }
      }
      
      if(GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton))
      {
        ptTrigBinMC = ptTrigBin+kmcPhoton*fNPtTrigBin;
        if( ptTrigBin >=0 )
        {
          fhPtTrigBinPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, GetEventWeight());
          fhPtTrigBinSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , GetEventWeight());
            
          if(fFillSSHisto)
          {
            fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, m02, GetEventWeight());
            fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , m02, GetEventWeight());
          }
        }
      } // photon MC
      
      // decays with lost pair
      if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
      {
        if     ( mcIndex == kmcPi0Decay ) ptTrigBinMC = ptTrigBin+kmcPi0DecayLostPair*fNPtTrigBin;
        else if( mcIndex == kmcEtaDecay ) ptTrigBinMC = ptTrigBin+kmcEtaDecayLostPair*fNPtTrigBin;
        
        if( ptTrigBin >=0 )
        {
          fhPtTrigBinPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, GetEventWeight());
          fhPtTrigBinSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , GetEventWeight());
          if(fFillSSHisto)
          {
            fhPtTrigBinLambda0vsPtLeadConeMC[ptTrigBinMC]->Fill(coneleadpt, m02, GetEventWeight());
            fhPtTrigBinLambda0vsSumPtConeMC [ptTrigBinMC]->Fill(coneptsum , m02, GetEventWeight());
          }
        }
      } // lost decays MC

    } // MC
  } // pT trigger bins
  
  // Shower shape dependent histograms
  if(fFillSSHisto)
  {
    fhELambda0 [isolated]->Fill(energy, m02, GetEventWeight());
    fhPtLambda0[isolated]->Fill(pt,     m02, GetEventWeight());
    //fhELambda1 [isolated]->Fill(energy, m20);
    
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

      
      Int_t noverlaps = 0;
      if(fFillOverlapHistograms)
      {
        const UInt_t nlabels = fCluster->GetNLabels();
        Int_t overpdg[nlabels];
        noverlaps = GetMCAnalysisUtils()->GetNOverlaps(fCluster->GetLabels(), nlabels,mcTag,-1,GetReader(),overpdg);
        
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
      }
      
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
      else if (noverlaps > 1 ) // More than 1 overlap
      {
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
          fhPtLambda0MCWithNOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
        {
          if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCWithNOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
          else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCWithNOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
        }
        
        fhPtLambda0MCWithNOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
        
        if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) )
        {
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) )
            fhPtLambda0MCConvWithNOverlap[kmcPhoton][isolated]->Fill(pt, m02, GetEventWeight());
          
          if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCDecayPairLost) )
          {
            if     ( mcIndex == kmcPi0Decay ) fhPtLambda0MCConvWithNOverlap[kmcPi0DecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
            else if( mcIndex == kmcEtaDecay ) fhPtLambda0MCConvWithNOverlap[kmcEtaDecayLostPair][isolated]->Fill(pt, m02, GetEventWeight());
          }
          
          fhPtLambda0MCConvWithNOverlap[mcIndex][isolated]->Fill(pt, m02, GetEventWeight());
        } // Conversion
      } // more than 1 overlap
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
  if(fFillTMHisto)
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
  snprintf(onePar, buffersize,"Isolation Cand Detector: %s;",fIsoDetectorString.Data()) ;
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
  
  TString parTitle = Form("#it{R} = %2.2f%s%s",GetIsolationCut()->GetConeSize(), sThreshold.Data(),sParticle.Data());
  
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
  
  fhEtaPhiNoIso  = new TH2F("hEtaPhiNoIso",
                            Form("Number of not isolated leading particles #eta vs #phi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiNoIso->SetXTitle("#eta");
  fhEtaPhiNoIso->SetYTitle("#phi");
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
                                  Form("#phi vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
                                  nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhPhiIsoMC[imc]->SetYTitle("#phi");
      fhPhiIsoMC[imc]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
      outputContainer->Add(fhPhiIsoMC[imc]) ;
      
      fhEtaIsoMC[imc]  = new TH2F(Form("hEtaMC%s",mcPartName[imc].Data()),
                                  Form("#phi vs #it{p}_{T} of isolated %s, %s",mcPartType[imc].Data(),parTitle.Data()),
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
                 Form("Number of %s leading Pi0 decay particles #eta vs #phi, bit %d, %s",isoTitle[iso].Data(),fDecayBits[ibit],parTitle.Data()),
                 netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiDecay[iso][ibit]->SetXTitle("#eta");
        fhEtaPhiDecay[iso][ibit]->SetYTitle("#phi");
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
    
    fhPhiIso  = new TH2F("hPhi",
                         Form("Number of isolated particles vs #phi, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhPhiIso->SetYTitle("#phi");
    fhPhiIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPhiIso) ;
    
    fhEtaIso  = new TH2F("hEta",
                         Form("Number of isolated particles vs #eta, %s",parTitle.Data()),
                         nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhEtaIso->SetYTitle("#eta");
    fhEtaIso->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhEtaIso) ;
    
    fhEtaPhiIso  = new TH2F("hEtaPhiIso",
                            Form("Number of isolated particles #eta vs #phi, %s",parTitle.Data()),
                            netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiIso->SetXTitle("#eta");
    fhEtaPhiIso->SetYTitle("#phi");
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
    
    fhConeSumPt  = new TH2F("hConePtSum",
                            Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                            nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma #it{p}_{T}");
    fhConeSumPt->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPt) ;
    
    fhConeSumPtTrigEtaPhi  = new TH2F("hConePtSumTrigEtaPhi",
                                      Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                      netabins,etamin,etamax,nphibins,phimin,phimax);
    fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
    fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
    fhConeSumPtTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
    outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
    
    fhPtInCone  = new TH2F("hPtInCone",
                           Form("#it{p}_{T} of clusters and tracks in isolation cone for #it{R} =  %2.2f",r),
                           nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInCone) ;
    
    if(fFillBackgroundBinHistograms)
    {
      fhPtLeadConeBin              = new TH1F*[fNBkgBin];
      fhSumPtConeBin               = new TH1F*[fNBkgBin];
      if(fFillSSHisto)
      {
        fhPtLeadConeBinLambda0     = new TH2F*[fNBkgBin];
        fhSumPtConeBinLambda0      = new TH2F*[fNBkgBin];
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
        
        if(fFillSSHisto)
        {
          fhPtLeadConeBinLambda0MC = new TH2F*[fNBkgBin*fgkNmcTypes];
          fhSumPtConeBinLambda0MC  = new TH2F*[fNBkgBin*fgkNmcTypes];
        }
      }
      
      for(Int_t ibin = 0; ibin < fNBkgBin; ibin++)
      {
        fhPtLeadConeBin[ibin]  = new TH1F
        (Form("hPtLeadCone_Bin%d",ibin),
         Form("cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
              fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
        fhPtLeadConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtLeadConeBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtLeadConeBin[ibin]) ;
        
        fhSumPtConeBin[ibin]  = new TH1F
        (Form("hSumPtCone_Bin%d",ibin),
         Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
              fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
        fhSumPtConeBin[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhSumPtConeBin[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhSumPtConeBin[ibin]) ;
        
        if(fFillTaggedDecayHistograms)
        {
          for(Int_t idecay = 0; idecay < fNDecayBits; idecay++)
          {
            Int_t bindecay = ibin+idecay*fNBkgBin;

            fhPtLeadConeBinDecay[bindecay]  = new TH1F
            (Form("hPtLeadCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
            fhPtLeadConeBinDecay[bindecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtLeadConeBinDecay[bindecay]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtLeadConeBinDecay[bindecay]) ;
            
            fhSumPtConeBinDecay[bindecay]  = new TH1F
            (Form("hSumPtCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c},  %s",
                  fDecayBits[idecay],fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
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
                  fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax);
            fhPtLeadConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtLeadConeBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhPtLeadConeBinMC[binmc]) ;
            
            fhSumPtConeBinMC[binmc]  = new TH1F
            (Form("hSumPtCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
             Form("in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                  fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax);
            fhSumPtConeBinMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhSumPtConeBinMC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            outputContainer->Add(fhSumPtConeBinMC[binmc]) ;
          } // MC particle loop
        }
        
        if(fFillSSHisto)
        {
          fhPtLeadConeBinLambda0[ibin]  = new TH2F
          (Form("hPtLeadConeLambda0_Bin%d",ibin),
           Form("#lambda_{0}, in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLeadConeBinLambda0[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtLeadConeBinLambda0[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLeadConeBinLambda0[ibin]) ;
          
          fhSumPtConeBinLambda0[ibin]  = new TH2F
          (Form("hSumPtConeLambda0_Bin%d",ibin),
           Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, %s",
                fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhSumPtConeBinLambda0[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhSumPtConeBinLambda0[ibin]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhSumPtConeBinLambda0[ibin]) ;
          
          if(IsDataMC())
          {
            for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
            {
              Int_t binmc = ibin+imc*fNBkgBin;
              fhPtLeadConeBinLambda0MC[binmc]  = new TH2F
              (Form("hPtLeadConeLambda0_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
               Form("#lambda_{0}, in cone %2.2f<#it{p}_{T}^{leading}<%2.2f GeV/#it{c}, MC %s, %s",
                    fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLeadConeBinLambda0MC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtLeadConeBinLambda0MC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              outputContainer->Add(fhPtLeadConeBinLambda0MC[binmc]) ;
              
              fhSumPtConeBinLambda0MC[binmc]  = new TH2F
              (Form("hSumPtConeLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#lambda_{0}, in cone %2.2f <#Sigma #it{p}_{T}< %2.2f GeV/#it{c}, MC %s, %s",
                    fBkgBinLimit[ibin],fBkgBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhSumPtConeBinLambda0MC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhSumPtConeBinLambda0MC[binmc]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
              outputContainer->Add(fhSumPtConeBinLambda0MC[binmc]) ;
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
      
      if(IsDataMC())
      {
        fhPtTrigBinPtLeadConeMC = new TH1F*[fNPtTrigBin*fgkNmcTypes];
        fhPtTrigBinSumPtConeMC  = new TH1F*[fNPtTrigBin*fgkNmcTypes];
      }

      if(fFillSSHisto)
      {
        fhPtTrigBinLambda0vsPtLeadCone = new TH2F*[fNPtTrigBin];
        fhPtTrigBinLambda0vsSumPtCone  = new TH2F*[fNPtTrigBin];
        
        if(IsDataMC())
        {
          fhPtTrigBinLambda0vsPtLeadConeMC = new TH2F*[fNPtTrigBin*fgkNmcTypes];
          fhPtTrigBinLambda0vsSumPtConeMC  = new TH2F*[fNPtTrigBin*fgkNmcTypes];
        }
      }
      
      for(Int_t ibin = 0; ibin < fNPtTrigBin; ibin++)
      {
        fhPtTrigBinPtLeadCone[ibin]  = new TH1F
        (Form("hPtTrigBin_PtLeadCone_Bin%d",ibin),
         Form("#it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
              fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
        fhPtTrigBinPtLeadCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtTrigBinPtLeadCone[ibin]->SetXTitle("#it{p}_{T}^{in cone} (GeV/#it{c})");
        outputContainer->Add(fhPtTrigBinPtLeadCone[ibin]) ;
        
        fhPtTrigBinSumPtCone[ibin]  = new TH1F
        (Form("hPtTrigBin_SumPtCone_Bin%d",ibin),
         Form("#Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
              fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptsumbins,ptsummin,ptsummax);
        fhPtTrigBinSumPtCone[ibin]->SetYTitle("d #it{N} / d #it{p}_{T}");
        fhPtTrigBinSumPtCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
        outputContainer->Add(fhPtTrigBinSumPtCone[ibin]) ;

        if(fFillTaggedDecayHistograms)
        {
          for(Int_t idecay = 0; idecay < fNDecayBits; idecay++)
          {
            Int_t binDecay = ibin+idecay*fNPtTrigBin;
            
            fhPtTrigBinPtLeadConeDecay[binDecay]  = new TH1F
            (Form("hPtTrigBin_PtLeadCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax);
            fhPtTrigBinPtLeadConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinPtLeadConeDecay[binDecay]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinPtLeadConeDecay[binDecay]) ;
            
            fhPtTrigBinSumPtConeDecay[binDecay]  = new TH1F
            (Form("hPtTrigBin_SumPtCone_Bin%d_DecayBit%d",ibin,fDecayBits[idecay]),
             Form("Decay bit %d, #Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                  fDecayBits[idecay],fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptsumbins,ptsummin,ptsummax);
            fhPtTrigBinSumPtConeDecay[binDecay]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinSumPtConeDecay[binDecay]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinSumPtConeDecay[binDecay]) ;
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
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax);
            fhPtTrigBinPtLeadConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinPtLeadConeMC[binmc]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinPtLeadConeMC[binmc]) ;
            
            fhPtTrigBinSumPtConeMC[binmc]  = new TH1F
            (Form("hPtTrigBin_SumPtCone_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
             Form("#Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                  fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptsumbins,ptsummin,ptsummax);
            fhPtTrigBinSumPtConeMC[binmc]->SetYTitle("d #it{N} / d #it{p}_{T}");
            fhPtTrigBinSumPtConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
            outputContainer->Add(fhPtTrigBinSumPtConeMC[binmc]) ;
          } // MC particle loop
        } // MC

        if(fFillSSHisto)
        {
          fhPtTrigBinLambda0vsPtLeadCone[ibin]  = new TH2F
          (Form("hPtTrigBin_PtLeadConeVSLambda0_Bin%d",ibin),
           Form("#lambda_{0} vs #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtTrigBinLambda0vsPtLeadCone[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtTrigBinLambda0vsPtLeadCone[ibin]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinLambda0vsPtLeadCone[ibin]) ;
          
          fhPtTrigBinLambda0vsSumPtCone[ibin]  = new TH2F
          (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d",ibin),
           Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone} %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, %s",
                fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], parTitle.Data()),nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
          fhPtTrigBinLambda0vsSumPtCone[ibin]->SetYTitle("#lambda_{0}^{2}");
          fhPtTrigBinLambda0vsSumPtCone[ibin]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
          outputContainer->Add(fhPtTrigBinLambda0vsSumPtCone[ibin]) ;
          
          if(IsDataMC())
          {
            for(Int_t imc = 0; imc < fgkNmcTypes; imc++)
            {
              Int_t binmc = ibin+imc*fNPtTrigBin;
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]  = new TH2F
              (Form("hPtTrigBin_PtLeadConeVSLambda0_Bin%d_MC%s",ibin, mcPartName[imc].Data()),
               Form("#lambda_{0} vs #it{p}_{T}^{lead. in cone}, %2.2f<#it{p}_{T}^{cand}<%2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtTrigBinLambda0vsPtLeadConeMC[binmc]->SetXTitle("#it{p}_{T}^{lead in cone} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinLambda0vsPtLeadConeMC[binmc]) ;
              
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]  = new TH2F
              (Form("hPtTrigBin_SumPtConeVSLambda0_Bin%d_MC%s",ibin,mcPartName[imc].Data()),
               Form("#lambda_{0} vs #Sigma #it{p}_{T}^{in cone}, %2.2f <#it{p}_{T}^{cand}< %2.2f GeV/#it{c}, MC %s, %s",
                    fPtTrigBinLimit[ibin],fPtTrigBinLimit[ibin+1], mcPartType[imc].Data(), parTitle.Data()),nptsumbins,ptsummin,ptsummax,ssbins,ssmin,ssmax);
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]->SetYTitle("#lambda_{0}^{2}");
              fhPtTrigBinLambda0vsSumPtConeMC[binmc]->SetXTitle("#Sigma #it{p}_{T}^{in cone} (GeV/#it{c})");
              outputContainer->Add(fhPtTrigBinLambda0vsSumPtConeMC[binmc]) ;
            } // MC particle loop
          } // MC
        } // SS histo
      } // pt trig bin loop
    } // pt trig bin histograms
    
    if(IsHighMultiplicityAnalysisOn())
    {
      fhPtInConeCent  = new TH2F("hPtInConeCent",
                                 Form("#it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                 100,0,100,nptinconebins,ptinconemin,ptinconemax);
      fhPtInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtInConeCent->SetXTitle("centrality");
      outputContainer->Add(fhPtInConeCent) ;
    }
    
    // Cluster only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyCharged)
    {
      fhConeSumPtCluster  = new TH2F("hConePtSumCluster",
                                     Form("Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtCluster->SetYTitle("#Sigma #it{p}_{T}");
      fhConeSumPtCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtCluster) ;
      
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
        fhConeSumPtCell->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUECluster  = new TH2F("hConePtSumEtaBandUECluster",
                                                "#Sigma cluster #it{p}_{T} in UE Eta Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUECluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUECluster) ;
        
        fhConeSumPtPhiBandUECluster  = new TH2F("hConePtSumPhiBandUECluster",
                                                "#Sigma cluster #it{p}_{T} UE Phi Band",
                                                nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtPhiBandUECluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUECluster) ;
        
        fhConeSumPtEtaBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumEtaBandUEClusterTrigEtaPhi",
                                                          "Trigger #eta vs #phi, #Sigma cluster #it{p}_{T} in UE Eta Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUEClusterTrigEtaPhi  = new TH2F("hConePtSumPhiBandUEClusterTrigEtaPhi",
                                                          "Trigger #eta vs #phi, #Sigma cluster #it{p}_{T} UE Phi Band",
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhi) ;
        if(fFillCellHistograms)
        {
          
          fhConeSumPtEtaBandUECell  = new TH2F("hConePtSumEtaBandUECell",
                                               "#Sigma cell #it{p}_{T} in UE Eta Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtEtaBandUECell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaBandUECell) ;
          
          fhConeSumPtPhiBandUECell  = new TH2F("hConePtSumPhiBandUECell",
                                               "#Sigma cell #it{p}_{T} UE Phi Band",
                                               nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPhiBandUECell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUECell->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiBandUECell) ;
          
          fhConeSumPtEtaBandUECellTrigEtaPhi  = new TH2F("hConePtSumEtaBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #phi, #Sigma cell #it{p}_{T} in UE Eta Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaBandUECellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaBandUECellTrigEtaPhi) ;
          
          fhConeSumPtPhiBandUECellTrigEtaPhi  = new TH2F("hConePtSumPhiBandUECellTrigEtaPhi",
                                                         "Trigger #eta vs #phi, #Sigma cell #it{p}_{T} UE Phi Band",
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiBandUECellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiBandUECellTrigEtaPhi) ;
        }
        
        fhEtaBandCluster  = new TH2F("hEtaBandCluster",
                                     Form("#eta vs #phi of clusters in #eta band isolation cone for #it{R} =  %2.2f",r),
                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandCluster->SetXTitle("#eta");
        fhEtaBandCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaBandCluster) ;
        
        fhPhiBandCluster  = new TH2F("hPhiBandCluster",
                                     Form("#eta vs #phi of clusters in #phi band isolation cone for #it{R} =  %2.2f",r),
                                     netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandCluster->SetXTitle("#eta");
        fhPhiBandCluster->SetYTitle("#phi");
        outputContainer->Add(fhPhiBandCluster) ;
        
        fhEtaPhiInConeCluster= new TH2F("hEtaPhiInConeCluster",
                                        Form("#eta vs #phi of clusters in cone for #it{R} =  %2.2f",r),
                                        netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInConeCluster->SetXTitle("#eta");
        fhEtaPhiInConeCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiInConeCluster) ;
        
        fhEtaPhiCluster= new TH2F("hEtaPhiCluster",
                                  Form("#eta vs #phi of all clusters"),
                                  netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiCluster->SetXTitle("#eta");
        fhEtaPhiCluster->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiCluster) ;
        
      }
      
      fhPtClusterInCone  = new TH2F("hPtClusterInCone",
                                    Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtClusterInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtClusterInCone) ;
      
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
                                  Form("#col vs #row of cells in #phi band isolation cone for #it{R} =  %2.2f",r),
                                  96,0,95,128,0,127);
        fhPhiBandCell->SetXTitle("#col");
        fhPhiBandCell->SetYTitle("#row");
        outputContainer->Add(fhPhiBandCell) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUESubCluster  = new TH2F("hConeSumPtEtaUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubCluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubCluster) ;
        
        fhConeSumPtPhiUESubCluster  = new TH2F("hConeSumPtPhiUESubCluster",
                                               Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                               nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubCluster->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubCluster) ;
        
        fhConeSumPtEtaUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubClusterTrigEtaPhi",
                                                         Form("Trigger #eta vs #phi, Clusters #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubClusterTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubClusterTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubClusterTrigEtaPhi",
                                                         Form("Trigger #eta vs #phi, Clusters #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubClusterTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubClusterTrigEtaPhi) ;
        
        if(fFillCellHistograms)
        {
          fhConeSumPtEtaUESubCell  = new TH2F("hConeSumPtEtaUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubCell) ;
          
          fhConeSumPtPhiUESubCell  = new TH2F("hConeSumPtPhiUESubCell",
                                              Form("Cells #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                              nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubCell) ;
          
          fhConeSumPtEtaUESubCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #phi, Cells #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubCellTrigEtaPhi",
                                                        Form("Trigger #eta vs #phi, Cells #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtPhiUESubCellTrigEtaPhi) ;
        }
        
        fhFractionClusterOutConeEta  = new TH2F("hFractionClusterOutConeEta",
                                                Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance",r),
                                                nptbins,ptmin,ptmax,100,0,1);
        fhFractionClusterOutConeEta->SetYTitle("#it{fraction}");
        fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionClusterOutConeEta) ;
        
        fhFractionClusterOutConeEtaTrigEtaPhi  = new TH2F("hFractionClusterOutConeEtaTrigEtaPhi",
                                                          Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance, in trigger #eta-#phi ",r),
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;
        
        fhFractionClusterOutConePhi  = new TH2F("hFractionClusterOutConePhi",
                                                Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #phi acceptance",r),
                                                nptbins,ptmin,ptmax,100,0,1);
        fhFractionClusterOutConePhi->SetYTitle("#it{fraction}");
        fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionClusterOutConePhi) ;
        
        fhFractionClusterOutConePhiTrigEtaPhi  = new TH2F("hFractionClusterOutConePhiTrigEtaPhi",
                                                          Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #phi acceptance, in trigger #eta-#phi ",r),
                                                          netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;
        
        fhConeSumPtSubvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCluster",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCluster);
        
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCluster",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCluster);
        
        fhConeSumPtSubvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCluster",
                                                          Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                          nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCluster);
        
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCluster",
                                                              Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                              nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCluster);
        
        fhConeSumPtVSUEClusterEtaBand  = new TH2F("hConeSumPtVSUEClusterEtaBand",
                                                  Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in eta band for cluster (before normalization), R=%2.2f",r),
                                                  nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
        fhConeSumPtVSUEClusterEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUEClusterEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUEClusterEtaBand);
        
        fhConeSumPtVSUEClusterPhiBand  = new TH2F("hConeSumPtVSUEClusterPhiBand",
                                                  Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in phi band for cluster (before normalization), R=%2.2f",r),
                                                  nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
        fhConeSumPtVSUEClusterPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUEClusterPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUEClusterPhiBand);
        
        if(fFillCellHistograms)
        {
          fhFractionCellOutConeEta  = new TH2F("hFractionCellOutConeEta",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConeEta->SetYTitle("#it{fraction}");
          fhFractionCellOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConeEta) ;
          
          fhFractionCellOutConeEtaTrigEtaPhi  = new TH2F("hFractionCellOutConeEtaTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #eta acceptance, in trigger #eta-#phi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConeEtaTrigEtaPhi) ;
          
          fhFractionCellOutConePhi  = new TH2F("hFractionCellOutConePhi",
                                               Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #phi acceptance",r),
                                               nptbins,ptmin,ptmax,100,0,1);
          fhFractionCellOutConePhi->SetYTitle("#it{fraction}");
          fhFractionCellOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
          outputContainer->Add(fhFractionCellOutConePhi) ;
          
          fhFractionCellOutConePhiTrigEtaPhi  = new TH2F("hFractionCellOutConePhiTrigEtaPhi",
                                                         Form("Fraction of the isolation cone #it{R} =  %2.2f, out of cells #phi acceptance, in trigger #eta-#phi ",r),
                                                         netabins,etamin,etamax,nphibins,phimin,phimax);
          fhFractionCellOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
          fhFractionCellOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhFractionCellOutConePhiTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhFractionCellOutConePhiTrigEtaPhi) ;
          
          
          fhConeSumPtSubvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubvsConeSumPtTotPhiCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubNormvsConeSumPtTotPhiCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotPhiCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCell);
          
          fhConeSumPtSubvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubvsConeSumPtTotEtaCell",
                                                         Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                         nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCell);
          
          fhConeSumPtSubNormvsConeSumPtTotEtaCell = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaCell",
                                                             Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                             nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
          fhConeSumPtSubNormvsConeSumPtTotEtaCell->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCell);
        }
      }
    }
    
    // Track only histograms
    if(GetIsolationCut()->GetParticleTypeInCone()!=AliIsolationCut::kOnlyNeutral)
    {
      fhConeSumPtTrack  = new TH2F("hConePtSumTrack",
                                   Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrack->SetYTitle("#Sigma #it{p}_{T}");
      fhConeSumPtTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrack) ;

      fhConePtLeadTrack  = new TH2F("hConeLeadPtTrack",
                                   Form("Track leading in isolation cone for #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhConePtLeadTrack->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLeadTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLeadTrack) ;
      
      fhPtTrackInCone  = new TH2F("hPtTrackInCone",
                                  Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f",r),
                                  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInCone) ;
      
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaBandUETrack  = new TH2F("hConePtSumEtaBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Eta Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtEtaBandUETrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaBandUETrack) ;
        
        fhConeSumPtPhiBandUETrack  = new TH2F("hConePtSumPhiBandUETrack",
                                              "#Sigma track #it{p}_{T} in UE Phi Band",
                                              nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax*8);
        fhConeSumPtPhiBandUETrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiBandUETrack) ;
        
        
        fhConeSumPtEtaBandUETrackTrigEtaPhi  = new TH2F("hConePtSumEtaBandUETrackTrigEtaPhi",
                                                        "Trigger #eta vs #phi, #Sigma track #it{p}_{T} in UE Eta Band",
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUETrackTrigEtaPhi  = new TH2F("hConePtSumPhiBandUETrackTrigEtaPhi",
                                                        "Trigger #eta vs #phi, #Sigma track #it{p}_{T} in UE Phi Band",
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhi) ;
        
        fhEtaBandTrack  = new TH2F("hEtaBandTrack",
                                   Form("#eta vs #phi of tracks in #eta band isolation cone for #it{R} =  %2.2f",r),
                                   netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandTrack->SetXTitle("#eta");
        fhEtaBandTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaBandTrack) ;
        
        fhPhiBandTrack  = new TH2F("hPhiBandTrack",
                                   Form("#eta vs #phi of tracks in #phi band isolation cone for #it{R} =  %2.2f",r),
                                   netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandTrack->SetXTitle("#eta");
        fhPhiBandTrack->SetYTitle("#phi");
        outputContainer->Add(fhPhiBandTrack) ;
        
        fhConeSumPtEtaUESubTrack  = new TH2F("hConeSumPtEtaUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubTrack) ;
        
        fhConeSumPtPhiUESubTrack  = new TH2F("hConeSumPtPhiUESubTrack",
                                             Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                             nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubTrack) ;
        
        fhConeSumPtEtaUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackTrigEtaPhi",
                                                       Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                       netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubTrackTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubTrackTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackTrigEtaPhi",
                                                       Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                       netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubTrackTrigEtaPhi) ;
        
        fhFractionTrackOutConeEta  = new TH2F("hFractionTrackOutConeEta",
                                              Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance",r),
                                              nptbins,ptmin,ptmax,100,0,1);
        fhFractionTrackOutConeEta->SetYTitle("#it{fraction}");
        fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
        outputContainer->Add(fhFractionTrackOutConeEta) ;
        
        fhFractionTrackOutConeEtaTrigEtaPhi  = new TH2F("hFractionTrackOutConeEtaTrigEtaPhi",
                                                        Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance, in trigger #eta-#phi ",r),
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
        fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
        fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
        
        fhConeSumPtSubvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubvsConeSumPtTotPhiTrack",
                                                        Form("#Sigma #it{p}_{T} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                        nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiTrack);
        
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotPhiTrack",
                                                            Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #phi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiTrack);
        
        fhConeSumPtSubvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubvsConeSumPtTotEtaTrack",
                                                        Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                        nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaTrack);
        
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack = new TH2F("hConeSumPtSubNormvsConeSumPtTotEtaTrack",
                                                            Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",r),
                                                            nptsumbins,ptsummin,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
        fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaTrack);
        
        
        // UE in perpendicular cone
        fhPerpConeSumPt  = new TH2F("hPerpConePtSum",
                                    Form("#Sigma #it{p}_{T} in isolation cone at #pm 45 degree phi from trigger particle, #it{R} =  %2.2f",r),
                                    nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPerpConeSumPt->SetYTitle("#Sigma #it{p}_{T}");
        fhPerpConeSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeSumPt) ;
        
        fhPtInPerpCone  = new TH2F("hPtInPerpCone",
                                   Form("#it{p}_{T} in isolation cone at #pm 45 degree phi from trigger particle, #it{R} =  %2.2f",r),
                                   nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInPerpCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
        fhPtInPerpCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInPerpCone) ;
        
        fhEtaPhiTrack= new TH2F("hEtaPhiTrack",
                                Form("#eta vs #phi of all Tracks"),
                                netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiTrack->SetXTitle("#eta");
        fhEtaPhiTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiTrack) ;
        
        fhEtaPhiInConeTrack= new TH2F("hEtaPhiInConeTrack",
                                      Form("#eta vs #phi of Tracks in cone for #it{R} =  %2.2f",r),
                                      netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInConeTrack->SetXTitle("#eta");
        fhEtaPhiInConeTrack->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiInConeTrack) ;
        
        fhConeSumPtVSUETracksEtaBand  = new TH2F("hConeSumPtVSUETracksEtaBand",
                                                 Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in eta band for tracks (before normalization), R=%2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
        fhConeSumPtVSUETracksEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUETracksEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUETracksEtaBand);
        
        fhConeSumPtVSUETracksPhiBand  = new TH2F("hConeSumPtVSUETracksPhiBand",
                                                 Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in phi band for tracks (before normalization), R=%2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
        fhConeSumPtVSUETracksPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUETracksPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUETracksPhiBand);
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
        fhConeSumPtCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
        fhConeSumPtCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtCellvsTrack) ;
        
        fhConeSumPtCellTrack = new TH2F("hConePtSumCellTrack",
                                        Form("Track and Cell #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtCellTrack->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCellTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtCellTrack) ;
        
        fhConeSumPtCellTrackTrigEtaPhi  = new TH2F("hConePtSumCellTrackTrigEtaPhi",
                                                   Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",r),
                                                   netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtCellTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtCellTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtCellTrackTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtCellTrackTrigEtaPhi) ;
      }
      
      if(fFillUEBandSubtractHistograms)
      {
        fhConeSumPtEtaUESub  = new TH2F("hConeSumPtEtaUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESub->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESub) ;
        
        fhConeSumPtPhiUESub  = new TH2F("hConeSumPtPhiUESub",
                                        Form("#Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                        nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESub->SetYTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESub) ;
        
        fhConeSumPtEtaUESubTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrigEtaPhi",
                                                  Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                  netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrigEtaPhi",
                                                  Form("Trigger #eta vs #phi, #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                  netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubTrigEtaPhi) ;
        
        fhConeSumPtEtaUESubClustervsTrack   = new TH2F("hConePtSumEtaUESubClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
        
        fhConeSumPtPhiUESubClustervsTrack   = new TH2F("hConePhiUESubPtSumClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
        
        fhEtaBandClustervsTrack   = new TH2F("hEtaBandClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhEtaBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhEtaBandClustervsTrack) ;
        
        fhPhiBandClustervsTrack   = new TH2F("hPhiBandClustervsTrack",
                                             Form("Track vs Cluster #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                             nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
        fhPhiBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhPhiBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhPhiBandClustervsTrack) ;
        
        fhEtaBandNormClustervsTrack   = new TH2F("hEtaBandNormClustervsTrack",
                                                 Form("Track vs Cluster #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhEtaBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhEtaBandNormClustervsTrack) ;
        
        fhPhiBandNormClustervsTrack   = new TH2F("hPhiBandNormClustervsTrack",
                                                 Form("Track vs Cluster #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                                 nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhPhiBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhPhiBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhPhiBandNormClustervsTrack) ;
        
        fhConeSumPtEtaUESubClustervsTrack   = new TH2F("hConePtSumEtaUESubClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
        
        fhConeSumPtPhiUESubClustervsTrack   = new TH2F("hConePhiUESubPtSumClustervsTrack",
                                                       Form("Track vs Cluster #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                       2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
        fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T} cluster");
        fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T} track");
        outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
        
        if(fFillCellHistograms)
        {
          
          fhConeSumPtEtaUESubCellvsTrack   = new TH2F("hConePtSumEtaUESubCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub eta band in isolation cone for #it{R} =  %2.2f",r),
                                                      2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhConeSumPtEtaUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhConeSumPtEtaUESubCellvsTrack) ;
          
          fhConeSumPtPhiUESubCellvsTrack   = new TH2F("hConePhiUESubPtSumCellvsTrack",
                                                      Form("Track vs Cell #Sigma #it{p}_{T} UE sub phi band in isolation cone for #it{R} =  %2.2f",r),
                                                      2*nptsumbins,-ptsummax,ptsummax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhConeSumPtPhiUESubCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhConeSumPtPhiUESubCellvsTrack) ;
          
          fhEtaBandCellvsTrack   = new TH2F("hEtaBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhEtaBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhEtaBandCellvsTrack) ;
          
          fhPhiBandCellvsTrack   = new TH2F("hPhiBandCellvsTrack",
                                            Form("Track vs Cell #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                            nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
          fhPhiBandCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhPhiBandCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhPhiBandCellvsTrack) ;
          
          fhEtaBandNormCellvsTrack   = new TH2F("hEtaBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in Eta band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhEtaBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhEtaBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhEtaBandNormCellvsTrack) ;
          
          fhPhiBandNormCellvsTrack   = new TH2F("hPhiBandNormCellvsTrack",
                                                Form("Track vs Cell #Sigma #it{p}_{T} in Phi band in isolation cone for #it{R} =  %2.2f",r),
                                                nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
          fhPhiBandNormCellvsTrack->SetXTitle("#Sigma #it{p}_{T} cell");
          fhPhiBandNormCellvsTrack->SetYTitle("#Sigma #it{p}_{T} track");
          outputContainer->Add(fhPhiBandNormCellvsTrack) ;
          
          fhConeSumPtEtaUESubTrackCell  = new TH2F("hConeSumPtEtaUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtEtaUESubTrackCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCell) ;
          
          fhConeSumPtPhiUESubTrackCell  = new TH2F("hConeSumPtPhiUESubTrackCell",
                                                   Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                   nptbins,ptmin,ptmax,2*nptsumbins,-ptsummax,ptsummax);
          fhConeSumPtPhiUESubTrackCell->SetYTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackCell->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiUESubTrackCell) ;
          
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtEtaUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from eta band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtEtaUESubTrackCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtEtaUESubTrackCellTrigEtaPhi) ;
          
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi  = new TH2F("hConeSumPtPhiUESubTrackCellTrigEtaPhi",
                                                             Form("Trigger #eta vs #phi, Tracks #Sigma #it{p}_{T} after bkg subtraction from phi band in the isolation cone for #it{R} =  %2.2f",r),
                                                             netabins,etamin,etamax,nphibins,phimin,phimax);
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtPhiUESubTrackCellTrigEtaPhi->SetYTitle("#phi_{trigger} (rad)");
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
         Form("%s - d#phi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhi[iso]->SetXTitle("E_{cluster} (GeV)");
        
        fhTrackMatchedDEtaDPhi[iso]  = new TH2F
        (Form("hTrackMatchedDEtaDPhi%s",isoName[iso].Data()),
         Form("%s - d#eta vs d#phi of cluster-track, %s",isoTitle[iso].Data(),parTitle.Data()),
         nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDEtaDPhi[iso]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDEtaDPhi[iso]->SetXTitle("d#eta");
          
        outputContainer->Add(fhTrackMatchedDEta[iso]) ;
        outputContainer->Add(fhTrackMatchedDPhi[iso]) ;
        outputContainer->Add(fhTrackMatchedDEtaDPhi[iso]) ;
         
        if(IsDataMC())
        {
          for(int imc = 0; imc < fgkNmcTypes; imc++)
          {
              fhTrackMatchedDEtaMC[imc][iso] = new TH2F(Form("hTrackMatchedDEta%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
                                                        Form("%s - d#eta of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
                                                        nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
              fhTrackMatchedDEtaMC[imc][iso]->SetYTitle("d#eta");
              fhTrackMatchedDEtaMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
              
              fhTrackMatchedDPhiMC[imc][iso] = new TH2F(Form("hTrackMatchedDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
                                                        Form("%s - d#phi of cluster-track vs cluster energy, %s",isoTitle[iso].Data(),parTitle.Data()),
                                                        nptbins,ptmin,ptmax,nresetabins,resphimin,resphimax);
              fhTrackMatchedDPhiMC[imc][iso]->SetYTitle("d#phi");
              fhTrackMatchedDPhiMC[imc][iso]->SetXTitle("E_{cluster} (GeV)");
              
              fhTrackMatchedDEtaDPhiMC[imc][iso]  = new TH2F
              (Form("hTrackMatchedDEtaDPhi%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
               Form("%s - d#eta vs d#phi of cluster-track, %s",isoTitle[iso].Data(),parTitle.Data()),
               nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
              fhTrackMatchedDEtaDPhiMC[imc][iso]->SetYTitle("d#phi (rad)");
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
            fhPtLambda0MC[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%s",isoName[iso].Data(),mcPartName[imc].Data()),
                                                Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MC[imc][iso]->SetYTitle("#lambda_{0}^{2}");
            fhPtLambda0MC[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MC[imc][iso]) ;

            fhPtLambda0MCConv[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%sConv",isoName[iso].Data(),mcPartName[imc].Data()),
                                                Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
            fhPtLambda0MCConv[imc][iso]->SetYTitle("#lambda_{0}^{2}");
            fhPtLambda0MCConv[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
            outputContainer->Add( fhPtLambda0MCConv[imc][iso]) ;
            
            if(fFillOverlapHistograms)
            {
              fhPtLambda0MCWith1Overlap[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                  Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                  nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCWith1Overlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCWith1Overlap[imc][iso]) ;
              
              fhPtLambda0MCConvWith1Overlap[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                      Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCConvWith1Overlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCConvWith1Overlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCConvWith1Overlap[imc][iso]) ;

              fhPtLambda0MCWithNOverlap[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%s_NOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                              Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCWithNOverlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCWithNOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCWithNOverlap[imc][iso]) ;
              
              fhPtLambda0MCConvWithNOverlap[imc][iso]  = new TH2F(Form("hPtLambda0%s_MC%sConv_NOverlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                                  Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                                  nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhPtLambda0MCConvWithNOverlap[imc][iso]->SetYTitle("#lambda_{0}^{2}");
              fhPtLambda0MCConvWithNOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtLambda0MCConvWithNOverlap[imc][iso]) ;
              
              
              fhPtNOverlap[imc][iso]  = new TH2F(Form("hPtNOverlaps%s_MC%s_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                              Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
                                                              nptbins,ptmin,ptmax,10,0,10);
              fhPtNOverlap[imc][iso]->SetYTitle("#it{N} overlaps");
              fhPtNOverlap[imc][iso]->SetXTitle("#it{p}_{T}(GeV/#it{c})");
              outputContainer->Add( fhPtNOverlap[imc][iso]) ;
              
              fhPtNOverlapConv[imc][iso]  = new TH2F(Form("hPtNOverlaps%s_MC%sConv_1Overlap",isoName[iso].Data(),mcPartName[imc].Data()),
                                                                  Form("%s cluster : #it{p}_{T} vs #lambda_{0}: %s %s, from conversion, 1 overlap",isoTitle[iso].Data(),mcPartType[imc].Data(),parTitle.Data()),
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
           Form("%s cluster: #it{p}_{T} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhPtLambda0TRD[iso]->SetYTitle("#lambda_{0}^{2}");
          fhPtLambda0TRD[iso]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          outputContainer->Add(fhPtLambda0TRD[iso]) ;
          
          fhELambda0TRD[iso]  = new TH2F
          (Form("hELambda0TRD%s",isoName[iso].Data()),
           Form("%s cluster: #it{E} vs #lambda_{0}, SM behind TRD, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
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
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax1[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax1[iso]) ;
          
          fhELambda1LocMax1[iso]  = new TH2F
          (Form("hELambda1LocMax1%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=1, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax1[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax1[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax1[iso]) ;
          
          fhELambda0LocMax2[iso]  = new TH2F
          (Form("hELambda0LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMax2[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMax2[iso]) ;
          
          fhELambda1LocMax2[iso]  = new TH2F
          (Form("hELambda1LocMax2%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}=2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMax2[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMax2[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMax2[iso]) ;
          
          fhELambda0LocMaxN[iso]  = new TH2F
          ( Form("hELambda0LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{0}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda0LocMaxN[iso]->SetYTitle("#lambda_{0}^{2}");
          fhELambda0LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda0LocMaxN[iso]) ;
          
          fhELambda1LocMaxN[iso]  = new TH2F
          (Form("hELambda1LocMaxN%s",isoName[iso].Data()),
           Form("%s cluster (#eta) pairs: #it{E} vs #lambda_{1}, #it{NLM}>2, %s",isoTitle[iso].Data(),parTitle.Data()),nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
          fhELambda1LocMaxN[iso]->SetYTitle("#lambda_{1}^{2}");
          fhELambda1LocMaxN[iso]->SetXTitle("#it{E} (GeV)");
          outputContainer->Add(fhELambda1LocMaxN[iso]) ;
        } // NLM
      } // SS histo
      
      
      if(GetCalorimeter() == kEMCAL && fFillEMCALRegionSSHistograms)
      {
        for(Int_t ieta = 0; ieta < 4; ieta++) 
        {  
          for(Int_t iphi = 0; iphi < 3; iphi++) 
          {
//            fhLam0EMCALRegion[iso][ieta][iphi] = 
//            new TH2F(Form("hLam0_%s_eta%d_phi%d",isoName[iso].Data(),ieta,iphi),
//                     Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, region eta %d, phi %d",
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
//                       Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, region eta %d, phi %d, SM covered by TRD",
//                             isoTitle[iso].Data(),ieta,iphi),
//                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
//              fhLam0EMCALRegionTRD[iso][ieta][iphi]->SetYTitle("#lambda_{0}^{2}");
//              fhLam0EMCALRegionTRD[iso][ieta][iphi]->SetXTitle("#it{p}_{T} (GeV)");
//              outputContainer->Add(fhLam0EMCALRegionTRD[iso][ieta][iphi]) ;
//            } // TRD
            
            for(Int_t ism = 0; ism < GetCaloUtils()->GetNumberOfSuperModulesUsed(); ism++) 
            {
              fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism] = 
              new TH2F(Form("hLam0_%s_eta%d_phi%d_sm%d",isoName[iso].Data(),ieta,iphi,ism),
                       Form("%s, cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, sm %d, region eta %d, phi %d",
                            isoTitle[iso].Data(),ism,ieta,iphi),
                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
              fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]->SetYTitle("#lambda_{0}^{2}");
              fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]->SetXTitle("#it{p}_{T} (GeV)");
              outputContainer->Add(fhLam0EMCALRegionPerSM[iso][ieta][iphi][ism]) ;
            } // ism
          } // iphi 
        } // ieta
        
        Float_t ptLimit[] = {2,3,4,5,6,8,10,12};
        for(Int_t ipt = 0; ipt < 7; ipt++)
        {
          fhEtaPhiLam0BinPtBin[iso][ipt]  = new TH2F
          (Form("hEtaPhiLam0BinPtBin%d%s",ipt,isoName[iso].Data()),
           Form("%s, #eta vs #phi in #it{p}_{T}=[%2.1f,%2.1f] GeV/#it{c} and #lambda^{2}_{0}=[0.3,0.4]",
                isoTitle[iso].Data(),ptLimit[ipt],ptLimit[ipt+1]),
           netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiLam0BinPtBin[iso][ipt]->SetYTitle("#phi (rad)");
          fhEtaPhiLam0BinPtBin[iso][ipt]->SetXTitle("#eta");
          outputContainer->Add(fhEtaPhiLam0BinPtBin[iso][ipt]) ;
        }
      } // regions in EMCal
      
      if(IsDataMC())
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
        
        if(fFillSSHisto)
        {
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
          
//          if(GetCalorimeter() == kEMCAL && fFillEMCALRegionSSHistograms)
//          {
//            for(Int_t ieta = 0; ieta < 4; ieta++) 
//            {  
//              for(Int_t iphi = 0; iphi < 3; iphi++) 
//              {
//                for(Int_t iReg = 0; iReg < 6; iReg++) 
//                {
//                  fhLam0EMCALRegionMCConvRcut[iso][ieta][iphi][iReg] = 
//                  new TH2F(Form("hMCPhotonConversionLambda0%s_R%d_eta%d_phi%d",isoName[iso].Data(),iReg,ieta,iphi),
//                           Form("%s,cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s, region eta %d, phi %d",
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
//                             Form("%s,cluster from converted photon, #it{p}_{T} vs #lambda_{0}^{2}, conversion in %s, region eta %d, phi %d, SM covered by TRD",
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

          
        } // Shower shape histograms
      }
    } // control histograms for isolated and non isolated objects
    
    
    if(IsPileUpAnalysisOn())
    {
      fhPtTrackInConeOtherBC  = new TH2F("hPtTrackInConeOtherBC",
                                         Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0",r),
                                         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeOtherBC->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeOtherBC->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeOtherBC) ;
      
      fhPtTrackInConeOtherBCPileUpSPD  = new TH2F("hPtTrackInConeOtherBCPileUpSPD",
                                                  Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC!=0, pile-up from SPD",r),
                                                  nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeOtherBCPileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeOtherBCPileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeOtherBCPileUpSPD) ;
      
      fhPtTrackInConeBC0  = new TH2F("hPtTrackInConeBC0",
                                     Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0",r),
                                     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeBC0) ;
      
      fhPtTrackInConeVtxBC0  = new TH2F("hPtTrackInConeVtxBC0",
                                        Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0",r),
                                        nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeVtxBC0->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeVtxBC0->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeVtxBC0) ;
      
      
      fhPtTrackInConeBC0PileUpSPD  = new TH2F("hPtTrackInConeBC0PileUpSPD",
                                              Form("#it{p}_{T} of tracks in isolation cone for #it{R} =  %2.2f, TOF from BC==0, pile-up from SPD",r),
                                              nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInConeBC0PileUpSPD->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeBC0PileUpSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInConeBC0PileUpSPD) ;
      
      
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
                                   Form("primary photon %s : #phi vs #it{p}_{T}, %s",pptype[i].Data(),parTitle.Data()),
                                   nptbins,ptmin,ptmax,200,0.,TMath::TwoPi());
        fhPhiPrimMC[i]->SetYTitle("#phi");
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
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
        fhEtaPhiPtThresIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtThresIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtThresIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtThresIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtFrac_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiPtFracIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtFracIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtFracIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtFracIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiPtSumIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiPtSumIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiSumDensity_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
        fhEtaPhiSumDensityIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiSumDensityIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiSumDensityIso[icone][ipt]->SetYTitle("#phi");
        outputContainer->Add(fhEtaPhiSumDensityIso[icone][ipt]) ;
        
        snprintf(name, buffersize,"hEtaPhiFracPtSum_Cone_%d_Pt%d",icone,ipt);
        snprintf(title, buffersize,"Isolated candidate #eta:#phi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
        fhEtaPhiFracPtSumIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
        fhEtaPhiFracPtSumIso[icone][ipt]->SetXTitle("#eta");
        fhEtaPhiFracPtSumIso[icone][ipt]->SetYTitle("#phi");
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
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{th} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtThresholds[ipt]);
          fhEtaPhiPtThresDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtThresDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtThresDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiPtFrac_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiPtFracDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtFracDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtFracDecayIso[icone][ipt]) ;
          
          
          snprintf(name, buffersize,"hEtaPhiPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiPtSumDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiSumDensity_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for density #it{R} =  %2.2f and #it{p}_{T}^{sum} = %2.2f GeV/#it{c}",fConeSizes[icone],fSumPtThresholds[ipt]);
          fhEtaPhiSumDensityDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiSumDensityDecayIso[icone][ipt]->SetYTitle("#phi");
          outputContainer->Add(fhEtaPhiSumDensityDecayIso[icone][ipt]) ;
          
          snprintf(name, buffersize,"hEtaPhiFracPtSum_Decay_Cone_%d_Pt%d",icone,ipt);
          snprintf(title, buffersize,"Isolated decay candidate #eta:#phi distribution for FracPtSum #it{R} =  %2.2f and #it{p}_{T}^{fr} = %2.2f GeV/#it{c}",fConeSizes[icone],fPtFractions[ipt]);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]  = new TH2F(name, title,netabins,etamin,etamax,nphibins,phimin,phimax);
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetXTitle("#eta");
          fhEtaPhiFracPtSumDecayIso[icone][ipt]->SetYTitle("#phi");
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

  fNPtTrigBin = 6;
  fPtTrigBinLimit[ 0] =  8; fPtTrigBinLimit[ 1] = 10; fPtTrigBinLimit[ 2] = 12; fPtTrigBinLimit[ 3] = 14; fPtTrigBinLimit[ 4] = 16;
  fPtTrigBinLimit[ 5] = 20; fPtTrigBinLimit[ 6] = 25; ;
  for(Int_t ibin = fNPtTrigBin+1; ibin < 20; ibin++) fPtTrigBinLimit[ibin] = 00.0;
  
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
      Int_t  trackID   = TMath::Abs(track->GetID()) ;
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
    Float_t phi        = aod->Phi();
    Float_t eta        = aod->Eta();
    Float_t m02        = aod->GetM02();
    Int_t   iSM        = aod->GetSModNumber();
    
    AliDebug(1,Form("pt %1.1f, eta %1.1f, phi %1.1f, Isolated %d",pt, eta, phi, isolated));
    
    //
    // EMCAL SM regions
    //
    AliVCluster *cluster = 0;
    if ( GetCalorimeter() == kEMCAL && fFillEMCALRegionSSHistograms && fFillSSHisto)
    {
      // Get original cluster, needed to feed the subregion selection method
      
      Int_t iclus = -1;
      cluster = FindCluster(GetEMCALClusters(),aod->GetCaloLabel(0),iclus);

      Int_t etaRegion = -1, phiRegion = -1;
      
      GetCaloUtils()->GetEMCALSubregion(cluster,GetReader()->GetEMCALCells(),etaRegion,phiRegion);
      
      if(etaRegion >= 0 && etaRegion < 4 && phiRegion >=0 && phiRegion < 3) 
      {
//        fhLam0EMCALRegion[isolated][etaRegion][phiRegion]->Fill(pt, m02, GetEventWeight());
//        
//        if(GetFirstSMCoveredByTRD() >= 0 && iSM >= GetFirstSMCoveredByTRD()  )
//        {
//          fhLam0EMCALRegionTRD[isolated][etaRegion][phiRegion]->Fill(pt, m02, GetEventWeight());
//        }
        
        fhLam0EMCALRegionPerSM[isolated][etaRegion][phiRegion][iSM]->Fill(pt, m02, GetEventWeight());
      }
      
      if(m02 >=0.3 && m02 <= 0.4)
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
        if(ptbin >= 0) fhEtaPhiLam0BinPtBin[isolated][ptbin]->Fill(eta, phi, GetEventWeight());
      }
    }
    
    // Conversion radius
    if( IsDataMC()                                                                 &&        
        GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
        GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)     &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)        &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)        
      )
    {
      Int_t pdg  = 0, status  = 0, momLabel  = -1;
      Int_t pdgD = 0, statusD = 0, daugLabel = -1;
      Bool_t ok = kFALSE, okD = kFALSE;
      
      //fPrimaryMom = 
      GetMCAnalysisUtils()->GetMother(aod->GetLabel(),GetReader(), pdg, status, ok, momLabel);     
      //fMomentum = 
      GetMCAnalysisUtils()->GetDaughter(0,momLabel,GetReader(),pdgD, statusD, okD, daugLabel, fProdVertex);
      
      if(okD)
      {
        Float_t prodR = TMath::Sqrt(fProdVertex.X()*fProdVertex.X()+fProdVertex.Y()*fProdVertex.Y());
        
        //printf("Iso %d, Conversion: mom pdg %d (stat %d), 1st daugher %d (stat %d), mom label %d, org label %d, daugh label %d, prodR %f\n", isolated, pdg,status, pdgD, statusD, 
        //       momLabel, aod->GetLabel(),daugLabel,prodR);
        
        fhMCConversionVertex[isolated]->Fill(pt,prodR,GetEventWeight());
        if(GetCalorimeter() == kEMCAL && GetFirstSMCoveredByTRD() >= 0 && iSM >= GetFirstSMCoveredByTRD() )
          fhMCConversionVertexTRD[isolated]->Fill(pt,prodR,GetEventWeight());

        if(fFillSSHisto)
        {
          
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
//            if ( GetCalorimeter() == kEMCAL && fFillEMCALRegionSSHistograms )
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
      } // ok MC daughter
    } // ok MC conversion
    
    //---------------------------------------------------------------
    // Fill pt/sum pT distribution of particles in cone or in UE band
    //---------------------------------------------------------------
    
    Float_t coneptLeadCluster= 0;
    Float_t coneptLeadTrack  = 0;
    Float_t coneptsumCluster = 0;
    Float_t coneptsumTrack   = 0;
    Float_t coneptsumCell    = 0;
    Float_t etaBandptsumClusterNorm = 0;
    Float_t etaBandptsumTrackNorm   = 0;
    
    CalculateTrackSignalInCone   (aod,coneptsumTrack  , coneptLeadTrack  );
    CalculateCaloSignalInCone    (aod,coneptsumCluster, coneptLeadCluster);
    if(fFillCellHistograms)
      CalculateCaloCellSignalInCone(aod,coneptsumCell   );
    
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
    
    Float_t coneptLead = coneptLeadTrack;
    if(coneptLeadCluster > coneptLeadTrack) coneptLead = coneptLeadCluster;
    fhConePtLead->Fill(pt, coneptLead, GetEventWeight());
    
    AliDebug(1,Form("Particle %d Energy Sum in Isolation Cone %2.2f, Leading pT in cone %2.2f",
             iaod, coneptsumTrack+coneptsumCluster, coneptLead));
    
    //normalize phi/eta band per area unit
    if(fFillUEBandSubtractHistograms)
      CalculateNormalizeUEBandPerUnitArea(aod, coneptsumCluster, coneptsumCell, coneptsumTrack, etaBandptsumTrackNorm, etaBandptsumClusterNorm) ;
    
    //  printf("Histograms analysis : cluster pt = %f, etaBandTrack = %f, etaBandCluster = %f, isolation = %d\n",aod->Pt(),etaBandptsumTrackNorm,etaBandptsumClusterNorm,aod->IsIsolated());
    
        
    //---------------------------------------------------------------
    // Fill Shower shape and track matching histograms
    //---------------------------------------------------------------
    
    FillTrackMatchingShowerShapeControlHistograms(aod, coneptsumTrack+coneptsumCluster, coneptLead, mcIndex);
    
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
  Int_t    nprim     = 0;
  
  TParticle        * primStack = 0;
  AliAODMCParticle * primAOD   = 0;
  
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

  for(Int_t i=0 ; i < nprim; i++)
  {
    if(GetReader()->AcceptOnlyHIJINGLabels() && !GetReader()->IsHIJINGLabel(i)) continue ;
    
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
      
      // Protection against floating point exception
      if ( primStack->Energy() == TMath::Abs(primStack->Pz()) || 
          (primStack->Energy() - primStack->Pz()) < 1e-3      ||
          (primStack->Energy() + primStack->Pz()) < 0           )  continue ; 
      
      //printf("i %d, %s %d  %s %d \n",i, stack->Particle(i)->GetName(), stack->Particle(i)->GetPdgCode(),
      //       primStack->GetName(), primStack->GetPdgCode());
      
      //photonY   = 0.5*TMath::Log((prim->Energy()+prim->Pz())/(prim->Energy()-prim->Pz())) ;
      
      //Photon kinematics
      primStack->Momentum(fMomentum);
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
      
      // Protection against floating point exception
      if ( primAOD->E() == TMath::Abs(primAOD->Pz()) || 
          (primAOD->E() - primAOD->Pz()) < 1e-3      || 
          (primAOD->E() + primAOD->Pz()) < 0           )  continue ; 
      
      //photonY   = 0.5*TMath::Log((prim->Energy()+prim->Pz())/(prim->Energy()-prim->Pz())) ;
      
      //Photon kinematics
      fMomentum.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());
    }
    
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
      if(GetReader()->ReadStack()          &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primStack)) continue ;
      if(GetReader()->ReadAODMCParticles() &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fIsoDetector, primAOD  )) continue ;
    }
    
    // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
    if(!GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),fIsoDetector)) continue ;
    
    // Get tag of this particle photon from fragmentation, decay, prompt ...
    // Set the origin of the photon.
    tag = GetMCAnalysisUtils()->CheckOrigin(i,GetReader(),fIsoDetector);
    
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
      Int_t ndaugh = GetMCAnalysisUtils()->GetNDaughters(i,GetReader(), okpi0);
     // printf("OK pi0 %d, ndaugh %d\n",okpi0,ndaugh);
      Int_t d1Pdg = 0, d1Status = 0; Bool_t ok1 = kFALSE;
      Int_t d2Pdg = 0, d2Status = 0; Bool_t ok2 = kFALSE;

      if ( ndaugh > 0 ) fMomDaugh1 = GetMCAnalysisUtils()->GetDaughter(0,i,GetReader(),d1Pdg, d1Status,ok1, pi0d1Label,fProdVertex);
      if ( ndaugh > 1 ) fMomDaugh2 = GetMCAnalysisUtils()->GetDaughter(1,i,GetReader(),d2Pdg, d2Status,ok2, pi0d2Label,fProdVertex);
      
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
    TParticle        * mcisopStack = 0;
    AliAODMCParticle * mcisopAOD   = 0;
    Int_t partInConeStatus = -1, partInConeMother = -1;
    Double_t partInConePt = 0, partInConeE = 0, partInConeEta = 0, partInConePhi = 0;
    Int_t partInConeCharge = 0, npart = 0;
    for(Int_t ip = 0; ip < nprim ; ip++)
    {
      if(ip==i) continue;
      
      if( (pdg==111 || pdg==221) && ( ip == pi0d1Label || ip == pi0d2Label ) )
      {
        //printf("Do not count pi0 decays in cone when isolating pi0 \n");
        continue;
      }
      
      if( GetReader()->ReadStack() )
      {
        mcisopStack = static_cast<TParticle*>(stack->Particle(ip));
        if( !mcisopStack ) continue;
        partInConeStatus = mcisopStack->GetStatusCode();
        
        // Consider only final state particles, but this depends on generator,
        // status 1 is the usual one, in case of not being ok, leave the possibility
        // to not consider this.
        if( partInConeStatus != 1 &&
            GetMCAnalysisUtils()->GetMCGenerator()!= AliMCAnalysisUtils::kBoxLike ) continue ;
        
        partInConeMother = mcisopStack->GetMother(0);
        partInConePt     = mcisopStack->Pt();
        partInConeE      = mcisopStack->Energy();
        partInConeEta    = mcisopStack->Eta();
        partInConePhi    = mcisopStack->Phi();
        partInConeCharge = TMath::Abs((Int_t) TDatabasePDG::Instance()->GetParticle(mcisopStack->GetPdgCode())->Charge());
        mcisopStack->Momentum(fMomIso);
      }
      else
      {
        mcisopAOD   = (AliAODMCParticle *) mcparticles->At(ip);
        if( !mcisopAOD )   continue;
        
        partInConeStatus = mcisopAOD->GetStatus();
        // Consider only final state particles, but this depends on generator,
        // status 1 is the usual one, in case of not being ok, leave the possibility
        // to not consider this.
        if( partInConeStatus != 1 &&
            GetMCAnalysisUtils()->GetMCGenerator() != AliMCAnalysisUtils::kBoxLike ) continue ;
        
        partInConeMother = mcisopAOD->GetMother();
        partInConePt     = mcisopAOD->Pt();
        partInConeE      = mcisopAOD->E();
        partInConeEta    = mcisopAOD->Eta();
        partInConePhi    = mcisopAOD->Phi();
        partInConeCharge = TMath::Abs(mcisopAOD->Charge());
        fMomIso.SetPxPyPzE(mcisopAOD->Px(),mcisopAOD->Py(),mcisopAOD->Pz(),mcisopAOD->E());
      }
      
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
            if(GetReader()->ReadStack()          &&
               !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), mcisopStack)) continue ;
            if(GetReader()->ReadAODMCParticles() &&
               !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(GetCalorimeter(), mcisopAOD  )) continue ;
          }
        }
      }
      //
      
      dR = GetIsolationCut()->Radius(photonEta, photonPhi, partInConeEta, partInConePhi);
      
      if(dR > GetIsolationCut()->GetConeSize())
        continue;
      
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
        if(mcIndex == kmcPrimPi0Decay) GetMCAnalysisUtils()->GetMotherWithPDG(i,111,GetReader(),okpi0, pi0label);
        else                           GetMCAnalysisUtils()->GetMotherWithPDG(i,221,GetReader(),okpi0, pi0label);
        
        if(okpi0)
        {
          ndaugh = GetMCAnalysisUtils()->GetNDaughters(pi0label,GetReader(), okpi0);
          if(ndaugh==2)
          {
            Int_t d1Pdg = 0, d1Status = 0;
            fMomDaugh1 = GetMCAnalysisUtils()->GetDaughter(0,pi0label,GetReader(),d1Pdg, d1Status,ok1, d1Label,fProdVertex);
            Int_t d2Pdg = 0, d2Status = 0;
            fMomDaugh2 = GetMCAnalysisUtils()->GetDaughter(1,pi0label,GetReader(),d2Pdg, d2Status,ok2, d2Label,fProdVertex);
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
/// Set the detrimeter for the analysis.
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
/// Set the detrimeter for the analysis.
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

