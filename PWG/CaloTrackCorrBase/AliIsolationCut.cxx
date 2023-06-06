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
#include <TCustomBinning.h>

// --- AliRoot system ---
#include "AliCaloTrackParticleCorrelation.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliMixedEvent.h"
#include "AliLog.h"
#include "AliRhoParameter.h"

// --- CaloTrackCorrelations --- 
#include "AliCaloTrackReader.h"
#include "AliCalorimeterUtils.h"
#include "AliCaloPID.h"
#include "AliFiducialCut.h"
#include "AliHistogramRanges.h"

#include "AliIsolationCut.h"

/// \cond CLASSIMP
ClassImp(AliIsolationCut) ;
/// \endcond

//____________________________________
/// Default constructor. Initialize parameters
//____________________________________
AliIsolationCut::AliIsolationCut() :
TObject(),
fFillHistograms(0),        fFillEtaPhiHistograms(0), fEtaPhiHistogramsMinPt(0), 
fFillLeadingHistograms(0), fFillLeadingVsUESubHisto(0), fFillHighMultHistograms(0),
fMakeConeExcessCorr(0),    fFillFractionExcessHistograms(0),
fConeSize(0.),       fConeSizeBandGap(0.),
fPtThreshold(0.),    fPtThresholdMax(10000.),
fSumPtThreshold(0.), fSumPtThresholdMax(10000.),    fSumPtThresholdGap(0.),
fPtFraction(0.),     fICMethod(0),                  fPartInCone(0),
fFracIsThresh(1),    fIsTMClusterInConeRejected(1), fDistMinToTrigger(-1.),
fUseLeadingPtUEFactor(0), fLeadingPtUEFactor(10000),
fUseMaxPtUE(0),      fMaxPtUE(1000),
fJetRhoTaskName(""),
fDebug(0),           fMomentum(),                   fTrackVector(),
fDMCEtaGap(-1),      fDMCPhiMin(-1),                fDMCPhiMax(-1),
fEMCEtaSize(-1),     fEMCPhiMin(-1),                fEMCPhiMax(-1),
fTPCEtaSize(-1),     fTPCPhiSize(-1),
// Histograms
fHistoRanges(0),                            fNCentBins(0),
fhPtInCone(0),       
fhPtClusterInCone(0),                       fhPtTrackInCone(0),
fhPtInConeCent(0),
fhPtClusterInConeCent(0),                   fhPtTrackInConeCent(0),
fhConeSumPt(0),      
fhConeSumPtCluster(0),                      fhConeSumPtClusterMatched(0),
fhConeSumPtClusterMatchedFraction(0),       fhConeSumPtClusterMatchedFraction3D(0),
fhConeSumPtTrack(0),
fhConeSumPtClustervsTrack(0),               fhConeSumPtClusterTrackFrac(0),            
fhConeSumPtTrigEtaPhi(0),
fhConeSumPtTrackTrigEtaPhi(0),              fhConeSumPtClusterTrigEtaPhi(0),
fhConeSumPtUESub(0),      
fhConeSumPtUESubCluster(0),                 fhConeSumPtUESubTrack(0),
fhConeSumPtUESubClusterCutMax(0),           fhConeSumPtUESubTrackCutMax(0),
fhConeSumPtUESubClusterCutLeadFactor(0),    fhConeSumPtUESubTrackCutLeadFactor(0),
fhConeSumPtUESubClustervsTrack(0),          fhConeSumPtUESubClusterTrackFrac(0),
fhConeSumPtUESubTrigEtaPhi(0),
fhConeSumPtUESubTrackTrigEtaPhi(0),         fhConeSumPtUESubClusterTrigEtaPhi(0),
fhConePtLead(0),     
fhConePtLeadCluster(0),                     fhConePtLeadTrack(0),
fhConePtLeadClustervsTrack(0),              fhConePtLeadClusterTrackFrac(0),
fhEtaPhiCluster(0),                         fhEtaPhiTrack(0),
fhEtaPhiInConeCluster(0),                   fhDeltaEtaPhiInConeCluster(0),   fhTriggerEtaPhiInConeClusterPt(0),
fhEtaPhiInConeTrack(0),                     fhDeltaEtaPhiInConeTrack(0),     fhTriggerEtaPhiInConeTrackPt(0),
fhPtInPerpCone(0),                          fhPerpConeSumPt(0),
fhPerpConeRho(0),                           fhPerpConeRhoCutMax(0),          fhPerpConeRhoCutLeadFactor(0),
fhEtaPhiInPerpCone(0),                      fhConeSumPtVSPerpCone(0),        fhPerpConeSumPtTrackSubVsNoSub(0),
fhPerpConeSumPtTrigEtaPhi(0),
fhDeltaEtaPhiInPerpCone(0),                 fhTriggerEtaPhiInPerpConeTrackPt(0),
fhPerpCone1SumPt(0),                        fhPerpCone1SumPtUESub(0),
fhJetRhoSumPt(0),                           fhJetRho(0),
fhConeSumPtVSJetRho(0),                     fhJetRhoSumPtTrackSubVsNoSub(0),
fhJetRhoSumPtTrigEtaPhi(0),
fhEtaBandClusterPt(0),                      fhPhiBandClusterPt(0),
fhEtaBandTrackPt(0),                        fhPhiBandTrackPt(0)           , fhPerpBandTrackPt(0),
fhConeSumPtEtaBandUECluster(0),             fhConeSumPtPhiBandUECluster(0),
fhConeSumPtEtaBandUETrack(0),               fhConeSumPtPhiBandUETrack(0)  , fhConeSumPtPerpBandUETrack(0),
fhEtaBandClusterEtaPhi(0),                  fhPhiBandClusterEtaPhi(0),
fhEtaBandTrackEtaPhi(0),                    fhPhiBandTrackEtaPhi(0)       , fhPerpBandTrackEtaPhi(0),
fhEtaBandClusterDeltaEtaPhi(0),             fhPhiBandClusterDeltaEtaPhi(0),
fhEtaBandTrackDeltaEtaPhi(0),               fhPhiBandTrackDeltaEtaPhi(0)  , fhPerpBandTrackDeltaEtaPhi(0),
fhEtaBandClusterPtTriggerEtaPhi(0),         fhPhiBandClusterPtTriggerEtaPhi(0),
fhEtaBandTrackPtTriggerEtaPhi(0),           fhPhiBandTrackPtTriggerEtaPhi(0),fhPerpBandTrackPtTriggerEtaPhi(0),
fhConeSumPtEtaBandUEClusterTrigEtaPhi(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhi(0),
fhConeSumPtEtaBandUETrackTrigEtaPhi(0),     fhConeSumPtPhiBandUETrackTrigEtaPhi(0), fhConeSumPtPerpBandUETrackTrigEtaPhi(0),
fhConeSumPtEtaBandUEClusterTrigEtaPhiCent(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhiCent(0),
fhConeSumPtEtaBandUETrackTrigEtaPhiCent(0),     fhConeSumPtPhiBandUETrackTrigEtaPhiCent(0), fhConeSumPtPerpBandUETrackTrigEtaPhiCent(0),
fhConeSumPtVSUETrackBand(0),                fhConeSumPtVSUEClusterBand(0),
fhConeSumPtUEBandNormCluster(0),            fhConeSumPtUEBandNormTrack(0), 
fhConeRhoUEBandCluster(0),                  fhConeRhoUEBandTrack(0),
fhConeRhoUEBandClusterCutMax(0),            fhConeRhoUEBandTrackCutMax(0),
fhConeRhoUEBandClusterCutLeadFactor(0),     fhConeRhoUEBandTrackCutLeadFactor(0),
fhFractionTrackOutConeEta(0),               fhFractionTrackOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConeEta(0),             fhFractionClusterOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConePhi(0),             fhFractionClusterOutConePhiTrigEtaPhi(0),
fhFractionClusterOutConeEtaPhi(0),          fhFractionClusterOutConeEtaPhiTrigEtaPhi(0),
fhConeSumPtUEBandSubClustervsTrack(0),       
fhBandClustervsTrack(0),                    fhBandNormClustervsTrack(0),             
fhConeSumPtTrackSubVsNoSub(0),              fhConeSumPtClusterSubVsNoSub(0),
fhConeSumPtCent(0),                         
fhConeSumPtClusterCent(0),                  fhConeSumPtTrackCent(0),
fhConeSumPtClustervsTrackCent(0),           fhConeSumPtClusterTrackFracCent(0),
fhConeSumPtTrigEtaPhiCent(0),
fhConeSumPtTrackTrigEtaPhiCent(0),          fhConeSumPtClusterTrigEtaPhiCent(0),
fhConeSumPtUESubCent(0),
fhConeSumPtUESubClusterCent (0),            fhConeSumPtUESubTrackCent(0),
fhConeSumPtUESubClusterCutMaxCent (0),      fhConeSumPtUESubTrackCutMaxCent(0),
fhConeSumPtUESubClusterCutLeadFactorCent(0),fhConeSumPtUESubTrackCutLeadFactorCent(0),
fhConeSumPtUESubClustervsTrackCent(0),      fhConeSumPtUESubClusterTrackFracCent(0),
fhConeSumPtUESubTrigEtaPhiCent(0),
fhConeSumPtUESubTrackTrigEtaPhiCent(0),     fhConeSumPtUESubClusterTrigEtaPhiCent(0),
fhPerpConeSumPtCent (0),
fhPerpConeRhoCent (0),                      fhPerpConeRhoCutMaxCent(0),        fhPerpConeRhoCutLeadFactorCent(0),
fhConeSumPtVSPerpConeCent(0),               fhPerpConeSumPtTrackSubVsNoSubCent(0),
fhPtInPerpConeCent(0),                      fhPerpConeSumPtTrigEtaPhiCent(0),
fhPerpCone1SumPtCent(0),                    fhPerpCone1SumPtUESubCent(0),
fhJetRhoSumPtCent(0),                       fhJetRhoCent(0),
fhConeSumPtVSJetRhoCent(0),                 fhJetRhoSumPtTrigEtaPhiCent(0),
fhConeSumPtUEBandNormClusterCent(0),        fhConeSumPtUEBandNormTrackCent(0),
fhConeRhoUEBandClusterCent(0),              fhConeRhoUEBandTrackCent(0),
fhConeRhoUEBandClusterCutMaxCent(0),        fhConeRhoUEBandTrackCutMaxCent(0),
fhConeRhoUEBandClusterCutLeadFactorCent(0), fhConeRhoUEBandTrackCutLeadFactorCent(0),
fhConeSumPtEtaBandUEClusterCent(0),         fhConeSumPtPhiBandUEClusterCent(0), 
fhConeSumPtEtaBandUETrackCent(0),           fhConeSumPtPhiBandUETrackCent(0),  fhConeSumPtPerpBandUETrackCent(0),
fhEtaBandClusterPtCent(0),                  fhPhiBandClusterPtCent(0),
fhEtaBandTrackPtCent(0),                    fhPhiBandTrackPtCent(0),           fhPerpBandTrackPtCent(0),
fhConeSumPtUEBandSubClustervsTrackCent(0),
fhBandClustervsTrackCent(0),                fhBandNormClustervsTrackCent(0),
fhConeSumPtVSUETrackBandCent(0),            fhConeSumPtVSUEClusterBandCent(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt(0)             , fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead(0)    , fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead(0),
fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt(0)           , fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt(0)         , fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent(0)         , fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent(0), fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent(0),
fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent(0)       , fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent(0),
fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent(0)     , fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent(0)
{
  InitParameters();
}

//_________________________________________________________________________________________________________________________________
/// Get the pt sum of the clusters inside the cone, the leading cluster pT and number of clusters 
///
/// \param aodParticle: Kinematics and + of candidate particle for isolation.
/// \param reader: pointer to AliCaloTrackReader. Needed to access event info.
/// \param bFillAOD: Indicate if particles in cone must be added to AOD particle object.
/// \param useRefs: Get the list of tracks or clusters in cone from references
/// \param aodArrayRefName: Name of array where list of tracks/clusters in cone is stored.
/// \param bgCls: List of clusters from mixed event background pool AliAnaParticleHadronCorrelation.
/// \param calorimeter: Which input trigger calorimeter used
/// \param pid: pointer to AliCaloPID. Needed to reject matched clusters in isolation cone.
/// \param nPart: number of tracks/clusters above threshold in cone, output.
/// \param nfrac: 1 if fraction pT cluster-track / pT trigger in cone avobe threshold, output.
/// \param coneptsumCluster: total momentum energy in cone (track+cluster), output.
/// \param coneptLeadCluster: momentum of leading cluster or track in cone, output.
/// \param etaBandPtSumCluster: sum of clusters in eta band, same phi as candidate
/// \param etaBandPtLeadCluster: momentum of leading cluster in eta band
/// \param etaBandPtSumClusterCutMax: sum of clusters in eta band, same phi as candidate, pt < max pt
/// \param etaBandPtSumClusterCutLeadFactor: sum of clusters in eta band, same phi as candidate, pt < lead in cone pt
/// \param phiBandPtSumCluster: sum of clusters in phi band, same eta as candidate (not useful need to restrict it)
/// \param phiBandPtLeadCluster: momentum of leading cluster in phi band
/// \param phiBandPtSumClusterCutMax: sum of clustersin phi band, same eta as candidate (not useful need to restrict it), pt <max pt
/// \param phiBandPtSumClusterCutLeadFactor: sum of clustersin phi band, same eta as candidate (not useful need to restrict it), pt < lead in cone pt
/// \param histoWeight: Histograms weight (event, pt depedent)
/// \param centrality:Centrality percentile
/// \param cenBin: Assigned centrality bin with defined range elsewhere
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateCaloSignalInCone 
(
 AliCaloTrackParticleCorrelation * pCandidate, AliCaloTrackReader * reader,
 Bool_t    bFillAOD           , Bool_t    useRefs, 
 TString   aodArrayRefName    , TObjArray  * bgCls,
 Int_t     calorimeter        , AliCaloPID * pid,
 Int_t   & nPart              , Int_t   & nfrac,
 Float_t & coneptsumCluster   , Float_t & coneptLeadCluster   ,
 Float_t & etaBandPtSumCluster, Float_t & etaBandPtLeadCluster, Float_t & etaBandPtSumClusterCutMax, Float_t & etaBandPtSumClusterCutLeadFactor,
 Float_t & phiBandPtSumCluster, Float_t & phiBandPtLeadCluster, Float_t & phiBandPtSumClusterCutMax, Float_t & phiBandPtSumClusterCutLeadFactor,
 Double_t  histoWeight        , 
 Float_t   centrality         , Int_t     cenBin
) 
{
  if ( fPartInCone == kOnlyCharged ) return ;
  
  // Get the array with clusters
  //
  TObjArray * plNe = 0x0; 
  
  if ( bgCls )
  {
    plNe = bgCls ;
  }
  else if ( !useRefs )
  {
    if      (calorimeter == AliFiducialCut::kPHOS )
      plNe = reader->GetPHOSClusters();
    else if (calorimeter == AliFiducialCut::kEMCAL)
      plNe = reader->GetEMCALClusters();
  }
  else
  {
    plNe = pCandidate->GetObjArray(aodArrayRefName+"Clusters"); 
  }
  
  if ( !plNe ) 
  {
    // In case of no entries, add 0 to histogram to keep number of trigger entries
    if ( fFillHistograms )
    {
      if( fICMethod != kSumBkgSubIC )
      {
        if ( fFillHighMultHistograms )
          fhConeSumPtClusterCent->Fill(pCandidate->Pt(), 0., centrality, histoWeight);
        else
          fhConeSumPtCluster->Fill(pCandidate->Pt(), 0., histoWeight);
      }
      
      if ( !fFillHighMultHistograms && fFillLeadingHistograms )
        fhConePtLeadCluster->Fill(pCandidate->Pt(), 0., histoWeight);
    }
    
    if ( coneptLeadCluster > 0  || coneptsumCluster > 0 ) 
      AliError(Form("No ref clusters!!! sum %f, lead %f",coneptsumCluster,coneptLeadCluster));
    
    return ;
  }
  
  // Init parameters
  //
  Float_t ptC   = pCandidate->Pt() ;
  Float_t phiC  = pCandidate->Phi() ;
  if ( phiC < 0 ) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  Float_t pt     = -100. ;
  Float_t eta    = -100. ;
  Float_t phi    = -100. ;
  Float_t rad    = -100. ;
  
  TObjArray * refclusters  = 0x0;
  Int_t       nclusterrefs = 0;
  
  //
  // Get the clusters in the cone
  //
  //printf("Loop calo\n");
  Float_t coneptsumClusterMatched = 0;
  for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ )
  {
    AliVCluster * calo = dynamic_cast<AliVCluster *>(plNe->At(ipr)) ;
    
    if ( calo )
    {
      // Get the index where the cluster comes, to retrieve the corresponding vertex
      Int_t evtIndex = 0 ;
      if ( reader->GetMixedEvent() )
        evtIndex=reader->GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ;
      
      // Do not count the candidate (photon or pi0) or the daughters of the candidate
      if ( calo->GetID() == pCandidate->GetCaloLabel(0) ||
           calo->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;
      
      // Assume that come from vertex in straight line
      calo->GetMomentum(fMomentum,reader->GetVertex(evtIndex)) ;
      
      pt  = fMomentum.Pt()  ;
      eta = fMomentum.Eta() ;
      phi = fMomentum.Phi() ;
    }
    else
    {// Mixed event stored in AliCaloTrackParticles
      AliCaloTrackParticle * calomix = dynamic_cast<AliCaloTrackParticle*>(plNe->At(ipr)) ;
      
      if ( !calomix )
      {
        AliWarning("Wrong calo data type, continue");
        continue;
      }
      
      pt  = calomix->Pt();
      eta = calomix->Eta();
      phi = calomix->Phi() ;
    }
    
    // ** Calculate distance between candidate and tracks **
    
    if( phi < 0 ) phi+=TMath::TwoPi();
    
    rad = Radius(etaC, phiC, eta, phi);
    
    // ** Exclude clusters too close to the candidate, inactive by default **
    
    if ( rad < fDistMinToTrigger ) continue ;
    
    //-------------------------------------------------------------
    // ** For the isolated particle **
    //-------------------------------------------------------------
    //
    //      // Only loop the particle at the same side of candidate
    //      if ( TMath::Abs(phi-phiC)>TMath::PiOver2() ) continue ;
    //
    //      // If at the same side has particle larger than candidate,
    //      // then candidate can not be the leading, skip such events
    //      if(pt > ptC)
    //      {
    //        n         = -1;
    //        nfrac     = -1;
    //        coneptsumCluster = -1;
    //        isolated  = kFALSE;
    //
    //        pCandidate->SetLeadingParticle(kFALSE);
    //
    //        if(bFillAOD)
    //        {
    //          if(reftracks)
    //          {
    //            reftracks  ->Clear();
    //            delete reftracks;
    //          }
    //
    //          if(refclusters)
    //          {
    //            refclusters->Clear();
    //            delete refclusters;
    //          }
    //        }
    //        return ;
    //      }
    //
    //AliDebug(2,Form("\t Cluster %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f", 
    //                ipr,pt,eta,phi,rad));
    
    //
    // Select calorimeter clusters inside the isolation radius 
    //
    if ( rad > fConeSize ) continue ;
    
    AliDebug(2,"Inside candidate cone");
    
    // Skip matched clusters with tracks in case of neutral+charged analysis
    if ( calo )
    {
      Bool_t bRes = kFALSE, bEoP = kFALSE;
      Bool_t matched = pid->IsTrackMatched(calo, reader->GetCaloUtils(),
                                           reader->GetInputEvent(),
                                           bEoP,bRes);

      Float_t dPhi = calo->GetTrackDx();

      AliVTrack *track = reader->GetCaloUtils()->GetMatchedTrack(calo, reader->GetInputEvent());

      if ( TMath::Abs(dPhi) < 999 && track )
      {
        Float_t ptTrack = track->Pt();
        //printf("MATCHED track pt %2.2f, cluster pt %2.2f, min %2.2f, max %2.2f\n",
        //       ptTrack,pt, reader->GetCTSPtMin(), reader->GetCTSPtMax() );
        if ( reader->GetCTSPtMin()  < ptTrack && reader->GetCTSPtMax()  > ptTrack )
          coneptsumClusterMatched+=pt;
      }
      if ( fIsTMClusterInConeRejected && fPartInCone == kNeutralAndCharged && matched ) continue ;
    }

    if ( bFillAOD )
    {
      nclusterrefs++;
      
      if ( nclusterrefs == 1 )
      {
        refclusters = new TObjArray(0);
        //refclusters->SetName(Form("Clusters%s",aodArrayRefName.Data()));
        TString tempo(aodArrayRefName)  ;
        tempo += "Clusters" ;
        refclusters->SetName(tempo);
        refclusters->SetOwner(kFALSE);
      }
      
      refclusters->Add(calo);
    }
    
    coneptsumCluster+=pt;
    
    if ( coneptLeadCluster < pt ) coneptLeadCluster = pt;
    
    if ( fFillHistograms )
    {
      if ( !fFillHighMultHistograms )
      {
        fhPtInCone->Fill(ptC, pt, histoWeight);

        if ( fPartInCone == kNeutralAndCharged )
          fhPtClusterInCone->Fill(ptC, pt, histoWeight);
      }
      else
      {
        fhPtInConeCent->Fill(ptC, pt, centrality, histoWeight);

        if ( fPartInCone == kNeutralAndCharged )
          fhPtClusterInConeCent->Fill(ptC, pt, centrality, histoWeight);
      }
      
      if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
      {
        fhEtaPhiInConeCluster->Fill(eta, phi, histoWeight);
        fhDeltaEtaPhiInConeCluster->Fill(etaC-eta, phiC-phi, histoWeight);
        fhTriggerEtaPhiInConeClusterPt->Fill(etaC, phiC, pt, histoWeight);
      }
    }
    
    // *Before*, count particles in cone
    //
    if ( pt > fPtThreshold && pt < fPtThresholdMax )  nPart++;
    
    //if fPtFraction*ptC<fPtThreshold then consider the fPtThreshold directly
    if ( fFracIsThresh )
    {
      if ( fPtFraction*ptC < fPtThreshold )
      {
        if ( pt > fPtThreshold )    nfrac++ ;
      }
      else
      {
        if ( pt > fPtFraction*ptC ) nfrac++;
      }
    }
    else
    {
      if ( pt > fPtFraction*ptC ) nfrac++;
    }
    
  } // in cone cluster loop

  //
  // Get the UE clusters out of the cone
  //

  for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ )
  {
    AliVCluster * calo = dynamic_cast<AliVCluster *>(plNe->At(ipr)) ;

    if ( calo )
    {
      // Get the index where the cluster comes, to retrieve the corresponding vertex
      Int_t evtIndex = 0 ;
      if ( reader->GetMixedEvent() )
        evtIndex=reader->GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ;

      // Do not count the candidate (photon or pi0) or the daughters of the candidate
      if ( calo->GetID() == pCandidate->GetCaloLabel(0) ||
           calo->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;

      // Skip matched clusters with tracks in case of neutral+charged analysis
      if ( fIsTMClusterInConeRejected )
      {
        Bool_t bRes = kFALSE, bEoP = kFALSE;
        Bool_t matched = pid->IsTrackMatched(calo, reader->GetCaloUtils(),
                                             reader->GetInputEvent(),
                                             bEoP,bRes);
        if ( fPartInCone == kNeutralAndCharged && matched ) continue ;
      }

      // Assume that come from vertex in straight line
      calo->GetMomentum(fMomentum,reader->GetVertex(evtIndex)) ;

      pt  = fMomentum.Pt()  ;
      eta = fMomentum.Eta() ;
      phi = fMomentum.Phi() ;
    }
    else
    {// Mixed event stored in AliCaloTrackParticles
      AliCaloTrackParticle * calomix = dynamic_cast<AliCaloTrackParticle*>(plNe->At(ipr)) ;

      if ( !calomix )
      {
        AliWarning("Wrong calo data type, continue");
        continue;
      }

      pt  = calomix->Pt();
      eta = calomix->Eta();
      phi = calomix->Phi() ;
    }

    // ** Calculate distance between candidate and tracks **

    if( phi < 0 ) phi+=TMath::TwoPi();

    rad = Radius(etaC, phiC, eta, phi);

    // ** For the background out of cone **

    // Angles between trigger and track
    Double_t dEta = etaC - eta;
    Double_t dPhi = phiC - phi;

    // Shift phi angle when trigger is close to 0 or 360
    if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
    if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();

    if ( fFillHistograms && fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
      fhEtaPhiCluster->Fill(eta, phi, histoWeight);

    if ( fICMethod > kSumBkgSubIC && rad > fConeSize+fConeSizeBandGap )
    {
      // phi band
      //
      Bool_t takeIt = kTRUE;

      // Look only at 90 degrees with respect candidate, avoid opposite side jet
      // In case of presence of DCal or combination EMCal+DCal+PHOS,
      // anyway Phi band should not be used for DCal/PHOS, there is no space.
      if ( TMath::Abs(dPhi) > TMath::PiOver2() ) takeIt = kFALSE ;

      // Within eta cone size
      if ( TMath::Abs(dEta) < (fConeSize+fConeSizeBandGap) && takeIt )
      {
        phiBandPtSumCluster += pt;
        if ( pt < fMaxPtUE )
          phiBandPtSumClusterCutMax += pt;
        if ( pt < coneptLeadCluster*fLeadingPtUEFactor )
          phiBandPtSumClusterCutLeadFactor += pt;

        if ( fFillHistograms && fICMethod == kSumBkgSubPhiBandIC )
        {
          if ( fFillHighMultHistograms )
            fhPhiBandClusterPtCent->Fill(ptC, pt, centrality, histoWeight);
          else
            fhPhiBandClusterPt->Fill(ptC, pt, histoWeight);

          if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
          {
            fhPhiBandClusterEtaPhi->Fill(eta, phi, histoWeight);
            fhPhiBandClusterDeltaEtaPhi->Fill(etaC-eta, phiC-phi, histoWeight);
            fhPhiBandClusterPtTriggerEtaPhi->Fill(etaC, phiC, pt, histoWeight);
          }
        }
      } // phi band

      // eta band
      //
      //takeIt = kTRUE;

      // Within phi cone size
      if ( TMath::Abs(dPhi) < (fConeSize+fConeSizeBandGap) )// && takeIt )
      {
        etaBandPtSumCluster += pt;
        if ( pt < fMaxPtUE )
          etaBandPtSumClusterCutMax += pt;
        if ( pt < coneptLeadCluster*fLeadingPtUEFactor )
          etaBandPtSumClusterCutLeadFactor += pt;

        if ( fFillHistograms && fICMethod == kSumBkgSubEtaBandIC )
        {
          if ( fFillHighMultHistograms )
            fhEtaBandClusterPtCent->Fill(ptC, pt, centrality, histoWeight);
          else
            fhEtaBandClusterPt->Fill(ptC, pt, histoWeight);

          if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
          {
            fhEtaBandClusterEtaPhi->Fill(eta, phi, histoWeight);
            fhEtaBandClusterDeltaEtaPhi->Fill(etaC-eta, phiC-phi, histoWeight);
            fhEtaBandClusterPtTriggerEtaPhi->Fill(etaC, phiC, pt, histoWeight);
          }
        }
      } // eta band
    } // out of cone
  } // UE neutral particle loop
  
  //printf("end loop, fill histo\n");
  
  if ( fFillHistograms ) // Do not fill in perp cones case
  {
    //printf("A\n");
    if( fPartInCone == kNeutralAndCharged )
    {
      if ( fFillHighMultHistograms )
      {
        fhConeSumPtClusterCent->Fill(ptC, coneptsumCluster, centrality, histoWeight);

        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
          fhConeSumPtClusterTrigEtaPhiCent[cenBin] ->Fill(etaC, phiC, coneptsumCluster , histoWeight);
      }
      else
      {
        fhConeSumPtCluster ->Fill(ptC, coneptsumCluster, histoWeight);
        fhConeSumPtClusterMatched->Fill(ptC, coneptsumClusterMatched, histoWeight);
        if ( coneptsumCluster > 0 )
        {
          //printf("pt %f, Sum pt %2.2f, matched %2.2f\n", ptC, coneptsumCluster, coneptsumClusterMatched);
          Float_t matchFraction = coneptsumClusterMatched / (coneptsumClusterMatched+coneptsumCluster);
          fhConeSumPtClusterMatchedFraction  ->Fill(ptC, matchFraction, histoWeight);
          fhConeSumPtClusterMatchedFraction3D->Fill(ptC, coneptsumCluster, matchFraction, histoWeight);
        }
        
        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
          fhConeSumPtClusterTrigEtaPhi ->Fill(etaC, phiC, coneptsumCluster , histoWeight);

        if ( fFillLeadingHistograms ) 
          fhConePtLeadCluster->Fill(ptC, coneptLeadCluster, histoWeight);
      }
    }  
    
    //printf("B\n");

    // UE substraction
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      if ( fFillHighMultHistograms )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
          fhConeSumPtEtaBandUEClusterCent->Fill(ptC, etaBandPtSumCluster, centrality, histoWeight);

        if ( fICMethod == kSumBkgSubPhiBandIC )
          fhConeSumPtPhiBandUEClusterCent->Fill(ptC, phiBandPtSumCluster, centrality, histoWeight);

        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
            fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, etaBandPtSumCluster, histoWeight);

          if ( fICMethod == kSumBkgSubPhiBandIC )
            fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, phiBandPtSumCluster, histoWeight);
        }
      }
      else
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUECluster->Fill(ptC, etaBandPtSumCluster, histoWeight);
        }

        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUECluster->Fill(ptC, phiBandPtSumCluster, histoWeight);
        }

        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
            fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaC, phiC, etaBandPtSumCluster, histoWeight); // Check

          if ( fICMethod == kSumBkgSubPhiBandIC )
            fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaC, phiC, phiBandPtSumCluster, histoWeight); // Check
        }
      }
    } // UE sub
  } // fill histo
  
  // Add reference clusters arrays to AOD when filling AODs only
  // Add selected clusters in cone to pCandidate, might be modified later 
  // if UE subtraction and normalization requested
  //
  if ( bFillAOD && refclusters ) pCandidate->AddObjArray(refclusters);  
}

//_________________________________________________________________________________________________________________________________
/// Get the pt sum of the tracks inside the cone and UE regions, the leading track pT and number of clusters.
/// Pass the calculated pT values, but also set them in pCandidate
///
/// \param pCandidate: Kinematics and + of candidate particle for isolation.
/// \param reader: pointer to AliCaloTrackReader. Needed to access event info.
/// \param bFillAOD: Indicate if particles in cone must be added to AOD particle object.
/// \param useRefs: Get the list of tracks or clusters in cone from references
/// \param aodArrayRefName: Name of array where list of tracks/clusters in cone is stored.
/// \param bgTrk: optional external array of tracks
/// \param nPart: number of tracks/clusters above threshold in cone, output.
/// \param nfrac: 1 if fraction pT cluster-track / pT trigger in cone avobe threshold, output.
/// \param coneptsumTrack: total momentum energy in cone track, output.
/// \param coneptLeadTrack: momentum of leading track in cone, output.
/// \param etaBandPtSumTrack: sum of tracks in eta band, same phi as candidate
/// \param etaBandPtLeadTrack: momentum of leading track in eta band
/// \param etaBandPtSumTrackCutMax: sum of tracks in eta band, same phi as candidate, pt < max pt
/// \param etaBandPtSumTrackCutLeadFactor: sum of tracks in eta band, same phi as candidate, pt < lead in cone pt
/// \param phiBandPtSumTrack: sum of tracks in phi band, same eta as candidate (not useful need to restrict it)
/// \param phiBandPtLeadTrack: momentum of leading track in phi band
/// \param phiBandPtSumTrackCutMax: sum of tracks in phi band, same eta as candidate (not useful need to restrict it), pt < max pt
/// \param phiBandPtSumTrackCutLeadFactor: sum of tracks in phi band, same eta as candidate (not useful need to restrict it), pt < lead in cone pt
/// \param perpConePtSumTrack: sum of tracks in perpendicular cones in phi, return result divided by 2
/// \param perpConePtLeadTrack: momentum of leading track in  perpendicular cones
/// \param perpConePtSumTrackCutMax: sum of tracks in perpendicular cones in phi, return result divided by 2, pt < max pt
/// \param perpConePtSumTrackCutLeadFactor: sum of tracks in perpendicular cones in phi, return result divided by 2, pt < lead in cone pt
/// \param perpBandPtSumTrack: sum of tracks in perpendicular  eta band, return result divided by 2
/// \param perpBandPtLeadTrack: momentum of leading track in perpendicular bands
/// \param perpBandPtSumTrackCutMax: sum of tracks in perpendicular  eta band, return result divided by 2, pt <max pt
/// \param perpBandPtSumTrackCutLeadFactor: sum of tracks in perpendicular  eta band, return result divided by 2, pt < lead in cone pt
/// \param perpCone1PtSumTrack: sum of tracks in 1 perpendicular cone in phi
/// \param histoWeight: Histograms weight (event, pt depedent)
/// \param centrality:Centrality percentile
/// \param cenBin: Assigned centrality bin with defined range elsewhere
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateTrackSignalInCone
(
 AliCaloTrackParticleCorrelation * pCandidate, AliCaloTrackReader * reader,
 Bool_t    bFillAOD         , Bool_t    useRefs, 
 TString   aodArrayRefName  , TObjArray * bgTrk,
 Int_t   & nPart            , Int_t   & nfrac,
 Float_t & coneptsumTrack   , Float_t & coneptLeadTrack,
 Float_t & etaBandPtSumTrack, Float_t &  etaBandPtLeadTrack, Float_t &  etaBandPtSumTrackCutMax, Float_t &  etaBandPtSumTrackCutLeadFactor,
 Float_t & phiBandPtSumTrack, Float_t &  phiBandPtLeadTrack, Float_t &  phiBandPtSumTrackCutMax, Float_t &  phiBandPtSumTrackCutLeadFactor,
 Float_t & perpConePtSumTrack,Float_t & perpConePtLeadTrack, Float_t & perpConePtSumTrackCutMax, Float_t & perpConePtSumTrackCutLeadFactor,
 Float_t & perpBandPtSumTrack,Float_t & perpBandPtLeadTrack, Float_t & perpBandPtSumTrackCutMax, Float_t & perpBandPtSumTrackCutLeadFactor,
 Float_t & perpCone1PtSumTrack,
 Double_t  histoWeight      ,
 Float_t   centrality       , Int_t    cenBin
) 
{  
  if ( fPartInCone == kOnlyNeutral ) return ;

  // Get the array with tracks
  //
  TObjArray * plCTS    = 0x0; 
  if ( bgTrk )
  {
    plCTS = bgTrk;
  }
  else if ( !useRefs )
  {
    plCTS = reader->GetCTSTracks();
  }
  else
  {
    plCTS = pCandidate->GetObjArray(aodArrayRefName+"Tracks");
  }
  
  if ( !plCTS ) 
  {
    // In case of no entries, add 0 to histogram to keep number of trigger entries
    if ( fFillHistograms  && fFillHighMultHistograms )
    {
      fhConeSumPtTrack      ->Fill(pCandidate->Pt(), 0., histoWeight);

      if ( fFillLeadingHistograms )
        fhConePtLeadTrack     ->Fill(pCandidate->Pt(), 0., histoWeight);
    } // histograms

    if ( coneptLeadTrack > 0  || coneptsumTrack > 0 )
      AliError(Form("No ref tracks!!! sum %f, lead %f",coneptsumTrack,coneptLeadTrack));

    return ;
  }

  //-----------------------------------------------------------
  // Init parameters
  //
  Float_t ptTrig   = pCandidate->Pt() ;
  Float_t phiTrig  = pCandidate->Phi() ;
  if ( phiTrig < 0 ) phiTrig+=TMath::TwoPi();
  Float_t etaTrig  = pCandidate->Eta() ;

  Float_t ptTrack  = -100. ;
  Float_t etaTrack = -100. ;
  Float_t phiTrack = -100. ;
  Float_t rad      = -100. ;

  TObjArray * reftracks  = 0x0;
  Int_t       ntrackrefs = 0;

  //-----------------------------------------------------------
  // Get the tracks in cone
  //-----------------------------------------------------------

  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ )
  {
    AliVTrack* track = dynamic_cast<AliVTrack*>(plCTS->At(ipr)) ;

    if ( track )
    {
      // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
      // do not count the candidate or the daughters of the candidate
      // in the isolation conte
      if ( pCandidate->GetDetectorTag() == AliFiducialCut::kCTS ) // make sure conversions are tagged as kCTS!!!
      {
        Int_t  trackID   = reader->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
        Bool_t contained = kFALSE;

        for(Int_t i = 0; i < 4; i++)
        {
          if( trackID == pCandidate->GetTrackLabel(i) ) contained = kTRUE;
        }

        if ( contained ) continue ;
      }

      fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
      ptTrack  = fTrackVector.Pt();
      etaTrack = fTrackVector.Eta();
      phiTrack = fTrackVector.Phi() ;
    }
    else
    {// Mixed event stored in AliCaloTrackParticles
      AliCaloTrackParticle * trackmix = dynamic_cast<AliCaloTrackParticle*>(plCTS->At(ipr)) ;
      if(!trackmix)
      {
        AliWarning("Wrong track data type, continue");
        continue;
      }

      ptTrack  = trackmix->Pt();
      etaTrack = trackmix->Eta();
      phiTrack = trackmix->Phi() ;
    }

    // ** Calculate distance between candidate and tracks **

    if ( phiTrack < 0 ) phiTrack+=TMath::TwoPi();

    rad = Radius(etaTrig, phiTrig, etaTrack, phiTrack);

    // ** Exclude tracks too close to the candidate, inactive by default **

    if ( rad < fDistMinToTrigger ) continue ;

    //-------------------------------------------------------------
    // ** For the isolated particle **
    //-------------------------------------------------------------
    //
    //      // Only loop the particle at the same side of candidate
    //      if ( TMath::Abs(phiTrack-phiTrig) > TMath::PiOver2() ) continue ;
    //
    //      // If at the same side has particle larger than candidate,
    //      // then candidate can not be the leading, skip such events
    //      if(pt > ptTrig)
    //      {
    //        n         = -1;
    //        nfrac     = -1;
    //        coneptsumTrack = -1;
    //        isolated  = kFALSE;
    //
    //        pCandidate->SetLeadingParticle(kFALSE);
    //
    //        if(bFillAOD && reftracks)
    //        {
    //          reftracks->Clear();
    //          delete reftracks;
    //        }
    //
    //        return ;
    //      }
    //
    //AliDebug(2,Form("\t Track %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f",
    //                ipr,ptTrack,etaTrack,phiTrack,rad));

    //
    // Select tracks inside the isolation radius
    //
    if ( rad > fConeSize ) continue ;

    AliDebug(2,"Inside candidate cone");

    if ( bFillAOD )
    {
      ntrackrefs++;

      if ( ntrackrefs == 1 )
      {
        reftracks = new TObjArray(0);
        //reftracks->SetName(Form("Tracks%s",aodArrayRefName.Data()));
        TString tempo(aodArrayRefName)  ;
        tempo += "Tracks" ;
        reftracks->SetName(tempo);
        reftracks->SetOwner(kFALSE);
      }

      reftracks->Add(track);
    }

    coneptsumTrack+=ptTrack;

    if ( coneptLeadTrack < ptTrack ) coneptLeadTrack = ptTrack;

    if ( fFillHistograms )
    {
      if ( !fFillHighMultHistograms )
      {
        fhPtInCone->Fill(ptTrig , ptTrack, histoWeight);

        if( fPartInCone == kNeutralAndCharged )
          fhPtTrackInCone->Fill(ptTrig , ptTrack, histoWeight);
      }
      else
      {
        fhPtInConeCent->Fill(ptTrig , ptTrack, centrality, histoWeight);

        if( fPartInCone == kNeutralAndCharged )
          fhPtTrackInConeCent->Fill(ptTrig , ptTrack, centrality, histoWeight);
      }

      if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
      {
        fhEtaPhiInConeTrack->Fill(etaTrack, phiTrack, histoWeight);
        fhDeltaEtaPhiInConeTrack->Fill(etaTrig-etaTrack, phiTrig-phiTrack, histoWeight);
        fhTriggerEtaPhiInConeTrackPt->Fill(etaTrig, phiTrig, ptTrack, histoWeight);
      }
    } // histograms

    // *Before*, count particles in cone
    if ( ptTrack > fPtThreshold && ptTrack < fPtThresholdMax )  nPart++;

    //if fPtFraction*ptTrig<fPtThreshold then consider the fPtThreshold directly
    if ( fFracIsThresh )
    {
      if ( fPtFraction*ptTrig < fPtThreshold )
      {
        if ( ptTrack > fPtThreshold )    nfrac++ ;
      }
      else
      {
        if ( ptTrack > fPtFraction*ptTrig ) nfrac++;
      }
    }
    else
    {
      if ( ptTrack > fPtFraction*ptTrig ) nfrac++;
    }
  } // track loop

  //-----------------------------------------------------------
  // Select the UE tracks
  //-----------------------------------------------------------

  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ )
  {
    AliVTrack* track = dynamic_cast<AliVTrack*>(plCTS->At(ipr)) ;

    if ( track )
    {
      // In case of isolation of single tracks or conversion photon (2 tracks) or pi0 (4 tracks),
      // do not count the candidate or the daughters of the candidate
      // in the isolation conte
      if ( pCandidate->GetDetectorTag() == AliFiducialCut::kCTS ) // make sure conversions are tagged as kCTS!!!
      {
        Int_t  trackID   = reader->GetTrackID(track) ; // needed instead of track->GetID() since AOD needs some manipulations
        Bool_t contained = kFALSE;

        for(Int_t i = 0; i < 4; i++)
        {
          if( trackID == pCandidate->GetTrackLabel(i) ) contained = kTRUE;
        }

        if ( contained ) continue ;
      }

      fTrackVector.SetXYZ(track->Px(),track->Py(),track->Pz());
      ptTrack  = fTrackVector.Pt();
      etaTrack = fTrackVector.Eta();
      phiTrack = fTrackVector.Phi() ;
    }
    else
    {// Mixed event stored in AliCaloTrackParticles
      AliCaloTrackParticle * trackmix = dynamic_cast<AliCaloTrackParticle*>(plCTS->At(ipr)) ;
      if(!trackmix)
      {
        AliWarning("Wrong track data type, continue");
        continue;
      }

      ptTrack  = trackmix->Pt();
      etaTrack = trackmix->Eta();
      phiTrack = trackmix->Phi() ;
    }

    // ** Calculate distance between candidate and tracks **

    if ( phiTrack < 0 ) phiTrack+=TMath::TwoPi();

    rad = Radius(etaTrig, phiTrig, etaTrack, phiTrack);

    if ( fFillHistograms && fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
      fhEtaPhiTrack->Fill(etaTrack, phiTrack, histoWeight);
    
    // Angles between trigger and track
    Double_t dEta = etaTrig - etaTrack;
    Double_t dPhi = phiTrig - phiTrack;

    // Shift phi angle when trigger is close to 0 or 360
    if ( dPhi >=  TMath::Pi() ) dPhi-=TMath::TwoPi();
    if ( dPhi <= -TMath::Pi() ) dPhi+=TMath::TwoPi();

    if ( fICMethod > kSumBkgSubIC && rad > fConeSize+fConeSizeBandGap )
    {
      // Phi band
      //
      Bool_t takeIt = kTRUE;
      
      // Look only half TPC with respect candidate, avoid opposite side jet 
      if ( TMath::Abs(dPhi) > TMath::PiOver2() )  takeIt = kFALSE;
      
      // Within eta cone size
      if ( TMath::Abs(dEta) < (fConeSize+fConeSizeBandGap) &&  takeIt )
      {
        phiBandPtSumTrack += ptTrack;
        if ( ptTrack < fMaxPtUE )
          phiBandPtSumTrackCutMax += ptTrack;
        if ( ptTrack < coneptLeadTrack*fLeadingPtUEFactor )
          phiBandPtSumTrackCutLeadFactor += ptTrack;
        if ( ptTrack > phiBandPtLeadTrack )
          phiBandPtLeadTrack = ptTrack ;

        if ( fFillHistograms && fICMethod == kSumBkgSubPhiBandIC )
        {
          if ( fFillHighMultHistograms )
            fhPhiBandTrackPtCent->Fill(ptTrig, ptTrack, centrality, histoWeight);
          else
            fhPhiBandTrackPt->Fill(ptTrig , ptTrack, histoWeight);

          if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
          {
            fhPhiBandTrackEtaPhi->Fill(etaTrack, phiTrack, histoWeight);
            fhPhiBandTrackDeltaEtaPhi->Fill(etaTrig-etaTrack, phiTrig-phiTrack, histoWeight);
            fhPhiBandTrackPtTriggerEtaPhi->Fill(etaTrig, phiTrig, ptTrack, histoWeight);
          }
        }
      } // phi band

      // Eta band
      //
      //takeIt = kTRUE;

      // Within phi cone size
      if ( TMath::Abs(dPhi) < (fConeSize+fConeSizeBandGap) )// && takeIt )
      {
        etaBandPtSumTrack += ptTrack;
        if ( ptTrack < fMaxPtUE )
          etaBandPtSumTrackCutMax += ptTrack;
        if ( ptTrack < coneptLeadTrack*fLeadingPtUEFactor )
          etaBandPtSumTrackCutLeadFactor += ptTrack;
        if ( ptTrack > etaBandPtLeadTrack )
          etaBandPtLeadTrack = ptTrack ;

        if ( fFillHistograms && fICMethod == kSumBkgSubEtaBandIC )
        {
          if ( fFillHighMultHistograms )
            fhEtaBandTrackPtCent->Fill(ptTrig, ptTrack, centrality, histoWeight);
          else
            fhEtaBandTrackPt->Fill(ptTrig , ptTrack, histoWeight);

          if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
          {
            fhEtaBandTrackEtaPhi->Fill(etaTrack, phiTrack, histoWeight);
            fhEtaBandTrackDeltaEtaPhi->Fill(etaTrig-etaTrack, phiTrig-phiTrack, histoWeight);
            fhEtaBandTrackPtTriggerEtaPhi->Fill(etaTrig, phiTrig, ptTrack, histoWeight);
          }
        }
      } // eta band
      
    } // out of cone

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fill the histograms at +-45 degrees in phi  
    // from trigger particle, perpedicular to trigger axis in phi
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( fICMethod == kSumBkgSubIC )
    {
      Double_t dPhiPlu = dPhi + TMath::PiOver2();
      Double_t dPhiMin = dPhi - TMath::PiOver2();
      
      Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
      Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;
      
      Bool_t fillPerp = kFALSE;
      if ( TMath::Sqrt(argPlu) < fConeSize ) { fillPerp = kTRUE ; perpCone1PtSumTrack+=ptTrack; }
      if ( TMath::Sqrt(argMin) < fConeSize ) fillPerp = kTRUE ;
      
      if ( fillPerp )
      {
        perpConePtSumTrack+=ptTrack;
        if ( ptTrack < fMaxPtUE )
          perpConePtSumTrackCutMax += ptTrack;
        if ( ptTrack < coneptLeadTrack*fLeadingPtUEFactor )
          perpConePtSumTrackCutLeadFactor += ptTrack;
        if ( ptTrack > perpConePtLeadTrack )
          perpConePtLeadTrack = ptTrack ;

        if ( fFillHistograms )
        {
          if ( !fFillHighMultHistograms )
            fhPtInPerpCone->Fill(ptTrig, ptTrack, histoWeight);
          else
            fhPtInPerpConeCent->Fill(ptTrig, ptTrack, centrality, histoWeight);
          
          if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
          {
            fhEtaPhiInPerpCone->Fill(etaTrack,phiTrack, histoWeight);
            fhDeltaEtaPhiInPerpCone->Fill(etaTrig-etaTrack, phiTrig-phiTrack, histoWeight);
            fhTriggerEtaPhiInPerpConeTrackPt->Fill(etaTrig, phiTrig, ptTrack, histoWeight);
          }
        }
      }
    }

    if ( fICMethod == kSumBkgSubPerpBandIC )
    {
      Double_t dPhiPlu = dPhi + TMath::PiOver2();
      Double_t dPhiMin = dPhi - TMath::PiOver2();

      if ( TMath::Abs(dPhiPlu) < fConeSize ||
           TMath::Abs(dPhiMin) < fConeSize   )
      {
        perpBandPtSumTrack += ptTrack;
        if ( ptTrack < fMaxPtUE )
          perpBandPtSumTrackCutMax += ptTrack;
        if ( ptTrack < coneptLeadTrack*fLeadingPtUEFactor )
          perpBandPtSumTrackCutLeadFactor += ptTrack;
        if ( ptTrack > perpBandPtLeadTrack )
          perpBandPtLeadTrack = ptTrack ;

        if ( fFillHistograms )
        {
          if ( fFillHighMultHistograms )
            fhPerpBandTrackPtCent->Fill(ptTrig, ptTrack, centrality, histoWeight);
          else
            fhPerpBandTrackPt->Fill(ptTrig , ptTrack, histoWeight);

          if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
          {
            fhPerpBandTrackEtaPhi->Fill(etaTrack, phiTrack, histoWeight);
            fhPerpBandTrackDeltaEtaPhi->Fill(etaTrig-etaTrack, phiTrig-phiTrack, histoWeight);
            fhPerpBandTrackPtTriggerEtaPhi->Fill(etaTrig, phiTrig, ptTrack, histoWeight);
          }
        }
      } // perp eta band
    }
  }// charged particle loop

  // 2 perpendicular cones added, divide by 2 total amount of energy.
  perpConePtSumTrack /= 2;
  perpConePtSumTrackCutMax /= 2;
  perpConePtSumTrackCutLeadFactor /= 2;

  if ( fFillHistograms )
  {
    if ( fPartInCone == kNeutralAndCharged )
    {
      if ( fFillHighMultHistograms )
      {
        fhConeSumPtTrackCent ->Fill(ptTrig, coneptsumTrack, centrality, histoWeight);
        
        if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
          fhConeSumPtTrackTrigEtaPhiCent[cenBin]->Fill(etaTrig, phiTrig, coneptsumTrack, histoWeight);
      }
      else
      {
        fhConeSumPtTrack ->Fill(ptTrig, coneptsumTrack, histoWeight);
        
        if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
          fhConeSumPtTrackTrigEtaPhi ->Fill(etaTrig, phiTrig, coneptsumTrack, histoWeight);

        if ( fFillLeadingHistograms ) 
          fhConePtLeadTrack->Fill(ptTrig, coneptLeadTrack, histoWeight);
      }
    }
    
    // UE subtraction
    if ( fICMethod > kSumBkgSubIC )
    {
      if ( fFillHighMultHistograms )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
          fhConeSumPtEtaBandUETrackCent->Fill(ptTrig, etaBandPtSumTrack , centrality, histoWeight);

        if ( fICMethod == kSumBkgSubPhiBandIC )
          fhConeSumPtPhiBandUETrackCent->Fill(ptTrig, phiBandPtSumTrack , centrality, histoWeight);

        if ( fICMethod == kSumBkgSubPerpBandIC )
          fhConeSumPtPerpBandUETrackCent->Fill(ptTrig,perpBandPtSumTrack , centrality, histoWeight);

        if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
            fhConeSumPtEtaBandUETrackTrigEtaPhiCent[cenBin]->Fill(etaTrig, phiTrig, etaBandPtSumTrack, histoWeight);

          if ( fICMethod == kSumBkgSubPhiBandIC )
            fhConeSumPtPhiBandUETrackTrigEtaPhiCent[cenBin]->Fill(etaTrig, phiTrig, phiBandPtSumTrack, histoWeight);

          if ( fICMethod == kSumBkgSubPerpBandIC )
            fhConeSumPtPerpBandUETrackTrigEtaPhiCent[cenBin]->Fill(etaTrig, phiTrig,perpBandPtSumTrack, histoWeight);
        }
      }
      else
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUETrack->Fill(ptTrig, etaBandPtSumTrack , histoWeight);
        }

        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUETrack->Fill(ptTrig, phiBandPtSumTrack , histoWeight);
        }

        if ( fICMethod == kSumBkgSubPerpBandIC )
        {
          fhConeSumPtPerpBandUETrack->Fill(ptTrig,perpBandPtSumTrack , histoWeight);
        }

        if ( fFillEtaPhiHistograms && ptTrig > fEtaPhiHistogramsMinPt )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
            fhConeSumPtEtaBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSumTrack, histoWeight);

          if ( fICMethod == kSumBkgSubPhiBandIC )
            fhConeSumPtPhiBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSumTrack, histoWeight);

          if ( fICMethod == kSumBkgSubPerpBandIC )
            fhConeSumPtPerpBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig,perpBandPtSumTrack, histoWeight);
        }
      }
    } // UE sub
  } // fill histograms

  // Add reference track arrays to AOD when filling AODs only
  // Add selected tracks in cone to pCandidate, might be modified later 
  // if UE subtraction and normalization requested
  //
  if ( bFillAOD && reftracks ) pCandidate->AddObjArray(reftracks);  
}

//_________________________________________________________________________________________________________________________________
/// Get normalization of cluster background band.
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandClusterNormalization
(
 Int_t     calorimeter,
 Float_t   etaC,                  Float_t   phiC,
 Float_t   excessEta,             Float_t   excessPhi,         
 Float_t   excessAreaEta,         Float_t   excessAreaPhi,         
 Float_t   etaUEptsumCluster,     Float_t   phiUEptsumCluster,
 Float_t & etaUEptsumClusterNorm, Float_t & phiUEptsumClusterNorm) const
{
  // Cone area
  // Area = pi R^2, isolation cone area
  Float_t coneArea = fConeSize*fConeSize*TMath::Pi();
  if ( fDistMinToTrigger > 0 ) coneArea -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  // Area = pi (R-gap)^2, isolation cone area
  Float_t coneAreaGap = (fConeSize + fConeSizeBandGap)*(fConeSize + fConeSizeBandGap)*TMath::Pi();
  if ( fDistMinToTrigger > 0 ) coneAreaGap -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  Float_t fEMCPhiSize = fEMCPhiMax-fEMCPhiMin;
  // DCal?
  if ( (phiC > 260*TMath::DegToRad() && phiC < 327*TMath::DegToRad()) &&  calorimeter != AliFiducialCut::kPHOS )
  {
    fEMCPhiSize = fDMCPhiMax-fDMCPhiMin;
  }
  
  Float_t excessEtaGap = excessEta,  excessPhiGap = excessPhi;
  Float_t excessAreaEtaGap = excessAreaEta,  excessAreaPhiGap = excessAreaPhi;
  if ( fConeSizeBandGap > 0 )
  {
    if ( TMath::Abs(etaC)+fConeSize+fConeSizeBandGap > fEMCEtaSize/2. )
      excessEtaGap = TMath::Abs(etaC) + fConeSize + fConeSizeBandGap - fEMCEtaSize/2.;

    if     ( TMath::Abs(phiC)+fConeSize+fConeSizeBandGap > fEMCPhiMax )
      excessPhiGap = (TMath::Abs(phiC) + fConeSize + fConeSizeBandGap) - fEMCPhiMax;
    else if( TMath::Abs(phiC)-fConeSize-fConeSizeBandGap < fEMCPhiMin )
      excessPhiGap = fEMCPhiMin - (TMath::Abs(phiC) - fConeSize - fConeSizeBandGap) ;

    excessAreaEtaGap = CalculateExcessAreaFraction(excessEtaGap,fConeSizeBandGap);
    excessAreaPhiGap = CalculateExcessAreaFraction(excessPhiGap,fConeSizeBandGap);
  }

  // UE band can also be out of acceptance, need to estimate corrected area
  //
  if ( excessEta != 0 ) coneArea /= excessAreaEta;
  if ( excessPhi != 0 ) coneArea /= excessAreaPhi;
  
  if ( excessEtaGap != 0 ) coneAreaGap /= excessAreaEtaGap;
  if ( excessPhiGap != 0 ) coneAreaGap /= excessAreaPhiGap;
  
  // Area of band, rectangle minus isolation region
  //
  Float_t etaBandArea = ( 2 * (fConeSize + fConeSizeBandGap) - excessPhiGap ) * fEMCEtaSize - coneAreaGap;
  Float_t phiBandArea = ( 2 * (fConeSize + fConeSizeBandGap) - excessEtaGap ) * fEMCPhiSize - coneAreaGap;
  
  // Calculate normalization factor and rescale UE band sum pT to cone area
  // 
  // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
  if ( phiBandArea > 0 ) 
    phiUEptsumClusterNorm = phiUEptsumCluster * ( coneArea / phiBandArea ); 
  
  // pi * R^2 / (2 R * 2*0.7)  -  trigger cone
  if ( etaBandArea > 0 ) 
    etaUEptsumClusterNorm = etaUEptsumCluster * ( coneArea / etaBandArea ); 
  
//  printf("clust band area: cone Area %2.2f, Gap Area eta %2.2f phi %2.2f, EMC eta size %2.2f, EMC phi size %2.2f\n"
//         "etaBandArea %2.2f phiBandArea %2.2f, norm factor phi %2.2f, norm factor eta %2.2f\n",
//         coneArea, coneAreaEta, coneAreaPhi, fEMCEtaSize,fEMCPhiSize,
//         etaBandArea,phiBandArea,coneArea/phiBandArea,coneArea/etaBandArea);
  
  if ( etaBandArea < 0 || phiBandArea < 0 )
    printf("Negative clust band area: cone Area %2.2f, Gap Area Gap %2.2f, EMC eta size %2.2f, EMC phi size %2.2f\n",
           coneArea, coneAreaGap, fEMCEtaSize, fEMCPhiSize);
}

//________________________________________________________________________________________________________________________________
/// Get normalization of track background band.
//________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandTrackNormalization  
(
 Float_t   etaC               ,  
 Float_t   excessEta          ,          
 Float_t   excessAreaEta      ,         
 Float_t   etaUEptsumTrack    ,  Float_t   phiUEptsumTrack,     Float_t   perpUEptsumTrack,
 Float_t & etaUEptsumTrackNorm,  Float_t & phiUEptsumTrackNorm, Float_t & perpUEptsumTrackNorm ) const
{
  // Cone area
  // Area = pi R^2, isolation cone area
  Float_t coneArea = fConeSize*fConeSize*TMath::Pi();
  if ( fDistMinToTrigger > 0 ) coneArea -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  // Area = pi ((R+gap)^2), isolation cone area
  Float_t coneAreaGap = (fConeSize + fConeSizeBandGap)*(fConeSize + fConeSizeBandGap)*TMath::Pi();
  if ( fDistMinToTrigger > 0 ) coneAreaGap -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  Float_t excessEtaGap = excessEta;
  Float_t excessAreaEtaGap = excessAreaEta;
  if ( fConeSizeBandGap > 0 )
  {
    if ( TMath::Abs(etaC)+fConeSize > fTPCEtaSize/2. )
      excessEtaGap = TMath::Abs(etaC) + fConeSize + fConeSizeBandGap - fTPCEtaSize/2.;

    excessAreaEtaGap = CalculateExcessAreaFraction(excessEtaGap,fConeSizeBandGap);
  }
  
  // UE band can also be out of acceptance, need to estimate corrected area
  //
  if ( excessEta != 0 ) coneArea /= excessAreaEta;

  if ( excessEtaGap != 0 ) coneAreaGap /= excessAreaEtaGap;
  
  // Area of band, rectangle minus isolation region
  //
  Float_t perpBandArea=   4* fConeSize                                    * fTPCEtaSize ;
  Float_t etaBandArea =   2*(fConeSize+fConeSizeBandGap)                  * fTPCEtaSize - coneAreaGap;
  Float_t phiBandArea = ( 2*(fConeSize+fConeSizeBandGap) - excessEtaGap ) * fTPCPhiSize - coneAreaGap;
  
  // Calculate normalization factor and rescale UE band sum pT to cone area
  //
  // pi * R^2 / (2 R * 2 pi) -  trigger cone
  if ( phiBandArea > 0 ) 
    phiUEptsumTrackNorm = phiUEptsumTrack * ( coneArea / phiBandArea ); 
  
  // pi * R^2 / (2 R * 1.6)  -  trigger cone
  if ( etaBandArea > 0 ) 
    etaUEptsumTrackNorm = etaUEptsumTrack * ( coneArea / etaBandArea ); 

  if ( perpBandArea > 0 )
    perpUEptsumTrackNorm = perpUEptsumTrack * ( coneArea / perpBandArea );
  
//  printf("track band area: cone Area %2.2f, Gap Area eta %2.2f phi %2.2f, EMC eta size %2.2f, EMC phi size %2.2f\n"
//          "etaBandArea %2.2f phiBandArea %2.2f, norm factor phi %2.2f, norm factor eta %2.2f\n",
//          coneArea, coneAreaEta, coneAreaPhi, fTPCEtaSize,fTPCPhiSize,
//          etaBandArea,phiBandArea,coneArea/phiBandArea,coneArea/etaBandArea);
  
  if ( etaBandArea < 0 || phiBandArea < 0 || perpBandArea < 0 )
    printf("Negative track band area: cone Area %2.2f, Gap Area eta %2.2f, TPC eta size %2.2f, TPC phi size %2.2f\n",
           coneArea, coneAreaGap, fTPCEtaSize, fTPCPhiSize);
}

//________________________________________________________________________
/// If isolation cone is outside a detector, calculate the area in excess.
/// \param excess: cone size minus acceptance of detector.
/// \param gap: fConeSizeBandGap, passed when correcting bkg excess
/// \return Area of a circunference segment 1/2 R^2 (angle-sin(angle)), angle = 2*ACos((R-excess)/R).
//________________________________________________________________________
Float_t AliIsolationCut::CalculateExcessAreaFraction(Float_t excess, Float_t gap) const
{
  Float_t coneSize = fConeSize + gap;
  Float_t angle   = 2*TMath::ACos( (coneSize-excess) / coneSize );

  Float_t coneA   = coneSize*coneSize*TMath::Pi(); // A = pi R^2, isolation cone area

  Float_t excessA = coneSize*coneSize / 2 * (angle-TMath::Sin(angle));

  if ( coneA > excessA )
  {
    return coneA / (coneA-excessA);
  }
  else
  {
    AliWarning(Form("Please Check : R = %2.2f + gap = %2.2f, excess %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f",
                    fConeSize,gap, excess,coneA, excessA, angle*TMath::RadToDeg()));
    return  1;
  }
}

//________________________________________________________________________
/// If isolation cone are is outside a detector, calculate the area in excess.
/// \param etaC: Candidate pseudorapidity.
/// \param phiC: Candidate azimuthal angle
/// \param det: Detector used
/// \param excessAreaTrkEta: cone excess out of TPC pseudo rapidity
/// \param excessAreaClsEta: cone excess out of EMC pseudo rapidity
/// \param excessAreaClsPhi: cone excess out of EMC azimuthal angle
//________________________________________________________________________
void AliIsolationCut::CalculateExcessAreaFractionForChargedAndNeutral
( Float_t etaC, Float_t phiC, Int_t det,
  Float_t & excessTrkEta, Float_t & excessAreaTrkEta, 
  Float_t & excessClsEta, Float_t & excessAreaClsEta, 
  Float_t & excessClsPhi, Float_t & excessAreaClsPhi  ) const
{
  excessTrkEta = 0;
  excessClsEta = 0;
  excessClsPhi = 0;
  Float_t excessClsEta0 = 0;

  if ( fPartInCone != kOnlyNeutral )
  {
    if ( TMath::Abs(etaC)+fConeSize > fTPCEtaSize/2. )
      excessTrkEta = TMath::Abs(etaC) + fConeSize - fTPCEtaSize/2.;
    
    if ( excessTrkEta < 0 )
      AliInfo(Form("Fix negative excess: trk eta %f",excessTrkEta));

    excessAreaTrkEta = CalculateExcessAreaFraction(excessTrkEta);
  }

  if ( fPartInCone != kOnlyCharged &&
      (det == AliCaloTrackReader::kEMCAL || det == AliCaloTrackReader::kDCAL || det == AliCaloTrackReader::kPHOS )  )
  {
    if ( TMath::Abs(etaC)+fConeSize > fEMCEtaSize/2. )
      excessClsEta = TMath::Abs(etaC) + fConeSize - fEMCEtaSize/2.;

    if ( phiC < 260*TMath::DegToRad() ) // EMCal
    {
      if     ( TMath::Abs(phiC)+fConeSize > fEMCPhiMax )
        excessClsPhi = (TMath::Abs(phiC) + fConeSize) - fEMCPhiMax;
      else if( TMath::Abs(phiC)-fConeSize < fEMCPhiMin )
        excessClsPhi = fEMCPhiMin - (TMath::Abs(phiC) - fConeSize) ;
    }
    else if ( TMath::Abs(etaC) >= fDMCEtaGap )// DCal
    {
      if ( TMath::Abs(etaC)-fDMCEtaGap < fConeSize )
        excessClsEta0 = fConeSize - (TMath::Abs(etaC)-fDMCEtaGap);
      //printf("DCal Eta C %f - R %f, excess %f\n",etaC, fConeSize, excessClsEta0);
      //printf("PhiC %f, MinPhi %f, MaxPhi %f\n",phiC,fDMCPhiMin,fDMCPhiMax);

      if     ( TMath::Abs(phiC)+fConeSize > fDMCPhiMax )
        excessClsPhi = (TMath::Abs(phiC) + fConeSize) - fDMCPhiMax;
      else if( TMath::Abs(phiC)-fConeSize < fDMCPhiMin )
        excessClsPhi = fDMCPhiMin - (TMath::Abs(phiC) - fConeSize) ;
    }

    if ( excessClsEta < 0 || excessClsEta0 < 0 || excessClsPhi < 0 )
      AliInfo(Form("Fix negative excess: cls eta %f, eta0 %f, cls phi %f",
                   excessClsEta, excessClsEta0, excessClsPhi));

    Float_t excessAreaClsEta0 = CalculateExcessAreaFraction(excessClsEta0);
    //printf("excessAreaClsEta0 %f\n",excessAreaClsEta0);
    excessAreaClsEta  = CalculateExcessAreaFraction(excessClsEta);
    //printf("excessAreaClsEta %f\n",excessAreaClsEta);
    excessAreaClsEta  *= excessAreaClsEta0;
    excessAreaClsPhi  = CalculateExcessAreaFraction(excessClsPhi);
    //printf("excessAreaClsPhi %f\n",excessAreaClsPhi);
  }
}

//_________________________________________________________________________________
/// Set TPC and EMCal angle limits. Do it once.
/// Get the hardcoded value set in the fiducial cut class.
/// \param reader: Access to reader class
//_________________________________________________________________________________
void AliIsolationCut::GetDetectorAngleLimits ( AliCaloTrackReader * reader, Int_t calorimeter )
{
  if ( fEMCEtaSize > 0 ) return ; // Already set
  
  fEMCEtaSize = reader->GetFiducialCut()->GetEMCALFidCutMaxEtaArray()->At(0) -
                reader->GetFiducialCut()->GetEMCALFidCutMinEtaArray()->At(0) ;
  fEMCPhiMin  = reader->GetFiducialCut()->GetEMCALFidCutMinPhiArray()->At(0) ;
  fEMCPhiMax  = reader->GetFiducialCut()->GetEMCALFidCutMaxPhiArray()->At(0) ;
  if ( fEMCPhiMax > 187 )
    fEMCPhiMax = 185.8; // Hack!
//  fDMCEtaGap  = TMath::Abs(reader->GetFiducialCut()->GetDCALFidCutMinEtaArray()->At(0));
//  fDMCPhiMin  = reader->GetFiducialCut()->GetDCALFidCutMinPhiArray()->At(0) ;
//  fDMCPhiMax  = reader->GetFiducialCut()->GetDCALFidCutMaxPhiArray()->At(0) ;

  if ( calorimeter == AliFiducialCut::kPHOS)
  {
    fEMCEtaSize = reader->GetFiducialCut()->GetPHOSFidCutMaxEtaArray()->At(0) -
                  reader->GetFiducialCut()->GetPHOSFidCutMinEtaArray()->At(0) ;
    fEMCPhiMin  = reader->GetFiducialCut()->GetPHOSFidCutMinPhiArray()->At(0) ;
    fEMCPhiMax  = reader->GetFiducialCut()->GetPHOSFidCutMaxPhiArray()->At(0) ;

    fDMCEtaGap  = 0;
    fDMCPhiMin  = 0;
    fDMCPhiMax  = 0;
  }
  
  // Info stored in degrees, put them on radians
  fEMCPhiMin *= TMath::DegToRad();
  fEMCPhiMax *= TMath::DegToRad();
  fDMCPhiMin *= TMath::DegToRad();
  fDMCPhiMax *= TMath::DegToRad();

  // Get the cut used for the TPC tracks in the reader, +-0.8, +-0.9 ...
  // Only valid in simple fidutial cut case and if the cut is applied, careful!
  fTPCEtaSize = reader->GetFiducialCut()->GetCTSFidCutMaxEtaArray()->At(0) -
                reader->GetFiducialCut()->GetCTSFidCutMinEtaArray()->At(0) ;
  fTPCPhiSize = TMath::Pi(); // Half TPC tracks with respect trigger candidate inspected

  AliDebug(1,Form("TPC: Deta %2.2f; Dphi %2.2f; Calo: Deta %2.2f, phi min %2.2f, max %2.2f, DCal eta gap %2.2f,  phi min %2.2f, max %2.2f",
                  fTPCEtaSize,fTPCPhiSize,fEMCEtaSize,
                  fEMCPhiMin*TMath::RadToDeg(),fEMCPhiMax*TMath::RadToDeg(),
                  fDMCEtaGap,
                  fDMCPhiMin*TMath::RadToDeg(),fDMCPhiMax*TMath::RadToDeg()));

//  printf("TPC: Deta %2.2f; Dphi %2.2f; Calo: Deta %2.2f, phi min %2.2f, max %2.2f, DCal eta gap %2.2f,  phi min %2.2f, max %2.2f\n",
//         fTPCEtaSize,fTPCPhiSize,fEMCEtaSize,
//         fEMCPhiMin*TMath::RadToDeg(),fEMCPhiMax*TMath::RadToDeg(),
//         fDMCEtaGap,
//         fDMCPhiMin*TMath::RadToDeg(),fDMCPhiMax*TMath::RadToDeg());
}

//_________________________________________________________________________________
/// Get good cell density (number of active cells over all cells in cone).
//_________________________________________________________________________________
Float_t AliIsolationCut::GetCellDensity(AliCaloTrackParticleCorrelation * pCandidate,
                                        AliCaloTrackReader * reader) const
{
  Double_t coneCells    = 0.; //number of cells in cone with radius fConeSize
  Double_t coneCellsBad = 0.; //number of bad cells in cone with radius fConeSize
  Double_t cellDensity  = 1.;

  Float_t phiC  = pCandidate->Phi() ;
  if(phiC<0) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;

  if(pCandidate->GetDetectorTag() == AliCaloTrackReader::kEMCAL)
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    AliCalorimeterUtils *cu = reader->GetCaloUtils();

    Int_t absId = -999;
    if (eGeom->GetAbsCellIdFromEtaPhi(etaC,phiC,absId))
    {
      // Get absolute (col,row) of candidate
      Int_t iEta=-1, iPhi=-1, iRCU = -1;
      Int_t nSupMod = cu->GetModuleNumberCellIndexes(absId, pCandidate->GetDetectorTag(), iEta, iPhi, iRCU);

      Int_t colC = iEta;
      if (nSupMod % 2) colC =  AliEMCALGeoParams::fgkEMCALCols + iEta ;
      Int_t rowC = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);

      Int_t sqrSize = int(fConeSize/0.0143) ; // Size of cell in radians
      Int_t status = 0;
      // Loop on cells in a square of side fConeSize to check cells in cone
      for(Int_t icol = colC-sqrSize; icol < colC+sqrSize;icol++)
      {
        for(Int_t irow = rowC-sqrSize; irow < rowC+sqrSize; irow++)
        {
          if (Radius(colC, rowC, icol, irow) < sqrSize)
          {
            coneCells += 1.;

            Int_t cellSM  = -999;
            Int_t cellEta = -999;
            Int_t cellPhi = -999;
            if(icol > AliEMCALGeoParams::fgkEMCALCols-1)
            {
              cellSM = 0+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
              cellEta = icol-AliEMCALGeoParams::fgkEMCALCols;
              cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
            }
            if(icol < AliEMCALGeoParams::fgkEMCALCols)
            {
              cellSM = 1+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
              cellEta = icol;
              cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
            }

            //Count as bad "cells" out of EMCAL acceptance
            if(icol < 0 || icol > AliEMCALGeoParams::fgkEMCALCols*2 ||
               irow < 0 || irow > AliEMCALGeoParams::fgkEMCALRows*16./3) //5*nRows+1/3*nRows
            {
              coneCellsBad += 1.;
            }
            // Count as bad "cells" marked as bad in the DataBase
            else if (cu->GetEMCALChannelStatus(cellSM,cellEta,cellPhi,status)==1)
            {
              coneCellsBad += 1. ;
            }
          }
        }
      }//end of cells loop
    }
    else AliWarning("Cluster with bad (eta,phi) in EMCal for energy density calculation");

    if (coneCells > 0.)
    {
      cellDensity = (coneCells-coneCellsBad)/coneCells;
      //printf("Energy density = %f\n", cellDensity);
    }
  }

  return cellDensity;
}

//___________________________________________________________________________________
/// Get good cell density (number of active cells over all cells in cone).
//___________________________________________________________________________________
void AliIsolationCut::GetCoeffNormBadCell(AliCaloTrackParticleCorrelation * pCandidate,
                                          AliCaloTrackReader * reader,
                                          Float_t &  coneBadCellsCoeff,
                                          Float_t &  etaBandBadCellsCoeff,
                                          Float_t & phiBandBadCellsCoeff)
{
  Double_t coneCells    = 0.; //number of cells in cone with radius fConeSize
  Double_t phiBandCells = 0.; //number of cells in band phi
  Double_t etaBandCells = 0.; //number of cells in band eta

  Float_t phiC  = pCandidate->Phi() ;
  if ( phiC < 0 ) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;

  if(pCandidate->GetDetectorTag() == AliCaloTrackReader::kEMCAL)
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    AliCalorimeterUtils *cu = reader->GetCaloUtils();

    Int_t absId = -999;
    if (eGeom->GetAbsCellIdFromEtaPhi(etaC,phiC,absId))
    {
      //Get absolute (col,row) of candidate
      Int_t iEta=-1, iPhi=-1, iRCU = -1;
      Int_t nSupMod = cu->GetModuleNumberCellIndexes(absId, pCandidate->GetDetectorTag(),
                                                     iEta, iPhi, iRCU);

      Int_t colC = iEta;
      if (nSupMod % 2) colC =  AliEMCALGeoParams::fgkEMCALCols + iEta ;
      Int_t rowC = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);

      Int_t sqrSize = int(fConeSize/0.0143) ; // Size of cell in radians
      Int_t status  = 0;
      for(Int_t icol = 0; icol < 2*AliEMCALGeoParams::fgkEMCALCols-1;icol++)
      {
        for(Int_t irow = 0; irow < 5*AliEMCALGeoParams::fgkEMCALRows -1; irow++)
        {
          //loop on cells in a square of side fConeSize to check cells in cone
          if     ( Radius(colC, rowC, icol, irow) < sqrSize ) { coneCells    += 1.; }
          else if( icol>colC-sqrSize  &&  icol<colC+sqrSize ) { phiBandCells += 1 ; }
          else if( irow>rowC-sqrSize  &&  irow<rowC+sqrSize ) { etaBandCells += 1 ; }

          Int_t cellSM  = -999;
          Int_t cellEta = -999;
          Int_t cellPhi = -999;
          if(icol > AliEMCALGeoParams::fgkEMCALCols-1)
          {
            cellSM = 0+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
            cellEta = icol-AliEMCALGeoParams::fgkEMCALCols;
            cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
          }
          if(icol < AliEMCALGeoParams::fgkEMCALCols)
          {
            cellSM = 1+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
            cellEta = icol;
            cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
          }

          if( (icol < 0 || icol > AliEMCALGeoParams::fgkEMCALCols*2-1 ||
               irow < 0 || irow > AliEMCALGeoParams::fgkEMCALRows*5 - 1) //5*nRows+1/3*nRows //Count as bad "cells" out of EMCAL acceptance
             || (cu->GetEMCALChannelStatus(cellSM,cellEta,cellPhi,status)==1))  //Count as bad "cells" marked as bad in the DataBase
          {
            if     ( Radius(colC, rowC, icol, irow) < sqrSize ) coneBadCellsCoeff    += 1.;
            else if( icol>colC-sqrSize  &&  icol<colC+sqrSize ) phiBandBadCellsCoeff += 1 ;
	          else if( irow>rowC-sqrSize  &&  irow<rowC+sqrSize ) etaBandBadCellsCoeff += 1 ;
          }
        }
      }//end of cells loop
    }
    else AliWarning("Cluster with bad (eta,phi) in EMCal for energy density coeff calculation");

    if (coneCells > 0.)
    {
      //   printf("Energy density coneBadCellsCoeff= %.2f coneCells%.2f\n", coneBadCellsCoeff,coneCells);
      coneBadCellsCoeff = (coneCells-coneBadCellsCoeff)/coneCells;
      //  printf("coneBadCellsCoeff= %.2f\n", coneBadCellsCoeff);
    }
    
    if (phiBandCells > 0.)
    {
      // printf("Energy density phiBandBadCellsCoeff = %.2f phiBandCells%.2f\n", phiBandBadCellsCoeff,phiBandCells);
      phiBandBadCellsCoeff = (phiBandCells-phiBandBadCellsCoeff)/phiBandCells;
      // printf("phiBandBadCellsCoeff = %.2f\n", phiBandBadCellsCoeff);
    }
    
    if (etaBandCells > 0.)
    {
      //printf("Energy density etaBandBadCellsCoeff = %.2f etaBandCells%.2f\n", etaBandBadCellsCoeff,etaBandCells);
      etaBandBadCellsCoeff = (etaBandCells-etaBandBadCellsCoeff)/etaBandCells;
      // printf("etaBandBadCellsCoeff = %.2f\n",etaBandBadCellsCoeff);
    }
  }
}


//_________________________________________________________
/// Create histograms to be saved in output file and 
/// store them in outputContainer of the analysis class that calls this class.
//_________________________________________________________
TList * AliIsolationCut::GetCreateOutputObjects()
{    
  if ( !fHistoRanges )
  {
    AliFatal("Histogram ranges data base not initialized!");
    return 0x0;
  }
  
  fFillHistograms       = kTRUE;
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("IsolationCutBase") ; 
  outputContainer->SetOwner(kFALSE);
  
  Int_t  nptbins        = fHistoRanges->GetHistoPtBins();
  Float_t ptmax         = fHistoRanges->GetHistoPtMax();
  Float_t ptmin         = fHistoRanges->GetHistoPtMin();

  Int_t   nphibins      = fHistoRanges->GetHistoPhiBins();
  Int_t   netabins      = fHistoRanges->GetHistoEtaBins();
//  Float_t phimax        = fHistoRanges->GetHistoPhiMax();
//  Float_t etamax        = fHistoRanges->GetHistoEtaMax();
//  Float_t phimin        = fHistoRanges->GetHistoPhiMin();
//  Float_t etamin        = fHistoRanges->GetHistoEtaMin();
  
  Int_t  nptsumbins     = fHistoRanges->GetHistoNPtSumBins();
  Float_t ptsummax      = fHistoRanges->GetHistoPtSumMax();
  Float_t ptsummin      = fHistoRanges->GetHistoPtSumMin();
  Int_t  nptinconebins  = fHistoRanges->GetHistoNPtInConeBins();
  Float_t ptinconemax   = fHistoRanges->GetHistoPtInConeMax();
  Float_t ptinconemin   = fHistoRanges->GetHistoPtInConeMin();
  Float_t ptsumminUESub  = fHistoRanges->GetHistoPtSumSubMin();
  Float_t ptsummaxUESub  = fHistoRanges->GetHistoPtSumSubMax();
  Int_t  nptsumbinsUESub = fHistoRanges->GetHistoNPtSumSubBins();
  
  Int_t nrhobins = 100;
  Float_t rhomin = 0;
  Float_t rhomax = 200;

  //
  // For TH3 histograms, more coarse and not constant binning
  //
  TArrayD  ptBinsArray = fHistoRanges->GetHistoPtArr();
  TArrayD ptCBinsArray = fHistoRanges->GetHistoPtInConeArr();
  TArrayD cenBinsArray = fHistoRanges->GetHistoCentralityArr();
  TArrayD sumBinsArray = fHistoRanges->GetHistoPtSumArr();
  TArrayD sueBinsArray = fHistoRanges->GetHistoPtSumSubArr();
  TArrayD etaBinsArray = fHistoRanges->GetHistoEtaArr();
  TArrayD phiBinsArray = fHistoRanges->GetHistoPhiArr();

  TCustomBinning excBinning;
  excBinning.SetMinimum(0.0);
  excBinning.AddStep(2, 0.05);
  TArrayD excBinsArray;
  excBinning.CreateBinEdges(excBinsArray);

  TCustomBinning exc2Binning;
  exc2Binning.SetMinimum(0.0);
  exc2Binning.AddStep(4, 0.05);
  TArrayD exc2BinsArray;
  exc2Binning.CreateBinEdges(exc2BinsArray);
  
  TCustomBinning lfrBinning;
  lfrBinning.SetMinimum(0.0);
  lfrBinning.AddStep(5, 0.1);
  TArrayD lfrBinsArray;
  lfrBinning.CreateBinEdges(lfrBinsArray);

  TCustomBinning fraBinning;
  fraBinning.SetMinimum(0.0);
  fraBinning.AddStep(1, 0.05);// 20
  fraBinning.AddStep(2, 0.10); // 10
  fraBinning.AddStep(4, 0.20); // 10
  fraBinning.AddStep(8, 0.40); // 10
  TArrayD fraBinsArray;
  fraBinning.CreateBinEdges(fraBinsArray);
  
  TCustomBinning fueBinning;
  fueBinning.SetMinimum(-8);
  fueBinning.AddStep(-4, 0.40); // 10
  fueBinning.AddStep(-2, 0.20); // 10
  fueBinning.AddStep(-1, 0.10); // 10
  fueBinning.AddStep( 1, 0.05);// 20
  fueBinning.AddStep( 2, 0.10); // 10
  fueBinning.AddStep( 4, 0.20); // 10
  fueBinning.AddStep( 8, 0.40); // 10
  TArrayD fueBinsArray;
  fueBinning.CreateBinEdges(fueBinsArray);
  
  //
  // Titles strings
  //
  TString sParticle = ", x^{ 0,#pm}";
  if      ( fPartInCone == kOnlyNeutral )  sParticle = ", x^{0}";
  else if ( fPartInCone == kOnlyCharged )  sParticle = ", x^{#pm}";
  
  TString parTitleR   = Form("#it{R} = %2.2f%s",fConeSize,sParticle.Data());
  
  //
  // Create histograms
  //
  if ( fPartInCone == kNeutralAndCharged )
  {
    if ( !fFillHighMultHistograms )
    {
      // Pt in cone
      fhPtTrackInCone  = new TH2F
      ("hPtTrackInCone",
       Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtTrackInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtTrackInCone) ;

      fhPtClusterInCone  = new TH2F
      ("hPtClusterInCone",
       Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      fhPtClusterInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtClusterInCone) ;

      // Leading in cone
      if ( fFillLeadingHistograms ) 
      {
        fhConePtLeadTrack  = new TH2F
        ("hConeLeadPtTrack",
         Form("Track leading in isolation cone for #it{R} = %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhConePtLeadTrack->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
        fhConePtLeadTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConePtLeadTrack) ;
        
        fhConePtLeadCluster  = new TH2F
        ("hConeLeadPtCluster",
         Form("Cluster leading in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhConePtLeadCluster->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
        fhConePtLeadCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConePtLeadCluster) ;
        
        fhConePtLeadClustervsTrack   = new TH2F
        ("hConePtLeadClustervsTrack",
         Form("Track vs Cluster lead #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
        fhConePtLeadClustervsTrack->SetXTitle("#it{p}^{leading cluster}_{T} (GeV/#it{c})");
        fhConePtLeadClustervsTrack->SetYTitle("#it{p}^{leading track}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConePtLeadClustervsTrack) ;
        
        fhConePtLeadClusterTrackFrac   = new TH2F
        ("hConePtLeadClusterTrackFraction",
         Form(" #it{p}^{leading cluster}_{T}/#it{p}^{leading track}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,200,0,5);
        fhConePtLeadClusterTrackFrac->SetYTitle("#it{p}^{leading cluster}_{T}/ #it{p}^{leading track}_{T}");
        fhConePtLeadClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConePtLeadClusterTrackFrac) ;
      }
      
      // Sum pt in cone
      fhConeSumPtTrack  = new TH2F
      ("hConePtSumTrack",
       Form("Track #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrack) ;

      fhConeSumPtCluster  = new TH2F
      ("hConePtSumCluster",
       Form("Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtCluster) ;

      fhConeSumPtClusterMatched  = new TH2F
      ("hConePtSumClusterMatched",
       Form("CPV Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtClusterMatched->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterMatched->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterMatched) ;

      fhConeSumPtClusterMatchedFraction  = new TH2F
      ("hConePtSumClusterMatchedFraction",
       Form("CPV Cluster #Sigma #it{p}_{T} fraction to total, #it{R}=%2.2f",fConeSize),
       nptbins,ptmin,ptmax,200,0,1);
      fhConeSumPtClusterMatchedFraction->SetYTitle("Fraction");
      fhConeSumPtClusterMatchedFraction->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterMatchedFraction) ;

      fhConeSumPtClusterMatchedFraction3D  = new TH3F
      ("hConePtSumClusterMatchedFraction3D",
       Form("CPV Cluster #Sigma #it{p}_{T} fraction to total, #it{R}=%2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax,200,0,1);
      fhConeSumPtClusterMatchedFraction3D->SetZTitle("Fraction");
      fhConeSumPtClusterMatchedFraction3D->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterMatchedFraction3D->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterMatchedFraction3D) ;

      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtTrackTrigEtaPhi = new TH3F
        ("hConePtSumTrackTrigEtaPhi",
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s",parTitleR.Data()),
         etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
        fhConeSumPtTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtTrackTrigEtaPhi) ;
        
        fhConeSumPtClusterTrigEtaPhi  = new TH3F
        (Form("hConePtSumClusterTrigEtaPhi"),
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} clusters in isolation cone for %s",parTitleR.Data()),
         etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
        fhConeSumPtClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtClusterTrigEtaPhi) ;
      }
      
      // No need with perpendicular cones filled with tracks
      if ( fICMethod != kSumBkgSubIC )
      {
        fhConeSumPtClustervsTrack   = new TH2F
        ("hConePtSumClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClustervsTrack) ;

        fhConeSumPtClusterTrackFrac   = new TH2F
        ("hConePtSumClusterTrackFraction",
         Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track}, #it{R}=%2.2f",fConeSize),
         nptbins,ptmin,ptmax,200,0,5);
        fhConeSumPtClusterTrackFrac->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
        fhConeSumPtClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterTrackFrac) ;
      }

      if ( fICMethod >= kSumBkgSubIC )
      {
        fhConeSumPtUESubTrack  = new TH2F
        ("hConePtSumUESubTrack",
         Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction",fConeSize),
         nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhConeSumPtUESubTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtUESubTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtUESubTrack) ;

        fhConeSumPtUESubCluster  = new TH2F
        ("hConePtSumUESubCluster",
         Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction",fConeSize),
         nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhConeSumPtUESubCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtUESubCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtUESubCluster) ;

        if ( fFillEtaPhiHistograms )
        {
          fhConeSumPtUESubTrackTrigEtaPhi = new TH3F
          ("hConePtSumUESubTrackTrigEtaPhi",
           Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s",parTitleR.Data()),
           etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
           phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
          fhConeSumPtUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
          fhConeSumPtUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtUESubTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtUESubTrackTrigEtaPhi) ;
          
          fhConeSumPtUESubClusterTrigEtaPhi  = new TH3F
          (Form("hConePtSumUESubClusterTrigEtaPhi"),
           Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} clusters in isolation cone for %s",parTitleR.Data()),
           etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
           phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
          fhConeSumPtUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
          fhConeSumPtUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhConeSumPtUESubClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtUESubClusterTrigEtaPhi) ;
        }
        
        // No need in perp cones analysis
        if ( fICMethod != kSumBkgSubIC )
        {
          fhConeSumPtUESubClustervsTrack   = new TH2F
          ("hConePtSumUESubClustervsTrack",
           Form("Track vs Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f, UE correction",fConeSize),
           nptsumbinsUESub,ptsumminUESub,ptsummaxUESub,
           nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
          fhConeSumPtUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhConeSumPtUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtUESubClustervsTrack) ;

          fhConeSumPtUESubClusterTrackFrac   = new TH2F
          ("hConePtSumUESubClusterTrackFraction",
           Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track}, #it{R}=%2.2f, UE correction",fConeSize),
           nptbins,ptmin,ptmax,200,0,5);
          fhConeSumPtUESubClusterTrackFrac->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
          fhConeSumPtUESubClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtUESubClusterTrackFrac) ;
        } // per cones
      } // UE sub
    }
    else
    {
      fhPtTrackInConeCent  = new TH3F
      ("hPtTrackInConeCent",
       Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f",fConeSize),
       //nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
       cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
      fhPtTrackInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtTrackInConeCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtTrackInConeCent->SetZTitle("Centrality (%)");
      outputContainer->Add(fhPtTrackInConeCent) ;

      fhPtClusterInConeCent  = new TH3F
      ("hPtClusterInConeCent",
       Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",fConeSize),
       //nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
       cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
      fhPtClusterInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
      fhPtClusterInConeCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fhPtClusterInConeCent->SetZTitle("Centrality (%)");
      outputContainer->Add(fhPtClusterInConeCent) ;
    }
  }

  if ( !fFillHighMultHistograms )
  {
    fhPtInCone  = new TH2F
    ("hPtInCone",
     Form("#it{p}_{T} of clusters and tracks in isolation cone for %s",parTitleR.Data()),
     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInCone) ;
    
    if ( fFillLeadingHistograms ) 
    {
      fhConePtLead  = new TH2F
      ("hConePtLead",
       Form("Track or Cluster  leading #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhConePtLead->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
      fhConePtLead->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConePtLead) ;
    }
 
    fhConeSumPt  = new TH2F
    ("hConePtSum",
     Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPt->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPt->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPt) ;
    
    if ( fICMethod >= kSumBkgSubIC )
    {
      fhConeSumPtUESub  = new TH2F
      ("hConePtSumUESub",
       Form("Track and/or Cluster #Sigma #it{p}_{T} in #it{R} = %2.2f, after UE correction",fConeSize),
       nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESub->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUESub) ;
    }
    
    if ( fFillEtaPhiHistograms )
    {
      fhConeSumPtTrigEtaPhi  = new TH3F
      ("hConePtSumTrigEtaPhi",
       Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} in isolation cone for %s",parTitleR.Data()),
       etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
       phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
       sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
      fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
      
      if ( fICMethod >= kSumBkgSubIC )
      {
        fhConeSumPtUESubTrigEtaPhi  = new TH3F
        ("hConePtSumUESubTrigEtaPhi",
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} for %s, after UE correction",parTitleR.Data()),
         etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1,  sueBinsArray.GetArray());
        fhConeSumPtUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtUESubTrigEtaPhi) ;
      }
    }
  }
  else
  {
    fhPtInConeCent  = new TH3F
    ("hPtInCone",
     Form("#it{p}_{T} of clusters and tracks in isolation cone for %s",parTitleR.Data()),
     //nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
      ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
     ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
     cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
    fhPtInConeCent->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInConeCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhPtInConeCent->SetZTitle("Centrality (%)");
    outputContainer->Add(fhPtInConeCent) ;
  }
  
  if ( fPartInCone != kOnlyCharged && fFillEtaPhiHistograms )
  {
    fhEtaPhiCluster = new TH2F
    ("hEtaPhiCluster",
     Form("#eta vs #varphi of all clusters"),
     netabins,-1,1,nphibins,0,TMath::TwoPi());
    fhEtaPhiCluster->SetXTitle("#eta");
    fhEtaPhiCluster->SetYTitle("#varphi (rad)");
    outputContainer->Add(fhEtaPhiCluster) ;
    
    fhEtaPhiInConeCluster = new TH2F
    ("hEtaPhiInConeCluster",
     Form("#eta vs #varphi of clusters in cone for #it{R} =  %2.2f",fConeSize),
     netabins,-1,1,nphibins,0,TMath::TwoPi());
    fhEtaPhiInConeCluster->SetXTitle("#eta");
    fhEtaPhiInConeCluster->SetYTitle("#varphi (rad)");
    outputContainer->Add(fhEtaPhiInConeCluster) ;

    fhDeltaEtaPhiInConeCluster = new TH2F
    ("hDeltaEtaPhiInConeCluster",
     Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-cluster} of Clusters in cone for #it{R} =  %2.2f",fConeSize),
     50,-fConeSize*1.2,fConeSize*1.2,60,-fConeSize*1.05,fConeSize*1.05);
    fhDeltaEtaPhiInConeCluster->SetXTitle("#Delta #eta_{trig-cluster}");
    fhDeltaEtaPhiInConeCluster->SetYTitle("#Delta #varphi_{trig-cluster}");
    outputContainer->Add(fhDeltaEtaPhiInConeCluster) ;

    fhTriggerEtaPhiInConeClusterPt = new TH3F
    ("hEtaPhiInConeClusterPt",
     Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in cone for #it{R} = %2.2f",fConeSize),
     etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
     phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
     ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
    fhTriggerEtaPhiInConeClusterPt->SetXTitle("#eta");
    fhTriggerEtaPhiInConeClusterPt->SetYTitle("#varphi");
    fhTriggerEtaPhiInConeClusterPt->SetZTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhTriggerEtaPhiInConeClusterPt) ;
  }
  
  // UE bands
  if ( fPartInCone != kOnlyCharged && fICMethod >= kSumBkgSubIC )
  {
    if ( !fFillHighMultHistograms )
    {
      if ( fFillLeadingVsUESubHisto )
      {
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadClusterInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt->SetYTitle("#it{p}^{lead}_{cluster} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt) ;

        fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt->SetYTitle("#it{p}^{lead UE}_{cluster} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt) ;

        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead->SetYTitle("#it{p}^{lead}_{cluster} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead) ;

        fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         lfrBinsArray.GetSize() - 1, lfrBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt->SetYTitle("#it{p}^{lead}_{track} Cone / UE ");
        fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt) ;
      }

      fhConeSumPtUESubClusterCutMax  = new TH2F
      ("hConePtSumUESubClusterCutMax",
       Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on UE leading pT",fConeSize),
       nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtUESubClusterCutMax->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESubClusterCutMax->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUESubClusterCutMax) ;

      fhConeSumPtUESubClusterCutLeadFactor  = new TH2F
      ("hConePtSumUESubClusterCutLeadFactor",
       Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on UE leading pT",fConeSize),
       nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtUESubClusterCutLeadFactor->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESubClusterCutLeadFactor->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUESubClusterCutLeadFactor) ;

      if ( fICMethod == kSumBkgSubEtaBandIC )
      {
        fhEtaBandClusterPt  = new TH2F
        ("hEtaBandClusterPt",
         Form("Clusters in #eta band out of cone for #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhEtaBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhEtaBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandClusterPt) ;
      }

      if ( fICMethod == kSumBkgSubPhiBandIC )
      {
        fhPhiBandClusterPt  = new TH2F
        ("hPhiBandClusterPt",
         Form("Clusters in #varphi band out of cone for #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPhiBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhPhiBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandClusterPt) ;
      }
    }
    
    if ( fFillEtaPhiHistograms )
    {
      if ( fICMethod == kSumBkgSubEtaBandIC )
      {
        fhEtaBandClusterEtaPhi  = new TH2F
        ("hEtaBandClusterEtaPhi",
         Form("#eta vs #varphi of clusters in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandClusterEtaPhi->SetXTitle("#eta");
        fhEtaBandClusterEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaBandClusterEtaPhi) ;

        fhEtaBandClusterDeltaEtaPhi = new TH2F
        ("hEtaBandClusterDeltaEtaPhi",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-cluster} of clusters in #varphi band for #it{R} =  %2.2f",fConeSize),
         70,-0.7,0.7, 60,-fConeSize*1.05,fConeSize*1.05);
        fhEtaBandClusterDeltaEtaPhi->SetXTitle("#Delta #eta_{trig-cluster}");
        fhEtaBandClusterDeltaEtaPhi->SetYTitle("#Delta #varphi_{trig-cluster}");
        outputContainer->Add(fhEtaBandClusterDeltaEtaPhi) ;

        fhEtaBandClusterPtTriggerEtaPhi = new TH3F
        ("hEtaBandClusterPtTriggerEtaPhi",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} clusters in #varphi band for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhEtaBandClusterPtTriggerEtaPhi->SetXTitle("#eta");
        fhEtaBandClusterPtTriggerEtaPhi->SetYTitle("#varphi");
        fhEtaBandClusterPtTriggerEtaPhi->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandClusterPtTriggerEtaPhi) ;
      }

      if ( fICMethod == kSumBkgSubPhiBandIC )
      {
        fhPhiBandClusterEtaPhi  = new TH2F
        ("hPhiBandClusterEtaPhi",
         Form("#eta vs #varphi of clusters in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandClusterEtaPhi->SetXTitle("#eta");
        fhPhiBandClusterEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhPhiBandClusterEtaPhi) ;

        fhPhiBandClusterDeltaEtaPhi = new TH2F
        ("hPhiBandClusterDeltaEtaPhi",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-cluster} of clusters in #varphi band for #it{R} =  %2.2f",fConeSize),
         60,-fConeSize*1.05,fConeSize*1.05, 120,-TMath::Pi()*1.1,TMath::Pi()*1.1);
        fhPhiBandClusterDeltaEtaPhi->SetXTitle("#Delta #eta_{trig-track}");
        fhPhiBandClusterDeltaEtaPhi->SetYTitle("#Delta #varphi_{trig-track}");
        outputContainer->Add(fhPhiBandClusterDeltaEtaPhi) ;

        fhPhiBandClusterPtTriggerEtaPhi = new TH3F
        ("hPhiBandClusterPtTriggerEtaPhi",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} clusters in #varphi band for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhPhiBandClusterPtTriggerEtaPhi->SetXTitle("#eta");
        fhPhiBandClusterPtTriggerEtaPhi->SetYTitle("#varphi");
        fhPhiBandClusterPtTriggerEtaPhi->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandClusterPtTriggerEtaPhi) ;
      }
    }
    
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      if ( !fFillHighMultHistograms )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUECluster  = new TH2F
          ("hConePtSumEtaBandUECluster",
           "#Sigma cluster #it{p}_{T} in UE Eta Band",
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtEtaBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaBandUECluster) ;
        }

        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUECluster  = new TH2F
          ("hConePtSumPhiBandUECluster",
           "#Sigma cluster #it{p}_{T} UE Phi Band",
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPhiBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiBandUECluster) ;
        }

        if ( fFillEtaPhiHistograms )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
          {
            fhConeSumPtEtaBandUEClusterTrigEtaPhi  = new TH3F
            ("hConePtSumEtaBandUEClusterTrigEtaPhi",
             "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE Eta Band",
             etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
             phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
            fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
            fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
            fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhi) ;
          }

          if ( fICMethod == kSumBkgSubPhiBandIC )
          {
            fhConeSumPtPhiBandUEClusterTrigEtaPhi  = new TH3F
            ("hConePtSumPhiBandUEClusterTrigEtaPhi",
             "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} UE Phi Band",
             etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
             phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
            fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
            fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
            fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhi) ;
          }
        }
        
        // Subtraction
      
        fhConeSumPtUEBandNormCluster  = new TH2F
        ("hConeSumPtUEBandNormCluster",
         Form("Clusters #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtUEBandNormCluster->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
        fhConeSumPtUEBandNormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtUEBandNormCluster) ;

        fhConeRhoUEBandCluster  = new TH2F
        ("hConeRhoUEBandCluster",
         Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandCluster->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandCluster) ;

        fhConeRhoUEBandClusterCutMax  = new TH2F
        ("hConeRhoUEBandClusterCutMax",
         Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandClusterCutMax->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandClusterCutMax->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandClusterCutMax) ;

        fhConeRhoUEBandClusterCutLeadFactor  = new TH2F
        ("hConeRhoUEBandClusterCutLeadFactor",
         Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandClusterCutLeadFactor->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandClusterCutLeadFactor->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandClusterCutLeadFactor) ;

        fhConeSumPtClusterSubVsNoSub = new TH2F
        ("hConeSumPtClusterSubVsNoSub",
         Form("Cluster #Sigma #it{p}_{T} in cone before vs after UE bkg sub, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhConeSumPtClusterSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtClusterSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtClusterSubVsNoSub);

        fhConeSumPtVSUEClusterBand  = new TH2F
        ("hConeSumPtVSUEClusterBand",
         Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in normalized UE band for cluster, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtVSUEClusterBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUEClusterBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUEClusterBand);
      }
    }
  } // clusters
  
  if ( fPartInCone != kOnlyNeutral && fFillEtaPhiHistograms )
  {
    fhEtaPhiTrack = new TH2F
    ("hEtaPhiTrack",
     Form("#eta vs #varphi of all Tracks"),
     netabins,-1,1,nphibins,0,TMath::TwoPi());
    fhEtaPhiTrack->SetXTitle("#eta");
    fhEtaPhiTrack->SetYTitle("#varphi (rad)");
    outputContainer->Add(fhEtaPhiTrack) ;
    
    fhEtaPhiInConeTrack = new TH2F
    ("hEtaPhiInConeTrack",
     Form("#eta vs #varphi of Tracks in cone for #it{R} =  %2.2f",fConeSize),
     netabins,-1,1,nphibins,0,TMath::TwoPi());
    fhEtaPhiInConeTrack->SetXTitle("#eta");
    fhEtaPhiInConeTrack->SetYTitle("#varphi");
    outputContainer->Add(fhEtaPhiInConeTrack) ;

    fhDeltaEtaPhiInConeTrack = new TH2F
    ("hDeltaEtaPhiInConeTrack",
     Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-track} of Tracks in cone for #it{R} =  %2.2f",fConeSize),
     50,-fConeSize*1.2,fConeSize*1.2,60,-fConeSize*1.05,fConeSize*1.05);
    fhDeltaEtaPhiInConeTrack->SetXTitle("#Delta #eta_{trig-track}");
    fhDeltaEtaPhiInConeTrack->SetYTitle("#Delta #varphi_{trig-track}");
    outputContainer->Add(fhDeltaEtaPhiInConeTrack) ;

    fhTriggerEtaPhiInConeTrackPt = new TH3F
    ("hEtaPhiInConeTrackPt",
     Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in cone for #it{R} = %2.2f",fConeSize),
     etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
     phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
     ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
    fhTriggerEtaPhiInConeTrackPt->SetXTitle("#eta");
    fhTriggerEtaPhiInConeTrackPt->SetYTitle("#varphi");
    fhTriggerEtaPhiInConeTrackPt->SetZTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhTriggerEtaPhiInConeTrackPt) ;
  }
  
  // UE subtraction
  if ( fPartInCone != kOnlyNeutral && fICMethod >= kSumBkgSubIC )
  {
    // UE in perpendicular cones
    //
    if ( fICMethod == kSumBkgSubIC )
    {
      if ( !fFillHighMultHistograms )
      {
        fhPtInPerpCone  = new TH2F
        ("hPtInPerpCone",
         Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPtInPerpCone->SetYTitle("#it{p}_{T in #perp cone} (GeV/#it{c})");
        fhPtInPerpCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPtInPerpCone) ;

        fhPerpConeSumPt  = new TH2F
        ("hPerpConePtSum",
         Form("#Sigma #it{p}_{T} in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPerpConeSumPt->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        fhPerpConeSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeSumPt) ;

        fhPerpCone1SumPt  = new TH2F
        ("hPerpCone1PtSum",
         Form("#Sigma #it{p}_{T} in 1 isolation cone at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhPerpCone1SumPt->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        fhPerpCone1SumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpCone1SumPt) ;

        fhPerpCone1SumPtUESub  = new TH2F
        ("hPerpCone1PtSumUESub",
         Form("#Sigma #it{p}_{T} - 1 #perp cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhPerpCone1SumPtUESub->SetYTitle("#Sigma #it{p}_{T, UE sub}^{iso} (GeV/#it{c})");
        fhPerpCone1SumPtUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpCone1SumPtUESub) ;

        fhPerpConeRho  = new TH2F
        ("hPerpConeRho",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhPerpConeRho->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRho->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeRho) ;

        fhPerpConeRhoCutMax  = new TH2F
        ("hPerpConeRhoCutMax",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhPerpConeRhoCutMax->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRhoCutMax->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeRhoCutMax) ;

        fhPerpConeRhoCutLeadFactor  = new TH2F
        ("hPerpConeRhoCutLeadFactor",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhPerpConeRhoCutLeadFactor->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRhoCutLeadFactor->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeRhoCutLeadFactor) ;

        fhConeSumPtVSPerpCone = new TH2F
        ("hConeSumPtVSPerpCone",
         Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in 2 isolation cones at #pm 45 degree #varphi, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtVSPerpCone->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtVSPerpCone->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSPerpCone);

        fhPerpConeSumPtTrackSubVsNoSub = new TH2F
        ("hPerpConeSumPtTrackSubVsNoSub",
         Form("Track #Sigma #it{p}_{T} in cone before and after UE bkg sub, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhPerpConeSumPtTrackSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPtTrackSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
        outputContainer->Add(fhPerpConeSumPtTrackSubVsNoSub);

        if ( fFillEtaPhiHistograms )
        {
          fhPerpConeSumPtTrigEtaPhi = new TH3F
          ("hPerpConeSumPtTrigEtaPhi",
           "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in #perp cones",
           etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
           phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
          fhPerpConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
          fhPerpConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
          fhPerpConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhPerpConeSumPtTrigEtaPhi) ;
        }
      }
      
      if ( fFillEtaPhiHistograms )
      {
        fhEtaPhiInPerpCone = new TH2F
        ("hEtaPhiInPerpCone",
         Form("#eta vs #varphi of all Tracks in #perp cones"),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInPerpCone->SetXTitle("#eta");
        fhEtaPhiInPerpCone->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiInPerpCone) ;

        fhDeltaEtaPhiInPerpCone = new TH2F
        ("hDeltaEtaPhiInPerpConeTrack",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-track} of Tracks in #perp cones for #it{R} =  %2.2f",fConeSize),
         50,-fConeSize*1.2,fConeSize*1.2,120,-TMath::Pi()-fConeSize*1.1,TMath::Pi()+fConeSize*1.1);
        fhDeltaEtaPhiInPerpCone->SetXTitle("#Delta #eta_{trig-track}");
        fhDeltaEtaPhiInPerpCone->SetYTitle("#Delta #varphi_{trig-track}");
        outputContainer->Add(fhDeltaEtaPhiInPerpCone) ;

        fhTriggerEtaPhiInPerpConeTrackPt = new TH3F
        ("hTriggerEtaPhiInPerpConeTrackPt",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in #perp cones for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhTriggerEtaPhiInPerpConeTrackPt->SetXTitle("#eta");
        fhTriggerEtaPhiInPerpConeTrackPt->SetYTitle("#varphi");
        fhTriggerEtaPhiInPerpConeTrackPt->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhTriggerEtaPhiInPerpConeTrackPt) ;
      }
    } // perpendicular
    
    if ( fICMethod == kSumBkgSubJetRhoIC && !fFillHighMultHistograms )
    {
      fhJetRhoSumPt  = new TH2F
      ("hJetRhoPtSum",
       Form("Jet #rho #it{R}^{2} #pi, #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhJetRhoSumPt->SetYTitle("#rho #it{R}^{2} #pi (GeV/#it{c})");
      fhJetRhoSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhJetRhoSumPt) ;

      fhJetRho  = new TH2F
      ("hJetRho",Form("Jet #rho, #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
      fhJetRho->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
      fhJetRho->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhJetRho) ;

      fhConeSumPtVSJetRho = new TH2F
      ("hConeSumPtVSJetRho",
       Form("Jet #rho #it{R}^{2} #pi versus #Sigma #it{p}_{T}, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtVSJetRho->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtVSJetRho->SetYTitle("Jet #rho #it{R}^{2} #pi (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSJetRho);

      fhJetRhoSumPtTrackSubVsNoSub = new TH2F
      ("hJetRhoSumPtTrackSubVsNoSub",
       Form("Track #Sigma #it{p}_{T} in cone before and after UE bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhJetRhoSumPtTrackSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhJetRhoSumPtTrackSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
      outputContainer->Add(fhJetRhoSumPtTrackSubVsNoSub);

      if ( fFillEtaPhiHistograms )
      {
        fhJetRhoSumPtTrigEtaPhi = new TH3F
        ("hJetRhoSumPtTrigEtaPhi",
         "Trigger #eta vs #varphi, Jet #rho #it{R}^{2} #pi",
         etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
        fhJetRhoSumPtTrigEtaPhi->SetZTitle("Jet #rho #it{R}^{2} #pi (GeV/#it{c})");
        fhJetRhoSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhJetRhoSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhJetRhoSumPtTrigEtaPhi) ;
      }
    }

    // UE bands
    if ( !fFillHighMultHistograms  )
    {
      if ( fFillLeadingVsUESubHisto )
      {
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadTrackInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt->SetYTitle("#it{p}^{lead}_{track} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt) ;

        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead->SetYTitle("#it{p}^{lead}_{track} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead) ;

        fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt->SetYTitle("#it{p}^{lead UE}_{track} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt) ;

        fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt = new TH3F
        ("hTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt",Form("%s",parTitleR.Data()),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         lfrBinsArray.GetSize() - 1, lfrBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
        fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt->SetYTitle("#it{p}^{lead}_{track} Cone / UE ");
        fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
        outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt) ;
      }

      fhConeSumPtUESubTrackCutMax  = new TH2F
      ("hConePtSumUESubTrackCutMax",
       Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on UE leading pT",fConeSize),
       nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtUESubTrackCutMax->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESubTrackCutMax->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUESubTrackCutMax) ;

      fhConeSumPtUESubTrackCutLeadFactor  = new TH2F
      ("hConePtSumUESubTrackCutLeadFactor",
       Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on UE leading pT",fConeSize),
       nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtUESubTrackCutLeadFactor->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESubTrackCutLeadFactor->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUESubTrackCutLeadFactor) ;

      if ( fICMethod == kSumBkgSubEtaBandIC )
      {
        fhEtaBandTrackPt  = new TH2F
        ("hEtaBandTrackPt",
         Form("Tracks in #eta band out of cone #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhEtaBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhEtaBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandTrackPt) ;
      }

      if ( fICMethod == kSumBkgSubPhiBandIC )
      {
        fhPhiBandTrackPt  = new TH2F
        ("hPhiBandTrackPt",
         Form("Tracks in #varphi band out of cone #it{R} = %2.2f and half TPC, #pm #pi",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPhiBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhPhiBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandTrackPt) ;
      }

      if ( fICMethod == kSumBkgSubPerpBandIC )
      {
        fhPerpBandTrackPt  = new TH2F
        ("hPerpBandTrackPt",
         Form("Tracks in #perp #eta band out of cone #it{R} = %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
        fhPerpBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
        fhPerpBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
        outputContainer->Add(fhPerpBandTrackPt) ;
      }
    }
    
    if ( fFillEtaPhiHistograms )
    {
      if ( fICMethod == kSumBkgSubEtaBandIC )
      {
        fhEtaBandTrackEtaPhi  = new TH2F
        ("hEtaBandTrackEtaPhi",
         Form("#eta vs #varphi of tracks in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaBandTrackEtaPhi->SetXTitle("#eta");
        fhEtaBandTrackEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaBandTrackEtaPhi) ;

        fhEtaBandTrackDeltaEtaPhi = new TH2F
        ("hEtaBandTrackDeltaEtaPhi",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-track} of tracks in #eta band for #it{R} =  %2.2f",fConeSize),
         100,-1,1, 60,-fConeSize*1.05,fConeSize*1.05);
        fhEtaBandTrackDeltaEtaPhi->SetXTitle("#Delta #eta_{trig-track}");
        fhEtaBandTrackDeltaEtaPhi->SetYTitle("#Delta #varphi_{trig-track}");
        outputContainer->Add(fhEtaBandTrackDeltaEtaPhi) ;

        fhEtaBandTrackPtTriggerEtaPhi = new TH3F
        ("hEtaBandTrackPtTriggerEtaPhi",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in #eta band for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhEtaBandTrackPtTriggerEtaPhi->SetXTitle("#eta");
        fhEtaBandTrackPtTriggerEtaPhi->SetYTitle("#varphi");
        fhEtaBandTrackPtTriggerEtaPhi->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandTrackPtTriggerEtaPhi) ;
      }

      if ( fICMethod == kSumBkgSubPhiBandIC )
      {
        fhPhiBandTrackEtaPhi  = new TH2F
        ("hPhiBandTrackEtaPhi",
         Form("#eta vs #varphi of tracks in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPhiBandTrackEtaPhi->SetXTitle("#eta");
        fhPhiBandTrackEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhPhiBandTrackEtaPhi) ;

        fhPhiBandTrackDeltaEtaPhi = new TH2F
        ("hPhiBandTrackDeltaEtaPhi",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-track} of tracks in #varphi band for #it{R} =  %2.2f",fConeSize),
         60,-fConeSize*1.05,fConeSize*1.05, 120,-TMath::Pi()*1.1,TMath::Pi()*1.1);
        fhPhiBandTrackDeltaEtaPhi->SetXTitle("#Delta #eta_{trig-track}");
        fhPhiBandTrackDeltaEtaPhi->SetYTitle("#Delta #varphi_{trig-track}");
        outputContainer->Add(fhPhiBandTrackDeltaEtaPhi) ;

        fhPhiBandTrackPtTriggerEtaPhi = new TH3F
        ("hPhiBandTrackPtTriggerEtaPhi",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in #varphi band for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhPhiBandTrackPtTriggerEtaPhi->SetXTitle("#eta");
        fhPhiBandTrackPtTriggerEtaPhi->SetYTitle("#varphi");
        fhPhiBandTrackPtTriggerEtaPhi->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandTrackPtTriggerEtaPhi) ;
      }

      if ( fICMethod == kSumBkgSubPerpBandIC )
      {
        fhPerpBandTrackEtaPhi  = new TH2F
        ("hPerpBandTrackEtaPhi",
         Form("#eta vs #varphi of tracks in #perp #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhPerpBandTrackEtaPhi->SetXTitle("#eta");
        fhPerpBandTrackEtaPhi->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhPerpBandTrackEtaPhi) ;

        fhPerpBandTrackDeltaEtaPhi = new TH2F
        ("hPerpBandTrackDeltaEtaPhi",
         Form("#Delta #eta_{trig-track} vs #Delta #varphi_{trig-track} of Tracks in #perp bands for #it{R} =  %2.2f",fConeSize),
         50,-1,1, 120,-TMath::Pi()-fConeSize*1.1,TMath::Pi()+fConeSize*1.1);
        fhPerpBandTrackDeltaEtaPhi->SetXTitle("#Delta #eta_{trig-track}");
        fhPerpBandTrackDeltaEtaPhi->SetYTitle("#Delta #varphi_{trig-track}");
        outputContainer->Add(fhPerpBandTrackDeltaEtaPhi) ;

        fhPerpBandTrackPtTriggerEtaPhi = new TH3F
        ("hPerpBandTrackPtTriggerEtaPhi",
         Form("#eta vs #varphi of trigger vs #it{p}_{T} tracks in #perp bands for #it{R} = %2.2f",fConeSize),
         etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
         ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray());
        fhPerpBandTrackPtTriggerEtaPhi->SetXTitle("#eta");
        fhPerpBandTrackPtTriggerEtaPhi->SetYTitle("#varphi");
        fhPerpBandTrackPtTriggerEtaPhi->SetZTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhPerpBandTrackPtTriggerEtaPhi) ;
      }
    }
    
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      if ( !fFillHighMultHistograms )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUETrack  = new TH2F
          ("hConePtSumEtaBandUETrack",
           "#Sigma track #it{p}_{T} in UE #eta Band",
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtEtaBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtEtaBandUETrack) ;
        }

        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUETrack  = new TH2F
          ("hConePtSumPhiBandUETrack",
           "#Sigma track #it{p}_{T} in UE #varphi Band",
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPhiBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPhiBandUETrack) ;
        }

        if ( fICMethod == kSumBkgSubPerpBandIC )
        {
          fhConeSumPtPerpBandUETrack  = new TH2F
          ("hConePtSumPerpBandUETrack",
           "#Sigma track #it{p}_{T} in UE #perp #eta Band",
           nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
          fhConeSumPtPerpBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPerpBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          outputContainer->Add(fhConeSumPtPerpBandUETrack) ;
        }

        if ( fFillEtaPhiHistograms )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
          {
            fhConeSumPtEtaBandUETrackTrigEtaPhi  = new TH3F
            ("hConePtSumEtaBandUETrackTrigEtaPhi",
             "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE #eta Band",
             etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
             phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
            fhConeSumPtEtaBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
            fhConeSumPtEtaBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
            fhConeSumPtEtaBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhi) ;
          }

          if ( fICMethod == kSumBkgSubPhiBandIC )
          {
            fhConeSumPtPhiBandUETrackTrigEtaPhi  = new TH3F
            ("hConePtSumPhiBandUETrackTrigEtaPhi",
             "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE #varphi Band",
             etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
             phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
            fhConeSumPtPhiBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
            fhConeSumPtPhiBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
            fhConeSumPtPhiBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhi) ;
          }

          if ( fICMethod == kSumBkgSubPerpBandIC )
          {
            fhConeSumPtPerpBandUETrackTrigEtaPhi  = new TH3F
            ("hConePtSumPerpBandUETrackTrigEtaPhi",
             "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE #perp #eta Band",
             etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
             phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
            fhConeSumPtPerpBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
            fhConeSumPtPerpBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
            fhConeSumPtPerpBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtPerpBandUETrackTrigEtaPhi) ;
          }
        }

        // Subtraction

        fhConeSumPtUEBandNormTrack  = new TH2F
        ("hConeSumPtUEBandNormTrack",
         Form("Tracks #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtUEBandNormTrack->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
        fhConeSumPtUEBandNormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtUEBandNormTrack) ;

        fhConeRhoUEBandTrack  = new TH2F
        ("hConeRhoUEBandTrack",
         Form("Tracks #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandTrack->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandTrack) ;

        fhConeRhoUEBandTrackCutMax  = new TH2F
        ("hConeRhoUEBandTrackCutMax",
         Form("Tracks #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandTrackCutMax->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandTrackCutMax->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandTrackCutMax) ;

        fhConeRhoUEBandTrackCutLeadFactor  = new TH2F
        ("hConeRhoUEBandTrackCutLeadFactor",
         Form("Tracks #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
         nptbins,ptmin,ptmax,nrhobins,rhomin,rhomax);
        fhConeRhoUEBandTrackCutLeadFactor->SetYTitle("#rho (GeV/#it{c} rad^{-2})");
        fhConeRhoUEBandTrackCutLeadFactor->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        outputContainer->Add(fhConeRhoUEBandTrackCutLeadFactor) ;

        fhConeSumPtTrackSubVsNoSub = new TH2F
        ("hConeSumPtTrackSubVsNoSub",
         Form("Track #Sigma #it{p}_{T} in cone before and after UE bkg sub, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
        fhConeSumPtTrackSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrackSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtTrackSubVsNoSub);

        fhConeSumPtVSUETrackBand  = new TH2F
        ("hConeSumPtVSUETrackBand",
         Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in normalized UE band for tracks, R=%2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhConeSumPtVSUETrackBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
        fhConeSumPtVSUETrackBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtVSUETrackBand);
      }
    }
  } // charged

  // UE subtraction
  
  if ( fPartInCone == kNeutralAndCharged && fICMethod >= kSumBkgSubEtaBandIC && !fFillHighMultHistograms )
  {
    fhConeSumPtUEBandSubClustervsTrack   = new TH2F
    ("hConePtSumUEBandSubClustervsTrack",
     Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptsumbinsUESub,ptsumminUESub,ptsummaxUESub,
     nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
    fhConeSumPtUEBandSubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
    fhConeSumPtUEBandSubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtUEBandSubClustervsTrack) ;
    
    fhBandClustervsTrack   = new TH2F
    ("hBandClustervsTrack",
     Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
    fhBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
    fhBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
    outputContainer->Add(fhBandClustervsTrack) ;
    
    fhBandNormClustervsTrack   = new TH2F
    ("hBandNormClustervsTrack",
     Form("Track vs Cluster Normalized #Sigma #it{p}_{T} in #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
    fhBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
    fhBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
    outputContainer->Add(fhBandNormClustervsTrack) ;
  } // UE, neutral + charged

  if ( fMakeConeExcessCorr && fFillFractionExcessHistograms )
  {
    if ( fPartInCone != kOnlyCharged )
    {
      fhFractionClusterOutConeEta  = new TH2F
      ("hFractionClusterOutConeEta",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta acceptance",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionClusterOutConeEta->SetYTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConeEta) ;

      fhFractionClusterOutConeEtaTrigEtaPhi  = new TH3F
      ("hFractionClusterOutConeEtaTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       //netabins,etamin,etamax,nphibins,phimin,phimax);
       etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
       phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;

      fhFractionClusterOutConePhi  = new TH2F
      ("hFractionClusterOutConePhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #varphi acceptance",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionClusterOutConePhi->SetYTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConePhi) ;

      fhFractionClusterOutConePhiTrigEtaPhi  = new TH3F
      ("hFractionClusterOutConePhiTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #varphi acceptance, in trigger #eta-#varphi ",fConeSize),
       etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
       phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;

      fhFractionClusterOutConeEtaPhi  = new TH2F
      ("hFractionClusterOutConeEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta x #varphi acceptance",fConeSize),
         ptBinsArray.GetSize() - 1,   ptBinsArray.GetArray(),
       exc2BinsArray.GetSize() - 1, exc2BinsArray.GetArray());
      fhFractionClusterOutConeEtaPhi->SetYTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaPhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConeEtaPhi) ;

      fhFractionClusterOutConeEtaPhiTrigEtaPhi  = new TH3F
      ("hFractionClusterOutConeEtaPhiTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta  x #varphi acceptance, in trigger #eta-#varphi ",fConeSize),
        etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
        phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
       exc2BinsArray.GetSize() - 1, exc2BinsArray.GetArray());
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetZTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConeEtaPhiTrigEtaPhi) ;
    }

    if ( fPartInCone != kOnlyNeutral )
    {
      fhFractionTrackOutConeEta  = new TH2F
      ("hFractionTrackOutConeEta",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of tracks #eta acceptance",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionTrackOutConeEta->SetYTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionTrackOutConeEta) ;

      fhFractionTrackOutConeEtaTrigEtaPhi  = new TH3F
      ("hFractionTrackOutConeEtaTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of tracks #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       etaBinsArray.GetSize() - 1, etaBinsArray.GetArray(),
       phiBinsArray.GetSize() - 1, phiBinsArray.GetArray(),
       excBinsArray.GetSize() - 1, excBinsArray.GetArray());
      fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
    }
  }

  if ( fFillHighMultHistograms )
  {
    fhConeSumPtCent  = new TH3F
    ("hConePtSumCent",
     Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",fConeSize),
      ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
     sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
     cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
    fhConeSumPtCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    fhConeSumPtCent->SetZTitle("Centrality (%)");
    outputContainer->Add(fhConeSumPtCent) ;
   
    if ( fFillEtaPhiHistograms )
    {
      fhConeSumPtTrigEtaPhiCent   = new TH3F*[fNCentBins] ;
      for(Int_t icen = 0; icen < fNCentBins; icen++)
      {
        fhConeSumPtTrigEtaPhiCent[icen]  = new TH3F
        (Form("hConePtSumTrigEtaPhi_Cen%d",icen),
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
         etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
         phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
        fhConeSumPtTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
        fhConeSumPtTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtTrigEtaPhiCent[icen]) ;
      }
    }
    
    if ( fPartInCone == kNeutralAndCharged )
    {
      fhConeSumPtTrackCent  = new TH3F
      ("hConePtSumTrackCent",
       Form("Track #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
       cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
      fhConeSumPtTrackCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      fhConeSumPtTrackCent->SetZTitle("Centrality (%)");
      outputContainer->Add(fhConeSumPtTrackCent) ;
      
      fhConeSumPtClusterCent  = new TH3F
      ("hConePtSumClusterCent",
       Form("Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
       cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
      fhConeSumPtClusterCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      fhConeSumPtClusterCent->SetZTitle("Centrality (%)");
      outputContainer->Add(fhConeSumPtClusterCent) ;
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtTrackTrigEtaPhiCent   = new TH3F*[fNCentBins] ;
        fhConeSumPtClusterTrigEtaPhiCent = new TH3F*[fNCentBins] ;
        for(Int_t icen = 0; icen < fNCentBins; icen++)
        {
          fhConeSumPtTrackTrigEtaPhiCent[icen]  = new TH3F
          (Form("hConePtSumTrackTrigEtaPhi_Cen%d",icen),
           Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
           etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
           phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
          fhConeSumPtTrackTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtTrackTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
          fhConeSumPtTrackTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtTrackTrigEtaPhiCent[icen]) ;
          
          fhConeSumPtClusterTrigEtaPhiCent[icen]  = new TH3F
          (Form("hConePtSumClusterTrigEtaPhi_Cen%d",icen),
           Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} clusters in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
           etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
           phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
          fhConeSumPtClusterTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtClusterTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
          fhConeSumPtClusterTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtClusterTrigEtaPhiCent[icen]) ;
        }
      }
      
      if ( fICMethod != kSumBkgSubIC )
      {
        fhConeSumPtClustervsTrackCent   = new TH3F
        ("hConePtSumClustervsTrackCent",
         Form("Track vs Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f",fConeSize),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtClustervsTrackCent->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtClustervsTrackCent->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhConeSumPtClustervsTrackCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtClustervsTrackCent) ;
        
        fhConeSumPtClusterTrackFracCent   = new TH3F
        ("hConePtSumClusterTrackFractionCent",
         Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track}, #it{R}=%2.2f, UE correction",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         fraBinsArray.GetSize() - 1, fraBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtClusterTrackFracCent->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
        fhConeSumPtClusterTrackFracCent->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
        fhConeSumPtClusterTrackFracCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtClusterTrackFracCent) ;
      }
    }
    
    if ( fICMethod >= kSumBkgSubIC )
    {
      fhConeSumPtUESubCent= new TH3F
      ("hConePtSumUESubCent",
       Form("Track and/or Cluster #Sigma #it{p}_{T} in #it{R} = %2.2f, after UE correction",fConeSize),
        ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
       sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
       cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
      fhConeSumPtUESubCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
      fhConeSumPtUESubCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      fhConeSumPtUESubCent->SetZTitle("Centrality (%)");
      outputContainer->Add(fhConeSumPtUESubCent) ;
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtUESubTrigEtaPhiCent   = new TH3F*[fNCentBins] ;
        for(Int_t icen = 0; icen < fNCentBins; icen++)
        {
          fhConeSumPtUESubTrigEtaPhiCent[icen]  = new TH3F
          (Form("hConePtSumUESubTrigEtaPhi_Cen%d",icen),
           Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
           etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
           phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
          fhConeSumPtUESubTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
          fhConeSumPtUESubTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
          fhConeSumPtUESubTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
          outputContainer->Add(fhConeSumPtUESubTrigEtaPhiCent[icen]) ;
        }
      }

      if ( fPartInCone == kNeutralAndCharged )
      {
        fhConeSumPtUESubTrackCent  = new TH3F
        ("hConePtSumUESubTrackCent",
         Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubTrackCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubTrackCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubTrackCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubTrackCent) ;   
        
        fhConeSumPtUESubClusterCent  = new TH3F
        ("hConePtSumUESubClusterCent",
         Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubClusterCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubClusterCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubClusterCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubClusterCent) ;
        
        if ( fFillEtaPhiHistograms )
        {
          fhConeSumPtUESubTrackTrigEtaPhiCent   = new TH3F*[fNCentBins] ;
          fhConeSumPtUESubClusterTrigEtaPhiCent = new TH3F*[fNCentBins] ;
          for(Int_t icen = 0; icen < fNCentBins; icen++)
          {
            fhConeSumPtUESubTrackTrigEtaPhiCent[icen]  = new TH3F
            (Form("hConePtSumUESubTrackTrigEtaPhi_Cen%d",icen),
             Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
             etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
             phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhConeSumPtUESubTrackTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
            fhConeSumPtUESubTrackTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
            fhConeSumPtUESubTrackTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtUESubTrackTrigEtaPhiCent[icen]) ;
            
            fhConeSumPtUESubClusterTrigEtaPhiCent[icen]  = new TH3F
            (Form("hConePtSumUESubClusterTrigEtaPhi_Cen%d",icen),
             Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} clusters in isolation cone for %s, cen bin %d",parTitleR.Data(),icen),
             etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
             phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhConeSumPtUESubClusterTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
            fhConeSumPtUESubClusterTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
            fhConeSumPtUESubClusterTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhConeSumPtUESubClusterTrigEtaPhiCent[icen]) ;
          }
        }
        
        if ( fICMethod != kSumBkgSubIC )
        {
          fhConeSumPtUESubClustervsTrackCent   = new TH3F
          ("hConePtSumUESubClustervsTrackCent",
           Form("Track vs Cluster #Sigma #it{p}_{T}, #it{R}=%2.2f, UE correction",fConeSize),
           sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
           sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtUESubClustervsTrackCent->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
          fhConeSumPtUESubClustervsTrackCent->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
          fhConeSumPtUESubClustervsTrackCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtUESubClustervsTrackCent) ;
          
          fhConeSumPtUESubClusterTrackFracCent   = new TH3F
          ("hConePtSumUESubClusterTrackFractionCent",
           Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track}, #it{R}=%2.2f, UE correction",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           fueBinsArray.GetSize() - 1, fueBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtUESubClusterTrackFracCent->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
          fhConeSumPtUESubClusterTrackFracCent->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
          fhConeSumPtUESubClusterTrackFracCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtUESubClusterTrackFracCent) ;
        }
      }
      
      if ( fICMethod == kSumBkgSubIC )
      {
        fhPtInPerpConeCent  = new TH3F
        ("hPtInPerpConeCent",
         Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f",fConeSize),
         ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
        ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
        cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhPtInPerpConeCent->SetYTitle("#it{p}_{T in #perp cone} (GeV/#it{c})");
        fhPtInPerpConeCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPtInPerpConeCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPtInPerpConeCent) ;

        fhPerpConeSumPtCent  = new TH3F
        ("hPerpConePtSumCent",
         Form("#Sigma #it{p}_{T} in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhPerpConeSumPtCent->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        fhPerpConeSumPtCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPtCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPerpConeSumPtCent) ;
        
        fhPerpCone1SumPtCent  = new TH3F
        ("hPerpCone1PtSumCent",
         Form("#Sigma #it{p}_{T} in 1 isolation cone at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
        sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
        cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhPerpCone1SumPtCent->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        fhPerpCone1SumPtCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPerpCone1SumPtCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPerpCone1SumPtCent) ;

        fhPerpCone1SumPtUESubCent  = new TH3F
        ("hPerpCone1PtSumUESubCent",
         Form("#Sigma #it{p}_{T} - 1 #perp cone, #it{R} =  %2.2f",fConeSize),
         ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
        sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
        cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhPerpCone1SumPtUESubCent->SetYTitle("#Sigma #it{p}_{T,UE sub} (GeV/#it{c})");
        fhPerpCone1SumPtUESubCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhPerpCone1SumPtUESubCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPerpCone1SumPtUESubCent) ;

        fhPerpConeRhoCent  = new TH2F
        ("hPerpConeRhoCent",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nrhobins,rhomin,rhomax, 100,0,100);
        fhPerpConeRhoCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRhoCent->SetYTitle("Centrality (%)");
        outputContainer->Add(fhPerpConeRhoCent) ;

        fhPerpConeRhoCutMaxCent  = new TH2F
        ("hPerpConeRhoCutMaxCent",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nrhobins,rhomin,rhomax, 100,0,100);
        fhPerpConeRhoCutMaxCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRhoCutMaxCent->SetYTitle("Centrality (%)");
        outputContainer->Add(fhPerpConeRhoCutMaxCent) ;

        fhPerpConeRhoCutLeadFactorCent  = new TH2F
        ("hPerpConeRhoCutLeadFactorCent",
         Form("#rho in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
         nrhobins,rhomin,rhomax, 100,0,100);
        fhPerpConeRhoCutLeadFactorCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
        fhPerpConeRhoCutLeadFactorCent->SetYTitle("Centrality (%)");
        outputContainer->Add(fhPerpConeRhoCutLeadFactorCent) ;

        fhConeSumPtVSPerpConeCent = new TH3F
        ("hConeSumPtVSPerpConeCent",
         Form("#perp cone #Sigma #it{p}_{T}, R=%2.2f",fConeSize),
        sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
        sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
        cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtVSPerpConeCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtVSPerpConeCent->SetYTitle("Jet #rho #it{R}^{2} #pi (GeV/#it{c})");
        fhConeSumPtVSPerpConeCent->SetZTitle("Centrality");
        outputContainer->Add(fhConeSumPtVSPerpConeCent);

        fhPerpConeSumPtTrackSubVsNoSubCent = new TH3F
        ("hPerpConeSumPtTrackSubVsNoSubCent",
         Form("Track #Sigma #it{p}_{T} in cone before and after UE bkg sub, R=%2.2f",fConeSize),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhPerpConeSumPtTrackSubVsNoSubCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhPerpConeSumPtTrackSubVsNoSubCent->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
        fhPerpConeSumPtTrackSubVsNoSubCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhPerpConeSumPtTrackSubVsNoSubCent);

        if ( fFillEtaPhiHistograms )
        {
          fhPerpConeSumPtTrigEtaPhiCent = new TH3F*[fNCentBins] ;
          for(Int_t icen = 0; icen < fNCentBins; icen++)
          {
            fhPerpConeSumPtTrigEtaPhiCent[icen]  = new TH3F
            (Form("hPerpConePtSumTrigEtaPhi_Cen%d",icen),
             Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} tracks in #perp cone for %s, cen bin %d",parTitleR.Data(),icen),
             etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
             phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
            fhPerpConeSumPtTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
            fhPerpConeSumPtTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
            fhPerpConeSumPtTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhPerpConeSumPtTrigEtaPhiCent[icen]) ;
          }
        }
      }
      
      if ( fICMethod == kSumBkgSubJetRhoIC )
      {
        fhJetRhoSumPtCent  = new TH3F
        ("hJetRhoPtSumCent",
         Form("Jet #rho #it{R}^{2} #pi, #it{R} =  %2.2f",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhJetRhoSumPtCent->SetYTitle("Jet #rho #it{R}^{2} #pi(GeV/#it{c})");
        fhJetRhoSumPtCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fhJetRhoSumPtCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhJetRhoSumPtCent) ;

        fhJetRhoCent  = new TH2F
        ("hJetRhoCent",
         Form("Jet #rho, #it{R} =  %2.2f",fConeSize),
         nrhobins,rhomin,rhomax,100,0,100);
        fhJetRhoCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
        fhJetRhoCent->SetYTitle("Centrality (%)");
        outputContainer->Add(fhJetRhoCent) ;

        fhConeSumPtVSJetRhoCent = new TH3F
        ("hConeSumPtVSJetRhoCent",
         Form("Jet #rho #it{R}^{2} #pi versus #Sigma #it{p}_{T}, R=%2.2f",fConeSize),
        sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
        sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
        cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtVSJetRhoCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtVSJetRhoCent->SetYTitle("Jet #rho #it{R}^{2} #pi (GeV/#it{c})");
        fhConeSumPtVSJetRhoCent->SetZTitle("Centrality");
        outputContainer->Add(fhConeSumPtVSJetRhoCent);

        fhJetRhoSumPtTrackSubVsNoSubCent = new TH3F
        ("hJetRhoSumPtTrackSubVsNoSubCent",
         Form("Track #Sigma #it{p}_{T} in cone before and after UE bkg sub, R=%2.2f",fConeSize),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhJetRhoSumPtTrackSubVsNoSubCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhJetRhoSumPtTrackSubVsNoSubCent->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
        fhJetRhoSumPtTrackSubVsNoSubCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhJetRhoSumPtTrackSubVsNoSubCent);

        if ( fFillEtaPhiHistograms )
        {
          fhJetRhoSumPtTrigEtaPhiCent = new TH3F*[fNCentBins] ;
          for(Int_t icen = 0; icen < fNCentBins; icen++)
          {
            fhJetRhoSumPtTrigEtaPhiCent[icen]  = new TH3F
            (Form("hJetRhoPtSumTrigEtaPhi_Cen%d",icen),
             Form("Jet #rho #it{R}^{2} #pi %s, cen bin %d",parTitleR.Data(),icen),
             etaBinsArray.GetSize()  -1, etaBinsArray.GetArray(),
             phiBinsArray.GetSize()  -1, phiBinsArray.GetArray(),
             sumBinsArray.GetSize() - 1, sumBinsArray.GetArray());
            fhJetRhoSumPtTrigEtaPhiCent[icen]->SetZTitle("Jet #rho #it{R}^{2} #pi (GeV/#it{c})");
            fhJetRhoSumPtTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
            fhJetRhoSumPtTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
            outputContainer->Add(fhJetRhoSumPtTrigEtaPhiCent[icen]) ;
          }
        }
      }

      if ( fICMethod >= kSumBkgSubEtaBandIC )
      {
        if ( fPartInCone != kOnlyCharged )
        {
          fhConeSumPtUEBandNormClusterCent  = new TH3F
          ("hConeSumPtUEBandNormClusterCent",
           Form("Clusters #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtUEBandNormClusterCent->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
          fhConeSumPtUEBandNormClusterCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtUEBandNormClusterCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtUEBandNormClusterCent) ;

          fhConeRhoUEBandClusterCent  = new TH2F
          ("hConeRhoUEBandClusterCent",
           Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandClusterCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandClusterCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandClusterCent) ;

          fhConeRhoUEBandClusterCutMaxCent  = new TH2F
          ("hConeRhoUEBandClusterCutMaxCent",
           Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandClusterCutMaxCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandClusterCutMaxCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandClusterCutMaxCent) ;

          fhConeRhoUEBandClusterCutLeadFactorCent  = new TH2F
          ("hConeRhoUEBandClusterCutLeadFactorCent",
           Form("Clusters #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandClusterCutLeadFactorCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandClusterCutLeadFactorCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandClusterCutLeadFactorCent) ;

          fhConeSumPtClusterSubVsNoSubCent = new TH3F
          ("hConeSumPtClusterSubVsNoSubCent",
           Form("Cluster #Sigma #it{p}_{T} in cone before vs after UE bkg sub, R=%2.2f",fConeSize),
          sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
          sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
          cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtClusterSubVsNoSubCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtClusterSubVsNoSubCent->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          fhConeSumPtClusterSubVsNoSubCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtClusterSubVsNoSubCent);

          fhConeSumPtVSUEClusterBandCent  = new TH3F
          ("hConeSumPtVSUEClusterBandCent",
           Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in normalized UE band for cluster, R=%2.2f",fConeSize),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtVSUEClusterBandCent->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUEClusterBandCent->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          fhConeSumPtVSUEClusterBandCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtVSUEClusterBandCent);
        }
        
        if ( fPartInCone != kOnlyNeutral )
        {
          fhConeSumPtUEBandNormTrackCent  = new TH3F
          ("hConeSumPtUEBandNormTrackCent",
           Form("Clusters #Sigma #it{p}_{T} in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtUEBandNormTrackCent->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
          fhConeSumPtUEBandNormTrackCent->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fhConeSumPtUEBandNormTrackCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtUEBandNormTrackCent) ;

          fhConeRhoUEBandTrackCent  = new TH2F
          ("hConeRhoUEBandTrackCent",
           Form("Tracks #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandTrackCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandTrackCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandTrackCent) ;

          fhConeRhoUEBandTrackCutMaxCent  = new TH2F
          ("hConeRhoUEBandTrackCutMaxCent",
           Form("Tracks #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandTrackCutMaxCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandTrackCutMaxCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandTrackCutMaxCent) ;

          fhConeRhoUEBandTrackCutLeadFactorCent  = new TH2F
          ("hConeRhoUEBandTrackCutLeadFactorCent",
           Form("Tracks #rho in normalized UE estimator, #it{R} =  %2.2f",fConeSize),
           nrhobins,rhomin,rhomax,100,0,100);
          fhConeRhoUEBandTrackCutLeadFactorCent->SetXTitle("#rho (GeV/#it{c} rad^{-2})");
          fhConeRhoUEBandTrackCutLeadFactorCent->SetYTitle("Centrality (%)");
          outputContainer->Add(fhConeRhoUEBandTrackCutLeadFactorCent) ;

          fhConeSumPtTrackSubVsNoSubCent = new TH3F
          ("hConeSumPtTrackSubVsNoSubCent",
           Form("Track #Sigma #it{p}_{T} in cone before vs after UE bkg sub, R=%2.2f",fConeSize),
          sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
          sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
          cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtTrackSubVsNoSubCent->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtTrackSubVsNoSubCent->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
          fhConeSumPtTrackSubVsNoSubCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtTrackSubVsNoSubCent);

          fhConeSumPtVSUETrackBandCent  = new TH3F
          ("hConeSumPtVSUETrackBandCent",
           Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in normalized UE band for tracks, R=%2.2f",fConeSize),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtVSUETrackBandCent->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
          fhConeSumPtVSUETrackBandCent->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
          fhConeSumPtVSUETrackBandCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtVSUETrackBandCent);
        }
      }
      
      if ( fFillEtaPhiHistograms )
      {
        if ( fPartInCone != kOnlyNeutral )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
          {
            fhConeSumPtEtaBandUETrackTrigEtaPhiCent = new TH3F*[fNCentBins] ;
            for(Int_t icen = 0; icen < fNCentBins; icen++)
            {
              fhConeSumPtEtaBandUETrackTrigEtaPhiCent[icen]  = new TH3F
              (Form("hConePtSumEtaBandUETrackTrigEtaPhi_Cen%d",icen),
               Form("Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE #eta Band, cen %d",icen),
               etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
               phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
              fhConeSumPtEtaBandUETrackTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}");
              fhConeSumPtEtaBandUETrackTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
              fhConeSumPtEtaBandUETrackTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
              outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhiCent[icen]) ;
            }
          }

          if ( fICMethod == kSumBkgSubPhiBandIC )
          {
            fhConeSumPtPhiBandUETrackTrigEtaPhiCent = new TH3F*[fNCentBins] ;
            for(Int_t icen = 0; icen < fNCentBins; icen++)
            {
              fhConeSumPtPhiBandUETrackTrigEtaPhiCent[icen]  = new TH3F
              (Form("hConePtSumPhiBandUETrackTrigEtaPhi_Cen%d",icen),
               Form("Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE #varphi Band, cen %d",icen),
               etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
               phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
              fhConeSumPtPhiBandUETrackTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}");
              fhConeSumPtPhiBandUETrackTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
              fhConeSumPtPhiBandUETrackTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
              outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhiCent[icen]) ;
            }
          }

          if ( fICMethod == kSumBkgSubPerpBandIC )
          {
            fhConeSumPtPerpBandUETrackTrigEtaPhiCent = new TH3F*[fNCentBins] ;
            for(Int_t icen = 0; icen < fNCentBins; icen++)
            {
              fhConeSumPtPerpBandUETrackTrigEtaPhiCent[icen]  = new TH3F
              (Form("hConePtSumPerpBandUETrackTrigEtaPhi_Cen%d",icen),
               Form("Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE #perp #eta Band, cen %d",icen),
               etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
               phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
              fhConeSumPtPerpBandUETrackTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}");
              fhConeSumPtPerpBandUETrackTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
              fhConeSumPtPerpBandUETrackTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
              outputContainer->Add(fhConeSumPtPerpBandUETrackTrigEtaPhiCent[icen]) ;
            }
          }
        } // tracks

        if ( fPartInCone != kOnlyCharged )
        {
          if ( fICMethod == kSumBkgSubEtaBandIC )
          {
            fhConeSumPtEtaBandUEClusterTrigEtaPhiCent = new TH3F*[fNCentBins] ;
            for(Int_t icen = 0; icen < fNCentBins; icen++)
            {
              fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[icen]  = new TH3F
              (Form("hConePtSumEtaBandUEClusterTrigEtaPhi_Cen%d",icen),
               Form("Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE #eta Band, cen %d",icen),
               etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
               phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
              fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}");
              fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
              fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
              outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhiCent[icen]) ;
            }
          }

          if ( fICMethod == kSumBkgSubPhiBandIC )
          {
            fhConeSumPtPhiBandUEClusterTrigEtaPhiCent = new TH3F*[fNCentBins] ;
            for(Int_t icen = 0; icen < fNCentBins; icen++)
            {
              fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[icen]  = new TH3F
              (Form("hConePtSumPhiBandUEClusterTrigEtaPhi_Cen%d",icen),
               Form("Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE #varphi Band, cen %d",icen),
               etaBinsArray.GetSize() - 1,  etaBinsArray.GetArray(),
               phiBinsArray.GetSize() - 1,  phiBinsArray.GetArray(),
               sumBinsArray.GetSize() - 1,  sumBinsArray.GetArray());
              fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[icen]->SetZTitle("#Sigma #it{p}_{T}");
              fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[icen]->SetXTitle("#eta_{trigger}");
              fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[icen]->SetYTitle("#varphi_{trigger} (rad)");
              outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhiCent[icen]) ;
            }
          }
        } // clusters
      }

      if ( fPartInCone != kOnlyCharged && fICMethod >= kSumBkgSubIC )
      {
        fhConeSumPtUESubClusterCutMaxCent  = new TH3F
        ("hConePtSumUESubClusterCutMaxCent",
         Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on leading UE pT",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubClusterCutMaxCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubClusterCutMaxCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubClusterCutMaxCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubClusterCutMaxCent) ;

        fhConeSumPtUESubClusterCutLeadFactorCent  = new TH3F
        ("hConePtSumUESubClusterCutLeadFactorCent",
         Form("Cluster #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on leading UE pT",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubClusterCutLeadFactorCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubClusterCutLeadFactorCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubClusterCutLeadFactorCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubClusterCutLeadFactorCent) ;

        if ( fFillLeadingVsUESubHisto )
        {
          fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent          = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent        = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent      = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent = new TH3F*[fNCentBins] ;

          for (Int_t icent = 0; icent < fNCentBins; icent++)
          {
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadClusterInConePt_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[icent]->SetYTitle("#it{p}^{lead}_{cluster} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[icent]->SetYTitle("#it{p}^{lead}_{cluster} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[icent]->SetYTitle("#it{p}^{lead}_{cluster} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt_Cent%d",icent),Form("%s",parTitleR.Data()),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             lfrBinsArray.GetSize() - 1, lfrBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[icent]->SetYTitle("#it{p}^{lead}_{track} Cone / UE ");
            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[icent]) ;
          }
        }

        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhEtaBandClusterPtCent  = new TH3F
          ("hEtaBandClusterPtCent",
           Form("Clusters in #eta band out of cone #it{R} =  %2.2f",fConeSize),
           ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhEtaBandClusterPtCent->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          fhEtaBandClusterPtCent->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
          fhEtaBandClusterPtCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhEtaBandClusterPtCent) ;
        }
        
        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhPhiBandClusterPtCent  = new TH3F
          ("hPhiBandClusterPtCent",
           Form("Clusters in #varphi band out of cone #it{R} =  %2.2f",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhPhiBandClusterPtCent->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          fhPhiBandClusterPtCent->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
          fhPhiBandClusterPtCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhPhiBandClusterPtCent) ;
        }
      }
      
      if ( fPartInCone != kOnlyCharged && fICMethod >= kSumBkgSubEtaBandIC )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUEClusterCent  = new TH3F
          ("hConePtSumEtaBandUEClusterCent",
           "#Sigma cluster #it{p}_{T} in UE Eta Band",
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtEtaBandUEClusterCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUEClusterCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          fhConeSumPtEtaBandUEClusterCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtEtaBandUEClusterCent) ;
        }
        
        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUEClusterCent  = new TH3F
          ("hConePtSumPhiBandUEClusterCent",
           "#Sigma cluster #it{p}_{T} UE Phi Band",
           ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtPhiBandUEClusterCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiBandUEClusterCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          fhConeSumPtPhiBandUEClusterCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtPhiBandUEClusterCent) ;
        }
      }
      
      if ( fPartInCone != kOnlyNeutral && fICMethod >= kSumBkgSubIC )
      {
        fhConeSumPtUESubTrackCutMaxCent  = new TH3F
        ("hConePtSumUESubTrackCutMaxCent",
         Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on UE max pT",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubTrackCutMaxCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubTrackCutMaxCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubTrackCutMaxCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubTrackCutMaxCent) ;

        fhConeSumPtUESubTrackCutLeadFactorCent  = new TH3F
        ("hConePtSumUESubTrackCutLeadFactorCent",
         Form("Track #Sigma #it{p}_{T},#it{R}=%2.2f, UE correction, cut on leading UE pT",fConeSize),
          ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUESubTrackCutLeadFactorCent->SetYTitle("#Sigma #it{p}_{T} - UE (GeV/#it{c})");
        fhConeSumPtUESubTrackCutLeadFactorCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
        fhConeSumPtUESubTrackCutLeadFactorCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUESubTrackCutLeadFactorCent) ;

        if ( fFillLeadingVsUESubHisto )
        {
          fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent          = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent        = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent      = new TH3F*[fNCentBins] ;
          fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent = new TH3F*[fNCentBins] ;

          for (Int_t icent = 0; icent < fNCentBins; icent++)
          {
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadTrackInConePt_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[icent]->SetYTitle("#it{p}^{lead}_{track} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[icent]->SetYTitle("#it{p}^{lead UE}_{track} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsTrackLeadInConePtLowUELead_Cent%d",icent),Form("%s, cent %d",parTitleR.Data(),icent),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[icent]->SetYTitle("#it{p}^{lead}_{track} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[icent]) ;

            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[icent] = new TH3F
            (Form("hTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt_Cent%d",icent),Form("%s",parTitleR.Data()),
              ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
             lfrBinsArray.GetSize() - 1, lfrBinsArray.GetArray(),
             sueBinsArray.GetSize() - 1, sueBinsArray.GetArray());
            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[icent]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[icent]->SetYTitle("#it{p}^{lead}_{track} Cone / UE ");
            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[icent]->SetZTitle("#it{p}_{T}^{iso} (GeV/#it{c})");
            outputContainer->Add(fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[icent]) ;
          }
        }

        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhEtaBandTrackPtCent  = new TH3F
          ("hEtaBandTrackPtCent",
           Form("Tracks in #eta band out of cone #it{R} =  %2.2f",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhEtaBandTrackPtCent->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          fhEtaBandTrackPtCent->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
          fhEtaBandTrackPtCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhEtaBandTrackPtCent) ;
        }
        
        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhPhiBandTrackPtCent  = new TH3F
          ("hPhiBandTrackPtCent",
           Form("Tracks in #varphi band out of cone #it{R} = %2.2f and half TPC, #pm #pi",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhPhiBandTrackPtCent->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          fhPhiBandTrackPtCent->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
          fhPhiBandTrackPtCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhPhiBandTrackPtCent) ;
        }

        if ( fICMethod == kSumBkgSubPerpBandIC )
        {
          fhPerpBandTrackPtCent  = new TH3F
          ("hPerpBandTrackPtCent",
           Form("Tracks in #perp #eta band out of cone #it{R} = %2.2f and half TPC, #pm #pi",fConeSize),
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           ptCBinsArray.GetSize() - 1, ptCBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhPerpBandTrackPtCent->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
          fhPerpBandTrackPtCent->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
          fhPerpBandTrackPtCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhPerpBandTrackPtCent) ;
        }
      }
      
      if ( fPartInCone != kOnlyNeutral && fICMethod >= kSumBkgSubEtaBandIC )
      {
        if ( fICMethod == kSumBkgSubEtaBandIC )
        {
          fhConeSumPtEtaBandUETrackCent  = new TH3F
          ("hConePtSumEtaBandUETrackCent",
           "#Sigma track #it{p}_{T} in UE #eta Band",
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtEtaBandUETrackCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtEtaBandUETrackCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          fhConeSumPtEtaBandUETrackCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtEtaBandUETrackCent) ;
        }
        
        if ( fICMethod == kSumBkgSubPhiBandIC )
        {
          fhConeSumPtPhiBandUETrackCent  = new TH3F
          ("hConePtSumPhiBandUETrackCent",
           "#Sigma track #it{p}_{T} UE #varphi Band",
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtPhiBandUETrackCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPhiBandUETrackCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          fhConeSumPtPhiBandUETrackCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtPhiBandUETrackCent) ;
        }

        if ( fICMethod == kSumBkgSubPerpBandIC )
        {
          fhConeSumPtPerpBandUETrackCent  = new TH3F
          ("hConePtSumPerpBandUETrackCent",
           "#Sigma track #it{p}_{T} UE #perp #eta Band",
            ptBinsArray.GetSize() - 1,  ptBinsArray.GetArray(),
           sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
           cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
          fhConeSumPtPerpBandUETrackCent->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
          fhConeSumPtPerpBandUETrackCent->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
          fhConeSumPtPerpBandUETrackCent->SetZTitle("Centrality (%)");
          outputContainer->Add(fhConeSumPtPerpBandUETrackCent) ;
        }
      }

      if ( fPartInCone == kNeutralAndCharged && fICMethod >= kSumBkgSubEtaBandIC )
      {
        fhConeSumPtUEBandSubClustervsTrackCent   = new TH3F
        ("hConePtSumUEBandSubClustervsTrackCent",
         Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         sueBinsArray.GetSize() - 1, sueBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhConeSumPtUEBandSubClustervsTrackCent->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtUEBandSubClustervsTrackCent->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhConeSumPtUEBandSubClustervsTrackCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhConeSumPtUEBandSubClustervsTrackCent) ;

        fhBandClustervsTrackCent   = new TH3F
        ("hBandClustervsTrackCent",
         Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhBandClustervsTrackCent->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhBandClustervsTrackCent->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhBandClustervsTrackCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhBandClustervsTrackCent) ;

        fhBandNormClustervsTrackCent   = new TH3F
        ("hBandNormClustervsTrackCent",
         Form("Track vs Cluster Normalized #Sigma #it{p}_{T} in #eta or #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         sumBinsArray.GetSize() - 1, sumBinsArray.GetArray(),
         cenBinsArray.GetSize()  -1, cenBinsArray.GetArray());
        fhBandNormClustervsTrackCent->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhBandNormClustervsTrackCent->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        fhBandNormClustervsTrackCent->SetZTitle("Centrality (%)");
        outputContainer->Add(fhBandNormClustervsTrackCent) ;
      } // UE, neutral + charged
    }
      
  }
  
  return outputContainer;
}
  
//____________________________________________
/// Put data member values in string to keep
/// in output container.
//____________________________________________
TString AliIsolationCut::GetICParametersList()
{
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliIsolationCut ---:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fConeSize=%1.2f; Gap %1.2f",fConeSize,fConeSizeBandGap) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtThreshold>%2.2f;<%2.2f;",fPtThreshold,fPtThresholdMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fSumPtThreshold<%2.2f, gap = %2.2f;<%2.2f;",fSumPtThreshold, fSumPtThresholdGap, fSumPtThresholdMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtFraction=%2.2f;",fPtFraction) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fICMethod=%d;",fICMethod) ;
  parList+=onePar ;
  if ( fICMethod == kSumBkgSubJetRhoIC )
  {
    snprintf(onePar,buffersize,"fJetRhoTaskName=%s",fJetRhoTaskName.Data());
  }
  snprintf(onePar,buffersize,"fPartInCone=%d;",fPartInCone) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fFracIsThresh=%i;",fFracIsThresh) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," fIsTMClusterInConeRejected=%i;", fIsTMClusterInConeRejected) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fDistMinToTrigger=%1.2f;",fDistMinToTrigger) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fFillHistograms=%d,fFillEtaPhiHistograms=%d; pt>%2.1f",fFillHistograms,fFillEtaPhiHistograms,fEtaPhiHistogramsMinPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMakeConeExcessCorr=%d;",fMakeConeExcessCorr) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNeutralOverChargedRatio={%1.2e,%1.2e,%1.2e,%1.2e};",
           fNeutralOverChargedRatio[0],fNeutralOverChargedRatio[1],fNeutralOverChargedRatio[2],fNeutralOverChargedRatio[3]) ;
  parList+=onePar ;
  
  return parList;
}

//____________________________________
// Initialize the parameters of the analysis.
//____________________________________
void AliIsolationCut::InitParameters()
{
  fFillHistograms       = kFALSE; // True in GetCreateOutputObjects();
  fFillEtaPhiHistograms = kFALSE; 
  fEtaPhiHistogramsMinPt= 10;
  fFillHighMultHistograms = kFALSE;
  fNCentBins            = 10 ;
  fConeSize             = 0.4 ;
  fConeSizeBandGap      = 0.0 ;
  fMakeConeExcessCorr   = kFALSE;
  fPtThreshold          = 0.5  ;
  fPtThresholdMax       = 10000.;
  fSumPtThreshold       = 2.0 ;
  fSumPtThresholdGap    = 0.5 ;
  fSumPtThresholdMax    = 10000. ;
  fPtFraction           = 0.1 ;
  fPartInCone           = kNeutralAndCharged;
  fICMethod             = kSumPtIC; // 0 pt threshol method, 1 cone pt sum method
  fFracIsThresh         = 1;
  fDistMinToTrigger     = -1.; // no effect
  fLeadingPtUEFactor    = 1.8;
  fMaxPtUE              = 3.0;
  fJetRhoTaskName       = "Rho";

  // Ratio charged to neutral
  // Based on pPb analysis, Erwann Masson Thesis 
  // Use eta band for charged and neutrals for estimation.
  fNeutralOverChargedRatio[0] = 0.363; 
  fNeutralOverChargedRatio[1] = 0.0;
  fNeutralOverChargedRatio[2] = 0.0;
  fNeutralOverChargedRatio[3] = 0.0;
  
  // Pb-Pb vs centrality (eta band, random trigger cone, LHC18qr)
//  fNeutralOverChargedRatio[0] = 1.49e-1; // 1.29e-1 phi band, 1.26e-1 trigger cone
//  fNeutralOverChargedRatio[1] =-2.54e-3; //-2.47e-3 phi band,-2.27E-3 trigger cone
//  fNeutralOverChargedRatio[2] = 0.0;
//  fNeutralOverChargedRatio[3] = 7.32e-7; // 7.38e-7 phi band, 6.91e-7 trigger cone
}

//________________________________________________________________________________
/// Declare a candidate particle isolated depending on the
/// cluster or track particle multiplicity and/or momentum.
///
/// \param pCandidate: Kinematics and + of candidate particle for isolation.
/// \param reader: pointer to AliCaloTrackReader. Needed to access event info.
/// \param bFillAOD: Indicate if particles in cone must be added to AOD particle object.
/// \param useRefs: Get the list of tracks or clusters in cone from references
/// \param aodArrayRefName: Name of array where list of tracks/clusters in cone is stored.
/// \param bgTrk: List of tracks from mixed event background pool AliAnaParticleHadronCorrelation.
/// \param bgCls: List of clusters from mixed event background pool AliAnaParticleHadronCorrelation.
/// \param calorimeter: Which input trigger calorimeter used
/// \param pid: pointer to AliCaloPID. Needed to reject matched clusters in isolation cone.
/// \param nPart: number of tracks/clusters above threshold in cone, output.
/// \param nfrac: 1 if fraction pT cluster-track / pT trigger in cone avobe threshold, output.
/// \param coneptsum: total momentum energy in cone (track+cluster), output.
/// \param ptLead: momentum of leading cluster or track in cone, output.
/// \param isolated: final bool with decission on isolation of candidate particle.
/// \param histoWeight: Histograms weight (event, pt depedent).
/// \param centrality:Centrality percentile
/// \param cenBin: Assigned centrality bin with defined range elsewhere
///
//________________________________________________________________________________
void  AliIsolationCut::MakeIsolationCut
(
 AliCaloTrackParticleCorrelation  *pCandidate, AliCaloTrackReader * reader,
 Bool_t    bFillAOD,    Bool_t    useRefs, TString aodArrayRefName,
 TObjArray * bgTrk, TObjArray * bgCls,
 Int_t     calorimeter, AliCaloPID * pid,
 Int_t   & nPart      , Int_t   & nfrac,
 Float_t & coneptsum  , Float_t & ptLead,
 Bool_t  & isolated   , Double_t histoWeight, 
 Float_t   centrality , Int_t cenBin
)
{
  if ( fPartInCone == kOnlyNeutral && fICMethod == kSumBkgSubIC )
  {
    AliFatal("Only neutral in cone and perpendicular cones combination not supported \n");
    return;
  }
  
  Float_t ptC   = pCandidate->Pt() ;
  Float_t phiC  = pCandidate->Phi() ;
  if ( phiC < 0 ) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  Float_t coneptsumCluster    = 0;
  Float_t coneptsumTrack      = 0;
  Float_t coneptLeadCluster   = 0;
  Float_t etaBandptLeadCluster= 0;
  Float_t phiBandptLeadCluster= 0;
  Float_t coneptLeadTrack     = 0;
  Float_t etaBandptLeadTrack  = 0;
  Float_t phiBandptLeadTrack  = 0;
  Float_t perpBandptLeadTrack = 0;
  Float_t perpConeptLeadTrack = 0;
  Float_t ueptLeadTrack       = 0;
  Float_t ueptLeadCluster     = 0;

  Float_t etaBandPtSumTrack   = 0;
  Float_t phiBandPtSumTrack   = 0;
  Float_t perpBandPtSumTrack  = 0;
  Float_t perpPtSumTrack      = 0;
  Float_t perp1PtSumTrack     = 0;

  Float_t etaBandPtSumTrackCutMax           = 0. ;
  Float_t phiBandPtSumTrackCutMax           = 0. ;
  Float_t perpBandPtSumTrackCutMax          = 0. ;
  Float_t perpPtSumTrackCutMax              = 0. ;

  Float_t etaBandPtSumTrackCutLeadFactor    = 0. ;
  Float_t phiBandPtSumTrackCutLeadFactor    = 0. ;
  Float_t perpBandPtSumTrackCutLeadFactor   = 0. ;
  Float_t perpPtSumTrackCutLeadFactor       = 0. ;

  Float_t etaBandPtSumCluster               = 0. ;
  Float_t phiBandPtSumCluster               = 0. ;
  Float_t etaBandPtSumClusterCutMax         = 0. ;
  Float_t phiBandPtSumClusterCutMax         = 0. ;
  Float_t etaBandPtSumClusterCutLeadFactor  = 0. ;
  Float_t phiBandPtSumClusterCutLeadFactor  = 0. ;
  Float_t rhoPtSumTrack                     = 0. ;

  Float_t  coneptsumTrackSub                = 0. ;
  Float_t  coneptsumClusterSub              = 0. ;
  Float_t  coneptsumTrackSubCutMax          = 0. ;
  Float_t  coneptsumClusterSubCutMax        = 0. ;
  Float_t  coneptsumTrackSubCutLeadFactor   = 0. ;
  Float_t  coneptsumClusterSubCutLeadFactor = 0. ;

  Float_t coneArea = fConeSize*fConeSize*TMath::Pi();
  // If there is a central hole in the trigger particle cone, recalculate cone area
  if ( fDistMinToTrigger > 0 ) coneArea -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  nPart     = 0 ;
  nfrac     = 0 ;
  isolated  = kFALSE;
  
  AliDebug(1,Form("Candidate pT %2.2f, eta %2.2f, phi %2.2f, cone %1.2f, thres %2.2f, Fill AOD? %d",
                  pCandidate->Pt(), pCandidate->Eta(), pCandidate->Phi()*TMath::RadToDeg(), 
                  fConeSize, fPtThreshold, bFillAOD));
  
  // ------------------------------------------
  // Get charged tracks and clusters in cone.
  // Add the pt or get the leading one
  // ------------------------------------------

  CalculateTrackSignalInCone   (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgTrk,
                                nPart              , nfrac,
                                coneptsumTrack     , coneptLeadTrack    ,
                                etaBandPtSumTrack  , etaBandptLeadTrack , etaBandPtSumTrackCutMax , etaBandPtSumTrackCutLeadFactor ,
                                phiBandPtSumTrack  , phiBandptLeadTrack , phiBandPtSumTrackCutMax , phiBandPtSumTrackCutLeadFactor ,
                                perpPtSumTrack     , perpConeptLeadTrack, perpPtSumTrackCutMax    , perpPtSumTrackCutLeadFactor    ,
                                perpBandPtSumTrack , perpBandptLeadTrack, perpBandPtSumTrackCutMax, perpBandPtSumTrackCutLeadFactor,
                                perp1PtSumTrack    , histoWeight,
                                centrality         , cenBin);

  CalculateCaloSignalInCone    (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgCls,
                                calorimeter        , pid, 
                                nPart              , nfrac,
                                coneptsumCluster   , coneptLeadCluster   ,
                                etaBandPtSumCluster, etaBandptLeadCluster, etaBandPtSumClusterCutMax, etaBandPtSumClusterCutLeadFactor,
                                phiBandPtSumCluster, phiBandptLeadCluster, phiBandPtSumClusterCutMax, phiBandPtSumClusterCutLeadFactor,
                                histoWeight        , 
                                centrality         , cenBin);

  // Add leading found information to candidate object
  pCandidate->SetNeutralLeadPtInCone(coneptLeadCluster);
  pCandidate->SetChargedLeadPtInCone(coneptLeadTrack);
  //printf("Lead pT track %2.2f, cluster %2.2f \n",coneptLeadTrack,coneptLeadCluster);

  // Get detectors acceptance
  // Do it once, needed for Band UE estimation and excess area determination
  //if ( fMakeConeExcessCorr || fICMethod >= kSumBkgSubEtaBandIC )  // comment not needed in all cases, but it does not take much, in case UE study done in AliAnaParticleIsolation
  GetDetectorAngleLimits(reader,calorimeter); 

  // Calculate how much of the cone got out the detectors acceptance
  //
  Float_t excessAreaTrkEta = 1;
  Float_t excessAreaClsEta = 1;
  Float_t excessAreaClsPhi = 1;
  Float_t excessTrkEta     = 0;
  Float_t excessClsEta     = 0;
  Float_t excessClsPhi     = 0;
  
  if ( fMakeConeExcessCorr )
  {
    CalculateExcessAreaFractionForChargedAndNeutral(etaC, phiC, 
                                                    pCandidate->GetDetectorTag(),
                                                    excessTrkEta, excessAreaTrkEta, 
                                                    excessClsEta, excessAreaClsEta, 
                                                    excessClsPhi, excessAreaClsPhi);
    
    // Store excess areas
    pCandidate->SetNeutralConeExcessAreaEta(excessAreaClsEta);
    pCandidate->SetNeutralConeExcessAreaPhi(excessAreaClsPhi);
    pCandidate->SetChargedConeExcessAreaEta(excessAreaTrkEta);
    pCandidate->SetChargedConeExcessAreaPhi(1);
    
    if ( excessAreaTrkEta > 3 || excessAreaClsPhi > 3 || excessAreaClsEta > 3 )
    {
      AliWarning(Form("Candidate excess area too large: pT %2.2f, eta %2.2f, phi %2.2f, cone %2.2f;"
                      " excess: trk eta %2.2f area %2.2f, cls eta %2.2f area %2.2f, cls phi %2.2f area %2.2f\n",
                      ptC, etaC, phiC, fConeSize,
                      excessTrkEta, excessAreaTrkEta,
                      excessClsEta, excessAreaClsEta,
                      excessClsPhi, excessAreaClsPhi));
    }

    if ( fFillHistograms && fFillFractionExcessHistograms )
    {
      if ( fPartInCone != kOnlyNeutral )
      {
        fhFractionTrackOutConeEta            ->Fill(ptC ,       excessAreaTrkEta, histoWeight);
        fhFractionTrackOutConeEtaTrigEtaPhi  ->Fill(etaC, phiC, excessAreaTrkEta, histoWeight);
      }

      if ( fPartInCone != kOnlyCharged )
      {
        fhFractionClusterOutConeEta          ->Fill(ptC ,       excessAreaClsEta, histoWeight);
        fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaC, phiC, excessAreaClsEta, histoWeight);
        fhFractionClusterOutConePhi          ->Fill(ptC ,       excessAreaClsPhi, histoWeight);
        fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaC, phiC, excessAreaClsPhi, histoWeight);
        fhFractionClusterOutConeEtaPhi       ->Fill(ptC ,       excessAreaClsPhi*excessAreaClsEta, histoWeight);
        fhFractionClusterOutConeEtaPhiTrigEtaPhi->Fill(etaC, phiC, excessAreaClsPhi*excessAreaClsEta, histoWeight);
      }
    }
  }
  
  // Calculate sum of cone energy, applying excess cone correction factor.
  //
  coneptsum = coneptsumCluster * excessAreaClsEta * excessAreaClsPhi + 
              coneptsumTrack   * excessAreaTrkEta;
  // coneptsum = coneptsumCluster + coneptsumTrack ;  
  
  //-------------------------------------------------------------------
  // Check isolation, depending on selected isolation criteria requested
  //-------------------------------------------------------------------
  Double_t coneptsumUESub        = coneptsumCluster+coneptsumTrack;
  Double_t coneptsumUESubCluster = coneptsumCluster;
  Double_t coneptsumUESubTrack   = coneptsumTrack;
  
  // *Now*, just check the leading particle in the cone if the threshold is passed
  if ( ptLead > fPtThreshold && ptLead < fPtThresholdMax)  nPart = 1;
  
  //if fPtFraction*ptC<fPtThreshold then consider the fPtThreshold directly
  if ( fFracIsThresh )
  {
    if ( fPtFraction*ptC < fPtThreshold )
    {
      if ( ptLead > fPtThreshold )    nfrac = 1 ;
    }
    else
    {
      if ( ptLead > fPtFraction*ptC ) nfrac = 1;
    }
  }
  else
  {
    if ( ptLead > fPtFraction*ptC )   nfrac = 1;
  }
  
  // Here start each of the selection methods
  //
  if      ( fICMethod == kPtThresIC )
  {
    if ( nPart == 0 ) isolated = kTRUE ;
    
    AliDebug(1,Form("pT Cand %2.2f, pT Lead %2.2f, %2.2f<pT Lead< %2.2f, isolated %d",
                    ptC,ptLead,fPtThreshold,fPtThresholdMax,isolated));
  }
  else if ( fICMethod == kSumPtIC )
  {
    if ( coneptsum > fSumPtThreshold &&
         coneptsum < fSumPtThresholdMax )
      isolated  =  kFALSE ;
    else
      isolated  =  kTRUE  ;
    
    AliDebug(1,Form("pT Cand %2.2f, SumPt %2.2f, %2.2f<Sum pT< %2.2f, isolated %d",
                    ptC,ptLead,fSumPtThreshold,fSumPtThresholdMax,isolated));
  }
  else if ( fICMethod == kPtFracIC )
  {
    if ( nfrac == 0 ) isolated = kTRUE ;
  }
  else if ( fICMethod == kSumPtFracIC )
  {
    //when the fPtFraction*ptC < fSumPtThreshold then consider the later case
    // printf("photon analysis IsDataMC() ?%i\n",IsDataMC());
    if ( fFracIsThresh )
    {
      if( fPtFraction*ptC < fSumPtThreshold  && coneptsum < fSumPtThreshold ) isolated  =  kTRUE ;
      if( fPtFraction*ptC > fSumPtThreshold  && coneptsum < fPtFraction*ptC ) isolated  =  kTRUE ;
    }
    else
    {
      if ( coneptsum < fPtFraction*ptC ) isolated  =  kTRUE ;
    }
  }
  else if ( fICMethod == kSumDensityIC )
  {
    // Get good cell density (number of active cells over all cells in cone)
    // and correct energy in cone.
    // Old method to study feasibility of analysis with cells as input and not clusters
    Float_t cellDensity = GetCellDensity(pCandidate,reader);
    
    if ( coneptsum < fSumPtThreshold*cellDensity )
      isolated = kTRUE;
  }
  else if ( fICMethod >= kSumBkgSubIC )
  {
    //printf("Do perp isolation\n");
    // First get the background to be subtracted to the cone content
    // Different ways: perpendicular cones, eta band, phi band
    //Double_t coneptsumBkg    = 0.;
    Double_t coneptsumBkgTrk = 0.;
    Double_t coneptsumBkgCls = 0.;
    
    Double_t coneptsumBkgTrkRaw = 0.;
    Double_t coneptsumBkgClsRaw = 0.;

    Double_t coneptsumBkgTrkCutMax = 0.;
    Double_t coneptsumBkgClsCutMax = 0.;
    Double_t coneptsumBkgTrkCutLeadFactor = 0.;
    Double_t coneptsumBkgClsCutLeadFactor = 0.;

    if ( fICMethod == kSumBkgSubIC )
    {
      ueptLeadTrack = perpConeptLeadTrack;

      // If there is a central hole in the trigger particle cone, scaledown the perp cone energy
      if ( fDistMinToTrigger > 0 )
      {
        perpPtSumTrack              *= (fConeSize*fConeSize - fDistMinToTrigger*fDistMinToTrigger) / (fConeSize*fConeSize);
        perpPtSumTrackCutMax        *= (fConeSize*fConeSize - fDistMinToTrigger*fDistMinToTrigger) / (fConeSize*fConeSize);
        perpPtSumTrackCutLeadFactor *= (fConeSize*fConeSize - fDistMinToTrigger*fDistMinToTrigger) / (fConeSize*fConeSize);
      }

      coneptsumBkgTrk                = perpPtSumTrack;
      coneptsumBkgTrkCutMax          = perpPtSumTrackCutMax;
      coneptsumBkgTrkCutLeadFactor   = perpPtSumTrackCutLeadFactor;
      coneptsumTrackSub              = coneptsumTrack - perpPtSumTrack;
      coneptsumTrackSubCutMax        = coneptsumTrack - perpPtSumTrackCutMax;
      coneptsumTrackSubCutLeadFactor = coneptsumTrack - perpPtSumTrackCutLeadFactor;

      if ( fPartInCone == kNeutralAndCharged || fPartInCone == kOnlyNeutral )
        coneptsumBkgCls = perpPtSumTrack*GetNeutralOverChargedRatio(centrality);
      //printf("centrality %f, neutral/charged %f\n",centrality,GetNeutralOverChargedRatio(centrality));
      
      // Add to candidate object
      if      ( fUseMaxPtUE )
        pCandidate->SetChargedPtSumInPerpCone(perpPtSumTrackCutMax*excessAreaTrkEta);
      else if ( fUseLeadingPtUEFactor )
        pCandidate->SetChargedPtSumInPerpCone(perpPtSumTrackCutLeadFactor*excessAreaTrkEta);
      else
        pCandidate->SetChargedPtSumInPerpCone(perpPtSumTrack*excessAreaTrkEta);
      
      if ( fFillHistograms )
      {
        Float_t rho           =  perpPtSumTrack / coneArea;
        Float_t rhoMax        =  perpPtSumTrackCutMax / coneArea;
        Float_t rhoLeadFactor =  perpPtSumTrackCutLeadFactor / coneArea;
        if ( fFillHighMultHistograms )
        {
          fhPerpConeSumPtCent->Fill(ptC, perpPtSumTrack * excessAreaTrkEta, centrality, histoWeight);

          fhPerpCone1SumPtCent->Fill(ptC, perp1PtSumTrack * excessAreaTrkEta, centrality, histoWeight);
          fhPerpCone1SumPtUESubCent->Fill(ptC, (coneptsumTrack-perp1PtSumTrack) * excessAreaTrkEta, centrality, histoWeight);

          if ( ptC > fEtaPhiHistogramsMinPt )
          {
            fhPerpConeRhoCent             ->Fill(rho           * excessAreaTrkEta, centrality, histoWeight);
            fhPerpConeRhoCutMaxCent       ->Fill(rhoMax        * excessAreaTrkEta, centrality, histoWeight);
            fhPerpConeRhoCutLeadFactorCent->Fill(rhoLeadFactor * excessAreaTrkEta, centrality, histoWeight);

            fhConeSumPtVSPerpConeCent->Fill(coneptsumTrack * excessAreaTrkEta,
                                            perpPtSumTrack * excessAreaTrkEta,
                                            centrality, histoWeight);

            fhPerpConeSumPtTrackSubVsNoSubCent->Fill(coneptsumTrack * excessAreaTrkEta,
                                                     (coneptsumTrack-perpPtSumTrack) * excessAreaTrkEta,
                                                     centrality, histoWeight);
            if ( fFillEtaPhiHistograms  && cenBin < fNCentBins && cenBin >= 0 )
              fhPerpConeSumPtTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, perpPtSumTrack * excessAreaTrkEta, histoWeight);
          }
        }
        else
        {
          fhPerpConeRho             ->Fill(ptC, rho           * excessAreaTrkEta, histoWeight);
          fhPerpConeRhoCutMax       ->Fill(ptC, rhoMax        * excessAreaTrkEta, histoWeight);
          fhPerpConeRhoCutLeadFactor->Fill(ptC, rhoLeadFactor * excessAreaTrkEta, histoWeight);

          fhPerpConeSumPt->Fill(ptC, perpPtSumTrack * excessAreaTrkEta, histoWeight);

          fhPerpCone1SumPt->Fill(ptC, perp1PtSumTrack * excessAreaTrkEta, histoWeight);
          fhPerpCone1SumPtUESub->Fill(ptC, (coneptsumTrack-perp1PtSumTrack) * excessAreaTrkEta, histoWeight);

          if ( ptC > fEtaPhiHistogramsMinPt )
          {
            fhConeSumPtVSPerpCone->Fill(coneptsumTrack * excessAreaTrkEta,
                                        perpPtSumTrack * excessAreaTrkEta,
                                        histoWeight);
            fhPerpConeSumPtTrackSubVsNoSub->Fill(coneptsumTrack * excessAreaTrkEta,
                                                 (coneptsumTrack-perpPtSumTrack) * excessAreaTrkEta,
                                                 histoWeight);
            if ( fFillEtaPhiHistograms )
              fhPerpConeSumPtTrigEtaPhi->Fill(etaC, phiC, perpPtSumTrack * excessAreaTrkEta, histoWeight);
          }
        }
      } // fill perp cone histograms
    } // UE subtraction by perpendicular cones
    else if ( fICMethod == kSumBkgSubJetRhoIC )
    {
      // Rely on Jet group methods to get the UE density.
      // https://alice-notes.web.cern.ch/system/files/notes/analysis/178/2015-Mar-02-analysis_note-fulltext.pdf
      //   "normal" rho approach is Sec. 5.1.1, sparse event 5.1.2

      AliRhoParameter * outrho = (AliRhoParameter*) reader->GetInputEvent()->FindListObject(fJetRhoTaskName);

      if ( !outrho )
        AliInfo(Form("Could not find rho container <%s>!",fJetRhoTaskName.Data()));
      else
      {
        rhoPtSumTrack = outrho->GetVal() * coneArea;

        coneptsumBkgTrk                = rhoPtSumTrack;
        coneptsumBkgTrkCutMax          = rhoPtSumTrack;
        coneptsumBkgTrkCutLeadFactor   = rhoPtSumTrack;
        coneptsumTrackSub              = coneptsumTrack - coneptsumBkgTrk;
        coneptsumTrackSubCutMax        = coneptsumTrack - coneptsumBkgTrkCutMax;
        coneptsumTrackSubCutLeadFactor = coneptsumTrack - coneptsumBkgTrkCutLeadFactor;

        if ( fPartInCone == kNeutralAndCharged || fPartInCone == kOnlyNeutral )
          coneptsumBkgCls = rhoPtSumTrack*GetNeutralOverChargedRatio(centrality);
        //printf("centrality %f, neutral/charged %f\n",centrality,GetNeutralOverChargedRatio(centrality));

        // Add to candidate object (put for the moment in the perp cone energy container)
        pCandidate->SetChargedPtSumInPerpCone(rhoPtSumTrack);

        if ( fFillHistograms )
        {
          if ( fFillHighMultHistograms )
          {
            fhJetRhoCent->Fill(outrho->GetVal(), centrality, histoWeight);
            fhJetRhoSumPtCent->Fill(ptC, rhoPtSumTrack, centrality, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeSumPtVSJetRhoCent->Fill(coneptsumTrack * excessAreaTrkEta,
                                            rhoPtSumTrack,
                                            centrality, histoWeight);
              fhJetRhoSumPtTrackSubVsNoSubCent->Fill(coneptsumTrack * excessAreaTrkEta,
                                                     coneptsumTrack * excessAreaTrkEta-rhoPtSumTrack,
                                                     centrality, histoWeight);

              if ( fFillEtaPhiHistograms && cenBin < fNCentBins && cenBin >= 0 )
                fhJetRhoSumPtTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, rhoPtSumTrack, histoWeight);
            }
          }
          else
          {
            fhJetRho->Fill(ptC, outrho->GetVal(), histoWeight);
            fhJetRhoSumPt->Fill(ptC, rhoPtSumTrack, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeSumPtVSJetRho->Fill(coneptsumTrack * excessAreaTrkEta,
                                        rhoPtSumTrack,
                                        histoWeight);
              fhJetRhoSumPtTrackSubVsNoSub->Fill(coneptsumTrack * excessAreaTrkEta,
                                                 coneptsumTrack * excessAreaTrkEta-rhoPtSumTrack,
                                                 histoWeight);
              if ( fFillEtaPhiHistograms )
                fhJetRhoSumPtTrigEtaPhi->Fill(etaC, phiC, rhoPtSumTrack, histoWeight);
            }
          }
        }
      }
    }
    else if ( fICMethod >= kSumBkgSubEtaBandIC ) // eta or phi band
    {
      //printf("UE band\n");
      Float_t  etaBandPtSumTrackNorm   = 0;
      Float_t  phiBandPtSumTrackNorm   = 0;
      Float_t  perpBandPtSumTrackNorm  = 0;
      Float_t  etaBandPtSumClusterNorm = 0;
      Float_t  phiBandPtSumClusterNorm = 0;
      
      Float_t  etaBandPtSumTrackCutMaxNorm   = 0;
      Float_t  phiBandPtSumTrackCutMaxNorm   = 0;
      Float_t  perpBandPtSumTrackCutMaxNorm  = 0;
      Float_t  etaBandPtSumClusterCutMaxNorm = 0;
      Float_t  phiBandPtSumClusterCutMaxNorm = 0;
      Float_t  etaBandPtSumClusterCutLeadFactorNorm = 0;
      Float_t  phiBandPtSumClusterCutLeadFactorNorm = 0;
      Float_t  etaBandPtSumTrackCutLeadFactorNorm   = 0;
      Float_t  phiBandPtSumTrackCutLeadFactorNorm   = 0;
      Float_t  perpBandPtSumTrackCutLeadFactorNorm  = 0;

      // Normalize background to cone area
      
      // Check clusters band
      //
      // If DCal or PHOS trigger clusters, it does not make sense, 
      // just use track based estimation later 
      Bool_t checkClustersBand = kTRUE;
      
      if ( fPartInCone == kNeutralAndCharged )
      {
        if      ( pCandidate->GetDetectorTag() != AliFiducialCut::kEMCAL ) 
          checkClustersBand = kFALSE;
        // if declared as emcal because analyze all calo together 
        // but within phos/dcal phi acceptance 
        else if ( phiC > 4.35 && phiC < 5.71 )  
          checkClustersBand = kFALSE; 
        
//        if ( !checkClustersBand ) printf("Cluster band not checked! det %d phi %f\n",
//                                         pCandidate->GetDetectorTag(),phiC);
      }
      
      if ( fPartInCone != kOnlyCharged && checkClustersBand )
      {
        CalculateUEBandClusterNormalization(calorimeter            ,
                                            etaC                   , phiC                  ,
                                            excessClsEta           , excessClsPhi          ,
                                            excessAreaClsEta       , excessAreaClsPhi      ,
                                            etaBandPtSumCluster    , phiBandPtSumCluster   ,
                                            etaBandPtSumClusterNorm, phiBandPtSumClusterNorm);
        
        CalculateUEBandClusterNormalization(calorimeter            ,
                                            etaC                   , phiC                  ,
                                            excessClsEta           , excessClsPhi          ,
                                            excessAreaClsEta       , excessAreaClsPhi      ,
                                            etaBandPtSumClusterCutMax    , phiBandPtSumClusterCutMax   ,
                                            etaBandPtSumClusterCutMaxNorm, phiBandPtSumClusterCutMaxNorm);

        CalculateUEBandClusterNormalization(calorimeter            ,
                                            etaC                   , phiC                  ,
                                            excessClsEta           , excessClsPhi          ,
                                            excessAreaClsEta       , excessAreaClsPhi      ,
                                            etaBandPtSumClusterCutLeadFactor    , phiBandPtSumClusterCutLeadFactor,
                                            etaBandPtSumClusterCutLeadFactorNorm, phiBandPtSumClusterCutLeadFactorNorm);

        // Add to candidate object the normalized band energy
        if      ( fUseMaxPtUE )
        {
          pCandidate->SetNeutralPtSumEtaBand(etaBandPtSumClusterCutMaxNorm);
          pCandidate->SetNeutralPtSumPhiBand(phiBandPtSumClusterCutMaxNorm);
        }
        else if ( fUseLeadingPtUEFactor )
        {
          pCandidate->SetNeutralPtSumEtaBand(etaBandPtSumClusterCutLeadFactorNorm);
          pCandidate->SetNeutralPtSumPhiBand(phiBandPtSumClusterCutLeadFactorNorm);
        }
        else
        {
          pCandidate->SetNeutralPtSumEtaBand(etaBandPtSumClusterNorm);
          pCandidate->SetNeutralPtSumPhiBand(phiBandPtSumClusterNorm);
        }

        if      ( fICMethod == kSumBkgSubEtaBandIC )
        {
          coneptsumBkgClsRaw  = etaBandPtSumCluster; 
          coneptsumBkgCls     = etaBandPtSumClusterNorm; 
          coneptsumBkgClsCutMax        = etaBandPtSumClusterCutMaxNorm;
          coneptsumBkgClsCutLeadFactor = etaBandPtSumClusterCutLeadFactorNorm;
        }
        else  if( fICMethod == kSumBkgSubPhiBandIC )
        {
          coneptsumBkgClsRaw  = phiBandPtSumCluster; 
          coneptsumBkgCls     = phiBandPtSumClusterNorm;
          coneptsumBkgClsCutMax         = phiBandPtSumClusterCutMaxNorm;
          coneptsumBkgClsCutLeadFactor  = phiBandPtSumClusterCutLeadFactorNorm;
        }

        coneptsumClusterSub              = coneptsumCluster - coneptsumBkgCls;
        coneptsumClusterSubCutMax        = coneptsumCluster - coneptsumBkgClsCutMax;
        coneptsumClusterSubCutLeadFactor = coneptsumCluster - coneptsumBkgClsCutLeadFactor;

//        printf("Cluster: sumpT %2.2f, \n \t phi: sum pT %2.2f, sumpT norm %2.2f;\n"
//               "\t eta: sum pT %2.2f, sumpT norm %2.2f;\n"
//               "\t subtracted:  %2.2f\n",
//               coneptsumCluster, 
//               phiBandPtSumCluster, phiBandPtSumClusterNorm, 
//               etaBandPtSumCluster, etaBandPtSumClusterNorm,
//               coneptsumClusterSub);
        
        if ( fFillHistograms )
        {
          Float_t rho           = coneptsumBkgCls / coneArea;
          Float_t rhoMax        = coneptsumBkgClsCutMax / coneArea;
          Float_t rhoLeadFactor = coneptsumBkgClsCutLeadFactor / coneArea;
          if ( fFillHighMultHistograms ) 
          {
            fhConeSumPtUEBandNormClusterCent->Fill(ptC, coneptsumBkgCls*excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeRhoUEBandClusterCent             ->Fill(rho          *excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);
              fhConeRhoUEBandClusterCutMaxCent       ->Fill(rhoMax       *excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);
              fhConeRhoUEBandClusterCutLeadFactorCent->Fill(rhoLeadFactor*excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);

              fhConeSumPtVSUEClusterBandCent->Fill(coneptsumCluster*excessAreaClsEta*excessAreaClsPhi,
                                                   coneptsumBkgCls *excessAreaClsEta*excessAreaClsPhi,
                                                   centrality, histoWeight);
              fhConeSumPtClusterSubVsNoSubCent->Fill(coneptsumCluster*excessAreaClsEta*excessAreaClsPhi,
                                                     coneptsumClusterSub*excessAreaClsEta*excessAreaClsPhi,
                                                     centrality, histoWeight);
            }
          }
          else
          {
            fhConeRhoUEBandCluster             ->Fill(ptC, rho          *excessAreaClsEta*excessAreaClsPhi, histoWeight);
            fhConeRhoUEBandClusterCutMax       ->Fill(ptC, rhoMax       *excessAreaClsEta*excessAreaClsPhi, histoWeight);
            fhConeRhoUEBandClusterCutLeadFactor->Fill(ptC, rhoLeadFactor*excessAreaClsEta*excessAreaClsPhi, histoWeight);

            fhConeSumPtUEBandNormCluster->Fill(ptC, coneptsumBkgCls*excessAreaClsEta*excessAreaClsPhi, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeSumPtVSUEClusterBand->Fill(coneptsumCluster*excessAreaClsEta*excessAreaClsPhi,
                                                   coneptsumBkgCls *excessAreaClsEta*excessAreaClsPhi,
                                                   histoWeight);
              fhConeSumPtClusterSubVsNoSub->Fill(coneptsumCluster*excessAreaClsEta*excessAreaClsPhi,
                                                 coneptsumClusterSub*excessAreaClsEta*excessAreaClsPhi,
                                                 histoWeight);
            }
          }
        } // histograms
      } // clusters in cone
      
      //printf("Pass cluster\n");

      if ( fPartInCone != kOnlyNeutral )
      {
        CalculateUEBandTrackNormalization
        (etaC                 ,
         excessTrkEta         ,
         excessAreaTrkEta     ,
         etaBandPtSumTrack    , phiBandPtSumTrack    , perpBandPtSumTrack,
         etaBandPtSumTrackNorm, phiBandPtSumTrackNorm, perpBandPtSumTrackNorm);
        
        CalculateUEBandTrackNormalization
        (etaC                 ,
         excessTrkEta         ,
         excessAreaTrkEta     ,
         etaBandPtSumTrackCutMax    , phiBandPtSumTrackCutMax    , perpBandPtSumTrackCutMax,
         etaBandPtSumTrackCutMaxNorm, phiBandPtSumTrackCutMaxNorm, perpBandPtSumTrackCutMaxNorm);

        CalculateUEBandTrackNormalization
        (etaC                 ,
         excessTrkEta         ,
         excessAreaTrkEta     ,
         etaBandPtSumTrackCutLeadFactor    , phiBandPtSumTrackCutLeadFactor    , perpBandPtSumTrackCutLeadFactor,
         etaBandPtSumTrackCutLeadFactorNorm, phiBandPtSumTrackCutLeadFactorNorm, perpBandPtSumTrackCutLeadFactorNorm);

        // Add to candidate object the normalized band energy
        if ( fUseMaxPtUE )
        {
          pCandidate->SetChargedPtSumEtaBand(etaBandPtSumTrackCutMaxNorm);
          pCandidate->SetChargedPtSumPhiBand(phiBandPtSumTrackCutMaxNorm);
          pCandidate->SetChargedPtSumInPerpCone(perpBandPtSumTrackCutMaxNorm); // ADD CORRESPONDING METHOD
        }
        else if ( fUseLeadingPtUEFactor )
        {
          pCandidate->SetChargedPtSumEtaBand(etaBandPtSumTrackCutLeadFactorNorm);
          pCandidate->SetChargedPtSumPhiBand(phiBandPtSumTrackCutLeadFactorNorm);
          pCandidate->SetChargedPtSumInPerpCone(perpBandPtSumTrackCutLeadFactorNorm); // ADD CORRESPONDING METHOD
        }
        else
        {
          pCandidate->SetChargedPtSumEtaBand(etaBandPtSumTrackNorm);
          pCandidate->SetChargedPtSumPhiBand(phiBandPtSumTrackNorm);
          pCandidate->SetChargedPtSumInPerpCone(perpBandPtSumTrackNorm); // ADD CORRESPONDING METHOD
        }

        if      ( fICMethod == kSumBkgSubEtaBandIC )
        {
          ueptLeadTrack      = etaBandptLeadTrack;
          coneptsumBkgTrkRaw = etaBandPtSumTrack; 
          coneptsumBkgTrk    = etaBandPtSumTrackNorm;
          coneptsumBkgTrkCutMax        = etaBandPtSumTrackCutMaxNorm;
          coneptsumBkgTrkCutLeadFactor = etaBandPtSumTrackCutLeadFactorNorm;
        }
        else  if( fICMethod == kSumBkgSubPhiBandIC )
        {
          ueptLeadTrack      = phiBandptLeadTrack;
          coneptsumBkgTrkRaw = phiBandPtSumTrack; 
          coneptsumBkgTrk    = phiBandPtSumTrackNorm;
          coneptsumBkgTrkCutMax        = phiBandPtSumTrackCutMaxNorm;
          coneptsumBkgTrkCutLeadFactor = phiBandPtSumTrackCutLeadFactorNorm;
        }
        else  if( fICMethod == kSumBkgSubPerpBandIC )
        {
          ueptLeadTrack      = perpBandptLeadTrack;
          coneptsumBkgTrkRaw = perpBandPtSumTrack;
          coneptsumBkgTrk    = perpBandPtSumTrackNorm;
          coneptsumBkgTrkCutMax        = perpBandPtSumTrackCutMaxNorm;
          coneptsumBkgTrkCutLeadFactor = perpBandPtSumTrackCutLeadFactorNorm;
        }
        
        coneptsumTrackSub              = coneptsumTrack - coneptsumBkgTrk;
        coneptsumTrackSubCutMax        = coneptsumTrack - coneptsumBkgTrkCutMax;
        coneptsumTrackSubCutLeadFactor = coneptsumTrack - coneptsumBkgTrkCutLeadFactor;

//        printf("Track: sumpT %2.2f, \n \t phi: sum pT %2.2f, sumpT norm %2.2f;\n"
//               "\t eta: sum pT %2.2f, sumpT norm %2.2f;\n"
//               "\t subtracted: %2.2f\n",
//               coneptsumTrack, 
//               phiBandPtSumTrack, phiBandPtSumTrackNorm, 
//               etaBandPtSumTrack, etaBandPtSumTrackNorm, 
//               coneptsumTrackSub);  
        
        if ( fFillHistograms )
        {          
          Float_t rho           = coneptsumBkgTrk / coneArea;
          Float_t rhoMax        = coneptsumBkgTrkCutMax / coneArea;
          Float_t rhoLeadFactor = coneptsumBkgTrkCutLeadFactor / coneArea;
          if ( fFillHighMultHistograms )
          {
            fhConeSumPtUEBandNormTrackCent->Fill(ptC, coneptsumBkgTrk*excessAreaTrkEta, centrality, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeRhoUEBandTrackCent             ->Fill(rho          *excessAreaTrkEta, centrality, histoWeight);
              fhConeRhoUEBandTrackCutMaxCent       ->Fill(rhoMax       *excessAreaTrkEta, centrality, histoWeight);
              fhConeRhoUEBandTrackCutLeadFactorCent->Fill(rhoLeadFactor*excessAreaTrkEta, centrality, histoWeight);

              fhConeSumPtVSUETrackBandCent->Fill(coneptsumTrack *excessAreaClsEta,
                                                 coneptsumBkgTrk*excessAreaClsEta,
                                                 centrality, histoWeight);
              fhConeSumPtTrackSubVsNoSubCent->Fill(coneptsumTrack*excessAreaTrkEta,
                                                   coneptsumTrackSub*excessAreaTrkEta,
                                                   centrality, histoWeight);
            }
          }
          else
          {
            fhConeRhoUEBandTrack             ->Fill(ptC, rho          *excessAreaTrkEta, histoWeight);
            fhConeRhoUEBandTrackCutMax       ->Fill(ptC, rhoMax       *excessAreaTrkEta, histoWeight);
            fhConeRhoUEBandTrackCutLeadFactor->Fill(ptC, rhoLeadFactor*excessAreaTrkEta, histoWeight);

            fhConeSumPtUEBandNormTrack->Fill(ptC, coneptsumBkgTrk*excessAreaTrkEta, histoWeight);
            if ( ptC > fEtaPhiHistogramsMinPt )
            {
              fhConeSumPtVSUETrackBand->Fill(coneptsumTrack *excessAreaClsEta,
                                             coneptsumBkgTrk*excessAreaClsEta,
                                             histoWeight);
              fhConeSumPtTrackSubVsNoSub->Fill(coneptsumTrack*excessAreaTrkEta,
                                               coneptsumTrackSub*excessAreaTrkEta,
                                               histoWeight);
            }
          }
        } // fill 
      } // tracks in cone

      // In case of DCal and PHOS, estimate the bkg in clusters from tracks
      //
      if ( fPartInCone == kNeutralAndCharged && !checkClustersBand )
      {
        coneptsumBkgCls     = coneptsumBkgTrk*GetNeutralOverChargedRatio(centrality); 
        coneptsumClusterSub = coneptsumCluster - coneptsumBkgCls;
      }
      
//     // Uncomment if perp cone is hacked to be calculated
//      printf("UE BKG per method:\n \t Perp Cone %2.2f (cl %2.2f, tr %2.2f)\n \t Phi Band %2.2f (cl %2.2f, tr %2.2f)\n \t Eta Band %2.2f (cl %2.2f, tr %2.2f)\n",
//             perpPtSumTrack+perpPtSumTrack*fNeutralOverChargedRatio, perpPtSumTrack*fNeutralOverChargedRatio,perpPtSumTrack,
//             phiBandPtSumTrackNorm+phiBandPtSumClusterNorm,phiBandPtSumClusterNorm,phiBandPtSumTrackNorm,
//             etaBandPtSumTrackNorm+etaBandPtSumClusterNorm,etaBandPtSumClusterNorm,etaBandPtSumTrackNorm);
      
      //printf("Pass track\n");

      if ( fPartInCone == AliIsolationCut::kNeutralAndCharged && 
           fFillHistograms )
      {
        if ( !fFillHighMultHistograms )
        {
          fhBandClustervsTrack    ->Fill(coneptsumBkgClsRaw, coneptsumBkgTrkRaw, histoWeight);
          fhBandNormClustervsTrack->Fill(coneptsumBkgCls   , coneptsumBkgTrk   , histoWeight);

          fhConeSumPtUEBandSubClustervsTrack->Fill(coneptsumClusterSub, coneptsumTrackSub, histoWeight);
        }
        else
        {
          fhBandClustervsTrackCent    ->Fill(coneptsumBkgClsRaw, coneptsumBkgTrkRaw, centrality, histoWeight);
          fhBandNormClustervsTrackCent->Fill(coneptsumBkgCls   , coneptsumBkgTrk   , centrality, histoWeight);

          fhConeSumPtUEBandSubClustervsTrackCent->Fill(coneptsumClusterSub, coneptsumTrackSub, centrality, histoWeight);
        }
      }
      //printf("Pass both\n");
      
    } // UE subtraction by Eta / Phi / Perp bands
    
    if ( fUseMaxPtUE )
    {
      coneptsumUESubCluster -= coneptsumBkgClsCutMax;
      coneptsumUESubTrack   -= coneptsumBkgTrkCutMax;
    }
    else if ( fUseLeadingPtUEFactor )
    {
      coneptsumUESubCluster -= coneptsumBkgClsCutLeadFactor;
      coneptsumUESubTrack   -= coneptsumBkgTrkCutLeadFactor;
    }
    else
    {
      coneptsumUESubCluster -= coneptsumBkgCls;
      coneptsumUESubTrack   -= coneptsumBkgTrk;
    }

//    printf("method %d, track cone %2.2f, no cut %2.2f, max cut %2.2f, lead cut %2.2f; sub %2.2f\n",
//           fICMethod, coneptsumTrack, coneptsumBkgTrk, coneptsumBkgTrkCutMax, coneptsumBkgTrkCutLeadFactor, coneptsumUESubTrack);
    
    Float_t fracLeadTrack = 0;
    if ( ueptLeadTrack > 0 )
      fracLeadTrack = coneptLeadTrack/ueptLeadTrack;
    Float_t fracLeadCluster = 0;
    if ( ueptLeadCluster > 0 )
      fracLeadCluster = coneptLeadTrack/ueptLeadCluster;

//    printf("method %d, Cen %2.0f; cone sum %2.2f, sub %2.2f; UE %2.2f Rho %2.2f; Cut sub %2.2f; UE %2.2f; Lead cone %2.2f, UE %2.2f, frac %2.2f\n",
//           fICMethod, centrality, coneptsumTrack,
//           coneptsumTrack - coneptsumBkgTrk, coneptsumBkgTrk,
//           coneptsumBkgTrk / ( TMath::Pi() * fConeSize * fConeSize),
//           coneptsumTrack - coneptsumBkgTrkCut, coneptsumBkgTrkCut,
//           coneptLeadTrack,ueptLeadTrack, fracLeadTrack);

    if ( fFillHistograms )
    {
      if ( fPartInCone != kOnlyNeutral  )
      {
        if ( !fFillHighMultHistograms )
        {
          if ( fFillLeadingVsUESubHisto )
          {
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePt->Fill(ptC, coneptLeadTrack,
                                                              coneptsumUESubTrack*excessAreaTrkEta,
                                                              histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePt->Fill(ptC, ueptLeadTrack,
                                                                coneptsumUESubTrack*excessAreaTrkEta,
                                                                histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePt->Fill(ptC, fracLeadTrack,
                                                                  coneptsumUESubTrack*excessAreaTrkEta,
                                                                  histoWeight);
            if ( ueptLeadTrack < coneptLeadTrack )
            {
              fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELead->Fill(ptC, coneptLeadTrack,
                                                                         coneptsumUESubTrack*excessAreaTrkEta,
                                                                         histoWeight);
            }
          }

          fhConeSumPtUESubTrackCutMax       ->Fill(ptC, coneptsumTrackSubCutMax        * excessAreaTrkEta, histoWeight);
          fhConeSumPtUESubTrackCutLeadFactor->Fill(ptC, coneptsumTrackSubCutLeadFactor * excessAreaTrkEta, histoWeight);
        }
        else if ( cenBin < fNCentBins && cenBin >= 0 )
        {
          if ( fFillLeadingVsUESubHisto )
          {
            fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtCent[cenBin]->Fill(ptC, coneptLeadTrack,
                                                                          coneptsumUESubTrack*excessAreaTrkEta,
                                                                          histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadUETrackInConePtCent[cenBin]->Fill(ptC, ueptLeadTrack,
                                                                            coneptsumUESubTrack*excessAreaTrkEta,
                                                                            histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadTrackFracInConePtCent[cenBin]->Fill(ptC, fracLeadTrack,
                                                                              coneptsumUESubTrack*excessAreaTrkEta,
                                                                              histoWeight);
            if ( ueptLeadTrack < coneptLeadTrack )
            {
              fhTrigPtVsSumPtUEBandSubVsLeadTrackInConePtLowUELeadCent[cenBin]->Fill(ptC, coneptLeadTrack,
                                                                                     coneptsumUESubTrack*excessAreaTrkEta,
                                                                                     histoWeight);
            }
          }

          fhConeSumPtUESubTrackCutMaxCent       ->Fill(ptC, coneptsumTrackSubCutMax        * excessAreaTrkEta, centrality, histoWeight);
          fhConeSumPtUESubTrackCutLeadFactorCent->Fill(ptC, coneptsumTrackSubCutLeadFactor * excessAreaTrkEta, centrality, histoWeight);
        }
      }

      if ( fPartInCone != kOnlyCharged )
      {
        if ( !fFillHighMultHistograms )
        {
          if ( fFillLeadingVsUESubHisto )
          {
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePt->Fill(ptC, coneptLeadCluster,
                                                                coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePt->Fill(ptC, ueptLeadCluster,
                                                                  coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                  histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePt->Fill(ptC, fracLeadCluster,
                                                                    coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                    histoWeight);
            if ( ueptLeadCluster < coneptLeadCluster )
            {
              fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELead->Fill(ptC, coneptLeadCluster,
                                                                           coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                           histoWeight);
            }
          }

          fhConeSumPtUESubClusterCutMax       ->Fill(ptC, coneptsumClusterSubCutMax        * excessAreaClsEta*excessAreaClsPhi, histoWeight);
          fhConeSumPtUESubClusterCutLeadFactor->Fill(ptC, coneptsumClusterSubCutLeadFactor * excessAreaClsEta*excessAreaClsPhi, histoWeight);
        }
        else if ( cenBin < fNCentBins && cenBin >= 0 )
        {
          if ( fFillLeadingVsUESubHisto )
          {
            fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtCent[cenBin]->Fill(ptC, coneptLeadCluster,
                                                                            coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                            histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadUEClusterInConePtCent[cenBin]->Fill(ptC, ueptLeadCluster,
                                                                              coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                              histoWeight);
            fhTrigPtVsSumPtUEBandSubVsLeadClusterFracInConePtCent[cenBin]->Fill(ptC, fracLeadCluster,
                                                                                coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                                histoWeight);
            if ( ueptLeadCluster < coneptLeadCluster )
            {
              fhTrigPtVsSumPtUEBandSubVsLeadClusterInConePtLowUELeadCent[cenBin]->Fill(ptC, coneptLeadCluster,
                                                                                       coneptsumUESubCluster*excessAreaClsEta*excessAreaClsPhi,
                                                                                       histoWeight);
            }
          }

          fhConeSumPtUESubClusterCutMaxCent       ->Fill(ptC, coneptsumClusterSubCutMax        * excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);
          fhConeSumPtUESubClusterCutLeadFactorCent->Fill(ptC, coneptsumClusterSubCutLeadFactor * excessAreaClsEta*excessAreaClsPhi, centrality, histoWeight);
        }
      }
    }

    // Calculated in case of fICMethod == kSumBkgSubEtaBandIC, 
    // reset excess areas to 1 if not used in final result.
    if ( !fMakeConeExcessCorr )
    {
      excessAreaTrkEta = 1;
      excessAreaClsEta = 1;
      excessAreaClsPhi = 1;
    }
    
    if      ( fPartInCone == kNeutralAndCharged )  coneptsumUESub  = coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta +
                                                                     coneptsumUESubTrack   * excessAreaTrkEta;
    else if ( fPartInCone == kOnlyCharged       )  coneptsumUESub  = coneptsumUESubTrack   * excessAreaTrkEta;
    else if ( fPartInCone == kOnlyNeutral       )  coneptsumUESub  = coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta;

    if ( coneptsumUESub > fSumPtThreshold && coneptsumUESub < fSumPtThresholdMax )
      isolated  =  kFALSE ;
    else
      isolated  =  kTRUE  ;
    //printf("isolated %d\n",isolated);

  } // UE subtraction, different options
  
  // Put final subtracted and corrected (if requested, if not it will be without the correction)
  pCandidate->SetNeutralPtSumInCone (coneptsumUESubCluster *excessAreaClsEta * excessAreaClsPhi);
  pCandidate->SetChargedPtSumInCone (coneptsumUESubTrack   *excessAreaTrkEta);
  
//  printf("MakeIso:: Method %d cone %0.2f: sumPt Ch %2.1f, Ne %2.1f; "
//         "eta band Ch %2.1f, Ne %2.1f; phi band Ch %2.1f, Ne %2.1f; perp Ch %2.1f \n ",
//         fICMethod, fConeSize,coneptsumTrack ,coneptsumCluster ,
//         etaBandPtSumTrack,etaBandPtSumCluster,
//         phiBandPtSumTrack,phiBandPtSumCluster,perpPtSumTrack);

  //-------------------------------------------------------------------
  // Fill histograms
  //-------------------------------------------------------------------

  if ( !fFillHistograms ) return;

  // Sum pt in cone
  //
  //printf("fill histo");
  if ( fFillHighMultHistograms )
  {
    fhConeSumPtCent->Fill(ptC, coneptsum, centrality, histoWeight);
    if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
      fhConeSumPtTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, coneptsum, histoWeight);
  }
  else
  {
    fhConeSumPt->Fill(ptC, coneptsum, histoWeight);
    if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
      fhConeSumPtTrigEtaPhi->Fill(etaC, phiC, coneptsum, histoWeight);
  }
  
  if ( fPartInCone == kNeutralAndCharged && fICMethod != kSumBkgSubIC ) // No need for perpendicular or charged/neutral only analysis
  {
    if ( fFillHighMultHistograms )
      fhConeSumPtClustervsTrackCent->Fill(coneptsumCluster, coneptsumTrack, centrality, histoWeight);
    else
      fhConeSumPtClustervsTrack ->Fill(coneptsumCluster, coneptsumTrack, histoWeight);

    if ( coneptsumTrack > 0) 
    {
      if ( fFillHighMultHistograms )
        fhConeSumPtClusterTrackFracCent->Fill(ptC, coneptsumCluster /coneptsumTrack, centrality, histoWeight);
      else
        fhConeSumPtClusterTrackFrac->Fill(ptC, coneptsumCluster /coneptsumTrack, histoWeight);
    }
  }
  
  // Here the sum in cone before subtraction, if done, to check the effect.
  //
  if ( fICMethod >= kSumBkgSubIC )
  {
    if ( fFillHighMultHistograms )
    {
      fhConeSumPtUESubCent->Fill(ptC, coneptsumUESub, centrality, histoWeight);
      if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
        fhConeSumPtUESubTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, coneptsumUESub, histoWeight);
    }
    else
    {
      fhConeSumPtUESub ->Fill(ptC, coneptsumUESub, histoWeight);
      if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
        fhConeSumPtUESubTrigEtaPhi->Fill(etaC, phiC, coneptsumUESub, histoWeight);
    }
    
    // No for charged/neutral only analysis
    if ( fPartInCone == kNeutralAndCharged )
    {
      if ( fFillHighMultHistograms ) 
      {
        fhConeSumPtUESubTrackCent  ->Fill(ptC, coneptsumUESubTrack * excessAreaTrkEta  , centrality, histoWeight);
        fhConeSumPtUESubClusterCent->Fill(ptC, coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, centrality, histoWeight);
        
        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt && cenBin < fNCentBins && cenBin >= 0 )
        {
          fhConeSumPtUESubTrackTrigEtaPhiCent  [cenBin]->Fill(etaC, phiC, coneptsumUESubTrack * excessAreaTrkEta  , histoWeight);
          fhConeSumPtUESubClusterTrigEtaPhiCent[cenBin]->Fill(etaC, phiC, coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, histoWeight);
        }
      }
      else
      {
        fhConeSumPtUESubTrack  ->Fill(ptC, coneptsumUESubTrack * excessAreaTrkEta, histoWeight);
        fhConeSumPtUESubCluster->Fill(ptC, coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, histoWeight);
        
        if ( fFillEtaPhiHistograms && ptC > fEtaPhiHistogramsMinPt )
        {
          fhConeSumPtUESubTrackTrigEtaPhi  ->Fill(etaC, phiC, coneptsumUESubTrack * excessAreaTrkEta, histoWeight);
          fhConeSumPtUESubClusterTrigEtaPhi->Fill(etaC, phiC, coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, histoWeight);
        }
      }
      
      // No need for perpendicular cones
      if ( fICMethod != kSumBkgSubIC )
      {
        if ( fFillHighMultHistograms )
          fhConeSumPtUESubClustervsTrackCent->Fill(coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, 
                                                   coneptsumUESubTrack * excessAreaTrkEta, centrality, histoWeight);
        else
          fhConeSumPtUESubClustervsTrack ->Fill(coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta, 
                                                coneptsumUESubTrack * excessAreaTrkEta, histoWeight);

        if ( TMath::Abs(coneptsumUESubTrack) > 0 ) 
        {
           if ( fFillHighMultHistograms ) 
             fhConeSumPtUESubClusterTrackFracCent->Fill(ptC, (coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta) / (coneptsumUESubTrack * excessAreaTrkEta), centrality, histoWeight);
          else
            fhConeSumPtUESubClusterTrackFrac->Fill(ptC, (coneptsumUESubCluster * excessAreaClsPhi*excessAreaClsEta) / (coneptsumUESubTrack * excessAreaTrkEta), histoWeight);
        }
      }
    }
  }

  // Leading in cone
  //
  Float_t coneptLead = coneptLeadTrack;
  if(coneptLeadCluster > coneptLeadTrack) 
    coneptLead = coneptLeadCluster;
  
  if ( !fFillHighMultHistograms && fFillLeadingHistograms )
  {
    if ( fPartInCone == kNeutralAndCharged )
    {
      fhConePtLeadClustervsTrack->Fill(coneptLeadCluster,coneptLeadTrack, histoWeight);

      if ( coneptLeadTrack > 0 )
        fhConePtLeadClusterTrackFrac->Fill(ptC, coneptLeadCluster/coneptLeadTrack, histoWeight);
    }

    fhConePtLead->Fill(ptC, coneptLead, histoWeight);
  }
}

//_____________________________________________________
/// Print some relevant parameters set for the analysis.
//_____________________________________________________
void AliIsolationCut::Print(const Option_t * opt) const
{
  if(! opt)
    return;

  printf("**** Print %s %s **** \n", GetName(), GetTitle() ) ;

  printf("IC method          =     %d\n",    fICMethod   ) ;
  if ( fICMethod == kSumBkgSubJetRhoIC )
    printf("Jet rho task name = %s\n",fJetRhoTaskName.Data());
  printf("Cone Size          =     %1.2f\n", fConeSize   ) ;
  printf("Cone Size UE Gap   =     %1.2f\n", fConeSizeBandGap ) ;
  printf("pT threshold       =     >%2.1f;<%2.1f\n", fPtThreshold, fPtThresholdMax) ;
  printf("Sum pT threshold   =     >%2.1f;<%2.1f, gap = %2.1f\n", fSumPtThreshold, fSumPtThresholdMax, fSumPtThresholdGap) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction ) ;
  printf("particle type in cone =  %d\n",    fPartInCone ) ;
  printf("using fraction for high pt leading instead of frac ? %i\n",fFracIsThresh);
  printf("minimum distance to candidate, R>%1.2f\n",fDistMinToTrigger);
  printf("leading pT UE factor = %1.2f, UE pt < %2.1f\n",fLeadingPtUEFactor,fMaxPtUE);
  printf("correct cone excess = %d \n",fMakeConeExcessCorr);
  printf("NeutralOverChargedRatio param={%1.2e,%1.2e,%1.2e,%1.2e} \n",
  fNeutralOverChargedRatio[0],fNeutralOverChargedRatio[1],fNeutralOverChargedRatio[2],fNeutralOverChargedRatio[3]) ;
  printf("    \n") ;
}

//______________________________________________________________
/// Calculate the distance to trigger from any particle.
/// \param etaC: pseudorapidity of candidate particle.
/// \param phiC: azimuthal angle of candidate particle.
/// \param eta: pseudorapidity of track/cluster to be considered in cone.
/// \param phi: azimuthal angle of track/cluster to be considered in cone.
//______________________________________________________________
Float_t AliIsolationCut::Radius(Float_t etaC, Float_t phiC,
                                Float_t eta , Float_t phi) const
{
  Float_t dEta = etaC-eta;
  Float_t dPhi = phiC-phi;
  
  if(TMath::Abs(dPhi) >= TMath::Pi())
    dPhi = TMath::TwoPi()-TMath::Abs(dPhi);

  return TMath::Sqrt( dEta*dEta + dPhi*dPhi );
}



