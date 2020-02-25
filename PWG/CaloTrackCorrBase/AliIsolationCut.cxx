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

// --- AliRoot system ---
#include "AliCaloTrackParticleCorrelation.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliMixedEvent.h"
#include "AliLog.h"

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
fFillHistograms(0),  fFillEtaPhiHistograms(0),      fMakeConeExcessCorr(0),
fConeSize(0.),       fPtThreshold(0.),              fPtThresholdMax(10000.),
fSumPtThreshold(0.), fSumPtThresholdMax(10000.),    fSumPtThresholdGap(0.),
fPtFraction(0.),     fICMethod(0),                  fPartInCone(0),
fFracIsThresh(1),    fIsTMClusterInConeRejected(1), fDistMinToTrigger(-1.),
fNeutralOverChargedRatio(0),
fDebug(0),           fMomentum(),                   fTrackVector(),
fEMCEtaSize(-1),     fEMCPhiMin(-1),                fEMCPhiMax(-1),
fTPCEtaSize(-1),     fTPCPhiSize(-1),
// Histograms
fHistoRanges(0),
fhPtInCone(0),       
fhPtClusterInCone(0),                       fhPtTrackInCone(0),  
fhConeSumPt(0),      
fhConeSumPtCluster(0),                      fhConeSumPtTrack(0),
fhConeSumPtClustervsTrack(0),               fhConeSumPtClusterTrackFrac(0), 
fhConeSumPtTrigEtaPhi(0),
fhConeSumPtUESub(0),      
fhConeSumPtUESubCluster(0),                 fhConeSumPtUESubTrack(0),
fhConeSumPtUESubClustervsTrack(0),          fhConeSumPtUESubClusterTrackFrac(0), 
fhConeSumPtUESubTrigEtaPhi(0),
fhConePtLead(0),     
fhConePtLeadCluster(0),                     fhConePtLeadTrack(0),
fhConePtLeadClustervsTrack(0),              fhConePtLeadClusterTrackFrac(0),
fhEtaPhiCluster(0),                         fhEtaPhiTrack(0),
fhEtaPhiInConeCluster(0),                   fhEtaPhiInConeTrack(0),   
fhPtInPerpCone(0),                          fhPerpConeSumPt(0),
fhEtaPhiInPerpCone(0),                      fhConeSumPtVSPerpCone(0), 
fhPerpConeSumPtTrigEtaPhi(0),
fhEtaBandClusterPt(0),                      fhPhiBandClusterPt(0),
fhEtaBandTrackPt(0),                        fhPhiBandTrackPt(0),
fhConeSumPtEtaBandUECluster(0),             fhConeSumPtPhiBandUECluster(0),
fhConeSumPtEtaBandUETrack(0),               fhConeSumPtPhiBandUETrack(0),
fhEtaBandClusterEtaPhi(0),                  fhPhiBandClusterEtaPhi(0),
fhEtaBandTrackEtaPhi(0),                    fhPhiBandTrackEtaPhi(0),
fhConeSumPtEtaBandUEClusterTrigEtaPhi(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhi(0),
fhConeSumPtEtaBandUETrackTrigEtaPhi(0),     fhConeSumPtPhiBandUETrackTrigEtaPhi(0),
fhConeSumPtVSUETracksEtaBand(0),            fhConeSumPtVSUETracksPhiBand(0),
fhConeSumPtVSUEClusterEtaBand(0),           fhConeSumPtVSUEClusterPhiBand(0),
fhConeSumPtUEBandNormCluster(0),            fhConeSumPtUEBandNormTrack(0), 
fhFractionTrackOutConeEta(0),               fhFractionTrackOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConeEta(0),             fhFractionClusterOutConeEtaTrigEtaPhi(0),
fhFractionClusterOutConePhi(0),             fhFractionClusterOutConePhiTrigEtaPhi(0),
fhFractionClusterOutConeEtaPhi(0),          fhFractionClusterOutConeEtaPhiTrigEtaPhi(0),
fhConeSumPtUEBandSubClustervsTrack(0),       
fhBandClustervsTrack(0),                    fhBandNormClustervsTrack(0),             
fhConeSumPtTrackSubVsNoSub(0),              fhConeSumPtClusterSubVsNoSub(0)  
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
/// \param histoWeight: Histograms weight (event, pt depedent)
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateCaloSignalInCone 
(
 AliCaloTrackParticleCorrelation * pCandidate, AliCaloTrackReader * reader,
 Bool_t    bFillAOD           , Bool_t    useRefs, 
 TString   aodArrayRefName    , TObjArray  * bgCls,
 Int_t     calorimeter        , AliCaloPID * pid,
 Int_t   & nPart              , Int_t   & nfrac,
 Float_t & coneptsumCluster   , Float_t & coneptLeadCluster,
 Float_t & etaBandPtSumCluster, Float_t & phiBandPtSumCluster, 
 Double_t histoWeight
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
        fhConeSumPtCluster ->Fill(pCandidate->Pt(), 0., histoWeight);
      
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
  
  // Get the clusters
  //
  //printf("Loop calo\n");
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
        if( fPartInCone == kNeutralAndCharged &&
           pid->IsTrackMatched(calo,reader->GetCaloUtils(),reader->GetInputEvent()) ) continue ;
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
    
    // ** Exclude clusters too close to the candidate, inactive by default **
    
    if ( rad < fDistMinToTrigger ) continue ;
    
    // ** For the background out of cone **
    
    if ( fFillHistograms && fFillEtaPhiHistograms )
      fhEtaPhiCluster->Fill(eta, phi, histoWeight);

    if ( fICMethod >= kSumBkgSubIC && rad > fConeSize )
    {
      // phi band
      if(eta > (etaC-fConeSize) && eta < (etaC+fConeSize))
      {
        phiBandPtSumCluster += pt;
        
        if ( fFillHistograms )
        {
          fhPhiBandClusterPt->Fill(ptC, pt, histoWeight);  
          
          if ( fFillEtaPhiHistograms )
            fhPhiBandClusterEtaPhi->Fill(eta, phi, histoWeight);
        }
      } // phi band
      
      // eta band
      if(phi > (phiC-fConeSize) && phi < (phiC+fConeSize))
      {
        etaBandPtSumCluster += pt;
        
        if ( fFillHistograms )
        {
          fhEtaBandClusterPt->Fill(ptC, pt, histoWeight);  
          
          if ( fFillEtaPhiHistograms )
            fhEtaBandClusterEtaPhi->Fill(eta, phi, histoWeight);
        }
      } // eta band
    } // out of cone
    
    //-------------------------------------------------------------
    // ** For the isolated particle **
    //-------------------------------------------------------------
    
    // Only loop the particle at the same side of candidate
    if ( TMath::Abs(phi-phiC)>TMath::PiOver2() ) continue ;
    
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
    
    AliDebug(2,Form("\t Cluster %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f", 
                    ipr,pt,eta,phi,rad));
    
    //
    // Select calorimeter clusters inside the isolation radius 
    //
    if ( rad > fConeSize ) continue ;
    
    AliDebug(2,"Inside candidate cone");
    
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
      fhPtInCone->Fill(ptC, pt, histoWeight);
      
      if ( fPartInCone == kNeutralAndCharged )
        fhPtClusterInCone->Fill(ptC, pt, histoWeight);
      
      if ( fFillEtaPhiHistograms )
        fhEtaPhiInConeCluster->Fill(eta, phi, histoWeight);
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
    
  }// neutral particle loop
  //printf("end loop, fill histo\n");
  if ( fFillHistograms ) // Do not fill in perp cones case
  {
    //printf("A\n");
    if( fPartInCone == kNeutralAndCharged )
    {
      fhConeSumPtCluster ->Fill(ptC, coneptsumCluster , histoWeight);
      fhConePtLeadCluster->Fill(ptC, coneptLeadCluster, histoWeight);
    }  
    
    //printf("B\n");

    // UE substraction
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      fhConeSumPtEtaBandUECluster->Fill(ptC, etaBandPtSumCluster, histoWeight);
      fhConeSumPtPhiBandUECluster->Fill(ptC, phiBandPtSumCluster, histoWeight);
      
      //printf("C\n");

      fhConeSumPtVSUEClusterEtaBand->Fill(coneptsumCluster, etaBandPtSumCluster, histoWeight);
      fhConeSumPtVSUEClusterPhiBand->Fill(coneptsumCluster, phiBandPtSumCluster, histoWeight);
      
      //printf("A\n");
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaC, phiC, etaBandPtSumCluster*histoWeight); // Check
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaC, phiC, phiBandPtSumCluster*histoWeight); // Check
      }
    } // UE sub
  } // fill histo
  
  // Add calculated values to pCandidate, might be modified later 
  // if UE subtraction and normalization requested
  
  // Add reference arrays to AOD when filling AODs only
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
/// \param coneptsumTrack: total momentum energy in cone (track+cluster), output.
/// \param coneptLeadTrack: momentum of leading cluster or track in cone, output.
/// \param etaBandPtSumTrack: sum of tracks in eta band, same phi as candidate
/// \param phiBandPtSumTrack: sum of tracks in phi band, same eta as candidate (not useful need to restrict it)
/// \param perpConePtSumTrack: sum of tracks in perpendicular cones in phi, return result divided by 2
/// \param histoWeight: Histograms weight (event, pt depedent)
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateTrackSignalInCone
(
 AliCaloTrackParticleCorrelation * pCandidate, AliCaloTrackReader * reader,
 Bool_t    bFillAOD         , Bool_t    useRefs, 
 TString   aodArrayRefName  , TObjArray * bgTrk,
 Int_t   & nPart            , Int_t   & nfrac,
 Float_t & coneptsumTrack   , Float_t & coneptLeadTrack,
 Float_t & etaBandPtSumTrack, Float_t & phiBandPtSumTrack,
 Float_t & perpConePtSumTrack, Double_t histoWeight
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
    if ( fFillHistograms )
    {
      fhConeSumPtTrack      ->Fill(pCandidate->Pt(), 0., histoWeight);
      
      fhConePtLeadTrack     ->Fill(pCandidate->Pt(), 0., histoWeight);
    } // histograms
    
    if(coneptLeadTrack > 0  || coneptsumTrack > 0) 
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
  //
  //-----------------------------------------------------------
  for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ )
  {
    AliVTrack* track = dynamic_cast<AliVTrack*>(plCTS->At(ipr)) ;
    
    if(track)
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
    
    // ** For the background out of cone **
    
    if ( fFillHistograms && fFillEtaPhiHistograms )
      fhEtaPhiTrack->Fill(etaTrack, phiTrack, histoWeight);
    
    if ( fICMethod >= kSumBkgSubIC && rad > fConeSize )
    {
      // Phi band
      if ( etaTrack > (etaTrig-fConeSize) && etaTrack < (etaTrig+fConeSize) &&
          TMath::Abs(phiTrig-phiTrack) < TMath::PiOver2() )  // Look only half TPC with respect candidate, avoid opposite side jet 
      {
        phiBandPtSumTrack += ptTrack;
        
        if ( fFillHistograms )
        {
          fhPhiBandTrackPt->Fill(ptTrig , ptTrack, histoWeight);
          
          if ( fFillEtaPhiHistograms )
            fhPhiBandTrackEtaPhi->Fill(etaTrig, phiTrig, histoWeight);
        }
      } // phi band
      
      // Eta band
      if ( phiTrack > (phiTrig-fConeSize) && phiTrack < (phiTrig+fConeSize) ) 
      {
        etaBandPtSumTrack += ptTrack;
        
        if ( fFillHistograms )
        {
          fhEtaBandTrackPt->Fill(ptTrig , ptTrack, histoWeight);
          
          if ( fFillEtaPhiHistograms )
            fhEtaBandTrackEtaPhi->Fill(etaTrig, phiTrig, histoWeight);
        }
      } // eta band
      
    } // out of cone
      
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fill the histograms at +-45 degrees in phi  
    // from trigger particle, perpedicular to trigger axis in phi
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( fICMethod == kSumBkgSubIC )
    {
      Double_t dEta    = etaTrig - etaTrack;
      
      Double_t dPhiPlu = phiTrig - phiTrack + TMath::PiOver2();
      Double_t dPhiMin = phiTrig - phiTrack - TMath::PiOver2();
      
      Double_t argPlu  = dPhiPlu*dPhiPlu + dEta*dEta;
      Double_t argMin  = dPhiMin*dPhiMin + dEta*dEta;
      
      Bool_t fillPerp = kFALSE;
      if ( TMath::Sqrt(argPlu) < fConeSize ) fillPerp = kTRUE ;
      if ( TMath::Sqrt(argMin) < fConeSize ) fillPerp = kTRUE ;
      
      if ( fillPerp )
      {
        perpConePtSumTrack+=ptTrack;
        
        if ( fFillHistograms )
        {
          fhPtInPerpCone->Fill(ptTrig, ptTrack, histoWeight);
          
          if ( fFillEtaPhiHistograms )
            fhEtaPhiInPerpCone->Fill(etaTrack,phiTrack, histoWeight);
        }
      }
    }
    //-------------------------------------------------------------
    // ** For the isolated particle **
    //-------------------------------------------------------------

    // Only loop the particle at the same side of candidate
    if ( TMath::Abs(phiTrack-phiTrig) > TMath::PiOver2() ) continue ;
    
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
    
    AliDebug(2,Form("\t Track %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f", 
                    ipr,ptTrack,etaTrack,phiTrack,rad));
    
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
      fhPtInCone->Fill(ptTrig , ptTrack, histoWeight);
      
      if( fPartInCone == kNeutralAndCharged )
        fhPtTrackInCone->Fill(ptTrig , ptTrack, histoWeight);
        
      if ( fFillEtaPhiHistograms )
        fhEtaPhiInConeTrack->Fill(etaTrack, phiTrack, histoWeight);
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
    
  }// charged particle loop

  // 2 perpendicular cones added, divide by 2 total amount of energy.
  perpConePtSumTrack /= 2;
  
  if ( fFillHistograms )
  {
    if ( fPartInCone == kNeutralAndCharged )
    {
      fhConeSumPtTrack ->Fill(ptTrig, coneptsumTrack    , histoWeight);
      fhConePtLeadTrack->Fill(ptTrig, coneptLeadTrack   , histoWeight);
    }
    
    // UE subtraction
    if ( fICMethod == kSumBkgSubIC )
    {
      fhPerpConeSumPt->Fill(ptTrig, perpConePtSumTrack, histoWeight);
      
      fhConeSumPtVSPerpCone->Fill(coneptsumTrack, perpConePtSumTrack, histoWeight);

      if ( fFillEtaPhiHistograms )
        fhPerpConeSumPtTrigEtaPhi->Fill(etaTrig, phiTrig, perpConePtSumTrack*histoWeight);
    }
    
    if ( fICMethod > kSumBkgSubIC )
    {
      fhConeSumPtEtaBandUETrack->Fill(ptTrig, etaBandPtSumTrack , histoWeight);
      fhConeSumPtPhiBandUETrack->Fill(ptTrig, phiBandPtSumTrack , histoWeight);
      
      fhConeSumPtVSUETracksEtaBand->Fill(coneptsumTrack, etaBandPtSumTrack, histoWeight);
      fhConeSumPtVSUETracksPhiBand->Fill(coneptsumTrack, phiBandPtSumTrack, histoWeight);
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, etaBandPtSumTrack *histoWeight); // check
        fhConeSumPtPhiBandUETrackTrigEtaPhi->Fill(etaTrig, phiTrig, phiBandPtSumTrack *histoWeight); // check
      }
    } // UE sub
  } // fill histograms

  // Add calculated values to pCandidate, might be modified later 
  // if UE subtraction and normalization requested
  
  // Add reference arrays to AOD when filling AODs only
  if ( bFillAOD && reftracks ) pCandidate->AddObjArray(reftracks);  
}

//_________________________________________________________________________________________________________________________________
/// Get normalization of cluster background band.
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandClusterNormalization( Float_t   etaC,                  Float_t   phiC,
                                                           Float_t   excessEta,             Float_t   excessPhi,         
                                                           Float_t   excessAreaEta,         Float_t   excessAreaPhi,         
                                                           Float_t   etaUEptsumCluster,     Float_t   phiUEptsumCluster,
                                                           Float_t & etaUEptsumClusterNorm, Float_t & phiUEptsumClusterNorm) const
{
  Float_t coneA = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  if ( fDistMinToTrigger > 0 ) coneA -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  Float_t fEMCPhiSize = fEMCPhiMax-fEMCPhiMin;
  
  /* //Catherine code
   if(((((2*fConeSize*fEMCPhiSize)-coneA))*phiBandBadCellsCoeff)!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*fConeSize*fEMCPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
   if(((((2*(fConeSize-excess)*fEMCPhiSize)-(coneA-excessFracEta))*etaBandBadCellsCoeff))!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA *coneBadCellsCoeff/ (((2*(fConeSize-excess)*fEMCPhiSize)-(coneA/excessFracEta))*etaBandBadCellsCoeff));
   if(((2*(fConeSize-excess)*fEMCEtaSize)-(coneA-excessFracPhi))*phiBandBadCellsCoeff!=0) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*(fConeSize-excess)*fEMCEtaSize)-(coneA/excessFracPhi))*phiBandBadCellsCoeff));
   */
  
  // UE band cans also be out of acceptance, need to estimate corrected area
  //
  if ( excessEta != 0 ) coneA /=  excessAreaEta;
  if ( excessPhi != 0 ) coneA /=  excessAreaPhi;
  if ( excessPhi != 0 ) coneA /=  excessAreaPhi;
  
  // Calculate normalization
  //    
  if(((2*fConeSize-excessEta)*fEMCPhiSize-coneA) != 0 ) phiUEptsumClusterNorm = phiUEptsumCluster*(coneA / ((((2*fConeSize-excessEta)*fEMCPhiSize)-coneA))); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
  if(((2*fConeSize-excessPhi)*fEMCEtaSize-coneA) != 0 ) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA / ((((2*fConeSize-excessPhi)*fEMCEtaSize)-coneA))); // pi * R^2 / (2 R * 2*0.7)  -  trigger cone
}

//________________________________________________________________________________________________________________________________
/// Get normalization of track background band.
//________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandTrackNormalization  (Float_t   etaC               ,  
                                                          Float_t   excessEta          ,          
                                                          Float_t   excessAreaEta      ,         
                                                          Float_t   etaUEptsumTrack    ,  Float_t   phiUEptsumTrack,
                                                          Float_t & etaUEptsumTrackNorm,  Float_t & phiUEptsumTrackNorm) const
{
  Float_t coneA = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  if ( fDistMinToTrigger > 0 ) coneA -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  /*//Catherine code
   //phiUEptsumTrackNorm = phiUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*tpcPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   //etaUEptsumTrackNorm = etaUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*fTPCEtaSize)-coneA))*etaBandBadCellsCoeff); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*fConeSize*tpcPhiSize-coneA)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*fConeSize*tpcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   if((2*fConeSize*fTPCEtaSize-coneA)!=0)etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / (((2*fConeSize*fTPCEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*(fConeSize-excess)*tpcPhiSize)-(coneA-excessFracEta)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*(fConeSize-excess)*tpcPhiSize)-(coneA/excessFracEta))));
   */ //end Catherine code

  // UE band cans also be out of acceptance, need to estimate corrected area
  //
  if ( excessEta != 0 ) coneA /=  excessAreaEta;
  
  // Calculate normalization
  //
  if(((2*fConeSize-excessEta)*fTPCPhiSize - coneA) !=0 ) phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / ((((2*fConeSize-excessEta)*fTPCPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
  if(( 2*fConeSize           *fTPCEtaSize - coneA) !=0 ) etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / ((( 2*fConeSize           *fTPCEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone

}

//________________________________________________________________________
/// If isolation cone are is outside a detector, calculate the area in excess.
/// \param excess: cone size minus acceptance of detector.
/// \return Area of a circunference segment 1/2 R^2 (angle-sin(angle)), angle = 2*ACos((R-excess)/R).
//________________________________________________________________________
Float_t AliIsolationCut::CalculateExcessAreaFraction(Float_t excess) const
{
  Float_t angle   = 2*TMath::ACos( (fConeSize-excess) / fConeSize );

  Float_t coneA   = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area

  Float_t excessA = fConeSize*fConeSize / 2 * (angle-TMath::Sin(angle));

  if ( coneA > excessA ) return coneA / (coneA-excessA);
  else
  {
    AliWarning(Form("Please Check : Excess %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f",
                    excess,coneA, excessA, angle*TMath::RadToDeg(), coneA / (coneA-excessA)));
    return  1;
  }
}

//________________________________________________________________________
/// If isolation cone are is outside a detector, calculate the area in excess.
/// \param etaC: Candidate pseudorapidity.
/// \param phiC: Candidate azimuthal angle
/// \param excessAreaTrkEta: cone excess out of TPC pseudo rapidity
/// \param excessAreaClsEta: cone excess out of EMC pseudo rapidity
/// \param excessAreaClsPhi: cone excess out of EMC azimuthal angle
//________________________________________________________________________
void AliIsolationCut::CalculateExcessAreaFractionForChargedAndNeutral
( Float_t etaC, Float_t phiC, 
  Float_t & excessTrkEta, Float_t & excessAreaTrkEta, 
  Float_t & excessClsEta, Float_t & excessAreaClsEta, 
  Float_t & excessClsPhi, Float_t & excessAreaClsPhi  ) const
{
  excessTrkEta = 0;
  excessClsEta = 0;
  excessClsPhi = 0;
  
  if ( TMath::Abs(etaC)+fConeSize > fTPCEtaSize/2. )
    excessTrkEta = TMath::Abs(etaC) + fConeSize - fTPCEtaSize/2.;
  
  if(TMath::Abs(etaC)+fConeSize > fEMCEtaSize/2.)
    excessClsEta = TMath::Abs(etaC) + fConeSize - fEMCEtaSize/2.;
    
  if     ( TMath::Abs(phiC)+fConeSize > fEMCPhiMax )
    excessClsPhi = (TMath::Abs(phiC) + fConeSize) - fEMCPhiMax;
  else if( TMath::Abs(phiC)-fConeSize < fEMCPhiMax )
    excessClsPhi = fEMCPhiMin - (TMath::Abs(phiC) - fConeSize) ;
  
  excessAreaTrkEta = CalculateExcessAreaFraction(excessTrkEta);
  excessAreaClsEta = CalculateExcessAreaFraction(excessClsEta);
  excessAreaClsPhi = CalculateExcessAreaFraction(excessClsPhi);
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
  
  if ( calorimeter == AliFiducialCut::kPHOS)
  {
    fEMCEtaSize = reader->GetFiducialCut()->GetPHOSFidCutMaxEtaArray()->At(0) -
                  reader->GetFiducialCut()->GetPHOSFidCutMinEtaArray()->At(0) ;
    fEMCPhiMin  = reader->GetFiducialCut()->GetPHOSFidCutMinPhiArray()->At(0) ;
    fEMCPhiMax  = reader->GetFiducialCut()->GetPHOSFidCutMaxPhiArray()->At(0) ;
  }
  
  // Info stored in degrees, put them on radians
  fEMCPhiMin *= TMath::DegToRad();
  fEMCPhiMax *= TMath::DegToRad();
  
  // Get the cut used for the TPC tracks in the reader, +-0.8, +-0.9 ...
  // Only valid in simple fidutial cut case and if the cut is applied, careful!
  fTPCEtaSize = reader->GetFiducialCut()->GetCTSFidCutMaxEtaArray()->At(0) -
                reader->GetFiducialCut()->GetCTSFidCutMinEtaArray()->At(0) ;
  fTPCPhiSize = TMath::Pi(); // Half TPC tracks with respect trigger candidate inspected

  AliDebug(1,Form("TPC: Deta %2.2f; Dphi %2.2f; Calo: Deta %2.2f, phi min %2.2f, max %2.2f",
                  fTPCEtaSize,fTPCPhiSize,fEMCEtaSize,
                  fEMCPhiMin*TMath::RadToDeg(),fEMCPhiMax*TMath::RadToDeg()));
  
//  printf("TPC: Deta %2.2f; Dphi %2.2f; Calo: Deta %2.2f, phi min %2.2f, max %2.2f\n",
//         fTPCEtaSize,fTPCPhiSize,fEMCEtaSize,
//         fEMCPhiMin*TMath::RadToDeg(),fEMCPhiMax*TMath::RadToDeg());
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

  Int_t   nptbins       = fHistoRanges->GetHistoPtBins();
  Int_t   nphibins      = fHistoRanges->GetHistoPhiBins();
  Int_t   netabins      = fHistoRanges->GetHistoEtaBins();
  Float_t ptmax         = fHistoRanges->GetHistoPtMax();
  Float_t phimax        = fHistoRanges->GetHistoPhiMax();
  Float_t etamax        = fHistoRanges->GetHistoEtaMax();
  Float_t ptmin         = fHistoRanges->GetHistoPtMin();
  Float_t phimin        = fHistoRanges->GetHistoPhiMin();
  Float_t etamin        = fHistoRanges->GetHistoEtaMin();
  
  Int_t   nptsumbins    = fHistoRanges->GetHistoNPtSumBins();
  Float_t ptsummax      = fHistoRanges->GetHistoPtSumMax();
  Float_t ptsummin      = fHistoRanges->GetHistoPtSumMin();
  Int_t   nptinconebins = fHistoRanges->GetHistoNPtInConeBins();
  Float_t ptinconemax   = fHistoRanges->GetHistoPtInConeMax();
  Float_t ptinconemin   = fHistoRanges->GetHistoPtInConeMin();
  
  // For UE subtracted histograms, shift it down by 20 GeV
  // keep same histogram binning.
  Float_t ptsumminUESub   = ptsummin-20;
  Float_t ptsummaxUESub   = ptsummax;
  Int_t   nptsumbinsUESub = nptsumbins*(1.+20./(ptsummax-ptsummin));
  
  TString sParticle = ", x^{0,#pm}";
  if      ( fPartInCone == kOnlyNeutral )  sParticle = ", x^{0}";
  else if ( fPartInCone == kOnlyCharged )  sParticle = ", x^{#pm}";
  
  TString parTitleR   = Form("#it{R} = %2.2f%s",fConeSize,sParticle.Data());
  
  if ( fPartInCone == kNeutralAndCharged )
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
  
  fhPtInCone  = new TH2F
  ("hPtInCone",
   Form("#it{p}_{T} of clusters and tracks in isolation cone for %s",parTitleR.Data()),
   nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
  fhPtInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
  fhPtInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtInCone) ;
  
  fhConePtLead  = new TH2F
  ("hConePtLead",
   Form("Track or Cluster  leading #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
   nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fhConePtLead->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
  fhConePtLead->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
  outputContainer->Add(fhConePtLead) ;
  
  fhConeSumPt  = new TH2F
  ("hConePtSum",
   Form("Track and Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",fConeSize),
   nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
  fhConeSumPt->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
  fhConeSumPt->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
  outputContainer->Add(fhConeSumPt) ;
  
  if ( fFillEtaPhiHistograms )
  {
    fhConeSumPtTrigEtaPhi  = new TH2F
    ("hConePtSumTrigEtaPhi",
     Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} in isolation cone for %s",parTitleR.Data()),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
    fhConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
    outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
  }
  
  if ( fICMethod >= kSumBkgSubIC )
  {
    fhConeSumPtUESub  = new TH2F
    ("hConePtSumUESub",
     Form("Track and/or Cluster #Sigma #it{p}_{T} in #it{R} = %2.2f, after UE correction",fConeSize),
     nptbins,ptmin,ptmax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
    fhConeSumPtUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtUESub->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtUESub) ;
    
    if ( fFillEtaPhiHistograms )
    {
      fhConeSumPtUESubTrigEtaPhi  = new TH2F
      ("hConePtSumUESubTrigEtaPhi",
       Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} for %s, after UE correction",parTitleR.Data()),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhConeSumPtUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtUESubTrigEtaPhi) ;
    }
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
  }
  
  // UE bands
  if ( fPartInCone != kOnlyCharged && fICMethod >= kSumBkgSubIC )
  {
    fhEtaBandClusterPt  = new TH2F
    ("hEtaBandClusterPt",
     Form("#it{p}_{T} of clusters in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhEtaBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhEtaBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
    outputContainer->Add(fhEtaBandClusterPt) ;
    
    fhPhiBandClusterPt  = new TH2F
    ("hPhiBandClusterPt",
     Form("#it{p}_{T} of clusters in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPhiBandClusterPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPhiBandClusterPt->SetYTitle("#it{p}_{T}^{cluster-band} (GeV/#it{c})");
    outputContainer->Add(fhPhiBandClusterPt) ;   
    
    if ( fFillEtaPhiHistograms )
    {
      fhEtaBandClusterEtaPhi  = new TH2F
      ("hEtaBandClusterEtaPhi",
       Form("#eta vs #varphi of clusters in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaBandClusterEtaPhi->SetXTitle("#eta");
      fhEtaBandClusterEtaPhi->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaBandClusterEtaPhi) ;
      
      fhPhiBandClusterEtaPhi  = new TH2F
      ("hPhiBandClusterEtaPhi",
       Form("#eta vs #varphi of clusters in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhPhiBandClusterEtaPhi->SetXTitle("#eta");
      fhPhiBandClusterEtaPhi->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhPhiBandClusterEtaPhi) ;
    }
    
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      fhConeSumPtEtaBandUECluster  = new TH2F
      ("hConePtSumEtaBandUECluster",
       "#Sigma cluster #it{p}_{T} in UE Eta Band",
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtEtaBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtEtaBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtEtaBandUECluster) ;
      
      fhConeSumPtPhiBandUECluster  = new TH2F
      ("hConePtSumPhiBandUECluster",
       "#Sigma cluster #it{p}_{T} UE Phi Band",
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtPhiBandUECluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtPhiBandUECluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtPhiBandUECluster) ;
      
      fhConeSumPtVSUEClusterEtaBand  = new TH2F
      ("hConeSumPtVSUEClusterEtaBand",
       Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #eta band for cluster (before normalization), R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
      fhConeSumPtVSUEClusterEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
      fhConeSumPtVSUEClusterEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSUEClusterEtaBand);
      
      fhConeSumPtVSUEClusterPhiBand  = new TH2F
      ("hConeSumPtVSUEClusterPhiBand",
       Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #varphi band for cluster (before normalization), R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
      fhConeSumPtVSUEClusterPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
      fhConeSumPtVSUEClusterPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSUEClusterPhiBand);    
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaBandUEClusterTrigEtaPhi  = new TH2F
        ("hConePtSumEtaBandUEClusterTrigEtaPhi",
         "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} in UE Eta Band",
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUEClusterTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUEClusterTrigEtaPhi  = new TH2F
        ("hConePtSumPhiBandUEClusterTrigEtaPhi",
         "Trigger #eta vs #varphi, #Sigma cluster #it{p}_{T} UE Phi Band",
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUEClusterTrigEtaPhi) ;
      }
      
      // Subtraction
      fhConeSumPtUEBandNormCluster  = new TH2F
      ("hConeSumPtUEBandNormCluster",
       Form("Clusters #Sigma #it{p}_{T} in normalized #eta or #varphi band, #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtUEBandNormCluster->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
      fhConeSumPtUEBandNormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUEBandNormCluster) ;
      
      fhConeSumPtClusterSubVsNoSub = new TH2F
      ("hConeSumPtClusterSubVsNoSub",
       Form("#Sigma #it{p}_{T} in cone before vs after UE bkg sub from #eta or #varphi band, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtClusterSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtClusterSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtClusterSubVsNoSub);
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
  }
  
  // UE subtraction
  if ( fPartInCone != kOnlyNeutral && fICMethod >= kSumBkgSubIC )
  {
    // UE in perpendicular cones
    //
    if ( fICMethod == kSumBkgSubIC )
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
      
      fhConeSumPtVSPerpCone = new TH2F
      ("hConeSumPtVSPerpCone",
       Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in 2 isolation cones at #pm 45 degree #varphi, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtVSPerpCone->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtVSPerpCone->SetYTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSPerpCone);
      
      if ( fFillEtaPhiHistograms )
      {
        fhEtaPhiInPerpCone = new TH2F
        ("hEtaPhiInPerpCone",
         Form("#eta vs #varphi of all Tracks in #perp cones"),
         netabins,-1,1,nphibins,0,TMath::TwoPi());
        fhEtaPhiInPerpCone->SetXTitle("#eta");
        fhEtaPhiInPerpCone->SetYTitle("#varphi (rad)");
        outputContainer->Add(fhEtaPhiInPerpCone) ;
        
        fhPerpConeSumPtTrigEtaPhi = new TH2F
        ("hPerpConeSumPtTrigEtaPhi",
         "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in #perp cones",
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhPerpConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}^{in #perp cone} (GeV/#it{c})");
        fhPerpConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhPerpConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhPerpConeSumPtTrigEtaPhi) ;
      }
    } // perpendicular
    
    // UE bands
    fhEtaBandTrackPt  = new TH2F
    ("hEtaBandTrackPt",
     Form("#it{p}_{T} of tracks in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhEtaBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhEtaBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
    outputContainer->Add(fhEtaBandTrackPt) ;
    
    fhPhiBandTrackPt  = new TH2F
    ("hPhiBandTrackPt",
     Form("#eta vs #varphi of tracks in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPhiBandTrackPt->SetXTitle("#it{p}_{T}^{trig} (GeV/#it{c})");
    fhPhiBandTrackPt->SetYTitle("#it{p}_{T}^{track-band} (GeV/#it{c})");
    outputContainer->Add(fhPhiBandTrackPt) ;
    
    if ( fFillEtaPhiHistograms )
    {
      fhEtaBandTrackEtaPhi  = new TH2F
      ("hEtaBandTrackEtaPhi",
       Form("#eta vs #varphi of tracks in #eta band isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaBandTrackEtaPhi->SetXTitle("#eta");
      fhEtaBandTrackEtaPhi->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaBandTrackEtaPhi) ;
      
      fhPhiBandTrackEtaPhi  = new TH2F
      ("hPhiBandTrackEtaPhi",
       Form("#eta vs #varphi of tracks in #varphi band isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhPhiBandTrackEtaPhi->SetXTitle("#eta");
      fhPhiBandTrackEtaPhi->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhPhiBandTrackEtaPhi) ;
    }
    
    if ( fICMethod >= kSumBkgSubEtaBandIC )
    {
      fhConeSumPtEtaBandUETrack  = new TH2F
      ("hConePtSumEtaBandUETrack",
       "#Sigma track #it{p}_{T} in UE Eta Band",
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtEtaBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtEtaBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtEtaBandUETrack) ;
      
      fhConeSumPtPhiBandUETrack  = new TH2F
      ("hConePtSumPhiBandUETrack",
       "#Sigma track #it{p}_{T} in UE Phi Band",
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtPhiBandUETrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtPhiBandUETrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtPhiBandUETrack) ;
      
      fhConeSumPtVSUETracksEtaBand  = new TH2F
      ("hConeSumPtVSUETracksEtaBand",
       Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #eta band for tracks (before normalization), R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,2*nptsumbins,ptsummin,2*ptsummax);
      fhConeSumPtVSUETracksEtaBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
      fhConeSumPtVSUETracksEtaBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSUETracksEtaBand);
      
      fhConeSumPtVSUETracksPhiBand  = new TH2F
      ("hConeSumPtVSUETracksPhiBand",
       Form("#Sigma #it{p}_{T} in cone versus #Sigma #it{p}_{T} in #varphi band for tracks (before normalization), R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,8*nptsumbins,ptsummin,8*ptsummax);
      fhConeSumPtVSUETracksPhiBand->SetXTitle("#Sigma #it{p}_{T} cone (GeV/#it{c})");
      fhConeSumPtVSUETracksPhiBand->SetYTitle("#Sigma #it{p}_{T} UE (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtVSUETracksPhiBand);
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaBandUETrackTrigEtaPhi  = new TH2F
        ("hConePtSumEtaBandUETrackTrigEtaPhi",
         "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE Eta Band",
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaBandUETrackTrigEtaPhi) ;
        
        fhConeSumPtPhiBandUETrackTrigEtaPhi  = new TH2F
        ("hConePtSumPhiBandUETrackTrigEtaPhi",
         "Trigger #eta vs #varphi, #Sigma track #it{p}_{T} in UE Phi Band",
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiBandUETrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiBandUETrackTrigEtaPhi) ;
      }
      
      // Subtraction
      fhConeSumPtUEBandNormTrack  = new TH2F
      ("hConeSumPtUEBandNormTrack",
       Form("Tracks #Sigma #it{p}_{T} in normalized #eta or #varphi band, #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
      fhConeSumPtUEBandNormTrack->SetYTitle("#Sigma #it{p}_{T}^{band-norm} (GeV/#it{c})");
      fhConeSumPtUEBandNormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtUEBandNormTrack) ;
      
      fhConeSumPtTrackSubVsNoSub = new TH2F
      ("hConeSumPtTrackSubVsNoSub",
       Form("#Sigma #it{p}_{T} in cone before and after UE bkg subtraction from #eta or #varphi band, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,nptsumbinsUESub,ptsumminUESub,ptsummaxUESub);
      fhConeSumPtTrackSubVsNoSub->SetXTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtTrackSubVsNoSub->SetYTitle("#Sigma #it{p}_{T, UE sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtTrackSubVsNoSub);
    }
  } // charged
    
  // UE subtraction
  
  if ( fPartInCone == kNeutralAndCharged && fICMethod >= kSumBkgSubEtaBandIC )
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
  
  if ( fFillHistograms)
  {
    if( fMakeConeExcessCorr || fICMethod >= kSumBkgSubEtaBandIC )
    {
      fhFractionClusterOutConeEta  = new TH2F
      ("hFractionClusterOutConeEta",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionClusterOutConeEta->SetYTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConeEta) ;
      
      fhFractionClusterOutConeEtaTrigEtaPhi  = new TH2F
      ("hFractionClusterOutConeEtaTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;
      
      fhFractionClusterOutConePhi  = new TH2F
      ("hFractionClusterOutConePhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #varphi acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionClusterOutConePhi->SetYTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConePhi) ;
      
      fhFractionClusterOutConePhiTrigEtaPhi  = new TH2F
      ("hFractionClusterOutConePhiTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #varphi acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;
      
      fhFractionClusterOutConeEtaPhi  = new TH2F
      ("hFractionClusterOutConeEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta x #varphi acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionClusterOutConeEtaPhi->SetYTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaPhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConeEtaPhi) ;
      
      fhFractionClusterOutConeEtaPhiTrigEtaPhi  = new TH2F
      ("hFractionClusterOutConeEtaPhiTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of clusters #eta  x #varphi acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetZTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConeEtaPhiTrigEtaPhi) ;
      
      fhFractionTrackOutConeEta  = new TH2F
      ("hFractionTrackOutConeEta",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of tracks #eta acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionTrackOutConeEta->SetYTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionTrackOutConeEta) ;
      
      fhFractionTrackOutConeEtaTrigEtaPhi  = new TH2F
      ("hFractionTrackOutConeEtaTrigEtaPhi",
       Form("Fraction of cone area (#it{R} =  %2.2f), out of tracks #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("1-(#it{A}_{cone}+#it{A}_{excess})/#it{A}_{cone}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
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
  snprintf(onePar,buffersize,"fConeSize=%1.2f;",fConeSize) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtThreshold>%2.2f;<%2.2f;",fPtThreshold,fPtThresholdMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fSumPtThreshold<%2.2f, gap = %2.2f;<%2.2f;",fSumPtThreshold, fSumPtThresholdGap, fSumPtThresholdMax) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtFraction=%2.2f;",fPtFraction) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fICMethod=%d;",fICMethod) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPartInCone=%d;",fPartInCone) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fFracIsThresh=%i;",fFracIsThresh) ;
  parList+=onePar ;
  snprintf(onePar,buffersize," fIsTMClusterInConeRejected=%i;", fIsTMClusterInConeRejected) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fDistMinToTrigger=%1.2f;",fDistMinToTrigger) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fFillHistograms=%d,fFillEtaPhiHistograms=%d;",fFillHistograms,fFillEtaPhiHistograms) ;
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
  fConeSize             = 0.4 ;
  fPtThreshold          = 0.5  ;
  fPtThresholdMax       = 10000.  ;
  fSumPtThreshold       = 2.0 ;
  fSumPtThresholdGap    = 0.5 ;
  fSumPtThresholdMax    = 10000. ;
  fPtFraction           = 0.1 ;
  fPartInCone           = kNeutralAndCharged;
  fICMethod             = kSumPtIC; // 0 pt threshol method, 1 cone pt sum method
  fFracIsThresh         = 1;
  fDistMinToTrigger     = -1.; // no effect
  fNeutralOverChargedRatio = 0.363; // Based on pPb analysis, to be confirmed on other systems. 
                                    // Use eta band for charged and neutrals for estimation.
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
 Bool_t  & isolated   , Double_t histoWeight   
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
  Float_t coneptLeadTrack     = 0;
  
  Float_t etaBandPtSumTrack   = 0;
  Float_t phiBandPtSumTrack   = 0;
  Float_t perpPtSumTrack      = 0;
  Float_t etaBandPtSumCluster = 0;
  Float_t phiBandPtSumCluster = 0;
  
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
  
  //printf("Get track signal\n");
  CalculateTrackSignalInCone   (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgTrk,
                                nPart              , nfrac,
                                coneptsumTrack     , coneptLeadTrack,
                                etaBandPtSumTrack  , phiBandPtSumTrack,
                                perpPtSumTrack     , histoWeight);
  
  //printf("Get calo signal\n");
  CalculateCaloSignalInCone    (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgCls,
                                calorimeter        , pid, 
                                nPart              , nfrac,
                                coneptsumCluster   , coneptLeadCluster,
                                etaBandPtSumCluster, phiBandPtSumCluster,
                                histoWeight);
  
  // Add leading found information to candidate object
  pCandidate->SetNeutralLeadPtInCone(coneptLeadCluster);
  pCandidate->SetChargedLeadPtInCone(coneptLeadTrack);
  
  // Calculate how much of the cone got out the detectors acceptance
  //
  Float_t excessAreaTrkEta = 1;
  Float_t excessAreaClsEta = 1;
  Float_t excessAreaClsPhi = 1;
  Float_t excessTrkEta     = 0;
  Float_t excessClsEta     = 0;
  Float_t excessClsPhi     = 0;
  
  if( fMakeConeExcessCorr || fICMethod >= kSumBkgSubEtaBandIC )
  {
    GetDetectorAngleLimits(reader,calorimeter); // Do it once, check done inside
    
    CalculateExcessAreaFractionForChargedAndNeutral(etaC, phiC, 
                                                    excessTrkEta, excessAreaTrkEta, 
                                                    excessClsEta, excessAreaClsEta, 
                                                    excessClsPhi, excessAreaClsPhi);
    
    // Store excess areas
    pCandidate->SetNeutralConeExcessAreaEta(excessAreaClsEta);
    pCandidate->SetNeutralConeExcessAreaPhi(excessAreaClsPhi);
    pCandidate->SetChargedConeExcessAreaEta(excessAreaTrkEta);
    pCandidate->SetChargedConeExcessAreaPhi(1);
    
    if ( fFillHistograms )
    {                   
      fhFractionTrackOutConeEta            ->Fill(ptC ,        excessAreaTrkEta-1, histoWeight);
      fhFractionTrackOutConeEtaTrigEtaPhi  ->Fill(etaC, phiC, (excessAreaTrkEta-1)*histoWeight); // check
      fhFractionClusterOutConeEta          ->Fill(ptC ,        excessAreaClsEta-1, histoWeight);
      fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaC, phiC, (excessAreaClsEta-1)*histoWeight); // check
      fhFractionClusterOutConePhi          ->Fill(ptC ,        excessAreaClsPhi-1, histoWeight);
      fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaC, phiC, (excessAreaClsPhi-1)*histoWeight); // check
      fhFractionClusterOutConeEtaPhi          ->Fill(ptC ,        excessAreaClsPhi*excessAreaClsEta-1, histoWeight);
      fhFractionClusterOutConeEtaPhiTrigEtaPhi->Fill(etaC, phiC, (excessAreaClsPhi*excessAreaClsEta-1)*histoWeight); // check
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
  Double_t coneptsumUESub        = coneptsum;
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
    
    if ( fICMethod == kSumBkgSubIC )
    {
      coneptsumBkgTrk = perpPtSumTrack;
      coneptsumBkgCls = perpPtSumTrack*fNeutralOverChargedRatio; 
      //coneptsumBkg    = perpPtSumTrack*(1+fNeutralOverChargedRatio);
      
      // Add to candidate object
      pCandidate->SetChargedPtSumInPerpCone(perpPtSumTrack*excessAreaTrkEta);
    } // UE subtraction by perpendicular cones
    else if ( fICMethod >= kSumBkgSubEtaBandIC ) // eta or phi band
    {
      //printf("UE band\n");
      Float_t  etaBandPtSumTrackNorm   = 0;
      Float_t  phiBandPtSumTrackNorm   = 0;
      Float_t  etaBandPtSumClusterNorm = 0;
      Float_t  phiBandPtSumClusterNorm = 0;
      
      Float_t  coneptsumTrackSub   = 0 ;
      Float_t  coneptsumClusterSub = 0 ;
  
      // Normalize background to cone area
      if     ( fPartInCone != kOnlyCharged )
      {
        CalculateUEBandClusterNormalization(etaC                   , phiC                  ,
                                            excessClsEta           , excessClsPhi          ,
                                            excessAreaClsEta       , excessAreaClsPhi      ,
                                            etaBandPtSumCluster    , phiBandPtSumCluster   ,
                                            etaBandPtSumClusterNorm, phiBandPtSumClusterNorm);
        
        if      ( fICMethod == kSumBkgSubEtaBandIC )
        {
          coneptsumBkgClsRaw  = etaBandPtSumCluster; 
          coneptsumBkgCls     = etaBandPtSumClusterNorm; 
          coneptsumClusterSub = coneptsumCluster - coneptsumBkgCls;

        }
        else  if( fICMethod == kSumBkgSubPhiBandIC )
        {
          coneptsumBkgClsRaw  = phiBandPtSumCluster; 
          coneptsumBkgCls     = phiBandPtSumClusterNorm;
          coneptsumClusterSub = coneptsumCluster - coneptsumBkgCls;
        }
        
        
        //      printf("Cluster: sumpT %2.2f, phi: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
        //             "\t eta: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
        //             "\t subtracted:  %2.2f\n",
        //             coneptsumCluster, phiBandPtSumCluster, coneptsumBkgCls, excessFracPhiCluster, 
        //             coneptsumBkgClsRaw, coneptsumBkgCls, excessFracEtaCluster,
        //             coneptsumClusterSubEta);
        
        if ( fFillHistograms )
        {
          fhConeSumPtUEBandNormCluster->Fill(ptC, coneptsumBkgCls, histoWeight);
          fhConeSumPtClusterSubVsNoSub->Fill(coneptsumCluster, coneptsumClusterSub, histoWeight);
        } // histograms
      } // clusters in cone
      
      //printf("Pass cluster\n");

      if     ( fPartInCone != kOnlyNeutral )
      {
        CalculateUEBandTrackNormalization(etaC                 ,
                                          excessTrkEta         ,    
                                          excessAreaTrkEta     ,    
                                          etaBandPtSumTrack    , phiBandPtSumTrack,
                                          etaBandPtSumTrackNorm, phiBandPtSumTrackNorm);
        
        if      ( fICMethod == kSumBkgSubEtaBandIC )
        {
          coneptsumBkgTrkRaw = etaBandPtSumTrack; 
          coneptsumBkgTrk    = etaBandPtSumTrackNorm; 
          coneptsumTrackSub  = coneptsumTrack - coneptsumBkgTrk;
        }
        else  if( fICMethod == kSumBkgSubPhiBandIC )
        {
          coneptsumBkgTrkRaw = phiBandPtSumTrack; 
          coneptsumBkgTrk    = phiBandPtSumTrackNorm; 
          coneptsumTrackSub  = coneptsumTrack - coneptsumBkgTrk;
        }
        
        //      printf("Track: sumpT %2.2f, phi: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
        //             "\t eta: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
        //             "\t subtracted: %2.2f\n",
        //             coneptsumTrack, phiBandPtSumTrack, coneptsumBkgTrk, excessFracPhiTrack, 
        //             etaBandPtSumTrack, coneptsumBkgTrk, excessFracEtaTrack,
        //             coneptsumTrackSub);  
        
        if ( fFillHistograms )
        {          
          fhConeSumPtUEBandNormTrack->Fill(ptC, coneptsumBkgTrk, histoWeight);
          fhConeSumPtTrackSubVsNoSub->Fill(coneptsumTrack, coneptsumTrackSub, histoWeight);
        } // fill 
      } // tracks in cone
      
      //printf("Pass track\n");

      if ( fPartInCone == AliIsolationCut::kNeutralAndCharged && 
           fFillHistograms )
      {
        fhBandClustervsTrack    ->Fill(coneptsumBkgClsRaw, coneptsumBkgTrkRaw, histoWeight);
        fhBandNormClustervsTrack->Fill(coneptsumBkgCls   , coneptsumBkgTrk   , histoWeight);
        
        fhConeSumPtUEBandSubClustervsTrack->Fill(coneptsumClusterSub, coneptsumTrackSub, histoWeight);
      }
        
      //printf("Pass both\n");

      // Add to candidate object
      pCandidate->SetNeutralPtSumEtaBand(coneptsumBkgCls);
      pCandidate->SetNeutralPtSumPhiBand(coneptsumBkgCls);
      
      pCandidate->SetChargedPtSumEtaBand(coneptsumBkgTrk);
      pCandidate->SetChargedPtSumPhiBand(coneptsumBkgTrk); 
      
    } // UE subtraction by Eta band
    
    coneptsumUESubCluster -= coneptsumBkgCls;
    coneptsumUESubTrack   -= coneptsumBkgTrk;
    
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
  
  //-------------------------------------------------------------------
  // Fill histograms
  //-------------------------------------------------------------------

  if ( !fFillHistograms ) return;
  
  // Sum pt in cone
  //
  //printf("fill histo");
  fhConeSumPt->Fill(ptC, coneptsum, histoWeight);
  
  if ( fFillEtaPhiHistograms )
    fhConeSumPtTrigEtaPhi->Fill(etaC, phiC, coneptsum*histoWeight); // check
  
  if ( fPartInCone == kNeutralAndCharged && fICMethod != kSumBkgSubIC ) // No need for perpendicular or charged/neutral only analysis
  {
    fhConeSumPtClustervsTrack ->Fill(coneptsumCluster, coneptsumTrack, histoWeight);
    if(coneptsumTrack > 0) fhConeSumPtClusterTrackFrac ->Fill(ptC, coneptsumCluster /coneptsumTrack, histoWeight);
  }
  
  // Here the sum in cone before subtraction, if done, to check the effect.
  //
  if ( fICMethod >= kSumBkgSubIC )
  {
    fhConeSumPtUESub ->Fill(ptC, coneptsumUESub, histoWeight);
    
    if ( fFillEtaPhiHistograms )
      fhConeSumPtUESubTrigEtaPhi->Fill(etaC, phiC, coneptsumUESub*histoWeight); // check
    
    // No for charged/neutral only analysis
    if ( fPartInCone == kNeutralAndCharged )
    {
      fhConeSumPtUESubTrack  ->Fill(ptC, coneptsumUESubTrack  , histoWeight);
      fhConeSumPtUESubCluster->Fill(ptC, coneptsumUESubCluster, histoWeight);

      // No need for perpendicular cones
      if( fICMethod != kSumBkgSubIC )
      {
        
        fhConeSumPtUESubClustervsTrack ->Fill(coneptsumUESubCluster, coneptsumUESubTrack, histoWeight);
        if ( TMath::Abs(coneptsumUESubTrack) > 0 ) 
          fhConeSumPtUESubClusterTrackFrac ->Fill(ptC, coneptsumUESubCluster /coneptsumUESubTrack, histoWeight);
      }
    }
  }
  
  // Leading in cone
  //
  Float_t coneptLead = coneptLeadTrack;
  if(coneptLeadCluster > coneptLeadTrack) 
    coneptLead = coneptLeadCluster;
  
  if ( fPartInCone == kNeutralAndCharged )
  {
    fhConePtLeadClustervsTrack->Fill(coneptLeadCluster,coneptLeadTrack, histoWeight);
   
    if ( coneptLeadTrack > 0 ) 
      fhConePtLeadClusterTrackFrac->Fill(ptC, coneptLeadCluster/coneptLeadTrack, histoWeight);
  }
  
  fhConePtLead->Fill(ptC, coneptLead, histoWeight);
  
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
  printf("Cone Size          =     %1.2f\n", fConeSize   ) ;
  printf("pT threshold       =     >%2.1f;<%2.1f\n", fPtThreshold, fPtThresholdMax) ;
  printf("Sum pT threshold   =     >%2.1f;<%2.1f, gap = %2.1f\n", fSumPtThreshold, fSumPtThresholdMax, fSumPtThresholdGap) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction ) ;
  printf("particle type in cone =  %d\n",    fPartInCone ) ;
  printf("using fraction for high pt leading instead of frac ? %i\n",fFracIsThresh);
  printf("minimum distance to candidate, R>%1.2f\n",fDistMinToTrigger);
  printf("correct cone excess = %d \n",fMakeConeExcessCorr);
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



