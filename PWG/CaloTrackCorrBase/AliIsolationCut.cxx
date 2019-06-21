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
fFillHistograms(0),  fFillEtaPhiHistograms(0),
fConeSize(0.),       fPtThreshold(0.),              fPtThresholdMax(10000.),
fSumPtThreshold(0.), fSumPtThresholdMax(10000.),    fPtFraction(0.),     
fICMethod(0),        fPartInCone(0),
fFracIsThresh(1),    fIsTMClusterInConeRejected(1), fDistMinToTrigger(-1.),
fDebug(0),           fMomentum(),                   fTrackVector(),
// Histograms
fHistoRanges(0),
fhPtInCone(0),       fhPtClusterInCone(0),          
fhPtTrackInCone(0),  fhPtInPerpCone(0),

fhConeSumPt(0),      fhConeSumPtCluster(0),         
fhConeSumPtTrack(0), fhPerpConeSumPt(0),
fhConePtLead(0),     fhConePtLeadCluster(0),        fhConePtLeadTrack(0),

fhConeSumPtClustervsTrack(0),     fhConeSumPtClusterTrackFrac(0),
fhConePtLeadClustervsTrack(0),    fhConePtLeadClusterTrackFrac(0),

fhConeSumPtTrigEtaPhi(0),

fhEtaPhiCluster(0),       fhEtaPhiTrack(0),
fhEtaPhiInConeCluster(0), fhEtaPhiInConeTrack(0),   fhEtaPhiInPerpCone(0),

fhEtaBandClusterPt(0), fhPhiBandClusterPt(0),
fhEtaBandTrackPt(0),   fhPhiBandTrackPt(0),

fhConeSumPtEtaBandUECluster(0),             fhConeSumPtPhiBandUECluster(0),
fhConeSumPtEtaBandUETrack(0),               fhConeSumPtPhiBandUETrack(0),

fhEtaBandClusterEtaPhi(0),        fhPhiBandClusterEtaPhi(0),
fhEtaBandTrackEtaPhi(0),          fhPhiBandTrackEtaPhi(0),

fhConeSumPtEtaBandUEClusterTrigEtaPhi(0),   fhConeSumPtPhiBandUEClusterTrigEtaPhi(0),
fhConeSumPtEtaBandUETrackTrigEtaPhi(0),     fhConeSumPtPhiBandUETrackTrigEtaPhi(0),
fhConeSumPtVSUETracksEtaBand(0),            fhConeSumPtVSUETracksPhiBand(0),
fhConeSumPtVSUEClusterEtaBand(0),           fhConeSumPtVSUEClusterPhiBand(0),

fhConeSumPtEtaUESub(0),                     fhConeSumPtPhiUESub(0),
fhConeSumPtEtaUESubTrigEtaPhi(0),           fhConeSumPtPhiUESubTrigEtaPhi(0),
fhConeSumPtEtaUENormCluster(0),             fhConeSumPtPhiUENormCluster(0),
fhConeSumPtEtaUESubCluster(0),              fhConeSumPtPhiUESubCluster(0),
fhConeSumPtEtaUESubClusterTrigEtaPhi(0),    fhConeSumPtPhiUESubClusterTrigEtaPhi(0),
fhConeSumPtEtaUENormTrack(0),               fhConeSumPtPhiUENormTrack(0),
fhConeSumPtEtaUESubTrack(0),                fhConeSumPtPhiUESubTrack(0),
fhConeSumPtEtaUESubTrackTrigEtaPhi(0),      fhConeSumPtPhiUESubTrackTrigEtaPhi(0),
fhConeSumPtEtaUESubClustervsTrack(0),       fhConeSumPtPhiUESubClustervsTrack(0),
fhEtaBandClustervsTrack(0),                 fhPhiBandClustervsTrack(0),
fhEtaBandNormClustervsTrack(0),             fhPhiBandNormClustervsTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiTrack(0),    fhConeSumPtSubNormvsConeSumPtTotPhiTrack(0),
fhConeSumPtSubvsConeSumPtTotEtaTrack(0),    fhConeSumPtSubNormvsConeSumPtTotEtaTrack(0),
fhConeSumPtSubvsConeSumPtTotPhiCluster(0),  fhConeSumPtSubNormvsConeSumPtTotPhiCluster(0),
fhConeSumPtSubvsConeSumPtTotEtaCluster(0),  fhConeSumPtSubNormvsConeSumPtTotEtaCluster(0)
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
    fhConeSumPtCluster ->Fill(pCandidate->Pt(), 0., histoWeight);
    
    fhConePtLeadCluster->Fill(pCandidate->Pt(), 0., histoWeight);
  
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
  
  if ( fFillHistograms )
  {
    if( fPartInCone == kNeutralAndCharged )
    {
      fhConeSumPtCluster ->Fill(ptC, coneptsumCluster , histoWeight);
      fhConePtLeadCluster->Fill(ptC, coneptLeadCluster, histoWeight);
    }  
    
    // UE substraction
    if ( fICMethod >= kSumBkgSubIC )
    {
      fhConeSumPtEtaBandUECluster->Fill(ptC, etaBandPtSumCluster, histoWeight);
      fhConeSumPtPhiBandUECluster->Fill(ptC, phiBandPtSumCluster, histoWeight);
      
      fhConeSumPtVSUEClusterEtaBand->Fill(coneptsumCluster, etaBandPtSumCluster, histoWeight);
      fhConeSumPtVSUEClusterPhiBand->Fill(coneptsumCluster, phiBandPtSumCluster, histoWeight);
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaBandUEClusterTrigEtaPhi->Fill(etaC, phiC, etaBandPtSumCluster*histoWeight); // Check
        fhConeSumPtPhiBandUEClusterTrigEtaPhi->Fill(etaC, phiC, phiBandPtSumCluster*histoWeight); // Check
      }
    } // UE sub
  } // fill histo
  
  // Add calculated values to pCandidate, might be modified later 
  // if UE subtraction and normalization requested

  pCandidate->SetNeutralLeadPtInCone(coneptLeadCluster);
  pCandidate->SetNeutralPtSumInCone (coneptsumCluster);
  
  pCandidate->SetNeutralPtSumEtaBand(etaBandPtSumCluster);
  pCandidate->SetNeutralPtSumPhiBand(phiBandPtSumCluster);
  
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
/// \param perpConePtSumTrack: sum of tracks in perpendicular cones in phi
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
      if ( etaTrack > (etaTrig-fConeSize) && etaTrack < (etaTrig+fConeSize) ) 
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
    if ( fICMethod >= kSumBkgSubIC )
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
        perpConePtSumTrack+=track->Pt();
        
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
    if ( fICMethod >= kSumBkgSubIC )
    {
      fhPerpConeSumPt->Fill(ptTrig, perpConePtSumTrack, histoWeight);
      
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

  pCandidate->SetChargedLeadPtInCone(coneptLeadTrack);
  pCandidate->SetChargedPtSumInCone (coneptsumTrack);
  
  pCandidate->SetChargedPtSumInPerpCone(perpConePtSumTrack);
  
  pCandidate->SetChargedPtSumEtaBand(etaBandPtSumTrack);
  pCandidate->SetChargedPtSumPhiBand(phiBandPtSumTrack);
  
  // Add reference arrays to AOD when filling AODs only
  if ( bFillAOD && reftracks ) pCandidate->AddObjArray(reftracks);  
}

//_________________________________________________________________________________________________________________________________
/// Get normalization of cluster background band.
//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandClusterNormalization(AliCaloTrackReader * /*reader*/, Float_t   etaC, Float_t /*phiC*/,
                                                           Float_t   phiUEptsumCluster,     Float_t   etaUEptsumCluster,
                                                           Float_t & phiUEptsumClusterNorm, Float_t & etaUEptsumClusterNorm,
                                                           Float_t & excessFracPhi,         Float_t & excessFracEta         ) const
{
  Float_t coneA = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  if ( fDistMinToTrigger > 0 ) coneA -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();

  //Careful here if EMCal limits changed .. 2010 (4 SM) to 2011-12 (10 SM), for the moment consider 100 deg in phi
  Float_t emcEtaSize = 0.7*2; // TO FIX
  Float_t emcPhiSize = TMath::DegToRad()*100.; // TO FIX

  /* //Catherine code
   if(((((2*fConeSize*emcPhiSize)-coneA))*phiBandBadCellsCoeff)!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*fConeSize*emcPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
   if(((((2*(fConeSize-excess)*emcPhiSize)-(coneA-excessFracEta))*etaBandBadCellsCoeff))!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA *coneBadCellsCoeff/ (((2*(fConeSize-excess)*emcPhiSize)-(coneA/excessFracEta))*etaBandBadCellsCoeff));
   if(((2*(fConeSize-excess)*emcEtaSize)-(coneA-excessFracPhi))*phiBandBadCellsCoeff!=0) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*(fConeSize-excess)*emcEtaSize)-(coneA/excessFracPhi))*phiBandBadCellsCoeff));
   */

  if((2*fConeSize*emcPhiSize-coneA)!=0) phiUEptsumClusterNorm = phiUEptsumCluster*(coneA / (((2*fConeSize*emcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
  if((2*fConeSize*emcEtaSize-coneA)!=0) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA / (((2*fConeSize*emcEtaSize)-coneA))); // pi * R^2 / (2 R * 2*0.7)  -  trigger cone

  //out of eta acceptance
  excessFracEta = 1;
  excessFracPhi = 1;

  if(TMath::Abs(etaC)+fConeSize > emcEtaSize/2.)
  {
    Float_t excess = TMath::Abs(etaC) + fConeSize - emcEtaSize/2.;
    excessFracEta  = CalculateExcessAreaFraction(excess);

    if ( excessFracEta != 0) coneA /=  excessFracEta;

    //UE band is also out of acceptance, need to estimate corrected area
    if(((2*fConeSize-excess)*emcPhiSize-coneA) != 0 ) phiUEptsumClusterNorm = phiUEptsumCluster*(coneA / ((((2*fConeSize-excess)*emcPhiSize)-coneA)));
    if(( 2*fConeSize        *emcEtaSize-coneA) != 0 ) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA / ((( 2*fConeSize        *emcEtaSize)-coneA)));
  }
}

//________________________________________________________________________________________________________________________________
/// Get normalization of track background band.
//________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandTrackNormalization  (AliCaloTrackReader * reader,    Float_t   etaC, Float_t /*phiC*/,
                                                          Float_t   phiUEptsumTrack,      Float_t   etaUEptsumTrack,
                                                          Float_t & phiUEptsumTrackNorm,  Float_t & etaUEptsumTrackNorm,
                                                          Float_t & excessFracPhi,        Float_t & excessFracEta         ) const
{
  Float_t coneA = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  if ( fDistMinToTrigger > 0 ) coneA -= fDistMinToTrigger*fDistMinToTrigger*TMath::Pi();
  
  // Get the cut used for the TPC tracks in the reader, +-0.8, +-0.9 ...
  // Only valid in simple fidutial cut case and if the cut is applied, careful!
  Float_t tpcEtaSize = reader->GetFiducialCut()->GetCTSFidCutMaxEtaArray()->At(0) -
  reader->GetFiducialCut()->GetCTSFidCutMinEtaArray()->At(0) ;
  Float_t tpcPhiSize = TMath::TwoPi();

  /*//Catherine code
   //phiUEptsumTrackNorm = phiUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*tpcPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   //etaUEptsumTrackNorm = etaUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*tpcEtaSize)-coneA))*etaBandBadCellsCoeff); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*fConeSize*tpcPhiSize-coneA)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*fConeSize*tpcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   if((2*fConeSize*tpcEtaSize-coneA)!=0)etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / (((2*fConeSize*tpcEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*(fConeSize-excess)*tpcPhiSize)-(coneA-excessFracEta)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*(fConeSize-excess)*tpcPhiSize)-(coneA/excessFracEta))));
   */ //end Catherine code

  //correct out of eta acceptance
  excessFracEta = 1;
  excessFracPhi = 1;

  if((2*fConeSize*tpcPhiSize-coneA)!=0) phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*fConeSize*tpcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
  if((2*fConeSize*tpcEtaSize-coneA)!=0) etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / (((2*fConeSize*tpcEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone

  if(TMath::Abs(etaC)+fConeSize > tpcEtaSize/2.)
  {
    Float_t excess = TMath::Abs(etaC) + fConeSize - tpcEtaSize/2.;
    excessFracEta  = CalculateExcessAreaFraction(excess);
    if (excessFracEta != 0) coneA /=  excessFracEta;

    //UE band is also out of acceptance, need to estimate corrected area
    if(((2*fConeSize-excess)*tpcPhiSize - coneA) !=0 ) phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / ((((2*fConeSize-excess)*tpcPhiSize)-coneA)));
    if(( 2*fConeSize        *tpcEtaSize - coneA) !=0 ) etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / ((( 2*fConeSize        *tpcEtaSize)-coneA)));
  }
  
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

  if(coneA > excessA) return coneA / (coneA-excessA);
  else
  {
    AliWarning(Form("Please Check : Excess Track %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f",
                    excess,coneA, excessA, angle*TMath::RadToDeg(), coneA / (coneA-excessA)));
    return  1;
  }
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
  
  TString sParticle = ", x^{0,#pm}";
  if      ( fPartInCone == kOnlyNeutral )  sParticle = ", x^{0}";
  else if ( fPartInCone == kOnlyCharged )  sParticle = ", x^{#pm}";
  
  TString parTitleR   = Form("#it{R} = %2.2f%s",fConeSize,sParticle.Data());

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
    
  fhConeSumPtTrigEtaPhi  = new TH2F
  ("hConePtSumTrigEtaPhi",
   Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} in isolation cone for %s",parTitleR.Data()),
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhConeSumPtTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
  fhConeSumPtTrigEtaPhi->SetXTitle("#eta_{trigger}");
  fhConeSumPtTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
  outputContainer->Add(fhConeSumPtTrigEtaPhi) ;
  
  if ( fPartInCone == kNeutralAndCharged )
  {
    fhPtClusterInCone  = new TH2F
    ("hPtClusterInCone",
     Form("#it{p}_{T} of clusters in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtClusterInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtClusterInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtClusterInCone) ;
    
    fhConePtLeadCluster  = new TH2F
    ("hConeLeadPtCluster",
     Form("Cluster leading in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhConePtLeadCluster->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
    fhConePtLeadCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConePtLeadCluster) ;
    
    fhConeSumPtCluster  = new TH2F
    ("hConePtSumCluster",
     Form("Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} = %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPtCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtCluster->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtCluster) ;
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
    fhConeSumPtEtaUENormCluster  = new TH2F
    ("hConeSumPtEtaUENormCluster",
     Form("Clusters #Sigma #it{p}_{T} in normalized #eta band, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtEtaUENormCluster->SetYTitle("#Sigma #it{p}_{T}^{#eta-band}_{norm} (GeV/#it{c})");
    fhConeSumPtEtaUENormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtEtaUENormCluster) ;
    
    fhConeSumPtPhiUENormCluster  = new TH2F
    ("hConeSumPtPhiUENormCluster",
     Form("Clusters #Sigma #it{p}_{T} in normalized #varphi band, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtPhiUENormCluster->SetYTitle("#Sigma #it{p}_{T}^{#varphi-band}_{norm} (GeV/#it{c})");
    fhConeSumPtPhiUENormCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtPhiUENormCluster) ;
    
    fhConeSumPtEtaUESubCluster  = new TH2F
    ("hConeSumPtEtaUESubCluster",
     Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtEtaUESubCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtEtaUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtEtaUESubCluster) ;
    
    fhConeSumPtPhiUESubCluster  = new TH2F
    ("hConeSumPtPhiUESubCluster",
     Form("Clusters #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtPhiUESubCluster->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtPhiUESubCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtPhiUESubCluster) ;
    
    if ( fFillEtaPhiHistograms )
    {
      fhConeSumPtEtaUESubClusterTrigEtaPhi  = new TH2F
      ("hConeSumPtEtaUESubClusterTrigEtaPhi",
       Form("Trigger #eta vs #varphi, Clusters #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhConeSumPtEtaUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
      fhConeSumPtEtaUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtEtaUESubClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtEtaUESubClusterTrigEtaPhi) ;
      
      fhConeSumPtPhiUESubClusterTrigEtaPhi  = new TH2F
      ("hConeSumPtPhiUESubClusterTrigEtaPhi",
       Form("Trigger #eta vs #varphi, Clusters #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhConeSumPtPhiUESubClusterTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
      fhConeSumPtPhiUESubClusterTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtPhiUESubClusterTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtPhiUESubClusterTrigEtaPhi) ;
      
      fhFractionClusterOutConeEta  = new TH2F
      ("hFractionClusterOutConeEta",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionClusterOutConeEta->SetYTitle("#it{fraction}");
      fhFractionClusterOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConeEta) ;
      
      fhFractionClusterOutConeEtaTrigEtaPhi  = new TH2F
      ("hFractionClusterOutConeEtaTrigEtaPhi",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionClusterOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConeEtaTrigEtaPhi) ;
      
      fhFractionClusterOutConePhi  = new TH2F
      ("hFractionClusterOutConePhi",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #varphi acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionClusterOutConePhi->SetYTitle("#it{fraction}");
      fhFractionClusterOutConePhi->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionClusterOutConePhi) ;
      
      fhFractionClusterOutConePhiTrigEtaPhi  = new TH2F
      ("hFractionClusterOutConePhiTrigEtaPhi",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of clusters #varphi acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionClusterOutConePhiTrigEtaPhi->SetZTitle("#it{fraction}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionClusterOutConePhiTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionClusterOutConePhiTrigEtaPhi) ;
      
      fhConeSumPtSubvsConeSumPtTotPhiCluster = new TH2F
      ("hConeSumPtSubvsConeSumPtTotPhiCluster",
       Form("#Sigma #it{p}_{T} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiCluster);
      
      fhConeSumPtSubNormvsConeSumPtTotPhiCluster = new TH2F
      ("hConeSumPtSubNormvsConeSumPtTotPhiCluster",
       Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubNormvsConeSumPtTotPhiCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiCluster);
      
      fhConeSumPtSubvsConeSumPtTotEtaCluster = new TH2F
      ("hConeSumPtSubvsConeSumPtTotEtaCluster",
       Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaCluster);
      
      fhConeSumPtSubNormvsConeSumPtTotEtaCluster = new TH2F
      ("hConeSumPtSubNormvsConeSumPtTotEtaCluster",
       Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubNormvsConeSumPtTotEtaCluster->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaCluster);
    }
    
  } // clusters
  
  if ( fPartInCone == kNeutralAndCharged )
  {
    fhPtTrackInCone  = new TH2F
    ("hPtTrackInCone",
     Form("#it{p}_{T} of tracks in isolation cone for #it{R} = %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtTrackInCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtTrackInCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtTrackInCone) ;
    
    fhConePtLeadTrack  = new TH2F
    ("hConeLeadPtTrack",
     Form("Track leading in isolation cone for #it{R} = %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhConePtLeadTrack->SetYTitle("#it{p}_{T, leading} (GeV/#it{c})");
    fhConePtLeadTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConePtLeadTrack) ;
    
    fhConeSumPtTrack  = new TH2F
    ("hConePtSumTrack",
     Form("Track #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPtTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtTrack->SetXTitle("#it{p}_{T, trigger} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtTrack) ;
  }
  
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
    
    fhPtInPerpCone  = new TH2F
    ("hPtInPerpCone",
     Form("#it{p}_{T} in isolation cone at #pm 45 degree #varphi from trigger particle, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptinconebins,ptinconemin,ptinconemax);
    fhPtInPerpCone->SetYTitle("#it{p}_{T in cone} (GeV/#it{c})");
    fhPtInPerpCone->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtInPerpCone) ;
    
    fhPerpConeSumPt  = new TH2F
    ("hPerpConePtSum",
     Form("#Sigma #it{p}_{T} in 2 isolation cones at #pm 45 degree #varphi from trigger particle, norm. to 1 cone, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,nptsumbins,ptsummin,ptsummax);
    fhPerpConeSumPt->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhPerpConeSumPt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPerpConeSumPt) ;
    
    if ( fFillEtaPhiHistograms )
    {
      fhEtaPhiInPerpCone = new TH2F
      ("hEtaPhiInPerpCone",
       Form("#eta vs #varphi of all Tracks"),
       netabins,-1,1,nphibins,0,TMath::TwoPi());
      fhEtaPhiInPerpCone->SetXTitle("#eta");
      fhEtaPhiInPerpCone->SetYTitle("#varphi (rad)");
      outputContainer->Add(fhEtaPhiInPerpCone) ;
    }
    
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
    fhConeSumPtEtaUENormTrack  = new TH2F
    ("hConeSumPtEtaUENormTrack",
     Form("Tracks #Sigma #it{p}_{T} in normalized #eta band, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtEtaUENormTrack->SetYTitle("#Sigma #it{p}_{T}^{#eta-band}_{norm} (GeV/#it{c})");
    fhConeSumPtEtaUENormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtEtaUENormTrack) ;
    
    fhConeSumPtPhiUENormTrack  = new TH2F
    ("hConeSumPtPhiUENormTrack",
     Form("Tracks #Sigma #it{p}_{T} in normalized #varphi band, #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtPhiUENormTrack->SetYTitle("#Sigma #it{p}_{T}^{#varphi-band}_{norm} (GeV/#it{c})");
    fhConeSumPtPhiUENormTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtPhiUENormTrack) ;
    
    fhConeSumPtEtaUESubTrack  = new TH2F
    ("hConeSumPtEtaUESubTrack",
     Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtEtaUESubTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtEtaUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtEtaUESubTrack) ;
    
    fhConeSumPtPhiUESubTrack  = new TH2F
    ("hConeSumPtPhiUESubTrack",
     Form("Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
    fhConeSumPtPhiUESubTrack->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
    fhConeSumPtPhiUESubTrack->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtPhiUESubTrack) ;
    
    if ( fFillEtaPhiHistograms )
    {
      fhConeSumPtEtaUESubTrackTrigEtaPhi  = new TH2F
      ("hConeSumPtEtaUESubTrackTrigEtaPhi",
       Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhConeSumPtEtaUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
      fhConeSumPtEtaUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtEtaUESubTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtEtaUESubTrackTrigEtaPhi) ;
      
      fhConeSumPtPhiUESubTrackTrigEtaPhi  = new TH2F
      ("hConeSumPtPhiUESubTrackTrigEtaPhi",
       Form("Trigger #eta vs #varphi, Tracks #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhConeSumPtPhiUESubTrackTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}");
      fhConeSumPtPhiUESubTrackTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhConeSumPtPhiUESubTrackTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhConeSumPtPhiUESubTrackTrigEtaPhi) ;
      
      fhFractionTrackOutConeEta  = new TH2F
      ("hFractionTrackOutConeEta",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance",fConeSize),
       nptbins,ptmin,ptmax,100,0,1);
      fhFractionTrackOutConeEta->SetYTitle("#it{fraction}");
      fhFractionTrackOutConeEta->SetXTitle("#it{p}_{T,trigger} (GeV/#it{c})");
      outputContainer->Add(fhFractionTrackOutConeEta) ;
      
      fhFractionTrackOutConeEtaTrigEtaPhi  = new TH2F
      ("hFractionTrackOutConeEtaTrigEtaPhi",
       Form("Fraction of the isolation cone #it{R} =  %2.2f, out of tracks #eta acceptance, in trigger #eta-#varphi ",fConeSize),
       netabins,etamin,etamax,nphibins,phimin,phimax);
      fhFractionTrackOutConeEtaTrigEtaPhi->SetZTitle("#it{fraction}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetXTitle("#eta_{trigger}");
      fhFractionTrackOutConeEtaTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
      outputContainer->Add(fhFractionTrackOutConeEtaTrigEtaPhi) ;
      
      fhConeSumPtSubvsConeSumPtTotPhiTrack = new TH2F
      ("hConeSumPtSubvsConeSumPtTotPhiTrack",
       Form("#Sigma #it{p}_{T} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubvsConeSumPtTotPhiTrack);
      
      fhConeSumPtSubNormvsConeSumPtTotPhiTrack = new TH2F
      ("hConeSumPtSubNormvsConeSumPtTotPhiTrack",
       Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #varphi band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubNormvsConeSumPtTotPhiTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotPhiTrack);
      
      fhConeSumPtSubvsConeSumPtTotEtaTrack = new TH2F
      ("hConeSumPtSubvsConeSumPtTotEtaTrack",
       Form("#Sigma #it{p}_{T} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubvsConeSumPtTotEtaTrack);
      
      fhConeSumPtSubNormvsConeSumPtTotEtaTrack = new TH2F
      ("hConeSumPtSubNormvsConeSumPtTotEtaTrack",
       Form("#Sigma #it{p}_{T, norm} in cone after bkg sub from #eta band vs #Sigma #it{p}_{T} in cone before bkg sub, R=%2.2f",fConeSize),
       nptsumbins,ptsummin,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetXTitle("#Sigma #it{p}_{T, tot} (GeV/#it{c})");
      fhConeSumPtSubNormvsConeSumPtTotEtaTrack->SetYTitle("#Sigma #it{p}_{T, sub norm} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtSubNormvsConeSumPtTotEtaTrack);
      
    }
  } // charged
  
  if ( fPartInCone == kNeutralAndCharged )
  {
    fhConeSumPtClustervsTrack   = new TH2F
    ("hConePtSumClustervsTrack",
     Form("Track vs Cluster #Sigma #it{p}_{T} in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
    fhConeSumPtClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
    fhConeSumPtClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtClustervsTrack) ;
    
    fhConeSumPtClusterTrackFrac   = new TH2F
    ("hConePtSumClusterTrackFraction",
     Form("#Sigma #it{p}_{T}^{cluster}/#Sigma #it{p}_{T}^{track} in isolation cone for #it{R} =  %2.2f",fConeSize),
     nptbins,ptmin,ptmax,200,0,5);
    fhConeSumPtClusterTrackFrac->SetYTitle("#Sigma #it{p}^{cluster}_{T} /#Sigma #it{p}_{T}^{track}");
    fhConeSumPtClusterTrackFrac->SetXTitle("#it{p}^{trigger}_{T} (GeV/#it{c})");
    outputContainer->Add(fhConeSumPtClusterTrackFrac) ;
    
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
    
    // UE subtraction
    if ( fICMethod >= kSumBkgSubIC )
    {
      fhConeSumPtEtaUESub  = new TH2F
      ("hConeSumPtEtaUESub",
       Form("#Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtEtaUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtEtaUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtEtaUESub) ;
      
      fhConeSumPtPhiUESub  = new TH2F
      ("hConeSumPtPhiUESub",
       Form("#Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
       nptbins,ptmin,ptmax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
      fhConeSumPtPhiUESub->SetYTitle("#Sigma #it{p}_{T} (GeV/#it{c})");
      fhConeSumPtPhiUESub->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhConeSumPtPhiUESub) ;
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtEtaUESubTrigEtaPhi  = new TH2F
        ("hConeSumPtEtaUESubTrigEtaPhi",
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} after bkg subtraction from #eta band in the isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtEtaUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}(GeV/#it{c})");
        fhConeSumPtEtaUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtEtaUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtEtaUESubTrigEtaPhi) ;
        
        fhConeSumPtPhiUESubTrigEtaPhi  = new TH2F
        ("hConeSumPtPhiUESubTrigEtaPhi",
         Form("Trigger #eta vs #varphi, #Sigma #it{p}_{T} after bkg subtraction from #varphi band in the isolation cone for #it{R} =  %2.2f",fConeSize),
         netabins,etamin,etamax,nphibins,phimin,phimax);
        fhConeSumPtPhiUESubTrigEtaPhi->SetZTitle("#Sigma #it{p}_{T}(GeV/#it{c})");
        fhConeSumPtPhiUESubTrigEtaPhi->SetXTitle("#eta_{trigger}");
        fhConeSumPtPhiUESubTrigEtaPhi->SetYTitle("#varphi_{trigger} (rad)");
        outputContainer->Add(fhConeSumPtPhiUESubTrigEtaPhi) ;
        
        fhConeSumPtEtaUESubClustervsTrack   = new TH2F
        ("hConePtSumEtaUESubClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #eta band in isolation cone for #it{R} =  %2.2f",fConeSize),
         1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtEtaUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtEtaUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtEtaUESubClustervsTrack) ;
        
        fhConeSumPtPhiUESubClustervsTrack   = new TH2F
        ("hConePhiUESubPtSumClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} UE sub #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         1.2*nptsumbins,-ptsummax*0.2,ptsummax,1.2*nptsumbins,-ptsummax*0.2,ptsummax);
        fhConeSumPtPhiUESubClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhConeSumPtPhiUESubClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhConeSumPtPhiUESubClustervsTrack) ;
        
        fhEtaBandClustervsTrack   = new TH2F
        ("hEtaBandClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhEtaBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandClustervsTrack) ;
        
        fhPhiBandClustervsTrack   = new TH2F
        ("hPhiBandClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax*4,nptsumbins,ptsummin,ptsummax*8);
        fhPhiBandClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPhiBandClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandClustervsTrack) ;
        
        fhEtaBandNormClustervsTrack   = new TH2F
        ("hEtaBandNormClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} in  #eta band in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhEtaBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhEtaBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhEtaBandNormClustervsTrack) ;
        
        fhPhiBandNormClustervsTrack   = new TH2F
        ("hPhiBandNormClustervsTrack",
         Form("Track vs Cluster #Sigma #it{p}_{T} in  #varphi band in isolation cone for #it{R} =  %2.2f",fConeSize),
         nptsumbins,ptsummin,ptsummax,nptsumbins,ptsummin,ptsummax);
        fhPhiBandNormClustervsTrack->SetXTitle("#Sigma #it{p}_{T}^{cluster} (GeV/#it{c})");
        fhPhiBandNormClustervsTrack->SetYTitle("#Sigma #it{p}_{T}^{track} (GeV/#it{c})");
        outputContainer->Add(fhPhiBandNormClustervsTrack) ;
      }
    } // UE
  } // neutral + charged
  
    return outputContainer;
  }
  
  //____________________________________________
  // Put data member values in string to keep
  // in output container.
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
  snprintf(onePar,buffersize,"fSumPtThreshold>%2.2f;<%2.2f;",fSumPtThreshold,fSumPtThresholdMax) ;
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
  fSumPtThresholdMax    = 10000. ;
  fPtFraction           = 0.1 ;
  fPartInCone           = kNeutralAndCharged;
  fICMethod             = kSumPtIC; // 0 pt threshol method, 1 cone pt sum method
  fFracIsThresh         = 1;
  fDistMinToTrigger     = -1.; // no effect
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
  Float_t ptC   = pCandidate->Pt() ;
  Float_t phiC  = pCandidate->Phi() ;
  if ( phiC < 0 ) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  Float_t coneptsumCluster = 0;
  Float_t coneptsumTrack   = 0;
  Float_t coneptLeadCluster= 0;
  Float_t coneptLeadTrack  = 0;
  
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
  
  CalculateTrackSignalInCone   (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgTrk,
                                nPart              , nfrac,
                                coneptsumTrack     , coneptLeadTrack,
                                etaBandPtSumTrack  , phiBandPtSumTrack,
                                perpPtSumTrack    , histoWeight);
  
  CalculateCaloSignalInCone    (pCandidate         , reader,
                                bFillAOD           , useRefs, 
                                aodArrayRefName    , bgCls,
                                calorimeter        , pid, 
                                nPart              , nfrac,
                                coneptsumCluster   , coneptLeadCluster,
                                etaBandPtSumCluster, phiBandPtSumCluster,
                                histoWeight);
  
  coneptsum = coneptsumCluster + coneptsumTrack;
  
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
  
  //-------------------------------------------------------------------
  // Check isolation, depending on selected isolation criteria requested
  //-------------------------------------------------------------------

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
  else if ( fICMethod == kSumBkgSubIC )
  {
    Double_t coneptsumBkg = 0.;
    
    Float_t  etaBandPtSumTrackNorm   = 0;
    Float_t  phiBandPtSumTrackNorm   = 0;
    Float_t  etaBandPtSumClusterNorm = 0;
    Float_t  phiBandPtSumClusterNorm = 0;
    
    Float_t  excessFracEtaTrack   = 1;
    Float_t  excessFracPhiTrack   = 1;
    Float_t  excessFracEtaCluster = 1;
    Float_t  excessFracPhiCluster = 1;
    
    Float_t  coneptsumTrackSubPhi   = 0 ;
    Float_t  coneptsumTrackSubEta   = 0 ;
    Float_t  coneptsumClusterSubPhi = 0 ;
    Float_t  coneptsumClusterSubEta = 0 ;
    
    Float_t  coneptsumTrackSubPhiNorm   = 0 ;
    Float_t  coneptsumTrackSubEtaNorm   = 0 ;
    Float_t  coneptsumClusterSubPhiNorm = 0 ;
    Float_t  coneptsumClusterSubEtaNorm = 0 ; 
    
    // Normalize background to cone area
    if     ( fPartInCone != kOnlyCharged )
    {
      CalculateUEBandClusterNormalization(reader, etaC, phiC,
                                          phiBandPtSumCluster    , etaBandPtSumCluster,
                                          phiBandPtSumClusterNorm, etaBandPtSumClusterNorm,
                                          excessFracPhiCluster   , excessFracEtaCluster    );
      
      pCandidate->SetNeutralConeExcessAreaEta(excessFracEtaCluster);
      pCandidate->SetNeutralConeExcessAreaPhi(excessFracPhiCluster);
      
      coneptsumClusterSubPhi = coneptsumCluster - phiBandPtSumClusterNorm;
      coneptsumClusterSubEta = coneptsumCluster - etaBandPtSumClusterNorm;
      
      //      printf("Cluster: sumpT %2.2f, phi: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
      //             "\t eta: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
      //             "\t subtracted: eta %2.2f, phi %2.2f\n",
      //             coneptsumCluster, phiBandPtSumCluster, phiBandPtSumClusterNorm, excessFracPhiCluster, 
      //             etaBandPtSumCluster, etaBandPtSumClusterNorm, excessFracEtaCluster,
      //             coneptsumClusterSubEta,coneptsumClusterSubPhi);
      
      if ( fFillHistograms )
      {
        fhConeSumPtPhiUENormCluster  ->Fill(ptC, phiBandPtSumClusterNorm, histoWeight);
        fhConeSumPtPhiUESubCluster   ->Fill(ptC, coneptsumClusterSubPhi , histoWeight);
        fhConeSumPtEtaUENormCluster  ->Fill(ptC, etaBandPtSumClusterNorm, histoWeight);
        fhConeSumPtEtaUESubCluster   ->Fill(ptC, coneptsumClusterSubEta , histoWeight);
        
        if ( coneptsumCluster != 0 )
        {
          coneptsumClusterSubPhiNorm = coneptsumClusterSubPhi/coneptsumCluster;
          coneptsumClusterSubEtaNorm = coneptsumClusterSubEta/coneptsumCluster;
        }
        
        if (fFillEtaPhiHistograms )
        {
          fhConeSumPtPhiUESubClusterTrigEtaPhi ->Fill(etaC, phiC, coneptsumClusterSubPhi *histoWeight); // check
          fhConeSumPtEtaUESubClusterTrigEtaPhi ->Fill(etaC, phiC, coneptsumClusterSubEta *histoWeight); // check
          
          fhFractionClusterOutConeEta          ->Fill(ptC ,        excessFracEtaCluster-1, histoWeight);
          fhFractionClusterOutConeEtaTrigEtaPhi->Fill(etaC, phiC, (excessFracEtaCluster-1)*histoWeight); // check
          fhFractionClusterOutConePhi          ->Fill(ptC ,        excessFracPhiCluster-1, histoWeight);
          fhFractionClusterOutConePhiTrigEtaPhi->Fill(etaC, phiC, (excessFracPhiCluster-1)*histoWeight); // check
          
          fhConeSumPtSubvsConeSumPtTotPhiCluster    ->Fill(coneptsumCluster,coneptsumClusterSubPhi    , histoWeight);
          fhConeSumPtSubNormvsConeSumPtTotPhiCluster->Fill(coneptsumCluster,coneptsumClusterSubPhiNorm, histoWeight);
          fhConeSumPtSubvsConeSumPtTotEtaCluster    ->Fill(coneptsumCluster,coneptsumClusterSubEta    , histoWeight);
          fhConeSumPtSubNormvsConeSumPtTotEtaCluster->Fill(coneptsumCluster,coneptsumClusterSubEtaNorm, histoWeight);
        }
      } // histograms
    } // clusters
    
    if     ( fPartInCone != kOnlyNeutral )
    {
      CalculateUEBandTrackNormalization(reader, etaC, phiC,
                                        phiBandPtSumTrack    , etaBandPtSumTrack  ,
                                        phiBandPtSumTrackNorm, etaBandPtSumTrackNorm,
                                        excessFracPhiTrack   , excessFracEtaTrack    );
      
      pCandidate->SetChargedConeExcessAreaEta(excessFracEtaTrack);
      pCandidate->SetChargedConeExcessAreaPhi(excessFracPhiTrack);
      
      coneptsumTrackSubPhi = coneptsumTrack - phiBandPtSumTrackNorm;
      coneptsumTrackSubEta = coneptsumTrack - etaBandPtSumTrackNorm;
      
//      printf("Track: sumpT %2.2f, phi: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
//             "\t eta: sum pT %2.2f, sumpT norm %2.2f, excess %2.2f;\n"
//             "\t subtracted: eta %2.2f, phi %2.2f\n",
//             coneptsumTrack, phiBandPtSumTrack, phiBandPtSumTrackNorm, excessFracPhiTrack, 
//             etaBandPtSumTrack, etaBandPtSumTrackNorm, excessFracEtaTrack,
//             coneptsumTrackSubEta,coneptsumTrackSubPhi);  
      
      if ( fFillHistograms )
      {
        fhConeSumPtPhiUESubTrack  ->Fill(ptC, coneptsumTrackSubPhi, histoWeight);
        fhConeSumPtEtaUESubTrack  ->Fill(ptC, coneptsumTrackSubEta, histoWeight);
        
        fhConeSumPtPhiUENormTrack ->Fill(ptC, phiBandPtSumTrackNorm , histoWeight);
        fhConeSumPtEtaUENormTrack ->Fill(ptC, etaBandPtSumTrackNorm , histoWeight);
        
        if ( coneptsumTrack > 0 )
        {
          coneptsumTrackSubPhiNorm = coneptsumTrackSubPhi/coneptsumTrack;
          coneptsumTrackSubEtaNorm = coneptsumTrackSubEta/coneptsumTrack;
        }
        
        if ( fFillEtaPhiHistograms )
        {       
          fhConeSumPtPhiUESubTrackTrigEtaPhi ->Fill(etaC, phiC, coneptsumTrackSubPhi *histoWeight); // check
          fhConeSumPtEtaUESubTrackTrigEtaPhi ->Fill(etaC, phiC, coneptsumTrackSubEta *histoWeight); // check
          
          fhFractionTrackOutConeEta          ->Fill(ptC ,      excessFracEtaTrack-1, histoWeight);
          fhFractionTrackOutConeEtaTrigEtaPhi->Fill(etaC, phiC,excessFracEtaTrack-1 *histoWeight); // check
          
          fhConeSumPtSubvsConeSumPtTotPhiTrack    ->Fill(coneptsumTrack, coneptsumTrackSubPhi    , histoWeight);
          fhConeSumPtSubNormvsConeSumPtTotPhiTrack->Fill(coneptsumTrack, coneptsumTrackSubPhiNorm, histoWeight);
          fhConeSumPtSubvsConeSumPtTotEtaTrack    ->Fill(coneptsumTrack, coneptsumTrackSubEta    , histoWeight);
          fhConeSumPtSubNormvsConeSumPtTotEtaTrack->Fill(coneptsumTrack, coneptsumTrackSubEtaNorm, histoWeight);
        }
      }
    }
    
    if ( fPartInCone == AliIsolationCut::kNeutralAndCharged )
    {
      Float_t sumPhiUESub = coneptsumClusterSubPhi + coneptsumTrackSubPhi;
      Float_t sumEtaUESub = coneptsumClusterSubEta + coneptsumTrackSubEta;
      
      fhConeSumPtPhiUESub ->Fill(ptC,  sumPhiUESub, histoWeight);
      fhConeSumPtEtaUESub ->Fill(ptC,  sumEtaUESub, histoWeight);
      
      if ( fFillEtaPhiHistograms )
      {
        fhConeSumPtPhiUESubTrigEtaPhi->Fill(etaC, phiC, sumPhiUESub *histoWeight); // check
        fhConeSumPtEtaUESubTrigEtaPhi->Fill(etaC, phiC, sumEtaUESub *histoWeight); // check
        
        fhEtaBandClustervsTrack    ->Fill(etaBandPtSumCluster    ,etaBandPtSumTrack    , histoWeight);
        fhPhiBandClustervsTrack    ->Fill(phiBandPtSumCluster    ,phiBandPtSumTrack    , histoWeight);
        fhEtaBandNormClustervsTrack->Fill(etaBandPtSumClusterNorm,etaBandPtSumTrackNorm, histoWeight);
        fhPhiBandNormClustervsTrack->Fill(phiBandPtSumClusterNorm,phiBandPtSumTrackNorm, histoWeight);
        
        fhConeSumPtEtaUESubClustervsTrack->Fill(coneptsumClusterSubEta, coneptsumTrackSubEta, histoWeight);
        fhConeSumPtPhiUESubClustervsTrack->Fill(coneptsumClusterSubPhi, coneptsumTrackSubPhi, histoWeight);
      }
    }
    
    if      ( fPartInCone == kOnlyCharged       ) coneptsumBkg = etaBandPtSumTrackNorm;
    else if ( fPartInCone == kOnlyNeutral       ) coneptsumBkg = etaBandPtSumClusterNorm;
    else if ( fPartInCone == kNeutralAndCharged ) coneptsumBkg = etaBandPtSumClusterNorm + etaBandPtSumTrackNorm;
    
    
    ///////////////// EXCESSS//////////////////
    //coneptsumCluster*=(coneBadCellsCoeff*excessFracEtaCluster*excessFracPhiCluster) ; // apply this correction earlier???
    // line commented out in last modif!!!
    //////////////////////////////////////////////////////////
    
    coneptsum = coneptsumCluster+coneptsumTrack;
    
    coneptsum -= coneptsumBkg;
    
    if ( coneptsum > fSumPtThreshold && coneptsum < fSumPtThresholdMax )
      isolated  =  kFALSE ;
    else
      isolated  =  kTRUE  ;
    
  }
  
  //-------------------------------------------------------------------
  // Fill histograms
  //-------------------------------------------------------------------

  if ( !fFillHistograms ) return;
  
  if ( fPartInCone == kNeutralAndCharged )
  {
    fhConeSumPtClustervsTrack ->Fill(coneptsumCluster, coneptsumTrack , histoWeight);
    fhConePtLeadClustervsTrack->Fill(coneptLeadCluster,coneptLeadTrack, histoWeight);
    
    if(coneptsumTrack  > 0) fhConeSumPtClusterTrackFrac ->Fill(ptC, coneptsumCluster /coneptsumTrack , histoWeight);
    if(coneptLeadTrack > 0) fhConePtLeadClusterTrackFrac->Fill(ptC, coneptLeadCluster/coneptLeadTrack, histoWeight);
  }
  
  fhConeSumPt              ->Fill(ptC,        coneptsum, histoWeight);
  fhConeSumPtTrigEtaPhi    ->Fill(etaC, phiC, coneptsum *histoWeight); // check

  Float_t coneptLead = coneptLeadTrack;
  if(coneptLeadCluster > coneptLeadTrack) 
    coneptLead = coneptLeadCluster;
  
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
  printf("pT threshold       =     >%2.1f;<%2.1f\n", fPtThreshold   ,   fPtThresholdMax) ;
  printf("Sum pT threshold   =     >%2.1f;<%2.1f\n", fSumPtThreshold,fSumPtThresholdMax) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction ) ;
  printf("particle type in cone =  %d\n",    fPartInCone ) ;
  printf("using fraction for high pt leading instead of frac ? %i\n",fFracIsThresh);
  printf("minimum distance to candidate, R>%1.2f\n",fDistMinToTrigger);
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



