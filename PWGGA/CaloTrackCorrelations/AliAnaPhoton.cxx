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

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
// Produces input for other analysis classes like AliAnaPi0,
// AliAnaParticleHadronCorrelation ...
//
// -- Author: Gustavo Conesa (LNF-INFN)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TH2F.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include "TParticle.h"
#include "TDatabasePDG.h"

// --- Analysis system ---
#include "AliAnaPhoton.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliCaloPID.h"
#include "AliMCAnalysisUtils.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODMCParticle.h"
#include "AliMixedEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"

// --- Detectors ---
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnaPhoton)

//____________________________
AliAnaPhoton::AliAnaPhoton() :
AliAnaCaloTrackCorrBaseClass(), fCalorimeter(""),
fMinDist(0.),                 fMinDist2(0.),                fMinDist3(0.),
fRejectTrackMatch(0),         fFillTMHisto(kFALSE),
fTimeCutMin(-10000),          fTimeCutMax(10000),
fNCellsCut(0),
fNLMCutMin(-1),               fNLMCutMax(10),
fFillSSHistograms(kFALSE),    fFillOnlySimpleSSHisto(1),
fFillPileUpHistograms(0),
fNOriginHistograms(8),        fNPrimaryHistograms(4),
// Histograms

// Control histograms
fhNCellsE(0),                 fhCellsE(0),
fhMaxCellDiffClusterE(0),     fhTimePt(0),                  fhEtaPhi(0),

fhEPhoton(0),                 fhPtPhoton(0),
fhPhiPhoton(0),               fhEtaPhoton(0),
fhEtaPhiPhoton(0),            fhEtaPhi05Photon(0),
fhPtCentralityPhoton(0),      fhPtEventPlanePhoton(0),

// Shower shape histograms
fhNLocMax(0),
fhDispE(0),                   fhLam0E(0),                   fhLam1E(0),
fhDispETRD(0),                fhLam0ETRD(0),                fhLam1ETRD(0),
fhDispETM(0),                 fhLam0ETM(0),                 fhLam1ETM(0),
fhDispETMTRD(0),              fhLam0ETMTRD(0),              fhLam1ETMTRD(0),

fhNCellsLam0LowE(0),          fhNCellsLam1LowE(0),          fhNCellsDispLowE(0),
fhNCellsLam0HighE(0),         fhNCellsLam1HighE(0),         fhNCellsDispHighE(0),

fhEtaLam0LowE(0),             fhPhiLam0LowE(0),
fhEtaLam0HighE(0),            fhPhiLam0HighE(0),
fhLam0DispLowE(0),            fhLam0DispHighE(0),
fhLam1Lam0LowE(0),            fhLam1Lam0HighE(0),
fhDispLam1LowE(0),            fhDispLam1HighE(0),
fhDispEtaE(0),                fhDispPhiE(0),
fhSumEtaE(0),                 fhSumPhiE(0),                 fhSumEtaPhiE(0),
fhDispEtaPhiDiffE(0),         fhSphericityE(0),
fhDispSumEtaDiffE(0),         fhDispSumPhiDiffE(0),

// MC histograms
fhMCPhotonELambda0NoOverlap(0),       fhMCPhotonELambda0TwoOverlap(0),      fhMCPhotonELambda0NOverlap(0),
// Embedding
fhEmbeddedSignalFractionEnergy(0),
fhEmbedPhotonELambda0FullSignal(0),   fhEmbedPhotonELambda0MostlySignal(0),
fhEmbedPhotonELambda0MostlyBkg(0),    fhEmbedPhotonELambda0FullBkg(0),
fhEmbedPi0ELambda0FullSignal(0),      fhEmbedPi0ELambda0MostlySignal(0),
fhEmbedPi0ELambda0MostlyBkg(0),       fhEmbedPi0ELambda0FullBkg(0),

fhTimePtPhotonNoCut(0),               fhTimePtPhotonSPD(0),
fhTimeNPileUpVertSPD(0),              fhTimeNPileUpVertTrack(0),
fhPtPhotonNPileUpSPDVtx(0),           fhPtPhotonNPileUpTrkVtx(0),
fhPtPhotonNPileUpSPDVtxTimeCut(0),    fhPtPhotonNPileUpTrkVtxTimeCut(0),
fhPtPhotonNPileUpSPDVtxTimeCut2(0),   fhPtPhotonNPileUpTrkVtxTimeCut2(0),

fhEClusterSM(0),                      fhEPhotonSM(0),
fhPtClusterSM(0),                     fhPtPhotonSM(0)
{
  //default ctor
  
  for(Int_t i = 0; i < 14; i++)
  {
    fhMCPt     [i] = 0;
    fhMCE      [i] = 0;
    fhMCPhi    [i] = 0;
    fhMCEta    [i] = 0;
    fhMCDeltaE [i] = 0;
    fhMCDeltaPt[i] = 0;
    fhMC2E     [i] = 0;
    fhMC2Pt    [i] = 0;
  }
  
  for(Int_t i = 0; i < 7; i++)
  {
    fhPtPrimMC [i] = 0;
    fhEPrimMC  [i] = 0;
    fhPhiPrimMC[i] = 0;
    fhEtaPrimMC[i] = 0;
    fhYPrimMC  [i] = 0;
    
    fhPtPrimMCAcc [i] = 0;
    fhEPrimMCAcc  [i] = 0;
    fhPhiPrimMCAcc[i] = 0;
    fhEtaPrimMCAcc[i] = 0;
    fhYPrimMCAcc  [i] = 0;
    
    fhDispEtaDispPhi[i] = 0;
    fhLambda0DispPhi[i] = 0;
    fhLambda0DispEta[i] = 0;

    fhPtPhotonPileUp[i] = 0;
    fhClusterTimeDiffPhotonPileUp [i] = 0;

    for(Int_t j = 0; j < 6; j++)
    {
      fhMCDispEtaDispPhi[i][j] = 0;
      fhMCLambda0DispEta[i][j] = 0;
      fhMCLambda0DispPhi[i][j] = 0;
    }
  }
  
  for(Int_t i = 0; i < 6; i++)
  {
    fhMCELambda0    [i]                  = 0;
    fhMCELambda1    [i]                  = 0;
    fhMCEDispersion [i]                  = 0;
    fhMCNCellsE     [i]                  = 0;
    fhMCMaxCellDiffClusterE[i]           = 0;
    fhLambda0DispEta[i]                  = 0;
    fhLambda0DispPhi[i]                  = 0;
    
    fhMCLambda0vsClusterMaxCellDiffE0[i] = 0;
    fhMCLambda0vsClusterMaxCellDiffE2[i] = 0;
    fhMCLambda0vsClusterMaxCellDiffE6[i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE0 [i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE2 [i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE6 [i] = 0;
    
    fhMCEDispEta       [i]               = 0;
    fhMCEDispPhi       [i]               = 0;
    fhMCESumEtaPhi     [i]               = 0;
    fhMCEDispEtaPhiDiff[i]               = 0;
    fhMCESphericity    [i]               = 0;
  }
  
  for(Int_t i = 0; i < 5; i++)
  {
    fhClusterCutsE [i] = 0;
    fhClusterCutsPt[i] = 0;
  }
  
  // Track matching residuals
  for(Int_t i = 0; i < 2; i++)
  {
    fhTrackMatchedDEta   [i] = 0;             fhTrackMatchedDPhi   [i] = 0;         fhTrackMatchedDEtaDPhi   [i] = 0;
    fhTrackMatchedDEtaNeg[i] = 0;             fhTrackMatchedDPhiNeg[i] = 0;         fhTrackMatchedDEtaDPhiNeg[i] = 0;
    fhTrackMatchedDEtaPos[i] = 0;             fhTrackMatchedDPhiPos[i] = 0;         fhTrackMatchedDEtaDPhiPos[i] = 0;
    fhTrackMatchedDEtaTRD[i] = 0;             fhTrackMatchedDPhiTRD[i] = 0;
    fhTrackMatchedDEtaMCOverlap[i] = 0;       fhTrackMatchedDPhiMCOverlap[i] = 0;
    fhTrackMatchedDEtaMCNoOverlap[i] = 0;     fhTrackMatchedDPhiMCNoOverlap[i] = 0;
    fhTrackMatchedDEtaMCConversion[i] = 0;    fhTrackMatchedDPhiMCConversion[i] = 0;
    fhTrackMatchedMCParticle[i] = 0;          fhTrackMatchedMCParticle[i] = 0;
    fhdEdx[i] = 0;                            fhEOverP[i] = 0;
    fhEOverPTRD[i] = 0;
  }
  
  //Initialize parameters
  InitParameters();
  
}

//_________________________________________________________________________________________
Bool_t  AliAnaPhoton::ClusterSelected(AliVCluster* calo, TLorentzVector mom, Int_t nMaxima)
{
  //Select clusters if they pass different cuts
  
  Float_t ptcluster  = mom.Pt();
  Float_t ecluster   = mom.E();
  Float_t etacluster = mom.Eta();
  Float_t phicluster = mom.Phi();

  if(phicluster < 0) phicluster+=TMath::TwoPi();
  
  Bool_t matched = IsTrackMatched(calo,GetReader()->GetInputEvent());
  
  if(GetDebug() > 2)
    printf("AliAnaPhoton::ClusterSelected() - Current Event %d; Before selection : E %2.2f, pT %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           ecluster,ptcluster, phicluster*TMath::RadToDeg(),etacluster);
  
  fhClusterCutsE [1]->Fill( ecluster);
  fhClusterCutsPt[1]->Fill(ptcluster);
  
  if(ecluster > 0.5) fhEtaPhi->Fill(etacluster, phicluster);
  
  Int_t   nSM  = GetModuleNumber(calo);
  if(nSM < GetCaloUtils()->GetNumberOfSuperModulesUsed() && nSM >=0)
  {
    fhEClusterSM ->Fill(ecluster ,nSM);
    fhPtClusterSM->Fill(ptcluster,nSM);
  }
  
  //.......................................
  //If too small or big energy, skip it
  if(ecluster < GetMinEnergy() || ecluster > GetMaxEnergy() ) return kFALSE ;
  
  if(GetDebug() > 2) printf("\t Cluster %d Pass E Cut \n",calo->GetID());
  
  fhClusterCutsE [2]->Fill( ecluster);
  fhClusterCutsPt[2]->Fill(ptcluster);
  
  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if(tof < fTimeCutMin || tof > fTimeCutMax) return kFALSE;
  
  if(GetDebug() > 2)  printf("\t Cluster %d Pass Time Cut \n",calo->GetID());
  
  fhClusterCutsE [3]->Fill( ecluster);
  fhClusterCutsPt[3]->Fill(ptcluster);
  
  //.......................................
  if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) return kFALSE;
  
  if(GetDebug() > 2) printf("\t Cluster %d Pass NCell Cut \n",calo->GetID());
  
  fhClusterCutsE [4]->Fill( ecluster);
  fhClusterCutsPt[4]->Fill(ptcluster);
  
  if(nMaxima < fNLMCutMin || nMaxima > fNLMCutMax) return kFALSE ;
  if(GetDebug() > 2) printf(" \t Cluster %d pass NLM %d of out of range \n",calo->GetID(), nMaxima);
  
  fhClusterCutsE [5]->Fill( ecluster);
  fhClusterCutsPt[5]->Fill(ptcluster);
  
  //.......................................
  //Check acceptance selection
  if(IsFiducialCutOn())
  {
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(! in ) return kFALSE ;
  }
  
  if(GetDebug() > 2) printf("\t Fiducial cut passed \n");
  
  fhClusterCutsE [6]->Fill( ecluster);
  fhClusterCutsPt[6]->Fill(ptcluster);
  
  //.......................................
  //Skip matched clusters with tracks
  
  // Fill matching residual histograms before PID cuts
  if(fFillTMHisto) FillTrackMatchingResidualHistograms(calo,0);
  
  if(fRejectTrackMatch)
  {
    if(matched)
    {
      if(GetDebug() > 2) printf("\t Reject track-matched clusters\n");
      return kFALSE ;
    }
    else
      if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  }// reject matched clusters
  
  fhClusterCutsE [7]->Fill( ecluster);
  fhClusterCutsPt[7]->Fill(ptcluster);
  
  //.......................................
  //Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  if(distBad < fMinDist)
  {//In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  }
  else if(GetDebug() > 2) printf("\t Bad channel cut passed %4.2f > %2.2f \n",distBad, fMinDist);
  
  fhClusterCutsE [8]->Fill( ecluster);
  fhClusterCutsPt[8]->Fill(ptcluster);
  
  if(GetDebug() > 0)
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           ecluster, ptcluster,mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //All checks passed, cluster selected
  return kTRUE;
  
}

//___________________________________________
void AliAnaPhoton::FillAcceptanceHistograms()
{
  //Fill acceptance histograms if MC data is available
  
  Double_t photonY   = -100 ;
  Double_t photonE   = -1 ;
  Double_t photonPt  = -1 ;
  Double_t photonPhi =  100 ;
  Double_t photonEta = -1 ;
  
  Int_t    pdg       =  0 ;
  Int_t    tag       =  0 ;
  Int_t    status    =  0 ;
  Int_t    mcIndex   =  0 ;
  Int_t    nprim     =  0 ;
  Bool_t   inacceptance = kFALSE ;
  
  TParticle        * primStack = 0;
  AliAODMCParticle * primAOD   = 0;
  TLorentzVector lv;
  
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
        printf("AliAnaPhoton::FillAcceptanceHistograms() - ESD primaries pointer not available!!\n");
        continue;
      }
      
      pdg    = primStack->GetPdgCode();
      status = primStack->GetStatusCode();
      
      if(primStack->Energy() == TMath::Abs(primStack->Pz()))  continue ; //Protection against floating point exception
      
      //printf("i %d, %s %d  %s %d \n",i, stack->Particle(i)->GetName(), stack->Particle(i)->GetPdgCode(),
      //       prim->GetName(), prim->GetPdgCode());
      
      //Photon kinematics
      primStack->Momentum(lv);
      
      photonY = 0.5*TMath::Log((primStack->Energy()+primStack->Pz())/(primStack->Energy()-primStack->Pz())) ;
    }
    else
    {
      primAOD = (AliAODMCParticle *) mcparticles->At(i);
      if(!primAOD)
      {
        printf("AliAnaPhoton::FillAcceptanceHistograms() - AOD primaries pointer not available!!\n");
        continue;
      }
      
      pdg    = primAOD->GetPdgCode();
      status = primAOD->GetStatus();
      
      if(primAOD->E() == TMath::Abs(primAOD->Pz()))  continue ; //Protection against floating point exception
      
      //Photon kinematics
      lv.SetPxPyPzE(primAOD->Px(),primAOD->Py(),primAOD->Pz(),primAOD->E());

      photonY = 0.5*TMath::Log((primAOD->E()+primAOD->Pz())/(primAOD->E()-primAOD->Pz())) ;
    }

    // Select only photons in the final state
    if(pdg != 22 ) continue ;
    
    // If too small or too large pt, skip, same cut as for data analysis
    photonPt  = lv.Pt () ;
    
    if(photonPt < GetMinPt() || photonPt > GetMaxPt() ) continue ;
    
    photonE   = lv.E  () ;
    photonEta = lv.Eta() ;
    photonPhi = lv.Phi() ;
    
    if(photonPhi < 0) photonPhi+=TMath::TwoPi();
    
    // Check if photons hit desired acceptance
    inacceptance = kTRUE;
    
    // Check same fidutial borders as in data analysis on top of real acceptance if real was requested.
    if( IsFiducialCutOn() && !GetFiducialCut()->IsInFiducialCut(lv,fCalorimeter)) inacceptance = kFALSE ;
    
    // Check if photons hit the Calorimeter acceptance
    if(IsRealCaloAcceptanceOn()) // defined on base class
    {
      if(GetReader()->ReadStack()          &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fCalorimeter, primStack)) inacceptance = kFALSE ;
      if(GetReader()->ReadAODMCParticles() &&
         !GetCaloUtils()->IsMCParticleInCalorimeterAcceptance(fCalorimeter, primAOD  )) inacceptance = kFALSE ;
    }
    
    // Get tag of this particle photon from fragmentation, decay, prompt ...
    // Set the origin of the photon.
    tag = GetMCAnalysisUtils()->CheckOrigin(i,GetReader(),fCalorimeter);
    
    if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton))
    {
      // A conversion photon from a hadron, skip this kind of photon
      // printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!\n ");
      // GetMCAnalysisUtils()->PrintMCTag(tag);
      
      continue;
    }
    
    // Consider only final state particles, but this depends on generator,
    // status 1 is the usual one, in case of not being ok, leave the possibility
    // to not consider this.
    if(status > 1) continue ; // Avoid "partonic" photons
    
    Bool_t takeIt  = kFALSE ;
    if(status == 1 && GetMCAnalysisUtils()->GetMCGenerator()!="" ) takeIt = kTRUE ;

    if      (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) continue;
    
    //Origin of photon
    if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))
    {
      mcIndex = kmcPPrompt;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
    {
      mcIndex = kmcPFragmentation ;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
    {
      mcIndex = kmcPISR;
    }
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay))
    {
      mcIndex = kmcPPi0Decay;
    }
    else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) ||
              GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay)))
    {
      mcIndex = kmcPOtherDecay;
    }
    else
    {
      // Other decay but from non final state particle
      mcIndex = kmcPOtherDecay;
    }//Other origin
    
    if(!takeIt &&  (mcIndex == kmcPPi0Decay || mcIndex == kmcPOtherDecay)) takeIt = kTRUE ;

    if(!takeIt) continue ;
      
    //Fill histograms for all photons
    fhYPrimMC[kmcPPhoton]->Fill(photonPt, photonY) ;
    if(TMath::Abs(photonY) < 1.0)
    {
      fhEPrimMC  [kmcPPhoton]->Fill(photonE ) ;
      fhPtPrimMC [kmcPPhoton]->Fill(photonPt) ;
      fhPhiPrimMC[kmcPPhoton]->Fill(photonE , photonPhi) ;
      fhEtaPrimMC[kmcPPhoton]->Fill(photonE , photonEta) ;
    }
    
    if(inacceptance)
    {
      fhEPrimMCAcc  [kmcPPhoton]->Fill(photonE ) ;
      fhPtPrimMCAcc [kmcPPhoton]->Fill(photonPt) ;
      fhPhiPrimMCAcc[kmcPPhoton]->Fill(photonE , photonPhi) ;
      fhEtaPrimMCAcc[kmcPPhoton]->Fill(photonE , photonEta) ;
      fhYPrimMCAcc  [kmcPPhoton]->Fill(photonE , photonY) ;
    }//Accepted
    
    //Fill histograms for photons origin
    if(mcIndex < fNPrimaryHistograms)
    {
      fhYPrimMC[mcIndex]->Fill(photonPt, photonY) ;
      if(TMath::Abs(photonY) < 1.0)
      {
        fhEPrimMC  [mcIndex]->Fill(photonE ) ;
        fhPtPrimMC [mcIndex]->Fill(photonPt) ;
        fhPhiPrimMC[mcIndex]->Fill(photonE , photonPhi) ;
        fhEtaPrimMC[mcIndex]->Fill(photonE , photonEta) ;
      }
      
      if(inacceptance)
      {
        fhEPrimMCAcc  [mcIndex]->Fill(photonE ) ;
        fhPtPrimMCAcc [mcIndex]->Fill(photonPt) ;
        fhPhiPrimMCAcc[mcIndex]->Fill(photonE , photonPhi) ;
        fhEtaPrimMCAcc[mcIndex]->Fill(photonE , photonEta) ;
        fhYPrimMCAcc  [mcIndex]->Fill(photonE , photonY) ;
      }//Accepted
    }
  }//loop on primaries

}

//________________________________________________________________________________
void AliAnaPhoton::FillPileUpHistograms(AliVCluster* cluster, AliVCaloCells *cells)
{
  // Fill some histograms to understand pile-up
  
  TLorentzVector mom;
  cluster->GetMomentum(mom,GetVertex(0));
  Float_t pt   = mom.Pt();
  Float_t time = cluster->GetTOF()*1.e9;
  
  AliVEvent * event = GetReader()->GetInputEvent();
  
  if(GetReader()->IsPileUpFromSPD())               fhPtPhotonPileUp[0]->Fill(pt);
  if(GetReader()->IsPileUpFromEMCal())             fhPtPhotonPileUp[1]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtPhotonPileUp[2]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtPhotonPileUp[3]->Fill(pt);
  if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtPhotonPileUp[4]->Fill(pt);
  if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtPhotonPileUp[5]->Fill(pt);
  if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtPhotonPileUp[6]->Fill(pt);
  
  fhTimePtPhotonNoCut->Fill(pt,time);
  if(GetReader()->IsPileUpFromSPD()) fhTimePtPhotonSPD->Fill(pt,time);
  
  // cells inside the cluster
  Float_t maxCellFraction = 0.;
  Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell( cells, cluster, maxCellFraction);
  
  //Loop on cells inside cluster, max cell must be over 100 MeV and time in BC=0
  if(cells->GetCellAmplitude(absIdMax) > 0.1 && TMath::Abs(time) < 30)
  {
    for (Int_t ipos = 0; ipos < cluster->GetNCells(); ipos++)
    {
      Int_t absId  = cluster->GetCellsAbsId()[ipos];
      
      if( absId == absIdMax ) continue ;
      
      Double_t tcell = cells->GetCellTime(absId);
      Float_t  amp   = cells->GetCellAmplitude(absId);
      Int_t    bc    = GetReader()->GetInputEvent()->GetBunchCrossNumber();
      
      GetCaloUtils()->GetEMCALRecoUtils()->AcceptCalibrateCell(absId,bc,amp,tcell,cells);
      tcell*=1e9;
      
      Float_t diff = (time-tcell);
      
      if( cells->GetCellAmplitude(absIdMax) < 0.1 ) continue ;
      
      if(GetReader()->IsPileUpFromSPD())               fhClusterTimeDiffPhotonPileUp[0]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromEMCal())             fhClusterTimeDiffPhotonPileUp[1]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhClusterTimeDiffPhotonPileUp[2]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhClusterTimeDiffPhotonPileUp[3]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhClusterTimeDiffPhotonPileUp[4]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhClusterTimeDiffPhotonPileUp[5]->Fill(pt, diff);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhClusterTimeDiffPhotonPileUp[6]->Fill(pt, diff);

    }//loop
  }
  
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
  
  // N pile up vertices
  Int_t nVtxSPD = -1;
  Int_t nVtxTrk = -1;
  
  if      (esdEv)
  {
    nVtxSPD = esdEv->GetNumberOfPileupVerticesSPD();
    nVtxTrk = esdEv->GetNumberOfPileupVerticesTracks();
    
  }//ESD
  else if (aodEv)
  {
    nVtxSPD = aodEv->GetNumberOfPileupVerticesSPD();
    nVtxTrk = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD
  
  if(pt < 8)
  {
    fhTimeNPileUpVertSPD  ->Fill(time,nVtxSPD);
    fhTimeNPileUpVertTrack->Fill(time,nVtxTrk);
  }
  
  fhPtPhotonNPileUpSPDVtx->Fill(pt,nVtxSPD);
  fhPtPhotonNPileUpTrkVtx->Fill(pt,nVtxTrk);
	
  if(TMath::Abs(time) < 25)
  {
    fhPtPhotonNPileUpSPDVtxTimeCut->Fill(pt,nVtxSPD);
    fhPtPhotonNPileUpTrkVtxTimeCut->Fill(pt,nVtxTrk);
  }
	
  if(time < 75 && time > -25)
  {
    fhPtPhotonNPileUpSPDVtxTimeCut2->Fill(pt,nVtxSPD);
    fhPtPhotonNPileUpTrkVtxTimeCut2->Fill(pt,nVtxTrk);
  }
  
}

//____________________________________________________________________________________
void  AliAnaPhoton::FillShowerShapeHistograms(AliVCluster* cluster, Int_t mcTag)
{
  //Fill cluster Shower Shape histograms
  
  if(!fFillSSHistograms || GetMixedEvent()) return;
  
  Float_t energy  = cluster->E();
  Int_t   ncells  = cluster->GetNCells();
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  Float_t disp    = cluster->GetDispersion()*cluster->GetDispersion();
  
  TLorentzVector mom;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
  {
    cluster->GetMomentum(mom,GetVertex(0)) ;
  }//Assume that come from vertex in straight line
  else
  {
    Double_t vertex[]={0,0,0};
    cluster->GetMomentum(mom,vertex) ;
  }
  
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  
  fhLam0E ->Fill(energy,lambda0);
  fhLam1E ->Fill(energy,lambda1);
  fhDispE ->Fill(energy,disp);
  
  if(fCalorimeter == "EMCAL" &&  GetFirstSMCoveredByTRD() >= 0 &&
     GetModuleNumber(cluster) >= GetFirstSMCoveredByTRD()  )
  {
    fhLam0ETRD->Fill(energy,lambda0);
    fhLam1ETRD->Fill(energy,lambda1);
    fhDispETRD->Fill(energy,disp);
  }
  
  Float_t l0   = 0., l1   = 0.;
  Float_t dispp= 0., dEta = 0., dPhi    = 0.;
  Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
  if(fCalorimeter == "EMCAL" && !fFillOnlySimpleSSHisto)
  {
    GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                                 l0, l1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);
    //printf("AliAnaPhoton::FillShowerShapeHistogram - l0 %2.6f, l1 %2.6f, disp %2.6f, dEta %2.6f, dPhi %2.6f, sEta %2.6f, sPhi %2.6f, sEtaPhi %2.6f \n",
    //       l0, l1, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi );
    //printf("AliAnaPhoton::FillShowerShapeHistogram - dispersion %f, dispersion eta+phi %f \n",
    //       disp, dPhi+dEta );
    fhDispEtaE        -> Fill(energy,dEta);
    fhDispPhiE        -> Fill(energy,dPhi);
    fhSumEtaE         -> Fill(energy,sEta);
    fhSumPhiE         -> Fill(energy,sPhi);
    fhSumEtaPhiE      -> Fill(energy,sEtaPhi);
    fhDispEtaPhiDiffE -> Fill(energy,dPhi-dEta);
    if(dEta+dPhi>0)fhSphericityE     -> Fill(energy,(dPhi-dEta)/(dEta+dPhi));
    if(dEta+sEta>0)fhDispSumEtaDiffE -> Fill(energy,(dEta-sEta)/((dEta+sEta)/2.));
    if(dPhi+sPhi>0)fhDispSumPhiDiffE -> Fill(energy,(dPhi-sPhi)/((dPhi+sPhi)/2.));
    
    Int_t ebin = -1;
    if      (energy < 2 ) ebin = 0;
    else if (energy < 4 ) ebin = 1;
    else if (energy < 6 ) ebin = 2;
    else if (energy < 10) ebin = 3;
    else if (energy < 15) ebin = 4;
    else if (energy < 20) ebin = 5;
    else                  ebin = 6;
    
    fhDispEtaDispPhi[ebin]->Fill(dEta   ,dPhi);
    fhLambda0DispEta[ebin]->Fill(lambda0,dEta);
    fhLambda0DispPhi[ebin]->Fill(lambda0,dPhi);
    
  }
  
  // if track-matching was of, check effect of track-matching residual cut
  
  if(!fRejectTrackMatch)
  {
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
    {
      dR = 2000., dZ = 2000.;
      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {
      fhLam0ETM ->Fill(energy,lambda0);
      fhLam1ETM ->Fill(energy,lambda1);
      fhDispETM ->Fill(energy,disp);
      
      if(fCalorimeter == "EMCAL" &&  GetFirstSMCoveredByTRD()   >= 0 &&
         GetModuleNumber(cluster) >= GetFirstSMCoveredByTRD()  )
      {
        fhLam0ETMTRD->Fill(energy,lambda0);
        fhLam1ETMTRD->Fill(energy,lambda1);
        fhDispETMTRD->Fill(energy,disp);
      }
    }
  }// if track-matching was of, check effect of matching residual cut
  
  
  if(!fFillOnlySimpleSSHisto)
  {
    if(energy < 2)
    {
      fhNCellsLam0LowE ->Fill(ncells,lambda0);
      fhNCellsLam1LowE ->Fill(ncells,lambda1);
      fhNCellsDispLowE ->Fill(ncells,disp);
      
      fhLam1Lam0LowE  ->Fill(lambda1,lambda0);
      fhLam0DispLowE  ->Fill(lambda0,disp);
      fhDispLam1LowE  ->Fill(disp,lambda1);
      fhEtaLam0LowE   ->Fill(eta,lambda0);
      fhPhiLam0LowE   ->Fill(phi,lambda0);
    }
    else
    {
      fhNCellsLam0HighE ->Fill(ncells,lambda0);
      fhNCellsLam1HighE ->Fill(ncells,lambda1);
      fhNCellsDispHighE ->Fill(ncells,disp);
      
      fhLam1Lam0HighE  ->Fill(lambda1,lambda0);
      fhLam0DispHighE  ->Fill(lambda0,disp);
      fhDispLam1HighE  ->Fill(disp,lambda1);
      fhEtaLam0HighE   ->Fill(eta, lambda0);
      fhPhiLam0HighE   ->Fill(phi, lambda0);
    }
  }
  
  if(IsDataMC())
  {
    AliVCaloCells* cells = 0;
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    //Fill histograms to check shape of embedded clusters
    Float_t fraction = 0;
    // printf("check embedding %i\n",GetReader()->IsEmbeddedClusterSelectionOn());
    
    if(GetReader()->IsEmbeddedClusterSelectionOn())
    {//Only working for EMCAL
      // 	printf("embedded\n");
      Float_t clusterE = 0; // recalculate in case corrections applied.
      Float_t cellE    = 0;
      for(Int_t icell  = 0; icell < cluster->GetNCells(); icell++)
      {
        cellE    = cells->GetCellAmplitude(cluster->GetCellAbsId(icell));
        clusterE+=cellE;
        fraction+=cellE*cluster->GetCellAmplitudeFraction(icell);
      }
      
      //Fraction of total energy due to the embedded signal
      fraction/=clusterE;
      
      if(GetDebug() > 1 )
        printf("AliAnaPhoton::FillShowerShapeHistogram() - Energy fraction of embedded signal %2.3f, Energy %2.3f\n",fraction, clusterE);
      
      fhEmbeddedSignalFractionEnergy->Fill(clusterE,fraction);
      
    }  // embedded fraction
    
    // Get the fraction of the cluster energy that carries the cell with highest energy
    Float_t maxCellFraction = 0.;
    Int_t absID = GetCaloUtils()->GetMaxEnergyCell(cells, cluster,maxCellFraction);
    
    if( absID < 0 ) AliFatal("Wrong absID");
      
    // Check the origin and fill histograms
    
    Int_t mcIndex = -1;
    
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)     &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)        &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta))
    {
      mcIndex = kmcssPhoton ;
      
      if(!GetReader()->IsEmbeddedClusterSelectionOn())
      {
        //Check particle overlaps in cluster
        
        // Compare the primary depositing more energy with the rest,
        // if no photon/electron as comon ancestor (conversions), count as other particle
        const UInt_t nlabels = cluster->GetNLabels();
        Int_t overpdg[nlabels];
        Int_t noverlaps = GetMCAnalysisUtils()->GetNOverlaps(cluster->GetLabels(), nlabels,mcTag,-1,GetReader(),overpdg);

        //printf("N overlaps %d \n",noverlaps);
        
        if(noverlaps == 0)
        {
          fhMCPhotonELambda0NoOverlap  ->Fill(energy, lambda0);
        }
        else if(noverlaps == 1)
        {
          fhMCPhotonELambda0TwoOverlap ->Fill(energy, lambda0);
        }
        else if(noverlaps > 1)
        {
          fhMCPhotonELambda0NOverlap   ->Fill(energy, lambda0);
        }
        else
        {
          printf("AliAnaPhoton::FillShowerShapeHistogram() - n overlaps = %d!!\n", noverlaps);
        }
      }//No embedding
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        if     (fraction > 0.9)
        {
          fhEmbedPhotonELambda0FullSignal   ->Fill(energy, lambda0);
        }
        else if(fraction > 0.5)
        {
          fhEmbedPhotonELambda0MostlySignal ->Fill(energy, lambda0);
        }
        else if(fraction > 0.1)
        {
          fhEmbedPhotonELambda0MostlyBkg    ->Fill(energy, lambda0);
        }
        else
        {
          fhEmbedPhotonELambda0FullBkg      ->Fill(energy, lambda0);
        }
      } // embedded
      
    }//photon   no conversion
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton)     &&
               GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
              !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)        &&
              !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta))
    {
      mcIndex = kmcssConversion ;
    }//conversion photon

    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron))
    {
      mcIndex = kmcssElectron ;
    }//electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  )
    {
      mcIndex = kmcssPi0 ;
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        if     (fraction > 0.9)
        {
          fhEmbedPi0ELambda0FullSignal   ->Fill(energy, lambda0);
        }
        else if(fraction > 0.5)
        {
          fhEmbedPi0ELambda0MostlySignal ->Fill(energy, lambda0);
        }
        else if(fraction > 0.1)
        {
          fhEmbedPi0ELambda0MostlyBkg    ->Fill(energy, lambda0);
        }
        else
        {
          fhEmbedPi0ELambda0FullBkg      ->Fill(energy, lambda0);
        }
      } // embedded
      
    }//pi0
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)  )
    {
      mcIndex = kmcssEta ;
    }//eta
    else
    {
      mcIndex = kmcssOther ;
    }//other particles
    
    fhMCELambda0           [mcIndex]->Fill(energy, lambda0);
    fhMCELambda1           [mcIndex]->Fill(energy, lambda1);
    fhMCEDispersion        [mcIndex]->Fill(energy, disp);
    fhMCNCellsE            [mcIndex]->Fill(energy, ncells);
    fhMCMaxCellDiffClusterE[mcIndex]->Fill(energy, maxCellFraction);
    
    if(!fFillOnlySimpleSSHisto)
    {
      if     (energy < 2.)
      {
        fhMCLambda0vsClusterMaxCellDiffE0[mcIndex]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0 [mcIndex]->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.)
      {
        fhMCLambda0vsClusterMaxCellDiffE2[mcIndex]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2 [mcIndex]->Fill(ncells,  maxCellFraction);
      }
      else
      {
        fhMCLambda0vsClusterMaxCellDiffE6[mcIndex]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6 [mcIndex]->Fill(ncells,  maxCellFraction);
      }
      
      if(fCalorimeter == "EMCAL")
      {
        fhMCEDispEta        [mcIndex]-> Fill(energy,dEta);
        fhMCEDispPhi        [mcIndex]-> Fill(energy,dPhi);
        fhMCESumEtaPhi      [mcIndex]-> Fill(energy,sEtaPhi);
        fhMCEDispEtaPhiDiff [mcIndex]-> Fill(energy,dPhi-dEta);
        if(dEta+dPhi>0)fhMCESphericity[mcIndex]-> Fill(energy,(dPhi-dEta)/(dEta+dPhi));
        
        Int_t ebin = -1;
        if      (energy < 2 ) ebin = 0;
        else if (energy < 4 ) ebin = 1;
        else if (energy < 6 ) ebin = 2;
        else if (energy < 10) ebin = 3;
        else if (energy < 15) ebin = 4;
        else if (energy < 20) ebin = 5;
        else                  ebin = 6;
        
        fhMCDispEtaDispPhi[ebin][mcIndex]->Fill(dEta   ,dPhi);
        fhMCLambda0DispEta[ebin][mcIndex]->Fill(lambda0,dEta);
        fhMCLambda0DispPhi[ebin][mcIndex]->Fill(lambda0,dPhi);
      }
    }
  }//MC data
  
}

//__________________________________________________________________________
void AliAnaPhoton::FillTrackMatchingResidualHistograms(AliVCluster* cluster,
                                                       Int_t cut)
{
  // If selected, fill histograms with residuals of matched clusters, help to define track matching cut
  // Residual filled for different cuts 0 (No cut), after 1 PID cut
  
  Float_t dZ = cluster->GetTrackDz();
  Float_t dR = cluster->GetTrackDx();
  
  if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
  {
    dR = 2000., dZ = 2000.;
    GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
  }
  
  AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
  
  Bool_t positive = kFALSE;
  if(track) positive = (track->Charge()>0);
  
  if(fhTrackMatchedDEta[cut] && TMath::Abs(dR) < 999)
  {
    fhTrackMatchedDEta[cut]->Fill(cluster->E(),dZ);
    fhTrackMatchedDPhi[cut]->Fill(cluster->E(),dR);
    if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhi[cut]->Fill(dZ,dR);

    if(track)
    {
      if(positive)
      {
        fhTrackMatchedDEtaPos[cut]->Fill(cluster->E(),dZ);
        fhTrackMatchedDPhiPos[cut]->Fill(cluster->E(),dR);
        if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhiPos[cut]->Fill(dZ,dR);
      }
      else
      {
        fhTrackMatchedDEtaNeg[cut]->Fill(cluster->E(),dZ);
        fhTrackMatchedDPhiNeg[cut]->Fill(cluster->E(),dR);
        if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhiNeg[cut]->Fill(dZ,dR);
      }
    }
    
    Int_t nSMod = GetModuleNumber(cluster);
    
    if(fCalorimeter=="EMCAL" &&   GetFirstSMCoveredByTRD() >= 0 &&
       nSMod >= GetFirstSMCoveredByTRD()   )
    {
      fhTrackMatchedDEtaTRD[cut]->Fill(cluster->E(),dZ);
      fhTrackMatchedDPhiTRD[cut]->Fill(cluster->E(),dR);
    }
    
    // Check dEdx and E/p of matched clusters

    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {
      if(track)
      {
        Float_t dEdx   = track->GetTPCsignal();
        Float_t eOverp = cluster->E()/track->P();
        
        fhdEdx[cut]  ->Fill(cluster->E(), dEdx);
        fhEOverP[cut]->Fill(cluster->E(), eOverp);
        
        if(fCalorimeter=="EMCAL" &&  GetFirstSMCoveredByTRD() >= 0 &&
           nSMod >= GetFirstSMCoveredByTRD()  )
          fhEOverPTRD[cut]->Fill(cluster->E(), eOverp);
        
        
      }
      else
        printf("AliAnaPhoton::FillTrackMatchingResidualHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
      
      
      
      if(IsDataMC())
      {
        
        Int_t tag = GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(),GetReader(),fCalorimeter);
        
        if  ( !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)  )
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)      ||
                    GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 2.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 0.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 1.5 );
          else                                                                                 fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 3.5 );
          
          // Check if several particles contributed to cluster and discard overlapped mesons
          if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) ||
             !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta))
          {
            if(cluster->GetNLabels()==1)
            {
              fhTrackMatchedDEtaMCNoOverlap[cut]->Fill(cluster->E(),dZ);
              fhTrackMatchedDPhiMCNoOverlap[cut]->Fill(cluster->E(),dR);
            }
            else
            {
              fhTrackMatchedDEtaMCOverlap[cut]->Fill(cluster->E(),dZ);
              fhTrackMatchedDPhiMCOverlap[cut]->Fill(cluster->E(),dR);
            }
            
          }// Check overlaps
          
        }
        else
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)      ||
                    GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 6.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 4.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 5.5 );
          else                                                                                 fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 7.5 );
          
          // Check if several particles contributed to cluster
          if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) ||
             !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta))
          {
            fhTrackMatchedDEtaMCConversion[cut]->Fill(cluster->E(),dZ);
            fhTrackMatchedDPhiMCConversion[cut]->Fill(cluster->E(),dR);
            
          }// Check overlaps
          
        }
        
      } // MC
      
    } // residuals window
    
  } // Small residual
  
}

//___________________________________________
TObjString *  AliAnaPhoton::GetAnalysisCuts()
{
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaPhoton ---\n") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist =%2.2f (Minimal distance to bad channel to accept cluster) \n",fMinDist) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist2=%2.2f (Cuts on Minimal distance to study acceptance evaluation) \n",fMinDist2) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinDist3=%2.2f (One more cut on distance used for acceptance-efficiency study) \n",fMinDist3) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRejectTrackMatch: %d\n",fRejectTrackMatch) ;
  parList+=onePar ;
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in PID class.
  parList += GetCaloPID()->GetPIDParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList()
  
  return new TObjString(parList) ;
}

//________________________________________________________________________
TList *  AliAnaPhoton::GetCreateOutputObjects()
{
  // Create histograms to be saved in output file and
  // store them in outputContainer
  TList * outputContainer = new TList() ;
  outputContainer->SetName("PhotonHistos") ;
	
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Float_t phimax = GetHistogramRanges()->GetHistoPhiMax(); Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins(); Float_t etamax = GetHistogramRanges()->GetHistoEtaMax(); Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  Int_t ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax   = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin   = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t nbins    = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   nmax    = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   nmin    = GetHistogramRanges()->GetHistoNClusterCellMin();
  Int_t ntimebins= GetHistogramRanges()->GetHistoTimeBins();         Float_t timemax = GetHistogramRanges()->GetHistoTimeMax();         Float_t timemin = GetHistogramRanges()->GetHistoTimeMin();
  
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
  
  Int_t bin[] = {0,2,4,6,10,15,20,100}; // energy bins for SS studies
  
  TString cut[] = {"Open","Reader","E","Time","NCells","NLM","Fidutial","Matching","Bad","PID"};
  for (Int_t i = 0; i < 10 ;  i++)
  {
    fhClusterCutsE[i] = new TH1F(Form("hE_Cut_%d_%s", i, cut[i].Data()),
                                Form("Number of clusters that pass cuts <= %d, %s", i, cut[i].Data()),
                                nptbins,ptmin,ptmax);
    fhClusterCutsE[i]->SetYTitle("d#it{N}/d#it{E} ");
    fhClusterCutsE[i]->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhClusterCutsE[i]) ;
    
    fhClusterCutsPt[i] = new TH1F(Form("hPt_Cut_%d_%s", i, cut[i].Data()),
                                Form("Number of clusters that pass cuts <= %d, %s", i, cut[i].Data()),
                                nptbins,ptmin,ptmax);
    fhClusterCutsPt[i]->SetYTitle("d#it{N}/d#it{E} ");
    fhClusterCutsPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhClusterCutsPt[i]) ;
  }
  
  fhEClusterSM = new TH2F("hEClusterSM","Raw clusters E and super-module number",
                          nptbins,ptmin,ptmax,
                          GetCaloUtils()->GetNumberOfSuperModulesUsed(),0,GetCaloUtils()->GetNumberOfSuperModulesUsed());
  fhEClusterSM->SetYTitle("SuperModule ");
  fhEClusterSM->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhEClusterSM) ;

  fhPtClusterSM = new TH2F("hPtClusterSM","Raw clusters #it{p}_{T} and super-module number",
                          nptbins,ptmin,ptmax,
                          GetCaloUtils()->GetNumberOfSuperModulesUsed(),0,GetCaloUtils()->GetNumberOfSuperModulesUsed());
  fhPtClusterSM->SetYTitle("SuperModule ");
  fhPtClusterSM->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhPtClusterSM) ;
  
  fhEPhotonSM = new TH2F("hEPhotonSM","Selected clusters E and super-module number",
                          nptbins,ptmin,ptmax,
                          GetCaloUtils()->GetNumberOfSuperModulesUsed(),0,GetCaloUtils()->GetNumberOfSuperModulesUsed());
  fhEPhotonSM->SetYTitle("SuperModule ");
  fhEPhotonSM->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhEPhotonSM) ;
  
  fhPtPhotonSM = new TH2F("hPtPhotonSM","Selected clusters #it{p}_{T} and super-module number",
                           nptbins,ptmin,ptmax,
                           GetCaloUtils()->GetNumberOfSuperModulesUsed(),0,GetCaloUtils()->GetNumberOfSuperModulesUsed());
  fhPtPhotonSM->SetYTitle("SuperModule ");
  fhPtPhotonSM->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhPtPhotonSM) ;
  
  fhNCellsE  = new TH2F ("hNCellsE","# of cells in cluster vs E of clusters", nptbins,ptmin,ptmax, nbins,nmin,nmax);
  fhNCellsE->SetXTitle("#it{E} (GeV)");
  fhNCellsE->SetYTitle("# of cells in cluster");
  outputContainer->Add(fhNCellsE);
  
  fhCellsE  = new TH2F ("hCellsE","energy of cells in cluster vs E of clusters", nptbins,ptmin,ptmax, nptbins*2,ptmin,ptmax);
  fhCellsE->SetXTitle("#it{E}_{cluster} (GeV)");
  fhCellsE->SetYTitle("#it{E}_{cell} (GeV)");
  outputContainer->Add(fhCellsE);
  
  fhTimePt  = new TH2F ("hTimePt","time of cluster vs pT of clusters", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimePt->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhTimePt->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimePt);
  
  fhMaxCellDiffClusterE  = new TH2F ("hMaxCellDiffClusterE","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                     nptbins,ptmin,ptmax, 500,0,1.);
  fhMaxCellDiffClusterE->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhMaxCellDiffClusterE->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
  outputContainer->Add(fhMaxCellDiffClusterE);
  
  fhEPhoton  = new TH1F("hEPhoton","Number of #gamma over calorimeter vs energy",nptbins,ptmin,ptmax);
  fhEPhoton->SetYTitle("#it{counts}");
  fhEPhoton->SetXTitle("#it{E}_{#gamma}(GeV)");
  outputContainer->Add(fhEPhoton) ;
  
  fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter vs #it{p}_{T}",nptbins,ptmin,ptmax);
  fhPtPhoton->SetYTitle("#it{counts}");
  fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/#it{c})");
  outputContainer->Add(fhPtPhoton) ;
  
  fhPtCentralityPhoton  = new TH2F("hPtCentralityPhoton","centrality vs #it{p}_{T}",nptbins,ptmin,ptmax, 100,0,100);
  fhPtCentralityPhoton->SetYTitle("Centrality");
  fhPtCentralityPhoton->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  outputContainer->Add(fhPtCentralityPhoton) ;
  
  fhPtEventPlanePhoton  = new TH2F("hPtEventPlanePhoton","centrality vs #it{p}_{T}",nptbins,ptmin,ptmax, 100,0,TMath::Pi());
  fhPtEventPlanePhoton->SetYTitle("Event plane angle (rad)");
  fhPtEventPlanePhoton->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtEventPlanePhoton) ;
  
  fhEtaPhi  = new TH2F
  ("hEtaPhi","cluster,#it{E} > 0.5 GeV, #eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhi->SetYTitle("#phi (rad)");
  fhEtaPhi->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhi) ;
  
  fhPhiPhoton  = new TH2F
  ("hPhiPhoton","#phi_{#gamma} vs #it{p}_{T}",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhPhiPhoton->SetYTitle("#phi (rad)");
  fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/#it{c})");
  outputContainer->Add(fhPhiPhoton) ;
  
  fhEtaPhoton  = new TH2F
  ("hEtaPhoton","#eta_{#gamma} vs #it{p}_{T}",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhEtaPhoton->SetYTitle("#eta");
  fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/#it{c})");
  outputContainer->Add(fhEtaPhoton) ;
  
  fhEtaPhiPhoton  = new TH2F
  ("hEtaPhiPhoton","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiPhoton->SetYTitle("#phi (rad)");
  fhEtaPhiPhoton->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiPhoton) ;
  if(GetMinPt() < 0.5)
  {
    fhEtaPhi05Photon  = new TH2F
    ("hEtaPhi05Photon","#eta vs #phi, E < 0.5",netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhi05Photon->SetYTitle("#phi (rad)");
    fhEtaPhi05Photon->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05Photon) ;
  }
  
  fhNLocMax = new TH2F("hNLocMax","Number of local maxima in cluster",
                       nptbins,ptmin,ptmax,10,0,10);
  fhNLocMax ->SetYTitle("N maxima");
  fhNLocMax ->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhNLocMax) ;
  
  //Shower shape
  if(fFillSSHistograms)
  {
    fhLam0E  = new TH2F ("hLam0E","#lambda_{0}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhLam0E->SetYTitle("#lambda_{0}^{2}");
    fhLam0E->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhLam0E);
    
    fhLam1E  = new TH2F ("hLam1E","#lambda_{1}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhLam1E->SetYTitle("#lambda_{1}^{2}");
    fhLam1E->SetXTitle("#it{E} (GeV)");
    outputContainer->Add(fhLam1E);
    
    fhDispE  = new TH2F ("hDispE"," dispersion^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
    fhDispE->SetYTitle("D^{2}");
    fhDispE->SetXTitle("#it{E} (GeV) ");
    outputContainer->Add(fhDispE);
    
    if(!fRejectTrackMatch)
    {
      fhLam0ETM  = new TH2F ("hLam0ETM","#lambda_{0}^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLam0ETM->SetYTitle("#lambda_{0}^{2}");
      fhLam0ETM->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam0ETM);
      
      fhLam1ETM  = new TH2F ("hLam1ETM","#lambda_{1}^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLam1ETM->SetYTitle("#lambda_{1}^{2}");
      fhLam1ETM->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam1ETM);
      
      fhDispETM  = new TH2F ("hDispETM"," dispersion^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhDispETM->SetYTitle("D^{2}");
      fhDispETM->SetXTitle("#it{E} (GeV) ");
      outputContainer->Add(fhDispETM);
    }
    
    if(fCalorimeter == "EMCAL" &&  GetFirstSMCoveredByTRD() >= 0)
    {
      fhLam0ETRD  = new TH2F ("hLam0ETRD","#lambda_{0}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLam0ETRD->SetYTitle("#lambda_{0}^{2}");
      fhLam0ETRD->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam0ETRD);
      
      fhLam1ETRD  = new TH2F ("hLam1ETRD","#lambda_{1}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhLam1ETRD->SetYTitle("#lambda_{1}^{2}");
      fhLam1ETRD->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhLam1ETRD);
      
      fhDispETRD  = new TH2F ("hDispETRD"," dispersion^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhDispETRD->SetYTitle("Dispersion^{2}");
      fhDispETRD->SetXTitle("#it{E} (GeV) ");
      outputContainer->Add(fhDispETRD);
      
      if(!fRejectTrackMatch &&  GetFirstSMCoveredByTRD() >=0 )
      {
        fhLam0ETMTRD  = new TH2F ("hLam0ETMTRD","#lambda_{0}^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhLam0ETMTRD->SetYTitle("#lambda_{0}^{2}");
        fhLam0ETMTRD->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhLam0ETMTRD);
        
        fhLam1ETMTRD  = new TH2F ("hLam1ETMTRD","#lambda_{1}^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhLam1ETMTRD->SetYTitle("#lambda_{1}^{2}");
        fhLam1ETMTRD->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhLam1ETMTRD);
        
        fhDispETMTRD  = new TH2F ("hDispETMTRD"," dispersion^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhDispETMTRD->SetYTitle("Dispersion^{2}");
        fhDispETMTRD->SetXTitle("#it{E} (GeV) ");
        outputContainer->Add(fhDispETMTRD);
      }
    }
    
    if(!fFillOnlySimpleSSHisto)
    {
      fhNCellsLam0LowE  = new TH2F ("hNCellsLam0LowE","N_{cells} in cluster vs #lambda_{0}^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsLam0LowE->SetXTitle("N_{cells}");
      fhNCellsLam0LowE->SetYTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsLam0LowE);
      
      fhNCellsLam0HighE  = new TH2F ("hNCellsLam0HighE","N_{cells} in cluster vs #lambda_{0}^{2}, #it{E} > 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsLam0HighE->SetXTitle("N_{cells}");
      fhNCellsLam0HighE->SetYTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsLam0HighE);
      
      fhNCellsLam1LowE  = new TH2F ("hNCellsLam1LowE","N_{cells} in cluster vs #lambda_{1}^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsLam1LowE->SetXTitle("N_{cells}");
      fhNCellsLam1LowE->SetYTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsLam1LowE);
      
      fhNCellsLam1HighE  = new TH2F ("hNCellsLam1HighE","N_{cells} in cluster vs #lambda_{1}^{2}, #it{E} > 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsLam1HighE->SetXTitle("N_{cells}");
      fhNCellsLam1HighE->SetYTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhNCellsLam1HighE);
      
      fhNCellsDispLowE  = new TH2F ("hNCellsDispLowE","N_{cells} in cluster vs dispersion^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsDispLowE->SetXTitle("N_{cells}");
      fhNCellsDispLowE->SetYTitle("D^{2}");
      outputContainer->Add(fhNCellsDispLowE);
      
      fhNCellsDispHighE  = new TH2F ("hNCellsDispHighE","N_{cells} in cluster vs dispersion^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax);
      fhNCellsDispHighE->SetXTitle("N_{cells}");
      fhNCellsDispHighE->SetYTitle("D^{2}");
      outputContainer->Add(fhNCellsDispHighE);
      
      fhEtaLam0LowE  = new TH2F ("hEtaLam0LowE","#eta vs #lambda_{0}^{2}, E < 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax);
      fhEtaLam0LowE->SetYTitle("#lambda_{0}^{2}");
      fhEtaLam0LowE->SetXTitle("#eta");
      outputContainer->Add(fhEtaLam0LowE);
      
      fhPhiLam0LowE  = new TH2F ("hPhiLam0LowE","#phi vs #lambda_{0}^{2}, E < 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax);
      fhPhiLam0LowE->SetYTitle("#lambda_{0}^{2}");
      fhPhiLam0LowE->SetXTitle("#phi");
      outputContainer->Add(fhPhiLam0LowE);
      
      fhEtaLam0HighE  = new TH2F ("hEtaLam0HighE","#eta vs #lambda_{0}^{2}, #it{E} > 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax);
      fhEtaLam0HighE->SetYTitle("#lambda_{0}^{2}");
      fhEtaLam0HighE->SetXTitle("#eta");
      outputContainer->Add(fhEtaLam0HighE);
      
      fhPhiLam0HighE  = new TH2F ("hPhiLam0HighE","#phi vs #lambda_{0}^{2}, #it{E} > 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax);
      fhPhiLam0HighE->SetYTitle("#lambda_{0}^{2}");
      fhPhiLam0HighE->SetXTitle("#phi");
      outputContainer->Add(fhPhiLam0HighE);
      
      fhLam1Lam0LowE  = new TH2F ("hLam1Lam0LowE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhLam1Lam0LowE->SetYTitle("#lambda_{0}^{2}");
      fhLam1Lam0LowE->SetXTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhLam1Lam0LowE);
      
      fhLam1Lam0HighE  = new TH2F ("hLam1Lam0HighE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of #it{E} > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhLam1Lam0HighE->SetYTitle("#lambda_{0}^{2}");
      fhLam1Lam0HighE->SetXTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhLam1Lam0HighE);
      
      fhLam0DispLowE  = new TH2F ("hLam0DispLowE","#lambda_{0}^{2} vs dispersion^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhLam0DispLowE->SetXTitle("#lambda_{0}^{2}");
      fhLam0DispLowE->SetYTitle("D^{2}");
      outputContainer->Add(fhLam0DispLowE);
      
      fhLam0DispHighE  = new TH2F ("hLam0DispHighE","#lambda_{0}^{2} vs dispersion^{2} in cluster of #it{E} > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhLam0DispHighE->SetXTitle("#lambda_{0}^{2}");
      fhLam0DispHighE->SetYTitle("D^{2}");
      outputContainer->Add(fhLam0DispHighE);
      
      fhDispLam1LowE  = new TH2F ("hDispLam1LowE","Dispersion^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhDispLam1LowE->SetXTitle("D^{2}");
      fhDispLam1LowE->SetYTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhDispLam1LowE);
      
      fhDispLam1HighE  = new TH2F ("hDispLam1HighE","Dispersion^{2} vs #lambda_{1^{2}} in cluster of #it{E} > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax);
      fhDispLam1HighE->SetXTitle("D^{2}");
      fhDispLam1HighE->SetYTitle("#lambda_{1}^{2}");
      outputContainer->Add(fhDispLam1HighE);
      
      if(fCalorimeter == "EMCAL")
      {
        fhDispEtaE  = new TH2F ("hDispEtaE","#sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhDispEtaE->SetXTitle("#it{E} (GeV)");
        fhDispEtaE->SetYTitle("#sigma^{2}_{#eta #eta}");
        outputContainer->Add(fhDispEtaE);
        
        fhDispPhiE  = new TH2F ("hDispPhiE","#sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhDispPhiE->SetXTitle("#it{E} (GeV)");
        fhDispPhiE->SetYTitle("#sigma^{2}_{#phi #phi}");
        outputContainer->Add(fhDispPhiE);
        
        fhSumEtaE  = new TH2F ("hSumEtaE","#delta^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i})^{2} / #Sigma w_{i} - <#eta>^{2} vs E",  nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhSumEtaE->SetXTitle("#it{E} (GeV)");
        fhSumEtaE->SetYTitle("#delta^{2}_{#eta #eta}");
        outputContainer->Add(fhSumEtaE);
        
        fhSumPhiE  = new TH2F ("hSumPhiE","#delta^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i})^{2}/ #Sigma w_{i} - <#phi>^{2} vs E",
                               nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
        fhSumPhiE->SetXTitle("#it{E} (GeV)");
        fhSumPhiE->SetYTitle("#delta^{2}_{#phi #phi}");
        outputContainer->Add(fhSumPhiE);
        
        fhSumEtaPhiE  = new TH2F ("hSumEtaPhiE","#delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",
                                  nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
        fhSumEtaPhiE->SetXTitle("#it{E} (GeV)");
        fhSumEtaPhiE->SetYTitle("#delta^{2}_{#eta #phi}");
        outputContainer->Add(fhSumEtaPhiE);
        
        fhDispEtaPhiDiffE  = new TH2F ("hDispEtaPhiDiffE","#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",
                                       nptbins,ptmin,ptmax,200, -10,10);
        fhDispEtaPhiDiffE->SetXTitle("#it{E} (GeV)");
        fhDispEtaPhiDiffE->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
        outputContainer->Add(fhDispEtaPhiDiffE);
        
        fhSphericityE  = new TH2F ("hSphericityE","(#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",
                                   nptbins,ptmin,ptmax, 200, -1,1);
        fhSphericityE->SetXTitle("#it{E} (GeV)");
        fhSphericityE->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
        outputContainer->Add(fhSphericityE);
        
        fhDispSumEtaDiffE  = new TH2F ("hDispSumEtaDiffE","#sigma^{2}_{#eta #eta} - #delta^{2}_{#eta #eta} / average vs E",  nptbins,ptmin,ptmax, 200,-0.01,0.01);
        fhDispSumEtaDiffE->SetXTitle("#it{E} (GeV)");
        fhDispSumEtaDiffE->SetYTitle("#sigma^{2}_{#eta #eta} - #delta^{2}_{#eta #eta} / average");
        outputContainer->Add(fhDispSumEtaDiffE);
        
        fhDispSumPhiDiffE  = new TH2F ("hDispSumPhiDiffE","#sigma^{2}_{#phi #phi} - #delta^{2}_{#phi #phi} / average vs E",  nptbins,ptmin,ptmax, 200,-0.01,0.01);
        fhDispSumPhiDiffE->SetXTitle("#it{E} (GeV)");
        fhDispSumPhiDiffE->SetYTitle("#sigma^{2}_{#phi #phi} - #delta^{2}_{#phi #phi} / average");
        outputContainer->Add(fhDispSumPhiDiffE);
        
        for(Int_t i = 0; i < 7; i++)
        {
          fhDispEtaDispPhi[i] = new TH2F (Form("hDispEtaDispPhi_EBin%d",i),Form("#sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",bin[i],bin[i+1]),
                                          ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
          fhDispEtaDispPhi[i]->SetXTitle("#sigma^{2}_{#eta #eta}");
          fhDispEtaDispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
          outputContainer->Add(fhDispEtaDispPhi[i]);
          
          fhLambda0DispEta[i] = new TH2F (Form("hLambda0DispEta_EBin%d",i),Form("#lambda^{2}_{0} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",bin[i],bin[i+1]),
                                          ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
          fhLambda0DispEta[i]->SetXTitle("#lambda^{2}_{0}");
          fhLambda0DispEta[i]->SetYTitle("#sigma^{2}_{#eta #eta}");
          outputContainer->Add(fhLambda0DispEta[i]);
          
          fhLambda0DispPhi[i] = new TH2F (Form("hLambda0DispPhi_EBin%d",i),Form("#lambda^{2}_{0}} vs #sigma^{2}_{#phi #phi} for %d < E < %d GeV",bin[i],bin[i+1]),
                                          ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
          fhLambda0DispPhi[i]->SetXTitle("#lambda^{2}_{0}");
          fhLambda0DispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
          outputContainer->Add(fhLambda0DispPhi[i]);
        }
      }
    }
  } // Shower shape
  
  // Track Matching
  
  if(fFillTMHisto)
  {
    TString cutTM [] = {"NoCut",""};
    
    for(Int_t i = 0; i < 2; i++)
    {
      fhTrackMatchedDEta[i]  = new TH2F
      (Form("hTrackMatchedDEta%s",cutTM[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEta[i]->SetYTitle("d#eta");
      fhTrackMatchedDEta[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhi[i]  = new TH2F
      (Form("hTrackMatchedDPhi%s",cutTM[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhi[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhi[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhi[i]  = new TH2F
      (Form("hTrackMatchedDEtaDPhi%s",cutTM[i].Data()),
       Form("d#eta vs d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhi[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDEtaDPhi[i]->SetXTitle("d#eta");
      
      fhTrackMatchedDEtaPos[i]  = new TH2F
      (Form("hTrackMatchedDEtaPos%s",cutTM[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaPos[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaPos[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiPos[i]  = new TH2F
      (Form("hTrackMatchedDPhiPos%s",cutTM[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiPos[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiPos[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiPos[i]  = new TH2F
      (Form("hTrackMatchedDEtaDPhiPos%s",cutTM[i].Data()),
       Form("d#eta vs d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiPos[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDEtaDPhiPos[i]->SetXTitle("d#eta");
      
      fhTrackMatchedDEtaNeg[i]  = new TH2F
      (Form("hTrackMatchedDEtaNeg%s",cutTM[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
      fhTrackMatchedDEtaNeg[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNeg[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNeg[i]  = new TH2F
      (Form("hTrackMatchedDPhiNeg%s",cutTM[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDPhiNeg[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNeg[i]->SetXTitle("#it{E}_{cluster} (GeV)");
      
      fhTrackMatchedDEtaDPhiNeg[i]  = new TH2F
      (Form("hTrackMatchedDEtaDPhiNeg%s",cutTM[i].Data()),
       Form("d#eta vs d#phi of cluster-track vs cluster energy, %s",cutTM[i].Data()),
       nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
      fhTrackMatchedDEtaDPhiNeg[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDEtaDPhiNeg[i]->SetXTitle("d#eta");
      
      fhdEdx[i]  = new TH2F (Form("hdEdx%s",cutTM[i].Data()),Form("matched track <dE/dx> vs cluster E, %s",cutTM[i].Data()),
                             nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax);
      fhdEdx[i]->SetXTitle("#it{E} (GeV)");
      fhdEdx[i]->SetYTitle("<dE/dx>");
      
      fhEOverP[i]  = new TH2F (Form("hEOverP%s",cutTM[i].Data()),Form("matched track E/p vs cluster E, %s",cutTM[i].Data()),
                               nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
      fhEOverP[i]->SetXTitle("#it{E} (GeV)");
      fhEOverP[i]->SetYTitle("E/p");
      
      outputContainer->Add(fhTrackMatchedDEta[i]) ;
      outputContainer->Add(fhTrackMatchedDPhi[i]) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhi[i]) ;
      outputContainer->Add(fhTrackMatchedDEtaPos[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiPos[i]) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhiPos[i]) ;
      outputContainer->Add(fhTrackMatchedDEtaNeg[i]) ;
      outputContainer->Add(fhTrackMatchedDPhiNeg[i]) ;
      outputContainer->Add(fhTrackMatchedDEtaDPhiNeg[i]) ;
      outputContainer->Add(fhdEdx[i]);
      outputContainer->Add(fhEOverP[i]);
      
      if(fCalorimeter=="EMCAL" &&  GetFirstSMCoveredByTRD() >=0 )
      {
        fhTrackMatchedDEtaTRD[i]  = new TH2F
        (Form("hTrackMatchedDEtaTRD%s",cutTM[i].Data()),
         Form("d#eta of cluster-track vs cluster energy, SM behind TRD, %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEtaTRD[i]->SetYTitle("d#eta");
        fhTrackMatchedDEtaTRD[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        fhTrackMatchedDPhiTRD[i]  = new TH2F
        (Form("hTrackMatchedDPhiTRD%s",cutTM[i].Data()),
         Form("d#phi of cluster-track vs cluster energy, SM behing TRD, %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhiTRD[i]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhiTRD[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        fhEOverPTRD[i]  = new TH2F
        (Form("hEOverPTRD%s",cutTM[i].Data()),
         Form("matched track E/p vs cluster E, behind TRD, %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax);
        fhEOverPTRD[i]->SetXTitle("#it{E} (GeV)");
        fhEOverPTRD[i]->SetYTitle("E/p");
        
        outputContainer->Add(fhTrackMatchedDEtaTRD[i]) ;
        outputContainer->Add(fhTrackMatchedDPhiTRD[i]) ;
        outputContainer->Add(fhEOverPTRD[i]);
      }
      
      if(IsDataMC())
      {
        fhTrackMatchedDEtaMCNoOverlap[i]  = new TH2F
        (Form("hTrackMatchedDEtaMCNoOverlap%s",cutTM[i].Data()),
         Form("d#eta of cluster-track vs cluster energy, no other MC particles overlap %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEtaMCNoOverlap[i]->SetYTitle("d#eta");
        fhTrackMatchedDEtaMCNoOverlap[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        fhTrackMatchedDPhiMCNoOverlap[i]  = new TH2F
        (Form("hTrackMatchedDPhiMCNoOverlap%s",cutTM[i].Data()),
         Form("d#phi of cluster-track vs cluster energy, no other MC particles overlap %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhiMCNoOverlap[i]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhiMCNoOverlap[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        outputContainer->Add(fhTrackMatchedDEtaMCNoOverlap[i]) ;
        outputContainer->Add(fhTrackMatchedDPhiMCNoOverlap[i]) ;
        fhTrackMatchedDEtaMCOverlap[i]  = new TH2F
        (Form("hTrackMatchedDEtaMCOverlap%s",cutTM[i].Data()),
         Form("d#eta of cluster-track vs cluster energy, several MC particles overlap %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEtaMCOverlap[i]->SetYTitle("d#eta");
        fhTrackMatchedDEtaMCOverlap[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        fhTrackMatchedDPhiMCOverlap[i]  = new TH2F
        (Form("hTrackMatchedDPhiMCOverlap%s",cutTM[i].Data()),
         Form("d#phi of cluster-track vs cluster energy, several MC particles overlap %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhiMCOverlap[i]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhiMCOverlap[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        outputContainer->Add(fhTrackMatchedDEtaMCOverlap[i]) ;
        outputContainer->Add(fhTrackMatchedDPhiMCOverlap[i]) ;
        
        fhTrackMatchedDEtaMCConversion[i]  = new TH2F
        (Form("hTrackMatchedDEtaMCConversion%s",cutTM[i].Data()),
         Form("d#eta of cluster-track vs cluster energy, no other MC particles overlap appart from conversions %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
        fhTrackMatchedDEtaMCConversion[i]->SetYTitle("d#eta");
        fhTrackMatchedDEtaMCConversion[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        fhTrackMatchedDPhiMCConversion[i]  = new TH2F
        (Form("hTrackMatchedDPhiMCConversion%s",cutTM[i].Data()),
         Form("d#phi of cluster-track vs cluster energy, no other MC particles overlap appart from conversions %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
        fhTrackMatchedDPhiMCConversion[i]->SetYTitle("d#phi (rad)");
        fhTrackMatchedDPhiMCConversion[i]->SetXTitle("#it{E}_{cluster} (GeV)");
        
        outputContainer->Add(fhTrackMatchedDEtaMCConversion[i]) ;
        outputContainer->Add(fhTrackMatchedDPhiMCConversion[i]) ;
        
        fhTrackMatchedMCParticle[i]  = new TH2F
        (Form("hTrackMatchedMCParticle%s",cutTM[i].Data()),
         Form("Origin of particle vs energy %s",cutTM[i].Data()),
         nptbins,ptmin,ptmax,8,0,8);
        fhTrackMatchedMCParticle[i]->SetXTitle("#it{E} (GeV)");
        //fhTrackMatchedMCParticle[i]->SetYTitle("Particle type");
        
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(1 ,"Photon");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(2 ,"Electron");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(4 ,"Rest");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
        fhTrackMatchedMCParticle[i]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
        
        outputContainer->Add(fhTrackMatchedMCParticle[i]);
      }
    }
  }
  
  if(fFillPileUpHistograms)
  {
    TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
    
    for(Int_t i = 0 ; i < 7 ; i++)
    {
      fhPtPhotonPileUp[i]  = new TH1F(Form("hPtPhotonPileUp%s",pileUpName[i].Data()),
                                      Form("Selected photon #it{p}_{T} distribution, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptmax);
      fhPtPhotonPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPhotonPileUp[i]);
      
      fhClusterTimeDiffPhotonPileUp[i]  = new TH2F(Form("hClusterTimeDiffPhotonPileUp%s",pileUpName[i].Data()),
                                             Form("Photon cluster E vs #it{t}_{max}-#it{t}_{cell} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                             nptbins,ptmin,ptmax,400,-200,200);
      fhClusterTimeDiffPhotonPileUp[i]->SetXTitle("#it{E} (GeV)");
      fhClusterTimeDiffPhotonPileUp[i]->SetYTitle("#it{t}_{max}-#it{t}_{cell} (ns)");
      outputContainer->Add(fhClusterTimeDiffPhotonPileUp[i]);
    }
    
    fhTimePtPhotonNoCut  = new TH2F ("hTimePtPhoton_NoCut","time of photon cluster vs pT of clusters, no cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimePtPhotonNoCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTimePtPhotonNoCut->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePtPhotonNoCut);
    
    fhTimePtPhotonSPD  = new TH2F ("hTimePtPhoton_SPD","time of  photon cluster vs pT of clusters, SPD cut", nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimePtPhotonSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fhTimePtPhotonSPD->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimePtPhotonSPD);
    
    fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,20,0,20);
    fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
    fhTimeNPileUpVertSPD->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertSPD);
    
    fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 20,0,20 );
    fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
    fhTimeNPileUpVertTrack->SetXTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNPileUpVertTrack);
    
    fhPtPhotonNPileUpSPDVtx  = new TH2F ("hPtPhoton_NPileUpVertSPD","pT of cluster vs N pile-up SPD vertex",
                                         nptbins,ptmin,ptmax,20,0,20);
    fhPtPhotonNPileUpSPDVtx->SetYTitle("# vertex ");
    fhPtPhotonNPileUpSPDVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpSPDVtx);
	  
    fhPtPhotonNPileUpTrkVtx  = new TH2F ("hPtPhoton_NPileUpVertTracks","pT of cluster vs N pile-up Tracks vertex",
                                         nptbins,ptmin,ptmax, 20,0,20 );
    fhPtPhotonNPileUpTrkVtx->SetYTitle("# vertex ");
    fhPtPhotonNPileUpTrkVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpTrkVtx);
	  
    fhPtPhotonNPileUpSPDVtxTimeCut  = new TH2F ("hPtPhoton_NPileUpVertSPD_TimeCut","pT of cluster vs N pile-up SPD vertex, |tof| < 25 ns",
                                                nptbins,ptmin,ptmax,20,0,20);
    fhPtPhotonNPileUpSPDVtxTimeCut->SetYTitle("# vertex ");
    fhPtPhotonNPileUpSPDVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpSPDVtxTimeCut);
	  
    fhPtPhotonNPileUpTrkVtxTimeCut  = new TH2F ("hPtPhoton_NPileUpVertTracks_TimeCut","pT of cluster vs N pile-up Tracks vertex, |tof| < 25 ns",
                                                nptbins,ptmin,ptmax, 20,0,20 );
    fhPtPhotonNPileUpTrkVtxTimeCut->SetYTitle("# vertex ");
    fhPtPhotonNPileUpTrkVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpTrkVtxTimeCut);
	  
    fhPtPhotonNPileUpSPDVtxTimeCut2  = new TH2F ("hPtPhoton_NPileUpVertSPD_TimeCut2","pT of cluster vs N pile-up SPD vertex, -25 < tof < 75 ns",
                                                 nptbins,ptmin,ptmax,20,0,20);
    fhPtPhotonNPileUpSPDVtxTimeCut2->SetYTitle("# vertex ");
    fhPtPhotonNPileUpSPDVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpSPDVtxTimeCut2);
	  
    fhPtPhotonNPileUpTrkVtxTimeCut2  = new TH2F ("hPtPhoton_NPileUpVertTracks_TimeCut2","pT of cluster vs N pile-up Tracks vertex,  -25 < tof < 75 ns",
                                                 nptbins,ptmin,ptmax, 20,0,20 );
    fhPtPhotonNPileUpTrkVtxTimeCut2->SetYTitle("# vertex ");
    fhPtPhotonNPileUpTrkVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPhotonNPileUpTrkVtxTimeCut2);
    
  }

  
  
  if(IsDataMC())
  {
    TString ptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}", "#pi^{0}","#eta",
      "e^{#pm}","#gamma->e^{#pm}","hadron?","Anti-N","Anti-P",
      "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}","String"                                } ;
    
    TString pname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Pi0","Eta","Electron",
      "Conversion", "Hadron", "AntiNeutron","AntiProton",
      "PhotonPrompt","PhotonFragmentation","PhotonISR","String" } ;
    
    for(Int_t i = 0; i < fNOriginHistograms; i++)
    {
      fhMCE[i]  = new TH1F(Form("hE_MC%s",pname[i].Data()),
                           Form("cluster from %s : E ",ptype[i].Data()),
                           nptbins,ptmin,ptmax);
      fhMCE[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMCE[i]) ;
      
      fhMCPt[i]  = new TH1F(Form("hPt_MC%s",pname[i].Data()),
                            Form("cluster from %s : #it{p}_{T} ",ptype[i].Data()),
                            nptbins,ptmin,ptmax);
      fhMCPt[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCPt[i]) ;
      
      fhMCEta[i]  = new TH2F(Form("hEta_MC%s",pname[i].Data()),
                             Form("cluster from %s : #eta ",ptype[i].Data()),
                             nptbins,ptmin,ptmax,netabins,etamin,etamax);
      fhMCEta[i]->SetYTitle("#eta");
      fhMCEta[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMCEta[i]) ;
      
      fhMCPhi[i]  = new TH2F(Form("hPhi_MC%s",pname[i].Data()),
                             Form("cluster from %s : #phi ",ptype[i].Data()),
                             nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhMCPhi[i]->SetYTitle("#phi (rad)");
      fhMCPhi[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMCPhi[i]) ;
      
      
      fhMCDeltaE[i]  = new TH2F (Form("hDeltaE_MC%s",pname[i].Data()),
                                 Form("MC - Reco E from %s",pname[i].Data()),
                                 nptbins,ptmin,ptmax, 200,-50,50);
      fhMCDeltaE[i]->SetYTitle("#Delta #it{E} (GeV)");
      fhMCDeltaE[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhMCDeltaE[i]);
      
      fhMCDeltaPt[i]  = new TH2F (Form("hDeltaPt_MC%s",pname[i].Data()),
                                  Form("MC - Reco #it{p}_{T} from %s",pname[i].Data()),
                                  nptbins,ptmin,ptmax, 200,-50,50);
      fhMCDeltaPt[i]->SetXTitle("p_{T,rec} (GeV/#it{c})");
      fhMCDeltaPt[i]->SetYTitle("#Delta #it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhMCDeltaPt[i]);
      
      fhMC2E[i]  = new TH2F (Form("h2E_MC%s",pname[i].Data()),
                             Form("E distribution, reconstructed vs generated from %s",pname[i].Data()),
                             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMC2E[i]->SetXTitle("#it{E}_{rec} (GeV)");
      fhMC2E[i]->SetYTitle("#it{E}_{gen} (GeV)");
      outputContainer->Add(fhMC2E[i]);
      
      fhMC2Pt[i]  = new TH2F (Form("h2Pt_MC%s",pname[i].Data()),
                              Form("p_T distribution, reconstructed vs generated from %s",pname[i].Data()),
                              nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
      fhMC2Pt[i]->SetXTitle("p_{T,rec} (GeV/#it{c})");
      fhMC2Pt[i]->SetYTitle("p_{T,gen} (GeV/#it{c})");
      outputContainer->Add(fhMC2Pt[i]);
      
      
    }
    
    TString pptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}",
      "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}"} ;
    
    TString ppname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay",
      "PhotonPrompt","PhotonFragmentation","PhotonISR"} ;
    
    for(Int_t i = 0; i < fNPrimaryHistograms; i++)
    {
      fhEPrimMC[i]  = new TH1F(Form("hEPrim_MC%s",ppname[i].Data()),
                               Form("primary photon %s : E ",pptype[i].Data()),
                               nptbins,ptmin,ptmax);
      fhEPrimMC[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEPrimMC[i]) ;
      
      fhPtPrimMC[i]  = new TH1F(Form("hPtPrim_MC%s",ppname[i].Data()),
                                Form("primary photon %s : #it{p}_{T} ",pptype[i].Data()),
                                nptbins,ptmin,ptmax);
      fhPtPrimMC[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMC[i]) ;
      
      fhYPrimMC[i]  = new TH2F(Form("hYPrim_MC%s",ppname[i].Data()),
                               Form("primary photon %s : Rapidity ",pptype[i].Data()),
                               nptbins,ptmin,ptmax,200,-2,2);
      fhYPrimMC[i]->SetYTitle("Rapidity");
      fhYPrimMC[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhYPrimMC[i]) ;
      
      fhEtaPrimMC[i]  = new TH2F(Form("hEtaPrim_MC%s",ppname[i].Data()),
                               Form("primary photon %s : #eta",pptype[i].Data()),
                               nptbins,ptmin,ptmax,200,-2,2);
      fhEtaPrimMC[i]->SetYTitle("#eta");
      fhEtaPrimMC[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEtaPrimMC[i]) ;
      
      fhPhiPrimMC[i]  = new TH2F(Form("hPhiPrim_MC%s",ppname[i].Data()),
                                 Form("primary photon %s : #phi ",pptype[i].Data()),
                                 nptbins,ptmin,ptmax,nphibins,0,TMath::TwoPi());
      fhPhiPrimMC[i]->SetYTitle("#phi (rad)");
      fhPhiPrimMC[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhPhiPrimMC[i]) ;
      
      
      fhEPrimMCAcc[i]  = new TH1F(Form("hEPrimAcc_MC%s",ppname[i].Data()),
                                  Form("primary photon %s in acceptance: E ",pptype[i].Data()),
                                  nptbins,ptmin,ptmax);
      fhEPrimMCAcc[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEPrimMCAcc[i]) ;
      
      fhPtPrimMCAcc[i]  = new TH1F(Form("hPtPrimAcc_MC%s",ppname[i].Data()),
                                   Form("primary photon %s in acceptance: #it{p}_{T} ",pptype[i].Data()),
                                   nptbins,ptmin,ptmax);
      fhPtPrimMCAcc[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      outputContainer->Add(fhPtPrimMCAcc[i]) ;
      
      fhYPrimMCAcc[i]  = new TH2F(Form("hYPrimAcc_MC%s",ppname[i].Data()),
                                  Form("primary photon %s in acceptance: Rapidity ",pptype[i].Data()),
                                  nptbins,ptmin,ptmax,100,-1,1);
      fhYPrimMCAcc[i]->SetYTitle("Rapidity");
      fhYPrimMCAcc[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhYPrimMCAcc[i]) ;

      fhEtaPrimMCAcc[i]  = new TH2F(Form("hEtaPrimAcc_MC%s",ppname[i].Data()),
                                  Form("primary photon %s in acceptance: #eta ",pptype[i].Data()),
                                  nptbins,ptmin,ptmax,netabins,etamin,etamax);
      fhEtaPrimMCAcc[i]->SetYTitle("#eta");
      fhEtaPrimMCAcc[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhEtaPrimMCAcc[i]) ;
      
      fhPhiPrimMCAcc[i]  = new TH2F(Form("hPhiPrimAcc_MC%s",ppname[i].Data()),
                                    Form("primary photon %s in acceptance: #phi ",pptype[i].Data()),
                                    nptbins,ptmin,ptmax,nphibins,phimin,phimax);
      fhPhiPrimMCAcc[i]->SetYTitle("#phi (rad)");
      fhPhiPrimMCAcc[i]->SetXTitle("#it{E} (GeV)");
      outputContainer->Add(fhPhiPrimMCAcc[i]) ;
      
    }
    
    if(fFillSSHistograms)
    {
      TString ptypess[] = { "#gamma","hadron?","#pi^{0}","#eta","#gamma->e^{#pm}","e^{#pm}"} ;
      
      TString pnamess[] = { "Photon","Hadron","Pi0","Eta","Conversion","Electron"} ;
      
      for(Int_t i = 0; i < 6; i++)
      {
        fhMCELambda0[i]  = new TH2F(Form("hELambda0_MC%s",pnamess[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{0}^{2}",ptypess[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCELambda0[i]->SetYTitle("#lambda_{0}^{2}");
        fhMCELambda0[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCELambda0[i]) ;
        
        fhMCELambda1[i]  = new TH2F(Form("hELambda1_MC%s",pnamess[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{1}^{2}",ptypess[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCELambda1[i]->SetYTitle("#lambda_{1}^{2}");
        fhMCELambda1[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCELambda1[i]) ;
        
        fhMCEDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pnamess[i].Data()),
                                       Form("cluster from %s : E vs dispersion^{2}",ptypess[i].Data()),
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCEDispersion[i]->SetYTitle("D^{2}");
        fhMCEDispersion[i]->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCEDispersion[i]) ;
        
        fhMCNCellsE[i]  = new TH2F (Form("hNCellsE_MC%s",pnamess[i].Data()),
                                    Form("# of cells in cluster from %s vs E of clusters",ptypess[i].Data()),
                                    nptbins,ptmin,ptmax, nbins,nmin,nmax);
        fhMCNCellsE[i]->SetXTitle("#it{E} (GeV)");
        fhMCNCellsE[i]->SetYTitle("# of cells in cluster");
        outputContainer->Add(fhMCNCellsE[i]);
        
        fhMCMaxCellDiffClusterE[i]  = new TH2F (Form("hMaxCellDiffClusterE_MC%s",pnamess[i].Data()),
                                                Form("energy vs difference of cluster energy from %s - max cell energy / cluster energy, good clusters",ptypess[i].Data()),
                                                nptbins,ptmin,ptmax, 500,0,1.);
        fhMCMaxCellDiffClusterE[i]->SetXTitle("#it{E}_{cluster} (GeV) ");
        fhMCMaxCellDiffClusterE[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
        outputContainer->Add(fhMCMaxCellDiffClusterE[i]);
        
        if(!fFillOnlySimpleSSHisto)
        {
          fhMCLambda0vsClusterMaxCellDiffE0[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE0_MC%s",pnamess[i].Data()),
                                                           Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, E < 2 GeV",ptypess[i].Data()),
                                                           ssbins,ssmin,ssmax,500,0,1.);
          fhMCLambda0vsClusterMaxCellDiffE0[i]->SetXTitle("#lambda_{0}^{2}");
          fhMCLambda0vsClusterMaxCellDiffE0[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
          outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE0[i]) ;
          
          fhMCLambda0vsClusterMaxCellDiffE2[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE2_MC%s",pnamess[i].Data()),
                                                           Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, 2< E < 6 GeV",ptypess[i].Data()),
                                                           ssbins,ssmin,ssmax,500,0,1.);
          fhMCLambda0vsClusterMaxCellDiffE2[i]->SetXTitle("#lambda_{0}^{2}");
          fhMCLambda0vsClusterMaxCellDiffE2[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
          outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE2[i]) ;
          
          fhMCLambda0vsClusterMaxCellDiffE6[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE6_MC%s",pnamess[i].Data()),
                                                           Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, #it{E} > 6 GeV",ptypess[i].Data()),
                                                           ssbins,ssmin,ssmax,500,0,1.);
          fhMCLambda0vsClusterMaxCellDiffE6[i]->SetXTitle("#lambda_{0}^{2}");
          fhMCLambda0vsClusterMaxCellDiffE6[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
          outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE6[i]) ;
          
          fhMCNCellsvsClusterMaxCellDiffE0[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE0_MC%s",pnamess[i].Data()),
                                                          Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, E < 2 GeV",ptypess[i].Data()),
                                                          nbins/5,nmin,nmax/5,500,0,1.);
          fhMCNCellsvsClusterMaxCellDiffE0[i]->SetXTitle("N cells in cluster");
          fhMCNCellsvsClusterMaxCellDiffE0[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
          outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE0[i]) ;
          
          fhMCNCellsvsClusterMaxCellDiffE2[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE2_MC%s",pnamess[i].Data()),
                                                          Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, 2< E < 6 GeV",ptypess[i].Data()),
                                                          nbins/5,nmin,nmax/5,500,0,1.);
          fhMCNCellsvsClusterMaxCellDiffE2[i]->SetXTitle("N cells in cluster");
          fhMCNCellsvsClusterMaxCellDiffE2[i]->SetYTitle("(#it{E}_{cluster} - #it{E}_{cell max})/ #it{E}_{cluster}");
          outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE2[i]) ;
          
          fhMCNCellsvsClusterMaxCellDiffE6[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE6_MC%s",pnamess[i].Data()),
                                                          Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, #it{E} > 6 GeV",ptypess[i].Data()),
                                                          nbins/5,nmin,nmax/5,500,0,1.);
          fhMCNCellsvsClusterMaxCellDiffE6[i]->SetXTitle("N cells in cluster");
          fhMCNCellsvsClusterMaxCellDiffE6[i]->SetYTitle("#it{E} (GeV)");
          outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE6[i]) ;
          
          if(fCalorimeter=="EMCAL")
          {
            fhMCEDispEta[i]  = new TH2F (Form("hEDispEtaE_MC%s",pnamess[i].Data()),
                                         Form("cluster from %s : #sigma^{2}_{#eta #eta} = #Sigma w_{i}(#eta_{i} - <#eta>)^{2}/ #Sigma w_{i} vs E",ptypess[i].Data()),
                                         nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
            fhMCEDispEta[i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispEta[i]->SetYTitle("#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEta[i]);
            
            fhMCEDispPhi[i]  = new TH2F (Form("hEDispPhiE_MC%s",pnamess[i].Data()),
                                         Form("cluster from %s : #sigma^{2}_{#phi #phi} = #Sigma w_{i}(#phi_{i} - <#phi>)^{2} / #Sigma w_{i} vs E",ptypess[i].Data()),
                                         nptbins,ptmin,ptmax, ssbins,ssmin,ssmax);
            fhMCEDispPhi[i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispPhi[i]->SetYTitle("#sigma^{2}_{#phi #phi}");
            outputContainer->Add(fhMCEDispPhi[i]);
            
            fhMCESumEtaPhi[i]  = new TH2F (Form("hESumEtaPhiE_MC%s",pnamess[i].Data()),
                                           Form("cluster from %s : #delta^{2}_{#eta #phi} = #Sigma w_{i}(#phi_{i} #eta_{i} ) / #Sigma w_{i} - <#phi><#eta> vs E",ptypess[i].Data()),
                                           nptbins,ptmin,ptmax, 2*ssbins,-ssmax,ssmax);
            fhMCESumEtaPhi[i]->SetXTitle("#it{E} (GeV)");
            fhMCESumEtaPhi[i]->SetYTitle("#delta^{2}_{#eta #phi}");
            outputContainer->Add(fhMCESumEtaPhi[i]);
            
            fhMCEDispEtaPhiDiff[i]  = new TH2F (Form("hEDispEtaPhiDiffE_MC%s",pnamess[i].Data()),
                                                Form("cluster from %s : #sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta} vs E",ptypess[i].Data()),
                                                nptbins,ptmin,ptmax,200,-10,10);
            fhMCEDispEtaPhiDiff[i]->SetXTitle("#it{E} (GeV)");
            fhMCEDispEtaPhiDiff[i]->SetYTitle("#sigma^{2}_{#phi #phi}-#sigma^{2}_{#eta #eta}");
            outputContainer->Add(fhMCEDispEtaPhiDiff[i]);
            
            fhMCESphericity[i]  = new TH2F (Form("hESphericity_MC%s",pnamess[i].Data()),
                                            Form("cluster from %s : (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi}) vs E",ptypess[i].Data()),
                                            nptbins,ptmin,ptmax, 200,-1,1);
            fhMCESphericity[i]->SetXTitle("#it{E} (GeV)");
            fhMCESphericity[i]->SetYTitle("s = (#sigma^{2}_{#phi #phi} - #sigma^{2}_{#eta #eta}) / (#sigma^{2}_{#eta #eta} + #sigma^{2}_{#phi #phi})");
            outputContainer->Add(fhMCESphericity[i]);
            
            for(Int_t ie = 0; ie < 7; ie++)
            {
              fhMCDispEtaDispPhi[ie][i] = new TH2F (Form("hMCDispEtaDispPhi_EBin%d_MC%s",ie,pnamess[i].Data()),
                                                    Form("cluster from %s : #sigma^{2}_{#phi #phi} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pnamess[i].Data(),bin[ie],bin[ie+1]),
                                                    ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
              fhMCDispEtaDispPhi[ie][i]->SetXTitle("#sigma^{2}_{#eta #eta}");
              fhMCDispEtaDispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
              outputContainer->Add(fhMCDispEtaDispPhi[ie][i]);
              
              fhMCLambda0DispEta[ie][i] = new TH2F (Form("hMCLambda0DispEta_EBin%d_MC%s",ie,pnamess[i].Data()),
                                                    Form("cluster from %s : #lambda^{2}_{0} vs #sigma^{2}_{#eta #eta} for %d < E < %d GeV",pnamess[i].Data(),bin[ie],bin[ie+1]),
                                                    ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
              fhMCLambda0DispEta[ie][i]->SetXTitle("#lambda^{2}_{0}");
              fhMCLambda0DispEta[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
              outputContainer->Add(fhMCLambda0DispEta[ie][i]);
              
              fhMCLambda0DispPhi[ie][i] = new TH2F (Form("hMCLambda0DispPhi_EBin%d_MC%s",ie,pnamess[i].Data()),
                                                    Form("cluster from %s :#lambda^{2}_{0} vs #sigma^{2}_{#phi #phi} for %d < E < %d GeV",pnamess[i].Data(),bin[ie],bin[ie+1]),
                                                    ssbins,ssmin,ssmax , ssbins,ssmin,ssmax);
              fhMCLambda0DispPhi[ie][i]->SetXTitle("#lambda^{2}_{0}");
              fhMCLambda0DispPhi[ie][i]->SetYTitle("#sigma^{2}_{#phi #phi}");
              outputContainer->Add(fhMCLambda0DispPhi[ie][i]);
            }
          }
        }
      }// loop
      
      if(!GetReader()->IsEmbeddedClusterSelectionOn())
      {
        fhMCPhotonELambda0NoOverlap  = new TH2F("hELambda0_MCPhoton_NoOverlap",
                                                "cluster from Photon : E vs #lambda_{0}^{2}",
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCPhotonELambda0NoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0NoOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCPhotonELambda0NoOverlap) ;
        
        fhMCPhotonELambda0TwoOverlap  = new TH2F("hELambda0_MCPhoton_TwoOverlap",
                                                 "cluster from Photon : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCPhotonELambda0TwoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0TwoOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCPhotonELambda0TwoOverlap) ;
        
        fhMCPhotonELambda0NOverlap  = new TH2F("hELambda0_MCPhoton_NOverlap",
                                               "cluster from Photon : E vs #lambda_{0}^{2}",
                                               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhMCPhotonELambda0NOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0NOverlap->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhMCPhotonELambda0NOverlap) ;
        
      } //No embedding
      
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        
        fhEmbeddedSignalFractionEnergy  = new TH2F("hEmbeddedSignal_FractionEnergy",
                                                   "Energy Fraction of embedded signal versus cluster energy",
                                                   nptbins,ptmin,ptmax,100,0.,1.);
        fhEmbeddedSignalFractionEnergy->SetYTitle("Fraction");
        fhEmbeddedSignalFractionEnergy->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbeddedSignalFractionEnergy) ;
        
        fhEmbedPhotonELambda0FullSignal  = new TH2F("hELambda0_EmbedPhoton_FullSignal",
                                                    "cluster from Photon embedded with more than 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPhotonELambda0FullSignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0FullSignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0FullSignal) ;
        
        fhEmbedPhotonELambda0MostlySignal  = new TH2F("hELambda0_EmbedPhoton_MostlySignal",
                                                      "cluster from Photon embedded with 50% to 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPhotonELambda0MostlySignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0MostlySignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0MostlySignal) ;
        
        fhEmbedPhotonELambda0MostlyBkg  = new TH2F("hELambda0_EmbedPhoton_MostlyBkg",
                                                   "cluster from Photon embedded with 10% to 50% energy in cluster : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPhotonELambda0MostlyBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0MostlyBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0MostlyBkg) ;
        
        fhEmbedPhotonELambda0FullBkg  = new TH2F("hELambda0_EmbedPhoton_FullBkg",
                                                 "cluster from Photonm embedded with 0% to 10% energy in cluster : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPhotonELambda0FullBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0FullBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0FullBkg) ;
        
        fhEmbedPi0ELambda0FullSignal  = new TH2F("hELambda0_EmbedPi0_FullSignal",
                                                 "cluster from Pi0 embedded with more than 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPi0ELambda0FullSignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0FullSignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0FullSignal) ;
        
        fhEmbedPi0ELambda0MostlySignal  = new TH2F("hELambda0_EmbedPi0_MostlySignal",
                                                   "cluster from Pi0 embedded with 50% to 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPi0ELambda0MostlySignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0MostlySignal->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0MostlySignal) ;
        
        fhEmbedPi0ELambda0MostlyBkg  = new TH2F("hELambda0_EmbedPi0_MostlyBkg",
                                                "cluster from Pi0 embedded with 10% to 50% energy in cluster : E vs #lambda_{0}^{2}",
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPi0ELambda0MostlyBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0MostlyBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0MostlyBkg) ;
        
        fhEmbedPi0ELambda0FullBkg  = new TH2F("hELambda0_EmbedPi0_FullBkg",
                                              "cluster from Pi0 embedded with 0% to 10% energy in cluster : E vs #lambda_{0}^{2}",
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhEmbedPi0ELambda0FullBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0FullBkg->SetXTitle("#it{E} (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0FullBkg) ;
        
      }// embedded histograms
      
      
    }// Fill SS MC histograms
    
  }//Histos with MC
  
  return outputContainer ;
  
}

//_______________________
void AliAnaPhoton::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD())
  {
    printf("AliAnaPhoton::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD())
  {
    printf("AliAnaPhoton::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
  if(GetReader()->GetDataType() == AliCaloTrackReader::kMC) GetCaloPID()->SwitchOnBayesian();
  
}

//____________________________________________________________________________
void AliAnaPhoton::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaPhoton_");
  
  fCalorimeter = "EMCAL" ;
  fMinDist     = 2.;
  fMinDist2    = 4.;
  fMinDist3    = 5.;
	
  fTimeCutMin  =-1000000;
  fTimeCutMax  = 1000000;
  fNCellsCut   = 0;
	
  fRejectTrackMatch       = kTRUE ;
	
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillAOD()
{
  //Do photon analysis and fill aods
  
  //Get the vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  
  //Select the Calorimeter of the photon
  TObjArray * pl = 0x0;
  AliVCaloCells* cells    = 0;
  if      (fCalorimeter == "PHOS" )
  {
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL")
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl)
  {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
  
  TLorentzVector mom;

  // Loop on raw clusters before filtering in the reader and fill control histogram
  if((GetReader()->GetEMCALClusterListName()=="" && fCalorimeter=="EMCAL") || fCalorimeter=="PHOS")
  {
    for(Int_t iclus = 0; iclus < GetReader()->GetInputEvent()->GetNumberOfCaloClusters(); iclus++ )
    {
      AliVCluster * clus = GetReader()->GetInputEvent()->GetCaloCluster(iclus);
      if     (fCalorimeter == "PHOS"  && clus->IsPHOS()  && clus->E() > GetReader()->GetPHOSPtMin() )
      {
        fhClusterCutsE [0]->Fill(clus->E());
        
        clus->GetMomentum(mom,GetVertex(0)) ;
        fhClusterCutsPt[0]->Fill(mom.Pt());
      }
      else if(fCalorimeter == "EMCAL" && clus->IsEMCAL() && clus->E() > GetReader()->GetEMCALPtMin())
      {
        fhClusterCutsE [0]->Fill(clus->E());
        
        clus->GetMomentum(mom,GetVertex(0)) ;
        fhClusterCutsPt[0]->Fill(mom.Pt());
      }
    }
  }
  else
  { // reclusterized
    TClonesArray * clusterList = 0;
    
    if(GetReader()->GetInputEvent()->FindListObject(GetReader()->GetEMCALClusterListName()))
      clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetInputEvent()->FindListObject(GetReader()->GetEMCALClusterListName()));
    else if(GetReader()->GetOutputEvent())
      clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetOutputEvent()->FindListObject(GetReader()->GetEMCALClusterListName()));
    
    if(clusterList)
    {
      Int_t nclusters = clusterList->GetEntriesFast();
      for (Int_t iclus =  0; iclus <  nclusters; iclus++)
      {
        AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
        if(clus)
        {
          fhClusterCutsE [0]->Fill(clus->E());
          
          clus->GetMomentum(mom,GetVertex(0)) ;
          fhClusterCutsPt[0]->Fill(mom.Pt());
        }
      }
    }
  }
  
  //Init arrays, variables, get number of clusters
  TLorentzVector mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  
  //----------------------------------------------------
  // Fill AOD with PHOS/EMCAL AliAODPWG4Particle objects
  //----------------------------------------------------
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));
    //printf("calo %d, %f\n",icalo,calo->E());
    
    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ;
    if (GetMixedEvent())
    {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ;
      //Get the vertex and check it is not too large in z
      if(TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut()) continue;
    }
    
    //Cluster selection, not charged, with photon id and in fiducial cut
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC)
    {
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;
    }//Assume that come from vertex in straight line
    else
    {
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
    
    //-----------------------------
    // Cluster selection
    //-----------------------------
    Int_t nMaxima = GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells); // NLM
    if(!ClusterSelected(calo,mom,nMaxima)) continue;
    
    //----------------------------
    //Create AOD for analysis
    //----------------------------
    AliAODPWG4Particle aodph = AliAODPWG4Particle(mom);
    
    //...............................................
    //Set the indeces of the original caloclusters (MC, ID), and calorimeter
    Int_t label = calo->GetLabel();
    aodph.SetLabel(label);
    aodph.SetCaloLabel(calo->GetID(),-1);
    aodph.SetDetector(fCalorimeter);
    //printf("Index %d, Id %d, iaod %d\n",icalo, calo->GetID(),GetOutputAODBranch()->GetEntriesFast());
    
    //...............................................
    //Set bad channel distance bit
    Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
    if     (distBad > fMinDist3) aodph.SetDistToBad(2) ;
    else if(distBad > fMinDist2) aodph.SetDistToBad(1) ;
    else                         aodph.SetDistToBad(0) ;
    //printf("DistBad %f Bit %d\n",distBad, aodph.DistToBad());
    
    //-------------------------------------
    // Play with the MC stack if available
    //-------------------------------------
    
    //Check origin of the candidates
    Int_t tag = -1;
    
    if(IsDataMC())
    {
      tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(),fCalorimeter);
      aodph.SetTag(tag);
      
      if(GetDebug() > 0)
        printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also
    
    //--------------------------------------------------------
    //Fill some shower shape histograms before PID is applied
    //--------------------------------------------------------
    
    FillShowerShapeHistograms(calo,tag);
    
    //-------------------------------------
    //PID selection or bit setting
    //-------------------------------------
    
    //...............................................
    // Data, PID check on
    if(IsCaloPIDOn())
    {
      // Get most probable PID, 2 options check bayesian PID weights or redo PID
      // By default, redo PID
      
      aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(calo));
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;
      
    }
    
    //...............................................
    // Data, PID check off
    else
    {
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetIdentifiedParticleType already called in SetPIDBits.
      
      GetCaloPID()->SetPIDBits(calo,&aodph, GetCaloUtils(),GetReader()->GetInputEvent());
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PID Bits set \n");
    }
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",
                              aodph.Pt(), aodph.GetIdentifiedParticleType());
    
    fhClusterCutsE [9]->Fill(calo->E());
    fhClusterCutsPt[9]->Fill(mom.Pt());
    
    Int_t   nSM  = GetModuleNumber(calo);
    if(nSM < GetCaloUtils()->GetNumberOfSuperModulesUsed() && nSM >=0)
    {
      fhEPhotonSM ->Fill(mom.E (),nSM);
      fhPtPhotonSM->Fill(mom.Pt(),nSM);
    }
    
    fhNLocMax->Fill(calo->E(),nMaxima);
    
    // Matching after cuts
    if( fFillTMHisto )          FillTrackMatchingResidualHistograms(calo,1);
    
    // Fill histograms to undertand pile-up before other cuts applied
    // Remember to relax time cuts in the reader
    if( fFillPileUpHistograms ) FillPileUpHistograms(calo,cells);
    
    // Add number of local maxima to AOD, method name in AOD to be FIXED
    aodph.SetFiducialArea(nMaxima);
    
    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());
  
}

//______________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms()
{
  //Fill histograms
  
  //In case of simulated data, fill acceptance histograms
  if(IsDataMC()) FillAcceptanceHistograms();
  
  // Get vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  //fhVertex->Fill(v[0],v[1],v[2]);
  if(TMath::Abs(v[2]) > GetZvertexCut()) return ; // done elsewhere for Single Event analysis, but there for mixed event
  
  //----------------------------------
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  Float_t cen = GetEventCentrality();
  // printf("++++++++++ GetEventCentrality() %f\n",cen);
  
  Float_t ep  = GetEventPlaneAngle();
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ph->GetIdentifiedParticleType();
    
    if(GetDebug() > 3)
      printf("AliAnaPhoton::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n",
             ph->GetIdentifiedParticleType(),ph->GetTag(), (ph->GetDetector()).Data()) ;
    
    //If PID used, fill histos with photons in Calorimeter fCalorimeter
    if(IsCaloPIDOn() && pdg != AliCaloPID::kPhoton) continue;
    
    if(ph->GetDetector() != fCalorimeter) continue;
    
    if(GetDebug() > 2)
      printf("AliAnaPhoton::MakeAnalysisFillHistograms() - ID Photon: pt %f, phi %f, eta %f\n", ph->Pt(),ph->Phi(),ph->Eta()) ;
    
    //................................
    //Fill photon histograms
    Float_t ptcluster  = ph->Pt();
    Float_t phicluster = ph->Phi();
    Float_t etacluster = ph->Eta();
    Float_t ecluster   = ph->E();
    
    fhEPhoton   ->Fill(ecluster);
    fhPtPhoton  ->Fill(ptcluster);
    fhPhiPhoton ->Fill(ptcluster,phicluster);
    fhEtaPhoton ->Fill(ptcluster,etacluster);
    if     (ecluster   > 0.5) fhEtaPhiPhoton  ->Fill(etacluster, phicluster);
    else if(GetMinPt() < 0.5) fhEtaPhi05Photon->Fill(etacluster, phicluster);
    
    fhPtCentralityPhoton ->Fill(ptcluster,cen) ;
    fhPtEventPlanePhoton ->Fill(ptcluster,ep ) ;
    
    //Get original cluster, to recover some information
    AliVCaloCells* cells    = 0;
    TObjArray * clusters    = 0;
    if(fCalorimeter == "EMCAL")
    {
      cells    = GetEMCALCells();
      clusters = GetEMCALClusters();
    }
    else
    {
      cells    = GetPHOSCells();
      clusters = GetPHOSClusters();
    }
    
    Int_t iclus = -1;
    AliVCluster *cluster = FindCluster(clusters,ph->GetCaloLabel(0),iclus);
    if(cluster)
    {
      Float_t maxCellFraction = 0;
      Int_t absID = GetCaloUtils()->GetMaxEnergyCell(cells, cluster, maxCellFraction);
      if( absID < 0 ) AliFatal("Wrong absID");

      // Control histograms
      fhMaxCellDiffClusterE->Fill(ph->E() ,maxCellFraction);
      fhNCellsE            ->Fill(ph->E() ,cluster->GetNCells());
      fhTimePt             ->Fill(ph->Pt(),cluster->GetTOF()*1.e9);

      if(cells)
      {
        for(Int_t icell = 0; icell <  cluster->GetNCells(); icell++)
          fhCellsE->Fill(ph->E(),cells->GetCellAmplitude(cluster->GetCellsAbsId()[icell]));
      }
    }
    
    //.......................................
    //Play with the MC data if available
    if(IsDataMC())
    {
      if(GetDebug()>0)
      {
        if(GetReader()->ReadStack() && !GetMCStack())
        {
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called?\n");
        }
        else if(GetReader()->ReadAODMCParticles() && !GetReader()->GetAODMCParticles())
        {
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
        }
      }
            
      //....................................................................
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      
      if(label < 0)
      {
        if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
        continue;
      }
      
      Float_t eprim   = 0;
      Float_t ptprim  = 0;
      Bool_t ok = kFALSE;
      TLorentzVector primary = GetMCAnalysisUtils()->GetMother(label,GetReader(),ok);
      if(ok)
      {
        eprim   = primary.Energy();
        ptprim  = primary.Pt();
      }
      
      Int_t tag =ph->GetTag();
      Int_t mcParticleTag = -1;
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[kmcPhoton])
      {
        fhMCE  [kmcPhoton] ->Fill(ecluster);
        fhMCPt [kmcPhoton] ->Fill(ptcluster);
        fhMCPhi[kmcPhoton] ->Fill(ecluster,phicluster);
        fhMCEta[kmcPhoton] ->Fill(ecluster,etacluster);
        
        fhMC2E     [kmcPhoton] ->Fill(ecluster, eprim);
        fhMC2Pt    [kmcPhoton] ->Fill(ptcluster, ptprim);
        fhMCDeltaE [kmcPhoton] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcPhoton] ->Fill(ptcluster,ptprim-ptcluster);
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) &&
           fhMCE[kmcConversion])
        {
          fhMCE  [kmcConversion] ->Fill(ecluster);
          fhMCPt [kmcConversion] ->Fill(ptcluster);
          fhMCPhi[kmcConversion] ->Fill(ecluster,phicluster);
          fhMCEta[kmcConversion] ->Fill(ecluster,etacluster);
          
          fhMC2E     [kmcConversion] ->Fill(ecluster, eprim);
          fhMC2Pt    [kmcConversion] ->Fill(ptcluster, ptprim);
          fhMCDeltaE [kmcConversion] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcConversion] ->Fill(ptcluster,ptprim-ptcluster);
        }
        
        if     (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt))
        {
          mcParticleTag = kmcPrompt;
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation))
        {
          mcParticleTag = kmcFragmentation;
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR))
        {
          mcParticleTag = kmcISR;
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) &&
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))
        {
          mcParticleTag = kmcPi0Decay;
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0))
        {
          mcParticleTag = kmcPi0;
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) &&
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay))
        {
          mcParticleTag = kmcEta;
        }
        else
        {
          mcParticleTag = kmcOtherDecay;
        }
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron))
      {
        mcParticleTag = kmcAntiNeutron;
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton))
      {
        mcParticleTag = kmcAntiProton;
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))
      {
        mcParticleTag = kmcElectron;
      }
      else if( fhMCE[kmcOther])
      {
        mcParticleTag = kmcOther;
        
        //		 printf(" AliAnaPhoton::MakeAnalysisFillHistograms() - Label %d, pT %2.3f Unknown, bits set: ",
        //					ph->GetLabel(),ph->Pt());
        //		  for(Int_t i = 0; i < 20; i++) {
        //			  if(GetMCAnalysisUtils()->CheckTagBit(tag,i)) printf(" %d, ",i);
        //		  }
        //		  printf("\n");
        
      }
      
      if(mcParticleTag >= 0 && fhMCE[mcParticleTag])
      {
        fhMCE  [mcParticleTag]->Fill(ecluster);
        fhMCPt [mcParticleTag]->Fill(ptcluster);
        fhMCPhi[mcParticleTag]->Fill(ecluster,phicluster);
        fhMCEta[mcParticleTag]->Fill(ecluster,etacluster);
        
        fhMC2E     [mcParticleTag]->Fill(ecluster, eprim);
        fhMC2Pt    [mcParticleTag]->Fill(ptcluster, ptprim);
        fhMCDeltaE [mcParticleTag]->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcParticleTag]->Fill(ptcluster,ptprim-ptcluster);
      }
    }//Histograms with MC
    
  }// aod loop
  
}


//__________________________________________________________________
void AliAnaPhoton::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
  printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
  printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  printf("Reject clusters with a track matched = %d\n",fRejectTrackMatch);
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
}
