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
#include <TH3D.h>
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
    fNCellsCut(0),                fFillSSHistograms(kFALSE),        
    fNOriginHistograms(8),        fNPrimaryHistograms(4),

    // Histograms
    fhNCellsE(0),                 fhCellsE(0),   // Control histograms            
    fhMaxCellDiffClusterE(0),     fhTimeE(0),    // Control histograms
    fhEPhoton(0),                 fhPtPhoton(0),  
    fhPhiPhoton(0),               fhEtaPhoton(0), 
    fhEtaPhiPhoton(0),            fhEtaPhi05Photon(0),

    // Shower shape histograms
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

    // MC histograms
    fhMCPhotonELambda0NoOverlap(0),       fhMCPhotonELambda0TwoOverlap(0),      fhMCPhotonELambda0NOverlap(0),
    // Embedding
    fhEmbeddedSignalFractionEnergy(0),     
    fhEmbedPhotonELambda0FullSignal(0),   fhEmbedPhotonELambda0MostlySignal(0),  
    fhEmbedPhotonELambda0MostlyBkg(0),    fhEmbedPhotonELambda0FullBkg(0),       
    fhEmbedPi0ELambda0FullSignal(0),      fhEmbedPi0ELambda0MostlySignal(0),    
    fhEmbedPi0ELambda0MostlyBkg(0),       fhEmbedPi0ELambda0FullBkg(0)

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
    fhYPrimMC  [i] = 0;
    
    fhPtPrimMCAcc [i] = 0;
    fhEPrimMCAcc  [i] = 0;
    fhPhiPrimMCAcc[i] = 0;
    fhYPrimMCAcc  [i] = 0;
  }  
  
  for(Int_t i = 0; i < 6; i++)
  {
    fhMCELambda0    [i]                  = 0;
    fhMCELambda1    [i]                  = 0;
    fhMCEDispersion [i]                  = 0;
    fhMCNCellsE     [i]                  = 0; 
    fhMCMaxCellDiffClusterE[i]           = 0; 
    fhMCLambda0vsClusterMaxCellDiffE0[i] = 0;
    fhMCLambda0vsClusterMaxCellDiffE2[i] = 0;
    fhMCLambda0vsClusterMaxCellDiffE6[i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE0 [i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE2 [i] = 0;
    fhMCNCellsvsClusterMaxCellDiffE6 [i] = 0;
  }
  
  for(Int_t i = 0; i < 5; i++) fhClusterCuts[i] = 0;
  
   // Track matching residuals
   for(Int_t i = 0; i < 2; i++)
   {
     fhTrackMatchedDEta[i] = 0;                fhTrackMatchedDPhi[i] = 0;         fhTrackMatchedDEtaDPhi[i] = 0;   
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

//__________________________________________________________________________
Bool_t  AliAnaPhoton::ClusterSelected(AliVCluster* calo, TLorentzVector mom) 
{
  //Select clusters if they pass different cuts
  if(GetDebug() > 2) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           calo->E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
    
  fhClusterCuts[1]->Fill(calo->E());
  
  //.......................................
  //If too small or big energy, skip it
  if(calo->E() < GetMinEnergy() || calo->E() > GetMaxEnergy() ) return kFALSE ; 
  
  if(GetDebug() > 2) printf("\t Cluster %d Pass E Cut \n",calo->GetID());
  
  fhClusterCuts[2]->Fill(calo->E());

  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if(tof < fTimeCutMin || tof > fTimeCutMax) return kFALSE;
  
  if(GetDebug() > 2)  printf("\t Cluster %d Pass Time Cut \n",calo->GetID());
  
  fhClusterCuts[3]->Fill(calo->E());

  //.......................................
  if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) return kFALSE;
  
  if(GetDebug() > 2) printf("\t Cluster %d Pass NCell Cut \n",calo->GetID());
  
  fhClusterCuts[4]->Fill(calo->E());

  //.......................................
  //Check acceptance selection
  if(IsFiducialCutOn()){
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(! in ) return kFALSE ;
  }
  
  if(GetDebug() > 2) printf("Fiducial cut passed \n");
  
  fhClusterCuts[5]->Fill(calo->E());

  //.......................................
  //Skip matched clusters with tracks
  
  // Fill matching residual histograms before PID cuts
  if(fFillTMHisto) FillTrackMatchingResidualHistograms(calo,0);
  
  if(fRejectTrackMatch){
    if(IsTrackMatched(calo,GetReader()->GetInputEvent())) {
      if(GetDebug() > 2) printf("\t Reject track-matched clusters\n");
      return kFALSE ;
    }
    else  
      if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  }// reject matched clusters
  
  fhClusterCuts[6]->Fill(calo->E());

  //.......................................
  //Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  if(distBad < fMinDist) {//In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  }
  else if(GetDebug() > 2) printf("\t Bad channel cut passed %4.2f > %2.2f \n",distBad, fMinDist);
  
  fhClusterCuts[7]->Fill(calo->E());
  
  if(GetDebug() > 0) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(), 
           calo->E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //All checks passed, cluster selected
  return kTRUE;
    
}

//_____________________________________________________________
void AliAnaPhoton::FillAcceptanceHistograms(){
  //Fill acceptance histograms if MC data is available
  
  if(GetReader()->ReadStack()){	
    AliStack * stack = GetMCStack();
    if(stack){
      for(Int_t i=0 ; i<stack->GetNtrack(); i++){
        TParticle * prim = stack->Particle(i) ;
        Int_t pdg = prim->GetPdgCode();
        //printf("i %d, %s %d  %s %d \n",i, stack->Particle(i)->GetName(), stack->Particle(i)->GetPdgCode(),
        //                             prim->GetName(), prim->GetPdgCode());
        
        if(pdg == 22){
          
          // Get tag of this particle photon from fragmentation, decay, prompt ...
          Int_t tag = GetMCAnalysisUtils()->CheckOrigin(i,GetReader(), 0);
          if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)){
            //A conversion photon from a hadron, skip this kind of photon
            //            printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!: tag %d, conv %d, pi0 %d, hadron %d, electron %d, unk %d, muon %d,pion %d, proton %d, neutron %d, kaon %d, antiproton %d, antineutron %d\n",tag,
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOther),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCUnknown),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCMuon), 
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPion),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCProton), 
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron),
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCKaon), 
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton), 
            //                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron));
            
            return;
          }
          
          //Get photon kinematics
          if(prim->Energy() == TMath::Abs(prim->Pz()))  continue ; //Protection against floating point exception	  
          
          Double_t photonY   = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
          Double_t photonE   = prim->Energy() ;
          Double_t photonPt  = prim->Pt() ;
          Double_t photonPhi = TMath::RadToDeg()*prim->Phi() ;
          if(photonPhi < 0) photonPhi+=TMath::TwoPi();
          Double_t photonEta = prim->Eta() ;
          
          
          //Check if photons hit the Calorimeter
          TLorentzVector lv;
          prim->Momentum(lv);
          Bool_t inacceptance = kFALSE;
          if     (fCalorimeter == "PHOS"){
            if(GetPHOSGeometry() && GetCaloUtils()->IsPHOSGeoMatrixSet()){
              Int_t mod ;
              Double_t x,z ;
              if(GetPHOSGeometry()->ImpactOnEmc(prim,mod,z,x)) 
                inacceptance = kTRUE;
              if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
            else{
              if(GetFiducialCut()->IsInFiducialCut(lv,fCalorimeter)) 
                inacceptance = kTRUE ;
              if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
          }	   
          else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
            if(GetEMCALGeometry()){
              
              Int_t absID=0;
              
              GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(prim->Eta(),prim->Phi(),absID);
              
              if( absID >= 0) 
                inacceptance = kTRUE;
              
              //                  if(GetEMCALGeometry()->Impact(phot1) && GetEMCALGeometry()->Impact(phot2)) 
              //                    inacceptance = kTRUE;
              if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
            else{
              if(GetFiducialCut()->IsInFiducialCut(lv,fCalorimeter)) 
                inacceptance = kTRUE ;
              if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
          }	  //In EMCAL
          
          //Fill histograms
          
          fhYPrimMC[kmcPPhoton]->Fill(photonPt, photonY) ;
          if(TMath::Abs(photonY) < 1.0)
          {
            fhEPrimMC  [kmcPPhoton]->Fill(photonE ) ;
            fhPtPrimMC [kmcPPhoton]->Fill(photonPt) ;
            fhPhiPrimMC[kmcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMC[kmcPPhoton]  ->Fill(photonE , photonEta) ;
          }
          if(inacceptance){
            fhEPrimMCAcc[kmcPPhoton]  ->Fill(photonE ) ;
            fhPtPrimMCAcc[kmcPPhoton] ->Fill(photonPt) ;
            fhPhiPrimMCAcc[kmcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMCAcc[kmcPPhoton]  ->Fill(photonE , photonY) ;
          }//Accepted
          
          //Origin of photon
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhEPrimMC[kmcPPrompt])
          {
            fhYPrimMC[kmcPPrompt]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPPrompt]->Fill(photonE ) ;
              fhPtPrimMC [kmcPPrompt]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPPrompt]  ->Fill(photonE , photonEta) ;
            }   
            if(inacceptance){
              fhEPrimMCAcc[kmcPPrompt]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPPrompt] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPPrompt]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) && fhEPrimMC[kmcPFragmentation])
          {
            fhYPrimMC[kmcPFragmentation]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPFragmentation]->Fill(photonE ) ;
              fhPtPrimMC [kmcPFragmentation]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPFragmentation]  ->Fill(photonE , photonEta) ;
            }  
            if(inacceptance){
              fhEPrimMCAcc[kmcPFragmentation]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPFragmentation] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPFragmentation]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) && fhEPrimMC[kmcPISR])
          {
            fhYPrimMC[kmcPISR]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPISR]->Fill(photonE ) ;
              fhPtPrimMC [kmcPISR]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPISR]->Fill(photonE , photonEta) ;
            }            
            if(inacceptance){
              fhEPrimMCAcc[kmcPISR]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPISR] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPISR]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)&& fhEPrimMC[kmcPPi0Decay])
          {
            fhYPrimMC[kmcPPi0Decay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPPi0Decay]->Fill(photonE ) ;
              fhPtPrimMC [kmcPPi0Decay]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPPi0Decay]  ->Fill(photonE , photonEta) ;
            }     
            if(inacceptance){
              fhEPrimMCAcc[kmcPPi0Decay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPPi0Decay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPPi0Decay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                    GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay)) && fhEPrimMC[kmcPOtherDecay])
          {
            fhYPrimMC[kmcPOtherDecay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPOtherDecay]->Fill(photonE ) ;
              fhPtPrimMC [kmcPOtherDecay]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPOtherDecay]  ->Fill(photonE , photonEta) ;
            } 
            if(inacceptance){
              fhEPrimMCAcc[kmcPOtherDecay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPOtherDecay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPOtherDecay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(fhEPrimMC[kmcPOther])
          {
            fhYPrimMC[kmcPOther]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPOther]->Fill(photonE ) ;
              fhPtPrimMC [kmcPOther]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPOther]  ->Fill(photonE , photonEta) ;
            }
            if(inacceptance){
              fhEPrimMCAcc[kmcPOther]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPOther] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPOther]  ->Fill(photonE , photonY) ;
            }//Accepted
          }//Other origin
        }// Primary photon 
      }//loop on primaries	
    }//stack exists and data is MC
  }//read stack
  else if(GetReader()->ReadAODMCParticles()){
    TClonesArray * mcparticles = GetReader()->GetAODMCParticles(0);
    if(mcparticles){
      Int_t nprim = mcparticles->GetEntriesFast();
      
      for(Int_t i=0; i < nprim; i++)
      {
        AliAODMCParticle * prim = (AliAODMCParticle *) mcparticles->At(i);   
        
        Int_t pdg = prim->GetPdgCode();
        
        if(pdg == 22){
          
          // Get tag of this particle photon from fragmentation, decay, prompt ...
          Int_t tag = GetMCAnalysisUtils()->CheckOrigin(i,GetReader(), 0);
          if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)){
            //A conversion photon from a hadron, skip this kind of photon
//            printf("AliAnaPhoton::FillAcceptanceHistograms() - not a photon, weird!: tag %d, conv %d, pi0 %d, hadron %d, electron %d, unk %d, muon %d,pion %d, proton %d, neutron %d, kaon %d, antiproton %d, antineutron %d\n",tag,
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOther),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCUnknown),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCMuon), 
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPion),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCProton), 
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron),
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCKaon), 
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton), 
//                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron));
            
            return;
          }
          
          //Get photon kinematics
          if(prim->E() == TMath::Abs(prim->Pz()))  continue ; //Protection against floating point exception	 
          
          Double_t photonY  = 0.5*TMath::Log((prim->E()-prim->Pz())/(prim->E()+prim->Pz())) ;
          Double_t photonE   = prim->E() ;
          Double_t photonPt  = prim->Pt() ;
          Double_t photonPhi = TMath::RadToDeg()*prim->Phi() ;
          if(photonPhi < 0) photonPhi+=TMath::TwoPi();
          Double_t photonEta = prim->Eta() ;
          
          //Check if photons hit the Calorimeter
          TLorentzVector lv;
          lv.SetPxPyPzE(prim->Px(),prim->Py(),prim->Pz(),prim->E());
          Bool_t inacceptance = kFALSE;
          if     (fCalorimeter == "PHOS"){
            if(GetPHOSGeometry() && GetCaloUtils()->IsPHOSGeoMatrixSet()){
              Int_t mod ;
              Double_t x,z ;
              Double_t vtx[]={prim->Xv(),prim->Yv(),prim->Zv()};
              if(GetPHOSGeometry()->ImpactOnEmc(vtx, prim->Theta(),prim->Phi(),mod,z,x))
                inacceptance = kTRUE;
              if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
            else{
              if(GetFiducialCut()->IsInFiducialCut(lv,fCalorimeter)) 
                inacceptance = kTRUE ;
              if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
          }	   
          else if(fCalorimeter == "EMCAL" && GetCaloUtils()->IsEMCALGeoMatrixSet()){
            if(GetEMCALGeometry()){
              
              Int_t absID=0;
              
              GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(prim->Eta(),prim->Phi(),absID);
              
              if( absID >= 0) 
                inacceptance = kTRUE;
              
              if(GetDebug() > 2) printf("In %s Real acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
            else{
              if(GetFiducialCut()->IsInFiducialCut(lv,fCalorimeter)) 
                inacceptance = kTRUE ;
              if(GetDebug() > 2) printf("In %s fiducial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
            }
          }	  //In EMCAL
          
          //Fill histograms
          
          fhYPrimMC[kmcPPhoton]->Fill(photonPt, photonY) ;
          if(TMath::Abs(photonY) < 1.0)
          {
            fhEPrimMC  [kmcPPhoton]->Fill(photonE ) ;
            fhPtPrimMC [kmcPPhoton]->Fill(photonPt) ;
            fhPhiPrimMC[kmcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMC[kmcPPhoton]->Fill(photonE , photonEta) ;
          }
          if(inacceptance){
            fhEPrimMCAcc[kmcPPhoton]  ->Fill(photonE ) ;
            fhPtPrimMCAcc[kmcPPhoton] ->Fill(photonPt) ;
            fhPhiPrimMCAcc[kmcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMCAcc[kmcPPhoton]  ->Fill(photonE , photonY) ;
          }//Accepted
          
          
          //Origin of photon
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhEPrimMC[kmcPPrompt])
          {
            fhYPrimMC[kmcPPrompt]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPPrompt]->Fill(photonE ) ;
              fhPtPrimMC [kmcPPrompt]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPPrompt]->Fill(photonE , photonEta) ;
            }   
            if(inacceptance){
              fhEPrimMCAcc[kmcPPrompt]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPPrompt] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPPrompt]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) && fhEPrimMC[kmcPFragmentation] )
          {
            fhYPrimMC[kmcPFragmentation]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPFragmentation]->Fill(photonE ) ;
              fhPtPrimMC [kmcPFragmentation]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPFragmentation]->Fill(photonE , photonEta) ;
            }  
            if(inacceptance){
              fhEPrimMCAcc[kmcPFragmentation]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPFragmentation] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPFragmentation]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) && fhEPrimMC[kmcPISR])
          {
            fhYPrimMC[kmcPISR]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPISR]->Fill(photonE ) ;
              fhPtPrimMC [kmcPISR]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPISR]->Fill(photonE , photonEta) ;
            }            
            if(inacceptance){
              fhEPrimMCAcc[kmcPISR]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPISR] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPISR]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)&& fhEPrimMC[kmcPPi0Decay])
          {
            fhYPrimMC[kmcPPi0Decay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPPi0Decay]->Fill(photonE ) ;
              fhPtPrimMC [kmcPPi0Decay]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPPi0Decay]->Fill(photonE , photonEta) ;
            }     
            if(inacceptance){
              fhEPrimMCAcc[kmcPPi0Decay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPPi0Decay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPPi0Decay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if((GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhEPrimMC[kmcPOtherDecay])
          {
            fhYPrimMC[kmcPOtherDecay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPOtherDecay]->Fill(photonE ) ;
              fhPtPrimMC [kmcPOtherDecay]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPOtherDecay]->Fill(photonE , photonEta) ;
            } 
            if(inacceptance){
              fhEPrimMCAcc[kmcPOtherDecay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPOtherDecay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPOtherDecay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(fhEPrimMC[kmcPOther])
          {
            fhYPrimMC[kmcPOther]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [kmcPOther]->Fill(photonE ) ;
              fhPtPrimMC [kmcPOther]->Fill(photonPt) ;
              fhPhiPrimMC[kmcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMC[kmcPOther]->Fill(photonE , photonEta) ;
            }
            if(inacceptance){
              fhEPrimMCAcc[kmcPOther]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[kmcPOther] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[kmcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[kmcPOther]  ->Fill(photonE , photonY) ;
            }//Accepted
          }//Other origin
        }// Primary photon 
      }//loop on primaries	
      
    }//kmc array exists and data is MC
  }	// read AOD MC
}

//__________________________________________________________________
void  AliAnaPhoton::FillShowerShapeHistograms(AliVCluster* cluster, const Int_t mcTag){
  
  //Fill cluster Shower Shape histograms
  
  if(!fFillSSHistograms || GetMixedEvent()) return;

  Float_t energy  = cluster->E();
  Int_t   ncells  = cluster->GetNCells();
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  Float_t disp    = cluster->GetDispersion()*cluster->GetDispersion();
  
  TLorentzVector mom;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
    cluster->GetMomentum(mom,GetVertex(0)) ;}//Assume that come from vertex in straight line
  else{
    Double_t vertex[]={0,0,0};
    cluster->GetMomentum(mom,vertex) ;
  }
  
  Float_t eta = mom.Eta();
  Float_t phi = mom.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  
  fhLam0E ->Fill(energy,lambda0);
  fhLam1E ->Fill(energy,lambda1);
  fhDispE ->Fill(energy,disp);
    
  if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5){
    fhLam0ETRD->Fill(energy,lambda0);
    fhLam1ETRD->Fill(energy,lambda1);
    fhDispETRD->Fill(energy,disp);
  }
  
  // if track-matching was of, check effect of track-matching residual cut 
  
  if(!fRejectTrackMatch)
  {
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn()){
      dR = 2000., dZ = 2000.;
      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }   
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {
      fhLam0ETM ->Fill(energy,lambda0);
      fhLam1ETM ->Fill(energy,lambda1);
      fhDispETM ->Fill(energy,disp);
      
      if(fCalorimeter == "EMCAL" && GetModuleNumber(cluster) > 5){
        fhLam0ETMTRD->Fill(energy,lambda0);
        fhLam1ETMTRD->Fill(energy,lambda1);
        fhDispETMTRD->Fill(energy,disp);
      }
    }
  }// if track-matching was of, check effect of matching residual cut  
  
  if(energy < 2){
    fhNCellsLam0LowE ->Fill(ncells,lambda0);
    fhNCellsLam1LowE ->Fill(ncells,lambda1);
    fhNCellsDispLowE ->Fill(ncells,disp);
    
    fhLam1Lam0LowE  ->Fill(lambda1,lambda0);
    fhLam0DispLowE  ->Fill(lambda0,disp);
    fhDispLam1LowE  ->Fill(disp,lambda1);
    fhEtaLam0LowE   ->Fill(eta,lambda0);
    fhPhiLam0LowE   ->Fill(phi,lambda0);  
    
  }
  else {
    fhNCellsLam0HighE ->Fill(ncells,lambda0);
    fhNCellsLam1HighE ->Fill(ncells,lambda1);
    fhNCellsDispHighE ->Fill(ncells,disp);

    fhLam1Lam0HighE  ->Fill(lambda1,lambda0);
    fhLam0DispHighE  ->Fill(lambda0,disp);
    fhDispLam1HighE  ->Fill(disp,lambda1);
    fhEtaLam0HighE   ->Fill(eta, lambda0);
    fhPhiLam0HighE   ->Fill(phi, lambda0);
  }
  
  if(IsDataMC()){

    AliVCaloCells* cells = 0;
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    //Fill histograms to check shape of embedded clusters
    Float_t fraction = 0;
    if(GetReader()->IsEmbeddedClusterSelectionOn()){//Only working for EMCAL

      Float_t clusterE = 0; // recalculate in case corrections applied.
      Float_t cellE    = 0;
      for(Int_t icell = 0; icell < cluster->GetNCells(); icell++){
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
    Int_t absID             =-1 ;
    Float_t maxCellFraction = 0.;
    
    absID = GetCaloUtils()->GetMaxEnergyCell(cells, cluster,maxCellFraction);
    
    // Check the origin and fill histograms
    if( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) && 
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0) &&
       !GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)){
      fhMCELambda0[kmcssPhoton]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssPhoton]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssPhoton] ->Fill(energy, disp);
      fhMCNCellsE[kmcssPhoton]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssPhoton]->Fill(energy,maxCellFraction);
            
      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      
      if(!GetReader()->IsEmbeddedClusterSelectionOn()){
        //Check particle overlaps in cluster
        
        // Compare the primary depositing more energy with the rest, 
        // if no photon/electron as comon ancestor (conversions), count as other particle
        Int_t ancPDG = 0, ancStatus = -1;
        TLorentzVector momentum; TVector3 prodVertex;
        Int_t ancLabel = 0;
        Int_t noverlaps = 1;      
        for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ ) {
          ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(cluster->GetLabels()[0],cluster->GetLabels()[ilab], GetReader(),ancPDG,ancStatus,momentum,prodVertex);
          if(ancPDG!=22 && TMath::Abs(ancPDG)!=11) noverlaps++;
        }
        //printf("N overlaps %d \n",noverlaps);
        
        if(noverlaps == 1){
          fhMCPhotonELambda0NoOverlap  ->Fill(energy, lambda0);
        }
        else if(noverlaps == 2){        
          fhMCPhotonELambda0TwoOverlap ->Fill(energy, lambda0);
        }
        else if(noverlaps > 2){          
          fhMCPhotonELambda0NOverlap   ->Fill(energy, lambda0);
        }
        else {
          printf("AliAnaPhoton::FillShowerShapeHistogram() - n overlaps = %d!!", noverlaps);
        }
      }//No embedding
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn()){
        
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
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCElectron)){
      fhMCELambda0[kmcssElectron]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssElectron]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssElectron] ->Fill(energy, disp);
      fhMCNCellsE[kmcssElectron]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssElectron]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssElectron] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssElectron] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssElectron] ->Fill(ncells,  maxCellFraction);
      }
    }//electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) && 
              GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) ){
      fhMCELambda0[kmcssConversion]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssConversion]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssConversion] ->Fill(energy, disp);
      fhMCNCellsE[kmcssConversion]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssConversion]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssConversion] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssConversion] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssConversion] ->Fill(ncells,  maxCellFraction);
      }      
      
    }//conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  ){
      fhMCELambda0[kmcssPi0]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssPi0]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssPi0] ->Fill(energy, disp);
      fhMCNCellsE[kmcssPi0]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssPi0]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssPi0] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssPi0] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssPi0] ->Fill(ncells,  maxCellFraction);
      }      
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn()){
        
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
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCEta)  ){
      fhMCELambda0[kmcssEta]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssEta]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssEta] ->Fill(energy, disp);
      fhMCNCellsE[kmcssEta]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssEta]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssEta] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssEta] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssEta] ->Fill(ncells,  maxCellFraction);
      }
      
    }//eta    
    else {
      fhMCELambda0[kmcssOther]    ->Fill(energy, lambda0);
      fhMCELambda1[kmcssOther]    ->Fill(energy, lambda1);
      fhMCEDispersion[kmcssOther] ->Fill(energy, disp);
      fhMCNCellsE[kmcssOther]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[kmcssOther]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[kmcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[kmcssOther] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[kmcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[kmcssOther] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[kmcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[kmcssOther] ->Fill(ncells,  maxCellFraction);
      }            
      
    }//other particles 
    
  }//MC data
  
}

//__________________________________________________________________________
void AliAnaPhoton::FillTrackMatchingResidualHistograms(AliVCluster* cluster, 
                                                       const Int_t cut)
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
    
  if(fhTrackMatchedDEta[cut] && TMath::Abs(dR) < 999)
  {
    fhTrackMatchedDEta[cut]->Fill(cluster->E(),dZ);
    fhTrackMatchedDPhi[cut]->Fill(cluster->E(),dR);
    
    if(cluster->E() > 0.5) fhTrackMatchedDEtaDPhi[cut]->Fill(dZ,dR);
    
    Int_t nSMod = GetModuleNumber(cluster);
    
    if(fCalorimeter=="EMCAL" &&  nSMod > 5)
    {
      fhTrackMatchedDEtaTRD[cut]->Fill(cluster->E(),dZ);
      fhTrackMatchedDPhiTRD[cut]->Fill(cluster->E(),dR);
    }
    
    // Check dEdx and E/p of matched clusters
    
    if(TMath::Abs(dZ) < 0.05 && TMath::Abs(dR) < 0.05)
    {      
      
      AliVTrack *track = GetCaloUtils()->GetMatchedTrack(cluster, GetReader()->GetInputEvent());
      
      if(track) 
      {
        
        Float_t dEdx   = track->GetTPCsignal();
        Float_t eOverp = cluster->E()/track->P();
        
        fhdEdx[cut]  ->Fill(cluster->E(), dEdx);
        fhEOverP[cut]->Fill(cluster->E(), eOverp);
        
        if(fCalorimeter=="EMCAL" && nSMod > 5)
          fhEOverPTRD[cut]->Fill(cluster->E(), eOverp);
        
        
      }
      else
          printf("AliAnaPhoton::FillTrackMatchingResidualHistograms() - Residual OK but (dR, dZ)= (%2.4f,%2.4f) no track associated WHAT? \n", dR,dZ);
      
      
      
      if(IsDataMC())
      {
        
        Int_t tag = GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(),GetReader(), 0);
        
        if  ( !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)  )
        {
          if       ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)      ||
                     GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)       ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 2.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton)    ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 0.5 );
          else if  ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)  ) fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 1.5 );
          else                                                                                 fhTrackMatchedMCParticle[cut]->Fill(cluster->E(), 3.5 );
          
          // Check if several particles contributed to cluster and discard overlapped mesons
          if(!GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) || 
             !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)){
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
             !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)){
            
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
  
  TString cut[] = {"Open","Reader","E","Time","NCells","Fidutial","Matching","Bad","PID"};
  for (Int_t i = 0; i < 9 ;  i++) 
  {
    fhClusterCuts[i] = new TH1F(Form("hCut_%d_%s", i, cut[i].Data()),
                                Form("Number of clusters that pass cuts <= %d, %s", i, cut[i].Data()),
                                nptbins,ptmin,ptmax); 
    fhClusterCuts[i]->SetYTitle("dN/dE ");
    fhClusterCuts[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhClusterCuts[i]) ;   
  }
  
  fhNCellsE  = new TH2F ("hNCellsE","# of cells in cluster vs E of clusters", nptbins,ptmin,ptmax, nbins,nmin,nmax); 
  fhNCellsE->SetXTitle("E (GeV)");
  fhNCellsE->SetYTitle("# of cells in cluster");
  outputContainer->Add(fhNCellsE);    
  
  fhCellsE  = new TH2F ("hCellsE","energy of cells in cluster vs E of clusters", nptbins,ptmin,ptmax, nptbins*2,ptmin,ptmax); 
  fhCellsE->SetXTitle("E_{cluster} (GeV)");
  fhCellsE->SetYTitle("E_{cell} (GeV)");
  outputContainer->Add(fhCellsE);    
  
  fhTimeE  = new TH2F ("hTimeE","time of cluster vs E of clusters", nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
  fhTimeE->SetXTitle("E (GeV)");
  fhTimeE->SetYTitle("time (ns)");
  outputContainer->Add(fhTimeE);  
  
  fhMaxCellDiffClusterE  = new TH2F ("hMaxCellDiffClusterE","energy vs difference of cluster energy - max cell energy / cluster energy, good clusters",
                                    nptbins,ptmin,ptmax, 500,0,1.); 
  fhMaxCellDiffClusterE->SetXTitle("E_{cluster} (GeV) ");
  fhMaxCellDiffClusterE->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
  outputContainer->Add(fhMaxCellDiffClusterE);  
  
  fhEPhoton  = new TH1F("hEPhoton","Number of #gamma over calorimeter vs energy",nptbins,ptmin,ptmax); 
  fhEPhoton->SetYTitle("N");
  fhEPhoton->SetXTitle("E_{#gamma}(GeV)");
  outputContainer->Add(fhEPhoton) ;   
  
  fhPtPhoton  = new TH1F("hPtPhoton","Number of #gamma over calorimeter vs p_{T}",nptbins,ptmin,ptmax); 
  fhPtPhoton->SetYTitle("N");
  fhPtPhoton->SetXTitle("p_{T #gamma}(GeV/c)");
  outputContainer->Add(fhPtPhoton) ; 
  
  fhPhiPhoton  = new TH2F
    ("hPhiPhoton","#phi_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
  fhPhiPhoton->SetYTitle("#phi (rad)");
  fhPhiPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhiPhoton) ; 
  
  fhEtaPhoton  = new TH2F
    ("hEtaPhoton","#eta_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,netabins,etamin,etamax); 
  fhEtaPhoton->SetYTitle("#eta");
  fhEtaPhoton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhEtaPhoton) ;
  
  fhEtaPhiPhoton  = new TH2F
  ("hEtaPhiPhoton","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
  fhEtaPhiPhoton->SetYTitle("#phi (rad)");
  fhEtaPhiPhoton->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiPhoton) ;
  if(GetMinPt() < 0.5){
    fhEtaPhi05Photon  = new TH2F
    ("hEtaPhi05Photon","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi05Photon->SetYTitle("#phi (rad)");
    fhEtaPhi05Photon->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhi05Photon) ;
  }
  
  //Shower shape
  if(fFillSSHistograms){
    
    fhLam0E  = new TH2F ("hLam0E","#lambda_{0}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLam0E->SetYTitle("#lambda_{0}^{2}");
    fhLam0E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam0E);  
    
    fhLam1E  = new TH2F ("hLam1E","#lambda_{1}^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhLam1E->SetYTitle("#lambda_{1}^{2}");
    fhLam1E->SetXTitle("E (GeV)");
    outputContainer->Add(fhLam1E);  
    
    fhDispE  = new TH2F ("hDispE"," dispersion^{2} vs E", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhDispE->SetYTitle("D^{2}");
    fhDispE->SetXTitle("E (GeV) ");
    outputContainer->Add(fhDispE);

    if(!fRejectTrackMatch)
    {
      fhLam0ETM  = new TH2F ("hLam0ETM","#lambda_{0}^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam0ETM->SetYTitle("#lambda_{0}^{2}");
      fhLam0ETM->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam0ETM);  
      
      fhLam1ETM  = new TH2F ("hLam1ETM","#lambda_{1}^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam1ETM->SetYTitle("#lambda_{1}^{2}");
      fhLam1ETM->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam1ETM);  
      
      fhDispETM  = new TH2F ("hDispETM"," dispersion^{2} vs E, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhDispETM->SetYTitle("D^{2}");
      fhDispETM->SetXTitle("E (GeV) ");
      outputContainer->Add(fhDispETM);
    }
    
    if(fCalorimeter == "EMCAL")
    {
      fhLam0ETRD  = new TH2F ("hLam0ETRD","#lambda_{0}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam0ETRD->SetYTitle("#lambda_{0}^{2}");
      fhLam0ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam0ETRD);  
      
      fhLam1ETRD  = new TH2F ("hLam1ETRD","#lambda_{1}^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhLam1ETRD->SetYTitle("#lambda_{1}^{2}");
      fhLam1ETRD->SetXTitle("E (GeV)");
      outputContainer->Add(fhLam1ETRD);  
      
      fhDispETRD  = new TH2F ("hDispETRD"," dispersion^{2} vs E, EMCAL SM covered by TRD", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhDispETRD->SetYTitle("Dispersion^{2}");
      fhDispETRD->SetXTitle("E (GeV) ");
      outputContainer->Add(fhDispETRD);
      
      if(!fRejectTrackMatch)
      {
        fhLam0ETMTRD  = new TH2F ("hLam0ETMTRD","#lambda_{0}^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam0ETMTRD->SetYTitle("#lambda_{0}^{2}");
        fhLam0ETMTRD->SetXTitle("E (GeV)");
        outputContainer->Add(fhLam0ETMTRD);  
        
        fhLam1ETMTRD  = new TH2F ("hLam1ETMTRD","#lambda_{1}^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhLam1ETMTRD->SetYTitle("#lambda_{1}^{2}");
        fhLam1ETMTRD->SetXTitle("E (GeV)");
        outputContainer->Add(fhLam1ETMTRD);  
        
        fhDispETMTRD  = new TH2F ("hDispETMTRD"," dispersion^{2} vs E, EMCAL SM covered by TRD, cut on track-matching residual |#Delta #eta| < 0.05,  |#Delta #phi| < 0.05", nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhDispETMTRD->SetYTitle("Dispersion^{2}");
        fhDispETMTRD->SetXTitle("E (GeV) ");
        outputContainer->Add(fhDispETMTRD);         
      } 
    } 
    
    fhNCellsLam0LowE  = new TH2F ("hNCellsLam0LowE","N_{cells} in cluster vs #lambda_{0}^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax); 
    fhNCellsLam0LowE->SetXTitle("N_{cells}");
    fhNCellsLam0LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam0LowE);  
    
    fhNCellsLam0HighE  = new TH2F ("hNCellsLam0HighE","N_{cells} in cluster vs #lambda_{0}^{2}, E > 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax); 
    fhNCellsLam0HighE->SetXTitle("N_{cells}");
    fhNCellsLam0HighE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam0HighE);  
    
    fhNCellsLam1LowE  = new TH2F ("hNCellsLam1LowE","N_{cells} in cluster vs #lambda_{1}^{2}, E < 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax); 
    fhNCellsLam1LowE->SetXTitle("N_{cells}");
    fhNCellsLam1LowE->SetYTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhNCellsLam1LowE);  
    
    fhNCellsLam1HighE  = new TH2F ("hNCellsLam1HighE","N_{cells} in cluster vs #lambda_{1}^{2}, E > 2 GeV", nbins,nmin, nmax, ssbins,ssmin,ssmax); 
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
    
    fhEtaLam0HighE  = new TH2F ("hEtaLam0HighE","#eta vs #lambda_{0}^{2}, E > 2 GeV", netabins,etamin,etamax, ssbins,ssmin,ssmax); 
    fhEtaLam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhEtaLam0HighE->SetXTitle("#eta");
    outputContainer->Add(fhEtaLam0HighE);  
    
    fhPhiLam0HighE  = new TH2F ("hPhiLam0HighE","#phi vs #lambda_{0}^{2}, E > 2 GeV", nphibins,phimin,phimax, ssbins,ssmin,ssmax); 
    fhPhiLam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhPhiLam0HighE->SetXTitle("#phi");
    outputContainer->Add(fhPhiLam0HighE);  
    
    fhLam1Lam0LowE  = new TH2F ("hLam1Lam0LowE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam1Lam0LowE->SetYTitle("#lambda_{0}^{2}");
    fhLam1Lam0LowE->SetXTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhLam1Lam0LowE);  
    
    fhLam1Lam0HighE  = new TH2F ("hLam1Lam0HighE","#lambda_{0}^{2} vs #lambda_{1}^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam1Lam0HighE->SetYTitle("#lambda_{0}^{2}");
    fhLam1Lam0HighE->SetXTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhLam1Lam0HighE);  
    
    fhLam0DispLowE  = new TH2F ("hLam0DispLowE","#lambda_{0}^{2} vs dispersion^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam0DispLowE->SetXTitle("#lambda_{0}^{2}");
    fhLam0DispLowE->SetYTitle("D^{2}");
    outputContainer->Add(fhLam0DispLowE);  
    
    fhLam0DispHighE  = new TH2F ("hLam0DispHighE","#lambda_{0}^{2} vs dispersion^{2} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhLam0DispHighE->SetXTitle("#lambda_{0}^{2}");
    fhLam0DispHighE->SetYTitle("D^{2}");
    outputContainer->Add(fhLam0DispHighE);  
    
    fhDispLam1LowE  = new TH2F ("hDispLam1LowE","Dispersion^{2} vs #lambda_{1}^{2} in cluster of E < 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhDispLam1LowE->SetXTitle("D^{2}");
    fhDispLam1LowE->SetYTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhDispLam1LowE);  
    
    fhDispLam1HighE  = new TH2F ("hDispLam1HighE","Dispersion^{2} vs #lambda_{1^{2}} in cluster of E > 2 GeV",  ssbins,ssmin,ssmax, ssbins,ssmin,ssmax); 
    fhDispLam1HighE->SetXTitle("D^{2}");
    fhDispLam1HighE->SetYTitle("#lambda_{1}^{2}");
    outputContainer->Add(fhDispLam1HighE);  
    
  } // Shower shape
  
  // Track Matching
  
  if(fFillTMHisto)
  {
        
    fhTrackMatchedDEta[0]  = new TH2F
    ("hTrackMatchedDEtaNoCut",
     "d#eta of cluster-track vs cluster energy, no photon cuts",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
    fhTrackMatchedDEta[0]->SetYTitle("d#eta");
    fhTrackMatchedDEta[0]->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDPhi[0]  = new TH2F
    ("hTrackMatchedDPhiNoCut",
     "d#phi of cluster-track vs cluster energy, no photon cuts",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
    fhTrackMatchedDPhi[0]->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhi[0]->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhi[0]  = new TH2F
    ("hTrackMatchedDEtaDPhiNoCut",
     "d#eta vs d#phi of cluster-track vs cluster energy, no photon cuts",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax); 
    fhTrackMatchedDEtaDPhi[0]->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhi[0]->SetXTitle("d#eta");   
        
    fhdEdx[0]  = new TH2F ("hdEdxNoCut","matched track <dE/dx> vs cluster E, no photon cuts ", 
                           nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
    fhdEdx[0]->SetXTitle("E (GeV)");
    fhdEdx[0]->SetYTitle("<dE/dx>");
    
    fhEOverP[0]  = new TH2F ("hEOverPNoCut","matched track E/p vs cluster E, no photon cuts ", 
                             nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
    fhEOverP[0]->SetXTitle("E (GeV)");
    fhEOverP[0]->SetYTitle("E/p");
    
    outputContainer->Add(fhTrackMatchedDEta[0]) ; 
    outputContainer->Add(fhTrackMatchedDPhi[0]) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhi[0]) ;
    outputContainer->Add(fhdEdx[0]);  
    outputContainer->Add(fhEOverP[0]);  

    fhTrackMatchedDEta[1]  = new TH2F
    ("hTrackMatchedDEta",
     "d#eta of cluster-track vs cluster energy, no photon cuts",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
    fhTrackMatchedDEta[1]->SetYTitle("d#eta");
    fhTrackMatchedDEta[1]->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDPhi[1]  = new TH2F
    ("hTrackMatchedDPhi",
     "d#phi of cluster-track vs cluster energy, no photon cuts",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
    fhTrackMatchedDPhi[1]->SetYTitle("d#phi (rad)");
    fhTrackMatchedDPhi[1]->SetXTitle("E_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhi[1]  = new TH2F
    ("hTrackMatchedDEtaDPhi",
     "d#eta vs d#phi of cluster-track vs cluster energy, no photon cuts",
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax); 
    fhTrackMatchedDEtaDPhi[1]->SetYTitle("d#phi (rad)");
    fhTrackMatchedDEtaDPhi[1]->SetXTitle("d#eta");   
    
    fhdEdx[1]  = new TH2F ("hdEdx","matched track <dE/dx> vs cluster E ", 
                           nptbins,ptmin,ptmax,ndedxbins, dedxmin, dedxmax); 
    fhdEdx[1]->SetXTitle("E (GeV)");
    fhdEdx[1]->SetYTitle("<dE/dx>");
    
    fhEOverP[1]  = new TH2F ("hEOverP","matched track E/p vs cluster E ", 
                             nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
    fhEOverP[1]->SetXTitle("E (GeV)");
    fhEOverP[1]->SetYTitle("E/p");
    
    outputContainer->Add(fhTrackMatchedDEta[1]) ; 
    outputContainer->Add(fhTrackMatchedDPhi[1]) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhi[1]) ;
    outputContainer->Add(fhdEdx[1]);  
    outputContainer->Add(fhEOverP[1]);      
    
    if(fCalorimeter=="EMCAL")
    {
      
      fhTrackMatchedDEtaTRD[0]  = new TH2F
      ("hTrackMatchedDEtaTRDNoCut",
       "d#eta of cluster-track vs cluster energy, SM behind TRD, no photon cuts",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaTRD[0]->SetYTitle("d#eta");
      fhTrackMatchedDEtaTRD[0]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiTRD[0]  = new TH2F
      ("hTrackMatchedDPhiTRDNoCut",
       "d#phi of cluster-track vs cluster energy, SM behing TRD, no photon cuts",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiTRD[0]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiTRD[0]->SetXTitle("E_{cluster} (GeV)");
            
      fhEOverPTRD[0]  = new TH2F ("hEOverPTRDNoCut","matched track E/p vs cluster E, behind TRD, no photon cuts ", 
                                  nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
      fhEOverPTRD[0]->SetXTitle("E (GeV)");
      fhEOverPTRD[0]->SetYTitle("E/p");
      
      outputContainer->Add(fhTrackMatchedDEtaTRD[0]) ; 
      outputContainer->Add(fhTrackMatchedDPhiTRD[0]) ;
      outputContainer->Add(fhEOverPTRD[0]);  
      
      fhTrackMatchedDEtaTRD[1]  = new TH2F
      ("hTrackMatchedDEtaTRD",
       "d#eta of cluster-track vs cluster energy, SM behind TRD",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaTRD[1]->SetYTitle("d#eta");
      fhTrackMatchedDEtaTRD[1]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiTRD[1]  = new TH2F
      ("hTrackMatchedDPhiTRD",
       "d#phi of cluster-track vs cluster energy, SM behing TRD",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiTRD[1]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiTRD[1]->SetXTitle("E_{cluster} (GeV)");
      
      fhEOverPTRD[1]  = new TH2F ("hEOverPTRD","matched track E/p vs cluster E, behind TRD ", 
                                  nptbins,ptmin,ptmax,nPoverEbins,pOverEmin,pOverEmax); 
      fhEOverPTRD[1]->SetXTitle("E (GeV)");
      fhEOverPTRD[1]->SetYTitle("E/p");
      
      outputContainer->Add(fhTrackMatchedDEtaTRD[1]) ; 
      outputContainer->Add(fhTrackMatchedDPhiTRD[1]) ;
      outputContainer->Add(fhEOverPTRD[1]);  
      
    }
    
    if(IsDataMC())
    {
      
      fhTrackMatchedDEtaMCNoOverlap[0]  = new TH2F
      ("hTrackMatchedDEtaMCNoOverlapNoCut",
       "d#eta of cluster-track vs cluster energy, no other MC particles overlap",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCNoOverlap[0]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCNoOverlap[0]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCNoOverlap[0]  = new TH2F
      ("hTrackMatchedDPhiMCNoOverlapNoCut",
       "d#phi of cluster-track vs cluster energy, no other MC particles overlap",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCNoOverlap[0]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCNoOverlap[0]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCNoOverlap[0]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCNoOverlap[0]) ;
      
      fhTrackMatchedDEtaMCNoOverlap[1]  = new TH2F
      ("hTrackMatchedDEtaMCNoOverlap",
       "d#eta of cluster-track vs cluster energy, no other MC particles overlap",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCNoOverlap[1]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCNoOverlap[1]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCNoOverlap[1]  = new TH2F
      ("hTrackMatchedDPhiMCNoOverlap",
       "d#phi of cluster-track vs cluster energy, no other MC particles overlap",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCNoOverlap[1]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCNoOverlap[1]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCNoOverlap[1]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCNoOverlap[1]) ;
      
      fhTrackMatchedDEtaMCOverlap[0]  = new TH2F
      ("hTrackMatchedDEtaMCOverlapNoCut",
       "d#eta of cluster-track vs cluster energy, several MC particles overlap",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCOverlap[0]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCOverlap[0]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCOverlap[0]  = new TH2F
      ("hTrackMatchedDPhiMCOverlapNoCut",
       "d#phi of cluster-track vs cluster energy, several MC particles overlap",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCOverlap[0]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCOverlap[0]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCOverlap[0]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCOverlap[0]) ;
      
      fhTrackMatchedDEtaMCOverlap[1]  = new TH2F
      ("hTrackMatchedDEtaMCOverlap",
       "d#eta of cluster-track vs cluster energy, several MC particles overlap",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCOverlap[1]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCOverlap[1]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCOverlap[1]  = new TH2F
      ("hTrackMatchedDPhiMCOverlap",
       "d#phi of cluster-track vs cluster energy, several MC particles overlap",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCOverlap[1]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCOverlap[1]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCOverlap[1]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCOverlap[1]) ;      
      
      fhTrackMatchedDEtaMCConversion[0]  = new TH2F
      ("hTrackMatchedDEtaMCConversionNoCut",
       "d#eta of cluster-track vs cluster energy, no other MC particles overlap appart from conversions",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCConversion[0]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCConversion[0]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCConversion[0]  = new TH2F
      ("hTrackMatchedDPhiMCConversionNoCut",
       "d#phi of cluster-track vs cluster energy, no other MC particles overlap appart from conversions",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCConversion[0]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCConversion[0]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCConversion[0]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCConversion[0]) ;
       
      
      fhTrackMatchedDEtaMCConversion[1]  = new TH2F
      ("hTrackMatchedDEtaMCConversion",
       "d#eta of cluster-track vs cluster energy, no other MC particles overlap appart from conversions",
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaMCConversion[1]->SetYTitle("d#eta");
      fhTrackMatchedDEtaMCConversion[1]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiMCConversion[1]  = new TH2F
      ("hTrackMatchedDPhiMCConversion",
       "d#phi of cluster-track vs cluster energy, no other MC particles overlap appart from conversions",
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiMCConversion[1]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiMCConversion[1]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaMCConversion[1]) ; 
      outputContainer->Add(fhTrackMatchedDPhiMCConversion[1]) ;
      
      
      fhTrackMatchedMCParticle[0]  = new TH2F
      ("hTrackMatchedMCParticleNoCut",
       "Origin of particle vs energy",
       nptbins,ptmin,ptmax,8,0,8); 
      fhTrackMatchedMCParticle[0]->SetXTitle("E (GeV)");   
      //fhTrackMatchedMCParticle[0]->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticle[0]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");
      
      fhTrackMatchedMCParticle[1]  = new TH2F
      ("hTrackMatchedMCParticle",
       "Origin of particle vs energy",
       nptbins,ptmin,ptmax,8,0,8); 
      fhTrackMatchedMCParticle[1]->SetXTitle("E (GeV)");   
      //fhTrackMatchedMCParticle[1]->SetYTitle("Particle type");
      
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(1 ,"Photon");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(2 ,"Electron");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(3 ,"Meson Merged");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(4 ,"Rest");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(5 ,"Conv. Photon");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(6 ,"Conv. Electron");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(7 ,"Conv. Merged");
      fhTrackMatchedMCParticle[1]->GetYaxis()->SetBinLabel(8 ,"Conv. Rest");      
      
      outputContainer->Add(fhTrackMatchedMCParticle[0]);            
      outputContainer->Add(fhTrackMatchedMCParticle[1]);       
      
    }
  }  
  
  
  if(IsDataMC()){
   
    TString ptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}", "#pi^{0}","#eta",
                        "e^{#pm}","#gamma->e^{#pm}","hadron?","Anti-N","Anti-P",
                        "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}","String"                                } ; 
    
    TString pname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Pi0","Eta","Electron",
                        "Conversion", "Hadron", "AntiNeutron","AntiProton",
                        "PhotonPrompt","PhotonFragmentation","PhotonISR","String" } ;
    
    for(Int_t i = 0; i < fNOriginHistograms; i++){ 
      
      fhMCE[i]  = new TH1F(Form("hE_MC%s",pname[i].Data()),
                                Form("cluster from %s : E ",ptype[i].Data()),
                                nptbins,ptmin,ptmax); 
      fhMCE[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCE[i]) ; 
      
      fhMCPt[i]  = new TH1F(Form("hPt_MC%s",pname[i].Data()),
                           Form("cluster from %s : p_{T} ",ptype[i].Data()),
                           nptbins,ptmin,ptmax); 
      fhMCPt[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhMCPt[i]) ;
      
      fhMCEta[i]  = new TH2F(Form("hEta_MC%s",pname[i].Data()),
                           Form("cluster from %s : #eta ",ptype[i].Data()),
                           nptbins,ptmin,ptmax,netabins,etamin,etamax); 
      fhMCEta[i]->SetYTitle("#eta");
      fhMCEta[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCEta[i]) ;
      
      fhMCPhi[i]  = new TH2F(Form("hPhi_MC%s",pname[i].Data()),
                           Form("cluster from %s : #phi ",ptype[i].Data()),
                           nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhMCPhi[i]->SetYTitle("#phi (rad)");
      fhMCPhi[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMCPhi[i]) ;
      
      
      fhMCDeltaE[i]  = new TH2F (Form("hDeltaE_MC%s",pname[i].Data()),
                                 Form("MC - Reco E from %s",pname[i].Data()), 
                                 nptbins,ptmin,ptmax, 200,-50,50); 
      fhMCDeltaE[i]->SetXTitle("#Delta E (GeV)");
      outputContainer->Add(fhMCDeltaE[i]);
      
      fhMCDeltaPt[i]  = new TH2F (Form("hDeltaPt_MC%s",pname[i].Data()),
                                  Form("MC - Reco p_{T} from %s",pname[i].Data()), 
                                  nptbins,ptmin,ptmax, 200,-50,50); 
      fhMCDeltaPt[i]->SetXTitle("#Delta p_{T} (GeV/c)");
      outputContainer->Add(fhMCDeltaPt[i]);
            
      fhMC2E[i]  = new TH2F (Form("h2E_MC%s",pname[i].Data()),
                             Form("E distribution, reconstructed vs generated from %s",pname[i].Data()), 
                             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
      fhMC2E[i]->SetXTitle("E_{rec} (GeV)");
      fhMC2E[i]->SetYTitle("E_{gen} (GeV)");
      outputContainer->Add(fhMC2E[i]);          
      
      fhMC2Pt[i]  = new TH2F (Form("h2Pt_MC%s",pname[i].Data()),
                              Form("p_T distribution, reconstructed vs generated from %s",pname[i].Data()), 
                              nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
      fhMC2Pt[i]->SetXTitle("p_{T,rec} (GeV/c)");
      fhMC2Pt[i]->SetYTitle("p_{T,gen} (GeV/c)");
      outputContainer->Add(fhMC2Pt[i]);
      
      
    }
    
    TString pptype[] = { "#gamma", "#gamma_{#pi decay}","#gamma_{other decay}","hadron?",
                         "#gamma_{prompt}","#gamma_{fragmentation}","#gamma_{ISR}"} ; 
    
    TString ppname[] = { "Photon","PhotonPi0Decay","PhotonOtherDecay","Hadron",
                         "PhotonPrompt","PhotonFragmentation","PhotonISR"} ;
    
    for(Int_t i = 0; i < fNPrimaryHistograms; i++){ 
      fhEPrimMC[i]  = new TH1F(Form("hEPrim_MC%s",ppname[i].Data()),
                           Form("primary photon %s : E ",pptype[i].Data()),
                           nptbins,ptmin,ptmax); 
      fhEPrimMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEPrimMC[i]) ; 
      
      fhPtPrimMC[i]  = new TH1F(Form("hPtPrim_MC%s",ppname[i].Data()),
                            Form("primary photon %s : p_{T} ",pptype[i].Data()),
                            nptbins,ptmin,ptmax); 
      fhPtPrimMC[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtPrimMC[i]) ;
      
      fhYPrimMC[i]  = new TH2F(Form("hYPrim_MC%s",ppname[i].Data()),
                             Form("primary photon %s : Rapidity ",pptype[i].Data()),
                             nptbins,ptmin,ptmax,800,-8,8); 
      fhYPrimMC[i]->SetYTitle("Rapidity");
      fhYPrimMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhYPrimMC[i]) ;
      
      fhPhiPrimMC[i]  = new TH2F(Form("hPhiPrim_MC%s",ppname[i].Data()),
                             Form("primary photon %s : #phi ",pptype[i].Data()),
                             nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiPrimMC[i]->SetYTitle("#phi (rad)");
      fhPhiPrimMC[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhPhiPrimMC[i]) ;
     
      
      fhEPrimMCAcc[i]  = new TH1F(Form("hEPrimAcc_MC%s",ppname[i].Data()),
                               Form("primary photon %s in acceptance: E ",pptype[i].Data()),
                               nptbins,ptmin,ptmax); 
      fhEPrimMCAcc[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhEPrimMCAcc[i]) ; 
      
      fhPtPrimMCAcc[i]  = new TH1F(Form("hPtPrimAcc_MC%s",ppname[i].Data()),
                                Form("primary photon %s in acceptance: p_{T} ",pptype[i].Data()),
                                nptbins,ptmin,ptmax); 
      fhPtPrimMCAcc[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtPrimMCAcc[i]) ;
      
      fhYPrimMCAcc[i]  = new TH2F(Form("hYPrimAcc_MC%s",ppname[i].Data()),
                                 Form("primary photon %s in acceptance: Rapidity ",pptype[i].Data()),
                                 nptbins,ptmin,ptmax,100,-1,1); 
      fhYPrimMCAcc[i]->SetYTitle("Rapidity");
      fhYPrimMCAcc[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhYPrimMCAcc[i]) ;
      
      fhPhiPrimMCAcc[i]  = new TH2F(Form("hPhiPrimAcc_MC%s",ppname[i].Data()),
                                 Form("primary photon %s in acceptance: #phi ",pptype[i].Data()),
                                 nptbins,ptmin,ptmax,nphibins,phimin,phimax); 
      fhPhiPrimMCAcc[i]->SetYTitle("#phi (rad)");
      fhPhiPrimMCAcc[i]->SetXTitle("E (GeV)");
      outputContainer->Add(fhPhiPrimMCAcc[i]) ;
      
    }
        	
    if(fFillSSHistograms){
      
      TString ptypess[] = { "#gamma","hadron?","#pi^{0}","#eta","#gamma->e^{#pm}","e^{#pm}"} ; 
      
      TString pnamess[] = { "Photon","Hadron","Pi0","Eta","Conversion","Electron"} ;
      
      for(Int_t i = 0; i < 6; i++){ 
        
        fhMCELambda0[i]  = new TH2F(Form("hELambda0_MC%s",pnamess[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{0}^{2}",ptypess[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCELambda0[i]->SetYTitle("#lambda_{0}^{2}");
        fhMCELambda0[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCELambda0[i]) ; 
        
        fhMCELambda1[i]  = new TH2F(Form("hELambda1_MC%s",pnamess[i].Data()),
                                    Form("cluster from %s : E vs #lambda_{1}^{2}",ptypess[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCELambda1[i]->SetYTitle("#lambda_{1}^{2}");
        fhMCELambda1[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCELambda1[i]) ; 
              
        fhMCEDispersion[i]  = new TH2F(Form("hEDispersion_MC%s",pnamess[i].Data()),
                                       Form("cluster from %s : E vs dispersion^{2}",ptypess[i].Data()),
                                       nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCEDispersion[i]->SetYTitle("D^{2}");
        fhMCEDispersion[i]->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCEDispersion[i]) ; 
                
        fhMCNCellsE[i]  = new TH2F (Form("hNCellsE_MC%s",pnamess[i].Data()),
                               Form("# of cells in cluster from %s vs E of clusters",ptypess[i].Data()), 
                                    nptbins,ptmin,ptmax, nbins,nmin,nmax); 
        fhMCNCellsE[i]->SetXTitle("E (GeV)");
        fhMCNCellsE[i]->SetYTitle("# of cells in cluster");
        outputContainer->Add(fhMCNCellsE[i]);  
        
        fhMCMaxCellDiffClusterE[i]  = new TH2F (Form("hMaxCellDiffClusterE_MC%s",pnamess[i].Data()),
                                             Form("energy vs difference of cluster energy from %s - max cell energy / cluster energy, good clusters",ptypess[i].Data()),
                                           nptbins,ptmin,ptmax, 500,0,1.); 
        fhMCMaxCellDiffClusterE[i]->SetXTitle("E_{cluster} (GeV) ");
        fhMCMaxCellDiffClusterE[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCMaxCellDiffClusterE[i]);  
        
        fhMCLambda0vsClusterMaxCellDiffE0[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE0_MC%s",pnamess[i].Data()),
                                    Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, E < 2 GeV",ptypess[i].Data()),
                                    ssbins,ssmin,ssmax,500,0,1.); 
        fhMCLambda0vsClusterMaxCellDiffE0[i]->SetXTitle("#lambda_{0}^{2}");
        fhMCLambda0vsClusterMaxCellDiffE0[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE0[i]) ; 
        
        fhMCLambda0vsClusterMaxCellDiffE2[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE2_MC%s",pnamess[i].Data()),
                                                         Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, 2< E < 6 GeV",ptypess[i].Data()),
                                                         ssbins,ssmin,ssmax,500,0,1.); 
        fhMCLambda0vsClusterMaxCellDiffE2[i]->SetXTitle("#lambda_{0}^{2}");
        fhMCLambda0vsClusterMaxCellDiffE2[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE2[i]) ; 
        
        fhMCLambda0vsClusterMaxCellDiffE6[i]  = new TH2F(Form("hLambda0vsClusterMaxCellDiffE6_MC%s",pnamess[i].Data()),
                                                         Form("cluster from %s : #lambda^{2}_{0} vs fraction of energy carried by max cell, E > 6 GeV",ptypess[i].Data()),
                                                         ssbins,ssmin,ssmax,500,0,1.); 
        fhMCLambda0vsClusterMaxCellDiffE6[i]->SetXTitle("#lambda_{0}^{2}");
        fhMCLambda0vsClusterMaxCellDiffE6[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCLambda0vsClusterMaxCellDiffE6[i]) ; 
        
        fhMCNCellsvsClusterMaxCellDiffE0[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE0_MC%s",pnamess[i].Data()),
                                                         Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, E < 2 GeV",ptypess[i].Data()),
                                                         nbins/5,nmin,nmax/5,500,0,1.); 
        fhMCNCellsvsClusterMaxCellDiffE0[i]->SetXTitle("N cells in cluster");
        fhMCNCellsvsClusterMaxCellDiffE0[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE0[i]) ; 
        
        fhMCNCellsvsClusterMaxCellDiffE2[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE2_MC%s",pnamess[i].Data()),
                                                         Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, 2< E < 6 GeV",ptypess[i].Data()),
                                                         nbins/5,nmin,nmax/5,500,0,1.); 
        fhMCNCellsvsClusterMaxCellDiffE2[i]->SetXTitle("N cells in cluster");
        fhMCNCellsvsClusterMaxCellDiffE2[i]->SetYTitle("(E_{cluster} - E_{cell max})/ E_{cluster}");
        outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE2[i]) ; 
        
        fhMCNCellsvsClusterMaxCellDiffE6[i]  = new TH2F(Form("hNCellsvsClusterMaxCellDiffE6_MC%s",pnamess[i].Data()),
                                                         Form("cluster from %s : N cells in cluster vs fraction of energy carried by max cell, E > 6 GeV",ptypess[i].Data()),
                                                         nbins/5,nmin,nmax/5,500,0,1.); 
        fhMCNCellsvsClusterMaxCellDiffE6[i]->SetXTitle("N cells in cluster");
        fhMCNCellsvsClusterMaxCellDiffE6[i]->SetYTitle("E (GeV)");
        outputContainer->Add(fhMCNCellsvsClusterMaxCellDiffE6[i]) ; 
        
       }// loop    
      
      if(!GetReader()->IsEmbeddedClusterSelectionOn())
      {
        fhMCPhotonELambda0NoOverlap  = new TH2F("hELambda0_MCPhoton_NoOverlap",
                                                "cluster from Photon : E vs #lambda_{0}^{2}",
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCPhotonELambda0NoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0NoOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCPhotonELambda0NoOverlap) ; 
        
        fhMCPhotonELambda0TwoOverlap  = new TH2F("hELambda0_MCPhoton_TwoOverlap",
                                                 "cluster from Photon : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCPhotonELambda0TwoOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0TwoOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCPhotonELambda0TwoOverlap) ; 
        
        fhMCPhotonELambda0NOverlap  = new TH2F("hELambda0_MCPhoton_NOverlap",
                                               "cluster from Photon : E vs #lambda_{0}^{2}",
                                               nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhMCPhotonELambda0NOverlap->SetYTitle("#lambda_{0}^{2}");
        fhMCPhotonELambda0NOverlap->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCPhotonELambda0NOverlap) ; 
        
      } //No embedding     
      
      //Fill histograms to check shape of embedded clusters
      if(GetReader()->IsEmbeddedClusterSelectionOn())
      {
        
        fhEmbeddedSignalFractionEnergy  = new TH2F("hEmbeddedSignal_FractionEnergy",
                                                    "Energy Fraction of embedded signal versus cluster energy",
                                                    nptbins,ptmin,ptmax,100,0.,1.); 
        fhEmbeddedSignalFractionEnergy->SetYTitle("Fraction");
        fhEmbeddedSignalFractionEnergy->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbeddedSignalFractionEnergy) ; 
        
        fhEmbedPhotonELambda0FullSignal  = new TH2F("hELambda0_EmbedPhoton_FullSignal",
                                                "cluster from Photon embedded with more than 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPhotonELambda0FullSignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0FullSignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0FullSignal) ; 
                
        fhEmbedPhotonELambda0MostlySignal  = new TH2F("hELambda0_EmbedPhoton_MostlySignal",
                                          "cluster from Photon embedded with 50% to 90% energy in cluster : E vs #lambda_{0}^{2}",
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPhotonELambda0MostlySignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0MostlySignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0MostlySignal) ; 
        
        fhEmbedPhotonELambda0MostlyBkg  = new TH2F("hELambda0_EmbedPhoton_MostlyBkg",
                                          "cluster from Photon embedded with 10% to 50% energy in cluster : E vs #lambda_{0}^{2}",
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPhotonELambda0MostlyBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0MostlyBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0MostlyBkg) ; 
                
        fhEmbedPhotonELambda0FullBkg  = new TH2F("hELambda0_EmbedPhoton_FullBkg",
                                                   "cluster from Photonm embedded with 0% to 10% energy in cluster : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPhotonELambda0FullBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPhotonELambda0FullBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPhotonELambda0FullBkg) ; 
        
        fhEmbedPi0ELambda0FullSignal  = new TH2F("hELambda0_EmbedPi0_FullSignal",
                                                    "cluster from Pi0 embedded with more than 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPi0ELambda0FullSignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0FullSignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0FullSignal) ; 
                
        fhEmbedPi0ELambda0MostlySignal  = new TH2F("hELambda0_EmbedPi0_MostlySignal",
                                                      "cluster from Pi0 embedded with 50% to 90% energy in cluster : E vs #lambda_{0}^{2}",
                                                      nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPi0ELambda0MostlySignal->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0MostlySignal->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0MostlySignal) ; 
        
        fhEmbedPi0ELambda0MostlyBkg  = new TH2F("hELambda0_EmbedPi0_MostlyBkg",
                                                   "cluster from Pi0 embedded with 10% to 50% energy in cluster : E vs #lambda_{0}^{2}",
                                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPi0ELambda0MostlyBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0MostlyBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0MostlyBkg) ; 
        
        fhEmbedPi0ELambda0FullBkg  = new TH2F("hELambda0_EmbedPi0_FullBkg",
                                                 "cluster from Pi0 embedded with 0% to 10% energy in cluster : E vs #lambda_{0}^{2}",
                                                 nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
        fhEmbedPi0ELambda0FullBkg->SetYTitle("#lambda_{0}^{2}");
        fhEmbedPi0ELambda0FullBkg->SetXTitle("E (GeV)");
        outputContainer->Add(fhEmbedPi0ELambda0FullBkg) ; 
  
      }// embedded histograms
      
      
    }// Fill SS MC histograms
    
  }//Histos with MC
      
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaPhoton::Init()
{
  
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaPhoton::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
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
  
  if(!pl) {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
  
  // Loop on raw clusters before filtering in the reader and fill control histogram
  if((GetReader()->GetEMCALClusterListName()=="" && fCalorimeter=="EMCAL") || fCalorimeter=="PHOS"){
    for(Int_t iclus = 0; iclus < GetReader()->GetInputEvent()->GetNumberOfCaloClusters(); iclus++ ){
      AliVCluster * clus = GetReader()->GetInputEvent()->GetCaloCluster(iclus);
      if     (fCalorimeter == "PHOS"  && clus->IsPHOS()  && clus->E() > GetReader()->GetPHOSPtMin() ) fhClusterCuts[0]->Fill(clus->E());
      else if(fCalorimeter == "EMCAL" && clus->IsEMCAL() && clus->E() > GetReader()->GetEMCALPtMin()) fhClusterCuts[0]->Fill(clus->E());
    }
  }
  else { // reclusterized
    TClonesArray * clusterList = 0;
    if(GetReader()->GetOutputEvent())
      clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetOutputEvent()->FindListObject(GetReader()->GetEMCALClusterListName()));
    if(clusterList){
      Int_t nclusters = clusterList->GetEntriesFast();
      for (Int_t iclus =  0; iclus <  nclusters; iclus++) {
        AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));        
        if(clus)fhClusterCuts[0]->Fill(clus->E());
      }  
    }
  }
  
  //Init arrays, variables, get number of clusters
  TLorentzVector mom, mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  
  //----------------------------------------------------
  // Fill AOD with PHOS/EMCAL AliAODPWG4Particle objects
  //----------------------------------------------------
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++){    
	  
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));	
    //printf("calo %d, %f\n",icalo,calo->E());
    
    //Get the index where the cluster comes, to retrieve the corresponding vertex
    Int_t evtIndex = 0 ; 
    if (GetMixedEvent()) {
      evtIndex=GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      //Get the vertex and check it is not too large in z
      if(TMath::Abs(GetVertex(evtIndex)[2])> GetZvertexCut()) continue;
    }
    
    //Cluster selection, not charged, with photon id and in fiducial cut	  
    if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
      calo->GetMomentum(mom,GetVertex(evtIndex)) ;}//Assume that come from vertex in straight line
    else{
      Double_t vertex[]={0,0,0};
      calo->GetMomentum(mom,vertex) ;
    }
    
    //--------------------------------------
    // Cluster selection
    //--------------------------------------
    if(!ClusterSelected(calo,mom)) continue;
    
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
    
    //--------------------------------------------------------------------------------------
    // Play with the MC stack if available
    //--------------------------------------------------------------------------------------
    
    //Check origin of the candidates
    Int_t tag = -1;
    
    if(IsDataMC()){
      
      tag = GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex());
      aodph.SetTag(tag);
      
      if(GetDebug() > 0)
        printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
        
    
    //--------------------------------------------------------------------------------------
    //Fill some shower shape histograms before PID is applied
    //--------------------------------------------------------------------------------------
    
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
    
    fhClusterCuts[8]->Fill(calo->E());
    
    // Matching after cuts
    if(fFillTMHisto) FillTrackMatchingResidualHistograms(calo,1);
    
    // Add number of local maxima to AOD, method name in AOD to be FIXED
    
    aodph.SetFiducialArea(GetCaloUtils()->GetNumberOfLocalMaxima(calo, cells));
    

    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  	
  if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD()  End fill AODs, with %d entries \n",GetOutputAODBranch()->GetEntriesFast());  
  
}

//__________________________________________________________________
void  AliAnaPhoton::MakeAnalysisFillHistograms() 
{
  //Fill histograms
  
  //-------------------------------------------------------------------
  // Access MC information in stack if requested, check that it exists.	
  AliStack         * stack       = 0x0;
  TParticle        * primary     = 0x0;   
  TClonesArray     * mcparticles = 0x0;
  AliAODMCParticle * aodprimary  = 0x0; 
  
  if(IsDataMC()){
    
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;
      if(!stack) {
        printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Stack not available, is the MC handler called? STOP\n");
        abort();
      }
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      
      //Get the list of MC particles
      mcparticles = GetReader()->GetAODMCParticles(0);
      if(!mcparticles && GetDebug() > 0) 	{
        printf("AliAnaPhoton::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
      }	
    }
  }// is data and MC
  
  
  // Get vertex
  Double_t v[3] = {0,0,0}; //vertex ;
  GetReader()->GetVertex(v);
  //fhVertex->Fill(v[0],v[1],v[2]);  
  if(TMath::Abs(v[2]) > GetZvertexCut()) return ; // done elsewhere for Single Event analysis, but there for mixed event
  
  //----------------------------------
  //Loop on stored AOD photons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaPhoton::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
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
    if     (ecluster > 0.5)   fhEtaPhiPhoton  ->Fill(etacluster, phicluster);
    else if(GetMinPt() < 0.5) fhEtaPhi05Photon->Fill(etacluster, phicluster);
    

    //Get original cluster, to recover some information
    Int_t absID             = 0; 
    Float_t maxCellFraction = 0;
    AliVCaloCells* cells    = 0;
    TObjArray * clusters    = 0; 
    if(fCalorimeter == "EMCAL"){
      cells    = GetEMCALCells();
      clusters = GetEMCALClusters();
    }
    else{
      cells    = GetPHOSCells();
      clusters = GetPHOSClusters();
    }
    
    Int_t iclus = -1;
    AliVCluster *cluster = FindCluster(clusters,ph->GetCaloLabel(0),iclus); 
    if(cluster)
    {
      absID = GetCaloUtils()->GetMaxEnergyCell(cells, cluster,maxCellFraction);
      
      // Control histograms
      fhMaxCellDiffClusterE->Fill(ph->E(),maxCellFraction);
      fhNCellsE            ->Fill(ph->E(),cluster->GetNCells());
      fhTimeE              ->Fill(ph->E(),cluster->GetTOF()*1.e9);
      if(cells)
      {
      for(Int_t icell = 0; icell <  cluster->GetNCells(); icell++)
        fhCellsE->Fill(ph->E(),cells->GetCellAmplitude(cluster->GetCellsAbsId()[icell]));
      }
    }
    
    //.......................................
    //Play with the MC data if available
    if(IsDataMC()){
      
      FillAcceptanceHistograms();
      
      //....................................................................
      // Access MC information in stack if requested, check that it exists.
      Int_t label =ph->GetLabel();
      if(label < 0) {
        if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** bad label ***:  label %d \n", label);
        continue;
      }
      
      Float_t eprim   = 0;
      Float_t ptprim  = 0;
      if(GetReader()->ReadStack()){
        
        if(label >=  stack->GetNtrack()) {
          if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
          continue ;
        }
        
        primary = stack->Particle(label);
        if(!primary){
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        eprim   = primary->Energy();
        ptprim  = primary->Pt();		
        
      }
      else if(GetReader()->ReadAODMCParticles()){
        //Check which is the input
        if(ph->GetInputFileIndex() == 0){
          if(!mcparticles) continue;
          if(label >=  mcparticles->GetEntriesFast()) {
            if(GetDebug() > 2)  printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** large label ***:  label %d, n tracks %d \n", 
                                       label, mcparticles->GetEntriesFast());
            continue ;
          }
          //Get the particle
          aodprimary = (AliAODMCParticle*) mcparticles->At(label);
          
        }
        
        if(!aodprimary){
          printf("AliAnaPhoton::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", label);
          continue;
        }
        
        eprim   = aodprimary->E();
        ptprim  = aodprimary->Pt();
        
      }
      
      Int_t tag =ph->GetTag();
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[kmcPhoton])
      {
        fhMCE  [kmcPhoton] ->Fill(ecluster);
        fhMCPt [kmcPhoton] ->Fill(ptcluster);
        fhMCPhi[kmcPhoton] ->Fill(ecluster,phicluster);
        fhMCEta[kmcPhoton] ->Fill(ecluster,etacluster);
        
        fhMC2E[kmcPhoton]     ->Fill(ecluster, eprim);
        fhMC2Pt[kmcPhoton]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[kmcPhoton] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcPhoton]->Fill(ptcluster,ptprim-ptcluster);     
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[kmcConversion])
        {
          fhMCE  [kmcConversion] ->Fill(ecluster);
          fhMCPt [kmcConversion] ->Fill(ptcluster);
          fhMCPhi[kmcConversion] ->Fill(ecluster,phicluster);
          fhMCEta[kmcConversion] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcConversion]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcConversion]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcConversion] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcConversion]->Fill(ptcluster,ptprim-ptcluster);     
        }			
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhMCE[kmcPrompt]){
          fhMCE  [kmcPrompt] ->Fill(ecluster);
          fhMCPt [kmcPrompt] ->Fill(ptcluster);
          fhMCPhi[kmcPrompt] ->Fill(ecluster,phicluster);
          fhMCEta[kmcPrompt] ->Fill(ecluster,etacluster);      
          
          fhMC2E[kmcPrompt]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcPrompt]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcPrompt] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcPrompt]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)&& fhMCE[kmcFragmentation])
        {
          fhMCE  [kmcFragmentation] ->Fill(ecluster);
          fhMCPt [kmcFragmentation] ->Fill(ptcluster);
          fhMCPhi[kmcFragmentation] ->Fill(ecluster,phicluster);
          fhMCEta[kmcFragmentation] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcFragmentation]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcFragmentation]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcFragmentation] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcFragmentation]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR)&& fhMCE[kmcISR])
        {
          fhMCE  [kmcISR] ->Fill(ecluster);
          fhMCPt [kmcISR] ->Fill(ptcluster);
          fhMCPhi[kmcISR] ->Fill(ecluster,phicluster);
          fhMCEta[kmcISR] ->Fill(ecluster,etacluster);    
          
          fhMC2E[kmcISR]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcISR]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcISR] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcISR]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[kmcPi0Decay])
        {
          fhMCE  [kmcPi0Decay] ->Fill(ecluster);
          fhMCPt [kmcPi0Decay] ->Fill(ptcluster);
          fhMCPhi[kmcPi0Decay] ->Fill(ecluster,phicluster);
          fhMCEta[kmcPi0Decay] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcPi0Decay]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcPi0Decay]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcPi0Decay] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcPi0Decay]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                  GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[kmcOtherDecay])
        {
          fhMCE  [kmcOtherDecay] ->Fill(ecluster);
          fhMCPt [kmcOtherDecay] ->Fill(ptcluster);
          fhMCPhi[kmcOtherDecay] ->Fill(ecluster,phicluster);
          fhMCEta[kmcOtherDecay] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcOtherDecay]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcOtherDecay]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcOtherDecay] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcOtherDecay]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [kmcPi0])
        {
          fhMCE  [kmcPi0] ->Fill(ecluster);
          fhMCPt [kmcPi0] ->Fill(ptcluster);
          fhMCPhi[kmcPi0] ->Fill(ecluster,phicluster);
          fhMCEta[kmcPi0] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcPi0]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcPi0]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcPi0] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcPi0]->Fill(ptcluster,ptprim-ptcluster);     
          
        } 
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[kmcEta])
        {
          fhMCE  [kmcEta] ->Fill(ecluster);
          fhMCPt [kmcEta] ->Fill(ptcluster);
          fhMCPhi[kmcEta] ->Fill(ecluster,phicluster);
          fhMCEta[kmcEta] ->Fill(ecluster,etacluster);
          
          fhMC2E[kmcEta]     ->Fill(ecluster, eprim);
          fhMC2Pt[kmcEta]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[kmcEta] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[kmcEta]->Fill(ptcluster,ptprim-ptcluster);     
          
        }      
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[kmcAntiNeutron])
      {
        fhMCE  [kmcAntiNeutron] ->Fill(ecluster);
        fhMCPt [kmcAntiNeutron] ->Fill(ptcluster);
        fhMCPhi[kmcAntiNeutron] ->Fill(ecluster,phicluster);
        fhMCEta[kmcAntiNeutron] ->Fill(ecluster,etacluster);
        
        fhMC2E[kmcAntiNeutron]     ->Fill(ecluster, eprim);
        fhMC2Pt[kmcAntiNeutron]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[kmcAntiNeutron] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcAntiNeutron]->Fill(ptcluster,ptprim-ptcluster);     
        
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[kmcAntiProton])
      {
        fhMCE  [kmcAntiProton] ->Fill(ecluster);
        fhMCPt [kmcAntiProton] ->Fill(ptcluster);
        fhMCPhi[kmcAntiProton] ->Fill(ecluster,phicluster);
        fhMCEta[kmcAntiProton] ->Fill(ecluster,etacluster);
        
        fhMC2E[kmcAntiProton]     ->Fill(ecluster, eprim);
        fhMC2Pt[kmcAntiProton]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[kmcAntiProton] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcAntiProton]->Fill(ecluster,ptprim-ptcluster);     
        
      } 
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[kmcElectron])
      {
        fhMCE  [kmcElectron] ->Fill(ecluster);
        fhMCPt [kmcElectron] ->Fill(ptcluster);
        fhMCPhi[kmcElectron] ->Fill(ecluster,phicluster);
        fhMCEta[kmcElectron] ->Fill(ecluster,etacluster);
        
        fhMC2E[kmcElectron]     ->Fill(ecluster, eprim);
        fhMC2Pt[kmcElectron]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[kmcElectron] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcElectron]->Fill(ecluster,ptprim-ptcluster);             
      }     
      else if( fhMCE[kmcOther]){
        fhMCE  [kmcOther] ->Fill(ecluster);
        fhMCPt [kmcOther] ->Fill(ptcluster);
        fhMCPhi[kmcOther] ->Fill(ecluster,phicluster);
        fhMCEta[kmcOther] ->Fill(ecluster,etacluster);
        
        fhMC2E[kmcOther]     ->Fill(ecluster, eprim);
        fhMC2Pt[kmcOther]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[kmcOther] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[kmcOther]->Fill(ecluster,ptprim-ptcluster);     
        
        //		 printf(" AliAnaPhoton::MakeAnalysisFillHistograms() - Label %d, pT %2.3f Unknown, bits set: ",
        //					ph->GetLabel(),ph->Pt());
        //		  for(Int_t i = 0; i < 20; i++) {
        //			  if(GetMCAnalysisUtils()->CheckTagBit(tag,i)) printf(" %d, ",i);
        //		  }
        //		  printf("\n");
        
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
