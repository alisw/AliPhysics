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
/* $Id: AliAnaPhoton.cxx 28688 2008-09-11 15:04:07Z gconesab $ */

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
//#include <Riostream.h>
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


ClassImp(AliAnaPhoton)
  
//____________________________________________________________________________
AliAnaPhoton::AliAnaPhoton() : 
    AliAnaPartCorrBaseClass(),    fCalorimeter(""), 
    fMinDist(0.),                 fMinDist2(0.),                fMinDist3(0.), 
    fRejectTrackMatch(0),         fTimeCutMin(-1),              fTimeCutMax(999999),         
    fNCellsCut(0),                fFillSSHistograms(kFALSE),    
    fNOriginHistograms(8),        fNPrimaryHistograms(4),
    fCheckConversion(kFALSE),     fRemoveConvertedPair(kFALSE), 
    fAddConvertedPairsToAOD(kFALSE), 
    fMassCut(0),                  fConvAsymCut(1.),             fConvDEtaCut(2.),
    fConvDPhiMinCut(-1.),         fConvDPhiMaxCut(7.), 

    // Histograms
    fhNCellsE(0),                 fhMaxCellDiffClusterE(0),
    fhEPhoton(0),                 fhPtPhoton(0),  
    fhPhiPhoton(0),               fhEtaPhoton(0), 
    fhEtaPhiPhoton(0),            fhEtaPhi05Photon(0),

    // Conversion histograms
    fhPtPhotonConv(0),            fhEtaPhiPhotonConv(0),        fhEtaPhi05PhotonConv(0),
    fhConvDeltaEta(0),            fhConvDeltaPhi(0),            fhConvDeltaEtaPhi(0), 
    fhConvAsym(0),                fhConvPt(0),
    fhConvDistEta(0),             fhConvDistEn(0),              fhConvDistMass(0),     
    fhConvDistEtaCutEta(0),       fhConvDistEnCutEta(0),        fhConvDistMassCutEta(0),
    fhConvDistEtaCutMass(0),      fhConvDistEnCutMass(0), 
    fhConvDistEtaCutAsy(0),       fhConvDistEnCutAsy(0),

    // Shower shape histograms
    fhDispE(0),                   fhLam0E(0),                   fhLam1E(0), 
    fhDispETRD(0),                fhLam0ETRD(0),                fhLam1ETRD(0),

    fhNCellsLam0LowE(0),          fhNCellsLam1LowE(0),          fhNCellsDispLowE(0),  
    fhNCellsLam0HighE(0),         fhNCellsLam1HighE(0),         fhNCellsDispHighE(0),

    fhEtaLam0LowE(0),             fhPhiLam0LowE(0), 
    fhEtaLam0HighE(0),            fhPhiLam0HighE(0), 
    fhLam0DispLowE(0),            fhLam0DispHighE(0), 
    fhLam1Lam0LowE(0),            fhLam1Lam0HighE(0), 
    fhDispLam1LowE(0),            fhDispLam1HighE(0),

    // MC histograms
    // Conversion MC histograms
    fhPtConversionTagged(0),           fhPtAntiNeutronTagged(0),       
    fhPtAntiProtonTagged(0),           fhPtUnknownTagged(0),
    fhEtaPhiConversion(0),             fhEtaPhi05Conversion(0),

    fhConvDeltaEtaMCConversion(0),     fhConvDeltaPhiMCConversion(0),  fhConvDeltaEtaPhiMCConversion(0),
    fhConvAsymMCConversion(0),         fhConvPtMCConversion(0),           
    fhConvDispersionMCConversion(0),   fhConvM02MCConversion(0),

    fhConvDeltaEtaMCAntiNeutron(0),    fhConvDeltaPhiMCAntiNeutron(0), fhConvDeltaEtaPhiMCAntiNeutron(0), 
    fhConvAsymMCAntiNeutron(0),        fhConvPtMCAntiNeutron(0), 
    fhConvDispersionMCAntiNeutron(0),  fhConvM02MCAntiNeutron(0),
    fhConvDeltaEtaMCAntiProton(0),     fhConvDeltaPhiMCAntiProton(0),  fhConvDeltaEtaPhiMCAntiProton(0),  
    fhConvAsymMCAntiProton(0),         fhConvPtMCAntiProton(0),  
    fhConvDispersionMCAntiProton(0),   fhConvM02MCAntiProton(0),
    fhConvDeltaEtaMCString(0),         fhConvDeltaPhiMCString(0),      fhConvDeltaEtaPhiMCString(0),      
    fhConvAsymMCString(0),             fhConvPtMCString(0),      
    fhConvDispersionMCString(0),       fhConvM02MCString(0),
    fhConvDistMCConversion(0),         fhConvDistMCConversionCuts(0),

    // Photon SS MC histograms
    fhMCPhotonELambda0NoOverlap(0),    fhMCPhotonELambda0TwoOverlap(0),  fhMCPhotonELambda0NOverlap(0),

    //Embedding
    fhEmbeddedSignalFractionEnergy(0),     
    fhEmbedPhotonELambda0FullSignal(0), fhEmbedPhotonELambda0MostlySignal(0),  
    fhEmbedPhotonELambda0MostlyBkg(0),  fhEmbedPhotonELambda0FullBkg(0),       
    fhEmbedPi0ELambda0FullSignal(0),    fhEmbedPi0ELambda0MostlySignal(0),    
    fhEmbedPi0ELambda0MostlyBkg(0),     fhEmbedPi0ELambda0FullBkg(0)         
{
  //default ctor
  
  for(Int_t i = 0; i < 14; i++){
    fhMCPt     [i] = 0;
    fhMCE      [i] = 0;
    fhMCPhi    [i] = 0;
    fhMCEta    [i] = 0;
    fhMCDeltaE [i] = 0;                
    fhMCDeltaPt[i] = 0;
    fhMC2E     [i] = 0;              
    fhMC2Pt    [i] = 0;
    
  }
  
  for(Int_t i = 0; i < 7; i++){
    fhPtPrimMC [i] = 0;
    fhEPrimMC  [i] = 0;
    fhPhiPrimMC[i] = 0;
    fhYPrimMC  [i] = 0;
    
    fhPtPrimMCAcc [i] = 0;
    fhEPrimMCAcc  [i] = 0;
    fhPhiPrimMCAcc[i] = 0;
    fhYPrimMCAcc  [i] = 0;
  }  
  
  for(Int_t i = 0; i < 6; i++){
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
  
  //Initialize parameters
  InitParameters();

}//____________________________________________________________________________
AliAnaPhoton::~AliAnaPhoton() 
{
  //dtor

}

//__________________________________________________________________
Bool_t  AliAnaPhoton::ClusterSelected(AliVCluster* calo, TLorentzVector mom) 
{
  //Select clusters if they pass different cuts
  if(GetDebug() > 2) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; Before selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(),
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
  //.......................................
  //If too small or big energy, skip it
  if(mom.E() < GetMinEnergy() || mom.E() > GetMaxEnergy() ) return kFALSE ; 
  if(GetDebug() > 2) printf("\t Cluster %d Pass E Cut \n",calo->GetID());
  
  //.......................................
  // TOF cut, BE CAREFUL WITH THIS CUT
  Double_t tof = calo->GetTOF()*1e9;
  if(tof < fTimeCutMin || tof > fTimeCutMax) return kFALSE;
  if(GetDebug() > 2)  printf("\t Cluster %d Pass Time Cut \n",calo->GetID());
  
  //.......................................
  if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) return kFALSE;
  if(GetDebug() > 2) printf("\t Cluster %d Pass NCell Cut \n",calo->GetID());
  
  //.......................................
  //Check acceptance selection
  if(IsFiducialCutOn()){
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
    if(! in ) return kFALSE ;
  }
  if(GetDebug() > 2) printf("Fiducial cut passed \n");
  
  //.......................................
  //Skip matched clusters with tracks
  if(fRejectTrackMatch){
    if(IsTrackMatched(calo)) {
      if(GetDebug() > 2) printf("\t Reject track-matched clusters\n");
      return kFALSE ;
    }
    else  
      if(GetDebug() > 2)  printf(" Track-matching cut passed \n");
  }// reject matched clusters
  
  //.......................................
  //Check Distance to Bad channel, set bit.
  Double_t distBad=calo->GetDistanceToBadChannel() ; //Distance to bad channel
  if(distBad < 0.) distBad=9999. ; //workout strange convension dist = -1. ;
  if(distBad < fMinDist) {//In bad channel (PHOS cristal size 2.2x2.2 cm), EMCAL ( cell units )
    return kFALSE ;
  }
  else if(GetDebug() > 2) printf("\t Bad channel cut passed %4.2f > %2.2f \n",distBad, fMinDist);
  //printf("Cluster %d Pass Bad Dist Cut \n",icalo);

  if(GetDebug() > 0) 
    printf("AliAnaPhoton::ClusterSelected() Current Event %d; After  selection : E %2.2f, pT %2.2f, Ecl %2.2f, phi %2.2f, eta %2.2f\n",
           GetReader()->GetEventNumber(), 
           mom.E(), mom.Pt(),calo->E(),mom.Phi()*TMath::RadToDeg(),mom.Eta());
  
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
          
          fhYPrimMC[mcPPhoton]->Fill(photonPt, photonY) ;
          if(TMath::Abs(photonY) < 1.0)
          {
            fhEPrimMC  [mcPPhoton]->Fill(photonE ) ;
            fhPtPrimMC [mcPPhoton]->Fill(photonPt) ;
            fhPhiPrimMC[mcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMC[mcPPhoton]  ->Fill(photonE , photonEta) ;
          }
          if(inacceptance){
            fhEPrimMCAcc[mcPPhoton]  ->Fill(photonE ) ;
            fhPtPrimMCAcc[mcPPhoton] ->Fill(photonPt) ;
            fhPhiPrimMCAcc[mcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMCAcc[mcPPhoton]  ->Fill(photonE , photonY) ;
          }//Accepted
          
          //Origin of photon
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhEPrimMC[mcPPrompt])
          {
            fhYPrimMC[mcPPrompt]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPPrompt]->Fill(photonE ) ;
              fhPtPrimMC [mcPPrompt]->Fill(photonPt) ;
              fhPhiPrimMC[mcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPPrompt]  ->Fill(photonE , photonEta) ;
            }   
            if(inacceptance){
              fhEPrimMCAcc[mcPPrompt]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPPrompt] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPPrompt]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) && fhEPrimMC[mcPFragmentation])
          {
            fhYPrimMC[mcPFragmentation]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPFragmentation]->Fill(photonE ) ;
              fhPtPrimMC [mcPFragmentation]->Fill(photonPt) ;
              fhPhiPrimMC[mcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPFragmentation]  ->Fill(photonE , photonEta) ;
            }  
            if(inacceptance){
              fhEPrimMCAcc[mcPFragmentation]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPFragmentation] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPFragmentation]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) && fhEPrimMC[mcPISR])
          {
            fhYPrimMC[mcPISR]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPISR]->Fill(photonE ) ;
              fhPtPrimMC [mcPISR]->Fill(photonPt) ;
              fhPhiPrimMC[mcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPISR]->Fill(photonE , photonEta) ;
            }            
            if(inacceptance){
              fhEPrimMCAcc[mcPISR]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPISR] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPISR]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)&& fhEPrimMC[mcPPi0Decay])
          {
            fhYPrimMC[mcPPi0Decay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPPi0Decay]->Fill(photonE ) ;
              fhPtPrimMC [mcPPi0Decay]->Fill(photonPt) ;
              fhPhiPrimMC[mcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPPi0Decay]  ->Fill(photonE , photonEta) ;
            }     
            if(inacceptance){
              fhEPrimMCAcc[mcPPi0Decay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPPi0Decay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPPi0Decay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                    GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay)) && fhEPrimMC[mcPOtherDecay])
          {
            fhYPrimMC[mcPOtherDecay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPOtherDecay]->Fill(photonE ) ;
              fhPtPrimMC [mcPOtherDecay]->Fill(photonPt) ;
              fhPhiPrimMC[mcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPOtherDecay]  ->Fill(photonE , photonEta) ;
            } 
            if(inacceptance){
              fhEPrimMCAcc[mcPOtherDecay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPOtherDecay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPOtherDecay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(fhEPrimMC[mcPOther])
          {
            fhYPrimMC[mcPOther]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPOther]->Fill(photonE ) ;
              fhPtPrimMC [mcPOther]->Fill(photonPt) ;
              fhPhiPrimMC[mcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPOther]  ->Fill(photonE , photonEta) ;
            }
            if(inacceptance){
              fhEPrimMCAcc[mcPOther]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPOther] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPOther]  ->Fill(photonE , photonY) ;
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
          
          fhYPrimMC[mcPPhoton]->Fill(photonPt, photonY) ;
          if(TMath::Abs(photonY) < 1.0)
          {
            fhEPrimMC  [mcPPhoton]->Fill(photonE ) ;
            fhPtPrimMC [mcPPhoton]->Fill(photonPt) ;
            fhPhiPrimMC[mcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMC[mcPPhoton]->Fill(photonE , photonEta) ;
          }
          if(inacceptance){
            fhEPrimMCAcc[mcPPhoton]  ->Fill(photonE ) ;
            fhPtPrimMCAcc[mcPPhoton] ->Fill(photonPt) ;
            fhPhiPrimMCAcc[mcPPhoton]->Fill(photonE , photonPhi) ;
            fhYPrimMCAcc[mcPPhoton]  ->Fill(photonE , photonY) ;
          }//Accepted
          
          
          //Origin of photon
          if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhEPrimMC[mcPPrompt])
          {
            fhYPrimMC[mcPPrompt]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPPrompt]->Fill(photonE ) ;
              fhPtPrimMC [mcPPrompt]->Fill(photonPt) ;
              fhPhiPrimMC[mcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPPrompt]->Fill(photonE , photonEta) ;
            }   
            if(inacceptance){
              fhEPrimMCAcc[mcPPrompt]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPPrompt] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPPrompt]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPPrompt]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation) && fhEPrimMC[mcPFragmentation] )
          {
            fhYPrimMC[mcPFragmentation]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPFragmentation]->Fill(photonE ) ;
              fhPtPrimMC [mcPFragmentation]->Fill(photonPt) ;
              fhPhiPrimMC[mcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPFragmentation]->Fill(photonE , photonEta) ;
            }  
            if(inacceptance){
              fhEPrimMCAcc[mcPFragmentation]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPFragmentation] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPFragmentation]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPFragmentation]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR) && fhEPrimMC[mcPISR])
          {
            fhYPrimMC[mcPISR]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPISR]->Fill(photonE ) ;
              fhPtPrimMC [mcPISR]->Fill(photonPt) ;
              fhPhiPrimMC[mcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPISR]->Fill(photonE , photonEta) ;
            }            
            if(inacceptance){
              fhEPrimMCAcc[mcPISR]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPISR] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPISR]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPISR]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay)&& fhEPrimMC[mcPPi0Decay])
          {
            fhYPrimMC[mcPPi0Decay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPPi0Decay]->Fill(photonE ) ;
              fhPtPrimMC [mcPPi0Decay]->Fill(photonPt) ;
              fhPhiPrimMC[mcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPPi0Decay]->Fill(photonE , photonEta) ;
            }     
            if(inacceptance){
              fhEPrimMCAcc[mcPPi0Decay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPPi0Decay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPPi0Decay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPPi0Decay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if((GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                   GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhEPrimMC[mcPOtherDecay])
          {
            fhYPrimMC[mcPOtherDecay]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPOtherDecay]->Fill(photonE ) ;
              fhPtPrimMC [mcPOtherDecay]->Fill(photonPt) ;
              fhPhiPrimMC[mcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPOtherDecay]->Fill(photonE , photonEta) ;
            } 
            if(inacceptance){
              fhEPrimMCAcc[mcPOtherDecay]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPOtherDecay] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPOtherDecay]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPOtherDecay]  ->Fill(photonE , photonY) ;
            }//Accepted
          }
          else if(fhEPrimMC[mcPOther])
          {
            fhYPrimMC[mcPOther]->Fill(photonPt, photonY) ;
            if(TMath::Abs(photonY) < 1.0){
              fhEPrimMC  [mcPOther]->Fill(photonE ) ;
              fhPtPrimMC [mcPOther]->Fill(photonPt) ;
              fhPhiPrimMC[mcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMC[mcPOther]->Fill(photonE , photonEta) ;
            }
            if(inacceptance){
              fhEPrimMCAcc[mcPOther]  ->Fill(photonE ) ;
              fhPtPrimMCAcc[mcPOther] ->Fill(photonPt) ;
              fhPhiPrimMCAcc[mcPOther]->Fill(photonE , photonPhi) ;
              fhYPrimMCAcc[mcPOther]  ->Fill(photonE , photonY) ;
            }//Accepted
          }//Other origin
        }// Primary photon 
      }//loop on primaries	
      
    }//mc array exists and data is MC
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
      
      if(GetDebug() > 1 ) printf("AliAnaPhoton::FillShowerShapeHistogram() - Energy fraction of embedded signal %2.3f, Energy %2.3f\n",fraction, clusterE);
      
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
      fhMCELambda0[mcssPhoton]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssPhoton]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssPhoton] ->Fill(energy, disp);
      fhMCNCellsE[mcssPhoton]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssPhoton]->Fill(energy,maxCellFraction);
      
      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssPhoton]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssPhoton] ->Fill(ncells,  maxCellFraction);
      }
      
      if(!GetReader()->IsEmbeddedClusterSelectionOn()){
        //Check particle overlaps in cluster
        
        //Compare the primary depositing more energy with the rest, if no photon/electron as comon ancestor (conversions), count as other particle
        Int_t ancPDG = 0, ancStatus = -1;
        TLorentzVector momentum; TVector3 prodVertex;
        Int_t ancLabel = 0;
        Int_t noverlaps = 1;      
        for (UInt_t ilab = 0; ilab < cluster->GetNLabels(); ilab++ ) {
          ancLabel = GetMCAnalysisUtils()->CheckCommonAncestor(cluster->GetLabels()[0],cluster->GetLabels()[ilab], GetReader(),ancPDG,ancStatus,momentum,prodVertex);
          if(ancPDG!=22 && TMath::Abs(ancPDG)!=11) noverlaps++;
        }
        
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
      fhMCELambda0[mcssElectron]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssElectron]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssElectron] ->Fill(energy, disp);
      fhMCNCellsE[mcssElectron]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssElectron]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssElectron] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssElectron] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssElectron]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssElectron] ->Fill(ncells,  maxCellFraction);
      }
    }//electron
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPhoton) && 
              GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCConversion) ){
      fhMCELambda0[mcssConversion]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssConversion]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssConversion] ->Fill(energy, disp);
      fhMCNCellsE[mcssConversion]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssConversion]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssConversion] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssConversion] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssConversion]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssConversion] ->Fill(ncells,  maxCellFraction);
      }      
      
    }//conversion photon
    else if  ( GetMCAnalysisUtils()->CheckTagBit(mcTag,AliMCAnalysisUtils::kMCPi0)  ){
      fhMCELambda0[mcssPi0]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssPi0]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssPi0] ->Fill(energy, disp);
      fhMCNCellsE[mcssPi0]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssPi0]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssPi0] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssPi0] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssPi0]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssPi0] ->Fill(ncells,  maxCellFraction);
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
      fhMCELambda0[mcssEta]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssEta]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssEta] ->Fill(energy, disp);
      fhMCNCellsE[mcssEta]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssEta]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssEta] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssEta] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssEta]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssEta] ->Fill(ncells,  maxCellFraction);
      }
      
    }//eta    
    else {
      fhMCELambda0[mcssOther]    ->Fill(energy, lambda0);
      fhMCELambda1[mcssOther]    ->Fill(energy, lambda1);
      fhMCEDispersion[mcssOther] ->Fill(energy, disp);
      fhMCNCellsE[mcssOther]     ->Fill(energy, ncells);
      fhMCMaxCellDiffClusterE[mcssOther]->Fill(energy,maxCellFraction);

      if     (energy < 2.){
        fhMCLambda0vsClusterMaxCellDiffE0[mcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE0[mcssOther] ->Fill(ncells,  maxCellFraction);
      }
      else if(energy < 6.){
        fhMCLambda0vsClusterMaxCellDiffE2[mcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE2[mcssOther] ->Fill(ncells,  maxCellFraction);
      }
      else{
        fhMCLambda0vsClusterMaxCellDiffE6[mcssOther]->Fill(lambda0, maxCellFraction);
        fhMCNCellsvsClusterMaxCellDiffE6[mcssOther] ->Fill(ncells,  maxCellFraction);
      }            
      
    }//other particles 
    
  }//MC data
  
}

//________________________________________________________________________
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
  snprintf(onePar,buffersize,"Conversion Selection: fConvAsymCut %1.2f, fConvDEtaCut %1.2f fConvDPhiCut (%1.2f,%1.2f)\n",
           fConvAsymCut, fConvDEtaCut, fConvDPhiMinCut, fConvDPhiMaxCut) ;
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
	
  Int_t nptbins  = GetHistoPtBins();  Float_t ptmax  = GetHistoPtMax();  Float_t ptmin  = GetHistoPtMin(); 
  Int_t nphibins = GetHistoPhiBins(); Float_t phimax = GetHistoPhiMax(); Float_t phimin = GetHistoPhiMin(); 
  Int_t netabins = GetHistoEtaBins(); Float_t etamax = GetHistoEtaMax(); Float_t etamin = GetHistoEtaMin();	
  Int_t ssbins   = GetHistoShowerShapeBins();  Float_t ssmax = GetHistoShowerShapeMax();  Float_t ssmin = GetHistoShowerShapeMin();
  Int_t nbins    = GetHistoNClusterCellBins(); Int_t nmax    = GetHistoNClusterCellMax(); Int_t nmin    = GetHistoNClusterCellMin(); 

  fhNCellsE  = new TH2F ("hNCellsE","# of cells in cluster vs E of clusters", nptbins,ptmin,ptmax, nbins,nmin,nmax); 
  fhNCellsE->SetXTitle("E (GeV)");
  fhNCellsE->SetYTitle("# of cells in cluster");
  outputContainer->Add(fhNCellsE);  
  
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
  
  //Conversion
  if(fCheckConversion){
    
    fhEtaPhiConversion  = new TH2F
    ("hEtaPhiConversion","#eta vs #phi for converted clusters",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiConversion->SetYTitle("#phi (rad)");
    fhEtaPhiConversion->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiConversion) ;
    if(GetMinPt() < 0.5){
      fhEtaPhi05Conversion  = new TH2F
      ("hEtaPhi05Conversion","#eta vs #phi, E > 0.5, for converted clusters",netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05Conversion->SetYTitle("#phi (rad)");
      fhEtaPhi05Conversion->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05Conversion) ;
    }    
    
    fhPtPhotonConv  = new TH1F("hPtPhotonConv","Number of #gamma over calorimeter, conversion",nptbins,ptmin,ptmax); 
    fhPtPhotonConv->SetYTitle("N");
    fhPtPhotonConv->SetXTitle("p_{T #gamma}(GeV/c)");
    outputContainer->Add(fhPtPhotonConv) ; 
    
    fhEtaPhiPhotonConv  = new TH2F
    ("hEtaPhiPhotonConv","#eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhiPhotonConv->SetYTitle("#phi (rad)");
    fhEtaPhiPhotonConv->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiPhotonConv) ;
    if(GetMinPt() < 0.5){
      fhEtaPhi05PhotonConv  = new TH2F
      ("hEtaPhi05PhotonConv","#eta vs #phi, E > 0.5",netabins,etamin,etamax,nphibins,phimin,phimax); 
      fhEtaPhi05PhotonConv->SetYTitle("#phi (rad)");
      fhEtaPhi05PhotonConv->SetXTitle("#eta");
      outputContainer->Add(fhEtaPhi05PhotonConv) ;
    }
    
    fhConvDeltaEta  = new TH2F
    ("hConvDeltaEta","#Delta #eta of selected conversion pairs",100,0,fMassCut,netabins*2,-0.5,0.5); 
    fhConvDeltaEta->SetYTitle("#Delta #eta");
    fhConvDeltaEta->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaEta) ;
    
    fhConvDeltaPhi  = new TH2F
    ("hConvDeltaPhi","#Delta #phi of selected conversion pairs",100,0,fMassCut,nphibins*2,-0.5,0.5); 
    fhConvDeltaPhi->SetYTitle("#Delta #phi");
    fhConvDeltaPhi->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvDeltaPhi) ;
    
    fhConvDeltaEtaPhi  = new TH2F
    ("hConvDeltaEtaPhi","#Delta #eta vs #Delta #phi of selected conversion pairs",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
    fhConvDeltaEtaPhi->SetYTitle("#Delta #phi");
    fhConvDeltaEtaPhi->SetXTitle("#Delta #eta");
    outputContainer->Add(fhConvDeltaEtaPhi) ;
    
    fhConvAsym  = new TH2F
    ("hConvAsym","Asymmetry of selected conversion pairs",100,0,fMassCut,100,0,1); 
    fhConvAsym->SetYTitle("Asymmetry");
    fhConvAsym->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvAsym) ;  
    
    fhConvPt  = new TH2F
    ("hConvPt","p_{T} of selected conversion pairs",100,0,fMassCut,100,0.,10.); 
    fhConvPt->SetYTitle("Pair p_{T} (GeV/c)");
    fhConvPt->SetXTitle("Pair Mass (GeV/c^2)");
    outputContainer->Add(fhConvPt) ;
    
    fhConvDistEta  = new TH2F
    ("hConvDistEta","distance to conversion vertex",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEta->SetXTitle("#eta");
    fhConvDistEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEta) ;
    
    fhConvDistEn  = new TH2F
    ("hConvDistEn","distance to conversion vertex",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEn->SetXTitle("E (GeV)");
    fhConvDistEn->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEn) ;

    fhConvDistMass  = new TH2F
    ("hConvDistMass","distance to conversion vertex",100,0,fMassCut,100,0.,5.); 
    fhConvDistMass->SetXTitle("m (GeV/c^2)");
    fhConvDistMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistMass) ;
    
    fhConvDistEtaCutEta  = new TH2F
    ("hConvDistEtaCutEta","distance to conversion vertex, dEta < 0.05",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutEta->SetXTitle("#eta");
    fhConvDistEtaCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutEta) ;
    
    fhConvDistEnCutEta  = new TH2F
    ("hConvDistEnCutEta","distance to conversion vertex, dEta < 0.05",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutEta->SetXTitle("E (GeV)");
    fhConvDistEnCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutEta) ;
    
    fhConvDistMassCutEta  = new TH2F
    ("hConvDistMassCutEta","distance to conversion vertex, dEta < 0.05",100,0,fMassCut,100,0.,5.); 
    fhConvDistMassCutEta->SetXTitle("m (GeV/c^2)");
    fhConvDistMassCutEta->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistMassCutEta) ;
    
    fhConvDistEtaCutMass  = new TH2F
    ("hConvDistEtaCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutMass->SetXTitle("#eta");
    fhConvDistEtaCutMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutMass) ;
    
    fhConvDistEnCutMass  = new TH2F
    ("hConvDistEnCutMass","distance to conversion vertex, dEta < 0.05, m < 10 MeV",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutMass->SetXTitle("E (GeV)");
    fhConvDistEnCutMass->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutMass) ;

    fhConvDistEtaCutAsy  = new TH2F
    ("hConvDistEtaCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",100,-0.7,0.7,100,0.,5.); 
    fhConvDistEtaCutAsy->SetXTitle("#eta");
    fhConvDistEtaCutAsy->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEtaCutAsy) ;
    
    fhConvDistEnCutAsy  = new TH2F
    ("hConvDistEnCutAsy","distance to conversion vertex, dEta < 0.05, m < 10 MeV, A < 0.1",nptbins,ptmin,ptmax,100,0.,5.); 
    fhConvDistEnCutAsy->SetXTitle("E (GeV)");
    fhConvDistEnCutAsy->SetYTitle(" distance (m)");
    outputContainer->Add(fhConvDistEnCutAsy) ;
    
  } // check conversion
  
  
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
    
    if(fCalorimeter == "EMCAL"){
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
        	
    if(fCheckConversion){  
      fhPtConversionTagged  = new TH1F("hPtMCConversionTagged","Number of converted #gamma over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtConversionTagged->SetYTitle("N");
      fhPtConversionTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtConversionTagged) ; 
      
      
      fhPtAntiNeutronTagged  = new TH1F("hPtMCAntiNeutronTagged","Number of AntiNeutron id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtAntiNeutronTagged->SetYTitle("N");
      fhPtAntiNeutronTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtAntiNeutronTagged) ; 
      
      fhPtAntiProtonTagged  = new TH1F("hPtMCAntiProtonTagged","Number of AntiProton id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtAntiProtonTagged->SetYTitle("N");
      fhPtAntiProtonTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtAntiProtonTagged) ; 
      
      fhPtUnknownTagged  = new TH1F("hPtMCUnknownTagged","Number of Unknown id as Photon over calorimeter, tagged as converted",nptbins,ptmin,ptmax); 
      fhPtUnknownTagged->SetYTitle("N");
      fhPtUnknownTagged->SetXTitle("p_{T #gamma}(GeV/c)");
      outputContainer->Add(fhPtUnknownTagged) ;     
      
      fhConvDeltaEtaMCConversion  = new TH2F
      ("hConvDeltaEtaMCConversion","#Delta #eta of selected conversion pairs from real conversions",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCConversion->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCConversion) ;
      
      fhConvDeltaPhiMCConversion  = new TH2F
      ("hConvDeltaPhiMCConversion","#Delta #phi of selected conversion pairs from real conversions",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCConversion->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCConversion) ;
      
      fhConvDeltaEtaPhiMCConversion  = new TH2F
      ("hConvDeltaEtaPhiMCConversion","#Delta #eta vs #Delta #phi of selected conversion pairs, from real conversions",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCConversion->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCConversion->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCConversion) ;
      
      fhConvAsymMCConversion  = new TH2F
      ("hConvAsymMCConversion","Asymmetry of selected conversion pairs from real conversions",100,0,fMassCut,100,0,1); 
      fhConvAsymMCConversion->SetYTitle("Asymmetry");
      fhConvAsymMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCConversion) ;
      
      fhConvPtMCConversion  = new TH2F
      ("hConvPtMCConversion","p_{T} of selected conversion pairs from real conversions",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCConversion->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCConversion->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCConversion) ;    
      
      fhConvDispersionMCConversion  = new TH2F
      ("hConvDispersionMCConversion","p_{T} of selected conversion pairs from real conversions",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCConversion->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCConversion->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCConversion) ;   
      
      fhConvM02MCConversion  = new TH2F
      ("hConvM02MCConversion","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCConversion->SetYTitle("M02 cluster 1");
      fhConvM02MCConversion->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCConversion) ;           
      
      fhConvDeltaEtaMCAntiNeutron  = new TH2F
      ("hConvDeltaEtaMCAntiNeutron","#Delta #eta of selected conversion pairs from anti-neutrons",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCAntiNeutron->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCAntiNeutron) ;
      
      fhConvDeltaPhiMCAntiNeutron  = new TH2F
      ("hConvDeltaPhiMCAntiNeutron","#Delta #phi of selected conversion pairs from anti-neutrons",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCAntiNeutron) ;
      
      fhConvDeltaEtaPhiMCAntiNeutron  = new TH2F
      ("hConvDeltaEtaPhiMCAntiNeutron","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-neutrons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCAntiNeutron->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCAntiNeutron->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCAntiNeutron) ;    
      
      fhConvAsymMCAntiNeutron  = new TH2F
      ("hConvAsymMCAntiNeutron","Asymmetry of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0,1); 
      fhConvAsymMCAntiNeutron->SetYTitle("Asymmetry");
      fhConvAsymMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCAntiNeutron) ;
      
      fhConvPtMCAntiNeutron  = new TH2F
      ("hConvPtMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCAntiNeutron->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCAntiNeutron->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCAntiNeutron) ;    
      
      fhConvDispersionMCAntiNeutron  = new TH2F
      ("hConvDispersionMCAntiNeutron","p_{T} of selected conversion pairs from anti-neutrons",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCAntiNeutron->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCAntiNeutron->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCAntiNeutron) ;       
      
      fhConvM02MCAntiNeutron  = new TH2F
      ("hConvM02MCAntiNeutron","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCAntiNeutron->SetYTitle("M02 cluster 1");
      fhConvM02MCAntiNeutron->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCAntiNeutron) ;  
      
      fhConvDeltaEtaMCAntiProton  = new TH2F
      ("hConvDeltaEtaMCAntiProton","#Delta #eta of selected conversion pairs from anti-protons",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCAntiProton->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCAntiProton) ;
      
      fhConvDeltaPhiMCAntiProton  = new TH2F
      ("hConvDeltaPhiMCAntiProton","#Delta #phi of selected conversion pairs from anti-protons",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCAntiProton->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCAntiProton) ;
      
      fhConvDeltaEtaPhiMCAntiProton  = new TH2F
      ("hConvDeltaEtaPhiMCAntiProton","#Delta #eta vs #Delta #phi of selected conversion pairs from anti-protons",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCAntiProton->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCAntiProton->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCAntiProton) ;    
      
      fhConvAsymMCAntiProton  = new TH2F
      ("hConvAsymMCAntiProton","Asymmetry of selected conversion pairs from anti-protons",100,0,fMassCut,100,0,1); 
      fhConvAsymMCAntiProton->SetYTitle("Asymmetry");
      fhConvAsymMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCAntiProton) ;
      
      fhConvPtMCAntiProton  = new TH2F
      ("hConvPtMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCAntiProton->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCAntiProton->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCAntiProton) ;
      
      fhConvDispersionMCAntiProton  = new TH2F
      ("hConvDispersionMCAntiProton","p_{T} of selected conversion pairs from anti-protons",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCAntiProton->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCAntiProton->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCAntiProton) ;       
      
      fhConvM02MCAntiProton  = new TH2F
      ("hConvM02MCAntiProton","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCAntiProton->SetYTitle("M02 cluster 1");
      fhConvM02MCAntiProton->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCAntiProton) ;       
      
      fhConvDeltaEtaMCString  = new TH2F
      ("hConvDeltaEtaMCString","#Delta #eta of selected conversion pairs from string",100,0,fMassCut,netabins,-0.5,0.5); 
      fhConvDeltaEtaMCString->SetYTitle("#Delta #eta");
      fhConvDeltaEtaMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaEtaMCString) ;
      
      fhConvDeltaPhiMCString  = new TH2F
      ("hConvDeltaPhiMCString","#Delta #phi of selected conversion pairs from string",100,0,fMassCut,nphibins,-0.5,0.5); 
      fhConvDeltaPhiMCString->SetYTitle("#Delta #phi");
      fhConvDeltaPhiMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvDeltaPhiMCString) ;
      
      fhConvDeltaEtaPhiMCString  = new TH2F
      ("hConvDeltaEtaPhiMCString","#Delta #eta vs #Delta #phi of selected conversion pairs from string",netabins,-0.5,0.5,nphibins,-0.5,0.5); 
      fhConvDeltaEtaPhiMCString->SetYTitle("#Delta #phi");
      fhConvDeltaEtaPhiMCString->SetXTitle("#Delta #eta");
      outputContainer->Add(fhConvDeltaEtaPhiMCString) ;    
      
      fhConvAsymMCString  = new TH2F
      ("hConvAsymMCString","Asymmetry of selected conversion pairs from string",100,0,fMassCut,100,0,1); 
      fhConvAsymMCString->SetYTitle("Asymmetry");
      fhConvAsymMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvAsymMCString) ;
      
      fhConvPtMCString  = new TH2F
      ("hConvPtMCString","p_{T} of selected conversion pairs from string",100,0,fMassCut,100,0.,10.); 
      fhConvPtMCString->SetYTitle("Pair p_{T} (GeV/c)");
      fhConvPtMCString->SetXTitle("Pair Mass (GeV/c^2)");
      outputContainer->Add(fhConvPtMCString) ;
      
      fhConvDispersionMCString  = new TH2F
      ("hConvDispersionMCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvDispersionMCString->SetYTitle("Dispersion cluster 1");
      fhConvDispersionMCString->SetXTitle("Dispersion cluster 2");
      outputContainer->Add(fhConvDispersionMCString) ;       
      
      fhConvM02MCString  = new TH2F
      ("hConvM02MCString","p_{T} of selected conversion pairs from string",100,0.,1.,100,0.,1.); 
      fhConvM02MCString->SetYTitle("M02 cluster 1");
      fhConvM02MCString->SetXTitle("M02 cluster 2");
      outputContainer->Add(fhConvM02MCString) ; 
      
      fhConvDistMCConversion  = new TH2F
      ("hConvDistMCConversion","calculated conversion distance vs real vertes for MC conversion",100,0.,5.,100,0.,5.); 
      fhConvDistMCConversion->SetYTitle("distance");
      fhConvDistMCConversion->SetXTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversion) ; 
      
      fhConvDistMCConversionCuts  = new TH2F
      ("hConvDistMCConversionCuts","calculated conversion distance vs real vertes for MC conversion, deta < 0.05, m < 10 MeV, asym < 0.1",100,0.,5.,100,0.,5.); 
      fhConvDistMCConversionCuts->SetYTitle("distance");
      fhConvDistMCConversionCuts->SetXTitle("vertex R");
      outputContainer->Add(fhConvDistMCConversionCuts) ; 
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
    
  //Store calo PID histograms
  if(fRejectTrackMatch && IsCaloPIDOn()){
    TList * caloPIDHistos = GetCaloPID()->GetCreateOutputObjects() ;
    for(Int_t i = 0; i < caloPIDHistos->GetEntries(); i++) {
      outputContainer->Add(caloPIDHistos->At(i)) ;
    }
    delete caloPIDHistos;
  }
  
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
  fMassCut     = 0.03; //30 MeV
	
  fTimeCutMin  = -1;
  fTimeCutMax  = 9999999;
  fNCellsCut   = 0;
	
  fRejectTrackMatch       = kTRUE ;
  fCheckConversion        = kFALSE;
  fRemoveConvertedPair    = kFALSE;
  fAddConvertedPairsToAOD = kFALSE;
	
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
  if(fCalorimeter == "PHOS")
    pl = GetPHOSClusters();
  else if (fCalorimeter == "EMCAL")
    pl = GetEMCALClusters();
  
  if(!pl) {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
  
  //Init arrays, variables, get number of clusters
  TLorentzVector mom, mom2 ;
  Int_t nCaloClusters = pl->GetEntriesFast();
  //List to be used in conversion analysis, to tag the cluster as candidate for conversion
  Bool_t * indexConverted = 0x0;
  if(fCheckConversion){
    indexConverted = new Bool_t[nCaloClusters];
    for (Int_t i = 0; i < nCaloClusters; i++) 
      indexConverted[i] = kFALSE;
	}
  
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
    //Play with the MC stack if available
    //--------------------------------------------------------------------------------------
    
    //Check origin of the candidates
    if(IsDataMC()){
      aodph.SetTag(GetMCAnalysisUtils()->CheckOrigin(calo->GetLabels(),calo->GetNLabels(),GetReader(), aodph.GetInputFileIndex()));

      if(GetDebug() > 0)
        printf("AliAnaPhoton::MakeAnalysisFillAOD() - Origin of candidate, bit map %d\n",aodph.GetTag());
    }//Work with stack also   
    
    //--------------------------------------------------------------------------------------
    //Fill some shower shape histograms before PID is applied
    //--------------------------------------------------------------------------------------
    
    FillShowerShapeHistograms(calo,aodph.GetTag());
    
    //-------------------------------------
    //PID selection or bit setting
    //-------------------------------------
    // MC
    if(GetReader()->GetDataType() == AliCaloTrackReader::kMC){
      //Get most probable PID, check PID weights (in MC this option is mandatory)
      aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());	 
      //If primary is not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;
    }	
    //...............................................
    // Data, PID check on
    else if(IsCaloPIDOn()){
      //Get most probable PID, 2 options check PID weights 
      //or redo PID, recommended option for MCEal.		
      if(!IsCaloPIDRecalculationOn())
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,calo->GetPID(),mom.E()));//PID with weights
      else
        aodph.SetIdentifiedParticleType(GetCaloPID()->GetIdentifiedParticleType(fCalorimeter,mom,calo));//PID recalculated
      
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PDG of identified particle %d\n",aodph.GetIdentifiedParticleType());
      
      //If cluster does not pass pid, not photon, skip it.
      if(aodph.GetIdentifiedParticleType() != AliCaloPID::kPhoton) continue ;			
      
    }
    //...............................................
    // Data, PID check off
    else{
      //Set PID bits for later selection (AliAnaPi0 for example)
      //GetPDG already called in SetPIDBits.
      GetCaloPID()->SetPIDBits(fCalorimeter,calo,&aodph, GetCaloUtils());
      if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - PID Bits set \n");		
    }
    
    if(GetDebug() > 1) printf("AliAnaPhoton::MakeAnalysisFillAOD() - Photon selection cuts passed: pT %3.2f, pdg %d\n",aodph.Pt(), aodph.GetIdentifiedParticleType());
    
    
    //--------------------------------------------------------------------------------------
    // Conversions pairs analysis
    // Check if cluster comes from a conversion in the material in front of the calorimeter
    // Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
    //--------------------------------------------------------------------------------------
    
    // Do analysis only if there are more than one cluster
    if( nCaloClusters > 1 && fCheckConversion){
      Bool_t bConverted = kFALSE;
      Int_t id2 = -1;
		  
      //Check if set previously as converted couple, if so skip its use.
      if (indexConverted[icalo]) continue;
		  
      // Second cluster loop
      for(Int_t jcalo = icalo + 1 ; jcalo < nCaloClusters ; jcalo++) {
        //Check if set previously as converted couple, if so skip its use.
        if (indexConverted[jcalo]) continue;
        //printf("Check Conversion indeces %d and %d\n",icalo,jcalo);
        AliVCluster * calo2 =  (AliVCluster*) (pl->At(jcalo));  //Get cluster kinematics
        
        
        //Mixed event, get index of event
        Int_t evtIndex2 = 0 ; 
        if (GetMixedEvent()) {
          evtIndex2=GetMixedEvent()->EventIndexForCaloCluster(calo2->GetID()) ; 
          
        }      
        
        //Get kinematics of second cluster
        if(GetReader()->GetDataType() != AliCaloTrackReader::kMC){
          calo2->GetMomentum(mom2,GetVertex(evtIndex2)) ;}//Assume that come from vertex in straight line
        else{
          Double_t vertex[]={0,0,0};
          calo2->GetMomentum(mom2,vertex) ;
        }
        
        //--------------------------------------
        // Cluster selection
        //--------------------------------------
        
        if(!ClusterSelected(calo2,mom2)) continue;  
        
        //................................................
        // Get TOF of each cluster in pair, calculate difference if small, 
        // take this pair. Way to reject clusters from hadrons (or pileup?)
        Double_t t12diff = calo2->GetTOF()-calo->GetTOF()*1e9;
        if(TMath::Abs(t12diff) > GetPairTimeCut()) continue;
        
        //................................................
        //Get mass of pair, if small, take this pair.
        Float_t pairM     = (mom+mom2).M();
        //printf("\t both in calo, mass %f, cut %f\n",pairM,fMassCut);
        if(pairM < fMassCut){  
          aodph.SetTagged(kFALSE);
          id2 = calo2->GetID();
          indexConverted[icalo]=kTRUE;
          indexConverted[jcalo]=kTRUE; 
          Float_t asymmetry = TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E());
          Float_t dPhi      = mom.Phi()-mom2.Phi();
          Float_t dEta      = mom.Eta()-mom2.Eta();  
          
          //...............................................
          //Fill few histograms with kinematics of the pair
          //FIXME, move all this to MakeAnalysisFillHistograms ...
          
          fhConvDeltaEta   ->Fill( pairM, dPhi      );
          fhConvDeltaPhi   ->Fill( pairM, dEta      );
          fhConvAsym       ->Fill( pairM, asymmetry );
          fhConvDeltaEtaPhi->Fill( dEta , dPhi      );
          fhConvPt         ->Fill( pairM, (mom+mom2).Pt());          
          
          //Estimate conversion distance, T. Awes, M. Ivanov
          //Under the assumption that the pair has zero mass, and that each electron 
          //of the pair has the same momentum, they will each have the same bend radius 
          //given by R=p/(qB) = p / (300 B) with p in [MeV/c], B in [Tesla] and R in [m]. 
          //With nominal ALICE magnet current of 30kA B=0.5T, and so with E_cluster=p,  
          //R = E/1.5 [cm].  Under these assumptions, the distance from the conversion 
          //point to the MCEal can be related to the separation distance, L=2y, on the MCEal 
          //as d = sqrt(R^2 -(R-y)^2) = sqrt(2Ry - y^2). And since R>>y we can write as 
          //d = sqrt(E*L/1.5) where E is the cluster energy and L is the distance in cm between 
          //the clusters.
          Float_t pos1[3];
          calo->GetPosition(pos1); 
          Float_t pos2[3];
          calo2->GetPosition(pos2); 
          Float_t clustDist = TMath::Sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+
                                          (pos1[1]-pos2[1])*(pos1[1]-pos2[1])+
                                          (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
          
          Float_t convDist  = TMath::Sqrt(mom.E() *clustDist*0.01/0.15);
          Float_t convDist2 = TMath::Sqrt(mom2.E()*clustDist*0.01/0.15);
          //printf("l = %f, e1 = %f, d1=%f, e2 = %f, d2=%f\n",clustDist,mom.E(),convDist,mom2.E(),convDist2);
          if(GetDebug() > 2)
            printf("AliAnaPhoton::MakeAnalysisFillAOD(): Pair with mass %2.3f < %2.3f, %1.2f < dPhi %2.2f < %2.2f, dEta %f < %2.2f, asymmetry %2.2f< %2.2f; \n    cluster1 id %d, e %2.3f  SM %d, eta %2.3f, phi %2.3f ; \n    cluster2 id %d, e %2.3f, SM %d,eta %2.3f, phi %2.3f\n",
                   pairM,fMassCut,fConvDPhiMinCut, dPhi, fConvDPhiMaxCut, dEta, fConvDEtaCut, asymmetry, fConvAsymCut,
                   calo->GetID(),calo->E(),GetCaloUtils()->GetModuleNumber(calo), mom.Eta(), mom.Phi(),
                   id2, calo2->E(), GetCaloUtils()->GetModuleNumber(calo2),mom2.Eta(), mom2.Phi());
          
          fhConvDistEta ->Fill(mom .Eta(), convDist );
          fhConvDistEta ->Fill(mom2.Eta(), convDist2);
          fhConvDistEn  ->Fill(mom .E(),   convDist );
          fhConvDistEn  ->Fill(mom2.E(),   convDist2);        
          fhConvDistMass->Fill((mom+mom2).M(), convDist );
          //dEta cut
          if(dEta<0.05){
            fhConvDistEtaCutEta ->Fill(mom .Eta(), convDist );
            fhConvDistEtaCutEta ->Fill(mom2.Eta(), convDist2);
            fhConvDistEnCutEta  ->Fill(mom .E(),   convDist );
            fhConvDistEnCutEta  ->Fill(mom2.E(),   convDist2);        
            fhConvDistMassCutEta->Fill((mom+mom2).M(), convDist );
            //mass cut
            if(pairM<0.01){//10 MeV
              fhConvDistEtaCutMass ->Fill(mom .Eta(), convDist );
              fhConvDistEtaCutMass ->Fill(mom2.Eta(), convDist2);
              fhConvDistEnCutMass  ->Fill(mom .E(),   convDist );
              fhConvDistEnCutMass  ->Fill(mom2.E(),   convDist2);        
              // asymmetry cut
              if(asymmetry<0.1){
                fhConvDistEtaCutAsy ->Fill(mom .Eta(), convDist );
                fhConvDistEtaCutAsy ->Fill(mom2.Eta(), convDist2);
                fhConvDistEnCutAsy  ->Fill(mom .E(),   convDist );
                fhConvDistEnCutAsy  ->Fill(mom2.E(),   convDist2); 
              }//asymmetry cut
            }//mass cut            
          }//dEta cut
          
          //...............................................
          //Select pairs in a eta-phi window
          if(TMath::Abs(dEta) < fConvDEtaCut    && 
             TMath::Abs(dPhi) < fConvDPhiMaxCut &&
             TMath::Abs(dPhi) > fConvDPhiMinCut && 
             asymmetry        < fConvAsymCut       ){
            bConverted = kTRUE;          
          }
          //printf("Accepted? %d\n",bConverted);
          //...........................................
          //Fill more histograms, simulated data
          //FIXME, move all this to MakeAnalysisFillHistograms ...
          if(IsDataMC()){
            
            //Check the origin of the pair, look for conversion, antinucleons or jet correlations (strings)
            Int_t ancPDG    = 0;
            Int_t ancStatus = 0;
            TLorentzVector momentum;
            TVector3 prodVertex;
            Int_t ancLabel  = GetMCAnalysisUtils()->CheckCommonAncestor(calo->GetLabel(), calo2->GetLabel(), 
                                                                        GetReader(), ancPDG, ancStatus, momentum, prodVertex);
            
            // printf("AliAnaPhoton::MakeAnalysisFillHistograms() - Common ancestor label %d, pdg %d, name %s, status %d; \n",
            //                          ancLabel,ancPDG,TDatabasePDG::Instance()->GetParticle(ancPDG)->GetName(),ancStatus);
            
            Int_t tag2 = GetMCAnalysisUtils()->CheckOrigin(calo2->GetLabels(),calo2->GetNLabels(),GetReader(), 0);
            if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCConversion)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCConversion) && (ancPDG==22 || TMath::Abs(ancPDG)==11) && ancLabel > -1){
                fhConvDeltaEtaMCConversion   ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCConversion   ->Fill( pairM, dPhi      );
                fhConvAsymMCConversion       ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCConversion->Fill( dEta , dPhi      );
                fhConvPtMCConversion         ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCConversion ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCConversion        ->Fill( calo->GetM02(), calo2->GetM02());
                fhConvDistMCConversion       ->Fill( convDist , prodVertex.Mag() );
                fhConvDistMCConversion       ->Fill( convDist2, prodVertex.Mag() );
                
                if(dEta<0.05 && pairM<0.01 && asymmetry<0.1){
                  fhConvDistMCConversionCuts->Fill( convDist , prodVertex.Mag() );
                  fhConvDistMCConversionCuts->Fill( convDist2, prodVertex.Mag() );
                }
                
              }              
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiNeutron)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiNeutron) && ancPDG==-2112 && ancLabel > -1){
                fhConvDeltaEtaMCAntiNeutron    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiNeutron    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiNeutron        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiNeutron ->Fill( dEta , dPhi      );
                fhConvPtMCAntiNeutron          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiNeutron  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiNeutron         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            else if(GetMCAnalysisUtils()->CheckTagBit(aodph.GetTag(),AliMCAnalysisUtils::kMCAntiProton)){
              if(GetMCAnalysisUtils()->CheckTagBit(tag2,AliMCAnalysisUtils::kMCAntiProton) && ancPDG==-2212 && ancLabel > -1){
                fhConvDeltaEtaMCAntiProton    ->Fill( pairM, dEta      );
                fhConvDeltaPhiMCAntiProton    ->Fill( pairM, dPhi      );
                fhConvAsymMCAntiProton        ->Fill( pairM, asymmetry );
                fhConvDeltaEtaPhiMCAntiProton ->Fill( dEta , dPhi      );
                fhConvPtMCAntiProton          ->Fill( pairM, (mom+mom2).Pt());
                fhConvDispersionMCAntiProton  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
                fhConvM02MCAntiProton         ->Fill( calo->GetM02(), calo2->GetM02());
              }
            }
            
            //Pairs coming from fragmenting pairs.
            if(ancPDG < 22 && ancLabel > 7 && (ancStatus == 11 || ancStatus == 12) ){
              fhConvDeltaEtaMCString    ->Fill( pairM, dPhi);
              fhConvDeltaPhiMCString    ->Fill( pairM, dPhi);
              fhConvAsymMCString        ->Fill( pairM, TMath::Abs(mom.E()-mom2.E())/(mom.E()+mom2.E()) );
              fhConvDeltaEtaPhiMCString ->Fill( dPhi, dPhi );
              fhConvPtMCString          ->Fill( pairM, (mom+mom2).Pt());
              fhConvDispersionMCString  ->Fill( calo->GetDispersion(), calo2->GetDispersion());
              fhConvM02MCString         ->Fill( calo->GetM02(), calo2->GetM02());
            }
            
          }// Data MC
          
          break;
        }
			  
      }//Mass loop
		  
      //..........................................................................................................
      //Pair selected as converted, remove both clusters or recombine them into a photon and put them in the AOD
      if(bConverted){ 
        //Add to AOD
        if(fAddConvertedPairsToAOD){
          //Create AOD of pair analysis
          TLorentzVector mpair = mom+mom2;
          AliAODPWG4Particle aodpair = AliAODPWG4Particle(mpair);
          aodpair.SetLabel(aodph.GetLabel());
          //aodpair.SetInputFileIndex(input);
          
          //printf("Index %d, Id %d\n",icalo, calo->GetID());
          //Set the indeces of the original caloclusters  
          aodpair.SetCaloLabel(calo->GetID(),id2);
          aodpair.SetDetector(fCalorimeter);
          aodpair.SetIdentifiedParticleType(aodph.GetIdentifiedParticleType());
          aodpair.SetTag(aodph.GetTag());
          aodpair.SetTagged(kTRUE);
          //Add AOD with pair object to aod branch
          AddAODParticle(aodpair);
          //printf("\t \t both added pair\n");
        }
        
        //Do not add the current calocluster
        if(fRemoveConvertedPair) continue;
        else {
          //printf("TAGGED\n");
          //Tag this cluster as likely conversion
          aodph.SetTagged(kTRUE);
        }
      }//converted pair
    }//check conversion
    //printf("\t \t added single cluster %d\n",icalo);
	  
    //FIXME, this to MakeAnalysisFillHistograms ...
    Int_t absID             = 0; 
    Float_t maxCellFraction = 0;
    AliVCaloCells* cells    = 0;
    
    if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
    else                        cells = GetPHOSCells();
    
    absID = GetCaloUtils()->GetMaxEnergyCell(cells, calo,maxCellFraction);
    
    fhMaxCellDiffClusterE->Fill(aodph.E(),maxCellFraction);
    fhNCellsE            ->Fill(aodph.E(),calo->GetNCells());

    //Add AOD with photon object to aod branch
    AddAODParticle(aodph);
    
  }//loop
  
  delete [] indexConverted;
	
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
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ph =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ph->GetIdentifiedParticleType();
    
    if(GetDebug() > 3) 
      printf("AliAnaPhoton::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ph->GetIdentifiedParticleType(),ph->GetTag(), (ph->GetDetector()).Data()) ;
    
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
    
    if(fCheckConversion &&ph->IsTagged()){
      fhPtPhotonConv->Fill(ptcluster);
      if(ecluster > 0.5)        fhEtaPhiPhotonConv  ->Fill(etacluster, phicluster);
      else if(GetMinPt() < 0.5) fhEtaPhi05PhotonConv->Fill(etacluster, phicluster);
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
      
      if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && fhMCE[mcPhoton])
      {
        fhMCE  [mcPhoton] ->Fill(ecluster);
        fhMCPt [mcPhoton] ->Fill(ptcluster);
        fhMCPhi[mcPhoton] ->Fill(ecluster,phicluster);
        fhMCEta[mcPhoton] ->Fill(ecluster,etacluster);
        
        fhMC2E[mcPhoton]     ->Fill(ecluster, eprim);
        fhMC2Pt[mcPhoton]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[mcPhoton] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcPhoton]->Fill(ptcluster,ptprim-ptcluster);     
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion) && fhMCE[mcConversion])
        {
          fhMCE  [mcConversion] ->Fill(ecluster);
          fhMCPt [mcConversion] ->Fill(ptcluster);
          fhMCPhi[mcConversion] ->Fill(ecluster,phicluster);
          fhMCEta[mcConversion] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcConversion]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcConversion]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcConversion] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcConversion]->Fill(ptcluster,ptprim-ptcluster);     
          
          if(fCheckConversion){
            if(ph->IsTagged()) fhPtConversionTagged ->Fill(ptcluster);
            if(ptcluster > 0.5)fhEtaPhiConversion   ->Fill(etacluster,phicluster);
            else               fhEtaPhi05Conversion ->Fill(etacluster,phicluster);
          }
        }			
        
        if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt) && fhMCE[mcPrompt]){
          fhMCE  [mcPrompt] ->Fill(ecluster);
          fhMCPt [mcPrompt] ->Fill(ptcluster);
          fhMCPhi[mcPrompt] ->Fill(ecluster,phicluster);
          fhMCEta[mcPrompt] ->Fill(ecluster,etacluster);      
          
          fhMC2E[mcPrompt]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcPrompt]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcPrompt] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcPrompt]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation)&& fhMCE[mcFragmentation])
        {
          fhMCE  [mcFragmentation] ->Fill(ecluster);
          fhMCPt [mcFragmentation] ->Fill(ptcluster);
          fhMCPhi[mcFragmentation] ->Fill(ecluster,phicluster);
          fhMCEta[mcFragmentation] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcFragmentation]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcFragmentation]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcFragmentation] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcFragmentation]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR)&& fhMCE[mcISR])
        {
          fhMCE  [mcISR] ->Fill(ecluster);
          fhMCPt [mcISR] ->Fill(ptcluster);
          fhMCPhi[mcISR] ->Fill(ecluster,phicluster);
          fhMCEta[mcISR] ->Fill(ecluster,etacluster);    
          
          fhMC2E[mcISR]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcISR]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcISR] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcISR]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) && 
                !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE[mcPi0Decay])
        {
          fhMCE  [mcPi0Decay] ->Fill(ecluster);
          fhMCPt [mcPi0Decay] ->Fill(ptcluster);
          fhMCPhi[mcPi0Decay] ->Fill(ecluster,phicluster);
          fhMCEta[mcPi0Decay] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcPi0Decay]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcPi0Decay]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcPi0Decay] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcPi0Decay]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if( (GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || 
                  GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay) ) && fhMCE[mcOtherDecay])
        {
          fhMCE  [mcOtherDecay] ->Fill(ecluster);
          fhMCPt [mcOtherDecay] ->Fill(ptcluster);
          fhMCPhi[mcOtherDecay] ->Fill(ecluster,phicluster);
          fhMCEta[mcOtherDecay] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcOtherDecay]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcOtherDecay]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcOtherDecay] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcOtherDecay]->Fill(ptcluster,ptprim-ptcluster);     
          
        }
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) && fhMCE  [mcPi0])
        {
          fhMCE  [mcPi0] ->Fill(ecluster);
          fhMCPt [mcPi0] ->Fill(ptcluster);
          fhMCPhi[mcPi0] ->Fill(ecluster,phicluster);
          fhMCEta[mcPi0] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcPi0]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcPi0]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcPi0] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcPi0]->Fill(ptcluster,ptprim-ptcluster);     
          
        } 
        else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta) && fhMCE[mcEta])
        {
          fhMCE  [mcEta] ->Fill(ecluster);
          fhMCPt [mcEta] ->Fill(ptcluster);
          fhMCPhi[mcEta] ->Fill(ecluster,phicluster);
          fhMCEta[mcEta] ->Fill(ecluster,etacluster);
          
          fhMC2E[mcEta]     ->Fill(ecluster, eprim);
          fhMC2Pt[mcEta]    ->Fill(ptcluster, ptprim);     
          fhMCDeltaE[mcEta] ->Fill(ecluster,eprim-ecluster);
          fhMCDeltaPt[mcEta]->Fill(ptcluster,ptprim-ptcluster);     
          
        }      
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron) && fhMCE[mcAntiNeutron])
      {
        fhMCE  [mcAntiNeutron] ->Fill(ecluster);
        fhMCPt [mcAntiNeutron] ->Fill(ptcluster);
        fhMCPhi[mcAntiNeutron] ->Fill(ecluster,phicluster);
        fhMCEta[mcAntiNeutron] ->Fill(ecluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtAntiNeutronTagged ->Fill(ptcluster);
        
        fhMC2E[mcAntiNeutron]     ->Fill(ecluster, eprim);
        fhMC2Pt[mcAntiNeutron]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[mcAntiNeutron] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcAntiNeutron]->Fill(ptcluster,ptprim-ptcluster);     
        
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton) && fhMCE[mcAntiProton])
      {
        fhMCE  [mcAntiProton] ->Fill(ecluster);
        fhMCPt [mcAntiProton] ->Fill(ptcluster);
        fhMCPhi[mcAntiProton] ->Fill(ecluster,phicluster);
        fhMCEta[mcAntiProton] ->Fill(ecluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtAntiProtonTagged ->Fill(ptcluster);
        
        fhMC2E[mcAntiProton]     ->Fill(ecluster, eprim);
        fhMC2Pt[mcAntiProton]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[mcAntiProton] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcAntiProton]->Fill(ecluster,ptprim-ptcluster);     
        
      } 
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron) && fhMCE[mcElectron])
      {
        fhMCE  [mcElectron] ->Fill(ecluster);
        fhMCPt [mcElectron] ->Fill(ptcluster);
        fhMCPhi[mcElectron] ->Fill(ecluster,phicluster);
        fhMCEta[mcElectron] ->Fill(ecluster,etacluster);
        
        fhMC2E[mcElectron]     ->Fill(ecluster, eprim);
        fhMC2Pt[mcElectron]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[mcElectron] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcElectron]->Fill(ecluster,ptprim-ptcluster);             
      }     
      else if( fhMCE[mcOther]){
        fhMCE  [mcOther] ->Fill(ecluster);
        fhMCPt [mcOther] ->Fill(ptcluster);
        fhMCPhi[mcOther] ->Fill(ecluster,phicluster);
        fhMCEta[mcOther] ->Fill(ecluster,etacluster);
        if(ph->IsTagged() && fCheckConversion) fhPtUnknownTagged ->Fill(ptcluster);
        
        fhMC2E[mcOther]     ->Fill(ecluster, eprim);
        fhMC2Pt[mcOther]    ->Fill(ptcluster, ptprim);     
        fhMCDeltaE[mcOther] ->Fill(ecluster,eprim-ecluster);
        fhMCDeltaPt[mcOther]->Fill(ecluster,ptprim-ptcluster);     
        
        
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
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("Min Distance to Bad Channel   = %2.1f\n",fMinDist);
  printf("Min Distance to Bad Channel 2 = %2.1f\n",fMinDist2);
  printf("Min Distance to Bad Channel 3 = %2.1f\n",fMinDist3);
  printf("Reject clusters with a track matched = %d\n",fRejectTrackMatch);
  printf("Check Pair Conversion                = %d\n",fCheckConversion);
  printf("Add conversion pair to AOD           = %d\n",fAddConvertedPairsToAOD);
  printf("Conversion pair mass cut             = %f\n",fMassCut);
  printf("Conversion selection cut : A < %1.2f; %1.3f < Dphi < %1.3f; Deta < %1.3f\n",
         fConvAsymCut,fConvDPhiMinCut, fConvDPhiMaxCut, fConvDEtaCut);
  printf("Time Cut: %3.1f < TOF  < %3.1f\n", fTimeCutMin, fTimeCutMax);
  printf("Number of cells in cluster is        > %d \n", fNCellsCut);
  printf("    \n") ;
	
} 
