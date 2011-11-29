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
/* $Id:  $ */

//_________________________________________________________________________
// Class containing methods for the isolation cut. 
// An AOD candidate (AliAODPWG4ParticleCorrelation type)
// is passed. Look in a cone around the candidate and study
// the hadronic activity inside to decide if the candidate is isolated
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 

//-Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)
//-Yaxian Mao (check the candidate particle is the leading particle or not at the same hemishere)

//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TLorentzVector.h>
#include <TObjArray.h>

// --- AliRoot system --- 
#include "AliIsolationCut.h" 
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliCaloTrackReader.h"
#include "AliMixedEvent.h"

ClassImp(AliIsolationCut)
  
//____________________________________
AliIsolationCut::AliIsolationCut() : 
TObject(),
fConeSize(0.),
fPtThreshold(0.), 
fSumPtThreshold(0.), 
fPtFraction(0.), 
fICMethod(0),
fPartInCone(0)

{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
}

//____________________________________________
TString AliIsolationCut::GetICParametersList()
{
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliIsolationCut ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fConeSize: (isolation cone size) %1.2f\n",fConeSize) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtThreshold =%1.2f (isolation pt threshold) \n",fPtThreshold) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtFraction=%1.2f (isolation pt threshold fraction ) \n",fPtFraction) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fICMethod=%d (isolation cut case) \n",fICMethod) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPartInCone=%d \n",fPartInCone) ;
  parList+=onePar ;
  
  return parList; 
}

//____________________________________________________________________________
void AliIsolationCut::InitParameters()
{
  //Initialize the parameters of the analysis.
  
  fConeSize       = 0.4 ; 
  fPtThreshold    = 1.  ; 
  fSumPtThreshold = 0.5 ; 
  fPtFraction     = 0.1 ; 
  fPartInCone     = kOnlyCharged;
  fICMethod       = kSumPtFracIC; // 0 pt threshol method, 1 cone pt sum method
  
}

//________________________________________________________________________________
void  AliIsolationCut::MakeIsolationCut(const TObjArray * plCTS, 
                                        const TObjArray * plNe, 
                                        const AliCaloTrackReader * reader, 
                                        const Bool_t bFillAOD, 
                                        AliAODPWG4ParticleCorrelation  *pCandidate, 
                                        const TString & aodArrayRefName,
                                        Int_t & n, 
                                        Int_t & nfrac, 
                                        Float_t &coneptsum,  
                                        Bool_t  &isolated) const
{  
  //Search in cone around a candidate particle if it is isolated 
  Float_t phiC  = pCandidate->Phi() ;
  if(phiC<0) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  Float_t ptC   = pCandidate->Pt() ;
  Float_t pt    = -100. ;
  Float_t eta   = -100. ;
  Float_t phi   = -100. ;
  Float_t rad   = -100. ;
  
  n         = 0 ;
  nfrac     = 0 ;
  coneptsum = 0.; 
  isolated  = kFALSE;
  
  //Initialize the array with refrences
  TObjArray * refclusters = 0x0;
  TObjArray * reftracks   = 0x0;
  Int_t ntrackrefs   = 0;
  Int_t nclusterrefs = 0;
  //Check charged particles in cone.
  if(plCTS && (fPartInCone==kOnlyCharged || fPartInCone==kNeutralAndCharged)){
    TVector3 p3;
    for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ ){
      AliAODTrack* track = (AliAODTrack *)(plCTS->At(ipr)) ; 
      //Do not count the candidate (pion, conversion photon) or the daughters of the candidate
      if(track->GetID() == pCandidate->GetTrackLabel(0) || track->GetID() == pCandidate->GetTrackLabel(1) 
         || track->GetID() == pCandidate->GetTrackLabel(2) || track->GetID() == pCandidate->GetTrackLabel(3) 
         ) continue ;
      p3.SetXYZ(track->Px(),track->Py(),track->Pz());
      pt   = p3.Pt();
      eta  = p3.Eta();
      phi  = p3.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      
      //only loop the particle at the same side of candidate
      if(TMath::Abs(phi-phiC)>TMath::PiOver2()) continue ;
      //if at the same side has particle larger than candidate, then candidate can not be the leading, skip such events
      if(pt > ptC){
        n         = -1;
        nfrac     = -1;
        coneptsum = -1;
        isolated  = kFALSE;
        if(bFillAOD && reftracks) {
          reftracks->Clear(); 
          delete reftracks;
        }
        return ;
      }
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      
      if(rad < fConeSize){
        if(bFillAOD) {
          ntrackrefs++;
          if(ntrackrefs == 1){
            reftracks = new TObjArray(0);
            //reftracks->SetName(Form("Tracks%s",aodArrayRefName.Data()));
            TString tempo(aodArrayRefName)  ; 
            tempo += "Tracks" ; 
            reftracks->SetName(tempo);
            reftracks->SetOwner(kFALSE);
          }
          reftracks->Add(track);
        }
        //printf("charged in isolation cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
        coneptsum+=pt;
        if(pt > fPtThreshold )    n++;
        if(pt > fPtFraction*ptC ) nfrac++;  
      } // Inside cone
    }// charged particle loop
  }//Tracks
  
  //Check neutral particles in cone.  
  if(plNe && (fPartInCone==kOnlyNeutral || fPartInCone==kNeutralAndCharged)){
	  
    
    TLorentzVector mom ;
    for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ ){
      AliVCluster * calo = (AliVCluster *)(plNe->At(ipr)) ;
      
      //Get the index where the cluster comes, to retrieve the corresponding vertex
      Int_t evtIndex = 0 ; 
      if (reader->GetMixedEvent()) {
        evtIndex=reader->GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
      }
      
      //Do not count the candidate (photon or pi0) or the daughters of the candidate
      if(calo->GetID() == pCandidate->GetCaloLabel(0) || calo->GetID() == pCandidate->GetCaloLabel(1)) continue ;      //Skip matched clusters with tracks
      
      if(calo->GetNTracksMatched() > 0) continue ; 
      
      calo->GetMomentum(mom,reader->GetVertex(evtIndex)) ;//Assume that come from vertex in straight line
      
      pt   = mom.Pt();
      eta  = mom.Eta();
      phi  = mom.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      //only loop the particle at the same side of candidate
      
      if(TMath::Abs(phi-phiC)>TMath::PiOver2()) continue ;
      //if at the same side has particle larger than candidate, then candidate can not be the leading, skip such events
      if(pt > ptC){
        n         = -1;
        nfrac     = -1;
        coneptsum = -1;
        isolated  = kFALSE;
        if(bFillAOD){
          if(reftracks){  
            reftracks  ->Clear();
            delete reftracks;
          }
          if(refclusters){
            refclusters->Clear(); 
            delete refclusters;
          }
        }
        return ;
      }
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      if(rad < fConeSize){
        if(bFillAOD) {
          nclusterrefs++;
          if(nclusterrefs==1){
            refclusters = new TObjArray(0);
            //refclusters->SetName(Form("Clusters%s",aodArrayRefName.Data()));
            TString tempo(aodArrayRefName)  ; 
            tempo += "Clusters" ; 
            refclusters->SetName(tempo);
            refclusters->SetOwner(kFALSE);
          }
          refclusters->Add(calo);
        }
        //printf("neutral in isolation cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
        coneptsum+=pt;
        if(pt > fPtThreshold )     n++;
        //if fPtFraction*ptC<fPtThreshold then consider the fPtThreshold directly
        if(fPtFraction*ptC<fPtThreshold) {
          if(pt>fPtThreshold)    nfrac++ ;
        }
        else {
          if(pt>fPtFraction*ptC) nfrac++; 
        }
      }//in cone
    }// neutral particle loop
  }//neutrals
  
  //printf("Isolation Cut: in cone with: pT>pTthres %d, pT > pTfrac*pTcandidate %d \n",n,nfrac);
  
  //Add reference arrays to AOD when filling AODs only
  if(bFillAOD) {
    if(refclusters)	pCandidate->AddObjArray(refclusters);
    if(reftracks)	  pCandidate->AddObjArray(reftracks);
  }
  //Check isolation, depending on method.
  if( fICMethod == kPtThresIC){
    if(n==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtIC){
    if(coneptsum < fSumPtThreshold)
      isolated  =  kTRUE ;
  }
  else if( fICMethod == kPtFracIC){
    if(nfrac==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtFracIC){
    //when the fPtFraction*ptC < fSumPtThreshold then consider the later case
    if(fPtFraction*ptC < fSumPtThreshold  && coneptsum < fSumPtThreshold) isolated  =  kTRUE ;
    if(fPtFraction*ptC > fSumPtThreshold  && coneptsum < fPtFraction*ptC) isolated  =  kTRUE ;
  }
  
}

//_____________________________________________________
void AliIsolationCut::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s **** \n", GetName(), GetTitle() ) ;
  
  printf("IC method          =     %d\n",    fICMethod   ) ; 
  printf("Cone Size          =     %1.2f\n", fConeSize   ) ; 
  printf("pT threshold       =     %2.1f\n", fPtThreshold) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction ) ;
  printf("particle type in cone =  %d\n",    fPartInCone ) ;
  printf("    \n") ;
  
} 
