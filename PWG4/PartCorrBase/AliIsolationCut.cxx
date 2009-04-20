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
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
//#include <Riostream.h>
#include <TLorentzVector.h>
#include <TRefArray.h>

// --- AliRoot system --- 
#include "AliIsolationCut.h" 
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODTrack.h"
#include "AliAODCaloCluster.h"

ClassImp(AliIsolationCut)
  
//____________________________________________________________________________
  AliIsolationCut::AliIsolationCut() : 
    TObject(),
    fConeSize(0.),fPtThreshold(0.), fPtFraction(0.), fICMethod(0)
 
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliIsolationCut::AliIsolationCut(const AliIsolationCut & g) : 
  TObject(g),
  fConeSize(g.fConeSize),
  fPtThreshold(g.fPtThreshold),
  fPtFraction(g.fPtFraction), 
  fICMethod(g.fICMethod)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliIsolationCut & AliIsolationCut::operator = (const AliIsolationCut & source)
{
  // assignment operator
  
  if(&source == this) return *this;
  
  fConeSize = source.fConeSize ;
  fPtThreshold = source.fPtThreshold ; 
  fICMethod = source.fICMethod ;
  fPtFraction = source.fPtFraction ;

  return *this;
  
}

//____________________________________________________________________________
TString AliIsolationCut::GetICParametersList()
{
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  char onePar[255] ;
  
  sprintf(onePar,"--- AliIsolationCut ---\n") ;
  parList+=onePar ;	
  sprintf(onePar,"fConeSize: (isolation cone size) %1.2f\n",fConeSize) ;
  parList+=onePar ;
  sprintf(onePar,"fPtThreshold =%1.2f (isolation pt threshold) \n",fPtThreshold) ;
  parList+=onePar ;
  sprintf(onePar,"fPtFraction=%1.2f (isolation pt threshold fraction ) \n",fPtFraction) ;
  parList+=onePar ;
  sprintf(onePar,"fICMethod=%d (isolation cut case) \n",fICMethod) ;
  parList+=onePar ;
  
  return parList; 
}

//____________________________________________________________________________
void AliIsolationCut::InitParameters()
{
  //Initialize the parameters of the analysis.
  
  fConeSize             = 0.4 ; 
  fPtThreshold         = 1. ; 
  fPtFraction        = 0.1 ; 
  
  fICMethod = kPtThresIC; // 0 pt threshol method, 1 cone pt sum method
  
}

//__________________________________________________________________
void  AliIsolationCut::MakeIsolationCut(TRefArray * plCTS,  TRefArray * plNe, Double_t * vertex, 
					const Bool_t fillAOD, AliAODPWG4ParticleCorrelation  *pCandidate, 
					Int_t & n, Int_t & nfrac, Float_t &coneptsum,  Bool_t  &isolated) const
{  
  //Search in cone around a candidate particle if it is isolated 
  Float_t phiC  = pCandidate->Phi() ;
  Float_t etaC = pCandidate->Eta() ;
  Float_t ptC = pCandidate->Pt() ;
  Float_t pt     = -100. ;
  Float_t eta   = -100.  ;
  Float_t phi    = -100.  ;
  Float_t rad   = -100 ;
  Bool_t first   = kTRUE;
  n = 0 ;
  coneptsum = 0.; 
  isolated = kFALSE;
  
  //Check charged particles in cone.
  if(plCTS){
    TVector3 p3;
    for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ ){
      AliAODTrack* track = (AliAODTrack *)(plCTS->At(ipr)) ; 
      //Do not count the candidate (pion, conversion photon) or the daughters of the candidate
      if(track->GetID() == pCandidate->GetTrackLabel(0) || track->GetID() == pCandidate->GetTrackLabel(1)) continue ;
      p3.SetXYZ(track->Px(),track->Py(),track->Pz());
      pt   = p3.Pt();
      eta  = p3.Eta();
      phi  = p3.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      
      if(rad < fConeSize){
	if(fillAOD) {
	  if(first) {
	    new (pCandidate->GetRefIsolationConeTracks()) TRefArray(TProcessID::GetProcessWithUID(track)); 
	    first = kFALSE;
	  }
	  pCandidate->AddIsolationConeTrack(track);
	}
	//printf("charged in isolation cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
	if(pt > fPtFraction*ptC ) nfrac++;  
      }
    }// charged particle loop
  }//Tracks
  
  //Check neutral particles in cone.  
  if(plNe){
    first= kTRUE;
    TLorentzVector mom ;
    for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ ){
      AliAODCaloCluster * calo = (AliAODCaloCluster *)(plNe->At(ipr)) ;
      
      //Do not count the candidate (photon or pi0) or the daughters of the candidate
      if(calo->GetID() == pCandidate->GetCaloLabel(0) || calo->GetID() == pCandidate->GetCaloLabel(1)) continue ;      //Skip matched clusters with tracks
      
      if(calo->GetNTracksMatched() > 0) continue ; 
      
      calo->GetMomentum(mom,vertex);//Assume that come from vertex in straight line
      pt   = mom.Pt();
      eta  = mom.Eta();
      phi  = mom.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      if(rad < fConeSize){
	if(fillAOD) {
	  if(first) {
	    new (pCandidate->GetRefIsolationConeClusters()) TRefArray(TProcessID::GetProcessWithUID(calo)); 
	    first = kFALSE;
	  }
	  pCandidate->AddIsolationConeCluster(calo);
	}
	//printf("neutral in isolation cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
	if(pt > fPtFraction*ptC ) nfrac++;
      }//in cone
    }// neutral particle loop
  }//neutrals
  
  //printf("Isolation Cut: in cone with: pT>pTthres %d, pT > pTfrac*pTcandidate %d \n",n,nfrac);
  
  //Check isolation, depending on method.
  if( fICMethod == kPtThresIC){
    if(n==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtIC){
    if(coneptsum < fPtThreshold)
      isolated  =  kTRUE ;
  }
  else if( fICMethod == kPtFracIC){
    if(nfrac==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtFracIC){
    if(coneptsum < fPtFraction*ptC)
      isolated  =  kTRUE ;
  }
}

//__________________________________________________________________
void AliIsolationCut::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s **** \n", GetName(), GetTitle() ) ;
  
  printf("IC method          =     %d\n", fICMethod) ; 
  printf("Cone Size          =     %1.2f\n", fConeSize) ; 
  printf("pT threshold       =     %2.1f\n", fPtThreshold) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction) ;

  printf("    \n") ;
  
} 
