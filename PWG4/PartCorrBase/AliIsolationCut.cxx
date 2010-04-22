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
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
//#include <Riostream.h>
#include <TLorentzVector.h>
#include <TObjArray.h>

// --- AliRoot system --- 
#include "AliIsolationCut.h" 
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODTrack.h"
#include "AliAODCaloCluster.h"
#include "AliCaloTrackReader.h"

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
void  AliIsolationCut::MakeIsolationCut(TObjArray * const plCTS,  TObjArray * const plNe, AliCaloTrackReader * const reader, 
					const Bool_t fillAOD, AliAODPWG4ParticleCorrelation  *pCandidate, 
					const TString aodArrayRefName,
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
  n = 0 ;
  coneptsum = 0.; 
  isolated = kFALSE;

  //Initialize the array with refrences
  TObjArray * refclusters = 0x0;
  TObjArray * reftracks    =0x0;
  Int_t ntrackrefs = 0;
  Int_t nclusterrefs = 0;
  
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
	  ntrackrefs++;
	  if(ntrackrefs==1){
	    reftracks    = new TObjArray(0);
	    reftracks->SetName(aodArrayRefName+"Tracks");
	    reftracks->SetOwner(kFALSE);
	  }
	  reftracks->Add(track);
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
	  
	//Get vertex for photon momentum calculation
	Double_t vertex[]  = {0,0,0} ; //vertex ;
	Double_t vertex2[] = {0,0,0} ; //vertex second AOD input ;
	if(!reader->GetDataType()== AliCaloTrackReader::kMC) 
	{
		reader->GetVertex(vertex);
		if(reader->GetSecondInputAODTree()) reader->GetSecondInputAODVertex(vertex2);
	}
    TLorentzVector mom ;
    for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ ){
      AliAODCaloCluster * calo = (AliAODCaloCluster *)(plNe->At(ipr)) ;
      
      //Do not count the candidate (photon or pi0) or the daughters of the candidate
      if(calo->GetID() == pCandidate->GetCaloLabel(0) || calo->GetID() == pCandidate->GetCaloLabel(1)) continue ;      //Skip matched clusters with tracks
      
      if(calo->GetNTracksMatched() > 0) continue ; 
      //Input from second AOD?
      Int_t input = 0;
      if     (pCandidate->GetDetector() == "EMCAL" && reader->GetAODEMCALNormalInputEntries() <= ipr) input = 1 ;
      else if(pCandidate->GetDetector() == "PHOS"  && reader->GetAODPHOSNormalInputEntries()  <= ipr) input = 1;
      
      //Get Momentum vector, 
      if     (input == 0) calo->GetMomentum(mom,vertex) ;//Assume that come from vertex in straight line
      else if(input == 1) calo->GetMomentum(mom,vertex2);//Assume that come from vertex in straight line  
      
      pt   = mom.Pt();
      eta  = mom.Eta();
      phi  = mom.Phi() ;
      if(phi<0) phi+=TMath::TwoPi();
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      rad = TMath::Sqrt((eta-etaC)*(eta-etaC)+ (phi-phiC)*(phi-phiC));
      if(rad < fConeSize){
	if(fillAOD) {
	  nclusterrefs++;
	  if(nclusterrefs==1){
	    refclusters    = new TObjArray(0);
	    refclusters->SetName(aodArrayRefName+"Clusters");
	    //refclusters->SetOwner(kFALSE);
	  }
	  refclusters->Add(calo);
	}
	//printf("neutral in isolation cone pt %f, phi %f, eta %f, R %f \n",pt,phi,eta,rad);
	coneptsum+=pt;
	if(pt > fPtThreshold ) n++;
	if(pt > fPtFraction*ptC ) nfrac++;
      }//in cone
    }// neutral particle loop
  }//neutrals
  
  //printf("Isolation Cut: in cone with: pT>pTthres %d, pT > pTfrac*pTcandidate %d \n",n,nfrac);
  
  //Add reference arrays to AOD when filling AODs only
  if(fillAOD) {
    if(refclusters)	pCandidate->AddObjArray(refclusters);
    if(reftracks)	pCandidate->AddObjArray(reftracks);
  }

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

  //  if(refclusters) delete refclusters;
  //if(reftracks)   delete reftracks;

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
