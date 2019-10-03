//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowKineMaker.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//        AliFlowKineMaker provides the method to create AliFlowEvent(s)
// creates AliFlowEvent from the KineTree . 
// TParticle(s) is translated into AliFlowTrack(s), with exact momentum, 
// P.Id., etc. Very basic track cuts are applyed (like primaries). 
// The present class can be used in a simple AliRoot macro or in a 
// more complex enviroment such as AliSelector or AliTask.
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWKINEMAKER_CXX
#define ALIFLOWKINEMAKER_CXX

// ROOT things
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TTree.h>
#include "TClonesArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
//#include "TDatabasePDG.h"

// AliRoot things (...not used here, but in the macro)
//#include "AliRun.h"
//#include "AliRunLoader.h"
//#include "AliStack.h"

// Flow things
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include "AliFlowKineMaker.h"

// ANSI things
#include <stdlib.h>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowKineMaker) 
//-----------------------------------------------------------------------
AliFlowKineMaker::AliFlowKineMaker():
  fEventNumber(0), fPartNumber(0), fGoodTracks(0), fGoodV0s(0),
  fGoodTracksEta(0), fPosiTracks(0), fNegaTracks(0), fUnconstrained(0),
  fSumAll(0), fCutEvts(0), fCutParts(0), fNewAli(kFALSE), fLoopParts(kTRUE), fCounter(0), 
  fKTree(0x0), fParticle(0x0), fParticlePDG(0x0), fCharge(0),
  fRunID(0), fNumberOfEvents(0), fNumberOfParticles(0), fMagField(0), 
  fFlowEvent(0x0), fFlowTrack(0x0), fFlowV0(0x0), 
  fAbsEta(2.1), fElow(0.001), fEup(1000.), fPrimary(kTRUE)
{
 // default constructor 
 // resets counters , sets defaults
 
 for(Int_t bb=0;bb<5;bb++) { fBayesianAll[bb] = 0 ; }
 for(Int_t vv=0;vv<3;vv++) { fVertex[vv] = 0. ; }

 // particle cut
 fLabel[0] = 0 ; fLabel[1] = -1 ;

//  // TGeant3::AddParticlesToPdgDataBase() ---  Stolen From TGeant3.cxx ----(
//  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
//  const Int_t kion=10000000;
//  const Int_t kspe=50000000;
//  const Double_t kAu2Gev=0.9314943228;
//  const Double_t khSlash = 1.0545726663e-27;
//  const Double_t kErg2Gev = 1/1.6021773349e-3;
//  const Double_t khShGev = khSlash*kErg2Gev;
//  const Double_t kYear2Sec = 3600*24*365.25;
//  // Ions
//  pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,0,3,"Ion",kion+10020);
//  pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,khShGev/(12.33*kYear2Sec),3,"Ion",kion+10030);
//  pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,khShGev/(12.33*kYear2Sec),6,"Ion",kion+20040);
//  pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,0,6,"Ion",kion+20030);
//  // Special particles
//  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,0,0,"Special",kspe+50);
//  pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,0,0,"Special",kspe+51);
//  // ----------------------------------------------------------------------)
}
//-----------------------------------------------------------------------
AliFlowKineMaker::~AliFlowKineMaker()
{
 // default destructor (no actions) 
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AliFlowEvent* AliFlowKineMaker::FillFlowEvent(TTree* fKTree)
{
 // From the MC KineTree (input) fills the AliFlowEvent (output) . 
 // It loops on the stored TParticles and calls the methods to fill the 
 // arrays in the AliFlowEvent (charged -> tracks , neutral -> v0s) . 

 fFlowEvent = new AliFlowEvent(10000) ; if(!fFlowEvent) { return 0 ; }
 //cout << " -evt- " << fFlowEvent << endl ;

 fRunID = -1 ;
 fEventNumber = -1 ;
 fNumberOfParticles = fKTree->GetEntries() ; 
 //
 cout << " *evt n. " << fEventNumber << " (run " << fRunID << ")  -  tracks/v0s : " << fNumberOfParticles << endl ;

 // Event id 
 fFlowEvent->SetRunID(fRunID) ;  	       
 fFlowEvent->SetEventID(fEventNumber) ;         
 fFlowEvent->SetOrigMult((UInt_t)fNumberOfParticles) ;

 // Run information (fixed - ???)
 fMagField = 4 ; // (?)
 fFlowEvent->SetMagneticField(fMagField) ;	
 fFlowEvent->SetCenterOfMassEnergy(AliFlowConstants::fgCenterOfMassEnergy) ;  
 fFlowEvent->SetBeamMassNumberEast(AliFlowConstants::fgBeamMassNumberEast) ;  
 fFlowEvent->SetBeamMassNumberWest(AliFlowConstants::fgBeamMassNumberWest) ;  

 // Trigger information (now is: ULon64_t - some trigger mask)
 fFlowEvent->SetL0TriggerWord(-1); 
 
 // Get primary vertex position
 fVertex[0] = 0. ; fVertex[1] = 0. ; fVertex[2] = 0. ;         // fVertex = // ?! how to get primary vertex !?
 fFlowEvent->SetVertexPos(fVertex[0],fVertex[1],fVertex[2]) ; 

 // Zero Degree Calorimeter information
 Int_t zdcp = (Int_t)(TMath::Sqrt(TMath::Sqrt(fNumberOfParticles))) ;
 Float_t zdce[3] ; zdce[0] = -1 ; zdce[1] = -1 ; zdce[2] = -1 ;
 fFlowEvent->SetZDCpart(zdcp);  			
 fFlowEvent->SetZDCenergy(zdce[0],zdce[1],zdce[2]);	

 fKTree->SetBranchAddress("Particles",&fParticle) ;

 // Track (& V0) loop
 if(fLoopParts)
 {
  Int_t badPart = 0 ;   
  for(fPartNumber=0;fPartNumber<fNumberOfParticles;fPartNumber++) 
  {
   fKTree->GetEntry(fPartNumber) ;
   if(CheckTrack(fParticle))
   {
    // fParticlePDG = fParticle->GetPDG() ; 
    fCharge = (Int_t)((fParticle->GetPDG()->Charge())/3) ; // cout << fCharge << endl ;	 
    // fCharge = (Int_t)(TMath::Sign(1,(fParticle->GetPdgCode()))) ;		 
  
    if(TMath::Abs(fCharge) > 0)
    {
     FillFlowTrack(fParticle) ;   
     fGoodTracks++ ;				       
    }
    else if(fCharge == 0)
    {
     FillFlowV0(fParticle) ;
     fGoodV0s++ ;			 
    }
   }
   else { badPart++ ; continue ; }
  }
  fCutParts += badPart ;
 } // cout << " -particle number- :  " << fPartNumber << endl ;

 // Evt setting stuff
 fFlowEvent->SetCentrality();	
 				
 fCounter++ ; 
 return fFlowEvent ;
}
//----------------------------------------------------------------------
AliFlowTrack* AliFlowKineMaker::FillFlowTrack(TParticle* fParticle)
{
 // From a charged TParticle (input) fills the AliFlowTrack (output) .

 TString name = "" ; name += fPartNumber ;
 Int_t idx = fFlowEvent->TrackCollection()->GetEntries() ;
 fFlowTrack = (AliFlowTrack*)(fFlowEvent->TrackCollection()->New(idx)) ;
 fFlowTrack->SetName(name.Data()) ;
 
 // cout << " -tr- " << name.Data() << "(" << idx << ")"  << endl ;

 // TParticle label (link: KineTree-ESD)
 Int_t label = TMath::Abs(fPartNumber);
 fFlowTrack->SetLabel(label) ; 			

 // signed DCA from ESDtrack
 Float_t x = fParticle->Vx() ; 
 Float_t y = fParticle->Vy() ; 
 Float_t z = fParticle->Vz() ; 
 Float_t xy = TMath::Sqrt(x*x + y*y) ; 
 fFlowTrack->SetDcaSigned(xy,z) ;     // cout << "DCA :   xy = " << xy << "   z = " << z << endl ;	    

 // error on the DCA = 0
 fFlowTrack->SetDcaError(0.,0.,0.) ; 		    

 // UnConstrained (global) first
 Double_t gD[3] ; 				
 gD[0] = fParticle->Px() ; gD[1] = fParticle->Py() ; gD[2] = fParticle->Pz() ;			
 // -
 Float_t phiGl = (Float_t)Phi(gD) ;  
 if(phiGl<0) { phiGl += 2*TMath::Pi() ; }
 fFlowTrack->SetPhiGlobal(phiGl) ;		
 Float_t ptGl = (Float_t)Pt(gD) ;  if(ptGl<=0) { cout << " !!! ptGlobal = " << ptGl << endl ; }
 fFlowTrack->SetPtGlobal(ptGl) ;		
 Float_t etaGl = (Float_t)Eta(gD) ; 
 fFlowTrack->SetEtaGlobal(etaGl) ;		

 // Constrained (same, if primary)
 if((fParticle->IsPrimary()) && (Norm(gD)!=0.))  
 {  							
  fFlowTrack->SetPhi(phiGl) ;                 		
  fFlowTrack->SetPt(ptGl) ;          			
  fFlowTrack->SetEta(etaGl) ; 				

  // number of constrainable tracks with |eta| < AliFlowConstants::fgEtaGood (0.9)
  if(TMath::Abs(etaGl) < AliFlowConstants::fgEtaMid)  { fGoodTracksEta++ ; }
 }
 else  // in case Constriction impossible for track, fill the UnConstrained (global)
 {
  fUnconstrained++ ; 	
  fFlowTrack->SetPhi(0.) ;
  fFlowTrack->SetPt(0.) ;   
  fFlowTrack->SetEta(0.) ; 
 }	     

 // positive - negative tracks
 
 //Int_t fCharge = (Int_t)(fParticle->GetPDG()->Charge()) ; 
 fFlowTrack->SetCharge(fCharge) ;		
 if(fCharge>0) 		{ fPosiTracks++ ; }
 else if(fCharge<0) 	{ fNegaTracks++ ; }
 else 			{ return 0 ; }

 // Track parametrization (p at, hits, clusters, dE/dx) 
 Double_t pVecAt[3] ; 
 for(Int_t gg=0;gg<3;gg++) { pVecAt[gg] = gD[gg] ; }
 Float_t pAt = (Float_t)Norm(pVecAt) ;
 
 Int_t nClus[9] ; Int_t nHits[9] ; 
 nClus[0] = 1 ;   nHits[0] = 1 ; 			    // ITS - pixel
 nClus[1] = 1 ;   nHits[1] = 1 ;  
 nClus[2] = 1 ;   nHits[2] = 1 ; 			    // ITS - drift
 nClus[3] = 1 ;   nHits[3] = 1 ;  
 nClus[4] = 1 ;   nHits[4] = 1 ; 			    // ITS - strips
 nClus[5] = 1 ;   nHits[5] = 1 ;  
 nClus[6] = 160 ; nHits[6] = 1 ; 			    // TPC
 nClus[7] = 130 ; nHits[7] = 1 ; 			    // TRD
 nClus[8] = 1 ;   nHits[8] = 1 ; 			    // TOF

 Int_t pdgcode = 0 ;
 Float_t dEdx[4] ; for(Int_t de=0;de<4;de++) { dEdx[de] = -1. ; }
 Float_t detResFun[6] ; for(Int_t de=0;de<6;de++) { detResFun[de] = 0. ; }
 Float_t zFirst = fVertex[2] ; 
 Float_t zLast  = 1000. ;
 Float_t rFirst = TMath::Sqrt((fVertex[0]*fVertex[0]) + (fVertex[1]*fVertex[1])) ; 
 Float_t rLast  = 1000. ;

// // Geometrical acceptance (calculated assuming straight tracks) 
//
//  Float_t rDet[9][2] ; Float_t zDet[9] ;      
//  rDet[0][0] = 3.9 ;   rDet[0][1] = 3.9 ;   zDet[0] = 14.1/2. ;     // ITS - pixel
//  rDet[1][0] = 7.6 ;   rDet[1][1] = 7.6 ;   zDet[1] = 14.1/2. ;    
//  rDet[2][0] = 15.0 ;  rDet[2][1] = 15.0 ;  zDet[2] = 22.2/2. ;     // ITS - drift
//  rDet[3][0] = 23.9 ;  rDet[3][1] = 23.9 ;  zDet[3] = 29.7/2. ;    
//  rDet[4][0] = 37.8 ;  rDet[4][1] = 38.4 ;  zDet[4] = 43.1/2. ;     // ITS - strips
//  rDet[5][0] = 42.8 ;  rDet[5][1] = 43.4 ;  zDet[5] = 48.9/2. ;    
//  rDet[6][0] = 84.5 ;  rDet[6][1] = 246.6 ; zDet[6] = 500./2. ;     // TPC
//  rDet[7][0] = 290. ;  rDet[7][1] = 370. ;  zDet[7] = 700./2. ;     // TRD
//  rDet[8][0] = 370. ;  rDet[8][1] = 399. ;  zDet[8] = 745./2. ;     // TOF
//
//  Float_t Atheta = fParticle->Pz()/fParticle->Pt() ; 
//  Float_t z0 = fParticle->Vz() ; 
//  Float_t r0 = TMath::Sqrt((fParticle->Vx()*fParticle->Vx())+(fParticle->Vy()*fParticle->Vy())) ;
//  if((fParticle->Vx()*fParticle->Px()+fParticle->Vy()*fParticle->Py())>0) 	 { r0 *= 1. ; }   // sign given basing on track direction in respect to position
//  else if((fParticle->Vx()*fParticle->Px()+fParticle->Vy()*fParticle->Py())<0) { r0 *= -1.; }
//  else 						 			 { r0  = 0. ; }
//
//  // rFirst = rDet[0][0] ; rLast  = rDet[0][0] ;
//  zFirst = z0 + Atheta * (rDet[0][0] - r0) ; 
//  zLast  = z0 + Atheta * (rDet[4][1] - r0) ;
//  Float_t Pin  = 0. ; Float_t Pout = 0. ; Float_t Rout = 0. ; 
//  for(int dd=0;dd<9;dd++)
//  {
//   Pin  = z0 + Atheta * (rDet[dd][0] - r0) ; 
//   if(Pin<zDet[dd]) 
//   {
//    Pout = z0 + Atheta * (rDet[dd][1] - r0) ; 
//    if(TMath::Abs(Pout<zDet[dd]))  			// track gets in and out inside acceptance -> full hits
//    { 
//     nHits[dd] = nClus[dd] ; 
//     Rout = rDet[dd][1] ; 
//     rLast = TMath::Abs(Rout) ;		
//    } 
//    else 						// track goes out from one side -> SOME hits (...)
//    {
//     Rout = r0 + ((TMath::Sign(zDet[dd],eta)-z0)/Atheta) ; 
//     rLast = TMath::Abs(Rout) ;	
//     Float_t proportion = TMath::Abs((rLast-rDet[dd][0])/(rDet[dd][1]-rDet[dd][0])) ; proportion *= nClus[dd] ;	
//     nHits[dd] = (Int_t)proportion ; 
//    }   			
//   }
//   else 		 				// track does not get in -> zero hits
//   {
//    nHits[0] = 0 ; rFirst = 0. ; //rLast = 0. ; 
//   } 	
//  }
//
//  if(nHits[7]) 		{ dEdx[0] = 1. ; }  // implement bethe-block for TPC
//  if(nHits[5] || nHits[6]) 	{ dEdx[1] = 1. ; }  // implement bethe-block for ITS
//  if(nHits[8]) 		{ dEdx[2] = 1. ; }  // implement transition-radiation for TRD
//  if(nHits[9]) 		{ dEdx[3] = 1. ; }  // implement time of flight for TOF
//

 // P.id. (basing on the pdg code from MC -> an exact P.Id (probability=1) is given for [e , mu , pi , K , p , d], others are 0) 
 pdgcode = fParticle->GetPdgCode() ;			     // cout << FlowDebug << "PDG code = " << pdgcode << endl ; 
 if(TMath::Abs(pdgcode) == 11)  	     { detResFun[0] = 1. ; fBayesianAll[0]++ ; }
 else if(TMath::Abs(pdgcode) == 13)	     { detResFun[1] = 1. ; fBayesianAll[1]++ ; }
 else if(TMath::Abs(pdgcode) == 211)	     { detResFun[2] = 1. ; fBayesianAll[2]++ ; } 
 else if(TMath::Abs(pdgcode) == 321)	     { detResFun[3] = 1. ; fBayesianAll[3]++ ; } 
 else if(TMath::Abs(pdgcode) == 2212)	     { detResFun[4] = 1. ; fBayesianAll[4]++ ; }
 else if(TMath::Abs(pdgcode) == 10010020)    { detResFun[5] = 1. ; fBayesianAll[5]++ ; }
 else 		{ for(Int_t de=0;de<6;de++)  { detResFun[de] = 0.2 ; } }
 fSumAll++ ;

 // Fill the (fake) track parameters (fit , P.Id. ...)
 fFlowTrack->SetMostLikelihoodPID(pdgcode); 

 fFlowTrack->SetElectronPositronProb(detResFun[0]);		  
 fFlowTrack->SetMuonPlusMinusProb(detResFun[1]);
 fFlowTrack->SetPionPlusMinusProb(detResFun[2]);
 fFlowTrack->SetKaonPlusMinusProb(detResFun[3]);
 fFlowTrack->SetProtonPbarProb(detResFun[4]);
 fFlowTrack->SetDeuteriumAntiDeuteriumProb(detResFun[5]);    // *!* implement P.Id. for Deuterium

 fFlowTrack->SetZFirstPoint(zFirst) ;  fFlowTrack->SetZLastPoint(zLast) ;	     
 fFlowTrack->SetTrackLength(TMath::Sqrt((zLast-zFirst)*(zLast-zFirst)+(rLast-rFirst)*(rLast-rFirst))) ;
 fFlowTrack->SetChi2(0.) ;

 // Fill the (fake) detector information (nHits, dE/dx, det, resp. func.) 
 fFlowTrack->SetFitPtsTPC(nHits[6]) ;		     // cout << FlowDebug << "nHits TPC = " << nHits[6] << endl ; 
 fFlowTrack->SetMaxPtsTPC(nClus[6]) ;		     // cout << FlowDebug << "nClus = " << nClus[6] << endl ;  
 fFlowTrack->SetChi2TPC(nHits[6]/nClus[6]) ;	     // cout << FlowDebug << "Chi2 = " << nHits[6]/nClus[6] << endl ;	     
 fFlowTrack->SetDedxTPC(dEdx[0]) ; 		     // cout << FlowDebug << "Dedx = " << dEdx << endl ; 
 fFlowTrack->SetPatTPC(pAt) ; 	 		     // cout << FlowDebug << "p = " << pAt << endl ; 				
 fFlowTrack->SetRespFunTPC(detResFun) ; 	     // cout << FlowDebug << "response function = " << detResFun << endl ;  
 // -
 Int_t nITShits = 0 ; for(int dd=0;dd<6;dd++) { nITShits += nHits[dd] ; } 
 Int_t nITSclus = 6 ; 
 fFlowTrack->SetFitPtsITS(nITShits) ;		     // cout << FlowDebug << "nHits ITS = " << nITShits << endl ;				    
 fFlowTrack->SetMaxPtsITS(nITSclus) ;		     // cout << FlowDebug << "nClus = " << nITSclus << endl ; 
 fFlowTrack->SetChi2ITS(nITShits/nITSclus) ;	     // cout << FlowDebug << "Chi2 = " << nITShits/nITSclus << endl ; 
 fFlowTrack->SetDedxITS(dEdx[1]) ; 		     // cout << FlowDebug << "Dedx = " << dEdx << endl ; 
 fFlowTrack->SetPatITS(pAt) ; 	 		     // cout << FlowDebug << "p = " << pAt << endl ;  					
 fFlowTrack->SetRespFunITS(detResFun) ; 	     // cout << FlowDebug << "response function = " << detResFun << endl ;  
 // -
 fFlowTrack->SetNhitsTRD(nHits[7]) ;		     // cout << FlowDebug << "nHits TOF = " << nHits[7] << endl ;						    
 fFlowTrack->SetMaxPtsTRD(nClus[7]) ;		     // cout << FlowDebug << "nClus = " << nClus[7] << endl ;  
 fFlowTrack->SetChi2TRD(nHits[7]/nClus[7]) ;	     // cout << FlowDebug << "Chi2 = " << nHits[7]/nClus[7] << endl ; 
 fFlowTrack->SetSigTRD(dEdx[2]) ;  		     // cout << FlowDebug << "Dedx = " << dEdx << endl ; 
 fFlowTrack->SetPatTRD(pAt) ; 	 		     // cout << FlowDebug << "p = " << pAt << endl ;  					
 fFlowTrack->SetRespFunTRD(detResFun) ; 	     // cout << FlowDebug << "response function = " << detResFun << endl ;  
 // -
 fFlowTrack->SetNhitsTOF(nHits[8]) ;		     // cout << FlowDebug << "nHits TOF = " << nHits[8] << endl ;						     
 fFlowTrack->SetMaxPtsTOF(nClus[8]) ;		     // cout << FlowDebug << "nClus = " << nClus[8] << endl ; 
 fFlowTrack->SetChi2TOF(nHits[8]/nClus[8]) ;	     // cout << FlowDebug << "Chi2 = " << nHits[8]/nClus[8] << endl ; 
 fFlowTrack->SetTofTOF(dEdx[3]) ;  		     // cout << FlowDebug << "Dedx = " << 1. << endl ; 
 fFlowTrack->SetPatTOF(pAt) ; 	 		     // cout << FlowDebug << "p = " << pAt << endl ;  					
 fFlowTrack->SetRespFunTOF(detResFun) ; 	     // cout << FlowDebug << "response function = " << detResFun << endl ;  
 
 return fFlowTrack ; 
}
//-----------------------------------------------------------------------
AliFlowV0* AliFlowKineMaker::FillFlowV0(TParticle* fParticle)
{
 // From a neutral TParticle (input) fills the AliFlowV0 (output) .

 TString name = "" ; name += fPartNumber ;
 Int_t idx = fFlowEvent->V0Collection()->GetEntries() ;
 fFlowV0 = (AliFlowV0*)(fFlowEvent->V0Collection()->New(idx)) ;
 fFlowV0->SetName(name.Data()) ;
 
 // cout << " -v0- " << name.Data() << "(" << idx << ")"  << endl ;

 // TParticle label (link: KineTree-ESD)
 Int_t label = TMath::Abs(fPartNumber);
 fFlowV0->SetLabel(label) ;

 // reconstructed position of the V0 
 Double_t xyz[3] ; 		
 xyz[0] = fParticle->Vx() ; 
 xyz[1] = fParticle->Vy() ; 
 xyz[2] = fParticle->Vz() ; 
 fFlowV0->SetCrossPoint(xyz[0],xyz[1],xyz[2]) ;

 // V0's impact parameter & error (chi2 , DCA , sigma , pointing angle)  
 fFlowV0->SetDca((Float_t)Norm(xyz)) ;
 fFlowV0->SetSigma(0.) ;
 fFlowV0->SetCosPointingAngle(1.) ;  
 fFlowV0->SetDaughtersDca(0.) ;
 fFlowV0->SetChi2(0.) ;

 // reconstructed momentum of the V0
 Double_t pxyz[3] ;
 pxyz[0] = fParticle->Px() ; pxyz[1] = fParticle->Py() ; pxyz[2] = fParticle->Pz() ;			
 Float_t phi = (Float_t)Phi(pxyz) ; if(phi<0) { phi += 2*TMath::Pi() ; }
 fFlowV0->SetPhi(phi) ;		
 Float_t pt = (Float_t)Pt(pxyz) ; 
 fFlowV0->SetPt(pt) ;		
 Float_t eta = (Float_t)Eta(pxyz) ; 
 fFlowV0->SetEta(eta) ; 	

 // P.id. 
 Int_t pdgCode = fParticle->GetPdgCode() ;
 fFlowV0->SetMostLikelihoodPID(pdgCode); 	

 // mass 
 fFlowV0->SetVmass((Float_t)fParticle->GetMass()) ; 
 
 // daughters (should be taken on puprose from the KineTree, and wrote into the flow event)
 //Int_t nDaughters = fParticle->GetNDaughters() ;
 //Int_t d1 = fParticle->GetDaughter(nDaughters-1) ;
 //Int_t d2 = fParticle->GetDaughter(nDaughters-2) ;
//
// AliFlowTrack* pos = (AliFlowTrack*)fFlowEvent->TrackCollection()->At(pN) ; 
// AliFlowTrack* neg = (AliFlowTrack*)fFlowEvent->TrackCollection()->At(nN) ; 
// fFlowV0->SetDaughters(pos,neg) ;
//
 // d1 + d2 ; // warning: statement with no effect :)

 return fFlowV0 ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowKineMaker::CheckTrack(TParticle* fParticle) const
{
 // applies track cuts (pE , eta , label)

 Float_t eta   = (Float_t)fParticle->Eta() ;
 Float_t pE    = (Float_t)fParticle->P() ;
 Int_t   label = -1 ;   // check how to assign this !
 Bool_t  prim  = fParticle->IsPrimary() ;
 
 if(fAbsEta && (eta > fAbsEta)) 					 { return kFALSE ; }
 if((fElow < fEup) && ((pE<fElow) || (pE>fEup)))			 { return kFALSE ; }
 if((fLabel[0] < fLabel[1]) && ((label<fLabel[0]) || (label>fLabel[1]))) { return kFALSE ; }
 if(fPrimary && !prim) 							 { return kFALSE ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowKineMaker::CheckEvent(TTree* fKTree) const
{
 // applies event cuts (dummy)
 
 if(!fKTree) { return kFALSE ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
void AliFlowKineMaker::PrintCutList()
{ 
 // Prints the list of Cuts 

 cout << " * ESD cuts list * " << endl ; 
 if(fLabel[0]<fLabel[1])
 { 
  cout << "  Track Label [ " << fLabel[0] << " , " << fLabel[1] << " ] " << endl ; 
 }
 if(fAbsEta)
 { 
  cout << "  |eta| < " << fAbsEta << endl ; 
 }
 if(fElow<fEup)
 { 
  cout << "  Track Energy (P_total) [ " << fElow << " , " << fEup << " ] " << endl ; 
 }
 if(fPrimary) 							 
 { 
  cout << "  Primary Particles " << endl ;
 }
 cout << " *               * " << endl ;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//*** USEFULL METHODS for a 3-array of double (~ TVector3) ***
//-----------------------------------------------------------------------
Double_t AliFlowKineMaker::Norm(Double_t nu[3])
{ 
 // returns the norm of a double[3] 

 Double_t norm2 = nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2] ;
 return TMath::Sqrt(norm2) ; 
}
//-----------------------------------------------------------------------
Double_t AliFlowKineMaker::Phi(Double_t nu[3])
{
 // returns the azimuthal angle of a double[3] 

 if(nu[0]==0 && nu[1]==0) { return 0. ; }
 else 			  { return TMath::ATan2(nu[1],nu[0]) ; }
}
//-----------------------------------------------------------------------
Double_t AliFlowKineMaker::Pt(Double_t nu[3])
{
 // returns the transvers momentum of a double[3] 

 Double_t trans = nu[0]*nu[0] + nu[1]*nu[1] ;
 return TMath::Sqrt(trans) ; 
}
//-----------------------------------------------------------------------
Double_t AliFlowKineMaker::Eta(Double_t nu[3])
{
 // returns the PseudoRapidity of a double[3] 
 // if transvers momentum = 0 --> returns +/- 1.000

 Double_t m = Norm(nu) ;
 if(nu[0]!=0 || nu[1]!=0) { return 0.5*TMath::Log((m+nu[2])/(m-nu[2])) ; }
 else     	 	  { return TMath::Sign((Double_t)1000.,nu[2]) ; }
}
//-----------------------------------------------------------------------

#endif

//////////////////////////////////////////////////////////////////////
// - one way to open the alice KineTree -
//////////////////////////////////////////////////////////////////////
//
//   TString fileName = "galice.root" ;
//   rl = AliRunLoader::Open(fileName.Data(),"MyEvent","read");
//   rl->LoadgAlice();
//   gAlice = rl->GetAliRun();
//   rl->LoadHeader();
//   rl->LoadKinematics();
//   fNumberOfEvents = rl->GetNumberOfEvents() ;
//
//   Int_t exitStatus = rl->GetEvent(getEv) ; if(exitStatus!=0) { return kFALSE ; }
//
//   TTree* pKTree = (TTree*)rl->TreeK();	  // Particles' TTree (KineTree)
//   AliStack* pStack = gAlice->Stack();  	  // Particles' Stack - "Label()" to get the number in the stack 
//
// // else if(rl)      // opens files one by one (unload and reload)
// // {
// //  rl->UnloadgAlice() ;
// //  rl->UnloadHeader() ;
// //  rl->UnloadKinematics() ;
// //  delete rl ; rl = 0 ;
// // }
//
//   fNumberOfParticles = pKTree->GetEntries() ;
//   nPart = pStack->GetNtrack() ;
//   cout << " Event : " << evtN << " :  particles : " << fNumberOfParticles << "  (stack: " << nPart << ") . " << endl ; }
//
//////////////////////////////////////////////////////////////////////

