//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowMaker.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//        AliFlowMaker provides the method to create AliFlowEvent(s) 
// from AliESD(s). Very basic track cuts are applyed.
// The present class can be used in a simple AliRoot macro or in a 
// more complex enviroment such as AliSelector or AliTask.
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWMAKER_CXX
#define ALIFLOWMAKER_CXX

// ROOT things
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include "TClonesArray.h"

// AliRoot things
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
//#include "AliKalmanTrack.h"

// Flow things
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include "AliFlowMaker.h"

// ANSI things
#include <stdlib.h>
#include <vector>

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowMaker) 
//-----------------------------------------------------------------------
AliFlowMaker::AliFlowMaker():
  fEventNumber(0), fTrackNumber(0), fV0Number(0), fGoodTracks(0), fGoodV0s(0),
  fGoodTracksEta(0), fPosiTracks(0), fNegaTracks(0), fUnconstrained(0),
  fSumAll(0), fCutEvts(0), fCutTrks(0), fCutV0s(0), fCounter(0), fMovedTr(0x0),
  fNewAli(kFALSE), fLoopTrks(kTRUE), fLoopV0s(kTRUE),
  fESD(0x0), fTrack(0x0), fV0(0x0), fVertex(0x0),
  fRunID(0), fNumberOfEvents(0), fNumberOfTracks(0), 
  fNumberOfV0s(0), fMagField(0), 
  fFlowEvent(0x0), fFlowTrack(0x0), fFlowV0(0x0),
  fNHits(0), fElow(0.001), fEup(1000.)
{
 // default constructor 
 // resets counters , sets defaults
 
 for(Int_t bb=0;bb<5;bb++) { fBayesianAll[bb] = 0 ; } ;

 // trak cuts
 fLabel[0] = 0 ; fLabel[1] = -1 ;
}
//-----------------------------------------------------------------------
AliFlowMaker::~AliFlowMaker()
{
 // default destructor (no actions) 
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AliFlowEvent* AliFlowMaker::FillFlowEvent(AliESD* fESD)
{
 // From the AliESD (input) fills the AliFlowEvent (output) . 
 // It loops on track & v0 and calls the methods to fill the arrays . 
 // ... . 

 fFlowEvent = new AliFlowEvent(10000) ; if(!fFlowEvent) { return 0 ; }
 //cout << " -evt- " << fFlowEvent << endl ;

 fRunID = fESD->GetRunNumber() ;
 fEventNumber = -1 ;
 // fEventNumber = fESD->GetEventNumber() ;
 fNumberOfTracks = fESD->GetNumberOfTracks() ;
 fNumberOfV0s = fESD->GetNumberOfV0s() ;
 //
 cout << " *evt n. " << fEventNumber << " (run " << fRunID << ")  -  tracks: " << fNumberOfTracks << " ,   v0s " << fNumberOfV0s << endl ;
 
 // Clean the vector of links esd-flowEvt
 fMovedTr.clear(); // Remove all elements

 // Event id 
 fFlowEvent->SetRunID(fRunID) ;  	       
 fFlowEvent->SetEventID(fEventNumber) ;         
 fFlowEvent->SetOrigMult((UInt_t)fNumberOfTracks) ;

 // Run information (fixed - ???)
 fMagField = fESD->GetMagneticField() ; // cout << " *fMagField " << fMagField << endl ;
 fFlowEvent->SetMagneticField(fMagField) ;	
 fFlowEvent->SetCenterOfMassEnergy(AliFlowConstants::fgCenterOfMassEnergy) ;  
 fFlowEvent->SetBeamMassNumberEast(AliFlowConstants::fgBeamMassNumberEast) ;  
 fFlowEvent->SetBeamMassNumberWest(AliFlowConstants::fgBeamMassNumberWest) ;  

 // Trigger information (now is: ULon64_t - some trigger mask)
 fFlowEvent->SetL0TriggerWord((Int_t)fESD->GetTriggerMask()); 
 
 // Get primary vertex position
 fVertex = (AliESDVertex*)fESD->GetVertex() ;
 Double_t position[3] ; 
 fVertex->GetXYZ(position) ;
 fFlowEvent->SetVertexPos((Float_t)position[0],(Float_t)position[1],(Float_t)position[2]) ; 

 // Zero Degree Calorimeter information
 Int_t zdcp = fESD->GetZDCParticipants() ; 
 Float_t zdce[3] ; 
 zdce[0] = fESD->GetZDCN1Energy() + fESD->GetZDCN2Energy(); 
 zdce[1] = fESD->GetZDCP1Energy() + fESD->GetZDCP2Energy() ; 
 zdce[2] = fESD->GetZDCEMEnergy() ;
 fFlowEvent->SetZDCpart(zdcp);  			
 fFlowEvent->SetZDCenergy(zdce[0],zdce[1],zdce[2]);	

 // Track loop
 if(fLoopTrks)
 {
  Int_t badTrks = 0 ;   
  for(fTrackNumber=0;fTrackNumber<fNumberOfTracks;fTrackNumber++) 
  {
   fTrack = fESD->GetTrack(fTrackNumber) ;
   if(CheckTrack(fTrack))
   {
    FillFlowTrack(fTrack) ;   
    fGoodTracks++ ;				      
    fMovedTr.push_back(fFlowEvent->TrackCollection()->GetLast()) ; 
   }
   else { fMovedTr.push_back(-1) ; badTrks++ ; continue ; }
  }
  fCutTrks += badTrks ;
 } // cout << " -track number- :  " << fTrackNumber << endl ;

 // V0 loop
 if(fLoopV0s)
 {
  Int_t badV0s = 0 ;
  for(fV0Number=0;fV0Number<fNumberOfV0s;fV0Number++) 
  {
   fV0 = fESD->GetV0(fV0Number) ;			
   if(CheckV0(fV0))
   {
    FillFlowV0(fV0) ;
    fGoodV0s++ ; 			
   }
   else { badV0s++ ; continue ; }
  }
  fCutV0s += badV0s ;
 } // cout << " -v0 number- :  " << fV0Number << endl ;

 // Evt setting stuff
 fFlowEvent->SetCentrality();	
 				
 fCounter++ ; 
 return fFlowEvent ;
}
//----------------------------------------------------------------------
AliFlowTrack* AliFlowMaker::FillFlowTrack(AliESDtrack* fTrack)
{
 // From the AliESDtrack (input) fills the AliFlowTrack (output) .

 TString name = "" ; name += fTrackNumber ;
 Int_t idx = fFlowEvent->TrackCollection()->GetEntries() ;
 fFlowTrack = (AliFlowTrack*)(fFlowEvent->TrackCollection()->New(idx)) ;
 fFlowTrack->SetName(name.Data()) ;
 
 // cout << " -tr- " << name.Data() << "(" << idx << ")"  << endl ;

 // ESD particle label (link: KineTree-ESD)
 Int_t label = TMath::Abs(fTrack->GetLabel());
 fFlowTrack->SetLabel(label) ; 			

 // signed DCA from ESDtrack
 Float_t xy = 0 ; Float_t z = 0 ; 
 fTrack->GetImpactParameters(xy,z) ; 
 fFlowTrack->SetDcaSigned(xy,z) ; 		    

 // error on the DCA
 Float_t dcaBis[2] ; Float_t dcaCov[3] ; 
 for(Int_t dd=0;dd<3;dd++) { dcaCov[dd] = 0. ; }
 fTrack->GetImpactParameters(dcaBis, dcaCov) ;
 fFlowTrack->SetDcaError(dcaCov[0],dcaCov[1],dcaCov[2]) ; 		    

 // UnConstrained (global) first
 Double_t gD[3] ; 				
 fTrack->GetPxPyPz(gD) ;			
 // -
 Float_t phiGl = (Float_t)Phi(gD) ;  
 if(phiGl<0) { phiGl += 2*TMath::Pi() ; }
 fFlowTrack->SetPhiGlobal(phiGl) ;		
 Float_t ptGl = (Float_t)Pt(gD) ;  if(ptGl<=0) { cout << " !!! ptGlobal = " << ptGl << endl ; }
 fFlowTrack->SetPtGlobal(ptGl) ;		
 Float_t etaGl = (Float_t)Eta(gD) ; 
 fFlowTrack->SetEtaGlobal(etaGl) ;		

 // Constrained (NEW)
 Double_t cD[3] ;
 Double_t par1 ; Double_t par2 ;  Double_t par3[3] ;
 if(fTrack->GetConstrainedExternalParameters(par1,par2,par3))
 {
  fTrack->GetConstrainedPxPyPz(cD) ;     
 }
 else { for(Int_t iii=0;iii<3;iii++) { cD[iii] =0 ; } }

 if(Norm(cD)!=0.)   // ConstrainedPxPyPz != 0 if ConstrainedChi2 < something ...
 {  							
  Float_t phi = (Float_t)Phi(cD) ; 
  if(phi<0) { phi += 2*TMath::Pi() ; }
  fFlowTrack->SetPhi(phi) ;                 		
  Float_t pt = (Float_t)Pt(cD) ;   if(pt<=0) { cout << " !!! pt = " << pt << endl ; }
  fFlowTrack->SetPt(pt) ;          			
  Float_t eta = (Float_t)Eta(cD) ; 
  fFlowTrack->SetEta(eta) ; 				
 
  // number of constrainable tracks with |eta| < AliFlowConstants::fgEtaGood (0.9)
  if(TMath::Abs(eta) < AliFlowConstants::fgEtaMid)  { fGoodTracksEta++ ; }
 }
 else  // in case Constriction impossible for track, fill the UnConstrained (global)
 {
  fUnconstrained++ ; 	
  fFlowTrack->SetPhi(0.) ;
  fFlowTrack->SetPt(0.) ;   
  fFlowTrack->SetEta(0.) ; 
 }	     

 // positive - negative tracks
 Int_t trkSign = (Int_t)fTrack->GetSign() ; 
 fFlowTrack->SetCharge(trkSign) ;		
 if(trkSign>0) 	{ fPosiTracks++ ; }
 else if(trkSign<0) 	{ fNegaTracks++ ; }
 else 			{ return 0 ; }

 // Tracking parameters (fit , TPC , ITS , dE/dx)
 fFlowTrack->SetChi2(fTrack->GetConstrainedChi2()) ;		
 fFlowTrack->SetTrackLength(fTrack->GetIntegratedLength()) ;	
 // -
 Int_t idXt[180] ; // used for Cluster Map ( see AliESDtrack::GetTPCclusters() )    // old:    Int
 Int_t idX[6] ;    // used for Cluster Map ( see AliESDtrack::GetITSclusters() )    // old:    UInt
 Int_t idxr[130] ; // used for Cluster Map ( see AliESDtrack::GetTRDclusters() )    // old:    UInt
 Int_t nClus = 0 ;	
 Int_t fNFound = 0 ;  					// *!* fNFoundable (in AliTPCtrack) ... added by M.Ianov 
 // -
 Double_t detPid[5] ;		
 Float_t  detPid6[AliFlowConstants::kPid] ; 

 Double_t pVecAt[3] ; 
 for(Int_t gg=0;gg<3;gg++) { pVecAt[gg] = gD[gg] ; }
 Bool_t boh ; Float_t pAt = 0 ; 			// to get p at each detector
 // -
 if(fNewAli) { boh = fTrack->GetPxPyPzAt(AliFlowConstants::fgTPCx, fMagField, pVecAt) ; }
 else 	     { boh = fTrack->GetInnerParam()->GetPxPyPzAt(AliFlowConstants::fgTPCx, fMagField, pVecAt) ; }
 pAt = (Float_t)Norm(pVecAt) ; if(!pAt) { pAt = (Float_t)Norm(gD) ; }
 nClus = fTrack->GetTPCclusters(idXt) ;
 fNFound = fTrack->GetTPCNclsF() ;  // was 160
 if( (fTrack->GetStatus() & AliESDtrack::kTPCpid) != 0 )
 {
  fTrack->GetTPCpid(detPid) ; 
  for(Int_t bb=0;bb<5;bb++) { detPid6[bb] = detPid[bb] ; } detPid6[5] = 0. ; 
  fFlowTrack->SetRespFunTPC(detPid6) ; 
 }
 fFlowTrack->SetMaxPtsTPC(fNFound) ; 	 			
 fFlowTrack->SetFitPtsTPC(nClus) ;                		
 fFlowTrack->SetDedxTPC(fTrack->GetTPCsignal()) ; 		
 fFlowTrack->SetChi2TPC((Float_t)(fTrack->GetTPCchi2())) ; 	
 fFlowTrack->SetPatTPC(pAt) ; 					
 // -
 if(fNewAli) { boh = fTrack->GetPxPyPzAt(AliFlowConstants::fgITSx, fMagField, pVecAt) ; }
 else 	     { boh = fTrack->GetInnerParam()->GetPxPyPzAt(AliFlowConstants::fgITSx, fMagField, pVecAt) ; }
 pAt = (Float_t)Norm(pVecAt) ; if(!pAt) { pAt = (Float_t)Norm(gD) ; }
 nClus = fTrack->GetITSclusters(idX) ;
 fNFound = 6 ; // ? fixed
 if( (fTrack->GetStatus() & AliESDtrack::kITSpid) != 0 )
 {
  fTrack->GetITSpid(detPid) ; 
  for(Int_t bb=0;bb<5;bb++) { detPid6[bb] = detPid[bb] ; } detPid6[5] = 0. ; 
  fFlowTrack->SetRespFunITS(detPid6) ; 
 } 
 fFlowTrack->SetMaxPtsITS(fNFound) ; 	 			
 fFlowTrack->SetFitPtsITS(nClus) ;				
 fFlowTrack->SetDedxITS(fTrack->GetITSsignal()) ;		
 fFlowTrack->SetChi2ITS((Float_t)(fTrack->GetITSchi2())) ; 	
 fFlowTrack->SetPatITS(pAt) ; 					
 // -
 if(fNewAli) { boh = fTrack->GetPxPyPzAt(AliFlowConstants::fgTRDx, fMagField, pVecAt) ; } 
 else 	     { boh = fTrack->GetInnerParam()->GetPxPyPzAt(AliFlowConstants::fgTRDx, fMagField, pVecAt) ; } 
 pAt = (Float_t)Norm(pVecAt) ; if(!pAt) { pAt = (Float_t)Norm(gD) ; }
 nClus = fTrack->GetTRDclusters(idxr) ;
 fNFound = fTrack->GetTRDncls() ;  // was 130
 if( (fTrack->GetStatus() & AliESDtrack::kTRDpid) != 0 )
 {
  fTrack->GetTRDpid(detPid) ; 
  for(Int_t bb=0;bb<5;bb++) { detPid6[bb] = detPid[bb] ; } detPid6[5] = 0. ; 
  fFlowTrack->SetRespFunTRD(detPid6) ; 
 }
 fFlowTrack->SetMaxPtsTRD(fNFound) ;	 			
 fFlowTrack->SetNhitsTRD(nClus) ;				
 fFlowTrack->SetSigTRD(fTrack->GetTRDsignal()) ;		
 fFlowTrack->SetChi2TRD((Float_t)fTrack->GetTRDchi2()) ; 	
 fFlowTrack->SetPatTRD(pAt) ; 					
 // -
 if(fNewAli) { boh = fTrack->GetPxPyPzAt(AliFlowConstants::fgTOFx, fMagField, pVecAt) ; }
 else 	     { boh = fTrack->GetInnerParam()->GetPxPyPzAt(AliFlowConstants::fgTOFx, fMagField, pVecAt) ; }
 pAt = (Float_t)Norm(pVecAt) ; if(!pAt) { pAt = (Float_t)Norm(gD) ; }
 nClus = fTrack->GetTOFcluster() ;
 fNFound = 0 ; if(fTrack->GetTOFCalChannel() > 0) { fNFound = 1 ; }
 if( (fTrack->GetStatus() & AliESDtrack::kTOFpid) != 0 )
 {
  fTrack->GetTOFpid(detPid) ; 
  for(Int_t bb=0;bb<5;bb++) { detPid6[bb] = detPid[bb] ; } detPid6[5] = 0. ; 
  fFlowTrack->SetRespFunTOF(detPid6) ; 
 }
 fFlowTrack->SetMaxPtsTOF(fNFound) ;				
 fFlowTrack->SetNhitsTOF(nClus) ;				
 fFlowTrack->SetTofTOF(fTrack->GetTOFsignal()) ;		
 fFlowTrack->SetChi2TOF(fTrack->GetTOFchi2()) ; 		
 fFlowTrack->SetPatTOF(pAt) ; 					
 // -
 Double_t rIn[3]  ; rIn[0] = 0.  ;  rIn[1] = 0. ;  rIn[2] = 0. ;
 Double_t rOut[3] ; rOut[0] = 0. ; rOut[1] = 0. ; rOut[2] = 0. ;
 // -
 fTrack->GetInnerXYZ(rIn) ;			 		
 fFlowTrack->SetZFirstPoint(rIn[2]) ; 
 //fTrack->GetXYZAt(AliFlowConstants::fgTPCx,fMagField,rOut) ;		
 fTrack->GetOuterXYZ(rOut) ;			 		
 fFlowTrack->SetZLastPoint(rOut[2]) ; 
 
 // ESD-P.Id. = 5-vector of Best detectors probabilities for [e , mu , pi , K , p] 
 Double_t trkPid[5] ; fTrack->GetESDpid(trkPid) ;		
 Double_t trkPid6[AliFlowConstants::kPid] ; 
 for(Int_t bb=0;bb<5;bb++) { trkPid6[bb] = trkPid[bb] ; } 
 trkPid6[5] = 0. ;						// *!* implement P.Id. for Deuterim

 // Bayesian P.Id. method (weighted probabilities for [e , mu , pi , K , p , d])
 Double_t bsum = 0 ; 
 Double_t bayePid[AliFlowConstants::kPid] ;    // normalized P.id
 Double_t storedPid[AliFlowConstants::kPid] ;  // stored P.id
 for(Int_t nB=0;nB<AliFlowConstants::kPid;nB++)  { bsum += trkPid6[nB]*AliFlowConstants::fgBayesian[nB] ; }
 if(bsum)
 {
  for(Int_t nB=0;nB<AliFlowConstants::kPid;nB++) 
  { 
   bayePid[nB] = trkPid6[nB]*AliFlowConstants::fgBayesian[nB] / bsum ; 
   storedPid[nB] = trkPid6[nB] ; 
  }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
 
 fFlowTrack->SetElectronPositronProb(storedPid[0]);               
 fFlowTrack->SetMuonPlusMinusProb(storedPid[1]);
 fFlowTrack->SetPionPlusMinusProb(storedPid[2]);
 fFlowTrack->SetKaonPlusMinusProb(storedPid[3]);
 fFlowTrack->SetProtonPbarProb(storedPid[4]);
 fFlowTrack->SetDeuteriumAntiDeuteriumProb(storedPid[5]); 	// *!* implement P.Id. for Deuterim

 // P.id. label given via the weighted prob.
 const Int_t kCode[]   =  {11,13,211,321,2212,10010020} ;
 Int_t kkk = 2 ; 			// if No id. -> then is a Pi
 Float_t pidMax = bayePid[2] ; 	// (if all equal, Pi probability get's the advantage to be the first)
 for(Int_t iii=0; iii<5; iii++) 
 {
  if(bayePid[iii]>pidMax) { kkk = iii ; pidMax = bayePid[iii] ; }  // !!! Bayesian as well !!!
 }
 fBayesianAll[kkk]++ ; fSumAll++ ; 	// goes on filling the vector of observed abundance 
 //-
 Int_t pdgCode = trkSign*kCode[kkk] ;
 fFlowTrack->SetMostLikelihoodPID(pdgCode);
 
 return fFlowTrack ; 
}
//-----------------------------------------------------------------------
AliFlowV0* AliFlowMaker::FillFlowV0(AliESDv0* fV0)
{
 // From the AliESDv0 (input) fills the AliFlowV0 (output) .

 TString name = "" ; name += fV0Number ;
 Int_t idx = fFlowEvent->V0Collection()->GetEntries() ;
 fFlowV0 = (AliFlowV0*)(fFlowEvent->V0Collection()->New(idx)) ;
 fFlowV0->SetName(name.Data()) ;
 
 // cout << " -v0- " << name.Data() << "(" << idx << ")"  << endl ;

 // ESD particle label (link: KineTree-ESD)
 Int_t label = -1 ; // TMath::Abs(fV0->GetLabel());
 fFlowV0->SetLabel(label) ; 			

 // reconstructed momentum of the V0
 Double_t pxyz[3] ;
 fV0->GetPxPyPz(pxyz[0],pxyz[1],pxyz[2]) ;

 Float_t phi = (Float_t)Phi(pxyz) ; if(phi<0) { phi += 2*TMath::Pi() ; }
 fFlowV0->SetPhi(phi) ;		
 Float_t pt = (Float_t)Pt(pxyz) ; 
 fFlowV0->SetPt(pt) ;		
 Float_t eta = (Float_t)Eta(pxyz) ; 
 fFlowV0->SetEta(eta) ; 	

 // reconstructed position of the V0 
 Double_t xyz[3] ; 		
 fV0->GetXYZ(xyz[0],xyz[1],xyz[2]) ;	        
 fFlowV0->SetCrossPoint(xyz[0],xyz[1],xyz[2]) ;

 // V0's impact parameter & error (chi2 , DCA , sigma , pointing angle)  
 //fFlowV0->SetDca((Float_t)fV0->GetD()) ;    // GetDistNorm 
 //fFlowV0->SetSigma((Float_t)fV0->GetDistSigma()) ;
 //fFlowV0->SetCosPointingAngle((Float_t)fV0->GetV0CosineOfPointingAngle()) ;  
 //fFlowV0->SetDaughtersDca(fV0->GetDcaV0Daughters()) ;
 //fFlowV0->SetChi2((Float_t)fV0->GetChi2V0()) ;     // AliRoot v4-04-Release (December 2006)  
 //fFlowV0->SetChi2((Float_t)fV0->GetChi2()) ;       // AliRoot v4-04-Release (old)
 // ...when they'll stop changing the methods I'll enable the above lines. For now:
 fFlowV0->SetDca(0.);
 fFlowV0->SetSigma(0.1); 
 fFlowV0->SetCosPointingAngle(1.) ;	
 fFlowV0->SetDaughtersDca(0.) ;
 fFlowV0->SetChi2(1.) ;

 // P.id. 
 Int_t pdgCode = fV0->GetPdgCode() ;
 fFlowV0->SetMostLikelihoodPID(pdgCode); 	

 // mass 
 fFlowV0->SetVmass((Float_t)fV0->GetEffMass()) ; 
 
 // daughters 
 Int_t pN = fV0->GetPindex() ;
 Int_t nN = fV0->GetNindex() ;
 fFlowV0->SetDaughters(fMovedTr[pN],fMovedTr[nN]) ;

 return fFlowV0 ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowMaker::CheckTrack(AliESDtrack* fTrack) const
{
 // applies track cuts (pE , nHits , label)

 if(!fTrack) { return kFALSE ; }

 Int_t idXt[180] ;  // used for Cluster Map ( see AliESDtrack::GetTPCclusters() )
 Int_t   nHits = fTrack->GetTPCclusters(idXt) ;
 Float_t pE    = fTrack->GetP() ;
 Int_t   label = fTrack->GetLabel() ;
 
 if(fNHits && (nHits<=fNHits)) 						 { return kFALSE ; }
 if((fElow < fEup) && ((pE<fElow) || (pE>fEup)))			 { return kFALSE ; }
 if((fLabel[0] < fLabel[1]) && ((label<fLabel[0]) || (label>fLabel[1]))) { return kFALSE ; }
 //if(fPrimary && ...) 							 { return kFALSE ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowMaker::CheckV0(AliESDv0* fV0) const
{
 // applies v0 cuts (mass>fElow, daughters included in the tracks)
 
 if(!fV0) { return kFALSE ; }
 
 Int_t pN = fV0->GetPindex() ;
 Int_t nN = fV0->GetNindex() ;
 Float_t iMass = (Float_t)fV0->GetEffMass() ;
 
 if(iMass < fElow)      { return kFALSE ; }
 if(fMovedTr[pN] < 0)   { return kFALSE ; }
 if(fMovedTr[nN] < 0)   { return kFALSE ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowMaker::CheckEvent(AliESD* fESD) const
{
 // applies event cuts (dummy)
 
 if(!fESD) { return kFALSE ; }

 return kTRUE ;
}
//-----------------------------------------------------------------------
void AliFlowMaker::PrintCutList()
{ 
 // Prints the list of Cuts 

 cout << " * ESD cuts list * " << endl ; 
 if(fLabel[0]<fLabel[1])
 { 
  cout << "  Track Label [ " << fLabel[0] << " , " << fLabel[1] << " ] " << endl ; 
 }
 if(fNHits)
 { 
  cout << "  TPC clusters > " << fNHits << endl ; 
 }
 if(fElow<fEup)
 { 
  cout << "  Track Energy (P_total) [ " << fElow << " , " << fEup << " ] " << endl ; 
 }
 if(fElow>0.)
 { 
  cout << "  V0 invariant mass > " << fElow << endl ; 
 }
 cout << " *               * " << endl ;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//*** USEFULL METHODS for a 3-array of double (~ TVector3) ***
//-----------------------------------------------------------------------
Double_t AliFlowMaker::Norm(Double_t nu[3])
{ 
 // returns the norm of a double[3] 

 Double_t norm2 = nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2] ;
 return TMath::Sqrt(norm2) ; 
}
//-----------------------------------------------------------------------
Double_t AliFlowMaker::Phi(Double_t nu[3])
{
 // returns the azimuthal angle of a double[3] 

 if(nu[0]==0 && nu[1]==0) { return 0. ; }
 else 			  { return TMath::ATan2(nu[1],nu[0]) ; }
}
//-----------------------------------------------------------------------
Double_t AliFlowMaker::Pt(Double_t nu[3])
{
 // returns the transvers momentum of a double[3] 

 Double_t trans = nu[0]*nu[0] + nu[1]*nu[1] ;
 return TMath::Sqrt(trans) ; 
}
//-----------------------------------------------------------------------
Double_t AliFlowMaker::Eta(Double_t nu[3])
{
 // returns the PseudoRapidity of a double[3] 
 // if transvers momentum = 0 --> returns +/- 1.000

 Double_t m = Norm(nu) ;
 if(nu[0]!=0 || nu[1]!=0) { return 0.5*TMath::Log((m+nu[2])/(m-nu[2])) ; }
 else     	 	  { return TMath::Sign((Double_t)1000.,nu[2]) ; }
}
//-----------------------------------------------------------------------

#endif
