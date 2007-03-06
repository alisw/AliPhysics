//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: track optimized for flow study, part of of AliFlowEvent. 
//                                                   Adapted from STAR 
// Original Authors:                 Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWTRACK_H
#define ALIFLOWTRACK_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TMath.h"
#include "TNamed.h"
#include "TString.h"
#include "TVector.h"

#include "AliFlowConstants.h"

class AliFlowTrack : public TNamed {

public:

                AliFlowTrack();
                AliFlowTrack(const Char_t* name);
  virtual       ~AliFlowTrack();

 // Gets - Track variables
  const Char_t* Pid()         	     const ;
  Float_t       Mass()        	     const ;
  Float_t  	InvMass()     	     const ; 
  Float_t       P()           	     const ;
  Float_t       Y()           	     const ;
  Float_t       PGlobal()     	     const ;
  Float_t       YGlobal()     	     const ;
  Bool_t        IsConstrainable()    const ;
  Float_t       Dca()	             const ;
  Float_t       DcaError() 	     const ;
  void          DcaError3(Float_t err[3]) const ;
  Int_t         Charge()      	     const { return TMath::Sign(1,MostLikelihoodPID()) ; } 
  Float_t       Phi()         	     const { if(IsConstrainable()) { return fPhi ; } else { return PhiGlobal() ; } }       
  Float_t       Eta()         	     const { if(IsConstrainable()) { return fEta ; } else { return EtaGlobal() ; } }
  Float_t       Pt()          	     const { if(IsConstrainable()) { return fPt  ; } else { return PtGlobal()  ; } }
  Float_t       PhiGlobal()   	     const { return fPhiGlobal ; }
  Float_t       EtaGlobal()   	     const { return fEtaGlobal ; }
  Float_t       PtGlobal()    	     const { return fPtGlobal ; }
  Float_t       Chi2()        	     const { return fChi2 ; }  
  Float_t       TrackLength() 	     const { return fTrackLength ; }
  Float_t       ZFirstPoint() 	     const { return fZFirstPoint ; }
  Float_t       ZLastPoint()  	     const { return fZLastPoint ; }
  Int_t         Label()       	     const { return fLabel ; }
  Float_t       TransDcaSigned()     const { return fDcaSigned[0] ; }
  Float_t       TransDca()    	     const { return TMath::Abs(fDcaSigned[0]) ; }  
  Float_t       ZedDca()    	     const { return TMath::Abs(fDcaSigned[1]) ; }  
  Float_t       TransDcaError()      const { return TMath::Abs(fDcaError[0]) ; }
  Float_t       ZedDcaError()        const { return TMath::Abs(fDcaError[1]) ; }

 // Gets - P.id.
  Float_t       MostLikelihoodProb()     	const;
  Float_t       PidProb(Int_t nn)   		const; 			 // normalized detector response (using AliFlowConstants::fgBayesian[] )
  TVector 	PidProbs() 			const; 
  void          PidProbs(Float_t pidN[AliFlowConstants::kPid])    const; 
  void          RawPidProbs(Float_t pidV[AliFlowConstants::kPid]) const; // raw ESD detector responses
  Float_t       ElectronPositronProb()  	const { return PidProb(0) ; } 
  Float_t       MuonPlusMinusProb()     	const { return PidProb(1) ; }
  Float_t       PionPlusMinusProb()     	const { return PidProb(2) ; }
  Float_t       KaonPlusMinusProb()     	const { return PidProb(3) ; }
  Float_t       ProtonPbarProb()        	const { return PidProb(4) ; }
  Float_t       DeuteriumAntiDeuteriumProb() 	const { return PidProb(5) ; }
  Int_t         MostLikelihoodPID()     	const { return fMostLikelihoodPID ; }

 // Gets - Detector response
  Float_t       Chi2TPC()     const { return fFitChi2[0] ; }
  Float_t       Chi2ITS()     const { return fFitChi2[1] ; }
  Float_t       Chi2TRD()     const { return fFitChi2[2] ; }
  Float_t       Chi2TOF()     const { return fFitChi2[3] ; }
  Int_t         FitPtsTPC()   const { return fFitPts[0]  ; }
  Int_t         FitPtsITS()   const { return fFitPts[1]  ; }
  Int_t         NhitsTRD()    const { return fFitPts[2]  ; }
  Int_t         NhitsTOF()    const { return fFitPts[3]  ; }
  Int_t         MaxPtsTPC()   const { return fMaxPts[0]  ; }
  Int_t         MaxPtsITS()   const { return fMaxPts[1]  ; }
  Int_t         MaxPtsTRD()   const { return fMaxPts[2]  ; }
  Int_t         MaxPtsTOF()   const { return fMaxPts[3]  ; }
  Float_t       DedxTPC()     const { return fDedx[0]	 ; }
  Float_t       DedxITS()     const { return fDedx[1]	 ; }
  Float_t       SigTRD()      const { return fDedx[2]	 ; }
  Float_t       TofTOF()      const { return fDedx[3]	 ; }
  Float_t       PatTPC()      const { return fMom[0] ; }
  Float_t       PatITS()      const { return fMom[1] ; }
  Float_t       PatTRD()      const { return fMom[2] ; }
  Float_t       PatTOF()      const { return fMom[3] ; }

 // Sets - Track variables
  void SetPid(const Char_t* pid);
  void SetCharge(Int_t charge) 			  { fMostLikelihoodPID = TMath::Sign(TMath::Abs(fMostLikelihoodPID),charge) ; }
  void SetPhi(Float_t phi)			  { fPhi = phi; }
  void SetEta(Float_t eta)			  { fEta = eta; }
  void SetPt(Float_t pt)  			  { fPt = pt; }
  void SetPhiGlobal(Float_t phi)		  { fPhiGlobal = phi; }
  void SetEtaGlobal(Float_t eta)		  { fEtaGlobal = eta; }
  void SetPtGlobal(Float_t pt)  		  { fPtGlobal = pt; }
  void SetDcaSigned(Float_t xy, Float_t z = 0.)   { fDcaSigned[0] = xy ; fDcaSigned[1] = z ; }
  void SetChi2(Float_t chi2)			  { fChi2 = chi2; }
  void SetTrackLength(Float_t tl)		  { fTrackLength = tl; }
  void SetZFirstPoint(Float_t zFirst)		  { fZFirstPoint = zFirst; }
  void SetZLastPoint(Float_t zLast)		  { fZLastPoint = zLast; }
  void SetLabel(Int_t label)			  { fLabel = label ; }
  void SetDcaError(Float_t xye, Float_t ze, Float_t xyze = 0.) { fDcaError[0] = xye ; fDcaError[1] = ze ; fDcaError[2] = xyze ; }

 // Sets - P.id.
  void SetMostLikelihoodPID(Int_t val)  	  { fMostLikelihoodPID = val ; }
  void SetElectronPositronProb(Float_t val)	  { fPidProb[0] = TMath::Abs(val) ; } 
  void SetMuonPlusMinusProb(Float_t val)	  { fPidProb[1] = TMath::Abs(val) ; } 
  void SetPionPlusMinusProb(Float_t val)	  { fPidProb[2] = TMath::Abs(val) ; } 
  void SetKaonPlusMinusProb(Float_t val)	  { fPidProb[3] = TMath::Abs(val) ; } 
  void SetProtonPbarProb(Float_t val)		  { fPidProb[4] = TMath::Abs(val) ; } 
  void SetDeuteriumAntiDeuteriumProb(Float_t val) { fPidProb[5] = TMath::Abs(val) ; }

 // Sets - Detector response
  void SetChi2TPC(Float_t chi2)			  { fFitChi2[0] = chi2   ; }
  void SetChi2ITS(Float_t chi2)			  { fFitChi2[1] = chi2   ; }
  void SetChi2TRD(Float_t chi2)			  { fFitChi2[2] = chi2   ; }
  void SetChi2TOF(Float_t chi2)			  { fFitChi2[3] = chi2   ; }
  void SetFitPtsTPC(Int_t fitPts)		  { fFitPts[0]  = fitPts ; }
  void SetFitPtsITS(Int_t nhits)		  { fFitPts[1]  = nhits ; }
  void SetNhitsTRD(Int_t fitPts)		  { fFitPts[2]  = fitPts ; }
  void SetNhitsTOF(Int_t hits)		  	  { fFitPts[3]  = hits ; }
  void SetMaxPtsTPC(Int_t maxPts)		  { fMaxPts[0]  = maxPts ; }
  void SetMaxPtsITS(Int_t maxPts)		  { fMaxPts[1]  = maxPts ; }
  void SetMaxPtsTRD(Int_t maxPts)		  { fMaxPts[2]  = maxPts ; }
  void SetMaxPtsTOF(Int_t maxPts)		  { fMaxPts[3]  = maxPts ; }
  void SetDedxTPC(Float_t dedx)	   	          { fDedx[0]    = dedx   ; }
  void SetDedxITS(Float_t dedx)	   	      	  { fDedx[1]    = dedx   ; }
  void SetSigTRD(Float_t trd)	   	          { fDedx[2]    = trd   ; }
  void SetTofTOF(Float_t tof)	   	          { fDedx[3]    = tof   ; }
  void SetPatTPC(Float_t p)			  { fMom[0] = p ; }
  void SetPatITS(Float_t p)			  { fMom[1] = p ; }
  void SetPatTRD(Float_t p)			  { fMom[2] = p ; }
  void SetPatTOF(Float_t p)			  { fMom[3] = p ; }

 // Selection methods
  void	 SetSelect(Int_t harmonic, Int_t selection) 		      { fSelection[harmonic][selection] = kTRUE ; }
  void	 SetSubevent(Int_t harmonic, Int_t selection, Int_t subevent) { fSubevent[harmonic][selection] = subevent ; }
  void	 ResetSelection() ;
  void	 PrintSelection() const ;
  Bool_t Select(Int_t harmonic, Int_t selection, Int_t subevent= -1) const ;

 // TRICKY
  void SetConstrainable();
  void SetUnConstrainable();


private:

 // to make the code checker happy
  AliFlowTrack &operator=(const AliFlowTrack &flowTrack) ; // Assignment Operator (dummy)

 // Data Members
  Float_t  fPhi;				// azimuthal angle of the constrained track
  Float_t  fEta;				// pseudorapidity of the constrained track
  Float_t  fPt;					// transverse momentum of the constrained track 
  Float_t  fZFirstPoint;  		       	// Z position at beginning of TPC
  Float_t  fZLastPoint;   		       	// Z position at end of PHOS
  Float_t  fChi2 ;				// constrained chi2
  Float_t  fTrackLength; 			// lenght of the track
  Int_t    fMostLikelihoodPID;  	       	// PDG code of the most probable P.Id.
  Double_t fDcaSigned[2] ; 			// positive or negative tracks -> P changes differently including vertex ...
  Double_t fDcaError[3] ; 			// error over the DCA
  Float_t  fPhiGlobal;				// azimuthal angle of the unconstrained track
  Float_t  fEtaGlobal;				// pseudorapidity of the unconstrained track
  Float_t  fPtGlobal;				// transverse momentum of the unconstrained track 
  Int_t    fLabel ;			       	// Label of the track (link: KineTree-ESD) 
 // -
  Float_t  fPidProb[AliFlowConstants::kPid] ;	// Array of probabilities to be (e,mu,pi,k,p,d)
  Int_t    fFitPts[4] ; 			// Fit Points (in: TPC,ITS,TRD,TOF)
  Int_t    fMaxPts[4] ; 			// Foundable Clusters (in: TPC,ITS,TRD,TOF)
  Float_t  fFitChi2[4] ;			// Chi2 (in: TPC,ITS,TRD,TOF)
  Float_t  fDedx[4] ;		 		// dE/dx from TPC and ITS , TRD signal , time of flight (TOF)
  Float_t  fMom[4] ;		 		// Track momentum at the entrance of TPC,ITS,TRD,TOF     

 // Selection's array for R.P. & sub-event selection (not stored)  
  Bool_t   fSelection[AliFlowConstants::kHars][AliFlowConstants::kSels]; //! 
  Short_t  fSubevent[AliFlowConstants::kHars][AliFlowConstants::kSels];	 //! 

  ClassDef(AliFlowTrack,5) ;                  	// macro for rootcint
};

#endif
//////////////////////////////////////////////////////////////////////
