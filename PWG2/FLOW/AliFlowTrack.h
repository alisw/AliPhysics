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

#ifndef AliFlowTrack_h
#define AliFlowTrack_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include "TMath.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TVector.h"

#include "AliFlowConstants.h"

class Flow ;

class AliFlowTrack : public TNamed {

public:

                AliFlowTrack();
                AliFlowTrack(const Char_t* name);
  virtual       ~AliFlowTrack();

 // Gets - Track variables
  const Char_t* Pid()         const;
  Int_t         Charge()      const;
  Float_t       Mass()        const;
  Float_t  	InvMass()     const; 
  Float_t       Phi()         const;
  Float_t       Eta()         const;
  Float_t       Pt()          const;
  Float_t       P()           const;
  Float_t       Y()           const;
  Float_t       PhiGlobal()   const;
  Float_t       EtaGlobal()   const;
  Float_t       PtGlobal()    const;
  Float_t       PGlobal()     const;
  Float_t       YGlobal()     const;
  Float_t       Chi2()        const;
  Float_t       TrackLength() const;
  Float_t       ZFirstPoint() const;
  Float_t       ZLastPoint()  const;
  Bool_t        IsConstrainable()    const;
  Float_t       TransDcaSigned()     const;
  Float_t       TransDca()    const;	    
  Float_t       Dca()	      const;
  Float_t	Dca2(Float_t dca[2]) const;
  Float_t       Dca3(Float_t dca[3]) const;
  Int_t         Label() 	        	const;
 // Gets - P.id.
  Int_t         MostLikelihoodPID()     	const;
  Float_t       MostLikelihoodProb()     	const;
  Float_t       ElectronPositronProb()  	const;
  Float_t       MuonPlusMinusProb()     	const;
  Float_t       PionPlusMinusProb()     	const;
  Float_t       KaonPlusMinusProb()     	const;
  Float_t       ProtonPbarProb()        	const;
  Float_t       DeuteriumAntiDeuteriumProb() 	const;
  Float_t       PidProb(Int_t nn)   		const;          // normalized detector response (using Flow::fBayesian[] )
  TVector 	PidProbs() 			const; 
  void          PidProbs(Float_t pidN[Flow::nPid])    const; 
  void          RawPidProbs(Float_t pidV[Flow::nPid]) const; 	// raw ESD detector responses
 // Gets - Detector response
  Float_t       Chi2TPC()     const;
  Int_t         FitPtsTPC()   const;
  Int_t         MaxPtsTPC()   const;
  Float_t       DedxTPC()     const;
  Float_t       Chi2ITS()     const;
  Int_t         FitPtsITS()   const;
  Int_t         MaxPtsITS()   const;
  Float_t       DedxITS()     const;
  Float_t       Chi2TRD()     const;
  Int_t         NhitsTRD()    const;
  Int_t         MaxPtsTRD()   const;
  Float_t       SigTRD()      const;
  Float_t       Chi2TOF()     const;
  Int_t         NhitsTOF()    const;
  Int_t         MaxPtsTOF()   const;
  Float_t       TofTOF()      const;  
  Float_t       PatTPC()      const;           // new
  Float_t       PatITS()      const;           // new
  Float_t       PatTRD()      const;           // new
  Float_t       PatTOF()      const;           // new


 // Sets - Track variables
  void SetPid(const Char_t* pid);
  void SetCharge(Int_t charge);
  void SetPhi(Float_t phi);
  void SetEta(Float_t eta);
  void SetPt(Float_t pt);
  void SetPhiGlobal(Float_t phi);
  void SetEtaGlobal(Float_t eta);
  void SetPtGlobal(Float_t pt);
  void SetDcaSigned(Float_t xy, Float_t z = 0.);
  void SetTransDcaSigned(Float_t xy);
  void SetChi2(Float_t chi2);
  void SetTrackLength(Float_t tl);
  void SetZFirstPoint(Float_t zFirst);
  void SetZLastPoint(Float_t zLast);
  void SetLabel(Int_t label);
 // Sets - P.id.
  void SetMostLikelihoodPID(Int_t val); 
  void SetElectronPositronProb(Float_t val);
  void SetMuonPlusMinusProb(Float_t val);
  void SetPionPlusMinusProb(Float_t val);
  void SetKaonPlusMinusProb(Float_t val);
  void SetProtonPbarProb(Float_t val);
  void SetDeuteriumAntiDeuteriumProb(Float_t val);
 // Sets - Detector response
  void SetChi2TPC(Float_t chi2);
  void SetFitPtsTPC(Int_t fitPts);
  void SetMaxPtsTPC(Int_t maxPts);
  void SetDedxTPC(Float_t dedx);
  void SetChi2ITS(Float_t chi2);
  void SetFitPtsITS(Int_t nhits);
  void SetMaxPtsITS(Int_t maxPts);
  void SetDedxITS(Float_t dedx);
  void SetChi2TRD(Float_t chi2);
  void SetNhitsTRD(Int_t fitPts);
  void SetMaxPtsTRD(Int_t maxPts);
  void SetSigTRD(Float_t dedx);
  void SetChi2TOF(Float_t chi2);
  void SetNhitsTOF(Int_t fitPts);
  void SetMaxPtsTOF(Int_t maxPts);
  void SetTofTOF(Float_t dedx);
  void SetPatTPC(Float_t p);           // new
  void SetPatITS(Float_t p);           // new
  void SetPatTRD(Float_t p);           // new
  void SetPatTOF(Float_t p);           // new
 // TRICKY
  void SetConstrainable();
  void SetUnConstrainable();

 // Selection's methods
  Bool_t Select(Int_t harmonic, Int_t selection, Int_t subevent= -1) const;
  void	 SetSelect(Int_t harmonic, Int_t selection);
  void	 SetSubevent(Int_t harmonic, Int_t selection, Int_t subevent);
  void	 PrintSelection() ;
  void	 ResetSelection() ;


private:

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
  //Double_t fDcaError[2] ; 			// error over the DCA (next to come)
  Float_t  fPhiGlobal;				// azimuthal angle of the unconstrained track
  Float_t  fEtaGlobal;				// pseudorapidity of the unconstrained track
  Float_t  fPtGlobal;				// transverse momentum of the unconstrained track 
  Int_t    fLabel ;			       	// Label of the track (link: KineTree-ESD) 
 // -
  Float_t  fPidProb[Flow::nPid] ;		// Array of probability to be   (e,mu,pi,k,p,d)
  Int_t    fFitPts[4] ; 			// Array of Fit Points 		(in: TPC,ITS,TRD,TOF)
  Int_t    fMaxPts[4] ; 			// Array of Foundable Clusters  (in: TPC,ITS,TRD,TOF)
  Float_t  fFitChi2[4] ;			// Array of Chi2                (in: TPC,ITS,TRD,TOF)
  Float_t  fDedx[4] ;		 		// Array of dE/dx               (in: TPC,ITS,TRD,TOF)
  Float_t  fMom[4] ;		 		// Array of momentum at the entrance of TPC,ITS,TRD,TOF     

 // Selection's array for R.P. & sub-event selection (not stored)  
  Bool_t   fSelection[Flow::nHars][Flow::nSels]; //! 
  Short_t  fSubevent[Flow::nHars][Flow::nSels];	 //! 

  ClassDef(AliFlowTrack,4) ;                  	// macro for rootcint
};

#endif
//////////////////////////////////////////////////////////////////////
