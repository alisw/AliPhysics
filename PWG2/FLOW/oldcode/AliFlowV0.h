//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowV0.h 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: an array of V0s is part of AliFlowEvent, 
//  this class allow flow study neutral secundary vertices (V0s)
//  as they are reconstructer by AliRoot at the ESD level .
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWV0_H
#define ALIFLOWV0_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include "TMath.h"
#include "TNamed.h"

#include "AliFlowTrack.h"

class AliFlowV0 : public TNamed {

public:

                AliFlowV0();
                AliFlowV0(const Char_t* name);
  virtual       ~AliFlowV0();

 // Gets
  const Char_t* Pid()          const ;
  Float_t       P()            const ;   
  Float_t       Y()            const ;
  Float_t       Phi()          const { return fPhi; }  
  Float_t       Eta()          const { return fEta; }  
  Float_t       Pt()           const { return fPt; }   
  Short_t       Charge()       const { return 0 ; }
  Float_t       Dca()          const { return fDca ; }
  Float_t       Chi2()         const { return fChi2; }
  Float_t       Mass()         const { return fMass ; }
  Int_t         DaughterP()    const { return fDaughterP ; } 
  Int_t         DaughterN()    const { return fDaughterN ; } 
  Float_t       V0Lenght()     const { return TMath::Sqrt(fCrossPoint[0]*fCrossPoint[0] + fCrossPoint[1]*fCrossPoint[1] + fCrossPoint[2]*fCrossPoint[2]) ; }
  Float_t       Sigma()        const { return fSigma ; }
  Float_t       CrossPointX()  const { return fCrossPoint[0] ; }
  Float_t       CrossPointY()  const { return fCrossPoint[1] ; } 
  Float_t       CrossPointZ()  const { return fCrossPoint[2] ; }
  Float_t       DaughtersDca() const { return fCrossDCA ; }
  Int_t         Label()        const { return fLabel ; }
  Float_t       GetCosPointingAngle() 	    const { return fPointAngle ; }
  Int_t 	MostLikelihoodPID() 	    const { return fMostLikelihoodPID; } 
  void          CrossPoint(Float_t Pxyz[3]) const { for(Int_t ii=0;ii<3;ii++) { Pxyz[ii] = fCrossPoint[ii] ; } }

 // Sets
  void SetPhi(Float_t phi) 		{ fPhi = phi; }   
  void SetEta(Float_t eta) 		{ fEta = eta; }   
  void SetPt(Float_t pt)   		{ fPt = pt; }     
  void SetChi2(Float_t chi2)	   	{ fChi2 = chi2; }   
  void SetVmass(Float_t mass) 	   	{ fMass = mass ; }
  void SetDca(Float_t dca)              { fDca = dca; }
  void SetSigma(Float_t sigma)	  	{ fSigma = sigma ; }	  
  void SetLabel(Int_t label)		{ fLabel = label ; }
  void SetDaughtersDca(Float_t dca)     { fCrossDCA = dca ; }
  void SetCosPointingAngle(Float_t cos) { fPointAngle = cos ; }
  void SetMostLikelihoodPID(Int_t pdgCode) 		  { fMostLikelihoodPID = pdgCode ; }
  void SetCrossPoint(Float_t pox,Float_t poy,Float_t poz) { fCrossPoint[0] = pox ; fCrossPoint[1] = poy ; fCrossPoint[2] = poz ; }
  void SetDaughters(Int_t pos, Int_t neg) 		  { fDaughterP = pos ; fDaughterN = neg ; }


private:

 // to make the code checker happy
  AliFlowV0(const AliFlowV0 &flowV0) ; 		  // Copy Constructor (dummy)
  AliFlowV0 &operator=(const AliFlowV0 &flowV0) ; // Assignment Operator (dummy)

 // Data Members
  Float_t   	fPhi;	     	    		// reconstructed azimuthal angle of the v0
  Float_t   	fPt;	     	    		// reconstructed transverse momentum of the v0 
  Float_t   	fEta;	     	    		// reconstructed pseudorapidity of the v0
  Float_t   	fChi2;       	    		// chi2 of the reconstructed v0
  Float_t   	fMass;       	    		// reconstructed v0 mass 
  Float_t  	fDca;		    		// distance of closest approach of the reconstructed v0 to the main vertex
  Float_t   	fCrossPoint[3] ;    		// crossing point coordinates of the two daughter tracks
  Float_t   	fCrossDCA ;  	    		// DCA between the 2 daughter tracks at the crossing point 
  Float_t   	fSigma ;     	    		// sigma of the DCA of the 2 daughter tracks at the crossing point 
  Int_t     	fLabel ;     	    		// Label of the V0 (link: KineTree-ESD) 
  Int_t    	fMostLikelihoodPID; 		// most probable P.Id. hypotesis
  Float_t   	fPointAngle;	    		// cosine of the pointing angle
  Int_t         fDaughterP ;		       	// positive daughter track (position in the TracksCollection())
  Int_t         fDaughterN ;		       	// negative daughter track (position in the TracksCollection())
  
  ClassDef(AliFlowV0,2) ;                  	// macro for rootcint
};

#endif
//////////////////////////////////////////////////////////////////////
