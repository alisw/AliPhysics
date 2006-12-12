//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: part of AliFlowEvent, allow flow study of v0s .
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowV0_h
#define AliFlowV0_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include "TMath.h"
#include "TObject.h"
#include "TNamed.h"

#include "AliFlowConstants.h"
#include "AliFlowTrack.h"

class Flow ;

class AliFlowV0 : public TNamed {

public:

                AliFlowV0();
                AliFlowV0(const Char_t* name);
  virtual       ~AliFlowV0();

 // Gets - V0's variables
  const Char_t* Pid()         const;
  Float_t       Phi()         const;
  Float_t       Eta()         const;
  Float_t       Pt()          const;
  Float_t       P()           const;
  Float_t       Y()           const;
  Short_t       Charge()      const;
  Float_t       Dca()         const;
  Float_t       Chi2()        const;
  Float_t       Mass()        const; 
  AliFlowTrack* DaughterP()   const; 
  AliFlowTrack* DaughterN()   const; 
  Float_t       V0Lenght()    const;
  Float_t       Sigma()       const;
  void          CrossPoint(Float_t Pxyz[3]) const;
  //TVector3      CrossPoint() const ;
  Float_t       CrossDca()    const;
  Int_t 	MostLikelihoodPID() const;
  Int_t         Label()        	const;
 // Sets - V0's variables
  void SetPid(const Char_t* pid);
  void SetPhi(Float_t phi); 	
  void SetEta(Float_t eta); 	
  void SetPt(Float_t pt); 	
  void SetChi2(Float_t chi2) ;	
  void SetVmass(Float_t mass);
  void SetCrossPoint(Float_t pox,Float_t poy,Float_t poz);
  void SetCrossDca(Float_t dca) ;
  void SetDca(Float_t dca) ; 
  void SetSigma(Float_t sigma) ;	  
  void SetDaughters(AliFlowTrack* pos, AliFlowTrack* neg);
  void SetMostLikelihoodPID(Int_t pdg_code);
  void SetLabel(Int_t label);

private:

 // Data Members
  Float_t  fPhi;				// reconstructed azimuthal angle of the v0
  Float_t  fPt;					// reconstructed transverse momentum of the v0 
  Float_t  fEta;				// reconstructed pseudorapidity of the v0
  Float_t  fChi2;				// chi2 of the reconstructed v0
  Float_t  fMass;    			       	// reconstructed v0 mass 
  Float_t  fDca;                             	// distance of closest approach of the reconstructed v0 to the main vertex
  Float_t  fCrossPoint[3] ;			// crossing point coordinates of the two daughter tracks
  Float_t  fCrossDCA ;				// DCA of the 2 tracks at the crossing point 
  Float_t  fSigma ;				//! DCA of the 2 tracks at the crossing point 
  Int_t    fMostLikelihoodPID;  	       	// most probable P.Id. hypotesis (...need to be implemented)
  Int_t    fLabel ;			       	// Label of the V0 (link: KineTree-ESD) 
  AliFlowTrack* fDaughters[2] ;		       	// daughter particles' pointers 
  
  ClassDef(AliFlowV0,1) ;                  	// macro for rootcint
};

#endif
//////////////////////////////////////////////////////////////////////
