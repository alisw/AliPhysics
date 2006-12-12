//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: 
//         an array of AliFlowV0 is part of the AliFlowEvent, 
// The object AliFlowV0 contains data members wich summarize the V0s 
// information most useful for flow study (Pt, eta, phi, ecc.). 
// AliFlowV0 also contains a references to the two daughter tracks
// in the TrackClollection() of the AliFlowEvent from wich the V0 has 
// been reconstructed.
//

#include "AliFlowV0.h"
#include "AliFlowConstants.h"
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowV0) ;
//////////////////////////////////////////////////////////////////////////////
AliFlowV0::AliFlowV0()
{
 // default constructor  
 
 fPhi = 0. ;
 fEta = 0. ;
 fPt = 0. ;
 fChi2 = 0. ;
 fMass = 0. ;			
 fDca = 0. ;
 fCrossDCA = 0. ;
 fSigma = 1. ;
 fLabel = 0 ;
 for(Int_t dd=0;dd<3;dd++) { fCrossPoint[dd] = 0. ; }
 for(Int_t dd=0;dd<2;dd++) { fDaughters[dd] = 0  ; }
}
//////////////////////////////////////////////////////////////////////////////
AliFlowV0::AliFlowV0(const Char_t* name) 
{
 // TNamed constructor 
 
 SetName(name) ;
 AliFlowV0() ;
}
//////////////////////////////////////////////////////////////////////////////
AliFlowV0::~AliFlowV0() 
{
 // default destructor (dummy)
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowV0::P() const 
{ 
 // Returns the reconstructed momentum of the v0

 float momentum = Pt()/TMath::Sqrt(1-(tanh(Eta())*tanh(Eta()))) ; 
 return momentum; 
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowV0::Y() const 
{
 // Rapidity of the v0

 if(TMath::Abs((Float_t)P()) == TMath::Abs((Float_t)Pt())) 	{ return 0. ; }
 else if(TMath::Abs((Float_t)P()) < TMath::Abs((Float_t)Pt()))	{ cout << "v0: " << GetName() << "has  Pt() > P() !!!" << endl ; }
 // -
 float M = Mass() ; 
 double Pz = TMath::Sqrt(P()*P() - Pt()*Pt()); 
 if (Eta() < 0) { Pz = -Pz ; }
 double E = TMath::Sqrt(P()*P() + M*M) ;
 float rapidity = 0.5*TMath::Log((E + Pz)/(E - Pz)) ;
 return rapidity ;
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowV0::SetDaughters(AliFlowTrack* pos, AliFlowTrack* neg) 
{ 
 // Sets positive and negative daughter tracks from the TrackCollection()

 fDaughters[0] = pos ; 
 fDaughters[1] = neg ; 
}
//////////////////////////////////////////////////////////////////////////////
void  AliFlowV0::CrossPoint(Float_t Pxyz[3]) const 
{
 // Coordinates of the v0

 for(Int_t ii=0;ii<3;ii++) { Pxyz[ii] = fCrossPoint[ii] ; }
}
//////////////////////////////////////////////////////////////////////////////
// TVector3 AliFlowV0::CrossPoint() const 
// {
//  // Coordinates of the v0
// 
//  TVector3 pxyz ;
//  pxyz.SetXYZ(fCrossPoint[0],fCrossPoint[1],fCrossPoint[2]) ;
//  return pxyz ; 
// }
//////////////////////////////////////////////////////////////////////////////
void AliFlowV0::SetCrossPoint(Float_t pox,Float_t poy,Float_t poz)
{ 
 // Sets the coordinates of the v0

 fCrossPoint[0] = pox ; fCrossPoint[1] = poy ; fCrossPoint[2] = poz ; 
}
//////////////////////////////////////////////////////////////////////////////
const char* AliFlowV0::Pid() const
{
 // Returns the P.Id. in characters (,...) basing on the stored pdg code 
 // (MostLikelihoodPID()) .
 
 const char *name[] = {"gamma","K0","K0s","K0l","Lambda0"} ;
 int pdg_code = TMath::Abs(MostLikelihoodPID()) ;

 TString p_id = "" ;
 if(pdg_code == 22)	   { p_id = name[0] ; }
 else if(pdg_code == 311)  { p_id = name[1] ; }
 else if(pdg_code == 310)  { p_id = name[2] ; }
 else if(pdg_code == 130)  { p_id = name[3] ; }
 else if(pdg_code == 3122) { p_id = name[4] ; }
 // ...
 else 			   { p_id = "0" ; }
 
 return p_id.Data() ; 
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowV0::SetPid(const Char_t* pid)	   	
{ 
 // Sets the P.Id. hypotesis of the track from a String imput ("K0","K0s",
 // "K0l","Lambda0"). The string itself is not stored, what is stored is the
 // PDG code.
 
 if(strstr(pid,"gamma"))        { fMostLikelihoodPID = 22   ; }
 else if(strstr(pid,"K0"))      { fMostLikelihoodPID = 311  ; }
 else if(strstr(pid,"K0s"))     { fMostLikelihoodPID = 310  ; }
 else if(strstr(pid,"K0l"))     { fMostLikelihoodPID = 130  ; }
 else if(strstr(pid,"Lambda0")) { fMostLikelihoodPID = 3122 ; }
 // ...
 else { fMostLikelihoodPID = 0 ; cout << "AliFlowV0 - !BAD IMPUT!" << endl ; }
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowV0::Mass() 		  const { return fMass ; }
Float_t  AliFlowV0::Phi()	      	  const { return fPhi; }    	     
Float_t  AliFlowV0::Eta()	      	  const { return fEta; }    	     
Float_t  AliFlowV0::Pt()	      	  const { return fPt; }     	     
Short_t  AliFlowV0::Charge() 	      	  const { return 0 ; } 	     
Float_t  AliFlowV0::Dca()		  const { return fDca ; }		  
Float_t  AliFlowV0::Chi2()		  const { return fChi2; }		  
Float_t  AliFlowV0::CrossDca()    	  const { return fCrossDCA ; }
Float_t  AliFlowV0::V0Lenght()	  	  const { return TMath::Sqrt(fCrossPoint[0]*fCrossPoint[0] + fCrossPoint[1]*fCrossPoint[1] + fCrossPoint[2]*fCrossPoint[2]) ; }	  
Float_t  AliFlowV0::Sigma()	  	  const { return fSigma ; }	  
Int_t	 AliFlowV0::MostLikelihoodPID()   const { return fMostLikelihoodPID; } 
Int_t    AliFlowV0::Label()               const { return fLabel ; }	      
AliFlowTrack* AliFlowV0::DaughterP()      const { return fDaughters[0] ; }
AliFlowTrack* AliFlowV0::DaughterN()      const { return fDaughters[1] ; }
//////////////////////////////////////////////////////////////////////////////
void AliFlowV0::SetMostLikelihoodPID(Int_t val) { fMostLikelihoodPID = val ; } 
void AliFlowV0::SetVmass(Float_t mass) 	   	{ fMass = mass ; }
void AliFlowV0::SetPhi(Float_t phi)		{ fPhi = phi; }
void AliFlowV0::SetEta(Float_t eta)		{ fEta = eta; }
void AliFlowV0::SetPt(Float_t pt)		{ fPt = pt; }
void AliFlowV0::SetChi2(Float_t chi2)	   	{ fChi2 = chi2; }
void AliFlowV0::SetDca(Float_t dca)             { fDca = dca; }
void AliFlowV0::SetCrossDca(Float_t dca)        { fCrossDCA = dca ; }
void AliFlowV0::SetSigma(Float_t sigma)	  	{ fSigma = sigma ; }	  
void AliFlowV0::SetLabel(Int_t label)		{ fLabel = label ; }
//////////////////////////////////////////////////////////////////////////////
