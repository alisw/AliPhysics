//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowV0.cxx 18618 2007-05-16 15:38:22Z snelling $
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
#include "AliFlowTrack.h"
#include "AliFlowEvent.h"

#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowV0) 
//////////////////////////////////////////////////////////////////////////////
AliFlowV0::AliFlowV0():
 fPhi(0.), fPt(0.), fEta(0.), fChi2(0.), fMass(0.), fDca(0.), fCrossDCA(0.), fSigma(1.), fLabel(0), fMostLikelihoodPID(0), fPointAngle(0.),
 fDaughterP(-1), fDaughterN(-1)
{
 // default constructor  
  for(Int_t dd=0;dd<3;dd++) { fCrossPoint[dd] = 0. ; }
  // fDaughterP = -1  ; fDaughterN = -1  ;
}
//////////////////////////////////////////////////////////////////////////////
AliFlowV0::AliFlowV0(const Char_t* name):
 fPhi(0.), fPt(0.), fEta(0.), fChi2(0.), fMass(0.), fDca(0.), fCrossDCA(0.), fSigma(1.), fLabel(0), fMostLikelihoodPID(0), fPointAngle(0.),
 fDaughterP(-1), fDaughterN(-1)
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

 if(Pt()<=0)  { return 0. ; }
 if(Eta()==0) { return 0. ; }

 float momentum = Pt()/TMath::Sqrt(1-(tanh(Eta())*tanh(Eta()))) ; 
 return momentum; 
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowV0::Y() const 
{
 // Rapidity of the v0

 if(TMath::Abs((Float_t)P()) == TMath::Abs((Float_t)Pt())) 	{ return 0. ; }
 else if(TMath::Abs((Float_t)P()) < TMath::Abs((Float_t)Pt()))	{ cout << "v0: " << GetName() << " has  Pt() > P() !!!" << endl ; return -1000. ; }
 // -
 Float_t mass = Mass() ; 
 Double_t pz = TMath::Sqrt(P()*P() - Pt()*Pt()); 
 if(Eta() < 0) { pz = -pz ; }
 Double_t e = TMath::Sqrt(P()*P() + mass*mass) ;
 Float_t rapidity = 0.5 * TMath::Log((e + pz)/(e - pz)) ;
 return rapidity ;
}
//////////////////////////////////////////////////////////////////////////////
const char* AliFlowV0::Pid() const
{
 // Returns the P.Id. in char* basing on the stored pdg code (MostLikelihoodPID()) .
 
 const char *name[] = {"gamma","K0","K0s","K0l","Lambda0"} ;
 int pdgCode = TMath::Abs(MostLikelihoodPID()) ;

 TString pId = "" ;
 if(pdgCode == 22)	  { pId = name[0] ; }
 else if(pdgCode == 311)  { pId = name[1] ; }
 else if(pdgCode == 310)  { pId = name[2] ; }
 else if(pdgCode == 130)  { pId = name[3] ; }
 else if(pdgCode == 3122) { pId = name[4] ; }
 // ...
 else 			  { pId = "0" ; }
 
 return pId.Data() ; 
}
//////////////////////////////////////////////////////////////////////////////
// void AliFlowV0::SetPid(const Char_t* pid)	   	
// { 
//  // Sets the P.Id. hypotesis of the track from a char* imput ("K0","Lambda0",...). 
//  // The string itself is not stored, what is stored is the PDG code.
//  
//  if(strstr(pid,"gamma"))        { fMostLikelihoodPID = 22   ; }
//  else if(strstr(pid,"K0"))      { fMostLikelihoodPID = 311  ; }
//  else if(strstr(pid,"K0s"))     { fMostLikelihoodPID = 310  ; }
//  else if(strstr(pid,"K0l"))     { fMostLikelihoodPID = 130  ; }
//  else if(strstr(pid,"Lambda0")) { fMostLikelihoodPID = 3122 ; }
//  // ...
//  else { fMostLikelihoodPID = 0 ; cout << "AliFlowV0 - !BAD IMPUT!" << endl ; }
// }
//////////////////////////////////////////////////////////////////////////////
