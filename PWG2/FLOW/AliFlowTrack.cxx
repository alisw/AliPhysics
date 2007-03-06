//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
//_____________________________________________________________
//
// Description: 
//         an array of AliFlowTrack is the core of the AliFlowEvent. 
// The object AliFlowTrack contains data members wich summarize the track 
// information most useful for flow study (such as Pt, eta, phi angle, 
// p.id hypothesis, some detector informations like fitpoints, chi2, 
// energy loss, t.o.f., ecc.). 
// This class is optimized for reaction plane calculation and sub-event 
// selection, through the appropriate methods in AliFlowEvent.  
// Two arrays of flags in the AliFlowTrack object (fSelection[][] and 
// fSubevent[][]) specify whether a track is included or not in a certain 
// Selection or Subevent. Those flags are switched on/off by the method 
// AliFlowEvent::SetSelection(), according to Eta, Pt, dca, constrainability  
// of the track, number of TPC hits, and p.id hypotesis (static member of the
// AliFlowEvent class (fEtaTpcCuts[][][], fPtTpcCuts[][][], fDcaGlobalCuts[], 
// fPid[], ecc.). 
//
// The AliFlowTrack Class is adapted from the original StFlowTrack,
// succesfully employed to study flow in the STAR experiment at RICH.
// Original Authors:                 Raimond Snellings & Art Poskanzer
//

#include "AliFlowTrack.h"
#include "AliFlowConstants.h"

#include "TMath.h"
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowTrack)
//-------------------------------------------------------------
AliFlowTrack::AliFlowTrack()
{
 // Default constructor
 
 fPhi = 0. ;
 fEta = 0. ;
 fPt  = 0. ;
 fPhiGlobal = 0. ;
 fEtaGlobal = 0. ;
 fPtGlobal  = 0. ;
 fChi2 = 0. ;
 fZFirstPoint = 0. ;  		      
 fZLastPoint = 0. ;			      
 fTrackLength = 0. ;
 fMostLikelihoodPID = 0 ;		      
 for(Int_t ii=0;ii<2;ii++) 	                { fDcaSigned[ii] = 0. ; }			      
 for(Int_t ii=0;ii<AliFlowConstants::kPid;ii++) { fPidProb[ii]   = 0. ; }
 for(Int_t ii=0;ii<4;ii++)
 {
  fFitPts[4]  = 0  ; fMaxPts[4]  = 0  ; 
  fFitChi2[4] = 0. ; fDedx[4]    = 0. ;
  fMom[4] = 0. ;
 }
 ResetSelection() ;			      
}
//-------------------------------------------------------------
AliFlowTrack::AliFlowTrack(const Char_t* name) 
{
 // TNamed constructor
  
 SetName(name) ;
 AliFlowTrack() ;
}
//-------------------------------------------------------------
AliFlowTrack::~AliFlowTrack() 
{
 // default destructor (dummy)
}
//-------------------------------------------------------------

//-------------------------------------------------------------
Float_t  AliFlowTrack::P() const 
{ 
 // Returns the total momentum of the constrained track (calculated from Pt & Eta).
 // If the track is not constrainable or eta is infinite, returns 0.
 
 if(!IsConstrainable()) { return 0. ; }		      

 Float_t momentum = 0. ; 
 Float_t conv = 1 - ( TMath::TanH(Eta()) * TMath::TanH(Eta()) ) ;
 if(conv>0) { momentum = Pt() / TMath::Sqrt(conv) ; } 

 return momentum; 
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::PGlobal() const 
{ 
 // Returns the total momentum of the unconstrained track (calculated from PtGlobal & EtaGlobal). 
 // If etaGlobal is infinite, returns 0.

 Float_t momentum = 0. ; 
 Float_t conv = 1 - ( TMath::TanH(EtaGlobal()) * TMath::TanH(EtaGlobal()) ) ;
 if(conv>0) { momentum = PtGlobal() / TMath::Sqrt(conv) ; }

 return momentum; 
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::Mass() const 
{ 
 // Returns the mass of the track, basing on its P.Id. hypotesis

 Float_t mass = 0.13957 ;  // pion mass 

 if(MostLikelihoodPID() == 0) 	     		      { mass = -1.    ; }
 else if(TMath::Abs(MostLikelihoodPID()) == 11)	      { mass = 0.00051; }
 else if(TMath::Abs(MostLikelihoodPID()) == 13)	      { mass = 0.10566; }
 else if(TMath::Abs(MostLikelihoodPID()) == 211)      { mass = 0.49368; }
 else if(TMath::Abs(MostLikelihoodPID()) == 321)      { mass = 0.13957; }
 else if(TMath::Abs(MostLikelihoodPID()) == 2212)     { mass = 0.93827; }
 else if(TMath::Abs(MostLikelihoodPID()) == 10010020) { mass = 1.87505; }

 return mass ;
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::InvMass() const 
{ 
 // Returns the invariant mass of the track, calculated from T.O.F. and P_tot

 Float_t mass = -1 ;       	    // if no calculation possible, it returns -1
 Float_t c = TMath::Ccgs()/1e+12 ;  // = 0.02998 ;  --> speed of light in cm/psec
 Float_t tof = TofTOF() ; 	    // in pico-seconds
 Float_t lenght = TrackLength() ;   // in cm
 Float_t ptot = 0 ;   		    // constrained parameters are used if there, otherwise unconstrained
 Float_t beta, gamma ;
 if(IsConstrainable()) { ptot = P() ; }
 else                  { ptot = PGlobal() ; }
 if(tof>0 && lenght>0)
 {
  beta = lenght / (tof * c) ;
  if(beta<1) 
  {
   gamma = 1 / TMath::Sqrt( 1 - (beta * beta) ) ;
   mass = TMath::Abs( ptot / (beta * gamma) ) ;  
  }
  //cout << " (yes) Mass = " << mass << " , t/l = " << (tof/lenght) << " , beta = " << beta << " , gamma = " << gamma << endl ; 
 }   
 //else { cout << " (no)  TOF = " << tof << "	, Lenght = " << lenght << endl ; }

 return mass ;
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::Y() const 
{
 // Rapidity of the constrained track.
 // If track is not constrainable, the unconstrained rapidity is returned.
 
 if(!IsConstrainable()) { return YGlobal() ; }
 if(TMath::Abs((Float_t)P()) == TMath::Abs((Float_t)Pt())) 	{ return 0. ; }
 else if(TMath::Abs((Float_t)P()) < TMath::Abs((Float_t)Pt()))	{ cout << "track: " << GetName() << "has  Pt() > P() !!!" << endl ; }
 // -
 Float_t mass = Mass() ; 
 Double_t pz = TMath::Sqrt(P()*P() - Pt()*Pt()); 
 if(Eta()<0) { pz = -pz ; }
 Double_t e = TMath::Sqrt(P()*P() + mass*mass) ;
 Float_t rapidity = 0.5*TMath::Log((e + pz)/(e - pz)) ;
 
 return rapidity ;
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::YGlobal() const 
{
 // Rapidity of the unconstrained track

 if(TMath::Abs((Float_t)PGlobal()) == TMath::Abs((Float_t)PtGlobal())) 	   { return 0. ; }
 else if(TMath::Abs((Float_t)PGlobal()) < TMath::Abs((Float_t)PtGlobal())) { cout << "track: " << GetName() << "has  Pt() > P() !!!" << endl ; }
 // -
 Float_t mass = Mass() ; 
 Double_t pz = TMath::Sqrt(PGlobal()*PGlobal() - PtGlobal()*PtGlobal()); 
 if(EtaGlobal()<0) { pz = -pz ; }
 Double_t e = TMath::Sqrt(PGlobal()*PGlobal() + mass*mass) ;
 Float_t rapidity = 0.5*TMath::Log((e + pz)/(e - pz)) ;
 
 return rapidity ;
}
//-------------------------------------------------------------
const Char_t* AliFlowTrack::Pid() const
{
 // Returns the P.Id. in characters (e,mu,pi,k,p,d) basing on the stored pdg code 
 // and its sign (both stored in MostLikelihoodPID() ) .
 
 const Char_t *name[] = {"e","mu","pi","k","pr","d"} ;
 Int_t pdgCode = TMath::Abs(MostLikelihoodPID()) ;
 Int_t charge = Charge() ;

 TString pId = "" ;
 if(pdgCode == 11)            { pId = name[0] ; }
 else if(pdgCode == 13)       { pId = name[1] ; }
 else if(pdgCode == 211)      { pId = name[2] ; }
 else if(pdgCode == 321)      { pId = name[3] ; }
 else if(pdgCode == 2212)     { pId = name[4] ; }
 else if(pdgCode == 10010020) { pId = name[5] ; }
 else			      { pId = "0" ; }

 if(charge>0) 		      { pId += "+" ; } 
 else if(charge<0) 	      { pId += "-" ; }
 else		 	      { pId += "" ; }
 
 return pId.Data() ; 
}
//-------------------------------------------------------------
void AliFlowTrack::PidProbs(Float_t pidN[AliFlowConstants::kPid]) const 
{ 
 // Returns the normalized probability for the given track to be [e,mu,pi,k,p,d] 
 // The detector response is weighted by the bayesian vector of particles 
 // abundances, stored in AliFlowConstants::fgBayesian[] .

 Double_t sum = 0 ; 
 for(Int_t n=0;n<AliFlowConstants::kPid;n++)  { sum += fPidProb[n] * AliFlowConstants::fgBayesian[n] ; }
 if(sum)
 {
  for(Int_t n=0;n<AliFlowConstants::kPid;n++) { pidN[n] = fPidProb[n] * AliFlowConstants::fgBayesian[n] / sum ; }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
} 
//-------------------------------------------------------------
Float_t  AliFlowTrack::PidProb(Int_t nn)	const
{
 // Returns the normalized probability of the track to be [nn] (e,mu,pi,k,pi,d).

 Float_t pidN[AliFlowConstants::kPid] ; 
 PidProbs(pidN) ;
 return pidN[nn] ;
}
//-------------------------------------------------------------
TVector AliFlowTrack::PidProbs()  		const
{
 // Returns the normalized probability for the given track to be [e,mu,pi,k,p,d] 
 // as a TVector.

 TVector pidNvec(AliFlowConstants::kPid) ;
 Float_t pidN[AliFlowConstants::kPid] ; 
 PidProbs(pidN) ;
 for(Int_t n=0;n<AliFlowConstants::kPid;n++)  { pidNvec[n] = pidN[n] ; }

 return pidNvec ;
}
//-------------------------------------------------------------
void AliFlowTrack::RawPidProbs(Float_t pidV[AliFlowConstants::kPid]) const 
{ 
 // Returns the array of probabilities for the track to be [e,mu,pi,k,pi,d].

 for(Int_t ii=0;ii<AliFlowConstants::kPid;ii++) { pidV[ii] = fPidProb[ii] ; }
} 
//-------------------------------------------------------------
Float_t AliFlowTrack::MostLikelihoodProb()   	    const 
{ 
 // Returns the probability of the most probable P.id.
 // (Warning: THIS IS NOT WEIGHTED IN THE BAYESIAN WAY...) 

 Int_t pdgCode = TMath::Abs(MostLikelihoodPID()) ;
 if(pdgCode == 11)            { return fPidProb[0] ; }
 else if(pdgCode == 13)       { return fPidProb[1] ; }
 else if(pdgCode == 211)      { return fPidProb[2] ; }
 else if(pdgCode == 321)      { return fPidProb[3] ; }
 else if(pdgCode == 2212)     { return fPidProb[4] ; }
 else if(pdgCode == 10010020) { return fPidProb[5] ; }
 else { return 0. ; }
} 
//-------------------------------------------------------------
void AliFlowTrack::SetPid(const Char_t* pid)	   	
{ 
 // Sets the P.Id. hypotesis of the track from a String imput (that should be 
 // ("e","mu","pi","k","pr","d"). Sign is allowed as well. The string itself 
 // is not stored, what is stored is the signed PDG code.
 
 if(strstr(pid,"e"))  	   { fMostLikelihoodPID = 11	; }
 else if(strstr(pid,"mu")) { fMostLikelihoodPID = 13	; }
 else if(strstr(pid,"pi")) { fMostLikelihoodPID = 211	; }
 else if(strstr(pid,"k"))  { fMostLikelihoodPID = 321	; }
 else if(strstr(pid,"pr")) { fMostLikelihoodPID = 2212  ; }
 else if(strstr(pid,"d"))  { fMostLikelihoodPID = 10010020 ; }
 else { fMostLikelihoodPID = 0 ; cout << "AliFlowTrack - !BAD IMPUT!" << endl ; }
 
 if(strchr(pid,'+'))  	   { fMostLikelihoodPID = TMath::Abs(fMostLikelihoodPID) * 1 ; }
 else if(strchr(pid,'-'))  { fMostLikelihoodPID = TMath::Abs(fMostLikelihoodPID) * -1 ; } 
}
//-------------------------------------------------------------
Bool_t AliFlowTrack::IsConstrainable()      	    const 
{ 
 // Returns kTRUE if the track is constrainable. 
 // If the track is constrainable -> the constrained parameters are stored (and the
 // unconstrained ones as well in ...Global). 
 // Otherwise (if track is not constrainable) -> just the unconstrained parameters 
 // are stored, and returned instead of the constrained ones .
 
 if(fPt==0 && fEta==0) { return kFALSE ; }
 else 	               { return kTRUE ; }
} 
//-------------------------------------------------------------
Float_t AliFlowTrack::Dca() const 
{
 // distance of closest approach (norm(dca[2]))

 Double_t mag = 0 ;
 for(Int_t ii=0;ii<2;ii++) { mag += (fDcaSigned[ii]*fDcaSigned[ii]) ; } 
 return TMath::Sqrt(mag) ; 
}     
//-------------------------------------------------------------
Float_t AliFlowTrack::DcaError() const 
{
 // error on the distance of closest approach

 Double_t err = 0 ;
 for(Int_t ii=0;ii<2;ii++) { err = (fDcaError[ii]*fDcaError[ii]) ; } 
 return TMath::Sqrt(err) ; 
}     
//-------------------------------------------------------------
void AliFlowTrack::DcaError3(Float_t err[3]) const 
{
 // Dca Error (xy, z, cov(xy,z)) ...

 for(Int_t ii=0;ii<2;ii++) { err[ii] = TMath::Abs(fDcaError[ii]) ; }		
}     
//-------------------------------------------------------------

//-------------------------------------------------------------
Bool_t AliFlowTrack::Select(Int_t harmonic,Int_t selection,Int_t subevent) const 
{
 // Returns the selection flag for [harmonic] [selection] [sub-event]

 if((subevent == -1) || (subevent == fSubevent[harmonic][selection])) 
 {
  return fSelection[harmonic][selection] ; 
 } 
 return kFALSE ; 	
}
//-------------------------------------------------------------
void AliFlowTrack::PrintSelection() const
{
 // Prints a short string of 0 & 1 to visualize the selections' flags
 // [har1][sel0],[har1][sel1],[har1][sel2],...,[har2][sel0],...
 
 for(Int_t i=0;i<AliFlowConstants::kHars;i++)
 {
  for(Int_t j=0;j<AliFlowConstants::kSels;j++)
  {
   if(Select(i,j)) { cout << 1 ; }
   else    	   { cout << 0 ; }
  }
 }
 cout << endl ;
}
//-------------------------------------------------------------
void AliFlowTrack::ResetSelection() 
{
 // Re-sets all the selection/sub-event flags to 0 
 // (track won't be used in any R.P. calculation)

 for(Int_t ii=0;ii<AliFlowConstants::kHars;ii++)
 {
  for(Int_t jj=0;jj<AliFlowConstants::kSels;jj++)
  {
   fSelection[ii][jj] = kFALSE ; 
   fSubevent[ii][jj] = -1 ; 
  }
 }
}
//-------------------------------------------------------------

//-------------------------------------------------------------
void AliFlowTrack::SetConstrainable()         	   	      
{ 
 // fills the constrained parameters with the unconstrained ones, making it a constrainable track.
 //                                                   !!! TRICKY METHOD !!!

 if(!IsConstrainable()) 
 { 
  fPhi = fPhiGlobal ; 
  fEta = fEtaGlobal ; 
  fPt  = fPtGlobal ; 
 } 
}
//-------------------------------------------------------------
void AliFlowTrack::SetUnConstrainable()        	   	      
{ 
 // deletes the constrained parameters making it  an unconstrainable track. 
 //                                                   !!! TRICKY METHOD !!!

 if(IsConstrainable()) 
 { 
  fPhi = 0. ; 
  fEta = 0. ; 
  fPt  = 0. ; 
 } 
}
//-------------------------------------------------------------
