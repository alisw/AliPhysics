//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowTrack.cxx 18618 2007-05-16 15:38:22Z snelling $
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
// information most useful for flow study (such as Pt, eta, phi angle), some 
// detector signals (fitpoints, chi2 of the track fit, energy loss, time of 
// flight, ecc.), the p.id is stored as the detector response function, for 
// each detector (with the possibility to costomly combine them), and the 
// combined one (see bayesian P.Id. - chap.5 of the ALICE PPR). 
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
AliFlowTrack::AliFlowTrack():
fPhi(0.), fEta(0.), fPt(0.), fZFirstPoint(0.), fZLastPoint(0.), fChi2(0.), fTrackLength(0.), fMostLikelihoodPID(0), fPhiGlobal(0.), fEtaGlobal(0.), fPtGlobal(0.), fLabel(0) {
 // Default constructor
 
 for(Int_t dd=0;dd<2;dd++) 	                { fDcaSigned[dd] = 0. ; }			      
 for(Int_t ii=0;ii<AliFlowConstants::kPid;ii++) { fCombRespFun[ii]   = 0. ; }
 for(Int_t det=0;det<4;det++)
 {
  fFitPts[det]  = 0  ; fMaxPts[det]   = 0  ; 
  fFitChi2[det] = 0. ; fDedx[det]     = 0. ; fMom[det] = 0. ; 
  for(Int_t ii=0;ii<AliFlowConstants::kPid;ii++) { fRespFun[det][ii] = -1. ; }
 }
 ResetSelection() ;			      
}
//-------------------------------------------------------------
AliFlowTrack::AliFlowTrack(const Char_t* name):
fPhi(0.), fEta(0.), fPt(0.), fZFirstPoint(0.), fZLastPoint(0.), fChi2(0.), fTrackLength(0.), fMostLikelihoodPID(0), fPhiGlobal(0.), fEtaGlobal(0.), fPtGlobal(0.), fLabel(0) {
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
void AliFlowTrack::PidProbs(Float_t *pidN) const 
{ 
 // Returns the normalized probability (the "bayesian weights") for the given track to be [e,mu,pi,k,p,d] .  
 // The COMBINED detector response function is scaled by the a priori probabilities 
 // (the normalised particle abundances is stored in AliFlowConstants::fgBayesian[]) .

 Double_t sum = 0 ; 
 for(Int_t n=0;n<AliFlowConstants::kPid;n++) { pidN[n] = PidProb(n) ; sum += pidN[n] ; }
 if(sum) 
 {  
  for(Int_t n=0;n<AliFlowConstants::kPid;n++)  { pidN[n] /=  sum ; }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
} 
//-------------------------------------------------------------
void AliFlowTrack::PidProbsC(Float_t *pidN) const 
{ 
 // Returns the normalized probability (the "bayesian weights") for the given track to be [e,mu,pi,k,p,d] .  
 // The CUSTOM detector response function (see AliFlowTrack) is scaled by the a priori probabilities 
 // (the normalised particle abundances is stored in AliFlowConstants::fgBayesian[]) .
 
 Double_t sum = 0 ; 
 for(Int_t n=0;n<AliFlowConstants::kPid;n++) { pidN[n] = PidProbC(n) ; sum += pidN[n] ; }
 if(sum) 
 {  
  for(Int_t n=0;n<AliFlowConstants::kPid;n++)  { pidN[n] /=  sum ; }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
} 
//-------------------------------------------------------------
Float_t  AliFlowTrack::PidProb(Int_t nPid) const
{
 // Returns the bayesian weight of the track to be [nPid = e,mu,pi,k,pi,d].
 // The detector response function in use is the combined one (from the ESD)

 if(nPid > AliFlowConstants::kPid) { return 0. ; } 

 return (fCombRespFun[nPid] * AliFlowConstants::fgBayesian[nPid]) ;
}
//-------------------------------------------------------------
Float_t  AliFlowTrack::PidProbC(Int_t nPid) const
{
 // Returns the bayesian weight of the track to be [nPid = e,mu,pi,k,pi,d].
 // The detector response function in use is the custom one ...

 if(nPid > AliFlowConstants::kPid) { return 0. ; } 

 Float_t  customRespFun[AliFlowConstants::kPid] ;
 GetCustomRespFun(customRespFun) ;

 return (customRespFun[nPid] * AliFlowConstants::fgBayesian[nPid]) ;
}
//-------------------------------------------------------------
Float_t AliFlowTrack::MostLikelihoodRespFunc() const 
{ 
 // Returns the detector response function for the most probable P.id. hypothesis
 // (Warning: THIS IS NOT THE BAYESIAN WEIGHT !) 

 Int_t pdgCode = TMath::Abs(MostLikelihoodPID()) ;
 if(pdgCode == 11)            { return fCombRespFun[0] ; }
 else if(pdgCode == 13)       { return fCombRespFun[1] ; }
 else if(pdgCode == 211)      { return fCombRespFun[2] ; }
 else if(pdgCode == 321)      { return fCombRespFun[3] ; }
 else if(pdgCode == 2212)     { return fCombRespFun[4] ; }
 else if(pdgCode == 10010020) { return fCombRespFun[5] ; }
 else { return 0. ; }
} 
//-------------------------------------------------------------
void AliFlowTrack::SetRespFun(Int_t det, Float_t *r)
{
 // This function fills "AliFlowConstants::kPid" PID weights 
 // (detector response functions) of the track .
 // The method is private, cause it is called by the public 
 // methods: SetRespFunITS, SetRespFunTPC, ...

 for(Int_t i=0;i<AliFlowConstants::kPid;i++)
 {
  fRespFun[det][i] = r[i] ;
 }
}
//-------------------------------------------------------------
void AliFlowTrack::GetRespFun(Int_t det, Float_t *r) const
{
 // This function returns the response functions of the 
 // detector [det] for the P.Id. of the track .
 // The method is private, and it is called by the public 
 // methods: GetRespFunITS, GetRespFunTPC, ...

 for(Int_t i=0;i<AliFlowConstants::kPid;i++)
 {
  r[i] = fRespFun[det][i] ;
 }
}
//-------------------------------------------------------------
void AliFlowTrack::GetCombinedRespFun(Float_t *r) const 
{
 // This function returns the combined response functions for  
 // [e,mu,pi,k,pi,d] as it is stored in the ESD .
 
 for(Int_t i=0;i<AliFlowConstants::kPid;i++)
 {
  r[i] = fCombRespFun[i] ;
 }
}
//------------------------------------------------------------- 
void AliFlowTrack::GetCustomRespFun(Float_t *r) const
{
 // This function returns a combined response functions as setted 
 // by the user ... at the moment just the sum of single responses ...  
 // for the track to be [e,mu,pi,k,pi,d] .
 
 Int_t sum ;
 for(Int_t i=0;i<AliFlowConstants::kPid;i++)
 {
  r[i] = 0. ; sum = 0 ;
  for(Int_t det=0;det<4;det++) 
  { 
   if(fRespFun[det][i]>0) { r[i] += fRespFun[det][i] ; sum += 1 ; } 
  }
  if(sum) { r[i] /= sum ; }
 }
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
void AliFlowTrack::SetPid(Int_t pdgCode)	   	
{ 
 // Sets the P.Id. hypotesis of the track from the PDG code. 
 // Sign can be given as well. 
 
 fMostLikelihoodPID = pdgCode ;
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
 // fills the constrained parameters with the unconstrained ones,
 // making the track a constrainable track.
 //                                       !!! TRICKY METHOD !!!

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
 // deletes the constrained parameters making the track an 
 // unconstrainable track. 
 //                                       !!! TRICKY METHOD !!!

 if(IsConstrainable()) 
 { 
  fPhi = 0. ; 
  fEta = 0. ; 
  fPt  = 0. ; 
 } 
}
//-------------------------------------------------------------
