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
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliFlowTrack) ;
//////////////////////////////////////////////////////////////////////////////
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
 for(Int_t ii=0;ii<2;ii++) 	    { fDcaSigned[ii] = 0 ; }			      
 for(Int_t ii=0;ii<Flow::nPid;ii++) { fPidProb[ii] = 0. ; }
 for(Int_t ii=0;ii<4;ii++)
 {
  fFitPts[4]  = 0  ; fMaxPts[4]  = 0  ; 
  fFitChi2[4] = 0. ; fDedx[4]    = 0. ;
  fMom[4] = 0. ;
 }
 ResetSelection() ;			      
}
//////////////////////////////////////////////////////////////////////////////
AliFlowTrack::AliFlowTrack(const Char_t* name) 
{
 // TNamed constructor
  
 SetName(name) ;
 AliFlowTrack() ;
}
//////////////////////////////////////////////////////////////////////////////
AliFlowTrack::~AliFlowTrack() 
{
 // default destructor (dummy)
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::P() const 
{ 
 // Returns the total momentum of the constrained track (calculated from Pt & Eta),
 // if the track is not constrainable returns 0 (used as flag for contrainable tracks).
 
 if(!IsConstrainable()) { return 0 ; }		      
 Float_t momentum = Pt()/TMath::Sqrt(1-(TMath::TanH(Eta())*TMath::TanH(Eta()))); 
 return momentum; 
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::PGlobal()	      	    const 
{ 
 // Returns the total momentum of the unconstrained track (always>0). 

 Float_t momentum = PtGlobal()/TMath::Sqrt(1-(TMath::TanH(EtaGlobal())*TMath::TanH(EtaGlobal()))); 
 return momentum; 
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::Mass() const 
{ 
 // Returns the mass of the track, basing on its P.Id. hypotesis

 Float_t M = 0.13957 ;  // pion mass 

 if(MostLikelihoodPID() == 0) 	     		      { M = -1.    ; }
 else if(TMath::Abs(MostLikelihoodPID()) == 11)	      { M = 0.00051; }
 else if(TMath::Abs(MostLikelihoodPID()) == 13)	      { M = 0.10566; }
 else if(TMath::Abs(MostLikelihoodPID()) == 211)      { M = 0.49368; }
 else if(TMath::Abs(MostLikelihoodPID()) == 321)      { M = 0.13957; }
 else if(TMath::Abs(MostLikelihoodPID()) == 2212)     { M = 0.93827; }
 else if(TMath::Abs(MostLikelihoodPID()) == 10010020) { M = 1.87505; }

 return M ;
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::InvMass() const 
{ 
 // Returns the invariant mass of the track, calculated from T.O.F. and P_tot

 Float_t M = -1 ;       	    // if no calculation possible, it returns -1
 Float_t c = TMath::Ccgs()/1e+12 ;  // = 0.02998 ;  --> speed of light in cm/psec
 Float_t tof = TofTOF() ; 	    // in pico-seconds
 Float_t lenght = TrackLength() ;   // in cm
 Float_t Ptot = 0 ;   		    // constrained parameters are used if there, otherwise unconstrained
 if(IsConstrainable()) { Ptot = P() ; }
 else                  { Ptot = PGlobal() ; }
 if(tof>0 && lenght>0)
 { 
  Float_t beta = lenght / (tof * c) ;
  if(beta<1) 
  {
   Float_t gamma = 1/TMath::Sqrt(1-(beta*beta)) ;
   M = TMath::Abs(Ptot/(beta*gamma)) ;  
  } // cout << ii << ": Mass = " << cmass << " , t/l = " << (tof/lenght) << " , beta = " << beta << " , gamma = " << gamma << endl ; }
 }  // else { cout << ii << ": TOF = " << tof << "	, Lenght = " << lenght << endl ; }

 return M ;
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::Y() const 
{
 // Rapidity of the constrained track.
 // If track is not constrainable, the unconstrained rapidity is returned.
 
 if(!IsConstrainable()) { return YGlobal() ; }
 if(TMath::Abs((Float_t)P()) == TMath::Abs((Float_t)Pt())) 	{ return 0. ; }
 else if(TMath::Abs((Float_t)P()) < TMath::Abs((Float_t)Pt()))	{ cout << "track: " << GetName() << "has  Pt() > P() !!!" << endl ; }
 // -
 float M = Mass() ; 
 double Pz = TMath::Sqrt(P()*P() - Pt()*Pt()); 
 if(Eta()<0) { Pz = -Pz ; }
 double E = TMath::Sqrt(P()*P() + M*M) ;
 float rapidity = 0.5*TMath::Log((E + Pz)/(E - Pz)) ;
 
 return rapidity ;
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::YGlobal() const 
{
 // Rapidity of the unconstrained track

 if(TMath::Abs((Float_t)PGlobal()) == TMath::Abs((Float_t)PtGlobal())) 	   { return 0. ; }
 else if(TMath::Abs((Float_t)PGlobal()) < TMath::Abs((Float_t)PtGlobal())) { cout << "track: " << GetName() << "has  Pt() > P() !!!" << endl ; }
 // -
 float M = Mass() ; 
 double Pz = TMath::Sqrt(PGlobal()*PGlobal() - PtGlobal()*PtGlobal()); 
 if(EtaGlobal()<0) { Pz = -Pz ; }
 double E = TMath::Sqrt(PGlobal()*PGlobal() + M*M) ;
 float rapidity = 0.5*TMath::Log((E + Pz)/(E - Pz)) ;
 
 return rapidity ;
}
//////////////////////////////////////////////////////////////////////////////
Bool_t AliFlowTrack::Select(Int_t harmonic,Int_t selection,Int_t subevent) const 
{
 // Returns the selection flag for [harmonic] [selection] [sub-event]

 if((subevent == -1) || (subevent == fSubevent[harmonic][selection])) 
 {
  if(fSelection[harmonic][selection]) { return kTRUE ; }
 } 
 return kFALSE ; 	
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::SetSelect(Int_t harmonic,Int_t selection) 
{ 
 fSelection[harmonic][selection] = kTRUE ;
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::SetSubevent(Int_t harmonic,Int_t selection,Int_t subevent) 
{ 
 fSubevent[harmonic][selection] = subevent ; 
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::PrintSelection() 
{
 // Prints a short string of 0 & 1 to visualize the selections' flags
 // [har1][sel0],[har1][sel1],[har1][sel2],...,[har2][sel0],...
 
 for(int i=0;i<Flow::nHars;i++)
 {
  for(int j=0;j<Flow::nSels;j++)
  {
   if(Select(i,j)) { cout << 1 ; }
   else    	   { cout << 0 ; }
  }
 }
 cout << endl ;
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::ResetSelection() 
{
 // Re-sets all the selection/sub-event flags to 0 
 // (track won't be used in any R.P. calculation)

 for(Int_t ii=0;ii<Flow::nHars;ii++)
 {
  for(Int_t jj=0;jj<Flow::nSels;jj++)
  {
   fSelection[ii][jj] = kFALSE ; 
   fSubevent[ii][jj] = -1 ; 
  }
 }
}
//////////////////////////////////////////////////////////////////////////////
const char* AliFlowTrack::Pid() const
{
 // Returns the P.Id. in characters (e,mu,pi,k,p,d) basing on the stored pdg code 
 // and its sign (both stored in MostLikelihoodPID() ) .
 
 const char *name[] = {"e","mu","pi","k","pr","d"} ;
 int pdg_code = TMath::Abs(MostLikelihoodPID()) ;
 int charge = Charge() ;

 TString p_id = "" ;
 if(pdg_code == 11)            { p_id = name[0] ; }
 else if(pdg_code == 13)       { p_id = name[1] ; }
 else if(pdg_code == 211)      { p_id = name[2] ; }
 else if(pdg_code == 321)      { p_id = name[3] ; }
 else if(pdg_code == 2212)     { p_id = name[4] ; }
 else if(pdg_code == 10010020) { p_id = name[5] ; }
 else			       { p_id = "0" ; }

 if(charge>0) 		{ p_id += "+" ; } 
 else if(charge<0) 	{ p_id += "-" ; }
 else		 	{ p_id += "" ; }
 
 return p_id.Data() ; 
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::PidProbs(Float_t pidN[Flow::nPid]) const 
{ 
 // Returns the normalized probability for the given track to be [e,mu,pi,k,p,d] 
 // The detector response is weighted by the bayesian vector of particles 
 // abundances, stored in Flow::fBayesian[] .

 Double_t sum = 0 ; 
 for(Int_t n=0;n<Flow::nPid;n++)  { sum += fPidProb[n] * Flow::fBayesian[n] ; }
 if(sum)
 {
  for(Int_t n=0;n<Flow::nPid;n++) { pidN[n] = fPidProb[n] * Flow::fBayesian[n] / sum ; }
 }
 else { cout << " ERROR - Empty Bayesian Vector !!! " << endl ; }
} 
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::PidProb(Int_t nn)	const
{
 // Returns the normalized probability of the track to be [nn] (e,mu,pi,k,pi,d).

 Float_t pidN[Flow::nPid] ; 
 PidProbs(pidN) ;
 return pidN[nn] ;
}
//////////////////////////////////////////////////////////////////////////////
TVector AliFlowTrack::PidProbs()  		const
{
 // Returns the normalized probability for the given track to be [e,mu,pi,k,p,d] 
 // as a TVector.

 TVector pidNvec(Flow::nPid) ;
 Float_t pidN[Flow::nPid] ; 
 PidProbs(pidN) ;
 for(Int_t n=0;n<Flow::nPid;n++)  { pidNvec[n] = pidN[n] ; }

 return pidNvec ;
}
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::RawPidProbs(Float_t pidV[Flow::nPid]) const 
{ 
 // Returns the array of probabilities for the track to be [e,mu,pi,k,pi,d].

 for(Int_t ii=0;ii<Flow::nPid;ii++) { pidV[ii] = fPidProb[ii] ; }
} 
//////////////////////////////////////////////////////////////////////////////
Float_t AliFlowTrack::MostLikelihoodProb()   	    const 
{ 
 // Returns the probability of the most probable P.id.
 // (Warning: THIS IS NOT WEIGHTED IN THE BAYESIAN WAY...) 

 int pdg_code = TMath::Abs(MostLikelihoodPID()) ;
 if(pdg_code == 11)            { return fPidProb[0] ; }
 else if(pdg_code == 13)       { return fPidProb[1] ; }
 else if(pdg_code == 211)      { return fPidProb[2] ; }
 else if(pdg_code == 321)      { return fPidProb[3] ; }
 else if(pdg_code == 2212)     { return fPidProb[4] ; }
 else if(pdg_code == 10010020) { return fPidProb[5] ; }
 else { return 0. ; }
} 
//////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////
Int_t AliFlowTrack::Charge() const 
{ 
 // Tracks charge (stored as sign of fMostLikelihoodPID). +/-2 if deuterium.
 
 if(MostLikelihoodPID() == 10010020) { return TMath::Sign(2,MostLikelihoodPID()) ; }
 else 				     { return TMath::Sign(1,MostLikelihoodPID()) ; }
} 	      
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::SetCharge(Int_t charge)
{ 
 // Sets the charge of the track (+/- sign of MostLikelihoodPID()).
 
 fMostLikelihoodPID = TMath::Sign(TMath::Abs(fMostLikelihoodPID),charge) ;
}
//////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////
Float_t AliFlowTrack::Dca()		  	    const 
{
 Double_t mag = 0 ;
 for(Int_t ii=0;ii<2;ii++) { mag += (fDcaSigned[ii]*fDcaSigned[ii]) ; } 
 return TMath::Sqrt(mag) ; 
}     
//////////////////////////////////////////////////////////////////////////////
Float_t AliFlowTrack::Dca2(Float_t dca[2])	    const 
{ 
 Double_t mag = 0 ;
 for(Int_t ii=0;ii<2;ii++) 
 { 
  dca[ii] = fDcaSigned[ii] ; 
  mag += (fDcaSigned[ii]*fDcaSigned[ii]) ; 
 } 
 return TMath::Sqrt(mag) ; 
}
//////////////////////////////////////////////////////////////////////////////
Float_t AliFlowTrack::Dca3(Float_t dca[3])	    const 
{
 dca[0] = TMath::Abs(fDcaSigned[0])* TMath::Cos(Phi()) ;
 dca[1] = TMath::Abs(fDcaSigned[0])* TMath::Sin(Phi()) ;
 dca[2] = TMath::Abs(fDcaSigned[1]) ; 		
 return Dca() ; 
}     
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::Phi()		      	    const 
{ 
 if(IsConstrainable()) { return fPhi ; }
 else 		{ return PhiGlobal() ; }
}		    
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::Pt()		      	    const 
{  
 if(IsConstrainable())  { return fPt ; }
 else 		 { return PtGlobal() ; }
}
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::Eta()		      	    const 
{  
 if(IsConstrainable())  { return fEta ; }
 else 		 { return EtaGlobal() ; }
}		      
//////////////////////////////////////////////////////////////////////////////
Float_t  AliFlowTrack::PhiGlobal()		    const { return fPhiGlobal ; }		    
Float_t  AliFlowTrack::PtGlobal()		    const { return fPtGlobal ; }		    
Float_t  AliFlowTrack::EtaGlobal()		    const { return fEtaGlobal ; }		      
//////////////////////////////////////////////////////////////////////////////
Float_t AliFlowTrack::TransDcaSigned() 		    const { return fDcaSigned[0] ; }		
Float_t AliFlowTrack::TransDca() 		    const { return TMath::Abs(TransDcaSigned()) ; }	
Float_t  AliFlowTrack::Chi2()		      	    const { return fChi2 ; }		    
Float_t  AliFlowTrack::TrackLength()	      	    const { return fTrackLength ; }	    
Float_t  AliFlowTrack::ZFirstPoint()	      	    const { return fZFirstPoint ; }	    
Float_t  AliFlowTrack::ZLastPoint()	      	    const { return fZLastPoint ; }	    

Int_t	 AliFlowTrack::MostLikelihoodPID()    	    const { return fMostLikelihoodPID; } 
Float_t  AliFlowTrack::ElectronPositronProb() 	    const { return PidProb(0) ; }
Float_t  AliFlowTrack::MuonPlusMinusProb()    	    const { return PidProb(1) ; }
Float_t  AliFlowTrack::PionPlusMinusProb()    	    const { return PidProb(2) ; }
Float_t  AliFlowTrack::KaonPlusMinusProb()    	    const { return PidProb(3) ; }
Float_t  AliFlowTrack::ProtonPbarProb()       	    const { return PidProb(4) ; }
Float_t  AliFlowTrack::DeuteriumAntiDeuteriumProb() const { return PidProb(5) ; }

Int_t	 AliFlowTrack::FitPtsTPC() 	      	    const { return fFitPts[0]  ; }  
Int_t	 AliFlowTrack::MaxPtsTPC() 	      	    const { return fMaxPts[0]  ; }  
Float_t  AliFlowTrack::Chi2TPC()	      	    const { return fFitChi2[0] ; }	    
Float_t  AliFlowTrack::DedxTPC()	      	    const { return fDedx[0]    ; }  
Int_t	 AliFlowTrack::FitPtsITS()  	            const { return fFitPts[1]  ; }	 
Int_t	 AliFlowTrack::MaxPtsITS() 	      	    const { return fMaxPts[1]  ; }
Float_t  AliFlowTrack::Chi2ITS()	      	    const { return fFitChi2[1] ; }		 
Float_t  AliFlowTrack::DedxITS()	      	    const { return fDedx[1]    ; }  
Int_t	 AliFlowTrack::NhitsTRD()  	            const { return fFitPts[2]  ; }	 
Int_t	 AliFlowTrack::MaxPtsTRD() 	      	    const { return fMaxPts[2]  ; }  
Float_t  AliFlowTrack::Chi2TRD()	      	    const { return fFitChi2[2] ; }		 
Float_t  AliFlowTrack::SigTRD()	      	    	    const { return fDedx[2]    ; }  
Int_t	 AliFlowTrack::NhitsTOF()  	            const { return fFitPts[3]  ; }	 
Int_t	 AliFlowTrack::MaxPtsTOF() 	      	    const { return fMaxPts[3]  ; }  
Float_t  AliFlowTrack::Chi2TOF()	      	    const { return fFitChi2[3] ; }		 
Float_t  AliFlowTrack::TofTOF()		      	    const { return fDedx[3]    ; }  

Float_t  AliFlowTrack::PatTPC()      		    const { return fMom[0] ; }
Float_t  AliFlowTrack::PatITS()      		    const { return fMom[1] ; }
Float_t  AliFlowTrack::PatTRD()      		    const { return fMom[2] ; }
Float_t  AliFlowTrack::PatTOF()      		    const { return fMom[3] ; }

Int_t    AliFlowTrack::Label()                	    const { return fLabel ; }	      
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::SetConstrainable()         	   	      
{ 
 // fills the constrained parameters with the unconstrained ones, making it
 // a constrainable track.                            !!! TRICKY METHOD !!!

 if(!IsConstrainable()) 
 { 
  fPhi = fPhiGlobal ; 
  fEta = fEtaGlobal ; 
  fPt  = fPtGlobal ; 
 } 
}
//////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////
void AliFlowTrack::SetPhi(Float_t phi)		   	      { fPhi = phi; }
void AliFlowTrack::SetEta(Float_t eta)		   	      { fEta = eta; }
void AliFlowTrack::SetPt(Float_t pt)		   	      { fPt = pt; }
void AliFlowTrack::SetPhiGlobal(Float_t phi)		      { fPhiGlobal = phi; }
void AliFlowTrack::SetEtaGlobal(Float_t eta)		      { fEtaGlobal = eta; }
void AliFlowTrack::SetPtGlobal(Float_t pt)		      { fPtGlobal = pt; }
void AliFlowTrack::SetDcaSigned(Float_t xy, Float_t z)        { fDcaSigned[0] = xy ; fDcaSigned[1] = z ; }
void AliFlowTrack::SetTransDcaSigned(Float_t xy)  	      { fDcaSigned[0] = xy ; }
void AliFlowTrack::SetChi2(Float_t chi2)	   	      { fChi2 = chi2; }
void AliFlowTrack::SetTrackLength(Float_t tl)	   	      { fTrackLength = tl; }
void AliFlowTrack::SetZFirstPoint(Float_t zFirst)  	      { fZFirstPoint = zFirst; }
void AliFlowTrack::SetZLastPoint(Float_t zLast)    	      { fZLastPoint = zLast; }

void AliFlowTrack::SetMostLikelihoodPID(Int_t val)            { fMostLikelihoodPID = val ; }   // feed here the Signed pdg_code 
void AliFlowTrack::SetElectronPositronProb(Float_t val)       { fPidProb[0] = TMath::Abs(val) ; } 
void AliFlowTrack::SetMuonPlusMinusProb(Float_t val)	      { fPidProb[1] = TMath::Abs(val) ; } 
void AliFlowTrack::SetPionPlusMinusProb(Float_t val)	      { fPidProb[2] = TMath::Abs(val) ; } 
void AliFlowTrack::SetKaonPlusMinusProb(Float_t val)	      { fPidProb[3] = TMath::Abs(val) ; } 
void AliFlowTrack::SetProtonPbarProb(Float_t val)  	      { fPidProb[4] = TMath::Abs(val) ; } 
void AliFlowTrack::SetDeuteriumAntiDeuteriumProb(Float_t val) { fPidProb[5] = TMath::Abs(val) ; }

void AliFlowTrack::SetFitPtsTPC(Int_t fitPts)	     	      { fFitPts[0]  = fitPts ; }
void AliFlowTrack::SetMaxPtsTPC(Int_t maxPts)		      { fMaxPts[0]  = maxPts ; }
void AliFlowTrack::SetChi2TPC(Float_t chi2)	  	      { fFitChi2[0] = chi2   ; }
void AliFlowTrack::SetDedxTPC(Float_t dedx)	   	      { fDedx[0]    = dedx   ; }
void AliFlowTrack::SetFitPtsITS(Int_t nhits)  	     	      { fFitPts[1]  = nhits  ; }
void AliFlowTrack::SetMaxPtsITS(Int_t maxPts)		      { fMaxPts[1]  = maxPts ; }
void AliFlowTrack::SetChi2ITS(Float_t chi2)		      { fFitChi2[1] = chi2   ; }
void AliFlowTrack::SetDedxITS(Float_t dedx)	   	      { fDedx[1]    = dedx   ; }
void AliFlowTrack::SetNhitsTRD(Int_t fitPts)	 	      { fFitPts[2]  = fitPts ; }
void AliFlowTrack::SetMaxPtsTRD(Int_t maxPts)		      { fMaxPts[2]  = maxPts ; }
void AliFlowTrack::SetChi2TRD(Float_t chi2)	   	      { fFitChi2[2] = chi2   ; }
void AliFlowTrack::SetSigTRD(Float_t dedx)	   	      { fDedx[2]    = dedx   ; }
void AliFlowTrack::SetNhitsTOF(Int_t fitPts)	     	      { fFitPts[3]  = fitPts ; }
void AliFlowTrack::SetMaxPtsTOF(Int_t maxPts)		      { fMaxPts[3]  = maxPts ; }
void AliFlowTrack::SetChi2TOF(Float_t chi2)	   	      { fFitChi2[3] = chi2   ; }
void AliFlowTrack::SetTofTOF(Float_t dedx)	   	      { fDedx[3]    = dedx   ; }

void AliFlowTrack::SetPatTPC(Float_t p)  	              { fMom[0] = p ; }
void AliFlowTrack::SetPatITS(Float_t p)  	              { fMom[1] = p ; }
void AliFlowTrack::SetPatTRD(Float_t p)  	              { fMom[2] = p ; }
void AliFlowTrack::SetPatTOF(Float_t p)  	              { fMom[3] = p ; }

void AliFlowTrack::SetLabel(Int_t label)	   	      { fLabel = label ; }
//////////////////////////////////////////////////////////////////////////////
