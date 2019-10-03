//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowSelection.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description:  
//       the AliFlowSelection Class selects tracks from the AliFlowEvent. 
// There are two kind of selection :
// - Selection of tracks for the R.P. determination is performed via the
//  method Select(AliFlowTrack*): it simply checks the labels array of 
//  the track. Those arrays must be previously filled by calling 
//  AliFlowEvent::SetSelection(), wich uses the static cuts defined here.
//  Set the flag to Select via SetHarmonic() and SetSelection().   
// - Selection of tracks for Correlation Analysis is performed via the 
//  method SelectPart(AliFlowEvent*): cuts are applied to tracks 
//  (such as p.id, pT, total P, rapidity, eta, number of fit points, 
//  ratio of fit points to max points, chi2, and global dca). 
//  Set the cuts to apply via the appropriate Set...() methods.
// 
// The AliFlowSelection Class is adapted from the original StFlowSelection,
// succesfully used to study flow in the STAR experiment at RICH.
// Original Authors:                 Raimond Snellings & Art Poskanzer
//

#include "AliFlowSelection.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
using namespace std; //required for resolving the 'cout' symbol

// - 1st selection (both harmonic) is Disabled ! - the 2 following are the same (basic cuts, no P.id)
Float_t  AliFlowSelection::fgEtaTpcCuts[2][AliFlowConstants::kHars][AliFlowConstants::kSels] = {{{1.0,0.0},{1.0,0.0}},{{0.0,1.1},{0.0,1.1}}} ;
Float_t  AliFlowSelection::fgPtTpcCuts[2][AliFlowConstants::kHars][AliFlowConstants::kSels]  = {{{1.0,0.1},{1.0,0.1}},{{0.0,10.0},{0.0,10.0}}} ;
Float_t  AliFlowSelection::fgDcaGlobalCuts[2] 	  = { 0. , 1. } ;
Char_t   AliFlowSelection::fgPid[10]              = { '\0' } ;
Bool_t   AliFlowSelection::fgConstrainable        = kTRUE ;
Int_t    AliFlowSelection::fgTPChits[AliFlowConstants::kSels] = { 0 , 1 } ;

ClassImp(AliFlowSelection)
//-----------------------------------------------------------------------
AliFlowSelection::AliFlowSelection():
  fHarmonic(-1), fSelection(-1), fSubevent(-1), fPtBinsPart(AliFlowConstants::kPtBinsPart), fCent(-1), fRun(-1), fV0SideBand(0.), fConstrainablePart(kFALSE)
{
 // Default constructor: when initialized all selection cuts are disabled (lo>hi).

 fPidPart[0] = '\0' ;
 fPidProbPart[0]       = 1 ;    	// = 0.  
 fPidProbPart[1]       = 0 ;    	// = 1.  
 fPtPart[0]	       = 1 ; 		// = 0.  
 fPtPart[1]	       = 0 ; 		// = 10. 
 fPPart[0]	       = 1 ; 		// = 0.  
 fPPart[1]	       = 0 ; 		// = 10. 
 fEtaPart[0]	       = 1 ; 		// = -0.9
 fEtaPart[1]	       = 0 ; 		// = 0.9 
 fEtaAbsPart[0]	       = 1 ; 		// = 0.
 fEtaAbsPart[1]	       = 0 ; 		// = 0.9 
 fFitPtsPart[0]        = 1 ; 		// = 0   
 fFitPtsPart[1]        = 0 ; 		// = 0   
 fDedxPtsPart[0]       = 1 ;    	// = 0   
 fDedxPtsPart[1]       = 0 ;    	// = 0   
 fFitOverMaxPtsPart[0] = 1 ; 		// = 0.  
 fFitOverMaxPtsPart[1] = 0 ; 		// = 0.  
 fChiSqPart[0]         = 1 ; 		// = 0.  
 fChiSqPart[1]         = 0 ; 		// = 100.
 fDcaGlobalPart[0]     = 1 ; 		// = 0.  
 fDcaGlobalPart[1]     = 0 ; 		// = 1.  
 fDcaOverSigma[0]      = 1 ; 		// = 0. ;
 fDcaOverSigma[1]      = 0 ; 		// = 1. ;
 fYPart[0]	       = 1 ; 		// = 0.  
 fYPart[1]	       = 0 ; 		// = 0.  
 // -
 fV0Pt[0]	       = 1 ;  
 fV0Pt[1]	       = 0 ;  
 fV0P[0]	       = 1 ;  
 fV0P[1]	       = 0 ;  
 fV0Eta[0]	       = 1 ;  
 fV0Eta[1]	       = 0 ;  
 fV0EtaAbs[0]          = 1 ;  
 fV0EtaAbs[1]          = 0 ;  
 fV0Y[0]	       = 1 ;  
 fV0Y[1]	       = 0 ;  
 fV0ChiSq[0]	       = 1 ;  
 fV0ChiSq[1]	       = 0 ;  
 fV0Lenght[0]          = 1 ;  
 fV0Lenght[1]          = 0 ;  
 fV0LenghtOverSigma[0] = 1 ;  
 fV0LenghtOverSigma[1] = 0 ;  
 fV0DcaCross[0]        = 1 ;  
 fV0DcaCross[1]        = 0 ;  
 fV0Mass[0]	       = 1 ;
 fV0Mass[1]	       = 0 ;
}
//-----------------------------------------------------------------------
AliFlowSelection::~AliFlowSelection() 
{
 // default destructor (dummy)
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowEvent* pFlowEvent) const 
{
 // Returns kTRUE if the event is selected. 

 if((fCent>=0) && ((Int_t)pFlowEvent->Centrality() != fCent)) { return kFALSE ; }
 //if(pFlowEvent->RunID() != fRun) { return kFALSE ; }
 //if(pFlowEvent->... != ...) { return kFALSE ; }
 
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowTrack* pFlowTrack) const 
{
 // Selects particles for event plane determination.
 // Returns kTRUE if the track is selected.
 
 if(!pFlowTrack->Select(fHarmonic, fSelection, fSubevent)) { return kFALSE ; }
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowV0* pFlowV0) const  
{
 // Returns kTRUE if the v0 is selected. (dummy) 

 if(!pFlowV0) { return kFALSE ; }
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::SelectPart(AliFlowTrack* pFlowTrack) const  
{
 // Make track selection for Correlation Analysis & Vn (charged particles 
 // to correlate with the event plane).  
 // Returns kTRUE if the track is selected.

 // PID
 if(fPidPart[0] != '\0') 
 {
  if(strstr(fPidPart, "h")!=0)					// if fPidPart contains the char "h" 
  {
   int charge = pFlowTrack->Charge();
   if(strcmp("h+", fPidPart)==0 && charge != 1)  return kFALSE;	// if fPidPart == "h+"
   if(strcmp("h-", fPidPart)==0 && charge != -1) return kFALSE;
  } 
  else 
  {
   const Char_t* pid = pFlowTrack->Pid() ;
   if(strstr(pid, fPidPart)==0) return kFALSE;			// if pFlowTrack->Pid() does not contain fPidPart (RP. selected parts.)
  }
 } 
 
 // PID probability
 float pidProb = pFlowTrack->MostLikelihoodRespFunc() ;
 if(fPidProbPart[1] > fPidProbPart[0] &&  (pidProb < fPidProbPart[0] || pidProb > fPidProbPart[1])) return kFALSE;
 
 // Constrainable
 bool constrainable = pFlowTrack->IsConstrainable() ;
 if(fConstrainablePart && !constrainable)  return kFALSE;
 
 // Pt
 float pt = pFlowTrack->Pt();
 if(fPtPart[1] > fPtPart[0] &&  (pt < fPtPart[0] || pt > fPtPart[1])) return kFALSE;
 
 // P
 float totalp = pFlowTrack->P();
 if(fPPart[1] > fPPart[0] && (totalp < fPPart[0] || totalp > fPPart[1])) return kFALSE;
 
 // Eta
 float eta = pFlowTrack->Eta();
 if(fEtaPart[1] > fEtaPart[0] && (eta < fEtaPart[0] || eta > fEtaPart[1])) return kFALSE;
 
 // |Eta|
 float absEta = TMath::Abs(pFlowTrack->Eta());
 if(fEtaAbsPart[1] > fEtaAbsPart[0] && (absEta < fEtaAbsPart[0] || absEta > fEtaAbsPart[1])) return kFALSE;
 
 // Fit Points (TPC)
 int fitPts = pFlowTrack->FitPtsTPC();
 if(fFitPtsPart[1] > fFitPtsPart[0] && (fitPts < fFitPtsPart[0] || fitPts > fFitPtsPart[1])) return kFALSE;
 
 // Fit Points over Max Points (TPC)
 int maxPts = pFlowTrack->MaxPtsTPC();
 if(maxPts) 
 { 
  float fitOverMaxPts = (float)(fitPts)/(float)maxPts ;
  if(fFitOverMaxPtsPart[1] > fFitOverMaxPtsPart[0] && (fitOverMaxPts < fFitOverMaxPtsPart[0] || fitOverMaxPts > fFitOverMaxPtsPart[1])) return kFALSE;
 }
 
 // Chi Squared (main Vertex)
 float chiSq = pFlowTrack->Chi2();
 if(fChiSqPart[1] > fChiSqPart[0] && (chiSq < fChiSqPart[0] || chiSq > fChiSqPart[1])) return kFALSE;
 
 // DCA Global
 float globdca = pFlowTrack->Dca();
 if(fDcaGlobalPart[1] > fDcaGlobalPart[0] && (globdca < fDcaGlobalPart[0] || globdca > fDcaGlobalPart[1])) return kFALSE;
 
 // DCA / error
 float dcaSigma = 1. ;
 if(pFlowTrack->TransDcaError() != 0) { dcaSigma = pFlowTrack->TransDca() / pFlowTrack->TransDcaError() ; }
 if(fDcaOverSigma[1] > fDcaOverSigma[0] && (dcaSigma < fDcaOverSigma[0] || dcaSigma > fDcaOverSigma[1])) return kFALSE;

 // Rapidity
 float y = pFlowTrack->Y();
 if(fYPart[1] > fYPart[0] && (y < fYPart[0] || y > fYPart[1])) return kFALSE;

 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::SelectPart(AliFlowV0* pFlowV0) const 
{
 // Make v0 selection for Correlation Analysis & Vn (neutral particles 
 // to correlate with the event plane).  
 // Returns kTRUE if the v0 is selected.

 // InvMass
 float mass = pFlowV0->Mass() ;
 if(fV0Mass[1]>fV0Mass[0])
 {
  float massMin = fV0Mass[0]-fV0SideBand ;
  float massMax = fV0Mass[1]+fV0SideBand ;
  if((mass < massMin) || (mass > massMax)) return kFALSE ;
 }

 // Pt
 float pt = pFlowV0->Pt();
 if(fV0Pt[1] > fV0Pt[0] &&  (pt < fV0Pt[0] || pt > fV0Pt[1])) return kFALSE;

 // P
 float totalp = pFlowV0->P();
 if(fV0P[1] > fV0P[0] && (totalp < fV0P[0] || totalp > fV0P[1])) return kFALSE;

 // Eta
 float eta = pFlowV0->Eta();
 if(fV0Eta[1] > fV0Eta[0] && (eta < fV0Eta[0] || eta > fV0Eta[1])) return kFALSE;
 
 // |Eta|
 float absEta = TMath::Abs(pFlowV0->Eta());
 if(fV0EtaAbs[1] > fV0EtaAbs[0] && (absEta < fV0EtaAbs[0] || absEta > fV0EtaAbs[1])) return kFALSE;
 
 // Chi Squared (main Vertex)
 float chiSq = pFlowV0->Chi2();
 if(fV0ChiSq[1] > fV0ChiSq[0] && (chiSq < fV0ChiSq[0] || chiSq > fV0ChiSq[1])) return kFALSE;

 // DCA Cross
 float cdca = pFlowV0->DaughtersDca() ;
 if(fV0DcaCross[1] > fV0DcaCross[0] && (cdca < fV0DcaCross[0] || cdca > fV0DcaCross[1])) return kFALSE;

 // V0 lenght
 float lenght = pFlowV0->V0Lenght() ;
 if(fV0Lenght[1] > fV0Lenght[0] && (lenght < fV0Lenght[0] || lenght > fV0Lenght[1])) return kFALSE;

 // V0 lenght
 float sigma  = pFlowV0->Sigma() ;
 if(sigma) 
 {
  float lenghtOverSigma = lenght/sigma ;
  if(fV0LenghtOverSigma[1] > fV0LenghtOverSigma[0] && (lenghtOverSigma < fV0LenghtOverSigma[0] || lenghtOverSigma > fV0LenghtOverSigma[1])) return kFALSE;
 }

 // Rapidity
 float y = pFlowV0->Y();
 if(fV0Y[1] > fV0Y[0] && (y < fV0Y[0] || y > fV0Y[1])) return kFALSE;

//  // PID (useless)
//  if(fV0Pid[0] != '\0')
//  {
//   int numPid = TMath::Abs(pFlowV0->MostLikelihoodPID()) ;
//   if((fV0Pid = "gamma") && (fMostLikelihoodPID != 22))          { return kFALSE ; }
//   else if((fV0Pid = "K0") && (fMostLikelihoodPID != 311))       { return kFALSE ; }
//   else if((fV0Pid = "K0s") && (fMostLikelihoodPID != 310))      { return kFALSE ; }
//   else if((fV0Pid = "K0l") && (fMostLikelihoodPID != 130))      { return kFALSE ; }
//   else if((fV0Pid = "Lambda0") && (fMostLikelihoodPID != 3122)) { return kFALSE ; }
//  }

 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0Part(AliFlowV0* pFlowV0) const 
{ 
 // selects v0s in the invariant Mass window
 
 float mass = pFlowV0->Mass() ;
 if(fV0Mass[1]>fV0Mass[0])
 {
  if(mass < fV0Mass[0] || mass > fV0Mass[1]) return kFALSE ; 
 } 
 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0Side(AliFlowV0* pFlowV0) const 
{ 
 // selects v0s in the sidebands of the Mass window
 
 float mass = pFlowV0->Mass() ;
 if(fV0Mass[1]>fV0Mass[0])
 { 
  float massMin = fV0Mass[0]-fV0SideBand ;
  float massMax = fV0Mass[1]+fV0SideBand ;
  if((mass < massMin) || (mass >= fV0Mass[0] && mass <= fV0Mass[1]) || (mass > massMax)) return kFALSE ;
 }
 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0sxSide(AliFlowV0* pFlowV0) const 
{ 
 // selects v0s in the left hand sideband
 
 float mass = pFlowV0->Mass() ;
 float massMin = fV0Mass[0]-fV0SideBand ;
 if((mass >= massMin) && (mass < fV0Mass[0])) 	{ return kTRUE ; }
 else 						{ return kFALSE ; }
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0dxSide(AliFlowV0* pFlowV0) const 
{ 
 // selects v0s in the right hand sideband
 
 float mass = pFlowV0->Mass() ;
 float massMax = fV0Mass[1]+fV0SideBand ;
 if((mass > fV0Mass[1]) && (mass <= massMax)) 	{ return kTRUE ; }
 else 						{ return kFALSE ; }
}
//-----------------------------------------------------------------------
void AliFlowSelection::PrintList() const 
{
 // Prints the selection criteria for correlation analysis

 cout << "#################################################################" << endl;
 cout << "# Selection List (particles correlated to the event plane):" << endl;
 cout << "# " << endl;
 if(fPidPart[0]!='\0') 
 { cout << "# P.id for particles correlated to the event plane: " << fPidPart << "  " << endl ; }
 if(fPtPart[1]>fPtPart[0]) 
 { cout << "# Pt for particles correlated to the event plane: " << fPtPart[0] << " to " << fPtPart[1] << " GeV/c" << endl ; }
 if(fPPart[1]>fPPart[0]) 
 { cout << "# P for particles correlated to the event plane: " << fPPart[0] << " to " << fPPart[1] << " GeV/c" << endl ; }
 if(fEtaPart[1]>fEtaPart[0]) 
 { cout << "# Eta for particles correlated to the event plane: " << fEtaPart[0] << " to " << fEtaPart[1] << endl ; }
 if(fEtaAbsPart[1]>fEtaAbsPart[0]) 
 { cout << "# |Eta| for V0s correlated to the event plane: " << fEtaAbsPart[0] << " to " << fEtaAbsPart[1] << endl ; }
 if(fYPart[1]>fYPart[0]) 
 { cout << "# Y for particles correlated to the event plane: " << fYPart[0] << " to " << fYPart[1] << endl ; }
 if(fFitPtsPart[1]>fFitPtsPart[0]) 
 { cout << "# Fit Points for particles correlated to the event plane: " << fFitPtsPart[0] << " to " << fFitPtsPart[1] << endl ; }
 if(fDedxPtsPart[1]>fDedxPtsPart[0]) 
 { cout << "# Dedx Points for particles correlated to the event plane: " << fDedxPtsPart[0] << " to " << fDedxPtsPart[1] << endl ; }
 if(fFitOverMaxPtsPart[1]>fFitOverMaxPtsPart[0]) 
 { cout << "# Fit/Max Points for particles correlated to the event plane: " << fFitOverMaxPtsPart[0] << " to " << fFitOverMaxPtsPart[1] << endl ; }
 if(fChiSqPart[1]>fChiSqPart[0]) 
 { cout << "# Chi2 for particles correlated to the event plane: " << fChiSqPart[0] << " to " << fChiSqPart[1] << endl ; }
 if(fDcaGlobalPart[1]>fDcaGlobalPart[0]) 
 { cout << "# Global Dca for particles correlated with the event plane: " << fDcaGlobalPart[0] << " to " << fDcaGlobalPart[1] << endl ; }
 if(fDcaOverSigma[1]>fDcaOverSigma[0])
 { cout << "# Transverse Dca / Sigma for particles correlated with the event plane: " << fDcaOverSigma[0] << " to " << fDcaOverSigma[1] << endl ; }
 if(fConstrainablePart) 	 { cout << "# Constrainability cut:   Constrained  Tracks " << endl ; }
 cout << "#################################################################" << endl;
}
//-----------------------------------------------------------------------
void AliFlowSelection::PrintV0List() const 
{
 // Prints the selection criteria for V0s

 cout << "#################################################################" << endl;
 cout << "# Selection List (V0s correlated with the event plane):" << endl;
 cout << "# " << endl;
 // if(fV0Pid[0]!='\0') { cout << "# P.id for V0s correlated to the event plane: " << fV0Pid << "  " << endl ; }
 if(fV0Mass[1]>fV0Mass[0]) 
 {
  if(!fV0SideBand)
  { cout << "# Invariant Mass for V0s correlated to the event plane: " << fV0Mass[0] << " to " << fV0Mass[1] << endl ; }
  else
  { cout << "# Invariant Mass for V0-SideBands correlated to the event plane: " << fV0Mass[0]-fV0SideBand << " to " << fV0Mass[0]  << " & " << fV0Mass[1] << " to " << fV0Mass[1]+fV0SideBand << endl ; }
 }
 if(fV0Pt[1]>fV0Pt[0]) 
 { cout << "# Pt for V0s correlated to the event plane: " << fV0Pt[0] << " to " << fV0Pt[1] << " GeV/c" << endl ; }
 if(fV0P[1]>fV0P[0]) 
 { cout << "# P for V0s correlated to the event plane: " << fV0P[0] << " to " << fV0P[1] << " GeV/c" << endl ; }
 if(fV0Eta[1]>fV0Eta[0]) 
 { cout << "# Eta for V0s correlated to the event plane: " << fV0Eta[0] << " to " << fV0Eta[1] << endl ; }
 if(fV0EtaAbs[1]>fV0EtaAbs[0]) 
 { cout << "# |Eta| for V0s correlated to the event plane: " << fV0EtaAbs[0] << " to " << fV0EtaAbs[1] << endl ; }
 if(fV0Y[1]>fV0Y[0]) 
 { cout << "# Y for V0s correlated to the event plane: " << fV0Y[0] << " to " << fV0Y[1] << endl ; }
 if(fV0ChiSq[1]>fV0ChiSq[0]) 
 { cout << "# Chi2 for V0s correlated to the event plane: " << fV0ChiSq[0] << " to " << fV0ChiSq[1] << endl ; }
 if(fV0DcaCross[1]>fV0DcaCross[0]) 
 { cout << "# Closest approach between the daughter tracks for V0s correlated to the event plane: " << fV0DcaCross[0] << " to " << fV0DcaCross[1] << endl ; }
 if(fV0Lenght[1]>fV0Lenght[0]) 
 { cout << "# ... for V0s correlated to the event plane: " << fV0Lenght[0] << " to " << fV0Lenght[1] << endl ; }
 if(fV0LenghtOverSigma[1]>fV0LenghtOverSigma[0]) 
 { cout << "# ... for V0s correlated to the event plane: " << fV0LenghtOverSigma[0] << " to " << fV0LenghtOverSigma[1] << endl ; }
 cout << "#################################################################" << endl;
}
//-----------------------------------------------------------------------
void AliFlowSelection::PrintSelectionList() const
{
 // Prints the list of selection cuts for RP determination
 
 if(fCent>=0)      
 { 
  cout << "#######################################################" << endl;
  cout << "# Event centrality: " << endl ; 
  float lowC, hiC ;
  if(fCent==0)      { lowC=0 ; hiC=AliFlowConstants::fgCentNorm[0] * AliFlowConstants::fgMaxMult ; }
  else if(fCent>8)  { lowC=AliFlowConstants::fgCentNorm[8] * AliFlowConstants::fgMaxMult ; hiC=99999 ; }
  else              { lowC=AliFlowConstants::fgCentNorm[fCent-1] * AliFlowConstants::fgMaxMult ; hiC=AliFlowConstants::fgCentNorm[fCent] * AliFlowConstants::fgMaxMult ; }
  cout << "#  - Centrality Class = " << fCent << " ( " << (int)lowC << " < mult < " << (int)hiC << " ) . " << endl ; 
 }
 else      
 { 
  cout << "#######################################################" << endl;
  cout << "# All centralities " << endl ; 
 }

 cout << "#######################################################" << endl;
 cout << "# Tracks used for the event plane: " << endl ; 
 cout << "#  - Selection[0]    (for all " << AliFlowConstants::kHars << " Harmonics) :  " << endl ; 
 cout << "#   NO CUTS " << endl ; 
 cout << "#  - Selection[1+]   (for all " << AliFlowConstants::kHars << " Harmonics) : " << endl ; 
 if(Pid()[0] != '\0') 			{ cout << "#   Particle ID =  " << Pid() << endl ; } 
 if(ConstrainCut())     		{ cout << "#   Constrainable Tracks " << endl ; }
 if(DcaGlobalCutHi()>DcaGlobalCutLo()) 	{ cout << "#   Global Dca Tpc cuts =  " << DcaGlobalCutLo() << " , " << DcaGlobalCutHi() << endl ; }
 for (int k = 1; k < AliFlowConstants::kSels; k++) 
 {
  for (int j = 0; j < AliFlowConstants::kHars ; j++) 
  {
   cout << "#  - Selection[" << k << "] , Harmonic[" << j+1 << "] :" << endl ;
   if(NhitsCut(k)) 			{ cout << "#   Minimum TPC hits =  " << NhitsCut(k) << endl ; }
   if(EtaCutHi(j,k)>EtaCutLo(j,k)) 	{ cout << "#   TMath::Abs(Eta) cuts =  " << EtaCutLo(j,k) << " , "  << EtaCutHi(j,k) << endl ; }
   if(PtCutHi(j,k)>PtCutLo(j,k))   	{ cout << "#   Pt cuts = " << PtCutLo(j,k) << " , "  << PtCutHi(j,k) << endl ; }
  }
 }
 cout << "#######################################################" << endl;  
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetHarmonic(const Int_t& harN) 
{
 // sets the Harmonic #

 if(harN < 0 || harN >= AliFlowConstants::kHars) 
 {
  cout << "### Harmonic " << harN << " not valid" << endl;
  fHarmonic = 0;
 } 
 else { fHarmonic = harN; } 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetSelection(const Int_t& selN) 
{
 // sets the Selection #

 if(selN < 0 || selN >= AliFlowConstants::kSels) 
 {
  cout << "### Selection " << selN << " not valid" << endl;
  fSelection = 0;
 } 
 else { fSelection = selN; } 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetSubevent(const Int_t& subN) 
{
 // sets the Sub-Event # (-1 for the full-event)

 if(subN < -1 || subN > AliFlowConstants::kSubs) 
 {
  cout << "### Subevent " << subN << " not valid" << endl;
  fSubevent = -1;
 } 
 else { fSubevent = subN; } 
}
//-----------------------------------------------------------------------
Float_t AliFlowSelection::PtMaxPart() const		       
{
 // Returns the upper pT cut for particle used in correlation analysis
 
 if(fPtPart[1]>fPtPart[0]) { return fPtPart[1] ; } 
 else { return 0. ; } 
}
//-----------------------------------------------------------------------
