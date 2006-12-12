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
#include "AliFlowConstants.h"
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

// - 1st selection (both harmonic) is Disabled ! - the 2 following are the same (basic cuts, no P.id)
Float_t  AliFlowSelection::fEtaTpcCuts[2][Flow::nHars][Flow::nSels] = {{{1.0,0.0},{1.0,0.0}},{{0.0,2.1},{0.0,2.1}}} ;
Float_t  AliFlowSelection::fPtTpcCuts[2][Flow::nHars][Flow::nSels]  = {{{1.0,0.1},{1.0,0.1}},{{0.0,9.0},{0.0,9.0}}} ;
Float_t  AliFlowSelection::fDcaGlobalCuts[2] 	 = { 0. , 1. } ;
Char_t   AliFlowSelection::fPid[10]              = { '\0' } ;
Bool_t   AliFlowSelection::fConstrainable        = kTRUE ;
Int_t    AliFlowSelection::fTPChits[Flow::nSels] = { 0 , 1 } ;

ClassImp(AliFlowSelection)
//-----------------------------------------------------------------------
AliFlowSelection::AliFlowSelection()
{
 // Default constructor: when initialized all selection cuts are disabled (lo>hi).

 fJustLoopConstrainable = kFALSE ;  	// = kTRUE
 // -
 fCent		       = -1 ; 
 fRun		       = -1 ; 
 // -
 fPidPart[0] = '\0' ;
 fConstrainablePart    = kFALSE ; 	// = kTRUE
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
 fYPart[0]	       = 1 ; 		// = 0.  
 fYPart[1]	       = 0 ; 		// = 0.  
 // -
 fV0Pid[10] = '\0' ;
 fV0SideBand 	       = 0. ;
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
 // -
 fPtBinsPart 	       = Flow::nPtBinsPart ;  
 // -		        
 fHarmonic  = -1 ;   // harmonic
 fSelection = -1 ;   // selection
 fSubevent  = -1 ;   // sub-event
}
//-----------------------------------------------------------------------
AliFlowSelection::~AliFlowSelection() 
{
 // default destructor (dummy)
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowEvent* pFlowEvent) 
{
 // Returns kTRUE if the event is selected. 

 if((fCent>=0) && ((Int_t)pFlowEvent->Centrality() != fCent)) { return kFALSE ; }
 //if(pFlowEvent->RunID() != fRun) { return kFALSE ; }
 //if(pFlowEvent->... != ...) { return kFALSE ; }
 
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowTrack* pFlowTrack)
{
 // Selects particles for event plane determination.
 // Returns kTRUE if the track is selected.
 
 if(!pFlowTrack->Select(fHarmonic, fSelection, fSubevent)) { return kFALSE ; }
 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::Select(AliFlowV0* pFlowV0) 
{
 // Returns kTRUE if the v0 is selected. (dummy) 

 return kTRUE ;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::SelectPart(AliFlowTrack* pFlowTrack) 
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
 float pidProb = pFlowTrack->MostLikelihoodProb() ;
 if (fPidProbPart[1] > fPidProbPart[0] &&  (pidProb < fPidProbPart[0] || pidProb > fPidProbPart[1])) return kFALSE;
 
 // Constrainable
 bool constrainable = pFlowTrack->IsConstrainable() ;
 if (fConstrainablePart && !constrainable)  return kFALSE;
 
 // Pt
 float pt = pFlowTrack->Pt();
 if (fPtPart[1] > fPtPart[0] &&  (pt < fPtPart[0] || pt > fPtPart[1])) return kFALSE;
 
 // P
 float totalp = pFlowTrack->P();
 if (fPPart[1] > fPPart[0] && (totalp < fPPart[0] || totalp > fPPart[1])) return kFALSE;
 
 // Eta
 float eta = pFlowTrack->Eta();
 if (fEtaPart[1] > fEtaPart[0] && (eta < fEtaPart[0] || eta > fEtaPart[1])) return kFALSE;
 
 // |Eta|
 float absEta = TMath::Abs(pFlowTrack->Eta());
 if (fEtaAbsPart[1] > fEtaAbsPart[0] && (absEta < fEtaAbsPart[0] || absEta > fEtaAbsPart[1])) return kFALSE;
 
 // Fit Points (TPC)
 int fitPts = pFlowTrack->FitPtsTPC();
 if (fFitPtsPart[1] > fFitPtsPart[0] && (fitPts < fFitPtsPart[0] || fitPts > fFitPtsPart[1])) return kFALSE;
 
 // Fit Points over Max Points (TPC)
 int maxPts = pFlowTrack->MaxPtsTPC();
 if(maxPts) 
 { 
  float fitOverMaxPts = (float)(fitPts)/(float)maxPts ;
  if (fFitOverMaxPtsPart[1] > fFitOverMaxPtsPart[0] && (fitOverMaxPts < fFitOverMaxPtsPart[0] || fitOverMaxPts > fFitOverMaxPtsPart[1])) return kFALSE;
 }
 
 // Chi Squared (main Vertex)
 float chiSq = pFlowTrack->Chi2();
 if (fChiSqPart[1] > fChiSqPart[0] && (chiSq < fChiSqPart[0] || chiSq > fChiSqPart[1])) return kFALSE;
 
 // DCA Global
 float globdca = pFlowTrack->Dca();
 if (fDcaGlobalPart[1] > fDcaGlobalPart[0] && (globdca < fDcaGlobalPart[0] || globdca > fDcaGlobalPart[1])) return kFALSE;

 // Rapidity
 float Y = pFlowTrack->Y();
 if (fYPart[1] > fYPart[0] && (Y < fYPart[0] || Y > fYPart[1])) return kFALSE;

 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t AliFlowSelection::SelectPart(AliFlowV0* pFlowV0) 
{
 // Make v0 selection for Correlation Analysis & Vn (neutral particles 
 // to correlate with the event plane).  
 // Returns kTRUE if the v0 is selected.

 // PID
 if(fV0Pid[0] != '\0') 
 {
  if(strstr(fV0Pid, '\0')!=0)
  {
   int charge = pFlowV0->Charge();
   if(strcmp("0", fV0Pid)==0 && charge != 0) return kFALSE;
  } 
  else 
  {
   const Char_t* pid = pFlowV0->Pid() ;
   if(strstr(pid, fV0Pid)==0) return kFALSE;
  }
 }

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
 if (fV0Pt[1] > fV0Pt[0] &&  (pt < fV0Pt[0] || pt > fV0Pt[1])) return kFALSE;

 // P
 float totalp = pFlowV0->P();
 if (fV0P[1] > fV0P[0] && (totalp < fV0P[0] || totalp > fV0P[1])) return kFALSE;

 // Eta
 float eta = pFlowV0->Eta();
 if (fV0Eta[1] > fV0Eta[0] && (eta < fV0Eta[0] || eta > fV0Eta[1])) return kFALSE;
 
 // |Eta|
 float absEta = TMath::Abs(pFlowV0->Eta());
 if (fV0EtaAbs[1] > fV0EtaAbs[0] && (absEta < fV0EtaAbs[0] || absEta > fV0EtaAbs[1])) return kFALSE;
 
 // Chi Squared (main Vertex)
 float chiSq = pFlowV0->Chi2();
 if (fV0ChiSq[1] > fV0ChiSq[0] && (chiSq < fV0ChiSq[0] || chiSq > fV0ChiSq[1])) return kFALSE;

 // DCA Cross
 float cdca = pFlowV0->CrossDca() ;
 if (fV0DcaCross[1] > fV0DcaCross[0] && (cdca < fV0DcaCross[0] || cdca > fV0DcaCross[1])) return kFALSE;

 // V0 lenght
 float lenght = pFlowV0->V0Lenght() ;
 if (fV0Lenght[1] > fV0Lenght[0] && (lenght < fV0Lenght[0] || lenght > fV0Lenght[1])) return kFALSE;

 // V0 lenght
 float sigma  = pFlowV0->Sigma() ;
 if(sigma) 
 {
  float lenghtOverSigma = lenght/sigma ;
  if (fV0LenghtOverSigma[1] > fV0LenghtOverSigma[0] && (lenghtOverSigma < fV0LenghtOverSigma[0] || lenghtOverSigma > fV0LenghtOverSigma[1])) return kFALSE;
 }

 // Rapidity
 float Y = pFlowV0->Y();
 if (fV0Y[1] > fV0Y[0] && (Y < fV0Y[0] || Y > fV0Y[1])) return kFALSE;

 return kTRUE;
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0Part(AliFlowV0* pFlowV0)
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
Bool_t  AliFlowSelection::SelectV0Side(AliFlowV0* pFlowV0)
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
Bool_t  AliFlowSelection::SelectV0sxSide(AliFlowV0* pFlowV0)
{ 
 // selects v0s in the left hand sideband
 float mass = pFlowV0->Mass() ;
 float massMin = fV0Mass[0]-fV0SideBand ;
 if((mass >= massMin) && (mass < fV0Mass[0])) 	{ return kTRUE ; }
 else 						{ return kFALSE ; }
}
//-----------------------------------------------------------------------
Bool_t  AliFlowSelection::SelectV0dxSide(AliFlowV0* pFlowV0)
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
 { cout << "# P.id for particles correlated to the event plane: " << fPidPart << "  " <<endl ; }
 if(fPtPart[1]>fPtPart[0]) 
 { cout << "# Pt for particles correlated to the event plane: " << fPtPart[0] << " to " << fPtPart[1] << " GeV/c" <<endl ; }
 if(fPPart[1]>fPPart[0]) 
 { cout << "# P for particles correlated to the event plane: " << fPPart[0] << " to " << fPPart[1] << " GeV/c" <<endl ; }
 if(fEtaPart[1]>fEtaPart[0]) 
 { cout << "# Eta for particles correlated to the event plane: " << fEtaPart[0] << " to " << fEtaPart[1] <<endl ; }
 if(fEtaAbsPart[1]>fEtaAbsPart[0]) 
 { cout << "# |Eta| for V0s correlated to the event plane: " << fEtaAbsPart[0] << " to " << fEtaAbsPart[1] <<endl ; }
 if(fYPart[1]>fYPart[0]) 
 { cout << "# Y for particles correlated to the event plane: " << fYPart[0] << " to " << fYPart[1] <<endl ; }
 if(fFitPtsPart[1]>fFitPtsPart[0]) 
 { cout << "# Fit Points for particles correlated to the event plane: " << fFitPtsPart[0] << " to " << fFitPtsPart[1] <<endl ; }
 if(fDedxPtsPart[1]>fDedxPtsPart[0]) 
 { cout << "# Dedx Points for particles correlated to the event plane: " << fDedxPtsPart[0] << " to " << fDedxPtsPart[1] << endl ; }
 if(fFitOverMaxPtsPart[1]>fFitOverMaxPtsPart[0]) 
 { cout << "# Fit/Max Points for particles correlated to the event plane: " << fFitOverMaxPtsPart[0] << " to " << fFitOverMaxPtsPart[1] <<endl ; }
 if(fChiSqPart[1]>fChiSqPart[0]) 
 { cout << "# Chi2 for particles correlated to the event plane: " << fChiSqPart[0] << " to " << fChiSqPart[1] <<endl ; }
 if(fDcaGlobalPart[1]>fDcaGlobalPart[0]) 
 { cout << "# Global Dca for particles correlated with the event plane: " << fDcaGlobalPart[0] << " to " << fDcaGlobalPart[1] <<endl ; }
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
 if(fV0Pid[0]!='\0') 
 { cout << "# P.id for V0s correlated to the event plane: " << fV0Pid << "  " <<endl ; }
 if(fV0Mass[1]>fV0Mass[0]) 
 {
  if(!fV0SideBand)
  { cout << "# Invariant Mass for V0s correlated to the event plane: " << fV0Mass[0] << " to " << fV0Mass[1] << endl ; }
  else
  { cout << "# Invariant Mass for V0-SideBands correlated to the event plane: " << fV0Mass[0]-fV0SideBand << " to " << fV0Mass[0]  << " & " << fV0Mass[1] << " to " << fV0Mass[1]+fV0SideBand << endl ; }
 }
 if(fV0Pt[1]>fV0Pt[0]) 
 { cout << "# Pt for V0s correlated to the event plane: " << fV0Pt[0] << " to " << fV0Pt[1] << " GeV/c" <<endl ; }
 if(fV0P[1]>fV0P[0]) 
 { cout << "# P for V0s correlated to the event plane: " << fV0P[0] << " to " << fV0P[1] << " GeV/c" <<endl ; }
 if(fV0Eta[1]>fV0Eta[0]) 
 { cout << "# Eta for V0s correlated to the event plane: " << fV0Eta[0] << " to " << fV0Eta[1] <<endl ; }
 if(fV0EtaAbs[1]>fV0EtaAbs[0]) 
 { cout << "# |Eta| for V0s correlated to the event plane: " << fV0EtaAbs[0] << " to " << fV0EtaAbs[1] <<endl ; }
 if(fV0Y[1]>fV0Y[0]) 
 { cout << "# Y for V0s correlated to the event plane: " << fV0Y[0] << " to " << fV0Y[1] <<endl ; }
 if(fV0ChiSq[1]>fV0ChiSq[0]) 
 { cout << "# Chi2 for V0s correlated to the event plane: " << fV0ChiSq[0] << " to " << fV0ChiSq[1] <<endl ; }
 if(fV0DcaCross[1]>fV0DcaCross[0]) 
 { cout << "# Closest approach between the daughter tracks for V0s correlated to the event plane: " << fV0DcaCross[0] << " to " << fV0DcaCross[1] <<endl ; }
 if(fV0Lenght[1]>fV0Lenght[0]) 
 { cout << "# ... for V0s correlated to the event plane: " << fV0Lenght[0] << " to " << fV0Lenght[1] <<endl ; }
 if(fV0LenghtOverSigma[1]>fV0LenghtOverSigma[0]) 
 { cout << "# ... for V0s correlated to the event plane: " << fV0LenghtOverSigma[0] << " to " << fV0LenghtOverSigma[1] <<endl ; }
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
  if(fCent==0)      { lowC=0 ; hiC=Flow::fCentNorm[0] * Flow::fMaxMult ; }
  else if(fCent>8)  { lowC=Flow::fCentNorm[8] * Flow::fMaxMult ; hiC=99999 ; }
  else              { lowC=Flow::fCentNorm[fCent-1] * Flow::fMaxMult ; hiC=Flow::fCentNorm[fCent] * Flow::fMaxMult ; }
  cout << "#  - Centrality Class = " << fCent << " ( " << (int)lowC << " < mult < " << (int)hiC << " ) . " << endl ; 
 }
 else      
 { 
  cout << "#######################################################" << endl;
  cout << "# All centralities " << endl ; 
 }

 cout << "#######################################################" << endl;
 cout << "# Tracks used for the event plane: " << endl ; 
 cout << "#  - Selection[0]    (for all " << Flow::nHars << " Harmonics) :  " << endl ; 
 cout << "#   NO CUTS " << endl ; 
 cout << "#  - Selection[1+]   (for all " << Flow::nHars << " Harmonics) : " << endl ; 
 if(Pid()[0] != '\0') 			{ cout << "#   Particle ID =  " << Pid() << endl ; } 
 if(ConstrainCut())     		{ cout << "#   Constrainable Tracks " << endl ; }
 if(DcaGlobalCutHi()>DcaGlobalCutLo()) 	{ cout << "#   Global Dca Tpc cuts =  " << DcaGlobalCutLo() << " , " << DcaGlobalCutHi() << endl ; }
 for (int k = 1; k < Flow::nSels; k++) 
 {
  for (int j = 0; j < Flow::nHars ; j++) 
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
 if (harN < 0 || harN >= Flow::nHars) 
 {
  cout << "### Harmonic " << harN << " not valid" << endl;
  fHarmonic = 0;
 } 
 else { fHarmonic = harN; } 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetSelection(const Int_t& selN) 
{
 if (selN < 0 || selN >= Flow::nSels) 
 {
  cout << "### Selection " << selN << " not valid" << endl;
  fSelection = 0;
 } 
 else { fSelection = selN; } 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetSubevent(const Int_t& subN) 
{
 if (subN < -1 || subN > Flow::nSubs) 
 {
  cout << "### Subevent " << subN << " not valid" << endl;
  fSubevent = -1;
 } 
 else { fSubevent = subN; } 
}
//-----------------------------------------------------------------------
Float_t AliFlowSelection::PtMaxPart() const			       
{ 
 if(fPtPart[1]>fPtPart[0]) { return fPtPart[1] ; } 
 else { return 0. ; } 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetEtaCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN) 
{ 
 fEtaTpcCuts[0][harN][selN] = lo ; 
 fEtaTpcCuts[1][harN][selN] = hi ; 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetPtCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN)  
{ 
 fPtTpcCuts[0][harN][selN] = lo  ; 
 fPtTpcCuts[1][harN][selN] = hi ; 
}
//-----------------------------------------------------------------------
void AliFlowSelection::SetV0SideBands() 
{  
 Float_t massInterval = fV0Mass[1] - fV0Mass[0] ;
 SetV0SideBands(TMath::Abs(massInterval/2)) ; 
}			
//-----------------------------------------------------------------------
void AliFlowSelection::SetV0SideBands(Float_t sb)           	       
{ 
 fV0SideBand = sb ; 
}			
//-----------------------------------------------------------------------
void    AliFlowSelection::SetCentralityCut(Int_t cent) 		       { fCent = cent ; }
void    AliFlowSelection::SetRunIdCut(Int_t run)  		       { fRun = run ; }
void 	AliFlowSelection::SetDcaGlobalCut(Float_t lo, Float_t hi)      { fDcaGlobalCuts[0] = lo ; fDcaGlobalCuts[1] = hi ; }
void 	AliFlowSelection::SetPidCut(const Char_t* pid)  	       { strncpy(fPid, pid, 9) ; fPid[9] = '\0' ; }
void 	AliFlowSelection::SetConstrainCut(Bool_t tf)		       { fConstrainable = tf ; }
void 	AliFlowSelection::SetNhitsCut(Int_t hits,Int_t selN)	       { fTPChits[selN] = hits; }
void	AliFlowSelection::SetPidPart(const Char_t* pid) 	       { strncpy(fPidPart, pid, 9); fPidPart[9] = '\0'; }
void	AliFlowSelection::SetPidProbPart(Float_t lo, Float_t hi)       { fPidProbPart[0] = lo ; fPidProbPart[1] = hi; }
void	AliFlowSelection::SetPtPart(Float_t lo, Float_t hi)	       { fPtPart[0] = lo; fPtPart[1] = hi; }
void	AliFlowSelection::SetPtBinsPart(Int_t bins)		       { fPtBinsPart = bins; }
void	AliFlowSelection::SetPPart(Float_t lo, Float_t hi)	       { fPPart[0] = lo; fPPart[1] = hi; }
void	AliFlowSelection::SetEtaPart(Float_t lo, Float_t hi)	       { fEtaPart[0] = lo; fEtaPart[1] = hi; }
void    AliFlowSelection::SetEtaAbsPart(Float_t lo, Float_t hi)        { fEtaAbsPart[0] = TMath::Abs(lo); fEtaAbsPart[1] = TMath::Abs(hi); }
void	AliFlowSelection::SetYPart(Float_t lo, Float_t hi)	       { fYPart[0] = lo; fYPart[1] = hi; }
void	AliFlowSelection::SetFitPtsPart(Int_t lo, Int_t hi)	       { fFitPtsPart[0] = lo; fFitPtsPart[1] = hi; }
void	AliFlowSelection::SetDedxPtsPart(Int_t lo, Int_t hi)	       { fDedxPtsPart[0] = lo; fDedxPtsPart[1] = hi; }
void	AliFlowSelection::SetFitOverMaxPtsPart(Float_t lo, Float_t hi) { fFitOverMaxPtsPart[0] = lo; fFitOverMaxPtsPart[1] = hi; }
void	AliFlowSelection::SetChiSqPart(Float_t lo, Float_t hi)         { fChiSqPart[0] = lo; fChiSqPart[1] = hi; }
void	AliFlowSelection::SetDcaGlobalPart(Float_t lo, Float_t hi)     { fDcaGlobalPart[0] = lo; fDcaGlobalPart[1] = hi; }
void    AliFlowSelection::SetConstrainablePart(Bool_t constr)          { fConstrainablePart = constr ; }
void	AliFlowSelection::SetV0Pid(const Char_t* pid)  	       	       { strncpy(fV0Pid, pid, 9) ; fV0Pid[9] = '\0' ; }			
void    AliFlowSelection::SetV0Mass(Float_t lo, Float_t hi)            { fV0Mass[0] = lo ; fV0Mass[1] = hi; }
void	AliFlowSelection::SetV0Pt(Float_t lo, Float_t hi)              { fV0Pt[0] = lo ; fV0Pt[1] = hi; }
void	AliFlowSelection::SetV0P(Float_t lo, Float_t hi)               { fV0P[0] = lo ; fV0P[1] = hi; }
void	AliFlowSelection::SetV0Eta(Float_t lo, Float_t hi)             { fV0Eta[0] = lo ; fV0Eta[1] = hi; }
void	AliFlowSelection::SetV0EtaAbs(Float_t lo, Float_t hi)          { fV0EtaAbs[0] = lo ; fV0EtaAbs[1] = hi; }
void	AliFlowSelection::SetV0Y(Float_t lo, Float_t hi)               { fV0Y[0] = lo ; fV0Y[1] = hi; }
void	AliFlowSelection::SetV0ChiSqPart(Float_t lo, Float_t hi)       { fV0ChiSq[0] = lo ; fV0ChiSq[1] = hi; }
void	AliFlowSelection::SetV0Lenght(Float_t lo, Float_t hi)          { fV0Lenght[0] = lo ; fV0Lenght[1] = hi; }
void	AliFlowSelection::SetV0DcaCross(Float_t lo, Float_t hi)        { fV0DcaCross[0] = lo ; fV0DcaCross[1] = hi; }
void	AliFlowSelection::SetV0LenghtOverSigma(Float_t lo, Float_t hi) { fV0LenghtOverSigma[0] = lo ; fV0LenghtOverSigma[1] = hi; }
//-----------------------------------------------------------------------
Char_t* AliFlowSelection::PidPart()	 			       { return fPidPart; }
Int_t	AliFlowSelection::PtBinsPart() const			       { return fPtBinsPart; }
Int_t	AliFlowSelection::Sel() const				       { return fSelection; }
Int_t	AliFlowSelection::Har() const				       { return fHarmonic; }
Int_t	AliFlowSelection::Sub() const				       { return fSubevent; }
Float_t AliFlowSelection::EtaCutLo(Int_t harN, Int_t selN) const       { return fEtaTpcCuts[0][harN][selN] ; }
Float_t AliFlowSelection::EtaCutHi(Int_t harN, Int_t selN) const       { return fEtaTpcCuts[1][harN][selN] ; }
Float_t AliFlowSelection::PtCutLo(Int_t harN, Int_t selN) const        { return fPtTpcCuts[0][harN][selN] ; }
Float_t AliFlowSelection::PtCutHi(Int_t harN, Int_t selN) const        { return fPtTpcCuts[1][harN][selN] ; }
Char_t* AliFlowSelection::Pid() const			               { return fPid; }
Float_t AliFlowSelection::DcaGlobalCutLo() const  	               { return fDcaGlobalCuts[0] ; }
Float_t AliFlowSelection::DcaGlobalCutHi() const  	               { return fDcaGlobalCuts[1] ; }
Bool_t  AliFlowSelection::ConstrainCut() const		   	       { return fConstrainable ; }
Int_t   AliFlowSelection::NhitsCut(Int_t selN) const		       { return fTPChits[selN] ; }
//-----------------------------------------------------------------------
Int_t   AliFlowSelection::CentralityCut() const 		       { return fCent ; }
Int_t   AliFlowSelection::RunIdCut() const 		       	       { return fRun ; }
//-----------------------------------------------------------------------
void    AliFlowSelection::SetJustLoopConstrainable()		       { fJustLoopConstrainable = kTRUE ; }
Bool_t  AliFlowSelection::JustLoopConstrainable() const		       { return fJustLoopConstrainable ; }
//-----------------------------------------------------------------------
