//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: class for selections in flow study, adapted from STAR 
// Original Authors:                Raimond Snellings & Art Poskanzer
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowSelection_h
#define AliFlowSelection_h

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "TObject.h"
#include "TVector.h"
#include "TMath.h"
#include <TROOT.h>

#include "AliFlowSelection.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowConstants.h"

class AliFlowTrack ;
class AliFlowEvent ;
class Flow ;

class AliFlowSelection : public TObject {


 public:

          AliFlowSelection();
  virtual ~AliFlowSelection();

 // Selection Methods for ... 
  Bool_t  Select(AliFlowEvent*);		       // (dummy)
  Bool_t  Select(AliFlowTrack*);		       // selection for R.P.[nSel][nHar]
  Bool_t  Select(AliFlowV0*);    		       // (dummy) 
  Bool_t  SelectPart(AliFlowTrack*);		       // track selection for Correlation Analysis
  Bool_t  SelectPart(AliFlowV0*);    		       // v0 selection for Correlation Analysis (mass window + sidebands)
  Bool_t  SelectV0Part(AliFlowV0*);    		       // v0 mass window for Correlation Analysis 
  Bool_t  SelectV0Side(AliFlowV0*);    		       // v0 sidebands for Correlation Analysis 
  Bool_t  SelectV0sxSide(AliFlowV0*);		       // selects v0s in the left hand sideband
  Bool_t  SelectV0dxSide(AliFlowV0*);		       // selects v0s in the right hand sideband

 // Gets (Harmonic, Selection, Sub-event)
  Int_t   Sel() const				       { return fSelection; }	  	       // Returns the Harmonic 
  Int_t   Har() const				       { return fHarmonic; }	  	       // Returns the Selection
  Int_t   Sub() const				       { return fSubevent; }	  	       // Returns the Sub-Event

 // Gets (R.P. cuts) and CutList
  Float_t EtaCutLo(Int_t harN, Int_t selN) const       { return fEtaTpcCuts[0][harN][selN] ; } // Returns lower eta cut for R.P.[harN][selN] calculation (absolute values)
  Float_t EtaCutHi(Int_t harN, Int_t selN) const       { return fEtaTpcCuts[1][harN][selN] ; } // Returns upper eta cut for R.P.[harN][selN] calculation (absolute values)
  Float_t PtCutLo(Int_t harN, Int_t selN) const        { return fPtTpcCuts[0][harN][selN] ; }  // Returns lower pT cut for R.P.[harN][selN] calculation 
  Float_t PtCutHi(Int_t harN, Int_t selN) const        { return fPtTpcCuts[1][harN][selN] ; }  // Returns upper pT cut for R.P.[harN][selN] calculation
  Float_t DcaGlobalCutLo() const		       { return fDcaGlobalCuts[0] ; }	       // Returns lower DCA cut for R.P. calculation
  Float_t DcaGlobalCutHi() const		       { return fDcaGlobalCuts[1] ; }	       // Returns upper DCA cut for R.P. calculation
  Bool_t  ConstrainCut() const  		       { return fConstrainable ; }	       // Returns kTRUE/kFalse if the cut over un-constrainable tracks is enabled
  Int_t   NhitsCut(Int_t selN) const		       { return fTPChits[selN] ; }	       // Returns the minimum number of TPC hits for R.P.[selN] calculation
  Char_t* Pid() const				       { return fPid; } 		       // Returns particle specie used in R.P. calculation

 // Gets (Event cuts)
  Int_t   CentralityCut() const			       { return fCent ; } 		       // Returns Event Centrality class
  Int_t   RunIdCut() const			       { return fRun ; }		       // Returns Run number 

 // Gets (Correlation analysis cuts of tracks & V0s)
  Char_t* PidPart()				       { return fPidPart; }		       // Returns selected particle species wrt Reaction Plane
  Int_t	  PtBinsPart() const			       { return fPtBinsPart; }		       // Returns N. of pT binning	       
  Float_t PtMaxPart() const ;								       // Returns the max pT for evt.plane calc.

 // Cuts list
  void    PrintList() const ;								       // Prints the tracks cut-list (for correlation analysis)
  void    PrintSelectionList() const ;							       // Prints a summary of the selection criteria (for RP determination)
  void    PrintV0List() const ;								       // Prints the v0s cut-list (for correlation analysis)

 // Harmonic & Selection set (R.P.)
  void    SetHarmonic(const Int_t&);							       // Sets the Harmonic
  void    SetSelection(const Int_t&);							       // Sets the Selection
  void    SetSubevent(const Int_t&);							       // Sets the Sub-Event
  
 // Cuts set (Reaction Plane)
  static void  SetPidCut(const Char_t* pid);				   		       // Sets the particle specie used in R.P. calculation
  static void  SetEtaCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN);  		       // Sets |eta| cut for R.P.[harN][selN] calculation
  static void  SetPtCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN);   		       // Sets pT cut for R.P.[harN][selN] calculation
  static void  SetDcaGlobalCut(Float_t lo, Float_t hi); 		   		       // Sets DCA cut for R.P. calculation
  static void  SetConstrainCut(Bool_t tf = kTRUE) ;			   		       // Sets the cut over un-constrainable tracks
  static void  SetNhitsCut(Int_t hits, Int_t selN) ; 			   		       // Sets the minimum number of TPC hits for R.P.[selN] calculation

 // Sets (Event cuts)
  void    SetCentralityCut(Int_t cent)		       { fCent = cent ; }		       // Sets Event Centrality class
  void    SetRunIdCut(Int_t run)  		       { fRun = run ; } 		       // Sets Run number 

 // Cuts set (correlation analysis cuts of tracks & V0s)
  void    SetPtBinsPart(Int_t bins)		       { fPtBinsPart = bins; }						     // Sets N. of bins from fPtPart[0] to fPtPart[1]		     

  void    SetPidPart(const Char_t* pid) 	       { strncpy(fPidPart, pid, 9); fPidPart[9] = '\0'; }		     // Sets PID for particles wrt Reaction plane 
  void    SetPidProbPart(Float_t lo, Float_t hi)       { fPidProbPart[0] = lo ; fPidProbPart[1] = hi; } 		     // Sets PID probability for particles wrt Reaction plane
  void    SetPtPart(Float_t lo, Float_t hi)	       { fPtPart[0] = lo; fPtPart[1] = hi; }				     // Sets pT for particles wrt Reaction plane 
  void    SetPPart(Float_t lo, Float_t hi)	       { fPPart[0] = lo; fPPart[1] = hi; }				     // Sets Momentum for particles wrt Reaction plane 
  void    SetEtaPart(Float_t lo, Float_t hi)	       { fEtaPart[0] = lo; fEtaPart[1] = hi; }				     // Sets Eta for particles wrt Reaction plane    
  void    SetEtaAbsPart(Float_t lo, Float_t hi)        { fEtaAbsPart[0] = TMath::Abs(lo); fEtaAbsPart[1] = TMath::Abs(hi); } // Sets |Eta| for particles wrt Reaction plane	     
  void    SetYPart(Float_t lo, Float_t hi)	       { fYPart[0] = lo; fYPart[1] = hi; }				     // Sets Rapidity for particles (with sign.) wrt Reaction plane 
  void    SetFitPtsPart(Int_t lo, Int_t hi)	       { fFitPtsPart[0] = lo; fFitPtsPart[1] = hi; }			     // Sets FitPoints for particles wrt Reaction plane  
  void    SetDedxPtsPart(Int_t lo, Int_t hi)	       { fDedxPtsPart[0] = lo; fDedxPtsPart[1] = hi; }			     // Sets dE/dx for particles wrt Reaction plane	     
  void    SetFitOverMaxPtsPart(Float_t lo, Float_t hi) { fFitOverMaxPtsPart[0] = lo; fFitOverMaxPtsPart[1] = hi; }	     // Sets FitPoints/MaxPoints for particles wrt Reaction plane 
  void    SetChiSqPart(Float_t lo, Float_t hi)         { fChiSqPart[0] = lo; fChiSqPart[1] = hi; }			     // Sets Chi^2 for particles wrt Reaction plane	     
  void    SetDcaGlobalPart(Float_t lo, Float_t hi)     { fDcaGlobalPart[0] = lo; fDcaGlobalPart[1] = hi; }		     // Sets d.c.a. for particles wrt Reaction plane	     
  void    SetConstrainablePart(Bool_t constr)	       { fConstrainablePart = constr ; }				     // Sets constrainability for particles wrt Reaction plane      

  void    SetV0Pid(const Char_t* pid)		       { strncpy(fV0Pid, pid, 9) ; fV0Pid[9] = '\0' ; } 		     // Sets PID for v0 wrt plane (...)
  void    SetV0Mass(Float_t lo, Float_t hi)	       { fV0Mass[0] = lo ; fV0Mass[1] = hi; }				     // Sets invariant mass cut for v0 wrt plane 
  void    SetV0Pt(Float_t lo, Float_t hi)	       { fV0Pt[0] = lo ; fV0Pt[1] = hi; }				     // Sets pT for v0 wrt plane 
  void    SetV0P(Float_t lo, Float_t hi)	       { fV0P[0] = lo ; fV0P[1] = hi; } 				     // Sets Momentum for v0 wrt plane 
  void    SetV0Eta(Float_t lo, Float_t hi)	       { fV0Eta[0] = lo ; fV0Eta[1] = hi; }				     // Sets Eta cut for v0 wrt plane 
  void    SetV0EtaAbs(Float_t lo, Float_t hi)	       { fV0EtaAbs[0] = lo ; fV0EtaAbs[1] = hi; }			     // Sets |Eta| cut (absolute value) for v0 wrt plane 
  void    SetV0Y(Float_t lo, Float_t hi)	       { fV0Y[0] = lo ; fV0Y[1] = hi; } 				     // Sets Rapidity for v0 wrt plane 
  void    SetV0ChiSqPart(Float_t lo, Float_t hi)       { fV0ChiSq[0] = lo ; fV0ChiSq[1] = hi; } 			     // Sets Chi^2 for v0 wrt plane 
  void    SetV0Lenght(Float_t lo, Float_t hi)	       { fV0Lenght[0] = lo ; fV0Lenght[1] = hi; }			     // Sets distance to the main vertex for v0 wrt plane
  void    SetV0DcaCross(Float_t lo, Float_t hi)        { fV0DcaCross[0] = lo ; fV0DcaCross[1] = hi; }			     // Sets distance to the main vertex in sigma units for v0 wrt plane
  void    SetV0LenghtOverSigma(Float_t lo, Float_t hi) { fV0LenghtOverSigma[0] = lo ; fV0LenghtOverSigma[1] = hi; }	     // Sets closest approach (between the 2 daughter tracks) for v0 wrt plane  
  void    SetV0SideBands()			       { SetV0SideBands(TMath::Abs((fV0Mass[1]-fV0Mass[0])/2)) ; }	     // Includes the v0 sideband analysis wrt plane	     
  void    SetV0SideBands(Float_t sb)		       { fV0SideBand = sb ; }						     // Includes the v0 sideband analysis and a width		     


 private:

 // These are just 3 integers - simple way to look at the [nHar][nSel] and [nSub] array
  Int_t   fHarmonic;                         			// harmonic
  Int_t   fSelection;                        			// selection
  Int_t   fSubevent;                         			// sub-event

  Int_t   fPtBinsPart;                       			// N. of bins in pT histograms (pT binning)

 // Event Cuts  (new)
  Int_t  fCent ;           	        			// Event Centrality class
  Int_t  fRun ;           	        			// Run number 

 // Cuts for V0 correlated to the Raction Plane (new)
  Char_t  fV0Pid[10];           	        		// PID for v0 wrt plane (...)
  Float_t fV0SideBand ;						// width of the sidebands (using the sidebands' candidates)
  Float_t fV0Mass[2] ;						// mass cut for v0 wrt plane
  Float_t fV0Pt[2];                        			// pT for v0 wrt plane 
  Float_t fV0P[2];                         			// Momentum for v0 wrt plane 
  Float_t fV0Eta[2];                       			// Eta cut for v0 wrt plane 
  Float_t fV0EtaAbs[2];                       			// |Eta| cut (absolute value) for v0 wrt plane 
  Float_t fV0Y[2];                         			// Rapidity for v0 wrt plane 
  Float_t fV0ChiSq[2];                     			// Chi^2 for v0 wrt plane 
  Float_t fV0Lenght[2];                 			// distance to the main vertex for v0 wrt plane
  Float_t fV0LenghtOverSigma[2];                 		// distance to the main vertex in sigma units for v0 wrt plane
  Float_t fV0DcaCross[2];                 			// closest approach (between the 2 daughter tracks) for v0 wrt plane

 // Cuts for Tracks that will be related to the Raction Plane (original strategy from STAR)
  Char_t  fPidPart[10];           	        		// PID for parts. wrt plane (h+, h-, pi-, pi+, pi, k+, k-, k, pr+, pr-, pr, d+, d-, d, e+, e-, e)
  Float_t fPidProbPart[2] ;                     		// probability of the most likelihood p.id. (you should specify also PidPart())
  Float_t fPtPart[2];                        			// pT for parts. wrt plane
  Float_t fPPart[2];                         			// Momentum for parts. wrt plane
  Float_t fEtaPart[2];                       			// Eta cut for parts. wrt plane
  Float_t fEtaAbsPart[2];                       		// |Eta| cut (absolute value) for parts. wrt plane
  Float_t fYPart[2];                         			// Rapidity for parts. wrt plane 
  Int_t   fFitPtsPart[2];                    			// FitPoints for parts. wrt plane
  Int_t   fDedxPtsPart[2];                   			// dE/dx for parts. wrt plane
  Float_t fFitOverMaxPtsPart[2];             			// FitPoints/MaxPoints for parts. wrt plane
  Float_t fChiSqPart[2];                     			// Chi^2 for parts. wrt plane
  Float_t fDcaGlobalPart[2];                 			// closest approach (to the main vertex) for parts. wrt plane
  Bool_t  fConstrainablePart;			    		// constrainability for parts. wrt plane 

 // Cuts for Tracks used in determining the Raction Plane (in STAR this selection was done inside the AliFlowEvent class)
  static Float_t  fEtaTpcCuts[2][Flow::nHars][Flow::nSels]; 	//! eta range (absolute values)
  static Float_t  fPtTpcCuts[2][Flow::nHars][Flow::nSels];  	//! pT range
  static Float_t  fDcaGlobalCuts[2];			    	//! DCA cuts
  static Char_t   fPid[10];				    	//! h+, h-, pi-, pi+, pi, k+, k-, k, pr+, pr-, pr, e+, e-, e, d+, d-, d
  static Int_t    fTPChits[Flow::nSels];		    	//! minimum number of TPC hits
  static Bool_t   fConstrainable;			    	//! cut un-constrainable tracks 


  ClassDef(AliFlowSelection,1)               			// macro for rootcint
}; 

#endif
