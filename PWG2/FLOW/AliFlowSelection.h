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
  Bool_t  Select(AliFlowEvent*);				// (dummy)
  Bool_t  Select(AliFlowTrack*);				// selection for R.P.[nSel][nHar]
  Bool_t  Select(AliFlowV0*);    				// (dummy) 
  Bool_t  SelectPart(AliFlowTrack*);				// track selection for Correlation Analysis
  Bool_t  SelectPart(AliFlowV0*);    				// v0 selection for Correlation Analysis (mass window + sidebands)
  Bool_t  SelectV0Part(AliFlowV0*);    				// v0 mass window for Correlation Analysis 
  Bool_t  SelectV0Side(AliFlowV0*);    				// v0 sidebands for Correlation Analysis 
  Bool_t  SelectV0sxSide(AliFlowV0*);				// selects v0s in the left hand sideband
  Bool_t  SelectV0dxSide(AliFlowV0*);				// selects v0s in the right hand sideband

 // Gets (Harmonic, Selection, Sub-event)
  Int_t   Sel() const;						// Returns the Harmonic 		      
  Int_t   Har() const;						// Returns the Selection
  Int_t   Sub() const;						// Returns the Sub-Event

 // Gets (Event cuts)
  Int_t   CentralityCut() const ;				// Returns Event Centrality class
  Int_t   RunIdCut() const ;					// Returns Run number 

 // Gets (R.P. cuts) and CutList
  Float_t EtaCutLo(Int_t harN, Int_t selN) const;		// Returns lower eta cut for R.P.[harN][selN] calculation (absolute values)
  Float_t EtaCutHi(Int_t harN, Int_t selN) const;		// Returns upper eta cut for R.P.[harN][selN] calculation (absolute values)				
  Float_t PtCutLo(Int_t harN, Int_t selN)  const;		// Returns lower pT cut for R.P.[harN][selN] calculation 
  Float_t PtCutHi(Int_t harN, Int_t selN) const ;		// Returns upper pT cut for R.P.[harN][selN] calculation
  Float_t DcaGlobalCutLo() const;				// Returns lower DCA cut for R.P. calculation
  Float_t DcaGlobalCutHi() const;				// Returns upper DCA cut for R.P. calculation
  Char_t* Pid() const;					        // Returns particle specie used in R.P. calculation
  Bool_t  ConstrainCut() const;				        // Returns kTRUE/kFalse if the cut over un-constrainable tracks is enabled
  Int_t   NhitsCut(Int_t selN) const;  			        // Returns the minimum number of TPC hits for R.P.[selN] calculation

 // Gets (correlation cuts)
  Char_t* PidPart() ;						// Returns selected particle species wrt Reaction Plane
  Float_t PtMaxPart() const ;					// Returns the max pT for evt.plane calc.
  Int_t   PtBinsPart() const ;					// Returns N. of pT binning
  void    SetPtBinsPart(const Int_t);				// Sets N. of bins from fPtPart[0] to fPtPart[1]
  void    PrintList() const ;					// Prints the tracks cut-list (for correlation analysis)
  void    PrintSelectionList() const ;				// Prints a summary of the selection criteria (for RP determination)
  void    PrintV0List() const ;					// Prints the v0s cut-list (for correlation analysis)

 // Harmonic & Selection set (R.P.)
  void    SetHarmonic(const Int_t&);				// Sets the Harmonic
  void    SetSelection(const Int_t&);				// Sets the Selection
  void    SetSubevent(const Int_t&);				// Sets the Sub-Event
  
 // Cuts set (Events)
  void    SetCentralityCut(Int_t cent) ;			// Sets Event Centrality class
  void    SetRunIdCut(Int_t run) ;				// Sets Run number 
  
 // Cuts set (Reaction Plane)
  static void  SetPidCut(const Char_t* pid);				   // Sets the particle specie used in R.P. calculation
  static void  SetEtaCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN);  // Sets |eta| cut for R.P.[harN][selN] calculation
  static void  SetPtCut(Float_t lo, Float_t hi, Int_t harN, Int_t selN);   // Sets pT cut for R.P.[harN][selN] calculation
  static void  SetDcaGlobalCut(Float_t lo, Float_t hi); 		   // Sets DCA cut for R.P. calculation
  static void  SetConstrainCut(Bool_t tf = kTRUE) ;			   // Sets the cut over un-constrainable tracks
  static void  SetNhitsCut(Int_t hits, Int_t selN) ; 			   // Sets the minimum number of TPC hits for R.P.[selN] calculation

 // Cuts set (Correlation Analysis)
  void    SetPidPart(const Char_t*);				// Sets PID for particles wrt Reaction plane 
  void	  SetPidProbPart(const Float_t, const Float_t);         // Sets PID probability for particles wrt Reaction plane
  void    SetPtPart(const Float_t, const Float_t);		// Sets pT for particles wrt Reaction plane 
  void    SetPPart(const Float_t, const Float_t);		// Sets Momentum for particles wrt Reaction plane 
  void    SetEtaPart(const Float_t, const Float_t);		// Sets Eta for particles wrt Reaction plane  	
  void    SetEtaAbsPart(const Float_t, const Float_t);		// Sets |Eta| for particles wrt Reaction plane  	
  void    SetYPart(const Float_t, const Float_t);		// Sets Rapidity for particles (with sign.) wrt Reaction plane 
  void    SetFitPtsPart(const Int_t, const Int_t);		// Sets FitPoints for particles wrt Reaction plane  
  void    SetDedxPtsPart(const Int_t, const Int_t);		// Sets dE/dx for particles wrt Reaction plane  	
  void    SetFitOverMaxPtsPart(const Float_t, const Float_t);	// Sets FitPoints/MaxPoints for particles wrt Reaction plane 
  void    SetChiSqPart(const Float_t, const Float_t);		// Sets Chi^2 for particles wrt Reaction plane  	
  void    SetDcaGlobalPart(const Float_t, const Float_t);	// Sets d.c.a. for particles wrt Reaction plane  	
  void    SetConstrainablePart(Bool_t constr = kTRUE);		// Sets constrainability for particles wrt Reaction plane 	

 // Cuts set (V0 Analysis)
  void    SetV0Pid(const Char_t*) ;				// Sets PID for v0 wrt plane (...)
  void    SetV0Mass(const Float_t, const Float_t) ;		// Sets invariant mass cut for v0 wrt plane 
  void    SetV0Pt(const Float_t, const Float_t) ;	  	// Sets pT for v0 wrt plane 
  void    SetV0P(const Float_t, const Float_t) ; 	   	// Sets Momentum for v0 wrt plane 
  void    SetV0Eta(const Float_t, const Float_t) ;	  	// Sets Eta cut for v0 wrt plane 
  void    SetV0EtaAbs(const Float_t, const Float_t) ;		// Sets |Eta| cut (absolute value) for v0 wrt plane 
  void    SetV0Y(const Float_t, const Float_t) ; 		// Sets Rapidity for v0 wrt plane 
  void    SetV0ChiSqPart(const Float_t, const Float_t) ; 	// Sets Chi^2 for v0 wrt plane 
  void    SetV0DcaCross(const Float_t, const Float_t) ;    	// Sets distance to the main vertex for v0 wrt plane
  void    SetV0Lenght(const Float_t, const Float_t) ;	    	// Sets distance to the main vertex in sigma units for v0 wrt plane
  void    SetV0LenghtOverSigma(const Float_t, const Float_t) ;	// Sets closest approach (between the 2 daughter tracks) for v0 wrt plane  
  void	  SetV0SideBands() ;           				// Includes the v0 sideband analysis wrt plane		
  void	  SetV0SideBands(const Float_t) ;           		// Includes the v0 sideband analysis and a width		

 // For just constrainable track analysis (main loop, R.P. excluded)
  void    SetJustLoopConstrainable() ;				// Sets for the analysis loop just over constrainable track
  Bool_t  JustLoopConstrainable() const ;			// kFALSE (default) or kTRUE (set it above) 

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

 // Cuts for all tracks entering the main loop (R.P. loop is excluded)
  Bool_t  fJustLoopConstrainable ;			    	// constrainability for tracks entering the main loop 

  ClassDef(AliFlowSelection,1)               			// macro for rootcint
}; 

#endif
