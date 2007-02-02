//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: generates phi-weights and counts particle abundances .
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWEIGHTER_H
#define ALIFLOWEIGHTER_H

#include <TVector2.h>
#include <TFile.h>
#include "AliFlowConstants.h"

class TH1F;
class TH1D;
class TOrdCollection;

class AliFlowTrack;
class AliFlowV0;
class AliFlowEvent;
class AliFlowSelection;
class Flow;

class AliFlowWeighter {


 public:

  AliFlowWeighter(const AliFlowSelection* flowSelect = 0); 	// Constructor with selection object (default selection if no one given)
  virtual  ~AliFlowWeighter(); 					// Default destructor (no actions)

 // Steps of the flow analysis
  Bool_t   Init() ;						// Books wgt histograms, opens output file
  Bool_t   Finish() ;						// Saves histograms, closes stuff

 // Analysis of 1 event (can be called from outside)
  Bool_t   WeightEvent(AliFlowEvent* fFlowEvent = 0) ; 		// Runs on the event

 // Output 
  void	   SetWgtFileName(TString name) 			{ fWgtFileName = name ; }  				   // Sets output file name
  void	   SetWgtFile(TFile* file) 				{ fWgtFile = file ; fWgtFileName = fWgtFile->GetName() ; }  // Sets output file
  TString  GetWgtFileName() const				{ return fWgtFileName ; }

  void 	   PrintBayesian(Int_t selN = 0) ; 			// Prints normalized particle abundance (selN)


 protected:
 
 // Internal methods to fill the histogram
  void     TracksLoop(TObjArray* fFlowTracks) ;  		// Fills Phi and PId histograms
  Bool_t   Weightening() ;  					// Calculates weights and fills PhiWgt histograms


 private:

 // enumerators etc.			    
  Int_t            fEventNumber ;	  		    	//! progressive enumeration of AliFlowEvents
  Int_t            fTrackNumber ;	  		    	//! progressive enumeration of AliFlowTracks
  Int_t            fNumberOfV0s ;	  		        //! total number of V0s in the current event
  Int_t 	   fNumberOfTracks ;			    	//! total number of tracks in the current event

  Int_t            fPhiBins ;     
  Float_t 	   fPhiMin ;
  Float_t  	   fPhiMax ; 

 // Internal pointers
  AliFlowEvent*     fFlowEvent ;      				//! pointer to AliFlowEvent
  AliFlowTrack*     fFlowTrack ;      				//! pointer to AliFlowTrack
  AliFlowSelection* fFlowSelect ;     				//  selection object
  TObjArray*        fFlowTracks ;     				//! pointer to the TrackCollection
  //Float_t           fVertex[3] ;				//! Event's Vertex position 

 // PhiWgt File
  TFile* 	    fWgtFile ;  			      	//! phi weight file 
  TString	    fWgtFileName ;				//! Wgt File Name (histograms for weight)

 // Histograms
  TOrdCollection*   fPhiWgtHistList ;     			//! Weights:  histogram list
  struct fHistFullHars 
  {
    TH1D*       fHistPhiPlus;
    TH1D*       fHistPhiMinus;
    TH1D*       fHistPhiAll;
    TH1D*       fHistPhi;
    TH1D*       fHistPhiWgtPlus;
    TH1D*       fHistPhiWgtMinus;
    TH1D*       fHistPhiWgtAll;
    TH1D*       fHistPhiWgt;
    TH1D*       fHistPhiFlatPlus;
    TH1D*       fHistPhiFlatMinus;
    TH1D*       fHistPhiFlatAll;
    TH1D*       fHistPhiFlat;
  };
  
  struct fHistFulls;	
  friend struct fHistFulls;
  struct fHistFulls 
  {
   TH1F*     	fHistBayPidMult;
   struct fHistFullHars fHistFullHar[Flow::nHars];		// wgt, evts, trks, v0s (as defined above)
  };
  struct fHistFulls fHistFull[Flow::nSels];                     //!

  ClassDef(AliFlowWeighter,0)              			// macro for rootcint
};

#endif



