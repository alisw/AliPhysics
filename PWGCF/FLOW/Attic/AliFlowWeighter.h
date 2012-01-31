//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowWeighter.h 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: generates phi-weights and counts particle abundances .
//  So, in fact, you should run this thing before the analysis,
//  if you want to use phi-weight.
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWWEIGHTER_H
#define ALIFLOWWEIGHTER_H

#include <TFile.h>
#include "AliFlowConstants.h"

class TH1F;
class TH1D;
class TOrdCollection;
class TVector2;

class AliFlowTrack;
class AliFlowV0;
class AliFlowEvent;
class AliFlowSelection;

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

  void 	   PrintBayesian(Int_t selN = 0) const ; 		// Prints normalized particle abundance (selN)


 protected:
 
 // Internal methods to fill the histogram
  void     TracksLoop(TObjArray* fFlowTracks) ;  		// Fills Phi and PId histograms
  Bool_t   Weightening() ;  					// Calculates weights and fills PhiWgt histograms


 private:

 // to make the code checker happy
  AliFlowWeighter(const AliFlowWeighter &flowWgt) ; 		// Copy Constructor (dummy)
  AliFlowWeighter &operator=(const AliFlowWeighter &flowAnal) ; // Assignment Operator

 // enumerators etc.			    
  Int_t            fEventNumber ;	  		    	//! progressive enumeration of AliFlowEvents
  Int_t            fTrackNumber ;	  		    	//! progressive enumeration of AliFlowTracks
  Int_t            fNumberOfV0s ;	  		        //! total number of V0s in the current event
  Int_t 	   fNumberOfTracks ;			    	//! total number of tracks in the current event

  Int_t            fPhiBins ;   			    	//! phi bins   
  Float_t 	   fPhiMin ;			    		//! i.e. 0
  Float_t  	   fPhiMax ; 			    		//! i.e. 2 Pi

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
  struct AliHistFullHars 
  {
    TH1D*       fHistPhiPlus;					//! histogram ...
    TH1D*       fHistPhiMinus;					//! histogram ...
    TH1D*       fHistPhiAll;					//! histogram ...
    TH1D*       fHistPhi;					//! histogram ...
    TH1D*       fHistPhiWgtPlus;				//! histogram ...
    TH1D*       fHistPhiWgtMinus;				//! histogram ...
    TH1D*       fHistPhiWgtAll; 				//! histogram ...
    TH1D*       fHistPhiWgt;					//! histogram ...
    TH1D*       fHistPhiFlatPlus;				//! histogram ...
    TH1D*       fHistPhiFlatMinus;				//! histogram ...
    TH1D*       fHistPhiFlatAll;				//! histogram ...
    TH1D*       fHistPhiFlat;					//! histogram ...
  };
  
  struct AliHistFulls 
  {
   TH1F*     	fHistBayPidMult ;				//! histogram ...
   struct AliHistFullHars fHistFullHar[AliFlowConstants::kHars];		//! structure wgt, evts, trks, v0s (as defined above)
  };
  struct AliHistFulls fHistFull[AliFlowConstants::kSels];                   //! structure array ...

  ClassDef(AliFlowWeighter,0)              			// macro for rootcint
};

#endif



