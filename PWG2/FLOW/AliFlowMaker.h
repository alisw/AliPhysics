//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: parser class from AliESD to AliFlowEvent . 
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWMAKER_H
#define ALIFLOWMAKER_H

class AliFlowEvent ;
class AliFlowTrack ;
class AliFlowV0 ;

class AliESD ;
class AliESDVertex ;
class AliESDtrack ;
class AliESDv0 ;

class AliFlowMaker  {

  public:
    AliFlowMaker();
    virtual ~AliFlowMaker();

  // FLOW SPECIFIC METHODS (to fill the flowEvents)
    AliFlowEvent*  FillFlowEvent(AliESD* fESD) ;	    // fills up evt quantities 
    AliFlowTrack*  FillFlowTrack(AliESDtrack* fTrack) ;     // fills up track quantities (global & constrained) ; p.id. & bayesian calculation ; 
    AliFlowV0*	   FillFlowV0(AliESDv0* fV0) ;  	    // fills up v0 quantities ; links v0 to tracks and back ;

  // USEFULL METHODS
    Double_t	   Norm(Double_t nu[3]) ;		    // norm of a non-vector 3 array	 
    Double_t	   Phi(Double_t nu[3]) ;		    // phi of a non-vector 3 array	 
    Double_t	   Pt(Double_t nu[3]) ; 		    // pt of a non-vector 3 array	
    Double_t	   Eta(Double_t nu[3]) ;		    // eta of a non-vector 3 array	 

  // Cut METHODS
    Bool_t         CheckTrack(AliESDtrack* fTrack) ; 	    // checks track (applies track cuts)
    Bool_t         CheckV0(AliESDv0* fV0) ; 	    	    // checks v0 (dummy)
    Bool_t         CheckEvent(AliESD* fESD) ; 	    	    // checks event (dummy)
    void 	   PrintCutList() ;			    // prints the list of cuts
    void 	   SetNHitsCut(Int_t nHits) 	            { fNHits = nHits ; }			// exclude tracks with less than .. TPC hits 
    void 	   SetECut(Float_t eLow, Float_t eUp)       { fElow = eLow ; fEup = eUp ; }		// exclude tracks below and above .. GeV 
    void 	   SetLabelCut(Int_t labLo, Int_t labHi)    { fLabel[0] = labLo ; fLabel[1] = labHi ; } // exclude tracks outside label interval

  // Get METHODS
    Int_t          GetNgoodTracks()         		    { return fGoodTracks ; }	  
    Int_t          GetNgoodV0s()            		    { return fGoodV0s ; }	  
    Int_t          GetNgoodTracksEta()      		    { return fGoodTracksEta ; } 
    Int_t          GetNposiTracks()         		    { return fPosiTracks ; }	
    Int_t          GetNnegaTracks() 	    		    { return fNegaTracks ; }	
    Int_t          GetNunconstrained()      		    { return fUnconstrained ; } 
    Int_t          GetBayesian(Int_t i = 2)     	    { return fBayesianAll[i] ; }
    Float_t        GetBayesianNorm(Int_t i = 2) 	    { return (Float_t)fBayesianAll[i] / (Float_t)fSumAll ; }

 protected:
 
  // enumerators 			    
    Int_t            fEventNumber ;	  		    //! progressive enumeration of ESD events
    Int_t            fTrackNumber ;	  		    //! progressive enumeration of ESD tracks
    Int_t            fV0Number ;	  		    //! progressive enumeration of ESD V0

    Int_t            fGoodTracks ;	  		    //! enumerator for good tracks
    Int_t            fGoodV0s ; 	  		    //! enumerator for good v0s
    Int_t            fGoodTracksEta ;	  		    //! enumerator for good tracks in the good eta range (-0.9..0.9)
    Int_t            fPosiTracks ;	  		    //! enumerator for positive tracks
    Int_t            fNegaTracks ;	  		    //! enumerator for negative tracks
    Int_t            fUnconstrained ;	  		    //! enumerator for tracks not constrainable
    Int_t            fBayesianAll[5] ;    		    //! final particles abundance -> AliFlowEvent (see Bayesian P.Id.)
    Int_t            fSumAll ;  	  		    //! total particles abundance (all kind)

    Int_t            fCutEvts ; 	  		    //! total enumerator for discarded events
    Int_t            fCutTrks ; 	  		    //! total enumerator for discarded tracks
    Int_t            fCutV0s ;  	  		    //! total enumerator for discarded V0s
 
  // Flags
    Bool_t           fNewAli ;	  	  		    //! enables the new ESD features (since AliRoot 12/2006) 
    Bool_t           fLoopTrks ;	  		    //! flag to loop over tracks 
    Bool_t           fLoopV0s ; 	  		    //! flag to loop over v0s 

    Int_t 	     fCounter ;        	  		    //! number of processed events

 private:

  // ESDs
    AliESD*          fESD;      	  		    //! "ESD" branch in fChain
    AliESDtrack*     fTrack;    	  		    //! "ESD track" 
    AliESDv0*        fV0;       	  		    //! "ESD v0" 
    AliESDVertex*    fVertex;       	  		    //! "ESD primary vertex"  

    Int_t            fRunID;		  		    //! last run ID
    Int_t            fNumberOfEvents ;    		    //! total number of ESD events in file
    Int_t            fNumberOfTracks ;    		    //! total number of tracks in the current event
    Int_t            fNumberOfV0s ;	  		    //! total number of v0s in the current event
    Float_t 	     fMagField ; 	  		    //! magnetic field from the ESD

  // Flow
    AliFlowEvent*    fFlowEvent ;         		    //! pointer to flow event
    AliFlowTrack*    fFlowTrack;          		    //! pointer to flow track
    AliFlowV0*       fFlowV0;          	  		    //! pointer to flow V0

  // Tracks cuts
    Int_t   fNHits;            	        		    // exclude tracks with less than .. TPC hits 
    Float_t fElow ;		        		    // exclude tracks below .. GeV (~total Momentum)
    Float_t fEup ;		        		    // exclude tracks above .. GeV (~total Momentum)
    Int_t   fLabel[2] ;	                		    // exclude tracks outside label interval
 
  ClassDef(AliFlowMaker,0);
};			

#endif
