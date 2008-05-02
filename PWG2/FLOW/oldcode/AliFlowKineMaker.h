//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowKineMaker.h 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: parser class from KineTree (pure MC simulation) to 
//  AliFlowEvent . It does not use the AliRunLoades, but simply 
//  gets the KineTree as an imput (a macro mast provide them).
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWKINEMAKER_H
#define ALIFLOWKINEMAKER_H

class AliFlowEvent ;
class AliFlowTrack ;
class AliFlowV0 ;
class TClonesArray ;

class TParticle ;
class TParticlePDG ;
class TTree ;
//class AliRun ;
//class AliRunLoader ;
//class AliStack ;

class AliFlowKineMaker { 

  public:
  
    AliFlowKineMaker();
    virtual ~AliFlowKineMaker();

  // FLOW SPECIFIC METHODS (to fill the flowEvents)
    AliFlowEvent*  FillFlowEvent(TTree* fKTree) ;	   // fills up evt quantities 
    AliFlowTrack*  FillFlowTrack(TParticle* fParticle) ;   // fills up track quantities ; 
    AliFlowV0*	   FillFlowV0(TParticle* fParticle) ;  	   // fills up v0 quantities ;

  // USEFULL METHODS
    Double_t	   Norm(Double_t nu[3]) ;		    // norm of a non-vector 3 array	 
    Double_t	   Phi(Double_t nu[3]) ;		    // phi of a non-vector 3 array	 
    Double_t	   Pt(Double_t nu[3]) ; 		    // pt of a non-vector 3 array	
    Double_t	   Eta(Double_t nu[3]) ;		    // eta of a non-vector 3 array	 

  // Cut METHODS
    Bool_t         CheckTrack(TParticle* fParticle) const ; // checks the particle (applies particle cuts, returns the particle charge)
    Bool_t         CheckEvent(TTree* fKTree) const ; 	    // checks the KineTree (dummy)
    void 	   PrintCutList() ;			    // prints the list of cuts

    void 	   SetAbsEtaCut(Float_t aEta) 	            { fAbsEta = aEta ; }			// exclude tracks with eta > aEta 
    void 	   SetECut(Float_t eLow, Float_t eUp)       { fElow = eLow ; fEup = eUp ; }		// exclude tracks below and above .. GeV 
    void 	   SetLabelCut(Int_t labLo, Int_t labHi)    { fLabel[0] = labLo ; fLabel[1] = labHi ; } // exclude tracks outside label interval
    void 	   SetPrimaryCut(Bool_t prim = kTRUE)       { fPrimary = prim ;	}                	// exclude secundaries 

  // Get METHODS
    Int_t          GetNgoodTracks() const         	    { return fGoodTracks ; }	  
    Int_t          GetNgoodV0s() const            	    { return fGoodV0s ; }	  
    Int_t          GetNgoodTracksEta() const      	    { return fGoodTracksEta ; } 
    Int_t          GetNposiTracks() const         	    { return fPosiTracks ; }	
    Int_t          GetNnegaTracks() const 	    	    { return fNegaTracks ; }	
    Int_t          GetNunconstrained()  const		    { return fUnconstrained ; } 
    Int_t          GetBayesian(Int_t i = 2) const     	    { return fBayesianAll[i] ; }
    Float_t        GetBayesianNorm(Int_t i = 2) const 	    { return (Float_t)fBayesianAll[i] / (Float_t)fSumAll ; }


 protected:
 
  // enumerators 			    
    Int_t            fEventNumber ;	  		    //! progressive enumeration of KineTree events
    Int_t            fPartNumber ;	  		    //! progressive enumeration of TParticle

    Int_t            fGoodTracks ;	  		    //! enumerator for good tracks
    Int_t            fGoodV0s ; 	  		    //! enumerator for good v0s
    Int_t            fGoodTracksEta ;	  		    //! enumerator for good tracks in the good eta range (-0.9..0.9)
    Int_t            fPosiTracks ;	  		    //! enumerator for positive tracks
    Int_t            fNegaTracks ;	  		    //! enumerator for negative tracks
    Int_t            fUnconstrained ;	  		    //! enumerator for tracks not constrainable
    Int_t            fBayesianAll[6] ;    		    //! final particles abundance -> AliFlowEvent (see Bayesian P.Id.)
    Int_t            fSumAll ;  	  		    //! total particles abundance (all kind)

    Int_t            fCutEvts ; 	  		    //! total enumerator for discarded events
    Int_t            fCutParts ; 	  		    //! total enumerator for discarded particles

  // Flags
    Bool_t           fNewAli ;	  	  		    //! enables the new features (since AliRoot 12/2006) 
    Bool_t           fLoopParts ;	  		    //! flag to loop over tracks 

    Int_t 	     fCounter ;        	  		    //! number of processed events


 private:

  // to make the code checker happy
   AliFlowKineMaker(const AliFlowKineMaker &flowMak) ; 		  // Copy Constructor (dummy)
   AliFlowKineMaker &operator=(const AliFlowKineMaker &flowMak) ; // Assignment Operator (dummy)

  // KineTree stuff
    TTree*     	     fKTree ;      	  		    //! KineTree
    TParticle* 	     fParticle ;      	  		    //! TParticle (momentum, decay, etc.)
    TParticlePDG*    fParticlePDG ;    	  		    //! TParticlePDG (type, charge, etc.)
    Int_t            fCharge ;		     		    //! charge of the TParticlePDG
    Float_t          fVertex[3] ;      	  		    //! primary vertex  
//    AliStack*        fStack ;  			    //! particle stack
//    AliRunLoader*    fRunLoader ;	     		    //! AliRunLoader
//    AliRun*          gAlice ;		     		    //! pointer to the AliRun (gAlice)

    Int_t            fRunID;		  		    //! last run ID
    Int_t            fNumberOfEvents ;    		    //! total number of KineTree events in file
    Int_t            fNumberOfParticles ;    		    //! total number of TParticles in the current event
    Float_t 	     fMagField ; 	  		    //! magnetic field from the ESD

  // Flow
    AliFlowEvent*    fFlowEvent ;         		    //! pointer to flow event
    AliFlowTrack*    fFlowTrack;          		    //! pointer to flow track
    AliFlowV0*       fFlowV0;          	  		    //! pointer to flow V0

  // Tracks cuts
    Float_t fAbsEta;            	        	    // exclude tracks with |eta| bigger than 
    Float_t fElow ;		        		    // exclude tracks below .. GeV (~total Momentum)
    Float_t fEup ;		        		    // exclude tracks above .. GeV (~total Momentum)
    Int_t   fLabel[2] ;	                		    // exclude tracks outside label interval
    Bool_t  fPrimary ;	                		    // exclude secundary tracks 
    
  ClassDef(AliFlowKineMaker,1);
};			

#endif

