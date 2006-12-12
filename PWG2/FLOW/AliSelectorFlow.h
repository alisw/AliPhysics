/* $Id$ */
/* derived from AliSelector.h,v 1.10 2006/08/15 jgrosseo Exp $ */

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL
#ifndef ALISELECTORFLOW_H
#define ALISELECTORFLOW_H

// ANSI things
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

// ROOT things
#include <TROOT.h>
#include <TSelector.h>
#include <TVector3.h>
#include <TVector.h>
#include <TObject.h>
#include <TString.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

// // AliRoot things
// #include "AliESD.h"
// #include "AliESDtrack.h"
// #include "AliKalmanTrack.h"
// #include "AliITStrackV2.h"
// #include "AliESDVertex.h"
// #include "AliESDv0.h"

// FLOW things
#include "AliFlowEvent.h"
#include "AliFlowConstants.h"

class TFile;
class TChain;
class TTree;
class TParticle;

class AliESD;
class AliESDtrack;
class AliESDv0;

class AliSelectorFlow : public TSelector {

  public:
    AliSelectorFlow();
    virtual ~AliSelectorFlow();

    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree*);
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

  // output file name
    //void    SetFlowEventFileName(TString name) ;
    //TString GetFlowEventFileName() ;

  // FLOW SPECIFIC METHODS (to fill the flowEvents)
    AliFlowEvent*	FillFlowEvent(AliESD* fESD) ;		 // fills up evt quantities 
    AliFlowTrack*	FillFlowTrack(AliESDtrack* fTrack) ;	 // fills up track quantities (global & constrained) ; p.id. & bayesian calculation ; 
    AliFlowV0*		FillFlowV0(AliESDv0* fV0) ;  		 // fills up v0 quantities ; links v0 to tracks and back ;

  // USEFULL METHODS (to use 3-arrays)
    Double_t	   Norm(Double_t nu[3]) ;			 // norm of a non-vector 3 array      
    Double_t	   Phi(Double_t nu[3]) ;			 // phi of a non-vector 3 array       
    Double_t	   Pt(Double_t nu[3]) ; 			 // pt of a non-vector 3 array	     
    Double_t	   Eta(Double_t nu[3]) ;			 // eta of a non-vector 3 array       

 protected:
    TTree*  GetKinematics();
    void CheckOptions();

    TTree	*fTree;     	  //! pointer to the TTree containing the events
    AliESD*          fESD;      	  //! "ESD" branch in fChain
    AliESDtrack*     fTrack;    	  //! "ESD track" 
    AliESDv0*        fV0;       	  //! "ESD v0" 
    Int_t 	     fCountFiles ;        // number of processed file
    TString 	     fFlowEventFileName ; //! output file name 

 private:

    TFile*        pFlowfile ;
    AliFlowEvent* pFlowEvent ;          //! pointer to flow event
    AliFlowTrack* pFlowTrack;           //! pointer to flow track
    AliFlowV0*    pFlowV0;          	//! pointer to flow V0
 
  // enumerators 			     
    Int_t       fRunID;                 //! last run ID
    Int_t 	fEventNumber ;          //! progressive enumeration of ESD events
    Int_t       fTrackNumber ;  	//! progressive enumeration of ESD tracks
    Int_t       fV0Number ;		//! progressive enumeration of ESD V0
    Int_t       fNumberOfEvents ;       //! total number of ESD events in file
    Int_t       fNumberOfTracks ;	//! total number of tracks in the current event
    Int_t       fNumberOfV0s ;  	//! total number of v0s in the current event
  // other enumerators 			     
    Int_t  	fGoodTracks ;		//! enumerator for good tracks
    Int_t  	fGoodV0s ;		//! enumerator for good v0s
    Int_t  	fGoodTracksEta ;	//! enumerator for good tracks in the good eta range (-0.9..0.9)
    Int_t  	fPosiTracks ;		//! enumerator for positive tracks
    Int_t  	fNegaTracks ;		//! enumerator for negative tracks
    Int_t   	fUnconstrained ;	//! enumerator for tracks not constrainable
    Int_t       fBayesianAll[5] ;       //! final particles abundance -> AliFlowEvent (see Bayesian P.Id.)
    Int_t       fSumAll ; 		//! total particles abundance (all kind)

    Int_t       fCutEvts ;		//! total enumerator for discarded events
    Int_t       fCutTrks ;		//! total enumerator for discarded tracks
    Int_t       fCutV0s ;		//! total enumerator for discarded V0s
 
  // Flags
    Bool_t      fLoopV0 ;               //! flag to loop over v0s 
    Bool_t      fDoNothing ;            //! flag for a dummy execution 
    Bool_t      fOnFlyAnalysis ;        //! flag for on-fly analysis 

  // Cuts
    Float_t fEtrkLow ;
    Float_t fEtrkHig ;
    Int_t   fHitsTrk ;

  // ...
    Float_t 	fMagField ; 		//! magnetic field from the ESD
 
 
    void DeleteKinematicsFile();

    TFile*        fKineFile;            //! pointer to Kinematics.root if the file was opened

    AliSelectorFlow(const AliSelectorFlow&);
    AliSelectorFlow& operator=(const AliSelectorFlow&);

  ClassDef(AliSelectorFlow,0);
};

#endif
