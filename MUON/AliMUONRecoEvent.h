#ifndef ALIMUONRECOEVENT_H
#define ALIMUONRECOEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/


// Authors : M.Gheata, A.Gheata 09/10/00

#include <TObject.h>
#include <TFile.h>
#include <TParticle.h>
#include <AliDetector.h>
#include "AliMUONHit.h"
class AliMUONEventReconstructor;

class AliMUONRecoTrack;

/////////////////////////////////////////////////////////////////////
//                                                                 //
// AliMUONRecoEvent                                                //
//                                                                 //
// This class handles an array of reconstructed tracks.            //
// It provides :                                                   //
//	- filling the tracks array according to the information    //
//        stored in AliMUONEventReconstructor class ;              //
//	- printing event and track informations : event number,    //
//	  number of tracks, hits positions, reconstr. momentum.    //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONRecoEvent : public TObject 
{
  public:
    AliMUONRecoEvent(Int_t eventNo = 0);
    virtual ~AliMUONRecoEvent();

    AliMUONRecoTrack* AddEmptyTrack();
    void              Clear(Option_t *option = "");
    void              EventInfo();
    Int_t             GetNoEvent()  const {return fNevr;}
    Int_t             GetNoTracks() const {return fNtracks;}
    Bool_t            MakeDumpTracks(Int_t muons, TClonesArray *tracksPtr, AliMUONEventReconstructor *MuonReco);
    void              SetNoEvent(Int_t event) 	{fNevr = event;}
    void              SetNoTracks(Int_t ntracks) {fNtracks = ntracks;} 

    void              SetNoMuons(Int_t muons) {fMuons = muons;} 

    TClonesArray*     TracksPtr() {return fTracks;}
    
 protected:    
    AliMUONRecoEvent(const AliMUONRecoEvent& rhs);
    AliMUONRecoEvent& operator=(const AliMUONRecoEvent& rhs);

 private:
   Int_t             fNevr;          // event number
   Int_t             fNtracks;       // number of tracks
   Int_t             fMuons;         // number of muons within acceptance
   TClonesArray      *fTracks;      //-> list of AliMUONRecoTracks
   
   ClassDef(AliMUONRecoEvent,1)	// Reconstructed event for MUON module
};

#endif
