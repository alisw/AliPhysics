// Authors : M.Gheata, A.Gheata 09/10/00
#ifndef MUON_RECEVENT
#define MUON_RECEVENT

////////////////////////////////////////////////////////////////////
//                                                                //
// AliMUONRecoEvent and AliMUONRecoTrack                          //
//                                                                //
// This header contains all classes needed for storing MUON       //
// reconstructed events in a tree                                 //
//                                                                //
////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TFile.h>
#include <TParticle.h>
#include <AliDetector.h>
#include "AliMUONHit.h"

class AliMUONRecoTrack;


/////////////////////////////////////////////////////////////////////
//                                                                 //
// AliMUONRecoEvent                                                //
//                                                                 //
// This class handles an array of reconstructed tracks.            //
// It provides :                                                   //
//	- filling the tracks array according the information stored//
//	in AliMUONEventReconstructor class ;                       //
//	- printing event and track informations : event number,    //
//	number of tracks, hits positions, reconstr. momentum.      //
//                                                                 //
/////////////////////////////////////////////////////////////////////


class AliMUONRecoEvent:public TObject {

private:
   Int_t             fNevr;               // event number
   Int_t             fNtracks;            // number of tracks
   TClonesArray      *fTracks;            // list of AliMUONRecoTracks
   
public:
                     AliMUONRecoEvent(Int_t eventNo = 0);
   virtual           ~AliMUONRecoEvent();
   AliMUONRecoTrack* AddEmptyTrack();
   void              Clear(Option_t *option = "");
   void              EventInfo();
   Int_t             GetNoEvent()  const {return fNevr;}
   Int_t             GetNoTracks() const {return fNtracks;}
   Bool_t            MakeDumpTracks(TClonesArray *tracksPtr);
   void              SetNoEvent(Int_t event) 	{fNevr = event;}
   void              SetNoTracks(Int_t ntracks) {fNtracks = ntracks;} 
   TClonesArray*     TracksPtr() {return fTracks;}

   ClassDef(AliMUONRecoEvent,1)	// Reconstructed event for MUON module
};

////////////////////////////////////////////////////////////////////
//                                                                //
// AliMUONRecoTrack                                               //
//                                                                //
// This class represents a reconstructed muon track.              //
//                                                                //
////////////////////////////////////////////////////////////////////

class AliMUONRecoTrack:public TObject {

private:
    Int_t       fSign;                  // charge sign
    Double_t 	fZvr;                   // z of track vertex point
    Double_t 	fChi2r;	                // chi squared for reco. track
    Double_t 	fPr[3];	                // reconstr. momentum (same as in vertex)
    Double_t 	fPosX[10];              // hit X position in all chambers
    Double_t 	fPosY[10];              // hit Y position in all chambers    
    Double_t 	fPosZ[10];              // hit Z position in all chambers

public:
   AliMUONRecoTrack() { }
    AliMUONRecoTrack(Bool_t active);
    virtual        ~AliMUONRecoTrack() { }	//desctructor
    const Double_t GetChi2r() const {return fChi2r;};
    const Double_t GetMomReconstr(Int_t axis) const {return fPr[axis];};
    const Int_t    GetSign() const {return fSign;};
    const Double_t GetPosX(Int_t chamber) const {return fPosX[chamber];};
    const Double_t GetPosY(Int_t chamber) const {return fPosY[chamber];};
    const Double_t GetPosZ(Int_t chamber) const {return fPosZ[chamber];};
    const Double_t GetVertexPos() { return fZvr;};
    const Double_t P() {return TMath::Sqrt(fPr[0]*fPr[0] + fPr[1]*fPr[1] + fPr[2]*fPr[2]);};
    const Double_t Phi();
    void           SetChi2r(Double_t chi) { fChi2r = chi;};
    void           SetHitPosition(Int_t chamber, Double_t x, Double_t y, Double_t z);
    void           SetMomReconstr(Double_t px, Double_t py, Double_t pz);
    void           SetSign(Int_t sign) {fSign = sign;};
    void           SetVertexPos(Double_t zvr) {fZvr = zvr;};
    const Double_t Theta();
    void           TrackInfo();
    ClassDef(AliMUONRecoTrack,1)	// A reconstructed muon track
};

#endif