//
//  AliEveESDTracks.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#ifndef __AliEveESDTracks__
#define __AliEveESDTracks__

#include <AliEveTrack.h>


#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>

#include <TEveManager.h>
#include <TEveTrackPropagator.h>


class AliEveESDTracks
{
public:
    AliEveESDTracks();
    ~AliEveESDTracks();
    
    TEveElementList* ByCategory();
    TEveElementList* ByType();
    
    void SetColors(Color_t colors[9]){
        for(int i=0;i<9;i++){fColors[i]=colors[i];}
    }
    void SetWidths(Width_t widths[9]){
        for(int i=0;i<9;i++){fWidths[i]=widths[i];}
    }
    void SetDashBad(bool dashBad){fDashBad=dashBad;}
    
private:
    bool fUseIPonFailedITSrefit;
    bool fTrueField;
    bool fRKstepper;
    AliESDtrackCuts* fAnalCuts;
    
    void SetupPropagator(TEveTrackPropagator* trkProp,Float_t magF, Float_t maxR);
    TString GetTitle(AliESDtrack* t);
    void AddParam(AliEveTrack* track, const AliExternalTrackParam* tp);
    AliEveTrack* MakeTPCtrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont);
    AliEveTrack* MakeITSstandaloneTrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont);

    AliEveTrack* MakeITStrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont);
    AliEveTrack* MakeTrack(AliESDtrack *at, TEveTrackList* cont);
    TEveTrackList* TPCtracks();
    TEveTrackList* ITStracks();
    TEveTrackList* ITSstandaloneTracks();
    TEveTrackList* Tracks();
    TEveTrackList* MItracks();
    TEveTrackList* TracksFromArray(TCollection* col, AliESDEvent* esd=0);
    void AliAnalCutsDemo();
    Float_t GetSigmaToVertex(AliESDtrack* esdTrack);

    TEveElementList* ByAnalCuts();
    Color_t fColors[15];
    Width_t fWidths[15];
    bool fDashBad;
    
    AliEveESDTracks(const AliEveESDTracks&);
    AliEveESDTracks& operator=(const AliEveESDTracks&);
};

#endif
