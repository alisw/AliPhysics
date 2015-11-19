//
//  AliEveAODTracks.h
//
//  Created by Jeremi Niedziela on 17/11/15.
//
//

#ifndef __AliEveAODTracks__
#define __AliEveAODTracks__

#include <AliEveTrack.h>


#include <AliAODtrack.h>

#include <TEveManager.h>
#include <TEveTrackPropagator.h>


class AliEveAODTracks
{
public:
    AliEveAODTracks();
    ~AliEveAODTracks();
    
    TEveElementList* ByPID();
    
    void SetColorsByCategory(Color_t colors[9]){
        for(int i=0;i<9;i++){fColorsByCategory[i]=colors[i];}
    }
    void SetWidth(Width_t width){fWidth=width;}
    void SetDashNoRefit(bool dashNoRefit){fDashNoRefit=dashNoRefit;}
    void SetDrawNoRefit(bool drawNoRefit){fDrawNoRefit=drawNoRefit;}
    
private:
    bool fUseIPonFailedITSrefit;
    bool fTrueField;
    bool fRKstepper;
    
    void SetupPropagator(TEveTrackPropagator* trkProp,Float_t magF, Float_t maxR);
    TString GetTitle(AliAODTrack* t);
    void AddParam(AliEveTrack* track, const AliExternalTrackParam* tp);
    

    AliEveTrack* MakeTrack(AliAODTrack *at, TEveTrackList* cont);
    
    Color_t fColorsByCategory[9];
    Width_t fWidth;
    bool fDashNoRefit;
    bool fDrawNoRefit;
    
    AliEveAODTracks(const AliEveAODTracks&);
    AliEveAODTracks& operator=(const AliEveAODTracks&);
};

#endif
