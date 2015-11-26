//
//  AliEveAODTracks.h
//
//  Created by Jeremi Niedziela on 17/11/15.
//
//

#ifndef __AliEveAODTracks__
#define __AliEveAODTracks__

#include <AliAODTrack.h>
#include <AliEveTrack.h>

class AliEveAODTracks
{
public:
    AliEveAODTracks();
    ~AliEveAODTracks();
        
    TEveElementList* ByPID();
private:
    TString         GetTitle(AliAODTrack* t);
    void            AddParam(AliEveTrack* track, const AliExternalTrackParam* tp);
    AliEveTrack*    MakeTrack(AliAODTrack *at, TEveTrackList* cont);
    
    bool fDrawNoRefit;
    
    AliEveAODTracks(const AliEveAODTracks&);
    AliEveAODTracks& operator=(const AliEveAODTracks&);
};

#endif
