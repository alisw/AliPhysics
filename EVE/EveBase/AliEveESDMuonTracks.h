//
//  AliEveESDMuonTracks.h
//
//  Created by Jeremi Niedziela on 3/12/15
/// Main authors: P. Pillot, L. Aphecetche; Subatech
//

#ifndef __AliEveESDMuonTracks__
#define __AliEveESDMuonTracks__

#include <AliEveTrack.h>

#include <AliESDEvent.h>
#include <AliMUONESDInterface.h>

#include <TEveTrackPropagator.h>

class AliEveESDMuonTracks
{
public:
    AliEveESDMuonTracks(){}
    ~AliEveESDMuonTracks(){}
    
    void Draw(Bool_t showClusters=false, Bool_t showDigits=false);
    
private:
    void SetupTrackPropagator(TEveTrackPropagator* trkProp, Bool_t tracker, Bool_t trigger);
    void AddMuonTracks(AliESDEvent* esd, AliMUONESDInterface* data,TEveTrackList* match, TEveTrackList* nomatch, TEveTrackList* ghost);
    
    AliEveESDMuonTracks(const AliEveESDMuonTracks&);
    AliEveESDMuonTracks& operator=(const AliEveESDMuonTracks&);
};

#endif
