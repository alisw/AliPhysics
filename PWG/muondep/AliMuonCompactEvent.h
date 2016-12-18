#ifndef ALIMUONCOMPACTEVENT_H
#define ALIMUONCOMPACTEVENT_H

#include <vector>
#include "AliMuonCompactTrack.h"

struct AliMuonCompactEvent
{
    AliMuonCompactEvent() : mTracks() {}
    std::vector<AliMuonCompactTrack> mTracks;
    double mPt;
    double mY; 

    void Clear() { mTracks.clear(); }

    void SetInput(double pt, double y) { mPt = pt; mY = y; }

    friend std::ostream& operator<<(std::ostream& out,
            const AliMuonCompactEvent& event);
};

#endif 
