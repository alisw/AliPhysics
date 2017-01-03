#ifndef ALIMUONCOMPACTEVENT_H
#define ALIMUONCOMPACTEVENT_H

#include <vector>
#include "AliMuonCompactTrack.h"

/**

  @ingroup pwg_muondep_compact

  @struct AliMuonCompactEvent

  @brief A compact event (list of tracks and initial basic kinematics)
*/

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
