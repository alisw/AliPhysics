#include "AliMuonCompactEvent.h"
#include <iostream>

/// \ingroup compact 
std::ostream& operator<<(std::ostream& out, const AliMuonCompactEvent& event)
{
    out << "AliMuonCompactEvent Pt=" << event.mPt << " Y=" << event.mY << " nTracks=" << event.mTracks.size() << std::endl;

    for ( std::vector<AliMuonCompactTrack>::size_type i = 0;
            i < event.mTracks.size(); ++i )
    {
        const AliMuonCompactTrack& track = event.mTracks[i];
        out << "   " << track << std::endl;
    }
    return out;
}
