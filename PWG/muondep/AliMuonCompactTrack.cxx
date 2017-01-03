#include "AliMuonCompactTrack.h"

/// \ingroup compact
std::ostream& operator<<(std::ostream& out, const AliMuonCompactTrack& track)
{
    out << "AliMuonCompactTrack ";
    for ( std::vector<AliMuonCompactCluster>::size_type i = 0;
            i < track.mClusters.size(); ++i )
    {
        const AliMuonCompactCluster& cl = track.mClusters[i];
        out << cl << " ";
    }
    return out;
}
