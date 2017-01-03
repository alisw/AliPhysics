#ifndef ALIMUONCOMPACTTRACK_H
#define ALIMUONCOMPACTTRACK_H

#include <vector>
#include <iostream>
#include "AliMuonCompactCluster.h"

#include "Rtypes.h"

/**

  @ingroup pwg_muondep_compact

  @struct AliMuonCompactTrack

  @brief A very compact muon track
*/

struct AliMuonCompactTrack
{
    AliMuonCompactTrack(Double_t x=0.0, Double_t y=0.0, Double_t z=0.0) 
        : mPx(x), mPy(y), mPz(z), mClusters() {}
    Double_t Px() const { return mPx; }
    Double_t Py() const { return mPy; }
    Double_t Pz() const { return mPz; }

    Double_t mPx,mPy,mPz;
    std::vector<AliMuonCompactCluster> mClusters;

    friend std::ostream& operator<<(std::ostream& out,
            const AliMuonCompactTrack& ct);

};
#endif

