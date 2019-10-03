#ifndef ALIMUONCOMPACTCLUSTER_H 
#define ALIMUONCOMPACTCLUSTER_H 

#include <iostream>

/**

  @ingroup pwg_muondep_compact

  @struct AliMuonCompactCluster

  @brief A very minimal cluster

  The only information that we keep are the Manu locations (bending
  and/or non-bending).

*/

struct AliMuonCompactCluster
{
    AliMuonCompactCluster(int b=0, int nb=0)
        : mBendingManuIx(b), mNonBendingManuIx(nb) 
    {}

    int DetElemId() const;

    int BendingManuIndex() const { return mBendingManuIx; }
    int NonBendingManuIndex() const { return mNonBendingManuIx; }

    friend std::ostream& operator<<(std::ostream& out,const AliMuonCompactCluster& cl);

private:

    int mBendingManuIx; /// Absolute Bending Manu Index this cluster is in
    int mNonBendingManuIx; /// Absolute Non-bending Manu Index this cluster is in
};

#endif

