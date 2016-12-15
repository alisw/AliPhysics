#ifndef ALIMUONCOMPACTCLUSTER_H 
#define ALIMUONCOMPACTCLUSTER_H 

#include <iostream>

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

    int mBendingManuIx;
    int mNonBendingManuIx;
};

#endif

