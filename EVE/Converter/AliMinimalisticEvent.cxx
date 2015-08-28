//
// Created by mgrochow on 7/23/15.
//

#include "AliMinimalisticEvent.h"

ClassImp(AliMinimalisticEvent);


void AliMinimalisticEvent::AddTrack(const AliMinimalisticTrack &track)
{
    fTracks.push_back(track);
}

void AliMinimalisticEvent::AddCluster(const AliMinimalisticCluster &cluster)
{
    fClusters.push_back(cluster);
}

AliMinimalisticEvent::AliMinimalisticEvent(
        Double_t energy,
        Int_t multiplicity,
        const std::string &collidingSystem,
        time_t timeStamp
) : TObject(),
    fEnergy(energy),
    fMultiplicity(multiplicity),
    fCollidingSystem(collidingSystem),
    fTimeStamp(timeStamp)
{

}
