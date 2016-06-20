#include "AliMinimalisticEvent.h"

ClassImp(AliMinimalisticEvent);

/// Adds minimalistic track (it has to be already extracted) and saved it inside minimalistic event
void AliMinimalisticEvent::AddTrack(const AliMinimalisticTrack &track)
{
    fTracks.push_back(track);
}

/// Adds minimalistic cluster (it has to be already extracted) and saved it inside minimalistic event
void AliMinimalisticEvent::AddCluster(const AliMinimalisticCluster &cluster)
{
    fClusters.push_back(cluster);
}

/// Adds minimalistic calo cluster (it has to be already extracted) and saved it inside minimalistic event
void AliMinimalisticEvent::AddCaloCluster(const AliMinimalisticCaloCluster &cluster)
{
    fCaloClusters.push_back(cluster);
}


/// Ctor -- set the minimalistic event up
AliMinimalisticEvent::AliMinimalisticEvent(
        Double_t energy,
        Int_t multiplicity,
        TString collidingSystem,
        time_t timeStamp
) : TObject(),
    fEnergy(energy),
    fMultiplicity(multiplicity),
    fCollidingSystem(collidingSystem),
    fTimeStamp(timeStamp)
{

}
