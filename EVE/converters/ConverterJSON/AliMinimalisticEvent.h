/// \class AliMinimalisticEvent
/// Class - container of indispensable data members for 3D reconstruction purposes.
/// It is used in the proces of conversion from ALICE data formats (e.g. ESD)
/// to external ones. SEE AliExternalFormatConverter
///
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

#ifndef ALIROOT_ALIMINIMALISTICEVENT_H
#define ALIROOT_ALIMINIMALISTICEVENT_H

#include <ctime>

#include <TObject.h>

#include <AliMinimalisticCluster.h>
#include <AliMinimalisticCaloCluster.h>
#include <AliMinimalisticTrack.h>


class AliMinimalisticEvent: public TObject  {
ClassDef(AliMinimalisticEvent, 1);
public:
    void AddTrack(const AliMinimalisticTrack& track);
    void AddCluster(const AliMinimalisticCluster& cluster);
    void AddCaloCluster(const AliMinimalisticCaloCluster& cluster);
    AliMinimalisticEvent(Double_t energy, Int_t multiplicity, TString collidingSystem, time_t timeStamp);
private:
    Int_t fEventID; /// ID of the event extracted from the run
    Double_t fEnergy; /// energy of the collision
    Int_t fMultiplicity; /// number of particles reconstructed
    std::string fCollidingSystem; /// information about type of event (e.g. proton-proton)
    std::time_t fTimeStamp; /// Time of the collision
    std::vector<AliMinimalisticTrack> fTracks; /// an array of minimalistic tracks
    std::vector<AliMinimalisticCluster> fClusters; /// an array of minimalistic clusters
    std::vector<AliMinimalisticCaloCluster> fCaloClusters; /// an array of minimalistic calo clusters
};


#endif //ALIROOT_ALIMINIMALISTICEVENT_H
