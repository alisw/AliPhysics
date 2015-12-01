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
    Int_t fEventID;
    Double_t fEnergy;
    Int_t fMultiplicity;
    std::string fCollidingSystem;
    std::time_t fTimeStamp;
    std::vector<AliMinimalisticTrack> fTracks;
    std::vector<AliMinimalisticCluster> fClusters;
    std::vector<AliMinimalisticCaloCluster> fCaloClusters;
};


#endif //ALIROOT_ALIMINIMALISTICEVENT_H
