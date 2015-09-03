//
// Created by mgrochow on 9/2/15.
//

#ifndef ALIROOT_ALICONVERTERPOLYLINESENGINE_H
#define ALIROOT_ALICONVERTERPOLYLINESENGINE_H

#include <TRint.h>

#include <AliMinimalisticEvent.h>
#include <AliMinimalisticTrack.h>
#include <AliESDEvent.h>
#include <TEveVector.h>


class AliConverterPolylinesEngine {
public:
    AliConverterPolylinesEngine();
    ~AliConverterPolylinesEngine();
    void AssertGeometry() const;
    void AssertMagField() const;
    void AddPolylinesToMuonTracks(
            Int_t trackNumber, AliMinimalisticTrack &minimalisticMuonTrack
    ) const;
    void AddPolylinesToMinimalisticTrack(
            Int_t trackID,
            AliMinimalisticTrack &minimalisticTrack
    ) const;
    void AddPolyLinesToKinkTrack(
            Int_t kinkID,
            AliMinimalisticTrack &mTrack,
            AliMinimalisticTrack &dTrack
    ) const;
    void AddPolyLinesToV0Track(
            Int_t v0ID,
            AliMinimalisticTrack &negativeTrack,
            AliMinimalisticTrack &positiveTrack
    ) const;
    void AddPolylinesToCascade(
            Int_t cascadeID,
            AliMinimalisticTrack &negativeTrack,
            AliMinimalisticTrack &positiveTrack,
            AliMinimalisticTrack &bachelorTrack
    ) const;
    void InsertPolyPoints(
            AliMinimalisticTrack &Track, std::vector<TEveVector4D> &Points
    ) const;
    void InitializeEngine(AliESDEvent *event);
private:
    TRint *fApp;
    AliESDEvent *fESDEvent;
};


#endif //ALIROOT_ALICONVERTERPOLYLINESENGINE_H
