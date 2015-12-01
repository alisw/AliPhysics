/// \class AliExternalFormatConverter
/// Class for computing polylines. Used by ALiExternalFormatConverter during conversion from ESD to JSON/XML
/// Usage:
///     1) Create converter instance and initialize it by passing :
///         AliESDEvent* event; // This should be initialized
///         AliConverterPolylinesEngine fPolylineEngine;
///         fPolylineEngine.InitializeEngine(event);
///     2) Use one of the availabe methods e.g.:
///         Int_t v0Entry = 1 // ID of the V0
///         AliMinimalisticTrack negative;
///         AliMinimalisticTrack positive;
///         fPolylineEngine.AddPolyLinesToV0Track(v0Entry, negative, positive);
//
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

#ifndef ALIROOT_ALICONVERTERPOLYLINESENGINE_H
#define ALIROOT_ALICONVERTERPOLYLINESENGINE_H

#include <TRint.h>

#include <ConversionConstants.h>
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
    AliConverterPolylinesEngine(const AliConverterPolylinesEngine&) {};
    AliConverterPolylinesEngine& operator=(const AliConverterPolylinesEngine&) {};
    TRint *fApp;
    AliESDEvent *fESDEvent;
};


#endif //ALIROOT_ALICONVERTERPOLYLINESENGINE_H
