/// \class AliExternalFormatConverter
/// Class for serializing ALICE data files (ESD, RecPoints, RawData) to
/// web-oriented lightweight formats: JSON, XML, CSV.
/// NOT all data from the file are serialized - only those which are indispensable for
/// 3D reconstruction purposes.
///
/// Usage:
///     1) Create converter instance:
///         const char* dirPath = "Path/To/Directory/With/ESD/Files";
///         AliExternalFormatConverter Converter(dirPath);
///     2) Convert it to external format:
///         Int_t entryNumber = #; //Number of entry in ESD file.
///         const char* jsonPath = "Path/To/Generated/file.json";
///         Converter.WriteJSONToFile(jsonPath, entryNumber); // JSON is an example
//
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

#ifndef ALIROOT_ALIEXTERNALFORMATCONVERTER_H
#define ALIROOT_ALIEXTERNALFORMATCONVERTER_H


#include <set>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TBufferJSON.h>
#include <TBufferXML.h>

#include <AliESDEvent.h>

#include <AliConverterPolylinesEngine.h>
#include <AliMinimalisticEvent.h>
#include <AliMinimalisticCluster.h>
#include <AliMinimalisticTrack.h>
#include <ConversionConstants.h>


class AliExternalFormatConverter {
public:
    virtual ~AliExternalFormatConverter();
    AliExternalFormatConverter();
    AliExternalFormatConverter(const TString dirPath);
    AliExternalFormatConverter(TFile* ES, TTree*);

    void LoadFiles(const TString dirPath);
    void LoadFiles(TFile *ESDFile, TFile *friendFile);

    void WriteJSONToFile(const Char_t *outputPath, Int_t entry);
    void WriteJSONToFile(const Char_t *outputPath, AliESDEvent *event);
    void WriteXMLToFile(const Char_t *outputPath, Int_t entry);
    void WriteXMLToFile(const Char_t *outputPath, AliESDEvent *event);

    TString GenerateJSON(Int_t eventEntry);
    TString GenerateJSON(AliESDEvent *event);
    TString GenerateXML(Int_t eventEntry);
    TString GenerateXML(AliESDEvent *event);

    void SerializeAllEvents(TString path);

private:
    AliExternalFormatConverter(const AliExternalFormatConverter&) {/*Converter cannot be copied*/};
    AliExternalFormatConverter& operator=(const AliExternalFormatConverter&) {/*Converter cannot be assigned*/};

    void LoadESDFile(const Char_t * ESDFilePath);
    void LoadESDFriends(const Char_t *friendFilePath);

    void WriteToFile(const Char_t *filePath, TString str) const;

    AliMinimalisticEvent GenerateMinimalisticEvent();
    void LoadEvent(Int_t entry);
    void LoadEvent(AliESDEvent *event);

    void PopulateEvent(AliMinimalisticEvent &event) const;
    void PopulateEventWithStandardTracks(AliMinimalisticEvent &event, std::set<Int_t> &usedTracks) const;
    void PopulateEventWithKinkTracks(AliMinimalisticEvent &event, std::set<Int_t> &kinkTracks) const;
    void PopulateEventWithV0Tracks(
            AliMinimalisticEvent &event, std::set<Int_t> &usedTracks, Int_t &specialID) const;
    void PopulateEventWithCascadeTracks(
            AliMinimalisticEvent &event, std::set<Int_t> &usedTracks, Int_t &specialID) const;
    void PopulateEventWithMuonTracks(AliMinimalisticEvent &event, Int_t &specialID) const;
    void PopulateEventWithCaloClusters(AliMinimalisticEvent &event) const;

    void AddContentToEvent(
            AliMinimalisticEvent &event,
            Int_t negative,
            Int_t parentID=AliMinimalisticTrack::fgkNoParent,
            Int_t childID=AliMinimalisticTrack::fgkNoChild
    ) const;
    AliMinimalisticTrack GenerateMinimalisticTrack(
            Int_t trackNumber, Int_t parentID, TrackType trackType
    ) const;
    AliMinimalisticTrack GenerateMinimalisticV0ParentTrack(
        AliESDv0 *V0,
        Int_t nChild,
        Int_t pChild,
        Int_t myID,
        Double_t startXYZ[3],
        Int_t type,
        Int_t parentID=AliMinimalisticTrack::fgkNoParent
    ) const;
    AliMinimalisticTrack GenerateMinimalisticCascadeParenTrack(
            AliESDcascade *cascade, Int_t v0ChildID, Int_t singleChildID, Int_t myID
    ) const;

    AliMinimalisticCluster GenerateMinimalisticCluster(Int_t trackID) const;
    void ExtractTrackPointArrays(AliMinimalisticCluster &cluster, Int_t trackNumber) const;

    void LoadEvent();
    void LoadESDFile(const TFile *ESDFile);
    void CheckEvent() const;

    TFile* fESDFile;
    TTree* fESDTree;
    AliESDEvent* fESDEvent;
    AliESDfriend* fESDFriend;
    AliConverterPolylinesEngine fPolylineEngine;
};


#endif //ALIROOT_ALIEXTERNALFORMATCONVERTER_H
