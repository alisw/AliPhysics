/// \class AliExternalFormatConverter
/// Class for serializing ALICE data files (ESD, RecPoints, RawData) to
/// web-oriented lightweight formats: JSON, XML, CSV.
/// NOT all data from the file are serialized - only those which are indispensable for
/// 3D reconstruction purposes.
///
/// Usage:
///     1) Create converter instance:
///         const char* filePath = "Path/To/Your/File";
///         const char* filenameFriends = "Path/To/Your/File/Friend";
///         AliExternalFormatConverter Converter(filePath, filenameFriends);
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
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveManager.h>
#include <TInterpreter.h>
#include <TRint.h>
#include <TEvePathMark.h>
#include <TEveVector.h>

#include <AliESDEvent.h>

#include <AliMinimalisticEvent.h>
#include <AliMinimalisticCluster.h>
#include <AliMinimalisticTrack.h>




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
    const static TString fgkDetector[23];

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
    void PopulateEventWithMuonTracks(AliMinimalisticEvent &event) const;

    void AddContentToEvent(
            AliMinimalisticEvent &event, Int_t negative, Int_t parentID=-1, Int_t childID=-1) const;

    AliMinimalisticTrack GenerateMinimalisticTrack(Int_t trackNumber, Int_t parentID, Int_t trackType) const;
    AliMinimalisticTrack GenerateMinimalisticV0ParentTrack(
        AliESDv0 *V0, Int_t nChild, Int_t pChild, Int_t myID, Double_t startXYZ[3], Int_t type, Int_t parentID=-1
    ) const;
    AliMinimalisticTrack GenerateMinimalisticCascadeParenTrack(
            AliESDcascade *cascade, Int_t v0ChildID, Int_t singleChildID, Int_t myID) const;

    AliMinimalisticCluster GenerateMinimalisticCluster(Int_t trackID) const;
    void ExtractTrackPointArrays(AliMinimalisticCluster &cluster, Int_t trackNumber) const;

    void LoadEvent();

    void LoadESDFile(const TFile *ESDFile);

    void AddPolylinesToMinimalisticTrack(Int_t trackID, AliMinimalisticTrack &minimalisticTrack) const;
    void AddPolyLinesToKinkTrack(Int_t kinkID, AliMinimalisticTrack &mTrack, AliMinimalisticTrack &dTrack) const;
    void AddPolyLinesToV0Track(
            Int_t v0ID, AliMinimalisticTrack &negativeTrack, AliMinimalisticTrack &positiveTrack) const;
   void AddPolylinesToCascade(
           Int_t cascadeID,
           AliMinimalisticTrack &negativeTrack,
           AliMinimalisticTrack &positiveTrack,
           AliMinimalisticTrack &bachelorTrack
   ) const;

    // TODO ADD PolyLines for muontracks
    void AddPolylinesToMuonTracks() const;
    void CheckEvent() const;

    void InsertPolyPoints(AliMinimalisticTrack &Track, std::vector<TEveVector4D> &Points) const;
    //void CalculateMagneticField();
    //Double_t fEventMagneticField;
    TFile* fESDFile;
    TTree* fESDTree;
    AliESDEvent* fESDEvent;
    AliESDfriend* fESDFriend;
    TRint *fApp;
};


#endif //ALIROOT_ALIEXTERNALFORMATCONVERTER_H
