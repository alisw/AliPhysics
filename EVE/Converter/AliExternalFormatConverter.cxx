//
// Created by Maciej Grochowicz on 7/21/15.
//

#include "AliExternalFormatConverter.h"

#include <fstream>

#include <AliEveTrack.h>
#include <AliEveMUONTrack.h>

#include <AliESDcascade.h>
#include <AliESDkink.h>
#include <AliESDMuonTrack.h>
#include <AliESDtrack.h>
#include <AliGeomManager.h>
#include <AliPDG.h>
#include <AliEveKink.h>
#include <AliEveV0.h>
#include <AliEveCascade.h>

#include <AliMUONGeometryTransformer.h>


AliExternalFormatConverter::AliExternalFormatConverter()
: fESDFile(nullptr), fESDFriend(nullptr), fESDEvent(nullptr), fESDTree(nullptr), fPolylineEngine()
{ }


AliExternalFormatConverter::AliExternalFormatConverter(const TString dirPath)
: fESDFile(nullptr), fESDFriend(nullptr), fESDEvent(nullptr), fESDTree(nullptr), fPolylineEngine()
{

    LoadFiles(dirPath);
}

AliExternalFormatConverter::AliExternalFormatConverter(TFile *ESDFile, TTree *ESDTree)
{
    throw "Not implemented.";
    //    if (fESDFile)
    //        delete fESDFile;
    //    fESDFile = ESDFile;
}

AliExternalFormatConverter::~AliExternalFormatConverter()
{
    if (fESDFile)
        delete fESDFile;
    if (fESDEvent)
        delete fESDEvent;
    //if (fApp)
    // delete fApp;
}

/// Loads files from given paths. friendPath is not obligatory
void AliExternalFormatConverter::LoadFiles(const TString dirPath)
{
    TString filePath = dirPath + "/AliESDs.root";
    TString friendPath = dirPath + "/AliESDfriends.root";
    LoadESDFile(filePath);
    LoadESDFriends(friendPath);
    //    if (!friendPath){
    //        LoadESDFriends(friendPath);
    //        return; // exit after loading Friend
    //    }
    //    if (fESDFriend)
    //        delete fESDFriend; // If friendPath is not given, but fESDFriend has already been allocated.
    
}
void AliExternalFormatConverter::LoadFiles(TFile *ESDFile, TFile *friendFile)
{
    LoadESDFile(ESDFile);
    fESDTree->AddFriend("esdFriendTree", friendFile);
    fESDTree->SetBranchStatus("ESDfriend");
    fESDFriend = dynamic_cast<AliESDfriend *>(fESDEvent->FindListObject("AliESDfriend"));
}

void AliExternalFormatConverter::LoadESDFile(const TFile *ESDFile)
{
    if (fESDFile)
        delete fESDFile;
    fESDFile = ESDFile->CurrentFile();
    fESDTree = dynamic_cast<TTree *>(fESDFile->Get("esdTree"));
    LoadEvent();
}

/// Loads ESD file from given path.
void AliExternalFormatConverter::LoadESDFile(const Char_t *filePath)
{
    if (fESDFile)
        delete fESDFile;
    fESDFile = new TFile(filePath);
    fESDTree = dynamic_cast<TTree *>(fESDFile->Get("esdTree"));
    LoadEvent();
}

void AliExternalFormatConverter::LoadEvent()
{
    if (fESDEvent)
        delete fESDEvent;
    fESDEvent = new AliESDEvent;
    fESDEvent->ReadFromTree(fESDTree);
    fPolylineEngine.InitializeEngine(fESDEvent);
}

///Loads ESDFriend file from given path
void AliExternalFormatConverter::LoadESDFriends(const Char_t *friendPath)
{
    fESDTree->AddFriend("esdFriendTree", friendPath);
    fESDTree->SetBranchStatus("ESDfriend");
    fESDFriend = dynamic_cast<AliESDfriend *>(fESDEvent->FindListObject("AliESDfriend"));
    
    if (fESDFriend) {
        fESDTree->SetBranchAddress("ESDfriend.", &fESDFriend);
    } else {
        std::cerr <<
        "An error ocurred. AliESDfriend has not been loaded. Please check the path" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void AliExternalFormatConverter::WriteJSONToFile(const Char_t *path, Int_t entry)
{
    TString json = GenerateJSON(entry);
    WriteToFile(path, json);
}

void AliExternalFormatConverter::WriteJSONToFile(const Char_t *path, AliESDEvent *event)
{
    TString json = GenerateJSON(event);
    WriteToFile(path, json);
}

void AliExternalFormatConverter::WriteXMLToFile(const Char_t *path, Int_t entry)
{
    TString xml = GenerateXML(entry);
    WriteToFile(path, xml);
}

void AliExternalFormatConverter::WriteXMLToFile(const Char_t *path, AliESDEvent *event)
{
    TString xml = GenerateXML(event);
    WriteToFile(path, xml);
}
TString AliExternalFormatConverter::GenerateJSON(Int_t entry)
{
    LoadEvent(entry);
    AliMinimalisticEvent event = GenerateMinimalisticEvent();
    return TBufferJSON::ConvertToJSON(&event);
}

TString AliExternalFormatConverter::GenerateXML(Int_t entry)
{
    LoadEvent(entry);
    AliMinimalisticEvent event = GenerateMinimalisticEvent();
    return TBufferXML::ConvertToXML(&event);
}
TString AliExternalFormatConverter::GenerateJSON(AliESDEvent *event)
{
    LoadEvent(event);
    AliMinimalisticEvent minimalisticEvent = GenerateMinimalisticEvent();
    return TBufferJSON::ConvertToJSON(&minimalisticEvent);
}

TString AliExternalFormatConverter::GenerateXML(AliESDEvent *event)
{
    LoadEvent(event);
    AliMinimalisticEvent minimalisticEvent = GenerateMinimalisticEvent();
    return TBufferXML::ConvertToXML(&minimalisticEvent);
}

///This method 'guts' given TObject (from loaded file) and returns compact
/// AliMinimalisticEvent that contains only necessary fields
AliMinimalisticEvent AliExternalFormatConverter::GenerateMinimalisticEvent()
{
    const char *beamType = fESDEvent->GetBeamType();
    time_t time_stamp = fESDEvent->GetTimeStamp();
    Float_t energy = fESDEvent->GetBeamEnergy();
    Int_t multiplicity = fESDEvent->GetMultiplicity()->GetNumberOfTracklets();
    AliMinimalisticEvent event(energy, multiplicity, beamType, time_stamp);
    PopulateEvent(event);
    return event;
}

void AliExternalFormatConverter::LoadEvent(Int_t entry)
{
    fESDTree->GetEntry(entry);
    CheckEvent();
}

void AliExternalFormatConverter::LoadEvent(AliESDEvent *event)
{
    fESDEvent = event;
    CheckEvent();
}

void AliExternalFormatConverter::CheckEvent() const
{
    if (!fESDEvent){
        std::cerr << "Event is corrupted" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void AliExternalFormatConverter::PopulateEvent(AliMinimalisticEvent &event) const
{
    std::set<Int_t> usedTracks;
    Int_t specialID = fESDEvent->GetNumberOfTracks();
    
    PopulateEventWithCascadeTracks(event, usedTracks, specialID);
    PopulateEventWithV0Tracks(event, usedTracks, specialID);
    PopulateEventWithKinkTracks(event, usedTracks);
    PopulateEventWithStandardTracks(event, usedTracks);
    PopulateEventWithMuonTracks(event, specialID);
}

void AliExternalFormatConverter::PopulateEventWithStandardTracks(
        AliMinimalisticEvent &event, std::set<Int_t> &usedTracks
        ) const
{
    for (Int_t standardEntry = 0; standardEntry < fESDEvent->GetNumberOfTracks(); standardEntry++){
        if (usedTracks.find(standardEntry) != usedTracks.end())
            continue;
        AddContentToEvent(event, standardEntry);
        usedTracks.insert(standardEntry);
    }
}

void AliExternalFormatConverter::PopulateEventWithV0Tracks(
        AliMinimalisticEvent &event, std::set<Int_t> &usedTracks, Int_t &specialID) const
{
    AliESDv0 *v0;
    for (Int_t v0Entry = 0; v0Entry < fESDEvent->GetNumberOfV0s(); v0Entry++){
        v0 = fESDEvent->GetV0(v0Entry);
        if (!v0)
            continue;
        Int_t positiveID = v0->GetPindex();
        Int_t negativeID = v0->GetNindex();
        if (usedTracks.find(positiveID) != usedTracks.end() || usedTracks.find(negativeID) != usedTracks.end())
            continue;
        Int_t v0ParentID = specialID++;
        AliMinimalisticTrack negative =
        GenerateMinimalisticTrack(negativeID, v0ParentID, kV0NegativeDaughter);
        AliMinimalisticTrack positive =
        GenerateMinimalisticTrack(positiveID, v0ParentID, kV0PositiveDaughter);
        if (fESDFriend){
            AliMinimalisticCluster motherCluster = GenerateMinimalisticCluster(negativeID);
            AliMinimalisticCluster daughterCluster = GenerateMinimalisticCluster(positiveID);
            event.AddCluster(motherCluster);
            event.AddCluster(daughterCluster);
        }
        fPolylineEngine.AddPolyLinesToV0Track(v0Entry, negative, positive);
        Double_t startCOORDS[] = {.0, .0, .0};
        AliMinimalisticTrack V0Parenttrack = GenerateMinimalisticV0ParentTrack(
                v0, negativeID, positiveID, v0ParentID, startCOORDS, kV0Mother
        );
        event.AddTrack(V0Parenttrack);
        event.AddTrack(positive);
        event.AddTrack(negative);
        usedTracks.insert(negativeID);
        usedTracks.insert(positiveID);
    }
}

void AliExternalFormatConverter::PopulateEventWithCascadeTracks(
        AliMinimalisticEvent &event, std::set<Int_t> &usedTracks, Int_t &specialID) const
{
    AliESDcascade *cascade;
    for (int cascadeEntry = 0; cascadeEntry < fESDEvent->GetNumberOfCascades(); cascadeEntry++){
        cascade = fESDEvent->GetCascade(cascadeEntry);
        Int_t positiveID = cascade->GetPindex();
        Int_t negativeID = cascade->GetNindex();
        Int_t bachelorID = cascade->GetIndex();
        std::set<Int_t>::iterator tracksEnd = usedTracks.end();
        if (usedTracks.find(positiveID) != tracksEnd
            || usedTracks.find(negativeID) != tracksEnd
            || usedTracks.find(bachelorID) != tracksEnd
            ) {
            continue;
        }
        
        Int_t v0ParentID = specialID++;
        Int_t cascadeParentID = specialID++;
        AliMinimalisticTrack negative =
        GenerateMinimalisticTrack(negativeID, v0ParentID, kCascadeNegativeDaughter);
        AliMinimalisticTrack positive =
        GenerateMinimalisticTrack(positiveID, v0ParentID, kCascadePositiveDaughter);
        AliMinimalisticTrack bachelor =
        GenerateMinimalisticTrack(bachelorID, cascadeParentID, kCascadePrimaryDaughter);
        fPolylineEngine.AddPolylinesToCascade(cascadeEntry, negative, positive, bachelor);
        event.AddTrack(negative);
        event.AddTrack(positive);
        event.AddTrack(bachelor);
        if (fESDFriend){
            AliMinimalisticCluster motherCluster = GenerateMinimalisticCluster(negativeID);
            AliMinimalisticCluster daughterCluster = GenerateMinimalisticCluster(positiveID);
            AliMinimalisticCluster bachelorCluster = GenerateMinimalisticCluster(bachelorID);
            event.AddCluster(motherCluster);
            event.AddCluster(daughterCluster);
            event.AddCluster(bachelorCluster);
        }
        Double_t v0StartCoords[3]; cascade->XvYvZv(v0StartCoords);
        AliMinimalisticTrack V0Parenttrack =
                GenerateMinimalisticV0ParentTrack(
                        (AliESDv0*)cascade,
                        negativeID,
                        positiveID,
                        v0ParentID,
                        v0StartCoords,
                        kCascadeSecondaryMother,
                        cascadeParentID
                );
        event.AddTrack(V0Parenttrack);
        
        AliMinimalisticTrack cascadeParentTrack = GenerateMinimalisticCascadeParenTrack(
                cascade, v0ParentID, bachelorID, cascadeParentID
        );
        
        event.AddTrack(cascadeParentTrack);
        
        usedTracks.insert(negativeID);
        usedTracks.insert(positiveID);
        usedTracks.insert(bachelorID);
    }
}

void AliExternalFormatConverter::PopulateEventWithKinkTracks(
        AliMinimalisticEvent &event, std::set<Int_t> &usedTracks) const
{
    AliESDkink *kink;
    for (Int_t kinkEntry = 0; kinkEntry < fESDEvent->GetNumberOfKinks(); kinkEntry++){
        kink = fESDEvent->GetKink(kinkEntry);
        Int_t motherID = kink->GetIndex(0);
        Int_t daughterID = kink->GetIndex(1);
        if (usedTracks.find(motherID) != usedTracks.end() || usedTracks.find(daughterID) != usedTracks.end())
            continue;
        
        AliMinimalisticTrack daughter = GenerateMinimalisticTrack(
                daughterID, motherID, kKinkDaughter
        );
        AliMinimalisticTrack mother = GenerateMinimalisticTrack(
                motherID, AliMinimalisticTrack::fgkNoParent, kKinkMother
        );
        if (fESDFriend){
            AliMinimalisticCluster motherCluster = GenerateMinimalisticCluster(motherID);
            AliMinimalisticCluster daughterCluster = GenerateMinimalisticCluster(daughterID);
            event.AddCluster(motherCluster);
            event.AddCluster(daughterCluster);
        }
        fPolylineEngine.AddPolyLinesToKinkTrack(kinkEntry, mother, daughter);
        event.AddTrack(daughter);
        event.AddTrack(mother);
        usedTracks.insert(motherID);
        usedTracks.insert(daughterID);
    }
}

void AliExternalFormatConverter::AddContentToEvent(
        AliMinimalisticEvent &event, Int_t trackID, Int_t parentID, Int_t childID
) const
{
    AliMinimalisticTrack track = GenerateMinimalisticTrack(trackID, parentID, kStandard);
    if (childID != AliMinimalisticTrack::fgkNoChild)
        track.AddChild(childID);
    if (fESDFriend){
        AliMinimalisticCluster cluster = GenerateMinimalisticCluster(trackID);
        event.AddCluster(cluster);
    }
    fPolylineEngine.AddPolylinesToMinimalisticTrack(trackID, track);
    event.AddTrack(track);
}


void AliExternalFormatConverter::WriteToFile(const char *path, TString fileString) const
{
    ofstream outfile;
    outfile.open(path, std::ios::binary | std::ios::out);
    if (!outfile){
        std::cerr<<"\n\nCouldn't create output file!\n\n"<<std::endl;
        return;
    }
    outfile << fileString << std::endl;
    outfile.close();
}

AliMinimalisticCluster AliExternalFormatConverter::GenerateMinimalisticCluster(Int_t trackID) const
{
    AliMinimalisticCluster cluster(trackID);
    ExtractTrackPointArrays(cluster, trackID);
    return cluster;
}

AliMinimalisticTrack AliExternalFormatConverter::GenerateMinimalisticTrack(
        Int_t trackNumber, Int_t parentID, TrackType trackType) const
{
    AliESDtrack *track = fESDEvent->GetTrack(trackNumber);
    
    Double_t mass = track->GetMass();
    Int_t PID = AliPID::ParticleCode(track->GetPID());
    Double_t energy = track->E();
    
    Int_t charge = track->Charge();
    Double_t signedPT = track->GetSignedPt();
    Double_t theta = track->Theta();
    Double_t phi = track->Phi();
    Double_t b = fESDEvent->GetMagneticField() / 10.0;
    Double_t helixCurvature = track->GetC(b);
    Double_t startXYZ[3];
    track->GetXYZ(startXYZ);
    Double_t endXYZ[] = {0, 0, 0};
    Double_t pxpypz[3];
    track->GetPxPyPz(pxpypz);
    
    AliMinimalisticTrack minimalisticTrack = AliMinimalisticTrack(
            charge, energy, trackNumber, PID, mass, signedPT,
            startXYZ, endXYZ, pxpypz, parentID, phi, theta,
            helixCurvature, trackType
    );
    
    return minimalisticTrack;
}

AliMinimalisticTrack AliExternalFormatConverter::GenerateMinimalisticCascadeParenTrack(
     AliESDcascade *cascade, Int_t v0ChildID, Int_t singleChildID, Int_t myID
) const
{
    Int_t charge = cascade->Charge();
    Double_t v0q = ((AliESDv0*)cascade)->ChangeMassHypothesis();
    cascade->ChangeMassHypothesis(v0q); //! This has to be called before calling V0->E(). Check AliESDV0 docs.
    Double_t energy = cascade->E();
    Int_t PID = cascade->GetPdgCodeXi();
    Double_t mass = cascade->GetEffMassXi();
    Double_t signedPt = cascade->Pt();
    Double_t startXYZ[3] = {.0, .0, .0};
    Double_t endXYZ[3]; cascade->XvYvZv(endXYZ);
    Double_t PxPyPz[3]; cascade->PxPyPz(PxPyPz);
    Double_t phi = cascade->Phi();
    Double_t theta = cascade->Theta();
    Double_t helixCurvature = 0.0;
    Int_t parentID = AliMinimalisticTrack::fgkNoParent;

    AliMinimalisticTrack cascdeParent(
        charge, energy, myID, PID, mass, signedPt, startXYZ, endXYZ, PxPyPz, parentID,
        phi, theta, helixCurvature, kCascadePrimaryMother
    );
    cascdeParent.AddChild(v0ChildID);
    cascdeParent.AddChild(singleChildID);
    cascdeParent.AddPolyPoint(startXYZ);
    cascdeParent.AddPolyPoint(endXYZ);
    return cascdeParent;
}

void AliExternalFormatConverter::ExtractTrackPointArrays(
        AliMinimalisticCluster &cluster, Int_t trackNumber) const
{
    AliESDfriendTrack *track = fESDFriend->GetTrack(trackNumber);
    if (!track){
        std::cerr << "Corrupted friend track: " << trackNumber << std::endl;
        return;
    }
    const AliTrackPointArray *array = track->GetTrackPointArray();
    int nPoints = array->GetNPoints();
    cluster.InsertXArray(array->GetX(), nPoints);
    cluster.InsertYArray(array->GetY(), nPoints);
    cluster.InsertZArray(array->GetZ(), nPoints);
    const Char_t *description[nPoints];
    const UShort_t *volumeID = array->GetVolumeID();
    for (int i = 0; i < nPoints; i++){
        UShort_t detectorID = AliGeomManager::VolUIDToLayer(volumeID[i]);
        description[i] = fgkDetector[detectorID];
    }
    cluster.InsertValueDescription(description, nPoints);
    cluster.SetSource("ESD");
}

AliMinimalisticTrack AliExternalFormatConverter::GenerateMinimalisticV0ParentTrack(
        AliESDv0 *V0,
        Int_t negativeChild,
        Int_t positiveChild,
        Int_t myID,
        Double_t startXYZ[3],
        Int_t type,
        Int_t parentID
) const
{
    Int_t charge = V0->Charge();
    V0->ChangeMassHypothesis(); //! This has to be called before calling V0->E(). Check AliESDV0 docs.
    Double_t energy = V0->E();
    Int_t PID = V0->PdgCode();
    Double_t mass = V0->GetEffMass();
    Double_t signedPt = V0->Pt(); //! The particle is neutral therefore Pt is always positive.
    Double_t endXYZ[3]; V0->XvYvZv(endXYZ); //! V0 position is the end of parent track.
    Double_t PxPyPz[3]; V0->PxPyPz(PxPyPz);
    Double_t phi = V0->Phi();
    Double_t theta = V0->Theta();
    Double_t helixCurvature = 0.0;
    AliMinimalisticTrack parent(
            charge, energy, myID, PID, mass, signedPt, startXYZ, endXYZ, PxPyPz, parentID, phi,
            theta, helixCurvature, type);
    parent.AddChild(negativeChild);
    parent.AddChild(positiveChild);
    parent.AddPolyPoint(startXYZ);
    parent.AddPolyPoint(endXYZ);
    return parent;
}

void AliExternalFormatConverter::PopulateEventWithMuonTracks(
        AliMinimalisticEvent &event, Int_t &specialID
) const
{
    AliESDMuonTrack *mTrack;
    for (Int_t muonTrack = 0; muonTrack < fESDEvent->GetNumberOfMuonTracks(); muonTrack++){
        mTrack = fESDEvent->GetMuonTrack(muonTrack);
        if (mTrack->GetNHit() == 0)
            continue;
        Int_t charge = mTrack->Charge();
        Double_t energy = mTrack->E();
        Int_t parentID = AliMinimalisticTrack::fgkNoParent;
        Double_t signedPT = mTrack->Pt()*charge;
        Double_t mass = mTrack->M();
        Double_t PxPyPz[3]; mTrack->PxPyPz(PxPyPz);
        Double_t startXYZ[3];
        startXYZ[0] = mTrack->GetNonBendingCoorAtDCA();
        startXYZ[1] = mTrack->GetBendingCoorAtDCA();
        startXYZ[2] = mTrack->GetZ();
        Double_t endXYZ[] = {.0, .0, .0};
        Double_t theta = mTrack->Theta();
        Double_t phi = mTrack->Phi();
        AliMinimalisticTrack minimalisticMuonTrack(
                charge, energy, specialID++, AliPID::kMuon, mass, signedPT, startXYZ, endXYZ, PxPyPz,
                AliMinimalisticTrack::fgkNoParent, phi, theta, .0, 0
        );
        fPolylineEngine.AddPolylinesToMuonTracks(muonTrack, minimalisticMuonTrack);
        event.AddTrack(minimalisticMuonTrack);
    }
}
