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


const TString AliExternalFormatConverter::fgkDetector[23] = {
        "Invalid Layer", "First Layer", "SPD1", "SPD2", "SDD1", "SDD2", "SSD1", "SSD2", "TPC1", "TPC2",
        "TRD1", "TRD2", "TRD3", "TRD4", "TRD5", "TRD6", "TOF", "PHOS1", "PHOS2", "HMPID", "MUON", "EMCAL", "LastLayer"
};

AliExternalFormatConverter::AliExternalFormatConverter()
        : fESDFile(nullptr), fESDFriend(nullptr), fESDEvent(nullptr), fESDTree(nullptr), fApp(nullptr)
{ }


AliExternalFormatConverter::AliExternalFormatConverter(const TString dirPath)
        : fESDFile(nullptr), fESDFriend(nullptr), fESDEvent(nullptr), fESDTree(nullptr)
{
    char *argv = "aaaaaaaaaa";
    int a = 0;
    fApp = new TRint("App", &a, &argv, 0 , 0, kTRUE);
    TEveManager::Create(kFALSE);
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
    //fApp->Run();
    //fApp->Terminate(1);
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

        AliMinimalisticTrack negative = GenerateMinimalisticTrack(negativeID, v0ParentID);
        AliMinimalisticTrack positive = GenerateMinimalisticTrack(positiveID, v0ParentID);
        if (fESDFriend){
            AliMinimalisticCluster motherCluster = GenerateMinimalisticCluster(negativeID);
            AliMinimalisticCluster daughterCluster = GenerateMinimalisticCluster(positiveID);
            event.AddCluster(motherCluster);
            event.AddCluster(daughterCluster);
        }
        AddPolyLinesToV0Track(v0Entry, negative, positive);
        Double_t startCOORDS[] = {.0, .0, .0};
        AliMinimalisticTrack V0Parenttrack = GenerateMinimalisticV0ParentTrack(
                v0, negativeID, positiveID, v0ParentID, startCOORDS
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
        AliMinimalisticTrack negative = GenerateMinimalisticTrack(negativeID, v0ParentID);
        AliMinimalisticTrack positive = GenerateMinimalisticTrack(positiveID, v0ParentID);
        AliMinimalisticTrack bachelor = GenerateMinimalisticTrack(bachelorID, cascadeParentID);
        AddPolylinesToCascade(cascadeEntry, negative, positive, bachelor);
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
        AliMinimalisticTrack V0Parenttrack = GenerateMinimalisticV0ParentTrack(
                (AliESDv0*)cascade, negativeID, positiveID, v0ParentID, v0StartCoords, cascadeParentID);
        event.AddTrack(V0Parenttrack);

        AliMinimalisticTrack cascadeParentTrack = GenerateMinimalisticCascadeParenTrack(
                cascade, v0ParentID, bachelorID, cascadeParentID);

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

        AliMinimalisticTrack daughter = GenerateMinimalisticTrack(daughterID, motherID);
        AliMinimalisticTrack mother = GenerateMinimalisticTrack(
                motherID,
                AliMinimalisticTrack::fgkImaginaryParent
        );
        if (fESDFriend){
            AliMinimalisticCluster motherCluster = GenerateMinimalisticCluster(motherID);
            AliMinimalisticCluster daughterCluster = GenerateMinimalisticCluster(daughterID);
            event.AddCluster(motherCluster);
            event.AddCluster(daughterCluster);
        }
        AddPolyLinesToKinkTrack(kinkEntry, mother, daughter);
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
    AliMinimalisticTrack track = GenerateMinimalisticTrack(trackID, parentID);
    if (childID!=AliMinimalisticTrack::fgkImaginaryParent)
        track.AddChild(childID);
    if (fESDFriend){
        AliMinimalisticCluster cluster = GenerateMinimalisticCluster(trackID);
        event.AddCluster(cluster);
    }
    AddPolylinesToMinimalisticTrack(trackID, track);
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
        Int_t trackNumber, Int_t parentID
) const
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
            charge, energy, trackNumber, PID, mass, signedPT, startXYZ, endXYZ, pxpypz, parentID, phi, theta,
            helixCurvature
    );

    return minimalisticTrack;
}

void AliExternalFormatConverter::AddPolylinesToMinimalisticTrack(
        Int_t trackID, AliMinimalisticTrack &minimalisticTrack) const
{
    Int_t maxTrackRadius = 520;
    TEveTrackList* tEveTrackList = new TEveTrackList("ESD-Tracks");
    TEveTrackPropagator *pTrackPropagator = tEveTrackList->GetPropagator();

    pTrackPropagator->SetMagField(-fESDEvent->GetMagneticField() / 10.0);
    pTrackPropagator->SetMaxR(maxTrackRadius);
    AliESDtrack *track = fESDEvent->GetTrack(trackID);
    AliEveTrack* evetrack = new AliEveTrack(track, pTrackPropagator);
    evetrack->SetSourceObject(track);
    // Add inner/outer track parameters as path-marks.
    if (track->IsOn(AliESDtrack::kTPCrefit))
    {
        Double_t pbuf[3], vbuf[3];
        if (track->GetInnerParam() != 0) {
            track->GetInnerParam()->GetXYZ(vbuf);
            track->GetInnerParam()->GetPxPyPz(pbuf);
        }
        if (track->GetOuterParam() != 0) {
            track->GetOuterParam()->GetXYZ(vbuf);
            track->GetOuterParam()->GetPxPyPz(pbuf);
        }
        TEvePathMark pm(TEvePathMark::kReference);
        pm.fV.Set(vbuf);
        pm.fP.Set(pbuf);
        evetrack->AddPathMark(pm);
    }
    tEveTrackList->SetChildClass(evetrack->Class());
    tEveTrackList->AddElement(evetrack);
    evetrack->MakeTrack();
    tEveTrackList->MakeTracks();

    std::vector<TEveVector4D> pointsVec = pTrackPropagator->GetPoints();
    for(std::vector<TEveVector4D>::iterator iter = pointsVec.begin(); iter != pointsVec.end(); ++iter){
        minimalisticTrack.AddPolyPoint(*iter);
    }
}

AliMinimalisticTrack AliExternalFormatConverter::GenerateMinimalisticV0ParentTrack(
        AliESDv0 *V0, Int_t negativeChild, Int_t positiveChild, Int_t myID, Double_t startXYZ[3], Int_t parentID
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
            charge, energy, myID, PID, mass, signedPt, startXYZ, endXYZ, PxPyPz, parentID, phi, theta, helixCurvature
    );
    parent.AddChild(negativeChild);
    parent.AddChild(positiveChild);
    parent.AddPolyPoint(startXYZ);
    parent.AddPolyPoint(endXYZ);
    return parent;
}

AliMinimalisticTrack AliExternalFormatConverter::GenerateMinimalisticCascadeParenTrack(
        AliESDcascade *cascade, Int_t v0ChildID, Int_t singleChildID, Int_t myID) const
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
    Int_t parentID = AliMinimalisticTrack::fgkImaginaryParent;

    AliMinimalisticTrack cascdeParent(
            charge, energy, myID, PID, mass, signedPt, startXYZ, endXYZ, PxPyPz, parentID, phi, theta, helixCurvature
    );
    cascdeParent.AddChild(v0ChildID);
    cascdeParent.AddChild(singleChildID);
    cascdeParent.AddPolyPoint(startXYZ);
    cascdeParent.AddPolyPoint(endXYZ);
    return cascdeParent;
}

void AliExternalFormatConverter::ExtractTrackPointArrays(
        AliMinimalisticCluster &cluster, Int_t trackNumber
) const
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

void AliExternalFormatConverter::PopulateEventWithMuonTracks(AliMinimalisticEvent &event) const
{
    AliEveMUONTrackList* lt = new AliEveMUONTrackList("ESD-Tracks");
    TEveRecTrack rt;
    AliESDMuonTrack *mTrack;
    for (Int_t muonTrack = 0; muonTrack < fESDEvent->GetNumberOfMuonTracks(); muonTrack++){
        mTrack = fESDEvent->GetMuonTrack(muonTrack);
        if (mTrack->GetNHit() == 0) continue;
        rt.fLabel = muonTrack;
        AliEveMUONTrack* track = new AliEveMUONTrack(&rt, lt->GetPropagator());
        track->MakeESDTrack(mTrack);
        lt->AddElement(track);

        Int_t charge = mTrack->Charge();
        Double_t energy = mTrack->E();
        Int_t parentID = AliMinimalisticTrack::fgkImaginaryParent;
        Double_t signedPT = mTrack->Pt()*charge;
        Double_t mass = mTrack->M();
        Double_t PxPyPz[3]; mTrack->PxPyPz(PxPyPz);
        Double_t startXYZ[3];
        startXYZ[0] = mTrack->GetNonBendingCoorAtDCA();
        startXYZ[1] = mTrack->GetBendingCoorAtDCA();
        startXYZ[2] = mTrack->GetZ();
        Double_t theta = mTrack->Theta();
        Double_t phi = mTrack->Phi();
    }
      lt->HackMomentumLimits();
}

void AliExternalFormatConverter::AddPolyLinesToKinkTrack(
        Int_t kinkID, AliMinimalisticTrack &mTrack, AliMinimalisticTrack &dTrack
) const
{
    AliEveKinkList* cont = new AliEveKinkList("ESD kink");

    TEveTrackPropagator* rnrStyleMoth = cont->GetPropagatorMoth();
    TEveTrackPropagator* rnrStyleDaugh = cont->GetPropagatorDaugh();

    rnrStyleMoth->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStyleDaugh->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStyleDaugh->SetMaxR(520);

    AliESDkink *kink = fESDEvent->GetKink(kinkID);
    AliESDtrack* moth = fESDEvent->GetTrack(kink->GetIndex(0));
    AliESDtrack* daug = fESDEvent->GetTrack(kink->GetIndex(1));

    TEveRecTrack rcMoth;
    TEveRecTrack rcDaug;
    TEveRecKink rcKink;

    rcKink.fPMother.Set(kink->GetMotherP());
    rcKink.fPDaughter.Set(kink->GetDaughterP());
    rcKink.fVKink.Set(kink->GetPosition());

    Double_t pbuf[3], vbuf[3];

    rcMoth.fSign = moth->GetTPCInnerParam()->GetSign();
    moth->GetTPCInnerParam()->GetXYZ(vbuf); rcMoth.fV.Set(vbuf);
    moth->GetTPCInnerParam()->GetPxPyPz(pbuf); rcMoth.fP.Set(pbuf);

    rcDaug.fSign = daug->GetOuterParam()->GetSign();
    rcDaug.fV.Set(rcKink.fVKink);
    rcDaug.fP.Set(rcKink.fPDaughter);

    AliEveKink* myKink = new AliEveKink(&rcMoth, &rcDaug, &rcKink, rnrStyleMoth,rnrStyleDaugh);
    myKink->SetESDKinkIndex(kinkID);
    gEve->AddElement(myKink, cont);
    cont->MakeKinks();

    std::vector<TEveVector4D> mPoints = rnrStyleMoth->GetPoints();
    std::vector<TEveVector4D> dPoints = rnrStyleDaugh->GetPoints();
    std::cout << "Przed funkcja wielkosc: " << mPoints.size() << std::endl;
    InsertPolyPoints(mTrack, mPoints);
    InsertPolyPoints(dTrack, dPoints);
}

void AliExternalFormatConverter::AddPolyLinesToV0Track(
        Int_t v0ID, AliMinimalisticTrack &negativeTrack, AliMinimalisticTrack &positiveTrack) const
{
    AliEveV0List* cont = new AliEveV0List("ESD v0");
    TEveTrackPropagator* positivePropagator = cont->GetPropagatorPos();
    TEveTrackPropagator* negativePropagator = cont->GetPropagatorNeg();
    positivePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );
    negativePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );

    gEve->AddElement(cont);
    AliESDv0 *v0 = fESDEvent->GetV0(v0ID);

    TEveRecTrack rcPos;
    TEveRecTrack rcNeg;
    TEveRecV0 rcV0;

    Double_t pbuf[3], vbuf[3];

    rcNeg.fSign = v0->GetParamN()->GetSign();
    v0->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
    v0->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);
    rcPos.fSign = v0->GetParamP()->GetSign();
    v0->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
    v0->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);

    rcNeg.fIndex = v0->GetNindex();
    rcPos.fIndex = v0->GetPindex();

    AliEveV0* myV0 = new AliEveV0(&rcNeg, &rcPos, &rcV0, negativePropagator, positivePropagator);

    myV0->SetESDIndex(v0ID);
    myV0->SetOnFlyStatus(v0->GetOnFlyStatus());
    myV0->SetDaughterDCA(v0->GetDcaV0Daughters());

    gEve->AddElement(myV0, cont);

    cont->MakeV0s();

    std::vector<TEveVector4D> negPoints = negativePropagator->GetPoints();
    std::vector<TEveVector4D> posPoints = positivePropagator->GetPoints();
    std::cout << "wielkosc neg: " << negPoints.size() << std::endl;
    InsertPolyPoints(negativeTrack, negPoints);
    InsertPolyPoints(positiveTrack, posPoints);
}

void AliExternalFormatConverter::AddPolylinesToCascade(
        Int_t cascadeID,
        AliMinimalisticTrack &negativeTrack,
        AliMinimalisticTrack &positiveTrack,
        AliMinimalisticTrack &bachelorTrack) const
{
    AliESDVertex* primVtx = (AliESDVertex*) fESDEvent->GetPrimaryVertex();

    AliEveCascadeList* cont = new AliEveCascadeList("ESD cascade");
    TEveTrackPropagator* rnrStyleBac = cont->GetPropagatorBac();
    TEveTrackPropagator* rnrStyleNeg = cont->GetPropagatorNeg();
    TEveTrackPropagator* rnrStylePos = cont->GetPropagatorPos();
    rnrStyleBac->SetMagField( 0.1 * fESDEvent->GetMagneticField() );
    rnrStyleNeg->SetMagField( 0.1 * fESDEvent->GetMagneticField() );
    rnrStylePos->SetMagField( 0.1 * fESDEvent->GetMagneticField() );

    gEve->AddElement(cont);

    AliESDcascade *cascade = fESDEvent->GetCascade(cascadeID);

    TEveRecTrack rcPos,rcNeg,rcBac;
    TEveRecV0 rcV0;
    TEveRecCascade rcCascade;

    Double_t v[3],pBac[3]={0.}, pNeg[3]={0.}, pPos[3]={0.}, cv[21]={0.},pbuf[3], vbuf[3];

    cascade->GetBPxPyPz(pBac[0], pBac[1], pBac[2]);
    cascade->GetNPxPyPz(pNeg[0], pNeg[1], pNeg[2]);
    cascade->GetPPxPyPz(pPos[0], pPos[1], pPos[2]);
    cascade->GetXYZcascade(v[0], v[1], v[2]);

    AliExternalTrackParam *bacParam = new AliExternalTrackParam(v,pBac,cv,cascade->Charge());

    rcBac.fSign = bacParam->GetSign();
    rcBac.fIndex = cascade->GetBindex();
    bacParam->GetXYZ(vbuf); rcBac.fV.Set(vbuf);
    bacParam->GetPxPyPz(pbuf); rcBac.fP.Set(pbuf);

    rcNeg.fSign = cascade->GetParamN()->GetSign();
    rcNeg.fIndex = cascade->GetNindex();
    cascade->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
    cascade->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);

    rcPos.fSign = cascade->GetParamP()->GetSign();
    rcPos.fIndex = cascade->GetPindex();
    cascade->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
    cascade->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);

    AliEveCascade* myCascade = new AliEveCascade(&rcBac, &rcNeg, &rcPos, &rcV0, &rcCascade, rnrStyleBac, rnrStyleNeg, rnrStylePos);

    myCascade->SetESDIndex(cascadeID);
    myCascade->SetDaughterDCA(cascade->GetDcaXiDaughters());
    myCascade->SetLambdaP( pNeg[0]+pPos[0], pNeg[1]+pPos[1], pNeg[2]+pPos[2] );
    myCascade->SetBachP( pBac[0], pBac[1], pBac[2]);
    
    gEve->AddElement(myCascade, cont);
    cont->MakeCascades();

    std::vector<TEveVector4D> negPoints = rnrStyleNeg->GetPoints();
    std::vector<TEveVector4D> posPoints = rnrStylePos->GetPoints();
    std::vector<TEveVector4D> bacPoints = rnrStyleBac->GetPoints();
    InsertPolyPoints(negativeTrack, negPoints);
    InsertPolyPoints(positiveTrack, posPoints);
    InsertPolyPoints(bachelorTrack, bacPoints);
}

void AliExternalFormatConverter::InsertPolyPoints(
        AliMinimalisticTrack &Track, std::vector<TEveVector4D> &Points) const
{
    std::cout <<"po funkcji: " << Points.size() <<std::endl;
    for(std::vector<TEveVector4D>::iterator iter = Points.begin(); iter != Points.end(); ++iter){
        Track.AddPolyPoint(*iter);
    }
}
