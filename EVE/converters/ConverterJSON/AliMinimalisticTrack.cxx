//
// Created by mgrochow on 7/31/15.
//

#include "AliMinimalisticTrack.h"

ClassImp(AliMinimalisticTrack)

AliMinimalisticTrack::AliMinimalisticTrack(
        Int_t charge,
        Double_t energy,
        Int_t ID,
        Int_t PID,
        Double_t mass,
        Double_t signedPT,
        Double_t startXYZ[],
        Double_t endXYZ[],
        Double_t pxpypz[],
        Int_t parentID,
        Double_t phi,
        Double_t theta,
        Double_t helixCurvature,
        Int_t type
)
        :
        fCharge(charge),
        fE(energy),
        fParentID(parentID),
        fMass(mass),
        fPID(PID),
        fSignedPT(signedPT),
        fPhi(phi),
        fTheta(theta),
        fHelixCurvature(helixCurvature),
        TObject()
{
    AddMomentum(pxpypz);
    AddStartCoordinates(startXYZ);
    AddEndCoordinates(endXYZ);
    SetUniqueID(ID);
    fType = fgkTrackTypes[type];
}

void AliMinimalisticTrack::AddChild(Int_t childID)
{
    fChildrenIDs.push_back(childID);
}

void AliMinimalisticTrack::AddMomentum(Double_t pxpypz[3])
{
    for (int i = 0; i < 3; i++)
        fMomentum[i] = pxpypz[i];
}

void AliMinimalisticTrack::AddStartCoordinates(Double_t xyz[3])
{
    for (int i = 0; i < 3; i++)
        fStartCoordinates[i] = xyz[i];
}

void AliMinimalisticTrack::AddEndCoordinates(Double_t xyz[3])
{
    for (int i = 0; i < 3; i++)
        fEndCoordinates[i] = xyz[i];
}

void AliMinimalisticTrack::AddPolyPoint(Double_t x, Double_t y, Double_t z)
{
    fPolyX.push_back(x);
    fPolyY.push_back(y);
    fPolyZ.push_back(z);
}

void AliMinimalisticTrack::AddPolyPoint(Double_t *xyz)
{
    fPolyX.push_back(xyz[0]);
    fPolyY.push_back(xyz[1]);
    fPolyZ.push_back(xyz[2]);
}

void AliMinimalisticTrack::SetTrackType(TrackType type)
{
    fType = fgkTrackTypes[type];
}
