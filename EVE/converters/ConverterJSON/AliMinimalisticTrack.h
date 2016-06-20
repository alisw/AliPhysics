/// \class AliMinimalisticTrack
/// Minimalistic track is a specialized container storing extracted data regarding
/// reconstructed tracks and relations between them. Relations between particles are
/// strictly hierarchical, hence a structure that allows storing data in ranked order is
/// required. A data type called a tree meets these requirements.
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

#ifndef ALIROOT_ALIMINIMALISTICTRACK_H
#define ALIROOT_ALIMINIMALISTICTRACK_H


#include <iostream>

#include <TString.h>
#include <TObject.h>

#include <ConversionConstants.h>


class AliMinimalisticTrack : public TObject {
ClassDef(AliMinimalisticTrack, 1);
public:
    AliMinimalisticTrack() : TObject()
    {  }

    AliMinimalisticTrack(
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
    );

    void AddChild(Int_t childID);
    void AddPolyPoint(Double_t x, Double_t y, Double_t z);
    void AddPolyPoint(Double_t xyz[3]);
    void SetTrackType(TrackType type);

    static const Int_t fgkNoParent = -1;
    static const Int_t fgkNoChild = -1;
private:
    void AddStartCoordinates(Double_t xyz[3]);
    void AddEndCoordinates(Double_t xyz[3]);
    void AddMomentum(Double_t pxpypz[3]);

    TString fType; /// Type (standard, V0 mother, daughter etc.)
    Int_t fCharge; /// Charge of the particle
    Double_t fE; /// Energy of the particle
    Int_t fParentID; /// Identification number of the parent-track -- it is '-1' (fgkNoParent) if there is no parent
    Int_t fPID; /// Type of the particle
    Double_t fSignedPT;
    Double_t fMass; /// Mass of the particle
    Double_t fMomentum[3]; /// Momentum represented as a vector [p_x , p_y , p_z]
    Double_t fStartCoordinates[3]; /// Start Coordinates represented as a vector [x, y, z] -- it is important concerning decaying particles!
    Double_t fEndCoordinates[3]; /// End Coordinates represented as a vector [x, y, z] -- it is important concerning decaying particles!
    std::vector<Int_t> fChildrenIDs; /// Identification numbers of the child-tracks
    /// Helix parameters – parameters describing trajectory of the particle
    Double_t fHelixCurvature; /// Helix curvature of the trajectory
    Double_t fTheta; /// An angle from Z-axis to the radius vector pointing to the particle
    Double_t fPhi; /// An angle from X-axis to the radius vector pointing to the particle
    /// Polylines -- array of points along the trajectory of the track
    std::vector<Double_t> fPolyX;
    std::vector<Double_t> fPolyY;
    std::vector<Double_t> fPolyZ;
};


#endif //ALIROOT_ALIMINIMALISTICTRACK_H
