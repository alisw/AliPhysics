//
// Created by mgrochow on 7/31/15.
//


#ifndef ALIROOT_ALIMINIMALISTICTRACK_H
#define ALIROOT_ALIMINIMALISTICTRACK_H


#include <iostream>

#include <TObject.h>
#include <TString.h>


class AliMinimalisticTrack : public TObject {

public:
    AliMinimalisticTrack() : TObject()
    {  }

    AliMinimalisticTrack(
            Int_t charge, Double_t energy, Int_t ID, Int_t PID, Double_t mass, Double_t signedPT,
            Double_t startXYZ[], Double_t endXYZ[], Double_t pxpypz[], Int_t parentID, Double_t phi,
            Double_t theta, Double_t helixCurvature, Int_t type);

    void AddChild(Int_t childID);
    const std::vector<Int_t> &GetChildrenIDs() const;
    UInt_t GetID();

    void AddPolyPoint(Double_t x, Double_t y, Double_t z);
    void AddPolyPoint(Double_t xyz[3]);
    static const Int_t fgkImaginaryParent = -1;
    enum {
        kStandard, kKinkMother, kKinkDaughter, kV0NegativeDaughter, kV0PositiveDaughter, kV0Mother,
        kCascadePrimaryMother, kCascadePrimaryDaughter, kCascadeSecondaryMother,
        kCascadeNegativeDaughter, kCascadePositiveDaughter, kMuonMatched, kMuonNotMatched, kMuonGhost
    };
    void SetTrackType(Int_t type);
private:
    static const TString fgkTrackTypes[14];

    void AddStartCoordinates(Double_t xyz[3]);
    void AddEndCoordinates(Double_t xyz[3]);
    void AddMomentum(Double_t pxpypz[3]);

    TString fType;
    Int_t fCharge;
    Double_t fE;
    Int_t fParentID;
    Int_t fPID;
    Double_t fSignedPT;
    Double_t fMass;
    Double_t fMomentum[3];
    Double_t fStartCoordinates[3];
    Double_t fEndCoordinates[3];
    std::vector<Int_t> fChildrenIDs;
    Double_t fHelixCurvature;
    Double_t fTheta;
    Double_t fPhi;
    std::vector<Double_t> fPolyX;
    std::vector<Double_t> fPolyY;
    std::vector<Double_t> fPolyZ;


    ClassDef(AliMinimalisticTrack, 1);
};


#endif //ALIROOT_ALIMINIMALISTICTRACK_H
