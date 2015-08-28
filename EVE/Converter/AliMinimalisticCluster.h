//
// Created by mgrochow on 7/31/15.
//

#include <iostream>

#include <TObject.h>
#include <TString.h>

#ifndef ALIROOT_ALIMINIMALISTICCLUSTER_H
#define ALIROOT_ALIMINIMALISTICCLUSTER_H


class AliMinimalisticCluster : public TObject {
ClassDef(AliMinimalisticCluster, 1)
public:
    AliMinimalisticCluster(Int_t trackID);

    AliMinimalisticCluster() : TObject() { }

    void InsertXArray(const Float_t *x, Int_t numberOfPoints);
    void InsertYArray(const Float_t *y, Int_t numberOfPoints);
    void InsertZArray(const Float_t *z, Int_t numberOfPoints);
    void InsertValueDescription(const Char_t *description[], Int_t numberOfPoints);
    void InsertValue(const Int_t *value, Int_t numberOfPoints);
    void SetSource(TString source);
private:
    std::vector<Float_t> fXArray;
    std::vector<Float_t> fYArray;
    std::vector<Float_t> fZArray;
    std::vector<Int_t> fValue;
    std::vector<TString> fValueDescritpion;
    Int_t fTrackID;
    TString fSource;
};


#endif //ALIROOT_ALIMINIMALISTICCLUSTER_H
