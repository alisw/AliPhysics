/// \class AliMinimalisticCluster
/// \Prepare data for saving as JSON/XML
/// This class stores all essential data regarding recorded particle trajectory in the
/// detector. Data is stored in structure of arrays
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

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
    std::vector<Float_t> fXArray; /// X position of registered signal in space
    std::vector<Float_t> fYArray; /// Y position of registered signal in space
    std::vector<Float_t> fZArray; /// Z position of registered signal in space
    std::vector<Int_t> fValue;
    std::vector<TString> fValueDescritpion; /// Name of the detector in which given point was registered (e.g. TOF)
    Int_t fTrackID; /// id of the track that the cluster corresponds to
    TString fSource; /// Type of ile from which it was extracted -- for the time being ESD
};


#endif //ALIROOT_ALIMINIMALISTICCLUSTER_H
