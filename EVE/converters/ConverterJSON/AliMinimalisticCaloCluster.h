/// \class AliMinimalisticCaloCluster
/// In this data containers basic information from calorimeters is stored. Calorime-
/// ter cluster is represented as a rectangular cuboid.
///
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology

#include <iostream>

#include <TObject.h>
#include <TString.h>

#ifndef ALIROOT_ALIMINIMALISTICCALOCLUSTER_H
#define ALIROOT_ALIMINIMALISTICCALOCLUSTER_H


class AliMinimalisticCaloCluster : public TObject
{
    ClassDef(AliMinimalisticCaloCluster, 1)
public:
    AliMinimalisticCaloCluster(float r,float phi, float eta, float phiHalfLength, float etaHalfLength, float energy);
    
    AliMinimalisticCaloCluster() : TObject() { }
private:
    /// Coordinates of the position in space in the spherical coordinate system:
    Float_t fR; /// radial distance to the middle of the cluster
    Float_t fPhi; /// azimuthal angle
    Float_t fEta; /// polar angle
    /// Dimension of faces:
    Float_t fPhiHalfLength; /// width
    Float_t fEtaHalfLength; /// length
    Float_t fEnergy; /// the height of the cuboid is proportional to the measured energy
};


#endif
