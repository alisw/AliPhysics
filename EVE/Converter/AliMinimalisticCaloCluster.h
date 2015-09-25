//
// Created by jniedzie on 9/21/15.
//

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
    Float_t fR;
    Float_t fPhi;
    Float_t fEta;
    Float_t fPhiHalfLength;
    Float_t fEtaHalfLength;
    Float_t fEnergy;
};


#endif
