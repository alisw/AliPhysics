//
// Created by jniedzie on 9/21/15.
//

#include "AliMinimalisticCaloCluster.h"
ClassImp(AliMinimalisticCaloCluster);

AliMinimalisticCaloCluster::AliMinimalisticCaloCluster(
        float r, float phi, float eta, float phiHalfLength, float etaHalfLength, float energy
) :
    TObject(),
    fR(r),
    fPhi(phi),
    fEta(eta),
    fPhiHalfLength(phiHalfLength),
    fEtaHalfLength(etaHalfLength),
    fEnergy(energy)
{ }