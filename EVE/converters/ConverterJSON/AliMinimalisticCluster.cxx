//
// Created by mgrochow on 7/31/15.
//

#include "AliMinimalisticCluster.h"
ClassImp(AliMinimalisticCluster);

AliMinimalisticCluster::AliMinimalisticCluster(Int_t trackID) :
        TObject(), fTrackID(trackID)
{ }

void AliMinimalisticCluster::InsertXArray(const Float_t *x, Int_t numberOfPoints)
{
    fXArray.insert(fXArray.end(), x, x + numberOfPoints);
}

void AliMinimalisticCluster::InsertYArray(const Float_t *y, Int_t numberOfPoints)
{
    fYArray.insert(fYArray.end(), y, y + numberOfPoints);
}

void AliMinimalisticCluster::InsertZArray(const Float_t *z, Int_t numberOfPoints)
{
    fZArray.insert(fZArray.end(), z, z + numberOfPoints);
}

void AliMinimalisticCluster::InsertValueDescription(const Char_t *description[], Int_t numberOfPoints)
{
    fValueDescritpion.insert(fValueDescritpion.end(), description, description + numberOfPoints);
}

void AliMinimalisticCluster::InsertValue(const Int_t *value, Int_t numberOfPoints)
{
    fValue.insert(fValue.end(), value, value + numberOfPoints);
}

void AliMinimalisticCluster::SetSource(TString source)
{
    fSource = source;
}
