//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cluster information                          //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDcluster.h"
#include <stdio.h>

ClassImp(AliPMDcluster)

AliPMDcluster::AliPMDcluster()
{
  // Default constructor
  for (Int_t i = 0; i < 5; i++)
    {
      fClusData[i] = 0.;
    }
}
AliPMDcluster::AliPMDcluster(Float_t *clusdata)
{
  // Constructor
  for (Int_t i = 0; i < 5; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
AliPMDcluster::AliPMDcluster(const AliPMDcluster &pmdcluster):TObject(pmdcluster)
{
  //Copy Constructor 
  if(&pmdcluster == this) return;
  for(Int_t i=0; i<5; i++)
    {
      this->fClusData[i] = pmdcluster.fClusData[i];
    }
  return;
}

AliPMDcluster & AliPMDcluster::operator=(const AliPMDcluster &pmdcluster)
{
  // Assignment operator 
  if(&pmdcluster == this) return *this;
  for(Int_t i=0; i<5; i++)
    {
      this->fClusData[i] = pmdcluster.fClusData[i];
    }
  return *this;
}

AliPMDcluster::~AliPMDcluster()
{
  // Destructor
}
Float_t AliPMDcluster::GetClusX() const
{
  return fClusData[0];
}
Float_t AliPMDcluster::GetClusY() const
{
  return fClusData[1];
}
Float_t AliPMDcluster::GetClusADC() const
{
  return fClusData[2];
}
Float_t AliPMDcluster::GetClusCells() const
{
  return fClusData[3];
}
Float_t AliPMDcluster::GetClusRadius() const
{
  return fClusData[4];
}


