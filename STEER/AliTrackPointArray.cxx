/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//                          Class AliTrackPointArray                        //
//   This class contains the ESD track space-points which are used during   //
//   the alignment procedures. Each space-point consist of 3 coordinates    //
//   (and their errors) and the index of the sub-detector which contains    //
//   the space-point.                                                       //
//   cvetan.cheshkov@cern.ch 3/11/2005                                      //
//////////////////////////////////////////////////////////////////////////////

#include "AliTrackPointArray.h"

ClassImp(AliTrackPointArray)

//______________________________________________________________________________
AliTrackPointArray::AliTrackPointArray()
{
  fNPoints = fSize = 0;
  fX = fY = fZ = 0;
  fVolumeID = 0;
  fCov = 0;
}

//______________________________________________________________________________
AliTrackPointArray::AliTrackPointArray(Int_t npoints):
  fNPoints(npoints)
{
  // Constructor
  //
  fSize = 6*npoints;
  fX = new Float_t[npoints];
  fY = new Float_t[npoints];
  fZ = new Float_t[npoints];
  fVolumeID = new UShort_t[npoints];
  fCov = new Float_t[fSize];
}

//______________________________________________________________________________
AliTrackPointArray::AliTrackPointArray(const AliTrackPointArray &array):
  TObject(array)
{
  // Copy constructor
  //
  fNPoints = array.fNPoints;
  fSize = array.fSize;
  fX = new Float_t[fNPoints];
  fY = new Float_t[fNPoints];
  fZ = new Float_t[fNPoints];
  fVolumeID = new UShort_t[fNPoints];
  fCov = new Float_t[fSize];
  memcpy(fX,array.fX,fNPoints*sizeof(Float_t));
  memcpy(fY,array.fY,fNPoints*sizeof(Float_t));
  memcpy(fZ,array.fZ,fNPoints*sizeof(Float_t));
  memcpy(fVolumeID,array.fVolumeID,fNPoints*sizeof(UShort_t));
  memcpy(fCov,array.fCov,fSize*sizeof(Float_t));
}

//_____________________________________________________________________________
AliTrackPointArray &AliTrackPointArray::operator =(const AliTrackPointArray& array)
{
  // assignment operator
  //
  if(this==&array) return *this;
  ((TObject *)this)->operator=(array);

  fNPoints = array.fNPoints;
  fSize = array.fSize;
  fX = new Float_t[fNPoints];
  fY = new Float_t[fNPoints];
  fZ = new Float_t[fNPoints];
  fVolumeID = new UShort_t[fNPoints];
  fCov = new Float_t[fSize];
  memcpy(fX,array.fX,fNPoints*sizeof(Float_t));
  memcpy(fY,array.fY,fNPoints*sizeof(Float_t));
  memcpy(fZ,array.fZ,fNPoints*sizeof(Float_t));
  memcpy(fVolumeID,array.fVolumeID,fNPoints*sizeof(UShort_t));
  memcpy(fCov,array.fCov,fSize*sizeof(Float_t));

  return *this;
}

//______________________________________________________________________________
AliTrackPointArray::~AliTrackPointArray()
{
  // Destructor
  //
  delete [] fX;
  delete [] fY;
  delete [] fZ;
  delete [] fVolumeID;
  delete [] fCov;
}


//______________________________________________________________________________
Bool_t AliTrackPointArray::AddPoint(Int_t i, const AliTrackPoint *p)
{
  // Add a point to the array at position i
  //
  if (i >= fNPoints) return kFALSE;
  fX[i] = p->GetX();
  fY[i] = p->GetY();
  fZ[i] = p->GetZ();
  fVolumeID[i] = p->GetVolumeID();
  memcpy(&fCov[6*i],p->GetCov(),6*sizeof(Float_t));
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTrackPointArray::GetPoint(AliTrackPoint &p, Int_t i) const
{
  // Get the point at position i
  //
  if (i >= fNPoints) return kFALSE;
  p.SetXYZ(fX[i],fY[i],fZ[i],&fCov[6*i]);
  p.SetVolumeID(fVolumeID[i]);
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTrackPointArray::HasVolumeID(UShort_t volid) const
{
  // This method checks if the array
  // has at least one hit in the detector
  // volume defined by volid
  Bool_t check = kFALSE;
  for (Int_t ipoint = 0; ipoint < fNPoints; ipoint++)
    if (fVolumeID[ipoint] == volid) check = kTRUE;

  return check;
}

ClassImp(AliTrackPoint)

//______________________________________________________________________________
AliTrackPoint::AliTrackPoint()
{
  // Default constructor
  //
  fX = fY = fZ = 0;
  fVolumeID = 0;
  memset(fCov,0,6*sizeof(Float_t));
}


//______________________________________________________________________________
AliTrackPoint::AliTrackPoint(Float_t x, Float_t y, Float_t z, const Float_t *cov, UShort_t volid)
{
  // Constructor
  //
  SetXYZ(x,y,z,cov);
  SetVolumeID(volid);
}

//______________________________________________________________________________
AliTrackPoint::AliTrackPoint(const Float_t *xyz, const Float_t *cov, UShort_t volid)
{
  // Constructor
  //
  SetXYZ(xyz[0],xyz[1],xyz[2],cov);
  SetVolumeID(volid);
}

//______________________________________________________________________________
AliTrackPoint::AliTrackPoint(const AliTrackPoint &p):
  TObject(p)
{
  // Copy constructor
  //
  SetXYZ(p.fX,p.fY,p.fZ,&(p.fCov[0]));
  SetVolumeID(p.fVolumeID);
}

//_____________________________________________________________________________
AliTrackPoint &AliTrackPoint::operator =(const AliTrackPoint& p)
{
  // assignment operator
  //
  if(this==&p) return *this;
  ((TObject *)this)->operator=(p);

  SetXYZ(p.fX,p.fY,p.fZ,&(p.fCov[0]));
  SetVolumeID(p.fVolumeID);

  return *this;
}

//______________________________________________________________________________
void AliTrackPoint::SetXYZ(Float_t x, Float_t y, Float_t z, const Float_t *cov)
{
  // Set XYZ coordinates and their cov matrix
  //
  fX = x;
  fY = y;
  fZ = z;
  if (cov)
    memcpy(fCov,cov,6*sizeof(Float_t));
}

//______________________________________________________________________________
void AliTrackPoint::SetXYZ(const Float_t *xyz, const Float_t *cov)
{
  // Set XYZ coordinates and their cov matrix
  //
  SetXYZ(xyz[0],xyz[1],xyz[2],cov);
}

//______________________________________________________________________________
void AliTrackPoint::GetXYZ(Float_t *xyz, Float_t *cov) const
{
  xyz[0] = fX;
  xyz[1] = fY;
  xyz[2] = fZ;
  if (cov)
    memcpy(cov,fCov,6*sizeof(Float_t));
}
