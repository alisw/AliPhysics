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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//  Reconstructed space point class for set:ITS   
//  Reconstructed points are expressed simultaneously in two different 
//  reference frames, both differing from the global system.
//  The first is referred to the sensor (see AliITSsegmentation for the
//  definition) and each point is represented by two coordinates: fXloc and
//  fZloc. This system in the code is referred to as "local"
//  The second is used for tracking (V2, SA and MI versions) and the X axis 
//  represents the radial coordinate (this system is, in the bending plane, 
//  a rotated system w.r.t. the global reference system). 
//  Each reaconstructed point is represented by two coordinates: fY and fZ, 
//  inherited from AliCluster. This system in the code is referred to as 
//  "trackingV2".
///////////////////////////////////////////////////////////////////////////////

#include <TGeoMatrix.h>
#include "AliITSRecPoint.h"
#include "AliAlignObj.h"

ClassImp(AliITSRecPoint)

//_____________________________________________________________
AliITSRecPoint::AliITSRecPoint(): 
AliCluster(),
fXloc(0),
fZloc(0),
fdEdX(0),
fIndex(0),
fQ(0),
fLayer(0),
fNz(0),
fNy(0),
fChargeRatio(0),
fType(0),
fDeltaProb(0),
fDriftTime(0.),
fDriftSide(0)
{
    // default constructor
}

//________________________________________________________________________
AliITSRecPoint::AliITSRecPoint(Int_t *lab,Float_t *hit, Int_t *info, Bool_t local):
AliCluster(AliGeomManager::LayerToVolUID((info[2]+AliGeomManager::kSPD1),lab[3]&0x3FF),hit,0,0,lab),
fXloc(0),
fZloc(0),
fdEdX(0),
fIndex(lab[3]),
fQ(hit[4]),
fLayer(info[2]),
fNz(info[1]),
fNy(info[0]),
fChargeRatio(0),
fType(0),
fDeltaProb(0),
fDriftTime(0.),
fDriftSide(0)
{
  //standard constructor used in AliITSClusterFinderV2

  if (!local) { // Cluster V2
    Double_t txyz[3] = {GetX(), GetY(), GetZ()};
    Double_t lxyz[3] = {0, 0, 0};
    GetTracking2LocalMatrix()->LocalToMaster(txyz,lxyz);
    fXloc = lxyz[0]; fZloc = lxyz[2];
    if(fLayer==4) hit[5]=-hit[5];
    if( (fLayer==4) || (fLayer==5) ) SetSigmaYZ(hit[5]); 
  }
  else {
    switch (fLayer) {
    case 0:
    case 1:
      fdEdX = 0;
      break;
    case 2:
    case 3:
      fdEdX=fQ*1e-6;
      break;
    case 4:
      fdEdX=fQ*2.16;
      SetSigmaYZ(hit[5]);
    case 5:
      fdEdX=fQ*2.16;
      hit[5]=-hit[5];
      SetSigmaYZ(hit[5]);
      break;
    default:
      AliError(Form("Wrong ITS layer %d (0 -> 5)",fLayer));
      break;
    }
    fXloc = hit[0];
    fZloc = hit[1];
    Double_t lxyz[3] = {fXloc, 0, fZloc};
    Double_t txyz[3] = {0, 0, 0};
    GetTracking2LocalMatrix()->MasterToLocal(lxyz,txyz);

    SetX(0.); SetY(txyz[1]); SetZ(txyz[2]);

  }

}

//_______________________________________________________________________
AliITSRecPoint::AliITSRecPoint(const AliITSRecPoint& pt):
  AliCluster(pt),
  fXloc(pt.fXloc),
  fZloc(pt.fZloc),
  fdEdX(pt.fdEdX),
  fIndex(pt.fIndex),
  fQ(pt.fQ),
  fLayer(pt.fLayer),
  fNz(pt.fNz),
  fNy(pt.fNy),
  fChargeRatio(pt.fChargeRatio),
  fType(pt.fType),
  fDeltaProb(pt.fDeltaProb),
  fDriftTime(pt.fDriftTime),
  fDriftSide(pt.fDriftSide)
{
  //Copy constructor

}

//______________________________________________________________________
AliITSRecPoint& AliITSRecPoint::operator=(const AliITSRecPoint& source){
  // Assignment operator

  this->~AliITSRecPoint();
  new(this) AliITSRecPoint(source);
  return *this;

}

//----------------------------------------------------------------------
void AliITSRecPoint::Print(ostream *os){
    ////////////////////////////////////////////////////////////////////////
    // Standard output format for this class.
    ////////////////////////////////////////////////////////////////////////
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
 
    fmt = os->setf(ios::fixed);  // set fixed floating point output
    *os << GetLabel(0) << " " << GetLabel(1) << " " << GetLabel(2) << " ";
    fmt = os->setf(ios::scientific); // set scientific for dEdX.
    *os << GetX() <<" " << GetY() << " " << GetZ() << " " ;
    *os << GetSigmaY2() << " " << GetSigmaZ2() << " " << GetSigmaYZ() << " ";
    fmt = os->setf(ios::fixed);
    *os << GetVolumeId() << " "<< Misalign() /*fIsMisaligned*/ << " ";
    fmt = os->setf(ios::scientific); // set scientific for dEdX.
    *os << fXloc << " " << fZloc << " " << fdEdX << " ";
    fmt = os->setf(ios::fixed); // every fixed
    *os << fIndex <<" " << fQ << " "<<fLayer <<" "<<fNz<<" "<<fNy<<" ";
    *os << fChargeRatio<<" " << fType << " " << fDeltaProb << " " << fDriftTime<< " " << fDriftSide;
    os->flags(fmt); // reset back to old formating.
    return;
}

//----------------------------------------------------------------------
Int_t AliITSRecPoint::GetNpixels() const {
//
// returns the number of pixels used for the SPD clusters
//

 if(fLayer > 1) return -1;
 else return fType;

}

//----------------------------------------------------------------------
Int_t AliITSRecPoint::GetSPDclusterType() const {
//
// returns an Int_t with encoded information on cluster size
//    type <= 16: cluster type identifier according to conventional numbering
//    type >  16: Npixels+1000*Ny+1000000*Nz
//

 Int_t type = -1;
 if(fLayer > 1) return type;
 else {

   switch (fType) {
     case 1 : type = 1 ;break;
     case 2 : if(fNy == 2) type = 2;
              else         type = 3;
              break;
     case 3 : if(fNy == 3)      type = 4;
              else if(fNz == 3) type = 6;
              else              type = 5;
              break;
     case 4 : if(fNz == 1)                  type = 7;
              else if(fNz == 2 && fNy == 2) type = 8;
              else if(fNy == 2 && fNz == 3) type = 11;
              else if(fNy == 3 && fNz == 2) type = 9;
              else                          type = 15;
              break;
     case 5 : if(fNy == 3 && fNz == 2)  type = 10;
              if(fNy == 2 && fNz == 3 ) type = 12;
              if(fNy == 5)              type = 16;
              else                      type = fType+1000*fNy+1000000*fNz;
              break; 
     case 6 : if(fNy ==3 && fNz == 2) type = 13;
              if(fNy ==2 && fNz == 3) type = 14;
              else                    type = fType+1000*fNy+1000000*fNz;
              break; 
     default: type = fType+1000*fNy+1000000*fNz;
              break; 
   }  

   return type;
 }
}

//----------------------------------------------------------------------
Int_t AliITSRecPoint::GetSDDclusterType() const {
// returns an Int_t with encoded information on cluster size
// Byte1 = fNz Byte0=fNy, other two bytes empty for extra information
// max. allowed cluster size = 255
  Int_t typ=(fNz&0xFF)<<8;
  typ+=fNy&0xFF;
  if(fDriftSide==1) typ+=1<<16;
  return typ;
}
//----------------------------------------------------------------------
void  AliITSRecPoint::DecodeSDDclusterType(Int_t cluType, Int_t &cluSizAn, Int_t& cluSizTb, Int_t &drSide){
// Extract cluster sizes and drift side from cluster type 
  cluSizTb=cluType&0xFF;
  cluSizAn=(cluType>>8)&0xFF;
  drSide=(cluType>>16);
  return;
}
//----------------------------------------------------------------------
Int_t AliITSRecPoint::GetSSDclusterType() const {
// returns an Int_t with encoded information on cluster size
// Byte1 = fNz Byte0=fNy, other two bytes empty for extra information
// max. allowed cluster size = 255
  Int_t typ=(fNz&0xFF)<<8;
  typ+=fNy&0xFF;
  return typ;
}

//----------------------------------------------------------------------
void AliITSRecPoint::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
    Bool_t mis;
    Int_t lab[4];
    Float_t hit[6];
    lab[3] = 0; // ??
    *is >> lab[0] >> lab[1] >> lab[2];
    SetLabel(lab[0],0); SetLabel(lab[1],1); SetLabel(lab[2],2);
    *is >> hit[0] >> hit[1] >> hit[2] >> hit[3] >> hit[4] >> hit[5];
    SetX(hit[0]);SetY(hit[1]);SetZ(hit[2]);SetSigmaY2(hit[3]);
    SetSigmaZ2(hit[4]);//fSigmaYZ=hit[5];
    *is >> lab[0] >> mis;
    SetVolumeId(lab[0]);// fIsMisalinged = mis;
    *is >> fXloc >> fZloc >> fdEdX;
    *is >> fIndex >> fQ >> fLayer >> fNz >> fNy >> fChargeRatio >> fType;
    *is >> fDeltaProb >> fDriftTime >> fDriftSide;

    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSRecPoint &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////
 
    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSRecPoint &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////
 
    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
