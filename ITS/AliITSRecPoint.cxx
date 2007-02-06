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


#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
ClassImp(AliITSRecPoint)

//_____________________________________________________________
AliITSRecPoint::AliITSRecPoint(): AliCluster(),
fX(0),
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
fGeom(0){
    // default creator
}

//_____________________________________________________________
AliITSRecPoint::AliITSRecPoint(AliITSgeom* geom): AliCluster(),
fX(0),
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
fGeom(geom) {
    // default creator

}

//________________________________________________________________________
AliITSRecPoint::AliITSRecPoint(Int_t module,AliITSgeom* geom,Int_t *lab,Float_t *hit, Int_t *info):AliCluster(lab,hit),
fX(0),
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
fGeom(geom)
{
  //standard constructor used in AliITSClusterFinderV2


  fType=0;
  fDeltaProb=0.;
  
  fGeom = geom;
  fGeom->TrackingV2ToDetL(module,fY,fZ,fXloc,fZloc);
  if(module<fGeom->GetStartSDD()) fdEdX=0.;
  if(module>=fGeom->GetStartSDD() && module<fGeom->GetStartSSD()){
    fdEdX=fQ*1e-6;
  }
  if(module>=fGeom->GetStartSSD()) fdEdX=fQ*2.16;
  
  
}
//_______________________________________________________________________
AliITSRecPoint::AliITSRecPoint(const AliITSRecPoint& pt):AliCluster(pt),
fX(pt.fX),
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
fGeom(pt.fGeom){
  //Copy constructor

}

//______________________________________________________________________
AliITSRecPoint& AliITSRecPoint::operator=(const AliITSRecPoint& source){
  // Assignment operator

  this->~AliITSRecPoint();
  new(this) AliITSRecPoint(source);
  return *this;

}

//________________________________________________________________________
AliITSRecPoint::AliITSRecPoint(Int_t *lab,Float_t *hit, Int_t *info):AliCluster(lab,hit),
fX(0),
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
fGeom(0){
  //standard constructor used in AliITSClusterFinderV2
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
    *os << fTracks[0]<< " " << fTracks[1] << " " << fTracks[2] << " ";
    *os << fXloc << " " << fZloc << " " << fQ << " ";
    fmt = os->setf(ios::scientific); // set scientific for dEdX.
    *os << fdEdX << " ";
    fmt = os->setf(ios::fixed); // every fixed
    *os << fSigmaY2 << " " << fSigmaZ2;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITSRecPoint::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
 

    *is >> fTracks[0] >> fTracks[1] >> fTracks[2] >> fXloc >> fZloc >> fQ;
    *is >> fdEdX >> fSigmaY2 >> fSigmaZ2;
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
