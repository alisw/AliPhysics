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

/*
$Log$
Revision 1.3  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.2  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.2.1  2000/09/22 14:47:52  cblume
Add the tracking code

*/                        
                                
#include <TObject.h>

#include "AliRun.h"

#include "AliTRD.h" 
#include "AliTRDgeometry.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtimeBin.h" 
#include "AliTRDtrackingSector.h" 

ClassImp(AliTRDtrackingSector) 

//_______________________________________________________

AliTRDtrackingSector::~AliTRDtrackingSector()
{
  //
  // Destructor
  //

  delete[] fTimeBin;

}

//_______________________________________________________
AliTRDtimeBin &AliTRDtrackingSector::operator[](Int_t i)
{
  //
  // Index operator 
  //

  return *(fTimeBin+i);

}

//_______________________________________________________

void AliTRDtrackingSector::SetUp()
{ 

  AliTRD *TRD = (AliTRD*) gAlice->GetDetector("TRD");
  fGeom = TRD->GetGeometry();

  fTimeBinSize = fGeom->GetTimeBinSize();

  fN = AliTRDgeometry::Nplan() * (Int_t(AliTRDgeometry::DrThick()
                                       /fTimeBinSize) + 1);

  fTimeBin = new AliTRDtimeBin[fN]; 

}

//______________________________________________________

Double_t AliTRDtrackingSector::GetX(Int_t l) const
{
  if( (l<0) || (l>fN-1)) {
    fprintf(stderr,"AliTRDtrackingSector::GetX: TimeBin index is out of range !\n");
    return -99999.;
  }
  else { 
    
    Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();
    Int_t plane = l/tb_per_plane;
    Int_t time_slice = l%(Int_t(AliTRDgeometry::DrThick()
                               /fTimeBinSize) + 1);

    Float_t t0 = fGeom->GetTime0(plane);
    Double_t x = t0 + time_slice * fTimeBinSize;

    //    cerr<<"plane, tb, x = "<<plane<<","<<time_slice<<","<<x<<endl;

    return x;
  }
}

//______________________________________________________

Double_t AliTRDtrackingSector::GetMaxY(Int_t l) const 
{ 

  if((l<(fN-1)) && (l>-1)) { 
    Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();
    Int_t plane = l/tb_per_plane;
    return fGeom->GetChamberWidth(plane); 
  }
  else {
    fprintf(stderr,
    "AliTRDtrackingSector::GetMaxY: TimeBin index is out of range !\n");
    if(l<0) return fGeom->GetChamberWidth(0);
    else return fGeom->GetChamberWidth(AliTRDgeometry::Nplan()-1);
  }
}

//______________________________________________________

Int_t AliTRDtrackingSector::GetTimeBinNumber(Double_t x) const
{
  //Float_t r_out = fGeom->GetTime0(AliTRDgeometry::Nplan()-1) 
  //              + AliTRDgeometry::DrThick(); 
  // Changed to new time0 (CBL)
  Float_t r_out = fGeom->GetTime0(AliTRDgeometry::Nplan()-1); 
  Float_t r_in = fGeom->GetTime0(0);
  //  cerr<<"GetTimeBinNumber: r_in,r_out = "<<r_in<<","<<r_out<<endl;

  if(x >= r_out) return fN-1;
  if(x <= r_in) return 0;

  Float_t gap = fGeom->GetTime0(1) - fGeom->GetTime0(0);
  //  cerr<<"GetTimeBinNumber: gap = "<<gap<<endl;

  Int_t plane = Int_t((x - r_in + fTimeBinSize/2)/gap);
  //  cerr<<"GetTimeBinNumber: plane="<<plane<<endl;

  Int_t local_tb = Int_t((x-fGeom->GetTime0(plane))/fTimeBinSize + 0.5);
  //  cerr<<"GetTimeBinNumber: local_tb="<<local_tb<<endl;

  Int_t time_bin = plane * (Int_t(AliTRDgeometry::DrThick()
                                 /fTimeBinSize) + 1) + local_tb;
  //  cerr<<"GetTimeBinNumber: time_bin = "<<time_bin<<endl;
  

  return time_bin;
}

//______________________________________________________

Int_t AliTRDtrackingSector::GetTimeBin(Int_t det, Int_t local_tb) const 
{
  Int_t plane = fGeom->GetPlane(det);

  Int_t time_bin = plane * (Int_t(AliTRDgeometry::DrThick()
                                 /fTimeBinSize) + 1) + local_tb;  
  return time_bin;
}

