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
Revision 1.5  2000/12/08 16:07:02  cblume
Update of the tracking by Sergei

Revision 1.4  2000/10/16 01:16:53  cblume
Changed timebin 0 to be the one closest to the readout

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

Double_t AliTRDtrackingSector::GetX(Int_t tb) const
{
  if( (tb<0) || (tb>fN-1)) {
    fprintf(stderr,"AliTRDtrackingSector::GetX: TimeBin index is out of range !\n");
    return -99999.;
  }
  else { 
    
    Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();
    Int_t local_tb = tb_per_plane - tb%tb_per_plane - 1;

    Int_t plane = tb/tb_per_plane;
    Float_t t0 = fGeom->GetTime0(plane);
    Double_t x = t0 - (local_tb + 0.5) * fTimeBinSize;

    return x;
  }
}

//______________________________________________________

Int_t AliTRDtrackingSector::GetTimeBinNumber(Double_t x) const
{
  Float_t r_out = fGeom->GetTime0(AliTRDgeometry::Nplan()-1); 
  Float_t r_in = fGeom->GetTime0(0) - AliTRDgeometry::DrThick();


  if(x >= r_out) return fN-1;
  if(x <= r_in) return 0;

  Int_t plane;
  for (plane = AliTRDgeometry::Nplan()-1; plane >= 0; plane--) {
    if(x > (fGeom->GetTime0(plane) - AliTRDgeometry::DrThick())) break;
  }  
 
  Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();
  Int_t local_tb = Int_t((fGeom->GetTime0(plane)-x)/fTimeBinSize);

  if((local_tb < 0) || (local_tb >= tb_per_plane)) {
    printf("AliTRDtrackingSector::GetTimeBinNumber: \n");
    printf("local time bin %d is out of bounds [0..%d]: x = %f \n",
	   local_tb,tb_per_plane-1,x);
    return -1;
  }
      
  Int_t time_bin = (plane + 1) * tb_per_plane - 1 - local_tb;

  return time_bin;
}

//______________________________________________________

Int_t AliTRDtrackingSector::GetTimeBin(Int_t det, Int_t local_tb) const 
{
  Int_t plane = fGeom->GetPlane(det);

  Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();

  Int_t time_bin = (plane + 1) * tb_per_plane - 1 - local_tb;

  return time_bin;
}


//______________________________________________________

Bool_t AliTRDtrackingSector::TECframe(Int_t tb, Double_t y, Double_t z) const
{
// 
// Returns <true> if point defined by <x(tb),y,z> is within 
// the TEC G10 frame, otherwise returns <false>  
//  

  if((tb > (fN-1)) || (tb < 0)) return kFALSE; 

  Int_t tb_per_plane = fN/AliTRDgeometry::Nplan();
  Int_t plane = tb/tb_per_plane;
  
  Double_t x = GetX(tb);
  y = TMath::Abs(y);

  if((y > fGeom->GetChamberWidth(plane)/2.) &&
     (y < x*TMath::Tan(0.5*AliTRDgeometry::GetAlpha()))) return kTRUE; 

  Double_t zmin, zmax;
  Float_t  fRowPadSize, fRow0;
  Int_t    nPadRows;

  for(Int_t iCha = 1; iCha < AliTRDgeometry::Ncham(); iCha++) {

    fRow0 = fGeom->GetRow0(plane,iCha-1,0);
    fRowPadSize = fGeom->GetRowPadSize(plane,iCha-1,0);
    nPadRows = fGeom->GetRowMax(plane,iCha-1,0);
    zmin = fRow0 - fRowPadSize/2 + fRowPadSize * nPadRows;

    fRow0 = fGeom->GetRow0(plane,iCha,0);
    fRowPadSize = fGeom->GetRowPadSize(plane,iCha,0);
    zmax = fRow0 - fRowPadSize/2;

    if((z > zmin) && (z < zmax)) return kTRUE;     
  }

  return kFALSE;
}
