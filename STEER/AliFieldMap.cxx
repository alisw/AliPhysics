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

//-----------------------------------------------------------------------
//
// Class to handle the field
// I/O and interpolation
// of the field map to return the value
// of the magnetic field in an arbitrary position
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//
//-----------------------------------------------------------------------

#include <TClass.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliFieldMap.h"

ClassImp(AliFieldMap)

//_______________________________________________________________________
AliFieldMap::AliFieldMap():
  fXbeg(0),
  fYbeg(0),
  fZbeg(0),
  fXend(0),
  fYend(0),
  fZend(0),
  fXdel(0),
  fYdel(0),
  fZdel(0),
  fXdeli(0),
  fYdeli(0),
  fZdeli(0),
  fXn(0),
  fYn(0),
  fZn(0),
  fWriteEnable(0),
  fB(0)
{
  //
  // Standard constructor
  //
  SetWriteEnable();
}

//_______________________________________________________________________
AliFieldMap::AliFieldMap(const char *name, const char *title):
  TNamed(name,title),
  fXbeg(0),
  fYbeg(0),
  fZbeg(0),
  fXend(0),
  fYend(0),
  fZend(0),
  fXdel(0),
  fYdel(0),
  fZdel(0),
  fXdeli(0),
  fYdeli(0),
  fZdeli(0),
  fXn(0),
  fYn(0),
  fZn(0),
  fWriteEnable(0),
  fB(0)
{
  //
  // Standard constructor
  //
  ReadField();
  SetWriteEnable();
}

//_______________________________________________________________________
AliFieldMap::~AliFieldMap()
{
  //
  // Destructor
  //  
  delete fB;
}

//_______________________________________________________________________
AliFieldMap::AliFieldMap(const AliFieldMap &map):
  TNamed(map),
  fXbeg(0),
  fYbeg(0),
  fZbeg(0),
  fXend(0),
  fYend(0),
  fZend(0),
  fXdel(0),
  fYdel(0),
  fZdel(0),
  fXdeli(0),
  fYdeli(0),
  fZdeli(0),
  fXn(0),
  fYn(0),
  fZn(0),
  fWriteEnable(0),
  fB(0)
{
  //
  // Copy constructor
  //
  map.Copy(*this);
}

//_______________________________________________________________________
void AliFieldMap::ReadField()
{
  // 
  // Method to read the magnetic field map from file
  //
  FILE* magfile;
  //  FILE* endf = fopen("end.table", "r");
  //  FILE* out  = fopen("out", "w");
  
  Int_t   ix, iy, iz, ipx, ipy, ipz;
  Float_t bx, by, bz;
  char *fname = 0;
  AliInfo(Form("Reading Magnetic Field Map %s from file %s",
	       fName.Data(),fTitle.Data()));

  fname   = gSystem->ExpandPathName(fTitle.Data());
  magfile = fopen(fname,"r");
  delete [] fname;

  fscanf(magfile,"%d %d %d %f %f %f %f %f %f", 
	 &fXn, &fYn, &fZn, &fXbeg, &fYbeg, &fZbeg, &fXdel, &fYdel, &fZdel);
 
  fXdeli = 1./fXdel;
  fYdeli = 1./fYdel;	
  fZdeli = 1./fZdel;	
  fXend  = fXbeg + (fXn-1)*fXdel;
  fYend  = fYbeg + (fYn-1)*fYdel;
  fZend  = fZbeg + (fZn-1)*fZdel;
  
  Int_t nDim   = fXn*fYn*fZn;

  //  Float_t x,y,z,b;

  fB = new TVector(3*nDim);
  if (magfile) {
      for (ix = 0; ix < fXn; ix++) {
	  ipx=ix*3*(fZn*fYn);
	  for (iy = 0; iy < fYn; iy++) {
	      ipy=ipx+iy*3*fZn;
	      for (iz = 0; iz < fZn; iz++) {
		  ipz=ipy+iz*3;

		  if (iz == -1) {
//		      fscanf(endf,"%f %f %f", &bx,&by,&bz);
		  } else if (iz > -1) {
		      fscanf(magfile," %f %f %f", &bx, &by, &bz);
		  } else {
		      continue;
		  }
		  
//		  fscanf(magfile,"%f %f %f %f %f %f %f ",
//			 &x, &y, &z, &bx,&by,&bz, &b);
//		  fprintf(out, "%15.8e %15.8e %15.8e \n", bx, by, bz);

		  (*fB)(ipz+2) = 10.*bz;
		  (*fB)(ipz+1) = 10.*by;
		  (*fB)(ipz  ) = 10.*bx;
	      } //iz
	  } // iy
      } // ix
/*
//
// this part for interpolation between z = 700 and 720 cm to get
// z = 710 cm
//      
      for (ix = 0; ix < fXn; ix++) {
	  ipx=ix*3*(fZn*fYn);
	  for (iy = 0; iy < fYn; iy++) {
	      ipy=ipx+iy*3*fZn;
	      Float_t bxx = (Bx(ix,iy,0) + Bx(ix,iy,2))/2.;
	      Float_t byy = (By(ix,iy,0) + By(ix,iy,2))/2.;
	      Float_t bzz = (Bz(ix,iy,0) + Bz(ix,iy,2))/2.;	      
	      ipz=ipy+3;
	      (*fB)(ipz+2) = bzz;
	      (*fB)(ipz+1) = byy;
	      (*fB)(ipz  ) = bxx;
	  } // iy
      } // ix
*/      
  } else { 
      printf("%s: File %s not found !\n",ClassName(),fTitle.Data());
      exit(1);
  } // if mafile
}

//_______________________________________________________________________
void AliFieldMap::Field(Float_t *x, Float_t *b) const
{
  //
  // Use simple interpolation to obtain field at point x
  //
    Double_t ratx, raty, ratz, hix, hiy, hiz, ratx1, raty1, ratz1, 
	bhyhz, bhylz, blyhz, blylz, bhz, blz, xl[3];
    const Double_t kone=1;
    Int_t ix, iy, iz;
    b[0]=b[1]=b[2]=0;
    //
    
    xl[0] = x[0] - fXbeg;
    xl[1] = x[1] - fYbeg;
    xl[2] = x[2] - fZbeg;
    
    hix=TMath::Max(0.,TMath::Min(xl[0]*fXdeli,fXn-1.0001));
    ratx=hix-int(hix);
    ix=int(hix);
    
    hiy=TMath::Max(0.,TMath::Min(xl[1]*fYdeli,fYn-1.0001));
    raty=hiy-int(hiy);
    iy=int(hiy);
    
    hiz=TMath::Max(0.,TMath::Min(xl[2]*fZdeli,fZn-1.0001));
    ratz=hiz-int(hiz);
    iz=int(hiz);

    ratx1=kone-ratx;
    raty1=kone-raty;
    ratz1=kone-ratz;

    if (!fB) return;

    bhyhz = Bx(ix  ,iy+1,iz+1)*ratx1+Bx(ix+1,iy+1,iz+1)*ratx;
    bhylz = Bx(ix  ,iy+1,iz  )*ratx1+Bx(ix+1,iy+1,iz  )*ratx;
    blyhz = Bx(ix  ,iy  ,iz+1)*ratx1+Bx(ix+1,iy  ,iz+1)*ratx;
    blylz = Bx(ix  ,iy  ,iz  )*ratx1+Bx(ix+1,iy  ,iz  )*ratx;
    bhz   = blyhz             *raty1+bhyhz             *raty;
    blz   = blylz             *raty1+bhylz             *raty;
    b[0]  = blz               *ratz1+bhz               *ratz;
    //
    bhyhz = By(ix  ,iy+1,iz+1)*ratx1+By(ix+1,iy+1,iz+1)*ratx;
    bhylz = By(ix  ,iy+1,iz  )*ratx1+By(ix+1,iy+1,iz  )*ratx;
    blyhz = By(ix  ,iy  ,iz+1)*ratx1+By(ix+1,iy  ,iz+1)*ratx;
    blylz = By(ix  ,iy  ,iz  )*ratx1+By(ix+1,iy  ,iz  )*ratx;
    bhz   = blyhz             *raty1+bhyhz             *raty;
    blz   = blylz             *raty1+bhylz             *raty;
    b[1]  = blz               *ratz1+bhz               *ratz;
    //
    bhyhz = Bz(ix  ,iy+1,iz+1)*ratx1+Bz(ix+1,iy+1,iz+1)*ratx;
    bhylz = Bz(ix  ,iy+1,iz  )*ratx1+Bz(ix+1,iy+1,iz  )*ratx;
    blyhz = Bz(ix  ,iy  ,iz+1)*ratx1+Bz(ix+1,iy  ,iz+1)*ratx;
    blylz = Bz(ix  ,iy  ,iz  )*ratx1+Bz(ix+1,iy  ,iz  )*ratx;
    bhz   = blyhz             *raty1+bhyhz             *raty;
    blz   = blylz             *raty1+bhylz             *raty;
    b[2]  = blz               *ratz1+bhz               *ratz;
}

//_______________________________________________________________________
void AliFieldMap::Copy(TObject & /* magf */) const
{
  //
  // Copy *this onto magf -- Not implemented
  //
  AliFatal("Not implemented!");
}

//_______________________________________________________________________
AliFieldMap & AliFieldMap::operator =(const AliFieldMap &magf)
{
  magf.Copy(*this);
  return *this;
}

//_______________________________________________________________________
void AliFieldMap::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliFieldMap.
    TVector* save = 0;
    
    if (R__b.IsReading()) {
	AliFieldMap::Class()->ReadBuffer(R__b, this);
    } else {
	if (!fWriteEnable) {
	    save = fB;
	    fB = 0;
	}
	AliFieldMap::Class()->WriteBuffer(R__b, this);
	if (!fWriteEnable) fB = save;
    }
}
