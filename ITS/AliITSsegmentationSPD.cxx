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

#include <TMath.h>

#include "AliITSsegmentationSPD.h"

ClassImp(AliITSsegmentationSPD)


Float_t ColFromZ300(Float_t z) {
// Get column number for each z-coordinate taking into account the 
// extra pixels in z direction assuming 300 micron sized pixels.
     Float_t col = 0.0;
     Float_t pitchz = 300.0;
     col = Float_t (z/pitchz);
     return col;
}
//_____________________________________________________________________________
Float_t ZFromCol300(Int_t col) {
// same comments as above
// Get z-coordinate for each colunm number
  Float_t pitchz = 300.0;
  Float_t z = 0.0;
  z = (col+0.5)*pitchz;
  return z;
}
//_____________________________________________________________________________
Float_t ZpitchFromCol300(Int_t col) {
  // returns Z pixel pitch for 300 micron pixels.
  return 300.0;
}
//_____________________________________________________________________________
Float_t ColFromZ(Float_t z) {
// hard-wired - keep it like this till we can parametrise 
// and get rid of AliITSgeomSPD425
// Get column number for each z-coordinate taking into account the 
// extra pixels in z direction 

  Float_t col = 0;
  Float_t pitchz = 425;
  if( z < 13175) {
    col = Float_t(z/pitchz);
  } else if( z < 14425) {  
    pitchz = 625;
    col = 31 + (z - 13175)/pitchz;
  } else if( z < 27175) {  
    col = 33 + (z - 14425)/pitchz;
  } else if( z < 28425) {  
    pitchz = 625;
    col = 63 + (z - 27175)/pitchz;
  } else if( z < 41175) {  
    col = 65 + (z - 28425)/pitchz;
  } else if( z < 42425) {  
    pitchz = 625;
    col = 95 + (z - 41175)/pitchz;
  } else if( z < 55175) {  
    col = 97 + (z - 42425)/pitchz;
  } else if( z < 56425) {  
    pitchz = 625;
    col = 127 + (z - 55175)/pitchz;
  } else if( z < 69175) {  
    col = 129 + (z - 56425)/pitchz;
  } else if( z < 70425) {  
    pitchz = 625;
    col = 159 + (z - 69175)/pitchz;
  } else if( z < 83600) {  
    col = 161 + (z - 70425)/pitchz;
  }   
  return col;
}

//_____________________________________________________________________________
Float_t ZFromCol(Int_t col) {
// same comments as above
// Get z-coordinate for each colunm number

  Float_t pitchz = 425;
  Float_t z = 0;
  if( col >=0 && col <= 30 ) {  
    z = (col + 0.5)*pitchz;    
  } else if( col >= 31 && col <= 32) {  
    pitchz = 625;
    z = 13175 + (col -31 + 0.5)*pitchz;    
  } else if( col >= 33 && col <= 62) {  
    z = 14425 + (col -33 + 0.5)*pitchz;    
  } else if( col >= 63 && col <= 64) {  
    pitchz = 625;
    z = 27175 + (col -63 + 0.5)*pitchz;    
  } else if( col >= 65 && col <= 94) {  
    z = 28425 + (col -65 + 0.5)*pitchz;    
  } else if( col >= 95 && col <= 96) {  
    pitchz = 625;
    z = 41175 + (col -95 + 0.5)*pitchz;    
  } else if( col >= 97 && col <= 126) {  
    z = 42425 + (col -97 + 0.5)*pitchz;    
  } else if( col >= 127 && col <= 128) {  
    pitchz = 625;
    z = 55175 + (col -127 + 0.5)*pitchz;    
  } else if( col >= 129 && col <= 158) {  
    z = 56425 + (col -129 + 0.5)*pitchz;    
  } else if( col >= 159 && col <= 160) {  
    pitchz = 625;
    z = 69175 + (col -159 + 0.5)*pitchz;    
  } else if( col >= 161 && col <= 191) {  
    z = 70425 + (col -161 + 0.5)*pitchz;    
  }   

  return z;
}

Float_t ZpitchFromCol(Int_t col) {
// Get pitch size in z direction for each colunm

  Float_t pitchz = 425;
  if( col >=32 && col <= 33 ) {  
    pitchz = 625;
  } else if( col >= 64 && col <= 65) {  
    pitchz = 625;
  } else if( col >= 96 && col <= 97) {  
    pitchz = 625;
  } else if( col >= 128 && col <= 129) {  
    pitchz = 625;
  } else if( col >= 160 && col <= 161) {  
    pitchz = 625;
  }   
  return pitchz;
}

AliITSsegmentationSPD::AliITSsegmentationSPD(){
  // Default constructor
   fNpx = 0;
   fNpz = 0;
   fCorr=0;
   fGeom = 0;

}
//____________________________________________________________________________
AliITSsegmentationSPD::AliITSsegmentationSPD(AliITSgeom *gm){
  // Constructor
   fCorr=0;
   fNpx = 0;
   fNpz = 0;
   Init(); 
   fGeom = gm;

}
//____________________________________________________________________________
AliITSsegmentationSPD& AliITSsegmentationSPD::operator=(AliITSsegmentationSPD &source){
   // = operator
   Int_t i;
   if(this==&source) return *this;
   this->fNpx  = source.fNpx;
   this->fNpz  = source.fNpz;
   this->fDx   = source.fDx;
   this->fDy   = source.fDy;
   for(i=0;i<256;i++) this->fCellSizeX[i] = source.fCellSizeX[i];
   for(i=0;i<280;i++) this->fCellSizeZ[i] = source.fCellSizeZ[i];
   this->fCorr = new TF1(*(source.fCorr));// make a proper copy of the function
   this->fGeom = source.fGeom; // copy only the pointers.
   return *this;
}
//____________________________________________________________________________
AliITSsegmentationSPD::AliITSsegmentationSPD(AliITSsegmentationSPD &source){
  // copy constructor
   *this = source;
}
//------------------------------
void AliITSsegmentationSPD::Init300(){
// Initialize infromation for 6 read out chip 300X50 micron pixel SPD 
// detectors. This chip is 150 microns thick by 1.28 cm in x by 8.37 cm
// long. It has 256  50 micron pixels in x and 279 300 micron size
// pixels in z.

    Int_t i;
    //const Float_t kconv=10000.;
    fNpx = 256; // The number of X pixel Cell same as in fCellSizeX array size
    fNpz = 279; // The number of Z pixel Cell same as in fCellSizeZ array size
    for(i=0;i<fNpx;i++) fCellSizeX[i] = 50.0; // microns all the same
    for(i=0;i<fNpz;i++) fCellSizeZ[i] = ZpitchFromCol300(i); // microns
    for(i=fNpz;i<280;i++) fCellSizeZ[i] = 0.0; // zero out rest of array
    fDx = 0;
    for(i=0;i<fNpx;i++) fDx += fCellSizeX[i];
    fDz = 0;
    for(i=0;i<fNpz;i++) fDz += fCellSizeZ[i];
    fDy = 300.0; //microns  SPD sensitive layer thickness
}

//------------------------------
void AliITSsegmentationSPD::Init(){
// Initialize infromation for 6 read out chip 425X50 micron pixel SPD 
// detectors. This chip is 150 microns thick by 1.28 cm in x by 8.375 cm
// long. It has 256  50 micron pixels in x and 197 mostly 425 micron size
// pixels in z. The two pixels between each readout chip are 625 microns long.

    Int_t i;
    //const Float_t kconv=10000.;
    fNpx = 256; // The number of X pixel Cell same as in fCellSizeX array size
    fNpz = 192; // The number of Z pixel Cell same as in fCellSizeZ array size
    for(i=0;i<fNpx;i++) fCellSizeX[i] = 50.0; // microns all the same
    for(i=0;i<fNpz;i++) fCellSizeZ[i] = ZpitchFromCol(i); // microns
    for(i=fNpz;i<280;i++) fCellSizeZ[i] = 0.0; // zero out rest of array
    fDx = 0;
    for(i=0;i<fNpx;i++) fDx += fCellSizeX[i];
    fDz = 0;
    for(i=0;i<fNpz;i++) fDz += fCellSizeZ[i];
    fDy = 300.0; //microns  SPD sensitive layer thickness
    printf(" AliITSsegmentationSPD - Init: fNpx fNpz fDx fDz %d %d %f %f\n",fNpx, fNpz, fDx, fDz);

}
//------------------------------
void AliITSsegmentationSPD::SetNCells(Int_t p1, Int_t p2){
  // for SPD this function should be used ONLY when a beam test setup 
  // configuration is studied

    fNpx=p1;
    fNpz=p2;

}
//------------------------------
void AliITSsegmentationSPD::SetDetSize(Float_t p1, Float_t p2, Float_t p3){
  // for SPD this function should be used ONLY when a beam test setup 
  // configuration is studied

    fDx=p1;
    fDz=p2;
    fDy=p3;

}
//------------------------------
Float_t AliITSsegmentationSPD::Dpx(Int_t i){
   //returs x pixel pitch for a give pixel
   return fCellSizeX[i];
}
//------------------------------
Float_t AliITSsegmentationSPD::Dpz(Int_t i){
  // returns z pixel pitch for a give pixel
   return ZpitchFromCol(i);
}
//------------------------------
void AliITSsegmentationSPD::GetCellIxz(Float_t &x,Float_t &z,Int_t &ix,Int_t &iz){
//  Returns pixel coordinates (ix,iz) for given real local coordinates (x,z)
//

    // expects x, z in microns

    // same segmentation on x
    Float_t dpx=Dpx(0);
    ix = (Int_t)(x/dpx + 1);
    // different segmentation on z
    iz = (Int_t)(ColFromZ(z) + 1);

    x /= dpx;
    z = ColFromZ(z);

    if (iz >  fNpz) iz= fNpz;
    if (ix >  fNpx) ix= fNpx;

    /*
    if (iz < -fNpz) iz= -fNpz;
    if (ix < -fNpx) ix=-fNpx;
    */
}

//------------------------------
void AliITSsegmentationSPD::GetCellCxz(Int_t ix,Int_t iz,Float_t &x,Float_t&z){
    // Transform from pixel to real local coordinates

    // returns x, z in microns

    Float_t dpx=Dpx(0);

    x = (ix>0) ? Float_t(ix*dpx)-dpx/2. : Float_t(ix*dpx)+dpx/2.;
    z = ZFromCol(iz);


}
//------------------------------
void AliITSsegmentationSPD::
Neighbours(Int_t iX, Int_t iZ, Int_t* Nlist, Int_t Xlist[8], Int_t Zlist[8]){
  // returns the neighbouring pixels for use in Cluster Finders and the like.
  /*
    *Nlist=4;Xlist[0]=Xlist[1]=iX;Xlist[2]=iX-1;Xlist[3]=iX+1;
    Zlist[0]=iZ-1;Zlist[1]=iZ+1;Zlist[2]=Zlist[3]=iZ;
  */


    *Nlist=8;
    Xlist[0]=Xlist[1]=iX;
    Xlist[2]=iX-1;
    Xlist[3]=iX+1;
    Zlist[0]=iZ-1;
    Zlist[1]=iZ+1;
    Zlist[2]=Zlist[3]=iZ;

   // Diagonal elements
    Xlist[4]=iX+1;
    Zlist[4]=iZ+1;

    Xlist[5]=iX-1;
    Zlist[5]=iZ-1;

    Xlist[6]=iX-1;
    Zlist[6]=iZ+1;

    Xlist[7]=iX+1;
    Zlist[7]=iZ-1;
}
