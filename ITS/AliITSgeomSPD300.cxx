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
Revision 1.1.2.8  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.6  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.4  2000/06/10 20:34:37  nilsen
Fixed compilation warning with HP unix.

Revision 1.3  2000/06/10 10:43:04  nilsen
Fixed bug in copy and operator =.


*/

#include "TBRIK.h"
#include "AliITSgeomSPD300.h"

ClassImp(AliITSgeomSPD300)

AliITSgeomSPD300::AliITSgeomSPD300(){
////////////////////////////////////////////////////////////////////////
//    default constructor, for ITS TDR geometry.
////////////////////////////////////////////////////////////////////////
const Float_t kdx=0.6400,kdy=0.0075,kdz=4.1890; // cm; Standard pixel detector
                                                // size is 2dx wide, 2dz long,
                                                // and 2dy thick. Geant 3.12
                                                // units.
const Float_t kbinx0 = 0.0050; // cm; Standard pixel size in x direction.
const Int_t   knbinx = 256;    // number of pixels along x direction.
const Float_t kbinz0 = 0.0300; // cm; Standard pixel size in z direction.
const Int_t   knbinz = 279;    // number of pixels along z direction.
    Int_t i;

    fDx = kdx;  // default value.
    fDy = kdy;  // default value.
    fDz = kdz;  // default value.
    fNbinx = knbinx; // default number of bins in x.
    fNbinz = knbinz; // default number of bins in z.

    fBinSizeX = new Float_t[fNbinx]; // array of bin sizes along x.
    for(i=0;i<fNbinx;i++) fBinSizeX[i] = kbinx0; // default x bin size.
    fBinSizeZ = new Float_t[fNbinz]; // array of bin sizes along z.
    for(i=0;i<fNbinz;i++) fBinSizeZ[i] = kbinz0; // default z bin size.

    // correct detector size for bin size.
    fDx = 0.0;
    for(i=0;i<fNbinx;i++) fDx +=fBinSizeX[i];
    fDx *= 0.5;
    fDz = 0.0;
    for(i=0;i<fNbinz;i++) fDz +=fBinSizeZ[i];
    fDz *= 0.5;

    fShapeSPD = new TBRIK("ActiveSPD","Active volume of SPD","SPD SI DET",
			  fDx,fDy,fDz);
}
//______________________________________________________________________
AliITSgeomSPD300::AliITSgeomSPD300(Float_t dy,Int_t nx,Float_t *bx,
				                Int_t nz,Float_t *bz){
////////////////////////////////////////////////////////////////////////
//    default constructor, for a User modified TDR based geometry.
////////////////////////////////////////////////////////////////////////
    Int_t i;
    fDx = 0.0;
    fDy =  dy;
    fDz = 0.0;
    fNbinx = nx; // new number of bins in x.
    fNbinz = nz; // new number of bins in z.

    fBinSizeX = new Float_t[fNbinx]; // array of bin sizes along x.
    for(i=0;i<fNbinx;i++) fBinSizeX[i] = bx[i]; // new x bin size.
    fBinSizeZ = new Float_t[fNbinz]; // array of bin sizes along z.
    for(i=0;i<fNbinz;i++) fBinSizeZ[i] = bz[i]; // new z bin size.

    // correct detector size for bin size.
    for(i=0;i<fNbinx;i++) fDx +=fBinSizeX[i];
    fDx *= 0.5;
    for(i=0;i<fNbinz;i++) fDz +=fBinSizeZ[i];
    fDz *= 0.5;

    fShapeSPD = new TBRIK("ActiveSPD","Active volume of SPD","SPD SI DET",
			  fDx,fDy,fDz);
}
//______________________________________________________________________
AliITSgeomSPD300::AliITSgeomSPD300(AliITSgeomSPD300 &source){
  // copy constructor
    Int_t i;
    if(&source == this) return;
    this->fShapeSPD = new TBRIK(*(source.fShapeSPD));
    this->fDx = source.fDx;
    this->fDy = source.fDy;
    this->fDz = source.fDz;
    if(this->fBinSizeX) delete[] this->fBinSizeX; 
    if(this->fBinSizeX) delete[] this->fBinSizeZ;
    this->fNbinx = source.fNbinx;
    this->fBinSizeX = new Float_t[this->fNbinx];
    this->fNbinz = source.fNbinz;
    this->fBinSizeZ = new Float_t[this->fNbinz];
    for(i=0;i<fNbinx;i++) this->fBinSizeX[i] = source.fBinSizeX[i];
    for(i=0;i<fNbinz;i++) this->fBinSizeZ[i] = source.fBinSizeZ[i];
}
//______________________________________________________________________
AliITSgeomSPD300& AliITSgeomSPD300::operator=(AliITSgeomSPD300 &source){
  // = operator
    Int_t i;
    if(&source == this) return *this;
    this->fShapeSPD = new TBRIK(*(source.fShapeSPD));
    this->fDx = source.fDx;
    this->fDy = source.fDy;
    this->fDz = source.fDz;
    if(this->fBinSizeX) delete[] this->fBinSizeX; 
    if(this->fBinSizeX) delete[] this->fBinSizeZ;
    this->fNbinx = source.fNbinx;
    this->fBinSizeX = new Float_t[this->fNbinx];
    this->fNbinz = source.fNbinz;
    this->fBinSizeZ = new Float_t[this->fNbinz];
    for(i=0;i<fNbinx;i++) this->fBinSizeX[i] = source.fBinSizeX[i];
    for(i=0;i<fNbinz;i++) this->fBinSizeZ[i] = source.fBinSizeZ[i];
    return *this;
}
//______________________________________________________________________
AliITSgeomSPD300::~AliITSgeomSPD300(){
  // destructor
    delete[] fBinSizeX;
    delete[] fBinSizeZ;
    delete   fShapeSPD;
}
//______________________________________________________________________
void AliITSgeomSPD300::ReSetBins(Float_t dy,Int_t nx,Float_t *bx,
                                        Int_t nz,Float_t *bz){
////////////////////////////////////////////////////////////////////////
//    default constructor, for a User modified TDR based geometry.
////////////////////////////////////////////////////////////////////////
    Int_t i;
    fDx = 0.0;
    fDy =  dy;
    fDz = 0.0;
    fNbinx = nx; // new number of bins in x.
    fNbinz = nz; // new number of bins in z.

    if(fBinSizeX!=0) delete[] fBinSizeX;
    fBinSizeX = new Float_t[fNbinx]; // array of bin sizes along x.
    for(i=0;i<fNbinx;i++) fBinSizeX[i] = bx[i]; // new x bin size.
    if(fBinSizeZ!=0) delete[] fBinSizeZ;
    fBinSizeZ = new Float_t[fNbinz]; // array of bin sizes along z.
    for(i=0;i<fNbinz;i++) fBinSizeZ[i] = bz[i]; // new z bin size.

    // correct detector size for bin size.
    for(i=0;i<fNbinx;i++) fDx +=fBinSizeX[i];
    fDx *= 0.5;
    for(i=0;i<fNbinz;i++) fDz +=fBinSizeZ[i];
    fDz *= 0.5;

    fShapeSPD = new TBRIK("ActiveSPD","Active volume of SPD","SPD SI DET",
			  fDx,fDy,fDz);
}
/*
//----------------------------------------------------------------------
void AliITSgeomSPD300::Streamer(TBuffer &R__b){
    // Streamer function for the class AliItSgeomSPD300.
    Int_t i;
    UInt_t R__s, R__c;

    if(R__b.IsReading()){
	Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
	if(R__v==1){
	    TObject::Streamer(R__b);
	    fShapeSPD->Streamer(R__b);
	    R__b >> fDx;
	    R__b >> fDy;
	    R__b >> fDz;
	}else if (R__v==2){
	    AliITSgeomSPD::Streamer(R__b);
	    fShapeSPD->Streamer(R__b);
	    R__b >> fDx;
	    R__b >> fDy;
	    R__b >> fDz;
	    R__b >> fNbinx;
	    if(fBinSizeX!=0) delete[] fBinSizeX;
	    fBinSizeX = new Float_t[fNbinx];
	    for(i=0;i<fNbinx;i++) R__b >> fBinSizeX[i];
	    R__b >> fNbinz;
	    if(fBinSizeZ!=0) delete[] fBinSizeZ;
	    fBinSizeZ = new Float_t[fNbinz];
	    for(i=0;i<fNbinz;i++) R__b >> fBinSizeZ[i];
	    R__b.CheckByteCount(R__s, R__c, AliITSgeomSPD300::IsA());
	} // end if R__v==1
    } else { // IsWriting.
	R__c = R__b.WriteVersion(AliITSgeomSPD300::IsA(), kTRUE);
	AliITSgeomSPD::Streamer(R__b);
	fShapeSPD->Streamer(R__b);
	R__b << fDx;
	R__b << fDy;
	R__b << fDz;
	R__b << fNbinx;
	for(i=0;i<fNbinx;i++) R__b << fBinSizeX[i];
	R__b << fNbinz;
	for(i=0;i<fNbinz;i++) R__b << fBinSizeZ[i];
	R__b.SetByteCount(R__c, kTRUE);
    } // end if R__b.IsReading()
}
//----------------------------------------------------------------------
*/
