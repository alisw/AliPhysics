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

////////////////////////////////////////////////////////////////////////
// This class is for the Silicon Pixel Detector, SPD, specific geometry.
// It is being replaced by AliITSsegmentationSPD class. This file also
// constains classes derived from AliITSgeomSPD which do nothing but
// initilize this one with predefined values.
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TMath.h>

#include "AliITSgeomSPD.h"

ClassImp(AliITSgeomSPD)

AliITSgeomSPD::AliITSgeomSPD():
TObject(),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0),
fNbinx(0),
fNbinz(0),
fLowBinEdgeX(0),
fLowBinEdgeZ(0){
// Default Constructor. Set everthing to null.
}
//______________________________________________________________________
AliITSgeomSPD::AliITSgeomSPD(Float_t dy,Int_t nx,Float_t *bx,
			     Int_t nz,Float_t *bz):
TObject(),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0),
fNbinx(0),
fNbinz(0),
fLowBinEdgeX(0),
fLowBinEdgeZ(0){
// Standard Constructor. Set everthing to null.

    ReSetBins(dy,nx,bx,nz,bz);
    return;
}
//______________________________________________________________________
void AliITSgeomSPD::ReSetBins(Float_t dy,Int_t nx,Float_t *bx,
			      Int_t nz,Float_t *bz){
// delets the contents of this and replaces it with the given values.
    Int_t i;
    Float_t dx = 0.0, dz = 0.0;

    // Compute size in x and z (based on bins).
    for(i=0;i<nx;i++) dx += bx[i];
    for(i=0;i<nz;i++) dz += bz[i];
    dx *= 0.5;
    dz *= 0.5;

    if(fLowBinEdgeX) delete[] fLowBinEdgeX; // delete existing
    if(fLowBinEdgeZ) delete[] fLowBinEdgeZ; // delete existing
    fLowBinEdgeX = 0;
    fLowBinEdgeZ = 0;

    SetNbinX(nx);
    SetNbinZ(nz);
    InitLowBinEdgeX();
    InitLowBinEdgeZ();
    fLowBinEdgeX[0] = -dx;
    fLowBinEdgeZ[0] = -dz;
    for(i=0;i<nx;i++) fLowBinEdgeX[i+1] = fLowBinEdgeX[i] + bx[i];
    for(i=0;i<nz;i++) fLowBinEdgeZ[i+1] = fLowBinEdgeZ[i] + bz[i];
    SetShape("ActiveSPD","Active volume of SPD","",dx,dy,dz);
    return;
}
//______________________________________________________________________
AliITSgeomSPD::AliITSgeomSPD(AliITSgeomSPD &source) : TObject(source),
fName(source.fName),
fTitle(source.fTitle),
fMat(source.fMat),
fDx(source.fDx),
fDy(source.fDy),
fDz(source.fDz),
fNbinx(source.fNbinx),
fNbinz(source.fNbinz),
fLowBinEdgeX(0),
fLowBinEdgeZ(0){
    // Copy constructor
  InitLowBinEdgeX();
  InitLowBinEdgeZ();
  for(Int_t i=0;i<fNbinx;i++) fLowBinEdgeX[i] = source.fLowBinEdgeX[i];
  for(Int_t i=0;i<fNbinz;i++) fLowBinEdgeZ[i] = source.fLowBinEdgeZ[i];

  
}
//______________________________________________________________________
AliITSgeomSPD& AliITSgeomSPD::operator=(AliITSgeomSPD &source){
    // = operator
  this->~AliITSgeomSPD();
  new(this) AliITSgeomSPD(source);
  return *this;
}
//______________________________________________________________________
AliITSgeomSPD::~AliITSgeomSPD(){
// Destructor

    if(fLowBinEdgeX) delete[] fLowBinEdgeX;
    if(fLowBinEdgeZ) delete[] fLowBinEdgeZ;
    fNbinx       = 0;
    fNbinz       = 0;
    fLowBinEdgeX = 0;
    fLowBinEdgeZ = 0;
}
//______________________________________________________________________
void AliITSgeomSPD::SetShape(const char *name,const char *title,
			     const char * /*mat*/,Float_t dx,Float_t dy,Float_t dz){
    // Delete any existing shape info and replace it with a new
    // shape information.
    // Inputs:
    //   char * name  Name of the shape
    //   char * title Title of the shape
    //   char * mat   Material name for the shape
    //   Float_t dx   half width of the shape [cm]
    //   Float_t dy   half thickness of the shape [cm]
    //   Float_t dz   half length of the shape [cm]
    // Outputs:
    //   none.
    // Return:
    //   none.

    fName  = name;
    fTitle = title;
    fDx    = dx;
    fDy    = dy;
    fDz    = dz;
    return;
}
//______________________________________________________________________
void AliITSgeomSPD::LToDet(Float_t xl,Float_t zl,Int_t &row,Int_t &col){
// Returns the row and column pixel numbers for a given local coordinate
// system. If they are outside then it will return -1 or fNbinx/z.
    Int_t i;

    if(xl<fLowBinEdgeX[0]) row = -1;
    else{
	for(i=0;i<fNbinx;i++) if(xl<=fLowBinEdgeX[i]) break;
	row = i;
    } //end if too low.
    if(zl<fLowBinEdgeX[0]) col = -1;
    else{
	for(i=0;i<fNbinz;i++) if(zl<=fLowBinEdgeZ[i]) break;
	col = i;
    } //end if too low.
    return;
}
//______________________________________________________________________
void AliITSgeomSPD::DetToL(Int_t row,Int_t col,Float_t &xl,Float_t &zl){
// returns the pixel center local coordinate system location for a given
// row and column number. It the row or column number is outside of the 
// defined range then it will return the nearest detector edge.

    if(row>=0||row<fNbinx-1) xl = 0.5*(fLowBinEdgeX[row]+fLowBinEdgeX[row+1]);
    else if(row<0) xl = fLowBinEdgeX[0];else xl = fLowBinEdgeX[fNbinx-1];
    if(col>=0||col<fNbinz-1) zl = 0.5*(fLowBinEdgeZ[col]+fLowBinEdgeZ[col+1]);
    else if(col<0) zl = fLowBinEdgeZ[0];else zl = fLowBinEdgeZ[fNbinz-1];
    return;
}
//______________________________________________________________________
void AliITSgeomSPD::Print(ostream *os) const {
// Standard output format for this class
    Int_t i;
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

    fmt = os->setf(ios::scientific); // set scientific floating point output
    *os << "TBRIK" << " ";
    *os << setprecision(16) << GetDx() << " ";
    *os << setprecision(16) << GetDy() << " ";
    *os << setprecision(16) << GetDz() << " ";
    *os << fNbinx-1 << " " << fNbinz-1 << " ";
    for(i=0;i<fNbinx;i++) *os << setprecision(16) << fLowBinEdgeX[i] << " ";
    for(i=0;i<fNbinz;i++) *os << setprecision(16) << fLowBinEdgeZ[i] << " ";
    *os << endl;
    os->flags(fmt);
    return;
}
//______________________________________________________________________
void AliITSgeomSPD::Read(istream *is){
// Standard input format for this class
    Int_t i,j;
    Float_t dx,dy,dz;
    char shape[20];

    for(i=0;i<20;i++) shape[i]='\0';
    *is >> shape;
    if(strcmp(shape,"TBRIK")) Warning("::Read","Shape not a TBRIK");
    *is >> dx >> dy >> dz;
    SetShape("ActiveSPD","Active volume of SPD","",dx,dy,dz);
    *is >> i >> j;
    SetNbinX(i);
    SetNbinZ(j);
    InitLowBinEdgeX();
    InitLowBinEdgeZ();
    for(i=0;i<fNbinx;i++) *is >> fLowBinEdgeX[i];
    for(i=0;i<fNbinz;i++) *is >> fLowBinEdgeZ[i];
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSPD &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSPD &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//=====================================================================

ClassImp(AliITSgeomSPD300)

AliITSgeomSPD300::AliITSgeomSPD300() : AliITSgeomSPD(){
////////////////////////////////////////////////////////////////////////
//    default constructor, for ITS TDR geometry. This only consists of
// a default constructor to construct the defalut TDR SPD detector geometry
// 256 X 279 300 by 50 micron pixels.
////////////////////////////////////////////////////////////////////////
const Float_t kdx=0.6400,kdy=0.0075,kdz=4.1900; // cm; Standard pixel detector
                                                // size is 2dx wide, 2dz long,
                                                // and 2dy thick. Geant 3.12
                                                // units.
const Float_t kbinx0 = 0.0050; // cm; Standard pixel size in x direction.
const Int_t   knbinx = 256;    // number of pixels along x direction.
const Float_t kbinz0 = 0.0300; // cm; Standard pixel size in z direction.
const Float_t kbinz1 = 0.0350; // cm; Edge pixel size in z direction.
const Int_t   knbinz = 279;    // number of pixels along z direction.
    Int_t i;
    Float_t dx=0.0,dz=0.0;

    SetNbinX(knbinx); // default number of bins in x.
    SetNbinZ(knbinz); // default number of bins in z.

    for(i=0;i<knbinx;i++) dx += kbinx0; // Compute size x.
    dx *= 0.5;
    for(i=0;i<knbinz;i++) dz += kbinz0; // Compute size z.
    dz += 2.0*(kbinz1-kbinz0);
    dz *= 0.5;
    InitLowBinEdgeX();
    InitLowBinEdgeZ();
    SetLowBinEdgeX(0,-dx); // Starting position X
    for(i=0;i<knbinx;i++) SetLowBinEdgeX(i+1,GetBinLowEdgeX(i)+kbinx0);
    SetLowBinEdgeZ(0,-dz); // Starting position z
    SetLowBinEdgeZ(1,GetBinLowEdgeZ(0)+kbinz1);
    for(i=1;i<knbinz;i++) SetLowBinEdgeZ(i+1,GetBinLowEdgeZ(i)+kbinz0);
    SetLowBinEdgeZ(knbinz,GetBinLowEdgeZ(knbinz-1)+kbinz1);

    if(TMath::Abs(dx-kdx)>1.0E-4 || TMath::Abs(dz-kdz)>1.0E-4) 
	Warning("Default Creator","Detector size may not be write.");
    SetShape("ActiveSPD","Active volume of SPD","",dx,kdy,dz);
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSPD300 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSPD300 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//=====================================================================

ClassImp(AliITSgeomSPD425Short)

AliITSgeomSPD425Short::AliITSgeomSPD425Short() : AliITSgeomSPD(){
////////////////////////////////////////////////////////////////////////
//    default constructor, for ITS post TDR geometry. This only consists of
// a default constructor to construct the defalut post TDR SPD detector 
// geometry 256 X 197 425 by 50 micron pixels with interleaved 625 by 50
// micron pixels (large detector).
////////////////////////////////////////////////////////////////////////
}
//----------------------------------------------------------------------
AliITSgeomSPD425Short::AliITSgeomSPD425Short(Int_t npar,Float_t *par) :
                                                              AliITSgeomSPD(){
////////////////////////////////////////////////////////////////////////
//    Standard constructor, for ITS post TDR geometry. This only consists of
// a default constructor to construct the defalut post TDR SPD detector 
// geometry 256 X 197 425 by 50 micron pixels with interleaved 625 by 50
// micron pixels (large detector).
////////////////////////////////////////////////////////////////////////

    const Float_t kdx=0.6400/*,kdy=0.015*/,kdz=3.480; // cm; Standard pixel
                                                      // detector size is 2dx
                                                      //  wide, 2dz long, and
                                                      //  2dy thick. Geant 3.12
                                                      // units.
    const Float_t kbinx0 = 0.0050; // cm; Standard pixel size in x direction.
    const Int_t   knbinx = 256;    // number of pixels along x direction.
    const Float_t kbinz0 = 0.0425; // cm; Standard pixel size in z direction.
    const Float_t kbinz1 = 0.0625; // cm; Special pixel size in z direction.
    const Int_t   knbinz = 160;    // number of pixels along z direction.
    Int_t i;
    Float_t dx,dz,*binSizeX,*binSizeZ;

    SetNbinX(knbinx); // default number of bins in x.
    SetNbinZ(knbinz); // default number of bins in z.

    binSizeX = new Float_t[knbinx]; // array of bin sizes along x.
    for(i=0;i<knbinx;i++) binSizeX[i] = kbinx0; // default x bin size.
    binSizeZ = new Float_t[knbinz]; // array of bin sizes along z.
    for(i=0;i<knbinz;i++) binSizeZ[i] = kbinz0; // default z bin size.
    binSizeZ[ 31] = kbinz1;
    binSizeZ[ 32] = kbinz1;

    binSizeZ[ 63] = kbinz1;
    binSizeZ[ 64] = kbinz1;

    binSizeZ[ 95] = kbinz1;
    binSizeZ[ 96] = kbinz1;

    binSizeZ[127] = kbinz1;
    binSizeZ[128] = kbinz1;

    // correct detector size for bin size.
    dx = 0.0;
    for(i=0;i<knbinx;i++) dx += binSizeX[i];
    dx *= 0.5;
    dz = 0.0;
    for(i=0;i<knbinz;i++) dz += binSizeZ[i];
    dz *= 0.5;

    if(npar<3){
	Error("AliITSgeomSPD425Short",
              "npar=%d<3 array par must be at least [3] or larger",npar);
	return;
    } // end if
    SetShape("ActiveSPD","Active volume of SPD","",
	     par[0],par[1],par[2]);
    if(TMath::Abs(dx-kdx)>1.0E-4 || TMath::Abs(dz-kdz)>1.0E-4) 
	Warning("Default Creator","Detector size may not be write.");

    InitLowBinEdgeX(); // array of bin sizes along x.
    InitLowBinEdgeZ(); // array of bin sizes along x.
    SetLowBinEdgeX(0,-dx);
    SetLowBinEdgeZ(0,-dz);
    for(i=0;i<knbinx;i++) SetLowBinEdgeX(i+1,GetBinLowEdgeX(i)+binSizeX[i]);
    for(i=0;i<knbinz;i++) SetLowBinEdgeZ(i+1,GetBinLowEdgeZ(i)+binSizeZ[i]);
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSPD425Short &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSPD425Short &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================

ClassImp(AliITSgeomSPD425Long)

AliITSgeomSPD425Long::AliITSgeomSPD425Long(){
////////////////////////////////////////////////////////////////////////
//    default constructor, for ITS post TDR geometry. This only consists of
// a default constructor to construct the defalut post TDR SPD detector 
// geometry 256 X 197 425 by 50 micron pixels with interleaved 625 by 50
// micron pixels (large detector).
////////////////////////////////////////////////////////////////////////

    const Float_t kdx=0.6400,kdy=0.0075,kdz=4.2650; // cm; Standard pixel
                                                    // detector size is 2dx
                                                    //  wide, 2dz long, and
                                                    //  2dy thick. Geant 3.12
                                                    // units.
    const Float_t kbinx0 = 0.0050; // cm; Standard pixel size in x direction.
    const Int_t   knbinx = 256;    // number of pixels along x direction.
    const Float_t kbinz0 = 0.0425; // cm; Standard pixel size in z direction.
    const Float_t kbinz1 = 0.0625; // cm; Special pixel size in z direction.
    const Int_t   knbinz = 192;    // number of pixels along z direction.
    Int_t i;
    Float_t dx,dz,*binSizeX,*binSizeZ;

    SetNbinX(knbinx); // default number of bins in x.
    SetNbinZ(knbinz); // default number of bins in z.

    binSizeX = new Float_t[knbinx]; // array of bin sizes along x.
    for(i=0;i<knbinx;i++) binSizeX[i] = kbinx0; // default x bin size.
    binSizeZ = new Float_t[knbinz]; // array of bin sizes along z.
    for(i=0;i<knbinz;i++) binSizeZ[i] = kbinz0; // default z bin size.
    binSizeZ[ 31] = kbinz1;
    binSizeZ[ 32] = kbinz1;

    binSizeZ[ 63] = kbinz1;
    binSizeZ[ 64] = kbinz1;

    binSizeZ[ 95] = kbinz1;
    binSizeZ[ 96] = kbinz1;

    binSizeZ[127] = kbinz1;
    binSizeZ[128] = kbinz1;

    binSizeZ[159] = kbinz1;
    binSizeZ[160] = kbinz1;

    // correct detector size for bin size.
    dx = 0.0;
    for(i=0;i<knbinx;i++) dx += binSizeX[i];
    dx *= 0.5;
    dz = 0.0;
    for(i=0;i<knbinz;i++) dz += binSizeZ[i];
    dz *= 0.5;

    SetShape("ActiveSPD","Active volume of SPD","",dx,kdy,dz);
    if(TMath::Abs(dx-kdx)>1.0E-4 || TMath::Abs(dz-kdz)>1.0E-4) 
	Warning("Default Creator","Detector size may not be write.");

    InitLowBinEdgeX(); // array of bin sizes along x.
    InitLowBinEdgeZ(); // array of bin sizes along x.
    SetLowBinEdgeX(0,-dx);
    SetLowBinEdgeZ(0,-dz);
    for(i=0;i<knbinx;i++) SetLowBinEdgeX(i+1,GetBinLowEdgeX(i)+binSizeX[i]);
    for(i=0;i<knbinz;i++) SetLowBinEdgeZ(i+1,GetBinLowEdgeZ(i)+binSizeZ[i]);
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSPD425Long &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSPD425Long &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================
