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
Revision 1.17  2002/10/22 14:45:41  alibrary
Introducing Riostream.h

Revision 1.16  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.14.6.1  2002/06/10 17:51:15  hristov
Merged with v3-08-02

Revision 1.15  2002/05/19 18:17:03  hristov
Changes needed by ICC/IFC compiler (Intel)

Revision 1.14  2001/10/19 21:32:35  nilsen
Minor changes to remove compliation warning on gcc 2.92.2 compiler, and
cleanded up a little bit of code.

Revision 1.13  2001/10/12 22:07:20  nilsen
A patch for C++ io manipulation functions so that they will work both
with GNU gcc 2.96 and GNU gcc 3.01 compilers. Needs to be tested with
other platforms.

Revision 1.12  2001/08/24 21:06:37  nilsen
Added more documentation, fixed up some coding violations, and some
forward declorations.

Revision 1.11  2001/05/16 08:17:49  hristov
Bug fixed in the StepManager to account for the difference in the geometry 
tree for the ITS pixels. This fixes both the funny distribution of pixel 
coordinates and the missing hits/digits/points in many sectors of the ITS 
pixel barrel. Also included is a patch to properly get and use the detector 
dimensions through out the ITS code. (B.Nilsen)

Revision 1.10  2001/04/26 22:44:34  nilsen
Bug fix.

Revision 1.9  2001/02/09 00:00:57  nilsen
Fixed compatibility problem with HP unix {ios::fmtflags -> Int_t}. Fixed
bugs in iostream based streamers used to read and write .det files. Fixed
some detector sizes. Fixed bugs in some default-special constructors.

Revision 1.8  2001/02/03 00:00:30  nilsen
New version of AliITSgeom and related files. Now uses automatic streamers,
set up for new formatted .det file which includes detector information.
Additional smaller modifications are still to come.

*/

////////////////////////////////////////////////////////////////////////
// This class is for the Silicon Pixel Detector, SPD, specific geometry.
// It is being replaced by AliITSsegmentationSPD class. This file also
// constains classes derived from AliITSgeomSPD which do nothing but
// initilize this one with predefined values.
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TShape.h>
#include <TMath.h>

#include "AliITSgeomSPD.h"

ClassImp(AliITSgeomSPD)

AliITSgeomSPD::AliITSgeomSPD(){
// Default Constructor. Set everthing to null.

    fShapeSPD    = 0;
    fNbinx       = 0;
    fNbinz       = 0;
    fLowBinEdgeX = 0;
    fLowBinEdgeZ = 0;
}
//______________________________________________________________________
AliITSgeomSPD::AliITSgeomSPD(Float_t dy,Int_t nx,Float_t *bx,
			     Int_t nz,Float_t *bz){
// Standard Constructor. Set everthing to null.

    fShapeSPD    = 0;
    fNbinx       = 0;
    fNbinz       = 0;
    fLowBinEdgeX = 0;
    fLowBinEdgeZ = 0;
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

    delete fShapeSPD; // delete existing shape
    if(this->fLowBinEdgeX) delete[] this->fLowBinEdgeX; // delete existing
    if(this->fLowBinEdgeZ) delete[] this->fLowBinEdgeZ; // delete existing

    SetNbinX(nx);
    SetNbinZ(nz);
    InitLowBinEdgeX();
    InitLowBinEdgeZ();
    fLowBinEdgeX[0] = -dx;
    fLowBinEdgeZ[0] = -dz;
    for(i=0;i<nx;i++) fLowBinEdgeX[i+1] = fLowBinEdgeX[i] + bx[i];
    for(i=0;i<nz;i++) fLowBinEdgeZ[i+1] = fLowBinEdgeZ[i] + bz[i];
    SetShape("ActiveSPD","Active volume of SPD","SPD SI DET",dx,dy,dz);
    return;
}
//______________________________________________________________________
AliITSgeomSPD::AliITSgeomSPD(AliITSgeomSPD &source){
    // Copy constructor

    *this = source; // just use the = operator for now.
    return;
}
//______________________________________________________________________
AliITSgeomSPD& AliITSgeomSPD::operator=(AliITSgeomSPD &source){
    // = operator
    Int_t i;

    if(&source == this) return *this;
    this->fShapeSPD = new TBRIK(*(source.fShapeSPD));
    if(this->fLowBinEdgeX) delete[] this->fLowBinEdgeX;
    if(this->fLowBinEdgeZ) delete[] this->fLowBinEdgeZ;
    this->fNbinx = source.fNbinx;
    this->fNbinz = source.fNbinz;
    this->InitLowBinEdgeX();
    this->InitLowBinEdgeZ();
    for(i=0;i<fNbinx;i++) this->fLowBinEdgeX[i] = source.fLowBinEdgeX[i];
    for(i=0;i<fNbinz;i++) this->fLowBinEdgeZ[i] = source.fLowBinEdgeZ[i];
    return *this;
}
//______________________________________________________________________
AliITSgeomSPD::~AliITSgeomSPD(){
// Destructor

    delete fShapeSPD;
    if(this->fLowBinEdgeX) delete[] this->fLowBinEdgeX;
    if(this->fLowBinEdgeZ) delete[] this->fLowBinEdgeZ;
    fShapeSPD    = 0;
    fNbinx       = 0;
    fNbinz       = 0;
    fLowBinEdgeX = 0;
    fLowBinEdgeZ = 0;
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
#if defined __ICC || defined __ECC
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
    if(fShapeSPD!=0) delete fShapeSPD;
    SetShape("ActiveSPD","Active volume of SPD","SPD SI DET",dx,dy,dz);
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

/*
$Log$
Revision 1.17  2002/10/22 14:45:41  alibrary
Introducing Riostream.h

Revision 1.16  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.14.6.1  2002/06/10 17:51:15  hristov
Merged with v3-08-02

Revision 1.15  2002/05/19 18:17:03  hristov
Changes needed by ICC/IFC compiler (Intel)

Revision 1.14  2001/10/19 21:32:35  nilsen
Minor changes to remove compliation warning on gcc 2.92.2 compiler, and
cleanded up a little bit of code.

Revision 1.13  2001/10/12 22:07:20  nilsen
A patch for C++ io manipulation functions so that they will work both
with GNU gcc 2.96 and GNU gcc 3.01 compilers. Needs to be tested with
other platforms.

Revision 1.12  2001/08/24 21:06:37  nilsen
Added more documentation, fixed up some coding violations, and some
forward declorations.

Revision 1.11  2001/05/16 08:17:49  hristov
Bug fixed in the StepManager to account for the difference in the geometry tree for the ITS pixels. This fixes both the funny distribution of pixel coordinates and the missing hits/digits/points in many sectors of the ITS pixel barrel. Also included is a patch to properly get and use the detector dimensions through out the ITS code. (B.Nilsen)

Revision 1.10  2001/04/26 22:44:34  nilsen
Bug fix.

Revision 1.9  2001/02/09 00:00:57  nilsen
Fixed compatibility problem with HP unix {ios::fmtflags -> Int_t}. Fixed
bugs in iostream based streamers used to read and write .det files. Fixed
some detector sizes. Fixed bugs in some default-special constructors.

Revision 1.8  2001/02/03 00:00:30  nilsen
New version of AliITSgeom and related files. Now uses automatic streamers,
set up for new formatted .det file which includes detector information.
Additional smaller modifications are still to come.

Revision 1.7  2000/10/02 16:32:35  barbera
Forward declaration added

Revision 1.1.2.8  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.6  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.4  2000/06/10 20:34:37  nilsen
Fixed compilation warning with HP unix.

Revision 1.3  2000/06/10 10:43:04  nilsen
Fixed bug in copy and operator =.

*/

//#include "AliITSgeomSPD300.h"

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

//    cout << "AliITSgeomSPD300 default creator called: start" << endl;

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
    SetShape("ActiveSPD","Active volume of SPD","SPD SI DET",dx,kdy,dz);
//    cout << "AliITSgeomSPD300 default creator called: end" << endl;
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
/*
$Log$
Revision 1.17  2002/10/22 14:45:41  alibrary
Introducing Riostream.h

Revision 1.16  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.14.6.1  2002/06/10 17:51:15  hristov
Merged with v3-08-02

Revision 1.15  2002/05/19 18:17:03  hristov
Changes needed by ICC/IFC compiler (Intel)

Revision 1.14  2001/10/19 21:32:35  nilsen
Minor changes to remove compliation warning on gcc 2.92.2 compiler, and
cleanded up a little bit of code.

Revision 1.13  2001/10/12 22:07:20  nilsen
A patch for C++ io manipulation functions so that they will work both
with GNU gcc 2.96 and GNU gcc 3.01 compilers. Needs to be tested with
other platforms.

Revision 1.12  2001/08/24 21:06:37  nilsen
Added more documentation, fixed up some coding violations, and some
forward declorations.

Revision 1.11  2001/05/16 08:17:49  hristov
Bug fixed in the StepManager to account for the difference in the geometry tree for the ITS pixels. This fixes both the funny distribution of pixel coordinates and the missing hits/digits/points in many sectors of the ITS pixel barrel. Also included is a patch to properly get and use the detector dimensions through out the ITS code. (B.Nilsen)

Revision 1.10  2001/04/26 22:44:34  nilsen
Bug fix.

Revision 1.9  2001/02/09 00:00:57  nilsen
Fixed compatibility problem with HP unix {ios::fmtflags -> Int_t}. Fixed
bugs in iostream based streamers used to read and write .det files. Fixed
some detector sizes. Fixed bugs in some default-special constructors.

Revision 1.8  2001/02/03 00:00:30  nilsen
New version of AliITSgeom and related files. Now uses automatic streamers,
set up for new formatted .det file which includes detector information.
Additional smaller modifications are still to come.

Revision 1.7  2000/10/02 16:32:35  barbera
Forward declaration added

Revision 1.1.2.8  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.6  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.4  2000/06/10 20:34:22  nilsen
Fixed compilation warning with HP unix.

Revision 1.3  2000/06/10 10:42:49  nilsen
Fixed bug in copy and operator =.


*/

//#include "AliITSgeomSPD425Short.h"

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

    SetShape("ActiveSPD","Active volume of SPD","SPD SI DET",
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

/*
$Log$
Revision 1.17  2002/10/22 14:45:41  alibrary
Introducing Riostream.h

Revision 1.16  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.14.6.1  2002/06/10 17:51:15  hristov
Merged with v3-08-02

Revision 1.15  2002/05/19 18:17:03  hristov
Changes needed by ICC/IFC compiler (Intel)

Revision 1.14  2001/10/19 21:32:35  nilsen
Minor changes to remove compliation warning on gcc 2.92.2 compiler, and
cleanded up a little bit of code.

Revision 1.13  2001/10/12 22:07:20  nilsen
A patch for C++ io manipulation functions so that they will work both
with GNU gcc 2.96 and GNU gcc 3.01 compilers. Needs to be tested with
other platforms.

Revision 1.12  2001/08/24 21:06:37  nilsen
Added more documentation, fixed up some coding violations, and some
forward declorations.

Revision 1.11  2001/05/16 08:17:49  hristov
Bug fixed in the StepManager to account for the difference in the geometry tree for the ITS pixels. This fixes both the funny distribution of pixel coordinates and the missing hits/digits/points in many sectors of the ITS pixel barrel. Also included is a patch to properly get and use the detector dimensions through out the ITS code. (B.Nilsen)

Revision 1.10  2001/04/26 22:44:34  nilsen
Bug fix.

Revision 1.9  2001/02/09 00:00:57  nilsen
Fixed compatibility problem with HP unix {ios::fmtflags -> Int_t}. Fixed
bugs in iostream based streamers used to read and write .det files. Fixed
some detector sizes. Fixed bugs in some default-special constructors.

Revision 1.8  2001/02/03 00:00:30  nilsen
New version of AliITSgeom and related files. Now uses automatic streamers,
set up for new formatted .det file which includes detector information.
Additional smaller modifications are still to come.

Revision 1.7  2000/10/02 16:32:35  barbera
Forward declaration added

Revision 1.1.2.8  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.6  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.4  2000/06/10 20:34:22  nilsen
Fixed compilation warning with HP unix.

Revision 1.3  2000/06/10 10:42:49  nilsen
Fixed bug in copy and operator =.


*/

//#include "AliITSgeomSPD425Long.h"

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

    SetShape("ActiveSPD","Active volume of SPD","SPD SI DET",dx,kdy,dz);
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
