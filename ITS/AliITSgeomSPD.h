#ifndef ALIITSGEOMSPD_H
#define ALIITSGEOMSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TBRIK.h>

class ostream;
class istream;

class TShape;

class AliITSgeomSPD : public TObject {

 public:
    AliITSgeomSPD();
    AliITSgeomSPD(Float_t dy,Int_t nx,Float_t *bx,Int_t nz,Float_t *bz);
    AliITSgeomSPD(AliITSgeomSPD &source);
    AliITSgeomSPD& operator=(AliITSgeomSPD &source);
    virtual ~AliITSgeomSPD();
    virtual void ReSetBins(Float_t dy,Int_t nx,Float_t *bx,
			   Int_t nz,Float_t *bz);
    virtual TShape *GetShape() const {return fShapeSPD;}
    virtual Float_t GetDx() const { // Get TBRIK Dx
        if(fShapeSPD!=0) return fShapeSPD->GetDx();
        else return 0.0;}
    virtual Float_t GetDy() const {// Get TBRIK Dy
        if(fShapeSPD!=0) return fShapeSPD->GetDy();
        else return 0.0;}
    virtual Float_t GetDz() const {// Get TBRIK Dz
        if(fShapeSPD!=0) return fShapeSPD->GetDz();
        else return 0.0;}
    virtual Int_t GetNbinxX() const {return fNbinx-1;} // returns the number of bins x
    virtual Int_t GetNbinxZ() const {return fNbinz-1;} // returns the number of bins z
    virtual Float_t GetBinSizeX(Int_t i) const
	{if(i<fNbinx-1&&i>=0) return fLowBinEdgeX[i+1]-fLowBinEdgeX[i];
	else return 0.0;}; // give size of bin i in x.
    virtual Float_t GetBinSizeZ(Int_t i) const
	{if(i<fNbinx-1&&i>=0) return fLowBinEdgeZ[i+1]-fLowBinEdgeZ[i];
	else return 0.0;}; // give size of bin i in z.
    virtual Float_t GetBinLowEdgeX(Int_t i) const
	{if(i<fNbinx-1&&i>=0) return fLowBinEdgeX[i];
	else return 0.0;}; // give size of bin i in x.
    virtual Float_t GetBinLowEdgeZ(Int_t i) const
	{if(i<fNbinz-1&&i>=0) return fLowBinEdgeZ[i];
	else return 0.0;}; // give size of bin i in z.
    virtual void InitLowBinEdgeX(){// allocate memory for fLowBinEdgeX.
	if(fLowBinEdgeX!=0) delete[] fLowBinEdgeX;
	if(fNbinx>0) fLowBinEdgeX = new Float_t[fNbinx];else fLowBinEdgeX = 0;}
    virtual void InitLowBinEdgeZ(){// allocate memory for fLowBinEdgeZ.
	if(fLowBinEdgeZ!=0) delete[] fLowBinEdgeZ;
	if(fNbinz>0) fLowBinEdgeZ = new Float_t[fNbinz];else fLowBinEdgeZ = 0;}
    virtual void SetShape(const char *name,const char *title,const char *mat,
			  Float_t dx,Float_t dy,Float_t dz)
	{fShapeSPD = new TBRIK(name,title,mat,dx,dy,dz);};
    virtual void SetNbinX(Int_t i){fNbinx = i+1;} // Set nubmer of pixels in x
    virtual void SetNbinZ(Int_t i){fNbinz = i+1;} // Set nubmer of pixels in z
    virtual void SetLowBinEdgeX(Int_t i,Float_t s){//puts value in fLowbinEdgeX
        if(i>=0&&i<fNbinx) fLowBinEdgeX[i] = s;}
    virtual void SetLowBinEdgeZ(Int_t i,Float_t s){//puts value in fLowbinEdgeZ
	if(i>=0&&i<fNbinz) fLowBinEdgeZ[i] = s;};
    virtual void LToDet(Float_t xl,Float_t zl,Int_t &row,Int_t &col);
    virtual void DetToL(Int_t row,Int_t col,Float_t &xl,Float_t &zl);
    virtual void Print(ostream *os) const; // output streamer to standard out.
    virtual void Read(istream *is); // input streamer from standard in.

 protected:
    TBRIK  *fShapeSPD; // SPD active area shape
    Int_t fNbinx;  // the number of elements in fLowBinEdgeX (#bins-1)
    Int_t fNbinz;  // the number of elements in fLowBinEdgeZ (#bins-1)
    Float_t *fLowBinEdgeX; //[fNbinx] Array of X lower bin edges for the pixels
    Float_t *fLowBinEdgeZ; //[fNbinz] Array of Z lower bin edges for the pixels

    ClassDef(AliITSgeomSPD,1) // ITS SPD detector geometry class..

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSPD &source);
istream &operator>>(istream &os,AliITSgeomSPD &source);
#endif
//======================================================================

#ifndef ALIITSGEOMSPD300_H
#define ALIITSGEOMSPD300_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//#include "AliITSgeomSPD.h"

class ostream;
class istream;

class AliITSgeomSPD300 : public AliITSgeomSPD {

 public:
    AliITSgeomSPD300();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SPD detector geometry.
    ClassDef(AliITSgeomSPD300,0) // ITS SPD detector geometry class for 300X50 micron pixel size.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSPD300 &source);
istream &operator>>(istream &os,AliITSgeomSPD300 &source);
#endif
//======================================================================

#ifndef ALIITSGEOMSPD425SHORT_H
#define ALIITSGEOMSPD425SHORT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//#include "AliITSgeomSPD.h"


// temporary - this will migrate into the segmentation class

class ostream;
class istream;

class AliITSgeomSPD425Short : public AliITSgeomSPD {

 public:
    AliITSgeomSPD425Short();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SPD detector geometry.
    ClassDef(AliITSgeomSPD425Short,0) // ITS SPD detector geometry class for 425X50 micron pixel size Short version.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSPD425Short &source);
istream &operator>>(istream &os,AliITSgeomSPD425Short &source);
#endif
//======================================================================

#ifndef ALIITSGEOMSPD425LONG_H
#define ALIITSGEOMSPD425LONG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//#include "AliITSgeomSPD.h"


// temporary - this will migrate into the segmentation class

class ostream;
class istream;

class AliITSgeomSPD425Long : public AliITSgeomSPD {

 public:
    AliITSgeomSPD425Long();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SPD detector geometry.
    ClassDef(AliITSgeomSPD425Long,0) // ITS SPD detector geometry class for 425X50 micron pixel size.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSPD425Long &source);
istream &operator>>(istream &os,AliITSgeomSPD425Long &source);
#endif
