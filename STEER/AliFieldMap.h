#ifndef ALIFIELDMAP_H
#define ALIFIELDMAP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//
#include <TNamed.h>
#include <TVector.h>

class AliFieldMap : public TNamed
{
  //Alice Magnetic Field with constant mesh

public:
    AliFieldMap();
    AliFieldMap(const char *name, const char *title);
    AliFieldMap(const AliFieldMap &mag);
    virtual ~AliFieldMap();
    void Copy(AliFieldMap &map) const;
    virtual AliFieldMap & operator=(const AliFieldMap &map);

    virtual void Field(Float_t *x, Float_t *b);
    Float_t Bx(const Int_t ix, const Int_t iy, const Int_t iz) {
	return (*fB)(3*(ix*(fZn*fYn)+iy*fZn+iz));
    }
    Float_t By(const Int_t ix, const Int_t iy, const Int_t iz) {
	return (*fB)(3*(ix*(fZn*fYn)+iy*fZn+iz)+1);
    }
    Float_t Bz(const Int_t ix, const Int_t iy, const Int_t iz) {
	return (*fB)(3*(ix*(fZn*fYn)+iy*fZn+iz)+2);
    }

    Bool_t Inside(Float_t x, Float_t y, Float_t z) 
	{ return (x > fXbeg && x <= fXend &&
		  y > fYbeg && y <= fYend &&
		  z > fZbeg && z <= fZend);
	}
    Float_t Xmin()  {return fXbeg;}
    Float_t Xmax()  {return fXend;}
    Float_t DelX()  {return fXdel;}
    Float_t DeliX() {return fXdeli;}
    
    Float_t Ymin()  {return fYbeg;}
    Float_t Ymax()  {return fYend;}
    Float_t DelY()  {return fYdel;}
    Float_t DeliY() {return fYdeli;}
    
    Float_t Zmin()  {return fZbeg;}
    Float_t Zmax()  {return fZend;}
    Float_t DelZ()  {return fZdel;}
    Float_t DeliZ() {return fZdeli;}
    void    SetLimits(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
		      Float_t zmin, Float_t zmax)
	{
	    fXbeg = xmin; fXend = xmax; fYbeg = ymin; fYend = ymax;
	    fZbeg = zmin; fZend = zmax;
	}
 private:
    void    ReadField();
 protected:

    Float_t    fXbeg;     // Start of mesh in x
    Float_t    fYbeg;     // Start of mesh in y
    Float_t    fZbeg;     // Start of mesh in z
    Float_t    fXend;     // End of mesh in x
    Float_t    fYend;     // End of mesh in y
    Float_t    fZend;     // End of mesh in z
    Float_t    fXdel;     // Mesh step in x
    Float_t    fYdel;     // Mesh step in y
    Float_t    fZdel;     // Mesh step in z
    Double_t   fXdeli;    // Inverse of Mesh step in x
    Double_t   fYdeli;    // Inverse of Mesh step in y
    Double_t   fZdeli;    // Inverse of Mesh step in z
    Int_t      fXn;       // Number of mesh points in x
    Int_t      fYn;       // Number of mesh points in y
    Int_t      fZn;       // Number of mesh points in z
    TVector*   fB;        //!Field map
    
    ClassDef(AliFieldMap,1)  //Class for Field Map
};

#endif
