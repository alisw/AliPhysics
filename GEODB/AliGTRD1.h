#ifndef AliGTRD1_H
#define AliGTRD1_H

/* -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGTRD1 Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------
//AliGTRD1 is a subclass of AliGShape. Its dimensions are:
//      - Dx    half-length of the bosx along X-axis
//      - Dy    hlf-length of the box along Y-axis
//      - Dz    half-length of the box along Z-axis */


#include "AliGShape.h"

class AliGTRD1: public AliGShape {
    
    protected:
        Float_t fDx2;  //half length in x at the high z surface
        Float_t fDx1;
        Float_t fDy;
        Float_t fDz;

    public:
        AliGTRD1( Text_t* name,Text_t* title, Float_t dx1,Float_t dx2, Float_t dy, Float_t  dz ); /* Constructor */
        AliGTRD1( ); /* Default Constructor */
	AliGTRD1( AliGTRD1* trd1 );
        virtual ~AliGTRD1() {} /* Destructor */
        	
	        Float_t GetDx1() {return fDx1;}
                Float_t GetDx2() {return fDx2;}
	        Float_t GetDy() {return fDy;}
	        Float_t GetDz() {return fDz;}
                void    Paint(Option_t *option);
	        void    SetDx1(Float_t dx1) {fDx1 = dx1;}
	        void    SetDx2(Float_t dx2) {fDx2 = dx2;}
	        void    SetDy(Float_t dy) {fDy = dy;}
	        void    SetDz(Float_t dz) {fDz = dz;}
                void    SetPoints( Float_t* buff );

    ClassDef(AliGTRD1,1) // Simple trapezoid class
};

#endif

