#ifndef ALIGBOX_H
#define ALIGBOX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGBox Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------
// AliGBox is a subclass of AliGShape. Its dimensions are:
//      - Dx    half-length of the bosx along X-axis
//      - Dy    hlf-length of the box along Y-axis
//      - Dz    half-length of the box along Z-axis */


#include "AliGShape.h"

class AliGBox: public AliGShape {

    protected:
        Float_t fDx; /* X half Dimension */
        Float_t fDy; /* Y half Dimension */
        Float_t fDz; /* Z half Dimension */

    public:
        AliGBox( Text_t* name,Text_t* title, Float_t Dx, Float_t Dy, Float_t Dz ); /* Constructor */
        AliGBox( AliGBox* box=NULL ); /* Copy or Default Constructor */
        AliGBox* operator=( const AliGBox *box ); /* Operator = */
        virtual ~AliGBox() {} /* Destructor */

        //virtual void Draw( Option_t* option );
        virtual void    DrawShape( Option_t* option ); // *MENU*
	        Float_t GetX() {return fDx;}
	        Float_t GetY() {return fDy;}
	        Float_t GetZ() {return fDz;}
        virtual void    Paint( Option_t* option );
        virtual void    PaintGLPoints( Float_t* vertex );
        virtual void    SetPoints( Float_t* buff );
	        void    SetX(Float_t Dx) {fDx = Dx;}
	        void    SetY(Float_t Dy) {fDy = Dy;}
	        void    SetZ(Float_t Dz) {fDz = Dz;}
        virtual void    Sizeof3D() const;
        //void Streamer(TBuffer &b);

    ClassDef(AliGBox,1) // Simple box class
};

#endif

