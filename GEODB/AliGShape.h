#ifndef ALIGSHAPE_H
#define ALIGSHAPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* -*- C++ -*- */
/*             */
/* 1998/10/22  */
/* --------------------------------------------------------------------------  */
/*                                                                             */
/* AliGShape Class                                                             */
/*                                                                             */
/* This file is part of the ALICE Geometry Database .                          */
/*                                                                             */
/* Author:  Joana E. Santo                                                     */
/*                                                                             */
/* --------------------------------------------------------------------------  */
/* The shape virtual class represents the basic elements used in this          */
/* geometry to build the constituents of a detector. Deriving from this        */
/* class there are simple shapes like box, sphere, cone and thorus,            */
/* supershapes that are obtained from the previous by Boolean operations       */
/* and also user defined shapes.                                               */

#ifndef ROOT_TNamed
#include <TNamed.h>

#endif

//#ifndef ROOT_TMaterial
//#include <TMaterial.h>
//#endif

#ifndef ROOT_TAttLine
#include <TAttLine.h>
#endif

#ifndef ROOT_TAttFill
#include <TAttFill.h>
#endif

#ifndef ROOT_X3DBuffer
#include <X3DBuffer.h>
#endif

#ifndef ROOT_TPolyLine3D
#include <TPolyLine3D.h>
#endif

#include "TVector.h"
class AliGNode;
extern "C" { void FillX3DBuffer (X3DBuffer *buff); }

class AliGShape: public TNamed, public TAttLine, public TAttFill { 
    protected:
        Int_t fColor;

    public:
        AliGShape(Text_t* name, Text_t* title);
        AliGShape();
        AliGShape( AliGShape* shape ); // Copy Constructor
        virtual ~AliGShape() {}

                Int_t   DistancetoPrimitive(Int_t, Int_t);
                Int_t   GetCol() { return fColor; }
        virtual Bool_t  Is3D () {return kTRUE;}
        virtual void    Paint(Option_t *option="");
        virtual void    PaintGLPoints(Float_t *vertex);
        virtual void    PaintShape(X3DBuffer *buff, Bool_t rangeView);
                void    SetCol( Int_t color ) { fColor = color; }
        virtual void    SetName(const Text_t *name);
        virtual void    SetPoints(Float_t *buffer);

    ClassDef(AliGShape,1) //Generic shape class
};

R__EXTERN AliGNode* gNode;
R__EXTERN TVector*  gMatrix;
R__EXTERN Size3D    gSize3D;

inline void AliGShape::PaintGLPoints(Float_t *) { }
inline void AliGShape::SetName(const Text_t *) { }

#endif
