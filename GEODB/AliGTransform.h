#ifndef ALIGTRANSFORM_H
#define ALIGTRANSFORM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


/* -*- C++ -*-                                                                  */
/*                                                                              */
/* 1999/01/13                                                                   */
/* ---------------------------------------------------------------------------  */
/*                                                                              */
/* AliGTransform Class                                                          */
/*                                                                              */
/* This file is part of the ALICE Geometry Database .                           */
/*                                                                              */
/* Author:  Joana E. Santo                                                      */
/*                                                                              */
/* ---------------------------------------------------------------------------  */
/* This class represents the transformations applied to regions to position them*/
/* in the node structure.                                                       */


#include <TNamed.h>
#include <TVector.h>
#include <TSystem.h>
#include <TString.h>

class AliGTransform: public TNamed {

    private:
        TString  fExpression;
        TVector* fMatrix;
        Float_t  fX;
        Float_t  fY;
        Float_t  fZ;
	Float_t  fTheta;
        Float_t  fPsi;
        Float_t  fPhi;

    public:
        AliGTransform(); /* Default Constructor */
	AliGTransform(AliGTransform *tra);
        AliGTransform( Text_t* name, Text_t* title );
        AliGTransform( Text_t* name, Text_t* title, Text_t *expression );
        AliGTransform( Text_t* name, Text_t* title, Text_t *axis, Float_t angle);
	AliGTransform( Text_t* name, Text_t* title, Float_t theta1,Float_t phi1, 
	                                            Float_t theta2,
						    Float_t phi2, 
						    Float_t theta3,Float_t
phi3);
	                                           
	AliGTransform( Text_t* name, Text_t* title, Float_t a1,Float_t a2,Float_t a3,Float_t b1,Float_t b2,
	Float_t b3,Float_t c1,Float_t c2,Float_t c3,Float_t Dx,Float_t Dy,Float_t Dz);
        virtual ~AliGTransform(); /* Destructor */

        void     BuildMatrix(Float_t Dx=0., Float_t Dy=0., Float_t Dz=0., Float_t
theta=0., Float_t psi=0.,Float_t phi=0. );
        void     CheckExpression();
        TVector* GetMatrix() {return fMatrix;}

    ClassDef(AliGTransform,1) //Transformation class (Rotation and Translation)
};
#endif
