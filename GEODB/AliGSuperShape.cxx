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
*/

// -*- C++ -*-
// 
// 1998/10/22
// ---------------------------------------------------------------------------
//
// AliGSuperShape Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//

#include "AliGSuperShape.h"
#include <TString.h>

ClassImp(AliGSuperShape)

//----------------------------------------------------------------------

AliGSuperShape::AliGSuperShape()
{
    /* Default Constructor */
    fExpression = "";
    fName       = "";
    fShapes     = NULL;
    fTitle      = "";
    fTransf     = NULL;
}

//----------------------------------------------------------------------

AliGSuperShape::AliGSuperShape( Text_t* name, Text_t* title, AliGShape* shapes, AliGTransform* trans, Text_t* expression ):AliGShape(name,title)
{
    /* Constructor */
    fExpression = expression;
    fShapes     = new TObjArray();
    fTransf     = new TObjArray();
}

//----------------------------------------------------------------------

AliGSuperShape::~AliGSuperShape(){
    /* Destructor */
    if(fShapes)     delete fShapes;
    if(fTransf)     delete fTransf;
}

//----------------------------------------------------------------------

void AliGSuperShape::Add( AliGShape *shape )
{

}

//----------------------------------------------------------------------

void AliGSuperShape::Add( AliGTransform *trans )
{

}

//----------------------------------------------------------------------

void AliGSuperShape::Add( Text_t *expression )
{

}

//----------------------------------------------------------------------  



















