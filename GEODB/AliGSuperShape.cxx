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



















