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
// 1998/11/25
// ---------------------------------------------------------------------------
//
// AliGConfig Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo/David Collados/Antonino Bajeli
//
// ---------------------------------------------------------------------------

#include <iostream.h>
#include "AliGConfig.h"

ClassImp(AliGConfig)

//-------------------------------------------------------------------------

AliGConfig::AliGConfig( Text_t* name, Text_t* title, TStringLong formula, Text_t* detail, const Text_t* shapetype, const Text_t* shapename, const Text_t* materialname, Int_t beg, Int_t end ) : TNamed(name, title)
{
    /* Constructor */
    fBeg.Set(beg,0000);           // Time validity of node version
    fDetail       = detail;       // Level of detail
    fEnd.Set(end,0000);           // Time validity of node version

    fFormula      = formula;      //Formula describing the node structure below
    fMaterialName = materialname;
    fShapeName    = shapename;
    fShapeType    = shapetype;
}

//-------------------------------------------------------------------------

AliGConfig::AliGConfig( AliGConfig* Config )
{
    if( Config ) {
        /* Copy Constructor */
        fBeg.Copy(Config->fBeg);                // Time validity of node version
        fDetail       = Config->fDetail.Copy(); // Level of detail
        fEnd.Copy(Config->fEnd);                // Time validity of node version

        fFormula      = Config->fFormula.Copy();
        fMaterialName = Config->fMaterialName.Copy();
        fName         = Config->fName.Copy();
        fShapeName    = Config->fShapeName.Copy();
        fShapeType    = Config->fShapeType.Copy();
        fTitle        = Config->fTitle.Copy();
    }
    else {
        /* Default Constructor */
        fBeg.Set();           // Time validity of node version
        fDetail       = "";   // Level of detail
        fEnd.Set();           // Time validity of node version
        fFormula      = "";   // Formula describing the node structure below         
        fMaterialName = "";
        fName         = "";
        fShapeName    = "";
        fShapeType    = "";
        fTitle        = "";
    }
}

//-------------------------------------------------------------------------

AliGConfig::~AliGConfig() 
{
    /* Destructor */
}

//-------------------------------------------------------------------------

AliGConfig* AliGConfig::operator=( AliGConfig* Config )
{
    /* Operator = */
    if( this == Config) return this; // special case.

    Config->fBeg.Copy(this->fBeg);     // Time validity of node version
    fDetail       = Config->fDetail;   // Level of detail
    Config->fEnd.Copy(this->fEnd);     // Time validity of node version

    fFormula      = Config->fFormula.Copy();
    fMaterialName = Config->fMaterialName.Copy();
    fName         = Config->fName.Copy();
    fShapeName    = Config->fShapeName.Copy();
    fShapeType    = Config->fShapeType.Copy();
    fTitle        = Config->fTitle.Copy();

    return this;
}

//-------------------------------------------------------------------------

