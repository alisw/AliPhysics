#ifndef ALIGSUPERSHAPE_H
#define ALIGSUPERSHAPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TObjArray.h>
#include <TNamed.h>
#include "AliGShape.h"
#include "AliGTransform.h"

class AliGSuperShape : public AliGShape {
    protected:
        TString    fExpression;
        TObjArray* fShapes;
        TObjArray* fTransf;

    public:
        AliGSuperShape(); /* Default Constructor */
        AliGSuperShape( Text_t* name, Text_t* title, AliGShape* shapes=NULL, AliGTransform* trans=NULL, Text_t* expression=NULL); /* Constructor */
        virtual ~AliGSuperShape(); /* Destructor */

        void Add( AliGShape* shape );
        void Add( AliGTransform* trans );
        void Add( Text_t* expression );

    ClassDef(AliGSuperShape,1) //SuperShape class
};
#endif




