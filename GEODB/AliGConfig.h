#ifndef ALIGConfig_H
#define ALIGConfig_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TString.h>
#include <TStringLong.h>
#include <TDatime.h>

class AliGConfig: public TNamed {

    protected:
        TDatime     fBeg;        // Time validity of node version
        TString     fDetail;     // Level of detail
        TDatime     fEnd;        // Time validity of node version
        //char*       fFormula;    // Formula describing the node structure below         
        TStringLong fFormula;    // Formula describing the node structure below         
        TString     fMaterialName;
        TString     fShapeName;
        TString     fShapeType;

    public:
        AliGConfig( Text_t* name, Text_t* title, TStringLong formula, Text_t* detail=NULL, const Text_t* shapetype=NULL, const Text_t* shapename=NULL, const Text_t* materialname=NULL, Int_t beg=0, Int_t end=0 );
        AliGConfig( AliGConfig* Config=NULL ); /* Copy Constructor */
        virtual ~AliGConfig(); /* Destructor */
        AliGConfig* operator=( AliGConfig* Config ); /* Operator = */

        TStringLong GetFormula()      const { return fFormula;      }
        TString     GetMaterialName() const { return fMaterialName; }
        TString     GetShapeName()    const { return fShapeName;    }
        TString     GetShapeType()    const { return fShapeType;    }

    ClassDef(AliGConfig,1)   //Config class
};
#endif
