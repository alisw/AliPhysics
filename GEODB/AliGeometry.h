#ifndef ALIGEOMETRY_H
#define ALIGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// -*- C++ -*-
// 
// 1999/01/05
// ---------------------------------------------------------------------------
//
// AliGeometry Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------
// The Geometry class holds the detector,s geometry. Essentially it has a pointer
// to the Top Level AliGNode and an array of pointers to rules that specify the
// complete desgn below that node.

#include <TFile.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TList.h>
#include "AliGNode.h"

const Int_t kMAXLEVELS  = 20;
const Int_t kMatrixSize = 16;
const Int_t kVectorSize =  4;

class AliGeometry: public TNamed {

    protected:
        Float_t  fBomb;     //Bomb factor for exploded geometry
        Int_t    fGeomLevel;
        Int_t    fMaxDepth;
        TList*   fRules;
        TString  fTopNode;
        TList*   fTransformation;

    public:
        AliGeometry( Text_t* name, Text_t* title, AliGNode* topNode, int maxDepth=0 ); // Constructor
        AliGeometry( AliGeometry* Geom=NULL ); // Copy or Default Constructor
        virtual ~AliGeometry(); // Destructor
        AliGeometry* operator=( const AliGeometry* Geom );

                TList*       GetfTransf() const {return fTransformation;}
                Float_t      GetBomb() const {return fBomb;}
                TString      GetTop() const {return fTopNode;}
                TList*       GetRules() const {return fRules;}
                AliGNode*    FileToMemTree( TFile* file );
                //AliGNode*   FileToMemTree( Text_t* root_file, Text_t* geom_file );
                void         Local2Master( Float_t *local, Float_t *master);
        virtual Int_t        PopLevel(){return fGeomLevel>0?fGeomLevel--:0;}
        virtual Int_t        PushLevel(){return fGeomLevel++;}
                TVector*     PopMatrix();
                void         PushMatrix(TVector* matrix);
                AliGNode*    RecurToMem( const char* node_name, TFile* file, AliGNode* father, char* subtransf );
                void         RulesList( TList *fRules, AliGNode* node, int position );
        virtual void         SetBomb(Float_t bomb=1.4) {fBomb = bomb;}
                void         UpdateMatrix(AliGTransform* trans);


    ClassDef(AliGeometry,1) //Generic Geometry class
};

R__EXTERN AliGeometry *gAliGeometry;
R__EXTERN TVector* gMatrix;

#endif
