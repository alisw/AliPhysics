#ifndef ALINODE_H
#define ALINODE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNode.h"

class AliNode : public TNode 
{
public:
    AliNode(){}
    AliNode(const char* name, const char* title, const char* shapename,
	    Double_t x = 0, Double_t y = 0, Double_t z = 0, const char* matrixname="",
	    Option_t* option="");

    AliNode(const char* name, const char* title, TShape* shape,
	    Double_t x = 0, Double_t y = 0, Double_t z = 0, TRotMatrix* matrix = 0,
	    Option_t* option="");
    AliNode(const AliNode &node, AliNode* parent);
    
    virtual ~AliNode(){}

    virtual void SetDivision(Int_t ndiv, Int_t axis, Float_t start, Float_t step);
    virtual void ExpandDivisions();
    virtual Int_t   Axis()   const {return fAxis;}
    virtual Int_t   Ndiv()   const {return fNDivision;}
    virtual Float_t Step()   const {return fStep;}
    virtual Float_t StartC() const {return fStartC;}
    virtual void    AddSons(TList* list);
    virtual void    AddSon(AliNode* node);
    
	    
	    
private:
    Int_t   fAxis;         // division axis
    Int_t   fNDivision;    // number of divisions
    Float_t fStep;         // number of steps
    Float_t fStartC;       // start coordinate

    AliNode &operator=(const AliNode &) {return *this;}

    ClassDef(AliNode,1) // Material Object for GUI 
};

#endif








