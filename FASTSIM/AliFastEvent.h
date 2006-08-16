#ifndef ALIFASTEVENT_H
#define ALIFASTEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TArrayF.h>

class AliFastEvent : public TObject {
 public:
    AliFastEvent();
    virtual ~AliFastEvent(){;}
    virtual void  SetMultiplicty(Int_t mul) 
	{fMultiplicity = mul;}
    virtual Int_t GetMultiplicty() 
	{return fMultiplicity;}
    virtual void SetVertex(const TArrayF &o) 
	{
	    fEventVertex[0] = o.At(0);
	    fEventVertex[1] = o.At(1);
	    fEventVertex[2] = o.At(2);
	}

    virtual void GetVertex(TArrayF &o) const
	{
	    o[0] = fEventVertex.At(0);
	    o[1] = fEventVertex.At(1);
	    o[2] = fEventVertex.At(2);
	}

 protected:
    Int_t     fMultiplicity;    // Event Multiplicity
    TArrayF   fEventVertex;     // Event primary vertex
    
    ClassDef(AliFastEvent,1) // Base class for fast event
};

#endif 



