#ifndef ALIGENCOCKTAILEVENTHEADER_H
#define ALIGENCOCKTAILEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenCocktailEventHeader : public AliGenEventHeader
{
 public:
    AliGenCocktailEventHeader();
    AliGenCocktailEventHeader(const char* name);
    virtual ~AliGenCocktailEventHeader() {}
    virtual void AddHeader(AliGenEventHeader* header);
    virtual TList* GetHeaders() {return fHeaders;}
    
protected:
    TList  *fHeaders;     // List of Headers
    ClassDef(AliGenCocktailEventHeader,1)  // Event header for Cocktail event
};

#endif
