#ifndef FLUKAVOLUME_H
#define FLUKAVOLUME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"

class FlukaVolume : public TNamed 
{
public:
    FlukaVolume() {;}
    FlukaVolume(const char* name, const Int_t med);
    virtual ~FlukaVolume(){;}
    Int_t GetMedium() {return fMedium;}
private:
    Int_t    fMedium;
    ClassDef(FlukaVolume,1) // Transient storage for volume information
};
#endif
