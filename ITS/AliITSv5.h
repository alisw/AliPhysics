#ifndef ALIITSV5_H
#define ALIITSV5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5
////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TBRIK.h"
#include "AliITS.h"
#include "AliITSgeom.h"

class AliITSv5 : public AliITS {

 public:
    AliITSv5();
    AliITSv5(const char *name, const char *title);
    AliITSv5(const AliITSv5 &source); // copy constructor
    AliITSv5& operator=(const AliITSv5 &source); // assignment operator	 
    virtual       ~AliITSv5() ;
    virtual void   BuildGeometry();
    virtual void  CreateGeometry();
    virtual void  CreateMaterials();
    virtual void  Init();   
    virtual Int_t IsVersion() const {
	                                  // returns the ITS version number
	                                  return 5;
												}
    virtual void  StepManager();

 private:
    Int_t fId5N; // The number of layers for geometry version 5
                 // The name of the layers as defined in the Geant tree.
    char  **fId5Name; // They are the names of the sensitive volumes



    ClassDef(AliITSv5,1)//Hits manager for ITS version 5 Official detailed geometry
};
 
#endif
