#ifndef ALIITSVTEST_H
#define ALIITSVTEST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//     Manager and hits classes for  ITS version 5
////////////////////////////////////////////////////////////////////////


#include "AliITS.h"

class AliITSvtest : public AliITS {

public:
           AliITSvtest();
			  AliITSvtest(const char *fileeuc,const char *filetme,
			  const char *name, const char *title);
           AliITSvtest(const AliITSvtest &source); // copy constructor
           AliITSvtest& operator=(const AliITSvtest &source); // assignment operator			  
           virtual       ~AliITSvtest() ;
           virtual void  CreateGeometry();
           virtual void  CreateMaterials();
           virtual void  Init();
           virtual Int_t IsVersion() const {
			                                   // returns the ITS version number 
			                                   return -1;
														 }
           virtual void  StepManager();

private:
           Int_t fIdTestN; // The number of layers for test geometry version
                           // The name of the layers as defined in the Geant tree.
           TString *fIdTestName; // They are the names of the sensitive volumes


  
  ClassDef(AliITSvtest,1)  //Hits manager for ITS test version, Private ITS class for different test geometries
};
 
#endif
