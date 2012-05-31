#ifndef ALIBCM_H
#define ALIBCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector BCM            //
//  andreas.morsch@cern.ch                    //
////////////////////////////////////////////////
 
#include "AliDetector.h"
class AliLoader;
 
class AliBCM : public AliDetector {
 
 public:
    AliBCM();
    AliBCM(const char *name, const char *title);
    virtual           ~AliBCM();
    virtual void       CreateGeometry();
    virtual void       CreateMaterials();
    virtual void       Init();
    virtual void       StepManager();
    virtual void       MakeBranch(Option_t* option);
    virtual void       SetTreeAddress();
    virtual AliLoader* MakeLoader(const char* topfoldername);
    
    virtual Int_t IsVersion() const {return 0;}
 private:
    Int_t fVolId;  // Volume Id of the sensor
   
   ClassDef(AliBCM,1)  // Manager for detector BCM
};

#endif
