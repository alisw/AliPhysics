#ifndef ALIITSVPPRCOARSEASYMM_H
#define ALIITSVPPRCOARSEASYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 6    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSvPPRcoarseasymm : public AliITS {

 public:
    AliITSvPPRcoarseasymm();
    AliITSvPPRcoarseasymm(const char *name, const char *title);
    AliITSvPPRcoarseasymm(const AliITSvPPRcoarseasymm &source); // copy constructor
    AliITSvPPRcoarseasymm& operator=(const AliITSvPPRcoarseasymm &source); // assignment operator
    virtual       ~AliITSvPPRcoarseasymm() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual Int_t  IsVersion() const {// returns the ITS version number 
	return 6;}
    virtual void   DrawModule();
    virtual void   StepManager();
    virtual Int_t GetMajorVersion(){// return Major Version Number
	return fMajorVersion;}
    virtual Int_t GetMinorVersion(){// return Major Version Number
	return fMinorVersion;}


 private:
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number

    ClassDef(AliITSvPPRcoarseasymm,1)  //Hits manager for set:ITS version 6 
                                       // PPR coarse Geometry asymmetric
};
 
#endif
