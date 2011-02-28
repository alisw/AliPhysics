#ifndef ALIITSV11_H
#define ALIITSV11_H
/* Copyright(c) 2007-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//************************************************************************
//
//                 Inner Traking System geometry v11
//
//  Based on ROOT geometrical modeler
//
// B. Nilsen, L. Gaudichet
//************************************************************************
#include "AliITSInitGeometry.h"
#include "AliITS.h"
class AliITSv11GeometrySPD;
class AliITSv11GeometrySDD;
class AliITSv11GeometrySSD;
class AliITSv11GeometrySupport;

class AliITSv11 : public AliITS {

 public:
    AliITSv11();
    AliITSv11(const char *name, const char *title);
    AliITSv11(Int_t /* debugITS */, Int_t debugSPD = 0, Int_t debugSDD = 0,
	     Int_t debugSSD = 0, Int_t debugSUP = 0);
    virtual       ~AliITSv11() ;
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   StepManager();
    //virtual AliITSv11GeometrySPD*     GetSPDGeometry(){return fSPDgeom;}
    virtual AliITSv11GeometrySDD*    GetSDDGeometry(){return fSDDgeom;}
    //virtual AliITSv11GeometrySupport* GetSupGeometry(){return fSupgeom;}
    virtual Int_t  IsVersion() const { return kv11;}  // ITS version number 
    virtual Int_t  GetMajorVersion() const {// return Major Version Number
                    return fMajorVersion;}
    virtual Int_t  GetMinorVersion() const {// return Major Version Number
                    return fMinorVersion;}
    virtual void   SetMinorVersion(Int_t v=0){ // Choose between existing minor versions
                   fMinorVersion = v;}

    virtual void SetDensityServicesByThickness(){// uses services density
	// calculation based on the thickness of the services.
	fByThick = kTRUE;}
    virtual void SetDensityServicesByMass(){// uses services density
	// calculation based on the Mass of the services.
	fByThick = kFALSE;}

 private:
    AliITSv11(const AliITSv11 &source);            // copy constructor
    AliITSv11& operator=(const AliITSv11 &source); // assignment operator

    Bool_t fByThick;          // Flag to use services materials by thickness
                                // ture, or mass false.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    AliITSv11GeometrySPD *fSPDgeom;      //SPD Geometry
    AliITSv11GeometrySDD *fSDDgeom;      //! SDD Geometry
    AliITSv11GeometrySSD *fSSDgeom;  //SSD Geometry
    AliITSv11GeometrySupport *fSupgeom;  //Support Geometry
    AliITSInitGeometry fIgm; //! Geometry initlization object

    ClassDef(AliITSv11,1)  // ITS version 11 
};
 
#endif
