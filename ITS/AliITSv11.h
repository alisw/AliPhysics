#ifndef ALIITSV11_H
#define ALIITSV11_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/////////////////////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 11, 2003 geometry    //
/////////////////////////////////////////////////////////////////////////
 
#include "AliITS.h"
class TGeoVolume;
class TGeoPcon;
class AliITSv11GeometrySPD;
class AliITSv11GeometrySDD;
class AliITSv11GeometrySupport;
 
class AliITSv11 : public AliITS {

 public:
    AliITSv11();
    AliITSv11(Int_t debugITS,Int_t debugSPD=0,Int_t debugSDD=0,
              Int_t debugSSD=0,Int_t debugSUP=0);
    AliITSv11(const char *title);
    AliITSv11(const AliITSv11 &source); // copy constructor
    AliITSv11& operator=(const AliITSv11 &source); // assignment operator
    virtual       ~AliITSv11();
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {return 11;} // ITS version number
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   DrawModule();
    virtual void   StepManager();
    // Getters for different parts of the geometry initilization
    virtual AliITSv11GeometrySPD*     GetSPDGeometry(){return fSPDgeom;} // SPD
    virtual AliITSv11GeometrySupport* GetSupGeometry(){return fSupgeom;} //Sup.
  private:
    void InitAliITSgeom();

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    Float_t  fDet1;	      // thickness of detector in SPD layer 1
    Float_t  fDet2;	      // thickness of detector in SPD layer 2
    Float_t  fChip1;	      // thickness of chip in SPD layer 1   
    Float_t  fChip2;	      // thickness of chip in SPD layer 2   
    Int_t    fRails;          // switch rails on (=1) and off (=0)
    Int_t    fFluid;          // switch between water(=1) and freon(=0)
    //
    AliITSv11GeometrySPD     *fSPDgeom;  //! Pointer to SPD Geometry
    AliITSv11GeometrySDD     *fSDDgeom;  //! Pointer to SDD Geometry
    AliITSv11GeometrySupport *fSupgeom;  //! Pointer to ITS support Geometry

    ClassDef(AliITSv11,1)  //Hits manager for set:ITS version 11
};
 
#endif
