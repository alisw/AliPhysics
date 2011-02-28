#ifndef ALIITSV11HYBRID_H
#define ALIITSV11HYBRID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//========================================================================
//
//            Geometry of the Inner Tracking System
//
//  This geometry is a mix between the old geometry (originally coded
//  in AliITSvPPRasymmFMD) and the new TGeo geometry (v11).
// 
// Ludovic Gaudichet  (gaudichet@to.infn.it)
//
//========================================================================


// $Id$

// $Log$
// Revision 1.7  2007/12/17 14:48:24  masera
// Thermal shield between SPD and SDD (M. Sitta)
//
// Revision 1.6  2007/10/21 19:22:53  masera
// Coding conventions
//
// Revision 1.5  2007/08/24 14:32:57  hristov
// Introduction of SPD half-stave volumes, cleaning and new code (in relation to new SPD geometry) in AliITSv11Hybrid (Ludovic)
//
// Revision 1.4  2007/06/28 10:17:25  masera
// Introduction of the new SSD geometry in simulation (AliITSv11Hybrid) and suppression of overlaps between old and new parts
//
// Revision 1.3  2007/05/08 16:57:42  masera
// Updates concerning the geometry: versioning system, new V11hybrid version, bug fixes (B.Nilsend and L. Gaudichet
//

 
#include "AliITS.h"
#include "AliITSInitGeometry.h"

class  AliITSv11GeometrySPD;
class  AliITSv11GeometrySDD;
class  AliITSv11GeometrySSD;
class  AliITSv11GeometrySupport;
class  TGeoVolume;
class TGeoVolumeAssembly;

class AliITSv11Hybrid : public AliITS {

 public:
    AliITSv11Hybrid();
    AliITSv11Hybrid(const char *title);
    AliITSv11Hybrid(const char *name, const char *title);
    virtual       ~AliITSv11Hybrid() ;
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 110;} 
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   StepManager();
    virtual void   AddAlignableVolumes() const;
    virtual void   SetMinorVersion(Int_t v=2){ // Choose between existing minor versions
	fMinorVersion = v;}
    virtual void   SetThicknessDet1(Float_t v=200.){ 
	 // Set detector thickness in layer 1
	 fDet1 = v;}
    virtual void   SetThicknessDet2(Float_t v=200.){ 
	 // Set detector thickness in layer 2
	 fDet2 = v;}
    virtual void   SetThicknessChip1(Float_t v=150.){ 
	 // Set chip thickness in layer 1
	 fChip1 = v;}	 	 
    virtual void   SetThicknessChip2(Float_t v=150.){ 
	 // Set chip thickness in layer 2
	 fChip2 = v;}
    virtual void   SetRails(Int_t v=0){ 
	 // Set flag for rails
	 fRails = v;}	 
    virtual void   SetCoolingFluid(Int_t v=1){
	 // Set flag for cooling fluid
	 fFluid = v;}
    virtual void SetDensityServicesByThickness(){// uses services density
	// calculation based on the thickness of the services.
	fByThick = kTRUE;}
    virtual void SetDensityServicesByMass(){// uses services density
	// calculation based on the Mass of the services.
	fByThick = kFALSE;}
    virtual Int_t GetMajorVersion() const {// return Major Version Number
	return fMajorVersion;}
    virtual Int_t GetMinorVersion() const {// return Major Version Number
	return fMinorVersion;}
    virtual Float_t GetThicknessDet1() const { 
	 // Get detector thickness in layer 1
	 return fDet1;}
    virtual Float_t GetThicknessDet2() const { 
	 // Get detector thickness in layer 2
	 return fDet2;}
    virtual Float_t GetThicknessChip1() const { 
	 // Get chip thickness in layer 1
	 return fChip1;}	 	 
    virtual Float_t GetThicknessChip2()const { 
	 // Get chip thickness in layer 2
	 return fChip2;}
    virtual Int_t GetRails() const {
	 // Get flag for rails
	 return fRails;}
    virtual Int_t GetCoolingFluid() const{ 
	 // Get flag for cooling fluid
	 return fFluid;}

 protected:
    void CreateOldGeometry();
    void SetT2Lmatrix(Int_t uid, Double_t yShift,
		      Bool_t yFlip, Bool_t yRot180=kFALSE) const; // Set T2L matrix in TGeoPNEntries
    void CreateSPDThermalShield(TGeoVolume *moth);
    TGeoVolumeAssembly *CreateSPDThermalShieldAssembly(const char *name,
		   Double_t innerA, Double_t innerB, Double_t innerRadius,
		   Double_t outerA, Double_t outerB, Double_t outerRadius,
		   Double_t halflength, Double_t thickness,
	           Double_t thicknessOmega, Double_t theta);
    void SPDThermalShape(Double_t a, Double_t b, Double_t r, Double_t d,
			 Double_t t, Double_t *x, Double_t *y);
    void SPDOmegaShape(Double_t ina, Double_t inb, Double_t inr,
		       Double_t oua, Double_t oub, Double_t our,
		       Double_t dou, Double_t d  , Double_t t  ,
		       Double_t *x, Double_t *y);
    void FillSPDXtruShape(Double_t a, Double_t  b, Double_t  r,
			  Double_t t, Double_t *x, Double_t *y);

 private:
    AliITSv11Hybrid(const AliITSv11Hybrid &source); // copy constructor
    AliITSv11Hybrid& operator=(const AliITSv11Hybrid &source); // assignment operator

    Bool_t fByThick;          // Flag to use services materials by thickness
                              // ture, or mass false.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    Float_t  fDet1;	      // thickness of detector in SPD layer 1
    Float_t  fDet2;	      // thickness of detector in SPD layer 2
    Float_t  fChip1;	      // thickness of chip in SPD layer 1   
    Float_t  fChip2;	      // thickness of chip in SPD layer 2   
    Int_t    fRails;          // flag to switch rails on (=1) and off (=0)
    Int_t    fFluid;          // flag to switch between water (=1) and freon (=0)
    Int_t fIDMother;          //! ITS Mother Volume id.

    AliITSInitGeometry fInitGeom;   //! Get access to decoding and AliITSgeom init functins
    AliITSv11GeometrySPD *fSPDgeom; //! SPD Geometry
    AliITSv11GeometrySDD *fSDDgeom; //! SDD Geometry
    AliITSv11GeometrySSD *fSSDgeom; //! SSD Geometry
    AliITSv11GeometrySupport *fSupgeom; //! Support Geometry

    ClassDef(AliITSv11Hybrid,0)                          
};
 
#endif
