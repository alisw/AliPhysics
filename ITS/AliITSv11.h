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
 
class AliITSv11 : public AliITS {

 public:
    AliITSv11();
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
    virtual void   SetWriteDet(Bool_t det=kTRUE){ // set .det write
	                                         fGeomDetOut = det;}
    virtual void   SetWriteDet(const char *f){ // set write file
	                         strncpy(fWrite,f,60);fGeomDetOut = kTRUE;}
    virtual void   SetReadDet(Bool_t det=kTRUE){ //set .det read
	                                        fGeomDetIn = det;}
    virtual void   SetReadDet(const char *f){ // set read file
	                           strncpy(fRead,f,60);fGeomDetIn = kTRUE;}
    virtual void   SetEUCLID(Bool_t euclid=kTRUE){ // set write Euclid file
	                                          fEuclidOut = euclid;}
    virtual void   SetEUCLIDFileName(const char *f){ // set write file
	                     fEuclidGeometry=f;fEuclidOut = kTRUE;}
    virtual void   SetMinorVersion(Int_t v=00){ // Set minor version
	fMinorVersion = v;}
    virtual void   SetThicknessDet1(Float_t v=200.){ 
	 // Set detector thickness in layer 1
	 fDet1 = v;}
    virtual void   SetThicknessDet2(Float_t v=200.){ 
	 // Set detector thickness in layer 2
	 fDet2 = v;}
    virtual void   SetThicknessChip1(Float_t v=300.){ 
	 // Set chip thickness in layer 1
	 fChip1 = v;}	 	 
    virtual void   SetThicknessChip2(Float_t v=200.){ 
	 // Set chip thickness in layer 2
	 fChip2 = v;}
    virtual void   SetRails(Int_t v=1){ 
	 // Set flag for rails
	 fRails = v;}	 
    virtual void   SetCoolingFluid(Int_t v=1){ 
	 // Set flag for cooling fluid
	 fFluid = v;}	 	 
    virtual Bool_t GetEUCLID(){return fEuclidOut;}// returns value Euclid flag
    virtual const char  *GetEULIIDFileName() const{ // return .euc file name
	                               return fEuclidGeometry.Data();}
    virtual Bool_t GetWriteDet() { // returns value GeomDetOut flag.
	                          return fGeomDetOut;}
    virtual Bool_t GetReadDet() { // returns value GeomDetIn flag.
	                         return fGeomDetIn;}
    virtual char  *GetReadDetFileName(){ // return .det read file name
	if(fRead[0]!='\0') return fRead; else return fEuclidGeomDet;}
    virtual char  *GetWriteDetFileName(){ // return .det write file name
	if(fWrite[0]!='\0') return fWrite; else return fEuclidGeomDet;}
    virtual Int_t GetMajorVersion(){// return Major Version Number
	return fMajorVersion;}
    virtual Int_t GetMinorVersion(){// return Major Version Number
	return fMinorVersion;}
    virtual Float_t GetThicknessDet1(){ 
	 // Get detector thickness in layer 1
	 return fDet1;}
    virtual Float_t GetThicknessDet2(){ 
	 // Get detector thickness in layer 2
	 return fDet2;}
    virtual Float_t GetThicknessChip1(){ 
	 // Get chip thickness in layer 1
	 return fChip1;}	 	 
    virtual Float_t GetThicknessChip2(){ 
	 // Get chip thickness in layer 2
	 return fChip2;}
    virtual Int_t GetRails(){ 
	 // Get flag for rails
	 return fRails;}	 
    virtual Int_t GetCoolingFluid(){ 
	 // Get flag for cooling fluid
	 return fFluid;}	 	 	 
	 	 
 private:
    void InitAliITSgeom();
    void SPDCone(TGeoVolume *Moth);
    void SPDThermalSheald(TGeoVolume *Moth);
    void SDDCone(TGeoVolume *Moth);
    void SSDCone(TGeoVolume *Moth);
    Double_t RmaxFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z);
    Double_t RminFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z);
    Double_t RFrom2Points(Double_t *p,Double_t *Z,Int_t i1,Int_t i2,Double_t z);
    Double_t Zfrom2MinPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r);
    Double_t Zfrom2MaxPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r);
    Double_t Zfrom2Points(Double_t *Z,Double_t *p,Int_t i1,Int_t i2,Double_t r);
    Double_t RmaxFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,Double_t th=0.0);
    Double_t RmaxFromZpCone(Double_t *Rmax,Double_t *Z,Double_t tc,Double_t z,Double_t th=0.0);
    Double_t RminFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,Double_t th=0.0);
    Double_t RminFromZpCone(Double_t *Rmin,Double_t *Z,Double_t tc,Double_t z,Double_t th=0.0);
    Double_t ZFromRmaxpCone(TGeoPcon *p,Double_t tc,Double_t r,Double_t th=0.0);
    Double_t ZFromRmaxpCone(Double_t *GetRmax,Double_t *GetZ,Double_t tc,Double_t r,Double_t th=0.0);
    Double_t ZFromRminpCone(TGeoPcon *p,Double_t tc,Double_t r,Double_t th=0.0);
    void RadiusOfCurvature(Double_t rc,Double_t theta0,Double_t z0,
			   Double_t r0,Double_t theta1,Double_t &z1,
			   Double_t &r1);

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t fEuclidOut;        // Flag to write geometry in euclid format
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file
    Float_t  fDet1;	      // thickness of detector in SPD layer 1
    Float_t  fDet2;	      // thickness of detector in SPD layer 2
    Float_t  fChip1;	      // thickness of chip in SPD layer 1   
    Float_t  fChip2;	      // thickness of chip in SPD layer 2   
    Int_t    fRails;          // switch rails on (=1) and off (=0)
    Int_t    fFluid;          // switch between water(=1) and freon(=0)

    ClassDef(AliITSv11,1)  //Hits manager for set:ITS version 11
};
 
#endif
