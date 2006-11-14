#ifndef ALIITSV11_H
#define ALIITSV11_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

 
#include "AliITS.h"
class AliITSv11GeometrySPD;
class AliITSv11GeometrySDD;
class AliITSv11GeometrySupport;

class AliITSv11 : public AliITS {

 public:
    AliITSv11();
    AliITSv11(const char *name, const char *title);
    AliITSv11(Int_t debugITS, Int_t debugSPD = 0, Int_t debugSDD = 0,
	     Int_t debugSSD = 0, Int_t debugSUP = 0);
    virtual       ~AliITSv11() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   DrawModule() const;
    virtual void   StepManager();
    //virtual AliITSv11GeometrySPD*     GetSPDGeometry(){return fSPDgeom;}
    virtual AliITSv11GeometrySDD*    GetSDDGeometry(){return fSDDgeom;}
    //virtual AliITSv11GeometrySupport* GetSupGeometry(){return fSupgeom;}
    virtual Int_t  IsVersion() const { return fMajorVersion;}  // ITS version number 
    virtual Int_t  GetMajorVersion() const {// return Major Version Number
                    return fMajorVersion;}
    virtual Int_t  GetMinorVersion() const {// return Major Version Number
                    return fMinorVersion;}


    virtual void   SetWriteDet(Bool_t det=kTRUE){ // set .det write
	                                         fGeomDetOut = det;}
    virtual void   SetWriteDet(const char *f){ // set write file
	                             strncpy(fWrite,f,60);fGeomDetOut = kTRUE;}
    virtual void   SetReadDet(Bool_t det=kTRUE){ //set .det read
	                                        fGeomDetIn = det;}
    virtual void   SetReadDet(const char *f){ // set read file
	                               strncpy(fRead,f,60);fGeomDetIn = kTRUE;}
    virtual void   SetEUCLIDFileName(const char *f){ // set write file
	                     fEuclidGeometry=f;fEuclidOut = kTRUE;}
    virtual void   SetMinorVersion(Int_t v=0){ // Choose between existing minor versions
                   fMinorVersion = v;}

    virtual void SetDensityServicesByThickness(){// uses services density
	// calculation based on the thickness of the services.
	fByThick = kTRUE;}
    virtual void SetDensityServicesByMass(){// uses services density
	// calculation based on the Mass of the services.
	fByThick = kFALSE;}
    virtual Bool_t GetEUCLID() const {return fEuclidOut;}// returns value Euclid flag.
    virtual const char  *GetEULIIDFileName() const{ // return .euc file name
	                               return fEuclidGeometry.Data();}
    virtual Bool_t GetWriteDet() const { // returns value GeomDetOut flag.
	                          return fGeomDetOut;}
    virtual Bool_t GetReadDet() const { // returns value GeomDetIn flag.
	                         return fGeomDetIn;}
    virtual char  *GetReadDetFileName(){ // return .det read file name
	          if(fRead[0]!='\0') return fRead; else return fEuclidGeomDet;}
    virtual char  *GetWriteDetFileName(){ // return .det write file name
	        if(fWrite[0]!='\0') return fWrite; else return fEuclidGeomDet;}


 private:
    AliITSv11(const AliITSv11 &source);            // copy constructor
    AliITSv11& operator=(const AliITSv11 &source); // assignment operator
    void InitAliITSgeom();

    Bool_t   fGeomDetOut;       // Flag to write .det file out
    Bool_t   fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Bool_t   fByThick;          // Flag to use services materials by thickness
                                // ture, or mass false.
    Int_t    fMajorVersion;     // Major version number == IsVersion
    Int_t    fMinorVersion;     // Minor version number
    char     fEuclidGeomDet[60];// file where detector transormation are define.
    char     fRead[60];         //! file name to read .det file
    char     fWrite[60];        //! file name to write .det file


    //AliITSv11GeometrySPD *fSPDgeom;      //SPD Geometry
    AliITSv11GeometrySDD *fSDDgeom;      //! SDD Geometry
    //AliITSv11GeometrySupport /fSupgeom;  //Support Geometry

    ClassDef(AliITSv11,1)  // ITS version 11 
};
 
#endif
