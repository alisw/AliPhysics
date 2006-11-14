#ifndef ALIITSVPPRASYMMFMD_H
#define ALIITSVPPRASYMMFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 10   //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSvPPRasymmFMD : public AliITS {

 public:
    AliITSvPPRasymmFMD();
    AliITSvPPRasymmFMD(const char *name, const char *title);
    virtual       ~AliITSvPPRasymmFMD() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 10;} 
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   DrawModule() const;
    virtual void   StepManager();
    virtual void   AddAlignableVolumes() const;
    virtual void   SetWriteDet(Bool_t det=kFALSE){ // set .det write
	                                         fGeomDetOut = det;}
    virtual void   SetWriteDet(const char *f){ // set write file
	                             strncpy(fWrite,f,60);fGeomDetOut = kTRUE;}
    virtual void   SetReadDet(Bool_t det=kFALSE){ //set .det read
	                                        fGeomDetIn = det;}
    virtual void   SetReadDet(const char *f){ // set read file
	                               strncpy(fRead,f,60);fGeomDetIn = kTRUE;}
    virtual void   SetEUCLIDFileName(const char *f){ // set write file
	                     fEuclidGeometry=f; SetEUCLID();}
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

 private:
    AliITSvPPRasymmFMD(const AliITSvPPRasymmFMD &source); // copy constructor
    AliITSvPPRasymmFMD& operator=(const AliITSvPPRasymmFMD &source); // assignment operator
    void InitAliITSgeom();

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Bool_t fByThick;          // Flag to use services materials by thickness
                              // ture, or mass false.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file
    Float_t  fDet1;	      // thickness of detector in SPD layer 1
    Float_t  fDet2;	      // thickness of detector in SPD layer 2
    Float_t  fChip1;	      // thickness of chip in SPD layer 1   
    Float_t  fChip2;	      // thickness of chip in SPD layer 2   
    Int_t    fRails;          // flag to switch rails on (=1) and off (=0)
    Int_t    fFluid;          // flag to switch between water (=1) and freon (=0)
    Int_t fIDMother;          //! ITS Mother Volume id.

    ClassDef(AliITSvPPRasymmFMD,3)  //Hits manager for set:ITS version 8 
                                 // PPR detailed Geometry asymmetric
};
 
#endif
