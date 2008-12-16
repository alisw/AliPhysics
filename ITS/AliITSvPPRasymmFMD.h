#ifndef ALIITSVPPRASYMMFMD_H
#define ALIITSVPPRASYMMFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 10   //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
#include "AliITSInitGeometry.h"

class AliITSvPPRasymmFMD : public AliITS {

 public:
    AliITSvPPRasymmFMD();
    AliITSvPPRasymmFMD(const Char_t *title); // Standard Constructor
    // Extended Standard constructor
    AliITSvPPRasymmFMD(const char *name, const char *title);
    virtual       ~AliITSvPPRasymmFMD() ;
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 10;} 
    virtual void   Init(); 
    virtual void   SetDefaults();
    virtual void   DrawModule() const;
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
    //
    // Print class in ascii form to stream
    virtual void PrintAscii(ostream *os)const;
    // Read in class in ascii form from stream
    virtual void ReadAscii(istream *is);

 private:
    // copy constructor. Will not copy class, just gives warning.
    AliITSvPPRasymmFMD(const AliITSvPPRasymmFMD &source);
    // assignment operator. Will not copy class, just gives warning.
    AliITSvPPRasymmFMD& operator=(const AliITSvPPRasymmFMD &source);
    void InitAliITSgeom();
    Bool_t IsDensityServicesByThickness()const {return fByThick;}
    // Return Mother volume ID
    Int_t GetMotherID()const {return fIDMother;}
    // Return AliITSInitGeometry object
    const AliITSInitGeometry & GetGeomInit()const{return fIgm;}

    Bool_t  fByThick;       // Flag to use services materials by thickness
                            // ture, or mass false.
    Int_t   fMajorVersion;  // Major version number == IsVersion
    Int_t   fMinorVersion;  // Minor version number
    Float_t fDet1;	    // thickness of detector in SPD layer 1
    Float_t fDet2;          // thickness of detector in SPD layer 2
    Float_t fChip1;         // thickness of chip in SPD layer 1   
    Float_t fChip2;         // thickness of chip in SPD layer 2   
    Int_t   fRails;         // flag to switch rails on (=1) and off (=0)
    Int_t   fFluid;         // flag to switch between water (=1) and freon (=0)
    Int_t   fIDMother;      //! ITS Mother Volume id.
    AliITSInitGeometry fIgm;//! Get access to decoding and AliITSgeom init functins

    ClassDef(AliITSvPPRasymmFMD,5) //Hits manager for set:ITS version 10
                                   // PPR detailed Geometry asymmetric
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,const AliITSvPPRasymmFMD &s);
istream &operator>>(istream &is,AliITSvPPRasymmFMD &s);

#endif
