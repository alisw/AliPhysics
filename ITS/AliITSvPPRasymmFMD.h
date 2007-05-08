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
    //set .det write
    virtual void  SetWriteDet(Bool_t det=kFALSE){fGeomDetOut=det;}
    // set write file
    virtual void   SetWriteDet(const Char_t *f){fWrite=f;fGeomDetOut=kTRUE;}
    //set .det read
    virtual void   SetReadDet(Bool_t det=kFALSE){fGeomDetIn = det;}
    // set read file
    virtual void   SetReadDet(const Char_t *f){fRead=f;fGeomDetIn=kTRUE;}
    // set write file
    virtual void   SetEUCLIDFileName(const Char_t *f){fEuclidGeometry=f;
                                                          SetEUCLID();}
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
    virtual const Char_t *GetEULIIDFileName() const{ // return .euc file name
        return fEuclidGeometry.Data();}
    // returns value GeomDetOut flag.
    virtual Bool_t GetWriteDet() const {return fGeomDetOut;}
    // returns value GeomDetIn flag.
    virtual Bool_t GetReadDet() const{return fGeomDetIn;} 
    virtual const Char_t *GetReadDetFileName()const{//return .det read file name
        if(fRead.IsNull()) return fRead.Data();
        else return GetEULIIDFileName();}
    virtual const Char_t *GetWriteDetFileName()const{//return .det write file name
        if(fWrite.IsNull()) return fWrite.Data();
        else return GetEULIIDFileName();}
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
    // returns Euclid file name with where Euclid geometry is kept.
    const TString & GetEuclidFile()const{return fEuclidGeomDet;}
    // Return Mother volume ID
    Int_t GetMotherID()const {return fIDMother;}
    // Return AliITSInitGeometry object
    const AliITSInitGeometry & GetGeomInit()const{return fIgm;}

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t  fGeomDetOut;    // Flag to write .det file out
    Bool_t  fGeomDetIn;     // Flag to read .det file or directly from Geat.
    Bool_t  fByThick;       // Flag to use services materials by thickness
                            // ture, or mass false.
    Int_t   fMajorVersion;  // Major version number == IsVersion
    Int_t   fMinorVersion;  // Minor version number
    TString fEuclidGeomDet; // file where detector transormation are define.
    TString fRead;          //! file name to read .det file
    TString fWrite;         //! file name to write .det file
    Float_t fDet1;	    // thickness of detector in SPD layer 1
    Float_t fDet2;          // thickness of detector in SPD layer 2
    Float_t fChip1;         // thickness of chip in SPD layer 1   
    Float_t fChip2;         // thickness of chip in SPD layer 2   
    Int_t   fRails;         // flag to switch rails on (=1) and off (=0)
    Int_t   fFluid;         // flag to switch between water (=1) and freon (=0)
    Int_t   fIDMother;      //! ITS Mother Volume id.
    AliITSInitGeometry fIgm;//! Get access to decoding and AliITSgeom init functins

    ClassDef(AliITSvPPRasymmFMD,4) //Hits manager for set:ITS version 10
                                   // PPR detailed Geometry asymmetric
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,const AliITSvPPRasymmFMD &s);
istream &operator>>(istream &is,AliITSvPPRasymmFMD &s);

#endif
