#ifndef ALIITSVPPRASYMM_H
#define ALIITSVPPRASYMM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set: ITS version 6    //
/////////////////////////////////////////////////////////
 
#include "AliITS.h"
 
class AliITSvPPRasymm : public AliITS {

 public:
    AliITSvPPRasymm();
    AliITSvPPRasymm(const char *name, const char *title);
    AliITSvPPRasymm(const AliITSvPPRasymm &source); // copy constructor
    AliITSvPPRasymm& operator=(const AliITSvPPRasymm &source); // assignment operator
    virtual       ~AliITSvPPRasymm() ;
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init(); 
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 8;} 
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
    virtual void   SetMinorVersion(Int_t v){ // Choose between existing minor versions
	fMinorVersion = 1;
	if(v==1) fMinorVersion = 1;
	else if(v==2) fMinorVersion = 2;
	else Warning("SetMinorVersion","Undefined Minor Version setting =1");}
    virtual Bool_t GetEUCLID(){return fEuclidOut;}// returns value Euclid flag.
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

 private:
    void InitAliITSgeom();

    // TString fEuclidGeomtery,fEuclidMaterial defined in AliModule.
    Bool_t fEuclidOut;        // Flag to write geometry in euclid format
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file

    ClassDef(AliITSvPPRasymm,1)  //Hits manager for set:ITS version 8 
                                 // PPR detailed Geometry asymmetric
};
 
#endif
