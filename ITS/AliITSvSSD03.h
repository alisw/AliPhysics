#ifndef ALIITSVSSD03_H
#define ALIITSVSSD03_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// ITS step manager and geometry class for the ITS SSD test beam geometry //
// of june 2003.                                                          //
////////////////////////////////////////////////////////////////////////////
#include "AliITS.h"

class AliITSvSSD03 : public AliITS{
 public:
    AliITSvSSD03(); // default constructor
    AliITSvSSD03(const char *title,Int_t geomnum=2003); // standard constructor
    virtual ~AliITSvSSD03(); // destructor
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return 1;} 
    virtual void   Init();
    virtual void   SetDefaults();
    virtual void   DrawModule() const;
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
    virtual void   SetMinorVersion(Int_t v=22){ // Choose between existing minor versions
        fMinorVersion = v;} 
    virtual void   SetThicknessDet1(Float_t v=300.){
         // Set detector thickness in layer 1
         fDet1 = v;}
    virtual void   SetThicknessDet2(Float_t v=300.){
         // Set detector thickness in layer 2
         fDet2 = v;}
    virtual void   SetThicknessChip1(Float_t v=300.){
         // Set chip thickness in layer 1
         fChip1 = v;}
    virtual void   SetThicknessChip2(Float_t v=300.){
         // Set chip thickness in layer 2
         fChip2 = v;}
    // Replacement default simulation initilization.
    virtual void SetDefaultSimulation();
    //
  private:
    void BuildGeometry2003();
    void CreateGeometry2003();
    void CreateMaterials2003();

 private:  
    AliITSvSSD03(const AliITSvSSD03 &source); // Copy constructor
    AliITSvSSD03& operator=(const AliITSvSSD03 &source); // = operator
    void InitAliITSgeom();
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number 
    Int_t  fGeomNumber;       // Geometry version number (year)
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file 
    Float_t  fDet1;           // thickness of detector in SSD layer 1
    Float_t  fDet2;           // thickness of detector in SSD layer 2
    Float_t  fChip1;          // thickness of chip in SSD layer 1
    Float_t  fChip2;          // thickness of chip in SSD layer 2 
    Int_t fIDMother;          //! ITS Mother Volume id.

    ClassDef(AliITSvSSD03,2) // Hits manager and geometry for SSD testbeam
};
#endif
