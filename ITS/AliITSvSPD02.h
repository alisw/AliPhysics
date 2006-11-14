#ifndef ALIITSVSPD02_H
#define ALIITSVSPD02_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
  $Id$
*/
// ITS step manager and geometry class for the ITS SPD test beam geometry
// of Summer 2002.
///////////////////////////////////////////////////////////////////////
// Step manager and 
// geometry class
// for the ITS 
// SPD test beam
// geometry of summer 2002
// 
///////////////////////////////////////////////////////////////////////

#include "AliITS.h"

class AliITSvSPD02 : public AliITS{
 public:
    AliITSvSPD02(); // default constructor
    AliITSvSPD02(const char *title,Int_t geomnum=2002); // standard constructor
    virtual ~AliITSvSPD02(); // destructor
    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  DecodeDetector(Int_t id,Int_t cpy,Int_t &lay,Int_t &lad,
				  Int_t &det)const;
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
    void BuildGeometry2002();
    void CreateGeometry2002();
    void CreateMaterials2002();

 private:  
    AliITSvSPD02(const AliITSvSPD02 &source); // Copy constructor
    AliITSvSPD02& operator=(const AliITSvSPD02 &source); // = operator
    void InitAliITSgeom();
    Bool_t fGeomDetOut;       // Flag to write .det file out
    Bool_t fGeomDetIn;        // Flag to read .det file or directly from Geat.
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number 
    Int_t  fGeomNumber;       // Geometry version number (year)
    char   fEuclidGeomDet[60];// file where detector transormation are define.
    char   fRead[60];         //! file name to read .det file
    char   fWrite[60];        //! file name to write .det file 
    Float_t  fDet1;           // thickness of detector in SPD layer 1
    Float_t  fDet2;           // thickness of detector in SPD layer 2
    Float_t  fChip1;          // thickness of chip in SPD layer 1
    Float_t  fChip2;          // thickness of chip in SPD layer 2 
    Int_t fIDMother;          //! ITS Mother Volume id.

    ClassDef(AliITSvSPD02,2) // Hits manager and geometry for SPD testbeam
};
#endif
