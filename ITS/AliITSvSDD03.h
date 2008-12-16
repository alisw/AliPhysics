#ifndef ALIITSVSDD03_H
#define ALIITSVSDD03_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
  $Id$
*/
// ITS step manager and geometry class for the ITS SDD test beam geometry
// of Summer 2003 and later.
// At present, the geometry and materials need to be checked against the
// proper geometry of the 2003 test beam. In addition, because the SSD
// used in the test beam are single sided, the second side needs be ignored.
// This can cause some problems with the SSD reconstruction code.
#include "AliITSInitGeometry.h"
#include "AliITS.h"

typedef enum {
  kNoTarg=0,kIron=4,kLead,kSteel,
  kCarbon=8,kAl,kBe,kTi,kSn,kCopper,kGe,kTungsten=19
} TargTyp_t;

class AliITSvSDD03 : public AliITS{
 public:
    AliITSvSDD03(); // default constructor
    AliITSvSDD03(const char *title,Int_t year); // standard constructor
    virtual ~AliITSvSDD03(); // destructor
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return (Int_t)kvSDD03;} 
    virtual void   Init(); 
    //virtual void   SetDefaults();
    virtual void   DrawModule() const;
    virtual void   StepManager(); 
    virtual void   SetTargMat(TargTyp_t tt=kNoTarg){ // set target material 
                                                 fTarg=tt;}    
    virtual void   SetTargThick(Float_t th=1){ // set target thickness 
                                                 fTargThick=th;}    
    virtual void   SetMinorVersion(Int_t v=22){ // Choose between existing minor versions
        fMinorVersion = v;}
    // Replacement default simulation initilization.
    //virtual void SetDefaultSimulation();
    TargTyp_t GetTargMat() const {return fTarg;}
    // Decodes the id and copy nuber to give the layer, ladder, and detector 
    // numbers . Returns the module number.
    //virtual Int_t DecodeDetector(Int_t id,Int_t cpy,Int_t &lay,
    //                             Int_t &lad,Int_t &det) const;
         //
 private:  
    AliITSvSDD03(const AliITSvSDD03 &source); // Copy constructor
    AliITSvSDD03& operator=(const AliITSvSDD03 &source); // = operator
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number 
    Int_t fIDMother;          //! ITS Mother Volume id.
    Int_t fYear;              // Year flag to select different geometries.
    TargTyp_t fTarg;          // Target material
    Float_t fTargThick;       // TargetThickness in mm
    AliITSInitGeometry fIgm;  //! Get Access to decoding an dAliITSgeom init functions

    ClassDef(AliITSvSDD03,4) // Hits manager and geometry for SDD testbeam
};
#endif
