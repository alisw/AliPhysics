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
#include "AliITSInitGeometry.h"
#include "AliITS.h"

class AliITSvSPD02 : public AliITS{
 public:
    AliITSvSPD02(); // default constructor
    AliITSvSPD02(const char *title,Int_t geomnum=2002); // standard constructor
    virtual ~AliITSvSPD02(); // destructor
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return (Int_t)kvSPD02;} 
    virtual void   Init();
    //virtual void   SetDefaults();
    virtual void   DrawModule() const;
    virtual void   StepManager(); 
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
    //virtual void SetDefaultSimulation();
    //
  private:
    void BuildGeometry2002();
    void CreateGeometry2002();
    void CreateMaterials2002();

 private:  
    AliITSvSPD02(const AliITSvSPD02 &source); // Copy constructor
    AliITSvSPD02& operator=(const AliITSvSPD02 &source); // = operator
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number 
    Int_t  fGeomNumber;       // Geometry version number (year)
    Float_t  fDet1;           // thickness of detector in SPD layer 1
    Float_t  fDet2;           // thickness of detector in SPD layer 2
    Float_t  fChip1;          // thickness of chip in SPD layer 1
    Float_t  fChip2;          // thickness of chip in SPD layer 2 
    Int_t fIDMother;          //! ITS Mother Volume id.
    AliITSInitGeometry fIgm;//! Get access to decoding and AliITSgeom init functins

    ClassDef(AliITSvSPD02,5) // Hits manager and geometry for SPD testbeam
};
#endif
