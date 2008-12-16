#ifndef ALIITSVSSD03_H
#define ALIITSVSSD03_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// ITS step manager and geometry class for the ITS SSD test beam geometry //
// of june 2003.                                                          //
////////////////////////////////////////////////////////////////////////////
#include "AliITSInitGeometry.h"
#include "AliITS.h"

class AliITSvSSD03 : public AliITS{
 public:
    AliITSvSSD03(); // default constructor
    AliITSvSSD03(const char *title,Int_t geomnum=2003); // standard constructor
    virtual ~AliITSvSSD03(); // destructor
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Int_t  IsVersion() const {// returns the ITS version number 
                                      return (Int_t)kvSSD03;} 
    virtual void   Init();
    //virtual void   SetDefaults();
    virtual void   DrawModule() const;
    virtual void   StepManager(); 
    virtual void   SetMinorVersion(Int_t v=22){ // Choose between existing minor versions
        fMinorVersion = v;}
    // Replacement default simulation initilization.
    //virtual void SetDefaultSimulation();
    //
  private:
    void BuildGeometry2003();
    void CreateGeometry2003();
    void CreateMaterials2003();

 private:  
    AliITSvSSD03(const AliITSvSSD03 &source); // Copy constructor
    AliITSvSSD03& operator=(const AliITSvSSD03 &source); // = operator
    Int_t  fMajorVersion;     // Major version number == IsVersion
    Int_t  fMinorVersion;     // Minor version number 
    Int_t  fGeomNumber;       // Geometry version number (year)
    Int_t  fIDMother;         //! ITS Mother Volume id.
    AliITSInitGeometry fIgm;  //! AliITSInitGeometry object

    ClassDef(AliITSvSSD03,3) // Hits manager and geometry for SSD testbeam
};
#endif
