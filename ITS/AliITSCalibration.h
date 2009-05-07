#ifndef ALIITSCALIBRATION_H
#define ALIITSCALIBRATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////
// Base ITS calibration class               //
//////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include "AliLog.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSDriftSpeedArraySDD.h"

class AliITSsegmentation;
class TF1;
class AliITSgeom;


/////////////////////////////////////////////
//                                         //
// ITS calibration virtual base class      //
/////////////////////////////////////////////
class AliITSCalibration : public TObject {
 public:
    // Default Constructor
    AliITSCalibration();

    // Destructor.
    virtual ~AliITSCalibration() {;}

    // Check for dead modules anche chips
    // Return 1 if the module/chip is dead, 0 if it is ok
    virtual Bool_t IsBad() const {AliError("This method must be implemented in a derived class"); return kFALSE;}
    virtual Bool_t IsChipBad(Int_t) const {AliError("This method must be implemented in a derived class"); return kFALSE;}
    //
    // Configuration methods
    //
    // Temperature in [degree K]
    virtual void    SetTemperature(Double_t t=300.0) {fT = t;}
    // Get temperature [degree K]
    virtual Double_t Temperature() const {return fT;}
 
    // Get data type
    virtual const char  *DataType() const {return fDataType.Data();}
    // Type of data - real or simulated
    virtual void    SetDataType(const char *data="simulated") {fDataType=data;}
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) = 0;
    // Set detector parameters: gain, coupling ...
    virtual  void   SetDetParam(Double_t *) = 0;

    virtual Int_t  NDetParam() const = 0;
    virtual void   GetDetParam(Double_t *) const = 0;
    virtual void   SetMapA(Int_t, AliITSCorrMapSDD*) {AliError("This method must be implemented in a derived class");}
    virtual void   SetMapT(Int_t, AliITSCorrMapSDD*) {AliError("This method must be implemented in a derived class");}
    virtual void   SetDriftSpeed(Int_t, AliITSDriftSpeedArraySDD*) {AliError("This method must be implemented in a derived class");}
    // Set sigmas of the charge spread function
    virtual void    SetSigmaSpread(Double_t, Double_t) = 0;
    // Get sigmas for the charge spread
    virtual void    SigmaSpread(Double_t &,Double_t &) const = 0;

    // Needed for SSD bad modules retrieval in the tracker
    virtual void    SetModuleIndex(Int_t /*modId*/) {};

    // Prints out the content of this class in ASCII format.
    virtual void Print(ostream *os) const;
    // Reads in the content of this class in the format of Print
    virtual void Read(istream *is);
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}


 protected:
    AliITSCalibration(const AliITSCalibration &ob); // copy constructor
    AliITSCalibration& operator=(const AliITSCalibration& source); // ass.
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}

    TString  fDataType;   // data type - real or simulated
    
    Float_t fT;   // The temperature of the Si in Degree K.
                 // deleted here but in AliITSDetTypeSim and AliITSDetTypeRec

    ClassDef(AliITSCalibration,3) // Detector type response virtual base class 
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSCalibration &source);
istream& operator>>(istream &os,AliITSCalibration &source);
#endif
