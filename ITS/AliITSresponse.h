#ifndef ALIITSRESPONSE_H
#define ALIITSRESPONSE_H


#include <TObject.h>

#include "AliITSsegmentation.h"

class TF1;
class TString;
class AliITSgeom;

//----------------------------------------------
//
// ITS response virtual base class
//
class AliITSresponse :
public TObject {
 public:


    virtual ~AliITSresponse() {}
    //
    // Configuration methods
    //

    // Set Electronics
    virtual void    SetElectronics(Int_t p1) {}
    // Get Electronics
    virtual Int_t Electronics()  {return 0;}                       

    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Float_t p1) {}
    // Get maximum Adc-count value
    virtual Float_t MaxAdc()  {return 0.;}                       

    // Set maximum Adc-top value
    virtual void    SetDynamicRange(Float_t p1) {}
    // Get maximum Adc-top value
    virtual Float_t DynamicRange()  {return 0.0;}                       

    // Set Charge Loss Linear Coefficient
    virtual void    SetChargeLoss(Float_t p1) {}
    // Get Charge Loss Linear Coefficient
    virtual Float_t ChargeLoss()  {return 0.0;}                       

    // Diffusion coefficient
    virtual void    SetDiffCoeff(Float_t, Float_t) {}
    // Get diffusion coefficients
    virtual void    DiffCoeff(Float_t &,Float_t &) {}

    // Temperature
    virtual void    SetTemperature(Float_t) {}
    // Get temperature
    virtual Float_t Temperature() {return 0.;}  
    // Type of data - real or simulated
    virtual void    SetDataType(const char *data) {}
    // Get data type
    virtual const char  *DataType() const {return "";} 


 
    // Set parameters options: "same" or read from "file" or "SetInvalid" or ...
    virtual void   SetParamOptions(const char* opt1,const char* opt2) {}
    // Set noise parameters 
    virtual void   SetNoiseParam(Float_t, Float_t) {}
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) {}
    // Set detector parameters: gain, coupling ...
    virtual  void   SetDetParam(Float_t *) {}

    // Parameters options
    virtual void   ParamOptions(char *,char*) {}
    virtual Int_t  NDetParam() {return 0;}
    virtual void   GetDetParam(Float_t *) {} 
    virtual void   GetNoiseParam(Float_t&, Float_t&) {}

    // Zero-suppression option - could be 1D, 2D or non-ZeroSuppressed 
    virtual void   SetZeroSupp(const char* opt) {}
    // Get zero-suppression option
    virtual const char *ZeroSuppOption() const {return "";}
     // Set thresholds
    virtual void   SetThresholds(Float_t, Float_t) {}
    virtual void   Thresholds(Float_t &, Float_t &) {}
    // Set min val
    virtual void   SetMinVal(Int_t) {};
    virtual Int_t  MinVal() {return 0;};

    // Set filenames
    virtual void   SetFilenames(const char *f1,const char *f2,const char *f3) {}
    // Filenames
    virtual void   Filenames(char*,char*,char*) {}


    virtual Float_t DriftSpeed() {return 0.;}
    virtual Bool_t  OutputOption() {return kFALSE;}
    virtual Bool_t  Do10to8() {return kTRUE;}
    virtual void    GiveCompressParam(Int_t *x) {}

    //  
    // Detector type response methods
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetNSigmaIntegration(Float_t p1) {}
    // Get number of sigmas over which cluster disintegration is performed
    virtual Float_t NSigmaIntegration() {return 0.;}
    // Set number of bins for the gaussian lookup table
    virtual void    SetNLookUp(Int_t p1) {}
    // Get number of bins for the gaussian lookup table
    virtual Int_t GausNLookUp() {return 0;}
    // Get scaling factor for bin i-th from the gaussian lookup table
    virtual Float_t GausLookUp(Int_t) {return 0.;}
    // Set sigmas of the charge spread function
    virtual void    SetSigmaSpread(Float_t p1, Float_t p2) {}
    // Get sigmas for the charge spread 
    virtual void    SigmaSpread(Float_t &s1, Float_t &s2) {}


    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss) {return 0.;}
    // Charge disintegration 
    virtual Float_t IntXZ(AliITSsegmentation *) {return 0.;}

    ClassDef(AliITSresponse,1) // Detector type response virtual base class 
};
#endif







