#ifndef ALIITSRESPONSE_H
#define ALIITSRESPONSE_H


#include <TObject.h>
#include <TF1.h>
#include "AliITSsegmentation.h"

class AliITSgeom;

//----------------------------------------------
//
// ITS response virtual base class
//
class AliITSresponse :
public TObject {
 public:
    //
    // Configuration methods
    //

    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Float_t p1) {}
    // Get maximum Adc-count value
    virtual Float_t MaxAdc()  {return 0.;}                       

    // Set maximum Adc-magic value
    virtual void    SetMagicValue(Float_t p1) {}
    // Get maximum Adc-magic value
    virtual Float_t MagicValue()  {return 0.0;}                       

    // Diffusion coefficient
    virtual void    SetDiffCoeff(Float_t)                     =0;
    // Get diffusion coefficient
    virtual Float_t DiffCoeff()                               =0;
    virtual Float_t Qref() {return 0.;}

    // Temperature
    virtual void    SetTemperature(Float_t) {}
    // Get temperature
    virtual Float_t Temperature() {return 0.;}  
    // Type of data - real or simulated
    virtual void    SetDataType(char *data)        =0;
    // Get data type
    virtual const char  *DataType()                =0; 


 
    // parameters: "same" or read from "file" or "SetInvalid" or ...
    virtual void    SetParamOptions(Option_t *opt1, Option_t *opt2) {}
    virtual  void   SetNoiseParam(Float_t, Float_t) {}
    // gain, coupling ...
    // Number of parameters to be set
    virtual  void   SetNDetParam(Int_t) {}
    virtual  void   SetDetParam(Float_t *) {}

    // Parameters options
    virtual void   ParamOptions(Option_t *&, Option_t *&) {}
    virtual Int_t  NDetParam() {return 0;}
    virtual void   GetDetParam(Float_t *) {} 
    virtual void   GetNoiseParam(Float_t&, Float_t&) {}

    // Zero-suppression option - could be 1D, 2D or non-ZS 
    virtual void   SetZeroSupp(Option_t *opt) {}
    // Get zero-suppression option
    virtual Option_t   *ZeroSuppOption() {return "";}
     // Set thresholds
    virtual void   SetThresholds(Float_t, Float_t) {}
    virtual void   Thresholds(Float_t &, Float_t &) {}
    // Set min val
    virtual void   SetMinVal(Int_t) {};
    virtual Int_t  MinVal() {return 0;};

    // Set filenames
    virtual void   SetFilenames(char *,char *,char *) {}
    // Filenames
    virtual void   Filenames(const char*&,const char*&,const char*&) {}


    virtual Float_t DriftSpeed() {return 0.;}
    virtual Bool_t  OutputOption() {return kFALSE;}
    virtual void    GiveCompressParam(Int_t *x) {}

    //  
    // Detector type response methods
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetNSigmaIntegration(Float_t p1) {}
    // Get number of sigmas over which cluster disintegration is performed
    virtual Float_t NSigmaIntegration() {return 0.;}
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







