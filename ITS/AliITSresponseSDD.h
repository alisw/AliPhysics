#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H

#include "AliITSresponse.h"

// response for SDD

class AliITSresponseSDD :
  public AliITSresponse {
public:
  //
  // Configuration methods
  //
  
  AliITSresponseSDD();
  virtual ~AliITSresponseSDD() { 
    // destructor
  }
  AliITSresponseSDD(const AliITSresponseSDD &source); // copy constructor
  AliITSresponseSDD& operator=(const AliITSresponseSDD &source); // ass. op.
  
  virtual void    SetMaxAdc(Float_t p1=1023) {
    // Adc-count saturation value
    fMaxAdc=p1;
  }
  virtual Float_t MaxAdc()  {
    // Get maximum Adc-count value
    return fMaxAdc;
  }                       
  
  virtual void    SetMagicValue(Float_t p1=450.) {
    // Set maximum Adc-magic value
    fTopValue=p1;
  }
  virtual Float_t MagicValue()  {
    // Get maximum Adc-magic value
    return fTopValue;
  }                       
  
  virtual void    SetDiffCoeff(Float_t p1=5.) {
    // Diffusion coefficient
    fDiffCoeff=p1;
  }
  virtual Float_t DiffCoeff() {
    // Get diffusion coefficient
    return fDiffCoeff;
  } 
  
  virtual void    SetQref(Float_t p1=120.) {
    // Coulomb repulsion
    fQref=p1;
  }
  virtual Float_t Qref() {
    // qref
    return fQref;
  }
  
  virtual void    SetDriftSpeed(Float_t p1=7.5) {
    // Drift velocity
    fDriftSpeed=p1;
  }
  virtual Float_t DriftSpeed() {
    // drift speed
    return fDriftSpeed;
  } 
  
  virtual void    SetTemperature(Float_t p1=23.) {
    // Temperature
    fTemperature=p1;
  }
  virtual Float_t Temperature() {
    // Get temperature
    return fTemperature;
  } 
  
  virtual void    SetDataType(char *data="simulated") {
    // Type of data - real or simulated
    fDataType=data;
  }
  virtual char  *DataType() {
    // Get data type
    return fDataType;
  } 
  
  virtual void SetParamOptions(Option_t *opt1="same",Option_t *opt2="same"){
    // Parameters: "same" or read from "file" 
    fParam1=opt1; fParam2=opt2;
  }
  virtual void   ParamOptions(Option_t *&opt1,Option_t *&opt2) {
    // options
    opt1=fParam1; opt2=fParam2;
  }
  
  virtual  void  SetNoiseParam(Float_t n=3., Float_t b=20.){
    // Noise and baseline
    fNoise=n; fBaseline=b;
  }   
  virtual  void  GetNoiseParam(Float_t &n, Float_t &b) {
    // get noise param
    n=fNoise; b=fBaseline;
  }   
  
  virtual void    SetZeroSupp(Option_t *opt="2D") {
    // Zero-suppression option - could be 1D, 2D or non-ZS 
    fOption=opt;
  }
  virtual Option_t *ZeroSuppOption() {
    // Get zero-suppression option
    return fOption;
  }
  virtual  void  SetMinVal(Int_t mv=4) {
    // Min value used in 2D - could be used as a threshold setting
    fMinVal = mv;
  }
  virtual Int_t  MinVal() {
    // min val
    return fMinVal;
  }
  
  virtual void   SetFilenames(char *f1=0,char *f2=0, char *f3=0) {
    // Set filenames - input, output, parameters ....
    fFileName1=f1; fFileName2=f2; fFileName3=f3;
  }
  virtual void   Filenames(const char*input,const char*baseline,const char*param) {
    // Filenames
    input=fFileName1; baseline=fFileName2; param=fFileName3;
  }     
  
  
  virtual  void  SetOutputOption(Bool_t write=kFALSE) {
    // set output option
    fWrite = write;
  }
  Bool_t OutputOption()  {
    // output option
    return fWrite;
  }
  // 
  // Compression parameters
  virtual  void  SetCompressParam(Int_t cp[8]); 
  void  GiveCompressParam(Int_t *x);
  
  //  
  // Detector type response methods
  virtual void    SetNSigmaIntegration(Float_t p1) {
    // Set number of sigmas over which cluster disintegration is performed
  }
  virtual Float_t NSigmaIntegration() {
    // Get number of sigmas over which cluster disintegration is performed
    return 0.;
  }
  virtual void    SetSigmaSpread(Float_t p1, Float_t p2) {
    // Set sigmas of the charge spread function
  }
  virtual void    SigmaSpread(Float_t &s1, Float_t &s2) {
    // Get sigmas for the charge spread 
  }
  
  virtual Float_t IntPH(Float_t eloss) {
    // Pulse height from scored quantity (eloss)
    return 0.;
  }
  virtual Float_t IntXZ(AliITSsegmentation *) {
    // Charge disintegration 
    return 0.;
  }
  
  
protected:
  
  Int_t     fCPar[8];        // Hardware compression parameters
  //Int_t     fNDetPar;        // Number of detector param 
  //Float_t   fDetPar[fNDetPar];  
  
  Float_t   fNoise;          // Noise
  Float_t   fBaseline;       // Baseline
  Float_t   fTopValue;       // still unclear to me 
  Float_t   fTemperature;    // Temperature 
  Float_t   fDriftSpeed;     // Drift velocity
  
  Float_t    fMaxAdc;        // Adc saturation value
  Float_t    fDiffCoeff;     // Diffusion Coefficient
  Float_t    fQref;          // Coulomb repulsion
  
  Int_t      fZeroSuppFlag;  // Zero-suppression flag
  Int_t      fMinVal;        // Min value used in 2D zero-suppression algo
  
  Bool_t     fWrite;         // Write option for the compression algorithms
  Option_t   *fOption;       //! 
                             // Zero-suppresion option (1D, 2D or none)
  Option_t   *fParam1;       //! 
                             //Read baselines from file option
  Option_t   *fParam2;       //! 
                             //Read compression algo thresholds from file 
  
  char*         fDataType;         //!
				   // input keys : run, module #
  char*         fFileName1;        //!
                                   // input keys : run, module #
  char*         fFileName2;        //!
                                   // baseline & noise val or output coded                                        // signal or monitored bgr.
  char*         fFileName3;        //!
                                   // param values or output coded signal 
  
  ClassDef(AliITSresponseSDD,1) // SDD response 
    
    };
#endif







