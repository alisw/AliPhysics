#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H

#include "TArrayF.h"
#include <TString.h>
#include <iostream.h>
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

  virtual void SetElectronics(Int_t p1=1) {
    // Electronics: Pascal or OLA
    fElectronics=p1;
  }
  
  virtual Int_t Electronics() {
    // Electronics: 1 = Pascal; 2 = OLA
    return fElectronics;
  }
  
  virtual void    SetMaxAdc(Float_t p1=1024.) {
    // Adc-count saturation value
    fMaxAdc=p1;
  }
  virtual Float_t MaxAdc()  {
    // Get maximum Adc-count value
    return fMaxAdc;
  }                       
  
  virtual void    SetChargeLoss(Float_t p1=0.01) {
    // Set Linear Charge Loss Steepness  
    fChargeLoss=p1;
  }
  virtual Float_t ChargeLoss()  {
    // Get Charge Loss Coefficient
    return fChargeLoss;
  }                       
  
  virtual void    SetDynamicRange(Float_t p1=132.) {
    // Set Dynamic Range
    fDynamicRange=p1;
  }
  virtual Float_t DynamicRange()  {
    // Get Dynamic Range
    return fDynamicRange;
  }                       
  
  virtual void    SetDiffCoeff(Float_t p1=3.23,Float_t p2=30.) {
    // Diffusion coefficients
    fDiffCoeff=p1;
    fDiffCoeff1=p2;
  }
  virtual void DiffCoeff(Float_t&diff,Float_t&diff1) {
    // Get diffusion coefficients
    diff = fDiffCoeff;
    diff1 = fDiffCoeff1;
  } 
  
  virtual void    SetDriftSpeed(Float_t p1=7.3) {
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
  
  virtual void    SetDataType(const char *data="simulated") {
    // Type of data - real or simulated
    fDataType=data;
  }
  virtual const char  *DataType() const {
    // Get data type
    return fDataType.Data();
  } 
  
  virtual void SetParamOptions(const char *opt1="same",const char *opt2="same"){
    // Parameters: "same" or read from "file" 
    fParam1=opt1; fParam2=opt2;
  }
  virtual void   ParamOptions(char *opt1,char *opt2) {
    // options
    strcpy(opt1,fParam1.Data()); strcpy(opt2,fParam2.Data());
  }
  
  virtual  void  SetNoiseParam(Float_t n=8.3, Float_t b=20.){
    // Noise and baseline
    fNoise=n; fBaseline=b;
  }   
  virtual  void  SetNoiseAfterElectronics(Float_t n=1.6){
    // Noise after electronics (ADC units)
    fNoiseAfterEl=n;
  }   
  virtual  void  GetNoiseParam(Float_t &n, Float_t &b) {
    // get noise param
    n=fNoise; b=fBaseline;
  }  
  virtual  Float_t  GetNoiseAfterElectronics(){
    // Noise after electronics (ADC units)
    return fNoiseAfterEl;
  }   

  virtual  void  SetDo10to8(Bool_t bitcomp=kTRUE) {
    // set the option for 10 to 8 bit compression
    fBitComp = bitcomp;
  }

  Bool_t Do10to8() {
    // get 10 to 8 compression option
    return fBitComp;
  }   
  
  virtual void    SetZeroSupp (const char *opt="1D") {
    // Zero-suppression option - could be 1D, 2D or non-ZS 
    fOption=opt;
  }
  virtual const char *ZeroSuppOption() const {
    // Get zero-suppression option
    return fOption.Data();
  }
  virtual  void  SetMinVal(Int_t mv=4) {
    // Min value used in 2D - could be used as a threshold setting
    fMinVal = mv;
  }
  virtual Int_t  MinVal() {
    // min val
    return fMinVal;
  }
  
  virtual void   SetFilenames(const char *f1="",const char *f2="",const char *f3="") {
    // Set filenames - input, output, parameters ....
    fFileName1=f1; fFileName2=f2; fFileName3=f3;
  }
  virtual void   Filenames(char *input,char *baseline,char *param) {
    // Filenames
   strcpy(input,fFileName1.Data());  strcpy(baseline,fFileName2.Data());  
   strcpy(param,fFileName3.Data());
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
  virtual void    SetNSigmaIntegration(Float_t p1=3.) {
    // Set number of sigmas over which cluster disintegration is performed
    fNsigmas=p1;
  }
  virtual Float_t NSigmaIntegration() {
    // Get number of sigmas over which cluster disintegration is performed
    return fNsigmas;
  }
  virtual void SetNLookUp(Int_t p1=121) {
    // Set number of sigmas over which cluster disintegration is performed
    fNcomps=p1;
    fGaus = new TArrayF(fNcomps+1);
    for(Int_t i=0; i<=fNcomps; i++) {
      Float_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
      (*fGaus)[i] = exp(-((x*x)/2));
 //     cout << "fGaus[" << i << "]: " << fGaus->At(i) << endl;
    }
  }
  virtual Int_t GausNLookUp() { 
    // Get number of intervals in which the gaussian lookup table is divided
    return fNcomps;
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
  virtual Float_t GausLookUp(Int_t i) {
    if(i<0 || i>=fNcomps) return 0.;
    return fGaus->At(i);
  }
  virtual void    Print();
    
protected:
  
  Int_t     fCPar[8];        // Hardware compression parameters
  Float_t   fNoise;          // Noise
  Float_t   fBaseline;       // Baseline
  Float_t   fNoiseAfterEl;   // Noise after electronics
  Float_t   fDynamicRange;   // Set Dynamic Range 
  Float_t   fChargeLoss;     // Set Linear Coefficient for Charge Loss 
  Float_t   fTemperature;    // Temperature 
  Float_t   fDriftSpeed;     // Drift velocity
  Int_t     fElectronics;    // Electronics
  
  Float_t    fMaxAdc;        // Adc saturation value
  Float_t    fDiffCoeff;     // Diffusion Coefficient (scaling the time)
  Float_t    fDiffCoeff1;    // Diffusion Coefficient (constant term)
  Float_t    fNsigmas;       // Number of sigmas over which charge disintegration 
                             // is performed 
  TArrayF   *fGaus;          // Gaussian lookup table for signal generation
  Int_t      fNcomps;        // Number of samplings along the gaussian
  
  Int_t      fZeroSuppFlag;  // Zero-suppression flag
  Int_t      fMinVal;        // Min value used in 2D zero-suppression algo
  
  Bool_t     fWrite;         // Write option for the compression algorithms
  Bool_t     fBitComp;       // 10 to 8 bit compression option

  TString    fOption;        // Zero-suppresion option (1D, 2D or none)
  TString    fParam1;        // Read baselines from file option
  TString    fParam2;        // Read compression algo thresholds from file 
  
  TString         fDataType;         // data type - real or simulated
  TString         fFileName1;        // input keys : run, module #
  TString         fFileName2;        // baseline & noise val or output coded                                                 // signal or monitored bgr.
  TString         fFileName3;        // param values or output coded signal 
  
  ClassDef(AliITSresponseSDD,1) // SDD response 
    
    };
#endif







