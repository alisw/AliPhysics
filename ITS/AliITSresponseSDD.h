#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
*/

#include "TArrayF.h"
#include "AliITSresponse.h"

// response for SDD

class AliITSresponseSDD : public AliITSresponse {
 public:
  //
  // Configuration methods
  //
  
    AliITSresponseSDD();
    AliITSresponseSDD(const char *dataType);
    virtual ~AliITSresponseSDD();

// Implementation of virtual member functions declared in AliITSresponse 
    virtual void SetElectronics(Int_t p1=1) 
      {fElectronics=p1;  /* Electronics: Pascal (1) or OLA (2) */ }
    virtual Int_t Electronics() const {// Electronics: 1 = Pascal; 2 = OLA
	return fElectronics;}
    virtual void    SetMaxAdc(Float_t p1) {// Adc-count saturation value
	fMaxAdc=p1;}
    virtual Float_t MaxAdc() const {// Get maximum Adc-count value
	return fMaxAdc;}
    virtual void    SetChargeLoss(Float_t p1) {
	// Set Linear Charge Loss Steepness  // 0.01 for 20%
	fChargeLoss=p1;}
    Float_t ChargeLoss() const {// Get Charge Loss Coefficient
	return fChargeLoss;}
    virtual void    SetDynamicRange(Float_t p1) {// Set Dynamic Range
	fDynamicRange = p1;}
    virtual Float_t DynamicRange() const {// Get Dynamic Range
	return fDynamicRange;}
    // Diffusion coefficients
    virtual void    SetDiffCoeff(Float_t p1, Float_t p2) 
	{fDiffCoeff=p1;fDiffCoeff1=p2;}
    // Get diffusion coefficients
    virtual void DiffCoeff(Float_t&diff,Float_t&diff1) const 
	{diff = fDiffCoeff;diff1 = fDiffCoeff1;}
    virtual void    SetDriftSpeed(Float_t p1) {// Drift velocity
	fDriftSpeed=p1;}
    virtual Float_t DriftSpeed() const {// drift speed
	return fDriftSpeed;}
    virtual void SetParamOptions(const char *opt1,const char *opt2){
	// Parameters: "same" or read from "file" 
	fParam1=opt1; fParam2=opt2;}
    virtual void   ParamOptions(char *opt1,char *opt2) const {// options
	strcpy(opt1,fParam1.Data()); strcpy(opt2,fParam2.Data());}
    virtual void  SetNoiseParam(Float_t n, Float_t b){
	// Noise and baseline  // 10 for ALICE with beam test measurements 2001
	fNoise=n; fBaseline=b;}
    virtual void  GetNoiseParam(Float_t &n, Float_t &b) const {// get noise param
	n=fNoise; b=fBaseline;}
    virtual Bool_t Do10to8() const {// get 10 to 8 compression option
	return fBitComp;}
    void    SetZeroSupp (const char *opt) {
	// Zero-suppression option - could be 1D, 2D or non-ZS 
	fOption=opt;}
    const char *ZeroSuppOption() const {// Get zero-suppression option
	return fOption.Data();}
    virtual void   SetThresholds(Float_t  mv, Float_t /* b */){
       // Min value used in 2D - could be used as a threshold setting
	fMinVal = mv;}
    virtual void   Thresholds(Float_t &  mv, Float_t & /* b */) const 
      {mv = fMinVal;}
    virtual void  GiveCompressParam(Int_t *x) const;
    // Detector type response methods
    virtual void    SetNSigmaIntegration(Float_t p1) {
	// Set number of sigmas over which cluster disintegration is performed
	fNsigmas=p1;}
    virtual Float_t NSigmaIntegration() const {
	// Get number of sigmas over which cluster disintegration is performed
	return fNsigmas;}
    virtual void SetNLookUp(Int_t p1);
    // Get number of intervals in which the gaussian lookup table is divided
    virtual Int_t GausNLookUp() const {return fNcomps;}
    virtual Float_t GausLookUp(Int_t i) const  {
	if(i<0 || i>=fNcomps) return 0.;return fGaus->At(i);}
   
//Declaration of member functions peculiar to this class
    Int_t Convert8to10(Int_t signal) const; //undo 10 to 8 bit SDD compresion
    void  SetNoiseAfterElectronics(Float_t n=2.38){
	// Noise after electronics (ADC units)
	// 2.36 for ALICE from beam test measurements 2001
	fNoiseAfterEl=n;}
    Float_t  GetNoiseAfterElectronics() const {
	// Noise after electronics (ADC units)
	return fNoiseAfterEl;}
    void  SetJitterError(Float_t jitter=20) {
	// set Jitter error (20 um for ALICE from beam test measurements 2001)
	fJitterError=jitter;}
    Float_t  JitterError() const {// set Jitter error
	return fJitterError;}
    void  SetDo10to8(Bool_t bitcomp=kTRUE) {
	// set the option for 10 to 8 bit compression
	fBitComp = bitcomp;}
    // Compression parameters
    void  SetCompressParam(Int_t cp[8]); 
    void SetDeadChannels(Int_t nmodules=0, Int_t nchips=0, Int_t nchannels=0);
    Int_t GetDeadModules() const { return fDeadModules; }
    Int_t GetDeadChips() const { return fDeadChips; }
    Int_t GetDeadChannels() const { return fDeadChannels; }
    Float_t Gain(Int_t mod,Int_t chip,Int_t ch)const {return fGain[mod][chip][ch]; }
    // these functions should be move to AliITSsegmentationSDD
    const Int_t Modules() const{return fgkModules;}// Total number of SDD modules
    const Int_t Chips() const{return fgkChips;} // Number of chips/module
    const Int_t Channels() const { return fgkChannels;}// Number of channels/chip
    //********
    void    PrintGains();
    void    Print();

    // not implemented virtual methods (devlared in the mother class
    virtual  void   SetDetParam(Float_t *)
      {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Float_t *) const 
      {NotImplemented("GetDetParam");}
    virtual  void   SetNDetParam(Int_t /* n */)
      {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam() const
      {NotImplemented("NDetParam"); return 0;}
    virtual void    SetSigmaSpread(Float_t, Float_t) 
      {NotImplemented("SetSigmaSpread");}
    virtual void    SigmaSpread(Float_t & /* p1 */,Float_t & /* p2 */) const 
      {NotImplemented("SigmaSpread");}

 protected:
    // these statis const should be move to AliITSsegmentationSDD
    static const Int_t fgkModules = 520;     // Total number of SDD modules
    static const Int_t fgkChips = 4;        // Number of chips/module
    static const Int_t fgkChannels = 64;    // Number of channels/chip
    //*******
    static const Int_t fgkMaxAdcDefault; // default for fMaxAdc
    static const Float_t fgkDynamicRangeDefault; // default for fDynamicRange
    static const Float_t fgkfChargeLossDefault; // default for fChargeLoss
    static const Float_t fgkDiffCoeffDefault; // default for fDiffCoeff
    static const Float_t fgkDiffCoeff1Default; // default for fDiffCoeff1 
    static const Float_t fgkTemperatureDefault; // default for fT (Kelvin)
    static const TString fgkParam1Default; // default for fParam1
    static const TString fgkParam2Default; // default for fParam2
    static const Float_t fgkNoiseDefault; // default for fNoise
    static const Float_t fgkBaselineDefault; // default for fBaseline
    static const TString fgkOptionDefault; // default for fOption
    static const Float_t fgkMinValDefault; // default for fMinVal
    static const Float_t fgkDriftSpeedDefault; // default for fDriftSpeed
    static const Float_t fgkNsigmasDefault; //default for fNsigmas
    static const Int_t fgkNcompsDefault; //default for fNcomps
    //

    Int_t fDeadModules;                   // Total number of dead SDD modules
    Int_t fDeadChips;                     // Number of dead chips
    Int_t fDeadChannels;                  // Number of dead channels
    Float_t   fGain[fgkModules][fgkChips][fgkChannels];   // Array for channel gains
    Int_t     fCPar[8];        // Hardware compression parameters
    Float_t   fNoise;          // Noise
    Float_t   fBaseline;       // Baseline
    Float_t   fNoiseAfterEl;   // Noise after electronics
    Float_t   fJitterError;    // jitter error
    Float_t   fDynamicRange;   // Set Dynamic Range 
    Float_t   fChargeLoss;     // Set Linear Coefficient for Charge Loss 
    Float_t   fDriftSpeed;     // Drift velocity
    Int_t     fElectronics;    // Electronics
    Float_t   fMaxAdc;         // Adc saturation value
    Float_t   fDiffCoeff;      // Diffusion Coefficient (scaling the time)
    Float_t   fDiffCoeff1;     // Diffusion Coefficient (constant term)
    Float_t   fNsigmas;    // Number of sigmas over which charge disintegration
                               // is performed
    TArrayF   *fGaus;          // Gaussian lookup table for signal generation
    Int_t      fNcomps;        // Number of samplings along the gaussian
    Float_t      fMinVal;        // Min value used in 2D zero-suppression algo
    Bool_t     fBitComp;       // 10 to 8 bit compression option
    TString    fOption;        // Zero-suppresion option (1D, 2D or none)
    TString    fParam1;        // Read baselines from file option
    TString    fParam2;        // Read compression algo thresholds from file

 private:
    AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
    AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.


    ClassDef(AliITSresponseSDD,3) // SDD response 
    
    };
#endif
