#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
*/

#include "TArrayF.h"
#include <TString.h>
#include <Riostream.h>
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

    void SetElectronics(Int_t p1=1) {// Electronics: Pascal (1) or OLA (2)
	fElectronics=p1;}
    Int_t Electronics() const {// Electronics: 1 = Pascal; 2 = OLA
	return fElectronics;}
    void    SetMaxAdc(Float_t p1=1024.) {// Adc-count saturation value
	fMaxAdc=p1;}
    Float_t MaxAdc() const {// Get maximum Adc-count value
	return fMaxAdc;}
    void    SetChargeLoss(Float_t p1=0.0) {
	// Set Linear Charge Loss Steepness  // 0.01 for 20%
	fChargeLoss=p1;}
    Float_t ChargeLoss() const {// Get Charge Loss Coefficient
	return fChargeLoss;}
    void    SetDynamicRange(Float_t p1=132.) {// Set Dynamic Range
	fDynamicRange=p1;}
    Float_t DynamicRange() const {// Get Dynamic Range
	return fDynamicRange;}
    void    SetDiffCoeff(Float_t p1=3.23,Float_t p2=30.) {
	// Diffusion coefficients
	fDiffCoeff=p1;fDiffCoeff1=p2;}
    void DiffCoeff(Float_t&diff,Float_t&diff1) {// Get diffusion coefficients
	diff = fDiffCoeff;diff1 = fDiffCoeff1;}
    void    SetDriftSpeed(Float_t p1=7.3) {// Drift velocity
	fDriftSpeed=p1;}
    Float_t DriftSpeed() const {// drift speed
	return fDriftSpeed;}
    void    SetTemperature(Float_t p1=23.) {// Temperature
	fTemperature=p1;}
    Float_t Temperature() const {// Get temperature
	return fTemperature;}
    void    SetDataType(const char *data="simulated") {
	// Type of data - real or simulated
	fDataType=data;}
    const char  *DataType() const {// Get data type
	return fDataType.Data();}
    void SetParamOptions(const char *opt1="same",const char *opt2="same"){
	// Parameters: "same" or read from "file" 
	fParam1=opt1; fParam2=opt2;}
    void   ParamOptions(char *opt1,char *opt2) {// options
	strcpy(opt1,fParam1.Data()); strcpy(opt2,fParam2.Data());}
    void  SetNoiseParam(Float_t n=10., Float_t b=20.){
	// Noise and baseline  // 10 for ALICE with beam test measurements 2001
	fNoise=n; fBaseline=b;}
    void  SetNoiseAfterElectronics(Float_t n=2.38){
	// Noise after electronics (ADC units)
	// 2.36 for ALICE from beam test measurements 2001
	fNoiseAfterEl=n;}
    void  GetNoiseParam(Float_t &n, Float_t &b) {// get noise param
	n=fNoise; b=fBaseline;}
    Float_t  GetNoiseAfterElectronics(){
	// Noise after electronics (ADC units)
	return fNoiseAfterEl;}
    void  SetJitterError(Float_t jitter=20) {
	// set Jitter error (20 um for ALICE from beam test measurements 2001)
	fJitterError=jitter;}
    Float_t  JitterError() {// set Jitter error
	return fJitterError;}
    void  SetDo10to8(Bool_t bitcomp=kTRUE) {
	// set the option for 10 to 8 bit compression
	fBitComp = bitcomp;}
    Bool_t Do10to8() const {// get 10 to 8 compression option
	return fBitComp;}
    void    SetZeroSupp (const char *opt="1D") {
	// Zero-suppression option - could be 1D, 2D or non-ZS 
	fOption=opt;}
    const char *ZeroSuppOption() const {// Get zero-suppression option
	return fOption.Data();}
    void  SetMinVal(Int_t mv=4) {
	// Min value used in 2D - could be used as a threshold setting
	fMinVal = mv;}
    Int_t  MinVal() const {// min val
	return fMinVal;}
    void   SetFilenames(const char *f1="",const char *f2="",const char *f3=""){
	// Set filenames - input, output, parameters ....
	fFileName1=f1; fFileName2=f2; fFileName3=f3;}
    void   Filenames(char *input,char *baseline,char *param) {// Filenames
	strcpy(input,fFileName1.Data());  strcpy(baseline,fFileName2.Data());  
	strcpy(param,fFileName3.Data());}
    void  SetOutputOption(Bool_t write=kFALSE) {// set output option
	fWrite = write;}
    Bool_t OutputOption() const {// output option
	return fWrite;}
    // 
    // Compression parameters
    void  SetCompressParam(Int_t cp[8]); 
    void  GiveCompressParam(Int_t *x);
    //
    // Detector type response methods
    void    SetNSigmaIntegration(Float_t p1=3.) {
	// Set number of sigmas over which cluster disintegration is performed
	fNsigmas=p1;}
    Float_t NSigmaIntegration() const {
	// Get number of sigmas over which cluster disintegration is performed
	return fNsigmas;}
    void SetNLookUp(Int_t p1=121) {
	// Set number of sigmas over which cluster disintegration is performed
	fNcomps=p1;
	fGaus = new TArrayF(fNcomps+1);
	for(Int_t i=0; i<=fNcomps; i++) {
	    Float_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
	    (*fGaus)[i] = exp(-((x*x)/2));
	    //     cout << "fGaus[" << i << "]: " << fGaus->At(i) << endl;
	}
    }
    // Get number of intervals in which the gaussian lookup table is divided
    Int_t GausNLookUp() const {return fNcomps;}
    Float_t IntPH(Float_t) const {// Pulse height from scored quantity (eloss)
	return 0.;}
    Float_t IntXZ(AliITSsegmentation *) const {// Charge disintegration 
	return 0.;}
    Float_t GausLookUp(Int_t i) const  {
	if(i<0 || i>=fNcomps) return 0.;return fGaus->At(i);}
    void SetDeadChannels(Int_t nmodules=0, Int_t nchips=0, Int_t nchannels=0);
    Int_t GetDeadModules() { return fDeadModules; }
    Int_t GetDeadChips() { return fDeadChips; }
    Int_t GetDeadChannels() { return fDeadChannels; }
    Float_t Gain(Int_t mod,Int_t chip,Int_t ch){return fGain[mod][chip][ch]; }
    // these functions should be move to AliITSsegmentationSDD
    const Int_t Modules() const{return fModules;}// Total number of SDD modules
    const Int_t Chips() const{return fChips;} // Number of chips/module
    const Int_t Channels() const { return fChannels;}// Number of channels/chip
    //********
    void    PrintGains();
    void    Print();

 private:
    AliITSresponseSDD(const AliITSresponseSDD &source); // copy constructor
    AliITSresponseSDD& operator=(const AliITSresponseSDD &source); // ass. op.

 protected:
    // these statis const should be move to AliITSsegmentationSDD
    static const Int_t fModules = 520;     // Total number of SDD modules
    static const Int_t fChips = 4;        // Number of chips/module
    static const Int_t fChannels = 64;    // Number of channels/chip
    //*******
    Int_t fDeadModules;                   // Total number of dead SDD modules
    Int_t fDeadChips;                     // Number of dead chips
    Int_t fDeadChannels;                  // Number of dead channels
    Float_t   fGain[fModules][fChips][fChannels];   // Array for channel gains
    Int_t     fCPar[8];        // Hardware compression parameters
    Float_t   fNoise;          // Noise
    Float_t   fBaseline;       // Baseline
    Float_t   fNoiseAfterEl;   // Noise after electronics
    Float_t   fJitterError;    // jitter error
    Float_t   fDynamicRange;   // Set Dynamic Range 
    Float_t   fChargeLoss;     // Set Linear Coefficient for Charge Loss 
    Float_t   fTemperature;    // Temperature 
    Float_t   fDriftSpeed;     // Drift velocity
    Int_t     fElectronics;    // Electronics
    Float_t    fMaxAdc;        // Adc saturation value
    Float_t    fDiffCoeff;     // Diffusion Coefficient (scaling the time)
    Float_t    fDiffCoeff1;    // Diffusion Coefficient (constant term)
    Float_t    fNsigmas;   // Number of sigmas over which charge disintegration
                               // is performed
    TArrayF   *fGaus;          // Gaussian lookup table for signal generation
    Int_t      fNcomps;        // Number of samplings along the gaussian
    Int_t      fMinVal;        // Min value used in 2D zero-suppression algo
    Bool_t     fWrite;         // Write option for the compression algorithms
    Bool_t     fBitComp;       // 10 to 8 bit compression option
    TString    fOption;        // Zero-suppresion option (1D, 2D or none)
    TString    fParam1;        // Read baselines from file option
    TString    fParam2;        // Read compression algo thresholds from file
    TString         fDataType;         // data type - real or simulated
    TString         fFileName1;        // input keys : run, module #
    TString         fFileName2;        // baseline & noise val or output code
                                       // signal or monitored bgr.
    TString         fFileName3;        // param values or output coded signal

    ClassDef(AliITSresponseSDD,3) // SDD response 
    
    };
#endif
