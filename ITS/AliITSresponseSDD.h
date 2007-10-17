#ifndef ALIITSRESPONSESDD_H
#define ALIITSRESPONSESDD_H
 
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 

#include <AliITSresponse.h>
#include <TArrayF.h>

/* $Id$ */

/////////////////////////////////////////////////////////////
//  Base settings for the ITS response classes.            //  
//  The data member of this class are static and set once  //
//  for all the modules.                                   //    
///////////////////////////////////////////////////////////// 

class AliITSresponseSDD : public AliITSresponse {
  public:

    AliITSresponseSDD();
    virtual ~AliITSresponseSDD();

    virtual void SetElectronics(Int_t p1=1) 
      {fElectronics=p1;  /* Electronics: Pascal (1) or OLA (2) */ }
    virtual Int_t Electronics()  const {// Electronics: 1 = Pascal; 2 = OLA
	return fElectronics;}
    virtual void    SetMaxAdc(Double_t p1) {// Adc-count saturation value
	fMaxAdc=p1;}
    virtual Float_t MaxAdc() const  {// Get maximum Adc-count value
	return fMaxAdc;}
    virtual void    SetChargeLoss(Double_t p1) {
	// Set Linear Charge Loss Steepness  // 0.01 for 20%
	fChargeLoss=p1;}
    virtual Float_t ChargeLoss() const {// Get Charge Loss Coefficient
	return fChargeLoss;}
    virtual void    SetDynamicRange(Double_t p1) {// Set Dynamic Range
	fDynamicRange = p1;}
    virtual Float_t DynamicRange() const {// Get Dynamic Range
	return fDynamicRange;}

    static Float_t DefaultDriftSpeed() {return fgkDriftSpeedDefault;}

    virtual void SetTimeOffset(Float_t to){fTimeOffset = to;}
    virtual Float_t TimeOffset()const {return fTimeOffset;}
    static Float_t DefaultTimeOffset() {return fgkTimeOffsetDefault;}

    virtual void SetADC2keV(Float_t conv){fADC2keV=conv;}
    virtual Float_t ADC2keV()const {return fADC2keV;}
    static Float_t DefaulttADC2keV() {return fgkADC2keVDefault;}

    virtual void SetParamOptions(const char *opt1,const char *opt2) {
	// Parameters: "same" or read from "file" 
	fParam1=opt1; fParam2=opt2;}
    virtual void   ParamOptions(char *opt1,char *opt2) const {// options
	strcpy(opt1,fParam1.Data()); strcpy(opt2,fParam2.Data());}
 
    virtual Bool_t Do10to8() const {// get 10 to 8 compression option
	return fBitComp;}
    void    SetZeroSupp (const char *opt) {
	// Zero-suppression option - could be 1D, 2D or non-ZS 
	fOption=opt;}
    const char *ZeroSuppOption() const {// Get zero-suppression option
	return fOption.Data();}
    // Detector type response methods
    virtual void    SetNSigmaIntegration(Double_t p1) {
	// Set number of sigmas over which cluster disintegration is performed
	fNsigmas=p1;}
    virtual Float_t NSigmaIntegration() const {
	// Get number of sigmas over which cluster disintegration is performed
	return fNsigmas;}
    virtual void SetNLookUp(Int_t p1);
    // Get number of intervals in which the gaussian lookup table is divided
    virtual Int_t GausNLookUp() const {return fNcomps;}
    virtual Float_t GausLookUp(Int_t i)  {
      if (!fGaus) SetNLookUp(fgkNcompsDefault);
      if(i<0 || i>=fNcomps) return 0.;return fGaus->At(i);
    }
   
    Int_t Convert8to10(Int_t signal) const; //undo 10 to 8 bit SDD compresion

    void  SetJitterError(Float_t jitter=20) {
	// set Jitter error (20 um for ALICE from beam test measurements 2001)
	fJitterError=jitter;}
    Float_t  JitterError() const {// set Jitter error
	return fJitterError;}
    void  SetDo10to8(Bool_t bitcomp=kTRUE) {
	// set the option for 10 to 8 bit compression
	fBitComp = bitcomp;}

 protected:

    static const Int_t fgkMaxAdcDefault; // default for fMaxAdc
    static const Float_t fgkDynamicRangeDefault; // default for fDynamicRange
    static const Float_t fgkfChargeLossDefault; // default for fChargeLoss
    static const Float_t fgkDiffCoeffDefault; // default for fDiffCoeff
    static const Float_t fgkDiffCoeff1Default; // default for fDiffCoeff1 
    static const TString fgkParam1Default; // default for fParam1
    static const TString fgkParam2Default; // default for fParam2
    static const TString fgkOptionDefault; // default for fOption
    static const Float_t fgkDriftSpeedDefault; // default for drift speed
    static const Float_t fgkTimeOffsetDefault; // default for fTimeOffset
    static const Float_t fgkADC2keVDefault; // default for fADC2keV
    static const Float_t fgkNsigmasDefault; //default for fNsigmas
    static const Int_t fgkNcompsDefault; //default for fNcomps
    

    Float_t  fJitterError;    // jitter error
    Float_t  fDynamicRange;   // Set Dynamic Range 
    Float_t  fChargeLoss;     // Set Linear Coefficient for Charge Loss 
    Float_t  fTimeOffset;     // Time offset due to electronic delays 
    Float_t  fADC2keV;        // Conversion factor from ADC to keV
    Int_t    fElectronics;    // Electronics
    Float_t  fMaxAdc;         // Adc saturation value
    Float_t  fNsigmas;   // Number of sigmas over which charge disintegration
                          // is performed
    TArrayF   *fGaus;          // Gaussian lookup table for signal generation
    Int_t      fNcomps;        // Number of samplings along the gaussian
    Bool_t     fBitComp;       // 10 to 8 bit compression option
    TString    fOption;        // Zero-suppresion option (1D, 2D or none)
    TString    fParam1;        // Read baselines from file option
    TString    fParam2;        // Read compression algo thresholds from file

 private:

   AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
   AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.

    ClassDef(AliITSresponseSDD,10) // Base response class 
    
    };
#endif
