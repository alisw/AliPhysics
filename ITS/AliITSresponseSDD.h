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
    virtual ~AliITSresponseSDD(){};

    static Float_t DefaultDriftSpeed() {return fgkDriftSpeedDefault;}

    virtual void SetTimeOffset(Float_t to){fTimeOffset = to;}
    virtual Float_t TimeOffset()const {return fTimeOffset;}
    static Float_t DefaultTimeOffset() {return fgkTimeOffsetDefault;}

    virtual void SetADC2keV(Float_t conv){fADC2keV=conv;}
    virtual Float_t ADC2keV()const {return fADC2keV;}
    static Float_t DefaulttADC2keV() {return fgkADC2keVDefault;}

 
    void    SetZeroSupp (const char *opt) {
	// Zero-suppression option - could be ZS or NOTSUPP
	fOption=opt;}
    const char *ZeroSuppOption() const {// Get zero-suppression option
	return fOption.Data();}
    // Detector type response methods

 protected:

    static const TString fgkOptionDefault; // default for fOption
    static const Float_t fgkDriftSpeedDefault; // default for drift speed
    static const Float_t fgkTimeOffsetDefault; // default for fTimeOffset
    static const Float_t fgkADC2keVDefault; // default for fADC2keV

    Float_t  fTimeOffset;     // Time offset due to electronic delays 
    Float_t  fADC2keV;        // Conversion factor from ADC to keV
    TString    fOption;        // Zero-suppresion option (ZS or non-ZS)

 private:

   AliITSresponseSDD(const AliITSresponseSDD &ob); // copy constructor
   AliITSresponseSDD& operator=(const AliITSresponseSDD & /* source */); // ass. op.

    ClassDef(AliITSresponseSDD,13) // Base response class 
    
    };
#endif
