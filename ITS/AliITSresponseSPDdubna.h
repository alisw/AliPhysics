#ifndef ALIITSRESPONSESPDDUBNA_H
#define ALIITSRESPONSESPDDUBNA_H

#include "AliITSresponse.h"
#include <TString.h>

//----------------------------------------------
//
// ITS response class for SPD
//
class AliITSresponseSPDdubna :  public AliITSresponse {
 public:
    AliITSresponseSPDdubna();
    virtual ~AliITSresponseSPDdubna() { 
	// destructror
    }
    //
    // Configuration methods
    //
    virtual void SetDiffCoeff(Float_t p1=7.877e-3/*0.00433*/,Float_t dummy=0.){
	// Diffusion coefficient
	fDiffCoeff = p1;dummy = 0.;}
    virtual Double_t DiffusionSigma(Double_t dy);
    virtual void DiffCoeff(Float_t &diffc,Float_t &dummy) {
	// Get diffusion coefficient
	diffc= fDiffCoeff;dummy = 0.0;
    }
    virtual  void   SetNoiseParam(Float_t n=200., Float_t b=0.) {
	// set noise and baseline
	fNoise=n; fBaseline=b;
    }   
    virtual  void   GetNoiseParam(Float_t &n, Float_t &b) {
	// get noise and baseline
	n=fNoise; b=fBaseline;
    }   
    virtual void     SetMinVal(Int_t p1=2000) {
	// Zero-suppression option threshold 
	fThreshold=p1;
    }
    virtual Int_t MinVal() {
	// Get zero-suppression threshold
	return fThreshold;
    }
    virtual void    SetDataType(const char *data="simulated") {
	// Type of data - real or simulated
	fDataType=data;
    }
    virtual const char  *DataType() const {
	// Get data typer
	return fDataType.Data();
    } 
    virtual void SetGeVToCharge(Float_t e = 2.778E+08) {
	// sets the conversion factor to go from Energy GeV to charge
	// (electrons).
	fGeVtoElec = e;
    }
    virtual Float_t GetGeVToCharge() {
	// Returns the conversion factor to go from Energy GeV to charge
	// (electrons).
	return fGeVtoElec;
    }
    virtual const Float_t GeVToCharge(Float_t e) const {
	// Converts deposited energy into electrons in Si.
	return e*fGeVtoElec;
    }

 protected:
    Float_t fDiffCoeff;       // Diffusion Coefficient
    Float_t fNoise;           // Noise value
    Float_t fBaseline;        // Baseline value
    Int_t   fThreshold;       // Zero-Suppression threshold
    Float_t fGeVtoElec;       // Conversion factor from GeV to electons
    TString fDataType;        // Type of data - real or simulated

    ClassDef(AliITSresponseSPDdubna,2) // SPD response
};
#endif
