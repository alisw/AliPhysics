#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H

#include "AliITSresponse.h"
#include <TString.h>

//----------------------------------------------
//
// ITS response class for SPD
//
class AliITSresponseSPD :
  public AliITSresponse {
public:
  
  AliITSresponseSPD();
  virtual ~AliITSresponseSPD() { 
    // destructror
  }
  //
  // Configuration methods
  //
  virtual void    SetDiffCoeff(Float_t p1=0.00433,Float_t dummy=0.) {
    // Diffusion coefficient
    fDiffCoeff=p1;
  }
  virtual void DiffCoeff(Float_t &diffc,Float_t &dummy) {
    // Get diffusion coefficient
    diffc= fDiffCoeff;
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
  
  ClassDef(AliITSresponseSPD,1) // SPD response
    
    protected:
  
  Float_t fDiffCoeff;       // Diffusion Coefficient
  Float_t fNoise;           // Noise value
  Float_t fBaseline;        // Baseline value
  Int_t fThreshold;         // Zero-Suppression threshold
  
  TString fDataType;        // Type of data - real or simulated
};

#endif
