#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H

#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
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
  virtual void    SetDiffCoeff(Float_t p1=0.00433) {
    // Diffusion coefficient
    fDiffCoeff=p1;
  }
  virtual Float_t DiffCoeff() {
    // Get diffusion coefficient
    return fDiffCoeff;
  }
  virtual  void   SetNoiseParam(Float_t n=0., Float_t b=0.) {
    // set noise
    fNoise=n; fBaseline=b;
  }   
  virtual  void   GetNoiseParam(Float_t &n, Float_t &b) {
    // get noise
    n=fNoise; b=fBaseline;
  }   
  virtual void     SetMinVal(Int_t p1=0) {
    // Zero-suppression option threshold 
    fThreshold=p1;
  }
  virtual Int_t MinVal() {
    // Get zero-suppression threshold
    return fThreshold;
  }
  virtual void    SetDataType(char *data="simulated") {
    // Type of data - real or simulated
    fDataType=data;
  }
  virtual char  *DataType() {
    // Get data typer
    return fDataType;
  } 
  
  ClassDef(AliITSresponseSPD,1) // SPD response
    
    protected:
  
  Float_t fDiffCoeff;       // Diffusion Coefficient
  Float_t fNoise;           // Noise value
  Float_t fBaseline;        // Baseline value
  Int_t fThreshold;         // Zero-Suppression threshold
  
  char* fDataType;          //!
                            // Type of data - real or simulated
};

#endif
