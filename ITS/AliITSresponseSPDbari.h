#ifndef ALIITSRESPONSESPDBARI_H
#define ALIITSRESPONSESPDBARI_H

#include "AliITSresponse.h"
#include <TString.h>

//----------------------------------------------
//
// ITS response class for SPD
//
class AliITSresponseSPDbari :
  public AliITSresponse {
public:
  
  AliITSresponseSPDbari();
  virtual ~AliITSresponseSPDbari() { 
    // destructror
  }
  //
  // Configuration methods
  //
  
  
  virtual  void   SetDiffCoeff(Float_t p1=0) {
    // 
    fDiffCoeff=p1;
  }
  virtual  Float_t   DiffCoeff() {
    // 
    return fDiffCoeff;
  }
  virtual  void   SetThresholds(Float_t thresh=7.2e-6, Float_t sigma=1.0e-6) {
    // Set Threshold and noise + threshold fluctuations parameter values
    fThresh=thresh; fSigma=sigma;
  }
  virtual  void   Thresholds(Float_t &thresh, Float_t &sigma) {
    // Get Threshold and noise + threshold fluctuations parameter values
    thresh=fThresh; sigma=fSigma;
  }
  virtual  void   SetNoiseParam(Float_t col=0., Float_t row=0.) {
    // set coupling parameters
    fCouplCol=col; fCouplRow=row;
  }   
  virtual  void   GetNoiseParam(Float_t &col, Float_t &row) {
    // get coupling parameters
    col=fCouplCol; row=fCouplRow;
  }       
  virtual void    SetDataType(char *data="simulated") {
    // Type of data - real or simulated
    fDataType=data;
  }
  virtual const char  *DataType() {
    // Get data typer
    return fDataType.Data();
  } 
  
  ClassDef(AliITSresponseSPDbari,1) // SPD response
    
    protected:
  
  Float_t fDiffCoeff;       // Sigma diffusion coefficient (not used) 
  Float_t fThresh;          // Threshold value
  Float_t fSigma;           // Noise + threshold fluctuations value
  Float_t fCouplCol;        // Coupling probability along a column
  Float_t fCouplRow;        // Coupling probability along a row

  TString fDataType;        // Type of data - real or simulated
};

#endif
