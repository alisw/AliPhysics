#ifndef ALIITSRESPONSESSD_H
#define ALIITSRESPONSESSD_H
 
#include "AliITSresponse.h"

// response for SSD

//-----------------------------
class AliITSresponseSSD :
public AliITSresponse {

public:
  
  AliITSresponseSSD();
  virtual ~AliITSresponseSSD();
  AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
  AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.
  
  virtual void    SetDiffCoeff(Float_t p1=0.) {
    // Diffusion coefficient
    fDiffCoeff=p1;
  }
  virtual Float_t DiffCoeff() {
    // Get diffusion coefficient
    return fDiffCoeff;
  }
  
  virtual  void   SetNoiseParam(Float_t np=420., Float_t nn=625.) {
    // set noise par
    fNoiseP=np; fNoiseN=nn;
  }
  virtual void    GetNoiseParam(Float_t &np, Float_t &nn) {
    // get noise par
    np=fNoiseP; nn=fNoiseN;
  }
  
  virtual void    SetParamOptions(Option_t *opt1="", Option_t *opt2="") {
    // parameters: "SetInvalid" to simulate the invalid strips
    fOption1=opt1; fOption2=opt2;
  }
  virtual void    ParamOptions(Option_t *&opt1,Option_t *&opt2) {
    // options
    opt1=fOption1; opt2=fOption2;
  }
  
  // Number of parameters to be set
  // 4 couplings, mean number of invalid strips, RMS of invalid strips
  virtual  void   SetNDetParam(Int_t npar=6) {
    // set number of param
    fNPar=npar;
  }

  virtual  void   SetDetParam(Float_t *par);
  
  // Parameters options
  virtual Int_t   NDetParam() {
    // number of param
    return fNPar;
  }
  virtual void    GetDetParam(Float_t *dpar); 

  virtual void    SetDataType(char *data="simulated") {
  // Type of data - real or simulated
    fDataType=data;
  }
  virtual char  *DataType() {
    // Get data type
    return fDataType;
  } 
  
  virtual void    SetSigmaSpread(Float_t p1=3., Float_t p2=2.) {
    // Set sigmas of the charge spread function: Pside-Nside
    // square of (microns)
    fSigmaP=p1; fSigmaN=p2;
  }
  virtual void    SigmaSpread(Float_t &sP, Float_t &sN) {
    // Get sigmas for the charge spread 
    sP=fSigmaP; sN=fSigmaN;
  }
  
protected:
  Int_t   fNPar;            // Number of detector param 
  Float_t *fDetPar;         //!
                            // Array of parameters 
  
  Float_t fNoiseP;          // Noise on Pside
  Float_t fNoiseN;          // Noise on Nside
  
  Float_t fSigmaP;          // Sigma charge spread on Pside
  Float_t fSigmaN;          // Sigma charge spread on Nside
  Float_t fDiffCoeff;       // Diffusion Coefficient
  
  Option_t *fOption1;       //!
                            // Simulate invalid strips option
  Option_t *fOption2;       //! 
                            // Not used for the moment
  
  char*  fDataType;         //!
                            // Type of data - real or simulated
  
  ClassDef(AliITSresponseSSD,1) //Response class for SSD 

    };



#endif







