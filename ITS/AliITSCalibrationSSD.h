#ifndef ALIITSCALIBRATIONSSD_H
#define ALIITSCALIBRATIONSSD_H
 
#include "AliITSCalibration.h"
#include "AliITSresponseSSD.h"

//////////////////////////////////////////////
// Response class for SSD                   //
//                                          //
//////////////////////////////////////////////
class AliITSCalibrationSSD : public AliITSCalibration {

 public:
    AliITSCalibrationSSD();
    AliITSCalibrationSSD(const char *dataType);
    virtual ~AliITSCalibrationSSD();
    virtual  void   SetNoiseParam(Double_t np, Double_t nn) {
	// set noise par
	fNoiseP=np; fNoiseN=nn;
    }
    virtual void    GetNoiseParam(Double_t &np, Double_t &nn) const {
	// get noise par
	np=fNoiseP; nn=fNoiseN;
    }

    virtual  void   SetNDetParam(Int_t npar) {
	// set number of param
	fNPar=npar;
    }
    virtual  void   SetDetParam(Double_t *par);

    // Parameters options
    virtual Int_t   NDetParam() const {
	// number of param
	return fNPar;
    }
    virtual void    GetDetParam(Double_t *dpar) const; 

    virtual void    SetSigmaSpread(Double_t p1, Double_t p2) {
	// Set sigmas of the charge spread function: Pside-Nside
	// square of (microns)
	fSigmaP=p1; fSigmaN=p2;
    }
    virtual void    SigmaSpread(Double_t &sP, Double_t &sN) const {
	// Get sigmas for the charge spread 
	sP=fSigmaP; sN=fSigmaN;
    }

    virtual void   SetThresholds(Double_t /* a */, Double_t /* b */)
      {NotImplemented("SetThresholds");}
    virtual void   Thresholds(Double_t & /* a */, Double_t & /* b */) const 
      {NotImplemented("Thresholds");}
  
    virtual void SetParamOptions(const char *opt1, const char *opt2) {((AliITSresponseSSD*)fResponse)->SetParamOptions(opt1,opt2);}
    virtual void GetParamOptions(char *opt1,char *opt2) const {((AliITSresponseSSD*)fResponse)->ParamOptions(opt1,opt2);}
    virtual void SetADCpereV(Double_t a=50./30000.0) {((AliITSresponseSSD*)fResponse)->SetADCpereV(a);}
    virtual Double_t GetDEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->DEvToADC(eV);}
    virtual Int_t IEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->IEvToADC(eV);} 

 protected:
    static const Double_t fgkNoiseNDefault; // default for fNoiseN
    static const Double_t fgkNoisePDefault; // default for fNoiseP
    static const Int_t fgkNParDefault; // default for fNPar
    static const Double_t fgkSigmaPDefault; //default for fSigmaP
    static const Double_t fgkSigmaNDefault; //default for fSigmaP
    Int_t   fNPar;            // Number of detector param 
    Double_t *fDetPar;         //[fNPar] Array of parameters

    Double_t fNoiseP;          // Noise on Pside
    Double_t fNoiseN;          // Noise on Nside

    Double_t fSigmaP;          // Sigma charge spread on Pside
    Double_t fSigmaN;          // Sigma charge spread on Nside
    
 private:
    AliITSCalibrationSSD(const AliITSCalibrationSSD &source); // copy constructor
    AliITSCalibrationSSD& operator=(const AliITSCalibrationSSD &source); // ass. op.

    ClassDef(AliITSCalibrationSSD,1) //Response class for SSD
};
#endif
