#ifndef ALIITSRESPONSESSD_H
#define ALIITSRESPONSESSD_H
 
#include "AliITSresponse.h"


// response for SSD

//-----------------------------
class AliITSresponseSSD : public AliITSresponse {

 public:
    AliITSresponseSSD();
    AliITSresponseSSD(const char *dataType);
    virtual ~AliITSresponseSSD();

// Implementation of virtual member functions declared in AliITSresponse 
    virtual void    SetDiffCoeff(Double_t p1,Double_t /* dummy */) {
	// Diffusion coefficient
      fDiffCoeff=p1; }
    virtual void    DiffCoeff(Double_t &diffc,Double_t & /* dummy */) const {
	// Get diffusion coefficient
	diffc= fDiffCoeff;
    }

    virtual  void   SetNoiseParam(Double_t np, Double_t nn) {
	// set noise par
	fNoiseP=np; fNoiseN=nn;
    }
    virtual void    GetNoiseParam(Double_t &np, Double_t &nn) const {
	// get noise par
	np=fNoiseP; nn=fNoiseN;
    }

    virtual void    SetParamOptions(const char *opt1, const char *opt2) {
	// parameters: "SetInvalid" to simulate the invalid strips
	fOption1=opt1; fOption2=opt2;
    }
    virtual void    ParamOptions(char *opt1,char *opt2) const {
	// options
	strcpy(opt1,fOption1.Data());  strcpy(opt2,fOption2.Data());
    }

    // Number of parameters to be set
    // 4 couplings, mean number of invalid strips, RMS of invalid strips
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

//Declaration of member functions peculiar to this class

// Sets electron-hole pairs to ADC value conversion factor
// are rather arbitrary (need tuning)
// minimum ionizing particle--> ~30000 pairs--> ADC channel 50
    void SetADCpereV(Double_t a=50./30000.0){fADCpereV = a;}
// Converts electron-hole pairs to
// ADC value conversion factor are rather arbitrary (need tuning)
// minimum ionizing particle--> ~30000 pairs--> ADC channel 50
    Double_t DEvToADC(Double_t eV){return eV*fADCpereV;}
    Int_t IEvToADC(Double_t eV){ // Converts electron-hole pairs to
	// ADC value
	return ((Int_t) DEvToADC(eV)); }

 //abstract methods in AliITSresponse not implemented in this class
    virtual void    SetDriftSpeed(Double_t /* p1 */)
      {NotImplemented("SetDriftSpeed");}
    virtual Double_t DriftSpeed() const 
      {NotImplemented("DrifSpeed"); return 0.;}
    virtual void    GiveCompressParam(Int_t *) const
      {NotImplemented("GiveCompressParam");}
    virtual Double_t ChargeLoss() const 
      {NotImplemented("ChargeLoss"); return 0.;}
    virtual void   SetThresholds(Double_t /* a */, Double_t /* b */)
      {NotImplemented("SetThresholds");}
    virtual void   Thresholds(Double_t & /* a */, Double_t & /* b */) const 
      {NotImplemented("Thresholds");}
    virtual void   SetZeroSupp(const char*)
      {NotImplemented("SetZeroSupp");}
    virtual const char *ZeroSuppOption() const 
      {NotImplemented("ZeroSuppression"); return "";}
    virtual void    SetNSigmaIntegration(Double_t)
      {NotImplemented("SetNSigmaIntegration");}
    virtual Double_t NSigmaIntegration() const
      {NotImplemented("NSigmaIntegration"); return 0.;}
    virtual void    SetNLookUp(Int_t) 
      {NotImplemented("SetNLookUp");}
  
 protected:
    static const Double_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const TString fgkOption1Default; // default for fOption1
    static const TString fgkOption2Default; // default for fOption2
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
    Double_t fDiffCoeff;       // Diffusion Coefficient

    Double_t fADCpereV;        // Constant to convert eV to ADC.

    TString fOption1;         // Simulate invalid strips option
    TString fOption2;         // Not used for the moment

 private:
    AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
    AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.

    ClassDef(AliITSresponseSSD,2) //Response class for SSD
};
#endif
