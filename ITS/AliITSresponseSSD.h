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
    virtual void    SetDiffCoeff(Float_t p1,Float_t /* dummy */) {
	// Diffusion coefficient
      fDiffCoeff=p1; }
    virtual void    DiffCoeff(Float_t &diffc,Float_t & /* dummy */) const {
	// Get diffusion coefficient
	diffc= fDiffCoeff;
    }

    virtual  void   SetNoiseParam(Float_t np, Float_t nn) {
	// set noise par
	fNoiseP=np; fNoiseN=nn;
    }
    virtual void    GetNoiseParam(Float_t &np, Float_t &nn) const {
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

    virtual  void   SetDetParam(Float_t *par);

    // Parameters options
    virtual Int_t   NDetParam() const {
	// number of param
	return fNPar;
    }
    virtual void    GetDetParam(Float_t *dpar) const; 

    virtual void    SetSigmaSpread(Float_t p1, Float_t p2) {
	// Set sigmas of the charge spread function: Pside-Nside
	// square of (microns)
	fSigmaP=p1; fSigmaN=p2;
    }
    virtual void    SigmaSpread(Float_t &sP, Float_t &sN) const {
	// Get sigmas for the charge spread 
	sP=fSigmaP; sN=fSigmaN;
    }

//Declaration of member functions peculiar to this class

// Sets electron-hole pairs to ADC value conversion factor
// are rather arbitrary (need tuning)
// minimum ionizing particle--> ~30000 pairs--> ADC channel 50
    void SetADCpereV(Float_t a=50./30000.0){fADCpereV = a;}
// Converts electron-hole pairs to
// ADC value conversion factor are rather arbitrary (need tuning)
// minimum ionizing particle--> ~30000 pairs--> ADC channel 50
    Double_t DEvToADC(Double_t eV){return eV*fADCpereV;}
    Int_t IEvToADC(Double_t eV){ // Converts electron-hole pairs to
	// ADC value
	return ((Int_t) DEvToADC(eV)); }

 //abstract methods in AliITSresponse not implemented in this class
    virtual void    SetDriftSpeed(Float_t /* p1 */)
      {NotImplemented("SetDriftSpeed");}
    virtual Float_t DriftSpeed() const 
      {NotImplemented("DrifSpeed"); return 0.;}
    virtual void   SetElectronics(Int_t /* i */) 
                    {NotImplemented("SetElectronics");}
    virtual Int_t Electronics() const {NotImplemented("Electronics"); return 0;}
    virtual void    GiveCompressParam(Int_t *) const
      {NotImplemented("GiveCompressParam");}
    virtual void SetMaxAdc(Float_t /* adc */) {NotImplemented("SetMaxAdc");}
    virtual Float_t MaxAdc() const {NotImplemented("MaxAdc"); return 0.;}
    virtual void    SetDynamicRange(Float_t /*dr */) 
      {NotImplemented("SetDynamicRange");}
    virtual Float_t DynamicRange() const 
      {NotImplemented("DynamicRange"); return 0.;}
    virtual void    SetChargeLoss(Float_t /* cl */)
      {NotImplemented("SetChargeLoss"); }
    virtual Float_t ChargeLoss() const 
      {NotImplemented("ChargeLoss"); return 0.;}
    virtual void   SetThresholds(Float_t /* a */, Float_t /* b */)
      {NotImplemented("SetThresholds");}
    virtual void   Thresholds(Float_t & /* a */, Float_t & /* b */) const 
      {NotImplemented("Thresholds");}
    virtual void   SetZeroSupp(const char*)
      {NotImplemented("SetZeroSupp");}
    virtual const char *ZeroSuppOption() const 
      {NotImplemented("ZeroSuppression"); return "";}
    virtual void    SetNSigmaIntegration(Float_t)
      {NotImplemented("SetNSigmaIntegration");}
    virtual Float_t NSigmaIntegration() const
      {NotImplemented("NSigmaIntegration"); return 0.;}
    virtual void    SetNLookUp(Int_t) 
      {NotImplemented("SetNLookUp");}
  
 protected:
    static const Float_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const TString fgkOption1Default; // default for fOption1
    static const TString fgkOption2Default; // default for fOption2
    static const Float_t fgkNoiseNDefault; // default for fNoiseN
    static const Float_t fgkNoisePDefault; // default for fNoiseP
    static const Int_t fgkNParDefault; // default for fNPar
    static const Float_t fgkSigmaPDefault; //default for fSigmaP
    static const Float_t fgkSigmaNDefault; //default for fSigmaP
    Int_t   fNPar;            // Number of detector param 
    Float_t *fDetPar;         //[fNPar] Array of parameters

    Float_t fNoiseP;          // Noise on Pside
    Float_t fNoiseN;          // Noise on Nside

    Float_t fSigmaP;          // Sigma charge spread on Pside
    Float_t fSigmaN;          // Sigma charge spread on Nside
    Float_t fDiffCoeff;       // Diffusion Coefficient

    Float_t fADCpereV;        // Constant to convert eV to ADC.

    TString fOption1;         // Simulate invalid strips option
    TString fOption2;         // Not used for the moment

 private:
    AliITSresponseSSD(const AliITSresponseSSD &source); // copy constructor
    AliITSresponseSSD& operator=(const AliITSresponseSSD &source); // ass. op.

    ClassDef(AliITSresponseSSD,1) //Response class for SSD
};
#endif
