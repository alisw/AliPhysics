#ifndef ALIITSCALIBRATIONSSD_H
#define ALIITSCALIBRATIONSSD_H
 
#include "AliITSCalibration.h"
#include "AliITSNoiseSSD.h"
#include "AliITSGainSSD.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSresponseSSD.h"
#include "TArrayF.h"
#include "TArrayI.h"

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

    void AddNoiseP(Int_t c, Float_t n) { fNoise->AddNoiseP(c,n);}       
    TArrayF GetNoiseP() {return fNoise->GetNoiseP(); }
    Float_t GetNoiseP(Int_t n) {return fNoise->GetNoiseP(n); }
    void AddNoiseN(Int_t c, Float_t n) { fNoise->AddNoiseN(c,n);}
    TArrayF GetNoiseN() {return fNoise->GetNoiseN(); }
    Float_t GetNoiseN(Int_t n) {return fNoise->GetNoiseN(n); }
    void SetNoise( AliITSNoiseSSD* noise) {fNoise=noise;}

    void AddGainP(Int_t c, Float_t n) { fGain->AddGainP(c,n);}       
    TArrayF GetGainP() {return fGain->GetGainP(); }
    Float_t GetGainP(Int_t n) {return fGain->GetGainP(n); }
    void AddGainN(Int_t c, Float_t n) { fGain->AddGainN(c,n);}
    TArrayF GetGainN() {return fGain->GetGainN(); }
    Float_t GetGainN(Int_t n) {return fGain->GetGainN(n); }
    void SetGain( AliITSGainSSD* gain) {fGain=gain;}

    TArrayI GetBadPChannelsList() { return fBadChannels->GetBadPChannelsList(); } const
    TArrayI GetBadNChannelsList() { return fBadChannels->GetBadNChannelsList(); } const
    void SetBadChannels( AliITSBadChannelsSSD* badchannels) {fBadChannels=badchannels;}

    void SetNoisePThreshold(Int_t threshold) { fNoisePThreshold = threshold;}
    void AddNoisyPChannel(Int_t c, Int_t n) { fNoisyPChannelsList.AddAt(n,c);}
    TArrayI GetNoisyPChannelsList() const {return fNoisyPChannelsList; }
    void SetNoiseNThreshold(Int_t threshold) { fNoiseNThreshold = threshold;}
    void AddNoisyNChannel(Int_t c, Int_t n) { fNoisyNChannelsList.AddAt(n,c);}
    TArrayI GetNoisyNChannelsList() const {return fNoisyNChannelsList; }

    void SetNDeadPChannelsList(Int_t n) { fDeadPChannelsList.Set(n); }
    void AddDeadPChannel(Int_t c, Int_t n) { fDeadPChannelsList.AddAt(n,c);}
    TArrayI GetDeadPChannelsList() const {return fDeadPChannelsList; }
    void SetNDeadNChannelsList(Int_t n) { fDeadNChannelsList.Set(n); }
    void AddDeadNChannel(Int_t c, Int_t n) { fDeadNChannelsList.AddAt(n,c);}
    TArrayI GetDeadNChannelsList() const {return fDeadNChannelsList; }
    //


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
    virtual void SetADCpereV(Double_t a=120./24888.9) {((AliITSresponseSSD*)fResponse)->SetADCpereV(a);}
    virtual Double_t GetDEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->DEvToADC(eV);}
    virtual Int_t IEvToADC(Double_t eV) const {return ((AliITSresponseSSD*)fResponse)->IEvToADC(eV);} 

    virtual void SetKeVperADC(Double_t a=86.4/120.0) {((AliITSresponseSSD*)fResponse)->SetKeVperADC(a);}
    virtual Double_t ADCToKeV(Double_t adc) const {return ((AliITSresponseSSD*)fResponse)->ADCToKeV(adc);}

  virtual Double_t GetCouplingPR() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingPR();}
  virtual Double_t GetCouplingPL() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingPL();}
  virtual Double_t GetCouplingNR() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingNR();}
  virtual Double_t GetCouplingNL() const {return ((AliITSresponseSSD*)fResponse)->GetCouplingNL();}
  virtual void SetCouplings(Double_t pr, Double_t pl, Double_t nr, Double_t nl) 
    { ((AliITSresponseSSD*)fResponse)->SetCouplings(pr,pl,nr,nl);}

  virtual Int_t GetZSThreshold() const {return ((AliITSresponseSSD*)fResponse)->GetZSThreshold();}
  virtual void SetZSThreshold(Int_t zsth) 
    { ((AliITSresponseSSD*)fResponse)->SetZSThreshold(zsth);}

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
    
    AliITSNoiseSSD *fNoise;
    AliITSGainSSD *fGain;
    AliITSBadChannelsSSD *fBadChannels;

    TArrayF fGainP;           // Gain for P side channels
    TArrayF fGainN;           // Gain for N side channels

    TArrayF fNoisP;           // Noise for P side channels
    TArrayF fNoisN;           // Noise for N side channels

    Float_t fNoisePThreshold;     // need to decide if channel is noisy  
    TArrayI  fNoisyPChannelsList; // list of P side noisy channels
    Float_t fNoiseNThreshold;     // need to decide if channel is noisy  
    TArrayI  fNoisyNChannelsList; // list of N side noisy channels

    TArrayI  fDeadNChannelsList;  // list of P side dead channels
    TArrayI  fDeadPChannelsList;  // list of N side dead channels

 private:
    AliITSCalibrationSSD(const AliITSCalibrationSSD &source); // copy constructor
    AliITSCalibrationSSD& operator=(const AliITSCalibrationSSD &source); // ass. op.

    ClassDef(AliITSCalibrationSSD,1) //Response class for SSD
};
#endif
