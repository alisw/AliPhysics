#ifndef ALIITSRESPONSESPDDUBNA_H
#define ALIITSRESPONSESPDDUBNA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TRandom.h>
class TString;


#include "AliITSresponse.h"

//----------------------------------------------
// ITS response class for SPD
class AliITSresponseSPDdubna : public AliITSresponse {
 public:
  AliITSresponseSPDdubna();
  virtual ~AliITSresponseSPDdubna(){}// destructror
 
//abstract methods in AliITSresponse not implemented in this class
    virtual void    GiveCompressParam(Int_t *) const
      {NotImplemented("GiveCompressParam");}
    virtual void    SetDriftSpeed(Float_t /* p1 */)
      {NotImplemented("SetDriftSpeed");}
    virtual Float_t DriftSpeed() const 
      {NotImplemented("DrifSpeed"); return 0.;}
    virtual void   SetElectronics(Int_t /* i */) 
                    {NotImplemented("SetElectronics");}
    virtual Int_t Electronics() const {NotImplemented("Electronics"); return 0;}
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
    virtual void    SetDiffCoeff(Float_t /* a */, Float_t /* b */)
      {NotImplemented("SetDiffCoeff");}
    virtual void    DiffCoeff(Float_t & /* a */,Float_t & /* b */) const
      {NotImplemented("DiffCoeff");}
    virtual  void   SetDetParam(Float_t *)
      {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Float_t *) const 
      {NotImplemented("GetDetParam");}
    virtual  void   SetNDetParam(Int_t /* n */)
      {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam() const
      {NotImplemented("NDetParam"); return 0;}
    virtual void   SetParamOptions(const char* /* a */,const char* /* b */)
      {NotImplemented("SetParamOptions");}
    virtual void   ParamOptions(char *,char*) const
      {NotImplemented("ParamOptions");} 
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
    virtual void    SetSigmaSpread(Float_t, Float_t) 
      {NotImplemented("SetSigmaSpread");}
    virtual void    SigmaSpread(Float_t & /* p1 */,Float_t & /* p2 */) const 
      {NotImplemented("SigmaSpread");}

// Implementation of virtual member functions declared in AliITSresponse 
  // set noise and baseline
  virtual void SetNoiseParam(Float_t n, Float_t b) {
    fNoise=n; fBaseline=b;}
  virtual void GetNoiseParam(Float_t &n, Float_t &b) const {
    n=fNoise; b=fBaseline;} // get noise and baseline
  // Zero-suppression option threshold
  virtual void SetThresholds(Float_t th, Float_t /* dum */) {fThreshold=th;}
  virtual void Thresholds(Float_t & thr, Float_t & /*dum */) const
    {thr = fThreshold;}

  // Prints out the content of this class in ASCII format.
  virtual void Print(ostream *os) const;

  // Reads in the content of this class in the format of Print
  virtual void Read(istream *is);

//Declaration of member functions peculiar to this class
  // Applies a random noise and addes the baseline
  Float_t ApplyBaselineAndNoise() const {return fBaseline+
                                                 fNoise*gRandom->Gaus();}
  //Returns just baseline value
  Float_t GetBaseline() const {return fBaseline;}
  // Returns just noise value
  Float_t GetNoise() const {return fNoise;} 
  // Get zero-suppression threshold
  Float_t GetThreshold() const {return fThreshold;}
  // Sets the coupling probabilities for columns and rows
  void SetCouplings(Float_t ccol=0.,Float_t crow=0.){fCouplCol=ccol;
                                                       fCouplRow=crow;}
  // Gets the coupling probabilities for columns and rows
  void GetCouplings(Float_t &ccol,Float_t &crow) const 
              {ccol=fCouplCol; crow=fCouplRow;}
  // Returns the fraction of dead pixels
  Float_t GetFractionDeadPixels() const {return fDeadPixels;}
  // Sets the fraction of dead pixels
  void SetFractionDeadPixels(Float_t f=0.01){fDeadPixels = f;}
  // Returns a logical if a pixels is dead or not
  Bool_t IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const ;

  ClassDef(AliITSresponseSPDdubna,3) // SPD response
    
 protected:
  static const Float_t fgkNoiseDefault; // default for fNoise
  static const Float_t fgkThresholdDefault; //default for fThreshold
  Float_t fNoise;           // Noise value
  Float_t fBaseline;        // Baseline value
  Float_t fCouplCol;        // Coupling Probability along columns
  Float_t fCouplRow;        // Coupling Probability along rows
  Float_t fThreshold;       // Zero-Suppression threshold
  Float_t fDeadPixels;      // Fraction of dead pixels.
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSresponseSPDdubna &source);
istream& operator>>(istream &os,AliITSresponseSPDdubna &source);

#endif
