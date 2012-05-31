#ifndef ALIPMDNOISECUT_H
#define ALIPMDNOISECUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;

class AliPMDNoiseCut: public TNamed
{
 public:
  AliPMDNoiseCut();
  AliPMDNoiseCut(const char* name);
  AliPMDNoiseCut(const AliPMDNoiseCut &noisecut);
  AliPMDNoiseCut& operator= (const AliPMDNoiseCut &noisecut);
  virtual ~AliPMDNoiseCut();

  Float_t GetNoiseCut(Int_t imod) const {return fNoiseCut[imod];}

  void  SetNoiseCut(Int_t imod, Float_t cut) {fNoiseCut[imod] = cut;}

  virtual void Print(Option_t *) const;
  
 protected:

  Float_t fNoiseCut[48];


  ClassDef(AliPMDNoiseCut,2) // calibration class for gainfactors
};
#endif
