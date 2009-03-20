#ifndef ALIFMDANACALIBENERGYDISTRIBUTION_H
#define ALIFMDANACALIBENERGYDISTRIBUTION_H

#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>

class AliFMDAnaCalibEnergyDistribution : public TObject
{
  
 public:
  
  AliFMDAnaCalibEnergyDistribution();
  void  SetNetaBins(Int_t nbins) {fNetaBins = nbins;}
  Int_t GetNetaBins() { return fNetaBins;}
  void  SetEtaLimits(Float_t eta_min, Float_t eta_max) {fEtaMin = eta_min; fEtaMax = eta_max;}
  void  SetEnergyDistribution(Int_t det, Char_t ring, Float_t eta, TH1F* edist);
  TH1F* GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta);

  
 protected:
  void      Init();
  TObjArray fArray;
  Bool_t    fIsInit;
  Int_t     fNetaBins;
  Float_t   fEtaMax;
  Float_t   fEtaMin;
  
  ClassDef(AliFMDAnaCalibEnergyDistribution,1);
};

#endif
