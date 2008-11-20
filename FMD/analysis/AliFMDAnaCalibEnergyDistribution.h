#ifndef ALIFMDANACALIBENERGYDISTRIBUTION_H
#define ALIFMDANACALIBENERGYDISTRIBUTION_H

#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>

class AliFMDAnaCalibEnergyDistribution : public TObject
{
  
 public:
  
  AliFMDAnaCalibEnergyDistribution();
  void  SetEnergyDistribution(Int_t det, Char_t ring, TH1F* edist);
  TH1F* GetEnergyDistribution(Int_t det, Char_t ring);
  
  
 protected:
  void      Init();
  TObjArray fArray;
  Bool_t    fIsInit;
  ClassDef(AliFMDAnaCalibEnergyDistribution,1);
};

#endif
