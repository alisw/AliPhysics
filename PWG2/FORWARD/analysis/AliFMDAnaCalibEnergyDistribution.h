#ifndef ALIFMDANACALIBENERGYDISTRIBUTION_H
#define ALIFMDANACALIBENERGYDISTRIBUTION_H

#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>

/**
 * @ingroup FMD_ana
 * @brief Find most-probable value of MIP peak for one or more
 * particles. 
 * 
 */
class AliFMDAnaCalibEnergyDistribution : public TObject
{
  
 public:
  
  AliFMDAnaCalibEnergyDistribution();
  void  SetNetaBins(Int_t nbins) {fNetaBins = nbins;}
  Int_t GetNetaBins() { return fNetaBins;}
  void  SetEtaLimits(Float_t eta_min, Float_t eta_max) {fEtaMin = eta_min; fEtaMax = eta_max;}
  void  SetEnergyDistribution(Int_t det, Char_t ring, Float_t eta, TH1F* edist);
  void  SetEnergyDistribution(Int_t det, Char_t ring, Int_t etabin, TH1F* edist);
  
  TH1F* GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta);
  void  SetEmptyEnergyDistribution(Int_t det, Char_t ring, TH1F* edist);
  TH1F* GetEmptyEnergyDistribution(Int_t det, Char_t ring);
  void  SetRingEnergyDistribution(Int_t det, Char_t ring, TH1F* edist);
  TH1F* GetRingEnergyDistribution(Int_t det, Char_t ring);
 protected:
  void      Init();
  TObjArray fArray;
  TObjArray fEmptyArray;
  TObjArray fRingArray;
  Bool_t    fIsInit;
  Int_t     fNetaBins;
  Float_t   fEtaMax;
  Float_t   fEtaMin;
  
  ClassDef(AliFMDAnaCalibEnergyDistribution,3);
};

#endif
// Local Variables:
//   mode: C++
// End:
