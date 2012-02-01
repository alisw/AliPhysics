#ifndef ALIFMDANACALIBENERGYDISTRIBUTION_H
#define ALIFMDANACALIBENERGYDISTRIBUTION_H
#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>
class TBrowser;
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
  Int_t GetNetaBins() const { return fNetaBins;}
  void  SetEtaLimits(Float_t eta_min, Float_t eta_max) {fEtaMin = eta_min; fEtaMax = eta_max;}
  void  SetEnergyDistributionUser(Int_t det, Char_t ring, Float_t eta, TH1F* edist);
  void  SetEnergyDistribution(Int_t det, Char_t ring, Int_t etabin, TH1F* edist);
  
  TH1F* GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta);
  void  SetEmptyEnergyDistribution(Int_t det, Char_t ring, TH1F* edist);
  TH1F* GetEmptyEnergyDistribution(Int_t det, Char_t ring);
  void  SetRingEnergyDistribution(Int_t det, Char_t ring, TH1F* edist);
  TH1F* GetRingEnergyDistribution(Int_t det, Char_t ring);
  Bool_t IsFolder() const { return kTRUE; }
  void Browse(TBrowser* b);
 protected:
  void      Init();
  TObjArray fArray;       // Overall array
  TObjArray fEmptyArray;  // Array of empty events energy dep
  TObjArray fRingArray;   // Array of rings
  Bool_t    fIsInit;      // Are we init?
  Int_t     fNetaBins;    // Number of eta bins (intervals)
  Float_t   fEtaMax;      // Max Eta
  Float_t   fEtaMin;      // Min Eta
  
  ClassDef(AliFMDAnaCalibEnergyDistribution,3);
};

#endif
// Local Variables:
//   mode: C++
// End:
