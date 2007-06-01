#ifndef ALIITSONLINESDDCMN_H
#define ALIITSONLINESDDCMN_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD noise corrected for common mode analysis   //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSOnlineSDD.h"

class TH2F;
class TGraph;
class AliITSOnlineSDDCMN : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDCMN();
  AliITSOnlineSDDCMN(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDDCMN();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();
  void ReadBaselines();

  void SetMinNoise(Float_t ns=0.001){fMinCorrNoise=ns;}
  void SetMaxNoise(Float_t ns=9.){fMaxCorrNoise=ns;}
  void SetNSigmaNoise(Float_t ns=4.){fNSigmaNoise=ns;}

  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Float_t GetAnodeBaseline(Int_t iAnode) const{ return fBaseline[iAnode];}
  Float_t GetAnodeRawNoise(Int_t iAnode) const{ return fRawNoise[iAnode];}
  Float_t GetAnodeCommonMode(Int_t iAnode) const{ return fCMN[iAnode];}
  Float_t GetAnodeCorrNoise(Int_t iAnode) const{
    if(fNEvents>0) return fSumCorrNoise[iAnode]/fNEvents;
    else return 0;
  }

  Float_t CalcMeanNoise();
  Int_t GetNEvents() const {return fNEvents;}
  void WriteToFXS();

 protected:

 private:
  Int_t fNEvents;
  Bool_t fGoodAnode[fgkNAnodes];
  Float_t fBaseline[fgkNAnodes];
  Float_t fRawNoise[fgkNAnodes];
  Float_t fSumCorrNoise[fgkNAnodes];
  Float_t fCMN[fgkNAnodes];
  Float_t fMinCorrNoise;
  Float_t fMaxCorrNoise;
  Float_t fNSigmaNoise;

  ClassDef(AliITSOnlineSDDCMN,1);
};
#endif
