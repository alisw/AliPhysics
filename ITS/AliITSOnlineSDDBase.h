#ifndef ALIITSONLINESDDBASE_H
#define ALIITSONLINESDDBASE_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD baseline and noise analysis                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include"AliITSOnlineSDD.h"
#include<TMath.h>

class TH2F;
class TGraph;
class AliITSOnlineSDDBase : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDBase();
  AliITSOnlineSDDBase(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDDBase();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();

  void SetMinBaseline(Float_t bas=10.){fMinBaseline=bas;}
  void SetMaxBaseline(Float_t bas=150.){fMaxBaseline=bas;}
  void SetMinRawNoise(Float_t ns=0.001){fMinRawNoise=ns;}
  void SetMaxRawNoise(Float_t ns=9.){fMaxRawNoise=ns;}
  void SetNSigmaNoise(Float_t ns=4.){fNSigmaNoise=ns;}

  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Float_t GetAnodeBaseline(Int_t iAnode) const{
    if(fNEvents>0) return fSumBaseline[iAnode]/fNEvents;
    else return 0;
  }
  Float_t GetAnodeRawNoise(Int_t iAnode) const{
    if(fNEvents>0) return TMath::Sqrt(fSumRawNoise[iAnode]/fNEvents-TMath::Power(GetAnodeBaseline(iAnode),2));
    
    else return 0;
  }

  Float_t CalcMeanRawNoise();
  Float_t GetAnodeCommonMode(Int_t iAnode) const{
    if(fNEvents>0) return fSumCMN[iAnode]/fNEvents;
    else return 0;
  }
  Int_t GetNEvents() const {return fNEvents;}
  void WriteToFXS();

 protected:

 private:
  Int_t fNEvents;
  Bool_t fGoodAnode[fgkNAnodes];
  Float_t fSumBaseline[fgkNAnodes];
  Float_t fSumRawNoise[fgkNAnodes];
  Float_t fSumCMN[fgkNAnodes];
  Float_t fMinBaseline;
  Float_t fMaxBaseline;
  Float_t fMinRawNoise;
  Float_t fMaxRawNoise;
  Float_t fNSigmaNoise;

  ClassDef(AliITSOnlineSDDBase,1);
};
#endif
