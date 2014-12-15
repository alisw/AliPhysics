#ifndef ALIITSONLINESDDBASE_H
#define ALIITSONLINESDDBASE_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD baseline and noise analysis                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSOnlineSDD.h"
#include <TMath.h>

class TH2F;
class TGraph;
class AliITSOnlineSDDBase : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDBase();
  AliITSOnlineSDDBase(Int_t nddl, Int_t ncarlos, Int_t sid);
  virtual ~AliITSOnlineSDDBase();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();

  void SetMinBaseline(Float_t bas=10.){fMinBaseline=bas;}
  void SetMaxBaseline(Float_t bas=150.){fMaxBaseline=bas;}
  void SetMinRawNoise(Float_t ns=0.5){fMinRawNoise=ns;}
  void SetMaxRawNoise(Float_t ns=9.){fMaxRawNoise=ns;}
  void SetNSigmaNoise(Float_t ns=4.){fNSigmaNoise=ns;}
  void SetGoldenBaselineValue(Float_t val=20.){fGoldenBaseline=val;}
  void SetZeroSuppThresholds(Float_t vall=2.2,Float_t valh=4.){
    fLowThrFact=vall;
    fHighThrFact=valh;
  }

  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Float_t GetAnodeBaseline(Int_t iAnode) const{
    if(fNEvents>0) return fSumBaseline[iAnode]/fNEvents;
    else return 0;
  }
  void GetMinAndMaxBaseline(Float_t &basMin, Float_t &basMax) const;
  Float_t GetMinimumBaseline() const;
  Float_t GetAnodeRawNoise(Int_t iAnode) const;

  Int_t CountGoodAnodes() const{
    Int_t nGdAn=0;
    for(Int_t ian=0;ian<fgkNAnodes;ian++) if(fGoodAnode[ian]) nGdAn++;  
    return nGdAn;
  }
  Float_t CalcMeanRawNoise() const;
  Float_t GetAnodeCommonMode(Int_t iAnode) const{
    if(fNEvents>0) return fSumCMN[iAnode]/fNEvents;
    else return 0;
  }
  Int_t GetNEvents() const {return fNEvents;}
  void WriteToASCII();

 protected:

 private:
  static const Int_t fgkMaxCorr;     // maximum baseline correction in AMBRA (=63)

  Int_t fNEvents;                    // number of events
  Bool_t fGoodAnode[fgkNAnodes];     // anode quality: good(1) - bad (0)
  Float_t fSumBaseline[fgkNAnodes];  // baseline summed over events
  Float_t fSumRawNoise[fgkNAnodes];  // noise summed over events
  Float_t fSumCMN[fgkNAnodes];       // common mode noise coeff.
  Float_t fMinBaseline;              // Cut value for minimum baseline
  Float_t fMaxBaseline;              // Cut value for maximum baseline
  Float_t fMinRawNoise;              // Cut value for minimum noise
  Float_t fMaxRawNoise;              // Cut value for maximum noise
  Float_t fNSigmaNoise;              // Cut value for noise (n*sigma)
  Float_t fGoldenBaseline;           // golden value for equalizing baselines
  Float_t fLowThrFact;               // factor for low threshold
  Float_t fHighThrFact;              // factor for high threshold

  ClassDef(AliITSOnlineSDDBase,2);
};

inline Float_t AliITSOnlineSDDBase::GetAnodeRawNoise(Int_t iAnode) const{
  // compute raw noise for given anode
  Float_t noise2=0.;
  if(fNEvents>0) noise2=fSumRawNoise[iAnode]/fNEvents-fSumBaseline[iAnode]*fSumBaseline[iAnode]/fNEvents/fNEvents;
  if(noise2>0.) return TMath::Sqrt(noise2);
  else return 0;
}

#endif
