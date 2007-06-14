#ifndef ALIITSONLINESDDBTP_H
#define ALIITSONLINESDDBTP_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD baseline, noise and gain analysis          //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSOnlineSDD.h"
#include "TMath.h"

class TH2F;
class TGraph;
class AliITSOnlineSDDBTP : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDBTP();
  AliITSOnlineSDDBTP(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDDBTP();
  void Reset();
  void AddBaseEvent(TH2F* hrawd);
  void AddTPEvent(TH2F* hrawd, Float_t xDAC);
  void ValidateAnodes();
  void SetMinBaseline(Float_t bas=10.){fMinBaseline=bas;}
  void SetMaxBaseline(Float_t bas=150.){fMaxBaseline=bas;}
  void SetMinRawNoise(Float_t ns=0.001){fMinRawNoise=ns;}
  void SetMaxRawNoise(Float_t ns=9.){fMaxRawNoise=ns;}
  void SetNSigmaNoise(Float_t ns=4.){fNSigmaNoise=ns;}
  void SetNSigmaGain(Float_t sig=3.){fNSigmaGain=sig;}
  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Float_t GetAnodeBaseline(Int_t iAnode) const{
    if(fNBaseEvents>0) return fSumBaseline[iAnode]/fNBaseEvents;
    else return 0;
  }
  Float_t GetAnodeRawNoise(Int_t iAnode) const{
    if(fNBaseEvents>0) return TMath::Sqrt(fSumRawNoise[iAnode]/fNBaseEvents-TMath::Power(GetAnodeBaseline(iAnode),2));
    
    else return 0;
  }
  Float_t CalcMeanRawNoise() const;
  void StatGain(Float_t &mean, Float_t  &rms);
  Float_t GetAnodeCommonMode(Int_t iAnode) const{
    if(fNBaseEvents>0) return fSumCMN[iAnode]/fNBaseEvents;
    else return 0;
  }
  Float_t GetChannelGain(Int_t iAnode)const{
    if(fNTPEvents>0) return fSumTPPeak[iAnode]/fNTPEvents;
    else return 0;
  }
  Int_t GetNBaseEvents() const {return fNBaseEvents;}
  Int_t GetNTPEvents() const {return fNTPEvents;}
  void WriteToASCII();

 protected:

 private:

  Int_t fNBaseEvents;                // number of "empty" events
  Int_t fNTPEvents;                  // number of "Test Pulse" events
  Bool_t fGoodAnode[fgkNAnodes];     // anode quality: good(1) - bad (0)
  Float_t fSumBaseline[fgkNAnodes];  // baseline summed over events
  Float_t fSumRawNoise[fgkNAnodes];  // noise summed over events
  Float_t fSumCMN[fgkNAnodes];       // common mode noise coeff.
  Float_t fSumTPPeak[fgkNAnodes];    // Test Pulse ampl. summed over events
  Float_t fTPPos[fgkNAnodes];        // Test pulse peak position
  Float_t fMinBaseline;              // Cut value for minimum baseline
  Float_t fMaxBaseline;              // Cut value for maximum baseline
  Float_t fMinRawNoise;              // Cut value for minimum noise
  Float_t fMaxRawNoise;              // Cut value for maximum noise
  Float_t fNSigmaNoise;              // Cut value for noise (n*sigma)
  Float_t fNSigmaGain;               // Cut value for gain (n*sigma)

  ClassDef(AliITSOnlineSDDBTP,1);
};
#endif
