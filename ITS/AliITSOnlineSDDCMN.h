#ifndef ALIITSONLINESDDCMN_H
#define ALIITSONLINESDDCMN_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD noise corrected for common mode analysis   //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSOnlineSDD.h"

class TH1F;
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

  Float_t CalcMeanNoise() const;
  Int_t GetNEvents() const {return fNEvents;}
  
  TH1F* GetBaselineAnodeHisto() const;
  TH1F* GetRawNoiseAnodeHisto() const;
  TH1F* GetCorrNoiseAnodeHisto() const;
  TH1F* GetBaselineHisto() const;
  TH1F* GetRawNoiseHisto() const;
  TH1F* GetCorrNoiseHisto() const;

  void WriteToASCII();
  Bool_t WriteToROOT(TFile *fil);

 protected:

 private:
  Int_t fNEvents;                    // number of events
  Bool_t fGoodAnode[fgkNAnodes];     // anode quality: good(1) - bad (0)
  Float_t fBaseline[fgkNAnodes];     // array of anode baselines
  Float_t fRawNoise[fgkNAnodes];     // array of anode raw noise
  Float_t fSumCorrNoise[fgkNAnodes]; // corrected noise summed over events
  Float_t fCMN[fgkNAnodes];          // common mode noise coeff.
  Float_t fMinCorrNoise;             // Cut value for minimum corrected noise
  Float_t fMaxCorrNoise;             // Cut value for maximum corrected noise
  Float_t fNSigmaNoise;              // Cut value for corrected noise (n*sigma)

  ClassDef(AliITSOnlineSDDCMN,1);
};
#endif
