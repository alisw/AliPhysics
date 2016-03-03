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
  AliITSOnlineSDDCMN(Int_t nddl, Int_t ncarlos, Int_t sid);
  virtual ~AliITSOnlineSDDCMN();
  void Reset();
  TH2F* GetCleanEvent(const TH2F* hrawd) const;
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();
  void ValidateModule();
  void ReadBaselines();

  void SetMinNoise(Float_t ns=0.001){fMinCorrNoise=ns;}
  void SetMaxNoise(Float_t ns=9.){fMaxCorrNoise=ns;}
  void SetNSigmaNoise(Float_t ns=4.){fNSigmaNoise=ns;}
  void SetMinNumGoodAnForGoodMod(Int_t cut=16){fMinNumGoodAnForGoodMod=cut;}
  void SetCutOnGoodAnodeClusterSize(Int_t clu=5){fMinClusterOfGoodAn=clu;}

  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Int_t CountGoodAnodes() const{
    Int_t nGdAn=0;
    for(Int_t ian=0;ian<fgkNAnodes;ian++) if(fGoodAnode[ian]) nGdAn++;  
    return nGdAn;
  }

  Float_t GetAnodeBaseline(Int_t iAnode) const{ return fBaseline[iAnode];}
  Int_t GetAnodeEqualizedBaseline(Int_t iAnode) const{ return fEqBaseline[iAnode];}
  Int_t GetAnodeBaselineOffset(Int_t iAnode) const{ return fOffsetBaseline[iAnode];}
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
  TH1F* GetCMNCoefAnodeHisto() const;
  TH1F* GetStatusAnodeHisto() const;
  TH1F* GetBaselineHisto() const;
  TH1F* GetRawNoiseHisto() const;
  TH1F* GetCorrNoiseHisto() const;

  void WriteToASCII();
  Bool_t WriteToROOT(TFile *fil);

 protected:

 private:
  Int_t fNEvents;                    // number of events
  Int_t fLowThreshold;             // low threshold for zero supp.
  Int_t fHighThreshold;            // high threshold for zero supp.
  Bool_t fGoodAnode[fgkNAnodes];     // anode quality: good(1) - bad (0)
  Float_t fBaseline[fgkNAnodes];     // array of anode baselines
  Int_t fEqBaseline[fgkNAnodes];     // array of anode baselines after equalization
  Int_t fOffsetBaseline[fgkNAnodes]; // array of offsets for baseline equal.
  Float_t fRawNoise[fgkNAnodes];     // array of anode raw noise
  Float_t fSumCorrNoise[fgkNAnodes]; // corrected noise summed over events
  Float_t fCMN[fgkNAnodes];          // common mode noise coeff.
  Float_t fMinCorrNoise;             // Cut value for minimum corrected noise
  Float_t fMaxCorrNoise;             // Cut value for maximum corrected noise
  Float_t fNSigmaNoise;              // Cut value for corrected noise (n*sigma)
  Int_t   fMinNumGoodAnForGoodMod;   // Min. n. good anodes to tag mod as good
  Int_t   fMinClusterOfGoodAn;       // Min. n. of adjacent good anodes

  ClassDef(AliITSOnlineSDDCMN,4);
};
#endif
