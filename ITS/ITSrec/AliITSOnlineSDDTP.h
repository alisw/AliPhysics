#ifndef ALIITSONLINESDDTP_H
#define ALIITSONLINESDDTP_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD Test Pulse analysis                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSOnlineSDD.h"

class TH1F;
class TH2F;
class AliITSOnlineSDDTP : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDTP();
  AliITSOnlineSDDTP(Int_t nddl, Int_t ncarlos, Int_t sid, Float_t xDAC);
  virtual ~AliITSOnlineSDDTP();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();
  void ReadBaselines();

  void SetNSigmaGain(Float_t sig=3.){fNSigmaGain=sig;}
  void SetNSigmaNoise(Float_t sig=10.){fNSigmaNoise=sig;}
  Bool_t IsModuleGood()const;
  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Float_t GetAnodeBaseline(Int_t iAnode) const{ return fBaseline[iAnode];}
  Int_t GetAnodeEqualizedBaseline(Int_t iAnode) const{ return fEqBaseline[iAnode];}
  Int_t GetAnodeBaselineOffset(Int_t iAnode) const{ return fOffsetBaseline[iAnode];}
  Float_t GetAnodeRawNoise(Int_t iAnode) const{ return fRawNoise[iAnode];}
  Float_t GetAnodeCommonMode(Int_t iAnode) const{ return fCMN[iAnode];}
  Float_t GetAnodeCorrNoise(Int_t iAnode) const{return fCorrNoise[iAnode];}
  Float_t GetTimeBinTPPeak(Int_t iAnode) const{    
    if(fNEvents[iAnode]>0) return fTPPos[iAnode]/fNEvents[iAnode];
    else return 0;
  }
  Int_t GetNEvents(Int_t iAnode) const {return fNEvents[iAnode];}
  Float_t GetChannelGain(Int_t iAnode)const{
    if(fNEvents[iAnode]>0) return fSumTPPeak[iAnode]/fNEvents[iAnode]/fDAC;
    else return 0;
  }
  void StatGain(Float_t &mean, Float_t  &rms) const;
  void WriteToASCII();
  Bool_t WriteToROOT(TFile *fil);

  TH1F* GetBaselineAnodeHisto() const;
  TH1F* GetRawNoiseAnodeHisto() const;
  TH1F* GetCorrNoiseAnodeHisto() const;
  TH1F* GetCMNCoefAnodeHisto() const;
  TH1F* GetGainAnodeHisto() const;
  TH1F* GetStatusAnodeHisto() const;

 protected:

 private:
  Int_t fNEvents[fgkNAnodes];      // number of TP events for given anode
  Float_t fDAC;                    // Pascal Test Pulse amplitude (DAC units)
  Int_t fLowThreshold;           // low threshold for zero supp.
  Int_t fHighThreshold;          // high threshold for zero supp.
  Bool_t fGoodAnode[fgkNAnodes];   // array of anode quality (1 good, 0 bad) 
  Float_t fBaseline[fgkNAnodes];   // array of anode baselines
  Int_t fEqBaseline[fgkNAnodes];   // array of anode baselines after equalization
  Int_t fOffsetBaseline[fgkNAnodes]; // array of offsets for baseline equal.
  Float_t fRawNoise[fgkNAnodes];   // array of anode raw noise
  Float_t fCMN[fgkNAnodes];        // common mode noise coeff.
  Float_t fCorrNoise[fgkNAnodes];  // array of anode corrected noise
  Float_t fSumTPPeak[fgkNAnodes];  // test pulse amplitude summed over events
  Float_t fTPPos[fgkNAnodes];      // test pulse position
  Float_t fNSigmaGain;             // Cut value for gain (n*sigma)
  Float_t fNSigmaNoise;            // Threshold for TP signal identification
  ClassDef(AliITSOnlineSDDTP,4);
};
#endif
