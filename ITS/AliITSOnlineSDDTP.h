#ifndef ALIITSONLINESDDTP_H
#define ALIITSONLINESDDTP_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD Test Pulse analysis                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSOnlineSDD.h"

class TH2F;
class AliITSOnlineSDDTP : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDTP();
  AliITSOnlineSDDTP(Int_t mod, Int_t sid,Float_t xDAC);
  virtual ~AliITSOnlineSDDTP();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();
  void ReadBaselines();

  void SetNSigmaGain(Float_t sig=3.){fNSigmaGain=sig;}
  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Int_t GetNEvents() const {return fNEvents;}
  Float_t GetChannelGain(Int_t iAnode)const{
    if(fNEvents>0) return fSumTPPeak[iAnode]/fNEvents/fDAC;
    else return 0;
  }
  void StatGain(Float_t &mean, Float_t  &rms);
  void WriteToFXS();

 protected:

 private:
  Int_t fNEvents;                  // number of events
  Float_t fDAC;                     // Pascal Test Pulse amplitude (DAC units)
  Bool_t fGoodAnode[fgkNAnodes]; // array of anode quality (1 good, 0 bad) 
  Float_t fBaseline[fgkNAnodes];   // array of anode baselines
  Float_t fSumTPPeak[fgkNAnodes];  // test pulse amplitude summed over events
  Float_t fTPPos[fgkNAnodes];      // test pulse position
  Float_t fNSigmaGain;             // Cut value for gain (n*sigma)

  ClassDef(AliITSOnlineSDDTP,1);
};
#endif
