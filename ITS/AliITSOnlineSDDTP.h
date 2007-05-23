#ifndef ALIITSONLINESDDTP_H
#define ALIITSONLINESDDTP_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD Test Pulse analysis                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include"AliITSOnlineSDD.h"

class TH2F;
class TGraph;
class AliITSOnlineSDDTP : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDTP();
  AliITSOnlineSDDTP(Int_t mod, Int_t sid,Float_t xDAQ);
  virtual ~AliITSOnlineSDDTP();
  void Reset();
  void AddEvent(TH2F* hrawd);
  void ValidateAnodes();
  void ReadBaselines();

  void SetNSigmaGain(Float_t sig=3.){fNSigmaGain=sig;}
  Bool_t IsAnodeGood(Int_t iAnode)const{ return fGoodAnode[iAnode];}
  Int_t GetNEvents() const {return fNEvents;}
  Float_t GetChannelGain(Int_t iAnode)const{
    if(fNEvents>0) return fSumTPPeak[iAnode]/fNEvents/fDAQ;
    else return 0;
  }
  void StatGain(Float_t &mean, Float_t  &rms);
  void WriteToFXS();

 protected:

 private:
  Int_t fNEvents;
  Float_t fDAQ;
  Bool_t fGoodAnode[fgkNAnodes];
  Float_t fBaseline[fgkNAnodes];
  Float_t fSumTPPeak[fgkNAnodes];
  Float_t fTPPos[fgkNAnodes];
  Float_t fNSigmaGain;

  ClassDef(AliITSOnlineSDDTP,1);
};
#endif
