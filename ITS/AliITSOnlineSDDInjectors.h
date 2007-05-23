#ifndef ALIITSONLINESDDINJECTORS_H
#define ALIITSONLINESDDINJECTORS_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD injector analysis                           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include"AliITSOnlineSDD.h"

class TH2F;
class TGraphErrors;
class AliITSOnlineSDDInjectors : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDInjectors();      
  AliITSOnlineSDDInjectors(Int_t mod, Int_t sid);
  virtual ~AliITSOnlineSDDInjectors();

  void SetSide(Int_t sid){fSide=sid;}
  void SetThreshold(Float_t thr=75.){fThreshold=thr;}
  void SetRangeLine1(Int_t tbmin=40, Int_t tbmax=90){
    fTbMin[0]=tbmin; fTbMax[0]=tbmax; 
  }
  void SetRangeLine2(Int_t tbmin=90, Int_t tbmax=140){
    fTbMin[1]=tbmin; fTbMax[1]=tbmax; 
  }
  void SetRangeLine3(Int_t tbmin=170, Int_t tbmax=220){
    fTbMin[2]=tbmin; fTbMax[2]=tbmax; 
  }
  void SetPolOrder(Int_t n=3){fPolOrder=n;}
  void SetMinDriftVel(Float_t vmin=4.){fMinDriftVel=vmin;}
  void SetMaxDriftVel(Float_t vmax=9.){fMaxDriftVel=vmax;}

  TGraphErrors* GetLineGraph(Int_t jlin);
  TGraphErrors* GetDriftVelocityGraph();
  Float_t* GetDriftVelFitParam()const{ return fParam;}
  Float_t GetDriftVelocity(Int_t jlin) const{return fDriftVel[jlin];}
  Float_t GetSigmaDriftVelocity(Int_t jlin) const{return fSigmaDriftVel[jlin];}
  Float_t GetTimeBinZero() const{return fTbZero;}
  Float_t GetDriftCoordinate(Float_t cAnode, Float_t cTimeBin);
  Int_t GetAnodeNumber(Int_t iInjLine);

  void PrintInjMap();
  void PrintCentroids();
  void WriteToFXS();

  void Reset();      
  void AnalyzeEvent(TH2F* his);      
  void FindGoodInjectors();
  void FindCentroids();
  void CalcDriftVelocity(Int_t jlin);
  void CalcTimeBinZero();
  void FitDriftVelocityVsAnode();

 protected:
  void SetPositions();
 private:

  enum {
    kNInjectors = 33
  };

  AliITSOnlineSDDInjectors(const AliITSOnlineSDDInjectors& source);
  AliITSOnlineSDDInjectors& operator = (const AliITSOnlineSDDInjectors& source);
  static const Float_t fgkSaturation;
  static const Float_t fgkJitterTB;

  TH2F* fHisto;
  Float_t fTbZero;
  Float_t fPosition[3];
  UShort_t fTbMin[3];
  UShort_t fTbMax[3];
  Bool_t fGoodInj[kNInjectors][3];
  Float_t fCentroid[kNInjectors][3];
  Float_t fRMSCentroid[kNInjectors][3];
  Float_t fDriftVel[kNInjectors];
  Float_t fSigmaDriftVel[kNInjectors];
  Float_t *fParam;
  Int_t fPolOrder;
  Float_t fMinDriftVel;
  Float_t fMaxDriftVel;
  Float_t fThreshold;

  ClassDef(AliITSOnlineSDDInjectors,1)
};
#endif
