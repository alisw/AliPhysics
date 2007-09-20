#ifndef ALIITSONLINESDDINJECTORS_H
#define ALIITSONLINESDDINJECTORS_H


///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD injector analysis                           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSOnlineSDD.h"

/* $Id$ */
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
  TGraphErrors* GetDriftVelocityGraph() const;
  Float_t* GetDriftVelFitParam()const{ return fParam;}
  Float_t GetDriftVelocity(Int_t jlin) const{return fDriftVel[jlin];}
  Float_t GetSigmaDriftVelocity(Int_t jlin) const{return fSigmaDriftVel[jlin];}
  Float_t GetTimeBinZero() const{return fTbZero;}
  Float_t GetDriftCoordinate(Float_t cAnode, Float_t cTimeBin);
  Int_t GetAnodeNumber(Int_t iInjLine) const;
  Float_t GetCentroid(Int_t injnumb, Int_t injline) const {
    if(injnumb<kNInjectors && injline<3) return fCentroid[injnumb][injline];
    else return -9999.;
  }
  Bool_t IsInjectorGood(Int_t injnumb, Int_t injline) const {
    if(injnumb<kNInjectors && injline<3) return fGoodInj[injnumb][injline];
    else return 0;
  }
  void PrintInjMap();
  void PrintCentroids();
  void WriteToASCII(Int_t evNumb, UInt_t timeStamp, Int_t optAppend=0);

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

  TH2F* fHisto;                         // histogram of module channel counts
  Float_t fTbZero;                      // Time zero for injector event
  Float_t fPosition[3];                 // Coordinates of injector lines
  UShort_t fTbMin[3];                   // Minimum time bin for each line
  UShort_t fTbMax[3];                   // Maximum time bin for each line
  Bool_t fGoodInj[kNInjectors][3];      // array of good injectors
  Float_t fCentroid[kNInjectors][3];    // array of time centroids of injectors
  Float_t fRMSCentroid[kNInjectors][3]; // array of time rms of injectors
  Float_t fDriftVel[kNInjectors];       // drift velocity
  Float_t fSigmaDriftVel[kNInjectors];  // error on drift velocity
  Float_t *fParam;                      // parameters of polinomial fit
                                        //  of drift vel. vs. anode number
  Int_t fPolOrder;                      // order of polinomial fit
  Float_t fMinDriftVel;                 // Cut value for minimum drift speed
  Float_t fMaxDriftVel;                 // Cut value for maximum drift speed
  Float_t fThreshold;                   // Threshold for injector signal

  ClassDef(AliITSOnlineSDDInjectors,1)
};
#endif
