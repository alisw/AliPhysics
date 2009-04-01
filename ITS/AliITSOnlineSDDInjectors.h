#ifndef ALIITSONLINESDDINJECTORS_H
#define ALIITSONLINESDDINJECTORS_H


/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used for SDD injector analysis                           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSOnlineSDD.h"


class TH2F;
class TGraphErrors;
class AliITSOnlineSDDInjectors : public AliITSOnlineSDD {

 public:
  AliITSOnlineSDDInjectors();      
  AliITSOnlineSDDInjectors(Int_t nddl, Int_t ncarlos, Int_t sid);
  virtual ~AliITSOnlineSDDInjectors();

  void SetThresholds(Float_t tl, Float_t th){
    fLowThreshold=tl;
    fHighThreshold=th;
  }
  void SetInjLineRange(Int_t jlin, Int_t tbmin, Int_t tbmax){
    fTbMin[jlin]=tbmin;
    fTbMax[jlin]=tbmax;
  }
  void SetPolOrder(Int_t n){fPolOrder=n;}
  void SetMinDriftSpeed(Float_t vmin){fMinDriftSpeed=vmin;}
  void SetMaxDriftSpeed(Float_t vmax){fMaxDriftSpeed=vmax;}
  void SetMaxDriftSpeedErr(Float_t maxval){
    fMaxDriftSpeedErr=maxval;
  }
  void SetFitLimits(Int_t firstpad,Int_t lastpad){
    fFirstPadForFit=firstpad;
    fLastPadForFit=lastpad;
  }
  void SetPadStatusCutForFit(Int_t cutval=1){
    fPadStatusCutForFit=cutval;
  }
  void SetDefaults();
  void SetTimeStep(Double_t tstep) {
    fTimeStep=tstep;
  }
  void SetUseTimeZeroSignal(Bool_t useTZ=kTRUE){
    fUseTimeZeroSignal=useTZ;
  }
  TGraphErrors* GetTimeVsDistGraph(Int_t jpad) const;
  TGraphErrors* GetDriftSpeedGraph() const;
  TGraphErrors* GetSelectedDriftSpeedGraph(Int_t minAcceptStatus) const;
  Double_t* GetDriftSpeedFitParam()const{ return fParam;}
  Double_t GetDriftSpeed(Int_t jpad) const{return fDriftSpeed[jpad];}
  Double_t GetDriftSpeedErr(Int_t jpad) const{return fDriftSpeedErr[jpad];}
  Double_t GetTimeBinZero() const{return fTbZero;}

  Double_t GetTimeStep() const{return fTimeStep;}
  Int_t GetAnodeNumber(Int_t iInjPad) const;
  Int_t GetInjPadNumberFromAnode(Int_t nAnode) const;
  Int_t GetInjPadStatus(Int_t jpad) const;  
  Int_t GetAnodeStatus(Int_t nAnode) const{
    Int_t jpad=GetInjPadNumberFromAnode(nAnode);
    return GetInjPadStatus(jpad);
  }  
  Double_t GetCentroid(Int_t jpad, Int_t jlin) const {
    if(jpad<kInjPads && jlin<kInjLines) return fCentroid[jpad][jlin];
    else return -9999.;
  }
  Bool_t IsInjectorGood(Int_t jpad, Int_t jlin) const {
    if(jpad<kInjPads && jlin<kInjLines) return fGoodInj[jpad][jlin];
    else return 0;
  }
  void PrintInjectorStatus();
  void PrintCentroids();
  void WriteToASCII(Int_t evNumb, UInt_t timeStamp, Int_t optAppend=0);

  void Reset();
  void AnalyzeEvent(TH2F* his);
  void AddEvent(TH2F* his);
  void FindGoodInjectors();
  void FindCentroids();
  void CalcDriftSpeed(Int_t jpad);
  void CalcTimeBinZero();
  void FitDriftSpeedVsAnode();
  Double_t GetMeanDriftSpeed(Int_t ipad) const{
    if(fNEvents==0) return 0.;
    return fSumDriftSpeed[ipad]/(Double_t)fNEvents;
  }
  Double_t GetRMSDriftSpeed(Int_t ipad) const;
  void FitMeanDriftSpeedVsAnode();

 protected:
  void SetPositions();
 private:

  enum {kInjPads  = 33};
  enum {kInjLines = 3};

  AliITSOnlineSDDInjectors(const AliITSOnlineSDDInjectors& source);
  AliITSOnlineSDDInjectors& operator = (const AliITSOnlineSDDInjectors& source);
  static const Float_t fgkSaturation;        // ADC saturation value (1008)
  static const Float_t fgkDefaultLThreshold;  // Default for fLowThreshold
  static const Float_t fgkDefaultHThreshold;  // Default for fHighThreshold
  static const Float_t fgkDefaultMinSpeed;   // Default for fMinDriftSpeed
  static const Float_t fgkDefaultMaxSpeed;   // Default for fMaxDriftSpeed
  static const Float_t fgkDefaultMaxErr;     // Default for fMaxDriftSpeedErr
  static const Int_t   fgkDefaultPolOrder;   // Default for fPolOrder
  static const Float_t fgkDefaultTimeStep;   // Default for fTimeStep
  static const UShort_t   fgkDefaultTbMin[kInjLines];  // Defaults for fTbMin
  static const UShort_t   fgkDefaultTbMax[kInjLines];  // Defaults for fTbMax


  TH2F* fHisto;                              // histogram of channel counts
  Double_t fTbZero;                           // Time zero for injector event
  Double_t fRMSTbZero;                        // Error on time zero 
  Double_t fPosition[kInjLines];              // Coordinates of injector lines
  UShort_t fTbMin[kInjLines];                // Minimum time bin for each line
  UShort_t fTbMax[kInjLines];                // Maximum time bin for each line
  Bool_t fGoodInj[kInjPads][kInjLines];      // array of good injectors
  Double_t fCentroid[kInjPads][kInjLines];    // array of time bin centroids
  Double_t fRMSCentroid[kInjPads][kInjLines]; // array of time rms of injectors
  Double_t fDriftSpeed[kInjPads];             // drift speed  
  Double_t fDriftSpeedErr[kInjPads];          // error on drift speed
  Int_t fNEvents;                            // number of events
  Double_t fSumDriftSpeed[kInjPads];          // drift speed summed over events  
  Double_t fSumSqDriftSpeed[kInjPads];        // drift speed^2 sum
  
  Double_t *fParam;                           // parameters of polinomial fit to
                                             // drift speed vs. anode number
  Int_t fPolOrder;                   // order of polinomial fit
  Float_t fMinDriftSpeed;            // Minimum value for drift speed
  Float_t fMaxDriftSpeed;            // Maximum value for drift speed
  Float_t fMaxDriftSpeedErr;         // Maximum value for error on drift speed
  Float_t fLowThreshold;             // Low threshold for injector signal
  Float_t fHighThreshold;            // High threshold for injector signal

  Int_t fFirstPadForFit;             // first injector pad used in fit
  Int_t fLastPadForFit;              // last injector pad used in fit
  Int_t fPadStatusCutForFit;         // minimum value of pad status for fit

  Double_t fTimeStep;                 // time bin value (25 or 50 ns)
  Bool_t fUseTimeZeroSignal;         // flag for usage of time zero signal
                                     // in drift speed calculation

  ClassDef(AliITSOnlineSDDInjectors,6)
};
#endif
