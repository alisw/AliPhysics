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
  void SetTimeStep(Float_t tstep) {
    fTimeStep=tstep;
  }

  TGraphErrors* GetTimeVsDistGraph(Int_t jpad) const;
  TGraphErrors* GetDriftSpeedGraph() const;
  TGraphErrors* GetSelectedDriftSpeedGraph(Int_t minAcceptStatus) const;
  Float_t* GetDriftSpeedFitParam()const{ return fParam;}
  Float_t GetDriftSpeed(Int_t jpad) const{return fDriftSpeed[jpad];}
  Float_t GetDriftSpeedErr(Int_t jpad) const{return fDriftSpeedErr[jpad];}
  Float_t GetTimeBinZero() const{return fTbZero;}

  Float_t GetTimeStep() const{return fTimeStep;}
  Int_t GetAnodeNumber(Int_t iInjPad) const;
  Int_t GetInjPadNumberFromAnode(Int_t nAnode) const;
  Int_t GetInjPadStatus(Int_t jpad) const;  
  Int_t GetAnodeStatus(Int_t nAnode) const{
    Int_t jpad=GetInjPadNumberFromAnode(nAnode);
    return GetInjPadStatus(jpad);
  }  
  Float_t GetCentroid(Int_t jpad, Int_t jlin) const {
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
  void FindGoodInjectors();
  void FindCentroids();
  void CalcDriftSpeed(Int_t jpad);
  void CalcTimeBinZero();
  void FitDriftSpeedVsAnode();

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
  Float_t fTbZero;                           // Time zero for injector event
  Float_t fPosition[kInjLines];              // Coordinates of injector lines
  UShort_t fTbMin[kInjLines];                // Minimum time bin for each line
  UShort_t fTbMax[kInjLines];                // Maximum time bin for each line
  Bool_t fGoodInj[kInjPads][kInjLines];      // array of good injectors
  Float_t fCentroid[kInjPads][kInjLines];    // array of time bin centroids
  Float_t fRMSCentroid[kInjPads][kInjLines]; // array of time rms of injectors
  Float_t fDriftSpeed[kInjPads];             // drift speed
  Float_t fDriftSpeedErr[kInjPads];          // error on drift speed
  Float_t *fParam;                           // parameters of polinomial fit to
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

  Float_t fTimeStep;                 // time bin value (25 or 50 ns)

  ClassDef(AliITSOnlineSDDInjectors,4)
};
#endif
