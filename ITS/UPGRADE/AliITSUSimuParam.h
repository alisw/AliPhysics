#ifndef ALIITSUSIMUPARAM_H
#define ALIITSUSIMUPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store the parameters used in the simulation ITS      //
//                                                               //
///////////////////////////////////////////////////////////////////
#include <TRandom.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TF1.h>
#include "AliMathBase.h"
class AliITSUParamList;

class AliITSUSimuParam : public TObject {

 public:
  enum {kNoCouplingPix,kOldCouplingPix,kNewCouplingPix,kMaxCouplingOptPix};
  //
  AliITSUSimuParam();
  AliITSUSimuParam(UInt_t nLayers,UInt_t nPix);
  AliITSUSimuParam(const AliITSUSimuParam& simpar);
  // assignment operator 
  AliITSUSimuParam& operator=(const AliITSUSimuParam& source);
  ~AliITSUSimuParam();
  //
  void     SetNPix(Int_t np);
  void     SetNLayers(Int_t nl);
  //
  Double_t ApplyPixBaselineAndNoise(UInt_t mod) const;
  Double_t CalcProbNoiseOverThreshold(UInt_t mod) const;
  //
  void     SetLrROCycleShift(Double_t v=0, Int_t lr=-1);
  Double_t GetLrROCycleShift(Int_t lr)                                    const  {return fLrROCycleShift ? fLrROCycleShift[lr]:0;}
  //
  void     SetPixThreshold(Double_t thresh, Double_t sigma, int mod=-1);
  void     SetPixMinElToAdd(Double_t nel)                                        {fPixMinElToAddDef = nel;}
  void     GetPixThreshold(UInt_t mod, Double_t& thresh, Double_t& sigma) const;
  Double_t GetPixThreshold(UInt_t mod)                                    const;
  Double_t GetPixMinElToAdd()                                             const  {return fPixMinElToAddDef;}
  //
  void     SetPixFakeRate(Double_t frate)                                        { fPixFakeRateDef = frate;}
  Double_t GetPixFakeRate()                                               const  {return fPixFakeRateDef;}
  //
  void     SetPixNoiseInAllMod(Bool_t ison)                                      {fPixNoiseInAllMod = ison;}
  Bool_t   GetPixNoiseInAllMod()                                          const  {return fPixNoiseInAllMod;}
  //
  void     SetPixNoise(Double_t noise, Double_t baseline, Int_t mod=-1);
  void     GetPixNoise(UInt_t mod,Double_t &noise, Double_t &baseline)    const;
  //
  void     SetPixBiasVoltage(Double_t bias=18.182,Int_t mod=-1);
  Double_t GetPixBiasVoltage(UInt_t mod)                                  const;
  // 
  void     SetGeVToCharge(Double_t gc=fgkNcompsDefault)                          {fGeVcharge = gc;}
  Double_t GetGeVToCharge()                                               const  {return fGeVcharge;}
  Double_t GeVToCharge(Double_t gev)                                      const  {return gev/fGeVcharge;}
  //
  void     SetPixCouplingOption(UInt_t opt);
  UInt_t   GetPixCouplingOption()                                         const  {return fPixCouplOpt;}

  void     SetPixCouplingParam(Double_t col, Double_t row)                       {fPixCouplCol = col; fPixCouplRow = row;}
  void     GetPixCouplingParam(Double_t &col, Double_t &row)              const  {col = fPixCouplCol; row = fPixCouplRow;}

  void     SetPixLorentzDrift(Bool_t ison)                                       {fPixLorentzDrift=ison;}
  Bool_t   GetPixLorentzDrift()                                           const  {return fPixLorentzDrift;}
  void     SetPixLorentzHoleWeight(Double_t weight)                               {fPixLorentzHoleWeight=weight;}
  Double_t GetPixLorentzHoleWeight()                                      const  {return fPixLorentzHoleWeight;}
  
  void     SetPixAddNoisyFlag(Bool_t value)                                      {fPixAddNoisyFlag = value;}
  Bool_t   GetPixAddNoisyFlag()                                           const  {return fPixAddNoisyFlag;}
  void     SetPixRemoveDeadFlag(Bool_t value)                                    {fPixRemoveDeadFlag = value;}
  Bool_t   GetPixRemoveDeadFlag()                                         const  {return fPixRemoveDeadFlag;}
  //
  Double_t LorentzAngleElectron(Double_t bz)                              const;
  Double_t LorentzAngleHole(Double_t bz)                                  const;
  //
  Int_t    GetNRespFunParams()                                            const  {return fRespFunParam.GetEntriesFast();}
  const    AliITSUParamList* GetRespFunParams(Int_t i)                    const  {return (const AliITSUParamList*)fRespFunParam[i];}
  const    AliITSUParamList* FindRespFunParams(Int_t detId)               const;
  void     AddRespFunParam(AliITSUParamList* pr);
  //
  virtual void Print(Option_t *opt = "")                                  const; 
  //
  static Double_t CalcProbNoiseOverThreshold(double base, double noise, double thresh);
  static Double_t GenerateNoiseQFunction(double prob, double mean, double sigma);
  //
 protected:

  static const Float_t fgkPixBiasVoltageDefault;//default for fPixBiasVoltage
  static const Float_t fgkPixThreshDefault; //default for fThresh
  static const Float_t fgkPixMinElToAddDefault; // default min number of electrons to add to sdigit
  static const Float_t fgkPixThrSigmaDefault; //default for fSigma
  static const Float_t fgkPixCouplColDefault; //default for fPixCouplCol
  static const Float_t fgkPixCouplRowDefault; //default for fPixCouplRow
  static const Float_t fgkPixEccDiffDefault;//default for fPixEccDiff
  static const Float_t fgkPixLorentzHoleWeightDefault;//default for fPixLorentzHoleWeight
  static const UInt_t  fgkPixCouplingOptDefault;  // type of pixel Coupling (old or new)
  static const Float_t fgkDOverVDefault;             // default distance over voltage 
  static const Float_t fgkGeVtoChargeDefault;        // default energy to ionize (free an electron) in GeV
  static const Float_t fgkTDefault;                  // default temperature
  static const Float_t fgkPixFakeRateDefault;        // default monopix fake rate
  static const Bool_t  fgkPixNoiseInAllMod;          // default switch to add digital noise to every module
  
  static const Float_t fgkNsigmasDefault; //default for fNsigmas
  static const Int_t   fgkNcompsDefault; //default for fNcomps

 private:
  //
  Float_t    fGeVcharge;          // Energy to ionize (free an electron) in GeV
  Float_t    fDOverV;             // The parameter d/v where d is the disance over which the the potential v is applied d/v [cm/volts]
  Float_t    fT;                  // The temperature of the Si in Degree K.
  //
  UInt_t     fNLayers;             // number of layers  
  UInt_t     fNPix;                // number of Pix type detectors
  UInt_t     fPixCouplOpt;         // Pix Coupling Option
  Float_t    fPixCouplCol;         // Pix Coupling parameter along the cols
  Float_t    fPixCouplRow;         // Pix Coupling parameter along the rows
  Bool_t     fPixLorentzDrift;     // Flag to decide whether to simulate the Lorentz Drift or not in Pix
  Float_t    fPixLorentzHoleWeight;// Lorentz Angle is computed for Pix as average of Hole and Electron
  //                                    this parameter gives the relative weights between the two
  Bool_t     fPixAddNoisyFlag;     // Flag saying whether noisy pixels should be added to digits
  Bool_t     fPixRemoveDeadFlag;   // Flag saying whether dead pixels should be removed from digits
  //
  Float_t    fPixThreshDef;        // Pix Threshold value
  Float_t    fPixThrSigmaDef;      // Pix Threshold fluctuation
  Float_t    fPixBiasVoltageDef;   // Bias Voltage for the Pix
  Float_t    fPixNoiseDef;         // Pix electronic noise: sigma
  Float_t    fPixBaselineDef;      // Pix electronic noise: baseline
  Float_t    fPixMinElToAddDef;    // min number of electrons to add
  Float_t    fPixFakeRateDef;      // Fake rate for the monopix
  Bool_t     fPixNoiseInAllMod;    // Add digital noise to every module
  //
  Float_t*   fLrROCycleShift; //[fNLayers] optional RO cycle shift for each layer (in fraction of RO cycle length if abs<1)
                              // if abs>1, then each module will have random phase
  Float_t*   fPixThresh;      //[fNPix] Pix Threshold value
  Float_t*   fPixThrSigma;    //[fNPix] Pix Threshold fluctuation
  Float_t*   fPixBiasVoltage; //[fNPix] Bias Voltage for the Pix
  Float_t*   fPixSigma;       //[fNPix] Pix threshold fluctuations spread
  Float_t*   fPixNoise;       //[fNPix] Pix electronic noise: sigma
  Float_t*   fPixBaseline;    //[fNPix] Pix electronic noise: baseline
  //
  TObjArray  fRespFunParam;   // set of parameterizations for response function (AliITSUParamList)

  ClassDef(AliITSUSimuParam,2);  // ITSU simulataion params
};


#endif
