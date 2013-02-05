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
#include "AliMathBase.h"
class AliParamList;

class AliITSUSimuParam : public TObject {

 public:
  enum {kOldCouplingPix,kNewCouplingPix,kMaxCouplingOptPix};
  //
  AliITSUSimuParam();
  AliITSUSimuParam(UInt_t nPix);
  AliITSUSimuParam(const AliITSUSimuParam& simpar);
  // assignment operator 
  AliITSUSimuParam& operator=(const AliITSUSimuParam& source);
  ~AliITSUSimuParam();

  Double_t ApplyPixBaselineAndNoise(UInt_t mod) const;
  Double_t CalcProbNoiseOverThreshold(UInt_t mod) const;
  //
  void     SetPixThreshold(Double_t thresh, Double_t sigma, int mod=-1);
  void     SetPixMinElToAdd(Double_t nel)                                        {fPixMinElToAddDef = nel;}
  void     GetPixThreshold(UInt_t mod, Double_t& thresh, Double_t& sigma) const;
  Double_t GetPixThreshold(UInt_t mod)                                    const;
  Double_t GetPixMinElToAdd()                                             const  {return fPixMinElToAddDef;}
  //
  void     SetPixNoise(Double_t noise, Double_t baseline, Int_t mod=-1);
  void     GetPixNoise(UInt_t mod,Double_t &noise, Double_t &baseline)    const;
  //
  void     SetPixBiasVoltage(Double_t bias=18.182,Int_t mod=-1);
  Double_t GetPixBiasVoltage(UInt_t mod)                                  const;

 
  void     SetGeVToCharge(Double_t gc=fgkNcompsDefault)                             {fGeVcharge = gc;}
  Double_t GetGeVToCharge()                                                  const  {return fGeVcharge;}
  Double_t GeVToCharge(Double_t gev)                                         const  {return gev/fGeVcharge;}
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
  Double_t LorentzAngleElectron(Double_t bz)                                 const;
  Double_t LorentzAngleHole(Double_t bz)                                     const;
  //
  Int_t    GetNRespFunParams()                                               const {return fRespFunParam.GetEntriesFast();}
  const AliParamList* GetRespFunParams(Int_t i)                              const {return (const AliParamList*)fRespFunParam[i];}
  const AliParamList* FindRespFunParams(Int_t detId)                         const;
  void     AddRespFunParam(AliParamList* pr);
  //
  virtual void Print(Option_t *opt = "")                                     const; 
  //
  static Double_t CalcProbNoiseOverThreshold(double base, double noise, double thresh);
  static Double_t GenerateNoiseQFunction(double prob, double mean, double sigma);
  //
 protected:

  static const Double_t fgkPixBiasVoltageDefault;//default for fPixBiasVoltage
  static const Double_t fgkPixThreshDefault; //default for fThresh
  static const Double_t fgkPixMinElToAddDefault; // default min number of electrons to add to sdigit
  static const Double_t fgkPixThrSigmaDefault; //default for fSigma
  static const Double_t fgkPixCouplColDefault; //default for fPixCouplCol
  static const Double_t fgkPixCouplRowDefault; //default for fPixCouplRow
  static const Double_t fgkPixEccDiffDefault;//default for fPixEccDiff
  static const Double_t fgkPixLorentzHoleWeightDefault;//default for fPixLorentzHoleWeight
  static const UInt_t   fgkPixCouplingOptDefault;  // type of pixel Coupling (old or new)
  static const Double_t fgkDOverVDefault;             // default distance over voltage 
  static const Double_t fgkGeVtoChargeDefault;        // default energy to ionize (free an electron) in GeV
  static const Double_t fgkTDefault;                  // default temperature

  static const Double_t fgkNsigmasDefault; //default for fNsigmas
  static const Int_t fgkNcompsDefault; //default for fNcomps

 private:
  //
  Double_t   fGeVcharge;          // Energy to ionize (free an electron) in GeV
  Double_t   fDOverV;             // The parameter d/v where d is the disance over which the the potential v is applied d/v [cm/volts]
  Double_t   fT;                  // The temperature of the Si in Degree K.
  //
  UInt_t     fNPix;            // number of Pix type detectors
  UInt_t     fPixCouplOpt;     // Pix Coupling Option
  Double_t   fPixCouplCol;     // Pix Coupling parameter along the cols
  Double_t   fPixCouplRow;     // Pix Coupling parameter along the rows
  Bool_t     fPixLorentzDrift;     // Flag to decide whether to simulate the Lorentz Drift or not in Pix
  Double_t   fPixLorentzHoleWeight;// Lorentz Angle is computed for Pix as average of Hole and Electron
  //                                    this parameter gives the relative weights between the two
  Bool_t     fPixAddNoisyFlag;     // Flag saying whether noisy pixels should be added to digits
  Bool_t     fPixRemoveDeadFlag;   // Flag saying whether dead pixels should be removed from digits
  //
  Double_t   fPixThreshDef;      // Pix Threshold value
  Double_t   fPixThrSigmaDef;    // Pix Threshold fluctuation
  Double_t   fPixBiasVoltageDef; // Bias Voltage for the Pix
  Double_t   fPixNoiseDef;       // Pix electronic noise: sigma
  Double_t   fPixBaselineDef;    // Pix electronic noise: baseline
  Double_t   fPixMinElToAddDef;  // min number of electrons to add
  //
  Double_t*  fPixThresh;      //[fNPix] Pix Threshold value
  Double_t*  fPixThrSigma;    //[fNPix] Pix Threshold fluctuation
  Double_t*  fPixBiasVoltage; //[fNPix] Bias Voltage for the Pix
  Double_t*  fPixSigma;       //[fNPix] Pix threshold fluctuations spread
  Double_t*  fPixNoise;       //[fNPix] Pix electronic noise: sigma
  Double_t*  fPixBaseline;    //[fNPix] Pix electronic noise: baseline
  //
  TObjArray  fRespFunParam;   // set of parameterizations for response function (AliParamList)

  ClassDef(AliITSUSimuParam,1);  // ITSU simulataion params
};

//_______________________________________________________________________
inline Double_t AliITSUSimuParam::CalcProbNoiseOverThreshold(double mean, double sigma, double thresh) 
{
  // calculate probability of gaussian noise exceeding the threshold
  if (mean+6*sigma<thresh) return 0;
  if (mean-6*sigma>thresh) return 1.;
  const double ksqrt2 = 1.41421356237309515e+00;
  return 0.5*AliMathBase::ErfcFast( (thresh-mean)/(sigma*ksqrt2));
}

//_______________________________________________________________________
inline Double_t AliITSUSimuParam::GenerateNoiseQFunction(double prob, double mean, double sigma) 
{
  // generate random noise exceeding threshold probability prob, i.e. find a random point in the right
  // tail of the gaussian(base,noise), provided that the tail integral = prob
  const double ksqrt2 = 1.41421356237309515e+00;
  return mean+sigma*ksqrt2*TMath::ErfcInverse(2*prob*(1.-gRandom->Rndm()));
}



#endif
