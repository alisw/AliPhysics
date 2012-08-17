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
#include <TMath.h>

class AliITSUSimuParam : public TObject {

 public:
  enum {kOldCouplingPixUpg,kNewCouplingPixUpg,kMaxCouplingOptPixUpg};
  //
  AliITSUSimuParam();
  AliITSUSimuParam(UInt_t nPixUpg);
  AliITSUSimuParam(const AliITSUSimuParam& simpar);
  // assignment operator 
  AliITSUSimuParam& operator=(const AliITSUSimuParam& source);
  ~AliITSUSimuParam();

  Double_t ApplyPixUpgBaselineAndNoise(UInt_t mod) const;
  Double_t CalcProbNoiseOverThreshold(UInt_t mod) const;
  //
  void     SetPixUpgThreshold(Double_t thresh, Double_t sigma, int mod=-1);
  void     GetPixUpgThreshold(UInt_t mod, Double_t& thresh, Double_t& sigma) const;
  Double_t GetPixUpgThreshold(UInt_t mod)                                    const;
  //
  void     SetPixUpgNoise(Double_t noise, Double_t baseline, Int_t mod=-1);
  void     GetPixUpgNoise(UInt_t mod,Double_t &noise, Double_t &baseline)    const;
  //
  void     SetPixUpgBiasVoltage(Double_t bias=18.182,Int_t mod=-1);
  Double_t GetPixUpgBiasVoltage(UInt_t mod)                                  const;

 
  void     SetGeVToCharge(Double_t gc=fgkNcompsDefault)                             {fGeVcharge = gc;}
  Double_t GetGeVToCharge()                                                  const  {return fGeVcharge;}
  Double_t GeVToCharge(Double_t gev)                                         const  {return gev/fGeVcharge;}
  //
  void     SetDistanceOverVoltage(Double_t d,Double_t v)                            {fDOverV = d/v;}
  void     SetDistanceOverVoltage(Double_t dv=fgkDOverVDefault)                     {fDOverV = dv;}
  Double_t GetDistanceOverVoltage()                                          const  {return fDOverV;}
  //
  void     SetPixUpgCouplingOption(UInt_t opt);
  UInt_t   GetPixUpgCouplingOption()                                         const  {return fPixUpgCouplOpt;}

  void     SetPixUpgCouplingParam(Double_t col, Double_t row)                       {fPixUpgCouplCol = col; fPixUpgCouplRow = row;}
  void     GetPixUpgCouplingParam(Double_t &col, Double_t &row)              const  {col = fPixUpgCouplCol; row = fPixUpgCouplRow;}

  void     SetPixUpgSigmaDiffusionAsymmetry(Double_t ecc)                           {fPixUpgEccDiff=ecc;}   
  void     GetPixUpgSigmaDiffusionAsymmetry(Double_t &ecc)                   const  {ecc=fPixUpgEccDiff;}

  void     SetPixUpgLorentzDrift(Bool_t ison)                                       {fPixUpgLorentzDrift=ison;}
  Bool_t   GetPixUpgLorentzDrift()                                           const  {return fPixUpgLorentzDrift;}
  void     SetPixUpgLorentzHoleWeight(Double_t weight)                               {fPixUpgLorentzHoleWeight=weight;}
  Double_t GetPixUpgLorentzHoleWeight()                                      const  {return fPixUpgLorentzHoleWeight;}
  
  void     SetPixUpgAddNoisyFlag(Bool_t value)                                      {fPixUpgAddNoisyFlag = value;}
  Bool_t   GetPixUpgAddNoisyFlag()                                           const  {return fPixUpgAddNoisyFlag;}
  void     SetPixUpgRemoveDeadFlag(Bool_t value)                                    {fPixUpgRemoveDeadFlag = value;}
  Bool_t   GetPixUpgRemoveDeadFlag()                                         const  {return fPixUpgRemoveDeadFlag;}
  //
  Double_t LorentzAngleElectron(Double_t bz)                                 const;
  Double_t LorentzAngleHole(Double_t bz)                                     const;
  //
  Double_t SigmaDiffusion3D(Double_t  l)                                     const;
  Double_t SigmaDiffusion2D(Double_t l)                                      const;
  Double_t SigmaDiffusion1D(Double_t l)                                      const;
  //
  virtual void Print(Option_t *opt = "")                                     const; 
  //
  static Double_t CalcProbNoiseOverThreshold(double base, double noise, double thresh);
  static Double_t GenerateNoiseQFunction(double prob, double mean, double sigma);
  //
 protected:

  static const Double_t fgkPixUpgBiasVoltageDefault;//default for fPixUpgBiasVoltage
  static const Double_t fgkPixUpgThreshDefault; //default for fThresh
  static const Double_t fgkPixUpgThrSigmaDefault; //default for fSigma
  static const Double_t fgkPixUpgCouplColDefault; //default for fPixUpgCouplCol
  static const Double_t fgkPixUpgCouplRowDefault; //default for fPixUpgCouplRow
  static const Double_t fgkPixUpgEccDiffDefault;//default for fPixUpgEccDiff
  static const Double_t fgkPixUpgLorentzHoleWeightDefault;//default for fPixUpgLorentzHoleWeight
  static const UInt_t   fgkPixUpgCouplingOptDefault;  // type of pixel Coupling (old or new)
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
  UInt_t     fNPixUpg;            // number of PixUpg type detectors
  UInt_t     fPixUpgCouplOpt;     // PixUpg Coupling Option
  Double_t   fPixUpgCouplCol;     // PixUpg Coupling parameter along the cols
  Double_t   fPixUpgCouplRow;     // PixUpg Coupling parameter along the rows
  Double_t   fPixUpgEccDiff;      // Eccentricity (i.e. asymmetry parameter) in the  Gaussian diffusion for PixUpg  
  Bool_t     fPixUpgLorentzDrift;     // Flag to decide whether to simulate the Lorentz Drift or not in PixUpg
  Double_t   fPixUpgLorentzHoleWeight;// Lorentz Angle is computed for PixUpg as average of Hole and Electron
  //                                    this parameter gives the relative weights between the two
  Bool_t     fPixUpgAddNoisyFlag;     // Flag saying whether noisy pixels should be added to digits
  Bool_t     fPixUpgRemoveDeadFlag;   // Flag saying whether dead pixels should be removed from digits
  //
  Double_t   fPixUpgThreshDef;      // PixUpg Threshold value
  Double_t   fPixUpgThrSigmaDef;    // PixUpg Threshold fluctuation
  Double_t   fPixUpgBiasVoltageDef; // Bias Voltage for the PixUpg
  Double_t   fPixUpgNoiseDef;       // PixUpg electronic noise: sigma
  Double_t   fPixUpgBaselineDef;    // PixUpg electronic noise: baseline
  //
  Double_t*  fPixUpgThresh;      //[fNPixUpg] PixUpg Threshold value
  Double_t*  fPixUpgThrSigma;    //[fNPixUpg] PixUpg Threshold fluctuation
  Double_t*  fPixUpgBiasVoltage; //[fNPixUpg] Bias Voltage for the PixUpg
  Double_t*  fPixUpgSigma;       //[fNPixUpg] PixUpg threshold fluctuations spread
  Double_t*  fPixUpgNoise;       //[fNPixUpg] PixUpg electronic noise: sigma
  Double_t*  fPixUpgBaseline;    //[fNPixUpg] PixUpg electronic noise: baseline
  //

  ClassDef(AliITSUSimuParam,1);
};

//_______________________________________________________________________
inline Double_t AliITSUSimuParam::CalcProbNoiseOverThreshold(double mean, double sigma, double thresh) 
{
  // calculate probability of gaussian noise exceeding the threshold
  const double ksqrt2 = 1.41421356237309515e+00;
  return 0.5*TMath::Erfc( (thresh-mean)/(sigma*ksqrt2));
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
