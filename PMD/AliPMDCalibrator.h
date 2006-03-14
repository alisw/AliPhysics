#ifndef ALIPMDCALIBRATOR_H
#define ALIPMDCALIBRATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TTask;
class TObjArray;
class TH1F;  

class AliPMDCalibData;

class AliPMDCalibrator
{
 public:
  AliPMDCalibrator() ;              // ctor
  virtual ~AliPMDCalibrator() ;     // dtor
  virtual void Exec();
  void CalculateIsoCell();          //calculates gains
  void Init();
  Bool_t Store();
  
 private:
  Float_t fGainFact[2][24][96][96];
  TH1F *fHsmIso[2][24];             //histos of isolated cell modulewise
  TH1F *fHadcIso[2][24][96][96];    // histos of isolated cells cellwise
  AliPMDCalibData *fCalibData;

ClassDef(AliPMDCalibrator,1)        // description 
};
#endif // AliPMDCALIBRATOR_H
