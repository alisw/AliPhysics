#ifndef ALIPMDCALIBDATA_H
#define ALIPMDCALIBDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;
class AliPMD;

class AliPMDCalibData: public TNamed
{
 public:
  AliPMDCalibData();
  AliPMDCalibData(const char* name);
  AliPMDCalibData(const AliPMDCalibData &calibda);
  AliPMDCalibData& operator= (const AliPMDCalibData &calibda);
  virtual ~AliPMDCalibData();
  void    Reset();
  void    SetGainFact(Int_t det, Int_t smn, Int_t row, Int_t col,
		      Float_t gain);
  Float_t GetGainFact(Int_t det, Int_t smn, Int_t row, Int_t col) const;
  virtual void Print(Option_t *) const;
  
 protected:

  enum
      {
	  kDet = 2,        // Number of plane
	  kModule = 24,    // Modules per plane
	  kRow    = 48,    // Maximum row
	  kCol    = 96     // Maximum Column
      };
  Float_t fGainFact[kDet][kModule][kRow][kCol];

  ClassDef(AliPMDCalibData,2) // calibration class for gainfactors
};
#endif
