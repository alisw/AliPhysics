#ifndef ALIPMDPEDESTAL_H
#define ALIPMDPEDESTAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;
class AliPMD;

class AliPMDPedestal: public TNamed
{
 public:
  AliPMDPedestal();
  AliPMDPedestal(const char* name);
  AliPMDPedestal(const AliPMDPedestal &pedestal);
  AliPMDPedestal& operator= (const AliPMDPedestal &pedestal);
  virtual ~AliPMDPedestal();
  void  Reset();
  void  SetPedMeanRms(Int_t det, Int_t smn, Int_t row, Int_t col,
		      Float_t pedmean, Float_t pedrms);
  Int_t GetPedMeanRms(Int_t det, Int_t smn, Int_t row, Int_t col) const;
  virtual void Print(Option_t *) const;
  
 protected:

  enum
      {
	  kDet    = 2,   // Number of planes
	  kModule = 24,  // Number of modules per plane
	  kRow    = 48,  // Row
          kCol    = 96   // Column
      };

  Int_t fPedMeanRms[kDet][kModule][kRow][kCol];

  ClassDef(AliPMDPedestal,1) // Pedestal class
};
#endif
