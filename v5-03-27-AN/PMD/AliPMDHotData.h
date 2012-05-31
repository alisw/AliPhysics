#ifndef ALIPMDHOTDATA_H
#define ALIPMDHOTDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;

class AliPMDHotData: public TNamed
{
 public:
  AliPMDHotData();
  AliPMDHotData(const char* name);
  AliPMDHotData(const AliPMDHotData &hotda);
  AliPMDHotData& operator= (const AliPMDHotData &hotda);
  virtual ~AliPMDHotData();
  void    Reset();
  void    SetHotChannel(Int_t det, Int_t smn, Int_t row, Int_t col,Float_t flag);
  Float_t GetHotChannel(Int_t det, Int_t smn, Int_t row, Int_t col) const;
  virtual void Print(Option_t *) const;
  
 protected:

  enum
      {
	  kDet = 2,        // Number of plane
	  kModule = 24,    // Modules per plane
	  kRow    = 48,    // Maximum row
	  kCol    = 96     // Maximum Column
      };
  Float_t fHotChannel[kDet][kModule][kRow][kCol];

  ClassDef(AliPMDHotData,1) // class for hot cells in PMD
};
#endif
