#ifndef ALIPMDMAPPINGDATA_H
#define ALIPMDMAPPINGDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;
class AliPMD;

class AliPMDMappingData: public TNamed
{
 public:
  AliPMDMappingData();
  AliPMDMappingData(const char* name);
  AliPMDMappingData(const AliPMDMappingData &mapda);
  AliPMDMappingData& operator= (const AliPMDMappingData &mapda);
  virtual ~AliPMDMappingData();
  void  Reset();

  Int_t GetBeginPatchBus(Int_t iddl, Int_t imod) const;
  Int_t GetEndPatchBus(Int_t iddl, Int_t imod) const;
  Int_t GetModuleNo(Int_t iddl, Int_t ibus) const;
  Int_t GetMcmperBus(Int_t iddl, Int_t ibus) const;
  Int_t GetStartRowBus(Int_t iddl, Int_t ibus) const;
  Int_t GetEndRowBus(Int_t iddl, Int_t ibus) const;
  Int_t GetStartColBus(Int_t iddl, Int_t ibus) const;
  Int_t GetEndColBus(Int_t iddl, Int_t ibus) const;

  
  void  SetPatchBus(Int_t iddl, Int_t imod, Int_t bpatchbus, Int_t epatchbus);
  void  SetModuleNo(Int_t iddl, Int_t ibus, Int_t modno);
  void  SetMcmperBus(Int_t iddl, Int_t ibus, Int_t totmcm);
  void  SetRowBus(Int_t iddl, Int_t ibus, Int_t rows, Int_t rowe);
  void  SetColBus(Int_t iddl, Int_t ibus, Int_t cols, Int_t cole);

  virtual void Print(Option_t *) const;
  
 protected:
  enum
      {
	kDdl = 6,     // Number of DDL
	kBus = 51    // Modules of patch bus
      };

  Int_t fBeginPatchBus[6][48];
  Int_t fEndPatchBus[6][48];
  Int_t fModuleNo[kDdl][kBus];
  Int_t fMcmperBus[kDdl][kBus];
  Int_t fStartRowBus[kDdl][kBus];
  Int_t fEndRowBus[kDdl][kBus];
  Int_t fStartColBus[kDdl][kBus];
  Int_t fEndColBus[kDdl][kBus];


  ClassDef(AliPMDMappingData,2) // calibration class for gainfactors
};
#endif
