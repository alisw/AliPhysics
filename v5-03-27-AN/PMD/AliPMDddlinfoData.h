#ifndef ALIPMDDDLINFODATA_H
#define ALIPMDDDLINFODATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TNamed;
class AliCDBEntry;

class AliPMDddlinfoData: public TNamed
{
 public:
  AliPMDddlinfoData();
  AliPMDddlinfoData(const char* name);
  AliPMDddlinfoData(const AliPMDddlinfoData &ddlinfoda);
  AliPMDddlinfoData& operator= (const AliPMDddlinfoData &ddlinfoda);
  virtual ~AliPMDddlinfoData();
  void  Reset();

  Int_t GetNoOfModulePerDdl(Int_t iddl) const;
  Int_t GetModulesPerDdl(Int_t iddl, Int_t imod) const;
  Int_t GetStartRowA(Int_t idet, Int_t ismn) const;
  Int_t GetStartRowB(Int_t idet, Int_t ismn) const;
  Int_t GetEndRowA(Int_t idet, Int_t ismn) const;
  Int_t GetEndRowB(Int_t idet, Int_t ismn) const;
  Int_t GetStartColA(Int_t idet, Int_t ismn) const;
  Int_t GetStartColB(Int_t idet, Int_t ismn) const;
  Int_t GetEndColA(Int_t idet, Int_t ismn) const;
  Int_t GetEndColB(Int_t idet, Int_t ismn) const;

  void SetNoOfModulePerDdl(Int_t iddl, Int_t nmod);
  void SetModuleNoPerDdl(Int_t iddl, Int_t mod[]);
  void SetStartRowA(Int_t srowa[][24]);
  void SetStartRowB(Int_t srowb[][24]);
  void SetEndRowA(Int_t erowa[][24]);
  void SetEndRowB(Int_t erowb[][24]);
  void SetStartColA(Int_t scola[][24]);
  void SetStartColB(Int_t scolb[][24]);
  void SetEndColA(Int_t ecola[][24]);
  void SetEndColB(Int_t ecolb[][24]);

  virtual void Print(Option_t *) const;
  
 protected:

  Int_t fModules[6];        // Total no. of modules per DDL
  Int_t fModuleNo[6][12];   // Serial Module nos. per DDL (12 nos)

  Int_t fStartRowA[2][24];  // removed from row A
  Int_t fStartRowB[2][24];  // removed from row B
  Int_t fEndRowA[2][24];    // removed upto row A
  Int_t fEndRowB[2][24];    // removed upto row B
  Int_t fStartColA[2][24];  // removed from col A
  Int_t fStartColB[2][24];  // removed from col B
  Int_t fEndColA[2][24];    // removed upto row A
  Int_t fEndColB[2][24];    // removed upto row B

  ClassDef(AliPMDddlinfoData,1) // ddlinfo database
};
#endif
