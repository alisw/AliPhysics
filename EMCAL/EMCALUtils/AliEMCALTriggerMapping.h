#ifndef ALIEMCALTRIGGERMAPPING_H
#define ALIEMCALTRIGGERMAPPING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
Mapping ABC 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/
#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliEMCALGeometry;

class AliEMCALTriggerMapping : public TObject 
{
public:
           AliEMCALTriggerMapping() : fNTRU(0), fGeometry(0x0) {}
           AliEMCALTriggerMapping(const Int_t ntru, const AliEMCALGeometry* geo) : fNTRU(ntru), fGeometry(geo) {}
  virtual ~AliEMCALTriggerMapping() {}
  
  virtual Bool_t GetAbsFastORIndexFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& id)                             const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInTRU(const Int_t iTRU, const Int_t iEta, const Int_t iPhi, Int_t& id) const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInSM(const Int_t  iSM, const Int_t iEta, const Int_t iPhi, Int_t& id)  const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInEMCAL(const Int_t iEta, const Int_t iPhi, Int_t& id)                 const = 0;
  virtual Bool_t GetTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iADC)                                  const = 0;
  virtual Bool_t GetPositionInTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi)           const = 0;
  virtual Bool_t GetPositionInSMFromAbsFastORIndex(const Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi)             const = 0;
  virtual Bool_t GetPositionInEMCALFromAbsFastORIndex(const Int_t id, Int_t& iEta, Int_t& iPhi)                      const = 0;
  virtual Bool_t GetFastORIndexFromCellIndex(const Int_t id, Int_t& idx)                                             const = 0;
  virtual Bool_t GetCellIndexFromFastORIndex(const Int_t id, Int_t idx[4])                                           const = 0;
  virtual Bool_t GetTRUIndexFromSTUIndex(const Int_t id, Int_t& idx)                                                 const = 0;
  virtual Bool_t GetTRUIndexFromOnlineIndex(const Int_t id, Int_t& idx)                                              const = 0;
  virtual Int_t  GetTRUIndexFromOnline(Int_t hwAdd, Int_t ddl, Int_t sm)                                             const = 0;
  virtual Bool_t GetOnlineIndexFromTRUIndex(const Int_t id, Int_t& idx)                                              const = 0;
  virtual Bool_t GetFastORIndexFromL0Index(const Int_t iTRU, const Int_t id, Int_t idx[], const Int_t size)          const = 0;

  virtual Int_t  GetTRUIndexFromSTUIndex(   const Int_t id)                                                          const = 0;
  virtual Int_t  GetTRUIndexFromOnlineIndex(const Int_t id)                                                          const = 0;
  virtual Int_t  GetOnlineIndexFromTRUIndex(const Int_t id)                                                          const = 0;
  
  virtual void  GetNTRU(Int_t& n) {n = fNTRU;}
  virtual Int_t GetNTRU() {return fNTRU;}
  
protected:  
  Int_t fNTRU;
  const AliEMCALGeometry* fGeometry;
  
private:
  AliEMCALTriggerMapping(const AliEMCALTriggerMapping& rhs);
  AliEMCALTriggerMapping& operator=(const AliEMCALTriggerMapping& rhs);
  
  ClassDef(AliEMCALTriggerMapping,1)
};
 
#endif

