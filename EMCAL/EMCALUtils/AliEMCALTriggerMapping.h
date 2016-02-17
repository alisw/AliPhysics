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
           AliEMCALTriggerMapping(Int_t ntru, const AliEMCALGeometry* geo) : fNTRU(ntru), fGeometry(geo) {}
  virtual ~AliEMCALTriggerMapping() {}
  
  virtual Bool_t GetAbsFastORIndexFromTRU(Int_t iTRU, Int_t iADC, Int_t& id)                             const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInTRU(Int_t iTRU, Int_t iEta, Int_t iPhi, Int_t& id)       const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInSM(Int_t  iSM, Int_t iEta, Int_t iPhi, Int_t& id)        const = 0;
  virtual Bool_t GetAbsFastORIndexFromPositionInEMCAL(Int_t iEta, Int_t iPhi, Int_t& id)                 const = 0;
  virtual Bool_t GetAbsFastORIndexFromPHOSSubregion( Int_t iPHOS, Int_t& id)                             const = 0;
  virtual Bool_t GetTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iADC)                            const = 0;
  virtual Bool_t GetPositionInTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi)     const = 0;
  virtual Bool_t GetPositionInSMFromAbsFastORIndex(Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi)       const = 0;
  virtual Bool_t GetPositionInEMCALFromAbsFastORIndex(Int_t id, Int_t& iEta, Int_t& iPhi)                const = 0;
  virtual Bool_t GetFastORIndexFromCellIndex(Int_t id, Int_t& idx)                                       const = 0;
  virtual Bool_t GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4])                                     const = 0;
  virtual Bool_t GetTRUIndexFromSTUIndex(Int_t id, Int_t& idx, Int_t detector)                           const = 0;
  virtual Bool_t GetTRUIndexFromOnlineIndex(Int_t id, Int_t& idx)                                        const = 0;
  virtual Int_t  GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm)                            const = 0;
  virtual Bool_t GetOnlineIndexFromTRUIndex(Int_t id, Int_t& idx)                                        const = 0;
  virtual Bool_t GetFastORIndexFromL0Index(Int_t iTRU, Int_t id, Int_t idx[], Int_t size)                const = 0;
  virtual Int_t  GetTRUIndexFromSTUIndex(   Int_t id, Int_t detector)                                    const = 0;
  virtual Int_t  GetTRUIndexFromOnlineIndex(Int_t id)                                                    const = 0;
  virtual Int_t  GetOnlineIndexFromTRUIndex(Int_t id)                                                    const = 0;


  virtual Bool_t  GetSTUIndexFromTRUIndex(    Int_t id, Int_t& idx                              ) const = 0 ;
  virtual Int_t   GetSTUIndexFromTRUIndex(    Int_t id                                          ) const = 0 ;
  virtual Bool_t  GetTRUFromSTU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC, Int_t detector) const = 0 ;
  virtual Bool_t  GetSTUFromTRU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC                ) const = 0 ;
  virtual Bool_t  GetTRUFromSTU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi, Int_t detector) const = 0 ;
  virtual Bool_t  GetSTUFromTRU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi                ) const = 0 ;


  virtual void  GetNTRU(Int_t& n) { n = fNTRU    ; }
  virtual Int_t GetNTRU()         { return fNTRU ; }
  
protected:  
  Int_t fNTRU;
  const AliEMCALGeometry* fGeometry;
  
private:
  AliEMCALTriggerMapping(           const AliEMCALTriggerMapping& rhs);
  AliEMCALTriggerMapping& operator=(const AliEMCALTriggerMapping& rhs);
  
  ClassDef(AliEMCALTriggerMapping,1)
};
 
#endif

