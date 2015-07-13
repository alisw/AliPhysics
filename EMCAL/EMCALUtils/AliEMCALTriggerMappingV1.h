#ifndef ALIEMCALTRIGGERMAPPINGV1_H
#define ALIEMCALTRIGGERMAPPINGV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerMapping.h"
class AliEMCALGeometry;
class AliEMCALTriggerMappingV1 : public AliEMCALTriggerMapping 
{
public:
	
	               AliEMCALTriggerMappingV1();
		       AliEMCALTriggerMappingV1(Int_t ntru, const AliEMCALGeometry* geo);
              virtual ~AliEMCALTriggerMappingV1() {}

  Bool_t   GetAbsFastORIndexFromTRU(Int_t iTRU, Int_t iADC, Int_t& id) const;
  Bool_t                    GetAbsFastORIndexFromPositionInTRU(Int_t iTRU, Int_t iEta, Int_t iPhi, Int_t& id) const;	
  Bool_t                    GetAbsFastORIndexFromPositionInSM( Int_t  iSM, Int_t iEta, Int_t iPhi, Int_t& id) const;	
  Bool_t                    GetAbsFastORIndexFromPositionInEMCAL(                Int_t iEta, Int_t iPhi, Int_t& id) const;
  Bool_t                                GetTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iADC) const;
  Bool_t   GetPositionInTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const;
  Bool_t    GetPositionInSMFromAbsFastORIndex(Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const;
  Bool_t GetPositionInEMCALFromAbsFastORIndex(Int_t id, Int_t& iEta, Int_t& iPhi) const;
  Bool_t          GetFastORIndexFromCellIndex(Int_t id, Int_t& idx) const;
  Bool_t          GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4]) const;
  Bool_t              GetTRUIndexFromSTUIndex(Int_t id, Int_t& idx, Int_t detector) const;
  Int_t               GetTRUIndexFromSTUIndex(Int_t id, Int_t detector) const;
  Bool_t           GetTRUIndexFromOnlineIndex(Int_t id, Int_t& idx) const;
  Int_t            GetTRUIndexFromOnlineIndex(Int_t id) const;
  Int_t            GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm) const;
  Bool_t           GetOnlineIndexFromTRUIndex(Int_t id, Int_t& idx) const;
  Int_t            GetOnlineIndexFromTRUIndex(Int_t id) const;
  Bool_t            GetFastORIndexFromL0Index(Int_t iTRU, Int_t id, Int_t idx[], Int_t size) const;
	
private:
	                    AliEMCALTriggerMappingV1(const AliEMCALTriggerMappingV1& rhs);
	         AliEMCALTriggerMappingV1& operator=(const AliEMCALTriggerMappingV1& rhs);
	
	ClassDef(AliEMCALTriggerMappingV1,1)
};
 
#endif

