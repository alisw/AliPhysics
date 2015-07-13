#ifndef ALIEMCALTRIGGERMAPPINGV2_H
#define ALIEMCALTRIGGERMAPPINGV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
 
 
Author: 
H. YOKOYAMA Tsukuba University
R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerMapping.h"
class AliEMCALGeometry;
class AliEMCALTriggerMappingV2 : public AliEMCALTriggerMapping 
{
public:

  //********************************************
  //static constant
  //********************************************
  static const Int_t  fNumberOfSuperModules = 20  ;//Total SM in EMCAL
  static const Int_t  fNTotalTRU            = 52  ;//Total TRU in EMCAL
  static const Int_t  fNTRU                 =  3  ;//#TRUs/SM  
  static const Int_t  fNTRUEta              =  3  ;//#TRUs/SM in Eta  
  static const Int_t  fNTRUPhi              =  1  ;//#TRUs/SM in Phi
  static const Int_t  fNEta                 = 24  ;//#FastOR/SM in Eta
  static const Int_t  fNPhi                 = 12  ;//#FastOR/SM in Phi 
  static const Int_t  fNModulesInTRU        = 96  ;//#FastOR/TRU
  static const Int_t  fNModulesInTRUEta     =  8  ;//#FastOR/TRU in Eta
  static const Int_t  fNModulesInTRUPhi     = 12  ;//#FastOR/TRU in Phi
  static const Int_t  fSTURegionNEta        = 2 * fNTRUEta * fNModulesInTRUEta;                //EMCAL+DCAL region eta size
  static const Int_t  fSTURegionNPhi        = 8 * fNTRUPhi * fNModulesInTRUPhi + 2 * fNPhi / 3;//EMCAL+DCAL region phi size
  
  //********************************************
  //SM type 
  //********************************************
  enum fEMCSMType { 
    kEMCAL_Standard = 0, 
    kEMCAL_Half     = 1, 
    kEMCAL_3rd      = 2, 
    kDCAL_Standard  = 3, 
    kDCAL_Ext       = 4 
  }; // possible SM Type
	

  //********************************************
  //constructor,destructor
  //********************************************
  AliEMCALTriggerMappingV2();
  AliEMCALTriggerMappingV2(Int_t ntru, const AliEMCALGeometry* geo);
  virtual ~AliEMCALTriggerMappingV2() {}

  //********************************************
  //Get FastOR index from TRU/SM/EMCAL Geometry
  //********************************************
  Bool_t  GetAbsFastORIndexFromTRU(          Int_t iTRU, Int_t iADC, Int_t& id) const;
  Bool_t  GetAbsFastORIndexFromPositionInTRU(Int_t iTRU, Int_t iEta, Int_t iPhi, Int_t& id) const;	
  Bool_t  GetAbsFastORIndexFromPositionInSM( Int_t  iSM, Int_t iEta, Int_t iPhi, Int_t& id) const;	
  Bool_t  GetAbsFastORIndexFromPositionInEMCAL(                Int_t iEta, Int_t iPhi, Int_t& id) const;

  //********************************************
  //Get TRU/SM/EMCAL Geometry from FastOR index
  //********************************************
  Bool_t  GetTRUFromAbsFastORIndex(             Int_t id, Int_t& iTRU , Int_t& iADC) const;
  Bool_t  GetPositionInTRUFromAbsFastORIndex(   Int_t id, Int_t& iTRU , Int_t& iEta, Int_t& iPhi) const;
  Bool_t  GetPositionInSMFromAbsFastORIndex(    Int_t id, Int_t& iSM  , Int_t& iEta, Int_t& iPhi) const;
  Bool_t  GetPositionInEMCALFromAbsFastORIndex( Int_t id,               Int_t& iEta, Int_t& iPhi) const;

  //********************************************
  //Cell Index
  //********************************************
  Bool_t  GetFastORIndexFromCellIndex(Int_t id, Int_t& idx) const;
  Bool_t  GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4]) const;

  //********************************************
  //TRU index
  //********************************************
  Bool_t  GetTRUIndexFromSTUIndex(    Int_t id, Int_t& idx) const;
  Int_t   GetTRUIndexFromSTUIndex(    Int_t id            ) const;
  Bool_t  GetTRUIndexFromOnlineIndex( Int_t id, Int_t& idx) const;
  Int_t   GetTRUIndexFromOnlineIndex( Int_t id            ) const;
  Int_t   GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm) const;
  Bool_t  GetOnlineIndexFromTRUIndex( Int_t id, Int_t& idx) const;
  Int_t   GetOnlineIndexFromTRUIndex( Int_t id            ) const;

  //********************************************
  //L0 Index
  //********************************************
  Bool_t  GetFastORIndexFromL0Index(Int_t iTRU, Int_t id, Int_t idx[], Int_t size) const;

  Bool_t GetInfoFromAbsFastORIndex(
    Int_t id, 
    Int_t& iTRU , Int_t& iADC , Int_t& iEta_TRU , Int_t& iPhi_TRU , 
    Int_t& iSM  ,               Int_t& iEta_SM  , Int_t& iPhi_SM  
    ) const;
	
private:
	         
  AliEMCALTriggerMappingV2(const AliEMCALTriggerMappingV2& rhs);
  AliEMCALTriggerMappingV2& operator=(const AliEMCALTriggerMappingV2& rhs);
	
  //********************************************
  //fastOR offset parameters
  //********************************************
  Int_t   fTRUFastOROffsetX[fNTotalTRU] ;//FastOR offset[#of TRU]
  Int_t   fTRUFastOROffsetY[fNTotalTRU] ;//
  Int_t   fnFastORInTRUPhi[ fNTotalTRU] ;//TRU size
  Int_t   fnFastORInTRUEta[ fNTotalTRU] ;//
  
  Int_t   fSMFastOROffsetX[ 20] ;//FastOR offset[#of SM ]
  Int_t   fSMFastOROffsetY[ 20] ;//
  Int_t   fnFastORInSMPhi[  20] ;//SM size
  Int_t   fnFastORInSMEta[  20] ;//

  Int_t   fnModuleInEMCALPhi    ;//#FastOR/EMCAL in Phi

  //********************************************
  //Initialization of FastOR index offset of each SM/TRU
  //********************************************
  Bool_t Init_TRU_offset()  ;
  Bool_t Init_SM_offset()   ;
  
  //********************************************
  //convert AbsFastORIndex from type-A(B) to type-B(A)
  //********************************************
  Int_t ConvAbsFastORIndexA2B(  Int_t idA) const  ;
  Int_t ConvAbsFastORIndexB2A(  Int_t idB) const  ;

  //********************************************
  //SM type
  //********************************************
  Int_t   GetSMType(Int_t iSM)      const {
    if( iSM<0 || iSM >= fNumberOfSuperModules)return -1  ;
    if( iSM < 10) return kEMCAL_Standard ;
    if( iSM < 12) return kEMCAL_3rd      ;
    if( iSM < 18) return kDCAL_Standard  ;
    if( iSM < 20) return kDCAL_Ext       ;
    return -1 ;
  }
    
  ClassDef(AliEMCALTriggerMappingV2,1)
};
 
#endif

