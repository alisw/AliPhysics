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
  static const Int_t  fNumberOfSuperModules = 20                                      ;//Total SM in EMCAL
  static const Int_t  fNTotalTRU            = 52                                      ;//Total TRU in EMCAL
  static const Int_t  fNModulesInTRUEta     =  8                                      ;//#FastOR/TRU in Eta
  static const Int_t  fNModulesInTRUPhi     = 12                                      ;//#FastOR/TRU in Phi
  static const Int_t  fNModulesInTRU        = fNModulesInTRUEta * fNModulesInTRUPhi   ;//#FastOR/TRU
  static const Int_t  fNTRUEta              =  3                                      ;//#TRUs/SM in Eta  
  static const Int_t  fNTRUPhi              =  1                                      ;//#TRUs/SM in Phi
  static const Int_t  fNTRU                 = fNTRUEta * fNTRUPhi                     ;//#TRUs/SM  
  static const Int_t  fNEta                 = fNModulesInTRUEta * fNTRUEta            ;//#FastOR/SM in Eta
  static const Int_t  fNPhi                 = fNModulesInTRUPhi * fNTRUPhi            ;//#FastOR/SM in Phi 
  static const Int_t  fSTURegionNEta        = 2/*Aside,Cside*/ * fNEta                ;//EMCAL+DCAL region eta size
  static const Int_t  fSTURegionNPhi        = (5 * fNPhi) + (1 * fNPhi/3)   /*EMCAL*/      
                                            + (3 * fNPhi) + (1 * fNPhi/3)   /*DCAL */ ;//#FastOR/EMCALs in Phi
  static const Int_t  fSTURegionN           = fSTURegionNEta * fSTURegionNPhi         ;//#FastOR/EMCALs
  
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
  enum fDetType { 
    kEMCAL  = 0, 
    kDCAL   = 1, 
    kPHOS   = 2
  }; 
	

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
  Bool_t  GetAbsFastORIndexFromPHOSSubregion( Int_t iPHOS, Int_t& id) const;

  //********************************************
  //Get TRU/SM/EMCAL Geometry from FastOR index
  //********************************************
  Bool_t  GetTRUFromAbsFastORIndex(             Int_t id, Int_t& iTRU , Int_t& iADC) const;
  Bool_t  GetPositionInTRUFromAbsFastORIndex(   Int_t id, Int_t& iTRU , Int_t& iEta, Int_t& iPhi) const;
  Bool_t  GetPositionInSMFromAbsFastORIndex(    Int_t id, Int_t& iSM  , Int_t& iEta, Int_t& iPhi) const;
  Bool_t  GetPositionInEMCALFromAbsFastORIndex( Int_t id,               Int_t& iEta, Int_t& iPhi) const;

  //********************************************
  //TRU vs. STU
  //********************************************
  Bool_t  GetTRUFromSTU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC, Int_t detector) const;
  Bool_t  GetSTUFromTRU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC                ) const;
  Bool_t  GetTRUFromSTU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi, Int_t detector) const;
  Bool_t  GetSTUFromTRU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi                ) const;

  //********************************************
  //Cell Index
  //********************************************
  Bool_t  GetFastORIndexFromCellIndex(Int_t id, Int_t& idx) const;
  Bool_t  GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4]) const;

  //********************************************
  //TRU index
  //********************************************
  Bool_t  GetTRUIndexFromSTUIndex(    Int_t id, Int_t& idx , Int_t detector ) const;
  Int_t   GetTRUIndexFromSTUIndex(    Int_t id             , Int_t detector ) const;
  Bool_t  GetSTUIndexFromTRUIndex(    Int_t id, Int_t& idx                  ) const;
  Int_t   GetSTUIndexFromTRUIndex(    Int_t id                              ) const;
  Bool_t  GetTRUIndexFromOnlineIndex( Int_t id, Int_t& idx) const{idx = id  ; return kTRUE  ;};
  Int_t   GetTRUIndexFromOnlineIndex( Int_t id            ) const{return id ;};
  Bool_t  GetOnlineIndexFromTRUIndex( Int_t id, Int_t& idx) const{idx = id  ; return kTRUE  ;};
  Int_t   GetOnlineIndexFromTRUIndex( Int_t id            ) const{return id ;};
  Int_t   GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm) const;

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
  Bool_t  fTRUIsCside[ fNTotalTRU] ;//
  
  Int_t   fSMFastOROffsetX[ fNumberOfSuperModules] ;//FastOR offset[#of SM ]
  Int_t   fSMFastOROffsetY[ fNumberOfSuperModules] ;//
  Int_t   fnFastORInSMPhi[  fNumberOfSuperModules] ;//SM size
  Int_t   fnFastORInSMEta[  fNumberOfSuperModules] ;//

  Int_t   fnModuleInEMCALPhi[5]    ;//#FastOR/EMCAL in Phi

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
  Bool_t  GetSMIsCside(Int_t iSM)   const {
    return (iSM%2 == 1)? kTRUE : kFALSE ;
  }
    
  ClassDef(AliEMCALTriggerMappingV2,1)
};
 
#endif

