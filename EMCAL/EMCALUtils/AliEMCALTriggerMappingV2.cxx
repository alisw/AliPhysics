/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerMappingV2.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

ClassImp(AliEMCALTriggerMappingV2)

//________________________________________________________________________________________________
AliEMCALTriggerMappingV2::AliEMCALTriggerMappingV2() : AliEMCALTriggerMapping()
{
  // Ctor
  SetUniqueID(2);
  
  for(Int_t iTRU=0; iTRU<fNTotalTRU; iTRU++){
    fTRUFastOROffsetX[iTRU] = 0 ;
    fTRUFastOROffsetY[iTRU] = 0 ;
    fnFastORInTRUPhi[ iTRU] = 0 ;
    fnFastORInTRUEta[ iTRU] = 0 ;
    fTRUIsCside[      iTRU] = kFALSE  ;
  }
  for(Int_t iSM=0; iSM<fNumberOfSuperModules; iSM++){
    fSMFastOROffsetX[ iSM]  = 0 ;
    fSMFastOROffsetY[ iSM]  = 0 ;
    fnFastORInSMPhi[  iSM]  = 0 ;
    fnFastORInSMEta[  iSM]  = 0 ;
  }
  for(Int_t iB=0; iB<5; iB++){
    fnModuleInEMCALPhi[iB]  = 0 ;
  }
}
//________________________________________________________________________________________________
AliEMCALTriggerMappingV2::AliEMCALTriggerMappingV2(Int_t ntru, const AliEMCALGeometry* geo) : AliEMCALTriggerMapping(ntru, geo)
{
  // Ctor
  SetUniqueID(2);
  
  for(Int_t iTRU=0; iTRU<fNTotalTRU; iTRU++){
    fTRUFastOROffsetX[iTRU] = 0;
    fTRUFastOROffsetY[iTRU] = 0;
    fnFastORInTRUPhi[ iTRU] = 0;
    fnFastORInTRUEta[ iTRU] = 0;
    fTRUIsCside[      iTRU] = kFALSE  ;
  }
  for(Int_t iSM=0; iSM<fNumberOfSuperModules; iSM++){
    fSMFastOROffsetX[ iSM]  = 0;
    fSMFastOROffsetY[ iSM]  = 0;
    fnFastORInSMPhi[  iSM]  = 0;
    fnFastORInSMEta[  iSM]  = 0;
  }
  for(Int_t iB=0; iB<5; iB++){
    fnModuleInEMCALPhi[iB]  = 0 ;
  }

  Init_TRU_offset()  ;
  Init_SM_offset()   ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetAbsFastORIndexFromTRU(Int_t iTRU, Int_t iADC, Int_t& id) const
{
  //Trigger mapping method, get  FastOr Index from TRU
  if (iTRU > fNTotalTRU-1     || iTRU < 0 || 
      iADC > fNModulesInTRU-1 || iADC < 0
  ){
    AliError(Form("Out of range! iTRU=%d, iADC=%d", iTRU, iADC));	
    return kFALSE;
  }
  
  Int_t iADCtmp = (fTRUIsCside[ iTRU])? (fNModulesInTRU - iADC - 1) : iADC   ;
  Int_t x = fTRUFastOROffsetX[iTRU]                              + int(iADCtmp / fnFastORInTRUPhi[iTRU])  ;
  Int_t y = fTRUFastOROffsetY[iTRU] + fnFastORInTRUPhi[iTRU] - 1 - int(iADCtmp % fnFastORInTRUPhi[iTRU])  ;
  id      = y * fSTURegionNEta + x          ;
  id      = ConvAbsFastORIndexA2B(id)       ;
  return kTRUE  ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetAbsFastORIndexFromPositionInTRU(Int_t iTRU, Int_t iEta, Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, get Index if FastOr from Position in TRU
  if (iTRU > fNTotalTRU-1               || iTRU < 0 || 
      iEta > fnFastORInTRUEta[iTRU] - 1 || iEta < 0 ||
      iPhi > fnFastORInTRUPhi[iTRU] - 1 || iPhi < 0
  ){
    AliError(Form("Out of range! iTRU=%d, iEta=%d, iPhi=%d", iTRU, iEta, iPhi));	
    return kFALSE;
  }
  Int_t iEtatmp = iEta  ;//XXX
  Int_t iPhitmp = iPhi  ;//XXX
  //Int_t iEtatmp = ( fTRUIsCside[ iTRU])? (fnFastORInTRUEta[iTRU] - 1 - iEta) : iEta  ;
  //Int_t iPhitmp = (!fTRUIsCside[ iTRU])? (fnFastORInTRUPhi[iTRU] - 1 - iPhi) : iPhi  ;
  Int_t x = fTRUFastOROffsetX[iTRU] + iEtatmp ;
  Int_t y = fTRUFastOROffsetY[iTRU] + iPhitmp ;
  id      = y * fSTURegionNEta + x            ;
  id      = ConvAbsFastORIndexA2B(id)         ;
  return kTRUE  ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetAbsFastORIndexFromPositionInSM(Int_t  iSM, Int_t iEta, Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, from position in SM Index get FastOR index 
  if (iSM  > fNumberOfSuperModules-1  || iSM  < 0 || 
      iEta > fnFastORInSMEta[iSM] -1  || iEta < 0 ||
      iPhi > fnFastORInSMPhi[iSM] -1  || iPhi < 0
  ){
    AliError(Form("Out of range! iSM=%d, iEta=%d, iPhi=%d", iSM, iEta, iPhi));	
    return kFALSE;
  }
  //Int_t iEtatmp = (GetSMIsCside(iSM) && GetSMType(iSM) == kDCAL_Standard)?(iEta + 8):iEta ;
  //Int_t x = fSMFastOROffsetX[iSM] + iEtatmp ;
  Int_t x = fSMFastOROffsetX[iSM] + iEta    ;
  Int_t y = fSMFastOROffsetY[iSM] + iPhi    ;
  id      = y * fSTURegionNEta + x          ;
  id      = ConvAbsFastORIndexA2B(id)       ;
  return kTRUE  ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetAbsFastORIndexFromPositionInEMCAL(Int_t iEta, Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, from position in EMCAL Index get FastOR index 
  if (
      iEta > fSTURegionNEta - 1 || iEta < 0 || 
      iPhi > fSTURegionNPhi - 1 || iPhi < 0  
  ){
    AliError(Form("Out of range! eta: %2d phi: %2d", iEta, iPhi));
    return kFALSE;
  }
  id      = iPhi * fSTURegionNEta + iEta  ;
  id      = ConvAbsFastORIndexA2B(id)     ;
  return kTRUE  ;
}

//________________________________________________________________________________________________
Bool_t   AliEMCALTriggerMappingV2::GetAbsFastORIndexFromPHOSSubregion( Int_t iPHOS, Int_t& id) const
{
  if(iPHOS > 35 || iPHOS < 0){
    AliError(Form("Out of range! phos subregion index: %2d ", iPHOS));
    return kFALSE;
  }
  Int_t iEta  = 16  + 4 * (Int_t)(iPHOS % 4) ;
  Int_t iPhi  = 64  + 4 * (Int_t)(iPHOS / 4) ;
  return GetAbsFastORIndexFromPositionInEMCAL(iEta,iPhi,id);
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iADC) const
{
  //Trigger mapping method, get TRU number from FastOr Index
  Int_t iEta_TRU , iPhi_TRU , iSM , iEta_SM , iPhi_SM ;
  Int_t id_tmp = ConvAbsFastORIndexB2A(id);
  return GetInfoFromAbsFastORIndex(
    id_tmp   , 
    iTRU , iADC , iEta_TRU , iPhi_TRU , 
    iSM  ,        iEta_SM  , iPhi_SM  
    );
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetPositionInTRUFromAbsFastORIndex(Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in TRU from FasOr Index
  Int_t iADC , iSM , iEta_SM , iPhi_SM ;
  Int_t id_tmp = ConvAbsFastORIndexB2A(id);
  return GetInfoFromAbsFastORIndex(
    id_tmp   , 
    iTRU , iADC , iEta     , iPhi     , 
    iSM  ,        iEta_SM  , iPhi_SM  
    );
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetPositionInSMFromAbsFastORIndex(Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in Super Module from FasOr Index
  Int_t iTRU , iADC , iEta_TRU , iPhi_TRU ;
  Int_t id_tmp = ConvAbsFastORIndexB2A(id);
  return GetInfoFromAbsFastORIndex(
    id_tmp   , 
    iTRU , iADC , iEta_TRU , iPhi_TRU , 
    iSM  ,        iEta     , iPhi  
    );
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetPositionInEMCALFromAbsFastORIndex(Int_t id, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in EMCAL from FastOR index
  if (id > fSTURegionN-1 || id < 0){
    AliError("Id out of range!");
    return kFALSE;
  }
  Int_t id_tmp = ConvAbsFastORIndexB2A(id);
  iEta  = id_tmp % fSTURegionNEta ;
  iPhi  = id_tmp / fSTURegionNEta ;
  return kTRUE;
}


//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetFastORIndexFromCellIndex(Int_t id, Int_t& idx) const
{
  // Trigger mapping method, from cell index get FastOR index 

  Int_t iSupMod, nModule, nIphi, nIeta, iphim, ietam;
  Bool_t isOK = fGeometry->GetCellIndex( id, iSupMod, nModule, nIphi, nIeta );
  
  //- XXX fGeometry->GetModulePhiEtaIndexInSModule( iSupMod, nModule, iphim, ietam );
 
  //+ => XXX  
  fGeometry->GetCellPhiEtaIndexInSModule(iSupMod, nModule, nIphi, nIeta, iphim, ietam);
  //ietam:0-31  for DCAL Cside
  if( GetSMType(iSupMod)==kDCAL_Standard && (iSupMod%2)==1)
    fGeometry->ShiftOfflineToOnlineCellIndexes(iSupMod, iphim, ietam);
  //ietam:16-47 for DCAL Cside
  iphim /= 2  ;
  ietam /= 2  ;
  //+ <= XXX
  
  if (isOK && GetAbsFastORIndexFromPositionInSM(iSupMod, ietam, iphim, idx)) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetCellIndexFromFastORIndex(Int_t id, Int_t idx[4]) const
{
  //Trigger mapping method, from FASTOR index get cell index 
  Int_t iSM=-1, iEta=-1, iPhi=-1;
  if (GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi))
  {
    Int_t ix = 2 * iEta;
    Int_t iy = 2 * iPhi;

    //+ => XXX  
    //ietam:16-47 for DCAL Cside
    if( GetSMType(iSM)==kDCAL_Standard ){
      if( iSM%2==1 )
        fGeometry->ShiftOnlineToOfflineCellIndexes(iSM, iy, ix);
      if(ix < 0 || ix > 31)
        return kFALSE ;
    }
    //ietam:0-31  for DCAL Cside
    //+ <= XXX
    
    idx[0] = fGeometry->GetAbsCellIdFromCellIndexes(iSM, iy    , ix    );
    idx[1] = fGeometry->GetAbsCellIdFromCellIndexes(iSM, iy    , ix + 1);
    idx[2] = fGeometry->GetAbsCellIdFromCellIndexes(iSM, iy + 1, ix    );
    idx[3] = fGeometry->GetAbsCellIdFromCellIndexes(iSM, iy + 1, ix + 1);
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetTRUIndexFromSTUIndex(Int_t id, Int_t& idx, Int_t detector ) const
{
  //Trigger mapping method, from STU index get TRU index 
  idx = GetTRUIndexFromSTUIndex(id,detector);
  return (idx>0)?kTRUE:kFALSE;
}

//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV2::GetTRUIndexFromSTUIndex(Int_t id, Int_t detector) const
{
  //Trigger mapping method, from STU index get TRU index 
  if ((id > 31 && detector == kEMCAL) || (id > 13 && detector == kDCAL) || id < 0){
    AliError(Form("TRU index out of range: %d",id));
  }
  if(detector == kEMCAL ){
    return id ;
  }
  else if(detector == kDCAL){
    return 32 + ((int)(id/4) * 6 ) + ((id%4 < 2)?(id%4):(id%4+2))  ;
  }
  return -1 ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetSTUIndexFromTRUIndex(Int_t id, Int_t& idx) const
{
  //Trigger mapping method, from STU index get TRU index 
  idx = GetSTUIndexFromTRUIndex(id);
  return (idx>0)?kTRUE:kFALSE;
}
//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV2::GetSTUIndexFromTRUIndex(Int_t id) const
{
  if(id < 32)return id  ;
  else{
    Int_t STUid = id ;
    if(STUid >= 48 ) STUid -= 2 ;
    if(STUid >= 42 ) STUid -= 2 ;
    if(STUid >= 36 ) STUid -= 2 ;
    STUid -= 32 ;
    return STUid  ;
  }
  return -1 ;
}

///
/// \return TRU  global offline number from:
/// \param hwAdd: hardware address
/// \param ddl number
/// \param sm: uper-module number
/// Used in AliEMCALTriggerRawDigitMaker::Add()
///
Int_t  AliEMCALTriggerMappingV2::GetTRUIndexFromOnlineHwAdd(Int_t hwAdd, Int_t ddl, Int_t sm) const
{    
  // 1/3 SMs
  
  if ( sm == 10 ) return 30;
  if ( sm == 11 ) return 31;
  if ( sm == 18 ) return 50;
  if ( sm == 19 ) return 51;

  // Full EMCal/DCal SMs
  
  UShort_t iBranch = ( hwAdd >> 11 ) & 0x1;  // 0/1
  
  Int_t iTRU = ( (ddl << 1) | iBranch ) - 1; // 0..2
  
  iTRU = (sm%2) ? 2-iTRU : iTRU;
    
  if(sm < 10) iTRU +=  3 * sm;    // EMCal
  else        iTRU += (3 * sm-4); // DCal
  
  if (iTRU > fNTotalTRU - 1 || iTRU < 0) 
  {
    AliError(Form("TRU index out of range: %d",iTRU));
    return -1;
  }
  
  return iTRU;
}



//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetFastORIndexFromL0Index(Int_t iTRU, Int_t id, Int_t idx[], Int_t size) const
{
  //Trigger mapping method, from L0 index get FastOR index 

  if (size <= 0 ||size > 4){
    AliError("Size not supported!");
    return kFALSE;
  }
		
  Int_t motif[4];
  motif[0] = 0;
  motif[2] = 1;
  motif[1] = fnFastORInTRUPhi[iTRU]     ;
  motif[3] = fnFastORInTRUPhi[iTRU] + 1 ;
 
  switch (size)
  {
    case 1: // Cosmic trigger
      if (!GetAbsFastORIndexFromTRU(iTRU, id, idx[1])) return kFALSE;
      break;
    case 4: // 4 x 4
      for (Int_t k = 0; k < 4; k++)
      {
        Int_t iADC = fnFastORInTRUPhi[iTRU] * int(id/(fnFastORInTRUPhi[iTRU]-1))
                   + (id%(fnFastORInTRUPhi[iTRU]-1))
                   + motif[k]   ;
				
        if (!GetAbsFastORIndexFromTRU(iTRU, iADC, idx[k])) return kFALSE;
      }
      break;
    default:
      break;
  }
	
  return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::Init_TRU_offset(){
  fTRUFastOROffsetX[0]   = 0 ;
  fTRUFastOROffsetY[0]   = 0 ;
  Int_t iTRU  = 0 ;

  for(int iSM=0; iSM<fNumberOfSuperModules; iSM++){
    Int_t       SM_type   = GetSMType(iSM);
    Bool_t      IsCside   = GetSMIsCside( iSM);
    Int_t       TRU_type  = 0             ;
    
    //===================
    //TRU ieta/iphi size 
    Int_t nTRU_inSM         = fNTRU                ;
    Int_t nTRU_inSM_phi     = fNTRUPhi             ;
    Int_t nTRU_inSM_eta     = fNTRUEta             ;
    Int_t nModule_inTRU_phi = fNModulesInTRUPhi    ;
    Int_t nModule_inTRU_eta = fNModulesInTRUEta    ;
    if(     SM_type == kEMCAL_3rd
         || SM_type == kDCAL_Ext       ){
      nTRU_inSM         = (Int_t)((Float_t)nTRU_inSM         / 3. );
      nTRU_inSM_eta     = (Int_t)((Float_t)nTRU_inSM_eta     / 3. );
      nModule_inTRU_phi = (Int_t)((Float_t)nModule_inTRU_phi / 3. );
      nModule_inTRU_eta = nModule_inTRU_eta                  * 3   ;
    }
    
    //===================
    //TRU ieta/iphi offset calculation 
    for(Int_t i=0; i<nTRU_inSM; i++){
      fnFastORInTRUPhi[iTRU]  = nModule_inTRU_phi ;
      fnFastORInTRUEta[iTRU]  = nModule_inTRU_eta ;
      fTRUIsCside[     iTRU]  = IsCside           ;
    
      if((iTRU+1) >= fNTotalTRU)break;
      
      TRU_type  = 0 ;
      if( i==nTRU_inSM-1 && IsCside ){//last TRU in SM
        TRU_type = 1  ;//right 
      }
      if(      TRU_type == 0){
        fTRUFastOROffsetX[iTRU+1]  = fTRUFastOROffsetX[iTRU] + nModule_inTRU_eta      ;
        fTRUFastOROffsetY[iTRU+1]  = fTRUFastOROffsetY[iTRU]                          ;
      }else if(TRU_type == 1){
        fTRUFastOROffsetX[iTRU+1]  = 0                                                ;
        fTRUFastOROffsetY[iTRU+1]  = fTRUFastOROffsetY[iTRU] + nModule_inTRU_phi      ;
      }
      iTRU++  ;
    }//TRU loop
  }//SM loop
  
  return kTRUE ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::Init_SM_offset(){

  fSMFastOROffsetX[0]  = 0 ;
  fSMFastOROffsetY[0]  = 0 ;
  fnModuleInEMCALPhi[0] = 0 ;
  Int_t iB=0  ;

  Int_t SM_type_buf = -1  ;
  for(int iSM=0; iSM<fNumberOfSuperModules; iSM++){
    Int_t   SM_type = GetSMType(iSM);
    Bool_t  IsCside = GetSMIsCside( iSM);
    
    Int_t nModule_inSM_phi  = fNPhi    ;
    Int_t nModule_inSM_eta  = fNEta    ;
    
    if(     SM_type == kEMCAL_3rd
         || SM_type == kDCAL_Ext       ){
      nModule_inSM_phi  = (Int_t)((Float_t)nModule_inSM_phi      / 3.);
    }

    fnFastORInSMPhi[iSM]  = nModule_inSM_phi ;
    fnFastORInSMEta[iSM]  = nModule_inSM_eta ;
     
    if(!IsCside){ 
      if(SM_type_buf == SM_type){
        fnModuleInEMCALPhi[iB]    += nModule_inSM_phi  ;
      }else{
        fnModuleInEMCALPhi[iB+1]  =  fnModuleInEMCALPhi[iB] + nModule_inSM_phi  ;
        iB++  ;
      }
      SM_type_buf = SM_type ;
    }
  
    if( (iSM+1) >= fNumberOfSuperModules)break  ;

    if(IsCside){//right SM
      fSMFastOROffsetX[iSM+1]  = 0                                        ;
      fSMFastOROffsetY[iSM+1]  = fSMFastOROffsetY[iSM] + nModule_inSM_phi ;
    }
    else{//left SM
      fSMFastOROffsetX[iSM+1]  = fSMFastOROffsetX[iSM] + nModule_inSM_eta     ;
      fSMFastOROffsetY[iSM+1]  = fSMFastOROffsetY[iSM]                        ;
    }
  }//SM loop

  return kTRUE ;
}

//________________________________________________________________________________________________
Bool_t AliEMCALTriggerMappingV2::GetInfoFromAbsFastORIndex(//conv from A
    Int_t id, 
    Int_t& iTRU , Int_t& iADC , Int_t& iEta_TRU , Int_t& iPhi_TRU , 
    Int_t& iSM  ,               Int_t& iEta_SM  , Int_t& iPhi_SM  
    ) const
{

  if (id > fSTURegionN-1 || id < 0){
    AliError("Id out of range!");
    return kFALSE;
  }
  Int_t idB = ConvAbsFastORIndexA2B(id) ;
  
  iTRU      = idB / fNModulesInTRU  ;
  iADC      = idB % fNModulesInTRU  ;
  if( iTRU > fNTotalTRU - 1 ) return kFALSE;

  iEta_TRU  = iADC / fnFastORInTRUPhi[iTRU] ;
  iPhi_TRU  = iADC % fnFastORInTRUPhi[iTRU] ;
  iADC      = fnFastORInTRUPhi[iTRU] * (( fTRUIsCside[iTRU])? (fnFastORInTRUEta[iTRU] - 1 - iEta_TRU  ) : iEta_TRU ) 
            +                          ((!fTRUIsCside[iTRU])? (fnFastORInTRUPhi[iTRU] - 1 - iPhi_TRU  ) : iPhi_TRU )
            ;

  Int_t x = id % fSTURegionNEta ;
  Int_t y = id / fSTURegionNEta ;
  Int_t idtmp = (y<fnModuleInEMCALPhi[2])? id : (id + fNModulesInTRU * 4) ;
  iSM = 2 * (int)(idtmp/(2 * fNEta * fNPhi)) + (int)(fTRUIsCside[iTRU]);
  if( iSM > fNumberOfSuperModules - 1 ) return kFALSE ;

  iEta_SM = x   % fnFastORInSMEta[iSM] ;
  iPhi_SM = idB % fnFastORInSMPhi[iSM] ;
  
  return kTRUE;
}

//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV2::ConvAbsFastORIndexA2B(Int_t idA) const
{
  Int_t det_phi           = int(idA/fSTURegionNEta) ;
  Int_t nModule_inSM_phi  = fNPhi                   ;
  Int_t idB               = 0                       ;
  for(int i=1;i<5;i++){
    if(det_phi < fnModuleInEMCALPhi[i]){
      idB = fSTURegionNEta * fnModuleInEMCALPhi[i-1] ;
      if(i == 2 || i == 4) nModule_inSM_phi  /= 3;
      break ;
    }
  }

  Int_t tmp0 = idA - idB   ;
  Int_t tmp1 = (int)(tmp0 / (fSTURegionNEta * nModule_inSM_phi))      ;
  Int_t tmp2 = (int)(tmp0 % (fSTURegionNEta * nModule_inSM_phi))      ;
  idB += tmp1 * (fSTURegionNEta * nModule_inSM_phi)     ;
  idB += (int)(tmp2 % fSTURegionNEta) *nModule_inSM_phi ;
  idB += (int)(tmp2 / fSTURegionNEta)     ;

  return idB ;
}
//________________________________________________________________________________________________
Int_t AliEMCALTriggerMappingV2::ConvAbsFastORIndexB2A(Int_t idB) const
{
  Int_t det_phi   = int(idB/fSTURegionNEta) ;
  Int_t idA  = 0 ;
  Int_t nModule_inSM_phi  = fNPhi ;
  for(int i=1;i<5;i++){
    if(det_phi < fnModuleInEMCALPhi[i]){
      idA = fSTURegionNEta * fnModuleInEMCALPhi[i-1] ;
      if(i == 2 || i == 4) nModule_inSM_phi  /= 3;
      break ;
    }
  }

  Int_t tmp0 = idB - idA   ;
  Int_t tmp1 = (int)(tmp0 / (fSTURegionNEta * nModule_inSM_phi))  ;
  Int_t tmp2 = (int)(tmp0 % (fSTURegionNEta * nModule_inSM_phi))  ;
  Int_t x = tmp2 / nModule_inSM_phi;
  Int_t y = tmp2 % nModule_inSM_phi;

  idA += tmp1 * (fSTURegionNEta * nModule_inSM_phi);
  idA += y*fSTURegionNEta + x;

  return idA  ;
}

//________________________________________________________________________________________________
Bool_t  AliEMCALTriggerMappingV2::GetTRUFromSTU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC, Int_t detector)const 
{
  Int_t ieta, iphi, oeta, ophi;
  oTRU  = GetTRUIndexFromSTUIndex(iTRU,detector);
  if (oTRU == -1) return kFALSE;
  ieta  = iADC % fnFastORInTRUEta[oTRU] ;
  iphi  = iADC / fnFastORInTRUEta[oTRU] ;
  oeta  = (fTRUIsCside[oTRU])? (fnFastORInTRUEta[oTRU] - ieta - 1) : ieta;
  ophi  = (fTRUIsCside[oTRU])? iphi : (fnFastORInTRUPhi[oTRU] - iphi - 1);
  oADC  = oeta * fnFastORInTRUPhi[oTRU] + ophi  ;
  return kTRUE  ;
}
//________________________________________________________________________________________________
Bool_t  AliEMCALTriggerMappingV2::GetTRUFromSTU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi, Int_t detector) const
{
  oTRU  = GetTRUIndexFromSTUIndex(iTRU,detector);
  if (oTRU == -1) return kFALSE;
  oeta  = (fTRUIsCside[oTRU])? (fnFastORInTRUEta[oTRU] - ieta - 1) : ieta;
  ophi  = (fTRUIsCside[oTRU])? iphi : (fnFastORInTRUPhi[oTRU] - iphi - 1);
  return kTRUE  ;
}
//________________________________________________________________________________________________
Bool_t  AliEMCALTriggerMappingV2::GetSTUFromTRU(Int_t iTRU, Int_t iADC, Int_t& oTRU, Int_t& oADC)const
{
  Int_t ieta, iphi, oeta, ophi;
  ieta  = iADC / fnFastORInTRUPhi[iTRU] ;
  iphi  = iADC % fnFastORInTRUPhi[iTRU] ;
  oeta  = (fTRUIsCside[iTRU])? (fnFastORInTRUEta[iTRU] - ieta - 1) : ieta;
  ophi  = (fTRUIsCside[iTRU])? iphi : (fnFastORInTRUPhi[iTRU] - iphi - 1);
  oTRU  = GetSTUIndexFromTRUIndex(iTRU);
  oADC  = ophi * fnFastORInTRUEta[iTRU] + oeta  ;
  return kTRUE  ;
}
//________________________________________________________________________________________________
Bool_t  AliEMCALTriggerMappingV2::GetSTUFromTRU(Int_t iTRU, Int_t ieta, Int_t iphi, Int_t& oTRU, Int_t& oeta, Int_t& ophi                ) const
{
  oTRU  = GetSTUIndexFromTRUIndex(iTRU);
  oeta  = (fTRUIsCside[iTRU])? (fnFastORInTRUEta[iTRU] - ieta - 1) : ieta;
  ophi  = (fTRUIsCside[iTRU])? iphi : (fnFastORInTRUPhi[iTRU] - iphi - 1);
  return kTRUE  ;
}


