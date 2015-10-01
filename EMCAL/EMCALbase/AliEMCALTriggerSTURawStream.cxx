
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


#include "AliEMCALTriggerSTURawStream.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TBits.h"

#include <cstdlib>
#include <iostream>


using std::cout;
using std::setw;
using std::endl;
using std::dec;
ClassImp(AliEMCALTriggerSTURawStream)

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream() : TObject(),
  fRawReader(0x0)         ,
  fDetector(kEMCAL)       ,
  fGetRawData(0)          ,
  fPayload(V0)            ,
  fNL0GammaPatch(0)       ,
  fL0GammaPatchIndex()    ,
  fNL1GammaPatch()        ,
  fG()                    ,
  fL1GammaThreshold()     ,
  fL1GammaPatchIndex()    ,
  fNL1JetPatch()          ,
  fJ()                    ,
  fL1JetThreshold()       ,
  fL1JetPatchIndex()      ,
  fADC()                  ,
  fV0A(0)                 ,
  fV0C(0)                 ,
  fS()                    ,
  fRho(0)                 ,
  fPatchSize(0)           ,
  fRegionEnable(0)        ,
  fFrameReceived(0)       ,
  fFwVersion(0)     
{
  for(int i=0;i<max_L0GammaPatchIndex;i++) 
    fL0GammaPatchIndex[i]     = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
    fNL1GammaPatch   [i]      = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
    fL1GammaThreshold[i]      = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
  for(int j=0;j<3;                    j++)
    fG[j][i]                  = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
  for(int j=0;j<max_L1GammaPatchIndex;j++)
    fL1GammaPatchIndex[j][i]  = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
    fNL1JetPatch[i]           = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
    fL1JetThreshold[i]        = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
  for(int j=0;j<3;                    j++)
    fJ[j][i]                  = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
  for(int j=0;j<max_L1JetPatchIndex;  j++)
    fL1JetPatchIndex[j][i]    = 0 ;
  for(int i=0;i<max_nmoduleInTRU;     i++)
  for(int j=0;j<max_nTRU;             j++)
    fADC[j][i]                = 0 ;
  for(int i=0;i<4;                    i++)
    fS[i] = 0 ;
}

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream(AliRawReader* rawReader) : TObject(),
  fRawReader(rawReader)   ,
  fDetector(kEMCAL)       ,
  fGetRawData(0)          ,
  fPayload(V0)            ,
  fNL0GammaPatch(0)       ,
  fL0GammaPatchIndex()    ,
  fNL1GammaPatch()        ,
  fG()                    ,
  fL1GammaThreshold()     ,
  fL1GammaPatchIndex()    ,
  fNL1JetPatch()          ,
  fJ()                    ,
  fL1JetThreshold()       ,
  fL1JetPatchIndex()      ,
  fADC()                  ,
  fV0A(0)                 ,
  fV0C(0)                 ,
  fS()                    ,
  fRho(0)                 ,
  fPatchSize(0)           ,
  fRegionEnable(0)        ,
  fFrameReceived(0)       ,
  fFwVersion(0)     
{
  //
  fRawReader->Reset();
  fRawReader->Select("EMCAL",AliDAQ::GetFirstSTUDDL());
  //
  for(int i=0;i<max_L0GammaPatchIndex;i++) 
    fL0GammaPatchIndex[i]     = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
    fNL1GammaPatch   [i]      = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
    fL1GammaThreshold[i]      = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
  for(int j=0;j<3;                    j++)
    fG[j][i]                  = 0 ;
  for(int i=0;i<max_L1Gamma;          i++)
  for(int j=0;j<max_L1GammaPatchIndex;j++)
    fL1GammaPatchIndex[j][i]  = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
    fNL1JetPatch[i]           = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
    fL1JetThreshold[i]        = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
  for(int j=0;j<3;                    j++)
    fJ[j][i]                  = 0 ;
  for(int i=0;i<max_L1Jet;            i++)
  for(int j=0;j<max_L1JetPatchIndex;  j++)
    fL1JetPatchIndex[j][i]    = 0 ;
  for(int i=0;i<max_nmoduleInTRU;     i++)
  for(int j=0;j<max_nTRU;             j++)
    fADC[j][i]                = 0 ;
  for(int i=0;i<4;                    i++)
    fS[i] = 0 ;
}

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::~AliEMCALTriggerSTURawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::Reset()
{
  // Reset
  if (fRawReader) fRawReader->Reset();
  fNL0GammaPatch = 0;
  for(int i=0;i<max_L1Gamma;i++)fNL1GammaPatch[i] = 0 ;
  for(int i=0;i<max_L1Jet;  i++)fNL1JetPatch  [i] = 0 ;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::ReadPayLoad()
{
  UInt_t word32[max_payload_size]; // 32b words
  for (Int_t i=0; i<max_payload_size; i++) word32[i] = 0;
  
  Int_t iword = 0;
  
  fNL0GammaPatch = 0;
  for(int i=0;i<max_L1Gamma ;i++)    fNL1GammaPatch  [i] = 0 ;
  for(int i=0;i<max_L1Jet   ;i++)    fNL1JetPatch    [i] = 0 ;
  
  Int_t eqId    = -1  ;
  Int_t eqSize  =  0  ;
  
  UInt_t w32;
 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // migrate raw data to local array
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while (fRawReader->ReadNextInt(w32)){
    if (!iword){
      eqId   = fRawReader->GetEquipmentId();
      eqSize = fRawReader->GetEquipmentSize();
    }
    word32[iword++] = w32;
    //if(iword<905)cout<<dec<<iword+9    <<" , "<<hex<<w32<<endl;
    //else         cout<<dec<<iword-905+9<<" , "<<hex<<w32<<endl;

  }

  Int_t   poffset = 0     ;
  
  //payload type selector
  if(false){}
  else if(iword==(kPayLoadSizeV0                                 )){
    poffset         = 0         ;
    fPayload        = V0        ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV0       + kPayLoadSizeV0_Raw      )){
    poffset         = 0         ;
    fPayload        = V0Raw     ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1                                 )){
    poffset         = 0         ;
    fPayload        = V1        ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1       + kPayLoadSizeV1_Raw      )){
    poffset         = 0         ;
    fPayload        = V1Raw     ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1_2                               )){
    poffset         = 0         ;
    fPayload        = V1_2      ;
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1_2     + kPayLoadSizeV1_2_Raw    )){
    poffset         = 0         ;
    fPayload        = V1_2Raw   ;
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV2_DCAL                            + kPayLoadSizeV2_EMCAL                           )){
    poffset   = (fDetector == kDCAL)? 0           : kPayLoadSizeV2_DCAL                           ;
    fPayload  = (fDetector == kDCAL)? V2DCAL      : V2EMCAL                                       ;
  }
  else if(iword==(kPayLoadSizeV2_DCAL                            + kPayLoadSizeV2_EMCAL + kPayLoadSizeV2_EMCAL_Raw)){
    poffset   = (fDetector == kDCAL)? 0           : kPayLoadSizeV2_DCAL                           ;
    fPayload  = (fDetector == kDCAL)? V2DCAL      : V2EMCALRaw                                    ;
  }
  else if(iword==(kPayLoadSizeV2_DCAL  + kPayLoadSizeV2_DCAL_Raw + kPayLoadSizeV2_EMCAL                           )){
    poffset   = (fDetector == kDCAL)? 0           : kPayLoadSizeV2_DCAL + kPayLoadSizeV2_DCAL_Raw ;
    fPayload  = (fDetector == kDCAL)? V2DCALRaw   : V2EMCAL                                       ;
  }
  else if(iword==(kPayLoadSizeV2_DCAL  + kPayLoadSizeV2_DCAL_Raw + kPayLoadSizeV2_EMCAL + kPayLoadSizeV2_EMCAL_Raw)){
    poffset   = (fDetector == kDCAL)? 0           : kPayLoadSizeV2_DCAL + kPayLoadSizeV2_DCAL_Raw ;
    fPayload  = (fDetector == kDCAL)? V2DCALRaw   : V2EMCALRaw                                    ;
  }
  else{
    AliError(Form("STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32", eqId, eqSize, iword));
    return kFALSE;
  }
  AliDebug(1, Form("STU (eqId: %d, eqSize: %d) payload size: %d word32", eqId, eqSize, iword));

  switch (fPayload){
    case V0         :    case V1         :    case V2EMCAL    :    case V2DCAL     :  case V1_2    :  {fGetRawData=0;  break;}
    case V0Raw      :    case V1Raw      :    case V2EMCALRaw :    case V2DCALRaw  :  case V1_2Raw :  {fGetRawData=1;  break;}
    default :    {}
  }//end case

  int index         = poffset;
  int offset        = 0;
  int nJetThresh    = 0;
  int nGammaThresh  = 0;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // read setting parameters
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (fPayload){
    //##################################################
    case  V0         :
    case  V0Raw      :
    {
      fL1JetThreshold   [0] = ((word32[index]>>16) & 0xFFFF);
      fL1GammaThreshold [0] = ( word32[index]      & 0xFFFF);
      index++ ;
      offset        = 1 ;
      nJetThresh    = 1 ;
      nGammaThresh  = 1 ;
      break;
    }
    //##################################################
    case  V1         :
    case  V1Raw      :
    {
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][0] = word32[index++] ;
      fRegionEnable  = word32[index++];
      fFrameReceived = word32[index++];
      fFwVersion     = word32[index++];
      offset = index;
      
      fL1GammaThreshold [0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
      fL1JetThreshold   [0] = GetThreshold(fJ[0][0], fJ[1][0], fJ[2][0], fV0A, fV0C);
      nJetThresh    = 1;
      nGammaThresh  = 1;
      break;
    }
    //##################################################
    case  V1_2    :
    case  V1_2Raw :
    {
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][1] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][1] = word32[index++] ;
      fRegionEnable   = word32[index++] ;
      fFrameReceived  = word32[index++] ;
      fFwVersion      = ( word32[index]        & 0x0000FFFF)  ;
      fPatchSize      = ((word32[index] >> 16) & 0x0000FFFF)  ;
      index++ ;
      offset = index  ;      
      
      fL1GammaThreshold [0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
      fL1JetThreshold   [0] = GetThreshold(fJ[0][0], fJ[1][0], fJ[2][0], fV0A, fV0C);
      fL1GammaThreshold [1] = GetThreshold(fG[0][1], fG[1][1], fG[2][1], fV0A, fV0C); 
      fL1JetThreshold   [1] = GetThreshold(fJ[0][1], fJ[1][1], fJ[2][1], fV0A, fV0C);
      nJetThresh    = 2;
      nGammaThresh  = 2;
      break;
    }

    //##################################################
    case  V2EMCAL    :
    case  V2EMCALRaw :
    {
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][1] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][1] = word32[index++] ;
      fRho            = word32[index++]        & 0x3FFFFFFF   ;
      fRegionEnable   = word32[index++] ;
      fFrameReceived  = word32[index++] ;
      fFwVersion      = ( word32[index]        & 0x0000FFFF)  ;
      fPatchSize      = ((word32[index] >> 16) & 0x0000FFFF)  ;
      index++ ;
      offset = index  ;      
      
      fL1GammaThreshold [0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
      fL1JetThreshold   [0] = GetThreshold(/*fJ[0][0], fJ[1][0],*/0,0, fJ[2][0], fV0A, fV0C);
      fL1GammaThreshold [1] = GetThreshold(fG[0][1], fG[1][1], fG[2][1], fV0A, fV0C); 
      fL1JetThreshold   [1] = GetThreshold(/*fJ[0][1], fJ[1][1],*/0,0, fJ[2][1], fV0A, fV0C);
      nJetThresh    = 2;
      nGammaThresh  = 2;
      break;
    }
    //##################################################
    case  V2DCAL     :
    case  V2DCALRaw  :
    {
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][1] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][1] = word32[index++] ;
      for(int i=0;i<4;i++) fS[i]    = word32[index++] ;
      fRho            = word32[index++] ; 
      fRegionEnable   = word32[index++] ;
      fFrameReceived  = word32[index++] ;
      fFwVersion      = ( word32[index]        & 0x0000FFFF)  ;
      fPatchSize      = ((word32[index] >> 16) & 0x0000FFFF)  ;
      index++ ;
      offset = index  ;
      
      fL1GammaThreshold [0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
      fL1JetThreshold   [0] = GetThreshold(/*fJ[0][0], fJ[1][0],*/0,0, fJ[2][0], fV0A, fV0C);
      fL1GammaThreshold [1] = GetThreshold(fG[0][1], fG[1][1], fG[2][1], fV0A, fV0C); 
      fL1JetThreshold   [1] = GetThreshold(/*fJ[0][1], fJ[1][1],*/0,0, fJ[2][1], fV0A, fV0C);
      nGammaThresh  = 2;
      nJetThresh    = 2;
      break;
    }
    default :
    {
      AliError(Form("STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32", eqId, eqSize, iword));
      return kFALSE;
    }
  }//end switch
  
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // start decoding
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  int L1JetIndexOffset   = (fDetector == kEMCAL)? 11 : 
                           (fDetector == kDCAL )? 11 :
                           0 ;
  int L0GammaIndexOffset = (fDetector == kEMCAL)? 6 * nTRU_EMCAL / 2 : 
                           (fDetector == kDCAL )? 6 * nTRU_DCAL  / 2 :
                           0 ;
  int L1GammaIndexOffset = (fDetector == kEMCAL)? 8 * nTRU_EMCAL / 2 : 
                           (fDetector == kDCAL )? 8 * nTRU_DCAL  / 2 :
                           0 ;

  //Jet patch index
  for (int i = 0; i < nJetThresh; i++) {
    DecodeL1JetPatchIndexes(i, word32, offset);
    offset += L1JetIndexOffset ;
  }
  
  //L0 Gamma index
  DecodeL0GammaPatchIndexes(word32, offset);
  offset += L0GammaIndexOffset  ;
  
  //L1 Gamma index
  for (int i = 0; i < nGammaThresh; i++) {
    DecodeL1GammaPatchIndexes(i, word32, offset);  
    offset += L1GammaIndexOffset;
  }
 
  //decide if Raw data included 
  if (fGetRawData == 0){
    return kTRUE;
  }
  else{
    DecodeTRUADC(word32, offset);
    //DCAL : PHOSsubregion data
    if(fDetector == kDCAL){
      offset  += nTRU_DCAL * nMod_DCAL / 2  ;
      DecodePHOSSubregion(word32, offset);
    }
    return kTRUE;
  }
  return kFALSE ;
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL0GammaPatchIndexes(UInt_t *word32, const int offset)
{
  if(   fPayload == V0         
    ||  fPayload == V0Raw      
    ||  fPayload == V1         
    ||  fPayload == V1Raw      
    ||  fPayload == V1_2         
    ||  fPayload == V1_2Raw      
    )
  {
    unsigned short  truL0indexes[max_nTRU][6];
    Int_t nTRU  = GetnTRU();

    // extraction from stream
    for (Int_t index=0;index<6;index++)
    for (Int_t tru_num=0;tru_num<nTRU/2;tru_num++){
      truL0indexes[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
      truL0indexes[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
    }
    for (Int_t tru_num=0;tru_num<nTRU;tru_num++)
    for (Int_t index=0;index<6;index++)
    for (Int_t bit_num=0;bit_num<12;bit_num++){
      if ((truL0indexes[tru_num][index] & (1 << bit_num))){
        Int_t idx = 12 * index + bit_num;
        fNL0GammaPatch++;
        fL0GammaPatchIndex[fNL0GammaPatch-1] = (((idx << 5) & 0x7E0) | (tru_num & 0x1F));
      }
    }
  //V2 : L0 index is filled by 0  => no mean
  }else{
    return ; 
  }
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1JetPatchIndexes(const int i, UInt_t *word32, const int offset)
{
  Int_t nSubregion_eta  = (fDetector == kEMCAL)? nSubregion_eta_EMCAL :
                          (fDetector == kDCAL )? nSubregion_eta_DCAL  :
                          0 ;
  Int_t nSubregion_phi  = (fDetector == kEMCAL)? nSubregion_phi_EMCAL :
                          (fDetector == kDCAL )? nSubregion_phi_DCAL  :
                          0 ;

  int jetSize = 2 + fPatchSize ; //0->2x2, 2->4x4
  for (Int_t ix = 0; ix < nSubregion_eta - (jetSize-1) ; ix++){
    UInt_t currentrow = word32[offset + ix];
    for (Int_t iy = 0; iy < nSubregion_phi - (jetSize-1); iy++){
      if (currentrow & (1 << iy)){
        fNL1JetPatch[i]                          = fNL1JetPatch[i] + 1;
        fL1JetPatchIndex[fNL1JetPatch[i] - 1][i] = ((ix << 8) & 0xFF00) | (iy & 0xFF);
      }
    }
  }
  return ;
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset)
{
  Int_t nTRU  = GetnTRU();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(   fPayload == V0         
    ||  fPayload == V0Raw      
    ||  fPayload == V1         
    ||  fPayload == V1Raw      
    ||  fPayload == V1_2         
    ||  fPayload == V1_2Raw      
    )
  {
    unsigned short  truL1indexes0[max_nTRU][8];
    // extraction from stream
    for (Int_t index  =0; index   <8      ; index++   )
    for (Int_t tru_num=0; tru_num <nTRU/2 ; tru_num++ ){
      truL1indexes0[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
      truL1indexes0[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
    }  
    // interpretation
    int gammacolnum;
    short indexcopy;
    for (Int_t tru_num=0; tru_num<nTRU ; tru_num++  )
    for (Int_t index  =0; index  <8    ; index++    )
    for (Int_t bit_num=0; bit_num<12   ; bit_num++  ){
      if ((truL1indexes0[tru_num][index] & (1<<bit_num)) != 0){
        if (index<4){ // Even
          gammacolnum = (2*bit_num  );
          indexcopy   = index;
        }
        else{         // Odd
          gammacolnum = (2*bit_num+1);
          indexcopy   = index-4;
        }            
        fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
        fL1GammaPatchIndex[fNL1GammaPatch[i] - 1][i] = (
              ((indexcopy   << 10) & 0xC00) 
            | ((gammacolnum <<  5) & 0x3E0) 
            | ( tru_num            &  0x1F)
            );
      }//if
    }//for
  }//fPayload  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else if(fPayload == V2EMCAL
      ||  fPayload == V2EMCALRaw
      ||  fPayload == V2DCAL
      ||  fPayload == V2DCALRaw
      )
  {
    Int_t thirdTRU  = (fDetector == kEMCAL)? 30 : 
                      (fDetector == kDCAL )? 12 :
                      0 ;
    const Int_t nModPhi = 12  ;
    const Int_t nModEta =  8  ;
    unsigned short  truL1indexes1[max_nTRU][nModPhi];
    for (Int_t iphi=0; iphi<nModPhi/2 ; iphi++)
    for (Int_t itru=0; itru<nTRU/2    ; itru++){
      truL1indexes1[2*itru  ][2*iphi   ] = ( word32[offset + iphi * nTRU/2 + itru] >> 0  & 0xFF);
      truL1indexes1[2*itru  ][2*iphi+1 ] = ( word32[offset + iphi * nTRU/2 + itru] >> 8  & 0xFF);
      truL1indexes1[2*itru+1][2*iphi   ] = ( word32[offset + iphi * nTRU/2 + itru] >> 16 & 0xFF);
      truL1indexes1[2*itru+1][2*iphi+1 ] = ( word32[offset + iphi * nTRU/2 + itru] >> 24 & 0xFF);
    }  
    for (Int_t itru=0; itru<nTRU    ; itru++)
    for (Int_t iphi=0; iphi<nModPhi ; iphi++)
    for (Int_t ieta=0; ieta<nModEta ; ieta++){
      if ( (truL1indexes1[itru][iphi] >> ieta) & 0x1 ){
        fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
        
        short iphi_tmp  = (itru >= thirdTRU)? ( iphi%2 + 2* (int)(iphi/6)       ) : iphi  ;
        short ieta_tmp  = (itru >= thirdTRU)? ( ieta   + 8*((int)(iphi/2) % 3 ) ) : ieta  ;

        fL1GammaPatchIndex[fNL1GammaPatch[i] - 1][i] = (
              ((iphi_tmp  << 10) & 0x7C00) 
            | ((ieta_tmp  <<  5) & 0x03E0) 
            | ((itru      <<  0) & 0x001F)
            );
        
      }//if
    }//for
  }//fPayload  
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeTRUADC(UInt_t *word32, const int offset)
{
  Int_t nTRU  = GetnTRU();
  Int_t nMod  = GetnMod();

  // extraction from stream
  for (Int_t index=0;index<nMod;index++){
    for (Int_t tru_num=0;tru_num<nTRU/2;tru_num++){
      fADC[2*tru_num  ][index] = ( word32[offset + index * nTRU/2 + tru_num]        & 0xFFFF);
      fADC[2*tru_num+1][index] = ((word32[offset + index * nTRU/2 + tru_num] >> 16) & 0xFFFF);
    }
  }  
 
  if(   fPayload == V0Raw         
    ||  fPayload == V1Raw      
    ||  fPayload == V1_2Raw      
    )
  {
    for (Int_t tru_num=16;tru_num<32;tru_num++){ // A side
      Int_t v[96];
      for (Int_t index=0;index<96;index++) v[index] = fADC[tru_num][95-index];
      for (Int_t index=0;index<96;index++) fADC[tru_num][index] = v[index];
    }
  }
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodePHOSSubregion(UInt_t *word32, const int offset){
  for (Int_t index=0;index<36;index++){
    fPHOSSubregion[index] = word32[offset + index] & 0x7FFFF  ;
  }
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL0GammaPatch(const Int_t i, Int_t& tru, Int_t& idx) const
{
  // L0 gamma patches sent to STU (original access to L0 patch indexes)
  
  if (i > fNL0GammaPatch) return kFALSE;
  
  tru =  fL0GammaPatchIndex[i] & 0x1F;
  idx = (fL0GammaPatchIndex[i] & 0x7E0) >> 5;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& tru, Int_t& col, Int_t& row) const
{
  // L1 gamma patch indexes
  if (j >= max_L1Gamma || i > fNL1GammaPatch[j]) return kFALSE;
  row = (fL1GammaPatchIndex[i][j] >> 10) & 0x1F ;
  col = (fL1GammaPatchIndex[i][j] >>  5) & 0x1F ;
  tru = (fL1GammaPatchIndex[i][j] >>  0) & 0x1F ;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1JetPatch(const Int_t i, const Int_t j, Int_t& col, Int_t& row) const
{
  // L1 jet patch indexes
  if (j >= max_L1Jet || i > fNL1JetPatch[j]) return kFALSE;
  col = (fL1JetPatchIndex[i][j] >> 0) & 0xFF  ;
  row = (fL1JetPatchIndex[i][j] >> 8) & 0xFF  ;
  return kTRUE;
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::GetADC(Int_t iTRU, UInt_t ADC[])
{
  // Time sums
  Int_t nMod  = GetnMod() ;
  for (Int_t i=0; i<nMod; i++) ADC[i] = fADC[iTRU][i];
}
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::GetPHOSSubregion(UInt_t PHOSSubregion[])
{
  for (Int_t i=0; i<36; i++) PHOSSubregion[i] = fPHOSSubregion[i];
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DumpPayLoad(const Option_t *option) const
{
  // Dump STU payload
  
  TString op = option;
  
  printf("V0A:             %d\n"  , fV0A              );
  printf("V0C:             %d\n"  , fV0C              );

  printf("Number of L0:    %d\n"  , fNL0GammaPatch    );
  for (int i = 0; i < max_L1Gamma; i++) {
    for (int j = 0; j < 3         ; j++) {
      printf("G[%d][%d]: %d\n", j, i, fG[j][i]);
    }
    printf("Gamma threshold[%d]: %d\n", i, fL1GammaThreshold[i]);
    printf("Number of L1-g[%d]: %d\n", i, fNL1GammaPatch[i]);
  }
  for (int i = 0; i < max_L1Jet; i++) {
    for (int j = 0; j < 3       ; j++) {
      printf("J[%d][%d]: %d\n", j, i, fJ[j][i]);
    }
    printf("Jet Threshold[%d]:   %d\n", i, fL1JetThreshold[i]);
    printf("Number of L1-j[%d]: %d\n", i, fNL1JetPatch[i]);
  }

  printf("Jet Rho:         %d\n"  , fRho              );
  printf("RegionEnable:    %8x\n" , fRegionEnable     ); 
  printf("FrameReceived:   %8x\n" , fFrameReceived    );
  printf("FwVersion:       %x\n"  , fFwVersion        );
  printf("PatchSize:       %d\n"  , fPatchSize        );
  printf("RawData:         %d\n"  , fGetRawData       );

  Int_t itru, col, row;
  // 
  if (op.Contains("L0") || op.Contains("ALL")){
    for (Int_t i = 0;i < fNL0GammaPatch; i++) {
      if (GetL0GammaPatch(i,itru,col)) 
        cout  << dec
              << "> Found L0 gamma in TRU #"  << setw(2) << itru 
              << " at idx: "                  << setw(2) << col 
              << endl;
    }
  }
  // 
  if (op.Contains("L1") || op.Contains("ALL")) {
    for (int j = 0; j < max_L1Gamma; j++) {
      for (Int_t i = 0; i < fNL1GammaPatch[j]; i++) {
        if (GetL1GammaPatch(i, j, itru, col, row)) 
          cout  << dec
                << "> Found L1 gamma "  << j 
                << " in TRU #"          << setw(2) << itru 
                << " at: ( col: "       << setw(2) << col 
                << " , row: "           << setw(2) << row 
                << " )" << endl;
      }
    }
    for (int j = 0; j < max_L1Jet; j++) {
      for (Int_t i = 0; i < fNL1JetPatch[j]; i++) {
        if (GetL1JetPatch(i, j, col, row)) 
          cout  << dec 
                << "> Found L1 jet "  << j 
                << " at: ( col: "     << setw(2) << col 
                << " , row: "         << setw(2) << row 
                << " )" << endl;
      }
    }
  }
  
  if ((op.Contains("ADC") || op.Contains("ALL")) && fGetRawData) {
    for (Int_t i = 0; i < GetnTRU(); i++) {
      cout << "--------\n";
      cout << "TRU #" << setw(2) << dec << i << ":";
      for (Int_t j = 0;j < GetnMod(); j++) { 
        //TBits xadc(12); xadc.Set(12,&fADC[i][j]); 
        if ((j % 4) == 0) cout << endl;
        //cout << setw(2) << j << ": " << xadc << " ";
        printf("%2d: %3x / ",j,fADC[i][j]); 
      }
      cout << "\n";
    }
  }
}

//_____________________________________________________________________________
UShort_t AliEMCALTriggerSTURawStream::GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const
{
  // Get threshold 
  ULong64_t v0sum = V0A + V0C;
  ULong64_t sqrV0 = v0sum * v0sum;          
  sqrV0 *= A;
  sqrV0 >>= 32;
  v0sum *= B;
  v0sum >>= 16;
  return (UShort_t)(sqrV0 + v0sum + C);
}
