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

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerSTURawStream) ;
/// \endcond

///
/// Default constructor
//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream() : TObject(),
  fRawReader(0x0)         ,
  fDetector(kEMCAL)       ,
  fGetRawData(0)          ,
  fPayload(v0)            ,
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
  for(int i=0;i<maxL0GammaPatchIndex;i++) 
    fL0GammaPatchIndex[i]     = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    fNL1GammaPatch   [i]      = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    fL1GammaThreshold[i]      = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    for(int j=0;j<3;                 j++)
      fG[j][i]                  = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    for(int j=0;j<maxL1GammaPatchIndex;j++)
      fL1GammaPatchIndex[j][i]  = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    fNL1JetPatch[i]           = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    fL1JetThreshold[i]        = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    for(int j=0;j<3;                 j++)
      fJ[j][i]                  = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    for(int j=0;j<maxL1JetPatchIndex;j++)
      fL1JetPatchIndex[j][i]    = 0 ;
  for(int i=0;i<maxnmoduleInTRU;     i++)
    for(int j=0;j<maxnTRU;           j++)
      fADC[j][i]                = 0 ;
  for(int i=0;i<4;                   i++)
    fS[i] = 0 ;
}

///
/// Constructor
//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream(AliRawReader* rawReader) : TObject(),
  fRawReader(rawReader)   ,
  fDetector(kEMCAL)       ,
  fGetRawData(0)          ,
  fPayload(v0)            ,
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
  for(int i=0;i<maxL0GammaPatchIndex;i++) 
    fL0GammaPatchIndex[i]     = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    fNL1GammaPatch   [i]      = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    fL1GammaThreshold[i]      = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    for(int j=0;j<3;                 j++)
      fG[j][i]                  = 0 ;
  for(int i=0;i<maxL1Gamma;          i++)
    for(int j=0;j<maxL1GammaPatchIndex;j++)
      fL1GammaPatchIndex[j][i]  = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    fNL1JetPatch[i]           = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    fL1JetThreshold[i]        = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    for(int j=0;j<3;                 j++)
      fJ[j][i]                  = 0 ;
  for(int i=0;i<maxL1Jet;            i++)
    for(int j=0;j<maxL1JetPatchIndex;j++)
      fL1JetPatchIndex[j][i]    = 0 ;
  for(int i=0;i<maxnmoduleInTRU;     i++)
    for(int j=0;j<maxnTRU;           j++)
      fADC[j][i]                = 0 ;
  for(int i=0;i<4;                   i++)
    fS[i] = 0 ;
}

///
/// Destructor
//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::~AliEMCALTriggerSTURawStream()
{ }

///
/// Reset arrays and raw stream
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::Reset()
{
  if (fRawReader) fRawReader->Reset();
  
  fNL0GammaPatch = 0;
  
  for(int i=0;i<maxL1Gamma;i++)fNL1GammaPatch[i] = 0 ;
  for(int i=0;i<maxL1Jet;  i++)fNL1JetPatch  [i] = 0 ;
}

///
/// Read pay load
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::ReadPayLoad()
{
  UInt_t word32[maxpayloadSize]; // 32b words
  for (Int_t i=0; i<maxpayloadSize; i++) word32[i] = 0;
  
  Int_t iword = 0;
  
  fNL0GammaPatch = 0;
  for(int i=0;i<maxL1Gamma ;i++)    fNL1GammaPatch  [i] = 0 ;
  for(int i=0;i<maxL1Jet   ;i++)    fNL1JetPatch    [i] = 0 ;
  
  Int_t eqId    = -1  ;
  Int_t eqSize  =  0  ;
  
  UInt_t w32;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // migrate raw data to local array
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while (fRawReader->ReadNextInt(w32))
  {
    if (!iword)
    {
      eqId   = fRawReader->GetEquipmentId();
      eqSize = fRawReader->GetEquipmentSize();
    }
    
    if(iword > maxpayloadSize -1 )return kFALSE ;
    word32[iword++] = w32;
    //cout<<dec<<setw(5)<<iword<<" , "<<hex<<setw(10)<<w32<<dec<<setw(10)<<w32<<endl;
  }
  
  Int_t   poffset = 0     ;
  
  //payload type selector
  if(false){}
  else if(iword==(kPayLoadSizeV0                               ))
  {
    poffset         = 0         ;
    fPayload        = v0        ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV0       + kPayLoadSizeV0Raw     ))
  {
    poffset         = 0         ;
    fPayload        = v0Raw     ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1                               ))
  {
    poffset         = 0         ;
    fPayload        = v1        ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV1       + kPayLoadSizeV1Raw     ))
  {
    poffset         = 0         ;
    fPayload        = v1Raw     ; 
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV12                              ))
  {
    poffset         = 0         ;
    fPayload        = v12       ;
    fDetector       = kEMCAL    ;
  }
  else if(iword==(kPayLoadSizeV12     + kPayLoadSizeV12Raw     ))
  {
    poffset         = 0         ;
    fPayload        = v12Raw    ;
    fDetector       = kEMCAL    ;
  }
  
  else if(iword==(kPayLoadSizeV2DCAL                           ))
  {
    poffset         = 0         ;
    fPayload        = v2DCAL    ;
    fDetector       = kDCAL     ;
  }
  
  else if(iword==(kPayLoadSizeV2EMCAL                          ))
  {
    poffset         = 0         ;
    fPayload        = v2EMCAL   ;
    fDetector       = kEMCAL    ;
  }
  
  else if(iword==(kPayLoadSizeV2DCAL  + kPayLoadSizeV2DCALRaw  ))
  {
    poffset         = 0           ;
    fPayload        = v2DCALRaw   ;
    fDetector       = kDCAL       ;
  }
  
  else if(iword==(kPayLoadSizeV2EMCAL + kPayLoadSizeV2EMCALRaw ))
  {
    poffset         = 0           ;
    fPayload        = v2EMCALRaw  ;
    fDetector       = kEMCAL      ;
  }
  
  else
  {
    AliError(Form("STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32", eqId, eqSize, iword));
    return kFALSE;
  }
 
  AliDebug(1, Form("STU (eqId: %d, eqSize: %d) payload size: %d word32", eqId, eqSize, iword));
  
  switch (fPayload)
  {
    case v0         :    case v1         :    case v2EMCAL    :    case v2DCAL     :  case v12    :  {fGetRawData=0;  break;}
    case v0Raw      :    case v1Raw      :    case v2EMCALRaw :    case v2DCALRaw  :  case v12Raw :  {fGetRawData=1;  break;}
    default :    {}
  }//end case
  
  int index         = poffset;
  int offset        = 0;
  int nJetThresh    = 0;
  int nGammaThresh  = 0;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // read setting parameters
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (fPayload)
  {
    //##################################################
    case  v0         :
    case  v0Raw      :
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
    case  v1         :
    case  v1Raw      :
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
    case  v12    :
    case  v12Raw :
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
    case  v2EMCAL    :
    case  v2EMCALRaw :
    {
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][1] = word32[index++] ;
      for(int i=0;i<3;i++) fJ[i][1] = word32[index++] ;
      fRho            = word32[index++]        & 0x000FFFFF   ;
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
    case  v2DCAL     :
    case  v2DCALRaw  :
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
      fL1JetThreshold   [0] = GetThreshold(fJ[0][0], fJ[1][0], fJ[2][0], fV0A, fV0C);
      fL1GammaThreshold [1] = GetThreshold(fG[0][1], fG[1][1], fG[2][1], fV0A, fV0C); 
      fL1JetThreshold   [1] = GetThreshold(fJ[0][1], fJ[1][1], fJ[2][1], fV0A, fV0C);
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
  
  //##############################
  //firmware version confirmation    
  if( 
     (fDetector==kEMCAL && ((fFwVersion & 0xf000)!=0xb000) ) ||
     (fDetector==kDCAL  && ((fFwVersion & 0xf000)!=0xd000) ) 
     )
  {
    return kFALSE;
  }
  //##############################
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // start decoding
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  int L1JetIndexOffset   = (fDetector == kEMCAL)? 11 : 
  (fDetector == kDCAL )? 11 :
  0 ;
  int L0GammaIndexOffset = (fDetector == kEMCAL)? 6 * nTRUEMCAL / 2 : 
  (fDetector == kDCAL )? 6 * nTRUDCAL  / 2 :
  0 ;
  int L1GammaIndexOffset = (fDetector == kEMCAL)? 8 * nTRUEMCAL / 2 : 
  (fDetector == kDCAL )? 8 * nTRUDCAL  / 2 :
  0 ;
  
  //Jet patch index
  for (int i = 0; i < nJetThresh; i++) 
  {
    DecodeL1JetPatchIndexes(i, word32, offset);
    offset += L1JetIndexOffset ;
  }
  
  //L0 Gamma index
  DecodeL0GammaPatchIndexes(word32, offset);
  offset += L0GammaIndexOffset  ;
  
  //L1 Gamma index
  for (int i = 0; i < nGammaThresh; i++) 
  {
    DecodeL1GammaPatchIndexes(i, word32, offset);  
    offset += L1GammaIndexOffset;
  }
  
  //decide if Raw data included 
  if (fGetRawData == 0)
  {
    return kTRUE;
  }
  else
  {
    DecodeTRUADC(word32, offset);
    
    //DCAL : PHOSsubregion data
    if(fDetector == kDCAL){
      offset  += nTRUDCAL * nModDCAL / 2  ;
      DecodePHOSSubregion(word32, offset);
    }
    return kTRUE;
  }
  return kFALSE ;
}

///
/// Decode L0 Gamma patch indeces
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL0GammaPatchIndexes(UInt_t *word32, const int offset)
{
  if(    fPayload == v0         
     ||  fPayload == v0Raw      
     ||  fPayload == v1         
     ||  fPayload == v1Raw      
     ||  fPayload == v12         
     ||  fPayload == v12Raw      
     )
  {
    unsigned short  truL0indexes[maxnTRU][6];
    Int_t nTRU  = GetnTRU();
    
    // extraction from stream
    for (Int_t index=0;index<6;index++) 
    {
      for (Int_t tru_num=0;tru_num<nTRU/2;tru_num++)
      {
        truL0indexes[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
        truL0indexes[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
      }
    }
    
    for (Int_t tru_num=0;tru_num<nTRU;tru_num++)
    {
      for (Int_t index=0;index<6;index++)
      {
        for (Int_t bit_num=0;bit_num<12;bit_num++)
        {
          if ((truL0indexes[tru_num][index] & (1 << bit_num)))
          {
            Int_t idx = 12 * index + bit_num;
            fNL0GammaPatch++;
            fL0GammaPatchIndex[fNL0GammaPatch-1] = (((idx << 5) & 0x7E0) | (tru_num & 0x1F));
          }
        }
      }
    }
    //v2 : L0 index is filled by 0  => no mean
  }
  else
  {
    return ; 
  }
}

///
/// Decode L1 Jet patch indeces
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1JetPatchIndexes(const int i, UInt_t *word32, const int offset)
{
  Int_t nSubregionEta  = (fDetector == kEMCAL)? nSubregionEtaEMCAL :
  (fDetector == kDCAL )? nSubregionEtaDCAL  :
  0 ;
  
  Int_t nSubregionPhi  = (fDetector == kEMCAL)? nSubregionPhiEMCAL :
  (fDetector == kDCAL )? nSubregionPhiDCAL  :
  0 ;
  
  int jetSize = 2 + fPatchSize ; //0->2x2, 2->4x4
  
  for (Int_t ix = 0; ix < nSubregionEta - (jetSize-1) ; ix++)
  {
    UInt_t currentrow = word32[offset + ix];
    for (Int_t iy = 0; iy < nSubregionPhi - (jetSize-1); iy++)
    {
      Int_t bit = (fDetector == kDCAL)? (currentrow & (1 << iy)) : (currentrow & (1<<(nSubregionPhi-jetSize-iy)));
      //if (currentrow & (1 << iy))
      
      if (bit)
      {
        fNL1JetPatch[i]                          = fNL1JetPatch[i] + 1;
        fL1JetPatchIndex[fNL1JetPatch[i] - 1][i] = ((ix << 8) & 0xFF00) | (iy & 0xFF);
      }
    }
  }
  
  return ;
}

///
/// Get max gamma patch 
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1GammaMaxPatch(Int_t& x, Int_t& y, Int_t& z) const
{
    z =   fV0C & 0xF;
    y = ( fV0C >> 4 ) & 0x1F;
    x = ( fV0C >> 9 ) & 0x1F;

    return kTRUE;
}

///
/// Decode L1 gamma patch indeces
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset)
{
  Int_t nTRU  = GetnTRU();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(    fPayload == v0         
     ||  fPayload == v0Raw      
     ||  fPayload == v1         
     ||  fPayload == v1Raw      
     ||  fPayload == v12         
     ||  fPayload == v12Raw      
     )
  {
    unsigned short  truL1indexes0[maxnTRU][8];
    
    // extraction from stream
    for (Int_t index  =0; index   <8      ; index++   )
    {
      for (Int_t tru_num=0; tru_num <nTRU/2 ; tru_num++ )
      {
        truL1indexes0[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
        truL1indexes0[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
      }  
    }
    
    // interpretation
    int gammacolnum;
    short indexcopy;
    for (Int_t tru_num=0; tru_num<nTRU ; tru_num++  )
    {
      for (Int_t index  =0; index  <8    ; index++    )
      {
        for (Int_t bit_num=0; bit_num<12   ; bit_num++  )
        {
          if ((truL1indexes0[tru_num][index] & (1<<bit_num)) != 0)
          {
            if (index<4)
            {
              // Even
              gammacolnum = (2*bit_num  );
              indexcopy   = index;
            }
            else
            { 
              // Odd
              gammacolnum = (2*bit_num+1);
              indexcopy   = index-4;
            }            
            fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
            fL1GammaPatchIndex[fNL1GammaPatch[i] - 1][i] = (
                                                              ((indexcopy   << 10) & 0xC00) 
                                                            | ((gammacolnum <<  5) & 0x3E0) 
                                                            | ( tru_num            &  0x1F)
                                                            );
          } // if
        } // for bit_num
      } // for index
    } // for tru_num
  }//fPayload  
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else if(    fPayload == v2EMCAL
          ||  fPayload == v2EMCALRaw
          ||  fPayload == v2DCAL
          ||  fPayload == v2DCALRaw
          )
  {
    Int_t thirdTRU  = (fDetector == kEMCAL)? 30 : 
    (fDetector == kDCAL )? 12 :
    0 ;
    const Int_t nModPhi = 12  ;
    const Int_t nModEta =  8  ;
    unsigned short  truL1indexes1[maxnTRU][nModPhi];
    for (Int_t iphi=0; iphi<nModPhi/2 ; iphi++)
    {
      for (Int_t itru=0; itru<nTRU/2    ; itru++)
      {
        truL1indexes1[2*itru  ][2*iphi   ] = ( word32[offset + iphi * nTRU/2 + itru] >> 0  & 0xFF);
        truL1indexes1[2*itru  ][2*iphi+1 ] = ( word32[offset + iphi * nTRU/2 + itru] >> 8  & 0xFF);
        truL1indexes1[2*itru+1][2*iphi   ] = ( word32[offset + iphi * nTRU/2 + itru] >> 16 & 0xFF);
        truL1indexes1[2*itru+1][2*iphi+1 ] = ( word32[offset + iphi * nTRU/2 + itru] >> 24 & 0xFF);
      }  
    }
    
    for (Int_t itru=0; itru<nTRU    ; itru++)
    {
      for (Int_t iphi=0; iphi<nModPhi ; iphi++)
      {
        for (Int_t ieta=0; ieta<nModEta ; ieta++)
        {
          if ( (truL1indexes1[itru][iphi] >> ieta) & 0x1 )
          {
            fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
            
            short iphi_tmp  = (itru >= thirdTRU)? ( iphi%2 + 2* (int)(iphi/6)       ) : iphi  ;
            short ieta_tmp  = (itru >= thirdTRU)? ( ieta   + 8*((int)(iphi/2) % 3 ) ) : ieta  ;
            
            fL1GammaPatchIndex[fNL1GammaPatch[i] - 1][i] = (
                                                              ((iphi_tmp  << 10) & 0x7C00) 
                                                            | ((ieta_tmp  <<  5) & 0x03E0) 
                                                            | ((itru      <<  0) & 0x001F)
                                                            );
          } // if
        } // for ieta
      } // for iphi
    } // itru
  } // fPayload  
}

///
/// Decode TRU ADC
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeTRUADC(UInt_t *word32, const int offset)
{
  Int_t nTRU  = GetnTRU();
  Int_t nMod  = GetnMod();
  
  // extraction from stream
  for (Int_t index=0;index<nMod;index++)
  {
    for (Int_t tru_num=0;tru_num<nTRU/2;tru_num++)
    {
      fADC[2*tru_num  ][index] = ( word32[offset + index * nTRU/2 + tru_num]        & 0xFFFF);
      fADC[2*tru_num+1][index] = ((word32[offset + index * nTRU/2 + tru_num] >> 16) & 0xFFFF);
    }
  }  
  
  if(    fPayload == v0Raw         
     ||  fPayload == v1Raw      
     ||  fPayload == v12Raw      
     )
  {
    for (Int_t tru_num=16;tru_num<32;tru_num++)
    {
      // A side
      Int_t v[96];
      for (Int_t index=0;index<96;index++) v[index] = fADC[tru_num][95-index];
      for (Int_t index=0;index<96;index++) fADC[tru_num][index] = v[index];
    }
  }
}

///
/// Decode PHOS sub-regions
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodePHOSSubregion(UInt_t *word32, const int offset)
{
  for (Int_t index=0;index<36;index++)
  {
    fPHOSSubregion[index] = word32[offset + index] & 0xFFFFFFFF  ;
  }
}

///
/// L0 gamma patches sent to STU (original access to L0 patch indexes)
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL0GammaPatch(const Int_t i, Int_t& tru, Int_t& idx) const
{  
  if (i > fNL0GammaPatch) return kFALSE;
  
  tru =  fL0GammaPatchIndex[i] & 0x1F;
  idx = (fL0GammaPatchIndex[i] & 0x7E0) >> 5;
  
  return kTRUE;
}

///
/// L1 gamma patch indexes
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1GammaPatch(const Int_t i, const Int_t j, 
                                                    Int_t& tru, Int_t& col, Int_t& row) const
{
  if (j >= maxL1Gamma || i > fNL1GammaPatch[j]) return kFALSE;
  
  row = (fL1GammaPatchIndex[i][j] >> 10) & 0x1F ;
  col = (fL1GammaPatchIndex[i][j] >>  5) & 0x1F ;
  tru = (fL1GammaPatchIndex[i][j] >>  0) & 0x1F ;
  
  return kTRUE;
}

///
/// L1 jet patch indexes
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1JetPatch(const Int_t i, const Int_t j,
                                                  Int_t& col, Int_t& row) const
{
  if (j >= maxL1Jet || i > fNL1JetPatch[j]) return kFALSE;
  
  col = (fL1JetPatchIndex[i][j] >> 0) & 0xFF  ;
  row = (fL1JetPatchIndex[i][j] >> 8) & 0xFF  ;
  
  return kTRUE;
}

///
/// Time sums
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::GetADC(Int_t iTRU, UInt_t ADC[])
{
  Int_t nMod  = GetnMod() ;
  for (Int_t i=0; i<nMod; i++) ADC[i] = fADC[iTRU][i];
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::GetPHOSSubregion(UInt_t PHOSSubregion[])
{
  for (Int_t i=0; i<36; i++) PHOSSubregion[i] = fPHOSSubregion[i];
}

///
/// Dump STU payload
//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DumpPayLoad(const Option_t *option) const
{  
  TString op = option;
  
  printf("V0A:             %d\n"  , fV0A              );
  printf("V0C:             %d\n"  , fV0C              );
  
  printf("Number of L0:    %d\n"  , fNL0GammaPatch    );
  for (int i = 0; i < maxL1Gamma; i++) 
  {
    for (int j = 0; j < 3         ; j++) 
    {
      printf("G[%d][%d]: %d\n", j, i, fG[j][i]);
    }
    
    printf("Gamma threshold[%d]: %d\n", i, fL1GammaThreshold[i]);
    printf("Number of L1-g[%d]: %d\n", i, fNL1GammaPatch[i]);
  }
  
  for (int i = 0; i < maxL1Jet; i++) 
  {
    for (int j = 0; j < 3       ; j++) 
    {
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
  if (op.Contains("L0") || op.Contains("ALL"))
  {
    for (Int_t i = 0;i < fNL0GammaPatch; i++) 
    {
      if (GetL0GammaPatch(i,itru,col)) 
        cout  << dec
        << "> Found L0 gamma in TRU #"  << setw(2) << itru 
        << " at idx: "                  << setw(2) << col 
        << endl;
    }
  }
  // 
  if (op.Contains("L1") || op.Contains("ALL")) 
  {
    for (int j = 0; j < maxL1Gamma; j++)
    {
      for (Int_t i = 0; i < fNL1GammaPatch[j]; i++) 
      {
        if (GetL1GammaPatch(i, j, itru, col, row)) 
          cout  << dec
          << "> Found L1 gamma "  << j 
          << " in TRU #"          << setw(2) << itru 
          << " at: ( col: "       << setw(2) << col 
          << " , row: "           << setw(2) << row 
          << " )" << endl;
      }
    }
    
    for (int j = 0; j < maxL1Jet; j++)
    {
      for (Int_t i = 0; i < fNL1JetPatch[j]; i++) 
      {
        if (GetL1JetPatch(i, j, col, row)) 
          cout  << dec 
          << "> Found L1 jet "  << j 
          << " at: ( col: "     << setw(2) << col 
          << " , row: "         << setw(2) << row 
          << " )" << endl;
      }
    }
  }
  
  if ((op.Contains("ADC") || op.Contains("ALL")) && fGetRawData) 
  {
    for (Int_t i = 0; i < GetnTRU(); i++) 
    {
      cout << "--------\n";
      cout << "TRU #" << setw(2) << dec << i << ":";
      for (Int_t j = 0;j < GetnMod(); j++) 
      { 
        //TBits xadc(12); xadc.Set(12,&fADC[i][j]); 
        if ((j % 4) == 0) cout << endl;
        //cout << setw(2) << j << ": " << xadc << " ";
        printf("%2d: %3x / ",j,fADC[i][j]); 
      }
      cout << "\n";
    }
  }
}

///
/// Get threshold 
//_____________________________________________________________________________
UShort_t AliEMCALTriggerSTURawStream::GetThreshold(Short_t a, Short_t b, Short_t c, 
                                                   UShort_t v0A, UShort_t v0C) const
{
  ULong64_t v0sum = v0A + v0C;
  ULong64_t sqrV0 = v0sum * v0sum;          
  
  sqrV0 *= a;
  sqrV0 >>= 32;
  
  v0sum *= b;
  v0sum >>= 16;
  
  return (UShort_t)(sqrV0 + v0sum + c);
}
