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

//_________________________________________________________________________
//  This class provides access to PHOS STU DDL raw data.
//--
//  Author      : Hiroki Yokoyama (Univ. of TSUKUBA / Univ. of Grenoble)
//  contact     : hiroki.yokoyama@cern.ch
//  Last Update : 23. Nov. 2015
//_________________________________________________________________________

#include "AliPHOSTriggerSTURawStream.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TBits.h"

#include <cstdlib>


using std::cout;
using std::setw;
using std::endl;
ClassImp(AliPHOSTriggerSTURawStream)

//_____________________________________________________________________________
AliPHOSTriggerSTURawStream::AliPHOSTriggerSTURawStream() : TObject(),
  fRawReader(0x0)         ,
  fGetRawData(0)          ,
  fPayload(V2PHOS)        ,
  fNL0GammaPatch(0)       ,
  fL0GammaPatchIndex()    ,
  fNL1GammaPatch()        ,
  fG()                    ,
  fL1GammaThreshold()     ,
  fL1GammaPatchIndex()    ,
  fADC()                  ,
  fV0A(0)                 ,
  fV0C(0)                 ,
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
  for(int i=0;i<max_nmoduleInTRU;     i++)
  for(int j=0;j<max_nTRU;             j++)
    fADC[j][i]                = 0 ;
}

//_____________________________________________________________________________
AliPHOSTriggerSTURawStream::AliPHOSTriggerSTURawStream(AliRawReader* rawReader) : TObject(),
  fRawReader(rawReader)   ,
  fGetRawData(0)          ,
  fPayload(V2PHOS)        ,
  fNL0GammaPatch(0)       ,
  fL0GammaPatchIndex()    ,
  fNL1GammaPatch()        ,
  fG()                    ,
  fL1GammaThreshold()     ,
  fL1GammaPatchIndex()    ,
  fADC()                  ,
  fV0A(0)                 ,
  fV0C(0)                 ,
  fRegionEnable(0)        ,
  fFrameReceived(0)       ,
  fFwVersion(0)     
{
  //
  fRawReader->Reset();
  fRawReader->Select("PHOS",20);
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
  for(int i=0;i<max_nmoduleInTRU;     i++)
  for(int j=0;j<max_nTRU;             j++)
    fADC[j][i]                = 0 ;
}

//_____________________________________________________________________________
AliPHOSTriggerSTURawStream::~AliPHOSTriggerSTURawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliPHOSTriggerSTURawStream::Reset()
{
  // Reset
  if (fRawReader) fRawReader->Reset();
  fNL0GammaPatch = 0;
  for(int i=0;i<max_L1Gamma;i++)fNL1GammaPatch[i] = 0 ;
}

//_____________________________________________________________________________
Bool_t AliPHOSTriggerSTURawStream::ReadPayLoad()
{
  UInt_t word32[max_payload_size]; // 32b words
  for (Int_t i=0; i<max_payload_size; i++) 
    word32[i] = 0;
  
  Int_t iword = 0;
  
  fNL0GammaPatch = 0;
  for(int i=0;i<max_L1Gamma;i++)
    fNL1GammaPatch  [i] = 0 ;
  
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
  }

  //payload type selector
  if(false){}
  else if(iword==(kPayLoadSizeV2_PHOS                            )){fPayload = V2PHOS     ;}
  else if(iword==(kPayLoadSizeV2_PHOS  + kPayLoadSizeV2_PHOS_Raw )){fPayload = V2PHOSRaw  ;}
  else{
    AliError(Form("STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32", eqId, eqSize, iword));
    return kFALSE;
  }
  AliDebug(1, Form("STU (eqId: %d, eqSize: %d) payload size: %d word32", eqId, eqSize, iword));

  switch (fPayload){
    case V2PHOS     :    {fGetRawData=0;  break;}
    case V2PHOSRaw  :    {fGetRawData=1;  break;}
    default :    {}
  }//end case

  int offset        = 0;
  int nGammaThresh  = 0;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // read setting parameters
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (fPayload){
    //##################################################
    case  V2PHOS     :
    case  V2PHOSRaw  :
    {
      int index = 0 ;
      fV0A = ((word32[index]>>16) & 0xFFFF);
      fV0C = ( word32[index]      & 0xFFFF);
      index++ ;
      for(int i=0;i<3;i++) fG[i][0] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][1] = word32[index++] ;
      for(int i=0;i<3;i++) fG[i][2] = word32[index++] ;
      fRegionEnable   = word32[index++] ;
      fFrameReceived  = word32[index++] ;
      fFwVersion      = word32[index++] ;
      offset = index  ;
      
      fL1GammaThreshold [0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
      fL1GammaThreshold [1] = GetThreshold(fG[0][1], fG[1][1], fG[2][1], fV0A, fV0C); 
      fL1GammaThreshold [2] = GetThreshold(fG[0][2], fG[1][2], fG[2][2], fV0A, fV0C); 
      nGammaThresh  = 3;
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

  int L0GammaIndexOffset = 6 * nTRU_PHOS  / 2 ;
  int L1GammaIndexOffset = 8 * nTRU_PHOS  / 2 ;

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
    return kTRUE;
  }
  return kFALSE ;
}

//_____________________________________________________________________________
void AliPHOSTriggerSTURawStream::DecodeL0GammaPatchIndexes(UInt_t *word32, const int offset)
{
  return ; 
}

//_____________________________________________________________________________
void AliPHOSTriggerSTURawStream::DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset)
{
  Int_t nTRU  = GetnTRU();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(     fPayload == V2PHOS
      ||  fPayload == V2PHOSRaw
      )
  {
    const Int_t nModPhi =  8   ;
    const Int_t nModEta =  14  ;
    unsigned short  truL1indexes2[max_nTRU][nModPhi];
    for (Int_t iphi=0; iphi<nModPhi   ; iphi++)
    for (Int_t itru=0; itru<nTRU/2    ; itru++){
      truL1indexes2[2*itru  ][iphi ] = ( word32[offset + iphi * nTRU/2 + itru] >> 0  & 0x3FFF);
      truL1indexes2[2*itru+1][iphi ] = ( word32[offset + iphi * nTRU/2 + itru] >> 16 & 0x3FFF);
    }  
    for (Int_t itru=0; itru<nTRU    ; itru++)
    for (Int_t iphi=0; iphi<nModPhi ; iphi++)
    for (Int_t ieta=0; ieta<nModEta ; ieta++){
      if ( (truL1indexes2[itru][iphi] >> ieta) & 0x1 ){
        fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
        
        short iphi_tmp  = iphi  ;
        short ieta_tmp  = ieta  ;

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
void AliPHOSTriggerSTURawStream::DecodeTRUADC(UInt_t *word32, const int offset)
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
}

//_____________________________________________________________________________
Bool_t AliPHOSTriggerSTURawStream::GetL0GammaPatch(const Int_t i, Int_t& tru, Int_t& idx) const
{
  // L0 gamma patches sent to STU (original access to L0 patch indexes)
  
  if (i > fNL0GammaPatch) return kFALSE;
  
  tru =  fL0GammaPatchIndex[i] & 0x1F;
  idx = (fL0GammaPatchIndex[i] & 0x7E0) >> 5;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPHOSTriggerSTURawStream::GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& tru, Int_t& col, Int_t& row) const
{
  // L1 gamma patch indexes
  if (j >= max_L1Gamma || i > fNL1GammaPatch[j]) return kFALSE;
  row = (fL1GammaPatchIndex[i][j] >> 10) & 0x1F ;
  col = (fL1GammaPatchIndex[i][j] >>  5) & 0x1F ;
  tru = (fL1GammaPatchIndex[i][j] >>  0) & 0x1F ;
  return kTRUE;
}

//_____________________________________________________________________________
void AliPHOSTriggerSTURawStream::GetADC(Int_t iTRU, UInt_t ADC[])
{
  // Time sums
  Int_t nMod  = GetnMod() ;
  for (Int_t i=0; i<nMod; i++) ADC[i] = fADC[iTRU][i];
}

//_____________________________________________________________________________
void AliPHOSTriggerSTURawStream::DumpPayLoad(const Option_t *option) const
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

  printf("RegionEnable:    %8x\n" , fRegionEnable     ); 
  printf("FrameReceived:   %8x\n" , fFrameReceived    );
  printf("FwVersion:       %x\n"  , fFwVersion        );
  printf("RawData:         %d\n"  , fGetRawData       );

  Int_t itru, col, row;
  // 
  if (op.Contains("L0") || op.Contains("ALL")){
    for (Int_t i = 0;i < fNL0GammaPatch; i++) {
      if (GetL0GammaPatch(i,itru,col)) 
        cout << "> Found L0 gamma in TRU #" << setw(2) << itru <<  " at idx: " << setw(2) << col << endl;
    }
  }
  // 
  if (op.Contains("L1") || op.Contains("ALL")) {
    for (int j = 0; j < max_L1Gamma; j++) {
      for (Int_t i = 0; i < fNL1GammaPatch[j]; i++) {
        if (GetL1GammaPatch(i, j, itru, col, row)) 
          cout << "> Found L1 gamma " << j << " in TRU #" << setw(2) << itru <<  " at: ( col: " << setw(2) << col << " , row: " << setw(2) << row << " )" << endl;
      }
    }
  }
  
  if ((op.Contains("ADC") || op.Contains("ALL")) && fGetRawData) {
    for (Int_t i = 0; i < GetnTRU(); i++) {
      cout << "--------\n";
      cout << "TRU #" << setw(2) << i << ":";
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
UShort_t AliPHOSTriggerSTURawStream::GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const
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
