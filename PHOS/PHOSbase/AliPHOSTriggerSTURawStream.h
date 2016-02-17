#ifndef ALIPHOSTRIGGERSTURAWSTREAM_H
#define ALIPHOSTRIGGERSTURAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
Author: H. Yokoyama (Univ. of TSUKUBA/ Univ. of Grenoble)
*/

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
//#include <map>

class AliRawReader;

class AliPHOSTriggerSTURawStream: public TObject 
{
  public:
    AliPHOSTriggerSTURawStream();
    AliPHOSTriggerSTURawStream(AliRawReader* rawReader);
    virtual ~AliPHOSTriggerSTURawStream();

    virtual void    Reset();
    virtual Bool_t  ReadPayLoad();

    virtual void    DumpPayLoad(const Option_t *option = "ALL") const;

    virtual void    GetADC(           Int_t iTRU, UInt_t ADC[]            );
    virtual UInt_t  GetL1GammaThreshold(const int i)  const {return fL1GammaThreshold[i];}

    virtual Int_t   GetNL0GammaPatch()                const {return fNL0GammaPatch      ;}
    virtual Int_t   GetNL1GammaPatch(   const int i)  const {return fNL1GammaPatch[i]   ;}

    virtual Bool_t  GetL0GammaPatch(const Int_t i,                Int_t& tru, Int_t& idx            ) const;
    virtual Bool_t  GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& tru, Int_t& col, Int_t& row) const;

    virtual UInt_t  GetV0A()            const {return fV0A           ;}
    virtual UInt_t  GetV0C()            const {return fV0C           ;}
    virtual UInt_t  GetG(int i, int j)  const {return fG[i][j]       ;}//[ABC][high/low]
    virtual UInt_t  GetRegionEnable()   const {return fRegionEnable  ;}
    virtual UInt_t  GetFrameReceived()  const {return fFrameReceived ;}
    virtual UInt_t  GetFwVersion()      const {return fFwVersion     ;}
    virtual Int_t   GetRawData()        const {return fGetRawData    ;}

  private:

    AliPHOSTriggerSTURawStream(const AliPHOSTriggerSTURawStream& rhs);
    AliPHOSTriggerSTURawStream& operator = (const AliPHOSTriggerSTURawStream& rhs);

    UShort_t GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const;

    AliRawReader* fRawReader;   // object for reading the raw data

    enum fPayloadType { 
      V2PHOS      =   8 , 
      V2PHOSRaw   =   9  
    }; 

    Int_t         fGetRawData ;//Set by word size
    fPayloadType  fPayload    ;//Set by word size

    static const Int_t kPayLoadSizeV2_PHOS        = 433   ;//13+ 0+ 0+84+112+112+112
    static const Int_t kPayLoadSizeV2_PHOS_Raw    = 1568  ;//112*28/2                  

    //static const Int_t max_payload_size = kPayLoadSizeV2_PHOS + kPayLoadSizeV2_PHOS_Raw ;
    static const Int_t max_payload_size = 2001 ;

    static const Int_t max_L0GammaPatchIndex  = 3100 ;  // (28-1)*(112-1) = 2997  (PHOS)
    static const Int_t max_L1Gamma            =    3 ;  // L1Gamma_low,mid,high   (PHOS)
    static const Int_t max_L1GammaPatchIndex  = 3100 ;  // (28-1)*(112-1) = 2997  (PHOS)

    static const Int_t max_nTRU               =   28 ;  // 28  (PHOS)
    static const Int_t max_nmoduleInTRU       =  112 ;  // 112 (PHOS)
    
    static const Int_t nTRU_PHOS              = 28  ;
    static const Int_t nMod_PHOS              = 112 ;

    Int_t     fNL0GammaPatch                                        ; // N L0 Gamma Patch
    UShort_t  fL0GammaPatchIndex[max_L0GammaPatchIndex]             ; // L0 Gamma Patch Index

    Int_t     fNL1GammaPatch                           [max_L1Gamma]; // N L1 Gamma Patch
    UInt_t    fG                                    [3][max_L1Gamma]; // Gamma threshold parameter:A,B,C
    UInt_t    fL1GammaThreshold                        [max_L1Gamma]; // L1 Gamma Threshold
    UShort_t  fL1GammaPatchIndex[max_L1GammaPatchIndex][max_L1Gamma]; // L1 Gamma Patch Index
    
    UInt_t    fADC[max_nTRU][max_nmoduleInTRU];

    UInt_t    fV0A        ;      // V0A
    UInt_t    fV0C        ;      // V0C

    UInt_t    fRegionEnable   ;  // Region Enable
    UInt_t    fFrameReceived  ;  // Frame Received
    UInt_t    fFwVersion      ;  // Fw Version

    Int_t     GetnTRU()const{ return  nTRU_PHOS;}
    Int_t     GetnMod()const{ return  nMod_PHOS;}
    
    virtual void    DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset);
    virtual void    DecodeL0GammaPatchIndexes(             UInt_t *word32, const int offset);
    virtual void    DecodeTRUADC(                          UInt_t *word32, const int offset);


    ClassDef(AliPHOSTriggerSTURawStream,2)   // class for reading PHOS STU DDL raw data
};

#endif
