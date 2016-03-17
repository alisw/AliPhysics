#ifndef ALIEMCALTRIGGERSTURAWSTREAM_H
#define ALIEMCALTRIGGERSTURAWSTREAM_H
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

class AliEMCALTriggerSTURawStream: public TObject 
{
  public:
    AliEMCALTriggerSTURawStream();
    AliEMCALTriggerSTURawStream(AliRawReader* rawReader);
    virtual ~AliEMCALTriggerSTURawStream();


    virtual void    SetDetector(int det){fDetector = static_cast<fDetType>(det);}
    virtual Int_t   GetDetector(){return static_cast<Int_t>(fDetector) ;}

    virtual void    Reset();
    virtual Bool_t  ReadPayLoad();
    
    virtual void    DumpPayLoad(const Option_t *option = "ALL") const;

    virtual void    GetADC(           Int_t iTRU, UInt_t ADC[]            );
    virtual void    GetPHOSSubregion(             UInt_t PHOSSubregion[]  );
    virtual UInt_t  GetL1JetThreshold(  const int i)  const {return fL1JetThreshold  [i];}
    virtual UInt_t  GetL1GammaThreshold(const int i)  const {return fL1GammaThreshold[i];}

    virtual Int_t   GetNL0GammaPatch()                const {return fNL0GammaPatch      ;}
    virtual Int_t   GetNL1GammaPatch(   const int i)  const {return fNL1GammaPatch[i]   ;}
    virtual Int_t   GetNL1JetPatch(     const int i)  const {return fNL1JetPatch  [i]   ;}

    virtual Bool_t  GetL0GammaPatch(const Int_t i,                Int_t& x, Int_t& y          ) const;
    virtual Bool_t  GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& x, Int_t& y, Int_t& z) const;
    virtual Bool_t  GetL1JetPatch(  const Int_t i, const Int_t j, Int_t& x, Int_t& y          ) const;

    virtual UInt_t  GetV0A()            const {return fV0A           ;}
    virtual UInt_t  GetV0C()            const {return fV0C           ;}
    virtual UInt_t  GetG(int i, int j)  const {return fG[i][j]       ;}//[ABC][high/low]
    virtual UInt_t  GetJ(int i, int j)  const {return fJ[i][j]       ;}//[ABC][high/low]
    virtual UInt_t  GetPHOSScale(int i) const {return fS[i]          ;}
    virtual UInt_t  GetRho()            const {return fRho           ;}
    virtual UInt_t  GetRegionEnable()   const {return fRegionEnable  ;}
    virtual UInt_t  GetFrameReceived()  const {return fFrameReceived ;}
    virtual UInt_t  GetFwVersion()      const {return fFwVersion     ;}
    virtual UInt_t  GetPatchSize()      const {return fPatchSize     ;}
    virtual Int_t   GetRawData()        const {return fGetRawData    ;}
    virtual Int_t   GetnTRU()           const {return (fDetector==kEMCAL)?nTRU_EMCAL:(fDetector==kDCAL)?nTRU_DCAL:0;}

  private:

    AliEMCALTriggerSTURawStream(const AliEMCALTriggerSTURawStream& rhs);
    AliEMCALTriggerSTURawStream& operator = (const AliEMCALTriggerSTURawStream& rhs);

    UShort_t GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const;

    AliRawReader* fRawReader;   // object for reading the raw data

    enum fDetType { 
      kEMCAL  = 0, 
      kDCAL   = 1 
    }; 
    
    enum fPayloadType { 
      V0          =   0 ,
      V0Raw       =   1 ,
      V1          =   2 ,
      V1Raw       =   3 ,
      V1_2        =   8 ,
      V1_2Raw     =   9 ,

      V2EMCAL     =   4 ,  
      V2EMCALRaw  =   5 ,
      V2DCAL      =   6 , 
      V2DCALRaw   =   7 ,
      def         =   8
    }; 

    fDetType      fDetector   ;//Set by function
    Int_t         fGetRawData ;//Set by word size
    fPayloadType  fPayload    ;//Set by word size

    static const Int_t kPayLoadSizeV0             = 236                         ;
    static const Int_t kPayLoadSizeV0_Raw         = 1536                        ;
    static const Int_t kPayLoadSizeV1             = 245                         ;
    static const Int_t kPayLoadSizeV1_Raw         = 1536                        ;
    //XXX
    static const Int_t kPayLoadSizeV1_2           = 390                         ;
    static const Int_t kPayLoadSizeV1_2_Raw       = 1536                        ;

    static const Int_t kPayLoadSizeV2_EMCAL       = 391   ;//17+11+11+96+128+128    
    static const Int_t kPayLoadSizeV2_EMCAL_Raw   = 1536  ;//96*32/2                   
    static const Int_t kPayLoadSizeV2_DCAL        = 197   ;//21+11+11+42+ 56+ 56    
    static const Int_t kPayLoadSizeV2_DCAL_Raw    = 708   ;//96*14/2 +36               

    static const Int_t max_payload_size = 1937 ;

    static const Int_t max_L0GammaPatchIndex  = 3100 ;  // (48-1)*(64-1)  = 2961  (EMCAL)
                                                        // (48-1)*(40-1)  = 1833  (DCAL)
    static const Int_t max_L1Gamma            =    2 ;  // L1Gamma_low,high       (EMCAL)
    static const Int_t max_L1GammaPatchIndex  = 3100 ;  // (48-1)*(64-1)  = 2961  (EMCAL)
                                                        // (48-1)*(40-1)  = 1833  (DCAL)
    static const Int_t max_L1Jet              =    2 ;  // L1Jet_low,high         (EMCAL,DCAL)
    static const Int_t max_L1JetPatchIndex    =  200 ;  // (12-1)*(16-1)  = 165   (EMCAL)
                                                        // (12-1)*(10-1)  =  99   (DCAL)

    static const Int_t max_nTRU               =   32 ;  // 32  (EMCAL)
                                                        // 14  (DCAL)
    static const Int_t max_nmoduleInTRU       =   96 ;  // 96  (EMCAL)
                                                        // 96  (DCAL)
    
    static const Int_t nTRU_EMCAL             = 32  ;
    static const Int_t nTRU_DCAL              = 14  ;

    static const Int_t nMod_EMCAL             = 96  ;
    static const Int_t nMod_DCAL              = 96  ;

    static const Int_t nSubregion_eta_EMCAL   = 12  ;
    static const Int_t nSubregion_phi_EMCAL   = 16  ;

    static const Int_t nSubregion_eta_DCAL    = 12  ;
    static const Int_t nSubregion_phi_DCAL    = 10  ;


    Int_t     fNL0GammaPatch                                        ; // N L0 Gamma Patch
    UShort_t  fL0GammaPatchIndex[max_L0GammaPatchIndex]             ; // L0 Gamma Patch Index

    Int_t     fNL1GammaPatch                           [max_L1Gamma]; // N L1 Gamma Patch
    UInt_t    fG                                    [3][max_L1Gamma]; // Gamma threshold parameter:A,B,C
    UInt_t    fL1GammaThreshold                        [max_L1Gamma]; // L1 Gamma Threshold
    UShort_t  fL1GammaPatchIndex[max_L1GammaPatchIndex][max_L1Gamma]; // L1 Gamma Patch Index
    
    Int_t     fNL1JetPatch                         [max_L1Jet]      ; // N L1 Jet Patch
    UInt_t    fJ                                [3][max_L1Jet]      ; // Jet threshold parameter:A,B,C
    UInt_t    fL1JetThreshold                      [max_L1Jet]      ; // L1 Jet Threshold
    UShort_t  fL1JetPatchIndex[max_L1JetPatchIndex][max_L1Jet]      ; // L1 Jet Patch Index

    UInt_t    fADC[max_nTRU][max_nmoduleInTRU];
    UInt_t    fPHOSSubregion[36];

    UInt_t    fV0A        ;      // V0A
    UInt_t    fV0C        ;      // V0C
    UInt_t    fS[4]       ;      // PHOS Scale parameter
    UInt_t    fRho        ;      // background Rho
    UInt_t    fPatchSize  ;      // jet patch size

    UInt_t    fRegionEnable   ;  // Region Enable
    UInt_t    fFrameReceived  ;  // Frame Received
    UInt_t    fFwVersion      ;  // Fw Version

    Int_t     GetnMod()const{ return  (fDetector == kEMCAL)? nMod_EMCAL :
                                      (fDetector == kDCAL )? nMod_DCAL  :
                                      0 ;
    }



    virtual void    DecodeL1JetPatchIndexes(  const int i, UInt_t *word32, const int offset);
    virtual void    DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset);
    virtual void    DecodeL0GammaPatchIndexes(             UInt_t *word32, const int offset);
    virtual void    DecodeTRUADC(                          UInt_t *word32, const int offset);
    virtual void    DecodePHOSSubregion(                   UInt_t *word32, const int offset);


    ClassDef(AliEMCALTriggerSTURawStream,2)   // class for reading EMCAL STU DDL raw data
};

#endif
