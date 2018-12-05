#ifndef ALIEMCALTRIGGERSTURAWSTREAM_H
#define ALIEMCALTRIGGERSTURAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//________________________________________________
/// \class AliEMCALTriggerSTURawStream
/// \ingroup EMCALbase
/// \brief provides access to EMCAL/DCAL STU DDL raw data.
///
/// This class provides access to EMCAL/DCAL STU DDL raw data.
///
/// \author: Hiroki Yokoyama, <hiroki.yokoyama@cern.ch>, (Univ. of TSUKUBA / Univ. of Grenoble)
//________________________________________________

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliRawReader;

class AliEMCALTriggerSTURawStream: public TObject 
{
  
public:
  
  AliEMCALTriggerSTURawStream();
  AliEMCALTriggerSTURawStream(AliRawReader* rawReader);
  virtual ~AliEMCALTriggerSTURawStream();
  
  virtual void    SetDetector(int det) { fDetector = static_cast<fDetType>(det) ; }
  virtual Int_t   GetDetector()        { return static_cast<Int_t>(fDetector)   ; }
  
  virtual void    Reset();
  virtual Bool_t  ReadPayLoad();
  
  virtual void    DumpPayLoad(const Option_t *option = "ALL") const;
  
  virtual void    GetADC(           Int_t iTRU, UInt_t ADC[]            );
  virtual void    GetPHOSSubregion(             UInt_t PHOSSubregion[]  );
  virtual UInt_t  GetL1JetThreshold(  const int i)  const { return fL1JetThreshold  [i] ; }
  virtual UInt_t  GetL1GammaThreshold(const int i)  const { return fL1GammaThreshold[i] ; }
  
  virtual Int_t   GetNL0GammaPatch()                const { return fNL0GammaPatch       ; }
  virtual Int_t   GetNL1GammaPatch(   const int i)  const { return fNL1GammaPatch[i]    ; }
  virtual Int_t   GetNL1JetPatch(     const int i)  const { return fNL1JetPatch  [i]    ; }
  
  virtual Bool_t  GetL0GammaPatch(const Int_t i,                Int_t& x, Int_t& y          ) const;
  virtual Bool_t  GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& x, Int_t& y, Int_t& z) const;
  virtual Bool_t  GetL1GammaMaxPatch(                           Int_t& x, Int_t& y, Int_t& z) const; // TRU, X, Y
  virtual Bool_t  GetL1JetPatch(  const Int_t i, const Int_t j, Int_t& x, Int_t& y          ) const;
  
  virtual UInt_t  GetV0A()            const { return fV0A            ; }
  virtual UInt_t  GetV0C()            const { return fV0C            ; }
  virtual UInt_t  GetG(int i, int j)  const { return fG[i][j]        ; }//[ABC][high/low]
  virtual UInt_t  GetJ(int i, int j)  const { return fJ[i][j]        ; }//[ABC][high/low]
  virtual UInt_t  GetPHOSScale(int i) const { return fS[i]           ; }
  virtual UInt_t  GetRho()            const { return (fRho & 0x3FFFF); }
  virtual UInt_t  GetRegionEnable()   const { return fRegionEnable   ; }
  virtual UInt_t  GetFrameReceived()  const { return fFrameReceived  ; }
  virtual UInt_t  GetFwVersion()      const { return fFwVersion      ; }
  virtual UInt_t  GetPatchSize()      const { return fPatchSize      ; }
  virtual Int_t   GetRawData()        const { return fGetRawData     ; }
  virtual Int_t   GetnTRU()           const { return (fDetector==kEMCAL)?nTRUEMCAL:(fDetector==kDCAL)?nTRUDCAL:0 ; }
  
  virtual Int_t   GetEMCALDCALFrameReceived() const { return (fRho >> 18  & 0x3); } // 1b per port
  
private:
  
  AliEMCALTriggerSTURawStream             (const AliEMCALTriggerSTURawStream& rhs);
  AliEMCALTriggerSTURawStream& operator = (const AliEMCALTriggerSTURawStream& rhs);
  
  UShort_t GetThreshold(Short_t a, Short_t b, Short_t c, UShort_t v0A, UShort_t v0C) const;
  
  AliRawReader* fRawReader;   ///< Object for reading the raw data
  
  enum fDetType 
  { 
    kEMCAL  = 0, 
    kDCAL   = 1 
  }; 
  
  enum fPayloadType 
  { 
    v0          =   0 ,
    v0Raw       =   1 ,
    v1          =   2 ,
    v1Raw       =   3 ,
    v12         =   8 ,
    v12Raw      =   9 ,
    
    v2EMCAL     =   4 ,  
    v2EMCALRaw  =   5 ,
    v2DCAL      =   6 , 
    v2DCALRaw   =   7 ,
    def         =   8
  }; 
  
  fDetType      fDetector   ; ///< Set by function
  Int_t         fGetRawData ; ///< Set by word size
  fPayloadType  fPayload    ; ///< Set by word size
  
  static const Int_t kPayLoadSizeV0          = 236  ;
  static const Int_t kPayLoadSizeV0Raw       = 1536 ;
  static const Int_t kPayLoadSizeV1          = 245  ;
  static const Int_t kPayLoadSizeV1Raw       = 1536 ;
  //XXX
  static const Int_t kPayLoadSizeV12         = 390  ;
  static const Int_t kPayLoadSizeV12Raw      = 1536 ;
  
  static const Int_t kPayLoadSizeV2EMCAL     = 391  ; ///< 17+11+11+96+128+128    
  static const Int_t kPayLoadSizeV2EMCALRaw  = 1536 ; ///< 96*32/2                   
  static const Int_t kPayLoadSizeV2DCAL      = 197  ; ///< 21+11+11+42+ 56+ 56    
  static const Int_t kPayLoadSizeV2DCALRaw   = 708  ; ///< 96*14/2 +36               
  
  static const Int_t maxpayloadSize          = 1937 ;
  
  static const Int_t maxL0GammaPatchIndex    = 3100 ; ///< (48-1)*(64-1)  = 2961  (EMCAL)
                                                      ///< (48-1)*(40-1)  = 1833  (DCAL)
  static const Int_t maxL1Gamma              =    2 ; ///< L1Gamma_low,high       (EMCAL)
  static const Int_t maxL1GammaPatchIndex    = 3100 ; ///< (48-1)*(64-1)  = 2961  (EMCAL)
                                                      ///< (48-1)*(40-1)  = 1833  (DCAL)
  static const Int_t maxL1Jet                =    2 ; ///< L1Jet_low,high         (EMCAL,DCAL)
  static const Int_t maxL1JetPatchIndex      =  200 ; ///< (12-1)*(16-1)  = 165   (EMCAL)
                                                      ///< (12-1)*(10-1)  =  99   (DCAL)
  
  static const Int_t maxnTRU                 =   32 ; ///< 32  (EMCAL)
                                                      ///< 14  (DCAL)
  static const Int_t maxnmoduleInTRU         =   96 ; ///< 96  (EMCAL)
                                                      ///< 96  (DCAL)
  
  static const Int_t nTRUEMCAL               = 32   ;
  static const Int_t nTRUDCAL                = 14   ;
  
  static const Int_t nModEMCAL               = 96   ;
  static const Int_t nModDCAL                = 96   ;
  
  static const Int_t nSubregionEtaEMCAL      = 12   ;
  static const Int_t nSubregionPhiEMCAL      = 16   ;
  
  static const Int_t nSubregionEtaDCAL       = 12   ;
  static const Int_t nSubregionPhiDCAL       = 10   ;
  
  
  Int_t     fNL0GammaPatch                                      ; ///< N L0 Gamma Patch
  UShort_t  fL0GammaPatchIndex[maxL0GammaPatchIndex]            ; ///< L0 Gamma Patch Index
  
  Int_t     fNL1GammaPatch                          [maxL1Gamma]; ///< N L1 Gamma Patch
  UInt_t    fG                                   [3][maxL1Gamma]; ///< Gamma threshold parameter:A,B,C
  UInt_t    fL1GammaThreshold                       [maxL1Gamma]; ///< L1 Gamma Threshold
  UShort_t  fL1GammaPatchIndex[maxL1GammaPatchIndex][maxL1Gamma]; ///< L1 Gamma Patch Index
  
  Int_t     fNL1JetPatch                        [maxL1Jet]      ; ///< N L1 Jet Patch
  UInt_t    fJ                               [3][maxL1Jet]      ; ///< Jet threshold parameter:A,B,C
  UInt_t    fL1JetThreshold                     [maxL1Jet]      ; ///< L1 Jet Threshold
  UShort_t  fL1JetPatchIndex[maxL1JetPatchIndex][maxL1Jet]      ; ///< L1 Jet Patch Index
  
  UInt_t    fADC[maxnTRU][maxnmoduleInTRU]; ///< ADCs
  UInt_t    fPHOSSubregion[36];             ///< PHOS subregions
  
  UInt_t    fV0A        ;      ///< V0A
  UInt_t    fV0C        ;      ///< V0C
  UInt_t    fS[4]       ;      ///< PHOS Scale parameter
  UInt_t    fRho        ;      ///< Background Rho
  UInt_t    fPatchSize  ;      ///< Jet patch size
  
  UInt_t    fRegionEnable   ;  ///< Region Enable
  UInt_t    fFrameReceived  ;  ///< Frame Received
  UInt_t    fFwVersion      ;  ///< Firmware Version
  
  Int_t     GetnMod() const                     { 
    return  (fDetector == kEMCAL)? nModEMCAL :
    (fDetector == kDCAL )? nModDCAL  :
    0 ;                                         }
  
  virtual void    DecodeL1JetPatchIndexes(  const int i, UInt_t *word32, const int offset);
  virtual void    DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset);
  virtual void    DecodeL0GammaPatchIndexes(             UInt_t *word32, const int offset);
  virtual void    DecodeTRUADC(                          UInt_t *word32, const int offset);
  virtual void    DecodePHOSSubregion(                   UInt_t *word32, const int offset);
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerSTURawStream,2) ;
  /// \endcond
  
};

#endif //ALIEMCALTRIGGERSTURAWSTREAM_H
