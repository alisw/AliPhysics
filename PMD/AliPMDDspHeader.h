#ifndef ALIPMDDSPHEADER_H
#define ALIPMDDSPHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Author - Basanta K. Nandi

#include <TObject.h>

class AliPMDDspHeader : public TObject {

public:
   AliPMDDspHeader();
   AliPMDDspHeader(const AliPMDDspHeader &dsph);
   AliPMDDspHeader& operator=(const AliPMDDspHeader &dsph);

   virtual ~AliPMDDspHeader();

   // dsp header
   void  SetDataKey(Int_t dkey)            {fDataKey = dkey;}
   void  SetTotalLength(Int_t totlength)   {fTotalLength = totlength;}
   void  SetRawDataLength(Int_t rawlength) {fRawDataLength = rawlength;}
   void  SetDspId(Int_t dspid)             {fDspId = dspid;}
   void  SetBlkL1ATrigger(Int_t trword1)   {fBlkL1ATrigger = trword1;}
   void  SetMiniEventId(Int_t trword2)     {fMiniEventId = trword2;}
   void  SetL1ATrigger(Int_t trword3)      {fL1ATrigger = trword3;}
   void  SetL1RTrigger(Int_t trword4)      {fL1RTrigger = trword4;}
   void  SetPaddingWord(UInt_t padword)    {fPaddingWord = padword;}
   void  SetErrorWord(Int_t errw)          {fErrorWord = errw;}

   void  SetHeader(Int_t *header);

   Int_t  GetHeaderLength()       const {return fgkHeaderLength;}
   UInt_t GetDefaultPaddingWord() const {return fgkDefaultPaddingWord;}

   Int_t GetDataKey()       const {return fDataKey;}
   Int_t GetTotalLength()   const {return fTotalLength;}
   Int_t GetRawDataLength() const {return fRawDataLength;}
   Int_t GetDspId()         const {return fDspId;}
   Int_t GetBlkL1Trigger()  const {return fBlkL1ATrigger;}  
   Int_t GetMiniEventId()   const {return fMiniEventId;}  
   Int_t GetL1ATrigger()    const {return fL1ATrigger;}  
   Int_t GetL1RTrigger()    const {return fL1RTrigger;}  
   Int_t GetPaddingWord()   const {return fPaddingWord;}  
   Int_t GetErrorWord()     const {return fErrorWord;}  

 private:

   Int_t     fDataKey;        // Data key word for FRT header
   Int_t     fTotalLength;    // total length of block structure
   Int_t     fRawDataLength;  // length of raw data
   Int_t     fDspId;          // Dsp id
   Int_t     fBlkL1ATrigger;  // 1st trigger word
   Int_t     fMiniEventId;    // 1st trigger word
   Int_t     fL1ATrigger;     // 1st trigger word
   Int_t     fL1RTrigger;     // 1st trigger word
   UInt_t    fPaddingWord;    // padding word (nb words odd:1, even:0)
   Int_t     fErrorWord;      // Error word (nb words odd:1, even:0)

   static const Int_t  fgkHeaderLength;       // header length in word
   static const UInt_t fgkDefaultPaddingWord; // Default padding word

   ClassDef(AliPMDDspHeader,1)  // PMD dsp Header
};
#endif
