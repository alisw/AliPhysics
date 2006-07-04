#ifndef ALIPMDBLOCKHEADER_H
#define ALIPMDBLOCKHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>

class AliPMDBlockHeader : public TObject {

public:
   AliPMDBlockHeader();
   AliPMDBlockHeader(const AliPMDBlockHeader &blockh);
   AliPMDBlockHeader& operator=(const AliPMDBlockHeader &blockh);

   virtual ~AliPMDBlockHeader();

   // Block header

   void  SetTotalLength(Int_t totlength)   {fTotalLength = totlength;}
   void  SetRawDataLength(Int_t rawlength) {fRawDataLength = rawlength;}
   void  SetDspId(Int_t dspid)             {fDspId = dspid;}
   void  SetTriggerWord1(Int_t trword1)    {fTrWord1 = trword1;}
   void  SetTriggerWord2(Int_t trword2)    {fTrWord2 = trword2;}
   void  SetTriggerWord3(Int_t trword3)    {fTrWord3 = trword3;}
   void  SetTriggerWord4(Int_t trword4)    {fTrWord4 = trword4;}
   void  SetPaddingWord(Int_t padword)     {fPadWord = padword;}
   void  SetHeader(Int_t *header);


   Int_t GetHeaderLength()  const {return fgkHeaderLength;}
   Int_t GetTotalLength()   const {return fTotalLength;}
   Int_t GetRawDataLength() const {return fRawDataLength;}
   Int_t GetDspId()         const {return fDspId;}
   Int_t GetTriggerWord1()  const {return fTrWord1;}  
   Int_t GetTriggerWord2()  const {return fTrWord2;}  
   Int_t GetTriggerWord3()  const {return fTrWord3;}  
   Int_t GetTriggerWord4()  const {return fTrWord4;}  
   Int_t GetPaddingWord()   const {return fPadWord;}  

 private:

   Int_t     fTotalLength;    // total length of block structure
   Int_t     fRawDataLength;  // length of raw data
   Int_t     fDspId;          // Dsp id
   Int_t     fTrWord1;        // 1st trigger word
   Int_t     fTrWord2;        // 1st trigger word
   Int_t     fTrWord3;        // 1st trigger word
   Int_t     fTrWord4;        // 1st trigger word
   Int_t     fPadWord;        // Padding word

   static const Int_t fgkHeaderLength;   // header length in word

   ClassDef(AliPMDBlockHeader,0)  // PMD Block Header
};
#endif
