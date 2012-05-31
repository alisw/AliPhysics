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

   void  SetDataKey(Int_t dkey)            {fDataKey = dkey;}
   void  SetTotalLength(Int_t totlength)   {fTotalLength = totlength;}
   void  SetRawDataLength(Int_t rawlength) {fRawDataLength = rawlength;}
   void  SetDspId(Int_t dspid)             {fDspId = dspid;}
   void  SetL0Trigger(Int_t trword1)       {fL0Trigger = trword1;}
   void  SetMiniEventId(Int_t trword2)     {fMiniEventId = trword2;}
   void  SetEventId1(Int_t trword3)        {fEventId1 = trword3;}
   void  SetEventId2(Int_t trword4)        {fEventId2 = trword4;}

   void  SetHeader(Int_t *header);


   Int_t GetHeaderLength()  const {return fgkHeaderLength;}
   Int_t GetDataKey()       const {return fDataKey;}
   Int_t GetTotalLength()   const {return fTotalLength;}
   Int_t GetRawDataLength() const {return fRawDataLength;}
   Int_t GetDspId()         const {return fDspId;}
   Int_t GetL0Trigger()     const {return fL0Trigger;}  
   Int_t GetMiniEventId()   const {return fMiniEventId;}  
   Int_t GetEventId1()      const {return fEventId1;}  
   Int_t GetEventId2()      const {return fEventId2;}  


 private:

   Int_t     fDataKey;        // Data key word for CRT header
   Int_t     fTotalLength;    // total length of block structure
   Int_t     fRawDataLength;  // length of raw data
   Int_t     fDspId;          // Dsp id
   Int_t     fL0Trigger;      // L0 trigger word
   Int_t     fMiniEventId;    // Bunch crossing for mini-event id
   Int_t     fEventId1;       // Event Id in bunch crossing
   Int_t     fEventId2;       // Event Id in orbit number

   static const Int_t fgkHeaderLength;   // header length in word

   ClassDef(AliPMDBlockHeader,1)  // PMD Block Header
};
#endif
