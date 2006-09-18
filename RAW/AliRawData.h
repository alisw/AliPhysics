#ifndef ALIRAWDATA_H
#define ALIRAWDATA_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawData                                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif


class AliRawData : public TObject {

public:
   AliRawData();
   virtual ~AliRawData() { if (fOwner) delete [] fRawData; }

   inline void SetSize(Int_t size);
   inline void SetBuffer(void *buf, Int_t size);
   Int_t       GetSize() const { return fSize; }
   void       *GetBuffer() { return fRawData; }

private:
   Int_t   fSize;         // number of raw data bytes
   Int_t   fBufSize;      //!actual size of fRawData
   char   *fRawData;      //[fSize] raw event data
   Bool_t  fOwner;        //!if true object owns fRawData buffer

   AliRawData(const AliRawData &);      // not implemented, usage causes
   AliRawData &operator=(const AliRawData &);  // link time error

   ClassDef(AliRawData,1)  // Alice raw event buffer
};

void AliRawData::SetSize(Int_t size)
{
   if (size > fBufSize) {
      if (fOwner) delete [] fRawData;
      fRawData = new char [size];
      fBufSize = size;
      fOwner   = kTRUE;
   }
   fSize = size;
}

void AliRawData::SetBuffer(void *buf, Int_t size)
{
   if (fOwner) delete [] fRawData;
   fRawData = (char *) buf;
   fBufSize = size;
   fSize    = size;
   fOwner   = kFALSE;
}

#endif
