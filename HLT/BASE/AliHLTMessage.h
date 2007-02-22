// $Id$

/** @file   AliHLTMessage.h
    @author Matthias Richter (customization of Root TMessage )
    @date   
    @brief  Serialization of Root objects in the ALICE HLT. */

// This is the original Root TMessage implementation with a few minor
// modifications, original revision:
// root/net: v5-14-00 $: TMessage.h,v 1.9 2005/12/09 15:12:19 rdm
// Author: Fons Rademakers   19/12/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ALIHLTMESSAGE_H
#define ALIHLTMESSAGE_H


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMessage                                                        //
//                                                                      //
// Message buffer class used for serializing objects and sending them   //
// over the network.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TBuffer
#include "TBuffer.h"
#endif
#ifndef ROOT_MessageTypes
#include "MessageTypes.h"
#endif

#include "AliHLTLogging.h"
/**
 * @class AliHLTMessage
 * Serialization of Root objects for transport in the Alice HLT analysis
 * chain.
 * This is the original Root TMessage implementation with a few minor
 * modifications.
 * - the AliHLTMessage(void *buf, Int_t bufsize) constructor has been made
 *   public in order to be used externally.
 */
class AliHLTMessage : public TBuffer, public AliHLTLogging {

public:
   AliHLTMessage(UInt_t what = kMESS_ANY);
   AliHLTMessage(void *buf, Int_t bufsize);
   virtual ~AliHLTMessage();

   void SetLength() const;

   void     Forward();
   TClass  *GetClass() const { return fClass; }
   void     Reset();
   void     Reset(UInt_t what) { SetWhat(what); Reset(); }
   UInt_t   What() const { return fWhat; }
   void     SetWhat(UInt_t what);

   void     SetCompressionLevel(Int_t level = 1);
   Int_t    GetCompressionLevel() const { return fCompress; }
   Int_t    Compress();
   Int_t    Uncompress();
   char    *CompBuffer() const { return fBufComp; }
   Int_t    CompLength() const { return (Int_t)(fBufCompCur - fBufComp); }

private:
   UInt_t   fWhat;        //Message type
   TClass  *fClass;       //If message is kMESS_OBJECT pointer to object's class
   Int_t    fCompress;    //Compression level from 0 (not compressed) to 9 (max compression)
   char    *fBufComp;     //Compressed buffer
   char    *fBufCompCur;  //Current position in compressed buffer
   char    *fCompPos;     //Position of fBufCur when message was compressed

   // AliHLTMessage objects cannot be copied or assigned
   AliHLTMessage(const AliHLTMessage &);           // not implemented
   void operator=(const AliHLTMessage &);     // not implemented

   ClassDef(AliHLTMessage,0)  // Message buffer class
};

#endif // ALIHLTMESSAGE_H
