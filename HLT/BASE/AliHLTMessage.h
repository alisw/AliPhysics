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
// TMessage                                                             //
//                                                                      //
// Message buffer class used for serializing objects and sending them   //
// over the network.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

// TBuffer has been made pure virtual in root version v5-15-02, this
// requires to inherit from TBufferFile instead of TBuffer.
// TMessage is not really used by this class but by including it we also get
// TBufferFile if this exists. The define ROOT_TBufferFile can than be used
// to differentiate between the usage of TBuffer or TBufferFile.
#include "TMessage.h"

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
 * - the <tt> AliHLTMessage(void *buf, Int_t bufsize)</tt> constructor has been
 *   made public in order to be used externally.
 *
 * The class can be used to extract an object from an input data block, or a
 * data block received via the HOMER interface or from the file writer.
 * <pre>
 *  AliHLTMessage msg(buffer, size);
 *  TObject* pObj=msg.ReadObject(msg.GetClass());
 * </pre>
 *
 * A simple test macro for a file can look like
 * <pre>
 *  const char* filename="myobject.dat";
 *  TString param=filename;
 *  param+="?filetype=raw";
 *  TFile file(param);
 *  if (file.IsZombie()) {
 *    cout << "can not open file " << filename << endl;
 *    return;
 *  }
 *  
 *  TArrayC buffer(file.GetSize());
 *  TArrayC tgtbuffer(file.GetSize());
 *  if (file.ReadBuffer(buffer.GetArray(), buffer.GetSize())) {
 *    cout << "error reading file " << filename << endl;
 *    return;
 *  }
 *
 *  AliHLTMessage msg(buffer.GetArray(), buffer.GetSize());
 *  TObject* pObj=msg.ReadObject(msg.GetClass());
 * </pre>
 *
 * @see AliHLTRootFileWriterComponent for an easy way to save objects
 * exported via AliHLTMessage in a ROOT file.
 *
 * To serialize an object into a buffer, the normal ROOT TMessage mechanism
 * can be used.
 * <pre>
 *    AliHLTMessage msg(kMESS_OBJECT);
 *    msg.WriteObject(pObject);
 *    Int_t iMsgLength=msg.Length();
 *    if (iMsgLength>0) {
 *      msg.SetLength(); // sets the length to the first (reserved) word
 *      char* pMsgBuffer msg.Buffer();
 *      // do something with pMsgBuffer and iMsgLenghth
 *    }
 * </pre>
 */
class AliHLTMessage 
:
# if defined(ROOT_TBufferFile)
public TBufferFile,
#else
public TBuffer,
#endif
public AliHLTLogging {

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

   /** the minimum size of a serialized TObject */
   static const Int_t fgkMinimumSize; //!transient

   /** a default buffer describing an empty message */
   static UInt_t fgkDefaultBuffer[2]; //!transient

   ClassDef(AliHLTMessage,0)  // Message buffer class
};

#endif // ALIHLTMESSAGE_H
