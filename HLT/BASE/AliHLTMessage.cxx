// $Id$

/** @file   AliHLTMessage.cxx
    @author Matthias Richter (customization of Root TMessage )
    @date   
    @brief  Serialization of Root objects in the ALICE HLT. */

// This is the original Root TMessage implementation with a few minor
// modifications, original revision:
// root/net: v5-14-00 $: TMessage.cxx,v 1.6 2004/05/07 09:51:58 brun
// Author: Fons Rademakers   19/12/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMessage                                                             //
//                                                                      //
// Message buffer class used for serializing objects and sending them   //
// over a network. This class inherits from TBuffer the basic I/O       //
// serializer.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliHLTMessage.h"
#include "Bytes.h"
#include "TFile.h"

extern "C" void R__zip (Int_t cxlevel, Int_t *nin, char *bufin, Int_t *lout, char *bufout, Int_t *nout);
extern "C" void R__unzip(Int_t *nin, UChar_t *bufin, Int_t *lout, char *bufout, Int_t *nout);
const Int_t kMAXBUF = 0xffffff;

ClassImp(AliHLTMessage)

//______________________________________________________________________________
AliHLTMessage::AliHLTMessage(UInt_t what) 
  :
# ifdef ROOT_TBufferFile
  TBufferFile(kWrite),
# else
  TBuffer(kWrite),
# endif
  AliHLTLogging(),
  fWhat(what),
  fClass(0),
  fCompress(0),
  fBufComp(0),
  fBufCompCur(0),
  fCompPos(0)
{
   // Create a AliHLTMessage object for storing objects. The "what" integer
   // describes the type of message. Predifined ROOT system message types
   // can be found in MessageTypes.h. Make sure your own message types are
   // unique from the ROOT defined message types (i.e. 0 - 10000 are
   // reserved by ROOT). In case you OR "what" with kMESS_ACK, the message
   // will wait for an acknowledgement from the remote side. This makes
   // the sending process synchronous. In case you OR "what" with kMESS_ZIP,
   // the message will be compressed in TSocket using the zip algorithm
   // (only if message is > 256 bytes).

   // space at the beginning of the message reserved for the message length
   UInt_t   reserved = 0;
   *this << reserved;

   *this << what;

}

const Int_t AliHLTMessage::fgkMinimumSize=30;
UInt_t AliHLTMessage::fgkDefaultBuffer[2]={0,0};

//______________________________________________________________________________
AliHLTMessage::AliHLTMessage(void *buf, Int_t bufsize)
  :
# if defined(ROOT_TBufferFile)
  TBufferFile(kRead, bufsize>fgkMinimumSize?bufsize:sizeof(fgkDefaultBuffer), bufsize>fgkMinimumSize?buf:&fgkDefaultBuffer, 0),
# else
  TBuffer(kRead, bufsize>fgkMinimumSize?bufsize:sizeof(fgkDefaultBuffer), bufsize>fgkMinimumSize?buf:&fgkDefaultBuffer, 0),
# endif
  AliHLTLogging(),
  fWhat(0),
  fClass(0),
  fCompress(0),
  fBufComp(0),
  fBufCompCur(0),
  fCompPos(0)
{
   // Create a AliHLTMessage object for reading objects. The objects will be
   // read from buf. Use the What() method to get the message type.

   // skip space at the beginning of the message reserved for the message length
   fBufCur += sizeof(UInt_t);

   *this >> fWhat;

   if (fWhat & kMESS_ZIP) {
      // if buffer has kMESS_ZIP set, move it to fBufComp and uncompress
      fBufComp    = fBuffer;
      fBufCompCur = fBuffer + bufsize;
      fBuffer     = 0;
      Uncompress();
   }

   if (fWhat == kMESS_OBJECT) {
      InitMap();
      fClass = ReadClass();     // get first the class stored in message
      SetBufferOffset(sizeof(UInt_t) + sizeof(fWhat));
      ResetMap();
   } else {
      fClass = 0;
   }
}

//______________________________________________________________________________
AliHLTMessage::~AliHLTMessage()
{
   // Clean up compression buffer.
   delete [] fBufComp;
}

//______________________________________________________________________________
void AliHLTMessage::Forward()
{
   // Change a buffer that was received into one that can be send, i.e.
   // forward a just received message.

   if (IsReading()) {
      SetWriteMode();
      SetBufferOffset(fBufSize);

      if (fBufComp) {
         fCompPos = fBufCur;
      }
   }
}

//______________________________________________________________________________
void AliHLTMessage::Reset()
{
   // Reset the message buffer so we can use (i.e. fill) it again.

   SetBufferOffset(sizeof(UInt_t) + sizeof(fWhat));
   ResetMap();

   if (fBufComp) {
      delete [] fBufComp;
      fBufComp    = 0;
      fBufCompCur = 0;
      fCompPos    = 0;
   }
}

//______________________________________________________________________________
void AliHLTMessage::SetLength() const
{
   // Set the message length at the beginning of the message buffer.

   if (IsWriting()) {
      char *buf = Buffer();
      *((UInt_t*)buf) = (UInt_t)(Length() - sizeof(UInt_t));

      if (fBufComp) {
         buf = fBufComp;
	 *((UInt_t*)buf) = (UInt_t)(Length() - sizeof(UInt_t));
      }
   }
}

//______________________________________________________________________________
void AliHLTMessage::SetWhat(UInt_t what)
{
   // Using this method one can change the message type a-posteriory.
   // In case you OR "what" with kMESS_ACK, the message will wait for
   // an acknowledgement from the remote side. This makes the sending
   // process synchronous.

   fWhat = what;

   char *buf = Buffer();
   buf += sizeof(UInt_t);   // skip reserved length space
   tobuf(buf, what);

   if (fBufComp) {
      buf = fBufComp;
      buf += sizeof(UInt_t);   // skip reserved length space
      tobuf(buf, what | kMESS_ZIP);
   }
}

//______________________________________________________________________________
void AliHLTMessage::SetCompressionLevel(Int_t level)
{
   // Set the message compression level. Can be between 0 and 9 with 0
   // being no compression and 9 maximum compression. In general the default
   // level of 1 is the best compromise between achieved compression and
   // cpu time. Compression will only happen when the message is > 256 bytes.

   if (level < 0) level = 0;
   if (level > 9) level = 9;

   if (level != fCompress && fBufComp) {
      delete [] fBufComp;
      fBufComp    = 0;
      fBufCompCur = 0;
      fCompPos    = 0;
   }
   fCompress = level;
}

//______________________________________________________________________________
Int_t AliHLTMessage::Compress()
{
   // Compress the message. The message will only be compressed if the
   // compression level > 0 and the if the message is > 256 bytes.
   // Returns -1 in case of error (when compression fails or
   // when the message increases in size in some pathological cases),
   // otherwise returns 0.

   if (fCompress == 0) {
      // no compression specified
      if (fBufComp) {
         delete [] fBufComp;
         fBufComp    = 0;
         fBufCompCur = 0;
         fCompPos    = 0;
      }
      return 0;
   }

   if (fBufComp && fCompPos == fBufCur) {
      // the message was already compressed
      return 0;
   }

   // remove any existing compressed buffer before compressing modified message
   if (fBufComp) {
      delete [] fBufComp;
      fBufComp    = 0;
      fBufCompCur = 0;
      fCompPos    = 0;
   }

   if (Length() <= (Int_t)(256 + 2*sizeof(UInt_t))) {
      // this message is too small to be compressed
      return 0;
   }

   Int_t hdrlen   = 2*sizeof(UInt_t);
   Int_t messlen  = Length() - hdrlen;
   Int_t nbuffers = messlen / kMAXBUF;
   Int_t chdrlen  = 3*sizeof(UInt_t);   // compressed buffer header length
   Int_t buflen   = TMath::Max(512, chdrlen + messlen + 9*nbuffers);
   fBufComp       = new char[buflen];
   char *messbuf  = Buffer() + hdrlen;
   char *bufcur   = fBufComp + chdrlen;
   Int_t noutot   = 0;
   Int_t nzip     = 0;
   Int_t nout, bufmax;
   for (Int_t i = 0; i <= nbuffers; i++) {
      if (i == nbuffers)
         bufmax = messlen - nzip;
      else
         bufmax = kMAXBUF;
      R__zip(fCompress, &bufmax, messbuf, &bufmax, bufcur, &nout);
      if (nout == 0 || nout >= messlen) {
         //this happens when the buffer cannot be compressed
         delete [] fBufComp;
         fBufComp    = 0;
         fBufCompCur = 0;
         fCompPos    = 0;
         return -1;
      }
      bufcur  += nout;
      noutot  += nout;
      messbuf += kMAXBUF;
      nzip    += kMAXBUF;
   }
   fBufCompCur = bufcur;
   fCompPos    = fBufCur;

   bufcur = fBufComp;
   tobuf(bufcur, (UInt_t)(CompLength() - sizeof(UInt_t)));
   Int_t what = fWhat | kMESS_ZIP;
   tobuf(bufcur, what);
   tobuf(bufcur, Length());    // original uncompressed buffer length

   return 0;
}

//______________________________________________________________________________
Int_t AliHLTMessage::Uncompress()
{
   // Uncompress the message. The message will only be uncompressed when
   // kMESS_ZIP is set. Returns -1 in case of error, 0 otherwise.

   if (!fBufComp || !(fWhat & kMESS_ZIP))
      return -1;

   Int_t buflen;
   Int_t hdrlen = 2*sizeof(UInt_t);
   UChar_t *bufcur = (UChar_t*)fBufComp + hdrlen;
   frombuf((char *&)bufcur, &buflen);
   fBuffer  = new char[buflen];
   fBufSize = buflen;
   fBufCur  = fBuffer + sizeof(UInt_t) + sizeof(fWhat);
   fBufMax  = fBuffer + fBufSize;
   char *messbuf = fBuffer + hdrlen;

   Int_t nin, nout, nbuf;
   Int_t noutot = 0;
   while (1) {
      nin  = 9 + ((Int_t)bufcur[3] | ((Int_t)bufcur[4] << 8) | ((Int_t)bufcur[5] << 16));
      nbuf = (Int_t)bufcur[6] | ((Int_t)bufcur[7] << 8) | ((Int_t)bufcur[8] << 16);
      R__unzip(&nin, bufcur, &nbuf, messbuf, &nout);
      if (!nout) break;
      noutot += nout;
      if (noutot >= buflen - hdrlen) break;
      bufcur  += nin;
      messbuf += nout;
   }

   fWhat &= ~kMESS_ZIP;
   fCompress = 1;

   return 0;
}

