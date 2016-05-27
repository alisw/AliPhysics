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
#include "TVirtualStreamerInfo.h"
#include "Bytes.h"
#include "TFile.h"
#include "TProcessID.h"
#include "TClass.h"

extern "C" void R__zip (Int_t cxlevel, Int_t *nin, char *bufin, Int_t *lout, char *bufout, Int_t *nout);
extern "C" void R__unzip(Int_t *nin, UChar_t *bufin, Int_t *lout, char *bufout, Int_t *nout);
const Int_t kMAXBUF = 0xffffff;

Bool_t AliHLTMessage::fgEvolution = kFALSE;

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
  , fBufUncompressed(0)
  , fBitsPIDs(0)
  , fInfos(NULL)
  , fEvolution(kFALSE)
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

   SetBit(kCannotHandleMemberWiseStreaming);
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
  , fBufUncompressed(0)
  , fBitsPIDs(0)
  , fInfos(NULL)
  , fEvolution(kFALSE)
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
      // Matthias Sep 2008
      // NOTE: this is not done in TMessage and will lead to the deletion
      // of the buffer. This is not allowed in case of HLT where the
      // buffer is handled by the framework. In general, I think this
      // is a very bad idea to do it like that in TMessage
      fBufComp    = NULL;
      fBufCompCur = 0;
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
  Reset();
}

//______________________________________________________________________________
void AliHLTMessage::EnableSchemaEvolutionForAll(Bool_t enable)
{
   // Static function enabling or disabling the automatic schema evolution.
   // By default schema evolution support is off.

   fgEvolution = enable;
}

//______________________________________________________________________________
Bool_t AliHLTMessage::UsesSchemaEvolutionForAll()
{
   // Static function returning status of global schema evolution.

   return fgEvolution;
}

//______________________________________________________________________________
void AliHLTMessage::ForceWriteInfo(TVirtualStreamerInfo *info, Bool_t /* force */)
{
   // Force writing the TStreamerInfo to the message.

   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TObjArray();
      if (fInfos->FindObject(info->GetName())==NULL) {
	fInfos->Add(info);
      }
   }
}

//______________________________________________________________________________
void AliHLTMessage::Forward()
{
   // Change a buffer that was received into one that can be send, i.e.
   // forward a just received message.

   if (IsReading()) {
      SetWriteMode();
      SetBufferOffset(fBufSize);
      SetBit(kCannotHandleMemberWiseStreaming);

      if (fBufComp) {
         fCompPos = fBufCur;
      }
   }
}

//______________________________________________________________________________
void AliHLTMessage::TagStreamerInfo(TVirtualStreamerInfo *info)
{
   // Remember that the StreamerInfo is being used in writing.

   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TObjArray();
      fInfos->Add(info);
   }
}

//______________________________________________________________________________
void AliHLTMessage::IncrementLevel(TVirtualStreamerInfo *info)
{
   // Increment level.

   TBufferFile::IncrementLevel(info);

   if (!info) return;
   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TObjArray();

      // add the streamer info, but only once
      // this assumes that there is only one version
      if (fInfos->FindObject(info->GetName())==NULL) {
	fInfos->Add(info);
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
   if (fBufUncompressed) {
     delete [] fBufUncompressed;
     fBufUncompressed=NULL;
   }

   delete fInfos; fInfos=NULL;
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
	 *((UInt_t*)buf) = (UInt_t)(CompLength() - sizeof(UInt_t));
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
   char *bufcur1 = fBufComp + hdrlen;
   frombuf(bufcur1, &buflen);
   UChar_t *bufcur = (UChar_t*)bufcur1;
   fBuffer  = new char[buflen];
   fBufUncompressed = fBuffer;
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

//______________________________________________________________________________
void AliHLTMessage::WriteObject(const TObject *obj)
{
   // Write object to message buffer.
   // When support for schema evolution is enabled the list of TStreamerInfo
   // used to stream this object is kept in fInfos. This information is used
   // by TSocket::Send that sends this list through the socket. This list is in
   // turn used by TSocket::Recv to store the TStreamerInfo objects in the
   // relevant TClass in case the TClass does not know yet about a particular
   // class version. This feature is implemented to support clients and servers
   // with either different ROOT versions or different user classes versions.

   if (fgEvolution || fEvolution) {
      if (fInfos)
         fInfos->Clear();
      else
         fInfos = new TObjArray();
   }

   fBitsPIDs.ResetAllBits();
   WriteObjectAny(obj, TObject::Class());
}

//______________________________________________________________________________
UShort_t AliHLTMessage::WriteProcessID(TProcessID *pid)
{
   // Check if the ProcessID pid is already in the message.
   // If not, then:
   //   - mark bit 0 of fBitsPIDs to indicate that a ProcessID has been found
   //   - mark bit uid+1 where uid id the uid of the ProcessID

   if (fBitsPIDs.TestBitNumber(0)) return 0;
   if (!pid)
      pid = TProcessID::GetPID();
   if (!pid) return 0;
   fBitsPIDs.SetBitNumber(0);
   UInt_t uid = pid->GetUniqueID();
   fBitsPIDs.SetBitNumber(uid+1);
   return 1;
}

AliHLTMessage* AliHLTMessage::Stream(TObject* pSrc, Int_t compression, unsigned verbosity, bool enableSchema)
{
  /// Helper function to stream an object into an AliHLTMessage
  /// The returned instance must be cleaned by the caller
  ///
  /// Get the data and data size from the message:
  ///  first check
  ///    pMsg->CompLength();
  ///    pMsg->CompBuffer();
  ///  if that is NULL
  ///    pMsg->Length();
  ///    pMsg->Buffer();
  ///
  /// Note: accessing scheme will be change din the future to just have the two
  ///       latter ones.
  if (!pSrc) return NULL;

  AliHLTLogging log;
  AliHLTMessage* pMsg=new AliHLTMessage(kMESS_OBJECT);
  if (!pMsg) {
    log.LoggingVarargs(kHLTLogError, "AliHLTMessage", "Stream" , __FILE__ , __LINE__ , "memory allocation failed");
    return NULL;
  }

  pMsg->EnableSchemaEvolution(enableSchema);
  pMsg->SetCompressionLevel(compression);
  pMsg->WriteObject(pSrc);
  if (pMsg->Length()>0) {
    // Matthias Sep 2008
    // NOTE: AliHLTMessage does implement it's own SetLength method
    // which is not architecture independent. The original SetLength
    // stores the size always in network byte order.
    // I'm trying to remember the rational for that, might be that
    // it was just some lack of knowledge. Want to change this, but
    // has to be done carefully to be backward compatible.
    pMsg->SetLength(); // sets the length to the first (reserved) word

    // does nothing if the level is 0
    pMsg->Compress();

    if (pMsg->CompBuffer()) {
      pMsg->SetLength(); // set once more to have the byte order
      if (verbosity>0) log.LoggingVarargs(kHLTLogInfo, "AliHLTMessage", "Stream" , __FILE__ , __LINE__ , "object %p type %s streamed: size %d", pSrc, pSrc->GetName(), pMsg->CompLength());
    } else {
      if (verbosity>0) log.LoggingVarargs(kHLTLogInfo, "AliHLTMessage", "Stream" , __FILE__ , __LINE__ , "object %p type %s streamed: size %d", pSrc, pSrc->GetName(), pMsg->Length());
    }
  }
  return pMsg;
}

TObject* AliHLTMessage::Extract(const void* pBuffer, unsigned bufferSize, unsigned verbosity)
{
   /// Helper function to extract an object from a buffer.
   /// The returned object must be cleaned by the caller
  AliHLTLogging log;
  if (!pBuffer || bufferSize<sizeof(AliHLTUInt32_t)) {
    if (verbosity>0) log.LoggingVarargs(kHLTLogWarning, "AliHLTMessage", "Extract" , __FILE__ , __LINE__ , "invalid input buffer %p %d", pBuffer, bufferSize);
    return NULL;
  }

  AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)pBuffer);
  if (firstWord==bufferSize-sizeof(AliHLTUInt32_t) &&
      firstWord>=34 /*thats the minimum size of a streamed TObject*/) {
    AliHLTMessage msg((AliHLTUInt8_t*)pBuffer, bufferSize);
    TClass* objclass=msg.GetClass();
    TObject* pObject=msg.ReadObject(objclass);
    if (pObject && objclass) {
      if (verbosity>0) log.LoggingVarargs(kHLTLogInfo, "AliHLTMessage", "Extract" , __FILE__ , __LINE__ , "object %p type %s created", pObject, objclass->GetName());
      return pObject;
    } else {
      if (verbosity>0) log.LoggingVarargs(kHLTLogWarning, "AliHLTMessage", "Extract" , __FILE__ , __LINE__ , "failed to create object from buffer of size %d", bufferSize);
    }
  } else {
    if (verbosity>0) log.LoggingVarargs(kHLTLogWarning, "AliHLTMessage", "Extract" , __FILE__ , __LINE__ , "not a streamed TObject: block size %d, indicated %d", bufferSize, firstWord+sizeof(AliHLTUInt32_t));
  }
  return NULL;
}

TObject* AliHLTMessage::Extract(const char* filename, unsigned verbosity)
{
   /// Helper function to extract an object from a file containing the streamed object.
   /// The returned object must be cleaned by the caller
  if (!filename) return NULL;
  
  AliHLTLogging log;
  TString input=filename;
  input+="?filetype=raw";
  TFile* pFile=new TFile(input);
  if (!pFile) return NULL;
  TObject* pObject=NULL;
  if (!pFile->IsZombie()) {
    pFile->Seek(0);
    TArrayC buffer;
    buffer.Set(pFile->GetSize());
    if (pFile->ReadBuffer(buffer.GetArray(), buffer.GetSize())==0) {
      pObject=Extract(buffer.GetArray(), buffer.GetSize(), verbosity);
    } else {
      log.LoggingVarargs(kHLTLogError, "AliHLTMessage", "Extract" , __FILE__ , __LINE__ , "failed reading %d byte(s) from file %s", pFile->GetSize(), filename);
    }
  }

  delete pFile;
  return pObject;
}
