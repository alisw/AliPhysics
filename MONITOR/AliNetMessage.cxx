#include <TVirtualStreamerInfo.h>
#include <Bytes.h>
#include <TFile.h>
#include <TClass.h>

#include "AliNetMessage.h"

Bool_t AliNetMessage::fgEvolution = kFALSE;

ClassImp(AliNetMessage)

//______________________________________________________________________________
AliNetMessage::AliNetMessage(UInt_t what) 
  :
  TBufferFile(kWrite),
  fWhat(what),
  fClass(0),
  fBufUncompressed(0), 
  fInfos(NULL), 
  fEvolution(kFALSE)
{
   // Create a AliNetMessage object for storing objects. The "what" integer
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


//______________________________________________________________________________
AliNetMessage::AliNetMessage(void *buf, Int_t bufsize)
  :
  TBufferFile(kRead, bufsize, buf),
  fWhat(0),
  fClass(0),
  fBufUncompressed(0), 
  fInfos(NULL), 
  fEvolution(kFALSE)
{
   // Create a AliNetMessage object for reading objects. The objects will be
   // read from buf. Use the What() method to get the message type.

   // skip space at the beginning of the message reserved for the message length
   fBufCur += sizeof(UInt_t);

   *this >> fWhat;

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
AliNetMessage::~AliNetMessage()
{
   // Clean up
  Reset();
}

//______________________________________________________________________________
void AliNetMessage::EnableSchemaEvolutionForAll(Bool_t enable)
{
   // Static function enabling or disabling the automatic schema evolution.
   // By default schema evolution support is off.

   fgEvolution = enable;
}

//______________________________________________________________________________
Bool_t AliNetMessage::UsesSchemaEvolutionForAll()
{
   // Static function returning status of global schema evolution.

   return fgEvolution;
}

//______________________________________________________________________________
void AliNetMessage::ForceWriteInfo(TVirtualStreamerInfo *info, Bool_t /* force */)
{
   // Force writing the TStreamerInfo to the message.

   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TList();
				fInfos->Add(info);
   }
}

//______________________________________________________________________________
void AliNetMessage::Forward()
{
   // Change a buffer that was received into one that can be send, i.e.
   // forward a just received message.

   if (IsReading()) {
      SetWriteMode();
      SetBufferOffset(fBufSize);
      SetBit(kCannotHandleMemberWiseStreaming);
   }
}

//______________________________________________________________________________
void AliNetMessage::TagStreamerInfo(TVirtualStreamerInfo *info)
{
   // Remember that the StreamerInfo is being used in writing.

   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TList();
      fInfos->Add(info);
   }
}

//______________________________________________________________________________
void AliNetMessage::IncrementLevel(TVirtualStreamerInfo *info)
{
   // Increment level.

   TBufferFile::IncrementLevel(info);

   if (!info) return;
   if (fgEvolution || fEvolution) {
      if (!fInfos) fInfos = new TList();

      // add the streamer info, but only once
      // this assumes that there is only one version
      if (fInfos->FindObject(info->GetName())==NULL) {
				fInfos->Add(info);
      }
   }
}

//______________________________________________________________________________
void AliNetMessage::Reset()
{
   // Reset the message buffer so we can use (i.e. fill) it again.

   SetBufferOffset(sizeof(UInt_t) + sizeof(fWhat));
   ResetMap();

   if (fBufUncompressed) {
     delete [] fBufUncompressed;
     fBufUncompressed=NULL;
   }
}

//______________________________________________________________________________
void AliNetMessage::SetLength() const
{
   // Set the message length at the beginning of the message buffer.

   if (IsWriting()) {
      char *buf = Buffer();
      *((UInt_t*)buf) = (UInt_t)(Length() - sizeof(UInt_t));
   }
}

//______________________________________________________________________________
void AliNetMessage::SetWhat(UInt_t what)
{
   // Using this method one can change the message type a-posteriory.
   // In case you OR "what" with kMESS_ACK, the message will wait for
   // an acknowledgement from the remote side. This makes the sending
   // process synchronous.

   fWhat = what;

   char *buf = Buffer();
   buf += sizeof(UInt_t);   // skip reserved length space
   tobuf(buf, what);
}

//______________________________________________________________________________
void AliNetMessage::WriteObject(const TObject *obj)
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
         fInfos = new TList();
   }

   WriteObjectAny(obj, TObject::Class());
}

