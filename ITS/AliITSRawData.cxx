////////////////////////////////////////////////
//  RawData classes for set:ITS               //
////////////////////////////////////////////////


#include "AliITSRawData.h"

ClassImp(AliITSRawData)

ClassImp(AliITSInStream)

//_____________________________________________________________________________

AliITSInStream::AliITSInStream()
{
  //default constructor
  fStreamLen=0;
  fInStream=0;
}
//_____________________________________________________________________________

AliITSInStream::AliITSInStream(ULong_t length)
{
  //
  // Creates a stream of unsigned chars
  //
  
  fStreamLen = length;
  fInStream = new UChar_t[length]; 
  
  ClearStream(); 
  
}

//_____________________________________________________________________________
AliITSInStream::~AliITSInStream()
{
  //destructor
  if (fInStream) delete[] fInStream;
}

//__________________________________________________________________________
AliITSInStream::AliITSInStream(const AliITSInStream &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fStreamLen = source.fStreamLen;
  this->fInStream = source.fInStream;
  return;
}

//_________________________________________________________________________
AliITSInStream& 
  AliITSInStream::operator=(const AliITSInStream &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fStreamLen = source.fStreamLen;
  this->fInStream = source.fInStream;
  return *this;
}

//_____________________________________________________________________________
void AliITSInStream::ClearStream()
{
  //clear the array
  memset(fInStream,0,sizeof(UChar_t)*fStreamLen);
}


//_____________________________________________________________________________
Bool_t AliITSInStream::CheckCount(ULong_t count) {
  //check boundaries
  if (count <= (ULong_t)fStreamLen) return kTRUE;
  else {
    Error("CheckCount", "actual size is %d, the necessary size is %d",fStreamLen,count);
    return kFALSE;
  }
}

//____________________________________________________________________________
void AliITSInStream::Streamer(TBuffer &R__b){
  // Stream an object of class AliITSInStream.
  
  static unsigned char *array;
  static Bool_t make=kTRUE;
  
  if (R__b.IsReading()) {
	  R__b >> fStreamLen;
	  //printf("Streamer: fStreamLen %d\n",fStreamLen);
          if (make) array=new unsigned char[fStreamLen];
          make=kFALSE;
          memset(array,0,sizeof(UChar_t)*fStreamLen);
          fInStream=array;
	  R__b.ReadFastArray(fInStream,fStreamLen);
	  
  } else {
    R__b << fStreamLen;
    R__b.WriteFastArray(fInStream,fStreamLen);
  }
}



ClassImp(AliITSOutStream)
  
  //_______________________________________________________________________
  
  AliITSOutStream::AliITSOutStream() {
  //default constructor
  fStreamLen=0;
  fOutStream=0;
}

//__________________________________________________________________________

AliITSOutStream::AliITSOutStream(ULong_t length) {
  //
  // Creates a stream of unsigned chars
  //
  
  fStreamLen = length;
  fOutStream = new ULong_t[length];  
  ClearStream(); 
  
}

//_____________________________________________________________________________
AliITSOutStream::~AliITSOutStream()
{
  //destructor
  if (fOutStream) delete[] fOutStream;
}

//__________________________________________________________________________
AliITSOutStream::AliITSOutStream(const AliITSOutStream &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fStreamLen = source.fStreamLen;
  this->fOutStream = source.fOutStream;
  return;
}

//_________________________________________________________________________
AliITSOutStream& 
  AliITSOutStream::operator=(const AliITSOutStream &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fStreamLen = source.fStreamLen;
  this->fOutStream = source.fOutStream;
  return *this;
}

//_____________________________________________________________________________
void AliITSOutStream::ClearStream()
{
  // clear stream
  memset(fOutStream,0,sizeof(ULong_t)*fStreamLen);
}

//_____________________________________________________________________________
Bool_t AliITSOutStream::CheckCount(ULong_t count)
{
  //check boundaries
  if (count < fStreamLen) return kTRUE;
  else {
    Error("CheckCount", "actual size is %d, the necessary size is %d",fStreamLen,count);
    return kFALSE;
  }
}

//____________________________________________________________________________
void AliITSOutStream::Streamer(TBuffer &R__b){
  
  // Stream an object of class AliITSOutStream.
  
  static unsigned long *array;
  static Bool_t make=kTRUE;
  
  if (R__b.IsReading()) {
	  R__b >> fStreamLen;
	  //printf("Streamer: fStreamLen %d\n",fStreamLen);
          if (make) array=new unsigned long[fStreamLen];
          make=kFALSE;
          memset(array,0,sizeof(ULong_t)*fStreamLen);
          fOutStream=array;
	  R__b.ReadFastArray(fOutStream,fStreamLen);
	  
  } else {
    R__b << fStreamLen;
    R__b.WriteFastArray(fOutStream,fStreamLen);
  }
}

