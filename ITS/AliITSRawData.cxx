////////////////////////////////////////////////
//  RawData classes for set:ITS               //
////////////////////////////////////////////////


#include "AliITSRawData.h"

ClassImp(AliITSRawData)

ClassImp(AliITSInStream)

//_____________________________________________________________________________

AliITSInStream::AliITSInStream():
fStreamLen(0),
fInStream(0){
  //default constructor
}
//_____________________________________________________________________________

AliITSInStream::AliITSInStream(UInt_t length):
fStreamLen(length),
fInStream(0){
  //
  // Creates a stream of unsigned chars
  //
  
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
AliITSInStream::AliITSInStream(const AliITSInStream &source) : TObject(source),
fStreamLen(source.fStreamLen),
fInStream(source.fInStream){
  //     Copy Constructor 

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
Bool_t AliITSInStream::CheckCount(UInt_t count) {
  //check boundaries
  if (count <= (UInt_t)fStreamLen) return kTRUE;
  else {
    Error("CheckCount", "actual size is %d, the necessary size is %d",fStreamLen,count);
    return kFALSE;
  }
}


ClassImp(AliITSOutStream)
  
  //_______________________________________________________________________
  
AliITSOutStream::AliITSOutStream():
fStreamLen(0),
fOutStream(0){
  //default constructor

}

//__________________________________________________________________________

AliITSOutStream::AliITSOutStream(UInt_t length):
fStreamLen(length),
fOutStream(0){
  //
  // Creates a stream of unsigned chars
  //
  
  fOutStream = new UInt_t[length];  
  ClearStream(); 
  
}

//_____________________________________________________________________________
AliITSOutStream::~AliITSOutStream()
{
  //destructor
  if (fOutStream) delete[] fOutStream;
}

//__________________________________________________________________________
AliITSOutStream::AliITSOutStream(const AliITSOutStream &source):TObject(source),
fStreamLen(source.fStreamLen),
fOutStream(source.fOutStream){
  //     Copy Constructor 

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
  memset(fOutStream,0,sizeof(UInt_t)*fStreamLen);
}

//_____________________________________________________________________________
Bool_t AliITSOutStream::CheckCount(UInt_t count)
{
  //check boundaries
  if (count < fStreamLen) return kTRUE;
  else {
    Error("CheckCount", "actual size is %d, the necessary size is %d",fStreamLen,count);
    return kFALSE;
  }
}
