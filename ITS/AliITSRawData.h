#ifndef ALIITSRAWDATA_H
#define ALIITSRAWDATA_H

////////////////////////////////////////////////
//  RawData classes for set:ITS               //
////////////////////////////////////////////////

#include <TObject.h>



class AliITSRawData: public TObject  {

  // most probably it should have a class AliITSHeaderEvent as data member
   
 public:
    AliITSRawData() {
      // constructor
    }
    virtual   ~AliITSRawData() {
      // destructor
    }

    ClassDef(AliITSRawData,1)     //RawData object for set:ITS

};

//___________________________________________
class AliITSInStream: public TObject{

 public:
  AliITSInStream();
  AliITSInStream(ULong_t length);
  virtual   ~AliITSInStream();
  AliITSInStream(const AliITSInStream &source); // copy constructor
  AliITSInStream& operator=(const AliITSInStream &source); // ass. operator
    
  void ClearStream();
  Bool_t CheckCount(ULong_t count);
  ULong_t  StreamLength() {
    // stream length
    return fStreamLen;
  }
  UChar_t *Stream() {
    // stream
    return fInStream;
  }
  
protected:
  
  // input stream of unsigned chars
  
  ULong_t     fStreamLen;       // Length of the array
  UChar_t    *fInStream;        // Pointer to an array of input unsigned chararacters
  
  
  
  ClassDef(AliITSInStream,1)     //Input Stream  object for set:ITS
    };

//___________________________________________
class AliITSOutStream: public TObject{
  
public:
  AliITSOutStream();
  
  AliITSOutStream(ULong_t length);
  virtual   ~AliITSOutStream();
  AliITSOutStream(const AliITSOutStream &source); // copy constructor
  AliITSOutStream& operator=(const AliITSOutStream &source); // assignment operator
    
  void ClearStream();
  Bool_t CheckCount(ULong_t count);
  ULong_t  StreamLength() {
    // stream length
    return fStreamLen;
  }
  ULong_t *Stream() {
    // stream
    return fOutStream;
  }
  
protected:
  
  // output stream of unsigned chars
  
  ULong_t     fStreamLen;        // Length of the array
  ULong_t    *fOutStream;        // Pointer to an array of unsigned long
  
  

    ClassDef(AliITSOutStream,1)     //Output Stream  object for set:ITS
};



#endif
