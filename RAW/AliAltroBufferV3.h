/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
// Class used for read-write the ALTRO data format //
/////////////////////////////////////////////////////

/*This class is an interface between the altro format file and the 
  user, and can be used in write or read mode
  In the write mode a new altro file is created and filled using the method FillBuffer().
  The name of the file is specified as parameter in the constructor as well as the type mode.
  In the Read mode the specified file is open and the values can be read using the
  methods GetNext() and GetNextBackWord().
  The first method is used to read the file forward while the second is used to read backward 
*/

#ifndef AliALTROBUFFERV3_H
#define AliALTROBUFFERV3_H

#include "AliAltroBuffer.h"

class AliAltroBufferV3: public AliAltroBuffer {
 public:
  AliAltroBufferV3(const char* fileName, AliAltroMapping *mapping = NULL);
  virtual ~AliAltroBufferV3();

  virtual void  FillBuffer(Int_t val);
  //this method stores a word into the buffer

  virtual void  WriteTrailer(Int_t wordsNumber, Short_t hwAddress); 
  //this method is used to write the trailer

  virtual void  WriteRCUTrailer(Int_t rcuId);
  //this method is used to write the RCU trailer

  enum { kMaxWords = 1024 };

 protected:
  void          ReverseAndWrite();
  //this method reverse the altro data order and write the buffer to the file

  AliAltroBufferV3(const AliAltroBufferV3& source);
  AliAltroBufferV3& operator = (const AliAltroBufferV3& source);

  UShort_t fArray[kMaxWords]; // Temporary array needed in reverting data order
  Int_t    fN;                // Size of the temporary array

  ClassDef(AliAltroBufferV3,0)  // Interface to the Altro format
};

#endif
