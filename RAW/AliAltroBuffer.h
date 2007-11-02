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

#ifndef AliALTROBUFFER_H
#define AliALTROBUFFER_H

#include <TObject.h>

class AliFstream;

class AliAltroMapping;

class AliAltroBuffer: public TObject {
 public:
  AliAltroBuffer(const char* fileName, AliAltroMapping *mapping = NULL);
  virtual ~AliAltroBuffer();

  void  FillBuffer(Int_t val);
  //this method stores a word into the buffer
  Int_t GetFreeCellNumber()const {return fFreeCellBuffer;}
  //this method returns the number of free cells of the internal buffer

  void  WriteTrailer(Int_t wordsNumber, Int_t padNumber, 
		     Int_t rowNumber, Int_t secNumber);
  //this method is used to write the trailer
  void  WriteTrailer(Int_t wordsNumber, Short_t hwAddress); 
  //this method is used to write the trailer
  void  WriteDummyTrailer(Int_t wordsNumber, Int_t padNumber, 
			  Int_t rowNumber, Int_t secNumber);
  //this method is used to write a dummy trailer

  void  WriteChannel(Int_t padNumber, Int_t rowNumber, Int_t secNumber,
		     Int_t nTimeBins, const Int_t* adcValues, 
		     Int_t threshold = 0);
  //this method is used to write all ADC values and the trailer of a channel
  void  WriteChannel(Short_t hwAddress,
		     Int_t nTimeBins, const Int_t* adcValues, 
		     Int_t threshold = 0);
  //this method is used to write all ADC values and the trailer of a channel
  Int_t WriteBunch(Int_t nTimeBins, const Int_t* adcValues,
		   Int_t threshold = 0);
  //this method is used to write all ADC values

  void  WriteDataHeader(Bool_t dummy, Bool_t compressed);
  //this method is used to write the data header

  void  WriteRCUTrailer(Int_t rcuId);
  //this method is used to write the RCU trailer

  void  SetVerbose(Int_t val) {fVerbose = val;}
  //this method is used to set the verbose level 
  //level  0 no output messages
  //level !=0 some messages are displayed during the run
  void  Flush();
  //this method is used to fill the buffer with 2AA hexadecimal value and save it into the output file

  void  SetMapping(AliAltroMapping *mapping) { fMapping = mapping; }

 protected:
  AliAltroBuffer(const AliAltroBuffer& source);
  AliAltroBuffer& operator = (const AliAltroBuffer& source);

  UInt_t fBuffer[5];    //Buffer dimension is 32*5=160 bits and it contains 16 values
                        //A value is never splitted in two Buffer


  Int_t fShift;         //This variable contains the number of free bits in the current cell of
                        //the Buffer after that the value Val is been inserted.
                        //size of Int_t is 32 bit that is the same size of a cell of Buffer so 
                        //the shift operation are performed only on value Val.
  Int_t fCurrentCell;   //This variable contains the cell number of the cell currently used 
  Int_t fFreeCellBuffer;//number of free cells of the buffer
  Int_t fVerbose;       //verbose level
  AliFstream* fFile;    //logical name of the I/O file
  UInt_t fDataHeaderPos;//Data header position

  // Now the parameters for the mapping
  AliAltroMapping*    fMapping;      // Pointer to the mapping handler

  ClassDef(AliAltroBuffer,0)  // Interface to the Altro format
};

#endif
