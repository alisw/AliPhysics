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
#ifdef __CINT__
class fstream;
#else
#include "Riostream.h"
#endif

class AliAltroMapping;

class AliAltroBuffer: public TObject {
 public:
  AliAltroBuffer(const char* fileName, Int_t flag, const AliAltroMapping *mapping = NULL);
  virtual ~AliAltroBuffer();

  void  FillBuffer(Int_t val);
  //this method stores a word into the buffer
  Int_t GetFreeCellNumber()const {return fFreeCellBuffer;}
  //this method returns the number of free cells of the internal buffer
  Int_t GetNextBackWord();
  //this method returns the next word of 10 bit (reading the file backward) if it exists -1 otherwise
  Int_t GetNext();
  //this method returns the next word of 10 bit (reading the file forward) if it exists -1 otherwise

  void  WriteTrailer(Int_t wordsNumber, Int_t padNumber, 
		     Int_t rowNumber, Int_t secNumber);
  //this method is used to write the trailer
  void  WriteTrailer(Int_t wordsNumber, Short_t hwAdress); 
  //this method is used to write the trailer
  void  WriteDummyTrailer(Int_t wordsNumber, Int_t padNumber, 
			  Int_t rowNumber, Int_t secNumber);
  //this method is used to write a dummy trailer
  Bool_t ReadTrailer(Int_t& wordsNumber, Int_t& padNumber, 
		     Int_t& rowNumber, Int_t &secNumber);
  //this method is used to read the trailer when the file is read forward
  Bool_t ReadTrailer(Int_t& wordsNumber, Short_t& hwAdress); 
  //this method is used to read the trailer when the file is read forward
  Bool_t ReadDummyTrailer(Int_t& wordsNumber, Int_t& padNumber, 
			  Int_t& rowNumber, Int_t &secNumber);
  //this method is used to read the trailer when the file is read forward
  Bool_t ReadTrailerBackward(Int_t& wordsNumber, Int_t& padNumber, 
			     Int_t& rowNumber, Int_t& secNumber);
  //this method is used to read the trailer when the file is read backward
  Bool_t ReadTrailerBackward(Int_t& wordsNumber, Short_t& hwAdress); 
  //this method is used to read the trailer when the file is read backward
  Bool_t ReadDummyTrailerBackward(Int_t& wordsNumber, Int_t& padNumber, 
				  Int_t& rowNumber, Int_t& secNumber);
  //this method is used to read the trailer when the file is read backward

  void  WriteChannel(Int_t padNumber, Int_t rowNumber, Int_t secNumber,
		     Int_t nTimeBins, const Int_t* adcValues, 
		     Int_t threshold = 0);
  //this method is used to write all ADC values and the trailer of a channel
  void  ReadChannel(Int_t padNumber, Int_t rowNumber,  Int_t secNumber,
		    Int_t& nTimeBins, Int_t* adcValues);
  //this method is used to read all ADC values and the trailer of a channel

  void  WriteDataHeader(Bool_t dummy, Bool_t compressed);
  //this method is used to write the data header
  Bool_t ReadDataHeader();
  //this method is used to read the data header
  void  SetVerbose(Int_t val) {fVerbose = val;}
  //this method is used to set the verbose level 
  //level  0 no output messages
  //level !=0 some messages are displayed during the run
  void  Flush();
  //this method is used to fill the buffer with 2AA hexadecimal value and save it into the output file
  Int_t GetFillWordsNum() const {return fEndingFillWords;}

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
  Int_t fFlag;          //0 read  1 write
  Int_t fVerbose;       //verbose level
  fstream* fFile;       //logical name of the I/O file
  Int_t fMaskBackward;  //bit mask for backward reading of a file
  UInt_t fFilePosition;//'pointer' to the actual position in the file
  UInt_t fFileEnd;     //position of the last element of the file (File dimension)
  UInt_t fDataHeaderPos;//Data header position
  Int_t  fEndingFillWords;//Few words at the end of the stream

  // Now the parameters for the mapping
  const AliAltroMapping*    fMapping;      // Pointer to the mapping handler

  ClassDef(AliAltroBuffer,0)  // Interface to the Altro format
};

#endif
