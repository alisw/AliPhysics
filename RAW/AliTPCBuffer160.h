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

#ifndef AliTPCBUFFER160_H
#define AliTPCBUFFER160_H

#include <TObject.h>
#ifdef __CINT__
class fstream;
#else
#include "Riostream.h"
#endif


class AliTPCBuffer160:public TObject{
public:
  AliTPCBuffer160(){}//default constructor
  AliTPCBuffer160(const char* fileName,Int_t flag);//constructor
  AliTPCBuffer160(fstream* file, Int_t size);//constructor for reading a file
  virtual ~AliTPCBuffer160();//destructor
  AliTPCBuffer160(const AliTPCBuffer160 &source); // copy constructor
  AliTPCBuffer160& operator=(const AliTPCBuffer160 &source); // ass. op.
  void  FillBuffer(Int_t Val);
  //this method store a word into the buffer
  Int_t GetFreeCellNumber()const{return fFreeCellBuffer;}
  //this method return the number of free cells of the internal buffer
  Int_t GetNextBackWord();
  //this method return the next word of 10 bit (reading the file backward) if it exists -1 otherwise
  Int_t GetNext();
  //this method return the next word of 10 bit (reading the file forward) if it exists -1 otherwise
  void  WriteTrailer(Int_t WordsNumber,Int_t PadNumber,Int_t RowNumber,Int_t SecNumber);
  //this method is used to write the trailer
  void  ReadTrailer(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber);
  //this method is used to read the trailer when the file is read forward
  Int_t ReadTrailerBackward(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber);
  //this method is used to read the trailer when the file is read backward
  void  WriteDataHeader(Bool_t dummy, Bool_t comressed);
  //this method is used to write the data header
  void  SetVerbose(Int_t val){fVerbose=val;}
  //this method is used to set the verbose level 
  //level  0 no output messages
  //level !=0 some messages are displayed during the run
  void  Flush();
  //this method is used to fill the buffer with 2AA hexadecimal value and save it into the output file
  Int_t GetFillWordsNum()const{return fEndingFillWords;}
private:
  void  PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit);
  //this method is used to pack bits into a word of 32 bits
  void  UnpackWord(UInt_t PackedWord, Int_t StartBit, Int_t StopBit, UInt_t &Word);
  //this method is used to read a precise number of bits from a word of 32 bits
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
  fstream* f;           //logical name of the I/O file
  Bool_t fCreated;      //true if f was created by the buffer
  Int_t fMaskBackward;  //bit mask for backward reading of a file
  UInt_t fFilePosition;//'pointer' to the actual position in the file
  UInt_t fFileEnd;     //position of the last element of the file (File dimension)
  UInt_t fDataHeaderPos;//Data header position
  Int_t  fEndingFillWords;//Few words at the end of the stream
  ClassDef(AliTPCBuffer160,1)
};

#endif
