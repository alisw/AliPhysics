#ifndef ALIDIGITS_H
#define ALIDIGITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment digits
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include   "AliArrayI.h"
#include   "AliArrayS.h"
#include   "AliSegmentID.h"
class AliH2F;

class AliDigits: public AliSegmentID{ 
public:
  AliDigits();
  ~AliDigits();
  inline Short_t GetDigitFast(Int_t row, Int_t column);  //return value at given row and collumn
  inline void  SetDigitFast(Short_t value,Int_t row, Int_t column);  //set value at given row and collumn
  inline Bool_t BoundsOK(const char *where, Int_t row, Int_t col) ;  //Check If Bound Ok
  Bool_t OutOfBoundsError(const char *where, Int_t row, Int_t column);
  virtual void Allocate(Int_t rows, Int_t columns);  //construct empty buffer fDigits with size rows x columns
  virtual Short_t GetDigit(Int_t row, Int_t column);
  virtual void ExpandBuffer();  //expand buffer to twodimensional array
  virtual void CompresBuffer(Int_t bufferType,Int_t threshold); //compres buffer according buffertype algorithm   
  virtual Bool_t First(); //adjust  first valid current digit
  virtual Bool_t Next();  //addjust next valid current digit
  void SetThreshold(Int_t th) {fThreshold = th;}
  Int_t  GetThreshold() {return fThreshold;}
  Int_t CurrentRow(){ return fCurrentRow;}  //return current row
  Int_t CurrentColumn(){ return fCurrentCol;} //return current column
  Int_t CurrentDigit() {return fElements->At(fCurrentIndex);} //return degit for current row and column
  void AcceptHisto(AliH2F * his);  //update buffer for - it will content histogram values
  AliH2F * GenerHisto();           //generate 2 dimensional histogram with digits
  AliH2F *  Draw( const char *option=0,Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); //draw digits
  Int_t GetSize() {return fNelems;} //return main buffer size
protected:
  virtual  void Invalidate();  
  void ExpandBuffer1(); //expand buffer of type to twodimensional array
  void CompresBuffer1(); //compres buffer according  algorithm 1
  Bool_t First0();  //first for the buffer type 0
  Bool_t Next0();  //next for the buffer type 0
  Bool_t First1(); //first for the buffer type 1
  Bool_t Next1();//next for the buffer type 1
  Short_t  GetDigit1(Int_t row, Int_t column); //return digit for given row and column
 
  Int_t     fNrows;   //number of rows in Segment
  Int_t     fNcols; //number of collumns in Segment 
private:
  AliArrayS *fElements;  //buffer of 2 bytes integers for digits
  AliArrayI *fIndex;  //index position of column
  Int_t     fBufType; //type of the buffer - define compression algorithm  
  Int_t     fThreshold; //treshold for zero suppresion
  Int_t     fNelems;  //total number of elements 
  Int_t fCurrentRow;   //!current row  iteration
  Int_t fCurrentCol;   //!current column iteration
  Int_t fCurrentIndex; //!current index in field
 
  ClassDef(AliDigits,1) 
};
 


inline Bool_t AliDigits::BoundsOK(const char *where, Int_t row, Int_t col) 
{
  //Check If Bound Ok
  if ( (col>=fNcols) || (col<0) ) return OutOfBoundsError(where,row,col);
  Int_t index =(*fIndex).At(col)+row;
  if ( (index<0) || (index>fNelems)) return OutOfBoundsError(where,row,col);
  return kTRUE;  
}

Short_t AliDigits::GetDigitFast(Int_t row, Int_t column)
{
  //
  //return digit from  fDigits array
  //if out of range return dummy value  ( value at row = 0, collumn = 0)
  //
  return fElements->At(fIndex->At(column)+row); 
}

void  AliDigits::SetDigitFast(Short_t value, Int_t row, Int_t column)
{
  //
  //set  digit 
  //
  if ( (row<0) || (row>=fNrows)  || (column<0) || (column>=fNrows) ) 
       ::Error("AliDigits::SetDigitFast", "row %d  col %d out of bounds (size: %d x %d, this: 0x%08x)", 
	   row, column, fNrows, fNcols, this);
  (*fElements)[fIndex->At(column)+row]=value; 
}

#endif

