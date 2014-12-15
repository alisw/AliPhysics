#ifndef ALIDIGITS_H
#define ALIDIGITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment digits
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include   <TArrayI.h>
#include   <TArrayS.h>
#include   "AliSegmentID.h"
class AliH2F;

class AliDigits: public AliSegmentID{ 
public:
  AliDigits();
  AliDigits(const AliDigits &digits); //copy constructor
  AliDigits &operator = (const AliDigits & digits); //assignment operator
  virtual ~AliDigits();
  Short_t * GetDigits(){return fElements->GetArray();}   //return row  pointer to the array digits
  Short_t GetDigitFast(Int_t row, Int_t column);  //return value at given row and collumn
  Short_t GetDigitUnchecked(Int_t row, Int_t column);  //return value at given row and collumn
  void  SetDigitFast(Short_t value,Int_t row, Int_t column);  //set value at given row and collumn
  Bool_t BoundsOK(const char *where, Int_t row, Int_t col) ;  //Check If Bound Ok
  Bool_t OutOfBoundsError(const char *where, Int_t row, Int_t column);
  virtual void Allocate(Int_t rows, Int_t columns);  //construct empty buffer fDigits with size rows x columns
  virtual Short_t GetDigit(Int_t row, Int_t column);
  virtual void ExpandBuffer();  //expand buffer to twodimensional array
  virtual void CompresBuffer(Int_t bufferType,Int_t threshold); //compres buffer according buffertype algorithm   
  virtual Bool_t First(); //adjust  first valid current digit
  virtual Bool_t Next();  //addjust next valid current digit
  void SetThreshold(Int_t th) {fThreshold = th;} //set threshold
  Int_t  GetThreshold() {return fThreshold;}  //return threshold    
  Int_t GetNRows(){return fNrows;}
  Int_t GetNCols(){return fNcols;}
  Int_t CurrentRow(){ return fCurrentRow;}  //return current row
  Int_t CurrentColumn(){ return fCurrentCol;} //return current column
  Int_t CurrentDigit() {return fElements->At(fCurrentIndex);} //return degit for current row and column
  void AcceptHisto(AliH2F * his);  //update buffer for - it will content histogram values
  AliH2F * GenerHisto();           //generate 2 dimensional histogram with digits
  AliH2F *DrawDigits( const char *option=0,Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); //draw digits
  
  Int_t GetSize();//return total size of object in bytes
  Int_t GetDigitSize(); //return total size of pure digits 
  Int_t GetOverTh(Float_t threshold,Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); //return number of digits over threshold 

  inline Short_t * GetDigitsColumn(Int_t row);                              //return row  pointer to the array digits

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
  TArrayS *fElements;  //buffer of 2 bytes integers for digits
  TArrayI *fIndex;  //index position of column
  Int_t     fBufType; //type of the buffer - define compression algorithm  
  Int_t     fThreshold; //treshold for zero suppresion
  Int_t     fNelems;  //total number of elements 
  Int_t fCurrentRow;   //!current row  iteration
  Int_t fCurrentCol;   //!current column iteration
  Int_t fCurrentIndex; //!current index in field
 
  ClassDef(AliDigits,2) 
};
 


inline Bool_t AliDigits::BoundsOK(const char *where, Int_t row, Int_t col) 
{
  //Check If Bound Ok
  if ( (col>=fNcols) || (col<0) ) return OutOfBoundsError(where,row,col);
  Int_t index =(*fIndex).At(col)+row;
  if ( (index<0) || (index>fNelems)) return OutOfBoundsError(where,row,col);
  return kTRUE;  
}

inline Short_t AliDigits::GetDigitFast(Int_t row, Int_t column)
{
  //
  //return digit from  fDigits array
  //if out of range return dummy value  ( value at row = 0, collumn = 0)
  //
  return fElements->At(fIndex->At(column)+row); 
}

inline Short_t AliDigits::GetDigitUnchecked(Int_t row, Int_t column)
{
  //
  //return digit from  fDigits array
  //if out of range return dummy value  ( value at row = 0, collumn = 0)
  //
  return fElements->fArray[fIndex->fArray[column]+row]; 
}

inline Short_t * AliDigits::GetDigitsColumn(Int_t column){
  //
  //return row  pointer to the array digits
  //
  return &(fElements->fArray[fIndex->fArray[column]]);
}


inline void  AliDigits::SetDigitFast(Short_t value, Int_t row, Int_t column)
{
  //
  //set  digit 
  //
  if ( (row<0) || (row>=fNrows)  || (column<0) || (column>=fNcols) ) 
       Error("AliDigits::SetDigitFast", "row %d  col %d out of bounds (size: %d x %d, this: 0x%08lx)", 
	     row, column, fNrows, fNcols, (ULong_t)this);
  (*fElements)[fIndex->At(column)+row]=value; 
}

#endif

