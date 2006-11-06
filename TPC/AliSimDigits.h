#ifndef ALISIMDIGITS_H
#define ALISIMDIGITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class generaol Alice segment digits
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////
#include <TError.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include   "AliDigits.h"

class AliH2F;


class AliSimDigits : public AliDigits{
public: 
  AliSimDigits();
  AliSimDigits(const AliSimDigits &param);
  AliSimDigits &operator = (const AliSimDigits & digits); 
  virtual ~AliSimDigits();
  void AllocateTrack(Int_t length);  //construct empty buffer fTracks with size rows x column x length (number of tracks for one digit)
  Int_t *GetTracks(){return fTracks->GetArray();}
  Int_t GetTrackIDFast(Int_t row, Int_t column,Int_t level);  //return track ID  at given row and collumn
  void  SetTrackIDFast(Int_t value,Int_t row, Int_t column,Int_t level);  //set ID track at given row and collumn
  virtual Int_t GetTrackID(Int_t row, Int_t column, Int_t level);
  virtual void ExpandTrackBuffer();  //expand buffer to twodimensional array
  virtual void CompresTrackBuffer(Int_t bufType); //compres buffer according buffertype algorithm 
  AliH2F *  DrawTracks( const char *option=0,Int_t level=0, 
		  Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); //draw tracks
  TClonesArray * GenerTPCClonesArray(TClonesArray * arr); //generate TClonnesArray of digits
  //only for demonstration purpose
private:
  void InvalidateTrack();
 
  Int_t GetTrackID1(Int_t row, Int_t column, Int_t level);  //returnb track ID of digits - for buffer compresion 1 
  void  ExpandTrackBuffer1(); //comress track according algorithm 1 (track ID comression independent to the digit compression) 
  void  CompresTrackBuffer1(); //comress track according algorithm 1 (track ID comression independent to the digit compression)
 
  Int_t GetTrackID2(Int_t row, Int_t column, Int_t level);  //returnb track ID of digits - for buffer compresion 2
  void  ExpandTrackBuffer2(); //comress track according algorithm 2 (track ID comression according  digit compression)
  void  CompresTrackBuffer2(); //comress track according algorithm 2 (track ID comression according  digit compression)

  TArrayI * fTracks;     //buffer of track index 
  TArrayI * fTrIndex;    //index position of column
  Int_t       fNlevel;   //number of tracks etries  for one digit
  Int_t       fTrBufType;  //buffer type of the tracks
  // Bool_t      ClassError( ); //signalize class error 
  ClassDef(AliSimDigits,2) 
};



inline Int_t AliSimDigits::GetTrackIDFast(Int_t row, Int_t column,Int_t level)
{
  //
  //return track ID  at given row and column
  //  return fTracks[level].At(fTrIndex[level][column]+row); 
  return fTracks->At(level*fNrows*fNcols+fNrows*column+row); 
}
 
inline void AliSimDigits::SetTrackIDFast(Int_t value,Int_t row, Int_t column,Int_t level)
{
  //
  value+=2;
  //set ID track at given row and collumn
  //  fTracks[level][fTrIndex[level][column]+row]=value; 
  if ( (row<0) || (row>=fNrows)  || (column<0) || (column>=fNcols) ) 
       ::Error("AliSimDigits::SetTrackIDFast", "row %d  col %d out of bounds (size: %d x %d, this: 0x%08x)", 
	   row, column, fNrows, fNcols, this);
  if ( (level<0) || (level>=fNlevel)) ::Error("AliSimDigits::SetTrackIDFast", "index %d out of bounds", level);
  (*fTracks)[level*fNrows*fNcols+fNrows*column+row]=value; 
}



#endif






