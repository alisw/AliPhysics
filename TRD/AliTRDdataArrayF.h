#ifndef TRDdataArrayF_H
#define TRDdataArrayF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdataArrayF.h,v */
 
#include   "AliTRDdataArray.h"

/////////////////////////////////////////////////////////////
//                                                         //
//  General container for float data from TRD detector     //
//  segments.                                              //
//  Adapted from AliDigits, origin M.Ivanov                //
//                                                         //
/////////////////////////////////////////////////////////////

class AliTRDdataArrayF : public AliTRDdataArray {

 public:

  AliTRDdataArrayF();
  AliTRDdataArrayF(Int_t nrow, Int_t ncol,Int_t ntime);
  ~AliTRDdataArrayF();

  virtual void    Allocate(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void    Compress(Int_t bufferType, Float_t threshold);
  virtual void    Compress(Int_t bufferType); 
  virtual void    Expand();
  virtual Bool_t  First();
  virtual Bool_t  Next(); 
  virtual void    Reset();

  inline  void    SetData(Int_t row, Int_t col, Int_t time, Float_t value);
  virtual void    SetThreshold(Float_t threshold) { fThreshold = threshold; };

  virtual Float_t GetData(Int_t row, Int_t col, Int_t time);
  virtual Float_t GetThreshold()                  { return fThreshold;  };

  virtual Int_t   GetSize();
  virtual Int_t   GetDataSize(); 
  virtual Int_t   GetOverThreshold(Float_t threshold);  

 protected:

  inline  void    SetDataFast(Int_t idx1, Int_t idx2, Float_t value); 
  inline  Float_t GetDataFast(Int_t idx1, Int_t idx2); 

  Float_t         GetData1(Int_t idx1, Int_t idx2); 
  void            Expand1(); 
  void            Compress1(); 
  void            Expand2();
  void            Compress2();
  Bool_t          First0();
  Bool_t          Next0(); 
  Bool_t          First1();
  Bool_t          Next1();

  AliTRDarrayF  *fElements;        // Buffer of 4 bytes floats for the array content
  Float_t        fThreshold;       // Threshold for zero suppression
 
  ClassDef(AliTRDdataArrayF,1)     // Container for float data of one TRD detector segment

};
 

//_____________________________________________________________________________
inline Float_t AliTRDdataArrayF::GetDataFast(Int_t idx1, Int_t idx2)
{
  //
  // Returns the value at a given position in the array
  //

  return fElements->At(fIndex->At(idx2) + idx1); 

}

//_____________________________________________________________________________
inline void AliTRDdataArrayF::SetData(Int_t row, Int_t col, Int_t time
                                                          , Float_t value)
{
  //
  // Sets the data value at a given position of the array
  // Includes boundary checking
  //

  if ((row >= 0) && (col >= 0) && (time >= 0)) {
    Int_t idx1 = GetIdx1(row,col);
    if ((idx1 >= 0) && (time < fNdim2)) {
      SetDataFast(idx1,time,value);
    }
    else {
      if (idx1 >= 0) {
        TObject::Error("SetData"
                      ,"time %d out of bounds (size: %d, this: 0x%08x)"
                      ,time,fNdim2,this);
      }
    }
  }

}

//_____________________________________________________________________________
inline void  AliTRDdataArrayF::SetDataFast(Int_t idx1, Int_t idx2, Float_t value)
{
  //
  // Set the value at a given position in the array
  //

  if ((idx1 < 0) || (idx1 >= fNdim1) || 
      (idx2 < 0) || (idx2 >= fNdim2)) { 
    TObject::Error("SetDataFast"
                  ,"idx1 %d  idx2 %d out of bounds (size: %d x %d, this: 0x%08x)"
                  ,idx1,idx2,fNdim1,fNdim2,this);
  }

  (*fElements)[fIndex->fArray[idx2] + idx1] = value; 

}

#endif

