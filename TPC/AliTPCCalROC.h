#ifndef ALITPCCALROC_H
#define ALITPCCALROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalROC.h,v */

//////////////////////////////////////////////////
//                                              //
//  TPC calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>

//_____________________________________________________________________________
class AliTPCCalROC : public TObject {

 public:
  
  AliTPCCalROC();
  AliTPCCalROC(Int_t sector);
  AliTPCCalROC(const AliTPCCalROC &c);
  virtual           ~AliTPCCalROC();  
  Int_t        GetNrows() const                  { return fgNRows[fIndex]; };
  Int_t        GetNchannels()       const       { return fgNChannels[fIndex];   };
  Float_t      GetValue(Int_t row, Int_t pad)  { return fData[fgRowPosIndex[fIndex][row]+pad]; };
  void         SetValue(Int_t row, Int_t pad, Float_t vd)
                                                {  fData[fgRowPosIndex[fIndex][row]+pad]= vd; };
  static void Init(); 
 public:
  Int_t     fSector;          // sector number
  Int_t     fIndex;           // 0- if inner 1- outer
  Float_t  *fData;            //[fNchannels] Data
  //
  static Int_t  fgNSectorsAll;     // number of sectors
  static Int_t  fgNSectors[2];     // number of sectors - inner outer
  static Int_t  fgNRows[2];        // number of row     - inner outer
  static Int_t  fgNChannels[2];    // total number of pads   - inner sector - outer sector
  static Int_t *fgNPads[2];        // number of pads in row  - inner - outer      
  static Int_t *fgRowPosIndex[2];  // index array            - inner - outer
  // 
  ClassDef(AliTPCCalROC,1)    //  TPC ROC calibration class

};

#endif
