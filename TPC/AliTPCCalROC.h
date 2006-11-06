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
#include <AliTPCROC.h>

//_____________________________________________________________________________
class AliTPCCalROC : public TObject {

 public:
  
  AliTPCCalROC();
  AliTPCCalROC(UInt_t sector);
  AliTPCCalROC(const AliTPCCalROC &c);
 AliTPCCalROC &operator = (const AliTPCCalROC & param);
  virtual           ~AliTPCCalROC();  
  UInt_t        GetNrows() const               { return fNRows;};
  UInt_t        GetNchannels()       const     { return fNChannels;};
  UInt_t        GetNPads(UInt_t row)  const     { return (row<fNRows)? AliTPCROC::Instance()->GetNPads(fSector,row):0;};
  Float_t      GetValue(UInt_t row, UInt_t pad) const { return ( (row<fNRows) && (fIndexes[row]+pad)<fNChannels)? fData[fIndexes[row]+pad]: 0; };
  Float_t      GetValue(UInt_t channel) const { return  fData[channel]; };
  void         SetValue(UInt_t row, UInt_t pad, Float_t vd) { if ( row<fNRows && (fIndexes[row]+pad)<fNChannels)fData[fIndexes[row]+pad]= vd; };
  void         SetValue(UInt_t channel, Float_t vd) {fData[channel]= vd; };
  virtual void Draw(Option_t* option = "");
  static void Test();
 protected:
  UInt_t     fSector;          // sector number
  UInt_t     fNChannels;       // number of channels
  UInt_t     fNRows;           // number of rows
  const UInt_t* fIndexes;      //!indexes
  Float_t  *fData;            //[fNChannels] Data
  ClassDef(AliTPCCalROC,1)    //  TPC ROC calibration class

};

#endif
