/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////
// Class used to store the CTP (trigger) DDL raw data      //
/////////////////////////////////////////////////////////////

#ifndef AliTPCDDLRAWDATA_H
#define AliTPCDDLRAWDATA_H


class AliCTPRawData:public TObject{
 public:
  AliCTPRawData();                                       // default constructor
  virtual ~AliCTPRawData() {;}                           // destructor
  AliCTPRawData(const AliCTPRawData &source);            // copy constructor
  AliCTPRawData& operator=(const AliCTPRawData &source); // assignment operator
  void  RawData();  //This method is used to create the slides (sequence of files)

  ClassDef(AliCTPRawData,1)
};
    
#endif
