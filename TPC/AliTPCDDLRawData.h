/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
// Class used for the fifth ALICE data challenge  //
/////////////////////////////////////////////////////

#ifndef AliTPCDDLRAWDATA_H
#define AliTPCDDLRAWDATA_H


class AliTPCDDLRawData:public TObject{
 public:
  AliTPCDDLRawData(){;}//default constructor
  virtual ~AliTPCDDLRawData(){;}//destructor
  AliTPCDDLRawData(const AliTPCDDLRawData &source); // copy constructor
  AliTPCDDLRawData& operator=(const AliTPCDDLRawData &source); // ass. op.
  //This method is used to create the slides (sequence of files)
  void  RawData(Int_t LDCsNumber);
  //This method is used to create the compressed slides starting from the uncompressed ones 
  //or it can be used to decompress a sequence of compressed slices
  Int_t RawDataCompDecompress(Int_t LDCsNumber,Int_t Comp=0);
  //This method is used to create the Altro format file from "AliTPCDDL.dat"
  void  RawDataAltro();
  //This method is used to Construct an Altro format file starting from
  //the slices compressed or uncompressed
  void RawDataAltroDecode(Int_t LDCsNumber,Int_t Comp=0);
 private:
  ClassDef(AliTPCDDLRawData,1)
};
    
#endif
