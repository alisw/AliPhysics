/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
// Class used for the ALICE data challenges        //
/////////////////////////////////////////////////////

#ifndef AliTPCDDLRAWDATA_H
#define AliTPCDDLRAWDATA_H


class AliTPCDDLRawData:public TObject{
 public:
  AliTPCDDLRawData(){fVerbose=0;}//default constructor
  virtual ~AliTPCDDLRawData(){;}//destructor
  AliTPCDDLRawData(const AliTPCDDLRawData &source); // copy constructor
  AliTPCDDLRawData& operator=(const AliTPCDDLRawData &source); // ass. op.
  void  RawData(Int_t LDCsNumber,Int_t EventNumber);
  //This method is used to create the slides (sequence of files)
  Int_t RawDataCompDecompress(Int_t LDCsNumber,Int_t EventNumber,Int_t Comp=0);
  //This method is used to create the compressed slides starting from the uncompressed ones 
  //or it can be used to decompress a sequence of compressed slices
  void  RawDataAltro()const;
  //This method is used to create the Altro format file from "AliTPCDDL.dat"
  void RawDataAltroDecode(Int_t LDCsNumber,Int_t EventNumber,Int_t Comp=0);
  //This method is used to construct an Altro format file starting from
  //the slices compressed or uncompressed
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
 private:
  Int_t fVerbose;         //Verbose level 0:Silent, 1: cout msg, 2:txt files for debugging
  ClassDef(AliTPCDDLRawData,1)
};
    
#endif
