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
  void  RawData(const char* inputFileName = "AliTPCDDL.dat");
  //This method is used to create the slides (sequence of files)
  Int_t RawDataCompDecompress(Bool_t compress = kTRUE);
  //This method is used to create the compressed slides starting from the uncompressed ones 
  //or it can be used to decompress a sequence of compressed slices
  void  RawDataAltro(const char* inputFileName = "AliTPCDDL.dat", const char* outputFileName = "AltroFormatDDL.dat")const;
  //This method is used to create the Altro format file from "AliTPCDDL.dat"
  void RawDataAltroDecode(const char* outputFileName);
  //This method is used to construct an Altro format file starting from
  //the slices compressed or uncompressed
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
 private:
  Int_t fVerbose;         //Verbose level 0:Silent, 1: cout msg, 2:txt files for debugging
  enum {kDDLOffset = 0};  //offset for DDL number
  ClassDef(AliTPCDDLRawData,1)
};
    
#endif
