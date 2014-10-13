/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
// Class used for the ALICE data challenges        //
/////////////////////////////////////////////////////

#ifndef AliTPCDDLRAWDATA_H
#define AliTPCDDLRAWDATA_H


class AliTPCDDLRawData:public TObject{
 public:
  AliTPCDDLRawData():TObject(),
    fVerbose(0){}//default constructor
  virtual ~AliTPCDDLRawData(){;}//destructor
  AliTPCDDLRawData(const AliTPCDDLRawData &source); // copy constructor
  AliTPCDDLRawData& operator=(const AliTPCDDLRawData &source); // ass. op.
  void  RawData(const char* inputFileName = "AliTPCDDL.dat");
  //This method is used to create the slides (sequence of files)
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
 private:
  Int_t fVerbose;         //Verbose level 0:Silent, 1: cout msg, 2:txt files for debugging
  enum {kDDLOffset = 0};  //offset for DDL number
  ClassDef(AliTPCDDLRawData,1)
};
    
#endif
