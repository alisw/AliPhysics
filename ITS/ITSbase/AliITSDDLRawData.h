#ifndef AliITSDDLRAWDATA_H
#define AliITSDDLRAWDATA_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////////
// Class used for generating the files containing raw data,                //
// required for  Data Challenge                                            //
/////////////////////////////////////////////////////////////////////////////

/* $Id$ */

class TTree;
class TBranch;
class AliITSDDLModuleMapSDD;
class AliITSFOSignalsSPD;

class AliITSDDLRawData:public TObject{
 public:
  AliITSDDLRawData();//default constructor
  virtual ~AliITSDDLRawData(){;}//destructor
  AliITSDDLRawData(const AliITSDDLRawData &source); // copy constructor
  AliITSDDLRawData& operator=(const AliITSDDLRawData &source); // ass. op.
  Int_t RawDataSPD(TBranch* branch, AliITSFOSignalsSPD* foSignals = NULL);
  // This method generates the files with the Silicon pixel detector data
  Int_t RawDataSDD(TBranch* branch, const AliITSDDLModuleMapSDD* ddlsdd); 
  // This method generates the files with the Silicon drift detector data
  Int_t RawDataSSD(TBranch* branch, UInt_t runno);
  // This method generates the files with the Silicon pixel detector data
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
  // To set the verbose level
  void SetSDDRawFormat(Char_t val=7){
    fSDDRawFormat=val;
  }
 private: 
  void  GetDigitsSPD(TClonesArray *ITSdigits, Int_t mod,Int_t ddl,UInt_t *buf, AliITSFOSignalsSPD* foSignals = NULL);
  //This method formats and stores in buf all the digits of a SPD module
  void  GetDigitsSDDCompressed(TClonesArray *ITSdigits, Int_t mod,UInt_t *buf);
  void  GetDigitsSDD(TClonesArray *ITSdigits, Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf);
  //This method formats and stores in buf all the digits of a SDD module
  void  GetDigitsSSD(TClonesArray *ITSdigits, Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf,UInt_t runno);
  //This method formats and stores in buf all the digits of a SSD module
  void  WriteChipHeader(Int_t ChipAddr,Int_t halfStave,UInt_t &BaseWord);
  void  WriteChipTrailer(UInt_t *buf, Int_t ChipHitCount, Bool_t foBit, UInt_t &BaseWord);
  void  WriteHit(UInt_t *buf,Int_t RowAddr,Int_t HitAddr,UInt_t &BaseWord);
  //The three previous  methods are used to store the data according to the 
  //Silicon pixel detector data format
  Int_t fVerbose;            //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  Int_t fIndex;             //number of 32 words to be stored into the output file
  Int_t fHalfStaveModule;     //first or second half of an Half Stave module
  Char_t fSDDRawFormat;          // index for SDD raw format

  ClassDef(AliITSDDLRawData,3)
};
    
#endif
