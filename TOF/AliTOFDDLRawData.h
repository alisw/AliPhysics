#ifndef ALITOFDDLRAWDATA_H
#define ALITOFDDLRAWDATA_H

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Class used for generating the files containing raw data,               //
// required for Data Challenge                                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
//#include "AliTOFCableLengthMap.h"

class TBranch;
class TClonesArray;

class AliTOFDigitMap;
//class AliTOFRawStream;

class AliTOFDDLRawData:public TObject {

 public:

  AliTOFDDLRawData();                               // default constructor
  virtual ~AliTOFDDLRawData();                              // destructor
  AliTOFDDLRawData(const AliTOFDDLRawData &source); // copy constructor
  AliTOFDDLRawData& operator=(const AliTOFDDLRawData &source); // ass. op.
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;} // To set the verbose level

  Int_t RawDataTOF(TBranch* branch); 

  void SetPackedAcquisitionMode(Bool_t mode) {fPackedAcquisition=mode;};
  void SetFakeOrphaneProduction(Bool_t flag) {fFakeOrphaneProduction=flag;};
  void SetMatchingWindow(Int_t matWin) {fMatchingWindow=matWin;}; // setter for fMatchingWindow [bin number]
  Bool_t GetPackedAcquisitionMode() const {return fPackedAcquisition;};
  Bool_t GetFakeOrphaneProduction() const {return fFakeOrphaneProduction;};
  Int_t  GetMatchingWindow() const {return fMatchingWindow;}; // getter for fMatchingWindow [bin number]

 private:

  void  GetDigits();

  void  ReverseArray(UInt_t a[], Int_t n) const;

  void  MakeDRMheader(Int_t nDDL, UInt_t *buf);
  void  MakeDRMtrailer(UInt_t *buf);
  void  MakeLTMheader(UInt_t *buf);
  void  MakeLTMdata(UInt_t *buf);
  void  MakeLTMtrailer(UInt_t *buf);
  void  MakeTRMheader(Int_t nTRM, UInt_t *buf);
  void  MakeTRMtrailer(UInt_t *buf);
  void  MakeTRMfiller(UInt_t *buf);
  void  MakeTRMchainHeader(Int_t iChain, Int_t nTRM, UInt_t *buf);
  void  MakeTRMchainTrailer(Int_t iChain, UInt_t *buf);
  void  MakeTDCdigits(Int_t nDDL, Int_t nTRM, Int_t iChain, UInt_t *buf);

  UInt_t MakeFiller() const;

  Bool_t HeadOrTail() const;

  Int_t fVerbose;                 //Verbose level (0:no msg, 1:msg,
				  //2:digits in txt files)
  Int_t fIndex;                   //number of 32-bit words to be
				  //stored into the output file
  Bool_t fPackedAcquisition;      //flag for packed/no packed acquisition
  Bool_t fFakeOrphaneProduction;  //flag to insert fake orphane
				  //(leading or trailing) time
				  //measurements
  Int_t fMatchingWindow;          //time window [bin number] where to
				  //search time-of-flight measurements
				  //for the current event

  AliTOFDigitMap *fTOFdigitMap;   //Pointer to the channel-TOF map

  TClonesArray *fTOFdigitArray;   //Pointer to the TOF digits

  Int_t fWordsPerDRM;
  Int_t fWordsPerTRM;
  Int_t fWordsPerChain;

  ClassDef(AliTOFDDLRawData,4)

};
    
#endif // ALITOFDDLRAWDATA_H
