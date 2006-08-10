#ifndef AliTOFDDLRAWDATA_H
#define AliTOFDDLRAWDATA_H

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Class used for generating the files containing raw data,               //
// required for Data Challenge                                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TBranch;
class TClonesArray;

class AliTOFDigitMap;
class AliTOFGeometry;
class AliTOFRawStream;

class AliTOFDDLRawData:public TObject {

 public:

  AliTOFDDLRawData();                               // default constructor
  AliTOFDDLRawData(AliTOFGeometry *tofGeom);        // constructor
  virtual ~AliTOFDDLRawData(){;}                    // destructor
  AliTOFDDLRawData(const AliTOFDDLRawData &source); // copy constructor
  AliTOFDDLRawData& operator=(const AliTOFDDLRawData &source); // ass. op.
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;} // To set the verbose level

  Int_t RawDataTOF(TBranch* branch); 

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
  void  MakeTRMfiller(UInt_t *buf, UInt_t nWordsPerTRM);
  void  MakeTRMchainHeader(Int_t iChain, Int_t nTRM, UInt_t *buf);
  void  MakeTRMchainTrailer(Int_t iChain, UInt_t *buf);
  void  MakeTDCdigits(Int_t nDDL, Int_t nTRM, Int_t iChain, UInt_t *buf, UInt_t &nWordsPerTRM);

  UInt_t  MakeFiller();

  Int_t fVerbose;                 //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  Int_t fIndex;                   //number of 32-bit words to be stored into the output file

  AliTOFGeometry *fTOFgeometry;   //Pointer to the TOF geometry

  AliTOFDigitMap *fTOFdigitMap;   //Pointer to the channel-TOF map

  TClonesArray *fTOFdigitArray;   //Pointer to the TOF digits

  AliTOFRawStream *fTOFrawStream; //Pointer to the AliTOFRawStream class

  ClassDef(AliTOFDDLRawData,1)

};
    
#endif
