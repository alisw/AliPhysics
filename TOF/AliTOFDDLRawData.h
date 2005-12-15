/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// Class used for generating the files containing raw data,               //
// required for  Data Challenge                                           //
////////////////////////////////////////////////////////////////////////////

#ifndef AliTOFDDLRAWDATA_H
#define AliTOFDDLRAWDATA_H

class AliTOF;
class AliTOFGeometry;
class TTree;

class AliTOFDDLRawData:public TObject{
 public:
  AliTOFDDLRawData();                               // default constructor
  AliTOFDDLRawData(AliTOFGeometry *tofGeom);        // constructor
  virtual ~AliTOFDDLRawData(){;}                    // destructor
  AliTOFDDLRawData(const AliTOFDDLRawData &source); // copy constructor
  AliTOFDDLRawData& operator=(const AliTOFDDLRawData &source); // ass. op.
  Int_t RawDataTOF(TBranch* branch); 
  // This method generates the files with the TOF detector data
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
  // To set the verbose level
 private:
  void  GetDigits(TClonesArray *TOFdigits, Int_t ddl,UInt_t *buf);
  //This method formats and stores in buf all the digits of a TOF module

  Int_t fVerbose;               //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  Int_t fIndex;                 //number of 32 words to be stored into the output file
  AliTOFGeometry *fTOFgeometry; //

  ClassDef(AliTOFDDLRawData,1)

};
    
#endif
