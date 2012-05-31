/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
// Class used for read-write the ALTRO data format //
/////////////////////////////////////////////////////

/*This class is an interface between the altro format file and the 
  user, and can be used in write or read mode
  In the write mode a new altro file is created and filled using the method FillBuffer().
  The name of the file is specified as parameter in the constructor as well as the type mode.
  In the Read mode the specified file is open and the values can be read using the
  methods GetNext() and GetNextBackWord().
  The first method is used to read the file forward while the second is used to read backward 
*/

#ifndef AliALTROBUFFERV3_H
#define AliALTROBUFFERV3_H

#include "AliAltroBuffer.h"

class AliAltroBufferV3: public AliAltroBuffer {
 public:
  AliAltroBufferV3(const char* fileName, AliAltroMapping *mapping = NULL);
  virtual ~AliAltroBufferV3();

  virtual void  FillBuffer(Int_t val);
  //this method stores a word into the buffer

  virtual void  WriteTrailer(Int_t wordsNumber, Short_t hwAddress); 
  //this method is used to write the trailer

  virtual UChar_t WriteRCUTrailer(Int_t rcuId);
  //this method is used to write the RCU trailer

  void  SetFECERRA(UInt_t v) { fFECERRA = v; }
  void  SetFECERRB(UInt_t v) { fFECERRB = v; }
  void  SetERRREG2(UInt_t v) { fERRREG2 = v; }
  void  SetERRREG3(UInt_t v) { fERRREG3 = v; }
  void  SetActiveFECsA(UShort_t m) { fActiveFECsA = m; }
  void  SetActiveFECsB(UShort_t m) { fActiveFECsB = m; }
  void  SetALTROCFG1(UInt_t cfg1) { fALTROCFG1 = cfg1; }
  void  SetALTROCFG2(UInt_t cfg2) { fALTROCFG2 = cfg2; }
  void  SetTSample(Double_t v) { fTSample = v; }
  void  SetL1Phase(Double_t v) { fL1Phase = v; }

  void SetNChAddrMismatch(UShort_t v) { SetField(fERRREG3, 0, 0xFFF, v); }
  void SetNChLengthMismatch(UShort_t v) { SetField(fERRREG3, 12, 0x1FFF, v); }
  
  void SetBaselineCorr(UChar_t v) { SetField(fALTROCFG1, 0, 0xF, v); }
  void SetPolarity(Bool_t v) { SetField(fALTROCFG1, 4, 0x1, v); }
  void SetNPresamples(UChar_t v) { SetField(fALTROCFG1, 5, 0x3, v); }
  void SetNPostsamples(UChar_t v) { SetField(fALTROCFG1, 7, 0xF, v); }
  void SetSecondBaselineCorr(Bool_t v) { SetField(fALTROCFG1, 11, 0x1, v); }
  void SetGlitchFilter(UChar_t v) { SetField(fALTROCFG1, 12, 0x3, v); }
  void SetNNonZSPostsamples(UChar_t v) { SetField(fALTROCFG1, 14, 0x7, v); }
  void SetNNonZSPresamples(UChar_t v) { SetField(fALTROCFG1, 17, 0x3, v); }
  void SetZeroSupp(Bool_t v) { SetField(fALTROCFG1, 19, 0x1, v); }
  void SetNAltroBuffers(Bool_t v) { SetField(fALTROCFG2, 24, 0x1, v); }
  void SetNPretriggerSamples(UChar_t  v) { SetField(fALTROCFG2, 20, 0xF, v); }
  void SetNSamplesPerCh(UShort_t v) { SetField(fALTROCFG2, 10, 0x3FF, v); }
  void SetSparseRO(Bool_t v) { SetField(fALTROCFG2, 9, 0x1, v); }


  enum { kMaxWords = 1024 };

 protected:
  UInt_t        SetField(UInt_t& input, UShort_t start, UInt_t mask, UInt_t val) const;
  void          ReverseAndWrite();
  //this method reverse the altro data order and write the buffer to the file

  AliAltroBufferV3(const AliAltroBufferV3& source);
  AliAltroBufferV3& operator = (const AliAltroBufferV3& source);

  UShort_t fArray[kMaxWords]; // Temporary array needed in reverting data order
  Int_t    fN;                // Size of the temporary array
  UInt_t   fFECERRA;	// FECERRA
  UInt_t   fFECERRB;	// FECERRB
  UInt_t   fERRREG2;	// ERRREG2
  UInt_t   fERRREG3;	// ERRREG3
  UShort_t fActiveFECsA;// ActiveFECsA
  UShort_t fActiveFECsB;// ActiveFECsB
  UInt_t   fALTROCFG1;	// ALTROCFG1
  UInt_t   fALTROCFG2;	// ALTROCFG2
  Double_t fTSample;	// TSample
  Double_t fL1Phase;	// L1Phase

  ClassDef(AliAltroBufferV3,0)  // Interface to the Altro format
};

#endif
