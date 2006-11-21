#ifndef ALIPMDRAWTOSDIGITS_H
#define ALIPMDRAWTOSDIGITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : October 06 2006                              //
//                                                     //
//-----------------------------------------------------//


class TClonesArray;
class TTree;

class AliLoader;
class AliRunLoader;
class AliRawReader;

class AliPMDsdigit;
class AliPMDdigit;

class AliPMDRawToSDigits : public TObject
{
 public:

  AliPMDRawToSDigits();
  AliPMDRawToSDigits(const AliPMDRawToSDigits & /* pmdr2sd */);     // copy constructor
  AliPMDRawToSDigits &operator=(const AliPMDRawToSDigits & /* pmdr2sd */); // assi op
  virtual ~AliPMDRawToSDigits();
  
  void Raw2SDigits(AliRunLoader *runLoader, AliRawReader *rawReader);
  void Raw2Digits(AliRunLoader *runLoader, AliRawReader *rawReader);
  void AdcToMeV(Int_t adc, Float_t &edep);
  void AddSDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
		 Int_t irow, Int_t icol, Float_t adc);
  void AddDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
		Int_t irow, Int_t icol, Float_t adc);

  void ResetSDigit();
  void ResetDigit();


 protected:
  TClonesArray *fSDigits;    //! List of digits
  TClonesArray *fDigits;     //! List of digits

  Int_t   fNsdigit;          // Digits counter
  Int_t   fNdigit;           // Digits counter

  ClassDef(AliPMDRawToSDigits,0)    // Coverts Raw to SDigits
};
#endif

