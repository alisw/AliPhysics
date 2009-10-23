#ifndef ALIMUONVQADATAMAKERREC_H
#define ALIMUONVQADATAMAKERREC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVQADataMakerRec
/// \brief Interface for a MUON QADataMakerRec
/// 
/// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ALIRECOPARAM_H
#  include "AliRecoParam.h"
#endif

class AliESDEvent;
class AliQADataMakerRec;
class AliMUONRecoParam;
class AliRawReader;
class TH1;
class TObjArray;
class TTree;

class AliMUONVQADataMakerRec : public TObject
{
public:
  AliMUONVQADataMakerRec(AliQADataMakerRec* master);
  virtual ~AliMUONVQADataMakerRec();
  
  virtual void InitDigits() = 0; 
  virtual void InitESDs() = 0; 
  virtual void InitRaws() = 0; 
  virtual void InitRecPoints() = 0; 
  
  virtual void MakeRaws(AliRawReader* rawReader) = 0; 
  virtual void MakeDigits(TTree* dig) = 0; 
  virtual void MakeRecPoints(TTree* recpo) = 0;
  virtual void MakeESDs(AliESDEvent* esd) = 0;
  
  virtual void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list) = 0;
  virtual void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list) = 0;
  virtual void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list) = 0;
  virtual void EndOfDetectorCycleDigits(Int_t specie, TObjArray** list) = 0;

protected:

  Int_t RunNumber() const;
  
  AliRecoParam::EventSpecie_t CurrentEventSpecie() const;
  
  const AliMUONRecoParam* GetRecoParam() const;
  
  TH1* GetDigitsData(Int_t index) const;
  TH1* GetESDsData(Int_t index) const;
  TH1* GetRecPointsData(Int_t index) const;
  TH1* GetRawsData(Int_t index) const;
                    
  Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);  
  Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);                      
  Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);                 
  Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE);  
  
private:
  AliMUONVQADataMakerRec(const AliMUONVQADataMakerRec& rhs);
  AliMUONVQADataMakerRec& operator=(const AliMUONVQADataMakerRec& rhs);
  
  AliQADataMakerRec* fMaster; ///< master to get access to its methods
  
  ClassDef(AliMUONVQADataMakerRec,1) // Interface for a MUON QADataMakerRec
};

#endif
