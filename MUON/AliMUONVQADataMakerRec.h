#ifndef ALIMUONVQADATAMAKERREC_H
#define ALIMUONVQADATAMAKERREC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVQADataMakerRec
/// \brief Interface for a MUON QADataMakerRec
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ALIRECOPARAM_H
#  include "AliRecoParam.h"
#endif

#ifndef ALIQAv1_H
#  include "AliQAv1.h"
#endif

#include "AliQADataMakerRec.h"

class AliESDEvent;
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
  
  /// Initialization for handling Digits
  virtual void InitDigits() = 0; 
  /// Initialization for handling ESD
  virtual void InitESDs() = 0; 
  /// Initialization for handling Raws
  virtual void InitRaws() = 0; 
  /// Initialization for handling RecPoints
  virtual void InitRecPoints() = 0; 
  
  /// Produces QA data for Raws
  virtual void MakeRaws(AliRawReader* rawReader) = 0; 
  /// Produces QA data for Digits
  virtual void MakeDigits(TTree* dig) = 0; 
  /// Produces QA data for RecPoints
  virtual void MakeRecPoints(TTree* recpo) = 0;
  /// Produces QA data for ESD
  virtual void MakeESDs(AliESDEvent* esd) = 0;
  
  /// Wrap up things at each cycle for Raws
  virtual void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list) = 0;
  /// Wrap up things at each cycle for RecPoints
  virtual void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list) = 0;
  /// Wrap up things at each cycle for ESD
  virtual void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list) = 0;
  /// Wrap up things at each cycle for Digits
  virtual void EndOfDetectorCycleDigits(Int_t specie, TObjArray** list) = 0;

  /// Reset anything that must be reset for Raws
  virtual void ResetDetectorRaws(TObjArray* list) { ResetDetector(list); }
  /// Reset anything that must be reset for RecPoints
  virtual void ResetDetectorRecPoints(TObjArray* list) { ResetDetector(list); }
  /// Reset anything that must be reset for ESD
  virtual void ResetDetectorESDs(TObjArray* list) { ResetDetector(list); }
  /// Reset anything that must be reset for Digits
  virtual void ResetDetectorDigits(TObjArray* list) { ResetDetector(list); }
  ///
  void   ForbidCloning(TH1* h, Bool_t v=kTRUE)  { if (fMaster) fMaster->ForbidCloning(h,v); }
  ///  
protected:

  void ResetDetector(const TObjArray* list);
  
  Int_t RunNumber() const;
  
  AliRecoParam::EventSpecie_t CurrentEventSpecie() const;
  
  const AliMUONRecoParam* GetRecoParam() const;
  
  TObject* GetDigitsData(Int_t index)    const  {return fMaster ? fMaster->GetDigitsData(index)    : 0x0;}
  TObject* GetESDsData(Int_t index)      const  {return fMaster ? fMaster->GetESDsData(index)      : 0x0;}
  TObject* GetRecPointsData(Int_t index) const  {return fMaster ? fMaster->GetRecPointsData(index) : 0x0;}
  TObject* GetRawsData(Int_t index)      const  {return fMaster ? fMaster->GetRawsData(index)      : 0x0;}
  
  TH1*     GetDigitsData(Int_t index, int trCl)    const  {return fMaster ? fMaster->GetDigitsData(index, trCl)    : 0x0;}
  TH1*     GetESDsData(Int_t index, int trCl)      const  {return fMaster ? fMaster->GetESDsData(index, trCl)      : 0x0;}
  TH1*     GetRecPointsData(Int_t index, int trCl) const  {return fMaster ? fMaster->GetRecPointsData(index, trCl) : 0x0;}
  TH1*     GetRawsData(Int_t index, int trCl)      const  {return fMaster ? fMaster->GetRawsData(index,trCl)       : 0x0;}
  
  TH1*     GetMatchingRawsHisto(Int_t index, Int_t trigId) const {return fMaster ? fMaster->GetMatchingRawsHisto(index, trigId) : 0x0;}
  
  Int_t Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);  
  Int_t Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);                      
  Int_t Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE);                 
  Int_t Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert = kFALSE, const Bool_t image = kFALSE, const Bool_t saveForCorr = kFALSE);  
  //
  Int_t  FillDigitsData(Int_t index, double x)                 {return fMaster ? fMaster->FillDigitsData(index,x) : 0;}
  Int_t  FillRawsData(Int_t index, double x)                   {return fMaster ? fMaster->FillRawsData(index,x) : 0;}
  Int_t  FillRecPointsData(Int_t index, double x)              {return fMaster ? fMaster->FillRecPointsData(index,x) : 0;}
  Int_t  FillESDsData(Int_t index, double x)                   {return fMaster ? fMaster->FillESDsData(index,x) : 0;}
  //
  Int_t  FillDigitsData(Int_t index, double x, double y)       {return fMaster ? fMaster->FillDigitsData(index,x,y) : 0;}
  Int_t  FillRawsData(Int_t index, double x, double y)         {return fMaster ? fMaster->FillRawsData(index,x,y) : 0;}
  Int_t  FillRecPointsData(Int_t index, double x, double y)    {return fMaster ? fMaster->FillRecPointsData(index,x,y) : 0;}
  Int_t  FillESDsData(Int_t index, double x, double y)         {return fMaster ? fMaster->FillESDsData(index,x,y) : 0;}
  //
  Int_t  FillDigitsData(Int_t index, double x, double y, double w)     {return fMaster ? fMaster->FillDigitsData(index,x,y,w) : 0;}
  Int_t  FillRawsData(Int_t index, double x, double y, double w)       {return fMaster ? fMaster->FillRawsData(index,x,y,w) : 0;}
  Int_t  FillRecPointsData(Int_t index, double x, double y, double w)  {return fMaster ? fMaster->FillRecPointsData(index,x,y,w) : 0;}
  Int_t  FillESDsData(Int_t index, double x, double y, double w)       {return fMaster ? fMaster->FillESDsData(index,x,y,w) : 0;}
  //
  void   ClonePerTrigClass(AliQAv1::TASKINDEX_t task);  

  AliQADataMakerRec* fMaster; ///< master to get access to its methods

private:
  /// Not implemented
  AliMUONVQADataMakerRec(const AliMUONVQADataMakerRec& rhs);
  /// Not implemented
  AliMUONVQADataMakerRec& operator=(const AliMUONVQADataMakerRec& rhs);
  
  ClassDef(AliMUONVQADataMakerRec,1) // Interface for a MUON QADataMakerRec
};

#endif
