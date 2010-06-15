/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONVQADataMakerRec.h"

///
/// \class AliMUONVQADataMakerRec
/// 
/// Interface for a MUON QADataMakerRec, common to MCH and MTR
/// 
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONVQADataMakerRec)
/// \endcond

#include "AliQADataMakerRec.h"
#include "AliMUONRecoParam.h"
#include "AliCDBManager.h"
#include "TH1.h"

//_____________________________________________________________________________
AliMUONVQADataMakerRec::AliMUONVQADataMakerRec(AliQADataMakerRec* master)
: fMaster(master)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVQADataMakerRec::~AliMUONVQADataMakerRec()
{
  /// dtor
}

//_____________________________________________________________________________
Int_t 
AliMUONVQADataMakerRec::Add2DigitsList(TH1 * hist, const Int_t index, const Bool_t expert , const Bool_t image )
{
  /// fwd
  return fMaster ? fMaster->Add2DigitsList(hist,index,expert,image) : -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONVQADataMakerRec::Add2ESDsList(TH1 * hist, const Int_t index, const Bool_t expert , const Bool_t image )
{
  /// fwd
  return fMaster ? fMaster->Add2ESDsList(hist,index,expert,image) : -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONVQADataMakerRec::Add2RecPointsList(TH1 * hist, const Int_t index, const Bool_t expert , const Bool_t image )
{
  /// fwd
  return fMaster ? fMaster->Add2RecPointsList(hist,index,expert,image) : -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONVQADataMakerRec::Add2RawsList(TH1 * hist, const Int_t index, const Bool_t expert , const Bool_t image , const Bool_t saveForCorr )
{
  /// fwd
  return fMaster ? fMaster->Add2RawsList(hist,index,expert,image,saveForCorr) : -1;
}

//_____________________________________________________________________________
AliRecoParam::EventSpecie_t 
AliMUONVQADataMakerRec::CurrentEventSpecie() const
{
  /// fwd
  return fMaster ? fMaster->GetEventSpecie() : AliRecoParam::kDefault;
}

//_____________________________________________________________________________
const AliMUONRecoParam* 
AliMUONVQADataMakerRec::GetRecoParam() const
{
  /// fwd
  return fMaster ? dynamic_cast<const AliMUONRecoParam*>(fMaster->GetRecoParam()) : 0x0;
}

//_____________________________________________________________________________
TH1* 
AliMUONVQADataMakerRec::GetDigitsData(Int_t index) const
{
  /// fwd
  return fMaster ? fMaster->GetDigitsData(index) : 0x0;
}

//_____________________________________________________________________________
TH1* 
AliMUONVQADataMakerRec::GetESDsData(Int_t index) const
{
  /// fwd
  return fMaster ? fMaster->GetESDsData(index) : 0x0;
}

//_____________________________________________________________________________
TH1* 
AliMUONVQADataMakerRec::GetRecPointsData(Int_t index) const
{
  /// fwd
  return fMaster ? fMaster->GetRecPointsData(index) : 0x0;
}

//_____________________________________________________________________________
TH1* 
AliMUONVQADataMakerRec::GetRawsData(Int_t index) const
{
  /// fwd
  return fMaster ? fMaster->GetRawsData(index) : 0x0;
}

//_____________________________________________________________________________
void 
AliMUONVQADataMakerRec::ResetDetector(const TObjArray* list)
{
  /// Reset all histograms found in list, that match either trigger or tracker

  TString cn(ClassName());
  TString pattern;
  
  if ( cn.Contains("Trigger") ) pattern = "Trigger";
  if ( cn.Contains("Tracker") ) pattern = "Tracker";
  
  TIter next(list); 
  TObject* o;
  while ( (o = next()) ) 
  {
    TH1* h = dynamic_cast<TH1*>(o);
    if ( h ) 
    {
      TString hcn(h->ClassName());
      if ( hcn.Contains(pattern) ) 
      {
        h->Reset();
      }
    }
  }
}

//_____________________________________________________________________________
Int_t 
AliMUONVQADataMakerRec::RunNumber() const
{
  /// fwd
  return fMaster ? fMaster->GetRun() : -1;
}
