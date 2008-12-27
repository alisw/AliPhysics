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

#include "AliMUONQAChecker.h"

/// \class AliMUONQAChecker
///
/// Implementation of AliQACheckerBase for MCH and MTR
///
/// For the moment we only implement the checking of raw data QA for the tracker
/// by looking at the occupancy at the manu level.
/// We count the number of manus above a given occupancy threshold, and
/// depending on that number, the resulting QA flag is warning, error or fatal.
/// (there's no "info" type in this case).
///
/// \author Laurent Aphecetche, Subatech

#include "AliLog.h"
#include "AliMUONVTrackerData.h"
#include "AliMpManuIterator.h"
#include "AliQA.h"
#include <TDirectory.h>
#include <TH1.h>

/// \cond CLASSIMP
ClassImp(AliMUONQAChecker)
/// \endcond

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker() : 
    AliQACheckerBase("MUON","MUON Quality Assurance Data Maker") 
{
	/// ctor
}          

//__________________________________________________________________
AliMUONQAChecker::~AliMUONQAChecker() 
{
	/// dtor
}

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker(const AliMUONQAChecker& qac) : 
    AliQACheckerBase(qac.GetName(), qac.GetTitle()) 
{
	/// copy ctor 
}   

//______________________________________________________________________________
Double_t *
AliMUONQAChecker::Check(AliQA::ALITASK_t /*index*/)
{
  /// Check data
  
  AliError(Form("This method is not implemented. Should it be ? fDataSubDir = %p (%s)",
                fDataSubDir, ( fDataSubDir ? fDataSubDir->GetPath() : "")));
  return NULL;
}

//______________________________________________________________________________
Double_t *
AliMUONQAChecker::Check(AliQA::ALITASK_t index, TObjArray ** list)
{
  /// Check objects in list
  
  if ( index == AliQA::kRAW ) 
  {
    return CheckRaws(list);
  }
  
  if ( index == AliQA::kREC)
  {
    return CheckRecPoints(list);
  }
  
  if ( index == AliQA::kESD )
  {
    return CheckESD(list);
  }
  
  AliWarning(Form("Checker for task %d not implement for the moment",index));
  return NULL;
}

//______________________________________________________________________________
TH1* 
AliMUONQAChecker::GetHisto(TObjArray* list, const char* hname) const
{
  /// Get a given histo from the list
  TH1* h = static_cast<TH1*>(list->FindObject(hname));
  if (!h)
  {
    AliError(Form("Did not find expected histo %s",hname));
  }
  return h;
}

//______________________________________________________________________________
Double_t *
AliMUONQAChecker::CheckRecPoints(TObjArray ** list)
{
  /// Check rec points
  /// Very binary check for the moment. 
  
  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 1.0 ; 
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    TH1* h = GetHisto(list[specie],"hTrackerNumberOfClustersPerDE");
  
    if ( !h ) rv[specie] =  0.75; // only a warning if histo not found, in order not to kill anything because QA failed...
  
    else if ( h->GetMean() == 0.0 ) rv[specie] =  MarkHisto(*h,0.0);
  }
  return rv;
}

//______________________________________________________________________________
Double_t 
AliMUONQAChecker::MarkHisto(TH1& histo, Double_t value) const
{
  /// Mark histo as originator of some QA error/warning
  
  if ( value != 1.0 )
  {
    histo.SetBit(AliQA::GetQABit());
  }
  
  return value;
}

//______________________________________________________________________________
Double_t *
AliMUONQAChecker::CheckESD(TObjArray ** list)
{
  /// Check ESD
  
  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 1.0 ; 
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    
    TH1* h = GetHisto(list[specie],"hESDnTracks");
  
    if (!h) rv[specie] = 0.75;
  
    else if ( h->GetMean() == 0.0 ) rv[specie] =  MarkHisto(*h,0.0); // no track -> fatal
  
    h = GetHisto(list[specie],"hESDMatchTrig");
  
    if (!h) rv[specie] =  0.75;
  
    else if (h->GetMean() == 0.0 ) rv[specie] = MarkHisto(*h,0.25); // no trigger matching -> error
  }
  return rv;
}

//______________________________________________________________________________
Double_t *
AliMUONQAChecker::CheckRaws(TObjArray ** list)
{
  /// Check raws

  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 1.0 ; 
 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    TIter next(list[specie]);
    TObject* object;
    AliMUONVTrackerData* data(0x0);
  
    while ( (object=next()) && !data )
      {
        if (object->InheritsFrom("AliMUONVTrackerData"))
          {
            data = static_cast<AliMUONVTrackerData*>(object);
          }
      }

    if ( !data ) 
      {
        AliError("Did not find TrackerData in the list !");
        return NULL;
      }
  
    AliMpManuIterator it;
    Int_t detElemId;
    Int_t manuId;
    Int_t n50(0); // number of manus with occupancy above 0.5
    Int_t n75(0); // number of manus with occupancy above 0.75
    Int_t n(0); // number of manus with some occupancy
  
    while ( it.Next(detElemId,manuId) )
      {
        Float_t occ = data->Manu(detElemId,manuId,2);
        if (occ > 0 ) ++n;
        if (occ >= 0.5 ) ++n50;
        if (occ >= 0.75 ) ++n75;    
      }

    AliInfo(Form("n %d n50 %d n75 %d",n,n50,n75));
  
    if ( n == 0 ) 
      {
        AliError("Oups. Got zero occupancy in all manus ?!");
        rv[specie] =  0.0;
      }

    if ( n75 ) 
      {
        AliError(Form("Got %d manus with occupancy above 0.75",n75));
        rv[specie] =  0.1;
      }
    
    if ( n50 ) 
      {
        AliWarning(Form("Got %d manus with occupancy above 0.5",n50));
        rv[specie] =  0.9;
      }
  }
  return rv;
}

//______________________________________________________________________________
void AliMUONQAChecker::Init(const AliQA::DETECTORINDEX_t det) 
{
  // intialises QA and QA checker settings
  AliQA::Instance(det) ; 
  Float_t * hiValue  = new Float_t[AliQA::kNBIT] ; 
  Float_t * lowValue = new Float_t[AliQA::kNBIT] ;
  lowValue[AliQA::kINFO]      = 0.999   ; 
  hiValue[AliQA::kINFO]       = 1.0 ; 
  hiValue[AliQA::kWARNING]    = 0.99 ; 
  lowValue[AliQA::kWARNING]   = 0.5 ; 
  lowValue[AliQA::kERROR]     = 0.0   ; 
  hiValue[AliQA::kERROR]      = 0.5 ; 
  lowValue[AliQA::kFATAL]     = -1.0   ; 
  hiValue[AliQA::kFATAL]      = 0.0 ; 
  SetHiLo(hiValue, lowValue) ; 
}

//______________________________________________________________________________
void 
AliMUONQAChecker::SetQA(AliQA::ALITASK_t index, const Double_t * value) const
{
  /// sets the QA according the return value of the Check

  AliQA * qa = AliQA::Instance(index);
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    qa->UnSet(AliQA::kFATAL, specie);
    qa->UnSet(AliQA::kWARNING, specie);
    qa->UnSet(AliQA::kERROR, specie);
    qa->UnSet(AliQA::kINFO, specie);

    if ( ! value ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQA::kFATAL, specie) ; 
    } else {
      if ( value[specie] >= fLowTestValue[AliQA::kFATAL] && value[specie] < fUpTestValue[AliQA::kFATAL] ) 
        qa->Set(AliQA::kFATAL, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQA::kERROR] && value[specie] <= fUpTestValue[AliQA::kERROR]  )
        qa->Set(AliQA::kERROR, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQA::kWARNING] && value[specie] <= fUpTestValue[AliQA::kWARNING]  )
        qa->Set(AliQA::kWARNING, specie) ;
      else if ( value[specie] > fLowTestValue[AliQA::kINFO] && value[specie] <= fUpTestValue[AliQA::kINFO] ) 
        qa->Set(AliQA::kINFO, specie) ; 	
    }
  }
}
