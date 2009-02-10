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

#include "AliMUONQADataMakerSim.h"
#include "AliMUONHit.h"  
#include "AliMUONDigit.h"  
#include "AliMUONVHitStore.h"
#include "AliMUONVDigitStore.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TTree.h>

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMakerSim
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerSim)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerSim::AliMUONQADataMakerSim() : 
    AliQADataMakerSim(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker"),
    fHitStore(0x0),
    fDigitStore(0x0)   
{
  /// Default constructor

  AliDebug(1,"");
}

//____________________________________________________________________________ 
AliMUONQADataMakerSim::AliMUONQADataMakerSim(const AliMUONQADataMakerSim& qadm) :
    AliQADataMakerSim(),
  fHitStore(0x0),
  fDigitStore(0x0)
{
  /// Copy constructor

  AliDebug(1,"");

    if ( qadm.fHitStore ) 
    {
      fHitStore = static_cast<AliMUONVHitStore*>(qadm.fHitStore->Clone());
    }
    if ( qadm.fDigitStore ) 
    {
      fDigitStore = static_cast<AliMUONVDigitStore*>(qadm.fDigitStore->Clone());
    }
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliMUONQADataMakerSim& AliMUONQADataMakerSim::operator = (const AliMUONQADataMakerSim& qadm )
{
  /// Assignment operator

    AliDebug(1,"");

    this->~AliMUONQADataMakerSim();
    new(this) AliMUONQADataMakerSim(qadm);
    return *this;
}

//__________________________________________________________________
AliMUONQADataMakerSim::~AliMUONQADataMakerSim()
{
  /// Destructor

  AliDebug(1,"");

  delete fHitStore;
  delete fDigitStore;
}

//__________________________________________________________________
void AliMUONQADataMakerSim::InitHits() 
{
  /// Initialized hit spectra
  TH1F* h0 = new TH1F("hHitDetElem", "DetElemId distribution in Hits", 1400, 100., 1500.); 
  Add2HitsList(h0, 0);

  TH1F* h1 = new TH1F("hHitPtot", "P distribution in Hits ", 300, 0., 300.); 
  Add2HitsList(h1, 1);
  return;
} 

//__________________________________________________________________
void AliMUONQADataMakerSim::InitSDigits() 
{
  /// Initialized SDigits spectra
  TH1I* h0 = new TH1I("hSDigitsDetElem", "Detection element distribution in SDigits",  1400, 100, 1500); 
  Add2SDigitsList(h0, 0);

  TH1F* h1 = new TH1F("hSDigitsCharge", "Charge distribution in SDigits", 4096, 0, 4095); 
  Add2SDigitsList(h1, 1);

}  

//__________________________________________________________________
void AliMUONQADataMakerSim::InitDigits() 
{
  /// Initialized Digits spectra 
  TH1I* h0 = new TH1I("hDigitsDetElem", "Detection element distribution in Digits",  1400, 100, 1500); 
  Add2DigitsList(h0, 0);

  TH1I* h1 = new TH1I("hDigitsADC", "ADC distribution in Digits", 4096, 0, 4095); 
  Add2DigitsList(h1, 1);  

} 

//__________________________________________________________________
void AliMUONQADataMakerSim::MakeHits(TTree* hitsTree)        
{
  /// makes data from Hits
  if (!fHitStore)
    fHitStore = AliMUONVHitStore::Create(*hitsTree);
  fHitStore->Connect(*hitsTree, false);
  hitsTree->GetEvent(0);
    
  TIter next(fHitStore->CreateIterator());

  AliMUONHit* hit = 0x0;

  while ( ( hit = static_cast<AliMUONHit*>(next()) ) )
  {
    GetHitsData(0)->Fill(hit->DetElemId());
    GetHitsData(1)->Fill(hit->Momentum());
  }

  
}

//__________________________________________________________________
void AliMUONQADataMakerSim::MakeSDigits(TTree* sdigitsTree)        
{
  /// makes data from SDigits
  if (!fDigitStore)
    fDigitStore = AliMUONVDigitStore::Create(*sdigitsTree);
  fDigitStore->Connect(*sdigitsTree, false);
  sdigitsTree->GetEvent(0);
    
  TIter next(fDigitStore->CreateIterator());

  AliMUONVDigit* dig = 0x0;

  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) )
  {
    GetSDigitsData(0)->Fill(dig->DetElemId());
    GetSDigitsData(1)->Fill(dig->Charge());
  }
} 

//__________________________________________________________________
void AliMUONQADataMakerSim::MakeDigits(TTree* digitsTree)         
{
   /// makes data from Digits
  if (!fDigitStore)
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
  fDigitStore->Connect(*digitsTree, false);
  digitsTree->GetEvent(0);
    
  TIter next(fDigitStore->CreateIterator());

  AliMUONVDigit* dig = 0x0;

  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) )
  {
    GetDigitsData(0)->Fill(dig->DetElemId());
    GetDigitsData(1)->Fill(dig->ADC());
  }
}
      
//____________________________________________________________________________ 
void AliMUONQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray** list)
{
    ///Detector specific actions at end of cycle
    // do the QA checking
    AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;  
}


//____________________________________________________________________________ 
void AliMUONQADataMakerSim::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle
  
}
