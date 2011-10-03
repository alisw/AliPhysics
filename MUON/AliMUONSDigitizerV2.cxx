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

#include <TClonesArray.h>
#include <TParticle.h>
#include <TTree.h>

#include "AliMUONSDigitizerV2.h"

#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONVDigit.h"
#include "AliMUONHit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVHitStore.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONConstants.h"

#include "AliMpCDB.h"
#include "AliMpDEManager.h"

#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliHeader.h"
#include "AliGenCocktailEventHeader.h"

//-----------------------------------------------------------------------------
/// The sdigitizer performs the transformation from hits (energy deposits by
/// the transport code) to sdigits (equivalent of charges on pad).
///
/// It does so by converting the energy deposit into a charge and then spreading
/// this charge over several pads, according to the response function (a 
/// Mathieson distribution, basically).
/// 
/// See also AliMUONResponseV0, which is doing the real job (in DisIntegrate
/// method), while this sdigitizer is just "steering" the process.
///
/// Please note that we do *not* merge sdigits after creation, which means
/// that after sdigitization, a given pad can have several sdigits. This
/// merging is taken care of later on by the digitizer(V3).
//-----------------------------------------------------------------------------

ClassImp(AliMUONSDigitizerV2)

Float_t  AliMUONSDigitizerV2::fgkMaxIntTime = 10.0;
Float_t  AliMUONSDigitizerV2::fgkMaxPosTimeDif = 1.22E-6;
Float_t  AliMUONSDigitizerV2::fgkMaxNegTimeDif = -3.5E-6;
Float_t  AliMUONSDigitizerV2::fgkMinTimeDif = 25E-9;

//_____________________________________________________________________________
AliMUONSDigitizerV2::AliMUONSDigitizerV2() 
: TTask("AliMUONSDigitizerV2","From Hits to SDigits for MUON")
{
  ///
  /// ctor.
  ///

  // Load mapping
  if ( ! AliMpCDB::LoadMpSegmentation() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
}

//_____________________________________________________________________________
AliMUONSDigitizerV2::~AliMUONSDigitizerV2()
{
  ///
  /// dtor.
  ///
}

//_____________________________________________________________________________
void
AliMUONSDigitizerV2::Exec(Option_t*)
{
  ///
  /// Go from hits to sdigits.
  ///
  /// In the code below, apart from the loop itself (which look complicated
  /// but is really only a loop on each hit in the input file) the main
  /// work is done in AliMUONResponse::DisIntegrate method, which converts
  /// a single hit in (possibly) several sdigits.
  ///
  
  AliDebug(1,"");
  
  AliRunLoader* runLoader = AliRunLoader::Instance();
  AliLoader* loader = runLoader->GetDetectorLoader("MUON");

  loader->LoadHits("READ");
  
  AliMUON* muon = static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  
  Int_t nofEvents(runLoader->GetNumberOfEvents());
  
  TString classname = muon->DigitStoreClassName();
  
  AliMUONVDigitStore* sDigitStore = AliMUONVDigitStore::Create(classname.Data());
  
  if (!sDigitStore)
  {
    AliFatal(Form("Could not create digitstore of class %s",classname.Data()));
  }
  
  AliDebug(1,Form("Will use digitStore of type %s",sDigitStore->ClassName()));
          
  // average arrival time to chambers, for pileup studies

  for ( Int_t iEvent = 0; iEvent < nofEvents; ++iEvent ) 
  {    
    // Loop over events.
    TObjArray tdlist;
    tdlist.SetOwner(kTRUE);
    
    AliDebug(1,Form("iEvent=%d",iEvent));
    runLoader->GetEvent(iEvent);
  
    // for pile up studies
    float t0=fgkMaxIntTime;  int aa=0;
    AliHeader* header = runLoader->GetHeader();   
    AliGenCocktailEventHeader* cocktailHeader =
      dynamic_cast<AliGenCocktailEventHeader*>(header->GenEventHeader());
    if (cocktailHeader) {
      AliGenCocktailEventHeader* genEventHeader = (AliGenCocktailEventHeader*) (header->GenEventHeader());
      TList* headers = genEventHeader->GetHeaders();
      TIter nextH(headers);
      AliGenEventHeader *entry; 
      while((entry = (AliGenEventHeader*)nextH())) {
	float t = entry->InteractionTime();	
	if (TMath::Abs(t)<TMath::Abs(t0)) t0 = t;      
	aa++;
      }
    } else {
      AliGenEventHeader* evtHeader = 
	(AliGenEventHeader*)(header->GenEventHeader());
      float t = evtHeader->InteractionTime();		
      if (TMath::Abs(t)<TMath::Abs(t0)) t0 = t;           
      aa++;
    }

    loader->MakeSDigitsContainer();

    TTree* treeS = loader->TreeS();

    if ( !treeS )
    {
      AliFatal("");
    }

    sDigitStore->Connect(*treeS);
    
    TTree* treeH = loader->TreeH();

    AliMUONVHitStore* hitStore = AliMUONVHitStore::Create(*treeH);
    hitStore->Connect(*treeH);
    
    Long64_t nofTracks = treeH->GetEntries();
    
    for ( Long64_t iTrack = 0; iTrack < nofTracks; ++iTrack )
    {
      // Loop over the tracks of this event.
      treeH->GetEvent(iTrack);

      AliMUONHit* hit;
      TIter next(hitStore->CreateIterator());
      Int_t ihit(0);
      
      while ( ( hit = static_cast<AliMUONHit*>(next()) ) )       
      {
	Int_t chamberId = hit->Chamber()-1;
 	Float_t age = hit->Age()-t0;

        AliMUONChamber& chamber = muon->Chamber(chamberId);
        AliMUONResponse* response = chamber.ResponseModel();
        
        // This is the heart of this method : the dis-integration
        TList digits;        
	if (aa>1){  // if there are pileup events
	  Float_t chamberTime = AliMUONConstants::AverageChamberT(chamberId);
	  Float_t timeDif=age-chamberTime;	  
	  if (timeDif>fgkMaxPosTimeDif || timeDif<fgkMaxNegTimeDif) {
	    continue;
	  }
	  if(TMath::Abs(timeDif)>fgkMinTimeDif){
	    response->DisIntegrate(*hit,digits,timeDif);
	  }
	  else{
	    response->DisIntegrate(*hit,digits,0.);
	  }
 	}
 	else{
	  response->DisIntegrate(*hit,digits,0.);
	}
        
        TIter nextd(&digits);
        AliMUONVDigit* d;
        while ( ( d = (AliMUONVDigit*)nextd() ) )
        {
          // Update some sdigit information that could not be known
          // by the DisIntegrate method
          d->SetHit(ihit);
	  d->SetTime(age);
          d->AddTrack(hit->GetTrack(),d->Charge());
          tdlist.Add(d);
        }
        ++ihit;
      }
      hitStore->Clear();
    } // end of loop on tracks within an event
    
    TIter next(&tdlist);
    AliMUONVDigit* d;
    
    while ( ( d = static_cast<AliMUONVDigit*>(next()) ) )
    {
      d->ChargeInFC(kTRUE);
      
      AliMUONVDigit* added = sDigitStore->Add(*d,AliMUONVDigitStore::kMerge);
      if (!added)
      {
        AliError("Could not add digit to digitStore");
      }
    }
  
    treeS->Fill();
    
    loader->WriteSDigits("OVERWRITE");
    
    sDigitStore->Clear();
    
    loader->UnloadSDigits();
    
    delete hitStore;
    
  } // loop on events
  
  loader->UnloadHits();  
  
  delete sDigitStore;

}
