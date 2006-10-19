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

#include "AliMUONSDigitizerV2.h"

#include "AliLog.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMpDEManager.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "TObjArray.h"


///
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
///

ClassImp(AliMUONSDigitizerV2)

//_____________________________________________________________________________
AliMUONSDigitizerV2::AliMUONSDigitizerV2() 
: TTask("AliMUONSDigitizerV2","From Hits to SDigits for MUON")
{
  //
  // ctor.
  //
}

//_____________________________________________________________________________
AliMUONSDigitizerV2::~AliMUONSDigitizerV2()
{
  //
  // dtor.
  //
}

//_____________________________________________________________________________
void
AliMUONSDigitizerV2::Exec(Option_t*)
{
  //
  // Go from hits to sdigits.
  //
  // In the code below, apart from the loop itself (which look complicated
  // but is really only a loop on each hit in the input file) the main
  // work is done in AliMUONResponse::DisIntegrate method, which converts
  // a single hit in (possibly) several sdigits.
  //
  
  AliDebug(1,"");
  
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader();
  AliLoader* fLoader = runLoader->GetLoader("MUONLoader");

  fLoader->LoadHits("READ");
  
  AliMUONData muonData(fLoader,"MUON","MUON");

  AliMUON* muon = static_cast<AliMUON*>(gAlice->GetModule("MUON"));
    
  Int_t nofEvents(runLoader->GetNumberOfEvents());
  for ( Int_t iEvent = 0; iEvent < nofEvents; ++iEvent ) 
  {    
    // Loop over events.
    TObjArray tdlist;
    tdlist.SetOwner(kTRUE);
    
    AliDebug(1,Form("iEvent=%d",iEvent));
    runLoader->GetEvent(iEvent);
    TTree* treeS = fLoader->TreeS();
    AliDebug(1,Form("TreeS=%p",treeS));
    if ( !treeS )
    {
      AliDebug(1,"MakeSDigitsContainer");
      fLoader->MakeSDigitsContainer();
      treeS = fLoader->TreeS();
    }
    AliDebug(1,Form("TreeS=%p",treeS));
    muonData.MakeBranch("S");
    muonData.SetTreeAddress("S");
    
    muonData.SetTreeAddress("H");
    TTree* treeH = fLoader->TreeH();
    AliDebug(1,Form("TreeH=%p",treeH));
             
    Long64_t nofTracks = treeH->GetEntries();
    for ( Long64_t iTrack = 0; iTrack < nofTracks; ++iTrack )
    {
      // Loop over the tracks of this event.
      treeH->GetEvent(iTrack);
      TClonesArray* hits = muonData.Hits();
      Int_t nofHits = hits->GetEntriesFast();
      for ( Int_t ihit = 0; ihit < nofHits; ++ihit )
      {
        // Loop over the hits of this track.
        AliMUONHit* hit = static_cast<AliMUONHit*>(hits->At(ihit)); 
        Int_t chamberId = hit->Chamber()-1;
        AliMUONChamber& chamber = muon->Chamber(chamberId);
        AliMUONResponse* response = chamber.ResponseModel();
        
        // This is the heart of this method : the dis-integration
        TList digits;        
        response->DisIntegrate(*hit,digits);
        
        TIter next(&digits);
        AliMUONDigit* d;
        while ( ( d = (AliMUONDigit*)next() ) )
        {
          // Update some sdigit information that could not be known
          // by the DisIntegrate method
          d->SetHit(ihit);
          d->AddTrack(iTrack,d->Signal());
          tdlist.Add(d);
        }
      }
      muonData.ResetHits();
    } // end of loop on tracks within an event
    
    for ( Int_t i = 0; i <= tdlist.GetLast(); ++i )
    {
      AliMUONDigit* d = (AliMUONDigit*)tdlist[i];
      StdoutToAliDebug(1,d->Print(););
      if ( d->Signal() > 0 ) // that check would be better in the disintegrate
        // method, but to compare with old sdigitizer, it has to be there.
      {
        muonData.AddSDigit(AliMpDEManager::GetChamberId(d->DetElemId()),*d);
      }
    }
    muonData.Fill("S");
    fLoader->WriteSDigits("OVERWRITE");
    
    muonData.ResetSDigits();
    fLoader->UnloadSDigits();
  } // loop on events
  
  fLoader->UnloadHits();  
}
