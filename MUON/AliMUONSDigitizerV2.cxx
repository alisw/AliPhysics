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
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "TObjArray.h"

ClassImp(AliMUONSDigitizerV2)

//_____________________________________________________________________________
AliMUONSDigitizerV2::AliMUONSDigitizerV2() 
: TTask("AliMUONSDigitizerV2","From Hits to SDigits for MUON")
{
}

//_____________________________________________________________________________
AliMUONSDigitizerV2::~AliMUONSDigitizerV2()
{
}

//_____________________________________________________________________________
void
AliMUONSDigitizerV2::Exec(Option_t*)
{
  //
  // Go from hits to sdigits.
  //
  
  AliDebug(1,"");
  
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader();
  AliLoader* fLoader = runLoader->GetLoader("MUONLoader");

  fLoader->LoadHits("READ");
  
  AliMUONData muonData(fLoader,"MUON","MUON");

  AliMUON* muon = static_cast<AliMUON*>(gAlice->GetModule("MUON"));
    
  const Int_t nofEvents(runLoader->GetNumberOfEvents());
  for ( Int_t iEvent = 0; iEvent < nofEvents; ++iEvent ) 
  {    
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
      treeH->GetEvent(iTrack);
      TClonesArray* hits = muonData.Hits();
      Int_t nofHits = hits->GetEntriesFast();
      for ( Int_t ihit = 0; ihit < nofHits; ++ihit )
      {
        AliMUONHit* hit = static_cast<AliMUONHit*>(hits->At(ihit)); 
        Int_t chamberId = hit->Chamber()-1;
        AliMUONChamber& chamber = muon->Chamber(chamberId);
        AliMUONResponse* response = chamber.ResponseModel();
        TList digits;
        response->DisIntegrate(*hit,digits);
        TIter next(&digits);
        AliMUONDigit* d;
        while ( ( d = (AliMUONDigit*)next() ) )
        {
          d->SetHit(ihit);
          tdlist.Add(d);
        }
      }
      muonData.ResetHits();
    } // loop on tracks within an event
    tdlist.Sort(); // not really needed, except for comparing with old sdigitizer
    for ( Int_t i = 0; i <= tdlist.GetLast(); ++i )
    {
      AliMUONDigit* d = (AliMUONDigit*)tdlist[i];
      StdoutToAliDebug(1,d->Print(););
      if ( d->Signal() > 0 ) // that check would be better in the disintegrate
        // method, but to compare with old sdigitizer, it has to be there.
      {
        muonData.AddSDigit(d->DetElemId()/100-1,*d);
      }
    }
    muonData.Fill("S");
    fLoader->WriteSDigits("OVERWRITE");
    muonData.ResetSDigits();
    fLoader->UnloadSDigits();
  } // loop on events
  
  fLoader->UnloadHits();  
}
