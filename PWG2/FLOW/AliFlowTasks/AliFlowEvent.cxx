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

/*****************************************************************
  AliFlowEvent: Event container for flow analysis

  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#include "Riostream.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCFManager.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowEvent.h"

ClassImp(AliFlowEvent)

//-----------------------------------------------------------------------

AliFlowEvent::AliFlowEvent():
  AliFlowEventSimple()
{
  //ctor
  cout << "AliFlowEvent: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------

AliFlowEvent::AliFlowEvent(const AliFlowEvent& event):
  AliFlowEventSimple(event)
{
  //cpy ctor
}

//-----------------------------------------------------------------------

AliFlowEvent& AliFlowEvent::operator=(const AliFlowEvent& event)
{
  //assignment operator
  AliFlowEventSimple::operator=(event);
  return *this;
}

//-----------------------------------------------------------------------
void AliFlowEvent::SetMCReactionPlaneAngle(const AliMCEvent* mcEvent)
{
  //sets the event plane angle from the proper header in the MC

  //COCKTAIL with HIJING
  if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Cocktail Header"))   //returns 0 if matches
  {
    AliGenCocktailEventHeader *headerC = dynamic_cast<AliGenCocktailEventHeader *> (mcEvent-> GenEventHeader());
    if (headerC)
    {
      TList *lhd = headerC->GetHeaders();
      if (lhd)
      {
        AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(0));
        if (hdh) AliFlowEventSimple::SetMCReactionPlaneAngle( hdh->ReactionPlaneAngle() );
      }
    }
  }
  //GEVSIM
  else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"GeVSim header"))   //returns 0 if matches
  {
    AliGenGeVSimEventHeader* headerG = dynamic_cast<AliGenGeVSimEventHeader*>(mcEvent->GenEventHeader());
    if (headerG) AliFlowEventSimple::SetMCReactionPlaneAngle( headerG->GetEventPlane() );
  }
  //HIJING
  else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Hijing"))   //returns 0 if matches
  {
    AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
    if (headerH) AliFlowEventSimple::SetMCReactionPlaneAngle( headerH->ReactionPlaneAngle() );
  }
  //EPOS
  else if (!strcmp(mcEvent->GenEventHeader()->GetName(),"EPOS"))
  {
    AliGenEposEventHeader* headerE = dynamic_cast<AliGenEposEventHeader*>(mcEvent->GenEventHeader());
    if (headerE) AliFlowEventSimple::SetMCReactionPlaneAngle( headerE->ReactionPlaneAngle() );
  }
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliMCEvent* anInput,
                            const AliCFManager* intCFManager,
                            const AliCFManager* diffCFManager):
  AliFlowEventSimple(20)
{
  //Fills the event from the MC kinematic information

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    //get input particle
    AliMCParticle* pParticle = dynamic_cast<AliMCParticle*>(anInput->GetTrack(itrkN));
    if (!pParticle) continue;

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (intCFManager && diffCFManager)
    {
      rpOK = intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle);
      poiOK = diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle);
    }
    if (!(rpOK||poiOK)) continue;

    //TODO maybe make a class AliFlowTrack with a constructor from AliVParticle
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetEta(pParticle->Eta());
    pTrack->SetPhi(pParticle->Phi());
    pTrack->SetPt(pParticle->Pt());

    if (rpOK && intCFManager)
    {
      pTrack->SetForRPSelection(kTRUE);
      fEventNSelTracksRP++;
    }
    if (poiOK && diffCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
    }

    AddTrack(pTrack) ;
  }//for all tracks
  SetMCReactionPlaneAngle(anInput);
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
                            const AliCFManager* intCFManager,
                            const AliCFManager* diffCFManager ):
  AliFlowEventSimple(20)
{
  //Fills the event from the ESD

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (intCFManager && diffCFManager)
    {
      rpOK = ( intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
               intCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
      poiOK = ( diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
                diffCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
    }
    if (!(rpOK || poiOK)) continue;

    //make new AliFLowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(pParticle->Pt() );
    pTrack->SetEta(pParticle->Eta() );
    pTrack->SetPhi(pParticle->Phi() );

    //marking the particles used for int. flow:
    if(rpOK && intCFManager)
    {
      pTrack->SetForRPSelection(kTRUE);
      fEventNSelTracksRP++;
    }
    //marking the particles used for diff. flow:
    if(poiOK && diffCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
    }

    AddTrack(pTrack);
  }//end of while (itrkN < iNumberOfInputTracks)
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliAODEvent* anInput,
                            const AliCFManager* intCFManager,
                            const AliCFManager* diffCFManager):
  AliFlowEventSimple(20)
{
  //Fills the event from the AOD
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    AliAODTrack* pParticle = anInput->GetTrack(itrkN);   //get input particle

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (intCFManager && diffCFManager)
    {
      rpOK = ( intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
               intCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
      poiOK = ( diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
                diffCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
    }
    if (!(rpOK || poiOK)) continue;

    //make new AliFlowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(pParticle->Pt() );
    pTrack->SetEta(pParticle->Eta() );
    pTrack->SetPhi(pParticle->Phi() );

    if (rpOK && intCFManager)
    {
      pTrack->SetForRPSelection(kTRUE);
      fEventNSelTracksRP++;
    }
    if (poiOK && diffCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
    }
    AddTrack(pTrack);
  }

  //  if (iSelParticlesRP >= fMinMult && iSelParticlesRP <= fMaxMult)
  //  {
  //    if ( (++fCount % 100) == 0)
  //    {
  //      if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  //      else cout<<" MC Reaction Plane Angle = unknown "<< endl;
  //      cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  //      cout<<" # of RP selected tracks = "<<iSelParticlesRP<<endl;
  //      cout<<" # of POI selected tracks = "<<iSelParticlesPOI<<endl;
  //      cout << "# " << fCount << " events processed" << endl;
  //    }
  //    return pEvent;
  //  }
  //  else
  //  {
  //    cout<<"Not enough tracks in the FlowEventSimple"<<endl;
  //    return 0;
  //  }
  //}
  //else
  //{
  //  cout<<"Event does not pass multiplicity cuts"<<endl;
  //  return 0;
  //}

}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
                            const AliMCEvent* anInputMc,
                            KineSource anOption,
                            const AliCFManager* intCFManager,
                            const AliCFManager* diffCFManager ):
  AliFlowEventSimple(20)
{
  //fills the event with tracks from the ESD and kinematics from the MC info via the track label
  if (anOption==kNoKine)
  {
    cout<<"WRONG OPTION IN AliFlowEventMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, KineSource anOption)"<<endl;
    exit(1);
  }

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  Int_t iNumberOfInputTracksMC = anInputMc->GetNumberOfTracks() ;
  if (iNumberOfInputTracksMC==-1)
  {
    cout<<"Skipping Event -- No MC information available for this event"<<endl;
    return;
  }

  //loop over ESD tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //get Label
    Int_t iLabel = pParticle->GetLabel();
    //match to mc particle
    AliMCParticle* pMcParticle = (AliMCParticle*) anInputMc->GetTrack(TMath::Abs(iLabel));

    //check
    if (TMath::Abs(pParticle->GetLabel())!=pMcParticle->Label()) cout<<"pParticle->GetLabel()!=pMcParticle->Label() "<<pParticle->GetLabel()<<"  "<<pMcParticle->Label()<<endl;

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (intCFManager && diffCFManager)
    {
      if(anOption == kESDkine)
      {
        //cout<<"take the PID from the MC & the kinematics from the ESD"<<endl;
        if (intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts1") &&
            intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle))
          rpOK=kTRUE;
        if (diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts2") &&
            diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle))
          poiOK=kTRUE;
      }
      else if (anOption == kMCkine)
      {
        //cout<<"take the PID and kinematics from the MC"<<endl;
        if (intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle))
          rpOK=kTRUE;
        if (diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle))
          poiOK=kTRUE;
      }
    }

    if (!(rpOK || poiOK)) continue;

    //make new AliFlowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    if(anOption == kESDkine)   //take the PID from the MC & the kinematics from the ESD
    {
      pTrack->SetPt(pParticle->Pt() );
      pTrack->SetEta(pParticle->Eta() );
      pTrack->SetPhi(pParticle->Phi() );
    }
    else if (anOption == kMCkine)   //take the PID and kinematics from the MC
    {
      pTrack->SetPt(pMcParticle->Pt() );
      pTrack->SetEta(pMcParticle->Eta() );
      pTrack->SetPhi(pMcParticle->Phi() );
    }
    else
    {
      cout<<"Not a valid option"<<endl;
    }

    if (rpOK && intCFManager)
    {
      fEventNSelTracksRP++;
      pTrack->SetForRPSelection();
    }
    if (poiOK && diffCFManager) pTrack->SetForPOISelection();

    AddTrack(pTrack);
  }
  SetMCReactionPlaneAngle(anInputMc);
}

