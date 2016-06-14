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
  mods:     Redmer A. Bertens (rbertens@cern.ch)
*****************************************************************/

#include "Riostream.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TArrayD.h"
#include "TProfile.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCFManager.h"
#include "AliESDtrack.h"
#include "AliESDPmdTrack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliOADBContainer.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrack.h"
#include "AliFlowVector.h"
#include "AliFlowEvent.h"
#include "AliLog.h"

using std::endl;
using std::cout;
ClassImp(AliFlowEvent)

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent():
  AliFlowEventSimple(), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
    // constructor
    for(Int_t i(0); i < 9; i++) {
        for(Int_t j(0); j < 2; j++) {
            for(Int_t k(0); k < 2; k++) {
                fMeanQ[i][j][k] = 0.; 
                fWidthQ[i][j][k] = 0.;  
                fMeanQv3[i][j][k] = 0.; 
                fWidthQv3[i][j][k] = 0.;
            }
        }
    }
    for(Int_t i(0); i < 5; i++) {
       fQxavsV0[i] = 0x0;
       fQyavsV0[i] = 0x0;  
       fQxcvsV0[i] = 0x0;
       fQycvsV0[i] = 0x0;
    }
    //ctor
  cout << "AliFlowEvent: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent(Int_t n):
  AliFlowEventSimple(n), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
    // constructor
    for(Int_t i(0); i < 9; i++) {
        for(Int_t j(0); j < 2; j++) {
            for(Int_t k(0); k < 2; k++) {
                fMeanQ[i][j][k] = 0.; 
                fWidthQ[i][j][k] = 0.;  
                fMeanQv3[i][j][k] = 0.; 
                fWidthQv3[i][j][k] = 0.;
            }
        }
    }
    for(Int_t i(0); i < 5; i++) {
       fQxavsV0[i] = 0x0;
       fQyavsV0[i] = 0x0;  
       fQxcvsV0[i] = 0x0;
       fQycvsV0[i] = 0x0;
    }
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent(const AliFlowEvent& event):
  AliFlowEventSimple(event), fApplyRecentering(event.fApplyRecentering), fDivSigma(event.fDivSigma), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // copy constructor 
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }
}

//-----------------------------------------------------------------------
AliFlowEvent& AliFlowEvent::operator=(const AliFlowEvent& event)
{
  //assignment operator
  if (&event==this) return *this;       // check self-assignment

  fApplyRecentering = event.fApplyRecentering;
  fCachedRun = event.fCachedRun;
  fVZEROcentralityBin = event.fVZEROcentralityBin;
  fDivSigma = event.fDivSigma;
  fEvent = 0x0; // should never be copied
  fChi2A = 0x0; // do not clone these; if 0x0 they will be retrieved from the rp cuts object
  fChi2C = 0x0;
  fChi3A = 0x0;
  fChi3C = 0x0;
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = event.fMeanQ[i][j][k]; 
              fWidthQ[i][j][k] = event.fWidthQ[i][j][k];  
              fMeanQv3[i][j][k] = event.fMeanQv3[i][j][k]; 
              fWidthQv3[i][j][k] = event.fWidthQv3[i][j][k];
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }

  AliFlowEventSimple::operator=(event);
  return *this;
}
//-----------------------------------------------------------------------
AliFlowTrack* AliFlowEvent::GetTrack(Int_t i)
{
  //get track i from collection
  if (i>=fNumberOfTracks) return NULL;
  AliFlowTrack* pTrack = static_cast<AliFlowTrack*>(fTrackCollection->At(i)) ;
  return pTrack;
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
  //THERMINATOR
  else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Therminator"))   //returns 0 if matches
  {
    AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
    if (headerH) AliFlowEventSimple::SetMCReactionPlaneAngle( headerH->ReactionPlaneAngle() );
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
  //AMPT
  else if (!strcmp(mcEvent-> GenEventHeader()->GetName(),"Ampt"))   //returns 0 if matches
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
  //Hydjet
  else
  {
    AliCollisionGeometry* header = dynamic_cast<AliCollisionGeometry*>(mcEvent->GenEventHeader());
    if (header) AliFlowEventSimple::SetMCReactionPlaneAngle( header->ReactionPlaneAngle() );
  }
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliMCEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }

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
    if (rpCFManager && poiCFManager)
    {
      rpOK = rpCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle);
      poiOK = poiCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle);
    }
    if (!(rpOK||poiOK)) continue;

    AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
    pTrack->SetSource(AliFlowTrack::kFromMC);

    if (rpOK && rpCFManager)
    {
      pTrack->SetForRPSelection(kTRUE);
      IncrementNumberOfPOIs(0);
    }
    if (poiOK && poiCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
      IncrementNumberOfPOIs(1);
    }

    AddTrack(pTrack) ;
  }//for all tracks
  SetMCReactionPlaneAngle(anInput);
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
   for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }
  //set run number
  if(anInput->GetRunNumber()) fRun = anInput->GetRunNumber();
 
  //Fills the event from the ESD

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (rpCFManager && poiCFManager)
    {
      rpOK = ( rpCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
               rpCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
      poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
                poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
    }
    if (!(rpOK || poiOK)) continue;

    //make new AliFLowTrack
    AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
    pTrack->SetSource(AliFlowTrack::kFromESD);

    //marking the particles used for int. flow:
    if(rpOK && rpCFManager)
    {
      pTrack->SetForRPSelection(kTRUE);
      IncrementNumberOfPOIs(0);
    }
    //marking the particles used for diff. flow:
    if(poiOK && poiCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
      IncrementNumberOfPOIs(1);
    }

    AddTrack(pTrack);
  }//end of while (itrkN < iNumberOfInputTracks)
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliAODEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }

  //set run number
  if(anInput->GetRunNumber()) fRun = anInput->GetRunNumber();
 
  //Fills the event from the AOD
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
  {
    AliAODTrack* pParticle = dynamic_cast<AliAODTrack*>(anInput->GetTrack(itrkN));
    if(!pParticle) AliFatal("Not a standard AOD");   //get input particle

    //check if pParticle passes the cuts
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
    if (rpCFManager && poiCFManager)
    {
      rpOK = ( rpCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
               rpCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
      poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
                poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
    }
    if (!(rpOK || poiOK)) continue;

    //make new AliFlowTrack
    AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
    pTrack->SetSource(AliFlowTrack::kFromAOD);

    if (rpOK /* && rpCFManager */ ) // to be fixed - with CF managers uncommented only empty events (NULL in header files)
    {
      pTrack->SetForRPSelection(kTRUE);
      IncrementNumberOfPOIs(0);
    }
    if (poiOK /* && poiCFManager*/ )
    {
      pTrack->SetForPOISelection(kTRUE);
      IncrementNumberOfPOIs(1);
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
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }
  //fills the event with tracks from the ESD and kinematics from the MC info via the track label
  if (anOption==kNoKine)
  {
    AliFatal("WRONG OPTION IN AliFlowEventMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, KineSource anOption)");
    exit(1);
  }

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  Int_t iNumberOfInputTracksMC = anInputMc->GetNumberOfTracks() ;
  if (iNumberOfInputTracksMC==-1)
  {
    AliError("Skipping Event -- No MC information available for this event");
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
    if (TMath::Abs(pParticle->GetLabel())!=pMcParticle->Label())
      AliWarning(Form("pParticle->GetLabel()!=pMcParticle->Label(), %i, %i", pParticle->GetLabel(), pMcParticle->Label()));

    //check if pParticle passes the cuts
    Bool_t rpOK = kFALSE;
    Bool_t poiOK = kFALSE;
    if (rpCFManager && poiCFManager)
    {
      if(anOption == kESDkine)
      {
        if (rpCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts1") &&
            rpCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle))
          rpOK=kTRUE;
        if (poiCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts2") &&
            poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle))
          poiOK=kTRUE;
      }
      else if (anOption == kMCkine)
      {
        if (rpCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle))
          rpOK=kTRUE;
        if (poiCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle))
          poiOK=kTRUE;
      }
    }

    if (!(rpOK || poiOK)) continue;

    //make new AliFlowTrack
    AliFlowTrack* pTrack = NULL;
    if(anOption == kESDkine)   //take the PID from the MC & the kinematics from the ESD
    {
      pTrack = new AliFlowTrack(pParticle);
    }
    else if (anOption == kMCkine)   //take the PID and kinematics from the MC
    {
      pTrack = new AliFlowTrack(pMcParticle);
    }

    if (rpOK && rpCFManager)
    {
      IncrementNumberOfPOIs(0);
      pTrack->SetForRPSelection();
    }
    if (poiOK && poiCFManager) 
    { 
      IncrementNumberOfPOIs(1);
      pTrack->SetForPOISelection();
    }

    AddTrack(pTrack);
  }
  SetMCReactionPlaneAngle(anInputMc);
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
			    const AliMultiplicity* anInputTracklets,
			    const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }

  //Select the particles of interest from the ESD
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
    {
      AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle

      //check if pParticle passes the cuts
      Bool_t poiOK = kTRUE;
      if (poiCFManager)
	{
	  poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
		    poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
	}
      if (!poiOK) continue;
      
      //make new AliFLowTrack
      AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
          
      //marking the particles used for the particle of interest (POI) selection:
      if(poiOK && poiCFManager)
	{
          IncrementNumberOfPOIs(1);
	  pTrack->SetForPOISelection(kTRUE);
	  pTrack->SetSource(AliFlowTrack::kFromESD);
	}

      AddTrack(pTrack);
    }//end of while (itrkN < iNumberOfInputTracks)

  //Select the reference particles from the SPD tracklets
  anInputTracklets = anInput->GetMultiplicity();
  Int_t multSPD = anInputTracklets->GetNumberOfTracklets();
  
  //loop over tracklets
  for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
    Float_t thetaTr= anInputTracklets->GetTheta(itracklet);
    Float_t phiTr= anInputTracklets->GetPhi(itracklet);
    // calculate eta
    Float_t etaTr = -TMath::Log(TMath::Tan(thetaTr/2.));
    
    //make new AliFLowTrackSimple
    AliFlowTrack* pTrack = new AliFlowTrack();
    pTrack->SetPt(0.0);
    pTrack->SetEta(etaTr);
    pTrack->SetPhi(phiTr);
    //marking the particles used for the reference particle (RP) selection:
    IncrementNumberOfPOIs(0);
    pTrack->SetForRPSelection(kTRUE);
    pTrack->SetSource(AliFlowTrack::kFromTracklet);

    //Add the track to the flowevent
    AddTrack(pTrack);
  }

}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* esd,
			    const AliCFManager* poiCFManager,
                            Bool_t hybrid):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }
  //Select the particles of interest from the ESD
  Int_t iNumberOfInputTracks = esd->GetNumberOfTracks() ;

  //Double_t gPt = 0.0, gP = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
//  Double_t dca3D = 0.0;       FIXME unused variable

  AliESDtrack trackTPC;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
    {

      if (!esd->GetTrack(itrkN)) continue;

      Bool_t useTPC = kFALSE;

      AliESDtrack* pParticle = esd->GetTrack(itrkN);   //get input particle

      //check if pParticle passes the cuts
      Bool_t poiOK = kTRUE;

      if (poiCFManager)
      {
        poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
                  poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
      }

      if (!(poiOK)) continue;

      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)pParticle->GetTPCInnerParam();

      if (tpcTrack)
      {

//      gPt = tpcTrack->Pt();
//      gP = tpcTrack->P();

        useTPC = kTRUE;

        const AliESDVertex *vertexSPD = esd->GetPrimaryVertexSPD();
        const AliESDVertex *vertexTPC = esd->GetPrimaryVertexTPC();

        AliExternalTrackParam copy(*tpcTrack);
        if(hybrid)
          copy.PropagateToDCA(vertexSPD,esd->GetMagneticField(),100.,dca,cov);
        else
          copy.PropagateToDCA(vertexTPC,esd->GetMagneticField(),100.,dca,cov);

//        dca3D = TMath::Sqrt(TMath::Power(dca[0],2)+TMath::Power(dca[1],2));   FIXME unused variable

      }

      //make new AliFLowTrack
      AliFlowTrack* pTrack = new AliFlowTrack(pParticle);

      pTrack->SetSource(AliFlowTrack::kFromESD);

      //marking the particles used for diff. flow:
      if(poiOK && poiCFManager)
      {
        pTrack->SetForPOISelection(kTRUE);
        IncrementNumberOfPOIs(1);
      }

      if(useTPC)
      {
        pTrack->SetForRPSelection(kTRUE);
        IncrementNumberOfPOIs(0);
      }

      AddTrack(pTrack);

    }//end of while (itrkN < iNumberOfInputTracks)

}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
			    const TH2F* anInputFMDhist,
			    const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20), fApplyRecentering(-1), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)
{
    // constructor
    for(Int_t i(0); i < 9; i++) {
        for(Int_t j(0); j < 2; j++) {
            for(Int_t k(0); k < 2; k++) {
                fMeanQ[i][j][k] = 0.; 
                fWidthQ[i][j][k] = 0.;  
                fMeanQv3[i][j][k] = 0.; 
                fWidthQv3[i][j][k] = 0.;
            }
        }
    }
    for(Int_t i(0); i < 5; i++) {
       fQxavsV0[i] = 0x0;
       fQyavsV0[i] = 0x0;  
       fQxcvsV0[i] = 0x0;
       fQycvsV0[i] = 0x0;
    }

  //Select the particles of interest from the ESD
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
    {
      AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle

      //check if pParticle passes the cuts
      Bool_t poiOK = kTRUE;
      if (poiCFManager)
	{
	  poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
		    poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
	}
      if (!poiOK) continue;
 
      //make new AliFLowTrack
      AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
          
      //marking the particles used for the particle of interest (POI) selection:
      if(poiOK && poiCFManager)
	{
    IncrementNumberOfPOIs(1);
	  pTrack->SetForPOISelection(kTRUE);
	  pTrack->SetSource(AliFlowTrack::kFromESD);
	}

      AddTrack(pTrack);
    }//end of while (itrkN < iNumberOfInputTracks)

  //Select the reference particles from the FMD hits
  //loop over FMD histogram
  Int_t iBinsEta = anInputFMDhist->GetNbinsX();
  Int_t iBinsPhi = anInputFMDhist->GetNbinsY();
  
  for (Int_t iEta = 1; iEta <= iBinsEta; iEta++){
    Double_t etaFMD = anInputFMDhist->GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= iBinsPhi; iPhi++){
      Double_t phiFMD = anInputFMDhist->GetYaxis()->GetBinCenter(iPhi);
      Double_t weightFMD = anInputFMDhist->GetBinContent(iEta,iPhi);
    
      if (weightFMD > 0.0) { //do not add empty bins
	//make new AliFLowTrackSimple
	AliFlowTrack* pTrack = new AliFlowTrack();
	pTrack->SetPt(0.0);
	pTrack->SetEta(etaFMD);
	pTrack->SetPhi(phiFMD);
	pTrack->SetWeight(weightFMD);
	//marking the particles used for the reference particle (RP) selection:
	pTrack->TagRP();
	IncrementNumberOfPOIs(0);
	pTrack->SetSource(AliFlowTrack::kFromFMD);

	//Add the track to the flowevent
	AddTrack(pTrack);
	
      }
    }
  }

}

//-----------------------------------------------------------------------
void AliFlowEvent::FindDaughters(Bool_t keepDaughtersInRPselection)
{
  //each flow track holds it's esd track index as well as its daughters esd index.
  //fill the array of daughters for every track with the pointers to flow tracks
  //to associate the mothers with daughters directly
  for (Int_t iTrack=0; iTrack<fMothersCollection->GetEntriesFast(); iTrack++)
  {
    AliFlowTrack* mother = static_cast<AliFlowTrack*>(fMothersCollection->At(iTrack));
    if (!mother) continue;
    if (mother->GetNDaughters()<1) continue;
    for (Int_t iDaughterCandidate=0; iDaughterCandidate<fNumberOfTracks; iDaughterCandidate++)
    {
      AliFlowTrack* daughterCandidate = static_cast<AliFlowTrack*>(fTrackCollection->At(iDaughterCandidate));
      Int_t esdIndexDaughterCandidate = daughterCandidate->GetID();
      for (Int_t iDaughter=0; iDaughter<mother->GetNDaughters(); iDaughter++)
      {
        Int_t esdIndexDaughter = mother->GetIDDaughter(iDaughter);
        if (esdIndexDaughter==esdIndexDaughterCandidate)
        {
          mother->SetDaughter(iDaughter,daughterCandidate);
          daughterCandidate->SetForRPSelection(keepDaughtersInRPselection);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------
void AliFlowEvent::Fill( AliFlowTrackCuts* rpCuts,
                         AliFlowTrackCuts* poiCuts )
{
  //Fills the event from a vevent: AliESDEvent,AliAODEvent,AliMCEvent
  //the input data needs to be attached to the cuts
  //we have two cases, if we're cutting the same collection of tracks
  //(same param type) then we can have tracks that are both rp and poi
  //in the other case we want to have two exclusive sets of rps and pois
  //e.g. one tracklets, the other PMD or global - USER IS RESPOSIBLE
  //FOR MAKING SURE THEY DONT OVERLAP OR ELSE THE SAME PARTICLE WILL BE
  //TAKEN TWICE

  ClearFast();
  if (!rpCuts || !poiCuts) return;
  AliFlowTrackCuts::trackParameterType sourceRP = rpCuts->GetParamType();
  AliFlowTrackCuts::trackParameterType sourcePOI = poiCuts->GetParamType();
  AliFlowTrack* pTrack=NULL;
 
  //set run
 if(rpCuts->GetRun()) fRun = rpCuts->GetRun();
 
  // if the source for rp's or poi's is the VZERO detector, get the calibration 
  // and set the calibration parameters
  if (sourceRP == AliFlowTrackCuts::kBetaVZERO) {
    SetBetaVZEROCalibrationForTrackCuts(rpCuts);
    fDivSigma = rpCuts->GetDivSigma();
    if(!rpCuts->GetApplyRecentering()) {
      // if the user does not want to recenter, switch the flag
      fApplyRecentering = -1;
    }
    // note: this flag is used in the overloaded implementation of Get2Qsub()
    // and tells the function to use as Qsub vectors the recentered Q-vectors
    // from the VZERO oadb file or from the event header
  }
  if (sourcePOI == AliFlowTrackCuts::kBetaVZERO) {
      // probably no-one will choose vzero tracks as poi's ...
      SetBetaVZEROCalibrationForTrackCuts(poiCuts);
      fDivSigma = poiCuts->GetDivSigma();
  }
 
 if (sourceRP == AliFlowTrackCuts::kDeltaVZERO) {
  SetDeltaVZEROCalibrationForTrackCuts(rpCuts);
  fDivSigma = rpCuts->GetDivSigma();
  if(!rpCuts->GetApplyRecentering()) {
   // if the user does not want to recenter, switch the flag
   fApplyRecentering = -1;
  }
 }
 if (sourcePOI == AliFlowTrackCuts::kDeltaVZERO) {
  // probably no-one will choose vzero tracks as poi's ...
  SetDeltaVZEROCalibrationForTrackCuts(poiCuts);
  fDivSigma = poiCuts->GetDivSigma();
 }

 if (sourceRP == AliFlowTrackCuts::kVZERO) {
      SetVZEROCalibrationForTrackCuts(rpCuts);
      if(!rpCuts->GetApplyRecentering()) {
          // if the user does not want to recenter, switch the flag
          fApplyRecentering = -1;
      }
      // note: this flag is used in the overloaded implementation of Get2Qsub()
      // and tells the function to use as Qsub vectors the recentered Q-vectors
      // from the VZERO oadb file or from the event header
  }
  if (sourcePOI == AliFlowTrackCuts::kVZERO) {
      // probably no-one will choose vzero tracks as poi's ...
      SetVZEROCalibrationForTrackCuts(poiCuts); 
  }

  if (sourceRP==sourcePOI)
  {
    //loop over tracks
    Int_t numberOfInputObjects = rpCuts->GetNumberOfInputObjects();
    for (Int_t i=0; i<numberOfInputObjects; i++)
    {
      //get input object (particle)
      TObject* particle = rpCuts->GetInputObject(i);

      Bool_t rp = rpCuts->IsSelected(particle,i);
      Bool_t poi = poiCuts->IsSelected(particle,i);

      if (!(rp||poi)) continue;

      //make new AliFlowTrack
      if (rp)
      {
        pTrack = rpCuts->FillFlowTrack(fTrackCollection,fNumberOfTracks);
        if (!pTrack) continue;
        pTrack->Tag(0); IncrementNumberOfPOIs(0);
        if (poi) {pTrack->Tag(1); IncrementNumberOfPOIs(1);}
        if (pTrack->GetNDaughters()>0) fMothersCollection->Add(pTrack);
      }
      else if (poi)
      {
        pTrack = poiCuts->FillFlowTrack(fTrackCollection,fNumberOfTracks);
        if (!pTrack) continue;
        pTrack->Tag(1); IncrementNumberOfPOIs(1);
        if (pTrack->GetNDaughters()>0) fMothersCollection->Add(pTrack);
      }
      fNumberOfTracks++;
    }//end of while (i < numberOfTracks)
  }
  else if (sourceRP!=sourcePOI)
  {
    //here we have two different sources of particles, so we fill
    //them independently
    //POI
    for (Int_t i=0; i<poiCuts->GetNumberOfInputObjects(); i++)
    {
      TObject* particle = poiCuts->GetInputObject(i);
      Bool_t poi = poiCuts->IsSelected(particle,i);
      if (!poi) continue;
      pTrack = poiCuts->FillFlowTrack(fTrackCollection,fNumberOfTracks);
      if (!pTrack) continue;
      pTrack->Tag(1);
      IncrementNumberOfPOIs(1);
      fNumberOfTracks++;
      if (pTrack->GetNDaughters()>0) fMothersCollection->Add(pTrack);
    }
    //RP
    Int_t numberOfInputObjects = rpCuts->GetNumberOfInputObjects();
    for (Int_t i=0; i<numberOfInputObjects; i++)
      {
      TObject* particle = rpCuts->GetInputObject(i);
      Bool_t rp = rpCuts->IsSelected(particle,i);
      if (!rp) continue;
      pTrack = rpCuts->FillFlowTrack(fTrackCollection,fNumberOfTracks);
      if (!pTrack) continue;
      pTrack->Tag(0);
      IncrementNumberOfPOIs(0);
      fNumberOfTracks++;
      if (pTrack->GetNDaughters()>0) fMothersCollection->Add(pTrack);
    }
  }
}

//-----------------------------------------------------------------------
void AliFlowEvent::InsertTrack(AliFlowTrack *track) {
  // adds a flow track at the end of the container
  AliFlowTrack *pTrack = ReuseTrack( fNumberOfTracks++ );
  *pTrack = *track;
  if (track->GetNDaughters()>0)
  {
    fMothersCollection->Add(pTrack);
  }
  return;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowEvent::ReuseTrack(Int_t i)
{
  //try to reuse an existing track, if empty, make new one
  AliFlowTrack* pTrack = static_cast<AliFlowTrack*>(fTrackCollection->At(i));
  if (pTrack)
  {
    pTrack->Clear();
  }
  else 
  {
    pTrack = new AliFlowTrack();
    fTrackCollection->AddAtAndExpand(pTrack,i);
  }
  return pTrack;
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( AliFlowTrackCuts* rpCuts,
                            AliFlowTrackCuts* poiCuts ):
  AliFlowEventSimple(20), fApplyRecentering(kFALSE), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)

{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
    fQxavsV0[i] = 0x0;
    fQyavsV0[i] = 0x0;  
    fQxcvsV0[i] = 0x0;
    fQycvsV0[i] = 0x0;
  }   
  //Fills the event from a vevent: AliESDEvent,AliAODEvent,AliMCEvent
  //the input data needs to be attached to the cuts
  Fill(rpCuts,poiCuts);
}

//-------------------------------------------------------------------//
//---- Including PMD tracks as RP --------------------------//

AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
			    const AliESDPmdTrack *pmdtracks,
			    const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20), fApplyRecentering(kFALSE), fDivSigma(kTRUE), fCachedRun(-1), fVZEROcentralityBin(-1), fEvent(0x0), fChi2A(0x0), fChi2C(0x0), fChi3A(0x0), fChi3C(0x0)

{
  // constructor
  for(Int_t i(0); i < 9; i++) {
      for(Int_t j(0); j < 2; j++) {
          for(Int_t k(0); k < 2; k++) {
              fMeanQ[i][j][k] = 0.; 
              fWidthQ[i][j][k] = 0.;  
              fMeanQv3[i][j][k] = 0.; 
              fWidthQv3[i][j][k] = 0.;
          }
      }
  }
  for(Int_t i(0); i < 5; i++) {
     fQxavsV0[i] = 0x0;
     fQyavsV0[i] = 0x0;  
     fQxcvsV0[i] = 0x0;
     fQycvsV0[i] = 0x0;
  }
  
  Float_t GetPmdEta(Float_t xPos, Float_t yPos, Float_t zPos);
  Float_t GetPmdPhi(Float_t xPos, Float_t yPos);
  //Select the particles of interest from the ESD
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  //loop over tracks
  for (Int_t itrkN=0; itrkN<iNumberOfInputTracks; itrkN++)
    {
      AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
      //check if pParticle passes the cuts
      Bool_t poiOK = kTRUE;
      if (poiCFManager)
	{
	  poiOK = ( poiCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
		    poiCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle));
	}
      if (!poiOK) continue;
      
      //make new AliFLowTrack
      AliFlowTrack* pTrack = new AliFlowTrack(pParticle);
      
      //marking the particles used for the particle of interest (POI) selection:
      if(poiOK && poiCFManager)
	{
          IncrementNumberOfPOIs(1);
	  pTrack->SetForPOISelection(kTRUE);
	  pTrack->SetSource(AliFlowTrack::kFromESD);
	}
      
      AddTrack(pTrack);
    }//end of while (itrkN < iNumberOfInputTracks)
  
  //Select the reference particles from the PMD tracks
  Int_t npmdcl = anInput->GetNumberOfPmdTracks();
  printf("======There are %d PMD tracks in this event\n-------",npmdcl);
  //loop over clusters 
  for(Int_t iclust=0; iclust < npmdcl; iclust++){
    //AliESDPmdTrack *pmdtr = anInput->GetPmdTrack(iclust);
    pmdtracks = anInput->GetPmdTrack(iclust);
    Int_t   det   = pmdtracks->GetDetector();
    //Int_t   smn   = pmdtracks->GetSmn();
    Float_t clsX  = pmdtracks->GetClusterX();
    Float_t clsY  = pmdtracks->GetClusterY();
    Float_t clsZ  = pmdtracks->GetClusterZ();
    Float_t ncell = pmdtracks->GetClusterCells();
    Float_t adc   = pmdtracks->GetClusterADC();
    //Float_t pid   = pmdtracks->GetClusterPID();
    Float_t etacls = GetPmdEta(clsX,clsY,clsZ);
    Float_t phicls = GetPmdPhi(clsX,clsY);
    //make new AliFLowTrackSimple
    AliFlowTrack* pTrack = new AliFlowTrack();
    //if(det == 0){ //selecting preshower plane only
    if(det == 0 && adc > 270 && ncell > 1){ //selecting preshower plane only
      //pTrack->SetPt(adc);//cluster adc
      pTrack->SetPt(0.0);
      pTrack->SetEta(etacls);
      pTrack->SetPhi(phicls);
      //marking the particles used for the reference particle (RP) selection:
      IncrementNumberOfPOIs(0);
      pTrack->SetForRPSelection(kTRUE);
      pTrack->SetSource(AliFlowTrack::kFromPMD);
      //Add the track to the flowevent
      AddTrack(pTrack);
    }//if det
  }
}
//----------------------------------------------------------------------------//
Float_t GetPmdEta(Float_t xPos, Float_t yPos, Float_t zPos)
{
  Float_t rpxpy, theta, eta;
  rpxpy  = TMath::Sqrt(xPos*xPos + yPos*yPos);
  theta  = TMath::ATan2(rpxpy,zPos);
  eta    = -TMath::Log(TMath::Tan(0.5*theta));
  return eta;
}
//--------------------------------------------------------------------------//
Float_t GetPmdPhi(Float_t xPos, Float_t yPos)
{
  Float_t pybypx, phi = 0., phi1;
  if(xPos==0)
    {
      if(yPos>0) phi = 90.;
      if(yPos<0) phi = 270.;
    }
  if(xPos != 0)
    {
      pybypx = yPos/xPos;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      
      if(xPos > 0 && yPos > 0) phi = phi1;        // 1st Quadrant
      if(xPos < 0 && yPos > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(xPos < 0 && yPos < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(xPos > 0 && yPos < 0) phi = 360 - phi1;  // 4th Quadrant
      
    }
  phi = phi*3.14159/180.;
  return   phi;
}
//---------------------------------------------------------------//
AliFlowVector AliFlowEvent::GetQ(Int_t n, TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights)
{
  // start with the uncalibrated Q-vector. if we're not looking at vzero data, this will be returned
  // so this is completely backward compatible
  AliFlowVector vQ = AliFlowEventSimple::GetQ(n, weightsList, usePhiWeights, usePhiWeights, useEtaWeights);
    // see if we want to re-center 2010 style
  if(fApplyRecentering == 2010) {
    AliFlowVector Qarray[2];
    AliFlowEvent::Get2Qsub(Qarray, n, weightsList, usePhiWeights, usePtWeights, useEtaWeights);
    AliFlowVector vA = Qarray[0];
    AliFlowVector vB = Qarray[1];

    // now we have the calibrated sub-event q-vectors, these will be merged using a resolution derived weight
    Double_t chiA(1.), chiC(1.), dQX(0.), dQY(0.);
    if(n==2) {
      chiA = fChi2A->At(fVZEROcentralityBin);
      chiC = fChi2C->At(fVZEROcentralityBin);
    } else if (n==3) {
      chiA = fChi3A->At(fVZEROcentralityBin);
      chiC = fChi3C->At(fVZEROcentralityBin);
    }

    // combine the vzera and vzeroc signal, note that vzeroa is stored in vB and vzeroc in vA
    dQX = chiA*chiA*vB.X()+chiC*chiC*vA.X();
    dQY = chiA*chiA*vB.Y()+chiC*chiC*vA.Y();
    vQ.Set(dQX, dQY);

  } else if (fApplyRecentering == 2011) { // 11h style recentering
    Double_t dQX = 0.;
    Double_t dQY = 0.;

    if(fEvent && fEvent->GetEventplane()) {
      fEvent->GetEventplane()->CalculateVZEROEventPlane(fEvent, 10, n, dQX, dQY);
      // set the new q-vectors (which in this case means REPLACING) 
      vQ.Set(dQX, dQY);
    } // if for some reason the info from the event header is not available, the AliFlowTrackCuts object
      // has provided the equalized multiplicity info so this should still be relatively safe
  }

  // return the Q-vector 
  return vQ;
}   
 //---------------------------------------------------------------//
void AliFlowEvent::Get2Qsub(AliFlowVector* Qarray, Int_t n, TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights)
{
  // get q vectors for the subevents. if no recentering is necessary, get the guy from the flow event simple
  AliFlowEventSimple::Get2Qsub(Qarray, n, weightsList, usePhiWeights, usePtWeights, useEtaWeights);
  AliFlowVector vA = Qarray[0];
  AliFlowVector vB = Qarray[1];
 
  // else get the recentering from the cached info
  if (fApplyRecentering == 2010)        // 10h style recentering, implemented for n=2 and n=3
  {     
    // first retrieve the q-vectors from the AliFlowEventSimple:: routine
    // extract the information form the current flow vectors
    Double_t Qxc(vA.X());       // IMPORTANT: user is responsible for the sign of eta
    Double_t Qyc(vA.Y());       // vzeroC has negative pseudorapidity and is taken as subevent A
    Double_t Qxa(vB.X());       // vzeroA has positive pseudorapidity and is taken as subevent B
    Double_t Qya(vB.Y());
    // init some values for the corrections
    
    // default values for vector a (VZEROA)
    Double_t Qxamean(0);
    Double_t Qxarms(1);
    Double_t Qyamean(0);
    Double_t Qyarms(1);
    // default values for vector b (VZEROC)
    Double_t Qxcmean(0);
    Double_t Qxcrms(1);
    Double_t Qycmean(0);
    Double_t Qycrms(1);	
    // note that defaults are chosen such that for n!=2||n!=3 re-centering is a null-operation
    
    if( n == 2) {       // second order symmetry
        Qxamean = fMeanQ[fVZEROcentralityBin][1][0];
        Qxarms  = fWidthQ[fVZEROcentralityBin][1][0];
        Qyamean = fMeanQ[fVZEROcentralityBin][1][1];
        Qyarms  = fWidthQ[fVZEROcentralityBin][1][1];

        Qxcmean = fMeanQ[fVZEROcentralityBin][0][0];
        Qxcrms  = fWidthQ[fVZEROcentralityBin][0][0];
        Qycmean = fMeanQ[fVZEROcentralityBin][0][1];
        Qycrms  = fWidthQ[fVZEROcentralityBin][0][1];	
    } else if (n == 3) {        // third order symmetry
        Qxamean = fMeanQv3[fVZEROcentralityBin][1][0];
        Qxarms  = fWidthQv3[fVZEROcentralityBin][1][0];
        Qyamean = fMeanQv3[fVZEROcentralityBin][1][1];
        Qyarms  = fWidthQv3[fVZEROcentralityBin][1][1];
  
        Qxcmean = fMeanQv3[fVZEROcentralityBin][0][0];
        Qxcrms  = fWidthQv3[fVZEROcentralityBin][0][0];
        Qycmean = fMeanQv3[fVZEROcentralityBin][0][1];
        Qycrms  = fWidthQv3[fVZEROcentralityBin][0][1];	
    }
    // do the correction    
    Double_t QxaCor = (Qxa - Qxamean)/Qxarms;
    Double_t QyaCor = (Qya - Qyamean)/Qyarms;
    Double_t QxcCor = (Qxc - Qxcmean)/Qxcrms;
    Double_t QycCor = (Qyc - Qycmean)/Qycrms;
    // update the vector
    vA.Set(QxcCor, QycCor);
    vB.Set(QxaCor, QyaCor);
  } else if (fApplyRecentering == 2011) { // 11h style recentering

    Double_t QxaCor = 0.;
    Double_t QyaCor = 0.;
    Double_t QxcCor = 0.;
    Double_t QycCor = 0.;

    if(fEvent && fEvent->GetEventplane()) {
      fEvent->GetEventplane()->CalculateVZEROEventPlane(fEvent, 8, n, QxaCor, QyaCor);
      fEvent->GetEventplane()->CalculateVZEROEventPlane(fEvent, 9, n, QxcCor, QycCor);
 
      // set the new q-vectors (which in this case means REPLACING) 
      vA.Set(QxcCor, QycCor);
      vB.Set(QxaCor, QyaCor);
    } // if for some reason the info from the event header is not available, the AliFlowTrackCuts object
      // has provided the equalized multiplicity info so this should still be relatively safe
  } else if (fApplyRecentering == 999) {
      // experimental VZERO recentering
    // first retrieve the q-vectors from the AliFlowEventSimple:: routine
    // extract the information form the current flow vectors
    Double_t Qxc(vA.X());       // IMPORTANT: user is responsible for the sign of eta
    Double_t Qyc(vA.Y());       // vzeroC has negative pseudorapidity and is taken as subevent A
    Double_t Qxa(vB.X());       // vzeroA has positive pseudorapidity and is taken as subevent B
    Double_t Qya(vB.Y());
    // init some values for the corrections
    
    Double_t Cen = fEvent->GetCentrality()->GetCentralityPercentile("V0M");
    // default values for vector a (VZEROA)
    Double_t Qxamean(fQxavsV0[n-1]->GetBinContent(fQxavsV0[n-1]->FindBin(Cen)));
    Double_t Qxarms(fQxavsV0[n-1]->GetBinError(fQxavsV0[n-1]->FindBin(Cen)));
    Double_t Qyamean(fQyavsV0[n-1]->GetBinContent(fQyavsV0[n-1]->FindBin(Cen)));
    Double_t Qyarms(fQyavsV0[n-1]->GetBinError(fQyavsV0[n-1]->FindBin(Cen)));
    // default values for vector b (VZEROC)
    Double_t Qxcmean(fQxcvsV0[n-1]->GetBinContent(fQxcvsV0[n-1]->FindBin(Cen)));
    Double_t Qxcrms(fQxcvsV0[n-1]->GetBinError(fQxcvsV0[n-1]->FindBin(Cen)));
    Double_t Qycmean(fQycvsV0[n-1]->GetBinContent(fQycvsV0[n-1]->FindBin(Cen)));
    Double_t Qycrms(fQycvsV0[n-1]->GetBinError(fQycvsV0[n-1]->FindBin(Cen)));
   
    // just a precaution ...
    if(n > 5) {
      Qxamean = 0;
      Qxarms = 1;
      Qyamean = 0;
      Qyarms = 1;
      Qxcmean = 0;
      Qxcrms = 1; 
      Qycmean = 0;    // this effectively disables the calibration
      Qycrms = 1;
    }
   
   Double_t QxaR = Qxa - Qxamean;
   Double_t QyaR = Qya - Qyamean;
   Double_t QxcR = Qxc - Qxcmean;
   Double_t QycR = Qyc - Qycmean;
   if(fDivSigma && Qxarms>0. && Qyarms>0. && Qxcrms>0. && Qycrms>0.) {
    QxaR /= Qxarms;
    QyaR /= Qyarms;
    QxcR /= Qxcrms;
    QycR /= Qycrms;
   }
   // update the vector
   vA.Set(QxcR, QycR);
   vB.Set(QxaR, QyaR);
   
  } else if (fApplyRecentering == 666) {
   // experimental VZERO recentering for 11h Full TPC Flow data, harmonics 1,2,3
   // first retrieve the q-vectors from the AliFlowEventSimple:: routine
   // extract the information form the current flow vectors
   if(n < 4) {
    Double_t Qxc(vA.X());       // IMPORTANT: user is responsible for the sign of eta
    Double_t Qyc(vA.Y());       // vzeroC has negative pseudorapidity and is taken as subevent A
    Double_t Qxa(vB.X());       // vzeroA has positive pseudorapidity and is taken as subevent B
    Double_t Qya(vB.Y());
    Double_t MultC = vA.GetMult();
    Double_t MultA = vB.GetMult();
    // init some values for the corrections
    Double_t Cen = - 999;
    //Added by Bernhard Hohlweger - bhohlweg@cern.ch
    if(fEvent->GetRunNumber() < 209122){
    	//For Run1 Data the Old Centrality Percentile Method is available whereas for Run2 a new method was implemented
    	//Cut was done for the first run of the LHC15a period
    	Cen = fEvent->GetCentrality()->GetCentralityPercentile("V0M");
    }else{
    	AliMultSelection *MultSelection = 0x0;
    	MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
    	if( !MultSelection) {
    		//If you get this warning (and EventCentrality -999) please check that the AliMultSelectionTask actually ran (before your task)
    		AliWarning("AliMultSelection not found, did you Run AliMultSelectionTask? \n");
    	}else{
    		Cen = MultSelection->GetMultiplicityPercentile("V0M");
    	}
    }
    //14062016 Bernhard Hohlweger
    // default values for vector a (VZEROA)
    Double_t Qxamean(fQxavsV0[n-1]->GetBinContent(fQxavsV0[n-1]->FindBin(Cen)));
    Double_t Qxarms(fQxavsV0[n-1]->GetBinError(fQxavsV0[n-1]->FindBin(Cen)));
    Double_t Qyamean(fQyavsV0[n-1]->GetBinContent(fQyavsV0[n-1]->FindBin(Cen)));
    Double_t Qyarms(fQyavsV0[n-1]->GetBinError(fQyavsV0[n-1]->FindBin(Cen)));
    // default values for vector b (VZEROC)
    Double_t Qxcmean(fQxcvsV0[n-1]->GetBinContent(fQxcvsV0[n-1]->FindBin(Cen)));
    Double_t Qxcrms(fQxcvsV0[n-1]->GetBinError(fQxcvsV0[n-1]->FindBin(Cen)));
    Double_t Qycmean(fQycvsV0[n-1]->GetBinContent(fQycvsV0[n-1]->FindBin(Cen)));
    Double_t Qycrms(fQycvsV0[n-1]->GetBinError(fQycvsV0[n-1]->FindBin(Cen)));
    
    Double_t QxaR = Qxa - Qxamean*MultA;
    Double_t QyaR = Qya - Qyamean*MultA;
    Double_t QxcR = Qxc - Qxcmean*MultC;
    Double_t QycR = Qyc - Qycmean*MultC;
    if(fDivSigma && Qxarms>0. && Qyarms>0. && Qxcrms>0. && Qycrms>0.) {
     QxaR /= Qxarms;
     QyaR /= Qyarms;
     QxcR /= Qxcrms;
     QycR /= Qycrms;
    }
    // update the vector
    vA.Set(QxcR, QycR);
    vB.Set(QxaR, QyaR);
    
   } else {
    cout << " WARNING: recentering not possible for harmonic " << n << " (delta calibration) " << endl;
   } // end of if(n < 4)
  }
 
 Qarray[0] = vA;
 Qarray[1] = vB;
}
//_____________________________________________________________________________
void AliFlowEvent::SetVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts) {
    // open calibration info, copied from AliAnalyisTaskVnV0.cxx
    fEvent = cuts->GetEvent();
    if(!fEvent) return; // coverity. we need to know the event to get the runnumber and centrlaity
    // get the vzero centrality percentile (cc dependent calibration)
    Float_t v0Centr(fEvent->GetCentrality()->GetCentralityPercentile("V0M"));
    if(v0Centr < 5) fVZEROcentralityBin = 0;
    else if(v0Centr < 10) fVZEROcentralityBin = 1;
    else if(v0Centr < 20) fVZEROcentralityBin = 2;
    else if(v0Centr < 30) fVZEROcentralityBin = 3;
    else if(v0Centr < 40) fVZEROcentralityBin = 4;
    else if(v0Centr < 50) fVZEROcentralityBin = 5;
    else if(v0Centr < 60) fVZEROcentralityBin = 6;
    else if(v0Centr < 70) fVZEROcentralityBin = 7;
    else fVZEROcentralityBin = 8;

    // if this event is from the same run as the previous event
    // we can use the cached calibration values, no need to re-open the 
    // aodb file, else cache the new run
    Int_t run(fEvent->GetRunNumber());
    if(fCachedRun == run) return;
    else fCachedRun = run;
    
    // relevant for 2010 data: check if the proper chi weights for merging vzero a and vzero c ep are present
    // if not, use sane defaults. centrality binning is equal to that given in the fVZEROcentralityBin snippet
    //
    // chi values can be calculated using the static helper function 
    // AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res) where res is the event plane
    // resolution in a given centrality bin
    //
    // the resolutions that were used for these defaults are
    // Double_t R2VZEROA[] = {.35, .40, .48, .50, .48, .45, .38, .26, .16};
    // Double_t R2VZEROC[] = {.45, .60, .70, .73, .68, .60, .40, .36, .17};
    //
    // Double_t R3VZEROA[] = {.22, .23, .22, .19, .15, .12, .08, .00, .00};
    // Double_t R3VZEROC[] = {.30, .30, .28, .25, .22, .17, .11, .00, .00};
    // this might need a bit of updating as they were read 'by-eye' from a performance plot ..
    Double_t chiC2[] = {0.771423, 1.10236, 1.38116, 1.48077, 1.31964, 1.10236, 0.674622, 0.600403, 0.273865};
    Double_t chiA2[] = {0.582214, 0.674622, 0.832214, 0.873962, 0.832214, 0.771423, 0.637146, 0.424255, 0.257385};
    Double_t chiC3[] = {0.493347, 0.493347, 0.458557, 0.407166, 0.356628, 0.273865, 0.176208, 6.10352e-05, 6.10352e-05};
    Double_t chiA3[] = {0.356628, 0.373474, 0.356628, 0.306702, 0.24115, 0.192322, 0.127869, 6.10352e-05, 6.10352e-05};

    // this may seem redundant but in this way the cuts object is owner of the arrays
    // even if they're created here (so we won't get into trouble with dtor, assigmnet and copying) 
    if(!cuts->GetChi2A()) cuts->SetChi2A(new TArrayD(9, chiA2));
    if(!cuts->GetChi2C()) cuts->SetChi2C(new TArrayD(9, chiC2));
    if(!cuts->GetChi3A()) cuts->SetChi3A(new TArrayD(9, chiA3));
    if(!cuts->GetChi3C()) cuts->SetChi3C(new TArrayD(9, chiC3));

    if(!fChi2A) fChi2A = cuts->GetChi2A();
    if(!fChi2C) fChi2C = cuts->GetChi2C();
    if(!fChi3A) fChi3A = cuts->GetChi3A();
    if(!fChi3C) fChi3C = cuts->GetChi3C();
 
   TFile *foadb = TFile::Open("$ALICE_PHYSICS/OADB/PWGCF/VZERO/VZEROcalibEP.root");
    if(!foadb){
	printf("OADB file $ALICE_PHYSICS/OADB/PWGCF/VZERO/VZEROcalibEP.root cannot be opened, CALIBRATION FAILED !");
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr");
    if(!cont){
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }
    if(!(cont->GetObject(run))){
        // if the multiplicity correction cannot be found for the specified run, 
        // loop over the 11h runs to see if it's 11h data
        Int_t runs11h[] = {170593, 170572, 170556, 170552, 170546, 170390, 170389, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170267, 170264, 170230, 170228, 170208, 170207, 170205, 170204, 170203, 170195, 170193, 170163, 170162, 170159, 170155, 170152, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170038, 170036, 170027, 169981, 169975, 169969, 169965, 169961, 169956, 169926, 169924, 169923, 169922, 169919, 169918, 169914, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169683, 169628, 169591, 169590, 169588, 169587, 169586, 169584, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169236, 169167, 169160, 169156, 169148, 169145, 169144, 169143, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168984, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168461, 168460, 168458, 168362, 168361, 168356, 168342, 168341, 168325, 168322, 168318, 168311, 168310, 168213, 168212, 168208, 168207, 168206, 168205, 168204, 168203, 168181, 168177, 168175, 168173, 168172, 168171, 168115, 168108, 168107, 168105, 168104, 168103, 168076, 168069, 168068, 168066, 167988, 167987, 167986, 167985, 167921, 167920, 167915, 167909, 167903, 167902, 167818, 167814, 167813, 167808, 167807, 167806, 167713, 167712, 167711, 167706, 167693};
        for(Int_t r(0); r < 176; r++) {
            if(run == runs11h[r]) {
                printf(" > run has been identified as 11h < \n");
                 if(cuts->GetVZEROgainEqualizationPerRing()) {
                    // enable or disable rings through the weights, weight 1. is enabled, 0. is disabled
                    // start with the vzero c rings (segments 0 through 31)
                    (cuts->GetUseVZERORing(0)) ? cuts->SetVZEROCpol(0, 1.) : cuts->SetVZEROCpol(0, 0.);
                    (cuts->GetUseVZERORing(1)) ? cuts->SetVZEROCpol(1, 1.) : cuts->SetVZEROCpol(1, 0.);
                    (cuts->GetUseVZERORing(2)) ? cuts->SetVZEROCpol(2, 1.) : cuts->SetVZEROCpol(2, 0.);
                    (cuts->GetUseVZERORing(3)) ? cuts->SetVZEROCpol(3, 1.) : cuts->SetVZEROCpol(3, 0.);
                    // same for vzero a
                    (cuts->GetUseVZERORing(4)) ? cuts->SetVZEROApol(0, 1.) : cuts->SetVZEROApol(0, 0.);
                    (cuts->GetUseVZERORing(5)) ? cuts->SetVZEROApol(1, 1.) : cuts->SetVZEROApol(1, 0.);
                    (cuts->GetUseVZERORing(6)) ? cuts->SetVZEROApol(2, 1.) : cuts->SetVZEROApol(2, 0.);
                    (cuts->GetUseVZERORing(7)) ? cuts->SetVZEROApol(3, 1.) : cuts->SetVZEROApol(3, 0.);
                } else {
                    // else enable all rings, which is also default
                    for(Int_t i(0); i < 4; i++) cuts->SetVZEROCpol(i, 1.);
                    for(Int_t i(0); i < 4; i++) cuts->SetVZEROApol(i, 1.);
                }
                // pass a NULL pointer to the track cuts object, the NULL pointer will identify 11h runs
                cuts->SetVZEROgainEqualisation(NULL);
                fApplyRecentering = 2011;
                return; // the rest of the steps are not necessary
            }
        }
        // the run has not been identified as lhc11h data, so we assume a template calibration
	printf("OADB object hMultVZEROBefCorr is not available for run %i (used default run 137366)\n",run);
	run = 137366;
    }
    printf(" > run has been identified as 10h < \n");
    // step 1) get the proper multiplicity weights from the vzero signal
    TProfile* fMultVZERO = ((TH2F *) cont->GetObject(run))->ProfileX();

    TF1 *fpol0 = new TF1("fpol0","pol0"); 
    if(cuts->GetVZEROgainEqualizationPerRing()) {
        // do the calibration per ring
        // start with the vzero c rings (segments 0 through 31)
        fMultVZERO->Fit(fpol0, "N0", "", 0, 8);
        (cuts->GetUseVZERORing(0)) ? cuts->SetVZEROCpol(0, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(0, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 8, 16);
        (cuts->GetUseVZERORing(1)) ? cuts->SetVZEROCpol(1, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(1, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 16, 24);
        (cuts->GetUseVZERORing(2)) ? cuts->SetVZEROCpol(2, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(2, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 24, 32);
        (cuts->GetUseVZERORing(3)) ? cuts->SetVZEROCpol(3, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(3, 0.);
        // same thing for vero A
        fMultVZERO->Fit(fpol0, "N0", "", 32, 40);
        (cuts->GetUseVZERORing(4)) ? cuts->SetVZEROApol(0, fpol0->GetParameter(0)) : cuts->SetVZEROApol(0, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 40, 48);
        (cuts->GetUseVZERORing(5)) ? cuts->SetVZEROApol(1, fpol0->GetParameter(0)) : cuts->SetVZEROApol(1, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 48, 56);
        (cuts->GetUseVZERORing(6)) ? cuts->SetVZEROApol(2, fpol0->GetParameter(0)) : cuts->SetVZEROApol(2, 0.);
        fMultVZERO->Fit(fpol0, "N0", "", 56, 64);
        (cuts->GetUseVZERORing(7)) ? cuts->SetVZEROApol(3, fpol0->GetParameter(0)) : cuts->SetVZEROApol(3, 0.);
    } else {
        // do the calibration in one go. the calibration will still be 
        // stored per ring, but each ring has the same weight now
       fMultVZERO->Fit(fpol0,"N0","",0,31);
       for(Int_t i(0); i < 4; i++) cuts->SetVZEROCpol(i, fpol0->GetParameter(0));
       fMultVZERO->Fit(fpol0,"N0","",32,64);
       for(Int_t i(0); i < 4; i++) cuts->SetVZEROApol(i, fpol0->GetParameter(0));
    }
    // the parameters to weigh the vzero track cuts have been extracted now, 
    // so we can pass them to the current track cuts obect
    cuts->SetVZEROgainEqualisation(fMultVZERO);       // passed as a TH1

    // step 2) reweight the q-vectors that will be  called by flow methods which use
    // subevents
    // underlying assumption is that subevent a uses VZEROA
    // and subevent b uses VZEROC
    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < 9;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
	
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}

                // after grabbing all the info, set the CORRECTION TERMS to
                // the 2nd and 3rd order qsub-vectors
                // we do this here for all centralities, so that subsequent events
                // can grab the correction from these cached values
                fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

		//for v3
		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya3_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

     	    }
	}
    }
    // set the recentering style (might be switched back to -1 if recentering is disabeled)
    fApplyRecentering = 2010;
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetBetaVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts) {
    // implementation of beta vzero calibration

    fEvent = cuts->GetEvent();
    if(!fEvent) return; // coverity. we need to know the event to get the runnumber and centrlaity
    // get the vzero centrality percentile (cc dependent calibration)
    

    // if this event is from the same run as the previous event
    // we can use the cached calibration values, no need to re-open the 
    // aodb file, else cache the new run
    Int_t run(fEvent->GetRunNumber());
    if(fCachedRun == run) return;
    else fCachedRun = run;
    
    // check if the proper chi weights for merging vzero a and vzero c ep are present
    // if not, use sane defaults. centrality binning is equal to that given in the fVZEROcentralityBin snippet
    //
    // chi values can be calculated using the static helper function 
    // AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res) where res is the event plane
    // resolution in a given centrality bin
    //
    // the resolutions that were used for these defaults are
    // Double_t R2VZEROA[] = {.35, .40, .48, .50, .48, .45, .38, .26, .16};
    // Double_t R2VZEROC[] = {.45, .60, .70, .73, .68, .60, .40, .36, .17};
    //
    // Double_t R3VZEROA[] = {.22, .23, .22, .19, .15, .12, .08, .00, .00};
    // Double_t R3VZEROC[] = {.30, .30, .28, .25, .22, .17, .11, .00, .00};
    Double_t chiC2[] = {0.771423, 1.10236, 1.38116, 1.48077, 1.31964, 1.10236, 0.674622, 0.600403, 0.273865};
    Double_t chiA2[] = {0.582214, 0.674622, 0.832214, 0.873962, 0.832214, 0.771423, 0.637146, 0.424255, 0.257385};
    Double_t chiC3[] = {0.493347, 0.493347, 0.458557, 0.407166, 0.356628, 0.273865, 0.176208, 6.10352e-05, 6.10352e-05};
    Double_t chiA3[] = {0.356628, 0.373474, 0.356628, 0.306702, 0.24115, 0.192322, 0.127869, 6.10352e-05, 6.10352e-05};

    // this may seem redundant but in this way the cuts object is owner of the arrays
    // even if they're created here (so we won't get into trouble with dtor, assigmnet and copying) 
    if(!cuts->GetChi2A()) cuts->SetChi2A(new TArrayD(9, chiA2));
    if(!cuts->GetChi2C()) cuts->SetChi2C(new TArrayD(9, chiC2));
    if(!cuts->GetChi3A()) cuts->SetChi3A(new TArrayD(9, chiA3));
    if(!cuts->GetChi3C()) cuts->SetChi3C(new TArrayD(9, chiC3));

    if(!fChi2A) fChi2A = cuts->GetChi2A();
    if(!fChi2C) fChi2C = cuts->GetChi2C();
    if(!fChi3A) fChi3A = cuts->GetChi3A();
    if(!fChi3C) fChi3C = cuts->GetChi3C();
 

    TFile *foadb = TFile::Open("$ALICE_PHYSICS/PWGCF/FLOW/database/calibV0_filtered.root");
    if(!foadb){
	printf("OADB file $ALICE_PHYSICS/PWGCF/FLOW/database/calibV0_filtered.root cannot be opened, CALIBRATION FAILED !");
	return;
    }

    // first get the mutiplicity of the vzero channels before gain equalization
    AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr_filtered");
    if(!cont){
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }
    if(!(cont->GetObject(run))) {
        // if the multiplicity correction cannot be found for the specified run, 
        // loop over the 11h runs to see if it's 11h data
        Int_t runs11h[] = {170593, 170572, 170556, 170552, 170546, 170390, 170389, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170267, 170264, 170230, 170228, 170208, 170207, 170205, 170204, 170203, 170195, 170193, 170163, 170162, 170159, 170155, 170152, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170038, 170036, 170027, 169981, 169975, 169969, 169965, 169961, 169956, 169926, 169924, 169923, 169922, 169919, 169918, 169914, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169683, 169628, 169591, 169590, 169588, 169587, 169586, 169584, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169236, 169167, 169160, 169156, 169148, 169145, 169144, 169143, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168984, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168461, 168460, 168458, 168362, 168361, 168356, 168342, 168341, 168325, 168322, 168318, 168311, 168310, 168213, 168212, 168208, 168207, 168206, 168205, 168204, 168203, 168181, 168177, 168175, 168173, 168172, 168171, 168115, 168108, 168107, 168105, 168104, 168103, 168076, 168069, 168068, 168066, 167988, 167987, 167986, 167985, 167921, 167920, 167915, 167909, 167903, 167902, 167818, 167814, 167813, 167808, 167807, 167806, 167713, 167712, 167711, 167706, 167693};
        for(Int_t r(0); r < 176; r++) {
            if(run == runs11h[r]) {
                printf(" > run has been identified as 11h < \n");
                if(cuts->GetVZEROgainEqualizationPerRing()) {
                    // enable or disable rings through the weights, weight 1. is enabled, 0. is disabled
                    // start with the vzero c rings (segments 0 through 31)
                    (cuts->GetUseVZERORing(0)) ? cuts->SetVZEROCpol(0, 1.) : cuts->SetVZEROCpol(0, 0.);
                    (cuts->GetUseVZERORing(1)) ? cuts->SetVZEROCpol(1, 1.) : cuts->SetVZEROCpol(1, 0.);
                    (cuts->GetUseVZERORing(2)) ? cuts->SetVZEROCpol(2, 1.) : cuts->SetVZEROCpol(2, 0.);
                    (cuts->GetUseVZERORing(3)) ? cuts->SetVZEROCpol(3, 1.) : cuts->SetVZEROCpol(3, 0.);
                    // same for vzero a
                    (cuts->GetUseVZERORing(4)) ? cuts->SetVZEROApol(0, 1.) : cuts->SetVZEROApol(0, 0.);
                    (cuts->GetUseVZERORing(5)) ? cuts->SetVZEROApol(1, 1.) : cuts->SetVZEROApol(1, 0.);
                    (cuts->GetUseVZERORing(6)) ? cuts->SetVZEROApol(2, 1.) : cuts->SetVZEROApol(2, 0.);
                    (cuts->GetUseVZERORing(7)) ? cuts->SetVZEROApol(3, 1.) : cuts->SetVZEROApol(3, 0.);
                } else {
                    // else enable all rings, which is also default
                    for(Int_t i(0); i < 4; i++) cuts->SetVZEROCpol(i, 1.);
                    for(Int_t i(0); i < 4; i++) cuts->SetVZEROApol(i, 1.);
                }
                // pass a NULL pointer to the track cuts object, the NULL pointer will identify 11h runs
                cuts->SetVZEROgainEqualisation(NULL);
                fApplyRecentering = 2011;
                return; // the rest of the steps are not necessary
            }
        }
        // the run has not been identified as lhc11h data, so we assume a template calibration
	printf("OADB object hMultVZEROBefCorr is not available for run %i (used default run 138275)\n",run);
	run = 138275;
    }
    printf(" > run has been identified as 10h < \n");
    // step 0) get the profile which contains average multiplicity per VZERO channel
    TProfile* fMultVZERO = static_cast<TProfile*>(cont->GetObject(run));

    TF1 *fpol0 = new TF1("fpol0","pol0"); 
    // step 1) extract the proper weights from the profile. Q-vector recentering relies on 
    // ring-by-ring gain equalization, so the terms are extracted ring-by-ring here

    // start with the vzero c rings (segments 0 through 31)
    fMultVZERO->Fit(fpol0, "N0", "", 0, 8);
    (cuts->GetUseVZERORing(0)) ? cuts->SetVZEROCpol(0, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(0, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 8, 16);
    (cuts->GetUseVZERORing(1)) ? cuts->SetVZEROCpol(1, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(1, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 16, 24);
    (cuts->GetUseVZERORing(2)) ? cuts->SetVZEROCpol(2, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(2, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 24, 32);
    (cuts->GetUseVZERORing(3)) ? cuts->SetVZEROCpol(3, fpol0->GetParameter(0)) : cuts->SetVZEROCpol(3, 0.);
    // same thing for vero A
    fMultVZERO->Fit(fpol0, "N0", "", 32, 40);
    (cuts->GetUseVZERORing(4)) ? cuts->SetVZEROApol(0, fpol0->GetParameter(0)) : cuts->SetVZEROApol(0, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 40, 48);
    (cuts->GetUseVZERORing(5)) ? cuts->SetVZEROApol(1, fpol0->GetParameter(0)) : cuts->SetVZEROApol(1, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 48, 56);
    (cuts->GetUseVZERORing(6)) ? cuts->SetVZEROApol(2, fpol0->GetParameter(0)) : cuts->SetVZEROApol(2, 0.);
    fMultVZERO->Fit(fpol0, "N0", "", 56, 64);
    (cuts->GetUseVZERORing(7)) ? cuts->SetVZEROApol(3, fpol0->GetParameter(0)) : cuts->SetVZEROApol(3, 0.);

    // the parameters to weigh the vzero track cuts have been extracted now, 
    // so we can pass them to the current track cuts obect
    cuts->SetVZEROgainEqualisation(fMultVZERO);       // passed as a TH1
 
    // step 2) extract the calibration histograms from the database and
    // pass them to the cuts object
    //
    // first index of the oadb array is the harmonic n, the second index is either qax, qay, qcx, qcy
    AliOADBContainer* h[5][4];
    for(Int_t i(0); i < 5; i++) {
      h[i][0] = (AliOADBContainer*)foadb->Get(Form("hQxa%i_filtered", i+1));
      if(h[i][0]) fQxavsV0[i] = static_cast<TH1F*>(h[i][0]->GetObject(run));
      h[i][1] = (AliOADBContainer*)foadb->Get(Form("hQya%i_filtered", i+1));
      if(h[i][1]) fQyavsV0[i] = static_cast<TH1F*>(h[i][1]->GetObject(run));
      h[i][2] = (AliOADBContainer*)foadb->Get(Form("hQxc%i_filtered", i+1));
      if(h[i][2]) fQxcvsV0[i] = static_cast<TH1F*>(h[i][2]->GetObject(run));
      h[i][3] = (AliOADBContainer*)foadb->Get(Form("hQyc%i_filtered", i+1));
      if(h[i][3]) fQycvsV0[i] = static_cast<TH1F*>(h[i][3]->GetObject(run));
    }
    
    // set the recentering style (might be switched back to -1 if recentering is disabeled)
    // FIXME as an ugly hack, for now I mark this as 999 to denote the experimental nature
    // of this and use it transparently without disrupting the existing calbiration
    fApplyRecentering = 999;
}

//-----------------------------------------------------------------------------

void AliFlowEvent::SetDeltaVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts) {
 // implementation of delta vzero calibration 2011
 
 fEvent = cuts->GetEvent();
 if(!fEvent) return; // coverity. we need to know the event to get the runnumber and centrlaity
 // get the vzero centrality percentile (cc dependent calibration)
 
 
 // if this event is from the same run as the previous event
 // we can use the cached calibration values, no need to re-open the
 // aodb file, else cache the new run
 Int_t run(fEvent->GetRunNumber());
 if(fCachedRun == run) return;
 else fCachedRun = run;
 
 // check if the proper chi weights for merging vzero a and vzero c ep are present
 // if not, use sane defaults. centrality binning is equal to that given in the fVZEROcentralityBin snippet
 //
 // chi values can be calculated using the static helper function
 // AliAnalysisTaskJetV2::CalculateEventPlaneChi(Double_t res) where res is the event plane
 // resolution in a given centrality bin
 //
 // the resolutions that were used for these defaults are
 // Double_t R2VZEROA[] = {.35, .40, .48, .50, .48, .45, .38, .26, .16};
 // Double_t R2VZEROC[] = {.45, .60, .70, .73, .68, .60, .40, .36, .17};
 //
 // Double_t R3VZEROA[] = {.22, .23, .22, .19, .15, .12, .08, .00, .00};
 // Double_t R3VZEROC[] = {.30, .30, .28, .25, .22, .17, .11, .00, .00};
 Double_t chiC2[] = {0.771423, 1.10236, 1.38116, 1.48077, 1.31964, 1.10236, 0.674622, 0.600403, 0.273865};
 Double_t chiA2[] = {0.582214, 0.674622, 0.832214, 0.873962, 0.832214, 0.771423, 0.637146, 0.424255, 0.257385};
 Double_t chiC3[] = {0.493347, 0.493347, 0.458557, 0.407166, 0.356628, 0.273865, 0.176208, 6.10352e-05, 6.10352e-05};
 Double_t chiA3[] = {0.356628, 0.373474, 0.356628, 0.306702, 0.24115, 0.192322, 0.127869, 6.10352e-05, 6.10352e-05};
 
 // this may seem redundant but in this way the cuts object is owner of the arrays
 // even if they're created here (so we won't get into trouble with dtor, assigmnet and copying)
 if(!cuts->GetChi2A()) cuts->SetChi2A(new TArrayD(9, chiA2));
 if(!cuts->GetChi2C()) cuts->SetChi2C(new TArrayD(9, chiC2));
 if(!cuts->GetChi3A()) cuts->SetChi3A(new TArrayD(9, chiA3));
 if(!cuts->GetChi3C()) cuts->SetChi3C(new TArrayD(9, chiC3));
 
 if(!fChi2A) fChi2A = cuts->GetChi2A();
 if(!fChi2C) fChi2C = cuts->GetChi2C();
 if(!fChi3A) fChi3A = cuts->GetChi3A();
 if(!fChi3C) fChi3C = cuts->GetChi3C();
 
 // get the calibration file for the gain equalization
 TFile *foadb = TFile::Open("alien:///alice/cern.ch/user/j/jmargutt/gainVZERO.LHC11h.root");
 if(!foadb){
  printf("file alien:///alice/cern.ch/user/j/jmargutt/gainVZERO.LHC11h.root cannot be opened, CALIBRATION FAILED !");
  return;
 }
 TH3F* Weights = dynamic_cast<TH3F*>(foadb->FindObjectAny("LHC11h")->FindObject("gHistVZEROChannelGainEqualizationMap"));
 if(!Weights){
  printf("gHistVZEROChannelGainEqualizationMap is not available in the file\n");
  return;
 }
 
 if(cuts->GetVZEROgainEqualizationPerRing()) {
  // enable or disable rings through the weights, weight 1. is enabled, 0. is disabled
  // start with the vzero c rings (segments 0 through 31)
  (cuts->GetUseVZERORing(0)) ? cuts->SetVZEROCpol(0, 1.) : cuts->SetVZEROCpol(0, 0.);
  (cuts->GetUseVZERORing(1)) ? cuts->SetVZEROCpol(1, 1.) : cuts->SetVZEROCpol(1, 0.);
  (cuts->GetUseVZERORing(2)) ? cuts->SetVZEROCpol(2, 1.) : cuts->SetVZEROCpol(2, 0.);
  (cuts->GetUseVZERORing(3)) ? cuts->SetVZEROCpol(3, 1.) : cuts->SetVZEROCpol(3, 0.);
  // same for vzero a
  (cuts->GetUseVZERORing(4)) ? cuts->SetVZEROApol(0, 1.) : cuts->SetVZEROApol(0, 0.);
  (cuts->GetUseVZERORing(5)) ? cuts->SetVZEROApol(1, 1.) : cuts->SetVZEROApol(1, 0.);
  (cuts->GetUseVZERORing(6)) ? cuts->SetVZEROApol(2, 1.) : cuts->SetVZEROApol(2, 0.);
  (cuts->GetUseVZERORing(7)) ? cuts->SetVZEROApol(3, 1.) : cuts->SetVZEROApol(3, 0.);
 } else {
  // else enable all rings, which is also default
  for(Int_t i(0); i < 4; i++) cuts->SetVZEROCpol(i, 1.);
  for(Int_t i(0); i < 4; i++) cuts->SetVZEROApol(i, 1.);
 }
 
  // loop over the 11h runs to see if it's 11h data (Full TPC Flow)
  Int_t runs11h[] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};
 Int_t RunBin = -1;
  for(Int_t r(0); r < 68; r++) {
   if(run == runs11h[r]) {
    RunBin = r+1;
   }
  }
  if(RunBin>0) {
   printf(" > run has been identified as 11h Full TPC Flow < \n");
  } else {
   printf(" > run has NOT been identified as 11h Full TPC Flow, use default for 11h < \n");
   // pass a NULL pointer to the track cuts object, the NULL pointer will identify 11h runs
   cuts->SetVZEROgainEqualisation(NULL);
   fApplyRecentering = 2011;
   return; // the rest of the steps are not necessary
  }
 
 // step 0) get the TH2D which contains the correction factor per VZERO channel per centrality bin
 //         and pass it to the current track cuts obect
 Weights->GetXaxis()->SetRange(RunBin,RunBin);
 TH2D* fMultVZEROCen = dynamic_cast<TH2D*>(Weights->Project3D("zy"));
 cuts->SetVZEROgainEqualisationCen(fMultVZEROCen);       // passed as a TH2
 
 // get the calibration file for the q-vector recentering
 TFile *fqvec = TFile::Open("alien:///alice/cern.ch/user/j/jmargutt/recenteringVZERO.LHC11h.root");
 if(!fqvec){
  printf("file alien:///alice/cern.ch/user/j/jmargutt/recenteringVZERO.LHC11h.root cannot be opened, CALIBRATION FAILED !");
  return;
 }
 
 // first index of the oadb array is the harmonic n, the second index is either qax, qay, qcx, qcy
 TH2D* h[3][4];
 for(Int_t i(0); i < 3; i++) {
  h[i][0] = (TH2D*)(fqvec->FindObjectAny("LHC11h")->FindObject(Form("gHistVZEROAQ%ixRecenteringMap",i+1)));
  if(h[i][0]) fQxavsV0[i] = static_cast<TH1D*>(h[i][0]->ProjectionY(Form("fQxavsV0[%i]",i),RunBin,RunBin));
  h[i][1] = (TH2D*)(fqvec->FindObjectAny("LHC11h")->FindObject(Form("gHistVZEROAQ%iyRecenteringMap",i+1)));
  if(h[i][1]) fQyavsV0[i] = static_cast<TH1D*>(h[i][1]->ProjectionY(Form("fQyavsV0[%i]",i),RunBin,RunBin));
  h[i][2] = (TH2D*)(fqvec->FindObjectAny("LHC11h")->FindObject(Form("gHistVZEROCQ%ixRecenteringMap",i+1)));
  if(h[i][2]) fQxcvsV0[i] = static_cast<TH1D*>(h[i][2]->ProjectionY(Form("fQxcvsV0[%i]",i),RunBin,RunBin));
  h[i][3] = (TH2D*)(fqvec->FindObjectAny("LHC11h")->FindObject(Form("gHistVZEROCQ%iyRecenteringMap",i+1)));
  if(h[i][3]) fQycvsV0[i] = static_cast<TH1D*>(h[i][3]->ProjectionY(Form("fQycvsV0[%i]",i),RunBin,RunBin));
 }
 
 // set the recentering style (might be switched back to -1 if recentering is disabeled)
 // FIXME as an ugly hack, for now I mark this as 666 to denote the experimental nature
 // of this and use it transparently without disrupting the existing calbiration
 fApplyRecentering = 666;
}

//-----------------------------------------------------------------------------

void AliFlowEvent::ClearFast()
{
  //clear the event without releasing any memory
  //note that cached run number of recentering settigns are not clear
  //(see AliFlowEvent::ClearCachedRun() )
  AliFlowEventSimple::ClearFast();
}

//_____________________________________________________________________________

void AliFlowEvent::ClearCachedRun()
{
    //clear the cached run (not in clear fast as cache needs to be persistent in most cases )
  fCachedRun=0;
  fApplyRecentering=0;
}
