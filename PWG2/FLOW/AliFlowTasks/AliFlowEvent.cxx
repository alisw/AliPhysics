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
#include "TH2F.h"
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
#include "AliMultiplicity.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrack.h"
#include "AliFlowEvent.h"
#include "AliLog.h"

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
  //EPOS
  else if (!strcmp(mcEvent->GenEventHeader()->GetName(),"EPOS"))
  {
    AliGenEposEventHeader* headerE = dynamic_cast<AliGenEposEventHeader*>(mcEvent->GenEventHeader());
    if (headerE) AliFlowEventSimple::SetMCReactionPlaneAngle( headerE->ReactionPlaneAngle() );
  }
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliMCEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager):
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
      fNumberOfRPs++;
    }
    if (poiOK && poiCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
    }

    AddTrack(pTrack) ;
  }//for all tracks
  SetMCReactionPlaneAngle(anInput);
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager ):
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
      fNumberOfRPs++;
    }
    //marking the particles used for diff. flow:
    if(poiOK && poiCFManager)
    {
      pTrack->SetForPOISelection(kTRUE);
    }

    AddTrack(pTrack);
  }//end of while (itrkN < iNumberOfInputTracks)
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliAODEvent* anInput,
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager):
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
      fNumberOfRPs++;
    }
    if (poiOK /* && poiCFManager*/ )
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
                            const AliCFManager* rpCFManager,
                            const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20)
{
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
    Bool_t rpOK = kTRUE;
    Bool_t poiOK = kTRUE;
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
      fNumberOfRPs++;
      pTrack->SetForRPSelection();
    }
    if (poiOK && poiCFManager) pTrack->SetForPOISelection();

    AddTrack(pTrack);
  }
  SetMCReactionPlaneAngle(anInputMc);
}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
			    const AliMultiplicity* anInputTracklets,
			    const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20)
{

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
    fNumberOfRPs++;
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
  AliFlowEventSimple(20)
{

  //Select the particles of interest from the ESD
  Int_t iNumberOfInputTracks = esd->GetNumberOfTracks() ;

  //Double_t gPt = 0.0, gP = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;

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

        if(hybrid)
          tpcTrack->PropagateToDCA(vertexSPD,esd->GetMagneticField(),100.,dca,cov);
        else
          tpcTrack->PropagateToDCA(vertexTPC,esd->GetMagneticField(),100.,dca,cov);

        dca3D = TMath::Sqrt(TMath::Power(dca[0],2)+TMath::Power(dca[1],2));

      }

      //make new AliFLowTrack
      AliFlowTrack* pTrack = new AliFlowTrack(pParticle);

      pTrack->SetSource(AliFlowTrack::kFromESD);

      //marking the particles used for diff. flow:
      if(poiOK && poiCFManager)
      {
        pTrack->SetForPOISelection(kTRUE);
      }

      if(useTPC)
      {
        pTrack->SetForRPSelection(kTRUE);
        fNumberOfRPs++;
      }

      AddTrack(pTrack);

    }//end of while (itrkN < iNumberOfInputTracks)

}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( const AliESDEvent* anInput,
			    const TH2F* anInputFMDhist,
			    const AliCFManager* poiCFManager ):
  AliFlowEventSimple(20)
{

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
	pTrack->SetForRPSelection(kTRUE);
	pTrack->SetSource(AliFlowTrack::kFromFMD);
	fNumberOfRPs++;

	//Add the track to the flowevent
	AddTrack(pTrack);
	
      }
    }
  }

}

//-----------------------------------------------------------------------
AliFlowEvent::AliFlowEvent( AliVEvent* event,
                            AliFlowTrackCuts* rpCuts,
                            AliFlowTrackCuts* poiCuts ):
  AliFlowEventSimple(20)
{
  //Fills the event from a vevent: AliESDEvent,AliAODEvent,AliMCEvent

  if (!rpCuts || !poiCuts) return;

  //if input event empty try to do MC analysis
  if (!event) event = rpCuts->GetMCevent();
  if (!event) return;

  Int_t numberOfTracks = event->GetNumberOfTracks() ;

  //loop over tracks
  for (Int_t i=0; i<numberOfTracks; i++)
  {
    AliVParticle* particle = event->GetTrack(i);   //get input particle

    Bool_t rp = rpCuts->IsSelected(particle);
    Bool_t poi = poiCuts->IsSelected(particle);
    
    if (!(rp||poi)) continue;

    //make new AliFLowTrack
    //here we need to be careful: if selected particle passes both rp and poi cuts
    //then both cuts should have been done on the same set of parameters, e.g. global
    //or TPConly. Otherwise we would have to introduce the same particle twice.
    //this means that in a sane scenario when we pass both rp and poi cuts we get our
    //parameters from any one of them (here rp).
    AliFlowTrack* pTrack = NULL;
    if (rp&&poi)
    {
      pTrack = rpCuts->MakeFlowTrack();
      pTrack->TagRP(); fNumberOfRPs++;
      pTrack->TagPOI();
    }
    else
    if (rp)
    {
      pTrack = rpCuts->MakeFlowTrack();
      pTrack->TagRP(); fNumberOfRPs++;
    }
    else
    if (poi)
    {
      pTrack = poiCuts->MakeFlowTrack();
      pTrack->TagPOI();
    }

    AddTrack(pTrack);
  }//end of while (i < numberOfTracks)
}

