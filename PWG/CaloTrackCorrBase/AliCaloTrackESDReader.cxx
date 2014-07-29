
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

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
// or other particle identification and correlations
//
//
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


//---- ANALYSIS system ----
#include "AliCaloTrackESDReader.h" 
#include "AliAODEvent.h"
#include "AliMultiEventInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

ClassImp(AliCaloTrackESDReader)

//______________________________________________
AliCaloTrackESDReader::AliCaloTrackESDReader() : 
AliCaloTrackReader(), fConstrainTrack(0),
fESDtrackCuts(0), fESDtrackComplementaryCuts(0)
{
  //Default Ctor
  
  //Initialize parameters
  fDataType           = kESD;
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
  fConstrainTrack     = kFALSE ; // constrain tracks to vertex

}

//_____________________________________________
AliCaloTrackESDReader::~AliCaloTrackESDReader()
{
  //Dtor
  
  AliCaloTrackReader::~AliCaloTrackReader();
  
  delete fESDtrackCuts;
  delete fESDtrackComplementaryCuts;
}

//_________________________________________________________
Bool_t AliCaloTrackESDReader::CheckForPrimaryVertex() const
{
  //Check if the vertex was well reconstructed, copy of conversion group
  
  AliESDEvent * esdevent = dynamic_cast<AliESDEvent*> (fInputEvent);
  if(!esdevent) return kFALSE;
  
  if(esdevent->GetPrimaryVertex()->GetNContributors() > 0)
  {
    return kTRUE;
  }
  
  if(esdevent->GetPrimaryVertex()->GetNContributors() < 1)
  {
    // SPD vertex
    if(esdevent->GetPrimaryVertexSPD()->GetNContributors() > 0)
    {
      return kTRUE;
      
    }
    if(esdevent->GetPrimaryVertexSPD()->GetNContributors() < 1)
    {
      return kFALSE;
    }
  }

  return kFALSE;

}

//________________________________
void AliCaloTrackESDReader::Init()
{
  //Init reader. Method to be called in AliAnaCaloTrackCorrMaker
  
  AliCaloTrackReader::Init();
  
  if(!fESDtrackCuts)
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //initialize with TPC only tracks
}

//______________________________________________________________________________
Bool_t AliCaloTrackESDReader::SelectTrack(AliVTrack* track, Double_t pTrack[3])
{
  // Select ESD track using the cuts declared in fESDtrackCuts
  // in case of hybrid tracks, 2 different sets of cuts defined.
  
  AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
  
  if(!esdTrack) return kFALSE;
  
  const AliExternalTrackParam* constrainParam = esdTrack->GetConstrainedParam();
  
  if(fESDtrackCuts->AcceptTrack(esdTrack))
  {
    track->GetPxPyPz(pTrack) ;
    
    if(fConstrainTrack)
    {
      if( !constrainParam ) return kFALSE;
      
      esdTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
      esdTrack->GetConstrainedPxPyPz(pTrack);
      
    } // use constrained tracks
    
    if(fSelectSPDHitTracks && !esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1))
      return kFALSE ; // Not much sense to use with TPC only or Hybrid tracks
  }
  
  // Complementary track to global : Hybrids (make sure that the previous selection is for Global)
  else if(fESDtrackComplementaryCuts && fESDtrackComplementaryCuts->AcceptTrack(esdTrack))
  {
    // constrain the track
    if( !constrainParam ) return kFALSE;
    
    esdTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
    esdTrack->GetConstrainedPxPyPz(pTrack);
    
  }
  else return kFALSE;
  
  return kTRUE;
}

//_______________________________________________________________
void  AliCaloTrackESDReader::SetTrackCuts(AliESDtrackCuts * cuts)
{
  // Set Track cuts
  
  if(fESDtrackCuts) delete fESDtrackCuts ;
  
  fESDtrackCuts = cuts ;
  
}

//____________________________________________________________________________
void  AliCaloTrackESDReader::SetTrackComplementaryCuts(AliESDtrackCuts * cuts)
{
  // Set Track cuts for complementary tracks (hybrids)
  
  if(fESDtrackComplementaryCuts) delete fESDtrackComplementaryCuts ;
  
  fESDtrackComplementaryCuts = cuts ;
  
}

//_________________________________________________________________
void AliCaloTrackESDReader::SetInputOutputMCEvent(AliVEvent* esd,
                                                  AliAODEvent* aod,
                                                  AliMCEvent* mc) 
{
  // Connect the data pointers
  
  Bool_t tesd = kFALSE ; 
  
  if ( strcmp(esd->GetName(), "AliMixedEvent") == 0 ) 
  {
    AliMultiEventInputHandler* multiEH = dynamic_cast<AliMultiEventInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(multiEH)
    {
      if (multiEH->GetFormat() == 0 ) 
      {
        tesd = kTRUE ; 
      }
    }
    else
    {
      printf("AliCaloTrackESDReader::SetInputOutputMCEvent() - MultiEventHandler is NULL");
      abort();
    }
  }
  if (strcmp(esd->GetName(),"AliESDEvent") == 0)
  {
    tesd = kTRUE ; 
  }
  
  if(!tesd)
  {
    AliFatal(Form("AliCaloTrackESDReader::SetInputOutputMCEvent() - STOP ::Wrong reader, here only ESDs. Input name: %s != AliESDEvent \n",esd->GetName()));
  }
  
  SetInputEvent(esd);
  SetOutputEvent(aod);
  SetMC(mc);
  
}



