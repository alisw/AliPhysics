
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

//---- ANALYSIS system ----
#include "AliCaloTrackESDReader.h" 
#include "AliAODEvent.h"
#include "AliMultiEventInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliCaloTrackESDReader) ;
/// \endcond 

//______________________________________________
/// Default constructor. Initialize parameters.
//______________________________________________
AliCaloTrackESDReader::AliCaloTrackESDReader() : 
AliCaloTrackReader(), fConstrainTrack(0),
fESDtrackCuts(0), fESDtrackComplementaryCuts(0)
{
  fDataType           = kESD;
  fConstrainTrack     = kFALSE ; // constrain tracks to vertex
}

//_____________________________________________
/// Default destructor.
//_____________________________________________
AliCaloTrackESDReader::~AliCaloTrackESDReader()
{  
  AliCaloTrackReader::DeletePointers();
  
  delete fESDtrackCuts;
  delete fESDtrackComplementaryCuts;
}

//_________________________________________________________
/// Check if the vertex was well reconstructed.
/// Copy method of PCM group.
//_________________________________________________________
Bool_t AliCaloTrackESDReader::CheckForPrimaryVertex() const
{  
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

//______________________________________________________________
/// \return pointer to Generated event header (AliGenEventHeader)
//______________________________________________________________
AliGenEventHeader* AliCaloTrackESDReader::GetGenEventHeader() const
{
  if ( !fMC ) return 0x0 ;
  
  AliGenEventHeader * eventHeader = fMC->GenEventHeader();
  
  if ( fMCGenerEventHeaderToAccept=="" ) return eventHeader ;
  
  AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
  
  if ( !cocktail ) return 0x0 ;
  
  TList *genHeaders = cocktail->GetHeaders();
  
  Int_t nGenerators = genHeaders->GetEntries();
  
  for(Int_t igen = 0; igen < nGenerators; igen++)
  {
    AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
    
    TString name = eventHeader2->GetName();
    //printf("ESD Event header %d %s\n",igen,name.Data());
    
    if(name.Contains(fMCGenerEventHeaderToAccept,TString::kIgnoreCase)) 
      return eventHeader2 ;
  }
  
  return 0x0;
}

//________________________________
/// Init reader. Method to be called in AliAnaCaloTrackCorrMaker.
//________________________________
void AliCaloTrackESDReader::Init()
{  
  AliCaloTrackReader::Init();
  
  if(!fESDtrackCuts)
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //initialize with TPC only tracks
}

//______________________________________________________________________________
/// Select ESD track using the cuts declared in *fESDtrackCuts*.
/// In case of hybrid tracks, 2 different sets of cuts defined.
//______________________________________________________________________________
Bool_t AliCaloTrackESDReader::SelectTrack(AliVTrack* track, Double_t pTrack[3])
{  
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
/// Set Track cuts.
//_______________________________________________________________
void  AliCaloTrackESDReader::SetTrackCuts(AliESDtrackCuts * cuts)
{  
  if(fESDtrackCuts) delete fESDtrackCuts ;
  
  fESDtrackCuts = cuts ;
}

//____________________________________________________________________________
/// Set Track cuts for complementary tracks (hybrids).
//____________________________________________________________________________
void  AliCaloTrackESDReader::SetTrackComplementaryCuts(AliESDtrackCuts * cuts)
{  
  if(fESDtrackComplementaryCuts) delete fESDtrackComplementaryCuts ;
  
  fESDtrackComplementaryCuts = cuts ;
}

//_________________________________________________________________
/// Connect the data pointers.
//_________________________________________________________________
void AliCaloTrackESDReader::SetInputOutputMCEvent(AliVEvent* esd,
                                                  AliAODEvent* aod,
                                                  AliMCEvent* mc) 
{  
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
      AliFatal("MultiEventHandler is NULL");
      return;
    }
  }
  if (strcmp(esd->GetName(),"AliESDEvent") == 0)
  {
    tesd = kTRUE ; 
  }
  
  if(!tesd)
  {
    AliFatal(Form("STOP ::Wrong reader, here only ESDs. Input name: %s != AliESDEvent",esd->GetName()));
  }
  
  SetInputEvent(esd);
  SetOutputEvent(aod);
  SetMC(mc);
}



