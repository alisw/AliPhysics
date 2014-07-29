
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
// Class for reading data (AODs) in order to do prompt gamma
//  or other particle identification and correlations.
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
// 
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////

//---- ANALYSIS system ----
#include "AliCaloTrackAODReader.h" 
#include "AliAODInputHandler.h"
#include "AliMultiEventInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"
#include "AliAODEvent.h"

ClassImp(AliCaloTrackAODReader)

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader() : 
  AliCaloTrackReader(),   fOrgInputEvent(0x0),
  fSelectHybridTracks(0), fSelectPrimaryTracks(0),
  fTrackFilterMask(0),    fTrackFilterMaskComplementary(0),
  fSelectFractionTPCSharedClusters(0), fCutTPCSharedClustersFraction(0)
{
  //Default Ctor
  
  //Initialize parameters
  fDataType = kAOD;
  
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
 
  fTrackFilterMask = 128;
  fTrackFilterMaskComplementary = 0; // in case of hybrid tracks, without using the standard method
  
  fSelectFractionTPCSharedClusters = kTRUE;
  fCutTPCSharedClustersFraction    = 0.4;
  
}

//_________________________________________________________
Bool_t AliCaloTrackAODReader::CheckForPrimaryVertex() const
{
  //Check if the vertex was well reconstructed, copy of conversion group
  
  AliAODEvent * aodevent = dynamic_cast<AliAODEvent*>(fInputEvent);
  if(!aodevent) return kFALSE;
  
  if (aodevent->GetPrimaryVertex() != NULL)
  {
    if(aodevent->GetPrimaryVertex()->GetNContributors() > 0)
    {
      return kTRUE;
    }
  }
  
  if(aodevent->GetPrimaryVertexSPD() != NULL)
  {
    if(aodevent->GetPrimaryVertexSPD()->GetNContributors() > 0)
    {
      return kTRUE;
    }
    else
    {
      AliWarning(Form("Number of contributors from bad vertex type:: %s",
                      aodevent->GetPrimaryVertex()->GetName()));
      return kFALSE;
    }
  }
  
  return kFALSE;
  
}

//_____________________________________________________________________________
Bool_t AliCaloTrackAODReader::SelectTrack(AliVTrack* track, Double_t pTrack[3])
{
  // Select AOD track using the AOD filter bits
  
  AliAODTrack *aodtrack = dynamic_cast <AliAODTrack*>(track);
  
  if(!aodtrack) return kFALSE;
  
  
  if(fDebug > 2 ) printf("AliCaloTrackAODReader::FillInputCTS():AOD track type: %d (primary %d), hybrid? %d \n",
                         aodtrack->GetType(),AliAODTrack::kPrimary,
                         aodtrack->IsHybridGlobalConstrainedGlobal());
  
  // Hybrid?
  if (fSelectHybridTracks && fTrackFilterMaskComplementary == 0)
  {
    if (!aodtrack->IsHybridGlobalConstrainedGlobal()) return kFALSE ;
  }
  else // Filter Bit?
  {
    Bool_t accept = aodtrack->TestFilterBit(fTrackFilterMask);
    
    if(!fSelectHybridTracks && !accept) return kFALSE ;
    
    if(fSelectHybridTracks) // Second filter bit for hybrids?
    {
      Bool_t acceptcomplement = aodtrack->TestFilterBit(fTrackFilterMaskComplementary);
      if (!accept && !acceptcomplement) return kFALSE ;
    }
  }
  
  //
  if(fSelectSPDHitTracks)
  { // Not much sense to use with TPC only or Hybrid tracks
    if(!aodtrack->HasPointOnITSLayer(0) && !aodtrack->HasPointOnITSLayer(1)) return kFALSE ;
  }
  
  //
  if ( fSelectFractionTPCSharedClusters )
  {
    Double_t frac = Double_t(aodtrack->GetTPCnclsS()) / Double_t(aodtrack->GetTPCncls());
    if (frac > fCutTPCSharedClustersFraction)
    {
      if (fDebug > 2 )printf("\t Reject track, shared cluster fraction %f > %f\n",frac, fCutTPCSharedClustersFraction);
      return kFALSE ;
    }
  }
  
  //
  if ( fSelectPrimaryTracks )
  {
    if ( aodtrack->GetType()!= AliAODTrack::kPrimary )
    {
      if (fDebug > 2 ) printf("\t Remove not primary track\n");
      return kFALSE ;
    }
  }
  
  if (fDebug > 2 ) printf("\t accepted track! \n");
  
  track->GetPxPyPz(pTrack) ;
  
  return kTRUE;
  
}


//_________________________________________________________________
void AliCaloTrackAODReader::SetInputOutputMCEvent(AliVEvent* input,
                                                  AliAODEvent* aod,
                                                  AliMCEvent* mc)
{
  // Connect the data pointers
  // If input is AOD, do analysis with input, if not, do analysis with the output aod.
  
  //printf("AODInputHandler %p, MergeEvents %d \n",aodIH, aodIH->GetMergeEvents());
  
  Bool_t tesd = kFALSE ; 
  Bool_t taod = kTRUE ; 
  if ( strcmp(input->GetName(), "AliMixedEvent") == 0 ) 
  {
    AliMultiEventInputHandler* multiEH = dynamic_cast<AliMultiEventInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if(multiEH){
      if (multiEH->GetFormat() == 0 ) 
      {
        tesd = kTRUE ; 
      } else if (multiEH->GetFormat() == 1) 
      {
        taod = kTRUE ; 
      }
    }
    else
    {
      printf("AliCaloTrackAODReader::SetInputOutputMCEvent() - MultiEventHandler is NULL");
      abort();
    }
  }
  if        (strcmp(input->GetName(),"AliESDEvent") == 0) 
  {
    tesd = kTRUE ; 
  } else if (strcmp(input->GetName(),"AliAODEvent") == 0) 
  {
    taod = kTRUE ; 
  }
  
  
  if(tesd)   
  {
    SetInputEvent(aod);
    SetOutputEvent(aod);
    fOrgInputEvent = input;
  }
  else if(taod)
  {
    AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    
	  if (aodIH && aodIH->GetMergeEvents()) 
    {
		  //Merged events, use output AOD.
		  SetInputEvent(aod);
		  SetOutputEvent(aod);
      fOrgInputEvent = input;
	  }
	  else
    {
		  SetInputEvent(input);
		  SetOutputEvent(aod);
	  }
  }
  else
  { 
    AliFatal(Form("AliCaloTrackAODReader::SetInputOutputMCEvent() - STOP : Wrong data format: %s\n",input->GetName()));
  }
  
  SetMC(mc);
  
}

