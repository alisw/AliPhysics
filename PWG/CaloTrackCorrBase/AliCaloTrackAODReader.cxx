
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
#include "AliCaloTrackAODReader.h" 
#include "AliAODInputHandler.h"
#include "AliMultiEventInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"
#include "AliAODEvent.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliCaloTrackAODReader) ;
/// \endcond

//______________________________________________
/// Default constructor. Initialize parameters
//______________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader() : 
  AliCaloTrackReader(),   fOrgInputEvent(0x0),
  fSelectHybridTracks(0), fSelectPrimaryTracks(0),
  fTrackFilterMask(0),    fTrackFilterMaskComplementary(0),
  fSelectFractionTPCSharedClusters(0), fCutTPCSharedClustersFraction(0)
{
  fDataType = kAOD;
  
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
 
  fTrackFilterMask = 128;
  fTrackFilterMaskComplementary = 0; // in case of hybrid tracks, without using the standard method
  
  fSelectFractionTPCSharedClusters = kTRUE;
  fCutTPCSharedClustersFraction    = 0.4;
}

//_________________________________________________________
/// Check if the vertex was well reconstructed.
/// Copy method of PCM group.
//_________________________________________________________
Bool_t AliCaloTrackAODReader::CheckForPrimaryVertex() const
{  
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

//____________________________________________________________
/// \return list of MC particles in AOD. Do it for the corresponding input event.
//____________________________________________________________
TClonesArray* AliCaloTrackAODReader::GetAODMCParticles() const
{  
  TClonesArray * particles = NULL ;

  AliAODEvent * aod = dynamic_cast<AliAODEvent*> (fInputEvent) ;
  if(aod) particles = (TClonesArray*) aod->FindListObject("mcparticles");
  
  return particles ;
}

//___________________________________________________________
/// \return MC header in AOD. Do it for the corresponding input event.
//___________________________________________________________
AliAODMCHeader* AliCaloTrackAODReader::GetAODMCHeader() const
{  
  AliAODMCHeader *mch = NULL;
  
  AliAODEvent * aod = dynamic_cast<AliAODEvent*> (fInputEvent);
  if(aod) mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
  
  return mch;
}

//_____________________________________________________________________________
/// Select AOD track using the AOD filter bits or predefined selection methods.
//_____________________________________________________________________________
Bool_t AliCaloTrackAODReader::SelectTrack(AliVTrack* track, Double_t pTrack[3])
{  
  AliAODTrack *aodtrack = dynamic_cast <AliAODTrack*>(track);
  
  if(!aodtrack) return kFALSE;
  
  AliDebug(2,Form("AOD track type: %d (primary %d), hybrid? %d",
                  aodtrack->GetType(),AliAODTrack::kPrimary,
                  aodtrack->IsHybridGlobalConstrainedGlobal()));
  
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
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if (frac > fCutTPCSharedClustersFraction)
    {
      AliDebug(2,Form("\t Reject track, shared cluster fraction %f > %f",frac, fCutTPCSharedClustersFraction));
      return kFALSE ;
    }
  }
  
  //
  if ( fSelectPrimaryTracks )
  {
    if ( aodtrack->GetType()!= AliAODTrack::kPrimary )
    {
      AliDebug(2,"\t Remove not primary track");
      return kFALSE ;
    }
  }
  
  AliDebug(2,"\t accepted track!");
  
  track->GetPxPyPz(pTrack) ;
  
  return kTRUE;
}

//_________________________________________________________________
/// Connect the data pointers
/// If input is AOD, do analysis with input, if not, do analysis with the output aod.
//_________________________________________________________________
void AliCaloTrackAODReader::SetInputOutputMCEvent(AliVEvent* input,
                                                  AliAODEvent* aod,
                                                  AliMCEvent* mc)
{
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
      AliFatal("MultiEventHandler is NULL");
      return;
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
    AliFatal(Form("STOP : Wrong data format: %s",input->GetName()));
  }
  
  SetMC(mc);
}

