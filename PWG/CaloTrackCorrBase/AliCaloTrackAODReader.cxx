
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
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include <TObjString.h>

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
 
  //fTrackFilterMask = 128;
  //fTrackFilterMaskComplementary = 0; // in case of hybrid tracks, without using the standard method
  
  //fSelectFractionTPCSharedClusters = kTRUE;
  fCutTPCSharedClustersFraction    = 0.4;
  
  for(Int_t i = 0; i < 8; i++)
  {
    fhCTSAODTrackCutsPt   [i] = 0;
    fhCTSAODTrackCutsPtCen[i] = 0;
    fhCTSAODTrackCutsPtSignal   [i] = 0;
    fhCTSAODTrackCutsPtCenSignal[i] = 0;
  }
}

//_________________________________________________________
/// Check if the vertex was well reconstructed.
/// Copy method of PCM group.
//_________________________________________________________
Bool_t AliCaloTrackAODReader::CheckForPrimaryVertex() const
{  
  AliAODEvent * aodevent = NULL;
  // In case of analysis of pure MC event used in embedding 
  // get bkg PbPb event since cuts for embedded event are based on 
  // data vertex and not on pp simu vertex
  if ( fEmbeddedEvent[1] ) // Input event is MC
    aodevent = dynamic_cast<AliAODEvent*>(((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetEvent());
  else 
    aodevent = dynamic_cast<AliAODEvent*>(fInputEvent);
  
  if ( !aodevent ) return kFALSE;
  
  if ( aodevent->GetPrimaryVertex() != NULL )
  {
    if ( aodevent->GetPrimaryVertex()->GetNContributors() > 0 )
    {
      return kTRUE;
    }
  }
  
  if ( aodevent->GetPrimaryVertexSPD() != NULL )
  {
    if ( aodevent->GetPrimaryVertexSPD()->GetNContributors() > 0 )
    {
      return kTRUE;
    }
    else
    {
      AliDebug(1,Form("Null number of contributors from bad vertex type:: %s",
                      aodevent->GetPrimaryVertex()->GetName()));
      return kFALSE;
    }
  }
  
  return kFALSE;
}


//___________________________________________________
/// Fill the output list of initialized control histograms.
/// Cluster or track spectra histograms, depending on different selection cuts.
/// First fill the histograms of the mother class, that are independent on ESD/AOD.
/// Then add the AOD specific ones.
//___________________________________________________
TList * AliCaloTrackAODReader::GetCreateControlHistograms()
{  
  AliCaloTrackReader::GetCreateControlHistograms();
    
  if ( fFillCTS )
  {
    TString names[] = {"FilterBit_Hybrid", "SPDHit", "NclustersITS", "Chi2ITS",
      "NclustersTPC", "Chi2TPC", "SharedCluster", "Primary"};
    
    for(Int_t i = 0; i < 8; i++)
    {
      if ( names[i].Contains("Filter") && !fTrackFilterMask && 
          !fTrackFilterMaskComplementary && !fSelectHybridTracks            ) continue;
      if ( names[i] == "SPDHit"        && !fSelectSPDHitTracks              ) continue;
      if ( names[i] == "NclustersITS"  && fSelectMinITSclusters <= 0        ) continue;
      if ( names[i] == "Chi2ITS"       && fSelectMaxChi2PerITScluster > 1000) continue;
      if ( names[i] == "NclustersTPC"  && fSelectMinTPCclusters <= 0        ) continue;
      if ( names[i] == "Chi2TPC"       && fSelectMaxChi2PerTPCcluster > 1000) continue;
      if ( names[i] == "SharedCluster" && !fSelectFractionTPCSharedClusters ) continue;
      if ( names[i] == "Primary"       && !fSelectPrimaryTracks             ) continue;
      
      if ( !IsHistoCentDependentOn() )
      {
        fhCTSAODTrackCutsPt[i] = new TH1F
        (Form("hCTSReaderAODTrackCuts_%d_%s",i,names[i].Data()),
         Form("AOD CTS Cut %d, %s",i,names[i].Data()), 
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
        fhCTSAODTrackCutsPt[i]->SetYTitle("# tracks");
        fhCTSAODTrackCutsPt[i]->SetXTitle("#it{p}_{T} (GeV)");
        fOutputContainer->Add(fhCTSAODTrackCutsPt[i]);

        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
        {
          fhCTSAODTrackCutsPtSignal[i] = new TH1F
          (Form("hCTSReaderAODTrackCutsSignal_%d_%s",i,names[i].Data()),
           Form("AOD CTS Cut %d, %s",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1]) ;
          fhCTSAODTrackCutsPtSignal[i]->SetYTitle("# tracks");
          fhCTSAODTrackCutsPtSignal[i]->SetXTitle("#it{p}_{T} (GeV)");
          fOutputContainer->Add(fhCTSAODTrackCutsPtSignal[i]);
        }
      }
      else
      {
        fhCTSAODTrackCutsPtCen[i] = new TH2F
        (Form("hCTSReaderAODTrackCutsCen_%d_%s",i,names[i].Data()),
         Form("AOD CTS Cut %d, %s",i,names[i].Data()), 
         fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
         100, 0, 100) ;
        fhCTSAODTrackCutsPtCen[i]->SetYTitle("# tracks");
        fhCTSAODTrackCutsPtCen[i]->SetXTitle("#it{p}_{T} (GeV)");
        fhCTSAODTrackCutsPtCen[i]->SetYTitle("Centrality (%)");
        fOutputContainer->Add(fhCTSAODTrackCutsPtCen[i]);

        if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
        {
          fhCTSAODTrackCutsPtCenSignal[i] = new TH2F
          (Form("hCTSReaderAODTrackCutsCenSignal_%d_%s",i,names[i].Data()),
           Form("AOD CTS Cut %d, %s",i,names[i].Data()),
           fEnergyHistogramNbins, fEnergyHistogramLimit[0], fEnergyHistogramLimit[1],
           100, 0, 100) ;
          fhCTSAODTrackCutsPtCenSignal[i]->SetYTitle("# tracks");
          fhCTSAODTrackCutsPtCenSignal[i]->SetXTitle("#it{p}_{T} (GeV)");
          fhCTSAODTrackCutsPtCenSignal[i]->SetYTitle("Centrality (%)");
          fOutputContainer->Add(fhCTSAODTrackCutsPtCenSignal[i]);
        }
      }
    }
  }
  
  return fOutputContainer ;
}

//________________________________________________________
/// Save parameters used for analysis in a string.
//________________________________________________________
TObjString *  AliCaloTrackAODReader::GetListOfParameters()
{
  // Recover the string from the mother class
  TString parList = (AliCaloTrackReader::GetListOfParameters())->GetString();
  
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"AOD Track: Hybrid %d, Filter bit %d, Complementary bit %d, Primary %d; ", 
           fSelectHybridTracks, (Int_t)fTrackFilterMask, (Int_t)fTrackFilterMaskComplementary, fSelectPrimaryTracks) ;
  parList+=onePar ;
  
  if ( fSelectFractionTPCSharedClusters )
  {
    snprintf(onePar,buffersize,"Fraction of TPC shared clusters ON: %2.2f ", fCutTPCSharedClustersFraction) ;
    parList+=onePar ;
  }
  
  return new TObjString(parList) ;

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
  
  AliAODEvent * aod =  NULL; 
  
  if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1] )
    aod =  dynamic_cast<AliAODEvent*> (AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent());
  else 
    aod = dynamic_cast<AliAODEvent*> (fInputEvent);
  
  if(aod) mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
  
  return mch;
}

//______________________________________________________________
/// \return pointer to Generated event header (AliGenEventHeader)
//______________________________________________________________
AliGenEventHeader* AliCaloTrackAODReader::GetGenEventHeader() const
{
  if ( !GetAODMCHeader() ) return 0x0;
  
  if ( fGenEventHeader ) return fGenEventHeader;
  
  Int_t nGenerators = GetAODMCHeader()->GetNCocktailHeaders();
  
  if ( nGenerators <= 0  ) return 0x0;
  
  if ( fMCGenerEventHeaderToAccept=="" ) 
    return GetAODMCHeader()->GetCocktailHeader(0);
  
  //if(nGenerators < 1) printf("No generators %d\n",nGenerators);
  
  for(Int_t igen = 0; igen < nGenerators; igen++)
  {
    AliGenEventHeader * eventHeader = GetAODMCHeader()->GetCocktailHeader(igen) ;
   
    TString name = eventHeader->GetName();
    
    AliDebug(2,Form("AOD Event header %d name %s event header class %s: Select if contains <%s>",
             igen,name.Data(),eventHeader->ClassName(),fMCGenerEventHeaderToAccept.Data()));
    
//   if(nGenerators < 2) 
//     printf("AOD Event header %d name %s event header class %s: Select if contains <%s>\n",
//           igen,name.Data(),eventHeader->ClassName(),fMCGenerEventHeaderToAccept.Data());
    
    if(name.Contains(fMCGenerEventHeaderToAccept,TString::kIgnoreCase))
      return eventHeader ;
  }
  
  return 0x0;
}

//_____________________________
/// 
/// Return track ID, different for ESD and AODs.
/// In AODs a track can correspond to several definitions.
/// If negative (hybrid constrained to vertex), it corresponds to global track 
/// with ID positive plus one.
///
/// See AliCaloTrackReader for ESD correspondent.
///
/// \return track ID
/// \param track: pointer to track
//_____________________________
Int_t AliCaloTrackAODReader::GetTrackID(AliVTrack* track) 
{
  Int_t id = track->GetID();
  
  if( id < 0 ) id = TMath::Abs(id) - 1;
  
  return id;
}


//________________________________________________________
/// Print parameters
//________________________________________________________
void AliCaloTrackAODReader::Print(const Option_t * opt) const
{  
  if(! opt)
    return;
  
  AliCaloTrackReader::Print(opt);
  
  printf("AOD Track: Hybrid %d, Filter bit %d, Complementary bit %d, Primary %d; \n", 
           fSelectHybridTracks, (Int_t)fTrackFilterMask, (Int_t)fTrackFilterMaskComplementary, fSelectPrimaryTracks) ;
  
  if ( fSelectFractionTPCSharedClusters )
  {
    printf("Fraction of TPC shared clusters ON: %2.2f ", fCutTPCSharedClustersFraction) ;
  }
}

//_____________________________________________________________________________
/// Select AOD track using the AOD filter bits or predefined selection methods.
//_____________________________________________________________________________
Bool_t AliCaloTrackAODReader::SelectTrack(AliVTrack* track, Double_t pTrack[3])
{  
  AliAODTrack *aodtrack = dynamic_cast <AliAODTrack*>(track);
  
  if(!aodtrack) return kFALSE;
  
  track->GetPxPyPz(pTrack) ;
  
  AliDebug(2,Form("AOD track type: %d (primary %d), hybrid? %d",
                  aodtrack->GetType(),AliAODTrack::kPrimary,
                  aodtrack->IsHybridGlobalConstrainedGlobal()));

  Int_t cen = GetEventCentrality();

  // Hybrid?
  if ( fSelectHybridTracks && fTrackFilterMaskComplementary == 0 )
  {
    if ( !aodtrack->IsHybridGlobalConstrainedGlobal() ) return kFALSE ;
  }
  else if ( fTrackFilterMask > 0 || fTrackFilterMaskComplementary > 0 ) // Filter Bit?
  {
    Bool_t accept = aodtrack->TestFilterBit(fTrackFilterMask);
    
    if ( !fSelectHybridTracks && !accept ) return kFALSE ;
    
    if ( fSelectHybridTracks ) // Second filter bit for hybrids?
    {
      Bool_t acceptcomplement = aodtrack->TestFilterBit(fTrackFilterMaskComplementary);
      if ( !accept && !acceptcomplement ) return kFALSE ;
    }
  }

  Bool_t fillEmbedSignalTrack = kFALSE;
  if ( fEmbeddedEvent[0] && !fEmbeddedEvent[1]  && aodtrack->GetLabel() >=0 )
    fillEmbedSignalTrack = kTRUE;

  if ( fSelectHybridTracks || fTrackFilterMaskComplementary || fTrackFilterMask )
  {
    AliDebug(2,"Pass mask cut");

    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [0]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[0]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [0]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[0]->Fill(track->Pt(),cen);
    }
  }
  
  //
  // ITS related cuts
  // Not much sense to use with TPC only or Hybrid tracks
  //
  
//  printf("SPD %d (%d,%d) - n clus %d >= %d - chi2 %f max chi2/n cls %f\n",
//         fSelectSPDHitTracks, aodtrack->HasPointOnITSLayer(0),!aodtrack->HasPointOnITSLayer(1),
//         aodtrack->GetITSNcls(),fSelectMinITSclusters,
//         aodtrack->GetITSchi2(),fSelectMaxChi2PerITScluster);
  
  if ( fSelectSPDHitTracks )
  { 
    if ( !aodtrack->HasPointOnITSLayer(0) && !aodtrack->HasPointOnITSLayer(1) ) 
      return kFALSE ;
    
    AliDebug(2,"Pass SPD layer cut");
    
    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [1]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[1]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [1]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[1]->Fill(track->Pt(),cen);
    }
  }
  
  Int_t nITScls = aodtrack->GetITSNcls();
  if ( fSelectMinITSclusters > 0 )
  {
    if ( nITScls < fSelectMinITSclusters  ) 
    {
      return kFALSE;    
    }
    
    AliDebug(2,"Pass n ITS cluster cut");
    
    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [2]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[2]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [2]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[2]->Fill(track->Pt(),cen);
    }
  }
  
  Float_t chi2PerITScluster = 0.;
  if ( nITScls > 0 ) chi2PerITScluster = aodtrack->GetITSchi2()/nITScls;
  
  if ( chi2PerITScluster > fSelectMaxChi2PerITScluster ) 
  {
    return kFALSE; 
  }
  
  if ( fSelectMaxChi2PerITScluster < 1000 )
  {
    AliDebug(2,"Pass ITS chi2/ncls cut");

    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [3]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[3]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [3]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[3]->Fill(track->Pt(),cen);
    }
  }
  
  //
  // TPC related cuts
  // Not much sense to use with ITS only
  //

//  printf("TPC - n clus %d >= %d - chi2 %f max chi2/ncls %f\n",
//         aodtrack->GetTPCNcls(),fSelectMinTPCclusters,
//         aodtrack->GetTPCchi2(),fSelectMaxChi2PerTPCcluster);

  Int_t nTPCcls = aodtrack->GetTPCNcls();
  if ( fSelectMinTPCclusters > 0 )
  {
    if ( nTPCcls < fSelectMinTPCclusters  )
    {
      return kFALSE;
    }

    AliDebug(2,"Pass n TPC cluster cut");

    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [4]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[4]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [4]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[4]->Fill(track->Pt(),cen);
    }
  }

  //
  Float_t chi2PerTPCcluster = 0.;
  if ( nTPCcls > 0 ) chi2PerTPCcluster = aodtrack->GetTPCchi2()/nTPCcls;

  if ( chi2PerTPCcluster > fSelectMaxChi2PerTPCcluster )
  {
    return kFALSE;
  }

  if ( fSelectMaxChi2PerTPCcluster < 1000 )
  {
    AliDebug(2,"Pass TPC chi2/ncls cut");

    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [5]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[5]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [5]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[5]->Fill(track->Pt(),cen);
    }
  }

  //
  if ( fSelectFractionTPCSharedClusters )
  {
    Double_t frac = 0;
    Float_t ncls  = Float_t(aodtrack->GetTPCncls ());
    Float_t nclsS = Float_t(aodtrack->GetTPCnclsS());
    if (  ncls> 0 )  frac =  nclsS / ncls ;
    
    if ( frac > fCutTPCSharedClustersFraction )
    {
      AliDebug(2,Form("\t Reject track, shared cluster fraction %f > %f",frac, fCutTPCSharedClustersFraction));
      
      return kFALSE ;
    }
    
    AliDebug(2,"Pass TPC shared cluster cut");
    
    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [6]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[6]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [6]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[6]->Fill(track->Pt(),cen);
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
    
    AliDebug(2,"\t accepted track!");
    
    if ( !IsHistoCentDependentOn() ) fhCTSAODTrackCutsPt   [7]->Fill(aodtrack->Pt());
    else                             fhCTSAODTrackCutsPtCen[7]->Fill(aodtrack->Pt(),cen);

    if ( fillEmbedSignalTrack )
    {
      if ( !fHistoCentDependent ) fhCTSAODTrackCutsPt         [7]->Fill(track->Pt());
      else                        fhCTSAODTrackCutsPtCenSignal[7]->Fill(track->Pt(),cen);
    }
  }
  
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

