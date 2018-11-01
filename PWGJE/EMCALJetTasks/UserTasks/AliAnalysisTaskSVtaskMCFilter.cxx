/**************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
//
// Task for filtering out EPOS part of the simulated detector level events   
//
//-----------------------------------------------------------------------
// Authors: Filip Krizek

//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TList.h>
#include <TROOT.h>
#include <TH3F.h>
#include <TRandom3.h>

#include "AliLog.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliStack.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSVtaskMCFilter.h"

ClassImp(AliAnalysisTaskSVtaskMCFilter)


AliAnalysisTaskSVtaskMCFilter::AliAnalysisTaskSVtaskMCFilter() :
fFilteredTracksArray(0x0),
fFilteredTracksName("mytracks"),
fInputTracksName("tracks"),
fAodEvent(0x0),
fMCHeader(0x0),
fMCPartArray(0x0)
{
   //constructor


}
//_______________________________________________________________________________
AliAnalysisTaskSVtaskMCFilter::AliAnalysisTaskSVtaskMCFilter(const char *name) :
fFilteredTracksArray(0x0),
fFilteredTracksName("mytracks"),
fInputTracksName("tracks"),
fAodEvent(0x0),
fMCHeader(0x0),
fMCPartArray(0x0)
{
   //constructor
}
//_______________________________________________________________________________
AliAnalysisTaskSVtaskMCFilter::~AliAnalysisTaskSVtaskMCFilter()
{
   //destructor
}
//_______________________________________________________________________________
void AliAnalysisTaskSVtaskMCFilter::UserCreateOutputObjects(){
 //dummy
  return;
}
//_______________________________________________________________________________
void AliAnalysisTaskSVtaskMCFilter::ExecOnce()
{

   //
   // To be executed only once, for the first event
   //

   AliDebug(2, "Entering ExecOnce()");

   // Load the event
   fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

   if(!fAodEvent){
       if (AODEvent() && IsStandardAOD()) {

          // In case there is an AOD handler writing a standard AOD, use the AOD
          // event in memory rather than the input (ESD) event.
          fAodEvent = dynamic_cast<AliAODEvent*>(AODEvent());
       }
   } 

   fMCPartArray = dynamic_cast<TClonesArray*>(fAodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
   if(!fMCPartArray) {
      AliError(Form("MC particles not found! Task '%s' will be disabled!", GetName()));
      return;
   }

   fMCHeader = (AliAODMCHeader*)fAodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
   if(!fMCHeader) {
      AliError(Form("MC header not found! Task '%s' will be disabled!", GetName()));
      return;
   }

  // add the filtered tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fFilteredTracksName))) {
     fFilteredTracksArray = new TClonesArray("AliAODTrack");
     fFilteredTracksArray->SetName(fFilteredTracksName);
     ::Info("AliAnalysisTaskSVtaskMCFilter::ExecOnce", "track collection with name '%s' has been added to the event.", fFilteredTracksName.Data());
     InputEvent()->AddObject(fFilteredTracksArray);
  }
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fFilteredTracksName.Data()));
    return;
  }
  
   AliAnalysisTaskEmcal::ExecOnce();
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskSVtaskMCFilter::Run()
{

   //
   // Analysis execution
   //

    //fFilteredTracksArray->Delete();
    fFilteredTracksArray->Clear();

   // fix for temporary bug in ESDfilter
   // the AODs with null vertex pointer didn't pass the PhysSel
   if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return kFALSE;

 

   GetTrackContainer(0)->ResetCurrentID();
   AliAODTrack* track = 0;
   TString nameGen;
   Int_t lab, mother;
   AliAODMCParticle *mcpart = NULL;
   Int_t n = fFilteredTracksArray->GetEntriesFast();

   AliTrackContainer * tracks = dynamic_cast<AliTrackContainer *>(GetParticleContainer(fInputTracksName.Data()));
   // Iterable approach (using C++11)  loop over detector level tracks
   for(auto trackIterator : tracks->accepted_momentum()){

      track = dynamic_cast<AliAODTrack *>(trackIterator.second);  // Get the full track

      lab     = TMath::Abs(track->GetLabel());
      nameGen = AliVertexingHFUtils::GetGenerator(lab,fMCHeader);

      while(nameGen.IsWhitespace()){    //if this particle does not have any generator name execute this loop
         mcpart = (AliAODMCParticle*) fMCPartArray->At(lab);  //find the corresponding MC particle
         if(!mcpart){
            printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
            break;
         }

         mother = mcpart->GetMother(); //find the index of the corresponding mother
         if(mother<0){   //the mother does not exist
            printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
            break;
         }

         lab = mother;   // change label to the mother   
         nameGen = AliVertexingHFUtils::GetGenerator(mother,fMCHeader); //change the name of generator to mother
      }

      if(!(nameGen.IsWhitespace() || nameGen.Contains("EPOS"))){
          new ((*fFilteredTracksArray)[n]) AliAODTrack(*track);
          n++;

      }
   }
}
