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
#include "AliGenEventHeader.h"

ClassImp(AliAnalysisTaskSVtaskMCFilter)


AliAnalysisTaskSVtaskMCFilter::AliAnalysisTaskSVtaskMCFilter() :
fFilteredTracksArray(0x0),
fFilteredClustersArray(0x0),
fFilteredTracksName("mytracks"),
fFilteredClustersName("myclusters"),
fInputTracksName("tracks"),
fInputClustersName("caloClusters"),
fFilterType(1),
fAodEvent(0x0),
fMCHeader(0x0),
fMCPartArray(0x0)
{
   //constructor


}
//_______________________________________________________________________________
AliAnalysisTaskSVtaskMCFilter::AliAnalysisTaskSVtaskMCFilter(const char *name) :
AliAnalysisTaskEmcal(name),
fFilteredTracksArray(0x0),
fFilteredClustersArray(0x0),
fFilteredTracksName("mytracks"),
fFilteredClustersName("myclusters"),
fInputTracksName("tracks"),
fInputClustersName("caloClusters"),
fFilterType(1),
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
//
//_______________________________________________________________________________
AliAnalysisTaskSVtaskMCFilter*  AliAnalysisTaskSVtaskMCFilter::AddTaskSVtaskMCFilter(const char* trkcontname, const char* outtrk,  Bool_t fFilterTracks, const char* clscontname, const char* outcls, Bool_t fFilterClusters){

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
     ::Error("AliAnalysisTaskSVtaskMCFilter.cxx", "No analysis manager to connect to.");
     return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()){
      ::Error("AddSVtaskMCFilter", "This task requires an input event handler");
      return NULL;
   }

   TString taskFiltername;
   if(fFilterClusters) taskFiltername = Form("SVFilterTask_%s-%s_to_%s-%s", trkcontname, clscontname, outtrk, outcls);
   else taskFiltername = Form("SVFilterTask_%s_to_%s", trkcontname, outtrk);

   AliAnalysisTaskSVtaskMCFilter* task = (AliAnalysisTaskSVtaskMCFilter*)mgr->GetTask(taskFiltername.Data());
   if(task){
      ::Info("AddTaskSVtaskMCFilter", "Task %s already exist, continue", taskFiltername.Data());
      return task;
   }else{
      ::Info("AddTaskSVtaskMCFilter", "Creating the task");

      // create the task
      task = new AliAnalysisTaskSVtaskMCFilter(taskFiltername.Data());
      AliTrackContainer*    trkCont;
      AliParticleContainer* mcpartCont;
      AliClusterContainer*  clsCont;
      if(fFilterTracks){
         trkCont = task->AddTrackContainer(trkcontname);//for data, and reco MC
         if(fFilterClusters) clsCont = task->AddClusterContainer(clscontname);
      }else{
         mcpartCont= task->AddMCParticleContainer(trkcontname);
      }
      mgr->AddTask(task);
   }

   Int_t filterparam;
   if(fFilterTracks){
     if(fFilterClusters) filterparam = 2;
     else filterparam = 1;
   }else filterparam = 0;

   task->SetInputTracksName(trkcontname);
   task->SetInputClustersName(clscontname);
   task->SetFilteredTracksName(outtrk);
   task->SetFilteredClustersName(outcls);
   task->SetFilterType(filterparam);
   // ------ input data ------
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   mgr->ConnectInput(task, 0, cinput);

   ::Info("AliAnalysisTaskSVtaskMCFilter", "Input and Output connected to the manager");
   return task;

}
//______________________________________________________________________________
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

  // add the filtered particles or tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fFilteredTracksName))) {
     if(fFilterType != 0){
        fFilteredTracksArray = new TClonesArray("AliAODTrack");
     }else{
        fFilteredTracksArray = new TClonesArray("AliAODMCParticle");
     }
     fFilteredTracksArray->SetName(fFilteredTracksName);
     ::Info("AliAnalysisTaskSVtaskMCFilter::ExecOnce", "particle or track collection with name '%s' has been added to the event.", fFilteredTracksName.Data());
     InputEvent()->AddObject(fFilteredTracksArray);
  }else{
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fFilteredTracksName.Data()));
    return;
  }

  // add the filtered clusters to event if not yet there and if needed
  if(fFilterType == 2){
    if (!(InputEvent()->FindListObject(fFilteredClustersName))) {
      fFilteredClustersArray = new TClonesArray("AliAODCaloCluster");
      fFilteredClustersArray->SetName(fFilteredClustersName);
      ::Info("AliAnalysisTaskSVtaskMCFilter::ExecOnce", "cluster collection with name '%s' has been added to the event.", fFilteredClustersName.Data());
      InputEvent()->AddObject(fFilteredClustersArray);
    }else{
      AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fFilteredClustersName.Data()));
      return;
    }
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
    if(fFilterType == 2) fFilteredClustersArray->Clear();

   // fix for temporary bug in ESDfilter
   // the AODs with null vertex pointer didn't pass the PhysSel
   if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return kFALSE;

   TString nameGen;
   Int_t lab, mother;
   AliAODMCParticle *mcpart = NULL;
   Int_t n = fFilteredTracksArray->GetEntriesFast();

   if (fFilterType != 0) { //detector level tracks
     AliAODTrack* track = 0;
     AliTrackContainer * tracks = dynamic_cast<AliTrackContainer *>(GetParticleContainer(fInputTracksName.Data()));
     // Iterable approach (using C++11)  loop over detector level tracks
     for(auto trackIterator : tracks->accepted_momentum()){

        track = dynamic_cast<AliAODTrack *>(trackIterator.second);  // Get the full track
        if(track){
           lab     = TMath::Abs(track->GetLabel());
           nameGen = GetGenerator(lab,fMCHeader);

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
              nameGen = GetGenerator(mother,fMCHeader); //change the name of generator to mother
           }
           if(!(nameGen.IsWhitespace() || nameGen.Contains("EPOS"))){
             new ((*fFilteredTracksArray)[n]) AliAODTrack(*track);
             n++;
           }
        }
     }

     // Filter clusters
     if(fFilterType == 2){
       TString nameGenClus;
       Int_t labClus, motherClus;
       AliAODMCParticle *mcpartClus = NULL;
       Int_t m = fFilteredClustersArray->GetEntriesFast();

       AliAODCaloCluster* cluster = 0;
       AliClusterContainer * clusters = dynamic_cast<AliClusterContainer *>(GetClusterContainer(fInputClustersName.Data()));
       // Iterable approach (using C++11)  loop over detector level clusters
       for(auto clusterIterator : clusters->accepted_momentum()){

          cluster = dynamic_cast<AliAODCaloCluster *>(clusterIterator.second);  // Get the full cluster
          if(cluster){
             labClus     = TMath::Abs(cluster->GetLabel());
             nameGenClus = AliVertexingHFUtils::GetGenerator(labClus,fMCHeader);

             while(nameGenClus.IsWhitespace()){    //if this cluster does not have any generator name execute this loop
                mcpartClus = (AliAODMCParticle*) fMCPartArray->At(labClus);  //find the corresponding MC particle
                if(!mcpartClus){
                   printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
                   break;
                }
                motherClus = mcpartClus->GetMother(); //find the index of the corresponding mother
                if(motherClus<0){   //the mother does not exist
                   printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
                   break;
                }
                labClus = motherClus;   // change label to the mother
                nameGenClus = AliVertexingHFUtils::GetGenerator(motherClus,fMCHeader); //change the name of generator to mother
             }
             if(!(nameGenClus.IsWhitespace() || nameGenClus.Contains("EPOS"))){
               new ((*fFilteredClustersArray)[m]) AliAODCaloCluster(*cluster);
               m++;
             }
          }
       }
     }

   }else{  //particle level particles
      AliAODMCParticle* mcParticle = NULL;
      AliAODMCParticle* initialmcParticle = NULL;
      AliMCParticleContainer *mcParticles = dynamic_cast<AliMCParticleContainer*>(GetMCParticleContainer(fInputTracksName.Data()));

      Int_t index=0;

      // Iterable approach (using C++11)  loop over detector level mcParticle
      for(auto trackIterator : mcParticles->all_momentum()){
	 initialmcParticle = dynamic_cast<AliAODMCParticle*>(trackIterator.second);  // Get the full track
         if(initialmcParticle){
            lab     = index;//TMath::Abs(initialmcParticle->GetLabel());
            nameGen = GetGenerator(lab,fMCHeader);

            while(nameGen.IsWhitespace()){    //if this particle does not have any generator name execute this loop
               mcParticle = (AliAODMCParticle*) fMCPartArray->At(lab);  //find the corresponding MC particle
               if(!mcParticle){
                  printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
                  break;
               }

               mother = mcParticle->GetMother(); //find the index of the corresponding mother
               if(mother<0){   //the mother does not exist
                  printf("AliAnalysisTaskMultCheck::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
                  break;
               }
               lab = mother;   // change label to the mother
               nameGen = GetGenerator(mother,fMCHeader); //change the name of generator to mother
            }

            if(!(nameGen.IsWhitespace() || nameGen.Contains("EPOS"))){
               new ((*fFilteredTracksArray)[n]) AliAODMCParticle(*initialmcParticle);
               n++;
            }
         }
         index++;
      }
   }


   return kTRUE;
}

//______________________________________________________________________
//Implementation originally from AliVertexingHFUtils. Moved to here due to library problems
TString AliAnalysisTaskSVtaskMCFilter::GetGenerator(Int_t label, AliAODMCHeader* header){
    /// get the name of the generator that produced a given particle

    Int_t nsumpart=0;
    TList *lh=header->GetCocktailHeaders();
    Int_t nh=lh->GetEntries();
    for(Int_t i=0;i<nh;i++){
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        Int_t npart=gh->NProduced();
        if(label>=nsumpart && label<(nsumpart+npart)) return genname;
        nsumpart+=npart;
    }
    TString empty="";
    return empty;
}
