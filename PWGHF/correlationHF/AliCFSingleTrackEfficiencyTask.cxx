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

/*__|______________________________________________________________________________|
 |                              -----Info(i)-----                                  |
 |  AliCFAnalysisTask which provides standard way of calculating single track      |
 |  efficiency between different steps of the procedure, ouptut of the task is a   |
 |  AliCFContainer from which the efficienciy can be calculated                    |
 |                                                                                 |
 |   ESDs<-->AODs (ON/OFF)                                                         |
 |                                                                                 |
 |                                                      Authors:                   |
 |                                                                                 |
 |_____________________________________________________________________________|___*/




#include "AliCFSingleTrackEfficiencyTask.h"
#include "AliSingleTrackEffCuts.h"

#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "TChain.h"

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"

#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"

#include "AliSingleTrackEffCuts.h"

ClassImp(AliCFSingleTrackEfficiencyTask)

//__________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask() :
  fReadTPCTracks(0),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(0x0),
  fTriggerMask(AliVEvent::kAny),
  fMCCuts(0x0),
  fSetFilterBit(kFALSE),
  fbit(0),
  fHistEventsProcessed(0x0)

{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const Char_t* name,AliESDtrackCuts *trackcuts, AliSingleTrackEffCuts * mccuts) :
  AliAnalysisTaskSE(name),
  fReadTPCTracks(0),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(trackcuts),
  fTriggerMask(AliVEvent::kAny),
  fMCCuts(mccuts),
  fSetFilterBit(kFALSE),
  fbit(0),
  fHistEventsProcessed(0x0)

{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFSingleTrackEfficiencyTask","Calling Constructor");

  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
  DefineOutput(3,TList::Class());
}

//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask& AliCFSingleTrackEfficiencyTask::operator=(const AliCFSingleTrackEfficiencyTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fReadTPCTracks = c.fReadTPCTracks ;
    fReadAODData = c.fReadAODData ;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList ;
    if(c.fTrackCuts) { delete fTrackCuts; fTrackCuts = new AliESDtrackCuts(*(c.fTrackCuts)); }
    fTriggerMask = c.fTriggerMask;
    if(c.fMCCuts) { delete fMCCuts; fMCCuts = new AliSingleTrackEffCuts(*(c.fMCCuts)); }

    fSetFilterBit  = c.fSetFilterBit;
    fbit = c.fbit ;
    fHistEventsProcessed = c.fHistEventsProcessed;



  }
  return *this;
}

//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const AliCFSingleTrackEfficiencyTask& c) :
  AliAnalysisTaskSE(c),
  fReadTPCTracks(c.fReadTPCTracks),
  fReadAODData(c.fReadAODData),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fTrackCuts(c.fTrackCuts),
  fTriggerMask(c.fTriggerMask),
  fMCCuts(c.fMCCuts),
  fSetFilterBit(c.fSetFilterBit),
  fbit(c.fbit),
  fHistEventsProcessed(c.fHistEventsProcessed)

{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::~AliCFSingleTrackEfficiencyTask() {
  //
  //destructor
  //
  Info("~AliCFSingleTrackEfficiencyTask","Calling Destructor");

  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
  if (fTrackCuts)           delete fTrackCuts;
  if (fMCCuts)              delete fMCCuts;


}

//_________________________________________________
void AliCFSingleTrackEfficiencyTask::Init()
{
  if(!fMCCuts) {
    AliFatal(" MC Cuts not defined");
    return;
  }
  if(!fTrackCuts) {
    AliFatal(" Track Cuts not defined");
    return;
  }
}

//_________________________________________________
void AliCFSingleTrackEfficiencyTask::UserExec(Option_t *)
{

       Info("UserExec","") ;

       AliVEvent*    fEvent = fInputEvent ;
       AliVParticle* track;
       
       if (!fInputEvent) {
	 AliFatal("NO EVENT FOUND!");
	 return;
       }

       if (!fMCEvent) {
	 AliFatal("NO MC INFO FOUND");
	 return ;
       }
       

       fHistEventsProcessed->Fill(0.5); // # of Event proceed        
       Bool_t IsEventMCSelected = kFALSE;
       Bool_t isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");
       

       // Step 0A. MC Gen Event Selection for ESDs/ AODs
       if(isAOD) {
	 if(!fEvent && AODEvent() && IsStandardAOD()) {
	   // In case there is an AOD handler writing a standard AOD, use the AOD 
	   // event in memory rather than the input (ESD) event.    
	   fEvent = dynamic_cast<AliAODEvent*> (AODEvent());
	 }
	 IsEventMCSelected = fMCCuts->IsMCEventSelected(fEvent);//AODs
       } else {
	 IsEventMCSelected = fMCCuts->IsMCEventSelected(fMCEvent);//ESDs
       }
       
       

       //pass the evt info to the cuts that need it 
       if(!isAOD) fCFManager->SetMCEventInfo (fMCEvent);
       else fCFManager->SetMCEventInfo (fEvent);
       fCFManager->SetRecEventInfo(fEvent);
       
       
       if(!IsEventMCSelected) {
	 AliDebug(3,"MC event not passing the quality criteria \n");
	 PostData(1,fHistEventsProcessed) ;
	 PostData(2,fCFManager->GetParticleContainer()) ;
	 PostData(3,fQAHistList) ;
	 return;
       }
       
       fHistEventsProcessed->Fill(1.5); // # of Event after passing MC cuts
       
       
       //Filling of MC generated particle -> below function 
       if(!isAOD) CheckESDParticles();
       else CheckAODParticles();
       
       // Step 0B. MC Reconstruction Event Selection for ESDs/ AODs
       Bool_t isRecoEventOk = fMCCuts->IsRecoEventSelected(fEvent);
     
       Double_t containerInput[5] ;
       Double_t containerInputMC[5] ; // for true pt
     
       if(isRecoEventOk) {
	 
	 fHistEventsProcessed->Fill(2.5); // # of Event after passing all cuts
	 const AliVVertex *vertex = fEvent->GetPrimaryVertex();
	 containerInput[4] = vertex->GetZ(); // Z Vertex of Event
	 containerInputMC[4]  =  containerInput[4];
	 //cout << "Z vtx of current event is @ Event " <<  containerInputMC[4] << endl; 
	 
	 // Reco tracks track loop
	 for (Int_t iTrack = 0; iTrack<fEvent->GetNumberOfTracks(); iTrack++) {

	   track = fEvent->GetTrack(iTrack);
	   
	   Double_t mom[3];
	   track->PxPyPz(mom);
	   Double_t pt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	   containerInput[0] = pt ;
	   containerInput[1] = track->Eta();
	   containerInput[2] = track->Phi() ;
	   containerInput[3] = track->Theta() ;
	   
	   // Step 4. Track that are recostructed and filling
	   fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;
	   
 
	   // Step 5. Track that are recostructed and +Kine acceptance filling
	   if (!fMCCuts->IsRecoParticleKineAcceptance(track)) continue;
	   fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoKineCuts) ;


	   // is track associated to particle ? if yes + implimenting the physical primary..
	   Int_t label = track->GetLabel();
	   if (label<0) {
	     AliDebug(3,"Particle not matching MC label \n");
	     continue;
	   }

	   // Step 6. Track that are recostructed + true pt and filling
	   AliVParticle *mcPart  = (AliVParticle*)fMCEvent->GetTrack(label);
	   containerInputMC[0] = (Float_t)mcPart->Pt();
	   containerInputMC[1] = mcPart->Eta() ;
	   containerInputMC[2] = mcPart->Phi() ;
	   containerInputMC[3] = mcPart->Theta() ;
 
	   if (!fMCCuts->IsMCParticleGenerated(mcPart)) continue;

	 	   
	    // for filter bit selection
	   AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(track);
	   if(isAOD && fSetFilterBit) if (!aodTrack->TestFilterMask(BIT(fbit))) continue;
	   Bool_t isESDtrack = track->IsA()->InheritsFrom("AliESDtrack");
	   AliESDtrack *tmptrack;
	   if(isESDtrack) {
	     tmptrack = dynamic_cast<AliESDtrack*>(track);
	   }
	   else {
	     tmptrack = ConvertTrack(aodTrack); //AODs
	   }
	   
	   // exclude global constrained and TPC only tracks (filter bits 128 and 512)
	   Int_t id = tmptrack->GetID();
	   if(isAOD && id<0) {
	     AliDebug(3,"Track removed bc corresponts to either filter bit 128 or 512 (TPC only tracks)\n");
	     continue;
	   }

         
	   // Step 7. Track that are recostructed + Quality + Kine criteria filling
	   if( fTrackCuts->IsSelected(tmptrack) ){
	     fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepReconstructedMC);
	     fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoQualityCuts);
	   }else AliDebug(3,"Reconstructed track not passing quality criteria\n");

	   if(isAOD) delete tmptrack;
	  
	 } 


       }
       else AliDebug(3,"Event not passing quality criteria\n");
       
  
       // PostData(0) is taken care of by AliAnalysisTaskSE 
       PostData(1,fHistEventsProcessed) ;
       PostData(2,fCFManager->GetParticleContainer()) ;
       PostData(3,fQAHistList) ;
       
  return;

}


//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::Terminate(Option_t*)
{

       Info("Terminate","");
       AliAnalysisTaskSE::Terminate();

       /*
       //draw some example plots....
       AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
       
       TH1D* h00 =   cont->ShowProjection(0,0) ;
       TH1D* h01 =   cont->ShowProjection(0,1) ;
       TH1D* h02 =   cont->ShowProjection(0,2) ;
       TH1D* h03 =   cont->ShowProjection(0,3) ;
       TH1D* h04 =   cont->ShowProjection(0,4) ;
       TH1D* h05 =   cont->ShowProjection(0,5) ;
       TH1D* h06 =   cont->ShowProjection(0,6) ;
       
       h00->SetMarkerStyle(23) ;
       h01->SetMarkerStyle(24) ;
       h02->SetMarkerStyle(25) ;
       h03->SetMarkerStyle(26) ;
       h04->SetMarkerStyle(27) ;
       h05->SetMarkerStyle(28) ;
       h06->SetMarkerStyle(28) ;
       
       
       
       TCanvas * c =new TCanvas("c","",1400,800);
       c->Divide(4,2);
       
       c->cd(1);
       h00->Draw("p");
       c->cd(2);
       h01->Draw("p");
       c->cd(3);
       h02->Draw("p");
       c->cd(4);
       h03->Draw("p");
       c->cd(5);
       h04->Draw("p");
       c->cd(6);
       h05->Draw("p");
       c->cd(7);
       h06->Draw("p");
       
       c->SaveAs("plots.eps");
       */
       
}


//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::UserCreateOutputObjects() {
  
       Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());

       //slot #1
       OpenFile(1);
       
       fHistEventsProcessed = new TH1I("fHistEventsProcessed","fHistEventsProcessed",3,0,3) ;
       fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"All events");
       fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Good MC events");
       fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Good Reconstructed events");
       
       PostData(1,fHistEventsProcessed) ;
       PostData(2,fCFManager->GetParticleContainer()) ;
       PostData(3,fQAHistList) ;
       
       return;
}


//__________| Function to fill MC Gen Particle in Diff Steps with EDSs |___________________
void AliCFSingleTrackEfficiencyTask::CheckESDParticles(){
  
       if (!fMCEvent) {
	 AliFatal("NO MC INFO FOUND");
	 return;
       }
       
       
       Double_t containerInput[5] ; //number of variables
       
       TArrayF vtxPos(3); // for Z vtx
       AliGenEventHeader *genHeader;
       genHeader = fMCEvent->GenEventHeader();
       genHeader->PrimaryVertex(vtxPos);
       containerInput[4]  = vtxPos[2];
       
       //loop on the MC Gen Partiles
       for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
	 
	 AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);  
	 containerInput[0] = (Float_t)mcPart->Pt(); 
	 containerInput[1] = mcPart->Eta() ;
	 containerInput[2] = mcPart->Phi() ;
	 containerInput[3] = mcPart->Theta() ;
	 
	 
	 // Step 1. Particle passing through Generation criteria and filling
	 if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
	   AliDebug(3,"MC Particle not passing through genetations criteria\n");
	   continue;
	 }
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);
	 
     
	 
	 // Step 2. Particle passing through Kinematic criteria and filling
	 if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
	   AliDebug(3,"MC Particle not in the kine acceptance\n");
	   continue;
	 }
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);
	 
	 // Step 3. Particle passing through Track ref criteria and filling
	 // did leave signal (enough clusters) on the detector
	 if( !fMCCuts->IsMCParticleInReconstructable(mcPart) ) {
	   AliDebug(3,"MC Particle not in the reconstructible\n");
	   continue;
	 }
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCAccpCut);
	 
	 
       }// end of particle loop
       

       return;
              
}

//__________| Function to fill MC Gen Particle in Diff Steps with AODs |___________________
void AliCFSingleTrackEfficiencyTask::CheckAODParticles(){

       if (!fInputEvent) {
	 AliFatal("NO EVENT FOUND!");
	 return;
       }

       AliVEvent* event = fInputEvent ;    
       if(!event && AODEvent() && IsStandardAOD()) {
	 // In case there is an AOD handler writing a standard AOD, use the AOD 
	 // event in memory rather than the input (ESD) event.    
	 event = dynamic_cast<AliAODEvent*> (AODEvent());
       }
       

       TClonesArray* mcArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
       if (!mcArray) {
	 AliError("Could not find Monte-Carlo in AOD");
	 return;
       }
       
       AliAODMCHeader *mcHeader;
       mcHeader = dynamic_cast<AliAODMCHeader*>(event->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
       if (!mcHeader) {
	 AliError("Could not find MC Header in AOD");
	 //return kFALSE;
       }
       
       
       Double_t containerInput[5] ;
       containerInput[4]  = mcHeader->GetVtxZ();
       // loop over AOD MC-Particles 
       for (Int_t ipart=0; ipart<mcArray->GetEntriesFast(); ipart++) { 
	 
	 AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(ipart));
	 containerInput[0] = (Float_t)mcPart->Pt();
	 containerInput[1] = mcPart->Eta() ;
	 containerInput[2] = mcPart->Phi() ;
	 containerInput[3] = mcPart->Theta() ;

	 if (!mcPart){
	   AliError("Failed casting particle from MC array!, Skipping particle");
	   continue;
	 }
	 

	
	 // Step 1. Particle passing through Generation criteria and filling
	 if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
	   AliDebug(3,"MC Particle not passing quality criteria\n");
	   continue;
	 }
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);
	 

	 // Step 2. Particle passing through Kinematic criteria and filling
	 if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
	   AliDebug(3,"MC Particle not in the acceptance\n");
	   continue;
	 }
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);
	 

	 // Step 3. Particle passing through Track Ref criteria and filling
	 // but no info available for Track ref in AOD fillng same as above
	 fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCAccpCut);
     

       }
       
       return;

}


//___________________________________________________________________________
AliESDtrack * AliCFSingleTrackEfficiencyTask::ConvertTrack(AliAODTrack *track)
{

       AliAODEvent *aodEvent = (AliAODEvent*)fInputEvent;
       const AliAODVertex *primary = aodEvent->GetPrimaryVertex();
       Double_t pos[3],cov[6];
       primary->GetXYZ(pos);
       primary->GetCovarianceMatrix(cov);
       const AliESDVertex vESD(pos,cov,100.,100);
       
       AliESDtrack *esdTrack =  new AliESDtrack(track);
       // set the TPC cluster info
       esdTrack->SetTPCClusterMap(track->GetTPCClusterMap());
       esdTrack->SetTPCSharedMap(track->GetTPCSharedMap());
       esdTrack->SetTPCPointsF(track->GetTPCNclsF());
       // needed to calculate the impact parameters
       esdTrack->RelateToVertex(&vESD,0.,3.); 

       //  std::cout << " primary vtx "<< primary << std::endl;
       //  std::cout << " esdvtx "<< vESD.GetName() << std::endl;

       //  std::cout<< " esdtrack pt "<< esdTrack.Pt() << " and status " << esdTrack.GetStatus() <<endl;
       //  std::cout << " aod track "<< track<< " and status " << track->GetStatus() << std::endl;
       //  std::cout << " esd track "<< esdTrack<< " and status " << esdTrack->GetStatus() << std::endl;
       
       return esdTrack;


}

