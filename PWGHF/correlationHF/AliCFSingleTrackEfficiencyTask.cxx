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

#include "AliSelectNonHFE.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelection.h"

#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "TChain.h"

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliHFEtools.h"


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
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliDxHFEToolsMC.h"

#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliVertexingHFUtils.h"

#include <memory>

ClassImp(AliCFSingleTrackEfficiencyTask)

//__________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const char* opt) :
fOption(opt),
  fReadTPCTracks(0),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(0x0),
  fTriggerMask(AliVEvent::kAny),
  fMCCuts(0x0),
  fSelNHFE(NULL),
  fElectrons(NULL),
  fElectronsKine(NULL),
  fVertUtil(0),
  fSetFilterBit(kFALSE),
  fbit(0),
  fMinNclsTPCPID(80),
  fRequireTOF(kTRUE),
  fMinRatioTPCcluster(0.6),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fTOFnSigma(0),
//  fParticleIDforPID(AliPID::kElectron),
  fUsePID(kTRUE),
  fUseTPCPID(kTRUE),
  fUseTOFPID(kTRUE),
  fMaxPtForTOFPID(999),
  fUseTOFonlyWhenPresent(kFALSE),
  fInvMassLow(0.10),
  fMaxRadius(999),
  fOriginMotherReco(-1),
  fOriginMotherKine(-1),
  fSelectElSource(kAll),
  fUseGenerator(kFALSE),
  fHistEventsProcessed(0x0),
  fElectronPt(NULL),
  fElectronPtStart(NULL),
  fHistInvMassLS(NULL),
  fHistInvMassULS(NULL),
  fhGenerator(NULL),
  fhOriginKine(NULL),
  fhOriginReco(NULL),
  fhElGenerator(NULL)


{
  //
  //Default ctor
  //

}
//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const char* opt,const Char_t* name,AliESDtrackCuts *trackcuts, AliSingleTrackEffCuts * mccuts) :
  AliAnalysisTaskSE(name),
  fOption(opt),
  fReadTPCTracks(0),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(trackcuts),
  fTriggerMask(AliVEvent::kAny),
  fMCCuts(mccuts),
  fSelNHFE(NULL),
  fElectrons(NULL),
  fElectronsKine(NULL),
  fVertUtil(0),
  fSetFilterBit(kFALSE),
  fbit(0),
  fMinNclsTPCPID(80),
  fRequireTOF(kTRUE),
  fMinRatioTPCcluster(0.6),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fTOFnSigma(0),
  //  fParticleIDforPID(AliPID::kElectron),
  fUsePID(kTRUE),
  fUseTPCPID(kTRUE),
  fUseTOFPID(kTRUE),
  fMaxPtForTOFPID(999),
  fUseTOFonlyWhenPresent(kFALSE),
  fInvMassLow(0.10),
  fMaxRadius(999),
  fOriginMotherReco(-1),
  fOriginMotherKine(-1),
  fSelectElSource(kAll),
  fUseGenerator(kFALSE),
  fHistEventsProcessed(0x0),
  fElectronPt(NULL),
  fElectronPtStart(NULL),
  fHistInvMassLS(NULL),
  fHistInvMassULS(NULL),
  fhGenerator(NULL),
  fhOriginKine(NULL),
  fhOriginReco(NULL),
  fhElGenerator(NULL)

{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFSingleTrackEfficiencyTask","Calling Constructor");
  fVertUtil=new AliVertexingHFUtils();

  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */

  DefineInput(0, TChain::Class());
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
    fSelNHFE=c.fSelNHFE;
    fElectrons=c.fElectrons;
    fElectronsKine=c.fElectronsKine;
    fVertUtil=c.fVertUtil;
    fSetFilterBit  = c.fSetFilterBit;
    fbit = c.fbit ;
    fMinNclsTPCPID=c.fMinNclsTPCPID;
    fRequireTOF=c.fRequireTOF;
    fMinRatioTPCcluster=c.fMinRatioTPCcluster;
    fMaxRadius=c.fMaxRadius;
    fOriginMotherReco=c.fOriginMotherReco;
    fOriginMotherKine=c.fOriginMotherKine;
    fSelectElSource=c.fSelectElSource;  
    fUseGenerator=c.fUseGenerator;

    fHistEventsProcessed = c.fHistEventsProcessed;
    fElectronPt=c.fElectronPt;
    fElectronPtStart=c.fElectronPtStart;
    fHistInvMassLS=c.fHistInvMassLS;
    fHistInvMassULS=c.fHistInvMassULS;
    fhGenerator=c.fhGenerator;
    fhOriginKine=c.fhOriginKine;
    fhOriginReco=c.fhOriginReco;
    fhElGenerator=c.fhElGenerator;

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
  fSelNHFE(c.fSelNHFE),
  fElectrons(c.fElectrons),
  fElectronsKine(c.fElectronsKine),
  fVertUtil(c.fVertUtil),
  fSetFilterBit(c.fSetFilterBit),
  fbit(c.fbit),
  fMinNclsTPCPID(c.fMinNclsTPCPID),
  fRequireTOF(c.fRequireTOF),
  fMinRatioTPCcluster(c.fMinRatioTPCcluster),
  fMaxRadius(c.fMaxRadius),
  fOriginMotherReco(c.fOriginMotherReco),
  fOriginMotherKine(c.fOriginMotherKine),
  fSelectElSource(c.fSelectElSource),
  fUseGenerator(c.fUseGenerator),

  fHistEventsProcessed(c.fHistEventsProcessed),
  fElectronPt(c.fElectronPt),
  fElectronPtStart(c.fElectronPtStart),
  fHistInvMassLS(c.fHistInvMassLS),
  fHistInvMassULS(c.fHistInvMassULS),
  fhGenerator(c.fhGenerator),
  fhOriginKine(c.fhOriginKine),
  fhOriginReco(c.fhOriginReco),
  fhElGenerator(c.fhElGenerator)

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
  if( fElectronPt)          delete fElectronPt;
  if( fElectronPtStart)     delete fElectronPtStart;
  if(fHistInvMassLS)        delete fHistInvMassLS;
  if(fHistInvMassULS)       delete fHistInvMassULS;
  if(fElectrons)            delete fElectrons;
  if(fElectronsKine)        delete fElectronsKine;
  if(fVertUtil)             delete fVertUtil;
  if(fhGenerator)           delete fhGenerator;
  if(fhOriginKine)          delete fhOriginKine;
  if(fhOriginReco)          delete fhOriginReco;
  if(fhElGenerator)         delete fhElGenerator;

  if(fSelNHFE){
    delete fSelNHFE;
    fSelNHFE=NULL;
  }



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
       
  //Info("UserExec","") ;
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
  
  AliPIDResponse *pidResponse=NULL;
  if(fUsePID){
    pidResponse= ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetPIDResponse();
    //    fpidResponse = fInputHandler->GetPIDResponse();
    if(!pidResponse){
      AliDebug(1, "Using default PID Response");
      pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
    }     
  }
  fSelNHFE->SetPIDresponse(pidResponse);

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
      fOriginMotherReco=-1;
      track = fEvent->GetTrack(iTrack);
      fElectronPtStart->Fill(track->Pt());

      AliAODTrack *aodtrack = NULL;
      if(isAOD) aodtrack= static_cast<AliAODTrack *>(track);
	   
      Double_t mom[3];
      track->PxPyPz(mom);
      Double_t pt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      containerInput[0] = pt ;
      containerInput[1] = track->Eta();
      containerInput[2] = track->Phi() ;
      containerInput[3] = track->Theta() ;

      /*
      int originvsGen=CheckBackgroundSource(track,const_cast<AliVEvent*>(fEvent));

      containerInput[5]=originvsGen;
      fhElGenerator->Fill(originvsGen);*/

      // Step 4. Track that are recostructed and filling
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;
	   
 
      // Step 5. Track that are recostructed and +Kine acceptance filling
      // Pt and eta
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
      containerInputMC[0] = mcPart->Pt();
      containerInputMC[1] = mcPart->Eta() ;
      containerInputMC[2] = mcPart->Phi() ;
      containerInputMC[3] = mcPart->Theta() ;
      
      if (!fMCCuts->IsMCParticleGenerated(mcPart)) continue;
      AliAODMCParticle *mcPart2=dynamic_cast<AliAODMCParticle*>(mcPart);
      
      if((aodtrack->Charge() != -1) && (aodtrack->Charge() != 1)) cout << "charge: " << aodtrack->Charge() << endl;

      Bool_t selected=kTRUE;


      int originvsGen=CheckBackgroundSource(track,const_cast<AliVEvent*>(fEvent));

      fhElGenerator->Fill(originvsGen);

      if(fUseGenerator){
	if(fSelectElSource==kHFPythia && originvsGen!=kHFPythia){
	  selected=kFALSE;
	}
	if(fSelectElSource==knonHFHijing && originvsGen!=knonHFHijing){
	  selected=kFALSE;
	}

	if(fSelectElSource==kConvElHijing && originvsGen!=kConvElHijing){
	  selected=kFALSE;
	}
      }
      else{

	if(fSelectElSource==kHF && originvsGen!=kHF){
	  selected=kFALSE;
	}
	if(fSelectElSource==knonHF && originvsGen!=knonHF){
	  selected=kFALSE;
	}

	if(fSelectElSource==kConvEl && originvsGen!=kConvEl){
	  selected=kFALSE;
	}
      }

      if(!selected) continue;
      
      fhOriginReco->Fill(fOriginMotherReco);
	
      Double_t x=mcPart2->Xv();
      Double_t y=mcPart2->Yv();
      double radius=TMath::Sqrt(x*x+y*y);
      if(radius > fMaxRadius) { continue;}
      
      
      //if(!selected) cout << "radius: " << radius << "    fMaxr: " << fMaxRadius << endl; 

	 	   
      // for filter bit selection
      AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(track);
      if(isAOD && fSetFilterBit) 
	if (!aodTrack->TestFilterMask(BIT(fbit))){
	  //	  cout << "cut due to filter bit " << fSetFilterBit << "  " <<fbit <<endl;
	  continue;
	}
      AliVTrack *vtrack=static_cast<AliVTrack *>(track);
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
	//cout << "Track removed" << endl;
	AliDebug(3,"Track removed bc corresponts to either filter bit 128 or 512 (TPC only tracks)\n");
	continue;
      }


      // Step 7. Track that are recostructed + Quality + Kine criteria filling
      if(! fTrackCuts->IsSelected(tmptrack) ){
	selected=kFALSE;
	continue;
      }
      if(selected){

	AliDebug(2,"Reconstructed track pass first quality criteria\n");
	//fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepReconstructedFirstTrackCutsMC);
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoFirstQualityCuts);
      }else {
	//cout <<"Not passing first " << endl;
	AliDebug(3,"Reconstructed track not passing first quality criteria\n");
      }


      Bool_t useTOFPID=kTRUE;
      if(track->Pt() > fMaxPtForTOFPID) useTOFPID=kFALSE;
      if(useTOFPID) 
	AliDebug(2,Form("Pt: %f, use CombinedPID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtForTOFPID));
      else 
	AliDebug(2,Form("Pt: %f, use only TPC PID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtForTOFPID));
      /*
	if(useTOFPID) printf(Form("Pt: %f, use CombinedPID (fMaxPtCombinedPID= %f)\n",track->Pt(),fMaxPtForTOFPID));
	if(fMinNclsTPCPID>0) cout <<"using TPC clusters for PID"<<endl;
	if(fMinRatioTPCcluster>0) cout <<"using ratio TPC clusters"<<endl;
	if(fRequireTOF) cout << "require hit in TOF" << endl;*/
      // DxHFE: Add here the extra cuts
      // 1. TPC PID clusters
      if(fMinNclsTPCPID>0 && selected){
	Int_t nclsTPCPID = tmptrack->GetTPCsignalN();
	if(nclsTPCPID<fMinNclsTPCPID){
	  //cout << "cut due to nr cls TPC PID " << endl;
	  AliDebug(2,Form("nlcTPCPID NOT selected - nrclusters=%d", nclsTPCPID));
	  selected=kFALSE;
	}
      }
 
      // 2. ratio TPC/findable
      if(fMinRatioTPCcluster>0 && selected){	     
	if(isAOD){
	  //AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
	  const TBits &clusterTPC = aodtrack->GetTPCClusterMap();
	  Int_t allclusters = clusterTPC.CountBits();
	  Double_t clusterRatio = vtrack->GetTPCNclsF() ? static_cast<Double_t>(allclusters)/static_cast<Double_t>(vtrack->GetTPCNclsF()) : 1.; // foundall/findable
	  AliDebug(2,Form("clusterRatio: %f  of all clusters: %d",clusterRatio,allclusters));
	  if(clusterRatio <= fMinRatioTPCcluster){
	    AliDebug(2,"clusterRatio NOT selected");
	    //cout << "cut due to clusterratio" << endl;
	    selected=kFALSE;
	  }
	}
      }
      
      // 3. TOF matching
      if(fRequireTOF && selected && useTOFPID){
	if(!(vtrack->GetStatus() & AliESDtrack::kTOFpid)){
	  if(fUseTOFonlyWhenPresent) 
	    useTOFPID=kFALSE;
	  else{
	    selected = kFALSE;
	    AliDebug(2,"Cut due to TOF requirement");
	  }
	}
      }
      if(selected){
	fElectronPt->Fill(tmptrack->Pt());
      }
	   
      if(selected){
	AliDebug(2,"Reconstructed track pass quality criteria\n");
	//fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepReconstructedMC);
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoQualityCuts);
      }else AliDebug(3,"Reconstructed track not passing quality criteria\n");
	   
      // PID requirement
      if(selected){
	//also check for pdg first????
	if(fUseTPCPID){
	  Float_t tpcNsigma = pidResponse->NumberOfSigmasTPC(vtrack, AliPID::kElectron); // change to particle
	  AliDebug(2, Form("Number of sigmas in TPC: %f", tpcNsigma));
	  //cout << "sigmaTPC " << tpcNsigma << "   " << fTPCnSigmaMin << " - " << fTPCnSigmaMax << endl;
	  if(tpcNsigma<fTPCnSigmaMin || tpcNsigma>fTPCnSigmaMax) { selected = false;}

	}

	if(useTOFPID && fUseTOFPID){

	  // if fMaxPtCombinedPID is set to lower than upper Ptlimit (10GeV/c), will separate
	  // PID into two regions: below fMaxptCombinedPID - both TPC and TOF, above only TPC
	  Float_t tofNsigma = pidResponse->NumberOfSigmasTOF(vtrack, AliPID::kElectron); //change to particle
	  AliDebug(2, Form("Number of sigmas in TOF: %f", tofNsigma));
	  // cout << "sigmaTOF " << tofNsigma << "   " << fTOFnSigma  << endl;
	  // for now: Assume symmetric cut for TOF
	  if(TMath::Abs(tofNsigma) > fTOFnSigma  && useTOFPID) {selected = false;}
	}

	if(selected){
	  //cout << "RECO: originvsGen: " << originvsGen << " kConvElHijing: " << kConvElHijing << endl;

	  // fill container for tracks
	  fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepRecoPIDMC);
	  fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPID);
	}

      }
 
      // invariant mass method - Not sure if needed here...
      if(selected){
	fSelNHFE->FindNonHFE(iTrack, aodtrack, const_cast<AliVEvent*>(fEvent));
	if(fSelNHFE->IsLS() || fSelNHFE->IsULS())
	  {
	    //Not selected
	    //	    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kINVMASS);
	    //cout <<"Not selected due to inv mass" << endl;
	    AliDebug(2,"Cut: Invmass");
	    selected=false;
	  }
	if(selected){

	  // fill container for tracks, with all electrons
	  //fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepRecoInvMassMC);
	  fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoInvMass);



	}
      }
      
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
int AliCFSingleTrackEfficiencyTask::CheckBackgroundSource(AliVParticle* track, const AliVEvent* pEvent,Bool_t useMCarray){

  // cout << "here " << endl;
  //cout <<  track->GetLabel() << endl;
  int res=0;
  if(useMCarray) res=fElectronsKine->CheckMC(track, const_cast<AliVEvent*>(pEvent));
  else res=fElectrons->CheckMC(track, const_cast<AliVEvent*>(pEvent));
  //cout << "here2 " << endl;
  int origin=-1;
  if(useMCarray){
    origin=fElectronsKine->GetOriginMother();
    fOriginMotherKine=origin;
  }
  else{
    origin=fElectrons->GetOriginMother();
    fOriginMotherReco=origin;
  }
  Bool_t isCharm=(origin==AliDxHFEToolsMC::kOriginCharm || 
		  origin==AliDxHFEToolsMC::kOriginGluonCharm);
  Bool_t isBeauty=(origin==AliDxHFEToolsMC::kOriginBeauty || 
		   origin==AliDxHFEToolsMC::kOriginGluonBeauty);

  //  if(origin==AliDxHFEToolsMC::kOriginGluonCharm || origin==AliDxHFEToolsMC::kOriginGluonBeauty) cout << "HELLO" << endl;

  Bool_t isConversion=(origin==AliDxHFEToolsMC::kNrOrginMother+2);
  //if(isConversion) cout << "Origin is conversion " << origin << endl;
  //else cout << "not conversion" << endl;
  TString nameGen;	
  int originvsGen=-2;
  Int_t generator=0;// adding per track: -1 hijing; 0-pythiaHF; 2-the rest  --> -2= pure hijing -> back ok (but check identity); 0= phytia,pythia-> ok for checking signal; reject the rest: -1 (Hij,pyt), 1 (hij, rest), 2 (pythia,rest) ,4 (rest,rest) 
      
  if(fUseGenerator){
    originvsGen=-1; 
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(pEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      //return kFALSE;
    }
    AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
    if(!aodtrack) {cout << "no AODtrack" << endl; return -1;}
	   
    if(useMCarray){
      TClonesArray* mcArray = dynamic_cast<TClonesArray*>(pEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
	AliError("Could not find Monte-Carlo in AOD");
	return -1;
      }
      //nameGen="";
      //fVertUtil->
      GetTrackPrimaryGenerator(aodtrack,mcHeader,mcArray,nameGen);
    }
    else{
      GetTrackPrimaryGenerator(aodtrack,mcHeader,fMCEvent,nameGen);
    }
  

    if(nameGen.Contains("ijing")){
      //cout << "hijing generator" << endl;
      generator=1;
    }
    else if(nameGen.Contains("ythia")){
      //cout << "pythia generator" << endl;
      generator=2;
    }
  }
  bool isHF=(isCharm || isBeauty);

  if(isHF){
    originvsGen=kHF;

    if(fUseGenerator){
      if(generator==0)
	originvsGen=kHFGen0;
      if(generator==1)
	originvsGen=kHFHijing;
      if(generator==2)
	originvsGen=kHFPythia;
    }
  }
  else{
    if(isConversion){
      originvsGen=kConvEl;
      if(fUseGenerator){
	if(generator==0)
	  originvsGen=kConvElGen0;
	if(generator==1)
	  originvsGen=kConvElHijing;
	if(generator==2)
	  originvsGen=kConvElPythia;
      }
      
    }
    else{
      originvsGen=knonHF;
      if(fUseGenerator){
	if(generator==0)
	  originvsGen=knonHFGen0;
	if(generator==1)
	  originvsGen=knonHFHijing;
	if(generator==2)
	  originvsGen=knonHFPythia;
      }      

    }
  }
  
  return originvsGen;
  
}

//_____________________________________________________________________
void AliCFSingleTrackEfficiencyTask::GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  // method to check if a track comes from a given generator

  Int_t lab=track->GetLabel();
  if(lab<0) return;
  nameGen=GetGenerator(lab,header);
  
  Int_t countControl=0;

  AliAODMCParticle *mcpart=(AliAODMCParticle*)track;

  while(nameGen.IsWhitespace()){
    if(countControl>0) {
      mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    }
    if(!mcpart){
      printf("IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=GetGenerator(mother,header);
    countControl++;
    /**
    if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
      printf("BREAK: Protection from infinite loop active\n");
      break;
      }*/
  }
  
  return;
}
//_____________________________________________________________________
void AliCFSingleTrackEfficiencyTask::GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,AliVEvent *arrayMC,TString &nameGen){

  // method to check if a track comes from a given generator

  Int_t lab=track->GetLabel();
  if(lab<0) return;
  nameGen=GetGenerator(lab,header);
  
  Int_t countControl=0;

  while(nameGen.IsWhitespace()){
    //    cout << "lab: " << lab << endl;
    AliAODMCParticle * mcpart= (AliAODMCParticle*)arrayMC->GetTrack(lab);
    if(!mcpart){
      printf("IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=GetGenerator(mother,header);
    countControl++;
    /*
    if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
      printf("BREAK: Protection from infinite loop active\n");
      break;
      }*/
  }
  
  return;
}
//______________________________________________________________________
TString AliCFSingleTrackEfficiencyTask::GetGenerator(Int_t label, AliAODMCHeader* header){
  // get the name of the generator that produced a given particle

  Int_t nsumpart=0;
  TString empty="";
  TList *lh=header->GetCocktailHeaders();
  if(!lh){
    cout << "no cocktail header" << endl;
    return empty;
  }
  Int_t nh=lh->GetEntries();
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
    TString genname=gh->GetName();
    Int_t npart=gh->NProduced();
    if(label>=nsumpart && label<(nsumpart+npart)) return genname;
    nsumpart+=npart;
  }

  return empty;
}

//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::Terminate(Option_t*)
{

  //  Info("Terminate","");
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

  ParseArguments(fOption.Data());

  fQAHistList = new TList;
  fQAHistList->SetOwner();
       

  //Initialization of invariant mass cut function (AliSelectNonHFE)
  fSelNHFE= new AliSelectNonHFE("IM","IM");

  // Cuts for associated track in Inv Mass cut
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts();
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetRequireSigmaToVertex(kTRUE);
  trackCuts->SetMaxChi2PerClusterTPC(4.0);
  trackCuts->SetMinNClustersTPC(80);
  trackCuts->SetPtRange(0.3,1e10);
  
  fSelNHFE->SetTrackCuts(-3, 3, trackCuts);
  fSelNHFE->SetInvariantMassCut(fInvMassLow);
  //  fSelNHFE->SetAlgorithm("KF");
  fSelNHFE->SetAODanalysis(kTRUE);

  // Invariant mass LS and ULS without cut
  fHistInvMassLS= new TH1F("fInvMassLS","Invariant mass LS",1000,0,0.5) ;
  fHistInvMassULS = new TH1F("fInvMassULS","Invariant mass ULS",1000,0,0.5) ;
  

  fSelNHFE->SetHistMass(fHistInvMassULS);
  fSelNHFE->SetHistMassBack(fHistInvMassLS);
  fQAHistList->Add(fHistInvMassLS);
  fQAHistList->Add(fHistInvMassULS);

  // Setting up the electron selection class for MC
  TString option="usekine"; 
  if(!fMCCuts->GetUseIsPhysicalPrimary()) option+=" notusePhysPrim";
  fElectronsKine= new AliDxHFEParticleSelectionMCEl(option);
  int result=fElectronsKine->Init();
  if (result<0) {
    AliFatal(Form("initialization of worker class instance fElectronsKine failed with error %d", result));
  }

  fElectrons= new AliDxHFEParticleSelectionMCEl();
  result=fElectrons->Init();
  if (result<0) {
    AliFatal(Form("initialization of worker class instance fElectrons failed with error %d", result));
  }

  fHistEventsProcessed = new TH1I("fHistEventsProcessed","fHistEventsProcessed",3,0,3) ;
  fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"All events");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Good MC events");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Good Reconstructed events");

  fElectronPt = new TH1F("fElectronPt","fElectronPt",100,0.,10);
  fElectronPt->GetXaxis()->SetTitle("electron Pt");
  fQAHistList->Add(fElectronPt);

  fElectronPtStart = new TH1F("fElectronPtStart","fElectronPtStart",100,0,10);
  fElectronPtStart->GetXaxis()->SetTitle("electron Pt");
  fQAHistList->Add(fElectronPtStart);

  fhGenerator = new TH1F("fhGenerator","fhGenerator",3,-0.5,2.5);
  fhGenerator->GetXaxis()->SetTitle("Which generator");
  fQAHistList->Add(fhGenerator);

  fhOriginKine = new TH1F("fhOriginKine","fhOriginKine",15,-1.5,13.5);
  fhOriginKine->GetXaxis()->SetTitle("Electron source origin Kine");
  fQAHistList->Add(fhOriginKine);

  fhOriginReco = new TH1F("fhOriginReco","fhOriginReco",15,-1.5,13.5);
  fhOriginReco->GetXaxis()->SetTitle("Electron source origin Reco");
  fQAHistList->Add(fhOriginReco);

  fhElGenerator = new TH1F("fhElGenerator","fhElGenerator",13,-0.5,12.5);
  fhElGenerator->GetXaxis()->SetTitle("Which generator + el source");
  fQAHistList->Add(fhElGenerator);
       
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fQAHistList) ;
  cout << "at the end of createing" << endl;
  return;
}


int AliCFSingleTrackEfficiencyTask::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;

  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if(argument.BeginsWith("filterbit=")){
      argument.ReplaceAll("filterbit=", "");
      fSetFilterBit=kTRUE;
      fbit=argument.Atoi();
      AliInfo(Form("Setting filterbit to %d",fbit));
      continue;
    }
    if (argument.BeginsWith("TPCratio=")){
      argument.ReplaceAll("TPCratio=", "");
      fMinRatioTPCcluster=argument.Atof();	    
      AliInfo(Form("tpc clusters /findable %f",fMinRatioTPCcluster));
      continue;
    }	  
    if (argument.BeginsWith("notrequireTOF")){
      AliInfo("No requirement on TOF");
      fRequireTOF=kFALSE;	    
      continue;
    }
    if (argument.BeginsWith("clustersTPCPID=")){
      argument.ReplaceAll("clustersTPCPID=", "");
      fMinNclsTPCPID=argument.Atoi();;
      AliInfo(Form("Requirement on nr clusters for TPC PID %d ",fMinNclsTPCPID));
      continue;
    }
    if(argument.BeginsWith("maxTOFpt=")){
      argument.ReplaceAll("maxTOFpt=", "");
      SetUseTOFPID(kTRUE,argument.Atof());
      AliInfo(Form("Setting max pt for TOF PID to %f",argument.Atof()));
      continue;
    }
    if(argument.BeginsWith("maxradius=")){
      argument.ReplaceAll("maxradius=", "");
      fMaxRadius=argument.Atof();
      fMCCuts->SetMaxRadius(fMaxRadius);
      AliInfo(Form("Setting max radius for particles %f",fMaxRadius));
      continue;
    }
    if(argument.BeginsWith("TOFwhenpresent")){
      SetUseTOFWhenPresent(kTRUE);
      AliInfo("Only use TOF when it's present");
      continue;
    }
    if (argument.BeginsWith("notuseTOFPID")){
      AliInfo("Not Use TOF PID");
      fUseTOFPID=kFALSE;	    
      continue;
    }
    if(argument.BeginsWith("elsource=")){
      argument.ReplaceAll("elsource=", "");
      
      if(argument.CompareTo("HFEPythia")==0){ cout << " setting HFE as source from Pythia " << endl; fUseGenerator=kTRUE; fSelectElSource=kHFPythia;}
      else if(argument.CompareTo("nonHFEHijing")==0) {cout << " setting nonHFE as source from Hijing" << endl; fUseGenerator=kTRUE; fSelectElSource=knonHFHijing;}
      else if(argument.CompareTo("convHijing")==0){ cout << " setting conv as source from Hijing" << endl; fUseGenerator=kTRUE; fSelectElSource=kConvElHijing; }
      
      if(argument.CompareTo("HFE")==0){ cout << " setting HFE as source " << endl; fSelectElSource=kHF;}
      else if(argument.CompareTo("nonHFE")==0) {cout << " setting nonHFE as source " << endl; fSelectElSource=knonHF;}
      else if(argument.CompareTo("conv")==0){ cout << " setting conv as source " << endl; fSelectElSource=kConvEl; }

      else AliFatal(Form("unknown argument '%s'", argument.Data()));
      AliInfo(Form("Selecting only source %d",fSelectElSource));
      continue;
    }
    if(argument.BeginsWith("usegenerator")){
      AliInfo("Select source also based on generator");
      fUseGenerator=kTRUE;	    
      continue;

    }
    if (argument.BeginsWith("notusePhysPrim")){
      AliInfo("Not Use IsPhysicalPrimary()");
      fMCCuts->SetUseIsPhysicalPrimary(kFALSE);
      continue;
    }
    if (argument.BeginsWith("notrejectPileup")){
      AliInfo("Not reject Pileup");
      fMCCuts->SetRejectPileup(kFALSE);
      continue;
    }
    if (argument.BeginsWith("notusePhysicsSelection")){
      AliInfo("Not use physics selection");
      fMCCuts->SetUsePhysicsSelection(kFALSE);
      continue;
    }
    if (argument.BeginsWith("notuseTPCpid")){
      AliInfo("Not Use TPC PID");
      fUseTPCPID=kFALSE;	    
      continue;
    }
    AliWarning(Form("unknown argument '%s'", argument.Data()));
  }
  return 0;
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
    containerInput[0] = mcPart->Pt(); 
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
  //cout << "CheckAODParticles " << endl;
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
  //cout << "nr entries: " << mcArray->GetEntriesFast() << endl;
  /*
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return;
  selectedTracks->SetOwner(kFALSE); // creating new track objects below
  TIter next(mcArray);
  TObject* pObj=NULL;
  Int_t ipart=0;
  while ((pObj=next())) {

    AliAODMCParticle* mcPart=dynamic_cast<AliAODMCParticle*>(pObj);
    //if (!track) continue;*/

  for (Int_t ipart=0; ipart<mcArray->GetEntriesFast(); ipart++) { 
    if (ipart > mcArray->GetEntriesFast()) cout << "should not be here" << endl;
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(ipart));
    //ipart++;
 

    if (!mcPart){
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    //    cout << "mcPartlabel: " <<  TMath::Abs(mcPart->GetLabel()) << endl;
    containerInput[0] = (Float_t)mcPart->Pt();
    containerInput[1] = mcPart->Eta() ;
    containerInput[2] = mcPart->Phi() ;
    containerInput[3] = mcPart->Theta() ;

    // Step 1. Particle passing through Generation criteria and filling
    if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
      AliDebug(3,"MC Particle not passing quality criteria\n");
      continue;
    }
    /*
    cout << "ipart: " << ipart << endl;
    cout << "mcPartlabel: " <<  mcPart->Label() << endl;
    Int_t mother = mcPart->GetMother();
    cout << "mother: " << mother << endl;
    AliAODMCParticle * mcMother= (AliAODMCParticle*)mcArray->At(mother);
    if(!mcMother) {
      cout << "mother not here" << endl;
      return;
    }
    int daug=mcMother->GetDaughter(1);
    int daug2=mcMother->GetDaughter(0);
    cout << "daug: " << daug <<"   daug2: " << daug2 << endl; 
*/
    //cout << "mcPartlabel: " <<  TMath::Abs(mcPart->GetLabel()) << endl;
    /*
      containerInput[5] = originvsGen; */
    int originvsGen= CheckBackgroundSource(mcPart,const_cast<AliVEvent*>(event),kTRUE);
    //    int originvsGen=CheckBackgroundSource(track,const_cast<AliVEvent*>(fEvent));
    fhElGenerator->Fill(originvsGen);
    if(fSelectElSource==kHFPythia && originvsGen!=kHFPythia){
      //fhOriginReco->Fill(fOriginMotherReco);
      continue;
    }
    if(fSelectElSource==knonHFHijing && originvsGen!=knonHFHijing){
      //      fhOriginReco->Fill(fOriginMotherReco);
      continue;
    }

    if(fSelectElSource==kConvElHijing && originvsGen!=kConvElHijing){
      //fhOriginReco->Fill(fOriginMotherReco);
      continue;
    }
    fhOriginKine->Fill(fOriginMotherKine);

    //cout << "originvsGen: " << originvsGen << " kConvElHijing: " << kConvElHijing << endl;
      //    cout << "

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

