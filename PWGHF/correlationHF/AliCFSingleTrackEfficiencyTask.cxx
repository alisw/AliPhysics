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
  fReducedMode(kFALSE),
  fUsePt(kTRUE)

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
  fReducedMode(kFALSE),
  fUsePt(kTRUE)
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
    fReducedMode=c.fReducedMode;
    fUsePt=c.fUsePt;

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
  fReducedMode(c.fReducedMode),
  fUsePt(c.fUsePt)

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
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
  if (fTrackCuts)           delete fTrackCuts;
  if (fMCCuts)              delete fMCCuts;


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
  ((TH1F*)fQAHistList->FindObject("fHistEventsProcessed"))->Fill(0); // # of Event proceed        
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
    //    PostData(1,fHistEventsProcessed) ;
    PostData(2,fCFManager->GetParticleContainer()) ;
    PostData(3,fQAHistList) ;
    return;
  }
  
  ((TH1F*)fQAHistList->FindObject("fHistEventsProcessed"))->Fill(1); // # of Event after passing MC cuts
       
  //Filling of MC generated particle -> below function 
  if(!isAOD) CheckESDParticles();
  else CheckAODParticles();
       
  // Step 0B. MC Reconstruction Event Selection for ESDs/ AODs
  Bool_t isRecoEventOk = fMCCuts->IsRecoEventSelected(fEvent);
     
  int nr_input=5;
  if(fReducedMode)
    nr_input=4;
  Double_t containerInput[nr_input] ;
  Double_t containerInputMC[nr_input] ; // for true pt
     
  if(isRecoEventOk) {
	 
    ((TH1F*)fQAHistList->FindObject("fHistEventsProcessed"))->Fill(2); // # of Event after passing all cuts
    const AliVVertex *vertex = fEvent->GetPrimaryVertex();
    containerInput[kZvt] = vertex->GetZ(); // Z Vertex of Event 
    containerInputMC[kZvt]  =  containerInput[kZvt];
    //cout << "Z vtx of current event is @ Event " <<  containerInputMC[4] << endl; 
	 

    // Reco tracks track loop
    for (Int_t iTrack = 0; iTrack<fEvent->GetNumberOfTracks(); iTrack++) {
      fOriginMotherReco=-1;
      track = fEvent->GetTrack(iTrack);
      ((TH1F*)fQAHistList->FindObject("fElectronPtStart"))->Fill(track->Pt());

      AliAODTrack *aodtrack = NULL;
      if(isAOD) aodtrack= static_cast<AliAODTrack *>(track);
	   
      Double_t mom[3];
      track->PxPyPz(mom);
      Double_t pt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      if(fUsePt)
	containerInput[kPt] = pt ;
      else
	containerInput[kPt]=track->P();
      containerInput[kEta] = track->Eta();
      containerInput[kPhi] = track->Phi() ;
      if(!fReducedMode) containerInput[kTheta] = track->Theta() ;

      /*
	int originvsGen=CheckBackgroundSource(track,const_cast<AliVEvent*>(fEvent));

	containerInput[5]=originvsGen;
	fhElGenerator->Fill(originvsGen);*/

      // Step 4. Track that are recostructed and filling
      if(!fReducedMode) fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;
	   
 
      // Step 5. Track that are recostructed and +Kine acceptance filling
      // Pt and eta
      if (!fMCCuts->IsRecoParticleKineAcceptance(track)) continue;
      if(!fReducedMode) fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoKineCuts) ;
      
      // is track associated to particle ? if yes + implimenting the physical primary..
      Int_t label = track->GetLabel();
      if (label<0) {
	AliDebug(3,"Particle not matching MC label \n");
	continue;
      }

      // Step 6. Track that are recostructed + true pt and filling
      AliVParticle *mcPart  = (AliVParticle*)fMCEvent->GetTrack(label);
      if(fUsePt)
	containerInputMC[kPt] = mcPart->Pt();
      else
	containerInputMC[kPt] = mcPart->P();
      containerInputMC[kEta] = mcPart->Eta() ;
      containerInputMC[kPhi] = mcPart->Phi() ;
      if(!fReducedMode)containerInputMC[kTheta] = mcPart->Theta() ;
      if (!fMCCuts->IsMCParticleGenerated(mcPart)) continue;
      AliAODMCParticle *mcPart2=dynamic_cast<AliAODMCParticle*>(mcPart);
      
      if((aodtrack->Charge() != -1) && (aodtrack->Charge() != 1)) cout << "charge: " << aodtrack->Charge() << endl;

      Bool_t selected=kTRUE;


      int originvsGen=CheckBackgroundSource(track,const_cast<AliVEvent*>(fEvent));

      ((TH1F*)fQAHistList->FindObject("fhElGenerator"))->Fill(originvsGen);

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
	if(fSelectElSource==kHadronHijing && originvsGen!=kHadronHijing){
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
	if(fSelectElSource==kHadron && originvsGen!=kHadron){
	  selected=kFALSE;
	}
      }

      if(!selected) continue;
      
      Double_t x=mcPart2->Xv();
      Double_t y=mcPart2->Yv();
      double radius=TMath::Sqrt(x*x+y*y);
      if(radius > fMaxRadius) { continue;}
      
      
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
	AliDebug(3,"Reconstructed track not passing first quality criteria\n");
	continue;
      }

      AliDebug(2,"Reconstructed track pass first quality criteria\n");
      //fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepReconstructedFirstTrackCutsMC);
      if(fReducedMode)
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedRecoFirstQualityCuts);
      else
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoFirstQualityCuts);
      

      Bool_t useTOFPID=kTRUE;
      if(track->Pt() > fMaxPtForTOFPID) useTOFPID=kFALSE;
      if(useTOFPID) 
	AliDebug(2,Form("Pt: %f, use CombinedPID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtForTOFPID));
      else 
	AliDebug(2,Form("Pt: %f, use only TPC PID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtForTOFPID));
 
      // DxHFE: Add here the extra cuts
      // 1. TPC PID clusters
      if(fMinNclsTPCPID>0 ){
	Int_t nclsTPCPID = tmptrack->GetTPCsignalN();
	if(nclsTPCPID<fMinNclsTPCPID){
	  //cout << "cut due to nr cls TPC PID " << endl;
	  AliDebug(2,Form("nlcTPCPID NOT selected - nrclusters=%d", nclsTPCPID));
	  continue;
	}
      }
 
      // 2. ratio TPC/findable
      if(fMinRatioTPCcluster>0){	     
	if(isAOD){
	  //AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
	  const TBits &clusterTPC = aodtrack->GetTPCClusterMap();
	  Int_t allclusters = clusterTPC.CountBits();
	  Double_t clusterRatio = vtrack->GetTPCNclsF() ? static_cast<Double_t>(allclusters)/static_cast<Double_t>(vtrack->GetTPCNclsF()) : 1.; // foundall/findable
	  AliDebug(2,Form("clusterRatio: %f  of all clusters: %d",clusterRatio,allclusters));
	  if(clusterRatio <= fMinRatioTPCcluster){
	    AliDebug(2,"clusterRatio NOT selected");
	    //cout << "cut due to clusterratio" << endl;
	    continue;
	  }
	}
      }
      
      // 3. TOF matching
      if(fRequireTOF && useTOFPID){
	if(!(vtrack->GetStatus() & AliESDtrack::kTOFpid)){
	  if(fUseTOFonlyWhenPresent) 
	    useTOFPID=kFALSE;
	  else{
	    AliDebug(2,"Cut due to TOF requirement");
	    continue;
	  }
	}
      }
      ((TH1F*)fQAHistList->FindObject("fElectronPt"))->Fill(tmptrack->Pt());
      
      AliDebug(2,"Reconstructed track pass quality criteria\n");
      //fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepReconstructedMC);
      if(fReducedMode)
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedRecoQualityCuts);
      else
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoQualityCuts);
	   
      // PID requirement
      ((TH2F*)fQAHistList->FindObject("fhdEdxSigmavsEta"))->Fill(vtrack->Eta(), pidResponse->NumberOfSigmasTPC(vtrack, AliPID::kElectron));
      ((TH2F*)fQAHistList->FindObject("fhdEdxvsEta"))->Fill(vtrack->Eta(), vtrack->GetTPCsignal());

      //also check for pdg first????
      if(fUseTPCPID){
	Float_t tpcNsigma = pidResponse->NumberOfSigmasTPC(vtrack, AliPID::kElectron); // change to particle
	AliDebug(2, Form("Number of sigmas in TPC: %f", tpcNsigma));
	//cout << "sigmaTPC " << tpcNsigma << "   " << fTPCnSigmaMin << " - " << fTPCnSigmaMax << endl;
	if(tpcNsigma<fTPCnSigmaMin || tpcNsigma>fTPCnSigmaMax) { continue;}

      }

      ((TH2F*)fQAHistList->FindObject("fhdEdxSigmavsEtaTPC"))->Fill(vtrack->Eta(), pidResponse->NumberOfSigmasTPC(vtrack, AliPID::kElectron));
      ((TH2F*)fQAHistList->FindObject("fhdEdxvsEtaTPC"))->Fill(vtrack->Eta(), vtrack->GetTPCsignal());
    
      if(useTOFPID && fUseTOFPID){

	// if fMaxPtCombinedPID is set to lower than upper Ptlimit (10GeV/c), will separate
	// PID into two regions: below fMaxptCombinedPID - both TPC and TOF, above only TPC
	Float_t tofNsigma = pidResponse->NumberOfSigmasTOF(vtrack, AliPID::kElectron); //change to particle
	AliDebug(2, Form("Number of sigmas in TOF: %f", tofNsigma));
	// for now: Assume symmetric cut for TOF
	if(TMath::Abs(tofNsigma) > fTOFnSigma) {continue;}
      }


      ((TH2F*)fQAHistList->FindObject("fhdEdxSigmavsEtaTPCTOF"))->Fill(vtrack->Eta(), pidResponse->NumberOfSigmasTPC(vtrack, AliPID::kElectron));
      ((TH2F*)fQAHistList->FindObject("fhdEdxvsEtaTPCTOF"))->Fill(vtrack->Eta(), vtrack->GetTPCsignal());

      // fill container for tracks
      if(fReducedMode)
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedRecoPID);
      else{
	fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepRecoPIDMC);
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPID);
      }

      // invariant mass method - Not sure if needed here...
      fSelNHFE->FindNonHFE(iTrack, aodtrack, const_cast<AliVEvent*>(fEvent));
      if(fSelNHFE->IsLS() || fSelNHFE->IsULS())
	{
	  //Not selected
	  AliDebug(2,"Cut: Invmass");
	  continue;
	}

      // fill container for tracks, with all electrons
      //fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepRecoInvMassMC);
      if(fReducedMode)
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedRecoInvMass);
      else
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoInvMass);

      ((TH1F*)fQAHistList->FindObject("fhOriginReco"))->Fill(fOriginMotherReco);
	
      if(fSelectElSource==kHadronHijing || fSelectElSource==kHadron){((TH1F*)fQAHistList->FindObject("fPDGHadron"))->Fill(TMath::Abs(mcPart->PdgCode()));}
    if(isAOD) delete tmptrack;

	
    }
      
	
  }
  else AliDebug(3,"Event not passing quality criteria\n");
       
  
  // PostData(0) is taken care of by AliAnalysisTaskSE 
  //  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fQAHistList) ;
       
  return;

}
int AliCFSingleTrackEfficiencyTask::CheckBackgroundSource(AliVParticle* track, const AliVEvent* pEvent,Bool_t useMCarray){

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
  if(!(fSelectElSource==kHadronHijing || fSelectElSource==kHadron)){
    if(useMCarray) fElectronsKine->CheckMC(track, (AliVEvent*)pEvent);
    else fElectrons->CheckMC(track, (AliVEvent*)pEvent);

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
  }
  else{
    originvsGen=kHadron;
    if(fUseGenerator){
      if(generator==0)
	originvsGen=kHadronGen0;
      if(generator==1)
	originvsGen=kHadronHijing;
      if(generator==2)
	originvsGen=kHadronPythia;
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
  
  fQAHistList->Add(CreateControlHistogram("fInvMassLS","Invariant mass LS",1000,0,0.5));
  fQAHistList->Add(CreateControlHistogram("fInvMassULS","Invariant mass ULS",1000,0,0.5));

  fSelNHFE->SetHistMass((TH1F*)fQAHistList->FindObject("fInvMassULS"));
  fSelNHFE->SetHistMassBack((TH1F*)fQAHistList->FindObject("fInvMassLS"));

  // Setting up the electron selection class for MC
  TString option="usekine"; 
  if(!fMCCuts->GetUseIsPhysicalPrimary()) option+="-keepsecondary";
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

  const char* cutBinNames[]={
    "All events",
    "Good MC events",
    "Good Reconstructed events"
  };
  
  fQAHistList->Add(CreateControlHistogram("fHistEventsProcessed","Event Info",3,cutBinNames));

  fQAHistList->Add(CreateControlHistogram("fElectronPt","electron pt",100,0.,10));
  fQAHistList->Add(CreateControlHistogram("fElectronPtStart","electron pt Start",100,0.,10));
  fQAHistList->Add(CreateControlHistogram("fhGenerator","Which Generator",3,-0.5,2.5));
  fQAHistList->Add(CreateControlHistogram("fhOriginKine","Electron source origin Kine",15,-1.5,13.5));
  fQAHistList->Add(CreateControlHistogram("fhOriginReco","Electron source origin Reco",15,-1.5,13.5));

  fQAHistList->Add(CreateControlHistogram("fhElGenerator","Which generator + el source",kNrSources,-0.5,kNrSources-0.5));

  double dEdxvseta[6]={100,-1.,1.,200,0., 200.};
  double sigmavseta[6]={100,-1.,1.,200,-10., 10.};
  fQAHistList->Add(CreateControl2DHistogram("fhdEdxvsEta", "dEdx vs eta",dEdxvseta ,"dE/dx","#eta"));
  fQAHistList->Add(CreateControl2DHistogram("fhdEdxSigmavsEta", "dEdx vs eta", sigmavseta,"dE/dx","#eta"));

  fQAHistList->Add(CreateControl2DHistogram("fhdEdxvsEtaTPC", "dEdx vs eta",dEdxvseta ,"dE/dx","#eta"));
  fQAHistList->Add(CreateControl2DHistogram("fhdEdxSigmavsEtaTPC", "dEdx vs eta", sigmavseta,"dE/dx","#eta"));

  fQAHistList->Add(CreateControl2DHistogram("fhdEdxvsEtaTPCTOF", "dEdx vs eta",dEdxvseta ,"dE/dx","#eta"));
  fQAHistList->Add(CreateControl2DHistogram("fhdEdxSigmavsEtaTPCTOF", "dEdx vs eta", sigmavseta,"dE/dx","#eta"));
  fQAHistList->Add(CreateControlHistogram("fPDGHadron","PDG of hadrons",5000));  
       
  //PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fQAHistList) ;

  return;
}


TH1* AliCFSingleTrackEfficiencyTask::CreateControlHistogram(const char* name,
							    const char* title,
							    int nBins,
							    double min,
							    double max,
							    const char** binLabels) const
{
  /// create control histogram
  std::auto_ptr<TH1> h(new TH1D(name, title, nBins, min, max));
  if (!h.get()) return NULL;
  if (binLabels) {
    for (int iLabel=0; iLabel<nBins; iLabel++) {
      h->GetXaxis()->SetBinLabel(iLabel+1, binLabels[iLabel]);    
    }
  }
  
  return h.release();
}


TH2* AliCFSingleTrackEfficiencyTask::CreateControl2DHistogram(const char* name,
							      const char* title,
							      double* nBins,
							      const char* xaxis,
							      const char* yaxis
							      ) const
{
  /// create control 2D histogram. Requires as input:
  // name = name of histogram 
  // title = title of histogram
  // nBins (array with 6 elements) containing apropriate binning and range for x and y axis
  // xaxis = title of x axis 
  // yaxis = title of y axis 

  std::auto_ptr<TH2> h(new TH2D(name, title, (Int_t)nBins[0], nBins[1], nBins[2], (Int_t)nBins[3], nBins[4],nBins[5]));
  if (!h.get()) return NULL;
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  
  return h.release();
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
    if (argument.BeginsWith("invmasscut=")){
      argument.ReplaceAll("invmasscut=", "");
      fInvMassLow=argument.Atof();	    
      AliInfo(Form("Invariant mass cut %f",fInvMassLow));
      continue;
    }	   
    if (argument.BeginsWith("notrequireTOF")){
      AliInfo("No requirement on TOF");
      fRequireTOF=kFALSE;	    
      continue;
    }
    if (argument.BeginsWith("notuseTOFPID")){
      AliInfo("Not Use TOF PID");
      fUseTOFPID=kFALSE;
      fRequireTOF=kFALSE;	    
      continue;
    }
    if(argument.BeginsWith("TOFwhenpresent")){
      SetUseTOFWhenPresent(kTRUE);
      AliInfo("Only use TOF when it's present");
      continue;
    }

    if (argument.BeginsWith("useglobalmomentum") || argument.BeginsWith("useglobalmom")){
      AliInfo("Use global momentum");
      fUsePt=kFALSE;	    
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

    if (argument.BeginsWith("reducedmode")){
      AliInfo("Running in reduced mode");
      fReducedMode=kTRUE;	    
      continue;
    }
    if(argument.BeginsWith("elsource=")){
      argument.ReplaceAll("elsource=", "");
      
      if(argument.CompareTo("HFEPythia")==0){ cout << " setting HFE as source from Pythia " << endl; fUseGenerator=kTRUE; fSelectElSource=kHFPythia;}
      else if(argument.CompareTo("nonHFEHijing")==0) {cout << " setting nonHFE as source from Hijing" << endl; fUseGenerator=kTRUE; fSelectElSource=knonHFHijing;}
      else if(argument.CompareTo("convHijing")==0){ cout << " setting conv as source from Hijing" << endl; fUseGenerator=kTRUE; fSelectElSource=kConvElHijing; }
      else if(argument.CompareTo("hadronHijing")==0){ cout << " setting hadron as source from Hijing" << endl; fMCCuts->SetSelectPdg(AliSingleTrackEffCuts::kPDGSelectNotPdg); fUseGenerator=kTRUE; fSelectElSource=kHadronHijing; }
      else if(argument.CompareTo("HFE")==0){ cout << " setting HFE as source " << endl; fSelectElSource=kHF;}
      else if(argument.CompareTo("nonHFE")==0) {cout << " setting nonHFE as source " << endl; fSelectElSource=knonHF;}
      else if(argument.CompareTo("conv")==0){ cout << " setting conv as source " << endl; fSelectElSource=kConvEl; }
      else if(argument.CompareTo("hadron")==0){ cout << " setting hadron as source " << endl; fMCCuts->SetSelectPdg(AliSingleTrackEffCuts::kPDGSelectNotPdg); fSelectElSource=kHadron; }

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
  
  int nr_input=5;
  if(fReducedMode)
    nr_input=4;
  Double_t containerInput[nr_input] ; //number of variables
       
  TArrayF vtxPos(3); // for Z vtx
  AliGenEventHeader *genHeader;
  genHeader = fMCEvent->GenEventHeader();
  genHeader->PrimaryVertex(vtxPos);
  containerInput[kZvt]  = vtxPos[2];
       
  //loop on the MC Gen Partiles
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
	 
    AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);  
    if(fUsePt)    containerInput[kPt] = mcPart->Pt(); 
    else
      containerInput[kPt] = mcPart->P(); 
    containerInput[kEta] = mcPart->Eta() ;
    containerInput[kPhi] = mcPart->Phi() ;
    if(!fReducedMode) containerInput[kTheta] = mcPart->Theta() ;
	 
	 
    // Step 1. Particle passing through Generation criteria and filling
    if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
      AliDebug(3,"MC Particle not passing through genetations criteria\n");
      continue;
    }
    if(!fReducedMode)
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);
	 
     
	 
    // Step 2. Particle passing through Kinematic criteria and filling
    if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
      AliDebug(3,"MC Particle not in the kine acceptance\n");
      continue;
    }
    if(fReducedMode)
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedMCKineCut);
    else
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);
	 
    // Step 3. Particle passing through Track ref criteria and filling
    // did leave signal (enough clusters) on the detector
    if( !fMCCuts->IsMCParticleInReconstructable(mcPart) ) {
      AliDebug(3,"MC Particle not in the reconstructible\n");
      continue;
    }
    if(!fReducedMode)
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
       
  int nr_input=5;
  if(fReducedMode)
    nr_input=4;    
  Double_t containerInput[nr_input] ;
  containerInput[kZvt]  = mcHeader->GetVtxZ();

  for (Int_t ipart=0; ipart<mcArray->GetEntriesFast(); ipart++) { 
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(ipart));

    if (!mcPart){
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    //    cout << "mcPartlabel: " <<  TMath::Abs(mcPart->GetLabel()) << endl;
    if(fUsePt)
      containerInput[kPt] = (Float_t)mcPart->Pt();
    else
      containerInput[kPt] = (Float_t)mcPart->P();
    containerInput[kEta] = mcPart->Eta() ;
    containerInput[kPhi] = mcPart->Phi() ;
    if(!fReducedMode) containerInput[kTheta] = mcPart->Theta() ;

    // Step 1. Particle passing through Generation criteria and filling
    if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
      AliDebug(3,"MC Particle not passing quality criteria\n");
      continue;
    }
 
    int originvsGen= CheckBackgroundSource(mcPart,const_cast<AliVEvent*>(event),kTRUE);
    ((TH1F*)fQAHistList->FindObject("fhElGenerator"))->Fill(originvsGen);

    if(fUseGenerator){
      if(fSelectElSource==kHFPythia && originvsGen!=kHFPythia){
	continue;
      }
      if(fSelectElSource==knonHFHijing && originvsGen!=knonHFHijing){
	continue;
      }

      if(fSelectElSource==kConvElHijing && originvsGen!=kConvElHijing){
	continue;
      }
      if(fSelectElSource==kHadronHijing && originvsGen!=kHadronHijing){
	continue;
      }
    }
    else{

      if(fSelectElSource==kHF && originvsGen!=kHF){
	continue;
      }
      if(fSelectElSource==knonHF && originvsGen!=knonHF){
	continue;
      }

      if(fSelectElSource==kConvEl && originvsGen!=kConvEl){
	continue;
      }
      if(fSelectElSource==kHadron && originvsGen!=kHadron){
	continue;
      }
    }

    if(!fReducedMode) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);
	 

    // Step 2. Particle passing through Kinematic criteria and filling
    if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
      AliDebug(3,"MC Particle not in the acceptance\n");
      continue;
    }
    ((TH1F*)fQAHistList->FindObject("fhOriginKine"))->Fill(fOriginMotherKine);

    if(fReducedMode)
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepRedMCKineCut);
    else
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);
	 

    // Step 3. Particle passing through Track Ref criteria and filling
    // but no info available for Track ref in AOD fillng same as above
    if(!fReducedMode) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCAccpCut);
     

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

