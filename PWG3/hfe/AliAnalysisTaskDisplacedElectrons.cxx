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
//
// The analysis task:
// study displaced electrons from beauty and charm 
// with cut on impact parameters in various pT bins
// 
// 
// Authors:
//  Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//  Carlo Bombonati <Carlo.Bombonati@cern.ch>
//

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include "AliLog.h"

#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>
#include <TCanvas.h>

#include "AliAODEvent.h"
#include "AliCFManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPIDResponse.h"
#include "AliStack.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"

#include "AliAnalysisManager.h"

#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliHFEtools.h"
#include "AliHFEdisplacedElectrons.h"
#include "AliAnalysisTaskDisplacedElectrons.h"


//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons():
  AliAnalysisTaskSE("Task for displaced electron study")
  , fDeDebugLevel(0)
  , fNminITSCluster(0)
  , fNminPrimVtxContrib(0)
  , fDePIDdetectors("")
  , fDePIDstrategy(0)
  , fDePlugins(0)
  , fDeCuts(0x0)
  , fDePID(0x0)
  , fDeCFM(0x0)
  , fDisplacedElectrons(0x0)
  , fDeNEvents(0x0)
  , fElectronsMcPt(0x0)
  , fElectronsEsdPt(0x0)
  , fElectronsDataPt(0x0)
  , fDeCorrection(0x0)
  , fDeQA(0x0)
  , fHistDisplacedElectrons(0x0)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons(const char * name):
  AliAnalysisTaskSE(name)
  , fDeDebugLevel(0)
  , fNminITSCluster(0)
  , fNminPrimVtxContrib(0)
  , fDePIDdetectors("")
  , fDePIDstrategy(0)
  , fDePlugins(0)
  , fDeCuts(0x0)
  , fDePID(0x0)
  , fDeCFM(0x0)
  , fDisplacedElectrons(0x0)
  , fDeNEvents(0x0)
  , fElectronsMcPt(0x0)
  , fElectronsEsdPt(0x0)
  , fElectronsDataPt(0x0)
  , fDeCorrection(0x0)
  , fDeQA(0x0)
  , fHistDisplacedElectrons(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());

  // Initialize pid
  fDePID = new AliHFEpid("DEPID");

}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons(const AliAnalysisTaskDisplacedElectrons &ref):
  AliAnalysisTaskSE(ref)
  , fDeDebugLevel(ref.fDeDebugLevel)
  , fNminITSCluster(ref.fNminITSCluster)
  , fNminPrimVtxContrib(ref.fNminPrimVtxContrib)
  , fDePIDdetectors(ref.fDePIDdetectors)
  , fDePIDstrategy(ref.fDePIDstrategy)
  , fDePlugins(ref.fDePlugins)
  , fDeCuts(ref.fDeCuts)
  , fDePID(ref.fDePID)
  , fDeCFM(ref.fDeCFM)
  , fDisplacedElectrons(ref.fDisplacedElectrons)
  , fDeNEvents(ref.fDeNEvents)
  , fElectronsMcPt(ref.fElectronsMcPt)
  , fElectronsEsdPt(ref.fElectronsEsdPt)
  , fElectronsDataPt(ref.fElectronsDataPt)
  , fDeCorrection(ref.fDeCorrection)
  , fDeQA(ref.fDeQA)
  , fHistDisplacedElectrons(ref.fHistDisplacedElectrons)
{
  //
  // Copy Constructor
  //
}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons &AliAnalysisTaskDisplacedElectrons::operator=(const AliAnalysisTaskDisplacedElectrons &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  AliAnalysisTask::operator=(ref);
  fDeDebugLevel = ref.fDeDebugLevel;
  fNminITSCluster = ref.fNminITSCluster;
  fNminPrimVtxContrib = ref.fNminPrimVtxContrib;
  fDePIDdetectors = ref.fDePIDdetectors;
  fDePIDstrategy = ref.fDePIDstrategy;
  fDePlugins = ref.fDePlugins;
  fDePID = ref.fDePID;
  fDeCuts = ref.fDeCuts;
  fDeCFM = ref.fDeCFM;
  fDisplacedElectrons = ref.fDisplacedElectrons;
  fDeNEvents = ref.fDeNEvents;
  fElectronsMcPt = ref.fElectronsMcPt;
  fElectronsEsdPt = ref.fElectronsEsdPt;
  fElectronsDataPt = ref.fElectronsDataPt;
  fDeCorrection = ref.fDeCorrection;
  fDeQA = ref.fDeQA;
  fHistDisplacedElectrons = ref.fHistDisplacedElectrons;

  return *this;
}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::~AliAnalysisTaskDisplacedElectrons(){
  //
  // Destructor
  //

  if(fDePID) delete fDePID;
  if(fDeCFM) delete fDeCFM;
  if(fDisplacedElectrons) delete fDisplacedElectrons;  
  if(fDeNEvents) delete fDeNEvents;
  if(fElectronsMcPt) delete fElectronsMcPt;
  if(fElectronsEsdPt) delete fElectronsEsdPt;
  if(fElectronsDataPt) delete fElectronsDataPt;
  if(fDeCorrection){
    fDeCorrection->Clear();
    delete fDeCorrection;
  }
  if(fDeQA){
    fDeQA->Clear();
    delete fDeQA;
  }
  if(fHistDisplacedElectrons){ 
    fHistDisplacedElectrons->Clear();  
    delete fHistDisplacedElectrons;  
  }
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::UserCreateOutputObjects(){
  // create output objects
  // fDeNEvents
  // MC and Data containers


  if(!fDeQA) fDeQA = new TList;
  fDeQA->SetName("variousQAhistograms");
  
  fDeNEvents = new TH1I("nDeEvents", "Number of Events in the DE Analysis", 2, 0, 2); 
  const Int_t nBins = 14;
  const Float_t ptBins[nBins] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,12.0,16.0,20.0};
  fElectronsMcPt = new TH1F("mcElectronPt", "MC: p_{T} distribution of identified electrons (mcpid);p_{T} (GeV/c);Counts;", nBins-1, ptBins); 
  fElectronsEsdPt = new TH1F("esdElectronPt", "ESD: p_{T} distribution of identified electrons (hfepid);p_{T} (GeV/c);Counts;", nBins-1, ptBins); 
  fElectronsDataPt = new TH1F("dataElectronPt", "DATA: p_{T} distribution of identified electrons (hfepid);p_{T} (GeV/c);Counts;", nBins-1, ptBins); 

  fDeQA->AddAt(fDeNEvents,0);
  if(HasMCData()){
    fDeQA->AddAt(fElectronsMcPt, 1);
    fDeQA->AddAt(fElectronsEsdPt, 2); 
  }
  else{
    fDeQA->AddAt(fElectronsDataPt, 1);  
  }
  // Initialize correction Framework and Cuts
  fDeCFM = new AliCFManager;
  MakeEventContainer();
  MakeParticleContainer();
 
  if(!fDeCorrection) fDeCorrection = new TList();
  fDeCorrection->SetName("deCorrections");
  fDeCorrection->AddAt(fDeCFM->GetEventContainer(), 0);
  fDeCorrection->AddAt(fDeCFM->GetParticleContainer(), 1);
  fDeCorrection->Print();

  for(Int_t istep = 0; istep < fDeCFM->GetEventContainer()->GetNStep(); istep++)
    fDeCFM->SetEventCutsList(istep, 0x0);
  for(Int_t istep = 0; istep < fDeCFM->GetParticleContainer()->GetNStep(); istep++)
    fDeCFM->SetParticleCutsList(istep, 0x0);

  if(!fDeCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fDeCuts = new AliHFEcuts;
    fDeCuts->CreateStandardCuts();    
  }
  
  fDeCuts->Initialize(fDeCFM);
  
  fDePID->SetHasMCData(HasMCData());
  fDePID->AddDetector("TPC", 0);
  fDePID->AddDetector("TOF", 1);
	fDePID->ConfigureTOF();
  fDePID->ConfigureTPCdefaultCut();
  fDePID->InitializePID();     // Only restrictions to TPC allowed   


  // displaced electron study----------------------------------
  if(GetPlugin(kDisplacedElectrons)){
    
    fDisplacedElectrons = new AliHFEdisplacedElectrons;
    fDisplacedElectrons->SetDebugLevel(fDeDebugLevel);
    fDisplacedElectrons->SetHasMCData(HasMCData());
    fDisplacedElectrons->SetMinPrimVtxContrib(fNminPrimVtxContrib);
    fDisplacedElectrons->SetNitsCluster(fNminITSCluster);
  
    if(!fHistDisplacedElectrons) fHistDisplacedElectrons = new TList();
    fDisplacedElectrons->CreateOutputs(fHistDisplacedElectrons);
  }
  
}



//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::UserExec(Option_t *){
  //
  // Run the analysis
  // 

  if(fDeDebugLevel>=10) AliInfo("analyse single event");

  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }
 
  //  
  AliESDInputHandler *inH = dynamic_cast<AliESDInputHandler *>(fInputHandler);
  
  if(!inH){
    AliError("AliESDInputHandler dynamic cast failed!");
    return;
  }
  
  // from now on, only ESD are analyzed
  // using HFE pid, using HFE cuts
  // using CORRFW
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliAODEvent::Class());
  }
  fDePID->SetPIDResponse(pidResponse);

  if(!fDeCuts){
    AliError("HFE cuts not available");
    return;
  }

  // ---- CHOICE OF ANALYSIS ----
  if(HasMCData()){
    // Protect against missing MC trees
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){
      AliError("MC handler cuts not available");
      return;
    }

    if(!mcH->InitOk()) return;
    if(!mcH->TreeK()) return;
    if(!mcH->TreeTR()) return;            
    
    AliDebug(4, Form("MC Event: %p", fMCEvent));
    if(!fMCEvent){
      AliError("No MC Event, but MC Data required");
      return;
    }
    
    ProcessMC(); // PURE MC
    
    if(IsESDanalysis()) ProcessESD(); // ESD WITH MC
    
  }else if(IsESDanalysis()) ProcessData(); // PURE ESD
  

  fDeNEvents->Fill(1);  
  PostData(1, fHistDisplacedElectrons);
  PostData(2, fDeCorrection);
  PostData(3, fDeQA);
  
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::ProcessMC(){
  //
  // handle pure MC analysis
  //

  Int_t nMCelectrons = 0;
  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("ESD event not available");
    return;
  }

  Double_t mcContainer[4];   // container for the output in THnSparse
  memset(mcContainer, 0, sizeof(Double_t) * 4);

  fDeCFM->SetMCEventInfo(fMCEvent);

  Double_t nContributor[1] = {0};
  const AliVVertex *mcPrimVtx = fMCEvent->GetPrimaryVertex();
  if(mcPrimVtx) nContributor[0] = mcPrimVtx->GetNContributors();
  
  // 
  // cut at MC event level
  //

  if(!fDeCFM->CheckEventCuts(AliHFEcuts::kEventStepGenerated, fMCEvent)) return;
  if(GetPlugin(kCorrection)) fDeCFM->GetEventContainer()->Fill(nContributor,AliHFEcuts::kEventStepGenerated);
  
  AliStack *stack = 0x0;
  
  if(!fMCEvent->Stack())return;
  stack = fMCEvent->Stack();
  Int_t nTracks = stack->GetNtrack();

  AliMCParticle *mcTrack = 0x0;

  for(Int_t itrack = 0; itrack<nTracks; itrack++){
    if(!(stack->Particle(itrack))) continue;
    if(mcTrack) mcTrack = 0x0;
    mcTrack = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(itrack));
    if(!mcTrack) continue;
    //TParticle *mcPart = stack->Particle(itrack);
        
    mcContainer[0] = mcTrack->Pt();
    mcContainer[1] = mcTrack->Eta();
    mcContainer[2] = mcTrack->Phi();
    mcContainer[3] = mcTrack->Charge();
    


    if (!stack->IsPhysicalPrimary(mcTrack->GetLabel())) continue;
    // no cut but require primary
    if(GetPlugin(kCorrection))fDeCFM->GetParticleContainer()->Fill(mcContainer, 0);    
        
    // all pions for reference
    if(TMath::Abs(mcTrack->Particle()->GetPdgCode())==AliHFEdisplacedElectrons::kPDGpion && GetPlugin(kCorrection))
      fDeCFM->GetParticleContainer()->Fill(mcContainer, 1);
    
    // cut for signal: all MC electrons
    if(TMath::Abs(mcTrack->Particle()->GetPdgCode())==AliHFEdisplacedElectrons::kPDGelectron && GetPlugin(kCorrection)) 
      fDeCFM->GetParticleContainer()->Fill(mcContainer, 2);

    // cut at track level kinematics: pt and eta
    if(TMath::Abs(mcContainer[1])>=0.8 || mcContainer[0]>20 || mcContainer[0]<0.1) continue;

    if(TMath::Abs(mcTrack->Particle()->GetPdgCode())==AliHFEdisplacedElectrons::kPDGelectron){
      nMCelectrons++;
      fElectronsMcPt->Fill(mcContainer[0]);
    }  

    if(GetPlugin(kCorrection))fDeCFM->GetParticleContainer()->Fill(mcContainer, 3);

    
    // fill MC THnSparse
    fDisplacedElectrons->FillMcOutput(fESD, fMCEvent, mcTrack);
    
  }  // mc track loop 

  if(fDeDebugLevel>=10) printf("there are %d electrons in this MC event", nMCelectrons);
  
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::ProcessESD(){
  
  // this is to handle ESD tracks with MC information
  // MC pid is only used when HFE pid is implemented, for comparison
  // corrections are taken into account
  
  // process data: ESD tracks with MC information
  
  const Int_t kStepPID = AliHFEcuts::kStepHFEcutsTRD + 1;

  Double_t esdContainer[4];   // container for the output in THnSparse
  memset(esdContainer, 0, sizeof(Double_t) * 4);
  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("No ESD event available");
    return;
  }

  fDeCFM->SetRecEventInfo(fESD);
  Double_t nContrib[1] = {fESD->GetPrimaryVertex()->GetNContributors()};

  Bool_t alreadyseen = kFALSE;
  AliLabelContainer cont(fESD->GetNumberOfTracks());
  
  Int_t nHFEelectrons = 0;  
  AliESDtrack *track = 0x0;    
  AliStack *stack = 0x0;

  if(!(stack = fMCEvent->Stack())) return;
  
  //
  // cut at ESD event level
  //
  if(!fDeCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;
  if(GetPlugin(kCorrection)) fDeCFM->GetEventContainer()->Fill(nContrib, AliHFEcuts::kEventStepReconstructed);
  
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    track = fESD->GetTrack(itrack);
    
    if(GetPlugin(kDisplacedElectrons)) {
      
      esdContainer[0] = track->Pt();
      esdContainer[1] = track->Eta();
      esdContainer[2] = track->Phi();
      esdContainer[3] = track->Charge();
      
      // before any cut
      alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));  
      cont.Append(TMath::Abs(track->GetLabel()));  // check double counting
      if(alreadyseen) continue;  // avoid double counting
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+AliHFEcuts::kStepRecNoCut + AliHFEcuts::kNcutStepsMCTrack);	
      
      // 1st track cut
      // RecKine: ITSTPC cuts : ITS & TPC refit, covmatrix: (2, 2, 0.5, 0.5, 2); min_tpccls: 50, chi2_tpccls: 3.5
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack);
      
      // 2nd track cut
      // RecPrim: cut on track quality : DCA to vertex max: 3cm and 10cm; reject kink daughters
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack);
      
      // 3rd track cut
      // HFEcuts: ITS layers cuts: ITS pixel layer: kFirst, kSecond, kBoth, kNone or kAny
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection))fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack);

      /*
      //  4th track cut
      // TRD: number of tracklets in TRD
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack);
      */
      
      // 5th track cut
      // track accepted, do PID 
      // --> only electron candidate will be processed
      AliHFEpidObject hfetrack;
      hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      hfetrack.SetRecTrack(track);
      //if(HasMCData())hfetrack.SetMCTrack(mctrack);
      
      if(!fDePID->IsSelected(&hfetrack)) continue;
      else if(fDeDebugLevel>=10)
	AliInfo("ESD info: this particle is identified as electron by HFEpid method \n");
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(esdContainer, 1+kStepPID + AliHFEcuts::kNcutStepsMCTrack);
    
      // Fill Containers
      nHFEelectrons++;
      fElectronsEsdPt->Fill(esdContainer[0]);
      fDisplacedElectrons->FillEsdOutput(fESD, track, stack);
    
    }  // displaced electron analysis on ESD with MC plugin
  } // track loop  
  
  if(fDeDebugLevel>=10) printf("there are %d HFE electrons in this ESD event", nHFEelectrons);
  
}



//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::ProcessData(){

  // this is a track loop over real data 
  // no MC information at all 
  // HFE pid is used
  
  const Int_t kNcutStepsESDtrack = AliHFEcuts::kNcutStepsRecTrack + 1;
  //const Int_t kNcutStepsTrack = AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack;
  const Int_t kStepPID = AliHFEcuts::kStepHFEcutsTRD + 1;

  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("No ESD event available");
    return;
  }

  Double_t dataContainer[4];   // container for the output in THnSparse
  memset(dataContainer, 0, sizeof(Double_t) * 4);
  
  Bool_t alreadyseen = kFALSE;
  AliLabelContainer cont(fESD->GetNumberOfTracks());


  AliESDtrack *track = 0x0;
  Int_t nHFEelectrons= 0;
  
  fDeCFM->SetRecEventInfo(fESD);
  Double_t nContrib[1] = {fESD->GetPrimaryVertex()->GetNContributors()};
  if(!fDeCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;
  if(GetPlugin(kCorrection)) fDeCFM->GetEventContainer()->Fill(nContrib, AliHFEcuts::kEventStepReconstructed);
  
  
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    track = fESD->GetTrack(itrack);
    
    if(GetPlugin(kDisplacedElectrons)) {
      
      dataContainer[0] = track->Pt();
      dataContainer[1] = track->Eta();
      dataContainer[2] = track->Phi();
      dataContainer[3] = track->Charge();
      
      alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));  // double counted track
      cont.Append(TMath::Abs(track->GetLabel()));
      if(alreadyseen) continue;  // avoid double counting
      if(GetPlugin(kCorrection))fDeCFM->GetParticleContainer()->Fill(&dataContainer[4], 
								     1+AliHFEcuts::kStepRecNoCut + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);

      // 1st track cut
      // RecKine: ITSTPC cuts : ITS & TPC refit, covmatrix: (2, 2, 0.5, 0.5, 2); min_tpccls: 50, chi2_tpccls: 3.5
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(&dataContainer[4], 
								      1+AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);
      
      // 2nd track cut
      // RecPrim: cut on track quality : DCA to vertex max: 3cm and 10cm; reject kink daughters
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, track))  continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(&dataContainer[4], 
								      1+AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);
	
      //  3rd track cut
      // HFEcuts: ITS layers cuts: ITS pixel layer: kFirst, kSecond, kBoth, kNone or kAny
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)) continue;
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(&dataContainer[4],
								      1+AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);
      
      /*
      //  4th track cut
      // TRD: number of tracklets in TRD0
      if(!fDeCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, track)) continue;
      if(GetPlugin(kCorrection)) if(HasMCData())fDeCFM->GetParticleContainer()->Fill(&dataContainer[4], 
										     1+AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);
      */


      // 5th track cut
      // track accepted, do PID --> only electron candidate will be processed

      AliHFEpidObject hfetrack;
      hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      hfetrack.SetRecTrack(track);
      
      if(!fDePID->IsSelected(&hfetrack)) continue;
      else if(fDeDebugLevel>=10)
	AliInfo("ESD info: this particle is identified as electron by HFEpid method \n");
      if(GetPlugin(kCorrection)) fDeCFM->GetParticleContainer()->Fill(dataContainer,  1+kStepPID + AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack);

      nHFEelectrons++;
      fElectronsDataPt->Fill(dataContainer[0]);
      fDisplacedElectrons->FillDataOutput(fESD, track);
    } // analyze displaced electrons plugin switched on
  } // track loop  

  if(fDeDebugLevel>=10) printf("there are %d HFE electrons in this DATA event", nHFEelectrons);
}


//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //  
  
  fHistDisplacedElectrons = dynamic_cast<TList *>(GetOutputData(1));
  fDeCorrection = dynamic_cast<TList *>(GetOutputData(2));
  fDeQA = dynamic_cast<TList *>(GetOutputData(3));
  if(!fDeCorrection) AliError("correction list not available\n");
  if(!fHistDisplacedElectrons) AliError("de list not available\n");
  if(!fDeQA) AliError("qa list is not available\n");
  
  fHistDisplacedElectrons->Print();
  fDeCorrection->Print();
  fDeQA->Print();

  AliInfo("analysis done!\n");
  
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::PrintStatus() const {
  
  //
  // Print Analysis status
  //
  printf("\n");
  printf("\t Analysis Settings\n\t========================================\n");
  printf("\t running over %s\n", HasMCData()?"MC data":"pp collision data");
  printf("\t displaced electrons' analysis is %s\n", GetPlugin(kDisplacedElectrons)?"ON":"OFF");
  printf("\t correction container is %s\n", GetPlugin(kCorrection)?"ON":"OFF");
  printf("\t hfe pid qa is %s\n", GetPlugin(kDePidQA)?"ON":"OFF");
  printf("\t post processing  is %s\n", GetPlugin(kPostProcess)?"ON":"OFF");
  printf("\t cuts: %s\n", (fDeCuts != NULL) ? "YES" : "NO");
  printf("\t ");
  printf("\n");
}

//__________________________________________                                                  
void AliAnalysisTaskDisplacedElectrons::SwitchOnPlugin(Int_t plug){
  //                                            
  // Switch on Plugin          
  // Available:                                  
  //  - analyze impact parameter
  //  - Post Processing                                                                      
  
  switch(plug)
    {
    case kDisplacedElectrons: 
      SETBIT(fDePlugins, plug); 
      break;
    case kCorrection:
      SETBIT(fDePlugins, plug); 
      break;
    case kDePidQA:
      SETBIT(fDePlugins, plug); 
      break;
    case kPostProcess: 
      SETBIT(fDePlugins, plug); 
      break;
    default: 
      AliError("Unknown Plugin");
    };
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNcutStepsESDtrack = AliHFEcuts::kNcutStepsRecTrack + 1;
  const Int_t kNcutStepsTrack = AliHFEcuts::kNcutStepsMCTrack + kNcutStepsESDtrack;
 
  const Int_t kNvar   = 4;
  //number of variables on the grid:pt,eta, phi, charge
  const Double_t kPtbound[2] = {0.1, 10.};
  const Double_t kEtabound[2] = {-0.8, 0.8};
  const Double_t kPhibound[2] = {0., 2. * TMath::Pi()};

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 40; // bins in pt
  iBin[1] =  8; // bins in eta 
  iBin[2] = 18; // bins in phi
  iBin[3] =  2; // bins in charge

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  binEdges[0] = AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);
  binEdges[1] = AliHFEtools::MakeLinearBinning(iBin[1], kEtabound[0], kEtabound[1]);
  binEdges[2] = AliHFEtools::MakeLinearBinning(iBin[2], kPhibound[0], kPhibound[1]);
  binEdges[3] = AliHFEtools::MakeLinearBinning(iBin[3], -1.1, 1.1); // Numeric precision

  //------------------------------------------------
  //     one "container" for MC+ESD+Data          
  //----------pure MC track-------------------------
  // 0: MC generated
  // 1: MC pion total  ---- be careful!!!!
  // 2: MC electrons total
  // 3: MC electrons in acceptance 
  //-------ESD track with MC info-------------------
  // 4: ESD track with MC: no cut
  // 5: ESD track with MC: cut on kine its tpc
  // 6: ESD track with MC: rec prim
  // 7: ESD track with MC: hfe cuts its
  // 8: ESD track with MC: hfe cuts trd
  // 9: ESD track with MC: hfe pid 
  //-----------data track---------------------------
  // 10: DATA track wo MC: no cut
  // 11: DATA track wo MC: cut on kine its tpc
  // 12: DATA track wo MC: rec prim
  // 13: DATA track wo MC: hfe cuts its
  // 14: DATA track wo MC: hfe cuts trd
  // 15: DATA track wo MC: hfe pid 
  //------------------------------------------------

  AliCFContainer* container = new AliCFContainer("deTrackContainer", "Container for tracks", 
						 (1 + kNcutStepsTrack + kNcutStepsESDtrack), kNvar, iBin);
  
  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++){
    container -> SetBinLimits(ivar, binEdges[ivar]);
  }
  fDeCFM->SetParticleContainer(container);
}


//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::MakeEventContainer(){
  //
  // Create the event container for the correction framework and link it
  //
  
  // event container
  // 0: MC event
  // 1: ESD event

  const Int_t kNvar = 1;  // number of variables on the grid: number of tracks per event
  const Double_t kNTrackBound[2] = {-0.5, 200.5};
  const Int_t kNBins = 201;

  AliCFContainer *evCont = new AliCFContainer("deEventContainer", "Container for DE events", AliHFEcuts::kNcutStepsEvent, kNvar, &kNBins);

  Double_t trackBins[kNBins];
  for(Int_t ibin = 0; ibin < kNBins; ibin++) trackBins[ibin] = kNTrackBound[0] + static_cast<Double_t>(ibin);
  evCont->SetBinLimits(0,trackBins);

  fDeCFM->SetEventContainer(evCont);

}



//__________________________________________________________
void AliAnalysisTaskDisplacedElectrons::AddPIDdetector(TString detector){
  //
  // Adding PID detector to the task
  //
  if(!fDePIDdetectors.Length()) 
    fDePIDdetectors = detector;
  else
    fDePIDdetectors += ":" + detector;
}



//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliLabelContainer::AliLabelContainer(Int_t capacity):
  fContainer(NULL),
  fBegin(NULL),
  fEnd(NULL),
  fLast(NULL),
  fCurrent(NULL)
{
  //
  // Default constructor
  //
  fContainer = new Int_t[capacity];
  fBegin = &fContainer[0];
  fEnd = &fContainer[capacity - 1];
  fLast = fCurrent = fBegin;
}

//____________________________________________________________
Bool_t AliAnalysisTaskDisplacedElectrons::AliLabelContainer::Append(Int_t label){
  //
  // Add Label to the container
  //
  if(fLast > fEnd) return kFALSE;
  *fLast++ = label;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliAnalysisTaskDisplacedElectrons::AliLabelContainer::Find(Int_t label) const {
  //
  // Find track in the list of labels
  //
  for(Int_t *entry = fBegin; entry <= fLast; entry++) 
    if(*entry == label) return kTRUE;
  return kFALSE;
}

//____________________________________________________________
Int_t AliAnalysisTaskDisplacedElectrons::AliLabelContainer::Next() { 
  //
  // Mimic iterator
  //
  if(fCurrent > fLast) return -1; 
  fCurrent++;
  return *fCurrent;
}
