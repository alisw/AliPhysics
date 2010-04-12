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

#include "AliCFManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliMCParticle.h"


#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliHFEdisplacedElectrons.h"

#include "AliAnalysisTaskDisplacedElectrons.h"


//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons():
  AliAnalysisTaskSE("Task for displaced electron study")
  , fDebugLevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCuts(0x0)
  , fPID(0x0)
  , fCFM(0x0)
  , fNEvents(0x0)
  , fElectronsPt(0x0)
  , fOutput(0x0)
  , fCorrection(0x0)
  , fDisplacedElectrons(0x0)
  , fHistDisplacedElectrons(0x0)
{
  //
  // Dummy constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(2, TList::Class());  // output
  DefineOutput(1, TList::Class());  // output
  
  SetHasMCData();

  // Initialize pid
  fPID = new AliHFEpid;
  

}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons(const char * name):
  AliAnalysisTaskSE(name)
  , fDebugLevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCuts(0x0)
  , fPID(0x0)
  , fCFM(0x0)
  , fNEvents(0x0)
  , fElectronsPt(0x0)
  , fOutput(0x0)
  , fCorrection(0x0)
  , fDisplacedElectrons(0x0)
  , fHistDisplacedElectrons(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(1, TList::Class());

  SetHasMCData();
}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::AliAnalysisTaskDisplacedElectrons(const AliAnalysisTaskDisplacedElectrons &ref):
  AliAnalysisTaskSE(ref)
  , fDebugLevel(ref.fDebugLevel)
  , fPIDdetectors(ref.fPIDdetectors)
  , fPIDstrategy(ref.fPIDstrategy)
  , fPlugins(ref.fPlugins)
  , fESD(ref.fESD)
  , fMC(ref.fMC)
  , fCuts(ref.fCuts)
  , fPID(ref.fPID)
  , fCFM(ref.fCFM)
  , fNEvents(ref.fNEvents)
  , fElectronsPt(ref.fElectronsPt)
  , fOutput(ref.fOutput)
  , fCorrection(ref.fCorrection)
  , fDisplacedElectrons(ref.fDisplacedElectrons)
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
  fDebugLevel = ref.fDebugLevel;
  fPIDdetectors = ref.fPIDdetectors;
  fPIDstrategy = ref.fPIDstrategy;
  fPlugins = ref.fPlugins;
  fESD = ref.fESD;
  fMC = ref.fMC;
  fPID = ref.fPID;
  fCuts = ref.fCuts;
  fCFM = ref.fCFM;
  fNEvents = ref.fNEvents;
  fOutput = ref.fOutput;
  fCorrection = ref.fCorrection;
  fDisplacedElectrons = ref.fDisplacedElectrons;
  fHistDisplacedElectrons = ref.fHistDisplacedElectrons;

  return *this;
}

//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::~AliAnalysisTaskDisplacedElectrons(){
  //
  // Destructor
  //

  if(fPID) delete fPID;
  
  if(fESD) delete fESD;
  if(fMC) delete fMC;

  if(fCFM) delete fCFM;

  if(fNEvents) delete fNEvents;
  
  if(fElectronsPt) delete fElectronsPt;

  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  

  if(fCorrection){
    fCorrection->Clear();
    delete fCorrection;
  }
 
  if(fDisplacedElectrons) delete fDisplacedElectrons;  

  if(fHistDisplacedElectrons){ 
    fHistDisplacedElectrons->Clear();  
    delete fHistDisplacedElectrons;  
  }
   
  
}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::UserCreateOutputObjects(){
  // create output objects
  // fNEvents
  // MC and Data containers

  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); 
  const Int_t nBins = 14;
  const Float_t ptBins[nBins] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,12.0,16.0,20.0};
  fElectronsPt = new TH1F("esdPt", "p_{T} distribution of identified electrons (HFEpid);p_{T} (GeV/c);dN/dp_{T};", nBins-1, ptBins); 
  if(!fOutput) fOutput = new TList;
  fOutput->SetName("results");
  fOutput->AddAt(fElectronsPt,0);
  fOutput->AddAt(fNEvents,1);
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeEventContainer();
  MakeParticleContainer();
 
  if(!fCorrection) fCorrection = new TList();
  fCorrection->SetName("corrections");
  fCorrection->AddAt(fCFM->GetEventContainer(), 0);
  fCorrection->AddAt(fCFM->GetParticleContainer(), 1);
  fCorrection->Print();

  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetEventContainer()->GetNStep(); istep++)
    fCFM->SetEventCutsList(istep, 0x0);
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();

  }
  
  fCuts->SetCutITSpixel(AliHFEextraCuts::kBoth);
  
  fCuts->Initialize(fCFM);
    
  fPID->SetHasMCData(HasMCData());
  if(!fPIDdetectors.Length() && ! fPIDstrategy) AddPIDdetector("TPC");
  if(fPIDstrategy)
    fPID->InitializePID(Form("Strategy%d", fPIDstrategy));
  else
    fPID->InitializePID(fPIDdetectors.Data());     // Only restrictions to TPC allowed 

  // displaced electron study----------------------------------
  if(GetPlugin(kDisplacedElectrons)){
    
    fDisplacedElectrons = new AliHFEdisplacedElectrons;
    fDisplacedElectrons->SetDebugLevel(fDebugLevel);
    fDisplacedElectrons->SetHasMCData(HasMCData());
    
    if(!fHistDisplacedElectrons) fHistDisplacedElectrons = new TList();
    
    fDisplacedElectrons->CreateOutputs(fHistDisplacedElectrons);
    
    fOutput->AddAt(fHistDisplacedElectrons, 2);
  }


}



//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::UserExec(Option_t *){
  //
  // Run the analysis
  // 

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdH){      
    AliError("No ESD input handler");
    return;
  } else {
    fESD = esdH->GetEvent();
  }
 
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!mcH){       
    AliError("No MC truth handler");
    return;
  } else {
    fMC = mcH->MCEvent();
  }
  
  if(fDebugLevel>=5)AliInfo("Processing ESD events");
  if(!fESD){
    AliError("No ESD Event");
    return;
  }
  
  if(HasMCData()){
    if(fDebugLevel>=5)AliInfo(Form("MC Event: %p", fMC));
    if(!fMC){
      AliError("No MC Event, but MC Data required");
      return;
    }
  }
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }
  
  // process data: ESD tracks with MC information
  Double_t container[8];   // container for the output in THnSparse
  memset(container, 0, sizeof(Double_t) * 8);
  
  Bool_t signal = kTRUE;
  Bool_t alreadyseen = kFALSE;
  LabelContainer cont(fESD->GetNumberOfTracks());

  Int_t nElectrons=0;
  
  AliESDtrack *track = 0x0;    
  
  Double_t nContrib = fESD->GetPrimaryVertex()->GetNContributors();
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;
  if(GetPlugin(kCorrection)){
    fCFM->GetEventContainer()->Fill(&nContrib, AliHFEcuts::kEventStepReconstructed);
  }
  fCFM->SetRecEventInfo(fESD);

  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    track = fESD->GetTrack(itrack);
    
    //
    //  track quality cut: 1st step: ITS and TPC cut; 2nd step: Rec prim; 3rd step: hfe cut on ITS pixel layer
    //
    // require within rapidity range +/-0.9

    //   if((TMath::Abs(track->Eta()))>0.9) continue;

	
    if(GetPlugin(kDisplacedElectrons)) {

      // 1st cut
      // RecKine: ITSTPC cuts : ITS & TPC refit, covmatrix: (2, 2, 0.5, 0.5, 2); min_tpccls: 50, chi2_tpccls: 3.5
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, track)){
	signal = kFALSE;
	continue;
      }
      
      container[0] = track->Pt();
      container[1] = track->Eta();
      container[2] = track->Phi();
      container[3] = track->Charge();
      
      AliStack *stack = fMC->Stack();
      if(!stack) continue;	
      
      AliMCParticle *mctrack = NULL;
      if(HasMCData() ){ 
	mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel())));
	if(!mctrack)  continue;			
	
 	container[4] = mctrack->Pt();
 	container[5] = mctrack->Eta();
 	container[6] = mctrack->Phi();
 	container[7] = mctrack->Charge();

	//	if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
	if(TMath::Abs(mctrack->Eta())>0.9) {
	  signal = kFALSE;
	  continue;
	}  // cut on kinematics

      }  // has MC data

      
      if(signal && GetPlugin(kCorrection)){
	alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));  // double counted track
	cont.Append(TMath::Abs(track->GetLabel()));
	fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecKineITSTPC);
	fCFM->GetParticleContainer()->Fill(&container[0], AliHFEcuts::kStepRecKineITSTPC + 2*AliHFEcuts::kNcutStepsESDtrack);
	if(alreadyseen) 
	  fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsESDtrack);
	
      }  // fill correction --- 1st
    
      // second cut
      // RecPrim: cut on track quality : DCA to vertex max: 3cm and 10cm; reject kink daughters
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, track)) { 
	signal = kFALSE;
	continue;
      }

      if(signal) {
	alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));
	cont.Append(TMath::Abs(track->GetLabel()));
	
	fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecPrim);
	fCFM->GetParticleContainer()->Fill(&container[0], AliHFEcuts::kStepRecPrim + 2*AliHFEcuts::kNcutStepsESDtrack);
	if(alreadyseen) {
	  fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsESDtrack);
	}
      }
      
      //  third cut
      // HFEcuts: ITS layers cuts: ITS pixel layer: kFirst, kSecond, kBoth, kNone or kAny
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)){
	signal = kFALSE;
	continue;
      }

      if(signal) {
	alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));
	cont.Append(TMath::Abs(track->GetLabel()));
	
	fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepHFEcutsITS);
	fCFM->GetParticleContainer()->Fill(&container[0], AliHFEcuts::kStepHFEcutsITS + 2*AliHFEcuts::kNcutStepsESDtrack);
	if(alreadyseen) {
	  fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsESDtrack);
	}
      }
      
      fDisplacedElectrons->FillMCOutput(fESD, track, stack);
      
      // track accepted, do PID --> only electron candidate will be processed
      AliHFEpidObject hfetrack;
      hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
      hfetrack.fRecTrack = track;
      if((fPID->IsSelected(&hfetrack)==kTRUE))
	if(fDebugLevel>=10)
	  AliInfo(Form("ESD info: this particle is %s identified as electron by HFEpid method \n", (fPID->IsSelected(&hfetrack)==kTRUE)?" ":" NOT "));
      if(!fPID->IsSelected(&hfetrack))   continue;
      
         // Fill Containers

      nElectrons++;

      fElectronsPt->Fill(track->Pt());

      if(signal) {
	fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID + 2*AliHFEcuts::kNcutStepsESDtrack);
	fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepPID);
	if(alreadyseen) {
	  fCFM->GetParticleContainer()->Fill(&container[4], (AliHFEcuts::kStepPID + (AliHFEcuts::kNcutStepsESDtrack)));
	}
      }
      
      fDisplacedElectrons->FillESDOutput(fESD, track);
      
    } // analyze displaced electrons plugin switched on
  } // track loop  

  if(fDebugLevel>=5)
    AliInfo(Form("ESD info: number of electrons found in this event: %d\n", nElectrons));


  fNEvents->Fill(1);
  
  PostData(1, fOutput);
  PostData(2, fCorrection);

}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //  
  
  if(GetPlugin(kDisplacedElectrons))
    {
      
      fOutput = dynamic_cast<TList *>(GetOutputData(1));
      fCorrection= dynamic_cast<TList *>(GetOutputData(2));

      if(!fOutput || !fCorrection){
	if(!fCorrection) AliError("correction list not available\n");
	if(!fOutput) AliError("output list not available\n");
	
	return;
      }

      fOutput->Print();
      fCorrection->Print();
      
      AliInfo("analysis done!\n");
      
    }
}






//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::PrintStatus() const {
  
  //
  // Print Analysis status
  //
  printf("\n\tAnalysis Settings\n\t========================================\n");
  printf("\tdisplaced electrons' analysis is %s\n", GetPlugin(kDisplacedElectrons)?"ON":"OFF");
  printf("\tcorrection container  is %s\n", GetPlugin(kCorrection)?"ON":"OFF");
  printf("\tpost processing  is %s\n", GetPlugin(kPostProcess)?"ON":"OFF");
  printf("\tcuts: %s\n", (fCuts != NULL) ? "YES" : "NO");
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
      SETBIT(fPlugins, plug); 
      break;
    case kPostProcess: 
      SETBIT(fPlugins, plug); 
      break;
    case kCorrection:
      SETBIT(fPlugins, plug); 
      break;
    default: 
      AliError("Unknown Plugin");
    };
}


//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::MakeParticleContainer(){
  //
  // Create the particle container (borrowed from AliAnalysisTaskHFE)
  //
  const Int_t kNvar   = 4; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.1, kPtmax = 20.;   // used for fixed pt bins
  //  const Float_t ptBins[14] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,12.0,16.0,20.0};
  const Double_t kEtamin = -0.9, kEtamax = 0.9;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 40; //bins in pt   // used for fixed pt bins
  //  iBin[0] = 13; //bins in pt
  iBin[1] =  8; //bins in eta 
  iBin[2] = 18; // bins in phi
  iBin[3] =  2; // bins in charge
  
  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  //  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=ptBins[i];  // using variable bins
  for(Int_t i=0; i<=iBin[0]; i++) 
     binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);   // fixed pt bin
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;
  for(Int_t i=0; i<=iBin[3]; i++) binEdges[3][i]=1.1*i-1.1; // Numeric precision
  
  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks", 
						 (AliHFEcuts::kNcutStepsTrack + 1 + 2*(AliHFEcuts::kNcutStepsESDtrack + 1)), 
						 kNvar, iBin);

  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    container -> SetBinLimits(ivar, binEdges[ivar]);
  fCFM->SetParticleContainer(container);
  
  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  // add more containers


}

//____________________________________________________________
void AliAnalysisTaskDisplacedElectrons::MakeEventContainer(){
  //
  // Create the event container for the correction framework and link it
  //
  const Int_t kNvar = 1;  // number of variables on the grid: number of tracks per event
  const Double_t kNTrackBound[2] = {-0.5, 200.5};
  const Int_t kNBins = 201;

  AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, &kNBins);

  Double_t trackBins[kNBins];
  for(Int_t ibin = 0; ibin < kNBins; ibin++) trackBins[ibin] = kNTrackBound[0] + static_cast<Double_t>(ibin);
  evCont->SetBinLimits(0,trackBins);

  fCFM->SetEventContainer(evCont);

}



//__________________________________________________________
void AliAnalysisTaskDisplacedElectrons::AddPIDdetector(TString detector){
  //
  // Adding PID detector to the task
  //
  if(!fPIDdetectors.Length()) 
    fPIDdetectors = detector;
  else
    fPIDdetectors += ":" + detector;
}



//____________________________________________________________
AliAnalysisTaskDisplacedElectrons::LabelContainer::LabelContainer(Int_t capacity):
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
Bool_t AliAnalysisTaskDisplacedElectrons::LabelContainer::Append(Int_t label){
  //
  // Add Label to the container
  //
  if(fLast > fEnd) return kFALSE;
  *fLast++ = label;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliAnalysisTaskDisplacedElectrons::LabelContainer::Find(Int_t label) const {
  //
  // Find track in the list of labels
  //
  for(Int_t *entry = fBegin; entry <= fLast; entry++) 
    if(*entry == label) return kTRUE;
  return kFALSE;
}

//____________________________________________________________
Int_t AliAnalysisTaskDisplacedElectrons::LabelContainer::Next()  { 
  //
  // Mimic iterator
  //
  if(fCurrent > fLast) return -1; 
  return *fCurrent++;
}
