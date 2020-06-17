/*
  Comments:
  * 19 aug 2015: trichert: changed mc primary method
  * 19 aug 2015: trichert: changed vtx cut
  * 08 sep 2015: trichert: removed rapidity v0data->y = aodV0->Y(pdgV0); from AOD MC (made it crash), new
  * 15 sep 2015: trichert: fixed bug: p_track->n_track
  * 15 sep 2015: trichert: flags on injected MC particles: 
  - MC track: pidCodeMC>10 = injected
  - track: pidCode>10 = injected
  - p_track and n_track: n/p_pidCode>10 = injected
  * 23 sep 2015: trichert: hardcoded trigger conditions, search for TRIGGER CONDITION
  * 21 oct 2015: trichert: uncommented hardcoded trigger conditions (search for TRIGGER CONDITION)
  * 14 dec 2015: trichert: changed aodTrack->TestFilterBit(1) to aodTrack->TestFilterBit(128) for TPC_only, and pTrack->TestFilterBit(128) and nTrack->TestFilterBit(128) (from 1)
  * 11 dec 2017: implemented a setting function for cosPA ("fCosPACut")
  * 13 dec 2017: cleanup + filtering to reduce output size: 
  - eta: 0.9->0.8
  - setter function for DecayRCut: if (lV0Radius < fDecayRCut ) continue;
  * 15 dec 2017: setter for trigger, removed trigger settings in addtask parameters (now set it in train config)
  * 18 dec 2017: several trigger settings
  * 04 jan 2018: changed back to old trigger method
  * 17 may 2018: last debugging&cleaning -> zvtx is now ok in local tests 
  * 26 jun 2018: change logic to use AODEvent() if InputEvent() didn't work
  * 05 jul 2018: change type to AOD
  * 10 jul 2018: clean-up
  * 10 jul 2018: removed analysisMode == kTPCTrk and TPCBranch
  * 10 jul 2018: found and corrected bug:  if(pPid){nncl = nPid->GetTPCsignalN();}->if(nPid){nncl = nPid->GetTPCsignalN();}
  * 13 aug 2018: fill tree with AliAODEvent::GetNumberOfTracks() so I can plot is vs centrality
                 and: switch to use AliMultSelection instead of AliCentrality: CentFramework(CentFrameworkAliCen) kTRUE: use AliCentrality, kFALSE: use AliMultSelection
  * 31 oct 2018: fill tree with fEvent->refMult
  * 17 jan 2019: add x-check histos to get quick overview of which: finder, calib, trig
                 fTriggerInt;         // 0 = kMB, 1 = kCent, 2 = kSemiCent
                 fV0Finder;           // 0 = oldFinder, 1 = newFinder
                 fCentFramework;      // 0 = AliCentrality, 1 = AliMultSelection
  * 23 jan 2019: add more filter bit flags for v0 daughters, and make sure TPCrefit, nCrossedRowsTPC, and findable cuts are always there 
  * 30 jan 2019: setting functions for: SetContributorsVtxCut, SetContributorsVtxSPDCut, SetZvsSPDvtxCorrCut, SetVtxR2Cut, SetCrossedRowsCut, SetCrossedOverFindableCut, and removing the eta cut limit of 0.8 on V0s (but keeping it on tracks & daughter tracks)
  * 06 feb 2019: implement reject kink daughters option
  * 17 jan 2020: removing esd part of the code since heavily outdated 
  * 29 jan 2020: oob pileup cut and AliEventCuts for run-2 analyses
  * 23 apr 2020: preparing for pp, renaming PileupCut, clean-up
  * 4 jun  2020: cleaning up track filter bit flags, storing n sigma dedx, plus some other small things 
  * 16 jun 2020: fixing INEL>0 flag in vtxstatus

  */

// Responsible:
// Tuva Richert (Lund)


#include "AliAnalysisTaskHighPtDeDx.h"
#include <bitset>

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h> 
#include <AliAODVZERO.h>
#include "AliPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "AliMultSelection.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>
//#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include "AliCentrality.h" 
#include <AliKFVertex.h>
#include <AliAODVertex.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 
#include <iostream>
#include "AliEventCuts.h"

// STL includes
#include <iostream>
using namespace std;



ClassImp(AliAnalysisTaskHighPtDeDx)

const Double_t AliAnalysisTaskHighPtDeDx::fgkClight = 2.99792458e-2;

//_____________________________________________________________________________
AliAnalysisTaskHighPtDeDx::AliAnalysisTaskHighPtDeDx():
  AliAnalysisTaskSE(),
  fEventCuts(0),
  fAOD(0x0),
  fPIDResponse(0x0),
  fMC(0x0),
  fMCArray(0x0),
  fTrackFilter(0x0),
  fTrackFilterGolden(0x0),
  fTrackFilterTPC(0x0),
  fAnalysisType("AOD"),
  fCentDetector("V0M"),
  fAnalysisMC(kFALSE),
  fCentFrameworkAliCen(kFALSE),
  fAnalysisRun2(kFALSE),
  fVZEROBranch(kFALSE),
  fRandom(0x0),
  fEvent(0x0),
  fTrackArrayGlobalPar(0x0),
  fV0ArrayGlobalPar(0x0), 
  fTrackArrayMC(0x0),
  fVZEROArray(0x0),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  // fTrigType1(AliVEvent::kMB),
  // fTrigType2(AliVEvent::kCentral),
  // fTrigType3(AliVEvent::kSemiCentral),
  fVtxCut(10.0),  
  fEtaCut(0.8),  
  fEtaCutStack(1.2),  
  fMinPt(0.1),
  fMinPtV0(0.1),
  fCosPACut(0.95),
  fDecayRCut(0.2),
  fContributorsVtxCut(2),
  fContributorsVtxSPDCut(2),
  fZvsSPDvtxCorrCut(0.5),          
  fVtxR2Cut(10.0),           
  fCrossedRowsCut(70.0),     
  fCrossedOverFindableCut(0.8),
  fNegTrackStatus(0),
  fPosTrackStatus(0),
  fNegTOFExpTDiff(99999), 
  fPosTOFExpTDiff(99999),
  fRejectKinks(kTRUE),
  fLowPtFraction(0.01),
  fMassCut(0.1),
  fMinCent(0.0),
  fMaxCent(100.0),
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fTriggerInt(-999),
  fV0Finder(-999),
  fCentFramework(-999),
  fListOfObjects(0), 
  fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
  fn1(0),
  fn2(0),
  fcent(0),
  fTree(0x0)
{
  // Default constructor (should not be used)

  //  fRandom = new TRandom(0); // 0 means random seed
}

//______________________________________________________________________________
AliAnalysisTaskHighPtDeDx::AliAnalysisTaskHighPtDeDx(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(0),
  fAOD(0x0),
  fPIDResponse(0x0),
  fMC(0x0),
  fMCArray(0x0),
  fTrackFilter(0x0),
  fTrackFilterGolden(0x0),
  fTrackFilterTPC(0x0),
  fAnalysisType("AOD"),
  fCentDetector("V0M"),
  fAnalysisMC(kFALSE),
  fCentFrameworkAliCen(kFALSE),
  fAnalysisRun2(kFALSE),
  fVZEROBranch(kFALSE),
  fRandom(0x0),
  fEvent(0x0),
  fTrackArrayGlobalPar(0x0),
  fV0ArrayGlobalPar(0x0), //
  fTrackArrayMC(0x0),
  fVZEROArray(0x0),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  // fTrigType1(AliVEvent::kMB),
  // fTrigType2(AliVEvent::kCentral),
  // fTrigType3(AliVEvent::kSemiCentral),
  fVtxCut(10.0),  
  fEtaCut(0.8),
  fEtaCutStack(1.2),    
  fMinPt(0.1),
  fMinPtV0(0.1),
  fCosPACut(0.95),
  fDecayRCut(0.2),
  fContributorsVtxCut(2),
  fContributorsVtxSPDCut(2),
  fZvsSPDvtxCorrCut(0.5),          
  fVtxR2Cut(10.0),           
  fCrossedRowsCut(70.0),     
  fCrossedOverFindableCut(0.8),
  fNegTrackStatus(0),
  fPosTrackStatus(0),
  fNegTOFExpTDiff(99999), 
  fPosTOFExpTDiff(99999),
  fRejectKinks(kTRUE),
  fLowPtFraction(0.01),
  fMassCut(0.1),
  fMinCent(0.0),
  fMaxCent(100.0),
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fTriggerInt(-999),
  fV0Finder(-999),
  fCentFramework(-999),
  fListOfObjects(0), 
  fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
  fn1(0),
  fn2(0),
  fcent(0),
  fTree(0x0)
{

  if(fTree)fTree=0;
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliAnalysisTaskHighPtDeDx::~AliAnalysisTaskHighPtDeDx()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fListOfObjects;
    fListOfObjects = 0;
  }
  if (fRandom) delete fRandom;
  fRandom=0;
    
}

//______________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::UserCreateOutputObjects()
{ 
 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 
  // We also create the random generator here so it might get different seeds...
  fRandom = new TRandom(0); // 0 means random seed

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //
  // Histograms
  //  
  fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 3, 0, 3);
  fListOfObjects->Add(fEvents);

  fn1=new TH1F("fn1","fn1: trigger condition",5001,-1,5000);
  fListOfObjects->Add(fn1);

  fn2=new TH1F("fn2","fn2: trigger condition",5001,-1,5000);
  fListOfObjects->Add(fn2);

  fcent=new TH1F("fcent","fcent",104,-2,102);
  fListOfObjects->Add(fcent);

  fVtx = new TH1I("fVtx","Event Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
  fListOfObjects->Add(fVtx);

  fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);
  
  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

  if (fAnalysisMC) {    
    fVtxMC = new TH1I("fVtxMC","Vtx info - no trigger cut (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
    fListOfObjects->Add(fVtxMC);
  }

  //
  // Tree
  //  
  fTree = new TTree("tree","Event data");
  fEvent = new DeDxEvent();
  fTree->Branch("event", &fEvent);
  //array with tracks
  fTrackArrayGlobalPar = new TClonesArray("DeDxTrack", 1000);
  fTree->Bronch("trackGlobalPar", "TClonesArray", &fTrackArrayGlobalPar);
  //array with v0's
  fV0ArrayGlobalPar = new TClonesArray("DeDxV0", 1000);
  fTree->Bronch("v0GlobalPar", "TClonesArray", &fV0ArrayGlobalPar);
  
  if(fVZEROBranch){
    fVZEROArray = new TClonesArray("VZEROCell", 1000);
    fTree->Bronch("cellVZERO", "TClonesArray", &fVZEROArray);
  }
  if(fAnalysisMC) {    
    fTrackArrayMC = new TClonesArray("DeDxTrackMC", 1000);
    fTree->Bronch("trackMC", "TClonesArray", &fTrackArrayMC);
  }
  
  fTree->SetDirectory(0);
  fListOfObjects->Add(fTree);
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::UserExec(Option_t *) 
{
  // Main loop

  //
  // First we make sure that we have valid input(s)!
  //

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }
  
  if(fAnalysisRun2 && !fEventCuts.AcceptEvent(event)) return;
  
  //changes 26jun2018
  // fAOD = dynamic_cast<AliAODEvent*>(event);
  // if(!fAOD){
  //   Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
  //   this->Dump();
  //   return;
  //   }
  // }
  fAOD = dynamic_cast<AliAODEvent*>( InputEvent() );
  if (!fAOD) {
    AliWarning("ERROR: fAOD not available from InputEvent() trying with AODEvent()");
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      AliWarning("ERROR: fAOD not available from AODEvent() Aborting event!");
      return;
    }  
  }

  
  
  if (fAnalysisMC) { 

    fMC = dynamic_cast<AliMCEvent*>(MCEvent());
    
    fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
    if(!fMCArray){
      Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }    
    
  }//MC

  
  // Get trigger decision
  fTriggeredEventMB = 0; 

  // Bool_t lIsDesiredTrigger1 = AliMultSelectionTask::IsSelectedTrigger(event, fTrigType1);
  // if(lIsDesiredTrigger1)fTriggeredEventMB = 1;
 
  // Bool_t lIsDesiredTrigger2 = AliMultSelectionTask::IsSelectedTrigger(event, fTrigType2);
  // if(lIsDesiredTrigger2)fTriggeredEventMB += 2; 
  
  // Bool_t lIsDesiredTrigger3 = AliMultSelectionTask::IsSelectedTrigger(event, fTrigType3);
  // if(lIsDesiredTrigger3)fTriggeredEventMB += 4; 
  

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  if(inputHandler) fPIDResponse = inputHandler->GetPIDResponse();

  if(inputHandler->IsEventSelected() & ftrigBit1 ) fTriggeredEventMB = 1;

  //Addition to trigger x-check: fTriggerInt; // 0 = kMB, 1 = kCent, 2 = kSemiCent, 3 = kINT7
  if(inputHandler->IsEventSelected() & AliVEvent::kMB)          fTriggerInt = 0;
  if(inputHandler->IsEventSelected() & AliVEvent::kCentral)     fTriggerInt = 1; 
  if(inputHandler->IsEventSelected() & AliVEvent::kSemiCentral) fTriggerInt = 2;
  if(inputHandler->IsEventSelected() & AliVEvent::kINT7)        fTriggerInt = 3; 
  

  // Get process type for MC
  fMcProcessType = 0; // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD

  // real data that are not triggered we skip
  if(!fAnalysisMC && !fTriggeredEventMB) return; 
  

  if (fAnalysisMC) {
    
    AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader")); 
    
    if(mcHeader) {
      fZvtxMC = mcHeader->GetVtxZ();
      
      // PYTHIA:
      if(strstr(mcHeader->GetGeneratorName(), "Pythia")) {
	fMcProcessType =  GetPythiaEventProcessType(mcHeader->GetEventType());
      } else {
	// PHOJET:
	fMcProcessType =  GetDPMjetEventProcessType(mcHeader->GetEventType());
      }
    }
  }//MC
  
    //vtxstatus:
    //-2=badSPD
    //-1=noVtx
    //1=goodVtx
    //2=goodBut Zvtx-zvSPD>fZvsSPDvtxCorrCut (if cut value set)
    //3=goodBut radiusSq > fVtxR2Cut (if cut value set)
    //4=goodBut AliMultSelectionTask::IsINELgtZERO(fAOD) == kFALSE (use kTRUE for pp vs mult analyses and kFALSE for min bis)

  // We have to be careful how we normalize because we probably want to normalize to:
  // Nevents=(No vertex + outside + inside)/(out + in)*in
  
  
  // Global primary vertex
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  Double_t xVertex = vtx->GetX();
  Double_t yVertex = vtx->GetY();
  Double_t zVertex = vtx->GetZ();
  Double_t radiusSq = yVertex * yVertex + xVertex * xVertex;
  
  // SPD primary vertex
  const AliAODVertex *vtxSPD = fAOD->GetPrimaryVertexSPD();
  Float_t zvSPD = vtxSPD->GetZ();


  if(!vtx){
    fVtxStatus = -999;
    fVtxStatus = -1;
  }
  
  if(vtx->GetNContributors() <= fContributorsVtxCut){
    fZvtx = -1599;
    fVtxStatus = -2;
  }
  else if(vtxSPD->GetNContributors() <= fContributorsVtxSPDCut){
    fZvtx = -1599;
    fVtxStatus = -2;
  }
  else if(TMath::Abs( fZvtx - zvSPD ) > fZvsSPDvtxCorrCut){
    fZvtx = zVertex;
    fVtxStatus = 2;
  }
  else if(radiusSq > fVtxR2Cut){
    fZvtx = zVertex;
    fVtxStatus = 3;
  }
  else{
    fZvtx = zVertex;
    fVtxStatus = 1;
  }
  
  
  Float_t centrality = -10;

  // only analyze triggered events
  if(fTriggeredEventMB){
    
    if(fCentFrameworkAliCen){ // CentFramework in AliCentrality
      fCentFramework = 0; // x-check histo: 0 = AliCentrality, 1 = AliMultSelection
      
      AliCentrality *centObject =  fAOD->GetCentrality();
      if(centObject) centrality = centObject->GetCentralityPercentile(fCentDetector); 
      
    }else{ // CentFramework in AliMultSelection
      fCentFramework = 1; // x-check histo: 0 = AliCentrality, 1 = AliMultSelection
      
      AliMultSelection* centObject = 0x0; 
      centObject = (AliMultSelection*)fAOD->FindListObject("MultSelection");

      if(AliMultSelectionTask::IsINELgtZERO(fAOD) == kFALSE) fVtxStatus = 4.;
      
      if(centObject){
	centrality = centObject->GetMultiplicityPercentile(fCentDetector);	
      }
      if(!centObject) cout<<"no centObject: please check that the AliMultSelectionTask actually ran (before your task) "<<endl; 
    }
    
    if((centrality>fMaxCent)||(centrality<fMinCent))return;	
    
    fcent->Fill(centrality);
    AnalyzeAOD(fAOD);
     
  }//if triggered

  
  // store MC event data
  if(fAnalysisMC) ProcessMCTruthAOD();
  
  fEvent->process       = fMcProcessType;
  fEvent->trig          = fTriggeredEventMB;
  fEvent->triggerInt    = fTriggerInt;
  fEvent->zvtxMC        = fZvtxMC;
  fEvent->cent          = centrality;
  fEvent->centFramework = fCentFramework;

  //fEvent->centV0A      = centralityV0A;
  //fEvent->centZNA      = centralityZNA;
  //fEvent->centCL1      = centralityCL1;
  
  fTree->Fill();
  
  if(fTrackArrayGlobalPar)fTrackArrayGlobalPar->Clear();
  if(fV0ArrayGlobalPar)fV0ArrayGlobalPar->Clear();
  if(fVZEROArray)fVZEROArray->Clear();
  if(fAnalysisMC)fTrackArrayMC->Clear();
  
  // Post output data.
  PostData(1, fListOfObjects);
}


//________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::AnalyzeAOD(AliAODEvent* aodEvent)
{
  fRun  = aodEvent->GetRunNumber();
  fEventId = 0;
  if(aodEvent->GetHeader())
    fEventId = GetEventIdAsLong(aodEvent->GetHeader());
  
  UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();
  Float_t   magf      = aodEvent->GetMagneticField();

  if(fTriggeredEventMB) {// Only MC case can we have not triggered events
    
    // accepted event
    fEvents->Fill(0);
    
    //if(fVtxStatus!=1) return; // accepted vertex
    //    Int_t nAODTracks = aodEvent->GetNumberOfTracks();

    ProduceArrayTrksAOD( aodEvent, kGlobalTrk );
    ProduceArrayV0AOD( aodEvent, kGlobalTrk );

    fEvents->Fill(1);

    if(fVZEROBranch){//change 10/07-18: esdV0->aodV0
      AliAODVZERO *aodV0 = aodEvent->GetVZEROData();
      for (Int_t iCh=0; iCh<64; ++iCh) { 
	Float_t multv=aodV0->GetMultiplicity(iCh); 
	Int_t intexv=iCh;
	VZEROCell* cellv0 = new((*fVZEROArray)[iCh]) VZEROCell();
	cellv0->cellindex=intexv;
	cellv0->cellmult= multv;
      }   
    }

  } // end if triggered
  
  fEvent->run       = fRun;
  fEvent->eventid   = fEventId;
  fEvent->time      = time;
  fEvent->mag       = magf;
  fEvent->zvtx      = fZvtx;
  fEvent->vtxstatus = fVtxStatus;

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskHighPtDeDx::GetVertex(const AliVEvent* event) const
{
  Float_t zvtx = -999;
  
  const AliVVertex* primaryVertex = event->GetPrimaryVertex(); 
  
  if(primaryVertex->GetNContributors()>0)
    zvtx = primaryVertex->GetZ();

  return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskHighPtDeDx::GetPidCode(Int_t pdgCode) const 
{
  // return our internal code for pions, kaons, and protons
  
  Short_t pidCode = 6;
  
  switch (TMath::Abs(pdgCode)) {
  case 211:
    pidCode = 1; // pion
    break;
  case 321:
    pidCode = 2; // kaon
    break;
  case 2212:
    pidCode = 3; // proton
    break;
  case 11:
    pidCode = 4; // electron
    break;
  case 13:
    pidCode = 5; // muon
    break;
  default:
    pidCode = 6;  // something else?
  };
  
  return pidCode;
}

//_____________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::ProcessMCTruthAOD() 
{

  Short_t trackmult = 0;
  Short_t nadded    = 0;
  const Int_t nTracksMC = fMCArray->GetEntriesFast();

  for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

    //Cuts
    if(!(trackMC->IsPhysicalPrimary())) continue;
    Double_t chargeMC = trackMC->Charge();

    
    if (TMath::Abs(trackMC->Eta()) > fEtaCutStack ) continue;
    trackmult++;
    
    Float_t ptMC      = trackMC->Pt();
    Float_t pMC       = trackMC->P();
    Float_t etaMC     = trackMC->Eta();
    Float_t phiMC     = trackMC->Phi();
    Float_t yMC       = trackMC->Y();
    Int_t pdgCode     = trackMC->PdgCode();
    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(pdgCode);
   
    
    TString genname;
    Bool_t yesno=fMC->GetCocktailGenerator(iTracks, genname);
    if(yesno) {
      if(!genname.Contains("Hijing")) 
	pidCodeMC +=10;
    }
  
 
    Bool_t lIsStrangeness = kFALSE; 
    if ( TMath::Abs(pdgCode)==310 || TMath::Abs(pdgCode)==3122 || TMath::Abs(pdgCode)==3312 || TMath::Abs(pdgCode)==3334 ) lIsStrangeness = kTRUE; 

    if (trackMC->Pt() < fMinPt && !lIsStrangeness) {
      // Keep small fraction of low pT tracks
      if(fRandom->Rndm() > fLowPtFraction) continue; 
    }     
    
    if (trackMC->Pt() < fMinPtV0 && lIsStrangeness) {
      // Keep small fraction of low pT tracks
      if(fRandom->Rndm() > fLowPtFraction) continue; 
    }
    

    DeDxTrackMC* track = new((*fTrackArrayMC)[nadded]) DeDxTrackMC();
    nadded++;
    track->pMC   = pMC;
    track->ptMC  = ptMC;
    track->etaMC = etaMC;
    track->phiMC = phiMC;
    track->yMC   = yMC;
    track->qMC   = Short_t(chargeMC);
    track->pidMC = pidCodeMC;
    track->pdgMC = pdgCode;
  }
      
  Sort(fTrackArrayMC, kTRUE);

  fEvent->trackmultMC = trackmult;
  fEvent->nMC         = nadded;

}

//_____________________________________________________________________________
Short_t AliAnalysisTaskHighPtDeDx::GetPythiaEventProcessType(Int_t pythiaType) {
  //
  // Get the process type of the event.  PYTHIA
  //
  // source PWG0   dNdpt 

  Short_t globalType = -1; //init
      
  if(pythiaType==92||pythiaType==93){
    globalType = 2; //single diffractive
  }
  else if(pythiaType==94){
    globalType = 3; //double diffractive
  }
  //else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
  else {
    globalType = 1;  //non diffractive
  }
  return globalType;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskHighPtDeDx::GetDPMjetEventProcessType(Int_t dpmJetType) {
  //
  // get the process type of the event.  PHOJET
  //
  //source PWG0   dNdpt 
  // can only read pythia headers, either directly or from cocktalil header
  Short_t globalType = -1;
  
  if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
    globalType = 1;
  }
  else if (dpmJetType==5 || dpmJetType==6) {
    globalType = 2;
  }
  else if (dpmJetType==7) {
    globalType = 3;
  }
  return globalType;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskHighPtDeDx::GetEventIdAsLong(AliVHeader* header) const
{
  // To have a unique id for each event in a run!
  // Modified from AliRawReader.h
  return ((ULong64_t)header->GetBunchCrossNumber()+
	  (ULong64_t)header->GetOrbitNumber()*3564+
	  (ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskHighPtDeDx::FindPrimaryMotherAOD(AliAODMCParticle* startParticle)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart)
    {
      
      if(mcPart->IsPrimary())
	return mcPart;
      
      Int_t mother = mcPart->GetMother();

      mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
    }

  return 0;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskHighPtDeDx::FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  nSteps = 0;

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart)
    {
      
      if(mcPart->IsPrimary())
	return mcPart;
      
      Int_t mother = mcPart->GetMother();

      mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
      nSteps++; // 1 level down
    }

  return 0;
}
 
//_____________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::Sort(TClonesArray* array, Bool_t isMC) 
{
  const Int_t n = array->GetEntriesFast(); 

  if(n==0) {

    if(isMC) 
      fEvent->ptmaxMC = 0;
    else
      fEvent->ptmax   = 0;
      
    return;
  }

  Float_t ptArray[n];
  Int_t   indexArray[n];
  
  for(Int_t i = 0; i < n; i++) {

    Float_t pt = 0;
    if(isMC) {
      DeDxTrackMC* track = (DeDxTrackMC*)array->At(i);
      pt = track->ptMC;
    } else {
      DeDxTrack* track = (DeDxTrack*)array->At(i);
      pt = track->pt;
    }
    ptArray[i]    = pt;
    indexArray[i] = i;
  }

  TMath::Sort(n, ptArray, indexArray);
  
  // set max pt
  if(isMC) {
    DeDxTrackMC* track = (DeDxTrackMC*)array->At(indexArray[0]);
    fEvent->ptmaxMC = track->ptMC;
  } else {
    DeDxTrack* track = (DeDxTrack*)array->At(indexArray[0]);
    fEvent->ptmax   = track->pt;
  }
  
  // set order of each track
  for(Int_t i = 0; i < n; i++) {
    
    if(isMC) {
      DeDxTrackMC* track = (DeDxTrackMC*)array->At(indexArray[i]);
      track->orderMC = i;
    } else {
      DeDxTrack* track = (DeDxTrack*)array->At(indexArray[i]);
      track->order = i;
    }
  }
}




// ################################
// ############ AOD - updated ##### 
// ################################
//__________________________________________________________________
void AliAnalysisTaskHighPtDeDx::ProduceArrayTrksAOD( AliAODEvent *AODevent, AnalysisMode analysisMode ){

  Int_t trackmult = 0; // no pt cuts
  Int_t nadded = 0;
  Int_t nAODTracks = AODevent->GetNumberOfTracks();

  //31/10-18
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(AODevent->GetHeader());
  UInt_t iNTracksRef = 0;
  if(!header) AliError("Cannot get AOD header");
  else iNTracksRef = header->GetRefMultiplicityComb08(); // get reference multiplicity

  
  if(fTrackArrayGlobalPar) fTrackArrayGlobalPar->Clear();
    
  //const AliAODVertex*	vertexSPD= (AliAODVertex*)AODevent->GetPrimaryVertexSPD();//GetPrimaryVertex()
  //if( vertexSPD->GetNContributors() < 1 || TMath::Abs( vertexSPD->GetZ()) > 10.0 ) return;
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AddTaskInfoAOD145
  
  const AliAODVertex* vertex = AODevent->GetPrimaryVertex();

  for(Int_t iT = 0; iT < nAODTracks; iT++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(AODevent->GetTrack(iT));
    if(!aodTrack) AliFatal("Not a standard AOD");
      
    //FOR AOD068, filterCut_Set2=KTRUE WHEN THE TRACK SATISFIES THE CUTS FROM JET ANALYSIS
      
    UShort_t filterFlag = 0;
    
    // "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
    // if (fTrackFilterGolden) {
    //   if(aodTrack->TestFilterBit(1024)) filterFlag +=1;
    // }
    
    // FILTER 128 ARE TPC ONLY TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX
    // if (fTrackFilterTPC) {
    //   if(aodTrack->TestFilterBit(128)) filterFlag +=2;
    // }
    
    // tighter cuts on primary particles for high pT tracks
    // take the standard cuts, which include already 
    // ITSrefit and use only primaries...
    // if (fTrackFilter) {
    //   if(aodTrack->TestFilterBit(32)) filterFlag +=4;
    // }
    
    //https://twiki.cern.ch/twiki/bin/view/ALICE/HybridTracks
    //PbPb LHC10h	160	768	 
    //PbPb LHC11h	145	768
    // if(aodTrack->TestFilterBit(768)) filterFlag +=8; //before, when used with the other filterbit options
    
    if(aodTrack->TestFilterBit(768)) filterFlag +=1; 

    if(!aodTrack->TestFilterBit(768)) continue;
    
    //new method 23 Jan 2019:
    
    if(!aodTrack->IsOn(AliAODTrack::kTPCrefit)) continue;

    if(aodTrack->Chi2perNDF()>4) continue;

    Float_t nCrossedRowsTPC = aodTrack->GetTPCClusterInfo(2,1);
    Int_t findable = aodTrack->GetTPCNclsF();
    if(nCrossedRowsTPC < fCrossedRowsCut) continue;
    if(findable <= 0) continue;
    if(nCrossedRowsTPC/findable < fCrossedOverFindableCut) continue;        
            
    
    Double_t b[2], cov[3];
    if (!(aodTrack->PropagateToDCA(vertex, AODevent->GetMagneticField(), kVeryBig, b, cov))) filterFlag = 32; // propagation failed!!!!!; only need to use if filterbit is 128
    
    Float_t dcaxy = b[0];
    Float_t dcaz = b[1];
      
    TBits sharedTPC = aodTrack->GetTPCSharedMap();
    Int_t sharedtpcclusters = sharedTPC.CountBits(0)-sharedTPC.CountBits(159);
 
    if(TMath::Abs(aodTrack->Eta()) > fEtaCut) continue;

    if (aodTrack->Pt() < fMinPt) {
      // Keep small fraction of low pT tracks
      if(fRandom->Rndm() > fLowPtFraction) continue;
    }
      
    Short_t charge      = aodTrack->Charge();
    Float_t pt          = aodTrack->Pt();
    Float_t p           = aodTrack->P(); 
    Float_t eta         = aodTrack->Eta();
    Float_t phi         = aodTrack->Phi();
    AliAODPid* aodPid   = aodTrack->GetDetPid();
    Short_t neff        = 0; // neff is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res -> see https://github.com/alisw/AliPhysics/blob/master/PWGLF/STRANGENESS/Correlationpp/AliAnalysisTaskCorrelationhhK0s.cxx for implementation 
    Short_t ncl         = -10;
    Float_t dedx        = -10;
    if(aodPid){
      ncl     = aodPid->GetTPCsignalN();
      dedx    = aodPid->GetTPCsignal();
    }

    Float_t ptMC        = 0;
    Short_t pidCode     = 0; // 0 = real data / no mc track!
    Short_t primaryFlag = 0; // 0 = real data / not primary mc track  
    Int_t pdgMother     = 0;

    
    if(fAnalysisMC) {
	
      const Int_t label = TMath::Abs(aodTrack->GetLabel());
	
      AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));

      if (mcTrack){

	if(mcTrack->IsPhysicalPrimary()) primaryFlag = 1;
	  
	Int_t pdgCode = mcTrack->GetPdgCode();
	pidCode = GetPidCode(pdgCode);
	  
	TString genname;
	Bool_t yesno = fMC->GetCocktailGenerator(label, genname);
	if(yesno){
	  if(!genname.Contains("Hijing")) pidCode +=10;
	}
	// if(!yesno) Printf("no cocktail header list was found for this event");
	
	ptMC = mcTrack->Pt();
	AliAODMCParticle* mother = FindPrimaryMotherAOD(mcTrack);
	pdgMother = mother->GetPdgCode();	    

      }//mcTrack

    }//MC
    
    DeDxTrack* track = new((*fTrackArrayGlobalPar)[nadded]) DeDxTrack();
    nadded++;
	
    track->p          = p;
    track->pt         = pt;
    track->eta        = eta;
    track->phi        = phi;
    track->q          = charge;
    track->filter     = filterFlag;
    track->ncl        = ncl;
    track->neff       = neff;
    track->dedx       = dedx;
    track->protNSigma = 0;
    track->pionNSigma = 0;
    track->dcaxy      = dcaxy;
    track->dcaz       = dcaz;
    track->pid        = pidCode;
    track->primary    = primaryFlag;
    track->pttrue     = ptMC;
    track->mother     = pdgMother;
    track->tpcnclS    = sharedtpcclusters;

  }//end of track loop

  Sort(fTrackArrayGlobalPar, kFALSE);
  fEvent->trackmult  = trackmult;
  fEvent->n          = nadded;
  fEvent->nTracks    = nAODTracks;//13/08-18
  fEvent->refMult    = iNTracksRef;//31/10-18
}

// ################################
// #### BEGINNING OF V0 CODE ######
// ################################
void AliAnalysisTaskHighPtDeDx::ProduceArrayV0AOD( AliAODEvent *AODevent, AnalysisMode analysisMode ){
  Int_t nv0s = AODevent->GetNumberOfV0s();
  if(nv0s<1)return;
  Int_t     nadded    = 0;
  
  AliAODVertex *myBestPrimaryVertex = AODevent->GetPrimaryVertex();
  if(!myBestPrimaryVertex) return;
  
  if(fV0ArrayGlobalPar) fV0ArrayGlobalPar->Clear();
  
  // This is the begining of the V0 loop  
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
    AliAODv0 *aodV0 = AODevent->GetV0(iV0);
    if (!aodV0) continue;
    
    AliAODVertex* vertex = aodV0->GetSecondaryVtx();
    if(!vertex){
      Printf("ERROR: Could not retrieve vertex");
      continue;
    }
      
    AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0); 
    AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
    
    if(!pTrack || !nTrack){
      Printf("ERROR: Could not retrieve one of the daughter track");
      continue;
    }
    
    // Remove like-sign daughters
    if(pTrack->Charge() == nTrack->Charge()) continue; 
    
    // Make sure charge ordering is ok
    if(pTrack->Charge() < 0) {
      AliAODTrack* helpTrack = pTrack;
      pTrack = nTrack;
      nTrack = helpTrack;
    } 
        
    //reject kink daughters //06 Feb 2019:
    AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*)(nTrack->GetProdVertex());
    Char_t cTypeVtxProdNeg = prodVtxDaughterNeg->GetType(); 
    AliAODVertex* prodVtxDaughterPos = (AliAODVertex*)(pTrack->GetProdVertex());
    Char_t cTypeVtxProdPos = prodVtxDaughterPos->GetType(); 
    
    if(fRejectKinks){
      if(cTypeVtxProdNeg == AliAODVertex::kKink) continue;
      if(cTypeVtxProdPos == AliAODVertex::kKink) continue;
    }


    // neff is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res -> see https://github.com/alisw/AliPhysics/blob/master/PWGLF/STRANGENESS/Correlationpp/AliAnalysisTaskCorrelationhhK0s.cxx for implementation 
 
    //#### check positive tracks ####
    
    UShort_t filterFlag_p = 0;
    // if(fTrackFilterGolden){
    //   if(pTrack->TestFilterBit(1024)) filterFlag_p +=1;
    // }
    // if(fTrackFilterTPC){
    //   if(pTrack->TestFilterBit(128)) filterFlag_p +=2;	  
    // }
    // if (fTrackFilter){
    //   if(pTrack->TestFilterBit(32)) filterFlag_p +=4;
    // }
    // if(pTrack->TestFilterBit(768)) filterFlag_p +=8; 


    if(!pTrack->IsOn(AliAODTrack::kTPCrefit)) continue;
    
    Float_t nCrossedRowsTPC_p = pTrack->GetTPCClusterInfo(2,1);
    Int_t findable_p = pTrack->GetTPCNclsF();
    if(nCrossedRowsTPC_p < fCrossedRowsCut) continue;
    if(findable_p <= 0) continue;
    if(nCrossedRowsTPC_p/findable_p < fCrossedOverFindableCut) continue;        
    
    if(pTrack->Chi2perNDF()>4) continue;
       
    //#### check negative tracks ####

    UShort_t filterFlag_n = 0;
    // if(fTrackFilterGolden){
    //   if(nTrack->TestFilterBit(1024)) filterFlag_n +=1;
    // }
    // if(fTrackFilterTPC){
    //   if(nTrack->TestFilterBit(128)) filterFlag_n +=2;	  
    // }
    // if (fTrackFilter){
    //   if(nTrack->TestFilterBit(32)) filterFlag_n +=4;
    // }
    // if(nTrack->TestFilterBit(768)) filterFlag_n +=8; 


    if(!nTrack->IsOn(AliAODTrack::kTPCrefit)) continue;
    
    Float_t nCrossedRowsTPC_n = nTrack->GetTPCClusterInfo(2,1);
    Int_t findable_n = nTrack->GetTPCNclsF();
    if(nCrossedRowsTPC_n < fCrossedRowsCut) continue;
    if(findable_n <= 0) continue;
    if(nCrossedRowsTPC_n/findable_n < fCrossedOverFindableCut) continue;        
    
    if(nTrack->Chi2perNDF()>4) continue;


    //pileup rejection:
    
    Short_t oobPileupFlag = 0; // 1 = oob pileup rejection cut applied
    
    fPosTrackStatus = pTrack->GetStatus();
    fNegTrackStatus = nTrack->GetStatus();
    fNegTOFExpTDiff = nTrack->GetTOFExpTDiff(AODevent->GetMagneticField());
    fPosTOFExpTDiff = pTrack->GetTOFExpTDiff(AODevent->GetMagneticField());

    //ITS||TOF requirement
    Bool_t lITSorTOFsatisfied = kFALSE; 
    if( (fNegTrackStatus & nTrack->IsOn(AliAODTrack::kITSrefit) ||
	 (fPosTrackStatus & pTrack->IsOn(AliAODTrack::kITSrefit)) ) )
      lITSorTOFsatisfied = kTRUE; 
    if( (TMath::Abs(fNegTOFExpTDiff+2500.) > 1e-6) || (TMath::Abs(fPosTOFExpTDiff+2500.) > 1e-6) )
      lITSorTOFsatisfied = kTRUE;
    
    if(lITSorTOFsatisfied) oobPileupFlag = 1;
    
    
    Double_t lV0Radius         = aodV0->RadiusV0();
    Float_t alpha              = aodV0->AlphaV0();
    Float_t ptarm              = aodV0->PtArmV0();
    Double_t lV0DecayLength    = aodV0->DecayLength(myBestPrimaryVertex);
    Double_t deltaInvMassG     = aodV0->InvMass2Prongs(0,1,11,11);
    Double_t deltaInvMassK0s   = aodV0->MassK0Short()-0.498;
    Double_t deltaInvMassL     = aodV0->MassLambda()-1.116;
    Double_t deltaInvMassAntiL = aodV0->MassAntiLambda()-1.116;
      
    TBits psharedTPC = pTrack->GetTPCSharedMap();
    Int_t psharedtpcclusters = psharedTPC.CountBits(0)-psharedTPC.CountBits(159);
    Short_t pcharge  = pTrack->Charge();
    Float_t ppt      = pTrack->Pt();
    Float_t pp       = pTrack->P(); 
    Float_t peta     = pTrack->Eta();
    Float_t pphi     = pTrack->Phi();
    AliAODPid* pPid  = pTrack->GetDetPid();
    Short_t pncl     = -10;
    Short_t pneff    = 0;
    Float_t pProtNSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton ));
    Float_t pPionNSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion ));
    Float_t pdedx    = -10;
    if(pPid){
      pncl     = pPid->GetTPCsignalN();
      pdedx    = pPid->GetTPCsignal(); 
    }
    
    TBits nsharedTPC = nTrack->GetTPCSharedMap();
    Int_t nsharedtpcclusters = nsharedTPC.CountBits(0)-nsharedTPC.CountBits(159);
    Short_t ncharge  = nTrack->Charge();
    Float_t npt      = nTrack->Pt();
    Float_t np       = nTrack->P(); 
    Float_t neta     = nTrack->Eta();
    Float_t nphi     = nTrack->Phi();
    AliAODPid* nPid  = nTrack->GetDetPid();
    Short_t nncl     = -10;
    Short_t nneff    = 0; 
    Float_t nProtNSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton ));
    Float_t nPionNSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion ));
    Float_t ndedx    = -10;
    if(nPid){
      nncl     = nPid->GetTPCsignalN();
      ndedx    = nPid->GetTPCsignal(); 
    }

    
    // ### Pre-selection to reduce output size ###
    if(TMath::Abs(peta) > fEtaCut || TMath::Abs(neta) > fEtaCut) continue;
    // if(TMath::Abs(aodV0->Eta()) > fEtaCut) continue;
    if(aodV0->GetOnFlyStatus() != 0 ) continue; //Reject on-the-fly tracks
    if(lV0Radius < fDecayRCut ) continue;
    if(aodV0->Pt() < fMinPtV0) continue;
    if(aodV0->CosPointingAngle(myBestPrimaryVertex) < fCosPACut ) continue;
    if(TMath::Abs(deltaInvMassK0s) > fMassCut &&
       TMath::Abs(deltaInvMassL) > fMassCut &&
       TMath::Abs(deltaInvMassAntiL) > fMassCut &&
       TMath::Abs(deltaInvMassG) > fMassCut )
      continue;

    Int_t   primaryV0     = 0; // 0 means that the tracks are not both daughters of a primary particle (1 means they are)
    Int_t   pdgV0         = 0; // 0 means that they don't have same origin for MC (1 means they have the same original mother)
    Int_t   pdgmotherV0   = 0; 

    Float_t p_ptMC        = 0;
    Short_t p_pidCode     = 0; // 0 = real data / no mc track!
    Short_t p_primaryFlag = 0; // 0 = real data / not primary mc track  
    Int_t   p_pdgMother   = 0;

    Float_t n_ptMC        = 0;
    Short_t n_pidCode     = 0; // 0 = real data / no mc track!
    Short_t n_primaryFlag = 0; // 0 = real data / not primary mc track  
    Int_t   n_pdgMother   = 0;

    if(fAnalysisMC) {
      // The following code was modified by Peter (18/8-2015)
      // This is the AOD MC analysis of the daughter tracks

      AliAODMCParticle* p_mother = 0;
      Int_t p_mother_steps       = 0;
      Int_t p_mother_label       = 0;     
      AliAODMCParticle* n_mother = 0;
      Int_t n_mother_steps       = 0;
      Int_t n_mother_label       = 0;
	
      // positive track
      const Int_t p_label = TMath::Abs(pTrack->GetLabel());
      AliAODMCParticle* p_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(p_label));

      if (p_mcTrack){
	
	if(p_mcTrack->IsPhysicalPrimary()) p_primaryFlag = 1;
	
	Int_t p_pdgCode = p_mcTrack->GetPdgCode();
	p_pidCode = GetPidCode(p_pdgCode);	  
	p_ptMC = p_mcTrack->Pt();

	p_mother_label = p_mcTrack->GetMother(); 	  
	if(p_mother_label >= 0) {
	  p_mother = dynamic_cast<AliAODMCParticle*>(fMCArray->At(p_mother_label));
	  p_pdgMother = p_mother->GetPdgCode();
	}
      }

      // negative track
      const Int_t n_label = TMath::Abs(nTrack->GetLabel());
      AliAODMCParticle* n_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(n_label));

      if (n_mcTrack){
	  
	if(n_mcTrack->IsPhysicalPrimary()) n_primaryFlag = 1;
	  
	Int_t n_pdgCode = n_mcTrack->GetPdgCode();
	n_pidCode = GetPidCode(n_pdgCode);	  
	n_ptMC = n_mcTrack->Pt();

	n_mother_label = n_mcTrack->GetMother(); 	  
	if(n_mother_label >= 0) {
	  n_mother = dynamic_cast<AliAODMCParticle*>(fMCArray->At(n_mother_label));
	  n_pdgMother = n_mother->GetPdgCode();
	}
      }
      
      // Check if V0 is primary = first and the same mother of both partciles
      if(p_mother && n_mother && p_mother == n_mother) {
	
	TString genname;
	Bool_t yesno=fMC->GetCocktailGenerator(p_mother_label, genname);
	//    if(!yesno) Printf("no cocktail header list was found for this event");
	if(yesno) {
	  if(!genname.Contains("Hijing")) {
	    p_pidCode +=10;
	    n_pidCode +=10;
	  }
	}
	  
	pdgV0 = p_pdgMother;
	if(p_mother->IsPhysicalPrimary()){
	  primaryV0 = 1;
	}
      }
    }//MC
    
      
    DeDxV0* v0data = new((*fV0ArrayGlobalPar)[nadded]) DeDxV0();
    nadded++;
	
    // v0 data
    v0data->p            = aodV0->P();
    v0data->pt           = aodV0->Pt();
    v0data->eta          = aodV0->Eta();
    v0data->phi          = aodV0->Phi();
    v0data->pdca         = aodV0->DcaPosToPrimVertex();
    v0data->ndca         = aodV0->DcaNegToPrimVertex();
    v0data->dmassG       = deltaInvMassG;
    v0data->dmassK0      = deltaInvMassK0s;
    v0data->dmassL       = deltaInvMassL;
    v0data->dmassAL      = deltaInvMassAntiL;
    v0data->alpha        = alpha;
    v0data->ptarm        = ptarm;
    v0data->decayr       = lV0Radius;
    v0data->decayl       = lV0DecayLength;
    // v0data->pdca      = TMath::Sqrt(pdcaxy*pdcaxy + pdcaz*pdcaz);
    // v0data->ndca      = TMath::Sqrt(ndcaxy*ndcaxy + ndcaz*ndcaz);
    v0data->status       = aodV0->GetOnFlyStatus();
    v0data->chi2         = aodV0->Chi2V0();
    v0data->cospt        = aodV0->CosPointingAngle(myBestPrimaryVertex);
    v0data->dcav0        = aodV0->DcaV0ToPrimVertex();
    v0data->dcadaughters = aodV0->DcaV0Daughters();
    v0data->primary      = primaryV0;
    v0data->pdg          = pdgV0;
    v0data->pdgmother    = pdgmotherV0;
    v0data->oobPileupFlag = oobPileupFlag;
	 
	
    // positive track
    v0data->ptrack.p       = pp;
    v0data->ptrack.pt      = ppt;
    v0data->ptrack.eta     = peta;
    v0data->ptrack.phi     = pphi;
    v0data->ptrack.q       = pcharge;
    v0data->ptrack.ncl     = pncl;
    v0data->ptrack.neff    = pneff;
    v0data->ptrack.protNSigma = pProtNSigma;
    v0data->ptrack.pionNSigma = pPionNSigma;
    v0data->ptrack.dedx    = pdedx;
    v0data->ptrack.pid     = p_pidCode;
    v0data->ptrack.primary = p_primaryFlag;
    v0data->ptrack.pttrue  = p_ptMC;
    v0data->ptrack.mother  = p_pdgMother;
    v0data->ptrack.filter  = filterFlag_p;
    v0data->ptrack.tpcnclS = psharedtpcclusters;
	
    // negative track
    v0data->ntrack.p       = np;
    v0data->ntrack.pt      = npt;
    v0data->ntrack.eta     = neta;
    v0data->ntrack.phi     = nphi;
    v0data->ntrack.q       = ncharge;
    v0data->ntrack.ncl     = nncl;
    v0data->ntrack.neff    = nneff;
    v0data->ntrack.protNSigma = nProtNSigma;
    v0data->ntrack.pionNSigma = nPionNSigma;
    v0data->ntrack.dedx    = ndedx;
    v0data->ntrack.pid     = n_pidCode;
    v0data->ntrack.primary = n_primaryFlag;
    v0data->ntrack.pttrue  = n_ptMC;
    v0data->ntrack.mother  = n_pdgMother;
    v0data->ntrack.filter  = filterFlag_n;
    v0data->ntrack.tpcnclS = nsharedtpcclusters;
	
  }//end loop over v0's
  
}// ###### END AOD #######



