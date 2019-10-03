#include "AliAnalysisTaskHighPtDeDxV0.h"

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
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include <AliAODVertex.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>

#include "AliCentrality.h" 
#include "AliESDtrackCuts.h"

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 

// STL includes
#include <iostream>
using namespace std;

//_____________________________________________________________________________
//
// Responsible:
// Peter Christiansen (Lund)
//
//_____________________________________________________________________________


ClassImp(AliAnalysisTaskHighPtDeDxV0)

const Double_t AliAnalysisTaskHighPtDeDxV0::fgkClight = 2.99792458e-2;

//_____________________________________________________________________________
AliAnalysisTaskHighPtDeDxV0::AliAnalysisTaskHighPtDeDxV0():
  AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
  fTrackFilter(0x0),
  fTrackFilterGolden(0x0),
  fTrackFilterTPC(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  fEvent(0x0),
  fV0ArrayGlobalPar(0x0),
  fV0ArrayTPCPar(0x0),
  fTrackArrayMC(0x0),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinPt(0.1),
  fMinCent(0.0),
  fMaxCent(100.0),
  fMassCut(0.1),
  fTreeOption(0),
  fRequireRecV0(kTRUE),
  fStoreMcIn(kFALSE),
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0), 
  fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
  fTree(0x0),
  fnv0(0)
{
  // Default constructor (should not be used)

}

//______________________________________________________________________________
AliAnalysisTaskHighPtDeDxV0::AliAnalysisTaskHighPtDeDxV0(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
  fTrackFilter(0x0),
  fTrackFilterGolden(0x0),
  fTrackFilterTPC(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  fEvent(0x0),
  fV0ArrayGlobalPar(0x0),
  fV0ArrayTPCPar(0x0),
  fTrackArrayMC(0x0),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinPt(0.1),
  fMinCent(0.0),
  fMaxCent(100.0),
  fMassCut(0.1),
  fTreeOption(0),
  fRequireRecV0(kTRUE),
  fStoreMcIn(kFALSE),
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0), 
  fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
  fTree(0x0),
  fnv0(0)
{
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskHighPtDeDxV0::~AliAnalysisTaskHighPtDeDxV0()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fListOfObjects;
    fListOfObjects = 0;
  }
  
  // //for proof running; I cannot create tree do to memory limitations -> work with THnSparse 
  // if (fListOfObjects  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
  
  
  
}

//______________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::UserCreateOutputObjects()
{ 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //
  // Histograms
  //  
  fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 1, 0, 1);
  fListOfObjects->Add(fEvents);
  
  fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
  fListOfObjects->Add(fVtx);

  fnv0 = new TH1F("fnv0","# V0's", 10000,0,10000);
  fListOfObjects->Add(fnv0);

  fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);
  
  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

  if (fAnalysisMC) {    
    fVtxMC = new TH1I("fVtxMC","Vtx info - no trigger cut (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
    fListOfObjects->Add(fVtxMC);
  
  }

  if (fTreeOption) {

    fTree = new TTree("tree","Event data");
    fEvent = new DeDxEvent();
    fTree->Branch("event", &fEvent);

    fV0ArrayGlobalPar = new TClonesArray("DeDxV0", 1000);
    fTree->Bronch("v0GlobalPar", "TClonesArray", &fV0ArrayGlobalPar);

    fV0ArrayTPCPar = new TClonesArray("DeDxV0", 1000);
    fTree->Bronch("v0TPCPar", "TClonesArray", &fV0ArrayTPCPar);

    if (fAnalysisMC && fStoreMcIn) {    
      fTrackArrayMC = new TClonesArray("DeDxTrackMC", 1000);
      fTree->Bronch("trackMC", "TClonesArray", &fTrackArrayMC);
    }

    fTree->SetDirectory(0);

    fListOfObjects->Add(fTree);

  }

  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::UserExec(Option_t *) 
{
  // Main loop
  
  //
  // Comment: This method matches completely the same method for the high pT
  // tracks
  //


  //
  // First we make sure that we have valid input(s)!
  //
  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }    
  } else {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }    
  }

  if (fAnalysisMC) {

    if (fAnalysisType == "ESD"){
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(!fMC){
	Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }    

      fMCStack = fMC->Stack();
      
      if(!fMCStack){
	Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }    
    } else { // AOD

      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(fMC)
	fMC->Dump();

      fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
      if(!fMCArray){
	Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
      }    
    }
  }

    // Get trigger decision
  fTriggeredEventMB = 0; //init
 


  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
     ->IsEventSelected() & ftrigBit1 ){
    fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
     ->IsEventSelected() & ftrigBit2 ){
    // From AliVEvent:
    //    kINT7         = BIT(1), // V0AND trigger, offline V0 selection
    fTriggeredEventMB += 2;  
  }

  
  // Get process type for MC
  fMcProcessType = 0; // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD

  // real data that are not triggered we skip
  if(!fAnalysisMC && !fTriggeredEventMB)
    return; 
  
  if (fAnalysisMC) {
    
    if (fAnalysisType == "ESD"){

      AliHeader* headerMC = fMC->Header();
      if (headerMC) {
	
	AliGenEventHeader* genHeader = headerMC->GenEventHeader();
	TArrayF vtxMC(3); // primary vertex  MC 
	vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
	if (genHeader) {
	  genHeader->PrimaryVertex(vtxMC);
	}
	fZvtxMC = vtxMC[2];
	
	// PYTHIA:
	AliGenPythiaEventHeader* pythiaGenHeader =
	  dynamic_cast<AliGenPythiaEventHeader*>(headerMC->GenEventHeader());
	if (pythiaGenHeader) {  //works only for pythia
	  fMcProcessType =  GetPythiaEventProcessType(pythiaGenHeader->ProcessType());
	}
	// PHOJET:
	AliGenDPMjetEventHeader* dpmJetGenHeader =
	  dynamic_cast<AliGenDPMjetEventHeader*>(headerMC->GenEventHeader());
	if (dpmJetGenHeader) {
	  fMcProcessType = GetDPMjetEventProcessType(dpmJetGenHeader->ProcessType());
	}
      }
    } else { // AOD
      
      AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader")); 
      if(mcHeader) {
	fZvtxMC = mcHeader->GetVtxZ();
	
	if(strstr(mcHeader->GetGeneratorName(), "Pythia")) {
	  fMcProcessType =  GetPythiaEventProcessType(mcHeader->GetEventType());
	} else {
	  fMcProcessType =  GetDPMjetEventProcessType(mcHeader->GetEventType());
	}
      }
    }
  }
  
  // There are 3 cases
  // Vertex: NO  - status -1
  // Vertex: YES : outside cut - status 0
  //             : inside cut  - status 1
  // We have to be careful how we normalize because we probably want to
  // normalize to:
  // Nevents=(No vertex + outside + inside)/(out + in)*in
  
  if (fAnalysisType == "ESD")
    fZvtx = GetVertex(fESD);
  else // AOD
    fZvtx = GetVertex(fAOD);
    
  fVtxStatus = -999;
  
  if(fZvtx<-990) {
    
    fVtxStatus = -1;
    if(fTriggeredEventMB)
      fVtx->Fill(0);
    if(fAnalysisMC)
      fVtxMC->Fill(0);
  } else {
    
    if(fTriggeredEventMB)
      fVtx->Fill(1);
    if(fAnalysisMC)
      fVtxMC->Fill(1);
    fVtxBeforeCuts->Fill(fZvtx);
    fVtxStatus = 0;
    if (TMath::Abs(fZvtx) < fVtxCut) {	
      fVtxAfterCuts->Fill(fZvtx);
      fVtxStatus = 1;
    }
  }
  
  // store MC event data no matter what
  if(fAnalysisMC && fStoreMcIn) {

    if (fAnalysisType == "ESD") {
      ProcessMCTruthESD();
    } else { // AOD
      ProcessMCTruthAOD();
    }
  }      
  Float_t centrality = -10;
  // only analyze triggered events
  if(fTriggeredEventMB) {
    
    if (fAnalysisType == "ESD"){
      if(fAnalysisPbPb){
	AliCentrality *centObject = fESD->GetCentrality();
	centrality = centObject->GetCentralityPercentile("V0M"); 
	if((centrality>fMaxCent)||(centrality<fMinCent))return;
      }
      Int_t nv0test=fESD->GetNumberOfV0s();
      cout<<"&&&&&&&&&&&&&&&&   hello world"<<endl;
      fnv0->Fill(nv0test);
      AnalyzeESD(fESD);
      
    } else { // AOD
      if(fAnalysisPbPb){
	AliCentrality *centObject = fAOD->GetCentrality();
	centrality = centObject->GetCentralityPercentile("V0M"); 
	if((centrality>fMaxCent)||(centrality<fMinCent))return;
      }
      Int_t nv0test=fAOD->GetNumberOfV0s();
      fnv0->Fill(nv0test);
      AnalyzeAOD(fAOD);
      
    }
  }
  


  if( fTreeOption) {

    fEvent->process = fMcProcessType;
    fEvent->trig    = fTriggeredEventMB;
    fEvent->zvtxMC  = fZvtxMC;
    fEvent->cent      = centrality;

    if(!fRequireRecV0 || fEvent->n > 0) // only fill if there are accepted V0s
      fTree->Fill();
    fV0ArrayGlobalPar->Clear();
    fV0ArrayTPCPar->Clear();
    if (fAnalysisMC && fStoreMcIn)
      fTrackArrayMC->Clear();
  }
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::AnalyzeESD(AliESDEvent* esdEvent)
{
  fRun  = esdEvent->GetRunNumber();
  fEventId = 0;
  if(esdEvent->GetHeader())
    fEventId = GetEventIdAsLong(esdEvent->GetHeader());

  //  Int_t     event     = esdEvent->GetEventNumberInFile();
  UInt_t    time      = esdEvent->GetTimeStamp();
  //  ULong64_t trigger   = esdEvent->GetTriggerMask();
  Float_t   magf      = esdEvent->GetMagneticField();



  Short_t isPileup = esdEvent->IsPileupFromSPD();

  // centrality
  Float_t centrality = 120;
  AliCentrality *centObject = esdEvent->GetCentrality();
  if(centObject)
    centrality = centObject->GetCentralityPercentile("V0M"); 

  if(fTriggeredEventMB) {// Only MC case can we have not triggered events
    
    if(fVtxStatus==1) { // accepted vertex

      // accepted event
      fEvents->Fill(0);
      
      ProduceArrayTrksESD( esdEvent, kGlobalTrk );//produce array with global track parameters
      ProduceArrayTrksESD( esdEvent, kTPCTrk );//array with tpc track parametes


    } // end if vtx status
  } // end if triggered
  
  if(fTreeOption) {

    fEvent->run       = fRun;
    fEvent->eventid   = fEventId;
    fEvent->time      = time;
    //fEvent->cent      = centrality;
    fEvent->mag       = magf;
    fEvent->zvtx      = fZvtx;
    fEvent->vtxstatus = fVtxStatus;
    //fEvent->trackmult = trackmult;
    //fEvent->n         = nadded;
    fEvent->pileup    = isPileup;
  }
}

//________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::AnalyzeAOD(AliAODEvent* aodEvent)
{
  fRun  = aodEvent->GetRunNumber();
  fEventId = 0;
  if(aodEvent->GetHeader())
    fEventId = GetEventIdAsLong(aodEvent->GetHeader());

  //  Int_t     event     = aodEvent->GetEventNumberInFile();
  UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();
  //  ULong64_t trigger   = aodEvent->GetTriggerMask();
  Float_t   magf      = aodEvent->GetMagneticField();

  Int_t     trackmult = 0; // no pt cuts
  Int_t     nadded    = 0;

  Short_t isPileup = aodEvent->IsPileupFromSPD();

  // centrality
  Float_t centrality = 120;
  AliCentrality *centObject = aodEvent->GetCentrality();
  if(centObject)
    centrality = centObject->GetCentralityPercentile("V0M"); 

  if(fTriggeredEventMB) {// Only MC case can we have not triggered events
    
    if(fVtxStatus==1) { // accepted vertex

      // accepted event
      fEvents->Fill(0);
      
      ProduceArrayTrksAOD( aodEvent, kGlobalTrk );//produce array with global track parameters
      ProduceArrayTrksAOD( aodEvent, kTPCTrk );//array with tpc track parametes


    } // end if vtx status
  } // end if triggered
  
  if(fTreeOption) {

    fEvent->run       = fRun;
    fEvent->eventid   = fEventId;
    fEvent->time      = time;
    //fEvent->cent      = centrality;
    fEvent->mag       = magf;
    fEvent->zvtx      = fZvtx;
    fEvent->vtxstatus = fVtxStatus;
    //fEvent->trackmult = trackmult;
    //fEvent->n         = nadded;
    fEvent->pileup    = isPileup;
  }







}

//_____________________________________________________________________________
Float_t AliAnalysisTaskHighPtDeDxV0::GetVertex(const AliVEvent* event) const
{
  Float_t zvtx = -999;
  
  const AliVVertex* primaryVertex = event->GetPrimaryVertex(); 
  
  if(primaryVertex->GetNContributors()>0)
    zvtx = primaryVertex->GetZ();

  return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskHighPtDeDxV0::GetPidCode(Int_t pdgCode) const 
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
    pidCode = 5; // proton
    break;
  default:
    pidCode = 6;  // something else?
  };
  
  return pidCode;
}


//_____________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::ProcessMCTruthESD() 
{
  // Fill the special MC histogram with the MC truth info
  
  Short_t trackmult = 0;
  Short_t nadded    = 0;
  const Int_t nTracksMC = fMCStack->GetNtrack();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    //Cuts
    if(!(fMCStack->IsPhysicalPrimary(iTracks)))
      continue;
    
    TParticle* trackMC = fMCStack->Particle(iTracks);
    
    Double_t chargeMC = trackMC->GetPDG()->Charge();
    if (chargeMC != 0) // select charge 0 particles!
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;
    
    trackmult++;
    
    //   "charge:pt:p:eta:phi:pidCode"
    Float_t ptMC      = trackMC->Pt();
    Float_t pMC       = trackMC->P();
    Float_t etaMC     = trackMC->Eta();
    Float_t phiMC     = trackMC->Phi();
    
    Int_t pdgCode = trackMC->GetPdgCode();
    if(!(pdgCode==310 || pdgCode == 3122 || pdgCode == -3122))
      continue;
    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(pdgCode);
    
    // Warning: In the data code we cut on the daughters.....
    if (trackMC->Pt() < fMinPt)
      continue;
    
    if(fTreeOption) {
      
      DeDxTrackMC* track = new((*fTrackArrayMC)[nadded]) DeDxTrackMC();
      nadded++;
      
      track->pMC   = pMC;
      track->ptMC  = ptMC;
      track->etaMC = etaMC;
      track->phiMC = phiMC;
      track->qMC   = Short_t(chargeMC);
      track->pidMC = pidCodeMC;
      track->pdgMC = pdgCode;
    }
    
  }//MC track loop
  
  if(fTreeOption) {
    
    Sort(fTrackArrayMC, kTRUE);

    fEvent->trackmultMC = trackmult;
    fEvent->nMC         = nadded;
  }
  
}

//_____________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::ProcessMCTruthAOD() 
{
  // Fill the special MC histogram with the MC truth info

  Short_t trackmult = 0;
  Short_t nadded    = 0;
  const Int_t nTracksMC = fMCArray->GetEntriesFast();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
    
    //Cuts
    if(!(trackMC->IsPhysicalPrimary()))
      continue;
    
    
    Double_t chargeMC = trackMC->Charge();
    if (chargeMC != 0) // Keep the neutral particles
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;
    
    trackmult++;
    
    //   "charge:pt:p:eta:phi:pidCode"
    Float_t ptMC      = trackMC->Pt();
    Float_t pMC       = trackMC->P();
    Float_t etaMC     = trackMC->Eta();
    Float_t phiMC     = trackMC->Phi();
    
    Int_t pdgCode = trackMC->PdgCode();
    if(!(pdgCode==310 || pdgCode == 3122 || pdgCode == -3122))
      continue;
    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(pdgCode);
    
    if (trackMC->Pt() < fMinPt)
      continue;
    
    if(fTreeOption) {
      
      DeDxTrackMC* track = new((*fTrackArrayMC)[nadded]) DeDxTrackMC();
      nadded++;
      
      track->pMC   = pMC;
      track->ptMC  = ptMC;
      track->etaMC = etaMC;
      track->phiMC = phiMC;
      track->qMC   = Short_t(chargeMC);
      track->pidMC = pidCodeMC;
      track->pdgMC = pdgCode;
    }
    
  }//MC track loop
  
  if(fTreeOption) {
    
    Sort(fTrackArrayMC, kTRUE);

    fEvent->trackmultMC = trackmult;
    fEvent->nMC         = nadded;
  }
  
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskHighPtDeDxV0::GetPythiaEventProcessType(Int_t pythiaType) {
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
Short_t AliAnalysisTaskHighPtDeDxV0::GetDPMjetEventProcessType(Int_t dpmJetType) {
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
ULong64_t AliAnalysisTaskHighPtDeDxV0::GetEventIdAsLong(AliVHeader* header) const
{
  // To have a unique id for each event in a run!
  // Modified from AliRawReader.h
  return ((ULong64_t)header->GetBunchCrossNumber()+
	  (ULong64_t)header->GetOrbitNumber()*3564+
	  (ULong64_t)header->GetPeriodNumber()*16777215*3564);
}


//____________________________________________________________________
TParticle* AliAnalysisTaskHighPtDeDxV0::FindPrimaryMother(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  Int_t nSteps = 0;

  Int_t motherLabel = FindPrimaryMotherLabel(stack, label, nSteps);
  if (motherLabel < 0)
    return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskHighPtDeDxV0::FindPrimaryMotherLabel(AliStack* stack, Int_t label, Int_t& nSteps)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // returns its label
  //
  // Taken from AliPWG0Helper class
  //
  nSteps = 0;
  const Int_t nPrim  = stack->GetNprimary();
  
  while (label >= nPrim) {

    //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

    nSteps++; // 1 level down
    
    TParticle* particle = stack->Particle(label);
    if (!particle) {
      
      AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
      return -1;
    }
    
    // find mother
    if (particle->GetMother(0) < 0) {

      AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
      return -1;
    }
    
    label = particle->GetMother(0);
  }
  
  return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskHighPtDeDxV0::FindPrimaryMotherAOD(AliAODMCParticle* startParticle, Int_t& nSteps)
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
void AliAnalysisTaskHighPtDeDxV0::Sort(TClonesArray* array, Bool_t isMC) 
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
//__________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::ProduceArrayTrksESD( AliESDEvent *ESDevent, AnalysisMode analysisMode ){
  Int_t nv0s = ESDevent->GetNumberOfV0s();
  if(nv0s<1)return;
  Int_t     trackmult = 0; // no pt cuts
  Int_t     nadded    = 0;
  const AliESDVertex *myBestPrimaryVertex = ESDevent->GetPrimaryVertex();
  if (!myBestPrimaryVertex) return;
  if (!(myBestPrimaryVertex->GetStatus())) return;
  Double_t  lPrimaryVtxPosition[3];
  myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
  
  Double_t  lPrimaryVtxCov[6];
  myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
  Double_t  lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
  
  AliAODVertex* myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);


  if( analysisMode == kGlobalTrk ){
    if(fV0ArrayGlobalPar)
      fV0ArrayGlobalPar->Clear();
  } else if( analysisMode == kTPCTrk ){
    if(fV0ArrayTPCPar)
      fV0ArrayTPCPar->Clear();
  }
  
  if( analysisMode==kGlobalTrk ){  



    
    // ################################
    // #### BEGINNING OF V0 CODE ######
    // ################################

    
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
      
      // This is the begining of the V0 loop  
      AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
      if (!esdV0) continue;
      
      // AliESDTrack (V0 Daughters)
      UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());
      
      AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
      AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
	Printf("ERROR: Could not retreive one of the daughter track");
	continue;
      }
      
      // Remove like-sign
      if (pTrack->GetSign() == nTrack->GetSign()) {
	//cout<< "like sign, continue"<< endl;
	continue;
      } 
      
      // Eta cut on decay products
      if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
	continue;
      
      // Pt cut on decay products
      if (esdV0->Pt() < fMinPt)
	//	if (pTrack->Pt() < fMinPt && nTrack->Pt() < fMinPt)
	continue;

      //filter for positive track
      UShort_t filterFlag_p = 0;
      Bool_t filterCut_Set1_p = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2_p = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3_p = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
      
      UInt_t selectDebug_p = 0;
      if (fTrackFilterGolden) {
	selectDebug_p = fTrackFilterGolden->IsSelected(pTrack);
	if (selectDebug_p) {
	  filterFlag_p +=1;
	  filterCut_Set3_p=kTRUE;
	}
      }
      
      if (fTrackFilterTPC) {
	
	selectDebug_p = fTrackFilterTPC->IsSelected(pTrack);
	if (selectDebug_p){//only tracks which pass the TPC-only track cuts
	  filterFlag_p +=2;
	  filterCut_Set1_p=kTRUE;
	  
	}
	
      }
      
      if (fTrackFilter) {
	selectDebug_p = fTrackFilter->IsSelected(pTrack);
	if (selectDebug_p) {
	  filterCut_Set2_p=kTRUE;
	}
      }
      
      if(filterFlag_p ==0 )
	continue;
      

      //filter for negative track
      UShort_t filterFlag_n = 0;
      Bool_t filterCut_Set1_n = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2_n = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3_n = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
      
      UInt_t selectDebug_n = 0;
      if (fTrackFilterGolden) {
	selectDebug_n = fTrackFilterGolden->IsSelected(nTrack);
	if (selectDebug_n) {
	  filterFlag_n +=1;
	  filterCut_Set3_n=kTRUE;
	}
      }
      
      if (fTrackFilterTPC) {
	
	selectDebug_n = fTrackFilterTPC->IsSelected(nTrack);
	if (selectDebug_n){//only tracks which pass the TPC-only track cuts
	  filterFlag_n +=2;
	  filterCut_Set1_n=kTRUE;
	  
	}
	
      }
      
      if (fTrackFilter) {
	selectDebug_n = fTrackFilter->IsSelected(nTrack);
	if (selectDebug_n) {
	  filterCut_Set2_n=kTRUE;
	}
      }


      // Check if switch does anything!
      Bool_t isSwitched = kFALSE;
      if (pTrack->GetSign() < 0) { // switch
	
	isSwitched = kTRUE;
	AliESDtrack* helpTrack = nTrack;
	nTrack = pTrack;
	pTrack = helpTrack;
      }	
      
      Float_t alpha = esdV0->AlphaV0();
      Float_t ptarm = esdV0->PtArmV0();
      // Double_t pVtxPos= v0->PrimaryVtxPosition();      
      
      Double_t  lV0Position[3];
      esdV0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
      
      Double_t lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
      Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
					    TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
					    TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));
      AliKFVertex primaryVtxKF( *myPrimaryVertex );
      AliKFParticle::SetField(ESDevent->GetMagneticField());
      
      // Also implement switch here!!!!!!
      AliKFParticle* negEKF  = 0; // e-
      AliKFParticle* posEKF  = 0; // e+
      AliKFParticle* negPiKF = 0; // pi -
      AliKFParticle* posPiKF = 0; // pi +
      AliKFParticle* posPKF  = 0; // p
      AliKFParticle* negAPKF = 0; // p-bar
      
      if(!isSwitched) {
	negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
	posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
	negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
	posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
	posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
	negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
      } else { // switch + and - 
	negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
	posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
	negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
	posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
	posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
	negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
      }
      
      AliKFParticle v0GKF;  // Gamma e.g. from pi0
      v0GKF+=(*negEKF);
      v0GKF+=(*posEKF);
      v0GKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0K0sKF; // K0 short
      v0K0sKF+=(*negPiKF);
      v0K0sKF+=(*posPiKF);
      v0K0sKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0LambdaKF; // Lambda
      v0LambdaKF+=(*negPiKF);
      v0LambdaKF+=(*posPKF);	
      v0LambdaKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0AntiLambdaKF; // Lambda-bar
      v0AntiLambdaKF+=(*posPiKF);
      v0AntiLambdaKF+=(*negAPKF);
      v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);
      
      Double_t deltaInvMassG     = v0GKF.GetMass();
      Double_t deltaInvMassK0s   = v0K0sKF.GetMass()-0.498;
      Double_t deltaInvMassL     = v0LambdaKF.GetMass()-1.116;
      Double_t deltaInvMassAntiL = v0AntiLambdaKF.GetMass()-1.116;
      
      if(TMath::Abs(deltaInvMassK0s) > fMassCut &&
	 TMath::Abs(deltaInvMassL) > fMassCut &&
	 TMath::Abs(deltaInvMassAntiL) > fMassCut)
	continue;
      
      // Extract track information
      
      Short_t pcharge  = pTrack->Charge();
      Float_t ppt      = pTrack->Pt();
      Float_t pp       = pTrack->P(); 
      Float_t peta     = pTrack->Eta();
      Float_t pphi     = pTrack->Phi();
      Short_t pncl     = pTrack->GetTPCsignalN();
      Short_t pneff    = Short_t(pTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t pdedx    = pTrack->GetTPCsignal();
      
      Float_t ptpcchi  = 0;
      if(pTrack->GetTPCNcls() > 0)
	ptpcchi = pTrack->GetTPCchi2()/Float_t(pTrack->GetTPCNcls());
      
      Short_t ncharge  = nTrack->Charge();
      Float_t npt      = nTrack->Pt();
      Float_t np       = nTrack->P(); 
      Float_t neta     = nTrack->Eta();
      Float_t nphi     = nTrack->Phi();
      Short_t nncl     = nTrack->GetTPCsignalN();
      Short_t nneff    = Short_t(nTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t ndedx    = nTrack->GetTPCsignal();
      Float_t ntpcchi  = 0;
      if(nTrack->GetTPCNcls() > 0)
	ntpcchi = nTrack->GetTPCchi2()/Float_t(nTrack->GetTPCNcls());
      
      Float_t b[2];
      Float_t bCov[3];
      pTrack->GetImpactParameters(b,bCov);
      Float_t pdcaxy   = b[0];
      Float_t pdcaz    = b[1];
      nTrack->GetImpactParameters(b,bCov);
      Float_t ndcaxy   = b[0];
      Float_t ndcaz    = b[1];
      
      Int_t   primaryV0     = 0; // 0 means that the tracks are not both daughters of a primary particle (1 means they are)
      Int_t   pdgV0         = 0; // 0 means that they don't have same origin for MC (1 means they have the same original mother)
      Float_t p_ptMC        = 0;
      Short_t p_pidCode     = 0; // 0 = real data / no mc track!
      Short_t p_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   p_pdgMother   = 0;
      Float_t n_ptMC        = 0;
      Short_t n_pidCode     = 0; // 0 = real data / no mc track!
      Short_t n_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   n_pdgMother   = 0;
      if(fAnalysisMC) {
	
	Int_t p_mother_label = 0;
	Int_t p_mother_steps = 0;
	Int_t n_mother_label = 0;
	Int_t n_mother_steps = 0;
	
	// positive track
	const Int_t p_label = TMath::Abs(pTrack->GetLabel());
	TParticle* p_mcTrack = fMCStack->Particle(p_label);	    
	if (p_mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(p_label))
	    p_primaryFlag = 1;
	  
	  Int_t p_pdgCode = p_mcTrack->GetPdgCode();
	  p_pidCode = GetPidCode(p_pdgCode);
	  
	  p_ptMC      = p_mcTrack->Pt();
	  
	  p_mother_label = FindPrimaryMotherLabel(fMCStack, p_label, 
						  p_mother_steps);
	  if(p_mother_label>0) {
	    TParticle* p_mother = fMCStack->Particle(p_mother_label);
	    p_pdgMother = p_mother->GetPdgCode();
	  }
	}
	
	// negative track
	const Int_t n_label = TMath::Abs(nTrack->GetLabel());
	TParticle* n_mcTrack = fMCStack->Particle(n_label);	    
	if (n_mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(n_label))
	    n_primaryFlag = 1;
	  
	  Int_t n_pdgCode = n_mcTrack->GetPdgCode();
	  n_pidCode = GetPidCode(n_pdgCode);
	  
	  n_ptMC      = n_mcTrack->Pt();
	  
	  n_mother_label = FindPrimaryMotherLabel(fMCStack, n_label, 
						  n_mother_steps);
	  if(n_mother_label>0) {
	    TParticle* n_mother = fMCStack->Particle(n_mother_label);
	    n_pdgMother = n_mother->GetPdgCode();
	  }
	}
	
	// Check if V0 is primary = first and the same mother of both partciles
	if(p_mother_label>0 && n_mother_label>0 && p_mother_label == n_mother_label) {
	  pdgV0 = p_pdgMother;
	  if(p_mother_steps == 1 && n_mother_steps == 1) {
	    primaryV0 = 1;
	  }
	}
      }

 
    
      if(fTreeOption) {
	


	DeDxV0* v0data = new((*fV0ArrayGlobalPar)[nadded]) DeDxV0();
	nadded++;
	
	// v0 data
	v0data->p       = esdV0->P();
	v0data->pt      = esdV0->Pt();
	v0data->eta     = esdV0->Eta();
	v0data->phi     = esdV0->Phi();
	v0data->pdca    = TMath::Sqrt(pdcaxy*pdcaxy + pdcaz*pdcaz);
	v0data->ndca    = TMath::Sqrt(ndcaxy*ndcaxy + ndcaz*ndcaz);
	v0data->dmassG  = deltaInvMassG;
	v0data->dmassK0 = deltaInvMassK0s;
	v0data->dmassL  = deltaInvMassL;
	v0data->dmassAL = deltaInvMassAntiL;
	v0data->alpha   = alpha;
	v0data->ptarm   = ptarm;
	v0data->decayr  = lV0Radius;
	v0data->decayl  = lV0DecayLength;
	
	// New parameters
	v0data->status  = esdV0->GetOnFlyStatus();
	v0data->chi2    = esdV0->GetChi2V0();
	v0data->cospt   = esdV0->GetV0CosineOfPointingAngle(); 
	// cospt: as I understand this means that the pointing to the vertex
	// is fine so I remove the dcaxy and dcaz for the V= class
	v0data->dcadaughters = esdV0->GetDcaV0Daughters();
	v0data->primary = primaryV0;
	v0data->pdg     = pdgV0;
	
	// positive track
	v0data->ptrack.p       = pp;
	v0data->ptrack.pt      = ppt;
	//	  v0data->ptrack.ptcon   = ppt_con;
	//	  v0data->ptrack.tpcchi  = ptpcchi;
	v0data->ptrack.eta     = peta;
	v0data->ptrack.phi     = pphi;
	v0data->ptrack.q       = pcharge;
	v0data->ptrack.ncl     = pncl;
	v0data->ptrack.neff    = pneff;
	v0data->ptrack.dedx    = pdedx;
	v0data->ptrack.dcaxy   = pdcaxy;
	v0data->ptrack.dcaz    = pdcaz;
	v0data->ptrack.pid     = p_pidCode;
	v0data->ptrack.primary = p_primaryFlag;
	v0data->ptrack.pttrue  = p_ptMC;
	v0data->ptrack.mother  = p_pdgMother;
	v0data->ptrack.filter  = filterFlag_p;
	v0data->ptrack.filterset1 = filterCut_Set1_p;
	v0data->ptrack.filterset2 = filterCut_Set2_p;
	v0data->ptrack.filterset3 = filterCut_Set3_p;
	
	// negative track
	v0data->ntrack.p       = np;
	v0data->ntrack.pt      = npt;
	//	  v0data->ntrack.ptcon   = npt_con;
	//	  v0data->ntrack.tpcchi  = ntpcchi;
	v0data->ntrack.eta     = neta;
	v0data->ntrack.phi     = nphi;
	v0data->ntrack.q       = ncharge;
	v0data->ntrack.ncl     = nncl;
	v0data->ntrack.neff    = nneff;
	v0data->ntrack.dedx    = ndedx;
	v0data->ntrack.dcaxy   = ndcaxy;
	v0data->ntrack.dcaz    = ndcaz;
	v0data->ntrack.pid     = n_pidCode;
	v0data->ntrack.primary = n_primaryFlag;
	v0data->ntrack.pttrue  = n_ptMC;
	v0data->ntrack.mother  = n_pdgMother;
	v0data->ntrack.filter  = filterFlag_n;
	v0data->ntrack.filterset1 = filterCut_Set1_n;
	v0data->ntrack.filterset2 = filterCut_Set2_n;
	v0data->ntrack.filterset3 = filterCut_Set3_n;


      }
      
      // clean up loop over v0

      delete negPiKF;
      delete posPiKF;
      delete posPKF;
      delete negAPKF;
  


    }
  
    // clean up event
    //delete myPrimaryVertex;
  }
  else if( analysisMode==kTPCTrk ){  
    cout<<"&&&&&&&&&&&&&&&&&&&&&                Hello world"<<endl;
    const AliESDVertex *vtxSPD = ESDevent->GetPrimaryVertexSPD();
    if( vtxSPD->GetNContributors() < 1 || TMath::Abs(vtxSPD->GetZ()) > 10.0 ) return;
    
    
    // ################################
    // #### BEGINNING OF V0 CODE ######
    // ################################

    
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
      
      // This is the begining of the V0 loop  
      AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
      if (!esdV0) continue;
      
      // AliESDTrack (V0 Daughters)
      UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());
      
      AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
      AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
	Printf("ERROR: Could not retreive one of the daughter track");
	continue;
      }

      // Remove like-sign
      if (pTrack->GetSign() == nTrack->GetSign()) {
	//cout<< "like sign, continue"<< endl;
	continue;
      } 

      AliESDtrack *pTrackTPC = AliESDtrackCuts::GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(ESDevent),pTrack->GetID());
      AliESDtrack *nTrackTPC = AliESDtrackCuts::GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(ESDevent),nTrack->GetID());

      if (!pTrackTPC || !nTrackTPC) {
	Printf("ERROR: Could not retreive one of the daughter TPC track");
	continue;
      }

      //filter for positive track
      UInt_t selectDebug_p = 0;
      UShort_t filterFlag_p = 0;
      selectDebug_p = fTrackFilterTPC->IsSelected(pTrackTPC);
      //if(selectDebug_p==0) continue;

      //filter for negative track
      UInt_t selectDebug_n = 0;
      UShort_t filterFlag_n = 0;
      selectDebug_n = fTrackFilterTPC->IsSelected(nTrackTPC);
      if(selectDebug_n==0 || selectDebug_p==0) continue;

      if(selectDebug_p)
	filterFlag_p += 1;

      if(pTrackTPC->Pt()>0.){
	// only constrain tracks above threshold
	AliExternalTrackParam exParamp;
	// take the B-field from the ESD, no 3D fieldMap available at this point
	Bool_t relate_p = false;
	relate_p = pTrackTPC->RelateToVertexTPC(vtxSPD,ESDevent->GetMagneticField(),
						kVeryBig,&exParamp);
	if(!relate_p){
	  delete pTrackTPC;
	  continue;
	}
	pTrackTPC->Set(exParamp.GetX(),exParamp.GetAlpha(),exParamp.GetParameter(),
		       exParamp.GetCovariance());
      }
      else continue;     
      

      //filter for negative track
 
      if(selectDebug_n)
	filterFlag_n += 1;

  
      if(nTrackTPC->Pt()>0.){
	// only constrain tracks above threshold
	AliExternalTrackParam exParamn;
	// take the B-field from the ESD, no 3D fieldMap available at this point
	Bool_t relate_n = false;
	relate_n = nTrackTPC->RelateToVertexTPC(vtxSPD,ESDevent->GetMagneticField(),
						kVeryBig,&exParamn);
	if(!relate_n){
	  delete nTrackTPC;
	  continue;
	}
	nTrackTPC->Set(exParamn.GetX(),exParamn.GetAlpha(),exParamn.GetParameter(),
		       exParamn.GetCovariance());
      }
      else continue;  
      
      
      // Eta cut on decay products
      if(TMath::Abs(pTrackTPC->Eta()) > fEtaCut || TMath::Abs(nTrackTPC->Eta()) > fEtaCut)
	continue;
      
      // Pt cut on decay products
      if (esdV0->Pt() < fMinPt)
	//	if (pTrack->Pt() < fMinPt && nTrack->Pt() < fMinPt)
	continue;
 
      
      // Check if switch does anything!
      Bool_t isSwitched = kFALSE;
      if (pTrackTPC->GetSign() < 0) { // switch
	
	isSwitched = kTRUE;
	AliESDtrack* helpTrack = nTrack;
	nTrackTPC = pTrackTPC;
	pTrackTPC = helpTrack;
      }	
      

      // Extract track information
      
      Short_t pcharge  = pTrackTPC->Charge();
      Float_t ppt      = pTrackTPC->Pt();
      Float_t pp       = pTrackTPC->P(); 
      Float_t peta     = pTrackTPC->Eta();
      Float_t pphi     = pTrackTPC->Phi();
      Short_t pncl     = pTrackTPC->GetTPCsignalN();
      Short_t pneff    = Short_t(pTrackTPC->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t pdedx    = pTrackTPC->GetTPCsignal();
      
      Float_t ptpcchi  = 0;
      if(pTrackTPC->GetTPCNcls() > 0)
	ptpcchi = pTrackTPC->GetTPCchi2()/Float_t(pTrackTPC->GetTPCNcls());
      
      Short_t ncharge  = nTrackTPC->Charge();
      Float_t npt      = nTrackTPC->Pt();
      Float_t np       = nTrackTPC->P(); 
      Float_t neta     = nTrackTPC->Eta();
      Float_t nphi     = nTrackTPC->Phi();
      Short_t nncl     = nTrackTPC->GetTPCsignalN();
      Short_t nneff    = Short_t(nTrackTPC->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t ndedx    = nTrackTPC->GetTPCsignal();
      Float_t ntpcchi  = 0;
      if(nTrackTPC->GetTPCNcls() > 0)
	ntpcchi = nTrackTPC->GetTPCchi2()/Float_t(nTrackTPC->GetTPCNcls());
      
      Float_t bp[2]={0,0};
      Float_t bCovp[3]={0,0,0};
      pTrackTPC->GetImpactParameters(bp,bCovp);
      Float_t pdcaxy   = bp[0];
      Float_t pdcaz    = bp[1];

      Float_t bn[2]={0,0};
      Float_t bCovn[3]={0,0,0};
      nTrackTPC->GetImpactParameters(bn,bCovn);
      Float_t ndcaxy   = bn[0];
      Float_t ndcaz    = bn[1];


      Float_t alpha = esdV0->AlphaV0();
      Float_t ptarm = esdV0->PtArmV0();
      // Double_t pVtxPos= v0->PrimaryVtxPosition();      
      
      Double_t  lV0Position[3];
      esdV0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
      
      Double_t lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
      Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
					    TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
					    TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));
      AliKFVertex primaryVtxKF( *myPrimaryVertex );
      AliKFParticle::SetField(ESDevent->GetMagneticField());
      
      // Also implement switch here!!!!!!
      AliKFParticle* negEKF  = 0; // e-
      AliKFParticle* posEKF  = 0; // e+
      AliKFParticle* negPiKF = 0; // pi -
      AliKFParticle* posPiKF = 0; // pi +
      AliKFParticle* posPKF  = 0; // p
      AliKFParticle* negAPKF = 0; // p-bar
      
      
      if(!isSwitched) {
	/*	
		negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
		posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
		negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
		posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
		posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
		negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
	*/

	negEKF  = new AliKFParticle( *(dynamic_cast<AliVTrack*>(nTrackTPC)) , 11);
	posEKF  = new AliKFParticle( *(dynamic_cast<AliVTrack*>(pTrackTPC)) ,-11);
	negPiKF = new AliKFParticle( *(dynamic_cast<AliVTrack*>(nTrackTPC)) ,-211);
	posPiKF = new AliKFParticle( *(dynamic_cast<AliVTrack*>(pTrackTPC)) , 211);
	posPKF  = new AliKFParticle( *(dynamic_cast<AliVTrack*>(pTrackTPC)) , 2212);
	negAPKF = new AliKFParticle( *(dynamic_cast<AliVTrack*>(nTrackTPC)) ,-2212);



      } else { // switch + and - 
	/*
	negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
	posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
	negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
	posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
	posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
	negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
	*/

	negEKF  = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(pTrackTPC)), 11);
	posEKF  = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(nTrackTPC)),-11);
	negPiKF = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(pTrackTPC)),-211);
	posPiKF = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(nTrackTPC)), 211);
	posPKF  = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(nTrackTPC)), 2212);
	negAPKF = new AliKFParticle(  *(dynamic_cast<AliVTrack*>(pTrackTPC)),-2212);


      }
   



      
      AliKFParticle v0GKF;  // Gamma e.g. from pi0
      v0GKF+=(*negEKF);
      v0GKF+=(*posEKF);
      v0GKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0K0sKF; // K0 short
      v0K0sKF+=(*negPiKF);
      v0K0sKF+=(*posPiKF);
      v0K0sKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0LambdaKF; // Lambda
      v0LambdaKF+=(*negPiKF);
      v0LambdaKF+=(*posPKF);	
      v0LambdaKF.SetProductionVertex(primaryVtxKF);
      
      AliKFParticle v0AntiLambdaKF; // Lambda-bar
      v0AntiLambdaKF+=(*posPiKF);
      v0AntiLambdaKF+=(*negAPKF);
      v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);
      
      Double_t deltaInvMassG     = v0GKF.GetMass();
      Double_t deltaInvMassK0s   = v0K0sKF.GetMass()-0.498;
      Double_t deltaInvMassL     = v0LambdaKF.GetMass()-1.116;
      Double_t deltaInvMassAntiL = v0AntiLambdaKF.GetMass()-1.116;
      
      if(TMath::Abs(deltaInvMassK0s) > fMassCut &&
	 TMath::Abs(deltaInvMassL) > fMassCut &&
	 TMath::Abs(deltaInvMassAntiL) > fMassCut)
	continue;
      
      // TODO: Whe should these be different? Different mass hypothesis = energy loss
      // This is not important for us as we focus on the decay products!
      // Double_t ptK0s        = v0K0sKF.GetPt(); 
      // Double_t ptL          = v0LambdaKF.GetPt();
      // Double_t ptAntiL      = v0AntiLambdaKF.GetPt();     
      
      Int_t   primaryV0     = 0; // 0 means that the tracks are not both daughters of a primary particle (1 means they are)
      Int_t   pdgV0         = 0; // 0 means that they don't have same origin for MC (1 means they have the same original mother)
      Float_t p_ptMC        = 0;
      Short_t p_pidCode     = 0; // 0 = real data / no mc track!
      Short_t p_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   p_pdgMother   = 0;
      Float_t n_ptMC        = 0;
      Short_t n_pidCode     = 0; // 0 = real data / no mc track!
      Short_t n_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   n_pdgMother   = 0;
      if(fAnalysisMC) {
	
	Int_t p_mother_label = 0;
	Int_t p_mother_steps = 0;
	Int_t n_mother_label = 0;
	Int_t n_mother_steps = 0;
	
	// positive track
	const Int_t p_label = TMath::Abs(pTrackTPC->GetLabel());
	TParticle* p_mcTrack = fMCStack->Particle(p_label);	    
	if (p_mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(p_label))
	    p_primaryFlag = 1;
	  
	  Int_t p_pdgCode = p_mcTrack->GetPdgCode();
	  p_pidCode = GetPidCode(p_pdgCode);
	  
	  p_ptMC      = p_mcTrack->Pt();
	  
	  p_mother_label = FindPrimaryMotherLabel(fMCStack, p_label, 
						  p_mother_steps);
	  if(p_mother_label>0) {
	    TParticle* p_mother = fMCStack->Particle(p_mother_label);
	    p_pdgMother = p_mother->GetPdgCode();
	  }
	}
	
	// negative track
	const Int_t n_label = TMath::Abs(nTrackTPC->GetLabel());
	TParticle* n_mcTrack = fMCStack->Particle(n_label);	    
	if (n_mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(n_label))
	    n_primaryFlag = 1;
	  
	  Int_t n_pdgCode = n_mcTrack->GetPdgCode();
	  n_pidCode = GetPidCode(n_pdgCode);
	  
	  n_ptMC      = n_mcTrack->Pt();
	  
	  n_mother_label = FindPrimaryMotherLabel(fMCStack, n_label, 
						  n_mother_steps);
	  if(n_mother_label>0) {
	    TParticle* n_mother = fMCStack->Particle(n_mother_label);
	    n_pdgMother = n_mother->GetPdgCode();
	  }
	}
	
	// Check if V0 is primary = first and the same mother of both partciles
	if(p_mother_label>0 && n_mother_label>0 && p_mother_label == n_mother_label) {
	  pdgV0 = p_pdgMother;
	  if(p_mother_steps == 1 && n_mother_steps == 1) {
	    primaryV0 = 1;
	  }
	}
      }
      

      if(fTreeOption) {
	
	DeDxV0* v0datatpc = new((*fV0ArrayTPCPar)[nadded]) DeDxV0();
	nadded++;
	
	// v0 data
	v0datatpc->p       = esdV0->P();
	v0datatpc->pt      = esdV0->Pt();
	v0datatpc->eta     = esdV0->Eta();
	v0datatpc->phi     = esdV0->Phi();
	v0datatpc->pdca    = TMath::Sqrt(pdcaxy*pdcaxy + pdcaz*pdcaz);
	v0datatpc->ndca    = TMath::Sqrt(ndcaxy*ndcaxy + ndcaz*ndcaz);
	v0datatpc->dmassG  = deltaInvMassG;
	v0datatpc->dmassK0 = deltaInvMassK0s;
	v0datatpc->dmassL  = deltaInvMassL;
	v0datatpc->dmassAL = deltaInvMassAntiL;
	v0datatpc->alpha   = alpha;
	v0datatpc->ptarm   = ptarm;
	v0datatpc->decayr  = lV0Radius;
	v0datatpc->decayl  = lV0DecayLength;
	
	// New parameters
	v0datatpc->status  = esdV0->GetOnFlyStatus();
	v0datatpc->chi2    = esdV0->GetChi2V0();
	v0datatpc->cospt   = esdV0->GetV0CosineOfPointingAngle(); 
	// cospt: as I understand this means that the pointing to the vertex
	// is fine so I remove the dcaxy and dcaz for the V= class
	v0datatpc->dcadaughters = esdV0->GetDcaV0Daughters();
	v0datatpc->primary = primaryV0;
	v0datatpc->pdg     = pdgV0;
	
	// positive track
	v0datatpc->ptrack.p       = pp;
	v0datatpc->ptrack.pt      = ppt;
	//	  v0data->ptrack.ptcon   = ppt_con;
	//	  v0data->ptrack.tpcchi  = ptpcchi;
	v0datatpc->ptrack.eta     = peta;
	v0datatpc->ptrack.phi     = pphi;
	v0datatpc->ptrack.q       = pcharge;
	v0datatpc->ptrack.ncl     = pncl;
	v0datatpc->ptrack.neff    = pneff;
	v0datatpc->ptrack.dedx    = pdedx;
	v0datatpc->ptrack.dcaxy   = pdcaxy;
	v0datatpc->ptrack.dcaz    = pdcaz;
	v0datatpc->ptrack.pid     = p_pidCode;
	v0datatpc->ptrack.primary = p_primaryFlag;
	v0datatpc->ptrack.pttrue  = p_ptMC;
	v0datatpc->ptrack.mother  = p_pdgMother;
	v0datatpc->ptrack.filter  = filterFlag_p;
	v0datatpc->ptrack.filterset1 = 0;
	v0datatpc->ptrack.filterset2 = 0;
	v0datatpc->ptrack.filterset3 = 0;
	
	// negative track
	v0datatpc->ntrack.p       = np;
	v0datatpc->ntrack.pt      = npt;
	//	  v0data->ntrack.ptcon   = npt_con;
	//	  v0data->ntrack.tpcchi  = ntpcchi;
	v0datatpc->ntrack.eta     = neta;
	v0datatpc->ntrack.phi     = nphi;
	v0datatpc->ntrack.q       = ncharge;
	v0datatpc->ntrack.ncl     = nncl;
	v0datatpc->ntrack.neff    = nneff;
	v0datatpc->ntrack.dedx    = ndedx;
	v0datatpc->ntrack.dcaxy   = ndcaxy;
	v0datatpc->ntrack.dcaz    = ndcaz;
	v0datatpc->ntrack.pid     = n_pidCode;
	v0datatpc->ntrack.primary = n_primaryFlag;
	v0datatpc->ntrack.pttrue  = n_ptMC;
	v0datatpc->ntrack.mother  = n_pdgMother;
	v0datatpc->ntrack.filter  = filterFlag_n;
	v0datatpc->ntrack.filterset1 = 0;
	v0datatpc->ntrack.filterset2 = 0;
	v0datatpc->ntrack.filterset3 = 0;
      }
      
      // clean up loop over v0
      
      delete negPiKF;
      delete posPiKF;
      delete posPKF;
      delete negAPKF;
 


    }

 
 


  }
  delete myPrimaryVertex;

  if(fTreeOption) {
    
    if( analysisMode==kGlobalTrk ){

      
      fEvent->trackmult = trackmult;
      fEvent->n         = nadded;


    }
 
    
  }


  
}
//_______________________________________________________________________________________________________________
void AliAnalysisTaskHighPtDeDxV0::ProduceArrayTrksAOD( AliAODEvent *AODevent, AnalysisMode analysisMode ){
  Int_t nv0s = AODevent->GetNumberOfV0s();
  if(nv0s<1)return;
  Int_t     trackmult = 0; // no pt cuts
  Int_t     nadded    = 0;
  
  AliAODVertex *myBestPrimaryVertex = AODevent->GetPrimaryVertex();
  if (!myBestPrimaryVertex) return;
  
  
  if( analysisMode == kGlobalTrk ){
    if(fV0ArrayGlobalPar)
      fV0ArrayGlobalPar->Clear();
  } else if( analysisMode == kTPCTrk ){
    if(fV0ArrayTPCPar)
      fV0ArrayTPCPar->Clear();
  }
  
  if( analysisMode==kGlobalTrk ){  
    
    
    // ################################
    // #### BEGINNING OF V0 CODE ######
    // ################################
    // This is the begining of the V0 loop  
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
      AliAODv0 *aodV0 = AODevent->GetV0(iV0);
      if (!aodV0) continue;
      
      // common part
      
      // AliAODTrack (V0 Daughters)
      AliAODVertex* vertex = aodV0->GetSecondaryVtx();
      if (!vertex) {
	Printf("ERROR: Could not retrieve vertex");
	continue;
      }
      
      AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0);
      AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
      if (!pTrack || !nTrack) {
	Printf("ERROR: Could not retrieve one of the daughter track");
	continue;
      }
      
      // Remove like-sign
      if (pTrack->Charge() == nTrack->Charge()) {
	//cout<< "like sign, continue"<< endl;
	continue;
      } 
      
      // Make sure charge ordering is ok
      if (pTrack->Charge() < 0) {
	AliAODTrack* helpTrack = pTrack;
	pTrack = nTrack;
	nTrack = helpTrack;
      } 
      
      // Eta cut on decay products
      if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
	continue;
      
      // Pt cut on decay products
      if (aodV0->Pt() < fMinPt)
	//	if (pTrack->Pt() < fMinPt && nTrack->Pt() < fMinPt)
	continue;
      
      //check positive tracks
      UShort_t filterFlag_p = 0;
      Bool_t filterCut_Set1_p = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2_p = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3_p = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
      
      if (fTrackFilterGolden) {  
	// ITSTPC2010 cuts is bit 32 according to above macro, new versions of aliroot includes the golden cuts
	
	if(pTrack->TestFilterBit(32)) {
	  filterFlag_p +=1;
	  filterCut_Set3_p = kTRUE;
	}
      }
      
      
      if (fTrackFilterTPC) {
	// TPC only cuts is bit 1 according to above macro
	// Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
	if(pTrack->TestFilterBit(1)){
	  filterFlag_p +=2;
	  filterCut_Set1_p = kTRUE;
	  
	}
      }
      
      if(filterFlag_p==0)
	continue;
      
      //check negative tracks
      UShort_t filterFlag_n = 0;
      Bool_t filterCut_Set1_n = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2_n = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3_n = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
      
      
      if (fTrackFilterGolden) {
	
	// ITSTPC2010 cuts is bit 32 according to above macro, new versions of aliroot includes the golden cuts
	if(nTrack->TestFilterBit(32)) {
	  filterFlag_n +=1;
	  filterCut_Set3_n = kTRUE;
	}
      }
      
      
      if (fTrackFilterTPC) {
	// TPC only cuts is bit 1 according to above macro
	// Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
	if(nTrack->TestFilterBit(1)){
	  filterFlag_n +=2;
	  filterCut_Set1_n = kTRUE;
	  
	}
      }
      
      if(filterFlag_n==0)
	continue;
      
      
      Float_t alpha = aodV0->AlphaV0();
      Float_t ptarm = aodV0->PtArmV0();
      // Double_t pVtxPos= v0->PrimaryVtxPosition();      
      
      Double_t lV0Radius      = aodV0->RadiusV0();
      Double_t lV0DecayLength = aodV0->DecayLength(myBestPrimaryVertex);
      
      Double_t deltaInvMassG     = aodV0->InvMass2Prongs(0,1,11,11);
      Double_t deltaInvMassK0s   = aodV0->MassK0Short()-0.498;
      Double_t deltaInvMassL     = aodV0->MassLambda()-1.116;
      Double_t deltaInvMassAntiL = aodV0->MassAntiLambda()-1.116;
      
      if(TMath::Abs(deltaInvMassK0s) > fMassCut &&
	 TMath::Abs(deltaInvMassL) > fMassCut &&
	 TMath::Abs(deltaInvMassAntiL) > fMassCut)
	continue;
      
      // TODO: Why should these be different? Different mass hypothesis = energy loss
      // This is not important for us as we focus on the decay products!
      // Double_t ptK0s        = v0K0sKF.GetPt(); 
      // Double_t ptL          = v0LambdaKF.GetPt();
      // Double_t ptAntiL      = v0AntiLambdaKF.GetPt();     
      
      // Extract track information
      
      Double_t b[2], cov[3];
      if(!pTrack->PropagateToDCA(myBestPrimaryVertex, AODevent->GetMagneticField(), kVeryBig, b, cov))
	filterFlag_p += 32; // propagation failed!!!!!
      
      Float_t pdcaxy   = b[0];
      Float_t pdcaz    = b[1];
      if(!nTrack->PropagateToDCA(myBestPrimaryVertex, AODevent->GetMagneticField(), kVeryBig, b, cov))
	filterFlag_n += 32; // propagation failed!!!!!
      Float_t ndcaxy   = b[0];
      Float_t ndcaz    = b[1];
      
      Short_t pcharge  = pTrack->Charge();
      Float_t ppt      = pTrack->Pt();
      Float_t pp       = pTrack->P(); 
      Float_t peta     = pTrack->Eta();
      Float_t pphi     = pTrack->Phi();
      //	Float_t ptpcchi  = pTrack->Chi2perNDF();
      
      AliAODPid* pPid = pTrack->GetDetPid();
      Short_t pncl     = -10;
      Short_t pneff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t pdedx    = -10;
      Float_t pbeta = -99;
      if(pPid) {
	pncl     = pPid->GetTPCsignalN();
	pdedx    = pPid->GetTPCsignal();
	//TOF
	if (pTrack->GetStatus()&AliESDtrack::kTOFpid){
	  Double_t tof[5];
	  pPid->GetIntegratedTimes(tof);
	  pbeta = tof[0]/pPid->GetTOFsignal();
	}
      }
      
      Short_t ncharge  = nTrack->Charge();
      Float_t npt      = nTrack->Pt();
      Float_t np       = nTrack->P(); 
      Float_t neta     = nTrack->Eta();
      Float_t nphi     = nTrack->Phi();
      //	Float_t ntpcchi  = nTrack->Chi2perNDF();
      
      AliAODPid* nPid = nTrack->GetDetPid();
      Short_t nncl     = -10;
      Short_t nneff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      Float_t ndedx    = -10;
      Float_t nbeta = -99;
      if(pPid) {
	nncl     = nPid->GetTPCsignalN();
	ndedx    = nPid->GetTPCsignal();
	//TOF
	if (nTrack->GetStatus()&AliESDtrack::kTOFpid){
	  Double_t tof[5];
	  nPid->GetIntegratedTimes(tof);
	  nbeta = tof[0]/nPid->GetTOFsignal();
	}
      }
      
      Int_t   primaryV0     = 0; // 0 means that the tracks are not both daughters of a primary particle (1 means they are)
      Int_t   pdgV0         = 0; // 0 means that they don't have same origin for MC (1 means they have the same original mother)
      Float_t p_ptMC        = 0;
      Short_t p_pidCode     = 0; // 0 = real data / no mc track!
      Short_t p_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   p_pdgMother   = 0;
      Float_t n_ptMC        = 0;
      Short_t n_pidCode     = 0; // 0 = real data / no mc track!
      Short_t n_primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   n_pdgMother   = 0;
      if(fAnalysisMC) {
	
	AliAODMCParticle* p_mother = 0;
	Int_t p_mother_steps = 0;
	AliAODMCParticle* n_mother = 0;
	Int_t n_mother_steps = 0;
	
	// positive track
	const Int_t p_label = TMath::Abs(pTrack->GetLabel());
	
	AliAODMCParticle* p_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(p_label));
	if (p_mcTrack){
	  
	  if(p_mcTrack->IsPhysicalPrimary())
	    p_primaryFlag = 1;
	  
	  Int_t p_pdgCode = p_mcTrack->GetPdgCode();
	  p_pidCode = GetPidCode(p_pdgCode);
	  
	  p_ptMC      = p_mcTrack->Pt();
	  
	  p_mother = FindPrimaryMotherAOD(p_mcTrack, p_mother_steps);
	  if(p_mother)
	    p_pdgMother = p_mother->GetPdgCode();
	}
	
	// negative track
	const Int_t n_label = TMath::Abs(pTrack->GetLabel());
	
	AliAODMCParticle* n_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(n_label));
	if (n_mcTrack){
	  
	  if(n_mcTrack->IsPhysicalPrimary())
	    n_primaryFlag = 1;
	  
	  Int_t n_pdgCode = n_mcTrack->GetPdgCode();
	  n_pidCode = GetPidCode(n_pdgCode);
	  
	  n_ptMC      = n_mcTrack->Pt();
	  
	  n_mother = FindPrimaryMotherAOD(n_mcTrack, n_mother_steps);
	  if(n_mother)
	    n_pdgMother = n_mother->GetPdgCode();
	}
	
	// Check if V0 is primary = first and the same mother of both partciles
	if(p_mother && n_mother && p_mother == n_mother) {
	  pdgV0 = p_pdgMother;
	  if(p_mother_steps == 1 && n_mother_steps == 1) {
	    primaryV0 = 1;
	  }
	}
      }
      
      if(fTreeOption) {
	
	DeDxV0* v0data = new((*fV0ArrayGlobalPar)[nadded]) DeDxV0();
	nadded++;
	
	// v0 data
	v0data->p       = aodV0->P();
	v0data->pt      = aodV0->Pt();
	v0data->eta     = aodV0->Eta();
	v0data->phi     = aodV0->Phi();
	v0data->pdca    = aodV0->DcaPosToPrimVertex();
	v0data->ndca    = aodV0->DcaNegToPrimVertex();
	v0data->dmassG  = deltaInvMassG;
	v0data->dmassK0 = deltaInvMassK0s;
	v0data->dmassL  = deltaInvMassL;
	v0data->dmassAL = deltaInvMassAntiL;
	v0data->alpha   = alpha;
	v0data->ptarm   = ptarm;
	v0data->decayr  = lV0Radius;
	v0data->decayl  = lV0DecayLength;
	// v0data->pdca    = TMath::Sqrt(pdcaxy*pdcaxy + pdcaz*pdcaz);
	// v0data->ndca    = TMath::Sqrt(ndcaxy*ndcaxy + ndcaz*ndcaz);
	
	// New parameters
	v0data->status  = aodV0->GetOnFlyStatus();
	v0data->chi2    = aodV0->Chi2V0();
	v0data->cospt   = aodV0->CosPointingAngle(myBestPrimaryVertex);
	// cospt: as I understand this means that the pointing to the vertex
	// is fine so I remove the dcaxy and dcaz for the V= class
	v0data->dcav0   = aodV0->DcaV0ToPrimVertex();
	v0data->dcadaughters = aodV0->DcaV0Daughters();
	v0data->primary = primaryV0;
	v0data->pdg     = pdgV0;
	
	// positive track
	v0data->ptrack.p       = pp;
	v0data->ptrack.pt      = ppt;
	//	  v0data->ptrack.ptcon   = ppt_con;
	//	  v0data->ptrack.tpcchi  = ptpcchi;
	v0data->ptrack.eta     = peta;
	v0data->ptrack.phi     = pphi;
	v0data->ptrack.q       = pcharge;
	v0data->ptrack.ncl     = pncl;
	v0data->ptrack.neff    = pneff;
	v0data->ptrack.dedx    = pdedx;
	v0data->ptrack.dcaxy   = pdcaxy;
	v0data->ptrack.dcaz    = pdcaz;
	v0data->ptrack.pid     = p_pidCode;
	v0data->ptrack.primary = p_primaryFlag;
	v0data->ptrack.pttrue  = p_ptMC;
	v0data->ptrack.mother  = p_pdgMother;
	v0data->ptrack.filter  = filterFlag_p;
	v0data->ptrack.filterset1 = filterCut_Set1_p;
	v0data->ptrack.filterset2 = filterCut_Set2_p;
	v0data->ptrack.filterset3 = filterCut_Set3_p;
	
	
	// negative track
	v0data->ntrack.p       = np;
	v0data->ntrack.pt      = npt;
	//	  v0data->ntrack.ptcon   = npt_con;
	//	  v0data->ntrack.tpcchi  = ntpcchi;
	v0data->ntrack.eta     = neta;
	v0data->ntrack.phi     = nphi;
	v0data->ntrack.q       = ncharge;
	v0data->ntrack.ncl     = nncl;
	v0data->ntrack.neff    = nneff;
	v0data->ntrack.dedx    = ndedx;
	v0data->ntrack.dcaxy   = ndcaxy;
	v0data->ntrack.dcaz    = ndcaz;
	v0data->ntrack.pid     = n_pidCode;
	v0data->ntrack.primary = n_primaryFlag;
	v0data->ntrack.pttrue  = n_ptMC;
	v0data->ntrack.mother  = n_pdgMother;
	v0data->ntrack.filter  = filterFlag_n;
	v0data->ntrack.filterset1 = filterCut_Set1_n;
	v0data->ntrack.filterset2 = filterCut_Set2_n;
	v0data->ntrack.filterset3 = filterCut_Set3_n;
	
	
      }
    }//end loop over v0's
    
    
  }else if( analysisMode==kTPCTrk ){  
    
    const AliAODVertex*	vertexSPD= (AliAODVertex*)AODevent->GetPrimaryVertexSPD();//GetPrimaryVertex()
    cout<<"&&&&&&&&&&&&&&&&&&&&&                Hello world  0"<<endl;
    if( vertexSPD->GetNContributors() < 1 || TMath::Abs(vertexSPD->GetZ()) > 10.0 ) return;
    
    cout<<"&&&&&&&&&&&&&&&&&&&&&                Hello world"<<endl;
    // ################################
    // #### BEGINNING OF V0 CODE ######
    // ################################
    // This is the begining of the V0 loop  



     for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
       AliAODv0 *aodV0 = AODevent->GetV0(iV0);
       if (!aodV0) continue;
       
       // common part
       
       // AliAODTrack (V0 Daughters)
       AliAODVertex* vertex = aodV0->GetSecondaryVtx();
       if (!vertex) {
	 Printf("ERROR: Could not retrieve vertex");
	 continue;
       }
       
       AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0);
       AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
       if (!pTrack || !nTrack) {
	 Printf("ERROR: Could not retrieve one of the daughter track");
	 continue;
       }
       
       
       // Remove like-sign
       if (pTrack->Charge() == nTrack->Charge()) {
	 //cout<< "like sign, continue"<< endl;
	 continue;
       } 
       
       // Make sure charge ordering is ok
       if (pTrack->Charge() < 0) {
	 AliAODTrack* helpTrack = pTrack;
	 pTrack = nTrack;
	 nTrack = helpTrack;
       } 
       
       // Eta cut on decay products
       if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
	 continue;

       cout<<"Eta positive track:"<<pTrack->Eta()<<endl;

       
       // Pt cut on decay products
       if (aodV0->Pt() < fMinPt)
	 //	if (pTrack->Pt() < fMinPt && nTrack->Pt() < fMinPt)
	 continue;
       
       //check positive tracks
       UShort_t filterFlag_p = 0;
       Bool_t filterCut_Set1_p = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
       Bool_t filterCut_Set2_p = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
       Bool_t filterCut_Set3_p = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
       
       // TPC only cuts is bit 1 according to above macro
       // Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
       if(pTrack->TestFilterBit(128)) {
	 cout<<"este track paso el corte bit 128"<<endl;
	 filterFlag_p +=1;
       }
       cout<<"filterFlag_p="<<filterFlag_p<<endl;      




       if(filterFlag_p==0)
	 continue;
 


       //check negative tracks
       UShort_t filterFlag_n = 0;
       Bool_t filterCut_Set1_n = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
       Bool_t filterCut_Set2_n = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
       Bool_t filterCut_Set3_n = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
       
       // TPC only cuts is bit 1 according to above macro
       // Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
       if(nTrack->TestFilterBit(128)) {
	 filterFlag_n +=1;
       }
       
       if(filterFlag_n==0)
	 continue;
       
       
       Float_t alpha = aodV0->AlphaV0();
       Float_t ptarm = aodV0->PtArmV0();
       // Double_t pVtxPos= v0->PrimaryVtxPosition();      
       
       Double_t lV0Radius      = aodV0->RadiusV0();
       Double_t lV0DecayLength = aodV0->DecayLength(myBestPrimaryVertex);
       
       Double_t deltaInvMassG     = aodV0->InvMass2Prongs(0,1,11,11);
       Double_t deltaInvMassK0s   = aodV0->MassK0Short()-0.498;
       Double_t deltaInvMassL     = aodV0->MassLambda()-1.116;
       Double_t deltaInvMassAntiL = aodV0->MassAntiLambda()-1.116;
       
       if(TMath::Abs(deltaInvMassK0s) > fMassCut &&
	  TMath::Abs(deltaInvMassL) > fMassCut &&
	  TMath::Abs(deltaInvMassAntiL) > fMassCut)
	 continue;
       
       // TODO: Why should these be different? Different mass hypothesis = energy loss
       // This is not important for us as we focus on the decay products!
       // Double_t ptK0s        = v0K0sKF.GetPt(); 
       // Double_t ptL          = v0LambdaKF.GetPt();
       // Double_t ptAntiL      = v0AntiLambdaKF.GetPt();     
       
       // Extract track information
       
       Double_t b[2], cov[3];
       if(!pTrack->PropagateToDCA(vertexSPD, AODevent->GetMagneticField(), kVeryBig, b, cov))
	 filterFlag_p += 32; // propagation failed!!!!!
       
       Float_t pdcaxy   = b[0];
       Float_t pdcaz    = b[1];
       if(!nTrack->PropagateToDCA(vertexSPD, AODevent->GetMagneticField(), kVeryBig, b, cov))
	 filterFlag_n += 32; // propagation failed!!!!!
       Float_t ndcaxy   = b[0];
       Float_t ndcaz    = b[1];
       
       Short_t pcharge  = pTrack->Charge();
       Float_t ppt      = pTrack->Pt();
       Float_t pp       = pTrack->P(); 
       Float_t peta     = pTrack->Eta();
       Float_t pphi     = pTrack->Phi();
       //	Float_t ptpcchi  = pTrack->Chi2perNDF();
       
       AliAODPid* pPid = pTrack->GetDetPid();
       Short_t pncl     = -10;
       Short_t pneff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
       Float_t pdedx    = -10;
       Float_t pbeta = -99;
       if(pPid) {
	 pncl     = pPid->GetTPCsignalN();
	 pdedx    = pPid->GetTPCsignal();
	 //TOF
	 if (pTrack->GetStatus()&AliESDtrack::kTOFpid){
	   Double_t tof[5];
	   pPid->GetIntegratedTimes(tof);
	   pbeta = tof[0]/pPid->GetTOFsignal();
	 }
       }
       
       Short_t ncharge  = nTrack->Charge();
       Float_t npt      = nTrack->Pt();
       Float_t np       = nTrack->P(); 
       Float_t neta     = nTrack->Eta();
       Float_t nphi     = nTrack->Phi();
       //	Float_t ntpcchi  = nTrack->Chi2perNDF();
       
       AliAODPid* nPid = nTrack->GetDetPid();
       Short_t nncl     = -10;
       Short_t nneff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
       Float_t ndedx    = -10;
       Float_t nbeta = -99;
       if(pPid) {
	 nncl     = nPid->GetTPCsignalN();
	 ndedx    = nPid->GetTPCsignal();
	 //TOF
	 if (nTrack->GetStatus()&AliESDtrack::kTOFpid){
	   Double_t tof[5];
	   nPid->GetIntegratedTimes(tof);
	   nbeta = tof[0]/nPid->GetTOFsignal();
	 }
       }
       
       Int_t   primaryV0     = 0; // 0 means that the tracks are not both daughters of a primary particle (1 means they are)
       Int_t   pdgV0         = 0; // 0 means that they don't have same origin for MC (1 means they have the same original mother)
       Float_t p_ptMC        = 0;
       Short_t p_pidCode     = 0; // 0 = real data / no mc track!
       Short_t p_primaryFlag = 0; // 0 = real data / not primary mc track  
       Int_t   p_pdgMother   = 0;
       Float_t n_ptMC        = 0;
       Short_t n_pidCode     = 0; // 0 = real data / no mc track!
       Short_t n_primaryFlag = 0; // 0 = real data / not primary mc track  
       Int_t   n_pdgMother   = 0;
       if(fAnalysisMC) {
	 
	 AliAODMCParticle* p_mother = 0;
	 Int_t p_mother_steps = 0;
	 AliAODMCParticle* n_mother = 0;
	 Int_t n_mother_steps = 0;
	 
	 // positive track
	 const Int_t p_label = TMath::Abs(pTrack->GetLabel());
	 
	 AliAODMCParticle* p_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(p_label));
	 if (p_mcTrack){
	   
	   if(p_mcTrack->IsPhysicalPrimary())
	     p_primaryFlag = 1;
	   
	   Int_t p_pdgCode = p_mcTrack->GetPdgCode();
	   p_pidCode = GetPidCode(p_pdgCode);
	   
	   p_ptMC      = p_mcTrack->Pt();
	   
	   p_mother = FindPrimaryMotherAOD(p_mcTrack, p_mother_steps);
	   if(p_mother)
	     p_pdgMother = p_mother->GetPdgCode();
	 }
	 
	 // negative track
	 const Int_t n_label = TMath::Abs(pTrack->GetLabel());
	 
	 AliAODMCParticle* n_mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(n_label));
	 if (n_mcTrack){
	   
	   if(n_mcTrack->IsPhysicalPrimary())
	     n_primaryFlag = 1;
	   
	   Int_t n_pdgCode = n_mcTrack->GetPdgCode();
	   n_pidCode = GetPidCode(n_pdgCode);
	   
	   n_ptMC      = n_mcTrack->Pt();
	   
	   n_mother = FindPrimaryMotherAOD(n_mcTrack, n_mother_steps);
	   if(n_mother)
	     n_pdgMother = n_mother->GetPdgCode();
	 }
	 
	 // Check if V0 is primary = first and the same mother of both partciles
	 if(p_mother && n_mother && p_mother == n_mother) {
	   pdgV0 = p_pdgMother;
	   if(p_mother_steps == 1 && n_mother_steps == 1) {
	     primaryV0 = 1;
	   }
	 }
       }
       
       if(fTreeOption) {
	 
	 DeDxV0* v0datatpc = new((*fV0ArrayTPCPar)[nadded]) DeDxV0();
	 nadded++;
	 
	 // v0 data
	 v0datatpc->p       = aodV0->P();
	 v0datatpc->pt      = aodV0->Pt();
	 v0datatpc->eta     = aodV0->Eta();
	 v0datatpc->phi     = aodV0->Phi();
	 v0datatpc->pdca    = aodV0->DcaPosToPrimVertex();
	 v0datatpc->ndca    = aodV0->DcaNegToPrimVertex();
	 v0datatpc->dmassG  = deltaInvMassG;
	 v0datatpc->dmassK0 = deltaInvMassK0s;
	 v0datatpc->dmassL  = deltaInvMassL;
	 v0datatpc->dmassAL = deltaInvMassAntiL;
	 v0datatpc->alpha   = alpha;
	 v0datatpc->ptarm   = ptarm;
	 v0datatpc->decayr  = lV0Radius;
	 v0datatpc->decayl  = lV0DecayLength;
	 // v0data->pdca    = TMath::Sqrt(pdcaxy*pdcaxy + pdcaz*pdcaz);
	 // v0data->ndca    = TMath::Sqrt(ndcaxy*ndcaxy + ndcaz*ndcaz);
	 
	 // New parameters
	 v0datatpc->status  = aodV0->GetOnFlyStatus();
	 v0datatpc->chi2    = aodV0->Chi2V0();
	 v0datatpc->cospt   = aodV0->CosPointingAngle(myBestPrimaryVertex);
	 // cospt: as I understand this means that the pointing to the vertex
	 // is fine so I remove the dcaxy and dcaz for the V= class
	 v0datatpc->dcav0   = aodV0->DcaV0ToPrimVertex();
	 v0datatpc->dcadaughters = aodV0->DcaV0Daughters();
	 v0datatpc->primary = primaryV0;
	 v0datatpc->pdg     = pdgV0;
	 
	 // positive track
	 v0datatpc->ptrack.p       = pp;
	 v0datatpc->ptrack.pt      = ppt;
	 //	  v0data->ptrack.ptcon   = ppt_con;
	 //	  v0data->ptrack.tpcchi  = ptpcchi;
	 v0datatpc->ptrack.eta     = peta;
	 v0datatpc->ptrack.phi     = pphi;
	 v0datatpc->ptrack.q       = pcharge;
	 v0datatpc->ptrack.ncl     = pncl;
	 v0datatpc->ptrack.neff    = pneff;
	 v0datatpc->ptrack.dedx    = pdedx;
	 v0datatpc->ptrack.dcaxy   = pdcaxy;
	 v0datatpc->ptrack.dcaz    = pdcaz;
	 v0datatpc->ptrack.pid     = p_pidCode;
	 v0datatpc->ptrack.primary = p_primaryFlag;
	 v0datatpc->ptrack.pttrue  = p_ptMC;
	 v0datatpc->ptrack.mother  = p_pdgMother;
	 v0datatpc->ptrack.filter  = filterFlag_p;
	 v0datatpc->ptrack.filterset1 = 0;
	 v0datatpc->ptrack.filterset2 = 0;
	 v0datatpc->ptrack.filterset3 = 0;
	 
	 
	 // negative track
	 v0datatpc->ntrack.p       = np;
	 v0datatpc->ntrack.pt      = npt;
	 //	  v0data->ntrack.ptcon   = npt_con;
	 //	  v0data->ntrack.tpcchi  = ntpcchi;
	 v0datatpc->ntrack.eta     = neta;
	 v0datatpc->ntrack.phi     = nphi;
	 v0datatpc->ntrack.q       = ncharge;
	 v0datatpc->ntrack.ncl     = nncl;
	 v0datatpc->ntrack.neff    = nneff;
	 v0datatpc->ntrack.dedx    = ndedx;
	 v0datatpc->ntrack.dcaxy   = ndcaxy;
	 v0datatpc->ntrack.dcaz    = ndcaz;
	 v0datatpc->ntrack.pid     = n_pidCode;
	 v0datatpc->ntrack.primary = n_primaryFlag;
	 v0datatpc->ntrack.pttrue  = n_ptMC;
	 v0datatpc->ntrack.mother  = n_pdgMother;
	 v0datatpc->ntrack.filter  = filterFlag_n;
	 v0datatpc->ntrack.filterset1 = 0;
	 v0datatpc->ntrack.filterset2 = 0;
	 v0datatpc->ntrack.filterset3 = 0;
	 
	 
       }
     }//end v0's loop
    
  }
  
  if(fTreeOption) {
    
    if( analysisMode==kGlobalTrk ){

      
      fEvent->trackmult = trackmult;
      fEvent->n         = nadded;
      
      
    }
    
    
  }
  
  
  
}



