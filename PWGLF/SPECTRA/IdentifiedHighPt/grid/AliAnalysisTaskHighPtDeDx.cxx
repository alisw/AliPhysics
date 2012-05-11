#include "AliAnalysisTaskHighPtDeDx.h"

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
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include "AliCentrality.h" 

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 


// STL includes
#include <iostream>
using namespace std;


//
// Responsible:
// Alexandru Dobrin (Wayne State) 
// Peter Christiansen (Lund)
//

/* 
To do:

Separate the code into two

*/




ClassImp(AliAnalysisTaskHighPtDeDx)

const Double_t AliAnalysisTaskHighPtDeDx::fgkClight = 2.99792458e-2;

//_____________________________________________________________________________
AliAnalysisTaskHighPtDeDx::AliAnalysisTaskHighPtDeDx():
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
  fTPCBranch(kFALSE),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  fRandom(0x0),
  fEvent(0x0),
  fTrackArrayGlobalPar(0x0),
  fTrackArrayTPCPar(0x0),
  fTrackArrayMC(0x0),
  fVZEROArray(0x0),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinPt(0.1),
  fMinCent(0.0),
  fMaxCent(100.0),
  fLowPtFraction(0.01),
  fTreeOption(0),
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
  fn1(0),
  fn2(0),
  fcent(0)


{
  // Default constructor (should not be used)

  //  fRandom = new TRandom(0); // 0 means random seed
}

//______________________________________________________________________________
AliAnalysisTaskHighPtDeDx::AliAnalysisTaskHighPtDeDx(const char *name):
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
  fTPCBranch(kFALSE),
  ftrigBit1(0x0),
  ftrigBit2(0x0),
  fRandom(0x0),
  fEvent(0x0),
  fTrackArrayGlobalPar(0x0),
  fTrackArrayMC(0x0),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinPt(0.1),
  fLowPtFraction(0.01),
  fTreeOption(0),
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
  fn1(0),
  fn2(0),
  fcent(0)

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
  
  // //for proof running; I cannot create tree do to memory limitations -> work with THnSparse 
  // if (fListOfObjects  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
  
  
  
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

  fn1=new TH1F("fn1","fn1",5001,-1,5000);
  fListOfObjects->Add(fn1);

  fn2=new TH1F("fn2","fn2",5001,-1,5000);
  fListOfObjects->Add(fn2);

  fcent=new TH1F("fcent","fcent",104,-2,102);
  fListOfObjects->Add(fcent);

  fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
  fListOfObjects->Add(fVtx);

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

    fTrackArrayGlobalPar = new TClonesArray("DeDxTrack", 1000);
    fTree->Bronch("trackGlobalPar", "TClonesArray", &fTrackArrayGlobalPar);
    if(fTPCBranch){
      fTrackArrayTPCPar = new TClonesArray("DeDxTrack", 1000);
      fTree->Bronch("trackTPCPar", "TClonesArray", &fTrackArrayTPCPar);
    }

    fVZEROArray = new TClonesArray("VZEROCell", 1000);
    fTree->Bronch("cellVZERO", "TClonesArray", &fVZEROArray);

    if (fAnalysisMC) {    
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


  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);
    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }    
  } else {
    fAOD = dynamic_cast<AliAODEvent*>(event);
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
    fn1->Fill(1);
    fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
     ->IsEventSelected() & ftrigBit2 ){
    // From AliVEvent:
    //    kINT7         = BIT(1), // V0AND trigger, offline V0 selection
    fTriggeredEventMB += 2;  
    fn2->Fill(1);
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
  if(fAnalysisMC) {

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
      fcent->Fill(centrality);
      AnalyzeESD(fESD);
    } else { // AOD
      if(fAnalysisPbPb){
	AliCentrality *centObject = fAOD->GetCentrality();
	if(centObject)
	  centrality = centObject->GetCentralityPercentile("V0M"); 
	if((centrality>fMaxCent)||(centrality<fMinCent))return;
      }
      fcent->Fill(centrality);
      AnalyzeAOD(fAOD);
    }
  }

  if( fTreeOption) {
    
    fEvent->process = fMcProcessType;
    fEvent->trig    = fTriggeredEventMB;
    fEvent->zvtxMC  = fZvtxMC;
    fEvent->cent      = centrality;

    fTree->Fill();
    fTrackArrayGlobalPar->Clear();
    if(fTPCBranch)
      fTrackArrayTPCPar->Clear();
    fVZEROArray->Clear();

    if (fAnalysisMC)    
      fTrackArrayMC->Clear();
  }
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskHighPtDeDx::AnalyzeESD(AliESDEvent* esdEvent)
{
  fRun  = esdEvent->GetRunNumber();
  fEventId = 0;
  if(esdEvent->GetHeader())
    fEventId = GetEventIdAsLong(esdEvent->GetHeader());
  
  Short_t isPileup = esdEvent->IsPileupFromSPD();
  
  //  Int_t     event     = esdEvent->GetEventNumberInFile();
  UInt_t    time      = esdEvent->GetTimeStamp();
  //  ULong64_t trigger   = esdEvent->GetTriggerMask();
  Float_t   magf      = esdEvent->GetMagneticField();





  if(fTriggeredEventMB) {// Only MC case can we have not triggered events
    
    // accepted event
    fEvents->Fill(0);
    
    
    if(fVtxStatus!=1) return; // accepted vertex
    Int_t nESDTracks = esdEvent->GetNumberOfTracks();
    
    if(nESDTracks<1)return;
    cout<<"Multiplicity="<<nESDTracks<<endl;
    
    ProduceArrayTrksESD( esdEvent, kGlobalTrk );//produce array with global track parameters
    if(fTPCBranch)
      ProduceArrayTrksESD( esdEvent, kTPCTrk );//array with tpc track parametes

    fEvents->Fill(1);
    AliESDVZERO *esdV0 = esdEvent->GetVZEROData();// loop sobre canales del V0 para obtener las multiplicidad

    for (Int_t iCh=0; iCh<64; ++iCh) { 
      Float_t multv=esdV0->GetMultiplicity(iCh); 
      Int_t intexv=iCh;
      VZEROCell* cellv0 = new((*fVZEROArray)[iCh]) VZEROCell();
      cellv0->cellindex=intexv;
      cellv0->cellmult= multv;
    }   



  } // end if triggered
  
  if(fTreeOption) {

    fEvent->run       = fRun;
    fEvent->eventid   = fEventId;
    fEvent->time      = time;
    //fEvent->cent      = centrality;
    fEvent->mag       = magf;
    fEvent->zvtx      = fZvtx;
    fEvent->vtxstatus = fVtxStatus;
    fEvent->pileup    = isPileup;

  }




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

  //Int_t     trackmult = 0; // no pt cuts
  //Int_t     nadded    = 0;

  Short_t isPileup = aodEvent->IsPileupFromSPD();




  if(fTriggeredEventMB) {// Only MC case can we have not triggered events
    
    // accepted event
    fEvents->Fill(0);
    
    if(fVtxStatus!=1) return; // accepted vertex
    Int_t nAODTracks = aodEvent->GetNumberOfTracks();
    if(nAODTracks<1)return;      

    ProduceArrayTrksAOD( aodEvent, kGlobalTrk );
    if(fTPCBranch)
      ProduceArrayTrksAOD( aodEvent, kTPCTrk );
    fEvents->Fill(1);

    AliAODVZERO *esdV0 = aodEvent->GetVZEROData();// loop sobre canales del V0 para obtener las multiplicidad

    for (Int_t iCh=0; iCh<64; ++iCh) { 
      Float_t multv=esdV0->GetMultiplicity(iCh); 
      Int_t intexv=iCh;
      VZEROCell* cellv0 = new((*fVZEROArray)[iCh]) VZEROCell();
      cellv0->cellindex=intexv;
      cellv0->cellmult= multv;
    }   




  } // end if triggered
  
  if(fTreeOption) {

    //Sort(fTrackArrayGlobalPar, kFALSE);

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
void AliAnalysisTaskHighPtDeDx::ProcessMCTruthESD() 
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
    if (chargeMC == 0)
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
    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(pdgCode);
    
    // Here we want to add some of the MC histograms!
    
    // And therefore we first cut here!
    if (trackMC->Pt() < fMinPt) {
      
      // Keep small fraction of low pT tracks
      if(fRandom->Rndm() > fLowPtFraction)
	continue;
    } // else {
    // Here we want to add the high pt part of the MC histograms!
    //    }
    
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
void AliAnalysisTaskHighPtDeDx::ProcessMCTruthAOD() 
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
    if (chargeMC == 0)
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
    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(pdgCode);
    
    // Here we want to add some of the MC histograms!
    
    // And therefore we first cut here!
    if (trackMC->Pt() < fMinPt) {
      
      // Keep small fraction of low pT tracks
      if(fRandom->Rndm() > fLowPtFraction)
	continue;
    } // else {
    // Here we want to add the high pt part of the MC histograms!
    //    }
    
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
TParticle* AliAnalysisTaskHighPtDeDx::FindPrimaryMother(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
  if (motherLabel < 0)
    return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskHighPtDeDx::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // returns its label
  //
  // Taken from AliPWG0Helper class
  //
  const Int_t nPrim  = stack->GetNprimary();
  
  while (label >= nPrim) {

    //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

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
//__________________________________________________________________
void AliAnalysisTaskHighPtDeDx::ProduceArrayTrksESD( AliESDEvent *ESDevent, AnalysisMode analysisMode ){
  
  const Int_t nESDTracks = ESDevent->GetNumberOfTracks();
  Int_t trackmult=0;
  Int_t nadded=0;
  if( analysisMode == kGlobalTrk ){
    if(fTrackArrayGlobalPar)
      fTrackArrayGlobalPar->Clear();
  } else if( analysisMode == kTPCTrk ){
    if(fTrackArrayTPCPar)
      fTrackArrayTPCPar->Clear();
  }
  
  if( analysisMode==kGlobalTrk ){  
    
    for(Int_t iT = 0; iT < nESDTracks; iT++) {
      
      AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
      
      
      if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
	continue;
      
      UShort_t filterFlag = 0;
      Bool_t filterCut_Set1 = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2 = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3 = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
      
      UInt_t selectDebug = 0;
      if (fTrackFilterGolden) {
	selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
	if (selectDebug) {
	  filterFlag +=1;
	  filterCut_Set3=kTRUE;
	}
      }
      
      if (fTrackFilterTPC) {
	
	selectDebug = fTrackFilterTPC->IsSelected(esdTrack);
	if (selectDebug){//only tracks which pass the TPC-only track cuts
	  trackmult++;
	  filterFlag +=2;
	  filterCut_Set1=kTRUE;
	  
	}
	
      }
      
      if (fTrackFilter) {
	selectDebug = fTrackFilter->IsSelected(esdTrack);
	if (selectDebug) {
	  filterCut_Set2=kTRUE;
	}
      }
      
     
      if(filterFlag==0)
	continue;
      
      
      //
      // Track was accepted
      //      
      
      // Here we want to add histograms!
      
      if (esdTrack->Pt() < fMinPt) {
	
	// Keep small fraction of low pT tracks
	if(fRandom->Rndm() > fLowPtFraction)
	  continue;
      } // else {
      // Here we want to add the high pt part of the histograms!
      //    }
    
      Short_t charge  = esdTrack->Charge();
      Float_t pt      = esdTrack->Pt();
      Float_t p       = esdTrack->P(); 
      Float_t eta     = esdTrack->Eta();
      Float_t phi     = esdTrack->Phi();
      Short_t ncl     = esdTrack->GetTPCsignalN();
      Short_t neff    = Short_t(esdTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      //	  Short_t nclf    = esdTrack->GetTPCNclsF();
      Float_t dedx    = esdTrack->GetTPCsignal();
      Float_t tpcchi  = 0;
      if(esdTrack->GetTPCNcls() > 0)
	tpcchi = esdTrack->GetTPCchi2()/Float_t(esdTrack->GetTPCNcls());
      Float_t b[2];
      Float_t bCov[3];
      esdTrack->GetImpactParameters(b,bCov);
      Float_t dcaxy   = b[0];
      Float_t dcaz    = b[1];
      Double_t p_con[3] = {0, 0, 0};
      esdTrack->GetConstrainedPxPyPz(p_con);
      
      
      //	Float_t pt_con = (Float_t)TMath::Sqrt(p_con[0]*p_con[0] + p_con[1]*p_con[1]);
      // const AliExternalTrackParam* tpcParam = esdTrack->GetTPCInnerParam();
      // Float_t pttpc   = tpcParam->Pt();
      // Float_t ptpc    = tpcParam->P();
      
      Float_t ptMC        = 0;
      Short_t pidCode     = 0; // 0 = real data / no mc track!
      Short_t primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   pdgMother   = 0;
      
      
      //with Globaltrack parameters????
      if(fAnalysisMC) {
	
	const Int_t label = TMath::Abs(esdTrack->GetLabel());
	TParticle* mcTrack = fMCStack->Particle(label);	    
	if (mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(label))
	    primaryFlag = 1;
	  
	  Int_t pdgCode = mcTrack->GetPdgCode();
	  pidCode = GetPidCode(pdgCode);
	  
	  ptMC      = mcTrack->Pt();
	  
	  TParticle* mother = FindPrimaryMother(fMCStack, label);
	  pdgMother = mother->GetPdgCode();
	}
      }
      
    
      //TOF
      Float_t beta = -99;
      if (esdTrack->GetStatus()&AliESDtrack::kTOFpid){
	if ((esdTrack->GetIntegratedLength() != 0) && (esdTrack->GetTOFsignal() != 0))
	  beta = esdTrack->GetIntegratedLength()/esdTrack->GetTOFsignal()/fgkClight;
      }
      
      if(fTreeOption) {
	
	DeDxTrack* track = new((*fTrackArrayGlobalPar)[nadded]) DeDxTrack();
	nadded++;
	
	track->p          = p;
	track->pt         = pt;
	//	  track->ptcon   = pt_con;
	//	  track->tpcchi  = tpcchi;
	track->eta        = eta;
	track->phi        = phi;
	track->q          = charge;
	track->filter     = filterFlag;
	track->ncl        = ncl;
	track->neff       = neff;
	track->dedx       = dedx;
	track->beta       = beta;
	track->dcaxy      = dcaxy;
	track->dcaz       = dcaz;
	track->pid        = pidCode;
	track->primary    = primaryFlag;
	track->pttrue     = ptMC;
	track->mother     = pdgMother;
	track->filterset1 = filterCut_Set1;
	track->filterset2 = filterCut_Set2;
	track->filterset3 = filterCut_Set3;
	
	
      }
    }//end of track loop

  }//end global
  
  else if( analysisMode==kTPCTrk ){  
    const AliESDVertex *vtxSPD = ESDevent->GetPrimaryVertexSPD();
    if( vtxSPD->GetNContributors() < 1 || TMath::Abs(vtxSPD->GetZ()) > 10.0 ) return;
 

    for(Int_t iT = 0; iT < nESDTracks; iT++) {
      
      AliESDtrack* esdTrack = ESDevent->GetTrack(iT);

      AliESDtrack *trackTPC = AliESDtrackCuts::GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(ESDevent),esdTrack->GetID());
      
      if(!trackTPC) continue;
    
      UInt_t selectDebug = 0;
      selectDebug = fTrackFilterTPC->IsSelected(trackTPC);
      if(selectDebug==0) continue;
    
      if(trackTPC->Pt()>0.){
	// only constrain tracks above threshold
	AliExternalTrackParam exParam;
	// take the B-field from the ESD, no 3D fieldMap available at this point
	Bool_t relate = false;
	relate = trackTPC->RelateToVertexTPC(vtxSPD,ESDevent->GetMagneticField(),
					   kVeryBig,&exParam);
	if(!relate){
	  delete trackTPC;
	  continue;
	}
	trackTPC->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),
		      exParam.GetCovariance());
      }
      else continue;     
      
 


      if(TMath::Abs(trackTPC->Eta()) > fEtaCut)
	continue;
      
      //
      // Track was accepted
      //      
      
      // Here we want to add histograms!
      
      if (trackTPC->Pt() < fMinPt) {
	
	// Keep small fraction of low pT tracks
	if(fRandom->Rndm() > fLowPtFraction)
	  continue;
      } // else {
      // Here we want to add the high pt part of the histograms!
      //    }
      
      Short_t charge  = trackTPC->Charge();
      Float_t pt      = trackTPC->Pt();
      Float_t p       = trackTPC->P(); 
      Float_t eta     = trackTPC->Eta();
      Float_t phi     = trackTPC->Phi();
      Short_t ncl     = trackTPC->GetTPCsignalN();
      Short_t neff    = Short_t(trackTPC->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      //	  Short_t nclf    = esdTrack->GetTPCNclsF();
      Float_t dedx    = trackTPC->GetTPCsignal();
      Float_t tpcchi  = 0;
      if(trackTPC->GetTPCNcls() > 0)
	tpcchi = trackTPC->GetTPCchi2()/Float_t(trackTPC->GetTPCNcls());
      //Float_t b[2];
      //Float_t bCov[3];
      //trackTPC->GetImpactParameters(b,bCov);
      //Float_t dcaxy   = b[0];
      //Float_t dcaz    = b[1];

      Float_t dcaxy = 0.;
      Float_t dcaz = 0.;
      trackTPC->GetImpactParameters(dcaxy,dcaz);


      Double_t p_con[3] = {0, 0, 0};
      trackTPC->GetConstrainedPxPyPz(p_con);
      
    
      // Float_t pt_con = (Float_t)TMath::Sqrt(p_con[0]*p_con[0] + p_con[1]*p_con[1]);
      // const AliExternalTrackParam* tpcParam = esdTrack->GetTPCInnerParam();
      // Float_t pttpc   = tpcParam->Pt();
      // Float_t ptpc    = tpcParam->P();
      
      Float_t ptMC        = 0;
      Short_t pidCode     = 0; // 0 = real data / no mc track!
      Short_t primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   pdgMother   = 0;
      
      
      //with Globaltrack parameters????
      if(fAnalysisMC) {
	
	const Int_t label = TMath::Abs(esdTrack->GetLabel());
	TParticle* mcTrack = fMCStack->Particle(label);	    
	if (mcTrack){
	  
	  if(fMCStack->IsPhysicalPrimary(label))
	    primaryFlag = 1;
	  
	  Int_t pdgCode = mcTrack->GetPdgCode();
	  pidCode = GetPidCode(pdgCode);
	  
	  ptMC      = mcTrack->Pt();
	  
	  TParticle* mother = FindPrimaryMother(fMCStack, label);
	  pdgMother = mother->GetPdgCode();
	}
      }
    
    
      //TOF
      Float_t beta = -99;
      if (esdTrack->GetStatus()&AliESDtrack::kTOFpid){
	if ((esdTrack->GetIntegratedLength() != 0) && (esdTrack->GetTOFsignal() != 0))
	  beta = esdTrack->GetIntegratedLength()/esdTrack->GetTOFsignal()/fgkClight;
      }
      
      if(fTreeOption) {
	
	DeDxTrack* tracktpc = new((*fTrackArrayTPCPar)[nadded]) DeDxTrack();
	nadded++;
	
	tracktpc->p          = p;
	tracktpc->pt         = pt;
	//	  track->ptcon   = pt_con;
	//	  track->tpcchi  = tpcchi;
	tracktpc->eta        = eta;
	tracktpc->phi        = phi;
	tracktpc->q          = charge;
	tracktpc->filter     = 1;
	tracktpc->ncl        = ncl;
	tracktpc->neff       = neff;
	tracktpc->dedx       = dedx;
	tracktpc->beta       = beta;//computed with Global tracks
	tracktpc->dcaxy      = dcaxy;
	tracktpc->dcaz       = dcaz;
	tracktpc->pid        = pidCode;
	tracktpc->primary    = primaryFlag;
	tracktpc->pttrue     = ptMC;
	tracktpc->mother     = pdgMother;
	tracktpc->filterset1 = 0;
	tracktpc->filterset2 = 0;
	tracktpc->filterset3 = 0;
	
      }


      if(trackTPC)
	delete trackTPC;


    }//end of track loop
 
  }//end of: if isglobal==kFALSE


  if(fTreeOption) {

    if( analysisMode==kGlobalTrk ){
      Sort(fTrackArrayGlobalPar, kFALSE);
      
      fEvent->trackmult = trackmult;
      fEvent->n         = nadded;
    }
    else if( analysisMode==kTPCTrk ){
      Sort(fTrackArrayTPCPar, kFALSE);
    }

  }


}
//__________________________________________________________________
void AliAnalysisTaskHighPtDeDx::ProduceArrayTrksAOD( AliAODEvent *AODevent, AnalysisMode analysisMode ){


  Int_t     trackmult = 0; // no pt cuts
  Int_t     nadded    = 0;
  Int_t nAODTracks = AODevent->GetNumberOfTracks();
  if( analysisMode == kGlobalTrk ){
    
    if(fTrackArrayGlobalPar)
      fTrackArrayGlobalPar->Clear();
    
  } 
  if( analysisMode == kTPCTrk ){
    if(fTrackArrayTPCPar)
      fTrackArrayTPCPar->Clear();
  }




  //const AliAODVertex*	vertexSPD= (AliAODVertex*)AODevent->GetPrimaryVertexSPD();//GetPrimaryVertex()
  //if( vertexSPD->GetNContributors() < 1 || TMath::Abs( vertexSPD->GetZ()) > 10.0 ) return;


  
  if( analysisMode==kGlobalTrk ){  
 

     const AliAODVertex* vertex = AODevent->GetPrimaryVertex();
    for(Int_t iT = 0; iT < nAODTracks; iT++) {
      
      AliAODTrack* aodTrack = AODevent->GetTrack(iT);
      
 
      
      UShort_t filterFlag = 0;
      Bool_t filterCut_Set1 = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2 = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3 = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
  
      
      
      if (fTrackFilterGolden) {
	
	// ITSTPC2010 cuts is bit 32 according to above macro, new versions of aliroot includes the golden cuts
	if(aodTrack->TestFilterBit(32)) {
	  filterFlag +=1;
	  filterCut_Set3 = kTRUE;
	}
      }
      
      
      if (fTrackFilterTPC) {
	// TPC only cuts is bit 1 according to above macro
	// Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
	if(aodTrack->TestFilterBit(1)){
	  filterFlag +=2;
	  filterCut_Set1 = kTRUE;
	  trackmult++;

	}
      }
      
      
      
      if(filterFlag==0)
	continue;
      
      
      Double_t b[2], cov[3];
      if (!(aodTrack->PropagateToDCA(vertex, AODevent->GetMagneticField(), kVeryBig, b, cov)))
	filterFlag = 32; // propagation failed!!!!!;
      Float_t dcaxy   = b[0];
      Float_t dcaz    = b[1];
      




      // As I understand this routine recalculates the momentum so it should be called first!
      //Double_t b[2], cov[3];
  
      
      //if(!aodTrack->PropagateToDCA(vertex, AODevent->GetMagneticField(), kVeryBig, b, cov))
      //	filterFlag = 32; // propagation failed!!!!!
      
      if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
	continue;
      if (aodTrack->Pt() < fMinPt) {
	
	// Keep small fraction of low pT tracks
	if(fRandom->Rndm() > fLowPtFraction)
	  continue;
      } // else {
      // Here we want to add the high pt part of the histograms!
      //    }
     
     

      
      Short_t charge  = aodTrack->Charge();
      Float_t pt      = aodTrack->Pt();
      Float_t p       = aodTrack->P(); 
      Float_t eta     = aodTrack->Eta();
      Float_t phi     = aodTrack->Phi();
      AliAODPid* aodPid = aodTrack->GetDetPid();
      Short_t ncl     = -10;
      Short_t neff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      //	  Short_t nclf    = aodTrack->GetTPCNclsF();
      Float_t dedx    = -10;
      Float_t beta = -99;
      if(aodPid) {
	ncl     = aodPid->GetTPCsignalN();
	dedx    = aodPid->GetTPCsignal();
	//TOF
	if (aodTrack->GetStatus()&AliESDtrack::kTOFpid){
	  Double_t tof[5];
	  aodPid->GetIntegratedTimes(tof);
	  beta = tof[0]/aodPid->GetTOFsignal();
	}
      }
      //	Float_t tpcchi = aodTrack->Chi2perNDF();
      
      // Double_t p_con[3] = {0, 0, 0};
      // aodTrack->GetConstrainedPxPyPz(p_con);
      //	Float_t pt_con = 0; // This is not there! (Float_t)TMath::Sqrt(p_con[0]*p_con[0] + p_con[1]*p_con[1]);
      // const AliExternalTrackParam* tpcParam = aodTrack->GetTPCInnerParam();
      // Float_t pttpc   = tpcParam->Pt();
      // Float_t ptpc    = tpcParam->P();
      
      Float_t ptMC        = 0;
      Short_t pidCode     = 0; // 0 = real data / no mc track!
      Short_t primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   pdgMother   = 0;
      
      if(fAnalysisMC) {
	
	const Int_t label = TMath::Abs(aodTrack->GetLabel());
	
	AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));
	if (mcTrack){
	  
	  if(mcTrack->IsPhysicalPrimary())
	    primaryFlag = 1;
	  
	  Int_t pdgCode = mcTrack->GetPdgCode();
	  pidCode = GetPidCode(pdgCode);
	  
	  ptMC      = mcTrack->Pt();
	  
	  AliAODMCParticle* mother = FindPrimaryMotherAOD(mcTrack);
	  pdgMother = mother->GetPdgCode();	    
	}
      }
    
    
      if(fTreeOption) {
	
	DeDxTrack* track = new((*fTrackArrayGlobalPar)[nadded]) DeDxTrack();
	nadded++;
	
	track->p          = p;
	track->pt         = pt;
	//	  track->ptcon   = pt_con;
	//	  track->tpcchi  = tpcchi;
	track->eta        = eta;
	track->phi        = phi;
	track->q          = charge;
	track->filter     = filterFlag;
	track->ncl        = ncl;
	track->neff       = neff;
	track->dedx       = dedx;
	track->beta       = beta;
	track->dcaxy      = dcaxy;
	track->dcaz       = dcaz;
	track->pid        = pidCode;
	track->primary    = primaryFlag;
	track->pttrue     = ptMC;
	track->mother     = pdgMother;
	track->filterset1 = filterCut_Set1;
	track->filterset2 = filterCut_Set2;
	track->filterset3 = filterCut_Set3;
      }
    }//end of track loop

 

  }//end of global track analysis
  


  else if( analysisMode==kTPCTrk ){  
 
    const AliAODVertex*	vertexSPD= (AliAODVertex*)AODevent->GetPrimaryVertexSPD();//GetPrimaryVertex()
    if( vertexSPD->GetNContributors() < 1 || TMath::Abs(vertexSPD->GetZ()) > 10.0 ) return;


    //const AliAODVertex* vertex = AODevent->GetPrimaryVertex();
    for(Int_t iT = 0; iT < nAODTracks; iT++) {
      
      AliAODTrack* aodTrack = AODevent->GetTrack(iT);
      
 
      
      UShort_t filterFlag = 0;
      Bool_t filterCut_Set1 = kFALSE;//parameters from global tracks, with TPC cuts (filter bit =1 in AOD)
      Bool_t filterCut_Set2 = kFALSE;//parameters from global tracks, cuts tpc+its 2010 W/O golden cuts
      Bool_t filterCut_Set3 = kFALSE;//parameters from global tracks, cuts its+tpc 2010 WITH golden cuts
  
      
     
      // TPC only cuts is bit 1 according to above macro
      // Alex always uses 128, NOTE: FILTER 128 ARE TPC TRACKS (TPC PARAMETERS) CONTRAINED TO THE SPD VERTEX, 
      if(aodTrack->TestFilterBit(128)){
	filterFlag +=2;
	trackmult++;
      }
      
      
      if(filterFlag==0)
	continue;
      
      
      Double_t b[2], cov[3];
      //AliAODVertex* vertex = AODevent->GetPrimaryVertex();
      //Double_t b[2] = {-99., -99.};
      //Double_t bCov[3] = {-99., -99., -99.};
      //aodTrack->PropagateToDCA(aod->GetPrimaryVertex(), aod->GetMagneticField(), 100., b, bCov);
      if (!(aodTrack->PropagateToDCA(vertexSPD, AODevent->GetMagneticField(), kVeryBig, b, cov)))
	filterFlag = 32; // propagation failed!!!!!;
      Float_t dcaxy   = b[0];
      Float_t dcaz    = b[1];




      // As I understand this routine recalculates the momentum so it should be called first!
      //Double_t b[2], cov[3];
  
      
      //if(!aodTrack->PropagateToDCA(vertex, AODevent->GetMagneticField(), kVeryBig, b, cov))
      //	filterFlag = 32; // propagation failed!!!!!
      
      if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
	continue;
      if (aodTrack->Pt() < fMinPt) {
	
	// Keep small fraction of low pT tracks
	if(fRandom->Rndm() > fLowPtFraction)
	  continue;
      } // else {
      // Here we want to add the high pt part of the histograms!
      //    }
     
     

      
      Short_t charge  = aodTrack->Charge();
      Float_t pt      = aodTrack->Pt();
      Float_t p       = aodTrack->P(); 
      Float_t eta     = aodTrack->Eta();
      Float_t phi     = aodTrack->Phi();
      AliAODPid* aodPid = aodTrack->GetDetPid();
      Short_t ncl     = -10;
      Short_t neff    = 0; // This is not yet there! Short_t(aodTrack->GetTPCClusterInfo(2, 1)); // effective track length for pT res
      //	  Short_t nclf    = aodTrack->GetTPCNclsF();
      Float_t dedx    = -10;
      Float_t beta = -99;
      if(aodPid) {
	ncl     = aodPid->GetTPCsignalN();
	dedx    = aodPid->GetTPCsignal();
	//TOF
	if (aodTrack->GetStatus()&AliESDtrack::kTOFpid){
	  Double_t tof[5];
	  aodPid->GetIntegratedTimes(tof);
	  beta = tof[0]/aodPid->GetTOFsignal();
	}
      }
      //	Float_t tpcchi = aodTrack->Chi2perNDF();
      
      // Double_t p_con[3] = {0, 0, 0};
      // aodTrack->GetConstrainedPxPyPz(p_con);
      //	Float_t pt_con = 0; // This is not there! (Float_t)TMath::Sqrt(p_con[0]*p_con[0] + p_con[1]*p_con[1]);
      // const AliExternalTrackParam* tpcParam = aodTrack->GetTPCInnerParam();
      // Float_t pttpc   = tpcParam->Pt();
      // Float_t ptpc    = tpcParam->P();
      
      Float_t ptMC        = 0;
      Short_t pidCode     = 0; // 0 = real data / no mc track!
      Short_t primaryFlag = 0; // 0 = real data / not primary mc track  
      Int_t   pdgMother   = 0;
      
      if(fAnalysisMC) {
	
	const Int_t label = TMath::Abs(aodTrack->GetLabel());
	
	AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));
	if (mcTrack){
	  
	  if(mcTrack->IsPhysicalPrimary())
	    primaryFlag = 1;
	  
	  Int_t pdgCode = mcTrack->GetPdgCode();
	  pidCode = GetPidCode(pdgCode);
	  
	  ptMC      = mcTrack->Pt();
	  
	  AliAODMCParticle* mother = FindPrimaryMotherAOD(mcTrack);
	  pdgMother = mother->GetPdgCode();	    
	}
      }
    
    
      if(fTreeOption) {
	
	DeDxTrack* tracktpc = new((*fTrackArrayTPCPar)[nadded]) DeDxTrack();
	nadded++;
	
	tracktpc->p          = p;
	tracktpc->pt         = pt;
	//	  track->ptcon   = pt_con;
	//	  track->tpcchi  = tpcchi;
	tracktpc->eta        = eta;
	tracktpc->phi        = phi;
	tracktpc->q          = charge;
	tracktpc->filter     = filterFlag;
	tracktpc->ncl        = ncl;
	tracktpc->neff       = neff;
	tracktpc->dedx       = dedx;
	tracktpc->beta       = beta;
	tracktpc->dcaxy      = dcaxy;
	tracktpc->dcaz       = dcaz;
	tracktpc->pid        = pidCode;
	tracktpc->primary    = primaryFlag;
	tracktpc->pttrue     = ptMC;
	tracktpc->mother     = pdgMother;
	tracktpc->filterset1 = filterCut_Set1;
	tracktpc->filterset2 = filterCut_Set2;
	tracktpc->filterset3 = filterCut_Set3;
      }
    }//end of track loop


  }//end of global track analysis


  if(fTreeOption) {
    
    if( analysisMode==kGlobalTrk ){
      Sort(fTrackArrayGlobalPar, kFALSE);
      
      fEvent->trackmult = trackmult;
      fEvent->n         = nadded;


    }
    if( analysisMode==kTPCTrk ){
      Sort(fTrackArrayTPCPar, kFALSE);
    }
    
  }
  
}
