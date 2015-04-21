/*************************************************************************
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

/* Class AliAnalysisTaskSEHFJets:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/* Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it, svallero@to.infn.it, s.lapointe@cern.ch */

#include <TH1F.h>
#include <TH2F.h>
#include <TAxis.h>
#include <TArrayI.h>
#include <TArrayD.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TProfile.h>

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenerator.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskSEHFJets.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNormalizationCounter.h"
#include "AliLog.h"
#include "AliHFJetsContainerVertex.h"
#include "AliAnalysisHelperJetTasks.h"

#include "AliAnalysisUtils.h"

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"


class TCanvas;
class TTree;
class TChain;
// class AliAnalysisTaskSE;
class AliCFContainer;

ClassImp(AliAnalysisTaskSEHFJets)

AliAnalysisTaskSEHFJets::AliAnalysisTaskSEHFJets() 
: AliAnalysisTaskEmcalJet("AliAnalysisTaskSEHFJets", kTRUE),
  fhJets(0),
  fhQaVtx(0),
  fhBJets(0),
  fhJetVtx(0),
  fhJetVtxData(0),
  fEvent(0x0), 
  fCorrMode(kTRUE),
  fJetIdMeth(),
  fSelectPtHard(kFALSE),
  fUseTriggerData(kFALSE),
  fPtHardMin(),
  fFlavor(),
  fTaggingRadius(),
  fFiredClass(),
  fNentries(0),
  fNtriggers(0),
  fTagger(0),
  fCutsHFjets(0),
  fAnalysisUtils(0x0),
  fbJetArray(0),
  fArrayMC(0),
  fJetArray(0x0),
  fMcJetArray(0x0),
  fJetContName(""),
  fTrackContName(""), 
  fMcJetContName(""),
  fMcTrackContName(""),
  fGTIp(0),fGTIn(0),fTrackBuffSize(19000),
  fHistTrials(0),
  fHistXsection(0),
  fHistEvents(0),
  fOutputList(0x0),fMCTracksCont(0)
   
{

  // default constructor

}

//________________________________________________________________________
AliAnalysisTaskSEHFJets::AliAnalysisTaskSEHFJets(const char *name) 
  : AliAnalysisTaskEmcalJet(name,kTRUE),
    fhJets(0),
    fhQaVtx(0),
    fhBJets(0),
    fhJetVtx(0),
    fhJetVtxData(0),
    fEvent(0x0), 
    fCorrMode(kTRUE),
    fJetIdMeth(),
    fSelectPtHard(kFALSE),
    fUseTriggerData(kFALSE),
    fPtHardMin(),
    fFlavor(),
    fTaggingRadius(),
    fFiredClass(),
    fNentries(0),
    fNtriggers(0),
    fTagger(0),
    fCutsHFjets(0),
    fAnalysisUtils(0x0),
    fbJetArray(0),
    fArrayMC(0),
    fJetArray(0x0),
    fMcJetArray(0x0),
    fJetContName(""),
    fTrackContName(""), 
    fMcJetContName(""),
    fMcTrackContName(""),
    fGTIp(0),fGTIn(0),fTrackBuffSize(19000),
    fHistTrials(0),
    fHistXsection(0),
    fHistEvents(0),
    fOutputList(0x0),fMCTracksCont(0)
{ 
  // standard constructor
  AliInfo("+++ Executing Constructor +++");

  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
AliAnalysisTaskSEHFJets::~AliAnalysisTaskSEHFJets(){

  // destructor
  AliInfo("+++ Executing Destructor +++");

  // Do not delete outputs in proof mode or merging will fail
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    if (fOutputList) delete fOutputList;
    if (fhJets) delete fhJets;
    if (fhQaVtx) delete fhQaVtx;
    if (fhBJets) delete fhBJets;
    if (fhJetVtx) delete fhJetVtx;
    if (fhJetVtxData) delete fhJetVtxData;
  }

  if (fTagger) delete fTagger;
  if (fCutsHFjets) delete fCutsHFjets;

  // Array, note the [] with the delete
  if (fGTIp)
    delete[] fGTIp;
  fGTIp=0;
  if (fGTIn)
    delete[] fGTIn;
  fGTIn=0;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFJets::Init()
{
  // Initialization
  AliInfo("+++ Executing Init +++");

  AliLog::SetGlobalDebugLevel(AliLog::kError);

}
//________________________________________________________________________
void AliAnalysisTaskSEHFJets::UserCreateOutputObjects(){

  // AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AliInfo("+++ Executing UserCreateOutputObjects +++");

 
  // else {        //no jets, just analysis tracks and clusters
  //   fTracksCont       = GetParticleContainer(0);
  //   fCaloClustersCont = GetClusterContainer(0);
  // }
  // fTracksCont->SetClassName("AliVTrack");
  // fCaloClustersCont->SetClassName("AliAODCaloCluster");

  // Initialize output list of containers
  if (fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if (!fOutputList){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }  
  
  // Control histogram
  fNentries=new TH1F("nentriesChFr", "Analyzed sample properties", 9,-0.5,8.5);
  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nEvSel");
  fNentries->GetXaxis()->SetBinLabel(3,"nEvGoodVtx");
  fNentries->GetXaxis()->SetBinLabel(4,"nEvPile-up spd");
  fNentries->GetXaxis()->SetBinLabel(5,"nEvPile-up mv");
  fNentries->GetXaxis()->SetBinLabel(6,"nTracksEv");
  fNentries->GetXaxis()->SetBinLabel(7,"nJetsCand");
  fNentries->GetXaxis()->SetBinLabel(8,"nJetsTagged");
  fNentries->GetXaxis()->SetBinLabel(9,"nUnexpError");
  fOutputList->Add(fNentries);

  // Control histogram
  fNtriggers=new TH1F("ntriggers", "Analyzed sample properties", 5,-0.5,5.5);
  fNtriggers->GetXaxis()->SetBinLabel(1,"nJ1");
  fNtriggers->GetXaxis()->SetBinLabel(3,"nJ2");
  fNtriggers->GetXaxis()->SetBinLabel(5,"notJ1orJ2");
  fOutputList->Add(fNtriggers);

  fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
  fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistTrials->GetYaxis()->SetTitle("trials");
  fOutputList->Add(fHistTrials);

  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
  fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistEvents->GetYaxis()->SetTitle("total events");
  fOutputList->Add(fHistEvents);

  fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
  fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistXsection->GetYaxis()->SetTitle("xsection");
  fOutputList->Add(fHistXsection);

  // const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  // const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
    
  // for (Int_t i = 1; i < 12; i++) {
  //   fHistTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  //   fHistXsection->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  //   fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  // }

  // Create the containers for jet and vertex properties
  // reconstructed jets
  fhJets= new AliHFJetsContainerVertex("kJets",AliHFJetsContainerVertex::kJets);
  fOutputList->Add(fhJets);
  // vertices QA
  fhQaVtx= new AliHFJetsContainerVertex("kQaVtx",AliHFJetsContainerVertex::kQaVtx);
  fOutputList->Add(fhQaVtx);
  // tagged jets
  fhBJets= new AliHFJetsContainerVertex("kBJets",AliHFJetsContainerVertex::kBJets);
  fOutputList->Add(fhBJets);
  // vertices within the jet
  fhJetVtx= new AliHFJetsContainerVertex("kJetVtx",AliHFJetsContainerVertex::kJetVtx);
  fOutputList->Add(fhJetVtx);
  // vertices within the jet - data
  fhJetVtxData= new AliHFJetsContainerVertex("kJetVtxData",AliHFJetsContainerVertex::kJetVtxData);
  fOutputList->Add(fhJetVtxData);

  fGTIp = new AliAODTrack *[fTrackBuffSize]; // Array of pointers
  fGTIn = new AliAODTrack *[fTrackBuffSize]; // Array of pointers


  
  PostData(1,fOutputList);

}

// //________________________________________________________________________
// void AliAnalysisTaskSEHFJets::ExecOnce() {

//   AliAnalysisTaskEmcalJet::ExecOnce();

//   if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
//   if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
//   if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

// }

//________________________________________________________________________
void AliAnalysisTaskSEHFJets::UserExec(Option_t */*option*/){

  AliInfo("+++ Executing UserExec +++");


  // Execute analysis for current event
  if (fCorrMode) AnalyseCorrectionsMode(); // must be MC, all steps are filled for container kBJets (only)
  else AnalyseDataMode(); // can also be MC, only step kCFStepRecoB is filled also for kBJets


}

//________________________________________________________________________
void AliAnalysisTaskSEHFJets::AnalyseDataMode(){

  if(!fbJetArray)fbJetArray=new TClonesArray("AliAODVertex",0);
  Double_t arrDispersion[10000];

  // AOD input event
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!aod) {
    AliInfo("Input AOD not available, trying with output handler...");
    if(AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD 
      // event in memory rather than the input event.    
      aod = dynamic_cast<AliAODEvent*> (AODEvent());
    } else AliError(RED"No AOD handler found or no standard AOD!" Bee);
  }

  fEvent = aod;

  if (!GetArrays()) return;

  // ALL EVENTS
  fNentries->Fill(0); // EventsAnal 
  
  AliAODVertex *vtx1;
  vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  if(!fCutsHFjets->IsEventSelected(aod)){
    AliDebug(AliLog::kDebug,"Event did not pass event selection from AliRDHFJetsCuts!");
    return;
  }

  fNentries->Fill(1); // event selected, pileup, trigger, etc...

  Int_t mult = 0;

  // Init steps
  AliHFJetsContainer::CFSteps step;
 
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO 
  Double_t primvtx[3],primcov[6];
  vtx1->GetXYZ(primvtx);
  vtx1->GetCovarianceMatrix(primcov);
  Int_t nPrimContr=vtx1->GetNContributors();
  Double_t chi2=vtx1->GetChi2();
  AliESDVertex* v1 = new AliESDVertex(primvtx,primcov,chi2,nPrimContr);
  Double_t magzkG = (Double_t)aod->GetMagneticField();


  // Reset the reference array to the global tracks..
  ResetTrackReference();
  // ..and set it
  for (Int_t iTrack=0;iTrack<aod->GetNumberOfTracks();iTrack++){
    AliAODTrack *track = static_cast <AliAODTrack*> (aod->GetTrack(iTrack));
    if (!track) continue;
    
    // Store the reference of the global tracks
    StoreTrackReference(track);
  }

  // Loop on jets (clusterized on RECO particles)
  // Int_t oldnJets=arrayJets->GetEntries();
  Int_t nJets=fJetArray->GetEntries();
  // cout << "number of jets new = " << nJets << endl;

  AliEmcalJet *jet;
  Int_t nvtx=0;
 

  for(Int_t jetcand=0;jetcand<nJets;jetcand++){
    nvtx=0;
    jet=(AliEmcalJet*)fJetArray->UncheckedAt(jetcand);
    if (jet->Pt()>10) cout << "new jet pt is = " << jet->Pt() << endl;
    TClonesArray *fTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject("PicoTracks"));
    if (!fTrackArray) cout << "no track array"<< endl;
    Int_t numtracks = fTrackArray->GetEntries();

    if(!fCutsHFjets->IsJetSelected(jet)){
      AliDebug(AliLog::kDebug,Form("Jet not selected: pT=%f, eta=%f!", jet->Pt(),jet->Eta()));
      continue;
    }
    step = AliHFJetsContainer::kCFStepRecoB;
    // // Run b-tagger
    nvtx=fTagger->FindVertices(jet,fGTIp,fGTIn,fTrackArray,aod,v1,magzkG,fbJetArray,arrDispersion);
    if(1){ //*** NEW BY SV!!! ***
      
      fhJetVtxData->FillStepJetVtxData(step,mult,jet,fbJetArray,nvtx,vtx1,arrDispersion);

      fbJetArray->Clear();
    } else AliDebug(AliLog::kDebug,"*** nvtx=0 !! ***");
    
  }
  PostData(1,fOutputList);

  // delete v1;

}
//________________________________________________________________________
void AliAnalysisTaskSEHFJets::AnalyseCorrectionsMode(){
  
  if(!fbJetArray)fbJetArray=new TClonesArray("AliAODVertex",0);
  Double_t arrDispersion[10000];

  // AOD input event
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!aod) {
    AliInfo("Input AOD not available, trying with output handler...");
    if(AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD 
      // event in memory rather than the input event.    
      aod = dynamic_cast<AliAODEvent*> (AODEvent());
    } else AliError(RED"No AOD handler found or no standard AOD!" Bee);
  }

  fEvent = aod;

  if (!GetArrays()) return;
  
  // ALL EVENTS
  fNentries->Fill(0); // EventsAnal 

  // Set flags for event selection
  Bool_t flagTriggered=kFALSE;
  Bool_t flagVertex=kFALSE;
  
  AliAODVertex *vtx1;
  vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();  
  
  if(!fCutsHFjets->IsEventSelected(aod)){
    AliDebug(AliLog::kDebug,"Event did not pass event selection from AliRDHFJetsCuts!");
    return;
  }
  flagTriggered=kTRUE;
  fNentries->Fill(1); // event selected, pileup, trigger, etc...

  // // VERTEX SELECTION
  // // AOD primary vertex // pp
  vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
  // require "VertexerTracks" and at least 1 contributor
  if(!(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)) {
    AliDebug(AliLog::kDebug,"Event did not pass primary vertex selection!");
    return;
  }
     
  // flagVertex=kTRUE;
  fNentries->Fill(2); // EvGoodVtx

// PILEUP REJECTION - settings for pp 3,0.6,3.,2.,10.
  if(aod->IsPileupFromSPD(3.,0.6,3.,2.,10.)) {
    AliDebug(AliLog::kDebug,"Event did not pass SPD pileup!");
    return;
  }
  fNentries->Fill(3);
  
  // fill ntrials and xsec, at the momemt for 1 pt hard bin
  AliAODMCHeader *aodmcHeader=0x0;

  // load MC header
  aodmcHeader =
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!aodmcHeader) AliError(RED"MC header branch not found!" Bee);

   if(fSelectPtHard){
    TString title = aodmcHeader->GetGeneratorName();
    // cout << "title = " << title << endl;
    if (!title.Contains(fPtHardMin) || !title.Contains(fFlavor)) return; 
  }
  
  fNentries->Fill(4); // is the pt-hard en
  flagVertex=kTRUE;
  
  // We should fill the container kBJets at each step for the correction procedure!
  
  // Choose method to tag MC jets 
  Int_t meth=fJetIdMeth;

  // Multiplicity MC (for vertex reco correction)
  // Get array of MC particles
  fArrayMC = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!fArrayMC) AliError(RED"MC particles branch not found!" Bee);
  // Count number of primary MC partcles
  Int_t nMC = fArrayMC->GetEntries();
  Int_t multMC = 0;
  for (Int_t i = 0; i<nMC; i++){
    AliAODMCParticle *part = (AliAODMCParticle*)fArrayMC->At(i);
    if (part->IsPhysicalPrimary()) multMC++;
  }
  AliInfo(Form(cy"MC particles %d primaries %d" Bee,nMC, multMC ));

  // Init steps
  AliHFJetsContainer::CFSteps step;
  // Loop on MC jets
  Int_t nMCJets=fMcJetArray->GetEntries();
    
  // also write jets on TList, needed for matching 
  TList *listMCJets=new TList();

  TArrayI mcBJets(nMCJets); // = new TArrayI(nMCJets);
  mcBJets.Reset(0);

  Int_t count=-1;
  AliEmcalJet *jetMC;
  for(Int_t jetcand=0;jetcand<nMCJets;jetcand++){
    jetMC=(AliEmcalJet*)fMcJetArray->UncheckedAt(jetcand);
    if (!jetMC) continue; 
    // restrict jet eta and pT ranges
    if(!fCutsHFjets->IsJetSelected(jetMC)){
      AliDebug(AliLog::kDebug,Form("JetMC not selected: pT=%f, eta=%f!", jetMC->Pt(),jetMC->Eta()));

      continue;
    }
    
    count++;
    listMCJets->AddAt(jetMC,count);

    TClonesArray *fTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject("MCParticlesSelected"));
	
    // Get jet flavour from 2 methods
    Double_t partonnatMC[2];
    Double_t ptpartMC[2];
    Double_t contributionMC=0; // pT weight of mother parton (only method 1)
   
    GetFlavour2Methods(jetMC, partonnatMC, ptpartMC, contributionMC, fTaggingRadius);

    // this is for g,ud,s,c,b
    mcBJets[count]=partonnatMC[meth]; 
    Printf(RED"MC mcBJets: %d  partonnatMC: %f position: %d" Bee, mcBJets.At(count), partonnatMC[meth], count);
    // Fill container tagger
    step=AliHFJetsContainer::kCFStepAll;
    //fhBJets->FillStepBJets(step,multMC,jetMC,0,partonnatMC,contributionMC,ptpartMC[0]);
    if (flagTriggered) {
      step=AliHFJetsContainer::kCFStepTriggered;
      //fhBJets->FillStepBJets(step,multMC,jetMC,0,partonnatMC,contributionMC,ptpartMC[0]);
      if (flagVertex) {
	step=AliHFJetsContainer::kCFStepVertex;
	// fhBJets->FillStepBJets(step,multMC,jetMC,0,partonnatMC,contributionMC,ptpartMC[0]);
	//   	//fhJets->FillStepJets(step,multMC,jetMC,partonnatMC,contributionMC,ptpartMC);
	fhJetVtx->FillStepJetVtx(step,multMC,jetMC,0,0,0,0,partonnatMC,arrDispersion);
      }
    }
  } // end loop on jets

  if (!flagVertex) return;
 
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO 
  Double_t primvtx[3],primcov[6];
  vtx1->GetXYZ(primvtx);
  vtx1->GetCovarianceMatrix(primcov);
  Int_t nPrimContr=vtx1->GetNContributors();
  Double_t chi2=vtx1->GetChi2();
  AliESDVertex* v1 = new AliESDVertex(primvtx,primcov,chi2,nPrimContr);
  Double_t magzkG = (Double_t)aod->GetMagneticField();

  // Reset the reference array to the global tracks..
  ResetTrackReference();
  // ..and set it
  for (Int_t iTrack=0;iTrack<aod->GetNumberOfTracks();iTrack++){
    AliAODTrack *track = static_cast <AliAODTrack*> (aod->GetTrack(iTrack));
    if (!track) continue;
    
    // Store the reference of the global tracks
    StoreTrackReference(track);
  }

  
  TClonesArray *arrayMC=0x0;
  Double_t vtxTrue[3];

  // load MC particles
  arrayMC = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC) AliError(RED"MC particles branch not found!" Bee);


  // MC primary vertex
  aodmcHeader->GetVertex(vtxTrue);

  // Loop on jets (clusterized on RECO particles)
   
  Int_t nJets=fJetArray->GetEntries();
  AliEmcalJet *jet;

  Int_t nvtx=0;
 
  for(Int_t jetcand=0;jetcand<nJets;jetcand++){
    nvtx=0;
    jet=(AliEmcalJet*)fJetArray->UncheckedAt(jetcand);

    if(!fCutsHFjets->IsJetSelected(jet)){
      AliDebug(AliLog::kDebug,Form("Jet not selected: pT=%f, eta=%f!", jet->Pt(),jet->Eta()));
      continue;
    }  // !!!!

    TClonesArray *fTrackArrayMc = dynamic_cast<TClonesArray*>(fEvent->FindListObject("MCParticlesSelected"));
    TClonesArray *fTrackArrayRec = dynamic_cast<TClonesArray*>(fEvent->FindListObject("PicoTracks"));
 
    // Get jet flavour from 3 methods
    Double_t partonnat[2];
    Double_t ptpart[2];
    Double_t contribution=0; // pT weight of mother parton (only method 1)
    GetFlavour2Methods(jet, partonnat, ptpart, contribution, fTaggingRadius);

    step = AliHFJetsContainer::kCFStepRecoB;
    // Fill container jets
    fhJets->FillStepJets(step,multMC,jet,partonnat,contribution,ptpart);
    // Run b-tagger
    nvtx=fTagger->FindVertices(jet,fGTIp,fGTIn,fTrackArrayRec,aod,v1,magzkG,fbJetArray,arrDispersion);
    printf(" %d vertices, %d array size\n",nvtx,fbJetArray->GetEntries());
    // if(nvtx>0){ // sll change back to 1 !!!!!!!!!!!!
      if(1){ //*** NEW BY SV!!! ***
      // QA primary vertex selection
      fhQaVtx->FillStepQaVtx(step,multMC,jet,fbJetArray,arrDispersion,nvtx,vtx1,arrayMC,partonnat);
      // Fill container vertices
      fhJetVtx->FillStepJetVtx(step,multMC,jet,fbJetArray,nvtx,vtx1,arrayMC,partonnat,arrDispersion);
      // Fill container tagger
      fhBJets->FillStepBJets(step,multMC,jet,nvtx,partonnat,contribution,ptpart[0]);

      // Jet Matching
      TList *listTaggedJets=new TList();
      listTaggedJets->SetOwner(kTRUE);
      listTaggedJets->AddAt(jet, 0);

      // Int_t genJets=100; // consider all generated jets
      // Int_t recJets=1; // consider only current reco jet
      // TArrayI matchIndex(1);
      // TArrayF pTFraction(1);
      // Int_t debug=100; // default is 0
      // Float_t maxDist=0.3; // default is 0.3
      // Int_t mode=2; // default is 1
	
      // GetJetMatching(listMCJets, genJets, fTrackArrayMc, listTaggedJets, recJets, fTrackArrayRec, matchIndex, pTFraction, debug, maxDist, mode);
      // Int_t index=matchIndex.At(0);
      // Printf(MAG"JET %d INDEX %d pT fraction %f" Bee, recJets, index, pTFraction.At(0));
      // AliEmcalJet * matchedJet;
      //       if (index >= 0){
      // 	matchedJet=(AliAODJet*)listMCJets->At(index);
      // 	// below is to avoid matching more RECO jets to the same MC
      // 	listMCJets->RemoveAt(index);

    AliEmcalJet *matchedjet = NULL;

      matchedjet = jet->ClosestJet();
      if(!matchedjet) continue;

      GetFlavour2Methods(matchedjet, partonnat, ptpart, contribution, fTaggingRadius); 

      	// for purity
      	step = AliHFJetsContainer::kCFStepMatchedAny;
      	// //fhJets->FillStepJets(step,multMC,matchedJet,partonnat,contribution,ptpart);
      	fhJetVtx->FillStepJetVtx(step,multMC,matchedjet,fbJetArray,nvtx,vtx1,arrayMC,partonnat,arrDispersion); 
      	// efficiency for different flavours
      	Int_t nature = partonnat[meth];
      	switch (nature)
      	  {
      	  case 1: // gluons
                  step = AliHFJetsContainer::kCFStepMatchedGluon;
                  Printf(MAG"Matcehd to gluon-jet!!!");
      	    break;
      	  case 2: // light quarks
                  step = AliHFJetsContainer::kCFStepMatchedLight;
                  Printf(MAG"Matcehd to light-jet!!!");
      	    break;
      	  case 3: // C 
                  step = AliHFJetsContainer::kCFStepMatchedC;
                  Printf(MAG"Matcehd to C-jet!!!");
      	    break;
      	  case 4: // B 
                  step = AliHFJetsContainer::kCFStepMatchedB;
                  Printf(MAG"Matcehd to B-jet!!!");
      	    break;
      	  default: // 0 is not defined, should never be the case...
                  AliInfo(RED"Matching jet flavour not defined!");
      	    return;
      	  }

      	fhJetVtx->FillStepJetVtx(step,multMC,matchedjet,fbJetArray,nvtx,vtx1,arrayMC,partonnat,arrDispersion); 

        //     }
            fbJetArray->Clear();
    } else AliDebug(AliLog::kDebug,"*** nvtx=0 !! ***");
  }

    PostData(1,fOutputList);

  delete v1;
  delete listMCJets;
  //delete fArrayMC; 
  //delete fCutsHFjets;
  //delete fTagger;
  //fbJetArray->Delete();
  //delete fbJetArray;
  //delete fhJets;
  //delete fhQaVtx;
  //delete fhBJets;
  //delete fhJetVtx;
  //delete fArrayMC;
  //delete aod;

}

void AliAnalysisTaskSEHFJets::Terminate(const Option_t*){

  //TERMINATE METHOD: NOTHING TO DO
  AliInfo("+++ Executing Terminate +++");

}

void AliAnalysisTaskSEHFJets::GetFlavour2Methods(AliEmcalJet *jet, Double_t (&partonnat)[2], Double_t (&ptpart)[2], Double_t &contribution, Double_t radius){

  // 3 methods  to associate jet to mother parton

  /* Nature of the parton (methods 1)    *
   * 1 = gluon (pdg 21)                    *
   * 2 = light quark (pdg < 4)             *
   * 3 = c (pdg 4)                         *
   * 4 = b (pdg 5)                         *
   * Nature of the meson/barion (method 2) *
   * 2 = light                             *
   * 3 = with c                            *
   * 4 = with b                            */

  // Initialize output values
  for (Int_t i=0; i<2; i++){
    partonnat[i]=0;
    ptpart[i] = -1.;
  }
  AliAODMCParticle *parton[2];

  parton[0]=(AliAODMCParticle*)fTagger->IsMCJetParton(fArrayMC,jet,radius); // method 2

  parton[1]=(AliAODMCParticle*)fTagger->IsMCJetMeson(fArrayMC,jet,radius); // method 3


  if(parton[0]!=0){
    Int_t pdg=TMath::Abs(parton[0]->PdgCode());
    //if(pdg==4 || pdg==5) 
    AliInfo(Form(cy"parton method -> pdg parton: %d" Bee,pdg));
    if(pdg==21)partonnat[0]=1;
    else if(pdg<4)partonnat[0]=2;
    else if(pdg==4)partonnat[0]=3;
    else if(pdg==5)partonnat[0]=4;
    ptpart[0]=parton[0]->Pt();
  }

  if(parton[1]!=0){
    Int_t pdg=TMath::Abs(parton[1]->PdgCode());
    //if((pdg>=400 && pdg<=600) || (pdg>=4000 && pdg<=6000))
    AliInfo(Form(cy"meson method -> pdg parton: %d" Bee,pdg));
    if((pdg>=400 && pdg<=500) || (pdg>=4000 && pdg<=5000))partonnat[1]=3;
    else{
      if((pdg>=500 && pdg<=600) || (pdg>=5000 && pdg<=6000))partonnat[1]=4;
    }
    ptpart[1]=parton[1]->Pt();
  }
  else partonnat[1]=2;

}
Bool_t AliAnalysisTaskSEHFJets::GetArrays(){

  // Get jet collection
  if (!fJetContName.IsNull())
    {
      AliInfo(Form("Retrieve jets %s!", fJetContName.Data()));
      // cout << "woot woot " << endl;
      fJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fJetContName));
      if (!fJetArray)
	{
	  AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetContName.Data()));
	  return kFALSE;
	}
    }
  // Get Mc jet collection
  
  if (!fMcJetContName.IsNull() && fCorrMode)
    {
      AliInfo(Form("Retrieve jets %s!", fMcJetContName.Data()));
      // cout << "woot woot " << endl;
      fMcJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMcJetContName));
      if (!fMcJetArray)
	{
	  AliError(Form("%s: Could not retrieve MC jets %s!", GetName(), fMcJetContName.Data()));
	  return kFALSE;
	}
    }

  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFJets::StoreTrackReference(AliAODTrack *track){
  // Stores the pointer to the global (positive) track
  
  // Check that the id is positive
  if(track->GetID()>-1){
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    // return;
    // Check id is not too big for buffer
    if(track->GetID()>=fTrackBuffSize){
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	     ,track->GetID(),fTrackBuffSize);
      return;
    }

    // Warn if we overwrite a track
    if(fGTIp[track->GetID()]){
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if( (!track->GetFilterMap()) &&
	  (!track->GetTPCNcls())   )
	return;

      // Imagine the other way around,
      // the zero map zero clusters track
      // is stored and the good one wants 
      // to be added. We ommit the warning
      // and just overwrite the 'bad' track
      if( fGTIp[track->GetID()]->GetFilterMap() ||
	  fGTIp[track->GetID()]->GetTPCNcls()   ){
	// If we come here, there's a problem
	printf("Warning! global track info already there!");
	printf("         TPCNcls track1 %u track2 %u",
	       (fGTIp[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
	printf("         FilterMap track1 %u track2 %u\n",
	       (fGTIp[track->GetID()])->GetFilterMap(),track->GetFilterMap());
      }
    } // Two tracks same id

    // Assign the pointer
    (fGTIp[track->GetID()]) = track;
  }


  // Check that the id is negative
  if(track->GetID()<0){
    // printf("Warning: track has negative ID: %d\n",track->GetID());
    // return;
    // Check id is not too big for buffer
    if(TMath::Abs(track->GetID())>=fTrackBuffSize){
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	     ,track->GetID(),fTrackBuffSize);
      return;
    }

    // Warn if we overwrite a track
    if(fGTIn[TMath::Abs(track->GetID())]){
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if( (!track->GetFilterMap()) &&
	  (!track->GetTPCNcls())   )
	return;

      // Imagine the other way around,
      // the zero map zero clusters track
      // is stored and the good one wants 
      // to be added. We omit the warning
      // and just overwrite the 'bad' track
      if( fGTIn[TMath::Abs(track->GetID())]->GetFilterMap() ||
	  fGTIn[TMath::Abs(track->GetID())]->GetTPCNcls()   ){
	// If we come here, there's a problem
	printf("Warning! global track info already there!");
	printf("         TPCNcls track1 %u track2 %u",
	       (fGTIn[TMath::Abs(track->GetID())])->GetTPCNcls(),track->GetTPCNcls());
	printf("         FilterMap track1 %u track2 %u\n",
	       (fGTIn[TMath::Abs(track->GetID())])->GetFilterMap(),track->GetFilterMap());
      }
    } // Two tracks same id

    // Assign the pointer
    (fGTIn[TMath::Abs(track->GetID())]) = track;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEHFJets::ResetTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTIp[i]=0;
    fGTIn[i]=0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskSEHFJets::GetJetMatching(const TList *genJetsList, const Int_t &kGenJets, TClonesArray *fTrackArraymc,
					     const TList *recJetsList, const Int_t &kRecJets, TClonesArray *fTrackArrayrec,
					     TArrayI &iMatchIndex, TArrayF &fPtFraction,
					     Int_t iDebug, Float_t maxDist, Int_t mode){

                                            
  // Matching jets from two lists
  // Start with highest energetic jet from first list (generated/embedded)
  // Calculate distance (\Delta R) to all jets from second list (reconstructed)
  // Select N closest jets = kClosestJetsN
  // Check energy fraction from jets from first list in jets from second list
  // Matched jets = jet with largest energy fraction
  // Store index of matched jet in TArrayI iMatchIndex
                                  
  // reset index
  iMatchIndex.Reset(-1);
  fPtFraction.Reset(-1.);
    
  // N closest jets: store list with index and \Delta R
  const Int_t kClosestJetsN = 4; 
  Double_t closestJets[kClosestJetsN][2]; //[][0] = index, [][1] = \Delta R
        
  const Int_t nGenJets = TMath::Min(genJetsList->GetEntries(),kGenJets);
  const Int_t nRecJets = TMath::Min(recJetsList->GetEntries(),kRecJets);
  if(nRecJets==0||nGenJets==0) {
    if(iDebug>10) Printf("No jets nRecJets %d nGenJets %d\n",nRecJets,nGenJets);
    return;
  }
  AliEmcalJet *genJet = 0x0;
  AliEmcalJet *recJet = 0x0;
    
  // loop over generated/embedded jets
  for(Int_t ig=0; ig<nGenJets; ++ig){

    for(Int_t i=0; i<kClosestJetsN; ++i){
      closestJets[i][0] = -1;   // index
      closestJets[i][1] = 1e6;  // delta R
    }

    genJet = (AliEmcalJet*)genJetsList->At(ig);
    //if(!genJet || !JetSelected(genJet)) continue;
    if(!genJet) {
      if(iDebug>10) Printf("genJet %d does not exist",ig);
      continue;
    }
        
    // find N closest reconstructed jets
    Double_t deltaR = 0.;
    for(Int_t ir=0; ir<nRecJets; ++ir){
      recJet = (AliEmcalJet*)recJetsList->At(ir);
      //if(!recJet || !JetSelected(recJet)) continue;
      if(!recJet) {
	if(iDebug>10) Printf("recJet %d doesnot exist",ir);
	continue;
      }
      deltaR = genJet->DeltaR(recJet);
      cout << "delta R = " << deltaR << endl;
      Int_t i=kClosestJetsN-1;
      if(deltaR<closestJets[i][1] && deltaR<maxDist){
	closestJets[i][0] = (Double_t) ir; // index
	closestJets[i][1] = deltaR;
                
	// sort array (closest at the beginning)
	while(i>=1 && closestJets[i][1]<closestJets[i-1][1]){
	  Double_t tmpArr[2];
	  for(Int_t j=0; j<2; j++){
	    tmpArr[j] = closestJets[i-1][j];
	    closestJets[i-1][j]   = closestJets[i][j];
	    closestJets[i][j] = tmpArr[j];
	  }
	  i--;
	}
      } 
    } // end: loop over reconstructed jets
        
    // calculate fraction for the N closest jets
    Double_t maxFraction = -1.; // maximum found fraction in one jets
    Double_t cumFraction = 0.; // cummulated fraction of closest jets (for break condition)
    Double_t fraction = 0.;
    Int_t ir = -1;  // index of close reconstruced jet
        
    for(Int_t irc=0; irc<kClosestJetsN; irc++){
      ir = (Int_t)(closestJets[irc][0]);
      if(ir<0 || ir>nRecJets-1) continue;
      recJet = (AliEmcalJet*)recJetsList->At(ir);
      if(!(recJet)) continue;
            
      fraction = GetFractionOfJet(recJet, genJet, fTrackArrayrec, fTrackArraymc, mode);
            
      cumFraction += fraction;
            
      // check if jet fulfills current matching condition
      if(fraction>maxFraction){
	// avoid multiple links
	for(Int_t ij=0; ij<ig; ++ij){
	  if(iMatchIndex[ij]==ir) continue;
	}
	// set index
	maxFraction = fraction;
	fPtFraction[ig] = fraction;                
	iMatchIndex[ig] = ir;
      }
      // break condition: less energy left as already found in one jet or
      // as required for positiv matching
      if(1-cumFraction<maxFraction) break;
    } // end: loop over closest jets
        
    if(iMatchIndex[ig]<0){
      if(iDebug) Printf("Matching failed for (gen) jet #%d", ig);
    }
  }
}

Double_t AliAnalysisTaskSEHFJets::GetFractionOfJet(const AliEmcalJet *recJet, const AliEmcalJet *genJet, TClonesArray *fTrackArrayrec, TClonesArray *fTrackArraymc, Int_t mode){
  //
  // get the fraction of the signal jet in the full jt
  //
  Double_t ptGen = genJet->Pt();
  if(ptGen==0.) return 999.;
    
  Double_t ptAssocTracks = 0.; // sum of pT of tracks found in both jets
    
  // look at tracks inside jet
  Int_t nTracksGenJet = genJet->GetNumberOfTracks();
  Int_t nTracksRecJet = recJet->GetNumberOfTracks();
  cout << "no. of gen jets = " << nTracksGenJet << ", and no. rec jets " << nTracksRecJet << endl;
    for(Int_t ir=0; ir<nTracksRecJet; ++ir){
      AliVParticle* recTrack=((AliPicoTrack*) recJet->TrackAt(ir,fTrackArrayrec))->GetTrack();
        if(!recTrack) continue;
        for(Int_t ig=0; ig<nTracksGenJet; ++ig){
	  AliVParticle* genTrack =  (AliAODMCParticle*)  genJet->TrackAt(ig,fTrackArraymc);
            if(!genTrack) continue;
	    
      // look if it points to the same track
      if( (mode&1)!=0 && genTrack==recTrack){
	cout << "where am " << endl;
	ptAssocTracks += genTrack->Pt();
	break;
      }
 
      if( (mode&2)!=0 
	  && genTrack->GetLabel()>-1
	  && recTrack->GetLabel()>-1
	  && genTrack->GetLabel()==recTrack->GetLabel()){
	cout << "i?" << endl;
	ptAssocTracks += genTrack->Pt();
	break; 
      }
    }
  }
    
  // calculate fraction
  Double_t fraction = ptAssocTracks/ptGen;
  cout << "pt assoc 2 = " << ptAssocTracks <<  endl;
  cout << "pt gen 2 = " << ptGen <<  endl;
  cout << "fraction = " << fraction <<  endl;
    
  return fraction;
}
//------------------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEHFJets::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials)
{
  //
  // Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // Get the pt hard bin from the file path
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // (Partially copied from AliAnalysisHelperJetTasks)

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1,Form("File name: %s",file.Data()));

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); 
  
  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if (!key) {
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSEHFJets::UserNotify()
{
  // Called when file changes.

  if (!fCorrMode)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  AliAODMCHeader *aodmcHeader=0x0;
  // load MC header
  aodmcHeader =
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!aodmcHeader) AliError(RED"MC header branch not found!" Bee);

    TString title = aodmcHeader->GetGeneratorName();
    if (!title.Contains(fPtHardMin) || !title.Contains(fFlavor)) return kFALSE; 

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  // get pthard from the title, where bin 1 corresponds to pt-hard min of 10, bin 2 to pt-hard min 11(lf), 18(hf-b), 20(hf-c) and so on
  pthardbin = fPtHardMin.Atoi()/10.;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials);

  fHistTrials->Fill(pthardbin, trials);
  fHistXsection->Fill(pthardbin, xsection);
  fHistEvents->Fill(pthardbin, nevents);

  return kTRUE;
}

