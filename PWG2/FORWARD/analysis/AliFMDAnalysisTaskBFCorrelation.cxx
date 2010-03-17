 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "AliFMDAnalysisTaskBFCorrelation.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "TMath.h"
#include "AliFMDAnaParameters.h"
//#include "AliFMDGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
//#include "TDatabasePDG.h"
//#include "TParticlePDG.h"
#include "AliESDInputHandler.h"
ClassImp(AliFMDAnalysisTaskBFCorrelation)


AliFMDAnalysisTaskBFCorrelation::AliFMDAnalysisTaskBFCorrelation()
: fDebug(0),
  fOutputList(0),
  fInputList(0),
  fVertexString(0x0),
  fStandalone(kTRUE)
{
  // Default constructor
  DefineInput (0, TList::Class());
  DefineOutput(0, TList::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskBFCorrelation::AliFMDAnalysisTaskBFCorrelation(const char* name, Bool_t SE):
  AliAnalysisTask(name,name),
    fDebug(0),
    fOutputList(0),
    fInputList(0),
    fVertexString(0x0),
    fStandalone(kTRUE)
{
  fStandalone = SE;
  if(fStandalone) {
    DefineInput (0, TList::Class());
    DefineInput(1, TObjString::Class());
    DefineOutput(0, TList::Class());
    
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::CreateOutputObjects()
{
  //AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  if(!fOutputList) {
    fOutputList = new TList();
    fOutputList->SetName("BackgroundCorrected");
  }
  
  TH1F* test = new TH1F("test","test",10,0,10);
  fOutputList->Add(test);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ConnectInputData(Option_t */*option*/)
{
  if(fStandalone) {
    fInputList   = (TList*)GetInputData(0);
    fVertexString = (TObjString*)GetInputData(1);
  }
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  fVertexString = (TObjString*)fInputList->At(0);
  
  Int_t vtxbin   = fVertexString->GetString().Atoi();
  //for(UShort_t det=1;det<=3;det++) {
  //  Int_t nRings = (det==1 ? 1 : 2);
  //  for (UShort_t ir = 0; ir < nRings; ir++) {
  //    Char_t ringChar = (ir == 0 ? 'I' : 'O');
  //TH2F* hMultTotal = (TH2F*)fOutputList->FindObject(Form("dNdeta_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
  // TH2F* hMultTotalTrVtx = (TH2F*)fOutputList->FindObject(Form("dNdetaTrVtx_FMD%d%c_vtxbin%d",det,ringChar,vtxbin));
  
  
  TH1F* test = (TH1F*)fOutputList->FindObject("test");
  test->Fill(vtxbin);
  if(pars->GetProcessPrimary())
    ProcessPrimary();
  
  if(fStandalone) {
    PostData(0, fOutputList); 
  }
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::Terminate(Option_t */*option*/) {
  
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskBFCorrelation::ProcessPrimary() {
  /*
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if(!mcEvent)
    return;
  
    
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcEvent->Stack();
  
  TH1F* hPrimary = (TH1F*)fOutputList->FindObject("hMultvsEta");
  AliHeader* header            = mcEvent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  if(TMath::Abs(vertex.At(2)) > pars->GetVtxCutZ())
    return;
  Double_t delta           = 2*pars->GetVtxCutZ()/pars->GetNvtxBins();
  Double_t vertexBinDouble = (vertex.At(2) + pars->GetVtxCutZ()) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
    
  Bool_t firstTrack = kTRUE;
  
  // we loop over the primaries only unless we need the hits (diagnostics running slowly)
  Int_t nTracks = stack->GetNprimary();
  if(pars->GetProcessHits())
    nTracks = stack->GetNtrack();
  
  for(Int_t i = 0 ;i<nTracks;i++) {
    particle = (AliMCParticle*) mcEvent->GetTrack(i);
    if(!particle)
      continue;
   
    if(stack->IsPhysicalPrimary(i) && particle->Charge() != 0) {
      hPrimary->Fill(particle->Eta());
      

      TH1F* hPrimVtxBin = (TH1F*)fOutputList->FindObject(Form("primmult_vtxbin%d",vertexBin));
      hPrimVtxBin->Fill(particle->Eta());
      if(firstTrack) {
	fNMCevents.Fill(vertexBin);
	firstTrack = kFALSE;
      }
    }
  }
  */
}
//_____________________________________________________________________
//
//
// EOF
