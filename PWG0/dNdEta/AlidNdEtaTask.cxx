/* $Id$ */

#include "AlidNdEtaTask.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TObjString.h>
#include <TF1.h>
#include <TGraph.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliMultiplicity.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDInputHandlerRP.h>
#include <AliESDHeader.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "AliCorrection.h"
#include "AliCorrectionMatrix3D.h"
#include "dNdEta/dNdEtaAnalysis.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"

//#define FULLALIROOT

#ifdef FULLALIROOT
  #include "../ITS/AliITSRecPoint.h"
  #include "AliCDBManager.h"
  #include "AliCDBEntry.h"
  #include "AliGeomManager.h"
  #include "TGeoManager.h"
#endif

ClassImp(AlidNdEtaTask)

AlidNdEtaTask::AlidNdEtaTask(const char* opt) :
  AliAnalysisTaskSE("AlidNdEtaTask"),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode((AliPWG0Helper::AnalysisMode) (AliPWG0Helper::kTPC | AliPWG0Helper::kFieldOn)),
  fTrigger(AliTriggerAnalysis::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(-1),
  fReadMC(kFALSE),
  fUseMCVertex(kFALSE),
  fOnlyPrimaries(kFALSE),
  fUseMCKine(kFALSE),
  fSymmetrize(kFALSE),
  fMultAxisEta1(kFALSE),
  fDiffTreatment(AliPWG0Helper::kMCFlags),
  fEsdTrackCuts(0),
  fdNdEtaAnalysisESD(0),
  fMult(0),
  fMultVtx(0),
  fEvents(0),
  fVertexResolution(0),
  fdNdEtaAnalysis(0),
  fdNdEtaAnalysisND(0),
  fdNdEtaAnalysisNSD(0),
  fdNdEtaAnalysisOnePart(0),
  fdNdEtaAnalysisTr(0),
  fdNdEtaAnalysisTrVtx(0),
  fdNdEtaAnalysisTracks(0),
  fPartPt(0),
  fVertex(0),
  fVertexVsMult(0),
  fPhi(0),
  fRawPt(0),
  fEtaPhi(0),
  fModuleMap(0),
  fDeltaPhi(0),
  fDeltaTheta(0),
  fFiredChips(0),
  fTrackletsVsClusters(0),
  fTrackletsVsUnassigned(0),
  fStats(0),
  fStats2(0),
  fPtMin(0.15),
  fEta(0x0),
  fEtaMC(0x0),
  fHistEvents(0),
  fHistEventsMC(0),
  fTrigEffNum(0),
  fTrigEffDen(0),
  fVtxEffNum(0),
  fVtxEffDen(0),
  fVtxTrigEffNum(0),
  fVtxTrigEffDen(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  //DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  
  fZPhi[0] = 0;
  fZPhi[1] = 0;

  AliLog::SetClassDebugLevel("AlidNdEtaTask", AliLog::kWarning);
}

AlidNdEtaTask::~AlidNdEtaTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}

Bool_t AlidNdEtaTask::UserNotify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file: %s", count, ((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  return kTRUE;
}

//________________________________________________________________________
void AlidNdEtaTask::ConnectInputData(Option_t *opt)
{
  // Connect ESD
  // Called once
  
  AliAnalysisTaskSE::ConnectInputData(opt);

  Printf("AlidNdEtaTask::ConnectInputData called");

  /*
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();
    
    TString branches("AliESDHeader Vertex ");

    if (fAnalysisMode & AliPWG0Helper::kSPD || fTrigger & AliTriggerAnalysis::kOfflineFlag)
      branches += "AliMultiplicity ";
      
    if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS)
      branches += "Tracks ";
  
    // Enable only the needed branches
    esdH->SetActiveBranches(branches);
  }
  */

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
  
  #ifdef FULLALIROOT
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
    //AliCDBManager::Instance()->SetDefaultStorage("MC", "Residual");
    AliCDBManager::Instance()->SetRun(0);
    
    AliCDBManager* mgr = AliCDBManager::Instance();
    AliCDBEntry* obj = mgr->Get(AliCDBPath("GRP", "Geometry", "Data"));
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    
    AliGeomManager::GetNalignable("ITS");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  #endif
}

void AlidNdEtaTask::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AlidNdEtaTask::CreateOutputObjects");
  //AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);

  if (fOnlyPrimaries)
    Printf("WARNING: Processing only primaries (MC information used). This is for systematical checks only.");

  if (fUseMCKine)
    Printf("WARNING: Using MC kine information. This is for systematical checks only.");

  if (fUseMCVertex)
    Printf("WARNING: Replacing vertex by MC vertex. This is for systematical checks only.");

  fOutput = new TList;
  fOutput->SetOwner();

  fdNdEtaAnalysisESD = new dNdEtaAnalysis("fdNdEtaAnalysisESD", "fdNdEtaAnalysisESD", fAnalysisMode);
  fOutput->Add(fdNdEtaAnalysisESD);

  fMult = new TH1F("fMult", "fMult;Ntracks;Count", 201, -0.5, 200.5);
  fOutput->Add(fMult);

  fMultVtx = new TH1F("fMultVtx", "fMultVtx;Ntracks;Count", 201, -0.5, 200.5);
  fOutput->Add(fMultVtx);

  for (Int_t i=0; i<3; ++i)
  {
    fPartEta[i] = new TH1F(Form("dndeta_check_%d", i), Form("dndeta_check_%d", i), 60, -3, 3);
    fPartEta[i]->Sumw2();
    fOutput->Add(fPartEta[i]);
  }

  fEvents = new TH1F("dndeta_check_vertex", "dndeta_check_vertex", 800, -40, 40);
  fOutput->Add(fEvents);

  fVertexResolution = new TH2F("dndeta_vertex_resolution_z", ";z resolution;#Delta #phi", 1000, 0, 2, 100, 0, 0.2);
  fOutput->Add(fVertexResolution);

  fPhi = new TH1F("fPhi", "fPhi;#phi in rad.;count", 720, 0, 2 * TMath::Pi());
  fOutput->Add(fPhi);

  fEtaPhi = new TH2F("fEtaPhi", "fEtaPhi;#eta;#phi in rad.;count", 80, -4, 4, 18*20, 0, 2 * TMath::Pi());
  fOutput->Add(fEtaPhi);
  
  fStats = new TH1F("fStats", "fStats", 5, 0.5, 5.5);
  fStats->GetXaxis()->SetBinLabel(1, "vertexer 3d");
  fStats->GetXaxis()->SetBinLabel(2, "vertexer z");
  fStats->GetXaxis()->SetBinLabel(3, "trigger");
  fStats->GetXaxis()->SetBinLabel(4, "physics events");
  fStats->GetXaxis()->SetBinLabel(5, "physics events after veto");
  fOutput->Add(fStats);
  
  fStats2 = new TH2F("fStats2", "fStats2", 7, -0.5, 6.5, 12, -0.5, 11.5);
  
  fStats2->GetXaxis()->SetBinLabel(1, "No trigger");
  fStats2->GetXaxis()->SetBinLabel(2, "No Vertex");
  fStats2->GetXaxis()->SetBinLabel(3, "Vtx res. too bad");
  fStats2->GetXaxis()->SetBinLabel(4, "|z-vtx| > 15");
  fStats2->GetXaxis()->SetBinLabel(5, "|z-vtx| > 10");
  fStats2->GetXaxis()->SetBinLabel(6, "0 tracklets");
  fStats2->GetXaxis()->SetBinLabel(7, "Selected");
  
  fStats2->GetYaxis()->SetBinLabel(1, "n/a");
  fStats2->GetYaxis()->SetBinLabel(2, "empty");
  fStats2->GetYaxis()->SetBinLabel(3, "BBA");
  fStats2->GetYaxis()->SetBinLabel(4, "BBC");
  fStats2->GetYaxis()->SetBinLabel(5, "BBA BBC");
  fStats2->GetYaxis()->SetBinLabel(6, "BGA");
  fStats2->GetYaxis()->SetBinLabel(7, "BGC");
  fStats2->GetYaxis()->SetBinLabel(8, "BGA BGC");
  fStats2->GetYaxis()->SetBinLabel(9, "BBA BGC");
  fStats2->GetYaxis()->SetBinLabel(10, "BGA BBC");
  fStats2->GetYaxis()->SetBinLabel(11, "GFO");
  fStats2->GetYaxis()->SetBinLabel(12, "NO GFO");
  fOutput->Add(fStats2);

  fTrackletsVsClusters = new TH2F("fTrackletsVsClusters", ";tracklets;clusters in ITS", 50, -0.5, 49.5, 1000, -0.5, 999.5);
  fOutput->Add(fTrackletsVsClusters);
  
  if (fAnalysisMode & AliPWG0Helper::kSPD)
  {
    fDeltaPhi = new TH1F("fDeltaPhi", "fDeltaPhi;#Delta #phi;Entries", 500, -0.16, 0.16);
    fOutput->Add(fDeltaPhi);
    fDeltaTheta = new TH1F("fDeltaTheta", "fDeltaTheta;#Delta #theta;Entries", 500, -0.05, 0.05);
    fOutput->Add(fDeltaTheta);
    fFiredChips = new TH2F("fFiredChips", "fFiredChips;Chips L1 + L2;tracklets", 1201, -0.5, 1201, 50, -0.5, 49.5);
    fOutput->Add(fFiredChips);
    fTrackletsVsUnassigned = new TH2F("fTrackletsVsUnassigned", ";tracklets;unassigned clusters in L0", 50, -0.5, 49.5, 200, -0.5, 199.5);
    fOutput->Add(fTrackletsVsUnassigned);
    for (Int_t i=0; i<2; i++)
    {
      fZPhi[i] = new TH2F(Form("fZPhi_%d", i), Form("fZPhi Layer %d;z (cm);#phi (rad.)", i), 200, -15, 15, 200, 0, TMath::Pi() * 2);
      fOutput->Add(fZPhi[i]);
    }
    
    fModuleMap = new TH1F("fModuleMap", "fModuleMap;module number;cluster count", 240, -0.5, 239.5);
    fOutput->Add(fModuleMap);
  }

  if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS)
  {
    fRawPt =  new TH1F("fRawPt", "raw pt;p_{T};Count", 2000, 0, 100);
    fOutput->Add(fRawPt);
  }

  fVertex = new TH3F("vertex_check", "vertex_check;x;y;z", 100, -1, 1, 100, -1, 1, 100, -30, 30);
  fOutput->Add(fVertex);
  
  fVertexVsMult = new TH3F("fVertexVsMult", "fVertexVsMult;x;y;multiplicity", 100, -1, 1, 100, -1, 1, 100, -0.5, 99.5);
  fOutput->Add(fVertexVsMult);

  if (fReadMC)
  {
    fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysis);

    fdNdEtaAnalysisND = new dNdEtaAnalysis("dndetaND", "dndetaND", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisND);

    fdNdEtaAnalysisNSD = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisNSD);

    fdNdEtaAnalysisOnePart = new dNdEtaAnalysis("dndetaOnePart", "dndetaOnePart", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisOnePart);

    fdNdEtaAnalysisTr = new dNdEtaAnalysis("dndetaTr", "dndetaTr", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTr);

    fdNdEtaAnalysisTrVtx = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTrVtx);

    fdNdEtaAnalysisTracks = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisTracks);

    fPartPt =  new TH1F("dndeta_check_pt", "dndeta_check_pt", 1000, 0, 10);
    fPartPt->Sumw2();
    fOutput->Add(fPartPt);
  }

  if (fEsdTrackCuts)
  {
    fEsdTrackCuts->SetName("fEsdTrackCuts");
    fOutput->Add(fEsdTrackCuts);
  }

  fEta = new TH1D("fEta", "Eta;#eta;count", 80, -4, 4);
  fOutput->Add(fEta);
  
  fEtaMC = new TH1D("fEtaMC", "Eta, MC;#eta;count", 80, -4, 4);
  fOutput->Add(fEtaMC);

  fHistEvents = new TH1D("fHistEvents", "N. of Events;accepted;count", 2, 0, 2);
  fOutput->Add(fHistEvents);
  
  fHistEventsMC = new TH1D("fHistEventsMC", "N. of MC Events;accepted;count", 2, 0, 2);
  fOutput->Add(fHistEventsMC);

  fTrigEffNum = new TH1D("fTrigEffNum", "N. of triggered events", 100,0,100); 
  fOutput->Add(fTrigEffNum);
  fTrigEffDen = new TH1D("fTrigEffDen", "N. of true events", 100,0,100); 
  fOutput->Add(fTrigEffDen);
  fVtxEffNum = new TH1D("fVtxEffNum", "N. of events with vtx", 100,0,100); 
  fOutput->Add(fVtxEffNum);
  fVtxEffDen = new TH1D("fVtxEffDen", "N. of true events", 100,0,100); 
  fOutput->Add(fVtxEffDen);
  fVtxTrigEffNum = new TH1D("fVtxTrigEffNum", "N. of triggered events with vtx", 100,0,100); 
  fOutput->Add(fVtxTrigEffNum);
  fVtxTrigEffDen = new TH1D("fVtxTrigEffDen", "N. of triggered true events", 100,0,100); 
  fOutput->Add(fVtxTrigEffDen);
  
  //AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);
}

Bool_t AlidNdEtaTask::IsEventInBinZero()
{
  // checks if the event goes to the 0 bin
  
  Bool_t isZeroBin = kTRUE;
  fESD = (AliESDEvent*) fInputEvent;
  
  AliInputEventHandler* inputHandler = static_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return kFALSE;
  }
    
  static AliTriggerAnalysis* triggerAnalysis = 0;
  if (!triggerAnalysis)
  {
    AliPhysicsSelection* physicsSelection = dynamic_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if (physicsSelection)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
  }
  
  if (!triggerAnalysis)
  {
    Printf("ERROR: Could not receive trigger analysis object");
    return kFALSE;
  }
  
  if (!triggerAnalysis->IsTriggerFired(fESD, fTrigger))
    return kFALSE;
  
  // 0 bin check - from Michele
  
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  if (!mult){
    Printf("AlidNdEtaTask::IsBinZero: Can't get multiplicity object from ESDs");
    return kFALSE;
  }
  Int_t ntracklet = mult->GetNumberOfTracklets();
  const AliESDVertex * vtxESD = fESD->GetPrimaryVertexSPD();
  if(vtxESD) {
	  // If there is a vertex from vertexer z with delta phi > 0.02 we
	  // don't consider it rec (we keep the event in bin0). If quality
	  // is good eneough we check the number of tracklets
	  // if the vertex is more than 15 cm away, this is autamatically bin0
	  if( TMath::Abs(vtxESD->GetZ()) <= 15 ) {
		  if (vtxESD->IsFromVertexerZ()) {
			  if (vtxESD->GetDispersion()<=0.02 ) {
				  if(ntracklet>0) isZeroBin = kFALSE;
			  }
		  } else if(ntracklet>0) isZeroBin = kFALSE; // if the event is not from Vz we chek the n of tracklets
	  } 
  }
  return isZeroBin;
}

void AlidNdEtaTask::UserExec(Option_t*)
{
  // process the event

  fESD = (AliESDEvent*) fInputEvent;
  
  // these variables are also used in the MC section below; however, if ESD is off they stay with their default values
  Bool_t eventTriggered = kFALSE;
  const AliESDVertex* vtxESD = 0;

  // post the data already here
  PostData(1, fOutput);

  // ESD analysis
  if (fESD)
  {
    AliESDHeader* esdHeader = fESD->GetHeader();
    if (!esdHeader)
    {
      Printf("ERROR: esdHeader could not be retrieved");
      return;
    }
    
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    if (!inputHandler)
    {
      Printf("ERROR: Could not receive input handler");
      return;
    }
    
    // TODO use flags here!
    eventTriggered = inputHandler->IsEventSelected();
        
    static AliTriggerAnalysis* triggerAnalysis = 0;
    if (!triggerAnalysis)
    {
      AliPhysicsSelection* physicsSelection = dynamic_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
      if (physicsSelection)
        triggerAnalysis = physicsSelection->GetTriggerAnalysis();
    }
      
    if (eventTriggered)
      eventTriggered = triggerAnalysis->IsTriggerFired(fESD, fTrigger);
    
    AliTriggerAnalysis::V0Decision v0A = triggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kFALSE);
    AliTriggerAnalysis::V0Decision v0C = triggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kFALSE);
    Bool_t fastORHW = (triggerAnalysis->SPDFiredChips(fESD, 1) > 0);
    
    Int_t vZero = 0;
    if (v0A != AliTriggerAnalysis::kV0Invalid && v0C != AliTriggerAnalysis::kV0Invalid)
    {
      if (v0A == AliTriggerAnalysis::kV0Empty && v0C == AliTriggerAnalysis::kV0Empty) vZero = 1;
      if (v0A == AliTriggerAnalysis::kV0BB    && v0C == AliTriggerAnalysis::kV0Empty) vZero = 2;
      if (v0A == AliTriggerAnalysis::kV0Empty && v0C == AliTriggerAnalysis::kV0BB)    vZero = 3;
      if (v0A == AliTriggerAnalysis::kV0BB    && v0C == AliTriggerAnalysis::kV0BB)    vZero = 4;
      if (v0A == AliTriggerAnalysis::kV0BG    && v0C == AliTriggerAnalysis::kV0Empty) vZero = 5;
      if (v0A == AliTriggerAnalysis::kV0Empty && v0C == AliTriggerAnalysis::kV0BG)    vZero = 6;
      if (v0A == AliTriggerAnalysis::kV0BG    && v0C == AliTriggerAnalysis::kV0BG)    vZero = 7;
      if (v0A == AliTriggerAnalysis::kV0BB    && v0C == AliTriggerAnalysis::kV0BG)    vZero = 8;
      if (v0A == AliTriggerAnalysis::kV0BG    && v0C == AliTriggerAnalysis::kV0BB)    vZero = 9;
    }

    if (vZero == 0)
      Printf("VZERO: %d %d %d %d", vZero, eventTriggered, v0A, v0C);
      
    Bool_t filled = kFALSE;
      
    if (!eventTriggered)
    {
      fStats2->Fill(0.0, vZero);
      fStats2->Fill(0.0, (fastORHW) ? 10 : 11);
      filled = kTRUE;
    }
    
    if (eventTriggered)
      fStats->Fill(3);
      
    /*
    // ITS cluster tree
    AliESDInputHandlerRP* handlerRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalys
isManager()->GetInputEventHandler());
    if (handlerRP)
    {
      TTree* itsClusterTree = handlerRP->GetTreeR("ITS");
      if (!itsClusterTree)
        return;

      TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
      TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

      itsClusterBranch->SetAddress(&itsClusters);

      Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();

      Int_t totalClusters = 0;

      // loop over the its subdetectors
      for (Int_t iIts=0; iIts < nItsSubs; iIts++) {

        if (!itsClusterTree->GetEvent(iIts))
          continue;

        Int_t nClusters = itsClusters->GetEntriesFast();
        totalClusters += nClusters;

        #ifdef FULLALIROOT
          if (fAnalysisMode & AliPWG0Helper::kSPD)
          {
            // loop over clusters
            while (nClusters--) {
              AliITSRecPoint* cluster = (AliITSRecPoint*) itsClusters->UncheckedAt(nClusters);

              Int_t layer = cluster->GetLayer();

              if (layer > 1)
                continue;

              Float_t xyz[3] = {0., 0., 0.};
              cluster->GetGlobalXYZ(xyz);

              Float_t phi = TMath::Pi() + TMath::ATan2(-xyz[1], -xyz[0]);
              Float_t z = xyz[2];

              fZPhi[layer]->Fill(z, phi);
              fModuleMap->Fill(layer * 80 + cluster->GetDetectorIndex());
            }
          }
        #endif
      }
    }
    */
      
    vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
    if (!vtxESD)
    {
      if (!filled)
      {
        fStats2->Fill(1, vZero);
        fStats2->Fill(1, (fastORHW) ? 10 : 11);
        filled = kTRUE;
      }
    }
    else
    {
      if (!AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
      {
        if (!filled)
        {
          fStats2->Fill(2, vZero);
          fStats2->Fill(2, (fastORHW) ? 10 : 11);
          filled = kTRUE;
        }
      }
      
      Double_t vtx[3];
      vtxESD->GetXYZ(vtx);
      
      // try to compare spd vertex and vertexer from tracks
      // remove vertices outside +- 15 cm
      if (TMath::Abs(vtx[2]) > 15)
      {
        if (!filled)
        {
          fStats2->Fill(3, vZero);
          fStats2->Fill(3, (fastORHW) ? 10 : 11);
          filled = kTRUE;
        }
      }
      
      if (TMath::Abs(vtx[2]) > 10)
      {
        if (!filled)
        {
          fStats2->Fill(4, vZero);
          fStats2->Fill(4, (fastORHW) ? 10 : 11);
          filled = kTRUE;
        }
      }
        
      const AliMultiplicity* mult = fESD->GetMultiplicity();
      if (!mult)
      {
	Printf("Returning, no Multiplicity found");
        return;
      }
      
      if (mult->GetNumberOfTracklets() == 0)
      {
        if (!filled)
        {
          fStats2->Fill(5, vZero);
          fStats2->Fill(5, (fastORHW) ? 10 : 11);
          filled = kTRUE;
        }
      }
    }
    
    if (!filled)
    {
      fStats2->Fill(6, vZero);
      fStats2->Fill(6, (fastORHW) ? 10 : 11);
      //Printf("File: %s, IEV: %d, TRG: ---, Orbit: 0x%x, Period: %d, BC: %d", ((TTree*) GetInputData(0))->GetCurrentFile()->GetName(), fESD->GetEventNumberInFile(), fESD->GetOrbitNumber(),fESD->GetPeriodNumber(),fESD->GetBunchCrossNumber());
    }
      
    if (eventTriggered)
      fStats->Fill(3);
    
    fStats->Fill(5);
    
    // get the ESD vertex
    vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
    
    Double_t vtx[3];

    // fill z vertex resolution
    if (vtxESD)
    {
      if (strcmp(vtxESD->GetTitle(), "vertexer: Z") == 0)
        fVertexResolution->Fill(vtxESD->GetZRes(), vtxESD->GetDispersion());
      else
        fVertexResolution->Fill(vtxESD->GetZRes(), 0);
      
      //if (strcmp(vtxESD->GetTitle(), "vertexer: 3D") == 0)
      {
        fVertex->Fill(vtxESD->GetXv(), vtxESD->GetYv(), vtxESD->GetZv());
      }
      
      if (AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
      {
          vtxESD->GetXYZ(vtx);

          // vertex stats
          if (strcmp(vtxESD->GetTitle(), "vertexer: 3D") == 0)
            fStats->Fill(1);
          
          if (strcmp(vtxESD->GetTitle(), "vertexer: Z") == 0)
            fStats->Fill(2);
          
          // remove vertices outside +- 15 cm
          if (TMath::Abs(vtx[2]) > 15)
            vtxESD = 0;
      }
      else
        vtxESD = 0;
    }

    // needed for syst. studies
    AliStack* stack = 0;
    TArrayF vtxMC(3);

    if (fUseMCVertex || fUseMCKine || fOnlyPrimaries || fReadMC) {
      if (!fReadMC) {
        Printf("ERROR: fUseMCVertex or fUseMCKine or fOnlyPrimaries set without fReadMC set!");
        return;
      }

      AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
        Printf("ERROR: Could not retrieve MC event handler");
        return;
      }

      AliMCEvent* mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
        Printf("ERROR: Could not retrieve MC event");
        return;
      }

      AliHeader* header = mcEvent->Header();
      if (!header)
      {
        AliDebug(AliLog::kError, "Header not available");
        return;
      }

      // get the MC vertex
      AliGenEventHeader* genHeader = header->GenEventHeader();
      if (!genHeader)
      {
        AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
        return;
      }
      genHeader->PrimaryVertex(vtxMC);

      if (fUseMCVertex)
        vtx[2] = vtxMC[2];

      stack = mcEvent->Stack();
      if (!stack)
      {
        AliDebug(AliLog::kError, "Stack not available");
        return;
      }
    }
    
    // create list of (label, eta, pt) tuples
    Int_t inputCount = 0;
    Int_t* labelArr = 0;
    Float_t* etaArr = 0;
    Float_t* thirdDimArr = 0;
    if (fAnalysisMode & AliPWG0Helper::kSPD)
    {
      if (vtxESD)
      {
        // get tracklets
        const AliMultiplicity* mult = fESD->GetMultiplicity();
        if (!mult)
        {
          AliDebug(AliLog::kError, "AliMultiplicity not available");
          return;
        }
  
        Int_t arrayLength = mult->GetNumberOfTracklets();
        if (fAnalysisMode & AliPWG0Helper::kSPDOnlyL0)
          arrayLength += mult->GetNumberOfSingleClusters();
          
        labelArr = new Int_t[arrayLength];
        etaArr = new Float_t[arrayLength];
        thirdDimArr = new Float_t[arrayLength];
  
        // get multiplicity from SPD tracklets
        for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
        {
          //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));
  
          if (fOnlyPrimaries)
            if (mult->GetLabel(i, 0) < 0 || mult->GetLabel(i, 0) != mult->GetLabel(i, 1) || !stack->IsPhysicalPrimary(mult->GetLabel(i, 0)))
              continue;
  
          Float_t deltaPhi = mult->GetDeltaPhi(i);
          // prevent values to be shifted by 2 Pi()
          if (deltaPhi < -TMath::Pi())
            deltaPhi += TMath::Pi() * 2;
          if (deltaPhi > TMath::Pi())
            deltaPhi -= TMath::Pi() * 2;
  
          if (TMath::Abs(deltaPhi) > 1)
            printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);
  
          Int_t label = mult->GetLabel(i, 0);
          Float_t eta = mult->GetEta(i);
          
          // control histograms
          Float_t phi = mult->GetPhi(i);
          if (phi < 0)
            phi += TMath::Pi() * 2;
            
          // TEST exclude probably inefficient phi region
          //if (phi > 5.70 || phi < 0.06)
          //  continue;
            
          fPhi->Fill(phi);
          
          if (vtxESD && TMath::Abs(vtx[2]) < 10)
          {
            fEtaPhi->Fill(eta, phi);
            fDeltaPhi->Fill(deltaPhi);
            fDeltaTheta->Fill(mult->GetDeltaTheta(i));
          }
          
          if (fDeltaPhiCut > 0 && (TMath::Abs(deltaPhi) > fDeltaPhiCut || TMath::Abs(mult->GetDeltaTheta(i)) > fDeltaPhiCut / 0.08 * 0.025))
            continue;
  
          if (fUseMCKine)
          {
            if (label > 0)
            {
              TParticle* particle = stack->Particle(label);
              eta = particle->Eta();
              phi = particle->Phi();
            }
            else
              Printf("WARNING: fUseMCKine set without fOnlyPrimaries and no label found");
          }
          
          if (fSymmetrize)
            eta = TMath::Abs(eta);
  
          etaArr[inputCount] = eta;
          labelArr[inputCount] = label;
          thirdDimArr[inputCount] = phi;
          ++inputCount;
        }
        
        if (fAnalysisMode & AliPWG0Helper::kSPDOnlyL0)
        {
          // get additional clusters from L0 
          for (Int_t i=0; i<mult->GetNumberOfSingleClusters(); ++i)
          {
            etaArr[inputCount] = -TMath::Log(TMath::Tan(mult->GetThetaSingle(i)/2.));
            labelArr[inputCount] = -1;
            thirdDimArr[inputCount] = mult->GetPhiSingle(i);
            
            ++inputCount;
          }
        }
        
        if (!fFillPhi)
        {
          // fill multiplicity in third axis
          for (Int_t i=0; i<inputCount; ++i)
            thirdDimArr[i] = inputCount;
        }
  
        Int_t firedChips = mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
        fFiredChips->Fill(firedChips, inputCount);
        Printf("Accepted %d tracklets (%d fired chips)", inputCount, firedChips);
        
        fTrackletsVsUnassigned->Fill(inputCount, mult->GetNumberOfSingleClusters());
      }
    }
    else if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS)
    {
      if (!fEsdTrackCuts)
      {
        AliDebug(AliLog::kError, "fESDTrackCuts not available");
        return;
      }

      Bool_t foundInEta10 = kFALSE;
      
      if (vtxESD)
      {
        // get multiplicity from ESD tracks
        TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD, fAnalysisMode & AliPWG0Helper::kTPC);
        Int_t nGoodTracks = list->GetEntries();
        Printf("Accepted %d tracks out of %d total ESD tracks", nGoodTracks, fESD->GetNumberOfTracks());
  
        labelArr = new Int_t[nGoodTracks];
        etaArr = new Float_t[nGoodTracks];
        thirdDimArr = new Float_t[nGoodTracks];
  
        // loop over esd tracks
        for (Int_t i=0; i<nGoodTracks; i++)
        {
          AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
          if (!esdTrack)
          {
            AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", i));
            continue;
          }
          
          Float_t phi = esdTrack->Phi();
          if (phi < 0)
            phi += TMath::Pi() * 2;
  
          Float_t eta = esdTrack->Eta();
          Int_t label = TMath::Abs(esdTrack->GetLabel());
          Float_t pT  = esdTrack->Pt();
          
          // force pT to fixed value without B field
          if (!(fAnalysisMode & AliPWG0Helper::kFieldOn))
            pT = 1;
  
          fPhi->Fill(phi);
          fEtaPhi->Fill(eta, phi);
          if (eventTriggered && vtxESD)
            fRawPt->Fill(pT);
  
          if (esdTrack->Pt() < fPtMin) // even if the pt cut is already in the esdTrackCuts....
            continue;
          
          if (fOnlyPrimaries)
          {
	   if (label == 0)
	     continue;
	   
	   if (stack->IsPhysicalPrimary(label) == kFALSE)
	     continue;
          }

	  // 2 types of INEL>0 trigger - choose one

	  // 1. HL 
          // INEL>0 trigger
          // if (TMath::Abs(esdTrack->Eta()) < 1 && esdTrack->Pt() > 0.15)
          // foundInEta10 = kTRUE;

          // 2. MB working group
	  if (TMath::Abs(esdTrack->Eta()) < 0.8 && esdTrack->Pt() > fPtMin){  // this should be in the trigger selection as well, so one particle should always be found
	    foundInEta10 = kTRUE;
	  }
  
          if (fUseMCKine)
          {
	    if (label > 0)
            {
              TParticle* particle = stack->Particle(label);
              eta = particle->Eta();
              pT = particle->Pt();
	      // check when using INEL>0, MB Working Group definition
	      if (TMath::Abs(eta) >=0.8 || pT <= fPtMin){
		AliDebug(2,Form("WARNING ************* USING KINE: eta = %f, pt = %f",eta,pT));
	      }
            }
            else
              Printf("WARNING: fUseMCKine set without fOnlyPrimaries and no label found");
	  }
  
          if (fSymmetrize)
            eta = TMath::Abs(eta);
          etaArr[inputCount] = eta;
          labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
          thirdDimArr[inputCount] = pT;
          ++inputCount;
        }
        
        //if (inputCount > 30)
        //  Printf("Event with %d accepted TPC tracks. Period number: %d Orbit number: %x Bunch crossing number: %d", inputCount, fESD->GetPeriodNumber(), fESD->GetOrbitNumber(), fESD->GetBunchCrossNumber());
        
        // TODO restrict inputCount used as measure for the multiplicity to |eta| < 1
  
        delete list;
      }
      
      if (!foundInEta10){
        eventTriggered = kFALSE;
	fHistEvents->Fill(0.1);
	AliDebug(3,"Event rejected");
      }
      else{
	if (eventTriggered){
	  fHistEvents->Fill(1.1);
	  AliDebug(3,Form("Event Accepted, with inputcount = %d",inputCount));
	}
	else{
	  AliDebug(3,"Event has foundInEta10 but was not triggered");
	}
      }
    }
    else
      return;

    // Processing of ESD information (always)
    if (eventTriggered)
    {
      // control hist
      fMult->Fill(inputCount);
      fdNdEtaAnalysisESD->FillTriggeredEvent(inputCount);

      if (vtxESD)
      {
        // control hist
        
        if (strcmp(vtxESD->GetTitle(), "vertexer: 3D") == 0)
          fVertexVsMult->Fill(vtxESD->GetXv(), vtxESD->GetYv(), inputCount);
      
        Int_t part05 = 0;
        Int_t part10 = 0;
        for (Int_t i=0; i<inputCount; ++i)
        {
          Float_t eta = etaArr[i];
          Float_t thirdDim = thirdDimArr[i];
	  fEta->Fill(eta);

          fdNdEtaAnalysisESD->FillTrack(vtx[2], eta, thirdDim);

          if (TMath::Abs(vtx[2]) < 10)
          {
            fPartEta[0]->Fill(eta);

            if (vtx[2] < 0)
              fPartEta[1]->Fill(eta);
            else
              fPartEta[2]->Fill(eta);
          }
          
          if (TMath::Abs(eta) < 0.5)
            part05++;
          if (TMath::Abs(eta) < 1.0) // in the INEL>0, MB WG definition, this is in any case equivalent to <0.8, due to the EsdTrackCuts
            part10++;
        }
        
        Int_t multAxis = inputCount;
        if (fMultAxisEta1)
          multAxis = part10;

        fMultVtx->Fill(multAxis);
        //if (TMath::Abs(vtx[2]) < 10)
        //  fMultVtx->Fill(part05);

        // for event count per vertex
        fdNdEtaAnalysisESD->FillEvent(vtx[2], multAxis);

        // control hist
        if (multAxis > 0)
          fEvents->Fill(vtx[2]);

        if (fReadMC)
        {
          // from tracks is only done for triggered and vertex reconstructed events
          for (Int_t i=0; i<inputCount; ++i)
          {
            Int_t label = labelArr[i];

            if (label < 0)
              continue;

            //Printf("Getting particle of track %d", label);
            TParticle* particle = stack->Particle(label);

            if (!particle)
            {
              AliDebug(AliLog::kError, Form("ERROR: Could not retrieve particle %d.", label));
              continue;
            }

            Float_t thirdDim = -1;
            if (fAnalysisMode & AliPWG0Helper::kSPD)
            {
              if (fFillPhi)
              {
                thirdDim = particle->Phi();
              }
              else
                thirdDim = multAxis;
            }
            else
              thirdDim = particle->Pt();

            Float_t eta = particle->Eta();
            if (fSymmetrize)
              eta = TMath::Abs(eta);
            fdNdEtaAnalysisTracks->FillTrack(vtxMC[2], eta, thirdDim);
          } // end of track loop

          // for event count per vertex
          fdNdEtaAnalysisTracks->FillEvent(vtxMC[2], multAxis);
        }
      }
    }
    if (etaArr)
      delete[] etaArr;
    if (labelArr)
      delete[] labelArr;
    if (thirdDimArr)
      delete[] thirdDimArr;
  }

  if (fReadMC)   // Processing of MC information (optional)
  {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    AliStack* stack = mcEvent->Stack();
    if (!stack)
    {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    AliHeader* header = mcEvent->Header();
    if (!header)
    {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    // get the MC vertex
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if (!genHeader)
    {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }

    TArrayF vtxMC(3);
    genHeader->PrimaryVertex(vtxMC);
    
    // get process type
    Int_t processType = AliPWG0Helper::GetEventProcessType(fESD, header, stack, fDiffTreatment);
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", processType));

    if (processType == AliPWG0Helper::kInvalidProcess)
      AliDebug(AliLog::kError, Form("Unknown process type %d.", processType));

    // loop over mc particles
    Int_t nPrim  = stack->GetNprimary();

    Int_t nAcceptedParticles = 0;
    Bool_t oneParticleEvent = kFALSE;

    Int_t nMCparticlesinRange = 0; // number of true particles in my range of interest

    // count particles first, then fill
    for (Int_t iMc = 0; iMc < nPrim; ++iMc)
    {
      if (!stack->IsPhysicalPrimary(iMc))
        continue;

      //Printf("Getting particle %d", iMc);
      TParticle* particle = stack->Particle(iMc);

      if (!particle)
        continue;

      if (particle->GetPDG()->Charge() == 0)
        continue;

      //AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", iMc, particle->GetUniqueID()));
      Float_t eta = particle->Eta();

      AliDebug(2,Form("particle %d: eta = %f, pt = %f",iMc,particle->Eta(),particle->Pt()));

      // INEL>0: choose one
      // 1.HL
      // if (TMath::Abs(eta) < 1.0){
      //   oneParticleEvent = kTRUE;
      //   nMCparticlesinRange++;
      //}
      // 2.MB Working Group definition
      if (TMath::Abs(eta) < 0.8 && particle->Pt() > fPtMin){
	oneParticleEvent = kTRUE;
	nMCparticlesinRange++;
      }

      // make a rough eta cut (so that nAcceptedParticles is not too far off the true measured value (NB: these histograms are only gathered for comparison))
      if (TMath::Abs(eta) < 1.5) // && pt > 0.3)
        nAcceptedParticles++;
    }

    if (TMath::Abs(vtxMC[2]) < 10.){
      // if MC vertex is inside 10
      if (vtxESD){
	fVtxEffNum->Fill(nMCparticlesinRange);
      }
      fVtxEffDen->Fill(nMCparticlesinRange);
      if (eventTriggered){
	if (vtxESD){
	   fVtxTrigEffNum->Fill(nMCparticlesinRange);
	}
	fVtxTrigEffDen->Fill(nMCparticlesinRange);
	fTrigEffNum->Fill(nMCparticlesinRange);
      }
      fTrigEffDen->Fill(nMCparticlesinRange);
    }

    if (oneParticleEvent){
	    fHistEventsMC->Fill(1.1);
    }
    else{
	    fHistEventsMC->Fill(0.1);
    }	

    for (Int_t iMc = 0; iMc < nPrim; ++iMc)
    {
      if (!stack->IsPhysicalPrimary(iMc))
        continue;

      //Printf("Getting particle %d", iMc);
      TParticle* particle = stack->Particle(iMc);

      if (!particle)
        continue;

      if (particle->GetPDG()->Charge() == 0)
        continue;

      Float_t eta = particle->Eta();
      if (fSymmetrize)
        eta = TMath::Abs(eta);

      Float_t thirdDim = -1;

      if (fAnalysisMode & AliPWG0Helper::kSPD)
      {
        if (fFillPhi)
        {
          thirdDim = particle->Phi();
        }
        else
          thirdDim = nAcceptedParticles;
      }
      else
        thirdDim = particle->Pt();

      fdNdEtaAnalysis->FillTrack(vtxMC[2], eta, thirdDim);

      if (processType != AliPWG0Helper::kSD)
        fdNdEtaAnalysisNSD->FillTrack(vtxMC[2], eta, thirdDim);

      if (processType == AliPWG0Helper::kND)
        fdNdEtaAnalysisND->FillTrack(vtxMC[2], eta, thirdDim);
      
      if (oneParticleEvent){
	      AliDebug(3,Form("filling dNdEtaAnalysis object:: vtx = %f, eta = %f, pt = %f",vtxMC[2], eta, thirdDim));
	      fdNdEtaAnalysisOnePart->FillTrack(vtxMC[2], eta, thirdDim);
      }

      if (eventTriggered)
      {
        fdNdEtaAnalysisTr->FillTrack(vtxMC[2], eta, thirdDim);
        if (vtxESD)
          fdNdEtaAnalysisTrVtx->FillTrack(vtxMC[2], eta, thirdDim);
      }

      if (TMath::Abs(eta) < 1.0 && particle->Pt() > 0 && particle->P() > 0)
      {
        //Float_t value = 1. / TMath::TwoPi() / particle->Pt() * particle->Energy() / particle->P();
        Float_t value = 1;
        fPartPt->Fill(particle->Pt(), value);
	if (TMath::Abs(eta) < 0.8 && particle->Pt() > fPtMin){
	  fEtaMC->Fill(eta);
	}
      }
    }

    fdNdEtaAnalysis->FillEvent(vtxMC[2], nAcceptedParticles);
    if (processType != AliPWG0Helper::kSD)
      fdNdEtaAnalysisNSD->FillEvent(vtxMC[2], nAcceptedParticles);
    if (processType == AliPWG0Helper::kND)
      fdNdEtaAnalysisND->FillEvent(vtxMC[2], nAcceptedParticles);
    if (oneParticleEvent)
      fdNdEtaAnalysisOnePart->FillEvent(vtxMC[2], nAcceptedParticles);

    if (eventTriggered)
    {
      fdNdEtaAnalysisTr->FillEvent(vtxMC[2], nAcceptedParticles);
      if (vtxESD)
        fdNdEtaAnalysisTrVtx->FillEvent(vtxMC[2], nAcceptedParticles);
    }
  }
}

void AlidNdEtaTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput)
    Printf("ERROR: fOutput not available");

  TH1D* dNdEta = new TH1D("dNdEta","#eta;counts",80, -4, 4);
  TH1D* dNdEtaMC = new TH1D("dNdEtaMC","#eta,MC;counts",80, -4, 4);
  if (fOutput)
  {
    fdNdEtaAnalysisESD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("fdNdEtaAnalysisESD"));
    fMult = dynamic_cast<TH1F*> (fOutput->FindObject("fMult"));
    fMultVtx = dynamic_cast<TH1F*> (fOutput->FindObject("fMultVtx"));
    for (Int_t i=0; i<3; ++i)
      fPartEta[i] = dynamic_cast<TH1F*> (fOutput->FindObject(Form("dndeta_check_%d", i)));
    fEvents = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_vertex"));
    fVertexResolution = dynamic_cast<TH2F*> (fOutput->FindObject("dndeta_vertex_resolution_z"));

    fVertex = dynamic_cast<TH3F*> (fOutput->FindObject("vertex_check"));
    fVertexVsMult = dynamic_cast<TH3F*> (fOutput->FindObject("fVertexVsMult"));
    fPhi = dynamic_cast<TH1F*> (fOutput->FindObject("fPhi"));
    fRawPt = dynamic_cast<TH1F*> (fOutput->FindObject("fRawPt"));
    fEtaPhi = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaPhi"));
    for (Int_t i=0; i<2; ++i)
      fZPhi[i] = dynamic_cast<TH2F*> (fOutput->FindObject(Form("fZPhi_%d", i)));
    fModuleMap = dynamic_cast<TH1F*> (fOutput->FindObject("fModuleMap"));
    fDeltaPhi = dynamic_cast<TH1F*> (fOutput->FindObject("fDeltaPhi"));
    fDeltaTheta = dynamic_cast<TH1F*> (fOutput->FindObject("fDeltaTheta"));
    fFiredChips = dynamic_cast<TH2F*> (fOutput->FindObject("fFiredChips"));
    fTrackletsVsClusters = dynamic_cast<TH2F*> (fOutput->FindObject("fTrackletsVsClusters"));
    fTrackletsVsUnassigned = dynamic_cast<TH2F*> (fOutput->FindObject("fTrackletsVsUnassigned"));
    fStats = dynamic_cast<TH1F*> (fOutput->FindObject("fStats"));
    fStats2 = dynamic_cast<TH2F*> (fOutput->FindObject("fStats2"));

    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCuts"));

    fEta = dynamic_cast<TH1D*>(fOutput->FindObject("fEta"));
    fEtaMC = dynamic_cast<TH1D*>(fOutput->FindObject("fEtaMC"));
    fHistEvents = dynamic_cast<TH1D*>(fOutput->FindObject("fHistEvents"));
    fHistEventsMC = dynamic_cast<TH1D*>(fOutput->FindObject("fHistEventsMC"));
    fVtxEffDen = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxEffDen"));
    fVtxEffNum = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxEffNum"));
    fTrigEffDen = dynamic_cast<TH1D*>(fOutput->FindObject("fTrigEffDen"));
    fTrigEffNum = dynamic_cast<TH1D*>(fOutput->FindObject("fTrigEffNum"));
    fVtxTrigEffDen = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxTrigEffDen"));
    fVtxTrigEffNum = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxTrigEffNum"));
    if (!fHistEvents){
	    AliError("fHistEvents not found, impossible to determine the corresponding dNdEta distribution for the control file");
    }		    
    else Printf("Events selected = %f", fHistEvents->GetBinContent(fHistEvents->GetXaxis()->FindBin(1.1)));
    if (!fHistEventsMC){
	    AliError("fHistEventsMC not found, impossible to determine the corresponding dNdEta distribution for the control file");
    }	
    else Printf("Events selected frm MC = %f", fHistEventsMC->GetBinContent(fHistEventsMC->GetXaxis()->FindBin(1.1)));
    if (!fEta){
	    AliError("fEta not found, impossible to determine the corresponding dNdEta distribution for the control file");
    }	
    if (!fEtaMC){
	    AliError("fEtaMC not found, impossible to determine the corresponding dNdEta distribution for the control file");
    }		    
    Float_t nevents = 0;
    if (fHistEvents) nevents = fHistEvents->GetBinContent(fHistEvents->GetXaxis()->FindBin(1.1));
    Float_t neventsMC = 0;
    if (fHistEventsMC) neventsMC = fHistEventsMC->GetBinContent(fHistEventsMC->GetXaxis()->FindBin(1.1));
    if (fEta && fEtaMC){
	    for (Int_t ibin = 1; ibin <= fEta->GetNbinsX(); ibin++){
		    Float_t eta =0;
		    Float_t etaerr =0;
		    Float_t etaMC =0;
		    Float_t etaerrMC =0;
		    if (fEta && nevents > 0) {
			    eta = fEta->GetBinContent(ibin)/nevents/fEta->GetBinWidth(ibin);
			    etaerr = fEta->GetBinError(ibin)/nevents/fEta->GetBinWidth(ibin);
			    dNdEta->SetBinContent(ibin,eta);
			    dNdEta->SetBinError(ibin,etaerr);
		    }
		    if (fEtaMC && neventsMC > 0) {
			    etaMC = fEtaMC->GetBinContent(ibin)/neventsMC/fEtaMC->GetBinWidth(ibin);
			    etaerrMC = fEtaMC->GetBinError(ibin)/neventsMC/fEtaMC->GetBinWidth(ibin);
			    dNdEtaMC->SetBinContent(ibin,etaMC);
			    dNdEtaMC->SetBinError(ibin,etaerrMC);
		    }
	    }
    }
    new TCanvas("eta", " eta ",50, 50, 550, 550) ;
    if (fEta) fEta->Draw();
    new TCanvas("etaMC", " etaMC ",50, 50, 550, 550) ;
    if (fEtaMC) fEtaMC->Draw();
    new TCanvas("dNdEta", "#eta;dNdEta ",50, 50, 550, 550) ;
    dNdEta->Draw();
    new TCanvas("dNdEtaMC", "#eta,MC;dNdEta ",50, 50, 550, 550) ;
    dNdEtaMC->Draw();
    new TCanvas("Events", "Events;Events ",50, 50, 550, 550) ;
    if (fHistEvents) fHistEvents->Draw();
    new TCanvas("Events, MC", "Events, MC;Events ",50, 50, 550, 550) ;
    if (fHistEventsMC) fHistEventsMC->Draw();
    
    TFile* outputFileCheck = new TFile("histogramsCheck.root", "RECREATE");
    if (fEta) fEta->Write();
    if (fEtaMC) fEtaMC->Write();
    dNdEta->Write();
    dNdEtaMC->Write();
    outputFileCheck->Write();
    outputFileCheck->Close();

  }

  if (!fdNdEtaAnalysisESD)
  {
    AliDebug(AliLog::kError, "ERROR: fdNdEtaAnalysisESD not available");
  }
  else
  {
    if (fMult && fMultVtx)
    {
      new TCanvas;
      fMult->Draw();
      fMultVtx->SetLineColor(2);
      fMultVtx->Draw("SAME");
    }

    if (fFiredChips)
    {
      new TCanvas;
      fFiredChips->Draw("COLZ");
    }

    if (fPartEta[0] && fPartEta[1] && fPartEta[2] && fEvents && fMultVtx && fMult)
    {
      Int_t events1 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(-9.9), fEvents->GetXaxis()->FindBin(-0.001));
      Int_t events2 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(0.001), fEvents->GetXaxis()->FindBin(9.9));

      Printf("%d events with vertex used", events1 + events2);
      Printf("%d total events with vertex, %d total triggered events, %d triggered events with 0 multiplicity", (Int_t) fMultVtx->Integral(), (Int_t) fMult->Integral(), (Int_t) fMult->GetBinContent(1));

      if (events1 > 0 && events2 > 0)
      {
        fPartEta[0]->Scale(1.0 / (events1 + events2));
        fPartEta[1]->Scale(1.0 / events1);
        fPartEta[2]->Scale(1.0 / events2);

        for (Int_t i=0; i<3; ++i)
          fPartEta[i]->Scale(1.0 / fPartEta[i]->GetBinWidth(1));

        new TCanvas("control", "control", 500, 500);
        for (Int_t i=0; i<3; ++i)
        {
          fPartEta[i]->SetLineColor(i+1);
          fPartEta[i]->Draw((i==0) ? "" : "SAME");
        }
      }
    }

    if (fEvents)
    {
        new TCanvas("control3", "control3", 500, 500);
        fEvents->Draw();
    }
    
    if (fStats2)
    {
      new TCanvas;
      fStats2->Draw("TEXT");
    }

    TFile* fout = new TFile("analysis_esd_raw.root", "RECREATE");

    if (fdNdEtaAnalysisESD)
      fdNdEtaAnalysisESD->SaveHistograms();

    if (fEsdTrackCuts)
      fEsdTrackCuts->SaveHistograms("esd_track_cuts");

    if (fMult)
      fMult->Write();

    if (fMultVtx)
      fMultVtx->Write();

    for (Int_t i=0; i<3; ++i)
      if (fPartEta[i])
        fPartEta[i]->Write();

    if (fEvents)
      fEvents->Write();

    if (fVertexResolution)
    {
      fVertexResolution->Write();
      fVertexResolution->ProjectionX();
      fVertexResolution->ProjectionY();
    }

    if (fDeltaPhi)
      fDeltaPhi->Write();

    if (fDeltaTheta)
      fDeltaTheta->Write();
    
    if (fPhi)
      fPhi->Write();

    if (fRawPt)
      fRawPt->Write();

    if (fEtaPhi)
      fEtaPhi->Write();

    for (Int_t i=0; i<2; ++i)
      if (fZPhi[i])
        fZPhi[i]->Write();
    
    if (fModuleMap)
      fModuleMap->Write();
    
    if (fFiredChips)
      fFiredChips->Write();

    if (fTrackletsVsClusters)
      fTrackletsVsClusters->Write();
    
    if (fTrackletsVsUnassigned)
      fTrackletsVsUnassigned->Write();
    
    if (fStats)
      fStats->Write();

    if (fStats2)
      fStats2->Write();
    
    if (fVertex)
      fVertex->Write();

    if (fVertexVsMult)
      fVertexVsMult->Write();
    
    if (fVtxEffDen) fVtxEffDen->Write();
    else Printf("fVtxEffDen not found");

    if (fVtxEffNum) fVtxEffNum->Write();
    else Printf("fVtxEffNum not found");

    if (fVtxTrigEffNum) fVtxTrigEffNum->Write();
    else Printf("fVtxTrigEffNum not found");

    if (fVtxTrigEffDen) fVtxTrigEffDen->Write();
    else Printf("fVtxTrigEffDen not found");

    if (fTrigEffNum) fTrigEffNum->Write();
    else Printf("fTrigEffNum not found");

    if (fTrigEffDen) fTrigEffDen->Write();
    else Printf("fTrigEffDen not found");

    fout->Write();
    fout->Close();

    Printf("Writting result to analysis_esd_raw.root");
  }

  if (fReadMC)
  {
    if (fOutput)
    {
      fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
      fdNdEtaAnalysisND = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaND"));
      fdNdEtaAnalysisNSD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaNSD"));
      fdNdEtaAnalysisOnePart = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaOnePart"));
      fdNdEtaAnalysisTr = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTr"));
      fdNdEtaAnalysisTrVtx = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTrVtx"));
      fdNdEtaAnalysisTracks = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTracks"));
      fPartPt = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_pt"));
    }

    if (!fdNdEtaAnalysis || !fdNdEtaAnalysisTr || !fdNdEtaAnalysisTrVtx || !fPartPt)
    {
      AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) fdNdEtaAnalysis, (void*) fPartPt));
      return;
    }

    fdNdEtaAnalysis->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisND->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisNSD->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisOnePart->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTr->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTrVtx->Finish(0, -1, AlidNdEtaCorrection::kNone);
    fdNdEtaAnalysisTracks->Finish(0, -1, AlidNdEtaCorrection::kNone);

    Int_t events = (Int_t) fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram()->Integral();
    fPartPt->Scale(1.0/events);
    fPartPt->Scale(1.0/fPartPt->GetBinWidth(1));

    TFile* fout = new TFile("analysis_mc.root","RECREATE");

    fdNdEtaAnalysis->SaveHistograms();
    fdNdEtaAnalysisND->SaveHistograms();
    fdNdEtaAnalysisNSD->SaveHistograms();
    fdNdEtaAnalysisOnePart->SaveHistograms();
    fdNdEtaAnalysisTr->SaveHistograms();
    fdNdEtaAnalysisTrVtx->SaveHistograms();
    fdNdEtaAnalysisTracks->SaveHistograms();

    if (fPartPt)
      fPartPt->Write();

    fout->Write();
    fout->Close();

    Printf("Writting result to analysis_mc.root");

    if (fPartPt)
    {
      new TCanvas("control2", "control2", 500, 500);
      fPartPt->Draw();
    }
  }
}
