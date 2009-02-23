/* $Id$ */

#include "AlidNdEtaTask.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
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
#include <AliESDHeader.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "AliCorrection.h"
#include "AliCorrectionMatrix3D.h"
#include "dNdEta/dNdEtaAnalysis.h"

ClassImp(AlidNdEtaTask)

AlidNdEtaTask::AlidNdEtaTask(const char* opt) :
  AliAnalysisTask("AlidNdEtaTask", ""),
  fESD(0),
  fOutput(0),
  fOption(opt),
  fAnalysisMode(AliPWG0Helper::kTPC),
  fTrigger(AliPWG0Helper::kMB1),
  fFillPhi(kFALSE),
  fDeltaPhiCut(-1),
  fReadMC(kFALSE),
  fUseMCVertex(kFALSE),
  fOnlyPrimaries(kFALSE),
  fUseMCKine(kFALSE),
  fEsdTrackCuts(0),
  fdNdEtaAnalysisESD(0),
  fMult(0),
  fMultVtx(0),
  fEvents(0),
  fVertexResolution(0),
  fdNdEtaAnalysis(0),
  fdNdEtaAnalysisND(0),
  fdNdEtaAnalysisNSD(0),
  fdNdEtaAnalysisTr(0),
  fdNdEtaAnalysisTrVtx(0),
  fdNdEtaAnalysisTracks(0),
  fPartPt(0),
  fVertex(0),
  fPhi(0),
  fRawPt(0),
  fEtaPhi(0),
  fDeltaPhi(0),
  fFiredChips(0),
  fTriggerVsTime(0),
  fStats(0)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
  
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

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

Bool_t AlidNdEtaTask::Notify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file: %s", count, ((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  return kTRUE;
}

//________________________________________________________________________
void AlidNdEtaTask::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  Printf("AlidNdEtaTask::ConnectInputData called");

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();

    // Enable only the needed branches
    esdH->SetActiveBranches("AliESDHeader Vertex");

    if (fAnalysisMode == AliPWG0Helper::kSPD)
      esdH->SetActiveBranches("AliESDHeader Vertex AliMultiplicity");

    if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS) {
      esdH->SetActiveBranches("AliESDHeader Vertex Tracks");
    }
  }

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
}

void AlidNdEtaTask::CreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AlidNdEtaTask::CreateOutputObjects");

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
    fPartEta[i] = new TH1F(Form("dndeta_check_%d", i), Form("dndeta_check_%d", i), 60, -6, 6);
    fPartEta[i]->Sumw2();
    fOutput->Add(fPartEta[i]);
  }

  fEvents = new TH1F("dndeta_check_vertex", "dndeta_check_vertex", 800, -40, 40);
  fOutput->Add(fEvents);

  Float_t resMax = (fAnalysisMode == AliPWG0Helper::kSPD) ? 0.2 : 2;
  fVertexResolution = new TH1F("dndeta_vertex_resolution_z", "dndeta_vertex_resolution_z", 1000, 0, resMax);
  fOutput->Add(fVertexResolution);

  fPhi = new TH1F("fPhi", "fPhi;#phi in rad.;count", 720, 0, 2 * TMath::Pi());
  fOutput->Add(fPhi);

  fEtaPhi = new TH2F("fEtaPhi", "fEtaPhi;#eta;#phi in rad.;count", 80, -4, 4, 18*5, 0, 2 * TMath::Pi());
  fOutput->Add(fEtaPhi);

  fTriggerVsTime = new TGraph; //TH1F("fTriggerVsTime", "fTriggerVsTime;trigger time;count", 100, 0, 100);
  fTriggerVsTime->SetName("fTriggerVsTime");
  fTriggerVsTime->GetXaxis()->SetTitle("trigger time");
  fTriggerVsTime->GetYaxis()->SetTitle("count");
  fOutput->Add(fTriggerVsTime);

  fStats = new TH1F("fStats", "fStats", 2, 0.5, 2.5);
  fStats->GetXaxis()->SetBinLabel(1, "vertexer 3d");
  fStats->GetXaxis()->SetBinLabel(2, "vertexer z");
  fOutput->Add(fStats);

  if (fAnalysisMode == AliPWG0Helper::kSPD)
  {
    fDeltaPhi = new TH1F("fDeltaPhi", "fDeltaPhi;#Delta #phi;Entries", 500, -0.2, 0.2);
    fOutput->Add(fDeltaPhi);
    fFiredChips = new TH2F("fFiredChips", "fFiredChips;Chips L1 + L2;tracklets", 1201, -0.5, 1201, 50, -0.5, 49.5);
    fOutput->Add(fFiredChips);
    for (Int_t i=0; i<2; i++)
    {
      fZPhi[i] = new TH2F(Form("fZPhi_%d", i), Form("fZPhi Layer %d;z (cm);#phi (rad.)", i), 200, -20, 20, 180, 0, TMath::Pi() * 2);
      fOutput->Add(fZPhi[i]);
    }
  }

  if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS)
  {
    fRawPt =  new TH1F("fRawPt", "raw pt;p_{T};Count", 2000, 0, 100);
    fOutput->Add(fRawPt);
  }

  fVertex = new TH3F("vertex_check", "vertex_check", 100, -1, 1, 100, -1, 1, 100, -30, 30);
  fOutput->Add(fVertex);

  if (fReadMC)
  {
    fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysis);

    fdNdEtaAnalysisND = new dNdEtaAnalysis("dndetaND", "dndetaND", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisND);

    fdNdEtaAnalysisNSD = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD", fAnalysisMode);
    fOutput->Add(fdNdEtaAnalysisNSD);

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
}

void AlidNdEtaTask::Exec(Option_t*)
{
  // process the event

  // these variables are also used in the MC section below; however, if ESD is off they stay with their default values
  Bool_t eventTriggered = kFALSE;
  const AliESDVertex* vtxESD = 0;

  // post the data already here
  PostData(0, fOutput);

  // ESD analysis
  if (fESD)
  {
    // check event type (should be PHYSICS = 7)
    AliESDHeader* esdHeader = fESD->GetHeader();
    if (!esdHeader)
    {
      Printf("ERROR: esdHeader could not be retrieved");
      return;
    }

    /*
    UInt_t eventType = esdHeader->GetEventType();
    if (eventType != 7)
    {
      Printf("Skipping event because it is of type %d", eventType);
      return;
    }
    */

    // trigger definition
    eventTriggered = AliPWG0Helper::IsEventTriggered(fESD, fTrigger);

    // get the ESD vertex
    vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);

    Double_t vtx[3];

    // fill z vertex resolution
    if (vtxESD)
    {
      fVertexResolution->Fill(vtxESD->GetZRes());
      fVertex->Fill(vtxESD->GetXv(), vtxESD->GetYv(), vtxESD->GetZv());

      if (AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
      {
          vtxESD->GetXYZ(vtx);

          // vertex stats
          if (strcmp(vtxESD->GetTitle(), "vertexer: 3D") == 0)
          {
            fStats->Fill(1);
          }
          else if (strcmp(vtxESD->GetTitle(), "vertexer: Z") == 0)
            fStats->Fill(2);
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
    if (fAnalysisMode == AliPWG0Helper::kSPD)
    {
      // get tracklets
      const AliMultiplicity* mult = fESD->GetMultiplicity();
      if (!mult)
      {
        AliDebug(AliLog::kError, "AliMultiplicity not available");
        return;
      }

      labelArr = new Int_t[mult->GetNumberOfTracklets()];
      etaArr = new Float_t[mult->GetNumberOfTracklets()];
      thirdDimArr = new Float_t[mult->GetNumberOfTracklets()];

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
        fPhi->Fill(phi);
        fEtaPhi->Fill(eta, phi);
        
        if (deltaPhi < 0.01)
        {
          // layer 0
          Float_t z = vtx[2] + 3.9 / TMath::Tan(2 * TMath::ATan(TMath::Exp(- eta)));
          fZPhi[0]->Fill(z, phi);
          // layer 1
          z = vtx[2] + 7.6 / TMath::Tan(2 * TMath::ATan(TMath::Exp(- eta)));
          fZPhi[1]->Fill(z, phi);
        }

        fDeltaPhi->Fill(deltaPhi);

        if (fDeltaPhiCut > 0 && TMath::Abs(deltaPhi) > fDeltaPhiCut)
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
        
        etaArr[inputCount] = eta;
        labelArr[inputCount] = label;
        thirdDimArr[inputCount] = phi;
        ++inputCount;
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
    }
    else if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS)
    {
      if (!fEsdTrackCuts)
      {
        AliDebug(AliLog::kError, "fESDTrackCuts not available");
        return;
      }

      if (vtxESD)
      {
        // get multiplicity from ESD tracks
        TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD, (fAnalysisMode == AliPWG0Helper::kTPC));
        Int_t nGoodTracks = list->GetEntries();
        Printf("Accepted %d tracks", nGoodTracks);
  
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
  
          fPhi->Fill(phi);
          fEtaPhi->Fill(eta, phi);
          if (eventTriggered && vtxESD)
            fRawPt->Fill(pT);
  
          if (fOnlyPrimaries && label == 0)
            continue;
  
          if (fUseMCKine)
          {
            if (label > 0)
            {
              TParticle* particle = stack->Particle(label);
              eta = particle->Eta();
              pT = particle->Pt();
            }
            else
              Printf("WARNING: fUseMCKine set without fOnlyPrimaries and no label found");
          }
  
          etaArr[inputCount] = eta;
          labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
          thirdDimArr[inputCount] = pT;
          ++inputCount;
        }
        
        // TODO restrict inputCount used as measure for the multiplicity to |eta| < 1
  
        delete list;
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

      fTriggerVsTime->SetPoint(fTriggerVsTime->GetN(), fESD->GetTimeStamp(), 1);

      if (vtxESD)
      {
        // control hist
        fMultVtx->Fill(inputCount);

        for (Int_t i=0; i<inputCount; ++i)
        {
          Float_t eta = etaArr[i];
          Float_t thirdDim = thirdDimArr[i];

          fdNdEtaAnalysisESD->FillTrack(vtx[2], eta, thirdDim);

          if (TMath::Abs(vtx[2]) < 10)
          {
            fPartEta[0]->Fill(eta);

            if (vtx[2] < 0)
              fPartEta[1]->Fill(eta);
            else
              fPartEta[2]->Fill(eta);
          }
        }

        // for event count per vertex
        fdNdEtaAnalysisESD->FillEvent(vtx[2], inputCount);

        // control hist
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
            if (fAnalysisMode == AliPWG0Helper::kSPD)
            {
              if (fFillPhi)
              {
                thirdDim = particle->Phi();
              }
              else
                thirdDim = inputCount;
            }
            else
              thirdDim = particle->Pt();

            fdNdEtaAnalysisTracks->FillTrack(vtxMC[2], particle->Eta(), thirdDim);
          } // end of track loop

          // for event count per vertex
          fdNdEtaAnalysisTracks->FillEvent(vtxMC[2], inputCount);
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
    Int_t processType = AliPWG0Helper::GetEventProcessType(header);
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", processType));

    if (processType==AliPWG0Helper::kInvalidProcess)
      AliDebug(AliLog::kError, Form("Unknown process type %d.", processType));

    // loop over mc particles
    Int_t nPrim  = stack->GetNprimary();

    Int_t nAcceptedParticles = 0;

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

       // make a rough eta cut (so that nAcceptedParticles is not too far off the true measured value (NB: this histograms are only gathered for comparison))
      if (TMath::Abs(eta) < 1.5) // && pt > 0.3)
        nAcceptedParticles++;
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
      Float_t thirdDim = -1;

      if (fAnalysisMode == AliPWG0Helper::kSPD)
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

      if (eventTriggered)
      {
        fdNdEtaAnalysisTr->FillTrack(vtxMC[2], eta, thirdDim);
        if (vtxESD)
          fdNdEtaAnalysisTrVtx->FillTrack(vtxMC[2], eta, thirdDim);
      }

      if (TMath::Abs(eta) < 1.0)
        fPartPt->Fill(particle->Pt());
    }

    fdNdEtaAnalysis->FillEvent(vtxMC[2], nAcceptedParticles);
    if (processType != AliPWG0Helper::kSD)
      fdNdEtaAnalysisNSD->FillEvent(vtxMC[2], nAcceptedParticles);
    if (processType == AliPWG0Helper::kND)
      fdNdEtaAnalysisND->FillEvent(vtxMC[2], nAcceptedParticles);

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

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput)
    Printf("ERROR: fOutput not available");

  if (fOutput)
  {
    fdNdEtaAnalysisESD = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("fdNdEtaAnalysisESD"));
    fMult = dynamic_cast<TH1F*> (fOutput->FindObject("fMult"));
    fMultVtx = dynamic_cast<TH1F*> (fOutput->FindObject("fMultVtx"));
    for (Int_t i=0; i<3; ++i)
      fPartEta[i] = dynamic_cast<TH1F*> (fOutput->FindObject(Form("dndeta_check_%d", i)));
    fEvents = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_vertex"));
    fVertexResolution = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_vertex_resolution_z"));

    fPhi = dynamic_cast<TH1F*> (fOutput->FindObject("fPhi"));
    fRawPt = dynamic_cast<TH1F*> (fOutput->FindObject("fRawPt"));
    fEtaPhi = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaPhi"));
    for (Int_t i=0; i<2; ++i)
      fZPhi[i] = dynamic_cast<TH2F*> (fOutput->FindObject(Form("fZPhi_%d", i)));
    fDeltaPhi = dynamic_cast<TH1F*> (fOutput->FindObject("fDeltaPhi"));
    fFiredChips = dynamic_cast<TH2F*> (fOutput->FindObject("fFiredChips"));
    fTriggerVsTime = dynamic_cast<TGraph*> (fOutput->FindObject("fTriggerVsTime"));
    fStats = dynamic_cast<TH1F*> (fOutput->FindObject("fStats"));

    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCuts"));
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

    if (fPartEta[0])
    {
      Int_t events1 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(-19.9), fEvents->GetXaxis()->FindBin(-0.001));
      Int_t events2 = (Int_t) fEvents->Integral(fEvents->GetXaxis()->FindBin(0.001), fEvents->GetXaxis()->FindBin(19.9));

      Printf("%d events with vertex used", events1 + events2);

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

    TFile* fout = new TFile("analysis_esd_raw.root", "RECREATE");

    if (fdNdEtaAnalysisESD)
      fdNdEtaAnalysisESD->SaveHistograms();

    if (fEsdTrackCuts)
      fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

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
      fVertexResolution->Write();

    if (fDeltaPhi)
      fDeltaPhi->Write();

    if (fPhi)
      fPhi->Write();

    if (fRawPt)
      fRawPt->Write();

    if (fEtaPhi)
      fEtaPhi->Write();

    for (Int_t i=0; i<2; ++i)
      if (fZPhi[i])
        fZPhi[i]->Write();
    
    if (fFiredChips)
      fFiredChips->Write();

    if (fTriggerVsTime)
      fTriggerVsTime->Write();

    if (fStats)
      fStats->Write();

    fVertex = dynamic_cast<TH3F*> (fOutput->FindObject("vertex_check"));
    if (fVertex)
      fVertex->Write();

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
